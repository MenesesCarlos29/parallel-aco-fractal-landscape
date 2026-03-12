#ifndef _PHERONOME_HPP_
#define _PHERONOME_HPP_
#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <omp.h>
#include <utility>
#include <vector>
#include "basic_types.hpp"
#include <mpi.h>

/**
 * @brief Carte des phéronomes
 * @details Gère une carte des phéronomes avec leurs mis à jour ( dont l'évaporation )
 *
 */
class pheronome {
public:
    using size_t      = unsigned long;
    using pheronome_t = std::array< double, 2 >;

    /**
     * @brief Construit une carte initiale des phéronomes
     * @details La carte des phéronomes est initialisées à zéro ( neutre )
     *          sauf pour les bords qui sont marqués comme indésirables
     *
     * @param dim Nombre de cellule dans chaque direction
     * @param alpha Paramètre de bruit
     * @param beta Paramêtre d'évaporation
     */
    pheronome( size_t dim, const position_t& pos_food, const position_t& pos_nest,
               double alpha = 0.7, double beta = 0.9999, int rank, int size )
        : m_dim(dim), m_stride(dim + 2), m_alpha(alpha), m_beta(beta),
          m_pos_nest(pos_nest), m_pos_food(pos_food)
          {
        size_t rows_per_proc = dim / size;
        size_t remainder = dim % size;
        m_row_start = rank * rows_per_proc + std::min((size_t)rank, remainder);
        m_local_height = rows_per_proc + ((size_t)rank < remainder ? 1 : 0);

        m_map_of_pheronome.assign(m_stride * (m_local_height + 2), {{0., 0.}});

        for (size_t i = 0; i < m_local_height + 2; ++i) {
            m_map_of_pheronome[i * m_stride] = {{-1., -1.}};
            m_map_of_pheronome[i * m_stride + m_dim + 1] = {{-1., -1.}};
        }

        if (rank == 0) 
            for(size_t j=0; j<m_stride; ++j) m_map_of_pheronome[j] = {{-1.,-1.}};
        if (rank == size - 1)
            for(size_t j=0; j<m_stride; ++j) m_map_of_pheronome[j + (m_local_height+1)*m_stride] = {{-1.,-1.}};

        if (is_local(pos_food.y)) (*this)(pos_food.x, pos_food.y)[0] = 1.0;
        if (is_local(pos_nest.y)) (*this)(pos_nest.x, pos_nest.y)[1] = 1.0;

        m_buffer_pheronome = m_map_of_pheronome;    
    }
    pheronome( const pheronome& ) = delete;
    pheronome( pheronome&& )      = delete;
    ~pheronome( )                 = default;

    pheronome_t& operator( )( size_t i_global, size_t j_global ) {
        size_t j_local = j_global - m_row_start + 1; 
        return m_map_of_pheronome[j_local * m_stride + (i_global + 1)];
    }

    const pheronome_t& operator( )( size_t i_global, size_t j_global ) const {
        size_t j_local = j_global - m_row_start + 1; 
        return m_map_of_pheronome[j_local * m_stride + (i_global + 1)];
    }

    bool is_local(size_t j_global) const {
        return (j_global >= m_row_start && j_global < m_row_start + m_local_height);
    }

    void update_halos(int rank, int size) {
        int up = (rank == 0) ? MPI_PROC_NULL : rank - 1;
        int down = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;

        MPI_Sendrecv(m_map_of_pheronome.data() + m_stride, m_stride * 2, MPI_DOUBLE, up, 0,
                     m_map_of_pheronome.data(), m_stride * 2, MPI_DOUBLE, up, 0, 
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(m_map_of_pheronome.data() + m_local_height * m_stride, m_stride * 2, MPI_DOUBLE, down, 1,
                     m_map_of_pheronome.data() + (m_local_height + 1) * m_stride, m_stride * 2, MPI_DOUBLE, down, 1, 
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    void do_evaporation() {
        #pragma omp parallel for collapse(2)
        for (size_t j = 1; j <= m_local_height; ++j) {
            for (size_t i = 1; i <= m_dim; ++i) {
                m_map_of_pheronome[j * m_stride + i][0] *= m_beta;
                m_map_of_pheronome[j * m_stride + i][1] *= m_beta;
            }
        }
    }


    void mark_pheronome_xy(size_t i, size_t j_global) {
        if (!is_local(j_global)) return;
        
        size_t j_local = j_global - m_row_start + 1;

        auto get_val = [&](int di, int dj) {
            return m_map_of_pheronome[(j_local + dj) * m_stride + (i + 1 + di)];
        };
        for (int c = 0; c < 2; ++c) {
            double v_max = std::max({get_val(-1,0)[c], get_val(1,0)[c], get_val(0,-1)[c], get_val(0,1)[c]});
            double v_avg = 0.25 * (get_val(-1,0)[c] + get_val(1,0)[c] + get_val(0,-1)[c] + get_val(0,1)[c]);
            m_buffer_pheronome[j_local * m_stride + (i + 1)][c] = m_alpha * v_max + (1.0 - m_alpha) * v_avg;
        }
    }

    void update( ) {
        m_map_of_pheronome.swap( m_buffer_pheronome );

        if (is_local(m_pos_food.y)) (*this)(m_pos_food.x, m_pos_food.y)[0] = 1.0;
        if (is_local(m_pos_nest.y)) (*this)(m_pos_nest.x, m_pos_nest.y)[1] = 1.0;

        m_buffer_pheronome = m_map_of_pheronome;
    }

private:
    /**
     * @brief Mets à jour les conditions limites sur les cellules fantômes
     * @details Mets à jour les conditions limites sur les cellules fantômes :
     *     pour l'instant, on se contente simplement de mettre ces cellules avec
     *     des valeurs à -1 pour être sûr que les fourmis évitent ces cellules
     */
    unsigned long              m_dim, m_stride;
    size_t m_row_start, m_local_height;
    double                     m_alpha, m_beta;
    std::vector< pheronome_t > m_map_of_pheronome, m_buffer_pheronome;
    position_t m_pos_nest, m_pos_food;
};

#endif
