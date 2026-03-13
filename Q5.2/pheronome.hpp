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
               int rank, int size,
               double alpha = 0.7, double beta = 0.9999 )
        : m_dim( dim ),
          m_stride( dim + 2 ),
          m_alpha(alpha), m_beta(beta),
          m_pos_nest(pos_nest), m_pos_food(pos_food), m_rank(rank), m_size(size) 
          {
        size_t rows_per_proc = dim / size;
        size_t remainder = dim % size;
        m_row_start = rank * rows_per_proc + std::min((size_t)rank, remainder);
        size_t row_end = m_row_start + rows_per_proc + (rank < (int)remainder ? 1 : 0);
        m_local_height = row_end - m_row_start;

        m_map_of_pheronome.assign((m_local_height + 2) * m_stride, {{0., 0.}});

        if (is_local(pos_food.y)) (*this)(pos_food.x, pos_food.y) = {{1.0, 0.0}};
        if (is_local(pos_nest.y)) (*this)(pos_nest.x, pos_nest.y) = {{0.0, 1.0}};

        cl_update();
        m_buffer_pheronome = m_map_of_pheronome;
    }
    pheronome( const pheronome& ) = delete;
    pheronome( pheronome&& )      = delete;
    ~pheronome( )                 = default;

    pheronome_t& operator( )( size_t i, size_t j_global ) {
        size_t j_local = j_global - m_row_start + 1;
        return m_map_of_pheronome[j_local * m_stride + (i + 1)];
    }

    const pheronome_t& operator( )( size_t i, size_t j_global ) const {
        size_t j_local = j_global - m_row_start + 1;
        return m_map_of_pheronome[j_local * m_stride + (i + 1)];
    }

    bool is_local(size_t j_global) const {
        return j_global >= m_row_start && j_global < m_row_start + m_local_height;
    }

    size_t local_height() const { return m_local_height; }
    size_t row_start() const { return m_row_start; }
    size_t stride() const { return m_stride; }
    pheronome_t* data() { return m_map_of_pheronome.data(); }

    void sync_halos() {
        int north = (m_rank > 0) ? m_rank - 1 : MPI_PROC_NULL;
        int south = (m_rank < m_size - 1) ? m_rank + 1 : MPI_PROC_NULL;

        MPI_Sendrecv(&m_map_of_pheronome[1 * m_stride], m_stride * 2, MPI_DOUBLE, north, 0,
                     &m_map_of_pheronome[0 * m_stride], m_stride * 2, MPI_DOUBLE, north, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&m_map_of_pheronome[m_local_height * m_stride], m_stride * 2, MPI_DOUBLE, south, 0,
                     &m_map_of_pheronome[(m_local_height + 1) * m_stride], m_stride * 2, MPI_DOUBLE, south, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    void do_evaporation( ) {
       #pragma omp parallel for collapse(2)
        for (size_t j = 1; j <= m_local_height; ++j) {
            for (size_t i = 1; i <= m_dim; ++i) {
                m_map_of_pheronome[j * m_stride + i][0] *= m_beta;
                m_map_of_pheronome[j * m_stride + i][1] *= m_beta;
            }
        }
    }

    void update( ) {
        m_map_of_pheronome.swap(m_buffer_pheronome);
        sync_halos(); 
        cl_update();
        if (is_local(m_pos_food.y)) (*this)(m_pos_food.x, m_pos_food.y)[0] = 1.0;
        if (is_local(m_pos_nest.y)) (*this)(m_pos_nest.x, m_pos_nest.y)[1] = 1.0;
    }

private:
    /**
     * @brief Mets à jour les conditions limites sur les cellules fantômes
     * @details Mets à jour les conditions limites sur les cellules fantômes :
     *     pour l'instant, on se contente simplement de mettre ces cellules avec
     *     des valeurs à -1 pour être sûr que les fourmis évitent ces cellules
     */
    void cl_update( ) {
        // On mets tous les bords à -1 pour les marquer comme indésirables :
        for (size_t j = 0; j < m_local_height + 2; ++j) {
            m_map_of_pheronome[j * m_stride] = {{-1., -1.}};
            m_map_of_pheronome[j * m_stride + m_dim + 1] = {{-1., -1.}};
        }

        if (m_rank == 0) {
            for (size_t i = 0; i < m_stride; ++i) m_map_of_pheronome[i] = {{-1., -1.}};
        }
        if (m_rank == m_size - 1) {
            for (size_t i = 0; i < m_stride; ++i) 
                m_map_of_pheronome[(m_local_height + 1) * m_stride + i] = {{-1., -1.}};
        }
    }
    size_t m_dim, m_stride, m_row_start, m_local_height;
    double m_alpha, m_beta;
    int m_rank, m_size;
    std::vector<pheronome_t> m_map_of_pheronome, m_buffer_pheronome;
    position_t m_pos_nest, m_pos_food;
};

#endif
