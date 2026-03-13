#ifndef _PHERONOME_HPP_
#define _PHERONOME_HPP_

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <mpi.h>

#include "basic_types.hpp"

class pheronome {
public:
    using size_t = unsigned long;
    using pheronome_t = std::array<double, 2>;

    pheronome(size_t dim,
              const position_t& pos_food,
              const position_t& pos_nest,
              int rank,
              int size,
              double alpha = 0.7,
              double beta = 0.9999)
        : m_dim(dim),
          m_stride(dim + 2),
          m_rank(rank),
          m_size(size),
          m_alpha(alpha),
          m_beta(beta),
          m_pos_nest(pos_nest),
          m_pos_food(pos_food) {
        if (size <= 0) {
            throw std::runtime_error("Invalid MPI size in pheronome");
        }
        if (dim < 3) {
            throw std::runtime_error("Invalid pheromone dimensions");
        }
        if (static_cast<size_t>(size) > dim) {
            throw std::runtime_error("More MPI ranks than grid rows in pheronome");
        }

        const size_t rows_per_rank = dim / static_cast<size_t>(size);
        const size_t remainder = dim % static_cast<size_t>(size);
        m_row_start = static_cast<size_t>(rank) * rows_per_rank +
                      std::min(static_cast<size_t>(rank), remainder);
        m_local_height = rows_per_rank + (static_cast<size_t>(rank) < remainder ? 1UL : 0UL);

        if (m_local_height == 0) {
            throw std::runtime_error("local_height is zero in pheronome");
        }

        m_map.assign((m_local_height + 2) * m_stride, {{0.0, 0.0}});
        m_buffer = m_map;
        set_boundaries(m_map);
        set_boundaries(m_buffer);

        if (is_local_row(m_pos_food.y)) {
            (*this)(static_cast<size_t>(m_pos_food.x), static_cast<size_t>(m_pos_food.y))[0] = 1.0;
        }
        if (is_local_row(m_pos_nest.y)) {
            (*this)(static_cast<size_t>(m_pos_nest.x), static_cast<size_t>(m_pos_nest.y))[1] = 1.0;
        }
        m_buffer = m_map;
    }

    pheronome(const pheronome&) = delete;
    pheronome(pheronome&&) = delete;
    ~pheronome() = default;

    size_t dimensions() const { return m_dim; }
    size_t stride() const { return m_stride; }
    size_t row_start() const { return m_row_start; }
    size_t row_end() const { return m_row_start + m_local_height; }
    size_t local_height() const { return m_local_height; }

    bool is_local_row(int y_global) const {
        return y_global >= static_cast<int>(m_row_start) &&
               y_global < static_cast<int>(m_row_start + m_local_height);
    }

    pheronome_t& operator()(size_t x_global, size_t y_global) {
        return m_map[local_index(x_global, y_global)];
    }

    const pheronome_t& operator()(size_t x_global, size_t y_global) const {
        return m_map[local_index(x_global, y_global)];
    }

    double pheromone_value(int x_global, int y_global, int channel) const {
        if (channel < 0 || channel > 1) {
            return -1.0;
        }
        if (x_global < 0 || y_global < 0) {
            return -1.0;
        }
        if (x_global >= static_cast<int>(m_dim) || y_global >= static_cast<int>(m_dim)) {
            return -1.0;
        }

        const long long y_ll = static_cast<long long>(y_global);
        const long long start_ll = static_cast<long long>(m_row_start);
        const long long stop_ll = static_cast<long long>(m_row_start + m_local_height - 1);
        if (y_ll < start_ll - 1 || y_ll > stop_ll + 1) {
            return -1.0;
        }

        const size_t y_local = static_cast<size_t>(y_ll - start_ll + 1);
        const size_t x_local = static_cast<size_t>(x_global + 1);
        return m_map[y_local * m_stride + x_local][channel];
    }

    void do_evaporation() {
        for (size_t y = 1; y <= m_local_height; ++y) {
            for (size_t x = 1; x <= m_dim; ++x) {
                m_buffer[y * m_stride + x][0] *= m_beta;
                m_buffer[y * m_stride + x][1] *= m_beta;
            }
        }
    }

    void mark_pheronome_xy(size_t x_global, size_t y_global) {
        if (x_global >= m_dim || y_global >= m_dim) {
            return;
        }
        if (!is_local_row(static_cast<int>(y_global))) {
            return;
        }

        const size_t y_local = static_cast<size_t>(static_cast<long long>(y_global) -
                                                   static_cast<long long>(m_row_start) + 1LL);
        const size_t x_local = x_global + 1;

        auto value_or_zero = [&](int dx, int dy, int channel) -> double {
            const size_t xx = static_cast<size_t>(static_cast<int>(x_local) + dx);
            const size_t yy = static_cast<size_t>(static_cast<int>(y_local) + dy);
            return std::max(0.0, m_map[yy * m_stride + xx][channel]);
        };

        for (int c = 0; c < 2; ++c) {
            const double left = value_or_zero(-1, 0, c);
            const double right = value_or_zero(1, 0, c);
            const double up = value_or_zero(0, -1, c);
            const double down = value_or_zero(0, 1, c);
            const double vmax = std::max({left, right, up, down});
            const double vavg = 0.25 * (left + right + up + down);
            m_buffer[y_local * m_stride + x_local][c] = m_alpha * vmax + (1.0 - m_alpha) * vavg;
        }
    }

    void update() {
        m_map.swap(m_buffer);
        set_boundaries(m_map);

        if (is_local_row(m_pos_food.y)) {
            (*this)(static_cast<size_t>(m_pos_food.x), static_cast<size_t>(m_pos_food.y))[0] = 1.0;
        }
        if (is_local_row(m_pos_nest.y)) {
            (*this)(static_cast<size_t>(m_pos_nest.x), static_cast<size_t>(m_pos_nest.y))[1] = 1.0;
        }
    }

    void exchange_halos(double& t_comm_halo_ms) {
        const int up = (m_rank == 0) ? MPI_PROC_NULL : m_rank - 1;
        const int down = (m_rank == m_size - 1) ? MPI_PROC_NULL : m_rank + 1;
        const int row_count = static_cast<int>(m_stride * 2);

        double t0 = MPI_Wtime();
        MPI_Sendrecv(reinterpret_cast<double*>(m_map.data() + m_stride),
                     row_count,
                     MPI_DOUBLE,
                     up,
                     100,
                     reinterpret_cast<double*>(m_map.data()),
                     row_count,
                     MPI_DOUBLE,
                     up,
                     101,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

        MPI_Sendrecv(reinterpret_cast<double*>(m_map.data() + m_local_height * m_stride),
                     row_count,
                     MPI_DOUBLE,
                     down,
                     101,
                     reinterpret_cast<double*>(m_map.data() + (m_local_height + 1) * m_stride),
                     row_count,
                     MPI_DOUBLE,
                     down,
                     100,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        t_comm_halo_ms += (MPI_Wtime() - t0) * 1000.0;

        set_vertical_boundaries(m_map);
        if (m_rank == 0) {
            set_top_halo_boundaries(m_map);
        }
        if (m_rank == m_size - 1) {
            set_bottom_halo_boundaries(m_map);
        }
    }

    void pack_interior(std::vector<double>& out) const {
        out.clear();
        out.reserve(static_cast<std::size_t>(m_local_height) * m_dim * 2ULL);
        for (size_t y_global = row_start(); y_global < row_end(); ++y_global) {
            for (size_t x = 0; x < m_dim; ++x) {
                const auto& cell = (*this)(x, y_global);
                out.push_back(cell[0]);
                out.push_back(cell[1]);
            }
        }
    }

private:
    size_t local_index(size_t x_global, size_t y_global) const {
        if (x_global >= m_dim || y_global >= m_dim) {
            throw std::out_of_range("pheronome index out of global bounds");
        }
        const long long y_ll = static_cast<long long>(y_global);
        const long long start_ll = static_cast<long long>(m_row_start);
        const long long stop_ll = static_cast<long long>(m_row_start + m_local_height - 1);
        if (y_ll < start_ll - 1 || y_ll > stop_ll + 1) {
            throw std::out_of_range("pheronome y out of local + halo bounds");
        }
        const size_t y_local = static_cast<size_t>(y_ll - start_ll + 1);
        return y_local * m_stride + (x_global + 1);
    }

    void set_vertical_boundaries(std::vector<pheronome_t>& map) const {
        for (size_t y = 0; y < m_local_height + 2; ++y) {
            map[y * m_stride] = {{-1.0, -1.0}};
            map[y * m_stride + (m_dim + 1)] = {{-1.0, -1.0}};
        }
    }

    void set_top_halo_boundaries(std::vector<pheronome_t>& map) const {
        for (size_t x = 0; x < m_stride; ++x) {
            map[x] = {{-1.0, -1.0}};
        }
    }

    void set_bottom_halo_boundaries(std::vector<pheronome_t>& map) const {
        const size_t base = (m_local_height + 1) * m_stride;
        for (size_t x = 0; x < m_stride; ++x) {
            map[base + x] = {{-1.0, -1.0}};
        }
    }

    void set_boundaries(std::vector<pheronome_t>& map) const {
        set_vertical_boundaries(map);
        if (m_rank == 0) {
            set_top_halo_boundaries(map);
        }
        if (m_rank == m_size - 1) {
            set_bottom_halo_boundaries(map);
        }
    }

    size_t m_dim = 0;
    size_t m_stride = 0;
    int m_rank = 0;
    int m_size = 1;
    size_t m_row_start = 0;
    size_t m_local_height = 0;
    double m_alpha = 0.7;
    double m_beta = 0.9999;
    std::vector<pheronome_t> m_map;
    std::vector<pheronome_t> m_buffer;
    position_t m_pos_nest{};
    position_t m_pos_food{};
};

#endif
