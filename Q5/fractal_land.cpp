#include "fractal_land.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "rand_generator.hpp"

namespace {

std::vector<double> generate_global_altitude(fractal_land::dim_t ln2_dim,
                                             unsigned long nbSeeds,
                                             double deviation,
                                             int seed,
                                             unsigned long dim) {
    std::vector<double> altitude(dim * dim, 0.0);

    const unsigned long dim_ss_grid = 1UL << ln2_dim;
    auto get_global = [&](unsigned long x, unsigned long y) -> double& {
        return altitude[x + y * dim];
    };

    RandomGenerator gen(seed, 0.0, dim_ss_grid * deviation);
    for (fractal_land::dim_t x = 0; x < dim; x += dim_ss_grid) {
        for (fractal_land::dim_t y = 0; y < dim; y += dim_ss_grid) {
            get_global(x, y) = gen(x, y);
        }
    }

    auto compute_subgrid = [&](int log_subgrid_dim, int iB, int jB, double dev, std::size_t local_seed) {
        RandomGenerator gen_step(local_seed, -dev, dev);
        const unsigned long sub_dim = 1UL << log_subgrid_dim;
        const unsigned long xBeg = static_cast<unsigned long>(iB) * sub_dim;
        const unsigned long yBeg = static_cast<unsigned long>(jB) * sub_dim;
        const int mid = static_cast<int>(sub_dim / 2UL);
        const int xMid = static_cast<int>(xBeg) + mid;
        const int yMid = static_cast<int>(yBeg) + mid;
        const int xEnd = static_cast<int>(xBeg + sub_dim);
        const int yEnd = static_cast<int>(yBeg + sub_dim);

        get_global(static_cast<unsigned long>(xMid), yBeg) =
            0.5 * (get_global(xBeg, yBeg) + get_global(static_cast<unsigned long>(xEnd), yBeg)) +
            mid * gen_step(xMid, static_cast<int>(yBeg));
        get_global(xBeg, static_cast<unsigned long>(yMid)) =
            0.5 * (get_global(xBeg, yBeg) + get_global(xBeg, static_cast<unsigned long>(yEnd))) +
            mid * gen_step(static_cast<int>(xBeg), yMid);
        get_global(static_cast<unsigned long>(xMid), static_cast<unsigned long>(yEnd)) =
            0.5 * (get_global(xBeg, static_cast<unsigned long>(yEnd)) +
                   get_global(static_cast<unsigned long>(xEnd), static_cast<unsigned long>(yEnd))) +
            mid * gen_step(xMid, yEnd);
        get_global(static_cast<unsigned long>(xEnd), static_cast<unsigned long>(yMid)) =
            0.5 * (get_global(static_cast<unsigned long>(xEnd), yBeg) +
                   get_global(static_cast<unsigned long>(xEnd), static_cast<unsigned long>(yEnd))) +
            mid * gen_step(xEnd, yMid);
        get_global(static_cast<unsigned long>(xMid), static_cast<unsigned long>(yMid)) =
            0.25 * (get_global(static_cast<unsigned long>(xMid), yBeg) +
                    get_global(xBeg, static_cast<unsigned long>(yMid)) +
                    get_global(static_cast<unsigned long>(xMid), static_cast<unsigned long>(yEnd)) +
                    get_global(static_cast<unsigned long>(xEnd), static_cast<unsigned long>(yMid))) +
            mid * gen_step(xMid, yMid);
    };

    fractal_land::dim_t level = ln2_dim;
    unsigned long current_nbSeeds = nbSeeds;
    while (level > 1) {
        level -= 1;
        current_nbSeeds *= 2;
        for (unsigned long iB = 0; iB < current_nbSeeds; ++iB) {
            for (unsigned long jB = 0; jB < current_nbSeeds; ++jB) {
                compute_subgrid(static_cast<int>(level),
                                static_cast<int>(iB),
                                static_cast<int>(jB),
                                deviation,
                                static_cast<std::size_t>(seed));
            }
        }
    }

    return altitude;
}

} // namespace

fractal_land::fractal_land(const dim_t& log_size,
                           unsigned long nbSeeds,
                           double deviation,
                           int seed,
                           int rank,
                           int size)
    : m_dimensions(0), m_row_start(0), m_local_height(0) {
    if (size <= 0) {
        throw std::runtime_error("Invalid MPI size in fractal_land");
    }

    const unsigned long dim_ss_grid = 1UL << log_size;
    m_dimensions = nbSeeds * dim_ss_grid + 1UL;
    if (m_dimensions < 3) {
        throw std::runtime_error("Invalid fractal dimensions");
    }
    if (static_cast<unsigned long>(size) > m_dimensions) {
        throw std::runtime_error("More MPI ranks than grid rows in fractal_land");
    }

    const unsigned long rows_per_rank = m_dimensions / static_cast<unsigned long>(size);
    const unsigned long remainder = m_dimensions % static_cast<unsigned long>(size);
    m_row_start = static_cast<unsigned long>(rank) * rows_per_rank +
                  std::min(static_cast<unsigned long>(rank), remainder);
    m_local_height = rows_per_rank + (static_cast<unsigned long>(rank) < remainder ? 1UL : 0UL);
    if (m_local_height == 0) {
        throw std::runtime_error("local_height is zero in fractal_land");
    }

    const auto global_altitude =
        generate_global_altitude(log_size, nbSeeds, deviation, seed, m_dimensions);
    fill_from_global(global_altitude);
}

fractal_land::fractal_land(const dim_t& log_size,
                           unsigned long nbSeeds,
                           double deviation,
                           int seed)
    : fractal_land(log_size, nbSeeds, deviation, seed, 0, 1) {}

std::size_t fractal_land::local_index(unsigned long x_global, unsigned long y_global) const {
    if (x_global >= m_dimensions || y_global >= m_dimensions) {
        throw std::out_of_range("fractal_land global index out of bounds");
    }
    const long long y_ll = static_cast<long long>(y_global);
    const long long start_ll = static_cast<long long>(m_row_start);
    const long long stop_ll = static_cast<long long>(m_row_start + m_local_height - 1);
    if (y_ll < start_ll - 1 || y_ll > stop_ll + 1) {
        throw std::out_of_range("fractal_land y out of local + halo bounds");
    }
    const std::size_t y_local = static_cast<std::size_t>(y_ll - start_ll + 1);
    return static_cast<std::size_t>(x_global) + y_local * static_cast<std::size_t>(m_dimensions);
}

double fractal_land::operator()(unsigned long x_global, unsigned long y_global) const {
    return m_altitude[local_index(x_global, y_global)];
}

double& fractal_land::operator()(unsigned long x_global, unsigned long y_global) {
    return m_altitude[local_index(x_global, y_global)];
}

double fractal_land::value_at(int x_global, int y_global) const {
    if (x_global < 0 || y_global < 0) {
        return 1.0;
    }
    if (x_global >= static_cast<int>(m_dimensions) || y_global >= static_cast<int>(m_dimensions)) {
        return 1.0;
    }

    const long long y_ll = static_cast<long long>(y_global);
    const long long start_ll = static_cast<long long>(m_row_start);
    const long long stop_ll = static_cast<long long>(m_row_start + m_local_height - 1);
    if (y_ll < start_ll - 1 || y_ll > stop_ll + 1) {
        return 1.0;
    }

    return m_altitude[local_index(static_cast<unsigned long>(x_global),
                                  static_cast<unsigned long>(y_global))];
}

void fractal_land::fill_from_global(const std::vector<double>& global_altitude) {
    m_altitude.assign(static_cast<std::size_t>(m_local_height + 2) * m_dimensions, 0.0);

    auto global_at = [&](unsigned long x, unsigned long y) -> double {
        return global_altitude[x + y * m_dimensions];
    };

    for (unsigned long j = 0; j < m_local_height; ++j) {
        const unsigned long y_global = m_row_start + j;
        for (unsigned long x = 0; x < m_dimensions; ++x) {
            m_altitude[x + (j + 1) * m_dimensions] = global_at(x, y_global);
        }
    }

    const unsigned long y_top = (m_row_start > 0) ? (m_row_start - 1) : m_row_start;
    const unsigned long y_bottom = (m_row_start + m_local_height < m_dimensions)
                                       ? (m_row_start + m_local_height)
                                       : (m_row_start + m_local_height - 1);

    for (unsigned long x = 0; x < m_dimensions; ++x) {
        m_altitude[x] = global_at(x, y_top);
        m_altitude[x + (m_local_height + 1) * m_dimensions] = global_at(x, y_bottom);
    }
}
