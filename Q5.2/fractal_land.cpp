# include "fractal_land.hpp"
# include "rand_generator.hpp"
#include <algorithm>

bool fractal_land::is_local(dim_t y_global) const {
    return (y_global >= m_row_start && y_global < m_row_start + m_local_height);
}

void 
fractal_land::compute_subgrid(int log_subgrid_dim, int iB, int jB, double deviation, std::size_t seed)
{
    // Génère des réels pseudo-aléatoires compris dans [-deviation;+deviation]
    RandomGenerator gen(seed, -deviation, deviation);
    unsigned long dim_ss_grid = 1UL << (log_subgrid_dim);
    unsigned long iBeg = iB * dim_ss_grid;
    unsigned long jBeg = jB * dim_ss_grid;
    int mid_ind = dim_ss_grid / 2;

    unsigned long i_mid = iBeg + mid_ind, j_mid = jBeg + mid_ind;
    unsigned long iEnd = iBeg + dim_ss_grid, jEnd = jBeg + dim_ss_grid;

    auto safe_set = [&](unsigned long i, unsigned long j, double val) {
        if (is_local(j)) (*this)(i, j) = val;
    };

    safe_set(i_mid, jBeg, 0.5 * (get_global(iBeg, jBeg) + get_global(iEnd, jBeg)) + mid_ind * gen(i_mid, jBeg));
    safe_set(iBeg, j_mid, 0.5 * (get_global(iBeg, jBeg) + get_global(iBeg, jEnd)) + mid_ind * gen(iBeg, j_mid));
    safe_set(i_mid, jEnd, 0.5 * (get_global(iBeg, jEnd) + get_global(iEnd, jEnd)) + mid_ind * gen(i_mid, jEnd));
    safe_set(iEnd, j_mid, 0.5 * (get_global(iEnd, jBeg) + get_global(iEnd, jEnd)) + mid_ind * gen(iEnd, j_mid));
    
    safe_set(i_mid, j_mid, 0.25 * (get_global(i_mid, jBeg) + get_global(iBeg, j_mid) + 
                                  get_global(i_mid, jEnd) + get_global(iEnd, j_mid)) + mid_ind * gen(i_mid, j_mid));
}

fractal_land::fractal_land(const dim_t& ln2_dim, unsigned long nbSeeds, double deviation, int seed, int rank, int size)
{
    unsigned long dim_total = (nbSeeds << ln2_dim) + 1;
    m_dimensions = dim_total;

    unsigned long rows_per_proc = dim_total / size;
    unsigned long remainder = dim_total % size;
    m_row_start = rank * rows_per_proc + std::min((unsigned long)rank, remainder);
    unsigned long row_end = m_row_start + rows_per_proc + (rank < (int)remainder ? 1 : 0);
    m_local_height = row_end - m_row_start;

    m_altitude.assign(m_local_height * m_dimensions, 0.0);

    RandomGenerator gen(seed, 0., (1UL << ln2_dim) * deviation);

    for (dim_t j = 0; j < m_dimensions; j += (1UL << ln2_dim)) {
        if (is_local(j)) {
            for (dim_t i = 0; i < m_dimensions; i += (1UL << ln2_dim)) {
                (*this)(i, j) = gen(i, j);
            }
        }
    }

    dim_t ldim = ln2_dim;
    unsigned long current_ss_grid = (1UL << ln2_dim);
    unsigned long current_nbSeeds = nbSeeds;

    while (ldim > 0) {
        for (unsigned long iB = 0; iB < current_nbSeeds; ++iB) {
            for (unsigned long jB = 0; jB < current_nbSeeds; ++jB) {
                compute_subgrid(ldim, iB, jB, deviation, seed);
            }
        }
        ldim -= 1;
        current_ss_grid /= 2;
        current_nbSeeds *= 2;
    }
}