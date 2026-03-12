# include "fractal_land.hpp"
# include "rand_generator.hpp"
#include <vector>
#include <algorithm>

fractal_land::fractal_land(const dim_t& ln2_dim, unsigned long nbSeeds, double deviation, int seed, int rank, int size) :
    m_dimensions(0)
{
    // dim_ss_grid = 2^{ln2_dim}
    unsigned long dim_ss_grid = 1UL<<(ln2_dim);
    m_dimensions = nbSeeds*dim_ss_grid+1;

    unsigned long rows_per_proc = m_dimensions / size;
    unsigned long remainder = m_dimensions % size;

    m_row_start = rank * rows_per_proc + std::min((unsigned long)rank, remainder);
    unsigned long row_end = m_row_start + rows_per_proc + ((unsigned long)rank < remainder ? 1 : 0);
    m_local_height = row_end - m_row_start;

    std::vector<double> temp_global_altitude(m_dimensions * m_dimensions, 0.0);

    auto set_global = [&](unsigned long i, unsigned long j, double val) {
        temp_global_altitude[i + j * m_dimensions] = val;
    };

    auto get_global = [&](unsigned long i, unsigned long j) -> double& {
        return temp_global_altitude[i + j * m_dimensions];
    };

    // Seed the engine with an unsigned int
    RandomGenerator gen(seed, 0., dim_ss_grid*deviation);

    // Calcul des points initiaux :
    for ( dim_t i = 0; i < m_dimensions; i += dim_ss_grid )
        for ( dim_t j = 0; j < m_dimensions; j += dim_ss_grid )
            get_global(i, j) = gen(i, j);
    // Puis on itère pour calculer le paysage fractal :
    dim_t ldim = ln2_dim;
    unsigned long current_dim_ss_grid = dim_ss_grid;
    unsigned long current_nbSeeds = nbSeeds;

    while (ldim >= 1)
    {
        RandomGenerator gen_step(seed + ldim, -deviation, deviation);
        unsigned long mid = current_dim_ss_grid / 2;

        for (unsigned long iB = 0; iB < current_nbSeeds; ++iB) {
            for (unsigned long jB = 0; jB < current_nbSeeds; ++jB) {
                unsigned long iBeg = iB * current_dim_ss_grid;
                unsigned long jBeg = jB * current_dim_ss_grid;
                unsigned long iEnd = iBeg + current_dim_ss_grid;
                unsigned long jEnd = jBeg + current_dim_ss_grid;
                unsigned long iMid = iBeg + mid;
                unsigned long jMid = jBeg + mid;

                
                set_global(iMid, jBeg, 0.5 * (get_global(iBeg, jBeg) + get_global(iEnd, jBeg)) + mid * gen_step(iMid, jBeg));
                set_global(iBeg, jMid, 0.5 * (get_global(iBeg, jBeg) + get_global(iBeg, jEnd)) + mid * gen_step(iBeg, jMid));
                set_global(iMid, jEnd, 0.5 * (get_global(iBeg, jEnd) + get_global(iEnd, jEnd)) + mid * gen_step(iMid, jEnd));
                set_global(iEnd, jMid, 0.5 * (get_global(iEnd, jBeg) + get_global(iEnd, jEnd)) + mid * gen_step(iEnd, jMid));
                
                
                set_global(iMid, jMid, 0.25 * (get_global(iMid, jBeg) + get_global(iBeg, jMid) + 
                                              get_global(iMid, jEnd) + get_global(iEnd, jMid)) + mid * gen_step(iMid, jMid));
            }
        }
        ldim--;
        current_dim_ss_grid /= 2;
        current_nbSeeds *= 2;
    }

    m_altitude.resize(m_local_height * m_dimensions);
    for (unsigned long j = 0; j < m_local_height; ++j) {
        for (unsigned long i = 0; i < m_dimensions; ++i) {
            m_altitude[i + j * m_dimensions] = get_global   (i, j + m_row_start);
        }
    }
}