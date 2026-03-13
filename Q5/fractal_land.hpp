#ifndef _FRACTAL_LAND_HPP_
#define _FRACTAL_LAND_HPP_

#include <cstddef>
#include <vector>

class fractal_land {
public:
    using container = std::vector<double>;
    using dim_t = unsigned long;

    fractal_land(const dim_t& log_size,
                 unsigned long nbSeeds,
                 double deviation,
                 int seed,
                 int rank,
                 int size);
    fractal_land(const dim_t& log_size,
                 unsigned long nbSeeds,
                 double deviation,
                 int seed = 0);

    fractal_land(const fractal_land&) = delete;
    fractal_land(fractal_land&&) = default;
    ~fractal_land() = default;

    double operator()(unsigned long x_global, unsigned long y_global) const;
    double& operator()(unsigned long x_global, unsigned long y_global);

    double value_at(int x_global, int y_global) const;

    dim_t dimensions() const { return m_dimensions; }
    dim_t local_height() const { return m_local_height; }
    dim_t row_start() const { return m_row_start; }
    dim_t row_end() const { return m_row_start + m_local_height; }

private:
    std::size_t local_index(unsigned long x_global, unsigned long y_global) const;
    void fill_from_global(const std::vector<double>& global_altitude);

    dim_t m_dimensions = 0;
    dim_t m_row_start = 0;
    dim_t m_local_height = 0;
    container m_altitude;
};

#endif
