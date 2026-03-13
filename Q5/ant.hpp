#ifndef _ANT_HPP_
#define _ANT_HPP_

#include <cstddef>
#include <cstdint>
#include <vector>

#include "basic_types.hpp"
#include "fractal_land.hpp"
#include "pheronome.hpp"

struct AntData {
    int x = 0;
    int y = 0;
    int is_loaded = 0;
    std::uint32_t seed = 0U;
    double consumed_time = 0.0;
    std::uint64_t id = 0ULL;
};

class AntSwarm {
public:
    explicit AntSwarm(std::size_t nb_ants = 0);
    AntSwarm(const AntSwarm&) = default;
    AntSwarm(AntSwarm&&) = default;
    AntSwarm& operator=(const AntSwarm&) = default;
    AntSwarm& operator=(AntSwarm&&) = default;
    ~AntSwarm() = default;

    void resize(std::size_t nb_ants);
    std::size_t size() const { return m_ants.size(); }

    static void set_exploration_coef(double eps) { m_eps = eps; }

    void initialiser_positions_aleatoires(const fractal_land& land,
                                          std::size_t& seed_global,
                                          std::uint64_t id_offset = 0ULL);
    void set_from_data(const std::vector<AntData>& ants);
    void export_data(std::vector<AntData>& out) const;

    void advance_one(pheronome& phen,
                     const fractal_land& land,
                     const position_t& pos_food,
                     const position_t& pos_nest,
                     std::size_t& food_local,
                     double& t_comm_migration_ms);

    int x_at(std::size_t ant_id) const { return m_ants[ant_id].x; }
    int y_at(std::size_t ant_id) const { return m_ants[ant_id].y; }
    bool is_loaded_at(std::size_t ant_id) const { return m_ants[ant_id].is_loaded != 0; }
    double consumed_time_at(std::size_t ant_id) const { return m_ants[ant_id].consumed_time; }
    std::uint32_t seed_at(std::size_t ant_id) const { return m_ants[ant_id].seed; }
    std::uint64_t id_at(std::size_t ant_id) const { return m_ants[ant_id].id; }

private:
    static double m_eps;
    std::vector<AntData> m_ants;
};

#endif
