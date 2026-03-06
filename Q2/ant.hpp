#ifndef _ANT_HPP_
# define _ANT_HPP_
# include <cstddef>
# include <vector>
#include <cstdint>
# include "pheronome.hpp"
# include "fractal_land.hpp"
# include "basic_types.hpp"

/**
 * Stockage SoA des fourmis :
 * - positions x/y en tableaux contigus ;
 * - état chargé/non chargé dans un tableau dédié ;
 * - graine RNG propre à chaque fourmi.
 */
class AntSwarm
{
public:
    explicit AntSwarm(std::size_t nb_ants = 0);
    AntSwarm(const AntSwarm&) = default;
    AntSwarm(AntSwarm&&) = default;
    AntSwarm& operator=(const AntSwarm&) = default;
    AntSwarm& operator=(AntSwarm&&) = default;
    ~AntSwarm() = default;

    void resize(std::size_t nb_ants);
    std::size_t size() const { return m_x.size(); }

    static void set_exploration_coef(double eps) { m_eps = eps; }

    void initialiser_positions_aleatoires(const fractal_land& land, std::size_t& seed_global);

    void advance_one(pheronome& phen, const fractal_land& land,
                     const position_t& pos_food, const position_t& pos_nest, std::size_t& cpteur_food);

    int x_at(std::size_t ant_id) const { return m_x[ant_id]; }
    int y_at(std::size_t ant_id) const { return m_y[ant_id]; }
    bool is_loaded_at(std::size_t ant_id) const { return m_is_loaded[ant_id] != 0; }

private:
    static double m_eps;
    std::vector<int> m_x;
    std::vector<int> m_y;
    std::vector<std::uint8_t> m_is_loaded;
    std::vector<std::uint32_t> m_seed;
    std::vector<double> consumed_time;
    std::vector<std::uint32_t> active_loaded_ids;
    std::vector<std::uint32_t> active_unloaded_ids;
    std::vector<std::uint32_t> next_loaded_ids;
    std::vector<std::uint32_t> next_unloaded_ids;
};

#endif
