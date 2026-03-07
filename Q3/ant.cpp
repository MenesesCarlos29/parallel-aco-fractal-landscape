#include "ant.hpp"
#include <algorithm>
#include <cstdint>
#include <omp.h>
#include "rand_generator.hpp"

double AntSwarm::m_eps = 0.;

AntSwarm::AntSwarm(std::size_t nb_ants)
{
    resize(nb_ants);
}

void AntSwarm::resize(std::size_t nb_ants)
{
    m_x.assign(nb_ants, 0);
    m_y.assign(nb_ants, 0);
    m_is_loaded.assign(nb_ants, 0);
    m_seed.assign(nb_ants, 0);
    consumed_time.resize(nb_ants);
    active_loaded_ids.resize(nb_ants);
    active_unloaded_ids.resize(nb_ants);
    next_loaded_ids.resize(nb_ants);
    next_unloaded_ids.resize(nb_ants);
}
static inline std::uint32_t next_seed(std::uint32_t& seed)
{
    seed = static_cast<std::uint32_t>((1664525ull * seed + 1013904223ull) % 0xFFFFFFFFu);
    return seed;
}

static inline double rand_choice_01(std::uint32_t& seed)
{
    return static_cast<double>(next_seed(seed) % 2u);
}

static inline int rand_dir_14(std::uint32_t& seed)
{
    return 1 + static_cast<int>(next_seed(seed) % 4u);
}

void AntSwarm::initialiser_positions_aleatoires(const fractal_land& land, std::size_t& seed_global)
{
    auto gen_ant_pos = [&land, &seed_global]() {
        return rand_int32(1, static_cast<std::int32_t>(land.dimensions() - 2), seed_global);
    };

    for (std::size_t i = 0; i < m_x.size(); ++i) {
        m_x[i] = gen_ant_pos();
        m_y[i] = gen_ant_pos();
        m_is_loaded[i] = 0;
        m_seed[i] = seed_global;
    }
}

void AntSwarm::advance_one(pheronome& phen, const fractal_land& land,
                           const position_t& pos_food, const position_t& pos_nest, std::size_t& cpteur_food)
{
    // Debut d'une iteration globale: toutes les fourmis ont un budget de deplacement nul consomme.
    std::fill(consumed_time.begin(), consumed_time.end(), 0.0);

    // Separation des fourmis actives par etat pour eviter un branchement dans le noyau chaud.
    std::size_t active_loaded_count = 0;
    std::size_t active_unloaded_count = 0;

    for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(m_x.size()); ++i) {
        if (m_is_loaded[i]) {
            active_loaded_ids[active_loaded_count++] = i;
        } else {
            active_unloaded_ids[active_unloaded_count++] = i;
        }
    }

    // Buffers locaux OpenMP
    const int nb_threads = omp_get_max_threads();
    std::vector<std::vector<std::uint32_t>> next_loaded_local(static_cast<std::size_t>(nb_threads));
    std::vector<std::vector<std::uint32_t>> next_unloaded_local(static_cast<std::size_t>(nb_threads));
    std::vector<std::vector<position_t>> marks_local(static_cast<std::size_t>(nb_threads));
    std::vector<std::size_t> food_local(static_cast<std::size_t>(nb_threads), 0);
    const std::size_t reserve_par_thread =
        (m_x.size() / static_cast<std::size_t>(nb_threads)) + 64;
    for (int t = 0; t < nb_threads; ++t) {
        next_loaded_local[static_cast<std::size_t>(t)].reserve(reserve_par_thread);
        next_unloaded_local[static_cast<std::size_t>(t)].reserve(reserve_par_thread);
        marks_local[static_cast<std::size_t>(t)].reserve(reserve_par_thread);
    }

    auto traiter_liste = [&](std::vector<std::uint32_t>& liste_active,
                             std::size_t nb_actives,
                             int ind_pher,
                             std::vector<std::uint32_t>& liste_next_loaded,
                             std::size_t& nb_next_loaded,
                             std::vector<std::uint32_t>& liste_next_unloaded,
                             std::size_t& nb_next_unloaded)
    {
        for (int t = 0; t < nb_threads; ++t) {
            next_loaded_local[static_cast<std::size_t>(t)].clear();
            next_unloaded_local[static_cast<std::size_t>(t)].clear();
            marks_local[static_cast<std::size_t>(t)].clear();
            food_local[static_cast<std::size_t>(t)] = 0;
        }

        #pragma omp parallel if (nb_actives > 1024)
        {
            #pragma omp for schedule(guided,64) nowait
            for (std::int64_t pos = 0; pos < static_cast<std::int64_t>(nb_actives); ++pos) {
                const std::uint32_t i = liste_active[static_cast<std::size_t>(pos)];
                const int tid = omp_get_thread_num();

                double choix = rand_choice_01(m_seed[i]);
                position_t old_pos_ant{m_x[i], m_y[i]};
                position_t new_pos_ant = old_pos_ant;

                const position_t pos_gauche{new_pos_ant.x - 1, new_pos_ant.y};
                const position_t pos_droite{new_pos_ant.x + 1, new_pos_ant.y};
                const position_t pos_haut{new_pos_ant.x, new_pos_ant.y - 1};
                const position_t pos_bas{new_pos_ant.x, new_pos_ant.y + 1};

                double max_phen = std::max({phen[pos_gauche][ind_pher],
                                            phen[pos_droite][ind_pher],
                                            phen[pos_haut][ind_pher],
                                            phen[pos_bas][ind_pher]});

                if ((choix > m_eps) || (max_phen <= 0.0)) {
                    do {
                        new_pos_ant = old_pos_ant;
                        int d = rand_dir_14(m_seed[i]);

                        if (d == 1) new_pos_ant.x -= 1;
                        if (d == 2) new_pos_ant.y -= 1;
                        if (d == 3) new_pos_ant.x += 1;
                        if (d == 4) new_pos_ant.y += 1;

                    } while (phen[new_pos_ant][ind_pher] == -1.0);
                } else {
                    if (phen[pos_gauche][ind_pher] == max_phen)
                        new_pos_ant.x -= 1;
                    else if (phen[pos_droite][ind_pher] == max_phen)
                        new_pos_ant.x += 1;
                    else if (phen[pos_haut][ind_pher] == max_phen)
                        new_pos_ant.y -= 1;
                    else
                        new_pos_ant.y += 1;
                }

                m_x[i] = new_pos_ant.x;
                m_y[i] = new_pos_ant.y;
                marks_local[static_cast<std::size_t>(tid)].push_back(new_pos_ant);

                if (new_pos_ant == pos_nest) {
                    if (m_is_loaded[i] != 0) {
                        ++food_local[static_cast<std::size_t>(tid)];
                    }
                    m_is_loaded[i] = 0;
                }

                if (new_pos_ant == pos_food) {
                    m_is_loaded[i] = 1;
                }

                consumed_time[i] += land(static_cast<unsigned long>(new_pos_ant.x),
                                         static_cast<unsigned long>(new_pos_ant.y));

                if (consumed_time[i] < 1.0) {
                    if (m_is_loaded[i] != 0) {
                        next_loaded_local[static_cast<std::size_t>(tid)].push_back(i);
                    } else {
                        next_unloaded_local[static_cast<std::size_t>(tid)].push_back(i);
                    }
                }
            }
        }

        // Fusion séquentielle: on évite la contention pendant le noyau chaud.
        for (int t = 0; t < nb_threads; ++t) {
            cpteur_food += food_local[static_cast<std::size_t>(t)];

            for (const position_t& pos : marks_local[static_cast<std::size_t>(t)]) {
                phen.mark_pheronome_xy(static_cast<std::size_t>(pos.x),
                                       static_cast<std::size_t>(pos.y));
            }

            for (std::uint32_t i : next_loaded_local[static_cast<std::size_t>(t)]) {
                liste_next_loaded[nb_next_loaded++] = i;
            }

            for (std::uint32_t i : next_unloaded_local[static_cast<std::size_t>(t)]) {
                liste_next_unloaded[nb_next_unloaded++] = i;
            }
        }
    };

    // Ordonnancement par rondes: chaque fourmi active effectue un micro-deplacement par tour.
    while ((active_loaded_count + active_unloaded_count) != 0) {
        std::size_t next_loaded_count = 0;
        std::size_t next_unloaded_count = 0;

        traiter_liste(active_unloaded_ids,
                      active_unloaded_count,
                      0,
                      next_loaded_ids,
                      next_loaded_count,
                      next_unloaded_ids,
                      next_unloaded_count);

        traiter_liste(active_loaded_ids,
                      active_loaded_count,
                      1,
                      next_loaded_ids,
                      next_loaded_count,
                      next_unloaded_ids,
                      next_unloaded_count);

        std::swap(active_loaded_ids, next_loaded_ids);
        std::swap(active_unloaded_ids, next_unloaded_ids);

        active_loaded_count = next_loaded_count;
        active_unloaded_count = next_unloaded_count;
    }
}
