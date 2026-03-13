#include "ant.hpp"
#include <algorithm>
#include <cstdint>
#include <omp.h>
#include "rand_generator.hpp"

double AntSwarm::m_eps = 0.;

AntSwarm::AntSwarm(std::size_t nb_ants) { resize(nb_ants); }

void AntSwarm::resize(std::size_t nb_ants) {
    m_x.assign(nb_ants, 0);
    m_y.assign(nb_ants, 0);
    m_is_loaded.assign(nb_ants, 0);
    m_seed.assign(nb_ants, 0);
    consumed_time.assign(nb_ants, 0.0);
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

void AntSwarm::add_ant(const AntData& data) {
    m_x.push_back(data.x);
    m_y.push_back(data.y);
    m_is_loaded.push_back(data.is_loaded);
    m_seed.push_back(data.seed);
    consumed_time.push_back(0.0); 
}

void AntSwarm::advance_one(pheronome& phen, const fractal_land& land,
                           const position_t& pos_food, const position_t& pos_nest, 
                           std::size_t& cpteur_food, int rank, int size)
{
    std::size_t nb_ants = m_x.size();
    if (nb_ants == 0) goto mpi_exchange;

    std::size_t active_loaded_count = 0;
    std::size_t active_unloaded_count = 0;

    for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(nb_ants); ++i) {
        consumed_time[i] = 0.0;
        if (m_is_loaded[i]) active_loaded_ids[active_loaded_count++] = i;
        else active_unloaded_ids[active_unloaded_count++] = i;
    }

    const int nb_threads = omp_get_max_threads();
    std::vector<std::vector<std::uint32_t>> next_loaded_local(nb_threads);
    std::vector<std::vector<std::uint32_t>> next_unloaded_local(nb_threads);
    std::vector<std::vector<position_t>> marks_local(nb_threads);
    std::vector<std::size_t> food_local(nb_threads, 0);

    for (int t = 0; t < nb_threads; ++t) {
        next_loaded_local[t].reserve(nb_ants / nb_threads + 64);
        next_unloaded_local[t].reserve(nb_ants / nb_threads + 64);
        marks_local[t].reserve(nb_ants / nb_threads + 64);
    }

    int north = (rank > 0) ? rank - 1 : MPI_PROC_NULL;
    int south = (rank < size - 1) ? rank + 1 : MPI_PROC_NULL;

    auto exchange = [&](std::vector<AntData>& send_buf, int dest, int source) {
        int send_count = send_buf.size();
        int recv_count = 0;
        MPI_Sendrecv(&send_count, 1, MPI_INT, dest, 0,
                     &recv_count, 1, MPI_INT, source, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        if (recv_count > 0 || send_count > 0) {
            std::vector<AntData> recv_buf(recv_count);
            MPI_Sendrecv(send_buf.data(), send_count * sizeof(AntData), MPI_BYTE, dest, 1,
                         recv_buf.data(), recv_count * sizeof(AntData), MPI_BYTE, source, 1,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            for (const auto& ant : recv_buf) add_ant(ant);
        }
    };

    auto traiter_liste = [&](std::vector<std::uint32_t>& liste_active, std::size_t nb_actives, int ind_pher,
                             std::vector<std::uint32_t>& liste_next_loaded, std::size_t& nb_next_loaded,
                             std::vector<std::uint32_t>& liste_next_unloaded, std::size_t& nb_next_unloaded)
    {
        for (int t = 0; t < nb_threads; ++t) {
            next_loaded_local[t].clear(); next_unloaded_local[t].clear();
            marks_local[t].clear(); food_local[t] = 0;
        }

        #pragma omp parallel if (nb_actives > 1024)
        {
            #pragma omp for schedule(guided, 64) nowait
            for (std::int64_t pos = 0; pos < static_cast<std::int64_t>(nb_actives); ++pos) {
                const std::uint32_t i = liste_active[static_cast<std::size_t>(pos)];
                const int tid = omp_get_thread_num();

                double choix = rand_choice_01(m_seed[i]);
                position_t old_pos = {m_x[i], m_y[i]};
                position_t new_pos = old_pos;

                const position_t pos_gauche{old_pos.x - 1, old_pos.y};
                const position_t pos_droite{old_pos.x + 1, old_pos.y};
                const position_t pos_haut{old_pos.x, old_pos.y - 1};
                const position_t pos_bas{old_pos.x, old_pos.y + 1};

                double max_phen = std::max({phen[pos_gauche][ind_pher],
                                            phen[pos_droite][ind_pher],
                                            phen[pos_haut][ind_pher],
                                            phen[pos_bas][ind_pher]});

                if ((choix > m_eps) || (max_phen <= 0.0)) {
                    do {
                        new_pos = old_pos;
                        int d = rand_dir_14(m_seed[i]);
                        if (d == 1) new_pos.x -= 1;
                        else if (d == 2) new_pos.y -= 1;
                        else if (d == 3) new_pos.x += 1;
                        else if (d == 4) new_pos.y += 1;
                    } while (phen[new_pos][ind_pher] == -1.0); 
                } else {
                    if (phen[pos_gauche][ind_pher] == max_phen)      new_pos = pos_gauche;
                    else if (phen[pos_droite][ind_pher] == max_phen) new_pos = pos_droite;
                    else if (phen[pos_haut][ind_pher] == max_phen)   new_pos = pos_haut;
                    else                                             new_pos = pos_bas;
                }

                m_x[i] = new_pos.x; 
                m_y[i] = new_pos.y;
                marks_local[static_cast<std::size_t>(tid)].push_back(new_pos);

                if (new_pos == pos_nest) {
                    if (m_is_loaded[i]) {
                        food_local[static_cast<std::size_t>(tid)]++;
                        m_is_loaded[i] = 0;
                    }
                }

                if (new_pos == pos_food) {
                    m_is_loaded[i] = 1;
                }

                consumed_time[i] += land(static_cast<unsigned long>(new_pos.x),
                                         static_cast<unsigned long>(new_pos.y));

                if (consumed_time[i] < 1.0) {
                    if (m_is_loaded[i]) 
                        next_loaded_local[static_cast<std::size_t>(tid)].push_back(i);
                    else 
                        next_unloaded_local[static_cast<std::size_t>(tid)].push_back(i);
                }
            }
        }

        // Fusion séquentielle: on évite la contention pendant le noyau chaud.
        for (int t = 0; t < nb_threads; ++t) {
            cpteur_food += food_local[t];
            for (auto& p : marks_local[t]) phen.mark_pheronome_xy(p.x, p.y);
            for (auto id : next_loaded_local[t]) liste_next_loaded[nb_next_loaded++] = id;
            for (auto id : next_unloaded_local[t]) liste_next_unloaded[nb_next_unloaded++] = id;
        }
    };

    // Ordonnancement par rondes: chaque fourmi active effectue un micro-deplacement par tour.
    while ((active_loaded_count + active_unloaded_count) > 0) {
        std::size_t n_l_c = 0, n_u_c = 0;
        traiter_liste(active_unloaded_ids, active_unloaded_count, 0, next_loaded_ids, n_l_c, next_unloaded_ids, n_u_c);
        traiter_liste(active_loaded_ids, active_loaded_count, 1, next_loaded_ids, n_l_c, next_unloaded_ids, n_u_c);
        std::swap(active_loaded_ids, next_loaded_ids); std::swap(active_unloaded_ids, next_unloaded_ids);
        active_loaded_count = n_l_c; active_unloaded_count = n_u_c;
    }
    std::vector<AntData> to_north, to_south;
    std::size_t write_idx = 0;

    mpi_exchange:
        int r_start = (int)land.row_start();
        int r_height = (int)land.local_height();

        for (std::size_t i = 0; i < m_x.size(); ++i) {
            if (m_y[i] < r_start) 
                to_north.push_back({m_x[i], m_y[i], m_is_loaded[i], m_seed[i]});
            else if (m_y[i] >= r_start + r_height)
                to_south.push_back({m_x[i], m_y[i], m_is_loaded[i], m_seed[i]});
            else {
                m_x[write_idx] = m_x[i]; m_y[write_idx] = m_y[i];
                m_is_loaded[write_idx] = m_is_loaded[i]; m_seed[write_idx] = m_seed[i];
                write_idx++;
            }
        }
        m_x.resize(write_idx); m_y.resize(write_idx); m_is_loaded.resize(write_idx); m_seed.resize(write_idx);
        consumed_time.resize(write_idx);

    
        int north = (rank > 0) ? rank - 1 : MPI_PROC_NULL;
        int south = (rank < size - 1) ? rank + 1 : MPI_PROC_NULL;

        exchange(to_north, north, south);
        exchange(to_south, south, north);

        unsigned long local_food = static_cast<unsigned long>(cpteur_food);
        unsigned long global_food = 0;
        MPI_Allreduce(&local_food, &global_food, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        cpteur_food = static_cast<std::size_t>(global_food);

        if (m_x.size() > active_loaded_ids.size()) {
            active_loaded_ids.resize(m_x.size());
            active_unloaded_ids.resize(m_x.size());
            next_loaded_ids.resize(m_x.size());
            next_unloaded_ids.resize(m_x.size());
        }
}
