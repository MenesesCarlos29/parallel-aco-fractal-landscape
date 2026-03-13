#include "ant.hpp"
#include <algorithm>
#include <cstdint>
#include <mpi.h>
#include <omp.h>
#include "rand_generator.hpp"

double AntSwarm::m_eps = 0.;

AntSwarm::AntSwarm(std::size_t nb_ants)
{
    if (nb_ants > 0) {
            resize(nb_ants);
        }
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
    auto gen_ant_pos_x = [&land, &seed_global]() {
        return rand_int32(1, static_cast<std::int32_t>(land.dimensions() - 2), seed_global);
    };

    auto gen_ant_pos_y = [&land, &seed_global]() {
        const auto row_start = static_cast<std::int32_t>(land.row_start());
        const auto local_height = static_cast<std::int32_t>(land.local_height());
        const auto top = row_start;
        const auto bottom = row_start + local_height - 1;
        return rand_int32(top, bottom, seed_global);
    };

    for (std::size_t i = 0; i < m_x.size(); ++i) {
        m_x[i] = gen_ant_pos_x();
        m_y[i] = gen_ant_pos_y();
        m_is_loaded[i] = 0;
        m_seed[i] = seed_global;
    }
}

void AntSwarm::advance_one(pheronome& phen, const fractal_land& land,
                           const position_t& pos_food, const position_t& pos_nest, std::size_t& cpteur_food)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    std::vector<AntData> send_up;
    std::vector<AntData> send_down;
    std::vector<AntData> recv_up;
    std::vector<AntData> recv_down;
    std::vector<char> is_sent(m_x.size(), 0);

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

    active_loaded_ids.resize(active_loaded_count);
    active_unloaded_ids.resize(active_unloaded_count);

    // Buffers locaux OpenMP
    next_loaded_ids.clear();
    next_unloaded_ids.clear();

    next_loaded_ids.reserve(m_x.size() * 2 + 16);
    next_unloaded_ids.reserve(m_x.size() * 2 + 16);
    const int nb_threads = omp_get_max_threads();
    std::vector<std::vector<std::uint32_t>> next_loaded_local(static_cast<std::size_t>(nb_threads));
    std::vector<std::vector<std::uint32_t>> next_unloaded_local(static_cast<std::size_t>(nb_threads));
    std::vector<std::vector<position_t>> marks_local(static_cast<std::size_t>(nb_threads));
    std::vector<std::vector<AntData>> send_up_local(static_cast<std::size_t>(nb_threads));
    std::vector<std::vector<AntData>> send_down_local(static_cast<std::size_t>(nb_threads));
    std::vector<std::size_t> food_local(static_cast<std::size_t>(nb_threads), 0);
    const std::size_t reserve_par_thread =
        (m_x.size() / static_cast<std::size_t>(nb_threads)) + 64;
    for (int t = 0; t < nb_threads; ++t) {
        next_loaded_local[static_cast<std::size_t>(t)].reserve(reserve_par_thread);
        next_unloaded_local[static_cast<std::size_t>(t)].reserve(reserve_par_thread);
        marks_local[static_cast<std::size_t>(t)].reserve(reserve_par_thread);
        send_up_local[static_cast<std::size_t>(t)].reserve(reserve_par_thread);
        send_down_local[static_cast<std::size_t>(t)].reserve(reserve_par_thread);
    }

    auto traiter_liste = [&](std::vector<std::uint32_t>& liste_active,
                             std::size_t nb_actives,
                             int ind_pher,
                             std::vector<std::uint32_t>& liste_next_loaded,
                             std::vector<std::uint32_t>& liste_next_unloaded)
    {
        for (int t = 0; t < nb_threads; ++t) {
            send_up_local[t].clear();
            send_down_local[t].clear();
            next_loaded_local[static_cast<std::size_t>(t)].clear();
            next_unloaded_local[static_cast<std::size_t>(t)].clear();
            marks_local[static_cast<std::size_t>(t)].clear();
            food_local[static_cast<std::size_t>(t)] = 0;
        }

        #pragma omp parallel if (nb_actives > 1024)
        {
            const int tid = omp_get_thread_num();
            #pragma omp for schedule(guided,64) nowait
            for (std::int64_t pos = 0; pos < static_cast<std::int64_t>(nb_actives); ++pos) {
                const std::uint32_t i = liste_active[static_cast<std::size_t>(pos)];

                double choix = rand_choice_01(m_seed[i]);
                position_t old_pos_ant{m_x[i], m_y[i]};

                if (old_pos_ant.y < 0 || static_cast<std::size_t>(old_pos_ant.y) >= land.dimensions() ||
                    old_pos_ant.x < 0 || static_cast<std::size_t>(old_pos_ant.x) >= land.dimensions()) {
                    // pos invalide, remonter / descendre vers le bon rang
                    AntData d{old_pos_ant.x, old_pos_ant.y, static_cast<int>(m_is_loaded[i]), m_seed[i]};
                    is_sent[i] = 1;
                    if (old_pos_ant.y < static_cast<std::int32_t>(phen.row_start())) {
                        send_up_local[static_cast<std::size_t>(tid)].push_back(d);
                    } else {
                        send_down_local[static_cast<std::size_t>(tid)].push_back(d);
                    }
                    continue;
                }

                if (!phen.is_local(static_cast<std::size_t>(old_pos_ant.y))) {
                    AntData d{old_pos_ant.x, old_pos_ant.y, static_cast<int>(m_is_loaded[i]), m_seed[i]};
                    is_sent[i] = 1;
                    if (old_pos_ant.y < static_cast<std::int32_t>(phen.row_start())) {
                        send_up_local[static_cast<std::size_t>(tid)].push_back(d);
                    } else {
                        send_down_local[static_cast<std::size_t>(tid)].push_back(d);
                    }
                    continue;
                }

                position_t new_pos_ant = old_pos_ant;

                const position_t pos_gauche{new_pos_ant.x - 1, new_pos_ant.y};
                const position_t pos_droite{new_pos_ant.x + 1, new_pos_ant.y};
                const position_t pos_haut{new_pos_ant.x, new_pos_ant.y - 1};
                const position_t pos_bas{new_pos_ant.x, new_pos_ant.y + 1};

                auto get_pher = [&](const position_t& p) {
                    const int dim = static_cast<int>(land.dimensions());
                    if (p.x < 0 || p.x >= dim || p.y < 0 || p.y >= dim) return -1.0;

                    const long row_start = static_cast<long>(phen.row_start());
                    const long row_end = static_cast<long>(phen.row_start() + phen.local_height() - 1);
                    const long y = static_cast<long>(p.y);
                    if (y < row_start - 1 || y > row_end + 1) return -1.0;

                    return phen(static_cast<std::size_t>(p.x), static_cast<std::size_t>(p.y))[ind_pher];
                };

                double max_phen = std::max({
                    get_pher(pos_gauche),
                    get_pher(pos_droite),
                    get_pher(pos_haut),
                    get_pher(pos_bas)
                });

                if ((choix > m_eps) || (max_phen <= 0.0)) {
                    do {
                        new_pos_ant = old_pos_ant;
                        int d = rand_dir_14(m_seed[i]);

                        if (d == 1) new_pos_ant.x -= 1;
                        if (d == 2) new_pos_ant.y -= 1;
                        if (d == 3) new_pos_ant.x += 1;
                        if (d == 4) new_pos_ant.y += 1;

                    } while (get_pher(new_pos_ant) == -1.0);
                } else {
                    if (get_pher(pos_gauche) == max_phen)
                        new_pos_ant.x -= 1;
                    else if (get_pher(pos_droite) == max_phen)
                        new_pos_ant.x += 1;
                    else if (get_pher(pos_haut) == max_phen)
                        new_pos_ant.y -= 1;
                    else
                        new_pos_ant.y += 1;
                }

                m_x[i] = new_pos_ant.x;
                m_y[i] = new_pos_ant.y;

                if (!phen.is_local(static_cast<std::size_t>(new_pos_ant.y))) {
                    AntData d{new_pos_ant.x, new_pos_ant.y, static_cast<int>(m_is_loaded[i]), m_seed[i]};
                    is_sent[i] = 1;
                    if (new_pos_ant.y < old_pos_ant.y) {
                        send_up_local[static_cast<std::size_t>(tid)].push_back(d);
                    } else {
                        send_down_local[static_cast<std::size_t>(tid)].push_back(d);
                    }
                    continue;
                }

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

            for (const AntData& d : send_up_local[static_cast<std::size_t>(t)]) {
                send_up.push_back(d);
            }
            for (const AntData& d : send_down_local[static_cast<std::size_t>(t)]) {
                send_down.push_back(d);
            }

            for (const position_t& pos : marks_local[static_cast<std::size_t>(t)]) {
                phen.mark_pheronome_xy(static_cast<std::size_t>(pos.x),
                                       static_cast<std::size_t>(pos.y));
            }

            for (std::uint32_t i : next_loaded_local[static_cast<std::size_t>(t)]) {
                liste_next_loaded.push_back(i);
            }

            for (std::uint32_t i : next_unloaded_local[static_cast<std::size_t>(t)]) {
                liste_next_unloaded.push_back(i);
            }
        }
    };

    MPI_Datatype MPI_AntData;
    {
        AntData dummy;
        MPI_Aint base;
        MPI_Get_address(&dummy, &base);
        MPI_Aint displs[4];
        MPI_Get_address(&dummy.x, &displs[0]);
        MPI_Get_address(&dummy.y, &displs[1]);
        MPI_Get_address(&dummy.is_loaded, &displs[2]);
        MPI_Get_address(&dummy.seed, &displs[3]);
        for (int i = 0; i < 4; ++i) displs[i] -= base;
        int blocklens[4] = {1, 1, 1, 1};
        MPI_Datatype types[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_UINT32_T};
        MPI_Type_create_struct(4, blocklens, displs, types, &MPI_AntData);
        MPI_Type_commit(&MPI_AntData);
    }

    int up = (rank == 0) ? MPI_PROC_NULL : rank - 1;
    int down = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;

    std::size_t move_round = 0;
    const std::size_t max_move_rounds = std::max<std::size_t>(static_cast<std::size_t>(land.dimensions()) * 3, 1024);

    while ((active_loaded_count + active_unloaded_count) != 0 && move_round < max_move_rounds) {
        ++move_round;

        send_up.clear();
        send_down.clear();

        next_loaded_ids.clear();
        next_unloaded_ids.clear();

        traiter_liste(active_unloaded_ids,
                      active_unloaded_count,
                      0,
                      next_loaded_ids,
                      next_unloaded_ids);

        traiter_liste(active_loaded_ids,
                      active_loaded_count,
                      1,
                      next_loaded_ids,
                      next_unloaded_ids);

        int send_up_count = static_cast<int>(send_up.size());
        int send_down_count = static_cast<int>(send_down.size());
        int recv_up_count = 0;
        int recv_down_count = 0;

        MPI_Sendrecv(&send_up_count, 1, MPI_INT, up, 10,
                     &recv_down_count, 1, MPI_INT, down, 10,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&send_down_count, 1, MPI_INT, down, 11,
                     &recv_up_count, 1, MPI_INT, up, 11,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        recv_up.resize(static_cast<std::size_t>(recv_up_count));
        recv_down.resize(static_cast<std::size_t>(recv_down_count));

        MPI_Sendrecv(send_up.data(), send_up_count, MPI_AntData, up, 20,
                     recv_down.data(), recv_down_count, MPI_AntData, down, 20,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(send_down.data(), send_down_count, MPI_AntData, down, 21,
                     recv_up.data(), recv_up_count, MPI_AntData, up, 21,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        auto add_received = [&](const AntData &ant) {
            m_x.push_back(ant.x);
            m_y.push_back(ant.y);
            m_is_loaded.push_back(static_cast<std::uint8_t>(ant.is_loaded));
            m_seed.push_back(ant.seed);
            consumed_time.push_back(0.0);
            is_sent.push_back(0);
            std::uint32_t new_index = static_cast<std::uint32_t>(m_x.size() - 1);
            if (ant.is_loaded != 0) {
                next_loaded_ids.push_back(new_index);
            } else {
                next_unloaded_ids.push_back(new_index);
            }
        };

        for (const AntData &ant : recv_up) add_received(ant);
        for (const AntData &ant : recv_down) add_received(ant);

        std::swap(active_loaded_ids, next_loaded_ids);
        std::swap(active_unloaded_ids, next_unloaded_ids);

        active_loaded_count = active_loaded_ids.size();
        active_unloaded_count = active_unloaded_ids.size();
    }

    if (move_round >= max_move_rounds) {
        // Sécurité: si trop de rondes à cause de terrain 0, on force la fin du mouvement
        std::size_t stuck = active_loaded_count + active_unloaded_count;
        if (stuck > 0) {
            std::cerr << "avertissement: max_move_rounds atteint (" << move_round << ")";
            std::cerr << " | restants=" << stuck << "\n";
        }
    }

    std::size_t compact_write = 0;
    for (std::size_t read = 0, total = m_x.size(); read < total; ++read) {
        if (is_sent[read]) continue;
        if (compact_write != read) {
            m_x[compact_write] = m_x[read];
            m_y[compact_write] = m_y[read];
            m_is_loaded[compact_write] = m_is_loaded[read];
            m_seed[compact_write] = m_seed[read];
            consumed_time[compact_write] = consumed_time[read];
        }
        ++compact_write;
    }

    m_x.resize(compact_write);
    m_y.resize(compact_write);
    m_is_loaded.resize(compact_write);
    m_seed.resize(compact_write);
    consumed_time.resize(compact_write);
    is_sent.resize(compact_write);

    std::size_t global_food = 0;
    MPI_Allreduce(&cpteur_food, &global_food, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    cpteur_food = global_food;

    MPI_Type_free(&MPI_AntData);
}
