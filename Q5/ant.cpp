#include "ant.hpp"

#include <algorithm>
#include <array>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

#include <mpi.h>

#include "rand_generator.hpp"

double AntSwarm::m_eps = 0.0;

namespace {

static inline std::uint32_t next_seed(std::uint32_t& seed) {
    seed = static_cast<std::uint32_t>((1664525ULL * seed + 1013904223ULL) % 0xFFFFFFFFU);
    return seed;
}

static inline double rand_choice_01(std::uint32_t& seed) {
    return static_cast<double>(next_seed(seed) % 2U);
}

static inline int rand_dir_14(std::uint32_t& seed) {
    return 1 + static_cast<int>(next_seed(seed) % 4U);
}

int to_int_count(std::size_t n, const char* what) {
    if (n > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error(std::string(what) + " exceeds MPI int count");
    }
    return static_cast<int>(n);
}

void append_ant_from_bytes(const std::vector<unsigned char>& buffer, std::vector<AntData>& out) {
    if (buffer.empty()) {
        return;
    }
    if ((buffer.size() % sizeof(AntData)) != 0U) {
        throw std::runtime_error("Invalid ant migration payload size");
    }
    const std::size_t count = buffer.size() / sizeof(AntData);
    const auto* ptr = reinterpret_cast<const AntData*>(buffer.data());
    out.insert(out.end(), ptr, ptr + count);
}

} // namespace

AntSwarm::AntSwarm(std::size_t nb_ants) {
    resize(nb_ants);
}

void AntSwarm::resize(std::size_t nb_ants) {
    m_ants.assign(nb_ants, AntData{});
}

void AntSwarm::initialiser_positions_aleatoires(const fractal_land& land,
                                                std::size_t& seed_global,
                                                std::uint64_t id_offset) {
    for (std::size_t i = 0; i < m_ants.size(); ++i) {
        auto& ant = m_ants[i];
        ant.x = rand_int32(1, static_cast<std::int32_t>(land.dimensions() - 2), seed_global);
        ant.y = rand_int32(1, static_cast<std::int32_t>(land.dimensions() - 2), seed_global);
        ant.is_loaded = 0;
        ant.seed = static_cast<std::uint32_t>(seed_global);
        ant.consumed_time = 0.0;
        ant.id = id_offset + i;
    }
}

void AntSwarm::set_from_data(const std::vector<AntData>& ants) {
    m_ants = ants;
}

void AntSwarm::export_data(std::vector<AntData>& out) const {
    out = m_ants;
}

void AntSwarm::advance_one(pheronome& phen,
                           const fractal_land& land,
                           const position_t& pos_food,
                           const position_t& pos_nest,
                           std::size_t& food_local,
                           double& t_comm_migration_ms) {
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int up = (rank == 0) ? MPI_PROC_NULL : rank - 1;
    const int down = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;

    for (auto& ant : m_ants) {
        ant.consumed_time = 0.0;
    }

    std::vector<std::size_t> active_loaded;
    std::vector<std::size_t> active_unloaded;
    active_loaded.reserve(m_ants.size());
    active_unloaded.reserve(m_ants.size());
    for (std::size_t i = 0; i < m_ants.size(); ++i) {
        if (m_ants[i].is_loaded != 0) {
            active_loaded.push_back(i);
        } else {
            active_unloaded.push_back(i);
        }
    }

    const int y_start = static_cast<int>(phen.row_start());
    const int y_end_exclusive = static_cast<int>(phen.row_end());
    while (true) {

        int local_active = static_cast<int>(active_loaded.size() + active_unloaded.size());
        int global_active = 0;
        double t0_comm = MPI_Wtime();
        MPI_Allreduce(&local_active, &global_active, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        t_comm_migration_ms += (MPI_Wtime() - t0_comm) * 1000.0;
        if (global_active == 0) {
            break;
        }

        std::vector<AntData> send_up;
        std::vector<AntData> send_down;
        std::vector<std::size_t> next_loaded;
        std::vector<std::size_t> next_unloaded;
        std::vector<char> migrated(m_ants.size(), 0);

        send_up.reserve(active_loaded.size() + active_unloaded.size());
        send_down.reserve(active_loaded.size() + active_unloaded.size());
        next_loaded.reserve(active_loaded.size() + active_unloaded.size());
        next_unloaded.reserve(active_loaded.size() + active_unloaded.size());

        auto traiter_liste = [&](const std::vector<std::size_t>& active_ids, int ind_pher) {
            for (const std::size_t idx : active_ids) {
                auto& ant = m_ants[idx];

                const position_t old_pos{ant.x, ant.y};
                position_t new_pos = old_pos;
                const double choix = rand_choice_01(ant.seed);

                const position_t pos_left{new_pos.x - 1, new_pos.y};
                const position_t pos_right{new_pos.x + 1, new_pos.y};
                const position_t pos_up{new_pos.x, new_pos.y - 1};
                const position_t pos_down{new_pos.x, new_pos.y + 1};

                const double p_left = phen.pheromone_value(pos_left.x, pos_left.y, ind_pher);
                const double p_right = phen.pheromone_value(pos_right.x, pos_right.y, ind_pher);
                const double p_up = phen.pheromone_value(pos_up.x, pos_up.y, ind_pher);
                const double p_down = phen.pheromone_value(pos_down.x, pos_down.y, ind_pher);
                const double max_phen = std::max({p_left, p_right, p_up, p_down});

                if ((choix > m_eps) || (max_phen <= 0.0)) {
                    do {
                        new_pos = old_pos;
                        const int d = rand_dir_14(ant.seed);
                        if (d == 1) new_pos.x -= 1;
                        if (d == 2) new_pos.y -= 1;
                        if (d == 3) new_pos.x += 1;
                        if (d == 4) new_pos.y += 1;
                    } while (phen.pheromone_value(new_pos.x, new_pos.y, ind_pher) < 0.0);
                } else {
                    if (p_left == max_phen) {
                        new_pos.x -= 1;
                    } else if (p_right == max_phen) {
                        new_pos.x += 1;
                    } else if (p_up == max_phen) {
                        new_pos.y -= 1;
                    } else {
                        new_pos.y += 1;
                    }
                }

                ant.x = new_pos.x;
                ant.y = new_pos.y;

                if (new_pos == pos_nest) {
                    if (ant.is_loaded != 0) {
                        ++food_local;
                    }
                    ant.is_loaded = 0;
                }
                if (new_pos == pos_food) {
                    ant.is_loaded = 1;
                }

                ant.consumed_time += land.value_at(ant.x, ant.y);

                if (ant.y < y_start || ant.y >= y_end_exclusive) {
                    migrated[idx] = 1;
                    if (ant.y < y_start) {
                        send_up.push_back(ant);
                    } else {
                        send_down.push_back(ant);
                    }
                    continue;
                }

                phen.mark_pheronome_xy(static_cast<std::size_t>(ant.x), static_cast<std::size_t>(ant.y));

                if (ant.consumed_time < 1.0) {
                    if (ant.is_loaded != 0) {
                        next_loaded.push_back(idx);
                    } else {
                        next_unloaded.push_back(idx);
                    }
                }
            }
        };

        traiter_liste(active_unloaded, 0);
        traiter_liste(active_loaded, 1);

        std::vector<unsigned char> send_up_bytes;
        std::vector<unsigned char> send_down_bytes;
        send_up_bytes.resize(send_up.size() * sizeof(AntData));
        send_down_bytes.resize(send_down.size() * sizeof(AntData));
        if (!send_up.empty()) {
            std::memcpy(send_up_bytes.data(), send_up.data(), send_up_bytes.size());
        }
        if (!send_down.empty()) {
            std::memcpy(send_down_bytes.data(), send_down.data(), send_down_bytes.size());
        }

        int send_up_count = to_int_count(send_up_bytes.size(), "send_up bytes");
        int send_down_count = to_int_count(send_down_bytes.size(), "send_down bytes");
        int recv_up_count = 0;
        int recv_down_count = 0;

        t0_comm = MPI_Wtime();
        MPI_Sendrecv(&send_up_count, 1, MPI_INT, up, 10,
                     &recv_down_count, 1, MPI_INT, down, 10,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&send_down_count, 1, MPI_INT, down, 11,
                     &recv_up_count, 1, MPI_INT, up, 11,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        std::vector<unsigned char> recv_up_bytes(static_cast<std::size_t>(recv_up_count), 0U);
        std::vector<unsigned char> recv_down_bytes(static_cast<std::size_t>(recv_down_count), 0U);

        MPI_Sendrecv(send_up_count > 0 ? send_up_bytes.data() : nullptr,
                     send_up_count,
                     MPI_BYTE,
                     up,
                     20,
                     recv_down_count > 0 ? recv_down_bytes.data() : nullptr,
                     recv_down_count,
                     MPI_BYTE,
                     down,
                     20,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

        MPI_Sendrecv(send_down_count > 0 ? send_down_bytes.data() : nullptr,
                     send_down_count,
                     MPI_BYTE,
                     down,
                     21,
                     recv_up_count > 0 ? recv_up_bytes.data() : nullptr,
                     recv_up_count,
                     MPI_BYTE,
                     up,
                     21,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        t_comm_migration_ms += (MPI_Wtime() - t0_comm) * 1000.0;

        std::vector<int> remap(m_ants.size(), -1);
        std::vector<AntData> survivors;
        survivors.reserve(m_ants.size());
        for (std::size_t i = 0; i < m_ants.size(); ++i) {
            if (migrated[i] == 0) {
                remap[i] = static_cast<int>(survivors.size());
                survivors.push_back(m_ants[i]);
            }
        }
        m_ants.swap(survivors);

        std::vector<std::size_t> remapped_loaded;
        std::vector<std::size_t> remapped_unloaded;
        remapped_loaded.reserve(next_loaded.size());
        remapped_unloaded.reserve(next_unloaded.size());

        for (const std::size_t idx : next_loaded) {
            if (remap[idx] >= 0) {
                remapped_loaded.push_back(static_cast<std::size_t>(remap[idx]));
            }
        }
        for (const std::size_t idx : next_unloaded) {
            if (remap[idx] >= 0) {
                remapped_unloaded.push_back(static_cast<std::size_t>(remap[idx]));
            }
        }

        std::vector<AntData> received;
        received.reserve(static_cast<std::size_t>(recv_up_count + recv_down_count) / sizeof(AntData));
        append_ant_from_bytes(recv_up_bytes, received);
        append_ant_from_bytes(recv_down_bytes, received);

        for (const auto& ant : received) {
            const std::size_t new_index = m_ants.size();
            m_ants.push_back(ant);
            if (ant.y >= y_start && ant.y < y_end_exclusive) {
                phen.mark_pheronome_xy(static_cast<std::size_t>(ant.x), static_cast<std::size_t>(ant.y));
            }
            if (ant.consumed_time < 1.0) {
                if (ant.is_loaded != 0) {
                    remapped_loaded.push_back(new_index);
                } else {
                    remapped_unloaded.push_back(new_index);
                }
            }
        }

        active_loaded.swap(remapped_loaded);
        active_unloaded.swap(remapped_unloaded);
    }
}
