#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <mpi.h>

#include "ant.hpp"
#include "fractal_land.hpp"
#include "pheronome.hpp"
#include "rand_generator.hpp"
#include "renderer.hpp"
#include "window.hpp"

struct Config {
    std::size_t grid_size = 513;
    std::size_t nb_ants = 5000;
    double alpha = 0.7;
    double beta = 0.999;
    double eps = 0.8;
    std::size_t max_iters = 2000;
    int repetitions = 5;
    std::size_t seed = 2026;
    bool gui = true;
};

struct ParametresFractal {
    fractal_land::dim_t log2_sous_grille = 9;
    unsigned long nb_graines = 1;
    std::size_t taille_effective = 513;
};

struct PartitionRows {
    std::size_t row_start = 0;
    std::size_t row_end = 0;
    std::size_t local_height = 0;
};

struct MesuresLocales {
    double t_compute_move_ms = 0.0;
    double t_compute_evap_ms = 0.0;
    double t_compute_update_ms = 0.0;
    double t_render_ms = 0.0;
    double t_compute_ms = 0.0;

    double t_comm_halo_ms = 0.0;
    double t_comm_migration_ms = 0.0;
    double t_comm_food_ms = 0.0;
    double t_comm_ms = 0.0;

    double t_total_ms = 0.0;
    std::size_t first_iteration_food = 0;
    std::size_t food_kpi = 0;
    std::uint64_t total_ants_global = 0;
    std::uint64_t checksum_global = 0;

    double ants_min = 0.0;
    double ants_avg = 0.0;
    double ants_max = 0.0;
    double ants_stddev = 0.0;
    double imbalance_ratio = 0.0;
};

struct ReferenceBaseline {
    bool disponible = false;
    std::string message;
    std::filesystem::path source;
    std::vector<std::size_t> food_kpi;
    std::vector<std::uint64_t> checksums;
};

namespace {

int to_int_count(std::size_t n, const char* what) {
    if (n > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error(std::string(what) + " exceeds MPI int count");
    }
    return static_cast<int>(n);
}

PartitionRows decouper_lignes(std::size_t total_rows, int rank, int size) {
    const std::size_t base = total_rows / static_cast<std::size_t>(size);
    const std::size_t rem = total_rows % static_cast<std::size_t>(size);
    PartitionRows part;
    part.row_start = static_cast<std::size_t>(rank) * base +
                     std::min(static_cast<std::size_t>(rank), rem);
    part.local_height = base + (static_cast<std::size_t>(rank) < rem ? 1ULL : 0ULL);
    part.row_end = part.row_start + part.local_height;
    return part;
}

int owner_rank_for_row(int y_global, std::size_t total_rows, int size) {
    for (int r = 0; r < size; ++r) {
        const auto part = decouper_lignes(total_rows, r, size);
        if (y_global >= static_cast<int>(part.row_start) &&
            y_global < static_cast<int>(part.row_end)) {
            return r;
        }
    }
    return size - 1;
}

void append_ants_from_bytes(const std::vector<unsigned char>& bytes, std::vector<AntData>& out) {
    if (bytes.empty()) {
        return;
    }
    if ((bytes.size() % sizeof(AntData)) != 0U) {
        throw std::runtime_error("Invalid AntData payload size");
    }
    const std::size_t n = bytes.size() / sizeof(AntData);
    const auto* ptr = reinterpret_cast<const AntData*>(bytes.data());
    out.insert(out.end(), ptr, ptr + n);
}

std::uint64_t bits_double(double valeur) {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &valeur, sizeof(double));
    return bits;
}

void melanger_checksum(std::uint64_t& checksum, std::uint64_t valeur) {
    checksum ^= valeur + 0x9e3779b97f4a7c15ULL + (checksum << 6U) + (checksum >> 2U);
}

double calculer_moyenne(const std::vector<double>& valeurs) {
    if (valeurs.empty()) {
        return 0.0;
    }
    const double somme = std::accumulate(valeurs.begin(), valeurs.end(), 0.0);
    return somme / static_cast<double>(valeurs.size());
}

double calculer_ecart_type(const std::vector<double>& valeurs, double moyenne) {
    if (valeurs.empty()) {
        return 0.0;
    }
    const double somme_carres =
        std::inner_product(valeurs.begin(), valeurs.end(), valeurs.begin(), 0.0);
    double variance = somme_carres / static_cast<double>(valeurs.size()) - moyenne * moyenne;
    if (variance < 0.0 && variance > -1e-12) {
        variance = 0.0;
    }
    return std::sqrt(std::max(0.0, variance));
}

void afficher_aide(const char* programme) {
    std::cout
        << "Usage: " << programme << " [options]\n"
        << "Options:\n"
        << "  --grid-size N   Grid size (default: 513)\n"
        << "  --ants N        Number of ants (default: 5000)\n"
        << "  --alpha X       Alpha in [0,1] (default: 0.7)\n"
        << "  --beta X        Beta in [0,1] (default: 0.999)\n"
        << "  --eps X         Exploration coefficient in [0,1] (default: 0.8)\n"
        << "  --iter N        Max iterations (default: 2000)\n"
        << "  --rep N         Repetitions (default: 5)\n"
        << "  --seed N        Seed (default: 2026)\n"
        << "  --gui 0|1       SDL display (default: 1)\n"
        << "  --help          Show this help\n";
}

ParametresFractal construire_parametres_fractal(std::size_t taille_demandee) {
    if (taille_demandee < 3) {
        throw std::invalid_argument("--grid-size must be >= 3");
    }
    std::size_t base = taille_demandee - 1;
    fractal_land::dim_t puissance = 0;
    while ((base % 2ULL) == 0ULL && base > 1ULL) {
        base /= 2ULL;
        ++puissance;
    }
    if (puissance >= std::numeric_limits<std::size_t>::digits) {
        throw std::invalid_argument("--grid-size too large");
    }
    const auto nb_graines = static_cast<unsigned long>(base);
    const auto taille_effective =
        static_cast<std::size_t>(nb_graines) * (std::size_t{1} << puissance) + 1ULL;
    return ParametresFractal{puissance, nb_graines, taille_effective};
}

Config parse_args(int argc, char* argv[]) {
    Config cfg;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        const auto next_value = [&](const std::string& option) -> const char* {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value for " + option);
            }
            return argv[++i];
        };

        if (arg == "--grid-size") {
            cfg.grid_size = std::stoull(next_value(arg));
        } else if (arg == "--ants") {
            cfg.nb_ants = std::stoull(next_value(arg));
        } else if (arg == "--alpha") {
            cfg.alpha = std::stod(next_value(arg));
        } else if (arg == "--beta") {
            cfg.beta = std::stod(next_value(arg));
        } else if (arg == "--eps") {
            cfg.eps = std::stod(next_value(arg));
        } else if (arg == "--iter") {
            cfg.max_iters = std::stoull(next_value(arg));
        } else if (arg == "--rep") {
            cfg.repetitions = std::stoi(next_value(arg));
        } else if (arg == "--seed") {
            cfg.seed = std::stoull(next_value(arg));
        } else if (arg == "--gui") {
            const int gui = std::stoi(next_value(arg));
            if (gui != 0 && gui != 1) {
                throw std::invalid_argument("--gui expects 0 or 1");
            }
            cfg.gui = (gui == 1);
        } else if (arg == "--help") {
            afficher_aide(argv[0]);
            std::exit(0);
        } else {
            throw std::invalid_argument("Unknown option: " + arg);
        }
    }

    if (cfg.nb_ants == 0) {
        throw std::invalid_argument("--ants must be > 0");
    }
    if (cfg.max_iters == 0) {
        throw std::invalid_argument("--iter must be > 0");
    }
    if (cfg.repetitions <= 0) {
        throw std::invalid_argument("--rep must be > 0");
    }
    if (cfg.alpha < 0.0 || cfg.alpha > 1.0 || cfg.beta < 0.0 || cfg.beta > 1.0 ||
        cfg.eps < 0.0 || cfg.eps > 1.0) {
        throw std::invalid_argument("--alpha, --beta and --eps must be in [0,1]");
    }
    return cfg;
}

void normaliser_terrain(fractal_land& land) {
    double local_max = std::numeric_limits<double>::lowest();
    double local_min = std::numeric_limits<double>::max();

    for (unsigned long y_global = land.row_start(); y_global < land.row_end(); ++y_global) {
        for (unsigned long x = 0; x < land.dimensions(); ++x) {
            local_max = std::max(local_max, land(x, y_global));
            local_min = std::min(local_min, land(x, y_global));
        }
    }

    double global_max = 0.0;
    double global_min = 0.0;
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    const double delta = global_max - global_min;
    const unsigned long y_begin = (land.row_start() > 0) ? (land.row_start() - 1) : land.row_start();
    const unsigned long y_end_inclusive =
        (land.row_end() < land.dimensions()) ? land.row_end() : (land.row_end() - 1);

    if (delta <= std::numeric_limits<double>::epsilon()) {
        for (unsigned long y_global = y_begin; y_global <= y_end_inclusive; ++y_global) {
            for (unsigned long x = 0; x < land.dimensions(); ++x) {
                land(x, y_global) = 0.0;
            }
        }
        return;
    }

    for (unsigned long y_global = y_begin; y_global <= y_end_inclusive; ++y_global) {
        for (unsigned long x = 0; x < land.dimensions(); ++x) {
            land(x, y_global) = (land(x, y_global) - global_min) / delta;
        }
    }
}

std::vector<AntData> distribuer_fourmis_initiales(const Config& cfg,
                                                  std::size_t dim,
                                                  int rank,
                                                  int size) {
    std::vector<int> send_counts_bytes;
    std::vector<int> send_displs_bytes;
    std::vector<unsigned char> send_buffer;

    if (rank == 0) {
        std::vector<std::vector<AntData>> buckets(static_cast<std::size_t>(size));
        std::size_t seed = cfg.seed;
        for (std::size_t id = 0; id < cfg.nb_ants; ++id) {
            AntData ant;
            ant.x = rand_int32(1, static_cast<std::int32_t>(dim - 2), seed);
            ant.y = rand_int32(1, static_cast<std::int32_t>(dim - 2), seed);
            ant.is_loaded = 0;
            ant.seed = static_cast<std::uint32_t>(seed);
            ant.consumed_time = 0.0;
            ant.id = id;
            const int owner = owner_rank_for_row(ant.y, dim, size);
            buckets[static_cast<std::size_t>(owner)].push_back(ant);
        }

        send_counts_bytes.assign(static_cast<std::size_t>(size), 0);
        send_displs_bytes.assign(static_cast<std::size_t>(size), 0);
        std::size_t total_bytes = 0;
        for (int r = 0; r < size; ++r) {
            const std::size_t bytes =
                buckets[static_cast<std::size_t>(r)].size() * sizeof(AntData);
            send_counts_bytes[static_cast<std::size_t>(r)] =
                to_int_count(bytes, "initial ant bytes per rank");
            total_bytes += bytes;
        }
        for (int r = 1; r < size; ++r) {
            send_displs_bytes[static_cast<std::size_t>(r)] =
                send_displs_bytes[static_cast<std::size_t>(r - 1)] +
                send_counts_bytes[static_cast<std::size_t>(r - 1)];
        }

        send_buffer.assign(total_bytes, 0U);
        for (int r = 0; r < size; ++r) {
            const auto& ants = buckets[static_cast<std::size_t>(r)];
            if (ants.empty()) {
                continue;
            }
            const std::size_t offset = static_cast<std::size_t>(send_displs_bytes[static_cast<std::size_t>(r)]);
            const std::size_t bytes = ants.size() * sizeof(AntData);
            std::memcpy(send_buffer.data() + offset, ants.data(), bytes);
        }
    }

    int recv_bytes = 0;
    MPI_Scatter(rank == 0 ? send_counts_bytes.data() : nullptr,
                1,
                MPI_INT,
                &recv_bytes,
                1,
                MPI_INT,
                0,
                MPI_COMM_WORLD);

    std::vector<unsigned char> recv_buffer(static_cast<std::size_t>(recv_bytes), 0U);
    MPI_Scatterv(rank == 0 ? send_buffer.data() : nullptr,
                 rank == 0 ? send_counts_bytes.data() : nullptr,
                 rank == 0 ? send_displs_bytes.data() : nullptr,
                 MPI_BYTE,
                 recv_bytes > 0 ? recv_buffer.data() : nullptr,
                 recv_bytes,
                 MPI_BYTE,
                 0,
                 MPI_COMM_WORLD);

    std::vector<AntData> local_ants;
    append_ants_from_bytes(recv_buffer, local_ants);
    return local_ants;
}

void calculer_stats_charge(const AntSwarm& ants,
                           int size,
                           std::uint64_t& total_ants_global,
                           double& ants_min,
                           double& ants_avg,
                           double& ants_max,
                           double& ants_stddev,
                           double& imbalance_ratio) {
    const double local_count = static_cast<double>(ants.size());
    const double local_sq = local_count * local_count;
    double sum_counts = 0.0;
    double sum_sq = 0.0;

    MPI_Allreduce(&local_count, &ants_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_count, &ants_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_count, &sum_counts, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_sq, &sum_sq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    unsigned long long local_u64 = static_cast<unsigned long long>(ants.size());
    unsigned long long total_u64 = 0ULL;
    MPI_Allreduce(&local_u64, &total_u64, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    total_ants_global = static_cast<std::uint64_t>(total_u64);

    ants_avg = sum_counts / static_cast<double>(size);
    double variance = sum_sq / static_cast<double>(size) - ants_avg * ants_avg;
    if (variance < 0.0 && variance > -1e-12) {
        variance = 0.0;
    }
    ants_stddev = std::sqrt(std::max(0.0, variance));
    imbalance_ratio = (ants_avg > 0.0) ? (ants_max / ants_avg) : 0.0;
}

std::uint64_t calculer_checksum_global(const pheronome& phen,
                                       const AntSwarm& ants,
                                       int rank,
                                       int size) {
    const std::size_t dim = static_cast<std::size_t>(phen.dimensions());

    std::vector<double> local_map;
    phen.pack_interior(local_map);
    const int local_map_count = to_int_count(local_map.size(), "local map count");

    std::vector<int> map_counts;
    std::vector<int> map_displs;
    if (rank == 0) {
        map_counts.assign(static_cast<std::size_t>(size), 0);
        map_displs.assign(static_cast<std::size_t>(size), 0);
    }
    MPI_Gather(&local_map_count,
               1,
               MPI_INT,
               rank == 0 ? map_counts.data() : nullptr,
               1,
               MPI_INT,
               0,
               MPI_COMM_WORLD);

    std::vector<double> global_map;
    if (rank == 0) {
        int total = 0;
        for (int r = 0; r < size; ++r) {
            if (r > 0) {
                map_displs[static_cast<std::size_t>(r)] =
                    map_displs[static_cast<std::size_t>(r - 1)] +
                    map_counts[static_cast<std::size_t>(r - 1)];
            }
            total += map_counts[static_cast<std::size_t>(r)];
        }
        global_map.assign(static_cast<std::size_t>(total), 0.0);
    }

    MPI_Gatherv(local_map_count > 0 ? local_map.data() : nullptr,
                local_map_count,
                MPI_DOUBLE,
                rank == 0 ? global_map.data() : nullptr,
                rank == 0 ? map_counts.data() : nullptr,
                rank == 0 ? map_displs.data() : nullptr,
                MPI_DOUBLE,
                0,
                MPI_COMM_WORLD);

    std::vector<AntData> local_ants;
    ants.export_data(local_ants);
    const std::size_t local_bytes_sz = local_ants.size() * sizeof(AntData);
    const int local_bytes = to_int_count(local_bytes_sz, "local ant bytes");

    std::vector<int> ant_counts_bytes;
    std::vector<int> ant_displs_bytes;
    if (rank == 0) {
        ant_counts_bytes.assign(static_cast<std::size_t>(size), 0);
        ant_displs_bytes.assign(static_cast<std::size_t>(size), 0);
    }
    MPI_Gather(&local_bytes,
               1,
               MPI_INT,
               rank == 0 ? ant_counts_bytes.data() : nullptr,
               1,
               MPI_INT,
               0,
               MPI_COMM_WORLD);

    std::vector<unsigned char> global_ant_bytes;
    if (rank == 0) {
        int total = 0;
        for (int r = 0; r < size; ++r) {
            if (r > 0) {
                ant_displs_bytes[static_cast<std::size_t>(r)] =
                    ant_displs_bytes[static_cast<std::size_t>(r - 1)] +
                    ant_counts_bytes[static_cast<std::size_t>(r - 1)];
            }
            total += ant_counts_bytes[static_cast<std::size_t>(r)];
        }
        global_ant_bytes.assign(static_cast<std::size_t>(total), 0U);
    }

    MPI_Gatherv(local_bytes > 0 ? reinterpret_cast<const unsigned char*>(local_ants.data()) : nullptr,
                local_bytes,
                MPI_BYTE,
                rank == 0 ? global_ant_bytes.data() : nullptr,
                rank == 0 ? ant_counts_bytes.data() : nullptr,
                rank == 0 ? ant_displs_bytes.data() : nullptr,
                MPI_BYTE,
                0,
                MPI_COMM_WORLD);

    if (rank != 0) {
        return 0ULL;
    }

    std::vector<double> global_map_0(dim * dim, 0.0);
    std::vector<double> global_map_1(dim * dim, 0.0);
    for (int r = 0; r < size; ++r) {
        const auto part = decouper_lignes(dim, r, size);
        std::size_t offset = static_cast<std::size_t>(map_displs[static_cast<std::size_t>(r)]);
        for (std::size_t y = part.row_start; y < part.row_end; ++y) {
            for (std::size_t x = 0; x < dim; ++x) {
                if (offset + 1ULL >= global_map.size()) {
                    throw std::runtime_error("Invalid gathered pheromone payload");
                }
                const std::size_t idx = x * dim + y;
                global_map_0[idx] = global_map[offset++];
                global_map_1[idx] = global_map[offset++];
            }
        }
    }

    std::vector<AntData> global_ants;
    append_ants_from_bytes(global_ant_bytes, global_ants);
    std::sort(global_ants.begin(), global_ants.end(), [](const AntData& a, const AntData& b) {
        return a.id < b.id;
    });

    std::uint64_t checksum = 1469598103934665603ULL;
    for (std::size_t x = 0; x < dim; ++x) {
        for (std::size_t y = 0; y < dim; ++y) {
            const std::size_t idx = x * dim + y;
            melanger_checksum(checksum, bits_double(global_map_0[idx]));
            melanger_checksum(checksum, bits_double(global_map_1[idx]));
        }
    }
    for (const auto& ant : global_ants) {
        melanger_checksum(checksum, static_cast<std::uint64_t>(static_cast<std::uint32_t>(ant.x)));
        melanger_checksum(checksum, static_cast<std::uint64_t>(static_cast<std::uint32_t>(ant.y)));
        melanger_checksum(checksum, ant.is_loaded != 0 ? 0xA5A5A5A5ULL : 0x5A5A5A5AULL);
    }
    return checksum;
}

ReferenceBaseline charger_reference_baseline_q2(const std::vector<std::filesystem::path>& chemins) {
    ReferenceBaseline ref;
    std::filesystem::path valid_path;
    for (const auto& p : chemins) {
        if (std::filesystem::exists(p)) {
            valid_path = p;
            break;
        }
    }
    if (valid_path.empty()) {
        ref.message = "Q2 reference CSV not found";
        return ref;
    }

    std::ifstream file(valid_path);
    if (!file.is_open()) {
        ref.message = "Cannot open Q2 reference CSV";
        return ref;
    }

    std::string line;
    if (!std::getline(file, line)) {
        ref.message = "Q2 reference CSV is empty";
        return ref;
    }

    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        std::stringstream ss(line);
        std::string rep, first_it, food, checksum;
        if (!std::getline(ss, rep, ',')) continue;
        if (rep == "MEAN" || rep == "STD_DEV") break;
        if (!std::getline(ss, first_it, ',')) continue;
        if (!std::getline(ss, food, ',')) continue;
        if (!std::getline(ss, checksum, ',')) continue;
        try {
            (void)std::stoull(rep);
            ref.food_kpi.push_back(static_cast<std::size_t>(std::stoull(food)));
            ref.checksums.push_back(static_cast<std::uint64_t>(std::stoull(checksum)));
        } catch (...) {
            continue;
        }
    }

    if (ref.food_kpi.empty() || ref.checksums.empty()) {
        ref.message = "No usable repetitions in Q2 reference CSV";
        return ref;
    }

    ref.disponible = true;
    ref.source = valid_path;
    ref.message = "Q2 reference loaded";
    return ref;
}

void afficher_verification_baseline_q2(int mpi_size,
                                       const std::vector<std::size_t>& food_kpi_q5,
                                       const std::vector<std::uint64_t>& checksums_q5) {
    if (mpi_size != 1) {
        std::cout << "Baseline check Q2: A_VERIFIER"
                  << " | reason=verification only defined for np=1\n";
        return;
    }

    const std::vector<std::filesystem::path> candidats = {
        "Q2/results/Q2_timings_breakdown.csv",
        "../Q2/results/Q2_timings_breakdown.csv",
        "results/Q2_timings_breakdown.csv",
        "../results/Q2_timings_breakdown.csv"};

    const auto ref = charger_reference_baseline_q2(candidats);
    if (!ref.disponible) {
        std::cout << "Baseline check Q2: A_VERIFIER | reason=" << ref.message << "\n";
        return;
    }

    const std::size_t n = std::min(food_kpi_q5.size(), std::min(ref.food_kpi.size(), ref.checksums.size()));
    if (n == 0) {
        std::cout << "Baseline check Q2: A_VERIFIER | reason=no common repetitions\n";
        return;
    }

    bool all_match = true;
    for (std::size_t i = 0; i < n; ++i) {
        if (food_kpi_q5[i] != ref.food_kpi[i] || checksums_q5[i] != ref.checksums[i]) {
            all_match = false;
            break;
        }
    }

    if (all_match) {
        std::cout << "Baseline check Q2: MATCH_VALUE_ONLY"
                  << " | source=" << ref.source.string()
                  << " | repetitions=" << n << "/" << food_kpi_q5.size()
                  << " | reason=reference_csv_without_full_param_metadata\n";
    } else {
        std::cout << "Baseline check Q2: MISMATCH"
                  << " | source=" << ref.source.string()
                  << " | repetitions=" << n << "/" << food_kpi_q5.size() << "\n";
        for (std::size_t i = 0; i < n; ++i) {
            if (food_kpi_q5[i] != ref.food_kpi[i] || checksums_q5[i] != ref.checksums[i]) {
                std::cout << "  Rep " << (i + 1)
                          << " | Q5(Food_KPI=" << food_kpi_q5[i]
                          << ", Checksum=" << checksums_q5[i]
                          << ") vs Q2(Food_KPI=" << ref.food_kpi[i]
                          << ", Checksum=" << ref.checksums[i] << ")\n";
            }
        }
    }
}

MesuresLocales executer_repetition(const Config& cfg,
                                   const ParametresFractal& params,
                                   bool activer_gui,
                                   int rank,
                                   int size,
                                   bool calculer_checksum) {
    MesuresLocales mesures;

    fractal_land land(params.log2_sous_grille,
                      params.nb_graines,
                      1.0,
                      static_cast<int>(cfg.seed & 0x7fffffffU),
                      rank,
                      size);
    normaliser_terrain(land);

    const int bound_min = 1;
    const int bound_max = static_cast<int>(land.dimensions()) - 2;
    position_t pos_nest{static_cast<int>(land.dimensions() / 2), static_cast<int>(land.dimensions() / 2)};
    position_t pos_food{
        std::max(bound_min, static_cast<int>((3 * land.dimensions()) / 4)),
        std::max(bound_min, static_cast<int>((3 * land.dimensions()) / 4))};
    pos_food.x = std::min(pos_food.x, bound_max);
    pos_food.y = std::min(pos_food.y, bound_max);
    if (pos_food == pos_nest) {
        pos_food.x = std::min(pos_food.x + 1, bound_max);
    }

    AntSwarm::set_exploration_coef(cfg.eps);
    AntSwarm ants;
    ants.set_from_data(distribuer_fourmis_initiales(cfg, land.dimensions(), rank, size));

    pheronome phen(land.dimensions(), pos_food, pos_nest, rank, size, cfg.alpha, cfg.beta);
    phen.exchange_halos(mesures.t_comm_halo_ms);

    const bool gui_local = activer_gui && (rank == 0) && (size == 1);
    std::unique_ptr<Window> win;
    std::unique_ptr<Renderer> renderer;
    if (gui_local) {
        win = std::make_unique<Window>("Q5 MPI domain decomposition",
                                       static_cast<int>(2 * land.dimensions() + 10),
                                       static_cast<int>(land.dimensions() + 266));
        renderer = std::make_unique<Renderer>(land, phen, pos_nest, pos_food, ants);
    }

    bool continue_loop = true;
    bool first_food_seen = false;
    std::size_t food_quantity = 0;
    std::size_t it = 0;

    const double t_total_start = MPI_Wtime();
    while (continue_loop && it < cfg.max_iters) {
        ++it;

        if (gui_local) {
            SDL_Event event;
            while (SDL_PollEvent(&event)) {
                if (event.type == SDL_QUIT) {
                    continue_loop = false;
                }
            }
        }

        bool global_continue = continue_loop;
        MPI_Bcast(&global_continue, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
        continue_loop = global_continue;
        if (!continue_loop) {
            break;
        }

        std::size_t food_local = 0;
        double comm_migration_ms = 0.0;

        double t0 = MPI_Wtime();
        ants.advance_one(phen, land, pos_food, pos_nest, food_local, comm_migration_ms);
        double t1 = MPI_Wtime();
        mesures.t_compute_move_ms += (t1 - t0) * 1000.0;
        mesures.t_comm_migration_ms += comm_migration_ms;

        unsigned long long global_food = 0ULL;
        unsigned long long local_food_u64 = static_cast<unsigned long long>(food_local);
        t0 = MPI_Wtime();
        MPI_Allreduce(&local_food_u64, &global_food, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        t1 = MPI_Wtime();
        mesures.t_comm_food_ms += (t1 - t0) * 1000.0;
        food_quantity += static_cast<std::size_t>(global_food);

        t0 = MPI_Wtime();
        phen.do_evaporation();
        t1 = MPI_Wtime();
        mesures.t_compute_evap_ms += (t1 - t0) * 1000.0;

        t0 = MPI_Wtime();
        phen.update();
        t1 = MPI_Wtime();
        mesures.t_compute_update_ms += (t1 - t0) * 1000.0;

        phen.exchange_halos(mesures.t_comm_halo_ms);

        if (gui_local) {
            t0 = MPI_Wtime();
            renderer->display(*win, food_quantity);
            win->blit();
            t1 = MPI_Wtime();
            mesures.t_render_ms += (t1 - t0) * 1000.0;
        }

        if (!first_food_seen && food_quantity > 0) {
            mesures.first_iteration_food = it;
            first_food_seen = true;
        }
    }
    const double t_total_end = MPI_Wtime();

    mesures.food_kpi = food_quantity;
    mesures.t_compute_ms = mesures.t_compute_move_ms + mesures.t_compute_evap_ms +
                           mesures.t_compute_update_ms + mesures.t_render_ms;
    mesures.t_comm_ms = mesures.t_comm_halo_ms + mesures.t_comm_migration_ms + mesures.t_comm_food_ms;
    mesures.t_total_ms = (t_total_end - t_total_start) * 1000.0;

    calculer_stats_charge(ants,
                          size,
                          mesures.total_ants_global,
                          mesures.ants_min,
                          mesures.ants_avg,
                          mesures.ants_max,
                          mesures.ants_stddev,
                          mesures.imbalance_ratio);

    if (calculer_checksum) {
        mesures.checksum_global = calculer_checksum_global(phen, ants, rank, size);
    }

    return mesures;
}

} // namespace

int main(int argc, char* argv[]) {
    int rank = 0;
    int size = 1;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    try {
        Config cfg = parse_args(argc, argv);
        const ParametresFractal params = construire_parametres_fractal(cfg.grid_size);

        if (size != 1) {
            cfg.gui = false;
        }
        if (rank != 0) {
            cfg.gui = false;
        }

        if (rank == 0) {
            if (params.taille_effective != cfg.grid_size) {
                std::cerr << "Warning: requested grid size " << cfg.grid_size
                          << " adjusted to " << params.taille_effective << "\n";
            }
            std::cout << "--- Q5 MPI domain decomposition (row-wise) ---\n";
            std::cout << "Grid=" << params.taille_effective
                      << " | Ants=" << cfg.nb_ants
                      << " | Iter=" << cfg.max_iters
                      << " | Rep=" << cfg.repetitions
                      << " | GUI=" << (cfg.gui ? "ON" : "OFF")
                      << " | MPI np=" << size << "\n";
            if (size != 1) {
                std::cout << "GUI forced OFF for np>1\n";
            }
            std::cout << "Warm-up...\n";
        }

        Config warmup_cfg = cfg;
        warmup_cfg.max_iters = std::min<std::size_t>(cfg.max_iters, 1ULL);
        (void)executer_repetition(warmup_cfg, params, false, rank, size, false);

        if (cfg.gui && rank == 0) {
            if (SDL_Init(SDL_INIT_VIDEO) != 0) {
                throw std::runtime_error("Failed to initialize SDL");
            }
        }

        std::vector<std::size_t> food_kpi;
        std::vector<std::uint64_t> checksums;
        std::vector<std::uint64_t> total_ants;
        std::vector<std::size_t> first_it;
        std::vector<double> t_move;
        std::vector<double> t_evap;
        std::vector<double> t_update;
        std::vector<double> t_compute;
        std::vector<double> t_comm_halo;
        std::vector<double> t_comm_migration;
        std::vector<double> t_comm_food;
        std::vector<double> t_comm;
        std::vector<double> t_total;
        std::vector<double> ants_min_v;
        std::vector<double> ants_avg_v;
        std::vector<double> ants_max_v;
        std::vector<double> ants_std_v;
        std::vector<double> imbalance_v;

        if (rank == 0) {
            const std::size_t reps = static_cast<std::size_t>(cfg.repetitions);
            food_kpi.reserve(reps);
            checksums.reserve(reps);
            total_ants.reserve(reps);
            first_it.reserve(reps);
            t_move.reserve(reps);
            t_evap.reserve(reps);
            t_update.reserve(reps);
            t_compute.reserve(reps);
            t_comm_halo.reserve(reps);
            t_comm_migration.reserve(reps);
            t_comm_food.reserve(reps);
            t_comm.reserve(reps);
            t_total.reserve(reps);
            ants_min_v.reserve(reps);
            ants_avg_v.reserve(reps);
            ants_max_v.reserve(reps);
            ants_std_v.reserve(reps);
            imbalance_v.reserve(reps);
        }

        for (int rep = 0; rep < cfg.repetitions; ++rep) {
            MesuresLocales local = executer_repetition(cfg, params, cfg.gui, rank, size, true);

            double g_move = 0.0;
            double g_evap = 0.0;
            double g_update = 0.0;
            double g_compute = 0.0;
            double g_comm_halo = 0.0;
            double g_comm_migration = 0.0;
            double g_comm_food = 0.0;
            double g_comm = 0.0;
            double g_total = 0.0;

            MPI_Reduce(&local.t_compute_move_ms, &g_move, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.t_compute_evap_ms, &g_evap, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.t_compute_update_ms, &g_update, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.t_compute_ms, &g_compute, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.t_comm_halo_ms, &g_comm_halo, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.t_comm_migration_ms, &g_comm_migration, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.t_comm_food_ms, &g_comm_food, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.t_comm_ms, &g_comm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.t_total_ms, &g_total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

            unsigned long long food_local = static_cast<unsigned long long>(local.food_kpi);
            unsigned long long food_global = 0ULL;
            MPI_Reduce(&food_local, &food_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

            unsigned long long ants_local = static_cast<unsigned long long>(local.total_ants_global);
            unsigned long long ants_global = 0ULL;
            MPI_Reduce(&ants_local, &ants_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

            double ants_min_global = 0.0;
            double ants_avg_global = 0.0;
            double ants_max_global = 0.0;
            double ants_std_global = 0.0;
            double imbalance_global = 0.0;
            MPI_Reduce(&local.ants_min, &ants_min_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.ants_avg, &ants_avg_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.ants_max, &ants_max_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.ants_stddev, &ants_std_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.imbalance_ratio, &imbalance_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

            unsigned long long first_local = static_cast<unsigned long long>(local.first_iteration_food);
            unsigned long long first_global = 0ULL;
            MPI_Reduce(&first_local, &first_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

            if (rank == 0) {
                food_kpi.push_back(static_cast<std::size_t>(food_global));
                checksums.push_back(local.checksum_global);
                total_ants.push_back(static_cast<std::uint64_t>(ants_global));
                first_it.push_back(static_cast<std::size_t>(first_global));
                t_move.push_back(g_move);
                t_evap.push_back(g_evap);
                t_update.push_back(g_update);
                t_compute.push_back(g_compute);
                t_comm_halo.push_back(g_comm_halo);
                t_comm_migration.push_back(g_comm_migration);
                t_comm_food.push_back(g_comm_food);
                t_comm.push_back(g_comm);
                t_total.push_back(g_total);
                ants_min_v.push_back(ants_min_global);
                ants_avg_v.push_back(ants_avg_global);
                ants_max_v.push_back(ants_max_global);
                ants_std_v.push_back(ants_std_global);
                imbalance_v.push_back(imbalance_global);

                std::cout << "Rep " << (rep + 1) << "/" << cfg.repetitions
                          << " | T_compute=" << std::fixed << std::setprecision(6) << g_compute << " ms"
                          << " | T_comm=" << g_comm << " ms"
                          << " | T_total=" << g_total << " ms"
                          << " | First_Iteration=" << first_global
                          << " | Food_KPI=" << food_global
                          << " | Checksum=" << local.checksum_global
                          << " | Total_Ants_Global=" << ants_global
                          << " | Ants(min/avg/max)=" << ants_min_global << "/" << ants_avg_global << "/" << ants_max_global
                          << " | Imbalance=" << imbalance_global
                          << "\n";
            }
        }

        if (cfg.gui && rank == 0) {
            SDL_Quit();
        }

        if (rank == 0) {
            std::filesystem::create_directories("results");
            std::ofstream csv("results/Q5_timings_breakdown.csv");
            if (!csv.is_open()) {
                throw std::runtime_error("Cannot write results/Q5_timings_breakdown.csv");
            }

            csv << "Repetition,MPI_Procs,First_Iteration,Food_KPI,Checksum,Total_Ants_Global,"
                   "T_compute_move_ms,T_compute_evap_ms,T_compute_update_ms,T_compute_ms,"
                   "T_comm_halo_ms,T_comm_migration_ms,T_comm_food_ms,T_comm_ms,T_total_ms,"
                   "Ants_Min,Ants_Avg,Ants_Max,Ants_Stddev,Imbalance_Ratio\n";

            for (std::size_t i = 0; i < t_total.size(); ++i) {
                csv << (i + 1) << ","
                    << size << ","
                    << first_it[i] << ","
                    << food_kpi[i] << ","
                    << checksums[i] << ","
                    << total_ants[i] << ","
                    << std::fixed << std::setprecision(6)
                    << t_move[i] << ","
                    << t_evap[i] << ","
                    << t_update[i] << ","
                    << t_compute[i] << ","
                    << t_comm_halo[i] << ","
                    << t_comm_migration[i] << ","
                    << t_comm_food[i] << ","
                    << t_comm[i] << ","
                    << t_total[i] << ","
                    << ants_min_v[i] << ","
                    << ants_avg_v[i] << ","
                    << ants_max_v[i] << ","
                    << ants_std_v[i] << ","
                    << imbalance_v[i] << "\n";
            }

            const double mean_move = calculer_moyenne(t_move);
            const double mean_evap = calculer_moyenne(t_evap);
            const double mean_update = calculer_moyenne(t_update);
            const double mean_compute = calculer_moyenne(t_compute);
            const double mean_comm_halo = calculer_moyenne(t_comm_halo);
            const double mean_comm_migration = calculer_moyenne(t_comm_migration);
            const double mean_comm_food = calculer_moyenne(t_comm_food);
            const double mean_comm = calculer_moyenne(t_comm);
            const double mean_total = calculer_moyenne(t_total);
            const double mean_ants_min = calculer_moyenne(ants_min_v);
            const double mean_ants_avg = calculer_moyenne(ants_avg_v);
            const double mean_ants_max = calculer_moyenne(ants_max_v);
            const double mean_ants_std = calculer_moyenne(ants_std_v);
            const double mean_imbalance = calculer_moyenne(imbalance_v);

            const double std_move = calculer_ecart_type(t_move, mean_move);
            const double std_evap = calculer_ecart_type(t_evap, mean_evap);
            const double std_update = calculer_ecart_type(t_update, mean_update);
            const double std_compute = calculer_ecart_type(t_compute, mean_compute);
            const double std_comm_halo = calculer_ecart_type(t_comm_halo, mean_comm_halo);
            const double std_comm_migration = calculer_ecart_type(t_comm_migration, mean_comm_migration);
            const double std_comm_food = calculer_ecart_type(t_comm_food, mean_comm_food);
            const double std_comm = calculer_ecart_type(t_comm, mean_comm);
            const double std_total = calculer_ecart_type(t_total, mean_total);
            const double std_ants_min = calculer_ecart_type(ants_min_v, mean_ants_min);
            const double std_ants_avg = calculer_ecart_type(ants_avg_v, mean_ants_avg);
            const double std_ants_max = calculer_ecart_type(ants_max_v, mean_ants_max);
            const double std_ants_std = calculer_ecart_type(ants_std_v, mean_ants_std);
            const double std_imbalance = calculer_ecart_type(imbalance_v, mean_imbalance);

            csv << "MEAN,"
                << size
                << ",,,,,"
                << std::fixed << std::setprecision(6)
                << mean_move << ","
                << mean_evap << ","
                << mean_update << ","
                << mean_compute << ","
                << mean_comm_halo << ","
                << mean_comm_migration << ","
                << mean_comm_food << ","
                << mean_comm << ","
                << mean_total << ","
                << mean_ants_min << ","
                << mean_ants_avg << ","
                << mean_ants_max << ","
                << mean_ants_std << ","
                << mean_imbalance << "\n";

            csv << "STD_DEV,"
                << size
                << ",,,,,"
                << std::fixed << std::setprecision(6)
                << std_move << ","
                << std_evap << ","
                << std_update << ","
                << std_compute << ","
                << std_comm_halo << ","
                << std_comm_migration << ","
                << std_comm_food << ","
                << std_comm << ","
                << std_total << ","
                << std_ants_min << ","
                << std_ants_avg << ","
                << std_ants_max << ","
                << std_ants_std << ","
                << std_imbalance << "\n";

            std::cout << "\n=== Q5 MPI domain decomposition results ===\n";
            std::cout << "Mean timings (ms): T_compute=" << mean_compute
                      << " | T_comm=" << mean_comm
                      << " | T_total=" << mean_total << "\n";
            std::cout << "Mean load balance: ants(min/avg/max)="
                      << mean_ants_min << "/" << mean_ants_avg << "/" << mean_ants_max
                      << " | imbalance=" << mean_imbalance << "\n";
            std::cout << "CSV generated: results/Q5_timings_breakdown.csv\n";

            afficher_verification_baseline_q2(size, food_kpi, checksums);
        }

        MPI_Finalize();
        return 0;
    } catch (const std::exception& e) {
        if (rank == 0) {
            std::cerr << "Error: " << e.what() << "\n";
            afficher_aide(argv[0]);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
}
