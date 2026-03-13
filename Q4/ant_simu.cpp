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

struct PartitionFourmis {
    std::size_t debut = 0;
    std::size_t fin = 0;
    std::size_t nb_ants = 0;
};

struct MesuresLocales {
    double t_compute_move_ms = 0.0;
    double t_compute_evap_ms = 0.0;
    double t_compute_update_ms = 0.0;
    double t_render_ms = 0.0;
    double t_compute_ms = 0.0;
    double t_comm_food_ms = 0.0;
    double t_comm_pher_ms = 0.0;
    double t_comm_ms = 0.0;
    double t_total_ms = 0.0;
    double ecart_sanity_ms = 0.0;
    std::size_t premiere_iteration_nourriture = 0;
    std::size_t food_kpi = 0;
    std::uint64_t checksum_global = 0;
};

struct ReferenceBaseline {
    bool disponible = false;
    std::string message;
    std::filesystem::path chemin_source;
    std::vector<std::size_t> food_kpi;
    std::vector<std::uint64_t> checksums;
};

void afficher_aide(const char* programme) {
    std::cout
        << "Usage: " << programme << " [options]\n"
        << "Options disponibles:\n"
        << "  --grid-size N   Taille de la grille (defaut: 513)\n"
        << "  --ants N        Nombre de fourmis (defaut: 5000)\n"
        << "  --alpha X       Coefficient alpha dans [0,1] (defaut: 0.7)\n"
        << "  --beta X        Coefficient beta dans [0,1] (defaut: 0.999)\n"
        << "  --eps X         Coefficient d'exploration dans [0,1] (defaut: 0.8)\n"
        << "  --iter N        Nombre d'iterations (defaut: 2000)\n"
        << "  --rep N         Nombre de repetitions (defaut: 5)\n"
        << "  --seed N        Graine initiale (defaut: 2026)\n"
        << "  --gui 0|1       Active/desactive SDL (defaut: 1)\n"
        << "  --help          Affiche cette aide\n";
}

ParametresFractal construire_parametres_fractal(std::size_t taille_demandee) {
    if (taille_demandee < 3) {
        throw std::invalid_argument("--grid-size doit etre >= 3.");
    }

    std::size_t base = taille_demandee - 1;
    fractal_land::dim_t puissance = 0;
    while ((base % 2ULL) == 0ULL && base > 1ULL) {
        base /= 2ULL;
        ++puissance;
    }

    if (puissance >= std::numeric_limits<std::size_t>::digits) {
        throw std::invalid_argument("--grid-size trop grand pour les calculs internes.");
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
        const auto valeur_suivante = [&](const std::string& option) -> const char* {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Valeur manquante pour " + option + ".");
            }
            return argv[++i];
        };

        if (arg == "--grid-size") {
            cfg.grid_size = std::stoull(valeur_suivante(arg));
        } else if (arg == "--ants") {
            cfg.nb_ants = std::stoull(valeur_suivante(arg));
        } else if (arg == "--alpha") {
            cfg.alpha = std::stod(valeur_suivante(arg));
        } else if (arg == "--beta") {
            cfg.beta = std::stod(valeur_suivante(arg));
        } else if (arg == "--eps") {
            cfg.eps = std::stod(valeur_suivante(arg));
        } else if (arg == "--iter") {
            cfg.max_iters = std::stoull(valeur_suivante(arg));
        } else if (arg == "--rep") {
            cfg.repetitions = std::stoi(valeur_suivante(arg));
        } else if (arg == "--seed") {
            cfg.seed = std::stoull(valeur_suivante(arg));
        } else if (arg == "--gui") {
            const int gui_int = std::stoi(valeur_suivante(arg));
            if (gui_int != 0 && gui_int != 1) {
                throw std::invalid_argument("--gui accepte uniquement 0 ou 1.");
            }
            cfg.gui = (gui_int == 1);
        } else if (arg == "--help") {
            afficher_aide(argv[0]);
            std::exit(0);
        } else {
            throw std::invalid_argument("Option inconnue: " + arg);
        }
    }

    if (cfg.nb_ants == 0) {
        throw std::invalid_argument("--ants doit etre > 0.");
    }
    if (cfg.max_iters == 0) {
        throw std::invalid_argument("--iter doit etre > 0.");
    }
    if (cfg.repetitions <= 0) {
        throw std::invalid_argument("--rep doit etre > 0.");
    }
    if (cfg.alpha < 0.0 || cfg.alpha > 1.0) {
        throw std::invalid_argument("--alpha doit etre dans [0,1].");
    }
    if (cfg.beta < 0.0 || cfg.beta > 1.0) {
        throw std::invalid_argument("--beta doit etre dans [0,1].");
    }
    if (cfg.eps < 0.0 || cfg.eps > 1.0) {
        throw std::invalid_argument("--eps doit etre dans [0,1].");
    }

    return cfg;
}

void normaliser_terrain(fractal_land& land) {
    double max_val = std::numeric_limits<double>::lowest();
    double min_val = std::numeric_limits<double>::max();

    for (fractal_land::dim_t i = 0; i < land.dimensions(); ++i) {
        for (fractal_land::dim_t j = 0; j < land.dimensions(); ++j) {
            max_val = std::max(max_val, land(i, j));
            min_val = std::min(min_val, land(i, j));
        }
    }

    const double delta = max_val - min_val;
    if (delta <= std::numeric_limits<double>::epsilon()) {
        for (fractal_land::dim_t i = 0; i < land.dimensions(); ++i) {
            for (fractal_land::dim_t j = 0; j < land.dimensions(); ++j) {
                land(i, j) = 0.0;
            }
        }
        return;
    }

    for (fractal_land::dim_t i = 0; i < land.dimensions(); ++i) {
        for (fractal_land::dim_t j = 0; j < land.dimensions(); ++j) {
            land(i, j) = (land(i, j) - min_val) / delta;
        }
    }
}

std::uint64_t bits_double(double valeur) {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &valeur, sizeof(double));
    return bits;
}

void melanger_checksum(std::uint64_t& checksum, std::uint64_t valeur) {
    checksum ^= valeur + 0x9e3779b97f4a7c15ULL + (checksum << 6U) + (checksum >> 2U);
}

std::uint64_t calculer_checksum_feromones(const pheronome& phen) {
    std::uint64_t checksum = 1469598103934665603ULL;
    const std::size_t dim = phen.dimensions();
    for (std::size_t i = 0; i < dim; ++i) {
        for (std::size_t j = 0; j < dim; ++j) {
            const auto& cellule = phen(i, j);
            melanger_checksum(checksum, bits_double(cellule[0]));
            melanger_checksum(checksum, bits_double(cellule[1]));
        }
    }
    return checksum;
}

double calculer_moyenne(const std::vector<double>& valeurs) {
    const double somme = std::accumulate(valeurs.begin(), valeurs.end(), 0.0);
    return somme / static_cast<double>(valeurs.size());
}

double calculer_ecart_type(const std::vector<double>& valeurs, double moyenne) {
    const double somme_carres =
        std::inner_product(valeurs.begin(), valeurs.end(), valeurs.begin(), 0.0);
    double variance = somme_carres / static_cast<double>(valeurs.size()) - moyenne * moyenne;
    if (variance < 0.0 && variance > -1e-12) {
        variance = 0.0;
    }
    return std::sqrt(std::max(0.0, variance));
}

PartitionFourmis construire_partition(std::size_t nb_ants_total, int rank, int size) {
    const std::size_t base = nb_ants_total / static_cast<std::size_t>(size);
    const std::size_t reste = nb_ants_total % static_cast<std::size_t>(size);

    PartitionFourmis part;
    if (static_cast<std::size_t>(rank) < reste) {
        part.nb_ants = base + 1ULL;
        part.debut = static_cast<std::size_t>(rank) * part.nb_ants;
    } else {
        part.nb_ants = base;
        part.debut = reste * (base + 1ULL) + (static_cast<std::size_t>(rank) - reste) * base;
    }
    part.fin = part.debut + part.nb_ants;
    return part;
}

std::vector<int> construire_counts(std::size_t nb_ants_total, int size) {
    std::vector<int> counts(static_cast<std::size_t>(size), 0);
    for (int r = 0; r < size; ++r) {
        const auto part = construire_partition(nb_ants_total, r, size);
        if (part.nb_ants > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::overflow_error("Trop de fourmis pour MPI_Gatherv.");
        }
        counts[static_cast<std::size_t>(r)] = static_cast<int>(part.nb_ants);
    }
    return counts;
}

std::vector<int> construire_displacements(const std::vector<int>& counts) {
    std::vector<int> displs(counts.size(), 0);
    for (std::size_t i = 1; i < counts.size(); ++i) {
        displs[i] = displs[i - 1] + counts[i - 1];
    }
    return displs;
}

void avancer_seed_initialisation(const fractal_land& land, std::size_t nb_ants_skip, std::size_t& seed) {
    for (std::size_t i = 0; i < nb_ants_skip; ++i) {
        (void)rand_int32(1, static_cast<std::int32_t>(land.dimensions() - 2), seed);
        (void)rand_int32(1, static_cast<std::int32_t>(land.dimensions() - 2), seed);
    }
}

void initialiser_partition_fourmis(const fractal_land& land,
                                   const Config& cfg,
                                   const PartitionFourmis& part,
                                   AntSwarm& ants) {
    ants.resize(part.nb_ants);
    std::size_t seed_local = cfg.seed;
    avancer_seed_initialisation(land, part.debut, seed_local);
    ants.initialiser_positions_aleatoires(land, seed_local);
}

std::uint64_t calculer_checksum_global(const pheronome& phen,
                                       const AntSwarm& ants_local,
                                       std::size_t nb_ants_total,
                                       int rank,
                                       const std::vector<int>& counts,
                                       const std::vector<int>& displs) {
    const int local_count = static_cast<int>(ants_local.size());

    std::vector<int> local_x(static_cast<std::size_t>(local_count), 0);
    std::vector<int> local_y(static_cast<std::size_t>(local_count), 0);
    std::vector<unsigned char> local_loaded(static_cast<std::size_t>(local_count), 0U);
    for (int i = 0; i < local_count; ++i) {
        const auto idx = static_cast<std::size_t>(i);
        local_x[idx] = ants_local.x_at(idx);
        local_y[idx] = ants_local.y_at(idx);
        local_loaded[idx] = static_cast<unsigned char>(ants_local.is_loaded_at(idx) ? 1U : 0U);
    }

    std::vector<int> global_x;
    std::vector<int> global_y;
    std::vector<unsigned char> global_loaded;
    if (rank == 0) {
        global_x.resize(nb_ants_total, 0);
        global_y.resize(nb_ants_total, 0);
        global_loaded.resize(nb_ants_total, 0U);
    }

    MPI_Gatherv(local_x.empty() ? nullptr : local_x.data(),
                local_count,
                MPI_INT,
                (rank == 0 && !global_x.empty()) ? global_x.data() : nullptr,
                rank == 0 ? counts.data() : nullptr,
                rank == 0 ? displs.data() : nullptr,
                MPI_INT,
                0,
                MPI_COMM_WORLD);

    MPI_Gatherv(local_y.empty() ? nullptr : local_y.data(),
                local_count,
                MPI_INT,
                (rank == 0 && !global_y.empty()) ? global_y.data() : nullptr,
                rank == 0 ? counts.data() : nullptr,
                rank == 0 ? displs.data() : nullptr,
                MPI_INT,
                0,
                MPI_COMM_WORLD);

    MPI_Gatherv(local_loaded.empty() ? nullptr : local_loaded.data(),
                local_count,
                MPI_UNSIGNED_CHAR,
                (rank == 0 && !global_loaded.empty()) ? global_loaded.data() : nullptr,
                rank == 0 ? counts.data() : nullptr,
                rank == 0 ? displs.data() : nullptr,
                MPI_UNSIGNED_CHAR,
                0,
                MPI_COMM_WORLD);

    if (rank != 0) {
        return 0ULL;
    }

    std::uint64_t checksum = calculer_checksum_feromones(phen);
    for (std::size_t i = 0; i < nb_ants_total; ++i) {
        melanger_checksum(checksum, static_cast<std::uint64_t>(static_cast<std::uint32_t>(global_x[i])));
        melanger_checksum(checksum, static_cast<std::uint64_t>(static_cast<std::uint32_t>(global_y[i])));
        melanger_checksum(checksum, global_loaded[i] != 0U ? 0xA5A5A5A5ULL : 0x5A5A5A5AULL);
    }
    return checksum;
}

ReferenceBaseline charger_reference_baseline_q2(const std::vector<std::filesystem::path>& chemins) {
    ReferenceBaseline reference;

    std::filesystem::path chemin_valide;
    for (const auto& candidat : chemins) {
        if (std::filesystem::exists(candidat)) {
            chemin_valide = candidat;
            break;
        }
    }

    if (chemin_valide.empty()) {
        reference.message = "fichier de reference Q2 absent";
        return reference;
    }

    std::ifstream fichier(chemin_valide);
    if (!fichier.is_open()) {
        reference.message = "impossible d'ouvrir la reference Q2";
        return reference;
    }

    std::string ligne;
    if (!std::getline(fichier, ligne)) {
        reference.message = "reference Q2 vide";
        return reference;
    }

    while (std::getline(fichier, ligne)) {
        if (ligne.empty()) {
            continue;
        }

        std::stringstream flux(ligne);
        std::string col_rep;
        std::string col_first_it;
        std::string col_food_kpi;
        std::string col_checksum;

        if (!std::getline(flux, col_rep, ',')) {
            continue;
        }
        if (col_rep == "MEAN" || col_rep == "STD_DEV") {
            break;
        }
        if (!std::getline(flux, col_first_it, ',')) {
            continue;
        }
        if (!std::getline(flux, col_food_kpi, ',')) {
            continue;
        }
        if (!std::getline(flux, col_checksum, ',')) {
            continue;
        }

        try {
            (void)std::stoull(col_rep);
            reference.food_kpi.push_back(static_cast<std::size_t>(std::stoull(col_food_kpi)));
            reference.checksums.push_back(static_cast<std::uint64_t>(std::stoull(col_checksum)));
        } catch (...) {
            continue;
        }
    }

    if (reference.food_kpi.empty() || reference.checksums.empty()) {
        reference.message = "aucune repetition exploitable dans la reference Q2";
        return reference;
    }

    reference.disponible = true;
    reference.chemin_source = chemin_valide;
    reference.message = "reference Q2 chargee";
    return reference;
}

void afficher_verification_baseline_q2(int mpi_size,
                                       const std::vector<std::size_t>& food_kpi_q4,
                                       const std::vector<std::uint64_t>& checksums_q4) {
    if (mpi_size != 1) {
        std::cout << "Verification baseline Q2: A_VERIFIER"
                  << " | raison=verification reservee a np=1\n";
        return;
    }

    const std::vector<std::filesystem::path> candidats = {
        "Q2/results/Q2_timings_breakdown.csv",
        "../Q2/results/Q2_timings_breakdown.csv",
        "results/Q2_timings_breakdown.csv",
        "../results/Q2_timings_breakdown.csv"
    };
    const auto reference = charger_reference_baseline_q2(candidats);

    if (!reference.disponible) {
        std::cout << "Verification baseline Q2: A_VERIFIER"
                  << " | raison=" << reference.message << "\n";
        return;
    }

    const std::size_t nb_compare = std::min(
        food_kpi_q4.size(),
        std::min(reference.food_kpi.size(), reference.checksums.size())
    );
    if (nb_compare == 0) {
        std::cout << "Verification baseline Q2: A_VERIFIER"
                  << " | raison=aucune repetition commune\n";
        return;
    }

    bool all_match = true;
    for (std::size_t i = 0; i < nb_compare; ++i) {
        if (food_kpi_q4[i] != reference.food_kpi[i] || checksums_q4[i] != reference.checksums[i]) {
            all_match = false;
            break;
        }
    }

    if (!all_match) {
        std::cout << "Verification baseline Q2: MISMATCH"
                  << " | source=" << reference.chemin_source.string()
                  << " | repetitions_comparees=" << nb_compare << "/" << food_kpi_q4.size() << "\n";
        for (std::size_t i = 0; i < nb_compare; ++i) {
            if (food_kpi_q4[i] != reference.food_kpi[i] || checksums_q4[i] != reference.checksums[i]) {
                std::cout << "  Rep " << (i + 1)
                          << " | Q4(Food_KPI=" << food_kpi_q4[i]
                          << ", Checksum=" << checksums_q4[i]
                          << ") vs Q2(Food_KPI=" << reference.food_kpi[i]
                          << ", Checksum=" << reference.checksums[i] << ")\n";
            }
        }
        return;
    }

    std::cout << "Verification baseline Q2: MATCH_VALUE_ONLY"
              << " | source=" << reference.chemin_source.string()
              << " | repetitions_comparees=" << nb_compare << "/" << food_kpi_q4.size()
              << " | raison=reference_csv_sans_metadonnees_parametres\n";
}

MesuresLocales executer_repetition(const Config& cfg,
                                   const ParametresFractal& params_fractal,
                                   bool activer_gui,
                                   int rank,
                                   int size,
                                   const PartitionFourmis& partition,
                                   const std::vector<int>& counts,
                                   const std::vector<int>& displs,
                                   bool calculer_checksum) {
    const bool gui_local = activer_gui && (rank == 0);

    fractal_land land(
        params_fractal.log2_sous_grille,
        params_fractal.nb_graines,
        1.0,
        static_cast<int>(cfg.seed & 0x7fffffffU)
    );
    normaliser_terrain(land);

    const int borne_min = 1;
    const int borne_max = static_cast<int>(land.dimensions()) - 2;
    position_t pos_nest{static_cast<int>(land.dimensions() / 2), static_cast<int>(land.dimensions() / 2)};
    position_t pos_food{
        std::max(borne_min, static_cast<int>((3 * land.dimensions()) / 4)),
        std::max(borne_min, static_cast<int>((3 * land.dimensions()) / 4))
    };
    pos_food.x = std::min(pos_food.x, borne_max);
    pos_food.y = std::min(pos_food.y, borne_max);
    if (pos_food == pos_nest) {
        pos_food.x = std::min(pos_food.x + 1, borne_max);
    }

    AntSwarm::set_exploration_coef(cfg.eps);
    AntSwarm ants(partition.nb_ants);
    initialiser_partition_fourmis(land, cfg, partition, ants);

    pheronome phen(land.dimensions(), pos_food, pos_nest, cfg.alpha, cfg.beta);

    std::unique_ptr<Window> win;
    std::unique_ptr<Renderer> renderer;
    if (gui_local) {
        win = std::make_unique<Window>(
            "Simulation ACO MPI",
            static_cast<int>(2 * land.dimensions() + 10),
            static_cast<int>(land.dimensions() + 266)
        );
        renderer = std::make_unique<Renderer>(land, phen, pos_nest, pos_food, ants);
    }

    MesuresLocales mesures;
    std::size_t food_quantity = 0;
    std::size_t it = 0;
    bool cont_loop = true;
    bool premiere_nourriture_detectee = false;

    const std::size_t phen_count = phen.stride() * phen.stride() * 2ULL;
    if (phen_count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("Carte de pheromones trop grande pour MPI_Allreduce.");
    }

    const auto t_total_start = MPI_Wtime();
    while (cont_loop && it < cfg.max_iters) {
        ++it;

        if (gui_local) {
            SDL_Event event;
            while (SDL_PollEvent(&event)) {
                if (event.type == SDL_QUIT) {
                    cont_loop = false;
                }
            }
        }

        std::size_t batch_food = 0;
        auto t0 = MPI_Wtime();
        ants.advance_one(phen, land, pos_food, pos_nest, batch_food);
        auto t1 = MPI_Wtime();
        mesures.t_compute_move_ms += (t1 - t0) * 1000.0;

        unsigned long local_food = static_cast<unsigned long>(batch_food);
        unsigned long global_food = 0UL;
        t0 = MPI_Wtime();
        MPI_Allreduce(&local_food, &global_food, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
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

        t0 = MPI_Wtime();
        MPI_Allreduce(
            MPI_IN_PLACE,
            phen.data(),
            static_cast<int>(phen_count),
            MPI_DOUBLE,
            MPI_MAX,
            MPI_COMM_WORLD
        );
        t1 = MPI_Wtime();
        mesures.t_comm_pher_ms += (t1 - t0) * 1000.0;

        if (gui_local) {
            t0 = MPI_Wtime();
            renderer->display(*win, food_quantity);
            win->blit();
            t1 = MPI_Wtime();
            mesures.t_render_ms += (t1 - t0) * 1000.0;
        }

        if (!premiere_nourriture_detectee && food_quantity > 0) {
            mesures.premiere_iteration_nourriture = it;
            premiere_nourriture_detectee = true;
        }
    }
    const auto t_total_end = MPI_Wtime();

    mesures.food_kpi = food_quantity;
    mesures.t_compute_ms =
        mesures.t_compute_move_ms + mesures.t_compute_evap_ms +
        mesures.t_compute_update_ms + mesures.t_render_ms;
    mesures.t_comm_ms = mesures.t_comm_food_ms + mesures.t_comm_pher_ms;
    mesures.t_total_ms = (t_total_end - t_total_start) * 1000.0;
    mesures.ecart_sanity_ms = std::abs(mesures.t_total_ms - (mesures.t_compute_ms + mesures.t_comm_ms));

    if (calculer_checksum) {
        mesures.checksum_global = calculer_checksum_global(
            phen,
            ants,
            cfg.nb_ants,
            rank,
            counts,
            displs
        );
    }

    return mesures;
}

int main(int argc, char* argv[]) {
    int rank = 0;
    int size = 1;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    try {
        Config cfg = parse_args(argc, argv);
        const ParametresFractal params_fractal = construire_parametres_fractal(cfg.grid_size);
        const auto partition = construire_partition(cfg.nb_ants, rank, size);
        const auto counts = construire_counts(cfg.nb_ants, size);
        const auto displs = construire_displacements(counts);

        if (size != 1 && cfg.gui && rank == 0) {
            std::cout << "GUI desactivee: mode MPI multi-processus.\n";
        }
        if (size != 1) {
            cfg.gui = false;
        }
        if (rank != 0) {
            cfg.gui = false;
        }

        if (rank == 0) {
            if (params_fractal.taille_effective != cfg.grid_size) {
                std::cerr << "Avertissement: la taille de grille demandee " << cfg.grid_size
                          << " est ajustee a " << params_fractal.taille_effective << ".\n";
            }

            std::cout << "--- Q4 MPI pur (environnement replique) ---\n";
            std::cout << "Grille: " << params_fractal.taille_effective
                      << " | Fourmis: " << cfg.nb_ants
                      << " | Iterations max: " << cfg.max_iters
                      << " | Repetitions: " << cfg.repetitions
                      << " | GUI: " << (cfg.gui ? "ON" : "OFF")
                      << " | MPI np: " << size << "\n";
            std::cout << "Warm-up en cours...\n";
        }

        (void)executer_repetition(
            cfg,
            params_fractal,
            false,
            rank,
            size,
            partition,
            counts,
            displs,
            false
        );

        bool sdl_initialisee = false;
        if (cfg.gui && rank == 0) {
            if (SDL_Init(SDL_INIT_VIDEO) != 0) {
                throw std::runtime_error("Impossible d'initialiser SDL en mode GUI.");
            }
            sdl_initialisee = true;
        }

        std::vector<std::size_t> premiere_iteration_nourriture;
        std::vector<std::size_t> food_kpi;
        std::vector<std::uint64_t> checksums;
        std::vector<double> t_compute_move;
        std::vector<double> t_compute_evap;
        std::vector<double> t_compute_update;
        std::vector<double> t_compute;
        std::vector<double> t_comm_food;
        std::vector<double> t_comm_pher;
        std::vector<double> t_comm;
        std::vector<double> t_total;
        std::vector<double> t_sanity_ecart;

        if (rank == 0) {
            const auto rep_count = static_cast<std::size_t>(cfg.repetitions);
            premiere_iteration_nourriture.reserve(rep_count);
            food_kpi.reserve(rep_count);
            checksums.reserve(rep_count);
            t_compute_move.reserve(rep_count);
            t_compute_evap.reserve(rep_count);
            t_compute_update.reserve(rep_count);
            t_compute.reserve(rep_count);
            t_comm_food.reserve(rep_count);
            t_comm_pher.reserve(rep_count);
            t_comm.reserve(rep_count);
            t_total.reserve(rep_count);
            t_sanity_ecart.reserve(rep_count);
        }

        for (int rep = 0; rep < cfg.repetitions; ++rep) {
            MesuresLocales local = executer_repetition(
                cfg,
                params_fractal,
                cfg.gui,
                rank,
                size,
                partition,
                counts,
                displs,
                true
            );

            double g_compute_move = 0.0;
            double g_compute_evap = 0.0;
            double g_compute_update = 0.0;
            double g_compute = 0.0;
            double g_comm_food = 0.0;
            double g_comm_pher = 0.0;
            double g_comm = 0.0;
            double g_total = 0.0;
            double g_ecart_sanity = 0.0;

            MPI_Reduce(&local.t_compute_move_ms, &g_compute_move, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.t_compute_evap_ms, &g_compute_evap, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.t_compute_update_ms, &g_compute_update, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.t_compute_ms, &g_compute, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.t_comm_food_ms, &g_comm_food, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.t_comm_pher_ms, &g_comm_pher, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.t_comm_ms, &g_comm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.t_total_ms, &g_total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&local.ecart_sanity_ms, &g_ecart_sanity, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

            unsigned long long food_local = static_cast<unsigned long long>(local.food_kpi);
            unsigned long long food_global = 0ULL;
            MPI_Reduce(&food_local, &food_global, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

            unsigned long long first_it_local =
                static_cast<unsigned long long>(local.premiere_iteration_nourriture);
            unsigned long long first_it_global = 0ULL;
            MPI_Reduce(
                &first_it_local,
                &first_it_global,
                1,
                MPI_UNSIGNED_LONG_LONG,
                MPI_MAX,
                0,
                MPI_COMM_WORLD
            );

            if (rank == 0) {
                premiere_iteration_nourriture.push_back(static_cast<std::size_t>(first_it_global));
                food_kpi.push_back(static_cast<std::size_t>(food_global));
                checksums.push_back(local.checksum_global);
                t_compute_move.push_back(g_compute_move);
                t_compute_evap.push_back(g_compute_evap);
                t_compute_update.push_back(g_compute_update);
                t_compute.push_back(g_compute);
                t_comm_food.push_back(g_comm_food);
                t_comm_pher.push_back(g_comm_pher);
                t_comm.push_back(g_comm);
                t_total.push_back(g_total);
                t_sanity_ecart.push_back(g_ecart_sanity);

                const double tolerance = std::max(1.0, 0.05 * g_total);
                const bool sanity_ok = g_ecart_sanity <= tolerance;

                std::cout << "Rep " << (rep + 1) << "/" << cfg.repetitions
                          << " | T_compute=" << std::fixed << std::setprecision(6) << g_compute << " ms"
                          << " | T_comm=" << g_comm << " ms"
                          << " | T_total=" << g_total << " ms"
                          << " | Food_KPI=" << food_global
                          << " | Checksum=" << local.checksum_global
                          << " | Sanity=" << (sanity_ok ? "OK" : "A_VERIFIER")
                          << " (ecart_max_ranks=" << g_ecart_sanity << " ms)\n";
            }
        }

        if (sdl_initialisee) {
            SDL_Quit();
        }

        if (rank == 0) {
            std::filesystem::create_directories("results");
            std::ofstream csv_file("results/Q4_timings_breakdown.csv");
            if (!csv_file.is_open()) {
                throw std::runtime_error("Impossible d'ecrire results/Q4_timings_breakdown.csv");
            }

            csv_file << "Repetition,MPI_Procs,First_Iteration,Food_KPI,Checksum,"
                        "T_compute_move_ms,T_compute_evap_ms,T_compute_update_ms,T_compute_ms,"
                        "T_comm_food_ms,T_comm_pher_ms,T_comm_ms,T_total_ms,Sanity_ecart_ms\n";

            for (std::size_t i = 0; i < t_total.size(); ++i) {
                csv_file << (i + 1) << ","
                         << size << ","
                         << premiere_iteration_nourriture[i] << ","
                         << food_kpi[i] << ","
                         << checksums[i] << ","
                         << std::fixed << std::setprecision(6)
                         << t_compute_move[i] << ","
                         << t_compute_evap[i] << ","
                         << t_compute_update[i] << ","
                         << t_compute[i] << ","
                         << t_comm_food[i] << ","
                         << t_comm_pher[i] << ","
                         << t_comm[i] << ","
                         << t_total[i] << ","
                         << t_sanity_ecart[i] << "\n";
            }

            const double mean_compute_move = calculer_moyenne(t_compute_move);
            const double mean_compute_evap = calculer_moyenne(t_compute_evap);
            const double mean_compute_update = calculer_moyenne(t_compute_update);
            const double mean_compute = calculer_moyenne(t_compute);
            const double mean_comm_food = calculer_moyenne(t_comm_food);
            const double mean_comm_pher = calculer_moyenne(t_comm_pher);
            const double mean_comm = calculer_moyenne(t_comm);
            const double mean_total = calculer_moyenne(t_total);
            const double mean_sanity = calculer_moyenne(t_sanity_ecart);

            const double std_compute_move = calculer_ecart_type(t_compute_move, mean_compute_move);
            const double std_compute_evap = calculer_ecart_type(t_compute_evap, mean_compute_evap);
            const double std_compute_update = calculer_ecart_type(t_compute_update, mean_compute_update);
            const double std_compute = calculer_ecart_type(t_compute, mean_compute);
            const double std_comm_food = calculer_ecart_type(t_comm_food, mean_comm_food);
            const double std_comm_pher = calculer_ecart_type(t_comm_pher, mean_comm_pher);
            const double std_comm = calculer_ecart_type(t_comm, mean_comm);
            const double std_total = calculer_ecart_type(t_total, mean_total);
            const double std_sanity = calculer_ecart_type(t_sanity_ecart, mean_sanity);

            csv_file << "MEAN,"
                     << size << ",,,,"
                     << std::fixed << std::setprecision(6)
                     << mean_compute_move << ","
                     << mean_compute_evap << ","
                     << mean_compute_update << ","
                     << mean_compute << ","
                     << mean_comm_food << ","
                     << mean_comm_pher << ","
                     << mean_comm << ","
                     << mean_total << ","
                     << mean_sanity << "\n";

            csv_file << "STD_DEV,"
                     << size << ",,,,"
                     << std::fixed << std::setprecision(6)
                     << std_compute_move << ","
                     << std_compute_evap << ","
                     << std_compute_update << ","
                     << std_compute << ","
                     << std_comm_food << ","
                     << std_comm_pher << ","
                     << std_comm << ","
                     << std_total << ","
                     << std_sanity << "\n";

            std::cout << "\n=== Resultats Q4 MPI ===\n";
            std::cout << "Moyennes (ms): T_compute=" << mean_compute
                      << " | T_comm=" << mean_comm
                      << " | T_total=" << mean_total << "\n";
            std::cout << "CSV genere: results/Q4_timings_breakdown.csv\n";

            afficher_verification_baseline_q2(size, food_kpi, checksums);
        }

        MPI_Finalize();
        return 0;
    } catch (const std::exception& e) {
        if (rank == 0) {
            std::cerr << "Erreur: " << e.what() << "\n";
            afficher_aide(argv[0]);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
}
