#include <algorithm>
#include <chrono>
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
#include <stdexcept>
#include <string>
#include <vector>

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

struct MesuresRepetition {
    double t_move_ants_ms = 0.0;
    double t_evap_ms = 0.0;
    double t_pher_update_ms = 0.0;
    double t_render_ms = 0.0;
    double t_total_ms = 0.0;
    std::size_t premiere_iteration_nourriture = 0;
    std::size_t food_kpi = 0;
    std::size_t iterations_executees = 0;
    std::uint64_t checksum = 0;
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
        << "  --iter N        Nombre d'iterations (defaut: 3000)\n"
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

std::uint64_t calculer_checksum(const pheronome& phen, const std::vector<ant>& ants) {
    std::uint64_t checksum = 1469598103934665603ULL;
    const std::size_t dim = phen.dimensions();

    for (std::size_t i = 0; i < dim; ++i) {
        for (std::size_t j = 0; j < dim; ++j) {
            const auto& cellule = phen(i, j);
            melanger_checksum(checksum, bits_double(cellule[0]));
            melanger_checksum(checksum, bits_double(cellule[1]));
        }
    }

    for (const auto& a : ants) {
        const auto& pos = a.get_position();
        melanger_checksum(checksum, static_cast<std::uint64_t>(static_cast<std::uint32_t>(pos.x)));
        melanger_checksum(checksum, static_cast<std::uint64_t>(static_cast<std::uint32_t>(pos.y)));
        melanger_checksum(checksum, a.is_loaded() ? 0xA5A5A5A5ULL : 0x5A5A5A5AULL);
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

MesuresRepetition executer_repetition(const Config& cfg,
                                      const ParametresFractal& params_fractal,
                                      bool activer_gui)
{
    using horloge = std::chrono::high_resolution_clock;

    std::size_t current_seed = cfg.seed;

    fractal_land land(
        params_fractal.log2_sous_grille,
        params_fractal.nb_graines,
        1.0,
        static_cast<int>(current_seed & 0x7fffffffU)
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

    ant::set_exploration_coef(cfg.eps);

    std::vector<ant> ants;
    ants.reserve(cfg.nb_ants);
    auto gen_ant_pos = [&land, &current_seed]() {
        return rand_int32(1, static_cast<std::int32_t>(land.dimensions() - 2), current_seed);
    };
    for (std::size_t i = 0; i < cfg.nb_ants; ++i) {
        ants.emplace_back(position_t{gen_ant_pos(), gen_ant_pos()}, current_seed);
    }

    pheronome phen(land.dimensions(), pos_food, pos_nest, cfg.alpha, cfg.beta);

    std::unique_ptr<Window> win;
    std::unique_ptr<Renderer> renderer;
    if (activer_gui) {
        win = std::make_unique<Window>(
            "Simulation ACO",
            static_cast<int>(2 * land.dimensions() + 10),
            static_cast<int>(land.dimensions() + 266)
        );
        renderer = std::make_unique<Renderer>(land, phen, pos_nest, pos_food, ants);
    }

    MesuresRepetition mesures;
    std::size_t food_quantity = 0;
    std::size_t it = 0;
    bool cont_loop = true;
    bool premiere_nourriture_detectee = false;

    const auto debut_total = horloge::now();
    while (cont_loop && it < cfg.max_iters) {
        ++it;

        if (activer_gui) {
            SDL_Event event;
            while (SDL_PollEvent(&event)) {
                if (event.type == SDL_QUIT) {
                    cont_loop = false;
                }
            }
        }

        const auto debut_move = horloge::now();
        for (std::size_t i = 0; i < ants.size(); ++i) {
            ants[i].advance(phen, land, pos_food, pos_nest, food_quantity);
        }
        const auto fin_move = horloge::now();
        mesures.t_move_ants_ms += std::chrono::duration<double, std::milli>(fin_move - debut_move).count();

        const auto debut_evap = horloge::now();
        phen.do_evaporation();
        const auto fin_evap = horloge::now();
        mesures.t_evap_ms += std::chrono::duration<double, std::milli>(fin_evap - debut_evap).count();

        const auto debut_update = horloge::now();
        phen.update();
        const auto fin_update = horloge::now();
        mesures.t_pher_update_ms +=
            std::chrono::duration<double, std::milli>(fin_update - debut_update).count();

        if (activer_gui) {
            const auto debut_render = horloge::now();
            renderer->display(*win, food_quantity);
            win->blit();
            const auto fin_render = horloge::now();
            mesures.t_render_ms +=
                std::chrono::duration<double, std::milli>(fin_render - debut_render).count();
        }

        if (!premiere_nourriture_detectee && food_quantity > 0) {
            mesures.premiere_iteration_nourriture = it;
            premiere_nourriture_detectee = true;
        }
    }
    const auto fin_total = horloge::now();

    mesures.t_total_ms = std::chrono::duration<double, std::milli>(fin_total - debut_total).count();
    mesures.food_kpi = food_quantity;
    mesures.iterations_executees = it;
    mesures.checksum = calculer_checksum(phen, ants);

    return mesures;
}

void afficher_sanity_check(const MesuresRepetition& mesures, bool mode_gui) {
    const double somme_etapes =
        mesures.t_move_ants_ms + mesures.t_evap_ms + mesures.t_pher_update_ms + mesures.t_render_ms;
    const double ecart = std::abs(mesures.t_total_ms - somme_etapes);
    const double tolerance = std::max(1.0, 0.05 * mesures.t_total_ms);
    const bool check_total = ecart <= tolerance;

    std::cout << "Sanity Check: T_total ~= T_move_ants + T_evap + T_pher_update + T_render"
              << " | ecart=" << std::fixed << std::setprecision(6) << ecart << " ms"
              << " | statut=" << (check_total ? "OK" : "A_VERIFIER") << "\n";

    if (!mode_gui) {
        const bool check_render = std::abs(mesures.t_render_ms) <= 1e-9;
        std::cout << "Sanity Check (--gui 0): T_render ~ 0"
                  << " | valeur=" << std::fixed << std::setprecision(6) << mesures.t_render_ms << " ms"
                  << " | statut=" << (check_render ? "OK" : "A_VERIFIER") << "\n";
    }
}

int main(int argc, char* argv[]) {
    try {
        const Config cfg = parse_args(argc, argv);
        const ParametresFractal params_fractal = construire_parametres_fractal(cfg.grid_size);

        if (params_fractal.taille_effective != cfg.grid_size) {
            std::cerr << "Avertissement: la taille de grille demandee " << cfg.grid_size
                      << " est ajustee a " << params_fractal.taille_effective << ".\n";
        }

        std::cout << "--- Q1 Instrumentation monocoeur (timings par etape) ---\n";
        std::cout << "Grille: " << params_fractal.taille_effective
                  << " | Fourmis: " << cfg.nb_ants
                  << " | Iterations max: " << cfg.max_iters
                  << " | Repetitions: " << cfg.repetitions
                  << " | GUI: " << (cfg.gui ? "ON" : "OFF") << "\n";

        std::cout << "Exécution du Warm-up pour chauffer le cache...\n";
        (void)executer_repetition(cfg, params_fractal, false);

        bool sdl_initialisee = false;
        if (cfg.gui) {
            if (SDL_Init(SDL_INIT_VIDEO) != 0) {
                throw std::runtime_error("Impossible d'initialiser SDL en mode GUI.");
            }
            sdl_initialisee = true;
        }

        std::vector<double> t_move_ants;
        std::vector<double> t_evap;
        std::vector<double> t_pher_update;
        std::vector<double> t_render;
        std::vector<double> t_total;
        std::vector<std::size_t> premiere_iteration_nourriture;
        std::vector<std::size_t> food_kpi;
        std::vector<std::uint64_t> checksums;

        t_move_ants.reserve(static_cast<std::size_t>(cfg.repetitions));
        t_evap.reserve(static_cast<std::size_t>(cfg.repetitions));
        t_pher_update.reserve(static_cast<std::size_t>(cfg.repetitions));
        t_render.reserve(static_cast<std::size_t>(cfg.repetitions));
        t_total.reserve(static_cast<std::size_t>(cfg.repetitions));
        premiere_iteration_nourriture.reserve(static_cast<std::size_t>(cfg.repetitions));
        food_kpi.reserve(static_cast<std::size_t>(cfg.repetitions));
        checksums.reserve(static_cast<std::size_t>(cfg.repetitions));

        for (int rep = 0; rep < cfg.repetitions; ++rep) {
            const MesuresRepetition mesures = executer_repetition(cfg, params_fractal, cfg.gui);

            t_move_ants.push_back(mesures.t_move_ants_ms);
            t_evap.push_back(mesures.t_evap_ms);
            t_pher_update.push_back(mesures.t_pher_update_ms);
            t_render.push_back(mesures.t_render_ms);
            t_total.push_back(mesures.t_total_ms);
            premiere_iteration_nourriture.push_back(mesures.premiere_iteration_nourriture);
            food_kpi.push_back(mesures.food_kpi);
            checksums.push_back(mesures.checksum);

            std::cout << "Rep " << (rep + 1) << "/" << cfg.repetitions
                      << " | T_move_ants=" << std::fixed << std::setprecision(6) << mesures.t_move_ants_ms << " ms"
                      << " | T_evap=" << mesures.t_evap_ms << " ms"
                      << " | T_pher_update=" << mesures.t_pher_update_ms << " ms"
                      << " | T_render=" << mesures.t_render_ms << " ms"
                      << " | T_total=" << mesures.t_total_ms << " ms"
                      << " | First_Iteration=" << mesures.premiere_iteration_nourriture
                      << " | Food_KPI=" << mesures.food_kpi
                      << " | Checksum=" << mesures.checksum
                      << "\n";

            afficher_sanity_check(mesures, cfg.gui);
        }

        if (sdl_initialisee) {
            SDL_Quit();
        }

        std::filesystem::create_directories("results");
        std::ofstream csv_file("results/Q1_timings_breakdown.csv");
        if (!csv_file.is_open()) {
            std::cerr << "Erreur: impossible d'ecrire results/Q1_timings_breakdown.csv\n";
            return 1;
        }

        csv_file << "Repetition,Time_ms,First_Iteration,Food_KPI,Checksum,T_move_ants,T_evap,T_pher_update,T_render,T_total\n";
        for (std::size_t i = 0; i < t_total.size(); ++i) {
            csv_file << (i + 1) << ","
                     << std::fixed << std::setprecision(6) << t_total[i] << ","
                     << premiere_iteration_nourriture[i] << ","
                     << food_kpi[i] << ","
                     << checksums[i] << ","
                     << t_move_ants[i] << ","
                     << t_evap[i] << ","
                     << t_pher_update[i] << ","
                     << t_render[i] << ","
                     << t_total[i] << "\n";
        }

        const double mean_move = calculer_moyenne(t_move_ants);
        const double mean_evap = calculer_moyenne(t_evap);
        const double mean_update = calculer_moyenne(t_pher_update);
        const double mean_render = calculer_moyenne(t_render);
        const double mean_total = calculer_moyenne(t_total);

        const double std_move = calculer_ecart_type(t_move_ants, mean_move);
        const double std_evap = calculer_ecart_type(t_evap, mean_evap);
        const double std_update = calculer_ecart_type(t_pher_update, mean_update);
        const double std_render = calculer_ecart_type(t_render, mean_render);
        const double std_total = calculer_ecart_type(t_total, mean_total);

        csv_file << "MEAN,"
                 << std::fixed << std::setprecision(6) << mean_total << ","
                 << ",,,"
                 << mean_move << ","
                 << mean_evap << ","
                 << mean_update << ","
                 << mean_render << ","
                 << mean_total << "\n";
        csv_file << "STD_DEV,"
                 << std::fixed << std::setprecision(6) << std_total << ","
                 << ",,,"
                 << std_move << ","
                 << std_evap << ","
                 << std_update << ","
                 << std_render << ","
                 << std_total << "\n";

        std::cout << "\n=== Resultats Q1 (timings par etape) ===\n";
        std::cout << "Moyennes (ms):"
                  << " move_ants=" << mean_move
                  << " | evap=" << mean_evap
                  << " | pher_update=" << mean_update
                  << " | render=" << mean_render
                  << " | total=" << mean_total << "\n";
        std::cout << "CSV genere: results/Q1_timings_breakdown.csv\n";

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Erreur: " << e.what() << "\n";
        afficher_aide(argv[0]);
        return 1;
    }
}
