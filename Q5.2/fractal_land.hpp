#ifndef _FRACTAL_LAND_HPP_
# define _FRACTAL_LAND_HPP_
// Génération d'un fractal pour le coût énergie ( altitude ) de déplacement d'une fourmie
// L'algorithme prend plusieurs paramètres :
// 0. Taille : Nombre de "cases" par direction ( pair d'entier )
// 1. Nombre de graînes : Nombre de points initiaux par direction dont on détermine l'altitude à l'initialisation de l'agorithme
// 2. Déviation : degré de variation de l'altitude en fonction de la distance
// 3. Graîne aléatoire : détermine le paysage à retrouver
# include <vector>
# include <utility>
# include <cassert>
#include "rand_generator.hpp"

/**
 * @brief Génère un paysage fractal à l'aide d'un algorithme pseudo-aléatoire 
 * @details 
 *     Génère de façon récursive un paysage fractal à l'aide de algorithme pseudo-aléatoire :
 *        1. Créée une grille de taille \f$nbSeeds*2^{log\_size}+1\f$ cases par directions
 *        2. Génère une altitude pour les cases ayant des indices i et j multiples de 
 *           \f$2^{log\_size}\f$ de telle sorte que le gradient d'altitude entre deux points ne dépasse pas la valeur deviation
 *           On considère alors les sous-grilles ayant pour indices mimimals : \f$Ib*2^{log\_size}\f$, \f$Jb*2^{log\_size}\f$ avec
 *           Ib et Jb compris entre 0 et nbSeeds ( compris ) et de tailles \f$2^{log\_size}\f$.
 *           On note n=log_size le niveau des sous--grilles initiales.
 *        3. Pour chaque sous--grille, on génère l'altitude des points d'indices locaux multiples de \f$2^{n-1}\f$ ormi pour les
 *           coins de la grille en respectant le gradient de déviation.
 *        4. Puis on considère les sous--grilles de niveau n-1 auxquelles on reapplique l'algorithme à partir de 3 et on s'arrête dès que
 *           le niveau de la grille atteint zéro.
 * @param log_size Le logarithme base 2 de la dimension de chaque sous-grille initiale
 * @param nbSeeds  Le nombre de sous-grilles initiales par direction
 * @param deviation La valeur maximale du gradient entre deux altitudes.
 * @param seed Graîne de génération aléatoire
 * @return Un tableau contenant la carte des altitudes en fonctions des indices i et j.
 */
class fractal_land
{
public:
    using dim_t = unsigned long;
    using container = std::vector<double>;
    fractal_land(const dim_t& ln2_dim, unsigned long nbSeeds, double deviation, int seed, int rank, int size);
    fractal_land(const fractal_land&) = delete;
    fractal_land(fractal_land&&) = delete;
    ~fractal_land() = default;

    double& operator()(dim_t i, dim_t j_global) {
        assert(is_local(j_global));
        dim_t j_local = j_global - m_row_start;
        return m_altitude[j_local * m_dimensions + i];
    }

    double operator()(dim_t i, dim_t j_global) const {
        assert(is_local(j_global));
        dim_t j_local = j_global - m_row_start;
        return m_altitude[j_local * m_dimensions + i];
    }

    double get_global(dim_t i, dim_t j) const {
        if (is_local(j)) return (*this)(i, j);
        return 0.0; 
    }

    bool is_local(dim_t y_global) const;
    dim_t dimensions() const { return m_dimensions; }
    dim_t local_height() const { return m_local_height; }
    dim_t row_start() const { return m_row_start; }
    
private:
    void compute_subgrid(int log_subgrid_dim, int iB, int jB, double deviation, std::size_t seed);

    dim_t m_dimensions;     
    dim_t m_row_start;      
    dim_t m_local_height;   
    container m_altitude;
};
#endif