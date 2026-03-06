#include "ant.hpp"
#include <algorithm>
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
    m_is_loaded.assign(nb_ants, false);
    m_seed.assign(nb_ants, 0);
}

void AntSwarm::initialiser_positions_aleatoires(const fractal_land& land, std::size_t& seed_global)
{
    auto gen_ant_pos = [&land, &seed_global]() {
        return rand_int32(1, static_cast<std::int32_t>(land.dimensions() - 2), seed_global);
    };

    for (std::size_t i = 0; i < size(); ++i) {
        m_x[i] = gen_ant_pos();
        m_y[i] = gen_ant_pos();
        m_is_loaded[i] = false;
        m_seed[i] = seed_global;
    }
}

void AntSwarm::advance_one(std::size_t ant_id, pheronome& phen, const fractal_land& land,
                           const position_t& pos_food, const position_t& pos_nest, std::size_t& cpteur_food)
{
    auto ant_choice = [this, ant_id]() mutable { return rand_double(0., 1., this->m_seed[ant_id]); };
    auto dir_choice = [this, ant_id]() mutable { return rand_int32(1, 4, this->m_seed[ant_id]); };
    double consumed_time = 0.;

    // Tant que la fourmi peut encore bouger dans le pas de temps imparti
    while ( consumed_time < 1. ) {
        // Si la fourmi est chargée, elle suit les phéromones de deuxième type, sinon ceux du premier.
        int        ind_pher    = ( m_is_loaded[ant_id] ? 1 : 0 );
        double     choix       = ant_choice( );
        position_t old_pos_ant{m_x[ant_id], m_y[ant_id]};
        position_t new_pos_ant = old_pos_ant;
        const position_t pos_gauche{new_pos_ant.x - 1, new_pos_ant.y};
        const position_t pos_droite{new_pos_ant.x + 1, new_pos_ant.y};
        const position_t pos_haut{new_pos_ant.x, new_pos_ant.y - 1};
        const position_t pos_bas{new_pos_ant.x, new_pos_ant.y + 1};
        double max_phen    = std::max( {phen[pos_gauche][ind_pher],
                                     phen[pos_droite][ind_pher],
                                     phen[pos_haut][ind_pher],
                                     phen[pos_bas][ind_pher]} );
        if ( ( choix > m_eps ) || ( max_phen <= 0. ) ) {
            do {
                new_pos_ant = old_pos_ant;
                int d = dir_choice();
                if ( d==1 ) new_pos_ant.x  -= 1;
                if ( d==2 ) new_pos_ant.y -= 1;
                if ( d==3 ) new_pos_ant.x  += 1;
                if ( d==4 ) new_pos_ant.y += 1;

            } while ( phen[new_pos_ant][ind_pher] == -1 );
        } else {
            // On choisit la case où le phéromone est le plus fort.
            if ( phen[pos_gauche][ind_pher] == max_phen )
                new_pos_ant.x -= 1;
            else if ( phen[pos_droite][ind_pher] == max_phen )
                new_pos_ant.x += 1;
            else if ( phen[pos_haut][ind_pher] == max_phen )
                new_pos_ant.y -= 1;
            else  // if (phen(new_pos_ant.first,new_pos_ant.second+1)[ind_pher] == max_phen)
                new_pos_ant.y += 1;
        }
        consumed_time += land( new_pos_ant.x, new_pos_ant.y);
        phen.mark_pheronome( new_pos_ant );
        m_x[ant_id] = new_pos_ant.x;
        m_y[ant_id] = new_pos_ant.y;

        if ( new_pos_ant == pos_nest ) {
            if ( m_is_loaded[ant_id] ) {
                cpteur_food += 1;
            }
            m_is_loaded[ant_id] = false;
        }
        if ( new_pos_ant == pos_food ) {
            m_is_loaded[ant_id] = true;
        }
    }
}
