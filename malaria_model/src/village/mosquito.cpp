#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include <string>

// #include <boost/random/mersenne_twister.hpp> // boost::mt19937
// #include <boost/random/normal_distribution.hpp> // boost::random::uniform_real_distribution
// #include <boost/proto/functional/std/utility.hpp> //make_pair
#include <utility>      // std::pair
// #include <boost/array.hpp>


#include "village/mosquito.h"

#include "common/common.h"

#include "village/village.h"
#include "util/randomness.h"
#include "util/graph.h"


namespace village {

MosquitoEckhoff::Fx::Fx(
        double aquatic_mortality,
        double maturation_days,
        double immat_life_days,
        double adult_life_days,
        double laying_days,
        double eggs_per_lay,
        double incubation_days,
        double feed_rate,
        double max_larval_capacity,
        double infection_prob
    ) :
        kAquatic_mortality(aquatic_mortality),
        kImmature_adult(1.0/maturation_days),
        kIadult_life(1.0/immat_life_days),
        kAdult_life(1.0/adult_life_days),
        kLay_rate(1.0/laying_days),
        kOviposition(eggs_per_lay),
        kInfected_to_infectious_rate(1.0/incubation_days),
        kFeed_rate(feed_rate),
        kMax_larval_capacity(max_larval_capacity),
        kInfection_prob(infection_prob) {

}


void MosquitoEckhoff::Fx::operator() (
            MosquitoEckhoff::ode_state_t &y,
            MosquitoEckhoff::ode_state_t &dydt,
            MosquitoEckhoff::ode_time_t t
        ) const {
    // float Tk_stoch = util::get_rand_normal(298.15, 7.0);

    // double Tk = 273.15+25.0+5.0*cos(t/kNum_days_per_year*2.0*kPi);
    // double Tk = 273.15+25.0+0.8*5.0*cos((t-60)/kNum_days_per_year*2.0*kPi);
    // double Tk = 273.15+25.0+0.8*5.0*cos((t+kNum_days_per_year-60)/kNum_days_per_year*2.0*kPi);

    // double Tk = 273.15+25.0+0.8*5.0*cos((t+2*kNum_days_per_year)/kNum_days_per_year*2.0*kPi);
    // double Tk = 273.15+25.0+5.0*cos((t+2*kNum_days_per_year)/kNum_days_per_year*2.0*kPi);
    // double Tk = 273.15+25.0+5.0*cos(t/kNum_days_per_year*2.0*kPi); 
    // double Tk = 273.15+25.0+0.8*5.0*cos((static_cast<float>(t)-70)/kNum_days_per_year*2.0*kPi);
    
    double Tk = 273.15+25.0+0.8*5.0*cos(t/common::kDays_per_year*2.0*common::kPI);
    // double Tk = 273.15+25.0+0.8*5.0*cos((t)/kNum_days_per_year*2.0*kPi);
    // double Tk = Tk_stoch+0.8*5.0*cos((t)/kNum_days_per_year*2.0*kPi);

    // if ((int(t)%10) == 0) {
    //     std::cout << "ODE t=" << t << ", Tk=" << Tk << "\n";
    //     std::cin.ignore();
    // }


    //the expressions below are based on temperature input in Kelvin
    double larva_to_immature = kAquatic_arrhenius_1*exp(-kAquatic_arrhenius_2/Tk);

    // std::cout << "larva_to_immature=" << larva_to_immature << "\n";
    // double infected_to_infectious = kInfected_arrhenius_1*exp(-kInfected_arrhenius_2/Tk);
    // double infected_to_infectious = kInfected_arrhenius_1*exp(-kInfected_arrhenius_2/Tk)-0.02;
    double infected_to_infectious = kInfected_to_infectious_rate;
    double eggs_to_larva = kEgg_arrhenius1*exp(-kEgg_arrhenius2/Tk);

    //Eggs
    dydt[0] = -y[0]*(eggs_to_larva+kAquatic_mortality)
                + kLay_rate*kOviposition*(y[4]+kInfected_egg_lay*(y[6]+y[8]));
    //Larva
    dydt[1] = y[0]*eggs_to_larva*(1.0-y[1]/kMax_larval_capacity)
                - y[1]*(kAquatic_mortality+larva_to_immature);
    //Immature
    dydt[2] = y[1]*larva_to_immature
                - y[2]*(kImmature_adult+kIadult_life);
                // - y[2]*(kImmature_adult+1.0/7.0);
                // - y[2]*(kImmature_adult+1.0/5.0);
                // - y[2]*(kImmature_adult+kAdult_life);

    //Adult uninfected biting
    dydt[3] = y[2]*kImmature_adult 
                - y[3]*(kAdult_life+kFeed_rate)
                + y[4]*kLay_rate;
    //Adult uninfected laying eggs
    dydt[4] = y[3]*kFeed_rate*(1.0-kInfection_prob)*(1.0-kHuman_feeding_mort)
                - y[4]*(kAdult_life+kLay_rate); 

    //Adult infected - biting
    dydt[5] = -y[5]*(kAdult_life+infected_to_infectious+kFeed_rate)
                + y[6]*kLay_rate;
    //Adult infected - laying eggs
    dydt[6] = -y[6]*(kAdult_life+kLay_rate+infected_to_infectious)
                + y[5]*kFeed_rate*(1.0-kHuman_feeding_mort)
                + y[3]*kFeed_rate*kInfection_prob*(1.0-kHuman_feeding_mort);

    //Adult infectious biting
    dydt[7] = y[5]*infected_to_infectious
                - y[7]*(kAdult_life+kFeed_rate)
                + y[8]*kLay_rate;
    //Adult infectious laying eggs
    dydt[8] = y[6]*infected_to_infectious
                - y[8]*(kAdult_life+kLay_rate)
                + y[7]*kFeed_rate*(1.0-kHuman_feeding_mort);
}


MosquitoEckhoff::MosquitoEckhoff(
        double aquatic_mortality,
        double maturation_days,
        double immat_life_days,
        double adult_life_days,
        double laying_days,
        double eggs_per_lay,
        double incubation_days,
        double feed_rate,
        double max_larval_capacity 
    ) :
        kAquatic_mortality(aquatic_mortality),
        kMaturation_days(maturation_days),
        kImmat_life_days(immat_life_days),
        kAdult_life_days(adult_life_days),
        kLaying_days(laying_days),
        kEggs_per_lay(eggs_per_lay),
        kIncubation_days(incubation_days),
        kFeed_rate(feed_rate),
        kMax_larval_capacity(max_larval_capacity) {

    // this->population.push_back(0.0);   // Eggs
    // this->population.push_back(0.0);   // Larva
    // this->population.push_back(200.0); // Immature

    this->population.push_back(150);   // Eggs
    this->population.push_back(150);   // Larva
    this->population.push_back(80); // Immature

    this->population.push_back(0);   // Uninfected biting
    this->population.push_back(0);   // Uninfected laying

    this->population.push_back(0);   // Infected biting
    this->population.push_back(0);   // Infected laying

    this->population.push_back(0);   // Infectious biting
    this->population.push_back(0);   // Infectious laying

    this->reset_time();
}

MosquitoEckhoff::~MosquitoEckhoff(){

}

void MosquitoEckhoff::reset_time() {
    this->current_time = 0.0;
}


double MosquitoEckhoff::step(double infection_prob) {

    boost::mt19937 rng;

// typedef odeint::controlled_runge_kutta< odeint::runge_kutta_dopri5< ode_state_t > > dopri_stepper_type;
// typedef odeint::dense_output_runge_kutta< dopri_stepper_type > dense_stepper_type;


    // // odeint::integrate_adaptive(
    // // int num_steps = odeint::integrate_const(
    odeint::integrate_const(
        // this->ode_stepper,
        // dense_stepper_type(),
        odeint::runge_kutta4<ode_state_t>(),
        // odeint::euler<ode_state_t>(),
        Fx(
            this->kAquatic_mortality,
            this->kMaturation_days,
            this->kImmat_life_days,
            this->kAdult_life_days,
            this->kLaying_days,
            this->kEggs_per_lay,
            this->kIncubation_days,
            this->kFeed_rate,
            this->kMax_larval_capacity,
            infection_prob
        ),
        this->population,
        this->current_time,
        this->current_time + this->kTime_length_per_integration,
        this->kDt
    );

    // std::cout << "kLay_rate: " << this->kLay_rate << "\n";
    // std::cin.ignore();

    // odeint::integrate_const(
    //     // this->ode_stepper,
    //     // dense_stepper_type(),
    //     stochastic_euler(),
    //     std::make_pair(
    //         Fx(
    //             this->kFeed_rate,
    //             this->kMax_larval_capacity,
    //             infection_prob
    //         ),
    //         Fx_stoch(
    //             rng,
    //             0.0
    //         )
    //     ),
    //     this->population,
    //     this->current_time,
    //     this->current_time + this->kTime_length_per_integration,
    //     this->kDt
    // );


    this->current_time += this->kTime_length_per_integration;

    return this->population[7]; //number of adult infectious biting mosquitos
}

double MosquitoEckhoff::step_for_one_year(double infection_prob) {
    double avg_num_infectious_mosq_per_day = 0.0;
    // for (int ii = 0; ii < Fx::kNum_days_per_year; ii++) {
    for (int ii = 0; ii < common::kDays_per_year; ii++) {
        avg_num_infectious_mosq_per_day += this->step(infection_prob);
    }
    // avg_num_infectious_mosq_per_day /= Fx::kNum_days_per_year;
    avg_num_infectious_mosq_per_day /= common::kDays_per_year;
    return avg_num_infectious_mosq_per_day;
}


MosquitoManager::MosquitoManager(
        const int* population_of,
        float aquatic_mortality,
        float maturation_days,
        float immat_life_days,
        float adult_life_days,
        float laying_days,
        float eggs_per_lay,
        float incubation_days,
        double feed_rate,
        double carrying_capacity_multiplier,
        int num_mosquitos,
        std::string max_larval_capacity_file_name
    ) :
        MosquitoManager(
            population_of,
            aquatic_mortality,
            maturation_days,
            immat_life_days,
            adult_life_days,
            laying_days,
            eggs_per_lay,
            incubation_days,
            feed_rate,
            carrying_capacity_multiplier,
            this->init_construct_max_larval_capacity_list_from_file(
                num_mosquitos,
                max_larval_capacity_file_name
            )
        ) {

}

MosquitoManager::MosquitoManager(
        const int* population_of,
        float aquatic_mortality,
        float maturation_days,
        float immat_life_days,
        float adult_life_days,
        float laying_days,
        float eggs_per_lay,
        float incubation_days,
        double feed_rate,
        double carrying_capacity_multiplier,
        std::vector<double> max_larval_capacity_list
    ) {

    double tot_larval_capacity = 0.0;

    for (uint vv = 0; vv < max_larval_capacity_list.size(); vv++) {

        this->mosquito_of.push_back(
            new MosquitoEckhoff(
                aquatic_mortality,
                maturation_days,
                immat_life_days,
                adult_life_days,
                laying_days,
                eggs_per_lay,
                incubation_days,
                feed_rate,
                max_larval_capacity_list[vv] * carrying_capacity_multiplier * population_of[vv]
            )
        );
        tot_larval_capacity += max_larval_capacity_list[vv] * carrying_capacity_multiplier * population_of[vv];
        this->num_infectious_mosq_of.push_back(-1.0);

    }

    std::cout << "ODE initialised tot_larval_capacity=" << tot_larval_capacity << "\n";

}


std::vector<double> MosquitoManager::init_construct_max_larval_capacity_list_from_file(
        int num_mosquitos,
        std::string max_larval_capacity_file_name
    ) {

    std::ifstream input_file_stream(max_larval_capacity_file_name);
    assert(input_file_stream);

    std::string str_line;

    std::vector<double> max_larval_capacity_list;

    for (int ii = 0; ii < num_mosquitos; ii++) {
        if (!std::getline(input_file_stream, str_line).eof()) {
            max_larval_capacity_list.push_back(std::stod(str_line));
        } else {
            std::cout << "Need " << num_mosquitos << " mosquitos. "
                        << "Only " << ii
                        << " found in " << max_larval_capacity_file_name
                        << ".\n";
            break;
        }
    }

    assert( static_cast<int>(max_larval_capacity_list.size()) == num_mosquitos);

    return max_larval_capacity_list;
}

MosquitoManager::~MosquitoManager() {
    for (uint vv = 0; vv < this->mosquito_of.size(); vv++) {
        delete this->mosquito_of[vv];
    }
}

void MosquitoManager::step(const std::vector<float>& infection_prob_list) {
    this->sum_num_larva_capacity = 0;
    for (uint vv = 0; vv < this->mosquito_of.size(); vv++) {
        this->num_infectious_mosq_of[vv] = 
            this->mosquito_of[vv]->step(static_cast<double>(infection_prob_list[vv]));
            if (this->mosquito_of[vv]->if_capacity_reached()){
                this->sum_num_larva_capacity++;
            }
        // if (this->num_infectious_mosq_of[vv] < 0.001) {
        //     std::cout << "vv:" << vv << ", num_infectious_mosq_of[vv]=" << num_infectious_mosq_of[vv] << std::endl;
        //     std::cin.ignore(); 
        // }
        // if (!std::isnormal(this->num_infectious_mosq_of[vv])) { // #include <cmath>
        //     this->num_infectious_mosq_of[vv] = 0.0;
        // }
    }
}


bool MosquitoManager::init_step_until_equilibrium(const std::vector<float>& infection_prob_list) {

    assert(infection_prob_list.size() == this->mosquito_of.size());

    std::vector<double> current_year(infection_prob_list.size(), 0.0);
    std::vector<double> last_year(infection_prob_list.size(), 0.0);

    if (this->verbose) {
        std::cout << "MM.init_step_until_equilibrium\n"
                    << "\tkEquilibrium_epsilon=" << this->kEquilibrium_epsilon << "\n"
                    << "\tkEquilibrium_max_num_years=" << this->kEquilibrium_max_num_years
                    << "\n";
    }

    bool equilibrium_reached = false;
    int yy = 0;
    do {

        equilibrium_reached = true;
        for (uint vv = 0; vv < this->mosquito_of.size(); vv++) {
            current_year[vv] = 
                this->mosquito_of[vv]->step_for_one_year(static_cast<double>(infection_prob_list[vv]));

            if (std::abs(current_year[vv] - last_year[vv]) > this->kEquilibrium_epsilon) {
                equilibrium_reached = false;
            }
            last_year[vv] = current_year[vv];
        }

        if (this->verbose) {
            std::cout << "Year # " << yy << "\n";
        }
 
        yy++;

    } while ((!equilibrium_reached) && (yy <= this->kEquilibrium_max_num_years));

    if (this->verbose) {
        if (equilibrium_reached) {
            std::cout << "Equilibrium reached after year " << yy-1 <<".\n";

            std::cout << current_year[0] << "\n";
        } else {
            std::cout << "Failed to reach equilibrium.\n";
        }
    }

    return equilibrium_reached;

}




}