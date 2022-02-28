
#ifdef __NVCC__
#include "third_party/catch2/catch.hpp"
#include <iostream>
// #include <limits> //std::numeric_limits

#include <iomanip>

#include <algorithm>

// #include <vector>

#include <cmath>

// #include "village/mosquito.h"
// #include "common/common.h"

// #include "util/util.h"
// #include "util/randomness.h"

// #include "common_test.h"

#include "village/mosquito_ode.cuh"



TEST_CASE( "53.2.1: ", "[mosquito_ode_cu]" ) {
    // std::vector<double> max_larval_capacity_list{1000, 1000, 1000};
    // std::vector<double> max_larval_capacity_list{1000, 1500, 3700};
    // std::vector<float> max_larval_capacity_list{1000};

    const int kNum_villages = 10000;


    std::vector<float> max_larval_capacity_list(kNum_villages,1000);
    // village::MosquitoEckhoff_multi mm(
    //     max_larval_capacity_list
    // );
    // std::vector<float> infection_prob_list{0.1, 0.2, 0.3};
    // std::vector<float> infection_prob_list{0.1, 0.45, 0.07};
    // std::vector<float> infection_prob_list{0.1};
    std::vector<float> infection_prob_list(kNum_villages,0.1);

    // MosquitoEckhoff_multi::stepper_type  stepper(9*max_larval_capacity_list.size());

    village::MosquitoManagerODEThrust mm(
        max_larval_capacity_list
    );


    const int kTmax = 500;

    for (int tt = 0; tt < kTmax; tt++) {

        mm.step(infection_prob_list);

        // mm.pre_integration_update(infection_prob_list);
        
        // odeint::integrate_const(village::MosquitoEckhoff_multi::stepper_type(), mm, mm.get_population(), tt+0.0, tt+1.0, 0.01);

        // std::cout << tt << ":";
        


        // for (uint vv = 0; vv < max_larval_capacity_list.size(); vv++) {
        //     std::cout << "\t" << mm.get_num_infectious_mosq(vv);
        //     for (int pp = 0; pp < 9; pp++) {
        //         std::cout << "\t" << mm.get_mosq_pop(vv,pp);
        //     }
        // }
        // std::cout << tt << "\n";
    }

}

// TEST_CASE( "53.1.2: MosquitoManager", "[mosquito::mm::init_for_equilibrium]") {
//     std::vector<double> max_larval_capacity_list{1000, 1500, 3700};
//     // std::vector<double> max_larval_capacity_list{1000, 1000};

//     village::MosquitoManager mm(
//         max_larval_capacity_list
//     );
//     std::vector<float> infection_prob_list{0.1, 0.45, 0.07};
//     // std::vector<float> infection_prob_list{0.1, 0.1};
//     mm.set_verbose_on();
//     // mm.set_verbose_off();
//     bool equilibrium_reached = mm.init_step_until_equilibrium(infection_prob_list);
//     REQUIRE(equilibrium_reached);
// }

// TEST_CASE( "53.2: Mosquito / ODE / Eckhoff", "[mosquito::eckhoff::step_for_one_year]") {
//     const double kEpsilon = 0.001;

//     const double kMax_larval_capacity = 1000;
//     const double kInfection_prob = 0.1;
//     const int kNum_years = 10;

//     village::MosquitoEckhoff m(kMax_larval_capacity);

//     std::cout << "Average daily number of infectious mosquitoes:\n";
    
//     double current_year = 0.0;
//     double last_year = 0.0;
//     for (int ii = 0; ii < kNum_years; ii++) {
//         current_year = m.step_for_one_year(kInfection_prob);
//         std::cout << "Year " << ii << ": ";
//         std::cout << current_year << "\n";
//         if (ii > 1) { // approx. equilibrium after two years
//             REQUIRE(std::abs(last_year - current_year) < kEpsilon);
//         }
//         last_year = current_year;
//     }

// }
#endif
