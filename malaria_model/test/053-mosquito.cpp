#include "third_party/catch2/catch.hpp"
#include <iostream>
// #include <limits> //std::numeric_limits

#include <iomanip>

#include <algorithm>

// #include <vector>

#include <cmath>

#include "village/mosquito.h"
// #include "common/common.h"

// #include "util/util.h"
// #include "util/randomness.h"

// #include "common_test.h"


TEST_CASE( "53.1.1: MosquitoManager : step", "[mosquito::mm::step]" ) {
    // std::vector<double> max_larval_capacity_list{1000, 1000, 1000};
    // std::vector<double> max_larval_capacity_list{1000, 1500, 3700};
    // std::vector<double> max_larval_capacity_list{1000};
    const int kNum_villages = 1;
    const int kPop_per_village = 200;
    // const double kFeed_rate = 0.183;
    const double kFeed_rate = 0.14;
    const float kPrevalence = 0.10;

    int* population_of = new int[kNum_villages];
    for (int ii = 0; ii < kNum_villages; ii++) {
        population_of[ii] = kPop_per_village;
    }

    std::vector<double> max_larval_capacity_list(kNum_villages,1000);

    village::MosquitoManager mm(
        population_of,
        0.1, 2.0, 7.0, 15.0, 2.0, 20, 14,
        kFeed_rate,
        1.0/kPop_per_village,
        max_larval_capacity_list
    );
    // std::vector<float> infection_prob_list{0.1, 0.2, 0.3};
    // std::vector<float> infection_prob_list{0.1, 0.45, 0.07};
    // std::vector<float> infection_prob_list{0.1};
    // std::vector<float> infection_prob_list(1000,0.1);
    std::vector<float> infection_prob_list(kNum_villages,kPrevalence);


    const int kTmax = 365*3;

    for (int tt = 0; tt < kTmax; tt++) {
        mm.step(infection_prob_list);
        
        std::cout << tt << ":";
        
        for (uint vv = 0; vv < max_larval_capacity_list.size(); vv++) {
            // std::cout << "\t" << mm.get_num_infectious_mosq(vv);
            for (int pp = 0; pp < 9; pp++) {
                std::cout << "\t" << mm.get_mosq_pop(vv,pp);
            }
        }
        std::cout << "\n";
    }

    delete[] population_of;

}

TEST_CASE( "53.1.2: MosquitoManager", "[mosquito::mm::init_for_equilibrium]") {
    const double kFeed_rate = 0.178;

    const int kNum_villages = 3;
    int* population_of = new int[kNum_villages];
    for (int ii = 0; ii < kNum_villages; ii++) {
        population_of[ii] = 200;
    }


    std::vector<double> max_larval_capacity_list{1000, 1500, 3700};
    // std::vector<double> max_larval_capacity_list{1000, 1000};

    REQUIRE(max_larval_capacity_list.size() == kNum_villages);

    village::MosquitoManager mm(
        population_of,
        0.1, 2.0, 7.0, 15.0, 2.0, 20, 14,
        kFeed_rate,
        1.0,
        max_larval_capacity_list
    );
    std::vector<float> infection_prob_list{0.1, 0.45, 0.07};
    // std::vector<float> infection_prob_list{0.1, 0.1};
    mm.set_verbose_on();
    // mm.set_verbose_off();
    bool equilibrium_reached = mm.init_step_until_equilibrium(infection_prob_list);
    REQUIRE(equilibrium_reached);

    delete[] population_of;

}

// TEST_CASE( "53.2: Mosquito / ODE / Eckhoff", "[mosquito::eckhoff::step_for_one_year]") {
//     const double kEpsilon = 0.001;

//     const double kMax_larval_capacity = 1000;
//     const double kInfection_prob = 0.1;
//     const int kNum_years = 10;

//     const double kFeed_rate = 0.178;

//     village::MosquitoEckhoff m(
//         kFeed_rate,
//         kMax_larval_capacity
//     );

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
