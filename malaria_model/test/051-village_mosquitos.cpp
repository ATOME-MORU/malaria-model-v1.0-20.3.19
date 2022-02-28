#include "third_party/catch2/catch.hpp"
#include <iostream>
#include <limits> //std::numeric_limits

#include <iomanip>

#include <algorithm>

// #include <vector>

#include <cmath>

#include "village/village.h"
#include "common/common.h"

// #include "util/util.h"
#include "util/randomness.h"

// #include "common_test.h"


TEST_CASE( "51.1: Mosquito Functions: mosquito_incubation_kernel", "[mosquito::incubation]" ) {

    const uint test_rounds = 50;
    const uint test_rounds_printed = 10;
    const float target_incubation_probability = 1.0/14.0;
    const float test_tolerence = 0.05;

    const int base_parasite_count = 150;

    float init_total_infected = 0;
    float init_total_infectious = 0;

    float result_total_infected = 0;
    float result_total_infectious = 0;

    float result_rate = 0;


    const uint num_villages = 4000;
    const uint num_parasite_types = 4;

    int** infected_array = new int* [num_villages];
    int** infectious_array = new int* [num_villages];

    for (uint vv = 0; vv < num_villages; vv++) {
        infected_array[vv] = new int[num_parasite_types];
        infectious_array[vv] = new int[num_parasite_types];

        infected_array[vv][0] = 0;
        infectious_array[vv][0] = 0;
        for (uint pp = 1; pp < num_parasite_types; pp++) {
            float rnd_nmbr = 0;

            util::set_random_numbers_uniform(&rnd_nmbr, 1);
            infected_array[vv][pp] = static_cast<int>(rnd_nmbr * base_parasite_count);
            // infected_array[vv][pp] = static_cast<int>(base_parasite_count);
            init_total_infected += infected_array[vv][pp];

            util::set_random_numbers_uniform(&rnd_nmbr, 1);
            infectious_array[vv][pp] = static_cast<int>(rnd_nmbr * base_parasite_count);
            // infectious_array[vv][pp] = static_cast<int>(base_parasite_count);
            init_total_infectious += infectious_array[vv][pp];
        }
    }

    int init_vll0_type0_infected = infected_array[0][0];
    int init_vll0_type0_infectious = infectious_array[0][0];

    std::cout << "Running " << test_rounds << " time steps, printing first " << test_rounds_printed << ":\n";

    int log_num_incubations = 0;

    for (uint tt = 0; tt < test_rounds; tt++ ){
        village::VillageManager::mosquito_incubation_step_kernel(
            infected_array,
            infectious_array,
            num_villages,
            num_parasite_types,
            target_incubation_probability,
            log_num_incubations
        );

        result_total_infected = 0;
        result_total_infectious = 0;
        result_rate = 0;

        for (uint vv = 0; vv < num_villages; vv++) {
            for (uint pp = 0; pp < num_parasite_types; pp++) {

                result_total_infected += infected_array[vv][pp];
                result_total_infectious += infectious_array[vv][pp];
            }
        }

        // overall result
        REQUIRE((init_total_infected + init_total_infectious) == (result_total_infected + result_total_infectious));
        REQUIRE((init_total_infected - result_total_infected) == (result_total_infectious - init_total_infectious));
        REQUIRE(init_total_infected >= result_total_infected);
        REQUIRE(init_total_infectious <= result_total_infectious);

        result_rate = float(init_total_infected - result_total_infected) / init_total_infected;

        if (tt < test_rounds_printed) {

            std::cout << "[step #" << tt << "]";

            std::cout << std::setprecision(4) << std::fixed;

            std::cout << "Number of incubations: expected (binomial mean) " << init_total_infected * target_incubation_probability;
            std::cout << ", actual " << init_total_infected - result_total_infected
                        << "(" << float(std::abs(init_total_infected - result_total_infected) - init_total_infected * target_incubation_probability)
                                    /( init_total_infected * target_incubation_probability)*100 << "%)";

            std::cout << "; Result incubation rate: " << result_rate << ", "
                        << float(result_rate-target_incubation_probability)/target_incubation_probability*100 << "%"
                        << " (target: " << target_incubation_probability
                        << " +/- " << std::setprecision(1) << test_tolerence*100 << "%)";
        }

        REQUIRE(std::abs(result_rate - target_incubation_probability) < target_incubation_probability * test_tolerence);

        if (tt < test_rounds_printed) {
            std::cout << " OK\n";
        }

        init_total_infected = result_total_infected;
        init_total_infectious = result_total_infectious;


        // check result of one village
        REQUIRE( (infected_array[0][0] + infectious_array[0][0]) == (init_vll0_type0_infected + init_vll0_type0_infectious) );
        REQUIRE( (init_vll0_type0_infected - infected_array[0][0]) == (infectious_array[0][0] - init_vll0_type0_infectious) );
        REQUIRE( init_vll0_type0_infected >= infected_array[0][0] );
        REQUIRE( init_vll0_type0_infectious <= infectious_array[0][0] );

        init_vll0_type0_infected = infected_array[0][0];
        init_vll0_type0_infectious = infectious_array[0][0];

    }

}