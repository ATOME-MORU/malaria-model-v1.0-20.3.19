#include "third_party/catch2/catch.hpp"
#include <iostream>
#include <limits> //std::numeric_limits

#include <algorithm>

#include <vector>

#include <cmath>

#include "village/village.h"
#include "common/common.h"

#include "util/randomness.h"

TEST_CASE( "52.1: Mobility Functions: step_population_movement_static_kernel", "[mobility:step_movement_static]" ) {
    
    const int num_villages = 3;
    int village_population[] = {100,0,10};
    int num_humans = 0;
    for (int vv = 0; vv < num_villages; vv++){
        num_humans += village_population[vv];
    }

    std::vector<std::vector<int>> villager_reg;
    villager_reg.resize(num_villages);
    int* human_home_village_array = new int[num_humans];
    common::HumanMobilityType* human_mobility_type_array = new common::HumanMobilityType[num_humans];
    int human_index = 0;

    for (int vv = 0; vv < num_villages; vv++){
        for (int hh = 0; hh < village_population[vv]; hh++) {
            villager_reg[vv].push_back(human_index);
            human_home_village_array[human_index] = vv;
            human_mobility_type_array[human_index] = common::HumanMobilityType::kStatic;
            human_index++;
        }
    }

    std::vector<std::vector<int>> mobility_network_static;
    mobility_network_static.resize(3);
    mobility_network_static[0].push_back(1);
    mobility_network_static[1].push_back(2);
    mobility_network_static[2].push_back(0);


    float static_population_move_out_rate = 1;
    float static_population_return_home_rate = 1;
    float non_static_population_move_out_rate = 0;
    float non_static_population_return_home_rate = 0;

    float* random_numbers_array = new float[num_humans*2];
    util::set_random_numbers_uniform(random_numbers_array, num_humans*2);

    village::VillageManager::step_population_movement_static_kernel(
        num_villages,
        villager_reg,
        mobility_network_static,
        num_humans,
        human_mobility_type_array,
        human_home_village_array,
        static_population_move_out_rate,    // G_mobility_static
        static_population_return_home_rate, // G_static_return
        non_static_population_move_out_rate,   // G_mobility_mobile
        non_static_population_return_home_rate, // G_mobile_return

        random_numbers_array,
        num_humans*2
    );

    // std::cout << villager_reg[0].size() << "\n";
    // std::cout << villager_reg[1].size() << "\n";
    // std::cout << villager_reg[2].size() << "\n";

    REQUIRE(villager_reg[0].size()==village_population[2]);
    REQUIRE(villager_reg[1].size()==village_population[0]);
    REQUIRE(villager_reg[2].size()==village_population[1]);

    village::VillageManager::step_population_movement_static_kernel(
        num_villages,
        villager_reg,
        mobility_network_static,
        num_humans,
        human_mobility_type_array,
        human_home_village_array,
        static_population_move_out_rate,    // G_mobility_static
        static_population_return_home_rate, // G_static_return
        non_static_population_move_out_rate,   // G_mobility_mobile
        non_static_population_return_home_rate, // G_mobile_return

        random_numbers_array,
        num_humans*2
    );

    REQUIRE(villager_reg[0].size()==village_population[0]);
    REQUIRE(villager_reg[1].size()==village_population[1]);
    REQUIRE(villager_reg[2].size()==village_population[2]);

    delete[] human_home_village_array;
    delete[] human_mobility_type_array;

}