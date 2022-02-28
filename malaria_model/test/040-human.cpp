#include <iostream>
#include <map>

#include <cmath>
#include <string>

#include "third_party/catch2/catch.hpp"
#include "common/common.h"
// #include "human/blood_system.h"
#include "human/human.h"

// #include "util/util.h"
#include "util/randomness.h"

// #include "common_test.h"

// TEST_CASE( "4.1: HumanManager: Age Distribution)", "[human:age]" ) {

//     int vector_size = 5;
//     float density_vector[5] = {100.0, 200.0, 300.0, 400.0, 500.0};
//     common::ParasiteType type_vector[5] = {common::ParasiteType::kRa,
//                                            common::ParasiteType::kRa,
//                                            common::ParasiteType::kRb,
//                                            common::ParasiteType::kRb,
//                                            common::ParasiteType::kR0};
//     int num_tests = 900000;
//     common::ParasiteType dominant_type_selected[900000] = {common::ParasiteType::kNoParasite};
//     std::map<common::ParasiteType, float> type_accumulator;



//     // REQUIRE(type_accumulator[common::ParasiteType::kRb]
//     //             /type_accumulator[common::ParasiteType::kR0]
//     //         - (density_vector[2]+density_vector[3])
//     //             /density_vector[4] < 0.01);

//     // human::HumanManager::print_distribution<common::ParasiteType>(dominant_type_selected,num_tests);

// }

TEST_CASE( "40.2: HumanManager: init_itn_kernel)", "[human:init_itn_kernel]" ) {

    const uint data_size = 900000;
    const float if_itn_probability = 0.2;
    const float test_tolerence = 0.001;

    // float* random_numbers = new float[data_size];
    float random_numbers[data_size] = {0.0};
    util::set_random_numbers_uniform(random_numbers, data_size);

    bool if_itn[data_size] = {false};

    human::HumanManager::init_itn_kernel(
            if_itn_probability,
            if_itn,
            random_numbers,
            data_size
        );

    std::map<bool, int> if_itn_accumulator;
    for (uint ii = 0; ii < data_size; ii++){
        if_itn_accumulator[if_itn[ii]]++;
    }

    float result_itn_probability = if_itn_accumulator[true] / (float) data_size;
    std::cout << "result ITN probability: " << result_itn_probability
                << " (target:" << if_itn_probability << " +/- " << test_tolerence << ")"<< std::endl;
    REQUIRE(std::abs(result_itn_probability - if_itn_probability) < test_tolerence);

    // REQUIRE(type_accumulator[common::ParasiteType::kRb]
    //             /type_accumulator[common::ParasiteType::kR0]
    //         - (density_vector[2]+density_vector[3])
    //             /density_vector[4] < 0.01);

    human::HumanManager::print_distribution<bool>(if_itn,(int)data_size);

}

TEST_CASE( "40.3: HumanManager: func_death_probability)", "[human:func_death_probability]" ) {
    const int kMax_age = 100;
    const int kMax_star = 20;

    std::cout << "Kernal 0: (karenjan17.cpp)\n";
    for (int ii = 0; ii < kMax_age; ii+=10) {
        float prob = human::HumanManager::func_death_probability_kernel_0(ii);
        std::cout << ii << " " << prob << " :"
                << std::string(prob*10000*kMax_star, '*')
                 << std::endl;
    }

    std::cout << "Kernal 1_0: (durgrestest.cpp / aguasSIR.h)\n";
    for (int ii = 0; ii < kMax_age; ii+=10) {
        float prob = human::HumanManager::func_death_probability_kernel_1_0(ii);
        std::cout << ii << " " << prob << " :"
                << std::string(prob*10000*kMax_star, '*')
                 << std::endl;
    }

    std::cout << "Kernal 1_1: (mmc_wp2, durgrestest.cpp / aguasSIR.h)\n";
    for (int ii = 0; ii < kMax_age; ii+=10) {
        float prob = human::HumanManager::func_death_probability_kernel_1_1(ii);
        std::cout << ii << " " << prob << " :"
                << std::string(prob*10000*kMax_star, '*')
                 << std::endl;
    }

    std::cout << "Kernal 2: (mda_coverage, 03 Dec 2018) 0.01 scaled in number of *s\n";
    for (int ii = 0; ii < kMax_age; ii+=10) {
        float prob = human::HumanManager::func_death_probability_kernel_2(ii);
        std::cout << ii << " " << prob << " :"
                << std::string(prob*100*kMax_star, '*')
                 << std::endl;
    }

}