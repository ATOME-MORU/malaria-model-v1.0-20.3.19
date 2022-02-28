#include <iostream>
#include <map>

#include <algorithm>

#include <cmath>

#include "third_party/catch2/catch.hpp"
#include "util/util.h"
#include "util/randomness.h"
#include "common/common.h"
#include "human/human.h"
#include "human/blood_system.h"

TEST_CASE( "30.1: Blood_System functions: parasite_dominance_density_based()", "[blood_system:dominance_density_based]" ) {

    const int vector_size = 5;
    float density_vector[vector_size] = {100.0, 200.0, 300.0, 400.0, 500.0};
    common::ParasiteType type_vector[5] = {common::ParasiteType::kRa,
                                           common::ParasiteType::kRa,
                                           common::ParasiteType::kRb,
                                           common::ParasiteType::kRb,
                                           common::ParasiteType::kR0};
    const int num_tests = 900000;
    common::ParasiteType dominant_type_selected[num_tests] = {common::ParasiteType::kNoParasite};
    std::map<common::ParasiteType, float> type_accumulator;

    std::cout << "Testing BloodSystem::parasite_dominance_density_based() with sample size("
                << num_tests << "). This could take a while..." << std::flush;
    for (int tt = 0; tt < num_tests; tt++) {
        dominant_type_selected[tt] = 
            human::BloodSystemManager::parasite_dominance_density_based<float>(
                0,vector_size-1,density_vector,type_vector
                );
            
        type_accumulator[dominant_type_selected[tt]]++;
    }

    std::cout << "done.\nSample data:\n";

    for (int ii = 0; ii < vector_size; ii++){
        std::cout << type_vector[ii] << ":" << density_vector[ii] << std::endl;
    }

    std::cout << "Result:\n";

    REQUIRE(std::abs(
            type_accumulator[common::ParasiteType::kRa]
                /type_accumulator[common::ParasiteType::kRb]
            - (density_vector[0]+density_vector[1])
                /(density_vector[2]+density_vector[3])
            ) < 0.01);

    REQUIRE(std::abs(
            type_accumulator[common::ParasiteType::kRb]
                /type_accumulator[common::ParasiteType::kR0]
            - (density_vector[2]+density_vector[3])
                /density_vector[4]
            ) < 0.01);

    human::HumanManager::print_distribution<common::ParasiteType>(dominant_type_selected,num_tests);

}

TEST_CASE("30.2 Blood_System functions: clinical_resolution_kernel", "[blood_system:clinical_resolution_kernel]") {
    const uint array_size = 900000;
    const float clinical_percentage = 0.7;
    const float given_mellow_rate = 0.2;
    const float test_tolerence = 0.001;

    bool clinical_status_array[array_size] = {};

    std::fill(clinical_status_array, clinical_status_array+array_size, false);
    std::fill(clinical_status_array, clinical_status_array+int(array_size*clinical_percentage), true);

    float target_asymptomatic_percentage = 1 - clinical_percentage + given_mellow_rate * clinical_percentage;

    human::BloodSystemManager::clinical_resolution_kernel(
            clinical_status_array,
            array_size,
            given_mellow_rate
        );
    int clinical_count = 0;
    for (uint ii = 0; ii < array_size; ii++){
        if (clinical_status_array[ii]) clinical_count++;
    }

    float result_percentage = float(array_size-clinical_count)/array_size;
    std::cout << "Result asymptomatic percentage:" << result_percentage
                << "(" << array_size-clinical_count << "/" << array_size
                <<"), target percentage: " << target_asymptomatic_percentage
                << " +/- " << test_tolerence << std::endl;
    REQUIRE(std::abs(result_percentage - target_asymptomatic_percentage ) < test_tolerence );
}

TEST_CASE("30.3 Blood_System drug_failure rate calculation", "[blood_system:drug_failure_calibration]") {

    const int kNum_hosts = 100000;
    const int kDrug_length = 28; // days

    const common::ParasiteType kParasite_type_to_test = common::ParasiteType::kR0;
    const common::DrugName kDrug_to_test = common::DrugName::kPiparte;

    float* random_numbers_drug_loss = new float[kNum_hosts];
    int random_numbers_drug_loss_size = kNum_hosts;
    float rnd_nmb = 0.0;

    const int kNum_tests = 100;
    std::vector<int> test_results(kNum_tests, 0.0);

    for (int test_index = 0; test_index < kNum_tests; test_index++){

        std::cout << "run#" << test_index << "\n";
        human::BloodSystemManager bsys_mgr(
            kNum_hosts,
            
            0.0, 0.0, 0.0,
            
            true, true,

            0.0,
            
            0.0, 0.0,
            
            0.0, 0.0,
            
            0.0,
            
            0, //dominance function

            2,0.0,0.0, //clinical probability function
            
            10.0, // mellow days

            5.0, 2.0, 15.0, 5.0, 160.0, 50.0, // stage progression

            // 160.0, // natural recovery days
            0.0, // resistance cost

            28, // drug period
            false, false, false,
            1.0, // base effect multiplier
            1.0,1.0,1.0,1.0,1.0,1.0,1.0,

            0.0,
            0.0,
            0.0,
            true,true,

            0, 6,

            0, 0, 0.0
        );

        bsys_mgr.init_for_drug_failure_test(kParasite_type_to_test, kDrug_to_test);

        for (int tt = 1; tt <= kDrug_length+2; tt++) {

            util::set_random_numbers_uniform(&rnd_nmb, 1);

            util::set_random_numbers_uniform(
                random_numbers_drug_loss,
                random_numbers_drug_loss_size
            );
            // if (rnd_nmb < 0.5) {
                bsys_mgr.drug_effect(tt);
                bsys_mgr.drug_loss(random_numbers_drug_loss,random_numbers_drug_loss_size);
            // } else {
            //     bsys_mgr.drug_loss(random_numbers_drug_loss);
            //     bsys_mgr.drug_effect(tt);
            // }



            bsys_mgr.step_report_drug_failure_by_clinical_status_on_day(tt);
            bsys_mgr.step_report_drug_composition();

            // std::cout << "drug_composition:";
            // for (const auto& dd: util::Enum_iterator<common::DrugName>()) {
            //     std::cout << bsys_mgr.daily_record_of_drug_current[static_cast<int>(dd)] << ";";
            // }
            // std::cout << "\n";

            // std::cout << "drug_composition_infected:";
            // for (const auto& dd: util::Enum_iterator<common::DrugName>()) {
            //     std::cout << bsys_mgr.daily_record_of_drug_current_in_infected_population[static_cast<int>(dd)] << ";";
            // }
            // std::cout << "\n";


            if (test_index == kNum_tests-1) {
                
                std::cout << "Day " << tt << ":\t" ;
                for (const auto& dd: util::Enum_iterator<common::DrugName>()) {
                    std::cout << dd << "-" << bsys_mgr.daily_record_of_drug_failure_by_effetive_period[static_cast<int>(dd)] << "; " ;
                }
                std::cout << "\n";
            }

            if (tt == kDrug_length) {
                test_results[test_index] = bsys_mgr.daily_record_of_drug_failure_by_effetive_period[static_cast<int>(kDrug_to_test)];
            }

        }

    }

    float result_average = 0.0;
    for (const auto& rr : test_results) {
        result_average += rr;
    }
    result_average /= kNum_tests;

    std::cout << "Average failure count from " << kNum_tests << " runs: " << result_average << "\n";

    delete[] random_numbers_drug_loss;

}

TEST_CASE("30.3 Blood_System functions: parasite_dominance_mmc_wp2_kernel", "[blood_system:parasite_dominance_mmc_wp2_kernel]") {
    const int kNum_systems = 40;
    const int kNum_pgs_per_system = 3;
    const int kNum_pgs = kNum_systems * kNum_pgs_per_system;

    int kHead_index[kNum_systems] = {};
    int kTail_index[kNum_systems] = {};

    for (int ss = 0; ss < kNum_systems; ss++) {
        kHead_index[ss] = 0;
        kTail_index[ss] = 2;
    }

    common::StageName most_advanced_stage_array[kNum_systems] = {};
    common::ParasiteType dominant_type_array[kNum_systems] = {};
    
    const common::DrugName kDrug_ab = common::DrugName::kPiparte;
    const common::DrugName kDrug_a = common::DrugName::kArtesunate;
    const common::DrugName kDrug_b = common::DrugName::kPiperaquine;
    const common::DrugName kDrug_other = common::DrugName::kSP;
    const common::DrugName kDrug_none = common::DrugName::kNoDrug;

    const common::StageName pg_stage_array[kNum_pgs] = {
        //0
        common::StageName::kLiver, common::StageName::kLiver, common::StageName::kLiver,
        common::StageName::kLiver, common::StageName::kBlood, common::StageName::kLiver,
        common::StageName::kInfectious, common::StageName::kBlood, common::StageName::kLiver,
        common::StageName::kBlood, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,

        common::StageName::kBlood, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kBlood, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kBlood, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,

        //10
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,

        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,

        //20
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,

        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kInfectious, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kInfectious, common::StageName::kLiver, common::StageName::kLiver,
        common::StageName::kLiver, common::StageName::kInfectious, common::StageName::kLiver,

        //30
        common::StageName::kLiver, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kLiver, common::StageName::kLiver, common::StageName::kInfectious,
        common::StageName::kLiver, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kLiver, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kLiver, common::StageName::kInfectious, common::StageName::kInfectious,

        common::StageName::kLiver, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kLiver, common::StageName::kLiver, common::StageName::kInfectious,
        common::StageName::kLiver, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kLiver, common::StageName::kInfectious, common::StageName::kInfectious,
        common::StageName::kLiver, common::StageName::kInfectious, common::StageName::kInfectious

    };

    const common::ParasiteType pg_type_array[kNum_pgs] = {
        common::ParasiteType::kRa, common::ParasiteType::kRab, common::ParasiteType::kR0,
        common::ParasiteType::kR0, common::ParasiteType::kRa, common::ParasiteType::kR0,
        common::ParasiteType::kR0, common::ParasiteType::kRab, common::ParasiteType::kRa,
        common::ParasiteType::kRab, common::ParasiteType::kR0, common::ParasiteType::kR0,
        common::ParasiteType::kR0, common::ParasiteType::kR0, common::ParasiteType::kR0,

        common::ParasiteType::kRab, common::ParasiteType::kRa, common::ParasiteType::kR0,
        common::ParasiteType::kRab, common::ParasiteType::kRb, common::ParasiteType::kR0,
        common::ParasiteType::kRab, common::ParasiteType::kR0, common::ParasiteType::kR0,
        common::ParasiteType::kRab, common::ParasiteType::kR0, common::ParasiteType::kR0,
        common::ParasiteType::kR0, common::ParasiteType::kRab, common::ParasiteType::kRb,

        //10 - kRab
        common::ParasiteType::kRab, common::ParasiteType::kR0, common::ParasiteType::kRb,
        common::ParasiteType::kRab, common::ParasiteType::kR0, common::ParasiteType::kRb,
        common::ParasiteType::kRab, common::ParasiteType::kR0, common::ParasiteType::kRa,
        common::ParasiteType::kRab, common::ParasiteType::kRab, common::ParasiteType::kRab,
        common::ParasiteType::kRab, common::ParasiteType::kRab, common::ParasiteType::kRab,

        // 15 - kRa
        common::ParasiteType::kRa, common::ParasiteType::kR0, common::ParasiteType::kRb,
        common::ParasiteType::kR0, common::ParasiteType::kRa, common::ParasiteType::kR0,
        common::ParasiteType::kRb, common::ParasiteType::kR0, common::ParasiteType::kRa,
        common::ParasiteType::kRa, common::ParasiteType::kRa, common::ParasiteType::kRa,
        common::ParasiteType::kRa, common::ParasiteType::kRa, common::ParasiteType::kRa,

        // 20 - kRb
        common::ParasiteType::kR0, common::ParasiteType::kR0, common::ParasiteType::kRb,
        common::ParasiteType::kR0, common::ParasiteType::kRa, common::ParasiteType::kRb,
        common::ParasiteType::kR0, common::ParasiteType::kR0, common::ParasiteType::kRb,
        common::ParasiteType::kRb, common::ParasiteType::kRb, common::ParasiteType::kRb,
        common::ParasiteType::kRb, common::ParasiteType::kRb, common::ParasiteType::kRb,

        // 25 - kR0
        common::ParasiteType::kR0, common::ParasiteType::kR0, common::ParasiteType::kR0,
        common::ParasiteType::kR0, common::ParasiteType::kR0, common::ParasiteType::kR0,
        common::ParasiteType::kR0, common::ParasiteType::kR0, common::ParasiteType::kR0,
        common::ParasiteType::kR0, common::ParasiteType::kRb, common::ParasiteType::kRab,
        common::ParasiteType::kRb, common::ParasiteType::kR0, common::ParasiteType::kRa,

        // 30 - not the one in the liver
        common::ParasiteType::kR0, common::ParasiteType::kRa, common::ParasiteType::kRb,
        common::ParasiteType::kR0, common::ParasiteType::kR0, common::ParasiteType::kRab,
        common::ParasiteType::kRa, common::ParasiteType::kR0, common::ParasiteType::kRb,
        common::ParasiteType::kRb, common::ParasiteType::kR0, common::ParasiteType::kR0,
        common::ParasiteType::kRab, common::ParasiteType::kR0, common::ParasiteType::kRa,

        common::ParasiteType::kR0, common::ParasiteType::kRa, common::ParasiteType::kRb,
        common::ParasiteType::kR0, common::ParasiteType::kR0, common::ParasiteType::kRab,
        common::ParasiteType::kRa, common::ParasiteType::kR0, common::ParasiteType::kRb,
        common::ParasiteType::kRb, common::ParasiteType::kR0, common::ParasiteType::kR0,
        common::ParasiteType::kRab, common::ParasiteType::kR0, common::ParasiteType::kRa

    };

    const common::DrugName drug_array[kNum_systems] = {
        kDrug_none, kDrug_ab, kDrug_a, kDrug_ab, kDrug_none,

        kDrug_ab, kDrug_ab, kDrug_ab, kDrug_ab, kDrug_ab,

        //10 - 29
        kDrug_ab, kDrug_b, kDrug_a, kDrug_none, kDrug_other, // kRab
        kDrug_ab, kDrug_b, kDrug_a, kDrug_none, kDrug_other, // kRa
        kDrug_ab, kDrug_b, kDrug_a, kDrug_none, kDrug_other, // kRb
        kDrug_ab, kDrug_b, kDrug_a, kDrug_none, kDrug_other, // kR0

        //30 - 39
        kDrug_none, kDrug_none, kDrug_none, kDrug_none, kDrug_none, // no drug
        kDrug_other, kDrug_other, kDrug_other, kDrug_other, kDrug_other // no drug


    };

    for (int ss = 0; ss < kNum_systems; ss++) {
        human::BloodSystemManager::parasite_dominance_mmc_wp2_kernel(
            ss*kNum_pgs_per_system + kHead_index[ss],
            ss*kNum_pgs_per_system + kTail_index[ss],
            ss,
            dominant_type_array,
            most_advanced_stage_array,
            pg_stage_array,
            pg_type_array,
            drug_array,
            kDrug_ab,
            kDrug_a,
            kDrug_b
        );
    }

    REQUIRE(most_advanced_stage_array[0] == common::StageName::kLiver);
    REQUIRE(most_advanced_stage_array[1] == common::StageName::kBlood);
    REQUIRE(most_advanced_stage_array[2] == common::StageName::kInfectious);
    REQUIRE(most_advanced_stage_array[3] == common::StageName::kInfectious);
    REQUIRE(most_advanced_stage_array[4] == common::StageName::kInfectious);

    for( int ss = 5; ss < kNum_systems; ss++) {
        REQUIRE(most_advanced_stage_array[ss] == common::StageName::kInfectious);
    }

    REQUIRE(dominant_type_array[0] == common::ParasiteType::kNoParasite);
    REQUIRE(dominant_type_array[1] == common::ParasiteType::kNoParasite);
    REQUIRE(dominant_type_array[2] == common::ParasiteType::kR0);
    REQUIRE(dominant_type_array[3] == common::ParasiteType::kR0);
    REQUIRE(dominant_type_array[4] == common::ParasiteType::kR0);

    REQUIRE(dominant_type_array[5] == common::ParasiteType::kRa);
    REQUIRE(dominant_type_array[6] == common::ParasiteType::kRb);
    REQUIRE(dominant_type_array[7] == common::ParasiteType::kR0);
    REQUIRE(dominant_type_array[8] == common::ParasiteType::kRab);
    REQUIRE(dominant_type_array[9] == common::ParasiteType::kRab);

    for( int ss = 10; ss < 15; ss++) {
        REQUIRE(dominant_type_array[ss] == common::ParasiteType::kRab);
    }

    for( int ss = 15; ss < 20; ss++) {
        REQUIRE(dominant_type_array[ss] == common::ParasiteType::kRa);
    }

    for( int ss = 20; ss < 25; ss++) {
        REQUIRE(dominant_type_array[ss] == common::ParasiteType::kRb);
    }

    for( int ss = 25; ss < 30; ss++) {
        REQUIRE(dominant_type_array[ss] == common::ParasiteType::kR0);
    }

    int num_tests_rnd_pick = 100;

    for (int ii  = 0; ii < num_tests_rnd_pick; ii++) {

        for (int ss = 30; ss < kNum_systems; ss++) {
            human::BloodSystemManager::parasite_dominance_mmc_wp2_kernel(
                ss*kNum_pgs_per_system + kHead_index[ss],
                ss*kNum_pgs_per_system + kTail_index[ss],
                ss,
                dominant_type_array,
                most_advanced_stage_array,
                pg_stage_array,
                pg_type_array,
                drug_array,
                kDrug_ab,
                kDrug_a,
                kDrug_b
            );

        }

        REQUIRE(dominant_type_array[30] != common::ParasiteType::kR0);
        REQUIRE(dominant_type_array[31] != common::ParasiteType::kR0);
        REQUIRE(dominant_type_array[32] != common::ParasiteType::kRa);
        REQUIRE(dominant_type_array[33] != common::ParasiteType::kRb);
        REQUIRE(dominant_type_array[34] != common::ParasiteType::kRab);

        REQUIRE(dominant_type_array[35] != common::ParasiteType::kR0);
        REQUIRE(dominant_type_array[36] != common::ParasiteType::kR0);
        REQUIRE(dominant_type_array[37] != common::ParasiteType::kRa);
        REQUIRE(dominant_type_array[38] != common::ParasiteType::kRb);
        REQUIRE(dominant_type_array[39] != common::ParasiteType::kRab);

    }

}

TEST_CASE("30.4 Blood_System functions: func_clinical_probability_blood_to_infectious", "[blood_system:clinical_probability_blood_to_infectious_kernel]") {
    
    uint8_t kImmunity_level_max = 10;
    uint8_t kCummulative_exposures_max = 10;


    std::cout << "BloodSystemManager::func_clinical_probability_blood_to_infectious_kernel_2" << "\n";
    // BloodSystemManager::func_clinical_probability_blood_to_infectious_kernel_2(
    //     uint8_t immunity_level,
    //     uint8_t cummulative_exposures
    // ) 
    std::cout << "immunity_level,cumulative_exposures,clinical_probability\n";
    for(uint8_t ii = 0; ii < kImmunity_level_max; ii++) {
        for(uint8_t cc = 0; cc < kCummulative_exposures_max; cc++) {
            std::cout << (int) ii << "," << (int) cc << "," << human::BloodSystemManager::func_clinical_probability_blood_to_infectious_kernel_2(ii, cc) << "\n";
        }
    }


}