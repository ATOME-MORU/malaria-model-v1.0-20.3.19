#include "third_party/catch2/catch.hpp"
#include <iostream>
#include <limits> //std::numeric_limits

#include <algorithm>

#include <map>
#include <vector>

#include <cmath>

#include "village/village.h"
#include "common/common.h"

// #include "util/util.h"
#include "util/randomness.h"

// #include "common_test.h"


TEST_CASE( "50.1: Treatment Functions: treatment_by_malaria_post_kernel", "[village:treatment_by_malaria_post_kernel]" ) {
    
    const common::DrugName target_drug = common::DrugName::kPiparte;

    const uint test_size = 900000;
    const float target_rate_clinical = 0.5;
    const float target_rate_asymptomatic = 0.3;
    const float test_tolerence = 0.001;

    float result_rate = 0;

    const uint num_villages = 1;

    common::DrugName existing_drug[test_size] = {};
    std::fill(existing_drug, existing_drug + test_size, common::DrugName::kNoDrug);

    common::DrugName given_drug[test_size] = {};
    std::fill(given_drug, given_drug + test_size, common::DrugName::kNoDrug);

    bool is_clinical[test_size] = {};
    std::fill(is_clinical, is_clinical + test_size, true);

    common::ParasiteType dominant_type[test_size] = {};
    std::fill(dominant_type, dominant_type + test_size, common::ParasiteType::kR0);

    int num_humans = test_size;

    std::vector<std::vector<int>> at_vll_reg;
    at_vll_reg.resize(num_villages);

    for (uint hh = 0; hh < test_size; hh++) {
        at_vll_reg[0].push_back(hh);
    }
    int at_vll_reg_size = num_villages;
    float rnd_nmbrs[test_size] = {};
    util::set_random_numbers_uniform(rnd_nmbrs, test_size);

    float rates_clinical[num_villages] = {};
    std::fill(rates_clinical, rates_clinical + num_villages, target_rate_clinical);

    float rates_asymptomatic[num_villages] = {};
    std::fill(rates_asymptomatic, rates_asymptomatic + num_villages, target_rate_asymptomatic);

    village::VillageManager::treatment_by_malaria_post_kernel(
            existing_drug,
            given_drug,
            is_clinical,
            dominant_type,
            num_humans,
            at_vll_reg,
            at_vll_reg_size,
            rnd_nmbrs,
            rates_clinical,
            rates_asymptomatic
        );

    std::map<common::DrugName, int> drug_accumulator_clinical;
    for (uint hh = 0; hh < test_size; hh++) {
        drug_accumulator_clinical[given_drug[hh]]++;
    }

    result_rate = drug_accumulator_clinical[target_drug] / (float) test_size;

    std::cout << "Result from all clinical: \n target_drug:test_size ratio is "
                << result_rate
                << "(" << drug_accumulator_clinical[target_drug] << ":" << test_size << "), test target is "
                << target_rate_clinical << " +/- " << test_tolerence
                << std::endl;

    REQUIRE(std::abs(result_rate-target_rate_clinical) < test_tolerence);


    std::fill(is_clinical, is_clinical + test_size, false);
    std::fill(given_drug, given_drug + test_size, common::DrugName::kNoDrug);

    village::VillageManager::treatment_by_malaria_post_kernel(
            existing_drug,
            given_drug,
            is_clinical,
            dominant_type,
            num_humans,
            at_vll_reg,
            at_vll_reg_size,
            rnd_nmbrs,
            rates_clinical,
            rates_asymptomatic
        );

    std::map<common::DrugName, int> drug_accumulator_asymptomatic;
    for (uint hh = 0; hh < test_size; hh++) {
        drug_accumulator_asymptomatic[given_drug[hh]]++;
    }

    result_rate = drug_accumulator_asymptomatic[target_drug] / (float) test_size;

    std::cout << "Result from all asymptomatic: \n target_drug:test_size ratio is "
                << result_rate
                << "(" << drug_accumulator_asymptomatic[target_drug] << ":" << test_size << "), test target is "
                << target_rate_asymptomatic << " +/- " << test_tolerence
                << std::endl;

    REQUIRE(std::abs(result_rate-target_rate_asymptomatic) < test_tolerence);





    std::map<common::DrugName, int> drug_accumulator;
    float target_rate = 0;

    const float percent_no_infections = 0.001;
    target_rate = target_rate_asymptomatic * (1-percent_no_infections);

    std::fill(given_drug, given_drug + test_size, common::DrugName::kNoDrug);
    std::fill(dominant_type, dominant_type + int(test_size*percent_no_infections), common::ParasiteType::kNoParasite);

    village::VillageManager::treatment_by_malaria_post_kernel(
            existing_drug,
            given_drug,
            is_clinical,
            dominant_type,
            num_humans,
            at_vll_reg,
            at_vll_reg_size,
            rnd_nmbrs,
            rates_clinical,
            rates_asymptomatic
        );

    for (uint hh = 0; hh < test_size; hh++) {
        drug_accumulator[given_drug[hh]]++;
    }

    result_rate = drug_accumulator[target_drug] / (float) test_size;

    std::cout << "Result from partial ("
                << percent_no_infections << ") non-infectious): \n target_drug:test_size ratio is "
                << result_rate
                << "(" << drug_accumulator[target_drug] << ":" << test_size << "), test target is "
                << target_rate << " +/- " << test_tolerence
                << std::endl;

    REQUIRE(std::abs(result_rate-target_rate) < test_tolerence);






    const float percent_with_drug = 0.11;
    target_rate = target_rate_asymptomatic * (1-percent_with_drug);

    std::fill(given_drug, given_drug + test_size, common::DrugName::kNoDrug);
    drug_accumulator.clear();

    std::fill(dominant_type, dominant_type + test_size, common::ParasiteType::kR0);

    std::fill(existing_drug, existing_drug + int(test_size*percent_with_drug), common::DrugName::kPiperaquine);

    village::VillageManager::treatment_by_malaria_post_kernel(
            existing_drug,
            given_drug,
            is_clinical,
            dominant_type,
            num_humans,
            at_vll_reg,
            at_vll_reg_size,
            rnd_nmbrs,
            rates_clinical,
            rates_asymptomatic
        );

    for (uint hh = 0; hh < test_size; hh++) {
        drug_accumulator[given_drug[hh]]++;
    }

    result_rate = drug_accumulator[target_drug] / (float) test_size;

    std::cout << "Result from partial ( "
                << percent_with_drug <<") with-drug: \n target_drug:test_size ratio is "
                << result_rate
                << "(" << drug_accumulator[target_drug] << ":" << test_size << "), test target is "
                << target_rate << " +/- " << test_tolerence
                << std::endl;

    REQUIRE(std::abs(result_rate-target_rate) < test_tolerence);

}


TEST_CASE( "50.2: Treatment Functions: treatment_by_mda_kernel", "[village:treatment_by_mda_kernel]" ) {
    
    const common::DrugName target_drug = common::DrugName::kPiparte;

    const uint test_size = 900000;
    const float target_rate = 0.2;
    const float test_tolerence = 0.001;

    float result_rate = 0;

    const uint num_villages = 1;

    common::DrugName given_drug[test_size] = {};
    std::fill(given_drug, given_drug + test_size, common::DrugName::kNoDrug);

    int num_humans = test_size;

    std::vector<std::vector<int>> at_vll_reg;
    at_vll_reg.resize(num_villages);
    for (uint hh = 0; hh < test_size; hh++) {
        at_vll_reg[0].push_back(hh);
    }
    float rnd_nmbrs[test_size] = {};
    util::set_random_numbers_uniform(rnd_nmbrs, test_size);

    village::VillageManager::treatment_by_mda_kernel(
            given_drug,
            num_humans,
            at_vll_reg[0],
            rnd_nmbrs,
            target_rate,
            target_drug
        );

    std::map<common::DrugName, int> result_accumulator;
    for (uint hh = 0; hh < test_size; hh++) {
        result_accumulator[given_drug[hh]]++;
    }

    result_rate = result_accumulator[target_drug] / (float) test_size;

    std::cout << "Result: \n target_drug(" << target_drug << "):test_size ratio is "
                << result_rate
                << "(" << result_accumulator[target_drug] << ":" << test_size << "), test target is "
                << target_rate << " +/- " << test_tolerence
                << std::endl;

    REQUIRE(std::abs(result_rate-target_rate) < test_tolerence);

}