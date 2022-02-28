#include <iostream>
#include <bitset>
// #include <iomanip>

// #include <string>
// #include <cmath>

#include <algorithm> // std::fill_n
#include <vector>

#include "third_party/catch2/catch.hpp"
#include "util/randomness.h"
#include "util/util.h"
#include "util/statistics.h"

#include "human/genotype.h"

TEST_CASE( "90.1: Genotype: bitwise operations for mutation definition", "[genotype:basic_operations]" ) {

    const int kNum_tests = 5000;

    const int kNum_bits_in_byte = 8;
    const int kNum_positions = sizeof(human::genotype_t) * kNum_bits_in_byte;

    // const float kMutation_probability = 0.1;
    const int kNum_probabilities = 4; // each probability is applied to all positions when tested
    const float kMutation_probability[kNum_probabilities] = {0.0, 0.2, 0.7, 1.0};

    int num_mutations = 0;

    human::genotype_t gt = 0;
    human::genotype_t gt_temp = 0;
    human::genotype_t mask = 0;

    float rnd_nmb = 0.0;

    int num_printed_mutations = 3;

    for (int pp = 0; pp < kNum_probabilities; pp++){
        num_mutations = 0;

        std::cout << "Test 90.1/" << pp
                <<":\nFor a genotype with " << kNum_positions
                << " positions, each with " << kMutation_probability[pp] * 100 << "%"
                << " chance to mutate:\n";
        
        for (int tt = 0; tt < kNum_tests; tt++) {
            rnd_nmb = 0.0;
            mask = 0;

            if (tt < num_printed_mutations) {
                std::cout << "\nTime step " << tt;
                std::cout << " Genotype (before): " << std::bitset<kNum_positions>(gt) << ". ";
                std::cout << "Mask construction:\n";
            }

            for (uint bb = 0; bb < kNum_positions; bb++) {
                util::set_random_numbers_uniform(&rnd_nmb, 1);
                if (tt < num_printed_mutations){
                    std::cout << "Mutation Mask: from " << std::bitset<kNum_positions>(mask);
                }
                mask <<= 1;
                mask += (rnd_nmb < kMutation_probability[pp]) ? 1 : 0;
                if (tt < num_printed_mutations){
                    std::cout << " to " << std::bitset<kNum_positions>(mask)
                                << ", rnd_nmb: " << rnd_nmb << "(" << kMutation_probability[pp] << ")\n";
                }
            }
            gt_temp = gt;
            gt ^= mask;

            gt_temp ^= gt; // the number of ones in gt_temp is the number of bits that were flipped by mask
            num_mutations += std::bitset< kNum_positions >(gt_temp).count();

            if (tt < num_printed_mutations){
                std::cout << "Genotype (after): " << std::bitset<kNum_positions>(gt)
                        << " with " << std::bitset< kNum_positions >(gt_temp).count() << " mutations. \n";
            }

        }


        int target_num_mutations = kNum_tests*kNum_positions*kMutation_probability[pp];
        float tolerance = 0.02; // 1%

        std::cout << "\nAfter " << kNum_tests << " time steps, "
                  << num_mutations << " mutations happend, expecting "
                  << target_num_mutations << " +/-" << tolerance * 100 << "% ... ";
        REQUIRE(static_cast<float>(std::abs(num_mutations - target_num_mutations)) <= target_num_mutations * tolerance);
        std::cout << "PASS\n\n";
        
    }

    std::cout << "\nExample: \nGenotype (before mutation): " << std::bitset<kNum_positions>(gt)
                << " (" << std::bitset<kNum_positions>(gt).count() << " ones)" << "\n";

    std::bitset<kNum_positions> mask_b(mask);
    std::cout << "Mutation Mask: " << mask_b
                << " (" << mask_b.count() << " ones)" << "\n";
    gt ^= mask;
    std::cout << "Genotype (after mutation): "
                << std::bitset<kNum_positions>(gt)
                << " (" << std::bitset<kNum_positions>(gt).count() << " ones)" "\n";

}

TEST_CASE( "90.1.1: Genotype: Initilisation", "[genotype:init]" ) {

    const int kMax_num_sites_in_combined_mutation = 2;
    const std::vector<float> kMutation_probability{0.1, 0.2, 0.3, 0.4};
    
    human::GenotypeManager gm(kMutation_probability, kMax_num_sites_in_combined_mutation);


    gm.print_all_combined_probabilities();

    std::vector<float> all_probabilities = gm.get_all_combined_probabilities();


    REQUIRE(
        all_probabilities.size() ==
            util::from_n_choose_up_to_k_get_num_combinations(
                kMutation_probability.size(), kMax_num_sites_in_combined_mutation
            )
    );

    float non_mutation_probability = 1.0;
    for (uint ss = 0; ss < kMutation_probability.size(); ss++) {
        non_mutation_probability *= (1-kMutation_probability.at(ss));
    }

    REQUIRE( all_probabilities.at(0) == non_mutation_probability );
}

TEST_CASE( "90.2: Genotype: step_mutation_independent_kernel function", "[genotype:step_mutation_independent_kernel]" ) {

    const int kNum_positions = 8;
    const float kMutation_probability_max = 1.0;
    float mutation_probability_array[kNum_positions] = {};
    int position_mask_array[kNum_positions] = {};
    for (int ii = 0; ii < kNum_positions; ii++) {
        mutation_probability_array[ii] = kMutation_probability_max 
                                            / static_cast<float>(kNum_positions)
                                            * static_cast<float>(ii);
        position_mask_array[ii] = 1<<(kNum_positions-ii-1);
    }
    mutation_probability_array[kNum_positions] = kMutation_probability_max;


    const int kNum_infections = 100;
    human::genotype_t gt_array[kNum_infections] = {};
    std::fill(gt_array, gt_array + kNum_infections, 0);


    const int kNum_time_steps = 3650;

    int num_mutations[kNum_positions];
    std::fill(num_mutations, num_mutations + kNum_positions, 0);
    human::genotype_t gt_array_temp[kNum_infections];

    for (int tt = 0; tt < kNum_time_steps; tt++) {

        std::copy(gt_array, gt_array + kNum_infections, gt_array_temp);

        human::GenotypeManager::step_mutation_independent_kernel(
            gt_array,
            kNum_infections,

            mutation_probability_array,
            kNum_positions
        );

        for (int gg = 0; gg < kNum_infections; gg++) {
            gt_array_temp[gg] ^= gt_array[gg];
            for (int ii = 0; ii < kNum_positions; ii++) {
                if ((gt_array_temp[gg] & position_mask_array[ii]) == position_mask_array[ii]) {
                    num_mutations[ii]++;
                }
            }
        }

    }

    float tolerance = 0.01; // 1%

    std::cout << "After " << kNum_time_steps << " steps of " << kNum_infections << " infections:\n";
    for (int ii = 0; ii < kNum_positions; ii++) {
        int num_expected_mutations = mutation_probability_array[ii] * kNum_time_steps * kNum_infections;
        std::cout << "Position " << ii 
                    << " with mutation probability " << mutation_probability_array[ii]
                    << ", number of actual mutations:" << num_mutations[ii]
                    << ", expecting " << num_expected_mutations
                    << " +/- " << tolerance*100 << "% ...";
        REQUIRE(std::abs(num_expected_mutations - num_mutations[ii]) <= num_expected_mutations * tolerance);
        std::cout << " PASS\n";
    }

}


TEST_CASE("90.3: Genotype: killing function performance test", "[genotype:kill_perf]") {
    const int kNum_tests = 10;

    const int kNum_genotypes = 128;
    const int kNum_majority_types = 1;
    const int kNum_minority_types = kNum_genotypes - kNum_majority_types;
    // const int kMajority_type_count = 2000000*10;
    const int64_t kMajority_type_count = 10000000000;
    const int64_t kMinority_type_count = 200;
    const int64_t kTotal_count = kMajority_type_count*kNum_majority_types
                            + kMinority_type_count*kNum_minority_types;
    // const int64_t kKill_count = kTotal_count * 0.9;
    const int64_t kKill_count_left = 1000000;
    const int64_t kKill_count = kTotal_count - kKill_count_left;

    int64_t parasite_count[kNum_genotypes] = {};
    int64_t parasite_count_survives[kNum_genotypes] = {};

    std::fill(parasite_count, parasite_count+kNum_majority_types, kMajority_type_count);
    std::fill(parasite_count+kNum_majority_types, parasite_count+kNum_genotypes, kMinority_type_count);

    std::fill(parasite_count_survives, parasite_count_survives+kNum_genotypes, 0);



    double wall_time_begin = 0.0;
    double wall_time_end = 0.0;


    wall_time_begin = util::get_wall_time();
    
    // float uniform_random_array[kKill_count] = {};
    float rnd_nmb= 0.0;
    for (int tt = 0; tt < kNum_tests; tt++ ){

        std::fill(parasite_count, parasite_count+kNum_majority_types, kMajority_type_count);
        std::fill(parasite_count+kNum_majority_types, parasite_count+kNum_genotypes, kMinority_type_count);

        std::fill(parasite_count_survives, parasite_count_survives+kNum_genotypes, 0);



        // util::set_random_numbers_uniform(uniform_random_array, kKill_count);

        int64_t total_cout = kTotal_count;
        // for (int64_t kk = 0; kk < kKill_count; kk++) {
        for (int64_t kk = 0; kk < kKill_count_left; kk++) {

            // int kill_at = uniform_random_array[kk] * total_cout;
            util::set_random_numbers_uniform(&rnd_nmb, 1);
            // int kill_at = rnd_nmb * total_cout;

            int survive_at = rnd_nmb * kTotal_count;

            int accumulator = 0;

            for (int tt = 0; tt < kNum_genotypes; tt++) {
                // if (parasite_count[tt] + accumulator >= kill_at) {
                if (parasite_count[tt] + accumulator >= survive_at) {
                    // parasite_count[tt]--;
                    // total_cout--;

                    parasite_count_survives[tt]++;
                    break;
                } else {
                    accumulator += parasite_count[tt];
                }
            }
        }



        int64_t remain_count = 0;
        
        for (int tt = 0; tt < kNum_genotypes; tt++) {
            // remain_count += parasite_count[tt];
            remain_count += parasite_count_survives[tt];
        }

        // REQUIRE(remain_count == (kTotal_count - kKill_count));
        REQUIRE(remain_count == kKill_count_left);

    }
    wall_time_end = util::get_wall_time();

    std::cout << "Per-kill sample with cascaded loop locator: " << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";




    // wall_time_begin = util::get_wall_time();
    
    // // float uniform_random_array[kKill_count] = {};
    // // std::vector<int> all_parasites(kTotal_count, 0);
    // // int kill_array[kKill_count/200] = {};
    // int* kill_array = new int[kKill_count];
    // std::vector<int> parasite_histogram(kNum_genotypes, 0);

    // for (int tt = 0; tt < kNum_tests; tt++ ){

    //     std::fill(parasite_count, parasite_count+kNum_majority_types, kMajority_type_count);
    //     std::fill(parasite_count+kNum_majority_types, parasite_count+kNum_genotypes, kMinority_type_count);


    //     for (int tt = 0; tt < kNum_genotypes; tt++) {
    //         parasite_histogram[tt] = parasite_count[tt];
    //     }

    //     util::set_random_numbers_discrete<int>(kill_array, kKill_count, parasite_histogram);

    //     for (int kk = 0; kk < kKill_count; kk++) {
    //         parasite_count[kill_array[kk]]--;
    //     }

    //     int remain_count = 0;
    //     for (int tt = 0; tt < kNum_genotypes; tt++) {
    //         remain_count += parasite_count[tt];
    //     }
    //     REQUIRE(remain_count == (kTotal_count - kKill_count));

    // }
    // delete[] kill_array;
    // wall_time_end = util::get_wall_time();

    // std::cout << "Per-kill sample from discrete distribution: " << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";







    // wall_time_begin = util::get_wall_time();
    
    // std::vector<int> all_parasites(kTotal_count, 0);
    // for (int tt = 0; tt < kNum_tests; tt++ ){

    //     std::fill(parasite_count, parasite_count+kNum_majority_types, kMajority_type_count);
    //     std::fill(parasite_count+kNum_majority_types, parasite_count+kNum_genotypes, kMinority_type_count);

    //     int accumulator = 0;
    //     for (int tt = 0; tt < kNum_genotypes; tt++) {
    //         std::fill(  all_parasites.begin() + accumulator, 
    //                     all_parasites.begin() + accumulator + parasite_count[tt],
    //                     tt
    //         );
    //         accumulator += parasite_count[tt];
    //     }
    //     std::random_shuffle(all_parasites.begin(), all_parasites.end());

    //     std::fill(parasite_count, parasite_count + kNum_genotypes, 0);
    //     for (int pp = kKill_count; pp < kTotal_count; pp++) {
    //         parasite_count[all_parasites.at(pp)]++;
    //     }

    //     int remain_count = 0;
    //     for (int tt = 0; tt < kNum_genotypes; tt++) {
    //         remain_count += parasite_count[tt];
    //     }
    //     REQUIRE(remain_count == (kTotal_count - kKill_count));

    // }
    // wall_time_end = util::get_wall_time();

    // std::cout << "Lineup with shuffle and truncate: " << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";





    // wall_time_begin = util::get_wall_time();

    // std::vector<bool> all_parasites_kill_flag(kTotal_count, false);
    // std::fill(all_parasites_kill_flag.begin(), all_parasites_kill_flag.end(), false);

    
    // for (int tt = 0; tt < kNum_tests; tt++ ){

    //     std::fill(parasite_count, parasite_count+kNum_majority_types, kMajority_type_count);
    //     std::fill(parasite_count+kNum_majority_types, parasite_count+kNum_genotypes, kMinority_type_count);


    //     // util::set_random_numbers_uniform(uniform_random_array, kKill_count);

    //     int accumulator = 0;
    //     for (int tt = 0; tt < kNum_genotypes; tt++) {
    //         std::fill(  all_parasites.begin() + accumulator, 
    //                     all_parasites.begin() + accumulator + parasite_count[tt],
    //                     tt
    //         );
    //         accumulator += parasite_count[tt];
    //     }
    //     std::fill(all_parasites_kill_flag.begin(), all_parasites_kill_flag.end(), false);

    //     for (int kk = 0; kk < kKill_count; kk++) {

    //         // int kill_at = uniform_random_array[kk] * total_cout;
    //         float rnd_nmb = 0.0;
    //         util::set_random_numbers_uniform(&rnd_nmb, 1);
    //         int kill_at = rnd_nmb * kTotal_count;

    //         while(all_parasites_kill_flag[kill_at]) {
    //             kill_at++;
    //             if (kill_at >= kTotal_count) {
    //                 kill_at = 0;
    //             }
    //         }

    //         all_parasites_kill_flag[kill_at] = true;
    //         parasite_count[all_parasites.at(kill_at)]--;

    //     }

    //     int remain_count = 0;
    //     for (int tt = 0; tt < kNum_genotypes; tt++) {
    //         remain_count += parasite_count[tt];
    //     }
    //     REQUIRE(remain_count == (kTotal_count - kKill_count));

    // }
    // wall_time_end = util::get_wall_time();

    // std::cout << "Per-kill sample from lineup with flag for removal: " << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";




    // Not a good idea (erase slow) >>

    // wall_time_begin = util::get_wall_time();
    
    // for (int tt = 0; tt < kNum_tests; tt++ ){

    //     std::fill(parasite_count, parasite_count+kNum_majority_types, kMajority_type_count);
    //     std::fill(parasite_count+kNum_majority_types, parasite_count+kNum_genotypes, kMinority_type_count);

    //     util::set_random_numbers_uniform(uniform_random_array, kKill_count);

    //     int accumulator = 0;
    //     all_parasites.resize(kTotal_count);
    //     for (int tt = 0; tt < kNum_genotypes; tt++) {
    //         std::fill(  all_parasites.begin() + accumulator, 
    //                     all_parasites.begin() + accumulator + parasite_count[tt],
    //                     tt
    //         );
    //         accumulator += parasite_count[tt];
    //     }

    //     for (int kk = 0; kk < kKill_count; kk++) {

    //         int kill_at = uniform_random_array[kk] * (kTotal_count-kk);

    //         parasite_count[all_parasites.at(kill_at)]--;
    //         all_parasites.erase(all_parasites.begin() + kill_at);

    //         std::cout << kk << "\n";

    //     }

    //     int remain_count = 0;
    //     for (int tt = 0; tt < kNum_genotypes; tt++) {
    //         remain_count += parasite_count[tt];
    //     }
    //     REQUIRE(remain_count == (kTotal_count - kKill_count));

    // }
    // wall_time_end = util::get_wall_time();

    // std::cout << "Per-kill sample from lineup with per-kill lineup reduction: " << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";

}