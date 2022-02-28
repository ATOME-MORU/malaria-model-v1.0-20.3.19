#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdint> //uint8_t

#include <bitset>

#include "util/randomness.h"
#include "util/statistics.h"

#include "human/genotype.h"

namespace human {

GenotypeManager::GenotypeManager(

        const std::vector<float>& mutation_probability_vector,
        const int max_num_sites_in_combined_mutation
    
    ) : kNum_sites_in_genotype_code(static_cast<int>(mutation_probability_vector.size())),
        kMax_num_sites_in_combined_mutation(max_num_sites_in_combined_mutation),
        kNum_mutation_combinations(
            util::from_n_choose_up_to_k_get_num_combinations(
                kNum_sites_in_genotype_code,
                kMax_num_sites_in_combined_mutation
            )
        ) {

    // Independent Mutation Probability
    this->kMutation_independent_probability = new float[this->kNum_sites_in_genotype_code];
    for (int pp = 0; pp < this->kNum_sites_in_genotype_code; pp++) {
        this->kMutation_independent_probability[pp] = mutation_probability_vector.at(pp);
    }

    // Combined Mutation Probability
    this->kMutation_combined_probability = new float[this->kNum_mutation_combinations];
    this->kMutation_combined_masks = new genotype_t[this->kNum_mutation_combinations];
    std::fill(this->kMutation_combined_probability,
        this->kMutation_combined_probability +kNum_mutation_combinations,
        1
    );
    std::fill(this->kMutation_combined_masks,
        this->kMutation_combined_masks +kNum_mutation_combinations,
        0
    );

    int num_combinations_accumulator = 0;
    // no mutation at all sites
    this->kMutation_combined_masks[num_combinations_accumulator] = 0;
    for (int ss = 0; ss < kNum_sites_in_genotype_code; ss++) {
        this->kMutation_combined_probability[num_combinations_accumulator]
            *= (1 - this->kMutation_independent_probability[ss]);
    }
    num_combinations_accumulator++;

    // mutation at least one site
    for(int kk = 1; kk <= kMax_num_sites_in_combined_mutation; kk++) {
        std::vector<std::vector<int>> choose_kk_combination_set = util::from_n_choose_k(kNum_sites_in_genotype_code, kk);

        for (auto& each_combination : choose_kk_combination_set) {
            
            for (int ss = 0; ss < kNum_sites_in_genotype_code; ss++) {
                if (std::find(each_combination.begin(), each_combination.end(), ss) !=  each_combination.end()) {
                    this->kMutation_combined_masks[num_combinations_accumulator] <<= 1;
                    this->kMutation_combined_masks[num_combinations_accumulator] += 1;

                    this->kMutation_combined_probability[num_combinations_accumulator]
                        *= this->kMutation_independent_probability[ss];

                } else {

                    this->kMutation_combined_masks[num_combinations_accumulator] <<= 1;

                    this->kMutation_combined_probability[num_combinations_accumulator]
                        *= (1 - this->kMutation_independent_probability[ss]);

                }
            }

            num_combinations_accumulator++;
        }
    }

}

GenotypeManager::~GenotypeManager() {
    if (kMutation_independent_probability) delete[] kMutation_independent_probability;
    if (kMutation_combined_probability) delete[] kMutation_combined_probability;
    if (kMutation_combined_masks) delete[] kMutation_combined_masks;
}

void GenotypeManager::print_all_combined_probabilities() const {

    const int kMax_num_sites_in_genotype_code = 8;  //this is for functions such as std::bitset<int>
                                                    //since the compiler needs to have that int
                                                    //at compile time

    std::cout << "Given that the *independent* mutation probabilities at each of the "
            << this->kNum_sites_in_genotype_code << " sites are as [";
    for (int ss = 0; ss < this->kNum_sites_in_genotype_code; ss++) {
        std::cout << " " << this->kMutation_independent_probability[ss];
    }
    std::cout << " ]\n";

    std::cout << "Followings are the probabilities of possible combinations of mutations:\n";
    std::cout << "(Max number of concurrent mutation events: " << kMax_num_sites_in_combined_mutation << ")\n";

    std::cout << "[combined mutation event: probability]\n";
    for (int pp = 0; pp < kNum_mutation_combinations; pp++) {
        std::cout << std::bitset<kMax_num_sites_in_genotype_code>(this->kMutation_combined_masks[pp])
                    << ": " << this->kMutation_combined_probability[pp] << "\n";
    }
    std::cout << "with [" << std::bitset<kMax_num_sites_in_genotype_code>(0)
        << "] being the situation where no mutation occurs at any sites.\n";
}

std::vector<float> GenotypeManager::get_all_combined_probabilities() const {
    std::vector<float> all_probabilities;
    for (int pp = 0; pp < kNum_mutation_combinations; pp++) {
        all_probabilities.push_back(this->kMutation_combined_probability[pp]);
    }
    return all_probabilities;
}

void GenotypeManager::step_mutation(
        genotype_t* genotype_array,
        const int genotype_array_size
    ){

    this->step_mutation_independent_kernel(
        genotype_array,
        genotype_array_size,

        this->kMutation_independent_probability,
        this->kNum_sites_in_genotype_code
    );

}

// void GenotypeManager::step_mutation_combined_kernel(
//         genotype_t* genotype_array,
//         const int genotype_array_size,
        

//     ) {


// }

void GenotypeManager::step_mutation_independent_kernel(
        genotype_t* genotype_array,
        const int genotype_array_size,

        const float* mutation_probability_independent_array,
        const int mutation_probability_array_size
    ) {

    genotype_t mask = 0;
    float* rnd_nmb_array = new float[mutation_probability_array_size];

    for (int gg = 0; gg < genotype_array_size; gg++) {

        mask = 0;
        util::set_random_numbers_uniform(rnd_nmb_array, mutation_probability_array_size);
        
        for (int bb = 0; bb < mutation_probability_array_size; bb++) {
            mask <<= 1;
            mask += (rnd_nmb_array[bb] < mutation_probability_independent_array[bb]) ? 1 : 0;
        }

        genotype_array[gg] ^= mask;

    }

    delete[] rnd_nmb_array;

    

}

}