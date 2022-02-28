#ifndef GENOTYPE_H
#define GENOTYPE_H


namespace human {

typedef uint8_t genotype_t;


class GenotypeManager {

    const int kNum_bits_in_a_byte = 8;

    const int kNum_sites_in_genotype_code; 

    // Probability of mutation at one site regardless of other events
    float* kMutation_independent_probability;

    // Probability of combination of mutation events at all sites
    float* kMutation_combined_probability;
    genotype_t* kMutation_combined_masks;
    const int kMax_num_sites_in_combined_mutation;
    const int kNum_mutation_combinations;

public:

    GenotypeManager(
        const std::vector<float>& mutation_probability_vector,
        const int max_num_sites_in_combined_mutation
    );
    ~GenotypeManager();

    void print_all_combined_probabilities() const;
    std::vector<float> get_all_combined_probabilities() const;

    void step_mutation(
        genotype_t* genotype_array,
        const int genotype_array_size
    );
    static void step_mutation_independent_kernel(
        genotype_t* genotype_array,
        const int genotype_array_size,

        const float* mutation_probability_independent_array,
        const int mutation_probability_array_size
    );

};



}
#endif