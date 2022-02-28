#ifndef TMDA_H
#define TMDA_H


namespace intervention {

// Targeted MDA:
// Target a percentage of all villages starting
// with the highest prevalence
class Tmda: public Smda{ 

    const float kTarget_size;
    // e.g. 0.1 -> MDA top 10% of all villages by prevalence

public:
    Tmda(
        // Smda
        const int num_teams,

        const int num_days_spent_in_village,
        const int num_drug_rounds,
        const int num_days_between_drug_rounds,

        // Tmda
        const float target_size

    );
    ~Tmda();

    void step(
        const int current_time_step,
        village::VillageManager* vll_mgr_ptr
    );

    void create_mda(
        const int current_time_step,
        const std::vector<float>& prevalence_per_village,
        const float* const* village_distances
    );

    void print_mda_map(
        const int num_villages,
        const float* longitude_of,
        const float* latitude_of,
        const int* population_of,
        const std::string& output_prefix
    );

};

}

#endif