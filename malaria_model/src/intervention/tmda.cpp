#include <iostream>
#include <vector>

#include "util/statistics.h"

#include "village/village.h"

#include "intervention/mda.h"
#include "intervention/smda.h"
#include "intervention/tmda.h"

namespace intervention {

Tmda::Tmda(
        // Smda
        const int num_teams,

        const int num_days_spent_in_village,
        const int num_drug_rounds,
        const int num_days_between_drug_rounds,

        // Tmda
        const float target_size

    ) : 
        Smda(
            num_teams,
            num_days_spent_in_village,
            num_drug_rounds,
            num_days_between_drug_rounds
        ),

        kTarget_size(target_size) {

}

Tmda::~Tmda(){

}

void Tmda::step(
        const int current_time_step,
        village::VillageManager* vll_mgr_ptr
    ){

    if (this->triggered_mdas.size() > 0) {
        for (auto &mda : this->triggered_mdas) {
            mda.step(current_time_step, vll_mgr_ptr);

        }
    }

}

void Tmda::create_mda(
        const int current_time_step,
        const std::vector<float>& prevalence_per_village,
        const float* const* village_distances
    ){

    assert(current_time_step == this->scheduled_start_times.front());
    this->scheduled_start_times.erase(this->scheduled_start_times.begin());

    int num_targets = this->kTarget_size * prevalence_per_village.size();

    std::vector<size_t> village_ids_in_decending_prevalence =
        util::sort_by_value_and_return_indices_descend(prevalence_per_village);

    std::vector<int> list_of_targeted_villages;

    std::cout << "tMDA is set for the following " << num_targets << " villages:\n";
    for (int cc = 0; cc < num_targets; cc++) {

        int v_id = static_cast<int>(village_ids_in_decending_prevalence[cc]);
        list_of_targeted_villages.push_back(v_id);
        std::cout << v_id << "(" << prevalence_per_village[v_id] << "), ";

    }

    float** target_distances = new float* [num_targets];

    for(int vv1 = 0; vv1 < num_targets; vv1++ ){
        
        target_distances[vv1] = new float[num_targets];
        std::fill_n(target_distances[vv1], num_targets, std::numeric_limits<float>::infinity());

        for(int vv2 = 0; vv2 < num_targets; vv2++ ){
            target_distances[vv1][vv2] =
                village_distances[list_of_targeted_villages[vv1]][list_of_targeted_villages[vv2]];
        }

    }

    this->triggered_mdas.push_back(Mda(
        this->kNum_teams,
        num_targets,

        this->kNum_days_spent_in_village,
        this->kNum_drug_rounds,
        this->kNum_days_between_drug_rounds,
        
        target_distances,

        list_of_targeted_villages,
        current_time_step
    ));


    for(int vv = 0; vv < num_targets; vv++ ){
        delete[] target_distances[vv];
    }
    delete[] target_distances;


}

void Tmda::print_mda_map(
        const int num_villages,
        const float* longitude_of,
        const float* latitude_of,
        const int* population_of,
        const std::string& output_prefix
    ){

    int num_targets = this->kTarget_size * num_villages;
    const std::vector<int>& list_of_targeted_villages = this->triggered_mdas.back().get_village_ids();
    assert(num_targets == static_cast<int>(list_of_targeted_villages.size()));

    float* lng = new float[num_targets];
    float* lat = new float[num_targets];
    int* pop = new int[num_targets];

    for(int vv = 0; vv < num_targets; vv++ ){
        lng[vv] = longitude_of[list_of_targeted_villages[vv]];
        lat[vv] = latitude_of[list_of_targeted_villages[vv]];
        pop[vv] = population_of[list_of_targeted_villages[vv]];
    }

    std::string sys_command = this->triggered_mdas.back().print_map(
        output_prefix,
        lng,
        lat,
        pop
    ) + "&";
    if(system(sys_command.c_str())){}

    delete[] lng;
    delete[] lat;
    delete[] pop;
}

}