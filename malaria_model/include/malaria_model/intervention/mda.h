#ifndef MDA_H
#define MDA_H

#include <string>
#include <map> // to use vector instead

#include "util/graph.h"

#include "village/village.h"
#include "human/human.h"
#include "common/common.h"

namespace intervention {

class Mda{

    const int kAssigned_to_team_default = -1;

    // const int kDays_per_village = 4;
    // const int kVisits_per_village = 3;
    // const int kDays_between_visits = 32;
    
    const int kDays_per_village;
    const int kVisits_per_village;
    const int kDays_between_visits;

    // const int kRoute_step_offset = kDays_between_visits/kDays_per_village;
    const int kRoute_step_offset;

    // const bool kStart_mda_from_sorted_village_list;
    // std::vector<int> village_ids_sorted;

    const common::DrugName kMda_drug = common::DrugName::kPiparte;

    util::VillageGraph vg;
    std::map<int, std::vector<int> > routes;

    std::vector<int> village_ids;

    // std::vector<int> current_village_of;

    std::vector<int> assigned_to_team; // size: number of villages

    std::vector<int> scheduled_start_times;
    uint num_executed_mdas = 0;

    bool in_execution = false;
    int execution_step = 0;
    int route_step = 0;

    // std::vector<bool> visited_in_execution;
    std::vector<int> visited_days_in_execution;

public:

    Mda(
        const int num_teams,
        const int num_villages,

        const int days_per_village,
        const int visits_per_village,
        const int days_between_visits,

        const float* const* distances
    );
    // Mda(
    //     const int num_teams,
    //     const int num_villages,
 
    //     const int days_per_village,
    //     const int visits_per_village,
    //     const int days_between_visits,

    //     const float* const* distances,

    //     const bool start_mda_from_sorted_village_list,
    //     const std::vector<int> village_id_sorted

    // );
    Mda(
        const int num_teams,
        const int num_villages,

        const int days_per_village,
        const int visits_per_village,
        const int days_between_visits,
        
        const float* const* distances,

        const std::vector<int> village_id_array,
        const int start_time_step
    );
    ~Mda();

    void init(
        const int num_teams,
        const int num_villages
    );
    void reroute_teams(
        const std::vector<int>& village_ids_sorted
    );

    const util::VillageGraph& get_graph();
    inline const std::vector<int>& get_village_ids() const{
        return this->village_ids;
    }


    std::string print_map( 
        const std::string file_prefix,
        const float* longitudes,
        const float* latitudes,
        const int* populations
    ) const;

    void print_all_routes() const;

    void register_start_time(int time_step);

    void step(int time_step, village::VillageManager* vll_mgr);

private:

    void coordinate_team_routes();

    void print_route(int route_id) const;

    int  get_next_start_time() const;
    void start();
    void execute(village::VillageManager* vll_mgr);
    void finish();

};

}

#endif