#include <iostream>
#include <vector>

#include <algorithm> //std::remove

#include <cassert>


#include "intervention/mda.h"
// #include "util/util.h"
#include "util/randomness.h"


namespace intervention {

// MDA
Mda::Mda(
        const int num_teams,
        const int num_villages,
 
        const int days_per_village,
        const int visits_per_village,
        const int days_between_visits,

        const float* const* distances

    ) : 

        kDays_per_village(days_per_village),
        kVisits_per_village(visits_per_village),
        kDays_between_visits(days_between_visits),
            kRoute_step_offset(kDays_between_visits/kDays_per_village),

        // kStart_mda_from_sorted_village_list(false),

        vg(num_villages, distances)

    {

    for (int vv = 0; vv < num_villages; vv++) {
        village_ids.push_back(vv);
    }

    assert(this->village_ids.size() == static_cast<uint>(num_villages));
    this->init(num_teams,num_villages);

}

// // MDA, sorted starting points
// Mda::Mda(
//         const int num_teams,
//         const int num_villages,
 
//         const int days_per_village,
//         const int visits_per_village,
//         const int days_between_visits,

//         const float* const* distances,

//         const bool start_mda_from_sorted_village_list,
//         const std::vector<int> village_id_sorted

//     ) : 

//         kDays_per_village(days_per_village),
//         kVisits_per_village(visits_per_village),
//         kDays_between_visits(days_between_visits),
//             kRoute_step_offset(kDays_between_visits/kDays_per_village),
        
//         kStart_mda_from_sorted_village_list(start_mda_from_sorted_village_list),
//         village_ids_sorted(village_id_sorted),

//         vg(num_villages, distances)

//     {

//     for (int vv = 0; vv < num_villages; vv++) {
//         village_ids.push_back(vv);
//     }

//     assert(this->village_ids.size() == static_cast<uint>(num_villages));
//     this->init(num_teams,num_villages);

// }

// MDA created by tMDA, fMDA, i.e. MDA in sub-set of villages
Mda::Mda(
        const int num_teams,
        const int num_villages,

        const int days_per_village,
        const int visits_per_village,
        const int days_between_visits,

        const float* const* distances,

        const std::vector<int> village_id_array,
        const int start_time_step

    ) : 

        kDays_per_village(days_per_village),
        kVisits_per_village(visits_per_village),
        kDays_between_visits(days_between_visits),
            kRoute_step_offset(kDays_between_visits/kDays_per_village),

        // kStart_mda_from_sorted_village_list(false),

        vg(num_villages, distances),

        village_ids(village_id_array) {

    assert(this->village_ids.size() == static_cast<uint>(num_villages));
    this->init(num_teams,num_villages);

    this->scheduled_start_times.push_back(start_time_step);
}

Mda::~Mda(){

}

void Mda::init(
        const int num_teams,
        const int num_villages
    ) {

    this->assigned_to_team.resize(num_villages, this->kAssigned_to_team_default);

    std::cout << "Mda: Constructing village map..." << std::flush;
    this->vg.mark_graph_with_min_spanning_tree();
    this->vg.reduce_graph_to_min_spanning_tree();
    std::cout << "[done]\n";

    float* random_numbers_teams = new float[num_teams];
    util::set_random_numbers_uniform(random_numbers_teams, num_teams);

    for (int tt = 0; tt < num_teams; tt++) {
        int start_village = (int) (random_numbers_teams[tt]*num_villages);
        // int start_village = tt;
        this->routes[tt] = this->vg.get_bfs(start_village);
    }

    delete[] random_numbers_teams;

    this->coordinate_team_routes();

    // this->visited_in_execution.resize(num_villages, false);
    this->visited_days_in_execution.resize(num_villages, 0);

    std::cout << "MDA created (teams: "<< num_teams <<", villages: " << num_villages << ")\n";
    this->print_all_routes();

    // std::cin.ignore();

}

void Mda::reroute_teams(
    const std::vector<int>& village_ids_sorted
    ) {
    // this->assigned_to_team.resize(this->routes.size(), this->kAssigned_to_team_default);
    for (auto& ii : this->assigned_to_team) {
        ii = this->kAssigned_to_team_default;
    }
    for (size_t tt = 0; tt < this->routes.size(); tt++) {
        this->routes[tt] = this->vg.get_bfs(village_ids_sorted.at(tt));
    }
    this->coordinate_team_routes();
    std::cout << "MDA rerouted:\n";
    this->print_all_routes();
}

const util::VillageGraph& Mda::get_graph(){
    return this->vg;
}

//TODO: std::list may be faster for removal
void Mda::coordinate_team_routes(){

    std::vector<uint> offset(this->routes.size(), 0);
    const int kRemoval_value = -1;
    
    for (size_t ii = 0; ii < this->vg.get_num_villages(); ii++){
        for (size_t tt = 0; tt < this->routes.size(); tt++) {

            if ((ii + offset[tt]) < this->routes[tt].size()){        

                while( assigned_to_team[this->routes[tt].at(ii + offset[tt])] != kAssigned_to_team_default ){
                    
                    this->routes[tt].at(ii + offset[tt]) = kRemoval_value;
                    offset[tt]++;

                    if ((ii + offset[tt]) >= this->routes[tt].size()) break;
                }

                if ((ii + offset[tt]) < this->routes[tt].size()){
                    assigned_to_team[ this->routes[tt].at(ii+offset[tt]) ] = tt;
                }
            }

        }
    }

    for (size_t tt = 0; tt < this->routes.size(); tt++) {
        this->routes[tt].erase(
            std::remove(this->routes[tt].begin(), this->routes[tt].end(), kRemoval_value),
            this->routes[tt].end());
    }
}

std::string Mda::print_map(
    const std::string file_prefix,
    const float* longitudes,
    const float* latitudes,
    const int* populations
    ) const{

    std::string file_name = "mda_plan";
    std::string plot_title = "MDA Plan";
    std::string group_name_prefix = "Team";
    std::vector<int> labels;
    return this->vg.print_graph_to_file_gnuplot<int>(
                file_prefix,
                file_name,
                plot_title,
                group_name_prefix,
                longitudes,
                latitudes,
                populations,
                this->assigned_to_team,
                labels
    );

}

void Mda::print_all_routes() const{
    for (size_t tt = 0; tt < this->routes.size(); tt++) {
        this->print_route(tt);
    }
    std::cout << "Villages (" << this->assigned_to_team.size() << ") assigned to teams: ";
    for (const auto& tti: this->assigned_to_team) {std::cout << tti << " ";}
    std::cout << std::endl;
}
void Mda::print_route(int route_id) const{
    std::cout << "MDA Team#" << route_id << " Visit Route (" << this->routes.at(route_id).size() << " villages): ";
    for (const auto& vvi: this->routes.at(route_id)) {std::cout << this->village_ids[vvi] << " ";}
    std::cout << std::endl;
}


void Mda::register_start_time(int time_step) {
    this->scheduled_start_times.push_back(time_step);
}
int Mda::get_next_start_time() const{
    if (num_executed_mdas >= this->scheduled_start_times.size()){
        return -1;
    } else {
        return this->scheduled_start_times.at(this->num_executed_mdas);
    }
}

void Mda::step(int time_step, village::VillageManager* vll_mgr){
    if(time_step == this->get_next_start_time()) {
        
        std::cout << "Start scheduled MDA #"
                  << this->num_executed_mdas+1 << "/"
                  << this->scheduled_start_times.size()
                  << " on day #" << time_step << std::endl;

        this->start();
        // std::cin.ignore();
    }

    if(this->in_execution) {
        this->execute(vll_mgr);
    }
}

void Mda::start(){
    this->in_execution = true;
    this->execution_step = 0;
    this->route_step = 0;

    // this->visited_in_execution.clear();
    // this->visited_in_execution.resize(this->assigned_to_team.size(), false);

    this->visited_days_in_execution.clear();
    this->visited_days_in_execution.resize(this->assigned_to_team.size(), 0);

    this->num_executed_mdas++; //TODO this should be called started_mdas
}

void Mda::execute(village::VillageManager* vll_mgr) {

    bool exist_unfinished_team = false;

    for (uint tt = 0; tt < this->routes.size(); tt++){

        // if ( (uint)this->route_step < this->routes[tt].size() ) {

        //     vll_mgr->treatment_by_mda(
        //         this->routes[tt].at(this->route_step),
        //         this->kMda_drug
        //     );

        //     visited_in_execution[this->routes[tt].at(this->route_step)] = true;
        //     exist_unfinished_team = true;
        // }

        for (int vs = 0; vs < this->kVisits_per_village; vs++) { // vs: visit
            int rs = this->route_step - vs*this->kRoute_step_offset;// rs: route-step
            if ( rs >= 0 && rs < static_cast<int>(this->routes[tt].size())){
                vll_mgr->treatment_by_mda(
                    // this->routes[tt].at(rs),
                    this->village_ids[this->routes[tt].at(rs)],
                    this->kMda_drug
                );

                std::cout << "Team " << tt
                        << " treating village" << this->village_ids[this->routes[tt].at(rs)]
                        << "\n";

                visited_days_in_execution[this->routes[tt].at(rs)]++;
                exist_unfinished_team = true;

            }
            if(rs < 0){exist_unfinished_team = true;}
        }

    }

    if (exist_unfinished_team) {
        this->execution_step++;

        if (this->execution_step % this->kDays_per_village == 0) {
            this->route_step++;
        }

    } else {
        this->finish();
    }
    
}

void Mda::finish(){
    // bool all_villages_were_visited_in_execution = true;
    // for (bool bb: this->visited_in_execution){
    //     if(!bb){
    //         all_villages_were_visited_in_execution = false;
    //         break;
    //     } 
    // }
    // assert(all_villages_were_visited_in_execution);

    int days_sum = 0;
    for (int dd: this->visited_days_in_execution){
        days_sum += dd;
    }
    assert(days_sum == this->kDays_per_village
                        *this->kVisits_per_village
                        *static_cast<int>(this->assigned_to_team.size())
        );

    std::cout << "MDA # "
                << this->num_executed_mdas << "/"
                << this->scheduled_start_times.size()
                << " finished.\n";
    // std::cin.ignore();


    this->in_execution = false;

}



}