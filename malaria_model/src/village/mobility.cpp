#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
// #include <chrono>

// #include "common/common.h" // village.h
#include "village/village.h"
#include "human/human.h"
#include "util/randomness.h"
#include "util/graph.h"

namespace village {

void VillageManager::init_shortest_path_distances(const util::VillageGraph &vg){
        
    for (int vv1 = 0; vv1 < this->sum_num_villages; vv1++) {
        std::vector<float> dists = vg.find_shortest_paths(vv1);
        for (int vv2 = 0; vv2 < this->sum_num_villages; vv2++) {
            this->distances_shortest_path[vv1][vv2] = dists[vv2];
            // std::cout << vv1 << " to " << vv2 << " is " << dists[vv2] << "\n";
        }
    }
}


std::string VillageManager::output_shortest_path_from_village_gnuplot(
        const std::string file_prefix,
        const int village_id,
        const util::VillageGraph &vg
    ) const {
    std::string file_name = "path_from_" + std::to_string(village_id);
    std::string plot_title = "Travel from Village " + std::to_string(village_id);
    std::string group_name_prefix = "Distances from village";
    std::vector<int> group_data(this->sum_num_villages, village_id);
    std::vector<int> labels;
    labels.push_back(village_id);

    std::cout << "labels size:" << labels.size() << "\n";

    return vg.print_graph_to_file_gnuplot<float>(
                file_prefix,
                file_name,
                plot_title,
                group_name_prefix,
                this->longitude_of,
                this->latitude_of,
                this->distances_shortest_path[village_id],
                group_data,
                labels
    );
}

void VillageManager::build_static_mobility_network(){
    float G_move = 0.95;
    float crit = 1.5;
    // float** prob = new float*[rows];
    // int** nprob = new int*[rows];
    int tot = 5000;

    // for (d=0; d<G_xregion; d++){
    for (int vv1 = 0; vv1 < this->sum_num_villages; vv1++){
         // village=G_villages[d];
         // lat=G_villages[d]->latitude;
         // longit=G_villages[d]->longitude;
         //cout<<"tvill: "<<G_villages[d]->latitude<<endl;
         // prob[d]=new float[rows];
         std::vector<float> prob_vv1(this->sum_num_villages, 0.0);
         // float sump[G_xregion];
         float sump = 0;
         float eps =0.05;
         // for (i=0; i<G_xregion; i++){
         for (int vv2 = 0; vv2 < this->sum_num_villages; vv2++){
            // village2=G_villages[i];
            // lat2=G_villages[i]->latitude;
            // longit2=G_villages[i]->longitude;

                  // if (d==i){
                  if (vv1 == vv2){
                    // prob[d][i]=0;
                    prob_vv1[vv2] = 0;
                  }
                  else{
                      // double dist = distanceEarth(lat,longit,lat2,longit2);
                      float dist = this->distances[vv1][vv2];
                      // prob[d][i]=float(village2->npeople)*float(G_villages[i]->npeople)/((1.0+(1.0/(1.0+float(exp(10*(float(dist-crit)))))))*(float(pow(dist,.5))+eps));
                      prob_vv1[vv2]=float(this->population_of[vv2]*this->population_of[vv1])/((1.0+(1.0/(1.0+float(exp(10*(float(dist-crit)))))))*(float(pow(dist,.5))+eps));
                  }
                  // if (i==0){
                  if (vv2==0){
                     // sump[i]=prob[d][i];
                     sump = prob_vv1[vv2];
                  }
                  else{
                    // sump[i]=sump[i-1]+prob[d][i];
                    sump += prob_vv1[vv2];
                  }
         }
         // prob[d][d]=(sump[G_xregion-1]*G_move);
         prob_vv1[vv1] = sump*G_move;
         // float nsump=sump[G_xregion-1]+prob[d][d];
         float nsump = sump + prob_vv1[vv1];

         int k=0;
         // nprob[d]=new int[rows];
         std::vector<int> nprob_vv1(this->sum_num_villages,0);
         // for (int w=0; w<G_xregion; w++){
         for (int vv3 = 0; vv3 < this->sum_num_villages; vv3++){
              // prob[d][w]=prob[d][w]/nsump;
              prob_vv1[vv3]=prob_vv1[vv3]/nsump;
              // nprob[d][w]=int(prob[d][w]*tot);
              nprob_vv1[vv3]=int(prob_vv1[vv3]*tot);
              // k=k+nprob[d][w];
              k = k + nprob_vv1[vv3];
              if (k>tot){
                 // nprob[d][w]=tot-nprob[d][w-1];
                 nprob_vv1[vv3] = tot - nprob_vv1[vv3-1];
              }
              // for (int q=0;q<nprob[d][w];q++){
              for (int q=0;q < nprob_vv1[vv3];q++){
                  // village->staticnet.push_back(w);
                  this->static_mobility_network[vv1].push_back(vv3);
              }
         }
        //cout<< "region: "<<d<<" size: "<<  village->staticnet.size()<<endl;
    }

    // for (int vv = 0; vv < this->sum_num_villages; vv++) {
    //   std::cout << "v" << vv << "(" << static_mobility_network[vv].size() << ")" << ": ";
    //   // for (uint ii = 0; ii < this->static_mobility_network[vv].size(); ii++) {
    //   //   std::cout << this->static_mobility_network[vv][ii] << ", ";
    //   // }
    //   std::cout << "\n";
    // }
    // std::cin.ignore();
}

void VillageManager::step_population_movement_static(
        const int human_array_size,
        const common::HumanMobilityType* human_mobility_type_array,
        const int* human_home_village_array,

        const float* random_numbers_array,
        const int random_numbers_array_size
    ){

    if (this->kMobility_enabled) {
        
        this->step_population_movement_static_kernel(
            this->sum_num_villages,
            this->at_village_register,
            this->static_mobility_network,

            human_array_size,
            human_mobility_type_array,
            human_home_village_array,

            this->kStatic_population_move_out_probability,
            this->kStatic_population_return_home_probability,
            this->kNon_static_population_move_out_probability,
            this->kNon_static_population_return_home_probability,

            random_numbers_array,
            random_numbers_array_size
        );
        // this->villager_reg_clean_up(this->sum_num_villages,this->at_village_register);
    
    }


    this->step_update_at_village_attractiveness();
}

void VillageManager::step_population_movement_static_kernel(
        const int village_array_size,
        // std::vector<std::vector<int>>& villager_reg,
        // boost::container::vector<boost::container::vector<int>>& villager_reg,
        Register_t& villager_reg,
        const std::vector<std::vector<int>>& mobility_network_static,

        const int human_array_size,
        const common::HumanMobilityType* human_mobility_type_array,
        const int* human_home_village_array,

        const float static_population_move_out_rate,    // G_mobility_static
        const float static_population_return_home_rate, // G_static_return
        const float non_static_population_move_out_rate,   // G_mobility_mobile
        const float non_static_population_return_home_rate, // G_mobile_return

        const float* random_numbers_array,
        const int random_numbers_array_size // 2 * human_array_size
    ){

    assert(random_numbers_array_size == human_array_size * 2);

    const int kValue_for_erase = -1;

    int num_moves = 0;

    const float kMove_out_rate = static_population_move_out_rate;
    const float kReturn_home_rate = static_population_return_home_rate;
    (void)non_static_population_move_out_rate;
    (void)non_static_population_return_home_rate;
    (void)human_mobility_type_array;
    // (void)mobility_network_static;

    const float kAny_movement_rate = std::max(kMove_out_rate,kReturn_home_rate);

    // std::vector<bool> human_already_moved_in_this_step(human_array_size, false);
    static std::vector<bool> human_already_moved_in_this_step(human_array_size, false);
    std::fill(human_already_moved_in_this_step.begin(), human_already_moved_in_this_step.end(), false);

    float rndm_nmbr_if_move = 0.0;
    float rndm_nmbr_move_to = 0.0;
    // bool erase_executed_with_iterator = false;

// static std::chrono::duration<double, std::milli> time_move_out;
// static std::chrono::duration<double, std::milli> time_return;
// static std::chrono::duration<double, std::milli> time_all;
// static std::chrono::duration<double, std::milli> time_target;
// static std::chrono::duration<double, std::milli> time_target_mid;
// static std::chrono::duration<double, std::milli> time_erase;

    // std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();

    for (int vv = 0; vv < village_array_size; vv++) {
        // std::cout << "village " << vv << std::endl;
        // for (auto hh_it : villager_reg[vv]){
        //     std::cout << hh_it << ", ";
        // }
        // std::cout << std::endl;

        // for (auto hh_it = villager_reg[vv].begin(); hh_it != villager_reg[vv].end(); ++hh_it;){
        // auto - std::vector<int>::iterator
        const std::vector<int>& vv_static_network = mobility_network_static[vv];
        const int vv_static_network_size = static_cast<int>(vv_static_network.size());

        for (auto& hh : villager_reg[vv]) {

            // erase_executed_with_iterator = false;
            
            // std::cout << "*hh_it:" << *hh_it << std::endl;
            // std::cin.ignore();

            assert(hh < human_array_size);
            assert(hh >= 0);

            if (human_already_moved_in_this_step[hh]) {
                // ++hh_it;
                continue;
            }
            
            // util::set_random_numbers_uniform(&rndm_nmbr, 1);
            rndm_nmbr_if_move = random_numbers_array[hh];

            if (rndm_nmbr_if_move >= kAny_movement_rate) {
                continue;
            }


            // std::cout << "home:" << human_home_village_array[*hh_it] << "\n";

            // 1. where is this human in relation to his home?
            if (human_home_village_array[hh] == vv) {
                
                // 1.1 at home >> move out?

                // if (human_mobility_type_array[hh] == common::HumanMobilityType::kStatic) {
                //     move_out_rate = static_population_move_out_rate;
                // } else {
                //     move_out_rate = non_static_population_move_out_rate;
                // }

                // move_out_rate = static_population_move_out_rate;

                if (rndm_nmbr_if_move < kMove_out_rate) {

                    // 1.1.1 execute move out


                    // util::set_random_numbers_uniform(&rndm_nmbr, 1);
                    rndm_nmbr_move_to = random_numbers_array[hh + human_array_size];

                    // (void)mobility_network_static;
    // std::chrono::high_resolution_clock::time_point t_start_target = std::chrono::high_resolution_clock::now();

                    // int target_village = mobility_network_static[vv]
                    //     [static_cast<int>(rndm_nmbr_move_to * mobility_network_static[vv].size())];
                    // int target_village = mobility_network_static.at(vv)
                    //     .at(static_cast<int>(rndm_nmbr_move_to * mobility_network_static[vv].size()));

                    // int target_village = rndm_nmbr_move_to * village_array_size;
                    // int target_index = rndm_nmbr_move_to * mobility_network_static[vv].size();
                    int target_index = rndm_nmbr_move_to * vv_static_network_size;
    // std::chrono::high_resolution_clock::time_point t_mid_target = std::chrono::high_resolution_clock::now();

                    // int target_village = mobility_network_static[vv][target_index];
                    int target_village = vv_static_network[target_index];
                    

    // std::chrono::high_resolution_clock::time_point t_end_target = std::chrono::high_resolution_clock::now();
    // time_target += t_end_target - t_start_target;
    // time_target_mid += t_mid_target - t_start_target;


                    // int target_village = static_cast<int>(rndm_nmbr_move_to * village_array_size);

                    // std::cout << "move out to " << target_village << "\n";

                    if (target_village == vv) {// moved to home village
                        continue;
                    } else {
    // std::chrono::high_resolution_clock::time_point t_start_move_out = std::chrono::high_resolution_clock::now();
                        villager_reg[target_village].push_back(hh);
                        this->v_hmn_mgr->set_human_at_village(hh, target_village);
    // std::chrono::high_resolution_clock::time_point t_end_move_out = std::chrono::high_resolution_clock::now();
                        human_already_moved_in_this_step[hh] = true;
                        hh = kValue_for_erase;

                        num_moves++;
    // time_move_out += t_end_move_out - t_start_move_out;
                    }

    // time_move_out += std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count();
                    // hh_it = villager_reg[vv].erase(hh_it);

                    // erase_executed_with_iterator = true;
                }

            } else {
                
                // 1.2 away from home >> return home ?

                // if (human_mobility_type_array[hh] == common::HumanMobilityType::kStatic) {
                //     return_home_rate = static_population_return_home_rate;
                // } else {
                //     return_home_rate = non_static_population_return_home_rate;
                // }

                // return_home_rate = static_population_return_home_rate;

                if (rndm_nmbr_if_move < kReturn_home_rate) {

                    // 1.2.1 execute return home

                    // int home_village = human_home_village_array[hh];
                    // std::cout << "return home to " << home_village << "\n";
    // std::chrono::high_resolution_clock::time_point t_start_return = std::chrono::high_resolution_clock::now();

                    villager_reg[human_home_village_array[hh]].push_back(hh);
                    this->v_hmn_mgr->set_human_return_home_village(hh);
    // std::chrono::high_resolution_clock::time_point t_end_return = std::chrono::high_resolution_clock::now();
                    human_already_moved_in_this_step[hh] = true;

                    hh = kValue_for_erase;

                    num_moves++;
    // time_return += t_end_return - t_start_return;


                    // hh_it = villager_reg[vv].erase(hh_it);

                    // erase_executed_with_iterator = true;

                }

            }

        }
    // std::chrono::high_resolution_clock::time_point t_start_erase = std::chrono::high_resolution_clock::now();

        villager_reg[vv].erase(
            std::remove(
                villager_reg[vv].begin(),
                villager_reg[vv].end(),
                kValue_for_erase
            ),
            villager_reg[vv].end()
        );
    // std::chrono::high_resolution_clock::time_point t_end_erase = std::chrono::high_resolution_clock::now();
    // time_erase += t_end_erase - t_start_erase;

        // VillageManager::remove_villager(villager_reg, vv);
    }

    // std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();

    // time_all += t_end - t_start;


    // int total_population_check = 0;
    // for (int vv = 0; vv < village_array_size; vv++) {
    //     total_population_check += villager_reg[vv].size();
    // }
    // assert(total_population_check == human_array_size);
    // std::cout << "\ttotal_population_check=" << total_population_check << "\n";

    std::cout << "\tnum_static_movements=" << num_moves << "\n";
    // std::cout << "\t\ttime_all=" << time_all.count() << "\n";
    // std::cout << "\t\time_move_out=" << time_move_out.count() << "\n";
    // std::cout << "\t\time_return=" << time_return.count() << "\n";
    // std::cout << "\t\time_target=" << time_target.count() << "\n";
    // std::cout << "\t\time_target_mid=" << time_target_mid.count() << "\n";
    // std::cout << "\t\time_erase=" << time_erase.count() << "\n";
}

// void VillageManager::remove_villager(Register_t& villager_reg, int vv){

//   const int kValue_for_erase = -1;
//           villager_reg[vv].erase(
//             std::remove(
//                 villager_reg[vv].begin(),
//                 villager_reg[vv].end(),
//                 kValue_for_erase
//             ),
//             villager_reg[vv].end()
//         );

// }

// void VillageManager::villager_reg_clean_up(
//         const int village_array_size,
//         // std::vector<std::vector<int>>& villager_reg
//         // boost::container::vector<boost::container::vector<int>>& villager_reg
//         Register_t& villager_reg
//     ){
//     const int kValue_for_erase = -1;
//     for (int vv = 0; vv < village_array_size; vv++) {
//         villager_reg[vv].erase(
//             std::remove(
//                 villager_reg[vv].begin(),
//                 villager_reg[vv].end(),
//                 kValue_for_erase
//             ),
//             villager_reg[vv].end()
//         );
//     }
// }

// void VillageManager::step_population_movement_static_kernel(
//         const int village_array_size,
//         std::vector<std::vector<int>>& villager_reg,
//         const std::vector<std::vector<int>>& mobility_network_static,

//         const int human_array_size,
//         const common::HumanMobilityType* human_mobility_type_array,
//         const int* human_home_village_array,

//         const float static_population_move_out_rate,    // G_mobility_static
//         const float static_population_return_home_rate, // G_static_return
//         const float non_static_population_move_out_rate,   // G_mobility_mobile
//         const float non_static_population_return_home_rate, // G_mobile_return

//         const float* random_numbers_array,
//         const int random_numbers_array_size // 2 * human_array_size
//     ){

//     assert(random_numbers_array_size == human_array_size * 2);

//     float move_out_rate = 0.0;
//     float return_home_rate = 0.0;

//     std::vector<bool> human_already_moved_in_this_step(human_array_size, false);

//     float rndm_nmbr_if_move = 0.0;
//     float rndm_nmbr_move_to = 0.0;
//     bool erase_executed_with_iterator = false;

//     for (int vv = 0; vv < village_array_size; vv++) {
//         // std::cout << "village " << vv << std::endl;
//         // for (auto hh_it : villager_reg[vv]){
//         //     std::cout << hh_it << ", ";
//         // }
//         // std::cout << std::endl;
//         for (auto hh_it = villager_reg[vv].begin(); hh_it != villager_reg[vv].end(); ){
//         // auto - std::vector<int>::iterator
//             erase_executed_with_iterator = false;
            
//             // std::cout << "*hh_it:" << *hh_it << std::endl;
//             // std::cin.ignore();

//             assert(*hh_it < human_array_size);
//             assert(*hh_it >= 0);

//             if (human_already_moved_in_this_step[*hh_it]) {
//                 ++hh_it;
//                 continue;
//             }
            
//             // util::set_random_numbers_uniform(&rndm_nmbr, 1);
//             rndm_nmbr_if_move = random_numbers_array[*hh_it];


//             // std::cout << "home:" << human_home_village_array[*hh_it] << "\n";

//             // 1. where is this human in relation to his home?
//             if (human_home_village_array[*hh_it] == vv) {
                
//                 // 1.1 at home >> move out?

//                 if (human_mobility_type_array[*hh_it] == common::HumanMobilityType::kStatic) {
//                     move_out_rate = static_population_move_out_rate;
//                 } else {
//                     move_out_rate = non_static_population_move_out_rate;
//                 }

//                 move_out_rate = static_population_move_out_rate;

//                 if (rndm_nmbr_if_move < move_out_rate) {

//                     // 1.1.1 execute move out


//                     // util::set_random_numbers_uniform(&rndm_nmbr, 1);
//                     rndm_nmbr_move_to = random_numbers_array[*hh_it + human_array_size];

//                     int target_village = mobility_network_static[vv]
//                         [static_cast<int>(rndm_nmbr_move_to * mobility_network_static[vv].size())];

//                     // std::cout << "move out to " << target_village << "\n";

//                     if (target_village == vv) {
//                         ++hh_it;
//                         continue;
//                     }

//                     villager_reg[target_village].push_back(*hh_it);
//                     human_already_moved_in_this_step[*hh_it] = true;

//                     hh_it = villager_reg[vv].erase(hh_it);

//                     erase_executed_with_iterator = true;
//                 }

//             } else {
                
//                 // 1.2 away from home >> return home ?

//                 if (human_mobility_type_array[*hh_it] == common::HumanMobilityType::kStatic) {
//                     return_home_rate = static_population_return_home_rate;
//                 } else {
//                     return_home_rate = non_static_population_return_home_rate;
//                 }

//                 return_home_rate = static_population_return_home_rate;

//                 if (rndm_nmbr_if_move < return_home_rate) {

//                     // 1.2.1 execute return home

//                     int home_village = human_home_village_array[*hh_it];
//                     // std::cout << "return home to " << home_village << "\n";

//                     villager_reg[home_village].push_back(*hh_it);
//                     human_already_moved_in_this_step[*hh_it] = true;

//                     hh_it = villager_reg[vv].erase(hh_it);

//                     erase_executed_with_iterator = true;

//                 }

//             }

//             // 2. iterator movement required if iterator has not been updated via erase()
//             if (!erase_executed_with_iterator) {
//                 ++hh_it;
//             }

//             erase_executed_with_iterator = false;

//         }
//     }

// }

}