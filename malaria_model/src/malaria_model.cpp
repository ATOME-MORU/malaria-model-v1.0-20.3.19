#define DO_QUOTE(X)       #X
#define QUOTE(X)          DO_QUOTE(X)

#include <iostream>
// #include <vector> // util::InputParser, not if -std=c++11
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <string> //std::stoi
#include <fstream> //std::ifstream, std::getline
#include <algorithm>

#include <experimental/filesystem> // -lstdc++fs in LDLIBS
namespace fs = std::experimental::filesystem;

#include <iomanip> //std::put_time
#include <ctime>

#include <map>

#include <cassert>

// #include <thread>

#include "rapidjson/document.h"

#include "util/util.h"
#include "util/randomness.h"
#include "village/village.h"
#include "village/mosquito.h"
#include "village/mosquito_ibm_ra.h"

#include "human/human.h"
#include "human/blood_system.h"


#include "simulation/reporter.h"
#include "intervention/mda.h"
#include "intervention/fmda.h"
#include "intervention/smda.h"
#include "intervention/tmda.h"
#include "intervention/vc.h"

#include "simulation/simulator.h"

// #include "util/graph.h"

#include <stdio.h>

int main(int argc, char** argv){

    // Parse user command line input

    util::InputParser input_parser(argc, argv);

    if(input_parser.cmdOptionExists("--version")){

        std::cout << "Build (" << QUOTE(BUILD_DATE) << ") :\n";
        std::cout   << "\tHEAD : " << QUOTE(GIT_COMMIT) << "\n"
                    << "\tDate : " << QUOTE(GIT_DATE) << "\n";

        return EXIT_SUCCESS;
    }

    const std::string& config_file_name = input_parser.getCmdOption("-c");
    if(!config_file_name.empty()) {

        // if(!checkFileExists(config_file_name)) {
        // std::cout << fs::current_path()/config_file_name.c_str() << std::endl;
        if(!fs::exists(config_file_name.c_str())){
            std::cout << "File (" << config_file_name << ") does not exist." << std::endl;
            return EXIT_FAILURE;
        }

        std::cout << "Reading configuration from ";
        std::cout << config_file_name << "" << std::endl;

    } else {

        std::cout << "Please give name of configuration file (JSON) using the \"-c [file]\" commandline option.";
        std::cout << " (e.g. -c data/config.json)" << std::endl;
        // std::cout << "(" << filename << ")" << std::endl;
        return EXIT_FAILURE;
    }

    // if(!fs::exists(simulation::kJson_schema_config_file_name.c_str())){
    //     std::cout << "File (" << simulation::kJson_schema_config_file_name << ") does not exist." << std::endl;
    //     return EXIT_FAILURE;
    // }


    const std::string& output_folder_name = input_parser.getCmdOption("-o");
    if(!output_folder_name.empty()) {

        if(!fs::exists(output_folder_name.c_str())){
            if(!fs::create_directory(output_folder_name.c_str())) {
                std::cout << "Could not create directory at (" << output_folder_name << ")." << std::endl;
                return EXIT_FAILURE;
            }
        }

    } else {
        // TODO: give a default value and proceed

        std::cout << "Please give name of output directory using the \"-o [path]\" commandline option.";
        std::cout << " (e.g. -o outputs)" << std::endl;
        
        return EXIT_FAILURE;
    }


    // -i [prefix], optional, if not given a timestamp will be used to prefix this run 
    const std::string& output_prefix_name = input_parser.getCmdOption("-i");


#ifdef __NVCC__
    const std::string& selected_device = input_parser.getCmdOption("-d");
    if(!selected_device.empty()){
        util::cuda_device = std::stoi(selected_device);
    }
#endif

    if(!simulation::Simulator::validate_config(config_file_name)){
        std::cout << "Error in configuration file ("<< config_file_name <<")\n";
        return EXIT_FAILURE;

    } else {

        std::cout << "Schema check passed\n";
        // return EXIT_SUCCESS;

        simulation::Simulator sim(
            config_file_name,
            output_folder_name,
            output_prefix_name
        );


        sim.run();
        sim.close();
        
    }



    return EXIT_SUCCESS;



}


//     rapidjson::Document rj_doc_config = util::get_json_from_file(config_file_name);
//     // rapidjson::Document rj_doc_config_schema = util::get_json_from_file(simulation::kJson_schema_config_file_name);

//     // if (!util::validate_json_against_schema(&rj_doc_config_schema, &rj_doc_config)) {return EXIT_FAILURE;}

//     rapidjson::Value::MemberIterator rj_find_it = rj_doc_config["simulation"].FindMember("total_steps");
//     if (rj_find_it != rj_doc_config.MemberEnd()) {
//         std::cout << "found, value: " << rj_find_it->value.GetInt() <<"\n";
//     } else {
//         std::cout << "not found\n";
//         // std::cout << "default:" << rj_doc_config_schema["definitions"]["simulation"]["properties"]["total_steps"]["default"].GetString() << std::endl;
//     }

//     assert(rj_doc_config["simulation"]["total_steps"].IsInt());
//     // std::cout << rj_doc_config["simulation"]["total_steps"].GetInt() << std::endl;
//     const int kSimulation_steps_max = rj_doc_config["simulation"]["total_steps"].GetInt();
//     assert(rj_doc_config["human"]["male_percentage"].IsNumber());
//     assert(rj_doc_config["human"]["male_percentage"].IsFloat());
//     std::cout << rj_doc_config["human"]["male_percentage"].GetFloat() << std::endl;

//     assert(rj_doc_config["village"]["count"].IsInt());
//     const int kNum_villages = rj_doc_config["village"]["count"].GetInt();

//     assert(rj_doc_config["intervention"]["mda"]["num_teams"].IsInt());
//     const int kNum_mda_teams = rj_doc_config["intervention"]["mda"]["num_teams"].GetInt();

//     const rapidjson::Value& kMda_schedule = rj_doc_config["intervention"]["mda"]["schedule"];
//     assert(kMda_schedule.IsArray());
//     for(rapidjson::SizeType ii = 0; ii < kMda_schedule.Size(); ii++){
//         std::cout << "mda scheduled start time #" << ii << " is " << kMda_schedule[ii].GetInt() << std::endl;
//     }

//     // return EXIT_SUCCESS;


//     // const std::string& output_folder_name = input_parser.getCmdOption("-o");
//     // if(!output_folder_name.empty()) {

//     //     if(!fs::exists(output_folder_name.c_str())){
//     //         if(!fs::create_directory(output_folder_name.c_str())) {
//     //             std::cout << "Could not create directory at (" << output_folder_name << ")." << std::endl;
//     //             return EXIT_FAILURE;
//     //         }
//     //     }

//     // } else {
//     //     // TODO: give a default value and proceed

//     //     std::cout << "Please give name of output directory using the \"-o [path]\" commandline option.";
//     //     std::cout << " (e.g. -o outputs)" << std::endl;
        
//     //     return EXIT_FAILURE;
//     // }

//     std::string output_prefix = util::get_output_prefix();
//     std::cout << output_prefix << std::endl;



//     // Read configurations

//     boost::property_tree::ptree pt;//TODO: replace property_tree with rapidjson
//     boost::property_tree::json_parser::read_json(config_file_name, pt);

//     char delimiter_used_in_file = 't';

//     try {
//         delimiter_used_in_file = pt.get<char>("village.data_file.delimiter");
//     } catch (std::exception & e) {
//         std::cout << "village.data_file.delimiter not given by user, using default"
//                   << "(" << delimiter_used_in_file << ")" <<std::endl;
//     }


//     switch (delimiter_used_in_file) {
//         case 's': 
//                     delimiter_used_in_file = ' '; //TODO
//                     break;
//         case 't': 
//                     delimiter_used_in_file = '\t';
//                     break;
//     }



//     // Initialise simulation data

//     village::VillageManager vll_mgr(
//             pt.get<std::string>("village.data_file.name"),
//             delimiter_used_in_file,
//             std::stoi(pt.get<std::string>("village.count")),
//             kSimulation_steps_max
//         );

//     assert(rj_doc_config["mosquito"]["seasonality"]["data_file"]["name"].IsString());
//     vll_mgr.init_mosquito_seasonality_from_file(rj_doc_config["mosquito"]["seasonality"]["data_file"]["name"].GetString());

//     vll_mgr.print_all();
//     vll_mgr.init_treatment_rates_malaria_post();
//     vll_mgr.print_summary();


//     intervention::Mda mda(kNum_mda_teams, kNum_villages, vll_mgr.get_village_distances());
//     std::string mda_print_map_sys_command = mda.print_map(
//         output_folder_name + util::kPathSeparator + output_prefix,
//         vll_mgr.longitude_of, vll_mgr.latitude_of, vll_mgr.population_of
//         ) + " &";
//     std::cout << mda_print_map_sys_command << std::endl;
//     if(system(mda_print_map_sys_command.c_str())){};

//     for(rapidjson::SizeType ii = 0; ii < kMda_schedule.Size(); ii++){
//         mda.register_start_time(kMda_schedule[ii].GetInt());
//     }
//     mda.print_all_routes();


//     vll_mgr.init_shortest_path_distances(mda.get_graph());
//     std::string print_distance_sys_command = vll_mgr.output_shortest_path_from_village_gnuplot(
//         output_folder_name + util::kPathSeparator + output_prefix,
//         21,
//         mda.get_graph()
//     );
//     std::cout << print_distance_sys_command << std::endl;
//     if(system(print_distance_sys_command.c_str())){};

//     // return EXIT_SUCCESS;

//     // util::VillageGraph vg(vll_mgr.sum_num_villages, vll_mgr.get_village_distances());
    
//     // // vg.update_edge_weight_stats(false);

//     // std::cout << "Reducing graph to mininum spanning tree..." << std::flush;
//     // vg.mark_graph_with_min_spanning_tree();
//     // vg.reduce_graph_to_min_spanning_tree();
//     // std::cout << "[done]\n";

//     // std::string village_connections_gnuplot_command =
//     //     vg.print_graph_to_file_gnuplot(
//     //         output_folder_name + util::kPathSeparator + output_prefix,
//     //         vll_mgr.longitude_of, vll_mgr.latitude_of
//     //     ) + " &";
//     // std::cout << village_connections_gnuplot_command << std::endl;
//     // system(village_connections_gnuplot_command.c_str());

//     // // std::string vg_file_name = vg.print_graph_to_file_dot(output_folder_name + util::kPathSeparator + output_prefix  + "_mst_");
//     // // std::string xdot_command = "xdot " + vg_file_name + " &";
//     // // system(xdot_command.c_str());

//     // float* random_numbers_teams = new float[kNum_mda_teams];
//     // util::set_random_numbers_uniform(random_numbers_teams, kNum_mda_teams);

//     // for (int tt = 0; tt < kNum_mda_teams; tt++){
//     //     int start_village = (int)(random_numbers_teams[tt]*vll_mgr.sum_num_villages);
//     //     std::cout << "Invervention Team " << tt << " route: ";
//     //     std::vector<int> bfs = vg.get_bfs(start_village);
//     //     for (const auto& vvi: bfs){std::cout << vvi << " ";} 
//     //     std::cout << std::endl;
//     // }

//     // return EXIT_SUCCESS;


//     human::HumanManager hmn_mgr(vll_mgr.sum_total_population);
//     // hmn_mgr.print_all();
    
//     vll_mgr.init_get_villagers(&hmn_mgr);

//     // vll_mgr.give_drug_to_village(0, common::DrugName::kArtesunate);


//     hmn_mgr.init_mobility(vll_mgr.get_migrant_population_percentages());
//     // hmn_mgr.print_all();

//     hmn_mgr.read_age_distribution(pt.get<std::string>("human.age_weight_vector_file"));

//     hmn_mgr.init_age();
//     // hmn_mgr.print_age_distribution();

//     hmn_mgr.init_gender(std::stof(pt.get<std::string>("human.male_percentage")));
//     hmn_mgr.print_gender_distribution();


//     assert(rj_doc_config["human"]["itn_probability"].IsNumber());
//     assert(rj_doc_config["human"]["itn_probability"].IsFloat());
//     hmn_mgr.init_itn(rj_doc_config["human"]["itn_probability"].GetFloat());

//     hmn_mgr.print_all();

//     human::BloodSystemManager bsys_mgr(hmn_mgr.sum_num_humans);
//     // if (true) {
//         assert(rj_doc_config["within_host"]["init_file"]["name"].IsString());

//         bsys_mgr.init_from_file(rj_doc_config["within_host"]["init_file"]["name"].GetString(),'\t');

//     // }

//     bsys_mgr.update_susceptibility(hmn_mgr.get_if_itn_of(), hmn_mgr.sum_num_humans);


    
//     vll_mgr.init_infections(hmn_mgr.get_bitten_by(), hmn_mgr.sum_num_humans);
//     bsys_mgr.parasite_intake(hmn_mgr.get_bitten_by(), common::StageName::kInfectious);

//     bsys_mgr.parasite_dominance();

//     bsys_mgr.print_all();

//     simulation::ReportManager rp_mgr(
//                                 output_folder_name + util::kPathSeparator + output_prefix,
//                                 kSimulation_steps_max,
//                                 &vll_mgr,
//                                 &hmn_mgr,
//                                 &bsys_mgr
//                             );
//     rp_mgr.open_ofstreams();
//     // rp_mgr.end_of_day_reports_print_headings();

//     std::cout << "parasite_count:";
//     printf("%lu\n",bsys_mgr.get_parasite_count());
//     std::cout << "parasite group count: " << bsys_mgr.sum_num_pgs << std::endl;
    
//     // Humans <-> BloodSystems 
//     bool* birth_death_reset_flags = new bool[hmn_mgr.sum_num_humans];
//     std::fill_n(birth_death_reset_flags, hmn_mgr.sum_num_humans, false);

//     float* random_numbers_movement = new float[hmn_mgr.sum_num_humans*2];
//     float* random_numbers_drug_loss = new float[bsys_mgr.sum_num_systems*2];

//     int random_number_repository_sys2_size = bsys_mgr.sum_num_systems*2;
//     float* random_number_repository_sys2 = new float[random_number_repository_sys2_size];

//     // Start Simulation

//     bool if_threading = false;
//     std::thread thrd_rnd_nmbrs_drug_loss;

//     double wall_time_begin = util::get_wall_time();

//     for (int tt = 0; tt <= kSimulation_steps_max; tt++){

//         std::cout << "[Step.begin]Day #" << tt << "\n";

//             hmn_mgr.clear_given_drug();

//         std::cout << "[V.1] Statistics ...\n";
//             // at the end
//         std::cout << "[V.2] Interventions (Static) ...\n";

//             mda.step(tt, &vll_mgr);


//         std::cout << "[V.3] Interventions (Dynamic) ...\n";

//             vll_mgr.open_malaria_post(tt);
//             vll_mgr.treatment_by_malaria_post(&bsys_mgr, &hmn_mgr);

//         std::cout << "[H.1.1] Birth/Death ...\n";
//             std::fill_n(birth_death_reset_flags, hmn_mgr.sum_num_humans, false);
//             int num_birth_and_death = hmn_mgr.birth_and_death(birth_death_reset_flags);

//             bsys_mgr.birth_and_death(birth_death_reset_flags, hmn_mgr.sum_num_humans);

//             std::cout << num_birth_and_death << std::endl;
//         // if (num_birth_and_death > 0) std::cin >> num_birth_and_death;

//         std::cout << "[H.1.2] Ageing ...\n";
//             if ( (tt+1) % common::kNum_days_in_one_year == 0 ){
//                 hmn_mgr.increase_all_ages_by_one();
//             }

//         std::cout << "[H.2] Mobility ...\n";

//             util::set_random_numbers_uniform(random_numbers_movement, hmn_mgr.sum_num_humans*2);

//             vll_mgr.step_population_movement_static(
//                 hmn_mgr.sum_num_humans,
//                 hmn_mgr.get_mobility_of(),
//                 hmn_mgr.get_home_village_of(),

//                 random_numbers_movement,
//                 hmn_mgr.sum_num_humans*2
//             );

//         std::cout << "[H.3] Susceptibility ...\n"; // H.7
//         std::cout << "[H.4] Immunity Loss ...\n";
//             util::set_random_numbers_uniform(
//                 random_number_repository_sys2, 
//                 random_number_repository_sys2_size/2
//             );
//             bsys_mgr.step_immunity_loss(
//                 random_number_repository_sys2,
//                 random_number_repository_sys2_size/2
//             );
        
//         std::cout << "[H.5] Drug Intake ...\n";

//             // hmn_mgr.print_distribution<common::DrugName>(hmn_mgr.get_given_drug(), hmn_mgr.sum_num_humans);
//             bsys_mgr.drug_intake(hmn_mgr.get_given_drug());


//         std::cout << "[P.1-3] Parasite Density (Replication/Decay/PD) ...\n";
//             bsys_mgr.drug_effect();

//         std::cout << "[P.4] Parasites Stage Progression ...\n";
//             bsys_mgr.parasite_stage_progression();

//             bsys_mgr.clinical_resolution();
            
//         std::cout << "[P.5] Parasite Dominance ...\n";
//             bsys_mgr.parasite_dominance();

//         std::cout << "[H.6] Drug Loss (PK) ...\n";
//             if (tt == 0 && if_threading){
//                 // util::set_random_numbers_uniform(random_numbers_drug_loss, bsys_mgr.sum_num_systems * 2);
//                 thrd_rnd_nmbrs_drug_loss = std::thread(
//                     util::set_random_numbers_uniform_boost_thread, 
//                     random_numbers_drug_loss,
//                     bsys_mgr.sum_num_systems * 2
//                 );
//             }
            
//             if (if_threading) {
//                 thrd_rnd_nmbrs_drug_loss.join();
//             } else {
//                 util::set_random_numbers_uniform(random_numbers_drug_loss, bsys_mgr.sum_num_systems * 2);
//             }

//             bsys_mgr.drug_loss(random_numbers_drug_loss);

//             if (if_threading) {
//                 thrd_rnd_nmbrs_drug_loss = std::thread(
//                     util::set_random_numbers_uniform_boost_thread, 
//                     random_numbers_drug_loss,
//                     bsys_mgr.sum_num_systems * 2
//                 );
//             }

//         std::cout << "[M] Mosquito Dynamics ...\n" ;
//             vll_mgr.step_infection_human_to_mosquito(
//                 bsys_mgr.get_infectiousness(),
//                 bsys_mgr.get_dominant_type(),
//                 bsys_mgr.sum_num_systems,
//                 tt
//             );
//             vll_mgr.mosquito_survival_step();
//             vll_mgr.mosquito_incubation_step();

//         std::cout << "[H.7] Parasite Intake (Infection from bites) ...\n";

//             bsys_mgr.update_susceptibility(hmn_mgr.get_if_itn_of(), hmn_mgr.sum_num_humans);

//             vll_mgr.infection_step(
//                 hmn_mgr.get_bitten_by(),
//                 hmn_mgr.sum_num_humans,
//                 bsys_mgr.get_susceptibility(),
//                 bsys_mgr.sum_num_systems,
//                 bsys_mgr.get_most_advanced_stage(),
//                 tt
//             );
//             // hmn_mgr.print_distribution<common::ParasiteType>(hmn_mgr.get_bitten_by(), hmn_mgr.sum_num_humans);
//             bsys_mgr.parasite_intake(hmn_mgr.get_bitten_by(), common::StageName::kLiver);

//         std::cout << "[Step.end] ...\n";

//         rp_mgr.end_of_day_reports(tt);

//         // std::cout << "Total parasite_count:";
//         // printf("%lu",bsys_mgr.get_parasite_count());
//         // std::cout << std::endl;

//         // std::cin.ignore();

//     }

//     double wall_time_end = util::get_wall_time();

//     std::cout << "Time cost: " << wall_time_end - wall_time_begin
//                 << " for " << kSimulation_steps_max << " steps, "
//                 << ( wall_time_end - wall_time_begin )/kSimulation_steps_max 
//                 << " per step." << std::endl;

//     delete[] random_numbers_drug_loss;
//     delete[] random_numbers_movement;
//     // bsys_mgr.print_all();

//     rp_mgr.close_ofstreams();
    




//     int pg_io_counter_size = *std::max_element(
//                                 bsys_mgr.debug_bites_accumulator_add.begin(),
//                                 bsys_mgr.debug_bites_accumulator_add.end()
//                             )+1;
//     int pg_io_counter_add_max_index = std::max_element(
//                                 bsys_mgr.debug_bites_accumulator_add.begin(),
//                                 bsys_mgr.debug_bites_accumulator_add.end()
//                             )-bsys_mgr.debug_bites_accumulator_add.begin();
//     std::cout << "max add:" << pg_io_counter_size-1
//                 << " (host " << pg_io_counter_add_max_index
//                 << ": c_exp=" << int(bsys_mgr.get_cummulative_exposures_of()[pg_io_counter_add_max_index])
//                 << ")\n";
//     int* pg_io_counter = new int[pg_io_counter_size];
//     std::fill_n(pg_io_counter, pg_io_counter_size, 0);
//     int pg_io_counter_sum = 0;


//     for (uint ii = 0; ii < bsys_mgr.debug_bites_accumulator_add.size(); ii++) {
//         pg_io_counter[bsys_mgr.debug_bites_accumulator_add.at(ii)]++;
//         pg_io_counter_sum += bsys_mgr.debug_bites_accumulator_add.at(ii);
//     }
//     std::cout << "sum add: " << pg_io_counter_sum << "\n";

//     for (int cc = 0; cc < pg_io_counter_size; cc++) {
//         std::cout << std::fixed << std::setprecision(1)
//                     << cc << ": "
//                     // << std::string(pg_io_counter[cc]*100/pg_io_counter_sum, '*')
//                     << " " << pg_io_counter[cc] << "\n";
//     }

//     delete[] pg_io_counter;









//     int pg_io_counter_less_size = -*std::min_element(
//                                 bsys_mgr.debug_bites_accumulator_less.begin(),
//                                 bsys_mgr.debug_bites_accumulator_less.end()
//                             )+1;
//     std::cout << "max less:-" << pg_io_counter_less_size-1 << "\n";
//     int* pg_io_counter_less = new int[pg_io_counter_less_size];
//     std::fill_n(pg_io_counter_less, pg_io_counter_less_size, 0);
//     int pg_io_counter_less_sum = 0;


//     for (uint ii = 0; ii < bsys_mgr.debug_bites_accumulator_less.size(); ii++) {
//         pg_io_counter_less[-bsys_mgr.debug_bites_accumulator_less.at(ii)]++;
//         pg_io_counter_less_sum += bsys_mgr.debug_bites_accumulator_less.at(ii);
//     }
//     std::cout << "sum less: " << pg_io_counter_less_sum << "\n";

//     for (int cc = 0; cc < pg_io_counter_less_size; cc++) {
//         std::cout << std::fixed << std::setprecision(1)
//                     << "-" << cc << ": "
//                     // << std::string(pg_io_counter_less[cc]*100/pg_io_counter_less_sum, '*')
//                     << " " << pg_io_counter_less[cc] << "\n";
//     }

//     delete[] pg_io_counter_less;










//     std::vector<int> bites_accumulator_both(bsys_mgr.debug_bites_accumulator_less.size(), 0);
//     for (uint ii = 0; ii < bites_accumulator_both.size(); ii++) {
//         bites_accumulator_both[ii] = bsys_mgr.debug_bites_accumulator_add[ii] + bsys_mgr.debug_bites_accumulator_less[ii];
//     }

//     int pg_io_counter_both_size = *std::max_element(
//                                 bites_accumulator_both.begin(),
//                                 bites_accumulator_both.end()
//                             )+1;
//     std::cout << "max:" << pg_io_counter_both_size-1 << "\n";
//     int* pg_io_counter_both = new int[pg_io_counter_both_size];
//     std::fill_n(pg_io_counter_both, pg_io_counter_both_size, 0);
//     int pg_io_counter_both_sum = 0;


//     for (uint ii = 0; ii < bites_accumulator_both.size(); ii++) {
//         pg_io_counter_both[bites_accumulator_both.at(ii)]++;
//         pg_io_counter_both_sum += bites_accumulator_both.at(ii);
//     }
//     std::cout << "sum: " << pg_io_counter_both_sum << "\n";

//     for (int cc = 0; cc < pg_io_counter_both_size; cc++) {
//         std::cout << std::fixed << std::setprecision(1)
//                     << cc << ": "
//                     // << std::string(pg_io_counter_both[cc]*100/pg_io_counter_both_sum, '*')
//                     << " " << pg_io_counter_both[cc] << "\n";
//     }

//     delete[] pg_io_counter_both;



//     // std::cout << "Number of PG_q intakes: " << bsys_mgr.sum_num_pg_intakes << "\n";
//     std::cout << "Number of ignored PG_q full incidents: " << bsys_mgr.sum_num_pg_q_is_full_incidents
//                 << "(" << float(bsys_mgr.sum_num_pg_q_is_full_incidents)/bsys_mgr.sum_num_pg_intakes*100
//                 << "% " << "of " << bsys_mgr.sum_num_pg_intakes << " total pg intakes)\n";


//     hmn_mgr.Foo();
//     vll_mgr.Foo();

//     //Some Soft Checks
//     std::cout << "\033[4mSoft Checks:\033[0m\n";

//     // 1. The number of people in all three groups are the same
//     std::cout << "1. Human.Population(" << hmn_mgr.sum_num_humans
//                 <<") == Village.Population(" << vll_mgr.sum_total_population
//                 <<") == BloodSystem.Count(" << bsys_mgr.sum_num_systems << ") ";
//     if (hmn_mgr.sum_num_humans == vll_mgr.sum_total_population && vll_mgr.sum_total_population == bsys_mgr.sum_num_systems) {
//         std::cout << "\033[32;47m[Pass]\033[0m";
//     } else {
//         std::cout << "\033[34;47m[Fail]\033[0m";
//     }
//     std::cout << "\n";

//     // 2. No one lives/is located at unknown village
//     int home_max = *std::max_element(hmn_mgr.home_village_of, hmn_mgr.home_village_of+hmn_mgr.sum_num_humans);
//     int at_max = *std::max_element(hmn_mgr.at_village_of, hmn_mgr.at_village_of+hmn_mgr.sum_num_humans);
//     std::cout   << "2. max{Human.home_village} (" << home_max 
//                 <<")  < Village.Count(" << vll_mgr.sum_num_villages
//                 <<") AND max{Human.at_village} (" << at_max
//                 <<") < Village.Count(" << vll_mgr.sum_num_villages << ") ";
//     if (home_max < vll_mgr.sum_num_villages && at_max < vll_mgr.sum_num_villages) {
//         std::cout << "\033[32;47m[Pass]\033[0m";
//     } else {
//         std::cout << "\033[34;47m[Fail]\033[0m";
//     }
//     std::cout << "\n";

//     return EXIT_SUCCESS;
