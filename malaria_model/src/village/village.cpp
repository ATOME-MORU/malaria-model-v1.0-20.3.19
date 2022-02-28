#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include <vector>

#include <iomanip>

#include <algorithm> // std::fill_n, std::random_shuffle
#include <limits> //std::numeric_limits

#include <random>

#include <math.h> // pow, sqrt, cos

#include <cassert>

#include "util/util.h"
#include "util/randomness.h"
#include "village/village.h"
#include "human/human.h"
#include "human/blood_system.h"

namespace village {

VillageManager::VillageManager(

        const bool if_batch_mode,
        const int simulation_step_max,

        const std::string input_file_name,
        const char input_file_delimiter,
        const int rows_to_read,
        const int rows_to_skip,
        const bool init_overwrite_ra_rate_if,
        const float init_overwrite_ra_rate_with,
        
        const std::string func_treatment_general_version,

        const float G_timecomp,
        const float G_fullcourse,
        const float G_covab,
        const float G_nomp,
        const float asymptomatic_to_clinical_ratio,
        const float establish_malaria_post_during_mda,
        const float establish_malaria_post_during_survey,
        const int establish_malaria_post_on_day,

        const float treatment_rate_mda,

        const float treatment_coverage_mmc_wp2,
        const int treatment_drug_mmc_wp2_id,


        const bool init_with_uniform_prevalence_if,
        const float init_with_uniform_prevalence_with,

        const float init_transmission_coefficient_scaling_factor_a,
        const float init_transmission_coefficient_scaling_factor_b,
        const float init_infected_mosquito_to_infected_human_ratio,
        const float init_infectious_mosquito_to_infected_human_ratio,
        const float transmission_coefficient_to_beta_scaling_factor,

        const bool mosquito_seasonality_switched_on,

        const bool mosquito_seasonality_init_if_from_file,
        const std::string mosquito_seasonality_init_file_name,

        const float seasonality_amplitude,
        const float seasonality_amplitude_multiplier,

        const float cos_pi_multiplier,
        const float cos_days_offset,

        const int func_infection_m2h_version,

        const float infection_susceptibility_multiplier,
        const float infection_infectiousness_multiplier,

        const float mosquito_biting_rate,

        const float mosquito_incubation_days,
        const float mosquito_infected_days,
        const float mosquito_infectious_days,

        const float mosquito_min_filter_threshold,
        const int mosquito_min_filter_start_on_day,

        // Mobility parameters
        const bool mobility_enabled, 
        const float static_population_move_out_probability,
        const float static_population_return_home_probability,
        const float non_static_population_move_out_probability,
        const float non_static_population_return_home_probability

    ) : 
        kIf_batch_mode(if_batch_mode),

        kInit_overwrite_ra_rate_if(init_overwrite_ra_rate_if),
        kInit_overwirte_ra_rate_with(init_overwrite_ra_rate_with),

        kFunc_treatment_general_version(func_treatment_general_version),

        G_timecomp(G_timecomp),
        G_fullcourse(G_fullcourse),
        G_covab(G_covab),
        G_nomp(G_nomp),
        asymptomatic_to_clinical_ratio(asymptomatic_to_clinical_ratio),
            treatment_rate_clinical_has_mp((1.0 / G_timecomp) * G_fullcourse * G_covab),
            treatment_rate_clinical_no_mp(treatment_rate_clinical_has_mp * G_nomp),
        kEstablish_malaria_post_during_mda(establish_malaria_post_during_mda),
        kEstablish_malaria_post_during_survey(establish_malaria_post_during_survey),
        kEstablish_malaria_post_on_day(establish_malaria_post_on_day),
        
        kTreatment_rate_mda(treatment_rate_mda),

        kTreatment_coverage_mmc_wp2(treatment_coverage_mmc_wp2),
        kTreatment_drug_mmc_wp2(static_cast<common::DrugName>(treatment_drug_mmc_wp2_id)),

        kInit_with_uniform_prevalence_if(init_with_uniform_prevalence_if),
        kInit_with_uniform_prevalence_with(init_with_uniform_prevalence_with),

        kInit_transmission_coefficient_scaling_factor(
            init_transmission_coefficient_scaling_factor_a /
            init_transmission_coefficient_scaling_factor_b
        ),
        kInit_infected_mosquito_to_infected_human_ratio(init_infected_mosquito_to_infected_human_ratio),
        kInit_infectious_mosquito_to_infected_human_ratio(init_infectious_mosquito_to_infected_human_ratio),
        kTransmission_coefficient_to_beta_scaling_factor(transmission_coefficient_to_beta_scaling_factor),

        kMosquito_seasonality_switched_on(mosquito_seasonality_switched_on),

        kMosquito_seasonality_init_if_from_file(mosquito_seasonality_init_if_from_file),
        kMosquito_seasonality_init_file_name(mosquito_seasonality_init_file_name),

        kSeasonality_amplitude(seasonality_amplitude),
        kSeasonality_amplitude_multiplier(seasonality_amplitude_multiplier),

        kSeasonality_cos_pi_multiplier(cos_pi_multiplier),
        kSeasonality_cos_days_offset(cos_days_offset),

        kFunc_infection_m2h_version(func_infection_m2h_version),

        kInfection_susceptibility_multiplier(infection_susceptibility_multiplier),
        kInfection_infectiousness_multiplier(infection_infectiousness_multiplier),

        kMosquito_biting_rate(mosquito_biting_rate),

        kMosquito_incubation_probability(1.0/mosquito_incubation_days),
        kMosquito_death_probability_infected(1.0/mosquito_infected_days),
        kMosquito_death_probability_infectious(1.0/mosquito_infectious_days),

        kMosquito_min_filter_threshold(mosquito_min_filter_threshold),
        kMosquito_min_filter_start_on_day(mosquito_min_filter_start_on_day),


        kMobility_enabled(mobility_enabled),
        kStatic_population_move_out_probability(static_population_move_out_probability),
        kStatic_population_return_home_probability(static_population_return_home_probability),
        kNon_static_population_move_out_probability(non_static_population_move_out_probability),
        kNon_static_population_return_home_probability(non_static_population_return_home_probability),

        daily_record_of_new_bites_annual(rows_to_read, 0),
        daily_record_of_new_infections_annual(rows_to_read, 0),
        daily_record_of_new_clinical_cases_annual(rows_to_read, 0),
        daily_record_of_new_asymptomatic_cases_annual(rows_to_read, 0)

    {

    // TODO: rows_to_read == 0 >> read all rows

    int num_columns = columns::kInit_last - columns::kInit_first + 1;

    std::ifstream input_files_stream(input_file_name);
    
    if (!input_files_stream) {

        std::cout << "Could not open input file (" << input_file_name << ")." << std::endl;

    } else {

        sum_num_villages = 0;
        sum_total_population = 0;

        // float* x = new float[num_columns*rows_to_read];
        id_of = new int[rows_to_read];
        population_of   = new int[rows_to_read];
        longitude_of    = new float[rows_to_read];
        latitude_of     = new float[rows_to_read];
        transmission_coefficient_of = new float[rows_to_read];
        num_malaria_post_of             = new int[rows_to_read];
        num_days_malaria_post_opened_of = new int[rows_to_read];
        migrant_population_percentage_of = new float[rows_to_read];
        // artemisinin_resistance_percentage_of = new float[rows_to_read];

        infection_rates_of = new float *[rows_to_read];

        mosquito_infected_composition_of = new int* [rows_to_read];
        mosquito_infectious_composition_of = new int* [rows_to_read];

        sum_prevalence_count_of = new int *[rows_to_read];

        for (int vv = 0; vv < rows_to_read; vv++){
            int num_parasite_types = static_cast<int>(common::ParasiteType::Last)+1;

            infection_rates_of[vv] = new float[num_parasite_types];
            std::fill_n(infection_rates_of[vv], num_parasite_types, 0.0);

            mosquito_infected_composition_of[vv] = new int[num_parasite_types];
            std::fill_n(mosquito_infected_composition_of[vv], num_parasite_types, 0);

            mosquito_infectious_composition_of[vv] = new int[num_parasite_types];
            std::fill_n(mosquito_infectious_composition_of[vv], num_parasite_types, 0);

            sum_prevalence_count_of[vv] = new int[num_parasite_types];
            std::fill_n(sum_prevalence_count_of[vv], num_parasite_types, 0);

        }


        malaria_post_in_operation_of = new bool[rows_to_read];
        std::fill_n(malaria_post_in_operation_of, rows_to_read, false);


        for (const auto& pt: util::Enum_iterator<common::ParasiteType>()) {
            parasite_type_lookup.push_back(pt);
        }

        
        mosquito_seasonality_on_length = simulation_step_max+1; // TODO: check relation bt. +1 and season offset
        mosquito_seasonality_on = new float[mosquito_seasonality_on_length];
        std::fill_n(mosquito_seasonality_on, mosquito_seasonality_on_length, 0.0);

        treatment_rate_mp_clinical_of = new float[rows_to_read];
        treatment_rate_mp_asymptomatic_of = new float[rows_to_read];

        distances = new float *[rows_to_read];
        distances_shortest_path = new float*[rows_to_read];
        for (int vv = 0; vv < rows_to_read; vv++) {
            distances[vv] = new float[rows_to_read];
            distances_shortest_path[vv] = new float[rows_to_read];
            std::fill_n(distances[vv], rows_to_read, std::numeric_limits<float>::infinity());
            std::fill_n(distances_shortest_path[vv], rows_to_read, std::numeric_limits<float>::infinity());
        }

        // mosquito_bite_probability_of_on.resize(rows_to_read);
        // for (int vv = 0; vv < rows_to_read; vv++) {
        //     mosquito_bite_probability_of_on[vv].resize(simulation_step_max+1);
        // }
        beta_season_of.resize(rows_to_read);

        village_reg = new Village[rows_to_read];

        std::string str_line;


        for (int ii = 0; ii < rows_to_skip; ii++) {
            assert(!std::getline(input_files_stream, str_line).eof());
        }

        for (int i = 0; i < rows_to_read; i++) {
            
            village_reg[i].id = id_of + i;

            village_reg[i].population = population_of + i;
            village_reg[i].longitude = longitude_of + i;
            village_reg[i].latitude = latitude_of + i;
            village_reg[i].transmission_coefficient = transmission_coefficient_of + i;
            village_reg[i].num_malaria_post = num_malaria_post_of + i;
            village_reg[i].num_days_malaria_post_opened = num_days_malaria_post_opened_of + i;
            village_reg[i].migrant_population_percentage = migrant_population_percentage_of + i;
            // village_reg[i].artemisinin_resistance_percentage = artemisinin_resistance_percentage_of + i;

            // std::getline(input_files_stream, str_line); // TODO: rows_to_read > rows in file
            
            if (!std::getline(input_files_stream, str_line).eof()) {

                // std::cout << "line is:" << str_line << std::endl;

                std::string item;
                std::stringstream ss(str_line);
                
                int j = 0; 
                id_of[i] = sum_num_villages;
                j++; // column 0 contains village ID

                while (std::getline(ss, item, input_file_delimiter)) {

                    if (j<=num_columns) {

                        switch (j) {

                            // case columns::kID:
                            //     id_of[i] = sum_num_villages; //TODO: bring this outside switch
                            //     break;

                            // the order of these following cases does not matter
                            // but each case must direct to the correct array
                            case columns::kPopulation:
                                population_of[i] = std::stoi(item);
                                break;
                            
                            case columns::kLongitude:
                                longitude_of[i] = std::stof(item);
                                break;

                            case columns::kLatitude:
                                latitude_of[i] = std::stof(item);
                                break;

                            case columns::kTransmission_coefficient:
                                transmission_coefficient_of[i] = std::stof(item);
                                break;

                            case columns::kNum_malaria_post:
                                num_malaria_post_of[i] = std::stoi(item);
                                break;

                            case columns::kNum_days_malaria_post_opened:
                                num_days_malaria_post_opened_of[i] = std::stoi(item);
                                break;

                            case columns::kMigrant_population_percentage:
                                migrant_population_percentage_of[i] = std::stof(item);
                                break;

                            case columns::kArtemisinin_resistance_percentage:
                                // artemisinin_resistance_percentage_of[i] = std::stof(item);
                                if (this->kInit_overwrite_ra_rate_if) {
                                    this->infection_rates_of[i][static_cast<int>(common::ParasiteType::kRa)] = this->kInit_overwirte_ra_rate_with;
                                } else {
                                    assert(std::stof(item)<=1.0);
                                    this->infection_rates_of[i][static_cast<int>(common::ParasiteType::kRa)] = std::stof(item);
                                }

                                this->infection_rates_of[i][static_cast<int>(common::ParasiteType::kRb)] = 0.0;
                                this->infection_rates_of[i][static_cast<int>(common::ParasiteType::kR0)] = 1.0
                                     - this->infection_rates_of[i][static_cast<int>(common::ParasiteType::kRa)]
                                     - this->infection_rates_of[i][static_cast<int>(common::ParasiteType::kRb)];
                                 
                                break;

                        }

                    }

                    j++;

                }//while columns

                j--;

                if (j!=num_columns && !this->kIf_batch_mode) {

                    std::cout << "\033[31;47mVillageManager: Unexpected width of input data."
                              << " (Expecting "
                              << num_columns << " got "
                              << j << ")"
                              << "\033[0m" << std::endl;

                }

                sum_total_population += population_of[sum_num_villages];
                sum_num_villages++;

            } else {
                std::cout << "\033[31;47mVillageManager: Data file contains less number of rows than requested."
                          << "\033[0m" << std::endl;

                break;
            } // for .eof()

        } // for rows_to_read

        input_files_stream.close();

        // std::cout << "\npop[0]=" << population_of[0] << std::endl;
        // std::cin.ignore();

        this->at_village_register.resize(this->sum_num_villages);

        // this->init_mosquito_bite_probability_of_on();

        this->init_village_distances();
        this->static_mobility_network.resize(this->sum_num_villages);
        this->build_static_mobility_network();

        this->daily_record_of_new_bites.resize(this->sum_num_villages);
        this->daily_record_of_new_infections.resize(this->sum_num_villages);
        this->daily_record_of_wasted_bites.resize(this->sum_num_villages);

        this->sum_prevalence_infected_of.resize(this->sum_num_villages);

        this->sum_prevalence_of.resize(this->sum_num_villages);

    }

}

VillageManager::~VillageManager() {
    if (id_of) {delete[] id_of;}
    if (population_of) {delete[] population_of;}
    if (longitude_of) {delete[] longitude_of;}
    if (latitude_of) {delete[] latitude_of;}
    if (transmission_coefficient_of) {delete[] transmission_coefficient_of;}
    if (num_malaria_post_of) {delete[] num_malaria_post_of;}
    if (num_days_malaria_post_opened_of) {delete[] num_days_malaria_post_opened_of;}
    if (migrant_population_percentage_of) {delete[] migrant_population_percentage_of;}
    // if (artemisinin_resistance_percentage_of) {delete[] artemisinin_resistance_percentage_of;;}

    if (village_reg) {delete[] village_reg;}
    // if (distances) {delete[] distances;}

    for (int vv = 0; vv < this->sum_num_villages; vv++){
        delete[] infection_rates_of[vv];

        delete[] mosquito_infected_composition_of[vv];
        delete[] mosquito_infectious_composition_of[vv];

        delete[] sum_prevalence_count_of[vv];

        delete[] distances[vv];
        delete[] distances_shortest_path[vv];
    }

    delete[] infection_rates_of;

    delete[] mosquito_infected_composition_of;
    delete[] mosquito_infectious_composition_of;

    delete[] sum_prevalence_count_of;

    delete[] distances;
    delete[] distances_shortest_path;

    delete[] malaria_post_in_operation_of;
    delete[] mosquito_seasonality_on;

    delete[] treatment_rate_mp_clinical_of;
    delete[] treatment_rate_mp_asymptomatic_of;
    
}


void VillageManager::update_summary() {
    // Updates summary data
    
    // TODO: update sum_num_villages

    // sum_total_population = 0;
    
    // for(int ii = 0; ii < sum_num_villages; ii++) {
    //     sum_total_population += population_of[ii];
    // }

    static bool first_run(true);

    int yesterday_infected_sum = this->daily_record_of_mosquito_infected_sum;
    int yesterday_infectious_sum = this->daily_record_of_mosquito_infectious_sum;

    this->daily_record_of_mosquito_infected_sum = 0;
    this->daily_record_of_mosquito_infectious_sum = 0;

    for (int vv = 0; vv < this->sum_num_villages; vv++) {
        for (int pp = 1; pp < static_cast<int>(this->parasite_type_lookup.size()); pp++) {
            this->daily_record_of_mosquito_infected_sum += this->mosquito_infected_composition_of[vv][pp];
            this->daily_record_of_mosquito_infectious_sum += this->mosquito_infectious_composition_of[vv][pp];
        }
    }

    if (!first_run){
        assert( yesterday_infected_sum
                + this->daily_record_of_mosquito_newly_infected
                - this->daily_record_of_mosquito_death_infected
                - this->daily_record_of_mosquito_incubation
                ==
                this->daily_record_of_mosquito_infected_sum
        );
        assert( yesterday_infectious_sum
                - this->daily_record_of_mosquito_death_infectious
                + this->daily_record_of_mosquito_incubation
                ==
                this->daily_record_of_mosquito_infectious_sum
        );
    } else {
        first_run = false;
    }

}


void VillageManager::update_prevalence_count_of(const common::ParasiteType* dominant_type_array){
    int num_parasite_types = static_cast<int>(common::ParasiteType::Last)+1;
    for (int vv = 0; vv < this->sum_num_villages; vv++) {
        std::fill_n(this->sum_prevalence_count_of[vv], num_parasite_types, 0);
        for (const auto& hh : this->at_village_register.at(vv)) {
            this->sum_prevalence_count_of[vv][static_cast<int>(dominant_type_array[hh])]++;
        }
        // this->sum_prevalence_infected_of[vv] = 1 -
        //     this->sum_prevalence_count_of[vv][static_cast<int>(common::ParasiteType::kNoParasite)]
        //     / static_cast<float>(this->at_village_register[vv].size());

        // this->sum_prevalence_infected_of[vv] *= this->kOde_k;
    }
}

// Depends: a call to update_prevalence_count_of
void VillageManager::update_prevalence_of() {
    for(int vv = 0; vv < this->sum_num_villages; vv++) {
        this->sum_prevalence_of[vv] = 1 -
            this->sum_prevalence_count_of[vv][static_cast<int>(common::ParasiteType::kNoParasite)]
                / static_cast<float>(this->at_village_register[vv].size());
        // std::cout << ">" << this->sum_prevalence_of[vv] << "\n";
    }
}

// int** VillageManager::get_prevalence_count_of() const {
//     return this->sum_prevalence_count_of;
// }
std::vector<int> VillageManager::get_prevalence_greater_than(float prevalence_threshold) const {
    std::vector<int> village_ids;
    float prev = 0.0;
    for (int vv = 0; vv < this->sum_num_villages; vv++) {
        prev = 1-
            this->sum_prevalence_count_of[vv][static_cast<int>(common::ParasiteType::kNoParasite)]
            / static_cast<float>(this->at_village_register[vv].size());
        if (prev >= prevalence_threshold) {
            village_ids.push_back(vv);
        }
    }
    return village_ids;
}
// std::vector<int> VillageManager::survey_prevalence_greater_than(
//         float prevalence_threshold,
//         int pcr_sample_size,
//         float pcr_sensitivity
//     ) const {

//     // int inf=RO+RA+RB;
//     // propinf=float(inf)/float(POP);
//     // int detectcases = int(poissonHighPrecision(propinf*SAMPLES));
//     // int seennumber = int((detectcases*G_sensitivity));
//     // detecprev = float(seennumber)/float(SAMPLES);

//     // if (detecprev >= G_trigger){

//     std::vector<int> village_ids;

//     int vv_pop = 0;
//     float vv_true_prev = 0.0;
//     int vv_sample_size = 0;
//     int vv_expected_cases = 0;
//     int vv_detectable_cases = 0;
//     int vv_detected_cases = 0;
//     float vv_detected_prev = 0.0;


//     for (int vv = 0; vv < this->sum_num_villages; vv++) {
//         vv_pop = static_cast<int>(this->at_village_register[vv].size());
//         vv_true_prev = 1-
//             this->sum_prevalence_count_of[vv][static_cast<int>(common::ParasiteType::kNoParasite)]
//             / static_cast<float>(vv_pop);
//         vv_sample_size = (vv_pop >= pcr_sample_size) ? pcr_sample_size : vv_pop;
//         vv_expected_cases = vv_true_prev*vv_sample_size;
//         // std::cout << vv_true_prev << ", " << vv_sample_size << std::endl;
//         if (vv_expected_cases > 0) {

//             util::set_random_numbers_poisson(&vv_detectable_cases, 1, vv_expected_cases);
//             vv_detected_cases = vv_detectable_cases * pcr_sensitivity;
//             vv_detected_prev = static_cast<float>(vv_detected_cases) / static_cast<float>(vv_sample_size);

//             if (vv_detected_prev >= prevalence_threshold) {
//                 village_ids.push_back(vv);
//             }

//         }
//     }

//     return village_ids;

// }

std::vector<int> VillageManager::survey_prevalence_greater_than(
        float prevalence_threshold,
        int pcr_sample_size,
        float pcr_sensitivity
    ) {

    std::vector<int> village_ids;
    for (int vv = 0; vv < this->sum_num_villages; vv++) {
        if ( this->get_prevalence_of_village_by_survey(vv,pcr_sample_size,pcr_sensitivity) >= prevalence_threshold) {
            village_ids.push_back(vv);
        }
    }
    return village_ids;
}

float VillageManager::get_prevalence_of_village_by_survey(
        int vv,
        int pcr_sample_size,
        float pcr_sensitivity
    ) {

    if (this->kEstablish_malaria_post_during_survey) {
        this->open_malaria_post(vv);
    }

    int vv_pop = static_cast<int>(this->at_village_register[vv].size());
    int vv_sample_size = (vv_pop >= pcr_sample_size) ? pcr_sample_size : vv_pop;
    float vv_true_prev = this->get_prevalence_of(vv);
    float vv_expected_cases = vv_true_prev*static_cast<float>(vv_sample_size);
    // std::cout << vv_true_prev << ", " << vv_sample_size << std::endl;

    if (vv_expected_cases > 0) {

        // util::set_random_numbers_poisson(&vv_detectable_cases, 1, vv_expected_cases);
        int vv_detectable_cases = util::get_rand_poisson(vv_expected_cases);

        int vv_detected_cases = vv_detectable_cases * pcr_sensitivity;
        return static_cast<float>(vv_detected_cases) / static_cast<float>(vv_sample_size);

    } else {
        return 0.0;
    }

}


const std::vector<float>& VillageManager::get_sum_prevalence_infected_of(
        const common::ParasiteType* dominant_type_array,
        const float* infectiousness_array
    ) {

    for (int vv = 0; vv < this->sum_num_villages; vv++) {
        float per_village_infectiousness = 0.0;
        for (const auto& hh : this->at_village_register.at(vv)) {
            if (dominant_type_array[hh] != common::ParasiteType::kNoParasite) {
                per_village_infectiousness += infectiousness_array[hh];
            }
        }
        per_village_infectiousness /= static_cast<float>(this->at_village_register[vv].size());
        // this->sum_prevalence_infected_of[vv] = per_village_infectiousness * this->kOde_k;
        this->sum_prevalence_infected_of[vv] = per_village_infectiousness;
    }


    return this->sum_prevalence_infected_of;

}


void VillageManager::print_summary() const{
    // Prints summary data

    std::cout << "\033[34;47;4mVillage Data Summary\n";
    std::cout << "\033[24mNumber of Villages: " << sum_num_villages << "\n";
    std::cout << "  Total Population: " << sum_total_population << "\n";
    std::cout << "\033[0m\n";

}

void VillageManager::print_all() const{
    util::print_data_table(this, this->sum_num_villages, 6 , 1);
}

void VillageManager::print_one(const int index_village) const{
        std::cout   << std::setw(5) << id_of[index_village]
                    << "/" << std::setw(3) << population_of[index_village]
                    << "@[" << std::setprecision(4) << std::setw(4) << longitude_of[index_village]
                    << ","  << std::setprecision(4) << std::setw(4) << latitude_of[index_village]
                    << "] " << std::setw(8) << transmission_coefficient_of[index_village]
                    << ", " << std::setprecision(1) << std::setw(2) << num_malaria_post_of[index_village]
                    << ", " << std::setprecision(3) << std::setw(4) << num_days_malaria_post_opened_of[index_village]
                    << ", " << std::setprecision(3) << std::setw(4) << migrant_population_percentage_of[index_village] * 100 << "%"
                    << ", " << std::setprecision(3) << std::setw(4)
                        << this->infection_rates_of[index_village][static_cast<int>(common::ParasiteType::kRa)] * 100 << "%"
                    << ", " << std::setprecision(3) << std::setw(4)
                        << this->infection_rates_of[index_village][static_cast<int>(common::ParasiteType::kRb)] * 100 << "%"
                    << ", " << std::setprecision(3) << std::setw(4)
                        << this->infection_rates_of[index_village][static_cast<int>(common::ParasiteType::kR0)] * 100 << "%"
                    << " ";

}

void VillageManager::init_treatment_rates_malaria_post(){
    for (int vv = 0; vv < this->sum_num_villages; vv++) {
        update_treatment_rates_malaria_post(vv);
    }
}
void VillageManager::update_treatment_rates_malaria_post(int village_id){
    // if (this->num_malaria_post_of[village_id] > 0) {
    if (this->malaria_post_in_operation_of[village_id]) {
        this->treatment_rate_mp_clinical_of[village_id] = this->treatment_rate_clinical_has_mp;
    } else {
        this->treatment_rate_mp_clinical_of[village_id] = this->treatment_rate_clinical_no_mp;
    }
    this->treatment_rate_mp_asymptomatic_of[village_id] =
            this->treatment_rate_mp_clinical_of[village_id] * this->asymptomatic_to_clinical_ratio;
}

void VillageManager::establish_malaria_posts() {
    for (int vv = 0; vv < this->sum_num_villages; vv++) {
        this->open_malaria_post(vv);
    }
}
void VillageManager::open_malaria_post(int village_id){
    this->malaria_post_in_operation_of[village_id] = true;
    this->update_treatment_rates_malaria_post(village_id);
    std::cout << "Malaria post opened at village " << village_id << "\n";
}

// void VillageManager::open_malaria_post(int time_step){

//     // /** open posts according to data **/
//     // if (timestep==mptime[i] && mp[i]>0){
//     //     if (village->post==0){
//     //         village->post = 1;
//     //         village->mptime=mptime[i];
//     //     }
//     // }
//     for (int vv = 0; vv < this->sum_num_villages; vv++) {
//         // if (this->num_malaria_post_of[vv] > 0) continue;
//         if (this->num_days_malaria_post_opened_of[vv] == time_step){
//             // this->num_malaria_post_of[vv]++;
//             if(this->num_malaria_post_of[vv] > 0){
//                 // std::cout << "opening mp at village" << vv
//                 //             << ", before: treat_rate_clinical:" << treatment_rate_mp_clinical_of[vv]
//                 //             << ", treat_rate_asymptonmatic:" << treatment_rate_mp_asymptomatic_of[vv];
//                 this->malaria_post_in_operation_of[vv] = true;
//                 this->update_treatment_rates_malaria_post(vv);
//                 // std::cout   << ", after: treat_rate_clinical:" << treatment_rate_mp_clinical_of[vv]
//                 //             << ", treat_rate_asymptonmatic:" << treatment_rate_mp_asymptomatic_of[vv];
//                 // std::cin.ignore();
//             }
//         }
//     }
// }
void VillageManager::treatment_by_malaria_post(
        human::BloodSystemManager* bld_mgr,
        human::HumanManager* hmn_mgr,
        int time_step
    ){

    assert(bld_mgr->sum_num_systems == hmn_mgr->sum_num_humans);

    float* random_numbers = new float[hmn_mgr->sum_num_humans];
    util::set_random_numbers_uniform(random_numbers, hmn_mgr->sum_num_humans);

    if (time_step == this->kEstablish_malaria_post_on_day) {
        this->establish_malaria_posts();
    }

    this->treatment_by_malaria_post_kernel(
            bld_mgr->get_drug_of(),
            hmn_mgr->get_given_drug(),
            bld_mgr->get_is_clinical_of(),
            bld_mgr->get_dominant_type_of(),
            hmn_mgr->sum_num_humans,
            this->at_village_register,
            this->sum_num_villages,
            random_numbers,
            this->treatment_rate_mp_clinical_of,
            this->treatment_rate_mp_asymptomatic_of
    );

    delete[] random_numbers;
}
void VillageManager::treatment_by_malaria_post_kernel(
        const common::DrugName* existing_drug,
        common::DrugName* given_drug,
        const bool* is_clinical,
        const common::ParasiteType* dominant_type,
        const int num_humans,
        // const std::vector<std::vector<int>>& at_vll_reg,
        // const boost::container::vector<boost::container::vector<int>>& at_vll_reg,
        const Register_t& at_vll_reg,
        const int at_vll_reg_size,
        const float* rnd_nmbrs,
        const float* rates_clinical,
        const float* rates_asymptomatic
    ) {
    
    // std::cout << "rates_clinical[0]:" << rates_clinical[0] << ", rates_asymptomatic[0]:" << rates_asymptomatic[0] << "\n";

    const common::DrugName kMalaria_post_drug = common::DrugName::kPiparte;

    int mp_treatment_counter = 0;

    float applicable_rate = 0.0;

    for (int vv = 0; vv < at_vll_reg_size; vv++) {
        for (const auto& hhi : at_vll_reg.at(vv)) {
            assert(hhi < num_humans);
            if (existing_drug[hhi] != common::DrugName::kPiperaquine
                && existing_drug[hhi] != common::DrugName::kPiparte
                && existing_drug[hhi] != common::DrugName::kArtesunate ) {
                
                    // std::cout << "hhi:" << hhi << " d:" << dominant_type[hhi] << std::endl;

                if (dominant_type[hhi] != common::ParasiteType::kNoParasite) {
                    if (is_clinical[hhi]){
                        applicable_rate = rates_clinical[vv];
                    } else {
                        applicable_rate = rates_asymptomatic[vv];
                    }
                    // std::cout << "rate:" << applicable_rate << std::endl;

                    if (rnd_nmbrs[hhi] < applicable_rate) {
                        given_drug[hhi] = kMalaria_post_drug;
                        mp_treatment_counter++;
                    }
                }
                
            }
        }
    }

    std::cout << "\tmp_treatment_counter=" << mp_treatment_counter << std::endl;

}

void VillageManager::treatment_by_mda(
        const int village_id,
        const common::DrugName mda_drug
    ){

    // const float kMda_treatment_rate = 1.0/1.5; // G_tauab = 1.0/1.5;

    // kMda_treatment_rate is defined here if it is village-bound
    // otherwise move to mda.cpp

    int village_pop = this->at_village_register[village_id].size();

    float* random_numbers = new float[village_pop];
    util::set_random_numbers_uniform(random_numbers, village_pop);

    if (this->kEstablish_malaria_post_during_mda) {    
        // this->malaria_post_in_operation_of[village_id] = true;
        // this->update_treatment_rates_malaria_post(village_id);
        this->open_malaria_post(village_id);
    }

    this->treatment_by_mda_kernel(
            this->v_hmn_mgr->get_given_drug(),
            this->v_hmn_mgr->sum_num_humans,
            this->at_village_register[village_id],
            random_numbers,
            // kMda_treatment_rate,
            this->kTreatment_rate_mda,
            mda_drug
        );

    delete[] random_numbers;
}
void VillageManager::treatment_by_mda_kernel(
        // const common::DrugName* existing_drug,
        common::DrugName* given_drug,
        const int human_index_max,
        // const std::vector<int>& list_of_humans,
        // const boost::container::vector<int>& list_of_humans,
        const Register_row_t& list_of_humans,
        const float* rnd_nmbrs,
        const float treatment_rate,
        const common::DrugName mda_drug
    ){
    int villager_index = 0;
    for (const auto& hhi : list_of_humans) {
        assert(hhi < human_index_max);
        // if (existing_drug[hhi] != common::DrugName::kPiperaquine
        //     && existing_drug[hhi] != common::DrugName::kPiparte
        //     && existing_drug[hhi] != common::DrugName::kArtesunate ) {
            
                // std::cout << "hhi:" << hhi << " d:" << dominant_type[hhi] << std::endl;

            // if (rnd_nmbrs[hhi] < treatment_rate) {
            if (rnd_nmbrs[villager_index] < treatment_rate) {
                given_drug[hhi] = mda_drug;
            }
            
        // }
            villager_index++;
    }

}

void VillageManager::step_treatment_general(
        human::BloodSystemManager* bld_mgr,
        human::HumanManager* hmn_mgr,
        int time_step
    ) {

    for (auto& ii: this->daily_record_of_drug_failure_by_new_treatment) {
        ii = 0;
    }
    for (auto& ii: this->daily_record_of_new_treatments_by_parasite_type) {
        ii = 0;
    }

    if (this->kFunc_treatment_general_version == "mp") {
        
        this->treatment_by_malaria_post(
            bld_mgr,
            hmn_mgr,
            time_step
        );

    } else if ( this->kFunc_treatment_general_version == "mmc_wp2") {

        this->func_general_treatment_mmc_wp2(
            hmn_mgr->get_given_drug(),
            bld_mgr->get_drug_of(),
            bld_mgr->get_drug_given_of(),
            bld_mgr->get_drug_given_day_of(),
            bld_mgr->get_is_clinical_of(),
            bld_mgr->get_mellow_rate(),
            bld_mgr->get_drug_failure_if_clinical_on_day(),
            bld_mgr->get_dominant_type_of(),
            this->at_village_register,
            time_step,
            this->daily_record_of_drug_failure_by_new_treatment,
            this->daily_record_of_new_treatments_by_parasite_type
        );

    } else if (this->kFunc_treatment_general_version == "mmc_wp2_rev1") {

        this->func_general_treatment_mmc_wp2_rev1(
            hmn_mgr->get_given_drug(),
            bld_mgr->get_drug_of(),
            bld_mgr->get_drug_given_of(),
            bld_mgr->get_drug_given_day_of(),
            bld_mgr->get_is_clinical_of(),
            bld_mgr->get_mellow_rate(),
            bld_mgr->get_drug_failure_if_clinical_on_day(),
            bld_mgr->get_dominant_type_of(),
            this->at_village_register,
            time_step,
            this->daily_record_of_drug_failure_by_new_treatment,
            this->daily_record_of_new_treatments_by_parasite_type
        );

    } else {

        std::cout << "No valid selection. Skipping general treatment.\n";
    }

}
void VillageManager::func_general_treatment_mmc_wp2(
        common::DrugName* given_drug,
        const common::DrugName* drug_of,
        const common::DrugName* drug_given_of,
        const int* drug_given_day_of,
        const bool* is_clinical_of,
        const float mellow_rate,
        const int drug_failure_if_clinical_on_day,
        const common::ParasiteType* dominant_type_of,
        const Register_t& at_vll_reg,
        const int time_step,
        std::vector<int>& daily_record_of_drug_failure_by_new_treatment,
        std::vector<int>& daily_record_of_new_treatments_by_parasite_type
    ) {


    // const float kTreatment_coverage_mmc_wp2 = 0.5;
    // (void) mellow_rate;
    const float kTreatment_rate_if_clinical = 1.0 - pow(1.0-this->kTreatment_coverage_mmc_wp2, mellow_rate);
    // const float kTreatment_rate_if_clinical = 1.0 - pow(kTreatment_coverage_mmc_wp2, 1.0/9.0);

    const common::DrugName kDrug1 = common::DrugName::kPiparte;
    const common::DrugName kDrug2 = common::DrugName::kSP;

    const float kDrug1_probability = 2.0 / 5.0;

    // std::cout << "mellow_rate=" << mellow_rate << "\n";
    std::cout << "kTreatment_rate_if_clinical=" << kTreatment_rate_if_clinical << "\n";
    std::cout << "kDrug1_probability=" << kDrug1_probability << "\n";

    float rnd_nmbr = 0.0;

    for (uint vv = 0; vv < at_vll_reg.size(); vv++) {
        for (const auto& hh : at_vll_reg.at(vv)) {
            if (is_clinical_of[hh]) {

                if (drug_of[hh] == common::DrugName::kNoDrug){
                    
                    util::set_random_numbers_uniform(&rnd_nmbr,1);
                    if (rnd_nmbr < kTreatment_rate_if_clinical) {

                        if (drug_given_day_of[hh] >= 0) {
                            if ((time_step - drug_given_day_of[hh]) < drug_failure_if_clinical_on_day) {

                                // std::cout << "drug_of[hh]=" << drug_of[hh]
                                //         << ", drug_given_of[hh]=" << drug_given_of[hh]
                                //         << ", drug_given_day_of[hh]=" << drug_given_day_of[hh] << "\n";

                                daily_record_of_drug_failure_by_new_treatment[static_cast<int>(drug_given_of[hh])]++;
                            }
                        }

                        util::set_random_numbers_uniform(&rnd_nmbr,1);
                        if (rnd_nmbr < kDrug1_probability) {
                            given_drug[hh] = kDrug1;
                        } else {
                            given_drug[hh] = kDrug2;
                        }

                        daily_record_of_new_treatments_by_parasite_type[static_cast<int>(dominant_type_of[hh])]++;
                    }
                
                }

            // } else {
            //     util::set_random_numbers_uniform(&rnd_nmbr,1);
            //     if (rnd_nmbr < 0.1 * kTreatment_rate_if_clinical) {
            //         given_drug[hh] = kDrug1;
            //     }
            }
        }
    }

}

void VillageManager::func_general_treatment_mmc_wp2_rev1(
        common::DrugName* given_drug,
        const common::DrugName* drug_of,
        const common::DrugName* drug_given_of,
        const int* drug_given_day_of,
        const bool* is_clinical_of,
        const float mellow_rate,
        const int drug_failure_if_clinical_on_day,
        const common::ParasiteType* dominant_type_of,
        const Register_t& at_vll_reg,
        const int time_step,
        std::vector<int>& daily_record_of_drug_failure_by_new_treatment,
        std::vector<int>& daily_record_of_new_treatments_by_parasite_type
    ) {


    // const float kTreatment_coverage_mmc_wp2 = 0.4;
    // (void) mellow_rate;
    const float kTreatment_rate_if_clinical = 1.0 - pow(1.0-this->kTreatment_coverage_mmc_wp2, mellow_rate);
    // const float kTreatment_rate_if_clinical = 1.0 - pow(kTreatment_coverage_mmc_wp2, 1.0/9.0);
    // const float kTreatment_rate_if_clinical = 0.4;


    // const common::DrugName kDrug1 = common::DrugName::kPiparte;
    // const common::DrugName kDrug2 = common::DrugName::kSP;

    // const float kDrug1_probability = 2.0 / 5.0;

    // std::cout << "mellow_rate=" << mellow_rate << "\n";
    // std::cout << "kTreatment_rate_if_clinical=" << kTreatment_rate_if_clinical << "\n";
    // std::cout << "kDrug1_probability=" << kDrug1_probability << "\n";

    float rnd_nmbr = 0.0;

    for (uint vv = 0; vv < at_vll_reg.size(); vv++) {
        for (const auto& hh : at_vll_reg.at(vv)) {
            if (is_clinical_of[hh]) {

                if (drug_of[hh] == common::DrugName::kNoDrug){
                    
                    util::set_random_numbers_uniform(&rnd_nmbr,1);
                    if (rnd_nmbr < kTreatment_rate_if_clinical) {

                        if (drug_given_day_of[hh] >= 0) {
                            if ((time_step - drug_given_day_of[hh]) < drug_failure_if_clinical_on_day) {

                                // std::cout << "drug_of[hh]=" << drug_of[hh]
                                //         << ", drug_given_of[hh]=" << drug_given_of[hh]
                                //         << ", drug_given_day_of[hh]=" << drug_given_day_of[hh] << "\n";

                                daily_record_of_drug_failure_by_new_treatment[static_cast<int>(drug_given_of[hh])]++;
                            }
                        }

                        // util::set_random_numbers_uniform(&rnd_nmbr,1);
                        // if (rnd_nmbr < kDrug1_probability) {
                            given_drug[hh] = this->kTreatment_drug_mmc_wp2;

                            // std::cout << "given drug: " << given_drug[hh] << "\n";
                            // std::cin.ignore();
                        // } else {
                        //     given_drug[hh] = kDrug2;
                        // }

                        daily_record_of_new_treatments_by_parasite_type[static_cast<int>(dominant_type_of[hh])]++;
                    }
                
                }

            // } else {
            //     util::set_random_numbers_uniform(&rnd_nmbr,1);
            //     if (rnd_nmbr < 0.1 * kTreatment_rate_if_clinical) {
            //         given_drug[hh] = kDrug1;
            //     }
            }
        }
    }

}

// void VillageManager::give_drug_to_village(human::HumanManager* human_manager,
//         const int index_village, const common::DrugName drug){
//     human_manager->update_given_drug_by_at_village(index_village, drug);
// }
// void VillageManager::give_drug_to_village(
//         human::HumanManager* human_manager,
//         const int index_village,
//         const common::DrugName drug
//         ){
//     // for(auto const& hhi : at_village_register[index_village]){
//     for(auto const& hhi : at_village_register.at(index_village)) {
//         // *hhi->given_drug = drug;
//         if (human_manager->home_village_of[hhi] == index_village) {
//             human_manager->given_drug_of[hhi] = drug;
//         }
//     }
//     // human_manager->update_given_drug_by_at_village(index_village, drug);
// }

void VillageManager::init_get_villagers(human::HumanManager* human_manager){
    int human_index = 0;
    for (int vv = 0; vv < this->sum_num_villages; vv++) {

        this->at_village_register[vv].reserve(this->population_of[vv]+this->kRegister_reserve_margin);
        
        for (int hh = 0; hh < this->population_of[vv]; hh++){
            // this->at_village_register[vv].push_back(human_manager->human_reg + human_index);
            this->at_village_register[vv].push_back(human_index);

            human_manager->set_human_at_village(human_index, vv);
            human_manager->set_human_home_village(human_index, vv);
            human_index++;
        }
    }
    this->v_hmn_mgr = human_manager;

    this->step_update_at_village_attractiveness();
}

// void VillageManager::init_infections(float* village_bite_rate, float* village_ra_rate) {
//     I=int((vregion[i]->v[3]*5.1/0.75)*vregion[i]->v[0])
//     I = ( I > number ? number : I );
//     i_init[i]=I;

//     float initresa = 1.0;
//     G_pa = initresa;   // initial artemisinin resistance level
//     float presa=float(initresa*(vregion[i]->v[7]));
//     pa[i]=presa;

//     ia = int (pa[i] * i_init[i]);
//     ib = int (pb[i] * i_init[i]);
//     io = i_init[i]-(ia+ib);

// }
void VillageManager::init_infections(
        common::ParasiteType* infections_array,
        int infections_array_size
    ) {

    this->init_infections_kernel(
            this->kInit_with_uniform_prevalence_if,
            this->kInit_with_uniform_prevalence_with,

            this->transmission_coefficient_of,
            this->infection_rates_of,
            this->at_village_register,
            this->sum_num_villages,

            infections_array,
            infections_array_size,
            
            this->mosquito_infected_composition_of,
            this->mosquito_infectious_composition_of,

            this->kInit_transmission_coefficient_scaling_factor,
            this->kInit_infected_mosquito_to_infected_human_ratio,
            this->kInit_infectious_mosquito_to_infected_human_ratio
            // this->init_infection_human_to_mosquito_infected_probability,
            // this->init_infection_human_to_mosquito_infectious_probability
        );
}
void VillageManager::init_infections_kernel(
        const bool init_with_uniform_prevalence_if,
        const float init_with_uniform_prevalence_with,

        const float* transmission_coefficient_array,
        const float* const* parasite_proportion_array,
        // std::vector<std::vector<int>>& villager_reg,
        // boost::container::vector<boost::container::vector<int>>& villager_reg,
        Register_t& villager_reg,
        const int vll_array_size, //number of villages

        common::ParasiteType* infections_array,
        const int infections_array_size, // number of humans

        int** infected_msq_composition_array,
        int** infectious_msq_composition_array,

        const float init_transmission_coefficient_scaling_factor,
        const float init_infected_mosquito_to_infected_human_ratio,
        const float init_infectious_mosquito_to_infected_human_ratio
        // const float infected_probability,
        // const float infectious_probability
    ) {

    assert(static_cast<int>(villager_reg.size())==vll_array_size);

    static std::random_device rd;
    static std::default_random_engine rnd_eng(rd());
    // static std::default_random_engine rnd_eng; // use this for reproduciable results

    // const float init_transmission_coefficient_scaling_factor = 4.1/0.75;

    for (int vv = 0; vv < vll_array_size; vv++) {

        int village_population = villager_reg[vv].size();
        int village_infected_population = 0;

        if (init_with_uniform_prevalence_if) {
            village_infected_population = village_population * init_with_uniform_prevalence_with;
        } else {
            village_infected_population = village_population * ( 
                                            transmission_coefficient_array[vv]
                                            * init_transmission_coefficient_scaling_factor
                                        );        
        }

        if (village_infected_population > village_population) {
            if (!this->kIf_batch_mode) {
                std::cout << "init_infections_kernel: Village "
                            << vv << " infected population ("
                            << village_infected_population
                            << ") greater than village population ("
                            << village_population
                            << "). Press ENTER to continue, or CTRL+C to terminate.\n";

                std::cin.ignore();
            }
            village_infected_population = village_population;
        }

        std::vector<int> population_composition( static_cast<int>(common::ParasiteType::Last)+1, 0 );
        // population_composition[static_cast<int>(common::ParasiteType::Ra)]

        std::shuffle(villager_reg[vv].begin(), villager_reg[vv].end(), rnd_eng);

        // std::vector<int>::const_iterator ihh_of_village = villager_reg[vv].begin();
        // boost::container::vector<int>::iterator ihh_of_village = villager_reg[vv].begin();
        Register_row_t::iterator ihh_of_village = villager_reg[vv].begin();

        for (auto pt: util::Enum_iterator<common::ParasiteType>()) {
            int pt_index = static_cast<int>(pt);
            if (pt_index == 0) { // common::ParasiteType::kNoParasite
                
                population_composition[pt_index] = 0; // to be corrected 

            } else {  

                population_composition[pt_index] =
                    village_infected_population * parasite_proportion_array[vv][pt_index];

                for (int hh_of_village = 0; hh_of_village < population_composition[pt_index]; hh_of_village++){

                    assert(*ihh_of_village < infections_array_size);

                    infections_array[*ihh_of_village] = pt;
                    // std::cout << "init_infections: infections_array[*ihh_of_village="<< *ihh_of_village <<"] = " << pt << "\n";
                    ihh_of_village++;

                    // infect mosquitos at the same time
                    float rnd_nmbr = 0.0;
                    util::set_random_numbers_uniform(&rnd_nmbr,1);

                    if (rnd_nmbr < init_infected_mosquito_to_infected_human_ratio){
                        infected_msq_composition_array[vv][pt_index]++;
                    }

                    if (rnd_nmbr < init_infectious_mosquito_to_infected_human_ratio){
                        infectious_msq_composition_array[vv][pt_index]++;
                    }

                    // if (rnd_nmbr < infectious_probability){
                    //     infectious_msq_composition_array[vv][pt_index]++;
                    // }

                    // if (rnd_nmbr < infected_probability){
                    //     infected_msq_composition_array[vv][pt_index]++;
                    // }
                }

            }
            // std::cout << "init_infections: population_composition[pt_index=" << pt_index << "] = "
            //             << population_composition[pt_index] << "\n";
            // std::cout << "init_infections: infected_msq_composition_array[vv=" << vv << "][pt_index=" << pt_index << "]"
            //             << infected_msq_composition_array[vv][pt_index] << "\n";
            // std::cout << "init_infections: infectious_msq_composition_array[vv=" << vv << "][pt_index=" << pt_index << "]"
            //             << infectious_msq_composition_array[vv][pt_index] << "\n";
        }

    }

}

// void VillageManager::init_mosquito_bite_probability_of_on(){
//     int village_count = static_cast<int>(this->mosquito_bite_probability_of_on.size());
//     int step_count = static_cast<int>(this->mosquito_bite_probability_of_on[0].size());

//     assert(village_count == this->sum_num_villages);

//     if (this->kMosquito_seasonality_init_if_from_file) {
//         this->init_mosquito_seasonality_from_file(this->kMosquito_seasonality_init_file_name);
//     }

//     for (int vv = 0; vv < village_count; vv++) {
//         for (int tt = 0; tt < step_count; tt++) {

//             float beta = this->transmission_coefficient_of[vv]
//                             * this->kTransmission_coefficient_to_beta_scaling_factor;

//             if (!this->kMosquito_seasonality_switched_on) {
//                 this->mosquito_bite_probability_of_on[vv][tt] = beta;
//                 continue;
//             }

//             if (this->kMosquito_seasonality_init_if_from_file) {
//                 this->mosquito_bite_probability_of_on[vv][tt] =
//                     this->func_beta_to_beta_season_given_multiplier(
//                         beta,
//                         this->mosquito_seasonality_on[tt + this->kSeasonality_step_offset],
//                         this->kSeasonality_amplitude,
//                         this->kSeasonality_amplitude_multiplier
//                     );
//             } else {
//                 this->mosquito_bite_probability_of_on[vv][tt] =
//                     this->func_beta_to_beta_season_cos(
//                         beta,
//                         this->kSeasonality_amplitude,
//                         this->kSeasonality_cos_pi_multiplier,
//                         this->kSeasonality_cos_days_offset,
//                         tt // TODO: offset here?
//                     );
//             }

//         }
//     }

// }

void VillageManager::step_update_beta_season_of(int at_time_step){
    for (int vv = 0; vv < this->sum_num_villages; vv++) {

            float beta = this->transmission_coefficient_of[vv]
                            * this->kTransmission_coefficient_to_beta_scaling_factor;

            if (!this->kMosquito_seasonality_switched_on) {
                this->beta_season_of[vv] = beta;
                continue;
            }

            if (this->kMosquito_seasonality_init_if_from_file) {
                this->beta_season_of[vv] =
                    this->func_beta_to_beta_season_given_multiplier(
                        beta,
                        this->mosquito_seasonality_on[at_time_step + this->kSeasonality_step_offset],
                        this->kSeasonality_amplitude,
                        this->kSeasonality_amplitude_multiplier
                    );
            } else {
                this->beta_season_of[vv] =
                    this->func_beta_to_beta_season_cos(
                        beta,
                        this->kSeasonality_amplitude,
                        this->kSeasonality_cos_pi_multiplier,
                        this->kSeasonality_cos_days_offset,
                        at_time_step
                    );
            }

    }
    this->beta_season_of_last_updated_on = at_time_step;
}

float VillageManager::func_beta_to_beta_season_cos(
        float beta,
        float seasonality_amplitude,
        float cos_pi_multiplier,
        float cos_days_offset,
        int time_step
    ) {
//betaseason[i] = betas[i]+(G_amp*betas[i]*cos(2.0*3.1416*((double(timestep)-90.0)/365.0)));
    // return beta + (
    //             seasonality_amplitude * beta * cos(
    //                 4.0 * 3.1416 * ( ( static_cast<float>(time_step) - 34.0 ) / 365.0))
    //         );
    // return beta + (
    //             seasonality_amplitude * beta * cos(
    //                 cos_pi_multiplier * 3.1416
    //                 * ( ( static_cast<float>(time_step) - cos_days_offset ) / 365.0)
    //             )
    //         );
    return beta + (
                seasonality_amplitude * beta * cos(
                    cos_pi_multiplier * common::kPI
                    * ( ( static_cast<float>(time_step) - cos_days_offset ) / common::kDays_per_year)
                )
            );
}

float VillageManager::func_beta_to_beta_season_given_multiplier(
        float beta,
        float seasonality_multiplier,    // from seasonality file
        float seasonality_amplitude,     // G_amp
        float seasonality_amplitude_multiplier // 0.5
    ) {
// betaseason[i] = (betas[i]*seasons)+(betas[i]*G_amp*0.5);
    return beta * seasonality_multiplier
             + beta * seasonality_amplitude * seasonality_amplitude_multiplier;
}

int VillageManager::func_beta_season_to_num_bites(
        float beta_season,
        int population
    ) {
    int num_bites = 0;
    float mean = beta_season * population;
    if (mean <= 0.0) {
        return 0;
    }
    util::set_random_numbers_poisson(&num_bites, 1, beta_season * population);

    return num_bites;
}

int VillageManager::func_num_infectious_mosq_to_num_bites(
        double num_infectious_mosq,
        float ode_k
    ) {
    int num_bites = 0;
    float mean = num_infectious_mosq * ode_k;

    if (!std::isnormal(mean)) {
        return 0;
    }

    if (mean <= 0.0) {
        return 0;
    }

    // std::cout << "num_infectious_mosq: " << num_infectious_mosq;
    // std::cout << ", ode_k: " << ode_k << "\n";
    util::set_random_numbers_poisson(&num_bites, 1, mean);
    return num_bites;
}

int VillageManager::func_num_infectious_mosq_to_num_bites_uniform(
        double num_infectious_mosq,
        float ode_k
    ) {
    const double kEpsilon = 0.001;
    if (num_infectious_mosq < kEpsilon) {
        return 0;
    } else {
        return func_num_infectious_mosq_to_num_bites_uniform(
                    static_cast<int>(ceil(num_infectious_mosq)),
                    ode_k
                );
    }
}
int VillageManager::func_num_infectious_mosq_to_num_bites_uniform(
        int num_infectious_mosq,
        float ode_k
    ) {
    int num_bites = 0;
    float rnd_nmbr = 0;
    for (int mm = 0; mm < num_infectious_mosq; mm++) {
        util::set_random_numbers_uniform(&rnd_nmbr, 1);
        if (rnd_nmbr < ode_k) {
            num_bites++;
        }
    }
    return num_bites;
}

// int VillageManager::calculate_number_of_new_bites_in_village(
//         float prv_dist_base_value,
//         float prv_dist_base_to_trns_cff_scaling_factor,
//         float seasonality_factor,
//         float seasonality_amplitude,
//         float seasonality_amplitude_factor,
//         uint population,
//         int time_step
//     ){

//     float beta = prv_dist_base_value
//                 * prv_dist_base_to_trns_cff_scaling_factor;

// // Method 1:
// //betaseason[i] = betas[i]+(G_amp*betas[i]*cos(2.0*3.1416*((double(timestep)-90.0)/365.0)));
//     float beta_season = beta + (
//                             seasonality_amplitude * beta * cos(
//                                 2.0 * 3.1416 * ( ( double(time_step) - 90.0 ) / 365.0))
//                         );
//     (void) seasonality_factor;
//     (void) seasonality_amplitude_factor;

// // Method 2:
// //betaseason[i] = (betas[i]*seasons)+(betas[i]*G_amp*0.5);

//     // float beta_season = beta * seasonality_factor
//     //                     + beta * seasonality_amplitude * seasonality_amplitude_factor;
//     // (void) time_step;

//     // float beta_season = beta * seasonality_factor;
//     // (void) seasonality_amplitude;
//     // (void) seasonality_amplitude_factor;
//                         // + beta * seasonality_amplitude * seasonality_amplitude_factor;
    
//     // std::cout << std::setprecision(8);
//     // std::cout << "prv_dist_base_value:" << prv_dist_base_value << "\n";
//     // std::cout << "prv_dist_base_to_trns_cff_scaling_factor:" << prv_dist_base_to_trns_cff_scaling_factor << "\n";
//     // std::cout << "beta:" << beta << "\n";
//     // std::cout << "beta_sseason:" << beta_season << "\n";

//     int num_new_bites_mean = beta_season * population;

//     if (num_new_bites_mean == 0) return 0;

//     int num_new_bites = 0;
//     util::set_random_numbers_poisson(&num_new_bites, 1, num_new_bites_mean);
//     return num_new_bites;
// }

void VillageManager::step_update_at_village_attractiveness() {

    if (this->v_hmn_mgr->get_if_attractiveness_enabled()) {

        if (this->kMobility_enabled || this->at_village_attractiveness.empty() ) {
            
            std::vector<std::vector<float>> new_at_village_attractiveness;
            new_at_village_attractiveness.resize(this->sum_num_villages);

            for (int vv = 0; vv < this->sum_num_villages; vv++) {
                for(const int & hh: this->at_village_register[vv]) {
                    new_at_village_attractiveness[vv].push_back(this->v_hmn_mgr->get_attractiveness(hh));
                }
            }
            this->at_village_attractiveness = new_at_village_attractiveness;
        
        }

    }

    this->at_village_attractiveness_updated = true;
    
}

std::vector<int> VillageManager::func_get_humans_to_bite(int num_new_bites, int village_id) {

    std::vector<int> humans_to_bite(num_new_bites, -1);

    if (this->v_hmn_mgr->get_if_attractiveness_enabled()) {
        assert(this->at_village_attractiveness_updated);

        int* new_bites_target_human_aux;
        new_bites_target_human_aux = new int[num_new_bites];
        util::set_random_numbers_discrete(new_bites_target_human_aux, num_new_bites, this->at_village_attractiveness[village_id]);
        for (int bb = 0; bb < num_new_bites; bb++) {
            humans_to_bite[bb] =  this->at_village_register[village_id][new_bites_target_human_aux[bb]];
        }
        delete[] new_bites_target_human_aux;

    } else {        

        // float* new_bites_target_human_aux;
        // new_bites_target_human_aux = new float[num_new_bites];
        // util::set_random_numbers_uniform(new_bites_target_human_aux, num_new_bites);
        // for (int bb = 0; bb < num_new_bites; bb++) {
        //     humans_to_bite[bb] = villager_reg[village_id][new_bites_target_human_aux[bb] * villager_reg[village_id].size()];
        // }
        // delete[] new_bites_target_human_aux;

        for (int& hh : humans_to_bite) {
            hh = this->get_random_human_id_at_village(village_id);
        }
    }

    return humans_to_bite;

}

void VillageManager::infection_m2h( // mosquito to human
        common::ParasiteType* infections_array,
        const int infections_array_size,
        const float* susceptibility_array,
        const int susceptibility_array_size,
        const common::StageName* stage_array,
        const int time_step,
        const std::vector<double>& num_infectious_mosq_from_ode,
        // const double ode_feed_rate,
        const bool if_using_ode
    ) {

    switch(this->kFunc_infection_m2h_version) {
        case 0 :
            assert(this->beta_season_of_last_updated_on == time_step);

            this->infection_m2h_kernel(
                // this->transmission_coefficient_of,
                // this->mosquito_seasonality_on[time_step + this->kSeasonality_step_offset],

                this->sum_num_villages,
                this->at_village_register,
                this->mosquito_infectious_composition_of,
                this->parasite_type_lookup,

                // this->kTransmission_coefficient_to_beta_scaling_factor,
                // this->kSeasonality_amplitude,
                // this->kSeasonality_amplitude_multiplier,
                this->kInfection_susceptibility_multiplier,
                // this->kParasite_mosquito_to_human_probability,

                susceptibility_array,
                susceptibility_array_size,
                stage_array,

                infections_array,
                infections_array_size,

                // time_step,
                // this->mosquito_bite_probability_of_on,
                this->beta_season_of,

                this->daily_record_of_new_bites_sum,
                this->daily_record_of_new_infections_sum,

                this->daily_record_of_new_bites,
                this->daily_record_of_new_infections,

                num_infectious_mosq_from_ode,
                this->kMosquito_biting_rate
                // ode_feed_rate
            );
            break;

        case 1 :
            this->infection_m2h_kernel_use_all_bites(

                this->sum_num_villages,
                this->at_village_register,
                this->mosquito_infectious_composition_of,
                this->parasite_type_lookup,

                this->kInfection_susceptibility_multiplier,

                susceptibility_array,
                susceptibility_array_size,
                stage_array,

                infections_array,
                infections_array_size,

                this->daily_record_of_new_bites_sum,
                this->daily_record_of_new_infections_sum,

                this->daily_record_of_new_bites,
                this->daily_record_of_new_infections
            );
            break;

        case 2 :

            this->infection_m2h_kernel_per_bite(

                time_step,
                this->kMosquito_min_filter_threshold,
                this->kMosquito_min_filter_start_on_day,

                this->sum_num_villages,
                // this->at_village_register,
                this->mosquito_infectious_composition_of,

                this->kInfection_susceptibility_multiplier,

                susceptibility_array,
                susceptibility_array_size,
                stage_array,

                infections_array,
                infections_array_size,

                this->daily_record_of_new_bites_sum,
                this->daily_record_of_new_infections_sum,
                this->daily_record_of_wasted_bites_sum,

                this->daily_record_of_new_bites,
                this->daily_record_of_new_infections,
                this->daily_record_of_wasted_bites,

                num_infectious_mosq_from_ode,
                this->kMosquito_biting_rate,
                // ode_feed_rate,
                if_using_ode
            );
            break;


    }

  //           if ((time_step % 10) == 0) {
  //           if ((time_step ) == 1460) {

  // std::ofstream myfile;
  // myfile.open ("waste.txt");
  //               for (int vv = 0; vv < this->sum_num_villages; vv++) {
  //                   myfile << vv << "," << this->daily_record_of_wasted_bites[vv]
  //                          << "," << this->sum_prevalence_infected_of[vv]
  //                          << "," << this->sum_prevalence_count_of[vv][static_cast<int>(common::ParasiteType::kNoParasite)]
  //                          << "," << this->at_village_register[vv].size()
  //                          << "," << num_infectious_mosq_from_ode.at(vv)
  //                          << "," << this->mosquito_infectious_composition_of[vv][static_cast<int>(common::ParasiteType::kR0)]
  //                          << "\n";
  //               }
  // myfile.close();
  // std::cin.ignore();
  //           }
            
}

void VillageManager::infection_m2h_kernel(
        // const float* transmission_coefficient_array,
        // const float seasonality_factor,
        
        const int vll_array_size,
        // const std::vector<std::vector<int>> villager_reg,
        // const boost::container::vector<boost::container::vector<int>> villager_reg,
        const Register_t& villager_reg,
        const int * const* infectious_msq_composition_array,
        const std::vector<common::ParasiteType> parasite_types,

        // const float transmission_coefficient_to_beta_scaling_factor,
        // const float seasonality_amplitude,
        // const float seasonality_amplitude_factor,
        const float bite_success_probability,

        const float* susceptibility_array,
        const int susceptibility_array_size,
        const common::StageName* stage_array,

        common::ParasiteType* infections_array,//[h]
        const int infections_array_size,

        // const int time_step,
        // const std::vector< std::vector<float> >& mosquito_bite_probability,
        const std::vector<float>& beta_season_of_array,

        int& daily_record_of_new_bites_sum,
        int& daily_record_of_new_infections_sum,

        std::vector<int>& daily_record_of_new_bites,
        std::vector<int>& daily_record_of_new_infections,

        const std::vector<double>& num_infectious_mosq_from_ode,
        const double ode_feed_rate
    ) {

    bool using_ode = (num_infectious_mosq_from_ode.at(0) > 0.0);

    assert(susceptibility_array_size == infections_array_size);

    static std::vector<int> village_infection_weights;


    daily_record_of_new_bites_sum = 0;
    daily_record_of_new_infections_sum = 0;

    // int daily_record_of_already_infectious_sum = 0;

    if (using_ode) {
        DEBUG_MSG( "infection_m2h_kernel: " << " Using ODE");
    } else {
        DEBUG_MSG( "infection_m2h_kernel: " << " Not using ODE");
    }

    // int num_new_bites_total = 0;
    // int num_inc_total = 0;

    for(int vv = 0; vv < vll_array_size; vv++){

        // 1. Get the number of new bites

        int num_new_bites = 0;

        // num_new_bites = VillageManager::calculate_number_of_new_bites_in_village(
        //     transmission_coefficient_array[vv],
        //     transmission_coefficient_to_beta_scaling_factor,
        //     seasonality_factor,
        //     seasonality_amplitude,
        //     seasonality_amplitude_factor,
        //     villager_reg[vv].size(),
        //     time_step
        // );
        // (void) transmission_coefficient_array;
        // (void) seasonality_factor;
        // (void) transmission_coefficient_to_beta_scaling_factor;
        // (void) seasonality_amplitude;
        // (void) seasonality_amplitude_factor;

        daily_record_of_new_bites[vv] = 0;
        daily_record_of_new_infections[vv] = 0;

        int per_village_mosq_pop = 0;
        for (uint pp = 0; pp < parasite_types.size(); pp++){
            per_village_mosq_pop += infectious_msq_composition_array[vv][pp];
        }

        if (using_ode) {
            // std::cout << vv << ":" << num_infectious_mosq_from_ode.at(vv) << "\n";
            // num_new_bites = VillageManager::func_num_infectious_mosq_to_num_bites(
            //     num_infectious_mosq_from_ode.at(vv),
            //     ode_k
            // );

            num_new_bites = VillageManager::func_num_infectious_mosq_to_num_bites_uniform(
                num_infectious_mosq_from_ode.at(vv),
                ode_feed_rate
            );

        } else {
            num_new_bites = VillageManager::func_beta_season_to_num_bites(
                // mosquito_bite_probability[vv][time_step],
                beta_season_of_array[vv],
                // per_village_mosq_pop * 400
                villager_reg[vv].size()
            );
            // std::cout << "villager_reg[vv].size()" << villager_reg[vv].size() << ", "
            //           << "per_village_mosq_pop " << per_village_mosq_pop << "\n";
        }



        // if (num_new_bites != 0) {
        //     std::cout << "v" << vv <<" new bites:" << num_new_bites << "\n";
        // }

        // std::cout << "num_new_bites:" << num_new_bites << "\n";
        if (num_new_bites == 0) continue;

        // 2. Prepare for 5.3.1
        // int new_bites_parasite_type_index[num_new_bites] = {};
        // std::vector<int> village_infection_weights;
        if (!using_ode) {
            
            village_infection_weights.clear();
            int total_weight = 0;
            for (uint pp = 0; pp < parasite_types.size(); pp++){
                village_infection_weights.push_back(infectious_msq_composition_array[vv][pp]);
                total_weight += infectious_msq_composition_array[vv][pp];
            }
            if (total_weight == 0) continue;
            // util::set_random_numbers_discrete<int>(new_bites_parasite_type_index, num_new_bites, village_infection_weights);
        }
        // num_new_bites_total += num_new_bites;

        daily_record_of_new_bites[vv] = num_new_bites;
        daily_record_of_new_bites_sum += num_new_bites;

        // 3. Get list of humans to be bitten

        std::vector<int> humans_to_bite = this->func_get_humans_to_bite(num_new_bites, vv);
        assert(num_new_bites == static_cast<int>(humans_to_bite.size()));


        // // 4. Get the random numbers for biting probability
        // float* new_bites_probabilities_aux;
        // new_bites_probabilities_aux = new float[num_new_bites];
        // util::set_random_numbers_uniform(new_bites_probabilities_aux, num_new_bites);

        // 5. Implement the bites
        int num_succ_bites = 0;
        // for (int bb = 0; bb < num_new_bites; bb++) {
        for (const int& human_index: humans_to_bite) {
            this->v_hmn_mgr->inc_log_bitten_times_m2h(human_index);
            // num_inc_total++;

            // // 5.1 Who to bite
            // int human_village_index = new_bites_target_human_aux[bb] * villager_reg[vv].size();
            // int human_index = villager_reg[vv][human_village_index];

            // std::cout << "human_village_index=" << human_village_index << ", human_index=" << human_index << "\n";

            assert(human_index < susceptibility_array_size);
            assert(human_index >=0);

            // 5.1.2 If host already has infectious stage parasite, then anti-bodies prevent new infection
            if (stage_array[human_index] == common::StageName::kInfectious) {
                // daily_record_of_already_infectious_sum++;
                continue;
            }

            // 5.2 Get infection threshold Mosquito to Human
            float bite_success_threshold = susceptibility_array[human_index];
            if (stage_array[human_index] == common::StageName::kNotInSystem) {
                bite_success_threshold *= bite_success_probability;
            }

            // 5.3 If successful bite

            // std::cout << "new_bites_probabilities_aux[bb]:" << new_bites_probabilities_aux[bb]
            //             << ", bite_success_threshold:" << bite_success_threshold << "\n";

            // if ( new_bites_probabilities_aux[bb] < bite_success_threshold ) {
            if ( util::get_rand_uniform() < bite_success_threshold ) {

                if (using_ode) {

                    infections_array[human_index] = common::ParasiteType::kR0;

                } else {

                    // 5.3.1 Draw parasite type at the last minute
                    int bite_parasite_index = 0;
                    util::set_random_numbers_discrete<int>(&bite_parasite_index, 1, village_infection_weights);
                    
                    // 5.3.2 Record bite
                    infections_array[human_index] = parasite_types.at(bite_parasite_index);

                    
                }
            
                num_succ_bites++;
            }

            //

        }

        // std::cout << "num_new_bites=" << num_new_bites << ", num_succ_bites=" << num_succ_bites << "\n";

        daily_record_of_new_infections[vv] = num_succ_bites;
        daily_record_of_new_infections_sum += num_succ_bites;

        // delete[] new_bites_target_human_aux;
        // delete[] new_bites_probabilities_aux;

    }

    // std::cout << "VM: num_new_bites_total=" << num_new_bites_total << ", num_inc_total=" << num_inc_total << "\n";
    // std::cout << "daily_record_of_new_bites_sum=" << daily_record_of_new_bites_sum
    //             << ", daily_record_of_new_infections_sum=" << daily_record_of_new_infections_sum
    //             << ", daily_record_of_already_infectious_sum" << daily_record_of_already_infectious_sum
    //             << "\n";

}

// drugrestest version of infection()
// Differences:
// 1. All infectious bites get transmitted
// 2. All parasite types get transmitted proportionally
void VillageManager::infection_m2h_kernel_use_all_bites(
        const int vll_array_size,
        const Register_t& villager_reg,
        const int * const* infectious_msq_composition_array,
        const std::vector<common::ParasiteType> parasite_types,

        const float bite_success_probability,

        const float* susceptibility_array,
        const int susceptibility_array_size,
        const common::StageName* stage_array,

        common::ParasiteType* infections_array,//[h]
        const int infections_array_size,

        int& daily_record_of_new_bites_sum,
        int& daily_record_of_new_infections_sum,

        std::vector<int>& daily_record_of_new_bites,
        std::vector<int>& daily_record_of_new_infections
    ) {



    assert(static_cast<int>(villager_reg.size())==vll_array_size);

    assert(susceptibility_array_size == infections_array_size);

    daily_record_of_new_bites_sum = 0;
    daily_record_of_new_infections_sum = 0;

    for(int vv = 0; vv < vll_array_size; vv++){

        daily_record_of_new_bites[vv] = 0;
        daily_record_of_new_infections[vv] = 0;

        assert(infectious_msq_composition_array[vv][0] == 0);

        for (uint pp = 0; pp < parasite_types.size(); pp++){
            // if (parasite_types.at(pp) != common::ParasiteType::kNoParasite) {

                for (int bb = 0; bb < infectious_msq_composition_array[vv][pp]; bb++) {

                    daily_record_of_new_bites[vv]++;
                    daily_record_of_new_bites_sum++;

                    float rnd_nmbr = 0.0;
                    util::set_random_numbers_uniform(&rnd_nmbr, 1);

                    int human_village_index = rnd_nmbr * villager_reg[vv].size();
                    int human_index = villager_reg[vv][human_village_index];

                    assert(human_index < susceptibility_array_size);
                    
                    if (stage_array[human_index] == common::StageName::kInfectious) {continue;}

                    float bite_success_threshold = susceptibility_array[human_index];
                    if (stage_array[human_index] == common::StageName::kNotInSystem) {
                        bite_success_threshold *= bite_success_probability;
                    }

                    util::set_random_numbers_uniform(&rnd_nmbr, 1);
                    if ( rnd_nmbr < bite_success_threshold ) {

                        infections_array[human_index] = parasite_types.at(pp);
                    
                        daily_record_of_new_infections[vv]++;
                        daily_record_of_new_infections_sum++;

                    }
                }
            // }
        }
    }

    // int new_infections_count = 0;
    // float susceptibility_sum = 0.0;
    // for (int hh = 0; hh < susceptibility_array_size; hh++) {
    //     if (infections_array[hh] != common::ParasiteType::kNoParasite) {
    //         new_infections_count++;
    //     }
    //     susceptibility_sum += susceptibility_array[hh];
    // }
    // std::cout << "new_infections_count=" << daily_record_of_new_infections_sum
    //             << "\nnew_infections_h_count=" << new_infections_count
    //             << "\nsusceptibility_sum="<< susceptibility_sum << "\n";
    // std::cin.ignore();

}

void VillageManager::infection_m2h_kernel_per_bite(

        const int time_step,
        const float mosquito_min_filter_threshold,
        const int mosquito_min_filter_start_on_day,
        
        const int vll_array_size,
        // const Register_t& villager_reg,
        const int * const* infectious_msq_composition_array,
        // const std::vector<common::ParasiteType> parasite_types,

        const float bite_success_probability, // susceptibility multiplier

        const float* susceptibility_array,
        const int susceptibility_array_size,
        const common::StageName* stage_array,

        common::ParasiteType* infections_array,//[h]
        const int infections_array_size,

        // const int time_step,
        // const std::vector< std::vector<float> >& mosquito_bite_probability,
        // const std::vector<float>& beta_season_of_array,

        int& daily_record_of_new_bites_sum,
        int& daily_record_of_new_infections_sum,
        int& daily_record_of_wasted_bites_sum,

        std::vector<int>& daily_record_of_new_bites,
        std::vector<int>& daily_record_of_new_infections,
        std::vector<int>& daily_record_of_wasted_bites,

        const std::vector<double>& num_infectious_mosq_from_ode,
        const double ode_feed_rate,
        const bool if_using_ode
    ) {

    // bool using_ode = ((num_infectious_mosq_from_ode.at(0)+num_infectious_mosq_from_ode.at(1)) > 0.0);
    bool using_ode = if_using_ode;

    if(using_ode){
        std::cout << "using ODE\n";

    } else {
        std::cout << "NOT using ODE\n";
        // std::cin.ignore();
    }

    assert(susceptibility_array_size == infections_array_size);


    daily_record_of_new_bites_sum = 0;
    daily_record_of_new_infections_sum = 0;
    daily_record_of_wasted_bites_sum = 0;

    // int daily_record_of_already_infectious_sum = 0;

    for(int vv = 0; vv < vll_array_size; vv++){

        daily_record_of_new_bites[vv] = 0;
        daily_record_of_new_infections[vv] = 0;
        daily_record_of_wasted_bites[vv] = 0;

        int num_succ_bites = 0;

        // int temp = 0;
        for (const auto& pt: util::Enum_iterator<common::ParasiteType>()) {
            if (pt != common::ParasiteType::kNoParasite) {

                if (using_ode && pt != common::ParasiteType::kR0) {

                    continue;
                }

                // int num_infectious_mosq = 0;

        // 1. # infectious mosquito to # infectious bites

                int num_new_bites = 0;

                if (using_ode) {
                    // num_infectious_mosq = num_infectious_mosq_from_ode.at(vv);
                    // num_new_bites = VillageManager::func_num_infectious_mosq_to_num_bites_uniform(
                    //     num_infectious_mosq_from_ode.at(vv),
                    //     ode_k
                    // );

                    float num_infectious_mosq = num_infectious_mosq_from_ode.at(vv);
                    if (time_step >= mosquito_min_filter_start_on_day
                         && num_infectious_mosq < mosquito_min_filter_threshold) {
                        num_infectious_mosq = 0.0;
                    }

                    num_new_bites = VillageManager::func_num_infectious_mosq_to_num_bites(
                        // num_infectious_mosq_from_ode.at(vv),
                        num_infectious_mosq,
                        ode_feed_rate
                    );
                    // std::cout << "ode in use\n";
                } else {
                    // num_infectious_mosq = infectious_msq_composition_array[vv][static_cast<int>(pt)];
                    // num_new_bites = VillageManager::func_num_infectious_mosq_to_num_bites_uniform(
                    //     infectious_msq_composition_array[vv][static_cast<int>(pt)],
                    //     ode_k
                    // );
                    num_new_bites = VillageManager::func_num_infectious_mosq_to_num_bites(
                        infectious_msq_composition_array[vv][static_cast<int>(pt)],
                        ode_feed_rate
                    );
                }

                // (void) num_infectious_mosq;

                // int num_new_bites = num_infectious_mosq * ode_k;

                // std::cout << "pt=" << pt <<"\n";
                // // std::cout << "num_infectious_mosq=" << num_infectious_mosq <<"\n";
                // std::cout << "num_new_bites=" << num_new_bites <<"\n";


                // if (num_infectious_mosq == 0) continue;

                std::vector<int> humans_to_bite = this->func_get_humans_to_bite(num_new_bites, vv);

                
                // for (int mm = 0; mm < num_infectious_mosq; mm++) {
                // for (int mm = 0; mm < num_new_bites; mm++) {
                for (const int& human_index : humans_to_bite) {
                    this->v_hmn_mgr->inc_log_bitten_times_m2h(human_index);


                    // float rnd_nmbr = 0;
                    // util::set_random_numbers_uniform(&rnd_nmbr, 1);
                    // if (rnd_nmbr < ode_k) {

                        daily_record_of_new_bites[vv]++;
                        daily_record_of_new_bites_sum++;

        // 2. The bitten human

                        // float rnd_nmbr2 = 0;
                        // util::set_random_numbers_uniform(&rnd_nmbr2, 1);

                        // int reg_index = rnd_nmbr2 * villager_reg[vv].size();
                        // int human_index = villager_reg[vv][reg_index];

                        // std::cout << "vv=" << vv << ";reg_index=" << reg_index
                        //           << ";human_index=" << human_index << "\n";

                        if (stage_array[human_index] == common::StageName::kInfectious) {
                            daily_record_of_wasted_bites[vv]++;
                            daily_record_of_wasted_bites_sum++;
                            continue;
                        }

                        float bite_success_threshold = susceptibility_array[human_index];
                        if (stage_array[human_index] == common::StageName::kNotInSystem) {
                            bite_success_threshold *= bite_success_probability;
                        }

                        // util::set_random_numbers_uniform(&rnd_nmbr2, 1);
                        // if (rnd_nmbr2 < bite_success_threshold) {
                        if (util::get_rand_uniform() < bite_success_threshold) {

        // 3. Successful bite

                            infections_array[human_index] = pt;
                            num_succ_bites++;
                        } else {

                            daily_record_of_wasted_bites[vv]++;
                            daily_record_of_wasted_bites_sum++;
                            
                        }

                    // }

                } // each num_new_bites
            }

        } // pt

        daily_record_of_new_infections[vv] = num_succ_bites;
        daily_record_of_new_infections_sum += num_succ_bites;

    } // vv

    assert(daily_record_of_new_bites_sum == (daily_record_of_new_infections_sum+daily_record_of_wasted_bites_sum));

    // std::cout << "daily_record_of_new_bites_sum=" << daily_record_of_new_bites_sum
    //             << ", daily_record_of_new_infections_sum=" << daily_record_of_new_infections_sum
    //             << ", daily_record_of_wasted_bites_sum" << daily_record_of_wasted_bites_sum
    //             << "\n";

    // std::cin.ignore();
}


void VillageManager::infection_h2m(
        const float* infectiousness_array,
        const common::ParasiteType* dominant_type_array,
        const int host_array_size,
        const int time_step
    ) {

    assert(this->beta_season_of_last_updated_on == time_step);

    this->infection_h2m_kernel(
        // this->transmission_coefficient_of,
        // this->mosquito_seasonality_on[time_step + this->kSeasonality_step_offset],

        this->sum_num_villages,
        this->at_village_register,
        this->mosquito_infected_composition_of,

        // this->kTransmission_coefficient_to_beta_scaling_factor,
        // this->kSeasonality_amplitude,
        // this->kSeasonality_amplitude_multiplier,
        this->kInfection_infectiousness_multiplier,
        // this->kParasite_human_to_mosquito_probability,

        infectiousness_array,
        dominant_type_array,
        host_array_size,

        // time_step,
        // this->mosquito_bite_probability_of_on,
        this->beta_season_of,

        this->daily_record_of_mosquito_newly_infected,
        this->daily_record_of_mosquito_new_rb
    );

}
void VillageManager::infection_h2m_kernel(
        // const float* transmission_coefficient_array,
        // const float seasonality_factor,
        
        const int vll_array_size,
        // const std::vector<std::vector<int>> villager_reg,
        // const boost::container::vector<boost::container::vector<int>> villager_reg,
        const Register_t& villager_reg,
        int** infected_msq_composition_array,
        // const std::vector<common::ParasiteType> parasite_types,

        // const float transmission_coefficient_to_beta_scaling_factor,
        // const float seasonality_amplitude,
        // const float seasonality_amplitude_factor,
        const float bite_success_probability_h2m,

        // host arrays
        const float* infectiousness_array,
        const common::ParasiteType* dominant_type_array,
        const int host_array_size,

        // const int time_step,
        // const std::vector< std::vector<float> > mosquito_bite_probability,
        const std::vector<float>& beta_season_of_array,

        // common::ParasiteType* infections_array,//[h]
        // const int infections_array_size
        int& daily_record_of_mosquito_newly_infected,
        int& daily_record_of_mosquito_new_rb
    ) {

    daily_record_of_mosquito_newly_infected = 0;
    daily_record_of_mosquito_new_rb = 0;

    // int num_new_bites_sum = 0;
    // float newbites_mean_sum = 0;
    // float beta_season_sum = 0.0;
    // int pop_sum = 0;

    for(int vv = 0; vv < vll_array_size; vv++){

        // 1. Get the number of new bites

        int num_new_bites = 0;

        // num_new_bites = VillageManager::calculate_number_of_new_bites_in_village(
        //     transmission_coefficient_array[vv],
        //     transmission_coefficient_to_beta_scaling_factor,
        //     seasonality_factor,
        //     seasonality_amplitude,
        //     seasonality_amplitude_factor,
        //     villager_reg[vv].size(),
        //     time_step
        // );
        // (void) transmission_coefficient_array;
        // (void) seasonality_factor;
        // (void) transmission_coefficient_to_beta_scaling_factor;
        // (void) seasonality_amplitude;
        // (void) seasonality_amplitude_factor;
        num_new_bites = VillageManager::func_beta_season_to_num_bites(
            // mosquito_bite_probability[vv][time_step],
            beta_season_of_array[vv],
            villager_reg[vv].size()
        );

        // pop_sum += villager_reg[vv].size();
        // num_new_bites_sum += num_new_bites;
        // newbites_mean_sum += mosquito_bite_probability[vv][time_step] * villager_reg[vv].size();
        // beta_season_sum += mosquito_bite_probability[vv][time_step];
        // std::cout << "num_new_bites h2m:" << num_new_bites << "\n";

        if (num_new_bites == 0) continue;

        // // 2. Get the aux variable for humans bitten by these new bites
        // float* new_bites_target_human_aux;
        // new_bites_target_human_aux = new float[num_new_bites];
        // util::set_random_numbers_uniform(new_bites_target_human_aux, num_new_bites);

        std::vector<int> humans_to_bite = this->func_get_humans_to_bite(num_new_bites, vv);

        // 3. Implement bites
        // for (int bb = 0; bb < num_new_bites; bb++) {
        for (const int& human_index: humans_to_bite) {
            this->v_hmn_mgr->inc_log_bitten_times_h2m(human_index);

            // // 4.1 Who to bite
            // int human_village_index = new_bites_target_human_aux[bb] * villager_reg[vv].size();
            // int human_index = villager_reg[vv][human_village_index];

            assert(human_index < host_array_size);

            if (dominant_type_array[human_index] == common::ParasiteType::kNoParasite) continue;

            // 4.2 Get infection threshold Human to Mosquito
            float bite_success_threshold_h2m = infectiousness_array[human_index] * bite_success_probability_h2m;

            if (bite_success_threshold_h2m == 0) continue;

            // 4.3 If bite successful
            float new_bite_probability_aux = 0;
            util::set_random_numbers_uniform(&new_bite_probability_aux, 1);

            // std::cout << "infectiousness_array[human_index]=" << infectiousness_array[human_index]
            //         << "bite_success_threshold_h2m=" << bite_success_threshold_h2m 
            //         << "new_bite_probability_aux=" << new_bite_probability_aux << "\n";
            // std::cin.ignore();
            
            // if (new_bite_probability_aux < bite_success_threshold_h2m) {
            if (util::get_rand_uniform() < bite_success_threshold_h2m) {

                // 4.3.1 Add infection to mosquito bites composition
                infected_msq_composition_array[vv][static_cast<int>(dominant_type_array[human_index])]++;
                daily_record_of_mosquito_newly_infected++;

                if (dominant_type_array[human_index] == common::ParasiteType::kRb) {
                    daily_record_of_mosquito_new_rb++;
                }

            }

        }

        // delete[] new_bites_target_human_aux;

    }
    
    // std::cout << "h2m_new_bites_sum=" << num_new_bites_sum 
    //         << "\nh2m_newbites_mean_sum=" << newbites_mean_sum
    //         << "\nh2m_beta_season_sum=" << beta_season_sum
    //         << "\nh2m_pop_sum=" << pop_sum
    //         << "\nh2m_newly_infected=" << daily_record_of_mosquito_newly_infected << "\n";
    // std::cin.ignore();

}


void VillageManager::init_mosquito_seasonality_from_file(std::string file_name){
    std::ifstream input_file_stream(file_name);
    assert(input_file_stream);

    std::string str_line;

    for (int ii = 0; ii < this->mosquito_seasonality_on_length; ii++) {
        if(!std::getline(input_file_stream, str_line).eof()) {
            this->mosquito_seasonality_on[ii] = std::stof(str_line);
        } else {
            std::cout << "Simulation length greater than available seasonality definition "
                        << "from file " << file_name << "\n";
        }
    }

    // for (int ii = 0; ii < this->mosquito_seasonality_on_length; ii++) {
    //     std::cout << "Initialised seasonality vector: "
    //                 << this->mosquito_seasonality_on[ii] << "\n";
    // }
}

void VillageManager::mosquito_survival_step(){

    this->mosquito_survival_step_kernel(
        this->mosquito_infected_composition_of,
        this->mosquito_infectious_composition_of,
        this->sum_num_villages,
        static_cast<int>(this->parasite_type_lookup.size()),
        this->kMosquito_death_probability_infected,
        this->kMosquito_death_probability_infectious,

        this->daily_record_of_mosquito_death_infected,
        this->daily_record_of_mosquito_death_infectious
    );

}
void VillageManager::mosquito_survival_step_kernel(
        int** infected_composition_array,
        int** infectious_composition_array,
        const int array_size_vv,
        const int array_size_pp,
        const float death_probability_infected,
        const float death_probability_infectious,

        int& daily_record_of_mosquito_death_infected,
        int& daily_record_of_mosquito_death_infectious
    ){

    float rnd_nmbr = 0.0;
    int temp_num_bites = 0;

    daily_record_of_mosquito_death_infected = 0;
    daily_record_of_mosquito_death_infectious = 0;

    for (int vv = 0; vv < array_size_vv; vv++) {
        for (int pp = 1; pp < array_size_pp; pp++) { // note that this start at 1, ignoring kNoParasite

            temp_num_bites = infected_composition_array[vv][pp];
            for (int ii = 0; ii < temp_num_bites; ii++ ) {
                util::set_random_numbers_uniform(&rnd_nmbr, 1);
                if (rnd_nmbr < death_probability_infected) {
                    infected_composition_array[vv][pp]--;

                    daily_record_of_mosquito_death_infected++;
                }
            }

            temp_num_bites = infectious_composition_array[vv][pp];
            for (int ii = 0; ii < temp_num_bites; ii++ ) {
                util::set_random_numbers_uniform(&rnd_nmbr, 1);
                if (rnd_nmbr < death_probability_infectious) {
                    infectious_composition_array[vv][pp]--;

                    daily_record_of_mosquito_death_infectious++;
                }
            }
        }
    }
}

void VillageManager::mosquito_incubation_step(){

    this->mosquito_incubation_step_kernel(
        this->mosquito_infected_composition_of,
        this->mosquito_infectious_composition_of,
        this->sum_num_villages,
        static_cast<int>(this->parasite_type_lookup.size()),
        this->kMosquito_incubation_probability,

        this->daily_record_of_mosquito_incubation
    );

}
void VillageManager::mosquito_incubation_step_kernel(
        int** infected_composition_array,
        int** infectious_composition_array,
        const int array_size_vv,
        const int array_size_pp,
        const float incubation_probability,

        int& daily_record_of_mosquito_incubation
    ){

    float rnd_nmbr = 0.0;
    int temp_num_bites_to_incubate = 0;

    daily_record_of_mosquito_incubation = 0;

    for (int vv = 0; vv < array_size_vv; vv++) {
        for (int pp = 1; pp < array_size_pp; pp++) { // note that this start at 1, ignoring kNoParasite
            temp_num_bites_to_incubate = infected_composition_array[vv][pp];
            for (int ii = 0; ii < temp_num_bites_to_incubate; ii++ ) {
                util::set_random_numbers_uniform(&rnd_nmbr, 1);
                if (rnd_nmbr < incubation_probability) {
                    infected_composition_array[vv][pp]--;
                    infectious_composition_array[vv][pp]++;

                    daily_record_of_mosquito_incubation++;
                }
            }
        }
    }
}

// void VillageManager::mosquito_ode_update_composition(
//     int v_index, double num_infected, double num_infectious
//     ) {
//     this->mosquito_infected_composition_of[v_index][static_cast<int>(common::ParasiteType::kR0)]
//         = static_cast<int>(num_infected);
//     this->mosquito_infectious_composition_of[v_index][static_cast<int>(common::ParasiteType::kR0)]
//         = static_cast<int>(num_infectious);
// }



float* VillageManager::get_migrant_population_percentages(){
    return this->migrant_population_percentage_of;
}

// void VillageManager::init_villager_mobility(){
// // float randomEvent2 = ran2();
// // if(randomEvent2 <= mobile_pop){
// //     float randomEvent3 = ran2();
// //         if(randomEvent3 <= G_temp_mobile){
// //             G_demes[ind]->mobile = TEMP;}
// //         else{G_demes[ind]->mobile = MOBILE;}
// // }
// // else { G_demes[ind]->mobile = STATIC;}

//     int human_index = 0;
//     for

// }

void VillageManager::init_village_distances(){
    this->init_village_distances_euclidean();
}
void VillageManager::init_village_distances_euclidean(){
    for (int vv1 = 0; vv1 < this->sum_num_villages; vv1++) {
        distances[vv1][vv1] = 0.0;
        for (int vv2 = vv1+1; vv2 < this->sum_num_villages; vv2++) {
            // distances[vv1][vv2] = distances[vv2][vv1] =  sqrt(
            //         pow(longitude_of[vv1]-longitude_of[vv2], 2)
            //         + pow(latitude_of[vv1]-latitude_of[vv2], 2)
            //     );
            distances[vv1][vv2] =  this->get_earth_euclidean_distance (
                    latitude_of[vv1], longitude_of[vv1],
                    latitude_of[vv2], longitude_of[vv2]
                );
            distances[vv2][vv1] = distances[vv1][vv2];
        }
    }
}
float VillageManager::get_earth_euclidean_distance(float lat1, float lng1, float lat2, float lng2) {
    double earthRadius = 6371000; //meters
    double PI=3.14159265359;
    // double dLat = toRadians(lat2-lat1);
    // double dLng = toRadians(lng2-lng1);
    double dLat = (lat2-lat1)*PI/180.0;
    double dLng = (lng2-lng1)*PI/180.0;
    double a = sin(dLat/2) * sin(dLat/2) +
               cos(lat1*PI/180.0) * cos(lat2*PI/180) *
               sin(dLng/2) * sin(dLng/2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    float dist = (float) (earthRadius * c);

    return dist;
}
void VillageManager::update_village_distances_euclidean_filter(const float threshold_ratio){
    float average_distance = 0;
    for (int vv1 = 0; vv1 < this->sum_num_villages; vv1++) {
        for (int vv2 = vv1+1; vv2 < this->sum_num_villages; vv2++) {
            average_distance += distances[vv1][vv2];
        }
    }
    average_distance /= (pow(this->sum_num_villages, 2) - this->sum_num_villages) / 2.0;
    std::cout << "avg. dis:" << std::setprecision(5) << average_distance << std::endl;
    for (int vv1 = 0; vv1 < this->sum_num_villages; vv1++) {
        for (int vv2 = vv1+1; vv2 < this->sum_num_villages; vv2++) {
            if(distances[vv1][vv2] > average_distance*threshold_ratio) {
                distances[vv1][vv2] = distances[vv2][vv1] = std::numeric_limits<float>::infinity();
            }
        }
    }
}
void VillageManager::print_village_distances() const{
    for (int vv1 = 0; vv1 < this->sum_num_villages; vv1++) {
        for (int vv2 = 0; vv2 < this->sum_num_villages; vv2++) {
            std::cout << "[" << std::setw(2) << vv1
                      << "][" << std::setw(2) << vv2
                      << "]" << std::setw(8) << std::setprecision(5)
                      << this->distances[vv1][vv2] << " ";
        }
        std::cout << std::endl;
    } 
}
const float* const* VillageManager::get_village_distances() const{
    return this->distances;
}


// int VillageManager::read_sum_total_population() const{return sum_total_population;}

void VillageManager::output_human_current_locations(
        human::HumanManager& hmn_mgr
    ) const {
    for (int vv = 0; vv < this->sum_num_villages; vv++) {
        for (const auto& hh : this->at_village_register[vv]) {
            hmn_mgr.set_human_at_village(hh, vv);
        }
    }
}


int VillageManager::get_random_human_id_at_village(int village_id) const{
    assert(village_id < this->sum_num_villages);
    int human_index_within_village = util::get_rand_uniform() * this->at_village_register[village_id].size();
    return this->at_village_register[village_id][human_index_within_village];
}

std::vector<int> VillageManager::get_village_ids_ordered_by_tc_desc() const {
    std::vector<int> village_ids;
    for (int vv = 0; vv < this->sum_num_villages; vv++) {
        village_ids.push_back(vv);
    }

    std::stable_sort(
        village_ids.begin(), village_ids.end(),
        [&](int a, int b) { return this->transmission_coefficient_of[a] > this->transmission_coefficient_of[b]; }
    );

    return village_ids;
}

std::vector<int> VillageManager::get_village_ids_ordered_by_tc_ascd() const {
    std::vector<int> village_ids;
    for (int vv = 0; vv < this->sum_num_villages; vv++) {
        village_ids.push_back(vv);
    }

    std::stable_sort(
        village_ids.begin(), village_ids.end(),
        [&](int a, int b) { return this->transmission_coefficient_of[a] < this->transmission_coefficient_of[b]; }
    );

    return village_ids;
}

}