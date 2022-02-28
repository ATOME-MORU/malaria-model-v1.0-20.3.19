#ifndef VILLAGE_H
#define VILLAGE_H
// #include <map>
// #include <list>

#include "common/common.h"
// #include <boost/container/vector.hpp>

namespace human{
    class Human;
    class HumanManager;

    class BloodSystemManager;
}

namespace util{
    class VillageGraph;
}

namespace village {

enum columns{
    kID = 0,
    
    // the order of the following enums must match that of the input data file
    kPopulation = 1, 
    kLatitude, // e.g. Thailand 15.87 N
    kLongitude,// e.g. Thailand 100.9925 E
    kTransmission_coefficient,
    kNum_malaria_post,
    kNum_days_malaria_post_opened,
    kMigrant_population_percentage,
    kArtemisinin_resistance_percentage,

    // columns configured at initialisation
    kInit_first = kPopulation,
    kInit_last = kArtemisinin_resistance_percentage,
};

class Village {

    public:

        const int* id;
        const int* population;
        const float* longitude;
        const float* latitude;
        const float* transmission_coefficient;
        const int* num_malaria_post;
        const int* num_days_malaria_post_opened;
        const float* migrant_population_percentage;
        // const float* artemisinin_resistance_percentage;

};


class VillageManager {

    // Treatment Parameters
    // const float G_timecomp = 1.5;
    // const float G_fullcourse = 0.9;
    // const float G_covab = 0.9;
    // const float G_nomp = 0.1;

    // const float asymptomatic_to_clinical_ratio = 0.0001;

    const int kRegister_reserve_margin = 100;

    const bool kIf_batch_mode;

    const bool kInit_overwrite_ra_rate_if;
    const float kInit_overwirte_ra_rate_with;


    //// Treatment 

    const std::string kFunc_treatment_general_version;
    
    const float G_timecomp;
    const float G_fullcourse;
    const float G_covab;
    const float G_nomp;

    const float asymptomatic_to_clinical_ratio;

    // const float treatment_rate_clinical_has_mp = (1.0 / G_timecomp) * G_fullcourse * G_covab;
    // const float treatment_rate_clinical_no_mp = treatment_rate_clinical_has_mp * G_nomp;

    const float treatment_rate_clinical_has_mp;
    const float treatment_rate_clinical_no_mp;

    const bool kEstablish_malaria_post_during_mda;
    const bool kEstablish_malaria_post_during_survey;

    const int kEstablish_malaria_post_on_day;

    const float kTreatment_rate_mda; // G_tauab = 1.0/1.5;

    const float kTreatment_coverage_mmc_wp2; // mmc_wp2_rev17
    const common::DrugName kTreatment_drug_mmc_wp2;


    //// Infection
    const bool kInit_with_uniform_prevalence_if;
    const float kInit_with_uniform_prevalence_with;

    const float kInit_transmission_coefficient_scaling_factor; // 5.1/0.75
    const float kInit_infected_mosquito_to_infected_human_ratio;
    const float kInit_infectious_mosquito_to_infected_human_ratio;
    // const float kInit_infection_human_to_mosquito_infected_probability = 0.04;
    // const float kInit_infection_human_to_mosquito_infectious_probability = 0.02;
    const float kTransmission_coefficient_to_beta_scaling_factor; // betaf 0.83
    // const float kTransmission_coefficient_to_beta_scaling_factor = 0.25; // betaf 0.83

    
    
    const bool kMosquito_seasonality_switched_on;

    const bool kMosquito_seasonality_init_if_from_file;
    const std::string kMosquito_seasonality_init_file_name; 

    const int kSeasonality_step_offset = 0; //-1, seasons=seasonvec[timestep-1]->vseas[0];


    const float kSeasonality_amplitude;
    const float kSeasonality_amplitude_multiplier;
    // const float kSeasonality_amplitude = 0.6; // G_amp = 0.9
    // const float kSeasonality_amplitude_multiplier = 0.5; // (betas[i]*G_amp*0.5)

    const float kSeasonality_cos_pi_multiplier;
    const float kSeasonality_cos_days_offset;

    const int kFunc_infection_m2h_version;

    const float kInfection_susceptibility_multiplier;
    const float kInfection_infectiousness_multiplier;
    // const float kParasite_mosquito_to_human_probability = 1.0; // G_c = 1.0
    // const float kParasite_human_to_mosquito_probability = 1.0; // G_c = 1.0

    const float kMosquito_biting_rate;

    const float kMosquito_incubation_probability;        //G_mosqinc = 1.0/14.0; // extrinsic incubation period
    const float kMosquito_death_probability_infected;     // G_mosqdeath2 = 1.0/20.0;  //life expectancy of an exposed mosquito
    const float kMosquito_death_probability_infectious;   // G_mosqdeath1 = 1.0/7.0;   //life expectancy of an infectious mosquito
    // const float kMosquito_incubation_probability = 1.0/14.0;        //G_mosqinc = 1.0/14.0; // extrinsic incubation period
    // const float kMosquito_death_probability_infected = 1.0/20.0;     // G_mosqdeath2 = 1.0/20.0;  //life expectancy of an exposed mosquito
    // const float kMosquito_death_probability_infectious = 1.0/7.0;   // G_mosqdeath1 = 1.0/7.0;   //life expectancy of an infectious mosquito


    const float kMosquito_min_filter_threshold;
    const int kMosquito_min_filter_start_on_day;

    std::vector<common::ParasiteType> parasite_type_lookup;
    
    // Mobility

    // G_mobility_static = 125.0/365.0;    // static pop probability of moving to another village
    // G_static_return = 1.0/3.0;         // static pop average length stay in another village
    // G_mobility_mobile = 125.0/365.0;     // seasonal migrant prob of moving to another village
    // G_mobile_return = 1.0/3.0;         // seasonal migrant average length stay in another village
    // G_migrate = 0.0;                   // seasonal migrant prob to migrate to another village other than home place - redefined in function
    // G_migrate_mobile = 1.0/90.0;

    // const float kStatic_population_move_out_probability = 125.0/365.0; //G_mobility_static
    // const float kStatic_population_return_home_probability = 1.0/3.0; //G_static_return
    // const float kNon_static_population_move_out_probability = 125.0/365.0; //G_mobility_mobile
    // const float kNon_static_population_return_home_probability = 1.0/3.0; //G_mobile_return

    const bool kMobility_enabled;
    const float kStatic_population_move_out_probability; //G_mobility_static
    const float kStatic_population_return_home_probability; //G_static_return
    const float kNon_static_population_move_out_probability; //G_mobility_mobile
    const float kNon_static_population_return_home_probability; //G_mobile_return

    std::vector<float> sum_prevalence_infected_of;

    // const float kOde_k;

    std::vector<float> sum_prevalence_of;

    std::vector<std::vector<float>> at_village_attractiveness;
    bool at_village_attractiveness_updated = false;

    void step_update_at_village_attractiveness(); //at_village_attractiveness
    std::vector<int> func_get_humans_to_bite(int num_new_bites, int village_id);

    public:

    // column/attribute-wise data silos
    // data read from input file
    int*    id_of;
    int*    population_of;
    float*  longitude_of;
    float*  latitude_of;
    float*  transmission_coefficient_of;
    int*    num_malaria_post_of;
    int*    num_days_malaria_post_opened_of; // time step mp opens
    float*  migrant_population_percentage_of;
    // float*  artemisinin_resistance_percentage_of;
    float** infection_rates_of; //[num_village][p_types]

    float* mosquito_seasonality_on; // [num_simulation_steps], seasons
    int mosquito_seasonality_on_length = 0;

    // std::vector< std::vector<float> > mosquito_bite_probability_of_on; // [num_village][num_simulation_steps]
    std::vector<float> beta_season_of; // [num_villages] beta_season per village, needs update every time step
    int beta_season_of_last_updated_on; // the last time_step in which beta_season_of was updated.

    int** sum_prevalence_count_of; //[num_village][p_types]

    // derived consts
    float** distances; // earth surface straight line distances

    float** distances_shortest_path;
    std::vector<std::vector<int>> static_mobility_network;

    // std::map<int, std::list<human::Human*>> at_village_register;
    typedef std::vector<int> Register_row_t;
    typedef std::vector< Register_row_t > Register_t;
    // std::vector<std::vector<int>> at_village_register;
    // boost::container::vector<boost::container::vector<int>> at_village_register;
    Register_t at_village_register;


    int** mosquito_infected_composition_of; //[num_village][p_types]
    int** mosquito_infectious_composition_of; //[num_village][p_types]

    // 
    bool*   malaria_post_in_operation_of;
    float*  treatment_rate_mp_clinical_of;
    float*  treatment_rate_mp_asymptomatic_of;



    // row/object-wise data points/registers
    Village* village_reg;


    human::HumanManager* v_hmn_mgr;

    // TODO: move public: to here

    // summary data

    int daily_record_of_new_bites_sum=0;
    int daily_record_of_new_infections_sum=0;
    int daily_record_of_wasted_bites_sum=0;

    std::vector<int> daily_record_of_new_bites;
    std::vector<int> daily_record_of_new_infections;
    std::vector<int> daily_record_of_wasted_bites;

    std::vector<int> daily_record_of_new_bites_annual;
    std::vector<int> daily_record_of_new_infections_annual;
    std::vector<int> daily_record_of_new_clinical_cases_annual;
    std::vector<int> daily_record_of_new_asymptomatic_cases_annual;

    int daily_record_of_mosquito_newly_infected = 0;
    int daily_record_of_mosquito_death_infected = 0;
    int daily_record_of_mosquito_death_infectious = 0;
    int daily_record_of_mosquito_incubation = 0;
    int daily_record_of_mosquito_infected_sum = 0;
    int daily_record_of_mosquito_infectious_sum = 0;
    int daily_record_of_mosquito_clean_sum = 0;
    int daily_record_of_mosquito_new_adults_sum = 0;
    int daily_record_of_mosquito_death_adults_sum = 0;
    int daily_record_of_mosquito_death_adults_infectious_count_sum = 0;

    int daily_record_of_eggs_sum = 0;
    int daily_record_of_larv_sum = 0;
    int daily_record_of_imat_sum = 0;
    int daily_record_of_aliv_sum = 0;
    int daily_record_of_aliv_biting_sum = 0;

    int daily_record_of_max_capacity_sum = 0;

    int daily_record_of_mosquito_new_rb = 0;

    std::vector<int> daily_record_of_drug_failure_by_new_treatment{std::vector<int>(static_cast<int>(common::DrugName::Last) + 1, 0)};

    std::vector<int> daily_record_of_new_treatments_by_parasite_type{std::vector<int>(static_cast<int>(common::ParasiteType::Last) + 1, 0)};

    int sum_num_villages = 0;
    int sum_total_population = 0;

    VillageManager(

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

        //Mobility
        const bool mobility_enabled,
        const float static_population_move_out_probability,
        const float static_population_return_home_probability,
        const float non_static_population_move_out_probability,
        const float non_static_population_return_home_probability

    );

    ~VillageManager();

    float* get_migrant_population_percentages();

    void update_prevalence_count_of(const common::ParasiteType* dominant_type_array);
    void update_prevalence_of();
    inline const std::vector<float>& get_prevalence_of() const {
        return this->sum_prevalence_of;
    }
    inline float get_prevalence_of(int village_id) const {
        return this->sum_prevalence_of[village_id];
    }
    // int** get_prevalence_count_of() const;
    std::vector<int> get_prevalence_greater_than(float prevalence_threshold) const;
    std::vector<int> survey_prevalence_greater_than(
        float prevalence_threshold,
        int pcr_sample_size,
        float pcr_sensitivity
    );
    float get_prevalence_of_village_by_survey(
        int village_id,
        int pcr_sample_size,
        float pcr_sensitivity
    );

    // inline const std::vector<float>& get_sum_prevalence_infected_of() const {
    //     return this->sum_prevalence_infected_of;
    // }
    const std::vector<float>& get_sum_prevalence_infected_of(
        const common::ParasiteType* dominant_type_array,
        const float* infectiousness_array
    );

    void update_summary();  // Updates summary data
    void print_summary() const;   // Prints summary data
    void print_all() const;       // Prints all village data
    void print_one(const int index_village) const;

    void init_treatment_rates_malaria_post();
    void update_treatment_rates_malaria_post(int village_id);
    void establish_malaria_posts();
    void open_malaria_post(int village_id);
    void treatment_by_malaria_post(
            human::BloodSystemManager* bld_mgr,
            human::HumanManager* hmn_mgr, // to remove
            int time_step
         );
    static void treatment_by_malaria_post_kernel(
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
        );

    void treatment_by_mda(
            const int village_id,
            const common::DrugName mda_drug
    );
    static void treatment_by_mda_kernel(
        // const common::DrugName* existing_drug,
        common::DrugName* given_drug,
        const int human_index_max,
        // const std::vector<int>& list_of_humans,
        // const boost::container::vector<int>& list_of_humans,
        const Register_row_t& list_of_humans,
        const float* rnd_nmbrs,
        const float treatment_rate,
        const common::DrugName mda_drug
    );

    void step_treatment_general(
        human::BloodSystemManager* bld_mgr,
        human::HumanManager* hmn_mgr,
        int time_step
    );
    void func_general_treatment_mmc_wp2(
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
        std::vector<int>& drug_failure_register,
        std::vector<int>& daily_record_of_new_treatments_by_parasite_type
    );
    void func_general_treatment_mmc_wp2_rev1(
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
        std::vector<int>& drug_failure_register,
        std::vector<int>& daily_record_of_new_treatments_by_parasite_type
    );

    void init_get_villagers(human::HumanManager* human_manager);

    void init_infections(
            common::ParasiteType* infections_array,
            int infections_array_size
    );
    void init_infections_kernel(
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
            const float init_infected_human_to_infected_mosquito_ratio,
            const float init_infected_human_to_infectious_mosquito_ratio
            // const float infected_probability,
            // const float infectious_probability
    );

    // void init_mosquito_bite_probability_of_on();
    void step_update_beta_season_of(int at_time_step);
    inline std::vector<float>& get_beta_season_of() {
        return this->beta_season_of;
    };
    // static int calculate_number_of_new_bites_in_village(
    //     float prv_dist_base_value,
    //     float prv_dist_base_to_trns_cff_scaling_factor,
    //     float seasonality_factor,
    //     float seasonality_amplitude,
    //     float seasonality_amplitude_factor,
    //     uint population,
    //     int time_step
    // );

    

    void infection_m2h(
        common::ParasiteType* infections_array,
        const int infections_array_size,
        const float* susceptibility_array,
        const int susceptibility_array_size,
        const common::StageName* stage_array,
        const int time_step,
        const std::vector<double>& num_infectious_mosq_from_ode,
        // const double ode_feed_rate,
        const bool if_using_ode
    );
    // static void infection_m2h_kernel(
    void infection_m2h_kernel(
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

        common::ParasiteType* infections_array,
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
    );
    static void infection_m2h_kernel_use_all_bites(
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
    );
    void infection_m2h_kernel_per_bite(

        const int time_step,
        const float mosquito_min_filter_threshold,
        const int mosquito_min_filter_start_on_day,
        
        const int vll_array_size,
        // const Register_t& villager_reg,
        const int * const* infectious_msq_composition_array,

        const float bite_success_probability, // susceptibility

        const float* susceptibility_array,
        const int susceptibility_array_size,
        const common::StageName* stage_array,

        common::ParasiteType* infections_array,//[h]
        const int infections_array_size,

        int& daily_record_of_new_bites_sum,
        int& daily_record_of_new_infections_sum,
        int& daily_record_of_wasted_bites_sum,

        std::vector<int>& daily_record_of_new_bites,
        std::vector<int>& daily_record_of_new_infections,
        std::vector<int>& daily_record_of_wasted_bites,

        const std::vector<double>& num_infectious_mosq_from_ode,
        const double ode_feed_rate,
        const bool if_using_ode
    );

    void infection_h2m(
        const float* infectiousness_array,
        const common::ParasiteType* dominant_type_array,
        const int host_array_size,
        const int time_step
    );
    void infection_h2m_kernel(
        // const float* transmission_coefficient_array,
        // const float seasonality_factor,
        
        const int vll_array_size,
        // const std::vector<std::vector<int>> villager_reg,
        // const boost::container::vector<boost::container::vector<int>> villager_reg,
        const Register_t& villager_reg,
        int** infected_msq_composition_array,

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

        int& daily_record_of_mosquito_newly_infected,
        int& daily_record_of_mosquito_new_rb
    );

    void init_mosquito_seasonality_from_file(std::string file_name);

    void mosquito_survival_step();
    void mosquito_survival_step_kernel(
        int** infected_composition_array,
        int** infectious_composition_array,
        const int array_size_vv,
        const int array_size_pp,
        const float death_probability_infected,
        const float death_probability_infectious,

        int& daily_record_of_mosquito_death_infected,
        int& daily_record_of_mosquito_death_infectious
    );

    void mosquito_incubation_step();
    static void mosquito_incubation_step_kernel(
        int** infected_composition_array,
        int** infectious_composition_array,
        const int array_size_vv,
        const int array_size_pp,
        const float incubation_probability,

        int& daily_record_of_mosquito_incubation
    );

    // void mosquito_ode_update_composition(
    //     int v_index, double num_infected, double num_infectious
    // );

    void init_village_distances();
    void init_village_distances_euclidean();
    float get_earth_euclidean_distance(float lat1, float lng1, float lat2, float lng2);

    void update_village_distances_euclidean_filter(const float threshold_ratio);

    void print_village_distances() const;
    const float* const* get_village_distances() const;


    // void give_drug_to_village(human::HumanManager* human_manager,
    //                          const int index_village,
    //                          const common::DrugName drug);
    // void give_drug_to_village(
    //             const int index_village,
    //             const common::DrugName drug
    //         );
    // int read_sum_total_population() const;


    // Mobility
    void init_shortest_path_distances(const util::VillageGraph &vg);
    std::string output_shortest_path_from_village_gnuplot(
        const std::string file_prefix,
        const int village_id,
        const util::VillageGraph &vg
    ) const;

    void build_static_mobility_network();

    void step_population_movement_static(
        const int human_array_size,
        const common::HumanMobilityType* human_mobility_type_array,
        const int* human_home_village_array,

        const float* random_numbers_array,
        const int random_numbers_array_size
    );
    void step_population_movement_static_kernel(
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
    );
    // static void remove_villager(Register_t& villager_reg, int vv);
    void villager_reg_clean_up(
        const int village_array_size,
        // std::vector<std::vector<int>>& villager_reg
        // boost::container::vector<boost::container::vector<int>>& villager_reg
        Register_t& villager_reg
    );


    // functions
    static float func_beta_to_beta_season_cos(
        float beta,
        float season_amplitude,
        float cos_pi_multiplier,
        float cos_days_offset,
        int time_step
    );
    static float func_beta_to_beta_season_given_multiplier(
        float beta,
        float season_multiplier,    // from seasonality file
        float season_amplitude,     // G_amp
        float season_amplitude_multiplier // 0.5
    );
    static int func_beta_season_to_num_bites(
        float beta_season,
        int population
    );
    static int func_num_infectious_mosq_to_num_bites(
        double num_infectious_mosq,
        float ode_k
    );
    static int func_num_infectious_mosq_to_num_bites_uniform(
        double num_infectious_mosq,
        float ode_k
    );
    static int func_num_infectious_mosq_to_num_bites_uniform(
        int num_infectious_mosq,
        float ode_k
    );
    
    void output_human_current_locations(
        human::HumanManager& hmn_mgr
    ) const;

    int get_random_human_id_at_village(int village_id) const;

    std::vector<int> get_village_ids_ordered_by_tc_desc() const;
    std::vector<int> get_village_ids_ordered_by_tc_ascd() const;

    ///////////////////////////////////////////
    // 

    inline void step_reset_updated_indicators() {
        at_village_attractiveness_updated = false; // updated after mobility step
    }

    ///////////////////////////////////////////
    // Access

    inline const float* get_transmission_coefficient_of() const {
        return this->transmission_coefficient_of;
    }

    inline const int* get_population_of() const {
        return this->population_of;
    }
    inline int get_home_population_of_village(int village_id) const {
        return this->population_of[village_id];
    }
    inline int get_current_population_of_village(int village_id) const {
        return static_cast<int>(this->at_village_register[village_id].size());
    }

};

}

#endif