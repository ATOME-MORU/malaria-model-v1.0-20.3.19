#ifndef BLOOD_SYSTEM_H
#define BLOOD_SYSTEM_H


#include "common/common.h"


namespace human{

// number of preallocated parasite groups per blood system
constexpr int kNum_pgs_per_system = 8; // TODO: put this inside BloodSystemManager

// class BloodSystem {

//     public:
//         bool* immunity;
//         int* immunity_days;
// };
    
class BloodSystemManager {

    float* susceptibility_of;
    float* infectiousness_of;

    bool* is_clinical_of;
    uint8_t* moi_of;    //multiplicity of infection
    common::ParasiteType* dominant_type_of;
    common::StageName* most_advanced_stage_of;

    int* immunity_days_of; // "-1" no immunity, ">=0" # of days immune
    uint8_t* immunity_level_of; // 0 - no immunity
    uint8_t* cummulative_exposures_of;

    common::DrugName* drug_of; // drug which is currently in host
    // int* drug_days_of;
    common::DrugName* drug_given_of; // drug last given to the host
    int* drug_given_day_of; // day on which the host last received a drug 

    int* last_recovery_day_of; // day on which the last recovered (natural/drug)

    uint8_t* num_of_clinical_pgs_of; // number of infections currently in the host that caused clinical-ness

    //TODO: bool* has_prim;

    // parasite group ("pg") related data
    common::ParasiteType* pg_type_of; // parasite type
    common::StageName* pg_stage_of;   // liver/blood/infectious
    // int* pg_stage_progression_due_day_of; // -1 if newly
    uint16_t* pg_aquired_on_day_of;   // time stamp of when the pg entered the b_system
    // uint16_t* pg_infectious_on_day_of;// time stamp of when the pg became infectious
    uint32_t* pg_num_parasites_of;    // number of parasites
    uint32_t* pg_num_gametocytes_of;  // number of gametocytes
    bool* pg_is_drug_target_of;
    bool* pg_has_caused_clinicalness_of;

    int8_t* head_pg_of; 
    int8_t* tail_pg_of;

    // int8_t range -128 ~ 127
    // uint8_t range 0 ~ 255
    // uint16_t range 0 ~ 65535
    // uint32_t range 0 ~ 4,294,967,295

     // G_antroph = 0.95;
     // G_indoorpref = 0.95;
     // G_xitn =(0.0)/365.0;   //change
     // G_itneffect = 0.5*G_antroph*G_indoorpref;
     // G_itnprop = 0.2;

    const bool kTreat_clean_hosts;
    const bool kTreat_asymptomatic_hosts;
    const float kInit_clinical_probability;
    // G_initclin = 0.02;

    const float kSusceptibility_default;
    const int kFunc_update_susceptibility_version;

    const float kInfectiousness_default;
    // const float kSusceptibility_default = 1.0;
    // const float kInfectiousness_default = 1.0;

    const float kItn_effect;
    // const float kItn_effect = 0.5 * 0.95 * 0.95;
     // G_antroph = 0.95;
     // G_indoorpref = 0.95;
     // G_itneffect = 0.5*G_antroph*G_indoorpref;
    // const float kItn_effect = 0.5 * 0.2 * 0.95;



    const float kPhi_a; //G_phia = 0.5
    const float kPhi_c; //G_phic = 1.0
    // const float kPhi_a = 1.0; //G_phia = 0.5
    // const float kPhi_a = 0.5; //G_phia = 0.5
    // const float kPhi_c = 1.0; //G_phic = 1.0

    const float kProphylaxis_effect;// G_artprophylaxis = 1.0
    // const float kProphylaxis_effect = 1.0;// G_artprophylaxis = 1.0

    const int kFunc_dominance_version;

    const int kFunc_clinical_probability_blood_to_infectious_version;
    const float kFunc_clinical_probability_blood_to_infectious_version3_para1;
    const float kFunc_clinical_probability_blood_to_infectious_version3_para2;

    // clinical resolution
    const float kMellow_rate; // G_mellow = 1.0/10.0;
    // const float kMellow_rate = 1.0/3.0; // G_mellow = 1.0/10.0;
    // const float kMellow_rate = 0.1; // G_mellow = 1.0/10.0;

    const float kStage_progression_l2b_days_mean;
    const float kStage_progression_l2b_days_sd;
    const float kStage_progression_b2i_days_mean;
    const float kStage_progression_b2i_days_sd;
    const float kStage_progression_i2r_days_mean;
    const float kStage_progression_i2r_days_sd;

    // const float kStage_progression_i2r_days_mean;
    // const float kStage_progression_i2r_days_sd;

    // const float kNatural_recovery_rate; // G_delta = 1.0/80.0
    // const float kNatural_recovery_rate = 1.0/160.0; // G_delta = 1.0/80.0
    // const float kNatural_recovery_rate = 1.0/80.0; // G_delta = 1.0/80.0
    const float kNatural_recovery_resistance_cost; // fitness cost

    const int kDrug_failure_if_clinical_on_day;
    // const float kDrug_b_resistance;
    // const float kDrug_ab_resistance;
    const bool kAllow_act_to_a;
    const bool kPossible_to_have_drug_a_only;

    const common::DrugName kDrug_ab = common::DrugName::kPiparte;
    const common::DrugName kDrug_a = common::DrugName::kArtesunate;
    const common::DrugName kDrug_b = common::DrugName::kPiperaquine;

    const float kMutation_prob_single_resistance;
    const float kMutation_prob_single_resistance_w_act;
    const float kMutation_prob_single_type_a_type_b_ratio;
    const bool kMutation_type_a_enabled;
    const bool kMutation_type_b_enabled;
    // const float kMutation_prob_double_resistance;

    const int kFunc_parasite_mutation_version;
    const int kFunc_parasite_mutation_fixation_after_day;

    // immunity-related

    const int kImmunity_days_until_loss_start;
    const int kImmunity_level_until_loss_start; // **
    const float kImmunity_loss_rate; //G_alpha = 1.0/60.0
    // const int kImmunity_days_until_loss_start = 40;
    // const int kImmunity_level_until_loss_start = 1; // **
    // const float kImmunity_loss_rate = 1.0/60.0; //G_alpha = 1.0/60.0


float kPds_kill_probability[static_cast<int>(common::ParasiteType::Last)+1]
                     [static_cast<int>(common::DrugName::Last)+1]
                     [static_cast<int>(common::StageName::Last)+1]= {
    // no parasite
    {
        //no drug
        {0,0,0,0},
        // Artesunate
        {0,0,0,0},
        // piperaquine
        {0,0,0,0},
        // piparte
        {0,0,0,0},
    // C. AQ
        {0,0,0,0}, // mmc_wp2_rev17
    // A+C. ASAQ
        {0,0,0,0}, // mmc_wp2_rev17
        // lumefantrine
        {0,0,0,0},
        // al
        {0,0,0,0},
        // SP
        {0,0,0,0},
        // chloroquine
        {0,0,0,0}
    },
    // parasite ra
    {
        //no drug
        {0,0,0,0},
        // Artesunate
        {0,0,1.0/5.0*0.27,              1.0/3.0*0.27},
        // piperaquine
        {0,0,1.0/3.0,                   1.0/15.0},
        // piparte
        {0,0,1.0/5.0*0.27+1.0/3.0*0.73, 1.0/3.0*0.27+1.0/15.0*0.73},
    // C. AQ
        {0,0,0,0}, // mmc_wp2_rev17
    // A+C. ASAQ
        {0,0,0,0}, // mmc_wp2_rev17
        // lumefantrine
        {0,0,0,0},
        // al
        {0,0,0,0},
        // SP
        {0,0,1.0/5.0,                   1.0/3.0},
        // chloroquine
        {0,0,1.0/5.0,                   1.0/3.0}
    },
    // parasite rb
    {
        //no drug
        {0,0,0,0},
        // Artesunate
        {0,0,1.0/5.0,                   1.0/3.0},
        // {0,0,0.0,                   0.0},
        // piperaquine
        {0,0,1.0/3.0*0.8,               1.0/15.0*0.8},
        // {0,0,1.0/3.0*0.0,               1.0/15.0*0.0},
        // piparte
        // {0,0,1.0/5.0*0.8+1.0/5.0*0.2,   1.0/3.0*0.8+1.0/3.0*0.2},
        {0,0,1.0/3.0*0.8+1.0/5.0*0.2,   1.0/15.0*0.8+1.0/3.0*0.2},
        // {0,0,0.0,   0.0},
    // C. AQ
        {0,0,0,0}, // mmc_wp2_rev17
    // A+C. ASAQ
        {0,0,0,0}, // mmc_wp2_rev17
        // lumefantrine
        {0,0,0,0},
        // al
        {0,0,0,0},
        // SP
        {0,0,1.0/5.0,                   1.0/3.0},
        // chloroquine
        {0,0,1.0/5.0,                   1.0/3.0}
    },
    // parasite r0
    {
        //no drug
        {0,0,0,0},
    // A. Artesunate  
        {0,0,1.0/5.0,       1.0/3.0},  
    // B. Piperaquine
        // {0,0,1.0/3.0,       1.0/15.0}, // mda_coverage, mosq_approx
        // {0,0,1.0/3.0,       1.0/10.0}, // mmc_wp2_rev5, 5562e42
        // {0,0,1.0/3.0,       1.0/7.0}, // mmc_wp2_rev14a
        // {0,0,1.0/7.0,       1.0/7.0}, // mmc_wp2_rev16a, A4
        // {0,0,1.0/8.0,       1.0/8.0}, // mmc_wp2_rev17, A5, AQ
        {0,0,1.0/9.0,       1.0/9.0}, // mmc_wp2_rev17, A6, L
    // A+B. Piparte
        // {0,0,1.0/5.0,       1.0/3.0}, // mda_coverage, mosq_approx
        // {0,0,1.0/5.0,       1.0/2.0}, // mmc_wp2_rev14a
        // {0,0,1.0/2.0,       1.0/2.0}, // mmc_wp2_rev16a, A4
        // {0,0,1.0/2.0,       1.0/2.0}, // mmc_wp2_rev17, A5, ASAQ
        {0,0,1.0/3.0,       1.0/3.0}, // mmc_wp2_rev17, A6, AL
    // C. AQ
        {0,0,1.0/8.0,       1.0/8.0}, // mmc_wp2_rev17
    // A+C. ASAQ
        {0,0,1.0/2.0,       1.0/2.0}, // mmc_wp2_rev17
    // D. Lumefantrine
        {0,0,1.0/9.0,       1.0/9.0}, // mmc_wp2_rev17
    // A+D. AL
        {0,0,1.0/3.0,       1.0/3.0}, // mmc_wp2_rev17
    // SP
        {0,0,1.0/5.0,       1.0/3.0},
    // Chloroquine
        {0,0,1.0/5.0,       1.0/3.0}
    },
    // parasite rab
    {
        //no drug
        {0,0,0,0},
        // Artesunate
        {0,0,1.0/5.0,                   1.0/3.0},
        // {0,0,0.0,                   0.0},
        // piperaquine
        {0,0,1.0/3.0*0.8,               1.0/15.0*0.8},
        // {0,0,1.0/3.0*0.0,               1.0/15.0*0.0},
        // piparte
        // {0,0,1.0/5.0*0.8+1.0/5.0*0.2,   1.0/3.0*0.8+1.0/3.0*0.2},
        {0,0,1.0/3.0*0.8+1.0/5.0*0.2,   1.0/15.0*0.8+1.0/3.0*0.2},
        // {0,0,0.0,   0.0},
    // C. AQ
        {0,0,0,0}, // mmc_wp2_rev17
    // A+C. ASAQ
        {0,0,0,0}, // mmc_wp2_rev17
        // lumefantrine
        {0,0,0,0},
        // al
        {0,0,0,0},
        // SP
        {0,0,1.0/5.0,                   1.0/3.0},
        // chloroquine
        {0,0,1.0/5.0,                   1.0/3.0}
    },
};

    // Export Format
    const std::string kExport_file_name = "_agestruct.csv";
    const char kData_sep_char = ';';
    const char kMoi_liver_char = 'L';
    const char kMoi_blood_char = 'B';
    // const char kMoi_infectious_target_char = 'D';
    const char kMoi_infectious_clinical_char = 'C';
    const char kMoi_infectious_asymptomatic_char = 'A';
    const char kMoi_clean_char = 'O';
    const char kMoi_unknown_char = 'U';

public:

    int sum_num_systems = 0;
    int sum_num_pgs = 0;

    int sum_num_pg_q_is_full_incidents = 0;
    int sum_num_pg_intakes = 0;

    const bool kDebug_assert = true;
    const bool kDebug_print = (kDebug_assert && false);
    std::vector<int> debug_bites_accumulator_add;
    std::vector<int> debug_bites_accumulator_less;

    // metrics that are reported daily
    int sum_prevalence_count[static_cast<int>(common::ParasiteType::Last)+1] = {};
    // tf-treatment failure
    int sum_tf_by_drug_loss_count[static_cast<int>(common::ParasiteType::Last)+1] = {};
    int sum_clinical_count = 0;
    int sum_asymptomatic_count = 0;

    std::vector<int> daily_record_of_drug_failure_by_effetive_period{std::vector<int>(static_cast<int>(common::DrugName::Last)+1,0)} ;
    std::vector<int> daily_record_of_drug_intake{std::vector<int>(static_cast<int>(common::DrugName::Last)+1,0)} ;

    std::vector<int> daily_record_of_drug_current{std::vector<int>(static_cast<int>(common::DrugName::Last)+1,0)} ;
    std::vector<int> daily_record_of_drug_current_in_infected_population{std::vector<int>(static_cast<int>(common::DrugName::Last)+1,0)} ;

    std::vector<float> daily_record_of_genotype_weights{std::vector<float>(static_cast<int>(common::ParasiteType::Last)+1,0)};
    int daily_record_of_parasite_positive_systems = 0;

    int daily_record_of_num_pgs = 0;
    float daily_record_of_avg_clinical_pgs = 0.0;

    std::vector<int> daily_record_of_new_cases_clinical{std::vector<int>(static_cast<int>(common::ParasiteType::Last)+1,0)};
    std::vector<int> daily_record_of_new_cases_asymptomatic {std::vector<int>(static_cast<int>(common::ParasiteType::Last)+1,0)};


    std::vector<int> record_of_cleared_pgs_r0{std::vector<int>(200, 0)};
    std::vector<int> record_of_cleared_pgs_rb{std::vector<int>(200, 0)};

    std::vector<int> log_new_fixation{std::vector<int>(static_cast<int>(common::ParasiteType::Last)+1,0)};
    // These are per pg logs (log_new_mutation is per host in parasite_mutation_mmc_wp2_drug_days_fixation)
    std::vector<int> log_new_mutation{std::vector<int>(static_cast<int>(common::ParasiteType::Last)+1,0)};
    std::vector<int> log_new_clearance_natural{std::vector<int>(static_cast<int>(common::ParasiteType::Last)+1,0)};
    std::vector<int> log_new_clearance_drug{std::vector<int>(static_cast<int>(common::ParasiteType::Last)+1,0)};

    int num_pip_treatments = 0;
    int num_pip_treatments_success = 0;


    //reporting items

    bool* report_stage_progression_flags;

    std::vector<bool> report_daily_new_clinical_case; // size H
    std::vector<bool> report_daily_new_asymptomatic_case; // size H
    
    // std::list<common::ParasiteGroup> parasite_reg;
    // int* parasite_offset;

    BloodSystemManager(
        int init_number_of_systems,
        
        float itn_mosq_preference_indoor,
        float itn_mosq_preference_human,
        float itn_effectiveness,

        bool treat_clean_hosts,
        bool treat_asymptomatic_hosts,

        float init_clinical_probability,

        float susceptibility_default,
        int func_update_susceptibility_version,

        float infectiousness_default,

        float phi_a,
        float phi_c,

        float prophylaxis_effect,

        int func_dominance_version,

        int func_clinical_probability_blood_to_infectious_version,
        float func_clinical_probability_blood_to_infectious_version3_para1,
        float func_clinical_probability_blood_to_infectious_version3_para2,

        float mellow_days,

        float stage_progression_l2b_days_mean,
        float stage_progression_l2b_days_sd,
        float stage_progression_b2i_days_mean,
        float stage_progression_b2i_days_sd,
        float stage_progression_i2r_days_mean,
        float stage_progression_i2r_days_sd,

        // float natural_recovery_days,
        float natural_recovery_resistance_cost,

        int drug_failure_if_clinical_on_day,
        bool allow_act_to_a,
        bool possible_to_have_drug_a_only,
        // float drug_b_resistance,
        // float drug_ab_resistance,

        float base_effect_multiplier,

        bool resistance_multipliers_enabled,
        float drug_b_effect_on_rb,
        float drug_ab_effect_on_rb,
        float drug_a_effect_on_ra,
        float drug_ab_effect_on_ra,
        float drug_a_effect_on_rab,
        float drug_b_effect_on_rab,
        float drug_ab_effect_on_rab,

        float mutation_prob_single_resistance,
        float mutation_prob_single_resistance_w_act,
        float mutation_prob_single_type_a_type_b_ratio,
        bool mutation_type_a_enabled,
        bool mutation_type_b_enabled,
        // float mutation_prob_double_resistance,

        int func_parasite_mutation_version,
        int func_parasite_mutation_fixation_after_day,

        int immunity_days_until_loss_start,
        int immunity_level_until_loss_start,
        float immunity_loss_days
    );

    void reset_system_at(int system_index);
    void birth_and_death(bool* flags, int num_systems);
    ~BloodSystemManager();

    void init_for_test();
    void init_for_drug_failure_test(
        common::ParasiteType parasite_type,
        common::DrugName drug_to_test
    );
    void init_from_file(
        std::string input_file_name,
        char input_file_delimiter,
        human::HumanManager* human_manager,
        bool if_use_age_column
    );

    void export_to_file(
        const std::string output_file_prefix,
        const human::HumanManager* human_manager
    );

    void import_from_file(
        const std::string file_name,
        human::HumanManager* human_manager,
        const bool if_use_age_column
    );

    int add_one_pg_with_defaults(
        const int sys_index,
        const common::ParasiteType pg_type,
        const common::StageName pg_stage,
        const uint16_t pg_aquired_on_day
        // const uint32_t pg_num_parasites,
        // const uint32_t pg_num_gametocytes,
        // const bool pg_is_drug_target,
        // const bool pg_has_caused_clinicalness
    );

    void update_susceptibility(
        const bool* if_itn,
        int if_itn_length,
        const int* age_of
    );
    static void update_susceptibility_kernel_1(
        const int* age_array,
        const bool* if_itn_array,
        const bool* if_clinical_array,
        const common::DrugName* drug_array,
        const int array_size,
        float* susceptibility_array,
        float* infectiousness_array,
        // const float default_susceptibility,
        const float default_infectiousness,
        const float phi_a,
        const float phi_c,
        const float prophylaxis_effect,
        const float itn_effect
    );
    static void update_susceptibility_kernel(
        const bool* if_itn_array,
        const bool* if_clinical_array,
        const common::DrugName* drug_array,
        const int array_size,
        float* susceptibility_array,
        float* infectiousness_array,
        const float default_susceptibility,
        const float default_infectiousness,
        const float phi_a,
        const float phi_c,
        const float prophylaxis_effect,
        const float itn_effect
    );

    void clinical_resolution();
    static void clinical_resolution_kernel(
            bool* is_clinical_array,
            int array_size,
            float mellow_rate
        );

    void parasite_intake(
        common::ParasiteType* intake_list,
        const common::StageName begin_at_stage,
        const int time_step,
        const bool if_intake_for_init
    );
    void parasite_injection(
        // common::ParasiteType* parasite_injection_tray,
        const common::ParasiteType parasite_type_to_inject,
        const int num_cases
    );
    void drug_intake(
        const common::DrugName* drug_list,
        const int time_step
    );
    void drug_loss(
        float* random_numbers,
        int random_numbers_array_size
    );

    void drug_effect(
        const int time_step
    );
    
    // void drugs_meet_parasites();
    void step_report_drug_failure_by_clinical_status_on_day(
        int current_time_step
    );

    void step_report_drug_composition();

    void parasite_mutation(int current_time_step);
    void parasite_mutation_mmc_wp2_drug_days_fixation_per_site(int current_time_step);
    void parasite_mutation_mmc_wp2_drug_days_fixation(int current_time_step);
    void parasite_mutation_mmc_wp2();
    void parasite_mutation_mmc_wp2_drug_trigger();
    void parasite_mutation_not_used(uint16_t time_step);

    void parasite_stage_progression(
        const int time_step
    );
    bool parasite_stage_progression_probability_based(
        const int p_index,
        const int s_index
        );
    float func_clinical_probability_blood_to_infectious(
        uint8_t immunity_level,
        uint8_t cummulative_exposures,
        uint8_t moi
    );
    static float func_clinical_probability_blood_to_infectious_kernel_0(
        uint8_t immunity_level,
        uint8_t cummulative_exposures,
        uint8_t moi
    );
    static float func_clinical_probability_blood_to_infectious_kernel_1(
        uint8_t immunity_level,
        uint8_t cummulative_exposures,
        uint8_t moi
    );
    static float func_clinical_probability_blood_to_infectious_kernel_2(
        uint8_t immunity_level,
        uint8_t cummulative_exposures
    );
    static float func_clinical_probability_blood_to_infectious_kernel_3(
        uint8_t immunity_level,
        uint8_t cummulative_exposures,
        float para1,
        float para2
    );
    static float func_clinical_probability_blood_to_infectious_kernel_4(
        uint8_t immunity_level,
        uint8_t cummulative_exposures,
        uint8_t moi
    );

    void parasite_dominance();
    static void parasite_dominance_mmc_wp2_kernel(
        const int head_pg_index,
        const int tail_pg_index,
        const int system_index,
        common::ParasiteType* dominant_type_of_array,
        common::StageName* most_advanced_stage_of_array,
        const common::StageName* pg_stage_of_array,
        const common::ParasiteType* pg_type_of_array,
        const common::DrugName* drug_of_array,
        const common::DrugName drug_ab,
        const common::DrugName drug_a,
        const common::DrugName drug_b
    );
    void parasite_dominance_probability_based(
        const int head_pg_index,
        const int tail_pg_index,
        const int system_index
        );
    void step_report_genotype_weights();

    template <typename T>
    static common::ParasiteType parasite_dominance_density_based(
        const int head_pg_index,
        const int tail_pg_index,
        const T* densities,
        const common::ParasiteType* types
        );

    void step_immunity_loss(
        const float* uniform_random_number_array,
        const int uniform_random_number_array_size
    );
    static void step_immunity_loss_kernel(
        uint8_t* immunity_level_array,
        int* immunity_days_array,
        const int immunity_array_size,

        const int immunity_days_until_loss_start,
        const int immunity_level_until_loss_start,
        const float immunity_loss_rate, 

        const float* uniform_random_number_array,
        const int uniform_random_number_array_size
    );

    // const get functions
    float* get_susceptibility() const;
    float* get_infectiousness() const;
    common::StageName* get_most_advanced_stage() const;
    common::DrugName* get_drug_of() const;
    common::DrugName* get_drug_given_of() const;
    int* get_drug_given_day_of() const;
    bool* get_is_clinical_of() const;
    common::ParasiteType* get_dominant_type_of() const;
    common::ParasiteType* get_pg_type() const;
    common::StageName* get_pg_stage() const;

    uint8_t* get_immunity_level_of() const;
    uint8_t* get_cummulative_exposures_of() const;

    uint64_t get_parasite_count() const;

    inline int get_system_id_of_pg(int pg_index) const{
        return pg_index / kNum_pgs_per_system;
    }

    float get_mellow_rate() const;
    int get_drug_failure_if_clinical_on_day() const;

    void print_all() const;
    void print_one(const int index_bsys) const;

    void print_drug_parameters() const;

};

}

#endif