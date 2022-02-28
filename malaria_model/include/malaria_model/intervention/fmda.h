#ifndef FMDA_H
#define FMDA_H


namespace intervention {

class Fmda{
public:
    // Mda parameters
    const int kNum_teams;
    const int kNum_villages;

    const int kNum_days_spent_in_village;
    const int kNum_drug_rounds;
    const int kNum_days_between_drug_rounds;
    

    // Fmda specific parameters
    const float* const* kVillage_distances;

    const int kSurvey_data_delivery_delay;
    const float kPcr_sensitivity;
    const int kPcr_sample_size;
    const float kPrevalence_threshold;

    const std::string kOutput_prefix;

    const bool kPrint_map;

    std::vector<bool> latest_survey_results;

    std::vector<int> scheduled_start_times;

    std::vector<std::vector<int>> list_of_survey_schedules_time;
    std::vector<std::vector<size_t>> list_of_survey_schedules_vid;
    int last_surveyed_index = 0;


    // int num_actioned_fmdas = 0;

    std::vector<std::string> survey_record; // "village_id,survey_time,survey_result,survey_id"


    std::vector<Mda> triggered_mdas;

    Fmda(
        const int num_teams,
        const int num_villages,

        const int num_days_spent_in_village,
        const int num_drug_rounds,
        const int num_days_between_drug_rounds,

        const float* const* distances,
        const int delay,
        const float pcr_sensitivity,
        const float pcr_sample_size,
        const float prevalence_threshold,

        const std::string& output_prefix,

        const bool print_map

    );


    void register_start_time(int time_step, bool use_given_survey_schedule);
    void register_survey_time(int time_step_mean, int time_step_sd, int delivery_by_time_step);

    bool if_surveying_today(const int current_time_step) const;

    void step(
        const int current_time_step,
        village::VillageManager* vll_mgr_ptr
        // const common::ParasiteType* dominant_type_array
        // const std::string& output_prefix
    );

    void create_mda(
        const int current_time_step,
        village::VillageManager* vll_mgr_ptr
        // const common::ParasiteType* dominant_type_array
        // const std::string& output_prefix
    );

};

}

#endif