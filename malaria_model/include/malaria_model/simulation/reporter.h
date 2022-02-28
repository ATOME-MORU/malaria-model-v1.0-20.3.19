#ifndef REPORTER_H
#define REPORTER_H


namespace simulation{

class ReportManager {

    const bool kIf_batch_mode;

    const std::string kOutput_title;

    std::string output_prefix;
    const int sim_time_max;

    const int kGnuplot_pause_interval = 1; // seconds
    const int kOstream_flush_interval = 10; // steps

    const char kData_sep = ';';
    const char kCSV_sep = ',';

    std::vector<bool> kReport_selector;

    const int kPrevalence_daily_y_max;

    village::VillageManager* vll_mgr;
    human::HumanManager* hmn_mgr;
    human::BloodSystemManager* bld_mgr;


    const std::vector<std::string> kReport_file_names = {
        "_prevalence_daily.csv",         // 0
        "_drug_daily.csv",               // 1
        "_infection_daily.csv",          // 2
        "_infection_daily_village.csv",  // 3
        "_mosquito_daily.csv",           // 4
        "_treatment_daily.csv",          // 5
        "_mutation_daily.csv",           // 6
        "_clearance_natural_daily.csv",  // 7
        "_clearance_drug_daily.csv",     // 8
        "_village_prev_annual.csv",      // 9
        "_village_infection_annual.csv"  // 10
        // "_village_daily.csv",
        // "_dominant_type_per_village.csv",
        // "_dominant_type_per_age.csv",
        // "_immunity_per_village.csv",
        // "_age_distribution.csv",                                
    }; // order is important, needs to match kReport_selector

    std::map<std::string, std::ofstream*> filename_stream_map;

public:

    ReportManager(

        const bool if_batch_mode,

        const std::string s_title,
        const std::string s,
        const int time_max,

        const bool if_output_prevalence_daily,
        const bool if_output_drug_daily,
        const bool if_output_infection_daily,
        const bool if_output_infection_daily_village,
        const bool if_output_mosquito_daily,
        const bool if_output_treatment_daily,
        const bool if_output_mutation_daily,
        const bool if_output_clearance_natural_daily,
        const bool if_output_clearance_drug_daily,
        const bool if_output_village_prev_annual,
        const bool if_output_village_infection_annual,

        const int prevalence_daily_y_max,

        village::VillageManager* vm,
        human::HumanManager* hm,
        human::BloodSystemManager* bm
    );


    void open_ofstreams();
    void close_ofstreams();

    void end_of_day_reports(int sim_time);

    void annual_report_village_infection(
        int sim_time,
        village::VillageManager* v_mgr,
        human::BloodSystemManager* b_mgr,
        human::HumanManager* h_mgr,
        std::ofstream* os
    );

    void annual_report_village_prev(
        int sim_time,
        village::VillageManager* v_mgr,
        std::ofstream* os
    );

    void daily_report_mutation(
        const int sim_time,
        std::ofstream* os
    );

    void daily_report_clearance_natural(
        const int sim_time,
        std::ofstream* os
    );

    void daily_report_clearance_drug(
        const int sim_time,
        std::ofstream* os
    );

    void daily_report_treatment(
        const int sim_time,
        village::VillageManager* v_mgr,
        human::BloodSystemManager* b_mgr,
        human::HumanManager* h_mgr,
        std::ofstream* os
    );

    void daily_report_mosquito(
        const int sim_time,
        village::VillageManager* v_mgr,
        std::ofstream* os
    );

    void daily_report_infection(
        const int sim_time,
        village::VillageManager* v_mgr,
        std::ofstream* os
    );

    void daily_report_infection_village(
        const int sim_time,
        village::VillageManager* v_mgr,
        std::ofstream* os
    );

    void daily_report_drug(
        const int sim_time,
        human::BloodSystemManager* b_mgr,
        village::VillageManager* v_mgr,
        std::ofstream* os
    );


    void daily_report_prevalence(
        const int sim_time,
        human::BloodSystemManager* b_mgr,
        std::ofstream* os
    );
    void daily_report_prevalence_start_plotting(std::string data_file_name);



    void daily_report_village(int sim_time);
    void daily_report_village_print_headings(std::ofstream* os);

    void daily_report_on_dominant_type(
        const int sim_time,
        const common::ParasiteType* dominant_type_array,
        const int* at_village_array,
        const int* age_array,
        const int num_humans,
        const int num_villages,
        const int num_age_bins
        );
    void daily_report_on_dominant_type_headings(
        const int num_villages,
        const int num_age_bins,
        std::ofstream* os_village,
        std::ofstream* os_age
        );

    void daily_report_on_immunity(
        const int sim_time,
        const uint8_t* immunity_level_array,
        const uint8_t* cummulative_exposures_array,
        const int* at_village_array,
        const int num_systems,
        const int num_villages
        );
    void daily_report_on_immunity_headings(
        std::ofstream* os_immunity_level,
        const int num_villages
        );

    void daily_report_on_age_distribution(
        const int sim_time,
        const int* age_array,
        const int num_humans,
        const int num_age_bins
        );
    void daily_report_on_age_distribution_headings(
        std::ofstream* os_age_distribution,
        const int num_age_bins
        );
    // void end_of_day_reports_print_headings();

    // void prevalence_report() const;
};

}
#endif