#include <iostream>
#include <vector>

#include <cassert>

#include "util/randomness.h"
#include "util/statistics.h"

#include "village/village.h"

#include "intervention/mda.h"
#include "intervention/fmda.h"


namespace intervention {

Fmda::Fmda(
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

    ) :
        kNum_teams(num_teams),
        kNum_villages(num_villages),

        kNum_days_spent_in_village(num_days_spent_in_village),
        kNum_drug_rounds(num_drug_rounds),
        kNum_days_between_drug_rounds(num_days_between_drug_rounds),

        kVillage_distances(distances),
        kSurvey_data_delivery_delay(delay),
        kPcr_sensitivity(pcr_sensitivity),
        kPcr_sample_size(pcr_sample_size),
        kPrevalence_threshold(prevalence_threshold),
        kOutput_prefix(output_prefix),
        kPrint_map(print_map),

        latest_survey_results(std::vector<bool>(num_villages, false))
    {
        this->survey_record.push_back("village,at_pop,time_step,survey_prev,true_prev,survey_group,selected_last_survey,selected_next_survey");
}

void Fmda::register_start_time(int time_step, bool use_given_survey_schedule) {
    this->scheduled_start_times.push_back(time_step);
    if (use_given_survey_schedule == false) {
        this->register_survey_time(time_step - this->kSurvey_data_delivery_delay, 0, time_step);
    }
}

void Fmda::register_survey_time(int time_step_mean, int time_step_sd, int delivery_by_time_step) {
    std::vector<int> new_schedule(this->kNum_villages,0);
    int kReloop_max = this->kNum_villages * 2;
    int reloop_count = 0;
    for (int vv = 0; vv < this->kNum_villages; vv++) {
        new_schedule[vv] = util::get_rand_normal(time_step_mean, time_step_sd);
        if (delivery_by_time_step > 0) {
            assert(delivery_by_time_step > (time_step_mean+time_step_sd+this->kSurvey_data_delivery_delay) );
            if (new_schedule[vv] > (delivery_by_time_step - this->kSurvey_data_delivery_delay)) {
                vv--;
                reloop_count++;
                assert(reloop_count < kReloop_max);
            }
        }
    }
    this->list_of_survey_schedules_time.push_back(new_schedule);
    this->list_of_survey_schedules_vid.push_back(util::sort_by_value_and_return_indices_ascend(new_schedule));

    // std::vector<size_t> idx(new_schedule.size());
    // std::iota(idx.begin(), idx.end(), 0);
    // this->list_of_survey_schedules_vid.push_back(idx);

    // for (int vv = 0; vv < this->kNum_villages; vv++) {
    //     std::cout << "v" << this->list_of_survey_schedules_vid.back()[vv]
    //             << " on d" << this->list_of_survey_schedules_time.back()[this->list_of_survey_schedules_vid.back()[vv]]
    //             << ", ";
    // }

    // std::cin.ignore();

}

bool Fmda::if_surveying_today(const int current_time_step) const {
    for (size_t ss = 0; ss < this->list_of_survey_schedules_vid.size(); ss++) {
        if (!this->list_of_survey_schedules_vid[ss].empty()) {
            int vv = this->list_of_survey_schedules_vid[ss].front();
            if (this->list_of_survey_schedules_time[ss][vv] == current_time_step) {
                return true;
            }
        }
    }
    return false;
}

void Fmda::step(
        const int current_time_step,
        village::VillageManager* vll_mgr_ptr
        // const common::ParasiteType* dominant_type_array
        // const std::string& output_prefix
    ){

    // do survey
    for (size_t ss = 0; ss < this->list_of_survey_schedules_vid.size(); ss++) {
        while(!this->list_of_survey_schedules_vid[ss].empty()){        
            int vv = this->list_of_survey_schedules_vid[ss].front();
            if (this->list_of_survey_schedules_time[ss][vv] == current_time_step) {

                std::string selection_history = "";
                if (this->latest_survey_results[vv]) {
                    selection_history = "y";
                } else {
                    selection_history = "n";
                }

                float prev_by_survey = vll_mgr_ptr->get_prevalence_of_village_by_survey(
                                            vv,
                                            this->kPcr_sample_size,
                                            this->kPcr_sensitivity
                                        );
                this->latest_survey_results[vv] = (prev_by_survey > this->kPrevalence_threshold);

                if (this->latest_survey_results[vv]) {
                    selection_history += ",y";
                } else {
                    selection_history += ",n";
                }

                this->list_of_survey_schedules_vid[ss].erase(this->list_of_survey_schedules_vid[ss].begin());

                this->survey_record.push_back(
                    std::to_string(vv)
                    + "," + std::to_string(vll_mgr_ptr->get_current_population_of_village(vv))
                    + "," + std::to_string(current_time_step)
                    + "," + std::to_string(prev_by_survey)
                    + "," + std::to_string(vll_mgr_ptr->get_prevalence_of(vv))
                    + "," + std::to_string(ss)
                    + "," + selection_history
                );
                std::cout << "Fmda-survey: village " << vv << " surveyed showing prevalence of " << prev_by_survey << "\n";

            } else {
                break;
            }
        }
    }

    // create MDA
    if (!this->scheduled_start_times.empty()){

        if (current_time_step == this->scheduled_start_times.front()) {

            // std::vector<int> village_ids;

            this->create_mda(
                current_time_step,
                vll_mgr_ptr
                // dominant_type_array
                // output_prefix
            );

            this->scheduled_start_times.erase(this->scheduled_start_times.begin());
        }

    }

    // execute MDA
    if (this->triggered_mdas.size() > 0) {
        for (auto &mda : this->triggered_mdas) {
            mda.step(current_time_step, vll_mgr_ptr);

        }
    }

}

void Fmda::create_mda(
        const int current_time_step,
        village::VillageManager* vll_mgr_ptr
        // const common::ParasiteType* dominant_type_array
        // const std::string& output_prefix
    ){

    // vll_mgr_ptr->update_prevalence_count_of(dominant_type_array);

    // std::vector<int> village_ids = vll_mgr_ptr->survey_prevalence_greater_than(
    //     this->kPrevalence_threshold,
    //     this->kPcr_sample_size,
    //     this->kPcr_sensitivity
    // );

    std::vector<int> village_ids;

    std::cout << "Fmda-create: scheduling MDA including villages";
    for (int vv = 0; vv < this->kNum_villages; vv++) {
        if (this->latest_survey_results[vv]) {
            village_ids.push_back(vv);
            std::cout << " " << vv;
        }
    }
    std::cout << "\n";


    int num_villages_over_threshold = static_cast<int>(village_ids.size());

    float** distances = new float* [num_villages_over_threshold];

    for(int vv1 = 0; vv1 < num_villages_over_threshold; vv1++ ){
        
        distances[vv1] = new float[num_villages_over_threshold];
        std::fill_n(distances[vv1], num_villages_over_threshold, std::numeric_limits<float>::infinity());

        for(int vv2 = 0; vv2 < num_villages_over_threshold; vv2++ ){
            distances[vv1][vv2] = this->kVillage_distances[village_ids[vv1]][village_ids[vv2]];
        }

    }


    this->triggered_mdas.push_back(Mda(
        this->kNum_teams,
        num_villages_over_threshold,

        this->kNum_days_spent_in_village,
        this->kNum_drug_rounds,
        this->kNum_days_between_drug_rounds,
        
        distances,

        village_ids,
        // current_time_step + this->survey_data_delivery_delay
        current_time_step
    ));

    if (this->kPrint_map) {

        float* lng = new float[num_villages_over_threshold];
        float* lat = new float[num_villages_over_threshold];
        int* pop = new int[num_villages_over_threshold];

        for(int vv = 0; vv < num_villages_over_threshold; vv++ ){
            lng[vv] = vll_mgr_ptr->longitude_of[village_ids[vv]];
            lat[vv] = vll_mgr_ptr->latitude_of[village_ids[vv]];
            pop[vv] = vll_mgr_ptr->population_of[village_ids[vv]];
        }

        std::string sys_command = this->triggered_mdas.back().print_map(
            this->kOutput_prefix + std::to_string(current_time_step) + "_",
            lng,
            lat,
            pop
        ) + "&";
        if(system(sys_command.c_str())){}


        // std::cin.ignore();

        for(int vv = 0; vv < num_villages_over_threshold; vv++ ){
            delete[] distances[vv];
        }
        delete[] distances;

        delete[] lng;
        delete[] lat;
        delete[] pop;
    }

}

}