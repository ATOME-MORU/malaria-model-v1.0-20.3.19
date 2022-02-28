#include <cassert>
#include <functional> // std::multiplies 
#include <algorithm> // std::transform

#include <iostream>

#include "intervention/vc.h"
#include "util/statistics.h"

namespace intervention {

Vc::Vc(
	float vc_efficacy,
	int num_villages
) : 
	kVc_efficacy(vc_efficacy) {

	this->efficacy_of.clear();
	this->efficacy_of.resize(num_villages, this->kVc_efficacy);
}

void Vc::init_add_to_schedule(int start_on, int duration) {
	assert(start_on >=0);
	assert(duration >0);

	if (!this->schedule.empty()) {
		assert((this->schedule.back().first + this->schedule.back().second - 1) < start_on);
	}

	this->schedule.push_back(std::make_pair(start_on, duration));
}

int Vc::get_next_start_time() const {
	if (this->schedule.empty()) {
		return -1;
	} else {
		return this->schedule[0].first;
	}
}

void Vc::set_target_on_top_villages(
		const std::vector<float>& prevalence_per_village,
		float target_size
	) {

	assert(prevalence_per_village.size() == this->efficacy_of.size());
	assert(target_size <= 1.0);
	assert(target_size >= 0.0);

	int num_villages_to_vc = static_cast<int>(prevalence_per_village.size() * target_size);
	assert(num_villages_to_vc <= static_cast<int>(this->efficacy_of.size()));
	assert(num_villages_to_vc >= 0);

	for (auto& ee : this->efficacy_of) {
		ee = 0.0;
	}

	std::vector<size_t> village_ids_in_decending_prevalence =
		util::sort_by_value_and_return_indices_descend(prevalence_per_village);

	std::cout << "set vc efficacy(" << this->kVc_efficacy
				<< ") to the following "
				<< num_villages_to_vc  << " villages:\n";
	for (int cc = 0; cc < num_villages_to_vc; cc++) {

		int village_id = static_cast<int>(village_ids_in_decending_prevalence[cc]);

		this->efficacy_of[village_id] = this->kVc_efficacy;
		std::cout << village_id << "("
					<< prevalence_per_village[village_id] << "), ";

	}
	std::cout << "\n";

}

void Vc::step_reduce_beta_season(std::vector<float>& beta_season_of, int at_time_step) {
	if (!this->schedule.empty()) {
		if ( at_time_step >= this->schedule[0].first) {
			if ( at_time_step == (this->schedule[0].first + this->schedule[0].second)) {
				this->schedule.erase(this->schedule.begin());
				std::cout << "day" << at_time_step << ":vc_finished\n";
			} else {
				assert(beta_season_of.size() == this->efficacy_of.size());
				std::transform(
					beta_season_of.begin(),
					beta_season_of.end(),
					this->efficacy_of.begin(),
					beta_season_of.begin(),
					op_apply_efficacy_to_beta_season
				);
			}
		}
	}
}

float Vc::op_apply_efficacy_to_beta_season( float beta_season, float vc_efficacy ) {
	return beta_season * (1.0 - vc_efficacy);
}

// void Vc::set_target_villages() {

// }


}