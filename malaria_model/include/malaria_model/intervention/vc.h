#ifndef VC_H
#define VC_H

namespace intervention {

class Vc {

	const float kVc_efficacy;
	std::vector<std::pair<int, int>> schedule; // <start_on, duration>

	std::vector<float> efficacy_of;

public:

	Vc(
		float vc_efficacy,
		int num_villages
	);

	void init_add_to_schedule(int start_on, int duration);

	int get_next_start_time() const;
	void set_target_on_top_villages(
		const std::vector<float>& prevalence_per_village,
		float focus_size
	);

	void step_reduce_beta_season(std::vector<float>& beta_season_of, int at_time_step);
	static float op_apply_efficacy_to_beta_season( float beta_season, float vc_efficacy);
};

}


#endif