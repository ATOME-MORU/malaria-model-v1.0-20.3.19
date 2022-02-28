#include <iostream>
#include <vector>

#include "third_party/catch2/catch.hpp"


#include "intervention/vc.h"

TEST_CASE( "100.1: Vector Control", "[vc:basics]" ) {
	const float kEfficacy = 0.4;
	const int kNum_villages = 2;
	const float kBeta_season = 0.15;

	const int kStep_a = 1;
	const int kDuration_a = 2;

	const int kStep_b = 5;
	const int kDuration_b = 3;

	std::cout << "Intervention::Vc / Basics ...";

	REQUIRE( intervention::Vc::op_apply_efficacy_to_beta_season(kBeta_season, kEfficacy)
			 == (kBeta_season * (1-kEfficacy)) );

	intervention::Vc vc(kEfficacy, kNum_villages);
	std::vector<float> beta_season_of;
	beta_season_of.resize(kNum_villages, kBeta_season);

	vc.init_add_to_schedule(kStep_a, kDuration_a);
	vc.init_add_to_schedule(kStep_b, kDuration_b);

	float beta_season_target = 0.0;

	for (int tt = 0; tt < 10; tt++) {
		beta_season_of.clear();
		beta_season_of.resize(kNum_villages, kBeta_season);

		vc.step_reduce_beta_season(beta_season_of, tt);

		if (tt < kStep_a) {
			beta_season_target = kBeta_season;
		}

		if (tt >= kStep_a && tt < (kStep_a + kDuration_a)) {
			beta_season_target = kBeta_season*(1-kEfficacy);
		}

		if (tt >= (kStep_a + kDuration_a) && tt < kStep_b) {
			beta_season_target = kBeta_season;
		}

		if (tt >= kStep_b && tt < (kStep_b + kDuration_b)) {
			beta_season_target = kBeta_season*(1-kEfficacy);
		}

		if (tt >= (kStep_b + kDuration_b)) {
			beta_season_target = kBeta_season;
		}

		for (auto bb : beta_season_of) {
			// std::cout << "tt:" << tt
			// 			<< ", beta_season=" << bb 
			// 			<< ", target=" << beta_season_target
			// 			<< "\n";
			REQUIRE(bb == beta_season_target);
		}
	}

	std::cout << " PASS\n";

}