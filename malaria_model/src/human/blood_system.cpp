#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <string>

#include <fstream>
#include <sstream>


#include <random>

#include "human/human.h"
#include "human/blood_system.h"
#include "common/common.h"
#include "util/util.h"
#include "util/randomness.h"

namespace human {

BloodSystemManager::BloodSystemManager (
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
        
    ) : 
        kTreat_clean_hosts(treat_clean_hosts),
        kTreat_asymptomatic_hosts(treat_asymptomatic_hosts),
        kInit_clinical_probability(init_clinical_probability),

        kSusceptibility_default(susceptibility_default),
        kFunc_update_susceptibility_version(func_update_susceptibility_version),

        kInfectiousness_default(infectiousness_default),

        kItn_effect(itn_mosq_preference_indoor * itn_mosq_preference_human * itn_effectiveness),

        kPhi_a(phi_a),
        kPhi_c(phi_c),

        kProphylaxis_effect(prophylaxis_effect),

        kFunc_dominance_version(func_dominance_version),

        kFunc_clinical_probability_blood_to_infectious_version(
            func_clinical_probability_blood_to_infectious_version
        ),
        kFunc_clinical_probability_blood_to_infectious_version3_para1(
            func_clinical_probability_blood_to_infectious_version3_para1
        ),
        kFunc_clinical_probability_blood_to_infectious_version3_para2(
            func_clinical_probability_blood_to_infectious_version3_para2
        ),

        kMellow_rate(1.0/mellow_days),

        kStage_progression_l2b_days_mean(stage_progression_l2b_days_mean),
        kStage_progression_l2b_days_sd(stage_progression_l2b_days_sd),
        kStage_progression_b2i_days_mean(stage_progression_b2i_days_mean),
        kStage_progression_b2i_days_sd(stage_progression_b2i_days_sd),
        kStage_progression_i2r_days_mean(stage_progression_i2r_days_mean), // recovery
        kStage_progression_i2r_days_sd(stage_progression_i2r_days_sd),


        // kNatural_recovery_rate(1.0/natural_recovery_days),
        kNatural_recovery_resistance_cost(natural_recovery_resistance_cost),

        kDrug_failure_if_clinical_on_day(drug_failure_if_clinical_on_day),
        kAllow_act_to_a(allow_act_to_a),
        kPossible_to_have_drug_a_only(possible_to_have_drug_a_only),
        // kDrug_b_resistance(drug_b_resistance),
        // kDrug_ab_resistance(drug_ab_resistance),

        kMutation_prob_single_resistance(mutation_prob_single_resistance),
        kMutation_prob_single_resistance_w_act(mutation_prob_single_resistance_w_act),
        kMutation_prob_single_type_a_type_b_ratio(mutation_prob_single_type_a_type_b_ratio),
        kMutation_type_a_enabled(mutation_type_a_enabled),
        kMutation_type_b_enabled(mutation_type_b_enabled),


        // kMutation_prob_double_resistance(mutation_prob_double_resistance),

        kFunc_parasite_mutation_version(func_parasite_mutation_version),
        kFunc_parasite_mutation_fixation_after_day(func_parasite_mutation_fixation_after_day),

        kImmunity_days_until_loss_start(immunity_days_until_loss_start),
        kImmunity_level_until_loss_start(immunity_level_until_loss_start),
        kImmunity_loss_rate(1.0/immunity_loss_days),
        report_daily_new_clinical_case(init_number_of_systems, false),
        report_daily_new_asymptomatic_case(init_number_of_systems, false)
    {

    // Allocate memory spaces, per system variables

    this->susceptibility_of = new float[init_number_of_systems];
    this->infectiousness_of = new float[init_number_of_systems];

    this->is_clinical_of = new bool[init_number_of_systems];
    this->moi_of = new uint8_t[init_number_of_systems];
    this->dominant_type_of = new common::ParasiteType[init_number_of_systems];
    this->most_advanced_stage_of = new common::StageName[init_number_of_systems];

    this->immunity_days_of = new int[init_number_of_systems];
    this->immunity_level_of = new uint8_t[init_number_of_systems];
    this->cummulative_exposures_of = new uint8_t[init_number_of_systems];

    this->drug_of = new common::DrugName[init_number_of_systems];
    this->drug_given_of = new common::DrugName[init_number_of_systems];
    // this->drug_days_of = new int[init_number_of_systems];
    this->drug_given_day_of = new int[init_number_of_systems];

    this->last_recovery_day_of = new int[init_number_of_systems];

    this->num_of_clinical_pgs_of = new uint8_t[init_number_of_systems];

    // Allocate memory spaces, per pg variables, "pg" - parasite group
    int num_parasite_groups = init_number_of_systems * kNum_pgs_per_system;

    this->pg_type_of            = new common::ParasiteType[num_parasite_groups];
    this->pg_stage_of           = new common::StageName[num_parasite_groups];
    this->pg_aquired_on_day_of        = new uint16_t[num_parasite_groups];
    this->pg_num_parasites_of   = new uint32_t[num_parasite_groups]; // TODO: change type to float
    this->pg_num_gametocytes_of = new uint32_t[num_parasite_groups];
    this->pg_is_drug_target_of = new bool[num_parasite_groups];
    this->pg_has_caused_clinicalness_of = new bool[num_parasite_groups];

    this->head_pg_of = new int8_t[num_parasite_groups];
    this->tail_pg_of = new int8_t[num_parasite_groups];

    this->report_stage_progression_flags = new bool[num_parasite_groups];
 
    this->sum_num_pgs = 0;
    this->sum_num_systems = 0;

    if(kDebug_assert){
        this->debug_bites_accumulator_add.resize(init_number_of_systems,0);
        this->debug_bites_accumulator_less.resize(init_number_of_systems,0);
    }

    for (int ss = 0; ss < init_number_of_systems; ss++) {

        this->reset_system_at(ss);
        this->sum_num_pgs += kNum_pgs_per_system;
        this->sum_num_systems++;

    }

    // this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kRb)]
    //     [static_cast<int>(common::DrugName::kPiperaquine)]
    //     [static_cast<int>(common::StageName::kBlood)] = 
    //         this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
    //             [static_cast<int>(common::DrugName::kPiperaquine)]
    //             [static_cast<int>(common::StageName::kBlood)]
    //         * this->kDrug_b_resistance;

    // this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kRb)]
    //     [static_cast<int>(common::DrugName::kPiperaquine)]
    //     [static_cast<int>(common::StageName::kInfectious)] = 
    //         this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
    //             [static_cast<int>(common::DrugName::kPiperaquine)]
    //             [static_cast<int>(common::StageName::kInfectious)]
    //         * this->kDrug_b_resistance;

    // this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kRb)]
    //     [static_cast<int>(common::DrugName::kPiparte)]
    //     [static_cast<int>(common::StageName::kBlood)] = 
    //         this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
    //             [static_cast<int>(common::DrugName::kPiparte)]
    //             [static_cast<int>(common::StageName::kBlood)]
    //         * this->kDrug_ab_resistance;

    // this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kRb)]
    //     [static_cast<int>(common::DrugName::kPiparte)]
    //     [static_cast<int>(common::StageName::kInfectious)] = 
    //         this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
    //             [static_cast<int>(common::DrugName::kPiparte)]
    //             [static_cast<int>(common::StageName::kInfectious)]
    //         * this->kDrug_ab_resistance;


    this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
        [static_cast<int>(common::DrugName::kArtesunate)]
        [static_cast<int>(common::StageName::kBlood)]
    *= base_effect_multiplier;

    this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
        [static_cast<int>(common::DrugName::kArtesunate)]
        [static_cast<int>(common::StageName::kInfectious)]
    *= base_effect_multiplier;

    this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
        [static_cast<int>(common::DrugName::kPiperaquine)]
        [static_cast<int>(common::StageName::kBlood)]
    *= base_effect_multiplier;

    this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
        [static_cast<int>(common::DrugName::kPiperaquine)]
        [static_cast<int>(common::StageName::kInfectious)]
    *= base_effect_multiplier;

    this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
        [static_cast<int>(common::DrugName::kPiparte)]
        [static_cast<int>(common::StageName::kBlood)]
    *= base_effect_multiplier;

    this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
        [static_cast<int>(common::DrugName::kPiparte)]
        [static_cast<int>(common::StageName::kInfectious)]
    *= base_effect_multiplier;



    if(resistance_multipliers_enabled) {

        this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kRb)]
            [static_cast<int>(common::DrugName::kPiperaquine)]
            [static_cast<int>(common::StageName::kBlood)] = 
                this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
                    [static_cast<int>(common::DrugName::kPiperaquine)]
                    [static_cast<int>(common::StageName::kBlood)]
                * drug_b_effect_on_rb;

        this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kRb)]
            [static_cast<int>(common::DrugName::kPiperaquine)]
            [static_cast<int>(common::StageName::kInfectious)] = 
                this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
                    [static_cast<int>(common::DrugName::kPiperaquine)]
                    [static_cast<int>(common::StageName::kInfectious)]
                * drug_b_effect_on_rb;

        this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kRb)]
            [static_cast<int>(common::DrugName::kPiparte)]
            [static_cast<int>(common::StageName::kBlood)] = 
                this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
                    [static_cast<int>(common::DrugName::kPiparte)]
                    [static_cast<int>(common::StageName::kBlood)]
                * drug_ab_effect_on_rb;

        this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kRb)]
            [static_cast<int>(common::DrugName::kPiparte)]
            [static_cast<int>(common::StageName::kInfectious)] = 
                this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
                    [static_cast<int>(common::DrugName::kPiparte)]
                    [static_cast<int>(common::StageName::kInfectious)]
                * drug_ab_effect_on_rb;




        this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kRa)]
            [static_cast<int>(common::DrugName::kArtesunate)]
            [static_cast<int>(common::StageName::kBlood)] = 
                this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
                    [static_cast<int>(common::DrugName::kArtesunate)]
                    [static_cast<int>(common::StageName::kBlood)]
                * drug_a_effect_on_ra;

        this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kRa)]
            [static_cast<int>(common::DrugName::kArtesunate)]
            [static_cast<int>(common::StageName::kInfectious)] = 
                this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
                    [static_cast<int>(common::DrugName::kArtesunate)]
                    [static_cast<int>(common::StageName::kInfectious)]
                * drug_a_effect_on_ra;

        this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kRa)]
            [static_cast<int>(common::DrugName::kPiparte)]
            [static_cast<int>(common::StageName::kBlood)] = 
                this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
                    [static_cast<int>(common::DrugName::kPiparte)]
                    [static_cast<int>(common::StageName::kBlood)]
                * drug_ab_effect_on_ra;

        this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kRa)]
            [static_cast<int>(common::DrugName::kPiparte)]
            [static_cast<int>(common::StageName::kInfectious)] = 
                this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
                    [static_cast<int>(common::DrugName::kPiparte)]
                    [static_cast<int>(common::StageName::kInfectious)]
                * drug_ab_effect_on_ra;





        this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kRab)]
            [static_cast<int>(common::DrugName::kArtesunate)]
            [static_cast<int>(common::StageName::kBlood)] = 
                this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
                    [static_cast<int>(common::DrugName::kArtesunate)]
                    [static_cast<int>(common::StageName::kBlood)]
                * drug_a_effect_on_rab;

        this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kRab)]
            [static_cast<int>(common::DrugName::kArtesunate)]
            [static_cast<int>(common::StageName::kInfectious)] = 
                this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
                    [static_cast<int>(common::DrugName::kArtesunate)]
                    [static_cast<int>(common::StageName::kInfectious)]
                * drug_a_effect_on_rab;

        this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kRab)]
            [static_cast<int>(common::DrugName::kPiperaquine)]
            [static_cast<int>(common::StageName::kBlood)] = 
                this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
                    [static_cast<int>(common::DrugName::kPiperaquine)]
                    [static_cast<int>(common::StageName::kBlood)]
                * drug_b_effect_on_rab;

        this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kRab)]
            [static_cast<int>(common::DrugName::kPiperaquine)]
            [static_cast<int>(common::StageName::kInfectious)] = 
                this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
                    [static_cast<int>(common::DrugName::kPiperaquine)]
                    [static_cast<int>(common::StageName::kInfectious)]
                * drug_b_effect_on_rab;

        this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kRab)]
            [static_cast<int>(common::DrugName::kPiparte)]
            [static_cast<int>(common::StageName::kBlood)] = 
                this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
                    [static_cast<int>(common::DrugName::kPiparte)]
                    [static_cast<int>(common::StageName::kBlood)]
                * drug_ab_effect_on_rab;

        this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kRab)]
            [static_cast<int>(common::DrugName::kPiparte)]
            [static_cast<int>(common::StageName::kInfectious)] = 
                this->kPds_kill_probability[static_cast<int>(common::ParasiteType::kR0)]
                    [static_cast<int>(common::DrugName::kPiparte)]
                    [static_cast<int>(common::StageName::kInfectious)]
                * drug_ab_effect_on_rab;

    }





    // this->sum_num_systems = init_number_of_systems;

    assert((this->sum_num_pgs == init_number_of_systems*kNum_pgs_per_system) && 
            "BloodSystemManager.constructor: this->sum_num_pgs != init_number_of_systems*kNum_pgs_per_system");
    assert((this->sum_num_systems == init_number_of_systems) && 
            "BloodSystemManager.constructor: this->sum_num_systems != init_number_of_systems");

}
void BloodSystemManager::reset_system_at(int system_index){
    int lbound = system_index*kNum_pgs_per_system;

    for (int pp = lbound; pp < lbound+kNum_pgs_per_system; pp++){
        this->pg_type_of[pp]           = common::ParasiteType::kNoParasite;
        this->pg_stage_of[pp]          = common::StageName::kNotInSystem;
        this->pg_aquired_on_day_of[pp] = 0;
        this->pg_num_parasites_of[pp]  = 0;
        this->pg_num_gametocytes_of[pp] = 0;
        this->pg_is_drug_target_of[pp] = false;
        this->pg_has_caused_clinicalness_of[pp] = false;

    }

    this->head_pg_of[system_index] = -1;
    this->tail_pg_of[system_index] = -1;

    this->drug_of[system_index] = common::DrugName::kNoDrug;
    this->drug_given_of[system_index] = common::DrugName::kNoDrug;
    // this->drug_days_of[system_index] = -1;
    this->drug_given_day_of[system_index] = -1;

    this->last_recovery_day_of[system_index] = -1;

    this->num_of_clinical_pgs_of[system_index] = 0;

    this->susceptibility_of[system_index] = this->kSusceptibility_default;
    this->infectiousness_of[system_index] = this->kInfectiousness_default;

    this->is_clinical_of[system_index] = false;
    this->moi_of[system_index] = 0;
    this->dominant_type_of[system_index] = common::ParasiteType::kNoParasite;
    this->most_advanced_stage_of[system_index] = common::StageName::kNotInSystem;

    this->immunity_days_of[system_index] = -1;
    this->immunity_level_of[system_index] = 1;
    this->cummulative_exposures_of[system_index] = 0;

    this->report_stage_progression_flags[system_index] = false;

    this->report_daily_new_clinical_case[system_index] = false;
    this->report_daily_new_asymptomatic_case[system_index] = false;

    if(kDebug_assert){
        this->debug_bites_accumulator_add[system_index] = 0;
        this->debug_bites_accumulator_less[system_index] = 0;
    }

}
void BloodSystemManager::birth_and_death(bool* flags, int num_systems){
    assert((num_systems == this->sum_num_systems) && 
            "BloodSystemManager::birth_and_death(bool* flags, int num_systems): num_systems!=this->sum_num_systems" );
    for (int ss = 0; ss < num_systems; ss++) {
        if (flags[ss]) this->reset_system_at(ss);
    }
}

BloodSystemManager::~BloodSystemManager (){
    if (drug_of) delete[] drug_of;
    if (drug_given_of) delete[] drug_given_of;
    // if (drug_days_of) delete[] drug_days_of;
    if (drug_given_day_of) delete[] drug_given_day_of;

    if (last_recovery_day_of) delete[] last_recovery_day_of;

    if (num_of_clinical_pgs_of) delete[] num_of_clinical_pgs_of;

    if (susceptibility_of) delete[] susceptibility_of;
    if (infectiousness_of) delete[] infectiousness_of;

    if (is_clinical_of) delete[] is_clinical_of;
    if (moi_of) delete[] moi_of;
    if (dominant_type_of) delete[] dominant_type_of;
    if (most_advanced_stage_of) delete[] most_advanced_stage_of;

    if (immunity_days_of) delete[] immunity_days_of;
    if (immunity_level_of) delete[] immunity_level_of;
    if (cummulative_exposures_of) delete[] cummulative_exposures_of;

    if (pg_type_of            ) delete[] pg_type_of            ;
    if (pg_stage_of           ) delete[] pg_stage_of           ;
    if (pg_aquired_on_day_of        ) delete[] pg_aquired_on_day_of        ;
    if (pg_num_parasites_of   ) delete[] pg_num_parasites_of   ;
    if (pg_num_gametocytes_of ) delete[] pg_num_gametocytes_of ;
    if (pg_is_drug_target_of) delete[] pg_is_drug_target_of;
    if (pg_has_caused_clinicalness_of) delete[] pg_has_caused_clinicalness_of;

    if (head_pg_of) delete[] head_pg_of;
    if (tail_pg_of) delete[] tail_pg_of;

    if (report_stage_progression_flags) delete[] report_stage_progression_flags;

}

void BloodSystemManager::init_for_test(){
    for (int ss = 0; ss < this->sum_num_systems; ss++) {
        this->drug_of[ss] = common::DrugName::kArtesunate;
        this->drug_given_of[ss] = common::DrugName::kArtesunate;
        // this->drug_days_of[ss] = 0;
        this->drug_given_day_of[ss] = 0;

        int lbound = ss*kNum_pgs_per_system;
        for (int ii = 0; ii < 4; ii++) {
            this->pg_type_of[lbound+this->tail_pg_of[ss]+1]            = common::ParasiteType::kRa;
            this->pg_stage_of[lbound+this->tail_pg_of[ss]+1]           = common::StageName::kBlood;
            this->pg_aquired_on_day_of[lbound+this->tail_pg_of[ss]+1]        = 0;
            this->pg_num_parasites_of[lbound+this->tail_pg_of[ss]+1]   = 1000/(ii+1);
            if (ii==0) this->pg_num_parasites_of[lbound+this->tail_pg_of[ss]+1] = 100;
            if (ii==3) this->pg_num_parasites_of[lbound+this->tail_pg_of[ss]+1] = 1000;
            this->pg_num_gametocytes_of[lbound+this->tail_pg_of[ss]+1] = 0;
    
            this->tail_pg_of[ss] ++;
            if (this->head_pg_of[ss]<0) this->head_pg_of[ss] = this->tail_pg_of[ss];
            
        }

    }
}

void BloodSystemManager::init_for_drug_failure_test(
        common::ParasiteType parasite_type,
        common::DrugName drug_to_test
    ) {
    for (int ss = 0; ss < this->sum_num_systems; ss++) {
        this->drug_of[ss] = drug_to_test;
        this->drug_given_of[ss] = drug_to_test;
        this->drug_given_day_of[ss] = 1;

        int lbound = ss*kNum_pgs_per_system;
        this->pg_type_of[lbound] = parasite_type;
        this->pg_stage_of[lbound] = common::StageName::kInfectious;
        this->pg_aquired_on_day_of[lbound] = 0;
        this->pg_num_parasites_of[lbound] = 10;
        this->pg_num_gametocytes_of[lbound] = 10;
        this->pg_has_caused_clinicalness_of[lbound] = true;
        this->pg_is_drug_target_of[lbound] = true;

        this->head_pg_of[ss] = 0;
        this->tail_pg_of[ss] = 0;

        this->is_clinical_of[ss] = true;
        this->most_advanced_stage_of[ss] = common::StageName::kInfectious;
        this->moi_of[ss] = 1;
        this->num_of_clinical_pgs_of[ss] = 1;
    }

}

void BloodSystemManager::init_from_file(
        std::string input_file_name,
        char input_file_delimiter,
        human::HumanManager* human_manager,
        bool if_use_age_column
    ){

    // const int kNum_columns = 5;

    std::cout << "BlookSystemManager::init_from_file(" << input_file_name << ")" << std::endl;
    std::ifstream input_file_stream(input_file_name);
    assert(input_file_stream);

    std::string str_line;

    for (int ss = 0; ss < this->sum_num_systems; ss++) {
        if(!std::getline(input_file_stream, str_line).eof()) {

            std::string item;
            std::stringstream s_stream(str_line);

            int col = 0;

            while (std::getline(s_stream, item, input_file_delimiter)) {

                switch (col) {
                    case 0 :
                        if (if_use_age_column) {
                            human_manager->set_human_age(ss, std::stoi(item));
                        }
                        break;
                    case 1 :
                        this->moi_of[ss] = uint8_t(std::stoi(item));
                        break;
                    case 2 :
                        // #define CLINICAL 5
                        // #define ASYMPTOMATIC 6
                        // #define NOINFECTION 7
                        // if (std::stoi(item)==5){
                        //     this->is_clinical_of[ss] = true;
                        // } else {
                        //     this->is_clinical_of[ss] = false;
                        // }
                        break;
                    case 3 :
                        this->immunity_level_of[ss] = uint8_t(std::stoi(item));
                        if (this->immunity_level_of[ss] > 0){
                            this->immunity_days_of[ss] = this->kImmunity_days_until_loss_start;
                        }
                        break;
                    case 4 :
                        this->cummulative_exposures_of[ss] = uint8_t(std::stoi(item));
                        break;
                    default:
                        std::cout << "BloodSystemManager::init_from_file : un-recognised column:"
                                  << col << std::endl;
                        break;
                }

                col++;

            } // while columns

        } else {
            std::cout << "Warning: Input file has "<< ss+1 <<" lines. Expecting "
                        << this->sum_num_systems << ". Press any key to continue.\n";
            std::cin.ignore();
        }
        
    }

    input_file_stream.close();

}

void BloodSystemManager::export_to_file(
        const std::string output_file_prefix,
        const human::HumanManager* human_manager
    ) {

    assert(this->sum_num_systems == human_manager->sum_num_humans);

    std::ofstream ofs(output_file_prefix + this->kExport_file_name);
    assert(ofs);
    
    const int* age_array = human_manager->get_age_of();

    ofs << "age" << this->kData_sep_char                // 0
        << "moi_string" << this->kData_sep_char         // 1
        << "clinical" << this->kData_sep_char           // 2
        << "immunity_level" << this->kData_sep_char     // 3
        << "immunity_days" << this->kData_sep_char      // 4
        << "cummulative_exposures" << this->kData_sep_char  //5
        << "home-village" << this->kData_sep_char       // 6
        << "at-village" << this->kData_sep_char         // 7
        << "attractiveness" << this->kData_sep_char     // 8
        << "bitten_times_m2h" << this->kData_sep_char   // 9
        << "bitten_times_h2m" << this->kData_sep_char   // 10
        << "\n";

    for (int ss = 0; ss < this->sum_num_systems; ss++) {
        ofs << age_array[ss] << this->kData_sep_char;

        if (this->tail_pg_of[ss] < 0) {
            ofs << this->kMoi_clean_char;
        } else {
            int lbound = ss*kNum_pgs_per_system;
            for (int pp = lbound + this->head_pg_of[ss]; pp <= lbound + this->tail_pg_of[ss]; pp++) {
                switch (this->pg_stage_of[pp]) {
                    case common::StageName::kLiver :
                        ofs << this->kMoi_liver_char;
                        break;
                    case common::StageName::kBlood :
                        ofs << this->kMoi_blood_char;
                        break;
                    case common::StageName::kInfectious :
                        // if (this->pg_is_drug_target_of[pp]) {
                        //     ofs << this->kMoi_infectious_target_char;
                        // } else
                        if (this->pg_has_caused_clinicalness_of[pp]) {
                            ofs << this->kMoi_infectious_clinical_char;
                        } else {
                            ofs << this->kMoi_infectious_asymptomatic_char;
                        }
                        break;
                    default :
                        ofs << this->kMoi_unknown_char;
                        break;
                }
            }
        }
        ofs << this->kData_sep_char;

        // if (this->is_clinical_of[ss]) {
        //     ofs << "1";
        // } else {
        //     ofs << "0";
        // }
        ofs << this->is_clinical_of[ss];
        ofs << this->kData_sep_char;

        ofs << static_cast<int>(this->immunity_level_of[ss]) << this->kData_sep_char;
        ofs << static_cast<int>(this->immunity_days_of[ss]) << this->kData_sep_char;

        ofs << static_cast<int>(this->cummulative_exposures_of[ss]) << this->kData_sep_char;

        ofs << human_manager->get_home_village(ss) << this->kData_sep_char;
        ofs << human_manager->get_at_village(ss) << this->kData_sep_char;
        ofs << human_manager->get_attractiveness(ss) << this->kData_sep_char;
        ofs << human_manager->get_log_bitten_times_m2h(ss) << this->kData_sep_char;
        ofs << human_manager->get_log_bitten_times_h2m(ss) << this->kData_sep_char;

        ofs << "\n";
    }

    ofs.close();
}

void BloodSystemManager::import_from_file(
        const std::string file_name,
        human::HumanManager* human_manager,
        const bool if_use_age_column
    ) {

    assert(human_manager->sum_num_humans == this->sum_num_systems);

    std::cout << "B_sys: importing from " << file_name << "\n";
    std::ifstream ifs(file_name);
    assert(ifs);

    std::string str_per_row;
    std::getline(ifs, str_per_row); // skip header at first row
    
    for (int ss = 0; ss < this->sum_num_systems; ss++) {
        if(!std::getline(ifs, str_per_row).eof()) {

            std::stringstream stream_per_row(str_per_row);
            std::string str_per_cell;

            int column_index = 0;

            while (std::getline(stream_per_row, str_per_cell, this->kData_sep_char)) {
                switch (column_index++) {
                    case 0 :
                        if (if_use_age_column) {
                            human_manager->set_human_age(ss, std::stoi(str_per_cell));
                        }
                        break;
                    case 1 : {
                        for (char& moi_char : str_per_cell) {
                            // std::cout << moi_char << "of" << str_per_cell << "\n";
                            // if (static_cast<int>(this->moi_of[ss]) > kNum_pgs_per_system/2) {
                            //     continue;
                            // }
                            int new_pg_index;

                            if (moi_char == this->kMoi_clean_char) {
                                continue;
                            } else if (moi_char == this->kMoi_blood_char) {
                                new_pg_index = this->add_one_pg_with_defaults(
                                    ss,
                                    common::ParasiteType::kR0,
                                    common::StageName::kBlood,
                                    0
                                );
                                this->moi_of[ss]++;
                            } else if (moi_char == this->kMoi_infectious_clinical_char) {
                                new_pg_index = this->add_one_pg_with_defaults(
                                    ss,
                                    common::ParasiteType::kR0,
                                    common::StageName::kInfectious,
                                    0
                                );
                                this->pg_has_caused_clinicalness_of[new_pg_index] = true;
                                this->num_of_clinical_pgs_of[ss]++;
                                this->moi_of[ss]++;
                            } else if (moi_char == this->kMoi_infectious_asymptomatic_char) {
                                new_pg_index = this->add_one_pg_with_defaults(
                                    ss,
                                    common::ParasiteType::kR0,
                                    common::StageName::kInfectious,
                                    0
                                );
                                this->moi_of[ss]++;
                            } else if (moi_char == this->kMoi_liver_char) {
                                new_pg_index = this->add_one_pg_with_defaults(
                                    ss,
                                    common::ParasiteType::kR0,
                                    common::StageName::kLiver,
                                    0
                                );
                            } else if (moi_char == this->kMoi_unknown_char) {

                            } else {
                                std::cout << "BloodSystemManager::import_from_file : unrecognised moi_char("
                                            << moi_char << ") at row" << ss << "+1\n";
                            }

                            // switch (moi_char) {
                            //     case this->kMoi_blood_char :
                            //         new_pg_index = this->add_one_pg_with_defaults(
                            //             ss,
                            //             common::ParasiteType::kR0,
                            //             common::StageName::kBlood,
                            //             0
                            //         );
                            //         this->moi_of[ss]++;
                            //         break;
                            //     case this->kMoi_infectious_clinical_char :
                            //         new_pg_index = this->add_one_pg_with_defaults(
                            //             ss,
                            //             common::ParasiteType::kR0,
                            //             common::StageName::kInfectious,
                            //             0
                            //         );
                            //         this->pg_has_caused_clinicalness_of[new_pg_index] = true;
                            //         this->num_of_clinical_pgs_of[ss]++;
                            //         this->moi_of[ss]++;
                            //         break;
                            //     case this->kMoi_infectious_asymptomatic_char : 
                            //         new_pg_index = this->add_one_pg_with_defaults(
                            //             ss,
                            //             common::ParasiteType::kR0,
                            //             common::StageName::kInfectious,
                            //             0
                            //         );
                            //         this->moi_of[ss]++;
                            //         break;
                            //     case this->kMoi_liver_char :
                            //         break;
                            //     case this->kMoi_clean_char :
                            //         break;
                            //     case this->kMoi_unknown_char :
                            //         break;
                            //     default:
                            //         std::cout << "BloodSystemManager::import_from_file : unrecognised moi_char("
                            //                     << moi_char << ") at row" << ss << "+1\n";
                            // }
                        }
                        break;
                    }
                    case 2 : 
                        this->is_clinical_of[ss] = (str_per_cell == "1");
                        break;
                    case 3 :
                        this->immunity_level_of[ss] = static_cast<uint8_t>(std::stoi(str_per_cell));
                        if (this->immunity_level_of[ss] == 0) {
                            this->immunity_level_of[ss] =1;
                        }
                        break;
                    case 4 :
                        this->immunity_days_of[ss] = std::stoi(str_per_cell);
                        break;
                    case 5 :
                        this->cummulative_exposures_of[ss] = static_cast<uint8_t>(std::stoi(str_per_cell));
                        break;
                    case 6 :
                        // ignoring the home village column
                        break;
                    case 7 :
                        // ignoring the at village column
                        break;
                    case 8 :
                        human_manager->set_human_attractiveness(ss,std::stof(str_per_cell));
                        break;
                    case 9 :
                        // ignore get_log_bitten_times_m2h
                        break;
                    case 10 :
                        // ignore get_log_bitten_times_h2m
                        break;
                    default:
                        std::cout << "BloodSystemManager::import_from_file : unrecognised column ("
                                    << column_index << ") at row " << ss << "+2\n";
                }
            }
        } else {
            std::cout << "BloodSystemManager::import_from_file : "
                        << this->sum_num_systems << " hosts in the system, only "
                        << ss << "-1 given in " << file_name;
            std::cin.ignore();
        }
    }
    this->parasite_dominance();
    // for (const auto& pt: util::Enum_iterator<common::ParasiteType>()) {
    //     std::cout << this->sum_prevalence_count[(int)pt] << "\n";
    // }
    // std::cin.ignore();
}

int BloodSystemManager::add_one_pg_with_defaults(
        const int sys_index,
        const common::ParasiteType pg_type,
        const common::StageName pg_stage,
        const uint16_t pg_aquired_on_day
        // const uint32_t pg_num_parasites,
        // const uint32_t pg_num_gametocytes,
        // const bool pg_is_drug_target,
        // const bool pg_has_caused_clinicalness
    ) {

    if (this->tail_pg_of[sys_index]==(kNum_pgs_per_system-1) && this->head_pg_of[sys_index]==0) {
        this->sum_num_pg_q_is_full_incidents++;
        return -1;
    }

    int lbound = sys_index*kNum_pgs_per_system;

    if (this->tail_pg_of[sys_index]==(kNum_pgs_per_system-1)) {
        // no space behind tail, push the queue to the front of the block
        // a faster implementation would be to copy from tail to 0 until (it meets with head or empty)
        // which is less readable

        int8_t target_space_offset = 0;
        for (int pp = lbound+this->head_pg_of[sys_index]; pp <= lbound+this->tail_pg_of[sys_index]; pp++){

            this->pg_type_of[lbound + target_space_offset] = this->pg_type_of[pp];
            this->pg_stage_of[lbound + target_space_offset] = this->pg_stage_of[pp];
            this->pg_aquired_on_day_of[lbound + target_space_offset] = this->pg_aquired_on_day_of[pp];
            this->pg_num_parasites_of[lbound + target_space_offset] = this->pg_num_parasites_of[pp];
            this->pg_num_gametocytes_of[lbound + target_space_offset] = this->pg_num_gametocytes_of[pp];
            this->pg_is_drug_target_of[lbound + target_space_offset] = this->pg_is_drug_target_of[pp];
            this->pg_has_caused_clinicalness_of[lbound + target_space_offset] = this->pg_has_caused_clinicalness_of[pp];

            target_space_offset++;
        }

        this->head_pg_of[sys_index] = 0;
        this->tail_pg_of[sys_index] = target_space_offset-1;
    }
    
    // now we can simply add the new pg to tail + 1
    int new_pg_index = lbound+this->tail_pg_of[sys_index]+1;

    this->pg_type_of[new_pg_index] = pg_type;
    // this->pg_stage_of[new_pg_index] = common::StageName::kLiver;
    this->pg_stage_of[new_pg_index] = pg_stage;
    this->pg_aquired_on_day_of[new_pg_index] = pg_aquired_on_day;
    this->pg_num_parasites_of[new_pg_index] = 100;
    this->pg_num_gametocytes_of[new_pg_index] = 0;
    this->pg_is_drug_target_of[new_pg_index] = false;
    this->pg_has_caused_clinicalness_of[new_pg_index] = false;

    this->tail_pg_of[sys_index]++;
    if (this->head_pg_of[sys_index] < 0) {
        this->head_pg_of[sys_index] = 0;
    }

    if(this->kDebug_assert){
        this->debug_bites_accumulator_add[sys_index]++;
    }

    return new_pg_index;

}

void BloodSystemManager::update_susceptibility(
        const bool* if_itn,
        int if_itn_length,
        const int* age_of
    ){

    assert(if_itn_length == this->sum_num_systems);

    switch(this->kFunc_update_susceptibility_version) {
        case 0 :
            this->update_susceptibility_kernel(
                if_itn,
                this->is_clinical_of,
                this->drug_of,
                this->sum_num_systems,
                this->susceptibility_of,
                this->infectiousness_of,
                this->kSusceptibility_default,
                this->kInfectiousness_default,
                this->kPhi_a,
                this->kPhi_c,
                this->kProphylaxis_effect,
                this->kItn_effect
            );
            break;

        case 1 :
            this->update_susceptibility_kernel_1(
                age_of,
                if_itn,
                this->is_clinical_of,
                this->drug_of,
                this->sum_num_systems,
                this->susceptibility_of,
                this->infectiousness_of,
                // this->kSusceptibility_default,
                this->kInfectiousness_default,
                this->kPhi_a,
                this->kPhi_c,
                this->kProphylaxis_effect,
                this->kItn_effect
            );

            break;

        default:
            std::cout << "Bsys: Unknown kFunc_update_susceptibility_version" << std::endl;
            break;    
    }

}

void BloodSystemManager::update_susceptibility_kernel_1(
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
    ) {

    for (int ss = 0; ss < array_size; ss++){


        // susceptibility_array[ss] = default_susceptibility;
        infectiousness_array[ss] = default_infectiousness;

        //person->susceptibility = (1-(0.99*exp(-0.14*person->age));
        susceptibility_array[ss] = 1-(0.99*exp(-0.14*age_array[ss]));

        if (if_itn_array[ss]) {
            // itn_counter++;
            if (if_clinical_array[ss]){
                susceptibility_array[ss] *= (1.0-itn_effect);
                infectiousness_array[ss] *= (1.0-itn_effect) * phi_c;
            } else {
                susceptibility_array[ss] *= (1.0-itn_effect);
                infectiousness_array[ss] *= (1.0-itn_effect) * phi_a;
            }
        } else {
            if (if_clinical_array[ss]){
                // susceptibility_array[ss] *= (1.0-kItn_effect);
                infectiousness_array[ss] *= phi_c;
            } else {
                // susceptibility_array[ss] *= (1.0-kItn_effect);
                infectiousness_array[ss] *= phi_a;
            }
        }

        // prophylaxis from artesunate uptake
        if ( drug_array[ss] == common::DrugName::kArtesunate
             || drug_array[ss] == common::DrugName::kPiparte) {

            susceptibility_array[ss] *= prophylaxis_effect;
        }

    }

}

void BloodSystemManager::update_susceptibility_kernel(
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
    ) {

    // float susceptibility_sum = 0.0;
    // int itn_counter = 0;

    for (int ss = 0; ss < array_size; ss++){


        susceptibility_array[ss] = default_susceptibility;
        infectiousness_array[ss] = default_infectiousness;

        if (if_itn_array[ss]) {
            // itn_counter++;
            if (if_clinical_array[ss]){
                susceptibility_array[ss] *= (1.0-itn_effect);
                infectiousness_array[ss] *= (1.0-itn_effect) * phi_c;
            } else {
                susceptibility_array[ss] *= (1.0-itn_effect);
                infectiousness_array[ss] *= (1.0-itn_effect) * phi_a;
            }
        } else {
            if (if_clinical_array[ss]){
                // susceptibility_array[ss] *= (1.0-kItn_effect);
                infectiousness_array[ss] *= phi_c;
            } else {
                // susceptibility_array[ss] *= (1.0-kItn_effect);
                infectiousness_array[ss] *= phi_a;
            }
        }

        // prophylaxis from artesunate uptake
        if ( drug_array[ss] == common::DrugName::kArtesunate
             || drug_array[ss] == common::DrugName::kPiparte) {

            susceptibility_array[ss] *= prophylaxis_effect;
        }

        // susceptibility_sum += susceptibility_array[ss];
    }

    // std::cout << "itn_counter="<< itn_counter << "\n";
    // std::cout << "itn_effect="<< itn_effect << "\n";
    // std::cout << "susceptibility_sum="<< susceptibility_sum << "\n";

}

void BloodSystemManager::clinical_resolution(){
    this->clinical_resolution_kernel(
        this->is_clinical_of,
        this->sum_num_systems,
        this->kMellow_rate
        );
}
void BloodSystemManager::clinical_resolution_kernel(
        bool* is_clinical_array,
        int array_size,
        float mellow_rate
    ){

    float rnd_nmb = 0.0;
    for (int ss = 0; ss < array_size; ss++){
        if (is_clinical_array[ss]){
            util::set_random_numbers_uniform(&rnd_nmb, 1);
            if ( rnd_nmb <= mellow_rate ){
                is_clinical_array[ss] = false;
            }
        }
    }

}


void BloodSystemManager::parasite_intake(
        common::ParasiteType* intake_list,
        const common::StageName begin_at_stage,
        const int time_step,
        const bool if_intake_for_init // If the intake action is at the beginning of the simulation
    ) {
    for (int ss = 0; ss < this->sum_num_systems; ss++) {
        if (intake_list[ss] != common::ParasiteType::kNoParasite) {
            // welcome new parasite group

            if(this->kDebug_print) {
                std::cout << "ss: start intake" << ss
                            << ", this->head_pg_of[ss]:" << int(this->head_pg_of[ss] )
                            << ", this->tail_pg_of[ss]:" << int(this->tail_pg_of[ss] )
                            << "debug_bites_accumulator_add[ss]" << debug_bites_accumulator_add[ss]
                            << ", debug_bites_accumulator_less[ss]" << debug_bites_accumulator_less[ss] << "\n";

                if ((this->tail_pg_of[ss]==(kNum_pgs_per_system-1) && this->head_pg_of[ss]==0)) {
                    std::cout << "ss: " << ss << ", moi: " << int(this->moi_of[ss])
                                << ", this->head_pg_of[ss]:" << int(this->head_pg_of[ss] )
                                << ", this->tail_pg_of[ss]:" << int(this->tail_pg_of[ss] )
                                << ", this->cummulative_exposures_of[ss]" << int( this->cummulative_exposures_of[ss])
                                << std::endl;
                    std::cout << "debug_bites_accumulator_add[ss]:" << debug_bites_accumulator_add[ss];
                    std::cout << ", debug_bites_accumulator_less[ss]:" << debug_bites_accumulator_less[ss] << std::endl;
                }
            }

            if(this->kDebug_assert){
                assert(this->tail_pg_of[ss] < int8_t(kNum_pgs_per_system));
                assert(this->tail_pg_of[ss] > -2);
                assert(this->head_pg_of[ss] < int8_t(kNum_pgs_per_system));
                assert(this->head_pg_of[ss] > -2);
                // assert(this->head_pg_of[ss] < 2);
            }

            // assert (!(this->tail_pg_of[ss]==(kNum_pgs_per_system-1) && this->head_pg_of[ss]==0) &&
            //         "one of the pg queues is full");

            this->sum_num_pg_intakes++;
            if (this->tail_pg_of[ss]==(kNum_pgs_per_system-1) && this->head_pg_of[ss]==0) {
                this->sum_num_pg_q_is_full_incidents++;
                continue;
            }

            int lbound = ss*kNum_pgs_per_system;
            if (this->tail_pg_of[ss]==(kNum_pgs_per_system-1)) {
                // no space behind tail, push the queue to the front of the block
                // a faster implementation would be to copy from tail to 0 until (it meets with head or empty)
                // which is less readable

                int8_t target_space_offset = 0;
                for (int pp = lbound+this->head_pg_of[ss]; pp <= lbound+this->tail_pg_of[ss]; pp++){

                    this->pg_type_of[lbound + target_space_offset] = this->pg_type_of[pp];
                    this->pg_stage_of[lbound + target_space_offset] = this->pg_stage_of[pp];
                    this->pg_aquired_on_day_of[lbound + target_space_offset] = this->pg_aquired_on_day_of[pp];
                    this->pg_num_parasites_of[lbound + target_space_offset] = this->pg_num_parasites_of[pp];
                    this->pg_num_gametocytes_of[lbound + target_space_offset] = this->pg_num_gametocytes_of[pp];
                    this->pg_is_drug_target_of[lbound + target_space_offset] = this->pg_is_drug_target_of[pp];
                    this->pg_has_caused_clinicalness_of[lbound + target_space_offset] = this->pg_has_caused_clinicalness_of[pp];

                    target_space_offset++;
                }
                if(this->kDebug_print) {
                    std::cout << "target_space_offset" << int(target_space_offset) << "\n";
                }

                this->head_pg_of[ss] = 0;
                this->tail_pg_of[ss] = target_space_offset-1;
            }
            
            // now we can simply add pg to tail + 1
            int tail = lbound + this->tail_pg_of[ss];
            this->pg_type_of[tail+1] = intake_list[ss];
            // this->pg_stage_of[tail+1] = common::StageName::kLiver;
            this->pg_stage_of[tail+1] = begin_at_stage;
            this->pg_aquired_on_day_of[tail+1] = time_step;
            this->pg_num_parasites_of[tail+1] = 100;
            this->pg_num_gametocytes_of[tail+1] = 0;
            this->pg_is_drug_target_of[tail+1] = false;
            this->pg_has_caused_clinicalness_of[tail+1] = false;


            // post intake status update

            if (this->pg_stage_of[tail+1] != common::StageName::kLiver){
                this->moi_of[ss]++;
            }


            if (this->pg_stage_of[tail+1] == common::StageName::kInfectious){

                float random_number = 0.0;
                util::set_random_numbers_uniform(&random_number, 1);
                bool to_cause_a_clinical_case = false;

                if (if_intake_for_init) {
                    to_cause_a_clinical_case = (random_number < this->kInit_clinical_probability);
                    // this->is_clinical_of[ss] = (random_number < 0.02); //G_initclin = 0.02                        
                } else {
                    to_cause_a_clinical_case = 
                        (random_number < this->func_clinical_probability_blood_to_infectious(
                                            this->immunity_level_of[ss],
                                            this->cummulative_exposures_of[ss],
                                            this->moi_of[ss]
                                        )
                        );
                }

                if (to_cause_a_clinical_case) {
                    this->pg_has_caused_clinicalness_of[tail+1] = true;
                    this->num_of_clinical_pgs_of[ss]++;
                    this->daily_record_of_new_cases_clinical[static_cast<int>(this->pg_type_of[tail+1])]++;
                } else {
                    this->daily_record_of_new_cases_asymptomatic[static_cast<int>(this->pg_type_of[tail+1])]++;
                }

                if (this->is_clinical_of[ss] == false) {
                    this->is_clinical_of[ss] =  to_cause_a_clinical_case;
                }
            }

            this->tail_pg_of[ss]++;
            if (this->head_pg_of[ss] < 0) {
                this->head_pg_of[ss] = 0;
            }
            
            if(this->kDebug_assert){
                this->debug_bites_accumulator_add[ss]++;
            }

            if(this->kDebug_print) {
                if((debug_bites_accumulator_add[ss]+debug_bites_accumulator_less[ss])!=(this->tail_pg_of[ss]-this->head_pg_of[ss]+1)){
                    std::cout << "ss: " << ss
                                << ", debug_bites_accumulator_add[ss]" << debug_bites_accumulator_add[ss]
                                << ", debug_bites_accumulator_less[ss]" << debug_bites_accumulator_less[ss]
                                << ", this->tail_pg_of[ss]" << int(this->tail_pg_of[ss])
                                << ", this->head_pg_of[ss]" << int(this->head_pg_of[ss]) << std::endl;
                }
            }

            if(this->kDebug_assert){
                assert ((debug_bites_accumulator_add[ss]+debug_bites_accumulator_less[ss])==(this->tail_pg_of[ss]-this->head_pg_of[ss]+1));

                assert(this->tail_pg_of[ss] < kNum_pgs_per_system);
                assert(this->head_pg_of[ss] < kNum_pgs_per_system);
            }




            // std::cout << "added to system " << ss << ", updated tail: " << (int) this->tail_pg_of[ss] << std::endl;
        }


    
        intake_list[ss] = common::ParasiteType::kNoParasite;

    }
}

void BloodSystemManager::parasite_injection(
        const common::ParasiteType parasite_type_to_inject,
        const int num_cases
    ) {

    std::vector<bool> injected_already(this->sum_num_systems, false);
    float rnd_nmb = 0.0;

    for (int ii = 0; ii < num_cases; ii++) {

        util::set_random_numbers_uniform(&rnd_nmb, 1);
        int system_index = rnd_nmb * this->sum_num_systems;

        // if (injected_already[system_index]) {
        // if (injected_already[system_index] || this->tail_pg_of[system_index] >= 0 ) {
        // if (injected_already[system_index] || this->tail_pg_of[system_index] < 0 ) {
        if (injected_already[system_index] || this->most_advanced_stage_of[system_index] != common::StageName::kInfectious ) {
            ii--;
        } else {
            int lbound = system_index * kNum_pgs_per_system;
            for (int pp = lbound + this->head_pg_of[system_index]; pp <= lbound + this->tail_pg_of[system_index]; pp++){
                // if (this->pg_stage_of[pp] == common::StageName::kInfectious) {
                    this->pg_type_of[pp] = parasite_type_to_inject;
                // }
            }
            this->dominant_type_of[system_index] = parasite_type_to_inject;
            injected_already[system_index] = true;
        }

    }

}

void BloodSystemManager::drug_intake(
        const common::DrugName* drug_list,
        const int time_step
    ){
    // int mda_drug_counter = 0;
    for (auto& ii : this->daily_record_of_drug_intake) {
        ii = 0;
    }

    for (int ss = 0; ss < this->sum_num_systems; ss++) {
        if (drug_list[ss] != common::DrugName::kNoDrug){
            // this->drug_days_of[ss] = 0;
            this->drug_of[ss] = drug_list[ss];
            this->drug_given_of[ss] = drug_list[ss];
            this->drug_given_day_of[ss] = time_step;

            this->daily_record_of_drug_intake[static_cast<int>(drug_list[ss])]++;

            if ( (this->tail_pg_of[ss] < 0) && this->kTreat_clean_hosts) {
                continue;
            }

            assert(this->tail_pg_of[ss] >= 0);

            int lbound = ss * kNum_pgs_per_system;
            int latest_infection = -1;
            int latest_infection_p_index = -1;
            for (int pp = lbound+this->head_pg_of[ss]; pp <= lbound+ this->tail_pg_of[ss]; pp++){
                // if (this->pg_stage_of[pp] == common::StageName::kInfectious && this->pg_aquired_on_day_of[pp] > latest_infection) {
                // if (this->pg_stage_of[pp] == common::StageName::kInfectious && this->pg_has_caused_clinicalness_of[pp]) {
                if (this->kTreat_asymptomatic_hosts || this->pg_has_caused_clinicalness_of[pp]) {
                    if (this->pg_aquired_on_day_of[pp] > latest_infection) {
                        latest_infection_p_index = pp;
                        latest_infection = this->pg_aquired_on_day_of[pp];
                    }
                }
            }

            // this->print_one(ss);
            assert (latest_infection_p_index >= 0);
            
            this->pg_is_drug_target_of[latest_infection_p_index] = true;

        }
    }
}

void BloodSystemManager::drug_loss(
        float* random_numbers,
        int random_numbers_array_size
    ){

    assert(random_numbers_array_size == this->sum_num_systems);

     // G_xa0 = (1.0/7.0);     // lose ART
     // G_xai = 1.0/3.0;       // lose DHA in ACT
     // G_xb = 1.0/30.0;
     // G_xab = (1.0/30.0);    // lose PIP

     // G_xprim = (1.0/2.0);   // lose primaquine
     // G_xala = 1.0/3.0;      // lose ART in AL
     // G_xall = 1.0/15.0;     // lose lumefantrine in AL
     // G_xlum = 1.0/15.0;     // lose LUM


// A4, DHA-PIP
    const float kDrug_loss_Artesunate = 1.0/7.0; // G_xa0
    // const float kDrug_loss_Piperaquine = 1.0/30.0; // G_xab
    //// const float kDrug_loss_Lumefantrine = 1.0/15.0; // G_xlum
    const float kDrug_loss_Lumefantrine = 1.0/14.0; // mmc_wp2_rev17
    const float kDrug_loss_AQ = 1.0/25.0; // mmc_wp2_rev17

    // const float kDrug_loss_Piparte_to_Piperaquine = 1.0/3.0; // G_xai // mda_coverage
    // const float kDrug_loss_Piparte_to_Piperaquine = 1.0/4.0; // G_xai // mmc_wp2_rev5, 5562e42
    const float kDrug_loss_Piparte_to_Artesunate = 1.0/30.0; // G_xab

    // const float kDrug_loss_AL_to_Lumefantrine = 1.0/3.0; // G_xala
    const float kDrug_loss_AL_to_Lumefantrine = 1.0/4.0; // mmc_wp2_rev17
    const float kDrug_loss_AL_to_Artesunate = 1.0/15.0; // G_xall

    const float kDrug_loss_ASAQ_to_AQ = 1.0/4.0; // mmc_wp2_rev17
    const float kDrug_loss_ASAQ_to_Artesunate = 1.0/100.0; // mmc_wp2_rev17, placeholder only

    const float kDrug_loss_SP = 1.0/7.0;
    const float kDrug_loss_Chloroquine = 1.0/14.0;

// A5, ASAQ
    // const float kDrug_loss_Piperaquine = kDrug_loss_AQ;
    // const float kDrug_loss_Piparte_to_Piperaquine = kDrug_loss_ASAQ_to_AQ;

// // A6, AL
    const float kDrug_loss_Piperaquine = kDrug_loss_Lumefantrine;
    const float kDrug_loss_Piparte_to_Piperaquine = kDrug_loss_AL_to_Lumefantrine;
    

    bool loss_flag = false;
    common::DrugName drug_to_change_to = common::DrugName::kNoDrug;

    for(const auto& pt: util::Enum_iterator<common::ParasiteType>()) {
        this->sum_tf_by_drug_loss_count[static_cast<int>(pt)] = 0;
    }

    for (int ss = 0; ss < this->sum_num_systems; ss++) {
        switch (this->drug_of[ss]) {
            case common::DrugName::kArtesunate :
                assert(this->kPossible_to_have_drug_a_only);
                if (random_numbers[ss] <= kDrug_loss_Artesunate) {
                    loss_flag = true;
                    drug_to_change_to = common::DrugName::kNoDrug;
                }
                break;
            case common::DrugName::kPiperaquine :
                if (random_numbers[ss] <= kDrug_loss_Piperaquine) {
                    loss_flag = true;
                    drug_to_change_to = common::DrugName::kNoDrug;
                }
                break;
            case common::DrugName::kPiparte :
                if (random_numbers[ss] <= kDrug_loss_Piparte_to_Piperaquine) {
                    loss_flag = true;
                    drug_to_change_to = common::DrugName::kPiperaquine;
                //     (void) kDrug_loss_Piparte_to_Artesunate;
                // }
                } else {
                    if (this->kAllow_act_to_a) {
                        float rnd_nmb_pa_to_a = 0.0;
                        util::set_random_numbers_uniform(&rnd_nmb_pa_to_a, 1);
                        if (rnd_nmb_pa_to_a <= kDrug_loss_Piparte_to_Artesunate) {
                            loss_flag = true;
                            drug_to_change_to = common::DrugName::kArtesunate;
                        }
                    }
                }
                break;
            case common::DrugName::kAQ :
                if (random_numbers[ss] <= kDrug_loss_AQ) {
                    loss_flag = true;
                    drug_to_change_to = common::DrugName::kNoDrug;
                }
                break;
            case common::DrugName::kASAQ :
                if (random_numbers[ss] <= kDrug_loss_ASAQ_to_AQ) {
                    loss_flag = true;
                    drug_to_change_to = common::DrugName::kAQ;
                //     (void) kDrug_loss_Piparte_to_Artesunate;
                // }
                } else {
                    if (this->kAllow_act_to_a) {
                        float rnd_nmb_pa_to_a = 0.0;
                        util::set_random_numbers_uniform(&rnd_nmb_pa_to_a, 1);
                        if (rnd_nmb_pa_to_a <= kDrug_loss_ASAQ_to_Artesunate) {
                            loss_flag = true;
                            drug_to_change_to = common::DrugName::kArtesunate;
                        }
                    }
                }
                break;
            case common::DrugName::kAL :
                if (random_numbers[ss] <= kDrug_loss_AL_to_Lumefantrine) {
                    loss_flag = true;
                    drug_to_change_to = common::DrugName::kLumefantrine;
                } else {
                    if (this->kAllow_act_to_a) {
                        float rnd_nmb_al_to_a = 0.0;
                        util::set_random_numbers_uniform(&rnd_nmb_al_to_a, 1);
                        if (rnd_nmb_al_to_a <= kDrug_loss_AL_to_Artesunate) {
                            loss_flag = true;
                            drug_to_change_to = common::DrugName::kArtesunate;
                        }
                    }
                }
                break;
            case common::DrugName::kLumefantrine :
                if (random_numbers[ss] <= kDrug_loss_Lumefantrine) {
                    loss_flag = true;
                    drug_to_change_to = common::DrugName::kNoDrug;
                }
                break;
            case common::DrugName::kSP :
                if (random_numbers[ss] <= kDrug_loss_SP) {
                    loss_flag = true;
                    drug_to_change_to = common::DrugName::kNoDrug;
                }
                break;
            case common::DrugName::kChloroquine :
                if (random_numbers[ss] <= kDrug_loss_Chloroquine) {
                    loss_flag = true;
                    drug_to_change_to = common::DrugName::kNoDrug;
                }
                break;
            case common::DrugName::kNoDrug:
                break;
        }

        // TODO
        // if (deme->prim == YES ){
        //     randomEvent=ran2();
        //     if (randomEvent<=G_xprim){
        //         deme->prim = NO;
        //         G_groups->SUM_PRIM = G_groups->SUM_PRIM - 1;
        //         G_villdata[ind]->SUM_PRIM = G_villdata[ind]->SUM_PRIM - 1;
        //     }
        // }

        if (loss_flag) {

            // std::cout << "Drug loss: from" << this->drug_of[ss] << " to " << drug_to_change_to << "\n";
            // std::cin.ignore();

            this->drug_of[ss] = drug_to_change_to;
            
            if (drug_to_change_to == common::DrugName::kNoDrug && this->tail_pg_of[ss] >= 0) {
                //TODO: record drug failure

                this->sum_tf_by_drug_loss_count[static_cast<int>(this->dominant_type_of[ss])] += 1;

            }

            
        }
            
        loss_flag = false;
        drug_to_change_to = common::DrugName::kNoDrug;

    }

}

// void BloodSystemManager::drug_effect(){
//     this->drug_effect_density_based();
// }

// void BloodSystemManager::drug_effect_probability_based(){

// }

// bool BloodSystemManager::drug_on_parasite_density_based(){

// }

void BloodSystemManager::drug_effect(
        const int time_step
    ){

    int drug_recovery_counter = 0;

    for (auto& ii : this->log_new_clearance_drug) {
        ii = 0;
    }


    const bool kIf_probability_based = true;
    // common::assign_kill_probabilities();

    float* random_number_drug_effect = new float[1];
    
    // std::cout << "begin with:" << this->pg_num_parasites_of[0] << ","
    //                            << this->pg_num_parasites_of[1] << ","
    //                            << this->pg_num_parasites_of[2] << ","
    //                            << this->pg_num_parasites_of[3] << ","
    //                            << this->pg_num_parasites_of[4] << ","
    //                            << this->pg_num_parasites_of[5] << ","
    //                            << this->pg_num_parasites_of[6] << ","
    //                            << this->pg_num_parasites_of[7] << ";"
    //                            << (int)this->head_pg_of[0] << "," << (int)this->tail_pg_of[0] << std::endl;
    for (int ss = 0; ss < this->sum_num_systems; ss++) {

        //TODO: check if no drug => continue
        if (this->drug_of[ss] == common::DrugName::kNoDrug) {
            continue;
        }
        
        int lbound = ss*kNum_pgs_per_system;
        // float drug_full_life = common::kDrug_half_life_day[(int)this->drug_of[ss]]*2;



        if (this->head_pg_of[ss] < 0){
            // no pg in this system
            // if (this->immunity_days_of[ss] >= 0) this->immunity_days_of[ss]++; // move this to immunity loss

        } else {

            for (int pp = lbound + this->head_pg_of[ss]; pp <= lbound+this->tail_pg_of[ss]; pp++){


                // apply drug

                bool if_cleared = false;

                uint32_t kill_size = 0;
                // uint32_t kill_size =
                //         common::kPds_kill_size[(int)this->pg_type_of[pp]]
                //                               [(int)this->drug_of[ss]]
                //                               [(int)this->pg_stage_of[pp]]
                //         * (
                //             (drug_full_life - this->drug_days_of[ss]) > 0 ?
                //             (drug_full_life - this->drug_days_of[ss]) : 0
                //           )
                //         /drug_full_life;


                if (kIf_probability_based) {
                    util::set_random_numbers_uniform(random_number_drug_effect, 1);
                    // if (random_number_drug_effect[0]
                    //     < common::kPds_kill_probability[(int)this->pg_type_of[pp]]
                    //                                    [(int)this->drug_of[ss]]
                    //                                    [(int)this->pg_stage_of[pp]]
                    //     ) {
                    //     if_cleared = true;
                    // }
                    // std::cout << "this->pg_type_of[pp]=" << this->pg_type_of[pp]
                    //             << ", this->drug_of[ss]=" << this->drug_of[ss]
                    //             << ", this->pg_stage_of[pp]=" << this->pg_stage_of[pp]
                    //             << ", effect:" << this->kPds_kill_probability[(int)this->pg_type_of[pp]]
                    //                                    [(int)this->drug_of[ss]]
                    //                                    [(int)this->pg_stage_of[pp]]
                    //             << "\n";
                    // if (this->pg_type_of[pp] == common::ParasiteType::kRb){
                    //     std::cin.ignore();
                    // }


                    if (random_number_drug_effect[0]
                        < this->kPds_kill_probability[(int)this->pg_type_of[pp]]
                                                       [(int)this->drug_of[ss]]
                                                       [(int)this->pg_stage_of[pp]]
                        ) {
                        // if (this->pg_type_of[pp] != common::ParasiteType::kR0) {
                        //     std::cout << "drug:" << this->drug_of[ss] << "\n";
                        //     std::cout << "kPds:" << this->kPds_kill_probability[(int)this->pg_type_of[pp]]
                        //                                    [(int)this->drug_of[ss]]
                        //                                    [(int)this->pg_stage_of[pp]] << "\n";
                        //     std::cout << "type:" << this->pg_type_of[pp] << "\n";
                        //     std::cout << "clinical:" << this->is_clinical_of[ss] << "\n";

                        //     // std::cin.ignore();
                            
                        // }

                        if_cleared = true;
                    }

                    // if ( this->pg_stage_of[pp] == common::StageName::kInfectious
                    //     && this->pg_type_of[pp] == common::ParasiteType::kRb
                    //     && this->drug_of[ss] == common::DrugName::kPiparte) {
                    //     this->num_pip_treatments++;
                    //     if (if_cleared) {
                    //         this->num_pip_treatments_success++;
                    //     }
                    // }
                    // std::cout << "random_number_drug_effect[0]: " << random_number_drug_effect[0];
                    // std::cout << ", kPds_kill_probability:" << common::kPds_kill_probability[(int)this->pg_type_of[pp]]
                    //                                    [(int)this->drug_of[ss]]
                    //                                    [(int)this->pg_stage_of[pp]];
                    // std::cout << " " << if_cleared << "\n";


                } else {
                    if (kill_size > this->pg_num_parasites_of[pp]){
                        if_cleared = true;
                    }
                }

                // if (kill_size > this->pg_num_parasites_of[pp]){ // pg is cleared
                if (if_cleared){ // pg is cleared
                    //TODO: find alternative reporting methods of p-clearance

                    // if ( this->head_pg_of[ss] != this->tail_pg_of[ss]) {
                    //     if ( pp == (lbound + this->tail_pg_of[ss])) {
                    //         num_clear_tail++;
                    //     } else if ( pp == (lbound + this->head_pg_of[ss]) ) {
                    //         num_clear_head++;
                            
                    //     }else{
                    //         num_clear_other++;
                    //     }
                    //     std::cout << "cleared_from_t: " <<  num_clear_tail << "\n";  
                    //     std::cout << "cleared_from_h: " <<  num_clear_head << "\n";  
                    //     std::cout << "cleared_from_o: " <<  num_clear_other << "\n";  
                    // }

                    this->log_new_clearance_drug[static_cast<int>(this->pg_type_of[pp])]++;

                    if (time_step - this->drug_given_day_of[ss] < 200) {
                        
                        if (this->pg_type_of[pp] == common::ParasiteType::kRb){
                            this->record_of_cleared_pgs_rb[time_step - this->drug_given_day_of[ss]]++;
                        } else {
                            this->record_of_cleared_pgs_r0[time_step - this->drug_given_day_of[ss]]++;
                        }
                    }
                    

                    if (this->kDebug_assert){
                        this->debug_bites_accumulator_less[ss]--;
                    }

                    // std::cout << "pg clearance: ss-" << ss << ", pp-" << pp;

                    this->pg_num_parasites_of[pp] = 0;
                    // this->pg_type_of[pp] = common::ParasiteType::kNoParasite;

                    if (this->pg_stage_of[pp]!=common::StageName::kLiver) {
                        assert(this->moi_of[ss] > 0);
                        this->moi_of[ss]--;
                    }

                    if (this->pg_has_caused_clinicalness_of[pp]) {
                        assert(this->num_of_clinical_pgs_of[ss] > 0);
                        this->num_of_clinical_pgs_of[ss]--;
                        if (this->num_of_clinical_pgs_of[ss] == 0){
                            this->is_clinical_of[ss] = false;
                        }
                    }

                    // queue compaction
                    if (pp == lbound+this->tail_pg_of[ss]){

                        this->tail_pg_of[ss]--;

                        if (pp == lbound+this->head_pg_of[ss]){
                            this->head_pg_of[ss] = -1;
                            this->tail_pg_of[ss] = -1;

                            this->immunity_days_of[ss] = 0; //gains immunity
                            this->immunity_level_of[ss]++;
                            this->cummulative_exposures_of[ss]++;

                            this->is_clinical_of[ss] = false;
                            this->dominant_type_of[ss] = common::ParasiteType::kNoParasite;
                            this->most_advanced_stage_of[ss] = common::StageName::kNotInSystem;

                            this->last_recovery_day_of[ss] = time_step;

                            drug_recovery_counter++;
                        }

                    } else if (pp == lbound+this->head_pg_of[ss]) {

                        this->head_pg_of[ss]++;

                    } else {

                        this->pg_type_of[pp] = this->pg_type_of[lbound+this->tail_pg_of[ss]];
                        this->pg_stage_of[pp] = this->pg_stage_of[lbound+this->tail_pg_of[ss]];
                        this->pg_aquired_on_day_of[pp] = this->pg_aquired_on_day_of[lbound+this->tail_pg_of[ss]];
                        this->pg_num_parasites_of[pp] = this->pg_num_parasites_of[lbound+this->tail_pg_of[ss]];
                        this->pg_num_gametocytes_of[pp] = this->pg_num_gametocytes_of[lbound+this->tail_pg_of[ss]];
                        this->pg_is_drug_target_of[pp] = this->pg_is_drug_target_of[lbound+this->tail_pg_of[ss]];
                        this->pg_has_caused_clinicalness_of[pp] = this->pg_has_caused_clinicalness_of[lbound+this->tail_pg_of[ss]];

                        // this->pg_num_parasites_of[lbound+this->tail_pg_of[ss]] = 0;

                        pp--;

                        this->tail_pg_of[ss]--;
                    }

                } else { //pg survives

                    // this->pg_num_parasites_of[pp] -= kill_size;

                // parasites replicate

                    // this->pg_num_parasites_of[pp] *=
                    //         common::kParasite_replication_rate
                    //                     [(int)this->pg_type_of[pp]]
                    //                     [this->pg_aquired_on_day_of[pp] % common::kParasite_replication_cycle_length];

                // end of day

                }

            }// for pp
            
        }//this->head_pg_of[ss] < 0)

        // this->drug_days_of[ss]++;
    }
    // std::cout << "end with:" << this->pg_num_parasites_of[0] << ","
    //                            << this->pg_num_parasites_of[1] << ","
    //                            << this->pg_num_parasites_of[2] << ","
    //                            << this->pg_num_parasites_of[3] << ","
    //                            << this->pg_num_parasites_of[4] << ","
    //                            << this->pg_num_parasites_of[5] << ","
    //                            << this->pg_num_parasites_of[6] << ","
    //                            << this->pg_num_parasites_of[7] << ";"
    //                            << (int)this->head_pg_of[0] << "," << (int)this->tail_pg_of[0] << std::endl;

    delete[] random_number_drug_effect;

    std::cout << "\tdrug_recovery_counter=" << drug_recovery_counter << "\n";

}

void BloodSystemManager::step_report_drug_failure_by_clinical_status_on_day(
        int current_time_step
    ) {

    for (auto& ii : this->daily_record_of_drug_failure_by_effetive_period) {
        ii = 0;
    }

    for (int ss = 0; ss < this->sum_num_systems; ss++) {
        if (this->drug_given_day_of[ss] > 0) {
            if (this->tail_pg_of[ss] >= 0){
                if (this->most_advanced_stage_of[ss] != common::StageName::kLiver){



                    int num_days_since_drug_given = current_time_step - this->drug_given_day_of[ss] + 1;
                    // + 1 because of the drug_intake >> drug_effect >> drug_failure_check workflow
                    if (num_days_since_drug_given == this->kDrug_failure_if_clinical_on_day) {

                        if (this->last_recovery_day_of[ss] < 0 ||
                            (current_time_step - this->last_recovery_day_of[ss] + 1) > this->kDrug_failure_if_clinical_on_day) {
                        // + 1 because of the drug_intake >> drug_effect >> drug_failure_check workflow

                            // bool at_least_one_current_infection_predates_last_drug = false;
                            bool at_least_one_remaining_drug_target = false;
                            int lbound = ss * kNum_pgs_per_system;
                            for (int pp = lbound + this->head_pg_of[ss]; pp <= lbound + this->tail_pg_of[ss]; pp++ ){
                                // if (this->pg_stage_of[pp] != common::StageName::kLiver) {
                                    
                                //     // if (static_cast<int>(this->pg_aquired_on_day_of[pp]) < this->drug_given_day_of[ss] ) {
                                //     //     std::cout << "ss=" << ss << "\n";
                                //     //     std::cout << "pp=" << pp << "\n";
                                //     //     std::cout << "this->pg_aquired_on_day_of[pp]=" << this->pg_aquired_on_day_of[pp] << "\n";
                                //     //     std::cout << "this->drug_given_day_of[ss]=" << this->drug_given_day_of[ss] << "\n";
                                //     //     std::cout << "this->pg_stage_of[pp]=" << this->pg_stage_of[pp] << "\n";
                                //     //     std::cout << "this->moi_of[ss]=" << unsigned(this->moi_of[ss]) << "\n";
                                //     //     std::cout << "this->cummulative_exposures_of[ss]=" << unsigned(this->cummulative_exposures_of[ss]) << "\n";

                                //     //     // std::cin.ignore();

                                //     //     at_least_one_current_infection_predates_last_drug = true;
                                //     //     continue;
                                //     // }
                                // }
                                if (static_cast<int>(this->pg_aquired_on_day_of[pp]) < this->drug_given_day_of[ss]) {
                                    if (this->pg_is_drug_target_of[pp]) {
                                        at_least_one_remaining_drug_target = true;
                                        // this->print_one(ss);
                                        continue;
                                    }
                                }
                            }

                            // if (at_least_one_current_infection_predates_last_drug) {
                            //     this->daily_record_of_drug_failure_by_effetive_period[static_cast<int>(this->drug_given_of[ss])]++;
                            // }
                            if (at_least_one_remaining_drug_target) {
                                this->daily_record_of_drug_failure_by_effetive_period[static_cast<int>(this->drug_given_of[ss])]++;
                            }
                        } else {
                            // std::cout << "current_time_step=" << current_time_step
                            //             << "last_recovery_day_of[ss]=" << last_recovery_day_of[ss];
                            // std::cin.ignore();
                        }

                    }
                
                }
            }
        }
    }

}

// void BloodSystemManager::pg_queue_remove_one(){

// }

void BloodSystemManager::parasite_mutation(int current_time_step) {
    if (this->kFunc_parasite_mutation_version == 0) {
        this->parasite_mutation_mmc_wp2();
    } else if (this->kFunc_parasite_mutation_version == 1) {
        this->parasite_mutation_mmc_wp2_drug_days_fixation(current_time_step);
    } else if (this->kFunc_parasite_mutation_version == 2) {
        this->parasite_mutation_mmc_wp2_drug_days_fixation_per_site(current_time_step);
    } else {
        std::cout << "Unknown kFunc_parasite_mutation_version : " << kFunc_parasite_mutation_version << std::endl;
        assert(false);
    }
}

void BloodSystemManager::parasite_mutation_mmc_wp2_drug_days_fixation_per_site(int current_time_step) {
    float applicable_mutation_rate = 0.0;

    for(auto& ii : this->log_new_fixation) {
        ii = 0;
    }
    for(auto& ii : this->log_new_mutation) {
        ii = 0;
    }

    for (int ss = 0; ss < this->sum_num_systems; ss++) {
        if (this->head_pg_of[ss] < 0) {
            // nothing in this host
            continue;
        }
        if (this->most_advanced_stage_of[ss] <= common::StageName::kLiver) {
            continue;
        }

        applicable_mutation_rate = this->kMutation_prob_single_resistance;

        if (this->drug_of[ss] == this->kDrug_ab
                || this->drug_of[ss] == this->kDrug_a
                || this->drug_of[ss] == this->kDrug_b
            ) {

            if ( (current_time_step - this->drug_given_day_of[ss]) >= this->kFunc_parasite_mutation_fixation_after_day) {
                applicable_mutation_rate = this->kMutation_prob_single_resistance_w_act;
            }

        }
        // std::cout << "applicable_mutation_rate:" << applicable_mutation_rate << std::endl;


        for (int ii = 0; ii < 2; ii++) {

            common::ParasiteType new_type_to_add = common::ParasiteType::kR0;
            
            if (ii == 0) {
                new_type_to_add = common::ParasiteType::kRa;
                if (!this->kMutation_type_a_enabled) {
                    continue;
                }
            } else if (ii == 1) {
                new_type_to_add = common::ParasiteType::kRb;
                if (!this->kMutation_type_b_enabled) {
                    continue;
                }
            } else {
                std::cout << "parasite_mutation_mmc_wp2_drug_days_fixation_per_site: unknown mutation site.";
                assert(false);
            }

            float rnd_nmb = 0.0;
            util::set_random_numbers_uniform(&rnd_nmb, 1);

            if (rnd_nmb < applicable_mutation_rate ) {
                // mutation happens

                this->log_new_mutation[static_cast<int>(new_type_to_add)]++;

                bool fixation_happened = false;
                int lbound = ss * kNum_pgs_per_system;
                for (int pp = lbound + this->head_pg_of[ss]; pp <= lbound + this->tail_pg_of[ss]; pp++) {


                    if ( (this->pg_stage_of[pp] == common::StageName::kBlood)
                            || (this->pg_stage_of[pp] == common::StageName::kInfectious)
                        ) {

                        if (this->pg_type_of[pp] == common::ParasiteType::kRab) {
                            // already Rab, no further mutation applied
                        } else if (this->pg_type_of[pp] == common::ParasiteType::kR0) {

                            this->pg_type_of[pp] = new_type_to_add;
                            fixation_happened = true;

                        } else if (this->pg_type_of[pp] == common::ParasiteType::kRa && new_type_to_add == common::ParasiteType::kRb) {

                            this->pg_type_of[pp] = common::ParasiteType::kRab;
                            fixation_happened = true;

                        } else if (this->pg_type_of[pp] == common::ParasiteType::kRb && new_type_to_add == common::ParasiteType::kRa) {

                            this->pg_type_of[pp] = common::ParasiteType::kRab;
                            fixation_happened = true;

                        }

                    } else {

                        assert(this->pg_stage_of[pp] == common::StageName::kLiver);

                    }

                } //pp

                if (fixation_happened) {
                    this->log_new_fixation[static_cast<int>(new_type_to_add)]++;
                }

            }

        }

    }
}

void BloodSystemManager::parasite_mutation_mmc_wp2_drug_days_fixation(int current_time_step) {
    float applicable_mutation_rate = 0.0;

    for(auto& ii : this->log_new_mutation) {
        ii = 0;
    }
    for(auto& ii : this->log_new_fixation) {
        ii = 0;
    }

    for (int ss = 0; ss < this->sum_num_systems; ss++) {
        if (this->head_pg_of[ss] < 0) {
            // nothing in this host
            continue;
        }
        if (this->most_advanced_stage_of[ss] <= common::StageName::kLiver) {
            continue;
        }

        applicable_mutation_rate = this->kMutation_prob_single_resistance;

        if (this->drug_of[ss] == this->kDrug_ab
                || this->drug_of[ss] == this->kDrug_a
                || this->drug_of[ss] == this->kDrug_b
            ) {

            if ( (current_time_step - this->drug_given_day_of[ss]) >= this->kFunc_parasite_mutation_fixation_after_day) {
                applicable_mutation_rate = this->kMutation_prob_single_resistance_w_act;
            }

        }

        // std::cout << "applicable_mutation_rate:" << applicable_mutation_rate << std::endl;

        float rnd_nmb = 0.0;
        util::set_random_numbers_uniform(&rnd_nmb, 1);

        // applicable_mutation_rate *= ((this->tail_pg_of[ss]-this->head_pg_of[ss]) + 1);

        if (rnd_nmb < applicable_mutation_rate ) {
            // mutation happend, need to work out the type

            common::ParasiteType new_type_to_add = common::ParasiteType::kR0;

            util::set_random_numbers_uniform(&rnd_nmb, 1);

            if (rnd_nmb < this->kMutation_prob_single_type_a_type_b_ratio/(this->kMutation_prob_single_type_a_type_b_ratio+1.0)) {
                new_type_to_add = common::ParasiteType::kRa;
            } else {
                new_type_to_add = common::ParasiteType::kRb;
            }

            this->log_new_mutation[static_cast<int>(new_type_to_add)]++;

            // std::cout << "new mutation " << new_type_to_add << " with rate " << applicable_mutation_rate << " stage: " << this->most_advanced_stage_of[ss] << "\n";

            bool fixation_happened = false;
            int lbound = ss * kNum_pgs_per_system;
            for (int pp = lbound + this->head_pg_of[ss]; pp <= lbound + this->tail_pg_of[ss]; pp++) {


                if ( (this->pg_stage_of[pp] == common::StageName::kBlood)
                        || (this->pg_stage_of[pp] == common::StageName::kInfectious)
                    ) {

                    // std::cout << "-" << this->pg_type_of[pp] << "-";

                    if (this->pg_type_of[pp] == common::ParasiteType::kRab) {
                        // already Rab, no further mutation applied
                    } else if (this->pg_type_of[pp] == common::ParasiteType::kR0) {

                        this->pg_type_of[pp] = new_type_to_add;
                        fixation_happened = true;
                        // this->log_new_mutation[static_cast<int>(new_type_to_add)]++;

                    } else if (this->pg_type_of[pp] == common::ParasiteType::kRa && new_type_to_add == common::ParasiteType::kRb) {

                        this->pg_type_of[pp] = common::ParasiteType::kRab;
                        fixation_happened = true;
                        // this->log_new_mutation[static_cast<int>(common::ParasiteType::kRab)]++;
                        // this->log_new_mutation[static_cast<int>(new_type_to_add)]++;

                    } else if (this->pg_type_of[pp] == common::ParasiteType::kRb && new_type_to_add == common::ParasiteType::kRa) {

                        this->pg_type_of[pp] = common::ParasiteType::kRab;
                        fixation_happened = true;
                        // this->log_new_mutation[static_cast<int>(common::ParasiteType::kRab)]++;
                        // this->log_new_mutation[static_cast<int>(new_type_to_add)]++;

                    }

                } else {
                    // std::cout << " L ";
                    assert(this->pg_stage_of[pp] == common::StageName::kLiver);

                }

                // std::cout << " | ";


            } //pp
            
            // std::cout << "\n";


            if (fixation_happened) {
                this->log_new_fixation[static_cast<int>(new_type_to_add)]++;
            }
            //     // std::cout << "fixed";
            // } else {
            //     // std::cout << "not fixed";
            //     // std::cin.ignore();

            // }

        }

    }
}

void BloodSystemManager::parasite_mutation_mmc_wp2() { // MMC_WP2, version rev9f
    float applicable_mutation_rate = 0.0;

    for(auto& ii : this->log_new_mutation) {
        ii = 0;
    }

    for (int ss = 0; ss < this->sum_num_systems; ss++) {
        if (this->head_pg_of[ss] < 0) {
            // nothing in this host
            continue;
        }
        
        int lbound = ss * kNum_pgs_per_system;
        
        for (int pp = lbound + this->head_pg_of[ss]; pp <= lbound + this->tail_pg_of[ss]; pp++) {

            applicable_mutation_rate = this->kMutation_prob_single_resistance;

            if (this->pg_has_caused_clinicalness_of[pp]
                ) {
                applicable_mutation_rate = this->kMutation_prob_single_resistance_w_act;
            }


            float rnd_nmb = 0.0;
            util::set_random_numbers_uniform(&rnd_nmb, 1);

            // applicable_mutation_rate *= ((this->tail_pg_of[ss]-this->head_pg_of[ss]) + 1);

            if (rnd_nmb < applicable_mutation_rate ) {
                // mutation happend, need to work out the type

                common::ParasiteType new_type_to_add = common::ParasiteType::kR0;

                util::set_random_numbers_uniform(&rnd_nmb, 1);

                if (rnd_nmb < this->kMutation_prob_single_type_a_type_b_ratio/(this->kMutation_prob_single_type_a_type_b_ratio+1.0)) {
                    new_type_to_add = common::ParasiteType::kRa;
                } else {
                    new_type_to_add = common::ParasiteType::kRb;
                }


                if ( (this->pg_stage_of[pp] == common::StageName::kBlood)
                        || (this->pg_stage_of[pp] == common::StageName::kInfectious)
                    ) {

                    if (this->pg_type_of[pp] == common::ParasiteType::kRab) {
                        // already Rab, no further mutation applied
                    } else if (this->pg_type_of[pp] == common::ParasiteType::kR0) {

                        this->pg_type_of[pp] = new_type_to_add;
                        this->log_new_mutation[static_cast<int>(this->pg_type_of[pp])]++;

                    } else if (this->pg_type_of[pp] == common::ParasiteType::kRa && new_type_to_add == common::ParasiteType::kRb) {

                        this->pg_type_of[pp] = common::ParasiteType::kRab;
                        this->log_new_mutation[static_cast<int>(this->pg_type_of[pp])]++;

                    } else if (this->pg_type_of[pp] == common::ParasiteType::kRb && new_type_to_add == common::ParasiteType::kRa) {

                        this->pg_type_of[pp] = common::ParasiteType::kRab;
                        this->log_new_mutation[static_cast<int>(this->pg_type_of[pp])]++;

                    }

                } else {

                    assert(this->pg_stage_of[pp] == common::StageName::kLiver);

                }


            }

        } // pp
    }
}

void BloodSystemManager::parasite_mutation_mmc_wp2_drug_trigger() { // MMC_WP2, pre rev9f
    float applicable_mutation_rate = 0.0;

    for(auto& ii : this->log_new_mutation) {
        ii = 0;
    }

    for (int ss = 0; ss < this->sum_num_systems; ss++) {
        if (this->head_pg_of[ss] < 0) {
            // nothing in this host
            continue;
        }

        applicable_mutation_rate = this->kMutation_prob_single_resistance;

        if (this->drug_of[ss] == this->kDrug_ab
                || this->drug_of[ss] == this->kDrug_a
                || this->drug_of[ss] == this->kDrug_b
            ) {

            applicable_mutation_rate = this->kMutation_prob_single_resistance_w_act;

        // std::cout << "applicable_mutation_rate:" << applicable_mutation_rate << std::endl;
        }


        float rnd_nmb = 0.0;
        util::set_random_numbers_uniform(&rnd_nmb, 1);

        // applicable_mutation_rate *= ((this->tail_pg_of[ss]-this->head_pg_of[ss]) + 1);

        if (rnd_nmb < applicable_mutation_rate ) {
            // mutation happend, need to work out the type

            common::ParasiteType new_type_to_add = common::ParasiteType::kR0;

            util::set_random_numbers_uniform(&rnd_nmb, 1);

            if (rnd_nmb < this->kMutation_prob_single_type_a_type_b_ratio/(this->kMutation_prob_single_type_a_type_b_ratio+1.0)) {
                new_type_to_add = common::ParasiteType::kRa;
            } else {
                new_type_to_add = common::ParasiteType::kRb;
            }

            int lbound = ss * kNum_pgs_per_system;
            for (int pp = lbound + this->head_pg_of[ss]; pp <= lbound + this->tail_pg_of[ss]; pp++) {

                if ( (this->pg_stage_of[pp] == common::StageName::kBlood)
                        || (this->pg_stage_of[pp] == common::StageName::kInfectious)
                    ) {

                    if (this->pg_type_of[pp] == common::ParasiteType::kRab) {
                        // already Rab, no further mutation applied
                    } else if (this->pg_type_of[pp] == common::ParasiteType::kR0) {

                        this->pg_type_of[pp] = new_type_to_add;
                        this->log_new_mutation[static_cast<int>(this->pg_type_of[pp])]++;

                    } else if (this->pg_type_of[pp] == common::ParasiteType::kRa && new_type_to_add == common::ParasiteType::kRb) {

                        this->pg_type_of[pp] = common::ParasiteType::kRab;
                        this->log_new_mutation[static_cast<int>(this->pg_type_of[pp])]++;

                    } else if (this->pg_type_of[pp] == common::ParasiteType::kRb && new_type_to_add == common::ParasiteType::kRa) {

                        this->pg_type_of[pp] = common::ParasiteType::kRab;
                        this->log_new_mutation[static_cast<int>(this->pg_type_of[pp])]++;

                    }

                } else {

                    assert(this->pg_stage_of[pp] == common::StageName::kLiver);

                }


            }

            // std::cout << "drug:" << this->drug_of[ss] << "\n"
            //           << "clinical:" << this->is_clinical_of[ss] << "\n"
            //           << "dominant:" << this->dominant_type_of[ss] << "\n" 
            //             << "rate:" << applicable_mutation_rate << "\n"
            //             << "type to add:" << new_type_to_add
            //             << std::endl;
            // std::cin.ignore();

        }

    }
}

void BloodSystemManager::parasite_mutation_not_used(uint16_t time_step) {
    // const uint32_t kMutation_until = 16382;
    const uint32_t kMutation_until = 1000000;

    int mutation_counter = 0;

    for (int ss = 0; ss < this->sum_num_systems; ss++) {
        if (this->head_pg_of[ss] < 0) {

        } else {
            int lbound = ss * kNum_pgs_per_system;
            
            for ( int pp = lbound + this->head_pg_of[ss]; pp <= lbound + this->tail_pg_of[ss]; pp++) {
                
                if (this->pg_type_of[pp] != common::ParasiteType::kR0) {

                } else if (this->pg_num_gametocytes_of[pp] > 0) {


                    if (this->pg_num_gametocytes_of[pp] > kMutation_until) {

                        std::cout << "new pg from mutation on day " << static_cast<int>(time_step) << "\n";
                        
                        common::ParasiteType pt = common::ParasiteType::kR0;
                        
                        if(this->drug_of[ss] == common::DrugName::kArtesunate) {
                            pt = common::ParasiteType::kRa;
                        } else if (this->drug_of[ss] == common::DrugName::kPiperaquine) {
                            pt = common::ParasiteType::kRb;
                        } else {
                            float rnd_nmb_rtype = 0.0;
                            util::set_random_numbers_uniform(&rnd_nmb_rtype, 1);
                            if (rnd_nmb_rtype < 0.5) {
                                pt = common::ParasiteType::kRa;
                            } else {
                                pt = common::ParasiteType::kRb;
                            }
                        }
                        
                        if (this->pg_stage_of[pp] != common::StageName::kLiver) {
                            this->moi_of[ss]++;
                        }
                        this->pg_num_gametocytes_of[pp] = 0;
                        int new_pg_index = this->add_one_pg_with_defaults(
                            ss,
                            pt,
                            this->pg_stage_of[pp],
                            time_step
                        );
                        (void) new_pg_index;
                        break;
                        
                    } else {
                        if (this->pg_num_gametocytes_of[pp] % 2 == 0) {
                            this->pg_num_gametocytes_of[pp] += 1;
                        } else {
                            this->pg_num_gametocytes_of[pp] *= 2;
                        }
                    }


                } else {

                    float rnd_nmb = 0.0;
                    util::set_random_numbers_uniform(&rnd_nmb, 1);
                    float mutation_rate = 0.0;
                    // if (this->pg_has_caused_clinicalness_of[pp]) {
                    if (this->is_clinical_of[ss] ) {
                        mutation_rate = 0.1;
                    } else {
                        mutation_rate = 0.001;
                    }

                    if (rnd_nmb < mutation_rate) {
                        this->pg_num_gametocytes_of[pp] = 1;

                        mutation_counter++;

                        // std::cout << "mutation happend\n";
                        // std::cin.ignore();
                    }

                }


            }
        }
    }

    std::cout << "# mutations: " << mutation_counter << "\n";
    // std::cin.ignore();
}

void BloodSystemManager::parasite_stage_progression(
        const int time_step
    ){
    int natural_recovery_counter = 0; // full recovery

    for (auto& ii: this->log_new_clearance_natural) {
        ii = 0;
    }
    for (auto& ii: this->daily_record_of_new_cases_clinical) {
        ii = 0;
    }
    for (auto& ii: this->daily_record_of_new_cases_asymptomatic) {
        ii = 0;
    }

    for (int ss = 0; ss < this->sum_num_systems; ss++) {
        this->report_daily_new_clinical_case[ss] = false;
        this->report_daily_new_asymptomatic_case[ss] = false;
        if (this->head_pg_of[ss] < 0){
            // no pg in this system
        } else {
            int lbound = ss*kNum_pgs_per_system;
            bool if_recovery_happend = false;
            for (int pp = lbound + this->head_pg_of[ss]; pp <= lbound+this->tail_pg_of[ss]; pp++){
                if(this->parasite_stage_progression_probability_based(pp, ss)) {
                    if_recovery_happend = true;
                }
            }

            // check for naturally-recovered ones and compact pg queue
            if (if_recovery_happend) {

                for (int pp = lbound +  this->head_pg_of[ss]; pp <= lbound + this->tail_pg_of[ss]; pp++){
                    if (this->pg_stage_of[pp] == common::StageName::kNotInSystem) {

                        this->log_new_clearance_natural[static_cast<int>(this->pg_type_of[pp])]++;

                        if (this->kDebug_assert){
                            this->debug_bites_accumulator_less[ss]--;
                        }

                        if (this->pg_stage_of[pp]!=common::StageName::kLiver) {
                            assert(this->moi_of[ss] > 0);
                            this->moi_of[ss]--;
                        }

                        if (this->pg_has_caused_clinicalness_of[pp]) {
                            assert(this->num_of_clinical_pgs_of[ss] > 0);
                            this->num_of_clinical_pgs_of[ss]--;
                            if (this->num_of_clinical_pgs_of[ss] == 0){
                                this->is_clinical_of[ss] = false;
                            }
                        }

                        // queue compaction - natural recovery
                        if (pp == lbound+this->tail_pg_of[ss]){

                            this->tail_pg_of[ss]--;

                            if (pp == lbound+this->head_pg_of[ss]){
                                this->head_pg_of[ss] = -1;
                                this->tail_pg_of[ss] = -1;

                                this->immunity_days_of[ss] = 0; //gains immunity
                                this->immunity_level_of[ss]++;
                                this->cummulative_exposures_of[ss]++;

                                this->is_clinical_of[ss] = false;
                                this->dominant_type_of[ss] = common::ParasiteType::kNoParasite;
                                this->most_advanced_stage_of[ss] = common::StageName::kNotInSystem;

                                this->last_recovery_day_of[ss] = time_step;
                                
                                natural_recovery_counter++;
                            }

                        } else if (pp == lbound+this->head_pg_of[ss]) {

                            this->head_pg_of[ss]++;

                        } else {

                            this->pg_type_of[pp] = this->pg_type_of[lbound+this->tail_pg_of[ss]];
                            this->pg_stage_of[pp] = this->pg_stage_of[lbound+this->tail_pg_of[ss]];
                            this->pg_aquired_on_day_of[pp] = this->pg_aquired_on_day_of[lbound+this->tail_pg_of[ss]];
                            this->pg_num_parasites_of[pp] = this->pg_num_parasites_of[lbound+this->tail_pg_of[ss]];
                            this->pg_num_gametocytes_of[pp] = this->pg_num_gametocytes_of[lbound+this->tail_pg_of[ss]];
                            this->pg_is_drug_target_of[pp] = this->pg_is_drug_target_of[lbound+this->tail_pg_of[ss]];
                            this->pg_has_caused_clinicalness_of[pp] = this->pg_has_caused_clinicalness_of[lbound+this->tail_pg_of[ss]];

                            // this->pg_num_parasites_of[lbound+this->tail_pg_of[ss]] = 0;

                            pp--;

                            this->tail_pg_of[ss]--;
                        }

                    }
                }
            }
        }
    }
    std::cout << "\tnatural_recovery_counter=" << natural_recovery_counter << "\n";
}
bool BloodSystemManager::parasite_stage_progression_probability_based(
        const int p_index,
        const int s_index
        ){

    // const float kPbl_liver_to_blood = 1.0/5.0;//G_gamma = 1.0/5.0;
    // const float kPbl_blood_to_infectious = 1.0/15.0;//G_sigma = 1.0/15.0;

    const float kPbl_liver_to_blood = 1.0 / this->kStage_progression_l2b_days_mean;
    const float kPbl_blood_to_infectious = 1.0 / this->kStage_progression_b2i_days_mean;



    bool if_recovery_happend = false;

    float random_number = 0.0;
    util::set_random_numbers_uniform(&random_number, 1);

    switch (this->pg_stage_of[p_index]) {
        case common::StageName::kLiver :
            if (random_number < kPbl_liver_to_blood) {
                this->pg_stage_of[p_index] = common::StageName::kBlood;
                this->report_stage_progression_flags[p_index] = true;
                this->moi_of[s_index]++;
            }
            break;
        case common::StageName::kBlood :
            if (random_number < kPbl_blood_to_infectious) {

                this->pg_stage_of[p_index] = common::StageName::kInfectious;
                this->report_stage_progression_flags[p_index] = true;


                util::set_random_numbers_uniform(&random_number, 1);
                bool to_cause_a_clinical_case = 
                        (random_number < this->func_clinical_probability_blood_to_infectious(
                                            this->immunity_level_of[s_index],
                                            this->cummulative_exposures_of[s_index],
                                            this->moi_of[s_index]
                                        )
                        );
                if (to_cause_a_clinical_case) {
                    this->pg_has_caused_clinicalness_of[p_index] = true;
                    this->num_of_clinical_pgs_of[s_index]++;
                    this->daily_record_of_new_cases_clinical[static_cast<int>(this->pg_type_of[p_index])]++;

                    this->report_daily_new_clinical_case[s_index] = true;
                } else {
                    this->daily_record_of_new_cases_asymptomatic[static_cast<int>(this->pg_type_of[p_index])]++;
                    this->report_daily_new_asymptomatic_case[s_index] = true;
                }

                if (this->is_clinical_of[s_index] == false){
                    this->is_clinical_of[s_index] = to_cause_a_clinical_case;
                }

            }
            break;
        case common::StageName::kInfectious :  {// natural recovery
                // float applicable_natural_recovery_rate = this->kNatural_recovery_rate;
                float applicable_natural_recovery_rate = 1.0 / this->kStage_progression_i2r_days_mean;
                if (this->pg_type_of[p_index] != common::ParasiteType::kR0) { // pg carries resistance
                    applicable_natural_recovery_rate *= (1.0 + this->kNatural_recovery_resistance_cost);
                }
                if (random_number < applicable_natural_recovery_rate) {
                    this->pg_stage_of[p_index] = common::StageName::kNotInSystem;
                    if_recovery_happend = true;
                }
            }
            break;
        default:
            break;
    }

    return if_recovery_happend;
}

float BloodSystemManager::func_clinical_probability_blood_to_infectious(
        uint8_t immunity_level,
        uint8_t cummulative_exposures,
        uint8_t moi
    ) {

    assert(moi >= 1);

    switch(this->kFunc_clinical_probability_blood_to_infectious_version) {
        case 0 :
                return this->func_clinical_probability_blood_to_infectious_kernel_0(
                    immunity_level,
                    cummulative_exposures,
                    moi
                );
                break;
        case 1 :
                return this->func_clinical_probability_blood_to_infectious_kernel_1(
                    immunity_level,
                    cummulative_exposures,
                    moi
                );
                break;
        case 2 :
                return this->func_clinical_probability_blood_to_infectious_kernel_2(
                    immunity_level,
                    cummulative_exposures
                );
                break;
        case 3 :
                return this->func_clinical_probability_blood_to_infectious_kernel_3(
                    immunity_level,
                    cummulative_exposures,
                    this->kFunc_clinical_probability_blood_to_infectious_version3_para1,
                    this->kFunc_clinical_probability_blood_to_infectious_version3_para2
                );
                break;
        case 4 :
                return this->func_clinical_probability_blood_to_infectious_kernel_4(
                    immunity_level,
                    cummulative_exposures,
                    moi
                );
                break;
        default:
                std::cout << "Warning: un-recognised kFunc_clinical_probability_blood_to_infectious_version" << std::endl;
                return this->func_clinical_probability_blood_to_infectious_kernel_0(
                    immunity_level,
                    cummulative_exposures,
                    moi
                );
                break;
    }

}

float BloodSystemManager::func_clinical_probability_blood_to_infectious_kernel_4(
        uint8_t immunity_level,
        uint8_t cummulative_exposures,
        uint8_t moi
    ) {
// mda_coverage, review, group g, revision 5

    float clinical_probability = (
            0.55 * exp( -(cummulative_exposures - 1) * 0.1)
            + exp( -0.9 * (cummulative_exposures + 1))
            )
            / sqrt(immunity_level);

    if (moi > 1) {
        return clinical_probability * exp( -0.15 * (moi - 1) );
    } else {
        return clinical_probability;
    }  

}



float BloodSystemManager::func_clinical_probability_blood_to_infectious_kernel_3(
        uint8_t immunity_level,
        uint8_t cummulative_exposures,
        float para1, //0.2
        float para2  //0.95
    ) {
// mmc_wp2, rev09g
// parameterised kernel_2

    float G_immunity = 0.7;
    float clinical_probability = (1-(para1*(immunity_level)))*exp(-G_immunity*cummulative_exposures*(immunity_level-para2)/10)+(1-(immunity_level/10))*0.05;
    if (clinical_probability < 0.05) {
        clinical_probability = 0.05;
    }

    // std::cout << "clinical_probability=" << clinical_probability << "\n";

    return clinical_probability;

}

float BloodSystemManager::func_clinical_probability_blood_to_infectious_kernel_2(
        uint8_t immunity_level,
        uint8_t cummulative_exposures
    ) {
// mmc_wp2, 627f39f43f8ab059e26918f3db1df23e49e79548
// deme->probclinicallb=(1-(0.4*(level)))*exp(-G_immunity*cumm*(level-0.15)/10)+(1-(level/10))*0.05;
//    if(deme->probclinicallb<0.05){deme->probclinicallb=0.05;}

    float G_immunity = 0.7;
    float clinical_probability = (1-(0.2*(immunity_level)))*exp(-G_immunity*cummulative_exposures*(immunity_level-0.95)/10)+(1-(immunity_level/10))*0.05;
    // float clinical_probability = (1-(0.2*(immunity_level)))*exp(-G_immunity*cummulative_exposures*(immunity_level-0.15)/10)+(1-(immunity_level/10))*0.05;
    // float clinical_probability = (1-(0.4*(immunity_level)))*exp(-G_immunity*cummulative_exposures*(immunity_level-0.15)/10)+(1-(immunity_level/10))*0.05;
    if (clinical_probability < 0.05) {
        clinical_probability = 0.05;
    }

    // std::cout << "clinical_probability=" << clinical_probability << "\n";

    return clinical_probability;



}

float BloodSystemManager::func_clinical_probability_blood_to_infectious_kernel_0(
        uint8_t immunity_level,
        uint8_t cummulative_exposures,
        uint8_t moi
    ) {
// karenjan17.cpp
// deme->probclinicallb=(0.1*exp(-(cumm-2)*0.1)+exp((-0.9*cumm)))/sqrt(level);
// deme->probclinicalbb=exp(-0.15*(moi-1))*deme->probclinicallb;

// deme->infections.push_back(aux2);
// if (deme->moi>1){
//     pci=ran2();
//     if (pci<= deme->probclinicalbb){
//         deme->clinical = CLINICAL;
//         deme->clintoday=1;
//     }
//     else{deme->clinical = ASYMPTOMATIC;
//     deme->clintoday=0;}
// }
// else if (deme->moi == 1){
//     pci=ran2();
//     if (pci<= deme->probclinicallb){
//         deme->clinical = CLINICAL;
//         deme->clintoday=1;
//     }
//     else{deme->clinical = ASYMPTOMATIC;
//     deme->clintoday=0;}
// }

    float clinical_probability = (
            0.1 * exp( -(cummulative_exposures - 2) * 0.1)
            + exp( -0.9 * cummulative_exposures)
            )
            / sqrt(immunity_level);

    if (moi > 1) {
        return clinical_probability * exp( -0.15 * (moi - 1) );
    } else {
        return clinical_probability;
    }

}

float BloodSystemManager::func_clinical_probability_blood_to_infectious_kernel_1(
        uint8_t immunity_level,
        uint8_t cummulative_exposures,
        uint8_t moi
    ) {
// drugrestest.cpp
// deme->probclinicallb=(0.25*exp(-(cumm-1))+exp((-0.6*cumm)))/level;
// deme->probclinicalbb=exp(-0.15*(moi-1))*deme->probclinicallb;

// deme->infections.push_back(aux2);
// if (deme->moi>1){
//     pci=ran2();
//     if (pci<= deme->probclinicalbb){
//         deme->clinical = CLINICAL;
//     }
//     else{deme->clinical = ASYMPTOMATIC;}
// }
// else if (deme->moi == 1){
//     pci=ran2();
//     if (pci<= deme->probclinicallb){
//         deme->clinical = CLINICAL;
//     }
//     else{deme->clinical = ASYMPTOMATIC;}
// }
    float clinical_probability = (
            0.25 * exp( -(cummulative_exposures - 1) )
            +exp( (-0.6 * cummulative_exposures) )
            )
            / immunity_level;

    if (moi > 1) {
        return clinical_probability * exp( -0.15 * (moi - 1) );
    } else {
        return clinical_probability;
    }

}

void BloodSystemManager::parasite_dominance(){

    for(const auto& pt: util::Enum_iterator<common::ParasiteType>()) {
        // this->sum_prevalence_count[(int) pt] = 0;
        this->sum_prevalence_count[static_cast<int>(pt)] = 0;
    }
    this->sum_clinical_count = 0;
    this->sum_asymptomatic_count = 0;

    for (int ss = 0; ss < this->sum_num_systems; ss++) {

        this->dominant_type_of[ss] = common::ParasiteType::kNoParasite;
        this->most_advanced_stage_of[ss] = common::StageName::kNotInSystem;

        if (this->head_pg_of[ss] < 0){
            // no pg in this system
            this->dominant_type_of[ss] = common::ParasiteType::kNoParasite;
            this->most_advanced_stage_of[ss] = common::StageName::kNotInSystem;
        } else {

            if (this->kFunc_dominance_version == 0) {

                this->parasite_dominance_probability_based( // mda_coverage
                    ss*kNum_pgs_per_system + this->head_pg_of[ss],
                    ss*kNum_pgs_per_system + this->tail_pg_of[ss],
                    ss);

            } else if (this->kFunc_dominance_version == 1) {
                
                this->parasite_dominance_mmc_wp2_kernel(
                    ss * kNum_pgs_per_system + this->head_pg_of[ss],
                    ss * kNum_pgs_per_system + this->tail_pg_of[ss],
                    ss,
                    this->dominant_type_of,
                    this->most_advanced_stage_of,
                    this->pg_stage_of,
                    this->pg_type_of,
                    this->drug_of,
                    this->kDrug_ab,
                    this->kDrug_a,
                    this->kDrug_b
                );

            } else {
                std::cout << "Unknown version of the dominance function: "
                            << this->kFunc_dominance_version << std::endl;
                assert(false);
            }
        }

        if (this->dominant_type_of[ss] == common::ParasiteType::kNoParasite) {
            this->is_clinical_of[ss] = false;
        }

        if (this->dominant_type_of[ss] != common::ParasiteType::kNoParasite) {
            if (this->is_clinical_of[ss]) {
                this->sum_clinical_count++;
            } else {
                this->sum_asymptomatic_count++;
            }
        }

        this->sum_prevalence_count[(int)this->dominant_type_of[ss]]+=1;
    }
}

void BloodSystemManager::step_report_genotype_weights() {
    for (auto& ii: this->daily_record_of_genotype_weights) {
        ii = 0.0;
    }
    this->daily_record_of_parasite_positive_systems = 0;
    this->daily_record_of_num_pgs = 0;
    this->daily_record_of_avg_clinical_pgs = 0.0;
    int num_non_zero_clinical_pg_hosts = 0;

    for (int ss = 0; ss < this->sum_num_systems; ss++) {
        std::vector<int> per_person_counter(static_cast<int>(common::ParasiteType::Last)+1,0);
        if (this->head_pg_of[ss] >= 0) {

            this->daily_record_of_parasite_positive_systems++;
            this->daily_record_of_num_pgs += (this->tail_pg_of[ss] - this->head_pg_of[ss] + 1);

            int head = ss*kNum_pgs_per_system + this->head_pg_of[ss];
            int tail = ss*kNum_pgs_per_system + this->tail_pg_of[ss];
            for (int pp = head; pp <= tail; pp++) {
                per_person_counter[static_cast<int>(this->pg_type_of[pp])]++;
            }
            for (const auto& pt: util::Enum_iterator<common::ParasiteType>()) {

                float per_person_weight = (float)per_person_counter[static_cast<int>(pt)]
                    / (float)((int)this->tail_pg_of[ss] - (int)this->head_pg_of[ss] + 1);
                
                this->daily_record_of_genotype_weights[static_cast<int>(pt)] += per_person_weight;
                
                // if (this->dominant_type_of[ss] == common::ParasiteType::kRb) {

                //     std::cout << pt << " weight:" << this->daily_record_of_genotype_weights[static_cast<int>(pt)] << "\n";
                    
                //     std::cout << "per_person_counter[" << pt << "]=" << per_person_counter[(int)pt] << "\n";
                //     std::cout << "per_person_total" << (float)((int)this->tail_pg_of[ss] - (int)this->head_pg_of[ss] + 1) << "\n";

                    
                //     std::cin.ignore();
                // }
                // std::cout << "per_person_weight=" << per_person_weight << "\n";
            }
        }
        if (this->num_of_clinical_pgs_of[ss] != 0) {
            num_non_zero_clinical_pg_hosts++;
            this->daily_record_of_avg_clinical_pgs += this->num_of_clinical_pgs_of[ss];
        }
    }
    this->daily_record_of_avg_clinical_pgs /= (float)num_non_zero_clinical_pg_hosts;
}

void BloodSystemManager::step_report_drug_composition() {
    for (auto& ii: this->daily_record_of_drug_current) {
        ii = 0;
    }
    for (auto& ii: this->daily_record_of_drug_current_in_infected_population) {
        ii = 0;
    }
    for (int ss = 0; ss < this->sum_num_systems; ss++) {
        this->daily_record_of_drug_current[static_cast<int>(this->drug_of[ss])]++;
        if (this->tail_pg_of[ss] >= 0) {
            this->daily_record_of_drug_current_in_infected_population[static_cast<int>(this->drug_of[ss])]++;
        }
    }
}

void BloodSystemManager::parasite_dominance_mmc_wp2_kernel(
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
    ){

    assert(head_pg_index <= tail_pg_index);
    assert(system_index <= head_pg_index);

    // reset current values
    dominant_type_of_array[system_index] = common::ParasiteType::kNoParasite;
    most_advanced_stage_of_array[system_index] = common::StageName::kNotInSystem;

    // bool any_infectious_resistance_type = false;
    int num_infectious_r0 = 0;
    bool any_infectious_ra = false;
    bool any_infectious_rb = false;
    bool any_infectious_rab = false;

    static std::vector<common::ParasiteType> list_of_infectious_resistant_pts;
    list_of_infectious_resistant_pts.clear();

    for (int pp = head_pg_index; pp <= tail_pg_index; pp++){
        
        assert(pg_stage_of_array[pp] != common::StageName::kNotInSystem);
        assert(pg_type_of_array[pp] != common::ParasiteType::kNoParasite);
        
        if (pg_stage_of_array[pp] > most_advanced_stage_of_array[system_index]){
            most_advanced_stage_of_array[system_index] = pg_stage_of_array[pp];
        }

        if (pg_stage_of_array[pp] == common::StageName::kInfectious){

            if (pg_type_of_array[pp] != common::ParasiteType::kR0) {
                // any_infectious_resistance_type = true;
                list_of_infectious_resistant_pts.push_back(pg_type_of_array[pp]);
                if (pg_type_of_array[pp] == common::ParasiteType::kRa) {
                    any_infectious_ra = true;
                } else if (pg_type_of_array[pp] == common::ParasiteType::kRb) {
                    any_infectious_rb = true;
                } else if (pg_type_of_array[pp] == common::ParasiteType::kRab) {
                    any_infectious_rab = true;
                } else {
                    assert(false && "unknown resistance type");
                }

            } else {
                num_infectious_r0 += 1;
            }
        }
    }

    if (most_advanced_stage_of_array[system_index] < common::StageName::kInfectious) {
        // no infectious stage parasites
        dominant_type_of_array[system_index] = common::ParasiteType::kNoParasite;

    } else if (list_of_infectious_resistant_pts.empty()) {
        // infectious, but all R0
        dominant_type_of_array[system_index] = common::ParasiteType::kR0;

    // infectious and has at least one resistance type, check drugs
    } else if (drug_of_array[system_index] == drug_ab) {
        // Drug_AB: AB > A > B > 0

        if (any_infectious_rab) {
            dominant_type_of_array[system_index] = common::ParasiteType::kRab;
        } else if (any_infectious_ra) {
            dominant_type_of_array[system_index] = common::ParasiteType::kRa;
        } else {
            assert(any_infectious_rb);
            dominant_type_of_array[system_index] = common::ParasiteType::kRb;
        }

    } else if (drug_of_array[system_index] == drug_b) {
        // Drug_B: AB > B > A > 0

        if (any_infectious_rab) {
            dominant_type_of_array[system_index] = common::ParasiteType::kRab;
        } else if (any_infectious_rb) {
            dominant_type_of_array[system_index] = common::ParasiteType::kRb;
        } else {
            assert(any_infectious_ra);
            dominant_type_of_array[system_index] = common::ParasiteType::kRa;
        }

    } else if (drug_of_array[system_index] == drug_a) {
        // Drug_A: AB > A > B > 0
        // assuming there is situations that drug_a alone exist in blood
        // check against kPossible_to_have_drug_a_only in other functions

        if (any_infectious_rab) {
            dominant_type_of_array[system_index] = common::ParasiteType::kRab;
        } else if (any_infectious_ra) {
            dominant_type_of_array[system_index] = common::ParasiteType::kRa;
        } else {
            assert(any_infectious_rb);
            dominant_type_of_array[system_index] = common::ParasiteType::kRb;
        }

    } else {
        // No drug / drug other than AB/A/B: R0 > A/B/AB

        // if (num_infectious_r0 > 0) {
        //     dominant_type_of_array[system_index] = common::ParasiteType::kR0;
        // } else {
        //     float random_number = 0.0;
        //     util::set_random_numbers_uniform(&random_number, 1);
        //     dominant_type_of_array[system_index] =
        //         list_of_infectious_resistant_pts.at(
        //             static_cast<int>(random_number*list_of_infectious_resistant_pts.size())
        //         );
        // }


        // No drug / drug other than AB/A/B: R0/A/B/AB
        // rev8c

        for (int ii = 0; ii < num_infectious_r0; ii++) {
            list_of_infectious_resistant_pts.push_back(common::ParasiteType::kR0);
        }

        float random_number = 0.0;
        util::set_random_numbers_uniform(&random_number, 1);

        dominant_type_of_array[system_index] = 
            list_of_infectious_resistant_pts.at(
                    static_cast<int>(random_number*list_of_infectious_resistant_pts.size())
            );

    }


}

void BloodSystemManager::parasite_dominance_probability_based(
    const int head_pg_index,
    const int tail_pg_index,
    const int system_index
    ){

    this->dominant_type_of[system_index] = common::ParasiteType::kNoParasite;
    this->most_advanced_stage_of[system_index] = common::StageName::kNotInSystem;

    bool there_is_an_infectious_resistance_type = false;
    // bool there_is_an_infectious_type = false;

    static std::vector<common::ParasiteType> pt;
    pt.clear();
    // std::vector<int> pt_counter(static_cast<int>(common::ParasiteType::Last)+1,0);

    for (int pp = head_pg_index; pp <= tail_pg_index; pp++){

        if (this->pg_stage_of[pp] > this->most_advanced_stage_of[system_index]){
            this->most_advanced_stage_of[system_index] = this->pg_stage_of[pp];
        }

        if (this->pg_stage_of[pp] == common::StageName::kInfectious){
            // there_is_an_infectious_type = true;
            pt.push_back(this->pg_type_of[pp]);

            if (this->pg_type_of[pp] != common::ParasiteType::kR0) {
                there_is_an_infectious_resistance_type = true;
            }
        }

        // pt_counter[static_cast<int>(this->pg_type_of[pp])]++;

        // if (this->pg_stage_of[pp] == common::StageName::kInfectious &&
        //     this->pg_type_of[pp] != common::ParasiteType::kR0) {

        //     there_is_an_infectious_resistance_type = true;
        //     // break;
        // }
    }

    // if (!there_is_an_infectious_type) {
    if (pt.empty()) {
        this->dominant_type_of[system_index] = common::ParasiteType::kNoParasite;
    } else if (there_is_an_infectious_resistance_type){

        if (this->drug_of[system_index]==common::DrugName::kArtesunate
            || this->drug_of[system_index]==common::DrugName::kPiparte // check drug_given_of?
            ){
            this->dominant_type_of[system_index] = common::ParasiteType::kRa;

        // if (this->drug_of[system_index]==common::DrugName::kPiperaquine
        //     || this->drug_of[system_index]==common::DrugName::kPiparte // check drug_given_of?
        // // if (this->drug_given_of[system_index]==common::DrugName::kPiperaquine
        // //     || this->drug_given_of[system_index]==common::DrugName::kPiparte // check drug_given_of?
        //     ){
        //     this->dominant_type_of[system_index] = common::ParasiteType::kRb;

        // if (this->drug_of[system_index]==common::DrugName::kArtesunate) {
        //     this->dominant_type_of[system_index] = common::ParasiteType::kRa;
        // } else if (this->drug_of[system_index]==common::DrugName::kPiperaquine) {
        //     this->dominant_type_of[system_index] = common::ParasiteType::kRb;
        } else {
            float random_number = 0.0;
            util::set_random_numbers_uniform(&random_number, 1);
            this->dominant_type_of[system_index]
                // =this->pg_type_of[head_pg_index + int(random_number*(tail_pg_index - head_pg_index + 1))];
                = pt.at(int(random_number*pt.size()));
            // if (pt.size() > 1) {
            //     std::cout << "this->dominant_type_of[system_index]=" << this->dominant_type_of[system_index] << "\n";
            //     std::cin.ignore();
            // }
        }

    } else {
        this->dominant_type_of[system_index] = common::ParasiteType::kR0;
    }

    // Check for natural recoveries

    // if (G_demes[d]->infectious_qeue.empty()==true && G_demes[d]->infections.empty()==true){
    //     G_demes[d]->moi=0;
    //     rec+=1;
    //     G_demes[d]->transmission=SUSCEPTIBLE;
    //     G_demes[d]->immunity_level= G_demes[d]->immunity_level+1;
    //     G_demes[d]->immunity_days=1;
    //     G_demes[d]->immunity=IMMUNE;
    //     G_demes[d]->resistance=RESISTRO;
    //     G_demes[d]->clinical = NOINFECTION;
    // }

}

template <typename T>
common::ParasiteType BloodSystemManager::parasite_dominance_density_based(
    const int head_pg_index,
    const int tail_pg_index,
    const T* densities,
    const common::ParasiteType* types
    ){

    //TODO check stage, only infectious pgs count
    std::vector<T> discrete_weights(tail_pg_index - head_pg_index + 1);

    for (int pp = head_pg_index; pp <= tail_pg_index; pp++){
        discrete_weights[pp]=densities[pp];
    }

    std::random_device rd; // TODO: update this
    std::mt19937 gen(rd());
    std::discrete_distribution<> d(discrete_weights.begin(), discrete_weights.end());

    return types[d(gen)];

}
template common::ParasiteType BloodSystemManager::parasite_dominance_density_based<float>(
    const int head_pg_index,
    const int tail_pg_index,
    const float* densities,
    const common::ParasiteType* types
);

void BloodSystemManager::step_immunity_loss(
        const float* uniform_random_number_array,
        const int uniform_random_number_array_size
    ) {

    this->step_immunity_loss_kernel(

        this->immunity_level_of,
        this->immunity_days_of,
        this->sum_num_systems,

        this->kImmunity_days_until_loss_start,
        this->kImmunity_level_until_loss_start,
        this->kImmunity_loss_rate,

        uniform_random_number_array,
        uniform_random_number_array_size
    );
}

void BloodSystemManager::step_immunity_loss_kernel(

        uint8_t* immunity_level_of_array,
        int* immunity_days_of_array,
        const int data_array_size_sum_num_systems,

        const int immunity_days_until_loss_start,
        const int immunity_level_until_loss_start,
        const float immunity_loss_rate, 

        const float* uniform_random_number_array,
        const int uniform_random_number_array_size

    ) {
    // /** immunity loss **/
    // if (deme->immunity==IMMUNE && deme->immunity_days>40 && deme->immunity_level>1 ){
    //     alpha=ran2();
    //     if (alpha<=G_alpha){
    //         deme->immunity_level=deme->immunity_level-1;
    //         if (deme->immunity_level==0){
    //             deme->immunity==NONIMMUNE;
    //             G_groups->SUM_IMM = G_groups->SUM_IMM - 1;
    //             G_groups->SUM_NIMM = G_groups->SUM_NIMM + 1;
    //             deme->immunity_days=0;
    //         }
    //     }
    //  }

    // /** immunity count **/
    // if (G_demes[d]->immunity==IMMUNE){
    //     G_demes[d]->immunity_days=G_demes[d]->immunity_days+1;
    // }

    assert(uniform_random_number_array_size == data_array_size_sum_num_systems);

    for (int ss = 0; ss < data_array_size_sum_num_systems; ss++) {

        if (immunity_level_of_array[ss] > 0) {

            immunity_days_of_array[ss]++;

            if ((immunity_level_of_array[ss] > immunity_level_until_loss_start)
                && (immunity_days_of_array[ss] > immunity_days_until_loss_start)
                ) {
                
                if (uniform_random_number_array[ss] < immunity_loss_rate) {

                    assert(immunity_level_of_array[ss] > 0);
                    immunity_level_of_array[ss]--;

                    if (immunity_level_of_array[ss] == 0){
                        immunity_days_of_array[ss] = 0;
                    }
                }

            }

        }
    }
}



float* BloodSystemManager::get_susceptibility() const {
    return this->susceptibility_of;
}
float* BloodSystemManager::get_infectiousness() const {
    return this->infectiousness_of ;
}
common::StageName* BloodSystemManager::get_most_advanced_stage() const {
    return this->most_advanced_stage_of;
}

common::DrugName* BloodSystemManager::get_drug_of() const {
    return this->drug_of;
}

common::DrugName* BloodSystemManager::get_drug_given_of() const {
    return this->drug_given_of;
}

int* BloodSystemManager::get_drug_given_day_of() const {
    return this->drug_given_day_of;
}

bool* BloodSystemManager::get_is_clinical_of() const {
    return this->is_clinical_of;
}

common::ParasiteType* BloodSystemManager::get_dominant_type_of() const {
    return this->dominant_type_of;
}

common::ParasiteType* BloodSystemManager::get_pg_type() const {
    return this->pg_type_of;
}

common::StageName* BloodSystemManager::get_pg_stage() const {
    return this->pg_stage_of;
}

uint8_t* BloodSystemManager::get_immunity_level_of() const {
    return this->immunity_level_of;
}

uint8_t* BloodSystemManager::get_cummulative_exposures_of() const{
    return this->cummulative_exposures_of;
}

uint64_t BloodSystemManager::get_parasite_count() const {
    uint64_t total=0;

    for(int pp = 0; pp < this->sum_num_pgs; pp++) {
        total += this->pg_num_parasites_of[pp];
    }

    return total;
}



float BloodSystemManager::get_mellow_rate() const {
    return this->kMellow_rate;
}
int BloodSystemManager::get_drug_failure_if_clinical_on_day() const {
    return this->kDrug_failure_if_clinical_on_day;
}


void BloodSystemManager::print_all() const{
    std::cout << "Table of Blood Systems, | #id drug'd'drug_days immunity_days'i'immunity_times head~tail dominant [parasite_group_type/stage/days:#p,#g]| ..." << std::endl;
    util::print_data_table(this, this->sum_num_systems, 6, 1);
}



void BloodSystemManager::print_one(const int index_bsys) const{
    std::cout << std::setw(8) <<  index_bsys ;
    std::cout << " " << this->drug_of[index_bsys] ;
    std::cout << " " << this->drug_given_of[index_bsys] ;
    std::cout << " " << this->drug_given_day_of[index_bsys] ;

    // std::cout << "d" << (int) this->drug_days_of[index_bsys] ;

    std::cout << " " << (int) this->immunity_days_of[index_bsys] ;
    std::cout << "i" << (int) this->cummulative_exposures_of[index_bsys] ;

    std::cout << " " << (int) this->head_pg_of[index_bsys] ;
    std::cout << "~" << (int) this->tail_pg_of[index_bsys] ;

    std::cout << " " << this->dominant_type_of[index_bsys] ;

    if (this->head_pg_of[index_bsys]>=0){    
        int lbound = index_bsys*kNum_pgs_per_system;
        for (int pp = lbound+this->head_pg_of[index_bsys]; pp <= lbound + this->tail_pg_of[index_bsys]; pp++){
            if (this->is_clinical_of[index_bsys]) {
                std::cout << " c[";
            } else {
                std::cout << "  [";
            }
            std::cout << this->pg_type_of[pp] << "/"
                      << this->pg_stage_of[pp] << "/"
                      << this->pg_aquired_on_day_of[pp] << ":"
                      << this->pg_num_parasites_of[pp] << ","
                      << this->pg_num_gametocytes_of[pp] << ", dt:"
                      << this->pg_is_drug_target_of[pp]<< ", cc:"
                      << this->pg_has_caused_clinicalness_of[pp]<< "]";
        }
    } else{
        std::cout << " [empty] ";
    }
}

void BloodSystemManager::print_drug_parameters() const {
    std::cout << "[B_Sys] Drug Parameters:\n";
    for(auto pp : util::Enum_iterator<common::ParasiteType>()) {
        
        if (pp == common::ParasiteType::kNoParasite) continue;
        
        for(auto dd : util::Enum_iterator<common::DrugName>()) {
            
            if (dd == common::DrugName::kNoDrug) continue;
            
            for(auto ss: util::Enum_iterator<common::StageName>()) {

                if (ss == common::StageName::kNotInSystem) continue;

                std::cout << "kPds["
                            << pp << "]["
                            << dd << "]["
                            << ss << "]="
                            << this->kPds_kill_probability[int(pp)][int(dd)][int(ss)]
                            << "\n";
            }
        }
    }
}
    
}