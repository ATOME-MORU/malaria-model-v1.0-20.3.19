#include <string>
#include <cassert>
#include <vector>
#include <iomanip>

#include <boost/chrono.hpp>
#include <boost/chrono/process_cpu_clocks.hpp>

// #include "rapidjson/document.h"

#include <string>
#include <fstream> //std::ofstream
// #include <thread>

#include <experimental/filesystem> // -lstdc++fs in LDLIBS
namespace fs = std::experimental::filesystem;

#include <boost/filesystem.hpp>


#include "rapidjson/document.h"

#include "util/util.h"
#include "util/randomness.h"

#include "village/village.h"
#include "village/mosquito.h"
#include "village/mosquito_ibm_ra.h"
#include "human/human.h"
#include "human/blood_system.h"

#include "intervention/mda.h"
#include "intervention/fmda.h"
#include "intervention/smda.h"
#include "intervention/tmda.h"
#include "intervention/vc.h"

#include "simulation/reporter.h"

#include "simulation/simulator.h"



namespace simulation{

    
// std::string kJson_schema_config_file_name = "schema.config.json";


Simulator::Simulator (
        
        const std::string& configuration_file_name,
        const std::string output_directory_name,
        const std::string& output_prefix_name

    ) : config(util::get_json_from_file(configuration_file_name)),
        // config_checked(this->check_config(config)),
        output_directory(output_directory_name),
        output_prefix(output_prefix_name.empty() ? util::get_output_prefix() : output_prefix_name),
        num_time_steps(config["simulation"]["total_steps"].GetInt()),
        kIf_batch_mode(config["simulation"]["batch_mode"].GetBool()),
        vll_mgr(
            kIf_batch_mode,
            num_time_steps,

            config["village"]["data_file"]["name"].GetString(),
            '\t',
            config["village"]["data_file"]["rows_to_read"].GetInt(),
            config["village"]["data_file"]["rows_to_skip"].GetInt(),
            config["village"]["data_file"]["overwrite_ra_rate_if"].GetBool(),
            config["village"]["data_file"]["overwrite_ra_rate_with"].GetFloat(),
            
            config["treatment"]["func_treatment_general_version"].GetString(),

            config["treatment"]["malaria_post"]["timecomp"].GetFloat(),
            config["treatment"]["malaria_post"]["fullcourse"].GetFloat(),
            config["treatment"]["malaria_post"]["covab"].GetFloat(),
            config["treatment"]["malaria_post"]["nomp"].GetFloat(),
            config["treatment"]["malaria_post"]["asymptomatic_to_clinical_ratio"].GetFloat(),
            config["treatment"]["malaria_post"]["establish_malaria_post_during_mda"].GetBool(),
            config["treatment"]["malaria_post"]["establish_malaria_post_during_survey"].GetBool(),
            config["treatment"]["malaria_post"]["establish_malaria_post_on_day"].GetInt(),

            config["treatment"]["treatment_rate_mda"].GetFloat(),
            
            config["treatment"]["mmc_wp2"]["coverage"].GetFloat(),
            config["treatment"]["mmc_wp2"]["drug"].GetInt(),

            config["infection"]["init_with_uniform_prevalence_if"].GetBool(),
            config["infection"]["init_with_uniform_prevalence_with"].GetFloat(),

            config["infection"]["init_transmission_coefficient_scaling_factor_a"].GetFloat(),
            config["infection"]["init_transmission_coefficient_scaling_factor_b"].GetFloat(),
            config["mosquito"]["init"]["infected_mosquito_to_infected_human_ratio"].GetFloat(),
            config["mosquito"]["init"]["infectious_mosquito_to_infected_human_ratio"].GetFloat(),
            config["infection"]["transmission_coefficient_to_beta_scaling_factor"].GetFloat(),

            config["mosquito"]["seasonality"]["switched_on"].GetBool(),

            config["mosquito"]["seasonality"]["use_data_file"].GetBool(),
            config["mosquito"]["seasonality"]["data_file"]["name"].GetString(),

            config["mosquito"]["seasonality"]["amplitude"].GetFloat(),
            config["mosquito"]["seasonality"]["amplitude_multiplier"].GetFloat(),

            config["mosquito"]["seasonality"]["cos_pi_multiplier"].GetFloat(),
            config["mosquito"]["seasonality"]["cos_days_offset"].GetFloat(),

            config["infection"]["func_m2h"]["version"].GetInt(),

            config["infection"]["susceptibility_multiplier"].GetFloat(),
            config["infection"]["infectiousness_multiplier"].GetFloat(),

            config["mosquito"]["biting_rate"].GetFloat(),

            config["mosquito"]["incubation_days"].GetFloat(),
            config["mosquito"]["infected_days"].GetFloat(),
            config["mosquito"]["infectious_days"].GetFloat(),

            config["mosquito"]["min_filter"]["threshold"].GetFloat(),
            config["mosquito"]["min_filter"]["start_on_day"].GetInt(),

            config["mobility"]["enabled"].GetBool(),
            config["mobility"]["static_population_move_out_probability"].GetFloat(),
            config["mobility"]["static_population_return_home_probability"].GetFloat(),
            config["mobility"]["non_static_population_move_out_probability"].GetFloat(),
            config["mobility"]["non_static_population_return_home_probability"].GetFloat()

            // config["mosquito"]["ode"]["k_value"].GetFloat()
        ),
        mda(
            config["intervention"]["mda"]["num_teams"].GetInt(),
            vll_mgr.sum_num_villages,

            config["intervention"]["mda"]["num_days_spent_in_village"].GetInt(),
            config["intervention"]["mda"]["num_drug_rounds"].GetInt(),
            config["intervention"]["mda"]["num_days_between_drug_rounds"].GetInt(),
            
            vll_mgr.get_village_distances()
        ),
        fmda(
            config["intervention"]["fmda"]["num_teams"].GetInt(),
            vll_mgr.sum_num_villages,

            config["intervention"]["fmda"]["num_days_spent_in_village"].GetInt(),
            config["intervention"]["fmda"]["num_drug_rounds"].GetInt(),
            config["intervention"]["fmda"]["num_days_between_drug_rounds"].GetInt(),

            vll_mgr.get_village_distances(),
            config["intervention"]["fmda"]["survey_data_delivery_delay"].GetInt(),
            config["intervention"]["fmda"]["pcr_sensitivity"].GetFloat(),
            config["intervention"]["fmda"]["pcr_sample_size"].GetInt(),
            config["intervention"]["fmda"]["prevalence_threshold"].GetFloat(),

            output_directory + util::kPathSeparator + output_prefix + "_fmda_",
            config["intervention"]["fmda"]["print_map"].GetBool()

        ),
        tmda(
            config["intervention"]["tmda"]["num_teams"].GetInt(),

            config["intervention"]["tmda"]["num_days_spent_in_village"].GetInt(),
            config["intervention"]["tmda"]["num_drug_rounds"].GetInt(),
            config["intervention"]["tmda"]["num_days_between_drug_rounds"].GetInt(),

            config["intervention"]["tmda"]["target_size"].GetFloat()
        ),
        vc(
            config["intervention"]["vector_control"]["efficacy"].GetFloat(),
            vll_mgr.sum_num_villages
        ),
        hmn_mgr(
            vll_mgr.sum_total_population,
            config["human"]["func_death_probability"]["version"].GetInt(),
            config["human"]["attractiveness"]["enabled"].GetBool(),
            config["human"]["attractiveness"]["mean"].GetFloat(),
            config["human"]["attractiveness"]["sd"].GetFloat(),
            config["human"]["attractiveness"]["reset_log_bitten_times_at_the_end_of_year"].GetFloat()
        ),
        bsys_mgr(
            vll_mgr.sum_total_population,

            config["intervention"]["itn"]["preference_indoor"].GetFloat(),
            config["intervention"]["itn"]["preference_human"].GetFloat(),
            config["intervention"]["itn"]["effectiveness"].GetFloat(),

            config["treatment"]["treat_clean_hosts"].GetBool(),
            config["treatment"]["treat_asymptomatic_hosts"].GetBool(),

            config["within_host"]["init_clinical_probability"].GetFloat(),

            config["within_host"]["susceptibility"]["default"].GetFloat(),
            config["within_host"]["susceptibility"]["func_update_susceptibility"]["version"].GetFloat(),

            config["within_host"]["infectiousness"]["default"].GetFloat(),
            
            config["within_host"]["infectiousness"]["phi_a"].GetFloat(),
            config["within_host"]["infectiousness"]["phi_c"].GetFloat(),

            config["within_host"]["prophylaxis_effect"].GetFloat(),

            config["within_host"]["func_dominance"]["version"].GetInt(),
            

            config["within_host"]["func_clinical_probability_blood_to_infectious"]["version"].GetInt(),
            config["within_host"]["func_clinical_probability_blood_to_infectious"]["version3_para1"].GetFloat(),
            config["within_host"]["func_clinical_probability_blood_to_infectious"]["version3_para2"].GetFloat(),

            config["within_host"]["mellow_days"].GetFloat(),

            config["within_host"]["stage_progression"]["l2b_days_mean"].GetFloat(),
            config["within_host"]["stage_progression"]["l2b_days_sd"].GetFloat(),
            config["within_host"]["stage_progression"]["b2i_days_mean"].GetFloat(),
            config["within_host"]["stage_progression"]["b2i_days_sd"].GetFloat(),
            config["within_host"]["stage_progression"]["i2r_days_mean"].GetFloat(),
            config["within_host"]["stage_progression"]["i2r_days_sd"].GetFloat(),

            // config["within_host"]["natural_recovery_days"].GetFloat(),
            config["within_host"]["natural_recovery_resistance_cost"].GetFloat(),

            config["drug"]["failure_if_clinical_on_day"].GetInt(),
            config["drug"]["allow_act_to_a"].GetBool(),
            config["drug"]["possible_to_have_drug_a_only"].GetBool(),
            // config["drug"]["drug_b_resistance"].GetFloat(),
            // config["drug"]["drug_ab_resistance"].GetFloat(),

            config["drug"]["base_effect_multiplier"].GetFloat(),

            config["drug"]["resistance_multipliers"]["enabled"].GetBool(),
            config["drug"]["resistance_multipliers"]["drug_b_effect_on_rb"].GetFloat(),
            config["drug"]["resistance_multipliers"]["drug_ab_effect_on_rb"].GetFloat(),
            config["drug"]["resistance_multipliers"]["drug_a_effect_on_ra"].GetFloat(),
            config["drug"]["resistance_multipliers"]["drug_ab_effect_on_ra"].GetFloat(),
            config["drug"]["resistance_multipliers"]["drug_a_effect_on_rab"].GetFloat(),
            config["drug"]["resistance_multipliers"]["drug_b_effect_on_rab"].GetFloat(),
            config["drug"]["resistance_multipliers"]["drug_ab_effect_on_rab"].GetFloat(),

            config["genotype"]["mutation"]["prob_single_resistance"].GetFloat(),
            config["genotype"]["mutation"]["prob_single_resistance_w_act"].GetFloat(),
            config["genotype"]["mutation"]["prob_single_type_a_type_b_ratio"].GetFloat(),
            config["genotype"]["mutation"]["type_a_enabled"].GetBool(),
            config["genotype"]["mutation"]["type_b_enabled"].GetBool(),
            // config["genotype"]["mutation_prob_double_resistance"].GetFloat(),

            config["genotype"]["mutation"]["func_parasite_mutation"]["version"].GetInt(),
            config["genotype"]["mutation"]["func_parasite_mutation"]["fixation_after_day"].GetInt(),


            config["within_host"]["immunity"]["days_until_loss_start"].GetInt(),
            config["within_host"]["immunity"]["level_until_loss_start"].GetInt(),
            config["within_host"]["immunity"]["loss_days"].GetFloat()
        ),
        rpt_mgr(
            kIf_batch_mode,
            config["simulation"]["title"].GetString(),
            output_directory + util::kPathSeparator + output_prefix,
            num_time_steps,

            config["simulation"]["reports"]["prevalence_daily"].GetBool(),
            config["simulation"]["reports"]["drug_daily"].GetBool(),
            config["simulation"]["reports"]["infection_daily"].GetBool(),
            config["simulation"]["reports"]["infection_daily_village"].GetBool(),
            config["simulation"]["reports"]["mosquito_daily"].GetBool(),
            config["simulation"]["reports"]["treatment_daily"].GetBool(),
            config["simulation"]["reports"]["mutation_daily"].GetBool(),
            config["simulation"]["reports"]["clearance_natural_daily"].GetBool(),
            config["simulation"]["reports"]["clearance_drug_daily"].GetBool(),
            config["simulation"]["reports"]["village_prev_annual"].GetBool(),
            config["simulation"]["reports"]["village_infections_annual"].GetBool(),


            config["simulation"]["plots"]["prevalence_daily_y_max"].GetInt(),

            &vll_mgr,
            &hmn_mgr,
            &bsys_mgr
        ) {

        std::cout << "got here" << std::endl;

        boost::filesystem::copy_file(
            configuration_file_name,
            this->output_directory
                + util::kPathSeparator
                + this->output_prefix
                + "_config.json"
        );

        assert(this->vll_mgr.sum_total_population == this->hmn_mgr.sum_num_humans);

        //////////////////////////////////////////////////////////////////////
        // Villages
        this->vll_mgr.init_treatment_rates_malaria_post();
        //this->vll_mgr.print_summary();


        //////////////////////////////////////////////////////////////////////
        // MDA

        if (this->config["intervention"]["mda"]["enabled"].GetBool()) {
            assert(this->config["intervention"]["mda"]["schedule"].IsArray());
            for(rapidjson::SizeType ss = 0; ss < config["intervention"]["mda"]["schedule"].Size(); ss++){
                assert(this->config["intervention"]["mda"]["schedule"][ss].IsInt());
                this->mda.register_start_time(this->config["intervention"]["mda"]["schedule"][ss].GetInt());
            }
            
            // this->output_mda();

            if (this->config["intervention"]["mda"]["start_teams_at_highest_tc"].GetBool()) {
                this->mda.reroute_teams(this->vll_mgr.get_village_ids_ordered_by_tc_desc());
                // this->output_mda();
            }

            if (this->config["intervention"]["mda"]["start_teams_at_lowest_tc"].GetBool()) {
                this->mda.reroute_teams(this->vll_mgr.get_village_ids_ordered_by_tc_ascd());
                // this->output_mda();
            }
            
        }

        // std::cin.ignore();


        //////////////////////////////////////////////////////////////////////
        // fMDA
        if (this->config["intervention"]["fmda"]["enabled"].GetBool()) {
            assert(this->config["intervention"]["fmda"]["schedule"].IsArray());
            for(rapidjson::SizeType ss = 0; ss < config["intervention"]["fmda"]["schedule"].Size(); ss++){
                assert(this->config["intervention"]["fmda"]["schedule"][ss].IsInt());
                this->fmda.register_start_time(
                    this->config["intervention"]["fmda"]["schedule"][ss].GetInt(),
                    this->config["intervention"]["fmda"]["survey_schedule"]["enabled"].GetBool()
                );
            }

            if (this->config["intervention"]["fmda"]["survey_schedule"]["enabled"].GetBool()) {
                for(rapidjson::SizeType ss = 0; ss < config["intervention"]["fmda"]["survey_schedule"]["mean"].Size(); ss++){
                    this->fmda.register_survey_time(
                        this->config["intervention"]["fmda"]["survey_schedule"]["mean"][ss].GetInt(),
                        this->config["intervention"]["fmda"]["survey_schedule"]["sd"][ss].GetInt(),
                        this->config["intervention"]["fmda"]["survey_schedule"]["delivery_by"][ss].GetInt()
                    );
                }
            }
            // this->mda.print_all_routes();
        }

        //////////////////////////////////////////////////////////////////////
        // tMDA
        if (this->config["intervention"]["tmda"]["enabled"].GetBool()) {
            assert(this->config["intervention"]["tmda"]["schedule"].IsArray());
            for (rapidjson::SizeType ss = 0; ss < this->config["intervention"]["tmda"]["schedule"].Size(); ss++) {
                assert(this->config["intervention"]["tmda"]["schedule"][ss].IsInt());
                this->tmda.register_start_time(this->config["intervention"]["tmda"]["schedule"][ss].GetInt());
            }
        }

        //////////////////////////////////////////////////////////////////////
        // Vector Control
        if (this->config["intervention"]["vector_control"]["enabled"].GetBool()) {
            assert(this->config["intervention"]["vector_control"]["schedule"]["starts"].IsArray());
            assert(this->config["intervention"]["vector_control"]["schedule"]["durations"].IsArray());
            assert(this->config["intervention"]["vector_control"]["schedule"]["starts"].Size()
                    == this->config["intervention"]["vector_control"]["schedule"]["durations"].Size());
            for (rapidjson::SizeType ss = 0; ss < this->config["intervention"]["vector_control"]["schedule"]["starts"].Size(); ss++) {
                this->vc.init_add_to_schedule(
                    this->config["intervention"]["vector_control"]["schedule"]["starts"][ss].GetInt(),
                    this->config["intervention"]["vector_control"]["schedule"]["durations"][ss].GetInt()
                );
            }
        }

        //////////////////////////////////////////////////////////////////////
        // Humans
        assert(this->config["human"]["age_weight_vector_file"].IsString());
        this->hmn_mgr.read_age_distribution(this->config["human"]["age_weight_vector_file"].GetString());
        if (config["human"]["func_death_probability"]["version"].GetInt() == 3) {
            this->hmn_mgr.read_death_prob_distribution(this->config["human"]["death_prob_vector_file"].GetString());
        }
        this->hmn_mgr.init_age();
        // this->hmn_mgr.print_age_distribution();

        assert(this->config["human"]["male_percentage"].IsNumber());
        // assert(this->config["human"]["male_percentage"].IsFloat());
        this->hmn_mgr.init_gender(this->config["human"]["male_percentage"].GetFloat());
        // this->hmn_mgr.print_gender_distribution();



        assert(this->config["intervention"]["itn"]["init_probability_per_human"].IsNumber());
        // assert(this->config["human"]["itn_probability"].IsFloat());
        this->hmn_mgr.init_itn(this->config["intervention"]["itn"]["init_probability_per_human"].GetFloat());


        //////////////////////////////////////////////////////////////////////
        // BSys
        if (this->config["simulation"]["init_infections_from_agestruct"].GetBool()) {

            
            if (this->config["within_host"]["init_file"]["use_old_format"].GetBool()) {

                assert(this->config["within_host"]["init_file"]["name_old_format"].IsString());

                //// old agestruct*.txt
                this->bsys_mgr.init_from_file(
                    this->config["within_host"]["init_file"]["name_old_format"].GetString(),
                    '\t',
                    &(this->hmn_mgr),
                    this->config["within_host"]["init_file"]["use_age_column"].GetBool()
                );

            } else {

                assert(this->config["within_host"]["init_file"]["name"].IsString());

                //// new agestruct*.csv from export
                this->bsys_mgr.import_from_file(
                    this->config["within_host"]["init_file"]["name"].GetString(),
                    &(this->hmn_mgr),
                    this->config["within_host"]["init_file"]["use_age_column"].GetBool()
                );

                if (this->config["within_host"]["init_file"]["ra_rate"].GetFloat() > 0.0) {

                    this->bsys_mgr.step_report_genotype_weights();
                    this->bsys_mgr.parasite_injection(
                        common::ParasiteType::kRa,
                        this->bsys_mgr.daily_record_of_parasite_positive_systems
                            * this->config["within_host"]["init_file"]["ra_rate"].GetFloat()
                    );
                    this->bsys_mgr.parasite_dominance();
                
                }
                
            }

            if (this->config["within_host"]["init_file"]["use_age_column"].GetBool()) {
                this->hmn_mgr.update_all_death_probability();
            }

        }

        this->bsys_mgr.print_drug_parameters();
        // std::cin.ignore();

        // this->bsys_mgr.update_susceptibility(
        //     this->hmn_mgr.get_if_itn_of(),
        //     this->hmn_mgr.sum_num_humans
        // );


        //////////////////////////////////////////////////////////////////////
        // Humans <--> Villages
        //////////////////////////////////////////////////////////////////////
        this->vll_mgr.init_get_villagers(&(this->hmn_mgr));
        this->hmn_mgr.init_mobility(this->vll_mgr.get_migrant_population_percentages());
        //this->hmn_mgr.print_all();

        //////////////////////////////////////////////////////////////////////
        // Bsys <--> Human 
        //////////////////////////////////////////////////////////////////////
        this->bsys_mgr.update_susceptibility(
            this->hmn_mgr.get_if_itn_of(),
            this->hmn_mgr.sum_num_humans,
            this->hmn_mgr.get_age_of()
        );


        //////////////////////////////////////////////////////////////////////
        // Bsys <--> Human <--> Villages
        //////////////////////////////////////////////////////////////////////

        if (this->config["simulation"]["init_infections_from_regions"].GetBool()){

            this->vll_mgr.init_infections(
                this->hmn_mgr.get_bitten_by(),
                this->hmn_mgr.sum_num_humans
            );
            this->bsys_mgr.parasite_intake(
                this->hmn_mgr.get_bitten_by(),
                common::StageName::kInfectious,
                0,
                true
            );
            this->bsys_mgr.parasite_dominance();

        }

        //this->bsys_mgr.print_all();


    //////////////////////////////////////////////////////////////////////
    // Mosquito Manager
        
        if (config["mosquito"]["ibm"]["in_use"].GetBool()) {
            this->msq_mgr_ibm = new village::MosquitoManagerIbmRa(
                config["mosquito"]["aquatic_mortality"].GetFloat(),
                config["mosquito"]["maturation_days"].GetFloat(),
                config["mosquito"]["immat_life_days"].GetFloat(),
                config["mosquito"]["adult_life_days"].GetFloat(),

                config["mosquito"]["laying_days"].GetFloat(),
                config["mosquito"]["biting_rate"].GetFloat(),
                config["mosquito"]["incubation_days"].GetFloat(),
                config["mosquito"]["eggs_per_lay"].GetFloat(),

                config["mosquito"]["ibm"]["max_larval_to_max_adult_ratio"].GetFloat(),
                config["mosquito"]["ibm"]["scaletomosqnumbers"].GetFloat(),
                config["mosquito"]["ibm"]["func_swamp_step_inc_born_death"]["version"].GetInt(),
                vll_mgr.sum_num_villages,
                vll_mgr.get_population_of(),
                vll_mgr.get_transmission_coefficient_of()
            );

        } else if (config["mosquito"]["ode"]["in_use"].GetBool()) {
            this->msq_mgr = new village::MosquitoManager(
                vll_mgr.get_population_of(),
                config["mosquito"]["aquatic_mortality"].GetFloat(),
                config["mosquito"]["maturation_days"].GetFloat(),
                config["mosquito"]["immat_life_days"].GetFloat(),
                config["mosquito"]["adult_life_days"].GetFloat(),
                config["mosquito"]["laying_days"].GetFloat(),
                config["mosquito"]["eggs_per_lay"].GetFloat(),
                config["mosquito"]["incubation_days"].GetFloat(),
                config["mosquito"]["biting_rate"].GetDouble(),
                config["mosquito"]["ode"]["carrying_capacity_multiplier"].GetDouble(),
                vll_mgr.sum_num_villages,
                config["mosquito"]["ode"]["init"]["max_larval_capacity_file"]["name"].GetString()
            );
            
            this->bsys_mgr.parasite_dominance();
            this->vll_mgr.update_prevalence_count_of(this->bsys_mgr.get_dominant_type_of());
            
            this->msq_mgr->set_verbose_on();
            // assert(msq_mgr->init_step_until_equilibrium(
            //     vll_mgr.get_sum_prevalence_infected_of(
            //         this->bsys_mgr.get_dominant_type_of(),
            //         this->bsys_mgr.get_infectiousness()
            //     )
            // ));
            this->msq_mgr->set_verbose_off();
            // std::cin.ignore();
        } else {
            this->msq_mgr = new village::MosquitoManager(
                vll_mgr.get_population_of(),
                config["mosquito"]["aquatic_mortality"].GetFloat(),
                config["mosquito"]["maturation_days"].GetFloat(),
                config["mosquito"]["immat_life_days"].GetFloat(),
                config["mosquito"]["adult_life_days"].GetFloat(),
                config["mosquito"]["laying_days"].GetFloat(),
                config["mosquito"]["eggs_per_lay"].GetFloat(),
                config["mosquito"]["incubation_days"].GetFloat(),
                config["mosquito"]["biting_rate"].GetDouble(),
                0.0,
                std::vector<double>(vll_mgr.sum_num_villages, 0.0)
            );
        }

        //////////////////////////////////////////////////////////////////////
        // Reporter
        this->rpt_mgr.open_ofstreams();


}

Simulator::~Simulator() {
    if (config["mosquito"]["ibm"]["in_use"].GetBool()) {
        delete this->msq_mgr_ibm;
    } else {
        delete this->msq_mgr;
    }
}

bool Simulator::validate_config(const std::string& configuration_file_name) {
    
    bool check_pass = false;
    rapidjson::Document config_doc = util::get_json_from_file(configuration_file_name);

    //////////////////////////////////////////////////////////////////////
    // JSON Validation
    //////////////////////////////////////////////////////////////////////
    assert(config_doc["simulation"]["schema_file_name"].IsString());
    rapidjson::Document rj_doc_config_schema =
        util::get_json_from_file(config_doc["simulation"]["schema_file_name"].GetString());

    assert(util::validate_json_against_schema(&rj_doc_config_schema, &config_doc));


    assert(config_doc["simulation"]["total_steps"].IsInt());


    assert(config_doc["village"]["data_file"]["name"].IsString());
    assert(config_doc["village"]["data_file"]["rows_to_read"].IsInt());

    assert(config_doc["treatment"]["malaria_post"]["timecomp"].IsNumber());
    assert(config_doc["treatment"]["malaria_post"]["fullcourse"].IsNumber());
    assert(config_doc["treatment"]["malaria_post"]["covab"].IsNumber());
    assert(config_doc["treatment"]["malaria_post"]["nomp"].IsNumber());
    assert(config_doc["treatment"]["malaria_post"]["asymptomatic_to_clinical_ratio"].IsNumber());

    // assert(config_doc["treatment"]["timecomp"].IsFloat());
    // assert(config_doc["treatment"]["fullcourse"].IsFloat());
    // assert(config_doc["treatment"]["covab"].IsFloat());
    // assert(config_doc["treatment"]["nomp"].IsFloat());
    // assert(config_doc["treatment"]["asymptomatic_to_clinical_ratio"].IsFloat());


    assert(config_doc["intervention"]["mda"]["num_teams"].IsInt());

    check_pass = true;
    return check_pass;

}

void Simulator::run(){

    // Humans <-> BloodSystems 
    bool* birth_death_reset_flags = new bool[hmn_mgr.sum_num_humans];
    std::fill_n(birth_death_reset_flags, hmn_mgr.sum_num_humans, false);

    int random_numbers_humans_size = hmn_mgr.sum_num_humans;
    float* random_numbers_humans = new float[random_numbers_humans_size];

    float* random_numbers_movement = new float[hmn_mgr.sum_num_humans*2];

    int random_numbers_drug_loss_size = bsys_mgr.sum_num_systems;
    float* random_numbers_drug_loss = new float[random_numbers_drug_loss_size];

    int random_number_repository_sys2_size = bsys_mgr.sum_num_systems*2;
    float* random_number_repository_sys2 = new float[random_number_repository_sys2_size];


    int performance_timing_start = this->config["simulation"]["performance"]["timing_start"].GetInt();
    int performance_timing_finish = this->config["simulation"]["performance"]["timing_finish"].GetInt();
    double performance_start_wall_time = 0.0;
    double performance_time_total = 0.0;

    boost::chrono::process_user_cpu_clock timer;

    boost::chrono::process_user_cpu_clock::time_point performance_start_cpu_time;
    boost::chrono::process_user_cpu_clock::duration performance_cpu_time_duration;

    // boost::chrono::process_user_cpu_clock::time_point step_start_cpu_time;

    // boost::chrono::process_cpu_clock performance_timer;
    // boost::chrono::process_cpu_clock::time_point performance_timer_start;
    // boost::chrono::process_cpu_clock::duration performance_timer_duration;

    bool kMosq_method_using_ode = config["mosquito"]["ode"]["in_use"].GetBool();
    bool kMosq_method_using_ibm = config["mosquito"]["ibm"]["in_use"].GetBool();

    // bool if_threading = false;
    // std::thread thrd_rnd_nmbrs_drug_loss;

    std::vector<double> stepper_time_cost;

    std::vector<double> time_cost_per_time_step;
    double wall_time_at_time_step_begin = 0.0;
    
    int death_total = 0;

    double wall_time_begin = util::get_wall_time();

    std::ofstream temp_file;

    // Run one end of day processes before simulation starts
    this->hmn_mgr.step_end_of_day(-1);

    for(int tt = 0; tt <= this->num_time_steps; tt++) {

        wall_time_at_time_step_begin = util::get_wall_time();

        // step_start_cpu_time = timer.now();
        // std::cout << step_start_cpu_time << std::endl;

        std::cout << tt << ".begin"  << "\n";

        if (tt == performance_timing_start) {
            performance_start_wall_time = util::get_wall_time();

            performance_start_cpu_time = timer.now();

            // performance_timer_start = performance_timer.now();
            // std::cout << "performance_start_wall_time: " << performance_start_wall_time << "\n";
            // std::cin.ignore();
        }

//SOD   ///// Mandatory start of day processes

            // Village
            this->vll_mgr.step_reset_updated_indicators();
            
            this->vll_mgr.step_update_beta_season_of(tt);

        ////// Ageing, Birth/Death
        // std::cout << "[1] Ageing, Birth/Death ...\n";
            // if ( (tt+1) % common::kNum_days_in_one_year == 0 ){
            if ( (tt>1) && (((tt-1) % common::kNum_days_in_one_year)==0) ) { //increase all ages on first day of year
                this->hmn_mgr.increase_all_ages_by_one();
                std::cout << "Day " << tt << ": all ages + 1\n";
            }

            util::set_random_numbers_uniform(
                random_numbers_humans,
                random_numbers_humans_size
            ); 

            std::fill_n(birth_death_reset_flags, hmn_mgr.sum_num_humans, false);
            int num_birth_and_death = this->hmn_mgr.birth_and_death(
                                                birth_death_reset_flags,
                                                random_numbers_humans,
                                                random_numbers_humans_size
                                        );

            this->bsys_mgr.birth_and_death(
                            birth_death_reset_flags,
                            this->hmn_mgr.sum_num_humans
            );

            death_total += num_birth_and_death;
            std::cout << "\tnum_birth_and_death=" << num_birth_and_death << "\n";

        ////// Interventions (Static + Dynamic)
        // std::cout << "[2] Interventions (Static + Dynamic) ...\n";

            if (this->config["intervention"]["mda"]["enabled"].GetBool()) {
                this->mda.step(
                    tt,
                    &(this->vll_mgr)
                );
            }

            if (this->config["intervention"]["fmda"]["enabled"].GetBool()) {

                this->vll_mgr.update_prevalence_count_of(this->bsys_mgr.get_dominant_type_of());
                this->vll_mgr.update_prevalence_of();

                this->fmda.step(
                    tt,
                    &(this->vll_mgr)
                    // this->bsys_mgr.get_dominant_type_of()
                    // this->output_directory + util::kPathSeparator
                    // + this->output_prefix + "_fmda_" + std::to_string(tt) + "_"
                );
            }

            if (this->config["intervention"]["tmda"]["enabled"].GetBool()) {
                if (this->tmda.get_next_start_time() == tt) {

                    this->vll_mgr.update_prevalence_count_of(this->bsys_mgr.get_dominant_type_of());
                    this->vll_mgr.update_prevalence_of();

                    this->tmda.create_mda(
                        tt,
                        this->vll_mgr.get_prevalence_of(),
                        this->vll_mgr.get_village_distances()
                    );

                    if (this->config["intervention"]["tmda"]["print_mda_map"].GetBool()) {
                        this->tmda.print_mda_map(
                            this->vll_mgr.sum_num_villages,
                            this->vll_mgr.longitude_of,
                            this->vll_mgr.latitude_of,
                            this->vll_mgr.population_of,
                            this->output_directory + util::kPathSeparator
                            + this->output_prefix + "_tmda_" + std::to_string(tt) + "_"
                        );
                    }

                }

                this->tmda.step(
                    tt,
                    &(this->vll_mgr)
                );
            }

            if (this->config["intervention"]["vector_control"]["enabled"].GetBool()) {
                // this->vll_mgr.update_prevalence_count_of(this->bsys_mgr.get_dominant_type_of());
                if (this->vc.get_next_start_time() == tt) {
                    
                    this->vll_mgr.update_prevalence_count_of(this->bsys_mgr.get_dominant_type_of());
                    this->vll_mgr.update_prevalence_of();
                    this->vc.set_target_on_top_villages(
                        this->vll_mgr.get_prevalence_of(),
                        this->config["intervention"]["vector_control"]["target_size"].GetFloat()
                    );
                    
                    std::cout << "day" << tt << ":vc_started\n";

                }

                this->vc.step_reduce_beta_season(
                    this->vll_mgr.get_beta_season_of(),
                    tt
                );
            }

            this->vll_mgr.step_treatment_general(
                &(this->bsys_mgr),
                &(this->hmn_mgr),
                tt
            );


            // this->vll_mgr.open_malaria_post(tt);
            // this->vll_mgr.treatment_by_malaria_post(
            //     &(this->bsys_mgr),
            //     &(this->hmn_mgr)
            // );

        ////// Mobility
        // std::cout << "[3] Mobility ...\n";

            if (this->config["mobility"]["enabled"].GetBool()) {

                util::set_random_numbers_uniform(random_numbers_movement, this->hmn_mgr.sum_num_humans*2);

            }
                this->vll_mgr.step_population_movement_static(
                    this->hmn_mgr.sum_num_humans,
                    this->hmn_mgr.get_mobility_of(),
                    this->hmn_mgr.get_home_village_of(),

                    random_numbers_movement,
                    this->hmn_mgr.sum_num_humans*2
                );
            
            // }

        ////// Immunity Loss
        // std::cout << "[4] Immunity Loss ...\n";
            util::set_random_numbers_uniform(
                random_number_repository_sys2, 
                random_number_repository_sys2_size/2
            );
            this->bsys_mgr.step_immunity_loss(
                random_number_repository_sys2,
                random_number_repository_sys2_size/2
            );
        
        ////// Drug Intake + Effect + C-resolution
        // std::cout << "[5]\n";

            // hmn_mgr.print_distribution<common::DrugName>(hmn_mgr.get_given_drug(), hmn_mgr.sum_num_humans);
            this->bsys_mgr.drug_intake(this->hmn_mgr.get_given_drug(), tt);
            this->bsys_mgr.step_report_drug_composition();
        // std::cout << "[P.1-3] Parasite Density (Replication/Decay/PD) ...\n";
            this->bsys_mgr.drug_effect(tt);

            this->bsys_mgr.clinical_resolution();


        ////// Parasites Stage Progression + Dominance
        // std::cout << "[6]\n";

            if (config["genotype"]["mutation"]["enabled"].GetBool() && tt >= config["genotype"]["mutation"]["start_at_time_step"].GetInt()) {
                this->bsys_mgr.parasite_mutation(tt);
                // this->bsys_mgr.parasite_mutation_mmc_wp2();
                // this->bsys_mgr.parasite_mutation_mmc_wp2_drug_trigger();
                
            }

            this->bsys_mgr.parasite_stage_progression(tt); // includes natural recovery 
            this->bsys_mgr.parasite_dominance(); // prevalence and clinical stats done here

            this->bsys_mgr.step_report_drug_failure_by_clinical_status_on_day(tt);


        ////// Drug Loss (PK)
        // std::cout << "[7]\n";
            // if (tt == 0 && if_threading){
            //     // util::set_random_numbers_uniform(random_numbers_drug_loss, random_numbers_drug_loss_size);
            //     thrd_rnd_nmbrs_drug_loss = std::thread(
            //         util::set_random_numbers_uniform_boost_thread, 
            //         random_numbers_drug_loss,
            //         random_numbers_drug_loss_size
            //     );
            // }
            
            // if (if_threading) {
            //     thrd_rnd_nmbrs_drug_loss.join();
            // } else {
                util::set_random_numbers_uniform(
                    random_numbers_drug_loss,
                    random_numbers_drug_loss_size
                );
            // }

            this->bsys_mgr.drug_loss(
                random_numbers_drug_loss,
                random_numbers_drug_loss_size
            );

            // if (if_threading) {
            //     thrd_rnd_nmbrs_drug_loss = std::thread(
            //         util::set_random_numbers_uniform_boost_thread, 
            //         random_numbers_drug_loss,
            //         this->bsys_mgr.sum_num_systems * 2
            //     );
            // }

        ////// H2M Infection & Mosquito Dynamics
        // std::cout << "[8]\n" ;

            this->bsys_mgr.update_susceptibility(
                this->hmn_mgr.get_if_itn_of(),
                this->hmn_mgr.sum_num_humans,
                this->hmn_mgr.get_age_of()
            );

            if (kMosq_method_using_ibm) {

                double before_step = util::get_wall_time();

                // this->msq_mgr_ibm->step(
                this->msq_mgr_ibm->step_matching_ode(
                    tt,
                    this->vll_mgr,

                    this->bsys_mgr.get_susceptibility(),
                    config["infection"]["susceptibility_multiplier"].GetFloat(),
                    this->bsys_mgr.get_most_advanced_stage(),
                    this->hmn_mgr.get_bitten_by(),
                    this->hmn_mgr.sum_num_humans,

                    this->bsys_mgr.get_infectiousness(),
                    config["infection"]["infectiousness_multiplier"].GetFloat(),
                    this->bsys_mgr.get_dominant_type_of(),
                    this->bsys_mgr.sum_num_systems,
                    this->hmn_mgr.get_home_village_of()
                );

                std::cout << "#IBM mosq si requests: " << this->msq_mgr_ibm->sum_num_superinfection_requests << "\n";

                stepper_time_cost.push_back(util::get_wall_time()-before_step);


                this->vll_mgr.daily_record_of_mosquito_infected_sum = this->msq_mgr_ibm->sum_num_infected_biting;
                this->vll_mgr.daily_record_of_mosquito_infectious_sum = this->msq_mgr_ibm->sum_num_infectious_biting;
                this->vll_mgr.daily_record_of_mosquito_clean_sum = this->msq_mgr_ibm->sum_num_clean_biting;
                this->vll_mgr.daily_record_of_mosquito_new_adults_sum = this->msq_mgr_ibm->sum_num_new_adults;
                this->vll_mgr.daily_record_of_mosquito_death_adults_sum = this->msq_mgr_ibm->sum_num_death_adults;
                this->vll_mgr.daily_record_of_mosquito_death_adults_infectious_count_sum = this->msq_mgr_ibm->sum_num_death_adults_infectious_count;


                this->vll_mgr.daily_record_of_mosquito_newly_infected = this->msq_mgr_ibm->sum_num_new_infected;
                this->vll_mgr.daily_record_of_mosquito_death_infected = this->msq_mgr_ibm->sum_num_death_infected;
                this->vll_mgr.daily_record_of_mosquito_death_infectious = this->msq_mgr_ibm->sum_num_death_infectious;
                this->vll_mgr.daily_record_of_mosquito_incubation = this->msq_mgr_ibm->sum_num_incubation;

                this->vll_mgr.daily_record_of_eggs_sum = static_cast<int>(this->msq_mgr_ibm->sum_num_eggs);
                this->vll_mgr.daily_record_of_larv_sum = static_cast<int>(this->msq_mgr_ibm->sum_num_larv);
                this->vll_mgr.daily_record_of_imat_sum = static_cast<int>(this->msq_mgr_ibm->sum_num_imat);
                this->vll_mgr.daily_record_of_aliv_sum = static_cast<int>(this->msq_mgr_ibm->sum_num_malive);
                this->vll_mgr.daily_record_of_aliv_biting_sum = static_cast<int>(this->msq_mgr_ibm->sum_num_malive_biting);

                this->vll_mgr.daily_record_of_max_capacity_sum = static_cast<int>(this->msq_mgr_ibm->sum_num_larva_capacity);

                // std::cout << "IBM/sum_num_inc_biting = " << this->msq_mgr_ibm->sum_num_infected_biting;
                // std::cout << "\nIBM/sum_num_inf_biting = " << this->msq_mgr_ibm->sum_num_infectious_biting;
                // std::cout << "\n";

                this->vll_mgr.daily_record_of_new_bites_sum = this->msq_mgr_ibm->sum_num_infectious_bites;
                this->vll_mgr.daily_record_of_new_infections_sum = this->msq_mgr_ibm->sum_num_successful_infectious_bites;

                this->bsys_mgr.parasite_intake(
                    this->hmn_mgr.get_bitten_by(),
                    common::StageName::kLiver,
                    // common::StageName::kBlood,
                    tt,
                    false
                );

// if (tt == (730+180)) {


//     std::ofstream mosq_pop_file_ibm;
//     mosq_pop_file_ibm.open ("temp_ibm_mosq_pop_730.csv");
//         // float sum = 0;
//     for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++) {
//         // this->msq_mgr->set_mosq_time(vv, this->msq_mgr->get_mosq_time(vv));
//         // if ( vv > 500) {
//             mosq_pop_file_ibm << vv << ": ";

//             // mosq_pop_file_ibm << v_temp_mosq.at(vv);
//             std::vector<int> mosq_pop = this->msq_mgr_ibm->get_mosq_pop_of_swamp(vv);

//             for (int pop : mosq_pop) {
//                 mosq_pop_file_ibm << " " << pop;

//             }
//             mosq_pop_file_ibm << "\n";
//         // }
//         // std::cout << ii << ":" << sum<< "\n";
//     }
//     mosq_pop_file_ibm.close();
//     std::cin.ignore();
// }



            } else if (kMosq_method_using_ode) {


                // this->bsys_mgr.parasite_dominance();



// if (tt == (7300 +180)) {
// // if (tt == 200) {
//     village::MosquitoManager* msq_mgr_copy = this->msq_mgr;

//     std::vector<float> v_temp_mosq = vll_mgr.get_sum_prevalence_infected_of(
//             this->bsys_mgr.get_dominant_type_of(),
//             this->bsys_mgr.get_infectiousness()
//         );

//             this->msq_mgr = new village::MosquitoManager(
//                 // config["mosquito"]["ode"]["feed_rate"].GetDouble(),
//                 config["mosquito"]["biting_rate"].GetDouble(),
//                 config["mosquito"]["ode"]["carrying_capacity_multiplier"].GetDouble(),
//                 vll_mgr.sum_num_villages,
//                 config["mosquito"]["ode"]["init"]["max_larval_capacity_file"]["name"].GetString()
//             );

//                     std::ofstream mosq_pop_file;
//                     mosq_pop_file.open ("mosq_pop_7300_180.csv");
//                         // float sum = 0;
//                     for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++) {
//                         this->msq_mgr->set_mosq_time(vv, msq_mgr_copy->get_mosq_time(vv));
//                         // if ( vv > 500) {
//                             mosq_pop_file << vv << ": ";
//                             mosq_pop_file << v_temp_mosq.at(vv);

//                             for (int ii = 0; ii < 9; ii++) {
//                                 // sum += this->msq_mgr->get_mosq_pop(vv, ii);
//                                 double mosq_num_temp = msq_mgr_copy->get_mosq_pop(vv,ii);
//                                 if (mosq_num_temp < 0.00000001) {
//                                     mosq_num_temp = 0.0;
//                                 }
//                                 this->msq_mgr->set_mosq_pop(vv,ii,mosq_num_temp);
//                                 mosq_pop_file << " " << msq_mgr_copy->get_mosq_pop(vv,ii);

//                             }
//                             mosq_pop_file << "\n";
//                         // }
//                         // std::cout << ii << ":" << sum<< "\n";
//                     }
//                     mosq_pop_file.close();

//     delete msq_mgr_copy;

// }


                std::vector<float> v_temp = vll_mgr.get_sum_prevalence_infected_of(
                        this->bsys_mgr.get_dominant_type_of(),
                        this->bsys_mgr.get_infectiousness()
                    );

                // std::cout << tt << "-v10" << v_temp[10] << "\n";

                // this->vll_mgr.update_prevalence_count_of(this->bsys_mgr.get_dominant_type_of());
                double before_step = util::get_wall_time();
                // this->msq_mgr->step(
                //     vll_mgr.get_sum_prevalence_infected_of(
                //         this->bsys_mgr.get_dominant_type_of(),
                //         this->bsys_mgr.get_infectiousness()
                //     )
                // );
                this->msq_mgr->step(
                    v_temp
                );
                stepper_time_cost.push_back(util::get_wall_time()-before_step);


                double num_infected = 0.0;
                double num_infectious = 0.0;

                double num_eggs = 0.0;
                double num_larv = 0.0;
                double num_imat = 0.0;

                double num_mosq = 0.0;
                double num_mosq_biting = 0.0;

                for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++) {

                    num_infected += this->msq_mgr->get_mosq_pop(vv, 5); ////Adult infected - biting
                    num_infectious += this->msq_mgr->get_mosq_pop(vv, 7); ////Adult infectious biting

                    num_eggs += this->msq_mgr->get_mosq_pop(vv, 0);
                    num_larv += this->msq_mgr->get_mosq_pop(vv, 1);
                    num_imat += this->msq_mgr->get_mosq_pop(vv, 2);

                    num_mosq += this->msq_mgr->get_mosq_pop(vv, 3);
                    num_mosq += this->msq_mgr->get_mosq_pop(vv, 4);
                    num_mosq += this->msq_mgr->get_mosq_pop(vv, 5);
                    num_mosq += this->msq_mgr->get_mosq_pop(vv, 6);
                    num_mosq += this->msq_mgr->get_mosq_pop(vv, 7);
                    num_mosq += this->msq_mgr->get_mosq_pop(vv, 8);

                    num_mosq_biting += this->msq_mgr->get_mosq_pop(vv, 3);
                    num_mosq_biting += this->msq_mgr->get_mosq_pop(vv, 5);
                    num_mosq_biting += this->msq_mgr->get_mosq_pop(vv, 7);


                    // std::cout << "this->msq_mgr->get_mosq_pop(vv, 5)" << this->msq_mgr->get_mosq_pop(vv, 5)
                    //             << ", this->msq_mgr->get_mosq_pop(vv, 7)" << this->msq_mgr->get_mosq_pop(vv, 7)
                    //             << "\n";
                    // // std::cin.ignore();

                }

                // std::cout << "ODE_num_eggs: " << num_eggs << "\n";
                // std::cout << "ODE_num_larv: " << num_larv << "\n";
                // std::cout << "ODE_num_imat: " << num_imat << "\n\n";
                // std::cout << "ODE_num_mosq: " << num_mosq << "\n";
                // std::cout << "ODE_num_mosq_b: " << num_mosq_biting << "\n";

                // if (tt == this->num_time_steps) {
                    // for (int ii = 0; ii < 9; ii++) {
                    //     float sum = 0;
                    //     for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++) {
                    //         sum += this->msq_mgr->get_mosq_pop(vv, ii);
                    //     }
                    //     std::cout << ii << ":" << sum<< "\n";
                    // }
                    // std::cin.ignore();
                // }

                // std::cout << "num_infected=" << num_infected
                //             << ", num_infectious=" << num_infectious << "\n";

                vll_mgr.daily_record_of_mosquito_infected_sum = static_cast<int>(num_infected);
                vll_mgr.daily_record_of_mosquito_infectious_sum = static_cast<int>(num_infectious);
                vll_mgr.daily_record_of_mosquito_clean_sum = static_cast<int>(num_mosq_biting-num_infected-num_infectious);

                vll_mgr.daily_record_of_eggs_sum = static_cast<int>(num_eggs);
                vll_mgr.daily_record_of_larv_sum = static_cast<int>(num_larv);
                vll_mgr.daily_record_of_imat_sum = static_cast<int>(num_imat);
                vll_mgr.daily_record_of_aliv_sum = static_cast<int>(num_mosq);
                vll_mgr.daily_record_of_aliv_biting_sum = static_cast<int>(num_mosq_biting);

                vll_mgr.daily_record_of_max_capacity_sum = static_cast<int>(this->msq_mgr->sum_num_larva_capacity);
                // vll_mgr.mosquito_ode_update_composition(
                //     0,
                //     num_infected,
                //     num_infectious
                // );

// if ((tt == 18250) || (tt == 21900)) {
//     std::ofstream temp_file;
//     temp_file.open("temp/capacity_"+std::to_string(tt)+".csv");
//     for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++) {
//         temp_file << vv << ",";
//         temp_file << this->msq_mgr->get_mosq_pop(vv,1) << ",";
//         temp_file << this->msq_mgr->get_mosq_capacity(vv) << "\n";
//     }
//     temp_file.close();
// }
                
// if ((tt >= 18250) && (tt <= (18250+370))) {
    
//     if (tt == 18250) {
//         temp_file.open("temp/capacity_one_year.csv");
//         temp_file << "cap,";
//         for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++) {
//             temp_file << this->msq_mgr->get_mosq_capacity(vv) << ",";
//         }
//         temp_file << "\n";
//     }

//     if ((tt % 10) == 0) {
//         temp_file << tt << ",";
//         for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++) {
//             temp_file << this->msq_mgr->get_mosq_pop(vv,1) << ",";
//         }
//         temp_file << "\n";
//     }

//     if (tt == (18250+370)) {
//         temp_file.close();
//     }
// }


// if (tt == (730)) {

//     std::vector<float> v_temp_mosq = vll_mgr.get_sum_prevalence_infected_of(
//             this->bsys_mgr.get_dominant_type_of(),
//             this->bsys_mgr.get_infectiousness()
//         );


//     std::ofstream mosq_pop_file_1;
//     mosq_pop_file_1.open ("temp_ode_mosq_pop_730.csv");
//         // float sum = 0;
//     for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++) {
//         // this->msq_mgr->set_mosq_time(vv, this->msq_mgr->get_mosq_time(vv));
//         // if ( vv > 500) {
//             mosq_pop_file_1 << vv << ": ";

//             // mosq_pop_file_1 << v_temp_mosq.at(vv);
//             for (int ii = 0; ii < 9; ii++) {
//                 mosq_pop_file_1 << " " << this->msq_mgr->get_mosq_pop(vv,ii);

//             }
//             mosq_pop_file_1 << "\n";
//         // }
//         // std::cout << ii << ":" << sum<< "\n";
//     }
//     mosq_pop_file_1.close();
//     std::cin.ignore();
// }

            } else {

                // this->vll_mgr.update_prevalence_count_of(this->bsys_mgr.get_dominant_type_of());
                
                double before_step = util::get_wall_time();
                
                this->vll_mgr.infection_h2m(
                    this->bsys_mgr.get_infectiousness(),
                    this->bsys_mgr.get_dominant_type_of(),
                    this->bsys_mgr.sum_num_systems,
                    tt
                );
                this->vll_mgr.mosquito_survival_step();
                this->vll_mgr.mosquito_incubation_step();

                stepper_time_cost.push_back(util::get_wall_time()-before_step);
                
            }


        ////// M2H Infection
        // std::cout << "[9]";



            if (!kMosq_method_using_ibm) {

                ///// temp start
                // std::vector<float> v_temp = vll_mgr.get_sum_prevalence_infected_of(
                //         this->bsys_mgr.get_dominant_type_of(),
                //         this->bsys_mgr.get_infectiousness()
                //     );
                ///// temp finish

                this->vll_mgr.infection_m2h(
                    this->hmn_mgr.get_bitten_by(),
                    this->hmn_mgr.sum_num_humans,
                    this->bsys_mgr.get_susceptibility(),
                    this->bsys_mgr.sum_num_systems,
                    this->bsys_mgr.get_most_advanced_stage(),
                    tt,
                    this->msq_mgr->get_num_infectious_mosq_of(),
                    // this->msq_mgr->get_feed_rate(),
                    kMosq_method_using_ode
                );
            
                // hmn_mgr.print_distribution<common::ParasiteType>(hmn_mgr.get_bitten_by(), hmn_mgr.sum_num_humans);
                this->bsys_mgr.parasite_intake(
                    this->hmn_mgr.get_bitten_by(),
                    common::StageName::kLiver,
                    // common::StageName::kBlood,
                    tt,
                    false
                );
            }

        ////// Infection Injection
        if (this->config["infection"]["injection"]["enabled"].GetBool()) {
            if (this->config["infection"]["injection"]["at_time_step"].GetInt() == tt ) {

                std::cout << "[*] Injection ...\n";

                this->bsys_mgr.step_report_genotype_weights();

                common::ParasiteType parasite_type_to_inject = common::ParasiteType::kRb;
                this->bsys_mgr.parasite_injection(
                    parasite_type_to_inject,
                    this->bsys_mgr.daily_record_of_parasite_positive_systems
                        * config["infection"]["injection"]["size_as_in_percentage_of_positive_population"].GetFloat()
                    // config["infection"]["injection"]["count"].GetInt()
                );

            }
        }
        
        if (this->config["simulation"]["reports"]["genotype_weights"].GetBool()) {
            this->bsys_mgr.step_report_genotype_weights();
        }

//EOD   ////// Mandatory end of day processes
        //std::cout << "[Step.end] ...\n";
            if ((!config["mosquito"]["ode"]["in_use"].GetBool())
                && (!config["mosquito"]["ibm"]["in_use"].GetBool())) {

                this->vll_mgr.update_summary();
            }

            if (config["simulation"]["reports"]["village_prev_annual"].GetBool()
                && (tt % 365 == 0)){
                    this->vll_mgr.update_prevalence_count_of(this->bsys_mgr.get_dominant_type_of());
                    this->vll_mgr.update_prevalence_of();
            }

            this->rpt_mgr.end_of_day_reports(tt);


            this->hmn_mgr.step_end_of_day(tt);


        // std::cin.ignore();

        if (tt == performance_timing_finish) {
            // std::cout << "performance_start_wall_time: " << performance_start_wall_time << "\n";

            performance_time_total = util::get_wall_time() - performance_start_wall_time;

            performance_cpu_time_duration = timer.now() - performance_start_cpu_time;


            // std::cout << "performance_time_total: " << performance_time_total << "\n";
            // std::cout << "util::get_wall_time(): " << util::get_wall_time() << "\n";
            // std::cin.ignore();
        }

        time_cost_per_time_step.push_back(util::get_wall_time()-wall_time_at_time_step_begin);
        // std::cout << "wall: " << util::get_wall_time()-wall_time_at_time_step_begin << std::endl;
        // boost::chrono::process_user_cpu_clock::duration step_cpu_duration= timer.now() - step_start_cpu_time;
        // double per_step_user_cpu_duration = boost::chrono::duration_cast<boost::chrono::milliseconds>(step_cpu_duration).count();
        // std::cout << step_cpu_duration << std::endl;
        
        // std::cout << per_step_user_cpu_duration << std::endl;
        // std::cin.ignore();

        // time_cost_per_time_step.push_back(per_step_user_cpu_duration);
        
        if (config["simulation"]["export_at_intervals"]["enabled"].GetBool()) {

            std::string out_dir = this->output_directory
                                    + util::kPathSeparator
                                    + this->output_prefix;

            if(!fs::exists(out_dir.c_str())){
                if(!fs::create_directory(out_dir.c_str())) {
                    std::cout << "Could not create directory at (" << out_dir << ")." << std::endl;
                    assert(false);
                }
            }



            if( (tt % config["simulation"]["export_at_intervals"]["interval_days"].GetInt()) == 0 ) {

                this->vll_mgr.output_human_current_locations(this->hmn_mgr);

                this->bsys_mgr.export_to_file(
                    out_dir
                        + util::kPathSeparator
                        + "on_day_" + std::to_string(tt),
                    &(this->hmn_mgr)
                );
            }
        }

    } // tt


    std::cout << std::fixed << std::setprecision(5)
                << "Measured performance timing: " << performance_time_total
                << " seconds, from " << performance_timing_start
                << " to " << performance_timing_finish << "\n";

    double d = boost::chrono::duration_cast<boost::chrono::milliseconds>(performance_cpu_time_duration).count();

    std::cout << std::fixed << std::setprecision(5)
                << "Measured performance User CPU timing: " << d
                << " ms, from " << performance_timing_start
                << " to " << performance_timing_finish << "\n";

    // std::cout << std::fixed << std::setprecision(5)
    //             << "Measured performance CPU timing: "
    //             << performance_timer_duration << ", "
    //             // << std::get<1>(performance_timer_duration)
    //             // << performance_timer_duration.get<0>()
    //             << boost::chrono::duration_cast<boost::chrono::milliseconds>(performance_timer_duration).count()
    //             << " ms, from " << performance_timing_start
    //             << " to " << performance_timing_finish << "\n";
    // std::cin.ignore();



    if (config["simulation"]["export"].GetBool()) {

        this->vll_mgr.output_human_current_locations(this->hmn_mgr);

        this->bsys_mgr.export_to_file(
            this->output_directory
                + util::kPathSeparator
                + this->output_prefix,
            &(this->hmn_mgr)
        );
    }

    std::cout << "record_of_cleared_pgs_rb:\n";
    for (auto ii : this->bsys_mgr.record_of_cleared_pgs_rb) {
        std::cout << ii << "\n";
    }
    std::cout << "record_of_cleared_pgs_r0:\n";
    for (auto ii : this->bsys_mgr.record_of_cleared_pgs_r0) {
        std::cout << ii << "\n";
    }

    std::cout << " pip treatment #:" << this->bsys_mgr.num_pip_treatments
                <<  "\npip treatment (success) #:" << this->bsys_mgr.num_pip_treatments_success
                << "\n";

    double wall_time_lapsed = util::get_wall_time() - wall_time_begin;

    std::cout << "death_total=" << death_total << std::endl;

    std::cout << std::fixed << std::setprecision(5)
                << "Time cost: " << wall_time_lapsed
                << " for " << this->num_time_steps << " steps, "
                << wall_time_lapsed / this->num_time_steps 
                << " per step.\n";


    if (this->config["intervention"]["fmda"]["enabled"].GetBool()) {
        std::ofstream output_file_survey(output_directory + util::kPathSeparator + output_prefix + "_fMDA_survey_records.csv");
        std::ostream_iterator<std::string> output_iterator(output_file_survey, "\n");
        std::copy(this->fmda.survey_record.begin(), this->fmda.survey_record.end(), output_iterator);
    }



    std::ofstream output_file(output_directory + util::kPathSeparator + output_prefix + "_stepper_time_cost.txt");
    std::ostream_iterator<double> output_iterator(output_file, "\n");
    std::copy(stepper_time_cost.begin(), stepper_time_cost.end(), output_iterator);


    std::ofstream output_file_time(output_directory + util::kPathSeparator + output_prefix + "_time_cost_per_step.txt");
    std::ostream_iterator<double> output_iterator_2(output_file_time, "\n");
    std::copy(time_cost_per_time_step.begin(), time_cost_per_time_step.end(), output_iterator_2);


    delete[] birth_death_reset_flags;

    delete[] random_numbers_humans;

    delete[] random_numbers_movement;
    delete[] random_numbers_drug_loss;
    delete[] random_number_repository_sys2;

}

void Simulator::close(){

    this->rpt_mgr.close_ofstreams();



    //////////////////////////////////////////////////////////////////////
    // Checks
    //////////////////////////////////////////////////////////////////////


    if (this->bsys_mgr.kDebug_assert) {

        int pg_io_counter_size = *std::max_element(
                                    this->bsys_mgr.debug_bites_accumulator_add.begin(),
                                    this->bsys_mgr.debug_bites_accumulator_add.end()
                                )+1;
        int pg_io_counter_add_max_index = std::max_element(
                                    this->bsys_mgr.debug_bites_accumulator_add.begin(),
                                    this->bsys_mgr.debug_bites_accumulator_add.end()
                                )-this->bsys_mgr.debug_bites_accumulator_add.begin();
        std::cout << "max add:" << pg_io_counter_size-1
                    << " (host " << pg_io_counter_add_max_index
                    << ": c_exp=" << int(this->bsys_mgr.get_cummulative_exposures_of()[pg_io_counter_add_max_index])
                    << ")\n";
        int* pg_io_counter = new int[pg_io_counter_size];
        std::fill_n(pg_io_counter, pg_io_counter_size, 0);
        int pg_io_counter_sum = 0;


        for (uint ii = 0; ii < this->bsys_mgr.debug_bites_accumulator_add.size(); ii++) {
            pg_io_counter[this->bsys_mgr.debug_bites_accumulator_add.at(ii)]++;
            pg_io_counter_sum += this->bsys_mgr.debug_bites_accumulator_add.at(ii);
        }
        std::cout << "sum add: " << pg_io_counter_sum << "\n";

        for (int cc = 0; cc < pg_io_counter_size; cc++) {
            std::cout << std::fixed << std::setprecision(1)
                        << cc << ": "
                        // << std::string(pg_io_counter[cc]*100/pg_io_counter_sum, '*')
                        << " " << pg_io_counter[cc] << "\n";
        }

        delete[] pg_io_counter;





        int pg_io_counter_less_size = -*std::min_element(
                                    this->bsys_mgr.debug_bites_accumulator_less.begin(),
                                    this->bsys_mgr.debug_bites_accumulator_less.end()
                                )+1;
        std::cout << "max less:-" << pg_io_counter_less_size-1 << "\n";
        int* pg_io_counter_less = new int[pg_io_counter_less_size];
        std::fill_n(pg_io_counter_less, pg_io_counter_less_size, 0);
        int pg_io_counter_less_sum = 0;


        for (uint ii = 0; ii < this->bsys_mgr.debug_bites_accumulator_less.size(); ii++) {
            pg_io_counter_less[-this->bsys_mgr.debug_bites_accumulator_less.at(ii)]++;
            pg_io_counter_less_sum += this->bsys_mgr.debug_bites_accumulator_less.at(ii);
        }
        std::cout << "sum less: " << pg_io_counter_less_sum << "\n";

        for (int cc = 0; cc < pg_io_counter_less_size; cc++) {
            std::cout << std::fixed << std::setprecision(1)
                        << "-" << cc << ": "
                        // << std::string(pg_io_counter_less[cc]*100/pg_io_counter_less_sum, '*')
                        << " " << pg_io_counter_less[cc] << "\n";
        }

        delete[] pg_io_counter_less;




        std::vector<int> bites_accumulator_both(this->bsys_mgr.debug_bites_accumulator_less.size(), 0);
        for (uint ii = 0; ii < bites_accumulator_both.size(); ii++) {
            bites_accumulator_both[ii] =
                this->bsys_mgr.debug_bites_accumulator_add[ii]
                + this->bsys_mgr.debug_bites_accumulator_less[ii];
        }

        int pg_io_counter_both_size = *std::max_element(
                                    bites_accumulator_both.begin(),
                                    bites_accumulator_both.end()
                                )+1;
        std::cout << "max:" << pg_io_counter_both_size-1 << "\n";
        int* pg_io_counter_both = new int[pg_io_counter_both_size];
        std::fill_n(pg_io_counter_both, pg_io_counter_both_size, 0);
        int pg_io_counter_both_sum = 0;


        for (uint ii = 0; ii < bites_accumulator_both.size(); ii++) {
            pg_io_counter_both[bites_accumulator_both.at(ii)]++;
            pg_io_counter_both_sum += bites_accumulator_both.at(ii);
        }
        std::cout << "sum: " << pg_io_counter_both_sum << "\n";

        for (int cc = 0; cc < pg_io_counter_both_size; cc++) {
            std::cout << std::fixed << std::setprecision(1)
                        << cc << ": "
                        // << std::string(pg_io_counter_both[cc]*100/pg_io_counter_both_sum, '*')
                        << " " << pg_io_counter_both[cc] << "\n";
        }

        delete[] pg_io_counter_both;

    }




    std::cout << "Number of ignored PG_q full incidents: "
            << this->bsys_mgr.sum_num_pg_q_is_full_incidents
            << "("
            << float(this->bsys_mgr.sum_num_pg_q_is_full_incidents)/this->bsys_mgr.sum_num_pg_intakes*100
            << "% "
            << "of "
            << this->bsys_mgr.sum_num_pg_intakes
            << " total pg intakes)\n";


    //Some Soft Checks
    std::cout << "\033[4mSoft Checks:\033[0m\n";

    // 1. The number of people in all three groups are the same
    std::cout << "1. Human.Population(" << this->hmn_mgr.sum_num_humans
                <<") == Village.Population(" << this->vll_mgr.sum_total_population
                <<") == BloodSystem.Count(" << this->bsys_mgr.sum_num_systems << ") ";

    if (this->hmn_mgr.sum_num_humans == this->vll_mgr.sum_total_population && this->vll_mgr.sum_total_population == this->bsys_mgr.sum_num_systems) {
        std::cout << "\033[32;47m[Pass]\033[0m";
    } else {
        std::cout << "\033[34;47m[Fail]\033[0m";
    }
    std::cout << "\n";

    // 2. No one lives/is located at unknown village
    int home_max = *std::max_element(this->hmn_mgr.home_village_of, this->hmn_mgr.home_village_of+this->hmn_mgr.sum_num_humans);
    int at_max = *std::max_element(this->hmn_mgr.at_village_of, this->hmn_mgr.at_village_of+this->hmn_mgr.sum_num_humans);
    std::cout   << "2. max{Human.home_village} (" << home_max 
                <<")  < Village.Count(" << this->vll_mgr.sum_num_villages
                <<") AND max{Human.at_village} (" << at_max
                <<") < Village.Count(" << this->vll_mgr.sum_num_villages << ") ";
    if (home_max < this->vll_mgr.sum_num_villages && at_max < this->vll_mgr.sum_num_villages) {
        std::cout << "\033[32;47m[Pass]\033[0m";
    } else {
        std::cout << "\033[34;47m[Fail]\033[0m";
    }
    std::cout << "\n";




}


void Simulator::output_mda(){
    this->mda.print_all_routes();
    if (!this->kIf_batch_mode){
        std::string sys_command = this->mda.print_map(
            this->output_directory + util::kPathSeparator + this->output_prefix,
            this->vll_mgr.longitude_of,
            this->vll_mgr.latitude_of,
            this->vll_mgr.population_of
        ) + " &";
        std::cout << sys_command << std::endl;
    
        if(system(sys_command.c_str())){}
    }

}

}