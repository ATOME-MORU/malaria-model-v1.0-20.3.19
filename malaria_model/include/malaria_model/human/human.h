#ifndef HUMAN_H
#define HUMAN_H

// #include <vector>
#include "common/common.h"

namespace village{class Village;}
// namespace common{struct Parasite;}

namespace human {

constexpr int kNum_age_bins = 100;


class Human {

    // public:

    //     int*  const age;
    //     char* const gender; // 'm' or 'f'
    //     int*  const home_village;
    //     int*  const at_village;
    //     common::DrugName* const given_drug;

    // Human(
    //     int*  const a,
    //     char* const g,
    //     int*  const ho,
    //     int*  const at,
    //     common::DrugName* const dg
    //     ) : age(a), gender(g), home_village(ho), at_village(at), given_drug(dg){}

    // note: a const pointer / by reference would be safer.
    //       but because of the const property, we'd have to use
    //       a vector to contain human_reg.
    //       a vector stl is not native on GPUs.

    public:

        int*  age;
        char* gender; // 'm' or 'f'
        int*  home_village;
        int*  at_village;
        common::DrugName *given_drug;

    // Human(
    //     int* a,
    //     char* g,
    //     int* ho,
    //     int* at,
    //     common::DrugName *dg
    //     ) : age(a), gender(g), home_village(ho), at_village(at), given_drug(dg){}

    // Human(){}
        // const int* age;
        // const char* gender; // 'm' or 'f'
        // const int* home_village;
        // const int* at_village;
        // common::DrugName* const given_drug;


        // const village::Village* village_home;
        // const village::Village* village_at;
};


class HumanManager {

    const int kFunc_death_probability_version;

    const bool kAttractiveness_enabled;
    const float kAttractiveness_mean;
    const float kAttractiveness_sd;
    std::vector<float> attractiveness_of;
    const int kReset_log_bitten_times_on_day;
    std::vector<int> log_bitten_times_m2h; // size H, accumulated value over time
    std::vector<int> log_bitten_times_h2m; // size H, accumulated value over time

    float func_gen_attractiveness() const;
    void init_attractiveness();


    public:

    // column/attribute-wise data silos
    int* age_of;
    char* gender_of;
    float* death_probability_of;
    // village::Village* village_home_of;
    // village::Village* village_at_of;
    int* home_village_of;
    int* at_village_of;

    bool* if_itn_of;

    common::HumanMobilityType* mobility_of; // TODO: move type from common to human

    // row/object-wise data points/registers
    Human* human_reg;
    // std::vector<Human> human_reg;

    // std::list<common::Parasite> parasite_reg;
    // int* parasite_offset;

    std::vector<int> age_weight_vector{std::vector<int>(kNum_age_bins,0)};
    std::vector<float> death_prob_vector{std::vector<float>(kNum_age_bins,0.0)};



    common::ParasiteType* bitten_by_of;
    common::DrugName* given_drug_of;

    // summary data
    int sum_num_humans = 0;

    HumanManager(
        int initial_population
    );
    HumanManager(
        int initial_population,
        int func_death_probability_version,

        bool attractiveness_enabled,
        float attractiveness_mean,
        float attractiveness_sd,
        float reset_log_bitten_times_at_the_end_of_year
    );
    ~HumanManager();

    void print_all() const;
    void print_one(const int index_human) const;

    void read_age_distribution(std::string input_file_name);
    void read_death_prob_distribution(std::string input_file_name);

    void init_age();
    void init_age_discrete_distribution();
    void init_age_mixed_normal();
    void print_age_distribution() const;

    int birth_and_death(
        bool* reset_flags,
        const float* random_numbers_array,
        const int random_numbers_array_size
    );
    // bool birth_and_death_probability_based(
    //     const int h_index,
    //     const int* ages,
    //     const float* random_numbers
    // );
    void update_all_death_probability();
    float func_death_probability (int age_in_years);
    static float func_death_probability_kernel_0 (int age_in_years);
    static float func_death_probability_kernel_1_0 (int age_in_years);
    static float func_death_probability_kernel_1_1 (int age_in_years);
    static float func_death_probability_kernel_2 (int age_in_years);
    static float func_death_probability_kernel_2_1 (int age_in_years);
    float func_death_probability_kernel_3 (int age_in_years);
    static float death_probability_age_in_years_kernel(
        const int age_in_years,
        const float G_sig,
        const float G_mu1,
        const float G_mu
    );
    void reset_human_to_birth(const int h_index);

    void increase_all_ages_by_one();

    void init_gender(const float male_probability);
    void print_gender_distribution() const;


    void init_mobility(const float* village_init_mobility_percentages);

    void init_itn(float itn_probability);
    static void init_itn_kernel( //TODO: worth templating this?
        float probability,
        bool* data,
        float* random_numbers,
        uint data_size
    );



    inline void set_human_home_village(const int human_id, const int village_id) {
        this->home_village_of[human_id] = village_id;
    }

    inline void set_human_at_village(const int human_id, const int village_id) {
        this->at_village_of[human_id] = village_id;
    }
    inline void set_human_return_home_village(const int human_id) {
        this->at_village_of[human_id] = this->home_village_of[human_id];
    }

    inline void set_human_age(int human_id, int age){
        this->age_of[human_id] = age;
    }

    inline void set_human_attractiveness(int human_id, float attractiveness){
        this->attractiveness_of[human_id] = attractiveness;
    }

    inline void inc_log_bitten_times_m2h(int human_id) {
        this->log_bitten_times_m2h[human_id]++;
    }
    inline void inc_log_bitten_times_h2m(int human_id) {
        this->log_bitten_times_h2m[human_id]++;
    }

    // void init_parasites();

    // void update_bitten_by(float* village_bite_rate, float* village_ra_rate);
    common::ParasiteType* get_bitten_by();
    void reset_bitten_by();

    // void update_given_drug_by_at_village(const int index_village,
    //                                      const common::DrugName drug);


    inline common::DrugName* get_given_drug() const {
        return this->given_drug_of;
    }
    void clear_given_drug();


    ////////////////////////
    //

    inline void reset_log_bitten_times() {
        std::fill(log_bitten_times_m2h.begin(), log_bitten_times_m2h.end(), 0);
        std::fill(log_bitten_times_h2m.begin(), log_bitten_times_h2m.end(), 0);
    }

    void step_end_of_day(int time_step);


    ////////////////////////
    // Access

    inline int get_log_bitten_times_m2h(int h_index) const {
        return this->log_bitten_times_m2h[h_index];
    }
    inline int get_log_bitten_times_h2m(int h_index) const {
        return this->log_bitten_times_h2m[h_index];
    }

    // At village
    inline int get_at_village(int h_index) const{
        return this->at_village_of[h_index];
    }
    inline const int* get_at_village_of() const{
        return this->at_village_of;
    }

    // Home village
    inline int get_home_village(int h_index) const {
        return this->home_village_of[h_index];
    }
    inline const int* get_home_village_of() const{
        return this->home_village_of;
    }

    // Age
    inline const int* get_age_of() const{
        return this->age_of;
    }
    inline const bool* get_if_itn_of() const{
        return this->if_itn_of;
    }

    //
    inline const common::HumanMobilityType* get_mobility_of() const{
        return this->mobility_of;
    }

    //
    inline bool get_if_attractiveness_enabled() const {
        return this->kAttractiveness_enabled;
    }
    inline float get_attractiveness(int h_index) const {
        return this->attractiveness_of[h_index];
    }

    template <class T>
    static void print_distribution(T* data_array, int length);


};

}

#endif