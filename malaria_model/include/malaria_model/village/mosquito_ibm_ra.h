#ifndef MOSQUITO_IBM_RA_H
#define MOSQUITO_IBM_RA_H

namespace village {

class MosquitoSwampRa;

class Dememosq {

    // const MosquitoSwampRa* kHome_swamp;
    MosquitoSwampRa* const kHome_swamp;
    int age = 0;
    int infectious_count = 0; // #times an infection was passed on to a human

    bool alive = false;
    bool pregnant = false;
    bool biting = false; //bitingstate: true - SEEKING, false - RESTING

    // enum class InfectionState : char {
    //     kSusceptible,
    //     kIncubating,
    //     kInfectious
    // };
    // InfectionState infstate;

    // int resistance;
    // int dominant;
public:
    struct Infection {
        common::ParasiteType pt;
        MosquitoSwampRa* origin;
    };
private:
    std::vector<Infection> incubation_list;
    std::vector<Infection> infectious_list;
    // std::vector<int> inforigin;
    // std::vector<int> incorigin;

public:
    Dememosq(
        MosquitoSwampRa* home_swamp,
        int init_age
    );
    ~Dememosq();

    inline void inc_age() {
        this->age++;
    }
    inline int get_age() {
        return this->age;
    }
    inline void set_age(int given_age) {
        this->age = given_age;
    }
    inline void inc_infectious_count() {
        this->infectious_count++;
    }
    // inline int get_infectious_count() {
    //     return this->infectious_count;
    // }
    inline void set_infectious_count(int given_count) {
        this->infectious_count = given_count;
    }
    inline void set_alive(bool if_alive) {
        this->alive = if_alive;
    }
    inline void set_pregnant(bool if_pregnant) {
        this->pregnant = if_pregnant;
    }
    inline void set_biting(bool if_biting) {
        this->biting = if_biting;
    }

    inline bool is_alive() const {
        return this->alive;
    }
    inline bool is_dead() const {
        return !this->alive;
    }
    inline bool is_biting() const {
        return this->biting;
    }
    inline bool is_pregnant() const {
        return this->pregnant;
    }
    inline MosquitoSwampRa* get_home_swamp() const {
        return this->kHome_swamp;
    }

    // inline void set_infstate_susceptible() {
    //     this->infstate = InfectionState::kSusceptible;
    // };
    // inline void set_infstate_incubating() {
    //     this->infstate = InfectionState::kIncubating;
    // };
    // inline void set_infstate_infectious() {
    //     this->infstate = InfectionState::kInfectious;
    // };


    inline void add_incubation(common::ParasiteType pt) {
        this->incubation_list.push_back(
            {pt, this->kHome_swamp}
        );
    }
    inline void add_incubation(common::ParasiteType pt, MosquitoSwampRa* orig_swamp) {
        this->incubation_list.push_back(
            {pt, orig_swamp}
        );
    }
    inline void add_infectious(common::ParasiteType pt) {
        this->infectious_list.push_back(
            {pt, this->kHome_swamp}
        );
    }
    inline void add_infectious(common::ParasiteType pt, MosquitoSwampRa* orig_swamp) {
        this->infectious_list.push_back(
            {pt, orig_swamp}
        );
    }
    inline bool has_incubation() const {
        return !(this->incubation_list.empty());
    }
    inline bool has_infectious() const {
        return !(this->infectious_list.empty());
    }
    inline int get_infectious_size() const {
        return static_cast<int>(this->infectious_list.size());
    }
    Infection select_rnd_infectious() const;

    void reset();
    void turn_dead();
    void turn_alive();

    void step_layeggs(float G_eggslayed);

    int step_incubate(float incubation_prob);
    void step_pregnancy(
        float G_mpregnancy,
        float G_mgestation,
        float G_eggslayed
    );
    void step_biting_resting_switch(
        float kBiting_resting_prob
    );

};

class MosquitoSwampRa {

    const float G_egg1=6.16E+07;      // Egg_Arrhenius1
    const float G_egg2=5754.03;     // Egg_Arrhenius2
    const float G_larva1=8.42E+10;   // aquatic_arrhenius_1
    const float G_larva2=8328.0;    // aquatic_arrhenius_2
    // const float G_pupa1=8.42e10;
    // const float G_pupa2=8328.0;

    const float G_aquatic_mortality;
    // const float G_aquatic_mortality=0.1;

    const float G_mature;
    // const float G_mature=1.0/2.0;
    const float G_iadult_life; // immature death rate was 1/5
    // const float G_iadult_life=1.0/7.0; // immature death rate was 1/5
    // const float G_iadult_life=1.0/5.0; // immature death rate was 1/5

    const float G_mosqdeath1;   //life expectancy of an infectious mosquito
    // const float G_mosqdeath1 = 1.0/24.0;   //life expectancy of an infectious mosquito
    // const float G_mosqdeath1 = 1.0/10.0;   //life expectancy of an infectious mosquito

    // startvillagestrigger
    // for (int i=0; i<G_xregion; i++){
    //     village=G_villages[i];
    //     village->mda = NO;
    //     village->mda_days = 0;
    //     village->mda_times = 0;
    //     village->counter = 0;
    //     village->post = 0;
    //     village->delay = 0;
    //     village->mda_start = 0;
    //     village->eggs=150;
    //     village->larva=150;
    //     village->pupa=150;
    //     village->iadult=80;
    // }
    
    // initialised in constructor
    const int kHost_village_id;
    const int nummosquito_i;
    const int kG_Dememosquitoes_start_index;
    const float Kmax_i; // Max_Larval_Capacity

    // int eggs=150;
    // int larva=150;
    // int pupa=150;
    // int iadult=80;

    int eggs=0;
    int larva=0;
    int pupa=0;
    int iadult=0;

    int malive=0;
    int mdead=0;

    int step_expected_new_adult = 0;
    float step_expected_new_adult_prob = 0.0;
    int step_actual_new_adult_counter = 0;

    int step_expected_death = 0;
    float step_expected_death_prob = 0.0;
    int step_actual_death_counter = 0;

    int step_actual_death_infectious_counter = 0;

    int step_death_infected = 0;
    int step_death_infectious = 0;


public:

    MosquitoSwampRa(
        float aquatic_mortality,
        float maturation_days,
        float immat_life_days,
        float adult_life_days,
        int host_village_id,
        int vv_nummosquito_i,
        int vv_G_Dememosquitoes_start_index,
        float vv_Kmax_i
    );
    ~MosquitoSwampRa();

private:
    void step_reset();

public:
    int sum_num_mosq_slots_not_enough_incidents = 0;

    void step(
        float K_i
    );

    bool step_inc_born_death(
        float K_i,
        std::vector<Dememosq*> G_Dememosquitoes
    );

    bool step_inc_born_death_rev(
        float K_i,
        std::vector<Dememosq*> G_Dememosquitoes
    );

    bool step_inc_born_death_ind(
        float K_i,
        std::vector<Dememosq*> G_Dememosquitoes
    );

    inline float get_step_expected_new_adult_prob() const {
        return this->step_expected_new_adult_prob;
    };
    inline bool step_if_expected_new_adult_num_reached() const {
        return (this->step_actual_new_adult_counter >= this->step_expected_new_adult);
    };
    inline void step_register_new_adult() {
        this->step_actual_new_adult_counter++;
    };
    inline float get_step_expected_death_prob() const {
        return this->step_expected_death_prob;
    };
    inline bool step_if_expected_death_num_reached() const {
        return (this->step_actual_death_counter >= this->step_expected_death);
    };
    inline void step_register_death(int infectious_count) {
        this->step_actual_death_counter++;
        this->step_actual_death_infectious_counter+=infectious_count;
    };
    inline int get_step_expected_new_adult() const {
        return this->step_expected_new_adult;
    }
    inline int get_step_expected_death() const {
        return this->step_expected_death;
    }
    

    inline int get_host_village_id() const {
        return this->kHost_village_id;
    };
    inline bool has_iadault() const {
        return (this->iadult>0);
    };


    inline void inc_alive() {
        this->malive++;
    };
    inline void inc_dead() {
        this->mdead++;
    };
    inline void inc_egg(int new_num_eggs) {
        this->eggs += new_num_eggs;
    };

    inline int get_status_eggs() const {
        return this->eggs;
    };
    inline int get_status_larva() const {
        return this->larva;
    };
    inline int get_status_pupa() const {
        return this->pupa;
    };
    inline int get_status_iadult() const {
        return this->iadult;
    };
    inline int get_status_malive() const {
        return this->malive;
    };
    inline int get_status_mdead() const {
        return this->mdead;
    };

    inline int get_step_actual_death_counter() const {
        return this->step_actual_death_counter;
    }
    inline int get_step_actual_death_infectious_counter() const {
        return this->step_actual_death_infectious_counter;
    }

    inline int get_step_death_infected() const {
        return this->step_death_infected;
    }
    inline int get_step_death_infectious() const {
        return this->step_death_infectious;
    }

    // inline int get_nummosquito() const{
    //     return this->nummosquito_i;
    // }

};

class MosquitoManagerIbmRa {

    // const float kResting_to_biting_rate = 0.9; // 0.8 lay_rate
    const float kResting_to_biting_rate; // 1/resting_days, step_matching_ode
    // const float G_bitingrate = 0.175;
    // const float G_bitingrate = 0.14; // [mosq][biting_rate]
    const float kBiting_rate; // [mosq][biting_rate]
    // const float G_deathbite = 1.0/10.0;
    // const float G_deathbite = 0.05;
    // const float G_deathbite_transmission = 1.0;
    const float G_deathbite = 0.00;
    const float G_deathbite_transmission = 1.0;


    const float G_mosqinc;     // extrinsic incubation period, 1.0/incubation_days
    // const float G_mosqinc = 1.0/14.0;     // extrinsic incubation period

    // const float G_mpregnancy = 1.0/2.0;
    // const float G_mgestation = 1.0/4.0;
    const float G_mpregnancy = 1.0/1.0;
    const float G_mgestation = 1.0/1.0;



    const float G_eggslayed; // eggs_per_lay
    // const float G_eggslayed = 20.0;
    // const float G_eggslayed = 10.0;
    // const float G_eggslayed = 4.0;

    // const float kBiting_resting_prob = 0.33;
    // const float kBiting_resting_prob = 0.5;
    const float kBiting_resting_prob = 0.5; // bite / rest switch, step


    // Initialisation

    // const float G_carryingcapacity = 1.5;
    const float kG_carryingcapacity;
    const float betaf = 0.095;
    // const float scaletomosqnumbers = 1300.0;
    const float kScaletomosqnumbers;

    // const float kInitProbAlive = 0.5;
    // const float kInitProbPregnant = 1.0; // RA 0.1
    // const float kInitProbBiting = 1.0; // RA 0.2
    // const float kInitProbInfectious = 0.05;

    const float kInitProbAlive = 0.0;
    const float kInitProbPregnant = 1.0; // RA 0.1
    const float kInitProbBiting = 1.0; // RA 0.2
    const float kInitProbInfectious = 0.05;

    const int kFunc_swamp_step_inc_born_death_version;

private:
    std::vector<MosquitoSwampRa*> mosquito_swamp_of;
    std::vector<Dememosq*> G_Dememosquitoes;

public:

    int sum_num_infected_biting = 0;
    int sum_num_infectious_biting = 0;
    int sum_num_clean_biting = 0;

    int sum_num_new_infected = 0;
    int sum_num_incubation = 0;
    int sum_num_death_infected = 0;
    int sum_num_death_infectious = 0;

    int sum_num_eggs = 0;
    int sum_num_larv = 0;
    int sum_num_imat = 0;
    int sum_num_malive = 0;
    int sum_num_malive_biting = 0;

    int sum_num_new_eggs = 0;

    int sum_num_new_adults = 0;
    int sum_num_death_adults = 0;

    int sum_num_death_adults_infectious_count = 0;

    int sum_num_infectious_bites = 0;
    // int sum_num_successful_infectoius_bites = 0;
    int sum_num_successful_infectious_bites = 0;

    int sum_num_larva_capacity = 0;

    // int sum_num_infectious_mosq = 0;

    int sum_num_superinfection_requests = 0;

    //////////////////////////////////////////////////////////////////////
    // Constructors ODE

    MosquitoManagerIbmRa(
        float aquatic_mortality,
        float maturation_days,
        float immat_life_days,
        float adult_life_days,

        float resting_days,
        float biting_rate,
        float incubation_days,
        float eggs_per_lay,
        
        float g_carryingcapacity,
        float scaletomosqnumbers,
        int func_swamp_step_inc_born_death_version,
        int num_villages,
        const int* population_of,
        const float* transmission_coefficient_of
    );


    ~MosquitoManagerIbmRa();

    //////////////////////////////////////////////////////////////////////
    // 

    void step(
        int timestep,
        const VillageManager& vll_mgr,
        // m2h
        const float* susceptibility_array,
        const float susceptibility_multiplier, // G_c
        const common::StageName* most_advanced_stage_array,
        common::ParasiteType* infections_array,
        const int infections_array_size,
        // h2m
        const float* infectiousness_array,
        const float infectiousness_multiplier, // G_c
        const common::ParasiteType* dominant_type_array,
        const int host_array_size,
        const int* human_home_village_array
    );

    void step_matching_ode(
        int timestep,
        const VillageManager& vll_mgr,
        // m2h
        const float* susceptibility_array,
        const float susceptibility_multiplier, // G_c
        const common::StageName* most_advanced_stage_array,
        common::ParasiteType* infections_array,
        const int infections_array_size,
        // h2m
        const float* infectiousness_array,
        const float infectiousness_multiplier, // G_c
        const common::ParasiteType* dominant_type_array,
        const int host_array_size,
        const int* human_home_village_array
    );

    std::vector<int> get_mosq_pop_of_swamp(int swamp_id);

    void step_print_stats() const;

    //////////////////////////////////////////////////////////////////////
    // Access

    inline MosquitoSwampRa* get_swamp_of_village(int village_id) const {
        assert(village_id < int(this->mosquito_swamp_of.size()));
        return this->mosquito_swamp_of[village_id];
    };

//     inline double get_num_infectious_mosq(int mosq_index) const {
//         return this->num_infectious_mosq_of[mosq_index];
//     }

//     inline const std::vector<double>& get_num_infectious_mosq_of() const {
//         return this->num_infectious_mosq_of;
//     }

//     inline double get_mosq_pop(int mosq_index, int pop_index) const {
//         return this->mosquito_of[mosq_index]->get_pop(pop_index);
//     }

//     inline void set_verbose_on() {
//         this->verbose = true;
//     }

//     inline void set_verbose_off() {
//         this->verbose = false;
//     }



};

}


#endif