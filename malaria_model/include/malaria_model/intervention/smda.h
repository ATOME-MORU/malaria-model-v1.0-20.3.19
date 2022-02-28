#ifndef SMDA_H
#define SMDA_H


namespace intervention {

// Selected MDA:
// Interface class for MDAs that are planned only
// on selected villages.
class Smda{ 
protected:
    // Mda parameters
    const int kNum_teams;

    const int kNum_days_spent_in_village;
    const int kNum_drug_rounds;
    const int kNum_days_between_drug_rounds;
    
    // Smda parameters

    std::vector<int> scheduled_start_times;
    std::vector<Mda> triggered_mdas;

public:
    Smda(
        // Smda
        const int num_teams,

        const int num_days_spent_in_village,
        const int num_drug_rounds,
        const int num_days_between_drug_rounds

    ) : 
        kNum_teams(num_teams),

        kNum_days_spent_in_village(num_days_spent_in_village),
        kNum_drug_rounds(num_drug_rounds),
        kNum_days_between_drug_rounds(num_days_between_drug_rounds) {

    }

    virtual ~Smda(){}
    inline void register_start_time(int time_step){
        this->scheduled_start_times.push_back(time_step);
    }

    inline int get_next_start_time() const {
        if (this->scheduled_start_times.empty()) {
            return -1;
        } else {
            return this->scheduled_start_times.front();
        }
    }

    virtual void step(
        const int current_time_step,
        village::VillageManager* vll_mgr_ptr
    ){
        (void) current_time_step;
        (void) vll_mgr_ptr;
    }

    virtual void create_mda(
        const int current_time_step,
        const std::vector<float>& prevalence_per_village
    ){
        (void) current_time_step;
        (void) prevalence_per_village;
    }

};

}

#endif