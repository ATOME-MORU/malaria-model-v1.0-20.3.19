#include <vector>
#include <iostream>
#include <math.h> // pow, sqrt, cos
#include <cassert>

#include "common/common.h"

#include "village/village.h"
#include "village/mosquito_ibm_ra.h"
#include "util/randomness.h"


namespace village {


Dememosq::Dememosq(
        MosquitoSwampRa* const home_swamp,
        int init_age
    ) : 
        kHome_swamp(home_swamp)
    {
        this->reset();
        this->set_age(init_age);
}

Dememosq::~Dememosq() {

}

Dememosq::Infection Dememosq::select_rnd_infectious() const {
    //         // std::cout << "OK here ...BB" << std::endl; 
    // int index = int(util::get_rand_uniform() * int(this->infectious_list.size()));
    //         // std::cout << "index=" << index << std::endl; 
    //         // std::cout << "this->infectious_list.size()=" << this->infectious_list.size() << std::endl; 

    // return this->infectious_list[index];

    return this->infectious_list[int(util::get_rand_uniform() * int(this->infectious_list.size()))];
}

void Dememosq::reset() {
    this->set_age(0);
    this->set_infectious_count(0);
    this->set_alive(false);
    this->set_pregnant(false);
    this->set_biting(false);
    this->incubation_list.clear();
    this->infectious_list.clear();
}

void Dememosq::turn_dead() {
    this->get_home_swamp()->step_register_death(this->infectious_count);
    // if (this->has_infectious()) {
    //     std::cout << "inft_q_size/age/inft_c;"
    //                 << this->get_infectious_size() << ";"
    //                 << this->age << ";"
    //                 << this->infectious_count << "\n";
    // }
    this->reset();
}

void Dememosq::turn_alive() {
    this->reset();
    this->set_alive(true);
    this->set_biting(true);
    this->get_home_swamp()->step_register_new_adult();
}

int Dememosq::step_incubate(float prob_incubation) {
    int c = 0;
    for (std::vector<Infection>::iterator inf_it = this->incubation_list.begin(); inf_it != this->incubation_list.end(); ++inf_it) {
        if (util::get_rand_uniform() < prob_incubation){
            this->infectious_list.push_back(*(inf_it));
            this->incubation_list.erase(inf_it--);
            c++;
            // break;
        }
    }
    return c;

    //// When no member is const in struct
    // for (std::vector<Infection>::iterator inf_it = this->incubation_list.begin(); inf_it != this->incubation_list.end(); ++inf_it) {
    //     if (util::get_rand_uniform() < prob_incubation){
    //         this->infectious_list.insert(
    //             this->infectious_list.end(),
    //             std::make_move_iterator(inf_it),
    //             std::make_move_iterator(inf_it+1)
    //         );
    //         this->incubation_list.erase(inf_it--);
    //     }
    // }
}

void Dememosq::step_layeggs(
        float G_eggslayed
    ) {
    // this->get_home_swamp()->inc_egg(util::get_rand_uniform() * G_eggslayed);
    this->get_home_swamp()->inc_egg(util::get_rand_normal(G_eggslayed, G_eggslayed*0.2));
}

void Dememosq::step_pregnancy(
        float G_mpregnancy,
        float G_mgestation,
        float G_eggslayed
    ) {
    if (!this->is_pregnant()) {
        if ( util::get_rand_uniform() < G_mpregnancy) {
            this->set_pregnant(true);
        }
    } else {
        // if (!this->is_biting()) {

        //     if (util::get_rand_uniform() < G_mgestation) {
        //         int new_egges = int(util::get_rand_normal(G_eggslayed, G_eggslayed*0.1));
        //         this->get_home_swamp()->inc_egg(new_egges);
        //         this->set_pregnant(false);
        //     }
        // }
        if (util::get_rand_uniform() < G_mgestation) {
            int new_eggs = G_eggslayed;
            this->get_home_swamp()->inc_egg(new_eggs);
            // this->set_pregnant(false);
        }
    }
}

void Dememosq::step_biting_resting_switch(
        float kBiting_resting_prob
    ) {
    if (util::get_rand_uniform() < kBiting_resting_prob ) { 
        this->set_biting(!this->is_biting());
    }

}



MosquitoSwampRa::MosquitoSwampRa(
        float aquatic_mortality,
        float maturation_days,
        float immat_life_days,
        float adult_life_days,

        int host_village_id,
        int vv_nummosquito_i,
        int vv_G_Dememosquitoes_start_index,
        float vv_Kmax_i
    ) : 
        G_aquatic_mortality(aquatic_mortality),
        G_mature(1.0/maturation_days),
        G_iadult_life(1.0/immat_life_days),
        G_mosqdeath1(1.0/adult_life_days),
        kHost_village_id(host_village_id),
        nummosquito_i(vv_nummosquito_i),
        kG_Dememosquitoes_start_index(vv_G_Dememosquitoes_start_index),
        Kmax_i(vv_Kmax_i)
    {

        this->eggs = Kmax_i * 0.5;
        this->larva = Kmax_i * 0.5;
        this->iadult = Kmax_i * 0.5;

        this->malive = 0;
        this->mdead = nummosquito_i - this->malive;
}

MosquitoSwampRa::~MosquitoSwampRa() {

}

void MosquitoSwampRa::step_reset() {
    this->step_expected_new_adult = 0;
    this->step_expected_new_adult_prob = 0.0;
    this->step_actual_new_adult_counter = 0;

    this->step_expected_death = 0;
    this->step_expected_death_prob = 0.0;
    this->step_actual_death_counter = 0;

    this->step_actual_death_infectious_counter = 0;

    this->step_death_infected = 0;
    this->step_death_infectious = 0;
}

void MosquitoSwampRa::step(
        float K_i
    ) {

    this->step_reset();

    // line 2148, villages(int timestep, float seasons)
    float egg2larva_i=G_egg1*exp(-(G_egg2/K_i));
    float larva2pupa_i=G_larva1*exp(-(G_larva2/K_i));
    // float pupa2adult_i=G_pupa1*exp(-(G_pupa2/K_i));

    // std::cout << "larva2pupa_i=" << larva2pupa_i << "\n";
    

    // std::cout << K_i << ", " << egg2larva_i << ", " << larva2pupa_i << "\n";

    // if (kHost_village_id < 3) {
    //     std::cout << "v" << kHost_village_id << ", Kmax_i=" << Kmax_i;
    //     std::cout << ", Pre-step:\n";
    //     std::cout << "eggs: " << this->eggs
    //                 << ", larva: " << this->larva
    //                 << ", pupa: " << this->pupa
    //                 << ", iadult: " << this->iadult
    //                 << ", malive: " << this->malive
    //                 << ", mdead: " << this->mdead << "\n";
    // }


    // line 2156
    if (this->eggs > 0) {
            /** binomial with 1-(not turning to larvae and not dying) is the total events **/
            float pnotegg=1-(1-egg2larva_i)*(1-G_aquatic_mortality);
            // boost::binomial_distribution<> my_binomial(village->eggs,pnotegg);  // binomial distribution
            // boost::variate_generator<boost::mt19937&, boost::binomial_distribution<> >next_value(rng, my_binomial);     // glues randomness with mapping
            // int neggs=next_value();
            // std::cout<<"pnotegg="<<pnotegg<<std::endl;
            // std::cout<<"this->eggs"<<this->eggs<<std::endl;

            int neggs=util::get_rand_binomial(this->eggs,pnotegg);
            // std::cout<<"neggs="<<neggs<<std::endl;

            //%creating normalized probability for each of the events
            float normprobegg=egg2larva_i/(egg2larva_i+G_aquatic_mortality);
            // boost::binomial_distribution<> my_binomial2(neggs,normprobegg);  // binomial distribution
            // boost::variate_generator<boost::mt19937&, boost::binomial_distribution<> >next_value2(rng, my_binomial2);     // glues randomness with mapping
            // int movedeggs=next_value2();
            int movedeggs=util::get_rand_binomial(neggs,normprobegg);
            // std::cout<<"movedeggs="<<movedeggs<<std::endl;

            //%keep all who didn't move
            // village->eggs=village->eggs-neggs;
            this->eggs=this->eggs-neggs;

            //%move to larva after using a logistic term!
            // boost::binomial_distribution<> my_binomial3(movedeggs,(1-village->larva/Kmax_i));  // binomial distribution
            // boost::variate_generator<boost::mt19937&, boost::binomial_distribution<> >next_value3(rng, my_binomial3);     // glues randomness with mapping
            // int surveggstolarva=next_value3();
            int surveggstolarva=util::get_rand_binomial(movedeggs,(1-this->larva/Kmax_i));
            // std::cout<<"surveggstolarva="<<surveggstolarva<<std::endl;

            // village->larva=village->larva+surveggstolarva;
            // if (village->larva>Kmax_i){village->larva=Kmax_i;}

            this->larva=this->larva+surveggstolarva;
            if (this->larva > Kmax_i){this->larva = Kmax_i;}
    }
    // std::cin.ignore();
    /** Larva to iadult **/
    if (this->larva > 0) {
        float pnotlarva=1-(1-G_aquatic_mortality)*(1-larva2pupa_i);
        // std::cout<<"pnotlarva="<<pnotlarva<<std::endl;

        // boost::binomial_distribution<> my_binomial(village->larva,pnotlarva);  // binomial distribution
        // boost::variate_generator<boost::mt19937&, boost::binomial_distribution<> >next_value(rng, my_binomial);     // glues randomness with mapping
        // int nlarva=next_value();
        int nlarva=util::get_rand_binomial(this->larva,pnotlarva);
        // std::cout<<"nlarva="<<nlarva<<std::endl;

        //%keep all who didn't move
        this->larva=this->larva-nlarva;

        float normproblarv=larva2pupa_i/(larva2pupa_i+G_aquatic_mortality);
        // std::cout<<"normproblarv="<<normproblarv<<std::endl;
        // boost::binomial_distribution<> my_binomial2(nlarva,normproblarv);  // binomial distribution
        // boost::variate_generator<boost::mt19937&, boost::binomial_distribution<> >next_value2(rng, my_binomial2);     // glues randomness with mapping
        // int movedlarva=next_value2();
        int movedlarva=util::get_rand_binomial(nlarva,normproblarv);
        // std::cout<<"movedlarva="<<movedlarva<<std::endl;

        //%move to larva after using a logistic term!
        // boost::binomial_distribution<> my_binomial3(movedlarva,(1-village->pupa/Kmax_i));  // binomial distribution
        // boost::variate_generator<boost::mt19937&, boost::binomial_distribution<> >next_value3(rng, my_binomial3);     // glues randomness with mapping
        // int survlarvatopupa=next_value3();
        int survlarvatopupa=util::get_rand_binomial(movedlarva,(1-this->iadult/Kmax_i));


        // village->iadult=village->iadult+survlarvatopupa;
        // if (village->iadult>Kmax_i){village->iadult=Kmax_i;}
        this->iadult=this->iadult+survlarvatopupa;

        if (this->iadult > Kmax_i){this->iadult = Kmax_i;}
           // std::cout<<"this->iadult:"<<this->iadult<<std::endl;
           // std::cout<<", Kmax[i]:"<<Kmax_i<<std::endl;
           // std::cin.ignore();

        // this->iadult+=movedlarva;
        // if (this->iadult > Kmax_i) {
        //     this->iadult = Kmax_i;
        //     std::cout << "iadult:"<< this->iadult << " == Kmax_i:" << Kmax_i << "\n";
        //     std::cin.ignore();
        // }
    }

    int prev_step_malive = this->malive;



    // line 1312, mosquitodynamicsnew(int timestep)
    /** immature to mature/death **/
    if(this->iadult>0){
        float pnotimmature=1-(1-G_mature)*(1-G_iadult_life); 
        // boost::binomial_distribution<> my_binomial(village->iadult,pnotimmature);  // binomial distribution
        // boost::variate_generator<boost::mt19937&, boost::binomial_distribution<> >next_value(rng, my_binomial);     // glues randomness with mapping
        // int nimmature=next_value();
        int nimmature=util::get_rand_binomial(this->iadult,pnotimmature);
        // village->iadult=village->iadult-nimmature;
        this->iadult=this->iadult-nimmature;

        float normprobiadult=G_mature/(G_mature+G_iadult_life);
        // boost::binomial_distribution<> my_binomial2(nimmature,normprobiadult);  // binomial distribution
        // boost::variate_generator<boost::mt19937&, boost::binomial_distribution<> >next_value2(rng, my_binomial2);     // glues randomness with mapping
        // int movediadult=next_value2();
        int movediadult=util::get_rand_binomial(nimmature, normprobiadult);
        // y[j]=movediadult;
        step_expected_new_adult = movediadult;
        // cout<<"moved inc: "<<movediadult<<endl;
        // village->malive=village->malive+movediadult;
        this->malive = this->malive + movediadult;
        // village->mdead=village->mdead-movediadult;
    }
        // cout<<"Second: "<<village->malive<<endl;

    if(this->malive>0){
        float pnotalive=1-(1-G_mosqdeath1);
        // boost::binomial_distribution<> my_binomial(village->malive,pnotalive);  // binomial distribution
        // boost::variate_generator<boost::mt19937&, boost::binomial_distribution<> >next_value(rng, my_binomial);     // glues randomness with mapping
        // int nalive=next_value();
        int nalive=util::get_rand_binomial(this->malive,pnotalive);

        // w[j]=nalive;
        step_expected_death = nalive;
        // cout<<"dead adult: "<<nalive<<endl;
        // village->malive=village->malive-nalive;
        // village->mdead=(nummosquito[j])-(village->malive);
        //      village->mdead=village->mdead+nalive;
        this->malive=this->malive-nalive;
        this->mdead=(this->nummosquito_i)-(this->malive);
    }

    this->step_expected_death_prob = float(step_expected_death) / float(prev_step_malive);
    this->step_expected_new_adult_prob = float(step_expected_new_adult) / float(this->mdead);

    this->step_actual_new_adult_counter = 0;
    this->step_actual_death_counter = 0;


    // if (kHost_village_id < 3) {
    //     std::cout << "v" << kHost_village_id;
    //     std::cout << ", Post-step:\n";
    //     std::cout << "eggs: " << this->eggs
    //                 << ", larva: " << this->larva
    //                 << ", pupa: " << this->pupa
    //                 << ", iadult: " << this->iadult
    //                 << ", malive: " << this->malive
    //                 << ", mdead: " << this->mdead << "\n";
    //     std::cin.ignore();
    // }

}

bool MosquitoSwampRa::step_inc_born_death(
        float K_i,
        std::vector<Dememosq*> G_Dememosquitoes
    ) {

    if (util::get_rand_uniform() < 0.5) {
        return this->step_inc_born_death_rev(K_i, G_Dememosquitoes);
    }

    const bool kFunc_verbose = false;

    bool if_max_larva_capacity = false;

    this->step_reset();

    float egg2larva_i=G_egg1*exp(-(G_egg2/K_i));
    float larva2iadult_i=G_larva1*exp(-(G_larva2/K_i));

    
    if (kFunc_verbose) {

        std::cout << "\nMosquitoSwampRa::step_inc_born_death >>";
        std::cout << "\nkHost_village_id=" << this->kHost_village_id << "\n";
        std::cout << "Kmax_i:" << this->Kmax_i << "\n";
        std::cout << "kG_Dememosquitoes_start_index:" << this->kG_Dememosquitoes_start_index << "\n";
        std::cout << "nummosquito_i:" << this->nummosquito_i << "\n";

        std::cout << "egg2larva_i:" << egg2larva_i << "\n";
        std::cout << "larva2iadult_i:" << larva2iadult_i << "\n";

        
        std::cout << "eggs:" << this->eggs << "\n";
        std::cout << "larvae:" << this->larva << "\n";
        std::cout << "iadult:" << this->iadult << "\n";
        std::cout << "malive:" << this->malive << "\n--\n";
        
    }


    if (this->eggs > 0) {
        float pnotegg=1-(1-egg2larva_i)*(1-G_aquatic_mortality);

        int neggs=util::get_rand_binomial(this->eggs,pnotegg);
        this->eggs=this->eggs-neggs;



        float normprobegg=egg2larva_i/(egg2larva_i+G_aquatic_mortality);

        int movedeggs=util::get_rand_binomial(neggs,normprobegg);

        int surveggstolarva=util::get_rand_binomial(movedeggs,(1-this->larva/Kmax_i));

        if (kFunc_verbose) {
            std::cout << "pnotegg:" << pnotegg << "\n";
            std::cout << "neggs:" << neggs << "\n";
            std::cout << "normprobegg:" << normprobegg << "\n";
            std::cout << "movedeggs:" << movedeggs << "\n";
            std::cout << "surveggstolarva:" << surveggstolarva << "\n--\n";
        }


        this->larva=this->larva+surveggstolarva;
        if (this->larva > Kmax_i) {
            this->larva = Kmax_i;
            if_max_larva_capacity = true;
        }
    }

    /** Larva to iadult **/
    if (this->larva > 0) {
        float pnotlarva=1-(1-G_aquatic_mortality)*(1-larva2iadult_i);

        int nlarva=util::get_rand_binomial(this->larva,pnotlarva);
        this->larva=this->larva-nlarva;

        float normproblarv=larva2iadult_i/(larva2iadult_i+G_aquatic_mortality);
        int movedlarva=util::get_rand_binomial(nlarva,normproblarv);

        if (kFunc_verbose) {
            std::cout << "pnotlarva:" << pnotlarva << "\n";
            std::cout << "nlarva:" << nlarva << "\n";
            std::cout << "normproblarv:" << normproblarv << "\n";
            std::cout << "movedlarva:" << movedlarva << "\n--\n";
        }

        this->iadult+=movedlarva;
        if (this->iadult > Kmax_i) {
            std::cout << "iadult:"<< this->iadult << " > Kmax_i:" << Kmax_i << "\n";
            // std::cin.ignore();
            this->iadult = Kmax_i;
        }
    }

    if (kFunc_verbose) {
        std::cout << "eggs:" << this->eggs << "\n";
        std::cout << "larvae:" << this->larva << "\n";
        std::cout << "iadult:" << this->iadult << "\n";
        std::cout << "malive:" << this->malive << "\n--\n";
    }



    const int kMosq_index_min = this->kG_Dememosquitoes_start_index; 
    const int kMosq_index_max = kMosq_index_min + this->nummosquito_i - 1;
    const int kMosq_slots_available = this->nummosquito_i;

    int prev_step_malive = this->malive;

    /*iadult to adult, new born*/
    if(this->iadult>0){
        float pnotimmature=1-(1-G_mature)*(1-G_iadult_life); 
        int nimmature=util::get_rand_binomial(this->iadult,pnotimmature);
        this->iadult=this->iadult-nimmature;

        float normprobiadult=G_mature/(G_mature+G_iadult_life);
        int movediadult=util::get_rand_binomial(nimmature, normprobiadult);

        step_expected_new_adult = movediadult;
        this->malive = this->malive + movediadult;
    }

    if (kFunc_verbose) {
        std::cout << "new born > " << "\n";

        std::cout << "this->larva=" << this->larva << "\n"; 
        std::cout << "this->iadult=" << this->iadult << "\n"; 
        std::cout << "prev_step_malive=" << prev_step_malive << "\n"; 
        std::cout << "step_expected_new_adult=" << step_expected_new_adult << "\n"; 
    }

    if(this->malive > kMosq_slots_available) {
        std::cout << "this->malive " << this->malive
                    << " > kMosq_slots_available" << kMosq_slots_available << std::endl;
        this->sum_num_mosq_slots_not_enough_incidents++;
        step_expected_new_adult -= (this->malive - kMosq_slots_available);
        this->malive -= (this->malive - kMosq_slots_available);
        // assert(this->malive <= kMosq_slots_available);
    }

    int num_turned_alive = 0;
    for (int mm = kMosq_index_min; mm <= kMosq_index_max; mm++) {
        if(num_turned_alive == step_expected_new_adult) {
            break;
        }
        if (G_Dememosquitoes[mm]->is_dead()) {
            G_Dememosquitoes[mm]->turn_alive();
            num_turned_alive++;
        }
    }
    assert(num_turned_alive == step_expected_new_adult);
    if (kFunc_verbose) {
        std::cout << "step_actual_new_adult_counter=" << step_actual_new_adult_counter << "\n--\n"; 
        std::cout << "death > " << "\n";
    }


    /*adult death*/
    if(this->malive>0){
        float pnotalive=1-(1-G_mosqdeath1);
        int nalive=util::get_rand_binomial(this->malive,pnotalive);
        
        if (kFunc_verbose) {
            std::cout << "this->malive=" << this->malive << "\n";
            std::cout << "pnotalive=" << pnotalive << "\n";
            std::cout << "nalive=" << nalive << std::endl;
        }

        step_expected_death = nalive;
        assert(step_expected_death <= this->malive);
        this->malive=this->malive-nalive;
        this->mdead=(this->nummosquito_i)-(this->malive);
        assert(this->mdead>=0);
    }


    int miss_fired = 0;
    for (int mm = 0; mm < step_expected_death; mm++) {
        int selected_slot_offset = util::get_rand_uniform() * kMosq_slots_available;
        if (G_Dememosquitoes[kMosq_index_min+selected_slot_offset]->is_alive()) {
            G_Dememosquitoes[kMosq_index_min+selected_slot_offset]->turn_dead();
        } else {
            mm--;
            miss_fired++;
        }
        if (miss_fired > kMosq_slots_available*10) {
            std::cout << "miss_fired > kMosq_slots_available*10" << std::endl;
            assert(false);
        }
    }

    if (kFunc_verbose) {
        std::cout << "step_expected_death=" << step_expected_death << "\n";
        std::cout << "step_actual_death_counter=" << step_actual_death_counter << "\n";
    }
    // this->step_expected_death_prob = float(step_expected_death) / float(prev_step_malive);
    // this->step_expected_new_adult_prob = float(step_expected_new_adult) / float(this->mdead);

    this->step_expected_death_prob = 0.0;
    this->step_expected_new_adult_prob = 0.0;


    return if_max_larva_capacity;

}

bool MosquitoSwampRa::step_inc_born_death_rev(
        float K_i,
        std::vector<Dememosq*> G_Dememosquitoes
    ) {

    const bool kFunc_verbose = false;

    bool if_max_larva_capacity = false;

    this->step_reset();

    float egg2larva_i=G_egg1*exp(-(G_egg2/K_i));
    float larva2iadult_i=G_larva1*exp(-(G_larva2/K_i));



    const int kMosq_index_min = this->kG_Dememosquitoes_start_index; 
    const int kMosq_index_max = kMosq_index_min + this->nummosquito_i - 1;
    const int kMosq_slots_available = this->nummosquito_i;

    
    if (kFunc_verbose) {

        std::cout << "\nMosquitoSwampRa::step_inc_born_death >>";
        std::cout << "\nkHost_village_id=" << this->kHost_village_id << "\n";
        std::cout << "Kmax_i:" << this->Kmax_i << "\n";
        std::cout << "kG_Dememosquitoes_start_index:" << this->kG_Dememosquitoes_start_index << "\n";
        std::cout << "nummosquito_i:" << this->nummosquito_i << "\n";

        std::cout << "egg2larva_i:" << egg2larva_i << "\n";
        std::cout << "larva2iadult_i:" << larva2iadult_i << "\n";

        
        std::cout << "eggs:" << this->eggs << "\n";
        std::cout << "larvae:" << this->larva << "\n";
        std::cout << "iadult:" << this->iadult << "\n";
        std::cout << "malive:" << this->malive << "\n--\n";
        
    }

    /*adult death*/
    if(this->malive>0){
        float pnotalive=1-(1-G_mosqdeath1);
        int nalive=util::get_rand_binomial(this->malive,pnotalive);
        
        if (kFunc_verbose) {
            std::cout << "death > " << "\n";
            std::cout << "this->malive=" << this->malive << "\n";
            std::cout << "pnotalive=" << pnotalive << "\n";
            std::cout << "nalive=" << nalive << "\n";
        }

        step_expected_death = nalive;
        assert(step_expected_death <= this->malive);
        this->malive=this->malive-nalive;
        this->mdead=(this->nummosquito_i)-(this->malive);
        assert(this->mdead>=0);
    }


    int miss_fired = 0;
    for (int mm = 0; mm < step_expected_death; mm++) {
        int selected_slot_offset = util::get_rand_uniform() * kMosq_slots_available;
        if (G_Dememosquitoes[kMosq_index_min+selected_slot_offset]->is_alive()) {
            G_Dememosquitoes[kMosq_index_min+selected_slot_offset]->turn_dead();
        } else {
            mm--;
            miss_fired++;
        }
        if (miss_fired > kMosq_slots_available*100) {
            std::cout << "miss_fired > kMosq_slots_available*100" << std::endl;
            assert(false);
        }
    }

    assert(step_expected_death == step_actual_death_counter);

    if (kFunc_verbose) {
        std::cout << "step_expected_death=" << step_expected_death << "\n";
        std::cout << "step_actual_death_counter=" << step_actual_death_counter << "\n";
    }
    // this->step_expected_death_prob = float(step_expected_death) / float(prev_step_malive);
    // this->step_expected_new_adult_prob = float(step_expected_new_adult) / float(this->mdead);

    this->step_expected_death_prob = 0.0;
    this->step_expected_new_adult_prob = 0.0;


    int prev_step_malive = this->malive;

    /*iadult to adult, new born*/
    if(this->iadult>0){
        float pnotimmature=1-(1-G_mature)*(1-G_iadult_life); 
        int nimmature=util::get_rand_binomial(this->iadult,pnotimmature);
        this->iadult=this->iadult-nimmature;

        float normprobiadult=G_mature/(G_mature+G_iadult_life);
        int movediadult=util::get_rand_binomial(nimmature, normprobiadult);

        step_expected_new_adult = movediadult;
        this->malive = this->malive + movediadult;
    }

    if (kFunc_verbose) {
        std::cout << "new born > " << "\n";

        std::cout << "this->larva=" << this->larva << "\n"; 
        std::cout << "this->iadult=" << this->iadult << "\n"; 
        std::cout << "prev_step_malive=" << prev_step_malive << "\n"; 
        std::cout << "step_expected_new_adult=" << step_expected_new_adult << "\n"; 
    }

    if(this->malive > kMosq_slots_available) {
        std::cout << "this->malive " << this->malive
                    << " > kMosq_slots_available" << kMosq_slots_available << std::endl;
        this->sum_num_mosq_slots_not_enough_incidents++;
        step_expected_new_adult -= (this->malive - kMosq_slots_available);
        this->malive -= (this->malive - kMosq_slots_available);
        // assert(this->malive <= kMosq_slots_available);
    }

    int num_turned_alive = 0;
    for (int mm = kMosq_index_min; mm <= kMosq_index_max; mm++) {
        if(num_turned_alive == step_expected_new_adult) {
            break;
        }
        if (G_Dememosquitoes[mm]->is_dead()) {
            G_Dememosquitoes[mm]->turn_alive();
            num_turned_alive++;
        }
    }
    assert(num_turned_alive == step_expected_new_adult);
    if (kFunc_verbose) {
        std::cout << "step_actual_new_adult_counter=" << step_actual_new_adult_counter << "\n--\n"; 
    }

    this->mdead=(this->nummosquito_i)-(this->malive);
    assert(this->mdead>=0);


    /** Larva to iadult **/
    if (this->larva > 0) {
        float pnotlarva=1-(1-G_aquatic_mortality)*(1-larva2iadult_i);

        int nlarva=util::get_rand_binomial(this->larva,pnotlarva);
        this->larva=this->larva-nlarva;

        float normproblarv=larva2iadult_i/(larva2iadult_i+G_aquatic_mortality);
        int movedlarva=util::get_rand_binomial(nlarva,normproblarv);

        if (kFunc_verbose) {
            std::cout << "pnotlarva:" << pnotlarva << "\n";
            std::cout << "nlarva:" << nlarva << "\n";
            std::cout << "normproblarv:" << normproblarv << "\n";
            std::cout << "movedlarva:" << movedlarva << "\n--\n";
        }

        this->iadult+=movedlarva;
        if (this->iadult > Kmax_i) {
            std::cout << "iadult:"<< this->iadult << " > Kmax_i:" << Kmax_i << "\n";
            // std::cin.ignore();
            this->iadult = Kmax_i;
        }
    }

    if (kFunc_verbose) {
        std::cout << "eggs:" << this->eggs << "\n";
        std::cout << "larvae:" << this->larva << "\n";
        std::cout << "iadult:" << this->iadult << "\n";
        std::cout << "malive:" << this->malive << "\n--\n";
    }



    if (this->eggs > 0) {
        float pnotegg=1-(1-egg2larva_i)*(1-G_aquatic_mortality);

        int neggs=util::get_rand_binomial(this->eggs,pnotegg);
        this->eggs=this->eggs-neggs;



        float normprobegg=egg2larva_i/(egg2larva_i+G_aquatic_mortality);

        int movedeggs=util::get_rand_binomial(neggs,normprobegg);

        int surveggstolarva=util::get_rand_binomial(movedeggs,(1-this->larva/Kmax_i));

        if (kFunc_verbose) {
            std::cout << "pnotegg:" << pnotegg << "\n";
            std::cout << "neggs:" << neggs << "\n";
            std::cout << "normprobegg:" << normprobegg << "\n";
            std::cout << "movedeggs:" << movedeggs << "\n";
            std::cout << "surveggstolarva:" << surveggstolarva << "\n--\n";
        }


        this->larva=this->larva+surveggstolarva;
        if (this->larva > Kmax_i) {
            this->larva = Kmax_i;
            if_max_larva_capacity = true;
        }
    }

    return if_max_larva_capacity;

}

bool MosquitoSwampRa::step_inc_born_death_ind(
        float K_i,
        std::vector<Dememosq*> G_Dememosquitoes
    ) {

    const bool kFunc_verbose = false;

    bool if_max_larva_capacity = false;

    this->step_reset();

    float egg2larva_i=G_egg1*exp(-(G_egg2/K_i));
    float larva2iadult_i=G_larva1*exp(-(G_larva2/K_i));



    const int kMosq_index_min = this->kG_Dememosquitoes_start_index; 
    const int kMosq_index_max = kMosq_index_min + this->nummosquito_i - 1;
    const int kMosq_slots_available = this->nummosquito_i;

    
    if (kFunc_verbose) {

        std::cout << "\nMosquitoSwampRa::step_inc_born_death_ind >>";
        std::cout << "\nkHost_village_id=" << this->kHost_village_id << "\n";
        std::cout << "Kmax_i:" << this->Kmax_i << "\n";
        std::cout << "kG_Dememosquitoes_start_index:" << this->kG_Dememosquitoes_start_index << "\n";
        std::cout << "nummosquito_i:" << this->nummosquito_i << "\n";

        std::cout << "egg2larva_i:" << egg2larva_i << "\n";
        std::cout << "larva2iadult_i:" << larva2iadult_i << "\n";

        
        std::cout << "eggs:" << this->eggs << "\n";
        std::cout << "larvae:" << this->larva << "\n";
        std::cout << "iadult:" << this->iadult << "\n";
        std::cout << "malive:" << this->malive << "\n--\n";
        
    }

    const int kPre_eggs = this->eggs;
    const int kPre_larva = this->larva;
    const int kPre_iadult = this->iadult;
    const int kPre_malive = this->malive;


    if (kPre_eggs > 0) {
        float pnotegg=1-(1-egg2larva_i)*(1-G_aquatic_mortality);

        int neggs=util::get_rand_binomial(kPre_eggs,pnotegg);
        this->eggs=this->eggs-neggs;



        float normprobegg=egg2larva_i/(egg2larva_i+G_aquatic_mortality);

        int movedeggs=util::get_rand_binomial(neggs,normprobegg);

        int surveggstolarva=util::get_rand_binomial(movedeggs,(1-kPre_larva/Kmax_i));

        if (kFunc_verbose) {
            std::cout << "pnotegg:" << pnotegg << "\n";
            std::cout << "neggs:" << neggs << "\n";
            std::cout << "normprobegg:" << normprobegg << "\n";
            std::cout << "movedeggs:" << movedeggs << "\n";
            std::cout << "surveggstolarva:" << surveggstolarva << "\n--\n";
        }


        this->larva=this->larva+surveggstolarva;
        if (this->larva > Kmax_i) {
            this->larva = Kmax_i;
            if_max_larva_capacity = true;
        }
    }


    /** Larva to iadult **/
    if (kPre_larva > 0) {
        float pnotlarva=1-(1-G_aquatic_mortality)*(1-larva2iadult_i);

        int nlarva=util::get_rand_binomial(kPre_larva,pnotlarva);
        this->larva=this->larva-nlarva;

        float normproblarv=larva2iadult_i/(larva2iadult_i+G_aquatic_mortality);
        int movedlarva=util::get_rand_binomial(nlarva,normproblarv);

        if (kFunc_verbose) {
            std::cout << "pnotlarva:" << pnotlarva << "\n";
            std::cout << "nlarva:" << nlarva << "\n";
            std::cout << "normproblarv:" << normproblarv << "\n";
            std::cout << "movedlarva:" << movedlarva << "\n--\n";
        }

        this->iadult+=movedlarva;
        if (this->iadult > Kmax_i) {
            std::cout << "iadult:"<< this->iadult << " > Kmax_i:" << Kmax_i << "\n";
            // std::cin.ignore();
            this->iadult = Kmax_i;
        }
    }

    if (kFunc_verbose) {
        std::cout << "eggs:" << this->eggs << "\n";
        std::cout << "larvae:" << this->larva << "\n";
        std::cout << "iadult:" << this->iadult << "\n";
        std::cout << "malive:" << this->malive << "\n--\n";
    }

    // int prev_step_malive = this->malive;



    /*iadult to adult, new born*/
    if(kPre_iadult>0){
        float pnotimmature=1-(1-G_mature)*(1-G_iadult_life); 
        int nimmature=util::get_rand_binomial(kPre_iadult,pnotimmature);
        this->iadult=this->iadult-nimmature;

        float normprobiadult=G_mature/(G_mature+0.75*G_iadult_life);
        // float normprobiadult=G_mature/(G_mature+G_iadult_life);
        int movediadult=util::get_rand_binomial(nimmature, normprobiadult);

        step_expected_new_adult = movediadult;
        this->malive = this->malive + movediadult;
    }

    if (kFunc_verbose) {
        std::cout << "new born > " << "\n";

        std::cout << "kPre_larva=" << kPre_larva << "\n"; 
        std::cout << "kPre_iadult=" << kPre_iadult << "\n"; 
        std::cout << "kPre_malive=" << kPre_malive << "\n"; 
        std::cout << "step_expected_new_adult=" << step_expected_new_adult << "\n"; 
    }

    if(this->malive > kMosq_slots_available) {
        std::cout << "this->malive " << this->malive
                    << " > kMosq_slots_available" << kMosq_slots_available << std::endl;
        this->sum_num_mosq_slots_not_enough_incidents++;
        step_expected_new_adult -= (this->malive - kMosq_slots_available);
        this->malive -= (this->malive - kMosq_slots_available);
        // assert(this->malive <= kMosq_slots_available);
    }

    int num_turned_alive = 0;
    for (int mm = kMosq_index_min; mm <= kMosq_index_max; mm++) {
        if(num_turned_alive == step_expected_new_adult) {
            break;
        }
        if (G_Dememosquitoes[mm]->is_dead()) {
            G_Dememosquitoes[mm]->turn_alive();
            num_turned_alive++;
        }
    }
    assert(num_turned_alive == step_expected_new_adult);
    if (kFunc_verbose) {
        std::cout << "step_actual_new_adult_counter=" << step_actual_new_adult_counter << "\n--\n"; 
    }

    this->mdead=(this->nummosquito_i)-(this->malive);
    assert(this->mdead>=0);





    /*adult death*/
    int num_over_aged_mosqs = 0;
    if(kPre_malive>0){

        for (int mm = kMosq_index_min; mm < kMosq_index_min+kMosq_slots_available; mm++) {
            if (G_Dememosquitoes[mm]->is_alive()
                && G_Dememosquitoes[mm]->get_age() >= 40) {
                if (G_Dememosquitoes[mm]->has_incubation()) {
                    this->step_death_infected++;
                }
                if (G_Dememosquitoes[mm]->has_infectious()) {
                    this->step_death_infectious++;
                }
                G_Dememosquitoes[mm]->turn_dead();
                num_over_aged_mosqs++;
            }
        }

        float pnotalive=1-(1-G_mosqdeath1);
        int nalive=util::get_rand_binomial(kPre_malive,pnotalive);

        if (nalive < num_over_aged_mosqs) {
            nalive = num_over_aged_mosqs;
        }
        
        if (kFunc_verbose) {
            std::cout << "death > " << "\n";
            std::cout << "kPre_malive=" << kPre_malive << "\n";
            std::cout << "pnotalive=" << pnotalive << "\n";
            std::cout << "nalive=" << nalive << "\n";
            std::cout << "kMosq_index_min=" << kMosq_index_min << "\n";
        }

        step_expected_death = nalive;
        assert(step_expected_death <= kPre_malive);
        assert(step_expected_death <= this->malive);
        this->malive -= nalive;
        this->mdead=(this->nummosquito_i)-(this->malive);
        assert(this->mdead >= 0);
    }



    int miss_fired = 0;
    // for (int mm = 0; mm < step_expected_death; mm++) {
    for (int mm = num_over_aged_mosqs; mm < step_expected_death; mm++) {
        int selected_slot_offset = util::get_rand_uniform() * kMosq_slots_available;
        int selected_index = kMosq_index_min+selected_slot_offset;
        // if (G_Dememosquitoes[selected_index]->is_alive()) {
        if (G_Dememosquitoes[selected_index]->is_alive()
            && G_Dememosquitoes[selected_index]->get_age()!=0
            ) {

            if (G_Dememosquitoes[selected_index]->has_incubation()) {
                this->step_death_infected++;
            }
            if (G_Dememosquitoes[selected_index]->has_infectious()) {
                this->step_death_infectious++;
            }
            G_Dememosquitoes[selected_index]->turn_dead();
        } else {
            mm--;
            miss_fired++;
        }
        if (miss_fired > kMosq_slots_available*100) {
            std::cout << "miss_fired > kMosq_slots_available*100" << std::endl;
            assert(false);
        }
    }

    if (kFunc_verbose) {
        std::cout << "step_expected_death=" << step_expected_death << "\n";
        std::cout << "step_actual_death_counter=" << step_actual_death_counter << "\n";
    }
    // if (num_over_aged_mosqs <= step_expected_death) {
        assert(step_expected_death == step_actual_death_counter);
    // } else {
    //     std::cout << "num_over_aged_mosqs("<< num_over_aged_mosqs <<") > step_expected_death("<< step_expected_death <<")\n";
    // }

    // this->step_expected_death_prob = float(step_expected_death) / float(prev_step_malive);
    // this->step_expected_new_adult_prob = float(step_expected_new_adult) / float(this->mdead);

    this->step_expected_death_prob = 0.0;
    this->step_expected_new_adult_prob = 0.0;


    return if_max_larva_capacity;

}

// void MosquitoSwampRa::step_population_alignment(
//         std::vector<Dememosq*> G_Dememosquitoes
//     ) {

//     assert(this->nummosquito_i > 0);

//     const int kMosq_index_min = this->kG_Dememosquitoes_start_index; 
//     const int kMosq_index_max = kMosq_index_min + this->nummosquito_i - 1;
//     const int kMosq_slots_available = this->nummosquito_i;

//     // const int kDeath_diff = this->step_expected_death - this->step_actual_death_counter;
//     // const int kNew_diff = this->step_expected_new_adult - this->step_actual_new_adult_counter;

//     if (this->step_actual_death_counter == this->step_expected_death ) {

//     } else if (this->step_actual_death_counter > this->step_expected_death) {

//         int num_alignments_death_to_alive = this->step_actual_death_counter - this->step_expected_death;
//         assert(num_alignments_death_to_alive <= kMosq_slots_available);

//         int num_turned_alive = 0;
//         for (int mm = kMosq_index_min; mm <= kMosq_index_max; mm++) {
//             if (G_Dememosquitoes[mm]->is_dead()) {
//                 G_Dememosquitoes[mm]->turn_alive();
//                 num_turned_alive++;
//                 if(num_turned_alive == num_alignments_death_to_alive) {
//                     break;
//                 }
//             }
//         }

//         assert(num_turned_alive == num_alignments_death_to_alive);

//     } else { //this->step_actual_death_counter < this->step_expected_death

//         int num_alignments_alive_to_death = this->step_expected_death - this->step_actual_death_counter;


//     }


//         for (int mm = 0; mm < num_alignments_death_to_alive; mm++) {
//             int selected_slot = util::get_rand_uniform() * kMosq_slots_available;
//             if (G_Dememosquitoes[kMosq_index_min + selected_slot]->is_dead()) {
//                 G_Dememosquitoes[kMosq_index_min + selected_slot]->turn_alive();
//             } else {
//                 mm--;
//             }
//         }

// }

MosquitoManagerIbmRa::MosquitoManagerIbmRa(
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
    ) : 
        kResting_to_biting_rate(1.0/resting_days),
        kBiting_rate(biting_rate),
        G_mosqinc(1.0/incubation_days),
        G_eggslayed(eggs_per_lay),
        kG_carryingcapacity(g_carryingcapacity),
        kScaletomosqnumbers(scaletomosqnumbers),
        kFunc_swamp_step_inc_born_death_version(func_swamp_step_inc_born_death_version)
        {

    // const float kInitProbAlive = 0.5;
    // const float kInitProbPregnant = 0.1;
    // const float kInitProbBiting = 0.2;

    // const float kInitProbInfectious = 0.05;

    this->mosquito_swamp_of.clear();
    this->G_Dememosquitoes.clear();

    int total_nummosquito = 0;
    float total_Kmax = 0.0;

    for (int vv = 0; vv < num_villages; vv++) {

        // float beta_values=betaf*vregion[i]->v[3];
        float beta_values=(this->betaf)*transmission_coefficient_of[vv];
        // int mpp = float(beta_values*(this->kScaletomosqnumbers));
        float mpp = float(beta_values*(this->kScaletomosqnumbers));

        // line 287
        // nummosquito[i]=int(float(number)*mpp);
        int vv_nummosquito_i=int(float(population_of[vv])*mpp);
        // Kmax[i]=nummosquito[i]*G_carryingcapacity;
        float vv_Kmax_i=vv_nummosquito_i*(this->kG_carryingcapacity);

    //     if (vv == 0) {
    // std::cout << "for vv == 0:\n";
    // std::cout << "vv_nummosquito_i = " << vv_nummosquito_i;
    // std::cout << ", vv_Kmax_i = " << vv_Kmax_i << "\n";
    //     }

        total_nummosquito += vv_nummosquito_i;
        total_Kmax += vv_Kmax_i;

        this->mosquito_swamp_of.push_back(
            new MosquitoSwampRa(
                aquatic_mortality,
                maturation_days,
                immat_life_days,
                adult_life_days,
                vv,
                vv_nummosquito_i,
                static_cast<int>(this->G_Dememosquitoes.size()),
                vv_Kmax_i
            )
        );

        MosquitoSwampRa* current_swamp_ptr = this->mosquito_swamp_of.back();

        // // line 762, fitnessinit, aguasSIRk.h
        for (int mm = 0; mm < vv_nummosquito_i; mm++) {
            this->G_Dememosquitoes.push_back(
                new Dememosq(
                    current_swamp_ptr,
                    0
                )
            );

            Dememosq* current_mosq_ptr = this->G_Dememosquitoes.back();

            if (util::get_rand_uniform() < this->kInitProbAlive) {
                current_mosq_ptr->set_alive(true);
                current_swamp_ptr->inc_alive();

                if (util::get_rand_uniform() < this->kInitProbPregnant) {
                    current_mosq_ptr->set_pregnant(true);
                }

                if (util::get_rand_uniform() < this->kInitProbBiting) {
                    current_mosq_ptr->set_biting(true);
                }

                if (util::get_rand_uniform() < this->kInitProbInfectious) {
                    current_mosq_ptr->add_infectious(common::ParasiteType::kR0);
                }

                //// not implemented
                // float randomEvent4 = ran2();
                // if(randomEvent4 <= 0.01){G_Dememosquitoes[mm]->resistance=RESISTRA;}
                // else{G_Dememosquitoes[mm]->resistance=RESISTRO;}

            } else {
                current_mosq_ptr->set_alive(false);
                current_swamp_ptr->inc_dead();
            }
        }
    }

    std::cout << "after init:\n";

    std::cout << "total_nummosquito = " << total_nummosquito;
    std::cout << ", total_Kmax = " << total_Kmax << "\n";

    this->step_print_stats();
    // std::cin.ignore();

}

MosquitoManagerIbmRa::~MosquitoManagerIbmRa() {
    for (uint vv = 0; vv < this->mosquito_swamp_of.size(); vv++) {
        delete this->mosquito_swamp_of[vv];
    }
    for (uint mm = 0; mm < this->G_Dememosquitoes.size(); mm++) {
        delete this->G_Dememosquitoes[mm];
    }
}

void MosquitoManagerIbmRa::step(
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

    ) { 

    // static constexpr double kNum_days_per_year = 365.0;
    // const double kPi = 3.14;

    assert(infections_array_size == host_array_size);
    double t_double = static_cast<double>(timestep);

    // K[i]=273.15+25.0+5.0*cos(timestep/365.0*2.0*3.1416);
    // float K_i=273.15+25.0+0.8*5.0*cos(timestep/365.0*2.0*3.1416);
    // float K_i=273.15+25.0+0.8*5.0*cos((timestep+2*kNum_days_per_year)/365.0*2.0*kPi);

    // float K_i=273.15+25.0+0.8*5.0*cos(timestep/kNum_days_per_year*2.0*kPi);
    float K_i=273.15+25.0+0.8*5.0*cos(t_double/common::kDays_per_year*2.0*common::kPI);

    int total_eggs = 0;
    int total_larva = 0;
    int total_pupa = 0;
    int total_iadult = 0;
    int total_malive = 0;
    int total_mdead = 0;


    // if (timestep % 10 == 0) {
    //     std::cout << "IBM timestep=" << timestep << ", K_i=" << K_i << "\n";
    //     std::cin.ignore();
    // }

    // 2119, villages; 1307-1344 mosquitodynamicsnew;
    for (MosquitoSwampRa* swamp : this->mosquito_swamp_of) {
        swamp->step(K_i);

        // total_eggs += swamp->get_status_eggs();
        // total_larva += swamp->get_status_larva();
        // total_pupa += swamp->get_status_pupa();
        // total_iadult += swamp->get_status_iadult();
        // total_malive += swamp->get_status_malive();
        // total_mdead += swamp->get_status_mdead();

    }

    // this->sum_num_eggs = total_eggs;
    // this->sum_num_larv = total_larva;
    // this->sum_num_imat = total_iadult;

    // std::cout << "eggs: " << total_eggs
    //             << ", larva: " << total_larva
    //             << ", pupa: " << total_pupa
    //             << ", iadult: " << total_iadult
    //             << ", malive: " << total_malive
    //             << ", mdead: " << total_mdead << "\n";
    // std::cin.ignore();



    // this->step_print_stats();


    // float incubation=G_mosqinc*exp(-(G_mosqinc/K));
    // float incubation_prob = this->G_mosqinc*exp(-(this->G_mosqinc/K_i));
    // float incubation_prob = this->G_mosqinc;
    // std::cout << "incubation_prob:" << incubation_prob << "\n";

    // double infected_to_infectious = kInfected_arrhenius_1*exp(-kInfected_arrhenius_2/Tk);
    // float incubation_prob = 1.17E+11*exp(-8340/K_i);
    // float incubation_prob = (1.0/14.0)*exp(-(1.0/14.0)/K_i);

    float incubation_prob = 1.0/14.0;


    this->sum_num_infected_biting = 0;
    this->sum_num_infectious_biting = 0;

    this->sum_num_infectious_bites = 0;
    this->sum_num_successful_infectious_bites = 0;

    int kia = 0;

    // 1387, mosquitodynamicsnew
    for (Dememosq* mosq: this-> G_Dememosquitoes) {

        if (mosq->is_alive()) {
            // alive

            mosq->inc_age();

            if (util::get_rand_uniform() < mosq->get_home_swamp()->get_step_expected_death_prob()) {
                mosq->turn_dead();
                continue;
            } 

            bool killed_in_action = false;

            if (mosq->is_biting()) { 

                if (mosq->has_incubation()) {
                    this->sum_num_infected_biting++;
                }
                if (mosq->has_infectious()) {
                    this->sum_num_infectious_biting++;
                }

                int num_new_bites = util::get_rand_poisson(this->kBiting_rate);
                for (int bb = 0; bb < num_new_bites; bb++) { 

                    float prob_transmission = 1.0;

                    // dies while trying to bite
                    if ( !( mosq->get_home_swamp()->step_if_expected_death_num_reached() )
                        && ( util::get_rand_uniform() < this->G_deathbite )
                    ) {
                        killed_in_action = true;
                        prob_transmission = this->G_deathbite_transmission;

                    // } else { // see issue #4, https://github.com/ATOME-MORU/Mosquito-model/issues/4
                    }

                    if( (prob_transmission == 1.0)
                        || (util::get_rand_uniform() < prob_transmission)
                    ) {

                        int human_id = vll_mgr.get_random_human_id_at_village(
                            mosq->get_home_swamp()->get_host_village_id()
                        );

                        assert(human_id < infections_array_size);

                        // h2m
                        if (dominant_type_array[human_id] != common::ParasiteType::kNoParasite) {    
                            if (util::get_rand_uniform() < (infectiousness_multiplier * infectiousness_array[human_id])) {
                                mosq->add_incubation(
                                    dominant_type_array[human_id],
                                    this->get_swamp_of_village(human_home_village_array[human_id])
                                );
                            }
                        }

                        // std::cout << " mosq->select_rnd_infectious().pt=" <<  mosq->select_rnd_infectious().pt << std::endl;
                        // m2h
                        if (mosq->has_infectious()) {

                            this->sum_num_infectious_bites++;

                            if (most_advanced_stage_array[human_id] == common::StageName::kInfectious) {
                                continue;
                            }

                            if (util::get_rand_uniform() < (susceptibility_multiplier * susceptibility_array[human_id])) {
                                Dememosq::Infection infection = mosq->select_rnd_infectious();
                                infections_array[human_id] = infection.pt;

                                mosq->inc_infectious_count();
                                this->sum_num_successful_infectious_bites++;
                            }
                        }

                        if (killed_in_action) {
                            kia++;
                            mosq->turn_dead(); 
                            break;
                        }
                    }
                    // }
                } // bb


                if (!killed_in_action) {

                    // incubation
                    static_cast<void>(mosq->step_incubate(incubation_prob));

                    // // pregnancy / lay eggs
                    // mosq->step_pregnancy(
                    //     this->G_mpregnancy,
                    //     this->G_mgestation,
                    //     this->G_eggslayed
                    // );

                    // biting / resting switch
                    mosq->step_biting_resting_switch(
                        this->kBiting_resting_prob
                    );
                    
                }


            // mosq->is_biting()
            } else {

                // incubation
                static_cast<void>(mosq->step_incubate(incubation_prob));

                // pregnancy / lay eggs
                mosq->step_pregnancy(
                    this->G_mpregnancy,
                    this->G_mgestation,
                    this->G_eggslayed
                );

                // biting / resting switch
                mosq->step_biting_resting_switch(
                    this->kBiting_resting_prob
                );

            }



        } else {

            // not alive
            if (util::get_rand_uniform() < mosq->get_home_swamp()->get_step_expected_new_adult_prob()) {
                mosq->turn_alive();
                // assert(mosq->get_home_swamp()->has_iadault());
            }

        }

    }//  G_Dememosquitoes


    for (MosquitoSwampRa* swamp : this->mosquito_swamp_of) {

        total_eggs += swamp->get_status_eggs();
        total_larva += swamp->get_status_larva();
        total_pupa += swamp->get_status_pupa();
        total_iadult += swamp->get_status_iadult();
        total_malive += swamp->get_status_malive();
        total_mdead += swamp->get_status_mdead();

    }

    this->sum_num_eggs = total_eggs;
    this->sum_num_larv = total_larva;
    this->sum_num_imat = total_iadult;

    std::cout << "eggs: " << total_eggs
                << ", larva: " << total_larva
                << ", pupa: " << total_pupa
                << ", iadult: " << total_iadult
                << ", malive: " << total_malive
                << ", mdead: " << total_mdead << "\n";


    // std::cout << "kia=" << kia << std::endl;
    this->step_print_stats();
}

void MosquitoManagerIbmRa::step_matching_ode(
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

    ) { 

    // static constexpr double kNum_days_per_year = 365.0;
    // const double kPi = 3.14;
    const bool kFunc_verbose = false;

    double t_double = static_cast<double>(timestep);

    assert(infections_array_size == host_array_size);

    // K[i]=273.15+25.0+5.0*cos(timestep/365.0*2.0*3.1416);
    // float K_i=273.15+25.0+0.8*5.0*cos(timestep/365.0*2.0*3.1416);
    // float K_i=273.15+25.0+0.8*5.0*cos((timestep+2*kNum_days_per_year)/365.0*2.0*kPi);

    // float K_i=273.15+25.0+0.8*5.0*cos(timestep/kNum_days_per_year*2.0*kPi);
    float K_i=273.15+25.0+0.8*5.0*cos(t_double/common::kDays_per_year*2.0*common::kPI);


    int num_mosq_layed_eggs = 0;


    // if (timestep % 10 == 0) {
    //     std::cout << "IBM timestep=" << timestep << ", K_i=" << K_i << "\n";
    //     std::cin.ignore();
    // }

    this->sum_num_larva_capacity = 0;

    this->sum_num_death_infected = 0;
    this->sum_num_death_infectious = 0;

    for (MosquitoSwampRa* swamp : this->mosquito_swamp_of) {

        if (kFunc_swamp_step_inc_born_death_version == 0) {
            if( swamp->step_inc_born_death(
                    K_i,
                    this->G_Dememosquitoes
                )
            ) {
                this->sum_num_larva_capacity++;
            }
        } else if (kFunc_swamp_step_inc_born_death_version == 1) {        
            if( swamp->step_inc_born_death_ind(
                    K_i,
                    this->G_Dememosquitoes
                )
            ) {
                this->sum_num_larva_capacity++;
            }
        } else {
            std::cout << "Unknown kFunc_swamp_step_inc_born_death_version: " <<  kFunc_swamp_step_inc_born_death_version << "\n";
            assert(false);
        }

        this->sum_num_death_infected += swamp->get_step_death_infected();
        this->sum_num_death_infectious += swamp->get_step_death_infectious();

    }

    // float incubation_prob = 1.0/14.0;


    this->sum_num_infected_biting = 0;
    this->sum_num_infectious_biting = 0;
    this->sum_num_clean_biting = 0;

    this->sum_num_new_infected = 0;
    this->sum_num_incubation = 0;


    this->sum_num_infectious_bites = 0;
    this->sum_num_successful_infectious_bites = 0;

    this->sum_num_superinfection_requests = 0;

    int kia = 0;

    int num_mosq_alive = 0;
    int num_mosq_alive_biting = 0;


    // 1387, mosquitodynamicsnew
    for (Dememosq* mosq: this-> G_Dememosquitoes) {

        if (mosq->is_alive()) {
            // alive
            num_mosq_alive++;

            mosq->inc_age();

            // if (!( mosq->get_home_swamp()->step_if_expected_death_num_reached())) {
            //     if (util::get_rand_uniform() < mosq->get_home_swamp()->get_step_expected_death_prob()) {
            //         mosq->turn_dead();
            //         continue;
            //     }
            // }

                // mosq->step_incubate(incubation_prob);


            bool incubation_done = false;
            if (util::get_rand_uniform() < 0.5) {
                // this->sum_num_incubation += mosq->step_incubate(this->G_mosqinc);
                this->sum_num_incubation += mosq->step_incubate(this->G_mosqinc);
                incubation_done = true;
            }




            bool killed_in_action = false;

            if (mosq->is_biting()) {
                num_mosq_alive_biting++;

                if (mosq->has_incubation()) {
                    this->sum_num_infected_biting++;
                }
                if (mosq->has_infectious()) {
                    this->sum_num_infectious_biting++;
                }
                if ((!mosq->has_incubation()) && (!mosq->has_infectious())) {
                    this->sum_num_clean_biting++;
                }

                // int num_new_bites = util::get_rand_poisson(this->kBiting_rate);
                // for (int bb = 0; bb < num_new_bites; bb++) { 

                if( util::get_rand_uniform() < this->kBiting_rate) {

                    float prob_transmission = 1.0;

                    // // dies while trying to bite
                    // if ( !( mosq->get_home_swamp()->step_if_expected_death_num_reached() )
                    //     && ( util::get_rand_uniform() < this->G_deathbite )
                    // ) {
                    //     killed_in_action = true;
                    //     prob_transmission = this->G_deathbite_transmission;
                    // }

                    if( (prob_transmission == 1.0)
                        || (util::get_rand_uniform() < prob_transmission)
                    ) {

                        int human_id = vll_mgr.get_random_human_id_at_village(
                            mosq->get_home_swamp()->get_host_village_id()
                        );

                        assert(human_id < infections_array_size);

                        // h2m
                        if (dominant_type_array[human_id] != common::ParasiteType::kNoParasite) {                        
                            if ((mosq->has_incubation()) || (mosq->has_infectious())) {
                                this->sum_num_superinfection_requests++;
                            }
                        }

                        // if ((!mosq->has_incubation()) && (!mosq->has_infectious())) { // superinfection filter

                            if (dominant_type_array[human_id] != common::ParasiteType::kNoParasite) {    
                                if (util::get_rand_uniform() < (infectiousness_multiplier * infectiousness_array[human_id])) {
                                    mosq->add_incubation(
                                        dominant_type_array[human_id],
                                        this->get_swamp_of_village(human_home_village_array[human_id])
                                    );
                                    this->sum_num_new_infected++;
                                }
                            }
                        // }


                        // std::cout << " mosq->select_rnd_infectious().pt=" <<  mosq->select_rnd_infectious().pt << std::endl;
                        // m2h
                        if (mosq->has_infectious()) {

                            this->sum_num_infectious_bites++;

                            if (most_advanced_stage_array[human_id] == common::StageName::kInfectious) {
                                // wasted infectious bite
                            } else {

                                if (util::get_rand_uniform() < (susceptibility_multiplier * susceptibility_array[human_id])) {
                                    Dememosq::Infection infection = mosq->select_rnd_infectious();
                                    infections_array[human_id] = infection.pt;

                                    mosq->inc_infectious_count();

                                    this->sum_num_successful_infectious_bites++;
                                }
                            }

                        }

                    }

                    mosq->set_biting(false);

                }// kBiting_rate




            } else { // not biting

                mosq->get_home_swamp()->inc_egg(this->G_eggslayed*this->kResting_to_biting_rate);
                // mosq->get_home_swamp()->inc_egg(this->G_eggslayed);
                num_mosq_layed_eggs += 1;

                if( util::get_rand_uniform() < this->kResting_to_biting_rate) {
                    mosq->set_biting(true);
                }

    // std::cout << "kResting_to_biting_rate: " << this->kResting_to_biting_rate << "\n";
    // std::cin.ignore();


            }

            if (killed_in_action) {
                kia++;
                mosq->turn_dead();
                assert(!killed_in_action);

            } else {
                
                // incubation
                // this->sum_num_incubation += mosq->step_incubate(incubation_prob);
                if (!incubation_done) {
                    this->sum_num_incubation += mosq->step_incubate(this->G_mosqinc);
                }
            }


        } else {

            // // not alive
            // if (util::get_rand_uniform() < mosq->get_home_swamp()->get_step_expected_new_adult_prob()) {
            //     mosq->turn_alive();
            //     // assert(mosq->get_home_swamp()->has_iadault());
            // }

        }

    }//  G_Dememosquitoes


    int total_eggs = 0;
    int total_larva = 0;
    // int total_pupa = 0;
    int total_iadult = 0;
    int total_malive = 0;
    int total_mdead = 0;

    int total_new_adult = 0;
    int total_death_adult = 0;

    int total_death_adult_infectious_count = 0;

    int total_incidents = 0;

    for (MosquitoSwampRa* swamp : this->mosquito_swamp_of) {

        total_eggs += swamp->get_status_eggs();
        total_larva += swamp->get_status_larva();
        // total_pupa += swamp->get_status_pupa();
        total_iadult += swamp->get_status_iadult();
        total_malive += swamp->get_status_malive();
        total_mdead += swamp->get_status_mdead();

        total_new_adult +=swamp->get_step_expected_new_adult();
        total_death_adult +=swamp->get_step_expected_death();

        total_death_adult_infectious_count += swamp->get_step_actual_death_infectious_counter();

        total_incidents += swamp->sum_num_mosq_slots_not_enough_incidents;

    }


    if (num_mosq_alive != total_malive) {
        std::cout << "num_mosq_alive=" << num_mosq_alive << "\n";
        std::cout << "total_malive=" << total_malive << "\n";
        assert(num_mosq_alive == total_malive);
    }

    this->sum_num_eggs = total_eggs;
    this->sum_num_larv = total_larva;
    this->sum_num_imat = total_iadult;
    this->sum_num_malive = total_malive;
    this->sum_num_malive_biting = num_mosq_alive_biting;

    this->sum_num_new_adults = total_new_adult;
    this->sum_num_death_adults = total_death_adult;

    this->sum_num_death_adults_infectious_count = total_death_adult_infectious_count;


    this->sum_num_new_eggs = num_mosq_layed_eggs;

    if (total_incidents > 0) {
        std::cout << "\033[34;47m[mIBM incidents]\033[0m" << " : " << total_incidents << "\n";
    } else {
        if (kFunc_verbose) {
            std::cout << "\033[32;47m[mIBM no incidents]\033[0m" << " : " << total_incidents << "\n";
        }
    }

    if (kFunc_verbose) {
        std::cout << "\nMosquitoManagerIbmRa::step_matching_ode >>";

        std::cout << "\neggs: " << total_eggs
                    << ", larva: " << total_larva
                    // << ", pupa: " << total_pupa
                    << ", iadult: " << total_iadult
                    << ", malive: " << total_malive
                    << ", mdead: " << total_mdead
                    << ", new_egg_layer:" << num_mosq_layed_eggs << "\n";
        this->step_print_stats();
    }

}

std::vector<int> MosquitoManagerIbmRa::get_mosq_pop_of_swamp(int swamp_id) {
    std::vector<int> swamp_pops;
    swamp_pops.push_back(this->mosquito_swamp_of.at(swamp_id)->get_status_eggs());
    swamp_pops.push_back(this->mosquito_swamp_of.at(swamp_id)->get_status_larva());
    swamp_pops.push_back(this->mosquito_swamp_of.at(swamp_id)->get_status_pupa());
    swamp_pops.push_back(this->mosquito_swamp_of.at(swamp_id)->get_status_iadult());
    swamp_pops.push_back(this->mosquito_swamp_of.at(swamp_id)->get_status_malive());
    swamp_pops.push_back(this->mosquito_swamp_of.at(swamp_id)->get_status_mdead());

    return swamp_pops;
}

void MosquitoManagerIbmRa::step_print_stats() const {
    // this->sum_num_infectious_mosq = 0;

     int tMOSQ = 0;
     int tTOT = 0;
     int deadMOSQ = 0;
     int Mosqinf = 0;
     int Mosqicb = 0;
     int Mosqseek = 0;
     int mfailinfectious=0;

     int tm=0;

     int biting_inf = 0;
     int biting_inc = 0;
     int biting_both = 0;

    for (Dememosq* mosq: this-> G_Dememosquitoes) {

        tm+=1;

        if (mosq->is_biting()){
            Mosqseek+=1;
            assert(mosq->is_alive());
            if(mosq->has_infectious()) {
                biting_inf +=1;
            }
            if(mosq->has_incubation()) {
                biting_inc +=1;
            }
            if(mosq->has_infectious() && mosq->has_incubation()) {
                biting_both +=1;
            }
        }

        if (mosq->is_alive()){
            tMOSQ+=1;
            tTOT+=1;
        }
        if (!(mosq->is_alive())){
            deadMOSQ+=1;
            tTOT+=1;
            if (mosq->has_infectious()){
                mfailinfectious+=1;
            }
        }
        if (mosq->has_infectious()){
            Mosqinf+=1;
        }
        if (mosq->has_incubation()){
            Mosqicb+=1;
        }
    }

    std::cout << "aliveMOSQ:" << tMOSQ  << "\tdeadMosqs:" << deadMOSQ <<"\ttotalMosqs:" << tTOT << std::endl
                << "Mosqinf:" << Mosqinf << "\tMosqicb:" << Mosqicb
                << "Mosqseek:" << Mosqseek << "\tbiting_inf:" << biting_inf << "\tbiting_inc:" << biting_inc << "\tbiting_both:" << biting_both << std::endl;
    std::cout << "MOSQ FAIL: "<<mfailinfectious << "\ttm: " << tm <<std::endl;

}

}