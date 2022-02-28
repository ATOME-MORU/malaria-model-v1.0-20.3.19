#include <iostream>

#include "village/mosquito_ode.cuh"

namespace village {

MosquitoEckhoff_multi::Fx::Fx(
        MosquitoEckhoff_multi::value_t larva_to_immature,
        MosquitoEckhoff_multi::value_t infected_to_infectious,
        MosquitoEckhoff_multi::value_t eggs_to_larva
    ) : 
        kLarva_to_immature(larva_to_immature),
        kInfected_to_infectious(infected_to_infectious),
        kEggs_to_larva(eggs_to_larva) {

}

MosquitoEckhoff_multi::MosquitoEckhoff_multi(
        std::vector<float>& host_max_larval_capacity_list
    ) : 
        kNum_mosqs(host_max_larval_capacity_list.size()),
        kMax_larval_capacity_list(host_max_larval_capacity_list),
        population(9*kNum_mosqs) {

    thrust::fill(population.begin(), population.end(), 0.0);
    thrust::fill(
        population.begin() + 2 * kNum_mosqs,
        population.begin() + 3 * kNum_mosqs,
        200.0
    );

    this->reset_time();
}

void MosquitoEckhoff_multi::reset_time() {
    this->current_time = 0.0;
}

void MosquitoEckhoff_multi::pre_integration_update(
        std::vector<float> host_infection_prob_list
    ) {

    assert(host_infection_prob_list.size() == this->kNum_mosqs);
    this->infection_prob_list = host_infection_prob_list;

    this->current_temperature =
        273.15 + 25.0 + 0.8 * 5.0 * cos(
                                (this->current_time-60)
                                / this->kNum_days_per_year
                                * 2.0 * this->kPi
                            );
    this->current_larva_to_immature =
        this->kAquatic_arrhenius_1
        * exp(-this->kAquatic_arrhenius_2 / this->current_temperature);
    this->current_infected_to_infectious = 
        this->kInfected_arrhenius_1
        * exp(-this->kInfected_arrhenius_2 / this->current_temperature);
    this->current_eggs_to_larva = 
        this->kEgg_arrhenius1
        * exp(-this->kEgg_arrhenius2 / this->current_temperature);


    this->current_time += this->kTime_length_per_integration;
}

void MosquitoEckhoff_multi::operator()(state_t& x, state_t& dxdt, value_t t) const {
    (void) t;
    thrust::for_each(
        thrust::make_zip_iterator( thrust::make_tuple(
            this->kMax_larval_capacity_list.begin(),
            this->infection_prob_list.begin(),
            thrust::make_zip_iterator( thrust::make_tuple(
                x.begin(),                  //[0]
                x.begin() + 1 * kNum_mosqs, //[1]
                x.begin() + 2 * kNum_mosqs, //[2]
                x.begin() + 3 * kNum_mosqs, //[3]
                x.begin() + 4 * kNum_mosqs, //[4]
                x.begin() + 5 * kNum_mosqs, //[5]
                x.begin() + 6 * kNum_mosqs, //[6]
                x.begin() + 7 * kNum_mosqs, //[7]
                x.begin() + 8 * kNum_mosqs  //[8]
            ) ),
            thrust::make_zip_iterator( thrust::make_tuple(
                dxdt.begin(),                  //[0]
                dxdt.begin() + 1 * kNum_mosqs, //[1]
                dxdt.begin() + 2 * kNum_mosqs, //[2]
                dxdt.begin() + 3 * kNum_mosqs, //[3]
                dxdt.begin() + 4 * kNum_mosqs, //[4]
                dxdt.begin() + 5 * kNum_mosqs, //[5]
                dxdt.begin() + 6 * kNum_mosqs, //[6]
                dxdt.begin() + 7 * kNum_mosqs, //[7]
                dxdt.begin() + 8 * kNum_mosqs  //[8]
            ) )
        ) ),
        thrust::make_zip_iterator( thrust::make_tuple(
            this->kMax_larval_capacity_list.end(),
            this->infection_prob_list.end(),
            thrust::make_zip_iterator( thrust::make_tuple(
                x.begin() + 1 * kNum_mosqs, //[0]
                x.begin() + 2 * kNum_mosqs, //[1]
                x.begin() + 3 * kNum_mosqs, //[2]
                x.begin() + 4 * kNum_mosqs, //[3]
                x.begin() + 5 * kNum_mosqs, //[4]
                x.begin() + 6 * kNum_mosqs, //[5]
                x.begin() + 7 * kNum_mosqs, //[6]
                x.begin() + 8 * kNum_mosqs, //[7]
                x.begin() + 9 * kNum_mosqs //[8]
            ) ),
            thrust::make_zip_iterator( thrust::make_tuple(
                dxdt.begin() + 1 * kNum_mosqs, //[0]
                dxdt.begin() + 2 * kNum_mosqs, //[1]
                dxdt.begin() + 3 * kNum_mosqs, //[2]
                dxdt.begin() + 4 * kNum_mosqs, //[3]
                dxdt.begin() + 5 * kNum_mosqs, //[4]
                dxdt.begin() + 6 * kNum_mosqs, //[5]
                dxdt.begin() + 7 * kNum_mosqs, //[6]
                dxdt.begin() + 8 * kNum_mosqs, //[7]
                dxdt.begin() + 9 * kNum_mosqs //[8]
            ) )
        ) ),
        Fx(
            this->current_larva_to_immature,
            this->current_infected_to_infectious,
            this->current_eggs_to_larva
        )
    );
}

MosquitoManagerODEThrust::MosquitoManagerODEThrust(
        std::vector<float> max_larval_capacity_list
    ) {
    this->mosq_sys = new MosquitoEckhoff_multi(max_larval_capacity_list);

}

MosquitoManagerODEThrust::~MosquitoManagerODEThrust() {
    delete this->mosq_sys;
}

void MosquitoManagerODEThrust::step(
        const std::vector<float>& infection_prob_list
    ) {

    this->mosq_sys->pre_integration_update(infection_prob_list);

    odeint::integrate_const(
        this->ode_stepper,
        *(this->mosq_sys),
        this->mosq_sys->get_population(),
        0.0, 1.0, 0.01
    );
}

}