#ifndef MOSQUITO_ODE_CUH
#define MOSQUITO_ODE_CUH

#include <thrust/device_vector.h>
#include <thrust/tuple.h>
#include <thrust/for_each.h>
#include <thrust/iterator/zip_iterator.h>

#include <boost/numeric/odeint/external/thrust/thrust.hpp>

#include <boost/numeric/odeint.hpp>

namespace odeint = boost::numeric::odeint;

namespace village{

class MosquitoEckhoff_multi {

public:
    typedef double value_t;
    typedef thrust::device_vector<value_t> state_t;
    typedef odeint::runge_kutta4<state_t, value_t, state_t, value_t, odeint::thrust_algebra, odeint::thrust_operations> stepper_type;
    // typedef odeint::runge_kutta4<state_t, value_t, state_t, value_t> stepper_type;



private:

    const value_t kTime_length_per_integration = 1.0;

    static constexpr value_t kNum_days_per_year = 365.0;
    const value_t kPi = 3.14159265358979;
    
    const value_t kAquatic_arrhenius_1 = 8.42E+10;
    const value_t kAquatic_arrhenius_2 = 8328.0;
    const value_t kInfected_arrhenius_1 = 1.17E+11;
    const value_t kInfected_arrhenius_2 = 8340.0;
    const value_t kEgg_arrhenius1 = 6.16E+07;
    const value_t kEgg_arrhenius2 = 5754.03;

public:
    class Fx {
    public:

        const value_t kAquatic_mortality = 0.1;
        const value_t kImmature_adult = 1/2.0;
        const value_t kAdult_life = 1/10.0;
        const value_t kHuman_feeding_mort = 0.1;
        const value_t kInfect_human_feed = 1.5; 
        const value_t kInfected_egg_lay = 0.8;
        const value_t kFeed_rate = 1/2.0;
        const value_t kLay_rate = 1/1.0;

        // const double oviposition=50.0;
        const value_t kOviposition = 4.0;

        const value_t kLarva_to_immature;
        const value_t kInfected_to_infectious;
        const value_t kEggs_to_larva;

        Fx(
            value_t larva_to_immature,
            value_t infected_to_infectious,
            value_t eggs_to_larva
        );

        template<class T>
        __host__ __device__
        void operator() (T tp) const { // tp - tuple

            value_t max_larval_capacity = thrust::get<0>(tp);
            value_t infection_prob = thrust::get<1>(tp);

            value_t x_0 = thrust::get<0>(thrust::get<2>(tp));
            value_t x_1 = thrust::get<1>(thrust::get<2>(tp));
            value_t x_2 = thrust::get<2>(thrust::get<2>(tp));
            value_t x_3 = thrust::get<3>(thrust::get<2>(tp));
            value_t x_4 = thrust::get<4>(thrust::get<2>(tp));
            value_t x_5 = thrust::get<5>(thrust::get<2>(tp));
            value_t x_6 = thrust::get<6>(thrust::get<2>(tp));
            value_t x_7 = thrust::get<7>(thrust::get<2>(tp));
            value_t x_8 = thrust::get<8>(thrust::get<2>(tp));


            //Eggs
            thrust::get<0>(thrust::get<3>(tp)) =
                -x_0*(kEggs_to_larva+kAquatic_mortality)
                + kLay_rate*kOviposition*(x_4+kInfected_egg_lay*(x_6+x_8));

            //Larva
            thrust::get<1>(thrust::get<3>(tp)) = 
                x_0*kEggs_to_larva*(1.0-x_1/max_larval_capacity)
                - x_1*(kAquatic_mortality+kEggs_to_larva);

            //Immature
            thrust::get<2>(thrust::get<3>(tp)) = 
                x_1*kLarva_to_immature
                - x_2*(kImmature_adult+kAdult_life);

            //Adult uninfected biting
            thrust::get<3>(thrust::get<3>(tp)) = 
                x_2*kImmature_adult
                - x_3*(kAdult_life+kFeed_rate)
                + x_4*kLay_rate;

            //Adult uninfected laying eggs
            thrust::get<4>(thrust::get<3>(tp)) = 
                x_3*kFeed_rate*(1.0-infection_prob)*(1.0-kHuman_feeding_mort)
                - x_4*(kAdult_life+kLay_rate);

            //Adult infected biting
            thrust::get<5>(thrust::get<3>(tp)) =
                -x_5*(kAdult_life+kInfected_to_infectious+kFeed_rate)
                + x_6*kLay_rate;

            //Adult infected laying eggs
            thrust::get<6>(thrust::get<3>(tp)) = 
                -x_6*(kAdult_life+kLay_rate+kInfected_to_infectious)
                + x_5*kFeed_rate*(1.0-kHuman_feeding_mort)
                + x_3*kFeed_rate*infection_prob*(1.0-kHuman_feeding_mort);

            //Adult infectious biting
            thrust::get<7>(thrust::get<3>(tp)) =
                x_5*kInfected_to_infectious
                - x_7*(kAdult_life+kFeed_rate)
                + x_8*kLay_rate;

            //Adult infectious laying eggs
            thrust::get<8>(thrust::get<3>(tp)) =
                x_6*kInfected_to_infectious
                - x_8*(kAdult_life+kLay_rate)
                + x_7*kFeed_rate*(1.0-kHuman_feeding_mort);

        }

    };

private:
    const size_t kNum_mosqs;
    const state_t kMax_larval_capacity_list;

    state_t population; // 9*kNum_mosqs element vector

    state_t infection_prob_list;

    value_t current_time;
    value_t current_temperature;
    value_t current_larva_to_immature;
    value_t current_infected_to_infectious;
    value_t current_eggs_to_larva;

public:

    MosquitoEckhoff_multi(
        std::vector<float>& host_max_larval_capacity_list
    );

    inline state_t& get_population(){
        return this->population;
    }

    void reset_time();
    void pre_integration_update(std::vector<float> host_infection_prob_list);

    void operator()(state_t& x, state_t& dxdt, value_t t) const;


};

class MosquitoManagerODEThrust {

    MosquitoEckhoff_multi* mosq_sys;

    MosquitoEckhoff_multi::stepper_type ode_stepper;


public:
    MosquitoManagerODEThrust(
        std::vector<float> max_larval_capacity_list
    );

    ~MosquitoManagerODEThrust();

    void step(const std::vector<float>& infection_prob_list);





};

}

#endif