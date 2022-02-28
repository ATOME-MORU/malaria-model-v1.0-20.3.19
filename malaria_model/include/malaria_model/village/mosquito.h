#ifndef MOSQUITO_H
#define MOSQUITO_H

#ifdef __NVCC__
#include "village/mosquito_ode.cuh"
#endif

#include <boost/numeric/odeint.hpp>
#include <boost/random/mersenne_twister.hpp> // boost::mt19937
#include <boost/random/normal_distribution.hpp> // boost::random::uniform_real_distribution

namespace odeint = boost::numeric::odeint;

namespace village {

class MosquitoODE {
public:
    virtual ~MosquitoODE(){}
    virtual double step(const double infection_prob) {return 0.0*infection_prob;}
    virtual double step_for_one_year(double infection_prob) {return 0.0*infection_prob;}
    virtual inline double get_pop(const int pop_index) const {return 0.0*pop_index;}
    virtual inline void set_pop(int pop_index, double pop) {(void)pop_index; (void)pop;}
    virtual inline double get_time() const{return 0.0;}
    virtual inline void set_time(const double new_time){(void)new_time;}
    virtual inline double get_capacity() const{return 0.0;}
    virtual inline bool if_capacity_reached(){return false;}
};

class MosquitoEckhoff: public MosquitoODE {
public:

    typedef std::vector< float > ode_state_t;
    typedef double ode_time_t;

    class Fx {

        // const double kPi = 3.14; // defined in common::

        // const double kAquatic_mortality = 0.1;
        const double kAquatic_mortality;
        const double kAquatic_arrhenius_1 = 8.42E+10;
        const double kAquatic_arrhenius_2 = 8328.0;

        // const double kImmature_adult = 1/2.0;
        // const double kImmature_adult = 1/2.0; //G_mature
        const double kImmature_adult; //G_mature, immature to adult rate, maturation

        const double kIadult_life;  // death rate of immature adults
        const double kAdult_life;   // death rate of adults
        // const double kAdult_life = 1/24.0;
        // const double kAdult_life = 1/10.0;
        // const double kAdult_life = 1/8.0;

        const double kHuman_feeding_mort = 0.00;
        // const double kInfect_human_feed = 1.5; 
        const double kInfected_arrhenius_1 = 1.17E+11;
        const double kInfected_arrhenius_2 = 8340.0;
        const double kEgg_arrhenius1 = 6.16E+07;
        const double kEgg_arrhenius2 = 5754.03;
        
        // const double kInfected_egg_lay = 0.8;
        const double kInfected_egg_lay = 1.0;


        // const double kLay_rate = 1/1.0; //lay/bite switch daily
        const double kLay_rate; //lay/bite switch daily

        // const double kOviposition=25.0;
        // const double kOviposition=10.0;
        // const double kOviposition=20.0;
        const double kOviposition; // eggs per lay
        // const double kOviposition = 8.0;

        const double kInfected_to_infectious_rate; // incubation

        // const double kMax_Larval_Capacity = 1E+3; // per village
        const double kFeed_rate; // Biting_rate
        const double kMax_larval_capacity; //multiplier applied on this
        const double kInfection_prob;
        
    public:
        
        // static constexpr double kNum_days_per_year = 365.0;

        Fx(
            double aquatic_mortality,
            double maturation_days,
            double immat_life_days,
            double adult_life_days,
            double laying_days,
            double eggs_per_lay,
            double incubation_days,
            double feed_rate,
            double max_larval_capacity,
            double infection_prob
        );

        void operator()(
            MosquitoEckhoff::ode_state_t &y,
            MosquitoEckhoff::ode_state_t &dydt,
            MosquitoEckhoff::ode_time_t t
        ) const;
    };




class stochastic_euler {
public:

    // typedef boost::array< double , N > state_type;
    // typedef boost::array< double , N > deriv_type;
    // typedef double value_type;
    // typedef double time_type;
    typedef unsigned short order_type;

    typedef boost::numeric::odeint::stepper_tag stepper_category;

    static order_type order( void ) { return 1; }

    template< class System >
    void do_step(
        System system,
        MosquitoEckhoff::ode_state_t &x,
        MosquitoEckhoff::ode_time_t t,
        MosquitoEckhoff::ode_time_t dt
    ) const {

        MosquitoEckhoff::ode_state_t det , stoch;
        for( size_t ii=0; ii < x.size(); ++ii ) {
            det.push_back(0.0);
            stoch.push_back(0.0);
        }

        // std::cout << "t=" << t << "\n";
        // std::cout << "dt=" << dt << "\n";

        system.first( x, det, t );
        system.second( x, stoch );

        for( size_t ii=0 ; ii<x.size() ; ++ii ) {
            // std::cout << "x[ii]=" << x[ii] << ", ";

            x[ii] += dt * det[ii] + sqrt( dt ) * stoch[ii];
            // std::cout << x[ii] << "\n";

            if (x[ii] < 0) {
                x[ii] =0;
            }
        }
    }
};



    class Fx_stoch {

        boost::mt19937 &m_rng;
        boost::normal_distribution<> m_dist;

    public:

        Fx_stoch(
            boost::mt19937 &rng,
            double sigma
        ):  m_rng(rng),
            m_dist(0.0, sigma){

        }

        void operator()(
            const MosquitoEckhoff::ode_state_t &x,
            MosquitoEckhoff::ode_state_t &dxdt
        ) {
            for( size_t ii=0 ; ii<x.size() ; ++ii ) {
                dxdt[ii] = m_dist( m_rng );
            }
        }
    };


    MosquitoEckhoff(
        double aquatic_mortality,
        double maturation_days,
        double immat_life_days,
        double adult_life_days,
        double laying_days,
        double eggs_per_lay,
        double incubation_days,
        double feed_rate,
        double max_larval_capacity
    );
    ~MosquitoEckhoff();

    void reset_time();
    double step(const double infection_prob);

    double step_for_one_year(double infection_prob);

    inline double get_pop(const int pop_index) const{
        return this->population.at(pop_index);
    }

    inline void set_pop(int pop_index, double pop) {
        this->population.at(pop_index) = pop;
    }

    inline double get_time() const{
        return this->current_time;
    }

    inline void set_time(const double new_time) {
        this->current_time = new_time;
    }

    inline double get_capacity() const{
        return this->kMax_larval_capacity;
    }

    inline bool if_capacity_reached() {
        // return (this->kMax_larval_capacity == this->get_pop(1));
        return (this->get_pop(1) > this->kMax_larval_capacity*0.9);
    }

private:

// typedef odeint::controlled_runge_kutta< runge_kutta_dopri5< ode_state_t > > dopri_stepper_type;
// typedef odeint::dense_output_runge_kutta< dopri_stepper_type > dense_stepper_type;

    odeint::runge_kutta4<ode_state_t> ode_stepper;
    const ode_time_t kDt = 0.05;
    // const ode_time_t kDt = 1;

    const double kAquatic_mortality;
    const double kMaturation_days;
    const double kImmat_life_days;
    const double kAdult_life_days;

    const double kLaying_days;
    const double kEggs_per_lay;
    const double kIncubation_days;
    const double kFeed_rate;
    const double kMax_larval_capacity;
    const ode_time_t kTime_length_per_integration = 1.0;

    ode_state_t population;
    ode_time_t current_time;

};



class MosquitoManager {

    bool verbose = false;

    const double kEquilibrium_epsilon = 1;
    const int kEquilibrium_max_num_years = 10;

    std::vector<MosquitoODE*> mosquito_of;
    std::vector<double> num_infectious_mosq_of;

public:

    int sum_num_larva_capacity = 0;

    //////////////////////////////////////////////////////////////////////
    //

    MosquitoManager(
        const int* population_of,
        float aquatic_mortality,
        float maturation_days,
        float immat_life_days,
        float adult_life_days,
        float laying_days,
        float eggs_per_lay,
        float incubation_days,
        double feed_rate,
        double carrying_capacity_multiplier,
        int num_mosquitos,
        std::string max_larval_capacity_file_name
    );

    MosquitoManager(
        const int* population_of,
        float aquatic_mortality,
        float maturation_days,
        float immat_life_days,
        float adult_life_days,
        float laying_days,
        float eggs_per_lay,
        float incubation_days,
        double feed_rate,
        double carrying_capacity_multiplier,
        std::vector<double> max_larval_capacity_list
    );

    std::vector<double> init_construct_max_larval_capacity_list_from_file(
        int num_mosquitos,
        std::string max_larval_capacity_file_name
    );

    ~MosquitoManager();

    //////////////////////////////////////////////////////////////////////
    // 

    void step(const std::vector<float>& infection_prob_list);

    bool init_step_until_equilibrium(const std::vector<float>& infection_prob_list);


    //////////////////////////////////////////////////////////////////////
    // Access
    inline double get_num_infectious_mosq(int mosq_index) const {
        return this->num_infectious_mosq_of[mosq_index];
    }

    inline const std::vector<double>& get_num_infectious_mosq_of() const {
        return this->num_infectious_mosq_of;
    }

    inline double get_mosq_pop(int mosq_index, int pop_index) const {
        return this->mosquito_of[mosq_index]->get_pop(pop_index);
    }

    inline void set_mosq_pop(int mosq_index, int pop_index, double pop) {
        this->mosquito_of[mosq_index]->set_pop(pop_index, pop);
    }

    inline double get_mosq_time(int mosq_index) const {
        return this->mosquito_of[mosq_index]->get_time();
    }

    inline void set_mosq_time(int mosq_index, double new_time) {
        this->mosquito_of[mosq_index]->set_time(new_time);
    }

    inline double get_mosq_capacity(int mosq_index) const {
        return this->mosquito_of[mosq_index]->get_capacity();
    }

    inline void set_verbose_on() {
        this->verbose = true;
    }

    inline void set_verbose_off() {
        this->verbose = false;
    }
};

}


#endif