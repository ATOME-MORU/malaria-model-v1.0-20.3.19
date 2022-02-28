#include "common/common.h"
#include <cmath>//pow()
#include <limits>

// #include <boost/math/special_functions/binomial.hpp>

#include <iostream>
namespace common{

std::ostream& operator<<(std::ostream& os, ParasiteType pt) {
    switch(pt) {
        case ParasiteType::kNoParasite  : os << "NP"; break;
        case ParasiteType::kRa          : os << "RA"; break;
        case ParasiteType::kRb          : os << "RB"; break;
        case ParasiteType::kR0          : os << "R0"; break;
        case ParasiteType::kRab         : os << "RAB"; break;
        default : os << "UnknownParasiteType";
        // default    : os.setstate(std::ios_base::failbit);
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, StageName pt) {
    switch(pt) {
        case StageName::kNotInSystem    : os << "NIS"; break;
        case StageName::kLiver          : os << "LVR"; break;
        case StageName::kBlood          : os << "BLD"; break;
        case StageName::kInfectious     : os << "INF"; break;
        default : os << "UnknownStageName";
        // default    : os.setstate(std::ios_base::failbit);
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, DrugName pt) {
    switch(pt) {
        case DrugName::kNoDrug      : os << "NDG"; break;
        case DrugName::kArtesunate  : os << "Art"; break;
        case DrugName::kPiperaquine : os << "Pip"; break;
        case DrugName::kPiparte     : os << "P+A"; break;
        case DrugName::kAQ          : os << "AQ"; break;
        case DrugName::kASAQ        : os << "ASAQ"; break;
        case DrugName::kLumefantrine: os << "LUM"; break;
        case DrugName::kAL          : os << "A+L"; break;
        case DrugName::kSP          : os << "S+P"; break;
        case DrugName::kChloroquine : os << "CLQ"; break;
        default : os << "UnknownDrugName";
        // default    : os.setstate(std::ios_base::failbit);
    }
    return os;
}

// float get_clinical_probability_single_infection(const uint8_t immunity_level,
//                                                 const uint8_t cummulative_exposures){
// // deme->probclinicallb=(0.1*exp(-(cumm-2)*0.1)+exp((-0.9*cumm)))/sqrt(level);
//     return  (
//             0.1*exp(-(cummulative_exposures-2)*0.1)
//             +exp(-0.9*cummulative_exposures)
//             )
//             / sqrt(immunity_level);
// }
// float get_clinical_probability_multiple_infection(const uint8_t immunity_level,
//                                                   const uint8_t cummulative_exposures,
//                                                   const uint8_t moi){
// // deme->probclinicalbb=exp(-0.15*(moi-1))*deme->probclinicallb;
//     return  exp(-0.15*(moi-1))*common::get_clinical_probability_single_infection(immunity_level, cummulative_exposures);
// }




float fractional_killing_function(float c, float gamma, float EC_50){
    const float kE_max = 100.0;
    return kE_max*pow(c,gamma)/(pow(c,gamma)+pow(EC_50,gamma));
}


void parasite_decay_and_replication(
                        uint64_t* prs_size,
                        const ParasiteType* prs_type,
                        const uint16_t* prs_age,
                        const uint64_t* imunity_size,
                        const float* immunity_rate,
                        const int index
                        ){
    const uint64_t kPrs_size_max = pow(10,18);

    if( prs_size[index] >= kPrs_size_max || prs_size[index] == 0){
        // given_size[index] = kPrs_size_max;
        // 1. number too large, leave it as it is
        // 2. no parasite, no action
    } else {
        prs_size[index]
            *= common::kParasite_replication_rate
                    [(int)prs_type[index]]
                    [prs_age[index] % common::kParasite_replication_cycle_length];
        
        // uint64_t decay_size = get_immunity_kill_size_constant(prs_type, prs_age,index);

        uint64_t decay_size = imunity_size[index] * pow( immunity_rate[index],prs_age[index]);

        if (prs_size[index] <= decay_size) {
            prs_size[index] = 0;
        } else {
            prs_size[index] -= decay_size;
        }
    }
}

uint64_t get_immunity_kill_size_constant(
                        const ParasiteType* prs_type,
                        const uint16_t* prs_age,
                        const int index
                        ){
    return common::kParasite_immunity_kill_constant
                [(int)prs_type[index]]
                [prs_age[index] % common::kParasite_decay_cycle_length];
}


// nothing below this line is used

void parasite_development_replication_flat(uint64_t* given_size,
                                            const ParasiteType* prs_type,
                                            const uint16_t* prs_age,
                                            const int index){
    const uint64_t kPrs_size_max = pow(10,18);

    if( given_size[index] >= kPrs_size_max){
        given_size[index] = kPrs_size_max;
    } else {
        given_size[index]
            *= common::kParasite_replication_rate
                    [(int)prs_type[index]]
                    [prs_age[index] % common::kParasite_replication_cycle_length];
    }
}
void parasite_development_decay_flat(uint64_t* given_size,
                                            const ParasiteType* prs_type,
                                            const uint16_t* prs_age,
                                            const int index){
    const uint64_t kPrs_size_max = pow(10,18);

    if( given_size[index] >= kPrs_size_max){
        given_size[index] = kPrs_size_max;
    } else {
        given_size[index]
            *= common::kParasite_decay_base_rate
                    [(int)prs_type[index]]
                    [prs_age[index] % common::kParasite_decay_cycle_length];
    }
}

float parasite_decay_phi(uint16_t days_since_last_exposure, ParasiteType prs_type){

    if(days_since_last_exposure == 0 ){
        //host has no immunity
        return common::kParasite_decay_Phi0;
    }

    const uint16_t kD_max = common::kParasite_immunity_days_max[(int)prs_type];

    if(kD_max == 0){
        //parasite type has no immunity period
        return common::kParasite_decay_Phi0;
    }

    if(days_since_last_exposure == 1){
        return kD_max+common::kParasite_decay_Phi0;
        // return kD_max;
    } else {
        // 140, 150, 160 // 40, 50, 60 with after_replication_size = prs_size[index] * ...
        return (parasite_decay_phi(days_since_last_exposure-1, prs_type)
                -common::kParasite_decay_Phi0 + kD_max) / days_since_last_exposure
                +common::kParasite_decay_Phi0;
        // return (parasite_decay_phi(days_since_last_exposure-1, prs_type)
        //         -common::kParasite_decay_Phi0 + kD_max) / pow(days_since_last_exposure,0.3);
                // +common::kParasite_decay_Phi0;

        // 40, 50, 60
        // return kD_max / pow(days_since_last_exposure,0.5);
                // +common::kParasite_decay_Phi0;
    }

}
void parasite_natural_development( const uint16_t* days_since_last_exposure,
                                    float* phi,
                                    // uint8_t imm_level,
                                    // uint8_t exp_count,
                                    ParasiteType* prs_type,
                                    const uint16_t* prs_age,
                                    uint64_t* prs_size,
                                    const int index,
                                    uint64_t* step_delta){

    // uint16_t kD_max = common::kParasite_immunity_days_max[(int)prs_type[index]];
    // float kPhi0 = common::kParasite_decay_Phi0;
    const uint64_t kPrs_size_max = pow(10,18);

    // if (days_since_last_exposure[index]>0) {
    //     phi[index] = (phi[index]-kPhi0+kD_max)/days_since_last_exposure[index];
    // } else {
    //     phi[index] = kPhi0;
    // }


    phi[index] = parasite_decay_phi(days_since_last_exposure[index], prs_type[index]);

    uint64_t after_decay_size = 0;
    uint64_t after_replication_size = 0;

    if( prs_size[index] >= kPrs_size_max){

        after_decay_size = kPrs_size_max;
        after_replication_size = kPrs_size_max;
        prs_size[index] = kPrs_size_max;

    } else {

        after_decay_size = prs_size[index]
                        * pow(common::kParasite_decay_base_rate
                                    [(int)prs_type[index]]
                                    [prs_age[index] % common::kParasite_decay_cycle_length]
                            , phi[index]);
        
        // after_decay_size = prs_size[index]
        //                 * (common::kParasite_decay_base_rate
        //                             [(int)prs_type[index]]
        //                             [prs_age[index] % common::kParasite_decay_cycle_length]
        //                     / phi[index]);

        // after_replication_size = prs_size[index]
        after_replication_size = after_decay_size
                        * common::kParasite_replication_rate
                                    [(int)prs_type[index]]
                                    [prs_age[index] % common::kParasite_replication_cycle_length];
        // after_replication_size = after_decay_size
        //                 * common::kParasite_replication_rate
        //                             [(int)prs_type[index]]
        //                             [prs_age[index] % common::kParasite_replication_cycle_length];
    }

    step_delta[0] = prs_size[index] - after_decay_size;
    step_delta[1] = after_replication_size - prs_size[index];
    // prs_size[index] = after_decay_size + after_replication_size - prs_size[index];
    prs_size[index] = after_replication_size;
}



void Bernstein_curve::calibrate(uint16_t T, uint16_t n, uint64_t y_max){
    this->t_scale = T;
    this->i = 1;
    this->n = n;
    this->t_i = (float)T/n;
    // this->nCi = (uint32_t) boost::math::binomial_coefficient<float>(this->n, this->i);
    this->nCi = this->n;
    this->y_scale = (float) y_max / (
            this->nCi
            * pow((float)this->t_i/this->t_scale, this->i)
            * pow(1-(float)this->t_i/this->t_scale, this->n-this->i)
        );

}

uint64_t Bernstein_curve::get_at(uint16_t t){
    if (t > this->t_scale) {
        return 0;
    }
    
   return this->y_scale
            * this->nCi
            * pow((float)t/this->t_scale,this->i)
            * pow(1-(float)t/this->t_scale, this->n-this->i);
}

uint64_t Bernstein_curve::get_at_level(uint16_t t, int level){
    if (t > this->t_scale) {
        return 0;
    }

    float level_y_scale = 1.1;
    float level_t_scale = 0.9;

    if (level == 0){
        if (t <= this->t_i) {
            return this->get_at(t) * level_y_scale;
        } else {
            return this->get_at(t_i) * level_y_scale + (this->get_at(t_i) -this->get_at(t_i-1)) * (t-this->t_i) ;
        }
    } else if (level == 1) {
        return this->get_at(t);
    } else {
        if (t >= this->t_scale*pow(level_t_scale,level-1)) {
            return 0;
        } else {
            return y_scale / pow(level_y_scale,level-1)
                * this->nCi
                * pow((float)t/(this->t_scale*pow(level_t_scale,level-1)),this->i)
                * pow(1-(float)t/(this->t_scale*pow(level_t_scale,level-1)), this->n-this->i);
        }
    }
    
    return y_scale
            * this->nCi
            * pow((float)t/this->t_scale,this->i)
            * pow(1-(float)t/this->t_scale, n-i);
}


}