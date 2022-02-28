#ifndef COMMON_H
#define COMMON_H

#include <cstdint> //uint8_t, etc.
#include <iostream> //std::ostream&

namespace common {

const float kPI = 3.1415927;
constexpr double kDays_per_year = 365.0;
const int kNum_days_in_one_year = 365;

enum class HumanMobilityType: char {
    kStatic,    // general population short term movement
    kTemporary, // temporary migration
    kMobile     // mobile population migration
};

// G_antroph = 0.95;
// G_indoorpref = 0.95;
// G_xitn =(0.0)/365.0;   //change
// G_itneffect = 0.5*G_antroph*G_indoorpref;

// Simulation

// Interventions
// const float kItn_effect = 0.5*0.95*0.95;

// enum ParasiteType { //temporary change for printing, need to fix this
enum class ParasiteType: char { 
    kNoParasite,
    kRa,
    kRb,
    kR0,
    kRab,
    First = kNoParasite,
    Last = kRab
};

std::ostream& operator<<(std::ostream& os, ParasiteType pt);

const int kParasite_replication_cycle_length = 2;
// e.g. a parasite which replicates 20 times every 48 hrs has rate profile {0,20}
const float kParasite_replication_rate
    [static_cast<int>(ParasiteType::Last)+1]
    [kParasite_replication_cycle_length] = {
        { 1.0,  1.0},
        { 1.0, 20.0},
        { 1.0, 20.0},
        { 1.0, 20.0},
        { 1.0, 20.0}
};
const int kParasite_decay_cycle_length = 2;
const uint64_t kParasite_immunity_kill_constant
    [static_cast<int>(ParasiteType::Last)+1]
    [kParasite_decay_cycle_length] = {
        { 0,  0},
        { 18000,  18000},
        { 18096,  18096},
        { 6,  6},
        { 18096,  18096}
};
const uint64_t immunity_kill_size_day_zero
    [static_cast<int>(ParasiteType::Last)+1] = {
         8000,
         8000,
         8000,
         8000,
         8000
};

const float kParasite_decay_base_rate
    [static_cast<int>(ParasiteType::Last)+1]
    [kParasite_decay_cycle_length] = {
        { 1.0,  1.0},
        { 0.4,  0.4},
        { 0.2,  0.2},
        { 0.9,  0.9},
        { 0.2,  0.2}
};
const uint16_t kParasite_immunity_days_max
// parasite immunity-based decay multiplier (# of days)
    [static_cast<int>(ParasiteType::Last)+1] ={
        1,
        // 40,
        // 50,
        // 60
        140,
        150,
        160,
        150
};
// const uint16_t kParasite_d_max
// // parasite immunity-based decay multiplier (# of days)
//     [static_cast<int>(ParasiteType::Last)+1] ={
//         1,
//         40,
//         50,
//         60
// };
const float kParasite_decay_Phi0 = 1;

// const int init_parasite_distribution [static_cast<int>(ParasiteType::Last)+1] = {
//     10,
    
// };

enum class DrugName: char {
    /*0*/    kNoDrug,
    /*1*/    kArtesunate,
    /*2*/    kPiperaquine,
    /*3*/    kPiparte,
    /*4*/    kAQ, // mmc_wp2_rev17
    /*5*/    kASAQ, // mmc_wp2_rev17
    /*6*/    kLumefantrine,
    /*7*/    kAL,
    /*8*/    kSP,
    /*9*/    kChloroquine,

    First = kNoDrug,
    Last = kChloroquine
};
std::ostream& operator<<(std::ostream& os, DrugName drg);


const float kDrug_half_life_day
    [static_cast<int>(DrugName::Last)+1] = {
        0.0, 1.5, 20.0, 10.0 
};

// const float kDrug_

enum class StageName: char {
    kNotInSystem,
    kLiver,
    kBlood,
    kInfectious,

    First = kNotInSystem,
    Last = kInfectious
};
std::ostream& operator<<(std::ostream& os, StageName stg);

// number of parasites killed per day, per parasite/drug/stage, at full-life
const int kPds_kill_size[static_cast<int>(ParasiteType::Last)+1]
                     [static_cast<int>(DrugName::Last)+1]
                     [static_cast<int>(StageName::Last)+1]= {
    // no parasite
    {
        //no drug
        {0,0,0,0},
        // argesunate
        {0,0,0,0},
        // piperaquine
        {0,0,0,0},
        // piparte
        {0,0,0,0}
    },
    // parasite ra
    {
        //no drug
        {0,0,0,0},
        // argesunate
        {0,0,600,0},
        // piperaquine
        {0,0,0,0},
        // piparte
        {0,0,0,0}
    },
    // parasite rb
    {
        //no drug
        {0,0,0,0},
        // argesunate
        {0,0,0,0},
        // piperaquine
        {0,0,0,0},
        // piparte
        {0,0,0,0}
    },
    // parasite r0
    {
        //no drug
        {0,0,0,0},
        // argesunate      
        {0,0,0,0},  
        // piperaquine
        {0,0,0,0},
        // piparte
        {0,0,0,0}
    }
};


// probability of parasites killed per day, per parasite/drug/stage
// TODO: should make this a map? map not available in CUDA?
// const float kPds_kill_probability[static_cast<int>(ParasiteType::Last)+1]
//                      [static_cast<int>(DrugName::Last)+1]
//                      [static_cast<int>(StageName::Last)+1]= {
//     // no parasite
//     {
//         //no drug
//         {0,0,0,0},
//         // Artesunate
//         {0,0,0,0},
//         // piperaquine
//         {0,0,0,0},
//         // piparte
//         {0,0,0,0},
//         // lumefantrine
//         {0,0,0,0},
//         // al
//         {0,0,0,0},
//         // SP
//         {0,0,0,0},
//         // chloroquine
//         {0,0,0,0}
//     },
//     // parasite ra
//     {
//         //no drug
//         {0,0,0,0},
//         // Artesunate
//         {0,0,1.0/5.0*0.27,              1.0/3.0*0.27},
//         // piperaquine
//         {0,0,1.0/3.0,                   1.0/15.0},
//         // piparte
//         {0,0,1.0/5.0*0.27+1.0/3.0*0.73, 1.0/3.0*0.27+1.0/15.0*0.73},
//         // lumefantrine
//         {0,0,0,0},
//         // al
//         {0,0,0,0},
//         // SP
//         {0,0,1.0/5.0,                   1.0/3.0},
//         // chloroquine
//         {0,0,1.0/5.0,                   1.0/3.0}
//     },
//     // parasite rb
//     {
//         //no drug
//         {0,0,0,0},
//         // Artesunate
//         {0,0,1.0/5.0,                   1.0/3.0},
//         // {0,0,0.0,                   0.0},
//         // piperaquine
//         {0,0,1.0/3.0*0.8,               1.0/15.0*0.8},
//         // {0,0,1.0/3.0*0.0,               1.0/15.0*0.0},
//         // piparte
//         // {0,0,1.0/5.0*0.8+1.0/5.0*0.2,   1.0/3.0*0.8+1.0/3.0*0.2},
//         {0,0,1.0/3.0*0.8+1.0/5.0*0.2,   1.0/15.0*0.8+1.0/3.0*0.2},
//         // {0,0,0.0,   0.0},
//         // lumefantrine
//         {0,0,0,0},
//         // al
//         {0,0,0,0},
//         // SP
//         {0,0,1.0/5.0,                   1.0/3.0},
//         // chloroquine
//         {0,0,1.0/5.0,                   1.0/3.0}
//     },
//     // parasite r0
//     {
//         //no drug
//         {0,0,0,0},
//         // Artesunate      
//         {0,0,1.0/5.0,       1.0/3.0},  
//         // piperaquine
//         {0,0,1.0/3.0,       1.0/15.0},
//         // piparte
//         {0,0,1.0/5.0,       1.0/3.0},
//         // lumefantrine
//         {0,0,0,0},
//         // al
//         {0,0,0,0},
//         // SP
//         {0,0,1.0/5.0,       1.0/3.0},
//         // chloroquine
//         {0,0,1.0/5.0,       1.0/3.0}
//     }
// };

// // Artesunate (a) + Blood

// kPds_kill_probability[static_cast<int>ParasiteType:kR0][static_cast<int>DrugName:kArtesunate][static_cast<int>StageName:kBlood] = 1.0/5.0;
// kPds_kill_probability[static_cast<int>ParasiteType:kRa][static_cast<int>DrugName:kArtesunate][static_cast<int>StageName:kBlood] = 1.0/5.0 * 0.27;
// kPds_kill_probability[static_cast<int>ParasiteType:kRb][static_cast<int>DrugName:kArtesunate][static_cast<int>StageName:kBlood] = 1.0/5.0;

// // Artesunate (a) + Infectious

// kPds_kill_probability[static_cast<int>ParasiteType:kR0][static_cast<int>DrugName:kArtesunate][static_cast<int>StageName:kInfectious] = 1.0/3.0;
// kPds_kill_probability[static_cast<int>ParasiteType:kRa][static_cast<int>DrugName:kArtesunate][static_cast<int>StageName:kInfectious] = 1.0/3.0 * 0.27;
// kPds_kill_probability[static_cast<int>ParasiteType:kRb][static_cast<int>DrugName:kArtesunate][static_cast<int>StageName:kInfectious] = 1.0/3.0;

// // Piperaquine (b) + Blood

// kPds_kill_probability[static_cast<int>ParasiteType:kR0][static_cast<int>DrugName:kPiperaquine][static_cast<int>StageName:kBlood] = 1.0/3.0;
// kPds_kill_probability[static_cast<int>ParasiteType:kRa][static_cast<int>DrugName:kPiperaquine][static_cast<int>StageName:kBlood] = 1.0/3.0;
// kPds_kill_probability[static_cast<int>ParasiteType:kRb][static_cast<int>DrugName:kPiperaquine][static_cast<int>StageName:kBlood] = 1.0/3.0*0.8;

// // Piperaquine (b) + Infectious

// kPds_kill_probability[static_cast<int>ParasiteType:kR0][static_cast<int>DrugName:kPiperaquine][static_cast<int>StageName:kInfectious] = 1.0/15.0;
// kPds_kill_probability[static_cast<int>ParasiteType:kRa][static_cast<int>DrugName:kPiperaquine][static_cast<int>StageName:kInfectious] = 1.0/15.0;
// kPds_kill_probability[static_cast<int>ParasiteType:kRb][static_cast<int>DrugName:kPiperaquine][static_cast<int>StageName:kInfectious] = 1.0/15.0*0.8;

// // Piparte (ab) + Blood

// kPds_kill_probability[static_cast<int>ParasiteType:kR0][static_cast<int>DrugName:kPiparte][static_cast<int>StageName:kBlood] = 1.0/5.0;
// kPds_kill_probability[static_cast<int>ParasiteType:kRa][static_cast<int>DrugName:kPiparte][static_cast<int>StageName:kBlood] = 
//     kPds_kill_probability[static_cast<int>ParasiteType:kR0][static_cast<int>DrugName:kPiparte][static_cast<int>StageName:kBlood] * 0.27
//     + kPds_kill_probability[static_cast<int>ParasiteType:kR0][static_cast<int>DrugName:kPiperaquine][static_cast<int>StageName:kBlood] * (1-0.27);
// kPds_kill_probability[static_cast<int>ParasiteType:kRb][static_cast<int>DrugName:kPiparte][static_cast<int>StageName:kBlood] = 
//     kPds_kill_probability[static_cast<int>ParasiteType:kR0][static_cast<int>DrugName:kPiparte][static_cast<int>StageName:kBlood] * 0.8
//     + kPds_kill_probability[static_cast<int>ParasiteType:kR0][static_cast<int>DrugName:kArtesunate][static_cast<int>StageName:kBlood] * (1-0.8);

// // Piparte (ab) + Infectious

// kPds_kill_probability[static_cast<int>ParasiteType:kR0][static_cast<int>DrugName:kPiparte][static_cast<int>StageName:kInfectious] = 1.0/3.0;
// kPds_kill_probability[static_cast<int>ParasiteType:kRa][static_cast<int>DrugName:kPiparte][static_cast<int>StageName:kInfectious] = 
//     kPds_kill_probability[static_cast<int>ParasiteType:kR0][static_cast<int>DrugName:kPiparte][static_cast<int>StageName:kInfectious] * 0.27
//     + kPds_kill_probability[static_cast<int>ParasiteType:kR0][static_cast<int>DrugName:kPiperaquine][static_cast<int>StageName:kInfectious] * (1-0.27);
// kPds_kill_probability[static_cast<int>ParasiteType:kRb][static_cast<int>DrugName:kPiparte][static_cast<int>StageName:kInfectious] = 
//     kPds_kill_probability[static_cast<int>ParasiteType:kR0][static_cast<int>DrugName:kPiparte][static_cast<int>StageName:kInfectious] * 0.8
//     + kPds_kill_probability[static_cast<int>ParasiteType:kR0][static_cast<int>DrugName:kArtesunate][static_cast<int>StageName:kInfectious] * (1-0.8);


// struct ParasiteGroup {
//     ParasiteType type;
//     StageName stage;
//     // DrugName drug;
//     uint32_t size_parasite;
//     uint32_t size_gametocyte;
//     // uint32_t system_id;
//     uint16_t num_days;
// };

// float get_clinical_probability_single_infection(const uint8_t immunity_level,
//                                                 const uint8_t cummulative_exposures);
// float get_clinical_probability_multiple_infection(const uint8_t immunity_level,
//                                                   const uint8_t cummulative_exposures,
//                                                   const uint8_t moi);

// equation S3 of www.pnas.org/cgi/content/short/1006113108
// c - given drug concentration
float fractional_killing_function(float c, float gamma, float EC_50);

// void parasite_decay_and_replication(
//                         uint64_t* prs_size,
//                         const ParasiteType* prs_type,
//                         const uint16_t* prs_age,
//                         const int index
//                         );
void parasite_decay_and_replication(
                        uint64_t* prs_size,
                        const ParasiteType* prs_type,
                        const uint16_t* prs_age,
                        const uint64_t* imunity_size,
                        const float* immunity_rate,
                        const int index
                        );
uint64_t get_immunity_kill_size_constant(
                        const ParasiteType* prs_type,
                        const uint16_t* prs_age,
                        const int index
                        );



// nothing below this line is used

void parasite_development_replication_flat(uint64_t* given_size,
                                            const ParasiteType* prs_type,
                                            const uint16_t* prs_age,
                                            const int index);

void parasite_development_decay_flat(uint64_t* given_size,
                                            const ParasiteType* prs_type,
                                            const uint16_t* prs_age,
                                            const int index);

float parasite_decay_phi(uint16_t days_since_last_exposure, ParasiteType prs_type);
void parasite_natural_development( const uint16_t* days_since_last_exposure,
                                    float* phi,
                                    // uint8_t imm_level,
                                    // uint8_t exp_count,
                                    ParasiteType* prs_type,
                                    const uint16_t* prs_age,
                                    uint64_t* prs_size,
                                    const int index,
                                    uint64_t* step_delta);


struct Bernstein_curve{
    uint16_t i = 1;
    uint16_t n = 1;
    uint32_t nCi = 1;

    float y_scale = 1;
    uint16_t t_scale = 1;
    float t_i = 1.0;

    void calibrate(uint16_t t_width, uint16_t yt_max, uint64_t y_max);

    uint64_t get_at(uint16_t t);
    uint64_t get_at_level(uint16_t t, int level);
};

}

#endif