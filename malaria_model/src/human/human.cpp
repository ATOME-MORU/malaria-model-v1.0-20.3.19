#include <iostream>

#include <random>
#include <array>

#include <cassert>

#include <map>
#include <string>
#include <iomanip>

#include <fstream>

#include <algorithm>    // std::fill
#include <vector>
#include <list>

// #include "common/common.h"
#include "human/human.h"
#include "village/village.h"
#include "util/util.h"
#include "util/randomness.h"

namespace human {

HumanManager::HumanManager (
        int initial_population
    ) :
        HumanManager(
            initial_population,
            0,
            false, 0.0, 0.0, 1.0
        )
     {
}

HumanManager::HumanManager (
        int initial_population,
        int func_death_probability_version,

        bool attractiveness_enabled,
        float attractiveness_mean,
        float attractiveness_sd,
        float reset_log_bitten_times_at_the_end_of_year
    ) : 
        kFunc_death_probability_version(func_death_probability_version),

        kAttractiveness_enabled(attractiveness_enabled),
        kAttractiveness_mean(attractiveness_mean),
        kAttractiveness_sd(attractiveness_sd),
        attractiveness_of(initial_population,1.0),
        kReset_log_bitten_times_on_day(static_cast<int>(reset_log_bitten_times_at_the_end_of_year*365)),
        log_bitten_times_m2h(initial_population,0),
        log_bitten_times_h2m(initial_population,0)
    {

    if (initial_population<1) {

        std::cout << "Can not initialise human population to "
                  << initial_population << "." << std::endl;

    } else {

        this->age_of    = new int[initial_population];
        this->gender_of = new char[initial_population];
        this->death_probability_of = new float[initial_population];
        // village_home_of = new village::Village[initial_population];
        // village_at_of   = new village::Village[initial_population];
        this->home_village_of = new int[initial_population];
        this->at_village_of = new int[initial_population];

        this->if_itn_of = new bool[initial_population];
        this->mobility_of = new common::HumanMobilityType[initial_population];

        this->human_reg = new Human[initial_population];

        this->bitten_by_of = new common::ParasiteType[initial_population];
        this->given_drug_of = new common::DrugName[initial_population];

        this->sum_num_humans = 0;

        for (int hh = 0; hh < initial_population; hh++) {

            this->age_of[hh]    = 40;
            this->gender_of[hh] = 'f';
            this->death_probability_of[hh] = 0.0;
            this->home_village_of[hh] = 0;
            this->at_village_of[hh]   = 0;

            this->if_itn_of[hh] = false;
            this->mobility_of[hh] = common::HumanMobilityType::kStatic;

            // this->human_reg[hh] = Human(
            //                         this->age_of + hh,
            //                         this->gender_of + hh,
            //                         this->home_village_of + hh,
            //                         this->at_village_of + hh,
            //                         this->given_drug_of + hh
            //                         );

            // this->human_reg.push_back( Human(
            //                         this->age_of + hh,
            //                         this->gender_of + hh,
            //                         this->home_village_of + hh,
            //                         this->at_village_of + hh,
            //                         this->given_drug_of + hh
            //                         ));

            this->human_reg[hh].age = this->age_of + hh;
            this->human_reg[hh].gender = this->gender_of + hh;
            this->human_reg[hh].home_village = this->home_village_of + hh;
            this->human_reg[hh].at_village = this->at_village_of + hh;
            this->human_reg[hh].given_drug = this->given_drug_of + hh;

            this->bitten_by_of[hh] = common::ParasiteType::kNoParasite;
            this->given_drug_of[hh] = common::DrugName::kNoDrug;

            this->sum_num_humans++;

        }
        
    }

    std::cout << "Population: " << this->sum_num_humans << std::endl;

    this->init_attractiveness();

}

HumanManager::~HumanManager() {
    if (age_of) {delete[] age_of;}
    if (gender_of) {delete[] gender_of;}
    if (death_probability_of) {delete[] death_probability_of;}
    if (home_village_of) {delete[] home_village_of;}
    if (at_village_of) {delete[] at_village_of;}

    if (if_itn_of) {delete[] if_itn_of;}
    if (mobility_of) {delete[] mobility_of;}

    if (human_reg) {delete[] human_reg;}
    if (bitten_by_of) {delete[] bitten_by_of;}
    if (given_drug_of) {delete[] given_drug_of;}
}

void HumanManager::print_all() const{
    // util::Printables<human::HumanManager> pt;
    // pt.print_data_table(this, this->sum_num_humans, 5, 2);
    util::print_data_table(this, this->sum_num_humans, 7, 2);
}

void HumanManager::print_one(const int index_human) const{
    std::cout   << std::setw(8) <<  index_human
                << " home:" << std::setw(4) << this->home_village_of[index_human]
                << ", at:" << std::setw(4) << this->at_village_of[index_human]
                << ", age:" << std::setw(2) << this->age_of[index_human]
                << ", gender:" << this->gender_of[index_human];
}

void HumanManager::read_age_distribution(std::string input_file_name) {
    
    std::ifstream input_files_stream(input_file_name);

    if (!input_files_stream) {

        std::cout << "Could not open input file (" << input_file_name << ")." << std::endl;
    
    } else {

        unsigned int rr = 0;
        std::string str_line;

        // while(!std::getline(input_files_stream, str_line).eof()){
        while (input_files_stream >> str_line) {

            // if (rr==this->age_weight_vector.size()){
            if (rr >= kNum_age_bins){
                std::cout << "\033[31;47mHumanManager: Age distribution file has more rows than the requested "
                          << rr << "."
                          << "\033[0m" << std::endl;
                break;
            }

            this->age_weight_vector[rr++] = std::stoi(str_line);

        }

    }

    input_files_stream.close();

// vector<int> vec (arr, arr + sizeof(arr) / sizeof(arr[0]) );

}

void HumanManager::read_death_prob_distribution(std::string input_file_name) {
    
    std::ifstream input_files_stream(input_file_name);

    if (!input_files_stream) {

        std::cout << "Could not open input file (" << input_file_name << ")." << std::endl;
        std::cin.ignore();
    
    } else {

        int rr = 0;
        std::string str_line;

        // while(!std::getline(input_files_stream, str_line).eof()){
        while (input_files_stream >> str_line) {
            // if (rr==this->age_weight_vector.size()){
            if (rr >= kNum_age_bins){
                std::cout << "\033[31;47mHumanManager: Death probability distribution file has more rows than the requested "
                          << rr << "."
                          << "\033[0m" << std::endl;
                break;
            }
            this->death_prob_vector[rr++] = std::stof(str_line)*0.4/365.0;
        }
    }
    input_files_stream.close();

    // for (auto pp:this->death_prob_vector) {
    //     std::cout << pp << std::endl;
    // }
}

// void HumanManager::generate_age_distribution() {
//     for (int& aw : this->age_weight_vector) {
//         aw = 0;

//     }
//     for (int aa = 0; aa < this->age_weight_vector.size(); aa++) {
//         this->age_weight_vector.at(aa) = 0;
//         this->age_weight_vector.at(aa) = func_age_weight_polynmial(aa);
//     }

// }

// int HumanManager::func_age_weight_polynmial(int age) {
//     float G_age1 = 6.802e-03;
//     float G_age2 = 6.755e-03;
//     float G_age3 = -7.683e-04;
//     float G_age4 = 3.663e-05;
//     float G_age5 = -8.972e-07;
//     float G_age6 = 1.181e-08;
//     float G_age7 = -7.965e-11;
//     float G_age8 = 2.164e-13;
//     //agep[k]= G_age1 + G_age2*k + G_age3*pow(k,2) + G_age4*pow(k,3) + G_age5*pow(k,4) + G_age6*pow(k,5) + G_age7*pow(k,6) + G_age8*pow(k,7);
//     return static_cast<int>( 
//         G_age1
//         + G_age2*age
//         + G_age3*pow(age,2)
//         + G_age4*pow(age,3)
//         + G_age5*pow(age,4)
//         + G_age6*pow(age,5)
//         + G_age7*pow(age,6)
//         + G_age8*pow(age,7)
//     );
// }

void HumanManager::init_age() {
    this->init_age_discrete_distribution();
    this->update_all_death_probability();
}
void HumanManager::init_age_discrete_distribution() {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::discrete_distribution<std::size_t> dist(this->age_weight_vector.begin(), this->age_weight_vector.end());

    for (int hh = 0; hh < this->sum_num_humans; hh++) {
            this->age_of[hh] = dist(gen);
    }
}
void HumanManager::init_age_mixed_normal() {
    std::random_device rd;
    std::mt19937 gen(rd());

    auto normal_dist = std::array<std::normal_distribution<>,3> {
        std::normal_distribution<>(10,5),
        std::normal_distribution<>(30,5),
        std::normal_distribution<>(50,5)
    };

    std::vector<int> weight_list{3, 4, 3};

    // auto weights = std::discrete_distribution<std::size_t>{
    //     0.3,
    //     0.4,
    //     0.3
    // };

    std::discrete_distribution<std::size_t> weights(weight_list.begin(), weight_list.end());

    for (int hh = 0; hh < this->sum_num_humans; hh++) {
        auto which_dist = weights(gen);
        int age = std::round(normal_dist[which_dist](gen));
        if (age > 0 && age < 100) {

            this->age_of[hh] = age;

        } else {

            hh--;

        }
    }

}
void HumanManager::print_age_distribution() const{
    print_distribution<int>(this->age_of, this->sum_num_humans);
}

void HumanManager::increase_all_ages_by_one() {
    for (int hh = 0; hh < this->sum_num_humans; hh++) {
        this->age_of[hh]++;
    }
    this->update_all_death_probability();
}

int HumanManager::birth_and_death(
        bool* reset_flags,
        const float* random_numbers_array,
        const int random_numbers_array_size
    ) {
    
    int num_resets = 0;

    assert(random_numbers_array_size == this->sum_num_humans);

    for (int hh = 0; hh < this->sum_num_humans; hh++){
        // if (this->birth_and_death_probability_based(hh, this->age_of, random_numbers)) {
        // if (this->birth_and_death_probability_based(hh, this->age_of, random_numbers_array)) {
        if (random_numbers_array[hh] < this->death_probability_of[hh]) {
            this->reset_human_to_birth(hh);
            reset_flags[hh] = true;
            num_resets++;
        }
    }
    return num_resets;
}
// bool HumanManager::birth_and_death_probability_based(
//     const int h_index,
//     const int* ages,
//     const float* random_numbers
//     ) {
// // deathprob = (0.5*(exp(-age*G_sig)))+(2.5*G_sig*exp((1/G_mu1*age)/(1/1.8)));
// // deme->death = deathprob*G_mu;
//     const float G_sig = 1.0/15.0;
//     const float G_mu1 = 80.0;
//     const float G_mu = 0.0222/365.0; // death rate 0.0222 - 45 yrs life expectancy

//     return (random_numbers[h_index]
//                 < G_mu * (
//                             (0.5*(exp(-ages[h_index]*G_sig)))
//                             + (2.5*G_sig*exp((1/G_mu1*ages[h_index])/(1/1.8)))
//                          )
//             );
// }
void HumanManager::update_all_death_probability() {
    for (int hh = 0; hh < this->sum_num_humans; hh++) {
        this->death_probability_of[hh] = this->func_death_probability(
            this->age_of[hh]
        );
        // this->death_probability_of[hh] =
        //     this->death_probability_age_in_years_kernel(
        //         this->age_of[hh],
        //         this->kDeath_func_sig,
        //         this->kDeath_func_mu1,
        //         this->kDeath_func_mu
        //     );
        // std::cout << "age=" << this->age_of[hh]
        //             << "dprob=" << this->death_probability_of[hh]
        //             << "\n";
        // std::cin.ignore();
    }
}
float HumanManager::func_death_probability (
        const int age_in_years
    ) {

    assert(age_in_years >= 0);

    switch(this->kFunc_death_probability_version) {
        case 0 : // karenjan17.cpp
            return this->func_death_probability_kernel_0(
                age_in_years
            );
            break;
        case 1 : // mmc_wp2
            return this->func_death_probability_kernel_1_1(
                age_in_years
            );
            break;
        case 2 : // mda_coverage, 03 Dec 2018
            return this->func_death_probability_kernel_2(
                age_in_years
            );
            break;
        case 21 : // mosq_apx, 03 Sep 2019
            return this->func_death_probability_kernel_2_1(
                age_in_years
            );
            break;
        case 3 : // death prob read from file
            return this->func_death_probability_kernel_3(
                age_in_years
            );
            break;
        default :
            assert(false && "Unknown HM.kFunc_death_probability_version");
            break;
    }

}
float HumanManager::func_death_probability_kernel_3(
        const int age_in_years
    ) {
    assert(age_in_years >= 0);
    if (!(age_in_years < static_cast<int>(this->death_prob_vector.size()))){
        std::cout << static_cast<int>(this->death_prob_vector.size());
        assert(age_in_years < static_cast<int>(this->death_prob_vector.size()));
    }
    return this->death_prob_vector[age_in_years];
}

float HumanManager::func_death_probability_kernel_2_1(
        const int age_in_years
    ) {
// mosq_apx, 03 Sep 2019
// added the `age>75` filter

    float death_prob = HumanManager::func_death_probability_kernel_2(age_in_years);

    if (age_in_years > 75) {
        return death_prob * 2;
    } else {
        return death_prob;
    }

}
float HumanManager::func_death_probability_kernel_2(
        const int age_in_years
    ) {
// mda_coverage, 03 Dec 2018
// deathprob = (0.65*(exp(-age*G_sig)))+(2.5*G_sig*exp((1/G_mu1*age)/(1/6.5)));
// G_sig = 1.0/15.0;
// G_mu1 = 80.0;
    const float kMu = 0.0222/365.0;
    const float kMu1 = 80.0;
    const float kSig = 1.0/15.0;

// The above setting gives Africa age distribution
// Revised on 03 Sep 2019
    // const float kMu = 0.0282/365.0;
    // const float kMu = 0.0152/365.0;
    // const float kMu1 = 180.0;

    // return (
    //     kMu * (
    //         (0.65*(exp(-age_in_years*kSig)))
    //         +(2.5*kSig*exp((1/kMu1*age_in_years)/(1/6.5)))
    //     )
    // );
    return (
        kMu * (
            (0.65*(exp(-age_in_years*kSig)))
            +(2.5*kSig*exp((1/kMu1*age_in_years)/(1/3.5)))
        )
    );
}
float HumanManager::func_death_probability_kernel_1_1(
        const int age_in_years
    ) {
// mmc_wp2
// same function as durgrestest.cpp but different parameters

    // // mmc_wp2, scenario 1
    const float G_mu = 0.0222/365.0; // 6.0822e-05
    const float G_mu1 = 80.0;
    const float G_sig = 1.0/15.0;    // 0.0666667

    return (
        ((1-(exp(-(age_in_years/50.0))))
        + ( 1/(G_sig*sqrt(2*3.14)) )
        * exp(-(pow(age_in_years-G_mu1,2))/(2*(pow(G_sig,2)))
        ))
        * G_mu / 0.8
    );

}
float HumanManager::func_death_probability_kernel_1_0(
        const int age_in_years
    ) {
// durgrestest.cpp / aguasSIR.h, mmc_wp2

// drugrestest.cpp which uses setprobdeath() in aguasDIR.h
// deathprob = (1-(exp(-k1)))+(1/(G_sig*sqrt(2*3.14)))*exp(-(pow(age-G_mu1,2))/(2*(pow(G_sig,2))));
// deme->death = deathprob*G_mu/0.8;
// G_mu = 0.0222/365.0; // 6.0822e-05
// G_mu1 = 35.0;
// G_sig = 1.0/10.0
// G_mu2 = 50.0;
// k1 = age/G_mu2;

    // durgrestest.cpp aguasSIR.h
    const float G_mu = 0.0222/365.0; // 6.0822e-05
    const float G_mu1 = 35.0;
    const float G_sig = 1.0/10.0;

    return (
        ((1-(exp(-(age_in_years/50.0))))
        + ( 1/(G_sig*sqrt(2*3.14)) )
        * exp(-(pow(age_in_years-G_mu1,2))/(2*(pow(G_sig,2)))
        ))
        * G_mu / 0.8
    );

}
float HumanManager::func_death_probability_kernel_0(
        const int age_in_years
    ) {
// karenjan17.cpp
// deathprob = (0.5*(exp(-age*G_sig)))+(2.5*G_sig*exp((1/G_mu1*age)/(1/1.8)));
// deme->death = deathprob*G_mu;

    const float G_mu = 0.0222/365.0; // death rate 0.0222 - 45 yrs life expectancy
    const float G_mu1 = 80.0;
    const float G_sig = 1.0/15.0;

    // (void) age_in_years;
    // return 1/(365*55);

    return (
        G_mu * (
            (0.15*(exp(-age_in_years*G_sig)))
            + (2.5*G_sig*exp((1/G_mu1*age_in_years)/(1/1.8)))
        )
    );

}
float HumanManager::death_probability_age_in_years_kernel(
    const int age_in_years,
    const float G_sig,
    const float G_mu1,
    const float G_mu
    ) {
// karenjan17.cpp
// deathprob = (0.5*(exp(-age*G_sig)))+(2.5*G_sig*exp((1/G_mu1*age)/(1/1.8)));
// deme->death = deathprob*G_mu;

    // const float G_sig = 1.0/15.0;
    // const float G_mu1 = 80.0;
    // const float G_mu = 0.0222/365.0; // death rate 0.0222 - 45 yrs life expectancy

    // (void) age_in_years;
    // return 1/(365*55);

    // return (
    //     G_mu * (
    //         (0.15*(exp(-age_in_years*G_sig)))
    //         + (2.5*G_sig*exp((1/G_mu1*age_in_years)/(1/1.8)))
    //     )
    // );

// durgrestest.cpp aguasSIR.h, mmc_wp2
// float k1;
// G_mu2 = 50.0;
// k1 = age/G_mu2;
// deathprob = (1-(exp(-k1)))+(1/(G_sig*sqrt(2*3.14)))*exp(-(pow(age-G_mu1,2))/(2*(pow(G_sig,2))));
// deme->death = deathprob*G_mu/0.8;

    // return (
    //     ((1-(exp(-(age_in_years/50.0))))
    //     + ( 1/(G_sig*sqrt(2*3.14)) )
    //     * exp(-(pow(age_in_years-G_mu1,2))/(2*(pow(G_sig,2)))
    //     ))
    //     * G_mu / 0.8
    // );

// mda_coverage
// deathprob = (0.65*(exp(-age*G_sig)))+(2.5*G_sig*exp((1/G_mu1*age)/(1/6.5)));
// G_sig = 1.0/15.0;
// G_mu1 = 80.0;
    (void) G_sig;
    (void) G_mu1;
    (void) G_mu;
    const float kSig = 1.0/15.0;
    const float kMu1 = 80.0;
    const float kMu = 0.0222/365.0;

    return (
        kMu * (
            (0.65*(exp(-age_in_years*kSig)))
            +(2.5*kSig*exp((1/kMu1*age_in_years)/(1/6.5)))
        )
    );

}
void HumanManager::reset_human_to_birth(const int h_index) {
    this->age_of[h_index] = 0;
    this->death_probability_of[h_index] = this->func_death_probability(
        this->age_of[h_index]
    );
    // this->gender_of[h_index]
    this->at_village_of[h_index] = this->home_village_of[h_index];

    this->bitten_by_of[h_index] = common::ParasiteType::kNoParasite;
    this->given_drug_of[h_index] = common::DrugName::kNoDrug;

    this->attractiveness_of[h_index] = this->func_gen_attractiveness();
    this->log_bitten_times_m2h[h_index] = 0;
    this->log_bitten_times_h2m[h_index] = 0;
}

void HumanManager::init_gender(const float male_probability) {

    float* random_numbers = new float[this->sum_num_humans];
    util::set_random_numbers_uniform(random_numbers, this->sum_num_humans);

    for (int hh = 0; hh < this->sum_num_humans; hh++) {
        this->gender_of[hh] = (random_numbers[hh] < male_probability ? 'm' : 'f');
    }

    delete[] random_numbers;

}

void HumanManager::print_gender_distribution() const{
    print_distribution<char>(this->gender_of, this->sum_num_humans);
}



void HumanManager::init_mobility(const float* village_init_mobility_percentages){
// // float mobpop = float(vregion[i]->v[6]);

// // mob[i]=mobpop;

// // float mobile_pop = mob[i];

// //    G_mobility_static = 125.0/365.0;    // static pop probability of moving to another village
// //    G_static_return = 1.0/3.0;         // static pop average length stay in another village
// //    G_mobility_mobile = 125.0/365.0;     // seasonal migrant prob of moving to another village
// //    G_mobile_return = 1.0/3.0;         // seasonal migrant average length stay in another village
    

// // float randomEvent2 = ran2();
// // if(randomEvent2 <= mobile_pop){
// //     float randomEvent3 = ran2();
// //         if(randomEvent3 <= G_temp_mobile){
// //             G_demes[ind]->mobile = TEMP;}
// //         else{G_demes[ind]->mobile = MOBILE;}
// // }
// // else { G_demes[ind]->mobile = STATIC;}

    const float G_temp_mobile = 0.2; // G_temp_mobile not initialised in legacy version

    float* random_numbers_1 = new float[this->sum_num_humans];
    util::set_random_numbers_uniform(random_numbers_1, this->sum_num_humans);
    float* random_numbers_2 = new float[this->sum_num_humans];
    util::set_random_numbers_uniform(random_numbers_2, this->sum_num_humans);
    
    for(int hh = 0; hh < this->sum_num_humans; hh++){
        if(random_numbers_1[hh] <= village_init_mobility_percentages[this->home_village_of[hh]]){
            if(random_numbers_2[hh] <= G_temp_mobile){
                this->mobility_of[hh] = common::HumanMobilityType::kTemporary;
            } else {
                this->mobility_of[hh] = common::HumanMobilityType::kMobile;
            }
        } else {
            this->mobility_of[hh] = common::HumanMobilityType::kStatic;
        }
    }
    delete[] random_numbers_1;
    delete[] random_numbers_2;
}

void HumanManager::init_itn(float itn_probability){

    if (itn_probability == 0.0) {
        for(int hh = 0; hh < this->sum_num_humans; hh++){
            this->if_itn_of[hh] = false;
        }
    } else {    
        float* random_numbers = new float[this->sum_num_humans];
        util::set_random_numbers_uniform(random_numbers, this->sum_num_humans);

        this->init_itn_kernel(
                itn_probability,
                this->if_itn_of,
                random_numbers,
                (uint) this->sum_num_humans
            );
        delete[] random_numbers;
    }

}
void HumanManager::init_itn_kernel( //TODO: template this?
        float probability,
        bool* data,
        float* random_numbers,
        uint data_size
    ){

    for(uint hh = 0; hh < data_size; hh++){
        data[hh] = (random_numbers[hh] < probability);
    }
}

float HumanManager::func_gen_attractiveness() const {
    float aa = 0.0;
    // while(true) {
        aa = util::get_rand_normal(this->kAttractiveness_mean, this->kAttractiveness_sd);
        // if(aa>0) {
        if(aa<0) {
            aa=0;
        }
    // }
        // if (util::get_rand_uniform() < 0.8) {
        //     aa = 1.0;
        // } else {
        //     aa = 0.0;
        // }
    return aa;
}

void HumanManager::init_attractiveness() {
    assert(this->sum_num_humans == static_cast<int>(this->attractiveness_of.size()));
    if (this->kAttractiveness_enabled) {
        std::cout << "HM: initialising attractiveness...";
        for ( float& aa : this->attractiveness_of) {
            aa = this->func_gen_attractiveness();
            // while(true) {
            //     aa = 
            //     if(aa>0) {
            //         break;
            //     }
            // }
        }
        std::cout << " done\n";
    }
}


// void HumanManager::init_parasites(){
//     // this->parasite_offset = new int[this->sum_num_humans];

//     // // std::mt19937 gen(std::random_device{}());

//     // // std::discrete_distribution<> disc_dist(weight_list.begin(), weight_list.end());

//     // for (int hh = 0; hh < this->sum_num_humans; hh++) {
//     //     this->parasite_reg.push_back({  common::ParasiteType::kRa,
//     //                                     common::StageName::kBlood,
//     //                                     0,0}
//     //         );
//     //     this->parasite_offset[hh] = hh;
//     // }

// }

// void HumanManager::update_bitten_by(float* village_bite_rate, float* village_ra_rate){
//     float* random_numbers_1 = new float[this->sum_num_humans];
//     float* random_numbers_2 = new float[this->sum_num_humans];
//     util::set_random_numbers_uniform(random_numbers_1, this->sum_num_humans);
//     util::set_random_numbers_uniform(random_numbers_2, this->sum_num_humans);

//     for (int hh = 0; hh < this->sum_num_humans; hh++){
//         if (random_numbers_1[hh] < village_bite_rate[this->at_village_of[hh]] &&
//             random_numbers_2[hh] < village_ra_rate[this->at_village_of[hh]] ){
//             this->bitten_by_of[hh] = common::ParasiteType::kRa;
//         } else {
//             this->bitten_by_of[hh] = common::ParasiteType::kNoParasite;
//         }

//     }
//     print_distribution<common::ParasiteType>(this->bitten_by_of, this->sum_num_humans);
// }

common::ParasiteType* HumanManager::get_bitten_by(){
    return this->bitten_by_of;
}

void HumanManager::reset_bitten_by(){
    for (int hh = 0; hh < this->sum_num_humans; hh++) {
        this->bitten_by_of[hh] = common::ParasiteType::kNoParasite;
    }
}

// void HumanManager::update_given_drug_by_at_village( const int index_village,
//                                                     const common::DrugName drug){
//     for (int hh = 0; hh < this->sum_num_humans; hh++){
//         if (this->at_village_of[hh] == index_village) {
//             this->given_drug_of[hh] = drug;
//         } else {
//             this->given_drug_of[hh] = common::DrugName::kNoDrug;
//         }
//     }
// }

// void HumanManager::give_drug_by_home_village(
//     const int village_index,
//     const common::DrugName drug,
//     const float probability
//     ) {
//     for (int hh = 0; hh < this->sum_num_humans; hh++) {
//         // if (this->home_village_of[])
//     }
// }


void HumanManager::clear_given_drug() {
    for (int hh = 0; hh < this->sum_num_humans; hh++){
        this->given_drug_of[hh] = common::DrugName::kNoDrug;
    }
}

void HumanManager::step_end_of_day(int time_step) {
    if (kReset_log_bitten_times_on_day == time_step) {
        this->reset_log_bitten_times();
    }
    this->clear_given_drug();
    // bitten_by is reset at the end of bsys.parasite_intake

    // int sum = 0;
    // for (auto& num : this->log_bitten_times_m2h) {
    //     sum+= num;
    // }
    // std::cout << "HM: sum=" << sum << "\n";
    // std::cin.ignore();
}


template <class T> // TODO: move this to util!
void HumanManager::print_distribution(T* data_array, int length){
    std::map<T, int> hist;

    for (int ii = 0; ii < length; ii++) {
        hist[data_array[ii]]++;
    }

    const int kScreen_width = 50;
    int population_per_star = length / hist.size() / kScreen_width;
    population_per_star = population_per_star > 0 ? population_per_star : 1;

    for (auto bin : hist) { //TODO auto const &bin
        std::cout << std::fixed
                  << std::setprecision(2) << std::setw(2) << bin.first << ' ' 
                  << std::setprecision(2) << std::setw(3) << (bin.second * 100.0 / length ) << "% "
                  << std::string(bin.second / population_per_star, '*') << ' ' << bin.second
                  << std::endl;
    }
}
template void HumanManager::print_distribution<bool>(bool* data_array, int length);
template void HumanManager::print_distribution<common::ParasiteType>(common::ParasiteType* data_array, int length);
template void HumanManager::print_distribution<common::DrugName>(common::DrugName* data_array, int length);


}