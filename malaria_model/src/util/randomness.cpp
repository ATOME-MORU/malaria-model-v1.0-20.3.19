#include <boost/random/random_device.hpp> //
#include <boost/random/mersenne_twister.hpp> // boost::mt19937
#include <boost/random/poisson_distribution.hpp>  // boost::poisson_distribution
#include <boost/random/uniform_real_distribution.hpp> // boost::random::uniform_real_distribution
#include <boost/random/normal_distribution.hpp> // boost::random::uniform_real_distribution
#include <boost/random/discrete_distribution.hpp> // boost::random::discrete_distribution
#include <boost/random/binomial_distribution.hpp>

#include <boost/random/variate_generator.hpp>  // boost::variate_generator

#include <iostream>

#include <random>

#include <cstring> // std::memcpy

#include "util/randomness.h"

namespace util {

bool randomness_uniform_buffer = true;

void config_uniform_buffer(bool if_buffer) {
    randomness_uniform_buffer = if_buffer;
}


#ifdef __NVCC__
int randomness_uniform_method = 2;
#else
int randomness_uniform_method = 0;
#endif
int randomness_normal_method = 0;
int randomness_poisson_method = 0;
int randomness_discrete_method = 0;
int randomness_binomial_method = 0;

void config_uniform_method(int method_id){
    randomness_uniform_method = method_id;
}
void config_normal_method(int method_id){
    randomness_normal_method = method_id;
}
void config_poisson_method(int method_id){
    randomness_poisson_method = method_id;
}
void config_discrete_method(int method_id){
    randomness_discrete_method = method_id;
}
void config_binomial_method(int method_id){
    randomness_binomial_method = method_id;
}


void set_random_numbers_uniform_boost_thread(float* result_array, int array_length){

    static boost::mt19937 b_gen;

    static bool seeded(false);
    if (!seeded) {
        b_gen.seed(time(NULL));
        seeded = true;
    }

    static boost::random::uniform_real_distribution<float> b_udist(0,1);
    static boost::variate_generator<boost::mt19937&, boost::random::uniform_real_distribution<float>> b_uniform_rvg(b_gen,b_udist);  
    
    float value = 0.0;
    for (int ii= 0; ii < array_length; ii++) {
        value = b_uniform_rvg();
        result_array[ii] = (value >= 1.0 ? 1.0-kEPS : value);
    }

}

float get_rand_uniform() {
    float rnd_number = 0.0;
    set_random_numbers_uniform(&rnd_number, 1);
    return rnd_number;
}

void set_random_numbers_uniform(float* random_numbers, int length){
    if (randomness_uniform_buffer) {
        set_random_numbers_uniform_buffered(random_numbers, length);
    } else {
        set_random_numbers_uniform_switch(random_numbers, length);
    }
}

void set_random_numbers_uniform_buffered(float* random_numbers, int length) {
    static int head(-1);
    static float random_number_buffer[kRandomness_buffer_length];

    if (head == -1) {
        set_random_numbers_uniform_switch(random_number_buffer, kRandomness_buffer_length);
        head = 0;
    }

    assert(length < kRandomness_buffer_length);
    assert(head <= kRandomness_buffer_length);

    if (length <= (kRandomness_buffer_length - head)) {

        std::memcpy(random_numbers, random_number_buffer + head, length * sizeof(float));
        head = head + length;
        
    } else {

        int section_1_length = kRandomness_buffer_length - head;
        int section_2_length = length - section_1_length;

        std::memcpy(random_numbers, random_number_buffer + head, section_1_length * sizeof(float));

        set_random_numbers_uniform_switch(random_number_buffer, kRandomness_buffer_length);
        head = 0;

        std::memcpy(random_numbers + section_1_length, random_number_buffer + head, section_2_length * sizeof(float));
        head = section_2_length;

    }

}

void set_random_numbers_uniform_switch(float* random_numbers, int length){
    switch(randomness_uniform_method) {
        case 0: set_random_numbers_uniform_boost(random_numbers, length); break;
        case 1: set_random_numbers_uniform_std(random_numbers, length); break;
#ifdef __NVCC__
        case 2: set_random_numbers_uniform_curand(random_numbers, length); break;
#endif
        default: set_random_numbers_uniform_boost(random_numbers, length); break;
    }
}
void set_random_numbers_uniform_boost(float* result_array, int array_length){

    // std::cout << "boost uniform\n";

    static boost::mt19937 b_gen;

    static bool seeded(false);
    if (!seeded) {
        b_gen.seed(time(NULL));
        seeded = true;
    }


    static boost::random::uniform_real_distribution<float> b_udist(0,1);
    static boost::variate_generator<boost::mt19937&, boost::random::uniform_real_distribution<float>> b_uniform_rvg(b_gen,b_udist);  
    
    float value = 0.0;
    for (int ii= 0; ii < array_length; ii++) {
        value = b_uniform_rvg();
        result_array[ii] = (value >= 1.0 ? 1.0-kEPS : value);
    }

}
void set_random_numbers_uniform_std(float* result_array, int array_length){

    // std::cout << "std uniform:\n";

    // static std::random_device rd;  //Will be used to obtain a seed for the random number engine
    static std::mt19937 gen(std::random_device{}()); //Standard mersenne_twister_engine seeded with rd()
    static std::uniform_real_distribution<> dis(0.0, 1.0);

    float value = 0.0;

    for( int ii= 0; ii < array_length; ii++ ){
        value = dis(gen);
        result_array[ii] = (value >= 1.0 ? 1.0-kEPS : value);        
    }

}

float get_rand_normal(float mean, float sd) {
    float rnd_number = 0;
    set_random_numbers_normal(&rnd_number, 1);
    return (rnd_number * sd + mean); // convertion from standard normal distribution
}

void set_random_numbers_normal(float* random_numbers, int length){
    set_random_numbers_normal_switch(random_numbers, length);
}
void set_random_numbers_normal_switch(float* random_numbers, int length) {
    switch(randomness_normal_method) {
        case 0: set_random_numbers_normal_boost(random_numbers, length); break;
        case 1: set_random_numbers_normal_std(random_numbers, length); break;
        default: set_random_numbers_normal_boost(random_numbers, length); break;
    }
}
void set_random_numbers_normal_boost(float* result_array, int array_length) {
    static boost::mt19937 b_n_gen(time(NULL));
    static boost::variate_generator< boost::mt19937&, boost::normal_distribution<float> > b_n_rvg(b_n_gen, boost::random::normal_distribution<float>(0,1));

    // std::cout << "boost normal\n";

    for (int ii = 0; ii < array_length; ii++) {
        result_array[ii] = b_n_rvg();
    }
}
void set_random_numbers_normal_std(float* result_array, int array_length) {
    static std::mt19937 gen(std::random_device{}());
    // static std::mt19937 gen(time(NULL));
    static std::normal_distribution<> dis(0,1);
    
    for( int ii= 0; ii < array_length; ii++ ){
        result_array[ii] = dis(gen);
    }
}


int get_rand_poisson(float mean) {
    int rnd_number = 0;
    set_random_numbers_poisson(&rnd_number, 1, mean);
    return rnd_number;
}

void set_random_numbers_poisson(int* random_numbers, int length, float mean){
    switch(randomness_poisson_method) {
        case 0: set_random_numbers_poisson_boost(random_numbers, length, mean); break;
        case 1: set_random_numbers_poisson_std(random_numbers, length, mean); break;
        default: set_random_numbers_poisson_boost(random_numbers, length, mean); break;
    }
}
void set_random_numbers_poisson_boost(int* result_array, int array_length, float poisson_mean){
    static boost::mt19937 b_p_gen;

    static bool seeded(false);
    if (!seeded) {
        b_p_gen.seed(time(NULL));
        seeded = true;
    }

    boost::random::poisson_distribution<int> b_pdist(poisson_mean);
    boost::variate_generator<boost::mt19937&, boost::poisson_distribution<int>> b_p_rvg(b_p_gen, b_pdist);

    for (int ii = 0; ii < array_length; ii++){
        result_array[ii] = b_p_rvg();
    }

}
void set_random_numbers_poisson_std(int* result_array, int array_length, float poisson_mean){
    // static std::default_random_engine gen;
    static std::mt19937 gen(std::random_device{}());
    std::poisson_distribution<int> dist(poisson_mean);

    for (int ii = 0; ii < array_length; ii++){
        result_array[ii] = dist(gen);
    }

}

template <typename Weight>
int get_rand_discrete(std::vector<Weight> weights){
    int rnd_number = 0;
    set_random_numbers_discrete<Weight>(&rnd_number, 1, weights);
    return rnd_number;
}
template int get_rand_discrete<int>(std::vector<int> weights);
template int get_rand_discrete<float>(std::vector<float> weights);

template <typename Weight>
void set_random_numbers_discrete(int* random_numbers, int length, std::vector<Weight> weights){
    switch (randomness_discrete_method) {
        case 0: set_random_numbers_discrete_boost<Weight>(random_numbers, length, weights); break;
        case 1: set_random_numbers_discrete_std<Weight>(random_numbers, length, weights); break;
        default: set_random_numbers_discrete_boost<Weight>(random_numbers, length, weights); break;
    }
}
template void set_random_numbers_discrete<int>(int* random_numbers, int length, std::vector<int> weights);
template void set_random_numbers_discrete<float>(int* random_numbers, int length, std::vector<float> weights);
template <typename Weight>
void set_random_numbers_discrete_boost(int* result_array, int array_length, std::vector<Weight> disc_weights){
    static boost::mt19937 b_d_gen;

    static bool seeded(false);
    if (!seeded) {
        b_d_gen.seed(time(NULL));
        seeded = true;
    }

    boost::random::discrete_distribution<int, Weight> b_d_dist(disc_weights.begin(), disc_weights.end());
    boost::variate_generator<boost::mt19937&, boost::random::discrete_distribution<int, Weight>> b_d_rvg(b_d_gen, b_d_dist);
    
    for (int ii = 0; ii < array_length; ii++) {
        result_array[ii] = b_d_rvg();
    }
}
template void set_random_numbers_discrete_boost<int>(int* result_array, int array_length, std::vector<int> disc_weights);
template void set_random_numbers_discrete_boost<float>(int* result_array, int array_length, std::vector<float> disc_weights);
template <typename Weight>
void set_random_numbers_discrete_std(int* result_array, int array_length, std::vector<Weight> disc_weights){
    static std::mt19937 gen(std::random_device{}());
    std::discrete_distribution<int> dist(disc_weights.begin(), disc_weights.end());

    for (int ii = 0; ii < array_length; ii++) {
        result_array[ii] = dist(gen);
    }
}
template void set_random_numbers_discrete_std<int>(int* result_array, int array_length, std::vector<int> disc_weights);
template void set_random_numbers_discrete_std<float>(int* result_array, int array_length, std::vector<float> disc_weights);

template <typename T>
void set_random_draw_from_queue(T* result_array, int array_length, std::vector<T> queue) {
    float* uniform_random_array;
    uniform_random_array = new float[array_length];
    set_random_numbers_uniform(uniform_random_array, array_length);

    for (int ii = 0; ii < array_length; ii++) {
        result_array[ii] = queue.at(static_cast<int>(uniform_random_array[ii] * queue.size()));
    }

    delete[] uniform_random_array;

}
template void set_random_draw_from_queue<int>(int* result_array, int array_length, std::vector<int> queue);


int get_rand_binomial(int num_events, float success_fraction) {
    int temp = 0;
    set_random_numbers_binomial(&temp, 1, num_events, success_fraction);
    return temp;
}

void set_random_numbers_binomial(int* random_numbers, int length, int num_events, float success_fraction) {
    
    switch(randomness_binomial_method) {
        case 0: set_random_numbers_binomial_boost(random_numbers, length, num_events, success_fraction); break;
        case 1: set_random_numbers_binomial_std(random_numbers, length, num_events, success_fraction); break;
        case 2: set_random_numbers_binomial_bernoulli_process(random_numbers, length, num_events, success_fraction); break;
        default: set_random_numbers_binomial_boost(random_numbers, length, num_events, success_fraction); break;
    }
}
void set_random_numbers_binomial_boost(int* result_array, int array_length, int binomial_num_events, float binomial_success_fraction) {

    static boost::mt19937 b_gen;

    static bool seeded(false);
    if (!seeded) {
        b_gen.seed(time(NULL));
        seeded = true; 
    }

    boost::random::binomial_distribution<> b_b_dist(binomial_num_events,binomial_success_fraction);
    boost::variate_generator<boost::mt19937&, boost::random::binomial_distribution<>> b_binomial_rvg(b_gen,b_b_dist);  

    for (int ii= 0; ii < array_length; ii++) {
        result_array[ii] = b_binomial_rvg();
    }
}
void set_random_numbers_binomial_std(int* result_array, int array_length, int binomial_num_events, float binomial_success_fraction){
    static std::mt19937 gen(std::random_device{}()); //Standard mersenne_twister_engine seeded with rd()
    std::binomial_distribution<> dis(binomial_num_events,binomial_success_fraction);

    for( int ii= 0; ii < array_length; ii++ ){
        result_array[ii] = dis(gen);        
    }
}
void set_random_numbers_binomial_bernoulli_process(int* result_array, int array_length, int binomial_num_events, float binomial_success_fraction){
    int result_ii = 0;
    float* uniform_random_array;
    uniform_random_array = new float[binomial_num_events];
    for( int ii= 0; ii < array_length; ii++ ){

        set_random_numbers_uniform(uniform_random_array, binomial_num_events);

        result_ii = 0;
        for (int ee = 0; ee < binomial_num_events; ee++) {
            if (uniform_random_array[ee] < binomial_success_fraction) {
                result_ii++;
            }
        }
        result_array[ii] = result_ii;

    }
    delete[] uniform_random_array;
}

}

