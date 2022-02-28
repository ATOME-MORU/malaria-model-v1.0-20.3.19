#ifndef RANDOMNESS_H
#define RANDOMNESS_H

//#include <vector>

namespace util {

const float kEPS = 1.2e-7;

#ifdef __NVCC__
extern int cuda_device;
#endif
// class RandomManager {

// public:
//     RandomManager

// }


const int kRandomness_buffer_length = 10000000; // 10000000 * float = 40MB

void config_uniform_buffer(bool if_buffer);

void config_uniform_method(int method_id);
void config_normal_method(int method_id);
void config_poisson_method(int method_id);
void config_discrete_method(int method_id);
void config_binomial_method(int method_id);



void set_random_numbers_uniform_boost_thread(float* result_array, int array_length);

float get_rand_uniform();

void set_random_numbers_uniform(float* random_numbers, int length);
void set_random_numbers_uniform_buffered(float* random_numbers, int length);
void set_random_numbers_uniform_switch(float* random_numbers, int length);

void set_random_numbers_uniform_boost(float* result_array, int array_length);
void set_random_numbers_uniform_std(float* result_array, int array_length);
#ifdef __NVCC__
void set_random_numbers_uniform_curand(float* result_array, int array_length);
#endif

float get_rand_normal(float mean, float sd);

void set_random_numbers_normal(float* random_numbers, int length);
// void set_random_numbers_normal_buffered(float* random_numbers, int length);
void set_random_numbers_normal_switch(float* random_numbers, int length);

void set_random_numbers_normal_boost(float* result_array, int array_length);
void set_random_numbers_normal_std(float* result_array, int array_length);


int get_rand_poisson(float mean);

void set_random_numbers_poisson(int* random_numbers, int length, float mean);
void set_random_numbers_poisson_boost(int* result_array, int array_length, float poisson_mean);
void set_random_numbers_poisson_std(int* result_array, int array_length, float poisson_mean);


template <typename Weight>
int get_rand_discrete(std::vector<Weight> weights);
template <typename Weight>
void set_random_numbers_discrete(int* random_numbers, int length, std::vector<Weight> weights);
template <typename Weight>
void set_random_numbers_discrete_boost(int* result_array, int array_length, std::vector<Weight> disc_weights);
template <typename Weight>
void set_random_numbers_discrete_std(int* result_array, int array_length, std::vector<Weight> disc_weights);

template <typename T>
void set_random_draw_from_queue(T* result_array, int array_length, std::vector<T> queue) ;



int get_rand_binomial(int num_events, float success_fraction);

void set_random_numbers_binomial(int* random_numbers, int length, int num_events, float success_fraction);
void set_random_numbers_binomial_boost(int* result_array, int array_length, int binomial_num_events, float binomial_success_fraction);
void set_random_numbers_binomial_std(int* result_array, int array_length, int binomial_num_events, float binomial_success_fraction);
void set_random_numbers_binomial_bernoulli_process(int* result_array, int array_length, int binomial_num_events, float binomial_success_fraction);

}


#endif