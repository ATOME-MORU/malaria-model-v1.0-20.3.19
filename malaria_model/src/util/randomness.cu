#include <stdio.h>
#include <stdlib.h>

#include <iostream>

#include <cuda.h>
#include <curand.h>

#include <vector>

#include "util/randomness.h"

namespace util {


#ifdef __NVCC__
int cuda_device = 0;
#endif

#define CUDA_CALL(x) do { if((x)!=cudaSuccess) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
    return;}} while(0)
    // return EXIT_FAILURE;}} while(0)
#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) { \
    printf("Error at %s:%d\n",__FILE__,__LINE__);\
    return;}} while(0)
    // return EXIT_FAILURE;}} while(0)

__global__ void filter_ones(float* number_array, int array_length) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    if (index < array_length) {
        if (number_array[index] >= 1) {
            number_array[index] -= (1.0-kEPS);
        }
    }
}


void set_random_numbers_uniform_curand(float* result_array, int array_length){
    static curandGenerator_t gen;
    static bool seeded(false);

    // CUDA_CALL(cudaSetDevice(1));

    static float* dev_random_number_array;
    static int dev_random_number_array_size(1);

    const int kThreads_per_block = 1024;

    CUDA_CALL(
        cudaSetDevice(cuda_device)
    );

    if (!seeded) {
        // CURAND_CALL(curandCreateGenerator(&gen,
        //  CURAND_RNG_PSEUDO_DEFAULT));

            // CURAND_CALL(curandCreateGenerator(&gen, 
            //     CURAND_RNG_PSEUDO_DEFAULT));
        CURAND_CALL(
            curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT)
        );
        // CURAND_CALL(curandSetPseudoRandomGeneratorSeed(gen, time(NULL)));
        CURAND_CALL(
            curandSetPseudoRandomGeneratorSeed(gen, time(NULL))
        );
        seeded = true;

        CUDA_CALL(
            cudaMalloc((void**) &dev_random_number_array, dev_random_number_array_size * sizeof(float))
        );
    }


    if (array_length > dev_random_number_array_size) {
        CUDA_CALL(
            cudaFree(dev_random_number_array)
        );
        dev_random_number_array_size = array_length;
        CUDA_CALL(
            cudaMalloc((void**) &dev_random_number_array, dev_random_number_array_size * sizeof(float))
        );
    }


    // std::cout << "GPU ramdom\n";
    // CUDA_CALL(cudaMalloc((void**) &dev_random_number_array, array_length * sizeof(float)));

    // CURAND_CALL(curandGenerateUniform(gen, dev_random_number_array, array_length));

    // CUDA_CALL(cudaMemcpy(random_numbers, dev_random_number_array, array_length * sizeof(float)), cudaMemcpyDeviceToHost);

    // CUDA_CALL(cudaFree(dev_random_number_array));

    CURAND_CALL(
        curandGenerateUniform(gen, dev_random_number_array, array_length)
    );

    filter_ones<<<(array_length+kThreads_per_block-1)/kThreads_per_block,kThreads_per_block>>>(dev_random_number_array, array_length);

    CUDA_CALL(
        cudaMemcpy(result_array, dev_random_number_array, array_length * sizeof(float), cudaMemcpyDeviceToHost)
    );
    // cudaFree(dev_random_number_array);

}

}