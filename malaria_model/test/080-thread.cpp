#include <thread>
#include "third_party/catch2/catch.hpp"
#include <iostream>
#include <iomanip>

#include <string>
#include <cmath>

#include <algorithm> // std::fill_n

#include "util/randomness.h"

/**
    Test function util::set_random_numbers_uniform_boost_thread for:
    1. Each call produces a different series of uniformly distributed
    random numbers [0-1).
    2. The mean value of the random numbers produced by this function 
    is within a tolerance of 0.5.
*/
TEST_CASE( "80.1: Thread: uniform random", "[thread:uniform_random]" ) {

    const int kNum_stars = 100;

    const int kNum_tests = 100;
    float test_mean_array[kNum_tests];
    std::fill_n(test_mean_array, kNum_tests, 0);
    float test_mean_array_mean = 0.0;
    const float kMean_tolerance = 0.01;
    const float kTarget_mean = 0.5;

    const int kArray_length = 100000;
    float rnd_nmbr_array[kArray_length];
    std::fill_n(rnd_nmbr_array, kArray_length, 0);

    const int kResult_counter_size = 10;
    int result_counter[kResult_counter_size];
    std::fill_n(result_counter, kResult_counter_size, 0);

    std::thread thrd_rnd_nmbrs;

    for (int tt = 0; tt < kNum_tests; tt++) {

        thrd_rnd_nmbrs = std::thread(
            util::set_random_numbers_uniform_boost_thread,
            rnd_nmbr_array,
            kArray_length
        );

        thrd_rnd_nmbrs.join();


        for (int ii = 0; ii < kArray_length; ii++) {
            result_counter[int(rnd_nmbr_array[ii]*10)]++;
            test_mean_array[tt] += rnd_nmbr_array[ii];
        }
        test_mean_array[tt] /= kArray_length;
        test_mean_array_mean += test_mean_array[tt];

        if (tt > 0){
            REQUIRE(test_mean_array[tt] != test_mean_array[tt-1]);
        }

        if (tt < 3) {
            std::cout << "\nExample results from test #" << tt << "/" << kNum_tests << ":\n";
            for(int cc = 0; cc < kResult_counter_size; cc++) {
                std::cout << std::fixed << std::setprecision(1)
                            << float(cc)/10 << ":"
                            << std::string(result_counter[cc]*kNum_stars/kArray_length, '*')
                            << " " << result_counter[cc] << "\n";
            }
        }

        std::fill_n(result_counter, kResult_counter_size, 0);

    }
    test_mean_array_mean /= kNum_tests;

    // std::cout << std::setprecision(5) << test_mean_array_mean << "\n";
    REQUIRE(std::abs(test_mean_array_mean-kTarget_mean) < kTarget_mean*kMean_tolerance);


    // The following won't work because each thread will need its own boost::mt19937 object

    // const int kNum_threads = 2;
    // float thread_mean[kNum_threads];
    // std::fill_n(thread_mean, kNum_threads, 0);
    // const int kNum_numbers_per_thread = kArray_length / kNum_threads;
    // std::thread thrd_rnd_nmbrs_array[kNum_threads];
    // for (int tt = 0; tt < kNum_threads; tt++) {
    //     thrd_rnd_nmbrs_array[0] = std::thread(
    //         util::set_random_numbers_uniform_boost_thread,
    //         rnd_nmbr_array + tt*kNum_numbers_per_thread,
    //         kNum_numbers_per_thread
    //     );
    // }

    // test_mean_array_mean = 0;
    // for (int tt = 0; tt < kNum_threads; tt++) {
    //     thrd_rnd_nmbrs_array[tt].join();
    //     for (int ii = 0; ii < kNum_numbers_per_thread; ii++) {
    //         thread_mean[tt] += rnd_nmbr_array[tt*kNum_numbers_per_thread + ii];
    //     }
    //     thread_mean[tt] /= kNum_numbers_per_thread;
    //     if (tt > 0) {
    //         REQUIRE(thread_mean[tt] != thread_mean[tt-1]);
    //     }
    //     test_mean_array_mean += thread_mean[tt];
    // }
    // test_mean_array_mean /= kNum_threads;

    // std::cout << test_mean_array_mean << "\n";


}