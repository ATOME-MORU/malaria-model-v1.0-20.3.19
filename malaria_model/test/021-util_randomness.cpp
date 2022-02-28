#include "third_party/catch2/catch.hpp"
#include <iostream>
#include <iomanip>

#include <algorithm> // std::fill_n, std::max

#include <cmath>

#include <string>
#include <vector>

#include "util/numerical_poisson.h"
#include "util/randomness.h"
#include "util/util.h"

// #include "common_test.h"


TEST_CASE( "21.1: Randomness: uniform random numbers", "[util:rand:uniform]" ) {

    std::vector<std::string> method_names;
    method_names.push_back("Boost::random::uniform_real_distribution");
    method_names.push_back("std::uniform_real_distribution");
#ifdef __NVCC__
    method_names.push_back("curand::hostAPI");
#endif

    const int kNum_tests = 10000;
    const int kNum_stars = 100;
    const int kArray_size = 10000;
    float uniform_random_array[kArray_size] = {0.0};
    const int kResult_counter_size = 10;
    int result_counter[kResult_counter_size] = {0};

    bool result_all_less_than_one = true;
    bool result_no_consequtive_identical_values = true;

    double wall_time_begin = 0.0;
    double wall_time_end = 0.0;

    for (int bb = 0; bb < 2; bb++){

        if (bb == 0) {
            util::config_uniform_buffer(false);
            std::cout << "No buffer:\n";
        } else {
            util::config_uniform_buffer(true);
            std::cout << "\nWith buffer ("
                        << util::kRandomness_buffer_length 
                        << " x float = "
                        << util::kRandomness_buffer_length/1000000 *sizeof(float)
                        << " MB ):\n";
        }

        for (int mm = 0; mm < static_cast<int>(method_names.size()); mm++){
            
            util::config_uniform_method (mm); //0-boost; 1-std
            util::set_random_numbers_uniform(uniform_random_array, kArray_size);
            
            std::cout << "Result using method " << method_names.at(mm) << ":\n";
            for (int ii = 0; ii < kArray_size; ii++) {
                result_counter[int(uniform_random_array[ii]*10)]++;
                if (uniform_random_array[ii] >=1 ){
                    result_all_less_than_one = false;
                }
            }
            for (int cc = 0; cc < kResult_counter_size; cc++) {
                std::cout << std::fixed << std::setprecision(1)
                            << float(cc)/10 << ": "
                            << std::string(result_counter[cc]*kNum_stars/kArray_size, '*')
                            << " " << result_counter[cc] << "\n";
            }
            REQUIRE(result_all_less_than_one);

            std::fill(result_counter, result_counter + kResult_counter_size, 0);
            result_all_less_than_one = true;

            // Randomness
            for (int ii = 0; ii < kArray_size; ii++) {
                util::set_random_numbers_uniform(uniform_random_array + ii, 1);
                if (ii > 0) {
                    if(uniform_random_array[ii] == uniform_random_array[ii-1]) {
                        result_no_consequtive_identical_values = false;
                    }
                }
            }
            REQUIRE(result_no_consequtive_identical_values);


            std::cout << "Time cost for "
                        << kNum_tests << " draws of "
                        << kArray_size << " numbers: " << std::flush;

            wall_time_begin = util::get_wall_time();
            for (int tt = 0; tt < kNum_tests; tt++ ){
                util::set_random_numbers_uniform(uniform_random_array, kArray_size);
            }
            wall_time_end = util::get_wall_time();

            std::cout << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";


            std::cout << "Time cost for "
                        << kNum_tests << " draws of "
                        << 1 << " number: " << std::flush;

            wall_time_begin = util::get_wall_time();
            for (int tt = 0; tt < kNum_tests; tt++ ){
                util::set_random_numbers_uniform(uniform_random_array, 1);
            }
            wall_time_end = util::get_wall_time();

            std::cout << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";


            std::cout << "Time cost for "
                        << 1 << " draws of "
                        << kArray_size << " number: " << std::flush;

            wall_time_begin = util::get_wall_time();
            // for (int tt = 0; tt < kNum_tests; tt++ ){
                util::set_random_numbers_uniform(uniform_random_array, kArray_size);
            // }
            wall_time_end = util::get_wall_time();

            std::cout << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";

        }

    } // if_buffer

}

TEST_CASE( "21.1.1: Randomness: get_rand_uniform", "[util:rand:get_rand_uniform]" ) {

    std::vector<std::string> method_names;
    method_names.push_back("Boost::random::uniform_real_distribution");
    method_names.push_back("std::uniform_real_distribution");
#ifdef __NVCC__
    method_names.push_back("curand::hostAPI");
#endif

    const int kNum_tests = 10000;
    const int kNum_stars = 100;
    const int kArray_size = 10000;
    // float uniform_random_array[kArray_size] = {0.0};
    // float rnd_nmb = 0.0;
    const int kResult_counter_size = 10;
    int result_counter[kResult_counter_size] = {0};

    bool result_all_less_than_one = true;
    bool result_no_consequtive_identical_values = true;

    double wall_time_begin = 0.0;
    double wall_time_end = 0.0;

    for (int bb = 0; bb < 2; bb++){

        if (bb == 0) {
            util::config_uniform_buffer(false);
            std::cout << "No buffer:\n";
        } else {
            util::config_uniform_buffer(true);
            std::cout << "\nWith buffer ("
                        << util::kRandomness_buffer_length 
                        << " x float = "
                        << util::kRandomness_buffer_length/1000000 *sizeof(float)
                        << " MB ):\n";
        }

        for (int mm = 0; mm < static_cast<int>(method_names.size()); mm++){
            
            util::config_uniform_method (mm); //0-boost; 1-std
            // util::set_random_numbers_uniform(uniform_random_array, kArray_size);
            
            std::cout << "Result using method " << method_names.at(mm) << ":\n";
            for (int ii = 0; ii < kArray_size; ii++) {
                float rnd_nmb = util::get_rand_uniform();
                result_counter[int(rnd_nmb*10)]++;
                if (rnd_nmb >=1 ){
                    result_all_less_than_one = false;
                }
            }
            for (int cc = 0; cc < kResult_counter_size; cc++) {
                std::cout << std::fixed << std::setprecision(1)
                            << float(cc)/10 << ": "
                            << std::string(result_counter[cc]*kNum_stars/kArray_size, '*')
                            << " " << result_counter[cc] << "\n";
            }
            REQUIRE(result_all_less_than_one);

            std::fill(result_counter, result_counter + kResult_counter_size, 0);
            result_all_less_than_one = true;

            // Randomness
            for (int ii = 0; ii < kArray_size; ii++) {
                float rnd_nmb = util::get_rand_uniform();
                if(rnd_nmb == util::get_rand_uniform()) {
                    result_no_consequtive_identical_values = false;
                }
            }
            REQUIRE(result_no_consequtive_identical_values);


            if (bb == 1) {            
                std::cout << "Time cost for "
                            << kNum_tests << " draws of "
                            << kArray_size << " numbers: " << std::flush;

                wall_time_begin = util::get_wall_time();
                for (int tt = 0; tt < kNum_tests; tt++ ){
                    for (int tt1 = 0; tt1 < kArray_size; tt1++ ){
                        util::get_rand_uniform();
                    }
                }
                wall_time_end = util::get_wall_time();

                std::cout << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";
            }


            std::cout << "Time cost for "
                        << kNum_tests << " draws of "
                        << 1 << " number: " << std::flush;

            wall_time_begin = util::get_wall_time();
            for (int tt = 0; tt < kNum_tests; tt++ ){
                util::get_rand_uniform();
            }
            wall_time_end = util::get_wall_time();

            std::cout << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";


            std::cout << "Time cost for "
                        << 1 << " draws of "
                        << kArray_size << " number: " << std::flush;

            wall_time_begin = util::get_wall_time();
            for (int tt = 0; tt < kArray_size; tt++ ){
                util::get_rand_uniform();
            }
            wall_time_end = util::get_wall_time();

            std::cout << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";

        }

    } // if_buffer

}

TEST_CASE( "21.2: Randomness: poisson random numbers", "[util:rand:poisson]" ) {

    std::vector<std::string> method_names;
    method_names.push_back("Boost::random::poisson_distribution");
    method_names.push_back("std::poisson_distribution");

    method_names.push_back("numerical_poissonHighPrecision (selected tests only)");

    const int kNum_tests = 10000;
    const int kNum_tests_randomness = 10;
    const int kNum_stars = 100;
    const int kArray_size = 10000;
    const float poisson_mean = 4.5;
    const float tolerence = poisson_mean/50;

    int poisson_random_array[kArray_size] = {0};
    const int kResult_counter_size = 20;
    int result_counter[kResult_counter_size] = {0};

    bool result_all_identical_values = true;

    float result_mean = 0.0;
    float result_mean_prev = 0.0;

    double wall_time_begin = 0.0;
    double wall_time_end = 0.0;

    for (int mm = 0; mm < 3; mm++){
        
        // Illustrative test
        if (mm == 2) {
            for (int ii = 0; ii < kArray_size; ii++) {
                poisson_random_array[ii] = int(poissonHighPrecision(poisson_mean));
            }
        } else {
            
            util::config_poisson_method (mm); //0-boost; 1-std
            util::set_random_numbers_poisson(poisson_random_array, kArray_size, poisson_mean);
        }
        
        std::cout << "Result using method " << method_names.at(mm) << ":\n";
        for (int ii = 0; ii < kArray_size; ii++) {
            if (poisson_random_array[ii]<kResult_counter_size){
                result_counter[poisson_random_array[ii]]++;
            }
            result_mean += float(poisson_random_array[ii])/float(kArray_size);
        }
        for (int cc = 0; cc < kResult_counter_size; cc++) {
            std::cout << std::fixed << std::setprecision(0)
                        << cc << ": "
                        << std::string(result_counter[cc]*kNum_stars/kArray_size, '*')
                        << " " << result_counter[cc] << "\n";
        }
        std::cout << std::fixed << std::setprecision(4) << "Mean: " << result_mean
                    << " (target: " << poisson_mean << " +/-" << tolerence << ")\n";

        REQUIRE(std::abs(result_mean - poisson_mean) < tolerence);


        // Randomness and Mean tests
        for (int ii = 0; ii < kArray_size; ii++) {
            util::set_random_numbers_poisson(poisson_random_array + ii, 1, poisson_mean);
            if (ii > 0) {
                if(poisson_random_array[ii] != poisson_random_array[ii-1]) {
                    result_all_identical_values = false;
                }
            }
        }
        REQUIRE(!result_all_identical_values);


        std::fill(result_counter, result_counter + kResult_counter_size, 0);
        result_mean = 0.0;

        std::cout << "Repeating for " << kNum_tests_randomness << " times to check randomness... " << std::flush;
        for (int tt = 0; tt < kNum_tests_randomness; tt++ ){
            util::set_random_numbers_poisson(poisson_random_array, kArray_size, poisson_mean);

            result_mean = 0.0;
            for (int ii = 0; ii < kArray_size; ii++) {
                result_mean += float(poisson_random_array[ii])/float(kArray_size);
            }

            REQUIRE(std::abs(result_mean - poisson_mean) < tolerence);

            if (tt == 0) {
                result_mean_prev = result_mean;
            } else {
                REQUIRE(result_mean != result_mean_prev);
                result_mean_prev = result_mean;
            }
        }
        std::cout << "done.\n";

        std::fill(result_counter, result_counter + kResult_counter_size, 0);
        result_mean = 0.0;
        result_mean_prev = 0.0;


        // Performance
        std::cout << "Time cost for "
                    << kNum_tests << " draws of "
                    << kArray_size << " numbers: " << std::flush;

        if (mm == 2) {
            wall_time_begin = util::get_wall_time();
            for (int tt = 0; tt < kNum_tests; tt++ ){
                for (int ii = 0; ii < kArray_size; ii++) {
                    poisson_random_array[ii] = int(poissonHighPrecision(poisson_mean));
                }
            }
            wall_time_end = util::get_wall_time();

        } else {

            wall_time_begin = util::get_wall_time();
            for (int tt = 0; tt < kNum_tests; tt++ ){
                util::set_random_numbers_poisson(poisson_random_array, kArray_size, poisson_mean);
            }
            wall_time_end = util::get_wall_time();
        }
        
        std::cout << wall_time_end - wall_time_begin << " seconds.\n";

    }

}

TEST_CASE( "21.3: Randomness: discrete random numbers", "[util:rand:discrete]" ) {

    std::vector<std::string> method_names;
    method_names.push_back("Boost::random::discrete_distribution");
    method_names.push_back("std::discrete_distribution");

    const int kNum_tests = 10000;
    const int kNum_tests_randomness = 10;

    const int kNum_stars = 100;
    const int kArray_size = 10000;

    const int kNum_weights = 3;
    const int kWeight_0 = 10;
    const int kWeight_1 = 20;
    const int kWeight_2 = 30;

    std::vector<int> weights;
    weights.push_back(kWeight_0);
    weights.push_back(kWeight_1);
    weights.push_back(kWeight_2);

    const float target_mean = float(kWeight_0*0 + kWeight_1 + kWeight_2*2)/(kWeight_0 + kWeight_1 + kWeight_2);

    const float tolerence = target_mean/50;

    int discrete_random_array[kArray_size] = {0};
    int result_counter[kNum_weights] = {0};

    bool result_all_identical_values = true;

    float result_mean = 0.0;
    float result_mean_prev = 0.0;

    double wall_time_begin = 0.0;
    double wall_time_end = 0.0;

    for (int mm = 0; mm < 2; mm++){
        
        util::config_discrete_method (mm); //0-boost; 1-std
        util::set_random_numbers_discrete<int>(discrete_random_array, kArray_size, weights);
        
        std::cout << "Result using method " << method_names.at(mm) << ":\n";
        for (int ii = 0; ii < kArray_size; ii++) {
            result_counter[discrete_random_array[ii]]++;
            result_mean += float(discrete_random_array[ii])/float(kArray_size);
        }
        for (int cc = 0; cc < kNum_weights; cc++) {
            std::cout << std::fixed << std::setprecision(0)
                        << cc << "(" << weights.at(cc) << ") : "
                        << std::string(result_counter[cc]*kNum_stars/kArray_size, '*')
                        << " " << result_counter[cc] << "\n";
        }
        std::cout << std::fixed << std::setprecision(4) << "Mean: " << result_mean
                    << " (target: " << target_mean
                    << " +/-" << tolerence << ")\n";

        REQUIRE(std::abs(result_mean - target_mean) < tolerence);

        std::fill(result_counter, result_counter + kNum_weights, 0);
        result_mean = 0.0;


        // Randomness and Mean tests
        std::cout << "Time cost for "
                    << kArray_size << " calls to draw "
                    << 1 << " number: " << std::flush;

        wall_time_begin = util::get_wall_time();
        for (int ii = 0; ii < kArray_size; ii++) {
            util::set_random_numbers_discrete<int>(discrete_random_array + ii, 1, weights);
            if (ii > 0) {
                if(discrete_random_array[ii] != discrete_random_array[ii-1]) {
                    result_all_identical_values = false;
                }
            }
        }
        wall_time_end = util::get_wall_time();
        std::cout << wall_time_end - wall_time_begin << " seconds.\n";
        
        REQUIRE(!result_all_identical_values);

        std::cout << "Repeating for " << kNum_tests_randomness << " times to check randomness... " << std::flush;
        for (int tt = 0; tt < kNum_tests_randomness; tt++ ){
            util::set_random_numbers_discrete<int>(discrete_random_array, kArray_size, weights);

            result_mean = 0.0;
            for (int ii = 0; ii < kArray_size; ii++) {
                result_mean += float(discrete_random_array[ii])/float(kArray_size);
            }

            REQUIRE(std::abs(result_mean - target_mean) < tolerence);

            if (tt == 0) {
                result_mean_prev = result_mean;
            } else {
                REQUIRE(result_mean != result_mean_prev);
                result_mean_prev = result_mean;
            }
        }
        std::cout << "done.\n";

        result_mean = 0.0;
        result_mean_prev = 0.0;


        // Performance
        std::cout << "Time cost for "
                    << kNum_tests << " draws of "
                    << kArray_size << " numbers: " << std::flush;

        wall_time_begin = util::get_wall_time();
        for (int tt = 0; tt < kNum_tests; tt++ ){
            util::set_random_numbers_discrete<int>(discrete_random_array, kArray_size, weights);
        }
        wall_time_end = util::get_wall_time();

        std::cout << wall_time_end - wall_time_begin << " seconds.\n";

    }

    std::vector<int> queue;
    for(uint jj = 0; jj < weights.size(); jj++){
        for (int ii = 0; ii < weights.at(jj); ii++) {
            queue.push_back(jj);
        }
    }

    util::set_random_draw_from_queue<int>(discrete_random_array, kArray_size, queue);

    std::cout << "Result using (custom) method util::set_random_draw_from_queue:\n";
    for (int ii = 0; ii < kArray_size; ii++) {
        result_counter[discrete_random_array[ii]]++;
        result_mean += float(discrete_random_array[ii])/float(kArray_size);
    }
    for (int cc = 0; cc < kNum_weights; cc++) {
        std::cout << std::fixed << std::setprecision(0)
                    << cc << "(" << weights.at(cc) << ") : "
                    << std::string(result_counter[cc]*kNum_stars/kArray_size, '*')
                    << " " << result_counter[cc] << "\n";
    }
    std::cout << std::fixed << std::setprecision(4) << "Mean: " << result_mean
                << " (target: " << target_mean
                << " +/-" << tolerence << ")\n";

    REQUIRE(std::abs(result_mean - target_mean) < tolerence);

    // Randomness and Mean tests
    std::cout << "Time cost for "
            << kArray_size << " calls to draw "
            << 1 << " number: " << std::flush;
    wall_time_begin = util::get_wall_time();
    for (int ii = 0; ii < kArray_size; ii++) {
        util::set_random_draw_from_queue<int>(discrete_random_array + ii, 1, queue);
        if (ii > 0) {
            if(discrete_random_array[ii] != discrete_random_array[ii-1]) {
                result_all_identical_values = false;
            }
        }
    }
    wall_time_end = util::get_wall_time();
    std::cout << wall_time_end - wall_time_begin << " seconds.\n";

    REQUIRE(!result_all_identical_values);

    std::cout << "Repeating for " << kNum_tests_randomness << " times to check randomness... " << std::flush;
    for (int tt = 0; tt < kNum_tests_randomness; tt++ ){
        util::set_random_draw_from_queue<int>(discrete_random_array, kArray_size, queue);

        result_mean = 0.0;
        for (int ii = 0; ii < kArray_size; ii++) {
            result_mean += float(discrete_random_array[ii])/float(kArray_size);
        }

        REQUIRE(std::abs(result_mean - target_mean) < tolerence);

        if (tt == 0) {
            result_mean_prev = result_mean;
        } else {
            REQUIRE(result_mean != result_mean_prev);
            result_mean_prev = result_mean;
        }
    }
    std::cout << "done.\n";

    result_mean = 0.0;
    result_mean_prev = 0.0;


    // Performance
    std::cout << "Time cost for "
                << kNum_tests << " draws of "
                << kArray_size << " numbers: " << std::flush;

    wall_time_begin = util::get_wall_time();
    for (int tt = 0; tt < kNum_tests; tt++ ){
        util::set_random_draw_from_queue<int>(discrete_random_array, kArray_size, queue);
    }
    wall_time_end = util::get_wall_time();

    std::cout << wall_time_end - wall_time_begin << " seconds.\n";


}


TEST_CASE( "21.4: Randomness: binomial random numbers", "[util:rand:binomial]" ) {

    int kNum_tests = 1000;
    const int kNum_stars = 100;

    int kNum_methods = 3;

    const int kArray_size = 10000;
    int result_random_array[kArray_size] = {};
    const int kResult_counter_size = 20;
    int result_counter[kResult_counter_size] = {0};

    const int kNum_events = 100;
    const float kSuccess_fraction = 1.0/14.0;

    bool result_all_values_identical = true;

    double wall_time_begin = 0.0;
    double wall_time_end = 0.0;

    for (int mm = 0; mm < kNum_methods; mm++){
        
        util::config_binomial_method(mm); //0-boost; 1-std

        if (mm == 0) {
            std::cout << "Using boost::random::binomial_distribution:\n";
        }
        if (mm == 1) {
            std::cout << "Using std::binomial_distribution:\n";
        }
        if (mm == 2) {
            std::cout << "Using bernoulli_process method:\n";
            std::cout << "NOTE: reduction of test size from " << kNum_tests << " to ";
            kNum_tests /= 10;
            std::cout << kNum_tests << "\n";
        }
        // Function usage demo

        std::fill(result_counter, result_counter + kResult_counter_size, 0);

        util::set_random_numbers_binomial(result_random_array, kArray_size, kNum_events, kSuccess_fraction);
        
        std::cout << "Result from one function call for " << kArray_size << " random numbers:\n";
        for (int ii = 0; ii < kArray_size; ii++) {
            if (result_random_array[ii] < kResult_counter_size){
                result_counter[result_random_array[ii]]++;
                // result_all_less_than_one = false;
            }
        }
        for (int cc = 0; cc < kResult_counter_size; cc++) {
            std::cout << std::fixed << std::setprecision(1)
                        << cc << ": "
                        << std::string(result_counter[cc]*kNum_stars/kArray_size, '*')
                        << " " << result_counter[cc] << "\n";
        }
        

        std::fill(result_counter, result_counter + kResult_counter_size, 0);
        result_all_values_identical = true;

        // Randomness from multiple calls
        std::cout << "Result from " << kArray_size << " function calls for 1 random number:\n";
        for (int ii = 0; ii < kArray_size; ii++) {
            util::set_random_numbers_binomial(result_random_array + ii, 1, kNum_events, kSuccess_fraction);
            if (ii > 0) {
                if(result_random_array[ii] != result_random_array[ii-1]) {
                    result_all_values_identical = false;
                }
            }
            if (result_random_array[ii] < kResult_counter_size){
                result_counter[result_random_array[ii]]++;
                // result_all_less_than_one = false;
            }
        }
        for (int cc = 0; cc < kResult_counter_size; cc++) {
            std::cout << std::fixed << std::setprecision(1)
                        << cc << ": "
                        << std::string(result_counter[cc]*kNum_stars/kArray_size, '*')
                        << " " << result_counter[cc] << "\n";
        }
        REQUIRE(!result_all_values_identical);

        std::fill(result_counter, result_counter + kResult_counter_size, 0);

        // Performance
        std::cout << "Time cost for "
                    << kNum_tests << " draws of "
                    << kArray_size << " numbers: " << std::flush;

        wall_time_begin = util::get_wall_time();
        for (int tt = 0; tt < kNum_tests; tt++ ){
            util::set_random_numbers_binomial(result_random_array, kArray_size, kNum_events, kSuccess_fraction);
        }
        wall_time_end = util::get_wall_time();

        std::cout << wall_time_end - wall_time_begin << " seconds.\n";


        std::cout << "Time cost for "
                    << kNum_tests*kArray_size << " draws of "
                    << 1 << " number: " << std::flush;

        wall_time_begin = util::get_wall_time();
        for (int tt = 0; tt < kNum_tests*kArray_size; tt++ ){
            util::set_random_numbers_binomial(result_random_array, 1, kNum_events, kSuccess_fraction);
        }
        wall_time_end = util::get_wall_time();

        std::cout << wall_time_end - wall_time_begin << " seconds.\n";




    }
}

TEST_CASE( "21.5: Randomness: binomial random numbers (mean)", "[util:rand:binomial_mean]" ) {
    
    const int kNum_mean_test_steps = 4;

    const int kNum_mean_test_runs = 1000;

    int num_draws_in_one_test = 0;

    int kNum_methods = 3;

    int result_random_number = 0;

    const int kNum_events = 100;
    const float kSuccess_fraction = 1.0/14.0;

    const float kTarge_mean = kNum_events * kSuccess_fraction;
    const float kTolerence = 0.01;

    bool result_all_values_identical = true;

    float result_mean = 0.0;
    float result_mean_avg = 0.0;

    double wall_time_begin = 0.0;
    double wall_time_end = 0.0;


    std::cout << "Theoretical mean (target): " << kTarge_mean;
    std::cout << "; Number of events: " << kNum_events;
    std::cout << "; Success fraction: " << kSuccess_fraction << "\n";

    

    for(int ss = 0; ss < kNum_mean_test_steps; ss++) {

        int num_draws_in_one_test = std::pow(10,(ss+1));

        for (int mm = 0; mm < kNum_methods; mm++){
        
            util::config_binomial_method(mm); //0-boost; 1-std

            std::cout << "Draw size: " << num_draws_in_one_test << " ; ";
            if (mm == 0) {
                std::cout << "    [boost]";
            }
            if (mm == 1) {
                std::cout << "      [std]";
            }
            if (mm == 2) {
                std::cout << "[bernoulli]";
            }
            
            for (int tt_mean = 0; tt_mean < kNum_mean_test_runs; tt_mean++){

                for (int tt = 0; tt < num_draws_in_one_test; tt++) {
                    util::set_random_numbers_binomial(&result_random_number, 1, kNum_events, kSuccess_fraction);
                    result_mean += float(result_random_number)/num_draws_in_one_test;
                }

                result_mean_avg += float(result_mean)/kNum_mean_test_runs;

                result_mean = 0.0;

            }

            std::cout << std::setprecision(3) << std::fixed;

            std::cout << ", result mean (avg. of " << kNum_mean_test_runs << " runs): " << result_mean_avg << ", "
                        << std::abs(result_mean_avg-kTarge_mean)/kTarge_mean*100 << '%';
            std::cout << std::setprecision(0);
            std::cout << "; within " << kTolerence * 100 << '%' << " tolerence?";
            REQUIRE(std::abs(result_mean_avg-kTarge_mean) < kTarge_mean*kTolerence);

            std::cout << " YES\n";

            // result_mean = 0.0;
            result_mean_avg = 0.0;
        }
        std::cout << "\n";
    }

}

TEST_CASE( "21.5.1: Randomness: binomial get_rand_binomial", "[util:rand:get_rand_binomial]" ) {
    
    const int kNum_mean_test_steps = 3;

    const int kNum_mean_test_runs = 1000;

    int num_draws_in_one_test = 0;

    int kNum_methods = 3;

    // int result_random_number = 0;

    const std::vector<int> kNum_events{100, 200, 300};
    const std::vector<float> kSuccess_fraction{1.0/14.0, 1.0/14.0, 1.0/24.0};

    std::cout << kSuccess_fraction.at(0) << "\n";
    std::cout << kSuccess_fraction.at(1) << "\n";
    std::cout << kSuccess_fraction.at(2) << "\n";

    // const float kTarge_mean = kNum_events * kSuccess_fraction;
    const float kTolerence = 0.01;

    bool result_all_values_identical = true;

    float result_mean = 0.0;
    float result_mean_avg = 0.0;

    // double wall_time_begin = 0.0;
    // double wall_time_end = 0.0;

    std::streamsize default_percision = std::cout.precision();

    for (uint ii = 0; ii < kNum_events.size(); ii++) {
        float kTarge_mean = kNum_events.at(ii) * kSuccess_fraction.at(ii);
        std::cout << "Theoretical mean (target): " << kNum_events.at(ii) * kSuccess_fraction.at(ii);
        std::cout << "; Number of events: " << kNum_events.at(ii);
        std::cout << "; Success fraction: " << kSuccess_fraction.at(ii) << std::endl;

        for(int ss = 0; ss < kNum_mean_test_steps; ss++) {

            int num_draws_in_one_test = std::pow(10,(ss+1));

            for (int mm = 0; mm < kNum_methods; mm++){
            
                util::config_binomial_method(mm); //0-boost; 1-std

                std::cout << "Draw size: " << num_draws_in_one_test << " ; ";
                if (mm == 0) {
                    std::cout << "    [boost]";
                }
                if (mm == 1) {
                    std::cout << "      [std]";
                }
                if (mm == 2) {
                    std::cout << "[bernoulli]";
                }
                
                for (int tt_mean = 0; tt_mean < kNum_mean_test_runs; tt_mean++){

                    for (int tt = 0; tt < num_draws_in_one_test; tt++) {
                        // util::set_random_numbers_binomial(&result_random_number, 1, kNum_events, kSuccess_fraction);
                        result_mean += float(util::get_rand_binomial(kNum_events.at(ii), kSuccess_fraction.at(ii)))/num_draws_in_one_test;
                    }

                    result_mean_avg += float(result_mean)/kNum_mean_test_runs;

                    result_mean = 0.0;

                }

                std::cout << std::setprecision(3) << std::fixed;

                std::cout << ", result mean (avg. of " << kNum_mean_test_runs << " runs): " << result_mean_avg << ", "
                            << std::abs(result_mean_avg-kTarge_mean)/kTarge_mean*100 << '%';
                std::cout << std::setprecision(0);
                std::cout << "; within " << kTolerence * 100 << '%' << " tolerence?";
                REQUIRE(std::abs(result_mean_avg-kTarge_mean) < kTarge_mean*kTolerence);

                std::cout << " YES\n";

                std::cout.precision(default_percision); 

                // result_mean = 0.0;
                result_mean_avg = 0.0;
            }
            std::cout << "\n";
        }
        
        
    }



    


}

TEST_CASE( "21.6: Randomness: normal random numbers", "[util:rand:normal]" ) {

    std::vector<std::string> method_names;
    method_names.push_back("Boost::random::normal_distribution");
    method_names.push_back("std::normal_real_distribution");
// #ifdef __NVCC__
//     method_names.push_back("curand::hostAPI");
// #endif

    const int kNum_tests = 10000;
    const int kNum_stars = 200;
    const int kArray_size = 10000;
    float random_array[kArray_size] = {0.0};
    const int kResult_counter_size = 40;
    int result_counter[kResult_counter_size] = {0};

    std::fill(result_counter, result_counter + kResult_counter_size, 0);

    // bool result_all_less_than_one = true;
    bool result_no_consequtive_identical_values = true;

    double wall_time_begin = 0.0;
    double wall_time_end = 0.0;

    const int kMean = 15;
    const int kSd = 5;

    // for (int bb = 0; bb < 2; bb++){

    //     if (bb == 0) {
    //         util::config_uniform_buffer(false);
    //         std::cout << "No buffer:\n";
    //     } else {
    //         util::config_uniform_buffer(true);
    //         std::cout << "\nWith buffer ("
    //                     << util::kRandomness_buffer_length 
    //                     << " x float = "
    //                     << util::kRandomness_buffer_length/1000000 *sizeof(float)
    //                     << " MB ):\n";
    //     }

        for (int mm = 0; mm < static_cast<int>(method_names.size()); mm++){
            
            // illustrative pdf
            util::config_normal_method(mm); //0-boost; 1-std
            util::set_random_numbers_normal(random_array, kArray_size);
            
            std::cout << "Result using method " << method_names.at(mm) << ":\n";
            for (int ii = 0; ii < kArray_size; ii++) {
                int bin = std::min( std::max(0,int(std::round(random_array[ii]*kSd+kMean))), kResult_counter_size-1 );
                result_counter[bin]++;
                // if (uniform_random_array[ii] >=1 ){
                //     result_all_less_than_one = false;
                // }
            }
            for (int cc = 0; cc < kResult_counter_size; cc++) {
                std::cout << std::fixed << std::setprecision(1)
                            << cc << ": "
                            << std::string(result_counter[cc]*kNum_stars/kArray_size, '*')
                            << " " << result_counter[cc] << "\n";
            }
            // REQUIRE(result_all_less_than_one);

            std::fill(result_counter, result_counter + kResult_counter_size, 0);
            // result_all_less_than_one = true;

            // continue;
            // Randomness
            for (int ii = 0; ii < kArray_size; ii++) {
                util::set_random_numbers_normal(random_array + ii, 1);
                if (ii > 0) {
                    if(random_array[ii] == random_array[ii-1]) {
                        result_no_consequtive_identical_values = false;
                    }
                }
            }
            REQUIRE(result_no_consequtive_identical_values);


            std::cout << "Time cost for "
                        << kNum_tests << " draws of "
                        << kArray_size << " numbers: " << std::flush;

            wall_time_begin = util::get_wall_time();
            for (int tt = 0; tt < kNum_tests; tt++ ){
                util::set_random_numbers_normal(random_array, kArray_size);
            }
            wall_time_end = util::get_wall_time();

            std::cout << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";


            std::cout << "Time cost for "
                        << kNum_tests << " draws of "
                        << 1 << " number: " << std::flush;

            wall_time_begin = util::get_wall_time();
            for (int tt = 0; tt < kNum_tests; tt++ ){
                util::set_random_numbers_normal(random_array, 1);
            }
            wall_time_end = util::get_wall_time();

            std::cout << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";


            std::cout << "Time cost for "
                        << 1 << " draws of "
                        << kArray_size << " number: " << std::flush;

            wall_time_begin = util::get_wall_time();
            // for (int tt = 0; tt < kNum_tests; tt++ ){
                util::set_random_numbers_normal(random_array, kArray_size);
            // }
            wall_time_end = util::get_wall_time();

            std::cout << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";

        }

    // } // if_buffer

}
