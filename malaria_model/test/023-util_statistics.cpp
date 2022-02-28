#include "third_party/catch2/catch.hpp"
#include <iostream>

#include <vector>

// #include <algorithm> // std::fill_n
// #include <limits> //std::numeric_limits

#include "util/statistics.h"


TEST_CASE( "23.1: Util / Statictics / from_n_choose_k", "[util:statistics:from_n_choose_k]" ) {
    const int kN_max = 8;
    std::vector<std::vector<int>> result_per_k;
    std::vector<std::vector<int>> result_all;

    for (int nn = 0; nn <= kN_max; nn++) {

        for (int kk = 0; kk <= nn; kk++) {

            result_per_k = util::from_n_choose_k(nn, kk);

            for (auto& each_combination : result_per_k) {
                for (auto& each_item: each_combination) {
                    std::cout << " " << each_item;
                }
                std::cout << "\n";
            }

            int size = result_per_k.size();
            int expected_size = 0;

            if ( (nn == 0) || (kk == 0) ) {
                expected_size = 0;
            } else {
                expected_size = util::factorial(nn) / (util::factorial(kk) * util::factorial(nn-kk));
            }

            std::cout << "From n-" << nn << " choose k-" << kk << ": "
                        <<" expecting " << expected_size << " [ n!/(k!*(n-k)!) ], got " << size;
            REQUIRE(size == expected_size);
            std::cout << " ... PASS\n";
        }
        
    }

    // REQUIRE(--count == static_cast<int>(common::DrugName::Last));

}

TEST_CASE( "23.2: Util / Statictics / sort_by_value_and_return_indices_descend", "[util:statistics:sort_by_value_and_return_indices_descend]" ) {

    std::vector<float> values{1.0, 3.2, 3.4, 2.0, 10.5, 9.6};
    std::vector<size_t> sorted_indices = util::sort_by_value_and_return_indices_descend(values);

    REQUIRE(sorted_indices.size() == values.size());

    REQUIRE(sorted_indices[0] == 4);
    REQUIRE(sorted_indices[1] == 5);
    REQUIRE(sorted_indices[2] == 2);
    REQUIRE(sorted_indices[3] == 1);
    REQUIRE(sorted_indices[4] == 3);
    REQUIRE(sorted_indices[5] == 0);

    std::cout << "In descending order: {";
    for (auto ii: sorted_indices) {
        std::cout << values[ii] << ", ";
    }
    std::cout << "}\n";

    sorted_indices = util::sort_by_value_and_return_indices_ascend(values);

    REQUIRE(sorted_indices.size() == values.size());

    REQUIRE(sorted_indices[0] == 0);
    REQUIRE(sorted_indices[1] == 3);
    REQUIRE(sorted_indices[2] == 1);
    REQUIRE(sorted_indices[3] == 2);
    REQUIRE(sorted_indices[4] == 5);
    REQUIRE(sorted_indices[5] == 4);

    std::cout << "In ascending order: {";
    for (auto ii: sorted_indices) {
        std::cout << values[ii] << ", ";
    }
    std::cout << "}\n";

    std::vector<float> values2(100,20);

    sorted_indices = util::sort_by_value_and_return_indices_ascend(values2);

    std::cout << "In ascending order: {";
    for (auto ii: sorted_indices) {
        std::cout << ii << "-" << values2[ii] << ", ";
    }
    std::cout << "}\n";

}
