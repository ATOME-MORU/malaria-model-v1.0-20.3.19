#include <algorithm> // std::sort
#include <numeric>   // std::iota
#include <iostream>
#include <string>
#include <vector>
#include <cassert>


#include "util/statistics.h"

namespace util {

std::vector<std::vector<int>> from_n_choose_k(const int n, const int k) {

    assert( n >= 0);
    assert( k >= 0);
    assert( k <= n);

    std::vector<std::vector<int>> all_combinations;
    if ( (n == 0) || (k == 0) ) {
        return all_combinations;
    }

    // std::string bitmask(k, 1);
    // bitmask.resize(n, 0);

    std::vector<bool> bitmask(n);
    std::fill(bitmask.end() - k, bitmask.end(), true);


    do {

        std::vector<int> one_combination;
        for (int ii = 0; ii < n; ++ii) {
            if (bitmask[ii]) {

                one_combination.push_back(ii);
                // std::cout << " " << ii;   
            }
        }

        all_combinations.push_back(one_combination);
        // std::cout << std::endl;


    } while (std::next_permutation(bitmask.begin(), bitmask.end()));

    return all_combinations;
}

int from_n_choose_k_get_num_combinations(const int n, const int k) {
    return factorial(n) / (factorial(k) * factorial(n-k));
}

int from_n_choose_up_to_k_get_num_combinations(const int n, const int k) {
    if (k == 0) {
        return 1;
    } else {
        return from_n_choose_k_get_num_combinations(n, k)
                + from_n_choose_up_to_k_get_num_combinations(n, k-1);
    }
}

int factorial(int x) {
    assert(x >= 0);
    if (x == 0) {return 1;}
    return (x == 1 ? x : x * factorial(x - 1));
}


template <typename T>
std::vector<size_t> sort_by_value_and_return_indices_descend(const std::vector<T>& values) {

  // initialize original index locations
  std::vector<size_t> idx(values.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in values
  std::sort(idx.begin(), idx.end(),
       [&values](size_t ii_1, size_t ii_2) {return values[ii_1] > values[ii_2];});

  return idx;
}
template std::vector<size_t> sort_by_value_and_return_indices_descend<int>(const std::vector<int>& values);
template std::vector<size_t> sort_by_value_and_return_indices_descend<float>(const std::vector<float>& values);


template <typename T>
std::vector<size_t> sort_by_value_and_return_indices_ascend(const std::vector<T>& values) {

  // initialize original index locations
  std::vector<size_t> idx(values.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in values
  // std::sort(idx.begin(), idx.end(),
  //      [&values](size_t ii_1, size_t ii_2) {return values[ii_1] < values[ii_2];});
  std::stable_sort(idx.begin(), idx.end(),
       [&values](size_t ii_1, size_t ii_2) {return values[ii_1] < values[ii_2];});

  return idx;
}
template std::vector<size_t> sort_by_value_and_return_indices_ascend<int>(const std::vector<int>& values);
template std::vector<size_t> sort_by_value_and_return_indices_ascend<float>(const std::vector<float>& values);


}
