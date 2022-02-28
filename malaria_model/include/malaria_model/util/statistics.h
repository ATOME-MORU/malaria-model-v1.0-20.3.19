#ifndef UTIL_STATISTICS_H
#define UTIL_STATISTICS_H

namespace util {

std::vector<std::vector<int>> from_n_choose_k(int n, int k);
int from_n_choose_k_get_num_combinations(const int n, const int k);
int from_n_choose_up_to_k_get_num_combinations(const int n, const int k);
int factorial(int x);

template <typename T>
std::vector<size_t> sort_by_value_and_return_indices_descend(const std::vector<T>& values);
template <typename T>
std::vector<size_t> sort_by_value_and_return_indices_ascend(const std::vector<T>& values);

}


#endif