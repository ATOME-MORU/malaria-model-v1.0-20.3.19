#include "third_party/catch2/catch.hpp"
#include <iostream>

#include <algorithm> // std::fill_n
#include <limits> //std::numeric_limits

#include "util/util.h"
#include "common/common.h"

#include "common_test.h"

TEST_CASE( "20.1: Utility functions: iterate and print name of scoped enums", "[util:enum_iterator]" ) {
    int count = 0;
    std::cout << "Available Parasite Types:\n";
    for(auto pt: util::Enum_iterator<common::ParasiteType>()) {
        std::cout << static_cast<int>(pt) << ": " << pt << std::endl;
        count++;
    }
    REQUIRE(--count == static_cast<int>(common::ParasiteType::Last));

    count = 0;
    std::cout << "Available Stages:\n";
    for(auto pt: util::Enum_iterator<common::StageName>()) {
        std::cout << static_cast<int>(pt) << ": " << pt << std::endl;
        count++;
    }
    REQUIRE(--count == static_cast<int>(common::StageName::Last));

    count = 0;
    std::cout << "Available Drugs:\n";
    for(auto pt: util::Enum_iterator<common::DrugName>()) {
        std::cout << static_cast<int>(pt) << ": " << pt << std::endl;
        count++;
    }
    REQUIRE(--count == static_cast<int>(common::DrugName::Last));

}
