#include "third_party/catch2/catch.hpp"
#include <iostream>

#include "common/common.h"
#include "util/util.h"

TEST_CASE( "70.1: Common definitions: comparison of StageName enums", "[common::stage_enum_comparison]" ) {
    int count = 0;
    common::StageName stg_last = common::StageName::First;
    std::cout << "Comparison of Parasite Types:\n";
    for(auto stg: util::Enum_iterator<common::StageName>()) {
        if (count > 0){
            std::cout << stg_last << " < " << stg << "\n";
            REQUIRE(stg_last < stg);
        }
        stg_last = stg;
        count++;
    }
    count = 0;
}

// TEST_CASE( "70.2: Common drug kill rate outputs", "[common::kPds_kill_probability]" ) {

//     // for(auto pp : util::Enum_iterator<common::ParasiteType>()) {
        
//     //     if (pp == common::ParasiteType::kNoParasite) continue;
        
//     //     for(auto dd : util::Enum_iterator<common::DrugName>()) {
            
//     //         if (dd == common::DrugName::kNoDrug) continue;
            
//     //         for(auto ss: util::Enum_iterator<common::StageName>()) {

//     //             if (ss == common::StageName::kNotInSystem) continue;

//     //             std::cout << "kPds["
//     //                         << pp << "]["
//     //                         << dd << "]["
//     //                         << ss << "]="
//     //                         << common::kPds_kill_probability[int(pp)][int(dd)][int(ss)]
//     //                         << "\n";
//     //         }
//     //     }
//     // }

//     //Bsys.print_drug_parameters

// }