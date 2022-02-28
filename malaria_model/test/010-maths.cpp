#include <iostream>

#include <fstream>
#include <string>

#include <cmath>//pow()

#include "third_party/catch2/catch.hpp"
#include "common/common.h"
#include "common_test.h"

TEST_CASE( "10.1: Sigmoid Function (Fractional Killing)", "[maths:sigmoid]" ) {
    const float kEC_50 = 60.0;
    const float kGamma = 10.0;

    REQUIRE(common::fractional_killing_function(0,kGamma,kEC_50)==0);
    REQUIRE(common::fractional_killing_function(1,kGamma,kEC_50)<0.001);
    REQUIRE(common::fractional_killing_function(kEC_50,kGamma,kEC_50)==50);
    REQUIRE(common::fractional_killing_function(kEC_50*2,kGamma,kEC_50)-100>-0.1);

    std::ofstream output_file;
    output_file.open(kTest_output_dir + "test_sigmoid.dat");
    for (float cc = 0.0; cc<=kEC_50*2; cc++){
        output_file << cc << " " << common::fractional_killing_function(cc,kGamma,kEC_50) << std::endl;
    }

    std::string gnuplot_command = "gnuplot -p -e \" set title 'Fractional Killing Function';"
           " set xlabel 'Concentration (ng/ml)';"
           // " set term wxt enhanced dashed 'arial,16';"
           " set key right center;"
           " set key box;"
           " set ylabel 'Killing Effect(%)';"
           " set label 'EC_{50}' at ";

    gnuplot_command.append(std::to_string(kEC_50+5));
    gnuplot_command.append(",5 center;");
    gnuplot_command.append("set arrow from ");
    gnuplot_command.append(std::to_string(kEC_50));
    gnuplot_command.append(",0 to ");
    gnuplot_command.append(std::to_string(kEC_50));
    gnuplot_command.append(",");
    gnuplot_command.append(std::to_string(common::fractional_killing_function(kEC_50,kGamma,kEC_50)));
    gnuplot_command.append("nohead dt 2;");
    gnuplot_command.append(" set label 'E_{max}' at -10,100;");
    gnuplot_command.append(" set xtics 10;");
    gnuplot_command.append(" set grid ytics lc rgb '#bbbbbb' lw 1 lt 0;");
    gnuplot_command.append(" set grid xtics lc rgb '#bbbbbb' lw 1 lt 0;");
    gnuplot_command.append(" plot '" + kTest_output_dir + "test_sigmoid.dat'\"");

    output_file.flush();
    output_file.close();
    if(system(gnuplot_command.c_str())){};
}

TEST_CASE( "10.2: Bernstein Function (Parasite Natural Development)", "[maths:bernstein]" ) {
    const uint64_t kClinical_parasite_count_min = 20*pow(10,6);

    const uint16_t kImmunity_level_1_clearence_day = 60;
    const uint16_t kImmunity_level_1_max_parasite_count_day = 10;

    const int kNum_tested_levels = 6;

    common::Bernstein_curve bc;
    bc.calibrate(kImmunity_level_1_clearence_day,
                kImmunity_level_1_clearence_day/kImmunity_level_1_max_parasite_count_day,
                kClinical_parasite_count_min);

    std::cout << "bc.i: " << bc.i
              << ", bc.n: " << bc.n
              << ", bc.nCi: " << bc.nCi
              << ", bc.y_scale: " << bc.y_scale
              << "\n";

    std::ofstream output_file;
    output_file.open(kTest_output_dir + "test_bernstein_curve.dat");

    // write column titles
    output_file << "time_step";
    for (int ll = 0; ll<kNum_tested_levels; ll++){
        output_file << " immunity-level_{"<< ll <<"}";
    }
    output_file << "\n";


    for (uint16_t tt = 0; tt < kImmunity_level_1_clearence_day; tt++){
        output_file << tt;
        for (int ll = 0; ll<kNum_tested_levels; ll++){
            output_file << " " << bc.get_at_level(tt, ll);
        }
        output_file << "\n";
    }

    std::string gnuplot_command = "gnuplot -p -e \" set title 'Parasite Natural Development (Bernstein Curve) ';";
    gnuplot_command.append("set xlabel 'Days After Infection';");
    gnuplot_command.append("set ylabel 'Parasite Count';");

    gnuplot_command.append("set key right center;");
    gnuplot_command.append("set key box;");

    gnuplot_command.append(" set label 'Clinical (above) / Asymptomatic (Below)' at ");
    gnuplot_command.append(std::to_string(kImmunity_level_1_clearence_day/2));
    gnuplot_command.append(",");
    gnuplot_command.append(std::to_string(kClinical_parasite_count_min));
    gnuplot_command.append(" left; set arrow from 0,");
    gnuplot_command.append(std::to_string(kClinical_parasite_count_min));
    gnuplot_command.append(" to ");
    gnuplot_command.append(std::to_string(kImmunity_level_1_clearence_day));
    gnuplot_command.append(",");
    gnuplot_command.append(std::to_string(kClinical_parasite_count_min));
    // gnuplot_command.append(std::to_string(common::fractional_killing_function(kEC_50,kGamma,kEC_50)));
    gnuplot_command.append(" nohead dt 2;");

    // gnuplot_command.append(" set yrange [0.1:*]; set log y 10;");
    gnuplot_command.append(" set yrange [0:*];");
    gnuplot_command.append(" set xrange[0:");
    gnuplot_command.append(std::to_string(kImmunity_level_1_clearence_day));
    gnuplot_command.append("];");

    gnuplot_command.append(" plot for [col=2:");
    gnuplot_command.append(std::to_string(kNum_tested_levels+1));
    gnuplot_command.append("] '" + kTest_output_dir + "test_bernstein_curve.dat' using 1:col ps .5 with linespoints title columnheader\"");
    // gnuplot_command.append("plot 'test_bernstein_curve.dat' using 1:2 ps .5 with linespoints title columnheader\"");
    
    std::cout << gnuplot_command << std::endl;

    output_file.flush();
    output_file.close();
    if(system(gnuplot_command.c_str())){};

}


// TEST_CASE( "1.2: Parasite Decay Phi Function ", "[maths:3]" ) {

//     // REQUIRE(common::parasite_decay_phi(0,common::ParasiteType::kNoParasite)==common::kParasite_decay_Phi0);
//     // REQUIRE(common::parasite_decay_phi(1,common::ParasiteType::kNoParasite)==common::kParasite_decay_Phi0);
//     // REQUIRE(common::parasite_decay_phi(10,common::ParasiteType::kNoParasite)==common::kParasite_decay_Phi0);
//     // REQUIRE(common::parasite_decay_phi(100,common::ParasiteType::kNoParasite)==common::kParasite_decay_Phi0);

//     const int kNum_tested_types = 3;
//     const common::ParasiteType types[]={
//         // common::ParasiteType::kNoParasite,
//         common::ParasiteType::kRa,
//         common::ParasiteType::kRb,
//         common::ParasiteType::kR0
//     };
//     for (int tt = 0; tt < kNum_tested_types; tt++){
//         REQUIRE(common::parasite_decay_phi(0,types[tt])==common::kParasite_decay_Phi0);
//         REQUIRE(common::parasite_decay_phi(1,types[tt])==common::kParasite_immunity_days_max[types[tt]]+common::kParasite_decay_Phi0);
//         if( types[tt] != common::ParasiteType::kNoParasite){
//             REQUIRE(common::parasite_decay_phi(5,types[tt])<common::kParasite_immunity_days_max[types[tt]]+common::kParasite_decay_Phi0);
//         }
//     }


//     const uint16_t kDays_max = 200;

//     std::ofstream output_file;
//     output_file.open("outputs/test/test_decay_phi.dat");
//     // write column headers
//     output_file << "time_step ";
//     for (int tt = 0; tt < kNum_tested_types; tt++){
//         output_file << "type_{" << tt << "} ";
//     }
//     output_file << "\n";

//     // write data
//     for (uint16_t dd = 0; dd <= kDays_max; dd++){
//         output_file << dd << " ";
//         for (int tt = 0; tt < kNum_tested_types; tt++){
//             output_file << common::parasite_decay_phi(dd,types[tt]) << " ";
//         }
//         output_file << "\n";
//     }

//     std::string gnuplot_command = "gnuplot -p -e \" set title 'Parasite Decay Phi';";
//     gnuplot_command.append(" set xrange[0:");
//     gnuplot_command.append(std::to_string(kDays_max));
//     gnuplot_command.append("];");

//     gnuplot_command.append(" plot for [col=2:");
//     gnuplot_command.append(std::to_string(kNum_tested_types+1));
//     gnuplot_command.append("] 'outputs/test/test_decay_phi.dat' using 1:col title columnheader\"");

//     output_file.flush();
//     output_file.close();
//     // system(gnuplot_command.c_str());
// }

// TEST_CASE( "1.3: Parasite Natural Development (Decay/Replication) Function ", "[maths:4]" ) {
//     const float kTime_max = 100;
//     const int kNum_parasites = 20;

//     uint16_t days[kNum_parasites];
//     common::ParasiteType types[kNum_parasites];
//     float phis[kNum_parasites];
//     uint16_t ages[kNum_parasites];
//     uint64_t sizes[kNum_parasites];
    
//     for(int pp=0; pp<kNum_parasites; pp++){
//         days[pp] = pp * 5;
//         types[pp] = common::ParasiteType::kRa;
//         phis[pp] = common::parasite_decay_phi(days[pp],types[pp]);
//         ages[pp] = 0;
//         sizes[pp] = 100;
//     }


//     uint64_t step_temp[] = {0,0};

//     std::ofstream output_file;
//     output_file.open("test_decay_and_replication.dat");
//     // write column titles
//     output_file << "time_step ";
//     for(int pp = 0; pp < kNum_parasites; pp++){
//         output_file << "parasite_{" << pp << "}-idays-" << days[pp] << "-type-" << (int)types[pp] << " ";
//     }
//     output_file << "\n";

//     // wirte time_step 0
//     output_file << "0 ";
//     for(int pp = 0; pp < kNum_parasites; pp++){
//         output_file << sizes[pp] << " ";
//     }
//     output_file << "\n";

//     // write data
//     for (int tt = 1; tt <= kTime_max; tt++){
//         output_file << tt << " ";
//         for(int pp = 0; pp < kNum_parasites; pp++){
//             common::parasite_natural_development(
//                 days, phis, types, ages, sizes, pp, step_temp
//             );

//             output_file << sizes[pp] << " ";

//             days[pp] =  days[pp]>0 ? days[pp]+1 : 0;
//             // days[pp]++;
//             ages[pp]++;

//         }
//         output_file << "\n";
//     }
//     std::string gnuplot_command = "gnuplot -p -e \" set title 'Parasite Natural Development (Decay/Replication) Function';";
//     gnuplot_command.append(" set yrange [0.1:*]; set log y 10;");
//     // gnuplot_command.append(" set yrange [0:*];");
//     gnuplot_command.append(" set xrange[0:100];");

//     gnuplot_command.append(" plot for [col=2:");
//     gnuplot_command.append(std::to_string(kNum_parasites+1));
//     gnuplot_command.append("] 'test_decay_and_replication.dat' using 1:col ps .5 with linespoints title columnheader\"");
    
//     std::cout << gnuplot_command << std::endl;

//     output_file.flush();
//     output_file.close();
//     // system(gnuplot_command.c_str());

// }


// TEST_CASE( "1.6: Parasite Decay/Replication Function ", "[maths:7]" ) {
//     const int kNum_tests = 5;

//     common::ParasiteType types[kNum_tests];
//     uint64_t sizes[kNum_tests];

//     uint64_t immunity_sizes[kNum_tests];
//     float immunity_rates[kNum_tests];


//     uint16_t ages[kNum_tests];
//     for (int tt = 0; tt < kNum_tests; tt++){
//         sizes[tt] = 100000;
//         types[tt] = common::ParasiteType::kRa;
//         ages[tt] = 0;
//         // immunity_sizes[tt] = common::immunity_kill_size_day_zero[(int)types[tt]];
//         immunity_sizes[tt] = 80000;
//         immunity_rates[tt] = 1.0 + (tt+1)*0.2;
//     }
//     // types[1] = common::ParasiteType::kRb;
//     // types[2] = common::ParasiteType::kR0;


//     std::ofstream output_file;
//     output_file.open("test_dandr.dat");


//     output_file << "time_step";
//     for (int tt = 0; tt < kNum_tests; tt++){
//         output_file << " parasite-type_{" << (int) types[tt]
//                     << "}-r="
//                     << common::kParasite_replication_rate
//                     [(int)types[tt]]
//                     [common::kParasite_replication_cycle_length-1]
//                     << "-i="
//                     << immunity_sizes[tt]
//                     << "-ir="
//                     << immunity_rates[tt];

//     }
//     output_file << "\n";

//     // wirte time_step 0
//     output_file << "-1 ";
//     for (int ii = 0; ii < kNum_tests; ii++) {
//         output_file << sizes[ii] << " ";
//     }
//     output_file << "\n";

//     for (int tt = 0; tt < 100; tt++) {
//         output_file << tt;

//         for (int ii = 0; ii < kNum_tests; ii++){
//             common::parasite_decay_and_replication(sizes,types,ages,immunity_sizes,immunity_rates, ii);
//             output_file << " " << sizes[ii];
//             ages[ii]++;
//         }
//         output_file << "\n";

//     }

//     std::string gnuplot_command = "gnuplot -p -e \" set title 'Parasite Natural Development (Flat) ';";
//     gnuplot_command.append(" set yrange [0.1:*]; set log y 10;");
//     // gnuplot_command.append(" set yrange [0:*];");
//     gnuplot_command.append(" plot for [col=2:");
//     gnuplot_command.append(std::to_string(kNum_tests+1));

//     gnuplot_command.append("] 'test_dandr.dat' using 1:col ps .5 with linespoints title columnheader\"");

//     std::cout << gnuplot_command << std::endl;

//     output_file.flush();
//     output_file.close();
//     system(gnuplot_command.c_str());

// }

// TEST_CASE( "1.5: Exponent Function ", "[maths:6]" ) {

//     const int kNum_tested_parasites = 4;


//     std::ofstream output_file;
//     output_file.open("test_decay_exponent.dat");

//     std::string gnuplot_command = "gnuplot -p -e \" set title 'Parasite Natural Development (Bernstein Curve) ';";



//     gnuplot_command.append("] 'test_bernstein_curve.dat' using 1:col ps .5 with linespoints title columnheader\"");

//     std::cout << gnuplot_command << std::endl;

//     output_file.flush();
//     output_file.close();
//     system(gnuplot_command.c_str());

// }
