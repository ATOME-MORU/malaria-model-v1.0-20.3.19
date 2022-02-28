#include "third_party/catch2/catch.hpp"
#include <iostream>

#include <algorithm> // std::fill_n
#include <limits> //std::numeric_limits
#include <iomanip>

// #include "util/util.h"
#include "util/graph.h"
#include "common/common.h"

#include "common_test.h"


TEST_CASE( "22.1: Graph-related algorithms", "[util:graph]"){

    // const int num_villages = 6;
    // const float dist[num_villages][num_villages] = {
    //     {0.0, 9.0, 9.0, 9.0, 9.0, 9.0 },
    //     {0.0, 0.0, 2.0, 1.0, 3.0, 3.0 },
    //     {0.0, 0.0, 0.0, 3.0, 7.0, 9.0 },
    //     {0.0, 0.0, 0.0, 0.0, 2.0, 9.0 },
    //     {0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },
    //     {0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
    // };

    // float** distances;

    // distances = new float *[num_villages];
    // for (int vv = 0; vv < num_villages; vv++) {
    //     distances[vv] = new float[num_villages];
    //     std::fill_n(distances[vv], num_villages, std::numeric_limits<float>::infinity());
    // }


    // std::cout << "Village distances:\n";
    // for (int vv1 = 0; vv1 < num_villages; vv1++) {
    //     for (int vv2 = 0; vv2 < num_villages; vv2++) {
    //         distances[vv1][vv2] = dist[vv1][vv2];
    //         std::cout << "[" << std::setw(1) << vv1
    //                   << "][" << std::setw(1) << vv2
    //                   << "]" << std::setw(2) << std::setprecision(1)
    //                   << distances[vv1][vv2] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // util::VillageGraph vg(num_villages, distances);


    // std::string outfile_name = vg.print_graph_to_file(kTest_output_dir + "test_graph_original");
    // std::string xdot_command = "xdot " + outfile_name + " &";
    // system(xdot_command.c_str());

    // std::vector<int> bfs;

    // std::cout << "Before min-span-tree reduction:\n";
    // for (int vv = 0; vv < num_villages; vv++){
    //     std::cout << "vg.find_dfs(" << vv << "):\n";
    //     bfs=vg.get_bfs(vv);

    //     for (const auto& vvi: bfs){std::cout << vvi << " ";} 
    //     std::cout << std::endl;
        
    //     REQUIRE(bfs.front()==vv);
    //     REQUIRE(bfs.size()==num_villages);       
    // }

    // vg.reduce_graph_to_min_spanning_tree();


    // outfile_name = vg.print_graph_to_file(kTest_output_dir + "test_graph_min_span_tree");
    // xdot_command = "xdot " + outfile_name + " &";
    // system(xdot_command.c_str());

    // std::cout << "After min-span-tree reduction:\n";
    // for (int vv = 0; vv < num_villages; vv++){
    //     std::cout << "vg.find_dfs(" << vv << "):\n";
    //     bfs=vg.get_bfs(vv);

    //     for (const auto& vvi: bfs){std::cout << vvi << " ";} 
    //     std::cout << std::endl;
        
    //     REQUIRE(bfs.front()==vv);
    //     REQUIRE(bfs.size()==num_villages);       
    // }


}

TEST_CASE( "22.2: Shortest Path", "[util:graph::shortest_path]"){

    const int num_villages = 6;
    const float dist[num_villages][num_villages] = {
        {0.0, 1.0, 9.0, 9.0, 9.0, 9.0 },
        {0.0, 0.0, 2.0, 1.0, 3.0, 3.0 },
        {0.0, 0.0, 0.0, 3.0, 7.0, 9.0 },
        {0.0, 0.0, 0.0, 0.0, 2.0, 9.0 },
        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
    };

    float** distances;

    distances = new float *[num_villages];
    for (int vv = 0; vv < num_villages; vv++) {
        distances[vv] = new float[num_villages];
        std::fill_n(distances[vv], num_villages, std::numeric_limits<float>::infinity());
    }


    std::cout << "Village distances:\n";
    for (int vv1 = 0; vv1 < num_villages; vv1++) {
        for (int vv2 = 0; vv2 < num_villages; vv2++) {
            distances[vv1][vv2] = dist[vv1][vv2];
            std::cout << "[" << std::setw(1) << vv1
                      << "][" << std::setw(1) << vv2
                      << "]" << std::setw(2) << std::setprecision(1)
                      << distances[vv1][vv2] << " ";
        }
        std::cout << std::endl;
    }

    util::VillageGraph vg(num_villages, distances);

    std::vector<float> shortest_dist(num_villages);
    int start_village = 0;

    shortest_dist = vg.find_shortest_paths(start_village);
    for (int ii = 0; ii < num_villages; ii++){
        std::cout << start_village << "-" << ii << ": " << shortest_dist[ii] << "\n";
        
    }

    start_village = 2;

    shortest_dist = vg.find_shortest_paths(start_village);
    for (int ii = 0; ii < num_villages; ii++){
        std::cout << start_village << "-" << ii << ": " << shortest_dist[ii] << "\n";
        
    }


    // std::string outfile_name = vg.print_graph_to_file(kTest_output_dir + "test_graph_original");
    // std::string xdot_command = "xdot " + outfile_name + " &";
    // system(xdot_command.c_str());

    // std::vector<int> bfs;

    // std::cout << "Before min-span-tree reduction:\n";
    // for (int vv = 0; vv < num_villages; vv++){
    //     std::cout << "vg.find_dfs(" << vv << "):\n";
    //     bfs=vg.get_bfs(vv);

    //     for (const auto& vvi: bfs){std::cout << vvi << " ";} 
    //     std::cout << std::endl;
        
    //     REQUIRE(bfs.front()==vv);
    //     REQUIRE(bfs.size()==num_villages);       
    // }

    // vg.reduce_graph_to_min_spanning_tree();


    // outfile_name = vg.print_graph_to_file(kTest_output_dir + "test_graph_min_span_tree");
    // xdot_command = "xdot " + outfile_name + " &";
    // system(xdot_command.c_str());

    // std::cout << "After min-span-tree reduction:\n";
    // for (int vv = 0; vv < num_villages; vv++){
    //     std::cout << "vg.find_dfs(" << vv << "):\n";
    //     bfs=vg.get_bfs(vv);

    //     for (const auto& vvi: bfs){std::cout << vvi << " ";} 
    //     std::cout << std::endl;
        
    //     REQUIRE(bfs.front()==vv);
    //     REQUIRE(bfs.size()==num_villages);       
    // }


}