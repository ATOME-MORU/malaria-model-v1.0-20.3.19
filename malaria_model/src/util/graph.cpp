#include <iostream>

// #include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/breadth_first_search.hpp>

#include <boost/range/irange.hpp>
#include <boost/pending/indirect_cmp.hpp>

#include <limits> //std::numeric_limits

#include <boost/graph/kruskal_min_spanning_tree.hpp>

#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <fstream>

#include <cassert>

#include <algorithm> // std::max_element

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
namespace ba = boost::accumulators;

#include "util/graph.h"

namespace util {

VillageGraph::VillageGraph(
    const int num_vll,
    const float* const* distances
): vg(num_vll){
    
    this->num_villages = num_vll; //TODO: replace this with boost::num_vertices(this->vg)
    vg_weight_map = boost::get(boost::edge_weight,vg);
    vg_ifmst_map = boost::get(Ifmst(),vg);
    vg_ifotl_map = boost::get(Ifotl(),vg);

    for (int vv=0; vv < this->num_villages; vv++) {
        vg[boost::vertex(vv,vg)].id = vv;
    }

    this->set_distances_undirected(distances);
}

void VillageGraph::set_distances_undirected(const float* const* distances) {
    Edge e;
    bool added;
    for (int vv1=0; vv1 < this->num_villages; vv1++) {
        for (int vv2=vv1+1; vv2 < this->num_villages; vv2++) {
            // if (distances[vv1][vv2] < std::numeric_limits<float>::max()){
                boost::tie(e, added) = boost::add_edge(
                    boost::vertex(vv1, vg),
                    boost::vertex(vv2, vg),
                    // distances[vv1][vv2],
                    vg
                );
                this->vg_weight_map[e] = distances[vv1][vv2];
                this->vg_ifmst_map[e] = false;
                this->vg_ifotl_map[e] = false;
            // }
        }
    }
}

void VillageGraph::update_edge_weight_stats(bool if_mst_only) {
    ba::accumulator_set<double, ba::stats<ba::tag::mean, ba::stats<ba::tag::variance> > > acc;
    
    boost::graph_traits<Graph>::edge_iterator eei, eei_end;
    for (boost::tie(eei, eei_end) = boost::edges(vg); eei != eei_end; ++eei){
        if(if_mst_only){
            if (this->vg_ifmst_map[*eei]){
                acc(this->vg_weight_map[*eei]);
            }
        } else {
            acc(this->vg_weight_map[*eei]);
        }
    }

    this->edge_set_weight_mean = ba::mean(acc);
    this->edge_set_weight_sd = sqrt(ba::variance(acc));
    

    std::cout << " Edge set weight mean:" << this->edge_set_weight_mean << std::endl;
    std::cout << " Edge set weight sd:" << this->edge_set_weight_sd << std::endl;
}
int VillageGraph::mark_edge_weight_outliers_sd(const int sd_count) {
    int total = 0;

    float threshold = this->edge_set_weight_mean + sd_count*this->edge_set_weight_sd;

    std::cout << " Edge set weight threshold:" << threshold << std::endl;

    boost::graph_traits<Graph>::edge_iterator eei, eei_end;
    for (boost::tie(eei, eei_end) = boost::edges(vg); eei != eei_end; ++eei){
        if ( this->vg_weight_map[*eei] > threshold ){
            this->vg_ifotl_map[*eei] = true;
            total++;
        } else {
            this->vg_ifotl_map[*eei] = false;
        }
    }
    std::cout << " Number of outlier edges marked:" << total << std::endl;
    return total;
}
void VillageGraph::remove_marked_outlier_edges(){
    int temp_counter=0;

    std::cout << "num edges: " << boost::num_edges(vg) << std::endl;

    std::vector<Edge> edges_to_remove;
    boost::graph_traits<Graph>::edge_iterator eei, eei_end;
    for (boost::tie(eei, eei_end) = boost::edges(vg); eei != eei_end; ++eei){
        if (this->vg_ifotl_map[*eei]){
            edges_to_remove.push_back(*eei);
        }
    }

    for (std::vector<Edge>::iterator eei = edges_to_remove.begin(); eei != edges_to_remove.end(); ++eei){
        boost::remove_edge(*eei, this->vg);
        temp_counter++;
    }

    std::cout << "num edges: " << boost::num_edges(vg) << std::endl;

    std::cout << "end of remove_marked_outlier_edges: removed " << temp_counter << std::endl;
}

void VillageGraph::mark_graph_with_min_spanning_tree() {
    std::vector<Edge> mst;
    boost::kruskal_minimum_spanning_tree(vg, std::back_inserter(mst));

    // std::vector<Edge> non_mst;
    boost::graph_traits<Graph>::edge_iterator eei, eei_end;
    //TODO: just iterate mst and mark true
    for (boost::tie(eei, eei_end) = boost::edges(vg); eei != eei_end; ++eei){
        if (std::find(mst.begin(), mst.end(), *eei) == mst.end()){
            // non_mst.push_back(*eei);
            this->vg_ifmst_map[*eei] = false;
        } else {
            this->vg_ifmst_map[*eei] = true;
        }
    }
    // for (std::vector<Edge>::iterator eei = mst.begin(); eei != mst.end(); ++eei){
    //     this->vg_ifmst_map[*eei] = true;
    // }
}
void VillageGraph::reduce_graph_to_min_spanning_tree() {
    std::vector<Edge> edges_to_remove;
    boost::graph_traits<Graph>::edge_iterator eei, eei_end;
    for (boost::tie(eei, eei_end) = boost::edges(vg); eei != eei_end; ++eei){
        if (!this->vg_ifmst_map[*eei]){
            edges_to_remove.push_back(*eei);
        }
    }

    for (std::vector<Edge>::iterator eei = edges_to_remove.begin(); eei != edges_to_remove.end(); ++eei){
        boost::remove_edge(*eei, this->vg);
    }
}

std::vector<float> VillageGraph::find_shortest_paths(int from_vll_id) const {

    // std::vector<Graph::vertex_descriptor> p(boost::num_vertices(this->vg));
    std::vector<float> d(boost::num_vertices(this->vg));

    boost::dijkstra_shortest_paths(this->vg, boost::vertex(from_vll_id, this->vg),
        // predecessor_map(boost::make_iterator_property_map(p.begin(), boost::get(boost::vertex_index, this->vg))).
        distance_map(boost::make_iterator_property_map(d.begin(), boost::get(boost::vertex_index, this->vg)))
    );

    // boost::graph_traits<Graph>::vertex_iterator vi, vend;
    // for (boost::tie(vi, vend) = boost::vertices(this->vg); vi != vend; ++vi) {
    //     std::cout << "distance(" << *vi << ") = " << d[*vi]
    //                 << ", ";
    //                 // << "\n";
    //     std::cout << "parent(" << *vi << ") = " << p[*vi] << std::endl;
    // }

    return d;
}

std::string VillageGraph::print_graph_to_file_dot(std::string file_prefix) {
    const std::string file_name = "_village_graph.dot";
    std::ofstream out_file(file_prefix + file_name);
    out_file << "graph A {\n"
            << " rankdir=LR\n"
            << " size=\"3,3\"\n"
            << " ratio=\"filled\"\n"
            << " edge[style=\"bold\"]\n" << " node[shape=\"circle\"]\n";

    boost::graph_traits<Graph>::edge_iterator eei, eei_end;
    for (boost::tie(eei, eei_end) = boost::edges(vg); eei != eei_end; ++eei){
        out_file << boost::source(*eei, this->vg) << " -- " << boost::target(*eei, this->vg);
        if (this->vg_ifotl_map[*eei]) {
            out_file << "[color=\"red\", label=\"" << boost::get(boost::edge_weight, this->vg, *eei)
                        << "\"];\n"; 
        } else if (this->vg_ifmst_map[*eei]) {
            out_file << "[color=\"black\", label=\"" << boost::get(boost::edge_weight, this->vg, *eei)
                        << "\"];\n";
        } else {
            out_file << "[color=\"grey\", label=\"" << boost::get(boost::edge_weight, this->vg, *eei)
                        << "\"];\n";      
        }
    }
    out_file << "}\n";
    out_file.close();

    return file_prefix + file_name;
}

std::string VillageGraph::print_traversal_path_to_file_gnuplot(
                std::string file_prefix,
                float* x,
                float* y) {
        
    const std::string data_file_name = file_prefix
                                        + "_village_locations_traversal.dat."
                                        + std::to_string(this->discover_orders[0]);
    std::ofstream out_file_data(data_file_name);
    out_file_data << "id longitude latitude visit_time\n";

    for (int vv = 0; vv < this->num_villages; vv++){
        out_file_data << vv << " " << x[vv] << " " << y[vv] << " "
                        << this->discover_times[vv] <<"\n";
    }
    out_file_data.close();

    std::string gnuplot_script_traversal = "set title 'Traversal Map'\n";

    gnuplot_script_traversal += "set label 1 at "
                    + std::to_string(x[this->discover_orders[0]]) + ", "
                    + std::to_string(y[this->discover_orders[0]])
                    + " \"{/:Bold Start}\" offset char -1,2 center\n";
    gnuplot_script_traversal += "set label 2 at "
                    + std::to_string(x[this->discover_orders[this->num_villages-1]]) + ","
                    + std::to_string(y[this->discover_orders[this->num_villages-1]])
                    + " \"{/:Bold Finish}\" offset char -1,2 center\n";

    for (int vv = 0; vv < this->num_villages-1; vv++){
        gnuplot_script_traversal += "set arrow from "
                            + std::to_string(x[this->discover_orders[vv]]) + ","
                            + std::to_string(y[this->discover_orders[vv]]) + " to " 
                            + std::to_string(x[this->discover_orders[vv+1]]) + ","
                            + std::to_string(y[this->discover_orders[vv+1]]) + "\n";
    }

    gnuplot_script_traversal += "plot input_data_file using 2:3 \n";
    const std::string script_file_name_traversal = file_prefix
                                        + "_map_traversal.gnuplot."
                                        + std::to_string(this->discover_orders[0]);
    std::ofstream out_file_script(script_file_name_traversal);
    out_file_script << gnuplot_script_traversal;
    out_file_script.close();

    return "gnuplot -p -e \"input_data_file='"
            + data_file_name + "'\" " + script_file_name_traversal;

}

template <typename T_weight>
std::string VillageGraph::print_graph_to_file_gnuplot(
                const std::string file_prefix,
                const std::string file_name,
                const std::string plot_title,
                const std::string group_name_prefix,
                const float* x,
                const float* y,
                const T_weight* weight,
                const std::vector<int> group,
                const std::vector<int>& labelled_villages
    ) const {

    assert(group.size() == boost::num_vertices(this->vg));

    const std::string data_file_name = file_prefix
                                        + "_" + file_name + ".dat";

    std::ofstream out_file_data(data_file_name);
    out_file_data << "id longitude latitude weight group\n";
    int kTeam_column = 5;

    for (int vv = 0; vv < this->num_villages; vv++){
        out_file_data << vv << " " << x[vv] << " " << y[vv] << " " << weight[vv] << " " << group.at(vv) << "\n";
    }

    out_file_data.close();

    T_weight weight_max = *std::max_element(weight, weight+boost::num_vertices(this->vg));

    std::string gnuplot_script_connections = "set title '" + plot_title + "'\n";

    for (uint ll = 0; ll < labelled_villages.size(); ll++) {
            gnuplot_script_connections += "set label " + std::to_string(ll+1) + " at "
                    + std::to_string(x[labelled_villages[ll]]) + ", "
                    + std::to_string(y[labelled_villages[ll]])
                    + " \"{/:Bold v_{" + std::to_string(labelled_villages[ll])
                    + "} }\" offset 0,1 center point pt 2 lc rgb \"blue\"\n";
    }

    boost::graph_traits<Graph>::edge_iterator eei, eei_end;
    for (boost::tie(eei, eei_end) = boost::edges(vg); eei != eei_end; ++eei){

        gnuplot_script_connections += "set arrow from "
                + std::to_string(x[(int)boost::source(*eei, this->vg)]) + ","
                + std::to_string(y[(int)boost::source(*eei, this->vg)]) + " to " 
                + std::to_string(x[(int)boost::target(*eei, this->vg)]) + ","
                + std::to_string(y[(int)boost::target(*eei, this->vg)]);

        if(this->vg_ifotl_map[*eei]) {
            gnuplot_script_connections += " linecolor rgb \"red\"";
        } else if (this->vg_ifmst_map[*eei]) {
            gnuplot_script_connections += " linecolor rgb \"black\"";
        } else {
            gnuplot_script_connections += " linecolor rgb \"grey\"";
        }

        gnuplot_script_connections += " nohead\n";

    }

    gnuplot_script_connections += "set xlabel \"Longitude\"\n";
    gnuplot_script_connections += "set ylabel \"Latitude\"\n";
    gnuplot_script_connections += "set xtics 0.5\n";
    gnuplot_script_connections += "set ytics 0.5\n";
    gnuplot_script_connections += "set size ratio -1 \n";


    gnuplot_script_connections += "unset colorbox\n";

    auto max_group_id = *std::max_element(std::begin(group), std::end(group));

    int num_groups = 0;

    for (int ii = 0; ii <= max_group_id; ii++) {
        if(std::find(group.begin(), group.end(), ii) != group.end()){
            num_groups++;
        }
    }

    gnuplot_script_connections += "set palette maxcolors " + std::to_string(num_groups+1) + "\n";

    gnuplot_script_connections += "set palette defined (";
    gnuplot_script_connections += this->get_palette_string(num_groups+1);
    gnuplot_script_connections += ")\n";

    gnuplot_script_connections += "set cbrange [-0.5:" + std::to_string(num_groups+0.5) + "] \n";


    // plotting in groups solution 1
    // // gnuplot_script_connections += "plot input_data_file using 2:3:(0.3+2*$4/"
    // gnuplot_script_connections += "plot " + data_file_name + " using 2:3:(0.3+2*$4/"
    //                                 + std::to_string(weight_max)
    //                                 + ".0):5 title \"villages\" pt 6 ps var palette\n";


    // plotting in groups solution 2
    bool plot_issued = false;
    for (int tt = 0; tt <= max_group_id; tt++){

        if(std::find(group.begin(), group.end(), tt) != group.end()){
            if (!plot_issued) {
                gnuplot_script_connections += "plot ";
            } else {
                gnuplot_script_connections += ", ";
            }

            gnuplot_script_connections += "\"< awk '{if($" + std::to_string(kTeam_column) + "=="
                                        + std::to_string(tt) + ") print}' "
                                        + data_file_name +" \" using 2:3:(0.3+2*$4/"
                                        + std::to_string(weight_max)
                                        + "):5 title \"" + group_name_prefix + " "
                                        // + ".0) title \"Team "
                                        + std::to_string(tt) + "\" pt 6 ps var palette";
                                        // + std::to_string(tt) + "\" pt 6 ps var";
            plot_issued = true;
        }
    }



    gnuplot_script_connections += "\n";

    const std::string script_file_name_connections = file_prefix
                                        + "_" + file_name + ".gnuplot";
    std::ofstream out_file_script(script_file_name_connections);
    out_file_script << gnuplot_script_connections;
    out_file_script.close();

    // return "gnuplot -p " + script_file_name_connections;
    return "gnuplot -e \" set terminal pdfcairo enhanced size 5in,5in; set output '"
            + file_prefix + "_" + file_name + ".pdf'\" " + script_file_name_connections;
}
template std::string VillageGraph::print_graph_to_file_gnuplot<int>(
                const std::string file_prefix,
                const std::string file_name,
                const std::string plot_title,
                const std::string group_name_prefix,
                const float* x,
                const float* y,
                const int* weight,
                const std::vector<int> group,
                const std::vector<int>& labelled_villages
) const;
template std::string VillageGraph::print_graph_to_file_gnuplot<float>(
                const std::string file_prefix,
                const std::string file_name,
                const std::string plot_title,
                const std::string group_name_prefix,
                const float* x,
                const float* y,
                const float* weight,
                const std::vector<int> group,
                const std::vector<int>& labelled_villages
) const;


std::string VillageGraph::get_palette_string(uint num_colours) const{

    assert(num_colours > 0);

    std::vector<std::string> colour_definitions;

    colour_definitions.push_back(std::to_string(0) + " '#e41a1c'");
    colour_definitions.push_back(std::to_string(1) + " '#377eb8'");
    colour_definitions.push_back(std::to_string(2) + " '#4daf4a'");
    colour_definitions.push_back(std::to_string(3) + " '#984ea3'");
    colour_definitions.push_back(std::to_string(4) + " '#ff7f00'");
    colour_definitions.push_back(std::to_string(5) + " '#a65628'");
    // colour_definitions.push_back(std::to_string(5) + " '#ffff33'");//yellow
    colour_definitions.push_back(std::to_string(6) + " '#f781bf'");
    colour_definitions.push_back(std::to_string(7) + " '#999999'");
    colour_definitions.push_back(std::to_string(8) + " '#ffffb3'");

    std::string colour_string = "";
    if (num_colours<=colour_definitions.size()){

        colour_string = colour_definitions.at(0);

        for(uint cc = 1; cc < num_colours; cc++){
            colour_string += ", " + colour_definitions.at(cc); 
        }

    } else {

        colour_string = "0 '#000fff', 1 '#90ff70', 2 '#ee0000'";
    }

    return colour_string;
}

std::vector<int> VillageGraph::get_bfs(int head_village){

    // typedef boost::graph_traits<Graph>::vertices_size_type VSize;

    // std::vector<VSize> discover_times(boost::num_vertices(vg));
    // this->discover_times.clear();
    // this->discover_orders.clear();

    // this->discover_times.reserve(boost::num_vertices(vg));
    // this->discover_orders.reserve(boost::num_vertices(vg));
    this->discover_times.resize(boost::num_vertices(vg));
    this->discover_orders.resize(boost::num_vertices(vg));

    typedef boost::iterator_property_map<
                std::vector<VSize>::iterator,
                boost::property_map<Graph,boost::vertex_index_t>::const_type
            > DTimePm;
    DTimePm discover_times_pm(this->discover_times.begin(), boost::get(boost::vertex_index, vg));

    VSize time = 0;

    TimedBfsVisitor<DTimePm> vst(discover_times_pm, time);
    boost::breadth_first_search(vg, boost::vertex(head_village, vg), boost::visitor(vst));

    // std::vector<VSize> discover_orders(boost::num_vertices(vg));

    boost::integer_range<int> index_range(0,boost::num_vertices(vg));
    std::copy(index_range.begin(), index_range.end(), this->discover_orders.begin());
    std::sort(
        this->discover_orders.begin(),
        this->discover_orders.end(),
        boost::indirect_cmp<DTimePm, std::less<VSize>>(discover_times_pm)
    );

    std::vector<int> bfs;

    for (int vv = 0; vv < this->num_villages; vv++ ){
        bfs.push_back((int)this->discover_orders[vv]);
        // std::cout << (int)this->discover_orders[vv] << "->";
    }
    // std::cout << "\n";

    return bfs;

}

uint VillageGraph::get_num_villages() const{
    return boost::num_vertices(this->vg);
}

}