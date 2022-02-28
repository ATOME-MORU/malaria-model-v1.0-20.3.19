#ifndef UTIL_GRAPH_H
#define UTIL_GRAPH_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include <string>
#include <vector>

namespace util {


class VillageGraph {

    struct Node_Village {
        int id;
        // std::string name;
    };

    struct Ifmst { // If included in min-spanning-tree
        typedef boost::edge_property_tag kind;
    };

    struct Ifotl { // If is outlier
        typedef boost::edge_property_tag kind;
    };

    /*
    adjacency_list<OutEdgeContainer, VertexContainer, Directedness,
                   VertexProperties, EdgeProperties,
                   GraphProperties, EdgeList>
    */
    typedef boost::adjacency_list<
                    boost::vecS, boost::vecS, boost::undirectedS,
                    Node_Village,
                    boost::property<boost::edge_weight_t, float,
                    boost::property<Ifmst, bool,
                    boost::property<Ifotl, bool> > >,
                    boost::no_property
                > Graph;
    typedef boost::graph_traits<Graph>::edge_descriptor Edge;


    template <typename TimeMap> class TimedBfsVisitor : public boost::default_bfs_visitor{

        typedef typename boost::property_traits <TimeMap>::value_type T;

        public:

        TimeMap time_map;
        T & ti;

        TimedBfsVisitor(TimeMap tm, T & t) : time_map(tm), ti(t) {}

        // void initialize_vertex(const Graph::vertex_descriptor &s, const Graph &g) const {
        //   std::cout << "Initialize: " << g[s].name << std::endl;
        // }
        void discover_vertex(const Graph::vertex_descriptor &s, const Graph &g) const {
            
            // std::cout << "Discover: " << g[s].id << std::endl;
            (void)g;

            boost::put(time_map, s, ti++);
        }
        // void examine_vertex(const Graph::vertex_descriptor &s, const Graph &g) const {
        //   std::cout << "Examine vertex: " << g[s].name << std::endl;
        // }
        // void examine_edge(const Graph::edge_descriptor &e, const Graph &g) const {
        //   std::cout << "Examine edge: " << g[e] << std::endl;
        // }
        // void tree_edge(const Graph::edge_descriptor &e, const Graph &g) const {
        //   std::cout << "Tree edge: " << g[e] << std::endl;
        // }
        // void non_tree_edge(const Graph::edge_descriptor &e, const Graph &g) const {
        //   std::cout << "Non-Tree edge: " << g[e] << std::endl;
        // }
        // void gray_target(const Graph::edge_descriptor &e, const Graph &g) const {
        //   std::cout << "Gray target: " << g[e] << std::endl;
        // }
        // void black_target(const Graph::edge_descriptor &e, const Graph &g) const {
        //   std::cout << "Black target: " << g[e] << std::endl;
        // }
        // void finish_vertex(const Graph::vertex_descriptor &s, const Graph &g) const {
        //   std::cout << "Finish vertex: " << g[s].name << std::endl;
        // }

    }; //TimedBfsVisitor

    Graph vg;

    int num_villages; //TODO: replace this with boost::num_vertices(this->vg)
    boost::property_map<Graph, boost::edge_weight_t>::type vg_weight_map;
    boost::property_map<Graph, Ifmst>::type vg_ifmst_map;
    boost::property_map<Graph, Ifotl>::type vg_ifotl_map;

    typedef boost::graph_traits<Graph>::vertices_size_type VSize;
    std::vector<VSize> discover_times;  // t_8, t_1, t_10, t_3, t_0, ...
    std::vector<VSize> discover_orders; // v_5, v_1, v_13, v_4, ...

    float edge_set_weight_mean = 0.0;
    float edge_set_weight_sd = 0.0;

    std::string get_palette_string(uint num_colours) const;

public:

    VillageGraph(const int num_vll, const float* const* distances);
    
    void set_distances_undirected(const float* const* distances);

    void update_edge_weight_stats(bool = false);

    int mark_edge_weight_outliers_sd(const int = 3);
    void remove_marked_outlier_edges();

    void mark_graph_with_min_spanning_tree();
    void reduce_graph_to_min_spanning_tree();

    std::vector<float> find_shortest_paths(int from_vll_id) const;

    std::string print_graph_to_file_dot(std::string file_prefix);
    std::string print_traversal_path_to_file_gnuplot(
                std::string file_prefix,
                float* x,
                float* y);

    template <typename T_weight>
    std::string print_graph_to_file_gnuplot(
                const std::string file_prefix,
                const std::string file_name,
                const std::string plot_title,
                const std::string group_name_prefix,
                const float* x,
                const float* y,
                const T_weight* weight,
                const std::vector<int> group,
                const std::vector<int>& labelled_villages
    ) const;

    std::vector<int> get_bfs(int head_village);
    uint get_num_villages() const;

};

}


#endif