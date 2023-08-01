#pragma once

#include <omp.h>
#include <stdint.h>

#include <limits>
#include <string>
#include <vector>

namespace olsp {
struct Edge {
    int m_target;  // corresponds to index of node in graph datastructure
    int m_cost;
    int m_child_1;  // if child edges aren't -1 that means this is a shortcut edge (see Contraction Hierarchies)
    int m_child_2;
    Edge(int target, int cost) : m_target(target), m_cost(cost), m_child_1(-1), m_child_2(-1) {}
    Edge(int target, int cost, int child_1, int child_2)
        : m_target(target), m_cost(cost), m_child_1(child_1), m_child_2(child_2) {}

    bool isShortcut() { return (m_child_1 != -1) && (m_child_2 != -1); }
};

/**
 * @brief Used for calculating the CH
 *
 */
struct ContractionData {
    std::vector<bool> m_outgoing;
    std::vector<int> m_reset_outgoing;
    std::vector<bool> m_visited;
    std::vector<int> m_reset_visited;
    std::vector<int> m_distances;
    std::vector<int> m_reset_distances;
    std::vector<int> m_num_contracted_neighbours;
    std::vector<std::pair<int, Edge>> m_shortcuts_fwd;
    std::vector<std::pair<int, Edge>> m_shortcuts_bwd;

    ContractionData(int num_nodes)
        : m_outgoing(num_nodes, false),
          m_visited(num_nodes, false),
          m_distances(num_nodes, std::numeric_limits<int>::max()) {}
    ContractionData() = default;
};

/**
 * @brief Mainly used for calculating the lower bound.
 * Also used for verification and reduceESC improvement (which doesn't produces valid ESC)
 *
 */
struct LowerBoundData {
    int m_start_node;
    int m_threshold;
    std::vector<bool> m_marked;
    std::vector<int> m_distances;
    std::vector<int> m_previous_node;
    std::vector<int> m_reset_previous_node;

    LowerBoundData(int num_nodes)
        : m_marked(num_nodes, false),
          m_distances(num_nodes, std::numeric_limits<int>::max()),
          m_previous_node(num_nodes, -1) {}
    LowerBoundData() = default;
};

/**
 * @brief Used for distances and path queries.
 * Pre defining all data allows for many queries at a time without allocating/deallocating a huge amount of memory
 *
 */
struct QueryData {
    int m_start;
    int m_end;
    int m_meeting_node;
    int m_distance;
    bool m_path_needed;           // when the path is not needed the path related variables don't need to be filled`
    std::vector<int> m_fwd_prev;  // fwd = forward
    std::vector<int> m_bwd_prev;  // bwd = backward
    std::vector<int> m_shortest_path;

    std::vector<bool> visited_fwd;
    std::vector<bool> visited_bwd;
    std::vector<int> distances_fwd;
    std::vector<int> distances_bwd;
    std::vector<int> reset_fwd;
    std::vector<int> reset_bwd;

    QueryData(int start, int end, int num_nodes, bool path_needed)
        : m_start(start),
          m_end(end),
          m_meeting_node(-1),
          m_distance(std::numeric_limits<int>::max()),
          m_path_needed(path_needed),
          m_fwd_prev(0),
          m_bwd_prev(0),
          m_shortest_path(0),
          visited_fwd(num_nodes, false),
          visited_bwd(num_nodes, false),
          distances_fwd(num_nodes, std::numeric_limits<int>::max()),
          distances_bwd(num_nodes, std::numeric_limits<int>::max()) {}
    QueryData() = default;
};

/**
 * @brief Defines whether the used graph file already has pre-calculated contraction hierarchies.
 *
 */
enum ReadMode { NORMAL = 0, CONTRACTION_HIERARCHY = 1 };

/**
 * @brief Defines which metric is uses for edge distances.
 * 1 travel time equals 1/100 s
 *
 */
enum DistanceMode { TRAVEL_TIME = 0, DISTANCE_METERS = 1 };

/**
 * @brief Defines which heuristic is used for calculating contraction hierarchies
 * IN_OUT and EDGE_DIFFERENCE are baseline heuristics.
 * WEIGHTED_COST is the preferred heuristic when not using IS.
 * MICROSOFT is preferred when using IS.
 */
enum Heuristic { IN_OUT = 0, EDGE_DIFFERENCE = 1, WEIGHTED_COST = 2, MICROSOFT = 3 };

/**
 * @brief Stores graph data structure and allows for various operations on graph data structure.
 *
 */
class Graph {
   public:
    /**
     * @brief Construct a new Graph object independent sets (IS)
     *
     * @param path          specifies location of graph
     * @param read_mode     specifies if CH is pre-calculated and stored in graph-file
     * @param ch_available  specifies if CH should be used
     * @param prune_graph   specifies if graph should be efficiently stored by using CH characteristics
     * @param ch_heuristic  specifies which heuristic to use
     * @param dist_mode     specifies whether to use travel time or distance in meter metric
     */
    Graph(const std::string& path, ReadMode read_mode, bool ch_available, bool prune_graph, Heuristic ch_heuristic,
          DistanceMode dist_mode = DistanceMode::TRAVEL_TIME);

    // constructor with independent sets (IS)
    /**
     * @brief Construct a new Graph object with independent sets (IS)
     *
     * @param path          specifies location of graph
     * @param read_mode     specifies if CH is pre-calculated and stored in graph-file
     * @param ch_available  specifies if CH should be used
     * @param prune_graph   specifies if graph should be efficiently stored by using CH characteristics
     * @param num_threads   specifies the number of threads to use in createHubLabelsWithIS and createCHwithIS functions
     * @param ch_heuristic  specifies which heuristic to use
     * @param dist_mode     specifies whether to use travel time or distance in meter metric
     */
    Graph(const std::string& path, ReadMode read_mode, bool ch_available, bool prune_graph, int num_threads,
          Heuristic ch_heuristic, DistanceMode dist_mode = DistanceMode::TRAVEL_TIME);

    Graph() = default;
    ~Graph() = default;

    /**
     * @brief Conventional implementation of the dijkstra algorithm
     * Should only be used to compare query times to other methods.
     *
     * @param graph since the function is static the graph vector must be passed
     * @param start node
     * @param end node
     * @return distances between start and end
     */
    static int dijkstraQuery(std::vector<std::vector<Edge>>& graph, int start, int end);

    /**
     * @brief Calculates distance by using bidirectional dijkstra.
     * Is about two times faster when compared to normal dijkstra.
     *
     * @warning Don't use this function with a graph based on contraction hierarchies.
     * It will sometimes produce wrong results depending on the hierarchy.
     * instead use contractionHierachyQuery function
     *
     * @param data
     */
    void bidirectionalDijkstraQuery(QueryData& data);

    /**
     * @brief Extracts path from QueryData.
     * dijstraQuery function with m_path_needed = true needs to be called beforehand.
     *
     * @param data
     */
    void bidirectionalDijkstraGetPath(QueryData& data);

    /**
     * @brief Works similarly to bidirectionalDijkstraQuery.
     *  Additionally uses CH information to reduce query time.
     *
     * @param data set m_path_needed = true to make path extractable
     */
    void contractionHierachyQuery(QueryData& data);

    /**
     * @brief Create a hub labeling when not using IS for CH.
     *
     * @param threshold used for pruned hub labeling
     */
    void createHubLabelsWithoutIS(int threshold = std::numeric_limits<int>::max());

    /**
     * @brief Create a hub labeling when using IS for CH.
     * Implementation uses OpenMP for parallelization with specified number of threads.
     *
     * @param threshold used for pruned hub labeling
     */
    void createHubLabelsWithIS(int threshold = std::numeric_limits<int>::max());

    /**
     * @brief Calculates distance using the created hub labeling.
     * Make sure hub labeling has been created before using this function.
     *
     * @param data
     */
    void hubLabelQuery(QueryData& data);

    double averageLabelSize();

    int maxLabelSize();

    /**
     * @brief Returns the number of label with distance in range of upper and lower.
     *
     * @param upper
     * @param lower
     * @return number of labels in range
     */
    int numHubLabelsInRange(int upper, int lower);

    /**
     * @brief Creates an ESC set based on the algorithm in the bachelor thesis.
     *
     * @param threshold corresponds with range of the electric vehicle
     * @return ESC set
     */
    std::vector<int> createESC(int threshold);

    /**
     * @brief Verifies an ESC based on the algorithm given in the bachelor thesis.
     *
     * @param esc_set which has been calculated prior
     * @param threshold that has been used for esc_set
     * @return true means set is valid (takes a long time even for small graphs and thresholds), false means it is
     * invalid
     */
    bool verifyShortestPathCover(std::vector<int>& esc_set, int threshold);

    /**
     * @brief This corresponds to the second not working improvement in the bachelor thesis.
     *
     * @param esc_set which has been calculated prior
     * @param threshold that has been used for esc_set
     * @return new smaller esc_set
     */
    std::vector<int> reduceESC(std::vector<int>& esc_set, int threshold);

    /**
     * @brief Calculates lower bound based on algorithm from S. Funke et. al. (also described in bachelor thesis).
     *
     * @param esc_set which has been calculated prior
     * @param threshold that has been used for esc_set
     * @return lower bound set (only the size of the set is relevant)
     */
    std::vector<int> lowerBound(std::vector<int>& esc_set, int threshold);

    std::vector<std::vector<Edge>>& getGraph() { return m_graph; }

    int getNumNodes() { return m_num_nodes; }

    void setNumThreads(int num_of_threads) { m_num_threads = num_of_threads; }

    /**
     * @brief Clears all the hub label data from memory.
     * Useful when ESC has been extracted from hub labeling and hub labeling is no longer required.
     *
     */
    void clearHubLabel() {
        std::vector<uint64_t>().swap(m_fwd_indices);
        std::vector<uint64_t>().swap(m_bwd_indices);
        std::vector<std::pair<int, int>>().swap(m_fwd_hub_labels);
        std::vector<std::pair<int, int>>().swap(m_bwd_hub_labels);
    }

   private:
    bool m_ch_available;  // ch/CH = Contraction Hierarchy
    int m_num_nodes;
    bool m_is;          // determines whether IS are used
    int m_num_threads;  // sets threads when using independent sets (IS)

    std::vector<std::vector<Edge>> m_graph;
    std::vector<std::vector<Edge>> m_reverse_graph;
    std::vector<int> m_node_level;

    std::vector<ContractionData> m_contr_data;  // only used to share state when calculating CH

    std::vector<int> m_level_indices_sorted;
    std::vector<int> m_node_indices;
    std::vector<uint64_t> m_fwd_indices;
    std::vector<uint64_t> m_bwd_indices;
    std::vector<std::pair<int, int>> m_fwd_hub_labels;
    std::vector<std::pair<int, int>> m_bwd_hub_labels;

    void readGraph(const std::string& path, ReadMode read_mode, DistanceMode dist_mode);

    void createReverseGraph(bool prune_graph);

    void createReverseGraphCH();

    void createReverseGraphNormal();

    /**
     * @brief Calculates the distance of two given points based on the Haversine Formula
     *
     * @param lat_1
     * @param lon_1
     * @param lat_2
     * @param lon_2
     * @return distance between the points in meter
     */
    int greatCircleDistance(double lat_1, double lon_1, double lat_2, double lon_2);

    /**
     * @brief Calculates hierarchy and adds shortcuts to graph without using independent sets.
     * This means that each node has a distinct level/hierarchy value.
     *
     * @param heuristic
     */
    void createCHwithoutIS(Heuristic heuristic);

    /**
     * @brief Calculates hierarchy and adds shortcuts to graph without using independent sets.
     * This means that many nodes have the same level/hierarchy value.
     *
     * @param heuristic
     */
    void createCHwithIS(Heuristic heuristic);

    int inOutProductHeuristic(std::vector<bool>& contracted, int node);

    int edgeDifferenceHeuristic(std::vector<bool>& contracted, int node);

    int weightedCostHeuristic(std::vector<bool>& contracted, int node);

    int altWeightedCostHeuristic(std::vector<bool>& contracted, int node, std::vector<int>& longest_path_fwd,
                                 std::vector<int>& longest_path_bwd);

    int deletedNeighboursHeuristic(std::vector<bool>& contracted, int node, std::vector<int>& num_deleted_neighbours);

    int microsoftHeuristic(std::vector<bool>& contracted, int node, int cur_level);

    void contractNode(std::vector<bool>& contracted, int contracted_node, int thread_num);

    void contractionDijkstra(int start, int contracted_node, std::vector<bool>& contracted, int num_outgoing,
                             int max_distance, int thread_num);

    void lowerBoundDijkstra(LowerBoundData& lb_data);

    // used for reduceESC function (second not working improvement)
    bool forwardDijkstraSearch(LowerBoundData& lb_data);

    // used for reduceESC function (second not working improvement)
    bool backwardDijkstraSearch(LowerBoundData& lb_data);

    bool pathCoverVerificationDijkstra(LowerBoundData& lb_data);

    int simplifiedHubLabelQuery(std::vector<std::pair<int, int>>& fwd_labels, int node);

    int simplifiedHubLabelQuery(int node, std::vector<std::pair<int, int>>& bwd_labels);
};

}  // namespace olsp
