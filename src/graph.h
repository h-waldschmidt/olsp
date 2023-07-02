#pragma once

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

struct ContractionData {
    std::vector<bool> m_outgoing;
    std::vector<int> m_reset_outgoing;
    std::vector<bool> m_visited;
    std::vector<int> m_reset_visited;
    std::vector<int> m_distances;
    std::vector<int> m_reset_distances;
    std::vector<std::pair<int, Edge>> m_shortcuts_fwd;
    std::vector<std::pair<int, Edge>> m_shortcuts_bwd;

    std::vector<int> m_num_contracted_neighbours;

    ContractionData(int num_nodes)
        : m_outgoing(num_nodes, false),
          m_visited(num_nodes, false),
          m_distances(num_nodes, std::numeric_limits<int>::max()) {}
    ContractionData() = default;
};

struct AdvancedHubLabelData {
    std::vector<bool> m_visited_fwd;
    std::vector<bool> m_visited_bwd;
    std::vector<int> m_distances_fwd;
    std::vector<int> m_distances_bwd;

    std::vector<int> m_reset_nodes_fwd;
    std::vector<int> m_reset_nodes_bwd;

    AdvancedHubLabelData(int num_nodes)
        : m_visited_fwd(num_nodes, false),
          m_visited_bwd(num_nodes, false),
          m_distances_fwd(num_nodes, std::numeric_limits<int>::max()),
          m_distances_bwd(num_nodes, std::numeric_limits<int>::max()) {}
    AdvancedHubLabelData() = default;
};

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

enum ReadMode { NORMAL = 0, CONTRACTION_HIERARCHY = 1 };

enum DistanceMode { TRAVEL_TIME = 0, DISTANCE_METERS = 1 };

enum Heuristic { IN_OUT = 0, EDGE_DIFFERENCE = 1, WEIGHTED_COST = 2, MICROSOFT = 3 };

class Graph {
   public:
    Graph(const std::string& path, ReadMode read_mode, bool ch_available, bool prune_graph, int num_threads,
          Heuristic ch_heuristic, DistanceMode dist_mode = DistanceMode::TRAVEL_TIME);
    Graph(std::vector<std::vector<Edge>> graph);  // TODO: Adjust constructors for normal mode and ch mode
    ~Graph() = default;

    void readGraph(const std::string& path, ReadMode read_mode, DistanceMode dist_mode);

    void createReverseGraph(bool prune_graph);

    int getNumNodes() { return m_num_nodes; }

    static int dijkstraQuery(std::vector<std::vector<Edge>>& graph, int start, int end);

    // Don't use this function with a graph based on contraction hierachies
    // It will sometimes produce wrong results depending on the hierachiy
    // instead use contractionHierachyQuery function
    void bidirectionalDijkstraQuery(QueryData& data);

    void bidirectionalDijkstraGetPath(QueryData& data);

    void contractionHierachyQuery(QueryData& data);

    void createHubLabels(int threshold = std::numeric_limits<int>::max());

    void hubLabelQuery(QueryData& data);

    double averageLabelSize();

    int maxLabelSize();

    int numHubLabelsInRange(int upper, int lower);

    std::vector<int> createShortestPathCover(int threshold);

    bool verifyShortestPathCover(std::vector<int>& shortest_path_cover, int threshold);

    std::vector<int> reducePathCover(std::vector<int>& path_cover, int threshold);
    bool forwardDijkstraSearch(LowerBoundData& lb_data);
    bool backwardDijkstraSearch(LowerBoundData& lb_data);

    std::vector<int> lowerBound(std::vector<int>& shortest_path_cover, int threshold);

    std::vector<std::vector<Edge>>& getGraph() { return m_graph; }

    void setNumThreads(int num_of_threads) { m_num_threads = num_of_threads; }

    void clearHubLabel() {
        std::vector<uint64_t>().swap(m_fwd_indices);
        std::vector<uint64_t>().swap(m_bwd_indices);
        std::vector<std::pair<int, int>>().swap(m_fwd_hub_labels);
        std::vector<std::pair<int, int>>().swap(m_bwd_hub_labels);
    }

   private:
    bool m_ch_available;  // ch = Contraction Hierarchy
    int m_num_nodes;
    int m_num_threads;
    std::vector<std::vector<Edge>> m_graph;
    std::vector<std::vector<Edge>> m_reverse_graph;
    std::vector<int> m_node_level;
    std::vector<ContractionData> m_contr_data;

    std::vector<int> m_level_indices_sorted;
    std::vector<int> m_node_indices;
    std::vector<uint64_t> m_fwd_indices;
    std::vector<uint64_t> m_bwd_indices;
    std::vector<std::pair<int, int>> m_fwd_hub_labels;
    std::vector<std::pair<int, int>> m_bwd_hub_labels;
    // std::vector<std::vector<std::pair<int, int>>> m_fwd_hub_labels;
    // std::vector<std::vector<std::pair<int, int>>> m_bwd_hub_labels;

    void createReverseGraphCH();

    void createReverseGraphNormal();

    int greatCircleDistance(double lat_1, double lon_1, double lat_2, double lon_2);

    void createCH(Heuristic heuristic);

    int inOutProductHeuristic(std::vector<bool>& contracted, int node);

    int edgeDifferenceHeuristic(std::vector<bool>& contracted, int node);

    int weightedCostHeuristic(std::vector<bool>& contracted, int node);

    int altWeightedCostHeuristic(std::vector<bool>& contracted, int node, std::vector<int>& longest_path_fwd,
                                 std::vector<int>& longest_path_bwd);

    int deletedNeighboursHeuristic(std::vector<bool>& contracted, int node, std::vector<int>& num_deleted_neighbours);

    int microsoftHeuristic(std::vector<bool>& contracted, int node, int cur_level);

    void contractNode(std::vector<bool>& contracted, int contracted_node);

    void contractionDijkstra(int start, int contracted_node, std::vector<bool>& contracted, int num_outgoing,
                             int max_distance);

    void lowerBoundDijkstra(LowerBoundData& lb_data);

    bool pathCoverVerificationDijkstra(LowerBoundData& lb_data);

    int simplifiedHubLabelQuery(std::vector<std::pair<int, int>>& fwd_labels, int node);

    int simplifiedHubLabelQuery(int node, std::vector<std::pair<int, int>>& bwd_labels);
};

}  // namespace olsp