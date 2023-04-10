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

struct QueryData {
    int m_start;
    int m_end;
    int m_meeting_node;
    int m_distance;
    bool m_path_needed;           // when the path is not needed the path related variables don't need to be filled`
    std::vector<int> m_fwd_prev;  // fwd = forward
    std::vector<int> m_bwd_prev;  // bwd = backward
    std::vector<int> m_shortest_path;

    QueryData(int start, int end, bool path_needed)
        : m_start(start),
          m_end(end),
          m_meeting_node(-1),
          m_distance(std::numeric_limits<int>::max()),
          m_path_needed(path_needed),
          m_fwd_prev(0),
          m_bwd_prev(0),
          m_shortest_path(0) {}
};

enum ReadMode { NORMAL = 0, CONTRACTION_HIERARCHY = 1 };

class Graph {
   public:
    Graph(const std::string& path, ReadMode read_mode, bool ch_available);
    Graph(std::vector<std::vector<Edge>> graph);  // TODO: Adjust constructors for normal mode and ch mode
    ~Graph() = default;

    void readGraph(const std::string& path, ReadMode read_mode);

    void createReverseGraph();

    static int dijkstraQuery(std::vector<std::vector<Edge>>& graph, int start, int end);

    void bidirectionalDijkstraQuery(QueryData& data);

    void bidirectionalDijkstraGetPath(QueryData& data);

    void contractionHierachyQuery(QueryData& data);

    void createHubLabels();

    void hubLabelQuery(QueryData& data);

    double averageLabelSize();

    int maxLabelSize();

    std::vector<int> createShortestPathCover(int threshold);

    std::vector<std::vector<Edge>>& getGraphVec() { return m_graph; }

   private:
    bool m_ch_available;  // ch = Contraction Hierarchy
    int m_num_nodes;
    std::vector<std::vector<Edge>> m_graph;
    std::vector<std::vector<Edge>> m_reverse_graph;
    std::vector<int> m_node_level;
    std::vector<std::vector<std::pair<int, int>>> m_fwd_hub_labels;
    std::vector<std::vector<std::pair<int, int>>> m_bwd_hub_labels;

    void createReverseGraphCH();

    void createReverseGraphNormal();

    void createCH();

    int inOutProductHeuristic(std::vector<bool>& contracted, int node);

    void contractNode(std::vector<bool>& contracted, int contracted_node);

    void contractionDijkstra(std::vector<int>& distances, int start, int contracted_node, std::vector<bool>& contracted,
                             std::vector<bool>& outgoing_nodes, int num_outgoing, int max_distance);

    int simplifiedHubLabelQuery(std::vector<std::pair<int, int>>& fwd_labels,
                                std::vector<std::pair<int, int>>& bwd_labels);
};

}  // namespace olsp
