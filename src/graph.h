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

struct BiDirectionalDijkstraData {
    int m_start;
    int m_end;
    int m_meeting_node;
    int m_distance;
    std::vector<int> m_fwd_prev;  // fwd = forward
    std::vector<int> m_bwd_prev;  // bwd = backward
    std::vector<int> m_shortest_path;

    BiDirectionalDijkstraData(int start, int end)
        : m_start(start),
          m_end(end),
          m_meeting_node(-1),
          m_distance(std::numeric_limits<int>::max()),
          m_fwd_prev(0),
          m_bwd_prev(0),
          m_shortest_path(0) {}
};

enum ReadMode { NORMAL = 0, CONTRACTION_HIERACHIES = 1 };

class Graph {
   public:
    Graph(const std::string& path, ReadMode read_mode);
    Graph(std::vector<std::vector<Edge>> graph);
    ~Graph() = default;

    static void readGraph(const std::string& path, std::vector<std::vector<Edge>>& graph, ReadMode read_mode);

    static void createReverseGraph(const std::vector<std::vector<Edge>>& graph,
                                   std::vector<std::vector<Edge>>& reverse_graph);

    void bidirectionalDijkstraCalculateDistance(BiDirectionalDijkstraData& data);

    void bidirectionalDijkstraGetPath(BiDirectionalDijkstraData& data);

   private:
    int m_num_nodes;
    std::vector<std::vector<Edge>> m_graph;
    std::vector<std::vector<Edge>> m_reverse_graph;
};

}  // namespace olsp
