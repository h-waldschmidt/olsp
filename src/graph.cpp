#include "graph.h"

#include <iostream>
#include <queue>

namespace olsp {

Graph::Graph(const std::string& path, ReadMode read_mode) {
    readGraph(path, m_graph, read_mode);
    createReverseGraph(m_graph, m_reverse_graph);
}

Graph::Graph(std::vector<std::vector<Edge>> graph) : m_graph(graph) { createReverseGraph(m_graph, m_reverse_graph); }

void Graph::readGraph(const std::string& path, std::vector<std::vector<Edge>>& graph, ReadMode read_mode) {
    // TODO:
}

void Graph::createReverseGraph(const std::vector<std::vector<Edge>>& graph,
                               std::vector<std::vector<Edge>>& reverse_graph) {
    if (graph.empty()) {
        std::cout << "Can't create reverse graph, because graph is empty" << std::endl;
        return;
    }

    reverse_graph.clear();
    reverse_graph.resize(graph.size());

    for (int i = 0; i < graph.size(); ++i) {
        for (Edge e : graph[i]) {
            int new_source = e.m_target;
            e.m_target = i;

            // TODO: I'm not sure about this (should shortcut order be changed)
            int child_1_copy = e.m_child_1;
            e.m_child_1 = e.m_child_2;
            e.m_child_2 = child_1_copy;

            reverse_graph[new_source].push_back(e);
        }
    }
}

void Graph::bidirectionalDijkstraCalculateDistance(BiDirectionalDijkstraData& data) {
    if (data.m_start < 0 || data.m_end < 0) {
        std::cout << "Invalid start or end nodes!" << std::endl;
        return;
    }

    // reset data for bidirectional dijkstra
    data.m_fwd_prev.clear();
    data.m_bwd_prev.clear();
    data.m_fwd_prev.resize(m_num_nodes, -1);
    data.m_bwd_prev.resize(m_num_nodes, -1);
    data.m_distance = std::numeric_limits<int>::max();
    data.m_meeting_node = -1;

    std::vector<bool> visited_fwd(m_num_nodes, false);
    std::vector<bool> visited_bwd(m_num_nodes, false);

    std::vector<int> distances_fwd(m_num_nodes, std::numeric_limits<int>::max());
    std::vector<int> distances_bwd(m_num_nodes, std::numeric_limits<int>::max());

    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>>
        fwd_pq;
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>>
        bwd_pq;

    distances_fwd[data.m_start] = 0;
    distances_bwd[data.m_end] = 0;

    // first corresponds to distance and second is node index
    fwd_pq.push(std::make_pair(0, data.m_start));
    bwd_pq.push(std::make_pair(0, data.m_end));

    while (!fwd_pq.empty() && !bwd_pq.empty()) {
        std::pair<int, int> fwd_node = fwd_pq.top();
        std::pair<int, int> bwd_node = bwd_pq.top();
        fwd_pq.pop();
        bwd_pq.pop();

        visited_fwd[fwd_node.second] = true;
        visited_bwd[bwd_node.second] = true;

        // forward step
        for (Edge& e : m_graph[fwd_node.second]) {
            // relax edge
            if (!visited_fwd[e.m_target] && distances_fwd[e.m_target] > distances_fwd[fwd_node.second] + e.m_cost) {
                distances_fwd[e.m_target] = distances_fwd[fwd_node.second] + e.m_cost;
                fwd_pq.push(std::make_pair(distances_fwd[e.m_target], e.m_target));
                data.m_fwd_prev[e.m_target] = fwd_node.second;
            }

            if (visited_bwd[e.m_target] &&
                distances_fwd[fwd_node.second] + e.m_cost + distances_bwd[e.m_target] < data.m_distance) {
                data.m_distance = distances_fwd[fwd_node.second] + e.m_cost + distances_bwd[e.m_target];
                data.m_meeting_node = e.m_target;
            }
        }

        // backward step
        for (Edge& e : m_reverse_graph[bwd_node.second]) {
            // relax edge
            if (!visited_bwd[e.m_target] && distances_bwd[e.m_target] > distances_bwd[bwd_node.second] + e.m_cost) {
                distances_bwd[e.m_target] = distances_bwd[bwd_node.second] + e.m_cost;
                bwd_pq.push(std::make_pair(distances_bwd[e.m_target], e.m_target));
                data.m_bwd_prev[e.m_target] = bwd_node.second;
            }

            if (visited_fwd[e.m_target] &&
                distances_bwd[bwd_node.second] + e.m_cost + distances_fwd[e.m_target] < data.m_distance) {
                data.m_distance = distances_bwd[bwd_node.second] + e.m_cost + distances_fwd[e.m_target];
                data.m_meeting_node = e.m_target;
            }
        }

        // termination condition
        if (distances_fwd[fwd_node.second] + distances_bwd[bwd_node.second] >= data.m_distance) return;
    }
}

void Graph::bidirectionalDijkstraGetPath(BiDirectionalDijkstraData& data) {
    if (!(data.m_distance < std::numeric_limits<int>::max() || data.m_distance > -1)) {
        std::cout << "Can't return path for invalid data!" << std::endl;
        return;
    }

    data.m_shortest_path.clear();

    std::vector<int> fwd_path;
    std::vector<int> bwd_path;

    int cur_node = data.m_meeting_node;
    while (cur_node != data.m_start) {
        cur_node = data.m_fwd_prev[cur_node];
        fwd_path.push_back(cur_node);
    }

    cur_node = data.m_meeting_node;
    while (cur_node != data.m_end) {
        cur_node = data.m_bwd_prev[cur_node];
        bwd_path.push_back(cur_node);
    }

    for (int i = fwd_path.size() - 1; i >= 0; --i) {
        data.m_shortest_path.push_back(fwd_path[i]);
    }

    data.m_shortest_path.push_back(data.m_meeting_node);

    for (int i = 0; i < bwd_path.size(); ++i) {
        data.m_shortest_path.push_back(bwd_path[i]);
    }
}

}  // namespace olsp
