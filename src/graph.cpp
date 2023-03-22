#include "graph.h"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>

namespace olsp {

Graph::Graph(const std::string& path, ReadMode read_mode, bool ch_available) : m_ch_available(ch_available) {
    readGraph(path, read_mode);
    m_num_nodes = m_graph.size();

    if (read_mode == ReadMode::NORMAL && m_ch_available) {
        createCH();  // TODO:
    }

    createReverseGraph();
}

Graph::Graph(std::vector<std::vector<Edge>> graph) : m_graph(graph) {
    createReverseGraph();
    m_num_nodes = m_graph.size();
}

void Graph::readGraph(const std::string& path, ReadMode read_mode) {
    std::cout << "Started reading graph file." << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    std::ifstream infile(path);
    std::string line = "#";

    // skip the metadata which begins with #
    while (line[0] == '#') getline(infile, line);

    getline(infile, line);
    int num_nodes = std::stoi(line);
    getline(infile, line);
    int num_edges = std::stoi(line);

    m_node_level.clear();

    // save node rank if CH Mode, else skip the node information
    if (read_mode == ReadMode::CONTRACTION_HIERARCHY) {
        m_node_level.resize(num_nodes);
        for (int i = 0; i < num_nodes; ++i) {
            getline(infile, line);
            std::stringstream ss(line);
            std::string s;

            // skip unused fields
            getline(ss, s, ' ');
            getline(ss, s, ' ');
            getline(ss, s, ' ');
            getline(ss, s, ' ');
            getline(ss, s, ' ');

            getline(ss, s, ' ');
            m_node_level[i] = std::stoi(s);
        }
    } else {
        for (int i = 0; i < num_nodes; ++i) getline(infile, line);
    }

    m_graph.clear();
    m_graph.resize(num_nodes);

    // read edge information
    for (int i = 0; i < num_edges; ++i) {
        getline(infile, line);
        std::stringstream ss(line);

        std::string s;
        getline(ss, s, ' ');
        int src = std::stoi(s);
        getline(ss, s, ' ');
        int target = std::stoi(s);
        getline(ss, s, ' ');
        int cost = std::stoi(s);

        if (read_mode == ReadMode::NORMAL) {
            m_graph[src].push_back(Edge{target, cost});
        } else {
            // skip unused fields
            getline(ss, s, ' ');
            getline(ss, s, ' ');

            getline(ss, s, ' ');
            int child_1 = std::stoi(s);
            getline(ss, s, ' ');
            int child_2 = std::stoi(s);
            m_graph[src].push_back(Edge{target, cost, child_1, child_2});
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Finished reading graph file. Took " << elapsed.count() << " milliseconds " << std::endl;
}

void Graph::createReverseGraph() {
    std::cout << "Started creating reverse graph." << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    if (m_graph.empty()) {
        std::cout << "Can't create reverse graph, because graph is empty" << std::endl;
        return;
    }

    m_reverse_graph.clear();
    m_reverse_graph.resize(m_num_nodes);

    if (m_ch_available)
        createReverseGraphCH();
    else
        createReverseGraphNormal();

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Finished creating reverse graph. Took " << elapsed.count() << " milliseconds " << std::endl;
}

int Graph::dijkstraQuery(std::vector<std::vector<Edge>>& graph, int start, int end) {
    if (start < 0 || end < 0) {
        std::cout << "Invalid start or end nodes!" << std::endl;
        return -1;
    }

    int distance = std::numeric_limits<int>::max();
    std::vector<bool> visited(graph.size(), false);
    std::vector<int> distances(graph.size(), std::numeric_limits<int>::max());
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;

    distances[start] = 0;
    pq.push(std::make_pair(0, start));

    while (!pq.empty()) {
        std::pair<int, int> cur_node = pq.top();
        pq.pop();
        visited[cur_node.second] = true;

        if (cur_node.second == end) return distances[end];

        for (Edge& e : graph[cur_node.second]) {
            if (!visited[e.m_target] && distances[e.m_target] > distances[cur_node.second] + e.m_cost) {
                distances[e.m_target] = distances[cur_node.second] + e.m_cost;
                pq.push(std::make_pair(distances[e.m_target], e.m_target));
            }
        }
    }
    return -1;
}

void Graph::bidirectionalDijkstraQuery(QueryData& data) {
    if (data.m_start < 0 || data.m_end < 0) {
        std::cout << "Invalid start or end nodes!" << std::endl;
        return;
    }

    std::cout << "Started bidirectional Dijkstra run." << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // reset data for bidirectional dijkstra
    data.m_fwd_prev.clear();
    data.m_bwd_prev.clear();

    if (data.m_path_needed) {
        data.m_fwd_prev.resize(m_num_nodes, -1);
        data.m_bwd_prev.resize(m_num_nodes, -1);
    }

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

                if (data.m_path_needed) data.m_fwd_prev[e.m_target] = fwd_node.second;
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

                if (data.m_path_needed) data.m_bwd_prev[e.m_target] = bwd_node.second;
            }

            if (visited_fwd[e.m_target] &&
                distances_bwd[bwd_node.second] + e.m_cost + distances_fwd[e.m_target] < data.m_distance) {
                data.m_distance = distances_bwd[bwd_node.second] + e.m_cost + distances_fwd[e.m_target];
                data.m_meeting_node = e.m_target;
            }
        }

        // termination condition
        if (distances_fwd[fwd_node.second] + distances_bwd[bwd_node.second] >= data.m_distance) {
            auto end = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
            std::cout << "Finished bidirectional Dijkstra run. Took " << elapsed.count() << " milliseconds "
                      << std::endl;
            return;
        }
    }
}

void Graph::bidirectionalDijkstraGetPath(QueryData& data) {
    if (!(data.m_distance < std::numeric_limits<int>::max() || data.m_distance > -1) || !data.m_path_needed) {
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
    std::cout << "Level of Meeting Node: " << m_node_level[data.m_meeting_node] << std::endl;
    for (size_t i = 0; i < data.m_shortest_path.size(); i++) {
        std::cout << "NodeID: " << data.m_shortest_path[i] << " Level: " << m_node_level[data.m_shortest_path[i]]
                  << std::endl;
    }
}

void Graph::contractionHierachyQuery(QueryData& data) {
    if (data.m_start < 0 || data.m_end < 0) {
        std::cout << "Invalid start or end nodes!" << std::endl;
        return;
    }

    std::cout << "Started CH Query." << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // reset data for bidirectional dijkstra
    data.m_fwd_prev.clear();
    data.m_bwd_prev.clear();

    if (data.m_path_needed) {
        data.m_fwd_prev.resize(m_num_nodes, -1);
        data.m_bwd_prev.resize(m_num_nodes, -1);
    }

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

    visited_fwd[data.m_start] = true;
    visited_bwd[data.m_end] = true;

    std::pair<int, int> fwd_node;
    std::pair<int, int> bwd_node;

    while (!fwd_pq.empty() || !bwd_pq.empty()) {
        while (!fwd_pq.empty()) {
            fwd_node = fwd_pq.top();
            fwd_pq.pop();
            if (visited_fwd[fwd_node.second] && fwd_node.first > distances_fwd[fwd_node.second]) continue;
            if (fwd_node.first > data.m_distance) break;

            // forward step
            for (Edge& e : m_graph[fwd_node.second]) {
                // relax edge
                if (m_node_level[e.m_target] <= m_node_level[fwd_node.second]) continue;

                if (!visited_fwd[e.m_target] || distances_fwd[e.m_target] > fwd_node.first + e.m_cost) {
                    distances_fwd[e.m_target] = fwd_node.first + e.m_cost;
                    fwd_pq.push(std::make_pair(distances_fwd[e.m_target], e.m_target));
                    visited_fwd[e.m_target] = true;

                    if (data.m_path_needed) data.m_fwd_prev[e.m_target] = fwd_node.second;
                }

                if (visited_bwd[e.m_target] &&
                    distances_fwd[fwd_node.second] + e.m_cost + distances_bwd[e.m_target] < data.m_distance) {
                    data.m_distance = distances_fwd[fwd_node.second] + e.m_cost + distances_bwd[e.m_target];
                    data.m_meeting_node = e.m_target;
                }
            }

            break;
        }

        while (!bwd_pq.empty()) {
            bwd_node = bwd_pq.top();
            bwd_pq.pop();
            if (visited_bwd[bwd_node.second] && bwd_node.first > distances_fwd[bwd_node.second]) continue;
            if (bwd_node.first > data.m_distance) break;

            // backward step
            for (Edge& e : m_reverse_graph[bwd_node.second]) {
                // relax edge
                if (m_node_level[e.m_target] <= m_node_level[bwd_node.second]) continue;

                if (!visited_bwd[e.m_target] || distances_bwd[e.m_target] > bwd_node.first + e.m_cost) {
                    distances_bwd[e.m_target] = bwd_node.first + e.m_cost;
                    bwd_pq.push(std::make_pair(distances_bwd[e.m_target], e.m_target));
                    visited_bwd[e.m_target] = true;

                    if (data.m_path_needed) data.m_bwd_prev[e.m_target] = bwd_node.second;
                }

                if (visited_fwd[e.m_target] &&
                    distances_bwd[bwd_node.second] + e.m_cost + distances_fwd[e.m_target] < data.m_distance) {
                    data.m_distance = distances_bwd[bwd_node.second] + e.m_cost + distances_fwd[e.m_target];
                    data.m_meeting_node = e.m_target;
                }
            }
            break;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
    std::cout << "Finished CH Query. Took " << elapsed.count() << " microseconds " << std::endl;
}

void Graph::createHubLabels() {
    std::cout << "Started creating hub labels." << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    if (m_graph.empty() || m_reverse_graph.empty()) {
        std::cout << "Can't create hub labels, because graph or reverse grpah is empty" << std::endl;
        return;
    }

    m_fwd_hub_labels.clear();
    m_fwd_hub_labels.resize(m_num_nodes);
    m_bwd_hub_labels.clear();
    m_bwd_hub_labels.resize(m_num_nodes);

    // sort the levels, but don't change the original vector
    std::vector<int> indices_sorted(m_num_nodes);
    std::iota(indices_sorted.begin(), indices_sorted.end(), 0);
    std::sort(indices_sorted.begin(), indices_sorted.end(),
              [&](int i, int j) { return m_node_level[i] > m_node_level[j]; });

    for (int i = 0; i < m_num_nodes; ++i) {
        int node = indices_sorted[i];
        if (i == 0) {
            m_fwd_hub_labels[node].push_back(std::make_pair(node, 0));
            m_bwd_hub_labels[node].push_back(std::make_pair(node, 0));
            continue;
        }

        // fwd lables
        for (Edge& e : m_graph[node]) {
            if (m_node_level[node] >= m_node_level[e.m_target]) continue;

            for (std::pair<int, int>& hub : m_fwd_hub_labels[e.m_target])
                m_fwd_hub_labels[node].push_back(std::make_pair(hub.first, hub.second + e.m_cost));
            m_fwd_hub_labels[node].push_back(std::make_pair(e.m_target, e.m_cost));
        }
        auto& fwd_labels = m_fwd_hub_labels[node];
        // remove duplicates
        std::sort(m_fwd_hub_labels[node].begin(), m_fwd_hub_labels[node].end(),
                  [](auto& left, auto& right) { return left.first < right.first; });
        for (auto iter = m_fwd_hub_labels[node].begin(); iter != m_fwd_hub_labels[node].end();) {
            auto iter_2 = m_fwd_hub_labels[node].end() - 1;
            if (std::distance(iter, iter_2) != 0 && iter->first == (iter + 1)->first) {
                if (iter->second > (iter + 1)->second) {
                    iter = m_fwd_hub_labels[node].erase(iter);
                } else {
                    iter = m_fwd_hub_labels[node].erase(iter + 1);
                    iter--;
                }
            } else {
                ++iter;
            }
        }

        // bwd labels
        for (Edge& e : m_reverse_graph[node]) {
            if (m_node_level[node] >= m_node_level[e.m_target]) continue;

            for (std::pair<int, int>& hub : m_bwd_hub_labels[e.m_target])
                m_bwd_hub_labels[node].push_back(std::make_pair(hub.first, hub.second + e.m_cost));
            m_bwd_hub_labels[node].push_back(std::make_pair(e.m_target, e.m_cost));
        }

        auto& bwd_labels = m_bwd_hub_labels[node];
        // remove duplicates
        std::sort(m_bwd_hub_labels[node].begin(), m_bwd_hub_labels[node].end(),
                  [](auto& left, auto& right) { return left.first < right.first; });
        for (auto iter = m_bwd_hub_labels[node].begin(); iter != m_bwd_hub_labels[node].end();) {
            auto iter_2 = m_bwd_hub_labels[node].end() - 1;
            if (std::distance(iter, iter_2) != 0 && iter->first == (iter + 1)->first) {
                if (iter->second > (iter + 1)->second) {
                    iter = m_bwd_hub_labels[node].erase(iter);
                } else {
                    iter = m_bwd_hub_labels[node].erase(iter + 1);
                    iter--;
                }
            } else {
                ++iter;
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Finished creating hub labels. Took " << elapsed.count() << " milliseconds " << std::endl;
}

void Graph::hubLabelQuery(QueryData& data) {
    if (data.m_start < 0 || data.m_end < 0) {
        std::cout << "Invalid start or end nodes!" << std::endl;
        return;
    }

    std::cout << "Started Hub Query." << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    data.m_distance = std::numeric_limits<int>::max();
    data.m_meeting_node = -1;
    for (auto& fwd_hub : m_fwd_hub_labels[data.m_start]) {
        for (auto& bwd_hub : m_bwd_hub_labels[data.m_end]) {
            if (fwd_hub.first == bwd_hub.first && fwd_hub.second + bwd_hub.second <= data.m_distance) {
                data.m_distance = fwd_hub.second + bwd_hub.second;
                data.m_meeting_node = fwd_hub.first;
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout << "Finished Hub Label query. Took " << elapsed.count() << " nanoseconds " << std::endl;
}

double Graph::averageLabelSize() {
    double avg_label_size = 0.0;

    for (int i = 0; i < m_num_nodes; i++) {
        avg_label_size += m_fwd_hub_labels[i].size();
        avg_label_size += m_bwd_hub_labels[i].size();
    }
    return avg_label_size / static_cast<double>(2 * m_num_nodes);
}

int Graph::maxLabelSize() {
    int max_label_size = 0;

    for (int i = 0; i < m_num_nodes; i++) {
        if (m_fwd_hub_labels[i].size() > max_label_size) max_label_size = m_fwd_hub_labels[i].size();
        if (m_bwd_hub_labels[i].size() > max_label_size) max_label_size = m_bwd_hub_labels[i].size();
    }
    return max_label_size;
}

void Graph::createReverseGraphNormal() {
    for (int i = 0; i < m_num_nodes; ++i) {
        for (Edge e : m_graph[i]) {
            int new_source = e.m_target;
            e.m_target = i;

            // TODO: I'm not sure about this (should shortcut order be changed)
            int child_1_copy = e.m_child_1;
            e.m_child_1 = e.m_child_2;
            e.m_child_2 = child_1_copy;

            m_reverse_graph[new_source].push_back(e);
        }
    }
}

void Graph::createReverseGraphCH() {
    for (int i = 0; i < m_num_nodes; ++i) {
        for (auto iter = m_graph[i].begin(); iter != m_graph[i].end();) {
            if (m_node_level[i] > m_node_level[iter->m_target]) {
                Edge edge = *iter;
                int new_source = edge.m_target;
                edge.m_target = i;
                // TODO: I'm not sure about this (should shortcut order be changed)
                int child_1_copy = edge.m_child_1;
                edge.m_child_1 = edge.m_child_2;
                edge.m_child_2 = child_1_copy;

                m_reverse_graph[new_source].push_back(edge);
                iter = m_graph[i].erase(iter);
            } else {
                ++iter;
            }
        }
    }
}

void Graph::createCH() {
    // TODO:

    std::cout << "Started creating CH." << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    if (m_graph.empty()) {
        std::cout << "Can't create CH, because graph is empty." << std::endl;
        return;
    }
    /*
        auto graph = m_graph;
        while(graph.size() != 0){
            std::vector<int> edge_diff(graph.size());
            for(int i = 0)
        }
    */
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Finished creating CH. Took " << elapsed.count() << " milliseconds " << std::endl;
}

}  // namespace olsp
