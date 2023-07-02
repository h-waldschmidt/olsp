#include "graph.h"

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <unordered_set>

#define _USE_MATH_DEFINES
#include <math.h>

namespace olsp {

Graph::Graph(const std::string& path, ReadMode read_mode, bool ch_available, bool prune_graph, Heuristic ch_heuristic,
             DistanceMode dist_mode)
    : m_ch_available(ch_available) {
    readGraph(path, read_mode, dist_mode);
    m_num_nodes = m_graph.size();

    if (read_mode == ReadMode::NORMAL && m_ch_available) {
        createCH(ch_heuristic);  // TODO:
    }

    createReverseGraph(prune_graph);
}

Graph::Graph(std::vector<std::vector<Edge>> graph) : m_graph(graph) {
    createReverseGraph(false);
    m_num_nodes = m_graph.size();
}

void Graph::readGraph(const std::string& path, ReadMode read_mode, DistanceMode dist_mode) {
    std::cout << "Started reading graph file." << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();
    std::ifstream infile;
    try {
        infile.open(path);
        if (!infile.good()) throw std::runtime_error("File doesn't exist!");
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        exit(1);
    }

    std::string line = "#";

    // skip the metadata which begins with #
    while (line[0] == '#') getline(infile, line);

    getline(infile, line);
    int num_nodes = std::stoi(line);
    getline(infile, line);
    int num_edges = std::stoi(line);

    m_node_level.clear();
    std::vector<std::pair<double, double>> node_coords;

    // save node rank if CH Mode, else skip the node information
    // or save coords if distances in meters are required
    if (read_mode == ReadMode::CONTRACTION_HIERARCHY || dist_mode == DistanceMode::DISTANCE_METERS) {
        if (read_mode == ReadMode::CONTRACTION_HIERARCHY) m_node_level.resize(num_nodes);
        if (dist_mode == DistanceMode::DISTANCE_METERS) node_coords.resize(num_nodes);

        for (int i = 0; i < num_nodes; ++i) {
            getline(infile, line);
            std::stringstream ss(line);
            std::string s;

            // skip unused fields
            getline(ss, s, ' ');
            getline(ss, s, ' ');

            if (dist_mode == DistanceMode::DISTANCE_METERS) {
                std::pair<double, double> coords;
                getline(ss, s, ' ');
                coords.first = stod(s);
                getline(ss, s, ' ');
                coords.second = stod(s);
                node_coords[i] = coords;
            } else {
                getline(ss, s, ' ');
                getline(ss, s, ' ');
            }

            if (read_mode == ReadMode::CONTRACTION_HIERARCHY) {
                getline(ss, s, ' ');
                getline(ss, s, ' ');
                m_node_level[i] = std::stoi(s);
            }
        }
    } else {
        for (int i = 0; i < num_nodes; ++i) {
            getline(infile, line);
        }
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

        if (dist_mode == DistanceMode::DISTANCE_METERS)
            cost = greatCircleDistance(node_coords[src].first, node_coords[src].second, node_coords[target].first,
                                       node_coords[target].second);

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

void Graph::createReverseGraph(bool prune_graph) {
    std::cout << "Started creating reverse graph." << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    if (m_graph.empty()) {
        std::cout << "Can't create reverse graph, because graph is empty" << std::endl;
        return;
    }

    m_reverse_graph.clear();
    m_reverse_graph.resize(m_num_nodes);

    if (m_ch_available && prune_graph)
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
        if (fwd_node.first + bwd_node.first >= data.m_distance) {
            break;
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
    std::cout << "Finished bidirectional Dijkstra run. Took " << elapsed.count() << " microseconds " << std::endl;
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

    // std::cout << "Started CH Query." << std::endl;

    // auto begin = std::chrono::high_resolution_clock::now();

    // reset data for bidirectional dijkstra
    data.m_fwd_prev.clear();
    data.m_bwd_prev.clear();

    if (data.m_path_needed) {
        data.m_fwd_prev.resize(m_num_nodes, -1);
        data.m_bwd_prev.resize(m_num_nodes, -1);
    }

    data.m_distance = std::numeric_limits<int>::max();
    data.m_meeting_node = -1;

    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>>
        fwd_pq;
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>>
        bwd_pq;

    data.distances_fwd[data.m_start] = 0;
    data.distances_bwd[data.m_end] = 0;

    // first corresponds to distance and second is node index
    fwd_pq.push(std::make_pair(0, data.m_start));
    bwd_pq.push(std::make_pair(0, data.m_end));

    data.visited_fwd[data.m_start] = true;
    data.visited_bwd[data.m_end] = true;

    data.reset_fwd.push_back(data.m_start);
    data.reset_bwd.push_back(data.m_end);

    std::pair<int, int> fwd_node;
    std::pair<int, int> bwd_node;

    while (!fwd_pq.empty() || !bwd_pq.empty()) {
        while (!fwd_pq.empty()) {
            fwd_node = fwd_pq.top();
            fwd_pq.pop();
            if (data.visited_fwd[fwd_node.second] && fwd_node.first > data.distances_fwd[fwd_node.second]) continue;
            if (fwd_node.first > data.m_distance) break;

            // forward step
            for (Edge& e : m_graph[fwd_node.second]) {
                // relax edge
                if (m_node_level[e.m_target] <= m_node_level[fwd_node.second]) continue;

                if (!data.visited_fwd[e.m_target] || data.distances_fwd[e.m_target] > fwd_node.first + e.m_cost) {
                    if (!data.visited_fwd[e.m_target]) data.reset_fwd.push_back(e.m_target);

                    data.distances_fwd[e.m_target] = fwd_node.first + e.m_cost;
                    fwd_pq.push(std::make_pair(data.distances_fwd[e.m_target], e.m_target));
                    data.visited_fwd[e.m_target] = true;

                    if (data.m_path_needed) data.m_fwd_prev[e.m_target] = fwd_node.second;
                }

                if (data.visited_bwd[e.m_target] &&
                    data.distances_fwd[fwd_node.second] + e.m_cost + data.distances_bwd[e.m_target] < data.m_distance) {
                    data.m_distance = data.distances_fwd[fwd_node.second] + e.m_cost + data.distances_bwd[e.m_target];
                    data.m_meeting_node = e.m_target;
                }
            }

            break;
        }

        while (!bwd_pq.empty()) {
            bwd_node = bwd_pq.top();
            bwd_pq.pop();
            if (data.visited_bwd[bwd_node.second] && bwd_node.first > data.distances_fwd[bwd_node.second]) continue;
            if (bwd_node.first > data.m_distance) break;

            // backward step
            for (Edge& e : m_reverse_graph[bwd_node.second]) {
                // relax edge
                if (m_node_level[e.m_target] <= m_node_level[bwd_node.second]) continue;

                if (!data.visited_bwd[e.m_target] || data.distances_bwd[e.m_target] > bwd_node.first + e.m_cost) {
                    if (!data.visited_fwd[e.m_target]) data.reset_bwd.push_back(e.m_target);

                    data.distances_bwd[e.m_target] = bwd_node.first + e.m_cost;
                    bwd_pq.push(std::make_pair(data.distances_bwd[e.m_target], e.m_target));
                    data.visited_bwd[e.m_target] = true;

                    if (data.m_path_needed) data.m_bwd_prev[e.m_target] = bwd_node.second;
                }

                if (data.visited_fwd[e.m_target] &&
                    data.distances_bwd[bwd_node.second] + e.m_cost + data.distances_fwd[e.m_target] < data.m_distance) {
                    data.m_distance = data.distances_bwd[bwd_node.second] + e.m_cost + data.distances_fwd[e.m_target];
                    data.m_meeting_node = e.m_target;
                }
            }
            break;
        }

        // termination condition
        if (fwd_node.first >= data.m_distance || bwd_node.first >= data.m_distance) {
            break;
        }
    }

    // reset_data
    for (int& index : data.reset_fwd) {
        data.distances_fwd[index] = std::numeric_limits<int>::max();
        data.visited_fwd[index] = false;
    }
    data.reset_fwd.clear();
    for (int& index : data.reset_bwd) {
        data.distances_bwd[index] = std::numeric_limits<int>::max();
        data.visited_bwd[index] = false;
    }
    data.reset_bwd.clear();

    // auto end = std::chrono::high_resolution_clock::now();
    // auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - begin);
    // std::cout << "Finished CH Query. Took " << elapsed.count() << " microseconds " << std::endl;
}

void Graph::createHubLabels(int threshold) {
    std::cout << "Started creating hub labels." << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    if (m_graph.empty() || m_reverse_graph.empty()) {
        std::cout << "Can't create hub labels, because graph or reverse grpah is empty" << std::endl;
        return;
    }

    m_fwd_indices.clear();
    m_fwd_indices.resize(m_num_nodes + 1, 0);
    m_bwd_indices.clear();
    m_bwd_indices.resize(m_num_nodes + 1, 0);

    m_fwd_hub_labels.clear();
    m_bwd_hub_labels.clear();

    // sort the levels, but don't change the original vector
    m_level_indices_sorted.clear();
    m_level_indices_sorted.resize(m_num_nodes);
    std::iota(m_level_indices_sorted.begin(), m_level_indices_sorted.end(), 0);
    std::sort(m_level_indices_sorted.begin(), m_level_indices_sorted.end(),
              [&](int i, int j) { return m_node_level[i] > m_node_level[j]; });

    m_node_indices.clear();
    m_node_indices.resize(m_num_nodes);
    // save for each node its correspnding index
    for (int i = 0; i < m_num_nodes; ++i) {
        m_node_indices[m_level_indices_sorted[i]] = i;
    }

    std::vector<std::pair<int, int>> fwd_labels;
    std::vector<std::pair<int, int>> bwd_labels;

    for (int i = 0; i < m_num_nodes; ++i) {
        int node = m_level_indices_sorted[i];

        if (i == 0) {
            m_fwd_hub_labels.push_back(std::make_pair(node, 0));
            m_bwd_hub_labels.push_back(std::make_pair(node, 0));
            m_fwd_indices[1] = 1;
            m_bwd_indices[1] = 1;
            continue;
        }

        fwd_labels.clear();
        bwd_labels.clear();

        fwd_labels.push_back(std::make_pair(node, 0));
        bwd_labels.push_back(std::make_pair(node, 0));

        int test_1 = m_node_indices[node];
        int test_2 = m_level_indices_sorted[i];

        // fwd lables
        for (Edge& e : m_graph[node]) {
            if (m_node_level[node] >= m_node_level[e.m_target]) continue;

            for (uint64_t j = m_fwd_indices[m_node_indices[e.m_target]];
                 j < m_fwd_indices[m_node_indices[e.m_target] + 1]; j++) {
                /*
                if (i >= 16503330){
                    std::cout << "Node ID" << node << std::endl;
                    std::cout << "Target Node ID: " << e.m_target << std::endl;
                    std::cout << "Node Index: " << i << std::endl;
                    std::cout << "Target Node Index: " << m_node_indices[e.m_target] << std::endl;
                    std::cout << "Target Node fwd Index: " << m_fwd_indices[m_node_indices[e.m_target]] << std::endl;
                    std::cout << "Target Node bwd Index: " << m_bwd_indices[m_node_indices[e.m_target]] <<std::endl;
                    std::cout << "Num Fwd Label: " << m_fwd_hub_labels.size() << std::endl;
                    std::cout << "Do Fwd Label Produce segfault: " << j << std::endl;
                    std::cout << "Fwd Label: " << m_fwd_hub_labels[j].first << std::endl;
                    std::cout << "End Debug Output" << std::endl;
                    }
                */
                if (m_fwd_hub_labels[j].second + e.m_cost <= threshold)
                    fwd_labels.push_back(
                        std::make_pair(m_fwd_hub_labels[j].first, m_fwd_hub_labels[j].second + e.m_cost));
            }
        }
        // remove duplicates
        std::sort(fwd_labels.begin(), fwd_labels.end(),
                  [](auto& left, auto& right) { return left.first < right.first; });
        for (auto iter = fwd_labels.begin(); iter != fwd_labels.end();) {
            auto iter_2 = fwd_labels.end() - 1;
            if (std::distance(iter, iter_2) != 0 && iter->first == (iter + 1)->first) {
                if (iter->second >= (iter + 1)->second) {
                    iter = fwd_labels.erase(iter);
                } else {
                    iter = fwd_labels.erase(iter + 1);
                    --iter;
                }
            } else {
                ++iter;
            }
        }

        for (auto iter = fwd_labels.begin(); iter != fwd_labels.end();) {
            int best_dist = std::numeric_limits<int>::max();
            if (iter->first != node) best_dist = simplifiedHubLabelQuery(fwd_labels, iter->first);
            if (best_dist < iter->second)
                iter = fwd_labels.erase(iter);
            else
                ++iter;
        }

        // bwd lables
        for (Edge& e : m_reverse_graph[node]) {
            if (m_node_level[node] >= m_node_level[e.m_target]) continue;

            for (uint64_t j = m_bwd_indices[m_node_indices[e.m_target]];
                 j < m_bwd_indices[m_node_indices[e.m_target] + 1]; j++) {
                if (m_bwd_hub_labels[j].second + e.m_cost <= threshold)
                    bwd_labels.push_back(
                        std::make_pair(m_bwd_hub_labels[j].first, m_bwd_hub_labels[j].second + e.m_cost));
            }
        }
        // remove duplicates
        std::sort(bwd_labels.begin(), bwd_labels.end(),
                  [](auto& left, auto& right) { return left.first < right.first; });
        for (auto iter = bwd_labels.begin(); iter != bwd_labels.end();) {
            auto iter_2 = bwd_labels.end() - 1;
            if (std::distance(iter, iter_2) != 0 && iter->first == (iter + 1)->first) {
                if (iter->second >= (iter + 1)->second) {
                    iter = bwd_labels.erase(iter);
                } else {
                    iter = bwd_labels.erase(iter + 1);
                    --iter;
                }
            } else {
                ++iter;
            }
        }
        for (auto iter = bwd_labels.begin(); iter != bwd_labels.end();) {
            int best_dist = std::numeric_limits<int>::max();
            if (iter->first != node) best_dist = simplifiedHubLabelQuery(iter->first, bwd_labels);
            if (best_dist < iter->second)
                iter = bwd_labels.erase(iter);
            else
                ++iter;
        }
        // update hub label data
        m_fwd_indices[i + 1] = m_fwd_indices[i] + static_cast<uint64_t>(fwd_labels.size());
        m_bwd_indices[i + 1] = m_bwd_indices[i] + static_cast<uint64_t>(bwd_labels.size());

        for (auto label : fwd_labels) {
            m_fwd_hub_labels.push_back(label);
        }
        for (auto label : bwd_labels) {
            m_bwd_hub_labels.push_back(label);
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

    int fwd_node_index = m_fwd_indices[m_node_indices[data.m_start]];
    int fwd_next_index = m_fwd_indices[m_node_indices[data.m_start] + 1];

    int bwd_node_index = m_bwd_indices[m_node_indices[data.m_end]];
    int bwd_next_index = m_bwd_indices[m_node_indices[data.m_end] + 1];

    while (fwd_node_index < fwd_next_index && bwd_node_index < bwd_next_index) {
        if (m_fwd_hub_labels[fwd_node_index].first == m_bwd_hub_labels[bwd_node_index].first) {
            if (m_fwd_hub_labels[fwd_node_index].second + m_bwd_hub_labels[bwd_node_index].second < data.m_distance) {
                data.m_meeting_node = m_fwd_hub_labels[fwd_node_index].first;
                data.m_distance = m_fwd_hub_labels[fwd_node_index].second + m_bwd_hub_labels[bwd_node_index].second;
            }
            ++fwd_node_index;
            ++bwd_node_index;
        } else if (m_fwd_hub_labels[fwd_node_index].first < m_bwd_hub_labels[bwd_node_index].first) {
            ++fwd_node_index;
        } else {
            ++bwd_node_index;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    std::cout << "Finished Hub Label query. Took " << elapsed.count() << " nanoseconds " << std::endl;
}

double Graph::averageLabelSize() {
    double avg_label_size = 0.0;

    for (int i = 0; i < m_num_nodes; i++) {
        avg_label_size += m_fwd_indices[i + 1] - m_fwd_indices[i];
        avg_label_size += m_bwd_indices[i + 1] - m_bwd_indices[i];
    }
    return avg_label_size / static_cast<double>(2 * m_num_nodes);
}

int Graph::maxLabelSize() {
    int max_label_size = 0;

    for (int i = 0; i < m_num_nodes; i++) {
        if (m_fwd_indices[i + 1] - m_fwd_indices[i] > max_label_size)
            max_label_size = m_fwd_indices[i + 1] - m_fwd_indices[i];
        if (m_bwd_indices[i + 1] - m_bwd_indices[i] > max_label_size)
            max_label_size = m_bwd_indices[i + 1] - m_bwd_indices[i];
    }
    return max_label_size;
}

int Graph::numHubLabelsInRange(int lower, int upper) {
    int count = 0;

    for (auto& fwd_labels : m_fwd_hub_labels) {
        if (fwd_labels.second > lower && fwd_labels.second < upper) ++count;
    }
    for (auto& bwd_labels : m_bwd_hub_labels) {
        if (bwd_labels.second > lower && bwd_labels.second < upper) ++count;
    }

    return count;
}

std::vector<int> Graph::createShortestPathCover(int threshold) {
    std::cout << "Started creating Path Cover." << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    std::unordered_set<int> path_cover_set;
    for (auto& fwd_labels : m_fwd_hub_labels) {
        if (fwd_labels.second > static_cast<double>(threshold) * (static_cast<double>(5) / static_cast<double>(10)) &&
            fwd_labels.second < threshold)
            path_cover_set.emplace(fwd_labels.first);
    }

    for (auto& bwd_labels : m_bwd_hub_labels) {
        if (bwd_labels.second > static_cast<double>(threshold) * (static_cast<double>(5) / static_cast<double>(10)) &&
            bwd_labels.second < threshold)
            path_cover_set.emplace(bwd_labels.first);
    }

    // convert to vector
    std::vector<int> path_cover_vec;
    path_cover_vec.insert(path_cover_vec.end(), path_cover_set.begin(), path_cover_set.end());

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Finished creating Path Cover. Took " << elapsed.count() << " milliseconds" << std::endl;

    return path_cover_vec;
}

bool Graph::verifyShortestPathCover(std::vector<int>& shortest_path_cover, int threshold) {
    // TODO:
    std::cout << "Started verifiying path cover." << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    LowerBoundData lb_data(m_num_nodes);
    lb_data.m_threshold = threshold;

    for (int& node : shortest_path_cover) lb_data.m_marked[node] = true;

    for (int i = 0; i < m_num_nodes; ++i) {
        lb_data.m_start_node = i;
        bool covered = pathCoverVerificationDijkstra(lb_data);
        if (!covered) return false;

        // reset data for next run
        for (int& node : lb_data.m_reset_previous_node) {
            lb_data.m_distances[node] = std::numeric_limits<int>::max();
        }
        lb_data.m_reset_previous_node.clear();

        if (i % 10000 == 0) {
            std::cout << "Finished " << i << "\n";
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
    std::cout << "Finished verifying path cover. Took " << elapsed.count() << " seconds" << std::endl;
    return true;
}

std::vector<int> Graph::reducePathCover(std::vector<int>& path_cover, int threshold) {
    // TODO:
    std::cout << "Started reducing path cover." << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    std::vector<int> new_path_cover;

    LowerBoundData lb_data(m_num_nodes);
    lb_data.m_threshold = threshold;

    for (int& node : path_cover) lb_data.m_marked[node] = true;

    for (int& node : path_cover) {
        lb_data.m_start_node = node;
        bool forward = forwardDijkstraSearch(lb_data);
        // reset data for next run
        for (int& node : lb_data.m_reset_previous_node) {
            lb_data.m_distances[node] = std::numeric_limits<int>::max();
        }
        lb_data.m_reset_previous_node.clear();

        bool backward = backwardDijkstraSearch(lb_data);
        // reset data for next run
        for (int& node : lb_data.m_reset_previous_node) {
            lb_data.m_distances[node] = std::numeric_limits<int>::max();
        }
        lb_data.m_reset_previous_node.clear();

        if (!forward && !backward)
            lb_data.m_marked[node] = false;
        else
            new_path_cover.push_back(node);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
    std::cout << "Finished reducing path cover. Took " << elapsed.count() << " seconds" << std::endl;
    return new_path_cover;
}

bool Graph::forwardDijkstraSearch(LowerBoundData& lb_data) {
    // TODO:
    std::priority_queue<std::pair<int, std::pair<int, bool>>, std::vector<std::pair<int, std::pair<int, bool>>>,
                        std::greater<std::pair<int, std::pair<int, bool>>>>
        pq;

    lb_data.m_distances[lb_data.m_start_node] = 0;
    lb_data.m_reset_previous_node.push_back(lb_data.m_start_node);
    pq.push(std::make_pair(0, std::make_pair(lb_data.m_start_node, false)));

    while (!pq.empty()) {
        auto cur_node = pq.top();
        pq.pop();

        if (lb_data.m_distances[cur_node.second.first] != cur_node.first) continue;

        if (lb_data.m_distances[cur_node.second.first] > static_cast<double>(lb_data.m_threshold) * 0.4 &&
            cur_node.second.second)
            return false;

        if (static_cast<double>(cur_node.first) > static_cast<double>(lb_data.m_threshold) * 0.5) break;

        for (Edge& e : m_graph[cur_node.second.first]) {
            if (m_node_level[cur_node.second.first] > m_node_level[e.m_target]) continue;

            if (lb_data.m_distances[e.m_target] > lb_data.m_distances[cur_node.second.first] + e.m_cost) {
                if (lb_data.m_distances[e.m_target] == std::numeric_limits<int>::max())
                    lb_data.m_reset_previous_node.push_back(e.m_target);
                lb_data.m_distances[e.m_target] = lb_data.m_distances[cur_node.second.first] + e.m_cost;
                lb_data.m_previous_node[e.m_target] = cur_node.second.first;
                bool covered;
                if (lb_data.m_marked[e.m_target])
                    covered = true;
                else
                    covered = cur_node.second.second;
                pq.push(std::make_pair(lb_data.m_distances[e.m_target], std::make_pair(e.m_target, covered)));
            }
        }
    }

    return true;
}

bool Graph::backwardDijkstraSearch(LowerBoundData& lb_data) {
    // TODO:
    std::priority_queue<std::pair<int, std::pair<int, bool>>, std::vector<std::pair<int, std::pair<int, bool>>>,
                        std::greater<std::pair<int, std::pair<int, bool>>>>
        pq;

    lb_data.m_distances[lb_data.m_start_node] = 0;
    lb_data.m_reset_previous_node.push_back(lb_data.m_start_node);
    pq.push(std::make_pair(0, std::make_pair(lb_data.m_start_node, false)));

    while (!pq.empty()) {
        auto cur_node = pq.top();
        pq.pop();

        if (lb_data.m_distances[cur_node.second.first] != cur_node.first) continue;

        if (lb_data.m_distances[cur_node.second.first] > static_cast<double>(lb_data.m_threshold) * 0.4 &&
            cur_node.second.second)
            return false;

        if (static_cast<double>(cur_node.first) > static_cast<double>(lb_data.m_threshold) * 0.5) break;

        for (Edge& e : m_reverse_graph[cur_node.second.first]) {
            if (m_node_level[cur_node.second.first] < m_node_level[e.m_target]) continue;

            if (lb_data.m_distances[e.m_target] > lb_data.m_distances[cur_node.second.first] + e.m_cost) {
                if (lb_data.m_distances[e.m_target] == std::numeric_limits<int>::max())
                    lb_data.m_reset_previous_node.push_back(e.m_target);
                lb_data.m_distances[e.m_target] = lb_data.m_distances[cur_node.second.first] + e.m_cost;
                lb_data.m_previous_node[e.m_target] = cur_node.second.first;
                bool covered;
                if (lb_data.m_marked[e.m_target])
                    covered = true;
                else
                    covered = cur_node.second.second;
                pq.push(std::make_pair(lb_data.m_distances[e.m_target], std::make_pair(e.m_target, covered)));
            }
        }
    }

    return true;
}

std::vector<int> intersection(std::vector<int> v1, std::vector<int> v2) {
    std::vector<int> v3;

    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v3));
    return v3;
}

std::vector<int> Graph::lowerBound(std::vector<int>& shortest_path_cover, int threshold) {
    std::cout << "Started calculating Lower Bound." << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    LowerBoundData lb_data(m_num_nodes);
    lb_data.m_threshold = threshold;

    // get random node out of shortest_path_cover
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(0, static_cast<int>(shortest_path_cover.size()) - 1);

    int spc_nodes_covered = 0;  // spc = shortest path cover
    std::vector<bool> spc_node_covered(shortest_path_cover.size(), false);
    std::vector<int> lower_bound;

    std::vector<bool> node_in_spc(m_num_nodes, false);
    for (int& node : shortest_path_cover) node_in_spc[node] = true;

    // std::vector<std::vector<int>> shortest_paths_backwards;

    // std::vector<int> shortest_path;
    std::vector<int> end_nodes;

    while (spc_nodes_covered != shortest_path_cover.size()) {
        int cur_spc_node = dist(rng);
        while (spc_node_covered[cur_spc_node]) cur_spc_node = dist(rng);
        spc_node_covered[cur_spc_node] = true;
        lb_data.m_start_node = shortest_path_cover[cur_spc_node];

        if (lb_data.m_marked[lb_data.m_start_node]) {
            spc_nodes_covered++;
            continue;
        }

        // run dijkstra until threshold reached
        lowerBoundDijkstra(lb_data);

        // look at paths with length greater than threshold
        // 1. first find all end nodes with path longer than threshold
        for (int i = 0; i < lb_data.m_reset_previous_node.size(); ++i) {
            int id = lb_data.m_reset_previous_node[i];
            if (lb_data.m_distances[id] < 2 * threshold && lb_data.m_distances[id] >= threshold) {
                end_nodes.push_back(id);

                /*
                if (lb_data.m_distances[id] !=
                    simplifiedHubLabelQuery(m_fwd_hub_labels[lb_data.m_start_node], m_bwd_hub_labels[id])) {
                    std::cout << "Different Distances" << std::endl;
                }
                */
            }
            lb_data.m_distances[id] = std::numeric_limits<int>::max();
        }

        // 2. backtrack paths through end node
        int selected_end_node = -1;
        int lowest_num_spc_nodes = std::numeric_limits<int>::max();
        for (int& node : end_nodes) {
            int cur_num_spc_nodes = 0;
            int cur_node = node;
            if (lb_data.m_marked[cur_node]) continue;

            while (cur_node != lb_data.m_start_node) {
                int previous = lb_data.m_previous_node[cur_node];
                if (lb_data.m_marked[previous]) break;

                if (node_in_spc[cur_node]) ++cur_num_spc_nodes;

                cur_node = previous;
            }

            if (cur_node == lb_data.m_start_node && !lb_data.m_marked[cur_node] &&
                cur_num_spc_nodes < lowest_num_spc_nodes) {
                lowest_num_spc_nodes = cur_num_spc_nodes;
                selected_end_node = node;
            }
        }

        end_nodes.clear();

        // 3. go through path of selected end node
        if (selected_end_node != -1) {
            int cur_node = selected_end_node;
            // shortest_path.push_back(cur_node);
            lb_data.m_marked[cur_node] = true;
            while (cur_node != lb_data.m_start_node) {
                int previous = lb_data.m_previous_node[cur_node];
                lb_data.m_marked[previous] = true;
                cur_node = previous;
                // shortest_path.push_back(cur_node);
            }
            // shortest_paths_backwards.push_back(shortest_path);
            // shortest_path.clear();
            lower_bound.push_back(lb_data.m_start_node);
        }

        // reset data for next run
        // for (int& node : lb_data.m_reset_previous_node) lb_data.m_previous_node[node] = -1;
        lb_data.m_reset_previous_node.clear();

        // std::cout << "Finished: " << spc_nodes_covered << "\n";

        ++spc_nodes_covered;
    }

    // check for duplicates in paths
    /*
    for (int i = 0; i < shortest_paths_backwards.size(); i++) {
        for (int j = i + 1; j < shortest_paths_backwards.size(); j++) {
            auto inters = intersection(shortest_paths_backwards[i], shortest_paths_backwards[j]);
            if (inters.size() != 0) {
                std::cout << "Intersection not empty: "
                          << "i: " << i << " j: " << j << " size: " << inters.size() << "\n";
                std::cout << "Start Node: " << shortest_paths_backwards[i][shortest_paths_backwards[i].size() - 1]
                          << " End Node: " << shortest_paths_backwards[i][0] << "\n";
                for (int k = 0; k < inters.size(); ++k) {
                    std::cout << "Inside Intersection " << k << " : " << inters[k] << "\n";
                }
            }
        }
    }
    */

    /*
    std::cout << "Hier noch die Pfade: "
              << "\n";
    // print paths
    for (int i = 0; i < shortest_paths_backwards.size(); i++) {
        for (int j = 0; j < shortest_paths_backwards[i].size(); j++) {
            std::cout << shortest_paths_backwards[i][j] << " ";
        }
        std::cout << "\n";
    }
    */

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
    std::cout << "Finished calculating Lower Bound. Took " << elapsed.count() << " seconds" << std::endl;

    return lower_bound;
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

int Graph::greatCircleDistance(double lat_1, double lon_1, double lat_2, double lon_2) {
    // degrees to radians
    lat_1 *= (M_PI / 180.0);
    lon_1 *= (M_PI / 180.0);
    lat_2 *= (M_PI / 180.0);
    lon_2 *= (M_PI / 180.0);

    // Using Haversine Distance: https://en.wikipedia.org/wiki/Haversine_formula
    // example in JS: https://github.com/njj/haversine/blob/develop/haversine.js
    double d_lat = lat_2 - lat_1;
    double d_lon = lon_2 - lon_1;
    double a = pow(sin(d_lat / 2.0), 2) + pow(sin(d_lon / 2.0), 2) * cos(lat_1) * cos(lat_2);
    double c = 2.0 * atan2(sqrt(a), sqrt(1 - a));

    return 6371000 * c;  // return in meters
}

void Graph::createCH(Heuristic heuristic) {
    std::cout << "Started creating CH." << std::endl;

    auto begin = std::chrono::high_resolution_clock::now();

    if (m_graph.empty()) {
        std::cout << "Can't create CH, because graph is empty." << std::endl;
        return;
    }
    m_node_level.clear();
    m_node_level.resize(m_num_nodes);

    // reverse graph needed to find incoming edges
    m_ch_available = false;
    createReverseGraph(false);
    m_ch_available = true;

    std::vector<bool> contracted(m_num_nodes, false);
    int num_contracted = 0;

    m_contr_data = ContractionData(m_num_nodes);
    m_contr_data.m_num_contracted_neighbours.clear();
    m_contr_data.m_num_contracted_neighbours.resize(m_num_nodes, 0);

    // std::vector<int> longest_path_fwd(m_num_nodes, 0);
    // std::vector<int> longest_path_bwd(m_num_nodes, 0);

    // initialize importance for all nodes
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>>
        importance_pq;
    for (int i = 0; i < m_num_nodes; ++i) {
        int importance = 0;
        switch (heuristic) {
            case Heuristic::IN_OUT:
                importance = inOutProductHeuristic(contracted, i);
                break;
            case Heuristic::EDGE_DIFFERENCE:
                importance = edgeDifferenceHeuristic(contracted, i);
                break;
            case Heuristic::WEIGHTED_COST:
                importance = weightedCostHeuristic(contracted, i);
                break;
            case Heuristic::MICROSOFT:
                importance = microsoftHeuristic(contracted, i, num_contracted);
                break;
            default:
                break;
        }
        importance_pq.emplace(std::make_pair(importance, i));
    }

    while (num_contracted != m_num_nodes) {
        auto contracted_node = importance_pq.top();
        importance_pq.pop();

        int new_importance = 0;
        switch (heuristic) {
            case Heuristic::IN_OUT:
                new_importance = inOutProductHeuristic(contracted, contracted_node.second);
                break;
            case Heuristic::EDGE_DIFFERENCE:
                new_importance = edgeDifferenceHeuristic(contracted, contracted_node.second);
                break;
            case Heuristic::WEIGHTED_COST:
                new_importance = weightedCostHeuristic(contracted, contracted_node.second);
                break;
            case Heuristic::MICROSOFT:
                new_importance = microsoftHeuristic(contracted, contracted_node.second, num_contracted);
                break;
            default:
                break;
        }

        while (new_importance > importance_pq.top().first) {
            importance_pq.emplace(std::make_pair(new_importance, contracted_node.second));
            contracted_node = importance_pq.top();
            importance_pq.pop();

            switch (heuristic) {
                case Heuristic::IN_OUT:
                    new_importance = inOutProductHeuristic(contracted, contracted_node.second);
                    break;
                case Heuristic::EDGE_DIFFERENCE:
                    new_importance = edgeDifferenceHeuristic(contracted, contracted_node.second);
                    break;
                case Heuristic::WEIGHTED_COST:
                    new_importance = weightedCostHeuristic(contracted, contracted_node.second);
                    break;
                case Heuristic::MICROSOFT:
                    new_importance = microsoftHeuristic(contracted, contracted_node.second, num_contracted);
                    break;
                default:
                    break;
            }
        }

        contractNode(contracted, contracted_node.second);
        contracted[contracted_node.second] = true;
        m_node_level[contracted_node.second] = num_contracted;
        ++num_contracted;

        // std::cout << "Finished contracting: " << num_contracted << "\n";
    }

    m_contr_data.m_distances.clear();
    m_contr_data.m_num_contracted_neighbours.clear();
    m_contr_data.m_outgoing.clear();
    m_contr_data.m_visited.clear();

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - begin);
    std::cout << "Finished creating CH. Took " << elapsed.count() << " seconds" << std::endl;
}

int Graph::inOutProductHeuristic(std::vector<bool>& contracted, int node) {
    int num_outgoing = 0;
    int num_incomming = 0;

    for (auto& outgoing : m_graph[node]) {
        if (!contracted[outgoing.m_target]) ++num_outgoing;
    }

    for (auto& incoming : m_reverse_graph[node]) {
        if (!contracted[incoming.m_target]) ++num_incomming;
    }
    return num_outgoing * num_incomming;
}

int Graph::edgeDifferenceHeuristic(std::vector<bool>& contracted, int node) {
    int num_edges_deleted = 0;
    for (auto& outgoing : m_graph[node]) {
        if (!contracted[outgoing.m_target]) ++num_edges_deleted;
    }

    for (auto& incoming : m_reverse_graph[node]) {
        if (!contracted[incoming.m_target]) ++num_edges_deleted;
    }

    int num_added_shortcuts = 0;

    // mark outgoing nodes
    // will be used for pruning in contractionDijkstra
    int num_outgoing = 0;
    int max_distance_out = -1;
    for (auto& outgoing : m_graph[node]) {
        if (contracted[outgoing.m_target]) continue;
        if (outgoing.m_cost > max_distance_out) max_distance_out = outgoing.m_cost;

        ++num_outgoing;
        m_contr_data.m_reset_outgoing.push_back(outgoing.m_target);
        m_contr_data.m_outgoing[outgoing.m_target] = true;
    }

    for (auto& incoming : m_reverse_graph[node]) {
        if (contracted[incoming.m_target]) continue;
        int max_distance = incoming.m_cost + max_distance_out;

        contractionDijkstra(incoming.m_target, node, contracted, num_outgoing, max_distance);

        for (auto& outgoing : m_graph[node]) {
            if (m_contr_data.m_visited[outgoing.m_target]) continue;
            if (contracted[outgoing.m_target]) continue;
            if (m_contr_data.m_distances[outgoing.m_target] != incoming.m_cost + outgoing.m_cost) continue;
            if (outgoing.m_target == incoming.m_target) continue;

            m_contr_data.m_visited[outgoing.m_target] = true;
            m_contr_data.m_reset_visited.push_back(outgoing.m_target);

            ++num_added_shortcuts;
        }

        // reset contraction data
        for (int& num : m_contr_data.m_reset_visited) m_contr_data.m_visited[num] = false;
        for (int& num : m_contr_data.m_reset_distances) m_contr_data.m_distances[num] = std::numeric_limits<int>::max();
        m_contr_data.m_reset_visited.clear();
        m_contr_data.m_reset_distances.clear();
    }

    // reset contraction data
    for (int& num : m_contr_data.m_reset_outgoing) m_contr_data.m_outgoing[num] = false;
    m_contr_data.m_reset_outgoing.clear();

    return num_added_shortcuts - num_edges_deleted;
}

int Graph::weightedCostHeuristic(std::vector<bool>& contracted, int node) {
    int num_outgoing = 0;
    int num_incomming = 0;

    for (auto& outgoing : m_graph[node]) {
        if (!contracted[outgoing.m_target]) ++num_outgoing;
    }

    for (auto& incoming : m_reverse_graph[node]) {
        if (!contracted[incoming.m_target]) ++num_incomming;
    }

    int max_cost = 0;
    int num_added_shortcuts = 0;

    // mark outgoing nodes
    // will be used for pruning in contractionDijkstra
    int max_distance_out = -1;
    for (auto& outgoing : m_graph[node]) {
        if (contracted[outgoing.m_target]) continue;
        if (outgoing.m_cost > max_distance_out) max_distance_out = outgoing.m_cost;

        m_contr_data.m_reset_outgoing.push_back(outgoing.m_target);
        m_contr_data.m_outgoing[outgoing.m_target] = true;
    }

    for (auto& incoming : m_reverse_graph[node]) {
        if (contracted[incoming.m_target]) continue;
        int max_distance = incoming.m_cost + max_distance_out;

        contractionDijkstra(incoming.m_target, node, contracted, num_outgoing, max_distance);

        for (auto& outgoing : m_graph[node]) {
            if (m_contr_data.m_visited[outgoing.m_target]) continue;
            if (contracted[outgoing.m_target]) continue;
            if (m_contr_data.m_distances[outgoing.m_target] != incoming.m_cost + outgoing.m_cost) continue;
            if (outgoing.m_target == incoming.m_target) continue;

            m_contr_data.m_visited[outgoing.m_target] = true;
            m_contr_data.m_reset_visited.push_back(outgoing.m_target);

            ++num_added_shortcuts;

            if (incoming.m_cost + outgoing.m_cost > max_cost) max_cost = incoming.m_cost + outgoing.m_cost;
        }

        // reset contraction data
        for (int& num : m_contr_data.m_reset_visited) m_contr_data.m_visited[num] = false;
        for (int& num : m_contr_data.m_reset_distances) m_contr_data.m_distances[num] = std::numeric_limits<int>::max();
        m_contr_data.m_reset_visited.clear();
        m_contr_data.m_reset_distances.clear();
    }

    // reset contraction data
    for (int& num : m_contr_data.m_reset_outgoing) m_contr_data.m_outgoing[num] = false;
    m_contr_data.m_reset_outgoing.clear();

    return 0.8 * max_cost + 0.2 * num_outgoing * num_incomming;
    // return max_cost;
}

int Graph::altWeightedCostHeuristic(std::vector<bool>& contracted, int node, std::vector<int>& longest_path_fwd,
                                    std::vector<int>& longest_path_bwd) {
    // TODO:
    int num_outgoing = 0;
    int num_incomming = 0;

    for (auto& outgoing : m_graph[node]) {
        if (!contracted[outgoing.m_target]) ++num_outgoing;
    }

    for (auto& incoming : m_reverse_graph[node]) {
        if (!contracted[incoming.m_target]) ++num_incomming;
    }

    int max_cost = 0;
    int max_path = 0;

    // mark outgoing nodes
    // will be used for pruning in contractionDijkstra
    int max_distance_out = -1;
    for (auto& outgoing : m_graph[node]) {
        if (contracted[outgoing.m_target]) continue;
        if (outgoing.m_cost > max_distance_out) max_distance_out = outgoing.m_cost;

        m_contr_data.m_reset_outgoing.push_back(outgoing.m_target);
        m_contr_data.m_outgoing[outgoing.m_target] = true;
    }

    for (auto& incoming : m_reverse_graph[node]) {
        if (contracted[incoming.m_target]) continue;
        int max_distance = incoming.m_cost + max_distance_out;

        contractionDijkstra(incoming.m_target, node, contracted, num_outgoing, max_distance);

        for (auto& outgoing : m_graph[node]) {
            if (contracted[outgoing.m_target]) continue;
            if (outgoing.m_target == incoming.m_target) continue;
            if (m_contr_data.m_visited[outgoing.m_target]) continue;
            if (m_contr_data.m_distances[outgoing.m_target] != incoming.m_cost + outgoing.m_cost) continue;

            // TODO: Das auch noch für incoming machen in einer extra liste
            // update longest_path
            if (longest_path_fwd[outgoing.m_target] < longest_path_fwd[node] + outgoing.m_cost)
                longest_path_fwd[outgoing.m_target] = longest_path_fwd[node] + outgoing.m_cost;
            if (longest_path_bwd[outgoing.m_target] < longest_path_bwd[node] + incoming.m_cost)
                longest_path_bwd[outgoing.m_target] = longest_path_bwd[node] + incoming.m_cost;

            m_contr_data.m_visited[outgoing.m_target] = true;
            m_contr_data.m_reset_visited.push_back(outgoing.m_target);

            if (longest_path_fwd[outgoing.m_target] > max_path) max_path = longest_path_fwd[node];
            if (longest_path_bwd[outgoing.m_target] > max_path) max_path = longest_path_bwd[node];
            if (outgoing.m_cost + incoming.m_cost > max_cost) max_cost = outgoing.m_cost + incoming.m_cost;
        }

        // reset contraction data
        for (int& num : m_contr_data.m_reset_visited) m_contr_data.m_visited[num] = false;
        for (int& num : m_contr_data.m_reset_distances) m_contr_data.m_distances[num] = std::numeric_limits<int>::max();
        m_contr_data.m_reset_visited.clear();
        m_contr_data.m_reset_distances.clear();
    }

    // reset contraction data
    for (int& num : m_contr_data.m_reset_outgoing) m_contr_data.m_outgoing[num] = false;
    m_contr_data.m_reset_outgoing.clear();

    return 0.799 * max_cost + 0.2 * num_outgoing * num_incomming + 0.001 * max_path;
}

int Graph::deletedNeighboursHeuristic(std::vector<bool>& contracted, int node,
                                      std::vector<int>& num_deleted_neighbours) {
    int num_outgoing = 0;
    int num_incomming = 0;

    for (auto& outgoing : m_graph[node]) {
        if (!contracted[outgoing.m_target]) ++num_outgoing;
    }

    for (auto& incoming : m_reverse_graph[node]) {
        if (!contracted[incoming.m_target]) ++num_incomming;
    }

    int max_cost = 0;

    // mark outgoing nodes
    // will be used for pruning in contractionDijkstra
    int max_distance_out = -1;
    for (auto& outgoing : m_graph[node]) {
        if (contracted[outgoing.m_target]) continue;
        if (outgoing.m_cost > max_distance_out) max_distance_out = outgoing.m_cost;

        m_contr_data.m_reset_outgoing.push_back(outgoing.m_target);
        m_contr_data.m_outgoing[outgoing.m_target] = true;
    }

    for (auto& incoming : m_reverse_graph[node]) {
        if (contracted[incoming.m_target]) continue;
        int max_distance = incoming.m_cost + max_distance_out;
        contractionDijkstra(incoming.m_target, node, contracted, num_outgoing, max_distance);

        for (auto& outgoing : m_graph[node]) {
            if (m_contr_data.m_visited[outgoing.m_target]) continue;
            if (contracted[outgoing.m_target]) continue;
            if (m_contr_data.m_distances[outgoing.m_target] != incoming.m_cost + outgoing.m_cost) continue;
            if (outgoing.m_target == incoming.m_target) continue;

            m_contr_data.m_visited[outgoing.m_target] = true;
            m_contr_data.m_reset_visited.push_back(outgoing.m_target);

            if (incoming.m_cost + outgoing.m_cost > max_cost) max_cost = incoming.m_cost + outgoing.m_cost;
        }

        // reset contraction data
        for (int& num : m_contr_data.m_reset_visited) m_contr_data.m_visited[num] = false;
        for (int& num : m_contr_data.m_reset_distances) m_contr_data.m_distances[num] = std::numeric_limits<int>::max();
        m_contr_data.m_reset_visited.clear();
        m_contr_data.m_reset_distances.clear();
    }

    // reset contraction data
    for (int& num : m_contr_data.m_reset_outgoing) m_contr_data.m_outgoing[num] = false;
    m_contr_data.m_reset_outgoing.clear();

    return 0.8 * max_cost + 0.2 * num_outgoing * num_incomming + 0.0 * num_deleted_neighbours[node];
}

int Graph::microsoftHeuristic(std::vector<bool>& contracted, int node, int cur_level) {
    int num_outgoing = 0;
    int num_incomming = 0;

    for (auto& outgoing : m_graph[node]) {
        if (!contracted[outgoing.m_target]) ++num_outgoing;
    }

    for (auto& incoming : m_reverse_graph[node]) {
        if (!contracted[incoming.m_target]) ++num_incomming;
    }

    int max_cost = 0;
    int num_added_shortcuts = 0;

    int max_neighbour_level = -1;

    // mark outgoing nodes
    // will be used for pruning in contractionDijkstra
    int max_distance_out = -1;
    for (auto& outgoing : m_graph[node]) {
        if (contracted[outgoing.m_target]) continue;
        if (outgoing.m_cost > max_distance_out) max_distance_out = outgoing.m_cost;

        m_contr_data.m_reset_outgoing.push_back(outgoing.m_target);
        m_contr_data.m_outgoing[outgoing.m_target] = true;
    }

    int underlying_shortcuts = 0;

    for (auto& incoming : m_reverse_graph[node]) {
        if (contracted[incoming.m_target]) {
            if (m_node_level[incoming.m_target] > max_neighbour_level)
                max_neighbour_level = m_node_level[incoming.m_target];
            continue;
        }
        int max_distance = incoming.m_cost + max_distance_out;
        contractionDijkstra(incoming.m_target, node, contracted, num_outgoing, max_distance);

        for (auto& outgoing : m_graph[node]) {
            if (m_contr_data.m_visited[outgoing.m_target]) continue;
            if (contracted[outgoing.m_target]) {
                if (m_node_level[incoming.m_target] > max_neighbour_level)
                    max_neighbour_level = m_node_level[incoming.m_target];
                continue;
            }
            if (m_contr_data.m_distances[outgoing.m_target] != incoming.m_cost + outgoing.m_cost) continue;
            if (outgoing.m_target == incoming.m_target) continue;

            m_contr_data.m_visited[outgoing.m_target] = true;
            m_contr_data.m_reset_visited.push_back(outgoing.m_target);

            if (incoming.m_child_1 != -1) underlying_shortcuts += incoming.m_child_1;
            if (outgoing.m_child_1 != -1) underlying_shortcuts += outgoing.m_child_1;

            ++num_added_shortcuts;

            if (incoming.m_cost + outgoing.m_cost > max_cost) max_cost = incoming.m_cost + outgoing.m_cost;
        }

        // reset contraction data
        for (int& num : m_contr_data.m_reset_visited) m_contr_data.m_visited[num] = false;
        for (int& num : m_contr_data.m_reset_distances) m_contr_data.m_distances[num] = std::numeric_limits<int>::max();
        m_contr_data.m_reset_visited.clear();
        m_contr_data.m_reset_distances.clear();
    }

    // reset contraction data
    for (int& num : m_contr_data.m_reset_outgoing) m_contr_data.m_outgoing[num] = false;
    m_contr_data.m_reset_outgoing.clear();

    if (max_neighbour_level > cur_level)
        max_neighbour_level = 0;
    else
        ++max_neighbour_level;

    return static_cast<int>(0.001 * static_cast<double>(max_cost)) +
           2 * (num_added_shortcuts - (num_incomming + num_outgoing)) +
           1 * m_contr_data.m_num_contracted_neighbours[node] + 5 * max_neighbour_level + underlying_shortcuts;
}

void Graph::contractNode(std::vector<bool>& contracted, int contracted_node) {
    // mark outgoing nodes
    // will be used for pruning in contractionDijkstra
    int num_outgoing = 0;
    int max_distance_out = -1;
    for (auto& outgoing : m_graph[contracted_node]) {
        if (contracted[outgoing.m_target]) continue;
        if (outgoing.m_cost > max_distance_out) max_distance_out = outgoing.m_cost;

        ++m_contr_data.m_num_contracted_neighbours[outgoing.m_target];

        ++num_outgoing;
        m_contr_data.m_reset_outgoing.push_back(outgoing.m_target);
        m_contr_data.m_outgoing[outgoing.m_target] = true;
    }

    std::vector<std::pair<int, Edge>> fwd_shortcuts;
    std::vector<std::pair<int, Edge>> bwd_shortcuts;

    for (auto& incoming : m_reverse_graph[contracted_node]) {
        if (contracted[incoming.m_target]) continue;
        int max_distance = incoming.m_cost + max_distance_out;

        ++m_contr_data.m_num_contracted_neighbours[incoming.m_target];

        contractionDijkstra(incoming.m_target, contracted_node, contracted, num_outgoing, max_distance);

        for (auto& outgoing : m_graph[contracted_node]) {
            if (m_contr_data.m_visited[outgoing.m_target]) continue;
            if (contracted[outgoing.m_target]) continue;
            if (m_contr_data.m_distances[outgoing.m_target] != incoming.m_cost + outgoing.m_cost) continue;
            if (outgoing.m_target == incoming.m_target) continue;

            m_contr_data.m_visited[outgoing.m_target] = true;
            m_contr_data.m_reset_visited.push_back(outgoing.m_target);

            int underlying_arcs = 0;
            if (incoming.m_child_1 != -1) underlying_arcs += incoming.m_child_1;
            if (outgoing.m_child_1 != -1) underlying_arcs += outgoing.m_child_1;

            // add shortcut
            Edge fwd_shortcut(outgoing.m_target, incoming.m_cost + outgoing.m_cost, underlying_arcs, -1);
            fwd_shortcuts.push_back(std::make_pair(incoming.m_target, fwd_shortcut));
            Edge bwd_shortcut(incoming.m_target, incoming.m_cost + outgoing.m_cost, underlying_arcs, -1);
            bwd_shortcuts.push_back(std::make_pair(outgoing.m_target, bwd_shortcut));
        }

        // reset contraction data
        for (int& num : m_contr_data.m_reset_visited) m_contr_data.m_visited[num] = false;
        for (int& num : m_contr_data.m_reset_distances) m_contr_data.m_distances[num] = std::numeric_limits<int>::max();
        m_contr_data.m_reset_visited.clear();
        m_contr_data.m_reset_distances.clear();
    }

    // reset contraction data
    for (int& num : m_contr_data.m_reset_outgoing) m_contr_data.m_outgoing[num] = false;
    m_contr_data.m_reset_outgoing.clear();

    for (auto& fwd_shortcut : fwd_shortcuts) m_graph[fwd_shortcut.first].push_back(fwd_shortcut.second);
    for (auto& bwd_shortcut : bwd_shortcuts) m_reverse_graph[bwd_shortcut.first].push_back(bwd_shortcut.second);
}

void Graph::contractionDijkstra(int start, int contracted_node, std::vector<bool>& contracted, int num_outgoing,
                                int max_distance) {
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;

    int num_visited_outgoing = 0;
    m_contr_data.m_distances[start] = 0;
    m_contr_data.m_reset_distances.push_back(start);
    pq.push(std::make_pair(0, start));

    while (!pq.empty()) {
        std::pair<int, int> cur_node = pq.top();
        pq.pop();

        if (m_contr_data.m_distances[cur_node.second] != cur_node.first) continue;

        if (m_contr_data.m_outgoing[cur_node.second]) ++num_visited_outgoing;
        if (cur_node.first > max_distance || num_visited_outgoing == num_outgoing) break;

        for (Edge& e : m_graph[cur_node.second]) {
            if (contracted[e.m_target]) continue;
            if (m_contr_data.m_distances[e.m_target] > m_contr_data.m_distances[cur_node.second] + e.m_cost) {
                if (m_contr_data.m_distances[e.m_target] == std::numeric_limits<int>::max())
                    m_contr_data.m_reset_distances.push_back(e.m_target);
                m_contr_data.m_distances[e.m_target] = m_contr_data.m_distances[cur_node.second] + e.m_cost;
                pq.push(std::make_pair(m_contr_data.m_distances[e.m_target], e.m_target));
            }
        }
    }
}

void Graph::lowerBoundDijkstra(LowerBoundData& lb_data) {
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;

    lb_data.m_distances[lb_data.m_start_node] = 0;
    lb_data.m_reset_previous_node.push_back(lb_data.m_start_node);
    pq.push(std::make_pair(0, lb_data.m_start_node));

    while (!pq.empty()) {
        std::pair<int, int> cur_node = pq.top();
        pq.pop();

        if (lb_data.m_distances[cur_node.second] != cur_node.first) continue;

        if (cur_node.first > lb_data.m_threshold * 2) break;

        for (Edge& e : m_graph[cur_node.second]) {
            // if (lb_data.m_marked[e.m_target]) continue;  // TODO: passt das?
            if (lb_data.m_distances[e.m_target] > lb_data.m_distances[cur_node.second] + e.m_cost) {
                if (lb_data.m_distances[e.m_target] == std::numeric_limits<int>::max())
                    lb_data.m_reset_previous_node.push_back(e.m_target);
                lb_data.m_distances[e.m_target] = lb_data.m_distances[cur_node.second] + e.m_cost;
                lb_data.m_previous_node[e.m_target] = cur_node.second;
                pq.push(std::make_pair(lb_data.m_distances[e.m_target], e.m_target));
            }
        }
    }
}

bool Graph::pathCoverVerificationDijkstra(LowerBoundData& lb_data) {
    // TODO:
    std::priority_queue<std::pair<int, std::pair<int, bool>>, std::vector<std::pair<int, std::pair<int, bool>>>,
                        std::greater<std::pair<int, std::pair<int, bool>>>>
        pq;

    lb_data.m_distances[lb_data.m_start_node] = 0;
    lb_data.m_reset_previous_node.push_back(lb_data.m_start_node);
    pq.push(std::make_pair(0, std::make_pair(lb_data.m_start_node, lb_data.m_marked[lb_data.m_start_node])));

    while (!pq.empty()) {
        auto cur_node = pq.top();
        pq.pop();

        if (lb_data.m_distances[cur_node.second.first] != cur_node.first) continue;

        if (lb_data.m_distances[cur_node.second.first] > lb_data.m_threshold && !cur_node.second.second) return false;

        if (static_cast<double>(cur_node.first) > static_cast<double>(lb_data.m_threshold) * 1.1) break;

        for (Edge& e : m_graph[cur_node.second.first]) {
            if (lb_data.m_distances[e.m_target] > lb_data.m_distances[cur_node.second.first] + e.m_cost) {
                if (lb_data.m_distances[e.m_target] == std::numeric_limits<int>::max())
                    lb_data.m_reset_previous_node.push_back(e.m_target);
                lb_data.m_distances[e.m_target] = lb_data.m_distances[cur_node.second.first] + e.m_cost;

                bool covered;
                if (lb_data.m_marked[e.m_target])
                    covered = true;
                else
                    covered = cur_node.second.second;
                pq.push(std::make_pair(lb_data.m_distances[e.m_target], std::make_pair(e.m_target, covered)));
            }
        }
    }

    return true;
}

int Graph::simplifiedHubLabelQuery(std::vector<std::pair<int, int>>& fwd_labels, int node) {
    int distance = std::numeric_limits<int>::max();
    auto fwd_iter = fwd_labels.begin();
    int bwd_node_index = m_bwd_indices[m_node_indices[node]];
    int bwd_next_index = m_bwd_indices[m_node_indices[node] + 1];

    while (fwd_iter != fwd_labels.end() && bwd_node_index < bwd_next_index) {
        if (fwd_iter->first == m_bwd_hub_labels[bwd_node_index].first) {
            if (fwd_iter->second + m_bwd_hub_labels[bwd_node_index].second < distance)
                distance = fwd_iter->second + m_bwd_hub_labels[bwd_node_index].second;
            ++fwd_iter;
            ++bwd_node_index;
        } else if (fwd_iter->first < m_bwd_hub_labels[bwd_node_index].first) {
            ++fwd_iter;
        } else {
            ++bwd_node_index;
        }
    }

    return distance;
}

int Graph::simplifiedHubLabelQuery(int node, std::vector<std::pair<int, int>>& bwd_labels) {
    int distance = std::numeric_limits<int>::max();
    int fwd_node_index = m_fwd_indices[m_node_indices[node]];
    int fwd_next_index = m_fwd_indices[m_node_indices[node] + 1];
    auto bwd_iter = bwd_labels.begin();

    while (fwd_node_index < fwd_next_index && bwd_iter != bwd_labels.end()) {
        if (bwd_iter->first == m_fwd_hub_labels[fwd_node_index].first) {
            if (bwd_iter->second + m_fwd_hub_labels[fwd_node_index].second < distance)
                distance = bwd_iter->second + m_fwd_hub_labels[fwd_node_index].second;
            ++bwd_iter;
            ++fwd_node_index;
        } else if (bwd_iter->first < m_fwd_hub_labels[fwd_node_index].first) {
            ++bwd_iter;
        } else {
            ++fwd_node_index;
        }
    }

    return distance;
}

}  // namespace olsp
