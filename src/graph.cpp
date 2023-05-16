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

Graph::Graph(const std::string& path, ReadMode read_mode, bool ch_available, bool prune_graph, int num_threads,
             DistanceMode dist_mode)
    : m_ch_available(ch_available), m_num_threads(num_threads) {
    readGraph(path, read_mode, dist_mode);
    m_num_nodes = m_graph.size();

    if (read_mode == ReadMode::NORMAL && m_ch_available) {
        createCH();  // TODO:
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
    std::unordered_set<int> different_levels;
    for (int i = 0; i < m_num_nodes; i++) {
        different_levels.emplace(m_node_level[i]);
    }
    std::vector<int> different_levels_vec;
    different_levels_vec.insert(different_levels_vec.end(), different_levels.begin(), different_levels.end());
    std::sort(different_levels_vec.begin(), different_levels_vec.end(),
              [](int& left, int& right) { return left > right; });
    std::vector<std::vector<int>> level_buckets(different_levels_vec.size());
    for (int i = 0; i < m_num_nodes; ++i) {
        auto iter = find(different_levels_vec.begin(), different_levels_vec.end(), m_node_level[i]);
        level_buckets[iter - different_levels_vec.begin()].push_back(i);
    }

    int num_calculated = 0;

    for (int i = 0; i < different_levels_vec.size(); ++i) {
#pragma omp parallel for
        for (int j = 0; j < level_buckets[i].size(); ++j) {
            int node = level_buckets[i][j];
            m_fwd_hub_labels[node].push_back(std::make_pair(node, 0));
            m_bwd_hub_labels[node].push_back(std::make_pair(node, 0));

            if (i == 0) continue;

            // fwd lables
            for (Edge& e : m_graph[node]) {
                if (m_node_level[node] >= m_node_level[e.m_target]) continue;

                for (std::pair<int, int>& hub : m_fwd_hub_labels[e.m_target])
                    m_fwd_hub_labels[node].push_back(std::make_pair(hub.first, hub.second + e.m_cost));
            }

            // remove duplicates
            std::sort(m_fwd_hub_labels[node].begin(), m_fwd_hub_labels[node].end(),
                      [](auto& left, auto& right) { return left.first < right.first; });
            for (auto iter = m_fwd_hub_labels[node].begin(); iter != m_fwd_hub_labels[node].end();) {
                auto iter_2 = m_fwd_hub_labels[node].end() - 1;
                if (std::distance(iter, iter_2) != 0 && iter->first == (iter + 1)->first) {
                    if (iter->second >= (iter + 1)->second) {
                        iter = m_fwd_hub_labels[node].erase(iter);
                    } else {
                        iter = m_fwd_hub_labels[node].erase(iter + 1);
                        --iter;
                    }
                } else {
                    ++iter;
                }
            }

            for (auto iter = m_fwd_hub_labels[node].begin(); iter != m_fwd_hub_labels[node].end();) {
                int best_dist = std::numeric_limits<int>::max();
                if (iter->first != node)
                    best_dist = simplifiedHubLabelQuery(m_fwd_hub_labels[node], m_bwd_hub_labels[iter->first]);
                if (best_dist < iter->second)
                    iter = m_fwd_hub_labels[node].erase(iter);
                else
                    ++iter;
            }

            // bwd labels
            for (Edge& e : m_reverse_graph[node]) {
                if (m_node_level[node] >= m_node_level[e.m_target]) continue;

                for (std::pair<int, int>& hub : m_bwd_hub_labels[e.m_target])
                    m_bwd_hub_labels[node].push_back(std::make_pair(hub.first, hub.second + e.m_cost));
            }

            // remove duplicates
            std::sort(m_bwd_hub_labels[node].begin(), m_bwd_hub_labels[node].end(),
                      [](auto& left, auto& right) { return left.first < right.first; });
            for (auto iter = m_bwd_hub_labels[node].begin(); iter != m_bwd_hub_labels[node].end();) {
                auto iter_2 = m_bwd_hub_labels[node].end() - 1;
                if (std::distance(iter, iter_2) != 0 && iter->first == (iter + 1)->first) {
                    if (iter->second >= (iter + 1)->second) {
                        iter = m_bwd_hub_labels[node].erase(iter);
                    } else {
                        iter = m_bwd_hub_labels[node].erase(iter + 1);
                        --iter;
                    }
                } else {
                    ++iter;
                }
            }

            for (auto iter = m_bwd_hub_labels[node].begin(); iter != m_bwd_hub_labels[node].end();) {
                int best_dist = std::numeric_limits<int>::max();
                if (iter->first != node)
                    best_dist = simplifiedHubLabelQuery(m_fwd_hub_labels[iter->first], m_bwd_hub_labels[node]);
                if (best_dist < iter->second)
                    iter = m_bwd_hub_labels[node].erase(iter);
                else
                    ++iter;
            }
        }
        num_calculated += level_buckets[i].size();
        std::cout << "Finished Hub Labels for " << num_calculated << " num of nodes." << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Finished creating hub labels. Took " << elapsed.count() << " milliseconds " << std::endl;
}

void Graph::advancedCreateHubLabels() {
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

    std::vector<AdvancedHubLabelData> hub_label_data(m_num_threads);
    for (int i = 0; i < m_num_threads; ++i) {
        hub_label_data[i] = AdvancedHubLabelData(m_num_nodes);
    }

    // sort the levels, but don't change the original vector
    std::unordered_set<int> different_levels;
    for (int i = 0; i < m_num_nodes; i++) {
        different_levels.emplace(m_node_level[i]);
    }
    std::vector<int> different_levels_vec;
    different_levels_vec.insert(different_levels_vec.end(), different_levels.begin(), different_levels.end());
    std::sort(different_levels_vec.begin(), different_levels_vec.end(),
              [](int& left, int& right) { return left > right; });
    std::vector<std::vector<int>> level_buckets(different_levels_vec.size());
    for (int i = 0; i < m_num_nodes; ++i) {
        auto iter = find(different_levels_vec.begin(), different_levels_vec.end(), m_node_level[i]);
        level_buckets[iter - different_levels_vec.begin()].push_back(i);
    }

    int num_calculated = 0;

    for (int i = 0; i < different_levels_vec.size(); ++i) {
#pragma omp parallel for
        for (int k = 0; k < level_buckets[i].size(); ++k) {
            int node = level_buckets[i][k];
            int thread_num = omp_get_thread_num();

            m_fwd_hub_labels[node].push_back(std::make_pair(node, 0));
            m_bwd_hub_labels[node].push_back(std::make_pair(node, 0));

            if (i == 0) continue;

            // fwd lables
            // do forward search and pruning
            forwardCHSearch(hub_label_data[thread_num], node);
            for (int j = 0; j < hub_label_data[thread_num].m_reset_nodes_fwd.size(); ++j) {
                if (hub_label_data[thread_num].m_reset_nodes_fwd[j] == node) continue;
                // if (m_node_level[node] >= m_node_level[hub_label_data[thread_num].m_reset_nodes_fwd[j]]) continue;
                bool should_be_added = true;
                for (Edge& e : m_reverse_graph[hub_label_data[thread_num].m_reset_nodes_fwd[j]]) {
                    // if (m_node_level[e.m_target] <= m_node_level[hub_label_data[thread_num].m_reset_nodes_fwd[j]])
                    // continue;

                    if (hub_label_data[thread_num].m_distances_fwd[e.m_target] + e.m_cost <
                        hub_label_data[thread_num].m_distances_fwd[hub_label_data[thread_num].m_reset_nodes_fwd[j]]) {
                        should_be_added = false;
                        break;
                    }
                }

                if (should_be_added)
                    m_fwd_hub_labels[node].push_back(std::make_pair(
                        hub_label_data[thread_num].m_reset_nodes_fwd[j],
                        hub_label_data[thread_num].m_distances_fwd[hub_label_data[thread_num].m_reset_nodes_fwd[j]]));
            }

            std::sort(m_fwd_hub_labels[node].begin(), m_fwd_hub_labels[node].end(),
                      [](auto& left, auto& right) { return left.first < right.first; });

            // bootstrapping
            for (auto iter = m_fwd_hub_labels[node].begin(); iter != m_fwd_hub_labels[node].end();) {
                int best_dist = std::numeric_limits<int>::max();
                if (iter->first != node)
                    best_dist = simplifiedHubLabelQuery(m_fwd_hub_labels[node], m_bwd_hub_labels[iter->first]);
                if (best_dist < iter->second)
                    iter = m_fwd_hub_labels[node].erase(iter);
                else
                    ++iter;
            }

            // bwd labels
            // do backward search and pruning
            backwardCHSearch(hub_label_data[thread_num], node);
            for (int j = 0; j < hub_label_data[thread_num].m_reset_nodes_bwd.size(); ++j) {
                if (hub_label_data[thread_num].m_reset_nodes_bwd[j] == node) continue;
                // if (m_node_level[node] >= m_node_level[hub_label_data[thread_num].m_reset_nodes_bwd[j]]) continue;

                bool should_be_added = true;
                for (Edge& e : m_graph[hub_label_data[thread_num].m_reset_nodes_bwd[j]]) {
                    // if (m_node_level[e.m_target] <= m_node_level[hub_label_data[thread_num].m_reset_nodes_bwd[j]])
                    // continue;

                    if (hub_label_data[thread_num].m_distances_bwd[e.m_target] + e.m_cost <
                        hub_label_data[thread_num].m_distances_bwd[hub_label_data[thread_num].m_reset_nodes_bwd[j]]) {
                        should_be_added = false;
                        break;
                    }
                }

                if (should_be_added)
                    m_bwd_hub_labels[node].push_back(std::make_pair(
                        hub_label_data[thread_num].m_reset_nodes_bwd[j],
                        hub_label_data[thread_num].m_distances_bwd[hub_label_data[thread_num].m_reset_nodes_bwd[j]]));
            }

            std::sort(m_bwd_hub_labels[node].begin(), m_bwd_hub_labels[node].end(),
                      [](auto& left, auto& right) { return left.first < right.first; });

            // bootstrapping
            for (auto iter = m_bwd_hub_labels[node].begin(); iter != m_bwd_hub_labels[node].end();) {
                int best_dist = std::numeric_limits<int>::max();
                if (iter->first != node)
                    best_dist = simplifiedHubLabelQuery(m_fwd_hub_labels[iter->first], m_bwd_hub_labels[node]);
                if (best_dist < iter->second)
                    iter = m_bwd_hub_labels[node].erase(iter);
                else
                    ++iter;
            }

            // reset data for next run
            for (int& index : hub_label_data[thread_num].m_reset_nodes_fwd) {
                hub_label_data[thread_num].m_distances_fwd[index] = std::numeric_limits<int>::max();
                hub_label_data[thread_num].m_visited_fwd[index] = false;
            }
            hub_label_data[thread_num].m_reset_nodes_fwd.clear();

            for (int& index : hub_label_data[thread_num].m_reset_nodes_bwd) {
                hub_label_data[thread_num].m_distances_bwd[index] = std::numeric_limits<int>::max();
                hub_label_data[thread_num].m_visited_bwd[index] = false;
            }
            hub_label_data[thread_num].m_reset_nodes_bwd.clear();
        }

        num_calculated += level_buckets[i].size();
        std::cout << "Finished Hub Labels for " << num_calculated << " num of nodes." << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Finished creating hub labels. Took " << elapsed.count() << " milliseconds " << std::endl;
}

void Graph::forwardCHSearch(AdvancedHubLabelData& data, int start_node) {
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>>
        fwd_pq;

    data.m_distances_fwd[start_node] = 0;
    data.m_visited_fwd[start_node] = true;
    data.m_reset_nodes_fwd.push_back(start_node);

    // first corresponds to distance and second is node index
    fwd_pq.push(std::make_pair(0, start_node));

    std::pair<int, int> fwd_node;

    while (!fwd_pq.empty()) {
        fwd_node = fwd_pq.top();
        fwd_pq.pop();
        if (data.m_visited_fwd[fwd_node.second] && fwd_node.first > data.m_distances_fwd[fwd_node.second]) continue;

        for (Edge& e : m_graph[fwd_node.second]) {
            // relax edge
            // if (m_node_level[e.m_target] <= m_node_level[fwd_node.second]) continue;

            if (!data.m_visited_fwd[e.m_target] || data.m_distances_fwd[e.m_target] > fwd_node.first + e.m_cost) {
                if (!data.m_visited_fwd[e.m_target]) data.m_reset_nodes_fwd.push_back(e.m_target);
                data.m_visited_fwd[e.m_target] = true;

                data.m_distances_fwd[e.m_target] = fwd_node.first + e.m_cost;
                fwd_pq.push(std::make_pair(data.m_distances_fwd[e.m_target], e.m_target));
            }
        }
    }
}

void Graph::backwardCHSearch(AdvancedHubLabelData& data, int start_node) {
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>>
        bwd_pq;

    data.m_distances_bwd[start_node] = 0;
    data.m_visited_bwd[start_node] = true;
    data.m_reset_nodes_bwd.push_back(start_node);

    // first corresponds to distance and second is node index
    bwd_pq.push(std::make_pair(0, start_node));

    std::pair<int, int> bwd_node;

    while (!bwd_pq.empty()) {
        bwd_node = bwd_pq.top();
        bwd_pq.pop();
        if (!data.m_visited_bwd[bwd_node.second] && bwd_node.first > data.m_distances_bwd[bwd_node.second]) continue;

        for (Edge& e : m_reverse_graph[bwd_node.second]) {
            // relax edge
            // if (m_node_level[e.m_target] <= m_node_level[bwd_node.second]) continue;

            if (!data.m_visited_bwd[e.m_target] || data.m_distances_bwd[e.m_target] > bwd_node.first + e.m_cost) {
                if (!data.m_visited_bwd[e.m_target]) data.m_reset_nodes_bwd.push_back(e.m_target);
                data.m_visited_bwd[e.m_target] = true;

                data.m_distances_bwd[e.m_target] = bwd_node.first + e.m_cost;
                bwd_pq.push(std::make_pair(data.m_distances_bwd[e.m_target], e.m_target));
            }
        }
    }
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
    auto fwd_iter = m_fwd_hub_labels[data.m_start].begin();
    auto bwd_iter = m_bwd_hub_labels[data.m_end].begin();

    while (fwd_iter != m_fwd_hub_labels[data.m_start].end() && bwd_iter != m_bwd_hub_labels[data.m_end].end()) {
        if (fwd_iter->first == bwd_iter->first) {
            if (fwd_iter->second + bwd_iter->second < data.m_distance) {
                data.m_meeting_node = fwd_iter->first;
                data.m_distance = fwd_iter->second + bwd_iter->second;
            }
            ++fwd_iter;
            ++bwd_iter;
        } else if (fwd_iter->first < bwd_iter->first) {
            ++fwd_iter;
        } else {
            ++bwd_iter;
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

int Graph::numHubLabelsInRange(int lower, int upper) {
    int count = 0;

    for (auto& fwd_labels : m_fwd_hub_labels) {
        for (auto& fwd_label : fwd_labels) {
            if (fwd_label.second > lower && fwd_label.second < upper) ++count;
        }
    }
    for (auto& bwd_labels : m_bwd_hub_labels) {
        for (auto& bwd_label : bwd_labels) {
            if (bwd_label.second > lower && bwd_label.second < upper) ++count;
        }
    }

    return count;
}

std::vector<int> Graph::createShortestPathCover(int threshold) {
    std::cout << "Started creating Path Cover." << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    std::unordered_set<int> path_cover_set;
    for (auto& fwd_labels : m_fwd_hub_labels) {
        for (auto& fwd_label : fwd_labels) {
            if (fwd_label.second > threshold / 2 && fwd_label.second < threshold)
                path_cover_set.emplace(fwd_label.first);
        }
    }
    for (auto& bwd_labels : m_bwd_hub_labels) {
        for (auto& bwd_label : bwd_labels) {
            if (bwd_label.second > threshold / 2 && bwd_label.second < threshold)
                path_cover_set.emplace(bwd_label.first);
        }
    }

    // convert to vector
    std::vector<int> path_cover_vec;
    path_cover_vec.insert(path_cover_vec.end(), path_cover_set.begin(), path_cover_set.end());

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "Finished creating Path Cover. Took " << elapsed.count() << " milliseconds" << std::endl;

    return path_cover_vec;
}

std::vector<int> Graph::lowerBound(std::vector<int>& shortest_path_cover, int threshold) {
    std::cout << "Started creating Path Cover." << std::endl;
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
    while (spc_nodes_covered != shortest_path_cover.size()) {
        int cur_spc_node = dist(rng);
        while (spc_node_covered[cur_spc_node]) cur_spc_node = dist(rng);

        lb_data.m_start_node = cur_spc_node;

        // run dijkstra until threshold reached
        lowerBoundDijkstra(lb_data);

        // look at paths with length greater than threshold
        // 1. first find all end nodes with path longer than threshold
        std::vector<int> end_nodes;
        for (int i = 0; i < lb_data.m_distances.size(); ++i) {
            if (lb_data.m_distances[i] < std::numeric_limits<int>::max() && lb_data.m_distances[i] >= threshold)
                end_nodes.push_back(i);
            lb_data.m_distances[i] = std::numeric_limits<int>::max();
        }

        // 2. backtrack paths through end node
        int selected_end_node = -1;
        int lowest_num_spc_nodes = std::numeric_limits<int>::max();
        for (int& node : end_nodes) {
            int cur_num_spc_nodes = 0;
            int cur_node = node;
            while (cur_node != lb_data.m_start_node) {
                int previous = lb_data.m_previous_node[cur_node];
                if (lb_data.m_marked[previous]) break;

                if (std::find(shortest_path_cover.begin(), shortest_path_cover.end(), cur_node) !=
                    shortest_path_cover.end())
                    ++cur_num_spc_nodes;

                cur_node = previous;
            }

            if (cur_node == lb_data.m_start_node && cur_num_spc_nodes < lowest_num_spc_nodes) {
                lowest_num_spc_nodes = cur_num_spc_nodes;
                selected_end_node = cur_node;
            }
        }

        // 3. go through path of selected end node
        if (selected_end_node != -1) {
            int cur_node = selected_end_node;
            while (cur_node != lb_data.m_start_node) {
                int previous = lb_data.m_previous_node[cur_node];
                lb_data.m_marked[previous] = true;
            }
            lower_bound.push_back(lb_data.m_start_node);
        }

        // reset data for next run
        for (int& node : lb_data.m_reset_previous_node) lb_data.m_previous_node[node] = -1;
        lb_data.m_reset_previous_node.clear();

        ++spc_nodes_covered;
    }

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

void Graph::createCH() {
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
    int cur_level = 0;

    m_contr_data.resize(m_num_threads);
    for (int i = 0; i < m_num_threads; ++i) {
        m_contr_data[i] = ContractionData(m_num_nodes);
        m_contr_data[i].m_num_contracted_neighbours.resize(m_num_nodes);
    }

    // std::vector<int> longest_path_fwd(m_num_nodes, 0);
    // std::vector<int> longest_path_bwd(m_num_nodes, 0);

    // initialize importance for all nodes
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>>
        importance_pq;
    for (int i = 0; i < m_num_nodes; ++i) {
        int importance = microsoftHeuristic(contracted, i, 0);
        importance_pq.emplace(std::make_pair(importance, i));
    }

    while (num_contracted != m_num_nodes) {
        // create independet set
        std::vector<bool> marked(m_num_nodes, false);
        std::vector<int> independent_set;
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>>
            new_pq;
        while (!importance_pq.empty()) {
            auto contracted_node = importance_pq.top();
            importance_pq.pop();

            if (marked[contracted_node.second]) {
                new_pq.emplace(contracted_node);
                continue;
            }

            if (cur_level != 0) {
                int new_importance = microsoftHeuristic(contracted, contracted_node.second, cur_level);
                while (new_importance > importance_pq.top().first) {
                    // std::cout << importance_pq.size() << std::endl;
                    if (marked[contracted_node.second])
                        new_pq.emplace(std::make_pair(new_importance, contracted_node.second));
                    else
                        importance_pq.emplace(std::make_pair(new_importance, contracted_node.second));
                    contracted_node = importance_pq.top();
                    importance_pq.pop();
                    new_importance = microsoftHeuristic(contracted, contracted_node.second, cur_level);
                }
            }

            if (marked[contracted_node.second])
                new_pq.emplace(contracted_node);
            else {
                marked[contracted_node.second] = true;
                independent_set.push_back(contracted_node.second);
                for (Edge& e : m_graph[contracted_node.second]) marked[e.m_target] = true;
                for (Edge& e : m_reverse_graph[contracted_node.second]) marked[e.m_target] = true;
            }
        }
        importance_pq = new_pq;

// contract nodes from independent set in parallel
#pragma omp parallel for
        for (int i = 0; i < independent_set.size(); ++i) {
            contractNode(contracted, independent_set[i]);
            contracted[independent_set[i]] = true;
            m_node_level[independent_set[i]] = cur_level;
        }

        for (int i = 0; i < m_num_threads; ++i) {
            for (auto& fwd_shortcut : m_contr_data[i].m_shortcuts_fwd)
                m_graph[fwd_shortcut.first].push_back(fwd_shortcut.second);
            for (auto& bwd_shortcut : m_contr_data[i].m_shortcuts_bwd)
                m_reverse_graph[bwd_shortcut.first].push_back(bwd_shortcut.second);

            m_contr_data[i].m_shortcuts_fwd.clear();
            m_contr_data[i].m_shortcuts_bwd.clear();

            if (i == 0) continue;

            for (int j = 0; j < m_num_nodes; ++j) {
                m_contr_data[0].m_num_contracted_neighbours[j] += m_contr_data[i].m_num_contracted_neighbours[j];
                m_contr_data[i].m_num_contracted_neighbours[j] = 0;
            }
        }

        num_contracted += static_cast<int>(independent_set.size());

        if (num_contracted >= m_num_nodes * 0.8) {
            omp_set_num_threads(1);
        }

        ++cur_level;
        std::cout << "Finished contracting: " << num_contracted << "\n";
    }

    m_contr_data.clear();
    omp_set_num_threads(m_num_threads);

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
        m_contr_data[0].m_reset_outgoing.push_back(outgoing.m_target);
        m_contr_data[0].m_outgoing[outgoing.m_target] = true;
    }

    for (auto& incoming : m_reverse_graph[node]) {
        if (contracted[incoming.m_target]) continue;
        int max_distance = incoming.m_cost + max_distance_out;

        contractionDijkstra(incoming.m_target, node, contracted, num_outgoing, max_distance);

        for (auto& outgoing : m_graph[node]) {
            if (m_contr_data[0].m_visited[outgoing.m_target]) continue;
            if (contracted[outgoing.m_target]) continue;
            if (m_contr_data[0].m_distances[outgoing.m_target] != incoming.m_cost + outgoing.m_cost) continue;
            if (outgoing.m_target == incoming.m_target) continue;

            m_contr_data[0].m_visited[outgoing.m_target] = true;
            m_contr_data[0].m_reset_visited.push_back(outgoing.m_target);

            ++num_added_shortcuts;
        }

        // reset contraction data
        for (int& num : m_contr_data[0].m_reset_visited) m_contr_data[0].m_visited[num] = false;
        for (int& num : m_contr_data[0].m_reset_distances)
            m_contr_data[0].m_distances[num] = std::numeric_limits<int>::max();
        m_contr_data[0].m_reset_visited.clear();
        m_contr_data[0].m_reset_distances.clear();
    }

    // reset contraction data
    for (int& num : m_contr_data[0].m_reset_outgoing) m_contr_data[0].m_outgoing[num] = false;
    m_contr_data[0].m_reset_outgoing.clear();

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

        m_contr_data.at(0).m_reset_outgoing.push_back(outgoing.m_target);
        m_contr_data[0].m_outgoing[outgoing.m_target] = true;
    }

    for (auto& incoming : m_reverse_graph[node]) {
        if (contracted[incoming.m_target]) continue;
        int max_distance = incoming.m_cost + max_distance_out;

        contractionDijkstra(incoming.m_target, node, contracted, num_outgoing, max_distance);

        for (auto& outgoing : m_graph[node]) {
            if (m_contr_data[0].m_visited[outgoing.m_target]) continue;
            if (contracted[outgoing.m_target]) continue;
            if (m_contr_data[0].m_distances[outgoing.m_target] != incoming.m_cost + outgoing.m_cost) continue;
            if (outgoing.m_target == incoming.m_target) continue;

            m_contr_data[0].m_visited[outgoing.m_target] = true;
            m_contr_data[0].m_reset_visited.push_back(outgoing.m_target);

            ++num_added_shortcuts;

            if (incoming.m_cost + outgoing.m_cost > max_cost) max_cost = incoming.m_cost + outgoing.m_cost;
        }

        // reset contraction data
        for (int& num : m_contr_data[0].m_reset_visited) m_contr_data[0].m_visited[num] = false;
        for (int& num : m_contr_data[0].m_reset_distances)
            m_contr_data[0].m_distances[num] = std::numeric_limits<int>::max();
        m_contr_data[0].m_reset_visited.clear();
        m_contr_data[0].m_reset_distances.clear();
    }

    // reset contraction data
    for (int& num : m_contr_data[0].m_reset_outgoing) m_contr_data[0].m_outgoing[num] = false;
    m_contr_data[0].m_reset_outgoing.clear();

    return 0.6 * max_cost + 0.4 * num_outgoing * num_incomming;
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

        m_contr_data[0].m_reset_outgoing.push_back(outgoing.m_target);
        m_contr_data[0].m_outgoing[outgoing.m_target] = true;
    }

    for (auto& incoming : m_reverse_graph[node]) {
        if (contracted[incoming.m_target]) continue;
        int max_distance = incoming.m_cost + max_distance_out;

        contractionDijkstra(incoming.m_target, node, contracted, num_outgoing, max_distance);

        for (auto& outgoing : m_graph[node]) {
            if (contracted[outgoing.m_target]) continue;
            if (outgoing.m_target == incoming.m_target) continue;
            if (m_contr_data[0].m_visited[outgoing.m_target]) continue;
            if (m_contr_data[0].m_distances[outgoing.m_target] != incoming.m_cost + outgoing.m_cost) continue;

            // TODO: Das auch noch für incoming machen in einer extra liste
            // update longest_path
            if (longest_path_fwd[outgoing.m_target] < longest_path_fwd[node] + outgoing.m_cost)
                longest_path_fwd[outgoing.m_target] = longest_path_fwd[node] + outgoing.m_cost;
            if (longest_path_bwd[outgoing.m_target] < longest_path_bwd[node] + incoming.m_cost)
                longest_path_bwd[outgoing.m_target] = longest_path_bwd[node] + incoming.m_cost;

            m_contr_data[0].m_visited[outgoing.m_target] = true;
            m_contr_data[0].m_reset_visited.push_back(outgoing.m_target);

            if (longest_path_fwd[outgoing.m_target] > max_path) max_path = longest_path_fwd[node];
            if (longest_path_bwd[outgoing.m_target] > max_path) max_path = longest_path_bwd[node];
            if (outgoing.m_cost + incoming.m_cost > max_cost) max_cost = outgoing.m_cost + incoming.m_cost;
        }

        // reset contraction data
        for (int& num : m_contr_data[0].m_reset_visited) m_contr_data[0].m_visited[num] = false;
        for (int& num : m_contr_data[0].m_reset_distances)
            m_contr_data[0].m_distances[num] = std::numeric_limits<int>::max();
        m_contr_data[0].m_reset_visited.clear();
        m_contr_data[0].m_reset_distances.clear();
    }

    // reset contraction data
    for (int& num : m_contr_data[0].m_reset_outgoing) m_contr_data[0].m_outgoing[num] = false;
    m_contr_data[0].m_reset_outgoing.clear();

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

        m_contr_data[0].m_reset_outgoing.push_back(outgoing.m_target);
        m_contr_data[0].m_outgoing[outgoing.m_target] = true;
    }

    for (auto& incoming : m_reverse_graph[node]) {
        if (contracted[incoming.m_target]) continue;
        int max_distance = incoming.m_cost + max_distance_out;
        contractionDijkstra(incoming.m_target, node, contracted, num_outgoing, max_distance);

        for (auto& outgoing : m_graph[node]) {
            if (m_contr_data[0].m_visited[outgoing.m_target]) continue;
            if (contracted[outgoing.m_target]) continue;
            if (m_contr_data[0].m_distances[outgoing.m_target] != incoming.m_cost + outgoing.m_cost) continue;
            if (outgoing.m_target == incoming.m_target) continue;

            m_contr_data[0].m_visited[outgoing.m_target] = true;
            m_contr_data[0].m_reset_visited.push_back(outgoing.m_target);

            if (incoming.m_cost + outgoing.m_cost > max_cost) max_cost = incoming.m_cost + outgoing.m_cost;
        }

        // reset contraction data
        for (int& num : m_contr_data[0].m_reset_visited) m_contr_data[0].m_visited[num] = false;
        for (int& num : m_contr_data[0].m_reset_distances)
            m_contr_data[0].m_distances[num] = std::numeric_limits<int>::max();
        m_contr_data[0].m_reset_visited.clear();
        m_contr_data[0].m_reset_distances.clear();
    }

    // reset contraction data
    for (int& num : m_contr_data[0].m_reset_outgoing) m_contr_data[0].m_outgoing[num] = false;
    m_contr_data[0].m_reset_outgoing.clear();

    return 0.8 * max_cost + 0.2 * num_outgoing * num_incomming + 0.0 * num_deleted_neighbours[node];
}

// TODO: finish this and test it
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

        m_contr_data[0].m_reset_outgoing.push_back(outgoing.m_target);
        m_contr_data[0].m_outgoing[outgoing.m_target] = true;
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
            if (m_contr_data[0].m_visited[outgoing.m_target]) continue;
            if (contracted[outgoing.m_target]) {
                if (m_node_level[incoming.m_target] > max_neighbour_level)
                    max_neighbour_level = m_node_level[incoming.m_target];
                continue;
            }
            if (m_contr_data[0].m_distances[outgoing.m_target] != incoming.m_cost + outgoing.m_cost) continue;
            if (outgoing.m_target == incoming.m_target) continue;

            m_contr_data[0].m_visited[outgoing.m_target] = true;
            m_contr_data[0].m_reset_visited.push_back(outgoing.m_target);

            if (incoming.m_child_1 != -1) underlying_shortcuts += incoming.m_child_1;
            if (outgoing.m_child_1 != -1) underlying_shortcuts += outgoing.m_child_1;

            ++num_added_shortcuts;

            if (incoming.m_cost + outgoing.m_cost > max_cost) max_cost = incoming.m_cost + outgoing.m_cost;
        }

        // reset contraction data
        for (int& num : m_contr_data[0].m_reset_visited) m_contr_data[0].m_visited[num] = false;
        for (int& num : m_contr_data[0].m_reset_distances)
            m_contr_data[0].m_distances[num] = std::numeric_limits<int>::max();
        m_contr_data[0].m_reset_visited.clear();
        m_contr_data[0].m_reset_distances.clear();
    }

    // reset contraction data
    for (int& num : m_contr_data[0].m_reset_outgoing) m_contr_data[0].m_outgoing[num] = false;
    m_contr_data[0].m_reset_outgoing.clear();

    if (max_neighbour_level > cur_level)
        max_neighbour_level = 0;
    else
        ++max_neighbour_level;

    return 0.001 * max_cost + m_contr_data[0].m_num_contracted_neighbours[node] + 1 * underlying_shortcuts +
           2 * (num_added_shortcuts - (num_incomming + num_outgoing)) + 5 * max_neighbour_level;
}

void Graph::contractNode(std::vector<bool>& contracted, int contracted_node) {
    // mark outgoing nodes
    // will be used for pruning in contractionDijkstra
    int num_outgoing = 0;
    int max_distance_out = -1;
    int thread_num = omp_get_thread_num();
    for (auto& outgoing : m_graph[contracted_node]) {
        if (contracted[outgoing.m_target]) continue;
        if (outgoing.m_cost > max_distance_out) max_distance_out = outgoing.m_cost;

        ++m_contr_data[thread_num].m_num_contracted_neighbours[outgoing.m_target];

        ++num_outgoing;
        m_contr_data[thread_num].m_reset_outgoing.push_back(outgoing.m_target);
        m_contr_data[thread_num].m_outgoing[outgoing.m_target] = true;
    }

    for (auto& incoming : m_reverse_graph[contracted_node]) {
        if (contracted[incoming.m_target]) continue;
        int max_distance = incoming.m_cost + max_distance_out;

        ++m_contr_data[thread_num].m_num_contracted_neighbours[incoming.m_target];

        contractionDijkstra(incoming.m_target, contracted_node, contracted, num_outgoing, max_distance);

        for (auto& outgoing : m_graph[contracted_node]) {
            if (m_contr_data[thread_num].m_visited[outgoing.m_target]) continue;
            if (contracted[outgoing.m_target]) continue;
            if (m_contr_data[thread_num].m_distances[outgoing.m_target] != incoming.m_cost + outgoing.m_cost) continue;
            if (outgoing.m_target == incoming.m_target) continue;

            m_contr_data[thread_num].m_visited[outgoing.m_target] = true;
            m_contr_data[thread_num].m_reset_visited.push_back(outgoing.m_target);

            // add shortcut
            int underlying_arcs = 0;
            if (incoming.m_child_1 != -1) underlying_arcs += incoming.m_child_1;
            if (outgoing.m_child_1 != -1) underlying_arcs += outgoing.m_child_1;

            Edge fwd_shortcut(outgoing.m_target, incoming.m_cost + outgoing.m_cost, underlying_arcs, -1);
            m_contr_data[thread_num].m_shortcuts_fwd.push_back(std::make_pair(incoming.m_target, fwd_shortcut));
            Edge bwd_shortcut(incoming.m_target, incoming.m_cost + outgoing.m_cost, underlying_arcs, -1);
            m_contr_data[thread_num].m_shortcuts_bwd.push_back(std::make_pair(outgoing.m_target, bwd_shortcut));
        }

        // reset contraction data
        for (int& num : m_contr_data[thread_num].m_reset_visited) m_contr_data[thread_num].m_visited[num] = false;
        for (int& num : m_contr_data[thread_num].m_reset_distances)
            m_contr_data[thread_num].m_distances[num] = std::numeric_limits<int>::max();
        m_contr_data[thread_num].m_reset_visited.clear();
        m_contr_data[thread_num].m_reset_distances.clear();
    }

    // reset contraction data
    for (int& num : m_contr_data[thread_num].m_reset_outgoing) m_contr_data[thread_num].m_outgoing[num] = false;
    m_contr_data[thread_num].m_reset_outgoing.clear();
}

void Graph::contractionDijkstra(int start, int contracted_node, std::vector<bool>& contracted, int num_outgoing,
                                int max_distance) {
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;

    int num_visited_outgoing = 0;
    int thread_num = omp_get_thread_num();
    m_contr_data[thread_num].m_distances[start] = 0;
    m_contr_data[thread_num].m_reset_distances.push_back(start);
    pq.push(std::make_pair(0, start));

    while (!pq.empty()) {
        std::pair<int, int> cur_node = pq.top();
        pq.pop();

        if (m_contr_data[thread_num].m_distances[cur_node.second] != cur_node.first) continue;

        if (m_contr_data[thread_num].m_outgoing[cur_node.second]) ++num_visited_outgoing;
        if (cur_node.first > max_distance || num_visited_outgoing == num_outgoing) break;

        for (Edge& e : m_graph[cur_node.second]) {
            if (contracted[e.m_target]) continue;
            if (m_contr_data[thread_num].m_distances[e.m_target] >
                m_contr_data[thread_num].m_distances[cur_node.second] + e.m_cost) {
                if (m_contr_data[thread_num].m_distances[e.m_target] == std::numeric_limits<int>::max())
                    m_contr_data[thread_num].m_reset_distances.push_back(e.m_target);
                m_contr_data[thread_num].m_distances[e.m_target] =
                    m_contr_data[thread_num].m_distances[cur_node.second] + e.m_cost;
                pq.push(std::make_pair(m_contr_data[thread_num].m_distances[e.m_target], e.m_target));
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

        if (cur_node.first > lb_data.m_threshold) break;

        for (Edge& e : m_graph[cur_node.second]) {
            if (lb_data.m_marked[e.m_target]) continue;  // TODO: passt das?
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

int Graph::simplifiedHubLabelQuery(std::vector<std::pair<int, int>>& fwd_labels,
                                   std::vector<std::pair<int, int>>& bwd_labels) {
    int distance = std::numeric_limits<int>::max();
    auto fwd_iter = fwd_labels.begin();
    auto bwd_iter = bwd_labels.begin();

    while (fwd_iter != fwd_labels.end() && bwd_iter != bwd_labels.end()) {
        if (fwd_iter->first == bwd_iter->first) {
            if (fwd_iter->second + bwd_iter->second < distance) distance = fwd_iter->second + bwd_iter->second;
            ++fwd_iter;
            ++bwd_iter;
        } else if (fwd_iter->first < bwd_iter->first) {
            ++fwd_iter;
        } else {
            ++bwd_iter;
        }
    }

    return distance;
}

}  // namespace olsp
