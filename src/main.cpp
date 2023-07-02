#include <omp.h>

#include <fstream>
#include <iostream>
#include <random>
#include <regex>

#include "graph.h"

int parseLine(std::string line) {
    return stoi(std::regex_replace(line, std::regex("[^0-9]*([0-9]+).*"), std::string("$1")));
}

// https://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
int getMemoryUsage() {  // Note: this value is in KB!
    int result = -1;
    std::ifstream file("/proc/self/status");
    std::string line;
    while (getline(file, line)) {
        if (line.find("VmSize") != std::string::npos) {
            result = parseLine(line);
            break;
        }
    }

    file.close();
    return result;
}

double convertTravelTimeToMeters(std::string graph_path, olsp::Heuristic heuristic) {
    olsp::Graph travel_time_graph(graph_path, olsp::ReadMode::NORMAL, true, true, 14, heuristic,
                                  olsp::DistanceMode::TRAVEL_TIME);
    olsp::Graph meter_graph(graph_path, olsp::ReadMode::NORMAL, true, true, 14, heuristic,
                            olsp::DistanceMode::DISTANCE_METERS);

    // travel_time_graph.createHubLabels();
    // meter_graph.createHubLabels();

    double avg_dist_meter = 0.0;
    double avg_dist_travel_time = 0.0;

    // randomly generate 100000 start/end pairs
    olsp::QueryData data(0, 0, travel_time_graph.getNumNodes(), false);
    std::vector<bool> marked_nodes(travel_time_graph.getNumNodes(), false);
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(0, travel_time_graph.getNumNodes() - 1);

    for (int i = 0; i < 100000; i++) {
        // calculate distance for all pairs in both graphs
        int start = dist(rng);
        while (marked_nodes[start]) start = dist(rng);
        int end = dist(rng);
        while (marked_nodes[end]) end = dist(rng);

        marked_nodes[start] = true;
        marked_nodes[end] = true;

        data.m_start = start;
        data.m_end = end;
        travel_time_graph.contractionHierachyQuery(data);
        // data.m_distance = olsp::Graph::dijkstraQuery(travel_time_graph.getGraph(), start, end);
        if (data.m_distance == std::numeric_limits<int>::max()) continue;
        // std::cout << "Distances Travel Time: " << data.m_distance << std::endl;

        avg_dist_travel_time += data.m_distance;

        // data.m_distance = olsp::Graph::dijkstraQuery(meter_graph.getGraph(), start, end);
        meter_graph.contractionHierachyQuery(data);
        if (data.m_distance == std::numeric_limits<int>::max()) continue;

        avg_dist_meter += data.m_distance;
        // std::cout << "Distances Meters: " << data.m_distance << std::endl;
    }

    avg_dist_meter /= 100000.0;
    avg_dist_travel_time /= 100000.0;

    // average how much a meter converts into travel_time
    return avg_dist_meter / avg_dist_travel_time;
}

void test() {
    // double stgt_conversion =
    //     convertTravelTimeToMeters("/home/helmut/Documents/BachelorArbeit/bachelorarbeit/data/stgtregbz.fmi");
    // double german_conversion =
    //    convertTravelTimeToMeters("/home/helmut/Documents/BachelorArbeit/bachelorarbeit/data/bw.fmi");

    // std::cout << "Stuttgart Conversion Rate: " << stgt_conversion << "\n";
    // std::cout << "Germany Conversion Rate: " << german_conversion << "\n";

    // olsp::Graph g("/Users/helmut/Documents/BachelorArbeit/bachelorarbeit/data/stgtregbz_ch.fmi",
    // olsp::ReadMode::CONTRACTION_HIERACHIES);

    // don't prune graph when using advanced hub label creation
    olsp::Graph g("/home/helmut/Documents/BachelorArbeit/bachelorarbeit/data/stgtregbz.fmi", olsp::ReadMode::NORMAL,
                  true, false, 14, olsp::Heuristic::WEIGHTED_COST, olsp::DistanceMode::TRAVEL_TIME);

    // int dist = olsp::Graph::dijkstraQuery(g.getGraphVec(), 377371, 754742);
    // std::cout << "Distance: " << dist << std::endl;
    {
        olsp::QueryData bd_data(377371, 754742, g.getNumNodes(), false);

        g.bidirectionalDijkstraQuery(bd_data);
        std::cout << "Distance: " << bd_data.m_distance << std::endl;
        std::vector<int> path = bd_data.m_shortest_path;
        int meeting_node = bd_data.m_meeting_node;
        std::cout << meeting_node << std::endl;
        // g.bidirectionalDijkstraGetPath(bd_data);
    }
    {
        olsp::QueryData bd_data(377371, 754742, g.getNumNodes(), false);
        g.contractionHierachyQuery(bd_data);
        std::cout << "Distance: " << bd_data.m_distance << std::endl;
        std::cout << bd_data.m_meeting_node << std::endl;
        // g.bidirectionalDijkstraGetPath(bd_data);
    }

    int threshold = 4000;

    g.createHubLabels();
    {
        olsp::QueryData bd_data(377371, 754742, g.getNumNodes(), false);
        g.hubLabelQuery(bd_data);
        std::cout << "Distance: " << bd_data.m_distance << std::endl;
        std::cout << bd_data.m_meeting_node << std::endl;
    }

    std::cout << "Average Label size: " << g.averageLabelSize() << std::endl;
    std::cout << "Max Label size: " << g.maxLabelSize() << std::endl;

    std::cout << "Num Labels with weight between: " << threshold / 2 << " and " << threshold << " : "
              << g.numHubLabelsInRange(threshold / 2, threshold) << std::endl;

    auto path_cover = g.createShortestPathCover(threshold);
    std::cout << "Path Cover Size: " << path_cover.size() << std::endl;

    // auto new_path_cover = g.reducePathCover(path_cover, threshold);
    // std::cout << "New Path Cover Size: " << new_path_cover.size() << std::endl;

    // auto lower_bound = g.verifyShortestPathCover(path_cover, threshold);
    // std::cout << "Path Cover Valid?: " << lower_bound;
}

void singleBenchmark(std::string graph_path, olsp::DistanceMode metric, olsp::Heuristic heuristic, int threshold) {
    if (metric == olsp::DistanceMode::TRAVEL_TIME) {
        double conversion = convertTravelTimeToMeters(graph_path, heuristic);
        std::cout << "Conversion: " << conversion << std::endl;
        threshold = static_cast<int>(static_cast<double>(threshold) / conversion);
    }

    std::cout << "Threshold: " << threshold << std::endl;
    olsp::Graph g(graph_path, olsp::ReadMode::NORMAL, true, false, 14, heuristic, metric);

    std::cout << "Graph Memory usage in kb: " << getMemoryUsage() << std::endl;

    // g.createHubLabels(threshold);
    g.createHubLabels();
    std::cout << "Graph and Label Memory usage in kb: " << getMemoryUsage() << std::endl;
    int avg_hub_label = g.averageLabelSize();
    std::cout << "Avg. Label Size: " << avg_hub_label << std::endl;
    int max_hub_label = g.maxLabelSize();
    std::cout << "Max Label Size: " << max_hub_label << std::endl;

    std::vector<int> path_cover = g.createShortestPathCover(threshold);
    std::cout << "Path cover size: " << path_cover.size() << std::endl;

    g.clearHubLabel();

    std::vector<int> lower_bound = g.lowerBound(path_cover, threshold);
    std::cout << "Lower bound size: " << lower_bound.size() << std::endl;
}

void benchmark() {
    // TODO:
    std::string stuttgart_path = "data/stgtregbz.fmi";
    std::string bw_path = "data/bw.fmi";
    std::string germany_path = "data/germany.fmi";

    int small_threshold = 40000;
    int big_threshold = 125000;

    // Stuttgart IN_OUT
    {
        std::cout << "Stuttgart Graph with IN_OUT and TravelTime." << std::endl;
        singleBenchmark(stuttgart_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::IN_OUT, small_threshold);
    }

    {
        std::cout << "Stuttgart Graph with IN_OUT and Meter-Metric." << std::endl;
        singleBenchmark(stuttgart_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::IN_OUT, small_threshold);
    }

    // Stuttgart Edge Difference
    {
        std::cout << "Stuttgart Graph with EDGE_DIFFERENCE and TravelTime." << std::endl;
        singleBenchmark(stuttgart_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::EDGE_DIFFERENCE,
                        small_threshold);
    }

    {
        std::cout << "Stuttgart Graph with EDGE_DIFFERENCE and Meter-Metric." << std::endl;
        singleBenchmark(stuttgart_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::EDGE_DIFFERENCE,
                        small_threshold);
    }

    // Stuttgart weighted Cost
    {
        std::cout << "Stuttgart Graph with WEIGHTED_COST and TravelTime." << std::endl;
        singleBenchmark(stuttgart_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::WEIGHTED_COST,
                        small_threshold);
    }

    {
        std::cout << "Stuttgart Graph with WEIGHTED_COST and Meter-Metric." << std::endl;
        singleBenchmark(stuttgart_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::WEIGHTED_COST,
                        small_threshold);
    }

    // Stuttgart microsoft
    {
        std::cout << "Stuttgart Graph with MICROSOFT and TravelTime." << std::endl;
        singleBenchmark(stuttgart_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::MICROSOFT, small_threshold);
    }

    {
        std::cout << "Stuttgart Graph with MICROSOFT and Meter-Metric." << std::endl;
        singleBenchmark(stuttgart_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::MICROSOFT,
                        small_threshold);
    }

    // BW IN_OUT
    {
        std::cout << "BW Graph with IN_OUT and TravelTime." << std::endl;
        singleBenchmark(bw_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::IN_OUT, big_threshold);
    }

    {
        std::cout << "BW Graph with IN_OUT and Meter-Metric." << std::endl;
        singleBenchmark(bw_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::IN_OUT, big_threshold);
    }

    // BW Edge Difference
    {
        std::cout << "BW Graph with edge difference and TravelTime." << std::endl;
        singleBenchmark(bw_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::EDGE_DIFFERENCE, big_threshold);
    }

    {
        std::cout << "BW Graph with edge difference and Meter-Metric." << std::endl;
        singleBenchmark(bw_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::EDGE_DIFFERENCE, big_threshold);
    }

    // BW weighted Cost
    {
        std::cout << "BW Graph with weighted cost and TravelTime." << std::endl;
        singleBenchmark(bw_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::WEIGHTED_COST, big_threshold);
    }

    {
        std::cout << "BW Graph with weighted cost and Meter-Metric." << std::endl;
        singleBenchmark(bw_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::WEIGHTED_COST, big_threshold);
    }

    // BW microsoft
    {
        std::cout << "BW Graph with Microsoft and TravelTime." << std::endl;
        singleBenchmark(bw_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::MICROSOFT, big_threshold);
    }

    {
        std::cout << "BW Graph with Microsoft and Meter-Metric." << std::endl;
        singleBenchmark(bw_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::MICROSOFT, big_threshold);
    }

    // Germnay weighted Cost
    {
        std::cout << "Germany Graph with weighted cost and TravelTime." << std::endl;
        singleBenchmark(germany_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::WEIGHTED_COST, big_threshold);
    }

    {
        std::cout << "Germany Graph with weighted cost and Meter-Metric." << std::endl;
        singleBenchmark(germany_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::WEIGHTED_COST,
                        big_threshold);
    }

    // Germany microsoft
    {
        std::cout << "Germany Graph with Microsoft and TravelTime." << std::endl;
        singleBenchmark(germany_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::MICROSOFT, big_threshold);
    }

    {
        std::cout << "Germany Graph with Microsoft and Meter-Metric." << std::endl;
        singleBenchmark(germany_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::MICROSOFT, big_threshold);
    }
}

void germanyBenchmark() {
    std::string germany_path = "data/germany.fmi";

    int big_threshold = 125000;

    // Germnay weighted Cost
    {
        std::cout << "Germany Graph with weighted cost and TravelTime." << std::endl;
        singleBenchmark(germany_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::WEIGHTED_COST, big_threshold);
    }

    {
        std::cout << "Germany Graph with weighted cost and Meter-Metric." << std::endl;
        singleBenchmark(germany_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::WEIGHTED_COST,
                        big_threshold);
    }

    // Germany microsoft
    {
        std::cout << "Germany Graph with Microsoft and TravelTime." << std::endl;
        singleBenchmark(germany_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::MICROSOFT, big_threshold);
    }

    {
        std::cout << "Germany Graph with Microsoft and Meter-Metric." << std::endl;
        singleBenchmark(germany_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::MICROSOFT, big_threshold);
    }
}

int main(int argc, char *argv[]) {
    omp_set_num_threads(14);

    germanyBenchmark();
    return 0;
}