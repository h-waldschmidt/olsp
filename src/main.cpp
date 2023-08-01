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

/**
 * @brief Returns conversion rate. (How many meters convert to one travel time (1/100 s))
 *
 * @param graph_path specifies location of path
 * @param heuristic  specifies the used heuristic
 * @param is         specifies whether IS are used or not
 * @return conversion rate
 */
double getConversionRate(std::string graph_path, olsp::Heuristic heuristic, bool is) {
    olsp::Graph travel_time_graph;
    olsp::Graph meter_graph;

    if (is) {
        travel_time_graph =
            olsp::Graph(graph_path, olsp::ReadMode::NORMAL, true, true, 14, heuristic, olsp::DistanceMode::TRAVEL_TIME);
        meter_graph = olsp::Graph(graph_path, olsp::ReadMode::NORMAL, true, true, 14, heuristic,
                                  olsp::DistanceMode::DISTANCE_METERS);
    } else {
        travel_time_graph =
            olsp::Graph(graph_path, olsp::ReadMode::NORMAL, true, true, heuristic, olsp::DistanceMode::TRAVEL_TIME);
        meter_graph =
            olsp::Graph(graph_path, olsp::ReadMode::NORMAL, true, true, heuristic, olsp::DistanceMode::DISTANCE_METERS);
    }

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
        if (data.m_distance == std::numeric_limits<int>::max()) continue;

        avg_dist_travel_time += data.m_distance;

        meter_graph.contractionHierachyQuery(data);
        if (data.m_distance == std::numeric_limits<int>::max()) continue;

        avg_dist_meter += data.m_distance;
    }

    avg_dist_meter /= 100000.0;
    avg_dist_travel_time /= 100000.0;

    // average how much a meter converts into travel_time
    return avg_dist_meter / avg_dist_travel_time;
}

void example() {
    olsp::Graph g("data/stgtregbz.fmi", olsp::ReadMode::NORMAL, true, false, olsp::Heuristic::IN_OUT,
                  olsp::DistanceMode::TRAVEL_TIME);

    int dist = olsp::Graph::dijkstraQuery(g.getGraph(), 377371, 754742);
    std::cout << "Distance: " << dist << std::endl;
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

    int threshold = 40000;

    g.createHubLabelsWithoutIS();
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

    auto lower_bound = g.verifyShortestPathCover(path_cover, threshold);
    std::cout << "Path Cover Valid?: " << lower_bound;
}

/**
 * @brief Run one single benchmark permutation.
 * Creates CH, Hub-Labeling, ESC Set and Lower-Bound for that set.
 *
 * @param graph_path specifies location of path
 * @param metric     specifies whether travel time or distance in meters is used
 * @param heuristic  specifies the used heuristic
 * @param is         specifies whether IS are used or not
 * @param threshold  specifies the range of the electric vehicle
 * @param conversion specifies the conversion rate from meters to travel time
 */
void singleBenchmark(std::string graph_path, olsp::DistanceMode metric, olsp::Heuristic heuristic, bool is,
                     int threshold, double conversion = 0) {
    if (metric == olsp::DistanceMode::TRAVEL_TIME && conversion == 0) {
        double new_conversion = getConversionRate(graph_path, heuristic, is);
        std::cout << "Conversion: " << new_conversion << std::endl;
        threshold = static_cast<int>(static_cast<double>(threshold) / new_conversion);
    } else if (metric == olsp::DistanceMode::TRAVEL_TIME && conversion != 0) {
        std::cout << "Conversion: " << conversion << std::endl;
        threshold = static_cast<int>(static_cast<double>(threshold) / conversion);
    }

    std::cout << "Threshold: " << threshold << std::endl;

    olsp::Graph g;
    if (is)
        g = olsp::Graph(graph_path, olsp::ReadMode::NORMAL, true, false, 14, heuristic, metric);
    else
        g = olsp::Graph(graph_path, olsp::ReadMode::NORMAL, true, false, heuristic, metric);

    std::cout << "Graph Memory usage in kb: " << getMemoryUsage() << std::endl;

    // g.createHubLabels(threshold); // pass threshold to use pruned hub-labeling
    if (is)
        g.createHubLabelsWithIS();
    else
        g.createHubLabelsWithoutIS();

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

/**
 * @brief Runs various permutations for all three graphs
 *
 * @param is specifies whether IS are used or not
 */
void benchmark(bool is = false) {
    std::string stuttgart_path = "data/stgtregbz.fmi";
    std::string bw_path = "data/bw.fmi";
    std::string germany_path = "data/germany.fmi";

    int small_threshold = 40000;
    int big_threshold = 125000;

    // average conversion_rate over 8 runs
    double stuttgart_conversion = 0.22345675;
    double bw_conversion = 0.221872;
    // average conversion_rate over 3 runs
    double germany_conversion = 0.229944333;
    // Stuttgart IN_OUT

    {
        std::cout << "Stuttgart Graph with IN_OUT and TravelTime." << std::endl;
        singleBenchmark(stuttgart_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::IN_OUT, is, small_threshold,
                        stuttgart_conversion);
    }

    {
        std::cout << "Stuttgart Graph with IN_OUT and Meter-Metric." << std::endl;
        singleBenchmark(stuttgart_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::IN_OUT, is,
                        small_threshold);
    }

    // Stuttgart Edge Difference
    {
        std::cout << "Stuttgart Graph with EDGE_DIFFERENCE and TravelTime." << std::endl;
        singleBenchmark(stuttgart_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::EDGE_DIFFERENCE, is,
                        small_threshold, stuttgart_conversion);
    }

    {
        std::cout << "Stuttgart Graph with EDGE_DIFFERENCE and Meter-Metric." << std::endl;
        singleBenchmark(stuttgart_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::EDGE_DIFFERENCE, is,
                        small_threshold);
    }

    // Stuttgart weighted Cost
    {
        std::cout << "Stuttgart Graph with WEIGHTED_COST and TravelTime." << std::endl;
        singleBenchmark(stuttgart_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::WEIGHTED_COST, is,
                        small_threshold, stuttgart_conversion);
    }

    {
        std::cout << "Stuttgart Graph with WEIGHTED_COST and Meter-Metric." << std::endl;
        singleBenchmark(stuttgart_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::WEIGHTED_COST, is,
                        small_threshold);
    }

    // Stuttgart microsoft
    {
        std::cout << "Stuttgart Graph with MICROSOFT and TravelTime." << std::endl;
        singleBenchmark(stuttgart_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::MICROSOFT, is,
                        small_threshold, stuttgart_conversion);
    }

    {
        std::cout << "Stuttgart Graph with MICROSOFT and Meter-Metric." << std::endl;
        singleBenchmark(stuttgart_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::MICROSOFT, is,
                        small_threshold);
    }

    // BW IN_OUT
    {
        std::cout << "BW Graph with IN_OUT and TravelTime." << std::endl;
        singleBenchmark(bw_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::IN_OUT, is, big_threshold,
                        bw_conversion);
    }

    {
        std::cout << "BW Graph with IN_OUT and Meter-Metric." << std::endl;
        singleBenchmark(bw_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::IN_OUT, is, big_threshold);
    }

    // BW Edge Difference
    {
        std::cout << "BW Graph with edge difference and TravelTime." << std::endl;
        singleBenchmark(bw_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::EDGE_DIFFERENCE, is, big_threshold,
                        bw_conversion);
    }

    {
        std::cout << "BW Graph with edge difference and Meter-Metric." << std::endl;
        singleBenchmark(bw_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::EDGE_DIFFERENCE, is,
                        big_threshold);
    }

    // BW weighted Cost
    {
        std::cout << "BW Graph with weighted cost and TravelTime." << std::endl;
        singleBenchmark(bw_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::WEIGHTED_COST, is, big_threshold,
                        bw_conversion);
    }

    {
        std::cout << "BW Graph with weighted cost and Meter-Metric." << std::endl;
        singleBenchmark(bw_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::WEIGHTED_COST, is,
                        big_threshold);
    }

    // BW microsoft
    {
        std::cout << "BW Graph with Microsoft and TravelTime." << std::endl;
        singleBenchmark(bw_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::MICROSOFT, is, big_threshold,
                        bw_conversion);
    }

    {
        std::cout << "BW Graph with Microsoft and Meter-Metric." << std::endl;
        singleBenchmark(bw_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::MICROSOFT, is, big_threshold);
    }

    // Germnay weighted Cost
    {
        std::cout << "Germany Graph with weighted cost and TravelTime." << std::endl;
        singleBenchmark(germany_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::WEIGHTED_COST, is,
                        big_threshold, germany_conversion);
    }

    {
        std::cout << "Germany Graph with weighted cost and Meter-Metric." << std::endl;
        singleBenchmark(germany_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::WEIGHTED_COST, is,
                        big_threshold);
    }

    // Germany microsoft
    {
        std::cout << "Germany Graph with Microsoft and TravelTime." << std::endl;
        singleBenchmark(germany_path, olsp::DistanceMode::TRAVEL_TIME, olsp::Heuristic::MICROSOFT, is, big_threshold,
                        germany_conversion);
    }

    {
        std::cout << "Germany Graph with Microsoft and Meter-Metric." << std::endl;
        singleBenchmark(germany_path, olsp::DistanceMode::DISTANCE_METERS, olsp::Heuristic::MICROSOFT, is,
                        big_threshold);
    }
}

int main(int argc, char* argv[]) {
    // example(); // call example benchmark
    benchmark();  // tests permutations without IS
    // benchmark(true); // tests permutations wit IS

    return 0;
}
