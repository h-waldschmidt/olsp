#include <iostream>

#include "graph.h"

int main(int argc, char *argv[]) {
    // olsp::Graph g("/Users/helmut/Documents/BachelorArbeit/bachelorarbeit/data/stgtregbz_ch.fmi",
    // olsp::ReadMode::CONTRACTION_HIERACHIES);

    olsp::Graph g("/home/helmut/Documents/BachelorArbeit/bachelorarbeit/data/stgtregbz_ch.fmi",
                  olsp::ReadMode::CONTRACTION_HIERARCHY, false);

    int dist = olsp::Graph::dijkstraQuery(g.getGraphVec(), 377371, 754742);
    std::cout << "Distance: " << dist << std::endl;

    olsp::QueryData bd_data(377371, 754742, false);
    g.bidirectionalDijkstraQuery(bd_data);
    std::cout << "Distance: " << bd_data.m_distance << std::endl;
    std::vector<int> path = bd_data.m_shortest_path;
    int meeting_node = bd_data.m_meeting_node;
    std::cout << meeting_node << std::endl;
    g.bidirectionalDijkstraGetPath(bd_data);

    g.contractionHierachyQuery(bd_data);
    std::cout << "Distance: " << bd_data.m_distance << std::endl;
    std::cout << bd_data.m_meeting_node << std::endl;
    g.bidirectionalDijkstraGetPath(bd_data);
    g.createHubLabels();
    g.hubLabelQuery(bd_data);
    std::cout << "Distance: " << bd_data.m_distance << std::endl;
    std::cout << meeting_node << std::endl;

    std::cout << "Average Label size: " << g.averageLabelSize() << std::endl;
    std::cout << "Max Label size: " << g.maxLabelSize() << std::endl;
    return 0;
}
