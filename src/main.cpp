#include <iostream>

#include "graph.h"

int main(int argc, char *argv[]) {
    // olsp::Graph g("/Users/helmut/Documents/BachelorArbeit/bachelorarbeit/data/stgtregbz_ch.fmi",
    // olsp::ReadMode::CONTRACTION_HIERACHIES);

    olsp::Graph g("/home/helmut/Documents/BachelorArbeit/bachelorarbeit/data/germany.fmi", olsp::ReadMode::NORMAL);

    olsp::QueryData bd_data(10444777, 21255267, false);
    g.bidirectionalDijkstraCalculateDistance(bd_data);
    std::cout << "Distance: " << bd_data.m_distance << std::endl;

    g.bidirectionalDijkstraGetPath(bd_data);
    std::cout << "Meeting node: " << bd_data.m_meeting_node << std::endl;

    // for (int i = 0; i < bd_data.m_shortest_path.size(); i++) {
    //     std::cout << bd_data.m_shortest_path[i] << std::endl;
    // }

    return 0;
}
