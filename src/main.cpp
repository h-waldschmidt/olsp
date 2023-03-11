#include <iostream>

#include "graph.h"

int main(int argc, char *argv[]) {
    // olsp::Graph g("/Users/helmut/Documents/BachelorArbeit/bachelorarbeit/data/stgtregbz_ch.fmi",
    // olsp::ReadMode::CONTRACTION_HIERACHIES);

    olsp::Graph g("/Users/helmut/Documents/BachelorArbeit/bachelorarbeit/data/bw.fmi", olsp::ReadMode::NORMAL);
    olsp::BiDirectionalDijkstraData bd_data(3600519, 3598035);
    g.bidirectionalDijkstraCalculateDistance(bd_data);
    std::cout << "Distance: " << bd_data.m_distance << std::endl;

    return 0;
}
