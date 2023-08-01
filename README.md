# Optimal Loading Station Placement (OLSP)

## How to run?

1. Get Stuttgart, Baden-Württemberg and Deutschland graph from [here](https://fmi.uni-stuttgart.de/alg/research/stuff/)
2. Create `data` folder and place graph files in there
3. Compile with `sh build.sh -r`
4. Run with `./build/olsp` (Default: runs predefined benchmarks for all graphs without independent sets)
5. (To define custom benchmarks edit `src/main.cpp` file (see `example` function). Available graph functions can be found in `src/graph.h`)

## Issues

-   When using the function `createHubLabelsWithIS()` it sometimes produces segfaults. It isn't really predictable either. One one run every permutation works and on others it randomly crashes. The correctness of the produced hub labeling isn't impacted (With no crash point to point queries result in shortest distances and the produced esc is valid). It has probably something to do with the parallelization.
