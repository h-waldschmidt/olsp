# Optimal Loading Station Placement (OLSP)

## How to run?

1. Get Stuttgart, Baden-Württemberg and Deutschland graph from [here](https://fmi.uni-stuttgart.de/alg/research/stuff/)
2. Create `data` folder and place graph files in there
3. Compile with `sh build.sh -r`
4. Run with `./build/olsp` (Default: runs predefined benchmarks for all graphs without independent sets)
5. (To define custom benchmarks edit `src/main.cpp` file (see `example` function). Available graph functions can be found in `src/graph.h`)
