# Optimal Loading Station Placement (OLSP)
This code is a part of my bachelor thesis [[1]](#1).
## How to run?

1. Get Stuttgart, Baden-Württemberg and Deutschland graph from [here](https://fmi.uni-stuttgart.de/alg/research/stuff/)
2. Create `data` folder and place graph files in there
3. Compile with `sh build.sh -r`
4. Run with `./build/olsp` (Default: runs predefined benchmarks for all graphs without independent sets)
5. (To define custom benchmarks edit `src/main.cpp` file (see `example` function). Available graph functions can be found in `src/graph.h`)

## References

<a id="1">[1]</a> 
Waldschmidt Helmut
"Optimized placement of charging stations for electric cars"
http://dx.doi.org/10.18419/opus-13829
