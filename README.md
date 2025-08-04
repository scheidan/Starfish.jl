# Starfish.jl: Compute Shortest Fish Trajectories

[![Build Status](https://github.com/scheidan/Starfish.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/scheidan/Starfish.jl/actions/workflows/CI.yml?query=branch%3Amain)


`Starfish` is a Julia package to find the shortest trajectory of a
fish between acoustic detections. It ensures that the path is always
compatible the depth observations. It is based on the [A* search
algorithm](https://en.wikipedia.org/wiki/A*_search_algorithm).

** THIS IS A PROTOTYPE! **


## Usage

See the docstring for ` find_shortest_trajectory` for details. The input format is
identical as for [`Wahoo.jl`](https://github.com/scheidan/Wahoo.jl).

### Example




## Acknowledgments

Thanks a lot to the authors of
[`AStarSearch.jl`](https://github.com/PaoloSarti/AStarSearch.jl) that
provide the backbone of this package.
