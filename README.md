# Starfish.jl: Compute Shortest Fish Trajectories 🌟🐠

[![Build Status](https://github.com/scheidan/Starfish.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/scheidan/Starfish.jl/actions/workflows/CI.yml?query=branch%3Amain)


`Starfish` is a Julia package to find the shortest trajectory of a
fish between acoustic detections. It ensures that the path is always
compatible the depth observations. It is based on the [A* search
algorithm](https://en.wikipedia.org/wiki/A*_search_algorithm).

** THIS IS A PROTOTYPE! **

## Installation

This package is not yet registered. Install it with:
```julia
] add https://github.com/scheidan/Starfish.jl
```

## Usage

See the docstring for `find_shortest_trajectory` for details. The input format is
identical as for [`Wahoo.jl`](https://github.com/scheidan/Wahoo.jl).

### Example

A minimal example with artificial data:

```julia
using Starfish

import GeoArrays
using DelimitedFiles: readdlm

# -----------
# Read data

# Read example data that come with the package
pathdata = joinpath(pkgdir(Starfish), "example_data")

# -- bathymetry
bathymetry_map = GeoArrays.read(joinpath(pathdata, "bathymetry_200m.tif"))


# -- depth observations
depth_signals = readdlm(joinpath(pathdata, "depth_observations.csv"), ',', header=true)[1][:,2]


# -- acoustic observations

# read signals (1: detections, 0: no detection, -1: inactive)
acoustic_allsignals = readdlm(joinpath(pathdata, "acoustic_observations.csv"), ',', header=true)[1][:,2:3]
acoustic_allsignals = Int.(acoustic_allsignals')

# vector of signals:
acoustic_signals = [acoustic_allsignals[1,:],
                    acoustic_allsignals[2,:]]

# read positions
moorings = readdlm(joinpath(pathdata, "acoustic_moorings.csv"), ',', header=true)[1]
acoustic_pos = tuple.(moorings[:,2], moorings[:,3])



# -----------
# Find path

res = find_shortest_trajectory(bathymetry_map,
                               acoustic_signals, acoustic_pos,
                               depth_signals,
                               goal_tol = 2,
                               abstol = 10, reltol=0.1);

res.path_length


# -----------
# Visualize path

using CairoMakie

fig = Figure();
ax  = Axis(fig[1, 1]);
hm = heatmap!(ax, bathymetry_map, colorrange = (0.0, 300),
              lowclip = :honeydew3, colormap=:dense);
scatter!(ax, acoustic_pos, color=:orange, markersize = 10);
xy = Point2.(res.path);
scatter!(ax, xy; markersize = 3,
         color = :red);
Colorbar(fig[:, end+1], hm);
fig
```


## References

The example bathymetry data is derived from the following survey:

Howe JA, Anderton R, Arosio R, et al. The seabed geomorphology and geological structure of the Firth of Lorn, western Scotland, UK, as revealed by multibeam echo-sounder survey. Earth and Environmental Science Transactions of the Royal Society of Edinburgh. 2014;105(4):273-284. https://doi.org/10.1017/S1755691015000146


## Acknowledgments

Many thanks to the authors of [`AStarSearch.jl`](https://github.com/PaoloSarti/AStarSearch.jl)!
