module Starfish

include("neighborhood.jl")
include("costs.jl")

export find_shortest_trajectory

import AStarSearch
import GeoArrays
using ProgressMeter: Progress, next!



"""
    find_shortest_trajectory(
        bathymetry::GeoArrays.GeoArray,
        acoustic_signals,
        acoustic_pos,
        depth_signals;
        goal_tol::Int = 0,
        reltol::Float64 = 0.0,
        abstol::Float64 = 0.0
    )

Compute the shortest path through a bathymetry grid that visits all
acoustic observation locations in order, while enforcing time-varying
depth constraints. The inputs are fomats are identical as for
Wahoo.jl.

Note, only a path from the first to the last accoustic observation is returned.

## Arguments

- `bathymetry::GeoArrays.GeoArray`: Raster of bathyemtry, see Wahoo.jl.
- `acoustic_signals`: Collection of acoustic signal measurements, see Wahoo.jl.
- `acoustic_pos`: Collection of spatial coordinates, see Wahoo.jl.
- `depth_signals`: Vector of depths measurements, see Wahoo.jl.

- `goal_tol::Int = 0`: tolerance for reaching each acoustic target location measured in pixels.
                       Controls the receiver range.
- `reltol::Float64 = 0.1`: Relative tolerance for depth measurements
- `abstol::Float64 = 0.0`: Absolute tolerance for depth measurements


## Return Values
- `path::Vector{Tuple{Float64, Float64}}`: coordinates defining the path
- `time_start`: time index of begining of path
- `time_end`: time index of end of path
- `path_length`: the (spatial) length of the whole path measured in "pixel steps"
"""
function find_shortest_trajectory(bathymetry::GeoArrays.GeoArray,
                                  acoustic_signals, acoustic_pos,
                                  depth_signals;
                                  goal_tol::Int = 0,
                                  reltol = 0.0, abstol = 0.0)
    # --
    # 1) get time and location of all accoustic detections
    # N.B. if we have multibe detection at the same time, only one is used!
    obs_points = []
    obs_time = []
    for s in 1:length(acoustic_signals)
        for t in 1:length(acoustic_signals[1])
            if acoustic_signals[s][t] == 1 && t ∉ obs_time
                push!(obs_points, acoustic_pos[s])
                push!(obs_time, t)
            end
        end
    end

    # sort over time
    perm = sortperm(obs_time)
    obs_time = obs_time[perm]
    obs_points = obs_points[perm]

    @info("Found $(length(obs_points)) acoustic observations ($(length(unique(obs_points))) locations) within times $(extrema(obs_time)).")

    # --
    # 2) convert in cartesian index
    ci = [GeoArrays.indices(bathymetry, c) for c in obs_points]

    obs_points_ci = [CartesianIndex(ci[i][1], ci[i][2], obs_time[i]) for i in 1:length(obs_points)]

    # --
    # 3) find paths
    _get_neighbor = make_neighbor_function(bathymetry.A,  depth_signals; reltol = reltol, abstol = abstol)

    _isgoal(p, g)::Bool = isgoal(p, g, goal_tol)

    all_path = []
    total_costs = 0
    p = Progress(length(obs_points_ci)-1; dt=0.1, desc="Find paths...")
    for i in 1:(length(obs_points_ci)-1)

        # the maximal costs per time step is 2
        maxcosts = 2 * (obs_points_ci[i+1] - obs_points_ci[i])[3]

        result =  AStarSearch.astar(_get_neighbor,
                                    obs_points_ci[i], obs_points_ci[i+1],
                                    heuristic = dist_heuristic,
                                    cost = cost,
                                    isgoal =  _isgoal,
                                    maxcost = maxcosts # limits search space
                                    )

        if result.status != :success
            println(result.status)
            @warn "no path found from time $(obs_time[i]) to $(obs_time[i+1])!",  obs_points_ci[i], obs_points_ci[i+1]
            break
        end
        push!(all_path, result.path)
        total_costs += result.cost
        next!(p)
    end

    # --
    # 4) concatenate to one path

    # convert to coordinates and concatenate
    all_path_coord = [[Tuple(GeoArrays.coords(bathymetry, (c[1], c[2]))) for c in p] for p in all_path]

    path = [all_path_coord[1][1]]
    for p in all_path_coord
        append!(path,  p[2:end])
    end


    (path = path,
     time_start = minimum(obs_time),
     time_end =  minimum(obs_time) + length(path) - 1,
     path_length = total_costs - length(path)
     )
end


end
