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
- `path::Vector{Tuple{Float64, Float64}}`: The coordinates defining the path for each time step.
                                           Time steps for which no path was found hold `(NaN, NaN)`.
- `path_length`: total (spatial) length of the found path.
- `costs`: total cost of the found path.
"""
function find_shortest_trajectory(bathymetry::GeoArrays.GeoArray,
                                  acoustic_signals, acoustic_pos,
                                  depth_signals;
                                  goal_tol::Int = 0,
                                  reltol = 0.0, abstol = 0.0)
    # --
    # 1) get time and location of all accoustic detections
    # N.B. if we have multiple detections at the same time, only one is used!
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
    total_length = 0
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

        if result.status == :success
            push!(all_path, result.path)
            total_costs += result.cost
            total_length += result.cost - (obs_points_ci[i+1] - obs_points_ci[i])[3]
        else
            @warn "No path found from $(obs_points_ci[i].I[1:2]) to $(obs_points_ci[i+1].I[1:2]), time = $(obs_time[i]):$(obs_time[i+1])!"
        end
        next!(p)
    end

    # --
    # 4) concatenate to one path

    # convert to coordinates and concatenate
    path = fill((NaN, NaN), length(depth_signals))
    for p in all_path
        idxs = p[1][3]:p[end][3]
        path[idxs] .= [Tuple(GeoArrays.coords(bathymetry, (c[1], c[2]))) for c in p]
    end

    # --
    (path = path,
     path_length = total_length,
     costs = total_costs
     )
end


end
