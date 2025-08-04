# -----
# defines function that returns the neighborhood, that is all legal moves

# Return Moore neighborhood with max distance `d`
moore_neighborhood(d::Int) = (CartesianIndex(i, j, 1) for i in -d:d, j in -d:d)

function get_neighbor(pos, bathymetry, depth_signals; reltol = 0.1, abstol = 0.0)

    res = CartesianIndex[]
    for d in (CartesianIndex(0, 0, 1),                             # do not move
              CartesianIndex(1, 0, 1), CartesianIndex(-1, 0, 1),   # up, down
              CartesianIndex(0, 1, 1), CartesianIndex(0, -1, 1),   # left, right
              CartesianIndex(1, 1, 1), CartesianIndex(1, -1, 1),   # diagonals
              CartesianIndex(-1, 1, 1), CartesianIndex(-1, -1, 1)
              )

        n = pos + d
        p = CartesianIndex(n[1], n[2])
        time = n[3]

        is_accessible = (1 ≤ p[1] ≤ size(bathymetry)[1]) &&
            (1 ≤ p[2] ≤ size(bathymetry)[2]) &&
            time <= length(depth_signals) &&
            (bathymetry[p] > 0) &&
            (max(abstol, bathymetry[p] * (1+reltol)) > depth_signals[time])

        if is_accessible
            push!(res, n)
        end
    end
    return res
end

# helper to make a typestable closure
function make_neighbor_function(bathymetry::AbstractMatrix{Tb}, depth::AbstractVector{Td};
                                reltol::Tr = 0.1, abstol::Ta = 10.0) where {Tb,Td,Tr,Ta}
    (pos::CartesianIndex{3}) -> get_neighbor(pos, bathymetry, depth; reltol=reltol, abstol=abstol)
end
