# -----
# defines function that returns the neighborhood, that is all legal moves

"""
Moore neighborhood with max distance `d`
"""
function moore_neighborhood(d::Int)
    vals = sort(-d:d, by = abs)   # having "do not move" first somehow speeds up the search
    (CartesianIndex(i, j, 1) for i in vals, j in vals)
end


function get_neighbor(pos, bathymetry, depth_signals; seabed_tol = 0.0, benthic_tol = Inf)

    d = 1
    res = CartesianIndex[]
    sizehint!(res, (2*d+1)^2)

    for d in moore_neighborhood(d)

        n = pos + d
        p = CartesianIndex(n[1], n[2])
        time = n[3]

        is_accessible = (1 ≤ p[1] ≤ size(bathymetry)[1]) &&
            (1 ≤ p[2] ≤ size(bathymetry)[2]) &&
            time <= length(depth_signals) &&
            (bathymetry[p] > 0) &&
            bathymetry[p] + seabed_tol > depth_signals[time] &&
            bathymetry[p] - depth_signals[time] < benthic_tol

        if is_accessible
            push!(res, n)
        end
    end
    return res
end

# helper to make a typestable closure
function make_neighbor_function(bathymetry::AbstractMatrix{Tb}, depth::AbstractVector{Td};
                                benthic_tol::Tr = 0.1, seabed_tol::Ta = 10.0) where {Tb,Td,Tr,Ta}
    (pos::CartesianIndex{3}) -> get_neighbor(pos, bathymetry, depth; seabed_tol = seabed_tol, benthic_tol = benthic_tol)
end
