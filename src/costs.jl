
"""
Compute the costs for moving from point `a` to a neighboring point `b`.
"""
function cost(a::CartesianIndex, b::CartesianIndex)
    Δ = b - a
    max(abs(Δ[1]), abs(Δ[2])) + abs(Δ[3])
end


"""
Heuristic to estimate the total cost for the path from `a` to `b`.

As long as the heuristic does not overestimate the true costs, A* finds an optimal path.
"""
function dist_heuristic(a::CartesianIndex, b::CartesianIndex)
    # spatial chebyshev distance beween points, plus time difference
    Δ = b - a
    max(abs(Δ[1]), abs(Δ[2])) + abs(Δ[3])
end


"""
Determine if point `p` is close enough to point`goal` to declare success.
"""
function isgoal(p::CartesianIndex, goal::CartesianIndex, tol::Int)
    Δ = goal - p
    p[3] == goal[3] &&                   # time must match exactly
        max(abs(Δ[1]), abs(Δ[2])) <= tol # we can be off by some pixel
end
