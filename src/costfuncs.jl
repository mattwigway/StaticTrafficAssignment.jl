module VDF
import ..RoadClass, ..EdgeData

abstract type DelayFunc end

"""
    BPR(α, β)

Construct a BPR cost function with the specified α and β values. The BPR cost function is
r = 1 + \alpha \frac{V}{C}^\beta, where r is the ratio of congested to freeflow travel time,
V is volume, C is capacity, and α and β are specified.
"""
struct BPR <: DelayFunc
    α::Float64
    β::Float64
end

get_delay(f::BPR, _, v, c) = 1 + f.α * (v / c) ^ f.β

abstract type CapacityCalculator end

"""
    PerLaneCapacityByTypeCalculator()

    Takes a dict of per-lane capacities by road type, and a default per lane capacity.
"""
struct PerLaneCapacityByTypeCalculator <: CapacityCalculator
    values::Dict{RoadClass.T, Int64}
    default::Int64
end

function get_capacity(c::PerLaneCapacityByTypeCalculator, e::EdgeData)
    lanes = coalesce(e.lanes, 1)

    if haskey(c.values, e.this_class)
        c.values[e.this_class] * lanes
    else
        c.default * lanes
    end
end

end