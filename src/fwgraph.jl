struct VertexID
    id::Int64
end

@enumx RoadClass motorway motorway_link trunk trunk_link primary primary_link secondary secondary_link tertiary tertiary_link unclassified centroid_connector

function get_road_class(cln::String)
    if cln == "motorway"
        return RoadClass.motorway
    elseif cln == "motorway_link"
        return RoadClass.motorway_link
    elseif cln == "trunk"
        return RoadClass.trunk
    elseif cln == "trunk_link"
        return RoadClass.trunk_link
    elseif cln == "primary"
        return RoadClass.primary
    elseif cln == "primary_link"
        return RoadClass.primary_link
    elseif cln == "secondary"
        return RoadClass.secondary
    elseif cln == "secondary_link"
        return RoadClass.secondary_link
    elseif cln == "tertiary"
        return RoadClass.tertiary
    elseif cln == "tertiary_link"
        return RoadClass.tertiary_link
    elseif cln == "unclassified"
        return RoadClass.unclassified
    elseif cln == "centroid_connector"
        return RoadClass.centroid_connector
    else
        error("unknown highway class $cln")
    end
end

EdgeData = @NamedTuple begin
    length_m::UInt16
    this_class::RoadClass.T
    next_class::RoadClass.T
    turn_angle::Int16
    traffic_signal::UInt8
    speed_kmh::Union{Missing, UInt8}
    lanes::Union{Missing, UInt8}
    oneway::Bool
    weight::Float64
    freeflow_traversal_time_secs::Float64
    turn_cost_secs::Float64

    # single numeric ID for this edge
    eidx::Int64
end

VertexData = @NamedTuple begin
    geom::LatLon{Float64}
end

mutable struct FWGraph
    G::MetaGraph{Int64, DiGraph{Int64}, VertexID, VertexData, EdgeData}
    spidx::Union{LibSpatialIndex.RTree, Nothing}
    next_centroid_connector::Int64
    center_lat::Float64
end

function FWGraph()
    FWGraph(MetaGraph(
        DiGraph();
        label_type=VertexID,
        edge_data_type=EdgeData,
        vertex_data_type=VertexData,
        weight_function=ed -> ed.weight
    ), nothing, -1, 0)
end