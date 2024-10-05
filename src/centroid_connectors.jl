
"""
    create_centroid_connectors!(G, ids, locations)

Create centroid connectors for locations (which should be a Vector of LatLons). Return a dict
mapping a centroid ID (from ids, parallel to locations) to VertexIDs in the graph.

It is safe to use this function more than once on the same graph.
"""
function create_centroid_connectors!(G, ids, locations::AbstractVector{<:LatLon{<:Any}})
    isnothing(G.spidx) && build_spatial_index!(G)

    connectors = Dict{eltype(ids), VertexID}()

    for (id, location) in zip(ids, locations)
        candidates = LibSpatialIndex.knn(G.spidx, [location.lon * cosd(G.center_lat), location.lat], 15)

        if isempty(candidates)
            @warn "TAZ $id at $location not linked"
            continue
        end

        distances = map(candidates) do c
            return euclidean_distance(location, G.G[VertexID(c)].geom)
        end

        thres_dist = minimum(distances) * 1.1

        if (thres_dist > typemax(UInt16))
            @warn "TAZ $id at $location is more than $(typemax(UInt16))m from nearest road, clamping centroid connector length at $(typemax(UInt16))m"
        end

        vid = VertexID(G.next_centroid_connector)
        G.next_centroid_connector -= 1 # centroid connectors have negative indices
        G.G[vid] = (geom=location,)
        connectors[id] = vid

        linked = false
        for (candidate, distance) in zip(candidates, distances)
            if distance < thres_dist
                linked = true

                # Centroid connectors are kinda weird here–while the rest of the graph is turn-based,
                # centroid connectors are not—just a vertex with edges to each of the nearby ways
                G.G[vid, VertexID(candidate)] = (
                    length_m=round(UInt16, min(distance, typemax(UInt16))),
                    this_class=RoadClass.centroid_connector,
                    next_class=RoadClass.centroid_connector, # TODO not correct
                    turn_angle=zero(Int16),
                    traffic_signal=zero(UInt8),
                    speed_kmh=missing,
                    lanes=missing,
                    oneway=false,
                    weight=NaN,
                    freeflow_traversal_time_secs=NaN,
                    turn_cost_secs=NaN,
                    eidx=0
                )

                # Cannot re-use edge data, edge ID is different
                G.G[VertexID(candidate), vid] = (
                    length_m=round(UInt16, min(distance, typemax(UInt16))),
                    this_class=RoadClass.centroid_connector,
                    next_class=RoadClass.centroid_connector, # TODO not correct
                    turn_angle=zero(Int16),
                    traffic_signal=zero(UInt8),
                    speed_kmh=missing,
                    lanes=missing,
                    oneway=false,
                    weight=NaN,
                    freeflow_traversal_time_secs=NaN,
                    turn_cost_secs=NaN,
                    eidx=0
                )
            end
        end

        @assert linked
    end

    renumber_edges!(G)

    return connectors
end


"""
    renumber_edges!(G)

Update edge indices in G to be 1:ne(G)
"""
function renumber_edges!(G)
    for (idx, (src, tgt)) in enumerate(edge_labels(G.G))
        G.G[src, tgt] = (
            G.G[src, tgt]...,
            eidx=idx
        )
    end
end