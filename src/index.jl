"""
    build_spatial_index!(G::FWGraph)

Build a spatial index for the vertices of graph G, and assign it in-place
"""
function build_spatial_index!(G::FWGraph)
    G.spidx = LibSpatialIndex.RTree(2)

    G.center_lat = find_center_lat(G)

    for v in labels(G.G)
        if v.id < 0
            # don't index centroid connectors
            continue
        end

        geom = G.G[v].geom
        # scale lon to be comparable to lat
        LibSpatialIndex.insert!(G.spidx, v.id, [geom.lon * cosd(G.center_lat), geom.lat])
    end
end

"""
   find_center_lat(G): find the latitude at the center of the graph
"""
find_center_lat(G) = median([G.G[node].geom.lat for node in labels(G.G)])