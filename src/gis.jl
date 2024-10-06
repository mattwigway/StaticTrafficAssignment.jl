"""
    graph_to_gis(output_filename, G, geom_file, names_file=nothing)

Given a graph and the sidecar geometry file (and, optionally, the sidecar names file), save a GIS
file (any format supported by GeoDataFrames.jl is fair game). Centroid connectors will be included if they 
exist in the graph.
"""
function graph_to_gis(output_filename, G::FWGraph, geom_file, names_file=nothing)
    # build the basic dataframe
    # read the geometries
    gdf = DataFrame(geom = map(deserialize(geom_file)) do nodes::Vector{NodeAndCoord}
            AG.createlinestring(map(nc -> [nc.lon, nc.lat], nodes))
    end)

    gdf.centroid_connector .= false

    # add centroid connectors
    ccs = DataFrame(geom=[
        AG.createlinestring([[G.G[src].geom.lon, G.G[src].geom.lat], [G.G[tgt].geom.lon, G.G[tgt].geom.lat]])
        for (src, tgt) in edge_labels(G.G)
        if G.G[src, tgt].this_class == RoadClass.centroid_connector
    ])

    ccs.centroid_connector .= true

    if !isnothing(names_file)
        gdf.name = open(names_file) do stream
            readlines(GzipDecompressorStream(stream))
        end

        ccs.name .= "centroid_connector"
    end

    gdf = vcat(gdf, ccs)

    metadata!(gdf, "geometrycolumns", (:geom,))

    GeoDataFrames.write(output_filename, gdf; crs=GFT.EPSG(4326))
end