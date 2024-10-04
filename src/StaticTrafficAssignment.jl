module StaticTrafficAssignment
import Compat: @compat
import LibSpatialIndex
import Graphs: DiGraph, ne, nv, dijkstra_shortest_paths, strongly_connected_components, rem_vertex!
import MetaGraphsNext: MetaGraph, labels, edge_labels, code_for, label_for
import EnumX: @enumx
import Geodesy: LatLon, euclidean_distance
import DataStructures: counter, inc!, DefaultDict
import OpenStreetMapPBF: scan_pbf
import Logging: @warn, @error, @info
import Missings: passmissing
import Serialization: serialize
import CodecZlib: GzipCompressorStream
import CSV
import DataFrames: DataFrame, groupby
import ForwardDiff

include("fwgraph.jl")
include("compute_heading.jl")
include("build_fw_graph.jl")
include("index.jl")
include("centroid_connectors.jl")
include("costfuncs.jl")
include("fw_weights.jl")
include("assignment.jl")
include("flows_from_dataframe.jl")

export build_graph, create_centroid_connectors!, flows_from_dataframe, assign_frankwolfe!, VDF
end