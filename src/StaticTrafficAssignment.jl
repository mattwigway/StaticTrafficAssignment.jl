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
import Serialization: serialize, deserialize
import CodecZlib: GzipCompressorStream, GzipDecompressorStream
import CSV
import DataFrames: DataFrame, groupby, metadata!, nrow
import Optim: optimize, Brent, minimizer, converged
import GeoDataFrames
import GeoFormatTypes as GFT
import ArchGDAL as AG
import StatsBase: median

include("fwgraph.jl")
include("compute_heading.jl")
include("build_fw_graph.jl")
include("index.jl")
include("centroid_connectors.jl")
include("costfuncs.jl")
include("fw_weights.jl")
include("assignment.jl")
include("flows_from_dataframe.jl")
include("gis.jl")

export build_graph, create_centroid_connectors!, flows_from_dataframe, assign_frankwolfe!, VDF, graph_to_gis
end