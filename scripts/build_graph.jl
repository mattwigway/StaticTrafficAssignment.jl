import StaticTrafficAssignment
import ArgParse: ArgParseSettings, @add_arg_table
import ArgParse

function parse_args()
    s = ArgParseSettings()
    @add_arg_table s begin
        "osm_pbf"
            help = "An OpenStreetMap .pbf file to process"
        "network"
            help = "An output file to write the network to (extension .fwgr)"
        "--save-names"
            help = "Save a sidecar file with street names"
            action = :store_true
        "--save-geometries"
            help = "Save a sidecar file with way geometries"
            action = :store_true
    end
    return ArgParse.parse_args(s)
end

function main()
    args = parse_args()
    pbf = args["osm_pbf"]::String
    outf = args["network"]::String
    save_names = args["save-names"]::Bool
    save_geom = args["save-geometries"]::Bool

    StaticTrafficAssignment.build_graph(pbf, outf, save_names=save_names, save_geoms=save_geom)
end

main()