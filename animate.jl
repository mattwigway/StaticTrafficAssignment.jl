# Create an animation of the convergence process

using CairoMakie, Serialization, ArgParse, GeoDataFrames, StaticTrafficAssignment, ArchGDAL, DataFrames,
    ProgressMeter, MetaGraphsNext, Graphs
import GeoFormatTypes as GFT

const COST_FUNC = VDF.BPR(0.6, 7)
const CAPACITY_FUNC = VDF.PerLaneCapacityByTypeCalculator(Dict(), 2100)
const COLORMAP = Reverse(:BrBG_4)

function get_congestion(G, flows)
    map(zip(1:nv(G.G), flows)) do (code, flow)
        # we need to get an edge (any edge) to calculate capacity and VDF
        src = label_for(G.G, code)
        tgt = first(outneighbor_labels(G.G, src))
        edg = G.G[src, tgt]

        capacity = VDF.get_capacity(CAPACITY_FUNC, edg)
        VDF.get_delay(COST_FUNC, edg, flow, capacity) 
    end
end

function plot_flows(G, geoms, flows, max_flow, figarea, title)
    linewidth=max.(flows[1:nrow(geoms)] ./ max_flow .* 16, 1)

    ax = Axis(figarea, aspect=DataAspect(), title=title, xticklabelsvisible=false, yticklabelsvisible=false, xgridvisible=false, ygridvisible=false)

    congestion = get_congestion(G, flows)[1:nrow(geoms)]

    # draw most congested last
    sorter = sortperm(congestion)

    for (g, w, c) in zip(geoms.geom[sorter], linewidth[sorter], congestion[sorter])
        lines!(ax, [g], linewidth=w, color=c, linecap=:round, colormap=COLORMAP, colorrange=(1, 10), highclip=:red)
    end
end

function draw_frame(G, geoms, states, λs, relgaps, iter, max_flow, outdir)
    fig = Figure(size=(3180, 1900), fontsize=fontsize=50)

    # show previous current assignment, and corresponding all-or-nothing assignment
    if iter > 1
        plot_flows(G, geoms, states[iter - 1].state.current_segment_flows, max_flow, fig[1:5, 1], "Current")
    end
    plot_flows(G, geoms, states[iter].state.all_or_nothing_segment_flows, max_flow, fig[1:5, 2], "All-or-nothing")
    Colorbar(fig[6, 1:2], limits=(1, 10), colormap=COLORMAP, highclip=:red, label="Congested to freeflow ratio", vertical=false)

    axλ = Axis(fig[7, 1], yscale=log10, xlabel="Iteration", ylabel="λ")
    ylims!(axλ, (1e-6, 1))
    axrg = Axis(fig[7, 2], yscale=log10, xlabel="Iteration", ylabel="Rel. gap")
    ylims!(axrg, (1e-5, 1))

    if iter > 1
        lines!(axλ, 1:iter, λs[1:iter])
        lines!(axrg, 1:iter, relgaps[1:iter])
    end

    save(joinpath(outdir, "frame_$(lpad(iter - 1, 5, "0")).png"), fig)
end

function read_geoms(graph, projection)
    # TODO do these still line up after island removal?!?!?!?!
    raw_geom = deserialize(graph * ".geoms")

    gdf = DataFrame(
        geom = map(raw_geom) do nodes
            ArchGDAL.createlinestring(map(nc -> [nc.lon, nc.lat], nodes))
        end
    )

    metadata!(gdf, "geometrycolumns", (:geom,))

    gdf.geom = reproject(gdf.geom, GFT.EPSG(4326), GFT.EPSG(projection), order=:trad)

    return gdf
end

# function read_states(iteration_dir, states_to_read)
#     states = map(states_to_read) do i
#         filename = joinpath(iteration_dir, "iter_$i.state")
#         push!(states, deserialize(filename))
#         i += 1
#     end
#     return states
# end


function read_states(iteration_dir)
    states_to_read = filter(x -> occursin(r"^iter_[0-9]+.state", x), readdir(iteration_dir))
    states = map(states_to_read) do i
        filename = joinpath(iteration_dir, i)
        s = deserialize(filename)
        (s..., iter=parse(Int64, match(r"[0-9]+", filename).match))
    end

    sort!(states, by=x->x.iter)

    return states
end

function main()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "graph"
            help = "graph file"
        "iteration_dir"
            help = "directory with iterations"
        "output_dir"
            help = "directory to output slides"
        "--projection"
            help = "EPSG code number for projection to use"
            arg_type = Int
            default = 3857 # web mercator
        "--range"
            help = "Frames to render, as a number, a Julia range syntax (or several separated by commas)"
    end

    args = parse_args(s)

    G::StaticTrafficAssignment.FWGraph = deserialize(args["graph"])
    geoms = read_geoms(args["graph"], args["projection"])

    range = if isnothing(args["range"])
        1:length(filter(x -> occursin(r"^iter_[0-9]+.state", x), readdir(args["iteration_dir"])))
    else
        vcat(eval(Meta.parse(args["range"]))...)
    end

    states = read_states(args["iteration_dir"])

    println("Read $(length(states)) states")

    CairoMakie.activate!()

    max_flow = maximum(states[1].state.all_or_nothing_segment_flows)

    λs = [state.λ for state in states]
    relgaps = [state.rg for state in states]

    @showprogress for i in range
        draw_frame(G, geoms, states, λs, relgaps, i, max_flow, args["output_dir"])
    end
end

main()