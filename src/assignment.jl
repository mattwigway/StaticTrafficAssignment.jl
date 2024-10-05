# load trips onto the network

MIN_λ = 1e-5

# NB b/c we are using a turn-based graph, "segment" refers to vertices herein

struct AssignmentState
    current_segment_flows::Vector{Float64}
    current_turn_flows::Vector{Float64}
    all_or_nothing_segment_flows::Vector{Float64}
    all_or_nothing_turn_flows::Vector{Float64}
end

AssignmentState(G::FWGraph) = AssignmentState(
    zeros(Float64, nv(G.G)),
    zeros(Float64, ne(G.G)),
    zeros(Float64, nv(G.G)),
    zeros(Float64, ne(G.G))
)

# G is graph, odmat is o-d matrix, dest_offset is the index of the first origin/destination, costfunc is a cost function
# this is a bit of an odd function, as it will mutate the graph costs, but also returns the vertex (i.e. segment, remember, edge
# /turn based graph) level flows

"""
    assign_frankwolfe!(G, odflows, costfunc; maxiter=100_000, rel_gap_tol=0.0001)

Assign flows represented by odflows to the network. ODFlows should be a vector of tuples (Source VertexID, [(Target, VertexID, flow)...]). See
create_odflows for help creating it.

Cost func should be a cost function, and capacity func should be a capacity function
"""
function assign_frankwolfe!(G, odflows, costfunc::VDF.DelayFunc, capacityfunc::VDF.CapacityCalculator; maxiter=100_000, rel_gap_tol = 0.0001,
        centroid_collector_speed_kmh=40, iteration_callback=nothing)
    state = AssignmentState(G)

    @info "calculating freeflow speeds"
    compute_initial_weights!(G, centroid_collector_speed_kmh)

    for iter in 1:maxiter
        @info "assignment: begin iteration $iter"
        @info "..all-or-nothing assignment"
        # In Ortuzar and Willumsen, they say to just compute the relative gap after updating the flows and
        # costs, but gloss over that you need all or nothing flows given current costs to be the lower bound
        # so do the all or nothing assignment here, and save it for the next iteration
        tm = @elapsed all_or_nothing!(G, odflows, state)
        @info "....all-or-nothing completed in $(round(tm, digits=3)) seconds"
        
        λ = if iter == 1
            1.0
        else
            find_optimal_λ(G, state, iter, capacityfunc, costfunc)
        end

        @info "optimal λ = $λ"

        update_flows!(state, λ)

        # From AEquilibraE documentation:
        # The relative gap is computed with the cost used to compute the All-or-Nothing portion of the iteration, and although the literature on this is obvious,
        # we took some time to realize that we should re-compute the travel costs only AFTER checking for convergence
        rel_gap = if iter == 1
            1.0
        else
            compute_rel_gap(G, state)
        end

        update_costs!(G, state, capacityfunc, costfunc)

        if !isnothing(iteration_callback)
            iteration_callback(iter, state, λ)
        end

        if rel_gap <= rel_gap_tol
            @info "assignment CONVERGED after $(iter + 1) iterations! relative gap $rel_gap <= $rel_gap_tol"
            break
        elseif iter == maxiter
            @error "assignment FAILED TO CONVERGE in maximum $maxiter iterations! final relative gap: $rel_gap"
            break
        elseif iter > 1
            @info "relative gap: $rel_gap > $rel_gap_tol"
            continue
        end
    end

    return state.current_segment_flows, state.current_turn_flows
end

function all_or_nothing!(G::FWGraph, odflows, state::AssignmentState)
    # perform all-or-nothing assignment with current costs
    # split-apply-combine - allocate an array per thread and sum them up at the end
    # todo could allocate these per-thread flows once, but let's not optimize prematurely
    # threads is the second dimension for ease of summation below
    # this is a little more complicated than a stardard Frank-Wolfe algorithm b/c we have to track segment/vertex
    # flows which are used to calculate costs for each turn, as well as turn/edge flows because they are
    # used in the relative gap metric
    @info "....using $(Threads.nthreads()) threads for all-or-nothing assignment"

    fill!(state.all_or_nothing_segment_flows, 0)
    fill!(state.all_or_nothing_turn_flows, 0)

    # Most of the time is spent in Dijkstra's. Lock when writing per-segment counts
    lck = ReentrantLock()

    Threads.@threads for (origin, destinations) in odflows        
        # run a Dijkstra search on the graph
        # nb this or the enumeration is creating excessive GC. rewrite to save less data (e.g. record
        # flow during path tracing, don't trave all paths and then record flow
        djstate = dijkstra_shortest_paths(G.G, [code_for(G.G, origin)], allpaths=true)

        @debug "Routed origin $origin"

        lck = ReentrantLock()

        # indexed not by vertex index, but by vertex position in dest_vertices
        # this is extremely similar to enumerate paths in
        # https://github.com/JuliaGraphs/LightGraphs.jl/blob/master/src/shortestpaths/bellman-ford.jl#L113
        # but with fewer allocations
        lock(lck) do
            for (destvx, n_trips) in destinations
                if destvx == origin
                    continue
                end

                current_vertex = code_for(G.G, destvx)

                @assert djstate.parents[current_vertex] != 0 && djstate.parents[current_vertex] != current_vertex """
                No path from origin $origin to destination $destvx
                $n_trips trips cannot be completed
                consider using --single-strong-component
                """

                while (djstate.parents[current_vertex] != 0 && djstate.parents[current_vertex] != current_vertex)
                    parent = djstate.parents[current_vertex]
                    eidx = G.G[label_for(G.G, parent), label_for(G.G, current_vertex)].eidx                        
                    state.all_or_nothing_turn_flows[eidx] += n_trips
                    state.all_or_nothing_segment_flows[parent] += n_trips
                    current_vertex = parent 
                end
            end
        end

        @debug "finished origin $origin"
    end
end

# https://sboyles.github.io/teaching/ce392c/5-beckmannmsafw.pdf is a very useful reference for this
function find_optimal_λ(G, state, iter, capacityfunc, costfunc)
    opt_res = optimize(λ -> aggregate_cost(G, state, λ, capacityfunc, costfunc) ^ 2, 1. / iter, 1, Brent())
    @assert converged(opt_res)
    minimizer(opt_res)
end

# return the FW objective function value for the current and all-or-nothing flows in state, and the given λ
function aggregate_cost(G, state, λ, capacityfunc, costfunc)
    obj = 0.0
    for (src, tgt) in edge_labels(G.G)
        edg = G.G[src, tgt]
        eidx = edg.eidx

        # x̂ in presentation
        aon_flow = state.all_or_nothing_segment_flows[code_for(G.G, src)]
        current_flow = state.current_segment_flows[code_for(G.G, src)]
        current_segment_flow = λ * aon_flow + (1 - λ) * current_flow
        capacity = VDF.get_capacity(capacityfunc, edg)
        congested_to_ff_ratio = #if edg.this_class != RoadClass.centroid_connector
            VDF.get_delay(costfunc, edg, current_segment_flow, capacity)
        # else
        #     # no delay on centroid connectors
        #     1.0
        # end
        total_cost = edg.freeflow_traversal_time_secs * congested_to_ff_ratio + edg.turn_cost_secs

        last_turn_flow = state.current_turn_flows[eidx]
        aon_turn_flow = state.all_or_nothing_turn_flows[eidx]
        obj += total_cost * (aon_turn_flow - last_turn_flow)
    end

    isfinite(obj) || error("F-W objective not finite")

    return obj
end

function update_flows!(state::AssignmentState, λ::Float64)

    # first update the flows
    # need to do this before costs since costs depend on flows not only on the same link
    # okay to update segment and turn flows separately. as long as they sum at the start, they will still
    # sum after a linear transformation
    for i in eachindex(state.current_segment_flows)
        state.current_segment_flows[i] = λ * state.all_or_nothing_segment_flows[i] + (1 - λ) * state.current_segment_flows[i]
    end

    for i in eachindex(state.current_turn_flows)
        state.current_turn_flows[i] = λ * state.all_or_nothing_turn_flows[i] + (1 - λ) * state.current_turn_flows[i]
    end
end

function update_costs!(G::FWGraph, state::AssignmentState, capacityfunc, costfunc)
    for (src, tgt) in edge_labels(G.G)
        edg = G.G[src, tgt]
        #if edg.this_class != RoadClass.centroid_connector
            # centroid connectors do not become congested
            # This actually doesn't work, I think because centroid connector flow gets applied to parent above!
            # So it affects the travel time on the other turns, but then is not costed in appropriately.
            # x̂ in presentation
            # already updated, above
            current_segment_flow = state.current_segment_flows[code_for(G.G, src)]
            capacity = VDF.get_capacity(capacityfunc, edg)
            congested_to_ff_ratio = VDF.get_delay(costfunc, edg, current_segment_flow, capacity)
            total_cost = edg.freeflow_traversal_time_secs * congested_to_ff_ratio + edg.turn_cost_secs

            # 0 cost could be okay, if the segment has no length (data issue), but negative costs are definitely a no-no
            @assert total_cost >= 0.0 "Total cost was $total_cost at $src -> $tgt (volume $current_segment_flow, capacity $capacity)"

            # update the edge costs and flows
            G.G[src, tgt] = (
                G.G[src, tgt]...,
                weight=total_cost
            )
        #end
    end
end

function compute_rel_gap(G, state)
    agg_current_flows::Float64 = 0.0
    agg_aon_flows::Float64 = 0.0

    for (src, tgt) in edge_labels(G.G)
        edat = G.G[src, tgt]
        current_turn_flow = state.current_turn_flows[edat.eidx]
        aon_turn_flow = state.all_or_nothing_turn_flows[edat.eidx]
        turn_cost = edat.weight

        agg_current_flows += current_turn_flow * turn_cost
        agg_aon_flows += aon_turn_flow * turn_cost
    end

    @assert agg_current_flows ≥ agg_aon_flows "current flow ($agg_current_flows) less than AON flow ($agg_aon_flows)"

    return (agg_current_flows - agg_aon_flows) / agg_current_flows
end