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

function copy_state!(dst::AssignmentState, src::AssignmentState)::AssignmentState
    copy!(dst.current_segment_flows, src.current_segment_flows)
    copy!(dst.current_turn_flows, src.current_turn_flows)
    copy!(dst.all_or_nothing_segment_flows, src.all_or_nothing_segment_flows)
    copy!(dst.all_or_nothing_turn_flows, src.all_or_nothing_turn_flows)
    return dst
end

# G is graph, odmat is o-d matrix, dest_offset is the index of the first origin/destination, costfunc is a cost function
# this is a bit of an odd function, as it will mutate the graph costs, but also returns the vertex (i.e. segment, remember, edge
# /turn based graph) level flows

"""
    assign_frankwolfe!(G, odflows, costfunc; maxiter=100_000, rel_gap_tol=0.0001)

Assign flows represented by odflows to the network. ODFlows should be a vector of tuples (Source VertexID, [(Target, VertexID, flow)...]). See
create_odflows for help creating it.

Cost func should be a cost function, and capacity func should be a capacity function
"""
function assign_frankwolfe!(G, odflows, costfunc::VDF.DelayFunc, capacityfunc::VDF.CapacityCalculator; maxiter=100_000, rel_gap_tol = 0.0001, centroid_collector_speed_kmh=40)
    @info "calculating freeflow speeds"
    compute_initial_weights!(G, centroid_collector_speed_kmh)

    state = AssignmentState(G)

    @info "initial all-or-nothing assignment"
    # in Ortuzar and Willumsen they don't have an initial assignment, but they also have to do the all
    # or nothing assignment twice, implicitly, once explicitly and once to compute relgap
    tm = @elapsed all_or_nothing!(G, odflows, state)
    @info "....all-or-nothing completed in $(round(tm, digits=3)) seconds"

    # And in order to make sure that all flows are preserved, I believe we need to set current flows to the
    # all-or-nothing assignments (rather than zero), because otherwise every successive step will still have a
    # little of the initial zero mixed in, and the total flows won't sum up.
    copy!(state.current_segment_flows, state.all_or_nothing_segment_flows)
    copy!(state.current_turn_flows, state.all_or_nothing_turn_flows)

    λ = 0.5

    for iter in 1:maxiter
        @info "assignment: begin iteration $iter"   

        λ = find_optimal_λ(G, state, λ, capacityfunc, costfunc)

        @info "optimal λ = $λ"

        @info "updating flows and costs"
        update_flows_and_costs!(G, λ, state, capacityfunc, costfunc)

        @info "..all-or-nothing assignment"
        # In Ortuzar and Willumsen, they say to just compute the relative gap after updating the flows and
        # costs, but gloss over that you need all or nothing flows given current costs to be the lower bound
        # so do the all or nothing assignment here, and save it for the next iteration
        tm = @elapsed all_or_nothing!(G, odflows, state)
        @info "....all-or-nothing completed in $(round(tm, digits=3)) seconds"

        rel_gap = compute_rel_gap(G, state)

        if rel_gap <= rel_gap_tol
            @info "assignment CONVERGED after $(iter + 1) iterations! relative gap $rel_gap <= $rel_gap_tol"
            break
        elseif iter == maxiter
            @error "assignment FAILED TO CONVERGE in maximum $maxiter iterations! final relative gap: $rel_gap"
            break
        else
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
                    if eidx > 0
                        # don't record centroid connectors (eidx 0)
                        
                        state.all_or_nothing_turn_flows[eidx] += n_trips
                        state.all_or_nothing_segment_flows[parent] += n_trips
                    end
                    current_vertex = parent 
                end
            end
        end

        @debug "finished origin $origin"
    end
end

# https://sboyles.github.io/teaching/ce392c/5-beckmannmsafw.pdf is a very useful reference for this
function find_optimal_λ(G, state, λ, capacityfunc, costfunc)
    # newton's method
    prev_λ = -1
    while true
        objval = aggregate_cost(G, state, λ, capacityfunc, costfunc)
        
        if abs(objval) < 1e-5 || abs(λ - prev_λ) < 1e-5
            return λ
        end
        
        drv_objval = ForwardDiff.derivative(x -> aggregate_cost(G, state, x, capacityfunc, costfunc), λ)

        prev_λ = λ
        λ = clamp(λ - objval / drv_objval, MIN_λ, 1)
    end
end

# return the FW objective function value for the current and all-or-nothing flows in state, and the given λ
function aggregate_cost(G, state, λ, capacityfunc, costfunc)
    obj = 0.0
    for (src, tgt) in edge_labels(G.G)
        edg = G.G[src, tgt]
        # x̂ in presentation
        current_segment_flow = λ * state.all_or_nothing_segment_flows[code_for(G.G, src)] + (1 - λ) * state.current_segment_flows[code_for(G.G, tgt)]
        capacity = VDF.get_capacity(capacityfunc, edg)
        congested_to_ff_ratio = VDF.get_delay(costfunc, edg, current_segment_flow, capacity)
        total_cost = edg.freeflow_traversal_time_secs * congested_to_ff_ratio + edg.turn_cost_secs

        eidx = edg.eidx
        if eidx > 0
            # don't include centroid connectors
            last_turn_flow = state.current_turn_flows[eidx]
            aon_turn_flow = state.all_or_nothing_turn_flows[eidx]
            obj += total_cost * (aon_turn_flow - last_turn_flow)
        end
    end

    return obj
end

function update_flows_and_costs!(G::FWGraph, λ::Float64, state::AssignmentState, capacityfunc, costfunc)

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

    for (src, tgt) in edge_labels(G.G)
        edg = G[src, tgt]
        # x̂ in presentation
        current_segment_flow = λ * state.all_or_nothing_segment_flows[code_for(G.G, src)] + (1 - λ) * state.current_segment_flows[code_for(G.G, tgt)]
        capacity = VDF.get_capacity(capacityfunc, edg)
        congested_to_ff_ratio = VDF.get_delay(costfunc, edg, current_segment_flow, capacity)
        total_cost = edg.freeflow_traversal_time_secs * congested_to_ff_ratio + edge.turn_cost_secs

        # 0 cost could be okay, if the segment has no length (data issue), but negative costs are definitely a no-no
        @assert total_cost >= 0.0

        # update the edge costs and flows
        G[src, tgt] = (
            G[src, tgt]...,
            weight=total_cost
        )
    end
end

function compute_rel_gap(G, state)
    agg_current_flows::Float64 = 0.0
    agg_aon_flows::Float64 = 0.0

    for (src, tgt) in edge_labels(G)
        edat = G[src, tgt].eidx
        if edat.eidx > 0
            # don't worry about centroid connectors
            current_turn_flow = state.current_turn_flows[edat.eidx]
            aon_turn_flow = state.all_or_nothing_turn_flows[edat.eidx]
            turn_cost = edat.weight

            agg_current_flows += current_turn_flow * turn_cost
            agg_aon_flows += aon_turn_flow * turn_cost
        end
    end

    return (agg_current_flows - agg_aon_flows) / agg_current_flows
end