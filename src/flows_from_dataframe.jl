"""
    flows_from_dataframe(df, src, tgt, count, connectors)

Create an O-D flow dataset from a dataframe, with columns for source TAZ, target TAZ, the number of trips, and a connectors
dict from create_centroid_connectors!
"""
function flows_from_dataframe(df, src, tgt, count, connectors)
    result = Tuple{VertexID, Vector{Tuple{VertexID, Float64}}}[]

    for grp in groupby(df, src)
        destinations = map(zip(grp[!, tgt], grp[!, count])) do (dest, n)
            (connectors[dest], n)
        end

        origin = connectors[first(grp[!, src])]

        push!(result, (origin, destinations))
    end
    
    return result
end