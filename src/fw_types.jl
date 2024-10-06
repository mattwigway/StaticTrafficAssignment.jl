@enum RoadClass motorway motorway_link trunk trunk_link primary primary_link secondary secondary_link tertiary tertiary_link unclassified centroid_connector

function get_road_class(cln::String)
    if cln == "motorway"
        return motorway
    elseif cln == "motorway_link"
        return motorway_link
    elseif cln == "trunk"
        return trunk
    elseif cln == "trunk_link"
        return trunk_link
    elseif cln == "primary"
        return primary
    elseif cln == "primary_link"
        return primary_link
    elseif cln == "secondary"
        return secondary
    elseif cln == "secondary_link"
        return secondary_link
    elseif cln == "tertiary"
        return tertiary
    elseif cln == "tertiary_link"
        return tertiary_link
    elseif cln == "unclassified"
        return unclassified
    elseif cln == "centroid_connector"
        return centroid_connector
    else
        error("unknown highway class $cln")
    end
end

@enum TurnType left right straight



