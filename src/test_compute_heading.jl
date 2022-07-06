using Geodesy

include("compute_heading.jl")

lat0 = parse(Float64, ARGS[1])
lon0 = parse(Float64, ARGS[2])
lat1 = parse(Float64, ARGS[3])
lon1 = parse(Float64, ARGS[4])

heading = compute_heading(LLA(lat0, lon0, 0.0), LLA(lat1, lon1, 0.0))

ns0 = lat0 > 0 ? "N" : "S"
ew0 = lon0 > 0 ? "E" : "W"

ns1 = lat1 > 0 ? "N" : "S"
ew1 = lon1 > 0 ? "E" : "W"

@info "Bearing from $(abs(lat0)) $ns0, $(abs(lon0)) $ew0 -> $(abs(lat1)) $ns1, $(abs(lon1)) $ew1: $heading"