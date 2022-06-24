"""
    horizontal_walk(lon::AbstractFloat, lat::AbstractFloat, distance::AbstractFloat, azimuth::AbstractFloat)

Compute the end location given a starting location `lon` and `lat` in degrees, a distance `distance` in meters
and an azimuth `azimuth` in degrees (the reference direction is North)

julia> horizontal_walk(50., 4., 11100., 90.)
2-element Vector{Float64}:
 50.09995485698309
  3.9999938922206835
"""
function horizontal_walk(start::Vector{<:AbstractFloat}, distance::AbstractFloat, azimuth::AbstractFloat)
    proj = DEFAULT_PROJ[]
    q2, _ = geod_direct(lonlat2xy(start, proj), azimuth, distance, proj)
    xy2lonlat(q2, proj)
end
horizontal_walk(lon::AbstractFloat, lat::AbstractFloat, distance::AbstractFloat, azimuth::AbstractFloat) = horizontal_walk([lon, lat], distance, azimuth)

