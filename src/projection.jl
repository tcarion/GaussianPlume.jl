
struct PointSource
    lon::Real
    lat::Real
    # TODO : use the same wind structures as with ATP45
    # wind::Union{WindAzimuth, Wind...}
    params::GaussianPlumeParams
end

"""
    concentration(x, y, z, source::PointSource)
Return the final longitude latitude point and the concentration a this location given a `source`.
"""
function concentration(x, y, z, source::PointSource)
    conc = concentration(x, y, z, source.params)

    # TODO : get the azimuth from the pointsource struct
    azimuth = 45
    downwind = horizontal_walk(source.lon, source.lat, x, azimuth)
    lateral_direction = 111
    final = horizontal_walk(downwind, y, azimuth, lateral_direction)
    [final, conc]
end