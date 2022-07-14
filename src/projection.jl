abstract type AbstractWind end
mutable struct WindAzimuth <: AbstractWind
    speed::Real
    azimuth::Real
end

struct PointSource
    lon::Real
    lat::Real
    wind::WindAzimuth
    params::GaussianPlumeParams
end

"""
    concentration(x, y, z, source::PointSource)
Return the final longitude latitude point and the concentration a this location given a `source`.
"""
function concentration(x, y, z, source::PointSource)
    conc = concentration(x, y, z, source.params)
    azimuth = source.wind.azimuth
    downwind = horizontal_walk(source.lon, source.lat, x, azimuth)
    lateral_direction = azimuth + 90
    final = horizontal_walk(downwind, y, azimuth, lateral_direction)
    [final, conc]
end