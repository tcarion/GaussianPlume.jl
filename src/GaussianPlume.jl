module GaussianPlume

using Proj4

export ReleaseParams, DispersionParams, GaussianPlumeParams
export concentration, StabilityClass, Terrain, PasquillGiffordCriteria, pasquill_gifford
export Rural, Urban
export A, B, C, D, E, F
export Strong, Moderate, Slight, Cloudy, Clear

const DEFAULT_PROJ = Ref{Any}(C_NULL)

function __init__()
    DEFAULT_PROJ[] = Projection("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
end

include("helpers.jl")
include("concentration.jl")
include("projection.jl")

end
