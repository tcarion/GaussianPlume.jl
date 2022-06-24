@enum Terrain Rural Urban
@enum StabilityClass A B C D E F

"""
    PasquillGiffordCriteria
Critera for determining the stability class according to sky conditions (Pasquill, 1961; Gifford, 1961)
`Strong`, `Moderate` and `Slight` are for incoming solar radiation (for daytime)
`Cloudy` and `Clear` are for cloudiness (for nighttime)
"""
@enum PasquillGiffordCriteria Strong Moderate Slight Cloudy Clear

"""
    pasquill_gifford(criteria::PasquillGiffordCriteria, windspeed::Real)
Return a list of stability classes according to the Pasquill Gifford criteria

# Example
julia> pasquill_gifford(Moderate, 5.5)
Set{StabilityClass} with 2 elements:
 C
 D
"""
function pasquill_gifford(criteria::PasquillGiffordCriteria, windspeed::Real)
    # TODO (see Table 2.1 of reference book)
    # using DelimitedFiles ??
    Set([A, B])
end

"""
    getcoeff(terrain::Terrain, stability::StabilityClass)

Given the `terrain` and the `stability`, return the three coefficients a, b, c for σ_y and σ_z respectively.

# Example
```jl
julia> getcoeff(Rural, A)
((a = 0.22, b = 0.0001, c = -0.5), (a = 0.2, b = 0, c = -0.5))
```
"""
function getcoeff(terrain::Terrain, stability::StabilityClass)
    # TODO
    # Get the a, b c coefficients from the Briggs_Dispersion_Parameters.xls file (maybe convert it to CSV)
    (
        (a = 0.22, b = 0.0001, c = -0.5), # coefficients for σ_y
        (a = 0.2, b = 0, c = -0.5) # coefficients for σ_z
    )
end

"""
    disp_function(a, b, c) 

Return the function for calculating the dispersion parameters with respect to x given the coefficients of the equation ax(1 + bx)^c
```
"""
disp_function(a, b, c) = x -> a * x * (1 + b*x)^c

"""
    ReleaseParams
Structure related to the release conditions
"""
Base.@kwdef mutable struct ReleaseParams
    # Effective source height [m]
    h::Real = 1.
    # Emission rate [g/s]
    Q::Real = 1.
    # Wind speed [m/s]
    u::Real = 5.
end

Base.@kwdef mutable struct DispersionParams
    # Horizontal dispersion parameter [m]
    σ_y::Real = 2.15
    # Vertical dispersion parameter [m]
    σ_z::Real = 2.15
end

"""
    DispersionParams(x::Real, terrain::Terrain, stability::StabilityClass)  

Return a `DispersionParams` struct given the downwind direction `x`, the terrain and stability class
# Example
julia> DispersionParams(100, Rural, A)
DispersionParams(21.89081818461976, 20.0)
```
"""
function DispersionParams(x::Real, terrain::Terrain, stability::StabilityClass)
    ycoef, zcoef = getcoeff(terrain, stability)
    σ_y = disp_function(ycoef...)(x)
    σ_z = disp_function(zcoef...)(x)
    DispersionParams(σ_y, σ_z)
end
Base.:+(d1::DispersionParams, d2::DispersionParams) = DispersionParams(d1.σ_y + d2.σ_y, d1.σ_z + d2.σ_z) 
Base.:/(d1::DispersionParams, x::Real) = DispersionParams(d1.σ_y / x, d1.σ_z / x) 

"""
    Params
"""
Base.@kwdef mutable struct GaussianPlumeParams
    # release parameters
    release::ReleaseParams = ReleaseParams()
    # release parameters
    terrain::Terrain = Rural
    # set of stability classes
    stabilities::Set{<:StabilityClass} = Set([A])
    # if ground reflection is considered
    reflection::Bool = true
end
GaussianPlumeParams(release::ReleaseParams, terrain::Terrain, criteria::PasquillGiffordCriteria) = GaussianPlumeParams(release, terrain, pasquill_gifford(criteria, release.u))

_yterm(y, sig) = exp(-0.5 * (y / sig)^2)
_zterm_noreflect(z, σ_z, h) = exp(-0.5 * ((z - h) / σ_z)^2)
_zterm_reflect(z, σ_z, h) = exp(-0.5 * ((z - h) / σ_z)^2) +  exp(-0.5 * ((z + h) / σ_z)^2)
function concentration(y, z, disp::DispersionParams = DispersionParams(), release::ReleaseParams = ReleaseParams(); reflection = true)
    h, Q, u = _destruct(release)
    σ_y = disp.σ_y
    σ_z = disp.σ_z
    coef = Q / (2pi * u * σ_y * σ_z)
    yt = _yterm(y, σ_y)
    zt = reflection ? _zterm_reflect(z, σ_z, h) : _zterm_noreflect(z, σ_z, h)
    coef * yt * zt
end
# concentration(ys::AbstractVector, z::Real, release::ReleaseParams; reflection = false) = concentration.(ys, z, Ref(release), reflection = reflection)
# concentration(y::Real, zs::AbstractVector, release::ReleaseParams; reflection = false) = concentration.(y, zs, Ref(release), reflection = reflection)
# concentration(ys::AbstractVector, zs::AbstractVector, release::ReleaseParams; reflection = false) = [concentration(y, z, release, reflection = reflection) for y in ys, z in zs]

function concentration(x::Real, y::Real, z::Real, release::ReleaseParams, terrain::Terrain, stabilities::Set{<:StabilityClass}; reflection = true)
    disps = DispersionParams.(x, terrain, stabilities)
    # If more than 1 stability, we take the average of the dispersion parameters for each class
    disp = sum(disps) / length(disps)
    concentration(y, z, disp, release, reflection = reflection)
end

"""
    concentration(x::Real, y::Real, z::Real, params::Params)
Return the concentration in [g m^{-3}] at some point given the conditions stated in `params`
`x` is the downwind distance from the point source.
`y` is the horizontal distance perpendicular to the wind.
`z` is the vertical distance from the point source.
"""
concentration(x::Real, y::Real, z::Real, params::GaussianPlumeParams) = 
    concentration(x, y, z, params.release, params.terrain, params.stabilities; reflection = params.reflection)

function _destruct(p::ReleaseParams)
    ntuple(fieldcount(ReleaseParams)) do i
        getfield(p, i)
    end
end
    