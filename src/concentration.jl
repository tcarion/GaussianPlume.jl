@enum Terrain Rural Urban
@enum StabilityClass A B C D E F

const StabilityClasses = AbstractSet{StabilityClass}

"""
    $(TYPEDSIGNATURES)
Define the stability classes for the model according to Pasquill and Gifford (Pasquill, 1961; Gifford, 1961).
If multiple classes are defined, the average result for each class is considered.

# Example
```jl
julia> Stabilities(:A, :B)
Set{StabilityClass} with 2 elements:
  GaussianPlume.A
  GaussianPlume.B
```
"""
function Stabilities(stabs::Vararg{Union{String, Symbol}}) :: StabilityClasses
    Set(tostab(stabs...))
end

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
    if windspeed < 2
        if criteria == Strong
            Set([A])
        elseif criteria == Moderate
            Set([A, B])
        elseif criteria == Slight
            Set([B])
        elseif criteria == Cloudy
            Set([E])
        elseif criteria == Clear
            Set([F])
        end
    elseif 2 <= windspeed < 3
        if criteria == Strong
            Set([A, B])
        elseif criteria == Moderate
            Set([B])
        elseif criteria == Slight
            Set([C])
        elseif criteria == Cloudy
            Set([E])
        elseif criteria == Clear
            Set([F])
        end
    elseif 3 <= windspeed < 5
        if criteria == Strong
            Set([B])
        elseif criteria == Moderate
            Set([B, C])
        elseif criteria == Slight
            Set([C])
        elseif criteria == Cloudy
            Set([D])
        elseif criteria == Clear
            Set([E])
        end
    elseif 5 <= windspeed <= 6
        if criteria == Strong
            Set([C])
        elseif criteria == Moderate
            Set([C, D])
        elseif criteria == Slight
            Set([D])
        elseif criteria == Cloudy
            Set([D])
        elseif criteria == Clear
            Set([D])
        end
    else
        if criteria == Strong
            Set([C])
        elseif criteria == Moderate
            Set([D])
        elseif criteria == Slight
            Set([D])
        elseif criteria == Cloudy
            Set([D])
        elseif criteria == Clear
            Set([D])
        end
    end
end

"""
    getcoeff(terrain::Terrain, stability::StabilityClass)

Given the `terrain` and the `stability`, return the three coefficients a, b, c for ??_y and ??_z respectively.

# Example
```jl
julia> getcoeff(Rural, A)
((a = 0.22, b = 0.0001, c = -0.5), (a = 0.2, b = 0, c = 1.0))
```
"""
function getcoeff(terrain::Terrain, stability::StabilityClass)
    if terrain == Rural
        if stability == A
            sigycol = 1
            sigzcol = 7
        elseif stability == B
            sigycol = 2
            sigzcol = 8
        elseif stability == C
            sigycol = 3
            sigzcol = 9
        elseif stability == D 
            sigycol = 4
            sigzcol = 10
        elseif stability == E
            sigycol = 5
            sigzcol = 11
        elseif stability == F
            sigycol = 6
            sigzcol = 12
        end
    elseif terrain == Urban
        if stability == A || stability == B
            sigycol = 13
            sigzcol = 17
        elseif stability == C
            sigycol = 14
            sigzcol = 18
        elseif stability == D
            sigycol = 15
            sigzcol = 19
        elseif stability == E || stability == F
            sigycol = 16
            sigzcol = 20
        end
    end
    nkeys = (:a, :b, :c)
    ((; zip(nkeys, BRIGGS_COEFS[:, sigycol])...), (; zip(nkeys, BRIGGS_COEFS[:, sigzcol])...))
end

"""
    disp_function(a, b, c) 

Return the function for calculating the dispersion parameters with respect to x given the coefficients of the equation ax(1 + bx)^c
```
"""
disp_function(a, b, c) = x -> a * x * (1 + b*x)^c

"""
    $(TYPEDEF)
Structure related to the release conditions

    $(FIELDS)
"""
Base.@kwdef mutable struct ReleaseParams
    "Effective source height [m]"
    h::Real = 1.
    "Emission rate [g/s]"
    Q::Real = 1.
    "Wind speed [m/s]"
    u::Real = 5.
end

Base.@kwdef mutable struct DispersionParams
    "Horizontal dispersion parameter [m]"
    y::Real = 2.15
    "Vertical dispersion parameter [m]"
    z::Real = 2.15
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
    sigma_y = disp_function(ycoef...)(x)
    sigma_z = disp_function(zcoef...)(x)
    DispersionParams(sigma_y, sigma_z)
end
Base.:+(d1::DispersionParams, d2::DispersionParams) = DispersionParams(d1.y + d2.y, d1.z + d2.z) 
Base.:/(d1::DispersionParams, x::Real) = DispersionParams(d1.y / x, d1.z / x) 

"""
    $(TYPEDEF)
Structure related to the release conditions

    $(FIELDS)
"""
Base.@kwdef mutable struct GaussianPlumeParams
    "release parameters"
    release::ReleaseParams = ReleaseParams()
    "terrain (Urban/Rural)"
    terrain::Terrain = Rural
    "set of stability classes"
    stabilities::StabilityClasses = Set([A])
    "if ground reflection is considered"
    reflection::Bool = true
    "height of the mixing layer"
    hmix::Union{Real, Nothing} = nothing
end
GaussianPlumeParams(release::ReleaseParams, terrain::Terrain, criteria::PasquillGiffordCriteria) = GaussianPlumeParams(release, terrain, pasquill_gifford(criteria, release.u))

"""
    $(TYPEDSIGNATURES)
Return the concentration in [g m^{-3}] at some point given the conditions stated in `params`.
- `x` is the downwind distance from the point source.
- `y` is the horizontal distance perpendicular to the wind.
- `z` is the vertical distance from the point source.
"""
(plume::GaussianPlumeParams)(x::Real, y::Real, z::Real) = concentration(x, y, z, plume.release, plume.terrain, plume.stabilities; reflection = plume.reflection)

"""
    $(SIGNATURES)
Buoyancy flux parameter according to Briggs (1968).
"""
buoyancy_flux(gas_density, air_density, stack_radius, gas_velocity) = (1 - gas_density / air_density) * GRAVITY * stack_radius^2 * gas_velocity

"""
    $(SIGNATURES)
Plume rise `??h` [m] according to Briggs (1975) parametrization. `x` is the downwind distance in meter and `u` the downwind velocity in meter/second.
"""
function plume_rise(flux_param, x, u) 
    if flux_param <= 55.
        xf = 49. * flux_param^(5/8)
    else
        xf = 119. * flux_param^(2/5)
    end
    1.6 * flux_param^(1/3) * min(x, xf)^(2/3) / u
end

downwind_mass(Q, u) = Q / u
gaussian_exp(x, p1, p2) = exp(-0.5 * ((x + p1) / p2)^2)
crosswind_conc(y, sigma_y) = 1 / sqrt(2pi) / sigma_y * gaussian_exp(y, 0, sigma_y)
vertical_conc(z, sigma_z, h) = 1 / sqrt(2pi) / sigma_z * sum(gaussian_exp.(z, h, sigma_z))
function concentration(y, z, disp::DispersionParams = DispersionParams(), release::ReleaseParams = ReleaseParams(); reflection = true)
    h, Q, u = _destruct(release)
    sigma_y = disp.y
    sigma_z = disp.z
    hs = reflection ? [h, -h] : -h
    downwind_mass(Q, u) * crosswind_conc(y, sigma_y) * vertical_conc(z, sigma_z, hs)
end
# concentration(ys::AbstractVector, z::Real, release::ReleaseParams; reflection = false) = concentration.(ys, z, Ref(release), reflection = reflection)
# concentration(y::Real, zs::AbstractVector, release::ReleaseParams; reflection = false) = concentration.(y, zs, Ref(release), reflection = reflection)
# concentration(ys::AbstractVector, zs::AbstractVector, release::ReleaseParams; reflection = false) = [concentration(y, z, release, reflection = reflection) for y in ys, z in zs]

function concentration(x::Real, y::Real, z::Real, release::ReleaseParams, terrain::Terrain, stabilities::StabilityClasses; reflection = true)
    disps = DispersionParams.(x, terrain, stabilities)
    # If more than 1 stability, we take the average of the dispersion parameters for each class
    disp = sum(disps) / length(disps)
    concentration(y, z, disp, release, reflection = reflection)
end

# """
#     $(TYPEDSIGNATURES)
# Return the concentration in [g m^{-3}] at some point given the conditions stated in `params`.
# - `x` is the downwind distance from the point source.
# - `y` is the horizontal distance perpendicular to the wind.
# - `z` is the vertical distance from the point source.
# """
# concentration(x::Real, y::Real, z::Real, params::GaussianPlumeParams) = 
#     concentration(x, y, z, params.release, params.terrain, params.stabilities; reflection = params.reflection)

function _destruct(p::ReleaseParams)
    ntuple(fieldcount(ReleaseParams)) do i
        getfield(p, i)
    end
end

function tostab(stabs::Vararg{Symbol})
    [getproperty(@__MODULE__, s) for s in collect(stabs)]
end

function tostab(stabs::Vararg{AbstractString})
    sstabs = Symbol.(stabs)
    tostab(sstabs...)
end