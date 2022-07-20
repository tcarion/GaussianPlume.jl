using GaussianPlume
using Test

@testset "Parameters" begin
    @test pasquill_gifford(Moderate, 5.5) == Stabilities(:C, :D)
    @test GaussianPlume.getcoeff(Urban, GaussianPlume.D) == (
        (a = 0.16, b = 0.0004, c = -0.5),
        (a = 0.14, b = 0.0003, c = -0.5)
    )
end

# Example from reference book "Air Dispersion Modeling - Foundations and Applications, Alex De Visscher, p.25"
@testset "concentration" begin
    h = 30
    relpar = ReleaseParams(h = h, Q = 5, u = 2)

    params = GaussianPlumeParams(release = relpar)
    
    params.stabilities = Stabilities(:D)
    params.reflection = true

    cground = concentration(1000, 0, 0, params)
    # caxes = concentration(1000, 0, h, params)
    @test cground ≈ 201.139895390322 * 1e-6
end

# Example from reference book "Air Dispersion Modeling - Foundations and Applications, Alex De Visscher, p.31"
@testset "Plume rise" begin
    Q = 20 # [m^3/s]
    ρₐ = 1.17
    ρₛ = 0.935
    rₛ = 1.
    wₛ = Q / (π*rₛ^2)
    u = 3

    Fb = GaussianPlume.buoyancy_flux(ρₛ, ρₐ, rₛ, wₛ)
    @test Fb ≈ 12.539 atol=1e-3
    
    Δh = GaussianPlume.plume_rise(Fb, 1000., u)
    @test Δh ≈ 47.589 atol=1e-3
end

@testset "concentration - multiple classes" begin
    params = GaussianPlumeParams()
    params.stabilities = Stabilities(:A, :B)
    cboth= concentration(1000, 0, 0, params)

    params.stabilities = Stabilities(:A)
    ca= concentration(1000, 0, 0, params)

    params.stabilities = Stabilities(:B)
    cb= concentration(1000, 0, 0, params)
end