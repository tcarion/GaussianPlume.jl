using GaussianPlume
using Test

@testset "Parameters" begin
    @test pasquill_gifford(Moderate, 5.5) == Stabilities(:C, :D)
    @test GaussianPlume.getcoeff(Urban, GaussianPlume.D) == (
        (a = 0.16, b = 0.0004, c = -0.5),
        (a = 0.14, b = 0.0003, c = -0.5)
    )
end

# Example from reference book "Air Dispersion Modeling - Foundations and Applications, Alex De Visscher, p.24-25"
@testset "concentration" begin
    # Example 2.1
    h = 75 + 15
    stabs = pasquill_gifford(Moderate, u)
    relpar = ReleaseParams(h = h, Q = 100, u = 7)
    params = GaussianPlumeParams(release = relpar)
    params.stabilities = stabs
    cground1 = concentration(1500, 0, 0, params) * 1e6
    ccross = concentration(1500, 100, 0, params) * 1e6

    @test cground1 ≈ 160.3 atol=1e-1
    @test ccross ≈ 107.5 atol=1e-1

    h = 30
    relpar = ReleaseParams(h = h, Q = 5, u = 2)

    params = GaussianPlumeParams(release = relpar)
    
    params.stabilities = Stabilities(:D)
    params.reflection = true

    cground2 = concentration(1000, 0, 0, params)
    @test cground2 ≈ 201.139895390322 * 1e-6
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