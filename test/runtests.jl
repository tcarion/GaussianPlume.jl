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

    plume = GaussianPlumeParams(release = relpar)
    
    plume.stabilities = Stabilities(:D)
    plume.reflection = true

    cground = plume(1000, 0, 0)
    # caxes = concentration(1000, 0, h)
    @test cground ≈ 201.139895390322 * 1e-6
end

@testset "concentration - multiple classes" begin
    plume = GaussianPlumeParams()
    plume.stabilities = Stabilities(:A, :B)
    cboth= plume(1000, 0, 0)

    plume.stabilities = Stabilities(:A)
    ca= plume(1000, 0, 0)

    plume.stabilities = Stabilities(:B)
    cb= plume(1000, 0, 0)
end