using GaussianPlume
using Test

release = ReleaseParams()
params = GaussianPlumeParams()
terrain = params.terrain
stabilities = params.stabilities
x = 100.
@testset "GaussianPlume.jl" begin
    @test pasquill_gifford(Moderate, 5.5) == Set([C, D])
    @test GaussianPlume.getcoeff(Urban, D) == (
        (a = 0.16, b = 0.0004, c = -0.5),
        (a = 0.14, b = 0.0003, c = -0.5)
    )
    concentration(x, 1, 2, params)
end
