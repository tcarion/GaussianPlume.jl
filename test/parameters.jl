using GaussianPlume
using GaussianPlume: HeatBalanceParams, pasquill_gifford
using GaussianPlume: Aermod, Calpuff
using Test

@test pasquill_gifford(Moderate, 5.5) == Stabilities(:C, :D)
@test GaussianPlume.briggs_dispersion(Urban, GaussianPlume.D) == (
    (a = 0.16, b = 0.0004, c = -0.5),
    (a = 0.14, b = 0.0003, c = -0.5)
    )

@test HeatBalanceParams(Aermod) isa HeatBalanceParams
@test HeatBalanceParams(C_g = 0.3) isa HeatBalanceParams