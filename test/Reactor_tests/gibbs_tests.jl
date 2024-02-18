using ProcessSimulator
using Test

y = Gibbs(2.0)

@test y == 4.0