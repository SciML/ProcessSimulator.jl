module EntropyScalingExt

using EntropyScaling
using ProcessSimulator

const PS = ProcessSimulator


function PS.viscosity(model::M, p, T, z) where M <: EntropyScaling.EoSModel
    return EntropyScaling.viscosity(model, p, T, z)
end


end