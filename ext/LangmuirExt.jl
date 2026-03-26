module LangmuirExt

using Langmuir
using ProcessSimulator

const PS = ProcessSimulator


function PS.loading(model::M, p, T, z) where M <: Langmuir.MultiComponentIsothermModel 
    return Langmuir.loading(model, p, T, z)
end

function PS.loading(model::M, p, T) where M <: Union{Langmuir.IsothermModel,Langmuir.MultiComponentIsothermModel}
    return Langmuir.loading(model, p, T)
end

function PS.isosteric_heat(model::M, p, T, z) where M <: Langmuir.MultiComponentIsothermModel 
    return Langmuir.isosteric_heat(model, p, T, z)
end

function PS.isosteric_heat(model::M, p, T) where M <: Union{Langmuir.IsothermModel,Langmuir.MultiComponentIsothermModel}
    return Langmuir.isosteric_heat(model, p, T)
end


end