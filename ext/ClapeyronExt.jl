module ClapeyronExt

import ProcessSimulator
import Clapeyron 
import Symbolics


function EosBasedGuesses(EoSModel::M, p::V, T::V, z::D) where {M <: Clapeyron.ActivityModel, V <: Real, D <: AbstractArray{ <: Real}}

    ## Bubble and dew pressure
    pᵇ = bubble_pressure(EoSModel, T, z)[1]
    pᵈ = dew_pressure(EoSModel, T, z)[1]
    pᵇᵈ = [pᵇ, pᵈ]

    Nc = length(z)
    nᵢⱼ = zeros(V, Nc, 2)
    x = zeros(V, Nc, 3)
    ρ = zeros(V, 3)
    h = zeros(V, 3)

    if p ≤ pᵇ && p ≥ pᵈ
        sol = TP_flash(EoSModel, p, T, z)
        nᵢⱼ .= sol[1]
        x .= sol[2]

        ## Enthalpy
        hₗ = enthalpy(EoSModel, p, T, x[:, 2], phase = "liquid")
        hᵥ = enthalpy(EoSModel, p, T, x[:, 3], phase = "vapor")
        ϕ = sum(nᵢⱼ, dims = 1)/sum(nᵢⱼ)
        hₒᵥ = hₗ*ϕ[1] + hᵥ*ϕ[2] 
        h .= [hₒᵥ, hₗ, hᵥ]

    elseif p > pᵇ
        nᵢⱼ[:, 1] .= z
        x[:, 1:2] .= z
        ϕ = sum(nᵢⱼ, dims = 1)/sum(nᵢⱼ)
        ρ[2] = PT_molar_density(EoSModel, p, T, z, phase = "liquid")
        ρ[3] = PT_molar_density(EoSModel, p, T, z, phase = "vapor")
        ρ[1] = 1.0/(ϕ[1]/ρ[2] + ϕ[2]/ρ[3])
        hₗ = enthalpy(EoSModel, p, T, z, phase = :liquid)
        hᵥ = enthalpy(EoSModel, p, T, z, phase = :vapor)
        hₒᵥ = hₗ*ϕ[1] + hᵥ*ϕ[2] 
        h .= [hₒᵥ, hₗ, hᵥ]

    elseif p < pᵈ
        nᵢⱼ[:, 2] .= z
        x[:, [1, 3]] .= z
        ϕ = sum(nᵢⱼ, dims = 1)/sum(nᵢⱼ)
        ρ[2] = PT_molar_density(EoSModel, p, T, z, phase = "liquid")
        ρ[3] = PT_molar_density(EoSModel, p, T, z, phase = "vapor")
        ρ[1] = 1.0/(ϕ[1]/ρ[2] + ϕ[2]/ρ[3])
        hₗ = enthalpy(EoSModel, p, T, z, phase = :liquid)
        hᵥ = enthalpy(EoSModel, p, T, z, phase = :vapor)
        hₒᵥ = hₗ*ϕ[1] + hᵥ*ϕ[2] 
        h .= [hₒᵥ, hₗ, hᵥ]
    end


    return EosBasedGuesses(EoSModel, p, T, ρ, x, h, pᵇᵈ)
end

function PT_molar_density(EoSModel::M, p, T, x; phase = "unknown") where M <: Clapeyron.EoSModel
    Clapeyron.molar_density(EoSModel, p, T, x, phase = phase)   
end

#Assumes only two phases
function TP_flash(EoSModel::M, p, T, x) where M <: Clapeyron.EoSModel

         if two_phase_check(EoSModel, p, T, x)
            
            sol = Clapeyron.tp_flash(EoSModel, p, T, x, RRTPFlash())
            xᵢⱼ = transpose(sol[1]) |> Array
            nᵢⱼ = transpose(sol[2]) |> Array
            ϕ = sum(nᵢⱼ, dims = 1)/sum(nᵢⱼ)

         elseif vapor_check(EoSModel, p, T, x)

            xᵢⱼ = [ones(length(x))/length(x) x]
            ϕ = [0.0, 1.0 - 1e-7]

         else
                
            xᵢⱼ = [x ones(length(x))/length(x)]
            ϕ = [1.0 - 1e-7, 0.0]

         end

         return (ϕ, [x xᵢⱼ])  
end

function two_phase_check(EoSModel::M, p, T, x) where M <: Clapeyron.EoSModel
    istwophase = (p ≥ dewP(EoSModel, T, x) + 1e-3) && (p ≤ bubbleP(EoSModel, T, x) - 1e-3) ? true : false
    return istwophase
end

function vapor_check(EoSModel::M, p, T, x) where M <: Clapeyron.EoSModel
    isvapor = (p < dewP(EoSModel, T, x)) ? true : false
    return isvapor
end

function flash_mol_fractions(EoSModel::M, p, T, x) where M <: Clapeyron.EoSModel
    z = TP_flash(EoSModel, p, T, x)[2]
    return z
end

function flash_mol_fractions_liquid(EoSModel::M, p, T, x) where M <: Clapeyron.EoSModel
         z = flash_mol_fractions(EoSModel, p, T, x)[:, 2] 
         return z
end

function flash_mol_fractions_vapor(EoSModel::M, p, T, x) where M <: Clapeyron.EoSModel
    z = flash_mol_fractions(EoSModel, p, T, x)[:, 3] 
    return z
end

function flash_vaporized_fraction(EoSModel::M, p, T, x) where M <: Clapeyron.EoSModel
         ϕ = TP_flash(EoSModel, p, T, x)[1] 
         return ϕ
end

function bubbleP(EoSModel::M, T, x) where M <: Clapeyron.EoSModel

    if length(EoSModel.components) ≥ 2

        Clapeyron.bubble_pressure(EoSModel, T, x)[1]
    else

        Clapeyron.saturation_pressure(EoSModel, T)[1]
    end
end

function dewP(EoSModel::M, T, x) where M <: Clapeyron.EoSModel

    if length(EoSModel.components) ≥ 2

        Clapeyron.dew_pressure(EoSModel, T, x)[1]
    else

        Clapeyron.saturation_pressure(EoSModel, T)[1]
    end
end

function ρT_enthalpy(EoSModel::M, ρ, T, x) where M <: Clapeyron.EoSModel
    return Clapeyron.VT_enthalpy(EoSModel, 1.0/ρ, T, x)
end

function ρT_internal_energy(EoSModel::M, ρ, T, x) where M <: Clapeyron.EoSModel
    return Clapeyron.VT_internal_energy(EoSModel, 1.0/ρ, T, x)
end

function sat_temperature(EoSModel::M, p) where M <: Clapeyron.EoSModel
    return Clapeyron.saturation_temperature(EoSModel, p)
end

Symbolics.@register_array_symbolic TP_flash(model::Clapeyron.EoSModel, p, T, arr::AbstractVector) begin
    size = (2,)
    eltype = eltype(arr)
end

Symbolics.@register_array_symbolic flash_mol_fractions_liquid(model::Clapeyron.EoSModel, p, T, arr::AbstractVector) begin
    size = (length(arr), )
    eltype = eltype(arr)
end

Symbolics.@register_array_symbolic flash_mol_fractions_vapor(model::Clapeyron.EoSModel, p, T, arr::AbstractVector) begin
    size = (length(arr), )
    eltype = eltype(arr)
end

Symbolics.@register_array_symbolic flash_vaporized_fraction(model::Clapeyron.EoSModel, p, T, arr::AbstractVector) begin
    size = (2,)
    eltype = eltype(arr)
end

@register_symbolic sat_temperature(model::Clapeyron.EoSModel, p)

@register_symbolic ρT_enthalpy(model::Clapeyron.EoSModel, ρ, T, arr::AbstractVector)

@register_symbolic ρT_internal_energy(model::Clapeyron.EoSModel, ρ, T, arr::AbstractVector)

@register_symbolic bubbleP(model::Clapeyron.EoSModel, T, arr::AbstractVector)

@register_symbolic dewP(model::Clapeyron.EoSModel, T, arr::AbstractVector)

end


