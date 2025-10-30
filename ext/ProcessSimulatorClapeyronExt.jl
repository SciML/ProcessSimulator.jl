module ProcessSimulatorClapeyronExt

using ProcessSimulator
using Clapeyron
using Symbolics

const PS = ProcessSimulator

function PS.EosBasedGuesses(EoSModel::M, p::V, T::V, z::D, ::Val{:Pressure}) where {M <: Clapeyron.EoSModel, V <: Real, D <: AbstractArray{ <: Real}}

    sol = PS.TP_flash(EoSModel, p, T, z)
    ϕ = sol[1]
    x = sol[2]
    ρ = zeros(V, 3)
    h = zeros(V, 3)

    ## Enthalpy
    hₗ = Clapeyron.enthalpy(EoSModel, p, T, x[:, 2], phase = "liquid")
    hᵥ = Clapeyron.enthalpy(EoSModel, p, T, x[:, 3], phase = "vapor")
    ρ[2] = PS.PT_molar_density(EoSModel, p, T, x[:, 2], phase = "liquid")
    ρ[3] = PS.PT_molar_density(EoSModel, p, T, x[:, 3], phase = "vapor")
    ρ[1] = 1.0/(ϕ[1]/ρ[2] + ϕ[2]/ρ[3])
    hₒᵥ = hₗ*ϕ[1] + hᵥ*ϕ[2] 
    h .= [hₒᵥ, hₗ, hᵥ]

    return PS.EosBasedGuesses(EoSModel, p, T, ρ, x, h, ϕ)
end

function PS.EosBasedGuesses(EoSModel::M, V::K, T::K, N::D, ::Val{:Volume}) where {M <: Clapeyron.EoSModel, K <: Real, D <: AbstractArray{ <: Real}}
    
    # initialize properties
    ρ = zeros(K, 3)
    h = zeros(K, 3)

    res = Clapeyron.vt_flash(EoSModel, V, T, N)
    nphases = length(res.compositions)
    _x = N./sum(N)
    p = pressure(res)

    if nphases ≤ 1
        v = first(res.volumes)
        z_xᵢⱼ = vcat(_x, _x, _x)
        h .= Clapeyron.VT_enthalpy(EoSModel, v, T, _x)
        ρ .= 1.0/v
        if p*v/(8.314*T) ≥ 0.5 #Very rought test for compressibility factor
            ϕ = [0.0, 1.0]
        else
            ϕ = [1.0, 0.0]
        end
    else
        ϕ = res.fractions/sum(res.fractions)
        xᵢⱼ = mapreduce(x -> x, hcat, res.compositions)
        vl = first(res.volumes)
        vv = last(res.volumes)
        
        #Properties per phase
        h[2] = Clapeyron.VT_enthalpy(EoSModel, vl, T, @view xᵢⱼ[:, 1])
        h[3] = Clapeyron.VT_enthalpy(EoSModel, vv, T, @view xᵢⱼ[:, 2])
        ρ[2] = 1.0/vl
        ρ[3] = 1.0/vv

        ρ[1] = 1.0/(ϕ[1]/ρ[2] + ϕ[2]/ρ[3])
        h[1] = h[2]*ϕ[1] + h[3]*ϕ[2]
        z_xᵢⱼ = hcat(_x, xᵢⱼ)
    end

    return PS.EosBasedGuesses(EoSModel, p, T, ρ, z_xᵢⱼ, h, ϕ)
end

function PS.PT_molar_density(EoSModel::M, p, T, x; phase = "unknown") where M <: Clapeyron.EoSModel
    Clapeyron.molar_density(EoSModel, p, T, x, phase = phase)   
end

function PS.TP_flash(EoSModel::M, p, T, x; nonvolatiles = nothing, noncondensables = nothing) where M <: Clapeyron.EoSModel

    _x = abs.(x)

    #abs(v - vv_ideal) ≤ 1e-3
    #p*v/(8.314*T) ≥ 0.52

    if PS.is_stable(EoSModel, p, T, _x)

        v = Clapeyron.volume(EoSModel, p, T, _x, phase = :unknown)
        #vv_ideal = Clapeyron.volume(Clapeyron.idealmodel(EoSModel), p, T, _x)

        if p*v/(8.314*T) ≥ 0.5 #Very rought test for compressibility factor
            xᵢⱼ = [_x _x]
            ϕ = [0.0, 1.0]

            else

            xᵢⱼ = [_x _x]
            ϕ = [1.0, 0.0]
        end


    else 
        sol = Clapeyron.tp_flash(EoSModel, p, T, _x, MichelsenTPFlash(equilibrium = :vle, nonvolatiles = nonvolatiles, noncondensables = noncondensables))
        xᵢⱼ = transpose(sol[1]) |> Array
        nᵢⱼ = transpose(sol[2]) |> Array
        ϕ = sum(nᵢⱼ, dims = 1)/sum(nᵢⱼ)
    end

        return (ϕ, [_x xᵢⱼ])  
end



function PS.TP_flash(EoSModel::M, p, T, x; nonvolatiles = nothing, noncondensables = nothing) where M <: Clapeyron.CompositeModel

    _x = abs.(x)

    if p > first(Clapeyron.bubble_pressure(EoSModel, T, _x))
        xᵢⱼ = [_x _x]
        ϕ = [1.0, 0.0]
    elseif p ≥ first(Clapeyron.dew_pressure(EoSModel, T, _x))
        sol = Clapeyron.tp_flash(EoSModel, p, T, _x, MichelsenTPFlash(equilibrium = :vle, nonvolatiles = nonvolatiles, noncondensables = noncondensables))
        xᵢⱼ = transpose(sol[1]) |> Array
        nᵢⱼ = transpose(sol[2]) |> Array
        ϕ = sum(nᵢⱼ, dims = 1)/sum(nᵢⱼ)
    else 
        xᵢⱼ = [_x _x]
        ϕ = [0.0, 1.0]
    end

        return (ϕ, [_x xᵢⱼ])  
end

function PS.is_stable(EoSModel::M, p, T, x) where M <: Clapeyron.EoSModel
    return Clapeyron.isstable(EoSModel, p, T, x)
end

function PS.is_stable(EoSModel::M, p, T, x) where M <: Clapeyron.CompositeModel
    pᵦ = Clapeyron.bubble_pressure(EoSModel, T, x)
    pᵢ = Clapeyron.dew_pressure(EoSModel, T, x)
    if p < pᵦ || p > pᵢ
        return false
    else
        return true
    end
end

function PS.is_stable(EoSModel::M, ρ, T, x) where M <: Clapeyron.IdealModel
    return true
end

function PS.is_VT_stable(EoSModel::M, v, T, x) where M <: Clapeyron.EoSModel
        return Clapeyron.VT_isstable(EoSModel, v, T, x)
end

function PS.is_VT_stable(EoSModel::M, v, T, x) where M <: Clapeyron.IdealModel
    return true
end

function PS.flash_mol_fractions(EoSModel::M, p, T, x) where M <: Clapeyron.EoSModel
    z = PS.TP_flash(EoSModel, p, T, x)[2]
    return z
end

function PS.flash_mol_fractions_liquid(EoSModel::M, p, T, x) where M <: Clapeyron.EoSModel
         z = PS.flash_mol_fractions(EoSModel, p, T, x)[:, 2] 
         return z
end

function PS.flash_mol_fractions_vapor(EoSModel::M, p, T, x) where M <: Clapeyron.EoSModel
    z = PS.flash_mol_fractions(EoSModel, p, T, x)[:, 3] 
    return z
end

function PS.flash_vaporized_fraction(EoSModel::M, p, T, x) where M <: Clapeyron.EoSModel
         ϕ = PS.TP_flash(EoSModel, p, T, x)[1] 
         return ϕ
end

function PS.ρT_enthalpy(EoSModel::M, ρ, T, x) where M <: Clapeyron.EoSModel
    return Clapeyron.VT_enthalpy(EoSModel, 1.0/ρ, T, x)
end

function PS.pT_enthalpy(EoSModel::M, p, T, x, phase = :unknown) where M <: Clapeyron.EoSModel
    return Clapeyron.enthalpy(EoSModel, p, T, x, phase = phase)
end

function PS.ρT_internal_energy(EoSModel::M, ρ, T, x) where M <: Clapeyron.EoSModel
    return Clapeyron.VT_internal_energy(EoSModel, 1.0/ρ, T, x)
end

function PS.molecular_weight(model::M, z::AbstractVector) where M <: Clapeyron.EoSModel
    return Clapeyron.molecular_weight(model, z)
end


Symbolics.@register_array_symbolic PS.TP_flash(model::Clapeyron.EoSModel, p, T, arr::AbstractVector) begin
    size = (2,)
    eltype = eltype(arr)
end

Symbolics.@register_array_symbolic PS.flash_mol_fractions_liquid(model::Clapeyron.EoSModel, p, T, arr::AbstractVector) begin
    size = (length(arr), )
    eltype = eltype(arr)
end

Symbolics.@register_array_symbolic PS.flash_mol_fractions_vapor(model::Clapeyron.EoSModel, p, T, arr::AbstractVector) begin
    size = (length(arr), )
    eltype = eltype(arr)
end

Symbolics.@register_array_symbolic PS.flash_vaporized_fraction(model::Clapeyron.EoSModel, p, T, arr::AbstractVector) begin
    size = (2,)
    eltype = eltype(arr)
end

@register_symbolic PS.ρT_enthalpy(model::Clapeyron.EoSModel, ρ, T, arr::AbstractVector)

@register_symbolic PS.ρT_internal_energy(model::Clapeyron.EoSModel, ρ, T, arr::AbstractVector)

@register_symbolic PS.pT_enthalpy(model::Clapeyron.EoSModel, p, T, arr::AbstractVector, sym::Union{Symbol, String})


end


