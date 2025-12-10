using ModelingToolkit

const R = 8.31446261815324 # J/(mol K)

abstract type AbstractMassTransferModel end

abstract type AbstractViscosityModel end

abstract type AbstractThermalConductivityModel end

abstract type AbstractDiffusivityModel end

abstract type AbstractHeatTransferModel end

struct HomogeneousDiffusivityCoeff{V <: AbstractVector{<: Real}} <: AbstractMassTransferModel
Dh::V
end

struct ConstantHeatTransferCoeff{V <: Real} <: AbstractHeatTransferModel
    k::V # Heat transfer coefficient
end

struct ConstantMassTransferCoeff{V <: AbstractVector{<: Real}} <: AbstractMassTransferModel
    k::V # Mass transfer coefficient
end

struct TransportModel{M <: Union{Nothing, AbstractMassTransferModel},
    H <: Union{Nothing, AbstractHeatTransferModel},
    V,
    K,
    D}

    MassTransferModel::M
    HeatTransferModel::H
    ViscosityModel::V
    ThermalConductivityModel::K
    DiffusivityModel::D
end

function TransportModel(masstransfermodel, heattransfermodel)
    return TransportModel(masstransfermodel, heattransfermodel, nothing, nothing, nothing)
end

function TransportModel(masstransfermodel, heattransfermodel, viscositymodel)
    return TransportModel(masstransfermodel, heattransfermodel, viscositymodel, nothing, nothing)
end

export HomogeneousDiffusivityCoeff, ConstantHeatTransferCoeff, ConstantMassTransferCoeff, TransportModel