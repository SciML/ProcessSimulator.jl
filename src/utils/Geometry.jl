abstract type AbstractTank end

mutable struct CylindricalTank{V <: Real} <: AbstractTank
    D::V #Diameter constant
    L::V #Overall height
    H::V #Height of each segment (defaults to 1 segment)
    porosity::V #Porosity
end

function CylindricalTank(;D, L, H = L, porosity = 0.5)
    return CylindricalTank(D, L, H, porosity)
end

function discretize!(tank::CylindricalTank, nseg::Int)
    tank.H = tank.L/nseg
end

function volume_(tank::CylindricalTank)
    return π*(tank.D/2)^2*tank.H
end

function surface_area_(tank::CylindricalTank)
    return 2*π*(tank.D/2)*tank.H + π*(tank.D/2)^2
end

function cross_section_area_(tank::CylindricalTank)
    return π*(tank.D/2)^2
end

export CylindricalTank, volume_, surface_area_, cross_section_area_, discretize!