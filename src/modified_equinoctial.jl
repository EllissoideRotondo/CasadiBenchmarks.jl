"""
Affine-form of equations of motions in modified equinoctial elements
"""
function mee_dynamics(x::T, μ = 1) where {T <: AbstractVector}
    mee_r2bp(x, μ), 
    mee_pert(x, μ)
end


function mee_r2bp(x::T, μ = 1) where {T <: AbstractVector}
    p, f, g, L = x[1], x[2], x[3], x[6]
    sinL, cosL = sincos(L)
    w = 1 + f * cosL + g * sinL
    SA[0, 0, 0, 0, 0, √(μ * p) * (w / p)^2]
end


function mee_pert(x::T, μ = 1) where {T <: AbstractVector}
    p, f, g, h, k, L = x
    sinL, cosL = sincos(L)
    w = 1 + f * cosL + g * sinL
    s² = 1 + h^2 + k^2
    SMatrix{6, 3}(0, √(p / μ) * sinL, -√(p / μ) * cosL, 0, 0, 0,
        2p / w * √(p / μ), √(p / μ) / w * ((w + 1) * cosL + f), √(p / μ) / w *
        ((w + 1) * sinL + g), 0, 0, 0,
        0, -√(p / μ) * g / w * (h * sinL - k * cosL),
        √(p / μ) * f / w * (h * sinL - k * cosL), √(p / μ) * s² * cosL / 2 / w,
        √(p / μ) * s² * sinL / 2w, √(p / μ) / w * (h * sinL - k * cosL)
    )
end