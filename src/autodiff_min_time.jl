function autodiff_min_time(s, thr, ceff, μ, T)
    # Inputs
    x = @views SA[s[1:6]...]
    m = s[7] 
    λ = @views SA[s[8:13]...]

    # State dynamics
    _, B = mee_dynamics(x, μ)
    û = -normalize(B'*λ) # optimal control (primer vector)
    ṁ = -thr/ceff
    f = @closure x -> ode_min_time(x, m, û, thr, μ)
    # Costate dynamics 
    if T === "ForwardDiff"
        λ̇ = -ForwardDiff.jacobian(f, x)'*λ
    elseif T == "Enzyme"
        λ̇ = -Enzyme.jacobian(Forward, f, x)[1]'*λ
    end

    λ̇ₘ = -thr/m^2*norm(B'*λ)

    # Return derivatives
    return SA[f(x)..., ṁ, λ̇..., λ̇ₘ]
end


function ode_min_time(x, m, û, thr, μ) 
    A, B = mee_dynamics(x, μ) 
    A + thr/m*B*û
end




