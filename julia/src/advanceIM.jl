
# IMEX scheme: explicit nonlinear advection + implicit viscous diffusion.
# utm, vtm: du/dt (explicit part) from previous step.
# Returns new ut, vt.
function advance_im!(t::Float64, dt::Float64,
                     un::Matrix{Float64}, vn::Matrix{Float64}, pn::Matrix{Float64},
                     unp1::Matrix{Float64}, vnp1::Matrix{Float64}, pnp1::Matrix{Float64},
                     utm::Matrix{Float64}, vtm::Matrix{Float64}, par::Par)
    dx = par.dx; dy = par.dy; nu = par.nu
    tnp1 = t + dt
    I1, I2 = get_index(par.gid)

    # Viscous terms at time n (reused in both predictor and corrector)
    uxx = (un[I1.+1,I2] .- 2.0.*un[I1,I2] .+ un[I1.-1,I2]) ./ dx^2
    uyy = (un[I1,I2.+1] .- 2.0.*un[I1,I2] .+ un[I1,I2.-1]) ./ dy^2
    vxx = (vn[I1.+1,I2] .- 2.0.*vn[I1,I2] .+ vn[I1.-1,I2]) ./ dx^2
    vyy = (vn[I1,I2.+1] .- 2.0.*vn[I1,I2] .+ vn[I1,I2.-1]) ./ dy^2

    # ===== Predictor: explicit advection, semi-implicit viscosity =====
    ut, vt = get_ut(t, un, vn, pn, 0.0, par)   # nu_scale=0: no viscous terms in explicit RHS
    unp1[I1,I2] .= un[I1,I2] .+ par.ab1.*ut[I1,I2] .+ par.ab2.*utm[I1,I2] .+ (0.5*dt*nu).*(uxx .+ uyy)
    vnp1[I1,I2] .= vn[I1,I2] .+ par.ab1.*vt[I1,I2] .+ par.ab2.*vtm[I1,I2] .+ (0.5*dt*nu).*(vxx .+ vyy)
    solve_implicit_time_step!(unp1, vnp1, tnp1, par)
    apply_boundary_conditions!(unp1, vnp1, tnp1, par)
    pnp1 .= pressure_equation!(tnp1, unp1, vnp1, dt, false, par)

    # ===== Corrector: Crank-Nicolson =====
    utp, vtp = get_ut(tnp1, unp1, vnp1, pnp1, 0.0, par)
    unp1[I1,I2] .= un[I1,I2] .+ (0.5*dt).*(utp[I1,I2] .+ ut[I1,I2]) .+ (0.5*dt*nu).*(uxx .+ uyy)
    vnp1[I1,I2] .= vn[I1,I2] .+ (0.5*dt).*(vtp[I1,I2] .+ vt[I1,I2]) .+ (0.5*dt*nu).*(vxx .+ vyy)
    solve_implicit_time_step!(unp1, vnp1, tnp1, par)
    apply_boundary_conditions!(unp1, vnp1, tnp1, par)
    pnp1 .= pressure_equation!(tnp1, unp1, vnp1, dt, false, par)

    return ut, vt
end
