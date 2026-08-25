
# AB2 + AM2 predictor-corrector 2nd order time step.
# utm, vtm: du/dt from previous step (in).
# Returns new ut, vt (du/dt at current t, needed by next step).
function advance_pc!(t::Float64, dt::Float64,
                     un::Matrix{Float64}, vn::Matrix{Float64}, pn::Matrix{Float64},
                     unp1::Matrix{Float64}, vnp1::Matrix{Float64}, pnp1::Matrix{Float64},
                     utm::Matrix{Float64}, vtm::Matrix{Float64}, par::Par)
    tnp1 = t + dt
    I1, I2 = get_index(par.gid)

    # ===== AB2 Predictor =====
    ut, vt = get_ut(t, un, vn, pn, par.nuScaleFactor, par)
    unp1[I1,I2] .= un[I1,I2] .+ par.ab1.*ut[I1,I2] .+ par.ab2.*utm[I1,I2]
    vnp1[I1,I2] .= vn[I1,I2] .+ par.ab1.*vt[I1,I2] .+ par.ab2.*vtm[I1,I2]
    apply_boundary_conditions!(unp1, vnp1, tnp1, par)
    pnp1 .= pressure_equation!(tnp1, unp1, vnp1, dt, false, par)

    # ===== AM2 Corrector =====
    utp, vtp = get_ut(tnp1, unp1, vnp1, pnp1, par.nuScaleFactor, par)
    unp1[I1,I2] .= un[I1,I2] .+ (0.5*dt).*(utp[I1,I2] .+ ut[I1,I2])
    vnp1[I1,I2] .= vn[I1,I2] .+ (0.5*dt).*(vtp[I1,I2] .+ vt[I1,I2])
    apply_boundary_conditions!(unp1, vnp1, tnp1, par)
    pnp1 .= pressure_equation!(tnp1, unp1, vnp1, dt, false, par)

    return ut, vt
end
