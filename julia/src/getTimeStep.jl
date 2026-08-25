
function get_time_step(step::Int, un::Matrix{Float64}, vn::Matrix{Float64}, par::Par)
    dx = par.dx;  dy = par.dy
    nu = par.nu

    ts = par.ts
    if ts == "fe"
        aStab = -2.0;  bStab = 1.0
    elseif ts == "rk2"
        aStab = -2.0;  bStab = 1.5
    elseif ts == "rk4"
        aStab = -2.7;  bStab = 2.7
    elseif ts == "ab2"
        aStab = -1.0;  bStab = 0.8
    elseif ts == "pc2"
        aStab = -1.75; bStab = 1.3
    elseif ts == "im2"
        aStab = -1e20; bStab = 1.7
    else
        error("get_time_step: unknown ts=$ts")
    end

    I1, I2 = get_index(par.gid)

    reLambda = 4.0*nu*(1.0/dx^2 + 1.0/dy^2)
    imLambda = maximum(abs.(un[I1,I2])./dx .+ abs.(vn[I1,I2])./dy)

    if par.ad != 0
        _, maxGradU, _ = get_max_divergence(un, vn, par)
        cd22 = par.ad22 / par.nd^2
        reLambda += 8.0*(par.ad21 + cd22*maxGradU)
    end

    dte = par.cfl / sqrt((reLambda/aStab)^2 + (imLambda/bStab)^2)
    dti = par.cfl / (abs(imLambda/bStab) + 1e-30)

    dt = (ts == "im2") ? dti : dte
    dt = min(dt, par.dtMax)

    if par.idebug > 0 && step == 0
        @printf("getTimeStep: dt = %9.3e (dte=%9.3e, dti=%9.3e, dti/dte=%8.2e, dtMax=%9.2e)\n",
                dt, dte, dti, dti/(dte+1e-100), par.dtMax)
        @printf("  imLambda=%9.3e, reLambda=%9.3e\n", imLambda, reLambda)
    end

    par.dtOld = dt
    return dt
end
