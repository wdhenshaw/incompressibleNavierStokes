
# Compute momentum RHS:  ut = -(u*ux + v*uy + px) + nu*scale*(uxx+uyy) + forcing
function get_ut(t::Float64, un::Matrix{Float64}, vn::Matrix{Float64}, pn::Matrix{Float64},
                nu_scale_factor::Float64, par::Par)

    cpu0 = time()
    dx = par.dx;  dy = par.dy
    nu = par.nu

    I1, I2 = get_index(par.gid)

    ut = zeros(par.Ngx, par.Ngy)
    vt = zeros(par.Ngx, par.Ngy)

    # precompute forcing
    uf = zeros(par.Ngx, par.Ngy)
    vf = zeros(par.Ngx, par.Ngy)
    if par.ms != "none"
        uf[I1,I2] .= par.ufe(par.x[I1,I2,1], par.x[I1,I2,2], t)
        vf[I1,I2] .= par.vfe(par.x[I1,I2,1], par.x[I1,I2,2], t)
    end

    # central difference operators on index ranges
    ux = (un[I1.+1, I2] .- un[I1.-1, I2]) ./ (2.0*dx)
    uy = (un[I1, I2.+1] .- un[I1, I2.-1]) ./ (2.0*dy)
    vx = (vn[I1.+1, I2] .- vn[I1.-1, I2]) ./ (2.0*dx)
    vy = (vn[I1, I2.+1] .- vn[I1, I2.-1]) ./ (2.0*dy)
    px = (pn[I1.+1, I2] .- pn[I1.-1, I2]) ./ (2.0*dx)
    py = (pn[I1, I2.+1] .- pn[I1, I2.-1]) ./ (2.0*dy)

    adv_u = un[I1,I2].*ux .+ vn[I1,I2].*uy .+ px
    adv_v = un[I1,I2].*vx .+ vn[I1,I2].*vy .+ py

    if nu_scale_factor != 0.0
        uxx = (un[I1.+1,I2] .- 2.0.*un[I1,I2] .+ un[I1.-1,I2]) ./ dx^2
        uyy = (un[I1,I2.+1] .- 2.0.*un[I1,I2] .+ un[I1,I2.-1]) ./ dy^2
        vxx = (vn[I1.+1,I2] .- 2.0.*vn[I1,I2] .+ vn[I1.-1,I2]) ./ dx^2
        vyy = (vn[I1,I2.+1] .- 2.0.*vn[I1,I2] .+ vn[I1,I2.-1]) ./ dy^2
        ut[I1,I2] .= -adv_u .+ nu.*(uxx .+ uyy) .+ uf[I1,I2]
        vt[I1,I2] .= -adv_v .+ nu.*(vxx .+ vyy) .+ vf[I1,I2]
    else
        ut[I1,I2] .= -adv_u .+ uf[I1,I2]
        vt[I1,I2] .= -adv_v .+ vf[I1,I2]
    end

    if par.ad != 0

        ad21 = par.ad21
        cd22 = par.ad22 / par.nd^2
        adc  = ad21 .+ cd22.*(abs.(ux) .+ abs.(uy) .+ abs.(vx) .+ abs.(vy))

        # @printf("Add artificial dissipation:t=%9.2e,  ad21=%9.2e ad22=%9.2e, cd22=%9.2e..\n",t,par.ad21,par.ad22,cd22);
        ut[I1,I2] .+= adc.*(un[I1.+1,I2] .- 4.0.*un[I1,I2] .+ un[I1.-1,I2] .+ un[I1,I2.+1] .+ un[I1,I2.-1])
        vt[I1,I2] .+= adc.*(vn[I1.+1,I2] .- 4.0.*vn[I1,I2] .+ vn[I1.-1,I2] .+ vn[I1,I2.+1] .+ vn[I1,I2.-1])
    end

    par.cpuGetUt += time() - cpu0
    return ut, vt
end
