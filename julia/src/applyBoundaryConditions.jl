
# 3rd-order extrapolation: u at ghost from three interior points inward
# is1,is2 = shift direction toward interior (+1 or -1 per axis)
function extrap3(u::Matrix{Float64}, i1::Int, i2::Int, is1::Int, is2::Int)
    return 3.0*u[i1+is1,i2+is2] - 3.0*u[i1+2*is1,i2+2*is2] + u[i1+3*is1,i2+3*is2]
end

# Vectorized extrap3 over index ranges
function extrap3!(u::Matrix{Float64}, I1g, I2g, is1::Int, is2::Int)
    for i2 in I2g, i1 in I1g
        u[i1,i2] = extrap3(u, i1, i2, is1, is2)
    end
end

function apply_boundary_conditions!(u::Matrix{Float64}, v::Matrix{Float64},
                                    t::Float64, par::Par)
    cpu0 = time()
    dx = par.dx;  dy = par.dy

    ms_on = par.ms != "none"

    # ---- STAGE 1: Dirichlet-type BCs on boundary points ----
    for side in 1:2, axis in 1:2
        isv = [0, 0]; isv[axis] = 1 - 2*(side-1)
        is1 = isv[1]; is2 = isv[2]
        mbc = side + 2*(axis-1)
        I1b, I2b = get_boundary_index(side, axis, par)
        bc = par.bc[side, axis]

        if bc == DIRICHLET || bc == NO_SLIP_WALL || bc == INFLOW
            u[I1b, I2b] .= par.gu[mbc](par.x[I1b, I2b, 1], par.x[I1b, I2b, 2], t)
            v[I1b, I2b] .= par.gv[mbc](par.x[I1b, I2b, 1], par.x[I1b, I2b, 2], t)

        elseif bc == SLIP_WALL
            n1 = Float64(-is1);  n2 = Float64(-is2)
            nDotU = n1.*u[I1b,I2b] .+ n2.*v[I1b,I2b]
            if ms_on
                nDotU .-= n1.*par.ue(par.x[I1b,I2b,1], par.x[I1b,I2b,2], t)
                nDotU .-= n2.*par.ve(par.x[I1b,I2b,1], par.x[I1b,I2b,2], t)
            end
            u[I1b,I2b] .-= nDotU.*n1
            v[I1b,I2b] .-= nDotU.*n2

        elseif bc == PRESSURE_INFLOW
            n1 = Float64(-is1);  n2 = Float64(-is2)
            nDotU = n1.*u[I1b,I2b] .+ n2.*v[I1b,I2b]
            u[I1b,I2b] .= nDotU.*n1
            v[I1b,I2b] .= nDotU.*n2
            if ms_on
                if axis == 1
                    v[I1b,I2b] .= par.ve(par.x[I1b,I2b,1], par.x[I1b,I2b,2], t)
                else
                    u[I1b,I2b] .= par.ue(par.x[I1b,I2b,1], par.x[I1b,I2b,2], t)
                end
            end

        elseif bc == PERIODIC || bc == OUTFLOW
            # handled below

        else
            error("apply_boundary_conditions: unknown BC=$(par.bc[side,axis])")
        end
    end

    # ---- STAGE 2: ghost point extrapolation (first pass: extrapolate all) ----
    for side in 1:2, axis in 1:2
        isv = [0,0]; isv[axis] = 1-2*(side-1)
        is1 = isv[1]; is2 = isv[2]
        I1b, I2b = get_boundary_index(side, axis, par)
        I1g = (I1b[1]-is1):(I1b[end]-is1)
        I2g = (I2b[1]-is2):(I2b[end]-is2)
        extrap3!(u, I1g, I2g, is1, is2)
        extrap3!(v, I1g, I2g, is1, is2)
    end

    # ---- STAGE 2b: override ghost for special BCs ----
    for side in 1:2, axis in 1:2
        isv = [0,0]; isv[axis] = 1-2*(side-1)
        is1 = isv[1]; is2 = isv[2]
        I1b, I2b = get_boundary_index(side, axis, par)
        I1g = (I1b[1]-is1):(I1b[end]-is1)
        I2g = (I2b[1]-is2):(I2b[end]-is2)
        bc = par.bc[side, axis]

        if bc == SLIP_WALL
            if axis == 1
                gvx = ms_on ? Float64(-is1)*par.vex(par.x[I1b,I2b,1], par.x[I1b,I2b,2], t) : zeros(length(I2b))
                v[I1g,I2g] .= v[I1b.+is1, I2b.+is2] .+ (2.0*dx).*gvx
            else
                guy = ms_on ? Float64(-is2)*par.uey(par.x[I1b,I2b,1], par.x[I1b,I2b,2], t) : zeros(length(I1b))
                u[I1g,I2g] .= u[I1b.+is1, I2b.+is2] .+ (2.0*dy).*guy
            end
        end
        # outflow and pressureInflow ghosts were already extrapolated in first pass
    end

    # ---- STAGE 3: divergence constraint on ghost points (Cartesian) ----
    for side in 1:2, axis in 1:2
        isv = [0,0]; isv[axis] = 1-2*(side-1)
        is1 = isv[1]; is2 = isv[2]
        I1b, I2b = get_boundary_index(side, axis, par)
        I1g = (I1b[1]-is1):(I1b[end]-is1)
        I2g = (I2b[1]-is2):(I2b[end]-is2)
        bc = par.bc[side, axis]

        if bc == NO_SLIP_WALL || bc == SLIP_WALL || bc == INFLOW ||
           bc == PRESSURE_INFLOW || bc == OUTFLOW
            if axis == 1
                # u.x = -v.y  =>  u_ghost from: u(ghost) = u(int) + 2*dx*is1*(-v.y at boundary)
                vy_bd = (v[I1b, I2b.+1] .- v[I1b, I2b.-1]) ./ (2.0*dy)
                u[I1g, I2g] .= u[I1b.+is1, I2b.+is2] .+ (2.0*dx*is1).*vy_bd
            else
                # v.y = -u.x  =>  v_ghost from: v(ghost) = v(int) + 2*dy*is2*(-u.x at boundary)
                ux_bd = (u[I1b.+1, I2b] .- u[I1b.-1, I2b]) ./ (2.0*dx)
                v[I1g, I2g] .= v[I1b.+is1, I2b.+is2] .+ (2.0*dy*is2).*ux_bd
            end
        end
    end

    # ---- STAGE 4: corner diagonal extrapolation ----
    for side1 in 1:2, side2 in 1:2
        is1 = 1 - 2*(side1-1)
        is2 = 1 - 2*(side2-1)
        ix = par.gid[side1,1] - is1
        iy = par.gid[side2,2] - is2
        u[ix,iy] = extrap3(u, ix, iy, is1, is2)
        v[ix,iy] = extrap3(v, ix, iy, is1, is2)
    end

    # ---- STAGE 5: periodic BCs ----
    for axis in 1:2
        if par.bc[1,axis] == PERIODIC || par.bc[2,axis] == PERIODIC
            if !(par.bc[1,axis] == PERIODIC && par.bc[2,axis] == PERIODIC)
                error("apply_boundary_conditions: periodic on one side but not the other for axis=$axis")
            end
            isv = [0,0]; isv[axis] = 1
            is1 = isv[1]; is2 = isv[2]

            I1a, I2a = get_boundary_index(1, axis, par)
            I1b_p, I2b_p = get_boundary_index(2, axis, par)
            # extend to full array in tangential direction
            if axis == 1;  I2a = 1:par.Ngy;  I2b_p = 1:par.Ngy
            else;          I1a = 1:par.Ngx;  I1b_p = 1:par.Ngx  end

            # right/top = left/bottom
            u[I1b_p, I2b_p] .= u[I1a, I2a]
            v[I1b_p, I2b_p] .= v[I1a, I2a]
            # right/top ghost = left/bottom ghost
            u[I1b_p.+is1, I2b_p.+is2] .= u[I1a.+is1, I2a.+is2]
            v[I1b_p.+is1, I2b_p.+is2] .= v[I1a.+is1, I2a.+is2]
            # left/bottom ghost = right/top ghost
            u[I1a.-is1, I2a.-is2] .= u[I1b_p.-is1, I2b_p.-is2]
            v[I1a.-is1, I2a.-is2] .= v[I1b_p.-is1, I2b_p.-is2]
        end
    end

    par.cpuBC += time() - cpu0
end
