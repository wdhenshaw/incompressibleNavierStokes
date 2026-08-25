
# Solve the implicit velocity system M*u^{n+1} = rhs
# unp1, vnp1 on input hold the RHS (explicit part already in-place).
# On output, unp1, vnp1 are overwritten with the implicit solution.
function solve_implicit_time_step!(unp1::Matrix{Float64}, vnp1::Matrix{Float64},
                                   tnp1::Float64, par::Par)
    cpu0 = time()
    Ngx = par.Ngx; Ngy = par.Ngy
    iax = par.gid[1,1]; ibx = par.gid[2,1]
    iay = par.gid[1,2]; iby = par.gid[2,2]
    numGhost = par.numGhost

    eqn(i1,i2) = 1 + (i1-1) + Ngx*(i2-1)
    Ng = Ngx*Ngy

    rhsu = zeros(Ng); rhsv = zeros(Ng)
    for i2 in iay:iby, i1 in iax:ibx
        ie = eqn(i1,i2)
        rhsu[ie] = unp1[i1,i2]
        rhsv[ie] = vnp1[i1,i2]
    end

    ms_on = (par.ms != "none" || par.knownSolution != "none")

    for side in 1:2, axis in 1:2
        isv = [0,0]; isv[axis] = 1-2*(side-1)
        is1 = isv[1]; is2 = isv[2]
        mbc = side + 2*(axis-1)
        I1b, I2b = get_boundary_index(side, axis, par)
        bc = par.bc[side,axis]

        if bc == DIRICHLET || bc == NO_SLIP_WALL || bc == INFLOW
            I1adj, I2adj = get_adjusted_boundary_index(side, axis, par)
            for i2 in I2adj, i1 in I1adj
                ie = eqn(i1,i2)
                if ms_on
                    rhsu[ie] = par.ue(par.x[i1,i2,1], par.x[i1,i2,2], tnp1)
                    rhsv[ie] = par.ve(par.x[i1,i2,1], par.x[i1,i2,2], tnp1)
                else
                    rhsu[ie] = 0.0; rhsv[ie] = 0.0
                end
            end

        elseif bc == SLIP_WALL
            n1 = Float64(-is1); n2 = Float64(-is2)
            I1adj, I2adj = get_adjusted_boundary_index(side, axis, par)
            for i2 in I2adj, i1 in I1adj
                ie = eqn(i1,i2)
                if axis == 1   # left/right: set u
                    rhsu[ie] = ms_on ? par.gu[mbc](par.x[i1,i2,1], par.x[i1,i2,2], tnp1) : 0.0
                else           # bottom/top: set v
                    rhsv[ie] = ms_on ? par.gv[mbc](par.x[i1,i2,1], par.x[i1,i2,2], tnp1) : 0.0
                end
            end
            for i2 in I2b, i1 in I1b
                ie = eqn(i1-is1, i2-is2)
                if axis == 1   # v.x = ...
                    rhsv[ie] = ms_on ? n1*par.vex(par.x[i1,i2,1], par.x[i1,i2,2], tnp1) : 0.0
                else           # u.y = ...
                    rhsu[ie] = ms_on ? n2*par.uey(par.x[i1,i2,1], par.x[i1,i2,2], tnp1) : 0.0
                end
            end

        elseif bc == PRESSURE_INFLOW
            I1adj, I2adj = get_adjusted_boundary_index(side, axis, par)
            for i2 in I2adj, i1 in I1adj
                ie = eqn(i1,i2)
                if axis == 1   # left/right: set v (tangential)
                    rhsv[ie] = ms_on ? par.gv[mbc](par.x[i1,i2,1], par.x[i1,i2,2], tnp1) : 0.0
                else           # bottom/top: set u (tangential)
                    rhsu[ie] = ms_on ? par.gu[mbc](par.x[i1,i2,1], par.x[i1,i2,2], tnp1) : 0.0
                end
            end
        # OUTFLOW, PERIODIC: nothing to do
        end
    end

    # Slip-slip corner: both components pinned (mbc=4 at end of loop above, matching MATLAB)
    for side1 in 1:2, side2 in 1:2
        if par.bc[side1,1]==SLIP_WALL && par.bc[side2,2]==SLIP_WALL
            i1=par.gid[side1,1]; i2=par.gid[side2,2]
            ie = eqn(i1,i2)
            rhsu[ie] = ms_on ? par.gu[4](par.x[i1,i2,1], par.x[i1,i2,2], tnp1) : 0.0
            rhsv[ie] = ms_on ? par.gv[4](par.x[i1,i2,1], par.x[i1,i2,2], tnp1) : 0.0
        end
    end

    rhsu = par.dAimp[1] \ rhsu
    rhsv = par.multipleImplicitSolversNeeded ? par.dAimp[2] \ rhsv : par.dAimp[1] \ rhsv

    I1g = (iax-numGhost):(ibx+numGhost)
    I2g = (iay-numGhost):(iby+numGhost)
    for i2 in I2g, i1 in I1g
        ie = eqn(i1,i2)
        unp1[i1,i2] = rhsu[ie]
        vnp1[i1,i2] = rhsv[ie]
    end

    par.cpuImplicit += time() - cpu0
end
