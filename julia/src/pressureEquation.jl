
using SparseArrays, LinearAlgebra

# Solve the pressure Poisson equation.
# factor_matrix=true: (re)build and factor the sparse matrix.
# factor_matrix=false: reuse existing par.dA factorization.
function pressure_equation!(t::Float64, u::Matrix{Float64}, v::Matrix{Float64},
                            dt::Float64, factor_matrix::Bool, par::Par)
    nu  = par.nu
    cdv = par.cdv
    dx  = par.dx;  dy  = par.dy
    iax = par.gid[1,1]; ibx = par.gid[2,1]
    iay = par.gid[1,2]; iby = par.gid[2,2]
    Ngx = par.Ngx;  Ngy = par.Ngy
    numGhost = par.numGhost

    eqn(i1,i2) = 1 + (i1-1) + Ngx*(i2-1)
    Ng  = Ngx*Ngy

    # detect singularity (all-Neumann pressure)
    is_singular = true
    for side in 1:2, axis in 1:2
        bc = par.bc[side,axis]
        if bc == DIRICHLET || bc == PRESSURE_INFLOW || bc == OUTFLOW
            is_singular = false
        end
    end
    Ngs = is_singular ? Ng+1 : Ng

    if factor_matrix
        cpu_factor = time()
        # Interior: 5 entries each; singular adds Ni extra for column + Ni for constraint row
        nzz_est = is_singular ? Ngs*8 : Ngs*6
        ia = zeros(Int, nzz_est);  ja = zeros(Int, nzz_est);  aa = zeros(nzz_est)
        nzz = Ref(0)
        sv! = (ii,jj,val) -> (nzz[]+=1; ia[nzz[]]=ii; ja[nzz[]]=jj; aa[nzz[]]=val)

        I1i, I2i = get_index_interior(par.gid, PC, par)
        for i2 in I2i, i1 in I1i
            ie = eqn(i1,i2)
            sv!(ie, eqn(i1,  i2-1),           1.0/dy^2)
            sv!(ie, eqn(i1-1,i2  ),  1.0/dx^2          )
            sv!(ie, eqn(i1,  i2  ), -2.0*(1.0/dx^2+1.0/dy^2))
            sv!(ie, eqn(i1+1,i2  ),  1.0/dx^2          )
            sv!(ie, eqn(i1,  i2+1),           1.0/dy^2)
            if is_singular; sv!(ie, Ngs, 1.0); end
        end

        if is_singular
            ie = Ngs
            for i2 in I2i, i1 in I1i; sv!(ie, eqn(i1,i2), 1.0); end
        end

        dxv = [dx, dy]
        for side in 1:2, axis in 1:2
            isv = [0,0]; isv[axis] = 1-2*(side-1)
            is1 = isv[1]; is2 = isv[2]
            I1b, I2b = get_boundary_index(side, axis, par)
            bc = par.bc[side,axis]

            if bc == DIRICHLET
                J1b = I1b; J2b = I2b
                if axis == 2
                    i1a_ = par.gid[1,1]; i1b_ = par.gid[2,1]
                    if par.bc[1,1]==DIRICHLET||par.bc[1,1]==PRESSURE_INFLOW; i1a_+=1; end
                    if par.bc[2,1]==DIRICHLET||par.bc[2,1]==PRESSURE_INFLOW; i1b_-=1; end
                    J1b = i1a_:i1b_
                end
                for i2 in J2b, i1 in J1b; ie=eqn(i1,i2); sv!(ie,ie,1.0); end
                for i2 in I2b, i1 in I1b
                    ie=eqn(i1-is1,i2-is2)
                    sv!(ie,ie,1.0); sv!(ie,eqn(i1,i2),-3.0)
                    sv!(ie,eqn(i1+is1,i2+is2),3.0); sv!(ie,eqn(i1+2*is1,i2+2*is2),-1.0)
                end

            elseif bc == PRESSURE_INFLOW
                for i2 in I2b, i1 in I1b
                    ie=eqn(i1,i2); sv!(ie,ie,1.0)
                    ie=eqn(i1-is1,i2-is2)
                    sv!(ie,ie,1.0); sv!(ie,eqn(i1,i2),-3.0)
                    sv!(ie,eqn(i1+is1,i2+is2),3.0); sv!(ie,eqn(i1+2*is1,i2+2*is2),-1.0)
                end

            elseif bc == NO_SLIP_WALL || bc == SLIP_WALL || bc == INFLOW
                h = dxv[axis]
                for i2 in I2b, i1 in I1b
                    ie=eqn(i1-is1,i2-is2)
                    sv!(ie,eqn(i1-is1,i2-is2), 1.0/(2.0*h))
                    sv!(ie,eqn(i1+is1,i2+is2),-1.0/(2.0*h))
                end

            elseif bc == OUTFLOW
                a0=par.outflowPressureCoeffp; a1=par.outflowPressureCoeffpn; h=dxv[axis]
                for i2 in I2b, i1 in I1b
                    ie=eqn(i1-is1,i2-is2)
                    sv!(ie,eqn(i1,i2),a0)
                    sv!(ie,eqn(i1-is1,i2-is2), a1/(2.0*h))
                    sv!(ie,eqn(i1+is1,i2+is2),-a1/(2.0*h))
                end

            elseif bc == PERIODIC
                for i2 in I2b, i1 in I1b
                    ie=eqn(i1-is1,i2-is2)
                    sv!(ie,ie,1.0)
                    sv!(ie,eqn(i1+(ibx-iax)*is1-is1,i2+(iby-iay)*is2-is2),-1.0)
                end
            end
        end

        # corner diagonal extrapolation
        for side1 in 0:1, side2 in 0:1
            is1=1-2*side1; is2=1-2*side2
            i1=par.gid[side1+1,1]-is1; i2=par.gid[side2+1,2]-is2
            ie=eqn(i1,i2)
            sv!(ie,ie,1.0); sv!(ie,eqn(i1+is1,i2+is2),-3.0)
            sv!(ie,eqn(i1+2*is1,i2+2*is2),3.0); sv!(ie,eqn(i1+3*is1,i2+3*is2),-1.0)
        end

        if par.idebug > 0
            @printf("Optimized fill of pressure matrix: nzzEst=%d, nzz=%d\n", nzz_est, nzz[])
        end
        A = sparse(ia[1:nzz[]], ja[1:nzz[]], aa[1:nzz[]], Ngs, Ngs)
        par.dA = lu(A)
    end

    cpu0 = time()  # measure solve time separately

    # --- build RHS (full-size arrays as in MATLAB) ---
    rhs = zeros(Ngs)
    I1, I2 = get_index(par.gid)

    pf = zeros(Ngx, Ngy)
    pf[I1,I2] .= par.pfe(par.x[I1,I2,1], par.x[I1,I2,2], t)

    ux = zeros(Ngx,Ngy); uy = zeros(Ngx,Ngy)
    vx = zeros(Ngx,Ngy); vy = zeros(Ngx,Ngy)
    ux[I1,I2] .= (u[I1.+1,I2] .- u[I1.-1,I2]) ./ (2.0*dx)
    uy[I1,I2] .= (u[I1,I2.+1] .- u[I1,I2.-1]) ./ (2.0*dy)
    vx[I1,I2] .= (v[I1.+1,I2] .- v[I1.-1,I2]) ./ (2.0*dx)
    vy[I1,I2] .= (v[I1,I2.+1] .- v[I1,I2.-1]) ./ (2.0*dy)

    I1i, I2i = get_index_interior(par.gid, PC, par)
    for i2 in I2i, i1 in I1i
        ie = eqn(i1,i2)
        div_damp = (cdv/dt)*(ux[i1,i2]+vy[i1,i2])
        rhs[ie] = -(ux[i1,i2]^2 + 2.0*uy[i1,i2]*vx[i1,i2] + vy[i1,i2]^2) + div_damp + pf[i1,i2]
    end

    add_forcing = (par.ms != "none" || par.knownSolution != "none")
    ms_on       = (par.ms != "none")

    for side in 1:2, axis in 1:2
        isv = [0,0]; isv[axis] = 1-2*(side-1)
        is1 = isv[1]; is2 = isv[2]
        I1b, I2b = get_boundary_index(side, axis, par)
        bc = par.bc[side,axis]

        for i2 in I2b, i1 in I1b
            x1 = par.x[i1,i2,1]; x2 = par.x[i1,i2,2]

            if bc == DIRICHLET
                rhs[eqn(i1,i2)] = par.pe(x1,x2,t)

            elseif bc == PRESSURE_INFLOW
                rhs[eqn(i1,i2)] = add_forcing ? par.pe(x1,x2,t) : par.pressureInflowValue

            elseif bc == NO_SLIP_WALL || bc == SLIP_WALL || bc == INFLOW
                ie = eqn(i1-is1, i2-is2)
                # p.n = nu*(u.xx+u.yy) on wall; curl-curl: axis==1 → -nu*v.xy, axis==2 → -nu*u.xy
                vxy = (v[i1+1,i2+1] - v[i1-1,i2+1] - v[i1+1,i2-1] + v[i1-1,i2-1]) / (4.0*dx*dy)
                uxy = (u[i1+1,i2+1] - u[i1-1,i2+1] - u[i1+1,i2-1] + u[i1-1,i2-1]) / (4.0*dx*dy)
                if axis == 1
                    rhs[ie] = -is1*nu*(-vxy)
                    if add_forcing
                        rhs[ie] += -is1*(par.pex(x1,x2,t) + nu*par.vexy(x1,x2,t))
                    end
                else
                    rhs[ie] = -is2*nu*(-uxy)
                    if add_forcing
                        rhs[ie] += -is2*(par.pey(x1,x2,t) + nu*par.uexy(x1,x2,t))
                    end
                end

            elseif bc == OUTFLOW
                a0=par.outflowPressureCoeffp; a1=par.outflowPressureCoeffpn
                ie = eqn(i1-is1,i2-is2)
                if ms_on
                    rhs[ie] = a0*par.pe(x1,x2,t) - isv[axis]*a1*par.pex(x1,x2,t)
                else
                    rhs[ie] = a0*par.pOutflow
                end

            elseif bc == PERIODIC
                rhs[eqn(i1-is1,i2-is2)] = 0.0
            end
        end
    end

    if is_singular && add_forcing
        # Constraint row sums only over interior pressure unknowns
        extra_val = 0.0
        for i2 in I2i, i1 in I1i
            extra_val += par.pe(par.x[i1,i2,1], par.x[i1,i2,2], t)
        end
        rhs[Ngs] = extra_val
    end

    sol = par.dA \ rhs

    p = zeros(Ngx, Ngy)
    I1g, I2g = get_index(par.gid, numGhost)
    for i2 in I2g, i1 in I1g
        p[i1,i2] = sol[eqn(i1,i2)]
    end

    par.cpuPressure += time() - cpu0
    return p
end
