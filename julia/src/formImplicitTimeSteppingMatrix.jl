
using SparseArrays, LinearAlgebra

# Build and factor M = I - (0.5*dt*nu)*L for the implicit velocity solve.
# Stores LU factorizations in par.dAimp[1] (and [2] if multipleImplicitSolversNeeded).
function form_implicit_time_stepping_matrix!(dt::Float64, par::Par)
    cpu0 = time()
    nu  = par.nu
    Ngx = par.Ngx; Ngy = par.Ngy
    iax = par.gid[1,1]; ibx = par.gid[2,1]
    iay = par.gid[1,2]; iby = par.gid[2,2]
    dx  = par.dx;  dy  = par.dy

    eqn(i1,i2) = 1 + (i1-1) + Ngx*(i2-1)
    Ng = Ngx*Ngy

    numImplicitSolvers = par.multipleImplicitSolversNeeded ? 2 : 1
    par.dAimp = Vector{Any}(undef, numImplicitSolvers)

    dxv = [dx, dy]

    for iuv in 1:numImplicitSolvers
        nzz_est = Ng*8
        ia = zeros(Int, nzz_est); ja = zeros(Int, nzz_est); aa = zeros(nzz_est)
        nzz = Ref(0)
        sv! = (ii,jj,val) -> (nzz[]+=1; ia[nzz[]]=ii; ja[nzz[]]=jj; aa[nzz[]]=val)

        component = (iuv == 1) ? UC : VC
        I1a, I2a = get_index_interior(par.gid, component, par)

        # Interior: M*u = (I - 0.5*dt*nu*L)*u
        for i2 in I2a, i1 in I1a
            ie = eqn(i1,i2)
            sv!(ie, eqn(i1,  i2-1),     -(0.5*nu*dt)*(          1.0/dy^2))
            sv!(ie, eqn(i1-1,i2  ),     -(0.5*nu*dt)*( 1.0/dx^2          ))
            sv!(ie, eqn(i1,  i2  ), 1.0-(0.5*nu*dt)*(-2.0*(1.0/dx^2+1.0/dy^2)))
            sv!(ie, eqn(i1+1,i2  ),     -(0.5*nu*dt)*( 1.0/dx^2          ))
            sv!(ie, eqn(i1,  i2+1),     -(0.5*nu*dt)*(          1.0/dy^2))
        end

        for side in 1:2, axis in 1:2
            isv = [0,0]; isv[axis] = 1-2*(side-1)
            is1 = isv[1]; is2 = isv[2]
            h = dxv[axis]
            I1b, I2b = get_boundary_index(side, axis, par)
            bc = par.bc[side,axis]

            if bc == DIRICHLET || bc == INFLOW || bc == NO_SLIP_WALL
                I1adj, I2adj = get_adjusted_boundary_index(side, axis, par)
                for i2 in I2adj, i1 in I1adj; ie=eqn(i1,i2); sv!(ie,ie,1.0); end
                for i2 in I2b, i1 in I1b
                    ie=eqn(i1-is1,i2-is2)
                    sv!(ie,ie,1.0); sv!(ie,eqn(i1,i2),-3.0)
                    sv!(ie,eqn(i1+is1,i2+is2),3.0); sv!(ie,eqn(i1+2*is1,i2+2*is2),-1.0)
                end

            elseif bc == SLIP_WALL
                I1adj, I2adj = get_adjusted_boundary_index(side, axis, par)
                for i2 in I2adj, i1 in I1adj
                    # Dirichlet only on the normal component for this solver
                    if (axis==1 && iuv==1) || (axis==2 && iuv==2)
                        ie=eqn(i1,i2); sv!(ie,ie,1.0)
                    end
                end
                for i2 in I2b, i1 in I1b
                    if (axis==1 && iuv==1) || (axis==2 && iuv==2)
                        ie=eqn(i1-is1,i2-is2)
                        sv!(ie,ie,1.0); sv!(ie,eqn(i1,i2),-3.0)
                        sv!(ie,eqn(i1+is1,i2+is2),3.0); sv!(ie,eqn(i1+2*is1,i2+2*is2),-1.0)
                    else
                        ie=eqn(i1-is1,i2-is2)
                        sv!(ie,eqn(i1-is1,i2-is2), 1.0/(2.0*h))
                        sv!(ie,eqn(i1+is1,i2+is2),-1.0/(2.0*h))
                    end
                end

            elseif bc == PRESSURE_INFLOW
                I1adj, I2adj = get_adjusted_boundary_index(side, axis, par)
                for i2 in I2adj, i1 in I1adj
                    # Dirichlet on the tangential component
                    if (axis==1 && iuv==2) || (axis==2 && iuv==1)
                        ie=eqn(i1,i2); sv!(ie,ie,1.0)
                    end
                end
                for i2 in I2b, i1 in I1b
                    # Extrapolate ghost (both components)
                    ie=eqn(i1-is1,i2-is2)
                    sv!(ie,ie,1.0); sv!(ie,eqn(i1,i2),-3.0)
                    sv!(ie,eqn(i1+is1,i2+is2),3.0); sv!(ie,eqn(i1+2*is1,i2+2*is2),-1.0)
                end

            elseif bc == OUTFLOW
                for i2 in I2b, i1 in I1b
                    ie=eqn(i1-is1,i2-is2)
                    sv!(ie,ie,1.0); sv!(ie,eqn(i1,i2),-3.0)
                    sv!(ie,eqn(i1+is1,i2+is2),3.0); sv!(ie,eqn(i1+2*is1,i2+2*is2),-1.0)
                end

            elseif bc == PERIODIC
                for i2 in I2b, i1 in I1b
                    ie=eqn(i1-is1,i2-is2)
                    sv!(ie,ie,1.0)
                    sv!(ie,eqn(i1+(ibx-iax)*is1-is1,i2+(iby-iay)*is2-is2),-1.0)
                end
            end
        end

        # Corner ghost extrapolation
        if !(par.bc[1,1]==PERIODIC && par.bc[1,2]==PERIODIC)
            for side1 in 0:1, side2 in 0:1
                is1=1-2*side1; is2=1-2*side2
                i1=par.gid[side1+1,1]-is1; i2=par.gid[side2+1,2]-is2
                ie=eqn(i1,i2)
                sv!(ie,ie,1.0); sv!(ie,eqn(i1+is1,i2+is2),-3.0)
                sv!(ie,eqn(i1+2*is1,i2+2*is2),3.0); sv!(ie,eqn(i1+3*is1,i2+3*is2),-1.0)
            end
        else
            for side1 in 0:1, side2 in 0:1
                is1=1-2*side1; is2=1-2*side2
                i1=par.gid[side1+1,1]-is1; i2=par.gid[side2+1,2]-is2
                ie=eqn(i1,i2)
                sv!(ie,ie,1.0)
                sv!(ie,eqn(i1+(ibx-iax)*is1,i2+(iby-iay)*is2),-1.0)
            end
        end

        # Slip-slip corner: both components pinned
        for side1 in 1:2, side2 in 1:2
            if par.bc[side1,1]==SLIP_WALL && par.bc[side2,2]==SLIP_WALL
                i1=par.gid[side1,1]; i2=par.gid[side2,2]
                ie=eqn(i1,i2); sv!(ie,ie,1.0)
            end
        end

        if par.idebug > 0
            @printf("form_implicit_time_stepping_matrix[%d]: nzz_est=%d, nzz=%d\n", iuv, nzz_est, nzz[])
        end

        AA = sparse(ia[1:nzz[]], ja[1:nzz[]], aa[1:nzz[]], Ng, Ng)
        par.dAimp[iuv] = lu(AA)
        par.cpuFactorImpMatrix += time() - cpu0
        cpu0 = time()
    end
end
