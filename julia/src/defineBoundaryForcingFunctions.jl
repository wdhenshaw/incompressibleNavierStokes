
function define_boundary_forcing_functions!(par::Par)
    if par.ic == "shear"
        par.guax = (x,y,t) ->  0.0
        par.gubx = (x,y,t) ->  0.0
        par.guay = (x,y,t) -> -1.0
        par.guby = (x,y,t) -> +1.0
        par.gvax = (x,y,t) ->  0.0
        par.gvbx = (x,y,t) ->  0.0
        par.gvay = (x,y,t) ->  0.0
        par.gvby = (x,y,t) ->  0.0
        par.pgax = (x,y,t) ->  0.0
        par.pgbx = (x,y,t) ->  0.0
        par.pgay = (x,y,t) ->  0.0
        par.pgby = (x,y,t) ->  0.0

    elseif par.ms != "none" || par.knownSolution != "none"
        par.guax = (x,y,t) -> par.ue(x,y,t)
        par.gubx = (x,y,t) -> par.ue(x,y,t)
        par.guay = (x,y,t) -> par.ue(x,y,t)
        par.guby = (x,y,t) -> par.ue(x,y,t)
        par.u0   = (x,y)   -> par.ue(x,y,0.0)

        par.gvax = (x,y,t) -> par.ve(x,y,t)
        par.gvbx = (x,y,t) -> par.ve(x,y,t)
        par.gvay = (x,y,t) -> par.ve(x,y,t)
        par.gvby = (x,y,t) -> par.ve(x,y,t)
        par.v0   = (x,y)   -> par.ve(x,y,0.0)

        # pressure BC functions depend on BC type
        bc11 = par.bc[1,1];  bc21 = par.bc[2,1]
        bc12 = par.bc[1,2];  bc22 = par.bc[2,2]

        if bc11 == DIRICHLET || bc11 == PERIODIC || bc11 == PRESSURE_INFLOW
            par.pgax = (x,y,t) -> par.pe(x,y,t)
        elseif bc11 == NO_SLIP_WALL || bc11 == INFLOW || bc11 == SLIP_WALL
            par.pgax = (x,y,t) -> -par.pex(x,y,t)
        else
            error("define_boundary_forcing_functions: pgax: finish me for bc=$bc11")
        end

        if bc21 == DIRICHLET || bc21 == PERIODIC || bc21 == PRESSURE_INFLOW
            par.pgbx = (x,y,t) -> par.pe(x,y,t)
        elseif bc21 == NO_SLIP_WALL || bc21 == INFLOW || bc21 == SLIP_WALL
            par.pgbx = (x,y,t) -> par.pex(x,y,t)
        elseif bc21 == OUTFLOW
            par.pgbx = (x,y,t) -> par.outflowPressureCoeffp.*par.pe(x,y,t) .+ par.outflowPressureCoeffpn.*par.pex(x,y,t)
        else
            error("define_boundary_forcing_functions: pgbx: finish me for bc=$bc21")
        end

        if bc12 == DIRICHLET || bc12 == PERIODIC || bc12 == PRESSURE_INFLOW
            par.pgay = (x,y,t) -> par.pe(x,y,t)
        elseif bc12 == NO_SLIP_WALL || bc12 == INFLOW || bc12 == SLIP_WALL
            par.pgay = (x,y,t) -> -par.pey(x,y,t)
        else
            error("define_boundary_forcing_functions: pgay: finish me for bc=$bc12")
        end

        if bc22 == DIRICHLET || bc22 == PERIODIC || bc22 == PRESSURE_INFLOW
            par.pgby = (x,y,t) -> par.pe(x,y,t)
        elseif bc22 == NO_SLIP_WALL || bc22 == INFLOW || bc22 == SLIP_WALL
            par.pgby = (x,y,t) -> par.pey(x,y,t)
        else
            error("define_boundary_forcing_functions: pgby: finish me for bc=$bc22")
        end

    else
        # default case
        par.guax = (x,y,t) -> par.uInflow
        par.gubx = (x,y,t) -> 0.0
        par.guay = (x,y,t) -> 0.0
        par.guby = (x,y,t) -> 0.0
        par.gvax = (x,y,t) -> 0.0
        par.gvbx = (x,y,t) -> 0.0
        par.gvay = (x,y,t) -> 0.0
        par.gvby = (x,y,t) -> 0.0
        par.pgax = (x,y,t) -> 0.0
        par.pgbx = (x,y,t) -> 0.0
        par.pgay = (x,y,t) -> 0.0
        par.pgby = (x,y,t) -> 0.0
    end

    par.gu = Any[par.guax, par.gubx, par.guay, par.guby]
    par.gv = Any[par.gvax, par.gvbx, par.gvay, par.gvby]
    par.gp = Any[par.pgax, par.pgbx, par.pgay, par.pgby]
end
