
# Return boundary indices adjusted to avoid double-counting corners
# at Dirichlet-type junctions.  Used for velocity implicit matrix and RHS assembly.
function get_adjusted_boundary_index(side::Int, axis::Int, par::Par)
    axis2 = mod(axis, par.nd) + 1   # adjacent axis: 1→2, 2→1
    bid = [par.gid[1,axis2], par.gid[2,axis2]]

    bc_here = par.bc[side, axis]
    if bc_here == DIRICHLET || bc_here == NO_SLIP_WALL || bc_here == INFLOW
        dirichlet_type = 1
    elseif bc_here == SLIP_WALL || bc_here == PRESSURE_INFLOW
        dirichlet_type = 2
    elseif bc_here == OUTFLOW
        dirichlet_type = 0
    else
        error("get_adjusted_boundary_index: unknown bc=$(par.bc[side,axis]) for (side,axis)=($side,$axis)")
    end

    if (dirichlet_type != 0 && axis == 2) || dirichlet_type != 1
        for side2 in 1:2
            is2 = 1 - 2*(side2-1)
            bc2 = par.bc[side2, axis2]
            if (bc2 == DIRICHLET || bc2 == NO_SLIP_WALL || bc2 == INFLOW) ||
               (bc_here == SLIP_WALL && bc2 == SLIP_WALL)
                bid[side2] = par.gid[side2, axis2] + is2
            end
        end
    end

    if axis == 1
        I1b = par.gid[side,axis]:par.gid[side,axis]
        I2b = bid[1]:bid[2]
    else
        I1b = bid[1]:bid[2]
        I2b = par.gid[side,axis]:par.gid[side,axis]
    end
    return I1b, I2b
end
