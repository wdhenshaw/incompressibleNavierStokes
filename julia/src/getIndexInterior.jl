
# Return index ranges for interior equation points for a given component (UC, VC, or PC).
# Excludes boundary points where Dirichlet-type conditions are imposed.
function get_index_interior(gid::Matrix{Int}, component::Int, par::Par)
    gidi = copy(gid)

    for axis in 1:par.nd
        for side in 1:2
            is = 1 - 2*(side-1)
            bc = par.bc[side, axis]

            if bc == DIRICHLET
                gidi[side, axis] += is

            elseif bc == PERIODIC || bc == OUTFLOW
                # include boundary points

            elseif bc == NO_SLIP_WALL || bc == INFLOW
                if component != PC
                    gidi[side, axis] += is
                end

            elseif bc == SLIP_WALL
                if (axis == 1 && component == UC) || (axis == 2 && component == VC)
                    gidi[side, axis] += is
                end

            elseif bc == PRESSURE_INFLOW
                if component == PC ||
                   (axis == 1 && component == VC) ||
                   (axis == 2 && component == UC)
                    gidi[side, axis] += is
                end

            else
                error("get_index_interior: unknown bc=$bc")
            end
        end
    end

    return get_index(gidi)
end
