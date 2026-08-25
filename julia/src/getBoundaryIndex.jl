
# Return index ranges for all points on boundary (side,axis)
function get_boundary_index(side::Int, axis::Int, par::Par)
    if axis == 1
        I1b = par.gid[side,axis]:par.gid[side,axis]
        I2b = par.gid[1,2]:par.gid[2,2]
    elseif axis == 2
        I1b = par.gid[1,1]:par.gid[2,1]
        I2b = par.gid[side,axis]:par.gid[side,axis]
    else
        error("get_boundary_index: invalid axis=$axis")
    end
    return I1b, I2b
end
