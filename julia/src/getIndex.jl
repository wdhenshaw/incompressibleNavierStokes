
# Return index ranges for all boundary-to-boundary points (optionally extended by `extra` ghost layers)
function get_index(gid::Matrix{Int}, extra::Int=0)
    I1 = (gid[1,1]-extra):(gid[2,1]+extra)
    I2 = (gid[1,2]-extra):(gid[2,2]+extra)
    return I1, I2
end
