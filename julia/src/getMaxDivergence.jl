
function get_max_divergence(un::Matrix{Float64}, vn::Matrix{Float64}, par::Par)
    dx = par.dx;  dy = par.dy
    I1, I2 = get_index(par.gid)

    ux = (un[I1.+1, I2] .- un[I1.-1, I2]) ./ (2.0*dx)
    vy = (vn[I1, I2.+1] .- vn[I1, I2.-1]) ./ (2.0*dy)
    uy = (un[I1, I2.+1] .- un[I1, I2.-1]) ./ (2.0*dy)
    vx = (vn[I1.+1, I2] .- vn[I1.-1, I2]) ./ (2.0*dx)

    maxDivU  = maximum(abs.(ux .+ vy))
    maxGradU = maximum(abs.(ux) .+ abs.(uy) .+ abs.(vx) .+ abs.(vy))

    return maxDivU, maxGradU, par
end
