
# Compute max-norm errors vs exact solution and divergence.
# Returns (maxErr[4], perr, uerr, verr, div)
#   maxErr = [pErrMax, uErrMax, vErrMax, divMax]
function get_errors(tn::Float64,
                    un::Matrix{Float64}, vn::Matrix{Float64}, pn::Matrix{Float64},
                    par::Par)
    dx = par.dx; dy = par.dy
    I1, I2 = get_index(par.gid)

    uTrue = par.ue(par.x[:,:,1], par.x[:,:,2], tn)
    vTrue = par.ve(par.x[:,:,1], par.x[:,:,2], tn)
    pTrue = par.pe(par.x[:,:,1], par.x[:,:,2], tn)

    uerr = un .- uTrue
    verr = vn .- vTrue
    perr = pn .- pTrue

    uErrMax = maximum(abs.(un[I1,I2] .- uTrue[I1,I2]))
    vErrMax = maximum(abs.(vn[I1,I2] .- vTrue[I1,I2]))
    pErrMax = maximum(abs.(pn[I1,I2] .- pTrue[I1,I2]))

    div = zeros(par.Ngx, par.Ngy)
    ux = (un[I1.+1,I2] .- un[I1.-1,I2]) ./ (2.0*dx)
    vy = (vn[I1,I2.+1] .- vn[I1,I2.-1]) ./ (2.0*dy)
    div[I1,I2] .= ux .+ vy

    divMax = maximum(abs.(div[I1,I2]))

    maxErr = [pErrMax, uErrMax, vErrMax, divMax]
    return maxErr, perr, uerr, verr, div
end
