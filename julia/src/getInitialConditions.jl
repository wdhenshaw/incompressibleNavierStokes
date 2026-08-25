
function get_initial_conditions!(t::Float64, un::Matrix{Float64}, vn::Matrix{Float64}, par::Par)
    if par.ic == "default"
        un .= par.u0(par.x[:,:,1], par.x[:,:,2])
        vn .= par.v0(par.x[:,:,1], par.x[:,:,2])

    elseif par.ic == "zero"
        fill!(un, 0.0)
        fill!(vn, 0.0)

    elseif par.ic == "constant"
        fill!(un, par.uic)
        fill!(vn, par.vic)

    elseif par.ic == "shear"
        par.plotErrors = 0
        beta = par.shearBeta
        ym   = 0.5*(par.ya + par.yb)
        un  .= tanh.(beta .* (par.x[:,:,2] .- ym))
        delta = par.shearDeltav
        vn  .= delta .* cos.(2.0.*pi.*(par.x[:,:,2] .- ym)) .* sin.(par.kx .* par.x[:,:,1])

    else
        error("get_initial_conditions: unknown ic=$(par.ic)")
    end
end
