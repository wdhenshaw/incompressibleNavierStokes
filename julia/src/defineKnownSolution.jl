
function define_known_solution!(par::Par)
    ks = par.knownSolution
    nu = par.nu
    kx = par.kx

    if ks == "TaylorGreen"
        par.computeErrors = 1
        ampu = 1.0;  ampp = ampu^2 / 4.0
        ftg  = t -> exp(-2.0*nu*kx^2*t)
        ftgt = t -> (-2.0*nu*kx^2)*exp(-2.0*nu*kx^2*t)

        par.ue  = (x,y,t) ->  ampu .* sin.(kx.*x) .* cos.(kx.*y) .* ftg(t)
        par.ve  = (x,y,t) -> -ampu .* cos.(kx.*x) .* sin.(kx.*y) .* ftg(t)
        par.pe  = (x,y,t) ->  ampp .* (cos.(2.0.*kx.*x) .+ cos.(2.0.*kx.*y)) .* ftg(t)^2

        par.uet = (x,y,t) ->  ampu .* sin.(kx.*x) .* cos.(kx.*y) .* ftgt(t)
        par.vet = (x,y,t) -> -ampu .* cos.(kx.*x) .* sin.(kx.*y) .* ftgt(t)

        par.pex = (x,y,t) ->  ampp .* (-2.0*kx) .* sin.(2.0.*kx.*x) .* ftg(t)^2
        par.pey = (x,y,t) ->  ampp .* (-2.0*kx) .* sin.(2.0.*kx.*y) .* ftg(t)^2

        par.uexy = (x,y,t) -> (-ampu*kx^2) .* cos.(kx.*x) .* sin.(kx.*y) .* ftg(t)
        par.vexy = (x,y,t) -> ( ampu*kx^2) .* sin.(kx.*x) .* cos.(kx.*y) .* ftg(t)

        par.ufe = (x,y,t) -> zeros(size(x))
        par.vfe = (x,y,t) -> zeros(size(x))
        par.pfe = (x,y,t) -> zeros(size(x))

    elseif ks == "Poiseuille"
        par.computeErrors = 1
        ampu = 1.0;  ampp = -nu*ampu*8.0

        par.ue  = (x,y,t) ->  (4.0.*ampu) .* (1.0 .- y) .* y
        par.ve  = (x,y,t) ->  zeros(size(x))
        par.pe  = (x,y,t) ->  ampp .* x

        par.uet = (x,y,t) -> zeros(size(x))
        par.vet = (x,y,t) -> zeros(size(x))

        par.pex = (x,y,t) ->  ampp .* ones(size(x))
        par.pey = (x,y,t) ->  zeros(size(x))

        par.uexy = (x,y,t) -> zeros(size(x))
        par.vexy = (x,y,t) -> zeros(size(x))

        par.ufe = (x,y,t) -> zeros(size(x))
        par.vfe = (x,y,t) -> zeros(size(x))
        par.pfe = (x,y,t) -> zeros(size(x))

    elseif ks == "none"
        if par.ms == "none"
            par.uet = (x,y,t) -> zeros(size(x))
            par.vet = (x,y,t) -> zeros(size(x))
            par.ufe = (x,y,t) -> zeros(size(x))
            par.vfe = (x,y,t) -> zeros(size(x))
            par.pfe = (x,y,t) -> zeros(size(x))
        end
    else
        error("define_known_solution: unknown knownSolution=$ks")
    end
end
