
function define_manufactured_solution!(par::Par)
    ms      = par.ms
    degreex = par.degreex
    degreet = par.degreet
    kx      = par.kx
    ky      = par.ky
    kt      = par.kt
    nu      = par.nu
    tzScale = par.tzScale

    if ms == "none"
        return
    end

    par.computeErrors = 1

    if ms == "trig"

        if tzScale == 0
            ampu = ky / sqrt(kx^2 + ky^2)
            ampv = kx / sqrt(kx^2 + ky^2)
        else
            ampu = ky / (kx^2 + ky^2)
            ampv = kx / (kx^2 + ky^2)
        end
        ampp = 1.0

        par.ue  = (x,y,t) -> ampu .* cos.(kx.*x) .* cos.(ky.*y) .* cos(kt*t)
        par.ve  = (x,y,t) -> ampv .* sin.(kx.*x) .* sin.(ky.*y) .* cos(kt*t)
        par.pe  = (x,y,t) -> ampp .* cos.(kx.*x) .* sin.(ky.*y) .* cos(kt*t)

        par.uex = (x,y,t) -> ampu * (-kx) .* sin.(kx.*x) .* cos.(ky.*y) .* cos(kt*t)
        par.vex = (x,y,t) -> ampv * ( kx) .* cos.(kx.*x) .* sin.(ky.*y) .* cos(kt*t)
        par.pex = (x,y,t) -> ampp * (-kx) .* sin.(kx.*x) .* sin.(ky.*y) .* cos(kt*t)

        par.uey = (x,y,t) -> ampu * (-ky) .* cos.(kx.*x) .* sin.(ky.*y) .* cos(kt*t)
        par.vey = (x,y,t) -> ampv * ( ky) .* sin.(kx.*x) .* cos.(ky.*y) .* cos(kt*t)
        par.pey = (x,y,t) -> ampp * ( ky) .* cos.(kx.*x) .* cos.(ky.*y) .* cos(kt*t)

        par.uet = (x,y,t) -> ampu * (-kt) .* cos.(kx.*x) .* cos.(ky.*y) .* sin(kt*t)
        par.vet = (x,y,t) -> ampv * (-kt) .* sin.(kx.*x) .* sin.(ky.*y) .* sin(kt*t)
        par.pet = (x,y,t) -> ampp * (-kt) .* cos.(kx.*x) .* sin.(ky.*y) .* sin(kt*t)

        par.uexx = (x,y,t) -> ampu * (-kx^2) .* cos.(kx.*x) .* cos.(ky.*y) .* cos(kt*t)
        par.vexx = (x,y,t) -> ampv * (-kx^2) .* sin.(kx.*x) .* sin.(ky.*y) .* cos(kt*t)
        par.pexx = (x,y,t) -> ampp * (-kx^2) .* cos.(kx.*x) .* sin.(ky.*y) .* cos(kt*t)

        par.uexy = (x,y,t) -> ampu * ( kx*ky) .* sin.(kx.*x) .* sin.(ky.*y) .* cos(kt*t)
        par.vexy = (x,y,t) -> ampv * ( kx*ky) .* cos.(kx.*x) .* cos.(ky.*y) .* cos(kt*t)
        par.pexy = (x,y,t) -> ampp * (-kx*ky) .* sin.(kx.*x) .* cos.(ky.*y) .* cos(kt*t)

        par.ueyy = (x,y,t) -> ampu * (-ky^2) .* cos.(kx.*x) .* cos.(ky.*y) .* cos(kt*t)
        par.veyy = (x,y,t) -> ampv * (-ky^2) .* sin.(kx.*x) .* sin.(ky.*y) .* cos(kt*t)
        par.peyy = (x,y,t) -> ampp * (-ky^2) .* cos.(kx.*x) .* sin.(ky.*y) .* cos(kt*t)

        par.uett = (x,y,t) -> ampu * (-kt^2) .* cos.(kx.*x) .* cos.(ky.*y) .* cos(kt*t)
        par.vett = (x,y,t) -> ampv * (-kt^2) .* sin.(kx.*x) .* sin.(ky.*y) .* cos(kt*t)
        par.pett = (x,y,t) -> ampp * (-kt^2) .* cos.(kx.*x) .* sin.(ky.*y) .* cos(kt*t)

    elseif ms == "poly"

        cu0 = 1.0; cu1 = 0.0; cu2 = 0.0
        if degreet > 0  cu1 = 0.5  end
        if degreet > 1  cu2 = 0.25 end
        time_poly  = t -> (cu0 + t*(cu1 + t*cu2))
        time_polyt = t -> (cu1 + 2.0*cu2*t)
        time_polytt = t -> 2.0*cu2

        cu0_  = 1.0; cux = 0.0; cuy = 0.0; cuxx = 0.0; cuyy = 0.0; cuxy = 0.0
        if degreex > 0  cux = 1.0;  cuy  = 0.5  end
        if degreex > 1  cuxx = 1.0; cuyy = 1.0;  cuxy = 2.0 end

        par.ue  = (x,y,t) -> (cu0_ .+ cux.*x .+ cuy.*y .+ cuxx.*x.*x .+ cuxy.*x.*y .+ cuyy.*y.*y) .* time_poly(t)
        par.uex = (x,y,t) -> (cux .+ cuxy.*y .+ 2.0.*cuxx.*x) .* time_poly(t)
        par.uey = (x,y,t) -> (cuy .+ cuxy.*x .+ 2.0.*cuyy.*y) .* time_poly(t)
        par.uexx = (x,y,t) -> (2.0.*cuxx) .* time_poly(t) .* one.(x)
        par.uexy = (x,y,t) -> (cuxy) .* time_poly(t) .* one.(x)
        par.ueyy = (x,y,t) -> (2.0.*cuyy) .* time_poly(t) .* one.(x)
        par.uet  = (x,y,t) -> (cu0_ .+ cux.*x .+ cuy.*y .+ cuxx.*x.*x .+ cuxy.*x.*y .+ cuyy.*y.*y) .* time_polyt(t)
        par.uett = (x,y,t) -> (cu0_ .+ cux.*x .+ cuy.*y .+ cuxx.*x.*x .+ cuxy.*x.*y .+ cuyy.*y.*y) .* time_polytt(t)

        cv0 = 0.5; cvx = 0.0; cvy = 0.0; cvxx = 0.0; cvyy = 0.0; cvxy = 0.0
        if degreex > 0  cvx = 0.75; cvy = -cux  end
        if degreex > 1  cvxx = 1.0; cvyy = -0.5*cuxy; cvxy = -2.0*cuxx  end

        par.ve  = (x,y,t) -> (cv0 .+ cvx.*x .+ cvy.*y .+ cvxx.*x.*x .+ cvxy.*x.*y .+ cvyy.*y.*y) .* time_poly(t)
        par.vex = (x,y,t) -> (cvx .+ cvxy.*y .+ 2.0.*cvxx.*x) .* time_poly(t)
        par.vey = (x,y,t) -> (cvy .+ cvxy.*x .+ 2.0.*cvyy.*y) .* time_poly(t)
        par.vexx = (x,y,t) -> (2.0.*cvxx) .* time_poly(t) .* one.(x)
        par.vexy = (x,y,t) -> (cvxy) .* time_poly(t) .* one.(x)
        par.veyy = (x,y,t) -> (2.0.*cvyy) .* time_poly(t) .* one.(x)
        par.vet  = (x,y,t) -> (cv0 .+ cvx.*x .+ cvy.*y .+ cvxx.*x.*x .+ cvxy.*x.*y .+ cvyy.*y.*y) .* time_polyt(t)
        par.vett = (x,y,t) -> (cv0 .+ cvx.*x .+ cvy.*y .+ cvxx.*x.*x .+ cvxy.*x.*y .+ cvyy.*y.*y) .* time_polytt(t)

        cp0 = 1.0; cpx = 0.0; cpy = 0.0; cpxx = 0.0; cpxy = 0.0; cpyy = 0.0
        if degreex > 0  cpx = 0.50;  cpy  = 0.25  end
        if degreex > 1  cpxx = 1.0;  cpyy = 1.5;   cpxy = -0.5  end

        par.pe  = (x,y,t) -> cp0 .+ cpx.*x .+ cpy.*y .+ cpxx.*x.*x .+ cpxy.*x.*y .+ cpyy.*y.*y
        par.pex = (x,y,t) -> cpx .+ 2.0.*cpxx.*x .+ cpxy.*y
        par.pey = (x,y,t) -> cpy .+ 2.0.*cpyy.*y .+ cpxy.*x
        par.pexx = (x,y,t) -> 2.0.*cpxx .* one.(x)
        par.peyy = (x,y,t) -> 2.0.*cpyy .* one.(x)
        par.pexy = (x,y,t) -> cpxy .* one.(x)

    else
        error("define_manufactured_solution: unknown ms=$ms")
    end

    # MS forcing (composition of closures)
    par.ufe = (x,y,t) -> par.uet(x,y,t) .+ par.ue(x,y,t).*par.uex(x,y,t) .+ par.ve(x,y,t).*par.uey(x,y,t) .+ par.pex(x,y,t) .- nu.*(par.uexx(x,y,t) .+ par.ueyy(x,y,t))
    par.vfe = (x,y,t) -> par.vet(x,y,t) .+ par.ue(x,y,t).*par.vex(x,y,t) .+ par.ve(x,y,t).*par.vey(x,y,t) .+ par.pey(x,y,t) .- nu.*(par.vexx(x,y,t) .+ par.veyy(x,y,t))
    par.pfe = (x,y,t) -> par.pexx(x,y,t) .+ par.peyy(x,y,t) .+ (par.uex(x,y,t).*par.uex(x,y,t) .+ 2.0.*par.uey(x,y,t).*par.vex(x,y,t) .+ par.vey(x,y,t).*par.vey(x,y,t))
end
