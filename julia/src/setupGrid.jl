
# Construct the Cartesian grid and fill par.x[i1,i2,1:2]
function setup_grid!(par::Par)
    if par.Nx < 0
        par.Nx = par.N0
    end
    if par.Ny < 0
        par.Ny = max(1, round(Int, par.Nx * (par.yb - par.ya) / (par.xb - par.xa)))
    end

    par.dx = (par.xb - par.xa) / par.Nx
    par.dy = (par.yb - par.ya) / par.Ny
    par.dr[1] = par.dx
    par.dr[2] = par.dy

    numGhost = par.orderInSpace ÷ 2
    par.numGhost = numGhost

    iax = 1 + numGhost;  iay = 1 + numGhost
    ibx = iax + par.Nx;  iby = iay + par.Ny
    Ngx = ibx + numGhost
    Ngy = iby + numGhost

    par.gid[1,1] = iax;  par.gid[2,1] = ibx
    par.gid[1,2] = iay;  par.gid[2,2] = iby

    par.dim[1,1] = 1;  par.dim[2,1] = Ngx
    par.dim[1,2] = 1;  par.dim[2,2] = Ngy

    par.Ngx = Ngx
    par.Ngy = Ngy

    par.isCartesian = true

    par.x = zeros(Ngx, Ngy, 2)
    for iy in 1:Ngy
        for ix in 1:Ngx
            par.x[ix, iy, 1] = par.xa + (ix - iax) * par.dx
            par.x[ix, iy, 2] = par.ya + (iy - iay) * par.dy
        end
    end

    # Cartesian metrics: not needed but keep rx as placeholder for code that may check it
    par.rx = zeros(0, 0, 2, 2)
end
