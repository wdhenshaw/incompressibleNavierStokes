
# Try to load Plots.jl at include time; plotting is silently skipped if unavailable.
# Install with:  ] add Plots
const _HAVE_PLOTS = try; @eval(using Plots); true; catch; false; end

# Plot u, v, p, vorticity in one figure (Figure 1), and optionally errors (Figure 2).
# plotOption >= 2 enables plotting; movieMode=0 pauses between frames, movieMode=1 runs continuously.
function plot_solution(tn::Float64, un::Matrix{Float64}, vn::Matrix{Float64},
                       pn::Matrix{Float64}, par::Par)
    if !_HAVE_PLOTS
        @printf("plotSolution: Plots.jl not installed — skipping. Install with:  ] add Plots\n")
        return
    end

    cpu0 = time()

    I1, I2 = get_index(par.gid)
    xv = par.x[I1, par.gid[1,2], 1]   # x-coordinates (1D)
    yv = par.x[par.gid[1,1], I2, 2]   # y-coordinates (1D)

    ttl = @sprintf("INS: ts=%s t=%.2f step=%d nu=%g N=[%d,%d] cfl=%.2f bc=%s ad=%d",
                   par.ts, tn, par.step, par.nu, par.Ngx, par.Ngy, par.cfl, par.bcLabel, par.ad)

    # u, v, p, vorticity as top-view heatmaps (1x4 layout)
    pu = heatmap(xv, yv, un[I1,I2]', title="u", xlabel="x", ylabel="y",
                 color=:turbo, aspect_ratio=:equal, colorbar=true)
    pv = heatmap(xv, yv, vn[I1,I2]', title="v", xlabel="x", ylabel="y",
                 color=:turbo, aspect_ratio=:equal, colorbar=true)
    pp = heatmap(xv, yv, pn[I1,I2]', title="p", xlabel="x", ylabel="y",
                 color=:turbo, aspect_ratio=:equal, colorbar=true)

    dx = par.dx;  dy = par.dy
    vor = zeros(par.Ngx, par.Ngy)
    vor[I1, I2] = (un[I1, I2.+1] .- un[I1, I2.-1]) ./ (2.0*dy) .-
                  (vn[I1.+1, I2] .- vn[I1.-1, I2]) ./ (2.0*dx)
    pw = heatmap(xv, yv, vor[I1,I2]', title="vorticity", xlabel="x", ylabel="y",
                 color=:turbo, aspect_ratio=:equal, colorbar=true)

    fig1 = plot(pu, pv, pp, pw, layout=(1,4), plot_title=ttl, size=(1400, 380))
    display(fig1)

    # Error plots (only when an exact solution is available)
    if par.computeErrors != 0
        _, perr, uerr, verr, divf = get_errors(tn, un, vn, pn, par)
        pe1 = heatmap(xv, yv, uerr[I1,I2]', title="u-err", xlabel="x", ylabel="y", color=:turbo)
        pe2 = heatmap(xv, yv, verr[I1,I2]', title="v-err", xlabel="x", ylabel="y", color=:turbo)
        pe3 = heatmap(xv, yv, perr[I1,I2]', title="p-err", xlabel="x", ylabel="y", color=:turbo)
        pe4 = heatmap(xv, yv, divf[I1,I2]', title="div",   xlabel="x", ylabel="y", color=:turbo)
        fig2 = plot(pe1, pe2, pe3, pe4, layout=(2,2), plot_title=ttl, size=(900, 640))
        display(fig2)
    end

    if par.savePlots != 0
        mkpath(par.figDir)
        tl = replace(@sprintf("t%.1f", tn), "." => "p")
        savefig(fig1, joinpath(par.figDir, @sprintf("%s%sUVPVorticity.png", par.plotName, tl)))
    end

    par.cpuPlot += time() - cpu0
end
