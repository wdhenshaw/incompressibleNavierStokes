
using Printf

# include all source files
include("src/Par.jl")
include("src/getIndex.jl")
include("src/getBoundaryIndex.jl")
include("src/getAdjustedBoundaryIndex.jl")
include("src/getIndexInterior.jl")
include("src/setupGrid.jl")
include("src/defineManufacturedSolution.jl")
include("src/defineKnownSolution.jl")
include("src/defineBoundaryForcingFunctions.jl")
include("src/getInitialConditions.jl")
include("src/getTimeStep.jl")
include("src/getMaxDivergence.jl")
include("src/getUt.jl")
include("src/applyBoundaryConditions.jl")
include("src/pressureEquation.jl")
include("src/formImplicitTimeSteppingMatrix.jl")
include("src/solveImplicitTimeStep.jl")
include("src/advanceAdams.jl")
include("src/advancePC.jl")
include("src/advanceIM.jl")
include("src/getErrors.jl")
include("src/plotSolution.jl")
include("src/writeCheckFile.jl")

# Parse a "-key=value" or "key=value" arg and set the matching Par field.
function assign_command_line_option!(line::String, par::Par)
    m = match(r"^-?([A-Za-z_][A-Za-z0-9_]*)=(.*)", strip(line))
    m === nothing && return
    key = m.captures[1]
    val_str = strip(m.captures[2], [';', ' '])
    sym = Symbol(key)
    hasproperty(par, sym) || return
    T = typeof(getproperty(par, sym))
    if T == Int
        setproperty!(par, sym, parse(Int, val_str))
    elseif T == Float64
        setproperty!(par, sym, parse(Float64, val_str))
    elseif T == String
        setproperty!(par, sym, String(val_str))
    elseif T == Bool
        setproperty!(par, sym, val_str ∈ ("1","true"))
    end
end

# Parse the bcs string and fill par.bc[side,axis].
# bcs character order: left(1,1), right(2,1), bottom(1,2), top(2,2)
function parse_bcs!(par::Par)
    bc_map = Dict('d'=>DIRICHLET,'n'=>NO_SLIP_WALL,'s'=>SLIP_WALL,
                  'I'=>PRESSURE_INFLOW,'i'=>INFLOW,'o'=>OUTFLOW,'p'=>PERIODIC)
    par.bcLabel = ""
    par.multipleImplicitSolversNeeded = false
    for axis in 1:2, side in 1:2
        m = side + 2*(axis-1)
        m > length(par.bcs) && continue
        c = par.bcs[m]
        haskey(bc_map, c) || error("Unknown BC character '$c' in bcs=$(par.bcs)")
        par.bc[side,axis] = bc_map[c]
        par.bcLabel *= string(c)
        if c == 's' || c == 'I'
            par.multipleImplicitSolversNeeded = true
        end
    end
end

# Main solver function.  Call as: res = ins("ts=ab2", "ms=trig", ...)
# Returns NamedTuple (maxErr, cpu).
function ins(args::String...)
    par = Par()
    for a in args
        assign_command_line_option!(a, par)
    end

    if length(args) == 0
        println("Usage: ins(\"ts=ab2\", \"tf=0.25\", \"ms=trig\", \"bcs=nnnn\", \"N0=10\", ...)")
        return (maxErr=zeros(4), cpu=0.0)
    end

    par.tFinal = par.tf

    # scale wave numbers by 2π
    par.kx *= 2π;  par.ky *= 2π;  par.kt *= 2π

    par.cdv = (par.ms != "none" && par.map != "Cartesian") ? 0.0 : par.cdv

    parse_bcs!(par)

    # reset timings
    par.cpuSetup           = 0.0
    par.cpuGetUt           = 0.0
    par.cpuPressure        = 0.0
    par.cpuBC              = 0.0
    par.cpuImplicit        = 0.0
    par.cpuFactorImpMatrix = 0.0
    par.cpuPlot            = 0.0

    define_manufactured_solution!(par)
    define_known_solution!(par)
    define_boundary_forcing_functions!(par)

    cpu_setup_start = time()

    setup_grid!(par)

    if par.idebug > 0
        @printf("-----------------Incompressible Navier-Stokes ---------------------------\n")
        @printf(" ts=%s, tFinal=%g, nu=%g, cfl=%g, cdv=%g, knownSolution=%s, N0=%d, idebug=%d\n",
                par.ts, par.tFinal, par.nu, par.cfl, par.cdv, par.knownSolution, par.N0, par.idebug)
        @printf(" ms=%s, degreex=%d, degreet=%d, [kx,ky,kt]=[%g,%g,%g]*2pi, tzScale=%d\n",
                par.ms, par.degreex, par.degreet,
                par.kx/(2π), par.ky/(2π), par.kt/(2π), par.tzScale)
        @printf(" bcLabel=%s, bc=[%d,%d,%d,%d], gid=[%d,%d,%d,%d]\n",
                par.bcLabel, par.bc[1,1], par.bc[2,1], par.bc[1,2], par.bc[2,2],
                par.gid[1,1], par.gid[2,1], par.gid[1,2], par.gid[2,2])
        @printf(" ad=%d, ad21=%g, ad22=%g\n", par.ad, par.ad21, par.ad22)
        @printf("--------------------------------------------------------------------------\n")
    end

    Ngx = par.Ngx; Ngy = par.Ngy
    un   = zeros(Ngx, Ngy)
    vn   = zeros(Ngx, Ngy)
    pn   = zeros(Ngx, Ngy)
    unp1 = zeros(Ngx, Ngy)
    vnp1 = zeros(Ngx, Ngy)
    pnp1 = zeros(Ngx, Ngy)
    ut   = zeros(Ngx, Ngy)
    vt   = zeros(Ngx, Ngy)

    cpu_total_start = time()
    par.cpuSetup = cpu_total_start

    t = 0.0
    get_initial_conditions!(t, un, vn, par)

    par.dtOld = -1.0
    dt = get_time_step(0, un, vn, par)
    Nt = max(2, ceil(Int, par.tFinal / dt))
    dt = par.tFinal / Nt

    dtOld = dt
    par.ab1 =  dt * (1.0 + dt/(2.0*dtOld))   # → 1.5*dt when dt==dtOld
    par.ab2 = -dt *       dt/(2.0*dtOld)       # → -0.5*dt

    # Factor pressure matrix and solve at t=0
    pn .= pressure_equation!(t, un, vn, dt, true, par)

    par.nuScaleFactor = 1.0

    # Bootstrap du/dt from the previous step (t = -dt)
    if par.ts == "ab2" || par.ts == "pc2"
        t_boot = -dt
        ut .= par.uet(par.x[:,:,1], par.x[:,:,2], t_boot)
        vt .= par.vet(par.x[:,:,1], par.x[:,:,2], t_boot)
    elseif par.ts == "im2"
        t_boot = -dt
        unp1 .= par.ue(par.x[:,:,1], par.x[:,:,2], t_boot)
        vnp1 .= par.ve(par.x[:,:,1], par.x[:,:,2], t_boot)
        pnp1 .= par.pe(par.x[:,:,1], par.x[:,:,2], t_boot)
        ut, vt = get_ut(t_boot, unp1, vnp1, pnp1, 0.0, par)
    end

    if par.ts == "im2"
        form_implicit_time_stepping_matrix!(dt, par)
    end

    par.cpuSetup           = time() - par.cpuSetup
    par.cpuGetUt           = 0.0
    par.cpuPressure        = 0.0
    par.cpuBC              = 0.0
    par.cpuImplicit        = 0.0
    par.cpuFactorImpMatrix = 0.0

    t = 0.0
    nextTimeToPrint = par.tp
    nextTimeToPlot  = par.tp
    tnp1 = 0.0
    n = 0
    maxNumberOfSteps = 1_000_000

    # Plot initial condition
    if par.plotOption >= 4
        plot_solution(t, un, vn, pn, par)
        par.movieMode == 0 && (print("Press Enter to continue..."); readline())
    end

    for n_step in 1:maxNumberOfSteps
        n = n_step
        t    = tnp1
        tnp1 = t + dt
        par.step = n

        if par.ts == "ab2"
            ut, vt = advance_adams!(t, dt, un, vn, pn, unp1, vnp1, pnp1, ut, vt, par)
        elseif par.ts == "pc2"
            ut, vt = advance_pc!(t, dt, un, vn, pn, unp1, vnp1, pnp1, ut, vt, par)
        elseif par.ts == "im2"
            ut, vt = advance_im!(t, dt, un, vn, pn, unp1, vnp1, pnp1, ut, vt, par)
        else
            error("ins: unknown ts=$(par.ts)")
        end

        # swap arrays
        un, unp1 = unp1, un
        vn, vnp1 = vnp1, vn
        pn, pnp1 = pnp1, pn

        if par.idebug > 0 && tnp1 >= nextTimeToPrint - 0.5*dt
            nextTimeToPrint += par.tp
            cpu_cur = time() - cpu_total_start
            maxDivU, maxGradU, _ = get_max_divergence(un, vn, par)
            divOverGrad = maxDivU / max(maxGradU, 1e-10)
            if par.computeErrors != 0
                maxErr, _, _, _, _ = get_errors(tnp1, un, vn, pn, par)
                @printf("%s: t=%9.3e step=%6d dt=%9.3e err-[u,v,p]=[%8.2e,%8.2e,%8.2e] div/grad=%9.2e cpu=%9.2e(s)\n",
                        par.ts, tnp1, n, dt, maxErr[2], maxErr[3], maxErr[1], divOverGrad, cpu_cur)
            else
                @printf("%s: t=%9.3e step=%6d dt=%9.3e div/grad=%9.2e cpu=%9.2e(s)\n",
                        par.ts, tnp1, n, dt, divOverGrad, cpu_cur)
            end
        end

        tnp1 > par.tFinal - 0.5*dt && break

        # Mid-run plots: every step OR at tp intervals
        at_plot_time = tnp1 >= nextTimeToPlot - 0.5*dt
        do_plot = par.plotEveryStep == 1 ||
                  (at_plot_time && (par.plotOption >= 2 || par.movieMode == 1))
        if do_plot
            at_plot_time && (nextTimeToPlot += par.tp)
            plot_solution(tnp1, un, vn, pn, par)
            par.movieMode == 0 && (print("Press Enter to continue..."); readline())
        end

        if mod(n, par.checkTimeStep) == 0
            dtNew = get_time_step(n, un, vn, par)
            numStepsRemaining = max(1, ceil(Int, (par.tFinal - tnp1) / dtNew))
            dtNew = (par.tFinal - tnp1) / numStepsRemaining
            dtDiff = abs(dtNew - dt) / dt
            if dtDiff > 0.1
                @printf("Change time-step: step=%d dt=%10.4e dtNew=%10.4e (rel-diff=%9.2e)\n",
                        n, dt, dtNew, dtDiff)
                dtOld = dt; dt = dtNew
                par.ab1 =  dt*(1.0 + dt/(2.0*dtOld))
                par.ab2 = -dt*      dt/(2.0*dtOld)
                if par.ts == "im2"
                    @printf("time-step changed — refactoring implicit matrix...\n")
                    form_implicit_time_stepping_matrix!(dt, par)
                end
            end
        else
            par.ab1 = 1.5*dt
            par.ab2 = -0.5*dt
        end
    end

    # Final plot
    if par.plotOption >= 2
        plot_solution(tnp1, un, vn, pn, par)
        (print("Press Enter to continue..."); readline())
    end

    par.cpuTotal = time() - cpu_total_start
    Nt_actual = n
    par.Nt = Nt_actual
    par.dt = dt

    maxErr = zeros(4)
    if par.computeErrors != 0
        maxErr, _, _, _, _ = get_errors(tnp1, un, vn, pn, par)
        if par.plotOption >= 0
            @printf("%s: t=%8.2e Nx=%3d Ny=%3d Nt=%5d dt=%8.2e max-Err(p,u,v)=(%8.2e,%8.2e,%8.2e) cpu=%8.2e(s)\n",
                    par.ts, tnp1, par.Nx, par.Ny, Nt_actual, dt,
                    maxErr[1], maxErr[2], maxErr[3], par.cpuTotal)
        end
    else
        if par.plotOption >= 0
            @printf("%s: t=%8.2e Nx=%3d Ny=%3d Nt=%5d dt=%8.2e cpu=%8.2e(s)\n",
                    par.ts, tnp1, par.Nx, par.Ny, Nt_actual, dt, par.cpuTotal)
        end
    end

    if par.plotOption != -1
        @printf("  ------------- TIMINGS %s Nx=%3d Ny=%3d Nt=%5d --------------------\n",
                par.ts, par.Nx, par.Ny, Nt_actual)
        @printf("                  cpu (s)     %%\n")
        @printf("Total             %8.2e  %5.1f\n", par.cpuTotal, 100.0)
        @printf("  setup           %8.2e  %5.1f\n", par.cpuSetup, par.cpuSetup/par.cpuTotal*100)
        @printf("  getUt           %8.2e  %5.1f\n", par.cpuGetUt, par.cpuGetUt/par.cpuTotal*100)
        @printf("  implicit solve  %8.2e  %5.1f\n", par.cpuImplicit, par.cpuImplicit/par.cpuTotal*100)
        @printf("  implicit factor %8.2e  %5.1f\n", par.cpuFactorImpMatrix, par.cpuFactorImpMatrix/par.cpuTotal*100)
        @printf("  pressure solve  %8.2e  %5.1f\n", par.cpuPressure, par.cpuPressure/par.cpuTotal*100)
        @printf("  bc              %8.2e  %5.1f\n", par.cpuBC, par.cpuBC/par.cpuTotal*100)
    end

    par.maxErr = maxErr
    par.uNorm[1] = maximum(abs.(un))
    par.uNorm[2] = maximum(abs.(vn))
    par.uNorm[3] = maximum(abs.(pn))
    maxDivU, maxGradU, _ = get_max_divergence(un, vn, par)
    par.maxDivU  = maxDivU
    par.maxGradU = maxGradU

    write_check_file(par)

    return (maxErr=maxErr, cpu=par.cpuTotal)
end

# Allow running as a script: julia ins.jl ts=ab2 ms=trig ...
if abspath(PROGRAM_FILE) == @__FILE__
    ins(ARGS...)
end
