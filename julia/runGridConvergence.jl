
include("ins.jl")

using Printf

# Run a grid convergence study.
# Example:
#   julia runGridConvergence.jl ts=ab2 tf=0.25 ms=trig numResolutions=3
#
function run_grid_convergence(args::String...)
    # defaults
    num_res = 3; N0 = 10; ts = "ab2"; ms = "trig"; known = "none"; tf = 0.25; idebug = 0

    for a in args
        m = match(r"^-?(\w+)=(.*)", a)
        m === nothing && continue
        k, v = m.captures[1], strip(m.captures[2], [';',' '])
        k == "numResolutions" && (num_res = parse(Int, v))
        k == "N0"             && (N0      = parse(Int, v))
        k == "ts"             && (ts      = v)
        k == "ms"             && (ms      = v)
        k == "knownSolution"  && (known   = v)
        k == "tf"             && (tf      = parse(Float64, v))
        k == "idebug"         && (idebug  = parse(Int, v))
    end

    maxErr = zeros(4, num_res)
    cpuv   = zeros(num_res)
    Nv     = zeros(Int, num_res)

    for ires in 1:num_res
        Nx = N0 * 2^(ires-1)
        Nv[ires] = Nx
        res = ins("ts=$ts", "tzScale=1", "tf=$tf", "ms=$ms", "knownSolution=$known",
                  "idebug=$idebug", "nu=0.1", "bcs=nnnn", "N0=$Nx",
                  "plotOption=-1", "computeErrors=1")
        maxErr[:,ires] .= res.maxErr
        cpuv[ires] = res.cpu
    end

    @printf("  N     p-err  ratio   u-err  ratio   v-err  ratio   div    ratio    cpu(s)  ratio\n")
    @printf(" ---------------------------------------------------------------------------------\n")
    for ires in 1:num_res
        if ires == 1
            @printf("%4d  %8.2e       %8.2e       %8.2e       %8.2e        %8.2e\n",
                    Nv[ires], maxErr[1,ires], maxErr[2,ires], maxErr[3,ires], maxErr[4,ires], cpuv[ires])
        else
            @printf("%4d  %8.2e %4.1f  %8.2e %4.1f  %8.2e %4.1f  %8.2e %4.1f   %8.2e %4.1f\n",
                    Nv[ires],
                    maxErr[1,ires], maxErr[1,ires-1]/maxErr[1,ires],
                    maxErr[2,ires], maxErr[2,ires-1]/maxErr[2,ires],
                    maxErr[3,ires], maxErr[3,ires-1]/maxErr[3,ires],
                    maxErr[4,ires], maxErr[4,ires-1]/maxErr[4,ires],
                    cpuv[ires], cpuv[ires]/cpuv[ires-1])
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_grid_convergence(ARGS...)
end
