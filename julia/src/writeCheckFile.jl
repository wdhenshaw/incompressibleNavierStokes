
using Printf

function write_check_file(par::Par)
    open(par.checkFileName, "w") do f
        @printf(f, "--------------------- INS Check File --------------------------\n")
        @printf(f, "ms=%s, ts=%s, bcs=%s, knownSolution=%s\n", par.ms, par.ts, par.bcs, par.knownSolution)
        @printf(f, "nu=%9.2e, cfl=%9.2e, tFinal=%9.2e, cdv=%9.2e\n", par.nu, par.cfl, par.tFinal, par.cdv)
        @printf(f, "ad=%d, ad21=%9.2e, ad22=%9.2e\n", par.ad, par.ad21, par.ad22)
        @printf(f, "Nx=%d Ny=%d Nt=%d dt=%8.2e\n", par.Nx, par.Ny, par.Nt, par.dt)
        @printf(f, "uNorm=%8.2e, vNorm=%8.2e, pNorm=%8.2e\n", par.uNorm[1], par.uNorm[2], par.uNorm[3])
        @printf(f, "maxDivU=%8.2e, maxGradU=%8.2e\n", par.maxDivU, par.maxGradU)
        if par.computeErrors != 0
            @printf(f, "max-Err(p,u,v)=(%8.2e,%8.2e,%8.2e)\n", par.maxErr[1], par.maxErr[2], par.maxErr[3])
        end
    end
    if par.idebug > 0
        @printf("wrote checkFile=[%s]\n", par.checkFileName)
    end
end
