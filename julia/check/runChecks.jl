
# Regression tests for ins.jl
# Usage:
#   julia runChecks.jl            -- compare against saved check files
#   julia runChecks.jl replace=1  -- run cases and save new check files

script_dir = @__DIR__
include(joinpath(script_dir, "..", "ins.jl"))

using Printf

function run_checks(args::String...)
    replace = false
    check_dir = script_dir

    for a in args
        m = match(r"^-?(\w+)=(.*)", a)
        m === nothing && continue
        k, v = m.captures[1], m.captures[2]
        k == "replace" && (replace = v ∈ ("1","true"))
    end

    cases = [
        ("insAB2PolyTZ.check",
         ["ts=ab2","tzScale=1","tf=0.1","ms=poly","knownSolution=none","nu=0.1",
          "degreex=2","degreet=2","bcs=nnnn","N0=10","idebug=0","plotOption=-1","computeErrors=1"]),
        ("insPC2PolyTZ.check",
         ["ts=pc2","tzScale=1","tf=0.1","ms=poly","knownSolution=none","nu=0.1",
          "degreex=2","degreet=2","bcs=nnnn","N0=10","idebug=0","plotOption=-1","computeErrors=1"]),
        ("insPC2PolyTZnsns.check",
         ["ts=pc2","tzScale=1","tf=0.1","ms=poly","knownSolution=none","nu=0.1",
          "degreex=2","degreet=2","bcs=nsns","N0=10","idebug=0","plotOption=-1","computeErrors=1"]),
        ("insPC2PolyTZionn.check",
         ["ts=pc2","tzScale=1","tf=0.1","ms=poly","knownSolution=none","nu=0.1",
          "degreex=2","degreet=2","bcs=ionn","N0=10","idebug=0","plotOption=-1","computeErrors=1"]),
        ("insIM2PolyTZ.check",
         ["ts=im2","tzScale=1","tf=0.1","ms=poly","knownSolution=none","nu=0.1",
          "degreex=2","degreet=2","bcs=nnnn","N0=10","idebug=0","plotOption=-1","computeErrors=1"]),
        ("insIM2PolyTZsnsn.check",
         ["ts=im2","tzScale=1","tf=0.1","ms=poly","knownSolution=none","nu=0.1",
          "degreex=2","degreet=2","bcs=snsn","N0=10","idebug=0","plotOption=-1","computeErrors=1"]),
        ("insIM2TaylorGreen.check",
         ["ts=im2","tf=0.2","ms=none","knownSolution=TaylorGreen","nu=0.1",
          "bcs=nnnn","N0=20","idebug=0","plotOption=-1","computeErrors=1"]),
    ]

    num_failed = 0
    rt = fill(-1, length(cases))

    for (icheck, (check_name, run_args)) in enumerate(cases)
        old_check = joinpath(check_dir, check_name)
        new_check = joinpath(pwd(), "ins.check")

        old_exists = isfile(old_check)
        if !old_exists || replace
            new_check = old_check
        end

        full_args = vcat(run_args, ["checkFileName=$new_check"])
        @printf("Running case %d: %s\n", icheck, check_name)
        ins(full_args...)

        if old_exists
            if read(new_check, String) == read(old_check, String)
                @printf("Test %d: %s : success.\n", icheck, check_name)
                rt[icheck] = 0
            else
                @printf("Test %d: %s : FAILED (files differ)\n", icheck, check_name)
                run(`diff $new_check $old_check`, wait=true)
                rt[icheck] = 1
                num_failed += 1
            end
        else
            @printf("Test %2d: %-35s : new file created.\n", icheck, check_name)
            rt[icheck] = 12345
        end
    end

    @printf("\n ----- check SUMMARY ----\n")
    for (icheck, (check_name, _)) in enumerate(cases)
        if rt[icheck] == 0
            @printf("Test %2d: %-35s : success.\n", icheck, check_name)
        elseif rt[icheck] == 12345
            @printf("Test %2d: %-35s : new file created.\n", icheck, check_name)
        else
            @printf("Test %2d: %-35s : FAILED.\n", icheck, check_name)
        end
    end
    if num_failed == 0
        @printf("===== SUCCESS. All regression tests passed. ======\n")
    else
        @printf("===== ERROR: %d regression tests FAILED. ======\n", num_failed)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_checks(ARGS...)
end
