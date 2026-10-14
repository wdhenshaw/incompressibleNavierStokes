ADD option to compute relative errors

Let's have a look at the GravityCapillaryWave:

```
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.
--- gravityCapillary: add PERIODIC AND TRACTION TO TOP  ERRORS ARE VERY LARGE!!


pc2: t=9.843e-04 step=     2 Nx= 40 dt=4.921e-04 maxRelErr-[p,u,v]=[0.00e+00,9.22e-02,4.01e-03] div/grad= 1.49e-03 (grad=4.22e-05), cpu= 1.50e+00(s)
pc2: t=9.843e-04 step=     8 Nx= 80 dt=1.230e-04 maxRelErr-[p,u,v]=[0.00e+00,1.32e-01,6.56e-03] div/grad= 3.88e-04 (grad=4.15e-05), cpu= 2.98e+00(s)


ins -N0=40 -idebug=1 -useTrueForP=1 -computeErrors=1 -computeRelativeErrors=1 -plotOption=3 -plotEveryStep=1 -ts=pc2 -cdv=0 -tf=.5 -tp=.01 -ms=none -knownSolution=GravityCapillaryWave -bcs=ppdt -gravity=-1 -ampfs=1e-4 -map=freeSurface -motion=freeSurfaceMotion 

USE TRUE FOR P: BIG ERRORS FOR (u,v) ON THE TRACTION BC

pc2: t=4.916e-04 step=     1 Nx= 40 dt=4.916e-04 maxRelErr-[p,u,v]=[0.00e+00,4.93e-02,1.26e-04] div/grad= 1.60e-03 (grad=4.33e-03), cpu= 1.50e+00(s)

pc2: t=9.833e-03 step=    20 Nx= 40 dt=4.916e-04 maxRelErr-[p,u,v]=[0.00e+00,4.18e-01,7.07e-02] div/grad= 1.26e-03 (grad=3.44e-03), cpu= 2.20e+00(s)


pc2: t=9.833e-03 step=    20 Nx= 40 dt=4.916e-04 maxRelErr-[p,u,v]=[3.41e-04,4.57e-01,1.23e-01] div/grad= 1.20e-02 (grad=3.24e-03), cpu= 2.59e+00(s)
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

--- gravityCapillary: DIRICHLET BC's : errors decrease BETTER when CDV==0 
ins -N0=40 -idebug=1 -computeErrors=1 -computeRelativeErrors=1 -plotOption=3 -plotEveryStep=0 -ts=pc2 -cdv=0 -tf=.5 -tp=.01 -ms=none -knownSolution=GravityCapillaryWave -bcs=dddd -gravity=-1 -ampfs=1e-4 -map=freeSurface -motion=freeSurfaceMotion 
cdv==0
pc2: t=9.833e-03 step=    20 Nx= 40 dt=4.916e-04 maxRelErr-[p,u,v]=[4.00e-07,1.33e-04,6.47e-05] div/grad= 1.57e-03 (grad=4.44e-03), cpu= 1.68e+00(s)
pc2: t=9.963e-03 step=    81 Nx= 80 dt=1.230e-04 maxRelErr-[p,u,v]=[1.09e-07,2.26e-05,1.87e-05] div/grad= 4.13e-04 (grad=4.45e-03), cpu= 1.09e+01(s)

cvd=1:
pc2: t=9.833e-03 step=    20 Nx= 40 dt=4.916e-04 maxRelErr-[p,u,v]=[9.64e-06,1.69e-03,1.07e-03] div/grad= 8.42e-03 (grad=4.46e-03), cpu= 1.82e+00(s)
pc2: t=9.963e-03 step=    81 Nx= 80 dt=1.230e-04 maxRelErr-[p,u,v]=[3.87e-06,5.73e-04,2.04e-04] div/grad= 2.85e-03 (grad=4.46e-03), cpu= 1.06e+01(s)


gravityCapillary: plot every step

ins -N0=40 -idebug=1 -computeErrors=1 -computeRelativeErrors=1 -plotOption=3 -plotEveryStep=1 -ts=pc2 -cdv=0 -tf=.5 -tp=.001 -ms=none -knownSolution=GravityCapillaryWave -bcs=dddd -gravity=-1 -ampfs=1e-4 -map=freeSurface -motion=freeSurfaceMotion

>>> N=80
pc2: t=4.916e-04 step=     1 Nx= 40 dt=4.916e-04 maxRelErr-[p,u,v]=[3.99e-07,1.86e-05,5.31e-06] div/grad= 1.67e-03 (grad=4.44e-03), cpu= 1.19e+00(s)
*** FORM THE PRESSURE MATRIX : Ng=1849 t=  9.8328e-04 cur=3 **
pc2: t=9.833e-04 step=     2 Nx= 40 dt=4.916e-04 maxRelErr-[p,u,v]=[3.99e-07,3.48e-05,1.11e-05] div/grad= 1.66e-03 (grad=4.44e-03), cpu= 1.83e+00(s)
pc2: t=1.475e-03 step=     3 Nx= 40 dt=4.916e-04 maxRelErr-[p,u,v]=[3.99e-07,4.85e-05,1.53e-05] div/grad= 1.66e-03 (grad=4.44e-03), cpu= 2.52e+00(s)
pc2: t=1.967e-03 step=     4 Nx= 40 dt=4.916e-04 maxRelErr-[p,u,v]=[3.99e-07,6.01e-05,1.80e-05] div/grad= 1.65e-03 (grad=4.44e-03), cpu= 3.25e+00(s)
pc2: t=2.458e-03 step=     5 Nx= 40 dt=4.916e-04 maxRelErr-[p,u,v]=[3.99e-07,7.00e-05,2.02e-05] div/grad= 1.64e-03 (grad=4.44e-03), cpu= 3.93e+00(s)
pc2: t=2.950e-03 step=     6 Nx= 40 dt=4.916e-04 maxRelErr-[p,u,v]=[3.99e-07,7.86e-05,2.38e-05] div/grad= 1.63e-03 (grad=4.44e-03), cpu= 4.48e+00(s)
pc2: t=3.441e-03 step=     7 Nx= 40 dt=4.916e-04 maxRelErr-[p,u,v]=[3.99e-07,8.60e-05,2.72e-05] div/grad= 1.63e-03 (grad=4.44e-03), cpu= 5.03e+00(s)
pc2: t=3.933e-03 step=     8 Nx= 40 dt=4.916e-04 maxRelErr-[p,u,v]=[3.99e-07,9.26e-05,3.05e-05] div/grad= 1.62e-03 (grad=4.44e-03), cpu= 5.56e+00(s)
pc2: t=4.425e-03 step=     9 Nx= 40 dt=4.916e-04 maxRelErr-[p,u,v]=[3.99e-07,9.84e-05,3.36e-05] div/grad= 1.62e-03 (grad=4.44e-03), cpu= 6.12e+00(s)
pc2: t=4.916e-03 step=    10 Nx= 40 dt=4.916e-04 maxRelErr-[p,u,v]=[4.00e-07,1.03e-04,3.66e-05] div/grad= 1.61e-03 (grad=4.44e-03), cpu= 6.86e+00(s)

>>>> N=80 
pc2: t=1.230e-04 step=     1 Nx= 80 dt=1.230e-04 maxRelErr-[p,u,v]=[1.09e-07,2.30e-06,5.65e-07] div/grad= 4.36e-04 (grad=4.45e-03), cpu= 1.52e+00(s)
*** FORM THE PRESSURE MATRIX : Ng=6889 t=  2.4600e-04 cur=3 **
pc2: t=2.460e-04 step=     2 Nx= 80 dt=1.230e-04 maxRelErr-[p,u,v]=[1.09e-07,4.30e-06,1.23e-06] div/grad= 4.35e-04 (grad=4.45e-03), cpu= 2.40e+00(s)
pc2: t=3.690e-04 step=     3 Nx= 80 dt=1.230e-04 maxRelErr-[p,u,v]=[1.09e-07,5.99e-06,1.71e-06] div/grad= 4.34e-04 (grad=4.45e-03), cpu= 3.24e+00(s)
pc2: t=4.920e-04 step=     4 Nx= 80 dt=1.230e-04 maxRelErr-[p,u,v]=[1.09e-07,7.44e-06,2.02e-06] div/grad= 4.33e-04 (grad=4.45e-03), cpu= 3.94e+00(s)
pc2: t=6.150e-04 step=     5 Nx= 80 dt=1.230e-04 maxRelErr-[p,u,v]=[1.09e-07,8.68e-06,2.21e-06] div/grad= 4.32e-04 (grad=4.45e-03), cpu= 4.62e+00(s)
pc2: t=7.380e-04 step=     6 Nx= 80 dt=1.230e-04 maxRelErr-[p,u,v]=[1.09e-07,9.76e-06,2.34e-06] div/grad= 4.32e-04 (grad=4.45e-03), cpu= 5.32e+00(s)
pc2: t=8.610e-04 step=     7 Nx= 80 dt=1.230e-04 maxRelErr-[p,u,v]=[1.09e-07,1.07e-05,2.45e-06] div/grad= 4.31e-04 (grad=4.45e-03), cpu= 6.05e+00(s)
pc2: t=9.840e-04 step=     8 Nx= 80 dt=1.230e-04 maxRelErr-[p,u,v]=[1.09e-07,1.15e-05,2.53e-06] div/grad= 4.30e-04 (grad=4.45e-03), cpu= 6.75e+00(s)
pc2: t=1.107e-03 step=     9 Nx= 80 dt=1.230e-04 maxRelErr-[p,u,v]=[1.09e-07,1.22e-05,2.80e-06] div/grad= 4.30e-04 (grad=4.45e-03), cpu= 7.46e+00(s)
pc2: t=1.230e-03 step=    10 Nx= 80 dt=1.230e-04 maxRelErr-[p,u,v]=[1.09e-07,1.29e-05,3.08e-06] div/grad= 4.29e-04 (grad=4.45e-03), cpu= 8.13e+00(s)
pc2: t=1.353e-03 step=    11 Nx= 80 dt=1.230e-04 maxRelErr-[p,u,v]=[1.09e-07,1.35e-05,3.34e-06] div/grad= 4.29e-04 (grad=4.45e-03), cpu= 8.83e+00(s)


>>>>>>>>>>>>>>
```

me❯ I am going to pull the latest version from the git repository. Changes were made by another person on the project. I would like you to help me go over the changes.

>claude:
I've fetched the four incoming commits and read the code diffs. Now I'm checking whether they change existing behavior, mainly defaults and the regression suite.

  Ran 1 shell command

I fetched origin/main but haven't pulled or merged anything. The other person is Gabriel Weredyk, who pushed 4 commits between Sep 24 and Oct 8:

┌─────────────────────────────────────────────────────────┬───────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│                         Commit                          │                                                   What it does                                                    │
├─────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ 59142a7 Starting to incorporate Gravity Capillary Waves │ Adds the matlabUtilities/ folder to the repo and a new GravityCapillaryWave known solution                        │
├─────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ 5dac3ac Added Perturbed Poiseuille Flow ICs             │ Generalizes the Poiseuille solution, adds -ic=perturbedPoiseuille, and adds doc figures and tables                │
├─────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ c8baefa Start of ins-parameters.md                      │ Reorganizes the parameter block in ins.m, renames README to README.md, adds doc/ins-parameters.md and doc/TODO.md │
├─────────────────────────────────────────────────────────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ 6867d34 Start of debugging Traction                     │ Adds the "auto-defaults" mechanism (src/assignDefaults.m) and changes how runGridConvergence builds its commands  │
└─────────────────────────────────────────────────────────┴───────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘

A large part of the diff is LaTeX build output (.aux, .log, .fls, .synctex.gz, the PDF) and generated tables. The real code changes are in ins.m, runGridConvergence.m and about 9 files in src/.

Before you pull

git pull will abort. These files are untracked on your side but tracked upstream, so git won't overwrite them:
- ins.check
- doc/insDoc.fdb_latexmk
- doc/insDoc.fls
- doc/insDoc.synctex.gz

They're all generated files, so you can move or delete them first. The src/*Old.m files and check/ins.check don't conflict.

The main changes

1. GravityCapillaryWave known solution (src/defineKnownSolution.m)
- It's linear potential-flow theory for a small free-surface wave over depth H, using the dispersion relation ω² = gk(1 + γk²/ρg)·tanh(kH).
- I checked the math: the kinematic condition, u = φₓ, the velocity derivatives and the linearized Bernoulli pressure are all consistent.
- The amplitude is about 1e-8, so the nonlinear terms are negligible.
- One limitation: a potential flow satisfies the momentum equation for any ν, but it doesn't satisfy the viscous stress condition at the free surface unless ν = 0. With nu=0.1 this isn't an exact solution. That may explain the reduced order the author notes ("linear convergence in the pressure"), and why the last commit is "debugging Traction".
- To make it run, pressureEquation.m now turns off boundary forcing for this case.

2. New parameter handling (ins.m, src/assignDefaults.m)
- Command-line arguments are now read into a separate temp_par.
- assignDefaults then applies per-case presets and merges the user's arguments on top. For example, -knownSolution=GravityCapillaryWave automatically sets the domain, bcs=ppnt, map=freeSurface, the free-surface motion, gravity=-1, and so on.
- The defaults block in ins.m was regrouped by category.

3. Poiseuille changes
- The exact solution is now driven by the pressure drop: ampu = (pIn − pOut)/(2μ), with p = pIn at x = 0.
- It now defaults to bcs=Ionn (pressure inflow, outflow).

4. runGridConvergence now passes along only the flags you type, instead of building a fixed ins command.

5. Some defaults now live inside the functions that use them: degreex, degreet, uic, vic, shearBeta and shearDeltav.

Problems I found

1. The default BC changed from nnnn to dddd in ins.m. Any run without -bcs= now behaves differently, including the examples in CLAUDE.md. All 18 regression cases in check/check.m set -bcs, so they aren't affected.
2. runGridConvergence no longer passes its own defaults on to ins. It still sets ms='trig', bcs='dddd', tzScale=1 and nu=0.1 locally, but those reach ins only if you type them. For example, runGridConvergence -ts=ab2 now runs with ms=none, which has no errors to measure. The LaTeX table is still named from the local values, so the table name can disagree with what was actually run.
3. The -aa flag doesn't do anything. It's described as auto-assign on/off, but assignDefaults never checks it, so the presets always apply.
4. A typo in assignDefaults.m: it sets par.outflowPressureCoeff=1, but the real parameter is outflowPressureCoeffp.
5. assignCommandLineOption.m:
   - length(n1)>0 && length(n1)>0 tests n1 twice; the second test should be on n2.
   - Unknown parameter names are now accepted silently. echo=0 hides the warning, and the value is assigned anyway.
6. -icfs=cos only works with this one solution. getGrid.m uses par.DispersionRelation, which only the GravityCapillaryWave case defines, so any other use of -icfs=cos will crash.
7. Plotting changes that affect everyone:
   - setAspectRatio.m now returns immediately.
   - plotSolution.m turns off figure resizing with if(1==0).

   These look like personal preferences, not fixes.
8. Smaller items:
   - ins.m now calls addpath('src/') with a relative path, so it only works when started from the repo root.
   - A top-level ins.check (just the output of the last run) and the LaTeX build files are now tracked. They will keep causing conflicts, so a .gitignore would help.
   - The README's last example now has -knownSolution=tone, which is a typo.

What's fine

- None of the regression baselines in check/ changed.
- The parameter reorganization doesn't change any values except the bcs default.
- par.mu is defined (ins.m:214) before the new code uses it.

MATLAB R2023b is installed, so after you pull I can run runChecks to confirm nothing broke. I can also draft notes for Gabriel about items 1–7, or fix them myself.


