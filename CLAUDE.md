# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

MATLAB research code solving the 2D incompressible Navier-Stokes equations (velocity-pressure formulation) using a pressure-Poisson projection method. Supports Cartesian and curvilinear (mapped) grids.

## Running the solver

All commands are MATLAB function calls. The solver is invoked as `ins` with keyword arguments:

```matlab
% Basic run (Adams-Bashforth, polynomial manufactured solution, no-slip walls)
ins -ts=ab2 -tf=.1 -ms=poly -nu=0.1 -bcs=nnnn -N0=10 -idebug=1

% Taylor-Green vortex (known solution for error checking)
ins -ts=pc2 -tf=.2 -ms=none -knownSolution=TaylorGreen -nu=0.1 -bcs=nnnn -N0=20 -idebug=1

% Kelvin-Helmholtz shear flow instability
ins -ts=pc2 -tf=2 -ms=none -knownSolution=none -nu=.0005 -ya=-.5 -yb=1.5 -xb=2 -bcs=ppnn -N0=500 -ic=shear -movieMode=1 -ad=1 -ad21=.1 -ad22=.1 -cfl=0.95

% Curvilinear grid (Annulus)
ins -ts=ab2 -tf=.1 -ms=poly -nu=0.05 -bcs=nnnn -N0=5 -map=Annulus -idebug=1
```

Key parameters:
- `-ts=` : time-stepping scheme — `ab2` (Adams-Bashforth 2), `pc2` (predictor-corrector 2), `im2` (implicit 2)
- `-ms=` : manufactured solution — `none`, `poly`, `trig`
- `-bcs=` : four-character BC string, one per boundary (left, right, bottom, top): `d`=Dirichlet, `n`=no-slip, `s`=slip, `i`=inflow, `I`=pressure inflow, `o`=outflow, `p`=periodic
- `-map=` : grid mapping — `Cartesian` (default), `Rectangle` (curvilinear identity, for testing), `Annulus`
- `-knownSolution=` : `none`, `TaylorGreen`, `Poiseuille`

## Regression tests

```matlab
% Run all regression checks (compare against saved .check files)
runChecks

% Re-generate check files after intentional changes
runChecks -replace=1
```

Check files live in `check/`. Each `.check` file stores norms and max errors for a specific configuration. `src/check.m` drives all cases; `src/writeCheckFile.m` writes them.

To add a new regression case: add an `ins` call inside `src/check.m`, run `runChecks -replace=1` to save the baseline, then commit the new `.check` file.

## Grid convergence verification

```matlab
runGridConvergence -ts=ab2 -tf=.25 -ms=trig -numResolutions=3
runGridConvergence -ts=im2 -tf=1.0 -ms=none -knownSolution=TaylorGreen -numResolutions=4
```

Prints a table of max errors (p, u, v, div) and convergence ratios across doubling resolutions. Expected convergence: ~2nd order for `im2`/`ab2`, ~4th order for `pc2` on smooth manufactured solutions.

## Architecture

### Central data structure: `par`

Everything flows through a single MATLAB struct `par`. It is initialized in `ins.m` and passed by value into every function, which returns an updated copy. It holds grid geometry, parameters, metric derivatives, timing accumulators, BC type integers, and factored matrices.

### Grid indexing

Arrays have ghost points. For a grid of `N0` interior cells:
- `par.gid(side,axis)` gives the boundary indices: `gid(1,1)=iax`, `gid(2,1)=ibx`, etc.
- Interior range: `iax:ibx` × `iay:iby`
- Ghost points extend one cell beyond each boundary (for `orderInSpace=2`, `numGhost=1`)
- `getIndex(par.gid)` returns `[I1,I2]` index vectors for boundary-to-boundary range
- `getIndexInterior(par.gid)` returns strictly interior indices

### Curvilinear grids

For non-Cartesian grids, `setupGrid` calls `evalMap` (in `src/evalMap.m`) which populates:
- `par.x(i1,i2,1:2)` — physical coordinates
- `par.rx(i1,i2,1:2,1:2)` — metric derivatives ∂rᵢ/∂xⱼ (inverse Jacobian)
- `par.dr(1:2)` — computational grid spacings (always 1/N in mapped space)

Physical derivatives are formed using the chain rule: e.g., `u.x = rx*u.r + sx*u.s`. `getDerivatives.m` handles this for second-order accurate first and second derivatives at a single point.

### Time-stepping loop (in `ins.m`)

Each time step calls one of:
- `advanceAdams` — explicit AB2, calls `getUt` then `pressureEquation`
- `advancePC` — predictor-corrector: two stages, each with `getUt` + `pressureEquation`
- `advanceIM` — implicit: solves `solveImplicitTimeStep` (sparse linear system per velocity component) then `pressureEquation`

After each advance, BCs are applied via `applyBoundaryConditions`.

### Pressure solve

`pressureEquation.m` builds and solves a sparse Poisson system. The matrix is factored once (`factorMatrix=1`) and reused across time steps when the time step is constant. The RHS is `-(u.x² + 2v.x·u.y + v.y²)` plus optional divergence damping (`cdv` coefficient). Ghost-point BCs for pressure are Neumann (normal gradient = forcing) on no-slip/slip walls, Dirichlet on inflow/outflow.

### Implicit velocity solve

`formImplicitTimeSteppingMatrix.m` + `solveImplicitTimeStep.m` build and solve `(I - dt*nu*L)u = rhs` as a sparse system. When boundary conditions mix types (slip walls, pressure inflow), `par.multipleImplicitSolversNeeded=1` and separate matrices are formed for each velocity component.

### Manufactured solutions and known solutions

- `defineManufacturedSolution.m` — sets up `par.fu`, `par.fv`, `par.fp` forcing functions as anonymous functions
- `defineKnownSolution.m` — sets `par.uExact`, `par.vExact`, `par.pExact` for TaylorGreen/Poiseuille
- `defineBoundaryForcingFunctions.m` — derives boundary values from the manufactured solution
- Errors are computed in `getErrors.m` and returned in `outPar.maxErr(1:4)` = [p, u, v, div]
