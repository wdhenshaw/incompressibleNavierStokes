Matlab codes to solve the incompressible Navier Stokes equations

ins.m                : main routine
runChecks.m          : run regression checks in the check folder, check/check.m
runGridConvergence.m : run a grid refinement study



# SAMPLE RUNS:

## TAYLOR-GREEN vortex
```
   ins -ts=im2 -tf=1 -tp=0.1 -ms=none -knownSolution=TaylorGreen -nu=0.1  -bcs=nnnn -N0=40 -idebug=1 -plotOption=1 -movieMode=1 -dtMax=1e-2;
   ins -ts=pc2 -tf=1 -tp=0.1 -ms=none -knownSolution=TaylorGreen -nu=0.1  -bcs=nnnn -N0=40 -idebug=1 -plotOption=1 -movieMode=1;  
```

## SHEAR-FLOW -- Kelvin-Helmholtz instability
```
  -- finer grid N=500
  ins -ts=pc2 -tzScale=1 -tf=2 -tp=0.1 -ms=none -knownSolution=none -idebug=1 -nu=.0005 -ya=-.5 -yb=1.5 -xb=2 -bcs=ppnn -N0=500 -ic=shear -movieMode=1 -plotOption=3 -plotErrors=0 -kx=1 -checkTimeStep=20 -ad=1 -ad21=.1 -ad22=.1 -cfl=0.95 -plotVorticity=1 -numThreads=1;

  -- coarser grid N=300
  ins -ts=pc2 -tzScale=1 -tf=2 -tp=.05 -ms=none -knownSolution=none -idebug=1 -nu=.0005 -ya=-.5 -yb=1.5 -xb=2 -bcs=ppnn -N0=300 -ic=shear -movieMode=1 -plotOption=3 -plotErrors=0 -kx=1 -checkTimeStep=20 -ad=1 -ad21=.1 -ad22=.1 -cfl=0.95 -plotVorticity=1 -numThreads=1;
```


## MOVING GRID EXAMPLES:
```
  -- translating square + poly TZ 
    ins -ts=ab2 -tf=0.5 -ms=poly -degreex=1 -degreet=1 -kx=0.5 -ky=0.5 -knownSolution=none -nu=0.05 -cdv=0 -bcs=dddd -idebug=1 -plotGrid=1 -motion=translate -plotOption=1 -N0=5;

  -- rotating square + trig TZ
     ins -ts=im2 -tf=1 -ms=trig -kx=0.5 -ky=0.5 -knownSolution=none -nu=0.05 -cdv=0 -bcs=nnnn -idebug=1 -plotEveryStep=0 -plotGrid=1 -motion=rotate -movieMode=1 -N0=20
```


## FREE-SURFACE 

```
 -- sin(kx*x) initial surface
   ins -ts=pc2 -tf=1 -tp=0.05 -ms=none -kx=2 -knownSolution=none -nu=0.01 -bcs=ppst -idebug=1 -movieMode=1 -motion=freeSurfaceMotion -map=freeSurface -ampfs=0.05 -plotEveryStep=0 -plotOption=3 -gravity=-10 -ic=zero -plotGrid=1 -plotAspectRatio=2.5 -checkTimeStep=5 -N0=40;

 -- Gaussian initial surface
   PC2: 
   ins -ts=pc2 -tf=1 -tp=0.05 -ms=none -kx=2 -knownSolution=none -nu=0.01 -cdv=0 -bcs=ppst -idebug=1 -movieMode=1 -motion=freeSurfaceMotion -map=freeSurface -icfs=gaussian -ampfs=0.15 -plotEveryStep=0 -plotOption=3 -gravity=-10 -ic=zero -plotGrid=1 -plotAspectRatio=2.5 -checkTimeStep=5 -N0=40;

   IM2: 
   ins -ts=im2 -tf=1 -tp=0.05 -ms=none -kx=2 -knownSolution=tone -nu=0.01 -cdv=0 -bcs=ppst -idebug=1 -movieMode=1 -motion=freeSurfaceMotion -map=freeSurface -icfs=gaussian -ampfs=0.15 -plotEveryStep=0 -plotOption=3 -gravity=-10 -ic=zero -plotGrid=1 -plotAspectRatio=2.5 -checkTimeStep=5 -dtMax=1e-2 -N0=40;
```
