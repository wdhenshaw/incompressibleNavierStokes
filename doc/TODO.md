
# TODO 

 - [x] Add Documentation to different `ins.m` parameters
    - Not added in `ins-parameter.md`, but detailed in comments of `ins.m`
 - [x] Condense certain runs of `ins.m` to minimal specified parameters
   - Added `assignDefaults.m` file which changes default parameters depending on certain specified parameters
- [x] Getting quadratic convergence for the `GravityCapillaryWave` solution
    - [x] Solving a motionless, cartesian problem with fake dirichlet conditions
    - [x] Solving a problem with free surface motion and geometry, but a fake dirichlet BC on the surface
    - [ ] Solving for free surface motion & geometry, with proper boundary conditions (`ppnt`)

### 10.08.26

- [ ] Fix the traction boundary condition
    - Note, Henshaw doesn't seem to have proof of traction BC converging quadratically with grid refinement
    - [ ] Investigate how it works for other free surface solutions
    - [ ] Investigate if pressure is the sole source of error
- [ ] Fixing the `im2` timestepping scheme
- [ ] Implementing a Travelling Wave Solution for an Annulus
