
# Specifying parameters for `ins.m`

The main routine used to solve the incompressible Navier-Stokes equations is executed by `ins.m`. Not only can `ins.m` be ran independently by a user, but it is also called by `runGridConvergence.m` and `runChecks.m`. In either case, there are a plethora of parameters that could be specified about `ins.m`, and this document attempts to detail all of them outside of the code.

[//]: <> ( # Metaparameters Before discussing the parameters pertaining to the function of `ins.m`, it is worth pointing out some parameters that change how the parameters in `ins.m` are interpreted, aptly named metaparameters.)

## Auto Assign

There are certain sets of parameters that would be nonsensical to specify one but not the others. Namely, in specifying solutions to the INS equations for the code to solve, it's important to specify certain parameters.

For instance, say you wanted to solve for a known solution `FreeSurfaceSolution` where the description of the free surface is vital for the problem. 

```
ins -knownSolution=FreeSurfaceSolution -map=freeSurface -motion=freeSurface
```


```
ins -knownSolution=FreeSurfaceSolution -aa=0
```

```
ins -knownSolution=FreeSurfaceSolution
```

```AUTO-ASSIGN
-map=freeSurface -motion=freeSurface
```


# Solution Description Parameters

To start, we need to tell `ins.m` what solution of the Incompressible Navier-Stokes it should be solving for. The code currently supports three mutually exclusive parameters of specifying a solution

- [`-knownSolution`](#known-solutions) : Specify a known solution to the INS equations (usually one with little forcing / well known forcing)
- [`-ms`](#manufactured-solutions) : Specify a manufactured solution to the INS equations (usually one where the forcing terms are easy to compute)
- [`-ic`](#initial-conditions) : Specify an initial velocity profile with no solution to compare the code against

These three ways of specifying a solution all work by specifying a string describing one of the solutions which the code will fill in information for later. If left unspecified, all three parameters will default to being `none`, a value that tells the code that the parameter is not being used.

## Known Solutions:

The `-knownSolution` parameter supports the following inputs:

- [`TaylorGreen`](#taylor-green-vortex) : Taylor-Green vortex exact solution
- [`Poiseuille`](#plane-poiseuille-flow) : Plane Poiseuille Flow
- [`GravityCapillaryWave`](#gravity-capillary-wave) : Gravity Capillary Travelling Wave


### Taylor Green Vortex

As detailed in `insDoc.pdf`, the exact solution for a Taylor Green vortex is given by
$$
\begin{align*}
    u &= + \sin (kx) \cos (ky) F(t) \\
    v &= - \cos (kx) \sin (ky) F(t) \\
    p &= \frac{\rho }4 ( \cos(2kx) + \cos(2ky) ) F^2(t) \\
    F(t) &\equiv \exp ( - 2 \nu k^2 t )
\end{align*}
$$
Note that while the `ins.m` code supports having different wave numbers in the horizontal & vertical length scale, the Taylor Green solution requires one wave number only, in which case the code will use `-kx`. Also note that $\rho$ and $\nu$ represent the density and diffusive coefficient of the fluid, which can also be specified by parameters (see [Physical Parameters](#physical-parameters)). Here is an example of running `ins.m` to solve for a Taylor-Green vortex solution:
```
ins -knownSolution=TaylorGreen -tf=1. -nu=0.1 -kx=1 -bcs=dddd -N0=100 -plotOption=1
```
### Plane Poiseuille Flow

As detailed in `insDoc.pdf`, Planar Poiseuille Flow is a steady-state solution to the INS equations that satisfies the following parabolic flow over the reference domain $\Omega = [0, 1]^2$
$$
\begin{align*}
    u &= \frac{ p_{\text{in}} - p_{\text{out}} }{2 \rho \nu} y (1 - y) \\
    v &= 0 \\
    p &= (p_{\text{out}} - p_{\text{in}}) x + p_{\text{in}}
\end{align*}
$$

```AUTO-ASSIGN
-xa=0 -xb=1 -ya=0 -yb=1 -bcs=Ionn -outflowPressureCoeff=1 -outflowPressureCoeffpn=0
```

### Gravity Capillary Wave

## Manufactured Solutions:

### Trigonometric

$$
\begin{align*}
    u &= \frac{k_y}{(k_x^2 + k_y^2)^q} \cos(k_x x) \cos(k_y y) \cos(k_t t) \\
    v &= \frac{k_x}{(k_x^2 + k_y^2)^q} \sin(k_x x) \sin(k_y y) \cos(k_t t) \\
    p &= \cos(k_x x) \sin(k_y y) \cos(k_t t)
\end{align*}
$$

### Polynomial

$$
\begin{align*}
    u &= q_u(x, y) r_u(t) \\
    v &= q_v(x, y) r_v(t) \\
    p &= q_p(x, y) r_p(t)
\end{align*}
$$

$$ q_* \in \mathcal P^2(\Omega) $$

## Initial Conditions:

### Zero

### Constant

### Shear Flow

### Perturbed Poiseuille Flow



# Boundary Conditions

# Physical Parameters

# Time & Timestepping Parameters

# Geometry & Grid Parameters

# Grid Motion Parameters

# Artificial Terms Parameters

# Optimization & Computation Parameters

# Printing Parameters

# Plotting Parameters




