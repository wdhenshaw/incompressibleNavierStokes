# Specifying parameters for `ins.m`

The main routine used to solve the incompressible Navier-Stokes equations is executed by `ins.m`. Not only can `ins.m` be ran independently by a user, but it is also called by `runGridConvergence.m` and `runChecks.m`. In either case, there are a plethora of parameters that could be specified about `ins.m`, and this document attempts to detail all of them outside of the code.

[//]: <> ( # Metaparameters Before discussing the parameters pertaining to the function of `ins.m`, it is worth pointing out some parameters that change how the parameters in `ins.m` are interpreted, aptly named metaparameters.)

## Auto Assign

There are certain sets of parameters that would be nonsensical to specify one but not the others. Namely, in specifying solutions to the INS equations for the code to solve, it's important to specify certain parameters. To address this, the code utilizes a "metaparameter", Auto Assign [`-aa`](#auto-assign) that is turned on by default to add necessary parameters.

For instance, say you wanted to solve for a known solution `FreeSurfaceSolution` where the description of the free surface is vital for the problem. 

```
ins -knownSolution=FreeSurfaceSolution -map=freeSurface -motion=freeSurface
```

Without Auto Assign turned on, the parameters `-map` and `-motion` would default to `Cartesian` and `none` respectively. But in regard to the `FreeSurfaceSolution` known solution, it would make more sense for the default values to be `freeSurface`.
```
ins -knownSolution=FreeSurfaceSolution -aa=0
```

So instead, by the code defaulting Auto Assign to be on, when ever the `FreeSurfaceSolution` known solution is ran, the code will automatically add the necessary parameters to make the grid have free motion
```
ins -knownSolution=FreeSurfaceSolution
```
throughout the documentation, the specific parameters that Auto Assign adds are specified by code blocks starting with `(AA)` as shown here:
```
(AA) -map=freeSurface -motion=freeSurface
```


# Solution Description Parameters

To start, we need to tell `ins.m` what solution of the Incompressible Navier-Stokes it should be solving for. The code currently supports three mutually exclusive parameters of specifying a solution

- [`-knownSolution`](#known-solutions) : Specify a known solution to the INS equations (usually one with little forcing / well known forcing)
- [`-ms`](#manufactured-solutions) : Specify a manufactured solution to the INS equations (usually one where the forcing terms are easy to compute)
- [`-ic`](#initial-conditions) : Specify an initial velocity profile with no solution to compare the code against

These three ways of specifying a solution all work by specifying a string describing one of the solutions which the code will fill in information for later. If left unspecified, all three parameters will default to being `none`, a value that tells the code that the parameter is not being used.

## Known Solutions

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

As detailed in `insDoc.pdf`, Planar Poiseuille Flow is a steady-state solution to the INS equations that satisfies the following boundary value problem over the reference domain $\Omega = [0, 1]^2$ for some specified inlet pressure $p_{\text{in}}$ and outlet pressure

$$
\begin{align*}
    \mathbf u(x, 0, t) = \mathbf u(x, 1, t) &= 0 \\
    p(0, y, t) = p_{\text{in}}, \qquad v(0, y, t) &\equiv 0\\
    p(1, y, t) = p_{\text{out}}, \qquad v(1, y, t) &\equiv 0
\end{align*}
$$

The velocity and pressure profile that solve the boundary value problem are given by the following parabolic flow, and when running the code to solve for the `Poiseuille` known solution, the code will reference these exact solutions.

$$
\begin{align*}
    u &= \frac{ p_{\text{in}} - p_{\text{out}} }{2 \rho \nu} y (1 - y) \\
    v &= 0 \\
    p &= (p_{\text{out}} - p_{\text{in}}) x + p_{\text{in}}
\end{align*}
$$

Note that the inlet and outlet pressure values can be changed by specifying `-pressureInflowValue` and `-pOutflow` parameters, regardless if the code is using a Pressure Inlet or Outlet boundary condition. 

Using Auto Assign, when solving for the Planar Poiseuille Flow, the following parameters are specified by default

```
(AA) -xa=0 -xb=1 -ya=0 -yb=1 -bcs=Ionn -outflowPressureCoeff=1 -outflowPressureCoeffpn=0
```

### Gravity Capillary Wave

As detailed in `insDoc.pdf`, given a wave number $k$, 

$$
\begin{align*}
    \eta(x, t) &= \eta_0 \cos(k x - \omega t)
\end{align*}
$$

$$
\begin{align*}
    \varphi (x, y, t) &= \frac{g}{\omega} ( 1 + k^2 \text{Bo}^{-1} ) \frac{\cosh(k (y + 1))}{\cosh(k)} \cdot - \eta_0 \sin(kx - \omega t)
\end{align*}
$$

$$
\begin{align*}
    \mathbf u &= \nabla \cdot \varphi \\
    p &= - \rho (\varphi_t + g y)
\end{align*}
$$

```
(AA) -xa=0 -xb=1 -ya=-1 -yb=0 -bcs=ppnt -map=freeSurface -motion=freeSurfaceMotion -icfs=cos -ampfs=1e-8 -gravity=-1
```

## Manufactured Solutions

The `-ms` parameters

- [`trig`](#trigonometric)
- [`poly`](#polynomial)

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

$$ \begin{align*} q_* \in \mathcal P^2(\Omega) \end{align*} $$

## Initial Conditions:

- `default`
- `zero`
- `constant`
- `perturbedPoiseuille`

### Constant (Zero)

$$
\begin{align*}
    u(x, y, 0) &\equiv u_{\text{ic}} \\
    v(x, y, 0) &\equiv v_{\text{ic}}
\end{align*}
$$

### Shear Flow

$$ 
\begin{align*}
    u(x, y, 0) &= \tanh(\beta(y - y_m)) \\
    v(x, y, 0) &= \delta \cos(2 \pi (y - y_m)) \sin (2 \pi k_x x)
\end{align*}
$$

### Perturbed Poiseuille Flow

$$
\begin{align*}
    U(y) &= \frac{ p_{\text{in}} - p_{\text{out}} }{2 \rho \nu} y (1 - y)
\end{align*}
$$


$$
\begin{align*}
    u(x, y, 0) &= U(y) + \varepsilon u_p(x, y) \\
    v(x, y, 0) &= v_p(x, y)
\end{align*}
$$

$$
\begin{align*}
    u_p(x, y) &=  \cos( k x ) \sin( k y) \\
    v_p(x, y) &= -\sin( k x) \cos( k y)
\end{align*}
$$

# Boundary Conditions

## No-Slip Wall and Inflow

## Slip Wall

## Outflow

## Pressure Inflow

## Traction (Free Surface)

# Physical Parameters

## Appearing in the PDE

- `-rho`
- `-nu`
- `-gravity`
- `-gamma`

## Wave Numbers

- `-kx`
- `-ky`
- `-kt`
- `-tzScale`

## `-perturbation`

# Time & Timestepping Parameters

- `-tf`
- `-dtMax`
- `-cfl`
- `-checkTimeStep`
- `-ts`

## Adams-Bashforth 2
## Predictor-Corrector 2
## Implicit-Explicit Scheme

# Geometry & Grid Parameters

## Cartesian Map

## Rectangle Map

## Annulus Map

## Transfinite-Interpolation Map

## Rotated Square Map

## Free Surface Map

# Grid Motion Parameters

## Translational Motion

## Rotational Motion

## Deforming Motion

## Free Surface Motion

# Artificial Terms Parameters

## Divergence Damping

## Artificial Dissipation

# Optimization & Computation Parameters


# Printing Parameters

## `-idebug`

# Plotting Parameters

