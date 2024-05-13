# Developers Card
Collection of showcase applications that I have developed over the years. 

## Partial differential equation (PDE) solvers
### Porous convection
Porous convection is an important physical process that describes the dynamics of fluid extraction from an essentially undeforming solid in a multiphase aggregate. Applications include industrial filtration systems, subsurface hydrogeological reservoirs, or melt extraction below a volcano. The movie below shows porous convection across a 2D domain. To calculate the fluid motion through the pore space, Darcy's equation is coupled to the Heat equation. The solution is obtained iteratively applying a pseudo-transient time stepping method.

| Attribute           | Value                      |
| :------------------ | :------------------------- |
| Physics             | Darcy's law + Heat equation|
| Method              | central finite differences |
| Time integration    | fully implicit             |
| Advection scheme    | upwind (first order)       |
| Programmed language | Julia                      |
| Parallelization     | no                         |
| Rayleigh-No         | 1000                       |
| Porosity            | 10 % (const.)              |


https://github.com/lcandiot/DevelopersCard/assets/50524459/24bbff85-92bd-4f90-a380-bb2f09f45e02

### Analytical solutions
For benchmarking purposes it is useful to compare the numerical solution to analytical solutions. For some equations it is possible to derive a solution analytically for others they do not exist. In this section, I showcase a 1D Diffusion-Convection equation for which I have derived an analytical solution assuming an initial normal (or Gaussian) distribution of a scalar quantity and homogeneous Dirichlet boundary conditions. As shown in the movie below the numerical fits the analytical solution up to an acceptable tolerance.

https://github.com/lcandiot/DevelopersCard/assets/50524459/f39d5d7a-c08c-4a26-8761-3e0ef47f0344


## Data Analysis
### Fitting data using neural networks
The following figure shows an example of fitting the temperature-dependent rock density using neural networks. The linear regression algorithm was developed using the Flux.jl package.

![linearRegression_density_Tdependent](https://github.com/lcandiot/DevelopersCard/assets/50524459/4b08ad9a-981c-40b9-aa65-2f376c35a3f2)

## Numerical Methods
### Newton-Raphson iterative solution procedure
This technique allows to find the root of any function $f(x)$. Solution procedure starts by defining an initial guess ($x_n$), ideal already close to the solution. Determining the derivate $f'(x_n)$ and its root gives a new guess $x_{n+1}$ which is closer to the solution than the guess before. Repeat this procedure either until $|x_{n+1} - x_n| < \varepsilon$ or $f(x_{n+1}) \approx 0$. An advantage of the Newton-Raphson iteration is that convergence is quadratically in the vicinity of the root, reducing the computational cost. The example below shows the Newton-Raphson iterative solution procedure to find the root of $f(x) = x^2 - 2$. The derivative $f'(x_n)$ is calculated (1) analytically, (2) numerically as difference between $x_n$ and a perturbation of $x_n$ approximation, (3) via the Julia package **Zygote.jl** which enables automatic differentiation of any function. Note the scatter plot in the right panel indeed shows a quadratic convergence behaviour. 

https://github.com/lcandiot/DevelopersCard/assets/50524459/f3b32855-6646-4716-81c1-adf1f37e4232

This and similar methods are widely used in optimization algorithms and PDE solvers. I have also used this technique to find the liquidus temperature of geomaterials as a function of pressure and chemical composition. 
