# Analytical solution of the diffusion-convection equation 1D
It is demonstrated how to derive an analytical solution for the diffusion-convection equation in 1D. This PDE is one of the most important equations in Fluid Mechanics, as it describes the evolution of chemical species in a fluid, the distribution of heat in a system, or the diffusion of a fluid through porous media. We will largely follow the derivation given in [1].

The transient diffusion of any quantity $u(x, t)$ in the absence of source terms can be described mathematically as 

$$
\frac{\partial u}{\partial t} = D  \ \frac{\partial^2 u}{\partial x^2} \ - V_x \frac{\partial u}{\partial x},
$$

where $t$ is time, $x$ is coordinate, $D$ is some diffusion coefficient, and $V_x$ is velocity of the medium. The number of boundary conditions necessary to solve the equation is defined by the highest-order derivative, in this case the second-order derivative of the diffusion term. Hence, two conditions are required to describe the interaction between the internal part and the spatial boundary of the system. Since the change in $u$ over time is described by a first-order derivate, one initial condition is required.

## Boundary conditions
Most commonly three different types of boundary conditions are used in fluid dynamics:

* **Robin** boundary condition
    $$
    BC =
    \begin{cases}
        -k\frac{\partial u}{\partial x} (0, t) = -h[u(0, t) - u_{0, surr}] \\
        -k\frac{\partial u}{\partial x} (L, t) = h[u(L, t) - u_{L, surr}]
    \end{cases} \ ,
    $$
* **Dirichlet** boundary condition
    $$
    BC =
        \begin{cases}
        u(0, t) = u_{0, surr} \\
        u(L, t) = u_{L, surr}
    \end{cases} \ ,
    $$
* **Neumann** boundary condition
    $$
    BC =
    \begin{cases}
        \frac{\partial u}{\partial x} (0, t) = 0 \\
        \frac{\partial u}{\partial x} (L, t) = 0
    \end{cases} \ ,
    $$

where $BC$ denotes the boundary condition, $k$ is related to the diffusion coefficient, and $h$ is a transfer coefficient.

**Robin** BC can be used to set the flux at the boundary according to the interaction of the system ($u$) with its surroundings ($u_{surr}$). Both **Dirichlet** and **Neumann** BC can be considered special cases of the **Robin** BC. 

Dividing the above **Robin** BC by $h$ and assuming that the diffusion time scale is much lower than the transfer time scale ($k << h$)

$$
    \lim_{h \to \infty} \frac{k}{h} = 0
$$

yields the **Dirichlet** BC. This means that transfer of $u$ from the surroundings into the system is fast enough to keep $u$ constant at the boundaries.

Similarly, one can divide the **Robin** BC by $k$ and assume that transfer happens at much slower rates than diffusion ($h << k$)

$$
\lim_{k \to \infty} \frac{h}{k} = 0
$$

yielding the **Neumann** BC. In this case, the diffusion process in the system happens faster than transfer between the system and its surrounding allowing for changes of $u$ at the boundary. In the special case 

$$
u_{0, surr} = u_{L, surr} = 0
$$

the BC are called homogeneous.

## Initial conditions
Finally, the first derivative of $u$ in time requires definition of an initial condition describing the distribution of $u$ at the beginning of observation

$$
u(x, 0) = f(x) \ .
$$

## Problem statement
Let us consider the diffusion-convection equation

\begin{align}
\frac{\partial u}{\partial t} & = D  \ \frac{\partial^2 u}{\partial x^2} \ - V_x \frac{\partial u}{\partial x}, & (1)
\end{align}


where $D$ and $V_x$ are constants. We apply the following initial and boundary conditions

\begin{align}
u(x, 0) & = f(x) & (2) \\
u(0, t) & = u(L, t) = 0 & (3) \\
\frac{\partial u}{\partial x} (0, t) - a_0 u(0, t) & = \frac{\partial u}{\partial x} (L, t) + a_L u(L, t)  = 0 & (4) \\
\frac{\partial u}{\partial x} (0, t) & = \frac{\partial u}{\partial x} (L, t) = 0 \ , &(5)
\end{align}

where $(3)$ - $(5)$ denote **Dirichlet**, **Robin**, and **Neumann** BC respectively, and $a_0 > 0$ and $a_L > 0$.

## Eliminate the convection term
Let's first transform the diffusion-convection equation into a diffusion equation

\begin{align}
u(x, t) & = A(x, t) v(x, t) \ . & (6)
\end{align}

The goal here is to find $A(x, t)$ such that the convection term $-V_x\frac{\partial u}{\partial x}$ disappears. Substitute this formulation back into the original equation

\begin{align}
\frac{\partial (Av)}{\partial t} & = D \frac{\partial}{\partial x} \Bigg( \frac{\partial (Av)}{\partial x} \Bigg) - V_x \frac{\partial (Av)}{\partial x} & (7)
\end{align}

and open the brackets using product rule

\begin{align}
v\frac{\partial A}{\partial t} + A\frac{\partial v}{\partial t} &= D \frac{\partial}{\partial x} \Bigg( v\frac{\partial A}{\partial x} + A\frac{\partial v}{\partial x} \Bigg) - V_x \Bigg(v\frac{\partial A}{\partial x} + A \frac{\partial v}{\partial x} \Bigg)\\
\Leftrightarrow v\frac{\partial A}{\partial t} + A\frac{\partial v}{\partial t} & = D\Bigg( \frac{\partial v}{\partial x}\frac{\partial A}{\partial x} +v\frac{\partial^2 A}{\partial x^2} + \frac{\partial A}{\partial x}\frac{\partial v}{\partial x} + A\frac{\partial^2 v}{\partial x^2}\Bigg) - V_x \Bigg(v\frac{\partial A}{\partial x} + A \frac{\partial v}{\partial x} \Bigg)\\
\Leftrightarrow v\frac{\partial A}{\partial t} + A\frac{\partial v}{\partial t} & = D\Bigg( v\frac{\partial^2 A}{\partial x^2} + 2\frac{\partial A}{\partial x}\frac{\partial v}{\partial x} + A\frac{\partial^2 v}{\partial x^2}\Bigg) - V_x \Bigg(v\frac{\partial A}{\partial x} + A \frac{\partial v}{\partial x} \Bigg) & (8)
\end{align}

now divide by A and collect similar terms in front of $v$

\begin{align}
\frac{v}{A}\frac{\partial A}{\partial t} + \frac{\partial v}{\partial t} & = D \Bigg(\frac{v}{A}\frac{\partial^2 A}{\partial x^2} + \frac{2}{A}\frac{\partial A}{\partial x}\frac{\partial v}{\partial x} + \frac{\partial^2 v}{\partial x^2} \Bigg) - V_x\Bigg( \frac{v}{A}\frac{\partial A}{\partial x} + \frac{\partial v}{\partial x} \Bigg) \\
\Leftrightarrow \frac{\partial v}{\partial t} & = D \frac{\partial^2 v}{\partial x^2} + \Bigg( \frac{2D}{A}\frac{\partial A}{\partial x} - V_x \Bigg) \frac{\partial v}{\partial x} + \Bigg( -\frac{1}{A}\frac{\partial A}{\partial t} + \frac{D}{A}\frac{\partial^2 A}{\partial x^2} - \frac{V_x}{A}\frac{\partial A}{\partial x} \Bigg)v \ . & (9)
\end{align}


To yield a standard diffusion equation for $v$, we must satisfy

\begin{align}
\frac{2D}{A}\frac{\partial A}{\partial x} - V_x & = 0 & (10)\\
-\frac{1}{A}\frac{\partial A}{\partial t} + \frac{D}{A}\frac{\partial^2 A}{\partial x^2} - \frac{V_x}{A}\frac{\partial A}{\partial x} & = 0 \ . & (11)
\end{align}

Let's reformulate $(10)$ as

\begin{align}
\frac{1}{A}\frac{\partial A}{\partial x} & = \frac{V_x}{2D} \ . & (12)
\end{align}

Note that (12) is actually an ordinary differential equation (ODE) which we can easily solve for

\begin{align}
 \frac{1}{A}\frac{d A}{d x} & =\frac{V_x}{2D} \\
\Leftrightarrow \int \frac{1}{A} dA & = \int \frac{V_x}{2D} dx \\
\Leftrightarrow \ln(A) & = \frac{V_x}{2D} x + C_1(t) \\
\Leftrightarrow A(x, t) & = C_2(t)\exp\Bigg(\frac{V_x}{2D}x \Bigg) \ , & (13)
\end{align}

where $C_2(t) = \exp(C_1(t))$ is an integration constant that does not change with respect to $x$, but is changing in $t$ because $A(x, t)$. Hence, we must find $C_2(t)$ by solving $(11)$. To this end, we compute the partial derivatives using $(13)$

\begin{align}
\frac{\partial A}{\partial t} = \frac{\partial C_2}{\partial t} \exp\Bigg( \frac{V_x}{2D} x\Bigg) & = \frac{1}{C_2} \frac{\partial C_2}{\partial t} A & (14) \\
\frac{\partial A}{\partial x} = C_2 \exp\Bigg( \frac{V_x}{2D} x \Bigg) \frac{V_x}{2D} & = A \frac{V_x}{2D} & (15) \\
\frac{\partial^2 A}{\partial x^2} = C_2 \exp\Bigg( \frac{V_x}{2D} x \Bigg) \frac{V_x^2}{4D^2} & = A \frac{V_x^2}{4D^2} \ . & (16) \\
\end{align}

Substitute $(14)-(16)$ into $(11), rearrange and simplify to yield

\begin{align}
-\frac{1}{C_2}\frac{\partial C_2}{\partial t} + \frac{V_x^2}{4D} - \frac{V_x^2}{2D} & = 0 \\
\frac{1}{C_2}\frac{\partial C_2}{\partial t} & = -\frac{V_x^2}{4D} \ . & (17)
\end{align}

Note that this is again an ODE which can be solved similar to $(12)$ to find

\begin{align}
C_2(t) & = C_3\exp\Bigg( -\frac{V_x^2}{4D} t \Bigg) \ . & (18) 
\end{align}

Next, we can substitute $(18)$ back into $(13) to retrieve

\begin{align}
A(x, t) & = C_3\exp\Bigg( -\frac{V_x^2}{4D} t \Bigg)\exp\Bigg(\frac{V_x}{2D}x \Bigg) \ . & (19)
\end{align}

Assuming that $C_3 = 1$ we can now get the general transformation of $u$

\begin{align}
u(x, t) & = \exp\Bigg( -\frac{V_x^2}{4D} t \Bigg)\exp\Bigg(\frac{V_x}{2D}x \Bigg) v(x,t) \ . & (20)
\end{align}

Substituting $(20)$ back into $(1)$, opening the brackets, simplifying and rearranging yields indeed a simpler form the PDE as

\begin{align}
\frac{\partial v}{\partial t} & = D \frac{\partial^2 v}{\partial x^2} \ . & (21)
\end{align}

Note that must adapt the initial condition according to the transformation we just did as

\begin{align}
u(x, 0) & = A(x, 0) v(x, 0) = \exp\Bigg( \frac{V_x}{2D} x\Bigg) v(x, 0) = f(x) \\
\Leftrightarrow v(x, 0) & = \exp\Bigg( -\frac{V_x}{2D} x\Bigg)f(x) \ . & (22)
\end{align}

Before considering the different types of boundary conditions, we split the partial differential equation (PDE) (22) into a set of two ODEs using the separation of variable method

\begin{align}
v(x, t) & = X(x) T(t) \ . & (23)
\end{align}

In $(23)$, $X$ is only a function of $x$ and $T$ is only a function of $t$. Substituting $(23)$ into $(21)$, opening the brackets using the product rule and rearranging yields

\begin{align}
\frac{\partial (XT)}{\partial t} & = D \frac{\partial^2(XT)}{\partial x^2} \\
\Leftrightarrow X\frac{\partial T}{\partial t} & = D T \frac{\partial^2 X}{\partial x^2} \\
\Leftrightarrow \frac{1}{X} \frac{\partial^2 X}{\partial x^2} & = \frac{1}{D T}\frac{\partial T}{\partial t} = \lambda \ . & (24)
\end{align}

Since equation $(24)$ must hold at any instance in space and time up to the constant $\lambda$, we can now formulate the system of ODE

\begin{cases}
\frac{d^2 X}{d x^2} = \lambda X \\
\frac{d T}{dt} = D \lambda T \ . & (25)
\end{cases}

## Solution for Dirichlet BC
First, we should adapt the Dirichlet BC to the transformed function

\begin{align}
u(x, t) & = A(x, t) v(x,t) = \exp\Bigg( -\frac{V_x^2}{4D} t \Bigg)\exp\Bigg(\frac{V_x}{2D}x \Bigg) v(x,t) & (26\\
u(0, t) & = \exp\Bigg( -\frac{V_x^2}{4D} t \Bigg)\exp\Bigg(\frac{V_x}{2D} 0 \Bigg) v(0,t) = \exp\Bigg( -\frac{V_x^2}{4D} t \Bigg) v(0, t) = 0 & (27) \\
u(L, t) & = \exp\Bigg( -\frac{V_x^2}{4D} t \Bigg)\exp\Bigg(\frac{V_x}{2D} L \Bigg) v(L,t) = 0 \ . & (28)
\end{align}

For $(27)$ & $(28)$ to hold, it follows that 

\begin{align}
v(0,t) = v(L,t) = 0 \ .
\end{align}

Second, we consider the two ODEs derived above subsequently, starting with the spatial in $(25)$ and its boundary conditions

\begin{align}
X(0) & = X(L) = 0 \ . & (29)
\end{align}

The following three cases need to be considered: $\lambda = 0, \ \lambda > 0, \ \lambda < 0$. Let's which of the cases results in a non-trivial solution.

**Case 1, $\lambda = 0$**
Here, we consider

\begin{align}
\frac{d^2 X}{d x^2} & = 0 X = 0 \ & (30)
\end{align}

and integrating twice yields the general solution

\begin{align}
X(x) & = ax + b \ . & (31) 
\end{align}

Using the boundary conditions we can determine the values of $a = b = 0$ which is indeed the trivial solution to the problem.

**Case 2, $\lambda = k^2 > 0$**
Assuming any solution of type

\begin{align}
X(x) = \exp(rx) \ & (32)
\end{align}

and plugging this into $(25)$ yields

\begin{align}
& \frac{d^2 (\exp(rx))}{d x2} & = & \lambda \exp(rx) & \\
\Leftrightarrow \  & r^2 \exp(rx) & = & \lambda \exp(rx) & \\
\Leftrightarrow \ &  r^2 - \lambda & = & 0 & \\  
\Leftrightarrow \ &  r^2 - k^2    & = & 0 & (33)
\end{align}

for which we can determine values for

\begin{align}
r & = \pm k \ . & (34)
\end{align}

This will give the following general solution

\begin{align}
X(x) & = a \exp(kx) + b \exp(-kx) & (35)
\end{align}

and using our Dirichlet BC 

\begin{align}
X(0) & = a + b & = 0 & \rightarrow a = -b \\
X(L) & = a (\exp(kL) - \exp(-kL)) & = 0 & \rightarrow a = 0 & (36)
\end{align}

because $(\exp(kx) - \exp(-kx)) > 0$. As we can see, this case also yields a trivial solution.

**Case 3 $\lambda = -k^2 < 0$**
Following the approach of **Case 2** we can find

\begin{align}
r & = \pm ik & (37)
\end{align}

yielding the general solution

\begin{align}
& X(x) & = \ & a \exp(ikx) + b \exp(-ikx) \\
&      & = \ & a (\cos(kx) + i\sin(kx)) + b (cos(kx) - i\sin(kx)) \\
&      & = \ & (a + b)\cos(kx) + i(a-b)\sin(kx) \\
& X(x) & = \ & e\cos(kx) + f\sin(kx) \ . & (38)  
\end{align}

We now have to find $e$ and $f$ to satisfy the Dirichlet BC

\begin{align}
X(0) & =  e\cos(0) + f\sin(0) =  0 & & \rightarrow e = 0 \\
X(L) & = f\sin(kL) = 0 \ . & (39)
\end{align}

Since we are seeking for nontrivial solutions, we enforce $f \ne 0$ and $\sin(kL)=0$ which is true for 

\begin{align}
k & = \frac{n\pi}{L} \ \Bigg| \  n = 1,2,3 ... & \ (40)
\end{align}

yielding the following general solution to $(38)$

\begin{align}
X_n(x) & = f\sin\Bigg( \frac{n\pi}{L}x \Bigg) \ \Bigg| \ n = 1, 2, 3 ... & \ . (41) 
\end{align}

We can choose $f=1$ as arbitrary constant.

Let's consider the ODE in time now

\begin{align}
& & \frac{dT}{dt} & =  D \lambda T \\
& & \frac{1}{T} dT & =  D \lambda dt \ | \int \\
& &\ln(T) & =  D \lambda t + C_4 \\ 
& \Leftrightarrow & T(t) & =  C_5 \exp(D \lambda t) = C_5 \exp\Bigg( -D \Bigg( \frac{n\pi}{L} \Bigg)^2 t \Bigg) \ & (42)
\end{align}

and again $C_5 = 1$ can be arbitrarily choosen.

Finally we can combine the obtain solutions $(41)$ and $(42)$ to get

\begin{align}
X_n(x)T_n(t) & = \sin \Bigg( \frac{n\pi}{L}x \Bigg) \exp \Bigg( -D \Bigg( \frac{n\pi}{L} \Bigg)^2 t \Bigg) \ & \Bigg| \ n = 1,2,3, ... & (43)
\end{align}

This gives us the general solution to equation $(21)$ in the form of

\begin{align}
v(x,t) & = \sum^{\infty}_{n=1} \zeta_n \sin\Bigg(\frac{n\pi}{L} x \Bigg) \exp \Bigg( -D \Bigg( \frac{n\pi}{L} \Bigg)^2 t \Bigg) \ & (44)
\end{align}

where the coefficients $\zeta_n$ can be determined from the intial condition $(22)

\begin{align}
v(x, 0) & = \exp\Bigg( -\frac{V_x}{2D} x\Bigg)f(x) = \sum^{\infty}_{n=1} \zeta_n \sin\Bigg(\frac{n\pi}{L} x \Bigg) \ . & (45)
\end{align}

Multiplying $(45)$ by $\sin(\frac{m\pi}{L}x)$ and integrating yields

\begin{align}
\int^{L}_{0} \exp\Bigg( -\frac{V_x}{2D} x\Bigg)f(x) \sin \Bigg( \frac{m\pi}{L}x \Bigg) dx & = \int^{L}_{0} \sum^{\infty}_{n=1} \zeta_n \sin\Bigg(\frac{n\pi}{L} x \Bigg) \sin \Bigg( \frac{m\pi}{L}x \Bigg) dx \ . & (46)
\end{align}

Note that $\sin(n\pi x)$ is [orthogonal](https://tutorial.math.lamar.edu/classes/de/PeriodicOrthogonal.aspx) which allows us to simplify $(46)$ to

\begin{align}
\zeta_n = \frac{2}{L} \int^{L}_{0} \exp\Bigg( -\frac{V_x}{2D} x\Bigg)f(x) \sin \Bigg( \frac{n\pi}{L}x \Bigg) dx \ . & (47)
\end{align}

Finally, we have everything needed to obtain a general analytical solution to $(1)$ for Dirichlet BC and the initial condition as

\begin{align}
 u(x,t) & =  A(x,t)v(x,t) \\
        & =  \exp \Bigg(- \frac{V_x^2}{4D} t \Bigg) \exp \Bigg( \frac{V_x}{2D} x \Bigg) \sum^{\infty}_{n=1} \zeta_n \sin\Bigg(\frac{n\pi}{L} x \Bigg) \exp \Bigg( -D \Bigg( \frac{n\pi}{L} \Bigg)^2 t \Bigg) \ .& (48) 
\end{align}

## References
| No. | Reference |
| --- | --------- |
| [1] | Glinowiecka-Cox, M. B. (2022). Analytic Solution of 1D Diffusion-Convection Equation with Varying Boundary Conditions. |