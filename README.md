# Scanning-Tunneling-Microscope-Uncertainty-Calculation
This project was undertaken in the summer of 2022. It consists of investigating the minimum possible uncertainty of a scanning tunnelling microscope. This was done by first solving the time-independent Schr&ouml;dinger equation using a linearly biased potential in Mathematica. Then, an expression that can be used in the Newton-Raphson method was obtained by integrating the probability density analytically and equating it to a theoretical transmittance. Next, a Monte Carlo method was used to sample the probability per unit length of an electron tunnelling through the potential barrier; this was performed in C++. Python was used to analyse and visualise the uncertainties.

## Theory
The system is formed by an incoming flux of particles from negative infinity which encounters a potential barrier of initial height $\phi$ at $x=0$ which linearly decreases to a height of $\phi-b$ at $x=w$. This gives a potential function of the system of
```math
V(x) = 
\begin{cases}
    0 & ;\, x<0 \\
    \phi - x\frac{b}{w} & ;\, x \in [0,w]\\
    -b & ;\, x>w
\end{cases}
```
The left and right of the barrier have the known solutions
```math
\psi_{\textrm{I}}(x) = e^{ik_1x} + Ae^{-ik_1x} \,\, ;\, k_1 = \frac{\sqrt{2mE}}{\hbar}
```
for the incoming flux and
```math
\psi_{\textrm{III}}(x) = De^{ik_3x} \,\, ;\, k_3 = \frac{\sqrt{2m(E+b)}}{\hbar}
```
for the outgoing flux, where $`A`$ and $`D`$ are complex constants. The second-order linear ordinary differential equation to solve for region II takes the form
```math
\frac{d^2\psi_{\textrm{II}}}{dx^2} + \frac{2m}{\hbar}\psi_{\textrm{II}}(E-\phi + x\frac{b}{w})=0 \,\, ,
```
which has the solution
```math
\begin{split}
    \psi_{\textrm{II}}(x) =& c_1 \textrm{Ai}(\alpha(x)) + c_2 \textrm{Bi}(\alpha(x))  \\
                    &;\, \alpha(x) = -\frac{1}{b}\Big(\frac{2bm}{w\hbar^2}\Big)^{\frac{1}{3}} ((E-\phi)w + bx)
\end{split}
\,\, ,
```
where $`c_1`$ and $`c_2`$ are complex coefficients and $`\textrm{Ai}(\alpha(x))`$ and $`\textrm{Bi}(\alpha(x))`$ are the Airy functions. To maintain the continuity of the system it must satisfy the relation
```math
\begin{pmatrix}
    1+A \\
    De^{ik_3w}
\end{pmatrix}
=\underline{M}
\begin{pmatrix}
    c_1 \\
    c_2
\end{pmatrix}
\,\, ,
```
where
```math
\underline{M} =
 \begin{pmatrix}
    \textrm{Ai}(\alpha(0)) & \textrm{Bi}(\alpha(0)) \\
    \textrm{Ai}(\alpha(w)) & \textrm{Bi}(\alpha(w)) 
\end{pmatrix}
\,\, .
```
Then to maintain continuity in the derivatives the system must also satisfy
```math
\zeta\underline{M}'
\begin{pmatrix}
    c_1 \\
    c_2
\end{pmatrix}
=\begin{pmatrix}
    ik_1(1 - A) \\
    Dik_3e^{ik_3w}
\end{pmatrix}
\,\, ,
```
where
```math
\underline{M}' = 
\begin{pmatrix}
    \textrm{Ai}'(\alpha(0)) & \textrm{Bi}'(\alpha(0)) \\
    \textrm{Ai}'(\alpha(w)) & \textrm{Bi}'(\alpha(w)) 
\end{pmatrix}
```
and $`\zeta = \frac{\partial \alpha}{\partial x} = -\Big(\frac{2bm}{w\hbar^2}\Big)^{\frac{1}{3}}`$. These can be used to solve for the coefficients $`A`$ and $`D`$ through
```math
\zeta\underline{M}'\underline{M}^{-1}
\begin{pmatrix}
    1 + A \\
    De^{ik_1w}
\end{pmatrix}
=\begin{pmatrix}
    ik_1(1 - A) \\
    Dik_3e^{ik_3w}
\end{pmatrix}
\,\, .
```
By defining $`\underline{\kappa} := \zeta\underline{M}'\underline{M}^{-1}`$ this matrix equation can be solved to give
```math
A = -\frac{\kappa_{12}\kappa_{21} + ik_1\kappa_{22} - \kappa_{11}\kappa_{22} + k_1k_3 + i\kappa_{11}k_3}{\kappa_{12}\kappa_{21} - ik_1\kappa_{22} - \kappa_{11}\kappa_{22} - k_1k_3 + i\kappa_{11}k_3}
```
and
```math
D = -\frac{2k_1\kappa_{21}e^{-ik_3w}}{i\kappa_{12}\kappa_{21} + k_1\kappa_{22} - i\kappa_{11}\kappa_{22} - ik_1k_3 - \kappa_{11}k_3} \,\, .
```
The coefficients in $`\psi_{II}(x)`$ can then be determined using
```math
\begin{split}
    \begin{pmatrix}
        c_1 \\
        c_2
    \end{pmatrix}
    =&
    \begin{pmatrix}
        (1+A)\textrm{Bi}(\alpha(w)) - \textrm{Bi}(\alpha(0))De^{ik_1w} \\
        \textrm{Ai}(\alpha(0))De^{ik_1w}-(1+A)\textrm{Ai}(\alpha(w))
    \end{pmatrix} \\
    &\cdot\frac{1}{\textrm{Ai}(\alpha(0))\textrm{Bi}(\alpha(w)) - \textrm{Ai}(\alpha(w))\textrm{Bi}(\alpha(0))}
\end{split}
\,\, .
```


The proportion of electrons that tunnel through is equivalent to the area under the probability density in region III over the total area under the whole probability density. Only fluxes should be considered so arbitrary lengths can be used in the three regions, although in region II we will use a length of $w$. This gives a tunnelling proportion/probability of
```math
r = \frac{L_3T(w)}{L_1\langle|\psi_I(x)|^2\rangle + w\langle|\psi_{II}(x)|^2\rangle + L_3T(w)} \,\, ,
```
where the average function values apply over the whole space that part of the wavefunction acts and $`T(w)`$ is the transmittance which is equivalent to $|D|^2$ from the equation defined previously. The first average function value is calculated from the integral
```math
\langle|\psi_I(x)|^2\rangle = \frac{1}{L_1}\int\limits_{-L_1}^0 |e^{ik_1x} + Ae^{-ik_1x}| dx \, \, ,
```
which evaluates to
```math
\begin{split}
    \langle|\psi_I(x)|^2\rangle = \Big(&1 + |A|^2 \\
                        &+ \frac{1}{2ik_1L_1}\big(A^*(1 - e^{-2ik_1L_1}) \\
                        &+ A(e^{2ik_1L_1} - 1)\big)\Big)
\end{split}
\,\, .
```
It can be seen that if $`L_1\rightarrow\infty`$ or when $`L_1=\frac{n\pi}{k_1} \, \, ;\, n\in\mathbb{N}`$ this expression simplifies to
```math
\langle|\psi_I(x)|^2\rangle = 1 + |A|^2
```
so it would be best to use any $`L_1=\frac{n\pi}{k_1} \,\, ;\, n\in\mathbb{Z}^+`$ (for computational ease). The second average function is more difficultly calculated by the integral
```math
I(w) = \langle|\psi_{II}(x)|^2\rangle = \int\limits_0^w |c_1\textrm{Ai}(\alpha(x)) + c_2\textrm{Bi}(\alpha(x))|^2 dx \,\, ,
```
which was calculated using Mathematica by noticing that the Airy functions must return real values since they have real inputs.


By choosing $`L_1=\frac{\pi}{k_1}`$ and $`L_3=(\gamma-1) w`$ where $`1<\gamma\in\mathbb{R}`$ we find that
```math
T(w)(r(\gamma-1) w - \frac{\pi}{k_1}r - (\gamma-1) w) + \frac{2\pi r}{k_1} + rI(w)=0 \,\, ,
```
which can be expressed as the function
```math
f(w) = T(w)(\frac{\pi}{k_1} +( \gamma-1) w(\frac{1}{r} -1)) - \frac{2\pi}{k_1} - I(w)=0 \,\, ,
```
where a root-finding method like Newton-Raphson can be used to find $`w`$.
An appropriate value of $\gamma$ can be selected such that the tunnelling proportion is defined to be $`0.5`$ when the width is at its average value. This can be shown to be
```math
\gamma = 1+\frac{\frac{\pi}{k_1\overline{w}}(2 - T(\overline{w})) + I(\overline{w})}{T(\overline{w})} \,\, .
```
