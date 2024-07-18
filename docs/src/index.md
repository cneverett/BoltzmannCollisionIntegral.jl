# BinaryInteractionSpectra.jl

`BinaryInteractionSpectra.jl` is a package for the evaluation of the collision integral from binary interactions of arbitrary momentum and anisotropic particle distributions using a method of discretisation and Monte-Carlo integration. 

The evaluation of the collision integral is of great use when studying kinetic system. The evolution of particle distributions ``f(x^a,\boldsymbol{p})`` within such a system is dictated by the collisional relativistic Boltzmann equation (neglecting external forcing):
```math
p_1^\nu\partial_\nu f(x^\mu,\boldsymbol{p}_1)=C(x^\mu,\boldsymbol{p}_1),
```
where ``p^\mu=(p^0,\boldsymbol{p})`` is the 4-Momentum and ``C(x^a,\boldsymbol{p})`` is the collision integral. The collision integral is given by 
```math
\begin{multline}
    C(\boldsymbol{p}_1)=\int\frac{\mathrm{d}^3\boldsymbol{p}_2}{p_2^0}\frac{\mathrm{d}^3\boldsymbol{p}_3}{p_3^0}\frac{\mathrm{d}^3\boldsymbol{p}_4}{p_4^0}\left[\frac{f(\boldsymbol{p}_3)f(\boldsymbol{p}_4)}{1+\updelta_{34}}W(p_3^\mu,p_4^\mu|p_1^\mu,p_2^\mu)\right. \\ \left.- \frac{f(\boldsymbol{p}_1)f(\boldsymbol{p}_2)}{1+\updelta_{34}}W(p_1^\mu,p_2^\mu|p_3^\mu,p_4^\mu)\right],
    \end{multline}
```
with the transition rate ``W`` being given by:
```math
\begin{equation}
W(p_1^\mu,p_2^\mu|p_3^\mu,p_4^\mu)=\frac{sp_{\text{in}}^{*2}}{\pi}\sigma_{12|34}(s,t)\updelta^{(4)}(p_1^\mu+p_2^\mu-p_3^\mu-p_4^\mu). 
\end{equation}
```

The Boltzmann equation can then be discretised in momentum space by integrating over a discrete volume element ``\Delta^3\boldsymbol{p}_1`` to yield: 
```math 
\int_{\Delta^3\boldsymbol{p}_1}\frac{\mathrm{d}^3\boldsymbol{p}_1}{p_1^0}p_1^\mu\partial_\mu f(x^\mu,\boldsymbol{p}_1)=\int_{\Delta^3\boldsymbol{p}_1}\mathrm{d}^3\boldsymbol{p}_1~S_{spe}(x^\mu,\boldsymbol{p}_1)-T_{spe}(x^\mu,\boldsymbol{p}_1),
```
where the collision integral has been split into two parts, an emission spectrum
$$
\begin{multline}
    S_{\text{spe}}(\boldsymbol{p}_1)=\frac{1}{1+\updelta_{34}}\frac{1}{p^0_1}\int\frac{\mathrm{d}^3\boldsymbol{p}_2}{p_2^0}\frac{\mathrm{d}^3\boldsymbol{p}_3}{p_3^0}\frac{\mathrm{d}^3\boldsymbol{p}_4}{p_4^0} \\ \times f(\boldsymbol{p}_3)f(\boldsymbol{p}_4)W(p_3^\mu,p_4^\mu|p_1^\mu,p_2^\mu),
\end{multline}
$$
dictating the rate of gain of particles of a specific type at a specific momentum from a binary interaction, and an absorption spectrum 
```math
\begin{multline}
T_{\text{spe}}(\boldsymbol{p}_1)=\frac{1}{1+\updelta_{34}}\frac{1}{p^0_1}\int\frac{\mathrm{d}^3\boldsymbol{p}_2}{p_2^0}\frac{\mathrm{d}^3\boldsymbol{p}_3}{p_3^0}\frac{\mathrm{d}^3\boldsymbol{p}_4}{p_4^0} \\ \times f(\boldsymbol{p}_1)f(\boldsymbol{p}_2)W(p_1^\mu,p_2^\mu|p_3^\mu,p_4^\mu),
\end{multline}
```
dictating the rate of loss of particles of a specific type at a specific momentum from the same interaction.

## Assumptions
In order to enable the evaluation of the emission and absorption terms, we make the following set of assumptions:  
- The system is axisymmetric in momentum space coordinates i.e. only dependant on ``p=|\boldsymbol{p}|`` and the cosine of the angle to the axis of symmetry ``\mu`` where ``\hat{z}\cdot\boldsymbol{p}=p\cos\theta=p\mu`` 
- The distribution functions are redefined in an axisymmetric form i.e. ``f(\boldsymbol{p})\equiv\frac{f(p,\mu)}{2\pi p^2}``.
	- This is done such that a single particle is described by the distribution ``f(p,\mu)=\delta(p-p_i)\delta(\mu-\mu_i)``.
- The axisymmetric distribution functions are averaged over phase space intervals ``\Delta p\Delta\mu``, to generate discrete values i.e. ``f(t,p_i,\mu_j)=\frac{1}{\Delta p_i\Delta\mu_j\Delta z_\alpha}\int_{\Delta p_i\Delta\mu_j\Delta z_\alpha} \mathrm{d}p \mathrm{d}\mu \mathrm{d}z~f(z,t,p,\mu)``, such that they act as if they are constant over that interval of phase space.  