# BinaryInteractionSpectra.jl

`BinaryInteractionSpectra.jl` is a package for the evaluation of the collision integral from binary interactions of arbitrary momentum and anisotropic particle distributions using a method of discretisation and Monte-Carlo integration. 

The evaluation of the collision integral is of great use when studying kinetic system. The evolution of particle distributions ``f(x^a,\vec{p})`` within such a system is dictated by the collisional relativistic Boltzmann equation (neglecting external forcing):
```math
p_1^\nu\partial_\nu f(x^\mu,\boldsymbol{p}_1)=C(x^\mu,\boldsymbol{p}_1),
```
where ``p^\mu=(p^0,\boldsymbol{p})`` is the 4-Momentum and ``C(x^a,\vec{p})`` is the collision integral. The collision integral is given by 
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

The Boltzmann equation can then be discretised in momentum space by integrating over a discrete volume element ``\Delta^3\boldsymbol{p}_1``

