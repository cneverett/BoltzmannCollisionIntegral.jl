# BinaryInteractionSpectra.jl

`BinaryInteractionSpectra.jl` is a package for the evaluation of the collision integral from binary interactions of arbitrary momentum and anisotropic particle distributions using a method of discretisation and Monte-Carlo integration. 

The evaluation of the collision integral is of great use when studying kinetic system. The evolution of particle distributions ``f(x^a,\boldsymbol{p})`` within such a system is dictated by the collisional relativistic Boltzmann equation (neglecting external forcing):
```math
p_1^\nu\partial_\nu f(x^\mu,\boldsymbol{p}_1)=C(x^\mu,\boldsymbol{p}_1),
```
where ``p^\mu=(p^0,\boldsymbol{p})`` is the 4-Momentum and ``C(x^a,\boldsymbol{p})`` is the collision integral. The collision integral is given by 
```math
    C(\boldsymbol{p}_1)=\int\frac{\mathrm{d}^3\boldsymbol{p}_2}{p_2^0}\frac{\mathrm{d}^3\boldsymbol{p}_3}{p_3^0}\frac{\mathrm{d}^3\boldsymbol{p}_4}{p_4^0}\left[\frac{f(\boldsymbol{p}_3)f(\boldsymbol{p}_4)}{1+\updelta_{34}}W(p_3^\mu,p_4^\mu\|p_1^\mu,p_2^\mu)- \frac{f(\boldsymbol{p}_1)f(\boldsymbol{p}_2)}{1+\updelta_{34}}W(p_1^\mu,p_2^\mu\|p_3^\mu,p_4^\mu)\right],
```
with the transition rate ``W`` being given by:
```math
W(p_1^\mu,p_2^\mu|p_3^\mu,p_4^\mu)=\frac{sp_{\text{in}}^{*2}}{\pi}\sigma_{12|34}(s,t)\updelta^{(4)}(p_1^\mu+p_2^\mu-p_3^\mu-p_4^\mu). 
```

The Boltzmann equation can then be discretised in momentum space by integrating over a discrete volume element ``\Delta^3\boldsymbol{p}_1`` to yield: 
```math 
\int_{\Delta^3\boldsymbol{p}_1}\frac{\mathrm{d}^3\boldsymbol{p}_1}{p_1^0}p_1^\mu\partial_\mu f(x^\mu,\boldsymbol{p}_1)=\int_{\Delta^3\boldsymbol{p}_1}\mathrm{d}^3\boldsymbol{p}_1~S_{spe}(x^\mu,\boldsymbol{p}_1)-T_{spe}(x^\mu,\boldsymbol{p}_1),
```
where the collision integral has been split into two parts, an emission spectrum
```math
S_{\text{spe}}(\boldsymbol{p}_1)=\frac{1}{1+\updelta_{34}}\frac{1}{p^0_1}\int\frac{\mathrm{d}^3\boldsymbol{p}_2}{p_2^0}\frac{\mathrm{d}^3\boldsymbol{p}_3}{p_3^0}\frac{\mathrm{d}^3\boldsymbol{p}_4}{p_4^0}f(\boldsymbol{p}_3)f(\boldsymbol{p}_4)W(p_3^\mu,p_4^\mu|p_1^\mu,p_2^\mu),
```math
dictating the rate of gain of particles of a specific type at a specific momentum from a binary interaction, and an absorption spectrum 
```math
T_{\text{spe}}(\boldsymbol{p}_1)=\frac{1}{1+\updelta_{34}}\frac{1}{p^0_1}\int\frac{\mathrm{d}^3\boldsymbol{p}_2}{p_2^0}\frac{\mathrm{d}^3\boldsymbol{p}_3}{p_3^0}\frac{\mathrm{d}^3\boldsymbol{p}_4}{p_4^0} f(\boldsymbol{p}_1)f(\boldsymbol{p}_2)W(p_1^\mu,p_2^\mu|p_3^\mu,p_4^\mu),
```
dictating the rate of loss of particles of a specific type at a specific momentum from the same interaction.

## Assumptions
In order to enable the evaluation of the emission and absorption terms, we make the following set of assumptions (dropping dependence on space-time for notational simplicity):  
- The system is axisymmetric in momentum space coordinates i.e. only dependant on ``p=|\boldsymbol{p}|`` and the cosine of the angle to the axis of symmetry (``\hat{z}``) ``\mu`` where ``\hat{z}\cdot\boldsymbol{p}=p\cos\theta=p\mu`` 
- The distribution functions are redefined in an axisymmetric form i.e. ``f(\boldsymbol{p})\equiv\frac{f(p,\mu)}{2\pi p^2}``. This is done such that a single particle is described by the distribution ``f(p,\mu)=\delta(p-p_i)\delta(\mu-\mu_i)``.
- The axisymmetric distribution functions are averaged over phase space intervals ``\Delta p\Delta\mu``, to generate discrete values i.e. ``f(t,p_i,\mu_j)=\frac{1}{\Delta p_i\Delta\mu_j\int_{\Delta p_i\Delta\mu_j} \mathrm{d}p \mathrm{d}\mu~f(p,\mu)``, such that they act as if they are constant over that interval of phase space.  

## Discrete Form of the Boltzmann Equation
Following the assumptions laid out in the previous section. The axisymmetric, momentum-discrete Boltzmann equation is given by 
```math
\left(\partial_t+\braket{\beta}_i\braket{\mu}_j\partial_z\right)f_{1,ij}=f_{3,kl}f_{4,mn}S_{34|12,ijklmn}-f_{1,ij}f_{2,kl}T_{12|34,ijkl}
```
The discrete emission spectrum term is a 6D array, given by
```math
f_{3,kl}f_{4,mn}S_{34|12,ijklmn}=f(z,t,p_{3,k},\mu_{3,l})f(z,t,p_{4,m},\mu_{4,n})\frac{\Delta p_{3,k}\Delta \mu_{3,l}\Delta p_{4,m}\Delta \mu_{4,n}}{\Delta p_i}\frac{1}{N}\sum^{N}_{a=1}\left[S_{val}\right](\{\boldsymbol{p}_1,\boldsymbol{p}_3,\boldsymbol{p}_4\}_a), 
```
with 
```math
S_{val}=\frac{1+\delta_{12}}{1+\delta_{34}}\sum_\pm\frac{2p_\pm^2\mathcal{F}_{34}^2}{ p_3^0p_4^0}\frac{d\sigma_{34|12}}{dt}(s,t_\pm)\frac{1}{p_3^0p_\pm-p_\pm^0p_3\cos\Theta_{\pm3}+p_4^0p_\pm-p_\pm^0p_4\cos\Theta_{\pm4}},
```
while the discrete absorption term is a 4D Array, given by
```math
f_{1,ij}f_{2,kl}T_{12|34,ijkl}=f(z,t,p_i,\mu_j)f(z,t,p_{2,k},\mu_{2,l})\Delta\boldsymbol{p}_{2,k}\Delta\mu_{2,l}\frac{1}{N}\sum^{N}_{a=1}\left[T_{val}\right](\{\boldsymbol{p}_1,\boldsymbol{p}_2\}_a),
```
with
```math
T_{val}=\frac{1}{1+\delta_{12}}\frac{\mathcal{F}_{12}(s)\sigma_{12|34}(s)}{p^0p_2^0}.
```

## Evaluation
Given a user input of interaction to evaluate and what discretisation to use, the discrete emission spectrum ``S_{12|34,ijklmn}`` (note the change of order of the interaction) and absorption spectrum ``T_{12|34,ijkl}`` are solved simultaneously by Monte-Carlo integration. The evaluation is performed by the functions `SpectraEvaluateSerial()` or `SpectralEvaluateMultiThread()` depending on if the code is to run in serially on a single core or multi-threaded on multiple cores of a single CPU. The code will then save the resulting multidimensional arrays.

!!! warning The evaluation does not guarantee convergence of the Monte-Carlo integration, instead it will sample only a user defined number of points. A rough measure of convergence is calculated and stored as `SConverge` and `TConverge` in the output file.

A quick setup guide is provided in the [Getting Started](@ref getting_started)section.  