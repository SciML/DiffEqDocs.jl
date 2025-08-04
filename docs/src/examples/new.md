# Parcel Model with Explicit Condensation and Supersaturation

This example demonstrates a numerical implementation of a classic cloud physics model:  
**"An Elementary Parcel Model with Explicit Condensation and Supersaturation"**  
by R.R. Rogers (1975), using the Julia language and the [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/) library.

The purpose of this example is to reproduce figures from the original paper by simulating the temporal evolution of supersaturation, droplet radius, temperature, and liquid water content in an ascending air parcel.

---

## Physical Background

```@example new
print(a)
```

Cloud formation is a complex thermodynamic process. To analyze it, we use simplified approaches, such as the  **parcel model**, in which a small "parcel" of air ascends adiabatically through the atmosphere.

In this simplified model:

- An **air parcel** rises at a constant velocity `U`.
- As pressure decreases, **supersaturation** is created.
- Water vapor condenses on pre-existing droplets (resulting in their growth), releasing latent heat.
- This process changes the **droplet radius**, **temperature**, and **liquid water content** over time.

Key simplifying assumptions:
- No new droplet activation (fixed number of droplets).
- No coalescence or sedimentation.
- No mixing with the environment (idealized adiabatic ascent).


---

##  Model Setup

Initial conditions for the simulation (based on the paper):

- Initial droplet radius: $$r_0 = 8 \mu m$$
- Droplet concentration: $$n_0 = 200 cm^{-3}$$
- Updraft speed: $$U = 10 \frac{m}{s}$$
- Temperature: $$T_0 = 280.15 K (7^\circ C)$$
- Pressure: $$p_0 = 800 hPa$$
- Initial supersaturation: $$S_0 = 1.00$$ (saturated)
- Gravitational acceleration $$g = 9.81 \frac{m}{s^2}$$

These parameters are typical of rapidly developing cumulonimbus clouds in the lower troposphere.

The simulation integrates the following variables over time, which are described by following equations:
- Pressure 'p(t)' 
<div align="center">

$$
\frac{p}{t}  = -\rho_0 g U 
$$

</div>
where $$\rho = \frac{p}{R\cdot T}$$

- Droplet radius 'r(t)'
<div align="center">

 $$
  \frac{dr}{dt} = \frac{\sigma(S, T, p)}{r}
  $$

</div>
where $$\sigma = \frac{S-1}{F_D + F_K}$$


- Liquid water mixing ratio 'x(t)'
- Pressure 'p(t)'
<div align="center">

$$
\frac{S}{t}  = -\rho_0 g U 
$$

</div>



---

## Results

### Figure 1 — Droplet Growth and Supersaturation

![Droplet radius and supersaturation](/docs/build/assets/rogers_fig1.svg)



---

### Figure 2 — Temperature and Liquid Water Content

![Temperature and liquid water mixing ratio](/docs/build/assets/rogers_fig1.svg)


---

## Reference

> R.R. Rogers (1975), *An Elementary Parcel Model with Explicit Condensation and Supersaturation*,  
> Atmosphere, Vol. 13, No. 4, pp. 192–204.  
> DOI: [10.1080/00046973.1975.9648397](https://doi.org/10.1080/00046973.1975.9648397)


