# Noise Processes

A `NoiseProcess` is a type defined as

```julia
type NoiseProcess{class,inplace,F}
  noise_func::F
end
```

Its class is used for classifying the noise, and can be used by algorithms to
throw errors. For example, white noise uses `:White`. `inplace` denotes whether
the noise generating function is an inplace function. Lastly, we have the
`noise_func`. This is the function which is actually called in order to generate
the noise.

The signature for `noise_func` is

```julia
noise_func(rand_vec,integrator)
```

for inplace functions, and

```julia
rand_vec = noise_func(integrator)
```

otherwise. For not inplace noise functions where the equation is on an `AbstractArray`,
the signature

```julia
rand_vec = noise_func(x::Tuple,integrator)
```

where `x` is the size of the noise vector to make, is required. But it's highly
recommended for performance that one uses inplace noise updates with equations
on `AbstractArray`.

### White Noise

The default noise is `WHITE_NOISE`. This is the noise process which uses `randn!`.
A special dispatch is added for complex numbers for `(randn()+im*randn())/sqrt(2)`.
This function is `DiffEqBase.wiener_randn` (or with `!` respectively). Thus
its noise function is essentially:

```julia
noise_func(integrator) = randn()
noise_func(x::Tuple,integrator) = randn(x)
noise_func(rand_vec,integrator) = randn!(rand_vec)
```
