# DiffEqBiological.jl API
```@meta
CurrentModule = DiffEqBiological
```

## Reaction Network Generation Macros
```@docs
@min_reaction_network
@reaction_network
```

## `@min_reaction_network` modifiers
```@docs
addodes!
addsdes!
addjumps!
```

## Basic properties
```@docs
speciesmap
paramsmap 
numspecies 
numparams
numreactions
```

## Reaction Properties
```@docs
substrates
products
dependents
dependants
ismassaction
substratestoich
productstoich
netstoich
```

## Generated Functions for Models
```@docs
oderhsfun
jacfun
paramjacfun
odefun
noisefun
sdefun
jumps
regularjumps
```

## Generated Expressions
```@docs
odeexprs
jacobianexprs
noiseexprs
jumpexprs
rateexpr
oderatelawexpr
ssaratelawexpr
```

## Dependency Graphs
```@docs
rxtospecies_depgraph
speciestorx_depgraph
rxtorx_depgraph
```