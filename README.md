# NEPA
[![Build Status](https://travis-ci.org/kepiej/NEPA.jl.svg?branch=master)](https://travis-ci.org/kepiej/NEPA.jl)

Nonparametric Efficiency and Productivity Analysis for Julia

A toolbox for doing efficiency analysis and computing productivity change using nonparametric techniques. Efficiency can be computed using
* Data Envelopment Analysis (DEA)
* Free Disposal Hulls (FDH)
* Luenberger's directional distance functions

for the various returns to scale assumptions (VRS, CRS, NIRS and NDRS).

Productivity change can be calculated using
* Hicks-Moorsteen index
* Luenberger productivity indicator

and decomposed in the various components (technical change, technical (in)efficiency change, scale (in)efficiency change) to get insight in the drivers of productivity growth.

## Usage
```julia
  using NEPA

  # Read in input data in matrix X and output data in matrix Y
  # with dim(X) = [K,N] and dim(Y) = [K,M] where K is the number of observations, N the number of inputs and M the number of outputs

  # Set input = false for output-oriented DEA model
  input = true
  # Initialize an input-oriented DEA model with variable returns to scale (VRS)
  D = DEA_VRS(X,Y,input)
  # Solve the model for all K observations. theta is a vector containing K efficiency scores.
  theta = D()
  # DEA models with other returns to scale can also be initialized:
  #  - constant return to scale (CRS): DEA_CRS,
  #  - nonincreasing returns to scale (NIRS): DEA_NIRS,
  #  - nondecreasing returns to scale (NDRS): DEA_NDRS
```

Analogeously, FDH models are implemented under the various returns to scale assumptions. The efficiency scores can be computed using the function call:

```julia
  # Solve FDH model under VRS
  theta = FDH_VRS(X,Y,input)
  # Solve FDH model under CRS
  theta = FDH_CRS(X,Y,input)
  # Solve FDH model under NIRS
  theta = FDH_NIRS(X,Y,input)
  # Solve FDH model under NDRS
  theta = FDH_NDRS(X,Y,input)
```

## TODO

* Add Malmquist index;
* Luenberger-Hicks-Moorsteen TFP;
* Add outlier detection;
* ...
