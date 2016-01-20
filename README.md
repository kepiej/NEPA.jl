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

  # Initialize an input-oriented DEA model with variable returns to scale (VRS)
  input = true # Set input = false for output-oriented DEA model
  D = DEA_VRS(X,Y,input)
  # Solve the model for all K observations. theta is a vector containing K efficiency scores.
  theta = D()
  # DEA models with other returns to scale can be initialized: DEA_CRS (CRS), DEA_NIRS (IRS), DEA_NDRS (NDRS)
```

## TODO

* Add Malmquist index;
* Luenberger-Hicks-Moorsteen TFP;
* Add outlier detection;
* ...
