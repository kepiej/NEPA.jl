# NEPA
[![Build Status](https://travis-ci.org/kepiej/NEPA.jl.svg?branch=master)](https://travis-ci.org/kepiej/NEPA.jl) [![codecov](https://codecov.io/gh/kepiej/NEPA.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/kepiej/NEPA.jl)

Nonparametric Efficiency and Productivity Analysis for Julia

A toolbox for doing efficiency analysis and computing productivity change using nonparametric techniques. Efficiency can be computed using
* Data Envelopment Analysis (DEA)
* Free Disposal Hulls (FDH)
* Luenberger's directional distance functions
* Slacks-based measure (SBM) of Tone (2001)

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

Analogeously, FDH models are implemented under the various returns to scale assumptions: the function names are FDH_VRS, FDH_CRS, FDH_NIRS and FDH_NDRS.

The directional distance function (DDF) for a given direction (gx,gy) under a convex VRS technology can be computed using:

```julia
  # Choose direction vectors. Here, we set the direction vectors equal to the observations.
  gx = X
  gy = Y
  # Initialize a convex directional distance function under VRS
  D = DDF{Convex,VRS}(X,Y,gx,gy)
  # Solve the model for all K observations. beta is a vector containing K efficiency scores.
  beta = D()
```

The {Convex,VRS} in DDF{Convex,VRS} specifies that the DDF is a convex VRS technology. Of course, any other return to scale assumption (CRS, NDRS and NIRS) can be specified in the same way. Finally, non-convex technologies are specified using {FreeDisposal,VRS}. Only VRS is currently implemented for the non-convex DDF.

## TODO

* Add Malmquist index;
* Luenberger-Hicks-Moorsteen TFP;
* Add more outlier detection methods (order-m is already available);
* ...
