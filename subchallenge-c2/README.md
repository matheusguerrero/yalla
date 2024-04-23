# Subchallenge C2

Note that, whilst the training data are not provided (because they are too large), we do provide the trained estimator (used in the paper) in `/NBE/intermediates/` alongside the bootstrap quantile estimates.

* `c2-1-simulate.R` - This script simulates the training data. Note that it requires that you first run `c1-1-fits.R` with quantile level `q0 <- 0.6` and number of bootstrap estimates `n.boots <- 750`. The full training data are quite big and the computations can be quite time consuming; with `K = 350000` sets of training data, expect an Rdata file of approximately 60GB.
* `c2-2-train.jl` - This is a Julia script for training the neural Bayes estimator and evaluating it on the observational data (as well as producing bootstrap estimates). Training is done using the [NeuralEstimators.jl](https://github.com/msainsburydale/NeuralEstimators.jl) package; please set the linked repo for details on installation. This code will not run on a standard laptop and requires a high-end GPU to run in a reasonable timeframe.
 * `c2-3-plot.R` - Plots the results and prints out bootstrap estimates. You can run this script without running the two proceeding scripts.
