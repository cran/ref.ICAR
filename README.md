# ref.ICAR R Package

Maintainer: Erica M. Porter <emporte@clemson.edu>

Implements an objective Bayes intrinsic conditional autoregressive prior. This model provides an objective Bayesian approach for modeling spatially correlated areal data using an intrinsic conditional autoregressive prior on a vector of spatial random effects.

## Version 2.0.1 changes

Version 2.0.2 of ref.ICAR includes minor changes and bug fixes:

* Changed how the spatial data for Example 2 of the vignette is read in, made necessary by changes to spData package.

* Changed how the plot for Example 2 of the vignette is generated, due to a conflict with ggplot2 and the updated sf package.

* Updated the references with correct years and volume numbers.
  
## References

Porter, E.M., Franck, C.T., and Ferreira, M.A.R. (2024), “Objective Bayesian model selection for spatial hierarchical models with intrinsic conditional autoregressive priors,” Bayesian Analysis, International Society for Bayesian Analysis, 19(4), 985-1011. https://doi.org/10.1214/23-BA1375.

Ferreira, M.A.R., Porter, E.M., and Franck, C.T. (2021), “Fast and scalable computations for Gaussian hierarchical models with intrinsic conditional autoregressive spatial random effects,” Computational Statistics and Data Analysis, 162, 107264. https://doi.org/10.1016/j.csda.2021.107264.

Keefe, M.J., Ferreira, M.A.R., and Franck, C.T. (2018), “On the formal specification of sum-zero constrained intrinsic conditional autoregressive models,” Spatial Statistics, Elsevier {BV}, 24, 54–65. https://doi.org/10.1016/j.spasta.2018.03.007.

Keefe, M.J., Ferreira, M.A.R., and Franck, C.T. (2019), “Objective Bayesian analysis for Gaussian hierarchical models with intrinsic conditional autoregressive priors,” Bayesian Analysis, International Society for Bayesian Analysis, 14, 181–209. https://doi.org/10.1214/18-BA1107.

