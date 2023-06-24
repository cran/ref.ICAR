# ref.ICAR R Package

Maintainer: Erica M. Porter <ericamp@vt.edu>

Implements an objective Bayes intrinsic conditional autoregressive prior. This model provides an objective Bayesian approach for modeling spatially correlated areal data using an intrinsic conditional autoregressive prior on a vector of spatial random effects.

## Version 2.0 changes

Version 2.0 of ref.ICAR includes several changes:

* Replaces the use of R package dependencies `rgdal` and `maptools` with functions from   `sf`, since `rgdal` and `maptools` will be archived soon     (https://rsbivand.github.io/csds_jan23/bivand_csds_ssg_230117.pdf).

* The function `ref.MCMC` implements faster sampling of parameters and random effects
  based on results from Ferreira et al. (2021).  See references below and in package     documentation.
  
* Now includes functions for objective Bayesian model selection based on the spatial
  ICAR model, as developed by Porter et al (2023).  See references below and in package
  documentation.
  
## References

Porter, E.M., Franck, C.T., and Ferreira, M.A.R. (2023), “Objective bayesian model selection for spatial hierarchical models with intrinsic conditional autoregressive priors,” Bayesian Analysis, International Society for Bayesian Analysis, 1, 1–27. https://doi.org/10.1214/23-BA1375.

Ferreira, M.A.R., Porter, E.M., and Franck, C.T. (2021), “Fast and scalable computations for Gaussian hierarchical models with intrinsic conditional autoregressive spatial random effects,” Computational Statistics and Data Analysis, 162. https://doi.org/10.1016/j.csda.2021.107264.

Keefe, M.J., Ferreira, M.A.R., and Franck, C.T. (2018), “On the formal specification of sum-zero constrained intrinsic conditional autoregressive models,” Spatial Statistics, Elsevier {BV}, 24, 54–65. https://doi.org/10.1016/j.spasta.2018.03.007.

Keefe, M.J., Ferreira, M.A.R., and Franck, C.T. (2019), “Objective Bayesian Analysis for Gaussian Hierarchical Models with Intrinsic Conditional Autoregressive Priors,” Bayesian Analysis, International Society for Bayesian Analysis, 14, 181–209. https://doi.org/10.1214/18-BA1107.

