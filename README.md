# SENLM

An R package for modelling ecological data that is sampled across environmental gradients.

For installation and usage details, please visit the package website:

  https://primer-e.github.io/senlm/
  
## Other non-linear modeling packages

  - [gnlm: Generalized Nonlinear Regression Models](https://cran.r-project.org/web/packages/gnlm/index.html)
  - [gamlss.nl: Fitting non linear parametric GAMLSS models](https://cran.r-project.org/web/packages/gamlss.nl/index.html)
  - [Stan](https://mc-stan.org/)

## Contributing

After making changes to the packages run `devtools::build_site()`
in order to update the package's website (`docs/`) and 
`devtools::build_vignette()` in order to rebuild the vignettes (`doc/`).

You should read the following excellent books about developing R packages
(both of which happen to be written by Hadley Wickham):

  - [R packages](http://r-pkgs.had.co.nz/)
  - [Advanced R](https://adv-r.hadley.nz/)
  
If you want to update the package's website, then read the [pkgdown documentation](https://pkgdown.r-lib.org/).
