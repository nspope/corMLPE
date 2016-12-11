corMLPE
=======

This package implements a correlation structure for the R package `nlme` for Clarke's maximum likelihood population effects model (Clarke et al. 2002). This is useful, e.g. to construct regressions on distance matrices with nonlinearity and multiple random effects. Essentially, this is a model for symmetric, relation data.

The `corStruct` object allows a single grouping factor; for example among several species of an organism, we could fit an isolation-by-distance model with a syntax such as

```{r}
lme(genetic.distance ~ geographic.distance, random = ~geographic.distance|species, 
    correlation=corMLPE(form=~pop1+pop2|species), data=my.data)
```

In this case `pop1` and `pop2` are numerical labels for the populations, such that each observation corresponds to a pair of populations.

The package currently interfaces with `lme`, `nlme`, `gls`, and `gamm` (from package `mgcv`). The motivation for implementing this model in `nlme` was to leverage pre-existing machinery for nonlinear models, heteroskedasticity, and additive models. Look at Bill Peterman's package `ResistanceGA`[https://github.com/wpeterman/ResistanceGA] for a nice, `lme4`-based implementation.

Note that unlike the results presented in Clarke's paper, `nlme`/`lme` will return GLS standard errors rather than the `OLS` standard errors. If OLS standard errors are desired for some reason, see function `MLPE()` in the package.

# References
Clarke et al. 2002. Confidence Limits for Regression Relationships between Distance Matrices: Estimating Gene Flow with Distance. Journal of Agricultural, Biological, and Environmental Statistics 7: 361-372.
