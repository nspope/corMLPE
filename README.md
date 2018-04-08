corMLPE
=======

This package implements a correlation structure for the R package `nlme` for Clarke's maximum likelihood population effects model (Clarke et al. 2002). This is useful, e.g. to construct regressions on distance matrices with nonlinearity and multiple random effects. More generally, this is a model for symmetric, relational data where the row and column effects are random and are integrated from the likelihood.

Initially, this repository housed a hacky script created for a single application (which can still be found [here](http://github.com/nspope/corMLPE_unsupported)). The old implementation is much slower, does not scale well to large datasets, and is only still extant for the sake of reproducibility. The easiest installation of the current version is via the package `devtools`;

```{r}
devtools::install_github("nspope/corMLPE")
```

Missing pairwise comparisons are allowed (e.g. one does not need a full set of pairwise measurements, as was the case for prior versions of this package), as are multiple observations from the same pair. However, observations that are self comparisons are not allowed.

The `corStruct` object allows a single grouping factor; for example to model isolation by distance within several species (where there are no pairwise measurements *between* species), an appropriate model might be

```{r}
lme(genetic.distance ~ geographic.distance, random = ~geographic.distance|species, 
    correlation=corMLPE(form=~pop1+pop2|species), data=my.data)
```

In this case `pop1` and `pop2` are numerical labels for the populations, such that each observation corresponds to a pair of populations.

With only a single species, use `gls`,

```{r}
gls(genetic.distance ~ geographic.distance, correlation=corMLPE(form=~pop1+pop2), data=my.data)
```

The package currently interfaces with `lme`, `nlme`, `gls`, and `gamm` (from package `mgcv`). The motivation for implementing this model in `nlme` was to leverage pre-existing machinery for nonlinear models, heteroskedasticity, and additive models. Look at Bill Peterman's package `ResistanceGA`[https://github.com/wpeterman/ResistanceGA] for a nice, `lme4`-based implementation.

Note that unlike the results presented in Clarke's paper, `nlme`/`lme` will return GLS standard errors rather than the OLS standard errors. If OLS standard errors are desired for some reason, see function `MLPE()` in the package.

You can reach me at nspope@utexas.edu if you have any questions/comments/complaints. Modify as you see fit.


# References
Clarke *et al*. 2002. Confidence Limits for Regression Relationships between Distance Matrices: Estimating Gene Flow with Distance. *Journal of Agricultural, Biological, and Environmental Statistics* 7: 361-372.

