corMLPE
=======

This package implements a correlation structure for the R package `nlme` for Clarke's maximum likelihood population effects model (Clarke et al. 2002). This is useful, e.g. to construct regressions on distance matrices with nonlinearity and multiple random effects. More generally, this is a model for symmetric, relational data where the row and column effects are random and are integrated from the likelihood.

Initially, this repository housed a rough script created for a single application. I've since ported the computationally intensive part of the code to C++, extended it for 'larger' datasets, and packaged it for a smoother installation. The easiest installation is via the package `devtools`;

```{r}
devtools:::install_github("nspope/corMLPE")
```

**An important current limitation** is that the correlation structure assumes a complete set of pairwise distances; i.e. all possible pairwise combinations among samples. This assumption is necessary for fast computation of the likelihood, which relies on closed form for the spectrum of the model's correlation matrix. I'm currently working on relaxing this assumption.

The `corStruct` object allows a single grouping factor; for example among several species of an organism, we could fit an isolation-by-distance model with a syntax such as

```{r}
lme(genetic.distance ~ geographic.distance, random = ~geographic.distance|species, 
    correlation=corMLPE(form=~pop1+pop2|species), data=my.data)
```

In this case `pop1` and `pop2` are numerical labels for the populations, such that each observation corresponds to a pair of populations.

The package currently interfaces with `lme`, `nlme`, `gls`, and `gamm` (from package `mgcv`). The motivation for implementing this model in `nlme` was to leverage pre-existing machinery for nonlinear models, heteroskedasticity, and additive models. Look at Bill Peterman's package `ResistanceGA`[https://github.com/wpeterman/ResistanceGA] for a nice, `lme4`-based implementation.

Note that unlike the results presented in Clarke's paper, `nlme`/`lme` will return GLS standard errors rather than the OLS standard errors. If OLS standard errors are desired for some reason, see function `MLPE()` in the package.

You can reach me at `nspope@utexas.edu` if you have any questions/comments/complaints. Modify as you see fit.

# References
Clarke *et al*. 2002. Confidence Limits for Regression Relationships between Distance Matrices: Estimating Gene Flow with Distance. *Journal of Agricultural, Biological, and Environmental Statistics* 7: 361-372.
