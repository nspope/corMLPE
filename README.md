corMLPE
=======

This script containts construction functions of a corStruct class for the R package nlme, which are an implementation of Clarke's maximum likelihood population effects model. This is useful, ie. to construct regressions on distance matrices with nonlinearity and multiple random effects.

The corStruct object allows a single grouping factor; for example among several species of an organism, we could fit an isolation-by-distance model with a syntax such as

lme(genetic.distance ~ geographic.distance, random = ~geographic.distance|species, correlation=corMLPE(form=~popid1+popid2|species), data=my.data)

In this case "popid1" and "popid2" are numerical labels for the populations, such that each observation corresponds to a pair of populations.

Currently works with lme, nlme, gls, and gamm (from mgcv). The motivation for implementing in nlme was to allow compatibility with nonlinear models, heteroskedasticity, and additive models.

Note that unlike the results presented in Clarke's paper, nlme/lme will return GLS standard errors rather than the OLS standard errors. If OLS standard errors are desired for some reason, see function Clarke() in the script.

Reference:
Clarke et al. 2002. Confidence Limits for Regression Relationships between Distance Matrices: Estimating Gene Flow with Distance. Journal of Agricultural, Biological, and Environmental Statistics 7: 361-372.
