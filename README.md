Understanding the stochastic partial differential equation approach to smoothing
================================================================================

David L. Miller, Richard Glennie & Andrew E. Seaton

Paper accepted at the Journal of Agricultural, Biological and Environmental Statistics.

The [paper](https://github.com/dill/SPDE-smoothing/blob/master/spde_paper.pdf) and [appendix](https://github.com/dill/SPDE-smoothing/blob/master/spde_paper_appendix.pdf) are here, code examples are in `supplementary/`.


## Abstract

Correlation and smoothness are terms used to describe a wide variety of random quantities. In time, space, and many other domains, they both imply the same idea: quantities that occur closer together are more similar than those further apart. Two popular statistical models that represent this idea are basis-penalty smoothers (Wood, 2017) and stochastic partial differential equations (SPDE) (Lindgren et al., 2011). In this paper, we discuss how the SPDE can be interpreted as a smoothing penalty and can be fitted using the R package mgcv, allowing practitioners with existing knowledge of smoothing penalties to better understand the implementation and theory behind the SPDE approach.
