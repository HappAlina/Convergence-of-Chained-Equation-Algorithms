# Convergence-of-Chained-Equation-Algorithms
The default setting for the the number of FCS runs in mice (van Buuren &amp; Groothuis-Oudshoorn 2011) is rather low (maxit=5). Convergence in distribution is difficult to examine, but it is an important assumption for chained equation MI algorithms. Objectives: Conduct a simulation study to explore the effect of different missing data patterns and different quantities of interest on the convergence of the chained equations algorithm using traceplots and the Gelman-Rubin diagnostic Rhat (Gelman &amp; Rubin 1992).

- simulate dependent data
- create different missing patterns & mechanisms
- impute missing values with the mice package using the default setting of maxit=5
- evaluate convergence of mean & variance with Gelman-Rubin diagnostic Rhat
- implement an Rhat for the Rsquared of a multiple lineare regression to have a multivariate measure that takes dependencies between variables into account

Result: In more than half our cases the mean and variance estimators have not converged which shows that 5 iterations is often not enough. 
However, the Rhat for Rsquared showed good convergence behaviour in about 90 percent of the cases. Thus, dependencies between variables are not lost when imputing with mice (and a maxit of 5).
