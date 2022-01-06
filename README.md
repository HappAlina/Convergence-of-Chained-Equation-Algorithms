# Convergence-of-Chained-Equation-Algorithms
The default setting for the the number of FCS runs in mice (van Buuren &amp; Groothuis-Oudshoorn 2011) is rather low (maxit=5). Convergence in distribution is difficult to examine, but it is an important assumption for chained equation MI algorithms. Objectives: Conduct a simulation study to explore the effect of different missing data patterns and different quantities of interest on the convergence of the chained equations algorithm using traceplots and the Gelman-Rubin di- agnostic Rˆ (Gelman &amp; Rubin 1992).

To Do: 
- Create dataframe
- Apply different missing data patterns (MCAR, MAR, MNAR)
- Estimate missing values with mice by trying different maxit
- When does the increase of maxit not improve the estimation any further? Test with traceplot and Gelman Rubin Diagnostic (find point of convergence)

- create LaTeX Doc
- Write motivation
- Convergence in distribution theory
- GRD theory
