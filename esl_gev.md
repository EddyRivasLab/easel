# esl_gev - Generalized extreme value distributions

The generalized extreme value distribution (GEV) includes all three
types of extreme value distributions: Type I (Gumbel), type II
(Fréchet), and type III (Weibull). Empirically, the scores of some
sequence alignment algorithms appear to follow GEV distributions. The
`gev` module is used in estimating the statistical significance
of such scores.

Most local sequence alignment scores follow the Gumbel distribution.
Easel's `gumbel` module applies specifically to the Gumbel. The
`gev` module is used for Type II or III extreme value
distributions, or for determining which of the three types of
distribution that a dataset best fits.

The generalized extreme value distribution depends on three
parameters, $\mu$, $\lambda$, and $\alpha$. When these parameters are
known, the statistical significance (P-value) of a single score $x$ is
$P(S>x)$, obtained by a call to `esl_gev_surv()`. The E-value for
obtaining that score or better in searching a database of $N$
sequences is just $NP(S>x)$.

When the parameters are unknown, they can be estimated from scores
obtained from comparisons of simulated random data. The
`esl_gev_FitComplete()` function performs maximum likelihood
parameter fitting [Coles01].

The `eslGEV_EXAMPLE` example at the end of `esl_gev.c` samples 10,000
data points from a Fréchet distribution with $\mu=-20$, $\lambda=0.4$,
$\alpha=0.1$; reports the min and max samples, and the probability
mass to the left of the min and to the right of the max (both of which
should be about $\frac{1}{10000}$, since we took 10,000 samples); and
then fits those simulated data to a Gumbel and reports the fitted
$\mu$ and $\lambda$.

## GEV densities

The probability density function (pdf) and the cumulative distribution
function (cdf) of the generalized extreme value distribution are
[Coles01]:

$$
\begin{aligned}
P(X=x) & = \lambda \left[ 1 + \alpha \lambda (x - \mu) \right]^{-\frac{\alpha+1}{\alpha}}
       \exp \left\{ - \left[ 1 + \alpha \lambda (x - \mu)
       \right]^{-\frac{1}{\alpha}} \right\} \\
P(X \geq x) & = \exp \left\{ - \left[ 1 +
     \alpha\lambda(x-\mu) \right]^{-\frac{1}{\alpha}} \right\}
\end{aligned}
$$

The parameters $\mu$, $\lambda$, and $\alpha$ are location, scale, and
shape parameters, respectively, with $-\infty < \mu < \infty$, $0 <
\lambda < \infty$, and $-\infty < \alpha < \infty$.

The Type II (Fréchet) distribution corresponds to $\alpha > 0$,
and the Type III (Weibull) distribution corresponds to $\alpha < 0$.
The Type I (Gumbel) distribution arises in the limit $\alpha
\rightarrow 0$. At values $\alpha \simeq 0$, Easel's GEV functions
revert to the Gumbel limit case, as opposed to dividing by zero and
failing.

Technically the GEV is only defined for values of $x$ such that $1 +
\alpha \lambda (x - \mu) > 0$. However, Easel's functions return
sensible values outside this domain, such as 0 for nonexistent
densities.

For more details, see the excellent description in [Coles01].
Easel's $\{ \mu, \lambda, \alpha \}$ notation differs from the $\{
\mu, \sigma, \xi \}$ parameterization used by Coles. Use $\lambda =
\frac{1}{\sigma}$ and $\alpha = \xi$ to translate.

## fitting GEV distributions to observed data

Easel fits GEVs by maximum likelihood estimation by numerically
optimizing the log likelihood function, using first derivative
information and conjugate gradient descent. See the `gumbel`
chapter for a more general introduction to maximum likelihood fitting.

### maximum likelihood estimation, complete data

The function `esl_gev_FitComplete()` uses gradient information
to find parameters that optimize the likelihood function, using the
conjugate gradient descent code in the `minimizer` module.

Given $n$ samples $x_1..x_n$, we want to estimate maximum likelihood
parameter estimates $\{ \hat{\mu}, \hat{\lambda}, \hat{\alpha} \}$
that maximize the log likelihood:

$$
\log L(\lambda, \mu, \alpha) = n \log \lambda 
       - \frac{\alpha+1}{\alpha} 
           \sum_{i=1}^{n} \log\left[1+ \alpha\lambda(x_i - \mu) \right]
       - \sum_{i=1}^{n} \left[ 1 + \alpha\lambda (x_i - \mu) \right]^{\frac{1}{\alpha}}
$$

The $\left[ 1 + \alpha\lambda (x_i - \mu) \right]^{\frac{1}{\alpha}}$
term can be rewritten in a more conveniently differentiable form as

$$
\exp \left\{ \frac{1}{\alpha} \log \left[ 1 + \alpha\lambda (x_i - \mu) \right] \right\}.
$$

Since the $\lambda$ parameter is constrained to $\lambda > 0$ but the
numerical optimizer expects unconstrained parameters, we use a change
of variables $\lambda = e^w$ and optimize an unconstrained value $w$.

The gradient of the log likelihood with respect to $\mu$, $w$, and
$\alpha$ is:

$$
\begin{aligned}
\frac{\partial \log L}{\partial \mu} & =
  \sum_{i=1}^n \frac{\lambda (\alpha+1)}{1+\alpha\lambda(x_i-\mu)} 
 -\sum_{i=1}^n \lambda \exp 
    \left\{ -\frac{\alpha+1}{\alpha} \log
          \left[1+\alpha\lambda(x_i-\mu)\right] \right\}
\\
\frac{\partial \log L}{\partial w} & =
  n - \sum_{i=1}^{n} \frac{\lambda (\alpha+1) (x_i - \mu)} 
                          {1 + \alpha \lambda (x_i - \mu)}
  + \sum_{i=1}^n \lambda (x_i - \mu) 
         \exp \left\{ -\frac{\alpha+1}{\alpha} \log
          \left[1+\alpha\lambda(x_i-\mu)\right] \right\}
\\
\frac{\partial \log L}{\partial \alpha} & =
   \sum_{i=1}^n \left\{
      \begin{array}{l}
      - \frac{\alpha+1}{\alpha} \frac{\lambda(x_i-\mu)}
                                  {1 +\alpha\lambda(x_i-\mu)}\\
      + \frac{1}{\alpha^2} \log \left[ 1 + \alpha\lambda(x_i - \mu) \right]\\
      + \frac{1}{\alpha} \frac{\lambda(x_i-\mu)}
                          {1 +\alpha\lambda(x_i-\mu)}
      e^{-\frac{1}{\alpha} \log\left[ 1 + \alpha\lambda(x_i - \mu) \right]}\\
     -  \frac{1}{\alpha^2} \log \left[ 1 + \alpha\lambda(x_i - \mu) \right]
      e^{-\frac{1}{\alpha} \log\left[ 1 + \alpha\lambda(x_i - \mu)
	 \right]} 
     \end{array}
     \right.
\end{aligned}
$$

When $|\alpha\lambda(x_i - \mu)|$ approaches $0$, the GEV approximates
a Gumbel distribution and these equations can be simplified using the
approximation $\log(1+a) \simeq a$.


## References

* **[Coles01]** S Coles, [_An Introduction to Statistical Modeling of
  Extreme Values_](https://search.worldcat.org/search?q=isbn:1852334592),
  Springer-Verlag (2001).
