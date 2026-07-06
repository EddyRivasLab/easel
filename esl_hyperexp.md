# esl_hyperexp - Hyperexponential distributions

The hyperexponential (mixture exponential) distribution may be useful
for fitting fat-tailed empirical distributions.

## hyperexponential densities

The hyperexponential distribution is a mixture of $K$ independent
exponentials with a common location $\mu$ and different decay
constants $\lambda_k$.

The probability density function (PDF) is:

$$
P(X=x) = \sum_k^{K} q_k \lambda_k e^{- \lambda_k (x - \mu)}
$$

The cumulative distribution function (CDF) is:

$$
P(X \leq x) = \sum_k^{K} q_k (1 - e^{- \lambda_k (x - \mu)})
$$

Variate $x$ ranges $\mu \leq x < \infty$.

Mixture coefficients $q_k$ specify the prior probability of each
component $k$; $0 \leq q_k \leq 1$ and $\sum_k q_k = 1$.

The single location parameter $\mu$ is unconstrained, $-\infty < \mu <
\infty$. (Exponential distributions are usually represented without an
explicit location parameter, implicitly assuming $\mu = 0$.)

The scale parameters $\lambda_k$ for each component are nonnegative,
$\lambda_k > 0$.
