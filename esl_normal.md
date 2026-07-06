# esl_normal - Normal (Gaussian) distributions

| Quantity | Symbol | Type | Range |
|----------|--------|------|-------|
| Variate | $x$ | `double` | $-\infty < x < \infty$ |
| Location | $\mu$ | `double` | $-\infty < \mu < \infty$ |
| Scale | $\sigma$ | `double` | $\sigma > 0$ |

The probability density function (PDF) is:

$$
\mathrm{PDF} = P(X=x) =  \frac{1}{\sigma \sqrt{2\pi}} e^{\frac{-(x-\mu)^2}{2\sigma^2}}.
$$

The cumulative distribution function (CDF) does not have a convenient
closed-form expression. It is derived numerically in terms of the
error function, $\mathrm{erf}()$:

$$
\mathrm{CDF} = P(X<x) =  \frac{1}{2} + \frac{1}{2} \mathrm{erf}\left(\frac{x - \mu}{\sigma \sqrt{2}}\right).
$$
