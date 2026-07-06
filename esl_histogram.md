# esl_histogram - collecting and displaying histograms

The `histogram` module is for collecting scores, fitting them to
expected distributions, and displaying them.

The histogram automatically reallocates its bins as data points
arrive, so the caller only needs to provide some initial guidance
about bin size and "phase" (offset of the bins relative to the real
number line). It accumulates counts in 64-bit unsigned integers, so it
can handle over $10^{19}$ total counts. Optionally (and provided that
the caller knows it has enough memory to support this), a "full"
histogram can be created and used to collect a sorted vector of raw
(unbinned) values.

Various different ways of fitting histogram data to different sorts of
expected distributions are supported, with interfaces to all of
Easel's statistical distribution modules. Data fitting is oriented
toward the case where the values are scores, with high scores being of
the most interest; for instance, routines for obtaining and fitting
the right (high-scoring) tail are provided, but not for the left tail.

Several of the output functions output data as XY data files suitable
for input into the freely available
[xmgrace](http://plasma-gate.weizmann.ac.il/Grace/) graphing program.


## example

The `eslHISTOGRAM_EXAMPLE` example at the end of `esl_histogram.c`
stores 10,000 samples from a Gumbel distribution in a histogram,
retrieves a vector containing the sorted samples, fits a Gumbel
distribution to that dataset, sets the expected counts in the
histogram, prints the observed and expected counts in an ASCII
histogram, and evaluates the goodness-of-fit.

Some points of interest:

* When the histogram is created, the arguments `(-100, 100, 0.5)` tell
  it to bin data into bins of width 0.5, initially starting at -100
  and ending at 100. This initialization is described below (see
  "Specifying binning of data values").

* Samples are collected one at a time with `esl_histogram_Add()`.

* After the data have been collected in a _full_ histogram, a vector
  of sorted raw data values can be retrieved using functions like
  `esl_histogram_GetData()`, and used to fit parameters of an expected
  distribution to the data.

* In addition to the observed binned counts, you can optionally set
  _expected_ binned counts in the histogram by calling
  `esl_histogram_SetExpect()` and providing pointers to an appropriate
  distribution function and its parameters.

* The `esl_histogram_Print()` function shows an ASCII text
  representation of the observed counts (and expected counts, if set)
  that looks a lot like FASTA's nice histogram output.

* The `esl_histogram_Goodness()` function compares the observed and
  expected binned counts, and calculates two goodness of fit tests: a
  G-test, and a $\chi^2$ test.


## specifying binning of data values

The histogram collects data values into bins. When the histogram is
created, the bin width and the relative offset of the bins is
permanently set, and an initial range is allocated.

For example, the call `esl_histogram_Create(-10, 10, 0.5)` creates 40
bins of width 0.5 from -10 to 10, with the first bin collecting scores
from $-10 < x \leq -9.5$, and the last bin collecting scores $9.5 < x
\leq 10.0$.

The lower bound of the initialization permanently sets the relative
offset of the bins. That is, `esl_histogram_Create(-10, 10, 0.5)`
makes the first bin $-10 < x \leq -9.5$, whereas
`esl_histogram_Create(-10.1, 9.9, 0.5)` makes the first bin $-10.1 < x
\leq -9.6$.

Aside from that, the initial range is only a suggestion. You can add
any real-valued $x$ to the histogram. The histogram will silently
reallocate itself to a wider range as needed. The ability of a
histogram to store data is effectively unlimited. Up to $2^{64}-1$
(more than $10^{19}$) counts can be collected. The histogram requires
16 bytes of storage per bin, and the number of bins it allocates
scales as $x_{\mathrm{max}} - x_{\mathrm{min}} / w$.


## optional collection of raw data values: full histograms

Normally a histogram would store only binned counts, so it can
efficiently summarize even very large numbers of samples.

In some cases it is useful to keep a list of the raw data values --
for instance, for more accurate parameter fitting to expected
distributions. This can be done by creating a "full" histogram with
`esl_histogram_CreateFull()` instead of `esl_histogram_Create()`. (The
example code above did this, because it did parameter fitting to the
raw data.) After data have been collected in a full histogram,
individual raw values or pointers to sorted arrays of raw values can
be retrieved using the `esl_histogram_Get*` functions.

A full histogram may require much more memory: about 4 bytes per data
point. You may not want to use full histograms if your problem
involves collecting many ($> 10^9$, say) data points.


## different parameter fitting scenarios

By default, the data you collect are assumed to be _complete_. You
observed all samples; if you fit to any expected distribution, the
expected distribution is assumed to describe the complete data; the
parameters of the expected distribution are to be fitted to an array
of the complete raw data samples; and any goodness of fit test is to
be applied to the complete data. This is the simplest, most obvious
case.

Other situations may arise. In addition to complete data, Easel is
designed to deal with four other cases:

1. The collected data are complete, and they are fit to a distribution
   that describes the complete data, but parameter fitting is done
   only in the right (highest-scoring) tail. This makes parameter
   fitting focus on the most important, high-scoring region of a score
   distribution, and ignore low-scoring outliers.

2. The collected data are complete, but they are fit to a distribution
   that only describes the right (highest scoring) tail, and the
   goodness-of-fit test is only performed on that tail. This case
   arises when we don't know the form of the expected distribution for
   the complete data, but the tail follows a predictable decay (an
   exponential tail, for example).

3. The collected data are left-censored such that no values $< \phi$
   were recorded in the histogram, but the data are fit to a complete
   distribution that predicts the probability even of the censored
   (unobserved) values. Goodness of fit is only evaluated in the
   observed data. (This case is what is actually meant by
   left-censored data.)

4. The high-scoring right tail of the collected data are fit as the
   _binned_ counts in the histogram (not raw sample values) to a
   distribution that describes the tail, such as an exponential. This
   case becomes useful when the raw data values have limited precision
   (because of rounding, for example), which can cause numerical
   problems with parameter fitting to tails. Another case where this
   is useful is when there are so many data points that the data must
   be binned just as a matter of practicality (not enough memory to
   hold a full histogram).

A variety of other situations can be dealt with by using different
combinations of the function calls that deal with these four cases.


### focusing parameter fitting on the highest scores

The `eslHISTOGRAM_EXAMPLE2` example at the end of `esl_histogram.c`
shows how to focus a Gumbel parameter fit on the right (high-scoring)
half of an observed distribution.

The key differences from the complete data case are:

* Only the high-scoring 50% of the data samples are retrieved, by
  calling `esl_histogram_GetTailByMass(h, 0.5, &xv, &n, &z)`. This
  returns `z`, the number of samples that were _censored_.

* These data are fit to a Gumbel distribution as a _left-censored_
  dataset by calling `esl_gumbel_FitCensored(xv, n, z, xv[0], &mu,
  &lambda)`.

The expected counts and the goodness of fit tests are still evaluated
for the complete data, even though the fit was performed only on the
highest scores.


### fitting to a tail distribution

The `eslHISTOGRAM_EXAMPLE3` example at the end of `esl_histogram.c`
fits an exponential tail to the high-scoring 10% of a
Gumbel-distributed dataset.

The differences to note are:

* The tail is fit as if it is _complete_ data as far as the
  exponential distribution is concerned.

* As a result, to use the exponential tail to predict expected data,
  we have to keep in mind how much probability mass the tail is
  supposed to predict (here, 10%), and that is provided to
  `esl_histogram_SetExpectedTail()`, which specifically calculates
  expected counts for a tail.


### fitting left-censored data

Fitting a Gumbel distribution to data that are _truly_ left censored
looks a lot like the case where we extracted the high scoring data for
a censored fit, as shown by the `eslHISTOGRAM_EXAMPLE4` example at the
end of `esl_histogram.c`.


### fitting binned data to a tail distribution

Normally, you want to fit parameters to the actual individual data
samples, not to binned data, because you'll get more accurate results.
An exception can arise when the data samples have limited precision
because they've been rounded off. Most distributions are not sensitive
to this, but some tail densities are, especially those with
singularities ($P(X=x) \rightarrow \infty$) at their origin. In such a
case, a fit to binned data may be superior, especially if you can
match the histogram's bins to the rounding procedure that was used.

The `eslHISTOGRAM_EXAMPLE5` example at the end of `esl_histogram.c`
shows fitting for samples that were already rounded up to the nearest
integer before adding them to the histogram.

Issues to note:

* The `esl_histogram_Create(-100, 100, 1.0)` call defined bins that
  exactly match the rounding procedure defined by `ceil(x)` -- all $x$
  that are rounded to the same value by `ceil(x)` would also go in the
  same bin of the histogram.

* The `esl_histogram_SetTailByMass()` function sets flags in the
  histogram to demarcate the desired tail. However, because the data
  have been binned, and we can only define the tail by a range of
  bins, it will generally be impossible to match the requested tail
  mass with adequate accuracy; the actual tail mass is $\geq$ the
  requested tail mass. It is returned to the caller, and it is the
  actual mass, not the requested mass, that should be used when
  setting expected counts.

* The `esl_histogram_SetRounding()` declaration sets a flag in the
  histogram that tells binned parameter fitting functions that the
  origin of the fitted density ($\mu$) should be set at the lower
  bound of the smallest bin, rather than the smallest raw data value
  observed in that bin.
