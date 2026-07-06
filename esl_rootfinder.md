# esl_rootfinder - One-dimensional rootfinding

To use the bisection method, you provide a function that computes
$f(x)$, and you provide an initial bracket $(x_L..x_R)$ in which the
root is known to lie. The `eslROOTFINDER_EXAMPLE` example at the end
of `esl_rootfinder.c` uses the bisection method to compute a root of
the quadratic function $ax^2 + bx + c = 0$ for $a=5$, $b=2$ and $c=-1$
(note that this function has two roots, one $<0$ and one $>0$; the
bracket $0..100$ makes the example find only the positive root).

To use the Newton/Raphson method instead of bisection, you provide a
function that computes $f(x)$ and its first derivative $df(x)/dx$, and
you provide an initial guess for the root $x$. The
`eslROOTFINDER_EXAMPLE2` example uses the Newton/Raphson method to
compute the root of the same function above (which has a derivative
$df(x)/dx = 2ax + b$). In this example, because the initial guess was
negative, the other root gets found.

Currently, just these two rootfinding algorithms are implemented. The
bisection method does not require derivative information, and it
requires the caller to provide an interval $(x_L..x_R)$ in which the
root lies ($f(x_L)$ and $f(x_R)$ have opposite signs). Newton/Raphson
uses derivative information, and it only needs an initial guess for
$x$, not an interval. Thus there are two different `_Create*()`
routines, `esl_rootfinder_CreateBracketer()` for initializing a
bisection method, and `esl_rootfinder_CreatePolisher()` for
initializing a Newton/Raphson method. The reason for the more general
names (`CreateBracketer()` and `CreatePolisher()`) is that I expect
other rootfinding algorithms (if we ever implement any) will group
similarly: bracketing methods without derivative information, and
"polishing" methods that use derivative information. But this may be
misguided, and may change in the future.
