# esl_stopwatch - timing parts of programs

The `stopwatch` module measures the elapsed (wall clock) time, CPU
time, and system time consumed by any part of a program.

The simple way to measure the CPU time consumption in an ANSI C program
is:

```c
    clock_t  t0, t1;
    t0 = clock();
    /* do_stuff */
    t1 = clock();
    printf("cpu time: %.2f\n", (double) (t1-t0)/(double) CLOCKS_PER_SEC);
```

The stopwatch module is just an elaboration of this. It tracks elapsed
and system time, in addition to cpu time; it hides the details of
converting a time difference in hardware clock ticks to a
human-interpretable time in seconds; and it provides a standard output
function for formatting times, similar to the output of the standard
UNIX `time` command line utility for timing processes.

Starting a stopwatch with `esl_stopwatch_Start()` initializes a base
time, t0. Stopping a stopwatch with `esl_stopwatch_Stop()` takes the
current time t1, and internally computes and stores elapsed, cpu, and
system time differences (t1-t0). These stored times can be displayed at
any time using `esl_stopwatch_Display()`, until the next time the watch
is stopped. A stopwatch can be stopped any number of times, measuring
increasing time from the same base. A stopwatch can also be started any
number of times, resetting the base each time it is set.

The `eslSTOPWATCH_EXAMPLE` example at the end of `esl_stopwatch.c`
measures a boring `sleep(5)` call, which should of course show an
elapsed wall time of 5 seconds. Change the `sleep(5)` call to
something cpu- or system-intensive to see a non-zero measurement of
cpu or system time.


## displaying and retrieving times

The `esl_stopwatch_Display()` function prints a line containing the cpu
time, system time, aggregated cpu+system time, and the elapsed (wall
clock) time. For example:

```c
CPU Time: 142.55u 7.17s 00:02:29.72 Elapsed: 00:02:35
```

If you want to access the times in seconds for your own purposes, the
relevant fields in a stopped `ESL_STOPWATCH` object are:

```c
  double elapsed;               /* elapsed time, seconds */
  double user;                  /* CPU time, seconds     */
  double sys;                   /* system time, seconds  */
```

## stopwatch precision and system dependency

Elapsed wall time is typically measured at low resolution, in units of
seconds (depending on the ANSI C `time_t` definition on your
system). It is displayed with a precision of 1 sec.

CPU time is typically measured in high resolution, in units of
microseconds (depending on the value of POSIX `_SC_CLK_TCK` or ANSI C
`CLOCKS_PER_SEC` on your system). It is displayed with a precision of
0.01 sec.

System time is only determined on systems that provide a POSIX
`times()` function. Like CPU time, it is typically measured at high
resolution, in units of microseconds (depending on the POSIX
`_SC_CLK_TCK` value on your system). It is displayed with a precision of
0.01 sec. On systems that do not provide a POSIX-compliant `times()`
function, system time is always reported as 0.

## aggregate times in parallelized code

In parallelized code, you may want to aggregate results from multiple
stopwatches into a single overall time measurement. Examples include
aggregating times from worker processes in PVM or MPI applications, or
aggregating times from multiple execution threads on systems where the
`times()` function does not correctly aggregate threads for you.

The `esl_stopwatch_Include()` function adds the cpu and system times in
a "client" stopwatch to a "master" stopwatch. Both the client and the
master stopwatch must be stopped. The elapsed time in the master
stopwatch is not affected; it is assumed to be keeping track of the
real (wall clock) time.
