#! /usr/bin/env python3

# Integration test for `easel histplot` 
#
# Usage: easel-histplot-itest.py <builddir> <srcdir> <tmppfx>
#   <builddir>: path to Easel build dir. `easel` miniapp is <builddir>/miniapps/easel
#   <srcdir>:   path to Easel src dir.
#   <tmppfx>:   prefix we're allowed to use to create tmp files in current working dir.
#
import filecmp
import glob
import os
import re
import subprocess
import sys
import esl_itest

import math         # math.log()
import random       # random.random(), etc
import struct       # struct.pack() for writing binary data

progs_used = [ 'miniapps/easel' ]
files_used = [ ]

(builddir, srcdir, tmppfx) = esl_itest.getargs(sys.argv)
esl_itest.check_files(srcdir,   files_used)
esl_itest.check_progs(builddir, progs_used)

easel = f'{builddir}/miniapps/easel'

# `-h` help 
r = esl_itest.run(f'{easel} histplot -h')


# Generation of test data
# Either Gumbel-distributed, or normally-distributed
# In Easel itests, we need to do this with only the Python standard library.
# sample_gumbel() is derived from easel's esl_gumbel.c::esl_gumbel_Sample()
# `lambda` is a reserved word in Python, so the variable is called `lam`
#
def sample_gumbel(mu=0.0, lam=0.693, n=1):
    X = []
    for _ in range(n):
        while True:
            if (u := random.random()) > 0.0: break    # 0 > u <= 1.0  uniform positive 
        X.append( mu - ( math.log(-math.log(u)) / lam ) )  
    return X

def sample_normal(mu=0.0, sigma=1.0, n=1):
    X = []
    for _ in range(n):
        X.append(random.normalvariate(mu=mu,sigma=sigma))
    return X

# Make a test text datafile {tmppfx}.dat with 10000 Gumbel samples.
#
X = sample_gumbel(n=10000)
with open(f'{tmppfx}.dat', 'w') as f:
    for x in X:
        f.write(f'{x:.5f}\n')

# ... and a binary file with doubles, {tmppfx}.b
with open(f'{tmppfx}.b', 'wb') as f:
        for x in X:
            f.write(struct.pack('d', x))

# ... and a faux tabular file with the data in field 3, {tmppfx}.tbl
with open(f'{tmppfx}.tbl', 'w') as f:
    for x in X:
        f.write(f'a    b    {x:.5f}\n')

# ... and a truncated Gumbel, truncated at 0.0
with open(f'{tmppfx}.trunc', 'w') as f:
    for x in X:
        if x > 0.:
            f.write(f'{x:.5f}\n')


# basic
r = esl_itest.run(f'{easel} histplot {tmppfx}.dat')
if (m := re.search(r'^&', r.stdout, flags=re.MULTILINE)) is None: esl_itest.fail()
n = sum(int(line.split()[1]) for line in r.stdout.splitlines() if not line.startswith('&'))
if n != 10000: esl_itest.fail()

# -o   output file for plot
r2 = esl_itest.run(f'{easel} histplot -o {tmppfx}.xy {tmppfx}.dat')
with open(f'{tmppfx}.xy') as f: output = f.read()
if output != r.stdout: esl_itest.fail()

# -f   which field to use as data
r2 = esl_itest.run(f'{easel} histplot -f 3 {tmppfx}.tbl')
if r2.stdout != r.stdout: esl_itest.fail()

# -b   input file is binary, array of doubles
#  (we can't compare output exactly to r.stdout, because of floating point roundoff error)
#
r = esl_itest.run(f'{easel} histplot -b {tmppfx}.b')
n = sum(int(line.split()[1]) for line in r.stdout.splitlines() if not line.startswith('&'))
if n != 10000: esl_itest.fail()

#  -w   sets bin size for histogram
r = esl_itest.run(f'{easel} histplot -w 2.0 {tmppfx}.dat')
n = sum(int(line.split()[1]) for line in r.stdout.splitlines() if not line.startswith('&'))
if n != 10000: esl_itest.fail()

#  --min  initial lower bound of histogram  [-100.]
#  --max  initial upper bound of histogram  [100.]
#
# These are just hints to the ESL_HISTOGRAM, which resizes as needed. --min will reset the baseline,
# but --max has no effect, really.
#
r = esl_itest.run(f'{easel} histplot --min 1.0 {tmppfx}.dat')
n = sum(int(line.split()[1]) for line in r.stdout.splitlines() if not line.startswith('&'))
if n != 10000: esl_itest.fail()

r = esl_itest.run(f'{easel} histplot --max 1.0 {tmppfx}.dat')
n = sum(int(line.split()[1]) for line in r.stdout.splitlines() if not line.startswith('&'))
if n != 10000: esl_itest.fail()

#  --surv    : output survival plot, not histogram
#
r = esl_itest.run(f'{easel} histplot --surv {tmppfx}.dat')
has_one = False
for line in r.stdout.splitlines():
    if (m := re.search(r'-?\d+\.\d+\s+1\s*', line)) is not None:
        has_one = True
        break
if not has_one: esl_itest.fail()    

# --gumbel      fit data to a Gumbel distribution
# --exptail     fit tail to an exponential distribution
# --gev         fit data to a generalized EVD (Frechet or Weibull)
# --normal      fit data to a normal (Gaussian) distribution
# --gumloc      fit data to Gumbel of known lambda
# --exptailloc  fit data to location of exponential tail of known lambda
#
# -t <x>       : set tail mass to fit to  [0.01]
#
r = esl_itest.run(f'{easel} histplot --gumbel {tmppfx}.dat')
n = sum(float(line.split()[1]) for line in r.stdout.splitlines() if not line.startswith('&') and not line.startswith('#'))
if abs(n - 20000) > 2.0: esl_itest.fail()   # 10000 exactly from the data, and ~10000 from the fit +/- fp roundoff error

r = esl_itest.run(f'{easel} histplot --gev {tmppfx}.dat')
n = sum(float(line.split()[1]) for line in r.stdout.splitlines() if not line.startswith('&') and not line.startswith('#'))
if abs(n - 20000) > 2.0: esl_itest.fail()   

r = esl_itest.run(f'{easel} histplot --normal {tmppfx}.dat')
n = sum(float(line.split()[1]) for line in r.stdout.splitlines() if not line.startswith('&') and not line.startswith('#'))
if abs(n - 20000) > 2.0: esl_itest.fail()   

r = esl_itest.run(f'{easel} histplot --exptail {tmppfx}.dat')
if (m := re.match(r'# Exponential', r.stdout)) is None: esl_itest.fail()

r = esl_itest.run(f'{easel} histplot --exptail -t 0.02 {tmppfx}.dat')
if (m := re.match(r'# Exponential', r.stdout)) is None: esl_itest.fail()

r = esl_itest.run(f'{easel} histplot --gumloc --lambda 0.690 {tmppfx}.dat')
n = sum(float(line.split()[1]) for line in r.stdout.splitlines() if not line.startswith('&') and not line.startswith('#'))
if abs(n - 20000) > 2.0: esl_itest.fail()   

r = esl_itest.run(f'{easel} histplot --exptailloc {tmppfx}.dat')
r = esl_itest.run(f'{easel} histplot --exptailloc -t 0.02 {tmppfx}.dat')
# just making sure it doesn't fail completely

# --trunc <x>   with --gumbel, specify data are truncated, min value is <x>
#
r = esl_itest.run(f'{easel} histplot --gumbel --trunc 0.0 {tmppfx}.trunc')

# --showgum     plot a known Gumbel for comparison
# --showexp     plot a known exponential tail for comparison
# --showgev     plot a known GEV for comparison
# --alpha <x>   set known alpha (GEV shape parameter)  [0.0]
# --lambda <x>  set known lambda  [0.693]
# --mu <x>      set known mu  [0.0]
#
r = esl_itest.run(f'{easel} histplot --showgum --mu 0.01 --lambda 0.690              {tmppfx}.dat')
r = esl_itest.run(f'{easel} histplot --showexp --mu 4.00 --lambda 0.690 -t 0.02      {tmppfx}.dat')
r = esl_itest.run(f'{easel} histplot --showgev --mu 4.00 --lambda 0.690 --alpha 0.04 {tmppfx}.dat')

# Cleanup
for tmpfile in glob.glob(f'{tmppfx}.*'):
    os.remove(tmpfile)

print('ok')

