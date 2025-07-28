#! /usr/bin/env python3

# Integration test for `easel filter` 
#
# Usage: easel-filter-itest.py <builddir> <srcdir> <tmppfx>
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

progs_used = [ 'miniapps/easel' ]
files_used = [ 'testsuite/example-rna.sto',
               'testsuite/example-stockholm.sto' ]

(builddir, srcdir, tmppfx) = esl_itest.getargs(sys.argv)
esl_itest.check_files(srcdir,   files_used)
esl_itest.check_progs(builddir, progs_used)

easel = f'{builddir}/miniapps/easel'

# `-h` help 
r = esl_itest.run(f'{easel} filter -h')

# basic
r  = esl_itest.run(f'{easel} filter 0.50 {srcdir}/testsuite/example-rna.sto')
r2 = subprocess.run(f'{easel} alipid -'.split(), check=True, encoding='utf-8', capture_output=True, input=r.stdout)
maxpid = max(float(line.split()[2]) for line in r2.stdout.splitlines() if not line.startswith('#'))
if maxpid > 50.: esl_itest.fail()

r3 = esl_itest.run(f'{easel} filter -o {tmppfx}.sto 0.50 {srcdir}/testsuite/example-rna.sto')
with open(f'{tmppfx}.sto') as f: s = f.read()
if s != r.stdout: esl_itest.fail()

# --seed
r3 = esl_itest.run(f'{easel} filter --seed 42 0.50 {srcdir}/testsuite/example-rna.sto')
if r3.stdout != r.stdout: esl_itest.fail()   # --seed only affects --randorder, and subsampling of deep MSAs

# --rna
# --dna
r3  = esl_itest.run(f'{easel} filter --rna 0.50 {srcdir}/testsuite/example-rna.sto')
if r3.stdout != r.stdout: esl_itest.fail()    # compare MSAs

r3  = esl_itest.run(f'{easel} filter --dna 0.50 {srcdir}/testsuite/example-rna.sto')
r4 = subprocess.run(f'{easel} alipid -'.split(), check=True, encoding='utf-8', capture_output=True, input=r3.stdout)
if r4.stdout != r2.stdout: esl_itest.fail()   # compare alipid outputs, not the MSAs; --dna changes T to U in MSA

# --informat
# --outformat
r3  = esl_itest.run(f'{easel} filter --informat stockholm 0.50 {srcdir}/testsuite/example-rna.sto')
if r3.stdout != r.stdout: esl_itest.fail()    

r3  = esl_itest.run(f'{easel} filter --outformat stockholm 0.50 {srcdir}/testsuite/example-rna.sto')
if r3.stdout != r.stdout: esl_itest.fail()    

# options for deriving consensus:
#  --ignore-rf      : ignore any RF line; always determine our own consensus
#  --fragthresh <x> : seq is fragment if aspan/alen < fragthresh  [0.5]
#  --symfrac <x>    : col is consensus if nres/(nres+ngap) >= symfrac  [0.5]
#
r3 = esl_itest.run(f'{easel} filter --ignore-rf --fragthresh 0.6 --symfrac 0.6 0.50 {srcdir}/testsuite/example-rna.sto')
r4 = subprocess.run(f'{easel} alipid -'.split(), check=True, encoding='utf-8', capture_output=True, input=r3.stdout)
maxpid = max(float(line.split()[2]) for line in r4.stdout.splitlines() if not line.startswith('#'))
if maxpid > 50.: esl_itest.fail()

# options for deriving consensus by sampling (on deep MSAs):
#  --no-sampling    : never use subsampling to determine consensus
#  --nsamp <n>      : number of seqs to sample (if using sampling)  [10000]
#  --sampthresh <n> : switch to using sampling when nseq > nsamp  [50000]
#  --maxfrag <n>    : if sample has > maxfrag fragments, don't use sample  [5000]
#
r3 = esl_itest.run(f'{easel} filter --ignore-rf --no-sampling 0.50 {srcdir}/testsuite/example-rna.sto')
if r3.stdout != r.stdout: esl_itest.fail()    # --no-sampling has no effect on shallow MSAs. This isn't a very good test. Need deep MSA to trigger sampling.

r3 = esl_itest.run(f'{easel} filter --ignore-rf --nsamp 50 --sampthresh 100 0.50 {srcdir}/testsuite/example-rna.sto')  # original #=RF line is preserved on output,
r4 = subprocess.run(f'{easel} alipid -'.split(), check=True, encoding='utf-8', capture_output=True, input=r3.stdout)   # even though a new consensus is calculated w/ --ignore-rf
maxpid = max(float(line.split()[2]) for line in r4.stdout.splitlines() if not line.startswith('#'))
if maxpid > 50.: esl_itest.fail()

r3 = esl_itest.run(f'{easel} filter --ignore-rf --maxfrag 1 0.50 {srcdir}/testsuite/example-rna.sto')
r4 = subprocess.run(f'{easel} alipid -'.split(), check=True, encoding='utf-8', capture_output=True, input=r3.stdout) 
maxpid = max(float(line.split()[2]) for line in r4.stdout.splitlines() if not line.startswith('#'))
if maxpid > 50.: esl_itest.fail()

# options for sequence preference:
#  --conscover : keep seq whose alispan has better consensus coverage  [default]
#  --randorder :  ... or with random preference
#  --origorder :  ... or prefer seq that comes first in order
#
r3 = esl_itest.run(f'{easel} filter --conscover 0.50 {srcdir}/testsuite/example-rna.sto')
if r3.stdout != r.stdout: esl_itest.fail()

r3 = esl_itest.run(f'{easel} filter --randorder 0.50 {srcdir}/testsuite/example-rna.sto')
r4 = subprocess.run(f'{easel} alipid -'.split(), check=True, encoding='utf-8', capture_output=True, input=r3.stdout) 
maxpid = max(float(line.split()[2]) for line in r4.stdout.splitlines() if not line.startswith('#'))
if maxpid > 50.: esl_itest.fail()

r3 = esl_itest.run(f'{easel} filter --origorder 0.50 {srcdir}/testsuite/example-rna.sto')
r4 = subprocess.run(f'{easel} alipid -'.split(), check=True, encoding='utf-8', capture_output=True, input=r3.stdout) 
maxpid = max(float(line.split()[2]) for line in r4.stdout.splitlines() if not line.startswith('#'))
if maxpid > 50.: esl_itest.fail()

# --amino           ... that input MSA is protein
r  = esl_itest.run(f'{easel} filter 0.30 {srcdir}/testsuite/example-stockholm.sto')
r2 = subprocess.run(f'{easel} alipid -'.split(), check=True, encoding='utf-8', capture_output=True, input=r.stdout)
maxpid = max(float(line.split()[2]) for line in r2.stdout.splitlines() if not line.startswith('#'))
if maxpid > 50.: esl_itest.fail()

r3  = esl_itest.run(f'{easel} filter --amino 0.30 {srcdir}/testsuite/example-stockholm.sto')
if r3.stdout != r.stdout: esl_itest.fail()

# Cleanup
for tmpfile in glob.glob(f'{tmppfx}.*'):
    os.remove(tmpfile)

print('ok')

                                                   
