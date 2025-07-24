#! /usr/bin/env python3

# Integration test for `easel kmer` 
#
# Usage: easel-kmer-itest.py <builddir> <srcdir> <tmppfx>
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
files_used = [ 'testsuite/example-genbank.gb',
               'testsuite/example-uniprot.dat', ]

(builddir, srcdir, tmppfx) = esl_itest.getargs(sys.argv)
esl_itest.check_files(srcdir,   files_used)
esl_itest.check_progs(builddir, progs_used)

easel = f'{builddir}/miniapps/easel'

# `-h` help 
#
r = esl_itest.run(f'{easel} kmer -h')

# basic
#
r  = esl_itest.run(f'{easel} kmer 3 {srcdir}/testsuite/example-genbank.gb')
nlines = len(r.stdout.splitlines())
nkmer  = sum(int(line.split()[1]) for line in r.stdout.splitlines() if not line.startswith('#'))
if nlines != 4**3 + 2: esl_itest.fail()   # k=3 => 4^k k-mers, +2 header lines
if nkmer  != 72202:    esl_itest.fail()

# --dna
#
r2  = esl_itest.run(f'{easel} kmer --dna 3 {srcdir}/testsuite/example-genbank.gb')
if r2.stdout != r.stdout: esl_itest.fail()

# --rna
#
r2  = esl_itest.run(f'{easel} kmer --rna 3 {srcdir}/testsuite/example-genbank.gb')
nlines = len(r2.stdout.splitlines())
nkmer  = sum(int(line.split()[1]) for line in r2.stdout.splitlines() if not line.startswith('#'))
if nlines != 4**3 + 2: esl_itest.fail()  
if nkmer  != 72202:    esl_itest.fail()

# --amino      
#
r  = esl_itest.run(f'{easel} kmer --amino 1 {srcdir}/testsuite/example-uniprot.dat')
nlines = len(r.stdout.splitlines())
if nlines != 22: esl_itest.fail()  

# --informat 
r2 = esl_itest.run(f'{easel} kmer --informat uniprot 1 {srcdir}/testsuite/example-uniprot.dat')
if r2.stdout != r.stdout: esl_itest.fail()

# -d         double-stranded
r =  esl_itest.run(f'{easel} kmer -d 3 {srcdir}/testsuite/example-genbank.gb')
nkmer = sum(int(line.split()[1]) for line in r.stdout.splitlines() if not line.startswith('#'))
if nkmer != 144404:  esl_itest.fail()    

# -q             : quiet: suppress column headers
r  = esl_itest.run(f'{easel} kmer -q 3 {srcdir}/testsuite/example-genbank.gb')
nlines = len(r.stdout.splitlines())
if nlines != 4**3: esl_itest.fail()  


  
# Cleanup
for tmpfile in glob.glob(f'{tmppfx}.*'):
    os.remove(tmpfile)

print('ok')
