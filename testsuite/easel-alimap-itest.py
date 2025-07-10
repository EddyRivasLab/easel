#! /usr/bin/env python3

# Integration test for `easel alimap` 
#
# Usage: easel-alimap-itest.py <builddir> <srcdir> <tmppfx>
#   <builddir>: path to Easel build dir. `easel` miniapp is <builddir>/miniapps/easel
#   <srcdir>:   path to Easel src dir.
#   <tmppfx>:   prefix we're allowed to use to create tmp files in current working dir.
#
import glob
import os
import re
import subprocess
import sys
import esl_itest

progs_used = [ 'miniapps/easel' ]

(builddir, srcdir, tmppfx) = esl_itest.getargs(sys.argv)
esl_itest.check_progs(builddir, progs_used)

easel = f'{builddir}/miniapps/easel'

# -h
r = esl_itest.run(f'{easel} alimap -h')

# create small MSA file test cases, amino and dna
with open(f'{tmppfx}.1', 'w') as f:
    f.write('# STOCKHOLM 1.0\n\n')
    f.write('seq1    AACDEFGHIKKKLMNPQRRSTVWYY\n')
    f.write('seq2    -ACDEFGHI--KLMNPQ-RSTVW-Y\n')
    f.write('seq3    -ACDEFGHI--KLMNPQ-RSTVW-Y\n')
    f.write('#=GC RF .xxxxxxxx..xxxxxx.xxxxx.x\n')
    f.write('//\n')

with open(f'{tmppfx}.2', 'w') as f:    
    f.write('# STOCKHOLM 1.0\n\n')
    f.write('seq1    AACDEFGHIKKKLMNPQRRSTVWYY\n')
    f.write('seq2    A-CDEFGHIK--LMNPQR-STVWY-\n')
    f.write('seq3    A-CDEFGHIK--LMNPQR-STVWY-\n')
    f.write('#=GC RF x.xxxxxxxx..xxxxxx.xxxxx.\n')
    f.write('//\n')

with open(f'{tmppfx}.3', 'w') as f:
    f.write('# STOCKHOLM 1.0\n\n')
    f.write('seq1 ACGTACGGGTACGTACGT\n')
    f.write('seq2 ACGTACG--TACGTACGT\n')
    f.write('seq3 ACGTACG--TACGTACGT\n')
    f.write('//\n')

with open(f'{tmppfx}.4', 'w') as f:
    f.write('# STOCKHOLM 1.0\n\n')
    f.write('seq1 ACGTACGGGTACGTACGT\n')
    f.write('seq2 ACGTAC--GTACGTACGT\n')
    f.write('seq3 ACGTAC--GTACGTACGT\n')
    f.write('//\n')

with open(f'{tmppfx}.sub', 'w') as f:
    f.write('# STOCKHOLM 1.0\n\n')
    f.write('seq1    ACDEFGHI\n')
    f.write('seq2    ACDEFGHI\n')
    f.write('seq3    ACDEFGHI\n')
    f.write('#=GC RF xxxxxxxx\n')
    f.write('//\n')

# basic
r = esl_itest.run(f'{easel} alimap {tmppfx}.1 {tmppfx}.2')
if re.search(r'^# RF coverage:\s+52\s+/\s+60', r.stdout, flags=re.MULTILINE) is None: esl_itest.fail()

# --amino
r2 = esl_itest.run(f'{easel} alimap --amino {tmppfx}.1 {tmppfx}.2')
if r.stdout != r2.stdout: esl_itest.fail()

# --dna
# --rna
r  = esl_itest.run(f'{easel} alimap       {tmppfx}.3 {tmppfx}.4')
r2 = esl_itest.run(f'{easel} alimap --dna {tmppfx}.3 {tmppfx}.4')
if r.stdout != r2.stdout: esl_itest.fail()

r2 = esl_itest.run(f'{easel} alimap --rna {tmppfx}.3 {tmppfx}.4')
if r.stdout != r2.stdout: esl_itest.fail()

# -q
r = esl_itest.run(f'{easel} alimap -q {tmppfx}.1 {tmppfx}.2')
if re.search(r'^[^#]', r.stdout, flags=re.MULTILINE) is not None: esl_itest.fail()   # with -q, all output lines start with #

# --mask-a2a
r = esl_itest.run(f'{easel} alimap --mask-a2a {tmppfx}.mask {tmppfx}.1 {tmppfx}.2')
if re.search(r'^#\s+\(Length:\s+25', r.stdout, flags=re.MULTILINE) is None: esl_itest.fail()

# --mask-a2rf
r = esl_itest.run(f'{easel} alimap --mask-a2rf {tmppfx}.mask {tmppfx}.1 {tmppfx}.2')
if re.search(r'^#\s+\(Length:\s+25', r.stdout, flags=re.MULTILINE) is None: esl_itest.fail()

# --mask-rf2a
r = esl_itest.run(f'{easel} alimap --mask-rf2a {tmppfx}.mask {tmppfx}.1 {tmppfx}.2')
if re.search(r'^#\s+\(Length:\s+20', r.stdout, flags=re.MULTILINE) is None: esl_itest.fail()

# --mask--rf2rf
r = esl_itest.run(f'{easel} alimap --mask-rf2rf {tmppfx}.mask {tmppfx}.1 {tmppfx}.2')
if re.search(r'^#\s+\(Length:\s+20', r.stdout, flags=re.MULTILINE) is None: esl_itest.fail()

# --submap
r = esl_itest.run(f'{easel} alimap --submap {tmppfx}.map {tmppfx}.1 {tmppfx}.sub')
with open(f'{tmppfx}.map') as f: map = f.read()
if re.fullmatch(r'0111111110000000000000000\s*\n', map) is None: esl_itest.fail()
    

# Cleanup
for tmpfile in glob.glob(f'{tmppfx}.*'):
    os.remove(tmpfile)


print('ok')


