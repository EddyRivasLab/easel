#! /usr/bin/env python3

# Integration test for `easel compalign` 
#
# Usage: easel-compalign-itest.py <builddir> <srcdir> <tmppfx>
#   <builddir>: path to Easel build dir. `easel` miniapp is <builddir>/miniapps/easel
#   <srcdir>:   path to Easel src dir.
#   <tmppfx>:   prefix we're allowed to use to create tmp files in current working dir.
#
# The tests here are regression tests, not nasty tests that explore edge cases.
#
import glob
import os
import re
import subprocess
import sys
import esl_itest

progs_used = [ 'miniapps/easel' ]
files_used = [ 'testsuite/example-rna.sto' ]

(builddir, srcdir, tmppfx) = esl_itest.getargs(sys.argv)
esl_itest.check_files(srcdir,   files_used)
esl_itest.check_progs(builddir, progs_used)

easel = f'{builddir}/miniapps/easel'

# `-h` help 
r = esl_itest.run(f'{easel} compalign -h')

# create small test case that works for both amino and dna/rna
testmsa1 = """\
# STOCKHOLM 1.0

seq1     GGGAAA.CCaC
seq2     GGGAAAACC.C
seq3     GGG.AA.CC.C
#=GC RF  xxxxxxxxx.x
//
"""

testmsa2 = """\
# STOCKHOLM 1.0

seq1          GGG.AAACCaC
#=GR seq1 PP  999.9999999
seq2          GGGAAAACC.C
#=GR seq2 PP  999999999.9
seq3          GGG..AACC.C
#=GR seq3 PP  999..9999.9
#=GC RF       xxxxxxxxx.x
//
"""

testmask = "1110000000"    # for consensus RF columns; change mask if you change the number of consensus cols
rflen = len(testmask)   
alen  = rflen+1            # +1 because there's one insert col in examples above. Change this if you change the example MSA length.

with open(f'{tmppfx}.1.sto', 'w') as f:  f.write(testmsa1)
with open(f'{tmppfx}.2.sto', 'w') as f:  f.write(testmsa2)
with open(f'{tmppfx}.mask',  'w') as f:  f.write(testmask)




# basic. file against itself
r = esl_itest.run(f'{easel} compalign {srcdir}/testsuite/example-rna.sto {srcdir}/testsuite/example-rna.sto')
if (m := re.search(r'36375\s+/\s+36375\s+\(1.000\)', r.stdout, flags=re.MULTILINE)) is None: esl_itest.fail()

# --amino
# --dna
# --rna
#    All give the same result. 
r  = esl_itest.run(f'{easel} compalign --amino {tmppfx}.1.sto {tmppfx}.2.sto')
r2 = esl_itest.run(f'{easel} compalign --dna   {tmppfx}.1.sto {tmppfx}.2.sto')
r3 = esl_itest.run(f'{easel} compalign --rna   {tmppfx}.1.sto {tmppfx}.2.sto')
if r.stdout != r2.stdout: esl_itest.fail()
if r.stdout != r3.stdout: esl_itest.fail()

#   -c            : print per column statistics instead of per sequence stats
r = esl_itest.run(f'{easel} compalign -c --dna   {tmppfx}.1.sto {tmppfx}.2.sto')
if len( r.stdout.splitlines()) != alen+2: esl_itest.fail()    # <alen> output lines + 2 header lines that start with '#'

#   -p            : print stats on accuracy versus posterior probability (PP)
r = esl_itest.run(f'{easel} compalign -p --dna   {tmppfx}.1.sto {tmppfx}.2.sto')
if (m := re.search(r'^\s*9\s+22\s*/\s+27\s+\(0\.81481\)', r.stdout, flags=re.MULTILINE)) is None: esl_itest.fail()   # Checks the PP=9 output line: 22/27

#   --p-mask <f>  : with -p, only consider columns within mask ('1' columns) in <f>
r = esl_itest.run(f'{easel} compalign -p --dna --p-mask {tmppfx}.mask  {tmppfx}.1.sto {tmppfx}.2.sto')
if (m := re.search(r'^\s*9\s+9\s*/\s+9\s+\(1\.0+\)', r.stdout, flags=re.MULTILINE)) is None: esl_itest.fail()   # Checks the PP=9 output line

#   --c2dfile <f> : print per consensus column stats to esl-ssdraw --dfile file <f>
r = esl_itest.run(f'{easel} compalign -c --dna --c2dfile {tmppfx}.dfile  {tmppfx}.1.sto {tmppfx}.2.sto')
if (m := re.search(r'^# Draw file of per-column stats saved to file', r.stdout, flags=re.MULTILINE)) is None: esl_itest.fail()   
nlines = sum(1 for _ in open(f'{tmppfx}.dfile'))
if nlines != (rflen+1)*2: esl_itest.fail()    


# Cleanup
for tmpfile in glob.glob(f'{tmppfx}.*'):
    os.remove(tmpfile)

print('ok')
