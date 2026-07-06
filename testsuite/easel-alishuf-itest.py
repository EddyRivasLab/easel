#! /usr/bin/env python3

# Integration test for `easel alishuf` 
#
# Usage: easel-alishuf-itest.py <builddir> <srcdir> <tmppfx>
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

def parse_alistat(output):
    if (m := re.search(r'^Number of sequences:\s*(\d+)', output, flags=re.MULTILINE)) is None: esl_itest.fail()
    nseq = int(m.group(1))
    if (m := re.search(r'^Alignment length:\s*(\d+)',    output, flags=re.MULTILINE)) is None: esl_itest.fail()
    alen = int(m.group(1))
    return (nseq, alen)


r = esl_itest.run(f'{easel} alistat {srcdir}/testsuite/example-rna.sto')
expected_nseq, expected_alen = parse_alistat(r.stdout)

# `-h` help 
r = esl_itest.run(f'{easel} alishuf -h')

# basic
r  = esl_itest.run(f'{easel} alishuf {srcdir}/testsuite/example-rna.sto')
r2 = subprocess.run(f'{easel} alistat -'.split(), check=True, encoding='utf-8', capture_output=True, input=r.stdout)
nseq, alen = parse_alistat(r2.stdout)
if nseq != expected_nseq or alen != expected_alen: esl_itest.fail()

#  -o           direct output data to file <f>
# --seed <n>    set random number generator seed to <n>  
#
# By using --seed, we can compare two shuffles for exact equality.
#
r   = esl_itest.run(f'{easel} alishuf --seed 42                 {srcdir}/testsuite/example-rna.sto')
r2  = esl_itest.run(f'{easel} alishuf --seed 42 -o {tmppfx}.sto {srcdir}/testsuite/example-rna.sto')
with open(f'{tmppfx}.sto') as f: s = f.read()
if s != r.stdout: esl_itest.fail()

# --rna         assert <msafile> is RNA
#
r2  = esl_itest.run(f'{easel} alishuf --seed 42 --rna {srcdir}/testsuite/example-rna.sto')
if r2.stdout != r.stdout: esl_itest.fail()

# --informat    assert <msafile> is in format <s>
r2  = esl_itest.run(f'{easel} alishuf --seed 42 --informat stockholm {srcdir}/testsuite/example-rna.sto')
if r2.stdout != r.stdout: esl_itest.fail()

# --dna         assert <msafile> is DNA 
#
# Now the shuffle will differ, using DNA T instead of RNA U, so don't compare exactly.
#
r2 = esl_itest.run(f'{easel} alishuf --dna {srcdir}/testsuite/example-rna.sto')
r3 = subprocess.run(f'{easel} alistat -'.split(), check=True, encoding='utf-8', capture_output=True, input=r2.stdout)
nseq, alen = parse_alistat(r3.stdout)
if nseq != expected_nseq or alen != expected_alen: esl_itest.fail()

#  -v     shuffle residues in each column independently
#
r  = esl_itest.run(f'{easel} alishuf -v {srcdir}/testsuite/example-rna.sto')
r2 = subprocess.run(f'{easel} alistat -'.split(), check=True, encoding='utf-8', capture_output=True, input=r.stdout)
nseq, alen = parse_alistat(r2.stdout)
if nseq != expected_nseq or alen != expected_alen: esl_itest.fail()

#  -b    take bootstrapping samples
#  -N    generate <n> samples per input msa (e.g. bootstraps)
#
r  = esl_itest.run(f'{easel} alishuf -b -N 10 {srcdir}/testsuite/example-rna.sto')
r2 = subprocess.run(f'{easel} alistat -1 -'.split(), check=True, encoding='utf-8', capture_output=True, input=r.stdout)
for line in r2.stdout.splitlines():
    if line.startswith('#'): continue
    if (m := re.match(r'(\d+)\s+(\S+)\s+-\s+Stockholm\s+(\d+)\s+(\d+)', line)) is None: esl_itest.fail()
    n       = int(m.group(1))
    msaname = m.group(2)
    nseq    = int(m.group(3))
    alen    = int(m.group(4))
    if not msaname.endswith(f'-bootsample-{n}'): esl_itest.fail()
    if nseq != expected_nseq or alen != expected_alen: esl_itest.fail()
    

# --amino       assert <msafile> is protein
#
# Switch to using a protein MSA; recalculate the expected nseq,alen.
#
r = esl_itest.run(f'{easel} alistat {srcdir}/testsuite/example-stockholm.sto')
expected_nseq, expected_alen = parse_alistat(r.stdout)

r  = esl_itest.run(f'{easel} alishuf --amino {srcdir}/testsuite/example-stockholm.sto')
r2 = subprocess.run(f'{easel} alistat -'.split(), check=True, encoding='utf-8', capture_output=True, input=r.stdout)
nseq, alen = parse_alistat(r2.stdout)
if nseq != expected_nseq or alen != expected_alen: esl_itest.fail()



# Cleanup
for tmpfile in glob.glob(f'{tmppfx}.*'):
    os.remove(tmpfile)

print('ok')
