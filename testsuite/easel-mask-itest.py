#! /usr/bin/env python3

# Integration test for `easel mask` 
#
# Usage: easel-mask-itest.py <builddir> <srcdir> <tmppfx>
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
files_used = [ ]

(builddir, srcdir, tmppfx) = esl_itest.getargs(sys.argv)
esl_itest.check_files(srcdir,   files_used)
esl_itest.check_progs(builddir, progs_used)

easel = f'{builddir}/miniapps/easel'

# Create test files.
#
testseqs = """\
>seq1
aaAAAAAAAABBBBBBBBbbCCCCCCCCcc
>seq2
ddDDD
EEEee
FFFff
"""

testmask = """\
seq1  11  20
seq2   6  10
"""

with open(f'{tmppfx}.fa',   'w') as f: f.write(testseqs)
with open(f'{tmppfx}.mask', 'w') as f: f.write(testmask)


# `-h` help 
r = esl_itest.run(f'{easel} mask -h')

# basic
r = esl_itest.run(f'{easel} mask {tmppfx}.fa {tmppfx}.mask')
if (m := re.search(r'aaAAAAAAAAXXXXXXXXXXCCCCCCCCcc', r.stdout)) is None: esl_itest.fail()
if (m := re.search(r'ddDDDXXXXXFFFff',                r.stdout)) is None: esl_itest.fail()

# -o
r2 = esl_itest.run(f'{easel} mask -o {tmppfx}.out {tmppfx}.fa {tmppfx}.mask')
with open(f'{tmppfx}.out') as f: s = f.read()
if s != r.stdout: esl_itest.fail()

# --informat  specify input file format
r2 = esl_itest.run(f'{easel} mask --informat fasta {tmppfx}.fa {tmppfx}.mask')
if r2.stdout != r.stdout: esl_itest.fail()

# -R   random access: fetch seqs from ssi-indexed <sqfile>
r2 = subprocess.run(f'{easel} sindex {tmppfx}.fa'.split(), check=True, encoding='utf-8', capture_output=True)
r2 = esl_itest.run(f'{easel} mask -R {tmppfx}.fa {tmppfx}.mask')
if r2.stdout != r.stdout: esl_itest.fail()

# -r   reverse mask
r = esl_itest.run(f'{easel} mask -r {tmppfx}.fa {tmppfx}.mask')
if (m := re.search(r'XXXXXXXXXXBBBBBBBBbbXXXXXXXXXX', r.stdout)) is None: esl_itest.fail()
if (m := re.search(r'XXXXXEEEeeXXXXX',                r.stdout)) is None: esl_itest.fail()

# -l   convert masked residues to lower case
r = esl_itest.run(f'{easel} mask -l {tmppfx}.fa {tmppfx}.mask')
if (m := re.search(r'AAAAAAAAAAbbbbbbbbbbCCCCCCCCC', r.stdout)) is None: esl_itest.fail()
if (m := re.search(r'DDDDDeeeeeFFFFF',               r.stdout)) is None: esl_itest.fail()

# -m   convert masked residues to character <c> 
r = esl_itest.run(f'{easel} mask -m N {tmppfx}.fa {tmppfx}.mask')
if (m := re.search(r'aaAAAAAAAANNNNNNNNNNCCCCCCCCcc', r.stdout)) is None: esl_itest.fail()
if (m := re.search(r'ddDDDNNNNNFFFff',                r.stdout)) is None: esl_itest.fail()

# -x   mask additional <n> residues beyond <start>,<end>
r = esl_itest.run(f'{easel} mask -x 7 {tmppfx}.fa {tmppfx}.mask')
if (m := re.search(r'aaAXXXXXXXXXXXXXXXXXXXXXXXXCcc', r.stdout)) is None: esl_itest.fail()
if (m := re.search(r'XXXXXXXXXXXXXXX',                r.stdout)) is None: esl_itest.fail()



# Cleanup
for tmpfile in glob.glob(f'{tmppfx}.*'):
    os.remove(tmpfile)

print('ok')

