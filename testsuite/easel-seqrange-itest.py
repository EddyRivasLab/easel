#! /usr/bin/env python3

# Integration test for `easel seqrange` 
#
# Usage: easel-seqrange-itest.py <builddir> <srcdir> <tmppfx>
#   <builddir>: path to Easel build dir. `easel` miniapp is <builddir>/miniapps/easel
#   <srcdir>:   path to Easel src dir.
#   <tmppfx>:   prefix we're allowed to use to create tmp files in current working dir.
#
import glob
import os
import re
import sys
import esl_itest

progs_used = [ 'miniapps/easel' ]
files_used = [ ]

(builddir, srcdir, tmppfx) = esl_itest.getargs(sys.argv)
esl_itest.check_files(srcdir,   files_used)
esl_itest.check_progs(builddir, progs_used)

easel = f'{builddir}/miniapps/easel'

# Create test seqfile with 31 (short) sequences, then SSI index it.
testseqs = """\
>random0
CUGCUUCGCA
>random1
GAGUACGUGG
>random2
GCUACCCUAA
>random3
GGUAACCUAA
>random4
AUUAGGGCAU
>random5
CCGACUUUAG
>random6
ACCAUUUACA
>random7
GACUAGAAAC
>random8
AUGUAGAGUA
>random9
CGCAGCCGGC
>random10
CAGAACUUCG
>random11
GAGGUCAGGC
>random12
UCACUUGUCG
>random13
ACCGGGGAUG
>random14
CGAUUUUCGG
>random15
CUGGUCCUGG
>random16
AUGUGAAGAC
>random17
AAUGAAGGUU
>random18
UGCUCCGGCG
>random19
CGCGACAUGG
>random20
AAAGCGACCG
>random21
UCCUGUAAGC
>random22
CGUUUCAUGG
>random23
GCAAAACGGC
>random24
GUCGCAUAUU
>random25
GGCUAACAUC
>random26
CUUGGUCUGC
>random27
CCCAGAGUGU
>random28
GUGUGCGCGU
>random29
ACGGCACCAA
>random30
GGGCAGUGCG
"""

with open(f'{tmppfx}.fa', 'w') as f: f.write(testseqs)
r = esl_itest.run(f'{easel} sindex {tmppfx}.fa')

# `-h` help 
r = esl_itest.run(f'{easel} seqrange -h')

# Eric's tests, from previous Perl integrated test
r = esl_itest.run(f'{easel} seqrange {tmppfx}.fa 1 31')
if (m := re.fullmatch(r'1-1\s*', r.stdout)) is None: esl_itest.fail()

r = esl_itest.run(f'{easel} seqrange {tmppfx}.fa 17 31')
if (m := re.fullmatch(r'17-17\s*', r.stdout)) is None: esl_itest.fail()

r = esl_itest.run(f'{easel} seqrange {tmppfx}.fa 1 13')
if (m := re.fullmatch(r'1-3\s*', r.stdout)) is None: esl_itest.fail()

r = esl_itest.run(f'{easel} seqrange {tmppfx}.fa 1 3')
if (m := re.fullmatch(r'1-11\s*', r.stdout)) is None: esl_itest.fail()

r = esl_itest.run(f'{easel} seqrange {tmppfx}.fa 2 3')
if (m := re.fullmatch(r'12-21\s*', r.stdout)) is None: esl_itest.fail()

r = esl_itest.run(f'{easel} seqrange {tmppfx}.fa 3 3')
if (m := re.fullmatch(r'22-31\s*', r.stdout)) is None: esl_itest.fail()

# --informat
r = esl_itest.run(f'{easel} seqrange --informat fasta {tmppfx}.fa 3 3')
if (m := re.fullmatch(r'22-31\s*', r.stdout)) is None: esl_itest.fail()

# Eric's failure cases
r = esl_itest.run(f'{easel} seqrange --informat fasta {tmppfx}.fa 1 32', expect_success=False)
r = esl_itest.run(f'{easel} seqrange --informat fasta {tmppfx}.fa 31 1', expect_success=False)
r = esl_itest.run(f'{easel} seqrange --informat fasta {tmppfx}.fa -1 1', expect_success=False)
r = esl_itest.run(f'{easel} seqrange --informat fasta {tmppfx}.fa 1 -1', expect_success=False)


# Cleanup
for tmpfile in glob.glob(f'{tmppfx}.*'):
    os.remove(tmpfile)

print('ok')

