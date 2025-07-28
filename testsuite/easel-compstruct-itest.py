#! /usr/bin/env python3

# Integration test for `easel compstruct` 
#
# Usage: easel-compstruct-itest.py <builddir> <srcdir> <tmppfx>
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
files_used = [ ]

(builddir, srcdir, tmppfx) = esl_itest.getargs(sys.argv)
esl_itest.check_files(srcdir,   files_used)
esl_itest.check_progs(builddir, progs_used)

easel = f'{builddir}/miniapps/easel'

# Test files
# seq2 structure slips by -1; counts as correct under Mathews' definition.
# seq1 is underpredicted
# seq3 is overpredicted
#
trusted_msa = """\
# STOCKHOLM 1.0

seq1         GGGGAAACCCCUUU
#=GR seq1 SS <<<<.AA>>>>.aa
seq2         GGGGAAACCCCUUU
#=GR seq2 SS <<<<.AA>>>>.aa
seq3         GGGGAAACCCCUUU
#=GR seq3 SS <<<<.AA>>>>.aa
//
"""

test_msa = """\
# STOCKHOLM 1.0

seq1         GGGGAAACCCCUUU
#=GR seq1 SS <<<..A..>>>.a.
seq2         GGGGAAACCCCUUU
#=GR seq2 SS <<<..AA>>>..aa
seq3         GGGGAAACCCCUUU
#=GR seq3 SS <<<<<.>>>>>...
//
"""

with open(f'{tmppfx}.trust', 'w') as f: f.write(trusted_msa)
with open(f'{tmppfx}.test',  'w') as f: f.write(test_msa)



# `-h` help 
r = esl_itest.run(f'{easel} compstruct -h')

# basic
r  = esl_itest.run(f'{easel} compstruct {tmppfx}.trust {tmppfx}.test')
if re.search(r'^seq1\s*==\s+3\s+4\s+75\.00%\s+3\s+3\s+100\.00%', r.stdout, flags=re.MULTILINE) is None: esl_itest.fail()
if re.search(r'^seq2\s*==\s+0\s+4\s+0\.00%\s+0\s+3\s+0\.00%',    r.stdout, flags=re.MULTILINE) is None: esl_itest.fail()
if re.search(r'^seq3\s*==\s+4\s+4\s+100\.00%\s+4\s+5\s+80\.00%', r.stdout, flags=re.MULTILINE) is None: esl_itest.fail()
if re.search(r'58\.33% sensitivity',                             r.stdout)                     is None: esl_itest.fail()
if re.search(r'63\.64% PPV',                                     r.stdout)                     is None: esl_itest.fail()


# -m    use Mathews' relaxed criterion, allowing +1/-1 slip
#
r  = esl_itest.run(f'{easel} compstruct -m {tmppfx}.trust {tmppfx}.test')
if re.search(r'^seq1\s*==\s+3\s+4\s+75\.00%\s+3\s+3\s+100\.00%',  r.stdout, flags=re.MULTILINE) is None: esl_itest.fail()
if re.search(r'^seq2\s*==\s+4\s+4\s+100\.00%\s+3\s+3\s+100\.00%', r.stdout, flags=re.MULTILINE) is None: esl_itest.fail()
if re.search(r'^seq3\s*==\s+4\s+4\s+100\.00%\s+4\s+5\s+80\.00%',  r.stdout, flags=re.MULTILINE) is None: esl_itest.fail()
if re.search(r'91\.67% sensitivity',                              r.stdout)                     is None: esl_itest.fail()
if re.search(r'90\.91% PPV',                                      r.stdout)                     is None: esl_itest.fail()

# -p    count pseudoknot annotation too
r  = esl_itest.run(f'{easel} compstruct -p {tmppfx}.trust {tmppfx}.test')
if re.search(r'^seq1\s*==\s+3\s+6\s+50\.00%\s+3\s+4\s+75\.00%', r.stdout, flags=re.MULTILINE) is None: esl_itest.fail()
if re.search(r'^seq2\s*==\s+2\s+6\s+33\.33%\s+2\s+5\s+40\.00%', r.stdout, flags=re.MULTILINE) is None: esl_itest.fail()
if re.search(r'^seq3\s*==\s+4\s+6\s+66\.67%\s+4\s+5\s+80.00%',  r.stdout, flags=re.MULTILINE) is None: esl_itest.fail()
if re.search(r'50\.00% sensitivity',                            r.stdout)                     is None: esl_itest.fail()
if re.search(r'64\.29% PPV',                                    r.stdout)                     is None: esl_itest.fail()

# Cleanup
for tmpfile in glob.glob(f'{tmppfx}.*'):
    os.remove(tmpfile)

print('ok')

