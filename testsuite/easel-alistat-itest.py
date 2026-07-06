#! /usr/bin/env python3

# Integration test for `easel alistat` 
#
# Usage: easel-alistat-itest.py <builddir> <srcdir> <tmppfx>
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
files_used = [ 'testsuite/example-rfam.sto',
               'testsuite/example-rna.sto',
               'testsuite/example-stockholm.sto' ]

(builddir, srcdir, tmppfx) = esl_itest.getargs(sys.argv)
esl_itest.check_files(srcdir,   files_used)
esl_itest.check_progs(builddir, progs_used)

easel = f'{builddir}/miniapps/easel'

# -h
r = esl_itest.run(f'{easel} alistat -h')

# -1             tabular output
r = esl_itest.run(f'{easel} alistat -1 {srcdir}/testsuite/example-rfam.sto')
if ( m := re.search(r'3\s+DUF3800-I\s+RF02969', r.stdout, flags=re.MULTILINE)) is None: esl_itest.fail()

# -q             (with -1) quieter; suppress tabular header 
r = esl_itest.run(f'{easel} alistat -1q {srcdir}/testsuite/example-rfam.sto')
if ( m := re.search(r'^3\s+DUF3800-I\s+RF02969', r.stdout, flags=re.MULTILINE)) is     None: esl_itest.fail()     
if ( m := re.search(r'^#',                       r.stdout, flags=re.MULTILINE)) is not None: esl_itest.fail()     

r = esl_itest.run(f'{easel} alistat -q {srcdir}/testsuite/example-rfam.sto', expect_success=False)

# --recsize      (with -1) include MSA record size (bytes) in tabular output
r = esl_itest.run(f'{easel} alistat -1 --recsize {srcdir}/testsuite/example-rfam.sto')
if ( m := re.search(r'^#.+recsize\s+size/nres', r.stdout, flags=re.MULTILINE)) is None: esl_itest.fail()     

r = esl_itest.run(f'{easel} alistat --recsize {srcdir}/testsuite/example-rfam.sto', expect_success=False)


# --informat     assert <msafile> format or alphabet
# --amino       
# --dna         
# --rna         
r  = esl_itest.run(f'{easel} alistat --informat stockholm {srcdir}/testsuite/example-stockholm.sto')
r2 = esl_itest.run(f'{easel} alistat                      {srcdir}/testsuite/example-stockholm.sto')
if r2.stdout != r.stdout: esl_itest.fail()

r  = esl_itest.run(f'{easel} alistat --amino              {srcdir}/testsuite/example-stockholm.sto')
r2 = esl_itest.run(f'{easel} alistat                      {srcdir}/testsuite/example-stockholm.sto')
if r2.stdout != r.stdout: esl_itest.fail()

r  = esl_itest.run      (f'{easel} alistat --dna         {srcdir}/testsuite/example-rna.sto')
r2 = esl_itest.run_piped(f'{easel} reformat -d stockholm {srcdir}/testsuite/example-rna.sto', f'{easel} alistat -')
if r2.stdout != r.stdout: esl_itest.fail()

r  = esl_itest.run(f'{easel} alistat --rna                {srcdir}/testsuite/example-rna.sto')
r2 = esl_itest.run(f'{easel} alistat                      {srcdir}/testsuite/example-rna.sto')
if r2.stdout != r.stdout: esl_itest.fail()


# Cleanup
for tmpfile in glob.glob(f'{tmppfx}.*'):
    os.remove(tmpfile)


print('ok')

