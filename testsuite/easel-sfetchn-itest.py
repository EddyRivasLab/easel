#! /usr/bin/env python3

# Integration test for `easel sfetchn` 
#
# Usage: easel-sfetchn-itest.py <builddir> <srcdir> <tmppfx>
#   <builddir>: path to Easel build dir. `easel` miniapp is <builddir>/miniapps/easel
#   <srcdir>:   path to Easel src dir.
#   <tmppfx>:   prefix we're allowed to use to create tmp files in current working dir.
#
import filecmp
import glob
import os
import re
import shutil
import sys
import esl_itest

files_used = [ 'testsuite/example-genbank.gb',     # 4 phage DNA seqs:    NC_047788 NC_055916 NC_007046 NC_049972, accessions same as names
               'testsuite/example-uniprot.dat',    # 4 protein seqs:      MNME_BEII9 DEF_RICCK GPMI_YERP3 FABZ_PROM2; accessions B2IJQ3 A8EXV2 A7FCU8 A8G6E7
               'testsuite/example-uniprot.fa' ]    # same 4 seqs as .dat: sp|B2IJQ3|MNME_BEII9, sp|A8EXV2|DEF_RICCK, sp|A7FCU8|GPMI_YERP3, sp|A8G6E7|FABZ_PROM2

progs_used = [ 'miniapps/easel' ]


(builddir, srcdir, tmppfx) = esl_itest.getargs(sys.argv)
esl_itest.check_files(srcdir,   files_used)
esl_itest.check_progs(builddir, progs_used)

easel = f'{builddir}/miniapps/easel'


# -h
r = esl_itest.run(f'{easel} sfetchn -h')

# Make copies of three example files from Easel testsuite directory.
shutil.copyfile(f'{srcdir}/testsuite/example-genbank.gb',  f'{tmppfx}.gb')
shutil.copyfile(f'{srcdir}/testsuite/example-uniprot.dat', f'{tmppfx}.dat')
shutil.copyfile(f'{srcdir}/testsuite/example-uniprot.fa',  f'{tmppfx}.fa')

# Index them. (SSI is mandatory for sfetchn)
r = esl_itest.run(f'{easel} sindex --noacc {tmppfx}.gb')    # GenBank name = accession, so don't bother indexing accession
r = esl_itest.run(f'{easel} sindex {tmppfx}.dat')
r = esl_itest.run(f'{easel} sindex -u {tmppfx}.fa')

## Complete sequence fetching
##

# Uniprot .fa can be fetched using <id> or <acc> in <db>|<acc>|<id> names, when sindexed with -u
#
with open(f'{tmppfx}.list', 'w') as f:
    f.write('# comment\n\nGPMI_YERP3 ignore\n\nMNME_BEII9 other fields\n')
r = esl_itest.run(f'{easel} sfetchn {tmppfx}.fa {tmppfx}.list')
if re.search(r'^>sp\|A7FCU8\|GPMI_YERP3 (?s:.+)^>sp\|B2IJQ3\|MNME_BEII9 ', r.stdout,  flags=re.MULTILINE) is None: esl_itest.fail()  # indexed fetch: seqs in order of <keyfile>

with open(f'{tmppfx}.list', 'w') as f:
    f.write('# comment\n\nA7FCU8\n\nB2IJQ3\n')
r = esl_itest.run(f'{easel} sfetchn {tmppfx}.fa {tmppfx}.list')
if re.search(r'^>sp\|A7FCU8\|GPMI_YERP3 (?s:.+)^>sp\|B2IJQ3\|MNME_BEII9 ', r.stdout,  flags=re.MULTILINE) is None: esl_itest.fail()  

with open(f'{tmppfx}.list', 'w') as f:
    f.write('sp|A7FCU8|GPMI_YERP3\nsp|B2IJQ3|MNME_BEII9')
r = esl_itest.run(f'{easel} sfetchn {tmppfx}.fa {tmppfx}.list')
if re.search(r'^>sp\|A7FCU8\|GPMI_YERP3 (?s:.+)^>sp\|B2IJQ3\|MNME_BEII9 ', r.stdout,  flags=re.MULTILINE) is None: esl_itest.fail()  

# Fetching complete, nonrevcomp sequences is verbatim in original format. This also tests -o
#
with open(f'{tmppfx}.list', 'w') as f:
    f.write('NC_047788\nNC_055916\nNC_007046\nNC_049972\n')
r = esl_itest.run(f'{easel} sfetchn -o {tmppfx}.out {tmppfx}.gb {tmppfx}.list')
if filecmp.cmp(f'{tmppfx}.gb', f'{tmppfx}.out', shallow=False) == False: esl_itest.fail()

# -o will refuse to overwrite: this fails:
r = esl_itest.run(f'{easel} sfetchn -o {tmppfx}.out {tmppfx}.gb {tmppfx}.list', expect_success=False)

# -f allows the overwrite
r = esl_itest.run(f'{easel} sfetchn -fo {tmppfx}.out {tmppfx}.gb {tmppfx}.list')

# --informat
with open(f'{tmppfx}.list', 'w') as f:
    f.write('DEF_RICCK\nGPMI_YERP3\n')
r = esl_itest.run(f'{easel} sfetchn -fo {tmppfx}.out --informat uniprot {tmppfx}.dat {tmppfx}.list')

# -r  : reverse complement comes out in FASTA, not original format
with open(f'{tmppfx}.list', 'w') as f:
    f.write('NC_049972\nNC_007046\n')
r = esl_itest.run(f'{easel} sfetchn -r -fo {tmppfx}.out {tmppfx}.gb {tmppfx}.list')
r = esl_itest.run(f'{easel} seqstat {tmppfx}.out')
if re.search(r'^Format:\s+FASTA(?s:.+)^Total # residues:\s+37131', r.stdout, flags=re.MULTILINE) is None: esl_itest.fail()


## Subsequence fetching (-C)
##

# end=0 means fetch suffix
with open(f'{tmppfx}.list', 'w') as f:
    f.write('subseq1 101 200 NC_049972\nsubseq2 18000 0 NC_007046\n')
r = esl_itest.run(f'{easel} sfetchn -C -fo {tmppfx}.out {tmppfx}.gb {tmppfx}.list')
r = esl_itest.run(f'{easel} seqstat -Aq {tmppfx}.out')                             # .out: 100, 199nt subseqs
if re.search(r'^subseq1\s+100 (?s:.+)^subseq2\s+200 ', r.stdout, flags=re.MULTILINE) is None: esl_itest.fail()

# revcomp by -r
r = esl_itest.run(f'{easel} sfetchn -C -r -fo {tmppfx}.out2 {tmppfx}.gb {tmppfx}.list')  # .out2: revcomp of the 2 subseqs
with open(f'{tmppfx}.list2', 'w') as f:                                                  # .list2: for complete sfetchn of the 2 subseqs
    f.write('subseq1\nsubseq2\n')
r = esl_itest.run(f'{easel} sindex -f {tmppfx}.out2')
r = esl_itest.run(f'{easel} sfetchn -r -fo {tmppfx}.out3 {tmppfx}.out2 {tmppfx}.list2')      # .out3 now == .out: revcomp of revcomp
if filecmp.cmp(f'{tmppfx}.out3', f'{tmppfx}.out', shallow=False) == False: esl_itest.fail()

# revcomp by coord
with open(f'{tmppfx}.list', 'w') as f:
    f.write('subseq1 200 101 NC_049972\nsubseq2 18199 18000 NC_007046\n')
r = esl_itest.run(f'{easel} sfetchn -C -fo {tmppfx}.out2 {tmppfx}.gb {tmppfx}.list')        # again .out2: revcomp of the two subseqs
r = esl_itest.run(f'{easel} sindex -f {tmppfx}.out2')
r = esl_itest.run(f'{easel} sfetchn -r -fo {tmppfx}.out3 {tmppfx}.out2 {tmppfx}.list2')     # .out3 now == .out: revcomp of revcomp
if filecmp.cmp(f'{tmppfx}.out3', f'{tmppfx}.out', shallow=False) == False: esl_itest.fail()

for tmpfile in glob.glob(f'{tmppfx}.*'):
    os.remove(tmpfile)

print('ok')
