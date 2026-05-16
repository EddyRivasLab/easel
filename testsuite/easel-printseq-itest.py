#! /usr/bin/env python3

# Integration test for `easel printseq`
#
# Usage: easel-printseq-itest.py <builddir> <srcdir> <tmppfx>
#   <builddir>: path to Easel build dir. `easel` miniapp is <builddir>/miniapps/easel
#   <srcdir>:   path to Easel src dir.
#   <tmppfx>:   prefix we're allowed to use to create tmp files in current working dir.
#
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

# Create a small DNA test file.
#
testseq = """\
>testseq
ATTATTATCCATTACACGAAGATTTAGTATCATCTTTTTCTTAACTCTCC
ACTTGCAAAGTTGCAGCCTTTGCGGTCGCTTCATCAATATTGCCTACCAA
ATAAAAAGCTTGTTCAGGAA
"""
with open(f'{tmppfx}.fa', 'w') as f: f.write(testseq)

# Same sequence as RNA, for autodetection / --rna / --dna tests.
#
testseq_rna = """\
>testseq
AUUAUUAUCCAUUACACGAAGAUUUAGUAUCAUCUUUUUCUUAACUCUCC
ACUUGCAAAGUUGCAGCCUUUGCGGUCGCUUCAUCAAUAUUGCCUACCAA
AUAAAAAGCUUGUUCAGGAA
"""
with open(f'{tmppfx}.rna.fa', 'w') as f: f.write(testseq_rna)

# A protein sequence, for testing that printseq rejects amino acid input.
#
testseq_aa = """\
>protein1
MKLAVLDQERWQNDPYSFGRTHGIACDEFKLMNPQSTYVW
"""
with open(f'{tmppfx}.aa.fa', 'w') as f: f.write(testseq_aa)

# A two-sequence file, for testing that printseq refuses multi-sequence input.
#
testseq_multi = """\
>seq1
ACGTACGTACGT
>seq2
TTTTAAAACCCC
"""
with open(f'{tmppfx}.multi.fa', 'w') as f: f.write(testseq_multi)

# Sequence with IUPAC ambiguity codes, for testing the complement table.
# Complements: A<->T, C<->G, N->N, R<->Y.
#
testseq_ambig = """\
>ambig
ACGTNNNNYRACGT
"""
with open(f'{tmppfx}.ambig.fa', 'w') as f: f.write(testseq_ambig)

# 11-nt sequence (3 codons + 2 trailing bases) for testing -3 with a
# sequence length that isn't a multiple of 3. The trailing partial codon
# must be printed as bases but not translated.
#
testseq_partial = """\
>partial
ATGAAACCCGT
"""
with open(f'{tmppfx}.partial.fa', 'w') as f: f.write(testseq_partial)


# `-h` help
r = esl_itest.run(f'{easel} printseq -h')

# basic
r = esl_itest.run(f'{easel} printseq {tmppfx}.fa')
if re.search(r'^\s+1: ATTATTATCC ATTACACGAA .* ACTTGCAAAG 60\s*$', r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()
if re.search(r'^\s+: TAATAATAGG TAATGTGCTT',                       r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()
if re.search(r'^\s+61: TTGCAGCCTT .* TGTTCAGGAA 120\s*$',          r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()

# -s : single-stranded
r = esl_itest.run(f'{easel} printseq -s {tmppfx}.fa')
if re.search(r'^\s+1: ATTATTATCC',  r.stdout, flags=re.MULTILINE)  == None: esl_itest.fail()
if re.search(r'^\s+: TAATAATAGG',   r.stdout, flags=re.MULTILINE)  != None: esl_itest.fail()

# -a : all three top-strand reading frames, single-letter code, no grouping.
# Verify both line blocks (the 120-nt seq spans two lines at -n 60).
r = esl_itest.run(f'{easel} printseq -a {tmppfx}.fa')
if re.search(r'^\s+1: ATTATTATCCATTACACGAAGATTTAGTATCATCTTTTTCTTAACTCTCCACTTGCAAAG 60', r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()
if re.search(r'^\s+0 : I  I  I  H  Y  T  K  I',  r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()
if re.search(r'^\s+1 :  L  L  S  I',             r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()
if re.search(r'^\s+2 :   Y  Y  P  L',            r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()
if re.search(r'^\s+61: TTGCAGCCTTTGCGGTCGCTTCATCAATATTGCCTACCAAATAAAAAGCTTGTTCAGGAA 120', r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()
if re.search(r'^\s+0 : L  Q  P  L',  r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()

# -3 : groups of 3, frame-0 translation in three-letter code.
r = esl_itest.run(f'{easel} printseq -3 {tmppfx}.fa')
if re.search(r'^\s+1: ATT ATT ATC CAT TAC ACG AAG',   r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()
if re.search(r'Ile Ile Ile His Tyr Thr Lys',          r.stdout)                     == None: esl_itest.fail()

# -3 on a sequence whose length is not a multiple of 3: the trailing
# partial codon ("GT") is printed as bases but not translated.
r = esl_itest.run(f'{easel} printseq -3 {tmppfx}.partial.fa')
if re.search(r'^\s+1: ATG AAA CCC GT',  r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()
if re.search(r'Met Lys Pro\s*$',        r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()

# -f <n> : number first residue as <n>, not 1.
r = esl_itest.run(f'{easel} printseq -f 100 {tmppfx}.fa')
if re.search(r'^\s*100: ATTATTATCC',                  r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()
if re.search(r'\s159\s',                              r.stdout)                     == None: esl_itest.fail()

# -g 0 : no grouping 
r = esl_itest.run(f'{easel} printseq -g 0 {tmppfx}.fa')
if re.search(r'^\s+1: ATTATTATCCATTACACGAAGATTTAGTATCATCTTTTTCTTAACTCTCCACTTGCAAAG 60', r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()

# -g 20 : larger groups
r = esl_itest.run(f'{easel} printseq -g 20 {tmppfx}.fa')
if re.search(r'^\s+1: ATTATTATCCATTACACGAA GATTTAGTATCATCTTTTTC TTAACTCTCCACTTGCAAAG', r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()

# -n <n> : line length.
r = esl_itest.run(f'{easel} printseq -n 30 {tmppfx}.fa')
for coord in ('1', '31', '61', '91'):
    if re.search(rf'^\s*{coord}: ', r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()

# -t : coords at top
r = esl_itest.run(f'{easel} printseq -t {tmppfx}.fa')
if re.search(r'^\s+10\s+20\s+30\s+40\s+50\s+60\s*$', r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()
if re.search(r'^ATTATTATCC ATTACACGAA',              r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()

# -3 -t : top coordinates with triplet grouping (coords mark codon ends).
r = esl_itest.run(f'{easel} printseq -3 -t {tmppfx}.fa')
if re.search(r'^\s+3\s+6\s+9\s+12\s+15', r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()
if re.search(r'^ATT ATT ATC CAT TAC',    r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()
if re.search(r'Ile Ile Ile His Tyr',     r.stdout)                     == None: esl_itest.fail()

# -f with -a : -f offset is applied to the left-margin coordinate.
r = esl_itest.run(f'{easel} printseq -f 100 -a {tmppfx}.fa')
if re.search(r'^\s*100: ATTATTATCCATTACACGAAGATTTAGTATCATCTTTTTCTTAACTCTCCACTTGCAAAG 159', r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()
if re.search(r'^\s+0 : I  I  I  H  Y',   r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()

# Ambiguous bases: complement table handles N (->N), R (->Y), Y (->R).
r = esl_itest.run(f'{easel} printseq --dna {tmppfx}.ambig.fa')
if re.search(r'^\s+1: ACGTNNNNYR ACGT 14', r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()
if re.search(r'^\s+: TGCANNNNRY TGCA',     r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()

# RNA input
r = esl_itest.run(f'{easel} printseq {tmppfx}.rna.fa')
if re.search(r'^\s+1: AUUAUUAUCC AUUACACGAA',  r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()
if re.search(r'^\s+: UAAUAAUAGG UAAUGUGCUU',   r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()

# --rna 
r = esl_itest.run(f'{easel} printseq --rna {tmppfx}.fa')
if re.search(r'^\s+1: AUUAUUAUCC AUUACACGAA',  r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()

# --dna 
r = esl_itest.run(f'{easel} printseq --dna {tmppfx}.rna.fa')
if re.search(r'^\s+1: ATTATTATCC ATTACACGAA',  r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()

# --informat
r2 = esl_itest.run(f'{easel} printseq --informat fasta {tmppfx}.fa')
r  = esl_itest.run(f'{easel} printseq                  {tmppfx}.fa')
if r.stdout != r2.stdout: esl_itest.fail()
r = esl_itest.run(f'{easel} printseq --informat bogus {tmppfx}.fa', expect_success=False)
r = esl_itest.run(f'{easel} printseq --informat genbank {tmppfx}.fa', expect_success=False)

# Protein input is rejected
r = esl_itest.run(f'{easel} printseq {tmppfx}.aa.fa', expect_success=False)
if re.search(r'protein sequence', r.stderr) == None: esl_itest.fail()

# Multi-sequence input: printseq formats the first sequence and warns to stderr.
r = esl_itest.run(f'{easel} printseq {tmppfx}.multi.fa')
if re.search(r'^\s+1: ACGTACGTAC GT', r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()
if re.search(r'more than one sequence', r.stderr)                   == None: esl_itest.fail()

# Option-conflict tests. printseq enforces:
#   -a and -3 are mutually exclusive
#   -a requires -g 0    (an explicit nonzero -g with -a is fatal)
#   -3 requires -g 3    (an explicit -g != 3 with -3 is fatal)
#   -t requires nonzero -g
#   --dna and --rna are a toggle group (mutually exclusive)
#
r = esl_itest.run(f'{easel} printseq -a -3      {tmppfx}.fa', expect_success=False)
r = esl_itest.run(f'{easel} printseq -a -g 5    {tmppfx}.fa', expect_success=False)
r = esl_itest.run(f'{easel} printseq -3 -g 5    {tmppfx}.fa', expect_success=False)
r = esl_itest.run(f'{easel} printseq -t -g 0    {tmppfx}.fa', expect_success=False)
r = esl_itest.run(f'{easel} printseq --dna --rna {tmppfx}.fa', expect_success=False)


# Cleanup
for tmpfile in glob.glob(f'{tmppfx}.*'):
    os.remove(tmpfile)

print('ok')
