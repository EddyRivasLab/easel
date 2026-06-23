#! /usr/bin/env python3

# Integration test for `easel alicol`
#
# Usage: easel-alicol-itest.py <builddir> <srcdir> <tmppfx>
#   <builddir>: path to Easel build dir. `easel` miniapp is <builddir>/miniapps/easel
#   <srcdir>:   path to Easel src dir.
#   <tmppfx>:   prefix we're allowed to use to create tmp files in current working dir.
#
# Also uses `easel reformat`.

import glob
import os
import re
import subprocess
import sys
import esl_itest

# seq1(r)
#
# Pull the aligned seq1 residues out of an `easel alicol` Stockholm
# (or Pfam) result <r>, and return them as a string. All our test MSAs
# have a sequence named seq1; checking what columns survive in it is
# the most compact way to verify a column selection.
#
def seq1(r):
    for line in r.stdout.splitlines():
        if (m := re.match(r'seq1\s+(\S+)\s*$', line)): return m.group(1)
    esl_itest.fail('no seq1 in output')


# ss_cons(r)
#
# Pull the #=GC SS_cons consensus structure annotation out of an
# `easel alicol` result <r>. Used to check that small mode fixes broken
# base pairs (deletes a pair when a column subset removes one partner)
# exactly the way the default path does.
#
def ss_cons(r):
    for line in r.stdout.splitlines():
        if (m := re.match(r'#=GC SS_cons\s+(\S+)\s*$', line)): return m.group(1)
    esl_itest.fail('no #=GC SS_cons in output')


# write_msa(outfile)
#
# A small test MSA with #=GC RF consensus annotation.
#
def write_msa(outfile):
    msa = """\
# STOCKHOLM 1.0
#=GF ID test
seq1 AC-DEF-G
seq2 AC-DE--G
seq3 AC-DEFMG
#=GC RF xx.xxx.x
//"""
    with open(outfile, 'w') as f:
        print(msa, file=f)


# write_pp_msa(outfile)
#
# A small test MSA with posterior probability annotation: per-sequence
# #=GR PP, consensus #=GC PP_cons, and #=GC RF. Discretized PP symbols
# 0-9* are interpreted by alicol as the min of their bin (for --ppcons,
# --ppfrac) or the mean (for --ppavg):
#
#                       min     mean
#    '0'  0.00-0.05     0.00    0.025
#    '5'  0.45-0.55     0.45    0.5
#    '7'  0.65-0.75     0.65    0.7
#    '8'  0.75-0.85     0.75    0.8
#    '9'  0.85-0.95     0.85    0.9
#    '*'  0.95-1.00     0.95    0.975
#
# Per-column unweighted PP means (over the 2 seqs), used in tests below:
#    col:    1      2      3       4      5
#    PP:    **     99     8*      57     03
#    mean: .975    .9    .8875    .6    .1625
#
def write_pp_msa(outfile, weighted=False):
    wgt = ""
    if weighted:
        wgt = "#=GS seq1 WT 3.0\n#=GS seq2 WT 1.0\n"
    msa = """\
# STOCKHOLM 1.0
#=GF ID pptest
{}seq1         ACDEF
#=GR seq1 PP *9850
seq2         ACDEF
#=GR seq2 PP *9*73
#=GC RF      xxxxx
#=GC PP_cons *9750
//""".format(wgt)
    with open(outfile, 'w') as f:
        print(msa, file=f)


# write_rna_ss_msa(outfile)
#
# A small RNA test MSA with #=GC SS_cons secondary structure annotation
# (and #=GC RF). The structure "<<..>>" pairs cols 1-6 and 2-5. Removing
# either partner of a pair must delete the whole pair from SS_cons (a
# "broken base pair"); `easel alicol` does this for RNA|DNA MSAs. Small
# mode fixes broken pairs inside esl_msafile2_RegurgitatePfam(); we test
# that it matches the default path.
#
def write_rna_ss_msa(outfile):
    msa = """\
# STOCKHOLM 1.0
#=GF ID rnatest
seq1         ACGUUA
seq2         ACGUUA
#=GC RF      xxxxxx
#=GC SS_cons <<..>>
//"""
    with open(outfile, 'w') as f:
        print(msa, file=f)



def main():
    progs_used = [ 'miniapps/easel' ]
    files_used = [ ]

    (builddir, srcdir, tmppfx) = esl_itest.getargs(sys.argv)
    esl_itest.check_files(srcdir,   files_used)
    esl_itest.check_progs(builddir, progs_used)

    easel = f'{builddir}/miniapps/easel'

    write_msa   (f'{tmppfx}.sto')
    write_pp_msa(f'{tmppfx}.pp.sto')
    write_pp_msa(f'{tmppfx}.ppw.sto', weighted=True)

    # -h   help
    r = esl_itest.run(f'{easel} alicol -h')

    # At least one selection option is required.
    r = esl_itest.run(f'{easel} alicol {tmppfx}.sto', expect_success=False)

    # Gap/residue fraction selections, and their --mingap/--nogap aliases.
    #   --mingap removes the all-gap column 3:  AC-DEF-G -> ACDEF-G
    #   --nogap  keeps only all-residue cols:   AC-DEF-G -> ACDEG
    #
    if seq1(esl_itest.run(f'{easel} alicol --mingap        {tmppfx}.sto')) != 'ACDEF-G': esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol --nogap         {tmppfx}.sto')) != 'ACDEG':   esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol --gapfrac 1.0   {tmppfx}.sto')) != 'ACDEF-G': esl_itest.fail()  # --gapfrac 1 == --mingap
    if seq1(esl_itest.run(f'{easel} alicol --symfrac 1.0   {tmppfx}.sto')) != 'ACDEG':   esl_itest.fail()  # --symfrac 1 == --nogap

    # --rfonly keeps the 6 consensus columns, removes the 2 inserts.
    if seq1(esl_itest.run(f'{easel} alicol --rfonly        {tmppfx}.sto')) != 'ACDEFG':  esl_itest.fail()

    # --keeprf: consensus columns are always kept, selection only applies to inserts.
    # --nogap alone gives ACDEG (drops consensus col 6=F); with --keeprf, col 6 is kept: ACDEFG.
    if seq1(esl_itest.run(f'{easel} alicol --nogap --keeprf {tmppfx}.sto')) != 'ACDEFG': esl_itest.fail()

    # --span, in default 1..alen coords. cols 4-6 are DEF.
    if seq1(esl_itest.run(f'{easel} alicol --span 4:6      {tmppfx}.sto')) != 'DEF':     esl_itest.fail()
    # --span suffix rule "4:" means cols 4..alen.
    if seq1(esl_itest.run(f'{easel} alicol --span 4:       {tmppfx}.sto')) != 'DEF-G':   esl_itest.fail()

    # --span in 1..C consensus coords, with -c.
    #   consensus cols 1-3 are MSA cols 1,2,4; span keeps everything between them inclusive: AC-D
    #   consensus suffix "4:" is consensus cols 4-6 = MSA cols 5,6,8, keeping inserts between: EF-G
    if seq1(esl_itest.run(f'{easel} alicol -c --span 1:3   {tmppfx}.sto')) != 'AC-D':    esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol -c --span 4:    {tmppfx}.sto')) != 'EF-G':    esl_itest.fail()

    # --span error cases: reversed coords, reverse-complement suffix, off-the-end.
    r = esl_itest.run(f'{easel} alicol --span 6:3  {tmppfx}.sto', expect_success=False)
    r = esl_itest.run(f'{easel} alicol --span :5   {tmppfx}.sto', expect_success=False)
    r = esl_itest.run(f'{easel} alicol --span 4:99 {tmppfx}.sto', expect_success=False)

    # -c requires a consensus-coord selection option (--span|--mask)...
    r = esl_itest.run(f'{easel} alicol -c --nogap  {tmppfx}.sto', expect_success=False)
    # ...and requires consensus annotation to be present in the MSA.
    with open(f'{tmppfx}.norf.sto', 'w') as f:
        print("# STOCKHOLM 1.0\nseq1 ACDEF\nseq2 ACDEG\n//", file=f)
    r = esl_itest.run(f'{easel} alicol -c --span 1:3 {tmppfx}.norf.sto', expect_success=False)

    # Posterior probability selections (need Stockholm w/ PP annotation).
    #
    # --ppcons keeps consensus cols with PP_cons min-of-bin >= x.
    #   PP_cons = *9750; mins = .95 .85 .65 .45 .0
    if seq1(esl_itest.run(f'{easel} alicol --ppcons 0.9 {tmppfx}.pp.sto')) != 'A':     esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol --ppcons 0.6 {tmppfx}.pp.sto')) != 'ACD':   esl_itest.fail()

    # --ppcons removes inserts by default; --keepins keeps them. (No inserts in
    # this all-consensus MSA, so we just check both run and select the same cols.)
    if seq1(esl_itest.run(f'{easel} alicol --ppcons 0.6 --keepins {tmppfx}.pp.sto')) != 'ACD': esl_itest.fail()

    # --ppavg keeps cols whose mean PP (mean-of-bin) >= x.
    #   col means: .975 .9 .8875 .6 .1625  (avoid exact .6 boundary; fp roundoff)
    if seq1(esl_itest.run(f'{easel} alicol --ppavg 0.95 {tmppfx}.pp.sto')) != 'A':     esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol --ppavg 0.7  {tmppfx}.pp.sto')) != 'ACD':   esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol --ppavg 0.1  {tmppfx}.pp.sto')) != 'ACDEF': esl_itest.fail()

    # --ppfrac keeps cols where a fraction >= x of residues have PP (min-of-bin) >= ppfracthresh.
    #   At default thresh 0.95, only '*' qualifies. Per-col '*' counts (of 2): c1=2 c2=0 c3=1 c4=0 c5=0.
    #   ppfrac 1.0 -> only col1 (frac 1.0):     A
    #   ppfrac 0.5 -> cols 1 and 3 (frac >=.5): AD
    if seq1(esl_itest.run(f'{easel} alicol --ppfrac 1.0 {tmppfx}.pp.sto')) != 'A':     esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol --ppfrac 0.5 {tmppfx}.pp.sto')) != 'AD':    esl_itest.fail()
    # Lower the threshold so more residues qualify: thresh 0.6 -> 5,7,8,9,* all qualify.
    #   per-col qualifying counts: c1=2 c2=2 c3=2 c4=1 c5=0; ppfrac 1.0 -> c1,c2,c3: ACD
    if seq1(esl_itest.run(f'{easel} alicol --ppfrac 1.0 --ppfracthresh 0.6 {tmppfx}.pp.sto')) != 'ACD': esl_itest.fail()

    # --ppcons|--ppavg|--ppfrac require Stockholm input.
    subprocess.run(f'{easel} reformat -o {tmppfx}.afa afa {tmppfx}.sto'.split(), check=True)
    r = esl_itest.run(f'{easel} alicol --ppcons 0.5 {tmppfx}.afa', expect_success=False)

    # -w: weighted fractions/averages.
    # seq1 has weight 3, seq2 weight 1. col4 PP = '5','7' (means .5,.7).
    #   unweighted mean = .6;  weighted mean = (3*.5 + 1*.7)/4 = .55
    # So at threshold 0.575, col4 is kept without -w (mean .6) but dropped with -w (.55).
    if seq1(esl_itest.run(f'{easel} alicol    --ppavg 0.575 {tmppfx}.ppw.sto')) != 'ACDE': esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol -w --ppavg 0.575 {tmppfx}.ppw.sto')) != 'ACD':  esl_itest.fail()
    # -w requires a fraction/average selection option.
    r = esl_itest.run(f'{easel} alicol -w --rfonly {tmppfx}.sto', expect_success=False)

    # -o sends the MSA to a file; stdout becomes a verbose diagnostic table.
    r = esl_itest.run(f'{easel} alicol --rfonly -o {tmppfx}.out.sto {tmppfx}.sto')
    if 'selection mode' not in r.stdout: esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol --rfonly {tmppfx}.out.sto')) != 'ACDEFG':   esl_itest.fail()  # output is a valid MSA

    # -q suppresses the verbose table (requires -o); stdout must be empty.
    r = esl_itest.run(f'{easel} alicol --rfonly -q -o {tmppfx}.out.sto {tmppfx}.sto')
    if r.stdout != '': esl_itest.fail()
    # -q requires -o.
    r = esl_itest.run(f'{easel} alicol --rfonly -q {tmppfx}.sto', expect_success=False)

    # --outmask writes the final all-column mask (length alen); round-trip it through --mask.
    esl_itest.run(f'{easel} alicol --rfonly --outmask {tmppfx}.mask -q -o {tmppfx}.out.sto {tmppfx}.sto')
    with open(f'{tmppfx}.mask') as f: mask = f.read().strip()
    if mask != '11011101': esl_itest.fail()   # inserts at cols 3,7 are 0
    if seq1(esl_itest.run(f'{easel} alicol --mask {tmppfx}.mask {tmppfx}.sto')) != 'ACDEFG': esl_itest.fail()

    # --outrfmask writes the final mask in consensus coords (length C); round-trip with -c --mask.
    #   --nogap keeps ACDEG; in consensus coords that drops consensus col 6 (=F): mask 111101
    esl_itest.run(f'{easel} alicol --nogap --outrfmask {tmppfx}.rfmask -q -o {tmppfx}.out.sto {tmppfx}.sto')
    with open(f'{tmppfx}.rfmask') as f: rfmask = f.read().strip()
    if rfmask != '111101': esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol -c --mask {tmppfx}.rfmask {tmppfx}.sto')) != 'ACDEG': esl_itest.fail()

    # Selection options AND together.
    #   --span 4:6 keeps cols 4,5,6 (=DEF). --nogap drops col 6, which is a gap in seq2 (AC-DE--G).
    #   The AND of the two keeps cols 4,5 = DE.
    if seq1(esl_itest.run(f'{easel} alicol --span 4:6 --nogap {tmppfx}.sto')) != 'DE': esl_itest.fail()

    # --informat / --outformat
    subprocess.run(f'{easel} reformat -o {tmppfx}.pfam pfam {tmppfx}.sto'.split(), check=True)
    if seq1(esl_itest.run(f'{easel} alicol --informat pfam --rfonly {tmppfx}.pfam')) != 'ACDEFG': esl_itest.fail()
    r = esl_itest.run(f'{easel} alicol --rfonly --outformat afa {tmppfx}.sto')
    if '>seq1' not in r.stdout: esl_itest.fail()

    # reading from stdin with '-'
    r = esl_itest.run_piped(f'cat {tmppfx}.sto', f'{easel} alicol --rfonly -')
    if seq1(r) != 'ACDEFG': esl_itest.fail()


    # --small : small-memory two-pass mode.
    #
    # It reaches the same column selections as the default path, but via the
    # legacy ESL_MSAFILE2 parser (Pfam in, Pfam out). Our test MSAs are already
    # one-block, 1-line/seq Stockholm (= Pfam), so they're fed in directly with
    # no reformatting. We don't diff small vs default output byte-for-byte
    # (margin spacing and GC-line ordering differ cosmetically); instead we
    # check the same seq1 column selections, masks, and SS_cons fixes.
    #
    # --small requires an asserted alphabet (count arrays need a known abc); our
    # write_msa MSAs are protein, so these use --amino.

    # alphabet must be asserted with --small.
    r = esl_itest.run(f'{easel} alicol --small --rfonly {tmppfx}.sto', expect_success=False)

    # --rfonly: same 6 consensus columns as the default path.
    if seq1(esl_itest.run(f'{easel} alicol --small --amino --rfonly {tmppfx}.sto')) != 'ACDEFG': esl_itest.fail()

    # Count-based selections (small_mask_by_gapfrac/symfrac, from the abc_ct array).
    if seq1(esl_itest.run(f'{easel} alicol --small --amino --mingap        {tmppfx}.sto')) != 'ACDEF-G': esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol --small --amino --nogap         {tmppfx}.sto')) != 'ACDEG':   esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol --small --amino --gapfrac 1.0   {tmppfx}.sto')) != 'ACDEF-G': esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol --small --amino --symfrac 1.0   {tmppfx}.sto')) != 'ACDEG':   esl_itest.fail()
    # --keeprf force-keeps consensus cols in small mode too: --nogap drops consensus col 6, --keeprf keeps it.
    if seq1(esl_itest.run(f'{easel} alicol --small --amino --nogap --keeprf {tmppfx}.sto')) != 'ACDEFG': esl_itest.fail()

    # PP-based selections (small_mask_by_ppavg/ppfrac from pp_ct; ppcons from msa->pp_cons).
    if seq1(esl_itest.run(f'{easel} alicol --small --amino --ppcons 0.9 {tmppfx}.pp.sto')) != 'A':     esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol --small --amino --ppcons 0.6 {tmppfx}.pp.sto')) != 'ACD':   esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol --small --amino --ppavg 0.95 {tmppfx}.pp.sto')) != 'A':     esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol --small --amino --ppavg 0.7  {tmppfx}.pp.sto')) != 'ACD':   esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol --small --amino --ppavg 0.1  {tmppfx}.pp.sto')) != 'ACDEF': esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol --small --amino --ppfrac 1.0 {tmppfx}.pp.sto')) != 'A':     esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol --small --amino --ppfrac 0.5 {tmppfx}.pp.sto')) != 'AD':    esl_itest.fail()
    if seq1(esl_itest.run(f'{easel} alicol --small --amino --ppfrac 1.0 --ppfracthresh 0.6 {tmppfx}.pp.sto')) != 'ACD': esl_itest.fail()

    # --small guards.
    #   -w is incompatible (the legacy regurgitation parser doesn't carry weights)
    r = esl_itest.run(f'{easel} alicol --small --amino -w --symfrac 1.0 {tmppfx}.sto', expect_success=False)
    #   input must be Pfam: an asserted non-Pfam --informat is rejected ({tmppfx}.afa written earlier)
    r = esl_itest.run(f'{easel} alicol --small --amino --informat afa --rfonly {tmppfx}.afa', expect_success=False)
    #   two passes can't rewind a stdin stream, so '-' is rejected
    r = esl_itest.run_piped(f'cat {tmppfx}.sto', f'{easel} alicol --small --amino --rfonly -', expect_success=False)

    # --small multi-MSA: two concatenated records in one file.
    with open(f'{tmppfx}.multi.sto', 'w') as fout:
        for f in (f'{tmppfx}.sto', f'{tmppfx}.sto'):
            with open(f) as fin: fout.write(fin.read())
    #   --rfonly applies per-MSA, so both records are emitted (2 terminating '//' lines)
    r = esl_itest.run(f'{easel} alicol --small --amino --rfonly {tmppfx}.multi.sto')
    if r.stdout.count('//') != 2: esl_itest.fail()
    #   single-MSA-only options error on a multi-MSA file (in pass 1, before any output)
    r = esl_itest.run(f'{easel} alicol --small --amino --span 1:3                 {tmppfx}.multi.sto', expect_success=False)
    r = esl_itest.run(f'{easel} alicol --small --amino --rfonly --outmask {tmppfx}.m {tmppfx}.multi.sto', expect_success=False)

    # --small --outmask round-trips to the same mask as the default path.
    esl_itest.run(f'{easel} alicol --small --amino --rfonly --outmask {tmppfx}.smask -q -o {tmppfx}.out.sto {tmppfx}.sto')
    with open(f'{tmppfx}.smask') as f: smask = f.read().strip()
    if smask != '11011101': esl_itest.fail()

    # --small base-pair fixing: removing one partner of a pair deletes the pair from
    # SS_cons, matching the default path. SS_cons "<<..>>" pairs cols 1-6 and 2-5;
    # --span 1:5 drops col 6, breaking the 1-6 pair, so col 1 becomes unpaired:  .<..>
    write_rna_ss_msa(f'{tmppfx}.rna.sto')
    ss_default = ss_cons(esl_itest.run(f'{easel} alicol         --rna --span 1:5 {tmppfx}.rna.sto'))
    ss_small   = ss_cons(esl_itest.run(f'{easel} alicol --small --rna --span 1:5 {tmppfx}.rna.sto'))
    if ss_default != '.<..>': esl_itest.fail()   # the 1-6 pair is gone; the 2-5 pair survives
    if ss_small   != ss_default: esl_itest.fail()


    # Cleanup
    for tmpfile in glob.glob(f'{tmppfx}.*'):
        os.remove(tmpfile)
    print('ok')



if __name__ == "__main__":
    main()
