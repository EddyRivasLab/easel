#! /usr/bin/env python3

# Integration test for `easel alimerge`
#
# Usage: easel-alimerge-itest.py <builddir> <srcdir> <tmppfx>
#   <builddir>: path to Easel build dir. `easel` miniapp is <builddir>/miniapps/easel
#   <srcdir>:   path to Easel src dir.
#   <tmppfx>:   prefix we're allowed to use to create tmp files in current working dir.
#
# Also uses `easel alicol --mingap`.

import filecmp
import glob
import os
import random
import re
import subprocess
import sys
import esl_itest


# msasplit(msafile, easel, n, tmppfx)
#
# Split the MSA in <msafile> into <K> new smaller MSAs, remove columns
# that are now all-gap using "{easel} alicol", and write the smaller
# MSAs to <n> new files <tmppfx>.1.sto, ..., <tmppfx>.<K>.sto.
#
# <msafile> must only contain 1 MSA.
#
# It's important to keep the seqs in their original order. `easel
# alistat` calculation of avg pid has a tiny/annoying order dependence
# because of fp roundoff error.
#
# Wrote this as a function that I can lift somewhere else - might be
# useful elsewhere someday.
#
def msasplit(msafile, easel, K, tmppfx):
    if K <= 1: sys.exit("<K> needs to be > 1")

    # Slurp in the whole file, so we can make multiple passes through it.
    with open(msafile) as f:
        msa = f.read()

    # Get a list of sequence names, <S>. Keep it in original order.
    # Stockholm can be multiblock; only keep the
    # first time we see each seqname.
    # 
    S            = []
    seqname_seen = set()
    n_msa        = 0
    for line in msa.splitlines():
        if re.match(r'//', line):
            n_msa += 1
        elif m := re.match(r'([^#]\S+)', line):
            if not m.group(1) in seqname_seen: S.append(m.group(1)) 
            seqname_seen.add(m.group(1))
    if n_msa > 1: sys.exit("<msafile> must contain only 1 MSA")
    
    # Assign each sequence to k=0..K-1 splits
    # This is fairly generic code for uniformly sampling
    # a "composition" of K splits from N things.
    #
    N   = len(S)

    assign = {}                                        # Split the 0..N-1 elems of list <S> by setting K endpoints 
    ends   = sorted(random.sample(range(0,N-1), K-1))  # ... first sample K-1 endpoints 0..N-2 uniformly without replacement. 
    ends  += [ N-1 ]                                   # ... then add sentinel at end. ends[k=0..K-1] works for all segments including K-1 now.
    k      = 0
    for i in range(N):
        assign[S[i]] = k+1
        if i == ends[k]:
            k += 1
    # Now assign[0..N-1] = 1..K is the assignment of each sequence to a split

    # Make <nfiles> passes through the MSA, one pass for each of the K new sets.
    # 
    for k in range(1,K+1):
        with open(f'{tmppfx}.{k}.tmp', 'w') as f:
            for line in msa.splitlines():
                if   (m := re.match(r'#=GR\s+(\S+)', line)):
                    if assign[m.group(1)] == k:              print(line, file=f)
                elif (m := re.match(r'#=GS\s+(\S+)', line)):
                    if assign[m.group(1)] == k:              print(line, file=f)
                elif       re.match(r'#',            line):  print(line, file=f)
                elif       re.match(r'//',           line):  print(line, file=f)
                elif       re.match(r'\s*$',         line):  print(line, file=f)
                elif (m := re.match(r'(\S+)',        line)):
                    if assign[m.group(1)] == k:              print(line, file=f)

    # Run each output file through easel alicol to remove all-gap columns
    # Need --keeprf, so #=GC RF lines remain compatible
    #
    for k in range(1,K+1):
        r   = esl_itest.run(f'{easel} alicol --informat pfam --mingap --keeprf -o {tmppfx}.{k}.sto {tmppfx}.{k}.tmp')
        os.remove(f'{tmppfx}.{k}.tmp')


# write_small_testmsa()
#
# A small test case that exercises the various sorts of Stockholm
# annotation, parsed and unparsed; and aseq fragments, with flanking
# ~'s.
#
# Splitting this (with msasplit) and merging it back with `easel
# alimerge` almost but not quite diffs clean against the original. The
# stray issue is the inconsistent use of - vs . as gap characters.  To
# work around this issue, we don't want to replace all gap chars or
# we'd wipe out the ~'s, which are part of the test.  Instead, the
# itest code will pass the `easel alimerge` output through `easel
# reformat --replace=-:.`.
#
# (Revisit this if we ever clean up gap character inconsistencies in
# Stockholm/Pfam format files.)
#
def write_small_testmsa(outfile):
    msa = """\
# STOCKHOLM 1.0
#=GF ID test
#=GF AC RF12345
#=GF XX unparsed GF annotation

#=GS seq1 DE seq1 test description
#=GS seq2 DE seq2 test description

#=GS seq1 XX unparsed GS annotation

seq1         .AGGGAA.AAAACCCUUUUAA
seq2         ~~~~GAAAAAAACCCUUUUA~
seq3         AAGGGAA..AAACCCUU~~~~
#=GR seq3 SS ..<<<_....__>>>::....
#=GR seq3 XX 1111111..11111111....
#=GC SS_cons ..<<<_~~~~__>>>::::..
#=GC RF      ..xxxx~~~~xxxxxxxxx..
#=GC XX      ..2222....222222222..
//"""
    with open(outfile, 'w') as f:
        print(msa, file=f)


    
    
def main():
    progs_used = [ 'miniapps/easel' ]
    files_used = [ 'testsuite/example-rna.sto' ]      # RNase P. 116 seqs. 845 cols. 45% av pid.

    (builddir, srcdir, tmppfx) = esl_itest.getargs(sys.argv)
    esl_itest.check_files(srcdir,   files_used)
    esl_itest.check_progs(builddir, progs_used)

    easel = f'{builddir}/miniapps/easel'

    # We test for identical `easel alistat` output for original vs. merged
    #
    r = esl_itest.run(f'{easel} alistat {srcdir}/testsuite/example-rna.sto')
    orig_alistat = r.stdout

    # A 2-way split
    #
    msasplit(f'{srcdir}/testsuite/example-rna.sto', easel, 2, tmppfx)
    r  = esl_itest.run(f'{easel} alimerge {tmppfx}.1.sto {tmppfx}.2.sto')
    r2 = subprocess.run(f'{easel} alistat -'.split(), check=True, encoding='utf-8', capture_output=True, input=r.stdout)
    if orig_alistat != r2.stdout: esl_itest.fail()

    # A 3-way split
    #
    msasplit(f'{srcdir}/testsuite/example-rna.sto', easel, 3, tmppfx)
    r  = esl_itest.run(f'{easel} alimerge {tmppfx}.1.sto {tmppfx}.2.sto {tmppfx}.3.sto')
    r2 = subprocess.run(f'{easel} alistat -'.split(), check=True, encoding='utf-8', capture_output=True, input=r.stdout)
    if orig_alistat != r2.stdout: esl_itest.fail()

    # One multi-MSA file with a 3-way split, through stdin
    #
    msasplit(f'{srcdir}/testsuite/example-rna.sto', easel, 3, tmppfx)
    r  = esl_itest.run_piped(f'cat {tmppfx}.1.sto {tmppfx}.2.sto {tmppfx}.3.sto', f'{easel} alimerge -')
    r2 = subprocess.run(f'{easel} alistat -'.split(), check=True, encoding='utf-8', capture_output=True, input=r.stdout)
    if orig_alistat != r2.stdout: esl_itest.fail()
    
    # -h   help 
    r = esl_itest.run(f'{easel} alimerge -h')

    # use the same test split for several option tests below...
    msasplit(f'{srcdir}/testsuite/example-rna.sto', easel, 2, tmppfx)

    # -o <outfile>
    r  = esl_itest.run(f'{easel} alimerge -o {tmppfx}.merged.sto {tmppfx}.1.sto {tmppfx}.2.sto')
    r2 = esl_itest.run(f'{easel} alistat {tmppfx}.merged.sto')
    if orig_alistat != r2.stdout: esl_itest.fail()

    # --small
    # requires one-block Pfam format on input
    r  = esl_itest.run(f'{easel} alimerge --small {tmppfx}.1.sto {tmppfx}.2.sto')
    r2 = subprocess.run(f'{easel} alistat -'.split(), check=True, encoding='utf-8', capture_output=True, input=r.stdout)
    if orig_alistat != r2.stdout: esl_itest.fail()

    # --rfonly
    # (with and without --small)
    # by comparison to `easel alicol --rfonly` on the original complete MSA
    #
    r  = esl_itest.run (f'{easel} alicol --rfonly {srcdir}/testsuite/example-rna.sto')
    r2 = subprocess.run(f'{easel} alistat -'.split(), check=True, encoding='utf-8', capture_output=True, input=r.stdout)
    r  = esl_itest.run (f'{easel} alimerge --rfonly {tmppfx}.1.sto {tmppfx}.2.sto')
    r3 = subprocess.run(f'{easel} alistat -'.split(), check=True, encoding='utf-8', capture_output=True, input=r.stdout)
    r  = esl_itest.run (f'{easel} alimerge --rfonly --small {tmppfx}.1.sto {tmppfx}.2.sto')
    r4 = subprocess.run(f'{easel} alistat -'.split(), check=True, encoding='utf-8', capture_output=True, input=r.stdout)
    if r3.stdout != r2.stdout: esl_itest.fail()
    if r4.stdout != r2.stdout: esl_itest.fail()

    # --outformat (incompatible with --small, which always outputs Pfam format)
    r  = esl_itest.run (f'{easel} alimerge --outformat afa {tmppfx}.1.sto {tmppfx}.2.sto')
    r2 = subprocess.run(f'{easel} alistat --informat afa -'.split(),                        check=True, encoding='utf-8', capture_output=True, input=r.stdout)
    r  = subprocess.run(f'{easel} reformat afa {srcdir}/testsuite/example-rna.sto'.split(), check=True, encoding='utf-8', capture_output=True)
    r3 = subprocess.run(f'{easel} alistat --informat afa -'.split(),                        check=True, encoding='utf-8', capture_output=True, input=r.stdout)
    if r3.stdout != r2.stdout: esl_itest.fail()

    r = esl_itest.run(f'{easel} alimerge --small --outformat afa {tmppfx}.1.sto {tmppfx}.2.sto', expect_success=False)

    # -v: verbose, requires -o
    r = esl_itest.run(f'{easel} alimerge -v -o {tmppfx}.sto {tmppfx}.1.sto {tmppfx}.2.sto')
    if 'file name' not in r.stdout: esl_itest.fail()
    r2 = esl_itest.run(f'{easel} alistat {tmppfx}.merged.sto')
    if orig_alistat != r2.stdout: esl_itest.fail()

    r = esl_itest.run(f'{easel} alimerge -v {tmppfx}.1.sto {tmppfx}.2.sto', expect_success=False)

    # An example MSA that exercises Stockholm parsed and unparsed annotation, and aseq fragment ~ marking.
    #
    write_small_testmsa(f'{tmppfx}.sto')
    msasplit(f'{tmppfx}.sto', easel, 2, tmppfx)
    r  = esl_itest.run(f'{easel} alimerge {tmppfx}.1.sto {tmppfx}.2.sto')
    r2 = subprocess.run(f'{easel} reformat -o {tmppfx}.out --replace=-:. pfam -'.split(), check=True, encoding='utf-8', input=r.stdout)
    if not filecmp.cmp(f'{tmppfx}.sto', f'{tmppfx}.out', shallow=False): esl_itest.fail()


    # Cleanup
    for tmpfile in glob.glob(f'{tmppfx}.*'):
        os.remove(tmpfile)
    print('ok')



if __name__ == "__main__":
    main()
