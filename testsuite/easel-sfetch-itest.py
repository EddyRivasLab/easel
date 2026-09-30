#! /usr/bin/env python3

# Integration test for `easel sfetch` miniapp
#
# Usage: easel-sfetch-itest.py <builddir> <srcdir> <tmppfx>
#   <builddir>: path to Easel build dir. `easel` miniapp is <builddir>/miniapps/easel
#   <srcdir>:   path to Easel src dir.
#   <tmppfx>:   prefix we're allowed to use to create tmp files in current working dir.
#
import glob
import os
import re
import shutil
import subprocess
import sys
import esl_itest

files_used = [ 'testsuite/example-genbank.gb',     # 4 phage DNA seqs:    NC_047788 NC_055916 NC_007046 NC_049972, accessions same as names
               'testsuite/example-uniprot.dat',    # 4 protein seqs:      MNME_BEII9 DEF_RICCK GPMI_YERP3 FABZ_PROM2; accessions B2IJQ3 A8EXV2 A7FCU8 A8G6E7
               'testsuite/example-uniprot.fa' ]    # same 4 seqs as .dat: sp|B2IJQ3|MNME_BEII9, sp|A8EXV2|DEF_RICCK, sp|A7FCU8|GPMI_YERP3, sp|A8G6E7|FABZ_PROM2

progs_used = [ 'miniapps/easel' ]

seq1 = 'ACGTACGTAG' * 12 + 'ACGTA'  # 125 residues, written 10 per line
seq2 = 'AGCTTAGCTT' * 7             #  70 residues, all on one line

(builddir, srcdir, tmppfx) = esl_itest.getargs(sys.argv)
esl_itest.check_files(srcdir,   files_used)
esl_itest.check_progs(builddir, progs_used)

# -h
r = esl_itest.run(f'{builddir}/miniapps/easel sfetch -h')

##
## Two passes. Pass 1: unindexed. Pass 2: indexed.
##
for p in [ 'unindexed', 'indexed']:

    # Make copies of three example files from Easel testsuite directory.
    shutil.copyfile('{}/testsuite/example-genbank.gb'.format(srcdir), '{}.gb'.format(tmppfx))
    shutil.copyfile('{}/testsuite/example-uniprot.dat'.format(srcdir), '{}.dat'.format(tmppfx))
    shutil.copyfile('{}/testsuite/example-uniprot.fa'.format(srcdir), '{}.fa'.format(tmppfx))

    # Also create two FASTA files whose records disagree on a line length:
    # seq1 has 10 residues per line, seq2 is all on one line. Such files
    # aren't "well-formatted" (see esl_ssi.md), so `easel sindex` must not
    # mark them for fast subseq lookup. If it does, indexed `-c` fetches come
    # back silently shifted, by a byte per line. The .noeol.fa version also
    # leaves the newline off its last line, a second place such a bug hides.
    for (sfx, eol) in [ ('mixed', '\n'), ('noeol', '') ]:
        with open('{}.{}.fa'.format(tmppfx, sfx), 'w') as f:
            f.write('>seq1\n')
            for i in range(0, len(seq1), 10): f.write(seq1[i:i+10] + '\n')
            f.write('>seq2\n' + seq2 + eol)

    # index them, on that pass
    if p == 'indexed':
        r = esl_itest.run('{}/miniapps/easel sindex {}.gb'.format(builddir,tmppfx))
        r = esl_itest.run('{}/miniapps/easel sindex {}.dat'.format(builddir,tmppfx))
        r = esl_itest.run('{}/miniapps/easel sindex {}.fa'.format(builddir,tmppfx))
        for sfx in [ 'mixed', 'noeol' ]:
            r = esl_itest.run('{}/miniapps/easel sindex {}.{}.fa'.format(builddir,tmppfx,sfx))

    # by name, verbatim fetch: output is GenBank format
    r = esl_itest.run('{0}/miniapps/easel sfetch {1}.gb NC_055916'.format(builddir, tmppfx))            
    if re.search(r'^LOCUS\s+NC_055916\s+17056 bp\s+DNA\s+linear', r.stdout) == None: esl_itest_fail()

    # by accession, verbatim fetch: output is EMBL/Uniprot format
    r = esl_itest.run('{0}/miniapps/easel sfetch {1}.dat A8EXV2'.format(builddir, tmppfx))              
    if re.search(r'^ID\s+DEF_RICCK\s+Reviewed;\s+175 AA\.', r.stdout) == None: esl_itest_fail()

    # from a stream is *not* verbatim; now you get FASTA
    r = esl_itest.run_piped('cat {}.dat'.format(tmppfx), '{}/miniapps/easel sfetch - A8EXV2'.format(builddir))  
    if re.search(r'^>DEF_RICCK\s+A8EXV2', r.stdout) == None: esl_itest.fail()

    # '.' gives you the first sequence in seqfile (verbatim, here)
    r = esl_itest.run('{0}/miniapps/easel sfetch {1}.gb .'.format(builddir, tmppfx))            
    if re.search(r'^LOCUS\s+NC_047788\s+18023 bp\s+DNA\s+linear', r.stdout) == None: esl_itest.fail()

    # -o
    r = esl_itest.run('{0}/miniapps/easel sfetch -o {1}.out {1}.gb NC_055916'.format(builddir, tmppfx))            
    if re.search(r'^Retrieved sequence NC_055916\.', r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()
    r = esl_itest.run('{0}/miniapps/easel seqstat {1}.out'.format(builddir, tmppfx))
    if re.search(r'^Format:\s+GenBank(?s:.+)^Total # residues:\s+17056', r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()

    # -o will not overwrite an existing file
    r = esl_itest.run('{0}/miniapps/easel sfetch -o {1}.out {1}.gb NC_055916'.format(builddir, tmppfx), expect_success=False)

    # -f will allow it to
    r = esl_itest.run('{0}/miniapps/easel sfetch -f -o {1}.out {1}.gb NC_055916'.format(builddir, tmppfx))

    # -n : output is in FASTA, can't be a verbatim fetch if you change the name
    #      Also, create a new ${tmppfx}.fa2 file with sequence named ${tmppfx}.1 so we can test -O below
    r = esl_itest.run('{0}/miniapps/easel sfetch -n {1}.1 -fo {1}.fa2 {1}.gb NC_055916'.format(builddir, tmppfx))            
    if re.search(r'^Retrieved sequence NC_055916\.', r.stdout, flags=re.MULTILINE) == None: esl_itest.fail()

    # -O   : to test this, we first used -n above to create a new file ${tmppfx}.fa2 with a seq named ${tmppfx}.1.
    r = esl_itest.run('{0}/miniapps/easel sfetch -O {1}.fa2 {1}.1'.format(builddir, tmppfx))

    # -O will not overwrite an existing file
    r = esl_itest.run('{0}/miniapps/easel sfetch -O {1}.fa2 {1}.1'.format(builddir, tmppfx), expect_success=False)

    # -f will allow it to
    r = esl_itest.run('{0}/miniapps/easel sfetch -fO {1}.fa2 {1}.1'.format(builddir, tmppfx))

    # -r reverse complements; and is not a verbatim fetch, so output is FASTA
    r  = esl_itest.run('{0}/miniapps/easel sfetch -r {1}.gb NC_007046'.format(builddir, tmppfx))
    r2 = subprocess.run('{}/miniapps/easel seqstat -'.format(builddir).split(), check=True, encoding='utf-8', capture_output=True, input=r.stdout)
    if re.search(r'^Format:\s+FASTA(?s:.+)^Total # residues:\s+18199', r2.stdout, flags=re.MULTILINE) == None: esl_itest.fail()

    # -r on an obviously not-DNA file is an error
    r  = esl_itest.run('{0}/miniapps/easel sfetch -r {1}.dat GPMI_YERP3'.format(builddir, tmppfx), expect_success=False)
    if re.search(r'^Failed to reverse complement', r.stderr) == None: esl_itest.fail()

    # --informat
    r  = esl_itest.run('{0}/miniapps/easel sfetch --informat genbank {1}.gb NC_007046'.format(builddir, tmppfx))

    # -c  subseq fetching
    r  = esl_itest.run('{0}/miniapps/easel sfetch -c 101..200 {1}.dat .'.format(builddir,tmppfx))
    if re.search(r'^>MNME_BEII9\/101-200', r.stdout) == None: esl_itest.fail()

    # -c  must return exactly the requested residues, indexed or not, even
    #     from the two not-well-formatted files we made above.
    for sfx in [ 'mixed', 'noeol' ]:
        for (name, seq) in [ ('seq1', seq1), ('seq2', seq2) ]:
            for (start, end) in [ (1,1), (1,10), (9,12), (11,20), (21,30), (61,70), (1,len(seq)) ]:
                r    = esl_itest.run('{0}/miniapps/easel sfetch -c {1}..{2} {3}.{4}.fa {5}'.format(builddir, start, end, tmppfx, sfx, name))
                got  = ''.join(r.stdout.splitlines()[1:])
                what = '{} {}.fa {} -c {}..{}'.format(p, sfx, name, start, end)
                if re.search(r'^>{0}\/{1}-{2}'.format(name, start, end), r.stdout) == None:
                    esl_itest.fail('{}: bad name/coords line'.format(what))
                if got != seq[start-1:end]:
                    esl_itest.fail('{}: got {}, not {}'.format(what, got, seq[start-1:end]))

    for tmpfile in glob.glob('{}.*'.format(tmppfx)): os.remove(tmpfile)



print('ok')
