#! /usr/bin/env python3

""" Integration test for `easel construct`

Usage: easel-construct-itest.py <builddir> <srcdir> <tmppfx>
  <builddir>: path to Easel build dir. `easel` miniapp is <builddir>/miniapps/easel
  <srcdir>:   path to Easel src dir.
  <tmppfx>:   prefix we're allowed to use to create tmp files in current working dir.

"""


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

# Test files

TESTMSA1 = """\
# STOCKHOLM 1.0
#=GF AU Infernal 0.1

simpex-1           --AGA-CUUCG-GUCGCUCG-UAACAG
#=GR simpex-1   SS ..:<<-<____->>>-<<-<.___>>>
simpex-2           aaAAUACGUCGGCUG-AAUACCCAGUA
#=GR simpex-2   SS ..::<<<____>->>--<-.<___>>:
simpex-3           --ACGUUUUG-GAACGGG-U-CCAACC
#=GR simpex-3   SS ..::<<<____>->>-<<-<.___>>>
#=GC SS_cons       ..::<<<____>->>-<<-<.___>>>
#=GC RF            ..AAgaCUUCGGAucgggCg.AcAccc
//
"""

TESTMSA2 = """\
# STOCKHOLM 1.0
#=GF AU Infernal 0.1

simpex-1           --AGA-CTTCG-GTCGCTCG-TAACAG
simpex-2           aaAATACGTCGGCTG-AATACCCAGTA
simpex-3           --ACGTTTTG-GAACGGG-T-CCAACC
#=GC SS_cons       ..::<<<____>->>-<<-<.___>>>
#=GC RF            ..AAgaCTTCGGAtcgggCg.AcAccc
//
"""

TESTMSA3 = """\
# STOCKHOLM 1.0
#=GF AU Infernal 0.1

simpex-1           --AGA-CUUCG-GUCGCUCG-UAACAG
#=GR simpex-1   SS ..:<<-<____->>>-<<-<.___>>>
simpex-2           aaAAUACGUCGGCUG-AAUACCCAGUA
#=GR simpex-2   SS ..::<<<____>->>--<-.<___>>:
simpex-3           --ACGUUUUG-GAACGGG-U-CCAACC
#=GR simpex-3   SS ..::<<<____>->>-<<-<.___>>>
//
"""

with open(f'{tmppfx}.1', 'w', encoding='utf-8') as f:
    f.write(TESTMSA1)
with open(f'{tmppfx}.2', 'w', encoding='utf-8') as f:
    f.write(TESTMSA2)
with open(f'{tmppfx}.3', 'w', encoding='utf-8') as f:
    f.write(TESTMSA3)


# `-h` help
r = esl_itest.run(f'{easel} construct -h')


# basic
r = esl_itest.run(f'{easel} construct {tmppfx}.1')
if re.search(r'^\s*simpex-2\s+5\s+4\s+4\s+1\s*$',
             r.stdout, flags=re.MULTILINE) is None:
    esl_itest.fail()
if re.search(r'^\s*SS_cons\(consensus\)\s+6\s+6\s+6\s+0\s*$',
             r.stdout, flags=re.MULTILINE) is None:
    esl_itest.fail()
if re.search(r'^#\s+13/\s*17\s+\(0.765\)\s+overlap\s*$',
             r.stdout, flags=re.MULTILINE) is None:
    esl_itest.fail()

r2 = esl_itest.run(f'{easel} construct {tmppfx}.2')
if re.search(r'^\s+SS_cons\(consensus\)\s+6\s+6\s+6\s+0\s*$',
             r2.stdout, flags=re.MULTILINE) is None:
    esl_itest.fail()

r3 = esl_itest.run(f'{easel} construct {tmppfx}.3')
if re.search(r'^\s*simpex-2\s+5\s*$',
             r3.stdout, flags=re.MULTILINE) is None:
    esl_itest.fail()

# --rna   tmppfx.1 is RNA
#
r3 = esl_itest.run(f'{easel} construct --rna {tmppfx}.1')    
if r3.stdout != r.stdout:
    esl_itest.fail()

# --dna   tmppfx.2 is DNA    
r3 = esl_itest.run(f'{easel} construct --dna {tmppfx}.2')    
if r3.stdout != r2.stdout:
    esl_itest.fail()

# -a   print info on all conflicting bps in individual structures
#
r = esl_itest.run(f'{easel} construct -a {tmppfx}.1')
if re.search(r'^More than 1 right mates for left  mate\s+7\s+7:\s+12 bp exists in\s+2/\s+3 seqs',
             r.stdout, flags=re.MULTILINE) is None:
    esl_itest.fail()
if re.search(r'^More than 1 right mates for left  mate\s+7\s+7:\s+13 bp exists in\s+1/\s+3 seqs',
             r.stdout, flags=re.MULTILINE) is None:
    esl_itest.fail()
if re.search(r'^More than 1 left  mates for right mate\s+25\s+20:\s+25 bp exists in\s+2/\s+3 seqs',
             r.stdout, flags=re.MULTILINE) is None:
    esl_itest.fail()
if re.search(r'^More than 1 left  mates for right mate\s+25\s+21:\s+25 bp exists in\s+1/\s+3 seqs',
             r.stdout, flags=re.MULTILINE) is None:
    esl_itest.fail()


#  -v   be verbose
#
r = esl_itest.run(f'{easel} construct -v {tmppfx}.1')
if re.search(r'^\s*simpex-2\s+5\s+4\s+4\s+1\s*$',
             r.stdout, flags=re.MULTILINE) is None:
    esl_itest.fail()
if re.search(r'^\s*SS_cons\(consensus\)\s+6\s+6\s+6\s+0\s*$',
             r.stdout, flags=re.MULTILINE) is None:
    esl_itest.fail()
if re.search(r'^#\s+13/\s*17\s+\(0.765\)\s+overlap\s*$',
             r.stdout, flags=re.MULTILINE) is None:
    esl_itest.fail()
if re.search(r'^ali:\s+1 seq\s+0 \(simpex-1\) bp\s+5:\s+14 conflicts with consensus bp\s+5:\s+15',
             r.stdout, flags=re.MULTILINE) is None:
    esl_itest.fail()


#  -x       set SS_cons as max set of non-conflicting bps from indi SSs
#  -o       output new MSA with SS_cons information (required by -x and other SS_cons opts)
#  --pfam   output alignment in Pfam (non-interleaved, 1 line/seq) format
#
r = esl_itest.run(f'{easel} construct -x --pfam -o {tmppfx}.sto {tmppfx}.1')
with open(f'{tmppfx}.sto', encoding='utf-8') as f:
    msa_output = f.read()
if re.search(r'^#=GC SS_cons     ::::::::::::::::<<_______>>',
             msa_output, flags=re.MULTILINE) is None:
    esl_itest.fail()

#  -r   remove SS_cons basepairs that conflicts with > 0 indi SS
#
r = esl_itest.run(f'{easel} construct -r -o {tmppfx}.sto {tmppfx}.1')
with open(f'{tmppfx}.sto', encoding='utf-8') as f:
    msa_output = f.read()
if re.search(r'^#=GC SS_cons\s+:::::<_______>::<<_______>>',
             msa_output, flags=re.MULTILINE) is None:
    esl_itest.fail()

#  -c   set SS_cons as indi SS with max bps consistent with SS_cons
#
r = esl_itest.run(f'{easel} construct -c -o {tmppfx}.sto {tmppfx}.1')
with open(f'{tmppfx}.sto', encoding='utf-8') as f:
    msa_output = f.read()
if re.search(r'^#=GC SS_cons\s+::::<<<____>->>:<<-<____>>>',
             msa_output, flags=re.MULTILINE) is None:
    esl_itest.fail()

#  --rfc  with -c, set RF annotation as seq SS_cons structure comes from
#
r = esl_itest.run(f'{easel} construct -c --rfc -o {tmppfx}.sto {tmppfx}.1')
with open(f'{tmppfx}.sto', encoding='utf-8') as f:
    msa_output = f.read()
if re.search(r'#=GC SS_cons\s+::::<<<____>->>:<<-<____>>>',
             msa_output, flags=re.MULTILINE) is None:
    esl_itest.fail()
if re.search(r'#=GC RF\s+--ACGUUUUG-GAACGGG-U-CCAAC',
             msa_output, flags=re.MULTILINE) is None:
    esl_itest.fail()

#  --indi <s>  define SS_cons as individual SS for sequence <x>
#
r = esl_itest.run(f'{easel} construct --indi simpex-2 -o {tmppfx}.sto {tmppfx}.1')
with open(f'{tmppfx}.sto', encoding='utf-8') as f:
    msa_output = f.read()
if re.search(r'^#=GC SS_cons     ::::<<<____>->>::<--<___>>:',
              msa_output, flags=re.MULTILINE) is None:
    esl_itest.fail()

#  --rfindi    with --indi <x>, define RF annotation as <x>
#
r = esl_itest.run(f'{easel} construct --indi simpex-2 --rfindi -o {tmppfx}.sto {tmppfx}.1')
with open(f'{tmppfx}.sto', encoding='utf-8') as f:
    msa_output = f.read()
if re.search(r'^#=GC SS_cons\s+::::<<<____>->>::<--<___>>:',
              msa_output, flags=re.MULTILINE) is None:
    esl_itest.fail()
if re.search(r'^#=GC RF\s+AAAAUACGUCGGCUG-AAUACCCAGUA',
              msa_output, flags=re.MULTILINE) is None:
    esl_itest.fail()

#  --ffreq <x>    aln cols i:j become SS_cons bps if paired in > <x> indi SS
#
r = esl_itest.run(f'{easel} construct --ffreq 0.6 -o {tmppfx}.sto {tmppfx}.1')
with open(f'{tmppfx}.sto', encoding='utf-8') as f:
    msa_output = f.read()
if re.search(r'^#=GC SS_cons\s+::::<<<____>->>:<<-<____>>>',
             msa_output, flags=re.MULTILINE) is None:
    esl_itest.fail()

r = esl_itest.run(f'{easel} construct --ffreq 0.7 -o {tmppfx}.sto {tmppfx}.1')
with open(f'{tmppfx}.sto', encoding='utf-8') as f:
    msa_output = f.read()
if re.search(r'^#=GC SS_cons\s+:::::::::::::::::<_______>:',
             msa_output, flags=re.MULTILINE) is None:
    esl_itest.fail()

#  --fmin      same as --ffreq but find min <x> that gives consistent SS_cons
#
r = esl_itest.run(f'{easel} construct --fmin -o {tmppfx}.sto {tmppfx}.1')
with open(f'{tmppfx}.sto', encoding='utf-8') as f:
    msa_output = f.read()
if re.search(r'^#=GC SS_cons\s+::::<<<____>->>:<<-<____>>>',
             msa_output, flags=re.MULTILINE) is None:
    esl_itest.fail()

#  -l <f>     list seqs w/> 0 indi bp that conflicts w/a SS_cons bp to file <f>
#
r = esl_itest.run(f'{easel} construct --fmin -l {tmppfx}.list -o {tmppfx}.sto {tmppfx}.1')
with open(f'{tmppfx}.list', encoding='utf-8') as f:
    list_output = f.read()
if re.search(r'^simpex-1', list_output, flags=re.MULTILINE) is None:
    esl_itest.fail()
if re.search(r'^simpex-2', list_output, flags=re.MULTILINE) is None:
    esl_itest.fail()

#  --lmax <n> with -l, change maximum allowed conflicts of 0 to <x>  [0]  (n>=0)
#
# simpex-2 doesn't appear in list now
#
r = esl_itest.run(f'{easel} construct --lmax 1 --fmin -l {tmppfx}.list -o {tmppfx}.sto {tmppfx}.1')
with open(f'{tmppfx}.list', encoding='utf-8') as f:
    list_output = f.read()
if re.search(r'^simpex-1', list_output, flags=re.MULTILINE) is None:
    esl_itest.fail()
if re.search(r'^simpex-2', list_output, flags=re.MULTILINE) is not None:
    esl_itest.fail()


# Cleanup
for tmpfile in glob.glob(f'{tmppfx}.*'):
    os.remove(tmpfile)

print('ok')
