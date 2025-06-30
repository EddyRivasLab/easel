#! /usr/bin/env python3

# Check all options in the .c file that they're documented in .man.in and tested in -itest.py.
#
# Usage:    miniapp-optcheck.py <.c file> <.man page> <itest.py script>
# Example:  devkit/miniapp-optcheck.py miniapps/cmd_afetch.c miniapps/easel-afetch.man.in testsuite/easel-afetch-itest.py
#
# In the .c file, parses our usual `static ESL_OPTIONS` structure.
#
# In the .man.in file, parses everything between `.SH .*OPTIONS` and
# `.SH SEE ALSO` for `.TP` indented paragraphs, which are always
# options.
#
# In the <miniapp>.itest.py file, parses for `esl_itest.run('easel
# <miniapp> <args>')` commands that match the <miniapp> we're testing,
# then checks for any arg that starts with a `-`. Allows concat of
# simple options (e.g. '-abc` is `-a -b -c`.). For options tested
# in a way that can't be parsed like this, also checks for formatted
# comment lines of form:
#    # optcheck assertion: <comma-sep list of options>'
#
import getopt
import os
import re
import sys

def process_cfile(cfile):
    with open(cfile) as f:
        in_opts = False
        optlist = []
        for line in f:
            if m := re.match(r'static ESL_OPTIONS ',   line):
                in_opts = True
            elif m := re.match(r'\s*\{\s*0\s*(,\s*0\s*){9}\},', line):
                in_opts = False
            elif in_opts:
                if m := re.match(r'\s+\{\s*"(-\S+)",\s*eslARG', line):
                    optlist.append(m.group(1))
    return set(optlist)

def process_manfile(manfile):
    with open(manfile) as f:
        in_opts        = False
        expect_optline = False
        optlist        = []
        for line in f:
            if m := re.match(r'.SH .*OPTIONS', line):    # Had to match .SH OPTIONS, .SH GENERAL OPTIONS, .SH BASIC OPTIONS
                in_opts = True
            elif m := re.match(r'.SH SEE ALSO', line):
                in_opts = False
            elif in_opts:
                if m := re.match(r'.TP', line):
                    expect_optline = True
                elif expect_optline:
                    if m := re.match(r'\.BI?\s+(\\-(?:\\-)?\S+)', line):
                        optlist.append(re.sub(r'\\-', '-', m.group(1)))
                        expect_optline = False
    return set(optlist)

def process_itestfile(itestfile):
    if m := re.search(r'easel-(\S+)-itest.py', itestfile):
        miniapp = m.group(1)
    else: sys.exit('Expect itestfile name to look like easel-*-itest.py')

    optlist = []
    with open(itestfile) as f:
        for line in f:
            # This pattern needs to work for a variety of ways the test scripts call esl_itest.run or esl_itest.run_piped,
            # on things like '{0}/miniapps/easel afetch' or f'{easel} afetch'. We look for 'esl_itest.run*' followed by
            # any quoted string containing 'easel <miniapp>' with or without f-string braces around the easel command. 
            if m := re.search(r"""esl_itest\.run.+(['"]).*[\{]?easel[\}]?\s+""" + re.escape(miniapp) + r"""\s+(.*?)\1""", line):
                argv = m.group(2).split()

                for a in argv:  # We don't know when to stop parsing options. It's conceivable that we could mistake a mandatory arg for an option, if it starts with '-'.
                    if    m := re.fullmatch(r'-[^-]', a):
                        optlist.append(a)                              # -a      simple option
                    elif  m := re.match(r'-(\w+)', a):
                        for c in m.group(1):
                            optlist.append(f'-{c}')                    # -abc    concatenated options
                    elif  m := re.match(r'--', a):
                        optlist.append(a)                              # --long  long form opts
            elif m := re.match(r'#\s+optcheck assertion:\s+(.+)', line):
                for opt in m.group(1).split(','):
                    optlist.append(opt)

    return set(optlist)
                    

def main():
    (opts, args) = getopt.getopt(sys.argv[1:], '1')

    if len(args) != 3:
        sys.exit('Usage: miniapp-optcheck.py <.c file> <.man page> <itest.py script>')

    do_oneline = False
    for opt, arg in opts:
        if opt == '-1':  do_oneline = True

    cfile       = args[0]
    manfile     = args[1]; 
    itestfile   = args[2]; 

    cfile_optset     = process_cfile(cfile)         if os.path.isfile(cfile)     else sys.exit(f'failed to find or open C source file {cfile}')
    manfile_optset   = process_manfile(manfile)     if os.path.isfile(manfile)   else None
    itestfile_optset = process_itestfile(itestfile) if os.path.isfile(itestfile) else None

    undoc_opts       = cfile_optset - manfile_optset    if manfile_optset   is not None else []
    untest_opts      = cfile_optset - itestfile_optset  if itestfile_optset is not None else []
    extradoc_opts    = manfile_optset - cfile_optset    if manfile_optset   is not None else []
    extratest_opts   = itestfile_optset - cfile_optset  if itestfile_optset is not None else []

    if do_oneline:
        if   manfile_optset is None:  col1 = '[no manpage]'
        elif len(undoc_opts) > 0:     col1 = f'[{len(undoc_opts)} undocumented]'
        else:                         col1 = 'ok'

        if   itestfile_optset is None: col2 = '[no itest]'
        elif len(untest_opts) > 0:     col2 = f'[{len(untest_opts)} undocumented]'
        else:                          col2 = 'ok'

        col3 = 'ok' if len(extradoc_opts) == 0 and len(extratest_opts) == 0 else '[+warnings]'

        print(f'{col1:18s} {col2:18s} {col3:12s}')

    else:
        is_ok = True
        if manfile_optset is None:
            print(f'No man page:\n    failed to find or open man page {manfile}\n')
            is_ok = False
        if itestfile_optset is None:
            print(f'No itest script:\n    failed to find or open integrated test script {itestfile}\n')
            is_ok = False

        if len(undoc_opts) > 0:
            print(f'{len(undoc_opts)} options not documented in man page:')
            for opt in undoc_opts: print(f'    {opt}')
            print('')
            is_ok = False

        if len(untest_opts) > 0:
            print(f'{len(untest_opts)} options not tested in integrated test script:')
            for opt in untest_opts: print(f'    {opt}')
            print('')
            is_ok = False

        if len(extradoc_opts) > 0:
            print(f"Warning: {len(extradoc_opts)} options documented that aren't in .c file:")
            for opt in extradoc_opts: print(f'    {opt}')
            print('')
            is_ok = False

        if len(extratest_opts) > 0:
            print(f"Warning: {len(extratest_opts)} options tested that aren't in .c file:")
            for opt in extratest_opts: print(f'    {opt}')
            print('')
            is_ok = False

        if is_ok: print("ok")

if __name__ == "__main__":
    main()
