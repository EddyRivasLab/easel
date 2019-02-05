#! /usr/bin/env python3

# autodoc.py
#
# Extract function documentation from a .c file and output it in a
# Markdown file. Assumes that the .c file uses function headers that
# conform to Easel code style. 
# 
# Usage:    ./autodoc.py <.c file>
# Example:  ./autodoc.py esl_module.c > esl_module_functions.md
#
#
# Replaces old Perl autodoc, which converted to LaTeX. I'm moving
# toward online developer documentation, and Markdown. This simplifies
# the conversion; we don't have to protect special LaTeX characters
# (#_$) for example.
#
# SRE, Sun 27 Jan 2019

import os
import sys
import re

def process(text):
    """
    Remove the leading " *     " prefixes from a Purpose, Returns, or Throws
    multiiline block of text that we've just pulled out of the function header.
    
    Convert <text> to `text` (Markdown code).
    """
    text = re.sub(r'^ \*[ \t]+(\S)', r'\1',   text, flags=re.MULTILINE)               # remove leading " *   "
    text = re.sub(r'^ \*',           r'',     text, flags=re.MULTILINE)               # remove " *" alone
    text = re.sub(r'<(\S|\S.*?\S)>', r'`\1`', text, flags=re.MULTILINE | re.DOTALL)   # convert <text> to `text`
    return text

def output_argtable(argtext):
    print("|  arg | description |")
    print("|------|-------------|")
    for m in re.finditer(r'^(\S+)\s*[:-]\s*(.+?)\s*$', argtext, flags=re.MULTILINE):
        print("| `{0}` | {1} |".format(m.group(1), m.group(2)))
    print("\n")

def main():
    cfile = sys.argv[1]
    if not os.path.isfile(cfile): exit(".c file {0} not found".format(cfile))
    fp   = open(cfile)
    text = fp.read()

    #                       /* Function:  ...    */ ... }   blank line. Grabs header + implementation(s).
    #                       vv    vvvvvvvv       vv   vvv   vvv   v...or, grabs subheading line in a comment
    for m in re.finditer(r'^(/\*\s+Function:.+?^ \*/)(.+?^\})\s*$^\s*$|^\s*\*#\s*\d+\..+?$\s*', text, flags=re.MULTILINE | re.DOTALL):
        if m.group(0).startswith('/*'):
            header = m.group(1)    # comment header "/* Function: ... */"
            impl   = m.group(2)    # implementation(s) 

            m = re.match(r'/\*\s+Function:\s*(.+)$', header, flags=re.MULTILINE)    # Usually one function name, but could also be comma-delimited list    
            funcnames = [ a.lstrip().rstrip() for a in m.group(1).split(',') ]    

            m = re.search(r'^\s+\*\s+Synopsis:\s+(.+)$', header, flags=re.MULTILINE)
            synopsis = m.group(1) if m else None

            m = re.search(r'^\s+\*\s+Args:\s+(.+?)(?:^ \*/|^ \* \S)', header, flags=re.MULTILINE | re.DOTALL)
            argtext = process(m.group(1)) if m else None

            m = re.search(r'^\s+\*\s+Purpose:\s+(.+?)(?:^ \*/|^ \* \S)', header, flags=re.MULTILINE | re.DOTALL)
            purpose = process(m.group(1)) if m else None

            m = re.search(r'^\s+\*\s+Returns:\s+(.+?)(?:^ \*/|^ \* \S)', header, flags=re.MULTILINE | re.DOTALL)
            returns = process(m.group(1)) if m else None

            m = re.search(r'^\s+\*\s+Throws:\s+(.+?)(?:^ \*/|^ \* \S)', header, flags=re.MULTILINE | re.DOTALL)
            throws = process(m.group(1)) if m else None

            # pull the syntax out of the C implementation.
            # nontrivial to do well with just regexps, without a real grammar parser,
            # because we're covering the less common case where there's >1 function
            # documented by a single header.
            syntax = []
            for fname in funcnames:                                 # list of names like "esl_foo_Function()", with the ()
                fname = fname.rstrip('()')                          # now just "esl_foo_Function"
                pattern1 = r'^(\S+.+?)\s+' + fname + '\s*(\((?s:.+?)\))\s*\{'
                m = re.search(pattern1, impl, flags=re.MULTILINE)
                if m: syntax.append(m.group(1) + ' ' + fname + m.group(2))
                else: exit("failed to parse out the syntax for {0}".format(fname))

            for a in funcnames: print("### `{0}`\n".format(a))
            if synopsis: print("**{0}**\n".format(synopsis.rstrip())) 
            for s in syntax:    print("`{0}`\n".format(s))
            if argtext: output_argtable(argtext)
            if purpose: print(purpose);
            if returns: print("Returns: {0}".format(returns))
            if throws:  print("Throws: {0}".format(throws))
            print("------")

        else:  # or, we're a section heading.
            m = re.match('^\s*\*#\s*(\d+\..+)', m.group(0))
            secheading = m.group(1)
            print("## {0}\n".format(secheading))


if __name__ == "__main__":
    main()







    




