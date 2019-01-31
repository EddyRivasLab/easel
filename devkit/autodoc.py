#! /usr/bin/env python3

# autodoc.py
#
# Extract function documentation from a .c file, and output it in a
# Markdown file. Assumes that the .c file uses function headers that
# conform to Easel code style.
# 
# Usage:    ./autodoc.py <C file>
# Example:  ./autodoc.py esl_module.c > esl_module_functions.md
#
#
# Replaces old Perl autodoc, which converted to LaTeX. I'm moving
# toward online developer documentation, and Markdown. This simplifies
# the conversion; we don't have to protect special LaTeX characters
# (#_$) for example.
#
# SRE, Sun 27 Jan 2019

import sys
import re

def process(text):
    """
    Remove the leading " *     " prefixes from a Purpose, Returns, or Throws
    multiiline block of text that we've just pulled out of the function header.
    
    Convert <text> to `text` (Markdown code).
    """
    text = re.sub(r'^ \*[ \t]+(\S)', r'\1',   text, flags=re.MULTILINE)
    text = re.sub(r'^ \*',           r'',     text, flags=re.MULTILINE)
    text = re.sub(r'<(\S|\S.*?\S)>', r'`\1`', text, flags=re.MULTILINE | re.DOTALL)
    return text

def output_argtable(argtext):
    print("|  arg | description |")
    print("|------|-------------|")
    for m in re.finditer(r'^(\S+)\s*[:-]\s*(.+?)\s*$', argtext, flags=re.MULTILINE):
        print("| `{0}` | {1} |".format(m.group(1), m.group(2)))
    print("\n")

def main():
    cfile = sys.argv[1]
    fp    = open(cfile)
    text  = fp.read()

    for m in re.finditer(r'^/\*\s+Function:.+?^ \*/.+?\{', text, flags=re.MULTILINE | re.DOTALL):
        header = m.group(0)      # m.group(0) is now: `/* Function: ... */ int funcname(...) {`

        m = re.search(r'/\*\s+Function:\s+(\S+)', header)
        funcname = m.group(1)

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

        m = re.search(r'\*/\s*(.+?)\s*\{$', header, flags=re.DOTALL)
        syntax = re.sub(r'\s+', r' ', m.group(1))     # collapse extra whitespace including newlines

        print("### `{0}`\n".format(funcname))
        if synopsis: print("**{0}**\n".format(synopsis)) 
        print("`{0}`\n".format(syntax))
        if argtext: output_argtable(argtext)
        if purpose: print(purpose);
        if returns: print("Returns: {0}".format(returns))
        if throws:  print("Throws: {0}".format(throws))
        print("------")



if __name__ == "__main__":
    main()







    




