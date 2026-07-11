#! /usr/bin/env python3

# pandoc_postprocess.py
#
# Massages output of pandoc man => latex conversion to fit Tufteian
# userguide style.
#
# Example:
#    pandoc -f man -t latex --columns=80 hmmbuild.man | pandoc_postprocess.py > manpage.tex
#
# - first \section{NAME} converted to \section{progname - description}
# - remaining \section's converted to \subsection* (unnumbered)
# - if SYNOPSIS section has multiple example commands, indent them equally
# - change {description} environment to our customized {wideitem}
# - convert \item[\textbf  to \item[\monob
# - convert \textbf{% ...} to \user{% ...}
# - convert \emph to \monoi
# - convert \textbf to \mono  : literals are .B in .man, mono in .tex
#
# \mono, \monoi, \monob, \user, and {wideitems} environment are
# defined in tufte-common-local.tex customization of user guides.
#
import re
import sys

if len(sys.argv) == 1: f = sys.stdin
else:                  f = open(sys.argv[1])

for line in f:
    line = line.rstrip('\n')

    # State flags (where are we in the document)
    if re.match(r'\\section{SYNOPSIS}', line):
        in_synopsis = True
    elif re.match(r'\\section{', line):
        in_synopsis = False

    # Replace \section{NAME} with \section{progname - description}, using next line too.
    if re.match(r'\\section\{NAME\}', line):
        for line in f:
            if not re.fullmatch(r'\s*', line):
                break
        m = re.match(r'((?:easel )?\S+)\s*\\?-\s*(.+)$', line) 
        if m:
            print(r'\section{{\texorpdfstring{{\monob{{{0}}}}}{{{0}}} - {1}}}'.format(m.group(1), m.group(2)))
        else:
            print("Error: no progname/description line found");
            sys.exit(1)
        continue

    # Remove everything after \section{See Also), and finish.
    if re.match(r'\\section\{SEE ALSO', line) or re.match(r'\\end\{document', line):
        print("\\newpage");
        break

    # In synopsis, put \noindent in front of each commandline, and preserve the .B's as bold.
    if in_synopsis and re.match(r'\s*\\textbf{', line):
        line = re.sub(r'\\textbf\{', r'\\monob{', line)
        print("\\noindent")

    #
    # Substitutions within a line.
    # The order of these replacements is important. (More specific ones first.)
    #
    line = re.sub(r'\\begin\{description\}',  r'\\begin{wideitem}', line)
    line = re.sub(r'\\end\{description\}',    r'\\end{wideitem}',   line)
    line = re.sub(r'\\section\{',             r'\\subsection*{',    line)   # \subsection* suppresses inclusion in TOC
    line = re.sub(r'\\item\s*\[\\textbf',     r'\\item [\\monob',   line)   # option names in .TP are emphasized bold
    line = re.sub(r'\\textbf\{\\% ',          r'\\user{\\% ',       line)   # example command lines are bold, on their own line
    line = re.sub(r'\\emph\{',                r'\\monoi{',          line)   # metavariables (options, args) are .I in man, mono italic in tex
    line = re.sub(r'\\textbf\{',              r'\\mono{',           line)   # literals (commands, etc) are .B in man, normal mono in tex

    print (line)
    

if f != sys.stdin:
    f.close()

