# Each module should have .tex documentation of its API,
# in addition to the .c and .h files. The list of modules
# that actually comply follows:
#
APITEX =  alphabet\
	  easel\
	  random\
	  regexp\
	  sqio\
	  stack


# Each module should be able to be autodocumented, converting
# its function headers to a .tex appendix file. The list of
# modules that actually comply follows:
#
AUTOTEX = alphabet\
	  dmatrix\
	  random\
	  regexp\
	  sqio\
	  stack

pdf:
	make autodoc
	pdflatex main
	pdflatex main

autodoc:
	for module in ${AUTOTEX}; do\
	   ./autodoc_functions ../$$module.c > autotext/$${module}_functions.tex;\
	done

clean:
	-rm *~ main.log main.out main.pdf main.aux main.toc
	-rm autotext/*.tex
