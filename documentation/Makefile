# Each module should have .tex documentation of its API,
# in addition to the .c and .h files. The list of modules
# that actually comply follows:
#
APITEX =  alphabet\
	  easel\
	  getopts\
	  keyhash\
	  msa\
	  random\
	  regexp\
	  sqio\
	  stack\
	  wuss


# Each module should be able to be autodocumented, converting
# its function headers to a .tex appendix file. The list of
# modules that actually comply follows:
#
AUTOTEX = alphabet\
	  dmatrix\
	  getopts\
	  keyhash\
	  msa\
	  random\
	  regexp\
	  sqio\
	  stack\
	  wuss

pdf:
	make autodoc
	pdflatex main
	bibtex   main
	pdflatex main
	pdflatex main

autodoc:
	for module in ${AUTOTEX}; do\
	   ./autodoc_functions ../$$module.c > autotext/$${module}_functions.tex;\
	done

clean:
	-rm *~ TAGS
	-rm main.log main.out main.pdf main.aux main.toc 
	-rm main.bbl main.blg main.brf
	-rm autotext/*.tex

distclean:
	make clean
