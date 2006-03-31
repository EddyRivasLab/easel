# Each module should have .tex documentation of its API,
# in addition to the .c and .h files. The list of modules
# that actually comply follows:
#
APITEX =  easel\
	  esl_alphabet\
	  esl_dirichlet\
	  esl_exponential\
	  esl_fileparser\
	  esl_gamma\
	  esl_getopts\
	  esl_gev\
	  esl_gumbel\
	  esl_histogram\
	  esl_hyperexp\
	  esl_keyhash\
#	  esl_mixdchlet\
	  esl_mixgev\
	  esl_msa\
	  esl_random\
	  esl_regexp\
	  esl_sqio\
	  esl_ssi\
	  esl_stack\
	  esl_stopwatch\
	  esl_stretchexp\
	  esl_vectorops\
	  esl_weibull\
	  esl_wuss

# Each module should be able to be autodocumented, converting
# its function headers to a .tex appendix file. The list of
# modules that actually comply follows:
#
AUTOTEX = easel\
	  esl_alphabet\
	  esl_dirichlet\
	  esl_exponential\
	  esl_fileparser\
	  esl_gamma\
	  esl_getopts\
	  esl_gev\
	  esl_gumbel\
	  esl_histogram\
	  esl_hyperexp\
	  esl_keyhash\
	  esl_mixgev\
	  esl_msa\
	  esl_random\
	  esl_regexp\
	  esl_sqio\
	  esl_ssi\
	  esl_stack\
	  esl_stopwatch\
	  esl_stretchexp\
	  esl_vectorops\
	  esl_weibull\
	  esl_wuss

pdf:
	make autodoc
	make cexcerpt
	pdflatex main
	bibtex   main
	pdflatex main
	pdflatex main

autodoc:
	mkdir -p autotext
	for module in ${AUTOTEX}; do\
	   ./autodoc_functions ../$$module.c >  autotext/$${module}_functions.tex;\
	done

cexcerpt:
	mkdir -p cexcerpts
	for module in ${AUTOTEX}; do\
	   ~/src/ssdk/cexcerpt ../$$module.c cexcerpts;\
	   ~/src/ssdk/cexcerpt ../$$module.h cexcerpts;\
	done

clean:
	-rm *~ TAGS
	-rm main.log main.out main.pdf main.aux main.toc 
	-rm main.bbl main.blg main.brf
	-rm autotext/*.tex

distclean:
	-rm -rf cexcerpts
	-rm -rf autotext
	make clean
