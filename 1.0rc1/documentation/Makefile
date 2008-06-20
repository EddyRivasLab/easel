AUTODOC  = ../devkit/autodoc
CEXCERPT = ../devkit/cexcerpt


# Each module in this list obeys three documentation conventions:
#
#  1. It has a .tex file documenting its API.
#     These are handwritten, in the source directory with
#     the .c and .h files.
#
#  2. It can be processed with autodoc_functions to
#     produce the function documentation appendix to 
#     each module. These appendices go to the 
#     autotext/ subdir.
#
#  3. It can be processed with cexcerpts to extract 
#     verbatim code fragments for inclusion in the 
#     documentation. These fragments go to the cexcerpts/
#     subdir.
#
MODULES = easel\
	  esl_alphabet\
	  esl_cluster\
	  esl_dirichlet\
	  esl_distance\
	  esl_dmatrix\
	  esl_exponential\
	  esl_fileparser\
	  esl_gamma\
	  esl_getopts\
	  esl_gev\
	  esl_gumbel\
	  esl_histogram\
	  esl_hyperexp\
	  esl_keyhash\
	  esl_minimizer\
	  esl_mixgev\
	  esl_mpi\
	  esl_msa\
	  esl_msacluster\
	  esl_msaweight\
	  esl_normal\
	  esl_paml\
	  esl_random\
	  esl_randomseq\
	  esl_ratematrix\
	  esl_regexp\
	  esl_rootfinder\
	  esl_scorematrix\
	  esl_sq\
	  esl_sqio\
	  esl_sse\
	  esl_ssi\
	  esl_stack\
	  esl_stats\
	  esl_stopwatch\
	  esl_stretchexp\
	  esl_swat\
	  esl_tree\
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
	mv main.pdf current.pdf

autodoc:
	mkdir -p autotext
	for module in ${MODULES}; do\
	   ${AUTODOC} ../$$module.c >  autotext/$${module}_functions.tex;\
	done

cexcerpt:
	mkdir -p cexcerpts
	for module in ${MODULES}; do\
	   ${CEXCERPT} ../$$module.c cexcerpts;\
	   ${CEXCERPT} ../$$module.h cexcerpts;\
	done

clean:
	-rm *~ TAGS
	-rm main.log main.out main.pdf main.aux main.toc 
	-rm main.bbl main.blg main.brf
	-rm typography.log typography.aux typography.out
	-rm autotext/*.tex

distclean:
	-rm -rf cexcerpts
	-rm -rf autotext
	make clean
