

pdf:
	make autodoc
	pdflatex main
	pdflatex main

autodoc:
	./autodoc_functions ../dmatrix.c > autotext/dmatrix_functions.tex

clean:
	\rm *~ main.log main.out main.pdf
