CC      = gcc
CFLAGS  = -g -Wall
AR      = ar rcv
LN      = ln
RANLIB  = ranlib

HEADERS = \
	alphabet.h\
	bioparse_paml.h\
	dirichlet.h\
	dmatrix.h\
	easel.h\
	gamma.h\
	getopts.h\
	interface_gsl.h\
	interface_lapack.h\
	keyhash.h\
	msa.h\
	parse.h\
	random.h\
	ratematrix.h\
	regexp.h\
	sqio.h\
	stack.h\
	vectorops.h

OBJS    = \
	alphabet.o\
	bioparse_paml.o\
	dirichlet.o\
	dmatrix.o\
	easel.o\
	gamma.o\
	getopts.o\
	interface_gsl.o\
	interface_lapack.o\
	keyhash.o\
	msa.o\
	parse.o\
	random.o\
	ratematrix.o\
	regexp.o\
	sqio.o\
	stack.o\
	vectorops.o

all: libeasel.a

.c.o:
	${CC} -I. ${CFLAGS} ${DEFS} -c $<		

libeasel.a: $(OBJS)
	$(AR) libeasel.a $(OBJS)
	$(RANLIB) libeasel.a
	chmod 644 libeasel.a

symlinks:
	mkdir -p easel
	for header in ${HEADERS}; do\
	   (cd easel; ${LN} -s ../$$header .);\
	done

clean:
	-rm -f ${OBJS} *~ libeasel.a


reallyclean:
	make clean
	(cd documentation; make clean)
	-rm easel/*.h

tags:
	etags *.[ch] Makefile



# magic SVN for setting keyword ID replacement on a new module foo:
# svn propset svn:keywords "Id" foo.[ch]