CC      = gcc
CFLAGS  = -g
AR      = ar rcv
LN      = ln
RANLIB  = ranlib

HEADERS = \
	bioparse_paml.h\
	dirichlet.h\
	dmatrix.h\
	easel.h\
	gamma.h\
	interface_gsl.h\
	interface_lapack.h\
	parse.h\
	random.h\
	ratematrix.h\
	vectorops.h

OBJS    = \
	bioparse_paml.o\
	dirichlet.o\
	dmatrix.o\
	easel.o\
	gamma.o\
	interface_gsl.o\
	interface_lapack.o\
	parse.o\
	random.o\
	ratematrix.o\
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
	-rm -f ${OBJS} *~