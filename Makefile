CC      = gcc
CFLAGS  = -g
AR      = ar rcv
RANLIB  = ranlib

HEADERS = \
	bioparse_paml.h\
	dmatrix.h\
	easel.h\
	interface_gsl.h\
	interface_lapack.h\
	parse.h\
	random.h\
	ratematrix.h\
	vectorops.h

OBJS    = \
	bioparse_paml.o\
	dmatrix.o\
	easel.o\
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

clean:
	-rm -f ${OBJS} *~