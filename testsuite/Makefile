
CC      = gcc
CFLAGS  = -g
AR      = ar rcv
RANLIB  = ranlib

EASELLIB = ../
EASELINC = ../

PROGS = parse-t1\
	dmatrix-t1\
	dmatrix-t2\

ALL:  ${PROGS}

.c.o:
	${CC} -I ${EASELLIB} ${CFLAGS} ${DEFS} -c $<		

$(PROGS): %: %.o ../libeasel.a
	cc ${CFLAGS} -o $@ -L ${EASELLIB} $@.o -leasel -lm

check:
	./sqc 2 testsuite.sqc .

clean:
	-rm -f testdrive-*
	-rm -f *.o *~ ${PROGS} testdriver

distclean:
	make clean
