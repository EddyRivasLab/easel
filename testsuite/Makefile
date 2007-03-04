
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

dcheck:
	@echo    Running developer testsuite...
	@echo -n Testing that all drivers compile...
	@./driver_report.pl       >  driver.report
	@echo    done. [See driver.report]
	@echo -n Measuring test code coverage...
	@./coverage_report.pl     >  coverage.report
	@echo    done. [See coverage.report]
	@echo -n Looking for memory leaks...
	@./valgrind_report.pl     >  valgrind.report
	@echo    done. [See valgrind.report]
	@echo -n Running unit tests...
	@./sqc 2 testsuite.sqc .  >  utest.report
	@echo    done. [See utest.report]
	@cat driver.report coverage.report valgrind.report utest.report > dcheck.report
	@echo
	@echo    Finished developer testsuite.
	@echo    For a complete report, see dcheck.report.

clean:
	-rm -f esl[^_]?????
	-rm -f *_utest
	-rm -f *.gcno *.gcda *.gcov
	-rm -f *.o *~ 
	-rm -f driver.report coverage.report valgrind.report utest.report
	-rm -f ${PROGS} 

distclean:
	make clean
