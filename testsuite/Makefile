################################################################
# Makefile for the Easel testsuite
# Makefile is generated from Makefile.in by autoconf
# DO NOT EDIT Makefile: only edit Makefile.in
# 
# SRE, Mon Jun 11 08:19:32 2007 [Janelia]
# SVN $Id$

CC       = gcc -std=gnu99
CFLAGS   = -g -O2
CPPFLAGS = 
LDFLAGS  = 
DEFS     = -DHAVE_CONFIG_H
LIBS     =  -lm
MPILIBS  = 

AR      = /usr/bin/ar rcv
RANLIB  = ranlib

ESLDIR = ..

PROGS = parse-t1\
	dmatrix-t1\
	dmatrix-t2\

UNIT_TESTS = easel_utest\
	alphabet_utest\
	dirichlet_utest\
	distance_utest\
	dmatrix_utest\
	exponential_utest\
	fileparser_utest\
	gamma_utest\
	getopts_utest\
	gumbel_utest\
	histogram_utest\
	hyperexp_utest\
	keyhash_utest\
	msa_utest\
	msaweight_utest\
	random_utest\
	ratematrix_utest\
	regexp_utest\
	rootfinder_utest\
	scorematrix_utest\
	sqio_utest\
	ssi_utest\
	stack_utest\
	stats_utest\
	stretchexp_utest\
	tree_utest\
	vectorops_utest\
	weibull_utest\
	wuss_utest
#	gev_utest\
#	minimizer_utest\
#	mixgev_utest\
#	mpi_utest\
#	msacluster_utest\
#	normal_utest\
#	paml_utest\
#	stopwatch_utest\

#	swat_utest\


all:  ${UNIT_TESTS}

.c.o:
	${CC} -I ${EASELLIB} ${CFLAGS} ${DEFS} -c $<		

$(PROGS): %: %.o ../libeasel.a
	${CC} ${CFLAGS} -o $@ -L ${EASELLIB} $@.o -leasel -lm

# Unit test compilation:
# 
# To compile unit tests, we need to construct a command line
# from module basename, module's testdriver flag, and module's
# C file. This is accomplished via a complex-looking shell scripture;
# xref J1/81.
#   $@           =  driver name            msa_utest
#   $${BASENAME} =  base module name       msa
#   $${DFLAG}    =  driver flag            eslMSA_TESTDRIVE
#   $${DFILE}    =  driver's C file        esl_msa.c
#
$(UNIT_TESTS): FORCE
	@echo Compiling $@ ...
	@BASENAME=`echo $@ | sed -e 's/_utest//'` ;\
	DFLAG=`echo $${BASENAME} | sed -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`;\
	DFLAG=esl$${DFLAG}_TESTDRIVE ;\
	if test $@ = "easel_utest"; then \
           DFILE=easel.c ;\
        else \
           DFILE=esl_$${BASENAME}.c ;\
        fi;\
	${CC} ${CFLAGS} ${DEFS} ${LDFLAGS} -o $@ -I${ESLDIR} -L${ESLDIR} -D$${DFLAG} ${ESLDIR}/$${DFILE} ${MPILIBS} -leasel ${LIBS}

FORCE:


check:  ${UNIT_TESTS}
	@echo 
	@echo Running Easel test suite...
	@echo
	@../devkit/sqc 2 testsuite.sqc .

dcheck: ${UNIT_TESTS}
	@echo    Running developer testsuite...
	@echo -n Running unit tests...
	@../devkit/sqc 2 testsuite.sqc .  >  utest.report
	@echo    done. [See utest.report]
	@echo -n Testing that all drivers compile...
	@env CC="${CC}" CFLAGS="${CFLAGS}" ./driver_report.pl       >  driver.report
	@echo    done. [See driver.report]
	@echo -n Measuring test code coverage...
	@env CC="${CC}" CFLAGS="${CFLAGS}" ./coverage_report.pl     >  coverage.report
	@echo    done. [See coverage.report]
	@echo -n Looking for memory leaks...
	@env CC="${CC}" CFLAGS="${CFLAGS}" ./valgrind_report.pl     >  valgrind.report
	@echo    done. [See valgrind.report]
	@cat driver.report coverage.report valgrind.report utest.report > dcheck.report
	@echo
	@echo    Finished developer testsuite.
	@echo    For a complete report, see dcheck.report.

clean:
	-rm -f esl[^_]?????
	-rm -f *.gcno *.gcda *.gcov
	-rm -f *.o *~ 
	-rm -f driver.report coverage.report valgrind.report utest.report
	-rm -f ${UNIT_TESTS}
	-rm -f ${PROGS} 

distclean:
	make clean
