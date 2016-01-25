#FC=ifort
FC=gfortran
#FFLAG=-traceback -g -check all
FFLAG=-Wall -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace -g

%.o: %.f90
	${FC} ${FFLAG} -c $< -o $@

test_read_aoints:  utest.o read_aoints.o test_read_aoints.o
	${FC} ${FFLAG} $^ -o $@
check_read_aoints: test_read_aoints
	cd test/out3 && ../../$< 

test_read_intin: test_read_intin.f90 read_intin.f90
	${FC} ${FFLAG} $^ -o $@
check_read_intin: test_read_intin
	cd test/out3 && ../../test_read_intin < int.in 

test_read_mocoef:  read_mocoef.o read_aoints.o test_read_mocoef.o
	${FC} ${FFLAG} $^ -o $@
check_read_mocoef: test_read_mocoef
	cd test/out && ../../test_read_mocoef

test_read_ci: utest.o read_ci.o test_read_ci.o 
	${FC} ${FFLAG} $^ -o $@
check_read_ci: test_read_ci
	cd test/out && ../../test_read_ci < ci2.in

check: test_read_aoints
	cd test/out3 && ../../test_read_aoints

clean:
	rm -f *.mod
	rm -f *.o
	rm -f test_read_ci test_read_intin test_read_aoints
