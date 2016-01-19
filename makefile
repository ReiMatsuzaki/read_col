FC=ifort
#FC=gfortran
FFLAG=-traceback -g -check all

%.o: %.f90
	${FC} ${FFLAG} -c $< -o $@

test_read_aoints:  read_aoints.o test_read_aoints.o
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

check: test_read_aoints
	cd test/out3 && ../../test_read_aoints

clean:
	rm *.mod
