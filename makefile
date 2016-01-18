FC=ifort
#FC=gfortran
FFLAG=-traceback -g -check all
test_read_aoints: test_read_aoints.f90 read_aoints.f90
	${FC} ${FFLAG} test_read_aoints.f90 read_aoints.f90 -o $@

test_read_intin: test_read_intin.f90 read_intin.f90
	${FC} ${FFLAG} $^ -o $@
check_read_intin: test_read_intin
	cd test/out3 && ../../test_read_intin < int.in 


check: test_read_aoints
	cd test/out3 && ../../test_read_aoints

clean:
	rm *.mod
