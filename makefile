FC=ifort
#FC=gfortran
FFLAG=-traceback -g -check all
read_aoints: test_read_aoints.f90 read_aoints.f90
	${FC} ${FFLAG} test_read_aoints.f90 read_aoints.f90 -o read_aoints

check: read_aoints
	cd test/out3 && ../../read_aoints

clean:
	rm *.mod
