#FC=ifort
FC=gfortran
#FFLAG=-traceback -g -check all
FFLAG=-Wall -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace -g
MODS=-I${HOME}/local/mod
LIBLAPACK= -lblas -llapack -llapack95

%.o: %.f90
	${FC} ${FFLAG} ${MODS} -c $< -o $@

# ==== Binary ====
solve_1e: utils.o read_aoints.o solve_1e.o
	${FC} ${FFLAG} $^ ${LIBLAPACK} -o $@
run_solve_1e: solve_1e
	cd test/out && ../../solve_1e

gto_vec: utils.o read_intin.o gto_vec.o 
	${FC} ${FFLAG} $^ ${LIBLAPACK} -o $@
run_gto_vec: gto_vec
	cd test/out && ../../gto_vec int.in gto_vec_delta.in  && cat gto_vec_delta.dat


# ==== Unit Tests ====
test_utils: utest.o utils.o test_utils.o
	${FC} ${FFLAG} $^ ${LIBLAPACK} -o $@
check_utils: test_utils
	./test_utils

test_read_aoints:  utest.o read_aoints.o test_read_aoints.o
	${FC} ${FFLAG} $^ ${LIBLAPACK} -o $@
check_read_aoints: test_read_aoints
	cd test/out3 && ../../$< 

test_read_intin: utest.o utils.o read_intin.o test_read_intin.o
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

# ==== other ====

check: test_read_aoints
	cd test/out3 && ../../test_read_aoints

clean:
	rm -f *.mod
	rm -f *.o
	rm -f test_read_ci test_read_intin test_read_aoints
