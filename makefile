#FC=ifort
FC=gfortran
#FFLAG=-traceback -g -check all
FFLAG=-Wall -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace -g
MODS=-I${HOME}/local/mod
LIBLAPACK= -lblas -llapack -llapack95

%.o: %.f90
	${FC} ${FFLAG} ${MODS} -c $< -o $@

vector.f90: vector.f90t
	cat vector.f90t| sed s/XXXTypeXXX/integer/g | sed s/XXXNameXXX/I/g >  vector.f90
	cat vector.f90t| sed s/XXXTypeXXX/complex*16/g | sed s/XXXNameXXX/Z/g >> vector.f90
	cat vector.f90t| sed s/XXXTypeXXX/real*8/g | sed s/XXXNameXXX/D/g >> vector.f90

# ==== Binary ====
proj_one: utils.o lalgebra.o read_aoints.o read_mocoef.o read_ci.o proj_one.o run_proj_one.o
	${FC} ${FFLAG} $^ ${LIBLAPACK} -o $@
check_proj_one: proj_one
	cd test/out/ && ../../proj_one proj_one.in

solve_1e: utils.o lalgebra.o read_aoints.o solve_1e.o
	${FC} ${FFLAG} $^ ${LIBLAPACK} -o $@
run_solve_1e: solve_1e
	cd test/out2 && ../../solve_1e solve_1e.in

gto_vec: utils.o read_intin.o gto_vec.o 
	${FC} ${FFLAG} $^ ${LIBLAPACK} -o $@
run_gto_vec: gto_vec
	cd test/out && ../../gto_vec int.in gto_vec_delta.in  && cat gto_vec_delta.dat


# ==== Unit Tests ====
# ---- general ----
test_block: utils.o utest.o block.o test_block.o
	${FC} ${FFLAG} $^ -o $@
check_block: test_block
	./test_block

test_vector: utest.o vector.o test_vector.o
	${FC} ${FFLAG} $^ -o $@
	./test_vector

test_utils: utest.o utils.o lalgebra.o test_utils.o
	${FC} ${FFLAG} $^ ${LIBLAPACK} -o $@
check_utils: test_utils
	./test_utils

# ---- reader ----
test_read_intin: utils.o utest.o read_intin.o test_read_intin.o
	${FC} ${FFLAG} $^ -o $@
check_read_intin: test_read_intin
	cd test/out3 && ../../test_read_intin < int.in 

test_read_aoints:  utest.o block.o read_aoints.o test_read_aoints.o
	${FC} ${FFLAG} $^ ${LIBLAPACK} -o $@
check_read_aoints: test_read_aoints
	cd test/out3 && ../../$< 

test_read_mocoef:  utest.o block.o utils.o read_intin.o read_aoints.o read_mocoef.o test_read_mocoef.o
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
