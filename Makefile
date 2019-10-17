# source files
SRCFILES = src/checkpoint.c src/dupio.c src/greens_func.c src/hash_table.c src/kinetic.c src/lattice_fourier.c src/linalg.c src/measurement.c src/monte_carlo.c src/param_parser.c src/profiler.c src/progress.c src/random.c src/stratonovich.c src/time_flow.c src/util.c
# test files
TSTFILES = test/matrix_exp_test.c test/blockcyclic_qr_test.c test/blockcyclic_tri_test.c test/blockcyclic_inv_test.c test/kinetic_test.c test/kinetic_test2.c test/kinetic_test3.c test/kinetic_test4.c test/lattice_fourier_test.c test/time_flow_test1.c test/time_flow_test2.c test/time_flow_test3.c test/stratonovich_test.c test/greens_func_flip_test.c test/greens_func_wrap_test.c test/greens_func_init_test1.c test/greens_func_init_test2.c test/greens_func_init_test3.c test/greens_func_init_test4.c test/monte_carlo_iter_test.c test/monte_carlo_phonon_block_test.c test/monte_carlo_iter_phonon_test.c test/green_unequal_time_test.c test/measurement_test.c

# This configuration uses gcc with generic CBLAS and LAPACKE libraries, and enables OpenMP parallelization.
CC = gcc
# compiler options
CCOPTS = -Wall -Iinclude -O2 -DVERSION=\"$(shell git describe --always)\"
CCOPTS += -fopenmp
# set these with appropriate libraries for your system
LIBRARIES = -lm -lblas -llapacke

## The following configuration selects the Intel compiler with MKL, and enables OpenMP parallelization.
#CC = icc
## compiler options
#CCOPTS = -Wall -Iinclude -O2 -xHost -restrict -DUSE_MKL -DMKL_DIRECT_CALL -mkl:sequential -DVERSION=\"$(shell git describe --always)\"
## replace by -qopenmp-stubs to disable OpenMP shared-memory parallelization
#CCOPTS += -qopenmp
## set these with appropriate libraries for your system
#LIBRARIES = -mkl:sequential -lrt

# static libraries, useful for gprof
#MKLROOT = /opt/intel/mkl
#LIBRARIES = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread

all: proj_main proj_test

proj_main: ${SRCFILES} src/main.c
	${CC} ${CCOPTS} -DNDEBUG -DPROFILE_ENABLE -o bin/hubbard_dqmc $? ${LIBRARIES}

proj_test: ${SRCFILES} ${TSTFILES} test/run_tests.c
	${CC} ${CCOPTS} -g -o test/run_tests $? ${LIBRARIES}
