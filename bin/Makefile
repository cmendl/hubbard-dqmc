SRCFILES = checkpoint.o dupio.o greens_func.o hash_table.o kinetic.o lattice_fourier.o linalg.o measurement.o monte_carlo.o param_parser.o profiler.o progress.o random.o stratonovich.o time_flow.o util.o
TSTFILES = matrix_exp_test.o blockcyclic_qr_test.o blockcyclic_tri_test.o blockcyclic_inv_test.o kinetic_test.o kinetic_test2.o lattice_fourier_test.o time_flow_test1.o time_flow_test2.o time_flow_test3.o stratonovich_test.o greens_func_flip_test.o greens_func_wrap_test.o greens_func_init_test1.o greens_func_init_test2.o greens_func_init_test3.o greens_func_init_test4.o monte_carlo_iter_test.o monte_carlo_phonon_block_test.o monte_carlo_iter_phonon_test.o green_unequal_time_test.o measurement_test.o

# use Intel compiler
CC = icc

# compiler options
CCOPTS = -Wall -O2 -xHost -restrict -DNDEBUG -DMEM_DATA_ALIGN=64 -DVERSION=\"$(shell git describe --always)\"
# uncomment to enable profiling; make sure to run 'make clean' before recompiling
#CCOPTS += -DPROFILE_ENABLE
# replace by -qopenmp to enable OpenMP shared-memory parallelization
CCOPTS += -qopenmp-stubs

# set these with appropriate libraries for your system
MKL = -mkl=sequential
LIBRARIES = ${MKL} -lrt

# static libraries, useful for gprof
#MKLROOT = /opt/intel/mkl
#LIBRARIES = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread

all: proj_main proj_test

proj_main: ${SRCFILES} main.o
	${CC} ${CCOPTS} -o hubbard_dqmc $? ${LIBRARIES}

proj_test: ${SRCFILES} ${TSTFILES} run_tests.o
	${CC} ${CCOPTS} -o run_tests $? ${LIBRARIES}


%.o: ../src/%.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<


run_tests.o: ../test/run_tests.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<


matrix_det_test.o: ../test/matrix_det_test.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

matrix_exp_test.o: ../test/matrix_exp_test.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

blockcyclic_qr_test.o: ../test/blockcyclic_qr_test.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

blockcyclic_tri_test.o: ../test/blockcyclic_tri_test.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

blockcyclic_inv_test.o: ../test/blockcyclic_inv_test.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

kinetic_test.o: ../test/kinetic_test.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

kinetic_test2.o: ../test/kinetic_test2.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

lattice_fourier_test.o: ../test/lattice_fourier_test.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

time_flow_test1.o: ../test/time_flow_test1.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

time_flow_test2.o: ../test/time_flow_test2.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

time_flow_test3.o: ../test/time_flow_test3.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

stratonovich_test.o: ../test/stratonovich_test.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

greens_func_flip_test.o: ../test/greens_func_flip_test.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

greens_func_wrap_test.o: ../test/greens_func_wrap_test.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

greens_func_init_test1.o: ../test/greens_func_init_test1.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

greens_func_init_test2.o: ../test/greens_func_init_test2.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

greens_func_init_test3.o: ../test/greens_func_init_test3.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

greens_func_init_test4.o: ../test/greens_func_init_test4.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

monte_carlo_iter_test.o: ../test/monte_carlo_iter_test.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

monte_carlo_phonon_block_test.o: ../test/monte_carlo_phonon_block_test.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

monte_carlo_iter_phonon_test.o: ../test/monte_carlo_iter_phonon_test.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

green_unequal_time_test.o: ../test/green_unequal_time_test.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<

measurement_test.o: ../test/measurement_test.c
	${CC} -c -I../include ${CCOPTS} ${MKL} $<


clean:
	rm *.o
