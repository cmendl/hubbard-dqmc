
SRCFILES = dupio.o greens_func.o kinetic.o lattice_fourier.o linalg.o measurement.o monte_carlo.o param_parser.o random.o stratonovich.o time_flow.o util.o
TSTFILES = matrix_det_test.o matrix_exp_test.o blockcyclic_qr_test.o blockcyclic_tri_test.o blockcyclic_inv_test.o kinetic_test.o lattice_fourier_test.o time_flow_test1.o time_flow_test2.o time_flow_test3.o stratonovich_test.o greens_func_flip_test.o greens_func_wrap_test.o greens_func_init_test1.o greens_func_init_test2.o greens_func_init_test3.o monte_carlo_iter_test.o monte_carlo_phonon_block_test.o monte_carlo_iter_phonon_test.o

# use Intel compiler
CC = icc

# set these with appropriate libraries for your system
LIBRARIES = -L/opt/intel/mkl/lib/intel64/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread

# compiler options
CCOPTS = -Wall -O2 -xHost -DNDEBUG -restrict -DMEM_DATA_ALIGN=64


all: proj_main proj_test

proj_main: ${SRCFILES} main.o
	${CC} ${CCOPTS} ${LIBRARIES} -o hubbard_dqmc $?

proj_test: ${SRCFILES} ${TSTFILES} run_tests.o
	${CC} ${CCOPTS} ${LIBRARIES} -o run_tests $?


%.o: ../src/%.c
	${CC} -c -I../include ${CCOPTS} $<


run_tests.o: ../test/run_tests.c
	${CC} -c -I../include ${CCOPTS} $<


matrix_det_test.o: ../test/matrix_det_test.c
	${CC} -c -I../include ${CCOPTS} $<

matrix_exp_test.o: ../test/matrix_exp_test.c
	${CC} -c -I../include ${CCOPTS} $<

blockcyclic_qr_test.o: ../test/blockcyclic_qr_test.c
	${CC} -c -I../include ${CCOPTS} $<

blockcyclic_tri_test.o: ../test/blockcyclic_tri_test.c
	${CC} -c -I../include ${CCOPTS} $<

blockcyclic_inv_test.o: ../test/blockcyclic_inv_test.c
	${CC} -c -I../include ${CCOPTS} $<

kinetic_test.o: ../test/kinetic_test.c
	${CC} -c -I../include ${CCOPTS} $<

lattice_fourier_test.o: ../test/lattice_fourier_test.c
	${CC} -c -I../include ${CCOPTS} $<

time_flow_test1.o: ../test/time_flow_test1.c
	${CC} -c -I../include ${CCOPTS} $<

time_flow_test2.o: ../test/time_flow_test2.c
	${CC} -c -I../include ${CCOPTS} $<

time_flow_test3.o: ../test/time_flow_test3.c
	${CC} -c -I../include ${CCOPTS} $<

stratonovich_test.o: ../test/stratonovich_test.c
	${CC} -c -I../include ${CCOPTS} $<

greens_func_flip_test.o: ../test/greens_func_flip_test.c
	${CC} -c -I../include ${CCOPTS} $<

greens_func_wrap_test.o: ../test/greens_func_wrap_test.c
	${CC} -c -I../include ${CCOPTS} $<

greens_func_init_test1.o: ../test/greens_func_init_test1.c
	${CC} -c -I../include ${CCOPTS} $<

greens_func_init_test2.o: ../test/greens_func_init_test2.c
	${CC} -c -I../include ${CCOPTS} $<

greens_func_init_test3.o: ../test/greens_func_init_test3.c
	${CC} -c -I../include ${CCOPTS} $<

monte_carlo_iter_test.o: ../test/monte_carlo_iter_test.c
	${CC} -c -I../include ${CCOPTS} $<

monte_carlo_phonon_block_test.o: ../test/monte_carlo_phonon_block_test.c
	${CC} -c -I../include ${CCOPTS} $<

monte_carlo_iter_phonon_test.o: ../test/monte_carlo_iter_phonon_test.c
	${CC} -c -I../include ${CCOPTS} $<


clean:
	rm *.o
