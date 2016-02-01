#include <mkl.h>
#include <stdbool.h>
#include <stdio.h>

#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif


typedef int (*test_function_t)();

// test function declarations
int MatrixExpTest();
int BlockCyclicQRTest();
int BlockCyclicTriTest();
int BlockCyclicInvTest();
int KineticTest();
int KineticTest2();
int KineticTest3();
int LatticeFourierTest();
int TimeFlowTest1();
int TimeFlowTest2();
int TimeFlowTest3();
int StratonovichTest();
int GreensFuncFlipTest();
int GreensFuncWrapTest();
int GreensFuncInitTest1();
int GreensFuncInitTest2();
int GreensFuncInitTest3();
int GreensFuncInitTest4();
int MonteCarloIterTest();
int MonteCarloPhononBlockTest();
int MonteCarloIterPhononTest();
int GreenUnequalTimeTest();
int MeasurementTest();


int main()
{
	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	test_function_t tests[] = { MatrixExpTest, BlockCyclicQRTest, BlockCyclicTriTest, BlockCyclicInvTest, KineticTest, KineticTest2, KineticTest3, LatticeFourierTest, TimeFlowTest1, TimeFlowTest2, TimeFlowTest3, StratonovichTest, GreensFuncFlipTest, GreensFuncWrapTest, GreensFuncInitTest1, GreensFuncInitTest2, GreensFuncInitTest3, GreensFuncInitTest4, MonteCarloIterTest, MonteCarloPhononBlockTest, MonteCarloIterPhononTest, GreenUnequalTimeTest, MeasurementTest };

	bool pass = true;

	int i;
	for (i = 0; i < (int)(sizeof(tests)/sizeof(tests[0])); i++)
	{
		int status = tests[i]();

		// print status message
		if (status < 0)
		{
			printf("Test crashed with internal problem!\n");
		}
		else if (status > 0)
		{
			printf("Test failed!\n");
		}
		else	// status == 0
		{
			printf("Test passed.\n");
		}

		printf("________________________________________________________________________________\n");

		if (status != 0) {
			pass = false;
		}
	}

	if (pass)
	{
		printf("All tests passed :-)\n");
	}
	else
	{
		printf("At least one test failed or crashed :-(\n");
	}

	// clean up MKL library
	MKL_Free_Buffers();

	// check for MKL memory leaks
	int nbuffers;
	MKL_INT64 nbytes_alloc;
	nbytes_alloc = MKL_Mem_Stat(&nbuffers);
	if (nbytes_alloc > 0)
	{
		printf("\nMKL memory leak detected! MKL still uses %lld bytes in %d buffer(s).\n", nbytes_alloc, nbuffers);
	}
	else
	{
		printf("\nMKL memory leak check appears to be fine.\n");
	}

	return 0;
}
