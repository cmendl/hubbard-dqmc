#include <stdbool.h>
#include <stdio.h>

#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
#include <crtdbg.h>
#endif

void mkl_free_buffers(void);


typedef int (*test_function_t)();

// test function declarations
int MatrixExpTest();
int BlockCyclicQRTest();
int BlockCyclicTriTest();
int BlockCyclicInvTest();
int KineticTest();
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


int main()
{
	// enable run-time memory check for debug builds
	#if defined(_WIN32) & (defined(DEBUG) | defined(_DEBUG))
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	#endif

	#define NUM_TESTS 20
	test_function_t tests[NUM_TESTS] = { MatrixExpTest, BlockCyclicQRTest, BlockCyclicTriTest, BlockCyclicInvTest, KineticTest, LatticeFourierTest, TimeFlowTest1, TimeFlowTest2, TimeFlowTest3, StratonovichTest, GreensFuncFlipTest, GreensFuncWrapTest, GreensFuncInitTest1, GreensFuncInitTest2, GreensFuncInitTest3, GreensFuncInitTest4, MonteCarloIterTest, MonteCarloPhononBlockTest, MonteCarloIterPhononTest, GreenUnequalTimeTest };

	bool pass = true;

	int i;
	for (i = 0; i < NUM_TESTS; i++)
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
	mkl_free_buffers();

	return 0;
}
