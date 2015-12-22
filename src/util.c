#include "util.h"
#include "dupio.h"
#include <time.h>
#ifdef _WIN32
#include <windows.h>
#endif


//________________________________________________________________________________________________________________________
///
/// \brief Read 'n' items of size 'size' from file 'filename', expecting the file size to be exactly n*size
///
int ReadData(const char *filename, void *data, const size_t size, const size_t n)
{
	FILE *fd = fopen(filename, "rb");
	if (fd == NULL)
	{
		duprintf("'fopen()' failed during call of 'ReadData()'.\n");
		return -1;
	}

	// obtain the file size
	fseek(fd, 0, SEEK_END);
	long filesize = ftell(fd);
	rewind(fd);
	// duprintf("file size: %d\n", filesize);
	if ((size_t)filesize != n*size)
	{
		duprintf("'ReadData()' failed: expected file size does not match.\n");
		return -2;
	}

	// copy the file into the data array
	if (fread(data, size, n, fd) != n)
	{
		duprintf("'fread()' failed during call of 'ReadData()'.\n");
		return -3;
	}

	fclose(fd);

	return 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Write 'n' items of size 'size' to file 'filename'
///
int WriteData(const char *filename, const void *data, const size_t size, const size_t n, const bool append)
{
	const char *mode = append ? "ab" : "wb";

	FILE *fd = fopen(filename, mode);
	if (fd == NULL)
	{
		duprintf("'fopen()' failed during call of 'WriteData()'.\n");
		return -1;
	}

	// write data array to file
	if (fwrite(data, size, n, fd) != n)
	{
		duprintf("'fwrite()' failed during call of 'WriteData()'.\n");
		return -3;
	}

	fclose(fd);

	return 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Get current time tick
///
uint64_t GetTicks()
{
	#ifdef _WIN32

	LARGE_INTEGER t;
	QueryPerformanceCounter(&t);
	return (uint64_t)(t.QuadPart);

	#else

	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return (uint64_t)(1000000000ULL * t.tv_sec + t.tv_nsec);

	#endif
}
