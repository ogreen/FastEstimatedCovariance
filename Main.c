/*
 * Authors: 
 *  Oded Green (ogreen@gatech.edu)
 *
 * Publications (please cite):
 * 
 * Optimizes both the number of multiplications and the number of additions:
 * O. Green, Y. Birk, "A Computationally Efficient Algorithm for the 2D Covariance Method", ACM/IEEE International Conference on High Performance Computing, Networking, Storage and Analysis, Denver, Colorado, 2013
 * 
 * Reduces and optimizes only the number of multiplications:
 * O. Green, L. David, A. Galperin, Y. Birk, "Efficient parallel computation of the estimated covariance matrix", arXiv, 2013
 * L. David, A. Galperin, O. Green, Y. Birk, "Efficient parallel computation of the estimated covariance matrix" , 26th IEEE Convention of Electrical and Electronics Engineers in Israel (IEEEI), 2010
 * 
 *
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 * 
 * - Redistributions of source code must retain the above copyright notice, 
 *   this list of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, 
 *   this list of conditions and the following disclaimer in the documentation 
 *   and/or other materials provided with the distribution.
 * - Neither the name of the Georgia Institute of Technology nor the names of 
 *   its contributors may be used to endorse or promote products derived from 
 *   this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF 
 * THE POSSIBILITY OF SUCH DAMAGE.
 */



#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <math.h>
#include <errno.h>

#include "BuildCovMat.h"
#include "Utilities.h"
#include <stdio.h>
#include "Config.h"

#include "Globals.h"

typedef enum BenchMarkTypeOn
{
	BM_STRAIGHT_FORWARD=1,
	BM_PARALLEL_SIMPLE=2,
	BM_PARALLEL_DIAGS_AS_ARRAY=4,
	BM_PARALLEL_DIAGS_AS_ARRAY_OPTIMIZED=8,
	BM_PARALLEL_DIAGS_AS_ARRAY_OPTIMIZED_PAR=16,
	BM_SEQUENTIAL_MUL_ADD=32,
	BM_PARALLEL_MUL_ADD=64,
} BenchMarkBinaryArray;



#define RUN_FUNC(benchmark_flag,i,iterations, timer, func) \
{\
	if(benchmark_flag){\
		for (i=1; i<=iterations; i++){\
			tmp = func();\
			if (i==1){\
				timer = tmp;\
			} else {\
				timer+=tmp;\
			}\
		}\
		printf("%lf,",  timer/(double)iterations);\
	}\
}\


#define RUN_FUNC_2PARS(benchmark_flag,i,iterations, timer, func,par1, par2) \
{\
	if(benchmark_flag){\
		for (i=1; i<=iterations; i++){\
			tmp = func(par1,par2);\
			if (i==1){\
				timer = tmp;\
			} else {\
				timer+=tmp;\
			}\
		}\
		printf("%lf,",  timer/(double)iterations);\
	}\
}\


void hostParseArgsDataStructureTesting(int argc, char** argv);
void CreateDummyInput(ComplexFloat32** inputMat, int INPUT_ROWS,int INPUT_COLS);

ComplexFloat32** createComplexArray(int rows,int cols);
void destroyComplexArray(ComplexFloat32** array,int rows);


// using function pointer
double RunAndTime(double (*)(void));

double RunAndTime(double (*funcPointer)(void))
{
	double timingResults = (*funcPointer)();
	return timingResults;
}


int main(int argc, char* argv[])
{
	THREADS=1;
	ITERATIONS=1;
	BENCHMARK_FLAGS=0;
	SORT_FLAG=0;
	hostParseArgsDataStructureTesting(argc,argv);
	if(INPUT_ROWS==0)
		exit(1);

	NR=INPUT_ROWS;
	NC=INPUT_COLS;
	PR=WINDOW_ROWS;
	PC=WINDOW_COLS;

	int combCount = COMBINATIONS_COUNT;
	int sortFlag = SORT_FLAG;

	windowVec=NULL;
	combinationsTable = (Combination*)malloc(sizeof(Combination)*combCount);
	windowVec = (ComplexFloat32*)malloc(sizeof(ComplexFloat32)*WINDOW_ROWS*WINDOW_COLS);

	inputMat=createComplexArray(INPUT_ROWS,INPUT_COLS);
	outputMat=NULL;
	outputMat=createComplexArray(WINDOW_ROWS*WINDOW_COLS,WINDOW_ROWS*WINDOW_COLS);
	outputMato=createComplexArray(WINDOW_ROWS*WINDOW_COLS,WINDOW_ROWS*WINDOW_COLS);



	// the window vector which is created by column-stacking each window in
	// the serial algorithm.
	// it is also used for the VectorWrites implementation


	int i=0;
	double naive=0, parallel=0, diag=0, diag_enhanced=0, tmp=0, diag_enhanced_par=0;
	double diag_add_mul=0,diag_add_mul_omp=0;


	CreateDummyInput(inputMat,INPUT_ROWS,INPUT_COLS);


	RUN_FUNC(BENCHMARK_FLAGS&BM_STRAIGHT_FORWARD,i,ITERATIONS, naive, BuildCovMatNaive);
	RUN_FUNC(BENCHMARK_FLAGS&BM_PARALLEL_SIMPLE,i,ITERATIONS, parallel, BuildCovMatParallel);

	RUN_FUNC(BENCHMARK_FLAGS&BM_PARALLEL_DIAGS_AS_ARRAY,i,ITERATIONS, diag, BuildCovMatDiagonalWrites);
	RUN_FUNC(BENCHMARK_FLAGS&BM_PARALLEL_DIAGS_AS_ARRAY_OPTIMIZED,i,ITERATIONS, diag_enhanced, BuildCovMatDiagonalWritesEnhanced);


	RUN_FUNC_2PARS(BENCHMARK_FLAGS&BM_PARALLEL_DIAGS_AS_ARRAY_OPTIMIZED_PAR,i,ITERATIONS, diag_enhanced_par, BuildCovMatDiagonalJustMultiplicatedEnhancedParallel,SORT_FLAG,0);

	RUN_FUNC(BENCHMARK_FLAGS&BM_SEQUENTIAL_MUL_ADD,i,ITERATIONS, diag_add_mul, BuildCovMatDiagonalAddMul);


	RUN_FUNC_2PARS(BENCHMARK_FLAGS&BM_PARALLEL_MUL_ADD,i,ITERATIONS, diag_add_mul_omp, BuildCovMatDiagonalAddMulOPENMP,SORT_FLAG,0);


	/*
	   Compare2Mats();
	   naive /= (double)ITERATIONS;
	   parallel /= (double)ITERATIONS;
	   diag /= (double)ITERATIONS;
	   diag_enhanced /= (double)ITERATIONS;
	   diag_add_mul /= (double)ITERATIONS;
	   diag_add_mul_omp /=(double)ITERATIONS;

	// the next two line are very important. don't remove or change their structure
	printf("** ALGORITHMS,Naive,Parallel,DiagonalWrites,DiagonalWritesEnhanced,DiagonalOded,DiagonalOdedOMP\n");
	printf("** RESULTS,%f,%f,%f,%f,%f,%f\n\n", naive, parallel, diag, diag_enhanced,diag_add_mul,diag_add_mul_omp);

	printf("Speedup is by a factor of sequential = %f , omp = %f\n", naive/diag_add_mul,naive/diag_add_mul_omp);
	printf("omp/sequential = %f\n", diag_add_mul/diag_add_mul_omp);
	 */
	//printf("*** AVERAGE, Naive/Parallel Ratio = %f\n", naive/parallel);



	free(combinationsTable);
	free(windowVec);

	destroyComplexArray(inputMat,INPUT_ROWS);
	destroyComplexArray(outputMat,WINDOW_ROWS*WINDOW_COLS);
	destroyComplexArray(outputMato,WINDOW_ROWS*WINDOW_COLS);





	return 0;
}

void hostParseArgsDataStructureTesting(int argc, char** argv) {
	static struct option long_options[] = {
		{"Input", required_argument, 0, 'I'},
		{"Window", required_argument, 0, 'W'},
		{"Runs", required_argument, 0, 'R'},
		{"ThreadCount",optional_argument,0,'T'},
		{"Benchmark",optional_argument,0,'B'},
		{"SortFlag",optional_argument,0,'S'},
		{0, 0, 0, 0}
	};


	while(1) {
		int option_index = 0;
		int c = getopt_long(argc, argv, "I:W:R:T:B:S:?h", long_options, &option_index);
		extern char * optarg;
		extern int    optind, opterr, optopt;
		int intout = 0;

		if(-1 == c)
			break;

		switch(c) {
			default:
				printf("Unrecognized option: %c\n\n", c);
			case '?':
			case 'h':
				printf("\nCreates two arrays A and B and merges them into array C in parallel on GPU.");
				printf("\n\nUsage"
						"\n====="
						"\n\n\t-I --Input <number>\n\t\tSpecify the input array size , I x I .\n"
						"\n\t-W --Window <number>\n\t\tSpecify the window size, W x W.\n"
						"\n\t-R --Runs \n\t\tNumber of iterations to complete\n"
						"\n\t-T --Threads \n\n\tNumber of threads\n"
						"\n\t-B --Benchmarks \n\t\tSelection of the benchmarks that need to be run.\n"
						"\n\t-S --Sort \n\t\tSort flag - relevant only for part of the functions.\n");

				exit(0);
				break;

			case 'I':
				errno = 0;
				intout = strtol(optarg, NULL, 10);
				if(errno || intout < 0) {
					printf("Error - Input size must be positive %s\n", optarg);
					exit(-1);
				}
				INPUT_ROWS=INPUT_COLS=intout;
				break;
			case 'W':
				errno = 0;
				intout = strtol(optarg, NULL, 10);
				if(errno || intout < 0) {
					printf("Error - Window size must be positive %s\n", optarg);
					exit(-1);
				}
				WINDOW_ROWS=WINDOW_COLS=intout;
				break;
			case 'R':
				errno = 0;
				intout = strtol(optarg, NULL, 10);
				if(errno || intout < 0) {
					printf("Number of runs must be positive%s\n", optarg);
					exit(-1);
				}
				ITERATIONS=intout;
				break;
			case 'T':
				errno = 0;
				intout = strtol(optarg, NULL, 10);
				if(errno || intout < 0) {
					printf("Number of threads must be positive %s\n", optarg);
					exit(-1);
				}
				THREADS = intout;
				break;
			case 'B':
				errno = 0;
				intout = strtol(optarg, NULL, 10);
				if(errno || intout < 0) {
					printf("The benchmark numbers should be positive %s\n", optarg);
					exit(-1);
				}
				BENCHMARK_FLAGS = intout;
				break;
			case 'S':
				errno = 0;
				intout = strtol(optarg, NULL, 10);
				if(errno || intout < 0) {
					printf("The benchmark numbers should be positive %s\n", optarg);
					exit(-1);
				}
				SORT_FLAG = intout;
				break;
		}
	}
}



void CreateDummyInput(ComplexFloat32** inputMat, int INPUT_ROWS,int INPUT_COLS)
{
	// use dummy data
	int r,c;
	for(r=0;r<INPUT_ROWS;r++)
	{
		for(c=0;c<INPUT_COLS;c++)
		{
			inputMat[r][c].imag=(float)(r*INPUT_ROWS+c);
			inputMat[r][c].real=(float)(c*INPUT_COLS+r);
		}
	}
}

ComplexFloat32** createComplexArray(int rows,int cols)
{
	ComplexFloat32** array = (ComplexFloat32**)malloc(rows*sizeof(ComplexFloat32*));

	for(int r=0; r<rows;r++)
	{
		array[r]=(ComplexFloat32*)malloc(cols*sizeof(ComplexFloat32));
	}

	return array;
}

void destroyComplexArray(ComplexFloat32** array,int rows)
{
	for(int r=0; r<rows;r++)
	{
		free(array[r]);
	}

	free(array);
}


void ZeroOutputMatrix(ComplexFloat32** matToZero, int rows, int cols)
{
	for (int i=0; i<rows; ++i)
	{
		for (int j=0; j<cols; ++j)
		{
			matToZero[i][j].real = 0.0;
			matToZero[i][j].imag = 0.0;
		}
	}
}

