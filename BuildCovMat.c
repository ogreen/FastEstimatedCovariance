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



#include <string.h>
#include <math.h>
#include "Utilities.h"
#include "Config.h"
#include "Globals.h"
#include "BuildCovMat.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

//#include "Timing.h"

#include "Utilities.h"

#include "timing_util.h"
#include <omp.h>
void SortCombinations(int numberOfPartitions);

double BuildCovMatNaive(void)
{
	int idx_row, idx_col; // top-left coordinates of current window
	int i, j; // coordinates inside the window
	int m;

	double t;

	omp_set_num_threads(THREADS);
	ZeroOutputMatrix(outputMat,WINDOW_ROWS,WINDOW_COLS);

	BuildCombinationsTable(WINDOW_ROWS,WINDOW_COLS);
	ComplexFloat32*** outputForAllThreads[THREADS];
	ComplexFloat32** windowVecForAllThreads[THREADS];
	for (int t=0; t<THREADS; t++)
	{
		outputForAllThreads[t]=createComplexArray(WINDOW_ROWS*WINDOW_COLS,WINDOW_ROWS*WINDOW_COLS);

		windowVecForAllThreads[t]=malloc(sizeof(ComplexFloat32)*WINDOW_ROWS*WINDOW_COLS);  

	}

	tic_reset();


	// printf("1");

#pragma omp parallel for 
	for (idx_row = 0; idx_row < INPUT_ROWS-WINDOW_ROWS+1; ++idx_row)
	{
		int thread=omp_get_thread_num();
		ComplexFloat32** outputMatT=outputForAllThreads[thread];

		ComplexFloat32* windowVec=windowVecForAllThreads[thread];  

		for (idx_col = 0; idx_col < INPUT_COLS-WINDOW_COLS+1; ++idx_col)
		{
			// this loops is where each window get processed

			// copy all elements in window to the vector, note that we first go over columns!
			m = 0;
			for (i=0; i<WINDOW_COLS; ++i)
			{
				for (j=0; j<WINDOW_ROWS; ++j)
				{
					windowVec[m].real =inputMat[j+idx_col][i+idx_row].real;
					windowVec[m].imag =inputMat[j+idx_col][i+idx_row].imag;
					m++;
				}
			}
			// do the calculation
			for (i=0; i<WINDOW_ROWS*WINDOW_COLS; ++i)
			{
				for (j=0; j<WINDOW_ROWS*WINDOW_COLS; ++j)
				{
					if (idx_row==0 && idx_col==0) //using short-circuit to minimize checks
					{
						// this is the first run, so assume the output is full of garbage
						MULT_CONJ(windowVec[i], windowVec[j], outputMatT[i][j]);
					}
					else
					{
						ADD_MULT_CONJ(windowVec[i], windowVec[j], outputMatT[i][j]);
					}
				}
			}
		}
	}
	for (int thread=0; thread<THREADS; thread++   )
	{
		ComplexFloat32** threadMat=outputForAllThreads[thread];

		for (idx_row = 0; idx_row < INPUT_ROWS-WINDOW_ROWS+1; ++idx_row)
		{ 
			for (idx_col = 0; idx_col < INPUT_COLS-WINDOW_COLS+1; ++idx_col)
			{
				outputMat[idx_row][idx_col].real += threadMat[idx_row][idx_col].real;
				outputMat[idx_row][idx_col].imag += threadMat[idx_row][idx_col].imag;
			}
		}      
	}

	t=tic_sinceReset();
	for (int t=0; t<THREADS; t++)
	{ 
		free(windowVecForAllThreads[t]);
		destroyComplexArray(outputForAllThreads[t],WINDOW_ROWS*WINDOW_COLS);


	}

	//	PRINT_BANNER
	return t;
}

double BuildCovMatParallel(void)
{
	double t;
	unsigned int idx;
	ZeroOutputMatrix(outputMat,WINDOW_ROWS,WINDOW_COLS);

	BuildCombinationsTable(WINDOW_ROWS,WINDOW_COLS);
	tic_reset();

	for (idx=0; idx<COMBINATIONS_COUNT; ++idx)
	{
		//		printf("%d \n",idx);
		RunCombination(idx,inputMat,outputMat,WINDOW_ROWS,WINDOW_COLS);
	}
	t=tic_sinceReset();

	PRINT_BANNER
		return t;

}

double BuildCovMatDiagonalWrites(void)
{
	double t;
	unsigned int idx;
	ZeroOutputMatrix(outputMat,WINDOW_ROWS,WINDOW_COLS);

	tic_reset();

	BuildCombinationsTable(WINDOW_ROWS,WINDOW_COLS);

	// i'm zero-ing the whole matrix just for ouputMat to look nice when
	// debugging, you can safely change it back when done.
	//ZeroOutputMatrix();

	for (idx=0; idx<COMBINATIONS_COUNT; ++idx)
	{
		RunCombinationDiagonalWrites(idx,inputMat,outputMat,WINDOW_ROWS,WINDOW_COLS);
	}

	t=tic_sinceReset();

	PRINT_BANNER
		return t;
}

double BuildCovMatDiagonalWritesEnhanced()
{
	double t;
	unsigned int idx;
	ZeroOutputMatrix(outputMat,WINDOW_ROWS,WINDOW_COLS);

	tic_reset();

	BuildCombinationsTable(WINDOW_ROWS,WINDOW_COLS);

	// i'm zero-ing the whole matrix just for ouputMat to look nice when
	// debugging, you can safely change it back when done.
	//ZeroOutputMatrix();

	ComplexFloat32 outputVector[WINDOW_ROWS*WINDOW_COLS];
	for (idx=0; idx<COMBINATIONS_COUNT; ++idx)
	{
		RunCombinationDiagonalWritesEnhanced(idx,inputMat,outputMat,WINDOW_ROWS,WINDOW_COLS,outputVector);
	}

	t=tic_sinceReset();

	PRINT_BANNER
		return t;

}

double BuildCovMatDiagonalJustMultiplicatedEnhancedParallel(int doSort, int doBuildTable)
{

	double t;
	// Get the number of processors in this system
	int iCPU;

	ZeroOutputMatrix(outputMato,WINDOW_ROWS,WINDOW_COLS);
	if(doBuildTable)
		BuildCombinationsTable(WINDOW_ROWS,WINDOW_COLS);

	iCPU = THREADS;

	if(doSort)
		SortCombinations(iCPU);

	// Now set the number of threads
	omp_set_num_threads(iCPU);


	/*   
	     unsigned int idx;
	     tic_reset();
	     BuildCombinationsTable(WINDOW_ROWS,WINDOW_COLS);

	// i'm zero-ing the whole matrix just for ouputMat to look nice when
	// debugging, you can safely change it back when done.
	//ZeroOutputMatrix();

	ComplexFloat32 outputVector[WINDOW_ROWS*WINDOW_COLS];
	for (idx=0; idx<COMBINATIONS_COUNT; ++idx)
	{
	//		printf("%d \n",idx);
	RunCombinationDiagonalWritesEnhanced(idx,inputMat,outputMat,WINDOW_ROWS,WINDOW_COLS,outputVector);
	}

	t=tic_sinceReset();

	 */   

	ComplexFloat32 outputVector[WINDOW_ROWS*WINDOW_COLS*iCPU];

	tic_reset();

	BuildCombinationsTable(WINDOW_ROWS,WINDOW_COLS);
	/*
	   unsigned int idx;
	//      for (idx=0; idx<COMBINATIONS_COUNT; idx++)
	for (idx=0; idx<COMBINATIONS_COUNT; ++idx)
	{
	//RunCombinationDiagonalOded(startComb++,inputMat,WINDOW_ROWS, WINDOW_COLS,
	//                           outputVector,finalMatPosR, finalMatPosC);
	//            if(thread_id==0)
	//                printf("%d,",idx) ;
	//     RunCombinationDiagonalWritesEnhanced(idx,inputMat,NULL,WINDOW_ROWS,WINDOW_COLS,outputVector);
	RunCombinationDiagonalWritesEnhanced(idx,inputMat,outputMat,WINDOW_ROWS,WINDOW_COLS,outputVector);
	}
	 */
	// i'm zero-ing the whole matrix just for ouputMat to look nice when
	// debugging, you can safely change it back when done.
	//ZeroOutputMatrix();

	int combsPerCPU = COMBINATIONS_COUNT/iCPU;
	int remainder = COMBINATIONS_COUNT%iCPU;
	if(remainder>0)
		combsPerCPU++;

	//   for(int s=0; s<1000;s++)
	//	printf("\n\n%d %d\n",COMBINATIONS_COUNT,combsPerCPU);
#pragma omp parallel
	{ 
		int idx;
		int combs;
		int thread_id= omp_get_thread_num();
		int startComb,stopComb;
		if (remainder > thread_id)
		{
			combs=combsPerCPU;
			startComb=(thread_id)*combsPerCPU;
		}
		else
		{
			combs=combsPerCPU-1;
			startComb=remainder*combsPerCPU+(thread_id-remainder)*(combsPerCPU-1);
		}
		stopComb=startComb+combs; 
		//        int my_counter = thread_id*combs;
		// printf("%d %d \n",startComb,combs	);

		//	   printf(" %d %d %d\n",thread_id,startComb,combs);
		int ptr = WINDOW_ROWS*WINDOW_ROWS*thread_id;
		//	int finalMatPosR[WINDOW_ROWS*WINDOW_COLS];
		//	int finalMatPosC[WINDOW_ROWS*WINDOW_COLS];
		for (idx=startComb; idx<stopComb; idx++)
		{
			//RunCombinationDiagonalOded(startComb++,inputMat,WINDOW_ROWS, WINDOW_COLS,
			//                           outputVector,finalMatPosR, finalMatPosC);
			//            if(thread_id==0)
			//                printf("%d,",idx) ;
			RunCombinationDiagonalWritesEnhanced(idx,inputMat,NULL,WINDOW_ROWS,WINDOW_COLS,outputVector+ptr);
		}
		//			printf("%d \n",idx);

		//         if(thread_id==0)
		//            printf("\n") ;
	}




	t=tic_sinceReset();

	//	t=tic_sincelast();
	//  t=tic_total();

	PRINT_BANNER
		return t;


}



double BuildCovMatDiagonalAddMul(void)
{
	double t;
	unsigned int idx;
	ZeroOutputMatrix(outputMato,WINDOW_ROWS,WINDOW_COLS);
	BuildCombinationsTable(WINDOW_ROWS,WINDOW_COLS);


	tic_reset();

	// i'm zero-ing the whole matrix just for ouputMat to look nice when
	// debugging, you can safely change it back when done.
	//ZeroOutputMatrix();

	ComplexFloat32 outputVector[WINDOW_ROWS*WINDOW_COLS];
	int finalMatPosR[WINDOW_ROWS*WINDOW_COLS];
	int finalMatPosC[WINDOW_ROWS*WINDOW_COLS];

	for (idx=0; idx<COMBINATIONS_COUNT; ++idx)
	{
		RunCombinationDiagonalAddMul(idx,inputMat,WINDOW_ROWS,WINDOW_COLS,outputVector,finalMatPosR, finalMatPosC);

	}
	t=tic_total();
	PRINT_BANNER
		return t;


}


double BuildCovMatDiagonalAddMulOPENMP(int doSort, int doBuildTable)
{

	double t;
	// Get the number of processors in this system
	int iCPU;

	ZeroOutputMatrix(outputMato,WINDOW_ROWS,WINDOW_COLS);
	if(doBuildTable)
		BuildCombinationsTable(WINDOW_ROWS,WINDOW_COLS);

	iCPU = THREADS;

	if(doSort)
		SortCombinations(iCPU);

	// Now set the number of threads
	omp_set_num_threads(iCPU);

	tic_reset();

	// i'm zero-ing the whole matrix just for ouputMat to look nice when
	// debugging, you can safely change it back when done.
	//ZeroOutputMatrix();

	int combsPerCPU = COMBINATIONS_COUNT/iCPU;
	int remainder = COMBINATIONS_COUNT%iCPU;
	if(remainder>0)
		combsPerCPU++;

#pragma omp parallel
	{
		int idx;
		int combs;
		int thread_id= omp_get_thread_num();
		int startComb,stopComb;
		if (remainder > thread_id)
		{
			combs=combsPerCPU;
			startComb=(thread_id)*combsPerCPU;
		}
		else
		{
			combs=combsPerCPU-1;
			startComb=remainder*combsPerCPU+(thread_id-remainder)*(combsPerCPU-1);
		}
		stopComb=startComb+combs; 
		//        int my_counter = thread_id*combs;
		// printf("%d %d \n",startComb,combs	);

		ComplexFloat32 outputVector[WINDOW_ROWS*WINDOW_COLS];
		int finalMatPosR[WINDOW_ROWS*WINDOW_COLS];
		int finalMatPosC[WINDOW_ROWS*WINDOW_COLS];
		for (idx=0; idx<combs; idx++)
		{
			RunCombinationDiagonalAddMul(startComb++,inputMat,WINDOW_ROWS, WINDOW_COLS,
					outputVector,finalMatPosR, finalMatPosC);

		}


	}

	t=tic_total();

	PRINT_BANNER
		return t;

}

void SortCombinations(int numberOfPartitions)
{
	Combination tempCombinationsTable[COMBINATIONS_COUNT];
	Combination sortedCombinationsTable[COMBINATIONS_COUNT];
	int c=0;

	for(c=0;c<COMBINATIONS_COUNT; c++)
	{
		tempCombinationsTable[c].blockSize	= combinationsTable[c].blockSize;
		tempCombinationsTable[c].mult1r		= combinationsTable[c].mult1r;
		tempCombinationsTable[c].mult1c		= combinationsTable[c].mult1c;
		tempCombinationsTable[c].mult2r		= combinationsTable[c].mult2r;
		tempCombinationsTable[c].mult2c		= combinationsTable[c].mult2c;
		tempCombinationsTable[c].refCoords	= combinationsTable[c].refCoords;
		tempCombinationsTable[c].type2		= combinationsTable[c].type2;
		tempCombinationsTable[c].id			= combinationsTable[c].id;

	}

	for(c=0;c<COMBINATIONS_COUNT; c++)
	{
		int max=0,max_c=0;
		int j;
		for(j=0; j<COMBINATIONS_COUNT;j++)
		{
			Combination* comb = &tempCombinationsTable[j];

			unsigned int mult1r = comb->mult1r ;
			unsigned int mult1c = comb->mult1c;
			unsigned int mult2r = comb->mult2r;
			unsigned int mult2c = comb->mult2c;
			unsigned int deltaR = abs(mult1r-mult2r);
			unsigned int deltaC = abs(mult1c-mult2c);
			int val;
			if(comb->mult1r==INPUT_ROWS && comb->mult2c == INPUT_COLS)
				continue;

			val = (INPUT_ROWS-deltaR)*(INPUT_COLS-deltaC);

			if(val>max)
			{
				max_c=j;
				max = val;
			}
		}
		sortedCombinationsTable[c].blockSize	= tempCombinationsTable[max_c].blockSize;
		sortedCombinationsTable[c].mult1r		= tempCombinationsTable[max_c].mult1r;
		sortedCombinationsTable[c].mult1c		= tempCombinationsTable[max_c].mult1c;
		sortedCombinationsTable[c].mult2r		= tempCombinationsTable[max_c].mult2r;
		sortedCombinationsTable[c].mult2c		= tempCombinationsTable[max_c].mult2c;
		sortedCombinationsTable[c].refCoords	= tempCombinationsTable[max_c].refCoords;
		sortedCombinationsTable[c].type2		= tempCombinationsTable[max_c].type2;
		sortedCombinationsTable[c].id			= tempCombinationsTable[max_c].id;

		tempCombinationsTable[max_c].mult1r=INPUT_ROWS;
		tempCombinationsTable[max_c].mult2c=INPUT_COLS;
		tempCombinationsTable[max_c].blockSize=0;
		tempCombinationsTable[max_c].refCoords=0;
	}

	{
		int parti;

		int f_cores = COMBINATIONS_COUNT-(COMBINATIONS_COUNT/numberOfPartitions)*numberOfPartitions;

		int w_first_partitions = COMBINATIONS_COUNT/numberOfPartitions;
		int w_last_partitions = COMBINATIONS_COUNT/numberOfPartitions;

		int wCounter = 0;

		if(f_cores>0)
		{
			w_first_partitions++;
		}

		for(parti=0; parti<numberOfPartitions;parti++)
		{
			int work,w;
			if(parti<f_cores)
			{
				work = w_first_partitions;
			}
			else
			{
				work = w_last_partitions;
			}
			for(w=0; w<work;w++)
			{
				int id=numberOfPartitions*w+parti;
				combinationsTable[wCounter].blockSize	= sortedCombinationsTable[id].blockSize;
				combinationsTable[wCounter].mult1r		= sortedCombinationsTable[id].mult1r;
				combinationsTable[wCounter].mult1c		= sortedCombinationsTable[id].mult1c;
				combinationsTable[wCounter].mult2r		= sortedCombinationsTable[id].mult2r;
				combinationsTable[wCounter].mult2c		= sortedCombinationsTable[id].mult2c;
				combinationsTable[wCounter].refCoords	= sortedCombinationsTable[id].refCoords;
				combinationsTable[wCounter].type2		= sortedCombinationsTable[id].type2;
				combinationsTable[wCounter].id			= sortedCombinationsTable[id].id;
				//				printf("%d ",id);
				wCounter++;
			}
			//			printf("\n");
		}
	}
}


void BuildCovMatZeroPadding(void)
{

	PRINT_BANNER
}

void BuildCovMatOnTheFlyIndicsCalculations(void)
{

	PRINT_BANNER
}
