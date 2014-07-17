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



#include <omp.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "Globals.h"
#include "Types.h"
#include "Utilities.h"

void multconj(const ComplexFloat32 *a, const ComplexFloat32 *b, ComplexFloat32 *dst)
{
	dst->real = (a->real * b->real) + (a->imag * b->imag);
	dst->imag = (b->real * a->imag) - (a->real * b->imag);
}


void RunCombination(unsigned int comb_id, ComplexFloat32** inputMat,ComplexFloat32** outputMat,int WINDOW_ROWS,int WINDOW_COLS)
{
	// extract data from the combinations table
	Combination* comb = &combinationsTable[comb_id];
	unsigned int mult1r = comb->mult1r;
	unsigned int mult1c = comb->mult1c;
	unsigned int mult2r = comb->mult2r;
	unsigned int mult2c = comb->mult2c;
	unsigned int blockR = comb->blockSize >> 8;
	unsigned int blockC = comb->blockSize & 0xFF;
	unsigned int refCoordR = comb->refCoords >> 16;
	unsigned int refCoordC = comb->refCoords & 0xFFFF;


	int refR, refC;

	unsigned int nw_r, nw_c;
	int shamt_r, shamt_c;
	ComplexFloat32 res; // the multiplication result
	int m1r, m1c, m2r, m2c; // coordinates of multipliers
	int min_rows, min_cols, max_rows, max_cols; // shifting boundries

	// iterate over all possible NW-positions for this combination
	for (nw_r=0; nw_r<INPUT_ROWS-blockR; ++nw_r)
	{
		for (nw_c=0; nw_c<INPUT_COLS-blockC; ++nw_c)
		{
			// calculate the initial multipliers' coordinates in the
			// input-image axis
			m1r = nw_r + mult1r;
			m1c = nw_c + mult1c;
			m2r = nw_r + mult2r;
			m2c = nw_c + mult2c;

			// make the multiplication
			multconj(&inputMat[m1r][m1c], &inputMat[m2r][m2c], &res);

			// calculate boundaries for shifts
			max_rows = MYMIN(WINDOW_ROWS - blockR - 1, nw_r);
			max_cols = MYMIN(WINDOW_COLS - blockC - 1, nw_c);
			min_rows = MYMAX(0, ((int)nw_r + WINDOW_ROWS - INPUT_ROWS));
			min_cols = MYMAX(0, ((int)nw_c + WINDOW_COLS - INPUT_COLS));

			// iterate over all possible shifts
			for (shamt_r=min_rows; shamt_r<=max_rows; ++shamt_r)
			{
				for (shamt_c=min_cols; shamt_c<=max_cols; ++shamt_c)
				{
					refR = refCoordR+shamt_r+shamt_c*WINDOW_ROWS;
					refC = refCoordC+shamt_r+shamt_c*WINDOW_ROWS;
					outputMat[refR][refC].real += res.real;
					outputMat[refR][refC].imag += res.imag;
				}
			}
		}
	}
	//	PRINT_BANNER
}

void RunCombinationDiagonalWrites(unsigned int comb_id, ComplexFloat32** inputMat,ComplexFloat32** outputMat,int WINDOW_ROWS,int WINDOW_COLS)
{
	// extract data from the combinations table, just for convenience
	Combination* comb = &combinationsTable[comb_id];
	unsigned int mult1r = comb->mult1r;
	unsigned int mult1c = comb->mult1c;
	unsigned int mult2r = comb->mult2r;
	unsigned int mult2c = comb->mult2c;
	unsigned int blockR = comb->blockSize >> 8;
	unsigned int blockC = comb->blockSize & 0xFF;

	unsigned int nw_r, nw_c;
	int shamt_r, shamt_c;
	ComplexFloat32 res; // the multiplication result
	int m1r, m1c, m2r, m2c; // coordinates of multipliers
	int min_rows, min_cols, max_rows, max_cols; // shifting boundries

	int oPos;

	ComplexFloat32* pLine;
	ComplexFloat32* pRes;
	int initialOffest = 0; // this will be changed only when comb_id>=N (a.k.a set2)

	// calculate the initial position in outputMat for writes
	if (comb_id<WINDOW_ROWS*WINDOW_COLS) {
		pLine = outputMat[comb_id];
	}
	else {
		// the number of the diagonal (its distance from the main diagonal)
		pLine = outputMat[(mult2c-1)*WINDOW_ROWS + WINDOW_ROWS-mult1r];
		initialOffest = mult1r;
	}

	// iterate over all possible NW-positions for this combination
	for (nw_r=0; nw_r<INPUT_ROWS-blockR; ++nw_r)
	{
		for (nw_c=0; nw_c<INPUT_COLS-blockC; ++nw_c)
		{
			// calculate the initial multipliers' coordinates in the
			// input-image axis
			m1r = nw_r + mult1r;
			m1c = nw_c + mult1c;
			m2r = nw_r + mult2r;
			m2c = nw_c + mult2c;

			// make the multiplication
			multconj(&inputMat[m1r][m1c], &inputMat[m2r][m2c], &res);

			// calculate boundaries for shifts
			max_rows = MYMIN(WINDOW_ROWS - blockR - 1, nw_r);
			max_cols = MYMIN(WINDOW_COLS - blockC - 1, nw_c);
			min_rows = MYMAX(0, (int)nw_r + WINDOW_ROWS - INPUT_ROWS);
			min_cols = MYMAX(0, (int)nw_c + WINDOW_COLS - INPUT_COLS);

			// iterate over all possible shifts
			for (shamt_r=min_rows; shamt_r<=max_rows; ++shamt_r)
			{
				for (shamt_c=min_cols; shamt_c<=max_cols; ++shamt_c)
				{
					// TODO: swap the order r<->c
					oPos = initialOffest + shamt_r+shamt_c*WINDOW_ROWS;
					pRes = &(pLine[oPos]);
					pRes->real += res.real;
					pRes->imag += res.imag;
				}
			}
		}
	}
	//	PRINT_BANNER
}

void RunCombinationDiagonalWritesEnhanced(unsigned int comb_id, ComplexFloat32** inputMat,ComplexFloat32** outputMat,int WINDOW_ROWS,int WINDOW_COLS, ComplexFloat32* outputVector)
{
	// extract data from the combinations table, just for convenience
	Combination* comb = &combinationsTable[comb_id];
	unsigned int mult1r = comb->mult1r;
	unsigned int mult1c = comb->mult1c;
	unsigned int mult2r = comb->mult2r;
	unsigned int mult2c = comb->mult2c;
	unsigned int blockR = comb->blockSize >> 8;
	unsigned int blockC = comb->blockSize & 0xFF;

	unsigned int nw_r, nw_c;
	int shamt_r, shamt_c;
	ComplexFloat32 res; // the multiplication result
	int m1r, m1c, m2r, m2c; // coordinates of multipliers
	int min_rows, min_cols, max_rows, max_cols; // shifting boundries

	int offset;

	int oPos;

	//    printf("%d %d %d\n", outputMat,WINDOW_ROWS,WINDOW_COLS);

	ComplexFloat32* pLine;
	ComplexFloat32* pRes;
	int initialOffest = 0; // this will be changed only when comb_id>=N (a.k.a set2)

	// calculate the initial position in outputMat for writes
	if (comb_id<WINDOW_ROWS*WINDOW_COLS) {
		//		pLine = outputMat[comb_id];
		pLine = outputVector;   
	}
	else {
		// the number of the diagonal (its distance from the main diagonal)
		//		pLine = outputMat[(mult2c-1)*WINDOW_ROWS + WINDOW_ROWS-mult1r];
		pLine = outputVector;   
		initialOffest = mult1r;
	}

	// iterate over all possible NW-positions for this combination
	for (nw_r=0; nw_r<INPUT_ROWS-blockR; ++nw_r)
	{
		for (nw_c=0; nw_c<INPUT_COLS-blockC; ++nw_c)
		{
			// calculate the initial multipliers' coordinates in the
			// input-image axis
			m1r = nw_r + mult1r;
			m1c = nw_c + mult1c;
			m2r = nw_r + mult2r;
			m2c = nw_c + mult2c;

			// make the multiplication
			multconj(&inputMat[m1r][m1c], &inputMat[m2r][m2c], &res);
			//            printf("%d %d %d %d\n",m1r,m1c,m2r,m2c);


			// calculate boundaries for shifts
			max_rows = MYMIN(WINDOW_ROWS - blockR - 1, nw_r);
			max_cols = MYMIN(WINDOW_COLS - blockC - 1, nw_c);
			min_rows = MYMAX(0, (int)nw_r + WINDOW_ROWS - INPUT_ROWS);
			min_cols = MYMAX(0, (int)nw_c + WINDOW_COLS - INPUT_COLS);

			// iterate over all possible shifts
			for (shamt_c=min_cols; shamt_c<=max_cols; ++shamt_c)
			{
				offset = initialOffest + shamt_c*WINDOW_ROWS;
				for (shamt_r=min_rows; shamt_r<=max_rows; ++shamt_r)
				{
					oPos = offset + shamt_r;
					pRes = &(pLine[oPos]);
					//					pRes = (pLine+oPos);

					pRes->real += res.real;
					pRes->imag += res.imag;
				}
			}
			/*
			//offset = initialOffest;
			for (shamt_c=min_cols; shamt_c<=max_cols; ++shamt_c)
			{
			//offset += shamt_c*WINDOW_ROWS;
			for (shamt_r=min_rows; shamt_r<=max_rows; ++shamt_r)
			{
			pRes = &(pLine[initialOffest + shamt_c+shamt_r*WINDOW_ROWS]);
			pRes->real += res.real;
			pRes->imag += res.imag;
			}
			}
			 */
		}
	}
	//	PRINT_BANNER
}


void RunCombinationDiagonalAddMul(unsigned int comb_id,ComplexFloat32** inputMat,int WINDOW_ROWS,int WINDOW_COLS,
		ComplexFloat32* outputVector,int* finalMatPosR, int* finalMatPosC)
{
	// extract data from the combinations table, just for convenience
	Combination* comb = &combinationsTable[comb_id];
	unsigned int mult1r = comb->mult1r;
	unsigned int mult1c = comb->mult1c;
	unsigned int mult2r = comb->mult2r;
	unsigned int mult2c = comb->mult2c;


	unsigned int deltaR = (int)abs((int)(mult1r-mult2r));
	unsigned int deltaC = (int)abs((int)(mult1c-mult2c));

	int type2 = comb->type2;
	//	unsigned int deltaR = comb->blockSize >> 8;
	//	unsigned int deltaC = comb->blockSize & 0xFF;

	unsigned int numElementsInBlock = PR-abs(deltaR);
	unsigned int numBlocks = PC-deltaC;

	int nr=NR,nc=NC,pr=PR,pc=PC;
	const unsigned int DR= nr-pr;
	const unsigned int DC= nc-pc;

	/*
	   ComplexFloat32 outputVector[WINDOW_ROWS*WINDOW_COLS];
	   int finalMatPosR[WINDOW_ROWS*WINDOW_COLS];
	   int finalMatPosC[WINDOW_ROWS*WINDOW_COLS];
	 */
	int elementC=0;
	ComplexFloat32 temp_res = {0,0}; // the multiplication result
	int i,j,r,c;

	memset(outputVector,0,WINDOW_ROWS*WINDOW_COLS*sizeof(ComplexFloat32));
	memset(finalMatPosR,0,WINDOW_ROWS*WINDOW_COLS*sizeof(int));
	memset(finalMatPosC,0,WINDOW_ROWS*WINDOW_COLS*sizeof(int));


	//	if(comb_id<WINDOW_ROWS*WINDOW_COLS)
	if (!type2)
	{
		// Computing the first index of the combination.
		// This index is made up
		for(i=0; i<=( DR); i++)
		{
			int iPdr=i+deltaR;
			for(j=0; j<= (DC); j++)
			{
				int jPdc = j+deltaC;
				ADD_MULT_CONJ(inputMat[i][j], inputMat[iPdr][jPdc], temp_res);
			}
		}
	}
	else
	{
		// Computing the first index of the combination.
		// This index is made up
		for(i=0; i<=( DR); i++)
		{
			int iPdr=i+deltaR;
			for(j=0; j<= (DC); j++)
			{
				int jPdc = j+deltaC;
				ADD_MULT_CONJ(inputMat[iPdr][j], inputMat[i][jPdc], temp_res);
			}
		}
	}

	//// Computing the first index of the combination.
	//// This index is made up
	//for(i=0; i<=( DR); i++)
	//{
	//	int iPdr=i+deltaR;
	//	for(j=0; j<= (DC); j++)
	//	{
	//		int jPdc = j+deltaC;
	//		ADD_MULT_CONJ(inputMat[i][j], inputMat[iPdr][jPdc], temp_res);
	//	}
	//}

	outputVector[0].real = temp_res.real;
	outputVector[0].imag = temp_res.imag;

	// Checking if the first element belongs to the first set of combinatons.
	// The combination that the first element is above the second.
	//	if(comb_id<WINDOW_ROWS*WINDOW_COLS)
	if (!type2)
	{
		finalMatPosR[0]=0;
		finalMatPosC[0]=PR*deltaC+deltaR;
		elementC++;
	}
	else
	{
		finalMatPosR[0]=deltaR;
		finalMatPosC[0]=PR*deltaC;
		elementC++;
	}

	for(r=1;r<(PR-deltaR);r++)
	{
		ComplexFloat32 newRowSum = {0,0},oldRowSum = {0,0};
		ComplexFloat32 addRows  = {0,0};

		int rM1 = r-1;
		int k = DR+1 + rM1;

		int kPdr = k+deltaR;
		int rM1Pdr = rM1 + deltaR;
		int cPdc;


		//		if(comb_id<WINDOW_ROWS*WINDOW_COLS)
		if (!type2)
		{
			for(c=0;c<=(DC); c++)
			{
				cPdc = c+deltaC;
				ADD_MULT_CONJ(inputMat[k][c],inputMat[kPdr][cPdc],newRowSum);
				ADD_MULT_CONJ(inputMat[rM1][c],inputMat[rM1Pdr][cPdc],oldRowSum);
			}
		}
		else
		{
			for(c=0;c<=(DC); c++)
			{
				cPdc = c+deltaC;
				ADD_MULT_CONJ(inputMat[kPdr][c],inputMat[k][cPdc],newRowSum);
				ADD_MULT_CONJ(inputMat[rM1Pdr][c],inputMat[rM1][cPdc],oldRowSum);
			}
		}

		ADD_REG(newRowSum,addRows)
			SUB_REG(oldRowSum,addRows)
			ADD_REG(outputVector[rM1],outputVector[r])
			ADD_REG(addRows,outputVector[r])

			finalMatPosR[elementC]=finalMatPosR[elementC-1]+1;;
		finalMatPosC[elementC]=finalMatPosC[elementC-1]+1;
		elementC++;
	}



	for(c=1; c<numBlocks; c++)
	{
		ComplexFloat32 newColSum = {0,0},oldColSum = {0,0};
		ComplexFloat32 addCols  = {0,0};

		int cM1 = c-1;
		int dcPc = DC+c;

		int w = DC+1 + cM1;

		int cM1PdeltaC = cM1 + deltaC;
		int dcPcPdeltaC = dcPc+deltaC;

		//		if(comb_id<WINDOW_ROWS*WINDOW_COLS)
		if (!type2)
		{
			for(r=0;r<=(DR); r++)
			{
				int rPdr = r+deltaR;
				ADD_MULT_CONJ(inputMat[r][w],inputMat[rPdr][w+deltaC],newColSum)
					ADD_MULT_CONJ(inputMat[r][cM1],inputMat[rPdr][cM1+deltaC],oldColSum)
			}
		}
		else
		{
			for(r=0;r<=(DR); r++)
			{
				int rPdr = r+deltaR;
				ADD_MULT_CONJ(inputMat[rPdr][w],inputMat[r][w+deltaC],newColSum)
					ADD_MULT_CONJ(inputMat[rPdr][cM1],inputMat[r][cM1+deltaC],oldColSum)
			}
		}

		ADD_REG(newColSum,addCols)
			SUB_REG(oldColSum,addCols)
			ADD_REG(outputVector[(c-1)*numElementsInBlock],outputVector[c*numElementsInBlock])
			ADD_REG(addCols,outputVector[c*numElementsInBlock])

			finalMatPosR[elementC]=finalMatPosR[elementC-numElementsInBlock]+PR;
		finalMatPosC[elementC]=finalMatPosC[elementC-numElementsInBlock]+PR;
		elementC++;


		for(r=1; r<numElementsInBlock; r++)
		{
			ComplexFloat32 w = {0},x = {0},y = {0},z = {0},deltaRowSum = {0};
			ComplexFloat32 tempRes = {0};
			int rM1 = r-1;
			int drPr = DR+r;

			int rM1PdeltaR = rM1 + deltaR;
			int drPrPdeltaR = drPr+deltaR;

			//			if(comb_id<WINDOW_ROWS*WINDOW_COLS)
			if (!type2)
			{
				MULT_CONJ(inputMat[rM1][cM1],inputMat[rM1PdeltaR][cM1PdeltaC],w);
				MULT_CONJ(inputMat[rM1][dcPc],inputMat[rM1PdeltaR][dcPcPdeltaC],x);
				MULT_CONJ(inputMat[drPr][cM1],inputMat[drPrPdeltaR][cM1PdeltaC],y);
				MULT_CONJ(inputMat[drPr][dcPc],inputMat[drPrPdeltaR][dcPcPdeltaC],z);
			}
			else
			{
				MULT_CONJ(inputMat[rM1PdeltaR][cM1],inputMat[rM1][cM1PdeltaC],w);
				MULT_CONJ(inputMat[rM1PdeltaR][dcPc],inputMat[rM1][dcPcPdeltaC],x);
				MULT_CONJ(inputMat[drPrPdeltaR][cM1],inputMat[drPr][cM1PdeltaC],y);
				MULT_CONJ(inputMat[drPrPdeltaR][dcPc],inputMat[drPr][dcPcPdeltaC],z);
			}

			ADD_REG(w,tempRes)
				SUB_REG(x,tempRes)
				SUB_REG(y,tempRes)
				ADD_REG(z,tempRes)

				ADD_REG(outputVector[(c-1)*numElementsInBlock + r],deltaRowSum)
				SUB_REG(outputVector[(c-1)*numElementsInBlock+ rM1],deltaRowSum)


				ADD_REG(deltaRowSum,tempRes)
				ADD_REG(outputVector[c*numElementsInBlock+rM1],tempRes)
				ADD_REG(tempRes,outputVector[c*numElementsInBlock+r])

				finalMatPosR[elementC]=finalMatPosR[elementC-1]+1;;
			finalMatPosC[elementC]=finalMatPosC[elementC-1]+1;
			elementC++;

		}
	}
	return;
	for(i=0; i<numElementsInBlock*numBlocks; i++)
	{
		outputMato[finalMatPosR[i]][finalMatPosC[i]].real=outputVector[i].real;
		outputMato[finalMatPosR[i]][finalMatPosC[i]].imag=outputVector[i].imag;
	}

}

void Compare2Mats()
{
	int i,j;
	double sum=0;
	double temp;
	double mat_sum=0;

	for (i=0; i<WINDOW_ROWS*WINDOW_COLS; i++) {
		for (j=i; j<WINDOW_ROWS*WINDOW_COLS; j++) {
			temp = (outputMat[i][j].real-outputMato[i][j].real)*(outputMat[i][j].real-outputMato[i][j].real)
				+ (outputMat[i][j].imag-outputMato[i][j].imag)*(outputMat[i][j].imag-outputMato[i][j].imag);
			sum+=temp;
			mat_sum +=(outputMat[i][j].real)*(outputMat[i][j].real)+ (outputMat[i][j].imag)*(outputMat[i][j].imag);
		}
	}
}
