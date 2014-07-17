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



#include "Config.h"
#include "Types.h"
#include "Globals.h"
#include <string.h>
#include <stdio.h>

// note that indices here starts at 0, while in Matlab they start at 1
void BuildCombinationsTable(int WINDOW_ROWS,int WINDOW_COLS)
{
	int idx_row,idx_col;
	int comb_idx = 0;
	Combination comb;

	// first part
	// the first element is [0,0] and the second one travels
	for (idx_col=0; idx_col<WINDOW_COLS; ++idx_col)
	{
		for (idx_row=0; idx_row<WINDOW_ROWS; ++idx_row)
		{

			comb.mult1r=0;
			comb.mult1c=0;
			comb.mult2r=idx_row;
			comb.mult2c=idx_col;

			comb.blockSize = idx_row<<8 | idx_col;
			comb.refCoords = 0<<16 | (idx_row+(idx_col*WINDOW_ROWS));
			comb.type2 = 0;
			comb.id = comb_idx;
			memcpy(&combinationsTable[comb_idx++], &comb, sizeof(Combination));
		}
	}

	// second part
	// first element is on the leftmost column, and the
	// second is on the topmost one
	for (idx_row=1; idx_row<WINDOW_ROWS; ++idx_row)
	{
		for (idx_col=1; idx_col<WINDOW_COLS; ++idx_col)
		{
			comb.mult1r=idx_row;
			comb.mult1c=0;
			comb.mult2r=0;
			comb.mult2c=idx_col;

			comb.blockSize = idx_row<<8 | idx_col;
			comb.refCoords = idx_row<<16 | (idx_col*WINDOW_ROWS);
			comb.type2 = 1;
			comb.id = comb_idx;

			memcpy(&combinationsTable[comb_idx++], &comb, sizeof(Combination));
		}
	}
//	PRINT_BANNER
}
