#ifndef _BuildCovMat_h_
#define _BuildCovMat_h_

#include <omp.h>

#include "Types.h"
//#include "configsizes.cfg"

/////////////////////////////////////////////////////////////////////////////
// NAME: BuildCovMatNaive
// DESCRIPTION:	The naive way of building the covariance estimation matrix.
//				Runs over the input matrix, and for each position of the
//				sliding window, uses column-stacking, inner-multiplying and
//				later summation of the results.
//
// INPUTS: none
// OUTPUTS: none
/////////////////////////////////////////////////////////////////////////////
double BuildCovMatNaive(void);

/////////////////////////////////////////////////////////////////////////////
// NAME: BuildCovMatParallel
// DESCRIPTION:	Our original efficiant and parallel algorithm for building
//				the	the covariance matrix. iterates over all NW coordinates
//				and in the input matrix, and calls RunCombiantion() for
//				each combination, on each NW.
//				This version uses the 2D combinationsTable, without any
//				optimizations for runtime speedup.
//				The upper triangle of the output matrix must be zero prior to
//				calling this function.
//
// INPUTS: none
// OUTPUTS: none
/////////////////////////////////////////////////////////////////////////////
double  BuildCovMatParallel(void);

/////////////////////////////////////////////////////////////////////////////
// NAME: RunCombination
// DESCRIPTION:	This is what each combination does in the parallel algorithm.
//				Essentially, it parses the line comb_id in the combinations
//				table, and run over all valid NW indices in the input
//				matrix. In each iteration, it makes the multiplication once,
//				and calculates all valid shifts and target indices to append
//				the result to.
//
// INPUTS: none
// OUTPUTS: none
/////////////////////////////////////////////////////////////////////////////
void RunCombination(unsigned int comb_id, ComplexFloat32** inputMat,ComplexFloat32** outputMat,int WINDOW_ROWS,int WINDOW_COLS);

/////////////////////////////////////////////////////////////////////////////
// NAME: BuildCovMatDiagonalWrites
// DESCRIPTION:	First optimization for the parallel algo. This function does
//				not writes its output into C (the algorithm's output matrix).
//				Instead, each RunCombination() gets as argument a vector,
//				which represents the diagonal that combination writes to, and
//				writes the output to that diagonal.
//				At the end of the run (which can also be done off-line), the
//				full "right" output image is reordered by another function.
//				Also refer to the documentation of the function
//				RunCombinationDiagonalWrites().
//
// INPUTS: none
// OUTPUTS: none
/////////////////////////////////////////////////////////////////////////////
double  BuildCovMatDiagonalWrites(void);

/////////////////////////////////////////////////////////////////////////////
// NAME: BuildCovMatDiagonalWritesEnhanced
// DESCRIPTION:	Exactly identical to BuildCovMatDiagonalWrites, but calls
//				RunCombinationDiagonalWritesEnhanced() instead of
//				RunCombinationDiagonalWrites().
//
// INPUTS: none
// OUTPUTS: none
/////////////////////////////////////////////////////////////////////////////
double BuildCovMatDiagonalWritesEnhanced();

double BuildCovMatDiagonalJustMultiplicatedEnhancedParallel(int doSort, int doBuildTable) ;


/////////////////////////////////////////////////////////////////////////////
// NAME: BuildCovMatDiagonalAddMul
// DESCRIPTION:	Exactly identical to BuildCovMatDiagonalWrites, but calls
//				RunCombinationDiagonalWritesEnhanced() instead of
//				RunCombinationDiagonalWrites().
//
// INPUTS: none
// OUTPUTS: none
/////////////////////////////////////////////////////////////////////////////
double BuildCovMatDiagonalAddMul(void);

double BuildCovMatDiagonalAddMulOPENMP(int doSort, int doBuildTable);




/////////////////////////////////////////////////////////////////////////////
// NAME: RunCombinationDiagonalWrites
// DESCRIPTION:	Each combination, when called by BuildCovMatDiagonalWrites,
//				is executed (via this function) in the same way described by
//				the RunCombination() function, but with one difference: the
//				combination does not write to the "right" indices in the
//				output matrix, but instead writes to a continuous area in it.
//				the rule is simple:
//				Say the output matrix (outputMat) is sized NxN.
//				combination i<N writes its outputs to line i of outputMat.
//				combination i>=N writes its output to line X, where X=[the
//				column-stack-distance between the combination's multipliers].
//
//				for more details, read the code.
//
//				TODO: maybe putting this in the combiantions table will be
//				faster?
//
// INPUTS: none
// OUTPUTS: none
/////////////////////////////////////////////////////////////////////////////
void RunCombinationDiagonalWrites(unsigned int comb_id, ComplexFloat32** inputMat,ComplexFloat32** outputMat,int WINDOW_ROWS,int WINDOW_COLS);

/////////////////////////////////////////////////////////////////////////////
// NAME: RunCombinationDiagonalWritesEnhanced
// DESCRIPTION:	Exactly identical to RunCombinationDiagonalWrites, with two
//				changes:
//				1. The inner nesting loop order is swapped, i.e. running on
//				rows first and then on cols.
//				2. In those loops, some calculation were factored out, so
//				that it would be calculated less time.
// INPUTS: none
// OUTPUTS: none
/////////////////////////////////////////////////////////////////////////////
void RunCombinationDiagonalWritesEnhanced(unsigned int, ComplexFloat32** inputMat,ComplexFloat32** outputMat,int WINDOW_ROWS,int WINDOW_COLS,ComplexFloat32* outputVector);




/////////////////////////////////////////////////////////////////////////////
// NAME: RunCombinationDiagonalWritesEnhanced
// DESCRIPTION:	Exactly identical to RunCombinationDiagonalWrites, with two
//				changes:
//				1. The inner nesting loop order is swapped, i.e. running on
//				rows first and then on cols.
//				2. In those loops, some calculation were factored out, so
//				that it would be calculated less time.
// INPUTS: none
// OUTPUTS: none
/////////////////////////////////////////////////////////////////////////////
void RunCombinationDiagonalAddMul(unsigned int,ComplexFloat32** inputMat,int WINDOW_ROWS,int WINDOW_COLS,
                                ComplexFloat32* outputVector,int* finalMatPosR, int* finalMatPosC);





/////////////////////////////////////////////////////////////////////////////
// NAME: RemapDiagonalWritesOutputMatrix
// DESCRIPTION:	As explained before, BuildCovMatDiagonalWrites() function
//				writes to outputMat, but not to the "right" places. This
//				function transforms outputMat to the "right" way it should be,
//				so comparision can be made between other outputMat's, as they
//				are composed by the other algorithms.
//
// INPUTS: none
// OUTPUTS: none
/////////////////////////////////////////////////////////////////////////////
void RemapDiagonalWritesOutputMatrix(void);

/////////////////////////////////////////////////////////////////////////////
// NAME: BuildCombinationsTable
// DESCRIPTION:
//
//
// INPUTS: none
// OUTPUTS: none
/////////////////////////////////////////////////////////////////////////////
void BuildCombinationsTable(int WINDOW_ROWS,int WINDOW_COLS);



	void Compare2Mats();

#endif // _BuildCov_h_


