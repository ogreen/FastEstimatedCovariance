

#ifndef _Globals_h_
#define _Globals_h_

#include "Types.h"
#include "Config.h"

int INPUT_ROWS;
int INPUT_COLS;

int WINDOW_ROWS;
int WINDOW_COLS;

int NR;
int NC;

int PR;
int PC;

int ITERATIONS;
int THREADS;

int BENCHMARK_FLAGS;
int SORT_FLAG;

#define COMBINATIONS_COUNT (WINDOW_ROWS*WINDOW_COLS+(WINDOW_ROWS-1)*(WINDOW_COLS-1))




// the input image
ComplexFloat32** inputMat;//[INPUT_ROWS][INPUT_COLS];

// the window vector which is created by column-stacking each window in
// the serial algorithm.
// it is also used for the VectorWrites implementation
ComplexFloat32* windowVec;//[WINDOW_ROWS*WINDOW_COLS];

// the output matrix
ComplexFloat32** outputMat;//[WINDOW_ROWS*WINDOW_COLS][WINDOW_ROWS*WINDOW_COLS];

// the output matrix2
ComplexFloat32** outputMato;//[WINDOW_ROWS*WINDOW_COLS][WINDOW_ROWS*WINDOW_COLS];


// combinations table
Combination* combinationsTable;//[COMBINATIONS_COUNT];

//Uint16 OffsetsTable[PERMUTATIONS_COUNT][SUBAPSIZER][SUBAPSIZEC][2];

#endif //_Globals_h_
