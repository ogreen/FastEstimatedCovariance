// Configuration file

#ifndef _Config_h_
#define _Config_h_

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// User settable parameters
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// print an info line each time a function ends, for debug purposes
#define PRINT_FINISHED_FUNC_NAMES 0

// print the result of the matrix-comparisons in case of success as well
// if set to 0, only failed comparisons will cause notification
#define PRINT_CONPARISONS 0

// number of iterations for the loop to run - more iterations means more
// precise results, since they are averaged


// note that the images should be interleaved (complex), single
// precision (float32), binary file i.e first 32 bytes are
// a[0][0].real, next 32 bytes are a[0][0].imag, next 32 bytes are
// a[0][1].real, and so on.
// mind the endianess! 'le' for x86, 'be' for the simulator
// if INPUT_IMAGE_FILENAME is not set, the algorithm will use dummy data
//#define INPUT_IMAGE_FILENAME "image_le.dat"
#define OUTPUT_FILENAME "out_le.dat"

// if this is defined, the results will be compared to the reference in this
// filename once the run is over. in case this isn't defined, the comparisons
// will always return true
//#define OUTPUT_REF_FILENAME "ref_le.dat"

// maximum allowed delta between elements - this will affect only the
// comparison result
// specify the maximum allowed delta in percents!
//#define MAX_DELTA 1e-3f
#define MAX_DEVIATION 5


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Don't change below this line
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// number of combinations
#define COMBINATIONS_COUNT (WINDOW_ROWS*WINDOW_COLS+(WINDOW_ROWS-1)*(WINDOW_COLS-1))

// useful macro
#define PRINT_BANNER \
		if (PRINT_FINISHED_FUNC_NAMES) printf("Finished %s.\n", __FUNCTION__);

#endif // _Config_h_
