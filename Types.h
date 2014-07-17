/**
  Types.h - types and data structures
**/

#ifndef _Types_h_
#define _Types_h_

typedef float Float32;
typedef unsigned char Uint8;
typedef unsigned short Uint16;
typedef unsigned int Uint32;

typedef struct ComplexFloat32_t {
	Float32 real;
	Float32 imag;
} ComplexFloat32;

typedef struct Size_t {
	int rows;
	int cols;
} Size;

typedef struct Coordinate_t {
	int row;
	int col;
} Coordinate;

/**
this structure supports window sizes up to 181x181
because refCoords is 16bits and should be able to store
integers of up to 65522 = 2*(181^2) < 2^16 = 65536

for each size stored in the struct, the first half of its size in
bytes (MSB) stores the row property and the last half (LSB) stores
the column data.
for instance, m1 = Combination.mult1 is sizes 16 byes, and represent
the coordiantes of the first multiplier, so (m1 >> 8) is the row
index of the first multiplier, and (m1 & 0xFF) is the column index.

also, disregard data alignment to preserve RAM (this is critical
on plurality, but might make things slower on x86)
*/
#ifdef _WIN32
#pragma pack(push,1)
#endif
typedef struct {
    int mult1r;
    int mult1c;
    int mult2r;
    int mult2c;
    short blockSize; // block size minus one, e.g. for mult1=(0,0) and mult2=(0,0), block size will be (0,0)
    short refCoords; // the result coordinates for the zero-iteration. i.e. without any shifts
	int type2; // 0 - for the first P*P. 1 - for the next (P-1)(P-1)
	int id;
} Combination;

#endif // _Types_h_
