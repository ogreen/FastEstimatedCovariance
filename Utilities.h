#ifndef _Utilities_h_
#define _Utilities_h_

#include "Types.h"

#define MYMIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MYMAX(X,Y) ((X) > (Y) ? (X) : (Y))

// dst = conjugate(src)
//void conj_(const ComplexFloat32 *src, ComplexFloat32 *dst);

// dst = a*b for complex numbers
//void mult_(const ComplexFloat32 *a, const ComplexFloat32 *b, ComplexFloat32* dst);

// dst = a * conjugate(b) for complex numbers
void multconj(const ComplexFloat32 *a, const ComplexFloat32 *b, ComplexFloat32 *dst);

// dst += a * conjugate(b) for complex numbers
//void addmultconj(const ComplexFloat32 *a, const ComplexFloat32 *b, ComplexFloat32 *dst);

/* complex implementation of dst += src */
//void add_(const ComplexFloat32 *src, ComplexFloat32 *dst);

// clears top triangle of the output matrix
//void ZeroOutputMatrixTopTriag(ComplexFloat32** matToZero);

// clears the whole output matrix
//void ZeroOutputMatrix(ComplexFloat32** matToZero, int rows, int cols);

// compares the result with the supplied reference
void CompareResults(void);

// completes the output matrix by transposing it. M = M + M'
// this is done so that the comparision will not fail
void CompleteOutputByTransposition(void);


#define SUB_REG(src,dst) dst.real -= src.real;dst.imag -= src.imag;

#define ADD_REG(src,dst) dst.real += src.real; dst.imag += src.imag;

#define ADD_MULT_CONJ(a,b,dst) dst.real += a.real * b.real + a.imag * b.imag; 	dst.imag += b.real * a.imag - a.real * b.imag;

#define MULT_CONJ(a,b,dst) dst.real = a.real * b.real + a.imag * b.imag; 	dst.imag = b.real * a.imag - a.real * b.imag;



ComplexFloat32** createComplexArray(int rows,int cols);
void destroyComplexArray(ComplexFloat32** array,int rows);



#endif //_Utilities_h_
