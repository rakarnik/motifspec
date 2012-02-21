#include "fastmath.h"

// Cawley, GC. "On a Fast Compact Approximation of the Exponential Function"
// Neural Computation 2000;12(9);2009-2012

#define EXP_A (1048576/M_LN2)
#define EXP_C 60801

double fastexp(double y) {
	union {
		double d;
		struct { int j, i; } n;
	} eco;

	eco.n.i = (int) (EXP_A * (y)) + (1072693248 - EXP_C);
	eco.n.j = 0;

	return eco.d;
}
