#include "fastmath.h"

// Cawley, GC. "On a Fast Compact Approximation of the Exponential Function"
// Neural Computation 2000;12(9);2009-2012

double fastexp(double y) {
	union expun eco;

	eco.n.i = (int) (1512775 * y + 1072632447);
	eco.n.j = 0;

	return eco.d;
}
