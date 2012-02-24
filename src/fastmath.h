#include <math.h>

union expun {
	double d;
	struct { int j, i; } n;
};


double fastexp(double y);

