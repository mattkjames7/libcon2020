#include "polyeval.h"

namespace con2020 {

double polyeval(double x, double *c, int d) {
	double y = c[0];
	int i;
	for (i=1;i<=d;i++) {
		y = y*x + c[i];
	}
	return y;
}

double pol1eval(double x, double *c, int d) {
	double y = x + c[0];
	int i;
	for (i=1;i<d;i++) {
		y = y*x + c[i];
	}
	return y;
}

} // end namespace con2020