#ifndef __SGN_H__
#define __SGN_H__
#include <stdio.h>
#include <stdlib.h>

namespace con2020 {

template <typename T> T sgn(T x) {
	return (x > 0) - (x < 0);
}

} // end namespace con2020
#endif