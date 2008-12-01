#ifndef CxDistMatrixNj_h
#define CxDistMatrixNj_h

#include "../Cx.h"
#include "CxDistMatrix.h"

#ifndef CxmUseInlines
int
CxDistMatrixNjDistCompare(CxtDMDist a, CxtDMDist b);
#endif

#if (defined(CxmUseInlines) || defined(CxDistMatrixNj_c))
// Compare two distances, and consider them equal if they are close enough.
// CxmDistMatrixNjMaxUlps (ULP: Units in Last Place) specifies the maximum
// number of ulps that a and b may differ by and still be considered equal.
// Ideally, we would keep track of the maximum possible accumulated error
// during tree construction, but the error grows exponentially, so that
// approach ends up not being practical.
//
// Provide 9 (base 10) digits of accuracy.
#define CxmDistMatrixNjMaxUlps 0x3fffffU
CxmInline int
CxDistMatrixNjDistCompare(CxtDMDist a, CxtDMDist b) {
    int ret;
    int64_t aInt, bInt;

    // Convert a and b to lexicographically ordered ints.
    aInt = *(int64_t *) &a;
    if (aInt < 0) {
	aInt = 0x8000000000000000ULL - aInt;
    }

    bInt = *(int64_t *) &b;
    if (bInt < 0) {
	bInt = 0x8000000000000000ULL - bInt;
    }

    // Check if a and b are within aMaxUlps of each other.
    if (llabs(aInt - bInt) <= CxmDistMatrixNjMaxUlps) {
	ret = 0;
    } else if (a < b) {
	ret = -1;
    } else {
	ret = 1;
    }

    return ret;
}
#endif

#endif // CxDistMatrixNj_h
