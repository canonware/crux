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
    union {
	CxtDMDist d;
	int64_t i;
    } aU, bU;

    // Convert a and b to lexicographically ordered ints.
    aU.d = a;
    if (aU.i < 0) {
	aU.i = 0x8000000000000000ULL - aU.i;
    }

    bU.d = b;
    if (bU.i < 0) {
	bU.i = 0x8000000000000000ULL - bU.i;
    }

    // Check if a and b are within aMaxUlps of each other.
    if (llabs(aU.i - bU.i) <= CxmDistMatrixNjMaxUlps) {
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
