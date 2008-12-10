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
// Provide 5 (base 10) digits of accuracy.
#define CxmDistMatrixNjMaxUlps 0x7fU
CxmInline int
CxDistMatrixNjDistCompare(CxtDMDist a, CxtDMDist b) {
    int ret;
    union {
	CxtDMDist d;
	int32_t i;
    } aU, bU;

    CxmAssert(sizeof(CxtDMDist) == sizeof(int32_t));

    // Convert a and b to lexicographically ordered ints.
    aU.d = a;
    if (aU.i < 0) {
	aU.i = 0x80000000U - aU.i;
    }

    bU.d = b;
    if (bU.i < 0) {
	bU.i = 0x80000000U - bU.i;
    }

    // Check if a and b are within aMaxUlps of each other.
    if (abs(aU.i - bU.i) <= CxmDistMatrixNjMaxUlps) {
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
