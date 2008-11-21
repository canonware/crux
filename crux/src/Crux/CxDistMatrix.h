#ifndef CxDistMatrix_h
#define CxDistMatrix_h

#include "Cx.h"

typedef double CxtDMDist;
typedef unsigned long CxtDMSize;

#ifndef CxmUseInlines
CxtDMSize
CxDistMatrixNxy2i(CxtDMSize n, CxtDMSize x, CxtDMSize y);
#endif

#if (defined(CxmUseInlines) || defined(CxDistMatrix_c))
CxmInline CxtDMSize
CxDistMatrixNxy2i(CxtDMSize n, CxtDMSize x, CxtDMSize y) {
    CxmAssert(x < n);
    CxmAssert(y < n);
    CxmAssert(x != y);

    if (x > y) {
	CxtDMSize t;

	t = x;
	x = y;
	y = t;
    }

    return n*x + y - (((x+3)*x) >> 1) - 1;
}
#endif

CxtDMDist *
CxDistMatrixNew(CxtDMSize ntaxa);

void
CxDistMatrixDelete(CxtDMDist *matrix);

#endif // CxDistMatrix_h
