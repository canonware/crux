#define CxDistMatrix_c
#include "CxDistMatrix.h"

CxtDMDist *
CxDistMatrixNew(CxtDMSize ntaxa) {
    return (CxtDMDist *) calloc(CxDistMatrixNxy2i(ntaxa, ntaxa - 2, ntaxa - 1),
      sizeof(CxtDMDist));
}

void
CxDistMatrixDelete(CxtDMDist *matrix) {
    free(matrix);
}
