// C-compatible prototypes for the Fortran-based functions in LAPACK.
#ifndef CxLapack_h
#define CxLapack_h

void
dgeev_(char *jobvl, char *jobvr, int *n, double *A, int *lda, double *wr,
  double *wi, double *vl, int *ldvl, double *vr, int *ldvr, double *work,
  int *lwork, int *info);

void
dgetrf_(int *m, int *n, double *A, int *lda, int *ipiv, int *info);

void
dgetri_(int *n, double *A, int *lda, int *ipiv, double *work, int *lwork,
  int *info);

#endif // CxLapack_h
