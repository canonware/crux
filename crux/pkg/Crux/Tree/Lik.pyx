from Crux.Character cimport Character

from libc cimport *
from CxMat cimport CxMatQDecomp, CxMatPt

cdef class QMat:
    """
        Partially initialize the internals of a Q matrix, which is conceptually
        a product of the following form:

            Q = R * Pi
                 _            _     _            _
                |  -  a  b  c  |   | pA  0  0  0  |
                |              |   |              |
                |  a  -  d  e  |   |  0 pC  0  0  |
              = |              | * |              |,
                |  b  d  -  f  |   |  0  0 pG  0  |
                |              |   |              |
                |_ c  e  f  - _|   |_ 0  0  0 pT _|

        where {a,b,c,d,e,f} are relative substitution rates, and {pA,pC,pG,pT}
        are frequencies.

        During Q matrix computation:
        * The elements on R's diagonal are set such that each row in R sums to
          0.
        * The relative rates are rescaled to fix one rate to 1.0 (f in the
          example above).
        * The frequencies are rescaled to sum to 1.

        By default:
        * The relative rates are all equal (1.0).
        * The frequencies are all equal (0.25 in the example above).

        Note that the QMat internals are lazily updated as necessary when used
        by PMat.  Such updates perform eigenvalue/eigenvector matrix
        decomposition of the Q matrix, which are used by PMat to compute
        substitution probabilities.
    """

    def __cinit__(self):
        self.qMat = NULL
        self.rMat = NULL
        self.piMat = NULL
        self.eigVecCube = NULL
        self.eigVals = NULL

    def __dealloc__(self):
        if self.qMat != NULL:
            free(self.qMat)
            self.qMat = NULL
        if self.rMat != NULL:
            free(self.rMat)
            self.rMat = NULL
        if self.piMat != NULL:
            free(self.piMat)
            self.piMat = NULL
        if self.eigVecCube != NULL:
            free(self.eigVecCube)
            self.eigVecCube = NULL
        if self.eigVals != NULL:
            free(self.eigVals)
            self.eigVals = NULL

    def __init__(self, type charType):
        cdef unsigned i, j

        self.up2date = False
        self.char_ = charType.get()
        self.dim = self.char_.nstates()

        self.rMat = <double *>calloc(self.dim * self.dim, sizeof(double))
        if self.rMat == NULL:
            raise MemoryError("Error allocating %dx%d rMat" % \
              (self.dim, self.dim))
        # Initialize to default equal relative rates.
        for 0 <= i < self.dim:
            self.rMat[i * self.dim + i] = 1.0

        self.piMat = <double *>calloc(self.dim * self.dim, sizeof(double))
        if self.piMat == NULL:
            raise MemoryError("Error allocating %dx%d piMat" % \
              (self.dim, self.dim))
        # Initialize to default equal frequencies.
        for 0 <= i < self.dim:
            for i < j < self.dim:
                self.piMat[i * self.dim + j] = 1.0 / <double>self.dim

        self.qMat = <double *>malloc(self.dim * self.dim * sizeof(double))
        if self.qMat == NULL:
            raise MemoryError("Error allocating %dx%d qMat" % \
              (self.dim, self.dim))

        self.eigVecCube = <double *>malloc(self.dim * self.dim * self.dim * \
          sizeof(double))
        if self.eigVecCube == NULL:
            raise MemoryError("Error allocating %d-element eigVecCube" % \
              self.dim)

        self.eigVals = <double *>malloc(self.dim * sizeof(double))
        if self.eigVals == NULL:
            raise MemoryError("Error allocating %d-element eigVals" % self.dim)

    cdef double getRate(self, unsigned i, unsigned j) except -1.0:
        assert i < self.dim
        assert j < self.dim
        assert i != j

        # Only the upper triangle of rMat is maintained, so swap i and j if
        # necessary.
        if i > j:
            i, j = j, i

        return self.rMat[i * self.dim + j]

    cdef void setRate(self, unsigned i, unsigned j, double rate) except *:
        assert i < self.dim
        assert j < self.dim
        assert i != j

        # Only the upper triangle of rMat is maintained, so swap i and j if
        # necessary.
        if i > j:
            i, j = j, i

        self.rMat[i * self.dim + j] = rate

    cdef double getFreq(self, unsigned i) except -1.0:
        assert i < self.dim

        return self.piMat[i * self.dim + i]

    cdef void setFreq(self, unsigned i, double freq) except *:
        assert i < self.dim

        self.piMat[i * self.dim + i] = freq

    cdef void _update(self) except *:
        cdef double fixedRate, piSum, elm
        cdef unsigned i, j

        assert not self.up2date

        # Rescale rMat, if necessary.
        fixedRate = self.rMat[(self.dim-1)*self.dim - 1]
        if fixedRate != 1.0:
            for 0 <= i < self.dim:
                for i < j < self.dim:
                    self.rMat[i*self.dim + j] /= fixedRate

        # Rescale piMat, if necessary.
        piSum = 0.0
        for 0 <= i < self.dim:
            piSum += self.piMat[i*self.dim + i]
        if piSum != 1.0:
            self.piMat[i*self.dim + i] /= piSum

        # Fill in the upper triangle of qMat.
        for 0 <= i < self.dim:
            self.qMat[i*self.dim + i] = 0.0
        for 0 <= i < self.dim:
            for i < j < self.dim:
                elm =  self.rMat[i*self.dim + j] * self.piMat[j*self.dim + j]
                self.qMat[i*self.dim + j] = elm
                self.qMat[i*self.dim + i] -= elm
                self.qMat[j*self.dim + j] -= elm

        # Decompose qMat.
        CxMatQDecomp(self.dim, self.qMat, self.eigVecCube, self.eigVals)
        self.up2date = True

cdef class PMat:
    """
        Partially initialize a substitution probability matrix.  The matrix is
        not useful until the compute() method is called.
    """

    def __cinit__(self):
        self.pMat = NULL
        self.eigValsExp = NULL

    def __dealloc__(self):
        if self.pMat != NULL:
            free(self.pMat)
            self.pMat = NULL

    def __init__(self, type charType):
        self.char_ = charType.get()
        self.dim = self.char_.nstates()

        self.pMat = <double *>malloc(self.dim * self.dim * sizeof(double))
        if self.pMat == NULL:
            raise MemoryError("Error allocating %dx%d pMat" % \
              (self.dim, self.dim))

    cdef void compute(self, QMat Q, double mu, double t) except *:
        """
            Compute substitution probabilities, based on the following formula:

                      Q*mu*t
              P(t) = e
        """
        if not Q.up2date:
            Q._update()
        CxMatPt(self.dim, self.pMat, Q.eigVecCube, Q.eigVals, mu * t);
