import Crux.Exception

class Exception(Crux.Exception.Exception):
    pass

import exceptions

class ValueError(Exception, exceptions.ValueError):
    def __init__(self, str):
        self._str = str

    def __str__(self):
        return self._str

cimport Parsing
from Crux.Taxa cimport Taxon
cimport Crux.Taxa as Taxa
cimport Crux.Fasta as Fasta
from Crux.Character cimport Character, Dna, Protein
from Crux.DistMatrix cimport DistMatrix

from libc cimport *
from libm cimport *
from atlas cimport *
from CxMat cimport CxMatLogDet
from CxMath cimport pop, gcd

cdef extern from "Python.h":
    cdef object PyString_FromStringAndSize(char *s, Py_ssize_t len)

import sys

global __name__

cdef class CTMatrix:
    def __init__(self, input=None, type charType=Dna, taxaMap=None):
        assert issubclass(charType, Character)

        self.charType = charType

        # This is used to manage the rows of the data matrix (_taxonData).
        if taxaMap == None:
            taxaMap = Taxa.Map()
        self.taxaMap = taxaMap

        # Initialize serial number.  0 is skipped here so that it can be used
        # as a special value by users of CTMatrix.
        self.sn = 1

        # Row-major character data.  Each element in _taxonData is a string of
        # characters that belong to the corresponding taxon in taxaMap.  The
        # keys are the integer indices, as reported by self.taxaMap.indGet().
        self._taxonData = {}

        if input is not None:
            self._fastaNew(input, charType)

    cpdef _fastaNew(self, input, type charType):
        cdef Fasta.Parser parser

        parser = Fasta.Parser(self, self.taxaMap)
        parser.parse(input, charType)

    cpdef str fastaPrint(self, file outFile=None):
        cdef list lines, taxa
        cdef Taxon taxon
        cdef str s, taxonData
        cdef int i

        if outFile is None:
            lines = []

        # Print.
        taxa = self.taxaMap.taxaGet()
        for taxon in taxa:
            if self.dataGet(taxon) != None:
                s = ">%s\n" % taxon.label.replace(' ', '_')
                if outFile is None:
                    lines.append(s)
                else:
                    outFile.write(s)

                # Break into lines of length 75.
                taxonData = self.dataGet(taxon)
                for 0 <= i < len(taxonData) by 75:
                    if i + 75 < len(taxonData):
                        s = "%s\n" % taxonData[i:i+75]
                    else:
                        s = "%s\n" % taxonData[i:]

                    if outFile is None:
                        lines.append(s)
                    else:
                        outFile.write(s)

        if outFile is None:
            if len(lines) > 0:
                # Chop off the last '\n'.
                lines[-1] = lines[-1][:-1]
            return "".join(lines)
        else:
            return None

    # Return the character data for a taxon.
    cpdef str dataGet(self, Taxon taxon):
        if not self._taxonData.has_key(self.taxaMap.indGet(taxon)):
            rVal = None
        else:
            rVal = self._taxonData[self.taxaMap.indGet(taxon)]

        return rVal

    # Set the character data for a taxon.
    cpdef dataSet(self, Taxon taxon, str data):
        if (self.taxaMap.indGet(taxon) == None):
            raise ValueError("Taxon %r not in taxa map" % taxon.label)

        self._taxonData[self.taxaMap.indGet(taxon)] = data

        self.sn += 1

cdef struct PctIdentRel:
    bint isValid # True iff the character pair constitutes a valid site pair.
    float fident # 0.0 for unequal sites, 1.0 for identical sites, something
                 # in between for ambiguous sites.

cdef class PctIdent:
    cdef PctIdentRel tab[128 * 128]

    def __init__(self, type charType, bint avgAmbigs, bint scoreGaps):
        cdef Character char_
        cdef list codes
        cdef str iCode, jCode
        cdef char *iS, *jS
        cdef int iC, jC
        cdef int ambiguous, missing, i, j, iVal, jVal

        assert issubclass(charType, Character)
        char_ = charType.get()

        codes = char_.codes()

        ambiguous = char_.any
        missing = 0
        for 0 <= i < len(codes):
            iCode = <str>codes[i]
            iVal = char_.code2val(iCode)
            if scoreGaps and iVal == missing:
                iVal = ambiguous
            iS = iCode
            iC = iS[0]
            for 0 <= j < len(codes):
                jCode = <str>codes[j]
                jVal = char_.code2val(jCode)
                if scoreGaps and jVal == missing:
                    jVal = ambiguous
                jS = jCode
                jC = <char>jS[0]
                if iVal == missing or jVal == missing:
                    self.tab[(iC << 7) + jC].isValid = False
                else:
                    self.tab[(iC << 7) + jC].isValid = True
                    if avgAmbigs:
                        self.tab[(iC << 7) + jC].fident = \
                          <float>pop(iVal & jVal) / <float>(pop(iVal) \
                          * pop(jVal))
                    else:
                        if iVal == jVal:
                            self.tab[(iC << 7) + jC].fident = 1.0
                        else:
                            self.tab[(iC << 7) + jC].fident = 0.0

    # Compute comparison statistics for a vs. b.
    cdef void stats(self, char *a, char *b, unsigned seqlen,
      float *fident, unsigned *nvalid) except *:
        cdef float ident
        cdef unsigned valid, j, elm
        cdef int aC, bC

        ident = 0.0
        valid = 0
        for 0 <= j < seqlen:
            aC = a[j]
            bC = b[j]
            elm = (aC << 7) + bC
            if self.tab[elm].isValid:
                ident += self.tab[elm].fident
                valid += 1
        fident[0] = ident
        nvalid[0] = valid

cdef class K2p:
    cdef bint scoreGaps
    cdef int code2val[128]

    def __init__(self, bint scoreGaps):
        cdef Dna dna
        cdef str code
        cdef char *sCode
        cdef int iCode

        self.scoreGaps = scoreGaps

        dna = Dna.get()
        for code in dna.codes():
            sCode = <char *>code
            iCode = sCode[0]
            self.code2val[iCode] = dna.code2val(code)

    cdef void stats(self, char *a, char *b, unsigned len, unsigned *rN,
      unsigned *rKs, unsigned *rKv) except *:
        cdef unsigned j
        cdef int aC, bC, aI, bI
        cdef int ambiguous, aVal, bVal, aOff, bOff
        cdef double A[16]
        cdef double p, n, ks, kv, k

        ambiguous = Dna.get().any
        memset(A, 0, sizeof(A))
        for 0 <= j < len:
            aC = a[j]
            bC = b[j]
            aVal = self.code2val[aC]
            aPop = pop(aVal)
            if aC == bC and aPop== 1:
                aOff = ffs(aVal) - 1
                # Identical non-ambiguous characters.
                A[aOff*4 + aOff] += 1.0
            else:
                bVal = self.code2val[bC]
                bPop = pop(bVal)

                if aPop== 1 and bPop == 1:
                    aOff = ffs(aVal) - 1
                    bOff = ffs(bVal) - 1
                    A[aOff*4 + bOff] += 1.0
                elif not self.scoreGaps and (aPop== 0 or bPop == 0):
                    # Ignore this character.
                    pass
                else:
                    # Ambiguous character(s).

                    # Translate gaps to mean "any character".
                    if aVal == 0:
                        aVal = ambiguous
                        aPop= pop(aVal)
                    if bVal == 0:
                        bVal = ambiguous
                        bPop = pop(bVal)

                    # Each possible state combination is assigned an equal
                    # probability.  All combinations sum to 1.0.
                    p = <double>1.0 / <double>(aPop * bPop)

                    # Iterate over every state combination, and add p to the
                    # combinations that are represented by the ambiguity.
                    for 0 <= aI < 4:
                        for 0 <= bI < 4:
                            if (aVal & (1 << aI)) and (bVal & (1 << bI)):
                                A[aI*4 + bI] += p

        # Compute stats needed by the k2p correction method.
        #
        #    A C G T
        #  /--------  /------------
        # A| - v s v  |  0  1  2  3
        # C| v - v s  |  4  5  6  7
        # G| s v - v  |  8  9 10 11
        # T| v s v -  | 12 13 14 15
        assert self.code2val[c'A'] == 8
        assert self.code2val[c'C'] == 4
        assert self.code2val[c'G'] == 2
        assert self.code2val[c'T'] == 1
        ks = A[2] + A[7] + A[8] + A[13]
        kv = A[1] + A[3] + A[4] + A[6] + A[9] + A[11] + A[12] + A[14]
        k = ks + kv
        n = k + A[0] + A[5] + A[10] + A[15]

        # Tajima's Taylor expansion method for computing corrected distances
        # requires integer values.  Assuming it were possible to complete the
        # computation with real numbers, this rounding would be a minor
        # (unbiased) source of error for all but the shortest of sequences.
        #
        # Note that care is taken to do even rounding for ks and kv such that
        # they sum to the rounded k.
        rN[0] = <unsigned>round(n)
        rKs[0] = <unsigned>round(ks)
        rKv[0] = <unsigned>round(k - ks)

    #                           ^
    # Compute the first term of d, which can be factored as:
    #
    #     k
    #  ------          max
    #  \               ------   (i-j)    j   j
    #   \      i!      \      kv        2  ks
    #    >   ------- *  >     ------- * ------
    #   /         i    /      (i-j)!    j!
    #  /     2 i n     ------
    #  ------          j=min
    #   i=1
    cdef double dist1(self, unsigned n, unsigned ks, unsigned kv) except -1.0:
        cdef double ret, t, u
        cdef unsigned h, i, j, k, min, max

        ret = 0.0
        k = ks + kv
        for 1 <= i <= k:
            t = 0.5 / <double>i
            for 1 <= h <= i:
                t *= <double>h / <double>n

            min = (0 if i < kv else i - kv)
            max = (i if i < ks else ks)
            # Compute first iteration.
            j = min
            for 1 <= h <= i-j:
                t *= <double>kv
                t /= <double>h
            u = 2.0 * <double>ks
            for 1 <= h <= j:
                t *= u
                t /= <double>h
            ret += t

            # Compute remaining iterations based on previous.
            u = (2.0 / <double>kv) * <double>ks
            for j + 1 <= j <= max:
                t *= u
                if i-j > 0:
                    t /= <double>(i-j)
                t /= <double>j
                ret += t

        return ret

    #                            ^
    # Compute the second term of d, which can be factored as:
    #
    #        kv
    #      ----   i   i
    #   1  \     2  kv
    #  ---  >    ------
    #   4  /        i
    #      ----  i n
    #      i=1
    #
    # Reduce the chances of numerical overflow by computing the numerator and
    # denominator separately, trying to cancel common factors at each step.
    cdef double dist2(self, unsigned n, unsigned kv) except -1.0:
        cdef double ret, t, u
        cdef unsigned h, i

        ret = 0.0
        # Compute first iteration.
        i = 1
        t = 0.5
        t *= <double>kv
        t /= <double>n
        ret += t

        # Compute remaining iterations based on previous.
        u = 2.0 / <double>n * <double>kv
        for i + 1 <= i < kv:
            t *= u
            ret += t / <double>i

        return ret

    cdef double dist(self, char *a, char *b, unsigned len) except -1.0:
        cdef unsigned n, ks, kv

        self.stats(a, b, len, &n, &ks, &kv)

        if n == 0:
            return NAN

        return self.dist1(n, ks, kv) + self.dist2(n, kv)

cdef class LogDet:
    cdef Character char_
    cdef bint scoreGaps
    cdef int code2val[128]
    cdef double *A
    cdef unsigned n

    def __cinit__(self):
        self.A = NULL

    def __dealloc__(self):
        if self.A != NULL:
            free(self.A)
            self.A = NULL

    def __init__(self, type charType, bint scoreGaps):
        cdef str code
        cdef char *sCode
        cdef int iCode

        assert issubclass(charType, Character)
        self.char_ = charType.get()

        for code in self.char_.codes():
            sCode = <char *>code
            iCode = sCode[0]
            self.code2val[iCode] = self.char_.code2val(code)

        self.scoreGaps = scoreGaps

        self.n = self.char_.nstates()
        self.A = <double *>calloc(self.n * self.n, sizeof(double))
        if self.A == NULL:
            raise MemoryError("Matrix allocation failed")

    # Compute LogDet/paralinear distance for a vs. b.
    cdef double dist(self, char *a, char *b, unsigned seqlen) except -1.0:
        cdef unsigned j, aPop, bPop
        cdef int aC, bC, aI, bI
        cdef int ambiguous, aVal, bVal, aOff, bOff
        cdef double p, sum

        ambiguous = self.char_.any
        memset(self.A, 0, self.n * self.n * sizeof(double))

        for 0 <= j < seqlen:
            aC = a[j]
            bC = b[j]
            aVal = self.code2val[aC]
            aPop= pop(aVal)
            if aC == bC and aPop== 1:
                aOff = ffs(aVal) - 1
                # Identical non-ambiguous characters.
                self.A[aOff*self.n + aOff] += 1.0
            else:
                bVal = self.code2val[bC]
                bPop = pop(bVal)

                if aPop== 1 and bPop == 1:
                    aOff = ffs(aVal) - 1
                    bOff = ffs(bVal) - 1
                    self.A[aOff*self.n + bOff] += 1.0
                elif not self.scoreGaps and (aPop== 0 or bPop == 0):
                    # Ignore this character.
                    pass
                else:
                    # Ambiguous character(s).

                    # Translate gaps to mean "any character".
                    if aVal == 0:
                        aVal = ambiguous
                        aPop= pop(aVal)
                    if bVal == 0:
                        bVal = ambiguous
                        bPop = pop(bVal)

                    # Each possible state combination is assigned an equal
                    # probability.  All combinations sum to 1.0.
                    p = <double>1.0 / <double>(aPop * bPop)

                    # Iterate over every state combination, and add p to the
                    # combinations that are represented by the ambiguity.
                    for 0 <= aI < self.n:
                        for 0 <= bI < self.n:
                            if (aVal & (1 << aI)) and (bVal & (1 << bI)):
                                self.A[aI*self.n + bI] += p

        # If the matrix is empty, there are no data on which to base a distance.
        sum = cblas_dasum(self.n * self.n, self.A, 1)
        if sum == 0.0:
            return NAN

        return CxMatLogDet(self.n, self.A) / <double>self.n

cdef class Alignment:
    def __cinit__(self):
        self.rows = NULL
        self.cols = NULL

    def __dealloc__(self):
        if self.rows != NULL:
            free(self.rows)
            self.rows = NULL
        if self.cols != NULL:
            free(self.cols)
            self.cols = NULL

    def __init__(self, CTMatrix matrix=None, str pad=None,
      Taxa.Map taxaMap=None, int nchars=-1, type charType=Dna,
      bint rowMajor=True, bint colMajor=False):
        """
Create an alignment, optionally from a CTMatrix.  Unless a pad symbol is
specified, all sequences must be of equal length.

If no CTMatrix is specified, create an empty alignment for the taxa in taxaMap, with space for nchars characters.

By default, create only a row-major matrix, but support row-major and/or
column-major format.  The get{Row,Seq}() methods are only supported if rowMajor
is enabled.  Likewise, the get{Col,Char}() methods are only supported if
colMajor is enabled.
"""
        cdef int ntaxa, i
        cdef list taxa
        cdef Taxon taxon

        assert issubclass(charType, Character)

        if not rowMajor and not colMajor:
            raise ValueError("Enable at least one of {rowMajor,colMajor}")

        if matrix is not None:
            charType = matrix.charType
            taxaMap = Taxa.Map(matrix.taxaMap.taxaGet())
            ntaxa = taxaMap.ntaxa
            if ntaxa == 0:
                raise ValueError("Invalid alignment size: 0x0")

            taxa = taxaMap.taxaGet()
            nchars = len(matrix.dataGet(taxa[0]))
            if nchars == 0:
                raise ValueError("Invalid alignment size: %dx0" % ntaxa)
            for taxon in taxa:
                if len(matrix.dataGet(taxon)) != nchars:
                    if pad is None:
                        raise ValueError(
                          "Inconsistent sequence lengths (%d vs. %d)" %
                          (nchars, len(matrix.dataGet(taxon))))
                    if len(matrix.dataGet(taxon)) > nchars:
                        nchars = len(matrix.dataGet(taxon))
        else:
            if taxaMap is None:
                raise ValueError("Specify matrix or taxaMap")
            taxaMap = Taxa.Map(taxaMap.taxaGet())
            ntaxa = taxaMap.ntaxa
            if ntaxa == 0 or nchars <= 0:
                raise ValueError("Invalid alignment size: %dx%d" %
                  (ntaxa, nchars))

        self.charType = charType
        self.taxaMap = taxaMap
        self.ntaxa = ntaxa
        self.nchars = nchars
        if rowMajor:
            self.rows = self._allocMatrix(ntaxa, nchars, pad)
        if colMajor:
            self.cols = self._allocMatrix(ntaxa, nchars, pad)

        if matrix is not None:
            # Fill in matrix.
            i = 0
            for taxon in taxa:
                self.setSeq(i, 0, matrix.dataGet(taxon))
                i += 1

    cdef char *_allocMatrix(self, int ntaxa, int nchars, str pad) except NULL:
        cdef char *ret
        cdef char *padS

        ret = <char *>malloc(ntaxa * nchars)
        if ret == NULL:
            raise MemoryError("Matrix allocation failed")

        if pad is not None:
            padS = pad
            memset(ret, padS[0], ntaxa * nchars)

        return ret

    cdef void setRow(self, int row, int col, char *chars, unsigned len):
        cdef unsigned j

        assert row < self.ntaxa
        assert col + len <= self.nchars

        if self.rows != NULL:
            memcpy(&self.rows[(self.nchars * row) + col], chars, len)

        if self.cols != NULL:
            for 0 <= j < len:
                self.cols[(self.ntaxa * (col + j)) + row] = chars[col + j]

    cdef char *getRow(self, int row):
        assert self.rows != NULL
        assert row < self.ntaxa

        return &self.rows[row * self.nchars]

    cdef char *getCol(self, int col):
        assert self.cols != NULL
        assert col < self.nchars

        return &self.cols[col * self.ntaxa]

    cpdef setSeq(self, int seq, int col, str chars):
        cdef char *s = chars
        self.setRow(seq, col, s, len(chars))

    cpdef str getSeq(self, int seq):
        assert self.rows != NULL
        assert seq < self.ntaxa

        return PyString_FromStringAndSize(&self.rows[seq * self.nchars],
          self.nchars)

    cpdef str getChar(self, int char_):
        assert self.cols != NULL
        assert char_ < self.nchars

        return PyString_FromStringAndSize(&self.cols[char_ * self.ntaxa],
          self.taxa)

    cpdef str fastaPrint(self, file outFile=None):
        cdef list lines, taxa
        cdef Taxon taxon
        cdef str s, taxonData
        cdef int i, j

        if outFile is None:
            lines = []

        # Print.
        taxa = self.taxaMap.taxaGet()
        for 0 <= i < len(taxa):
            taxon = taxa[i]
            s = ">%s\n" % taxon.label.replace(' ', '_')
            if outFile is None:
                lines.append(s)
            else:
                outFile.write(s)

            # Break into lines of length 75.
            taxonData = self.getSeq(i)
            for 0 <= j < len(taxonData) by 75:
                if j + 75 < len(taxonData):
                    s = "%s\n" % taxonData[j:j+75]
                else:
                    s = "%s\n" % taxonData[j:]

                if outFile is None:
                    lines.append(s)
                else:
                    outFile.write(s)

        if outFile is None:
            if len(lines) > 0:
                # Chop off the last '\n'.
                lines[-1] = lines[-1][:-1]
            return "".join(lines)
        else:
            return None

    cpdef DistMatrix dists(self, bint avgAmbigs=True, bint scoreGaps=True):
        """
Calculate uncorrected pairwise distances.

If scoreGaps is enabled, gaps are treated as being ambiguous, rather than as
missing.

If avgAmbigs is enabled, ambiguity codes are handled by equal-weight averaging
of all the cases that the ambiguities could possibly resolve to.
"""
        cdef DistMatrix ret
        cdef PctIdent tab
        cdef int i, j
        cdef char *iRow, *jRow
        cdef unsigned nvalid
        cdef float fident, dist

        assert self.rows != NULL

        ret = DistMatrix(self.taxaMap)

        if self.ntaxa > 1:
            tab = PctIdent(self.charType, avgAmbigs, scoreGaps)
            for 0 <= i < self.ntaxa:
                iRow = self.getRow(i)
                for i + 1 <= j < self.ntaxa:
                    jRow = self.getRow(j)
                    tab.stats(iRow, jRow, self.nchars, &fident, &nvalid)
                    if nvalid > 0:
                        dist = 1.0 - (fident / <float>nvalid)
                    else:
                        dist = 1.0
                    ret.distanceSet(i, j, dist)

        return ret

    cpdef DistMatrix jukesDists(self, bint avgAmbigs=True,
      bint scoreGaps=False):
        """
Calculate pairwise distances, corrected for multiple hits using the
Jukes-Cantor method, and the computational methods described by:

  Tajima, F. (1993) Unbiased estimation of evolutionary distance between
  nucleotide sequences.  Mol. Biol. Evol. 10(3):677-688.

If avgAmbigs is enabled, ambiguity codes are handled by equal-weight averaging
of all the cases that the ambiguities could possibly resolve to.

If scoreGaps is enabled, gaps are treated as being ambiguous, rather than as
missing.
"""
        cdef DistMatrix ret
        cdef PctIdent tab
        cdef int nstates, i, j
        cdef char *iRow, *jRow
        cdef unsigned n, x
        cdef float fident, k, b, p, d, t

        assert self.rows != NULL

        ret = DistMatrix(self.taxaMap)

        if self.ntaxa > 1:
            tab = PctIdent(self.charType, avgAmbigs, scoreGaps)
            nstates = self.charType.get().nstates()
            b = <float>(nstates - 1) / <float>nstates

            for 0 <= i < self.ntaxa:
                iRow = self.getRow(i)
                for i + 1 <= j < self.ntaxa:
                    jRow = self.getRow(j)
                    tab.stats(iRow, jRow, self.nchars, &fident, &n)
                    if n > 0:
                        k = <float>n - fident
                        p = k / <float>n

                        d = 0.0
                        t = b
                        for 1 <= x <= <unsigned>k:
                            t *= ((k / <float>(x * n))) / b
                            d += t
                            if t < FLT_EPSILON:
                                # Terms in the summation monotonically decrease,
                                # so there's no point in continuing.
                                break
                    else:
                        d = NAN
                    ret.distanceSet(i, j, d)

        return ret

    cdef void _kimuraDistsDNA(self, DistMatrix m, bint scoreGaps):
        cdef K2p k2p
        cdef int a, b
        cdef char *aRow, *bRow
        cdef double d

        k2p = K2p(scoreGaps)
        for 0 <= a < self.ntaxa:
            aRow = self.getRow(a)
            for a + 1 <= b < self.ntaxa:
                bRow = self.getRow(b)
                d = k2p.dist(aRow, bRow, self.nchars)
                m.distanceSet(a, b, d)

    cdef void _kimuraDistsProtein(self, DistMatrix m, bint avgAmbigs,
      bint scoreGaps):
        cdef PctIdent tab
        cdef int a, b
        cdef char *aRow, *bRow
        cdef unsigned n, iud
        cdef float fident, ud, d
        cdef list dayhoff_pams = [ # PAMs, [75.0% .. 93.0%] in 0.1% increments.
          195, 196, 197, 198, 199, 200, 200, 201, 202, 203,
          204, 205, 206, 207, 208, 209, 209, 210, 211, 212,
          213, 214, 215, 216, 217, 218, 219, 220, 221, 222,
          223, 224, 226, 227, 228, 229, 230, 231, 232, 233,
          234, 236, 237, 238, 239, 240, 241, 243, 244, 245,
          246, 248, 249, 250, 252, 253, 254, 255, 257, 258,
          260, 261, 262, 264, 265, 267, 268, 270, 271, 273,
          274, 276, 277, 279, 281, 282, 284, 285, 287, 289,
          291, 292, 294, 296, 298, 299, 301, 303, 305, 307,
          309, 311, 313, 315, 317, 319, 321, 323, 325, 328,
          330, 332, 335, 337, 339, 342, 344, 347, 349, 352,
          354, 357, 360, 362, 365, 368, 371, 374, 377, 380,
          383, 386, 389, 393, 396, 399, 403, 407, 410, 414,
          418, 422, 426, 430, 434, 438, 442, 447, 451, 456,
          461, 466, 471, 476, 482, 487, 493, 498, 504, 511,
          517, 524, 531, 538, 545, 553, 560, 569, 577, 586,
          595, 605, 615, 626, 637, 649, 661, 675, 688, 703,
          719, 736, 754, 775, 796, 819, 845, 874, 907, 945,
          988]

        tab = PctIdent(self.charType, avgAmbigs, scoreGaps)
        for 0 <= a < self.ntaxa:
            aRow = self.getRow(a)
            for a + 1 <= b < self.ntaxa:
                bRow = self.getRow(b)
                tab.stats(aRow, bRow, self.nchars, &fident, &n)
                if n == 0:
                    m.distanceSet(a, b, NAN)
                    continue

                ud = 1.0 - (fident / <float>n)
                # Convert ud to iud/1000, in order to be able to quantize
                # uncorrected distance in tenths of a percent.
                iud = <unsigned>roundf(ud * 1000.0)
                if iud == 0:
                    d = 0.0
                elif iud < 750:
                    # Use Kimura's empirical formula.
                    d = -logf(1.0 - ud - ((ud / 5.0) * ud))
                elif iud <= 930:
                    # Correct distance using Dayhoff PAMs.
                    d = <float>(<int>dayhoff_pams[iud - 750]) / 100.0
                else:
                    d = 10.0 # Arbitrary excessively large value.

                m.distanceSet(a, b, d)

    cpdef DistMatrix kimuraDists(self, bint avgAmbigs=True,
      bint scoreGaps=False):
        """
Calculate corrected pairwise distances.

DNA:
  Correct for multiple hits using the Kimura two-parameter method, and the
  computational methods described by:

    Tajima, F. (1993) Unbiased estimation of evolutionary distance between
    nucleotide sequences.  Mol. Biol. Evol. 10(3):677-688.

Protein:
  Correct for multiple hits using the empirical formula presented as equation
  (4.8) in:

    Kimura, M. (1983) The neutral theory of molecular evolution.  p. 75,
    Cambridge University Press, Cambridge, England.

  Starting at ~75% divergence, this formula is inadequate, so use a method
  based on Dayhoff PAM matrices, as developed in CLUSTAL W.

    Dayhoff, M.O., R.M. Schwartz, B.C. Orcutt (1978) "A model of evolutionary
    change in proteins." In "Atlas of Protein Sequence and Structure, vol. 5,
    suppl. 3." M.O. Dayhoff (ed.), pp. 345-352, Natl.  Biomed. Res. Found.,
    Washington, DC.

    Thompson, J.D., D.G. Higgins, T.J. Gibson (1994) CLUSTAL W: improving the
    sensitivity of progressive multiple sequence alignment through sequence
    weighting, position specific gap penalties and weight matrix choice.
    Nucleic Acids Res. 22:4673-4680.

For DNA, ambiguity codes are handled by equal-weight averaging of all the cases
that the ambiguities could possibly resolve to.  Unlike for some of the other
distance correction methods (such as jukesDists()), it would be completely
non-sensical to treat DNA ambiguities as unique states.  For protein, if
avgAmbigs is enabled, ambiguity codes are handled by equal-weight averaging of
all the cases that the ambiguities could possibly resolve to.

If scoreGaps is enabled, gaps are treated as being ambiguous, rather than as
missing.  Ambiguity is interpreted to mean that all states are equally likely,
so enabling this option tends to substantially increase distance, due to
unlikely mutations being highly represented.
"""
        cdef DistMatrix ret

        assert self.rows != NULL

        ret = DistMatrix(self.taxaMap)

        if self.ntaxa > 1:
            if self.charType is Dna:
                self._kimuraDistsDNA(ret, scoreGaps)
            elif self.charType is Protein:
                self._kimuraDistsProtein(ret, avgAmbigs, scoreGaps)
            else:
                raise ValueError(
                  "Unsupported character type for Kimura distance correction")

        return ret

    cpdef DistMatrix logdetDists(self, bint scoreGaps=False):
        """
Calculate pairwise distances, corrected for unequal state frequencies using the
LogDet/paralinear method described by:

  Lake, J.A. (1994) Reconstructing evolutionary trees from DNA and protein
  sequences: Paralinear distances.  Proc. Natl. Acad. Sci. 91:1455-1459.

The distances are scaled such that they represent the mean number of
substitutions per site.

Note that the Lake (1994) method differs from the Lockhart et al. (1994) method
primarily in that it converges on additive distances for consistent data.

  Lockhart, P.J., M.A. Steel, M.D. Hendy, and D. Penny (1994) Recovering
  evolutionary trees under a more Realistic model of sequence evolution.  Mol.
  Biol. Evol. 11(4):605-612.

If scoreGaps is enabled, gaps are treated as being ambiguous, rather than as
missing.

Ambiguity codes are handled by equal-weight averaging of all the cases that the
ambiguities could possibly resolve to.  Unlike for some of the other distance
correction methods (such as jukesDists()), it would be completely non-sensical
to treat ambiguities as unique states.
"""
        cdef DistMatrix ret
        cdef LogDet logDet
        cdef int i, j
        cdef char *iRow, *jRow
        cdef double d

        assert self.rows != NULL

        ret = DistMatrix(self.taxaMap)

        if self.ntaxa > 1:
            logDet = LogDet(self.charType, scoreGaps)

            for 0 <= i < self.ntaxa:
                iRow = self.getRow(i)
                for i + 1 <= j < self.ntaxa:
                    jRow = self.getRow(j)
                    d = logDet.dist(iRow, jRow, self.nchars)
                    ret.distanceSet(i, j, d)

        return ret
