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
from Crux.Character cimport Character, Dna
from Crux.DistMatrix cimport DistMatrix

from libc cimport *
from libm cimport *
from atlas cimport *

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

# Compute the number of 1 bits in x.
cdef inline unsigned pop(unsigned x):
    assert sizeof(unsigned) == 4

    x = x - ((x >> 1) & 0x55555555)
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333)
    x = (x + (x >> 4)) & 0x0f0f0f0f
    x = x + (x >> 8)
    x = x + (x >> 16)
    return x & 0x0000003f

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
      float *fident, unsigned *nvalid):
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

    # Compute determinant of the pairwise frequency matrix for a vs. b.
    cdef double det(self, char *a, char *b, unsigned seqlen):
        cdef unsigned j, popA, popB
        cdef int aC, bC, aI, bI
        cdef int ambiguous, aVal, bVal, aOff, bOff
        cdef double sum, p

        ambiguous = self.char_.any
        memset(self.A, 0, self.n * self.n * sizeof(double))

        for 0 <= j < seqlen:
            aC = a[j]
            bC = b[j]
            aVal = self.code2val[aC]
            popA = pop(aVal)
            if aC == bC and popA == 1:
                aOff = ffs(aVal) - 1
                # Identical non-ambiguous characters.
                self.A[aOff*self.n + aOff] += 1.0
            else:
                bVal = self.code2val[bC]
                popB = pop(bVal)

                if popA == 1 and popB == 1:
                    aOff = ffs(aVal) - 1
                    bOff = ffs(bVal) - 1
                    self.A[aOff*self.n + bOff] += 1.0
                elif not self.scoreGaps and (popA == 0 or popB == 0):
                    # Ignore this character.
                    pass
                else:
                    # Ambiguous character(s).

                    # Translate gaps to mean "any character".
                    if aVal == 0:
                        aVal = ambiguous
                        popA = pop(aVal)
                    if bVal == 0:
                        bVal = ambiguous
                        popB = pop(bVal)

                    # Each possible state combination is assigned an equal
                    # probability.  All combinations sum to 1.0.
                    p = <double>1.0 / <double>(popA * popB)

                    # Iterate over every state combination, and add p to the
                    # combinations that are represented by the ambiguity.
                    for 0 <= aI < self.n:
                        for 0 <= bI < self.n:
                            if (aVal & (1 << aI)) and (bVal & (1 << bI)):
                                self.A[aI*self.n + bI] += p

        # Normalize matrix, such that the elements sum to 1.0.
        sum = cblas_dasum(self.n * self.n, self.A, 1)
        if sum > 0.0:
            cblas_dscal(self.n * self.n, <double>1.0 / sum, self.A, 1)

        return CxMatDdet(self.n, self.A)

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

    cpdef DistMatrix jukesDists(self, bint avgAmbigs=True, bint scoreGaps=True):
        """
Calculate pairwise distances, corrected for multiple hits using the
Jukes-Cantor method, and the computational methods described by:

  Tajima, F. (1993) Unbiased estimation of evolutionary distance between
  nucleotide sequences.  Mol. Biol. Evol. 10(3):677-688.

If scoreGaps is enabled, gaps are treated as being ambiguous, rather than as
missing.

If avgAmbigs is enabled, ambiguity codes are handled by equal-weight averaging
of all the cases that the ambiguities could possibly resolve to.
"""
        cdef DistMatrix ret
        cdef PctIdent tab
        cdef int nstates, i, j
        cdef char *iRow, *jRow
        cdef unsigned n, x
        cdef float NaN, fident, k, b, p, d, t

        assert self.rows != NULL

        ret = DistMatrix(self.taxaMap)

        if self.ntaxa > 1:
            tab = PctIdent(self.charType, avgAmbigs, scoreGaps)
            NaN = 1e10000 / 1e10000
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
                        d = NaN
                    ret.distanceSet(i, j, d)

        return ret

    cpdef DistMatrix logdetDists(self, bint scoreGaps=True):
        """
Calculate pairwise distances, corrected for unequal state frequencies using the
LogDet method described by:

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
        cdef int nstates, i, j
        cdef char *iRow, *jRow
        cdef unsigned n, x
        cdef float NaN

        assert self.rows != NULL


        ret = DistMatrix(self.taxaMap)
        logDet = LogDet(self.charType, scoreGaps)

        if self.ntaxa > 1:
            NaN = 1e10000 / 1e10000

            for 0 <= i < self.ntaxa:
                iRow = self.getRow(i)
                for i + 1 <= j < self.ntaxa:
                    jRow = self.getRow(j)
                    det = logDet.det(iRow, jRow, self.nchars)
                    if det > 0.0:
                        d = -log(det)
                    elif det == 0.0:
                        d = NaN
                    ret.distanceSet(i, j, d)

        return ret
