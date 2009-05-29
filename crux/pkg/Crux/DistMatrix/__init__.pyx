"""
    Pairwise distance matrices.
"""
import Crux.Exception

class Exception(Crux.Exception.Exception):
    pass

import exceptions

class SyntaxError(Exception, exceptions.SyntaxError):
    def __init__(self, line, str):
        self.line = line
        self.str = str

    def __str__(self):
        return "Line %d: %s" % (self.line, self.str)

class ValueError(Exception, exceptions.ValueError):
    def __init__(self, str):
        self._str = str

    def __str__(self):
        return self._str

class MemoryError(Exception, exceptions.MemoryError):
    def __init__(self, str):
        self._str = str

    def __str__(self):
        return self._str

from Crux.Taxa cimport Taxon
cimport Crux.Taxa as Taxa
cimport Crux.Tree as Tree

from libc cimport *

from CxDistMatrixLexer cimport *
cimport Crux.DistMatrix.Nj as Nj

import random
import re
import sys

# Forward declaration.
cdef class DistMatrix

#===============================================================================
#
# This parser implements a superset of the PHYLIP distance matrix format.  There
# are several rather byzantine limitations in the PHYLIP distance matrix format
# that are not worth faithfully enforcing.
#
# All tokens are separated by whitespace ([ \n\r\t]+).  The first token
# (integer) specifies the number of taxa.  The next token must be a taxon label,
# followed by zero or more distances.  Distances must be in decimal (4.2000),
# integer (42), or exponential format (4.2e-42, 4.2E42, etc.).  Taxon labels
# must not be interpretable as distances, but otherwise can be composed of any
# printable characters.
#
# Distance matrices may be specified as full matrices (these must be
# symmetrical):
#
#   5
#   Taxon_A 0.0 1.0 2.0 3.0 4.0
#   Taxon_B 1.0 0.0 1.5 2.5 3.5
#   Taxon_C 2.0 1.5 0.0 2.2 3.2
#   Taxon_D 3.0 2.5 2.2 0.0 3.1
#   Taxon_E 4.0 3.5 3.2 3.1 0.0
#
# or as upper triangle symmetric matrices:
#
#   5
#   Taxon_A     1.0 2.0 3.0 4.0
#   Taxon_B         1.5 2.5 3.5
#   Taxon_C             2.2 3.2
#   Taxon_D                 3.1
#   Taxon_E
#
# or as lower triangle symmetric matrices:
#
#   5
#   Taxon_A
#   Taxon_B 1.0
#   Taxon_C 2.0 1.5
#   Taxon_D 3.0 2.5 2.2
#   Taxon_E 4.0 3.5 3.2 3.1
#
#===============================================================================

cdef enum Format:
    FormatUninitialized,
    FormatUnknown,
    FormatFull,
    FormatUpper,
    FormatLower

cdef class Parser:
    cdef int line
    cdef DistMatrix matrix
    cdef Taxa.Map taxaMap
    cdef Format format
    cdef size_t i, j

    def __init__(self, DistMatrix matrix, Taxa.Map taxaMap):
        self.matrix = matrix
        self.taxaMap = taxaMap
        self.format = FormatUninitialized

    cdef void tokenDist(self, float dist, unsigned line) except *:
        cdef DistMatrix m
        cdef size_t i, j, k

        m = self.matrix
        i = self.i
        j = self.j

        if self.format == FormatUninitialized:
            raise SyntaxError(line, "Unexpected distance")
        elif self.format == FormatUnknown:
            if j < m.ntaxa:
                m.dists[nxy2i(m.ntaxa, i, j)] = dist
            else:
                # This is the last distance for the first row of a full matrix.
                # Up to now we have been working under the assumption that the
                # matrix would turn out to be in upper format, so now we need
                # to discard the diagonal and shift all the distances down one
                # element.
                self.format = FormatFull
                if m.dists[nxy2i(m.ntaxa, 0, 1)] != <float>0.0:
                    raise SyntaxError(line, "Non-zero distance on diagonal")
                # Shift the entire first row back one element.
                for 1 <= k < j - 1:
                    m.dists[nxy2i(m.ntaxa, 0, k)] = \
                      m.dists[nxy2i(m.ntaxa, 0, k + 1)]
                j -= 1
                m.dists[nxy2i(m.ntaxa, i, j)] = dist
        elif self.format == FormatFull:
            if j < i:
                if dist != m.dists[nxy2i(m.ntaxa, j, i)]:
                    raise SyntaxError(line,
                      "Non-symmetric distance for %r <--> %r" %
                      (self.taxaMap.taxonGet(j).label,
                      self.taxaMap.taxonGet(i).label))
            elif j == i:
                if dist != 0.0:
                    raise SyntaxError(line, "Non-zero distance on diagonal")
            elif j == m.ntaxa:
                raise SyntaxError(line, "Too many distances")
            else:
                m.dists[nxy2i(m.ntaxa, i, j)] = dist
        elif self.format == FormatUpper:
            if j == m.ntaxa:
                raise SyntaxError(line, "Too many distances")
            m.dists[nxy2i(m.ntaxa, i, j)] = dist
        elif self.format == FormatLower:
            if j == i:
                raise SyntaxError(line, "Too many distances")
            m.dists[nxy2i(m.ntaxa, i, j)] = dist
        else:
            assert False

        j += 1
        self.j = j

    cdef void nextRow(self, unsigned line) except *:
        cdef size_t i, j

        i = self.i
        j = self.j

        if self.format == FormatUninitialized:
            self.format = FormatUnknown
            i = 0
            j = 1
        elif self.format == FormatUnknown:
            if i == 0:
                if j == self.matrix.ntaxa:
                    self.format = FormatUpper
                    i = 1
                    j = 2
                elif j == 1:
                    self.format = FormatLower
                    i = 1
                    j = 0
                else:
                    raise SyntaxError(line, "Incorrect number of distances")
            else:
                assert False
        elif self.format == FormatFull:
            if j != self.matrix.ntaxa:
                raise SyntaxError(line, "Incorrect number of distances")
            i += 1
            j = 0
        elif self.format == FormatUpper:
            if j != self.matrix.ntaxa:
                raise SyntaxError(line, "Incorrect number of distances")
            i += 1
            j = i + 1
        elif self.format == FormatLower:
            if j != i:
                raise SyntaxError(line, "Incorrect number of distances")
            i += 1
            j = 0
        else:
            assert False

        self.i = i
        self.j = j

    cdef void tokenLabel(self, str label, unsigned line) except *:
        cdef Taxon taxon

        # Update row state.
        self.nextRow(line)

        # Create taxon mapping.
        taxon = Taxa.get(label)
        if self.taxaMap.indGet(taxon) != -1:
            raise SyntaxError(line, "Duplicate taxon label: %r" % label)
        if self.taxaMap.ntaxa == self.matrix.ntaxa:
            raise SyntaxError(line, "Excess taxon: %r" % label)
        self.taxaMap.map(taxon, self.taxaMap.ntaxa)

    cdef void parse(self, input, int line=1) except *:
        cdef yyscan_t scanner
        cdef CxtDistMatrixLexerExtra extra
        cdef int status, ntaxa, tokType

        assert type(input) in (file, str)

        extra.line = line
        extra.ioerror = 0
        if type(input) == file:
            extra.inputMode = CxeDistMatrixLexerInputModeFd
            extra.input.fd = input.fileno()
        elif type(input) == str:
            extra.inputMode = CxeDistMatrixLexerInputModeStr
            extra.input.s.s = input
            extra.input.s.len = len(input)
        else:
            assert False

        if CxDistMatrixLexer_lex_init_extra(&extra, &scanner) != 0:
            raise MemoryError("Error initializing lexer")

        try:
            # Get number of taxa.
            status = CxDistMatrixLexer_lex(scanner)
            if status > 0:
                if extra.tokType != CxeDistMatrixLexerTokTypeInt:
                    raise SyntaxError(extra.line, "Number of taxa expected")
                ntaxa = int(CxDistMatrixLexer_get_text(scanner))
                self.matrix._allocDists(ntaxa)
                self.matrix.ntaxa = ntaxa
                status = CxDistMatrixLexer_lex(scanner)

            # Get rows.
            while status > 0:
                tokType = extra.tokType
                if tokType in (CxeDistMatrixLexerTokTypeInt, \
                  CxeDistMatrixLexerTokTypeDist):
                    self.tokenDist(extra.tok.dist, extra.line)
                elif tokType == CxeDistMatrixLexerTokTypeLabel:
                    self.tokenLabel(extra.tok.label, extra.line)
                else:
                    assert False

                status = CxDistMatrixLexer_lex(scanner)
            if status == -1:
                raise SyntaxError(extra.line, "Invalid character")
            if extra.ioerror != 0:
                raise IOError("I/O error for %r" % input)
            self.nextRow(extra.line)
            if self.taxaMap.ntaxa != self.matrix.ntaxa:
                raise SyntaxError(extra.line, "Insufficient taxa")
        finally:
            CxDistMatrixLexer_lex_destroy(scanner);

#
#===============================================================================

cdef class DistMatrix:
    """
        Distance matrix.

        Construct a symmetric DistMatrix from one of the following inputs:

          file/str : Parse the input file/string as a distance matrix.

          Taxa.Map : Create an uninitialized distance matrix of the appropriate
                     size, given the number of taxa in the Taxa.Map.

          DistMatrix : Duplicate or sample from the input DistMatrix, depending
                       on the value of the 'sampleSize' parameter.
    """
    def __cinit__(self):
        self.dists = NULL

    def __dealloc__(self):
        if self.dists != NULL:
            free(self.dists)
            self.dists = NULL

    def __init__(self, input=None, int sampleSize=-1):
        if type(input) in (file, str):
            self._parse(input)
        elif type(input) is Taxa.Map:
            self.taxaMap = <Taxa.Map>input
            self._allocDists(input.ntaxa)
        elif type(input) is DistMatrix:
            self._dup(input, sampleSize)
        else:
            raise ValueError("File, string, Taxa.Map, or DistMatrix expected")

    cdef void _allocDists(self, size_t ntaxa) except *:
        assert self.dists == NULL
        assert ntaxa > 1
        self.dists = <float *>calloc(nxy2i(ntaxa, ntaxa - 2,
          ntaxa - 1) + 1, sizeof(float))
        if self.dists == NULL:
            raise MemoryError("Failed to allocate %d-taxon distance matrix" %
              <int>ntaxa)
        self.ntaxa = ntaxa

    cdef void _dup(self, DistMatrix other, int sampleSize) except *:
        cdef Taxa.Map taxaMap, otherMap
        cdef list sample
        cdef size_t i, j

        if sampleSize == -1:
            # Simple duplication.
            self.taxaMap = Taxa.Map(other.taxaMap.taxaGet())
            self._allocDists(other.ntaxa)
            memcpy(self.dists, other.dists, sizeof(float) * \
              (nxy2i(self.ntaxa, self.ntaxa - 2, self.ntaxa - 1) \
              + 1))
        else:
            if sampleSize < 2 or sampleSize > other.ntaxa:
                raise ValueError("sampleSize out of range (%d not in 2..%d)" \
                  % (sampleSize, other.ntaxa))

            # Create a sample of rows and construct a Taxa.Map for the new
            # DistMatrix.
            sample = random.sample(xrange(other.ntaxa), sampleSize)
            sample.sort()
            self.taxaMap = Taxa.Map([other.taxaMap.taxonGet(i) for i in sample])

            # Copy matrix data.
            self._allocDists(sampleSize)
            for 0 <= i < sampleSize:
                for i + 1 <= j < sampleSize:
                    self.dists[nxy2i(sampleSize, i, j)] = \
                      other.dists[nxy2i(other.ntaxa, sample[i], \
                      sample[j])]
        self.additive = other.additive

    cdef void _parse(self, input) except *:
        cdef Parser parser

        self.taxaMap = Taxa.Map()
        parser = Parser(self, self.taxaMap)
        parser.parse(input)

    cpdef float distanceGet(self, size_t x, size_t y):
        """
            Get distance at (x,y) in matrix.
        """
        assert x < self.ntaxa
        assert y < self.ntaxa

        if x == y:
            return 0.0

        return self.dists[nxy2i(self.ntaxa, x, y)]

    cpdef distanceSet(self, size_t x, size_t y, float distance):
        """
            Set distance at (x,y) in matrix.
        """
        assert x != y or distance == 0.0
        assert x < self.ntaxa
        assert y < self.ntaxa

        if x == y:
            return

        self.dists[nxy2i(self.ntaxa, x, y)] = distance

    cdef Tree _nj(self, bint random):
        cdef Nj.Nj nj

        nj = Nj.Nj()
        nj.prepare(self.dists, self.ntaxa, self.taxaMap)
        # The NJ code owns the matrix now.
        self.dists = NULL
        self.ntaxa = 0

        return nj.nj(random)

    cdef Tree _rnj(self, bint random, bint additive):
        cdef Tree ret
        cdef Nj.Rnj rnj

        rnj = Nj.Rnj()
        rnj.prepare(self.dists, self.ntaxa, self.taxaMap)
        # The RNJ code owns the matrix now.
        self.dists = NULL
        self.ntaxa = 0

        (ret, self.additive) = rnj.rnj(random, additive)
        return ret

    cdef void _rowsSwap(self, size_t a, size_t b):
        cdef size_t i
        cdef float distA, distB

        for 0 <= i < self.ntaxa:
            if i != a and i != b:
                distA = self.dists[nxy2i(self.ntaxa, i, a)]
                distB = self.dists[nxy2i(self.ntaxa, i, b)]

                self.dists[nxy2i(self.ntaxa, i, a)] = distB
                self.dists[nxy2i(self.ntaxa, i, b)] = distA

    cdef void _matrixShuffle(self, list order) except *:
        cdef list rowTab, curOrder
        cdef size_t i, a, b, frRow, t

        # Create a lookup table that maps original row to the current row of
        # the matrix (rowTab), as well as a table that represents the current
        # row order of the matrix (curOrder).  This is needed to keep track of
        # where rows end up as repeated row swaps are done.
        rowTab = [i for i in xrange(self.ntaxa)]
        curOrder = [i for i in xrange(self.ntaxa)]

        # Iteratively swap the correct row into row i.  The last row need not be
        # swapped with itself.
        for 0 <= i < self.ntaxa - 1:
            frRow = <size_t>order[i]
            a = i
            b = <size_t>rowTab[frRow]

            self._rowsSwap(a, b)

            # Update curOrder.
            t = curOrder[a]
            curOrder[a] = curOrder[b]
            curOrder[b] = t

            # Update rowTab.
            t = <size_t>rowTab[<size_t>curOrder[a]]
            rowTab[<size_t>curOrder[a]] = \
              <size_t>rowTab[<size_t>curOrder[b]]
            rowTab[<size_t>curOrder[b]] = t

    cpdef shuffle(self):
        """
            Randomly shuffle the order of rows/columns in the distance matrix,
            and make corresponding changes to the Taxa.Map.
        """
        cdef list order, taxa
        cdef Taxa.Map taxaMap
        cdef int i

        # Create a random shuffle order.
        order = [<size_t>i for i in xrange(self.ntaxa)]
        random.shuffle(order)

        # Shuffle the matrix.
        self._matrixShuffle(order)

        # Shuffle the Taxa.Map.
        taxaMap = self.taxaMap
        taxa = taxaMap.taxaGet()
        for i in xrange(self.ntaxa):
            taxaMap.map(<Taxon>taxa[<int>order[i]], i, True)

    cpdef Tree nj(self, bint joinRandom=False, bint destructive=False):
        """
            Construct a tree using the neighbor joining (NJ) algorithm.  If
            destructive=True, the matrix contents will be discarded, thus
            rendering the matrix useless for any further operations.
        """
        cdef DistMatrix m

        if (destructive):
            return self._nj(joinRandom)
        else:
            m = DistMatrix(self)
            return m._nj(joinRandom)

    cpdef Tree rnj(self, bint joinRandom=False, bint tryAdditive=True,
      bint destructive=False):
        """
            Construct a tree using the relaxed neighbor joining (RNJ)
            algorithm.  If destructive=True, the matrix contents will be
            discarded, thus rendering the matrix useless for any further
            operations.
        """
        cdef DistMatrix m

        if (destructive):
            return self._rnj(joinRandom, tryAdditive)
        else:
            m = DistMatrix(self)
            return m._rnj(joinRandom, tryAdditive)

    cpdef render(self, str format=None, str distFormat="%.7e",
      file outFile=None):
        """
            Print the matrix to a string in 'full', 'upper', or 'lower' format.
            If outFile is unspecified, print to stdout.
        """
        cdef size_t i, j
        cdef str s

        if format is None:
            format = 'lower'
        assert format in ('full', 'upper', 'lower')
        if outFile is None:
            outFile = sys.stdout

        distFormat = " " + distFormat

        outFile.write("%d\n" % self.ntaxa)
        if format == 'full':
            for 0 <= i < self.ntaxa:
                outFile.write("%-10s" % self.taxaMap.taxonGet(i).label)
                for 0 <= j < self.ntaxa:
                    s = distFormat % self.distanceGet(i, j)
                    outFile.write(s)
                outFile.write("\n")
        elif format == 'upper':
            for 0 <= i < self.ntaxa - 1:
                outFile.write("%-10s" % self.taxaMap.taxonGet(i).label)
                for 0 <= j < self.ntaxa:
                    s = distFormat % self.distanceGet(i, j)
                    if i < j:
                        outFile.write(s)
                    else:
                        outFile.write(" " * len(s))
                outFile.write("\n")
            outFile.write("%s\n" % self.taxaMap.taxonGet(self.ntaxa - 1).label)
        elif format == 'lower':
            outFile.write("%s\n" % self.taxaMap.taxonGet(0).label)
            for 1 <= i < self.ntaxa:
                outFile.write("%-10s" % self.taxaMap.taxonGet(i).label)
                for 0 <= j < i:
                    s = distFormat % self.distanceGet(i, j)
                    outFile.write(s)
                outFile.write("\n")
        else:
            assert False
