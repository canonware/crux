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

cimport Parsing
from Crux.Taxa cimport Taxon
cimport Crux.Taxa as Taxa
cimport Crux.Tree as Tree

cdef extern from "sys/types.h":
    ctypedef unsigned long size_t

cdef extern from "stdlib.h":
    cdef void *calloc(size_t nmemb, size_t size)
    cdef void free(void *ptr)

cdef extern from "string.h":
    cdef void *memcpy(void *dest, void *src, size_t n)

from CxDistMatrix cimport *
cimport Crux.DistMatrix.Nj as Nj

import random
import re
import sys

# Forward declarations.
cdef class pInt(Parsing.Precedence)
cdef class pNum(Parsing.Precedence)
cdef class pDists(Parsing.Precedence)
cdef class Token(Parsing.Token)
cdef class TokenInt(Token)
cdef class TokenNum(Token)
cdef class TokenLabel(Token)
cdef class Nonterm(Parsing.Nonterm)
cdef class Matrix(Nonterm)
cdef class Ntaxa(Nonterm)
cdef class Rows(Nonterm)
cdef class Row(Nonterm)
cdef class Label(Nonterm)
cdef class Dists(Nonterm)
cdef class Dist(Nonterm)
cdef class Parser(Parsing.Lr)
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
# Begin Precedence.
#

cdef class pInt(Parsing.Precedence):
    "%fail"
cdef class pNum(Parsing.Precedence):
    "%fail"
cdef class pDists(Parsing.Precedence):
    "%fail <pInt <pNum"

#
# End Precedence.
#===============================================================================
# Begin Token.
#

cdef class Token(Parsing.Token):
    cdef str raw
    cdef int line

    def __init__(self, Parser parser, str raw, int line):
        Parsing.Token.__init__(self, parser)

        self.raw = raw
        self.line = line

    def __repr__(self):
        return Parsing.Token.__repr__(self) + (" (%r)" % self.raw)

cdef class TokenInt(Token):
    "%token int [pInt]"
    cdef CxtDMSize i

    def __init__(self, Parser parser, str raw, int line):
        Token.__init__(self, parser, raw, line)

        self.i = int(raw)

cdef class TokenNum(Token):
    "%token num [pNum]"
    cdef CxtDMDist n

    def __init__(self, Parser parser, str raw, int line):
        Token.__init__(self, parser, raw, line)

        self.n = float(raw)

cdef class TokenLabel(Token):
    "%token label"
    cdef str label

    def __init__(self, Parser parser, str raw, int line):
        Token.__init__(self, parser, raw, line)

        self.label = raw

#
# End Token.
#===============================================================================
# Begin Nonterm.
#

cdef class Nonterm(Parsing.Nonterm):
    def __init__(self, Parser parser):
        Parsing.Nonterm.__init__(self, parser)

cdef class Matrix(Nonterm):
    "%start"
    cpdef reduce(self, Ntaxa Ntaxa, Rows Rows):
        "%reduce Ntaxa Rows"
        cdef Parser parser

        parser = <Parser>self.parser

        if parser.taxaMap.ntaxa != parser.matrix.ntaxa:
            raise SyntaxError(parser.line, "Insufficient taxa")

cdef class Ntaxa(Nonterm):
    "%nonterm"
    cpdef reduce(self, TokenInt int_):
        "%reduce int"
        cdef DistMatrix m

        m = (<Parser>self.parser).matrix
        m._allocDists(int_.i)
        m.ntaxa = int_.i

cdef class Rows(Nonterm):
    "%nonterm"
    cpdef reduceOne(self, Row row):
        "%reduce Row"
    cpdef reduceExtend(self, Rows rows, Row row):
        "%reduce Rows Row"

cdef class Row(Nonterm):
    "%nonterm"
    cpdef reduce(self, Label Label, Dists Dists):
        "%reduce Label Dists"
        cdef Parser parser
        cdef CxtDMSize i, j

        parser = <Parser>self.parser
        i = parser.i
        j = parser.j

        if parser.format == FormatUnknown:
            assert i == 0
            if j == parser.matrix.ntaxa:
                parser.format = FormatUpper
                i = 1
                j = 2
            elif j == 1:
                parser.format = FormatLower
                i = 1
                j = 0
            else:
                raise SyntaxError(Label.label.line,
                  "Incorrect number of distances")
        elif parser.format == FormatFull:
            if j != parser.matrix.ntaxa:
                raise SyntaxError(Label.label.line,
                  "Incorrect number of distances")
            i += 1
            j = 0
        elif parser.format == FormatUpper:
            if j != parser.matrix.ntaxa:
                raise SyntaxError(Label.label.line,
                  "Incorrect number of distances")
            i += 1
            j = i + 1
        elif parser.format == FormatLower:
            if j != i:
                raise SyntaxError(Label.label.line,
                  "Incorrect number of distances")
            i += 1
            j = 0
        else:
            assert False

        parser.i = i
        parser.j = j

cdef class Label(Nonterm):
    "%nonterm"
    cdef TokenLabel label

    cpdef reduce(self, TokenLabel label):
        "%reduce label"
        cdef Parser parser
        cdef Taxon taxon
        cdef Taxa.Map taxaMap

        self.label = label

        parser = <Parser>self.parser

        # Create taxon mapping.
        taxon = Taxa.get(label.label)
        taxaMap = parser.taxaMap
        if taxaMap.indGet(taxon) != -1:
            raise SyntaxError(label.line, "Duplicate taxon label: %r" %
              label.label)

        if taxaMap.ntaxa == parser.matrix.ntaxa:
            raise SyntaxError(label.line, "Excess taxon: %r" % label.label)

        taxaMap.map(taxon, taxaMap.ntaxa)

cdef class Dists(Nonterm):
    "%nonterm"
    cpdef reduceEmpty(self):
        "%reduce [pDists]"
    cpdef reduceOne(self, Dist Dist):
        "%reduce Dist"
    cpdef reduceExtend(self, Dists Dists, Dist Dist):
        "%reduce Dists Dist"

cdef class Dist(Nonterm):
    "%nonterm"
    cpdef reduceInt(self, TokenInt int_):
        "%reduce int"
        self._setDist(int_.i, int_.line)

    cpdef reduceNum(self, TokenNum num):
        "%reduce num"
        self._setDist(num.n, num.line)

    cdef void _setDist(self, CxtDMDist dist, int line) except *:
        cdef Parser parser
        cdef DistMatrix m
        cdef CxtDMSize i, j, k

        parser = <Parser>self.parser
        m = parser.matrix
        i = parser.i
        j = parser.j

        if parser.format == FormatUnknown:
            if j < m.ntaxa:
                m.dists[CxDistMatrixNxy2i(m.ntaxa, i, j)] = dist
            else:
                # This is the last distance for the first row of a full matrix.
                # Up to now we have been working under the assumption that the
                # matrix would turn out to be in upper format, so now we need
                # to discard the diagonal and shift all the distances down one
                # element.
                parser.format = FormatFull
                if m.dists[CxDistMatrixNxy2i(m.ntaxa, 0, 1)] != 0.0:
                    raise SyntaxError(line, "Non-zero distance on diagonal")
                # Shift the entire first row back one element.
                for 1 <= k < j - 1:
                    m.dists[CxDistMatrixNxy2i(m.ntaxa, 0, k)] = \
                      m.dists[CxDistMatrixNxy2i(m.ntaxa, 0, k + 1)]
                j -= 1
                m.dists[CxDistMatrixNxy2i(m.ntaxa, i, j)] = dist
        elif parser.format == FormatFull:
            if j < i:
                if dist != m.dists[CxDistMatrixNxy2i(m.ntaxa, j, i)]:
                    raise SyntaxError(line,
                      "Non-symmetric distance for %r <--> %r" %
                      (parser.taxaMap.taxonGet(j).label,
                      parser.taxaMap.taxonGet(i).label))
            elif j == i:
                if dist != 0.0:
                    raise SyntaxError(line, "Non-zero distance on diagonal")
            elif j == m.ntaxa:
                raise SyntaxError(line, "Too many distances")
            else:
                m.dists[CxDistMatrixNxy2i(m.ntaxa, i, j)] = dist
        elif parser.format == FormatUpper:
            if j == m.ntaxa:
                raise SyntaxError(line, "Too many distances")
            m.dists[CxDistMatrixNxy2i(m.ntaxa, i, j)] = dist
        elif parser.format == FormatLower:
            if j == i:
                raise SyntaxError(line, "Too many distances")
            m.dists[CxDistMatrixNxy2i(m.ntaxa, i, j)] = dist

        j += 1
        parser.j = j

#
# End Nonterm.
#===============================================================================

# Regex used to recognize tokens.
cdef _re

cdef Parsing.Spec _spec

cdef enum Format:
    FormatUnknown,
    FormatFull,
    FormatUpper,
    FormatLower

cdef class Parser(Parsing.Lr):
    cdef int line
    cdef DistMatrix matrix
    cdef Taxa.Map taxaMap
    cdef Format format
    cdef CxtDMSize i, j

    def __init__(self, DistMatrix matrix, Taxa.Map taxaMap):
        global _re, _spec

        if _spec is None:
            _spec = self._initSpec()
        Parsing.Lr.__init__(self, _spec)

        if _re is None:
            _re = self._initRe()

        self.matrix = matrix
        self.taxaMap = taxaMap
        self.format = FormatUnknown
        self.i = 0
        self.j = 1

    cdef Parsing.Spec _initSpec(self):
        return Parsing.Spec([sys.modules[__name__]],
          pickleFile="%s/Crux/parsers/DistMatrix.pickle" % Crux.Config.datadir,
          verbose=(False if (not __debug__ or Crux.Config.quiet) else True),
          skinny=(False if __debug__ else True),
          logFile="%s/Crux/parsers/DistMatrix.log" % Crux.Config.datadir)

    cdef _initRe(self):
        return re.compile(r"""
    ([1-9][0-9]*(?=\s))                          # integer (#taxa or dist)
  | ([-+]?
     [0-9]+(?:[.][0-9]+)?
     (?:[eE][-+]?[0-9]+)?
     (?=\s|$))                                   # distance
  | (^[A-Za-z_][A-Za-z0-9_\-.]*
     (?:[ \t\r\f\v]+[A-Za-z_][A-Za-z0-9_\-.]*)*) # label
  | ([ \t\r\f\v]+)                               # whitespace
  | ([\n])                                       # whitespace (newline)
""", re.M | re.X)

    cpdef parse(self, lines, int line=1, bint verbose=False):
        cdef str l
        cdef int pos

        assert getattr3(lines, '__iter__', None) is not None

        self.verbose = verbose
        self.line = line

        for l in lines:
            pos = 0
            while pos < len(l):
                m = _re.match(l, pos)
                if m is None:
                    raise SyntaxError(self.line, "Invalid token")
                idx = m.lastindex
                start = m.start(idx)
                end = m.end(idx)
                if idx == 1:   # integer (#taxa or distance)
                    self.token(TokenInt(self, l[start:end], self.line))
                elif idx == 2: # distance
                    self.token(TokenNum(self, l[start:end], self.line))
                elif idx == 3: # label
                    self.token(TokenLabel(self, l[start:end], self.line))
                elif idx == 4: # whitespace
                    pass
                elif idx == 5: # whitespace (newline)
                    self.line += 1
                else:
                    assert False
                pos = end

        self.eoi()

#
#===============================================================================

cdef class DistMatrix:
    def __cinit__(self):
        self.dists = NULL

    def __dealloc__(self):
        if self.dists != NULL:
            free(self.dists)
            self.dists = NULL

    # Construct a symmetric DistMatrix from one of the following inputs:
    #
    #   lines : Parse the iterable lines as a distance matrix.
    #
    #   Taxa.Map : Create an uninitialized distance matrix of the appropriate
    #              size, given the number of taxa in the Taxa.Map.
    #
    #   DistMatrix : Duplicate or sample from the input DistMatrix, depending
    #                on the value of the 'sampleSize' parameter.
    def __init__(self, input=None, int sampleSize=-1):
        if getattr3(input, '__iter__', None) is not None:
            self._parse(input)
        elif type(input) is Taxa.Map:
            self.taxaMap = <Taxa.Map>input
            self._allocDists(input.ntaxa)
        elif type(input) is DistMatrix:
            self._dup(input, sampleSize)
        else:
            raise ValueError("Iterable lines, Taxa.Map, or DistMatrix expected")

    cdef void _allocDists(self, CxtDMSize ntaxa) except *:
        assert self.dists == NULL
        self.dists = <CxtDMDist *>calloc(CxDistMatrixNxy2i(ntaxa, ntaxa - 2,
          ntaxa - 1) + 1, sizeof(CxtDMDist))
        if self.dists == NULL:
            raise MemoryError("Failed to allocate %d-taxon distance matrix" %
              <int>ntaxa)
        self.ntaxa = ntaxa

    cdef void _dup(self, DistMatrix other, int sampleSize) except *:
        cdef Taxa.Map taxaMap, otherMap
        cdef list sample
        cdef CxtDMSize i, j

        if sampleSize == -1:
            # Simple duplication.
            self.taxaMap = Taxa.Map(other.taxaMap.taxaGet())
            self._allocDists(other.ntaxa)
            memcpy(self.dists, other.dists, sizeof(CxtDMDist) * \
              (CxDistMatrixNxy2i(self.ntaxa, self.ntaxa - 2, self.ntaxa - 1) \
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
                    self.dists[CxDistMatrixNxy2i(sampleSize, i, j)] = \
                      other.dists[CxDistMatrixNxy2i(other.ntaxa, sample[i], \
                      sample[j])]

    cdef void _parse(self, input) except *:
        cdef Parser parser

        self.taxaMap = Taxa.Map()
        parser = Parser(self, self.taxaMap)
        parser.parse(input)

    cdef CxtDMDist distanceGet(self, CxtDMSize x, CxtDMSize y):
        if x == y:
            return 0.0

        return self.dists[CxDistMatrixNxy2i(self.ntaxa, x, y)]

    cdef void distanceSet(self, CxtDMSize x, CxtDMSize y, CxtDMDist distance):
        assert x != y or distance == 0.0
        if x == y:
            return

        self.dists[CxDistMatrixNxy2i(self.ntaxa, x, y)] = distance

    cdef Tree _nj(self, bint random):
        """
Construct a tree using the neighbor joining (NJ) algorithm.  This operation
discards the matrix contents, so maintain a duplicate matrix if needed.
"""
        cdef Nj.Nj nj

        nj = Nj.Nj()
        nj.prepare(self.dists, self.ntaxa, self.taxaMap)
        # The NJ code owns the matrix now.
        self.dists = NULL
        self.ntaxa = 0

        return nj.nj(random)

    cdef Tree _rnj(self, bint random, bint additive):
        """
Construct a tree using the relaxed neighbor joining (RNJ) algorithm.  This
operation discards the matrix contents, so maintain a duplicate matrix if
needed.
"""
        cdef Nj.Rnj rnj

        rnj = Nj.Rnj()
        rnj.prepare(self.dists, self.ntaxa, self.taxaMap)
        # The RNJ code owns the matrix now.
        self.dists = NULL
        self.ntaxa = 0

        return rnj.rnj(random, additive)

    cdef void _rowsSwap(self, CxtDMSize a, CxtDMSize b):
        cdef CxtDMSize i
        cdef CxtDMDist distA, distB

        for 0 <= i < self.ntaxa:
            if i != a and i != b:
                distA = self.dists[CxDistMatrixNxy2i(self.ntaxa, i, a)]
                distB = self.dists[CxDistMatrixNxy2i(self.ntaxa, i, b)]

                self.dists[CxDistMatrixNxy2i(self.ntaxa, i, a)] = distB
                self.dists[CxDistMatrixNxy2i(self.ntaxa, i, b)] = distA

    cdef void _matrixShuffle(self, list order):
        cdef list rowTab, curOrder
        cdef CxtDMSize i, a, b, frRow, t

        # Create a lookup table that maps original row to the current row of
        # the matrix (rowTab), as well as a table that represents the current
        # row order of the matrix (curOrder).  This is needed to keep track of
        # where rows end up as repeated row swaps are done.
        rowTab = range(self.ntaxa)
        curOrder = range(self.ntaxa)

        # Iteratively swap the correct row into row i.  The last row need not be
        # swapped with itself.
        for 0 <= i < self.ntaxa - 1:
            frRow = order[i]
            a = i
            b = rowTab[frRow]

            self._rowsSwap(a, b)

            # Update curOrder.
            t = curOrder[a]
            curOrder[a] = curOrder[b]
            curOrder[b] = t

            # Update rowTab.
            t = rowTab[curOrder[a]]
            rowTab[curOrder[a]] = rowTab[curOrder[b]]
            rowTab[curOrder[b]] = t

    # Randomly shuffle the order of rows/columns in the distance matrix, and
    # make corresponding changes to the Taxa.Map.
    cpdef shuffle(self):
        cdef list order, taxa
        cdef Taxa.Map taxaMap
        cdef int i

        # Create a random shuffle order.
        order = random.sample(xrange(self.ntaxa), self.ntaxa)

        # Shuffle the matrix.
        self._matrixShuffle(order)

        # Shuffle the Taxa.Map.
        taxaMap = self.taxaMap
        taxa = taxaMap.taxaGet()
        for i in xrange(self.ntaxa):
            taxaMap.map(<Taxon>taxa[<int>order[i]], i, True)

    # Construct a neighbor joining (NJ) tree from the distance matrix.
    cpdef Tree nj(self, bint joinRandom=False, bint destructive=False):
        cdef DistMatrix m

        if (destructive):
            return self._nj(joinRandom)
        else:
            m = DistMatrix(self)
            return m._nj(joinRandom)

    # Construct a relaxed neighbor joining (RNJ) tree from the distance matrix.
    cpdef Tree rnj(self, bint joinRandom=False, bint tryAdditive=True,
      bint destructive=False):
        cdef DistMatrix m

        if (destructive):
            return self._rnj(joinRandom, tryAdditive)
        else:
            m = DistMatrix(self)
            return m._rnj(joinRandom, tryAdditive)

    # Print the matrix to a string in 'full', 'upper', or 'lower' format.
    cpdef render(self, str format=None, str distFormat="%.5e",
      file outFile=None):
        cdef CxtDMSize i, j
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
