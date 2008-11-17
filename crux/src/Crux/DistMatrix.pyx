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

cimport Parsing
cimport Taxa
from Tree cimport Tree

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

    def __init__(self, Parser parser, str raw):
        Parsing.Token.__init__(self, parser)

        self.raw = raw

    def __repr__(self):
        return Parsing.Token.__repr__(self) + (" (%r)" % self.raw)

cdef class TokenInt(Token):
    "%token int [pInt]"
    cdef int i

    def __init__(self, Parser parser, str raw):
        Token.__init__(self, parser, raw)

        self.i = int(raw)

cdef class TokenNum(Token):
    "%token num [pNum]"
    cdef float n

    def __init__(self, Parser parser, str raw):
        Token.__init__(self, parser, raw)

        self.n = float(raw)

cdef class TokenLabel(Token):
    "%token label"
    cdef str label

    def __init__(self, Parser parser, str raw):
        Token.__init__(self, parser, raw)

        self.label = raw

#
# End Token.
#===============================================================================
# Begin Nonterm.
#

cdef class Nonterm(Parsing.Nonterm):
    def __init__(self, Parsing.Lr parser):
        Parsing.Nonterm.__init__(self, parser)

cdef class Matrix(Nonterm):
    "%start"
    cpdef reduce(self, Ntaxa Ntaxa, Rows Rows):
        "%reduce Ntaxa Rows"

cdef class Ntaxa(Nonterm):
    "%nonterm"
    cpdef reduce(self, TokenInt int_):
        "%reduce int"

cdef class Rows(Nonterm):
    "%nonterm"
    cpdef reduceOne(self, Row row):
        "%reduce Row"
    cpdef reduceExtend(self, Rows rows, Row row):
        "%reduce Rows Row"

cdef class Row(Nonterm):
    "%nonterm"
    cpdef reduce(self, TokenLabel label, Dists Dists):
        "%reduce label Dists"

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
    cpdef reduceNum(self, TokenNum num):
        "%reduce num"

#
# End Nonterm.
#===============================================================================

# Regex used to recognize tokens.
cdef _re

cdef Parsing.Spec _spec

cdef class Parser(Parsing.Lr):
    cdef DistMatrix matrix
    cdef Taxa.Map taxaMap

    def __init__(self, DistMatrix matrix, Taxa.Map taxaMap):
        global _re, _spec

        if _spec is None:
            _spec = self._initSpec()
        Parsing.Lr.__init__(self, _spec)

        if _re is None:
            _re = self._initRe()

        self.matrix = matrix
        self.taxaMap = taxaMap

    cdef Parsing.Spec _initSpec(self):
        return Parsing.Spec([sys.modules[__name__]],
          pickleFile="%s/share/Crux-%s/DistMatrix.pickle" %
          (Crux.Config.prefix, Crux.Config.version),
          verbose=(False if (not __debug__ or Crux.Config.quiet) else True),
          skinny=(False if __debug__ else True),
          logFile="%s/share/Crux-%s/DistMatrix.log" %
          (Crux.Config.prefix, Crux.Config.version))

    cdef _initRe(self):
        return re.compile(r"""
    ([1-9][0-9]*(?=\s))                       # integer (#taxa or distance)
  | ([-+]?
     [0-9]+(?:[.][0-9]+)?
     (?:[eE][-+]?[0-9]+)?
     (?=\s|$))                                # distance
  | (^[A-Za-z_][A-Za-z0-9_]*
     (?:[ \t\r\f\v]+[A-Za-z_][A-Za-z0-9_]*)*) # label
  | ([ \t\r\f\v]+)                            # whitespace
  | ([\n])                                    # whitespace (newline)
""", re.M | re.X)

    cpdef parse(self, lines, int line=1, bint verbose=False):
        cdef str l
        cdef int pos

        assert getattr3(lines, '__iter__', None) is not None

        self.verbose = verbose

        for l in lines:
            pos = 0
            while pos < len(l):
                m = _re.match(l, pos)
                if m is None:
                    raise SyntaxError(line, "Invalid token")
                idx = m.lastindex
                start = m.start(idx)
                end = m.end(idx)
                if idx == 1:   # integer (#taxa or distance)
                    self.token(TokenInt(self, l[start:end]))
                elif idx == 2: # distance
                    self.token(TokenNum(self, l[start:end]))
                elif idx == 3: # label
                    self.token(TokenLabel(self, l[start:end]))
                elif idx == 4: # whitespace
                    pass
                elif idx == 5: # whitespace (newline)
                    line += 1
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
            pass # XXX

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
            # XXX Create empty matrix.
            pass
        elif type(input) is DistMatrix:
            self._dup(input, sampleSize)
        else:
            raise ValueError("Iterable lines, Taxa.Map, or DistMatrix expected")

    cdef void _dup(self, DistMatrix other, int sampleSize) except *:
        cdef Taxa.Map taxaMap, otherMap
        cdef list sample

        if sampleSize == -1:
            taxaMap = Taxa.Map(other.taxaMap.taxaGet())
            # XXX Duplicate.
        else:
            if sampleSize < 2 or sampleSize > other.ntaxa:
                raise ValueError("sampleSize out of range (%d not in 2..%d)" \
                  % (sampleSize, other.ntaxa))

            # Create a sample of rows and construct a Taxa.Map for the new
            # DistMatrix.
            sample = random.sample(other.taxaGet(), sampleSize)
            sample.sort()
            taxaMap = Taxa.Map(sample)

            # Copy matrix data.
            # XXX

    cdef void _parse(self, input) except *:
        cdef Parser parser

        self.taxaMap = Taxa.Map()
        parser = Parser(self, self.taxaMap)
        parser.parse(input)

    property ntaxa:
        def __get__(self):
            pass # XXX

    property taxaMap:
        def __get__(self):
            pass # XXX

    cdef float distanceGet(self, int x, int y):
        pass # XXX
    cdef void distanceSet(self, int x, int y, float distance):
        pass # XXX

    cdef void _matrixShuffle(self, list order):
        pass # XXX

    cdef Tree _nj(self, bint random):
        rVal = Tree(taxaMap=self.taxaMapGet())
        pass # XXX

    cdef Tree _rnj(self, bint random, bint additive):
        rVal = Tree(taxaMap=self.taxaMapGet())
        pass # XXX

    cdef _render(self, format, distFormat, file file_=None):
        pass # XXX

    # Randomly shuffle the order of rows/columns in the distance matrix, and
    # make corresponding changes to the Taxa.Map.
    def shuffle(self):
        # Create a random shuffle order.
        order = random.sample(range(self.ntaxa), self.ntaxa)

        # Shuffle the matrix.
        self._matrixShuffle(order)

        # Shuffle the Taxa.Map.
        taxaMap = self.taxaMapGet()
        labels = taxaMap.taxaGet()
        for i in xrange(self.ntaxa):
            taxaMap.map(labels[order[i]], i, replace=True)

    # Construct a neighbor joining (NJ) tree from the distance matrix.
    def nj(self, joinRandom=False, destructive=False):
        if (destructive):
            return self._nj(joinRandom)
        else:
            return DistMatrix(self)._nj(joinRandom)

    # Construct a relaxed neighbor joining (RNJ) tree from the distance matrix.
    def rnj(self, joinRandom=False, tryAdditive=True, destructive=False):
        if (destructive):
            return self._rnj(joinRandom, tryAdditive)
        else:
            return DistMatrix(self)._rnj(joinRandom, tryAdditive)

    # Print the matrix to a string in 'full', 'upper', or 'lower' format.
    def render(self, format=None, distFormat="%.5e", outFile=None):
        if format is None:
            format = 'lower'

        # Make sure that distFormat isn't going to cause problems.
        distFormat % self.distanceGet(0, 1)

        if outFile == None:
            rVal = self._render(format, " " + distFormat)
        else:
            rVal = self._render(format, " " + distFormat, outFile)

        return rVal
