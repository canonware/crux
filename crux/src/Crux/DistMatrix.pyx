################################################################################
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
# printable non-whitespace characters.
#
# Distance matrices may be specified as full matrices (these need not be
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
################################################################################

import Crux.Exception

class Exception(Crux.Exception.Exception):
    pass

import exceptions

class ValueError(Exception, exceptions.ValueError):
    def __init__(self, str):
        self._str = str

    def __str__(self):
        return self._str

cimport Taxa
from Tree cimport Tree

import random

cdef class DistMatrix:
    # Construct a DistMatrix from one of the following inputs:
    #
    #   str : Parse the string as a distance matrix.
    #
    #   file : Parse the file as a distance matrix.
    #
    #   Taxa.Map : Create an uninitialized distance matrix of the appropriate
    #               size, given the number of taxa in the Taxa.Map.
    #
    #   DistMatrix : Duplicate or sample from the input DistMatrix, depending
    #                on the value of the 'sample' parameter.
    def __init__(self, input=None, symmetric=False, sample=None):
        # Validate input before calling the C_DistMatrix constructor.  I'm not
        # aware of a simple way to check the type of a Python-created class
        # in C code, which is why the check is done here.
        #
        # Also, make sure to pass in a Taxa.Map.
        if type(input) == file or type(input) == str:
            self._parse(input, Taxa.Map(), symmetric)
        elif type(input) == Taxa.Map:
            self._parse(None, input, symmetric)
        elif type(input) == DistMatrix:
            if sample == None:
                taxonMap = Taxa.Map(input.taxonMapGet().taxaGet())
                self._dup(input, taxonMap)
            else:
                if sample < 2 or sample > input.ntaxaGet():
                    raise ValueError("sample: Out of range (%d not in 2..%d)" \
                                      % (sample, input.ntaxaGet()))
                if symmetric:
                    raise ValueError("symmetric: Automatic for sampling")

                # Create a sample of rows.
                rows = random.sample(range(input.ntaxaGet()), sample)

                # Construct a TaxoMap for the new DistMatrix.
                inputMap = input.taxonMapGet()
                taxonMap = Taxa.Map()
                for row in rows:
                    taxonMap.map(inputMap.labelGet(row), taxonMap.ntaxaGet())

                self._sample(input, taxonMap, rows)
        else:
            raise ValueError(
                "input: File, string, Taxa.Map, or DistMatrix expected")

    cdef void _parse(self, input, Taxa.Map taxonMap, bint symmetric):
        pass # XXX

    # XXX Make a property.
    cdef int ntaxaGet(self):
        pass # XXX

    # XXX Make a property.
    cdef bint isSymmetric(self):
        pass # XXX

    cdef DistMatrix _dup(self, input, Taxa.Map taxonMap):
        pass # XXX

    cdef DistMatrix _sample(self, input, Taxa.Map taxonMap, list rows):
        pass # XXX

    # XXX Make a property.
    cdef Taxa.Map taxonMapGet(self):
        pass # XXX

    # XXX Make a property.
    cdef float distanceGet(self, int x, int y):
        pass # XXX
    cdef void distanceSet(self, int x, int y, float distance):
        pass # XXX

    cdef void _matrixShuffle(self, list order):
        pass # XXX

    cdef Tree _nj(self, bint random):
        rVal = Tree(taxonMap=self.taxonMapGet())
        pass # XXX

    cdef Tree _rnj(self, bint random, bint additive):
        rVal = Tree(taxonMap=self.taxonMapGet())
        pass # XXX

    cdef _render(self, format, distFormat, file file_=None):
        pass # XXX

    # Randomly shuffle the order of rows/columns in the distance matrix, and
    # make corresponding changes to the Taxa.Map.
    def shuffle(self):
        # Create a random shuffle order.
        order = random.sample(range(self.ntaxaGet()), self.ntaxaGet())

        # Shuffle the matrix.
        self._matrixShuffle(order)

        # Shuffle the Taxa.Map.
        taxonMap = self.taxonMapGet()
        labels = taxonMap.taxaGet()
        for i in xrange(self.ntaxaGet()):
            taxonMap.map(labels[order[i]], i, replace=True)

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
        if format == None:
            if self.isSymmetric():
                format = 'lower'
            else:
                format = 'full'

        # Make sure that distFormat isn't going to cause problems.
        distFormat % self.distanceGet(0, 1)

        if outFile == None:
            rVal = self._render(format, " " + distFormat)
        else:
            rVal = self._render(format, " " + distFormat, outFile)

        return rVal
