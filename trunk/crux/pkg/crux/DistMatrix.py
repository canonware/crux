################################################################################
#
# <Copyright = jasone>
# <License>
#
################################################################################
#
# Version: Crux <Version = crux>
#
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
# Note: The parser reads the input one line at a time.  In order to avoid
#       silently losing trailing text, the last line of the distance matrix
#       must either be the end of input (end of file or end of string), or there
#       must be no non-whitespace characters before the newline at the end of
#       the last line of the distance matrix.
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

from C_DistMatrix import *
import TaxonMap
import crux

import random

class DistMatrix(C_DistMatrix):
    # Construct a DistMatrix from one of the following inputs:
    #
    #   str : Parse the string as a distance matrix.
    #
    #   file : Parse the file as a distance matrix.
    #
    #   TaxonMap : Create an uninitialized distance matrix of the appropriate
    #              size, given the number of taxa in the TaxonMap.
    #
    #   DistMatrix : Duplicate or sample from the input DistMatrix, depending
    #                on the value of the 'sample' parameter.
    def __init__(self, input=None, symmetric=False, sample=None):
        # Validate input before calling the C_DistMatrix constructor.  I'm not
        # aware of a simple way to check the type of a Python-created class
        # in C code, which is why the check is done here.
        #
        # Also, make sure to pass in a TaxonMap.
        if type(input) == file or type(input) == str:
            C_DistMatrix._parse(self, input, TaxonMap.TaxonMap(), symmetric)
        elif type(input) == TaxonMap.TaxonMap:
            C_DistMatrix._parse(self, None, input, symmetric)
        elif type(input) == DistMatrix:
            if sample == None:
                map = TaxonMap.TaxonMap(input.taxonMapGet().taxaGet())
                C_DistMatrix._dup(self, input, map)
            else:
                if sample < 2 or sample > input.ntaxaGet():
                    raise crux.DistMatrix\
                          .ValueError("sample: Out of range (%d not in 2..%d)" \
                                      % (sample, input.ntaxaGet()))
                if symmetric:
                    raise crux.DistMatrix\
                          .ValueError("symmetric: Automatic for sampling")

                # Create a sample of rows.
                rows = random.sample(range(input.ntaxaGet()), sample)

                # Construct a TaxonMap for the new DistMatrix.
                inputMap = input.taxonMapGet()
                map = TaxonMap.TaxonMap()
                for row in rows:
                    map.map(inputMap.labelGet(row), map.ntaxaGet())

                C_DistMatrix._sample(self, input, map, rows)
        else:
            raise crux.DistMatrix\
                  .ValueError(
                "input: File, string, TaxonMap, or DistMatrix expected")

    # Randomly shuffle the order of rows/columns in the distance matrix, and
    # make corresponding changes to the TaxonMap.
    def shuffle(self):
        # Create a random shuffle order.
        order = random.sample(range(self.ntaxaGet()), self.ntaxaGet())

        # Shuffle the matrix.
        self._matrixShuffle(order)

        # Shuffle the TaxonMap.
        map = self.taxonMapGet()
        labels = map.taxaGet()
        for i in forints(self.ntaxaGet()):
            map.map(labels[order[i]], i, replace=True)

    # Print the matrix to a string in 'full', 'upper', or 'lower' format.
    def render(self, format='full', distFormat="%.5e", outFile=None):
        if outFile == None:
            retval = self._stringRender(format, distFormat)
        else:
            # Make sure that distFormat isn't going to cause problems.
            distFormat % self.distanceGet(0, 1)

            retval = self._fileRender(format, " " + distFormat, outFile)

        return retval

    def _stringRender(self, format, distFormat):
        retval = "%d\n" % self.taxonMapGet().ntaxaGet()
        if format == 'full':
            for x in forints(self.taxonMapGet().ntaxaGet()):
                retval += "%-10s" % self.taxonMapGet().labelGet(x)
                for y in forints(self.taxonMapGet().ntaxaGet()):
                    retval += (" " + distFormat) % self.distanceGet(x, y)
                retval += "\n"
        elif format == 'upper':
            for x in forints(self.taxonMapGet().ntaxaGet()):
                retval += "%-10s" % self.taxonMapGet().labelGet(x)
                for y in forints(x + 1):
                    retval += "%8s" % ""
                for y in forints(self.taxonMapGet().ntaxaGet(), start=x+1):
                    retval += (" " + distFormat) % self.distanceGet(x, y)
                retval += "\n"
        elif format == 'lower':
            for x in forints(self.taxonMapGet().ntaxaGet()):
                retval += "%-10s" % self.taxonMapGet().labelGet(x)
                for y in forints(x):
                    retval += (" " + distFormat) % self.distanceGet(x, y)
                retval += "\n"
        else:
            raise crux.DistMatrix\
                  .ValueError("Format must be 'full', 'upper', or 'lower'")

        return retval
#EOF
