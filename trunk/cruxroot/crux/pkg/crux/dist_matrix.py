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

import taxon_map

import re

class dist_matrix(object):
    def __init__(self):
        pass

    # Parse input (file or string) and return a tuple, where the first element
    # in the tuple is a taxon_map, and the second element is a
    # row-major matrix.
    #
    # Example input:
    #
    #   5
    #   Taxon_A     1.0 2.0 3.0 4.0
    #   Taxon_B         1.5 2.5 3.5
    #   Taxon_C             2.2 3.2
    #   Taxon_D                 3.1
    #   Taxon_E
    #
    # Corresponding return value:
    #
    # (<taxon_map: ['Taxon_A', 'Taxon_B', 'Taxon_C', 'Taxon_D', 'Taxon_E']>,
    #  [0.0, 1.0, 2.0, 3.0, 4.0,
    #   1.0, 0.0, 1.5, 2.5, 3.5,
    #   2.0, 1.5, 0.0, 2.2, 3.2,
    #   3.0, 2.5, 2.2, 0.0, 3.1,
    #   4.0, 3.5, 3.2, 3.1, 0.0])
    #
    def parse(self, input):
        self._input = input
        self._i = 0
        self._matrix_format = 'unknown' # 'unknown', 'full', 'upper', 'lower'

        # Get the number of taxa.
        token = self._get_token()
        self._ntaxa = int(token)

        # Create an empty taxon_map.
        self._map = taxon_map.taxon_map()

        # Create an empty distance matrix (fourth quadrant coordinates).
        self._matrix = [None] * (self._ntaxa * self._ntaxa);

        # Get the first taxon label.
        token = self._get_token()
        distance = self._token_to_distance(token)
        if distance != None:
            raise ValueError
        self._map.map(token, self._map.ntaxa_get())
        distances = []

        x = 0
        while True:
            #sys.stderr.write("A")
            token = self._get_token()
            #sys.stderr.write("B")
            distance = self._token_to_distance(token)
            #sys.stderr.write("C")
            if distance != None:
                # Get distance for taxon.
                #sys.stderr.write("D")
                distances.append(distance)
            else:
                # Merge distances into the matrix.
                self._merge_distances(distances, x)
                distances = []

                # Get taxon label, unless all taxon labels have already been
                # read.
                if self._map.ntaxa_get() == self._ntaxa:
                    break
                else:
                    self._map.map(token, self._map.ntaxa_get())
                    x += 1

        return (self._map, self._matrix)

    # Return the next token.
    #
    # XXX Add support for input files.
    def _get_token(self):
        retval = ""
        start = self._i
        while self._i < len(self._input):
            if self._input[self._i] == " " \
                   or self._input[self._i] == "\n" \
                   or self._input[self._i] == "\t":
                if self._i == start:
                    start = self._i + 1
                else:
                    retval = self._input[start:self._i]
                    break

            self._i += 1

        return retval

    # Return distance (float), or None if the token cannot be converted to a
    # distance.
    def _token_to_distance(self, token):
        try:
            retval = float(token)
        except ValueError:
            retval = None

        return retval

    # Merge distances into self._matrix.
    def _merge_distances(self, distances, row):
        if self._matrix_format == 'unknown':
            if len(distances) == self._ntaxa:
                self._matrix_format = 'full'
            elif len(distances) == row:
                self._matrix_format = 'lower'
            elif len(distances) == self._ntaxa - row - 1:
                self._matrix_format = 'upper'

        if self._matrix_format == 'full':
            if len(distances) != self._ntaxa:
                raise ValueError
            i = 0
            while i < len(distances):
                self._matrix[row * self._ntaxa + i] = distances[i]
                i += 1
        elif self._matrix_format == 'lower':
            if len(distances) != row:
                # XXX Add informative error messages?
                print "XXX row: %d, len: %d" % (row, len(distances))
                raise ValueError
            self._matrix[row * self._ntaxa + row] = 0.0
            i = 0
            while i < len(distances):
                self._matrix[row * self._ntaxa + i] = distances[i]
                self._matrix[i * self._ntaxa + row] = distances[i]
                i += 1
        elif self._matrix_format == 'upper':
            if len(distances) != self._ntaxa - row - 1:
                raise ValueError
            self._matrix[row * self._ntaxa + row] = 0.0
            i = 0
            while i < len(distances):
                self._matrix[(i + row + 1) * self._ntaxa + row] = distances[i]
                self._matrix[row * self._ntaxa + (i + row + 1)] = distances[i]
                i += 1
        else:
            # Invalid matrix format.
            raise ValueError
#EOF
