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

import string

import TaxonMap
import crux.Exception

class Exception(crux.Exception):
    pass

class SyntaxError(Exception, SyntaxError):
    def __init__(self, line, message, format=None):
        self._line = line
        self._format = format
        self._message = message

    def __str__(self):
        if self._format != None:
            retval = "Line %d (%s matrix format): %s" \
                     % (self._line, self._format, self._message)
        else:
            retval = "Line %d: %s" \
                     % (self._line, self._message)

        return retval

class ValueError(Exception, ValueError):
    def __init__(self, str):
        self._str = str

    def __str__(self):
        return self._str

class DistMatrix(object):
    # Construct a DistMatrix from one of the following inputs:
    #
    #   str : Parse the string as a distance matrix.
    #
    #   file : Parse the file as a distance matrix.
    #
    #   TaxonMap : Create an uninitialized distance matrix of the appropriate
    #              size, given the number of taxa in the TaxonMap.
    def __init__(self, input=None):
        self._ntaxa = 0
        self._map = None
        self._matrix = None

        if type(input) == crux.TaxonMap.TaxonMap:
            self._taxonMapNew(input)
        elif type(input) == str:
            self._strNew(input)
        elif type(input) == file:
            self._fileNew(input)
        elif input != None:
            raise crux.DistMatrix.ValueError("Input must be string or file")

        return (self._map, self._matrix)

    def ntaxaGet(self):
        return self._ntaxa

    def taxonMapGet(self):
        return self._map

    def distanceGet(self, fr, to):
        if fr >= self._ntaxa or to >= self._ntaxa:
            raise crux.DistMatrix.ValueError("Out of bounds matrix access")

        return self._distanceGet(fr, to)

    def _distanceGet(self, fr, to):
        return self._matrix[fr * self._ntaxa + to]

    def distanceSet(self, fr, to, distance):
        if fr >= self._ntaxa or to >= self._ntaxa:
            raise crux.DistMatrix.ValueError("Out of bounds matrix access")

        if type(distance) != float and type(distance) != int:
            raise crux.DistMatrix.ValueError("Distance must be a number")

        self._distanceSet(fr, to, distance)

    def _distanceSet(self, fr, to, distance):
        self._matrix[fr * self._ntaxa + to] = distance

    # Print the matrix to a string in 'full', 'upper', or 'lower' format.
    def prints(self, format='full'):
        retval = "%d\n" % self._ntaxa
        if format == 'full':
            for x in forints(self._ntaxa):
                retval += "%-10s" % self._map.labelGet(x)
                for y in forints(self._ntaxa):
                    retval += " %1.5f" % self._distanceGet(x, y)
                retval += "\n"
        elif format == 'upper':
            for x in forints(self._ntaxa):
                retval += "%-10s" % self._map.labelGet(x)
                for y in forints(x + 1):
                    retval += "%8s" % ""
                for y in forints(self._ntaxa, start=x+1):
                    retval += " %1.5f" % self._distanceGet(x, y)
                retval += "\n"
        elif format == 'lower':
            for x in forints(self._ntaxa):
                retval += "%-10s" % self._map.labelGet(x)
                for y in forints(x):
                    retval += " %1.5f" % self._distanceGet(x, y)
                retval += "\n"
        else:
            raise crux.DistMatrix\
                  .ValueError("Format must be 'full', 'upper', or 'lower'")

        return retval

    # Create an empty distance matrix that is the right size for the TaxonMap
    # that was passed in to the constructor.
    def _taxonMapNew(self, input):
        # The size of the distance matrix corresponds to the number of taxa in
        # the TaxonMap that was passed in.
        self._ntaxa = input.ntaxaGet()

        # Use the TaxonMap that was passed in as the map.
        self._map = input

        # Create an uninitialized distance matrix.
        self._matrix = [None] * (self._ntaxa * self._ntaxa)

    # Parse input (string) and return a tuple, where the first element in the
    # tuple is a TaxonMmap, and the second element is a row-major matrix.
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
    # (<TaxonMap: ['Taxon_A', 'Taxon_B', 'Taxon_C', 'Taxon_D', 'Taxon_E']>,
    #  [0.0, 1.0, 2.0, 3.0, 4.0,
    #   1.0, 0.0, 1.5, 2.5, 3.5,
    #   2.0, 1.5, 0.0, 2.2, 3.2,
    #   3.0, 2.5, 2.2, 0.0, 3.1,
    #   4.0, 3.5, 3.2, 3.1, 0.0])
    #
    def _strNew(self, input):
        self._tokenGet = self._strTokenGet

        self._input = input
        self._i = 0

        self._parse()

    def _fileNew(self, input):
        self._tokenGet = self._fileTokenGet

        self._input = input
        self._tokenBuf = [None] * 32 # This is iteratively doubled as necessary.

        self._parse()

    def _parse(self):
        import __builtin__

        self._line = 1
        self._matrixFormat = 'unknown' # 'unknown', 'full', 'upper', 'lower'

        # Get the number of taxa.
        (token, line) = self._tokenGet()
        try:
            self._ntaxa = int(token)
        except __builtin__.ValueError:
            raise crux.DistMatrix.SyntaxError(line,
                                              "Unspecified number of taxa")
        if self._ntaxa < 2:
            raise crux.DistMatrix.SyntaxError(line,
                                              "Too few taxa")

        # Create an empty TaxonMap.
        self._map = TaxonMap.TaxonMap()

        # Create an empty distance matrix.
        self._matrix = [None] * (self._ntaxa * self._ntaxa)

        # Get the first taxon label.
        (token, line) = self._tokenGet()
        distance = self._tokenToDistance(token)
        if distance != None:
            raise crux.DistMatrix.SyntaxError(line, "Missing taxon label")
        self._map.map(token, self._map.ntaxaGet())

        # Get the next token; if it is a taxon label, then this matrix is in
        # lower triangle format.
        (token, line) = self._tokenGet()
        distance = self._tokenToDistance(token)
        if distance == None:
            # This is a lower-triangle matrix.
            self._map.map(token, self._map.ntaxaGet())

            # Get second row of distances.
            for y in forints(1):
                (token, line) = self._tokenGet()
                distance = self._tokenToDistance(token)
                if distance == None:
                    raise crux.DistMatrix\
                          .SyntaxError(line,
                                       "Missing distance (%d, %d)" % (1, y),
                                       'lower')
                self._distanceSet(1, y, distance)

            # Get remaining rows.
            for x in forints(self._ntaxa, start=2):
                # Get taxon label.
                (token, line) = self._tokenGet()
                distance = self._tokenToDistance(token)
                if distance != None:
                    raise crux.DistMatrix.SyntaxError(line,
                                                      "Missing taxon label",
                                                      'lower')
                self._map.map(token, self._map.ntaxaGet())

                # Get distances.
                for y in forints(x):
                    (token, line) = self._tokenGet()
                    distance = self._tokenToDistance(token)
                    if distance == None:
                        raise crux.DistMatrix\
                              .SyntaxError(line,
                                           "Missing distance (%d, %d)" % (x, y),
                                           'lower')
                    self._distanceSet(x, y, distance)

            # Reflect matrix contents.
            for x in forints(self._ntaxa):
                for y in forints(self._ntaxa, x + 1):
                    self._distanceSet(x, y, self._distanceGet(y, x))
            # Initialize diagonal.
            for x in forints(self._ntaxa):
                self._distanceSet(x, x, 0.0)
        else:
            # Get the first row of distances, and insert them into the matrix
            # as though parsing an upper-triangle matrix.
            self._distanceSet(0, 1, distance)

            for y in forints(self._ntaxa, start=2):
                (token, line) = self._tokenGet()
                distance = self._tokenToDistance(token)
                if distance == None:
                    raise crux.DistMatrix\
                          .SyntaxError(line,
                                       "Missing distance (%d, %d)" % (0, y))
                self._distanceSet(0, y, distance)

            # Determine whether this is a full or upper-triangle matrix.
            (token, line) = self._tokenGet()
            distance = self._tokenToDistance(token)
            if distance == None:
                # This is an upper-triangle matrix.
                self._map.map(token, self._map.ntaxaGet())

                # Get second row of distances.
                for y in forints(self._ntaxa, start=2):
                    (token, line) = self._tokenGet()
                    distance = self._tokenToDistance(token)
                    if distance == None:
                        raise crux.DistMatrix\
                              .SyntaxError(line,
                                           "Missing distance (%d, %d)" % (1, y),
                                           'upper')
                    self._distanceSet(1, y, distance)

                # Get remaining rows.
                for x in forints(self._ntaxa, start=2):
                    # Get taxon label.
                    (token, line) = self._tokenGet()
                    distance = self._tokenToDistance(token)
                    if distance != None:
                        raise crux.DistMatrix.SyntaxError(line,
                                                          "Missing taxon label",
                                                          'upper')
                    self._map.map(token, self._map.ntaxaGet())

                    # Get distances.
                    for y in forints(self._ntaxa, start=x+1):
                        (token, line) = self._tokenGet()
                        distance = self._tokenToDistance(token)
                        if distance == None:
                            raise crux.DistMatrix\
                                  .SyntaxError(line,
                                               "Missing distance (%d, %d)"
                                               % (x, y),
                                               'upper')
                        self._distanceSet(x, y, distance)

                # Reflect matrix contents.
                for x in forints(self._ntaxa):
                    for y in forints(self._ntaxa, x + 1):
                        self._distanceSet(y, x, self._distanceGet(x, y))
                # Initialize diagonal.
                for x in forints(self._ntaxa):
                    self._distanceSet(x, x, 0.0)
            else:
                # This is a full matrix.

                # Shift the contents of the first row back one position.
                for y in forints(self._ntaxa - 1):
                    self._distanceSet(0, y, self._distanceGet(0, y + 1))

                # Set last distance on first row.
                self._distanceSet(0, self._ntaxa - 1, distance)

                # Get remaining rows.
                for x in forints(self._ntaxa, start=1):
                    # Get taxon label.
                    (token, line) = self._tokenGet()
                    distance = self._tokenToDistance(token)
                    if distance != None:
                        raise crux.DistMatrix.SyntaxError(line,
                                                          "Missing taxon label",
                                                          'full')
                    self._map.map(token, self._map.ntaxaGet())

                    # Get distances.
                    for y in forints(self._ntaxa):
                        (token, line) = self._tokenGet()
                        distance = self._tokenToDistance(token)
                        if distance == None:
                            raise crux.DistMatrix\
                                  .SyntaxError(line,
                                               "Missing distance (%d, %d)"
                                               % (x, y),
                                               'full')
                        self._distanceSet(x, y, distance)

    # Return the next token.
    def _strTokenGet(self):
        token = ""
        start = self._i
        line = self._line
        while self._i < len(self._input):
            c = self._input[self._i]

            if c == " " or c == "\n" or c == "\t":
                if c == "\n":
                    self._line += 1

                if self._i == start:
                    # Nothing but whitespace so far.
                    start = self._i + 1
                    line = self._line
                else:
                    token = self._input[start:self._i]
                    self._i += 1
                    break

            self._i += 1

        return (token, line)

    # Return the next token.
    def _fileTokenGet(self):
        token = ""
        tokenLen = 0
        line = self._line

        c = self._input.read(1)
        while len(c) == 1:
            if c == " " or c == "\n" or c == "\t":
                if c == "\n":
                    self._line += 1

                if tokenLen == 0:
                    # Nothing but whitespace so far.
                    line = self._line
                else:
                    token = string.join(self._tokenBuf[:tokenLen], "")
                    tokenLen = 0
                    break
            else:
                # Enlarge the buffer, if necessary.
                if len(self._tokenBuf) == tokenLen:
                    self._tokenBuf.extend([None] * tokenLen)

                self._tokenBuf[tokenLen] = c
                tokenLen += 1

            c = self._input.read(1)

        return (token, line)

    # Return distance (float), or None if the token cannot be converted to a
    # distance.
    def _tokenToDistance(self, token):
        import __builtin__

        try:
            retval = float(token)
        except __builtin__.ValueError:
            retval = None

        return retval
#EOF
