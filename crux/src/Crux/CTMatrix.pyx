import Crux.Exception

class Exception(Crux.Exception.Exception):
    pass

import exceptions

class ValueError(Exception, exceptions.ValueError):
    def __init__(self, str):
        self._str = str

    def __str__(self):
        return self._str

from TaxonMap cimport TaxonMap
from Fasta cimport Fasta
from Character cimport Dna, Protein

class _Fasta(Fasta):
    def __init__(self, matrix, taxonMap):
        self._matrix = matrix
        self._taxonMap = taxonMap
        self._lastLabel = None

    # Overridden method.
    def labelAccept(self):
        self._lastLabel = self.token()

    # Overridden method.
    def charsAccept(self):
        # Set the character data type if it hasn't already been done.
        if len(self._matrix.charsGet()) == 0:
            if self.charType() == 'DNA':
                self._matrix.charsAppend([Dna()]
                  * len(self.token()))
            else: # Protein data.
                self._matrix.charsAppend([Protein()]
                  * len(self.token()))

        if self._taxonMap.indGet(self._lastLabel) == None:
            # Define a taxon mapping for this label.
            self._taxonMap.map(self._lastLabel, self._taxonMap.ntaxaGet())

        # Set the character data for this taxon.
        self._matrix.dataSet(self._lastLabel, self.token())

cdef class CTMatrix:
    def __init__(self, taxonMap=None):
        # This is an array of Character objects, and stores the character
        # types for the columns in _taxonData.
        self._chars = []

        # This is used to manage the rows of the data matrix (_taxonData).
        if taxonMap == None:
            taxonMap = TaxonMap()
        self._taxonMap = taxonMap

        # Row-major character data.  Each element in _taxonData is a string of
        # characters that belong to the corresponding taxon in _taxonMap.  The
        # keys are the integer indices, as reported by self._taxonMap.indGet().
        self._taxonData = {}

        # Initialize sequence number.  0 is skipped here so that it can be used
        # as a special value by users of CTMatrix.
        self._seq = 1

    # Return the current sequence number.  This number is incremented every
    # time that the matrix is modified.
    def seqGet(self):
        return self._seq

    def fastaParse(self, input, chartype='DNA'):
        parser = _Fasta(self, self._taxonMap)
        parser.parse(input, chartype)

    def fastaPrint(self, outFile=None):
        # Set up for sending output to either a string or a file.
        if outFile == None:
            callback = self._stringPrintCallback
            self._string = ""
        else:
            callback = self._filePrintCallback
            self._outFile = outFile

        # Print.
        taxa = self._taxonMap.taxaGet()
        for taxon in taxa:
            if self.dataGet(taxon) != None:
                callback(">%s\n" % taxon)
                # Break into lines of length 75.
                for i in xrange(0, len(self.dataGet(taxon)), 75):
                    if i + 75 < len(self.dataGet(taxon)):
                        callback("%s\n" % self.dataGet(taxon)[i:i+75])
                    else:
                        callback("%s\n" % self.dataGet(taxon)[i:])

        # Clean up and set rVal according to where the output was sent.
        if outFile == None:
            rVal = self._printString
            self._printString = None
        else:
            rVal = None
            self._printOutFile = None
        return rVal

    # Append Character objects to _chars.
    def charsAppend(self, chars):
        # XXX Only allow this to succeed if no character data have been set
        # yet.
        self._chars += chars
        self._seq += 1

    # Return the array of character types.
    def charsGet(self):
        return self._chars

    # Return the taxon map.
    def taxonMapGet(self):
        return self._taxonMap

    # Return the character data for a taxon.
    def dataGet(self, taxon):
        if not self._taxonData.has_key(self._taxonMap.indGet(taxon)):
            rVal = None
        else:
            rVal = self._taxonData[self._taxonMap.indGet(taxon)]

        return rVal

    # Set the character data for a taxon.
    def dataSet(self, taxon, data):
        if len(data) != len(self._chars):
            raise ValueError(
              "Wrong number of characters (%d, expected %d)"
                              % (len(data), len(self._chars)))
        if (self._taxonMap.indGet(taxon) == None):
            raise ValueError("Taxon %r not in taxon map" % taxon)

        self._taxonData[self._taxonMap.indGet(taxon)] = data

        self._seq += 1

    # Callback method that is used by the print method for printing of the
    # CTMatrix in FASTA format.
    def _stringPrintCallback(self, string):
        # Append string to previous strings that were passed to this callback.
        self._string = "%s%s" % (self._string, string)

    # Callback method that is used by the print method for printing of the
    # CTMatrix in FASTA format.
    def _filePrintCallback(self, string):
        # Print string to self._outFile.
        self._outFile.write(string)
