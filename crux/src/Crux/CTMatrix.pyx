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
from Taxa cimport Taxon
cimport Taxa
cimport Fasta
from Character cimport Dna, Protein

import sys

global __name__

cdef class CTMatrix:
    def __init__(self, input=None, type charType=Dna, taxaMap=None):
        self.charType = charType

        # This is used to manage the rows of the data matrix (_taxonData).
        if taxaMap == None:
            taxaMap = Taxa.Map()
        self.taxaMap = taxaMap

        # Initialize sequence number.  0 is skipped here so that it can be used
        # as a special value by users of CTMatrix.
        self.seq = 1

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

    cpdef str fastaPrint(self):
        cdef list lines, taxa
        cdef Taxon taxon
        cdef str taxonData
        cdef int i

        lines = []

        # Print.
        taxa = self.taxaMap.taxaGet()
        for taxon in taxa:
            if self.dataGet(taxon) != None:
                lines.append(">%s" % taxon.label.replace(' ', '_'))
                # Break into lines of length 75.
                taxonData = self.dataGet(taxon)
                for 0 <= i < len(taxonData) by 75:
                    if i + 75 < len(taxonData):
                        lines.append("%s" % taxonData[i:i+75])
                    else:
                        lines.append("%s" % taxonData[i:])

        return "\n".join(lines)

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

        self.seq += 1
