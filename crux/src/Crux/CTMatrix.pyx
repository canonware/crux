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

# Forward declarations.
cdef class Row(Fasta.Row)
cdef class _FastaParser(Fasta.Parser)
cdef class CTMatrix

#===============================================================================
# Begin FASTA construction support.
#

cdef class Row(Fasta.Row):
    "%extend Row"
    cpdef reduce(self, Fasta.TokenDescr descr, Fasta.TokenChars chars):
        "%accept"
        self.parser._addTaxon(Taxa.get(descr.label), chars.chars)

cdef Parsing.Spec _FastaSpec

cdef class _FastaParser(Fasta.Parser):
    cdef CTMatrix matrix
    cdef Taxa.Map taxaMap

    def __init__(self, CTMatrix matrix, Taxa.Map taxaMap):
        global _FastaSpec

        if _FastaSpec is None:
            _FastaSpec = self._initSpec()
        Fasta.Parser.__init__(self, _FastaSpec)

        self.matrix = matrix
        self.taxaMap = taxaMap

    cdef Parsing.Spec _initSpec(self):
        return Parsing.Spec([sys.modules[__name__], Crux.Fasta],
          pickleFile="%s/share/Crux-%s/CTMatrix.pickle" %
          (Crux.Config.prefix, Crux.Config.version),
          verbose=(False if (not __debug__ or Crux.opts.quiet) else True),
          skinny=(False if __debug__ else True),
          logFile="%s/share/Crux-%s/CTMatrix.log" %
          (Crux.Config.prefix, Crux.Config.version))

    cpdef _addTaxon(self, Taxon taxon, str chars):
        if self.taxaMap.indGet(taxon) == -1:
            # Define a taxon mapping for this label.
            self.taxaMap.map(taxon, self.taxaMap.ntaxa)

        # Set the character data for this taxon.
        self.matrix.dataSet(taxon, chars)

#
# Begin FASTA construction support.
#===============================================================================

cdef class CTMatrix:
    def __init__(self, str input=None, type charType=Dna, taxaMap=None):
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

    cpdef _fastaNew(self, str input, type charType):
        cdef _FastaParser parser

        parser = _FastaParser(self, self.taxaMap)
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
