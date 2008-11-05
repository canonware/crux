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
from TaxonMap cimport TaxonMap
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
        self.parser._addTaxon(descr.label, chars.chars)

cdef Parsing.Spec _FastaSpec

cdef class _FastaParser(Fasta.Parser):
    cdef CTMatrix matrix
    cdef TaxonMap taxonMap

    def __init__(self, CTMatrix matrix, TaxonMap taxonMap):
        global _FastaSpec

        if _FastaSpec is None:
            _FastaSpec = self._initSpec()
        Fasta.Parser.__init__(self, _FastaSpec)

        self.matrix = matrix
        self.taxonMap = taxonMap

    cdef Parsing.Spec _initSpec(self):
        return Parsing.Spec([sys.modules[__name__], Crux.Fasta],
          pickleFile="%s/share/Crux-%s/CTMatrix.pickle" %
          (Crux.Config.prefix, Crux.Config.version),
          verbose=(False if (not __debug__ or Crux.opts.quiet) else True),
          skinny=(False if __debug__ else True),
          logFile="%s/share/Crux-%s/CTMatrix.log" %
          (Crux.Config.prefix, Crux.Config.version))

    cpdef _addTaxon(self, str label, str chars):
        if self.taxonMap.indGet(label) == -1:
            # Define a taxon mapping for this label.
            self.taxonMap.map(label, self.taxonMap.ntaxaGet())

        # Set the character data for this taxon.
        self.matrix.dataSet(label, chars)

#
# Begin FASTA construction support.
#===============================================================================

cdef class CTMatrix:
    def __init__(self, str input=None, type charType=Dna, taxonMap=None):
        self.charType = charType

        # This is used to manage the rows of the data matrix (_taxonData).
        if taxonMap == None:
            taxonMap = TaxonMap()
        self.taxonMap = taxonMap

        # Initialize sequence number.  0 is skipped here so that it can be used
        # as a special value by users of CTMatrix.
        self.seq = 1

        # Row-major character data.  Each element in _taxonData is a string of
        # characters that belong to the corresponding taxon in taxonMap.  The
        # keys are the integer indices, as reported by self.taxonMap.indGet().
        self._taxonData = {}

        if input is not None:
            self._fastaNew(input, charType)

    cpdef _fastaNew(self, str input, type charType):
        cdef _FastaParser parser

        parser = _FastaParser(self, self.taxonMap)
        parser.parse(input, charType)

    cpdef str fastaPrint(self):
        cdef list lines
        cdef str taxonData
        cdef int i

        lines = []

        # Print.
        taxa = self.taxonMap.taxaGet()
        for taxon in taxa:
            if self.dataGet(taxon) != None:
                lines.append(">%s" % taxon.replace(' ', '_'))
                # Break into lines of length 75.
                taxonData = self.dataGet(taxon)
                for 0 <= i < len(taxonData) by 75:
                    if i + 75 < len(taxonData):
                        lines.append("%s" % taxonData[i:i+75])
                    else:
                        lines.append("%s" % taxonData[i:])

        return "\n".join(lines)

    # Return the character data for a taxon.
    cpdef str dataGet(self, str taxon):
        if not self._taxonData.has_key(self.taxonMap.indGet(taxon)):
            rVal = None
        else:
            rVal = self._taxonData[self.taxonMap.indGet(taxon)]

        return rVal

    # Set the character data for a taxon.
    cpdef dataSet(self, str taxon, str data):
        if (self.taxonMap.indGet(taxon) == None):
            raise ValueError("Taxon %r not in taxon map" % taxon)

        self._taxonData[self.taxonMap.indGet(taxon)] = data

        self.seq += 1
