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
from Crux.Taxa cimport Taxon
cimport Crux.Taxa as Taxa
cimport Crux.Fasta as Fasta
from Crux.Character cimport Dna, Protein
from Crux.DistMatrix cimport DistMatrix

from CxDistMatrix cimport *

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

    cpdef str fastaPrint(self, file outFile=None):
        cdef list lines, taxa
        cdef Taxon taxon
        cdef str s, taxonData
        cdef int i

        if outFile is None:
            lines = []

        # Print.
        taxa = self.taxaMap.taxaGet()
        for taxon in taxa:
            if self.dataGet(taxon) != None:
                s = ">%s\n" % taxon.label.replace(' ', '_')
                if outFile is None:
                    lines.append(s)
                else:
                    outFile.write(s)

                # Break into lines of length 75.
                taxonData = self.dataGet(taxon)
                for 0 <= i < len(taxonData) by 75:
                    if i + 75 < len(taxonData):
                        s = "%s\n" % taxonData[i:i+75]
                    else:
                        s = "%s\n" % taxonData[i:]

                    if outFile is None:
                        lines.append(s)
                    else:
                        outFile.write(s)

        if outFile is None:
            if len(lines) > 0:
                # Chop off the last '\n'.
                lines[-1] = lines[-1][:-1]
            return "".join(lines)
        else:
            return None

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

    # Compute the uncorrected distance between two taxa.
    cdef CxtDMDist _dist(self, char *a, char *b, int sLen):
        cdef int i, nchars, diff
        cdef char aC, bC

        nchars = 0
        diff = 0
        for 0 <= i < sLen:
            aC = a[i]
            bC = b[i]
            if (aC != c'-' and aC != c'?') or (bC != c'-' and bC != c'?'):
                nchars += 1
                if aC != bC:
                    diff += 1

        if nchars > 0:
            return <CxtDMDist>diff / <CxtDMDist>nchars
        else:
            return 1.0

    cpdef DistMatrix distances(self, str correction=None):
        """
Calculate pairwise distances using an optional correction for multiple hits.
The following correction algorithms are supported:

  "jukes" : XXX
  "kimura" : XXX
"""
        cdef DistMatrix ret
        cdef int ntaxa, i, j, iSeqLen, jSeqLen
        cdef Taxon iTaxon, jTaxon
        cdef str iSeq, jSeq
        cdef char *iStr, *jStr
        cdef CxtDMDist dist

        assert correction in (None, "jukes", "kimura")

        # XXX
        if correction in ("jukes", "kimura"):
            print "XXX Not implemented"
            sys.exit(1)

        ret = DistMatrix(self.taxaMap)

        ntaxa = self.taxaMap.ntaxa
        if ntaxa > 1:
            for 0 <= i < ntaxa:
                iTaxon = self.taxaMap.taxonGet(i)
                iSeq = self.dataGet(iTaxon)
                iSeqLen = len(iSeq)
                iStr = <char *>iSeq
                for i + 1 <= j < ntaxa:
                    jTaxon = self.taxaMap.taxonGet(j)
                    jSeq = self.dataGet(jTaxon)
                    jSeqLen = len(jSeq)
                    if iSeqLen != jSeqLen:
                        raise ValueError("Sequence lengths differ (%d vs. %d)" %
                          (iSeqLen, jSeqLen))
                    jStr = <char *>jSeq
                    dist = self._dist(iStr, jStr, iSeqLen)
                    ret.distanceSet(i, j, dist)

        return ret
