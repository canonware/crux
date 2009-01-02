# Forward declarations.
cdef class CTMatrix
cdef class Alignment

from Crux.Taxa cimport Taxon
cimport Crux.Taxa as Taxa
from Crux.DistMatrix cimport DistMatrix

cdef class CTMatrix:
    cdef readonly type charType
    cdef readonly Taxa.Map taxaMap
    cdef readonly int sn
    cdef dict _taxonData

    cpdef _fastaNew(self, input, type charType)
    cpdef str fastaPrint(self, file outFile=*)
    cpdef str dataGet(self, Taxon taxon)
    cpdef dataSet(self, Taxon taxon, str data)

cdef class Alignment:
    cdef readonly type charType
    cdef readonly Taxa.Map taxaMap
    cdef readonly int ntaxa
    cdef readonly int nchars
    cdef char *rows
    cdef char *cols

    cdef char *_allocMatrix(self, int ntaxa, int nchars, str pad) except NULL

    cdef void setRow(self, int row, int col, char *chars, unsigned len)
    cdef char *getRow(self, int row)
    cdef char *getCol(self, int col)

    cpdef setSeq(self, int seq, int col, str chars)
    cpdef str getSeq(self, int seq)
    cpdef str getChar(self, int char_)

    cpdef str fastaPrint(self, file outFile=*)

    cpdef DistMatrix dists(self, bint avgAmbigs=*, bint scoreGaps=*)
    cpdef DistMatrix jukesDists(self, bint avgAmbigs=*, bint scoreGaps=*)
    cdef void _kimuraDistsDNA(self, DistMatrix m, bint scoreGaps)
    cdef void _kimuraDistsProtein(self, DistMatrix m, bint avgAmbigs,
      bint scoreGaps)
    cpdef DistMatrix kimuraDists(self, bint avgAmbigs=*, bint scoreGaps=*)
    cpdef DistMatrix logdetDists(self, bint scoreGaps=*)
