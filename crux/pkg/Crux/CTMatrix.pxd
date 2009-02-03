# Forward declarations.
cdef class CTMatrix
cdef class Alignment

from Crux.Character cimport Character
from Crux.Taxa cimport Taxon
cimport Crux.Taxa as Taxa
from Crux.DistMatrix cimport DistMatrix

cdef class CTMatrix:
    cdef readonly type charType
    cdef readonly Taxa.Map taxaMap
    cdef readonly int sn
    cdef dict _taxonData

    cpdef _fastaNew(self, input, type charType)
    cdef void _renderLine(self, str line, list lines, file outFile) except *
    cpdef str fastaPrint(self, file outFile=*)
    cpdef str dataGet(self, Taxon taxon)
    cpdef dataSet(self, Taxon taxon, str data)

cdef class Alignment:
    cdef readonly type charType
    cdef readonly Taxa.Map taxaMap
    cdef readonly int ntaxa
    cdef readonly int nchars
    cdef readonly bint rowMajor
    cdef readonly bint colMajor
    cdef char *rows
    cdef char *cols
    cdef unsigned *freqs

    cdef char *_allocMatrix(self, int ntaxa, int nchars, str pad) except NULL
    cdef void _initRowMajor(self) except *
    cdef void _initColMajor(self) except *

    cpdef pad(self, str pad, int npad)

    cpdef deRow(self)
    cpdef deCol(self)

    cdef void setRow(self, int row, int col, char *chars, unsigned len) except *
    cdef void setCol(self, int row, int col, char *chars, unsigned len) except *
    cdef char *getRow(self, int row) except NULL
    cdef char *getCol(self, int col) except NULL

    cpdef setSeq(self, int seq, int col, str chars)
    cpdef setChar(self, int seq, int char_, str chars)
    cpdef str getSeq(self, int seq)
    cpdef str getChar(self, int char_)

    cpdef unsigned getFreq(self, int col)
    cpdef setFreq(self, int col, unsigned freq)

    cdef void _fitchCanonize(self, Character char_) except *
    cpdef canonize(self, bint fitch=*)
    cdef int _fitchCompact(self, list deco) except -1
    cpdef int compact(self, bint fitch=*)

    cdef void _renderLine(self, str line, list lines, file outFile) except *
    cpdef str render(self, unsigned interleave=*, file outFile=*)
    cpdef str fastaPrint(self, file outFile=*)

    cpdef DistMatrix dists(self, bint scoreGaps=*)
    cpdef DistMatrix jukesDists(self, bint scoreGaps=*)
    cdef void _kimuraDistsDNA(self, DistMatrix m, bint scoreGaps)
    cdef void _kimuraDistsProtein(self, DistMatrix m, bint scoreGaps)
    cpdef DistMatrix kimuraDists(self, bint scoreGaps=*)
    cpdef DistMatrix logdetDists(self, bint scoreGaps=*, bint allowNan=*)
