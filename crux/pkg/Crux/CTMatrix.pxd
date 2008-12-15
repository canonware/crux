# Forward declaration.
cdef class CTMatrix

from Crux.Taxa cimport Taxon
cimport Crux.Taxa as Taxa
from Crux.DistMatrix cimport DistMatrix

cdef class CTMatrix:
    cdef readonly charType
    cdef readonly Taxa.Map taxaMap
    cdef readonly int seq
    cdef dict _taxonData

    cpdef _fastaNew(self, input, type charType)
    cpdef str fastaPrint(self, file outFile=*)
    cpdef str dataGet(self, Taxon taxon)
    cpdef dataSet(self, Taxon taxon, str data)
    cdef float _dist(self, char *a, char *b, int sLen)
    cpdef DistMatrix distances(self, str correction=*)
