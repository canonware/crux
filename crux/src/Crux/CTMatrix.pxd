from Taxa cimport Taxon
cimport Taxa

cdef class CTMatrix:
    cdef readonly charType
    cdef readonly Taxa.Map taxaMap
    cdef readonly int seq
    cdef dict _taxonData

    cpdef _fastaNew(self, str input, type charType)
    cpdef str fastaPrint(self)
    cpdef str dataGet(self, Taxon taxon)
    cpdef dataSet(self, Taxon taxon, str data)
