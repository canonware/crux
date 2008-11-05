from TaxonMap cimport TaxonMap

cdef class CTMatrix:
    cdef readonly charType
    cdef readonly TaxonMap taxonMap
    cdef readonly int seq
    cdef dict _taxonData

    cpdef _fastaNew(self, str input, type charType)
    cpdef str fastaPrint(self)
    cpdef str dataGet(self, str taxon)
    cpdef dataSet(self, str taxon, str data)
