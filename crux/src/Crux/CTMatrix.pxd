from TaxonMap cimport TaxonMap

cdef class CTMatrix:
    cdef list _chars
    cdef TaxonMap _taxonMap
    cdef dict _taxonData
    cdef int _seq
