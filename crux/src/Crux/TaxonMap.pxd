cdef class TaxonMap:
    cdef dict _label2ind, _ind2label
    cpdef int ntaxaGet(self)
    cpdef labelGet(self, int ind)
    cpdef indGet(self, label)
    cpdef list taxaGet(self)
    cpdef bint equal(self, other)
    cpdef map(self, label, int ind, bint replace=?)
