cdef class TaxonMap:
    cdef dict _label2ind, _ind2label
    cpdef int ntaxaGet(self)
    cpdef str labelGet(self, int ind)
    cpdef int indGet(self, str label)
    cpdef list taxaGet(self)
    cpdef bint equal(self, TaxonMap other)
    cpdef map(self, str label, int ind, bint replace=?)
