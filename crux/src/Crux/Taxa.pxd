# Do not directly construct Taxon instances; instead call the taxon() function.
cdef class Taxon:
    cdef object __weakref__
    cdef readonly str label

cpdef Taxon get(str label)

cdef class Map:
    cdef dict _taxon2ind, _ind2taxon

    cpdef int ntaxaGet(self)
    cpdef Taxon taxonGet(self, int ind)
    cpdef int indGet(self, Taxon taxon)
    cpdef list taxaGet(self)
    cpdef bint equal(self, Map other) # XXX Use __eq__ or similar.
    cpdef map(self, Taxon taxon, int ind, bint replace=?)
