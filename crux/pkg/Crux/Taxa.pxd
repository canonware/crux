# Do not directly construct Taxon instances; instead call the taxon() function.
cdef class Taxon:
    cdef object __weakref__
    cdef readonly str label

cpdef Taxon get(str label)

cdef class Map:
    cdef dict _taxon2ind, _ind2taxon
    cdef readonly int ntaxa

    cpdef Taxon taxonGet(self, int ind)
    cpdef int indGet(self, Taxon taxon)
    cdef list taxaGet(self)
    # property taxa
    cpdef map(self, Taxon taxon, int ind, bint replace=*)
