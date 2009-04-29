cdef class Character:
    cdef dict _pStates
    cdef dict _aStates
    cdef dict _vals
    cdef readonly int any
    cdef readonly int nstates

    cpdef list pcodes(self)
    cpdef list codes(self)
    cdef void _stateCodeAdd(self, str code) except *
    cdef void _ambiguityCodeAdd(self, str code, list oStates) except *
    cdef void _aliasCodeAdd(self, str code, str oState) except *
    cpdef int code2val(self, str code) except -1
    cpdef str val2code(self, int val)

cdef class Dna(Character):
    pass
    # get() returns a singleton.

cdef class Protein(Character):
    pass
    # get() returns a singleton.
