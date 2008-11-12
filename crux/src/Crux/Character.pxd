cdef class Character:
    cdef dict _pStates
    cdef dict _aStates
    cdef dict _vals

    cpdef int nstates(self)
    cpdef list codes(self)
    cpdef stateCodeAdd(self, str code)
    cpdef ambiguityCodeAdd(self, str code, list oStates=*)
    cpdef aliasCodeAdd(self, str code, str oState)
    cpdef int code2val(self, str code)
    cpdef str val2code(self, int val)

cdef class Dna(Character): pass

cdef class Protein(Character): pass
