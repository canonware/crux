# Forward declarations.
cdef class Token
cdef class Nonterm
cdef class Parser

cimport Parsing

cdef class Token(Parsing.Token):
    cdef readonly Token prev, next
    cdef str _input
    cdef readonly int begPos, endPos
    cdef readonly int line, col

    # property raw

cdef class Nonterm(Parsing.Nonterm):
    cdef readonly int begPos, endPos
    cdef readonly str variant

cdef class Parser(Parsing.Lr):
    cdef readonly Token first, last

    cdef void _appendToken(self, Token token) except *
    cdef str expandInput(self, str input, int pos, int line, int col)
    cpdef parse(self, str input, int begPos=?, int line=?, int col=?, \
      bint verbose=?)
