cimport Parsing

# Forward declarations.
cdef class Token(Parsing.Token)
cdef class TokenDescr(Token)
cdef class TokenChars(Token)
cdef class Nonterm(Parsing.Nonterm)
cdef class Matrix(Nonterm)
cdef class RowList(Nonterm)
cdef class Row(Nonterm)
cdef class Parser(Parsing.Lr)

#===============================================================================
# Begin Token.
#
cdef class Token(Parsing.Token):
    cdef readonly Token prev, next
    cdef str _input
    cdef readonly int begPos, endPos
    cdef readonly int line
    # property raw

cdef class TokenDescr(Token):
    "%token descr"
    # property label
    # property comment

cdef class TokenChars(Token):
    "%token chars"
    # property chars

cdef class TokenWhitespace(Token):
    "%token whitespace"

#
# End Token.
#===============================================================================
# Begin Nonterm.
#

cdef class Nonterm(Parsing.Nonterm):
    cdef readonly int begPos, endPos

cdef class Matrix(Nonterm):
    "%start"
    cpdef reduce(self, RowList RowList)
    #   "%reduce RowList"

cdef class RowList(Nonterm):
    "%nonterm"
    cpdef reduceExtend(self, RowList RowList, Row Row)
    #   "%reduce RowList Row"

    cpdef reduceOne(self, Row Row)
    #   "%reduce Row"

cdef class Row(Nonterm):
    "%nonterm"
    cpdef reduce(self, TokenDescr descr, TokenChars chars)
    #   "%reduce descr chars"

#
# End Nonterm.
#===============================================================================

cdef class Parser(Parsing.Lr):
    cdef readonly Token first, last

    cdef Parsing.Spec _initSpec(self)
    cdef _initReDna(self)
    cdef _initReProtein(self)
    cdef void _appendToken(self, Token token) except *
    cpdef parse(self, str input, type charType=?, int begPos=?, int line=?,
      bint verbose=?)
