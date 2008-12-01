from Crux.CTMatrix cimport CTMatrix
cimport Parsing
from Crux.Taxa cimport Taxon
cimport Crux.Taxa as Taxa

# Forward declarations.
cdef class Token(Parsing.Token)
cdef class TokenDescr(Token)
cdef class TokenChars(Token)
cdef class Nonterm(Parsing.Nonterm)
cdef class Matrix(Nonterm)
cdef class RowList(Nonterm)
cdef class Row(Nonterm)
cdef class Chars(Nonterm)
cdef class Parser(Parsing.Lr)

#===============================================================================
# Begin Token.
#

cdef class Token(Parsing.Token):
    cdef str raw

cdef class TokenDescr(Token):
    "%token descr"
    cdef str label
    cdef str comment

cdef class TokenChars(Token):
    "%token chars"
    cdef str chars

#
# End Token.
#===============================================================================
# Begin Nonterm.
#

cdef class Nonterm(Parsing.Nonterm): pass

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
    cpdef reduce(self, TokenDescr descr, Chars Chars)
    #   "%reduce descr chars"

cdef class Chars(Nonterm):
    "%nonterm"
    cdef list chars

    cpdef reduceOne(self, TokenChars chars)
    #   "%reduce chars"

    cpdef reduceExtend(self, Chars Chars, TokenChars chars)
    #   "%reduce Chars chars"

#
# End Nonterm.
#===============================================================================

cdef class Parser(Parsing.Lr):
    cdef CTMatrix matrix
    cdef Taxa.Map taxaMap

    cdef Parsing.Spec _initSpec(self)
    cdef _initReDna(self)
    cdef _initReProtein(self)
    cdef void _addTaxon(self, Taxon taxon, str chars) except *
    cpdef parse(self, lines, type charType=*, int line=*, bint verbose=*)
