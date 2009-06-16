from CxNewickLexer cimport *
cimport Crux.Taxa as Taxa
from Crux.Tree cimport Tree, Node

cdef class Parser:
    cdef int status
    cdef bint bufferedToken
    cdef yyscan_t scanner
    cdef CxtNewickLexerExtra extra
    cdef Tree _tree
    cdef Taxa.Map _taxaMap

    cdef void _labelNode(self, Node node, str label) except *
    cdef void nextToken(self) except *
    cdef getToken(self)
    cdef ungetToken(self)
    cdef str parseLabel(self)
    cdef parseColonBranchLength(self)
    cdef Node parseSubtree(self, double *brlen)
    cdef Node parseSubtreeList(self)
    cdef Node parseDescendantList(self)
    cdef void parseTree(self) except *

    cpdef parse(self, input, int line=*)
