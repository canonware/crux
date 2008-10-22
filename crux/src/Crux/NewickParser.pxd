cimport Parsing

cdef class Node(Parsing.Nonterm):
    cdef int begPos, endPos
    cdef object variant

cdef class NewickParser(Parsing.Lr):
    cdef parse(self, input)
    cdef acceptXXX(self, XXX)
