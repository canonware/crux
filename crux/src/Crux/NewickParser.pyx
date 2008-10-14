import sys
cimport Parsing

# XXX The basic idea is to write a Newick parser that calls NewickParser
# virtual methods upon accepting various productions.  This allows a subclass
# to provide methods that build a tree (or do whatever).

#===============================================================================
# Begin Precedence.
#


#
# End Precedence.
#===============================================================================
# Begin Token.
#

cdef class Node(Parsing.Nonterm):
    cdef int begPos, endPos
    cdef object variant

    def __init__(self, parser):
        Parsing.Nonterm.__init__(self, parser)

        self.begPos = -1
        self.endPos = -1
        self.variant = None

class Token(Parsing.Token, Node):
    pass # XXX

#
# End Token.
#===============================================================================
# Begin Nonterm.
#

#
# End Nonterm.
#===============================================================================

spec = Parsing.Spec(sys.modules[__name__], "NewickParser.pickle", \
  verbose=True, skinny=False, logFile="NewickParser.log")

cdef class NewickParser(Parsing.Lr):
    cdef parse(self, input):
        pass # XXX

    cdef acceptXXX(self, XXX):
        # Virtual method.
        pass
