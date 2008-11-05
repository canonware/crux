cimport Parsing

# Forward declarations.
cdef class pColon(Parsing.Precedence)
cdef class pSemicolon(Parsing.Precedence)
cdef class pSubtree(Parsing.Precedence)
cdef class pLabel(Parsing.Precedence)

cdef class Token(Parsing.Token)
cdef class TokenLparen(Token)
cdef class TokenRparen(Token)
cdef class TokenComma(Token)
cdef class TokenColon(Token)
cdef class TokenSemicolon(Token)
cdef class TokenBranchLength(Token)
cdef class TokenUnquotedLabel(Token)
cdef class TokenQuotedLabel(Token)
cdef class TokenComment(Token)
cdef class TokenWhitespace(Token)

cdef class Nonterm(Parsing.Nonterm)
cdef class Tree(Nonterm)
cdef class DescendantList(Nonterm)
cdef class SubtreeList(Nonterm)
cdef class Subtree(Nonterm)
cdef class Label(Nonterm)

cdef class Parser(Parsing.Lr)

#===============================================================================
# Begin Precedence.
#
cdef class pColon(Parsing.Precedence):
    "%fail"
cdef class pSemicolon(Parsing.Precedence):
    "%fail"
cdef class pSubtree(Parsing.Precedence):
    "%fail"
cdef class pLabel(Parsing.Precedence):
    "%fail <pColon <pSemicolon <pSubtree"

#
# End Precedence.
#===============================================================================
# Begin Token.
#

cdef class Token(Parsing.Token):
    cdef readonly Token prev, next
    cdef str _input
    cdef readonly int begPos, endPos
    cdef readonly int line, col
    # property raw

cdef class TokenLparen(Token):
    "%token lparen"
cdef class TokenRparen(Token):
    "%token rparen"
cdef class TokenComma(Token):
    "%token comma"
cdef class TokenColon(Token):
    "%token colon [pColon]"
cdef class TokenSemicolon(Token):
    "%token semicolon [pSemicolon]"
cdef class TokenBranchLength(Token):
    "%token branchLength"
cdef class TokenUnquotedLabel(Token):
    "%token unquotedLabel"
    # property label
cdef class TokenQuotedLabel(Token):
    "%token quotedLabel"
    # property label
cdef class TokenComment(Token):
    "%token comment"
cdef class TokenWhitespace(Token):
    "%token whitespace"

#
# End Token.
#===============================================================================
# Begin Nonterm.
#

cdef class Nonterm(Parsing.Nonterm):
    cdef readonly int begPos, endPos

cdef class Tree(Nonterm):
    "%start Tree"

    cpdef reduceDRB(self, DescendantList DescendantList, Label Label,
      TokenColon colon, TokenBranchLength branchLength,
      TokenSemicolon semicolon)
    #   "%reduce DescendantList Label colon branchLength semicolon"

    cpdef reduceDR(self, DescendantList DescendantList, Label Label,
      TokenSemicolon semicolon)
    #   "%reduce DescendantList Label semicolon"

    cpdef reduceDB(self, DescendantList DescendantList, TokenColon colon,
      TokenBranchLength branchLength, TokenSemicolon semicolon)
    #   "%reduce DescendantList colon branchLength semicolon"

    cpdef reduceRB(self, Label Label, TokenColon colon,
      TokenBranchLength branchLength, TokenSemicolon semicolon)
    #   "%reduce Label colon branchLength semicolon"

    cpdef reduceD(self, DescendantList DescendantList, TokenSemicolon semicolon)
    #   "%reduce DescendantList semicolon"

    cpdef reduceR(self, Label Label, TokenSemicolon semicolon)
    #   "%reduce Label semicolon"

    cpdef reduceB(self, TokenColon colon, TokenBranchLength branchLength,
      TokenSemicolon semicolon)
    #   "%reduce colon branchLength semicolon"

    cpdef reduce(self, TokenSemicolon semicolon)
    #   "%reduce semicolon"

cdef class DescendantList(Nonterm):
    "%nonterm"

    cpdef reduce(self, TokenLparen lparen, SubtreeList SubtreeList,
      TokenRparen rparen)
    #   "%reduce lparen SubtreeList rparen"

cdef class SubtreeList(Nonterm):
    "%nonterm"

    cpdef reduceOne(self, Subtree Subtree)
    #   "%reduce Subtree"

    cpdef reduceExtend(self, SubtreeList SubtreeList, TokenComma comma,
      Subtree Subtree)
    #   "%reduce SubtreeList comma Subtree"

cdef class Subtree(Nonterm):
    "%nonterm"

    cpdef reduceDIB(self, DescendantList DescendantList,
      Label Label, TokenColon colon,
      TokenBranchLength branchLength)
    #   "%reduce DescendantList Label colon branchLength"

    cpdef reduceDI(self, DescendantList DescendantList,
      Label Label)
    #   "%reduce DescendantList Label"

    cpdef reduceDB(self, DescendantList DescendantList,
      TokenBranchLength branchLength)
    #   "%reduce DescendantList branchLength"

    cpdef reduceLB(self, Label Label, TokenColon colon,
      TokenBranchLength branchLength)
    #   "%reduce Label colon branchLength"

    cpdef reduceL(self, Label Label)
    #   "%reduce Label"

cdef class Label(Nonterm):
    "%nonterm"

    cpdef reduceU(self, TokenUnquotedLabel unquotedLabel)
    #   "%reduce unquotedLabel"

    cpdef reduceQ(self, TokenQuotedLabel quotedLabel)
    #   "%reduce quotedLabel"

    cpdef reduceB(self, TokenBranchLength branchLength)
    #   "%reduce branchLength"

    cpdef reduceE(self)
    #   "%reduce [pLabel]"

#
# End Nonterm.
#===============================================================================

cdef class Parser(Parsing.Lr):
    cdef readonly Token first, last

    cdef Parsing.Spec _initSpec(self)
    cdef _initReMain(self)
    cdef _initReComment(self)
    cdef void _appendToken(self, Token token) except *
    cdef str expandInput(self, str input, int pos, int line, int col)

    cdef Parsing.Token newTokenComment(self, str input, int start, int end,
      int tokLine, int tokCol)
    cdef Parsing.Token newTokenLparen(self, str input, int start, int end,
      int tokLine, int tokCol)
    cdef Parsing.Token newTokenRparen(self, str input, int start, int end,
      int tokLine, int tokCol)
    cdef Parsing.Token newTokenComma(self, str input, int start, int end,
      int tokLine, int tokCol)
    cdef Parsing.Token newTokenColon(self, str input, int start, int end,
      int tokLine, int tokCol)
    cdef Parsing.Token newTokenSemicolon(self, str input, int start, int end,
      int tokLine, int tokCol)
    cdef Parsing.Token newTokenBranchLength(self, str input, int start, int end,
      int tokLine, int tokCol)
    cdef Parsing.Token newTokenUnquotedLabel(self, str input, int start,
      int end, int tokLine, int tokCol)
    cdef Parsing.Token newTokenQuotedLabel(self, str input, int start, int end,
      int tokLine, int tokCol)
    cdef Parsing.Token newTokenWhitespace(self, str input, int start, int end,
      int tokLine, int tokCol)

    cpdef parse(self, str input, int begPos=?, int line=?, int col=?, \
      bint verbose=?)
