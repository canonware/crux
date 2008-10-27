#===============================================================================
# Newick Tree parser.  This implementation expands from Gary Olsen's
# interpretation of the Newick Tree Format Standard, available at:
#
#   http://evolution.genetics.washington.edu/phylip/newick_doc.html
#
# Joe Felsenstein's description of the Newick Tree Format is at:
#
#   http://evolution.genetics.washington.edu/phylip/newicktree.html
#
# The two descriptions directly conflict in their definitions of the top level
# production.  Felsenstein claims that the following is a valid tree:
#
#   A;
#
# However, Olsen claims that the descendant_list component of the tree
# production is mandatory, and that the descendant_list production must not be
# empty.  This implementation follows Felsenstein's interpretation, since it is
# a superset of Olsen's interpretation.
#
# Similarly, Felsenstein claims that zero length labels are allowable, whereas
# Olsen claims that a label is a "string of printing characters".  This
# implementation allows "zero length strings of printing characters".
#===============================================================================

# Forward declarations.
cdef class Token
cdef class Nonterm
cdef class Parser

cimport Parsing

import Crux.Exception

class Exception(Crux.Exception.Exception):
    pass

import exceptions

class SyntaxError(Exception, exceptions.SyntaxError):
    def __init__(self, str):
        self._str = str

    def __str__(self):
        return self._str

import re
import sys

global __name__

#===============================================================================
# Begin Precedence.
#

class pColon(Parsing.Precedence):
    "%fail"
class pSemicolon(Parsing.Precedence):
    "%fail"
class pLabel(Parsing.Precedence):
    "%fail <pColon <pSemicolon"

#
# End Precedence.
#===============================================================================
# Begin Token.
#

cdef class Token(Parsing.Token):
    def __init__(self, Parser parser, str input, int begPos, int endPos,
      int line, int col):
        Parsing.Token.__init__(self, parser)

        self.prev = None
        self.next = None

        self._input = input
        self.begPos = begPos
        self.endPos = endPos
        self.line = line
        self.col = col

        parser._appendToken(self)

    def __repr__(self):
        return Parsing.Token.__repr__(self) + (" (%r)" % self.raw)

    property raw:
        def __get__(self):
            return self._input[self.begPos:self.endPos]

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

    property raw:
        def __get__(self):
            # Convert '_' to ' '.
            return self._input[self.begPos:self.endPos].replace('_', ' ')

cdef class TokenQuotedLabel(Token):
    "%token quotedLabel"

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
    def __init__(self, Parsing.Lr parser):
        Parsing.Nonterm.__init__(self, parser)

        self.begPos = -1
        self.endPos = -1
        self.variant = None

cdef class Tree(Nonterm):
    "%start Tree"

    def reduceDRB(self, DescendantList, RootLabel, colon, branchLength,
      semicolon):
        "%reduce DescendantList RootLabel colon branchLength semicolon"

    def reduceDR(self, DescendantList, RootLabel, semicolon):
        "%reduce DescendantList RootLabel semicolon"

    def reduceDB(self, DescendantList, colon, branchLength, semicolon):
        "%reduce DescendantList colon branchLength semicolon"

    def reduceRB(self, RootLabel, colon, branchLength, semicolon):
        "%reduce RootLabel colon branchLength semicolon"

    def reduceD(self, DescendantList, semicolon):
        "%reduce DescendantList semicolon"

    def reduceR(self, RootLabel, semicolon):
        "%reduce RootLabel semicolon"

    def reduceB(self, colon, branchLength, semicolon):
        "%reduce colon branchLength semicolon"

    def reduce(self, semicolon):
        "%reduce semicolon"

cdef class DescendantList(Nonterm):
    "%nonterm"

    def reduce(self, lparen, SubtreeList, rparen):
        "%reduce lparen SubtreeList rparen"

cdef class SubtreeList(Nonterm):
    "%nonterm"

    def reduceOne(self, Subtree):
        "%reduce Subtree"

    def reduceExtend(self, SubtreeList, comma, Subtree):
        "%reduce SubtreeList comma Subtree"

cdef class Subtree(Nonterm):
    "%nonterm"

    def reduceDIB(self, DescendantList, InternalLabel, colon, branchLength):
        "%reduce DescendantList InternalLabel colon branchLength"

    def reduceDI(self, DescendantList, InternalLabel):
        "%reduce DescendantList InternalLabel"

    def reduceDB(self, DescendantList, branchLength):
        "%reduce DescendantList branchLength"

    def reduceLB(self, LeafLabel, colon, branchLength):
        "%reduce LeafLabel colon branchLength"

    def reduceL(self, LeafLabel):
        "%reduce LeafLabel"

cdef class RootLabel(Nonterm):
    "%nonterm"

    def reduce(self, Label):
        "%reduce Label"

cdef class InternalLabel(Nonterm):
    "%nonterm"

    def reduce(self, Label):
        "%reduce Label"

cdef class LeafLabel(Nonterm):
    "%nonterm"

    def reduce(self, Label):
        "%reduce Label"

cdef class Label(Nonterm):
    "%nonterm"

    def reduceU(self, unquotedLabel):
        "%reduce unquotedLabel"

    def reduceQ(self, quotedLabel):
        "%reduce quotedLabel"

    def reduceE(self):
        "%reduce [pLabel]"

#
# End Nonterm.
#===============================================================================

cdef Parsing.Spec spec
spec = Parsing.Spec(sys.modules[__name__], "Newick.pickle", \
  verbose=True, skinny=False, logFile="Newick.log")

# Regex used to recognize all tokens except complex comments.
cdef _reMain
_reMain = re.compile(r"""
    (\[[^[\]\n]*\])        # simple comment (non-nested, single line)
  | (\[[^[\]\n]*\[)        # complex comment prefix: [...[
  | (\[[^[\]\n]*\n)        # complex comment prefix: [...\n
  | ([(])                  # (
  | ([)])                  # )
  | (,)                    # ,
  | (:)                    # :
  | (;)                    # ;
  | ([-+]?
     [0-9]+(?:[.][0-9]+)?
     (?:[eE][-+]?[0-9]+)?) # branch length
  | ([^ \t\r\n()[\]':;,]+) # unquoted label
  | (''(?!'))              # quoted label
  | ([ \t\r]+)             # whitespace
  | ([\n])                 # whitespace (newline)
""", re.X)

# Regex used to process complex comments.
cdef _reComment
_reComment = re.compile(r"""
    (\[)         # start comment
  | (\])         # end comment
  | (\n)         # newline
  | ([^\[\]\n]+) # text
""", re.X)

cdef class Parser(Parsing.Lr):
    def __init__(self):
        Parsing.Lr.__init__(self, spec)
        self.first = None
        self.last = None

    cdef void _appendToken(self, Token token) except *:
        if self.first is None:
            self.first = token
            self.last = token
        else:
            self.last.next = token
            token.prev = self.last
            self.last = token

    cdef str expandInput(self, str input, int pos, int line, int col):
        """
Called when end of input is reached.  By default a SyntaxError is raised.
"""
        raise SyntaxError("%d:%d: Invalid token or end of input reached" \
          % (line, col))

    cpdef parse(self, str input, int begPos=0, int line=1, int col=0,
      bint verbose=False):
        cdef int pos = begPos
        cdef object m
        cdef int idx, start, end
        cdef Token token
        cdef int nesting, tokPos

        self.verbose = verbose

        while True:
            tokLine = line
            tokCol = col

            m = _reMain.match(input, pos)
            while m is None:
                input = self.expandInput(input, pos, tokLine, tokCol)
                m = _reMain.match(input, pos)
            idx = m.lastindex
            start = m.start(idx)
            end = m.end(idx)
            col += end - start
            if idx == 1:    # simple comment
                token = TokenComment(self, input, start,
                  end, tokLine, tokCol)
            elif idx == 2 or idx == 3:  # complex comment prefix
                nesting = (2 if idx == 2 else 1)
                tokPos = pos
                line += (0 if idx == 2 else 1)
                pos += end - start
                while nesting > 0:
                    m = _reComment.match(input, pos)
                    while m is None:
                        input = self.expandInput(input, tokPos, tokLine, tokCol)
                        m = _reComment.match(input, pos)
                    idx = m.lastindex
                    start = m.start(idx)
                    end = m.end(idx)
                    col += end - start
                    if idx == 1: # start comment
                        nesting += 1
                    elif idx == 2: # end comment
                        nesting -= 1
                    elif idx == 3: # newline
                        line += 1
                        col = 0
                    elif idx == 4: # text
                        pass
                    else:
                        assert False
                    pos += end - start
                token = TokenComment(self, input, tokPos, end, tokLine, tokCol)
            elif idx == 4:  # (
                token = TokenLparen(self, input, start, end, tokLine, tokCol)
                self.token(token)
            elif idx == 5:  # )
                token = TokenRparen(self, input, start, end, tokLine, tokCol)
                self.token(token)
            elif idx == 6:  # ,
                token = TokenComma(self, input, start, end, tokLine, tokCol)
                self.token(token)
            elif idx == 7:  # :
                token = TokenColon(self, input, start, end, tokLine, tokCol)
                self.token(token)
            elif idx == 8:  # ;
                token = TokenSemicolon(self, input, start, end, tokLine, tokCol)
                self.token(token)

                # Finish.
                pos = end
                self.eoi()
                return (pos, line, col)
            elif idx == 9:  # branch length
                token = TokenBranchLength(self, input, start, end, tokLine,
                  tokCol)
                self.token(token)
            elif idx == 10:  # unquoted label
                token = TokenUnquotedLabel(self, input, start, end, tokLine,
                  tokCol)
                self.token(token)
            elif idx == 11: # quoted label
                token = TokenQuotedLabel(self, input, start, end, tokLine,
                  tokCol)
                self.token(token)
            elif idx == 12: # whitespace
                token = TokenWhitespace(self, input, start, end, tokLine,
                  tokCol)
            elif idx == 13: # whitespace (newline)
                token = TokenWhitespace(self, input, start, end, tokLine,
                  tokCol)
                line += 1
                col = 0
            else:
                assert False

            pos = end
