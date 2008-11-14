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

import Crux.Exception

class Exception(Crux.Exception.Exception):
    pass

import exceptions

class SyntaxError(Exception, exceptions.SyntaxError):
    def __init__(self, line, col, str):
        self.line = line
        self.col = col
        self.str = str

    def __str__(self):
        return "%d:%d: %s" % (self.line, self.col, self.str)

class Malformed(Exception, exceptions.SyntaxError):
    def __init__(self, str):
        self._str = str

    def __str__(self):
        return self._str

import re
import sys

cimport Parsing
from Tree cimport Tree, Node, Edge
cimport Taxa

global __name__

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
cdef class NewickTree(Nonterm)
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
    property label:
        def __get__(self):
            cdef str ret

            ret = self._input[self.begPos:self.endPos]
            # Convert '_' to ' '.
            ret = ret.replace('_', ' ')
            return ret
cdef class TokenQuotedLabel(Token):
    "%token quotedLabel"
    property label:
        def __get__(self):
            cdef str ret

            # Strip the enclosing '...'.
            ret = self._input[self.begPos + 1:self.endPos - 1]
            # Convert '' to '.
            ret = ret.replace("''", "'")
            return ret
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

cdef class NewickTree(Nonterm):
    "%start Tree"

    def __init__(self, Parsing.Lr parser):
        Nonterm.__init__(self, parser)

        self.root = Node((<Parser>self.parser)._tree)

    # (A,B)C:4.2;
    #
    #   .
    #   |4.2
    #   C
    #  / \
    # A   B
    cpdef reduceDRB(self, DescendantList DescendantList, Label Label,
      TokenColon colon, TokenBranchLength branchLength,
      TokenSemicolon semicolon):
        "%reduce DescendantList Label colon branchLength semicolon"
        cdef Tree tree
        cdef Edge edge

        tree = (<Parser>self.parser)._tree
        (<Parser>self.parser)._labelNode(DescendantList.node, Label)
        edge = Edge(tree)
        edge.length = float(branchLength.raw)
        edge.attach(self.root, DescendantList.node)
        tree.base = self.root

    # (A,B)C;
    #
    #   .
    #   |
    #   C
    #  / \
    # A   B
    cpdef reduceDR(self, DescendantList DescendantList, Label Label,
      TokenSemicolon semicolon):
        "%reduce DescendantList Label semicolon"
        cdef Tree tree
        cdef Edge edge

        tree = (<Parser>self.parser)._tree
        (<Parser>self.parser)._labelNode(DescendantList.node, Label)
        edge = Edge(tree)
        edge.attach(self.root, DescendantList.node)
        tree.base = self.root

    # (A,B):4.2;
    #
    #   .
    #   |4.2
    #   .
    #  / \
    # A   B
    cpdef reduceDB(self, DescendantList DescendantList, TokenColon colon,
      TokenBranchLength branchLength, TokenSemicolon semicolon):
        "%reduce DescendantList colon branchLength semicolon"
        cdef Tree tree
        cdef Edge edge

        tree = (<Parser>self.parser)._tree
        edge = Edge(tree)
        edge.length = float(branchLength.raw)
        edge.attach(self.root, DescendantList.node)
        tree.base = self.root

    # A:4.2;
    #
    #   .
    #   |4.2
    #   A
    cpdef reduceRB(self, Label Label, TokenColon colon,
      TokenBranchLength branchLength, TokenSemicolon semicolon):
        "%reduce Label colon branchLength semicolon"
        cdef Tree tree
        cdef Node node
        cdef Edge edge

        tree = (<Parser>self.parser)._tree
        node = Node(tree)
        (<Parser>self.parser)._labelNode(node, Label)
        edge = Edge(tree)
        edge.length = float(branchLength.raw)
        edge.attach(self.root, node)
        tree.base = self.root

    # (A,B);
    #
    #   .
    #   |
    #   .
    #  / \
    # A   B
    cpdef reduceD(self, DescendantList DescendantList,
      TokenSemicolon semicolon):
        "%reduce DescendantList semicolon"
        cdef Tree tree
        cdef Edge edge

        tree = (<Parser>self.parser)._tree
        edge = Edge(tree)
        edge.attach(self.root, DescendantList.node)
        tree.base = self.root

    # A;
    #
    #   .
    #   |
    #   A
    cpdef reduceR(self, Label Label, TokenSemicolon semicolon):
        "%reduce Label semicolon"
        cdef Tree tree
        cdef Node node
        cdef Edge edge

        tree = (<Parser>self.parser)._tree
        node = Node(tree)
        (<Parser>self.parser)._labelNode(node, Label)
        edge = Edge(tree)
        edge.attach(self.root, node)
        tree.base = self.root

    # :4.2;
    #
    #   .
    #   |4.2
    #   .
    cpdef reduceB(self, TokenColon colon, TokenBranchLength branchLength,
      TokenSemicolon semicolon):
        "%reduce colon branchLength semicolon"
        cdef Tree tree
        cdef Node node
        cdef Edge edge

        tree = (<Parser>self.parser)._tree
        node = Node(tree)
        edge = Edge(tree)
        edge.length = float(branchLength.raw)
        edge.attach(self.root, node)
        tree.base = self.root

    # ;
    #
    #  .
    cpdef reduce(self, TokenSemicolon semicolon):
        "%reduce semicolon"
        cdef Tree tree

        tree = (<Parser>self.parser)._tree
        tree.base = self.root

cdef class DescendantList(Nonterm):
    "%nonterm"

    cpdef reduce(self, TokenLparen lparen, SubtreeList SubtreeList,
      TokenRparen rparen):
        "%reduce lparen SubtreeList rparen"
        cdef Subtree subtree
        cdef Edge edge

        self.node = Node((<Parser>self.parser)._tree)
        subtree = SubtreeList.last
        while subtree is not None:
            edge = Edge((<Parser>self.parser)._tree)
            edge.length = subtree.len
            edge.attach(self.node, subtree.node)
            subtree = subtree.prev

cdef class SubtreeList(Nonterm):
    "%nonterm"

    cpdef reduceOne(self, Subtree Subtree):
        "%reduce Subtree"
        Subtree.prev = None
        self.last = Subtree

    cpdef reduceExtend(self, SubtreeList SubtreeList, TokenComma comma,
      Subtree Subtree):
        "%reduce SubtreeList comma Subtree"
        Subtree.prev = SubtreeList.last
        self.last = Subtree

cdef class Subtree(Nonterm):
    "%nonterm"

    cpdef reduceDIB(self, DescendantList DescendantList,
      Label Label, TokenColon colon,
      TokenBranchLength branchLength):
        "%reduce DescendantList Label colon branchLength"
        self.node = DescendantList.node
        (<Parser>self.parser)._labelNode(self.node, Label)
        self.len = float(branchLength.raw)

    cpdef reduceDI(self, DescendantList DescendantList, Label Label):
        "%reduce DescendantList Label"
        self.node = DescendantList.node
        (<Parser>self.parser)._labelNode(self.node, Label)
        self.len = 0.0

    cpdef reduceDB(self, DescendantList DescendantList,
      TokenBranchLength branchLength):
        "%reduce DescendantList branchLength [pSubtree]"
        self.node = DescendantList.node
        self.len = float(branchLength.raw)

    cpdef reduceLB(self, Label Label, TokenColon colon,
      TokenBranchLength branchLength):
        "%reduce Label colon branchLength"
        self.node = Node((<Parser>self.parser)._tree)
        (<Parser>self.parser)._labelNode(self.node, Label)
        self.len = float(branchLength.raw)

    cpdef reduceL(self, Label Label):
        "%reduce Label"
        self.node = Node((<Parser>self.parser)._tree)
        (<Parser>self.parser)._labelNode(self.node, Label)
        self.len = 0.0

cdef class Label(Nonterm):
    "%nonterm"

    cpdef reduceU(self, TokenUnquotedLabel unquotedLabel):
        "%reduce unquotedLabel"
        self.label = unquotedLabel.label

    cpdef reduceQ(self, TokenQuotedLabel quotedLabel):
        "%reduce quotedLabel"
        self.label = quotedLabel.label

    cpdef reduceB(self, TokenBranchLength branchLength):
        "%reduce branchLength [pLabel]"
        self.label = branchLength.raw

    cpdef reduceE(self):
        "%reduce [pLabel]"
        self.label = None

#
# End Nonterm.
#===============================================================================

# Regex used to recognize all tokens except complex comments.
cdef _reMain
# Regex used to process complex comments.
cdef _reComment

cdef Parsing.Spec _spec

cdef class Parser(Parsing.Lr):
    def __init__(self, Tree tree, Taxa.Map taxaMap=None,
      Parsing.Spec spec=None):
        global _reMain, _reComment

        if spec is None:
            spec = self._initSpec()

        Parsing.Lr.__init__(self, spec)
        self.first = None
        self.last = None
        self._tree = tree
        self._taxaMap = taxaMap

        if _reMain is None:
            _reMain = self._initReMain()
        if _reComment is None:
            _reComment = self._initReComment()

    # Check whether spec has been initialized here, rather than doing
    # initialization during module initialization.  Lazy initialization is
    # necessary because Crux.Newick isn't defined until the Newick submodule
    # has been fully initialized.
    #
    # Lazy initialization here has the serendipitous advantage of reducing
    # Crux startup time and memory usage.  Indeed, under normal
    # circumstances, Crux.Newick._spec is never even used, since the useful
    # Newick parsers are subclasses of this one.
    cdef Parsing.Spec _initSpec(self):
        global _spec

        if _spec is None:
            _spec = Parsing.Spec(sys.modules[__name__],
              pickleFile="%s/share/Crux-%s/Newick.pickle" %
              (Crux.Config.prefix, Crux.Config.version),
              verbose=(False if (not __debug__ or Crux.Config.quiet) else True),
              skinny=(False if __debug__ else True),
              logFile="%s/share/Crux-%s/Newick.log" %
              (Crux.Config.prefix, Crux.Config.version))
        return _spec

    cdef _initReMain(self):
        return re.compile(r"""
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
     (?:[eE][-+]?[0-9]+)?
     (?!_))                # branch length
  | ([^ \t\r\n()[\]':;,]+) # unquoted label
  | ('(?:''|[^'])*'(?!'))  # quoted label
  | ([ \t\r\f\v]+)         # whitespace
  | ([\n])                 # whitespace (newline)
""", re.X)

    cdef _initReComment(self):
        return re.compile(r"""
    (\[)         # start comment
  | (\])         # end comment
  | (\n)         # newline
  | ([^\[\]\n]+) # text
""", re.X)

    cdef void _appendToken(self, Token token) except *:
        if self.first is None:
            self.first = token
            self.last = token
        else:
            self.last.next = token
            token.prev = self.last
            self.last = token

    cdef void _labelNode(self, Node node, Label label) except *:
        cdef int ind

        if label.label is None:
            return

        if self._taxaMap is not None:
            try:
                taxon = self._ind2taxon(int(label.label))
            except:
                raise Malformed("No Taxa.Map index entry for %r" % label.label)
        else:
            taxon = Taxa.get(label.label)

        node.taxon = taxon

    cdef str expandInput(self, str input, int pos, int line, int col):
        """
Called when end of input is reached.  By default a SyntaxError is raised.
"""
        raise SyntaxError(line, col, "Invalid token or end of input reached")

    cpdef parse(self, str input, int begPos=0, int line=1, int col=0,
      bint verbose=False):
        cdef int pos = begPos
        cdef object m
        cdef int idx, start, end
        cdef Token token
        cdef int tokLine, tokCol, nesting, tokPos

        self.verbose = verbose

        # Iteratively tokenize the input and feed the tokens to the parser.
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
                token = TokenComment(self, input, start, end, tokLine, tokCol)
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
                token = TokenComment(self, input, tokPos, end, tokLine,
                  tokCol)
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
                token = TokenSemicolon(self, input, start, end, tokLine,
                  tokCol)
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
