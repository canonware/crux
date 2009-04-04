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
    def __init__(self, line, str):
        self.line = line
        self.str = str

    def __str__(self):
        return "Line %d: %s" % (self.line, self.str)

class Malformed(Exception, exceptions.SyntaxError):
    def __init__(self, str):
        self._str = str

    def __str__(self):
        return self._str

import re
import sys

cimport Parsing
from CxNewickLexer cimport *
from Crux.Tree cimport Tree, Node, Edge
cimport Crux.Taxa as Taxa

cdef extern from "Python.h":
    cdef object PyString_FromStringAndSize(char *s, Py_ssize_t len)

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
    def __init__(self, Parser parser, str input, int line):
        Parsing.Token.__init__(self, parser)

        self.prev = None
        self.next = None

        self.raw = input
        self.line = line

        parser._appendToken(self)

    def __repr__(self):
        return Parsing.Token.__repr__(self) + (" (%r)" % self.raw)

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
            # Convert '_' to ' '.
            return self.raw.replace('_', ' ')
cdef class TokenQuotedLabel(Token):
    "%token quotedLabel"
    property label:
        def __get__(self):
            cdef str ret

            # Strip the enclosing '...'.
            ret = self.raw[1:-1]
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

cdef Parsing.Spec _spec

cdef class Parser(Parsing.Lr):
    def __init__(self, Tree tree, Taxa.Map taxaMap=None,
      Parsing.Spec spec=None):
        if spec is None:
            spec = self._initSpec()

        Parsing.Lr.__init__(self, spec)
        self.first = None
        self.last = None
        self._tree = tree
        self._taxaMap = taxaMap

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
              pickleFile="%s/Crux/parsers/Newick.pickle" % Crux.Config.datadir,
              verbose=(True if (__debug__ and Crux.Config.verbose) else False),
              skinny=(False if __debug__ else True),
              logFile="%s/Crux/parsers/Newick.log" % Crux.Config.datadir)
        return _spec

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

    cpdef parse(self, input, int line=1, bint verbose=False):
        cdef yyscan_t scanner
        cdef CxtNewickLexerExtra extra
        cdef int status, tokType

        assert type(input) in (file, str)

        self.verbose = verbose

        extra.line = line
        extra.ioerror = 0
        if type(input) == file:
            extra.inputMode = CxeNewickLexerInputModeFd
            extra.input.fd = input.fileno()
        elif type(input) == str:
            extra.inputMode = CxeNewickLexerInputModeStr
            extra.input.s.s = input
            extra.input.s.len = len(input)
        else:
            assert False

        if CxNewickLexer_lex_init_extra(&extra, &scanner) != 0:
            raise MemoryError("Error initializing lexer")

        try:
            # Iteratively tokenize the input and feed the tokens to the parser.
            status = CxNewickLexer_lex(scanner)
            while status > 0:
                tokType = extra.tokType
                if tokType == CxeNewickLexerTokTypeComment: # [...]
                    token = TokenComment(self, \
                      PyString_FromStringAndSize(extra.tok.s, extra.tok.len), \
                      extra.line)
                elif tokType == CxeNewickLexerTokTypeLparen: # (
                    token = TokenLparen(self, "(", extra.line)
                    self.token(token)
                elif tokType == CxeNewickLexerTokTypeRparen: # )
                    token = TokenRparen(self, ")", extra.line)
                    self.token(token)
                elif tokType == CxeNewickLexerTokTypeComma: # ,
                    token = TokenComma(self, ",", extra.line)
                    self.token(token)
                elif tokType == CxeNewickLexerTokTypeColon: # :
                    token = TokenColon(self, ":", extra.line)
                    self.token(token)
                elif tokType == CxeNewickLexerTokTypeSemicolon: # ;
                    token = TokenSemicolon(self, ";", extra.line)
                    self.token(token)

                    # Finish.
                    self.eoi()
                    return

                elif tokType == CxeNewickLexerTokTypeBranchLength:
                    token = TokenBranchLength(self, \
                      PyString_FromStringAndSize(extra.tok.s, extra.tok.len), \
                      extra.line)
                    self.token(token)
                elif tokType == CxeNewickLexerTokTypeUnquotedLabel:
                    token = TokenUnquotedLabel(self, \
                      PyString_FromStringAndSize(extra.tok.s, extra.tok.len), \
                      extra.line)
                    self.token(token)
                elif tokType == CxeNewickLexerTokTypeQuotedLabel:
                    token = TokenQuotedLabel(self, \
                      PyString_FromStringAndSize(extra.tok.s, extra.tok.len), \
                      extra.line)
                    self.token(token)
                elif tokType == CxeNewickLexerTokTypeWhitespace: # whitespace
                    token = TokenWhitespace(self, \
                      PyString_FromStringAndSize(extra.tok.s, extra.tok.len), \
                      extra.line)
                else:
                    assert False

                status = CxNewickLexer_lex(scanner)

            if status == -1:
                raise SyntaxError(extra.line, "Invalid token")
            raise SyntaxError(extra.line, "End of input reached")
        finally:
            CxNewickLexer_lex_destroy(scanner);
