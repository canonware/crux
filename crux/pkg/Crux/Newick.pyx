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

import exceptions

class SyntaxError(exceptions.SyntaxError):
    def __init__(self, line, str):
        self.line = line
        self.str = str

    def __str__(self):
        return "Line %d: %s" % (self.line, self.str)

import re
import sys
import Crux.Config

from libm cimport *
from CxNewickLexer cimport *
from Crux.Tree cimport Tree, Node, Edge
cimport Crux.Taxa as Taxa

cdef extern from "Python.h":
    cdef object PyString_FromStringAndSize(char *s, Py_ssize_t len)

cdef class Parser:
    """
        Newick tree parser.
    """
    def __init__(self, Tree tree, Taxa.Map taxaMap=None):
        self._tree = tree
        self._taxaMap = taxaMap

    cdef void _labelNode(self, Node node, str label) except *:
        cdef int ind

        if self._taxaMap is not None:
            try:
                taxon = self._ind2taxon(int(label))
            except:
                raise exceptions.SyntaxError("No Taxa.Map index entry for %r" \
                  % label)
        else:
            taxon = Taxa.get(label)

        node.taxon = taxon

    cdef void nextToken(self) except *:
        cdef int tokType

        if self.status == 0:
            raise SyntaxError(self.extra.line, "End of input reached")
        elif self.status < 0:
            raise SyntaxError(self.extra.line, "Invalid token")

        while True:
            self.status = CxNewickLexer_lex(self.scanner)
            if self.status == 0:
                raise SyntaxError(self.extra.line, "End of input reached")
            elif self.status < 0:
                raise SyntaxError(self.extra.line, "Invalid token")
            tokType = self.extra.tokType
            if tokType == CxeNewickLexerTokTypeComment: # [...]
                # Swallow comments.
                pass
            elif tokType == CxeNewickLexerTokTypeWhitespace: # whitespace
                # Swallow whitespace.
                pass
            elif tokType == CxeNewickLexerTokTypeSemicolon: # ;
                # Make note that no more tokens should be parsed.
                self.status = -1
                return
            else:
                # self.extra contains info about the token.
                return

    cdef getToken(self):
        # Use buffered token if it exists.
        if self.bufferedToken:
            self.bufferedToken = False
        else:
            self.nextToken()
        s = PyString_FromStringAndSize(self.extra.tok.s, self.extra.tok.len) # XXX Remove.
#        print "Get: %s\t%r" % \
#          (["comment", "(", ")", ",", ":", ";", "brlen", "ulabel", "qlabel",
#          "whitespace"][self.extra.tokType],
#          [s, "(", ")", ",", ":", ";", s, s, s, s][self.extra.tokType])

    cdef ungetToken(self):
        assert not self.bufferedToken
        self.bufferedToken = True
#        print "Unget"

    cdef str parseLabel(self):
        cdef int tokType
        cdef str s

        self.getToken()
        tokType = self.extra.tokType
        if tokType == CxeNewickLexerTokTypeUnquotedLabel:
            s = PyString_FromStringAndSize(self.extra.tok.s, self.extra.tok.len)
            # Convert '_' to ' '.
            s = s.replace('_', ' ')
            return s
        elif tokType == CxeNewickLexerTokTypeQuotedLabel:
            s = PyString_FromStringAndSize(self.extra.tok.s, self.extra.tok.len)
            # Strip the enclosing '...'.
            s = s[1:-1]
            # Convert '' to '.
            s = s.replace("''", "'")
            return s
        elif tokType == CxeNewickLexerTokTypeBranchLength:
            s = PyString_FromStringAndSize(self.extra.tok.s, self.extra.tok.len)
            return s
        else:
            self.ungetToken()
            return None

    cdef parseColonBranchLength(self):
        cdef str s
        cdef double brlen

        self.getToken()
        if self.extra.tokType == CxeNewickLexerTokTypeColon:
            self.getToken()
            s = PyString_FromStringAndSize(self.extra.tok.s, self.extra.tok.len)
            if self.extra.tokType != CxeNewickLexerTokTypeBranchLength:
                raise SyntaxError(self.extra.line, "Unexpected token: %r" % s)
            brlen = float(s)
        else:
            self.ungetToken()
            brlen = NAN
        return brlen

    cdef Node parseSubtree(self, double *brlen):
        cdef Node node, descendantList

        descendantList = self.parseDescendantList()
        label = self.parseLabel()
        brlen[0] = self.parseColonBranchLength()
        if isnan(brlen[0]):
            brlen[0] = 0.0

        if descendantList is not None:
            node = descendantList
        else:
            node = Node(self._tree)

        if label is not None:
            self._labelNode(node, label)

        return node

    cdef Node parseSubtreeList(self):
        cdef Node ret, subtree
        cdef double brlen

        ret = Node(self._tree)
        while True:
            subtree = self.parseSubtree(&brlen)
            edge = Edge(self._tree)
            if not isnan(brlen):
                edge.length = brlen
            edge.attach(ret, subtree)
            ret.ring = ret.ring.next # Re-order ring.

            self.getToken()
            if self.extra.tokType != CxeNewickLexerTokTypeComma:
                self.ungetToken()
                return ret

    cdef Node parseDescendantList(self):
        cdef Node ret
        cdef list subtreeList
        cdef int tokType
        cdef str s

        self.getToken()
        if self.extra.tokType != CxeNewickLexerTokTypeLparen:
            self.ungetToken()
            return None
        ret = self.parseSubtreeList()
        self.getToken()
        if self.extra.tokType != CxeNewickLexerTokTypeRparen:
            s = PyString_FromStringAndSize(self.extra.tok.s, self.extra.tok.len)
            raise SyntaxError(self.extra.line, "Unexpected token: %r" % s)
        return ret

    cdef void parseTree(self) except *:
        cdef Node descendantList, root
        cdef str label, s
        cdef double brlen
        cdef Tree tree
        cdef Edge edge

        descendantList = self.parseDescendantList()
        label = self.parseLabel()
        brlen = self.parseColonBranchLength()

        self.getToken()
        if self.extra.tokType != CxeNewickLexerTokTypeSemicolon:
            s = PyString_FromStringAndSize(self.extra.tok.s, self.extra.tok.len)
            raise SyntaxError(self.extra.line, "Unexpected token: %r" % s)

        tree = self._tree
        root = Node(tree)
        self._tree.base = root

        if descendantList is not None:
            node = descendantList
            if label is not None:
                self._labelNode(node, label)
        elif label is not None:
            node = Node(tree)
            self._labelNode(node, label)
        elif not isnan(brlen):
            node = Node(tree)
        else:
            return

        edge = Edge(tree)
        if not isnan(brlen):
            edge.length = brlen
        edge.attach(root, node)

    cpdef parse(self, input, int line=1):
        """
            Parse Newick input from a file or string and store the results in
            the Tree that was passed to the constructor.
        """
        assert type(input) in (file, str)

        self.status = 1
        self.bufferedToken = False

        self.extra.line = line
        self.extra.ioerror = 0
        if type(input) == file:
            self.extra.inputMode = CxeNewickLexerInputModeFd
            self.extra.input.fd = input.fileno()
        elif type(input) == str:
            self.extra.inputMode = CxeNewickLexerInputModeStr
            self.extra.input.s.s = input
            self.extra.input.s.len = len(input)
        else:
            assert False

        if CxNewickLexer_lex_init_extra(&self.extra, &self.scanner) != 0:
            raise MemoryError("Error initializing lexer")

        try:
            # Use recursive descent parsing to parse the input.
            self.parseTree()
        finally:
            CxNewickLexer_lex_destroy(self.scanner);
