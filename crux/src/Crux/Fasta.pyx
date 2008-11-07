################################################################################
#
# FASTA parser.  This parser can parse a superset of the standard FASTA format.
# Ordinarily, a FASTA file looks something like:
#
#   >Taxon_A
#   ACGTACGT--ACGT
#   >Taxon_B
#   ACGTACGTTT--GT
#
# Taxon labels ordinarily should not exceed 30 characters, nor should they have
# any whitespace in them, but this parser does not enforce these limitations.
#
# The character data can be either DNA or protein sequences, and need not be
# aligned.  Due to the large overlap of legal letters in protein and DNA
# sequences, the parser must be manually told whether the input data are protein
# or DNA.
#
################################################################################

cimport Parsing
from Character cimport Dna, Protein

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

import re
import sys

global __name__

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
    def __init__(self, Parser parser, str input, int begPos, int endPos,
      int line):
        Parsing.Token.__init__(self, parser)

        self.prev = None
        self.next = None

        self._input = input
        self.begPos = begPos
        self.endPos = endPos
        self.line = line

        parser._appendToken(self)

    def __repr__(self):
        return Parsing.Token.__repr__(self) + (" (%r)" % self.raw)

    property raw:
        def __get__(self):
            return self._input[self.begPos:self.endPos]

cdef class TokenDescr(Token):
    "%token descr"

    property label:
        def __get__(self):
            cdef str ret

            # Strip the leading '>' and discard comments.
            ret = re.match(r">(\S+)", self.raw).group(1)
            # Convert '_' to ' '.
            ret = ret.replace('_', ' ')
            return ret

    property comment:
        def __get__(self):
            return re.match(r">\S+\s+(.*)$", self.raw).group(1)

cdef class TokenChars(Token):
    "%token chars"

    property chars:
        def __get__(self):
            cdef str ret

            # Discard whitespace.
            ret = re.sub(r'[ \t\r\n]+', r'', self.raw)
            return ret

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

cdef class Matrix(Nonterm):
    "%start"
    cpdef reduce(self, RowList RowList):
        "%reduce RowList"

cdef class RowList(Nonterm):
    "%nonterm"
    cpdef reduceExtend(self, RowList RowList, Row Row):
        "%reduce RowList Row"

    cpdef reduceOne(self, Row Row):
        "%reduce Row"

cdef class Row(Nonterm):
    "%nonterm"
    cpdef reduce(self, TokenDescr descr, TokenChars chars):
        "%reduce descr chars"

#
# End Nonterm.
#===============================================================================

# Regexs used to recognize all tokens.
cdef _reDna
cdef _reProtein

cdef Parsing.Spec _spec

cdef class Parser(Parsing.Lr):
    def __init__(self, Parsing.Spec spec=None):
        global _spec

        if spec is None:
            if _spec is None:
                _spec = self._initSpec()
            spec = _spec
        Parsing.Lr.__init__(self, spec)
        self.first = None
        self.last = None

    cdef Parsing.Spec _initSpec(self):
        return Parsing.Spec([sys.modules[__name__]],
          pickleFile="%s/share/Crux-%s/Fasta.pickle" %
          (Crux.Config.prefix, Crux.Config.version),
          verbose=(False if (not __debug__ or Crux.opts.quiet) else True),
          skinny=(False if __debug__ else True),
          logFile="%s/share/Crux-%s/Fasta.log" %
          (Crux.Config.prefix, Crux.Config.version))

    cdef _initReDna(self):
        return re.compile(r"""
    (^>\S+.*$)                         # description
  | ([ABCDGHKMRSTVWY\-NX]
     [ABCDGHKMRSTVWY\-NX \t\n\r\f\v]*) # characters
  | ([ \t\n\r\f\v]+)                   # whitespace

""", re.I | re.M | re.X)

    cdef _initReProtein(self):
        return re.compile(r"""
    (^>\S+.*$)                                  # description
  | ([ABCDEFGHIKLMNPQRSTUVWXYZ*\-]
     [ABCDEFGHIKLMNPQRSTUVWXYZ*\- \t\n\r\f\v]*) # characters
  | ([ \t\n\r\f\v]+)                            # whitespace
""", re.I | re.M | re.X)

    cdef void _appendToken(self, Token token) except *:
        if self.first is None:
            self.first = token
            self.last = token
        else:
            self.last.next = token
            token.prev = self.last
            self.last = token

    cpdef parse(self, str input, type charType=Dna, int begPos=0, int line=1,
      bint verbose=False):
        cdef int pos = begPos
        cdef object regex, m
        cdef int idx, start, end
        cdef Token token
        cdef int nesting

        self.verbose = verbose

        if charType is Dna:
            if _reDna is None:
                _reDna = self._initReDna()
            regex = _reDna
        else:
            assert charType is Protein
            if _reProtein is None:
                _reProtein = self._initReProtein()
            regex = _reProtein

        # Iteratively tokenize the input and feed the tokens to the parser.
        while pos < len(input):
            m = regex.match(input, pos)
            if m is None:
                raise SyntaxError(line, "Invalid token")
            idx = m.lastindex
            start = m.start(idx)
            end = m.end(idx)
            if idx == 1:
                self.token(TokenDescr(self, input, start, end, line))
            elif idx == 2:
                self.token(TokenChars(self, input, start, end, line))
                line += input.count('\n', start, end)
            elif idx == 3:
                TokenWhitespace(self, input, start, end, line)
                line += input.count('\n', start, end)
            else:
                assert False

            pos = end

        self.eoi()
