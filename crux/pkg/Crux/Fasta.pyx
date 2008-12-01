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

cimport Parsing
from Crux.Character cimport Dna, Protein
from Crux.Taxa cimport Taxon
cimport Crux.Taxa as Taxa

global __name__

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
    def __init__(self, Parser parser, str raw):
        Parsing.Token.__init__(self, parser)

        self.raw = raw

    def __repr__(self):
        return Parsing.Token.__repr__(self) + (" (%r)" % self.raw)

cdef class TokenDescr(Token):
    "%token descr"
    def __init__(self, Parser parser, str raw):
        cdef object m

        Token.__init__(self, parser, raw)

        # Strip the leading '>' and discard comments.
        self.label = re.match(r">(\S+)", raw).group(1)
        # Convert '_' to ' '.
        self.label = self.label.replace('_', ' ')

        m = re.match(r">\S+\s+(.*)$", raw)
        if m is not None:
            self.comment = m.group(1)
        else:
            self.comment = None

cdef class TokenChars(Token):
    "%token chars"
    def __init__(self, Parser parser, str raw):
        Token.__init__(self, parser, raw)

        # Discard whitespace.
        self.chars = re.sub(r'[ \t\r\n]+', r'', raw)

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
    cpdef reduce(self, TokenDescr descr, Chars Chars):
        "%reduce descr Chars"
        (<Parser>self.parser)._addTaxon(Taxa.get(descr.label),
          "".join(Chars.chars))

cdef class Chars(Nonterm):
    "%nonterm"
    cpdef reduceOne(self, TokenChars chars):
        "%reduce chars"
        self.chars = [chars.chars]

    cpdef reduceExtend(self, Chars Chars, TokenChars chars):
        "%reduce Chars chars"
        self.chars = Chars.chars
        self.chars.append(chars.chars)

#
# End Nonterm.
#===============================================================================

# Regexs used to recognize all tokens.
cdef _reDna
cdef _reProtein

cdef Parsing.Spec _spec

cdef class Parser(Parsing.Lr):
    def __init__(self, CTMatrix matrix, Taxa.Map taxaMap=None,
      Parsing.Spec spec=None):
        global _spec

        if spec is None:
            if _spec is None:
                _spec = self._initSpec()
            spec = _spec
        Parsing.Lr.__init__(self, spec)

        self.matrix = matrix
        self.taxaMap = taxaMap

    cdef Parsing.Spec _initSpec(self):
        return Parsing.Spec([sys.modules[__name__]],
          pickleFile="%s/Crux/parsers/Fasta.pickle" % Crux.Config.datadir,
          verbose=(False if (not __debug__ or Crux.Config.quiet) else True),
          skinny=(False if __debug__ else True),
          logFile="%s/Crux/parsers/Fasta.log" % Crux.Config.datadir)

    cdef _initReDna(self):
        return re.compile(r"""
    (^>\S+.*$)                          # description
  | (^[ \t\n\r\f\v]*
     [ABCDGHKMRSTVWY\-NX]
     [ABCDGHKMRSTVWY\-NX \t\n\r\f\v]*$) # characters
  | (^[ \t\n\r\f\v]*$)                  # whitespace

""", re.I | re.X)

    cdef _initReProtein(self):
        return re.compile(r"""
    (^>\S+.*$)                                   # description
  | (^[ \t\n\r\f\v]*
     [ABCDEFGHIKLMNPQRSTUVWXYZ*\-]
     [ABCDEFGHIKLMNPQRSTUVWXYZ*\- \t\n\r\f\v]*$) # characters
  | (^[ \t\n\r\f\v]*$)                           # whitespace
""", re.I | re.X)

    cdef void _addTaxon(self, Taxon taxon, str chars) except *:
        if self.taxaMap.indGet(taxon) == -1:
            # Define a taxon mapping for this label.
            self.taxaMap.map(taxon, self.taxaMap.ntaxa)

        # Set the character data for this taxon.
        self.matrix.dataSet(taxon, chars)

    cpdef parse(self, lines, type charType=Dna, int line=1, bint verbose=False):
        cdef object regex, m
        cdef int idx, start, end
        cdef str l

        assert getattr3(lines, '__iter__', None) is not None

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
        for l in lines:
            m = regex.match(l)
            if m is None:
                raise SyntaxError(line, "Invalid token")
            idx = m.lastindex
            start = m.start(idx)
            end = m.end(idx)
            if idx == 1:
                self.token(TokenDescr(self, l))
            elif idx == 2:
                self.token(TokenChars(self, l))
            elif idx == 3:
                # Whitespace.
                pass
            else:
                assert False
            line += 1

        self.eoi()
