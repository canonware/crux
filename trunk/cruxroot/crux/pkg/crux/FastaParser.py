################################################################################
#
# <Copyright = jasone>
# <License>
#
################################################################################
#
# Version: Crux <Version = crux>
#
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

import _FastaParser

class Exception(_FastaParser.Exception):
    pass

class ValueError(Exception, _FastaParser.ValueError):
    pass

class TypeError(Exception, _FastaParser.TypeError):
    pass

class SyntaxError(Exception, _FastaParser.SyntaxError):
    pass

class FastaParser(_FastaParser.FastaParser):
    def __init__(self):
        pass

    # Parse input, which has either 'DNA' or 'protein' character data.
    def parse(self, input, charType='DNA'):
        self._charType = charType
        _FastaParser.FastaParser.parse(self, input, charType)

    # Return the character type being parsed.
    def charType(self):
        return self._charType

    def labelAccept(self):
        # Virtual method.
        pass

    def commentAccept(self):
        # Virtual method.
        pass

    def charsAccept(self):
        # Virtual method.
        pass
#EOF
