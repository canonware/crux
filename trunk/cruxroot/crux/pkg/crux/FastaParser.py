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

import re
import string

import crux.Exception

class Exception(crux.Exception):
    pass

class SyntaxError(Exception, SyntaxError):
    def __init__(self, message, line=None, column=0, char='', token=[]):
        self._message = message
        self._line = line
        self._column = column
        self._char = char
        self._token = token

    def __str__(self):
        if self._line != None:
            retval = "At %d:%d (token %r, char %r): %s" \
                     % (self._line, self._column,
                        string.join(self._token, ""), self._char,
                        self._message)
        else:
            retval = self._message

        return retval

class FastaParser(object):
    def __init__(self):
        pass

    # Parse input, which has either 'DNA' or 'protein' character data.
    def parse(self, input, chartype='DNA'):
        self._src = input
        self._srcOffset = 0
        self._line = 1
        self._column = 0
        self._state = None # None, 'label', 'chars'
        self._token = []

        # Initialize regular expressions.
        if chartype == 'DNA':
            # Match DNA character data.
            reChars = re.compile(r'[NXVHMDRWABSYCKGT-]', re.I)
        else:
            # Match protein character data.
            reChars = re.compile(r'[ABCDEFGHIKLMNPQRSTUVWXYZ-]', re.I)
        # Match non-whitespace.
        reNotWhitespace = re.compile(r'[^\r\n\t ]')

        # Read one character at a time and handle it using a DFA with three
        # states (None (starting state), 'label', and 'chars').
        (c, line, column) = self._getc()
        while (c != ""):
            if self._state == 'label':
                if c == "\n":
                    if len(self._token) == 0:
                        raise crux.FastaParser\
                              .SyntaxError("Empty label",
                                           line, column,
                                           c, self._token)
                        
                    self._tokenAccept(self.labelAccept, 'chars')
                else:
                    self._token.append(c)
            elif self._state == 'chars':
                if c == ">" and column == 0:
                    if len(self._token) == 0:
                        raise crux.FastaParser\
                              .SyntaxError("Missing character data",
                                           line, column,
                                           c, self._token)
                    self._tokenAccept(self.charsAccept, 'label')
                else:
                    if reChars.match(c):
                        # Append to token.
                        self._token.append(c)
                    elif reNotWhitespace.match(c):
                        raise crux.FastaParser\
                              .SyntaxError("Invalid character data",
                                           line, column,
                                           c, self._token)
            else: # Starting state.
                if c == ">":
                    self._state = 'label'

            # Get next character.
            (c, line, column) = self._getc()

        # Make sure that the input ends with character data.
        if self._state != 'chars':
            raise crux.FastaParser\
                  .SyntaxError("Input ended while reading label",
                               line, column,
                               c, self._token)

        # Accept the last token.
        if len(self._token) == 0:
            raise crux.FastaParser\
                  .SyntaxError("Missing character data",
                               line, column,
                               c, self._token)
        self._tokenAccept(self.charsAccept, None)

    # Return the current token (string).
    def token(self):
        return string.join(self._token, "")

    # Return the current line within the input.
    def line(self):
        return self._line

    def labelAccept(self):
        # Virtual method.
        pass

    def charsAccept(self):
        # Virtual method.
        pass

    # Get a character.
    def _getc(self):
        line = self._line
        column = self._column

        # Get a character.
        if type(self._src) == str:
            # String input.
            if self._srcOffset < len(self._src):
                c = self._src[self._srcOffset]
                self._srcOffset += 1
            else:
                c = ""
        else:
            # File input.
            c = self._src.read(1)

        # Update line and column info.
        if c == '\n':
            self._line += 1
            self._column = 0
        else:
            self._column += 1

        # Return the character, plus the line and column that c was read from.
        return (c, line, column)

    def _tokenAccept(self, method, newState):
        method()

        self._token = []
        self._state = newState
#EOF
