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
# Recursive descent parser that implements Newick Tree Format parsing.  This
# implementation expands from Gary Olsen's interpretation of the Newick Tree
# Format Standard, available at:
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
# implementation allows "zero length strings of printing characters".  This has
# the implication that the label productions always accept.
#
################################################################################
#
# Notation:
#
#   Production: production_name ::= value
#
#   Production reference: <production>
#
#   Optional (0-1 occurrences): [optional]
#
#   Alternation: <this> | <that>
#
#   Regular expression: /regex/
#
#   Special character: "[", "]", etc.
#
#   Whitespace between components of a production alternative indicates optional
#   whitespace (<ws>), which consists of whitespace and/or comments.
#
# Notes:
#
#   "_" in <unquoted_label> is converted to " ".
#
# tree ::= [ <descendant_list> ][ <root_label> ][ : <branch_length> ] ;
#
# descendant_list ::= ( <subtree> [ , <subtree> ] )
#
# subtree ::= <descendant_list> [ <internal_label> ][ : <branch_length> ]
#           | <leaf_label> [ : <branch_length> ]
#
# root_label ::= <label>
#
# internal_label ::= <label>
#
# leaf_label ::= <label>
#
# label ::= <unquoted_label>
#         | <quoted_label>
#         |
#
# unquoted_label ::= <ulabel_char>[<unquoted_label>]
#
# ulabel_char ::= /[^ ()[\]':;,]/
#
# quoted_label ::= '<qlabel_char>[<quoted_label>]'
#
# qlabel_char ::= /[^']|('')/
#
# branch_length ::= <int>[.<int>]
#                 | +<int>[.<int>]
#                 | -<int>[.<int>]
#
# int ::= <digit>[<int>]
#
# digit ::= /[0-9]/
#
# ws ::= [<comment>][<whitespace>][<ws>]
#      | <e>
#
# whitespace ::= / \t\r\n/
#
# comment ::= \[ [<comment_chars>] [<comment>] [<comment_chars>] \]
#
# comment_chars ::= /[^[]/
#                 | 
#
# e ::= <epsilon (empty production)>
#
################################################################################

import sys
import re

import crux.Exception

class Exception(crux.Exception):
    pass

class SyntaxError(Exception, SyntaxError):
    def __init__(self, message, offset=None, char=None, token=None):
        self._message = message
        self._offset = offset
        self._char = char
        self._token = token

    def __str__(self):
        if self._offset != None:
            retval = "At offset %d (char '%s') (token \"%s\"): %s" \
                     % (self._offset, self._char, self._token, self._message)
        else:
            retval = self._message

        return retval

class NewickParser(object):
    def __init__(self):
        pass

    # Parse input and call the *Accept methods for each token that is accepted.
    def parse(self, input):
        self._src = input
        self._srcOffset = -1
        self._c = "?"
        self._lookaheadC = None
        self._token = ""

        self._treeProduction()

    # Return the current token (string).
    def token(self):
        return self._token

    # Return the current offset within the input.
    def offset(self):
        return self._srcOffset

    def openParenAccept(self):
        # Virtual method.
        pass

    def closeParenAccept(self):
        # Virtual method.
        pass

    def rootLabelAccept(self):
        # Virtual method.
        pass

    def internalLabelAccept(self):
        # Virtual method.
        pass

    def leafLabelAccept(self):
        # Virtual method.
        pass

    def colonAccept(self):
        # Virtual method.
        pass

    def lengthAccept(self):
        # Virtual method.
        pass

    def commaAccept(self):
        # Virtual method.
        pass

    def semicolonAccept(self):
        # Virtual method.
        pass

    def commentAccept(self):
        # Virtual method.
        pass

    def whitespaceAccept(self):
        # Virtual method.
        pass

    #
    # Instance data documentation.
    #
    # _src : Source, either a file or a string.
    #
    # _srcOffset : Byte offset from the beginning of the source, starting at 0.
    #
    # _c : Current character, stored as a one byte string.  Rather than
    #      re-defining this for every character, the current character is put
    #      into one string which is used over and over.
    #
    # _lookaheadC : Lookahead character, stored as a one byte string.  This
    #                supports :p_ungetc.
    #
    # _token : Current token (or partial token, if not yet accepted).

    #
    # Private methods.
    #

    # Get a character.
    def _getc(self):
        if self._lookaheadC != None:
	    # Lookahead character defined.
            retval = self._lookaheadC
            self._lookaheadC = None
        else:
            # Increment the offset.
            self._srcOffset += 1

            if type(self._src) == str:
                # String input.
                if self._srcOffset >= len(self._src):
                    raise crux.NewickParser.SyntaxError("End of input reached",
                                                        self._srcOffset,
                                                        self._c,
                                                        self._token)

                self._c = self._src[self._srcOffset]
            else:
                # File input.
                self._c = read(self._src, 1)
                if self._c == "":
                    raise crux.NewickParser.SyntaxError("End of input reached",
                                                        self._srcOffset,
                                                        self._c,
                                                        self._token)

            retval = self._c

        # Append character to token.
        self._token += self._c

        return retval

    # Set lookahead_c so that it will be used by the next _getc call.
    def _ungetc(self):
        self._lookaheadC = self._c

        # Remove character from token.
        self._token = self._token[:-1]

    # Call the *Accept specified method, then reset the internal state in
    # preparation for the next token.
    def _tokenAccept(self, method):
        method()

        self._token = ""

    # Top level production.
    def _treeProduction(self):
        self._wsProduction()
        self._descendantListProduction()
        self._wsProduction()
        self._labelProduction(self.rootLabelAccept)
        self._wsProduction()

        if self._getc() == ":":
            self._tokenAccept(self.colonAccept)
            self._wsProduction()
            self._branchLengthProduction()
            self._wsProduction()
        else:
            self._ungetc()

        if self._getc() != ";":
            raise crux.NewickParser.SyntaxError("';' expected",
                                                self._srcOffset,
                                                self._c,
                                                self._token)

        self._tokenAccept(self.semicolonAccept)

    # Return True if a descendant_list is accepted, False otherwise.
    #
    # accepted _descendantListProduction()
    def _descendantListProduction(self):
        if self._getc() == "(":
            self._tokenAccept(self.openParenAccept)

            self._wsProduction()
            self._subtree()
            self._wsProduction()

            while True:
                if self._getc() == ",":
                    self._tokenAccept(self.commaAccept)
                    self._wsProduction()
                    self._subtree()
                    self._wsProduction()
                else:
                    self._ungetc()
                    break

            if self._getc() != ")":
                raise crux.NewickParser.SyntaxError("',' or ')' expected",
                                                    self._srcOffset,
                                                    self._c,
                                                    self._token)

            self._tokenAccept(self.closeParenAccept)

            retval = True
        else:
            self._ungetc()
            retval = False

        return retval

    def _subtree(self):
        self._wsProduction()
        accepted = self._descendantListProduction()
        self._wsProduction()
        if accepted:
            self._labelProduction(self.internalLabelAccept)
        else:
            self._labelProduction(self.leafLabelAccept)
        self._wsProduction()

        if self._getc() == ":":
            self._tokenAccept(self.colonAccept)
            self._wsProduction()
            self._branchLengthProduction()
        else:
            self._ungetc()

    def _labelProduction(self, accept_method):
        if not self._quotedLabelProduction(accept_method):
            if not self._unquotedLabelProduction(accept_method):
                # Accept a zero length label.
                self._tokenAccept(accept_method)

    # If an unquoted label is accepted, return True, otherwise return False.
    #
    # accepted _unquotedLabelProduction(accept_method)
    def _unquotedLabelProduction(self, accept_method):
        if self._ulabelCharProduction():
            # Consume all ulabel_char characters.
            while True:
                if not self._ulabelCharProduction():
                    break

            # Convert '_' to ' ' before accepting the token.
            uscore = re.compile(r'_')
            self._token = uscore.sub(' ', self._token)

            self._tokenAccept(accept_method)
            retval = True
        else:
            retval = False

        return retval

    # If an unquoted label character is accepted, return True, otherwise return
    # False.
    #
    # accepted _ulabelCharProduction()
    def _ulabelCharProduction(self):
        ulabel_char = re.compile(r'[^ ()[\]\':;,]')
        if ulabel_char.match(self._getc()):
            retval = True
        else:
            self._ungetc()
            retval = False

        return retval

    # If a quoted label is accepted, return True, otherwise return False.
    #
    # accepted _quotedLabelProduction(accept_method)
    def _quotedLabelProduction(self, accept_method):
        if self._getc() == "'":
            while True:
                if not self._qlabelCharProduction():
                    break

            # Remove quotes before accepting the token.
            self._token = self._token[1:-1]

            # Collapse '' to ' before accepting the token.
            qcollapse = re.compile("''")
            self._token = qcollapse.sub("'", self._token)

            self._tokenAccept(accept_method)

            retval = True
        else:
            self._ungetc()
            retval = False

        return retval

    # Since there is only one lookahead character, this method is responsible
    # for reading the ' that terminates the quoted label.
    #
    # accepted _qlabelCharProduction()
    def _qlabelCharProduction(self):
        if self._getc() == "'":
            if self._getc() == "'":
                retval = True
            else:
                self._ungetc()
                retval = False
        else:
            retval = True

        return retval

    def _branchLengthProduction(self):
        sign = re.compile(r'[+-]')
        if not sign.match(self._getc()):
            self._ungetc()

        self._intProduction()

        if self._getc() == ".":
            self._intProduction()
        else:
            self._ungetc()

        self._tokenAccept(self.lengthAccept)

    def _intProduction(self):
        if not self._digitProduction():
            raise crux.NewickParser.SyntaxError("Expected digit",
                                                self._srcOffset,
                                                self._c,
                                                self._token)

        while True:
            if not self._digitProduction():
                break

    # accepted _digitProduction()
    def _digitProduction(self):
        digit = re.compile(r'\d')
        if digit.match(self._getc()):
            retval = True
        else:
            self._ungetc()
            retval = False

        return retval

    # accepted _commentProduction()
    def _commentProduction(self):
        if self._getc() == "[":
            while True:
                if self._getc() == "]":
                    break
                else:
                    self._ungetc()
                    if not self._commentProduction():
                        self._getc()
            retval = True
        else:
            self._ungetc()
            retval = False

        return retval

    # accepted _whitespaceProduction()
    def _whitespaceProduction(self):
        whitespace = re.compile(r'\s')
        if whitespace.match(self._getc()):
            while True:
                if not whitespace.match(self._getc()):
                    self._ungetc()
                    break
            retval = True
        else:
            self._ungetc()
            retval = False

        return retval

    def _wsProduction(self):
        again = True
        while again:
            again = False

            # Match a comment.
            if self._commentProduction():
                self._tokenAccept(self.commentAccept)
                again = True

            # Match whitespace.
            if self._whitespaceProduction():
                self._tokenAccept(self.whitespaceAccept)
                again = True
#EOF
