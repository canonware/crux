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
# root_label ::= label
#
# internal_label ::= label
#
# leaf_label ::= label
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

class _newick_error: pass

class newick(object):
    def __init__(self):
        pass

    # Parse input and call the *_accept methods for each token that is accepted.
    def parse(self, input):
        self._src = input
        self._src_offset = -1
        self._c = "?"
        self._lookahead_c = None
        self._token = ""

        try:
            self._p_tree()
            retval = False
        except _newick_error:
            retval = True

        return retval

    # Return the current token (string).
    def token(self):
        return self._token

    # Return the current offset within the input.
    def offset(self):
        return self._src_offset

    # Handle an error, and unwind to the parse method.
    def error_raise(self, str):
        self.error_print(str)

        X = _newick_error()
        raise X

    # Print an error message.
    def error_print(self, str):
        print >> sys.stderr, "At offset %d (char '%s') (token \"%s\"): %s" \
              % (self._src_offset, self._c, self._token, str)

    def open_paren_accept(self):
        # Virtual method.
        pass

    def close_paren_accept(self):
        # Virtual method.
        pass

    def root_label_accept(self):
        # Virtual method.
        pass

    def internal_label_accept(self):
        # Virtual method.
        pass

    def leaf_label_accept(self):
        # Virtual method.
        pass

    def colon_accept(self):
        # Virtual method.
        pass

    def length_accept(self):
        # Virtual method.
        pass

    def comma_accept(self):
        # Virtual method.
        pass

    def semicolon_accept(self):
        # Virtual method.
        pass

    def comment_accept(self):
        # Virtual method.
        pass

    def whitespace_accept(self):
        # Virtual method.
        pass

    #
    # Instance data documentation.
    #
    # _src : Source, either a file or a string.
    #
    # _src_offset : Byte offset from the beginning of the source, starting at 0.
    #
    # _c : Current character, stored as a one byte string.  Rather than
    #      re-defining this for every character, the current character is put
    #      into one string which is used over and over.
    #
    # _lookahead_c : Lookahead character, stored as a one byte string.  This
    #                supports :p_ungetc.
    #
    # _token : Current token (or partial token, if not yet accepted).

    #
    # Private methods.
    #

    # Get a character.
    def _p_getc(self):
        if self._lookahead_c != None:
	    # Lookahead character defined.
            retval = self._lookahead_c
            self._lookahead_c = None
        else:
            # Increment the offset.
            self._src_offset += 1

            if type(self._src) == str:
                # String input.
                if self._src_offset >= len(self._src):
                    self.error_raise("End of input reached")

                self._c = self._src[self._src_offset]
            else:
                # File input.
                self._c = read(self._src, 1)
                if self._c == "":
                    self.error_raise("End of input reached")

            retval = self._c

        # Append character to token.
        self._token += self._c

        return retval

    # Set lookahead_c so that it will be used by the next _p_getc call.
    def _p_ungetc(self):
        self._lookahead_c = self._c

        # Remove character from token.
        self._token = self._token[:-1]

    # Call the *_accept specified method, then reset the internal state in
    # preparation for the next token.
    def _p_token_accept(self, method):
        method()

        self._token = ""

    # Top level production.
    def _p_tree(self):
        self._p_ws()
        self._p_descendant_list()
        self._p_ws()
        self._p_label(self.root_label_accept)
        self._p_ws()

        if self._p_getc() == ":":
            self._p_token_accept(self.colon_accept)
            self._p_ws()
            self._p_branch_length()
            self._p_ws()
        else:
            self._p_ungetc()

        if self._p_getc() != ";":
            self.error_raise("';' expected")

        self._p_token_accept(self.semicolon_accept)

    # Return True if a descendant_list is accepted, False otherwise.
    #
    # accepted _p_descendant_list()
    def _p_descendant_list(self):
        if self._p_getc() == "(":
            self._p_token_accept(self.open_paren_accept)

            self._p_ws()
            self._p_subtree()
            self._p_ws()

            while True:
                if self._p_getc() == ",":
                    self._p_token_accept(self.comma_accept)
                    self._p_ws()
                    self._p_subtree()
                    self._p_ws()
                else:
                    self._p_ungetc()
                    break

            if self._p_getc() != ")":
                self.error_raise("',' or ')' expected")

            self._p_token_accept(self.close_paren_accept)

            retval = True
        else:
            self._p_ungetc()
            retval = False

        return retval

    def _p_subtree(self):
        self._p_ws()
        accepted = self._p_descendant_list()
        self._p_ws()
        if accepted:
            self._p_label(self.internal_label_accept)
        else:
            self._p_label(self.leaf_label_accept)
        self._p_ws()

        if self._p_getc() == ":":
            self._p_token_accept(self.colon_accept)
            self._p_ws()
            self._p_branch_length()
        else:
            self._p_ungetc()

    def _p_label(self, accept_method):
        if not self._p_quoted_label(accept_method):
            if not self._p_unquoted_label(accept_method):
                # Accept a zero length label.
                self._p_token_accept(accept_method)

    # If an unquoted label is accepted, return True, otherwise return False.
    #
    # accepted _p_unquoted_label(accept_method)
    def _p_unquoted_label(self, accept_method):
        if self._p_ulabel_char():
            # Consume all ulabel_char characters.
            while True:
                if not self._p_ulabel_char():
                    break

            # Convert '_' to ' ' before accepting the token.
            uscore = re.compile(r'_')
            self._token = uscore.sub(' ', self._token)

            self._p_token_accept(accept_method)
            retval = True
        else:
            retval = False

        return retval

    # If an unquoted label character is accepted, return True, otherwise return
    # False.
    #
    # accepted _p_ulabel_char()
    def _p_ulabel_char(self):
        ulabel_char = re.compile(r'[^ ()[\]\':;,]')
        if ulabel_char.match(self._p_getc()):
            retval = True
        else:
            self._p_ungetc()
            retval = False

        return retval

    # If a quoted label is accepted, return True, otherwise return False.
    #
    # accepted _p_quoted_label(accept_method)
    def _p_quoted_label(self, accept_method):
        if self._p_getc() == "'":
            while True:
                if not self._p_qlabel_char():
                    break

            # Remove quotes before accepting the token.
            self._token = self._token[1:-1]

            # Collapse '' to ' before accepting the token.
            qcollapse = re.compile("''")
            self._token = qcollapse.sub("'", self._token)

            self._p_token_accept(accept_method)

            retval = True
        else:
            self._p_ungetc()
            retval = False

        return retval

    # Since there is only one lookahead character, this method is responsible
    # for reading the ' that terminates the quoted label.
    #
    # accepted _p_qlabel_char()
    def _p_qlabel_char(self):
        if self._p_getc() == "'":
            if self._p_getc() == "'":
                retval = True
            else:
                self._p_ungetc()
                retval = False
        else:
            retval = True

        return retval

    def _p_branch_length(self):
        sign = re.compile(r'[+-]')
        if not sign.match(self._p_getc()):
            self._p_ungetc()

        self._p_int()

        if self._p_getc() == ".":
            self._p_int()
        else:
            self._p_ungetc()

        self._p_token_accept(self.length_accept)

    def _p_int(self):
        if not self._p_digit():
            self.error_raise("Expected digit")

        while True:
            if not self._p_digit():
                break

    # accepted _p_digit()
    def _p_digit(self):
        digit = re.compile(r'\d')
        if digit.match(self._p_getc()):
            retval = True
        else:
            self._p_ungetc()
            retval = False

        return retval

    # accepted _p_comment()
    def _p_comment(self):
        if self._p_getc() == "[":
            while True:
                if self._p_getc() == "]":
                    break
                else:
                    self._p_ungetc()
                    if not self._p_comment():
                        self._p_getc()
            retval = True
        else:
            self._p_ungetc()
            retval = False

        return retval

    # accepted _p_whitespace()
    def _p_whitespace(self):
        whitespace = re.compile(r'\s')
        if whitespace.match(self._p_getc()):
            while True:
                if not whitespace.match(self._p_getc()):
                    self._p_ungetc()
                    break
            retval = True
        else:
            self._p_ungetc()
            retval = False

        return retval

    def _p_ws(self):
        again = True
        while again:
            again = False

            # Match a comment.
            if self._p_comment():
                self._p_token_accept(self.comment_accept)
                again = True

            # Match whitespace.
            if self._p_whitespace():
                self._p_token_accept(self.whitespace_accept)
                again = True
#EOF
