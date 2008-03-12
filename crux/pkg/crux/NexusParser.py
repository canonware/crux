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

from C_NexusParser import *

import NexusParserBase
import NexusParserSyntaxError
import NexusParserFile
import NexusParserBlock
import NexusParserStatement
import NexusParserWord
import NexusParserSeparator

import crux

class Exception(crux.Exception):
    pass

class SyntaxError(Exception, SyntaxError):
    def __init__(self, inputName, line, column, offset, message):
        self._inputName = inputName
        self._line = line
        self._column = column
        self._offset = offset
        self._message = message

    def __str__(self):
        if self._line != None:
            line = "%d" % self._line
        else:
            line = "?"

        if self._column != None:
            column = "%d" % self._column
        else:
            column = "?"

        if self._offset != None:
            offset = "%d" % self._offset
        else:
            offset = "?"

        rVal = "At %s:%d:%d(%d): %s" \
               % (self._inputName, line, column, offset, self._message)

        return rVal

class NexusParser(C_NexusParser):
    def __init__(self, input, inputName="<?>"):
        pass

    NexusParserFile = NexusParserFile.NexusParserFile
