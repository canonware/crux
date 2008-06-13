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

import NexusParser
import NexusParserBase

import Crux

class NexusParserSyntaxError(NexusParserBase.NexusParserBase):
    def __init__(self, nexusParser, parent=None, siblingIndex=None,
                 offset=None, length=None, line=None, column=None,
                 message=None):
        NexusParserBase.NexusParserBase.__init__(self, nexusParser, parent,
                                                 siblingIndex, offset, length,
                                                 line, column)
        self._message = message

    def accessCallback(self):
        raise Crux.NexusParser.SyntaxError(
            self.nexusParser().inputName(),
            self.line(), self.column(), self.offset(),
            self._message)
