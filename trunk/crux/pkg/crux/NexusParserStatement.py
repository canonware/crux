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

import NexusParserBase

# A statement is composed of the following children:
#
#  <command> <word> ... <word> ;
class NexusParserStatement(NexusParserBase.NexusParserBase):
    def __init__(self, nexusParser, parent=None, siblingIndex=None,
                 offset=None, length=None, line=None, column=None):
        NexusParserBase.NexusParserBase.__init__(self, nexusParser, parent,
                                                 siblingIndex, offset, length,
                                                 line, column)

    def commandName(self):
        if self.finished():
            rVal = self.child(0).cText
        else:
            rVal = None

        return rVal

    def nWords(self):
        if self.finished():
            rVal = self.nChildren() - 2
        else:
            rVal = 0

        return rVal

    def wordGet(self, index):
        if self.finished():
            if index < self.nChildren() - 2:
                rVal = self.child(index + 1)
            else:
                rVal = None
        else:
            rVal = None

        return rVal
