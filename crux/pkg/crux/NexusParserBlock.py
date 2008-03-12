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
import NexusParserStatement

# The default statement handler for the default block handler does not need to
# be as verbose as the generic default statement handler.
class _NexusParserStatement(NexusParserStatement.NexusParserStatement):
    # XXX Implement.
    pass

# A block is composed of the following children:
#
#  begin <blockName> <statement> ... <statement> end ;
class NexusParserBlock(NexusParserBase.NexusParserBase):
    def __init__(self, nexusParser, parent=None, siblingIndex=None,
                 offset=None, length=None, line=None, column=None):
        NexusParserBase.NexusParserBase.__init__(self, nexusParser, parent,
                                                 siblingIndex, offset, length,
                                                 line, column)

    def blockName(self):
        if self.finished():
            rVal = self.child(1).cText
        else:
            rVal = None

        return rVal

    def statementSearch(self, command):
        if self.finished():
            rVal = None
            nChildren = self.nChildren()
            for i in forints(nChildren - 2, 2):
                child = self.child()
                if child.commandName() == command:
                    rVal = child
                    break
        else:
            rVal = None

        return rVal

    def nStatements(self):
        if self.finished():
            rVal = self.nChildren() - 4
        else:
            rVal = 0

        return rVal

    def statementGet(self, index):
        if self.finished():
            if index < self.nChildren() - 4:
                rVal = self.child(index + 2)
            else:
                rVal = None
        else:
            rVal = None

        return rVal

    NexusParserStatement = _NexusParserStatement
