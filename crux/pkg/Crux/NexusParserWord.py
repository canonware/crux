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

class NexusParserWord(NexusParserBase.NexusParserBase):
    def __init__(self, nexusParser, parent=None, siblingIndex=None,
                 offset=None, length=None, line=None, column=None):
        NexusParserBase.NexusParserBase.__init__(self, nexusParser, parent,
                                                 siblingIndex, offset, length,
                                                 line, column)
        self._cText = None

    def cText(self):
        if self._cText == None:
            # XXX Canonize self.text, and store it as self._cText.
            pass

        return self._cText

    def nSeparators(self):
        if self.finished():
            rVal = self.nChildren()
        else:
            rVal = 0

        return rVal

    def separatorGet(self, index):
        if self.finished():
            if index < self.nChildren():
                rVal = self.child(index)
            else:
                rVal = None
        else:
            rVal = None

        return rVal
