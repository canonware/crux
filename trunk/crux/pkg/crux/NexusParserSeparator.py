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

class NexusParserSeparator(NexusParserBase.NexusParserBase):
    def __init__(self, nexusParser, parent=None, siblingIndex=None,
                 offset=None, length=None, line=None, column=None):
        NexusParserBase.NexusParserBase.__init__(self, nexusParser, parent,
                                                 siblingIndex, offset, length,
                                                 line, column)
        self._cText = None

    def cText(self):
        if self._cText == None:
            # XXX Canonize self.text, and store it as self._cText.  Only special
            # comments such as [!...] and [\b] have a non-empty canonized form.
            pass

        return self._cText
