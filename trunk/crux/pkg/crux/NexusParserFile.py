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
import NexusParserBlock

# XXX Implement standard block handler classes.
class _NexusParserBlockTaxa(NexusParserBlock.NexusParserBlock):
    def __init__(self, nexusParser, parent=None, siblingIndex=None,
                 offset=None, length=None, line=None, column=None):
        NexusParserBase.NexusParserBase.__init__(self, nexusParser, parent,
                                                 siblingIndex, offset, length,
                                                 line, column)

class _NexusParserBlockCharacters(NexusParserBlock.NexusParserBlock):
    def __init__(self, nexusParser, parent=None, siblingIndex=None,
                 offset=None, length=None, line=None, column=None):
        NexusParserBase.NexusParserBase.__init__(self, nexusParser, parent,
                                                 siblingIndex, offset, length,
                                                 line, column)

class _NexusParserBlockUnaligned(NexusParserBlock.NexusParserBlock):
    def __init__(self, nexusParser, parent=None, siblingIndex=None,
                 offset=None, length=None, line=None, column=None):
        NexusParserBase.NexusParserBase.__init__(self, nexusParser, parent,
                                                 siblingIndex, offset, length,
                                                 line, column)

class _NexusParserBlockDistances(NexusParserBlock.NexusParserBlock):
    def __init__(self, nexusParser, parent=None, siblingIndex=None,
                 offset=None, length=None, line=None, column=None):
        NexusParserBase.NexusParserBase.__init__(self, nexusParser, parent,
                                                 siblingIndex, offset, length,
                                                 line, column)

class _NexusParserBlockSets(NexusParserBlock.NexusParserBlock):
    def __init__(self, nexusParser, parent=None, siblingIndex=None,
                 offset=None, length=None, line=None, column=None):
        NexusParserBase.NexusParserBase.__init__(self, nexusParser, parent,
                                                 siblingIndex, offset, length,
                                                 line, column)

class _NexusParserBlockAssumptions(NexusParserBlock.NexusParserBlock):
    def __init__(self, nexusParser, parent=None, siblingIndex=None,
                 offset=None, length=None, line=None, column=None):
        NexusParserBase.NexusParserBase.__init__(self, nexusParser, parent,
                                                 siblingIndex, offset, length,
                                                 line, column)

class _NexusParserBlockCodons(NexusParserBlock.NexusParserBlock):
    def __init__(self, nexusParser, parent=None, siblingIndex=None,
                 offset=None, length=None, line=None, column=None):
        NexusParserBase.NexusParserBase.__init__(self, nexusParser, parent,
                                                 siblingIndex, offset, length,
                                                 line, column)

class _NexusParserBlockTrees(NexusParserBlock.NexusParserBlock):
    def __init__(self, nexusParser, parent=None, siblingIndex=None,
                 offset=None, length=None, line=None, column=None):
        NexusParserBase.NexusParserBase.__init__(self, nexusParser, parent,
                                                 siblingIndex, offset, length,
                                                 line, column)

class _NexusParserBlockNotes(NexusParserBlock.NexusParserBlock):
    def __init__(self, nexusParser, parent=None, siblingIndex=None,
                 offset=None, length=None, line=None, column=None):
        NexusParserBase.NexusParserBase.__init__(self, nexusParser, parent,
                                                 siblingIndex, offset, length,
                                                 line, column)

# A file is composed of the following children:
#
#  #nexus <block> ... <block>
class NexusParserFile(NexusParserBase.NexusParserBase):
    def __init__(self, nexusParser, parent=None, siblingIndex=None,
                 offset=None, length=None, line=None, column=None):
        NexusParserBase.NexusParserBase.__init__(self, nexusParser, parent,
                                                 siblingIndex, offset, length,
                                                 line, column)

    def nBlocks(self):
        if self.finished():
            rVal = self.nChildren() - 1
        else:
            rVal = 0

        return rVal

    def blockGet(self, index):
        if self.finished():
            if index < self.nChildren - 1:
                rVal = self.child(index + 1)
            else:
                rVal = None
        else:
            rVal = None

        return rVal

    def blockSearch(self, name):
        if self.finished():
            rVal = None
            nChildren = self.nChildren()
            for i in forints(nChildren, 1):
                child = self.child(i)
                if child.blockName() == name:
                    rVal = child
                    break
        else:
            rVal = None

        return rVal
        
    NexusParserBlock = NexusParserBlock.NexusParserBlock
    NexusParserBlockTaxa = _NexusParserBlockTaxa
    NexusParserBlockCharacters = _NexusParserBlockCharacters
    NexusParserBlockUnaligned = _NexusParserBlockUnaligned
    NexusParserBlockDistances = _NexusParserBlockDistances
    NexusParserBlockSets = _NexusParserBlockSets
    NexusParserBlockAssumptions = _NexusParserBlockAssumptions
    NexusParserBlockCodons = _NexusParserBlockCodons
    NexusParserBlockTrees = _NexusParserBlockTrees
    NexusParserBlockNotes = _NexusParserBlockNotes
