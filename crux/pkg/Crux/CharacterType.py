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
# The CharacterType class encapsulates the functionality that is necessary to
# map character codes to bit vectors (as well as the reverse mapping).  Under
# normal circumstances, there are only a few CharacterType instances, though
# many characters may refer to them.
#
# The DnaCharacterType and ProteinCharacterType classes are pre-configured to
# support the standard codes.
#
################################################################################

import DistMatrix
import Crux

class Exception(Crux.Exception):
    pass

class ValueError(Exception, ValueError):
    def __init__(self, str):
        self._str = str

    def __str__(self):
        return self._str

class CharacterType(object):
    def __init__(self):
        self._pStates = {} # Primary states.
        self._aStates = {} # All states (including primary states).
        self._vals = {} # Reverse lookup of states.

    def nstates(self):
        return len(self._pStates)

    def codes(self):
        return self._aStates.keys()

    # Define state code.
    def stateCodeAdd(self, code):
        if self._aStates.has_key(code):
            raise Crux.CharacterType\
                  .ValueError("State already defined: %r" % code)

        # Define a key in _pStates and _aStates that makes it possible to get
        # the state's value in constant time.
        pval = len(self._pStates) + 1
        self._pStates[code] = pval
        val = 1 << (pval - 1)
        self._aStates[code] = val

        # Create a reverse lookup (val --> code).
        self._vals[val] = code

    def ambiguityCodeAdd(self, code, oStates=None):
        if self._aStates.has_key(code):
            raise Crux.CharacterType\
                  .ValueError("State already defined: %r" % code)

        # Define a key in _ambiguities that makes it possible to get the state's
        # value in constant time.
        val = 0
        if oStates != None:
            for oState in oStates:
                val |= self._aStates[oState]
        self._aStates[code] = val

        # Create a reverse lookup (val --> code) if one doesn't already exist.
        if not self._vals.has_key(val):
            self._vals[val] = code

    def aliasCodeAdd(self, code, oState):
        if self._aStates.has_key(code):
            raise Crux.CharacterType\
                  .ValueError("State already defined: %r" % code)
        if not self._aStates.has_key(oState):
            raise Crux.CharacterType.ValueError("State not defined: %r"
                                                % oState)

        # Define a key in _aStates that makes it possible to get the state's
        # value in constant time.
        val = self._aStates[oState]
        self._aStates[code] = val

        # Create a reverse lookup (val --> code) if one doesn't already exist.
        if not self._vals.has_key(val):
            self._vals[val] = code

    # Given a state code, return the associated value.
    def code2val(self, code):
        if not self._aStates.has_key(code):
            raise Crux.CharacterType.ValueError("State not defined: %r" % code)

        return self._aStates[code]

    # Given an index value, return the primary associated state key.
    def val2code(self, val):
        if not self._vals.has_key(val):
            raise Crux.CharacterType.ValueError("Value not defined: %r" % val)

        return self._vals[val]

class DnaCharacterType(CharacterType):
    def __init__(self):
        CharacterType.__init__(self)

        states = ['T', 'G', 'C', 'A']
        for state in states:
            self.stateCodeAdd(state)

        ambiguities = {'K' : ['T', 'G'          ],
                       'Y' : ['T',      'C'     ],
                       'S' : [     'G', 'C'     ],
                       'B' : ['T', 'G', 'C'     ],
                       'W' : ['T',           'A'],
                       'R' : [     'G',      'A'],
                       'D' : ['T', 'G',      'A'],
                       'M' : [          'C', 'A'],
                       'H' : ['T',      'C', 'A'],
                       'V' : [     'G', 'C', 'A'],
                       '-' : ['T', 'G', 'C', 'A']}
        for ambiguity in ambiguities:
            self.ambiguityCodeAdd(ambiguity, ambiguities[ambiguity])

        aliases = {'t': 'T',
                   'g': 'G',
                   'k': 'K',
                   'c': 'C',
                   'y': 'Y',
                   's': 'S',
                   'b': 'B',
                   'a': 'A',
                   'w': 'W',
                   'r': 'R',
                   'd': 'D',
                   'm': 'M',
                   'h': 'H',
                   'v': 'V',
                   'N': '-',
                   'n': '-',
                   'X': '-',
                   'x': '-'}
        for alias in aliases:
            self.aliasCodeAdd(alias, aliases[alias])

class ProteinCharacterType(CharacterType):        
    def __init__(self):
        CharacterType.__init__(self)

        states = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
                  'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
        for state in states:
            self.stateCodeAdd(state)

        self.ambiguityCodeAdd('-', states)
