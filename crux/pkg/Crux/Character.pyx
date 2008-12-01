################################################################################
#
# The Character class encapsulates the functionality that is necessary to
# map character codes to bit vectors (as well as the reverse mapping).  Under
# normal circumstances, there are only a few Character instances, though
# many characters may refer to them.
#
# The Dna and Protein classes are pre-configured to
# support the standard codes.
#
################################################################################

import Crux.Exception

class Exception(Crux.Exception.Exception):
    pass

import exceptions

class ValueError(Exception, exceptions.ValueError):
    def __init__(self, str):
        self._str = str

    def __str__(self):
        return self._str

cimport Crux.DistMatrix as DistMatrix

cdef class Character:
    def __init__(self):
        self._pStates = {} # Primary states.
        self._aStates = {} # All states (including primary states).
        self._vals = {} # Reverse lookup of states.

    cpdef int nstates(self):
        return len(self._pStates)

    cpdef list codes(self):
        return self._aStates.keys()

    # Define state code.
    cpdef stateCodeAdd(self, str code):
        cdef int pval, val

        if self._aStates.has_key(code):
            raise ValueError("State already defined: %r" % code)

        # Define a key in _pStates and _aStates that makes it possible to get
        # the state's value in constant time.
        pval = len(self._pStates) + 1
        self._pStates[code] = pval
        val = 1 << (pval - 1)
        self._aStates[code] = val

        # Create a reverse lookup (val --> code).
        self._vals[val] = code

    cpdef ambiguityCodeAdd(self, str code, list oStates=None):
        cdef int val
        cdef str oState

        if self._aStates.has_key(code):
            raise ValueError("State already defined: %r" % code)

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

    cpdef aliasCodeAdd(self, str code, str oState):
        cdef int val

        if self._aStates.has_key(code):
            raise ValueError("State already defined: %r" % code)
        if not self._aStates.has_key(oState):
            raise ValueError("State not defined: %r" % oState)

        # Define a key in _aStates that makes it possible to get the state's
        # value in constant time.
        val = self._aStates[oState]
        self._aStates[code] = val

        # Create a reverse lookup (val --> code) if one doesn't already exist.
        if not self._vals.has_key(val):
            self._vals[val] = code

    # Given a state code, return the associated value.
    cpdef int code2val(self, str code):
        if not self._aStates.has_key(code):
            raise ValueError("State not defined: %r" % code)

        return self._aStates[code]

    # Given an index value, return the primary associated state key.
    cpdef str val2code(self, int val):
        if not self._vals.has_key(val):
            raise ValueError("Value not defined: %r" % val)

        return self._vals[val]

cdef class Dna(Character):
    def __init__(self):
        Character.__init__(self)

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

cdef class Protein(Character):
    def __init__(self):
        Character.__init__(self)

        states = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
                  'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
        for state in states:
            self.stateCodeAdd(state)

        self.ambiguityCodeAdd('-', states)
