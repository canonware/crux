#===============================================================================
#
# The Character class encapsulates the functionality that is necessary to
# map character codes to bit vectors (as well as the reverse mapping).  Under
# normal circumstances, there are only a few Character instances, though
# many characters may refer to them.
#
# The dna and protein instances are pre-configured to support the standard
# codes.
#
#===============================================================================

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
    def __init__(self, list states=[], dict ambiguities={}, dict aliases={}):
        cdef str state, ambiguity, alias

        self.any = 0 # Union of all primary states.
        self._pStates = {} # Primary states.
        self._aStates = {} # All states (including primary states).
        self._vals = {} # Reverse lookup of states.

        for state in states:
            self._stateCodeAdd(state)

        for ambiguity in ambiguities:
            self._ambiguityCodeAdd(ambiguity, ambiguities[ambiguity])

        for alias in aliases:
            self._aliasCodeAdd(alias, aliases[alias])

    cpdef int nstates(self):
        return len(self._pStates)

    cpdef list codes(self):
        return self._aStates.keys()

    # Define state code.
    cdef void _stateCodeAdd(self, str code) except *:
        cdef int pval, val

        if self._aStates.has_key(code):
            raise ValueError("State already defined: %r" % code)

        # Define a key in _pStates and _aStates that makes it possible to get
        # the state's value in constant time.
        pval = len(self._pStates) + 1
        self._pStates[code] = pval
        val = 1 << (pval - 1)
        self._aStates[code] = val

        # Merge into any.
        self.any |= val

        # Create a reverse lookup (val --> code).
        self._vals[val] = code

    cdef void _ambiguityCodeAdd(self, str code, list oStates) except *:
        cdef int val
        cdef str oState

        if self._aStates.has_key(code):
            raise ValueError("State already defined: %r" % code)

        # Define a key in _ambiguities that makes it possible to get the state's
        # value in constant time.
        val = 0
        if oStates is not None:
            for oState in oStates:
                val |= self._aStates[oState]
        self._aStates[code] = val

        # Create a reverse lookup (val --> code) if one doesn't already exist.
        if not self._vals.has_key(val):
            self._vals[val] = code

    cdef void _aliasCodeAdd(self, str code, str oState) except *:
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
    cpdef int code2val(self, str code) except -1:
        if not self._aStates.has_key(code):
            raise ValueError("State not defined: %r" % code)

        return self._aStates[code]

    # Given an index value, return the primary associated state key.
    cpdef str val2code(self, int val):
        if not self._vals.has_key(val):
            raise ValueError("Value not defined: %r" % val)

        return self._vals[val]

#===============================================================================

# Forward declaration.
cdef class Dna(Character)

cdef Dna dna

cdef class Dna(Character):
    """
        Character codes for DNA/RNA nucleotides.
    """
    def __init__(self):
        global dna

        assert dna is None
        Character.__init__(self,
          states=['A', 'C', 'G', 'T'],
          ambiguities={
            'K' : [          'G', 'T'],
            'Y' : [     'C',      'T'],
            'S' : [     'C', 'G'     ],
            'B' : [     'C', 'G', 'T'],
            'W' : ['A',           'T'],
            'R' : ['A',      'G'     ],
            'D' : ['A',      'G', 'T'],
            'M' : ['A', 'C'          ],
            'H' : ['A', 'C',      'T'],
            'V' : ['A', 'C', 'G'     ],
            'N' : ['A', 'C', 'G', 'T'],
            '-' : [                  ]
          },
          aliases={
            'U': 'T', # RNA.
            'u': 'T', # RNA.
            't': 'T',
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
            'n': 'N',
            '?': 'N',
            '.': '-'
          }
        )

    def get(cls):
        global dna

        if dna is None:
            dna = Dna()
        return dna
    get = classmethod(get)

#===============================================================================

# Forward declaration.
cdef class Protein(Character)

cdef Protein protein

cdef class Protein(Character):
    """
        Character codes for protein amino acids.
    """
    def __init__(self):
        global protein

        assert protein is None
        Character.__init__(self,
          states=[
            'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
            'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y', '*'
          ],
          ambiguities={
            'B': ['D', 'N'],
            'Z': ['E', 'Q'],
            'X': ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
                  'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y', '*'],
            '-': []
          },
          aliases={
            'a': 'A',
            'b': 'B',
            'c': 'C',
            'd': 'D',
            'e': 'E',
            'f': 'F',
            'g': 'G',
            'h': 'H',
            'i': 'I',
            'k': 'K',
            'l': 'L',
            'm': 'M',
            'n': 'N',
            'p': 'P',
            'q': 'Q',
            'r': 'R',
            's': 'S',
            't': 'T',
            'u': 'U',
            'v': 'V',
            'w': 'W',
            'x': 'X',
            'y': 'Y',
            'z': 'Z',
            '?': 'X',
            '.': '-'
          }
        )

    def get(cls):
        global protein

        if protein is None:
            protein = Protein()
        return protein
    get = classmethod(get)
