cdef class Precedence:
    cdef readonly object name
    cdef readonly object assoc
    cdef readonly dict relationships
    cdef readonly list equiv
    cdef readonly list dominators

cdef class SymbolSpec:
    cdef readonly object name
    cdef public object prec
    cdef readonly list firstSet
    cdef readonly list followSet
    cdef readonly int seq
    cdef bint firstSetMerge(self, SymbolSpec sym)
    cdef bint followSetMerge(self, list set)

cdef class String:
    cdef list rhs
    cdef int dotPos
    cdef SymbolSpec lookahead
    cdef int hash
    cdef int _hash(String self)

cdef class Symbol:
    cdef object __symSpec
    cdef object __parser

cdef class Nonterm(Symbol):
    cpdef merge(self, Nonterm other)

cdef class NontermSpec(SymbolSpec):
    cdef readonly object qualified
    cdef readonly object nontermType
    cdef readonly list productions

cdef class Token(Symbol): pass

cdef class TokenSpec(SymbolSpec):
    cdef readonly object qualified
    cdef readonly object tokenType

cdef class EndOfInput(Token): pass

cdef class Epsilon(Token): pass

cdef class Production:
    cdef readonly object method
    cdef readonly object qualified
    cdef readonly Precedence prec
    cdef readonly NontermSpec lhs
    cdef readonly list rhs
    cdef readonly int seq

cdef class Start(Production): pass

cdef class Item:
    cdef Production __production
    cdef readonly int dotPos
    cdef readonly dict lookahead
    cdef readonly int hash
    cdef lr0__repr__(self)
    cdef void lookaheadInsert(self, SymbolSpec sym)
    cdef bint lookaheadDisjoint(self, Item other)

cdef class ItemSet:
    cdef readonly dict _items
    cdef dict _added
    cdef void append(self, Item item)
    cdef bint addedAppend(self, Item item)
    cdef void _closeItems(self, list items)
    cdef void closure(self)
    cdef ItemSet xgoto(self, SymbolSpec sym)
    cdef bint merge(self, ItemSet other)
    cdef bint weakCompat(self, ItemSet other)

cdef class Action: pass

cdef class ShiftAction(Action):
    cdef readonly int nextState

cdef class ReduceAction(Action):
    cdef readonly Production production

cdef class Spec:
    cdef object _skinny
    cdef object _verbose
    cdef Precedence _none
    cdef Precedence _split
    cdef dict _precedences
    cdef dict _nonterms
    cdef dict _tokens
    cdef readonly dict _sym2spec
    cdef list _productions
    cdef readonly NontermSpec _userStartSym
    cdef NontermSpec _startSym
    cdef Production _startProd
    cdef list _itemSets
    cdef dict _itemSetsHash
    cdef readonly list _action
    cdef readonly list _goto
    cdef int _startState
    cdef int _nActions
    cdef int _nConflicts
    cdef int _nImpure
    cdef void _prepare(self, modules, pickleFile, pickleMode, logFile,
      graphFile) except *
    cdef void _introspect(self, modules) except *
    cdef void _references(self, logFile, graphFile) except *
    cdef void _resolvePrec(self, graphFile)
    cdef void _pickle(self, file_, mode)
    cdef _unpickle(self, file_, mode)
    cdef _compatible(self, Spec other)
    cdef void _validate(self, logFile) except *
    cdef void _firstSets(self)
    cdef void _followSets(self)
    cdef void _items(self)
    cdef void _lr(self)
    cdef void _actionAppend(self, state, sym, action)
    cdef void _disambiguate(self)
    cdef _resolve(self, sym, oldAct, newAct)

cdef class Lr:
    cdef readonly Spec _spec
    cdef bint _verbose
    cdef list _start
    cdef readonly list _stack

cdef class Glr(Lr):
    cdef void _reduce(self, workQ, epsilons, path, production, symSpec) \
      except *
