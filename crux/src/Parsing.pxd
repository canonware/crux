# Forward declarations.
cdef class Precedence
cdef class SymbolSpec
cdef class String
cdef class StringFirstSetCache
cdef class Symbol
cdef class Nonterm(Symbol)
cdef class NontermSpec(SymbolSpec)
cdef class Token(Symbol)
cdef class TokenSpec(SymbolSpec)
cdef class EndOfInput(Token)
cdef class Epsilon(Token)
cdef class Production
cdef class Start(Production)
cdef class Item
cdef class ItemSet
cdef class Action
cdef class ShiftAction(Action)
cdef class ReduceAction(Action)
cdef class Spec
cdef class Lr
cdef class Gsse
cdef class Gssn
cdef class Glr(Lr)

cdef class Precedence:
    cdef str name
    cdef str assoc
    cdef dict relationships
    cdef list equiv
    cdef list dominators

    cdef void specify(self, str assoc, dict relationships) except *

cdef class SymbolSpec:
    cdef str name
    cdef Precedence prec
    cdef list firstSet
    cdef list followSet
    cdef int seq

    cdef bint firstSetMerge(self, SymbolSpec sym)
    cdef bint followSetMerge(self, list set, TokenSpec epsilon)

cdef class String:
    cdef list rhs
    cdef int dotPos
    cdef SymbolSpec lookahead
    cdef int hash

    cdef int _hash(String self)

cdef class StringFirstSetCache:
    cdef TokenSpec _epsilon
    cdef dict _cache

    cdef dict getFirstSet(self, list rhs, int dotPos, SymbolSpec lookahead)

cdef class Symbol:
    cdef object __symSpec
    cdef Lr parser

cdef class Nonterm(Symbol):
    cpdef merge(self, Nonterm other)

cdef class NontermSpec(SymbolSpec):
    cdef str qualified
    cdef type nontermType
    cdef list chain # [%nonterm, %extends, ..., %extends]
    cdef list productions

cdef class Token(Symbol): pass

cdef class TokenSpec(SymbolSpec):
    cdef str qualified
    cdef type tokenType

cdef class EndOfInput(Token): pass

cdef class Epsilon(Token): pass

cdef class Production:
    cdef object method
    cdef str qualified
    cdef Precedence prec
    cdef NontermSpec lhs
    cdef list rhs
    cdef int seq

cdef class Start(Production): pass

cdef class Item:
    cdef Production production
    cdef int dotPos
    cdef dict lookahead
    cdef int hash

    cdef lr0__repr__(self)
    cdef void lookaheadInsert(self, SymbolSpec sym)
    cdef bint lookaheadDisjoint(self, Item other)

cdef class ItemSet:
    cdef dict _items
    cdef dict _added

    cdef void append(self, Item item)
    cdef bint addedAppend(self, Item item)
    cdef void _closeItems(self, list items, StringFirstSetCache cache)
    cdef void closure(self, StringFirstSetCache cache)
    cdef ItemSet xgoto(self, SymbolSpec sym)
    cdef bint merge(self, ItemSet other, StringFirstSetCache cache)
    cdef bint weakCompat(self, ItemSet other)

cdef class Action: pass

cdef class ShiftAction(Action):
    cdef int nextState

cdef class ReduceAction(Action):
    cdef Production production

cdef class Spec:
    cdef object _skinny
    cdef object _verbose
    cdef Precedence _none
    cdef Precedence _split
    cdef TokenSpec _eoi
    cdef TokenSpec _epsilon
    cdef dict _precedences
    cdef dict _nonterms
    cdef dict _tokens
    cdef dict _sym2spec
    cdef list _productions
    cdef NontermSpec _userStartSym
    cdef type _explicitStartSym
    cdef NontermSpec _startSym
    cdef Production _startProd
    cdef list _itemSets
    cdef dict _itemSetsHash
    cdef list _action
    cdef list _goto
    cdef int _startState
    cdef int _nActions
    cdef int _nConflicts
    cdef int _nImpure

    cdef void _prepare(self, list modules, str pickleFile, str pickleMode,
      str logFile, str graphFile) except *
    cdef void _introspect(self, list modules) except *
    cdef void _references(self, str logFile, str graphFile) except *
    cdef void _resolvePrec(self, str graphFile) except *
    cdef void _pickle(self, file_, mode)
    cdef _unpickle(self, file_, mode)
    cdef _compatible(self, Spec other)
    cdef void _validate(self, str logFile) except *
    cdef void _firstSets(self)
    cdef void _followSets(self)
    cdef void _items(self) except *
    cdef void _lr(self)
    cdef void _actionAppend(self, state, sym, action)
    cdef void _disambiguate(self)
    cdef str _resolve(self, SymbolSpec sym, Action oldAct, Action newAct)

cdef class Lr:
    cdef Spec _spec
    cdef bint _verbose
    cdef list _start
    cdef list _stack

    cpdef reset(self)
    cpdef token(self, Token token)
    cpdef eoi(self)
    cdef void _act(self, Symbol sym, SymbolSpec symSpec) except *
    cdef void _printStack(self)
    cdef void _reduce(self, Production production) except *
    cdef Symbol _production(self, Production production, list rhs)

cdef class Gsse:
    cdef Gssn node
    cdef Symbol value

cdef class Gssn:
    cdef list _edges
    cdef int nextState

cdef class Glr(Lr):
    cpdef reset(self)
    cpdef token(self, Token token)
    cpdef eoi(self)
    cdef void _act(self, Symbol sym, SymbolSpec symSpec) except *
    cdef void _reductions(self, Symbol sym, SymbolSpec symSpec) except *
    cdef void _reduceOne(self, list workQ, dict epsilons, object path, \
      Production production, SymbolSpec symSpec) except *
    cdef void _enqueueLimitedReductions(self, list workQ, dict epsilons, \
      Gsse edge, SymbolSpec symSpec) except *
    cdef void _shifts(self, Symbol sym, SymbolSpec symSpec) except *
    cdef void _printStack(self)
