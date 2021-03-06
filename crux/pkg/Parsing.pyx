#===============================================================================
# Copyright (c) 2007 Jason Evans <jasone@canonware.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#===============================================================================
#
# Release history:
#
# 1.3 (8 August 2007): Retroactively number public releases.
#
#                      Back-port to Python 2.4.
#
#                      Remove some magic surrounding epsilon, in order to
#                      generalize/simplify.
#
# 1.2 (6 May 2007): Fix some off-by-one errors in production count reporting.
#
#                   Add some missing code that helps detect which definitions
#                   are used/unused when building the parser.
#
# 1.1 (22 March 2007): Optimize/generalize Lr._production() by using argument
#                      list expansion.
#
# 1.0 (19 March 2007): Initial public release.
#
#===============================================================================
"""
    The Parsing module implements an LR(1) parser generator, as well as the
    runtime support for using a generated parser, via the Lr and Glr parser
    drivers.  There is no special parser generator input file format, but the
    parser generator still needs to know what classes/methods correspond to
    various aspects of the parser.  This information is specified via
    docstrings, which the parser generator introspects in order to generate a
    parser.  It is simplest to embed only one parser specification in each
    module, but it is possible to embed multiple parsers, split one or more
    parsers across multiple modules, and even nest parsers.

    The parsing tables are LR(1), but they are generated using a fast algorithm
    that avoids creating duplicate states that result when using the generic
    LR(1) algorithm.  Creation time and table size are on par with the LALR(1)
    algorithm.  However, LALR(1) can create reduce/reduce conflicts that don't
    exist in a true LR(1) parser.  For more information on the algorithm, see:

        A Practical General Method for Constructing LR(k) Parsers
        David Pager
        Acta Informatica 7, 249-268 (1977)

    Parsing table generation requires non-trivial amounts of time for large
    grammars.  Internal pickling support makes it possible to cache the most
    recent version of the parsing table on disk, and use the table if the
    current parser specification is still compatible with the one that was used
    to generate the pickled parsing table.  Since the compatibility checking is
    quite fast, even for large grammars, this removes the need to use the
    standard code generation method that is used by most parser generators.

    Parser specifications are encapsulated by the Spec class.  Parser instances
    use Spec instances, but are themselves based on separate classes.  This
    allows multiple parser instances to exist simultaneously, without requiring
    multiple copies of the parsing tables.  There are two separate parser
    driver classes:

      Lr : Standard Characteristic Finite State Machine (CFSM) driver, based on
           unambiguous LR(1) parsing tables.  This driver is faster than the
           Glr driver, but it cannot deal with all parsing tables that the Glr
           driver can.

      Glr : Generalized LR driver, capable of tracking multiple parse trees
            simultaneously, if the %split precedence is used to mark ambiguous
            actions.  This driver is closely based on Elkhound's design, which
            is described in a technical report:

                Elkhound: A Fast, Practical GLR Parser Generator
                Scott McPeak
                Report No. UCB/CSD-2-1214 (December 2002)
                http://www.cs.berkeley.edu/~smcpeak/elkhound/

    Parser generator directives are embedded in docstrings, and must begin with
    a '%' character, followed immediately by one of several keywords:

        Precedence : %fail %nonassoc %left %right %split
             Token : %token %extend
      Non-terminal : %start %nonterm %extend
        Production : %reduce %accept %amend %suppress

    All of these directives are associated with classes except for the
    production directives.  Productions are associated with methods within
    non-terminal classes.  The Parsing module provides base classes from which
    precedences, tokens, and non-terminals must be derived.  This is not as
    restrictive as it sounds, since there is nothing preventing, for example, a
    master Token class that subclasses Parsing.Token, which all of the actual
    token types then subclass.  Also, nothing prevents using multiple
    inheritance.

    Following are the base classes to be subclassed by parser specifications:

      * Precedence
      * Token
      * Nonterm

    The Parsing module implements the following exception classes:

      * Exception
      * AttributeError
      * SpecError
      * SyntaxError
"""
__all__ = ["Exception", "SpecError", "SyntaxError", "AttributeError",
           "Nonterm", "Parser", "Precedence", "Spec", "Token", "Lr", "Glr"]

import cPickle
import exceptions
import inspect
import re
import sys
import types

global __name__

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
cdef class _GssnEdgesIterHelper
cdef class _GssnNodesIterHelper
cdef class _GssnPathsIterHelper
cdef class Gssn
cdef class Glr(Lr)

#===============================================================================
# Begin exceptions.
#

class Exception(exceptions.Exception):
    """
        Top level Parsing exception class, from which all other Parsing
        exception classes inherit.
    """

class AttributeError(Exception, exceptions.AttributeError):
    """
        Attribute error, no different from the builtin exception, except that
        it also derives from Parsing.Exception.
    """
    def __init__(self, str):
        self._str = str

    def __str__(self):
        return self._str

class SpecError(Exception):
    """
        Specification error exception.  SpecError arises when the Spec
        introspection machinery detects an error either during docstring
        parsing or parser specification generation.
    """
    def __init__(self, str):
        self._str = str

    def __str__(self):
        return self._str

class SyntaxError(Exception, exceptions.SyntaxError):
    """
        Parser syntax error.  SyntaxError arises when a Parser instance detects
        a syntax error according to the Spec it is using, for the input being
        fed to it.
    """
    def __init__(self, str):
        self._str = str

    def __str__(self):
        return self._str

#
# End exceptions.
#===============================================================================

cdef class Precedence:
    """
        Precedences can be associated with tokens, non-terminals, and
        productions.  Precedence isn't as important for GLR parsers as for LR
        parsers, since GLR parsing allows for parse-time resolution of
        ambiguity.  Still, precedence can be useful for reducing the volume of
        ambiguities that must be dealt with at run-time.

        There are five precedence types: %fail, %nonassoc, %left, %right, and
        %split.  Each precedence can have relationships with other precedences:
        <, >, or =.  These relationships specify a directed acyclic graph
        (DAG), which is used to compute the transitive closures of
        relationships among precedences.  If no path exists between two
        precedences that are compared during conflict resolution, parser
        generation fails.  < and > are reflexive; it does not matter which is
        used.  Conceptually, the = relationship causes precedences to share a
        node in the DAG.

        During conflict resolution, an error results if no path exists in the
        DAG between the precedences under consideration.  When such a path
        exists, the highest precedence non-terminal or production takes
        precedence.  Associativity only comes into play for shift/reduce
        conflicts, where the terminal and the production have equivalent
        precedences (= relationship).  In this case, the non-terminal's
        associativity determines how the conflict is resolved.

        The %fail and %split associativities are special because they can be
        mixed with other associativities.  During conflict resolution, if
        another action has non-%fail associativity, then the %fail (lack of)
        associativity is overridden.  Similarly, %split associativity overrides
        any other associativity.  In contrast, any mixture of associativity
        between %nonassoc/%left/%right causes an unresolvable conflict.

               %fail : Any conflict is a parser-generation-time error.

                       A pre-defined precedence, [none], is provided.  It has
                       %fail associativity, and has no pre-defined precedence
                       relationships.

           %nonassoc : Resolve shift/reduce conflicts by removing both
                       possibilities, thus making conflicts a parse-time error.

               %left : Resolve shift/reduce conflicts by reducing.

              %right : Resolve shift/reduce conflicts by shifting.

              %split : Do not resolve conflicts; the GLR algorithm will split
                       the parse stack when necessary.

                       A pre-defined precedence, [split], is provided.  It has
                       %split associativity, and has no pre-defined precedence
                       relationships.

        By default, all symbols have [none] precedence.  Each production
        inherits the precedence of its left-hand-side nonterminal's precedence
        unless a precedence is manually specified for the production.

        Following are some examples of how to specify precedence classes:

          class P1(Parsing.Precedence):
              "%split p1"

          class p2(Parsing.Precedence):
              "%left" # Name implicitly same as class name.

          class P3(Parsing.Precedence):
              "%left p3 >p2" # No whitespace is allowed between > and p2.

          class P4(Parsing.Precedence):
              "%left p4 =p3" # No whitespace is allowed between = and p3.
    """
    def __init__(self, str name=None, str assoc=None, dict relationships=None):
        if name is None:
            return

        self.name = name
        self.assoc = None

        if assoc is not None:
            self.specify(assoc, relationships)

    def __reduce__(self):
        return (type(self), (), self.__getstate__())

    def __getstate__(self):
        return (self.name, self.assoc, self.relationships, self.equiv,
          self.dominators)

    def __setstate__(self, data):
        (self.name, self.assoc, self.relationships, self.equiv,
          self.dominators) = data

    def __repr__(self):
        cdef list equiv, domin
        cdef Precedence prec

        equiv = [prec.name for prec in self.equiv]
        equiv.sort()
        domin = [prec.name for prec in self.dominators]
        domin.sort()
        return "[%%%s %s ={%s} <{%s}]" % (self.assoc, self.name, \
          ",".join(equiv), ",".join(domin))

    # Important for pickling/unpickling.
    def __richcmp__(Precedence self, Precedence other, int op):
        assert op == 2
        return self is other

    cdef void specify(self, str assoc, dict relationships) except *:
        # We potentially create Precedence instances before their specifications
        # have been introspected, in which case only the name is known.
        # Furthermore, after introspection is complete, we need a way to
        # recognize instances that are referred to but never specified.  For
        # such instances, (self.assoc is None).
        assert self.assoc is None

        assert assoc in ("fail", "nonassoc", "left", "right", "split")

        self.assoc = assoc
        self.relationships = relationships # Raw relationships specification.

        self.equiv = [self] # Set.  Precedences that have equivalent precedence.
        self.dominators = [] # Set.  Precedences that have higher precedence.

cdef int _SymbolSpecSeq
_SymbolSpecSeq = 0

cdef class SymbolSpec:
    def __init__(self, str name=None, Precedence prec=None):
        global _SymbolSpecSeq

        if name is None:
            return

        self.name = name
        self.prec = prec
        self.chain = [self]
        self.firstSet = [] # Set.
        self.followSet = [] # Set.

        # Used for ordering symbols and hashing.
        self.seq = _SymbolSpecSeq
        _SymbolSpecSeq += 1

    def __reduce__(self):
        return (type(self), (), self.__getstate__())

    def __getstate__(self):
        return (self.name, self.prec, self.chain, self.firstSet,
          self.followSet, self.seq)

    def __setstate__(self, data):
        (self.name, self.prec, self.chain, self.firstSet, self.followSet,
          self.seq) = data

    def __repr__(self):
        return "%s" % self.name

    def __richcmp__(SymbolSpec self, SymbolSpec other, int op):
        cdef bint cmp

        cmp = self.seq == other.seq
        if op == 2: # ==
            return cmp
        elif op == 3: # !=
            return not cmp
        else:
            return NotImplemented

    cdef bint firstSetMerge(self, SymbolSpec sym):
        for 0 <= i < len(self.firstSet):
            elm = self.firstSet[i]
            if sym == elm:
                return True
        self.firstSet.append(sym)
        return False

    cdef bint followSetMerge(self, list set, TokenSpec epsilon):
        ret = True
        for 0 <= i < len(set):
            sym = set[i]
            if sym is not epsilon and sym not in self.followSet:
                self.followSet.append(sym)
                ret = False
        return ret

cdef class String:
    def __init__(self, list rhs, int dotPos, SymbolSpec lookahead):
        self.rhs = rhs
        self.dotPos = dotPos
        self.lookahead = lookahead

        if __debug__:
            for dotPos+1 <= i < len(self.rhs):
                sym = rhs[i]
                assert(isinstance(sym, SymbolSpec))

        self.hash = self._hash()

    # Conceptually, a String is represented as constructed in the syms property.
    # However, it is possible to avoid directly contstructing the syms list,
    # thus avoiding significant object construction/copying overhead.
    property syms:
        def __get__(self): return self.rhs[self.dotPos+1:] + [self.lookahead]

    def __richcmp__(String self, String other, int op):
        cdef int i, lenS
        cdef SymbolSpec symS, symO

        assert op == 2

        if self.hash != other.hash:
            return False

        lenS = len(self.rhs) - self.dotPos
        if lenS != len(other.rhs) - other.dotPos:
            return False

        for 1 <= i < lenS:
            symS = self.rhs[i + self.dotPos]
            symO = other.rhs[i + other.dotPos]
            if symS != symO:
                return False
        if self.lookahead != other.lookahead:
            return False

        return True

    def __hash__(self):
        return self.hash

    cdef int _hash(String self):
        cdef int ret, i
        cdef SymbolSpec sym

        ret = 5381
        for self.dotPos+1 <= i < len(self.rhs):
            sym = self.rhs[i]
            ret = hash(((ret << 5) + ret) + sym.seq)
        ret = hash(((ret << 5) + ret) + self.lookahead.seq)
        return ret

cdef class StringFirstSetCache:
    def __init__(self, TokenSpec epsilon):
        self._epsilon = epsilon
        self._cache = {}

    # Maintain a cache of String-->firstSet mappings, so that each followSet is
    # constructed only once.
    cdef dict getFirstSet(self, list rhs, int dotPos, SymbolSpec lookahead):
        cdef String s
        cdef dict firstSet
        cdef bint mergeEpsilon, hasEpsilon
        cdef SymbolSpec sym, elm

        s = String(rhs, dotPos, lookahead)
        if s in self._cache:
            return <dict>self._cache[s]
        else:
            # Calculate the first set for the string encoded by the s vector.
            firstSet = {} # Use dict rather than list during computation.
            mergeEpsilon = True
            for sym in s.syms:
                hasEpsilon = False
                for elm in sym.firstSet:
                    if elm is self._epsilon:
                        hasEpsilon = True
                    else:
                        firstSet[elm] = elm
                if not hasEpsilon:
                    mergeEpsilon = False
                    break
            # Merge epsilon if it was in the first set of every symbol.
            if mergeEpsilon:
                firstSet[self._epsilon] = self._epsilon

            # Cache result.
            self._cache[s] = firstSet
            return firstSet

cdef class Symbol:
    def __init__(self, symSpec, parser):
        self.__symSpec = symSpec
        self.parser = parser

    def __repr__(self):
        return "%r" % self.symSpec

    property symSpec:
        def __get__(self): return self.__symSpec

cdef class Nonterm(Symbol):
    """
        Non-terminal symbols have sets of productions associated with them.
        The productions induce a parse forest on an input token stream.  There
        is typically one special non-terminal, which is denoted via the %start
        directive.  If there is more than one %start directive, the
        Parsing.Spec constructor must be explicitly told which one to use as
        the start symbol.

        All other non-terminals are denoted via the %nonterm and %extend
        directives.  %nonterm is used to declare non-terminals with associated
        %reduce production methods.  %extend is used to modify existing
        non-terminals (including those declared via %start), with associated
        %reduce/%accept/%amend/%suppress production methods.  %extend can be
        used to create arbitrarily long chains of non-terminal modifications,
        but no forks are allowed.

        In addition to production methods, the merge() method may be called
        during resolution of ambiguous parses.  See the merge() documentation
        for further details.

        Following are examples of how to specify non-terminal classes and their
        associated productions:

          class E(Parsing.Nonterm):
              "%start E"
              def __init__(self):
                  Parsing.Nonterm.__init__(self)
                  # ...

              # Productions.
              def reduceA(self, E, plus, T):
                  "%reduce E plus T [split]"
                  print "%r ::= %r %r %r." % (self, E, plus, T)

              def reduceB(self, T):
                  "%reduce T"

          class T(Parsing.Nonterm):
              "%nonterm" # Name implicitly same as class name.
              def reduceA(self, T, star, F):
                  "%reduce T star F"

              def reduceB(self, F):
                  "%reduce F [p1]"

          class F(Parsing.Nonterm):
              "%nonterm F [p2]"
              def reduceA(self, lparen, E, rparen):
                  "%reduce lparen E rparen"

              def reduceB(self, id):
                  "%reduce id"

              def reduceC(self, x):
                  "%reduce x"

          # Extend F.  Both subclassing F and the %extend directive are
          # required.
          class Fsub(F):
              "%extend F [p3]"
              def reduceA(self, lparen, E, rparen):
                  "%accept"
                  # Leave the production associated with F.reduceA intact, but
                  # provide a different method to call when the production
                  # accepts.

              def reduceB(self, x):
                  "%amend x"
                  # Replace F.reduceB.

          # Extend F again.
          class Fsubsub(Fsub):
              "%extend F"
              def reduceB(self, id):
                  "%reduce id"

              def reduceC(self, x):
                  "%suppress"
                  # Do not use F.reduceC.
    """
    def __init__(self, Lr parser):
        Symbol.__init__(self, parser.spec._sym2spec[type(self)], parser)

    cpdef merge(self, Nonterm other):
        """
            Merging happens when there is an ambiguity in the input that allows
            non-terminals to be part of multiple overlapping series of
            reductions.  If no merge() method is specified, the parser will
            raise a SyntaxError upon encountering an ambiguity that confounds
            reduction processing.  However, it may be useful to either discard
            one of the possible parses, or to explicitly record the ambiguity
            in the data structures being created during parsing.  In both of
            these cases, the non-terminal-specific merge() is the place to do
            the work; merge() returns an object that is stored by the parser
            onto the parse stack.  In the case where merge() discards one of
            the possible parses, it need only return the parse that is to be
            preserved (self or other).

            If multiple merges are necessary, they cause a series of merge()
            calls.  The first alternative (self) may be the result of a
            previous merge() call, whereas other will not have not been merged
            yet (unless as the result of merging further down in the parse
            forest).

            The alternative that is discarded is never touched by the parser
            again, so if any immediate cleanup is necessary, it should be done
            in merge().
        """
        raise SyntaxError("No merge() for %r; merging %r <--> %r" % \
          (type(self), self, other))

cdef class NontermSpec(SymbolSpec):
    def __init__(self, nontermType=None, name=None, qualified=None,
      Precedence prec=None):
        if nontermType is None:
            return

        SymbolSpec.__init__(self, name, prec)

        self.qualified = qualified
        self.nontermType = nontermType
        self.productions = [] # Set.

    def __reduce__(self):
        return (type(self), (), self.__getstate__())

    def __getstate__(self):
        return (SymbolSpec.__getstate__(self), self.qualified, self.productions)

    def __setstate__(self, data):
        cdef SymbolSpec_data
        cdef str qualified
        cdef list productions

        (SymbolSpec_data, qualified, productions) = data

        SymbolSpec.__setstate__(self, SymbolSpec_data)

        # Convert qualified name to a type reference.
        elms = qualified.split(".")
        nontermType = sys.modules[elms[0]]
        for elm in elms[1:]:
            nontermType = nontermType.__dict__[elm]

        (self.qualified, self.nontermType, self.productions) = \
          (qualified, nontermType, productions)

    def __hash__(self):
        return hash(self.qualified)

cdef class Token(Symbol):
    """
        Tokens are terminal symbols.  The parser is fed Token instances, which
        is what drives parsing.  Typically, the user will define a class that
        subclasses Parsing.Token and implement parser-specific machinery there,
        then derive all actual token types from that class.

        The %extend directive can be used to modify existing tokens, in a
        similar fashion to how non-terminals can be modified.  The primary
        reason to extend a token is to add some behavior during token creation
        in an extended parser.  Note, however, that the extended token type
        must be used when feeding tokens to the parser, so the original
        tokenizer must either be replaced or somehow modified for extended
        tokens to work properly.

          class Token(Parsing.Token):
              def __init__(self, Lr parser):
                  Parsing.Token.__init__(self, parser)
                  # ...

          class Plus(Token):
              "%token plus [p1]"

          class star(Token):
              "%token star [p2]" # Name implicitly same as class name.

          class lparen(Token):
              "%token [split]"

          class rparen(Token):
              "%token [none]" # [none] not necessary, since it's the default.

          class id(Token):
              "%token"
    """
    def __init__(self, Lr parser):
        Symbol.__init__(self, parser.spec._sym2spec[type(self)], parser)

# AKA terminal symbol.
cdef class TokenSpec(SymbolSpec):
    def __init__(self, type tokenType=None, str name=None, str qualified=None,
      Precedence prec=None):
        if tokenType is None:
            return

        assert issubclass(tokenType, Token)

        SymbolSpec.__init__(self, name, prec)

        self.qualified = qualified
        self.tokenType = tokenType

    def __reduce__(self):
        return (type(self), (), self.__getstate__())

    def __getstate__(self):
        return (SymbolSpec.__getstate__(self), self.qualified)

    def __setstate__(self, data):
        cdef SymbolSpec_data
        cdef str qualified

        (SymbolSpec_data, qualified) = data

        SymbolSpec.__setstate__(self, SymbolSpec_data)

        # Convert qualified name to a type reference.
        elms = qualified.split(".")
        tokenType = sys.modules[elms[0]]
        for elm in elms[1:]:
            tokenType = tokenType.__dict__[elm]

        self.tokenType = tokenType

    def __hash__(self):
        return hash(self.qualified)

# <$>.
cdef class EndOfInput(Token): pass

# <e>.
cdef class Epsilon(Token): pass

cdef int _ProductionSeq
_ProductionSeq = 0

cdef class Production:
    def __init__(self, method=None, qualified=None, Precedence prec=None,
      lhs=None, rhs=None):
        global _ProductionSeq

        if (method is None):
            return

        if __debug__:
            for elm in rhs:
                assert isinstance(elm, SymbolSpec)

        self.method = method
        self.qualified = qualified
        self.prec = prec
        self.lhs = lhs
        self.rhs = rhs

        # Used for hashing.
        self.seq = _ProductionSeq
        _ProductionSeq += 1

    def __reduce__(self):
        return (type(self), (), self.__getstate__())

    def __getstate__(self):
        return (self.qualified, self.prec, self.lhs, self.rhs, self.seq)

    def __setstate__(self, data):
        cdef str qualified
        cdef Precedence prec
        cdef NontermSpec lhs
        cdef list rhs
        cdef int seq

        # Convert qualified name to a function reference.
        (qualified, prec, lhs, rhs, seq) = data
        elms = qualified.split(".")
        method = sys.modules[elms[0]]
        for elm in elms[1:]:
            method = method.__dict__[elm]

        # Set state.
        self.method = method
        self.qualified = qualified
        self.prec = prec
        self.lhs = lhs
        self.rhs = rhs
        self.seq = seq

    def __repr__(self):
        cdef SymbolSpec elm

        return "%r ::= %s. [%s]" % \
          (self.lhs, " ".join(["%r" % elm for elm in self.rhs]), self.prec.name)

    # Optional callback method.
    #
    # Called when a production is reduced.
    def reduce(self, lhs, *rhs): pass

class NontermStart(Nonterm):
    def reduce(self, Nonterm userStartSym, EndOfInput eoi):
        pass

cdef class Start(Production):
    def __init__(self, startSym, userStartSym):
        Production.__init__(self, None, startSym, userStartSym)

cdef class Item:
    def __init__(self, Production production=None, int dotPos=0,
      dict lookahead=None):
        if production is None:
            return

        assert dotPos >= 0
        assert dotPos <= len(production.rhs)
        if __debug__:
            for elm in lookahead:
                assert isinstance(elm, SymbolSpec)
                assert lookahead[elm] == elm

        self.production = production
        self.dotPos = dotPos
        self.lookahead = dict(lookahead)

        self.hash = (dotPos * _ProductionSeq) + production.seq

    def __reduce__(self):
        return (type(self), (), self.__getstate__())

    def __getstate__(self):
        return (self.production, self.dotPos, self.lookahead, self.hash)

    def __setstate__(self, data):
        (self.production, self.dotPos, self.lookahead, self.hash) = data

    def __hash__(self):
        return self.hash

    def __richcmp__(Item self, Item other, int op):
        if op == 2: # ==
            return self.hash == other.hash
        else:
            assert op == 0 # <
            return self.hash < other.hash

    def __repr__(self):
        cdef list strs, deco
        cdef object elm
        cdef int i
        cdef SymbolSpec sym

        strs = []
        strs.append("[%r ::=" % self.production.lhs)
        assert self.dotPos <= len(self.production.rhs)
        i = 0
        while i < self.dotPos:
            strs.append(" %r" % self.production.rhs[i])
            i += 1
        strs.append(" *")
        while i < len(self.production.rhs):
            strs.append(" %r" % self.production.rhs[i])
            i += 1
        deco = [(sym.name, sym) for sym in self.lookahead.iterkeys()]
        deco.sort()
        strs.append("., %s] [%s]" % \
          ("/".join(["%r" % elm[1] for elm in deco]), \
          self.production.prec.name))

        return "".join(strs)

    cdef lr0__repr__(self):
        cdef list strs
        cdef int i

        strs = []
        strs.append("%r ::=" % self.production.lhs)
        assert self.dotPos <= len(self.production.rhs)
        i = 0
        while i < self.dotPos:
            strs.append(" %r" % self.production.rhs[i])
            i += 1
        strs.append(" *")
        while i < len(self.production.rhs):
            strs.append(" %r" % self.production.rhs[i])
            i += 1
        strs.append(". [%s]" % self.production.prec.name)

        return "".join(strs)

    cdef void lookaheadInsert(self, SymbolSpec sym):
        self.lookahead[sym] = sym

    cdef bint lookaheadDisjoint(self, Item other):
        cdef dict sLookahead, oLookahead
        cdef SymbolSpec sSym, oSym
        cdef list keys
        cdef int i

        sLookahead = self.lookahead
        oLookahead = other.lookahead

        keys = sLookahead.keys()
        for 0 <= i < len(keys):
            sSym = keys[i]
            if sSym in oLookahead:
                return False

        keys = oLookahead.keys()
        for 0 <= i < len(keys):
            oSym = keys[i]
            if oSym in sLookahead:
                return False

        return True

cdef class _ItemSetIterHelper:
    cdef _items
    cdef int _index

    def __init__(self, dict items, dict added):
        cdef Item item

        if __debug__:
            for item in items.iterkeys():
                assert item.production.lhs.name == "<S>" or item.dotPos != 0
            for item in added.iterkeys():
                assert item.dotPos == 0
                assert item.production.lhs.name != "<S>"
        self._items = items.keys() + added.keys()
        self._index = 0

    def __iter__(self):
        return self

    def __next__(self):
        cdef Item ret

        if self._index >= len(self._items):
            raise StopIteration
        ret = self._items[self._index]
        self._index += 1
        return ret

cdef class ItemSet:
    def __init__(self):
        self._items = {}
        self._added = {}

    def __reduce__(self):
        return (type(self), (), self.__getstate__())

    def __getstate__(self):
        return (self._items, self._added)

    def __setstate__(self, data):
        (self._items, self._added) = data

    def __len__(self):
        return len(self._items)

    def __repr__(self):
        cdef list kernel, added
        cdef Item item

        kernel = [item for item in self._items.iterkeys()]
        kernel.sort()
        added = [item for item in self._added.iterkeys()]
        added.sort()
        return "ItemSet(kernel: %s, added: %r)" % \
          (", ".join(["%r" % item for item in kernel]), \
          ", ".join(["%r" % item for item in added]))

    def __hash__(self):
        cdef int ret
        cdef Item item

        # This works because addition is transitive.
        ret = 0
        for item in self._items.iterkeys():
            ret += item.hash
        ret = hash(ret)
        return ret

    def __richcmp__(ItemSet self, ItemSet other, int op):
        cdef Item sItem

        assert op == 2

        if len(self._items) != len(other._items):
            return False

        for sItem in self._items.iterkeys():
            if sItem not in other._items and sItem not in other._added:
                return False

        return True

    def __iter__(self):
        return _ItemSetIterHelper(self._items, self._added)

    # Merge a kernel item.
    cdef void append(self, Item item):
        cdef Item tItem

        assert item.production.lhs.name == "<S>" or item.dotPos != 0

        if item in self._items or item in self._added:
            (<Item>self._items[item]).lookahead.update(item.lookahead)
        else:
            tItem = Item(item.production, item.dotPos, item.lookahead)
            self._items[tItem] = tItem

    # Merge an added item.
    cdef bint addedAppend(self, Item item):
        cdef dict lookahead
        cdef int oldLen

        assert item.dotPos == 0
        assert item.production.lhs.name != "<S>"

        if item in self._added:
            lookahead = (<Item>self._added[item]).lookahead
            oldLen = len(lookahead)
            lookahead.update(item.lookahead)
            return (oldLen != len(lookahead))
        else:
            self._added[item] = item
            return True

    # Given a list of items, compute their closure and merge the results into
    # the set of added items.
    cdef void _closeItems(self, list items, StringFirstSetCache cache):
        cdef int i, dotPos
        cdef list rhs
        cdef Item item, tItem
        cdef SymbolSpec lookahead
        cdef NontermSpec lhs
        cdef dict firstSet

        # Iterate over the items until no more can be added to the closure.
        i = 0
        while i < len(items):
            item = items[i]
            rhs = item.production.rhs
            dotPos = item.dotPos
            if dotPos < len(rhs) \
              and isinstance(rhs[dotPos], NontermSpec):
                for lookahead in item.lookahead.keys():
                    firstSet = cache.getFirstSet(rhs, dotPos, lookahead)
                    lhs = rhs[dotPos]
                    for prod in lhs.productions:
                        tItem = Item(prod, 0, firstSet)
                        if self.addedAppend(tItem):
                            items.append(tItem)
            i += 1

    # Calculate and merge the kernel's transitive closure.
    cdef void closure(self, StringFirstSetCache cache):
        cdef list items, rhs
        cdef Item item, tItem
        cdef int dotPos
        cdef SymbolSpec lookahead
        cdef NontermSpec lhs
        cdef dict firstSet
        cdef int i

        items = []
        for item in self._items:
            rhs = item.production.rhs
            dotPos = item.dotPos
            if dotPos < len(rhs) and isinstance(rhs[dotPos], NontermSpec):
                for lookahead in item.lookahead:
                    firstSet = cache.getFirstSet(rhs, dotPos, lookahead)
                    lhs = rhs[dotPos]
                    for 0 <= i < len(lhs.productions):
                        prod = lhs.productions[i]
                        tItem = Item(prod, 0, firstSet)
                        if self.addedAppend(tItem):
                            items.append(tItem)
        self._closeItems(items, cache)

    # Calculate the kernel of the goto set, given a particular symbol.
    cdef ItemSet xgoto(self, SymbolSpec sym):
        cdef ItemSet ret
        cdef Item item, tItem
        cdef list items, rhs
        cdef Production production
        cdef int dotPos, i, iLim

        ret = ItemSet()
        items = self._items.keys()
        items.extend(self._added.keys())
        iLim = len(items)
        for 0 <= i < iLim:
            item = <Item>items[i]
            production = item.production
            rhs = production.rhs
            dotPos = item.dotPos
            if dotPos < len(rhs) and <SymbolSpec>rhs[dotPos] == sym:
                tItem = Item(item.production, dotPos + 1, item.lookahead)
                ret.append(tItem)
        return ret

    # Merge the kernel of other into this ItemSet, then update the closure.
    # It is not sufficient to copy other's added items, since other has not
    # computed its closure.
    cdef bint merge(self, ItemSet other, StringFirstSetCache cache):
        cdef list items
        cdef Item item, tItem
        cdef dict lookahead, tLookahead
        cdef SymbolSpec sym

        items = []
        for item in other._items.iterkeys():
            if item in self._items or item in self._added:
                lookahead = (<Item>self._items[item]).lookahead
                tLookahead = {}
                for sym in item.lookahead.iterkeys():
                    if sym not in lookahead:
                        lookahead[sym] = sym
                        tLookahead[sym] = sym
                if len(tLookahead) > 0:
                    tItem = Item(item.production, item.dotPos, tLookahead)
                    items.append(tItem)
            else:
                tItem = Item(item.production, item.dotPos, item.lookahead)
                self._items[tItem] = tItem
                items.append(tItem)

        if len(items) > 0:
            self._closeItems(items, cache)
            return True
        else:
            return False

    # Determine if self and other are weakly compatible, as defined by the
    # Pager(1977) algorithm.
    cdef bint weakCompat(self, ItemSet other):
        cdef int i, j
        cdef Item isItem, ioItem, jsItem, joItem

        # Check for identical kernel LR(0) items, and pair items, for later use.
        if len(self) != len(other):
            return False
        pairs = []
        for sItem in self._items.iterkeys():
            if sItem not in other._items:
                return False
            oItem = other._items[sItem]
            pairs.append((sItem, oItem))

        # Check for lookahead compatibility.
        for 0 <= i < len(pairs)-1:
            iPair = pairs[i]
            isItem = iPair[0]
            ioItem = iPair[1]
            for i+1 <= j < len(pairs):
                jPair = pairs[j]
                jsItem = jPair[0]
                joItem = jPair[1]

                if isItem.lookaheadDisjoint(joItem) \
                  and ioItem.lookaheadDisjoint(jsItem):
                    pass
                elif not isItem.lookaheadDisjoint(jsItem):
                    pass
                elif not ioItem.lookaheadDisjoint(joItem):
                    pass
                else:
                    return False
        return True

cdef class Action:
    """
        Abstract base class, subclassed by {Shift,Reduce}Action.
    """
    def __init__(self): pass

    def __reduce__(self):
        return (type(self), (), self.__getstate__())

    def __getstate__(self):
        return ()

    def __setstate__(self, data):
        pass

cdef class ShiftAction(Action):
    """
        Shift action, with assocated nextState.
    """
    def __init__(self, nextState=None):
        if nextState is None:
            return

        Action.__init__(self)
        self.nextState = nextState

    def __reduce__(self):
        return (type(self), (), self.__getstate__())

    def __getstate__(self):
        return (Action.__getstate__(self), self.nextState)

    def __setstate__(self, data):
        cdef Action_data
        cdef int nextState

        (Action_data, nextState) = data
        Action.__setstate__(self, Action_data)
        self.nextState = nextState

    def __repr__(self):
        return "[shift %r]" % self.nextState

    def __richcmp__(ShiftAction self, other, int op):
        cdef bint equal

        if not isinstance(other, ShiftAction) or \
          self.nextState != (<ShiftAction>other).nextState:
            equal = False
        else:
            equal = True

        if op == 2: # ==
            return equal
        else:
            assert op == 3 # !=
            return not equal

cdef class ReduceAction(Action):
    """
        Reduce action, with associated production.
    """
    def __init__(self, production=None):
        if production is None:
            return

        Action.__init__(self)
        self.production = production

    def __reduce__(self):
        return (type(self), (), self.__getstate__())

    def __getstate__(self):
        return (Action.__getstate__(self), self.production)

    def __setstate__(self, data):
        cdef Action_data
        cdef Production production

        (Action_data, production) = data
        Action.__setstate__(self, Action_data)
        self.production = production

    def __repr__(self):
        return "[reduce %r]" % self.production

    def __richcmp__(ReduceAction self, other, int op):
        cdef bint equal

        if not isinstance(other, ReduceAction) or \
          self.production != (<ReduceAction>other).production:
            equal = False
        else:
            equal = True

        if op == 2: # ==
            return equal
        else:
            assert op == 3 # !=
            return not equal

cdef class Spec:
    """
        The Spec class contains the read-only data structures that the Parser
        class needs in order to parse input.  Parser generation results in a
        Spec instance, which can then be shared by multiple Parser instances.

        Constructor argument documentation:
        =====================================================================
        modules : Either a single module, or a list of modules, wherein to
                  look for parser generator directives in docstrings.

        startSym : Start symbol.  This argument must be specified when there
                   are multiple non-terminals that use the %start directive.
                   This can happen if there are multiple parser
                   specifications embedded in the same set of modules, or if
                   parser specifications are nested.

        pickleFile : The path of a file to use for Spec pickling/unpickling.
                     The file will be created if necessary, but the
                     containing directory must already exist.

        pickleMode :  "r" : Unpickle from pickleFile.
                      "w" : Pickle to pickleFile.
                      "rw" : Unpickle/pickle from/to pickleFile.

        skinny : If true, discard all data that are only strictly necessary
                 while constructing the parsing tables.  This reduces
                 available debugging context, but substantially reduces
                 pickle size.

        logFile : The path of a file to store a human-readable copy of the
                  parsing tables in.  The file will be created if necessary,
                  but the containing directory must already exist.

        graphFile : The path of a file to store a graphviz representation
                    (dot format) of the precedence relationship graph.  The
                    file will be created if necessary, but the containing
                    directory must already exist.

        verbose : If true, print progress information while generating the
                  parsing tables.
        =====================================================================
    """
    def __init__(self, modules=None, type startSym=None,
      str pickleFile=None, str pickleMode="rw", bint skinny=True,
      str logFile=None, str graphFile=None, bint verbose=False):
        if modules is None:
            return

        assert pickleMode in ("rw", "r", "w")
        assert startSym is None or issubclass(startSym, Nonterm)

        if startSym is not None and issubclass(startSym, Nonterm):
            if type(startSym.__doc__) != str \
              or startSym.__doc__.split("\n")[-1].split(" ")[0] != "%start":
                raise SpecError("Invalid start symbol %r: %r" % \
                  (startSym, startSym.__doc__))

        self._skinny = skinny
        self._verbose = verbose

        # Default (no) precedence.
        self._none = Precedence("none", "fail", {})
        self._split = Precedence("split", "split", {})

        # Special tokens.
        self._eoi = TokenSpec(EndOfInput, "<$>", "%s.EndOfInput" % __name__, \
          self._none)
        self._epsilon = TokenSpec(Epsilon, "<e>", "%s.Epsilon" % __name__, \
          self._none)

        # Symbols are maintained as two separate sets so that non-terminals and
        # terminals (tokens) can be operated on separately where needed.
        self._precedences = {self._none.name: self._none,
                             self._split.name: self._split}
        self._nonterms = {}

        self._tokens = {self._eoi.name: self._eoi, \
          self._epsilon.name: self._epsilon}
        self._sym2spec = {EndOfInput: self._eoi, Epsilon: self._epsilon}
        self._productions = []

        self._userStartSym = None
        self._explicitStartSym = startSym
        self._startSym = None
        self._startProd = None

        # Everything below this point is computed from the above (once
        # introspection is complete).

        self._itemSets = [] # Each element corresponds to an element in _action.
        self._itemSetsHash = None
        # LR parsing tables.  The tables conceptually contain one state per row,
        # where each row contains one element per symbol.  The table is
        # conceptually in row-major form, but each row is actually a dictionary.
        # If no entry for a symbol exists for a particular state, then input of
        # that symbol is an error for that state.
        self._action = []
        self._goto = []
        self._startState = -1
        self._nActions = 0
        self._nConflicts = 0
        self._nImpure = 0 # Number of LR impurities (does not affect GLR).

        # Introspect modules and generate parse tables.
        if type(modules) == types.ModuleType:
            # Wrap single module in a list.
            modules = [modules]
        self._prepare(modules, pickleFile, pickleMode, logFile, graphFile)

    def __reduce__(self):
        return (type(self), (), self.__getstate__())

    def __getstate__(self):
        return (self._skinny, self._verbose, self._none, self._split,
          self._eoi, self._epsilon, self._precedences, self._nonterms,
          self._tokens, self._sym2spec, self._productions, self._userStartSym,
          self._startSym, self._startProd, self._itemSets, self._itemSetsHash,
          self._action, self._goto, self._startState, self._nActions,
          self._nConflicts, self._nImpure)

    def __setstate__(self, data):
        (self._skinny, self._verbose, self._none, self._split, self._eoi,
          self._epsilon, self._precedences, self._nonterms, self._tokens,
          self._sym2spec, self._productions, self._userStartSym,
          self._startSym, self._startProd, self._itemSets, self._itemSetsHash,
          self._action, self._goto, self._startState, self._nActions,
          self._nConflicts, self._nImpure) = data

    property pureLR:
        def __get__(self): return (self._nConflicts + self._nImpure == 0)

    property conflicts:
        def __get__(self): return self._nConflicts

    def __repr__(self):
        cdef str ret, conflict, resolve
        cdef list lines, deco, deco2, items
        cdef object elm, elm2
        cdef Precedence prec
        cdef SymbolSpec sym
        cdef TokenSpec tokenSpec
        cdef NontermSpec nontermSpec
        cdef Production prod
        cdef int i, ntokens, nnonterms, nproductions, nstates
        cdef Item item
        cdef Action action, other

        if self._skinny:
            # Print a very reduced summary, since most info has been discarded.
            return "Parsing.Spec: %d states, %d actions (%d split)" % \
              (len(self._action), self._nActions, self._nImpure)

        lines = []

        #=======================================================================
        lines.append("Precedences:")
        deco = [(prec.name, prec) for prec in self._precedences.itervalues()]
        deco.sort()
        for elm in deco:
            prec = elm[1]
            lines.append("  %r" % prec)

        lines.append("Tokens:")
        deco = [(sym.name, sym) for sym in self._tokens.itervalues()]
        deco.sort()
        for elm in deco:
            tokenSpec = <TokenSpec>elm[1]
            lines.append("  %r %r" % (tokenSpec, tokenSpec.prec))
            deco2 = [(sym.name, sym) for sym in tokenSpec.firstSet]
            deco2.sort()
            lines.append("    First set: %r" % [elm2[1] for elm2 in deco2])
            deco2 = [(sym.name, sym) for sym in tokenSpec.followSet]
            deco2.sort()
            lines.append("    Follow set: %r" % [elm2[1] for elm2 in deco2])

        lines.append("Non-terminals:")
        deco = [(sym.name, sym) for sym in self._nonterms.itervalues()]
        deco.sort()
        for elm in deco:
            nontermSpec = <NontermSpec>elm[1]
            lines.append("  %r %r" % (nontermSpec, nontermSpec.prec))
            deco2 = [(sym.name, sym) for sym in nontermSpec.firstSet]
            deco2.sort()
            lines.append("    First set: %r" % [elm2[1] for elm2 in deco2])
            deco2 = [(sym.name, sym) for sym in nontermSpec.followSet]
            deco2.sort()
            lines.append("    Follow set: %r" % [elm2[1] for elm2 in deco2])
            lines.append("    Productions:")
            deco2 = [([sym.name for sym in prod.rhs], prod) \
              for prod in nontermSpec.productions]
            deco2.sort()
            for elm2 in deco2:
                prod = elm2[1]
                lines.append("      %r" % prod)

        lines.append("Item sets:")
        for 0 <= i < len(self._itemSets):
            lines.append("  %d: %r" % (i, self._itemSets[i]))
        #=======================================================================

        ntokens = len(self._tokens) - 2
        nnonterms = len(self._nonterms) - 1
        nproductions = len(self._productions) - 1
        nstates = len(self._action)
        lines.append(("Parsing.Spec: %d token%s, %d non-terminal%s, " + \
          "%d production%s, %d state%s, %d action%s (%d split):") % \
          (ntokens, ("s", "")[ntokens == 1], \
          nnonterms, ("s", "")[nnonterms == 1], \
          nproductions, ("s", "")[nproductions == 1], \
          nstates, ("s", "")[nstates == 1], \
          self._nActions, ("s", "")[self._nActions == 1], \
          self._nImpure))
        if self.pureLR:
            lines.append("Algorithm compatibility: GLR, LR")
        elif self._nConflicts == 0:
            lines.append("Algorithm compatibility: GLR")
        else:
            lines.append("Algorithm compatibility: None, due to ambiguity")
        lines.append("Parsing tables:")
        for 0 <= i < len(self._action):
            lines.append("  %s" % ("=" * 78))
            lines.append("  State %d:%s" % \
              (i, ("", " (start state)")[self._startState == i]))
            items = [item for item in self._itemSets[i]]
            items.sort()
            for item in items:
                lines.append(" %s%s" % (" " * (len("%d" % i) + 9),
                  item.lr0__repr__()))
            lines.append("    Goto:")
            deco = [(sym.name, sym) for sym in self._goto[i]]
            deco.sort()
            for elm in deco:
                sym = elm[1]
                lines.append("    %15r : %r" % (sym, self._goto[i][sym]))
            lines.append("    Action:")
            deco = [(sym.name, sym) for sym in self._action[i]]
            deco.sort()
            for elm in deco:
                sym = elm[1]
                for action in self._action[i][sym]:
                    conflict = "   "
                    for other in self._action[i][sym]:
                        if action != other:
                            resolution = self._resolve(sym, other, action)
                            if resolution == "err":
                                conflict = "XXX"
                                break

                    if type(action) == ShiftAction:
                        lines.append("%s %15r : %-6s %d [%s]" % \
                          (conflict, sym, "shift", \
                          (<ShiftAction>action).nextState, sym.prec.name))
                    else:
                        assert type(action) == ReduceAction
                        lines.append("%s %15r : %-6s %r" % \
                          (conflict, sym, "reduce", \
                          (<ReduceAction>action).production))

        ret = "\n".join(lines)
        return ret

    cdef void _prepare(self, list modules, str pickleFile, str pickleMode,
      str logFile, str graphFile) except *:
        """
            Compile the specification into data structures that can be used by
            the Parser class for parsing.
        """
        cdef int ntokens, nnonterms, nproductions
        cdef str compat

        # Get the grammar specification.
        self._introspect(modules)

        # Resolve references in the grammar specification.
        self._references(logFile, graphFile)

        # Augment grammar with a special start symbol and production:
        #
        #   <S> ::= S <$>.
        assert self._startSym is None
        assert isinstance(self._userStartSym, NontermSpec)
        self._startSym = NontermSpec(NontermStart, "<S>",
          "%s.NontermStart" % __name__, self._none)
        self._startProd = Production(NontermStart.reduce.im_func,
                                     "%s.NontermStart.reduce" % __name__,
                                     self._none, self._startSym, \
                                     [self._userStartSym, self._eoi])
        self._startSym.productions.append(self._startProd)
        self._nonterms["<S>"] = self._startSym
        self._productions.append(self._startProd)

        if self._verbose:
            ntokens = len(self._tokens) - 2
            nnonterms = len(self._nonterms) - 1
            nproductions = len(self._productions) - 1
            print \
              "Parsing.Spec: %d token%s, %d non-terminal%s, %d production%s" % \
              (ntokens, ("s", "")[ntokens == 1], \
              nnonterms, ("s", "")[nnonterms == 1], \
              nproductions, ("s", "")[nproductions == 1])

        # Check for a compatible pickle.
        compat = self._unpickle(pickleFile, pickleMode)

        if compat == "incompatible":
            # Create the collection of sets of LR(1) items.
            self._firstSets()
            self._followSets()
            self._items()

        if compat == "compatible":
            # Just because the pickle was compatible does not mean that it is
            # valid for parsing.
            if self._nConflicts != 0:
                raise SpecError( \
                  "Compatible pickle is invalid due to conflicts (%d)" % \
                  self._nConflicts)
        if compat in ["itemsets", "incompatible"]:
            # Generate LR(1) parsing tables.
            self._lr()

            # Disambiguate actions.
            self._disambiguate()

            # Check for unused or ambiguous definitions, as well as reporting
            # ambiguities.
            try:
                self._validate(logFile)
            finally:
                # Pickle the spec, if method parameters so dictate, even if
                # there were validation errors, so that the pickle might be
                # used in part during later runs.
                self._pickle(pickleFile, pickleMode)
        elif compat == "repickle":
            # Pickle the spec, if method parameters so dictate.
            self._pickle(pickleFile, pickleMode)

        if self._skinny:
            # Discard bulky data that are not needed during parsing.  Note that
            # _pickle() also discarded data that don't even need to be pickled.
            #self._none = ...
            #self._split = ...
            self._precedences = {}
            self._nonterms = {}
            self._tokens = {}
            self._productions = []

    # Introspect modules and find special parser declarations.  In order to be
    # a special class, the class must both 1) be subclassed from Token or
    # Nonterm, and 2) contain the appropriate %foo docstring.
    cdef void _introspect(self, list modules) except *:
        cdef list deferTokens, deferNonterms
        cdef object module, d, v
        cdef str k
        cdef list dirtoks
        cdef str name
        cdef dict relationships
        cdef int i
        cdef str tok
        cdef object m
        cdef Precedence prec
        cdef str precName
        cdef TokenSpec token, pToken
        cdef NontermSpec nonterm, pNonterm

        if self._verbose:
            print ("Parsing.Spec: Introspecting module%s to acquire formal" + \
            " grammar specification...") % ("s", "")[len(modules) == 1]

        deferTokens = []
        deferNonterms = []
        for module in modules:
            d = module.__dict__
            for k in d:
                v = d[k]
                if type(v) is types.TypeType and type(v.__doc__) is str:
                    dirtoks = v.__doc__.split("\n")[-1].split(" ")

                    #===========================================================
                    # Precedence.
                    #
                    if issubclass(v, Precedence) and dirtoks[0] in \
                      ["%fail", "%nonassoc", "%left", "%right", "%split"]:
                        name = k
                        relationships = {}
                        i = 1
                        while i < len(dirtoks):
                            tok = dirtoks[i]
                            m = re.compile(r'([<>=])([A-Za-z]\w*)').match(tok)
                            if m:
                                # Precedence relationship.
                                if m.group(2) in relationships:
                                    raise SpecError("Duplicate precedence " \
                                      "relationship: %s" \
                                      % v.__doc__)
                                relationships[m.group(2)] = m.group(1)
                            else:
                                m = re.compile(r'([A-Za-z]\w*)').match(tok)
                                if m:
                                    if i != 1:
                                        raise SpecError("Precedence name " \
                                          "must come before relationships: %s" \
                                          % v.__doc__)
                                    name = m.group(1)
                                else:
                                    raise SpecError( \
                                      "Invalid precedence specification: %s" % \
                                      v.__doc__)
                            i += 1

                        if name in self._precedences \
                          and (<Precedence>self._precedences[name]).assoc \
                          is not None:
                            raise SpecError("Duplicate precedence name: %s" % \
                              v.__doc__)
                        if name in self._tokens:
                            raise SpecError("Identical token/precedence " \
                              "names: %s" % v.__doc__)
                        if name in self._nonterms:
                            raise SpecError("Identical nonterm/precedence " \
                              "names: %s" % v.__doc__)
                        if name not in self._precedences:
                            prec = Precedence(name, dirtoks[0][1:],
                              relationships)
                            self._precedences[name] = prec
                        else:
                            prec = <Precedence>self._precedences[name]
                            prec.specify(dirtoks[0][1:], relationships)
                    #===========================================================
                    # Token.
                    #
                    elif issubclass(v, Token) and \
                      dirtoks[0] in ("%token", "%extend"):
                        name = None
                        precName = None
                        i = 1
                        while i < len(dirtoks):
                            tok = dirtoks[i]
                            m = re.compile(r'\[([A-Za-z]\w*)\]').match(tok)
                            if m:
                                if i < len(dirtoks) - 1:
                                    raise SpecError("Precedence must come " \
                                      "last in token specification: %s" % \
                                      v.__doc__)
                                precName = m.group(1)
                            else:
                                m = re.compile(r'([A-Za-z]\w*)').match(tok)
                                if m:
                                    name = m.group(1)
                                else:
                                    raise SpecError("Invalid token " \
                                      "specification: %s" % v.__doc__)
                            i += 1
                        if name is None:
                            name = k
                        if precName is None:
                            precName = "none"
                        if name in self._precedences:
                            raise SpecError("Identical precedence/token " \
                              "names: %s" % v.__doc__)
                        if name in self._nonterms:
                            raise SpecError("Identical nonterm/token names: " \
                              "%s" % v.__doc__)
                        if precName not in self._precedences:
                            prec = Precedence(precName)
                            self._precedences[precName] = prec
                        else:
                            prec = self._precedences[precName]
                        if dirtoks[0] == "%extend":
                            # Defer until all other introspection is complete.
                            token = TokenSpec(v, name,
                              "%s.%s" % (module.__name__, k), prec)
                            deferTokens.append(token)
                            continue
                        if name in self._tokens:
                            raise SpecError("Duplicate token name: %s" % \
                              v.__doc__)
                        token = TokenSpec(v, name,
                          "%s.%s" % (module.__name__, k), prec)
                        self._tokens[name] = token
                        self._sym2spec[v] = token
                    #===========================================================
                    # Nonterm.
                    #
                    elif issubclass(v, Nonterm) and \
                      dirtoks[0] in ("%start", "%nonterm", "%extend"):
                        name = None
                        precName = None
                        i = 1
                        while i < len(dirtoks):
                            tok = dirtoks[i]
                            m = re.compile(r'\[([A-Za-z]\w*)\]').match(tok)
                            if m:
                                if i < len(dirtoks) - 1:
                                    raise SpecError("Precedence must come " \
                                      "last in non-terminal specification: %s" \
                                      % v.__doc__)
                                precName = m.group(1)
                            else:
                                m = re.compile(r'([A-Za-z]\w*)').match(tok)
                                if m:
                                    name = m.group(1)
                                else:
                                    raise SpecError("Invalid non-terminal " \
                                      "specification: %s" % v.__doc__)
                            i += 1
                        if name is None:
                            name = k
                        if precName is None:
                            precName = "none"
                        if name in self._precedences:
                            raise SpecError("Identical precedence/nonterm " \
                              "names: %s" % v.__doc__)
                        if name in self._tokens:
                            raise SpecError( "Identical token/nonterm names: " \
                              "%s" % v.__doc__)
                        if precName not in self._precedences:
                            prec = Precedence(precName)
                            self._precedences[precName] = prec
                        else:
                            prec = self._precedences[precName]
                        if dirtoks[0] == "%extend":
                            # Defer until all other introspection is complete.
                            nonterm = NontermSpec(v, name,
                              "%s.%s" % (module.__name__, k), prec)
                            deferNonterms.append(nonterm)
                            continue
                        if name in self._nonterms:
                            raise SpecError("Duplicate nonterm name: %s" % \
                              v.__doc__)
                        nonterm = NontermSpec(v, name,
                          "%s.%s" % (module.__name__, k), prec)
                        self._nonterms[name] = nonterm
                        self._sym2spec[v] = nonterm

                        if dirtoks[0] == "%start":
                            # Start symbol.
                            if self._explicitStartSym is None:
                                if self._userStartSym is not None:
                                    raise SpecError("Only one start " \
                                      "non-terminal allowed: %s" % v.__doc__)
                                self._userStartSym = nonterm
                            elif self._explicitStartSym == v:
                                self._userStartSym = nonterm
                    #===========================================================

        # Link deferred %extend directives together with their parents.
        while len(deferTokens) > 0:
            i = 0
            while i < len(deferTokens):
                token = <TokenSpec>deferTokens[i]
                if token.name not in self._tokens:
                    raise SpecError("Missing parent for extended " \
                      "token: %s" % token.tokenType)
                pToken = <TokenSpec>self._tokens[token.name]
                if token.tokenType.__base__ is pToken.tokenType:
                    deferTokens.pop(i)
                    token.chain = pToken.chain
                    token.chain.append(token)
                    # Replace parent's entries in _tokens and _sym2spec.
                    self._tokens[token.name] = token
                    self._sym2spec.pop(pToken.tokenType)
                    self._sym2spec[token.tokenType] = token
                else:
                    if not issubclass(token.tokenType, pToken.tokenType):
                        raise SpecError("Extended token fork: %r: %r" % \
                          (token.tokenType, token.tokenType.__doc__))
                    i += 1

        while len(deferNonterms) > 0:
            i = 0
            while i < len(deferNonterms):
                nonterm = <NontermSpec>deferNonterms[i]
                if nonterm.name not in self._nonterms:
                    raise SpecError("Missing parent for extended " \
                      "non-terminal: %s" % nonterm.nontermType)
                pNonterm = <NontermSpec>self._nonterms[nonterm.name]
                if nonterm.nontermType.__base__ is pNonterm.nontermType:
                    deferNonterms.pop(i)
                    nonterm.chain = pNonterm.chain
                    nonterm.chain.append(nonterm)
                    # Replace parent's entries in _nonterms and _sym2spec.
                    self._nonterms[nonterm.name] = nonterm
                    self._sym2spec.pop(pNonterm.nontermType)
                    self._sym2spec[nonterm.nontermType] = nonterm
                else:
                    if not issubclass(nonterm.nontermType, \
                      pNonterm.nontermType):
                        raise SpecError("Extended non-nonterminal fork: %r: " \
                          "%r" % (nonterm.nontermType, \
                          nonterm.nontermType.__doc__))
                    i += 1

        if not isinstance(self._userStartSym, NontermSpec):
            raise SpecError("No start symbol specified")

    # Resolve all symbolic (named) references.
    cdef void _references(self, str logFile, str graphFile) except *:
        cdef TokenSpec token
        cdef NontermSpec nonterm, link
        cdef dict meth2prod
        cdef str k
        cdef object d # dictproxy.
        cdef object v
        cdef list dirtoks, rhs
        cdef str precName
        cdef Precedence prec
        cdef object m
        cdef Production prod
        cdef int i
        cdef str tok

        # Build the graph of Precedence relationships.
        self._resolvePrec(graphFile)

        # Verify that Token-->Precedence references are valid.
        for token in self._tokens.itervalues():
            if token.prec.assoc is None:
                raise SpecError("Unknown precedence in Token specification " \
                  " for %r: %s" % (token.tokenType, token.tokenType.__doc__))

        # Resolve Nonterm-->{Nonterm,Token,Precedence} references.
        for nonterm in self._nonterms.itervalues():
            # Verify that Nonterm-->Prececence references are valid.
            if nonterm.prec.assoc is None:
                raise SpecError("Unknown precedence in Nonterm " \
                  "specification for %r: %s" % (nonterm.nontermType, \
                  nonterm.nontermType.__doc__))

            # Update the start symbol if it is in this chain.
            if nonterm.chain[0] == self._userStartSym:
                self._userStartSym = nonterm

            meth2prod = {}

            # Iterate over the chain and merge productions into meth2prod.
            for link in nonterm.chain:
                # Resolve Nonterm-->Precedence references.
                if type(link.prec) is str:
                    link.prec = self._precedences[nonterm.prec]

                d = link.nontermType.__dict__
                for k in d:
                    v = d[k]
                    if inspect.isroutine(v) and type(v.__doc__) is str:
                        dirtoks = v.__doc__.split("\n")[-1].split(" ")
                        if dirtoks[0] == "%reduce":
                            if k in meth2prod:
                                raise SpecError("Production exists in " \
                                  "ancestor: %r.%r: %r" % \
                                  (link.nontermType.__name__, k, v.__doc__))
                        elif dirtoks[0] in ("%accept", "%amend", "%suppress"):
                            if k not in meth2prod:
                                raise SpecError("Production does not exist " \
                                  "in ancestor: %r.%r: %r" % \
                                  (link.nontermType.__name__, k, v.__doc__))

                        if dirtoks[0] in ("%reduce", "%amend"):
                            rhs = []
                            prec = None
                            for 1 <= i < len(dirtoks):
                                tok = dirtoks[i]
                                m = re.compile(r'([A-Za-z]\w*)').match(tok)
                                if m:
                                    # Symbolic reference.
                                    if tok in self._tokens:
                                        rhs.append(self._tokens[tok])
                                    elif tok in self._nonterms:
                                        rhs.append(self._nonterms[tok])
                                    else:
                                        raise SpecError("Unknown symbol '%s' " \
                                          "in reduction specification: %s" % \
                                          (tok, v.__doc__))
                                else:
                                    m = re.compile(r'\[([A-Za-z]\w*)\]'). \
                                      match(tok)
                                    if m:
                                        # Precedence.
                                        if i < len(dirtoks) - 1:
                                            raise SpecError("Precedence must " \
                                              "come last in reduction " \
                                              "specification: %s" % \
                                              v.__doc__)
                                        if (<Precedence>self._precedences[ \
                                          m.group(1)]).assoc is None:
                                            raise SpecError("Unknown " \
                                              "precedence in reduction " \
                                              "specification: %s" % v.__doc__)
                                        prec = self._precedences[m.group(1)]

                            if prec is None:
                                # Inherit the non-terminal's precedence.
                                prec = link.prec

                            prod = Production(v, "%s.%s" % \
                              (link.qualified, k), prec, nonterm, rhs)

                            meth2prod[k] = prod
                        elif dirtoks[0] == "%accept":
                            # Everything but the associated method stays the
                            # same.
                            (<Production>meth2prod[k]).method = v
                        elif dirtoks[0] == "%suppress":
                            # White out existing entry, rather than removing it,
                            # so that it is possible to properly process
                            # associated directives further down the chain.
                            meth2prod[k] = None

            # Now that the chain has been merged, store the resulting
            # productions.
            for prod in meth2prod.itervalues():
                if prod is not None:
                    assert prod not in link.productions
                    nonterm.productions.append(prod)
                    self._productions.append(prod)

        # Make sure all referenced Precedence names are defined.  All of the
        # above validation should have already found all references to
        # undefined Precedence names.
        if __debug__:
            for prec in self._precedences.itervalues():
                assert prec.assoc is not None

    # Build the graph of Precedence relationships.
    cdef void _resolvePrec(self, str graphFile) except *:
        cdef Precedence precA, precB, precC, precD, prec, p
        cdef str precBName, rel
        cdef file f
        cdef bint done
        cdef list cycles

        # Resolve symbolic references and populate equiv/dominators.
        for precA in self._precedences.itervalues():
            for precBName in precA.relationships:
                if precBName not in self._precedences:
                    raise SpecError("Precedence '%s' specifies a " \
                      "relationship with unknown Precedence '%s'" % \
                      (precA, precBName))
                precB = self._precedences[precBName]
                rel = precA.relationships[precBName]
                if rel == "=":
                    precA.equiv.append(precB)
                elif rel == "<":
                    if precB not in precA.dominators:
                        precA.dominators.append(precB)
                elif rel == ">":
                    if precA not in precB.dominators:
                        precB.dominators.append(precA)
                else:
                    assert False

        # Create equivalence classes for all Precedence classes.  Since the
        # Precedence classes are equivalent, they also share dominator sets.
        for precA in self._precedences.itervalues():
            for precB in precA.equiv[:]:
                if not precB.equiv is precA.equiv:
                    # Merge the sets of equivalent Precedence classes.
                    for prec in precB.equiv:
                        if prec not in precA.equiv:
                            precA.equiv.append(prec)
                    # Share the equiv set.
                    for prec in precA.equiv:
                        prec.equiv = precA.equiv

        # Use the equivalence classes to merge dominator sets and share them.
        for precA in self._precedences.itervalues():
            for precB in precA.equiv[1:]:
                # Merge the sets of dominator Precedence classes.
                for prec in precB.dominators:
                    if prec not in precA.dominators:
                        precA.dominators.append(prec)
                # Share the dominator set.
                precB.dominators = precA.dominators

        # Write graphviz precedence graph to graphFile, if graphFile was
        # specified.
        if graphFile is not None:
            f = open(graphFile, "w+")
            if self._verbose:
                print \
                  "Parsing.Spec: Writing graphviz precedence graph to '%s'..." \
                  % graphFile
            f.write('digraph Precedence {\n')
            f.write('    graph [bgcolor=black, labeljust="l"]\n')
            f.write(\
              ('    node [shape=record, style=filled, color=black, ' + \
              'fillcolor=gray, fontname=Helvetica, fontsize=10.0]\n'))
            f.write('    edge [color=gray]\n')
            for precA in self._precedences.itervalues():
                if precA == precA.equiv[0]:
                    f.write(\
                      ('    Precedence_%s [label="{%s}"]\n') % (precA.name, \
                      "\\n".join(["%s (%s)" % (p.name, p.assoc) \
                      for p in precA.equiv])))
                    for precB in precA.dominators:
                        f.write('    Precedence_%s -> Precedence_%s\n' % \
                          ((<Precedence>precB.equiv[0]).name, \
                          (<Precedence>precA.equiv[0]).name))
            f.write('}\n')
            f.close()

        # Iteratively build dominator sets until no more work can be done.
        done = False
        while not done:
            done = True
            for precA in self._precedences.itervalues():
                if precA == precA.equiv[0]: # No need to do more than this.
                    for precB in precA.dominators[:]:
                        for precC in precB.equiv:
                            if precC not in precA.dominators:
                                precA.dominators.append(precC)
                                done = False
                        for precC in precB.dominators:
                            for precD in precC.equiv:
                                if precD not in precA.dominators:
                                    precA.dominators.append(precD)
                                    done = False

        # Check for cycles in the graph.
        cycles = []
        for precA in self._precedences.itervalues():
            for precB in [precA] + precA.equiv:
                if precB in precA.dominators:
                    cycles.append( \
                      "Precedence relationship cycle involving '%s'" % \
                      precA.name)
        if len(cycles) > 0:
            raise SpecError("\n".join(cycles))

    # Store state to a pickle file, if requested.
    cdef void _pickle(self, object file, mode):
        if self._skinny:
            # Discard bulky data that don't need to be pickled.
            #self._startSym = ...
            #self._startProd = ...
            self._itemSets = []
            self._itemSetsHash = {}
            #self._startState = ...

        if file is not None and "w" in mode:
            if self._verbose:
                print "Parsing.Spec: Creating %s Spec pickle in %s..." % \
                  (("fat", "skinny")[self._skinny], file)
            f = open(file, "w")
            cPickle.dump(self, f, protocol=cPickle.HIGHEST_PROTOCOL)
            f.close()

    # Restore state from a pickle file, if a compatible one is provided.  This
    # method uses the same set of return values as does _compatible().
    cdef _unpickle(self, object file, mode):
        cdef Spec spec
        cdef NontermSpec nontermSpec, nontermSelf
        cdef Production prodSpec, prodSelf
        cdef TokenSpec tokenSpec, tokenSelf

        if file is not None and "r" in mode:
            if self._verbose:
                print \
                  "Parsing.Spec: Attempting to use pickle from file \"%s\"..." \
                  % file
            try:
                f = open(file, "r")
            except IOError:
                if self._verbose:
                    error = sys.exc_info()
                    print "Parsing.Spec: Pickle open failed: Exception %s: %s" \
                      % (error[0], error[1])
                return "incompatible"

            # Any exception at all in unpickling can be assumed to be due to
            # an incompatible pickle.
            try:
                spec = cPickle.load(f)
            except:
                if self._verbose:
                    error = sys.exc_info()
                    print "Parsing.Spec: Pickle load failed: Exception %s: %s" \
                      % (error[0], error[1])
                return "incompatible"

            compat = self._compatible(spec)
            if compat == "incompatible":
                if self._verbose:
                    print "Parsing.Spec: Pickle in \"%s\" is incompatible." % \
                      file
                return compat

            if self._verbose:
                print \
                  "Parsing.Spec: Using %s pickle in \"%s\" (%s)..." \
                  % (("fat", "skinny")[spec._skinny], file, compat)

            if compat in ["compatible", "repickle"]:
                # Copy spec's data structures.
                self._precedences = spec._precedences
                self._action = spec._action
                self._goto = spec._goto
                self._startState = spec._startState
                self._nActions = spec._nActions
                self._nConflicts = spec._nConflicts
                self._nImpure = spec._nImpure
            elif compat == "itemsets":
                # Precedences are incompatible, so great care has to be taken
                # when copying from the pickle.  Overwrite all precedence
                # specifications in spec with the new ones, then copy over all
                # of the new symbols/productions (but not the new precedences,
                # of course).  This still leaves table generation, which is
                # done by the _prepare() method later.

                # Nonterminals.
                for key in self._nonterms:
                    nontermSelf = self._nonterms[key]
                    nontermSpec = spec._nonterms[key]
                    nontermSpec.prec = nontermSelf.prec
                    # Productions.
                    for prodSelf in nontermSelf.productions:
                        for prodSpec in nontermSpec.productions:
                            if prodSelf.qualified == prodSpec.qualified:
                                prodSpec.prec = prodSelf.prec
                                break
                        assert prodSelf.qualified == prodSpec.qualified
                # Tokens.
                for key in self._tokens:
                    tokenSelf = self._tokens[key]
                    tokenSpec = spec._tokens[key]
                    tokenSpec.prec = tokenSelf.prec
            else:
                assert False

            # Copy spec data structures that are usable regardless of whether
            # the parsing tables need to be rebuilt.
            self._nonterms = spec._nonterms
            self._tokens = spec._tokens
            self._sym2spec = spec._sym2spec
            self._productions = spec._productions
            self._userStartSym = spec._userStartSym
            self._startSym = spec._startSym
            self._startProd = spec._startProd
            self._itemSets = spec._itemSets
            self._itemSetsHash = spec._itemSetsHash

            return compat
        else:
            return "incompatible"

    # Determine whether other is compatible with self.  Note that self is not
    # completely initialized; the idea here is to determine whether other's
    # data structures can be copied *before* doing the work of building parsing
    # tables.
    #
    # Itemsets and precedences are not directly related, other than that
    # symbols have precedences associated with them.  Therefore, we check for
    # the following cases:
    #
    #   "compatible" : Completely compatible.
    #
    #   "repickle" : Compatible, but pickle needs to be regenerated.
    #
    #   "itemsets" : Itemsets are compatible, but precedence specifications are
    #                not.
    #
    #   "incompatible" : No useful compatibility.
    cdef _compatible(self, Spec other):
        cdef object ret
        cdef object key
        cdef object prec
        cdef Precedence precA, precB
        cdef NontermSpec nontermA, nontermB
        cdef Production prodA, prodB
        cdef TokenSpec tokenA, tokenB

        ret = "compatible"

        if (not self._skinny) and other._skinny:
            return "incompatible"
        elif self._skinny != other._skinny:
            ret = "repickle"

        # Precedences.
        if len(self._precedences) != len(other._precedences):
            if self._verbose:
                print "Parsing.Spec: Unequal number of precedences (%d vs %d)" \
                  % (len(self._precedences), len(other._precedences))
            ret = "itemsets"
        for key in self._precedences:
            if key not in other._precedences:
                if self._verbose:
                    print "Parsing.Spec: Missing precedence: %s" % key
                ret = "itemsets"
                continue
            precA = self._precedences[key]
            precB = other._precedences[key]
            if precA.name != precB.name \
              or precA.assoc != precB.assoc \
              or len(precA.relationships) != len(precB.relationships):
                if self._verbose:
                    print "Parsing.Spec: Incompatible precedences: %r vs. %r" \
                      % (precA, precB)
                ret = "itemsets"
                continue
            for prec in precA.relationships:
                rel = precA.relationships[prec]
                if prec not in precB.relationships \
                  or precB.relationships[prec] != rel:
                    if self._verbose:
                        print \
                          "Parsing.Spec: Incompatible precedences: %r vs. %r" \
                          % (precA, precB)
                    ret = "itemsets"
                    break

        # Nonterminals.
        if len(self._nonterms) != len(other._nonterms):
            if self._verbose:
                print \
                  "Parsing.Spec: Unequal number of non-terminals (%d vs %d)" \
                  % (len(self._nonterms), len(other._nonterms))
            return "incompatible"
        for key in self._nonterms:
            if key not in other._nonterms:
                if self._verbose:
                    print "Parsing.Spec: Missing non-terminal: %s" % key
                return "incompatible"
            nontermA = self._nonterms[key]
            nontermB = other._nonterms[key]
            if nontermA.name != nontermB.name \
              or nontermA.qualified != nontermB.qualified \
              or nontermA.nontermType != nontermB.nontermType:
                if self._verbose:
                    print \
                      "Parsing.Spec: Incompatible non-terminals: %r vs. %r" \
                      % (nontermA, nontermB)
                return "incompatible"
            if nontermA.prec.name != nontermB.prec.name:
                if self._verbose:
                    print \
                      ("Parsing.Spec: Differing precedences for " + \
                      "non-terminal: %r") % nontermA
                ret = "itemsets"

            # Productions.
            if len(nontermA.productions) != len(nontermB.productions):
                if self._verbose:
                    print \
                      "Parsing.Spec: Unequal number of productions (%d vs %d)" \
                      % (len(self._productions) - 1, \
                      len(other._productions) - 1)
                return "incompatible"
            for prodA in nontermA.productions:
                match = False
                for prodB in nontermB.productions:
                    if prodA.qualified == prodB.qualified \
                      and prodA.lhs.name == prodB.lhs.name \
                      and len(prodA.rhs) == len(prodB.rhs):
                        match = True
                        for 0 <= i < len(prodA.rhs):
                            if (<SymbolSpec>prodA.rhs[i]).name != \
                              (<SymbolSpec>prodB.rhs[i]).name:
                                match = False
                                if self._verbose:
                                    print \
                                      ("Parsing.Spec: Incompatible" + \
                                      " productions: %r vs. %r") \
                                      % (prodA, prodB)
                                break
                        if prodA.prec.name != prodB.prec.name:
                            if self._verbose:
                                print \
                                  ("Parsing.Spec: Differing precedences " + \
                                  "for production: %r") % prodA
                            ret = "itemsets"
                if not match:
                    return "incompatible"

        # Tokens.
        if len(self._tokens) != len(other._tokens):
            if self._verbose:
                print "Parsing.Spec: Unequal number of tokens (%d vs %d)" \
                  % (len(self._tokens), len(other._tokens))
            return "incompatible"
        for key in self._tokens:
            if key not in other._tokens:
                if self._verbose:
                    print "Parsing.Spec: Missing token: %s" % key
                return "incompatible"
            tokenA = self._tokens[key]
            tokenB = other._tokens[key]
            if tokenA.name != tokenB.name \
              or tokenA.tokenType != tokenB.tokenType:
                if self._verbose:
                    print \
                      "Parsing.Spec: Incompatible tokens: %r vs. %r" \
                      % (tokenA, tokenB)
                return "incompatible"
            if tokenA.prec.name != tokenB.prec.name:
                if self._verbose:
                    print \
                      "Parsing.Spec: Differing precedences for token: %r" \
                      % tokenA
                ret = "itemsets"

        # User start symbol.
        if self._userStartSym.name != other._userStartSym.name:
            if self._verbose:
                print "Parsing.Spec: Differing start symbols: %s vs. %s" \
                  % (self._userStartSym.name, other._userStartSym.name)
            return "incompatible"

        if other._skinny and ret == "itemsets":
            # The itemsets have to be regenerated, since they weren't pickled.
            ret = "incompatible"
        return ret

    # Check for unused prececence/token/nonterm/reduce specifications, then
    # throw a SpecError if any ambiguities exist in the grammar.
    cdef void _validate(self, str logFile) except *:
        cdef list lines, productions
        cdef dict used
        cdef ItemSet itemSet
        cdef Item item
        cdef int nUnused
        cdef SymbolSpec sym
        cdef TokenSpec token
        cdef str name
        cdef Production production
        cdef int ntokens, nnonterms, nproductions

        if self._verbose:
            print "Parsing.Spec: Validating grammar..."

        lines = []
        if self._nConflicts > 0:
            lines.append("Parsing.Spec: %d unresolvable conflict%s" % \
              (self._nConflicts, ("s", "")[self._nConflicts == 1]))

        # Previous code guarantees that all precedence/token/nonterm names are
        # unique.  Therefore, we can build a single dictionary here that keys on
        # names.
        used = {}
        productions = []
        for itemSet in self._itemSets:
            for item in itemSet:
                productions.append(item.production)
                used[item.production.prec.name] = item.production.prec
                for sym in [item.production.lhs] + item.production.rhs:
                    used[sym.name] = sym
                    used[sym.prec.name] = sym.prec

                for token in item.lookahead.iterkeys():
                    used[token.prec.name] = token.prec

        nUnused = 0

        # Precedences.
        for name in self._precedences:
            if name not in (self._none.name, self._split.name):
                if name not in used:
                    nUnused += 1
                    lines.append("Parsing.Spec: Unused precedence: %r" % \
                      self._precedences[name])

        # Tokens.
        for name in self._tokens:
            if name not in (self._eoi.name, self._epsilon.name):
                if name not in used:
                    nUnused += 1
                    lines.append("Parsing.Spec: Unused token: %s" % \
                      self._tokens[name])

        # Nonterms.
        for name in self._nonterms:
            if name not in (self._startSym.name,):
                if name not in used:
                    nUnused += 1
                    lines.append("Parsing.Spec: Unused nonterm: %s" % \
                      self._nonterms[name])

        # Productions.
        for production in self._productions:
            if production not in productions:
                nUnused += 1
                lines.append("Parsing.Spec: Unused production: %r" % production)

        if nUnused > 0:
            lines.insert((1, 0)[self._nConflicts == 0], \
              "Parsing.Spec: %d unused definition%s" % \
              (nUnused, ("s", "")[nUnused == 1]))

        # Write to logFile, if one was specified.
        if logFile is not None:
            f = open(logFile, "w+")
            if self._verbose:
                print "Parsing.Spec: Writing log to '%s'..." % logFile
            f.write("%s" % "\n".join(lines + ["%r" % self]))
            f.close()

        # Conflicts are fatal.
        if self._nConflicts > 0:
            raise SpecError("%s" % ("\n".join(lines)))

        # Make sure to let the user know about unused symbols if verbosity is
        # enabled, and there weren't any conflicts to cause notification via an
        # exception.
        if self._verbose:
            ntokens = len(self._tokens) - 2
            nnonterms = len(self._nonterms) - 1
            nproductions = len(self._productions) - 1
            lines.append(
              "Parsing.Spec: %d token%s, %d non-terminal%s, %d production%s" \
              % (ntokens, ("s", "")[ntokens == 1], \
              nnonterms, ("s", "")[nnonterms == 1], \
              nproductions, ("s", "")[nproductions == 1]))
            sys.stdout.write("%s\n" % "\n".join(lines))

    # Compute the first sets for all symbols.
    cdef void _firstSets(self):
        cdef SymbolSpec elm, elmSym
        cdef TokenSpec tokenSpec
        cdef NontermSpec nontermSpec
        cdef bint done, containsEpsilon
        cdef Production prod

        # Terminals.
        # first(X) is X for terminals.
        for tokenSpec in self._tokens.itervalues():
            tokenSpec.firstSetMerge(tokenSpec)

        # Non-terminals.
        #
        # Repeat the following loop until no more symbols can be added to any
        # first set.
        done = False
        while not done:
            done = True
            for name in self._nonterms:
                nontermSpec = self._nonterms[name]
                for prod in nontermSpec.productions:
                    # Merge epsilon if there is an empty production.
                    if len(prod.rhs) == 0:
                        if not nontermSpec.firstSetMerge(self._epsilon):
                            done = False

                    # Iterate through the RHS and merge the first sets into
                    # this symbol's, until a preceding symbol's first set does
                    # not contain epsilon.
                    for elm in prod.rhs:
                        containsEpsilon = False
                        for elmSym in elm.firstSet:
                            if not nontermSpec.firstSetMerge(elmSym):
                                done = False
                            if elmSym is self._epsilon:
                                containsEpsilon = True
                        if not containsEpsilon:
                            break

    # Compute the follow sets for all symbols.
    cdef void _followSets(self):
        cdef int i, j
        cdef bint done
        cdef object name
        cdef SymbolSpec rhsSym
        cdef NontermSpec nontermSpec
        cdef Production prod

        self._startSym.followSet = [self._epsilon]

        # Repeat the following loop until no more symbols can be added to any
        # follow set.
        done = False
        while not done:
            done = True
            for name in self._nonterms:
                nontermSpec = self._nonterms[name]
                for prod in nontermSpec.productions:
                    # For all A ::= aBb, merge first(b) into follow(B).
                    for 0 <= i < len(prod.rhs) - 1:
                        for i+1 <= j < len(prod.rhs):
                            rhsSym = prod.rhs[i]
                            if not rhsSym.followSetMerge( \
                              (<SymbolSpec>prod.rhs[j]).firstSet, \
                              self._epsilon):
                                done = False
                            if self._epsilon not in \
                              (<SymbolSpec>prod.rhs[j]).firstSet:
                                break

                    # For A ::= ab, or A ::= aBb where first(b) contains <e>,
                    # merge follow(A) into follow(B).
                    for len(prod.rhs)-1 >= i > -1:
                        rhsSym = prod.rhs[i]
                        if not rhsSym.followSetMerge(prod.lhs.followSet, \
                          self._epsilon):
                            done = False
                        if self._epsilon not in \
                          (<SymbolSpec>prod.rhs[i]).firstSet:
                            break

    # Compute the collection of sets of LR(1) items.
    cdef void _items(self):
        cdef ItemSet tItemSet, itemSet, gotoSet, mergeSet
        cdef Item tItem
        cdef list worklist, syms
        cdef dict itemSetsHash
        cdef int nwork, i, j, k
        cdef bint merged
        cdef SymbolSpec sym
        cdef StringFirstSetCache cache

        cache = StringFirstSetCache(self._epsilon)

        # Add {[S' ::= * S $., <e>]} to _itemSets.
        tItemSet = ItemSet()
        tItem = Item(self._startProd, 0, {self._epsilon: self._epsilon})
        tItemSet.append(tItem)
        tItemSet.closure(cache)
        self._itemSets.append(tItemSet)

        # List of state numbers that need to be processed.
        worklist = [0]
        if self._verbose:
            nwork = len(worklist)
            print "Parsing.Spec: Generating LR(1) itemset collection... ",
            sys.stdout.write("+")
            sys.stdout.flush()
        else:
            nwork = 0 # Silence compiler warning.

        # itemSetsHash uses itemsets as keys.  A value is a list of _itemSets
        # indices; these itemsets are the ones referred to the key itemset.
        itemSetsHash = {tItemSet: [0]}

        syms = self._tokens.values() + self._nonterms.values()
        while len(worklist) > 0:
            if self._verbose:
                if abs(len(worklist) - nwork) >= 10:
                    nwork = len(worklist)
                    sys.stdout.write("[%d/%d]" % \
                      (len(worklist), len(self._itemSets)))
                    sys.stdout.flush()

            i = worklist.pop(0)
            itemSet = self._itemSets[i]
            for 0 <= j < len(syms):
                sym = syms[j]
                gotoSet = itemSet.xgoto(sym)
                if len(gotoSet) > 0:
                    merged = False
                    if gotoSet in itemSetsHash:
                        for k in itemSetsHash[gotoSet]:
                            mergeSet = self._itemSets[k]
                            if mergeSet.weakCompat(gotoSet):
                                merged = True
                                if mergeSet.merge(gotoSet, cache):
                                    if k not in worklist:
                                        worklist.insert(0, k)
                                        if self._verbose:
                                            sys.stdout.write(".")
                                            sys.stdout.flush()
                                break
                    if not merged:
                        gotoSet.closure(cache)
                        worklist.append(len(self._itemSets))
                        if gotoSet not in itemSetsHash:
                            itemSetsHash[gotoSet] = [len(self._itemSets)]
                        else:
                            itemSetsHash[gotoSet].append(len(self._itemSets))
                        self._itemSets.append(gotoSet)
                        if self._verbose:
                            sys.stdout.write("+")
                            sys.stdout.flush()

        if self._verbose:
            sys.stdout.write("\n")
            sys.stdout.flush()
        self._itemSetsHash = itemSetsHash

    # Compute LR parsing tables.
    cdef void _lr(self):
        cdef dict itemSetsHash
        cdef ItemSet itemSet, itemSetB, itemSetC
        cdef dict state
        cdef Item item
        cdef SymbolSpec sym
        cdef int i

        # The collection of sets of LR(1) items already exists.
        assert len(self._itemSets) > 0
        assert len(self._action) == 0
        assert len(self._goto) == 0
        assert self._startState == -1
        assert self._nConflicts == 0

        if self._verbose:
            print \
              "Parsing.Spec: Generating LR(1) parsing tables (%d state%s)... " \
              % (len(self._itemSets), \
              ("s", "")[len(self._itemSets) == 1]),
            sys.stdout.flush()

        itemSetsHash = self._itemSetsHash

        for itemSet in self._itemSets:
            if self._verbose:
                sys.stdout.write(".")
                sys.stdout.flush()
            #===================================================================
            # _action.
            state = {}
            self._action.append(state)
            for item in itemSet:
                # X ::= a*Ab
                if item.dotPos < len(item.production.rhs):
                    sym = item.production.rhs[item.dotPos]
                    if isinstance(sym, TokenSpec):
                        itemSetB = itemSet.xgoto(sym)
                        for i in itemSetsHash[itemSetB]:
                            itemSetC = self._itemSets[i]
                            if itemSetC.weakCompat(itemSetB):
                                self._actionAppend(state, sym, ShiftAction(i))

                    # Check if this is the start state.
                    if self._startState == -1 \
                      and item.production.lhs == self._startSym \
                      and item.dotPos == 0:
                        assert len(item.production.rhs) == 2
                        self._startState = len(self._action) - 1
                # X ::= a*
                elif item.dotPos == len(item.production.rhs):
                    for lookahead in item.lookahead.iterkeys():
                        self._actionAppend(state, lookahead, \
                          ReduceAction(item.production))
                else:
                    assert False
            #===================================================================
            # _goto.
            state = {}
            self._goto.append(state)
            for nonterm in self._nonterms.itervalues():
                itemSetB = itemSet.xgoto(nonterm)
                if itemSetB in itemSetsHash:
                    for i in itemSetsHash[itemSetB]:
                        itemSetC = self._itemSets[i]
                        if itemSetC.weakCompat(itemSetB):
                            assert nonterm not in state
                            state[nonterm] = i

        if self._verbose:
            sys.stdout.write("\n")
            sys.stdout.flush()

    # Add a symbol action to state, if the action doesn't already exist.
    cdef void _actionAppend(self, state, sym, action):
        assert type(state) == dict
        assert isinstance(sym, SymbolSpec)
        assert isinstance(action, Action)

        if sym not in state:
            state[sym] = [action]
        else:
            actions = state[sym]
            if action not in actions:
                state[sym].append(action)

    # Look for action ambiguities and resolve them if possible.
    cdef void _disambiguate(self):
        cdef unsigned stateInd, vNConflicts, nConflicts, nActions, nImpure
        cdef dict state
        cdef str vRes, res
        cdef list actStats, acts, newActs
        cdef Action act, actI, actJ

        assert self._nActions == 0
        assert self._nConflicts == 0
        assert self._nImpure == 0

        if self._verbose:
            print "Parsing.Spec: Disambiguating LR(1) parsing tables... ",
            sys.stdout.flush()

        for 0 <= stateInd < len(self._action):
            state = self._action[stateInd]
            if self._verbose:
                vRes = "."
                vNConflicts = 0
            for sym in state:
                nConflicts = 0
                acts = [act for act in state[sym]]
                # Construct a list that corresponds to acts; each element
                # indicates whether to preserve the action.
                actStats = [True] * len(acts)

                # Fill in the cells of actStats.
                for i in xrange(len(acts)):
                    actI = acts[i]
                    for j in xrange(i+1, len(acts)):
                        actJ = acts[j]
                        res = self._resolve(sym, actI, actJ)
                        if res == "neither":
                            actStats[i] = False
                            actStats[j] = False
                        elif res == "old":
                            actStats[j] = False
                        elif res == "both":
                            pass
                        elif res == "new":
                            actStats[i] = False
                        elif res == "err":
                            actStats[i] = False
                            actStats[j] = False
                            nConflicts += 1
                        else:
                            assert False

                # Look for actions that can coexist or dominate all other
                # actions.
                newActs = []
                for j in xrange(len(acts)):
                    if actStats[j]:
                        newActs.append(acts[j])
                # Replace the action set if there exists a valid resolution
                # among the actions.
                if len(newActs) > 0 or nConflicts == 0:
                    if self._verbose:
                        if len(newActs) != len(acts):
                            vRes = "_"
                    state[sym] = newActs
                    nConflicts = 0
                elif self._verbose:
                    vNConflicts += nConflicts

                nActions = len(state[sym])
                if nActions > 1:
                    nImpure = nActions
                else:
                    nImpure = 0

                # Update summary stats.
                self._nActions += nActions
                self._nConflicts += nConflicts
                self._nImpure += nImpure

            if self._verbose:
                if vNConflicts == 0:
                    sys.stdout.write("%s" % vRes)
                else:
                    sys.stdout.write("[%d:%d]" % (stateInd, vNConflicts))
                sys.stdout.flush()

        if self._verbose:
            sys.stdout.write("\n")
            sys.stdout.flush()

    # Compute how to resolve an action conflict.
    #
    # ret: "neither" : Discard both.
    #      "old"     : Keep old.
    #      "both"    : Keep both.
    #      "new"     : Keep new.
    #      "err"     : Unresolvable conflict.
    cdef str _resolve(self, SymbolSpec sym, Action oldAct, Action newAct):
        cdef str ret, assoc
        cdef Precedence oldPrec, newPrec

        if type(oldAct) == ShiftAction:
            oldPrec = sym.prec
        else:
            assert type(oldAct) == ReduceAction
            oldPrec = (<ReduceAction>oldAct).production.prec

        if type(newAct) == ShiftAction:
            newPrec = sym.prec
        else:
            assert type(newAct) == ReduceAction
            newPrec = (<ReduceAction>newAct).production.prec

        if oldPrec in newPrec.dominators:
            # Discard new action.
            ret = "old"
        elif newPrec in oldPrec.dominators:
            # Discard old action.
            ret = "new"
        elif oldPrec in newPrec.equiv:
            assert newPrec in oldPrec.equiv

            if oldPrec.assoc == "split" or newPrec.assoc == "split":
                ret = "both"
            elif type(newAct) == type(oldAct):
                assert type(newAct) == ReduceAction
                assert type(oldAct) == ReduceAction
                # Fatal reduce/reduce conflict.
                ret = "err"
            else:
                if oldPrec.assoc != "fail" and newPrec.assoc != "fail" \
                  and oldPrec.assoc != newPrec.assoc:
                    # Conflicting associativity.
                    ret = "err"
                else:
                    # Determine associativity.  If only one of the actions has
                    # %fail associativity, it is overridden by the other.
                    if oldPrec.assoc == "fail":
                        assoc = newPrec.assoc
                    else:
                        assoc = oldPrec.assoc
                    assert assoc in ["fail", "nonassoc", "left", "right"]

                    if assoc == "fail":
                        ret = "err"
                    elif assoc == "left":
                        if type(oldAct) == ShiftAction:
                            ret = "new"
                        else:
                            assert type(newAct) == ShiftAction
                            ret = "old"
                    elif assoc == "right":
                        if type(oldAct) == ShiftAction:
                            ret = "old"
                        else:
                            assert type(newAct) == ShiftAction
                            ret = "new"
                    elif assoc == "nonassoc":
                        ret = "neither"
                    else:
                        assert False
        else:
            if newPrec in oldPrec.equiv:
                print "%r <--> %r" % (oldPrec, newPrec)
            assert newPrec not in oldPrec.equiv
            # No specified relationship between precedences.
            ret = "err"

        return ret

cdef class Lr:
    """
        LR(1) parser.  The Lr class uses a Spec instance in order to parse
        input that is fed to it via the token() method, and terminated via the
        eoi() method.
    """
    def __init__(self, Spec spec):
        if __debug__:
            if type(self) == Lr:
                assert spec.pureLR
        assert spec.conflicts == 0
        self.spec = spec
        self.reset()
        self.verbose = False

    cpdef reset(self):
        """
            Reset the parser in preparation for parsing new input.
        """
        self.start = None
        self._stack = [(Epsilon(self), 0)]

    cpdef token(self, Token token):
        """
            Feed a token to the parser.
        """
        cdef TokenSpec tokenSpec

        tokenSpec = self.spec._sym2spec[type(token)]
        self._act(token, tokenSpec)

    cpdef eoi(self):
        """
            Signal end-of-input to the parser.
        """
        cdef Token token

        token = EndOfInput(self)
        self.token(token)

        assert (<Symbol>self._stack[-1][0]) == token # <$>.
        if self.verbose:
            self._printStack()
            print "   --> accept"
        self._stack.pop()

        self.start = [<Symbol>self._stack[1][0]]
        assert self.start[0].symSpec == self.spec._userStartSym

    cdef void _act(self, Symbol sym, SymbolSpec symSpec) except *:
        cdef int topState
        cdef list actions
        cdef Action action

        if self.verbose:
            self._printStack()
            print "INPUT: %r" % sym

        while True:
            topState = <int>self._stack[-1][1]
            if symSpec not in self.spec._action[topState]:
                raise SyntaxError("Unexpected token: %r" % sym)

            actions = self.spec._action[topState][symSpec]
            assert len(actions) == 1
            action = <Action>actions[0]

            if self.verbose:
                print "   --> %r" % action
            if type(action) == ShiftAction:
                self._stack.append((sym, (<ShiftAction>action).nextState))
                break
            else:
                assert type(action) == ReduceAction
                self._reduce((<ReduceAction>action).production)

            if self.verbose:
                self._printStack()

    cdef void _printStack(self):
        cdef object node
        cdef Symbol sym
        cdef int state

        print "STACK:",
        for node in self._stack:
            sym = <Symbol>node[0]
            print "%r" % sym,
        print
        print "      ",
        for node in self._stack:
            sym = <Symbol>node[0]
            state = <int>node[1]
            print "%r%s" % (state,
              (" " * (len("%r" % sym) - len("%r" % state)))),
        print

    cdef void _reduce(self, Production production) except *:
        cdef int nRhs, i, sLen, topState
        cdef list rhs
        cdef Symbol r

        nRhs = len(production.rhs)
        rhs = []
        sLen = len(self._stack)
        for sLen - nRhs <= i < sLen:
            rhs.append(self._stack[i][0])

        for 0 <= i < nRhs:
            self._stack.pop()

        r = self._production(production, rhs)

        topState = <int>self._stack[-1][1]
        self._stack.append((r, self.spec._goto[topState][production.lhs]))

    cdef Symbol _production(self, Production production, list rhs):
        cdef Symbol sym

        sym = production.lhs.nontermType(self)
        production.method(sym, *rhs)

        return sym

#===============================================================================
# Begin graph-structured stack (GSS) classes.
#

cdef class Gsse:
    """
        Graph-structured stack edge.
    """
    def __init__(self, Gssn below, Gssn above, Symbol value):
        self.node = below
        above._edges.append(self)
        self.value = value

    def __repr__(self):
        return "{%r}" % self.value

    def __richcmp__(Gsse self, Gsse other, int op):
        cdef bint cmp

        cmp = (self.node == other.node and self.value == other.value)
        if op == 2: # ==
            return cmp
        elif op == 3: # !=
            return not cmp
        else:
            return NotImplemented

cdef class _GssnEdgesIterHelper:
    cdef list _edges
    cdef int _index

    def __init__(self, Gssn gssn):
        self._edges = gssn._edges[:]
        self._index = 0

    def __iter__(self):
        return self

    def __next__(self):
        cdef Gsse ret

        if self._index >= len(self._edges):
            raise StopIteration
        ret = self._edges[self._index]
        self._index += 1
        return ret

cdef class _GssnNodesIterHelper:
    cdef list _nodes
    cdef int _index

    def __init__(self, Gssn gssn):
        self._nodes = [edge.node for edge in gssn._edges]
        self._index = 0

    def __iter__(self):
        return self

    def __next__(self):
        cdef Gssn ret

        if self._index >= len(self._nodes):
            raise StopIteration
        ret = self._nodes[self._index]
        self._index += 1
        return ret

cdef class _GssnPathsIterHelper:
    cdef list _paths
    cdef int _index

    def __init__(self, Gssn gssn, int pathLen):
        self._paths = []
        self._pathsRecurse(gssn, pathLen, [])
        self._index = 0

    def __iter__(self):
        return self

    def __next__(self):
        cdef list ret

        if self._index >= len(self._paths):
            raise StopIteration
        ret = self._paths[self._index]
        self._index += 1
        return ret

    cdef _pathsRecurse(self, Gssn gssn, int pathLen, list path):
        cdef Gsse edge

        path.insert(0, gssn)
        if pathLen == -1 and len(gssn._edges) == 0:
            self._paths.append(path[:])
        elif pathLen != -1 and len(path) - 1 == pathLen * 2:
            self._paths.append(path[:])
        else:
            for edge in gssn.edges():
                # Avoid infinite recursion due to <e>-production cycles.
                if len(path) < 3 or edge != path[1]:
                    path.insert(0, edge)
                    self._pathsRecurse(edge.node, pathLen, path)
                    path.pop(0)
        path.pop(0)

cdef class Gssn:
    """
        Graph-structured stack node.
    """
    def __init__(self, Gssn below, Symbol value, int nextState):
        assert isinstance(below, Gssn) or below is None

        self._edges = []
        if below is not None:
            Gsse(below, self, value)
        self.nextState = nextState

    def __repr__(self):
        return "[%d]" % self.nextState

    property edge:
        def __get__(self):
            assert len(self._edges) == 1
            return self._edges[0]

    def edges(self):
        return _GssnEdgesIterHelper(self)

    def nodes(self):
        return _GssnNodesIterHelper(self)

    # Iterate over all paths of length pathLen.  Path length is measured as the
    # number of edges in the path, so a path of length 0 still consists of a
    # single node.
    #
    # Each path is encoded as a list that alternates between nodes and edges,
    # where the first and last elements are always nodes.
    #
    # <e>-grammars can cause cycles, which requires that we avoid infinite
    # recursion.
    def paths(self, int pathLen=-1):
        return _GssnPathsIterHelper(self, pathLen)

#
# End graph-structured stack (GSS) classes.
#===============================================================================

cdef class Glr(Lr):
    """
        GLR parser.  The Glr class uses a Spec instance in order to parse input
        that is fed to it via the token() method, and terminated via the eoi()
        method.
    """
    def __init__(self, spec):
        Lr.__init__(self, spec)

    cpdef reset(self):
        """
            Reset the parser in preparation for parsing new input.
        """
        cdef Gssn top

        self.start = None

        # Initialize with a stack that is in the start state.
        self._gss = []
        top = Gssn(None, None, 0)
        self._gss.append(top)

        self._paths = []

    cpdef token(self, Token token):
        """
            Feed a token to the parser.
        """
        cdef TokenSpec tokenSpec

        if self.verbose:
            print "%s" % ("-" * 80)
            print "INPUT: %r" % token
        tokenSpec = self.spec._sym2spec[type(token)]
        self._act(token, tokenSpec)
        if len(self._gss) == 0:
            raise SyntaxError("Unexpected token: %r" % token)

    cpdef eoi(self):
        """
            Signal end-of-input to the parser.
        """
        cdef Token token
        cdef Gssn top
        cdef list path
        cdef Gsse edge

        token = EndOfInput(self)
        self.token(token)

        # Gather the start symbols from the stacks.
        self.start = []
        for top in self._gss:
            for path in top.paths():
                assert len(path) == 5
                if self.verbose:
                    print "   --> accept %r" % path
                edge = path[1]
                assert isinstance(edge.value, Nonterm)
                assert edge.value.symSpec == self.spec._userStartSym
                self.start.append(edge.value)

        if len(self.start) == 0:
            raise SyntaxError("Unexpected end of input")

        if self.verbose:
            print "Start: %r" % self.start
            print "%s" % ("-" * 80)

    cdef void _act(self, Symbol sym, SymbolSpec symSpec) except *:
        self._reductions(sym, symSpec)
        self._shifts(sym, symSpec)

    cdef void _reductions(self, Symbol sym, SymbolSpec symSpec) except *:
        cdef dict epsilons
        cdef int nReduces, i
        cdef Gssn top
        cdef list workQ, path
        cdef Action action
        cdef ReduceAction rAction
        cdef Production production

        # epsilons is a dictionary that maps production-->[tops].  The purpose
        # is to avoid repeating the same epsilon production on a particular
        # stack top.  Ordinary productions do not require this care because we
        # can notice when a path has already been used for a production.
        epsilons = {}

        if self.verbose:
            nReduces = 0

        # Enqueue work.
        workQ = []
        i = 0
        while i < len(self._gss):
            top = self._gss[i]
            if symSpec not in self.spec._action[top.nextState]:
                # Unexpected token for this stack.
                self._gss.pop(i)
            else:
                for action in self.spec._action[top.nextState][symSpec]:
                    if type(action) == ReduceAction:
                        rAction = <ReduceAction>action
                        if len(rAction.production.rhs) == 0:
                            if rAction.production not in epsilons:
                                assert len([path for path in top.paths(0)]) == 1
                                path = [p for p in top.paths(0)][0]
                                epsilons[rAction.production] = [top]
                                workQ.append((path, rAction.production))
                                if self.verbose:
                                    print "   --> enqueue(a) %r" % \
                                      rAction.production
                                    print "                  %r" % path
                            elif top not in epsilons[rAction.production]:
                                assert len([path for path in top.paths(0)]) == 1
                                path = [p for p in top.paths(0)][0]
                                epsilons[rAction.production].append(top)
                                workQ.append((path, rAction.production))
                                if self.verbose:
                                    print "   --> enqueue(b) %r" % \
                                      rAction.production
                                    print "                  %r" % path
                        else:
                            # Iterate over all reduction paths through stack and
                            # enqueue them.
                            for path in top.paths(len(rAction.production.rhs)):
                                workQ.append((path, rAction.production))
                                if self.verbose:
                                    print "   --> enqueue(c) %r" % \
                                      rAction.production
                                    print "                  %r" % path
                i += 1

        # Process the work queue.
        while len(workQ) > 0:
            (path, production) = workQ.pop(0)

            if self.verbose:
                print "   --> reduce %r" % production
                print "              %r" % path
                nReduces += 1

            self._reduceOne(workQ, epsilons, path, production, symSpec)

        if self.verbose:
            if nReduces > 0:
                self._printStack()

    cdef void _reduceOne(self, list workQ, dict epsilons, object path, \
      Production production, SymbolSpec symSpec) except *:
        cdef list rhs
        cdef Symbol r, value
        cdef Gssn below, top
        cdef bint done
        cdef Gsse edge

        assert len(path[1::2]) == len(production.rhs)

        # Build the list of RHS semantic values to pass to the reduction action.
        rhs = [edge.value for edge in path[1::2]]

        # Call the user reduction method.
        r = self._production(production, rhs)

        below = path[0]
        done = False
        for top in self._gss:
            if top.nextState == \
              self.spec._goto[below.nextState][production.lhs]:
                # top is compatible with the reduction result we want to add to
                # the set of stack tops.
                for edge in top.edges():
                    if edge.node == below:
                        # There is already a below<--top link, so merge
                        # competing interpretations.
                        if self.verbose:
                            print "   --> merge %r <--> %r" % (edge.value, r)
                        value = production.lhs.nontermType.merge(edge.value, r)
                        if self.verbose:
                            if value == edge.value:
                                print "             %s" % \
                                  ("-" * len("%r" % edge.value))
                            else:
                                print "             %s      %s" % \
                                  ((" " * len("%r" % edge.value)), \
                                  "-" * len("%r" % r))
                        edge.value = value
                        done = True
                        break
                if not done:
                    # Create a new below<--top link.
                    edge = Gsse(below, top, r)
                    if self.verbose:
                        print "   --> shift(b) %r" % top

                    # Enqueue reduction paths that were created as a result of
                    # the new link.
                    self._enqueueLimitedReductions(workQ, epsilons, edge, \
                      symSpec)
                    done = True
                break
        if not done:
            # There is no compatible stack top, so create a new one.
            top = Gssn(below, r, \
              self.spec._goto[below.nextState][production.lhs])
            self._gss.append(top)
            if self.verbose:
                print "   --> shift(c) %r" % \
                  self.spec._goto[below.nextState][production.lhs]
            self._enqueueLimitedReductions(workQ, epsilons, top.edge, symSpec)

    # Enqueue paths that incorporate edge.
    cdef void _enqueueLimitedReductions(self, list workQ, dict epsilons, \
      Gsse edge, SymbolSpec symSpec) except *:
        cdef Gssn top
        cdef Action action
        cdef ReduceAction rAction
        cdef list path

        for top in self._gss:
            if symSpec in self.spec._action[top.nextState]:
                for action in self.spec._action[top.nextState][symSpec]:
                    if type(action) == ReduceAction:
                        rAction = <ReduceAction>action
                        if len(rAction.production.rhs) == 0:
                            if self.spec._goto[top.nextState] \
                              [rAction.production.lhs] == top.nextState:
                                # Do nothing, since enqueueing a reduction
                                # would result in performing the same reduction
                                # twice.
                                pass
                            elif rAction.production not in epsilons:
                                path = [top]
                                epsilons[rAction.production] = [top]
                                workQ.append((path, rAction.production))
                                if self.verbose:
                                    print "   --> enqueue(d) %r" % \
                                      rAction.production
                                    print "                  %r" % path
                            elif top not in epsilons[rAction.production]:
                                path = [top]
                                epsilons[rAction.production].append(top)
                                workQ.append((path, rAction.production))
                                if self.verbose:
                                    print "   --> enqueue(e) %r" % \
                                      rAction.production
                                    print "                  %r" % path
                        else:
                            # Iterate over all reduction paths through stack and
                            # enqueue them if they incorporate edge.
                            for path in top.paths(len(rAction.production.rhs)):
                                if edge in path[1::2]:
                                    workQ.append((path, rAction.production))
                                    if self.verbose:
                                        print "   --> enqueue(f) %r" % \
                                          rAction.production
                                        print "                  %r" % path

    cdef void _shifts(self, Symbol sym, SymbolSpec symSpec) except *:
        cdef list prevGss
        cdef int nShifts
        cdef Gssn topA, topB
        cdef Action action
        cdef ShiftAction sAction
        cdef bint merged

        prevGss = self._gss
        self._gss = []

        if self.verbose:
            nShifts = 0

        for topA in prevGss:
            if symSpec in self.spec._action[topA.nextState]:
                for action in self.spec._action[topA.nextState][symSpec]:
                    if type(action) == ShiftAction:
                        sAction = <ShiftAction>action
                        merged = False
                        for topB in self._gss:
                            if topB.nextState == topA.nextState:
                                Gsse(topA, topB, sym)
                                merged = True
                                break
                        if not merged:
                            top = Gssn(topA, sym, sAction.nextState)
                            self._gss.append(top)
                            if self.verbose:
                                print "   --> shift(a) %d" % sAction.nextState
                                nShifts += 1
        if self.verbose:
            if nShifts > 0:
                self._printStack()

    cdef void _printStack(self):
        cdef int i
        cdef Gssn top
        cdef list path

        i = 0
        for top in self._gss:
            for path in top.paths():
                if i == 0:
                    print "STK 0:",
                else:
                    print "    %d:" % i,
                for elm in path:
                    print "%r" % elm,
                print
                i += 1
