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

print "Test begin"

#
# Generic CharacterType tests.
#

print "=== CharacterType ==="
c = crux.CharacterType.CharacterType()

for code in ['A', 'B', 'C', 'D']:
    print "==> stateCodeAdd(%r)" % code
    c.stateCodeAdd(code)
    val =  c.code2val(code)
    print "%s --> %d" % (code, val)
    print "%s <-- %d" % (c.val2code(val), val)
    print "nstates: %d" % c.nstates()
codes = c.codes()
codes.sort()
print "Codes: %r" % codes

aCodes = {'W': ['A', 'B'],
          'X': ['A', 'B', 'C'],
          'Y': ['C', 'D'],
          'Z': ['A', 'B', 'C', 'D']}
for code in aCodes.keys():
    print "==> ambiguityCodeAdd(%r, %r)" % (code, aCodes[code])
    c.ambiguityCodeAdd(code, aCodes[code])
    val =  c.code2val(code)
    print "%s --> %d" % (code, val)
    print "%s <-- %d" % (c.val2code(val), val)
codes = c.codes()
codes.sort()
print "Codes: %r" % codes

xCodes = {'a': 'A',
          'b': 'B',
          'c': 'C',
          'd': 'D'}
for code in xCodes.keys():
    print "==> aliasCodeAdd(%r, %r)" % (code, xCodes[code])
    c.aliasCodeAdd(code, xCodes[code])
    val =  c.code2val(code)
    print "%s --> %d" % (code, val)
    print "%s <-- %d" % (c.val2code(val), val)
codes = c.codes()
codes.sort()
print "Codes: %r" % codes
print "nstates: %d" % c.nstates()

# Error conditions.
try:
    c.stateCodeAdd('A')
except:
    import sys

    error = sys.exc_info()
    print "Exception %s: %s" % (error[0], error[1])

try:
    c.ambiguityCodeAdd('A', ['B', 'C'])
except:
    import sys

    error = sys.exc_info()
    print "Exception %s: %s" % (error[0], error[1])

try:
    c.aliasCodeAdd('A', 'B')
except:
    import sys

    error = sys.exc_info()
    print "Exception %s: %s" % (error[0], error[1])

try:
    c.aliasCodeAdd('j', 'J')
except:
    import sys

    error = sys.exc_info()
    print "Exception %s: %s" % (error[0], error[1])

try:
    c.code2val('a')
except:
    import sys

    error = sys.exc_info()
    print "Exception %s: %s" % (error[0], error[1])

try:
    c.val2code(42)
except:
    import sys

    error = sys.exc_info()
    print "Exception %s: %s" % (error[0], error[1])

#
# DnaCharacterType tests.
#

print "=== DnaCharacterType ==="
c = crux.CharacterType.DnaCharacterType()
codes = c.codes()
codes.sort()
print "Codes: %r" % codes
print "nstates: %d" % c.nstates()
for code in codes:
    val =  c.code2val(code)
    print "%s --> 0x%x --> %s" % (code, val, c.val2code(val))

#
# ProteinCharacterType tests.
#
print "=== ProteinCharacterType ==="
c = crux.CharacterType.ProteinCharacterType()
codes = c.codes()
codes.sort()
print "Codes: %r" % codes
print "nstates: %d" % c.nstates()
for code in codes:
    val =  c.code2val(code)
    print "%s --> 0x%08x --> %s" % (code, val, c.val2code(val))

print "Test end"
