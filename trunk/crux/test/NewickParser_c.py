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

import sys

print "Test begin"

teststrs = (
    "(A);",
    "(A,B);",
    "(A,B,C);",

    "(A,(B,C));",
    "(A:3,(B:1,C:2):3);",
    "(A:2,B:3);",
    "(A:3,B:1,C:2);"
    )

for str in teststrs:
    print "\n=== %r ===" % str
    try:
        t = crux.Tree.Tree(str, autoMap=True)
        t.canonize()
        print t.render(labels=True, lengths=True)
        t.render(labels=True, lengths=True, outFile=sys.stdout)
    except:
        error = sys.exc_info()
        print "Exception %s: %s" % (error[0], error[1])
        
print "Test end"
