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

trees = [';',

         'a;',

         '(a:1,b:1);',
         '(a:0,b:0);',

         '(a:1,b:1,c:1);',
         '(a:0,b:1,c:1);',
         '(a:1,b:0,c:1);',
         '(a:1,b:1,c:0);',

         '(a:1,b:1,(c:1,d:1):1);',
         '(a:0,b:1,(c:1,d:1):1);',
         '(a:1,b:0,(c:1,d:1):1);',
         '(a:1,b:1,(c:0,d:1):1);',
         '(a:1,b:1,(c:1,d:0):1);',
         '(a:1,b:1,(c:1,d:1):0);',
         '(a:1,b:1,(c:1,d:1):-1);',

         '(a:1,b:1,(c:1,(d:1,e:1):1):1);',
         '(a:1,b:1,(c:1,(d:1,e:1):0):1);',
         '(a:1,b:1,(c:1,(d:1,e:1):1):0);',
         '(a:1,b:1,(c:1,(d:1,e:1):0):0);',

         '(a:1,b:1,(c:1,d:1,e:1):1);',
         '(a:1,b:1,(c:1,d:1,e:1):0);',
         ]

for tree in trees:
    print "========================================================="
    t = crux.Tree.Tree(tree, newickAutoMap=True)
    t.canonize()
    print t.render(labels=True, lengths=True, lengthFormat="%.0f")
    nCollapsed = t.collapse()
    t.canonize()
    print t.render(labels=True, lengths=True, lengthFormat="%.0f")
    print "Number of collapsed edges: %d" % nCollapsed

print "========================================================="
print "Test end"
