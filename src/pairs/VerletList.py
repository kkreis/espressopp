from espresso import pmi
from espresso.esutil import cxxinit

from espresso.pairs.Set import *

from _espresso import pairs_VerletList
class VerletListLocal(SetLocal, pairs_VerletList):
    'Stores a list of pairs generated by Verlet for a pairs::Computer to be applied.'
    def __init__(self, *args, **kywds):

        tuplen = len(args)
        dctlen = len(kywds)
        sumlen = tuplen + dctlen

        if sumlen == 3:
            if tuplen == 0:
                bc=kywds['bc']; storage=kywds['storage']; posProperty=kywds['posProperty']
            elif tuplen == 1:
                bc=args[0]; storage=kywds['storage']; posProperty=kywds['posProperty']
            elif tuplen == 2:
                bc=args[0]; storage=args[1]; posProperty=kywds['posProperty']
            else:
                bc=args[0]; storage=args[1]; posProperty=args[2]
            cxxinit(self, pairs_VerletList, bc, storage, posProperty)

        elif sumlen == 5:
            if tuplen == 0:
                bc=kywds['bc']; storage1=kywds['storage1']; storage2=kywds['storage2']
                posProperty1=kywds['posProperty1']; posProperty2=kywds['posProperty2']
            elif tuplen == 1:
                bc=args[0]; storage1=kywds['storage1']; storage2=kywds['storage2']
                posProperty1=kywds['posProperty1']; posProperty2=kywds['posProperty2']
            elif tuplen == 2:
                bc=args[0]; storage1=args[1]; storage2=kywds['storage2']
                posProperty1=kywds['posProperty1']; posProperty2=kywds['posProperty2']
            elif tuplen == 3:
                bc=args[0]; storage1=args[1]; storage2=args[2]
                posProperty1=kywds['posProperty1']; posProperty2=kywds['posProperty2']
            elif tuplen == 4:
                bc=args[0]; storage1=args[1]; storage2=args[2]
                posProperty1=args[3]; posProperty2=kywds['posProperty2']
            else:
                bc=args[0]; storage1=args[1]; storage2=args[2]
                posProperty1=args[3]; posProperty2=args[4]
            cxxinit(self, pairs_VerletList, bc, storage1, storage2, posProperty1, posProperty2)
        else:
            raise ValueError('Number of arguments to constructor of pairs::List is invalid.')


if pmi.IS_CONTROLLER:
    class VerletList(Set):
        'PMI class of a list of pairs'
        pmiproxydefs = \
            dict(cls = 'espresso.pairs.VerletListLocal', 
                 pmicall = [ 'deletePair', 'size', 'findPair' ])

