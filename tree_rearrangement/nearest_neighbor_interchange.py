from random import seed

from ete4 import Tree

def iter_NNIs(tree, copies=False):
    """
    Generator yielding all neighbor trees using Nearest Neighbor Interchange search.
    
    WARNING: only works with binary trees.
    
    :param False copies: if True yield copies of input Tree. Otherwise (default) it 
       yields alwaysthe same object in different configurations (if stacked in a 
       list, this will finally result in a list of all identical trees)
    """
    for n in tree.iter_descendants():
        if n.is_leaf():
            continue
        # yielder function out of the loop for preformance
        if copies:
            yielder = lambda x: x.copy()
        else:
            yielder = lambda x: x
        # we only need to play with 3 nodes to get all combinations, the fourth can be forgotten
        # find the 3 nodes
        a, b = n.get_children()
        if n.up.is_root():
            if n is tree.get_children()[1]:  # root does not exist here:it's just one edge
                continue
            c = n.get_sisters()[0].get_children()[0]
        else:
            c = n.get_sisters()[0]
        # swap nodes for first neighbor tree
        au = a.up
        bu = b.up
        cu = c.up
        a.detach()
        c.detach()
        au.add_child(c)
        cu.add_child(a)
        yield yielder(tree)
        # rebuild tree
        a.detach()
        c.detach()
        cu.add_child(c)
        au.add_child(a)
        # swap nodes for second neighbor tree
        b.detach()
        c.detach()
        bu.add_child(c)
        cu.add_child(b)
        yield yielder(tree)
        # rebuild tree
        b.detach()
        c.detach()
        cu.add_child(c)
        bu.add_child(b)

seed(2)
t = Tree()
t.populate(5, names_library=map(str, range(1, 11)))
print(t)

for tt in iter_NNIs(t):
    print(tt)

seed(2)
t = Tree()
tlen = 10000
t.populate(tlen, names_library=map(str, range(1, tlen + 1)))

print('Test if the number of neighbor trees found is equal to 2x(N-3), with N the number of leaves:')
print(f'{len(list(iter_NNIs(t)))} == {2 * (tlen - 3)}')
print(len(list(iter_NNIs(t))) == 2 * (tlen - 3))