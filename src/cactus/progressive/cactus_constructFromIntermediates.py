"""
Utility to create a HAL from intermediates left by intermediateResultsUrl.
Requires pypi library "newick".
"""

import sys
from newick import loads, dumps
from argparse import ArgumentParser
from subprocess import check_output, check_call
from copy import deepcopy

def postorder_create(node, prefix, hal):
    if node.name is None:
        raise RuntimeError("Requires a tree with all ancestors labeled.")
    sys.stderr.write("working on node %r\n" % (node.name))
    c2h = prefix + '-' + node.name + '.c2h'
    hal_fa = prefix + '-' + node.name + '.hal.fa'
    # get outgroup list (everything in c2h except children / anc)
    outgroups = []
    for species_line in check_output("grep -E '^s' %s | cut -f 2 | uniq" % c2h, shell=True).splitlines():
        # strip ' marks on either side
        species = species_line[1:-1]
        if species != node.name and species not in [n.name for n in node.descendants]:
            outgroups.append(species)
    # get local newick string
    subtree = deepcopy(node)
    for child in subtree.descendants:
        child.descendants = []
    newick = dumps(subtree)
    # actually perform the addition
    cmd = ['halAppendCactusSubtree', c2h, hal_fa, newick, hal]
    if len(outgroups) > 0:
        cmd.extend(['--outgroups', ",".join(outgroups)])
    sys.stderr.write('Running command %r\n' % cmd)
    check_call(cmd)
    # recurse
    for child in node.descendants:
        if len(child.descendants) == 0:
            # Leaf
            continue
        postorder_create(child, prefix, hal)

if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('tree', help='newick tree (*WITH ANCESTORS LABELED*)')
    parser.add_argument('prefix', help='intermediate prefix, i.e. prefix shared'
                        ' by c2h / .hal.fa files')
    parser.add_argument('output', help='output HAL file')
    opts = parser.parse_args()
    tree = loads(opts.tree)[0]
    postorder_create(tree, opts.prefix, opts.output)
