"""
Moving the progressive decomposition interface here as part of big refactor to get rid of experiment XML and all the various project files that contain them
- parse seqfiles
- compute outgroups
- compute schedule
- compute subtree

These functions should replace all direct access to:
- Project
- MultiCactusProject
- ExperimentWrapper
- Outgroup

They should also reduce a lot of access (especially creation) of:
- SeqFile
- Schedule
- MultiCactusTree

The main goal of this refactor is to get rid of XML files like the cactus Experiment.  They add a terrible amount of cruft, especially when they contain Toil file IDs as strings and things like that.  The various project and wrapper classes are just as bad and need to go too.  It'd be nice to rewrite some of what's left like the MultiCactusTree, Outgroup and Schedule, but not essential especially since we know what's there works.  So instead, this module tries to wrap it all up behind a more usable interface.  
"""

#!/usr/bin/env python3
import logging
import os
import xml.etree.ElementTree as ET
import copy

from cactus.shared.configWrapper import ConfigWrapper
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.progressive.outgroup import GreedyOutgroup, DynamicOutgroup
from cactus.progressive.seqFile import SeqFile
from sonLib.nxnewick import NXNewick

def parse_seqfile(seqfile_path, config_wrapper, root_name = None):
    """
    parse the seqfile
    returns (tree, event->path map, og list (from *'s in seqfile)
    """
    seq_file = SeqFile(seqfile_path)

    mc_tree = MultiCactusTree(seq_file.tree)
    mc_tree.nameUnlabeledInternalNodes(config_wrapper.getDefaultInternalNodePrefix())
    mc_tree.computeSubtreeRoots()

    return (mc_tree, seq_file.pathMap, seq_file.outgroups)

def compute_outgroups(mc_tree, config_wrapper, outgroup_candidates = set(), root_name = None, include_dists = False):
    """
    computes the outgroups
    returns event->outgroups map
    """
    assert isinstance(mc_tree, MultiCactusTree)
    outgroup_candidates = set(outgroup_candidates)

    if not root_name:
        root_name = mc_tree.getRootName()
    
    outgroup = GreedyOutgroup()
    outgroup.importTree(mc_tree, mc_tree.getNodeId(root_name))
    
    if config_wrapper.getOutgroupStrategy() == 'greedy':
        # use the provided outgroup candidates, or use all outgroups
        # as candidates if none are given
        outgroup.greedy(threshold=config_wrapper.getOutgroupThreshold(),
                        candidateSet=outgroup_candidates,
                        candidateChildFrac=config_wrapper.getOutgroupAncestorQualityFraction(),
                        maxNumOutgroups=config_wrapper.getMaxNumOutgroups())
    elif config_wrapper.getOutgroupStrategy() == 'greedyLeaves':
        # use all leaves as outgroups, unless outgroup candidates are given
        outgroup.greedy(threshold=config_wrapper.getOutgroupThreshold(),
                        candidateSet=outgroup_candidates,
                        candidateChildFrac=2.0,
                        maxNumOutgroups=config_wrapper.getMaxNumOutgroups())
    elif config_wrapper.getOutgroupStrategy() == 'greedyPreference':
        # prefer the provided outgroup candidates, if any, but use
        # other nodes as "filler" if we can't find enough.
        outgroup.greedy(threshold=config_wrapper.getOutgroupThreshold(),
                        candidateSet=outgroup_candidates,
                        candidateChildFrac=config_wrapper.getOutgroupAncestorQualityFraction(),
                        maxNumOutgroups=config_wrapper.getMaxNumOutgroups())
        outgroup.greedy(threshold=config_wrapper.getOutgroupThreshold(),
                        candidateSet=None,
                        candidateChildFrac=config_wrapper.getOutgroupAncestorQualityFraction(),
                        maxNumOutgroups=config_wrapper.getMaxNumOutgroups())
    elif config_wrapper.getOutgroupStrategy() != 'none':
        raise RuntimeError('Outgroup strategy "{}" not supported. Must be one of [greedy, greedyLeaves, greedyPreference, none]'.format(config_wrapper.getOutgroupStrategy))

    if not include_dists:
        for k, v in outgroup.ogMap.items():
            outgroup.ogMap[k] = [x[0] for x in v]
            
    return outgroup.ogMap

def get_subtree(mc_tree, root_name, config_wrapper, outgroup_map, include_outgroups=True):
    """
    get the subtree for a given internal node -- this will contain the events and outgroups for a single cactus job
    returns the (multicactus) tree
    
    note: this is used for handling the --root option on the cli to clip the input seqfile at the outset
     (or for getting ingroup and outgroup names of a given sub problem)
    """
    assert isinstance(mc_tree, MultiCactusTree)

    # get the root id
    root_id = mc_tree.getNodeId(root_name)
    
    # compute the subtrees
    # ugly hack to keep track of external outgroups for root experiment (yuck)
    externalOutgroupNames = set()

    # dig out every outgroup
    outgroup_names = set()
    if root_name in outgroup_map:
        for og in outgroup_map[root_name]:
            outgroup_names.add(og)

    # find outgroups we want to extract
    kept_nodes = set(mc_tree.postOrderTraversal(root_id))
    dead_nodes = set()
    external_outgroup_names = set()
    if include_outgroups:
        for node in mc_tree.postOrderTraversal():
            if node not in kept_nodes:
                dead_nodes.add(node)
                name = mc_tree.getName(node)
                if name in outgroup_names:
                    external_outgroup_names.add(name)

    # reroot the tree!
    sub_tree = copy.deepcopy(mc_tree)
    orig_parent = sub_tree.getName(sub_tree.getParent(root_id)) if sub_tree.getParent(root_id) is not None else None
    sub_tree.reroot(root_id)
    
    # add the outgroups to the tree
    # computing distance to new root for each one
    for og_name in external_outgroup_names:
        og_id = sub_tree.getNodeId(og_name)
        dist = 0.
        x = og_id
        while sub_tree.hasParent(x):
            d = sub_tree.getWeight(sub_tree.getParent(x), x)
            if d is None:
                dist = None
                break
            else:
                dist += d
            x = sub_tree.getParent(x)
        sub_tree.addOutgroup(og_name, dist)

    # strip out everything but the immediate children of root
    for node in sub_tree.postOrderTraversal():
        parent = sub_tree.getParent(node)
        if node != root_id and parent != root_id:
            dead_nodes.add(node)

    if orig_parent:
        # after rerooting, the original parent is a child so won't be caught above
        dead_nodes.add(sub_tree.getNodeId(orig_parent))

    # flush out all unused nodes, set the new root, and update the
    # experiment template tree to match the new project tree
    for node in dead_nodes:
        assert sub_tree.hasParent(node)
        sub_tree.removeEdge(sub_tree.getParent(node), node)

    sub_tree.computeSubtreeRoots()    
    return sub_tree

def get_spanning_subtree(mc_tree, root_name, config_wrapper, outgroup_map):
    """
    get the tree that contains all leaf jobs
    """
    # make the multicactus tree
    assert isinstance(mc_tree, MultiCactusTree)

    # get the root id
    root_id = mc_tree.getNodeId(root_name)

    # get the outgroups
    node_id_set = set()
    if root_name in outgroup_map:
        for outgroup_name in outgroup_map[root_name]:
            node_id_set.add(mc_tree.getNodeId(outgroup_name))

    # get the ingroups
    for node_id in mc_tree.getChildren(root_id):
        node_id_set.add(node_id)

    # get the spanning tree
    spanning_tree = mc_tree.extractSpanningTree([mc_tree.getName(node) for node in node_id_set])

    spanning_tree.computeSubtreeRoots()
    return spanning_tree

def get_event_set(mc_tree, config_wrapper, outgroup_map, root_name):
    """
    compute all events we need on hand (ingroups and outgroups) for a given problem or subproblem
    (used to narrow down import to releavnt nodes)
    """
    event_set = set([mc_tree.getName(node) for node in mc_tree.postOrderTraversal()]).union(set(outgroup_map.keys()))
    if root_name:
        # make sure we don't download anything we don't need
        sub_tree = get_subtree(mc_tree, root_name, config_wrapper, outgroup_map)
        tree_events = set([sub_tree.getName(node) for node in sub_tree.postOrderTraversal()])
        event_set = event_set.intersection(tree_events)
        event_set.remove(root_name)
    return event_set
