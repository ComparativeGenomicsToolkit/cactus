
## Functions for manipulating sonLib BinaryTree trees

def get_subtree_nodes(tree):
    """
    Gets all the nodes in the tree including tree
    :param tree: BinaryTree
    :return: list of BinaryTree nodes
    """
    return [ tree ] + (get_subtree_nodes(tree.left) if tree.left is not None else []) + \
           (get_subtree_nodes(tree.right) if tree.right is not None else [])

def get_leaves(tree):
    """
    Gets all leaf nodes in the given tree.
    :param tree: BinaryTree
    :return: List of BinaryTree nodes
    """
    return [ i for i in get_subtree_nodes(tree) if i.left is None and i.right is None ]

def get_node(tree, node_name):
    """
    Get a node with the given id==node_name.
    :param tree: BinaryTree
    :param node_name: String
    :return: BinaryTree
    """
    return [ i for i in get_subtree_nodes(tree) if i.iD == node_name ][0]

def get_distances(root, distances={}):
    """
    Get the distances between all pairs of nodes in the tree rooted at "root"
    :param root: The BinaryTree
    :param distances: A dictionary to hold the distances
    :return: distances
    """
    def add_distance(node1, node2, distance):
        distances[(node1, node2)] = distance
        distances[(node2, node1)] = distance

    add_distance(root, root, 0.0) # Add the self distance (corner case)

    def add_distances_for_child_subtree(child):
        get_distances(child, distances) # Recursively add the distances for the child subtree
        for i in get_subtree_nodes(child): # Add the distance from nodes in a child's subtree to the root
            add_distance(i, root, distances[(i, child)] + child.distance)

    if root.left != None:
        add_distances_for_child_subtree(root.left)
        if root.right != None:
            add_distances_for_child_subtree(root.right)
            # Add the distances between the nodes in root.left subtree and root.right subtree
            for i in get_subtree_nodes(root.left):
                for j in get_subtree_nodes(root.right):
                    add_distance(i, j, distances[(i, root.left)] + distances[(j, root.right)] + root.left.distance + root.right.distance)
    elif root.right != None:
        add_distances_for_child_subtree(root.right)

    return distances

def get_event_pairs(tree, events):
    """
    Generates the sequence of event pairs and their pairwise distance in the tree
    :param tree: BinaryTree
    :return: Yields a sequence of (BinaryTree, BinaryTree, distance) tuples
    """
    # Returns all pairs of events in given list and the distance in the the given tree between them
    distances = get_distances(tree)
    for i in range(len(events)):
        for j in range(i+1, len(events)):
            yield events[i], events[j], distances[(events[i], events[j])]
