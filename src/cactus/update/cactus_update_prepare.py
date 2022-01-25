from io import StringIO
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from cactus.shared.common import cactus_call

# manipulate Newick-format tree
from sonLib.nxnewick import NXNewick
from sonLib.nxtree import NXTree


# def halStats_call_get_children(hal_file, genome):
#    """Specific call function to get the genome's children using the halStats tool"""
#    children = cactus_call(check_output=True, parameters=[
#                           "halStats", "--children", genome, hal_file])
#    return children.split()


def halStats_call_get_tree(hal_file):
    """Specific call function to get the tree of a HAL-format file using the halStats tool"""
    return cactus_call(check_output=True, parameters=["halStats", "--tree", hal_file])


def hal2fasta_call_get_fasta(hal_file, fasta_file, genome, options="--hdf5InMemory"):
    """Specific call function to transform a HAL-format file into a FASTA-format file using the hal2fasta tool"""
    cactus_call(parameters=["hal2fasta", hal_file,
                genome] + options.split() + ["> ", fasta_file])


@DeprecationWarning
def get_clade(hal_file, target_genome):
    """Get the clade of a target genome given its name"""

    # get the existing tree from the HAL-format file
    tree = Phylo.read(
        StringIO(halStats_call_get_tree(hal_file)), format="newick")

    # get the target clade
    filter = tree.find_elements(target_genome)
    try:
        target_clade = next(filter)
    except StopIteration:
        raise RuntimeError(
            "Genome \'{}\' not found in the tree given by the HAL-format file {}".format(target_genome, hal_file))

    # return
    return target_clade

def extract_newick_tree(hal_file):
    """Extract the newick tree from a HAL-format file

    Args:
        hal_file (str): The path for the HAL-format file

    Raises:
        RuntimeError: Failed to parse the newick tree extract from the HAL-format file

    Returns:
        sonLib.nxtree.NXTree: The tree employed in the HAL-format file
    """

    # extract the tree in a string format
    string_tree = halStats_call_get_tree(hal_file)

    # parse and sanity check
    newickParser = NXNewick()
    try:
        tree = newickParser.parseString(string_tree)
    except:
        raise RuntimeError(
            "Failed to parse newick tree: {}}".format(string_tree))

    return tree

def get_tailored_context(tree, node_name):
    """Get tailored information from from a given node name in order to apply the update procedure

    Args:
        tree (sonLib.nxtree.NXTree): A NXTree-format tree
        node_name (str): The name of a node in the tree

    Returns:
        Dict[str, int]: [description]
    """
    context = {
        'parent':{},
        'children': {}
    }

    for node_id in tree.breadthFirstTraversal():
        if tree.getName(node_id) == node_name:
            
            #store the current node info
            context['parent'][tree.getName(node_id)] = node_id
            
            # also stores the children info
            for child_id in tree.getChildren(node_id):
                context['children'][tree.getName(child_id)] = child_id

    return context


def adding2node_prepare(hal_file, target_genome, new_genomes):
    """Adding new genomes to a target genome using an existing alingment given by a HAL-format file"""

    # get the tree
    tree = extract_newick_tree(hal_file)

    # creating the seqFile as the "tree-patch" for the new alignment due to this update to a node
    # Exemple:
    #   - target_genome = genome_3
    #   - new_genomes = [genome_4]
    #   - string-format patch tree: (genome_1, genome_2, genome_4):genome_3;
    
    current_context = get_tailored_context(tree, target_genome)
    patch_tree = '('



    patch_tree += '{});'.format(target_genome)

# def adding2branch():

# def cactus_update(options):


def main():
    #print("children: {}".format(halStats_call_get_children("/Users/thiagogenez/Documents/buckets/adding-genome/lepidoptera/runs/0/steps/Anc05.hal", "Anc05")))

    # Same main parser as usual
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    # Same subparsers as usual
    subparsers = parser.add_subparsers(
        help='Desired action to perform the alignment update', dest='action')

    # Usual arguments which are applicable for the whole script / top-level args
    parser.add_argument('--verbose', help='Common top-level parameter',
                        action='store_true', required=False)

    # Create parent subparser. Note `add_help=False` and creation via `argparse.`
    parent_parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter, add_help=False)
    parent_parser.add_argument(
        "inHal", help="The input HAL-format file containing the existing alignment")
    parent_parser.add_argument(
        "genomeName", help="The name of the new genome that is being added")
    parent_parser.add_argument(
        "fastaFile", help="The FASTA file of the new genome")

    # taken from cactus-prepare
    parent_parser.add_argument(
        "--halOptions", type=str, default="--hdf5InMemory", help="options for every hal command")
    parent_parser.add_argument(
        "--outDir", help='Directory where the processed leaf sequence and ancestral sequences will be placed.')
    parent_parser.add_argument(
        "--outSeqFile", help="Path for annotated Seq file output [default: outDir/seqFile]")

    # Subparsers based on parent & Add some arguments exclusively for parser_create
    parser_create = subparsers.add_parser("node",  aliases=['no'], parents=[
                                          parent_parser], help='Adding a new genome to a node (aka, the update-node approach)')
    requiredNamed = parser_create.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument(
        "--genome", "-g", help="Name of the genome in the existing alignment", required=True, metavar='GENOME_NAME')

    parser_update = subparsers.add_parser("branch", aliases=['br'], parents=[
                                          parent_parser], help='Add a new genome to a branch (aka, the update-branch approach)')
    requiredNamed = parser_update.add_argument_group(
        'required named arguments')
    requiredNamed.add_argument("--topGenome", "-tg", help="Name of the genome in the existing alignment",
                               dest='parent_genome', required=True, metavar='GENOME_NAME')
    requiredNamed.add_argument("--bottomGenome", "-bg", help="Name of the genome in the existing alignment",
                               dest='child_genome', required=True, metavar='GENOME_NAME')

    options = parser.parse_args()

    print(options)


if __name__ == '__main__':
    main()
