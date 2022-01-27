
# Native library
from os import path
from io import StringIO
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

# Cactus library
from cactus.shared.common import cactus_call

# Sonlib library
from sonLib.nxnewick import NXNewick
from sonLib.nxtree import NXTree

def call_cactus_prepare():
    return

def extract_newick_tree(hal_filename):
    """Extract the newick tree from a HAL-format file using the halStats binary tool

    Args:
        hal_filename (str): The path for the HAL-format file

    Raises:
        RuntimeError: Failed to parse the newick tree extract from the HAL-format file

    Returns:
        sonLib.nxtree.NXTree: The tree employed in the HAL-format file
    """

    # extract the tree in a string format
    string_tree = cactus_call(check_output=True, parameters=[
                              'halStats', '-tree', hal_filename])

    # parse and sanity check
    newickParser = NXNewick()
    try:
        tree = newickParser.parseString(string_tree)
    except:
        raise RuntimeError(
            "Failed to parse newick tree: {}}".format(string_tree))

    return tree


def extract_fasta_files(hal_filename, genomes, hal_options, output_dir):
    """Function to extract a genome in FASTA-format from a HAL-format file using the hal2fasta binary tool

    Args:
        hal_filename (str): The path of the HAL-format file
        genomes (List[str]): List of genomes to be extracted from the HAL-format file
        hal_options (str): The options for customising customise the extraction 
        output_dir (str): The directory path where the FASTA-file will be stored

    Returns:
        [Dict[str,str]]: A dict containing the filename and path for each extracted FASTA-format file
    Raises:
        RuntimeError: [description]
    """
    fasta_paths = {}

    for genome in genomes:
        fasta_filename = path.join(output_dir, genome + '.fa')

        return_code = cactus_call(check_result=True, checkparameters=[
                                  "hal2fasta", hal_filename + hal_options.split(), + ["> ",  fasta_filename]])
        if return_code:
            raise RuntimeError("hal exception caught: Genome {} not found".format(genome))
        
        fasta_paths[genome] = fasta_filename

    return fasta_paths

def get_node_id(tree, node_name):
    """Get the internal node ID in the tree

    Args:
        tree (NXTree): A NXTree-format tree
        node_name (str): The name of a node in the tree

    Returns:
        [int]: The internal node ID in tree if the node_name is found
        [None]: Otherwise
    """
    for node_id in tree.breadthFirstTraversal():
        if tree.getName(node_id) == node_name:
            return node_id

    return None


def get_tree_patch_adding2node(tree, target_genome_id, new_children):
    """[summary]

    Args:
        tree ([type]): [description]
        target_genome_id ([type]): [description]
        new_children ([type]): [description]

    Returns:
        [type]: [description]
    """

   # the tree patch for adding the new genomes as children of the target_genome
    patch = '('

    # re-adding current target_genome children to the patch
    for child_id in tree.getChildren(target_genome_id):
        child_name = tree.getName(child_id)
        patch += '{}:{},'.format(child_name,
                                 tree.getWeight(target_genome_id, child_id))

    # adding the new target_getnome child(ren) to the patch
    for name, weight in new_children.items():
        patch += '{}:{},'.format(name, weight)

    # remove the last ','
    patch = patch.rstrip(patch[-1])

    # close the path and add the target_genome_name
    patch += '){}'.format(tree.getName(target_genome_id))
    if tree.hasParent(target_genome_id):
        # target_genome is not the root of the tree
        patch += ':{};'.format(tree.getWeight(
            tree.getParent(target_genome_id), target_genome_id))
    else:
        # target_genome is the root of the tree
        patch += ';'

    return patch


def adding2node_prepare(target_genome_name, new_genomes, options):
    """A funciton to prepare the addition of new genomes to a node (genome).

    Args:
        hal_filename ([type]): [description]
        target_genome_name ([type]): [description]
        new_genomes ([type]): [description]

    Raises:
        RuntimeError: [description]
    """

    # get the tree
    tree = extract_newick_tree(options.halFile)

    # creating the seqFile as the "tree-patch" for the new alignment due to this update to a node
    # Exemple:
    #   - target_genome_name = genome_3
    #   - new_children = [genome_4]
    #   - string-format patch tree: (genome_1, genome_2, genome_4):genome_3;

    # get node_id of the target_genome_name
    target_genome_id = get_node_id(tree, target_genome_name)

    if target_genome_id is None:
        raise RuntimeError(
            "Genome \'{}\' not found in the tree extract from the HAL-format file {}".format(target_genome_name, options.halFile))

    # get tree patch in newick format
    patch = get_tree_patch_adding2node(tree, target_genome_id, new_genomes)

    # get the target_genome's current children
    children = [ tree.getName(i) for i in tree.getChildren(target_genome_id)]
    fasta_filenames = extract_fasta_files(options.halFile, children, options.halOptions, options.outputDir )

    # write the output
    if options.outSeqFile:
        with open(options.outSeqFile, 'w') as out_sf:
            out_sf.write(patch)
            for name, path in fasta_filenames.items():
                out_sf.write('{}\t{}\n'.format(name, path))


    commands = call_cactus_prepare()
    
    
    

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
