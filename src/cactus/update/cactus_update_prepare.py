
# Native library
import sys
import os
import re
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, SUPPRESS
from tempfile import mkdtemp

# Cactus library
from cactus.shared.common import cactus_call

# Sonlib library
from sonLib.nxnewick import NXNewick
from sonLib.nxtree import NXTree


def call_cactus_prepare(seq_file, out_dir, cactus_prepare_options):

    # if --cactusOptions embedded, it has to been removed and pass as a single string
    pattern = "\s?--cactusOptions [\"'].*?[\"']"
    regex = re.compile(pattern, re.IGNORECASE)
    cactus_options = re.search(pattern, cactus_prepare_options, re.IGNORECASE)
    
    if cactus_options:
        # grab the --cactusOptions content within the given string --cactus-prepare-options & remove quotes
        cactus_options = ['--cactusOptions', cactus_options.group(0).split('--cactusOptions ')[-1].strip("\"'")]
        
        # remove it from --cactus-prepare-options
        cactus_prepare_options = regex.sub(r'', cactus_prepare_options)
    else:

        # nothing found, nothing to pass thought
        cactus_options = []

    cmd = 'cactus-prepare {} --outDir {}/steps --outSeqFile {}/steps/steps.txt --jobStore {}/jobstore {}'.format(
        seq_file, out_dir, out_dir, out_dir,  cactus_prepare_options)

    # call the built-in cactus_call function
    return cactus_call(check_output=True, parameters=cmd.split() + cactus_options)


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
                              'halStats', '--tree', hal_filename])

    # parse and sanity check
    newickParser = NXNewick()
    try:
        tree = newickParser.parseString(string_tree)
    except:
        raise RuntimeError(
            "Failed to parse newick tree: {}}".format(string_tree))

    return tree


def extract_fasta_files(hal_file, genome_names, hal_options, assembly_dir):
    """Function to extract a genome in FASTA-format from a HAL-format file using the hal2fasta binary tool

    Args:
        hal_file (str): The path of the HAL-format file
        genome_names (List[str]): List of names of genome to be extracted from the HAL-format file
        hal_options (str): The options for customising customise the extraction 
        assembly_dir (str): The directory path where the FASTA-file will be stored

    Returns:
        [Dict[str,str]]: A dict containing the filename and path for each extracted FASTA-format file
    Raises:
        RuntimeError: [description]
    """
    fastas = {}

    for genome in genome_names:

        # FASTA-format filename to the correspondent genome
        filename = os.path.join(assembly_dir, genome + '.fa')

        # create the command-line
        # example: hal2fasta steps/Anc78.hal Anc78 --hdf5InMemory > steps/Anc78.fa
        cmd = "hal2fasta {} {} {}".format(hal_file, genome, hal_options)

        # call the built-in cactus_call function
        cactus_call(parameters=cmd.split(), outfile=filename)

        fastas[genome] = os.path.abspath(filename)

    return fastas


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


def get_tree_patch(node_name, parent_weight, children):
    """[summary]

    Args:
        node_name ([type]): [description]
        parent_weight ([type]): [description]
        children ([type]): [description]

    Returns:
        [type]: [description]
    """

   # the tree patch for adding the new genomes as children of the target_genome
    patch = '('

    # adding children to the patch
    for name, weight in children.items():
        patch += '{}:{},'.format(name, weight)

    # remove the last ','
    patch = patch.rstrip(patch[-1])

    # close the path and add the node_name
    patch += '){}'.format(node_name)

    # add the target_node's weight
    if parent_weight:
        # target_genome is not the root of the tree
        patch += ':{}'.format(parent_weight)

    # target_genome is the root of the tree
    patch += ';'

    return patch


def adding2node_prepare(options):
    """A funciton to prepare the addition of new genomes to a node (genome).

    Args:
        hal_filename ([type]): [description]
        target_genome_name ([type]): [description]
        new_genomes ([type]): [description]

    Raises:
        RuntimeError: [description]
    """

    # get the tree
    nxtree = extract_newick_tree(options.inHal)

    # creating the seqFile as the "tree-patch" for the new alignment due to this update to a node
    # Exemple:
    #   - target_genome_name = genome_3
    #   - new_children = [genome_4]
    #   - string-format patch tree: (genome_1, genome_2, genome_4):genome_3;

    # get node_id of the target_genome_name
    genome_internal_id = get_node_id(nxtree, options.genome)

    if genome_internal_id is None:
        raise RuntimeError(
            "Genome \'{}\' not found in the tree extracted from the HAL-format file {}".format(options.genome, options.halFile))

    # get current children of the target genome
    children = {nxtree.getName(child_id): nxtree.getWeight(
        genome_internal_id, child_id) for child_id in nxtree.getChildren(genome_internal_id)}

    # extract FASTA files for each target_genome's children
    sequences = extract_fasta_files(
        options.inHal, children.keys(), options.halOptions, options.work_dir['assemblies_dir'])

    # also add the new children (genomes) and their correspondent assemblies (fasta-format files)
    for genome_name, value in options.inFasta.items():

        # extract data
        (weight, fasta_path) = value

        # update children
        if genome_name in children:
            raise RuntimeError('The target genome \'{}\' already has a child genome named \'{}\'. The new genome being added \'{}\' must have a different name'.format(
                options.genome, genome_name, genome_name))

        # adding the new genome and its weight between itself and its parent
        children[genome_name] = weight

        # update the dict with the fasta-format files of the new genomes
        sequences[genome_name] = fasta_path

    # get the weight between the node and its parent
    parent_weight = None
    if nxtree.hasParent(genome_internal_id):
        weight = nxtree.getWeight(nxtree.getParent(
            genome_internal_id), genome_internal_id)

    # finanly get the newick-formt tree "patch" to embed the update
    patch = get_tree_patch(options.genome, parent_weight, children)

    # write the cactus seqFile for the update aligment
    if options.seq_file:
        with open(options.seq_file, 'w') as out_sf:
            out_sf.write(patch + '\n')
            for name, path in sequences.items():
                out_sf.write('{} {}\n'.format(name, path))

    plan = call_cactus_prepare(
        options.seq_file, options.work_dir['root_dir'], options.cactus_prepare_options)
    
    print(plan)


def cactus_aligment_update(options):
    print(options)
    if 'node' == options.action:
        adding2node_prepare(options)


def add_subcommand_options(subparser, parent_parser, subcommand):

    if 'node' in subcommand:
        parser_node_approach = subparser.add_parser("node", parents=[
                                                    parent_parser], help='Adding a new genome to a node (aka, the update-node approach)')
        requiredNamed = parser_node_approach.add_argument_group(
            'required named arguments')

        # required args for subcommand "node"
        requiredNamed.add_argument(
            "--genome", "-g", help="Name of the genome in the existing alignment", required=True, metavar='GENOME_NAME')

    elif 'branch' in subcommand:

        parser_branch_approach = subparser.add_parser("branch", parents=[
                                                      parent_parser], help='Add a new genome to a branch (aka, the update-branch approach)')
        requiredNamed = parser_branch_approach.add_argument_group(
            'required named arguments')

        # required args for subcommand "branch"
        requiredNamed.add_argument("--topGenome", "-tg", help="Name of the genome in the existing alignment",
                                   dest='parent_genome', required=True, metavar='GENOME_NAME')
        requiredNamed.add_argument("--bottomGenome", "-bg", help="Name of the genome in the existing alignment",
                                   dest='child_genome', required=True, metavar='GENOME_NAME')


def inHal_sanity_check(filename):
    filename = os.path.abspath(filename)

    return filename


def work_dir_sanity_check(dir):

    # defince base directories
    work_dirs = {
        'root_dir': dir,
        'assemblies_dir': os.path.join(dir, 'assemblies')
    }

    # create them
    for sub_dir in work_dirs.values():
        if sub_dir != dir:
            os.makedirs(sub_dir, exist_ok=True)

    return work_dirs


def inFasta_sanity_check(filename):
    inFasta = {}

    with open(filename) as f:
        for line in f:
            line = line.split()

            # TODO: check if len(line) !=3
            if len(line) == 2:
                line.append(1.0)

            # TODO: validate name and weight
            (name, path, weight) = line[0], line[1], line[2]

            # better sanity check needed for name and weight
            if path.startswith('s3://') or path.startswith('http://'):
                raise RuntimeError(
                    "http:// and s3:// path are not supported yet for sequence name \'{}\' ".format(name))

            if not os.path.isfile(path):
                raise RuntimeError(
                    "Invalid path for the sequence name \'{}\': {}".format(name, path))

            inFasta[name] = (weight, os.path.abspath(path))

    return inFasta


def main():
    #print("children: {}".format(halStats_call_get_children("/Users/thiagogenez/Documents/buckets/adding-genome/lepidoptera/runs/0/steps/Anc05.hal", "Anc05")))

    # Same main parser as usual
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter,
                            usage='%(prog)s {node,branch} [-h] [Options]', add_help=False)
    parser._positionals.title = "subcommand"

    # Same subparsers as usual
    subparser = parser.add_subparsers(
        help='Desired alignment update approach', dest='action')

    # Usual arguments which are applicable for the whole script / top-level args
    # parser.add_argument('--verbose', help='Common top-level parameter',
    #                    action='store_true', required=False)

    # hack to add --help
    parser.add_argument('--help', '-h', action="store_true",
                        dest='help', help=SUPPRESS)

    # Create parent subparser. Note `add_help=False` and creation via `argparse.`
    parent_parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter, add_help=False)
    parent_parser.add_argument(
        "inHal", help="The input HAL-format file containing the existing alignment")
    parent_parser.add_argument(
        "inFasta", help="Tab-separated file. First column: genome name, second column: FASTA-format file path")

    # taken from cactus-prepare
    parent_parser.add_argument(
        "--halOptions", type=str, default="--hdf5InMemory", help="options for every hal command")
    parent_parser.add_argument(
        "--outDir",  default='.', help='Directory where assemblies and cactus-call dependencies will be placed.', dest='work_dir')
    parent_parser.add_argument(
        "--skip-halValidate", help="Skip the validation of the given HAL file", action="store_false", dest='skip_halValidate')
    parent_parser.add_argument("--cactus-prepare-options", dest='cactus_prepare_options',
                               type=str, default="--preprocessBatchSize 1", help="Options to override the default behaviour of cactus-prepare")

    # add subcommands options
    add_subcommand_options(subparser, parent_parser, 'node')
    add_subcommand_options(subparser, parent_parser, 'branch')

    options = parser.parse_args()

    if not options.action:
        parser.print_help()
        sys.exit(1)

    # validate the given inFasta text file
    options.inFasta = inFasta_sanity_check(options.inFasta)

    # validate the given inHal file
    options.inHal = inHal_sanity_check(options.inHal)

    # validate --outDir
    options.work_dir = work_dir_sanity_check(options.work_dir)

    # cactus input seqFile is defined automatically
    options.seq_file = os.path.join(
        options.work_dir['root_dir'], 'seqFile.txt')

    cactus_aligment_update(options)


if __name__ == '__main__':
    main()
