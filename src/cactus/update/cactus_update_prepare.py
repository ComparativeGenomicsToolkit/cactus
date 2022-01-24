from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from cactus.shared.common import cactus_call


def halStats_call_get_children(inHal_filepath, genome):
    """Specific call function to get the genome's children using the halStats tool"""

    children = cactus_call(check_output=True, parameters=["halStats", "--children", genome, inHal_filepath])
    return children.split()

def hal2fasta_call_get_fasta(inHal_filepath, outFasta_filepath, genome, options="--hdf5InMemory"):
    """Specific call function to transform a HAL-format file into a FASTA-format file using the hal2fasta tool"""
    
    cactus_call(parameters=["hal2fasta", inHal_filepath, genome] + options.split() + ["> ", outFasta_filepath])


#def adding2node():

#def adding2branch():

#def cactus_update(options):

def main(): 
    #print("children: {}".format(halStats_call_get_children("/Users/thiagogenez/Documents/buckets/adding-genome/lepidoptera/runs/0/steps/Anc05.hal", "Anc05")))

    # Same main parser as usual
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    
    # Same subparsers as usual
    subparsers = parser.add_subparsers(help='Desired action to perform the alignment update', dest='action')

    # Usual arguments which are applicable for the whole script / top-level args
    parser.add_argument('--verbose', help='Common top-level parameter', action='store_true', required=False)

    # Create parent subparser. Note `add_help=False` and creation via `argparse.`
    parent_parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter, add_help=False)
    parent_parser.add_argument("inHal", help="The input HAL-format file containing the existing alignment")

    parent_parser.add_argument("--halOptions", type=str, default="--hdf5InMemory", help="options for every hal command")
    parent_parser.add_argument("--outDir", help='Directory where the processed leaf sequence and ancestral sequences will be placed.')
    parent_parser.add_argument("--outSeqFile", help="Path for annotated Seq file output [default: outDir/seqFile]")
    
    # Subparsers based on parent & Add some arguments exclusively for parser_create
    parser_create = subparsers.add_parser("node",  aliases=['no'], parents=[parent_parser], help='Adding a new genome to a node (aka, update-node approach)')
    requiredNamed = parser_create.add_argument_group('required named arguments')
    requiredNamed.add_argument("--genome", help="name of the environment", required=True)

    parser_update = subparsers.add_parser("branch", aliases=['br'], parents=[parent_parser], help='Add a new genome to a branch (aka, update-branch approach)')
    requiredNamed = parser_update.add_argument_group('required named arguments')
    requiredNamed.add_argument("--parentGenome", help="name of the environment", required=True)
    requiredNamed.add_argument("--childGenome", help="name of the environment", required=True)

    options = parser.parse_args()

    print(options)
    
if __name__ == '__main__':
    main()
