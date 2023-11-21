# Native library
import sys
import os
import re
import shutil
from argparse import (
    ArgumentParser,
    ArgumentDefaultsHelpFormatter,
    SUPPRESS,
    ArgumentTypeError,
)

# Toil library
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options, add_logging_options

# Sonlib library
from sonLib.nxnewick import NXNewick

# Cactus library
from cactus.shared.common import cactus_call
from cactus.shared.common import setupBinaries
from cactus.progressive.cactus_prepare import get_generation_info


def check_positive_float(value):
    """Helper function for argparse to check if a number is positive number

    Args:
        value (float): A given number parameter

    Raises:
        ArgumentTypeError: if the given number is negative

    Returns:
        float: float-format of the given number
    """
    fValue = float(value)
    if fValue <= 0:
        raise ArgumentTypeError(f"{value} is an invalid positive number")
    return fValue


def create_seq_file(seq_file, tree, assemblies):
    """Creates the seq_file for cactus-prepare

    Args:
        seq_file (str): seq_file filename
        tree (str): a Newick-format tree
        assemblies (Dict[str, str]): A dictionary containing the name and path for each sequence in the tree
    """

    # sanity check
    assert seq_file

    # write the cactus seqFile for the update alignment
    with open(seq_file, "w") as out_sf:
        out_sf.write(tree + "\n")
        for name, path in assemblies.items():
            out_sf.write(f"{name} {path}\n")


def call_cactus_prepare(seq_file, outDir, jobStore, outSeqFile, cactus_prepare_options):
    """Helper function to call cactus-prepare using built-in cactus_call

    Args:
        seq_file (str): the path for the seq_file file used by cactus-prepare
        outDir (str): the path for the output directory
        jobStore (str): the path for the jobstore directory
        outSeqFile (str): the path for the outSeqFile cactus-prepare parameter
        cactus_prepare_options (str): extra options to input into the cactus-prepare all

    Returns:
        str: the output of cactus-prepare
    """

    # if --cactusOptions embedded, it has to been removed to be deliveried as
    # a single string
    pattern = "\\s{0,}-{2}cactusOptions\\s{1,}[\"'].*?[\"']"
    regex = re.compile(pattern, re.IGNORECASE)
    cactus_options = re.search(pattern, cactus_prepare_options, re.IGNORECASE)

    if cactus_options:
        # grab the --cactusOptions content within the given string
        # --cactus-prepare-options & remove quotes
        cactus_options = [
            "--cactusOptions",
            " ".join(
                cactus_options.group(0)
                .strip()
                .split("--cactusOptions ")[-1]
                .strip()  # remove extra spaces
                .strip("\"'")  # break using quotes as separator
                .strip()  # remove extra spaces
                .split()  # final list
            ),
        ]

        # remove it from --cactus-prepare-options
        cactus_prepare_options = regex.sub(r"", cactus_prepare_options)
    else:

        # nothing found, nothing to pass thought
        cactus_options = []

    # cactus_prepare_options is append at the tail of cmd to overide cactus
    # options we pre-defined below
    cmd = (
        f"cactus-prepare {seq_file} "
        f"--outDir {outDir} "
        f"--outSeqFile {outSeqFile} "
        f"--jobStore {jobStore} "
        f"{cactus_prepare_options}"
    )

    # call the built-in cactus_call function
    return cactus_call(check_output=True, parameters=cmd.split() + cactus_options)


def extract_newick_tree(hal_file):
    """Extracts the newick tree from a HAL-format file using the halStats binary tool

    Args:
        hal_file (str): The path for the HAL-format file

    Raises:
        RuntimeError: Failed to parse the newick tree extract from the HAL-format file

    Returns:
        sonLib.nxtree.NXTree: The tree employed in the HAL-format file
    """

    # extract the tree in a string format
    string_tree = cactus_call(
        check_output=True, parameters=["halStats", "--tree", hal_file]
    )

    # parse and sanity check
    newickParser = NXNewick()
    try:
        tree = newickParser.parseString(string_tree)
    except BaseException as error:
        raise RuntimeError(f"Failed to parse newick tree: {string_tree}") from error

    return tree


def extract_fasta_files(hal_file, genome_names, halOptions, assembly_dir):
    """Extracts a genome in FASTA-format from a HAL-format file using the hal2fasta binary tool

    Args:
        hal_file (str): The path of the HAL-format file
        genome_names (List[str]): List of names of genome to be extracted from the HAL-format file
        halOptions (str): The options for customizing  the extraction
        assembly_dir (str): The directory path where the FASTA-file will be stored

    Returns:
        [Dict[str,str]]: A dict containing the filename and path for each extracted FASTA-format file
    """
    fastas = {}

    for genome in genome_names:

        # FASTA-format filename to the correspondent genome
        filename = os.path.join(assembly_dir, genome + ".fa")

        # create the command-line
        # example: hal2fasta steps/Anc78.hal Anc78 --hdf5InMemory >
        # steps/Anc78.fa
        cmd = f"hal2fasta {hal_file} {genome} {halOptions}"

        # call the built-in cactus_call function
        cactus_call(parameters=cmd.split(), outfile=filename, check_output=False)

        fastas[genome] = os.path.relpath(filename)

    return fastas


def get_node_id(tree, node_name):
    """Gets the internal node ID in the tree

    Args:
        tree (NXTree): A NXTree-format tree
        node_name (str): The name of a node in the tree

    Returns:
        [int]: The internal node ID in tree if the node_name is found
        [None]: Otherwise

    Raises:
        RuntimeError: if a node_name is not found in the tree
    """
    for node_id in tree.breadthFirstTraversal():
        if tree.getName(node_id) == node_name:
            return node_id

    # sanity check
    raise RuntimeError(
        f"Genome {node_name} not found in the tree extracted from the given HAL-format file"
    )


def get_tree_patch(node_name, top_weight, children_length_map, close=True):
    """ "Creates a tree in a NEWICK format

    Args:
        node_name (str): the name of the current node
        top_weight (float): Edge weight connecting the current node with its parent
        children_length_map (Dict[str,float]): Node's children
        close (bool, optional): Indicator to close the tree string with `;`. Defaults to True.

    Returns:
        str: the newick-format tree
    """

    # the tree patch for adding the new genomes as children of the target_genome
    patch = "("

    # adding children to the patch
    for name, weight in children_length_map.items():
        patch += f"{name}:{weight},"

    # remove the last ','
    patch = patch.rstrip(patch[-1])

    # close the path and add the node_name
    patch += f"){node_name}"

    # add the target_node's weight
    if top_weight:
        # target_genome is not the root of the tree
        patch += f":{top_weight}"

    # target_genome is the root of the tree
    if close:
        patch += ";"

    return patch


def remove_unnecessary_cactus_preprocess(plan, input_names):
    """Removes cactus-preprocess jobs from the given plan

    Args:
        plan (str): a plan given by cactus-prepare
        input_names (List[str]): A list of sequence names that doesn't to be preprocessed

    Returns:
        str: plan without some cactus-preprocess jobs
    """
    plan = plan.split("\n")
    edited_plan = []

    idx = 0
    for idx, line in enumerate(plan):

        # job done as the rest of the plan is the from this point
        if "## Alignment" in line:
            break

        input_name_pattern = "-{2}inputNames\\s{0,}.*?\\s{0,1}(?=-{2}|$)"

        # search for --inputNames where unnecessary cactus_preprocess must be removed
        search = re.search(input_name_pattern, line, re.IGNORECASE)

        # if there is a hit
        if search:
            # get the inputNames list
            existing_input_names = (
                search.group(0).strip("--").strip("inputNames").strip().split()
            )

            # remove the names that already exists
            cleaned = set(existing_input_names) - set(input_names)

            # it becomes an useless cactus-preprocess job and it must be removed (no copy)
            if len(cleaned) == 0:
                continue

            # modify the current line
            line = re.sub(
                input_name_pattern,
                "--inputNames " + " ".join(cleaned) + " ",
                line,
            ).strip()  # remove extra spaces

        # just copy as normal
        edited_plan.append(line)

    # copy the rest of the plan
    edited_plan = edited_plan + plan[idx:-1]

    # job done
    return "\n".join(edited_plan)


def make_plan_amendments(plan, update_cmds, validation_cmds):
    """Amends the given plan with updates and validation command lines

    Args:
        plan (str): original parsed plan from cactus-prepare
        update_cmds (List[str]): list of halGenomeReplace or halAddToBranch command lines
        validation_cmds (List[str]): list of halValidade command lines

    Returns:
        str: the plan with update and validation command lines included
    """

    # add update commands
    plan += "\n## Alignment update\n"
    plan += "\n".join(update_cmds)

    # add validation commands
    plan += "\n\n## Alignment validation\n"
    plan += "\n".join(validation_cmds)

    # add a last empty line
    plan += "\n"

    # amendments done
    return plan


def make_plan(
    seq_file,
    outDir,
    jobStore,
    outSeqFile,
    cactus_prepare_options,
    patch,
    assemblies,
    preprocess_to_remove,
):
    """Creates the plan from cactus-prepare and performs needed amendments removing unnecessary command lines

    Args:
        seq_file (str): path for the seqFile used by cactus-prepare
        outDir (str): path for the output directory
        jobStore (str): path for the jobstore directory
        outSeqFile (str): path for outSeqFile used by cactus-prepare
        cactus_prepare_options (_type_): _description_
        patch (str): the newick tree to be included in the seq_file
        assemblies (Dict[float,str]): assemblies/sequencies information to be included in the seq_file
        preprocess_to_remove (List[str]): list of sequences that do not need preprocessing step

    Returns:
        str: the stepwise execution plan to execute Cactus
    """

    # create the seqFile used as an input in cactus-prepare
    create_seq_file(seq_file, patch, assemblies)

    # created new header
    header = get_generation_info()

    # run cactus_prepare
    plan = call_cactus_prepare(
        seq_file, outDir, jobStore, outSeqFile, cactus_prepare_options
    )

    # remove old header
    plan = re.sub("## generated by", "## wrapping", plan, re.IGNORECASE | re.MULTILINE)
    plan = re.sub("## date .*?\n", "", plan, re.IGNORECASE | re.MULTILINE)
    plan = re.sub("## cactus commit .*?\n", "", plan, re.IGNORECASE | re.MULTILINE)

    # attach new header
    plan = header + plan

    # clean up cactus-preprocess jobs for existing children
    plan = remove_unnecessary_cactus_preprocess(plan, preprocess_to_remove)

    # add --includeRoot option into Round #1's cactus-blast and cactus-align jobs only
    # the --includeRoot option includes the root's sequence in the alignment
    for cactus_job in re.findall(
        r"cactus-(?:blast|align) .*?(?=\s{0,}\n)",
        plan,
        re.IGNORECASE | re.MULTILINE,
    )[-2:]:
        plan = re.sub(
            cactus_job,
            cactus_job + " --includeRoot",
            plan,
            re.IGNORECASE | re.MULTILINE,
        )

    # removing "## HAL merging" as there is no merging while performing alignment updates
    plan = re.sub("\n## HAL merging\n", "", plan, re.IGNORECASE | re.MULTILINE)
    plan = re.sub("halAppendSubtree .*?\n", "", plan, re.IGNORECASE | re.MULTILINE)

    # remove the last hal2fasta as it will
    hal2fasta_job = re.findall(
        r"hal2fasta .*?(?=\s{0,})\n", plan, re.IGNORECASE | re.MULTILINE
    )[-1]
    plan = re.sub(hal2fasta_job, "", plan, re.IGNORECASE | re.MULTILINE)

    return plan


def get_plan_adding2node(
    genome,
    in_hal,
    in_fasta_map,
    seq_file,
    halOptions,
    outDir,
    jobStore,
    outSeqFile,
    cactus_prepare_options="",
    remove_genome=None
):
    """
    A function that automatises the instructions at
    https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/updating-alignments.md#adding-to-a-node

    """

    # get the tree
    nxtree = extract_newick_tree(in_hal)

    # get node_id of the target_genome_name
    genome_internal_id = get_node_id(nxtree, genome)

    # get current children of the target genome
    length_map = {
        nxtree.getName(child_id): nxtree.getWeight(genome_internal_id, child_id)
        for child_id in nxtree.getChildren(genome_internal_id)
    }

    # extract FASTA files for each target_genome's children and the target_genome itself
    assemblies = extract_fasta_files(
        in_hal, [*length_map] + [genome], halOptions, outDir
    )

    # the remains new ancestor's children are given by the in_fasta_map file
    for genome_name, value in in_fasta_map.items():
        # extract weight and path
        (length_map[genome_name], assemblies[genome_name]) = value

    # The weight of the genome's first ancestor must be set to None here;
    # otherwise, cactus-prepare will "wrap" the current genome with a
    # brand-new ancestor. Also, the real weight between the current genome
    # and its first ancestor doesn't need to be catch now. The new genome
    # alignment will be replace anyway during the halReplaceGenome step,
    # and this weight is already in the given hal_file
    top_weight = None

    # finally get the newick-format tree "patch" to embed the update
    patch = get_tree_patch(genome, top_weight, length_map, True)

    # filename of the last hal file
    out_hal = f"{os.path.join(outDir,genome)}.hal"

    # create execution plan using cactus_prepare
    plan = make_plan(
        seq_file,
        outDir,
        jobStore,
        outSeqFile,
        cactus_prepare_options + " --outHal " + out_hal,
        patch,
        assemblies,
        list(set([*length_map]) - set([*in_fasta_map])),
    )

    # Amendment commands
    update_cmds = [
        (
            # existing tree
            f"halReplaceGenome "
            # hal file containing an alignment of the genome and its children.
            f"--bottomAlignmentFile {out_hal} "
            # hal file containing an alignment of the genome, its parent, and its siblings.
            f"--topAlignmentFile {in_hal} "
            # existing tree
            f"{in_hal} "
            # name of genome to be replaced
            f"{genome} "
            f"{halOptions} "
        )
    ]
    if remove_genome:
        update_cmds.insert(0, f"halRemoveGenome {out_hal} {remove_genome}")

    # adding HAL-file validating instructions
    validation_cmds = [f"halValidate --genome {genome} {in_hal} {halOptions}"]

    # finally create the plan
    plan = make_plan_amendments(plan, update_cmds, validation_cmds)

    # alignment-patch done!
    return plan


def get_plan_adding2branch(
    parent_genome,
    child_genome,
    new_ancestor_name,
    in_hal,
    in_fasta_map,
    seq_file,
    halOptions,
    outDir,
    jobStore,
    outSeqFile,
    top_length,
    bottom_length,
    cactus_prepare_options="",
):
    """A function that automatises the instructions at
    https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/updating-alignments.md#adding-to-a-branch
    """

    # get the tree
    nxtree = extract_newick_tree(in_hal)

    # get the node_id of parent_genome
    parent_genome_internal_id = get_node_id(nxtree, parent_genome)

    # get the node_id of child_genome
    child_genome_internal_id = get_node_id(nxtree, child_genome)

    # sanity check
    if not child_genome_internal_id in nxtree.getChildren(
        parent_genome_internal_id
    ) or not parent_genome_internal_id == nxtree.getParent(child_genome_internal_id):
        raise RuntimeError(
            f'Not a valid branch. Top genome "{parent_genome}" is not a parent of the bottom genome"{child_genome}"'
        )

    #########
    # PREAMBLE:  Length check
    #########
    original_length = nxtree.getWeight(
        parent_genome_internal_id, child_genome_internal_id
    )

    if not bottom_length:
        # raise runtime error indicating to use --forceBottomBranchLength option
        if top_length > original_length:
            raise RuntimeError(
                f'The given value "{top_length}" via --topBranchLength '
                f"desired for the branch length between {parent_genome} "
                f"(top genome) and the new ancestor ({new_ancestor_name}) "
                f"is higher than the original branch length "
                f"({original_length}) between {parent_genome} (top genome) "
                f"and {child_genome} (bottom genome). "
                f"This would imply in a negative length "
                f"({original_length} - {top_length} = {original_length - top_length}) to be used between "
                f"the new ancestor ({new_ancestor_name}) and {child_genome} "
                f"(bottom genome). Use --forceBottomBranchLength to make a "
                f"compulsory positive length value."
            )

        # calculates the bottom length as a difference between lengths
        # https://github.com/ComparativeGenomicsToolkit/hal/blob/ab7d889a31c99d74d477132eb136a2f89e859654/api/hdf5_impl/hdf5Alignment.cpp#L272
        if top_length <= original_length:
            bottom_length = original_length - top_length

    elif bottom_length:
        logger.warning('The use "--forceBottomBranchLength" has been enforced')

        if bottom_length + top_length != original_length:
            logger.warning(
                f"{bottom_length} ({child_genome}-{new_ancestor_name} "
                f"branch length) + {top_length} ({new_ancestor_name}-{parent_genome} "
                f"branch length) != {original_length} (original {child_genome}-{parent_genome} "
                f"branch length). Please, after the HAL-format {in_hal} file have been updated, "
                f'run "halUpdateBranchLengths" to '
                f"update the branch lengths for the correct NEWICK tree display by HAL."
            )

    #########
    # STEP 1) Bottom half step:  inferring a new ancestor
    #########

    # one of new ancestor's children is the given bottom genome that is already presented
    # in the given alignment (HAL file)
    length_map = {nxtree.getName(child_genome_internal_id): bottom_length}

    # extract FASTA files for each children capture above
    assemblies = extract_fasta_files(in_hal, length_map.keys(), halOptions, outDir)

    # others new ancestor's children are given by the in_fasta_map (given file)
    for genome_name, value in in_fasta_map.items():
        # extract weight and the path for the fasta file
        (length_map[genome_name], assemblies[genome_name]) = value

    # infer the tree patch for the new ancestor
    patch = get_tree_patch(new_ancestor_name, None, length_map, False)

    #####
    # STEP 2) TOP half step:  addressing new ancestor's parent
    #########

    # Get the parent_genome's children (with the except of the bottom genome addressed above) because they will:
    #   - a) be the new ancestor's siblings
    #   - b) provide outgroup information to build the new ancestor more accurately
    length_map = {
        nxtree.getName(child_id): nxtree.getWeight(parent_genome_internal_id, child_id)
        for child_id in nxtree.getChildren(parent_genome_internal_id)
        if child_id != child_genome_internal_id
    }

    # extract FASTA files for each parent_genome's children (except the bottom genome) + the parent_genome itself
    assemblies = {
        **assemblies,
        **extract_fasta_files(
            in_hal, [*length_map] + [parent_genome], halOptions, outDir
        ),
    }

    # update the map adding the length connecting the patch (previously created) with its parent (i.e. new ancestor)
    length_map[patch] = top_length

    # update the tree's patch to address the new ancestor's parent
    patch = get_tree_patch(parent_genome, None, length_map, True)

    # filename of the last hal file
    top_half_hal = f"{os.path.join(outDir,parent_genome)}.hal"

    # create execution plan
    plan = make_plan(
        seq_file,
        outDir,
        jobStore,
        outSeqFile,
        cactus_prepare_options + " --outHal " + top_half_hal,
        patch,
        assemblies,
        [child_genome] + list(set([*length_map]) - set([patch])),
    )

    # Amendment commands

    # HACK: in case a list of genomes is delivered pick random one in the case of existing more than one
    new_leaf = next(iter(in_fasta_map))
    update_cmds = [
        (
            f"halAddToBranch "
            f"{in_hal} "  # existing tree
            # tree containing insert, its proper bottom segments, and the new leaf genome
            f"{os.path.join(outDir,new_ancestor_name)}.hal "
            f"{top_half_hal} "  # tree containing insert, its parent, and its proper top segments
            f"{parent_genome} "  # insert's future parent
            f"{new_ancestor_name} "  # insert name
            f"{child_genome } "  # insert's future child
            f"{new_leaf} "  # name of new leaf genome
            f"{top_length } "  # length of branch from parent to insert
            f"{in_fasta_map[new_leaf][0]} "  # leaf branch length
            f"{halOptions} "
        )
    ]
    # HACK: to tackle the above issue
    if len(in_fasta_map) > 1:
        update_cmds.append(
            (
                f"halReplaceGenome "  # existing tree
                # hal file containing an alignment of the genome and its children.
                f"--bottomAlignmentFile {os.path.join(outDir,new_ancestor_name)}.hal "
                # hal file containing an alignment of the genome, its parent, and its siblings.
                f"--topAlignmentFile {in_hal} "
                f"{in_hal} "  # existing tree
                f"{new_ancestor_name} "  # name of genome to be replaced
                f"{halOptions} "
            )
        )

    # adding HAL-file validating instructions
    hal_validate_cmds = [
        f"halValidate --genome {i} {in_hal} {halOptions}"
        for i in [parent_genome, new_ancestor_name, child_genome]
    ]
    hal_validate_cmds.extend(
        [f"halValidate --genome {i} {in_hal} {halOptions}" for i in in_fasta_map.keys()]
    )

    # finally create the plan
    plan = make_plan_amendments(plan, update_cmds, hal_validate_cmds)

    # alignment-patch done!
    return plan


def cactus_alignment_update(options):
    """Calls the proper update approach: node or branch

    Args:
        options (Namespace): Namespace object containing the parameter values

    Raises:
        RuntimeError: if action is not node nor branch
    """

    # A genome can be replaced (for example to update an assembly version) by
    # removing it and then following the "add-to-node" procedure to add the
    # new version back as a child of its parent.
    remove_genome = None
    if "replace" == options.action:

        # grab parent_genome info before genome removal
        parent_genome = cactus_call(
            check_output=True,
            parameters=["halStats", "--parent", options.genome, options.in_hal],
        ).strip()

        # grab branch_length info before genome removal
        branch_length = cactus_call(
            check_output=True,
            parameters=[
                "halStats",
                "--branchLength",
                options.genome,
                options.in_hal,
            ],
        ).strip()

        # keep the same length
        options.in_fasta_map[options.genome] = (
            branch_length,
            *options.in_fasta_map[options.genome][1:],
        )

        # remove the genome for replacement
        remove_genome = options.genome

        # following the "add-to-node" procedure below to add it again
        options.action = "add"
        options.subcommand = "node"
        options.genome = parent_genome

    if "add" == options.action:
        if "node" == options.subcommand:
            print(
                get_plan_adding2node(
                    options.genome,
                    options.in_hal,
                    options.in_fasta_map,
                    options.inSeqFile,
                    options.halOptions,
                    options.outDir,
                    options.jobStore,
                    options.outSeqFile,
                    options.cactus_prepare_options,
                    remove_genome=remove_genome
                )
            )

        elif "branch" == options.subcommand:

            # create a name for the new ancestor if not given one
            if not options.ancestor_name:
                options.ancestor_name = (
                    f"{options.child_genome}-Patch-{options.parent_genome}"
                )

            print(
                get_plan_adding2branch(
                    options.parent_genome,
                    options.child_genome,
                    options.ancestor_name,
                    options.in_hal,
                    options.in_fasta_map,
                    options.inSeqFile,
                    options.halOptions,
                    options.outDir,
                    options.jobStore,
                    options.outSeqFile,
                    options.top_length,
                    options.bottom_length,
                    options.cactus_prepare_options,
                )
            )


def subcommand_options(subparser, parent_parser, subcommand, action):
    """Adds the node and branch action to the main argparse

    Args:
        subparser (ArgumentParser): subparser object from the main ArgumentParser object
        parent_parser (ArgumentParser object): the parent of subparser
        action (str): action name
    """

    if "replace" in subcommand:
        parser_node_approach = subparser.add_parser(
            "replace",
            parents=[parent_parser],
            help="Replacing an existing genome",
        )
        requiredNamed = parser_node_approach.add_argument_group(
            "Replace genome options"
        )
        # required args for action "node"
        requiredNamed.add_argument(
            "--genome",
            help="Name of the genome in the existing alignment",
            required=True,
            metavar="GENOME_NAME",
        )

    elif "add" in subcommand:
        if "node" in action:
            parser_node_approach = subparser.add_parser(
                "node",
                parents=[parent_parser],
                help="Adding a new genome to a node (aka, the adding-to-a-node approach)",
            )
            requiredNamed = parser_node_approach.add_argument_group(
                "Adding-to-a-node options"
            )

            # required args for action "node"
            requiredNamed.add_argument(
                "--genome",
                help="Name of the genome in the existing alignment",
                required=True,
                metavar="GENOME_NAME",
            )

        elif "branch" in action:
            parser_branch_approach = subparser.add_parser(
                "branch",
                parents=[parent_parser],
                help="Adding a new genome to a branch (aka, the adding-to-a-branch approach)",
            )
            requiredNamed = parser_branch_approach.add_argument_group(
                "Adding-to-a-branch options"
            )

            # required args for action "branch"
            requiredNamed.add_argument(
                "--parentGenome",
                help="Name of the genome in the existing alignment",
                dest="parent_genome",
                required=True,
                metavar="GENOME_NAME",
            )
            requiredNamed.add_argument(
                "--childGenome",
                help="Name of the genome in the existing alignment",
                dest="child_genome",
                required=True,
                metavar="GENOME_NAME",
            )
            requiredNamed.add_argument(
                "--ancestorName",
                help="Name of the new inferred ancestor",
                dest="ancestor_name",
                metavar="GENOME_NAME",
            )
            requiredNamed.add_argument(
                "--topBranchLength",
                help="Length of the branch between the new ancestor and the top genome",
                dest="top_length",
                type=check_positive_float,
                metavar="BRANCH_LENGTH",
                default=1.0,
            )
            requiredNamed.add_argument(
                "--forceBottomBranchLength",
                help="Forcing the branch length between the new ancestor and the bottom genome",
                dest="bottom_length",
                type=check_positive_float,
                metavar="BRANCH_LENGTH",
            )


def sanity_checks(options):
    def in_hal_sanity_check():
        """Makes backups (if requested) and validates the given HAL-format file

        Raises:
            RuntimeError: if the given path doesn't point to a file
        """

        # is it a valid path?
        if not os.path.isfile(options.in_hal):
            raise RuntimeError(f"Invalid file: '{options.in_hal}'")

        # skip halValidate?
        if not options.skip_halValidate:
            cactus_call(check_output=True, parameters=["halValidate", options.in_hal])

        # skip backup?
        if not options.skip_backup:
            # make a copy
            shutil.copy2(options.in_hal, f"{options.in_hal}.bak", follow_symlinks=True)

    def cactus_prepare_options_sanity_check():
        """Fixes the embedded cactus-prepare option removing the parameters
        outDir, outSeqFile and jobStore to avoid conflicts with these parameters
        used by this run (cactus-update-prepare)
        """

        # Strip some options embedded in cactus_prepare_options to avoid confusions
        keywords = ["outDir", "outSeqFile", "jobStore"]

        for kw in keywords:
            options.cactus_prepare_options = re.sub(
                f"\\s{{0,}}-{{2}}{kw}\\s{{1,}}.*?\\s{{0,}}(?=-{{2}}|$)",
                " ",
                options.cactus_prepare_options,
                re.IGNORECASE,
            )

    def work_dir_sanity_check():
        """Creates the jobstore and output directories"""

        # just create them
        os.makedirs(options.jobStore, exist_ok=True)
        os.makedirs(options.outDir, exist_ok=True)

    def in_fasta_map_sanity_check():
        """Checks the given tab-separated text file"""
        options.in_fasta_map = {}

        with open(options.in_fasta) as f:
            for line in f:
                line = line.split()

                if len(line) < 2 or len(line) > 3:
                    raise RuntimeError(
                        f"Tab-separated text file must have 2 or 3 (tab-separated) columns per row"
                    )

                # no length give, add the default value of 1.0
                if len(line) == 2:
                    line.append(1.0)

                # TODO: validate name and weight
                (name, path, weight) = line[0], line[1], line[2]

                # check if path is not a url
                if not re.search("^(https?|s3)", path, re.IGNORECASE):
                    path = os.path.relpath(path)

                options.in_fasta_map[name] = (weight, path)

    def overall_genome_sanity_check():
        """Checks if the given genome name of parameters 'genome', 'parent_genomes', and 'child_genome'
        exist in the HAL file"""
        keywords = ["genome", "parent_genomes", "child_genome"]

        for kw in keywords:
            if kw in options:
                genome = eval(f"options.{kw}")
                if genome not in existing_genomes:
                    raise RuntimeError(
                        f"Genome '{genome}' given by parameter '--{kw}' "
                        f"not found in the HAL file '{options.in_hal}'. "
                        f"Action '{options.action}' requires it."
                    )

    def replace_action_sanity_check():
        """Checks several things for replace action
        a) if the replacement genome occurs for leaves only
        b) if the given file mention one genome for replacement only
        c) if the genome asked for replacement matches with the one in the input file
        """
        # a) leaf-genome sanity check
        children_genomes = (
            cactus_call(
                check_output=True,
                parameters=[
                    "halStats",
                    "--children",
                    options.genome,
                    options.in_hal,
                ],
            )
            .strip()
            .split()
        )
        if len(children_genomes) > 0:
            raise RuntimeError(
                f"Node {options.genome} has children {children_genomes}. "
                f"Genome replacement works for leaf genomes only."
            )

        # b) one replecement at time
        if len(options.in_fasta_map) != 1:
            raise RuntimeError(
                f"One genome replacement genome at time. "
                f"File file {options.in_fasta} states {len(options.in_fasta_map)} "
                f"genomes for replacement."
            )

        # c) replacement sanity-check 1: if the genome asked for replacement exist in the given HAL file
        if options.genome not in existing_genomes:
            raise RuntimeError(
                f"Genome asked for replacement '{options.genome}' "
                f"does not exist in the given HAL File  "
                f"file {options.in_hal}"
            )

        # c) if the genome asked for replacement matches with the
        # name given in the input file
        if options.genome not in options.in_fasta_map:
            raise RuntimeError(
                f"Genome asked for replacement '{options.genome}' "
                f"does not match with the one mentioned in the input "
                f"file {options.in_fasta}: {list(options.in_fasta_map.keys())}"
            )

    def add_action_sanity_check():
        """Checkes if the genomes stated in the given file doesn't exist in the HAL file"""

        found = list()
        for genome in options.in_fasta_map.keys():
            if genome in existing_genomes:
                found.append(genome)

        if len(found) > 0:
            raise RuntimeError(
                f"HAL file already contains {found}, therefore "
                f"the action of {options.action}ing {found} to a {options.subcommand} "
                f"can't continue."
            )

    # get existing genomes in the hal file
    existing_genomes = (
        cactus_call(
            check_output=True,
            parameters=["halStats", "--genomes", options.in_hal],
        )
        .strip()
        .split()
    )

    # validate the given in_fasta_map text file
    in_fasta_map_sanity_check()

    # validate the given in_hal file and create a bkp if not requested
    in_hal_sanity_check()

    # validate --cactus-prepare-options
    cactus_prepare_options_sanity_check()

    # validate --outDir and --jobStore, and --outSeqFile
    work_dir_sanity_check()

    # sanity check for --genome --parent_genomes and --child_genome parameters
    overall_genome_sanity_check()

    # sanity check depending on the action
    if "add" == options.action:
        # sanity check for add action
        add_action_sanity_check()

    elif "replace" == options.action:
        # sanity check for replace action
        replace_action_sanity_check()
    else:
        raise RuntimeError(f"Unknown subcommand '{options.action}'")

    return options


def main():
    """Main cactus-update-prepare function"""

    # Same main parser as usual
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        usage="%(prog)s {add,replace} [-h] [Options]",
        add_help=False,
    )
    parser._positionals.title = "Update subcommand"

    # hack to add --help
    parser.add_argument("--help", "-h", action="store_true", dest="help", help=SUPPRESS)

    # Create parent subparser. Note `add_help=False` and creation via argparse
    parent_parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter, add_help=False
    )
    parent_parser.add_argument(
        "in_hal",
        metavar="inHal",
        help="The input HAL-format file containing the existing alignment",
    )
    parent_parser.add_argument(
        "in_fasta",
        metavar="inFasta",
        help="Tab-separated file. First column: genome name, second column: FASTA-format file path",
    )

    # taken from cactus-prepare
    parent_parser.add_argument(
        "--halOptions",
        dest="halOptions",
        type=str,
        default="--hdf5InMemory",
        help="options for every hal command",
    )
    parent_parser.add_argument(
        "--outDir",
        type=str,
        default="./steps",
        help="Directory where assemblies and cactus-call dependencies will be placed.",
        dest="outDir",
    )
    parent_parser.add_argument(
        "--outSeqFile",
        type=str,
        help="Path for annotated Seq file output [default: outDir/seqFile]",
        dest="outSeqFile",
    )
    parent_parser.add_argument(
        "--jobStore",
        type=str,
        default="./jobstore",
        dest="jobStore",
        help="base directory of jobStores to use in suggested commands",
    )
    parent_parser.add_argument(
        "--latest",
        dest="latest",
        action="store_true",
        help="Use the latest version of the docker container "
        "rather than pulling one matching this version of cactus",
    )
    parent_parser.add_argument(
        "--containerImage",
        dest="containerImage",
        default=None,
        help="Use the the specified pre-built container image "
        "rather than pulling one from quay.io",
    )
    parent_parser.add_argument(
        "--binariesMode",
        choices=["docker", "local", "singularity"],
        help="The way to run the Cactus binaries (at top level; use --cactusOpts to set it in nested calls)",
        default=None,
    )
    # Taken from https://github.com/DataBiosphere/toil/blob/559c313b8003b31b73cb0083fcab0ba30cf47a77/src/toil/common.py#L508
    # because `workDir` is now used by setupBinaries
    # at https://github.com/ComparativeGenomicsToolkit/cactus/blob/12a49f8b535060216d53134658709ea62efd125f/src/cactus/shared/common.py#L348-L413
    # since commit 6d94a879173924c8ae38c243435725071b79e9aa
    parent_parser.add_argument(
        "--workDir",
        dest="workDir",
        default=None,
        help="Absolute path to directory where temporary files generated during the Toil "
        "run should be placed. Standard output and error from batch system jobs "
        "(unless --noStdOutErr) will be placed in this directory. A cache directory "
        "may be placed in this directory. Temp files and folders will be placed in a "
        "directory toil-<workflowID> within workDir. The workflowID is generated by "
        "Toil and will be reported in the workflow logs. Default is determined by the "
        "variables (TMPDIR, TEMP, TMP) via mkdtemp. This directory needs to exist on "
        "all machines running jobs; if capturing standard output and error from batch "
        "system jobs is desired, it will generally need to be on a shared file system. "
        "When sharing a cache between containers on a host, this directory must be "
        "shared between the containers.",
    )

    # new for cactus-update-prepare
    parent_parser.add_argument(
        "--skip-backup",
        action="store_true",
        help="Skip the backup of the given HAL file",
        dest="skip_backup",
        default=False,
    )
    parent_parser.add_argument(
        "--skip-halValidate",
        help="Skip the validation of the given HAL file",
        action="store_true",
        dest="skip_halValidate",
        default=False,
    )
    parent_parser.add_argument(
        "--cactus-prepare-options",
        dest="cactus_prepare_options",
        type=str,
        default="--preprocessBatchSize 1 --cactusOptions '--realTimeLogging --logInfo --retryCount 0'",
        help="Options to bypass local configuration for cactus-prepare",
    )

    # add logging option
    add_logging_options(parent_parser)

    # overall subparser
    subparsers = parser.add_subparsers(
        help="Desired update for the given alignment", dest="action"
    )

    # adding "add" subcommand and its "node" and "branch" actions
    add_parser = subparsers.add_parser("add", help="Adding a new genome")
    add_subparser = add_parser.add_subparsers(
        help="Methods to add a new genome", dest="subcommand"
    )
    subcommand_options(add_subparser, parent_parser, "add", "node")
    subcommand_options(add_subparser, parent_parser, "add", "branch")

    # adding "replace" subcommand as an action
    subcommand_options(subparsers, parent_parser, "replace", None)

    options = parser.parse_args()

    # there is None subcommnd for action "replace"
    if options.action == "replace":
        options.subcommand = None

    if not options.action:
        parser.print_help()
        sys.exit(1)

    if options.action != "replace" and not options.subcommand:
        add_parser.print_help()
        sys.exit(1)

    # From cactus-prepare
    setupBinaries(options)
    set_logging_from_options(options)

    # update options and sanity checks before move on
    options = sanity_checks(options)

    # cactus input seqFile is defined automatically
    options.inSeqFile = os.path.join(options.outDir, "seq_file.in")

    # cactus output seqFile is defined automatically
    options.outSeqFile = os.path.join(options.outDir, "seq_file.out")

    # cactus-prepare parser wrapper
    cactus_alignment_update(options)


if __name__ == "__main__":
    main()
