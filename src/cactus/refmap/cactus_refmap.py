#!/usr/bin/env python3

"""Feature Wishlist
Top priority:
Internal contig name system (see prependUniqueIDs in cactus/src/cactus/pipeline/cactus_workflow.py, line 465)
Called in line 529:
    renamedInputSeqDir = fileStore.getLocalTempDir()
    uniqueFas = prependUniqueIDs(sequences, renamedInputSeqDir)
    uniqueFaIDs = [fileStore.writeGlobalFile(seq, cleanup=True) for seq in uniqueFas]

    I'm currently calling prependUniqueIDs outside of the workflow/fileStore. I should
        check to make sure that Cactus also ultimately exports the uniqueID files,
        like I'm doing. Otherwise, I should delete the prepended IDs before finishing 
        the pipeline.



Implement options.pathOverrides
Implement "Progressive Cactus Options" (see cactus_blast ArgumentParser)
    This includes s3 compatibility, I think?
logger.info timer (see line 91 of cactus_blast)

    """


from toil.common import Toil
from toil.job import Job
from toil.statsAndLogging import logger

import os
from argparse import ArgumentParser
import collections as col
import xml.etree.ElementTree as ET

from cactus.refmap import paf_to_lastz
from cactus.refmap import fasta_preprocessing
from cactus.refmap import apply_dipcall_bed_filter

from cactus.shared.common import setupBinaries, importSingularityImage
from cactus.shared.common import makeURL
from cactus.shared.common import cactus_call
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import cactusRootPath
from cactus.progressive.progressive_decomposition import compute_outgroups, parse_seqfile, get_subtree, get_spanning_subtree, get_event_set
from cactus.preprocessor.checkUniqueHeaders import sanitize_fasta_headers

## utilitary fxns:

def variation_length(lastz_cig):
    """Determines how long the mapping is (sum of insertion/deletion/match), based on the 
    lastz cig.
    Args:
        lastz_cig (string): a lastz cigar, e.g. "M50D3I5M30"

    Returns:
        int: sum of I+D+M in cigar.
    """
    # Parsing cigars:
    # indices 0, 2, 4,... are the type of variation.
    # indices 1, 3, 5,... are the length of the variation of that type.
    # there should be an even number of entries in the cig list: pairs of type, value.
    var_len = 0
    for i in range(0, len(lastz_cig), 2):
        print(i, lastz_cig[i], lastz_cig[i+1])
        if lastz_cig[i] in "IDM":
            var_len += int(lastz_cig[i+1])
    return var_len

def unpack_promise(job, iterable, i):
    """
    passed an iterable and a location i, returns ith item. Useful for accessing the 
    contents of a toil promise.
    """
    return iterable[i]

def consolidate_mappings(job, mapping_files):
    """
    Warning: discards headers of all mapping files.
    Given a list of mapping files, consolidates the contents (not counting headers) into a
    single file.
    """
    consolidated_mappings = job.fileStore.getLocalTempFile()
    with open(consolidated_mappings,"w") as outfile:
        for mapping_file in mapping_files.values():
            with open(job.fileStore.readGlobalFile(mapping_file)) as inf:
                for line in inf:
                    if not line.startswith("@"):
                        outfile.write(line)
    return job.fileStore.writeGlobalFile(consolidated_mappings)

def empty(job):
    """
    An empty job, for easier toil job organization.
    """
    return

## dipcall filter functions

def apply_dipcall_vcf_filter(job, infile, min_var_len=50000, min_mapq=5):
    """Filters out all mappings below min_var_len and min_mapq from a lastz file.
    NOTE: Assumes all secondary mappings are already removed. 
    Also: lastz cigars need to have <score> field filled with mapq, not raw score.
    Args:
        infile (lastz cigar format): [description]
        outfile ([type]): [description]
        min_var_len ([type]): [description]
        min_mapq ([type]): [description]
    """
    with open(job.fileStore.readGlobalFile(infile)) as inf:
        filtered = job.fileStore.getLocalTempFile()
        with open(filtered, "w+") as outf:
            for line in inf:
                parsed = line.split()
                if variation_length(parsed[10:]) >= min_var_len:
                    if int(parsed[9]) >= min_mapq:
                        outf.write(line)

    return job.fileStore.writeGlobalFile(filtered)

def filter_out_secondaries_from_paf(job, paf):
    """
    Removes all secondary mappings from paf file.
    """
    primary_paf = job.fileStore.getLocalTempFile()
    with open(primary_paf, "w") as outf:
        with open(job.fileStore.readGlobalFile(paf)) as inf:
            for line in inf:
                parsed = line.split()
                for i in parsed[11:]:
                    if i[:6] == "tp:A:P" or i[:6] == "tp:A:I":
                        outf.write(line)
                        break
    return job.fileStore.writeGlobalFile(primary_paf)

## mapping fxns:

def run_cactus_reference_align(job, assembly_files, reference, debug_export=False, dipcall_bed_filter=False, dipcall_vcf_filter=False):
    """
    Preprocesses assemblies, then runs mappings.
    """
    sanitize_job = job.addChildJobFn(sanitize_fasta_headers, assembly_files)
    mappings = sanitize_job.addFollowOnJobFn(map_all_to_ref, sanitize_job.rv(), reference, debug_export, dipcall_bed_filter, dipcall_vcf_filter).rv()
    return mappings

def map_all_to_ref(job, assembly_files, reference, debug_export=False, dipcall_bed_filter=False, dipcall_vcf_filter=False):
    """The meat of cactus-reference-align. Performs all mappings; applies dipcall 
    filter_bed_filter if necessary, then converts to lastz cigars.

    Args:
        assembly_files (orderedDict): key: asm_name; value: asm_file
        reference (string): asm_name of the reference.
        debug_export (bool): Export some intermediate files of the workflow, for debugging. Defaults to False.
        dipcall_bed_filter (bool): Apply the dipcall bed filter to the mappings. This will:
                                    Guarantee that there will be no overlapping mappings.
                                    * include mappings >=min_var_len=50kb in size.
                                    * BUT: will exclude regions of these mappings which overlap any other mappings >=min_size_mapping=10kb in size. (This includes other >=50kb mappings).
                                    * all mappings considered for inclusion or overlap must have >= 5 mapQ.
                                        Defaults to False.
        dipcall_vcf_filter (bool): Applies the preliminary requirements for the less-stringent vcf-filter. Ultimately, vcf-filter:
                                    * removes all secondary mappings
                                    * Filters out all mappings below min_var_len=50k and min_mapq=5 from a lastz file
                                     Defaults to False.
    """
    lead_job = job.addChildJobFn(empty)

    # map all assemblies to the reference. Don't map reference to reference, though.
    ref_mappings = dict()
    secondary_mappings = dict()
    primary_mappings = dict()
    for assembly, assembly_file in assembly_files.items():
        if assembly != reference:
            # map to a to b
            print("about to run map a to b. a:", assembly, job.fileStore.readGlobalFile(assembly_file), "b (ref):", reference, job.fileStore.readGlobalFile(assembly_files[reference]))
            map_job = lead_job.addChildJobFn(map_a_to_b, assembly_file, assembly_files[reference], (dipcall_bed_filter or dipcall_vcf_filter))
            ref_mappings[assembly] = map_job.rv()

            if dipcall_bed_filter:
                secondaries_filter_job = map_job.addFollowOnJobFn(filter_out_secondaries_from_paf, ref_mappings[assembly])
                primary_paf = secondaries_filter_job.rv()

                dipcall_bed_filter_job = secondaries_filter_job.addFollowOnJobFn(apply_dipcall_bed_filter.apply_dipcall_bed_filter, primary_paf)
                bed_filtered_primary_mappings = dipcall_bed_filter_job.rv()
                paf_mappings = bed_filtered_primary_mappings
            else:
                # convert mapping to lastz (and filter into primary and secondary mappings)
                paf_mappings = ref_mappings[assembly]

            # extract the primary and secondary mappings.
            primary_mappings[assembly] = paf_mappings
            secondary_mappings[assembly] = None

    # consolidate the primary mappings into a single file; same for secondary mappings.
    all_primary = lead_job.addFollowOnJobFn(consolidate_mappings, primary_mappings).rv()
    return all_primary

def map_a_to_b(job, a, b, dipcall_filter):
    """Maps fasta a to fasta b.

    Args:
        a (global file): fasta file a. In map_all_to_ref, a is an assembly fasta.
        b (global file): fasta file b. In map_all_to_ref, b is the reference.

    Returns:
        [type]: [description]
    """
    
    print("in map a to b. a:", a, "b:", b)
    # map_to_ref_paf = job.fileStore.writeGlobalFile(job.fileStore.getLocalTempFile())
    tmp = job.fileStore.getLocalTempFile()
    map_to_ref_paf = job.fileStore.writeGlobalFile(tmp)

    if dipcall_filter:
        # note: in dipcall, they include argument "--paf-no-hit". 
        # I don't see why they would include these "mappings", only to be filtered out 
        # later. I have not included the argument.
        cactus_call(parameters=["minimap2", "-c", "-xasm5", "--cs", "-r2k", "-o", job.fileStore.readGlobalFile(map_to_ref_paf),
                        job.fileStore.readGlobalFile(b), job.fileStore.readGlobalFile(a)])
    else:
        cactus_call(parameters=["minimap2", "-cx", "asm5", "-o", job.fileStore.readGlobalFile(map_to_ref_paf),
                        job.fileStore.readGlobalFile(b), job.fileStore.readGlobalFile(a)])
    

    return map_to_ref_paf


## main fxn and interface:

def get_options():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    # addCactusWorkflowOptions(parser)
    
    # ### For quick debugging of apply_dipcall_bed_filter:
    # parser.add_argument('paf', type=str,
    #                     help='For quick debugging of apply_dipcall_bed_filter.')

    
    # options for basic input/output
    parser.add_argument('seqFile', type=str,
                        help='A file containing all the information specified by cactus in construction. This aligner ignores the newick tree.')
    parser.add_argument('reference', type=str, 
                        help='Specifies which asm in seqFile should be treated as the reference.')
    parser.add_argument("outputFile", type=str, help = "Output pairwise alignment file")
    parser.add_argument("--pathOverrides", nargs="*", help="paths (multiple allowd) to override from seqFile")
    parser.add_argument("--pathOverrideNames", nargs="*", help="names (must be same number as --paths) of path overrides")

    # dipcall-like filters
    parser.add_argument('--dipcall_bed_filter', action='store_true', 
                        help="Applies filters & minimap2 arguments used to make the bedfile in dipcall. Only affects the primary mappings file. Secondary mappings aren't used in dipcall.")
    parser.add_argument('--dipcall_vcf_filter', action='store_true', 
                        help="Applies filters & minimap2 arguments used to make the vcf in dipcall. Only affects the primary mappings file. Secondary mappings aren't used in dipcall.")

    # Progressive Cactus Options:
    parser.add_argument("--configFile", dest="configFile",
                        help="Specify cactus configuration file",
                        default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))
    parser.add_argument("--latest", dest="latest", action="store_true",
                        help="Use the latest version of the docker container "
                        "rather than pulling one matching this version of cactus")
    parser.add_argument("--binariesMode", choices=["docker", "local", "singularity"],
                        help="The way to run the Cactus binaries", default=None)      
    parser.add_argument("--containerImage", dest="containerImage", default=None,
                        help="Use the the specified pre-built containter image "
                        "rather than pulling one from quay.io")

    ## options for importing assemblies:
    # following arguments are only useful under --non_blast_output
    # parser.add_argument('--non_blast_output', action='store_true', 
    #                 help="Instead of using cactus-blast-style prepended ids, use an alternative import method that only alters contig ids if absolutely necessary.")
    # parser.add_argument('--all_unique_ids', action='store_true', 
    #                     help="Only take effect when called with --non_blast_output. Prevents the program from touching the assembly files; the user promises that they don't contain any duplicate contig ids. In reality, there should never be contig renamings if there are no duplicate fasta ids.")
    # parser.add_argument('--overwrite_assemblies', action='store_true', 
    #                     help="When cleaning the assembly files to make sure there are no duplicate contig ids, overwrite the assembly files. Copy them to a neigboring folder with the affix '_edited_for_duplicate_contig_ids' instead.")

    # # Useful in normal asms import
    # parser.add_argument('--assembly_save_dir', type=str, default='./unique_id_assemblies/',
    #                     help='While deduplicating contig ids in the input fastas, save the assemblies in this directory. Ignored when used in conjunction with --overwrite_assemblies.')
                        
    # for debugging:
    parser.add_argument('--debug_export', action='store_true',
                        help='Export several other files for debugging inspection.')
    parser.add_argument('--debug_export_dir', type=str, default='./debug_export_dir/',
                        help='Location of the exported debug files.')
    options = parser.parse_args()
    return options

def main():
    options = get_options()

    with Toil(options) as toil:
        setupBinaries(options)

        importSingularityImage(options)

        # load up the seqfile and figure out the outgroups and schedule
        config_node = ET.parse(options.configFile).getroot()
        config_wrapper = ConfigWrapper(config_node)
        config_wrapper.substituteAllPredefinedConstantsWithLiterals(options)
        mc_tree, input_seq_map, og_candidates = parse_seqfile(options.seqFile, config_wrapper)
        og_map = compute_outgroups(mc_tree, config_wrapper, set(og_candidates))
        event_set = get_event_set(mc_tree, config_wrapper, og_map, mc_tree.getRootName())

        # apply path overrides.  this was necessary for wdl which doesn't take kindly to
        # text files of local paths (ie seqfile).  one way to fix would be to add support
        # for s3 paths and force wdl to use it.  a better way would be a more fundamental
        # interface shift away from files of paths throughout all of cactus
        if options.pathOverrides:
            for name, override in zip(options.pathOverrideNames, options.pathOverrides):
                input_seq_map[name] = override

        # check --reference input
        if options.reference:
            leaves = [mc_tree.getName(leaf) for leaf in mc_tree.getLeaves()]
            if options.reference not in leaves:
                raise RuntimeError("Genome specified with --reference, {}, not found in tree leaves".format(options.reference))
                
        #import the sequences
        input_seq_id_map = {}
        for (genome, seq) in input_seq_map.items():
            if genome in event_set:
                if os.path.isdir(seq):
                    tmpSeq = getTempFile()
                    catFiles([os.path.join(seq, subSeq) for subSeq in os.listdir(seq)], tmpSeq)
                    seq = tmpSeq
                seq = makeURL(seq)
                logger.info("Importing {}".format(seq))
                input_seq_id_map[genome] = toil.importFile(seq)
        
        ## Perform alignments:
        if not toil.options.restart:
            alignments = toil.start(Job.wrapJobFn(run_cactus_reference_align, input_seq_id_map, options.reference, options.debug_export, options.dipcall_bed_filter, options.dipcall_vcf_filter))

        else:
            alignments = toil.restart()

        ## Save alignments:
        if options.dipcall_vcf_filter: # this is substantially less restrictive than the dipcall_bed_filter. 
            dipcall_filtered = toil.start(Job.wrapJobFn(apply_dipcall_vcf_filter, alignments))
            toil.exportFile(dipcall_filtered, makeURL(options.outputFile))
        else:
            toil.exportFile(alignments, makeURL(options.outputFile))

if __name__ == "__main__":
    main()
