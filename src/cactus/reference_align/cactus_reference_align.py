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

import subprocess
import os
from argparse import ArgumentParser
import collections as col

from cactus.reference_align import paf_to_lastz
from cactus.reference_align import fasta_preprocessing

from cactus.shared.common import makeURL
from cactus.pipeline.cactus_workflow import prependUniqueIDs

## utilitary fxns:

def apply_dipcall_flt(job, infile, min_var_len=50000, min_mapq=5):
    """Filters out all mappings below min_var_len and min_mapq.
    NOTE: Assumes all secondary mappings are already removed. 
    Also: lastz cigars need to have <score> field filled with mapq, not raw score.
    Args:
        infile (lastz cigar format): [description]
        outfile ([type]): [description]
        min_var_len ([type]): [description]
        min_mapq ([type]): [description]
    """
    # debug = 0
    with open(job.fileStore.readGlobalFile(infile)) as inf:
        filtered = job.fileStore.getLocalTempFile()
        with open(filtered, "w+") as outf:
            for line in inf:
                parsed = line.split()
                # print(parsed, "\t", variation_length(parsed[10:]), "\t", parsed[9])
                # debug += 1
                # if debug == 2:
                #     break
                # if (parsed[3] - parsed[2]) >= min_var_len: This measures the length of the SV from the query point of view. Dipcall also includes deletions in measuring variation length.
                if variation_length(parsed[10:]) >= min_var_len:
                    if int(parsed[9]) >= min_mapq:
                        outf.write(line)

    return job.fileStore.writeGlobalFile(filtered)

def variation_length(lastz_cig_list):
    # 0, 2, 4,... are the type of variation.
    # 1, 3, 5,... are the length of the variation of that type.
    # there should be an even number of entries in the cig list: pairs of type, value.
    var_len = 0
    for i in range(0, len(lastz_cig_list), 2):
        print(i, lastz_cig_list[i], lastz_cig_list[i+1])
        if lastz_cig_list[i] in "IDM":
            var_len += int(lastz_cig_list[i+1])
    return var_len

def unpack_promise(job, iterable, i):
    """
    passed an iterable and a location i, returns ith item.
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

def get_asms_from_seqfile(seqfile):
    asm_files = col.OrderedDict()

    with open(seqfile) as inf:
        #skip the Newick tree on the first line:
        next(inf)
        
        for line in inf:
            parsed = line.split()
            if len(parsed) >= 2:
                #note: fastas can be in a directory containing single consecutive spaces. Two spaces in a row breaks my file parsing system. 
                if parsed[0][0] == "*":
                    asm_files[parsed[0][1:]] = " ".join(parsed[1:])
                else:
                    asm_files[parsed[0]] = " ".join(parsed[1:])

    return asm_files

def import_asms_non_blast_output(options, workflow):
    """Import asms; deduplicating contig ids if not --all_unique_ids

    Args:
        seqfile ([type]): [description]
        workflow ([type]): [description]

    Returns:
        [type]: [description]
    """
    # asms is dictionary of all asms (not counting reference) with key: asm_name, value: imported global toil file.
    asms = get_asms_from_seqfile(options.seqFile)

    if not options.all_unique_ids:
        # deduplicate contig id names, if user hasn't guaranteed unique contig ids.
        # new_fastas is the location of the asms with unique ids.
        if options.overwrite_assemblies:
            # overwrite the original assemblies. Note that the reference is never overwritten, as it is never altered. 
            # (Code assumes that the reference is internally free of duplicate ids, and just ensures other asms don't use reference ids.)
            asms = fasta_preprocessing.rename_duplicate_contig_ids(asms, options.refID, asms)
        else:
            # don't overwrite the original assemblies.
            # first, determine the new asm save locations.
            if not os.path.isdir(options.assembly_save_dir):
                os.mkdir(options.assembly_save_dir)

            new_asms = dict()
            for asm_id, asm in asms.items():
                if asm_id != options.refID:
                    new_asms[asm_id] = options.assembly_save_dir + asm.split("/")[-1]
                else:
                    # reference file never needs deduplication of contig ids, since ref contig ids are counted before all other asms.
                    new_asms[asm_id] = asm
                    
            asms = fasta_preprocessing.rename_duplicate_contig_ids(asms, options.refID, new_asms)

    # Import asms.
    for asm_id, asm in asms.items():
        asms[asm_id] = workflow.importFile('file://' + os.path.abspath(asm))

    return asms

def import_asms(options, workflow):
    """Import asms; deduplicating contig ids if not --all_unique_ids

    Args:
        options ([type]):
        workflow ([type]): [description]

    Returns:
        [type]: [description]
    """
    # asms is orderedDict of all asms (not counting reference) with key: asm_name, value: imported global toil file.

    asms = get_asms_from_seqfile(options.seqFile)

    # Prepend unique IDs.
    uniqueFas = prependUniqueIDs(asms.values(), options.assembly_save_dir)

    i = 0
    for asm_id in asms.keys():
        asms[asm_id] = uniqueFas[i]
        i += 1
        
    # Import asms.
    for asm_id, asm in asms.items():
        asms[asm_id] = workflow.importFile(makeURL(asm))

    return asms

def empty(job):
    """
    An empty job, for easier toil job organization.
    """
    return

## mapping fxns:

def map_all_to_ref(job, assembly_files, reference, debug_export, dipcall_filter=False):
    """
    Primarily for use with option_all_to_ref_only. Otherwise, use map_all_to_ref_and_get_poor_mappings.
    """
    lead_job = job.addChildJobFn(empty)

    # map all assemblies to the reference. Don't map reference to reference, though.
    ref_mappings = dict()
    for assembly, assembly_file in assembly_files.items():
        if assembly != reference:
            ref_mappings[assembly] = lead_job.addChildJobFn(map_a_to_b, assembly_file, assembly_files[reference], dipcall_filter).rv()
    
    consolidate_job = lead_job.addFollowOnJobFn(consolidate_mappings, ref_mappings)
    paf_mappings = consolidate_job.rv()

    conversion_job = consolidate_job.addFollowOnJobFn(paf_to_lastz.paf_to_lastz, paf_mappings)
    lastz_mappings = conversion_job.rv()

    primary_mappings = conversion_job.addChildJobFn(unpack_promise, lastz_mappings, 0).rv()
    secondary_mappings = conversion_job.addChildJobFn(unpack_promise, lastz_mappings, 1).rv()

    if debug_export:
        return (primary_mappings, secondary_mappings, ref_mappings, paf_mappings )
    else:
        return (primary_mappings, secondary_mappings)

def map_a_to_b(job, a, b, dipcall_filter):
    """Maps fasta a to fasta b.

    Args:
        a (global file): fasta file a. In map_all_to_ref, a is an assembly fasta.
        b (global file): fasta file b. In map_all_to_ref, b is the reference.

    Returns:
        [type]: [description]
    """
    
    # map_to_ref_paf = job.fileStore.writeGlobalFile(job.fileStore.getLocalTempFile())
    tmp = job.fileStore.getLocalTempFile()
    map_to_ref_paf = job.fileStore.writeGlobalFile(tmp)

    if dipcall_filter:
        subprocess.call(["minimap2", "-cx", "asm5", "-r2k", "-o", job.fileStore.readGlobalFile(map_to_ref_paf),
                        job.fileStore.readGlobalFile(b), job.fileStore.readGlobalFile(a)])
    else:
        subprocess.call(["minimap2", "-cx", "asm5", "-o", job.fileStore.readGlobalFile(map_to_ref_paf),
                        job.fileStore.readGlobalFile(b), job.fileStore.readGlobalFile(a)])
    

    return map_to_ref_paf


## main fxn and interface:

def get_options():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    
    # options for basic input/output
    parser.add_argument('seqFile', type=str,
                        help='A file containing all the information specified by cactus in construction. This aligner ignores the newick tree.')
    parser.add_argument('refID', type=str, 
                        help='Specifies which asm in seqFile should be treated as the reference.')
    parser.add_argument("outputFile", type=str, help = "Output pairwise alignment file")
    parser.add_argument('--dipcall_filter', action='store_true', 
                        help="Applies filters & minimap2 arguments used in dipcall. Only affects the primary mappings file. Secondary mappings aren't used in dipcall.")
    # parser.add_argument('--secondary', default="secondary.cigar", type=str, 
    #                     help='Filename for where to write lastz cigar output for secondary mappings.')
                        
    ## options for importing assemblies:
    # following arguments are only useful under --non_blast_output
    parser.add_argument('--non_blast_output', action='store_true', 
                    help="Instead of using cactus-blast-style prepended ids, use an alternative import method that only alters contig ids if absolutely necessary.")
    parser.add_argument('--all_unique_ids', action='store_true', 
                        help="Prevents the program from touching the assembly files; the user promises that they don't contain any duplicate contig ids. In reality, there should never be contig renamings if there are no duplicate fasta ids.")
    parser.add_argument('--overwrite_assemblies', action='store_true', 
                        help="When cleaning the assembly files to make sure there are no duplicate contig ids, overwrite the assembly files. Copy them to a neigboring folder with the affix '_edited_for_duplicate_contig_ids' instead.")

    # Useful in normal asms import
    parser.add_argument('--assembly_save_dir', type=str, default='./unique_id_assemblies/',
                        help='While deduplicating contig ids in the input fastas, save the assemblies in this directory. Ignored when used in conjunction with --overwrite_assemblies.')
                        
    # for debugging:
    parser.add_argument('--debug_export', action='store_true',
                        help='Export several other files for debugging inspection.')


    options = parser.parse_args()
    return options

def main():
    options = get_options()

    with Toil(options) as workflow:
        ## Preprocessing:
        # Import asms; by default, prepends unique IDs in the technique used in cactus-blast.
        if options.non_blast_output:
            asms = import_asms_non_blast_output(options, workflow)
        else:
            asms = import_asms(options, workflow)
            
        ## Perform alignments:
        if not workflow.options.restart:
            alignments = workflow.start(Job.wrapJobFn(map_all_to_ref, asms, options.refID, options.debug_export, options.dipcall_filter))

        else:
            alignments = workflow.restart()

        if options.debug_export:
            print(alignments)
            # Then return value is: (primary_mappings, secondary_mappings, ref_mappings, paf_mappings )
            for asm, mapping_file in alignments[2].items():
                workflow.exportFile(mapping_file, 'file://' + os.path.abspath("debug_" + asm + "_mapping_to_ref.txt"))
                break
            workflow.exportFile(alignments[3], 'file://' + os.path.abspath("debug_paf_mappings.txt"))

        ## Save alignments:
        if not options.dipcall_filter:
            workflow.exportFile(alignments[0], makeURL(options.outputFile))
            workflow.exportFile(alignments[1], makeURL(options.outputFile + ".secondary"))
        else:
            dipcall_filtered = workflow.start(Job.wrapJobFn(apply_dipcall_flt, alignments[0]))
            workflow.exportFile(dipcall_filtered, makeURL(options.outputFile))
            workflow.exportFile(alignments[1], makeURL(options.outputFile + ".secondary"))

        
        # workflow.exportFile(alignments[0], 'file://' + os.path.abspath(options.primary))
        # workflow.exportFile(alignments[1], 'file://' + os.path.abspath(options.secondary))

if __name__ == "__main__":
    main()