#!/usr/bin/env python3

"""This sciprt will take the cactus-graphmap output and input and split it up into chromosomes.  This is done after a whole-genome graphmap as we need the whole-genome minigraph alignments to do the chromosome splitting in the first place.  For each chromosome, it will make a PAF, GFA and Seqfile (pointing to chrosome fastas)
"""
import os
from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import copy
import timeit

from operator import itemgetter

from cactus.progressive.seqFile import SeqFile
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.common import setupBinaries, importSingularityImage
from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.progressive.schedule import Schedule
from cactus.progressive.projectWrapper import ProjectWrapper
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.pipeline.cactus_workflow import CactusWorkflowArguments
from cactus.pipeline.cactus_workflow import addCactusWorkflowOptions
from cactus.pipeline.cactus_workflow import CactusTrimmingBlastPhase
from cactus.pipeline.cactus_workflow import prependUniqueIDs
from cactus.shared.common import makeURL, catFiles
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.common import cactus_call
from cactus.shared.common import getOptionalAttrib, findRequiredNode
from cactus.shared.common import unzip_gz
from toil.job import Job
from toil.common import Toil
from toil.lib.bioio import logger
from toil.lib.bioio import setLoggingFromOptions
from toil.realtimeLogger import RealtimeLogger
from toil.lib.threading import cpu_count

from sonLib.nxnewick import NXNewick
from sonLib.bioio import getTempDirectory

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    addCactusWorkflowOptions(parser)

    parser.add_argument("seqFile", help = "Seq file (gzipped fastas supported)")
    parser.add_argument("minigraphGFA", help = "Minigraph-compatible reference graph in GFA format (can be gzipped)")
    parser.add_argument("graphmapPAF", type=str, help = "Output pairwise alignment file in PAF format (can be gzipped)")
    parser.add_argument("--outDir", required=True, type=str, help = "Output directory")
    parser.add_argument("--refContigs", nargs="*", help = "Subset to these reference contigs (multiple allowed)", default=[])
    parser.add_argument("--refContigsFile", type=str, help = "Subset to (newline-separated) reference contigs in this file")
    parser.add_argument("--otherContig", type=str, help = "Lump all reference contigs unselected by above options into single one with this name")
    
    #Progressive Cactus Options
    parser.add_argument("--configFile", dest="configFile",
                        help="Specify cactus configuration file",
                        default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))
    parser.add_argument("--latest", dest="latest", action="store_true",
                        help="Use the latest version of the docker container "
                        "rather than pulling one matching this version of cactus")
    parser.add_argument("--containerImage", dest="containerImage", default=None,
                        help="Use the the specified pre-built containter image "
                        "rather than pulling one from quay.io")
    parser.add_argument("--binariesMode", choices=["docker", "local", "singularity"],
                        help="The way to run the Cactus binaries", default=None)

    options = parser.parse_args()

    setupBinaries(options)
    setLoggingFromOptions(options)
    enableDumpStack()

    # todo: would be very nice to support s3 here
    if options.outDir:
        if not os.path.isdir(options.outDir):
            os.makedirs(options.outDir)
        
    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    start_time = timeit.default_timer()
    runCactusGraphMapSplit(options)
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-graphmap-split has finished after {} seconds".format(run_time))

def runCactusGraphMapSplit(options):
    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            split_id_map = toil.restart()
        else:
            options.cactusDir = getTempDirectory()

            #load cactus config
            configNode = ET.parse(options.configFile).getroot()
            config = ConfigWrapper(configNode)
            config.substituteAllPredefinedConstantsWithLiterals()

            # load up the contigs if any
            ref_contigs = set(options.refContigs)
            # todo: use import?
            if options.refContigsFile:
                with open(options.refContigsFile, 'r') as rc_file:
                    for line in rc_file:
                        if len(line.strip()):
                            ref_contigs.add(line.strip().split()[0])

            if options.otherContig:
                assert options.otherContig not in ref_contigs

            # get the minigraph "virutal" assembly name
            graph_event = getOptionalAttrib(findRequiredNode(configNode, "graphmap"), "assemblyName", default="__MINIGRAPH_SEQUENCES__")

            # load the seqfile
            seqFile = SeqFile(options.seqFile)
            
            #import the graph
            gfa_id = toil.importFile(makeURL(options.minigraphGFA))

            #import the paf
            paf_id = toil.importFile(makeURL(options.graphmapPAF))

            #import the sequences (that we need to align for the given event, ie leaves and outgroups)
            seqIDMap = {}
            leaves = set([seqFile.tree.getName(node) for node in seqFile.tree.getLeaves()])
            for genome, seq in seqFile.pathMap.items():
                if genome in leaves:
                    if os.path.isdir(seq):
                        tmpSeq = getTempFile()
                        catFiles([os.path.join(seq, subSeq) for subSeq in os.listdir(seq)], tmpSeq)
                        seq = tmpSeq
                    seq = makeURL(seq)
                    logger.info("Importing {}".format(seq))
                    seqIDMap[genome] = (seq, toil.importFile(seq))

            # todo: better error -- its easy to make this mistake
            assert graph_event in seqIDMap

            # run the workflow
            split_id_map = toil.start(Job.wrapJobFn(graphmap_split_workflow, options, config, seqIDMap,
                                                    gfa_id, options.minigraphGFA,
                                                    paf_id, options.graphmapPAF, ref_contigs, options.otherContig))

        #export the split data
        export_split_data(toil, split_id_map, options.outDir)

def graphmap_split_workflow(job, options, config, seqIDMap, gfa_id, gfa_path, paf_id, paf_path, ref_contigs, other_contig):

    root_job = Job()
    job.addChild(root_job)

    # get the sizes before we overwrite below
    gfa_size = gfa_id.size
    paf_size = paf_id.size
    
    # use file extension to sniff out compressed input
    if gfa_path.endswith(".gz"):
        gfa_id = root_job.addChildJobFn(unzip_gz, gfa_path, gfa_id, disk=gfa_id.size * 10).rv()
        gfa_size *= 10
    if paf_path.endswith(".gz"):
        paf_id = root_job.addChildJobFn(unzip_gz, paf_path, paf_id, disk=paf_id.size * 10).rv()
        paf_size *= 10

    # cactus-ids use alphabetical ordering.  we need this as our paf file will have them
    cactus_id_map = {}
    for i, event in enumerate(sorted(set(list(seqIDMap.keys())))):
        cactus_id_map[event] = i
        
    # use rgfa-split to split the gfa and paf up by contig
    split_gfa_job = root_job.addFollowOnJobFn(split_gfa, config, gfa_id, paf_id, ref_contigs, cactus_id_map,
                                              other_contig,
                                              disk=(gfa_size + paf_size) * 5)

    # use the output of the above splitting to do the fasta splitting
    split_fas_job = split_gfa_job.addFollowOnJobFn(split_fas, seqIDMap, cactus_id_map, split_gfa_job.rv())

    # gather everythign up into a table
    gather_fas_job = split_fas_job.addFollowOnJobFn(gather_fas, seqIDMap, split_gfa_job.rv(), split_fas_job.rv())

    # return all the files
    return gather_fas_job.rv()

def split_gfa(job, config, gfa_id, paf_id, ref_contigs, cactus_id_map, other_contig):
    """ Use rgfa-split to divide a GFA and PAF into chromosomes.  The GFA must be in minigraph RGFA output using
    the desired reference. """

    work_dir = job.fileStore.getLocalTempDir()
    gfa_path = os.path.join(work_dir, "mg.gfa")
    paf_path = os.path.join(work_dir, "mg.paf")
    out_prefix = "split_"

    job.fileStore.readGlobalFile(gfa_id, gfa_path)
    job.fileStore.readGlobalFile(paf_id, paf_path)

    # get the minigraph "virutal" assembly name
    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="__MINIGRAPH_SEQUENCES__")
    # and look up its unique id prefix.  this will be needed to pick its contigs out of the list
    mg_id = cactus_id_map[graph_event]

    cmd = ['rgfa-split', '-i', 'id={}|'.format(mg_id), '-G',
           '-g', os.path.basename(gfa_path),
           '-p', os.path.basename(paf_path),
           '-b', out_prefix]    
    if other_contig:
        cmd += ['-o', other_contig]
        
    for contig in ref_contigs:
        cmd += ['-c', contig]

    cactus_call(parameters=cmd, work_dir=work_dir)

    output_id_map = {}
    for out_name in os.listdir(work_dir):
        file_name, ext = os.path.splitext(out_name)
        if file_name.startswith(out_prefix) and ext in [".gfa", ".paf", ".fa_contigs"]:
            name = file_name[len(out_prefix):]
            if name not in output_id_map:
                output_id_map[name] = {}
            output_id_map[name][ext[1:]] = job.fileStore.writeGlobalFile(os.path.join(work_dir, out_name))
            
    return output_id_map

def split_fas(job, seq_id_map, cactus_id_map, split_id_map):
    """ Use samtools to split a bunch of fasta files into reference contigs, using the output of rgfa-split as a guide"""

    root_job = Job()
    job.addChild(root_job)
    
    # map event name to dict of contgs.  ex fa_contigs["CHM13"]["chr13"] = file_id
    fa_contigs = {}
    # we do each fasta in parallel
    for event in seq_id_map.keys():
        fa_path, fa_id = seq_id_map[event]
        cactus_id = cactus_id_map[event]
        fa_contigs[event] = root_job.addChildJobFn(split_fa_into_contigs, event, fa_id, fa_path, cactus_id, split_id_map,
                                                   disk=fa_id.size * 3).rv()

    return fa_contigs

def split_fa_into_contigs(job, event, fa_id, fa_path, cactus_id, split_id_map):
    """ Use samtools turn on fasta into one for each contig. this relies on the informatino in .fa_contigs
    files made by rgfa-split """

    # download the fasta
    work_dir = job.fileStore.getLocalTempDir()
    fa_path = os.path.join(work_dir, os.path.basename(fa_path))
    is_gz = fa_path.endswith(".gz")
    job.fileStore.readGlobalFile(fa_id, fa_path, mutable=is_gz)
    if is_gz:
        # samtools can only work on bgzipped files.  so we uncompress here to be able to support gzipped too
        cactus_call(parameters=['gzip', '-fd', fa_path])
        fa_path = fa_path[:-3]

    unique_id = 'id={}|'.format(cactus_id)
                        
    contig_fa_dict = {}

    for ref_contig in split_id_map.keys():
        query_contig_list_id = split_id_map[ref_contig]['fa_contigs']
        list_path = os.path.join(work_dir, '{}.fa_contigs'.format(ref_contig))
        job.fileStore.readGlobalFile(query_contig_list_id, list_path)
        faidx_input_path = os.path.join(work_dir, '{}.fa_contigs.clean'.format(ref_contig))
        contig_count = 0
        with open(list_path, 'r') as list_file, open(faidx_input_path, 'w') as clean_file:
            for line in list_file:
                query_contig = line.strip()
                if query_contig.startswith(unique_id):
                    assert query_contig.startswith(unique_id)
                    query_contig = query_contig[len(unique_id):]
                    clean_file.write('{}\n'.format(query_contig))
                    contig_count += 1
        contig_fasta_path = os.path.join(work_dir, '{}_{}.fa'.format(event, ref_contig))
        if contig_count > 0:
            cmd = ['samtools', 'faidx', fa_path, '--region-file', faidx_input_path] 
            cactus_call(parameters=cmd, outfile=contig_fasta_path)
        else:
            # TODO: review how cases like this are handled
            with open(contig_fasta_path, 'w') as empty_file:
                empty_file.write("")
        contig_fa_dict[ref_contig] = job.fileStore.writeGlobalFile(contig_fasta_path)

    return contig_fa_dict

def gather_fas(job, seq_id_map, output_id_map, contig_fa_map):
    """ take the split_fas output which has everything sorted by event, and move into the ref-contig-based table
    from split_gfa.  return the updated table, which can then be exported into the chromosome projects """

    for ref_contig in output_id_map.keys():
        output_id_map[ref_contig]['fa'] = {}
        for event, fa_id in contig_fa_map.items():
            output_id_map[ref_contig]['fa'][event] = fa_id[ref_contig]

    return output_id_map

def export_split_data(toil, output_id_map, output_dir):
    """ download all the split data locally """

    chrom_file_map = {}
    
    for ref_contig in output_id_map.keys():
        ref_contig_path = os.path.join(output_dir, ref_contig)
        if not os.path.isdir(ref_contig_path):
            os.makedirs(ref_contig_path)

        # GFA: <output_dir>/<contig>/<contig>.gfa
        toil.exportFile(output_id_map[ref_contig]['gfa'], makeURL(os.path.join(ref_contig_path, '{}.gfa'.format(ref_contig))))

        # PAF: <output_dir>/<contig>/<contig>.paf
        toil.exportFile(output_id_map[ref_contig]['paf'], makeURL(os.path.join(ref_contig_path, '{}.paf'.format(ref_contig))))

        # Fasta: <output_dir>/<contig>/fasta/<event>_<contig>.fa ..
        seq_file_map = {}
        for event, ref_contig_fa_id in output_id_map[ref_contig]['fa'].items():
            fa_base = os.path.join(ref_contig_path, 'fasta')
            if not os.path.isdir(fa_base):
                os.makedirs(fa_base)
            fa_path = makeURL(os.path.join(fa_base, '{}_{}.fa'.format(event, ref_contig)))
            seq_file_map[event] = fa_path
            toil.exportFile(ref_contig_fa_id, fa_path)

        # Seqfile: <output_dir>/<contig>/<contig>.seqfile
        seq_file_path = os.path.join(ref_contig_path, '{}.seqfile'.format(ref_contig))
        with open(seq_file_path, 'w') as seq_file:
            for event, fa_path in seq_file_map.items():
                seq_file.write('{}\t{}\n'.format(event, fa_path))

        # Top-level seqfile
        chrom_file_map[ref_contig] = makeURL(seq_file_path)
        
        
    with open(os.path.join(output_dir, 'chromfile.txt'), 'w') as chromfile:
        for ref_contig, seqfile_path in chrom_file_map.items():
            chromfile.write('{}\t{}\n'.format(ref_contig, makeURL(seqfile_path)))
    
if __name__ == "__main__":
    main()
