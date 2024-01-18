#!/usr/bin/env python3

#Released under the MIT license, see LICENSE.txt

"""Set up a given seqfile to include proprocessed/ancestral sequences, as well as to

"""
import os, sys
from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import copy
import hashlib
from sonLib.bioio import getTempDirectory
from datetime import datetime
import subprocess
import timeit
import shutil
import shlex
import math

from operator import itemgetter

from cactus.progressive.seqFile import SeqFile
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.progressive.progressive_decomposition import compute_outgroups, get_subtree, check_branch_lengths
from cactus.shared.common import cactusRootPath
from cactus.shared.common import enableDumpStack, setupBinaries
from cactus.shared.common import getDockerImage
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.version import cactus_commit
from cactus.shared.common import findRequiredNode
from cactus.shared.common import makeURL, cactus_call, RoundedJob
from cactus.shared.common import write_s3, has_s3, get_aws_region

from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from cactus.shared.common import cactus_cpu_count
from toil.realtimeLogger import RealtimeLogger
from toil.lib.conversions import human2bytes, bytes2human
from toil.lib.accelerators import count_nvidia_gpus

def main_toil():
    return main(toil_mode=True)

def main(toil_mode=False):
    parser = ArgumentParser()
    if toil_mode:
        Job.Runner.addToilOptions(parser)
        parser.add_argument("--latest", dest="latest", action="store_true",
                            help="Use the latest version of the docker container "
                            "rather than pulling one matching this version of cactus")
        parser.add_argument("--containerImage", dest="containerImage", default=None,
                            help="Use the the specified pre-built containter image "
                            "rather than pulling one from quay.io")
        parser.add_argument("--binariesMode", choices=["docker", "local", "singularity"],
                            help="The way to run the Cactus binaries (at top level; use --cactusOpts to set it in nested calls)",
                            default=None)
    parser.add_argument("seqFile", help = "Seq file")
    parser.add_argument("--outDir", help='Directory where the processed leaf sequence and ancestral sequences will be placed.'
                        ' Required when not using --wdl')
    parser.add_argument("--outSeqFile", help="Path for annotated Seq file output [default: outDir/seqFile]")
    parser.add_argument("--outHal", help="Output HAL file [default: outDir/out.hal]", required=toil_mode)
    parser.add_argument("--chromInfo",
                        help="Two-column file mapping genome (col 1) to comma-separated list of sex chromosomes. This information "
                        "will be used to guide outgroup selection so that, where possible, all chromosomes are present in"
                        " at least one outgroup.")    
            
    if not toil_mode:
        parser.add_argument("--wdl", action="store_true", help="output wdl workflow instead of list of commands")
        parser.add_argument("--noLocalInputs", action="store_true", help="dont embed local input paths in WDL script (as they will need"
                            " to be respecified when running on Terra")
        parser.add_argument("--jobStore", type=str, default="./jobstore", help="base directory of jobStores to use in suggested commands")
        parser.add_argument("--seqFileOnly", action="store_true", help="Only create output SeqFile (with no ancestors); do not make plan")
    parser.add_argument("--configFile", default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))
    parser.add_argument("--preprocessBatchSize", type=int, default=10, help="size (number of genomes) of preprocessing jobs")
    parser.add_argument("--halAppendBatchSize", type=int, default=100, help="size (number of genomes) of halAppendSubtree jobs (WDL-only)")
    parser.add_argument("--halOptions", type=str, default="--hdf5InMemory", help="options for every hal command")
    parser.add_argument("--cactusOptions", type=str, default="", help="options for every cactus command")
    parser.add_argument("--preprocessOnly", action="store_true", help="only decompose into preprocessor and cactus jobs")
    parser.add_argument("--dockerImage", type=str, help="docker image to use as wdl runtime")
    if not toil_mode:
        # no good reason this isn't supported in toil_mode, just time to implement
        parser.add_argument("--includeRoot", action='store_true', help="Root node's sequence is set from input: --includeRoot will be specified for its align and blast jobs")
    
    parser.add_argument("--gpu", nargs='?', const='all', default=None, help="toggle on GPU-enabled lastz, and specify number of GPUs (all available if no value provided)")
    parser.add_argument("--gpuType", default="nvidia-tesla-v100", help="GPU type (to set in WDL runtime parameters, use only with --wdl)")
    parser.add_argument("--nvidiaDriver", default="470.82.01", help="Nvidia driver version")
    parser.add_argument("--gpuZone", default="us-central1-a", help="zone used for gpu task")
    parser.add_argument("--zone", default="us-west2-a", help="zone used for all but gpu tasks")

    if not toil_mode:
        parser.add_argument("--defaultCores", type=int, help="Number of cores for each job unless otherwise specified")
    parser.add_argument("--preprocessCores", type=int, help="Number of cores for each cactus-preprocess job")
    parser.add_argument("--blastCores", type=int, help="Number of cores for each cactus-blast job")
    parser.add_argument("--alignCores", type=int, help="Number of cores for each cactus-align job")

    if not toil_mode:
        parser.add_argument("--defaultMemory", type=human2bytesN, help="Memory for each job unless otherwise specified. "
                            "Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes)")
    parser.add_argument("--preprocessMemory", type=human2bytesN, help="Memory for each cactus-preprocess job. "
                        "Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes)")
    parser.add_argument("--blastMemory", type=human2bytesN, help="Memory for each cactus-blast job. "
                        "Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes)")                                
    parser.add_argument("--alignMemory", type=human2bytesN, help="Memory for each cactus-align job. "
                        "Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes)")

    if not toil_mode:
        parser.add_argument("--defaultDisk", type=human2bytesN, help="Disk for each job unless otherwise specified. "
                            "Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes)")
    parser.add_argument("--preprocessDisk", type=human2bytesN, help="Disk for each cactus-preprocess job. "
                        "Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes)")
    parser.add_argument("--blastDisk", type=human2bytesN, help="Disk for each cactus-blast job. "
                        "Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes)")
    parser.add_argument("--alignDisk", type=human2bytesN, help="Disk for each cactus-align job. "
                        "Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes)")
    parser.add_argument("--halAppendDisk", type=human2bytesN, help="Disk for each halAppendSubtree job. "
                        "Standard suffixes like K, Ki, M, Mi, G or Gi are supported (default=bytes)")

    parser.add_argument("--preprocessPreemptible", type=int, help="Preemptible attempt count for each cactus-preprocess job [default=1]", default=1)
    parser.add_argument("--blastPreemptible", type=int, help="Preemptible attempt count for each cactus-blast job [default=1]", default=1)
    parser.add_argument("--alignPreemptible", type=int, help="Preemptible attempt count for each cactus-align job [default=1]", default=1)
    parser.add_argument("--halAppendPreemptible", type=int, help="Preemptible attempt count for each halAppendSubtree job [default=1]", default=1)

    parser.add_argument("--preprocessRetries", type=int, help="Retry attempts for each cactus-preprocess job [default=1]", default=1)
    parser.add_argument("--blastRetries", type=int, help="Retry attempts for each cactus-blast job [default=1]", default=1)
    parser.add_argument("--alignRetries", type=int, help="Retry attempts for each cactus-align job [default=1]", default=0)
    parser.add_argument("--halAppendRetries", type=int, help="Retry attempts for each halAppendSubtree job [default=1]", default=0)

    options = parser.parse_args()
    #todo support root option
    options.root = None

    if toil_mode:
        options.wdl = False
        options.noLocalInputs = False
        options.includeRoot = False
        options.outDir = '.'
        setupBinaries(options)
        # need to avoid nested container calls, so set toil-inside-toil jobs to local by default
        if "--binariesMode" not in options.cactusOptions:
            options.cactusOptions += " --binariesMode local"
        if options.jobStore.startswith('aws'):
            if not options.outHal.startswith('s3://'):
                raise RuntimeError("--outHal must be s3:// address when using s3 job store")
            if not has_s3:
                raise RuntimeError("S3 support requires toil to be installed with [aws]")
        if options.gpu == 'all':
            if options.batchSystem.lower() in ['single_machine', 'singlemachine']:
                options.gpu = count_nvidia_gpus()
                if not options.gpu:
                    raise RuntimeError('Unable to automatically determine number of GPUs: Please set with --gpu N')
            else:
                raise RuntimeError('--gpu N required in order to use GPUs on non single_machine batch systems')

    options.toil = toil_mode

    if not options.wdl and not options.toil:
        if not options.outDir:
            raise RuntimeError("--outDir option required when not using --wdl")
        if not options.outSeqFile:
            options.outSeqFile = os.path.join(options.outDir, os.path.basename(options.seqFile))
            if os.path.abspath(options.seqFile) == os.path.abspath(options.outSeqFile):
                options.outSeqFile += '.1'
                
    if not options.wdl and options.gpuType != "nvidia-tesla-v100":
        raise RuntimeError("--gpuType can only be used with --wdl ")

    if options.wdl and options.gpu == 'all':
        raise RuntimeError("Number of gpus N must be specified with --gpu N when using --wdl")

    if not options.outHal:
        options.outHal = os.path.join(options.outDir if options.outDir else '', 'out.hal')

    if not options.toil:
        if options.wdl and options.seqFileOnly:
            raise RuntimeError("--wdl cannot be used with --seqFileOnly")
    else:
        options.seqFileOnly = False            

    if options.wdl:
        # wdl handles output file structure
        if options.outDir:
            sys.stderr.write("Warning: --outDir option ignored with --wdl\n")
        options.outDir = "."
        if options.outSeqFile:
            sys.stderr.write("Warning: --outSeqFile option ignored with --wdl\n")
            options.outSeqFile = None
        if options.preprocessOnly:
            raise RuntimeError('--preprocessOnly cannot be used in conjunction with --wdl')
    if not options.dockerImage:
        options.dockerImage = getDockerImage()
    # apply defaults
    if options.defaultCores:
        if not options.preprocessCores:
            options.preprocessCores = options.defaultCores
        if not options.blastCores:
            options.blastCores = options.defaultCores
        if not options.alignCores:
            options.alignCores = options.defaultCores
    if options.defaultMemory:
        if not options.preprocessMemory:
            options.preprocessMemory = options.defaultMemory
        if not options.blastMemory:
            options.blastMemory = options.defaultMemory
        if not options.alignMemory:
            options.alignMemory = options.defaultMemory
    if options.defaultDisk:
        if not options.preprocessDisk:
            options.preprocessDisk = options.defaultDisk
        if not options.blastDisk:
            options.blastDisk = options.defaultDisk
        if not options.alignDisk:
            options.alignDisk = options.defaultDisk
        if not options.halAppendDisk:
            options.halAppendDisk = options.defaultDisk

    # todo: no reason not to support non-1 batch size, but mirror wdl logic for now
    if options.toil:
        if options.preprocessBatchSize != 1:
            if options.preprocessBatchSize != 3:
                # hacky way to only warn for non-default
                sys.stderr.write("Warning: --preprocessBatchSize reset to 1 for --wdl support\n")
            options.preprocessBatchSize = 1
        # todo: could also support this
        assert not options.preprocessOnly

    # https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#gpucount-gputype-and-nvidiadriverversion
    # note: k80 not included as WGA_GPU doesn't run on it.  
    acceptable_gpus = ['nvidia-tesla-v100', 'nvidia-tesla-p100', 'nvidia-tesla-p4', 'nvidia-tesla-t4']
    if options.gpuType not in acceptable_gpus:
        raise RuntimeError('--gpuType {} not supported by Terra.  Acceptable types are {}'.format(
            options.gpuType, acceptable_gpus))

    # used to unique jobstore
    options.jobStoreCount = 0

    enableDumpStack()
    cactusPrepare(options)

def get_jobstore(options, task=None):
    if (options.wdl or options.toil) and '://' not in options.jobStore:
        # if using local jobstore in WDL, it must be relative
        if options.wdl and task:
            prefix = wdl_disk(options, task)[1]
        else:
            prefix = '.'
        return os.path.join(prefix, os.path.basename(options.jobStore.rstrip('/')))
    if '://' not in options.jobStore and options.jobStoreCount == 0:
        # toil won't be able to make a subdir for us if the base dir doesn't exist
        try:
            os.makedirs(options.jobStore)
        except:
            pass            
    # otherwise, make sure the path is unique, since many can be used at once
    js = os.path.join(options.jobStore, str(options.jobStoreCount))
    options.jobStoreCount += 1
    return js

def get_log_options(options, task, event):
    """ get a name of a logfile: this is only relevant for printed commands (ie not toil/wdl) """
    if options.outDir:
        log_base = os.path.join(options.outDir, 'logs')
        try:
            os.makedirs(log_base)
        except:
            pass
        return ' --logFile {}'.format(os.path.join(log_base, '{}-{}.log'.format(task, event)))
    return ''
    
def human2bytesN(s):
    if s is not None:
        sb = human2bytes(s)
        if sb < 10000000:
            raise RuntimeError("Suspiciously small disk or memory specification detected: {}.  Did you forget to add a \"Gi\" for gigabytes?".format(sb))
        return sb
    else:
        return None

def bytes2humanN(s):
    return bytes2human(s).replace(' ', '') if s else s

def bytes2gigs(n):
    return int(int(n)/pow(2,30))

def get_toil_resource_opts(options, task):
    if task == 'preprocess':
        cores = options.preprocessCores
        mem = options.preprocessMemory
    elif task == 'blast':
        cores = options.blastCores
        mem = options.blastMemory
    elif task == 'align':
        cores = options.alignCores
        mem = options.alignMemory
    elif task == 'halAppend':
        cores = 1
        mem = options.alignMemory
    else:
        cores = None
        mem = None
    s = ''
    if cores:
        s += '--maxCores {}'.format(cores)
    if mem and not options.wdl:
        if s:
            s += ' '
        s += '--maxMemory {}'.format(bytes2humanN(mem))
    return s

def wdl_disk(options, task, local_scale=375, local_multiples=[0, 1, 2, 3, 4, 5, 6, 7, 8, 16, 24]):
    """ get the wdl disk options and toil workdir"""
    if task == 'preprocess':
        disk = options.preprocessDisk
    elif task == 'blast':
        disk = options.blastDisk
    elif task == 'align':
        disk = options.alignDisk
    elif task == 'halAppend':
        disk = options.halAppendDisk        
    if not disk:
        return "", "."
    disk = bytes2gigs(disk)
    if disk > local_scale:
        local_multiples = sorted(local_multiples)
        # round up to nearlest local_scale
        for mult in local_multiples:
            if mult * local_scale >= disk:
                disk = mult * local_scale
                break
        if local_multiples and disk > local_multiples[-1] * local_scale:
            logger.critical('Specified disk {} > {}*{}: this probably won\'t run on Terra'.format(disk, local_scale, local_multiples[-1]))            
        disk = int(math.ceil(float(disk) / float(local_scale))) * local_scale
    wdl_disk_options = "local-disk {} LOCAL".format(disk)
    cactus_opts = "."
    return wdl_disk_options, cactus_opts

def get_leaves_and_outgroups(options, mc_tree, og_map, root):
    """ fish the leaves and outgroups out of the experiment xml """
    node = mc_tree.getNodeId(root)
    igs = mc_tree.getChildNames(root)
    ogs = og_map[root] if root in og_map else []
    return igs, ogs

def get_generation_info():
    """ print a comment describing version and command line """
    header = '## generated by : {}\n'.format(" ".join(map(shlex.quote, sys.argv)))
    header += '## date : {}\n'.format(datetime.now())
    header += '## cactus commit : {}\n'.format(cactus_commit)
    return header

def cactusPrepare(options):
    """ annotate a SeqFile with ancestral names as well as paths for output sequences."""

    # read the input
    seqFile = SeqFile(options.seqFile)
    configNode = ET.parse(options.configFile).getroot()
    config = ConfigWrapper(configNode)

    if not options.wdl and not options.toil:
        # prepare output sequence directory
        # todo: support remote (ie s3) output directory
        try:
            os.makedirs(options.outDir)
        except:
            pass
        if not os.path.isdir(options.outDir):
            raise RuntimeError('Unable to create output sequence directory \'{}\''.format(options.outDir))
        if not os.access(options.outDir, os.W_OK):
            logger.warning('Output sequence directory is not writeable: \'{}\''.format(options.outDir))

    if options.preprocessOnly:
        if options.preprocessOnly:
            # hack the configfile to skip preprocessing and write it to the output dir
            config.removePreprocessors()
        options.configFile = os.path.join(options.outDir, 'config-prepared.xml')
        sys.stderr.write("configuration saved in {}\n".format(options.configFile))
        config.writeXML(options.configFile)
        
    # pass through the config file to the options
    # todo (don't like second hard-code check of .xml path)
    if options.configFile != os.path.join(cactusRootPath(), "cactus_progressive_config.xml") and not options.wdl:
        options.cactusOptions += ' --configFile {}'.format(options.configFile)

    # get the ancestor names
    tree = MultiCactusTree(seqFile.tree)
    if not options.seqFileOnly:
        tree.nameUnlabeledInternalNodes(prefix = config.getDefaultInternalNodePrefix())

    if options.includeRoot and not tree.getRootName() in seqFile.pathMap:
        raise RuntimeError('Root genome ({}) must be present in sequfile when --includeRoot is set.'.format(tree.getRootName()))
        
    # make the output
    outSeqFile = SeqFile()
    outSeqFile.tree= tree
    outSeqFile.pathMap = copy.deepcopy(seqFile.pathMap)
    outSeqFile.outgroups = copy.deepcopy(seqFile.outgroups)

    # update paths for preprocessed leaves or inferred ancestors
    for node in outSeqFile.tree.breadthFirstTraversal():
        name = outSeqFile.tree.getName(node)
        leaf = outSeqFile.tree.isLeaf(node)
        if leaf or (not leaf and name not in seqFile.pathMap and not options.preprocessOnly):
            if options.seqFileOnly and name not in seqFile.pathMap:
                continue
            out_basename = seqFile.pathMap[name] if name in seqFile.pathMap else '{}.fa'.format(name)
            outSeqFile.pathMap[name] = os.path.join(options.outDir, os.path.basename(out_basename))
            if options.wdl:
                # uniquify name in wdl to prevent collisions
                outSeqFile.pathMap[name] += '.pp'

    # write the output
    if options.outSeqFile:
        with open(options.outSeqFile, 'w') as out_sf:
            out_sf.write(str(outSeqFile))

    # write the instructions
    if options.toil:
        with Toil(options) as toil:
            if options.restart:
                toil.restart()
            else:
                get_plan(options, seqFile, outSeqFile, config, toil=toil)
    elif not options.seqFileOnly:
        print(get_plan(options, seqFile, outSeqFile, config, toil=None))            

def get_plan(options, inSeqFile, outSeqFile, configWrapper, toil):

    plan = get_generation_info() + '\n'

    if options.wdl:
        plan += wdl_workflow_start(options, inSeqFile)
        options.pp_map = {}

    if options.toil:
        # kick things off with an empty job which we will hook subsequent jobs onto
        # (using RoundedJob because root job must be sublcass of Job,
        #  https://github.com/ComparativeGenomicsToolkit/cactus/pull/284#issuecomment-684125478)
        start_job = RoundedJob()
        parent_job = start_job
        job_idx = {}
    
    # preprocessing
    plan += '\n## Preprocessor\n'
    leaves = [outSeqFile.tree.getName(leaf) for leaf in outSeqFile.tree.getLeaves()]
    for i in range(0, len(leaves), options.preprocessBatchSize):
        pre_batch = leaves[i:i+options.preprocessBatchSize]
        if options.wdl:
            plan += wdl_call_preprocess(options, inSeqFile, outSeqFile, pre_batch)
        elif options.toil:
            job_idx[("preprocess", leaves[i])] = parent_job.addChildJobFn(toil_call_preprocess, options, inSeqFile, outSeqFile, leaves[i],
                                                                          cores=options.preprocessCores,
                                                                          memory=options.preprocessMemory,
                                                                          disk=options.preprocessDisk)
        else:
            plan += 'cactus-preprocess {} {} {} --inputNames {} {} {}{}{}\n'.format(
                get_jobstore(options), options.seqFile, options.outSeqFile, ' '.join(pre_batch),
                options.cactusOptions, get_toil_resource_opts(options, 'preprocess'),
                ' --gpu {}'.format(options.gpu) if options.gpu else '',
                get_log_options(options, 'preprocess', leaves[i]))

    if options.preprocessOnly:
        plan += '\n## Cactus\n'
        plan += 'cactus {} {} {} {}{}\n'.format(get_jobstore(options), options.outSeqFile,
                                                options.outHal, options.cactusOptions,
                                                get_log_options(options, 'preprocess', 'all'))
        return plan

    # shedule up the alignments
    tree = inSeqFile.tree
    mc_tree = MultiCactusTree(tree)
    mc_tree.nameUnlabeledInternalNodes(configWrapper.getDefaultInternalNodePrefix())
    mc_tree.computeSubtreeRoots()
    og_map = compute_outgroups(mc_tree, configWrapper, inSeqFile.outgroups, chrom_info_file = options.chromInfo)
    check_branch_lengths(mc_tree)

    # set of all jobs, as genome names from the (fully resolved, output) seqfile
    events = set(outSeqFile.pathMap.keys()) - set(leaves)
    resolved = set(leaves)

    def get_deps(event):
        subtree = get_subtree(outSeqFile.tree, event, configWrapper, og_map)
        return set([subtree.getName(leaf) for leaf in subtree.getLeaves()])

    # group jobs into rounds.  where all jobs of round i can be run in parallel
    groups = []
    while len(events) > 0:
        group = []
        to_remove = []
        added = 0
        for event in events:
            if all([dep in resolved for dep in get_deps(event)]):
                group.append(event)
                to_remove.append(event)
                added += 1
        if added == 0:
            sys.stderr.write("schedule deadlock:\n")
            for event in events:
                sys.stderr.write("{} has deps {}\n".format(event, get_deps(event)))
            sys.exit(1)
        for tr in to_remove:
            resolved.add(tr)
            events.remove(tr)
        # sort the group so wdl consistent between runs            
        groups.append(sorted(group))

    def halPath(event):
        if event == mc_tree.getRootName():
            return options.outHal
        else:
            return os.path.join(options.outDir, event + '.hal')
    def cigarPath(event):
        return os.path.join(options.outDir, event + '.paf')

    # alignment groups
    plan += '\n## Alignment\n'
    for i, group in enumerate(groups):
        plan += '\n### Round {}'.format(i)
        if options.toil:
            # advance toil phase
            # todo: recapitulate exact dependencies
            parent_job = parent_job.addFollowOn(Job())
        for event in sorted(group):
            plan += '\n'
            if options.wdl:
                plan += wdl_call_blast(options, inSeqFile, mc_tree, og_map, event, cigarPath(event))
                plan += wdl_call_align(options, inSeqFile, mc_tree, og_map, event, cigarPath(event), halPath(event), outSeqFile.pathMap[event])
            elif options.toil:
                # promises only get fulfilleed if they are passed directly as arguments to the toil job, so we pull out the ones we need here
                leaf_deps, anc_deps = get_dep_names(options, mc_tree, og_map, event)
                fa_promises = [job_idx[("preprocess", dep)].rv() for dep in leaf_deps] + [job_idx[("align", dep)].rv(0) for dep in anc_deps]
                job_idx[("blast", event)] = parent_job.addChildJobFn(toil_call_blast,
                                                                     options,
                                                                     outSeqFile,
                                                                     mc_tree,
                                                                     og_map,
                                                                     event,
                                                                     cigarPath(event),
                                                                     leaf_deps + anc_deps,
                                                                     *fa_promises,
                                                                     cores=options.blastCores,
                                                                     memory=options.blastMemory,
                                                                     disk=options.preprocessDisk)
                job_idx[("align", event)] = job_idx[("blast", event)].addFollowOnJobFn(toil_call_align,
                                                                                       options, outSeqFile,
                                                                                       mc_tree,
                                                                                       og_map,
                                                                                       event,
                                                                                       cigarPath(event),
                                                                                       halPath(event),
                                                                                       outSeqFile.pathMap[event],
                                                                                       job_idx[("blast", event)].rv(),
                                                                                       leaf_deps + anc_deps, *fa_promises,
                                                                                       cores=options.alignCores,
                                                                                       memory=options.alignMemory,
                                                                                       disk=options.alignDisk)
            else:
                # todo: support cactus interface (it's easy enough here, but cactus_progressive.py needs changes to handle)
                cactus_options = options.cactusOptions
                if options.includeRoot and event == mc_tree.getRootName():
                    cactus_options += ' --includeRoot'
                if options.chromInfo:
                    cactus_options += ' --chromInfo {}'.format(options.chromInfo)
                plan += 'cactus-blast {} {} {} --root {} {} {}{}{}\n'.format(
                    get_jobstore(options), options.outSeqFile, cigarPath(event), event,
                    cactus_options, get_toil_resource_opts(options, 'blast'),
                    ' --gpu {}'.format(options.gpu) if options.gpu else '',
                    get_log_options(options, 'blast', event))
                plan += 'cactus-align {} {} {} {} --root {} {} {}{}\n'.format(
                    get_jobstore(options), options.outSeqFile, cigarPath(event), halPath(event), event,
                    cactus_options, get_toil_resource_opts(options, 'align'),
                    get_log_options(options, 'align', event))
                # todo: just output the fasta in cactus-align.
                plan += 'hal2fasta {} {} {} > {}\n'.format(halPath(event), event, options.halOptions, outSeqFile.pathMap[event])

    # advance toil phase
    if options.toil:
        parent_job = parent_job.addFollowOn(Job())
                
    # stitch together the final tree
    plan += '\n## HAL merging\n'
    root = mc_tree.getRootName()
    prev_event = None
    append_count = 0
    event_list = []
    for group in reversed(groups):
        for event in group:
            if event != root:
                if not options.toil and not options.wdl:
                    plan += 'halAppendSubtree {} {} {} {} --merge {}\n'.format(
                        halPath(root), halPath(event), event, event, options.halOptions)
                append_count += 1
                event_list.append(event)
            prev_event = event

    if options.toil:
        job_idx['hal_append'] = parent_job.addChildJobFn(toil_call_hal_append_subtrees,
                                                         options,
                                                         mc_tree,
                                                         og_map,
                                                         root,
                                                         job_idx[('align', root)].rv(1),
                                                         event_list,
                                                         *[job_idx[('align', e)].rv(1) for e in event_list],
                                                         cores=1,
                                                         memory=options.alignMemory,
                                                         disk=options.halAppendDisk)

    if options.wdl:
        prev_event = mc_tree.getRootName()
        idx = 0
        while idx < len(event_list):
            event_subset = event_list[idx:idx+options.halAppendBatchSize]
            plan += wdl_call_hal_append(options, mc_tree, og_map, event_subset, prev_event)
            prev_event = event_subset[-1]
            idx += options.halAppendBatchSize
        plan += wdl_workflow_end(options, prev_event, append_count > 1)

    if options.toil:
        start_time = timeit.default_timer()
        toil.start(start_job)
        end_time = timeit.default_timer()
        run_time = end_time - start_time
        logger.info("cactus-prepare-toil has finished after {} seconds".format(run_time))
        
    return plan

def input_fa_name(name):
    """ map event name to input fasta wdl variable name """
    return '{}_fa'.format(name)

def preprocess_call_name(names):
    """ map an event name list to a preprocess call name """
    tag = '_'.join(names)
    # replace all non alphanumeric chars with underscore
    tag = "".join([ c if c.isalnum() else "_" for c in tag ])
    if len(tag) > 80:
        # avoid giant names
        cs = '_' + hashlib.md5(tag.encode()).hexdigest()
        tag = tag[:80-len(cs)] + cs
    return 'preprocess_{}'.format(tag)

def preprocess_output(options, name):
    """ map an event name to its preprocessed output file """
    return '{}.out_files[{}]'.format(options.pp_map[name][0], options.pp_map[name][1])

def blast_call_name(name):
    """ map an event name to a blast call name """
    return 'blast_{}'.format(name)

def align_call_name(name):
    """ map an event name to an align call name """
    return 'align_{}'.format(name)

def hal_append_call_name(name):
    """ map an event name to a hal append call name """
    return 'hal_append_{}'.format(name)

def get_dep_names(options, mc_tree, og_map, event):
    """ return all the ingroup and outgroup event names for which
    the given event needs sequence.  returned in two lists
    (leaves, internals) """    
    leaves, outgroups = get_leaves_and_outgroups(options, mc_tree, og_map, event)
    project_leaves = set([mc_tree.getName(leaf_node) for leaf_node in mc_tree.getLeaves()])
    input_names = list(set(leaves + outgroups))
    leaf_names = []
    anc_names = []
    for input_name in input_names:
        if input_name in project_leaves:
            leaf_names.append(input_name)
        else:
            anc_names.append(input_name)
    return leaf_names, anc_names
    
def wdl_workflow_start(options, in_seq_file):

    s = 'version 1.0\n\n'
    s += wdl_task_preprocess(options) + '\n'
    s += wdl_task_blast(options) + '\n'
    s += wdl_task_align(options) + '\n'
    s += wdl_task_hal_append(options) + '\n'
    
    s += 'workflow cactus_prepared {\n\n'

    # we need to explicitly import local files
    s += '    input {\n'

    s += '        File seq_file'
    if not options.noLocalInputs:
        s += '=\"{}\"'.format(os.path.abspath(options.seqFile))
    s += '\n'
    
    s += '        File? config_file'
    if not options.noLocalInputs and options.configFile != os.path.join(cactusRootPath(), "cactus_progressive_config.xml"):
        s += '=\"{}\"'.format(os.path.abspath(options.configFile))
    s += '\n'

    s += '        File? chrom_info_file'
    if not options.noLocalInputs and options.chromInfo:
        s += '=\"{}\"'.format(os.path.abspath(options.chromInfo))
    s += '\n'
    
    for name, fa_path in in_seq_file.pathMap.items():
        # todo: replace with check from toil
        if '://' not in fa_path:
            s += '        File {}'.format(input_fa_name(name))
            if not options.noLocalInputs:
                s += '=\"{}\"'.format(os.path.abspath(fa_path))
            s += '\n'
    s += '    }\n'
    return s

def wdl_workflow_end(options, event, was_appended):
    s = '\n'
    s += '    output {\n'
    if was_appended:
        s += '        File out_hal = {}.out_file\n'.format(hal_append_call_name(event))
    else:
        s += '        File out_hal = {}.out_hal_file\n'.format(align_call_name(event))
    s += '    }\n'
    s += '}\n'
    return s

def wdl_task_preprocess(options):
    """ preprocess a set of genomes. input either a list of files or a list of urls.  toil can take them both
    as inputs, but they are kept separate due to the distinction between File and String types
    """
    
    s = 'task cactus_preprocess {\n'
    s += '    input {\n'
    s += '        Array[File]? in_files\n'
    s += '        Array[String]? in_names\n'
    s += '        Array[String]? in_urls\n'
    s += '        File? in_config_file\n'
    s += '        Array[String] out_names\n'
    s += '    }\n'
    s += '    command {\n        '
    s += 'cactus-preprocess {} --inPaths ${{sep=\" \" default=\"\" in_files}} ${{sep=\" \" default=\"\" in_urls}}'.format(get_jobstore(options, 'preprocess'))
    s += ' --inputNames ${sep=\" \" default=\"\" in_names}'
    s += ' --outPaths ${{sep=\" \" out_names}} {} {} {}'.format(options.cactusOptions, get_toil_resource_opts(options, 'preprocess'),
                                                                ' --gpu {}'.format(options.gpu) if options.gpu else '')
    s += ' ${\"--configFile \" + in_config_file}'
    s += '\n'
    s += '    }\n'
    
    s += '    runtime {\n'
    s += '        preemptible: {}\n'.format(options.preprocessPreemptible)
    s += '        maxRetries: {}\n'.format(options.preprocessRetries)    
    if options.preprocessCores:
        s += '        cpu: {}\n'.format(options.preprocessCores)
    if options.preprocessMemory:
        s += '        memory: \"{}GB\"\n'.format(bytes2gigs(options.preprocessMemory))
    if options.preprocessDisk:
        s += '        disks: \"{}\"\n'.format(wdl_disk(options, 'preprocess')[0])
    if options.gpu:
        s += '        gpuType: \"{}\"\n'.format(options.gpuType)
        s += '        gpuCount: {}\n'.format(options.gpu)
        s += '        bootDiskSizeGb: 20\n'
        s += '        nvidiaDriverVersion: \"{}\"\n'.format(options.nvidiaDriver)
        s += '        docker: \"{}\"\n'.format(getDockerImage(gpu=True))
        s += '        zones: \"{}\"\n'.format(options.gpuZone)
    else:
        s += '        docker: \"{}\"\n'.format(options.dockerImage)
        s += '        zones: \"{}\"\n'.format(options.zone)
    s += '    }\n'
    s += '    output {\n        Array[File] out_files=out_names\n    }\n'
    s += '}\n'
    return s

def wdl_call_preprocess(options, in_seq_file, out_seq_file, names):

    in_paths = [in_seq_file.pathMap[name] for name in names]
    out_names = [os.path.basename(out_seq_file.pathMap[name]) for name in names]
    # save a link back from the name to its output
    for i, name in enumerate(names):
        assert name not in options.pp_map
        options.pp_map[name] = (preprocess_call_name(names), i)

    s = '    call cactus_preprocess as {} {{\n'.format(preprocess_call_name(names))
    if '://' in in_paths[0]:
        s += '        input: in_urls=[{}],'.format(', '.join(['\"{}\"'.format(in_path) for in_path in in_paths]))
    else:
        s += '        input: in_file=[{}],'.format(', '.join([input_fa_name(name) for name in names]))
    s += ' in_names=[{}],'.format(', '.join(['\"{}\"'.format(name) for name in names]))
    s += ' in_config_file=config_file,'
    s += ' out_names=[{}]'.format(', '.join(['\"{}\"'.format(out_name) for out_name in out_names]))
    s += '\n    }\n'

    return s

def toil_call_preprocess(job, options, in_seq_file, out_seq_file, name):

    work_dir = job.fileStore.getLocalTempDir()
    
    in_path = in_seq_file.pathMap[name]
    out_name = os.path.basename(out_seq_file.pathMap[name])

    cmd = ['cactus-preprocess', os.path.join(work_dir, 'js'), '--inPaths', in_path,
           '--outPaths', out_name, '--inputNames', name, '--workDir', work_dir,
           '--maxCores', str(int(job.cores)), '--maxDisk', bytes2humanN(job.disk), '--maxMemory', bytes2humanN(job.memory)] + options.cactusOptions.strip().split(' ')
    if options.gpu:
        cmd += ['--gpu', options.gpu]
    
    cactus_call(parameters=cmd)

    out_fa_id = job.fileStore.writeGlobalFile(out_name)

    return out_fa_id
    
def wdl_task_blast(options):
    s = 'task cactus_blast {\n'
    s += '    input {\n'
    s += '        File in_seq_file\n'
    s += '        Array[String] in_fa_names\n'
    s += '        Array[File] in_fa_files\n'
    s += '        Array[String]? in_fa_urls\n'
    s += '        String in_root\n'
    s += '        File? in_config_file\n'
    s += '        File? in_chrom_info_file\n'    
    s += '        String? in_options\n'
    s += '        String out_name\n'
    s += '    }\n'
    s += '    command {\n        '
    s += 'cactus-blast {} ${{in_seq_file}} ${{out_name}} --root ${{in_root}}'.format(get_jobstore(options, 'blast'))
    s += ' --pathOverrides ${sep=\" \" in_fa_files} ${sep=\" \" in_fa_urls} --pathOverrideNames ${sep=\" \" in_fa_names}'
    s += ' {} {} {} ${{\"--configFile \" + in_config_file}} ${{in_options}}'.format(options.cactusOptions, get_toil_resource_opts(options, 'blast'),
                                                                                    ' --gpu {}'.format(options.gpu) if options.gpu else '')
    s += ' ${\"--chromInfo \" + in_chrom_info_file}'
    s += '\n    }\n'
    s += '    runtime {\n'
    s += '        preemptible: {}\n'.format(options.blastPreemptible)
    s += '        maxRetries: {}\n'.format(options.blastRetries)
    if options.blastCores:
        s += '        cpu: {}\n'.format(options.blastCores)
    if options.blastMemory:
        s += '        memory: \"{}GB\"\n'.format(bytes2gigs(options.blastMemory))
    if options.blastDisk:
        s += '        disks: \"{}\"\n'.format(wdl_disk(options, 'blast')[0])    
    if options.gpu:
        s += '        gpuType: \"{}\"\n'.format(options.gpuType)
        s += '        gpuCount: {}\n'.format(options.gpu)
        s += '        bootDiskSizeGb: 20\n'
        s += '        nvidiaDriverVersion: \"{}\"\n'.format(options.nvidiaDriver)
        s += '        docker: \"{}\"\n'.format(getDockerImage(gpu=True))
        s += '        zones: \"{}\"\n'.format(options.gpuZone)
    else:
        s += '        docker: \"{}\"\n'.format(options.dockerImage)
        s += '        zones: \"{}\"\n'.format(options.zone)
    s += '    }\n'
    s += '    output {\n        Array[File] out_files=glob(\"${out_name}*\")\n    }\n'
    s += '}\n'
    
    return s

def wdl_call_blast(options, in_seq_file, mc_tree, og_map, event, cigar_name):

    leaf_deps, anc_deps = get_dep_names(options, mc_tree, og_map, event)
    input_fas = []
    input_names = []
    input_urls = []
    for input_name in leaf_deps:
        # take fasta from cactus-preprocess for leaves
        input_fas.append(preprocess_output(options, input_name))
        input_names.append(input_name)
    for input_name in anc_deps:
        # take fasta from cactus-align for internal nodes
        input_fas.append('{}.out_fa_file'.format(align_call_name(input_name)))
        input_names.append(input_name)
    if options.includeRoot and event == mc_tree.getRootName():
        input_urls.append('\"{}\"'.format(in_seq_file.pathMap[event]))
        input_names.append(event)
        extra_options = '--includeRoot'
    else:
        extra_options = ''
        
    s = '    call cactus_blast as {} {{\n'.format(blast_call_name(event))
    s += '        input:'
    s += ' in_seq_file=seq_file,'
    s += ' in_fa_names=[{}],'.format(', '.join(['\"{}\"'.format(name) for name in input_names]))
    s += ' in_fa_files=[{}],'.format(', '.join(input_fas))
    if input_urls:
        s += ' in_fa_urls=[{}],'.format(', '.join(input_urls))
    s += ' in_root=\"{}\",'.format(event)
    s += ' in_config_file=config_file,'
    s += ' in_chrom_info_file=chrom_info_file,'
    if extra_options:
        s += ' in_options=\"{}\",'.format(extra_options)
    s += ' out_name=\"{}\"'.format(os.path.basename(cigar_name))
    s += '\n    }\n'

    return s

def toil_call_blast(job, options, seq_file, mc_tree, og_map, event, cigar_name, dep_names, *dep_fa_ids):

    work_dir = job.fileStore.getLocalTempDir()

    # serialize the seqfile so cactus-blast can use it 
    seq_file_path = os.path.join(work_dir, 'seqfile.txt')
    with open(seq_file_path, 'w') as sf:
        sf.write(str(seq_file))

    # read the fasta files
    assert len(dep_names) == len(dep_fa_ids)
    fa_paths = [os.path.join(work_dir, "{}.pp.fa".format(name)) for name in dep_names]
    for fa_path, fa_id in zip(fa_paths, dep_fa_ids):
        job.fileStore.readGlobalFile(fa_id, fa_path)

    blast_cmd = ['cactus-blast', os.path.join(work_dir, 'js'), seq_file_path, os.path.join(work_dir, os.path.basename(cigar_name)),
                 '--root', event, '--pathOverrides'] + fa_paths+ ['--pathOverrideNames'] + dep_names + \
                ['--workDir', work_dir, '--maxCores', str(int(job.cores)), '--maxDisk', bytes2humanN(job.disk), '--maxMemory', bytes2humanN(job.memory)] + \
                options.cactusOptions.strip().split(' ')
    if options.gpu:
        blast_cmd += ['--gpu', options.gpu]
    if options.chromInfo:
        #todo this won't support distributed execution
        blast_cmd += ['--chromInfo', options.chromInfo]
        
    cactus_call(parameters=blast_cmd)

    # scrape the output files out of the workdir
    out_nameids = []
    for out_file in [f for f in os.listdir(work_dir) if os.path.isfile(os.path.join(work_dir, f))]:
        if out_file.startswith(os.path.basename(cigar_name)):
            out_nameids.append((os.path.basename(out_file), job.fileStore.writeGlobalFile(os.path.join(work_dir, out_file))))
            
    return out_nameids

def wdl_task_align(options):
    s = 'task cactus_align {\n'
    s += '    input {\n'
    s += '        File in_seq_file\n'
    s += '        Array[String] in_fa_names\n'
    s += '        Array[File] in_fa_files\n'
    s += '        Array[String]? in_fa_urls\n'    
    s += '        Array[File] in_blast_files\n'
    s += '        String in_root\n'
    s += '        File? in_config_file\n'
    s += '        File? in_chrom_info_file\n'
    s += '        String? in_options\n'    
    s += '        String out_hal_name\n'
    s += '        String out_fa_name\n'
    s += '    }\n'
    s += '    command {\n        '
    s += 'cactus-align {} ${{in_seq_file}} ${{sep=\" \" in_blast_files}} ${{out_hal_name}} --root ${{in_root}}'.format(get_jobstore(options, 'align'))
    s += ' --pathOverrides ${{sep=\" \" in_fa_files}} ${{sep=\" \" in_fa_urls}} --pathOverrideNames ${{sep=\" \" in_fa_names}} {}'.format(options.cactusOptions)
    s += ' ${\"--chromInfo \" + in_chrom_info_file}'
    s += ' {} ${{\"--configFile \" + in_config_file}} ${{in_options}}'.format(get_toil_resource_opts(options, 'align'))
    s += '\n        '
    s += 'hal2fasta ${{out_hal_name}} ${{in_root}} {} > ${{out_fa_name}}'.format(options.halOptions)
    s += '\n    }\n'
    s += '    runtime {\n'
    s += '        docker: \"{}\"\n'.format(options.dockerImage)
    s += '        preemptible: {}\n'.format(options.alignPreemptible)
    s += '        maxRetries: {}\n'.format(options.alignRetries)    
    if options.alignCores:
        s += '        cpu: {}\n'.format(options.alignCores)
    if options.alignMemory:
        s += '        memory: \"{}GB\"\n'.format(bytes2gigs(options.alignMemory))
    if options.alignDisk:
        s += '        disks: \"{}\"\n'.format(wdl_disk(options, 'align')[0])
    s += '        zones: \"{}\"\n'.format(options.zone)
    s += '    }\n'
    s += '    output {\n'
    s += '        File out_hal_file=\"${out_hal_name}\"\n'
    s += '        File out_fa_file=\"${out_fa_name}\"\n'
    s += '    }\n'
    s += '}\n'
    
    return s

def wdl_call_align(options, in_seq_file, mc_tree, og_map, event, cigar_name, hal_path, fa_path):

    leaf_deps, anc_deps = get_dep_names(options, mc_tree, og_map, event)
    input_fas = []
    input_names = []
    input_urls = []
    for input_name in leaf_deps:
        # take fasta from cactus-preprocess for leaves
        input_fas.append(preprocess_output(options, input_name))
        input_names.append(input_name)
    for input_name in anc_deps:
        # take fasta from cactus-align for internal nodes
        input_fas.append('{}.out_fa_file'.format(align_call_name(input_name)))
        input_names.append(input_name)
    if options.includeRoot and event == mc_tree.getRootName():
        input_urls.append('\"{}\"'.format(in_seq_file.pathMap[event]))        
        input_names.append(event)
        extra_options = '--includeRoot'
    else:
        extra_options = ''
        
    s = '    call cactus_align as {} {{\n'.format(align_call_name(event))
    s += '        input:'
    s += ' in_seq_file=seq_file,'
    s += ' in_fa_names=[{}],'.format(', '.join(['\"{}\"'.format(name) for name in input_names]))
    s += ' in_fa_files=[{}],'.format(', '.join(input_fas))
    if input_urls:
        s += ' in_fa_urls=[{}],'.format(', '.join(input_urls))
    s += ' in_blast_files={}.out_files,'.format(blast_call_name(event))
    s += ' in_root=\"{}\",'.format(event)
    s += ' in_config_file=config_file,'
    s += ' in_chrom_info_file=chrom_info_file,'
    if extra_options:
        s += ' in_options=\"{}\",'.format(extra_options)    
    s += ' out_hal_name=\"{}\",'.format(os.path.basename(hal_path))
    s += ' out_fa_name=\"{}\"'.format(os.path.basename(fa_path))
    s += '\n    }\n'

    return s

def toil_call_align(job, options, seq_file, mc_tree, og_map, event, cigar_name, hal_path, fa_path, blast_output, dep_names, *dep_fa_ids):

    work_dir = job.fileStore.getLocalTempDir()

    # serialize the seqfile so cactus-blast can use it 
    seq_file_path = os.path.join(work_dir, 'seqfile.txt')
    with open(seq_file_path, 'w') as sf:
        sf.write(str(seq_file))

    # download the blast output from the file store
    blast_files = []
    for blast_file_name, blast_file_id in blast_output:
        blast_files.append(os.path.join(work_dir, blast_file_name))
        job.fileStore.readGlobalFile(blast_file_id, blast_files[-1])

    # read the fasta files
    assert len(dep_names) == len(dep_fa_ids)
    fa_paths = [os.path.join(work_dir, "{}.pp.fa".format(name)) for name in dep_names]
    for fa_path, fa_id in zip(fa_paths, dep_fa_ids):
        job.fileStore.readGlobalFile(fa_id, fa_path)

    # call cactus-align
    out_hal_path = os.path.join(work_dir, os.path.basename(hal_path))
    cactus_call(parameters=['cactus-align', os.path.join(work_dir, 'js'), seq_file_path] + blast_files +
                [out_hal_path, '--root', event,
                 '--pathOverrides'] + fa_paths + ['--pathOverrideNames'] + dep_names +
                ['--workDir', work_dir, '--maxCores', str(int(job.cores)), '--maxDisk', bytes2humanN(job.disk), '--maxMemory', bytes2humanN(job.memory)] + options.cactusOptions.strip().split(' ') + (['--chromInfo', options.chromInfo] if options.chromInfo else []))

    out_hal_id = job.fileStore.writeGlobalFile(out_hal_path)

    # export the fasta while we're at it
    out_fa_path = os.path.join(work_dir, '{}.fa'.format(event))
    cactus_call(parameters=['hal2fasta', out_hal_path, event] + options.halOptions.strip().split(' '),
                outfile=out_fa_path)
    out_fa_id = job.fileStore.writeGlobalFile(out_fa_path)

    return out_fa_id, out_hal_id    
    
def wdl_task_hal_append(options):
    s = 'task hal_append_subtree {\n'
    s += '    input {\n'
    s += '        File in_hal_parent\n'
    s += '        Array[File] in_hal_childs\n'
    s += '        Array[String] in_names\n'
    s += '    }\n'
    s += '    String parent_name = basename("${in_hal_parent}")\n'
    s += '    command <<<\n'
    
    # note: I've been unable to modify an input file then return it as an output
    #       so we explicitly copy it into a local string here first
    s += '        cp ~{in_hal_parent} ./~{parent_name}\n'
    # convert WDL arrays to bash arrays (don't think spaces in names will be supported, but that's true almost everywhere else too)
    s += '        HA_CHILDS=(~{sep=" " in_hal_childs})\n'
    s += '        HA_NAMES=(~{sep=" " in_names})\n'    
    s += '        for i in "${!HA_NAMES[@]}"; do\n'
    s += '             halAppendSubtree ./~{{parent_name}} ${{HA_CHILDS[$i]}} ${{HA_NAMES[$i]}} ${{HA_NAMES[$i]}} --merge --hdf5InMemory {}\n'.format(
options.halOptions)
    s += '        done\n'
    s += '    >>>\n'
    s += '    runtime {\n'
    s += '        docker: \"{}\"\n'.format(options.dockerImage)
    s += '        preemptible: {}\n'.format(options.halAppendPreemptible)
    s += '        maxRetries: {}\n'.format(options.halAppendRetries)        
    s += '        cpu: 1\n'
    if options.alignMemory:
        s += '        memory: \"{}GB\"\n'.format(bytes2gigs(options.alignMemory))
    if options.halAppendDisk:
        s += '        disks: \"{}\"\n'.format(wdl_disk(options, 'halAppend')[0])
    s += '        zones: \"{}\"\n'.format(options.zone)
    s += '    }\n'
    s += '    output {\n'
    s += '        File out_file=\"./~{parent_name}\"\n'
    s += '    }\n'
    s += '}\n'

    return s

def wdl_call_hal_append(options, mc_tree, og_map, events, prev_event):

    if prev_event == mc_tree.getRootName():
        # first time: the parent comes out of align
        parent_hal = '{}.out_hal_file'.format(align_call_name(prev_event))
    else:
        # otherwise, it comes from the previous append
        parent_hal = '{}.out_file'.format(hal_append_call_name(prev_event))
    child_hals = ['{}.out_hal_file'.format(align_call_name(event)) for event in events]
    s = '    call hal_append_subtree as {} {{\n'.format(hal_append_call_name(events[-1]))
    s += '        input:'
    s += ' in_hal_parent={},'.format(parent_hal)
    s += ' in_hal_childs=[{}],'.format(",".join(child_hals))
    s += ' in_names=[{}]'.format(",".join(['\"{}\"'.format(event) for event in events]))
    s += '\n    }\n'
    return s

def toil_call_hal_append_subtrees(job, options, mc_tree, og_map, root_name, root_hal_id, event_names, *event_ids):

    work_dir = job.fileStore.getLocalTempDir()

    # donload the root hal file
    root_file = os.path.join(work_dir, '{}.hal'.format(root_name))
    job.fileStore.readGlobalFile(root_hal_id, root_file, mutable=True)
    
    # download the hal files from the file store
    hal_files = []
    for event_name, event_id in zip(event_names, event_ids):
        hal_files.append(os.path.join(work_dir, '{}.hal'.format(event_name)))
        job.fileStore.readGlobalFile(event_id, hal_files[-1])

        # append to the root
        cactus_call(parameters=['halAppendSubtree', root_file, hal_files[-1], event_name, event_name, '--merge', '--hdf5InMemory'] +
                    options.halOptions.strip().split(' '))

    # bypassing toil.exportFile for now as it only works on promises returned by the
    # start job, which isn't how this is set up. also in practice it's often more convenient
    # to output to s3
    # todo: can we just use job.fileStore?
    if options.outHal.startswith('s3://'):
        # write it directly to s3
        write_s3(root_file, options.outHal, region=get_aws_region(options.jobStore))
    else:
        # write the output to disk
        shutil.copy2(root_file,  options.outHal)

    return job.fileStore.writeGlobalFile(root_file)

if __name__ == '__main__':
    main()
