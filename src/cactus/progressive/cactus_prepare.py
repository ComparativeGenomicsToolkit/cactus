#!/usr/bin/env python3

#Released under the MIT license, see LICENSE.txt

"""Set up a given seqfile to include proprocessed/ancestral sequences, as well as to

"""
import os, sys
from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import copy
from sonLib.bioio import getTempDirectory
from datetime import datetime

from operator import itemgetter

from cactus.progressive.seqFile import SeqFile
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.common import cactusRootPath
from cactus.shared.common import enableDumpStack
from cactus.shared.common import getDockerImage
from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.shared.configWrapper import ConfigWrapper
from cactus.progressive.schedule import Schedule
from cactus.progressive.projectWrapper import ProjectWrapper
from cactus.shared.version import cactus_commit

def main():
    parser = ArgumentParser()
    parser.add_argument("seqFile", help = "Seq file")
    parser.add_argument("outSeqDir", help='Directory where the processed leaf sequence and ancestral sequences will be placed')
    parser.add_argument("outSeqFile", help = "Path for annotated Seq file output")
    parser.add_argument("outputHal", type=str, help = "Output HAL file")
    parser.add_argument("--configFile", default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))
    parser.add_argument("--preprocessBatchSize", type=int, default=3, help="size (number of genomes) of suggested preprocessing jobs")
    parser.add_argument("--jobStore", type=str, default="./jobstore", help="base directory of jobStores to use in suggested commands")
    parser.add_argument("--halOptions", type=str, default="--hdf5InMemory", help="options for every hal command")
    parser.add_argument("--cactusOptions", type=str, default="--realTimeLogging --logInfo", help="options for every cactus command")
    parser.add_argument("--preprocessOnly", action="store_true", help="only decompose into preprocessor and cactus jobs")
    parser.add_argument("--wdl", action="store_true", help="output wdl workflow instead of list of commands")
    parser.add_argument("--dockerImage", type=str, help="docker image to use as wdl runtime")

    parser.add_argument("--preprocessCores", type=int, help="Number of cores for each cactus-preprocess job")
    parser.add_argument("--blastCores", type=int, help="Number of cores for each cactus-blast job")
    parser.add_argument("--alignCores", type=int, help="Number of cores for each cactus-align job")

    parser.add_argument("--preprocessMem", type=float, help="Memory in GB for each cactus-preprocess job")
    parser.add_argument("--blastMem", type=float, help="Memory in GB for each cactus-blast job")
    parser.add_argument("--alignMem", type=float, help="Memory in GB for each cactus-align job")
    
    options = parser.parse_args()
    options.database = 'kyoto_tycoon'
    #todo support root option
    options.root = None

    if options.wdl:
        if options.preprocessBatchSize != 1:
            sys.stderr.write("Warning: --preprocessBatchSize reset to 1 for --wdl support\n")
            options.preprocessBatchSize = 1
        if options.outSeqDir != ".":
            # todo: generalize
            options.outSeqDir = "."
        if options.preprocessOnly:
            raise RuntimeError('--preprocessOnly cannot be used in conjunction with --wdl')
        # keep toil scratch directory inside cromwell's working directory
        if '--workDir' not in options.cactusOptions:
            options.cactusOptions += ' --workDir .'
    if not options.dockerImage:
        options.dockerImage = getDockerImage()

    # need to go through this garbage (copied from the main() in progressive_cactus) to
    # come up with the project
    options.cactusDir = getTempDirectory()
    #Create the progressive cactus project
    projWrapper = ProjectWrapper(options, options.configFile)
    projWrapper.writeXml()
    # used to unique jobstore
    options.jobStoreCount = 0

    pjPath = os.path.join(options.cactusDir, ProjectWrapper.alignmentDirName,
                          '%s_project.xml' % ProjectWrapper.alignmentDirName)
    assert os.path.exists(pjPath)

    project = MultiCactusProject()

    if not os.path.isdir(options.cactusDir):
        os.makedirs(options.cactusDir)

    project.readXML(pjPath)

    enableDumpStack()
    cactusPrepare(options, project)

def get_jobstore(options):
    if options.wdl and '://' not in options.jobStore:
        # if using local jobstore in WDL, it must be relative
        return os.path.join('.', os.path.basename(options.jobStore.rstrip('/')))
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

def get_toil_resource_opts(options, task):
    if task == 'preprocess':
        cores = options.preprocessCores
        mem = options.preprocessMem
    elif task == 'blast':
        cores = options.blastCores
        mem = options.blastMem
    elif task == 'align':
        cores = options.alignCores
        mem = options.alignMem
    elif task == 'halAppend':
        cores = 1
        mem = options.alignMem
    else:
        cores = None
        mem = None
    s = ''
    if cores:
        s += '--maxCores {}'.format(cores)
    if mem:
        if s:
            s += ' '
        s += '--maxMemory {}G'.format(mem)
    return s

def get_leaves_and_outgroups(options, project, root):
    """ fish the leaves and outgroups out of the experiment xml """
    # open up the experiment (as we do in ProgressiveUp.run)
    experimentFile = project.expMap[root]
    expXml = ET.parse(experimentFile).getroot()
    experiment = ExperimentWrapper(expXml)
    tree = MultiCactusTree(experiment.getTree()).extractSubTree(root)
    leaves = tree.getChildNames(tree.getRootName())
    outgroups = experiment.getOutgroupGenomes()
    return leaves, outgroups

def get_generation_info():
    """ print a comment describing version and command line """
    header = '## generated by : {}\n'.format(' '.join(sys.argv))
    header += '## date : {}\n'.format(datetime.now())
    header += '## cactus commit : {}\n'.format(cactus_commit)
    return header

def cactusPrepare(options, project):
    """ annotate a SeqFile with ancestral names as well as paths for output sequences."""

    # read the input
    seqFile = SeqFile(options.seqFile)
    configNode = ET.parse(options.configFile).getroot()
    config = ConfigWrapper(configNode)

    # prepare output sequence directory
    # todo: support remote (ie s3) output directory
    try:
        os.makedirs(options.outSeqDir)
    except:
        pass
    if not os.path.isdir(options.outSeqDir):
        raise RuntimeError('Unable to create output sequence directory \'{}\''.format(options.outSeqDir))
    if not os.access(options.outSeqDir, os.W_OK):
        logger.warning('Output sequence directory is not writeable: \'{}\''.format(options.outSeqDir))

    # hack the configfile to skip preprocessing and write it to the output dir
    if options.preprocessOnly:
        config.removePreprocessors()
        options.configFile = os.path.join(options.outSeqDir, 'config.xml')
        config.writeXML(options.configFile)
        
    # pass through the config file to the options
    # todo (don't like second hard-code check of .xml path)
    if options.configFile != os.path.join(cactusRootPath(), "cactus_progressive_config.xml"):
        options.cactusOptions += ' --configFile {}'.format(options.configFile)

    # get the ancestor names
    tree = MultiCactusTree(seqFile.tree)
    tree.nameUnlabeledInternalNodes(prefix = config.getDefaultInternalNodePrefix())

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
            out_basename = seqFile.pathMap[name] if name in seqFile.pathMap else '{}.fa'.format(name)
            outSeqFile.pathMap[name] = os.path.join(options.outSeqDir, os.path.basename(out_basename))
            if options.wdl:
                # uniquify name in wdl to prevent collisions
                outSeqFile.pathMap[name] += '.pp'

    # write the output
    with open(options.outSeqFile, 'w') as out_sf:
        out_sf.write(str(outSeqFile))

    # write the instructions
    print(get_plan(options, project, seqFile, outSeqFile))

def get_plan(options, project, inSeqFile, outSeqFile):

    plan = get_generation_info() + '\n'

    if options.wdl:
        plan += wdl_workflow_start(options, inSeqFile)
    
    # preprocessing
    plan += '\n## Preprocessor\n'
    leaves = [outSeqFile.tree.getName(leaf) for leaf in outSeqFile.tree.getLeaves()]
    for i in range(0, len(leaves), options.preprocessBatchSize):
        pre_batch = leaves[i:i+options.preprocessBatchSize]
        if options.wdl:
            assert len(pre_batch) == 1
            plan += wdl_call_preprocess(options, inSeqFile, outSeqFile, leaves[i])
        else:
            plan += 'cactus-preprocess {} {} {} --inputNames {} {} {}\n'.format(
                get_jobstore(options), options.seqFile, options.outSeqFile, ' '.join(pre_batch),
                options.cactusOptions, get_toil_resource_opts(options, 'preprocess'))

    if options.preprocessOnly:
        # if we're only making preprocess jobs, we can end early
        plan += '\n## Cactus\n'
        plan += 'cactus {} {} {} {}\n'.format(get_jobstore(options), options.outSeqFile,
                                              options.outputHal, options.cactusOptions)
        return plan

    # shedule up the alignments
    schedule = Schedule()
    schedule.loadProject(project)
    schedule.compute()

    # set of all jobs, as genome names from the (fully resolved, output) seqfile
    events = set(outSeqFile.pathMap.keys()) - set(leaves)
    resolved = set(leaves)

    # convert follow-ons to dependencies
    follow_on_deps = {}
    for event in events:
        fo = schedule.followOn(event)
        if fo:
            follow_on_deps[fo] = event

    def get_deps(event):
        deps = set(schedule.deps(event))
        if event in follow_on_deps:
            deps = deps.union(set(follow_on_deps[event]))
        # I don't know why the schedule doesn't always give the children
        # todo: understand!
        try:
            has_name = outSeqFile.tree.getNodeId(event) is not None
        except:
            has_name = False
        if has_name:
            for node in outSeqFile.tree.getChildren(outSeqFile.tree.getNodeId(event)):
                if not outSeqFile.tree.isLeaf(node):
                    deps.add(outSeqFile.tree.getName(node))
        return deps

    events_and_virtuals = set()
    for event in events:
        events_and_virtuals.add(event)
        if get_deps(event):
            events_and_virtuals = events_and_virtuals.union(get_deps(event))

    # group jobs into rounds.  where all jobs of round i can be run in parallel
    groups = []
    while len(events_and_virtuals) > 0:
        group = []
        to_remove = []
        added = 0
        for event in events_and_virtuals:
            if all([dep in resolved for dep in get_deps(event)]):
                if not schedule.isVirtual(event):
                    group.append(event)
                to_remove.append(event)
                added += 1
        if added == 0:
            sys.stderr.write("schedule deadlock:\n")
            for event in events_and_virtuals:
                sys.stderr.write("{} has deps {}\b".format(event, get_deps(event)))
            sys.exit(1)
        for tr in to_remove:
            resolved.add(tr)
            events_and_virtuals.remove(tr)
        groups.append(group)

    def halPath(event):
        if event == project.mcTree.getRootName():
            return options.outputHal
        else:
            return os.path.join(options.outSeqDir, event + '.hal')
    def cigarPath(event):
        return os.path.join(options.outSeqDir, event + '.cigar')

    # alignment groups
    plan += '\n## Alignment\n'
    for i, group in enumerate(groups):
        plan += '\n### Round {}'.format(i)
        for event in sorted(group):
            plan += '\n'
            if options.wdl:
                plan += wdl_call_blast(options, project, event, cigarPath(event))
                plan += wdl_call_align(options, project, event, cigarPath(event), halPath(event), outSeqFile.pathMap[event])
            else:
                # todo: support cactus interface (it's easy enough here, but cactus_progressive.py needs changes to handle)
                plan += 'cactus-blast {} {} {} --root {} {} {}\n'.format(
                    get_jobstore(options), options.outSeqFile, cigarPath(event), event,
                    options.cactusOptions, get_toil_resource_opts(options, 'blast'))
                plan += 'cactus-align {} {} {} {} --root {} {} {}\n'.format(
                    get_jobstore(options), options.outSeqFile, cigarPath(event), halPath(event), event,
                    options.cactusOptions, get_toil_resource_opts(options, 'align'))
                # todo: just output the fasta in cactus-align.
                plan += 'hal2fasta {} {} {} > {}\n'.format(halPath(event), event, options.halOptions, outSeqFile.pathMap[event])

    # stitch together the final tree
    plan += '\n## HAL merging\n'
    root = project.mcTree.getRootName()
    prev_event = None
    for group in reversed(groups):
        for event in group:
            if event != root:
                if options.wdl:
                    plan += wdl_call_hal_append(options, project, event, prev_event)
                else:
                    plan += 'halAppendSubtree {} {} {} {} --merge {}\n'.format(
                        halPath(root), halPath(event), event, event, options.halOptions)
            prev_event = event

    if options.wdl:
        plan += wdl_workflow_end(options, prev_event)

    return plan

def input_fa_name(name):
    """ map event name to input fasta wdl variable name """
    return '{}_fa'.format(name)

def preprocess_call_name(name):
    """ map an event name to a preprocess call name """
    return 'preprocess_{}'.format(name)

def blast_call_name(name):
    """ map an event name to a blast call name """
    return 'blast_{}'.format(name)

def align_call_name(name):
    """ map an event name to an align call name """
    return 'align_{}'.format(name)

def hal_append_call_name(name):
    """ map an event name to a hal append call name """
    return 'hal_append_{}'.format(name)

def wdl_workflow_start(options, in_seq_file):

    s = 'version 1.0\n\n'
    s += wdl_task_preprocess(options) + '\n'
    s += wdl_task_blast(options) + '\n'
    s += wdl_task_align(options) + '\n'
    s += wdl_task_hal_append(options) + '\n'
    
    s += 'workflow cactus_prepared {\n\n'

    # we need to explicitly import local files
    s += '    input {\n'
    s += '        File prep_seq_file=\"{}\"\n'.format(os.path.abspath(options.outSeqFile))
    for name, fa_path in in_seq_file.pathMap.items():
        # todo: replace with check from toil
        if '://' not in fa_path:
            s += '        File {}=\"{}\"\n'.format(input_fa_name(name), os.path.abspath(fa_path))
    s += '    }\n'
    return s

def wdl_workflow_end(options, event):
    s = '\n'
    s += '    output {\n'
    s += '        File out_hal = {}.out_file\n'.format(hal_append_call_name(event))
    s += '    }\n'
    s += '}\n'
    return s

def wdl_task_preprocess(options):
    """ preprocess a single genome. expects one of in_file or in_url.  toil can take them both
    as inputs, but they are kept separate due to the distinction between File and String types
    """
    
    s = 'task cactus_preprocess {\n'
    s += '    input {\n'
    s += '        File? in_file\n'
    s += '        String? in_url\n'
    s += '        String out_name\n'
    s += '    }\n'
    s += '    command {\n        '
    s += 'cactus-preprocess {} --inPaths ${{default=\"\" in_file}} ${{default=\"\" in_url}}'.format(get_jobstore(options))
    s += ' --outPaths ${{out_name}} {} {}\n'.format(options.cactusOptions, get_toil_resource_opts(options, 'preprocess'))
    s += '    }\n'
    
    s += '    runtime {\n'
    s += '        docker: \"{}\"\n'.format(options.dockerImage)
    if options.preprocessCores:
        s+= '        cpu: \"{}\"\n'.format(options.preprocessCores)
    if options.preprocessMem:
        s+= '        memory: \"{}GB\"\n'.format(options.preprocessMem)
    s += '    }\n'
    s += '    output {\n        File out_file="${out_name}"\n    }\n'
    s += '}\n'
    return s

def wdl_call_preprocess(options, in_seq_file, out_seq_file, name):

    in_path = in_seq_file.pathMap[name]
    out_name = os.path.basename(out_seq_file.pathMap[name])

    s = '    call cactus_preprocess as {} {{\n'.format(preprocess_call_name(name))
    if '://' in in_path:
        s += '        input: in_url=\"{}\"'.format(in_path)
    else:
        s += '        input: in_file={}'.format(input_fa_name(name))
    s +=', out_name=\"{}\"'.format(out_name)
    s += '\n    }\n'

    return s
    
def wdl_task_blast(options):
    s = 'task cactus_blast {\n'
    s += '    input {\n'
    s += '        File in_seq_file\n'
    s += '        Array[String] in_fa_names\n'
    s += '        Array[File] in_fa_files\n'
    s += '        String in_root\n'
    s += '        String out_name\n'
    s += '    }\n'
    s += '    command {\n        '
    s += 'cactus-blast {} ${{in_seq_file}} ${{out_name}} --root ${{in_root}}'.format(get_jobstore(options))
    s += ' --pathOverrides ${sep=\" \" in_fa_files} --pathOverrideNames ${sep=\" \" in_fa_names}'
    s += ' {} {}'.format(options.cactusOptions, get_toil_resource_opts(options, 'blast'))
    s += '\n    }\n'
    s += '    runtime {\n'
    s += '        docker: \"{}\"\n'.format(options.dockerImage)
    if options.blastCores:
        s+= '        cpu: \"{}\"\n'.format(options.blastCores)
    if options.blastMem:
        s+= '        memory: \"{}GB\"\n'.format(options.blastMem)
    s += '    }\n'
    s += '    output {\n        Array[File] out_files=glob(\"${out_name}*\")\n    }\n'
    s += '}\n'
    
    return s

def wdl_call_blast(options, project, event, cigar_name):

    leaves, outgroups = get_leaves_and_outgroups(options, project, event)
    project_leaves = set([project.mcTree.getName(leaf_node) for leaf_node in project.mcTree.getLeaves()])
    input_names = list(set(leaves + outgroups))
    input_fas = []
    for input_name in input_names:
        # todo: support option to not preprocess
        if input_name in project_leaves:
            # take fasta from cactus-preprocess for leaves
            input_fas.append('{}.out_file'.format(preprocess_call_name(input_name)))
        else:
            # take fasta from cactus-align for internal nodes
            input_fas.append('{}.out_fa_file'.format(align_call_name(input_name)))
            
    s = '    call cactus_blast as {} {{\n'.format(blast_call_name(event))
    s += '        input:'
    s += ' in_seq_file=prep_seq_file,'
    s += ' in_fa_names=[{}],'.format(', '.join(['\"{}\"'.format(name) for name in input_names]))
    s += ' in_fa_files=[{}],'.format(', '.join(input_fas))
    s += ' in_root=\"{}\",'.format(event)
    s += ' out_name=\"{}\"'.format(os.path.basename(cigar_name))
    s += '\n    }\n'

    return s

def wdl_task_align(options):
    s = 'task cactus_align {\n'
    s += '    input {\n'
    s += '        File in_seq_file\n'
    s += '        Array[String] in_fa_names\n'
    s += '        Array[File] in_fa_files\n'
    s += '        Array[File] in_blast_files\n'
    s += '        String in_root\n'    
    s += '        String out_hal_name\n'
    s += '        String out_fa_name\n'
    s += '    }\n'
    s += '    command {\n        '
    s += 'cactus-align {} ${{in_seq_file}} ${{sep=\" \" in_blast_files}} ${{out_hal_name}} --root ${{in_root}}'.format(get_jobstore(options))
    s += ' --pathOverrides ${{sep=\" \" in_fa_files}} --pathOverrideNames ${{sep=\" \" in_fa_names}} {}'.format(options.cactusOptions)
    s += ' {}'.format(get_toil_resource_opts(options, 'align'))
    s += '\n        '
    s += 'hal2fasta ${{out_hal_name}} ${{in_root}} {} > ${{out_fa_name}}'.format(options.halOptions)
    s += '\n    }\n'
    s += '    runtime {\n'
    s += '        docker: \"{}\"\n'.format(options.dockerImage)
    if options.alignCores:
        s+= '        cpu: \"{}\"\n'.format(options.alignCores)
    if options.alignMem:
        s+= '        memory: \"{}GB\"\n'.format(options.alignMem)        
    s += '    }\n'
    s += '    output {\n'
    s += '        File out_hal_file=\"${out_hal_name}\"\n'
    s += '        File out_fa_file=\"${out_fa_name}\"\n'
    s += '    }\n'
    s += '}\n'
    
    return s

def wdl_call_align(options, project, event, cigar_name, hal_path, fa_path):

    leaves, outgroups = get_leaves_and_outgroups(options, project, event)
    project_leaves = set([project.mcTree.getName(leaf_node) for leaf_node in project.mcTree.getLeaves()])
    input_names = list(set(leaves + outgroups))
    input_fas = []
    for input_name in input_names:
        # todo: support option to not preprocess
        if input_name in project_leaves:
            # take fasta from cactus-preprocess for leaves
            input_fas.append('{}.out_file'.format(preprocess_call_name(input_name)))
        else:
            # take fasta from cactus-align for internal nodes
            input_fas.append('{}.out_fa_file'.format(align_call_name(input_name)))

    s = '    call cactus_align as {} {{\n'.format(align_call_name(event))
    s += '        input:'
    s += ' in_seq_file=prep_seq_file,'
    s += ' in_fa_names=[{}],'.format(', '.join(['\"{}\"'.format(name) for name in input_names]))
    s += ' in_fa_files=[{}],'.format(', '.join(input_fas))    
    s += ' in_blast_files={}.out_files,'.format(blast_call_name(event))
    s += ' in_root=\"{}\",'.format(event)
    s += ' out_hal_name=\"{}\",'.format(os.path.basename(hal_path))
    s += ' out_fa_name=\"{}\"'.format(os.path.basename(fa_path))
    s += '\n    }\n'

    return s

def wdl_task_hal_append(options):
    s = 'task hal_append_subtree {\n'
    s += '    input {\n'
    s += '        File in_hal_parent\n'
    s += '        File in_hal_child\n'
    s += '        String in_name\n'
    s += '    }\n'
    s += '    String parent_name = basename("${in_hal_parent}")\n'
    s += '    command {\n'
    
    # note: I've been unable to modify an input file then return it as an output
    #       so we explicitly copy it into a local string here first
    s += '        cp ${in_hal_parent} ${parent_name}\n'
    
    s += '        halAppendSubtree ${{parent_name}} ${{in_hal_child}} ${{in_name}} ${{in_name}} --merge {}'.format(options.halOptions)
    s += '\n    }\n'
    s += '    runtime {\n'
    s += '        docker: \"{}\"\n'.format(options.dockerImage)
    if options.alignMem:
        s+= '        memory: \"{}GB\"\n'.format(options.alignMem)        
    s += '    }\n'
    s += '    output {\n        File out_file=\"${parent_name}\"\n    }\n'
    s += '}\n'

    return s

def wdl_call_hal_append(options, project, event, prev_event):

    if prev_event == project.mcTree.getRootName():
        # first time: the parent comes out of align
        parent_hal = '{}.out_hal_file'.format(align_call_name(prev_event))
    else:
        # otherwise, it comes from the previous append
        parent_hal = '{}.out_file'.format(hal_append_call_name(prev_event))
    child_hal = '{}.out_hal_file'.format(align_call_name(event))
    s = '    call hal_append_subtree as {} {{\n'.format(hal_append_call_name(event))
    s += '        input:'
    s += ' in_hal_parent={},'.format(parent_hal)
    s += ' in_hal_child={},'.format(child_hal)
    s += ' in_name=\"{}\"'.format(event)
    s += '\n    }\n'
    return s
    
if __name__ == '__main__':
    main()
