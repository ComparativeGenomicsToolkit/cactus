#!/usr/bin/env python3

# Released under the MIT license, see LICENSE.txt

"""Export a Cactus seqfile (plus gzipped per-genome FASTAs) for the subtree rooted at a
given ancestor of a HAL alignment.

For the chosen --root ancestor and every genome below it in the tree, this runs
``hal2fasta | bgzip`` (one Toil job per genome, in parallel) to write
``<outDir>/<genome>.fa.gz``, and writes a seqfile (newick subtree on the first line,
then ``<genome>\\t<path>`` lines) to <seqFile>.  The ancestor's own sequence is included,
so the seqfile is ready for re-aligning the subclade and pinning it back with e.g.
``cactus-blast``/``cactus-align --includeRoot``.

The input HAL is symlinked (not copied) into the jobstore and read by each per-genome
job via symlink, so re-running this on a large HAL doesn't replicate it once per genome.

Example:
    cactus-hal2seqfile ./jobstore mammals.hal ./anc1-fastas ./anc1.seqfile --root Anc1
"""

import os
import sys
import timeit

from sonLib.nxnewick import NXNewick

from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options

from cactus.shared.common import setupBinaries, importSingularityImage, cactus_fast_walltime
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options, add_cactus_toil_options
from cactus.shared.common import makeURL, cactus_call, cactus_clamp_memory
from cactus.shared.version import cactus_commit


def hal_tree_and_lengths(hal_path):
    """Return (newick tree string, {genome: sequence length}) from a single `halStats` call."""
    out = cactus_call(parameters=['halStats', hal_path], check_output=True)
    tree, lengths, in_table = None, {}, False
    for line in out.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('(') and line.endswith(';'):
            tree = line
        elif line.startswith('GenomeName'):
            in_table = True
        elif in_table:
            toks = [t.strip() for t in line.split(',')]
            if len(toks) >= 3:
                lengths[toks[0]] = int(toks[2])
    if tree is None:
        tree = cactus_call(parameters=['halStats', hal_path, '--tree'], check_output=True).strip()
    return tree, lengths


def get_subtree(tree_string, root):
    """Parse the newick tree and return (subtree_newick, [genome names]) for the subtree
    rooted at ``root`` (the ancestor itself plus everything below it, in pre-order)."""
    tree = NXNewick().parseString(tree_string)
    name_to_id = {tree.getName(node): node for node in tree.breadthFirstTraversal()}
    if root not in name_to_id:
        raise RuntimeError('--root "{}" not found in HAL tree.  Available genomes: {}'.format(
            root, ' '.join(sorted(n for n in name_to_id if n))))
    root_id = name_to_id[root]

    def genome_names(node_id):
        names = [tree.getName(node_id)]
        for child in tree.getChildren(node_id):
            names += genome_names(child)
        return names

    def to_newick(node_id):
        children = list(tree.getChildren(node_id))
        if not children:
            return tree.getName(node_id)
        parts = []
        for child in children:
            sub = to_newick(child)
            weight = tree.getWeight(node_id, child, None)
            if weight is not None:
                sub += ':{:g}'.format(weight)
            parts.append(sub)
        return '({}){}'.format(','.join(parts), tree.getName(node_id))

    return to_newick(root_id) + ';', genome_names(root_id)


def export_subtree_fastas(job, hal_id, hal_name, genomes, lengths):
    """Fan out one hal2fasta|bgzip job per genome (run in parallel by Toil).
    Returns a {genome: fasta file id} map."""
    fa_ids = {}
    for genome in genomes:
        # the HAL is read by symlink (below), so each job only needs disk for its own FASTA
        length = lengths.get(genome, hal_id.size)
        fa_ids[genome] = job.addChildJobFn(hal2fasta_gz, hal_id, hal_name, genome,
                                            memory=cactus_clamp_memory(3000000000),
                                            disk=max(int(length * 2), 2**20)).rv()
    return fa_ids


def hal2fasta_gz(job, hal_id, hal_name, genome):
    """hal2fasta <genome> | bgzip, reading the HAL by symlink so that parallel jobs share a
    single copy of the (possibly very large) HAL rather than each copying it out."""
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = os.path.join(work_dir, hal_name)
    job.fileStore.readGlobalFile(hal_id, hal_path, mutable=False, symlink=True)
    fa_path = os.path.join(work_dir, '{}.fa.gz'.format(genome))
    cactus_call(parameters=[['hal2fasta', hal_path, genome], ['bgzip', '--threads', str(job.cores)]],
                outfile=fa_path)
    return job.fileStore.writeGlobalFile(fa_path)


def main():
    parser = Job.Runner.getDefaultArgumentParser()
    add_cactus_toil_options(parser)

    parser.add_argument("halFile", help="input HAL alignment")
    parser.add_argument("outDir", help="output directory for the gzipped per-genome FASTAs (created if needed)")
    parser.add_argument("seqFile", help="output seqfile to write (newick subtree followed by <genome>\\t<path> lines)")
    parser.add_argument("--root", default=None,
                        help="ancestor (internal node) whose subtree -- the node and everything below it -- is "
                        "exported.  Its own sequence is included so the seqfile works with --includeRoot. "
                        "[default: the root of the HAL tree]")

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
    set_logging_from_options(options)
    enableDumpStack()

    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    start_time = timeit.default_timer()

    # Work out the subtree (genomes + newick) and per-genome sizes up front: halStats only reads
    # tree metadata, so this is cheap even for big HALs, and we need the genome list to fan out the
    # per-genome export jobs and the sizes to right-size their disk requirements.
    tree_string, lengths = hal_tree_and_lengths(options.halFile)
    root = options.root if options.root else cactus_call(
        parameters=['halStats', options.halFile, '--root'], check_output=True).strip()
    subtree_newick, genomes = get_subtree(tree_string, root)

    if not os.path.isdir(options.outDir):
        os.makedirs(options.outDir)
    out_dir = os.path.abspath(options.outDir)

    paths = {}
    with Toil(options) as toil:
        importSingularityImage(options)
        if options.restart:
            fa_ids = toil.restart()
        else:
            # symlink the HAL into the jobstore (rather than copying it) -- import_file already
            # defaults to symlink=True in Toil, but we set it explicitly so this holds within the tool
            hal_id = toil.importFile(makeURL(options.halFile), symlink=True)
            fa_ids = toil.start(Job.wrapJobFn(export_subtree_fastas, hal_id,
                                              os.path.basename(options.halFile), genomes, lengths, walltime=cactus_fast_walltime()))

        # export each gzipped fasta to <outDir>/<genome>.fa.gz
        for genome, fa_id in fa_ids.items():
            out_path = os.path.join(out_dir, '{}.fa.gz'.format(genome))
            toil.exportFile(fa_id, makeURL(out_path))
            paths[genome] = out_path

    # write the seqfile: newick subtree, a blank line, then one <genome>\t<path> per line
    with open(options.seqFile, 'w') as seqfile_handle:
        seqfile_handle.write(subtree_newick + '\n\n')
        for genome in genomes:
            seqfile_handle.write('{}\t{}\n'.format(genome, paths[genome]))
    logger.info("Wrote seqfile for subtree rooted at {} ({} genomes) to {}".format(
        root, len(genomes), options.seqFile))

    run_time = timeit.default_timer() - start_time
    logger.info("cactus-hal2seqfile has finished after {} seconds".format(run_time))


if __name__ == '__main__':
    main()
