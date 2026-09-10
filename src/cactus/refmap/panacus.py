#!/usr/bin/env python3
"""Run panacus (pangenome coverage/growth statistics with native interactive-HTML plots) on one or
more whole-genome pangenome GFAs, and return the source table data plus a single combined report.

panacus renders its plots directly from the (Rust) binary as a single self-contained HTML file, so
no python plotting stack is needed.  A panacus report is a list of graph sections, so we put every
requested graph type (clip/full/filter) and count type (bp/node/...) into one report -- each as its
own !Gfa section with a distinct graph path (panacus keys sections by graph path, so they must
differ) and a clean human-readable name.

This function is reusable: it takes the merged-GFA dicts and is callable from cactus-graphmap-join,
cactus-pangenome, or a future standalone entry point.

panacus: https://github.com/codialab/panacus  --  please cite panacus if you use these outputs.
"""
import os
from toil.realtimeLogger import RealtimeLogger
from cactus.shared.common import cactus_call, getOptionalAttrib, findRequiredNode

def run_panacus(job, options, config, phase_out_dicts):
    """ Build one panacus report (plus per-graph/per-count-type source TSVs) covering every requested
    graph type.

    :param phase_out_dicts: {phase: make_vg_indexes-output-dict} for each requested phase, where each
           dict holds the merged GFA under the key '<phase>.gfa.gz' (phase is 'clip'/'full'/'filter').
    :returns: a dict of FileIDs keyed for automatic export by graphmap-join (a single
              'panacus.report.html' plus '<phase>.panacus.histgrowth.<count>.tsv' tables).
    """
    work_dir = job.fileStore.getLocalTempDir()

    # config knobs (see the <panacus> element in the config XML)
    pnode = findRequiredNode(config.xmlRoot, "panacus")
    count_types = [c.strip().lower() for c in getOptionalAttrib(pnode, "countType", typeFn=str, default="bp,node").split(',') if c.strip()]
    # dedupe (identical count types would collide on the same graph path) and never end up empty
    count_types = list(dict.fromkeys(count_types)) or ['bp']
    group_by = getOptionalAttrib(pnode, "groupBy", typeFn=str, default="haplotype").lower()
    quorum = getOptionalAttrib(pnode, "quorum", typeFn=str, default="0,0.1,0.5,1")
    cluster_method = getOptionalAttrib(pnode, "clusterMethod", typeFn=str, default="Average")

    # panacus PanSN grouping flag for the histgrowth CLI (-H haplotype / -S sample); the report YAML
    # uses the capitalized form.  count types are lowercase for the CLI, capitalized for the YAML.
    group_flag = {'haplotype': ['-H'], 'sample': ['-S']}.get(group_by, [])
    # the report's !Growth needs coverage and quorum as equal-length paired lists; hold coverage at 1
    # so the quorum values alone select pangenome (0) / shell / core (1) growth curves.
    coverage = ','.join(['1'] * len(quorum.split(',')))
    grouping = group_by.capitalize()
    out_label = getattr(options, 'outName', None) or 'graph'

    # match cactus's exported-filename convention for the graph label: clip is the default and has no
    # label (its '.clip' is stripped on export), full stays 'full', filter becomes 'd<filter>'.
    def phase_label(phase):
        if phase == 'clip':
            return ''
        if phase == 'filter':
            return 'd{}'.format(getattr(options, 'filter', '') or '')
        return phase

    out = {}
    sections = []  # (label, count_type, graph_name) for the report YAML, in order

    # process each requested graph type: decompress its merged GFA (bgzipped) once, then for each
    # count type make a distinctly-named hardlink (panacus keys sections by graph path) and a TSV.
    for phase in sorted(phase_out_dicts):
        label = phase_label(phase)
        gfa_gz = os.path.join(work_dir, '{}.{}.gfa.gz'.format(out_label, phase))
        job.fileStore.readGlobalFile(phase_out_dicts[phase]['{}.gfa.gz'.format(phase)], gfa_gz)
        base_gfa = os.path.join(work_dir, '{}.{}.gfa'.format(out_label, phase))
        cactus_call(parameters=['bgzip', '-d', '-c', gfa_gz], outfile=base_gfa)
        RealtimeLogger.info('panacus: prepared {} graph {}'.format(phase, os.path.basename(base_gfa)))

        for count_type in count_types:
            # a distinct path per (phase, count type) so panacus keeps the sections separate; the name
            # mirrors the exported filenames -- <outName>[.<label>].<count>.gfa (clip has no label)
            graph_name = '.'.join([out_label] + ([label] if label else []) + [count_type, 'gfa'])
            link = os.path.join(work_dir, graph_name)
            if not os.path.exists(link):
                os.link(base_gfa, link)
            tsv = os.path.join(work_dir, '{}.panacus.histgrowth.{}.tsv'.format(phase, count_type))
            cactus_call(parameters=['panacus', 'histgrowth', '-t', str(int(job.cores)),
                                    '-c', count_type, '-q', quorum, '-a'] + group_flag + [graph_name],
                        outfile=tsv, work_dir=work_dir, job_memory=job.memory)
            out['{}.panacus.histgrowth.{}.tsv'.format(phase, count_type)] = job.fileStore.writeGlobalFile(tsv)
            sections.append((label, count_type, graph_name))

    # one self-contained HTML report with a section per (graph type, count type).  panacus ties an
    # analysis' count type to its !Gfa section (only !Hist/!OrderedGrowth/!Similarity can override it,
    # !Growth cannot), so a section per count type is needed.  !Info and !NodeDistribution are
    # count-type independent and describe the graph, so they go in each graph's first section only.
    # panacus otherwise titles a section by its graph path (ugly, and identical across count types);
    # give each an explicit human-readable 'name'.
    yaml_name = 'panacus.report.yaml'
    with open(os.path.join(work_dir, yaml_name), 'w') as yf:
        for label, count_type, graph_name in sections:
            ct = count_type.capitalize()
            suffix = ' ({})'.format(count_type) if len(count_types) > 1 else ''
            name = out_label + (' ' + label if label else '') + suffix
            yf.write('- !Gfa\n')
            yf.write('  graph: {}\n'.format(graph_name))
            yf.write('  name: "{}"\n'.format(name))
            yf.write('  count_type: {}\n'.format(ct))
            yf.write('  grouping: {}\n'.format(grouping))
            yf.write('  analyses:\n')
            if count_type == count_types[0]:
                yf.write('    - !Info\n')
            yf.write('    - !Hist\n      count_type: {}\n'.format(ct))
            yf.write('    - !Growth\n      coverage: {}\n      quorum: {}\n'.format(coverage, quorum))
            yf.write('    - !OrderedGrowth\n      count_type: {}\n      coverage: {}\n      quorum: {}\n'.format(ct, coverage, quorum))
            yf.write('    - !Similarity\n      count_type: {}\n      cluster_method: {}\n'.format(ct, cluster_method))
            if count_type == count_types[0]:
                yf.write('    - !NodeDistribution\n')
    html = os.path.join(work_dir, 'panacus.report.html')
    cactus_call(parameters=['panacus', 'report', '-t', str(int(job.cores)), yaml_name],
                outfile=html, work_dir=work_dir, job_memory=job.memory)
    out['panacus.report.html'] = job.fileStore.writeGlobalFile(html)

    return out
