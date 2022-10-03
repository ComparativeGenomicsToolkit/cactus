#!/usr/bin/env python3
import re
import argparse
from os.path import basename, splitext
from bdsg.bdsg import HashGraph, ODGI, PackedGraph
from cyvcf2 import VCF, Writer

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", required=True, help="output VCF file name")
parser.add_argument("graph", help="pangenome graph in HashGraph/ODGI/PackedGraph format")
parser.add_argument("vcffile", help="input biallelic VCF file")
args = parser.parse_args()

def get_steps(path):
    return re.findall("[><][^><]+", path)

def get_at_steps(path):
    return re.findall("[><]\d+", path)

def reverse_path(path):
    d = str.maketrans("><", "<>")
    rev_path = []
    for step in get_steps(path)[::-1]:
        rev_path.append(step.translate(d))
    return "".join(rev_path)

def decompose_traversal(query, target, coords):
    query_steps = get_steps(query)
    target_steps = get_steps(target)
    if not (target_steps[0] == query_steps[0] and target_steps[-1] == query_steps[-1]):
        print(query_steps)
        print(target_steps)
    assert(target_steps[0] == query_steps[0] and target_steps[-1] == query_steps[-1])
    # Remove tailing consecutive identical steps
    for i, (query_step, target_step) in enumerate(zip(query_steps[::-1], target_steps[::-1])):
        if query_step != target_step:
            query_end = len(query_steps) - i + 1
            target_end = len(target_steps) - i + 1
            query_steps = query_steps[:query_end]
            target_steps = target_steps[:target_end]
            coords = coords[:target_end]
            break
    source = (0, 0)
    for i, target_step in enumerate(target_steps[:-1]):
        if i > 0 and target_step in query_steps[source[1]+1:-1]:
            j = query_steps[:-1].index(target_step, source[1]+1)
            target_interval = i - source[0] - 1
            query_interval = j - source[1] - 1
            if target_interval > 0 or query_interval > 0:
                sink = (i, j)
                target_path = "".join(target_steps[source[0]:sink[0]+1])
                target_coord = coords[source[0]:sink[0]+1]
                query_decomposed_steps = query_steps[source[1]:sink[1]+1]
                query_path = "".join(query_decomposed_steps)
                rev_target_steps = get_at_steps(reverse_path(target_path))
                source = (i, j)
                has_inv = False
                # INV should only happen between source and sink steps (exclusive)
                if set(rev_target_steps[1:-1]) & set(get_at_steps(query_path)[1:-1]):
                    has_inv = True
                yield query_path, target_path, target_coord, has_inv
            else:
                source = (i, j)

    # Process the tailing variant
    sink = (len(target_steps) - 1, len(query_steps) - 1)
    target_interval = sink[0] - source[0] - 1
    query_interval = sink[1] - source[1] - 1
    if target_interval > 0 or query_interval > 0:
        target_path = "".join(target_steps[source[0]:sink[0]+1])
        target_coord = coords[source[0]:sink[0]+1]
        query_decomposed_steps = query_steps[source[1]:sink[1]+1]
        query_path = "".join(query_decomposed_steps)
        rev_target_steps = get_at_steps(reverse_path(target_path))
        has_inv = False
        # INV should only happen between source and sink steps (exclusive)
        if set(rev_target_steps[1:-1]) & set(get_at_steps(query_path)[1:-1]):
            has_inv = True
        yield query_path, target_path, target_coord, has_inv

def get_allele_seq(path, graph):
    seq = ""
    # Exclude the first and last steps because they are the same as REF
    for step in get_at_steps(path)[1:-1]:
        strand = step[0]
        node = int(step[1:])
        if strand == ">":
            seq += graph.get_sequence(graph.get_handle(node, False))
        else:
            seq += graph.get_sequence(graph.get_handle(node, True))
    return seq

graph_format = splitext(basename(args.graph))[1]
assert(graph_format in [".hg", ".og", ".pg"])
if graph_format == ".hg":
    graph = HashGraph()
elif graph_format == ".og":
    graph = ODGI()
else:
    graph = PackedGraph()

graph.deserialize(args.graph)

use_ut = False
vcf = VCF(args.vcffile)
w = Writer(args.output, vcf)
w.add_info_to_header({"ID": "SVTYPE", "Number": "1", "Type": "String", "Description": "Type of variant"})
w.add_info_to_header({"ID": "END", "Number": "1", "Type": "Integer", "Description": "End position of the variant described in this record"})
w.add_info_to_header({"ID": "SVLEN", "Number": "1", "Type": "Integer", "Description": "Length of variant"})
w.add_info_to_header({"ID": "SS", "Number": "1", "Type": "String", "Description": "ID of source snarl"})
for variant in vcf:
    chrom = variant.CHROM
    pos = variant.POS
    qual = variant.QUAL
    filt = variant.FILTER
    if filt is None:
        filt = '.'
    level = variant.INFO.get("LV")
    source_snarl = variant.ID
    gt_str = []
    for samp_i in range(len(variant.genotypes)):
        gt_a = []
        for al in variant.genotypes[samp_i][:-1]:
            if al == -1:
                al = '.'
            gt_a.append(str(al))
        if variant.genotypes[samp_i][-1]:
            gt_str.append("|".join(gt_a))
        else:
            gt_str.append("/".join(gt_a))
    gt_str = '\t'.join(gt_str)
    # gt = set(variant.genotypes[0][:-1])
    if variant.INFO.get("UT"):
        use_ut = True
        uts = variant.INFO.get("UT").split(",")
        ats = []
        for ut in uts:
            at = "".join(re.findall("[><]\d+", ut))
            ats.append(at)
    else:
        ats = variant.INFO.get("AT").split(",")
    # There should be only two allele traversals because input VCF is biallelic
    assert(len(ats) == 2)
    ref_at, alt_at = ats

    # Check if this record has been added one preceding base
    ref_at_seq = get_allele_seq(ref_at, graph)
    diff_len = len(variant.REF) - len(ref_at_seq)
    # if not (diff_len == 0 or diff_len == 1):
    #     variant.REF = ref_at_seq
    #     diff_len = 0
    added_onebase = True if diff_len == 1 else False

    ref_steps = get_steps(ref_at)

    # Create an untangled reference traversal
    ref_coords = []
    start = pos
    if added_onebase:
        start += 1
    # Process the first step
    first_step = ref_steps[0]
    strand = first_step[0]
    node = int(first_step[1:])
    if strand == ">":
        seq = graph.get_sequence(graph.get_handle(node, False))
    else:
        seq = graph.get_sequence(graph.get_handle(node, True))
    seq_len = len(seq)
    ref_coords.append((chrom, start - seq_len, seq_len, seq))

    # Process the remaining steps
    for step in ref_steps[1:]:
        strand = step[0]
        node = int(step[1:])
        if strand == ">":
            seq = graph.get_sequence(graph.get_handle(node, False))
        else:
            seq = graph.get_sequence(graph.get_handle(node, True))
        seq_len = len(seq)
        ref_coords.append((chrom, start, seq_len, seq))
        start += seq_len

    # If using INFO/UT field, check if derived coordinates match those in UT
    if use_ut:
        ref_ut, alt_ut = uts
        ut_steps = get_steps(ref_ut)
        for ut_step, coord in zip(ut_steps, ref_coords):
            m = re.search("[><]\d+_(\d+)_(\d+)", ut_step)
            start = int(m[1])
            end = int(m[2])
            length = end - start
            assert(start == coord[1] and length == coord[2])
        alt_traversal = alt_ut
        ref_traversal = ref_ut
    else:
        alt_traversal = alt_at
        ref_traversal = ref_at

    for i, (alt_path, ref_path, ref_coord, has_inv) in enumerate(decompose_traversal(alt_traversal, ref_traversal, ref_coords)):
        ref_seq = get_allele_seq(ref_path, graph)
        alt_seq = get_allele_seq(alt_path, graph)
        chrom, source_start, step_size, source_seq = ref_coord[0]
        variant_start = source_start + step_size
        if len(ref_seq) == len(alt_seq) and not has_inv:
            anchor_start = source_start + step_size
            svlen = len(ref_seq)
            svtype = "SNP" if svlen == 1 else "MNP"
        else:
            if has_inv:
                svlen = len(ref_seq)
                svtype = "INV"
            else:
                svlen = len(alt_seq) - len(ref_seq)
                svtype = "DEL" if svlen < 0 else "INS"
            anchor_start = source_start + step_size - 1
            anchor_end = anchor_start + len(ref_seq)
            ref_seq = source_seq[-1] + ref_seq
            alt_seq = source_seq[-1] + alt_seq

        if svtype in ["SNP", "MNP"]:
            id = f"{chrom}-{variant_start}-{svtype}-{ref_seq}-{alt_seq}"
            if use_ut:
                variant_str = f"{chrom}\t{anchor_start}\t{id}\t{ref_seq}\t{alt_seq}\t{qual:.0f}\t{filt}\tSVTYPE={svtype};UT={ref_path},{alt_path};LV={level};SS={source_snarl}\tGT\t{gt_str}"
            else:
                variant_str = f"{chrom}\t{anchor_start}\t{id}\t{ref_seq}\t{alt_seq}\t{qual:.0f}\t{filt}\tSVTYPE={svtype};AT={ref_path},{alt_path};LV={level};SS={source_snarl}\tGT\t{gt_str}"
        else:
            id = f"{chrom}-{variant_start}-{svtype}-{abs(svlen)}"
            if use_ut:
                variant_str = f"{chrom}\t{anchor_start}\t{id}\t{ref_seq}\t{alt_seq}\t{qual:.0f}\t{filt}\tSVTYPE={svtype};END={anchor_end};SVLEN={svlen};UT={ref_path},{alt_path};LV={level};SS={source_snarl}\tGT\t{gt_str}"
            else:
                variant_str = f"{chrom}\t{anchor_start}\t{id}\t{ref_seq}\t{alt_seq}\t{qual:.0f}\t{filt}\tSVTYPE={svtype};END={anchor_end};SVLEN={svlen};AT={ref_path},{alt_path};LV={level};SS={source_snarl}\tGT\t{gt_str}"

        v = w.variant_from_string(variant_str)
        w.write_record(v)

vcf.close()
w.close()
