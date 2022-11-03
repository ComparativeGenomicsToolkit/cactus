#!/usr/bin/env python3
"""Pass through an input BED (or PAF) file and mask directly with that, bypassing all usual
preprocessor functionality.  It's a bit ugly, but currently needed to the pangenome pipeline
in order to mask coverage gaps (todo: clipping them in the PAF itself would be much cleaner)
"""

import os
import re
import sys
import shutil
import xml.etree.ElementTree as ET
from cactus.shared.common import cactus_cpu_count

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from sonLib.bioio import catFiles

from cactus.shared.common import cactus_call
from cactus.shared.common import RoundedJob
from cactus.shared.common import getOptionalAttrib
from cactus.shared.common import makeURL
from cactus.shared.common import get_faidx_subpath_rename_cmd

from toil.realtimeLogger import RealtimeLogger

class FileMaskingJob(RoundedJob):
    def __init__(self, fastaID, minLength=None, action=None, inputBedID=None, eventName=None):
        disk = 2*(fastaID.size)
        memory = fastaID.size
        RoundedJob.__init__(self, disk=disk, memory=memory, preemptable=True)
        self.fastaID = fastaID
        self.minLength = minLength
        self.action = action
        self.inputBedID = inputBedID
        self.eventName = eventName

    def run(self, fileStore):
        """
        extract any existing masking, merge it with the input bed, then apply it to the fasta
        """
        work_dir = fileStore.getLocalTempDir()
        fastaFile = os.path.join(work_dir, 'seq.fa')
        fileStore.readGlobalFile(self.fastaID, fastaFile)

        # download the input bed (which we'll merge in with the bed we compute here)
        inputBedFile = os.path.join(work_dir, 'input-regions.bed')
        fileStore.readGlobalFile(self.inputBedID, inputBedFile)

        if self.minLength is None:
            self.minLength = 0

        # extract the existing masked regions to merge in
        bedFile = get_mask_bed_from_fasta(self, self.eventName, None, fastaFile, self.minLength, work_dir)

        # load the fasta sequence information (needed for below clipping and/or bed merging)
        cactus_call(parameters=['samtools', 'faidx', fastaFile])
        # load the contig lengths
        contig_lengths = {}
        with open(fastaFile + '.fai', 'r') as fai:
            for line in fai:
                toks = line.strip().split('\t')
                contig_lengths[toks[0]] = int(toks[1])

        # merge in the input bed file
        if inputBedFile:
            input_line_count = 0
            with open(inputBedFile, 'r') as inputBedStream, open(bedFile, 'a') as bedStream:
                if self.eventName:
                    eventPrefix = 'id={}|'.format(self.eventName)
                else:
                    eventPrefix = ''
                for line in inputBedStream:
                    toks = line.split('\t')
                    if toks:
                        # our PAF file probably has prefixes like id=EVENT| which won't match up to the fasta
                        # so we strip here:
                        from_event = toks[0].startswith(eventPrefix)
                        if from_event:
                            toks[0] = toks[0][len(eventPrefix):]
                        # we may have given a whole-genome paf, so filter down our bed to the
                        # relevant contigs for this fasta (won't change masking output, but bed output will be cleaner)
                        if toks[0] in contig_lengths:
                            assert from_event
                            bedStream.write('\t'.join(toks))
                            input_line_count += 1
            RealtimeLogger.info("Merged in {} bed lines from input bed file for eventPrefix=\"{}\"".format(input_line_count, eventPrefix))
                            
        # merge up the intervals into a new bed file
        mergedBedFile = os.path.join(work_dir, 'filtered.bed')
        merge_cmd = []
        merge_cmd.append(['awk', '{{if($3-$2 > {}) print $1\"\\t\"$2\"\\t\"$3}}'.format(self.minLength), bedFile])
        merge_cmd.append(['bedtools', 'sort', '-i', '-'])
        merge_cmd.append(['bedtools', 'merge', '-i', '-', '-d', str(self.minLength)])
        if self.eventName:
            merge_cmd.append(['sed', '-e', 's/id={}|//g'.format(self.eventName)])
        cactus_call(outfile=mergedBedFile, parameters=merge_cmd)

        maskedFile = os.path.join(work_dir, 'masked.fa')
        
        if self.action in ('softmask', 'hardmask'):
            mask_cmd = ['cactus_fasta_softmask_intervals.py', '--origin=zero', mergedBedFile]
            if self.minLength:
                mask_cmd += ['--minLength={}'.format(self.minLength)]
            if self.action == 'hardmask':
                mask_cmd += ['--mask=N']
            # do the softmasking
            cactus_call(infile=fastaFile, outfile=maskedFile, parameters=mask_cmd)
        else:
            assert self.action == "clip"
            # to clip, we need a bed of the regions we want to *keep*.  We'll start with the whole thing
            allRegionsFile = os.path.join(work_dir, 'chroms.bed')
            cactus_call(outfile=allRegionsFile, parameters=['awk', '{print $1 "\\t0\\t" $2}', fastaFile + '.fai'])
            # now we cut out the regions
            clippedRegionsFile = os.path.join(work_dir, 'clipped.bed')
            cactus_call(outfile=clippedRegionsFile, parameters=['bedtools', 'subtract', '-a', allRegionsFile, '-b', mergedBedFile])
            # now we make a fiadx input regions
            faidxRegionsFile = os.path.join(work_dir, 'faidx_regions.txt')
            with open(clippedRegionsFile, 'r') as clipFile, open(mergedBedFile, 'a') as mergeFile, open(faidxRegionsFile, 'w') as listFile:
                for line in clipFile:
                    toks = line.strip().split("\t")
                    if len(toks) > 2:
                        seq, start, end = toks[0], int(toks[1]), int(toks[2])
                        if end - start > self.minLength or contig_lengths[seq] <= self.minLength:
                            region = seq
                            if end - start < contig_lengths[seq]:
                                # go from 0-based end exlusive to 1-based end inclusive when
                                # converting from BED to samtools region
                                region += ':{}-{}'.format(start + 1, end)
                            else:
                                assert start == 0 and end == contig_lengths[seq]
                            listFile.write('{}\n'.format(region))
                        else:
                            # the region was too small, we remember it in our filtered bed file
                            mergeFile.write(line)
            # and cut the fasta apart with samtools
            cmd = [['samtools', 'faidx', fastaFile, '-r', faidxRegionsFile]]
            # transform chr1:10-15 (1-based inclusive) into chr1_sub_9_15 (0-based end open)
            # this is a format that contains no special characters in order to make assembly hubs
            # happy.  But it does require conversion going into vg which wants chr[9-15] and
            # hal2vg can is updated to do this autmatically
            cmd.append(get_faidx_subpath_rename_cmd())
            cactus_call(outfile=maskedFile, parameters=cmd)
        
        return fileStore.writeGlobalFile(maskedFile), fileStore.writeGlobalFile(bedFile), fileStore.writeGlobalFile(mergedBedFile)


def maskJobOverride(job, config_node, mask_file_path, mask_file_id, mask_file_action, min_length):
    """ return a hijacked config file that does just one preprocessing job: mask each fasta sequence with 
    the given bed file.  if paf_length is specified, the file is treated as a PAF file, and a BED is extracted
    from it using coverage gaps of at least the given length. 
    """
    # this was unzipped upstream
    if mask_file_path.endswith('.gz'):
        mask_file_path = mask_file_path[:-3]
        
    if mask_file_path.endswith('.paf'):
        # convert the PAF to BED
        paf_file = job.fileStore.readGlobalFile(mask_file_id)
        bed_file = job.fileStore.getLocalTempFile()

        if not min_length:
            min_length = 1

        cactus_call(parameters=['pafcoverage', paf_file, '-g', '-m', str(min_length)], outfile=bed_file)

        mask_file_id = job.fileStore.writeGlobalFile(bed_file)

    # rewrite the config
    for node in config_node.findall("preprocessor"):
        config_node.remove(node)
        
    mask_node = ET.SubElement(config_node, 'preprocessor')
    mask_node.attrib['preprocessJob'] = 'maskFile'
    mask_node.attrib['inputBedID'] = mask_file_id
    mask_node.attrib['action'] = mask_file_action
    mask_node.attrib['minLength'] = min_length

    return config_node

def get_mask_bed_from_fasta(job, event, fa_id, fa_path, min_length, work_dir = None):
    """ make a bed file from one fasta"""
    return_id = False # hack in a toggle (work_dir) that lets this be called as a job or a function
    if not work_dir:
        work_dir = job.fileStore.getLocalTempDir()
        return_id = True
    bed_path = os.path.join(work_dir, os.path.basename(fa_path) + '.mask.bed')
    fa_path = os.path.join(work_dir, os.path.basename(fa_path))
    is_gz = fa_path.endswith(".gz")
    if return_id:
        job.fileStore.readGlobalFile(fa_id, fa_path, mutable=is_gz)
    if is_gz:
        cactus_call(parameters=['gzip', '-fd', fa_path])
        fa_path = fa_path[:-3]
    cactus_call(parameters=['cactus_softmask2hardmask', fa_path, '-b', '-m', str(min_length)], outfile=bed_path)
    if return_id:
        return job.fileStore.writeGlobalFile(bed_path)
    else:
        return bed_path
