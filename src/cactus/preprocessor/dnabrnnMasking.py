#!/usr/bin/env python3
"""Uses dna-brnn to mask alpha satellites with a given length threshold
"""

import os
import re
import sys
import shutil

from cactus.shared.common import cactus_cpu_count

from sonLib.bioio import catFiles

from cactus.shared.common import cactus_call
from cactus.shared.common import RoundedJob
from cactus.shared.common import cactusRootPath
from cactus.shared.common import getOptionalAttrib
from cactus.shared.common import makeURL
from cactus.shared.common import get_faidx_subpath_rename_cmd

from toil.realtimeLogger import RealtimeLogger

def loadDnaBrnnModel(toil, configNode, maskAlpha = False):
    """ store the model in a toil file id so it can be used in any workflow """
    for prepXml in configNode.findall("preprocessor"):
        if prepXml.attrib["preprocessJob"] == "dna-brnn":
            if maskAlpha or getOptionalAttrib(prepXml, "active", typeFn=bool, default=False):
                dnabrnnOpts = getOptionalAttrib(prepXml, "dna-brnnOpts", default="")
                if '-i' in dnabrnnOpts:
                    model_path = dnabrnnOpts[dnabrnnOpts.index('-i') + 1]
                else:
                    model_path = os.path.join(cactusRootPath(), 'attcc-alpha.knm')
                os.environ["CACTUS_DNA_BRNN_MODEL_ID"] = toil.importFile(makeURL(model_path))

class DnabrnnMaskJob(RoundedJob):
    def __init__(self, fastaID, dnabrnnOpts, cpu, minLength=None, action=None, inputBedID=None, eventName=None):
        memory = 4*1024*1024*1024
        disk = 2*(fastaID.size)
        cores = min(cactus_cpu_count(), cpu)
        RoundedJob.__init__(self, memory=memory, disk=disk, cores=cores, preemptable=True)
        self.fastaID = fastaID
        self.minLength = minLength
        self.action = action
        self.dnabrnnOpts = dnabrnnOpts        
        self.inputBedID = inputBedID #todo: moved to fileMasking --> remove from here
        self.eventName = eventName

    def run(self, fileStore):
        """
        mask alpha satellites with dna-brnn. returns (masked fasta, dna-brnn's raw output bed, filtered bed used for masking)
        where the filter bed has the minLength filters applied.  When clip is the selected action, suffixes
        get added to the contig names in the format of :<start>-<end> (one-based, inclusive)
        """
        work_dir = fileStore.getLocalTempDir()
        fastaFile = os.path.join(work_dir, 'seq.fa')
        fileStore.readGlobalFile(self.fastaID, fastaFile)

        # download the model
        modelFile = os.path.join(work_dir, 'model.knm')
        assert os.environ.get("CACTUS_DNA_BRNN_MODEL_ID") is not None        
        modelID = os.environ.get("CACTUS_DNA_BRNN_MODEL_ID")
        fileStore.readGlobalFile(modelID, modelFile)

        # download the input bed (which we'll merge in with the bed we compute here)
        if self.inputBedID:
            inputBedFile = os.path.join(work_dir, 'input-regions.bed')
            fileStore.readGlobalFile(self.inputBedID, inputBedFile)
        else:
            inputBedFile = None

        # ignore existing model flag
        if '-i' in self.dnabrnnOpts:
            i = self.dnabrnnOpts.index('-i')
            del self.dnabrnnOpts[i]
            del self.dnabrnnOpts[i]

        cmd = ['dna-brnn', fastaFile] + self.dnabrnnOpts.split() + ['-i', modelFile]
        
        if self.cores:
            cmd += ['-t', str(self.cores)]

        bedFile = os.path.join(work_dir, 'regions.bed')

        # run dna-brnn to make a bed file
        cactus_call(outfile=bedFile, parameters=cmd)

        if self.minLength is None:
            self.minLength = 0

        # load the fasta sequence information (needed for below clipping and/or bed merging)
        if inputBedFile or self.action == "clip":
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
        merge_cmd.append(['awk', '{{if($3-$2 > {}) print}}'.format(self.minLength), bedFile])
        merge_cmd.append(['bedtools', 'sort', '-i', '-'])
        merge_cmd.append(['bedtools', 'merge', '-i', '-', '-d', str(self.minLength)])            
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


def computePAFCoverage(job, config_node, paf_id):
    """ compute the gaps in PAF coverage, store them as a bed file, and add the bed file's filestore id
    into the config's dna-brnn xml element """
    paf_file = job.fileStore.readGlobalFile(paf_id)
    bed_file = job.fileStore.getLocalTempFile()

    dnabrnn_node = None
    for node in config_node.findall("preprocessor"):
        if getOptionalAttrib(node, "preprocessJob") == 'dna-brnn':
            dnabrnn_node = node
            break

    assert dnabrnn_node is not None

    min_length = max(1, getOptionalAttrib(dnabrnn_node, 'minLength', typeFn=int, default=0))

    cactus_call(parameters=['pafcoverage', paf_file, '-g', '-m', str(min_length)], outfile=bed_file)

    dnabrnn_node.attrib["inputBedID"] = job.fileStore.writeGlobalFile(bed_file)

    return config_node
