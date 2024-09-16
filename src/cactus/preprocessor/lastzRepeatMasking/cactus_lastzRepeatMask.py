#!/usr/bin/env python3

## USE LASTZ TO SOFTMASK REPEATS OF A GIVEN FASTA SEQUENCE FILE.

import os
import re
import sys
import shutil

from cactus.shared.common import cactus_cpu_count

from sonLib.bioio import catFiles

from cactus.shared.common import cactus_call
from cactus.shared.common import RoundedJob
from toil.realtimeLogger import RealtimeLogger

class RepeatMaskOptions:
    def __init__(self,
            fragment=200,
            minPeriod=10,
            lastzOpts="",
            unmaskInput=False,
            unmaskOutput=False,
            proportionSampled=1.0,
            gpu=0,
            cpu=None,
            lastz_memory=None,
            gpuLastzInterval=3000000,
            eventName='seq'):
        self.fragment = fragment
        self.minPeriod = minPeriod
        self.lastzOpts = lastzOpts
        self.unmaskInput = unmaskInput
        self.unmaskOutput = unmaskOutput
        self.proportionSampled = proportionSampled
        self.gpu = gpu
        self.cpu = cpu
        self.lastz_memory = lastz_memory
        self.gpuLastzInterval = gpuLastzInterval
        self.eventName = eventName

        self.period = max(1, round(self.proportionSampled * self.minPeriod))

        # make sure fragment size is even so they can overlap by exactly one half.
        if self.fragment % 2:
            self.fragment += 1


class LastzRepeatMaskJob(RoundedJob):
    def __init__(self, repeatMaskOptions, queryID, targetIDs):
        targetsSize = sum(targetID.size for targetID in targetIDs)
        if repeatMaskOptions.lastz_memory:
            memory = repeatMaskOptions.lastz_memory
        elif repeatMaskOptions.gpu:
            memory = min(40 * targetsSize, 512e9)
        else:
            memory = 4*1024*1024*1024
        disk = max(4*(queryID.size + targetsSize), memory)
        cores = repeatMaskOptions.cpu
        accelerators = ['cuda:{}'.format(repeatMaskOptions.gpu)] if repeatMaskOptions.gpu else None            
        RoundedJob.__init__(self, memory=memory, disk=disk, cores=cores, accelerators=accelerators, preemptable=True)
        self.repeatMaskOptions = repeatMaskOptions
        self.queryID = queryID
        self.targetIDs = targetIDs

    def getFragments(self, fileStore, queryFile):
        """
        Chop up the query fasta into fragments of a certain size, overlapping by half their length.
        """
        fragments = os.path.join(self.work_dir, self.repeatMaskOptions.eventName + '_frag')        
        cactus_call(infile=queryFile, outfile=fragments,
                    parameters=["cactus_fasta_fragments.py",
                                "--fragment=%s" % str(self.repeatMaskOptions.fragment),
                                "--step=%s" % (str(self.repeatMaskOptions.fragment // 2)),
                                "--origin=zero"])
        return fragments

    def alignFastaFragments(self, fileStore, targetFiles, fragments):
        """
        Align each query fragment against all the target chunks, stopping
        early to avoid exponential blowup if too many alignments are found.
        """
        target = os.path.join(self.work_dir, self.repeatMaskOptions.eventName + '.fa')
        catFiles(targetFiles, target)
        lastZSequenceHandling  = ['%s[multiple][nameparse=darkspace]' % os.path.basename(target), '%s[nameparse=darkspace]' % os.path.basename(fragments)]
        if self.repeatMaskOptions.unmaskInput:
            lastZSequenceHandling  = ['%s[multiple,unmask][nameparse=darkspace]' % os.path.basename(target), '%s[unmask][nameparse=darkspace]' % os.path.basename(fragments)]
        if self.repeatMaskOptions.gpu:
            assert not self.repeatMaskOptions.unmaskInput
            lastZSequenceHandling = ['%s' % os.path.basename(target), '%s' % os.path.basename(fragments)]
        alignment = os.path.join(self.work_dir, self.repeatMaskOptions.eventName + '.cigar')
        # Each time a fragment aligns to a base in the sequence, that
        # base's match count is incremented.  the plus three for the
        # period parameter is a fudge to ensure sufficient alignments
        # are found
        tool = 'run_kegalign' if self.repeatMaskOptions.gpu else 'lastz'
        # Note that --querydepth has no effect when --ungapped is passed (which is by default)
        lastz_cmd = [tool] + lastZSequenceHandling + \
            self.repeatMaskOptions.lastzOpts.split() + \
            ["--querydepth=keep,nowarn:%i" % (self.repeatMaskOptions.period+3),
             "--format=general:name1,zstart1,end1,name2,zstart2+,end2+",
             "--markend"]
        kegalign_messages = cactus_call(outfile=alignment,
                                        work_dir=self.work_dir,
                                        parameters=lastz_cmd,
                                        job_memory=self.memory,
                                        returnStdErr=self.repeatMaskOptions.gpu,
                                        gpus=self.repeatMaskOptions.gpu)
        if self.repeatMaskOptions.gpu:
            # run_kegalign can crash and still exit 0, so it's worth taking a moment to check the log for errors
            kegalign_messages = kegalign_messages.lower()
            for line in kegalign_messages.split("\n"):
                if not line.startswith("signals delivered"):
                    for keyword in ['terminate', 'error', 'fail', 'assert', 'signal', 'abort', 'segmentation', 'sigsegv', 'kill']:
                        if keyword in line and 'signals' not in line:
                            job.fileStore.logToMaster("KegAlign offending line: " + line)  # Log the messages
                            raise RuntimeError('{} exited 0 but keyword "{}" found in stderr'.format(lastz_cmd, keyword))

            # kegalign will write an invalid file if the output is empty.  correct it here!
            empty_output = False
            with open(alignment, 'r') as alignment_file:
                first_line = alignment_file.readline()
                if first_line.lower().startswith('no alignment'):
                    empty_output = True
                    assert os.path.getsize(alignment) < 32
            if empty_output:
                with open(alignment, 'w') as alignment_file:
                    pass
                
        return alignment

    def maskCoveredIntervals(self, fileStore, queryFile, alignment):
        """
        Mask the query fasta using the alignments to the target. Anything with more alignments than the period gets masked.
        """
        #This runs Bob's covered intervals program, which combines the lastz alignment info into intervals of the query.
        maskInfo = os.path.join(self.work_dir, self.repeatMaskOptions.eventName + '.maskinfo')

        # covered_intervals is part of segalign, so only run if not in gpu mode
        # * 2 takes into account the effect of the overlap
        scale_period = 2

        covered_call_cmd = ["cactus_covered_intervals",
                            "--origin=one",
                            "M=%s" % (int(self.repeatMaskOptions.period * scale_period))]

        covered_call_cmd += ["--queryoffsets"]
        cactus_call(infile=alignment, outfile=maskInfo, parameters=covered_call_cmd, job_memory=self.memory)

        # the previous lastz command outputs a file of intervals (denoted with indices) to softmask.
        # we finish by applying these intervals to the input file, to produce the final, softmasked output.
        args = ["--origin=one"]
        if self.repeatMaskOptions.unmaskOutput:
            args.append("--unmask")
        args.append(os.path.basename(maskInfo))
        maskedQuery = os.path.join(self.work_dir, self.repeatMaskOptions.eventName + '.maskedQeury')
        cactus_call(infile=queryFile, outfile=maskedQuery, work_dir=self.work_dir,
                    parameters=["cactus_fasta_softmask_intervals.py"] + args)
        return maskedQuery

    def run(self, fileStore):
        """
        Using sampled target fragments, mask repetitive regions of the query.
        """
        assert len(self.targetIDs) >= 1
        assert self.repeatMaskOptions.fragment > 1
        self.work_dir = fileStore.getLocalTempDir()
        queryFile = os.path.join(self.work_dir, self.repeatMaskOptions.eventName + '.query')
        fileStore.readGlobalFile(self.queryID, queryFile)
        targetFiles = [os.path.join(self.work_dir, '{}_{}.tgt'.format(self.repeatMaskOptions.eventName, i)) for i in range(len(self.targetIDs))]
        for targetFile, fileID in zip(targetFiles, self.targetIDs):
            fileStore.readGlobalFile(fileID, targetFile)

        fragments = self.getFragments(fileStore, queryFile)
        alignment = self.alignFastaFragments(fileStore, targetFiles, fragments)
        maskedQuery = self.maskCoveredIntervals(fileStore, queryFile, alignment)
        return fileStore.writeGlobalFile(maskedQuery)
