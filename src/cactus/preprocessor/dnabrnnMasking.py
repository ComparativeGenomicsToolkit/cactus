#!/usr/bin/env python3
"""Uses dna-brnn to mask alpha satellites with a given length threshold
"""

import os
import re
import sys
import shutil

from toil.lib.threading import cpu_count

from sonLib.bioio import catFiles

from cactus.shared.common import cactus_call
from cactus.shared.common import RoundedJob
from cactus.shared.common import cactusRootPath
from toil.realtimeLogger import RealtimeLogger

class DnabrnnMaskJob(RoundedJob):
    def __init__(self, fastaID, minLength, dnabrnnOpts, hardmask):
        memory = 4*1024*1024*1024
        disk = 2*(fastaID.size)
        # todo: clean up
        cores = cpu_count()
        RoundedJob.__init__(self, memory=memory, disk=disk, cores=cores, preemptable=True)
        self.fastaID = fastaID
        self.minLength = minLength
        self.dnabrnnOpts = dnabrnnOpts
        self.hardmask = hardmask

    def run(self, fileStore):
        """
        mask alpha satellites with dna-brnn
        """
        fastaFile = fileStore.readGlobalFile(self.fastaID)

        cmd = ['dna-brnn', fastaFile] + self.dnabrnnOpts.split()
        if '-i' not in self.dnabrnnOpts:
            # pull up the model
            # todo: is there are more robust way?
            cmd += ['-i', os.path.join(cactusRootPath(), 'attcc-alpha.knm')]
        
        if self.cores:
            cmd += ['-t', str(self.cores)]

        bedFile = fileStore.getLocalTempFile()

        # run dna-brnn to make a bed file
        cactus_call(outfile=bedFile, parameters=cmd)

        maskedFile = fileStore.getLocalTempFile()

        mask_cmd = ['cactus_fasta_softmask_intervals.py', '--origin=zero', '--minLength={}'.format(self.minLength), bedFile]

        if self.hardmask:
            mask_cmd += ['--mask=N']

        # do the softmasking
        cactus_call(infile=fastaFile, outfile=maskedFile, parameters=mask_cmd)

        return fileStore.writeGlobalFile(maskedFile)


