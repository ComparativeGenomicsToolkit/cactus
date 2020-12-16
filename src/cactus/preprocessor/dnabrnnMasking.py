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
from cactus.shared.common import getOptionalAttrib
from cactus.shared.common import makeURL

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
    def __init__(self, fastaID, dnabrnnOpts, hardmask, cpu, minLength=None):
        memory = 4*1024*1024*1024
        disk = 2*(fastaID.size)
        cores = min(cpu_count(), cpu)
        RoundedJob.__init__(self, memory=memory, disk=disk, cores=cores, preemptable=True)
        self.fastaID = fastaID
        self.minLength = minLength
        self.dnabrnnOpts = dnabrnnOpts
        self.hardmask = hardmask

    def run(self, fileStore):
        """
        mask alpha satellites with dna-brnn
        """
        work_dir = fileStore.getLocalTempDir()
        fastaFile = os.path.join(work_dir, 'seq.fa')
        fileStore.readGlobalFile(self.fastaID, fastaFile)

        # download the model
        modelFile = os.path.join(work_dir, 'model.knm')
        assert os.environ.get("CACTUS_DNA_BRNN_MODEL_ID") is not None        
        modelID = os.environ.get("CACTUS_DNA_BRNN_MODEL_ID")
        fileStore.readGlobalFile(modelID, modelFile)

        # ignore existing model flag
        if '-i' in self.dnabrnnOpts:
            i = self.dnabrnnOpts.index('-i')
            del self.dnabrnnOpts[i]
            del self.dnabrnnOpts[i]

        cmd = ['dna-brnn', fastaFile] + self.dnabrnnOpts.split() + ['-i', modelFile]
        
        if self.cores:
            cmd += ['-t', str(self.cores)]

        bedFile = fileStore.getLocalTempFile()

        # run dna-brnn to make a bed file
        cactus_call(outfile=bedFile, parameters=cmd)

        maskedFile = fileStore.getLocalTempFile()

        mask_cmd = ['cactus_fasta_softmask_intervals.py', '--origin=zero', bedFile]
        if self.minLength:
            mask_cmd += '--minLength={}'.format(self.minLength)

        if self.hardmask:
            mask_cmd += ['--mask=N']

        # do the softmasking
        cactus_call(infile=fastaFile, outfile=maskedFile, parameters=mask_cmd)

        return fileStore.writeGlobalFile(maskedFile)


