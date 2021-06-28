#!/usr/bin/env python3
"""Cut prefixes or suffixes from fasta headeres
"""


import os
import re
import sys
import shutil

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from cactus.shared.common import RoundedJob
from toil.realtimeLogger import RealtimeLogger

class CutHeadersJob(RoundedJob):
    def __init__(self, fastaID, cutBefore, cutBeforeOcc, cutAfter):
        disk = 2*(fastaID.size)
        RoundedJob.__init__(self, disk=disk, preemptable=True)
        self.fastaID = fastaID
        self.cutBefore = cutBefore
        self.cutBeforeOcc = cutBeforeOcc
        self.cutAfter = cutAfter

    def run(self, fileStore):
        """
        Cut before cutBefore and after cutAfter
        
        If cutBefore is # then something like
        >HG02055#1#h1tg000001l
        would become 
        >h1tg000001l

        If cutBeforeOcc is specified, it will only cut up to cutBeforeOcc characters
        so cutBefore # with cutBeforeOcc = 2 would change
        >HG02055#1#h1tg000001l#EBV
        into
        >h1tg000001l#EBV

        If cutAfter is a whitespace, then something like 
        >chr1  AC:CM000663.2  gi:568336023  LN:248956422  rl:Chromosome  M5:6aef897c3d6ff0c78aff06
        would become
        >chr1

        """
        work_dir = fileStore.getLocalTempDir()
        input_path = os.path.join(work_dir, 'seq.fa')
        fileStore.readGlobalFile(self.fastaID, input_path)
        output_path = os.path.join(work_dir, 'seq.cut.fa')

        with open(input_path, 'r') as in_file, open(output_path, 'w') as out_file:
            for seq_record in SeqIO.parse(in_file, 'fasta'):
                header = seq_record.description
                if self.cutBefore:
                    occs = [i for i, c in enumerate(header) if c in self.cutBefore]
                    if occs:
                        if not self.cutBeforeOcc:
                            pos = occs[-1]
                        else:
                            pos = occs[min(len(occs), self.cutBeforeOcc) - 1]
                        if pos >= 0:
                            if pos < len(header) - 1:
                                header = header[pos + 1:]
                            else:
                                header = ""
                if self.cutAfter:
                    pos_list = [header.find(c) for c in self.cutAfter if header.find(c) >= 0]
                    if pos_list:
                        pos = min(pos_list)
                        header = header[0:pos]

                if not header:
                    raise RuntimeError("Error: applying cutHeaders preprocessor removes entire header: {}".format(seq_record.description))
                seq_record.description = header
                seq_record.id = header
                SeqIO.write(seq_record, out_file, 'fasta')

        return fileStore.writeGlobalFile(output_path)
