import os
from toil.common import Toil
from toil.job import Job
from cactus.shared.common import cactus_call
from toil.statsAndLogging import logger
from toil.realtimeLogger import RealtimeLogger

def paf_to_lastz(job, paf_file, sort_secondaries=True, mask_bed_id=None, paf_to_stable=False):
    """
    Makes lastz output using paf2lastz. Also splits the input paf_file into two files
    in the output, one for the primary and the other for secondary.

    sort_secondaries bool, if true, will cause fxn to return two files instead of one.
    
    """

    work_dir = job.fileStore.getLocalTempDir()
    paf_path = os.path.join(work_dir, "alignments.paf")
    lastz_path = os.path.join(work_dir, "alignments.cigar")
    secondary_lastz_path = os.path.join(work_dir, "secondary_alignments.cigar")

    job.fileStore.readGlobalFile(paf_file, paf_path)

    cmd = []
    if paf_to_stable:
        cmd.append(['paf2stable', paf_path])
        
    if mask_bed_id:
        mask_bed_path = os.path.join(work_dir, "mask.bed")
        job.fileStore.readGlobalFile(mask_bed_id, mask_bed_path)
        cmd.append(['pafmask', paf_path if not cmd else '-', mask_bed_path])

    paf2lastz_cmd = ['paf2lastz', paf_path if not cmd else '-', '-q']
    if sort_secondaries:
        paf2lastz_cmd += ['-s', secondary_lastz_path]
    cmd.append(paf2lastz_cmd)

    if len(cmd) == 1:
        cmd = cmd[0]

    cactus_call(parameters=cmd, outfile=lastz_path)

    lastz_id = job.fileStore.writeGlobalFile(lastz_path)

    if sort_secondaries:
        secondary_id = job.fileStore.writeGlobalFile(secondary_lastz_path)
        return [lastz_id, secondary_id]
    else:
        return lastz_id
