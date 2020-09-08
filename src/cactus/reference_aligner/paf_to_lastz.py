from toil.common import Toil
from toil.job import Job

import subprocess

def paf_to_lastz(job, paf_file):
    """
    Makes lastz output using paftools.js. Also splits the input paf_file into two files
    in the output, one for the primary and the other for secondary.
    """
    primary = list()
    secondary = list()
    other = list()
    
    print("in paf_to_Lastz - looking for the cg tag.")
    with open(job.fileStore.readGlobalFile(paf_file)) as inf:
        for line in inf:
            print(line)
            if "tp:A:P" in line or "tp:A:I" in line:
                #then the line is a primary output file.
                primary.append(line)
            # elif "tp:A:S" in line:
            else:
                #then the line is a secondary output file.
                secondary.append(line)

    # write output to files; convert to lastz:
    lines = [primary, secondary]
    sort_files = [job.fileStore.getLocalTempFile() for i in range(len(lines))]
    paftool_files = [job.fileStore.getLocalTempFile() for i in range(len(lines))]
    out_files = [job.fileStore.writeGlobalFile(job.fileStore.getLocalTempFile()) for i in range(len(lines))]

    stderr_debug = job.fileStore.getLocalTempFile()
    with open(stderr_debug, "w") as debugf:
        for i in range(len(lines)):
            with open(sort_files[i], "w") as sortf:
                sortf.writelines(lines[i])
            with open(paftool_files[i], "w") as outf:
                subprocess.run(["paftools.js", "view", "-f", "lastz-cigar", sort_files[i]], stdout=outf, stderr=debugf)
            fix_negative_strand_mappings(paftool_files[i], job.fileStore.readGlobalFile(out_files[i]))

    with open(stderr_debug) as debugf:
        for line in debugf:
            print("stderr_debug\t", line)

            
    return out_files

def fix_negative_strand_mappings(infile, outfile):
    """
    paftools.js outputs lastz files with a problem. On "-" strand mappings, the start, stop
    coordinates of the mapping are shown with [smaller coord], [larger coord], when in fact
    the larger coord should be first. This fixes that.
    Note: requires that outfile != infile.
    """
    with open(infile) as inf:
        with open(outfile, "w") as outf:
            for line in inf:
                parsed = line.split()
                if parsed[4] == "-" or parsed[8] == "-":
                    if parsed[4] == "-":
                        first = parsed[3]
                        parsed[3] = parsed[2]
                        parsed[2] = first
                    if parsed[8] == "-":
                        first = parsed[6]
                        parsed[6] = parsed[5]
                        parsed[5] = first
                    outf.write(" ".join(parsed) + "\n")
                else:
                    outf.write(line)
