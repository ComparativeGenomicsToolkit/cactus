from toil.common import Toil
from toil.job import Job

import subprocess

def paf_to_lastz(job, paf_file, sort_secondaries=True):
    """
    Makes lastz output using paftools.js. Also splits the input paf_file into two files
    in the output, one for the primary and the other for secondary.

    sort_secondaries bool, if true, will cause fxn to return two files instead of one.
    
    """
    primary = list()
    primary_mapqs = list()
    secondary = list()
    secondary_mapqs = list()
    
    if not sort_secondaries:
        print("putting all mappings into primary")
        with open(job.fileStore.readGlobalFile(paf_file)) as inf:
            for line in inf:
                primary.append(line)
                primary_mapqs.append(line.split()[11])
    else:
        # print("in paf_to_Lastz - looking for the cg tag.")
        with open(job.fileStore.readGlobalFile(paf_file)) as inf:
            for line in inf:
                print("checking line")
                if "tp:A:S" in line:
                    print("line of interest: ", line)
                if "tp:A:P" in line or "tp:A:I" in line:
                    print("sent to primary")
                    #then the line is a primary output file.
                    primary.append(line)
                    primary_mapqs.append(line.split()[11])
                # elif "tp:A:S" in line:
                else:
                    print("sent to secondary")
                    #then the line is a secondary output file.
                    secondary.append(line)
                    secondary_mapqs.append(line.split()[11])

    # write output to files; convert to lastz:
    lines = [primary, secondary]
    mapqs = [primary_mapqs, secondary_mapqs]
    sort_files = [job.fileStore.getLocalTempFile() for i in range(len(lines))]
    paftool_files = [job.fileStore.getLocalTempFile() for i in range(len(lines))]
    fixed_paftool_files = [job.fileStore.getLocalTempFile() for i in range(len(lines))]
    out_files = [job.fileStore.writeGlobalFile(job.fileStore.getLocalTempFile()) for i in range(len(lines))]

    print("lines in primary:", len(lines[0]))
    print("len(lines in secondary:", len(lines[1]))

    stderr_debug = job.fileStore.getLocalTempFile()
    with open(stderr_debug, "w") as debugf:
        for i in range(len(lines)):
            with open(sort_files[i], "w") as sortf:
                sortf.writelines(lines[i])
            with open(paftool_files[i], "w") as outf:
                subprocess.run(["paftools.js", "view", "-f", "lastz-cigar", sort_files[i]], stdout=outf, stderr=debugf)
            fix_negative_strand_mappings(paftool_files[i], fixed_paftool_files[i])
            add_original_mapqs( mapqs[i], fixed_paftool_files[i], job.fileStore.readGlobalFile(out_files[i]))

    # print("debug info!")
    # with open(stderr_debug) as debugf:
    #     for line in debugf:
    #         print(line)

        
    # check that the lines going into paftools.js are in same order as lines going out.
    with open(job.fileStore.readGlobalFile(out_files[0])) as inf:
        i = 0
        for line in inf:
            #comparing primary from paf to final lastz output.
            paf_parsed = lines[0][i].split()
            lastz_parsed = line.split()
            if (lastz_parsed[3] == "+" and paf_parsed[2] != lastz_parsed[1]) or (lastz_parsed[3] == "-" and paf_parsed[2] != lastz_parsed[2]):
                # print("Lines differ between paf and paftools.js lastz output! Paftools.js may be acting in an unexpected manner. paf line: " + lines[0][i] + " lastz line " + line)
                raise ValueError("Lines differ between paf and paftools.js lastz output! Paftools.js may be acting in an unexpected manner. paf line: " + lines[0][i] + " lastz line " + line)
            i += 1

    # with open(stderr_debug) as debugf:
    #     for line in debugf:
    #         print("stderr_debug\t", line)

    if not sort_secondaries:
        return out_files[0]
    else:
        return out_files

def add_original_mapqs(mapqs, infile, outfile):
    """
    paftools.js uses "raw mapping score" instead of mapq. This replaces those fields with mapq again.
    """
    with open (outfile, "w") as outf:
        with open(infile) as inf:
            i = 0
            for line in inf:
                # this code assumes that paftools.js never drops any mappings (or reorder them) from paf when converting to lastz. Preliminary investigation suggests this is the case. 
                parsed = line.split()
                parsed[9] = mapqs[i]
                i += 1
                outf.write(" ".join(parsed) + "\n")
    if len(mapqs) != i:
        raise ValueError("len(mapqs) != len(lines). Conversion from paf to lastz caused some lines to drop from alignment. len(mapqs): " + str(len(mapqs)) + " len(lines): " + str(len(lines)))


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
