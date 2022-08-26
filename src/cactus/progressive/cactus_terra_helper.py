#!/usr/bin/env python3

#Released under the MIT license, see LICENSE.txt

"""Terra's workflow resume functionality absolutely stinks (I may be a little bitter atm):  Say you are running a large Cactus WDL job, and it's gone for many weeks and thousands of $$$, but something happens causing a job to fail.  You update the WDL to increase the memory, or switch out the Docker image (or even just rerun with the same image that you updated on quay) etc and, 9 times out of 10, Terra will restart the whole thing from scratch.  

This is a quick and dirty script to take the WDL you want to rerun and, given a google bucket prefix, fish out all the intermediate results and inline them in the WDL.

Example:  

1) run primates.wdl on Terra
2) it crashes in blast_Anc07
3) make primates.2.wdl that uses an updated docker image
4) fish out the google bucket prefix, which would look like gs://fc-8d260e46-b01a-4bc0-8710-b2e95124dece/d3ec9c38-3794-4590-8452-0f367b265495/cactus_prepared/969c5d70-c0ee-4238-b6f7-5b9df3d1f041/
5) gsutil ls -r <bucket> | cactus-terra-helper resume primates.wdl > primates.2.wdl  

and it will print out a new WDL where things like preprocess_Piliocolobus_tephrosceles_Pan_paniscus_Pan_trog_4b5e36051c04efdc81c34c6206d44ee1.out_files[2] are replaced with the full path gs://fc-8d260e46-b01a-4bc0-8710-b2e95124dece/d3ec9c38-3794-4590-8452-0f367b265495/cactus_prepared/969c5d70-c0ee-4238-b6f7-5b9df3d1f041/call-preprocess_Piliocolobus_tephrosceles_Pan_paniscus_Pan_trog_4b5e36051c04efdc81c34c6206d44ee1/attempt-2/PD_0802_megahit_final_contigs.min1k.assembly.masked.fa.pp etc.

Warning:

The parsing is pretty quick and hacky and really expects the WDL to look a lot like what comes out of cactus-prepare --wdl (at time of writing), and the terra bucket structure to be what I'm currently seeing.  Also, it assumes that if an output file is found in the job's bucket, then that job was a success and the file should be used (which I think is fine)

Scraping logs:

Use the scrape-logs subcommand to download all the logs
"""

import os, sys
from datetime import datetime
import subprocess

def scrape_logs(dirtree):
    """ download all the logs (stderr).  input must come from gsutil ls -l -r """    
    job_to_log = {}
    for line in dirtree:
        line_toks = line.strip().split()
        if len(line_toks) != 3:
            continue
        full_path = line_toks[-1]        
        if not full_path.startswith("gs://"):
            continue
        if os.path.basename(full_path) != 'stderr':
            continue
        toks = full_path.split("/")
        pi = toks.index("cactus_prepared")
        if pi + 3 >= len(toks):
            continue
        job_name = toks[pi + 2]
        assert job_name.startswith("call-")
        job_name = job_name[5:]
        
        assert len(line_toks) == 3
        size = int(line_toks[0])
        date = datetime.strptime(line_toks[1], "%Y-%m-%dT%H:%M:%SZ")

        # note we accept older logs if they are much bigger to filter out empty logs due to restart shenanigans
        if job_name not in job_to_log or \
           (date > job_to_log[job_name][0] and size * 10 > job_to_log[job_name][1]) or \
           (date <= job_to_log[job_name][0] and size > 10 * job_to_log[job_name][1]):
            job_to_log[job_name] = (date, size, full_path)

    for k,v in job_to_log.items():
        subprocess.check_call(['gsutil', 'cp', v[2], './{}.log'.format(k)])

def load_dirtree(dirtree):
    """ load up the results of gsutil ls -r"""
    # scab the gsutil ls -r results to make a mapping of job name to (full path) of output files
    pp_files = {} 
    blast_files = {}
    align_files = {}
    append_files = {}
    for line in dirtree:
        if not line.startswith("gs://"):
            continue
        full_path = line.strip()
        file_name = os.path.basename(full_path)
        toks = full_path.split("/")
        pi = toks.index("cactus_prepared")
        if pi + 3 >= len(toks):
            continue
        job_name = toks[pi + 2]
        assert job_name.startswith("call-")
        
        if job_name.startswith('call-preprocess'):
            job_name = job_name[5:]
            if file_name.endswith('.pp'):
                if job_name not in pp_files:
                    pp_files[job_name] = []
                pp_files[job_name].append(full_path)

        elif job_name.startswith('call-blast'):
            job_name = job_name[5:]
            if full_path.find('.cigar') > 0:
                if job_name not in blast_files:
                    blast_files[job_name] = []
                blast_files[job_name].append(full_path)

        elif job_name.startswith('call-align'):
            job_name = job_name[5:]
            if full_path.endswith('.hal') or full_path.endswith('.pp') or full_path.endswith('.fa'):
                if job_name not in align_files:
                    align_files[job_name] = []
                align_files[job_name].append(full_path)

        elif job_name.startswith('call-hal_append'):
            job_name = job_name[5:]
            if full_path.endswith('.hal'):
                if job_name not in append_files:
                    append_files[job_name] = []
                append_files[job_name].append(full_path)

    return pp_files, blast_files, align_files, append_files

def fix_pp_order(pp_files, wdl_lines):
    """ make sure that pp_files arrays are in the right order """
    for i, wdl_line in enumerate(wdl_lines):
        wdl_line = wdl_line.strip()
        if wdl_line.startswith('call cactus_preprocess as'):
            # get the name of the called WDL method
            pp_job_name = wdl_line.split()[3]
            assert pp_job_name.startswith('preprocess')
            if pp_job_name in pp_files:
                # we have some directory results, go to the next line and get the inputs in order
                next_line = wdl_lines[i+1]
                p = next_line.find('[')
                q = next_line.find(']')
                assert p > 0 and q > 1
                wdl_array = next_line[p+1:q].split(", ")
                pp_array = pp_files[pp_job_name]
                assert len(wdl_array) == len(pp_array)
                # now use the ordering to fix up our cached files
                fixed_array = []
                for wdl_file in wdl_array:
                    name = os.path.basename(wdl_file).strip('"')
                    cache_name = None
                    for pp_file in pp_array:
                        if os.path.basename(pp_file) == name + '.pp':
                            assert not cache_name
                            cache_name = pp_file
                    assert cache_name
                    fixed_array.append(cache_name)
                pp_files[pp_job_name] = fixed_array
    return pp_files
                    
def resolve_pp_files(pp_files, wdl_lines):
    """ replace all possible references to preprocessed files with their cached counterparts """
    out_lines = []

    # brute force replace rules
    replace_map = {}
    for job_name, cache_array in pp_files.items():
        for i, cache_file in enumerate(cache_array):
            replace_map['{}.out_files[{}]'.format(job_name, i)] = '"{}"'.format(cache_file)

    for wdl_line in wdl_lines:
        fixed_line = wdl_line
        for query, target in replace_map.items():
            fixed_line = fixed_line.replace(query, target)
        out_lines.append(fixed_line)

    return out_lines

def resolve_blast_files(blast_files, wdl_lines):
    """ replace all references to blast output arrays with their cached counterparts """
    out_lines = []

    # brute force replace rules
    replace_map = {}
    for job_name, cache_array in blast_files.items():
        replace_map['{}.out_files'.format(job_name)] = '[{}]'.format(','.join(['"{}"'.format(f) for f in cache_array]))

    for wdl_line in wdl_lines:
        fixed_line = wdl_line
        for query, target in replace_map.items():
            fixed_line = fixed_line.replace(query, target)
        out_lines.append(fixed_line)

    return out_lines

def resolve_align_files(align_files, wdl_lines):
    """ replace all references to hal and ancestor fa files with their cached counterparts """
    out_lines = []

    # brute force replace rules
    replace_map = {}
    for job_name, cache_array in align_files.items():
        assert len(cache_array) == 2
        if cache_array[0].endswith('.hal'):
            cache_array = cache_array[1], cache_array[0]
        assert cache_array[0].endswith('.fa') or cache_array[0].endswith('.pp')
        assert cache_array[1].endswith('.hal')
        replace_map['{}.out_fa_file'.format(job_name)] = '"{}"'.format(cache_array[0])
        replace_map['{}.out_hal_file'.format(job_name)] = '"{}"'.format(cache_array[1])

    for wdl_line in wdl_lines:
        fixed_line = wdl_line
        for query, target in replace_map.items():
            fixed_line = fixed_line.replace(query, target)
        out_lines.append(fixed_line)

    return out_lines

def resolve_append_files(append_files, wdl_lines):
    """ replace all references to hal files with their cached counterparts """
    out_lines = []

    # brute force replace rules
    replace_map = {}
    for job_name, cache_array in append_files.items():
        assert len(cache_array) == 1 and cache_array[0].endswith('.hal')
        cache_file = cache_array[0]
        replace_map['{}.out_file'.format(job_name)] = '"{}"'.format(cache_file)

    for wdl_line in wdl_lines:
        fixed_line = wdl_line
        for query, target in replace_map.items():
            fixed_line = fixed_line.replace(query, target)
        out_lines.append(fixed_line)

    return out_lines
    

def remove_jobs(job_names, wdl_lines):
    """ if a job's in our cache, we remove its call from the wdl """
    out_lines = []

    hide_line = False
    for i, wdl_line in enumerate(wdl_lines):
        wdl_line = wdl_line.strip()
        if wdl_line.startswith('call '):
            assert not hide_line
            job_name = wdl_line.split()[3]
            if job_name in job_names:
                hide_line = True
        if not hide_line:
            out_lines.append(wdl_lines[i])
                                
        if hide_line and wdl_line.startswith('}'):
            hide_line = False

    return out_lines    


def main():

    if not ((len(sys.argv) == 2 and sys.argv[1] == 'scrape-logs') or (len(sys.argv) == 3 and sys.argv[1] == 'resume')):
        sys.stderr.write("{}: Extract logs or edit WDL to use temporary outputs cached in a Google Bucket.\n".format(os.path.basename(sys.argv[0])))
        sys.stderr.write("Usage:\n  To resume:   gsutil ls -r <gs://bucket/prefix> | {} resume aln.wdl > aln.resume.wdl\n".format(os.path.basename(sys.argv[0])))
        sys.stderr.write("  To get logs: gsutil ls -lr <gs://bucket/prefix> | {} scrape-logs\n".format(os.path.basename(sys.argv[0])))
        return 1

    if sys.argv[1] == 'scrape-logs':
        scrape_logs(sys.stdin)
        return 0
    assert sys.argv[1] == 'resume'

    # load the directory tree for the bucket
    pp_files, blast_files, align_files, append_files = load_dirtree(sys.stdin)

    # load the wdl file into list of lines
    wdl_lines = []
    with open(sys.argv[2]) as wdl_file:
        for line in wdl_file:
            wdl_lines.append(line)

    # fix up the pp_files order (very important)
    pp_files = fix_pp_order(pp_files, wdl_lines)

    # update the wdl to resolve the preprocess substitutions
    wdl_lines = resolve_pp_files(pp_files, wdl_lines)
    wdl_lines = remove_jobs(pp_files.keys(), wdl_lines)

    # do the same for blast arrays
    wdl_lines = resolve_blast_files(blast_files, wdl_lines)
    wdl_lines = remove_jobs(blast_files.keys(), wdl_lines)

    # and the hal/fa align outputs
    wdl_lines = resolve_align_files(align_files, wdl_lines)
    wdl_lines = remove_jobs(align_files.keys(), wdl_lines)

    # and the hal files from hal_append_cactus_subtree
    wdl_lines = resolve_append_files(append_files, wdl_lines)
    wdl_lines = remove_jobs(append_files.keys(), wdl_lines)

    for line in wdl_lines:
        sys.stdout.write(line)    

    return 0

if __name__ == '__main__':
    main()


