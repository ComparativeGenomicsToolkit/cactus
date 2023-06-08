#!/usr/bin/env python3
"""Checks headers are all unique.
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
from cactus.shared.common import cactus_call

def checkUniqueHeaders(inputFile, outputFile, eventName, checkAlphaNumeric=False, checkUCSC=False, checkAssemblyHub=True):
    """Check that headers are unique and meet certain requirements."""
    seen = set()
    for seq_record in SeqIO.parse(inputFile, 'fasta'):
        header = seq_record.description
        seq = seq_record.seq
        if " " in header or "\t" in header:
            raise RuntimeError("The fasta header '{}' for event '{}' contains spaces or tabs. These characters will cause issues in space-separated formats like MAF, and may not function properly when viewed in a browser. Please remove these characters from the input headers and try again.".format(header, eventName))
        mungedHeader = header.split()[0]
        if mungedHeader.startswith('id=') and mungedHeader.find('|') > 3:
            # ignore unique prefix in this check -- it won't end up in final hal
            mungedHeader = mungedHeader[mungedHeader.find('|') + 1:]            
        if checkAlphaNumeric and "".join([ i for i in mungedHeader if str.isalnum(i) ]) != mungedHeader: #Check is only alpha numeric
            raise RuntimeError("We found a non-alpha numeric character in the fasta header, and the config file (checkAlphaNumeric option) demands that all fasta headers be alpha numeric: '{}' in '{}'".format(header, eventName))
        if checkUCSC:
            mungedHeader = mungedHeader.split('.')[-1]
            if "".join([ i for i in mungedHeader if (str.isalnum(i) or i == '_' or i == '-' or i == ':') ]) != mungedHeader:
                raise RuntimeError("We found a non-alpha numeric, '-', ':' or '_' prefix in the fasta header (UCSC Names option), please modify the first word after the '>' and after the last '.' in every fasta header to only contain alpha-numeric, '_', ':' or '-' characters, or consider using a more lenient option like --checkForAssemblyHub. The offending header: '{}' in '{}'".format(header, eventName))
        if checkAssemblyHub:
            if "".join([ i for i in mungedHeader if (str.isalnum(i) or i == '_' or i == '-' or i == ':' or i == ".") ]) != mungedHeader:
                raise RuntimeError("An invalid character was found in the first word of a fasta header. Acceptable characters for headers in an assembly hub include alphanumeric characters plus '_', '-', ':', and '.'. Please modify your headers to eliminate other characters. The offending header: '{}' in '{}'".format(header, eventName))
        if mungedHeader in seen:
            raise RuntimeError("We found a duplicated fasta header, the first word of each fasta header should be unique within each genome, as this is a requirement for the output HAL file or any MAF file subsequently created. Please modify the input fasta file. Offending duplicate header: '{}' in '{}'".format(header, eventName))
        seen.add(mungedHeader)

        SeqIO.write(seq_record, outputFile, 'fasta')

    
def sanitize_fasta_headers(job, fasta_id_map, pangenome=False, log_stats=True):
    """ input must be map of event -> fasta id"""
    out_fasta_id_map = {}
    for event, fasta_id in fasta_id_map.items():
        out_fasta_id_map[event] = job.addChildJobFn(sanitize_fasta_header, fasta_id, event, pangenome, log_stats,
                                                    disk=fasta_id.size*7).rv()
    return out_fasta_id_map

def sanitize_fasta_header(job, fasta_id, event, pangenome, log_stats):
    """ run the fasta through cactus_sanitizeFastaHeaders and ungzip if necessary.
    This doesn't do the full check above (though it could), but will catch serious errors as well as
    make sure everything has a id=EVENT| prefix

    if pangenome is set to true, it will strip everything up to and including last #
    so HG002#0#chr2 would get changed to just chr2 (and then id=HG002| would be prefixed as above)

    this keeps redundant information out of the sequence names, otherwise it can be duplicated
    in the final output.  it also keeps # symbols out of sequence names which can still be
    and issue for, say, hal2assmblyhub.py
    """
    
    work_dir = job.fileStore.getLocalTempDir()
    in_fa_path = os.path.join(work_dir, '{}.fa'.format(event))
    out_fa_path = os.path.join(work_dir, '{}.sanitized.fa'.format(event))
    job.fileStore.readGlobalFile(fasta_id, in_fa_path)
    with open(in_fa_path, 'rb') as in_fa_file:
        # https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed
        is_gzipped = in_fa_file.read(2) == b'\x1f\x8b'
    pg_opts = ['-p'] if pangenome else []
    if is_gzipped:
        cmd = [['gzip', '-dc', in_fa_path], ['cactus_sanitizeFastaHeaders', '-', event] + pg_opts]
    else:
        cmd = ['cactus_sanitizeFastaHeaders', in_fa_path, event] + pg_opts
    cactus_call(parameters=cmd, outfile=out_fa_path)
    job.fileStore.deleteGlobalFile(fasta_id)

    if log_stats:
        analysisString = cactus_call(parameters=["cactus_analyseAssembly", os.path.basename(out_fa_path)],
                                     work_dir=work_dir, check_output=True)
        job.fileStore.logToMaster("Assembly stats for %s: %s" % (event, analysisString))
    
    return job.fileStore.writeGlobalFile(out_fa_path)
    

    
    
