#!/usr/bin/env python3
"""Checks headers are all unique.
"""
from sonLib.bioio import fastaRead

def checkUniqueHeaders(inputFile, checkAlphaNumeric=False, checkUCSC=False, checkAssemblyHub=True):
    """Check that headers are unique and meet certain requirements."""
    seen = set()
    for header, seq in fastaRead(inputFile):
        if " " in header or "\t" in header:
            raise RuntimeError("The fasta header '%s' contains spaces or tabs. These characters will cause issues in space-separated formats like MAF, and may not function properly when viewed in a browser. Please remove these characters from the input headers and try again." % header)
        mungedHeader = header.split()[0]
        if checkAlphaNumeric and "".join([ i for i in mungedHeader if str.isalnum(i) ]) != mungedHeader: #Check is only alpha numeric
            raise RuntimeError("We found a non-alpha numeric character in the fasta header, and the config file (checkAlphaNumeric option) demands that all fasta headers be alpha numeric: %s" % header)
        if checkUCSC:
            mungedHeader = mungedHeader.split('.')[-1]
            if "".join([ i for i in mungedHeader if (str.isalnum(i) or i == '_' or i == '-' or i == ':') ]) != mungedHeader:
                raise RuntimeError("We found a non-alpha numeric, '-', ':' or '_' prefix in the fasta header (UCSC Names option), please modify the first word after the '>' and after the last '.' in every fasta header to only contain alpha-numeric, '_', ':' or '-' characters, or consider using a more lenient option like --checkForAssemblyHub. The offending header: %s" % header)
        if checkAssemblyHub:
            if "".join([ i for i in mungedHeader if (str.isalnum(i) or i == '_' or i == '-' or i == ':' or i == ".") ]) != mungedHeader:
                raise RuntimeError("An invalid character was found in the first word of a fasta header. Acceptable characters for headers in an assembly hub include alphanumeric characters plus '_', '-', ':', and '.'. Please modify your headers to eliminate other characters. The offending header: %s" % header)
        if mungedHeader in seen:
            raise RuntimeError("We found a duplicated fasta header, the first word of each fasta header should be unique within each genome, as this is a requirement for the output HAL file or any MAF file subsequently created. Please modify the input fasta file. Offending duplicate header: %s" % header)
        seen.add(mungedHeader)
