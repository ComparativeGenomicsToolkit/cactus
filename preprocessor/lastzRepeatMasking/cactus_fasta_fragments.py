#!/usr/bin/env python3
"""
Break a fasta file into fragments.

$$$ todo: spread out the fragment starts so that the last fragment ends at the
$$$       .. end of a sequence, if possible

$$$ todo: find runs of N and reset the fragment start position to skip past
$$$       .. such runs

:Author: Bob Harris (rsharris@bx.psu.edu)
"""

from sys    import argv,stdin,stderr,exit
from random import seed as random_seed,shuffle


def usage(s=None):
        message = """fasta_fragments [options] < fasta_file > fasta_file
  Split a fasta file into overlapping fragments.

  options:
    --fragment=<length>  length of each fragment
                         (default is 100)
    --step=<length>      distance between the start of each fragment
                         (default is 50)
    --shuffle[=<seed>]   randomly shuffle the order that fragments are output;
                         this can be very memory intensive, as all fragments
                         are collected in a list before any are output
                         (by default, fragments are output in sequence order)
    --origin=one         output positions are origin-one
                         (surprisingly, this is the default)
    --origin=zero        output positions are origin-zero
    --head=<number>      limit the number of fragments emitted"""

        if (s == None): exit (message)
        else:           exit ("%s\n%s" % (s,message))


def main():

        fragmentLength = 100
        stepLength     = 50
        shuffleEm      = False
        origin         = "one"
        headLimit      = None

        for arg in argv[1:]:
                if ("=" in arg):
                        argVal = arg.split("=",1)[1]

                if (arg.startswith("--fragment=")):
                        fragmentLength = int(argVal)
                elif (arg.startswith("--step=")):
                        stepLength = int(argVal)
                elif (arg == "--shuffle"):
                        shuffleEm = True
                elif (arg.startswith("--shuffle=")):
                        shuffleEm = True
                        random_seed(argVal)
                elif (arg.startswith("--origin=")):
                        origin = argVal
                        if (origin == "0"): origin = "zero"
                        if (origin == "1"): origin = "one"
                        assert (origin in ["zero","one"]), "can't understand %s" % arg
                elif (arg.startswith("--head=")):
                        headLimit = int_with_unit(argVal)
                elif (arg.startswith("--")):
                        usage("can't understand %s" % arg)
                else:
                        usage("can't understand %s" % arg)

        allN = "N" * fragmentLength

        # process the sequences

        if (shuffleEm):
                fragments = []

        fragNum = 0
        for (name,seq) in fasta_sequences(stdin):
                if (headLimit != None) and (fragNum > headLimit): break

                seq = seq.upper()
                for ix in range(0,len(seq),stepLength):
                        end = min(ix + fragmentLength, len(seq))
                        frag = seq[ix:end]
                        if (frag == allN): continue

                        fragNum += 1
                        if (headLimit != None) and (fragNum > headLimit):
                                print("limit of %d emitted fragments reached" % headLimit, file=stderr)
                                break

                        if (origin == "zero"): header = ">%s_%d" % (name,ix)
                        else:                  header = ">%s_%d" % (name,ix+1)
                        if (shuffleEm):
                                fragments += [(header,frag)]
                        else:
                                print(header)
                                print(frag)

        if (shuffleEm):
                shuffle(fragments)
                for (header,frag) in fragments:
                        print(header)
                        print(frag)


# fasta_sequences--
#       Read the fasta sequences from a file

def fasta_sequences(f):
        seqName = None
        seqNucs = None

        for line in f:
                line = line.strip()

                if (line.startswith(">")):
                        if (seqName != None):
                                yield (seqName,"".join(seqNucs))
                        seqName = line[1:].strip().split()[0]
                        seqNucs = []
                elif (seqName == None):
                        assert (False), "first sequence has no header"
                else:
                        seqNucs += [line]

        if (seqName != None):
                yield (seqName,"".join(seqNucs))


# int_with_unit--
#       Parse a string as an integer, allowing unit suffixes

def int_with_unit(s):
        if (s.endswith("K")):
                multiplier = 1000
                s = s[:-1]
        elif (s.endswith("M")):
                multiplier = 1000 * 1000
                s = s[:-1]
        elif (s.endswith("G")):
                multiplier = 1000 * 1000 * 1000
                s = s[:-1]
        else:
                multiplier = 1

        try:               return               int(s)   * multiplier
        except ValueError: return int(math.ceil(float(s) * multiplier))


if __name__ == "__main__": main()
