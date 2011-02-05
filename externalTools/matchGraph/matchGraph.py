#!/usr/bin/env python

#Copyright (C) 2011 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

"""matchGraph.py reports the maximum weight matching for a graph.
The default is to report the max weight matching. Alternatively the user can
request that the max cardinality matching be reported.

The input file format is:

Line 0 contains the number of vertices V and the number of edges E.
Lines 1 to E contain the edges referenced by vertex index i, j with weight w.
  Vertices i and j, as well as weight w are integers and the indices start  
  with 0.

The output file format is identical to the input format except that the
weights are not reported.
"""

import getopt
import mwmatching
import sys

def usage():
    print "matchGraph.py -e inputFile -w outputFile [-c]"
    print "\t-e inputFile is the input filename"
    print "\t-w outputFile is the output filename"
    print ""
    print "Optional:"
    print "\t-c Report the max cardinality matching"   

def main():
    # Process the cmd-line arguments
    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'e:w:ch')
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(-1)

    maxCard = False
    for o, a in optlist:
       if o == '-e':
           inputFile = a
       elif o == '-w':
           outputFile = a
       elif o == '-c':
           maxCard = True
       elif o == '-h':
           usage()
           sys.exit()
       else:   
           assert False, "unhandled option"

    # Must specify input and output filenames
    try:
       inputFile
    except:
       print "Need to specify input filename"
       usage()
       sys.exit()
    try:
       outputFile
    except:
       print "Need to specify output filename"
       usage()
       sys.exit()

    # Read the input graph file
    graphArray = []
    try:
        f = open(inputFile, 'r')
    except IOError:
        print "Cannot open %s" % (inputFile)
        sys.exit()
    lineIdx = 0
    for line in f:
        line = line.rstrip()
        if lineIdx == 0:
            (vertexNum, edgeNum) = line.split()
            vertexNum = int(vertexNum)
            edgeNum = int(edgeNum)
        else:
            (vertexI, vertexJ, weight) = line.split()
            vertexI = int(vertexI)
            vertexJ = int(vertexJ)
            weight = int(weight)
            if weight < 0:
                print "Error: Weights must be >= 0"
                sys.exit(-1)
            graphArray.append((vertexI, vertexJ, weight))
        lineIdx += 1

    # Call the maximum weight matching API function
    matchArray = mwmatching.maxWeightMatching(graphArray, maxCard)

    # Convert raw output
    vertexCount = len(matchArray)
    matchVertexCount = 0
    idx = 0
    matchGraphHash = {}
    for val in matchArray:
        if val != -1:
             key = " ".join([str(num) for num in sorted([idx, val])])
             if not key in matchGraphHash:
                 matchGraphHash[key] = 1
             else:
                 matchGraphHash[key] += 1
             matchVertexCount += 1
        idx += 1

    # Write output to file
    try:
        f = open(outputFile, 'w')
    except IOError:
        print "Cannot open %s" % (outputFile)
        sys.exit()
       
    header = "%d %d" % (vertexNum, matchVertexCount/2)
    f.write(header + "\n")
    for edge in matchGraphHash:
        f.write(edge + "\n")
    f.close()

    return

if __name__ == "__main__":
    main()
