#include <ctype.h>
#include <assert.h>

#include "cactus.h"
#include "sonLib.h"

#include "recursiveThreadBuilder.h"

static void writeMafHeaderLine(FILE *fileHandle, Block *block) {
    /*
     * Write the header for a MAF file.
     */
    if (block_getRootInstance(block) != NULL) {
        /* Get newick tree string with internal labels and no unary events */
        char *newickTreeString = block_makeNewickString(block, 1, 0);
        assert(newickTreeString != NULL);
        fprintf(fileHandle, "a score=%i tree='%s'\n",
                block_getLength(block) * block_getInstanceNumber(block),
                newickTreeString);
        free(newickTreeString);
    } else {
        fprintf(fileHandle, "a score=%i\n",
                block_getLength(block) * block_getInstanceNumber(block));
    }
}

static char *formatSequenceHeader(Sequence *sequence) {
    /*
     * Format the header of a sequence.
     */
    const char *sequenceHeader = sequence_getHeader(sequence);
    if (strlen(sequenceHeader) > 0) {
        char *cA = st_malloc(sizeof(char) * (1 + strlen(sequenceHeader)));
        sscanf(sequenceHeader, "%s", cA);
        return cA;
    } else {
        return cactusMisc_nameToString(sequence_getName(sequence));
    }
}

static void writeMafSequenceLine2(FILE *fileHandle, Segment *segment,
        char *segmentString) {
    Sequence *sequence = segment_getSequence(segment);
    assert(sequence != NULL);
    char *sequenceHeader = formatSequenceHeader(sequence);
    int32_t start;
    if (segment_getStrand(segment)) {
        start = segment_getStart(segment) - sequence_getStart(sequence);
    } else { //start with respect to the start of the reverse complement sequence
        start
                = (sequence_getStart(sequence) + sequence_getLength(sequence)
                        - 1) - segment_getStart(segment);
    }
    int32_t length = segment_getLength(segment);
    char *strand = segment_getStrand(segment) ? "+" : "-";
    int32_t sequenceLength = sequence_getLength(sequence);
    fprintf(fileHandle, "s\t%s\t%i\t%i\t%s\t%i\t%s\n", sequenceHeader, start,
            length, strand, sequenceLength, segmentString);
    free(sequenceHeader);
}

static void writeMafSequenceLine(FILE *fileHandle, Segment *segment) {
    /*
     * Write a maf sequence line for the given segment.
     */
    char *segmentString = segment_getString(segment);
    writeMafSequenceLine2(fileHandle, segment, segmentString);
    free(segmentString);
}

static void writeMafSequenceLineShowingOnlyReferenceSubstitutions(
        FILE *fileHandle, Segment *segment, const char *referenceString) {
    /*
     * Write a maf sequence line for the given segment, but showing only differences from reference sequence.
     */
    char *segmentString = segment_getString(segment);
    for (int32_t i = 0; i < segment_getLength(segment); i++) {
        if (toupper(segmentString[i]) == toupper(referenceString[i])) {
            segmentString[i] = '*';
        }
    }
    writeMafSequenceLine2(fileHandle, segment, segmentString);
    free(segmentString);
}

static bool writeMafBlock_showOnlySubstitutionsWithRespectToReference = 0;

static void writeMafBlock(FILE *fileHandle, Segment *referenceSegment) {
    /*
     * Writes a maf block.
     */
    Block *block = segment_getBlock(referenceSegment);
    writeMafHeaderLine(fileHandle, block);
    writeMafSequenceLine(fileHandle, referenceSegment);
    char *referenceString = segment_getString(referenceSegment);
    Block_InstanceIterator *iterator = block_getInstanceIterator(block);
    Segment *segment;
    while ((segment = block_getNext(iterator)) != NULL) {
        if (segment != referenceSegment && segment_getSequence(segment) != NULL) {
            if (writeMafBlock_showOnlySubstitutionsWithRespectToReference) {
                writeMafSequenceLineShowingOnlyReferenceSubstitutions(
                        fileHandle, segment, referenceString);
            } else {
                writeMafSequenceLine(fileHandle, segment);
            }
        }
    }
    fprintf(fileHandle, "\n");
    block_destructInstanceIterator(iterator);
    free(referenceString);
}

void makeMAFHeader(Flower *flower, FILE *fileHandle) {
    fprintf(fileHandle, "##maf version=1 scoring=N/A\n");
    char *cA = eventTree_makeNewickString(flower_getEventTree(flower));
    fprintf(fileHandle, "# cactus %s\n\n", cA);
    free(cA);
}

void makeMaf(Flower *flower, RecursiveFileBuilder *recursiveFileBuilder, Event *referenceEvent,
        FILE *parentFileHandle,
        bool showOnlySubstitutionsWithRespectToReference, bool hasParent) {
    /*
     * Makes mafs for the given flower. If outputFile != NULL then it makes a single maf
     * in the given file, else it writes a maf for each adjacency in the parent directory.
     * Child mafs are located in the child directory and added into the maf file as they are
     * needed.
     */
    //Cheeky global variable
    writeMafBlock_showOnlySubstitutionsWithRespectToReference
            = showOnlySubstitutionsWithRespectToReference;
    End *end;
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    if(!hasParent) {
        makeMAFHeader(flower, parentFileHandle);
    }
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if (end_isStubEnd(end) && end_isAttached(end)) {
            Cap *cap = end_getCapForEvent(end, event_getName(referenceEvent));
            assert(cap != NULL);
            assert(cap_getSequence(cap) != NULL);
            cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
            if (!cap_getSide(cap)) {
                recursiveFileBuilder_writeThread(recursiveFileBuilder, cap, writeMafBlock, NULL);
            }
        }
    }
    flower_destructEndIterator(endIt);
}
