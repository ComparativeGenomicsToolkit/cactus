#include <ctype.h>
#include <assert.h>

#include "cactus.h"
#include "sonLib.h"

#include "recursiveFileBuilder.h"

static Cap *getCapForReferenceEvent(End *end, Name referenceEventName) {
    /*
     * Get the cap for a given event.
     */
    End_InstanceIterator *it = end_getInstanceIterator(end);
    Cap *cap;
    while ((cap = end_getNext(it)) != NULL) {
        if (event_getName(cap_getEvent(cap)) == referenceEventName) {
            end_destructInstanceIterator(it);
            return cap;
        }
    }
    end_destructInstanceIterator(it);
    assert(0);
    return NULL;
}

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

static void writeMafBlock(FILE *fileHandle, Segment *referenceSegment,
        bool showOnlySubstitutionsWithRespectToReference) {
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
            if (showOnlySubstitutionsWithRespectToReference) {
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

static char *getMafFileName(Cap *cap, const char *directory) {
    return stString_print("%s/%s.maf", directory,
            cactusMisc_nameToStringStatic(cap_getName(cap)));
}

static void addChildMafToFile(FILE *fileHandle, Cap *cap,
        const char *childDirectory) {
    /*
     * Adds the maf from the child problems to this directory.
     */
    char *mafFileName = getMafFileName(cap, childDirectory);
    FILE *childFileHandle = fopen(mafFileName, "r");
    char c;
    while ((c = getc(childFileHandle)) != EOF) {
        putc(c, fileHandle);
    }
    fclose(childFileHandle);
    free(mafFileName);
}

static void makeMafForAdjacency(RecursiveFileBuilder *recursiveFileBuilder, Cap *cap,
        bool showOnlySubstitutionsWithRespectToReference) {
    /*
     * Iterate along thread and build maf file.
     */

    //Maf file name
    while (1) {
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(cap_getCoordinate(adjacentCap) - cap_getCoordinate(cap) >= 1);
        if (cap_getCoordinate(adjacentCap) - cap_getCoordinate(cap) > 1) {
            recursiveFileBuilder_writeAdjacency(recursiveFileBuilder, cap);
        }
        if ((cap = cap_getOtherSegmentCap(adjacentCap)) == NULL) {
            break;
        }
        recursiveFileBuilder_writeSegment(recursiveFileBuilder, cap_getSegment(adjacentCap),
                writeMafBlock);
    }
}

void makeMAFHeader(Flower *flower, FILE *fileHandle) {
    fprintf(fileHandle, "##maf version=1 scoring=N/A\n");
    char *cA = eventTree_makeNewickString(flower_getEventTree(flower));
    fprintf(fileHandle, "# cactus %s\n\n", cA);
    free(cA);
}

void makeMaf(Flower *flower, const char *referenceEventString,
        const char *childDirectory, const char *parentDirectory,
        bool showOnlySubstitutionsWithRespectToReference,
        const char *outputFile) {
    /*
     * Makes mafs for the given flower. If outputFile != NULL then it makes a single maf
     * in the given file, else it writes a maf for each adjacency in the parent directory.
     * Child mafs are located in the child directory and added into the maf file as they are
     * needed.
     */
    Event *referenceEvent = eventTree_getEventByHeader(
            flower_getEventTree(flower), referenceEventString);
    assert(referenceEvent != NULL);
    End *end;
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    FILE *fileHandle = NULL;
    if (outputFile != NULL) {
        fileHandle = fopen(outputFile, "w");
        makeMAFHeader(flower, fileHandle);
    }
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if (end_isStubEnd(end) && end_isAttached(end)) {
            Cap *cap = getCapForReferenceEvent(end,
                    event_getName(referenceEvent)); //The cap in the reference
            assert(cap != NULL);
            assert(cap_getSequence(cap) != NULL);
            cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
            if (!cap_getSide(cap)) {
                if (outputFile == NULL) {
                    char *mafFileName = getMafFileName(cap, parentDirectory);
                    fileHandle = fopen(mafFileName, "w");
                    free(mafFileName);
                }
                makeMafForAdjacency(fileHandle, cap, referenceEvent,
                        childDirectory, parentDirectory,
                        showOnlySubstitutionsWithRespectToReference);
                if (outputFile == NULL) {
                    fclose(fileHandle);
                }
            }
        }
    }
    if (outputFile != NULL) {
        fclose(fileHandle);
    }
    flower_destructEndIterator(endIt);
}
