#include "cactus.h"
#include "sonLib.h"


int32_t writeFlowerSequencesInFile(Flower *flower, const char *tempFile1, int32_t minimumSequenceLength) {
    FILE *fileHandle = NULL;
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    End *end;
    int32_t sequencesWritten = 0;
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
        Cap *cap;
        while ((cap = end_getNext(instanceIterator)) != NULL) {
            cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
            Cap *cap2 = cap_getAdjacency(cap);
            assert(cap2 != NULL);
            assert(cap_getStrand(cap2));

            if (!cap_getSide(cap)) {
                assert(cap_getSide(cap2));
                int32_t length = cap_getCoordinate(cap2) - cap_getCoordinate(cap) - 1;
                assert(length >= 0);
                if (length >= minimumSequenceLength) {
                    Sequence *sequence = cap_getSequence(cap);
                    assert(sequence != NULL);
                    if(fileHandle == NULL) {
                        fileHandle = fopen(tempFile1, "w");
                    }
                    char *string = sequence_getString(sequence, cap_getCoordinate(cap) + 1, length, 1);
                    fprintf(fileHandle, ">%s|1|%i\n%s\n", cactusMisc_nameToStringStatic(cap_getName(cap)), //sequence)),
                            cap_getCoordinate(cap) + 1, string);
                    free(string);
                    sequencesWritten++;
                }
            }
        }
        end_destructInstanceIterator(instanceIterator);
    }
    flower_destructEndIterator(endIterator);
    if(fileHandle != NULL) {
        fclose(fileHandle);
    }
    return sequencesWritten;
}

void printSequence(FILE *fH, Sequence *sequence, int64_t start, int64_t length) {
    if(start + length > sequence_getLength(sequence)) {
        length = sequence_getLength(sequence) - start;
    }
    fprintf("%lli|%i\n", sequence_getName(sequence), start);
    char *cA = sequence_getString(sequence, start, length, 1);
    fprintf("%s\n", cA);
    free(cA);
}

FILE *getTempChunkFile(char *tempDirRoot, int64_t *tempFileIndex) {
    char *fileName = stString_print("%s/%lli.fa" % (tempDirRoot, (*tempFileIndex)++));
    fprintf(stdout, "%s\n", fileName);
    FILE *fH = fopen(fileName, "w");
    free(fileName);
    return fH;
}

void makeChunksP(stList *chunkSeqs,
        int64_t chunkSize, int64_t overlapSize,
        char *tempDirRoot, int64_t *tempFileIndex) {
    FILE *fH = getTempChunkFile(rootdir, tempDirIndex);
    int64_t j=0;
    for(int32_t i=0; i<stList_length(chunkSeqs)-1; i++) {
        Sequence *seq = stList_get(chunkSeqs, i);
        printSequence(fH, seq, 0, INT64_MAX);
        j += sequence_getLength(sequence);
    }
    Sequence *seq = stList_pop(chunkSeqs);
    int64_t i = chunkSize-j;
    printSequence(fH, seq, 0, i);
    while(i < sequence_getLength(sequence)) {
        fclose(fH);
        fH = getTempChunkFile(rootdir, tempDirIndex);
        printSequence(fH, seq, i, chunkSize);
        i += chunkSize - overlapSize;
    }
    fclose(fH);
}

void makeChunks(Flower *flower,
        int64_t minimumSequenceLength,
        int64_t chunkSize, int64_t overlapSize,
        char *tempDirRoot, int64_t *tempFileIndex) {
    Sequence *seq;
    Flower_SequenceIterator *seqIt = flower_getSequenceIterator(flower);
    int64_t i = 0, tempDirIndex = 0;
    stList *chunkSeqs = stList_construct();
    while((seq = flower_getNextSequence(seqIt)) != NULL) {
        if(sequence_getLength(minimumSequenceLength) > minimumSequenceLength) {
            i += sequence_getLength(seq);
            stList_append(chunkSeqs, seq);
            if(i > chunkSize) {
                makeChunksP(chunkSeqs, chunkSize, overlapSize, tempDirRoot, tempFileIndex);
                while(stList_length(chunkSeqs) > 0) {
                    stList_pop(chunkSeqs);
                }
                i=0;
            }
        }
    }
    if(stList_length(chunkSeqs) > 0) {
        makeChunksP(chunkSeqs, chunkSize, overlapSize, tempDirRoot, tempFileIndex);
    }
    stList_destruct(chunkSeqs);
    flower_destructSequenceIterator(seqIt);
}

int64_t st_parse64BitInt(const char *cA) {
    int64_t j;
    int32_t i = sscanf(argv[4], "%lli", &j);
        (void)i;
        assert(i == 1);
        return j;
}

int main(int argc, char *argv[]) {
    CactusDisk *cactusDisk;
    Flower *flower;
    st_setLogLevelFromString(argv[1]);
    assert(argc == 8);
    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(argv[2]);
    cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(argv[3]));
    assert(flower != NULL);
    st_logInfo("Read the flower\n");

    //Parse the other inputs
    int64_t minimumSequenceLength = st_parse64BitInt(argv[4]);
    assert(minimumSequenceLength >= 0);
    int64_t chunkSize = st_parse64BitInt(argv[5]);
    assert(chunkSize > 0);
    int64_t overlapSize = st_parse64BitInt(argv[6]);
    assert(overlapSize >= 0);
    char *tempDirRoot = argv[7];
    int64_t tempFileIndex = 0;
    st_logInfo("Going to chunk up the sequences with %lli chunk-size,"
            "%lli overlap-size, %lli minimum-sequence-length and %s temporary directory"
            "for files \n", chunkSize, overlapSize, minimumSequenceLength, tempDirRoot);

    makeChunks(flower, minimumSequenceLength,
           chunkSize, overlapSize,
           tempDirRoot, &tempFileIndex);

    st_logInfo("Finished making chunks\n");
    return 0;
}
