/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusEndsTestShared.h"

static bool nestedTest = 0;

void cactusCapTestSetup(CuTest* testCase) {
    if (!nestedTest) {
        cactusEndsTestSharedSetup(testCase->name);
    }
}

void cactusCapTestTeardown(CuTest* testCase) {
    if (!nestedTest) {
        cactusEndsTestSharedTeardown(testCase->name);
    }
}

void testCap_construct(CuTest* testCase) {
    nestedTest = 0;
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, rootCap != NULL);
    cactusCapTestTeardown(testCase);
}

void testCap_getName(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getName(rootCap) != NULL_NAME);
    CuAssertTrue(testCase, end_getInstance(end, cap_getName(rootCap)) == cap_getReverse(rootCap));
    CuAssertTrue(testCase, cap_getName(leaf2Cap) != NULL_NAME);
    CuAssertTrue(testCase, end_getInstance(end, cap_getName(leaf2Cap)) == leaf2Cap);

    cactusCapTestTeardown(testCase);
}

void testCap_getOrientation(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getOrientation(rootCap) == end_getOrientation(cap_getEnd(rootCap)));
    CuAssertTrue(testCase, cap_getOrientation(leaf1Cap) == end_getOrientation(cap_getEnd(leaf1Cap)));
    CuAssertTrue(testCase, cap_getOrientation(leaf2Cap) == end_getOrientation(cap_getEnd(leaf2Cap)));

    CuAssertTrue(testCase, cap_getOrientation(cap_getReverse(rootCap)) == end_getOrientation(end_getReverse(cap_getEnd(rootCap))));
    CuAssertTrue(testCase, cap_getOrientation(cap_getReverse(leaf1Cap)) == end_getOrientation(end_getReverse(cap_getEnd(leaf1Cap))));
    CuAssertTrue(testCase, cap_getOrientation(cap_getReverse(leaf2Cap)) == end_getOrientation(end_getReverse(cap_getEnd(leaf2Cap))));

    CuAssertTrue(testCase, cap_getOrientation(leaf1Cap) == cap_getOrientation(rootCap));
    CuAssertTrue(testCase, cap_getOrientation(leaf1Cap) != cap_getOrientation(leaf2Cap));

    cactusCapTestTeardown(testCase);
}

void testCap_getReverse(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getReverse(rootCap) != NULL);
    CuAssertTrue(testCase, cap_getReverse(cap_getReverse(rootCap)) == rootCap);
    cactusCapTestTeardown(testCase);
}

void testCap_getEvent(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getEvent(rootCap) == rootEvent);
    CuAssertTrue(testCase, cap_getEvent(cap_getReverse(rootCap)) == rootEvent);
    cactusCapTestTeardown(testCase);
}

void testCap_getEnd(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getEnd(rootCap) == end_getReverse(end));
    CuAssertTrue(testCase, cap_getEnd(cap_getReverse(rootCap)) == end);
    CuAssertTrue(testCase, cap_getEnd(leaf2Cap) == end);
    CuAssertTrue(testCase, cap_getEnd(cap_getReverse(leaf2Cap)) == end_getReverse(end));
    cactusCapTestTeardown(testCase);
}

void testCap_getSegment(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    Block *block = block_construct(2, flower);
    Segment *segment = segment_construct(block, rootEvent);
    CuAssertTrue(testCase, cap_getSegment(segment_get5Cap(segment)) == segment);
    CuAssertTrue(testCase, cap_getSegment(segment_get3Cap(segment)) == segment);
    CuAssertTrue(testCase, cap_getOrientation(segment_get5Cap(segment)) == segment_getOrientation(segment));
    CuAssertTrue(testCase, cap_getOrientation(segment_get3Cap(segment)) == segment_getOrientation(segment));
    CuAssertTrue(testCase, cap_getSegment(cap_getReverse(segment_get5Cap(segment))) == segment_getReverse(segment));
    CuAssertTrue(testCase, cap_getSegment(cap_getReverse(segment_get3Cap(segment))) == segment_getReverse(segment));
    cactusCapTestTeardown(testCase);
}

void testCap_getOtherSegmentCap(CuTest *testCase) {
    cactusCapTestSetup(testCase);

    Block *block = block_construct(3, flower);
    Segment *segment = segment_construct2(block, 2, 1, sequence);
    Cap *_5Cap = segment_get5Cap(segment);
    Cap *_3Cap = segment_get3Cap(segment);

    CuAssertTrue(testCase, cap_getOtherSegmentCap(leaf1Cap) == NULL);
    CuAssertTrue(testCase, cap_getOtherSegmentCap(cap_getReverse(leaf1Cap)) == NULL);

    CuAssertTrue(testCase, cap_getOtherSegmentCap(_5Cap) == _3Cap);
    CuAssertTrue(testCase, cap_getOtherSegmentCap(_3Cap) == _5Cap);

    CuAssertTrue(testCase, cap_getOtherSegmentCap(cap_getReverse(_5Cap)) == cap_getReverse(_3Cap));
    CuAssertTrue(testCase, cap_getOtherSegmentCap(cap_getReverse(_3Cap)) == cap_getReverse(_5Cap));

    cactusCapTestTeardown(testCase);
}

void testCap_segmentCoordinates(CuTest* testCase) {
    /*
     * Tests the coordinates of an segment and its 5 and 3 prime caps.
     */
    cactusCapTestSetup(testCase);

    Block *block = block_construct(3, flower);
    Segment *segment = segment_construct2(block, 2, 1, sequence);
    Cap *_5Cap = segment_get5Cap(segment);
    Cap *_3Cap = segment_get3Cap(segment);

    CuAssertTrue(testCase, cap_getSide(_5Cap));
    CuAssertTrue(testCase, !cap_getSide(_3Cap));
    CuAssertTrue(testCase, cap_getStrand(_5Cap));
    CuAssertTrue(testCase, cap_getStrand(_3Cap));
    CuAssertIntEquals(testCase, 2, cap_getCoordinate(_5Cap));
    CuAssertIntEquals(testCase, 4, cap_getCoordinate(_3Cap));
    CuAssertTrue(testCase, segment_getStrand(segment));
    CuAssertIntEquals(testCase, 2, segment_getStart(segment));
    CuAssertIntEquals(testCase, 3, segment_getLength(segment));

    CuAssertTrue(testCase, !cap_getSide(cap_getReverse(_5Cap)));
    CuAssertTrue(testCase, cap_getSide(cap_getReverse(_3Cap)));
    CuAssertTrue(testCase, !cap_getStrand(cap_getReverse(_5Cap)));
    CuAssertTrue(testCase, !cap_getStrand(cap_getReverse(_3Cap)));
    CuAssertIntEquals(testCase, 2, cap_getCoordinate(cap_getReverse(_5Cap)));
    CuAssertIntEquals(testCase, 4, cap_getCoordinate(cap_getReverse(_3Cap)));
    CuAssertTrue(testCase, !segment_getStrand(segment_getReverse(segment)));
    CuAssertIntEquals(testCase, 4, segment_getStart(segment_getReverse(segment)));
    CuAssertIntEquals(testCase, 3, segment_getLength(segment_getReverse(segment)));

    cactusCapTestTeardown(testCase);
}

void testCap_segmentCoordinatesReverseStrand(CuTest* testCase) {
    /*
     * Tests the coordinates of an segment and its 5 and 3 prime caps.
     */
    cactusCapTestSetup(testCase);

    Block *block = block_construct(3, flower);
    Segment *segment = segment_construct2(block, 2, 0, sequence);
    Cap *_5Cap = segment_get5Cap(segment);
    Cap *_3Cap = segment_get3Cap(segment);

    CuAssertTrue(testCase, cap_getSide(_5Cap));
    CuAssertTrue(testCase, !cap_getSide(_3Cap));
    CuAssertTrue(testCase, !cap_getStrand(_5Cap));
    CuAssertTrue(testCase, !cap_getStrand(_3Cap));
    CuAssertIntEquals(testCase, 4, cap_getCoordinate(_5Cap));
    CuAssertIntEquals(testCase, 2, cap_getCoordinate(_3Cap));
    CuAssertTrue(testCase, !segment_getStrand(segment));
    CuAssertIntEquals(testCase, 4, segment_getStart(segment));
    CuAssertIntEquals(testCase, 3, segment_getLength(segment));

    CuAssertTrue(testCase, !cap_getSide(cap_getReverse(_5Cap)));
    CuAssertTrue(testCase, cap_getSide(cap_getReverse(_3Cap)));
    CuAssertTrue(testCase, cap_getStrand(cap_getReverse(_5Cap)));
    CuAssertTrue(testCase, cap_getStrand(cap_getReverse(_3Cap)));
    CuAssertIntEquals(testCase, 4, cap_getCoordinate(cap_getReverse(_5Cap)));
    CuAssertIntEquals(testCase, 2, cap_getCoordinate(cap_getReverse(_3Cap)));
    CuAssertTrue(testCase, segment_getStrand(segment_getReverse(segment)));
    CuAssertIntEquals(testCase, 2, segment_getStart(segment_getReverse(segment)));
    CuAssertIntEquals(testCase, 3, segment_getLength(segment_getReverse(segment)));

    cactusCapTestTeardown(testCase);
}

void testCap_getCoordinate(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getCoordinate(rootCap) == INT64_MAX);
    CuAssertTrue(testCase, cap_getCoordinate(cap_getReverse(rootCap)) == INT64_MAX);
    CuAssertTrue(testCase, cap_getCoordinate(leaf1Cap) == 4);
    CuAssertTrue(testCase, cap_getCoordinate(cap_getReverse(leaf1Cap)) == 4);
    cactusCapTestTeardown(testCase);
}

void testCap_setCoordinate(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getCoordinate(rootCap) == INT64_MAX);
    CuAssertTrue(testCase, cap_getStrand(rootCap));
    CuAssertTrue(testCase, cap_getSequence(rootCap) == NULL);
    cap_setCoordinates(rootCap, 5, 0, NULL);
    CuAssertTrue(testCase, cap_getCoordinate(rootCap) == 5);
    CuAssertTrue(testCase, !cap_getStrand(rootCap));
    CuAssertTrue(testCase, cap_getSequence(rootCap) == NULL);
    cap_setCoordinates(rootCap, INT64_MAX, 1, NULL);
    CuAssertTrue(testCase, cap_getCoordinate(rootCap) == INT64_MAX);
    CuAssertTrue(testCase, cap_getStrand(rootCap));
    CuAssertTrue(testCase, cap_getSequence(rootCap) == NULL);
    cactusCapTestTeardown(testCase);
}

void testCap_getStrand(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getStrand(rootCap));
    CuAssertTrue(testCase, !cap_getStrand(cap_getReverse(rootCap)));
    CuAssertTrue(testCase, cap_getStrand(leaf1Cap));
    CuAssertTrue(testCase, !cap_getStrand(cap_getReverse(leaf1Cap)));
    CuAssertTrue(testCase, !cap_getStrand(leaf2Cap));
    CuAssertTrue(testCase, cap_getStrand(cap_getReverse(leaf2Cap)));
    cactusCapTestTeardown(testCase);
}

void testCap_getSide(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, !cap_getSide(rootCap));
    CuAssertTrue(testCase, cap_getSide(cap_getReverse(rootCap)));
    CuAssertTrue(testCase, !cap_getSide(leaf1Cap));
    CuAssertTrue(testCase, cap_getSide(cap_getReverse(leaf1Cap)));
    CuAssertTrue(testCase, cap_getSide(leaf2Cap));
    CuAssertTrue(testCase, !cap_getSide(cap_getReverse(leaf2Cap)));
    cactusCapTestTeardown(testCase);
}

void testCap_getSequence(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getSequence(rootCap) == NULL);
    CuAssertTrue(testCase, cap_getSequence(cap_getReverse(rootCap)) == NULL);
    CuAssertTrue(testCase, cap_getSequence(leaf1Cap) == sequence);
    CuAssertTrue(testCase, cap_getSequence(cap_getReverse(leaf1Cap)) == sequence);
    cactusCapTestTeardown(testCase);
}

void testCap_adjacent(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getAdjacency(leaf1Cap) == NULL);
    CuAssertTrue(testCase, cap_getAdjacency(leaf3Cap) == NULL);
    cap_makeAdjacent(leaf1Cap, leaf3Cap);
    CuAssertTrue(testCase, cap_getAdjacency(leaf1Cap) == cap_getReverse(leaf3Cap));
    CuAssertTrue(testCase, cap_getAdjacency(leaf3Cap) == cap_getReverse(leaf1Cap));
    CuAssertTrue(testCase, cap_getAdjacency(cap_getReverse(leaf1Cap)) == leaf3Cap);
    CuAssertTrue(testCase, cap_getAdjacency(cap_getReverse(leaf3Cap)) == leaf1Cap);
    cactusCapTestTeardown(testCase);
}

void testCap_getTopCap(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    End *end1 = end_construct(0, flower);
    End *end2 = end_construct(0, flower);
    End *end3 = end_construct(0, flower);
    Event *intermediateEvent = event_construct4(NULL, 0.0, rootEvent, leafEvent, eventTree);

    Cap *cap1T = cap_construct(end1, rootEvent);
    Cap *cap1I = cap_construct(end1, intermediateEvent);
    Cap *cap1L1 = cap_construct(end1, leafEvent);
    Cap *cap1L2 = cap_construct(end1, leafEvent);
    cap_makeParentAndChild(cap1I, cap1L1);
    cap_makeParentAndChild(cap1I, cap1L2);
    cap_makeParentAndChild(cap1T, cap1I);
    end_setRootInstance(end1, cap1T);
    assert(end_getRootInstance(end1) == cap1T);

    CuAssertTrue(testCase, cap_getTopCap(cap1L1) == NULL);
    CuAssertTrue(testCase, cap_getTopCap(cap_getReverse(cap1L1)) == NULL);
    CuAssertTrue(testCase, cap_getTopCap(cap1L2) == NULL);
    CuAssertTrue(testCase, cap_getTopCap(cap1I) == NULL);

    Cap *cap2T = cap_construct(end2, rootEvent);
    Cap *cap2L = cap_construct(end2, leafEvent);
    cap_makeParentAndChild(cap2T, cap2L);
    end_setRootInstance(end2, cap2T);
    cap_makeAdjacent(cap1L1, cap2L);

    CuAssertTrue(testCase, cap_getTopCap(cap1L1) == cap1T);
    CuAssertTrue(testCase, cap_getTopCap(cap_getReverse(cap1L1)) == cap_getReverse(cap1T));
    CuAssertTrue(testCase, cap_getTopCap(cap1I) == NULL);

    Cap *cap3T = cap_construct(end3, rootEvent);
    Cap *cap3I = cap_construct(end3, intermediateEvent);
    cap_makeParentAndChild(cap3T, cap3I);
    end_setRootInstance(end3, cap3T);
    cap_makeAdjacent(cap1I, cap3I);
    cap_makeAdjacent(cap1T, cap3T);

    CuAssertTrue(testCase, cap_getTopCap(cap1L1) == cap1I);
    CuAssertTrue(testCase, cap_getTopCap(cap_getReverse(cap1L1)) == cap_getReverse(cap1I));
    CuAssertTrue(testCase, cap_getTopCap(cap1I) == cap1T);
    CuAssertTrue(testCase, cap_getTopCap(cap_getReverse(cap1I)) == cap_getReverse(cap1T));

    CuAssertTrue(testCase, cap_getTopCap(cap1T) == NULL);

    cactusCapTestTeardown(testCase);
}

void testCap_getTopFace(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    Face *face = face_construct(flower);
    cap_setTopFace(rootCap, face);
    CuAssertTrue(testCase, cap_getTopFace(rootCap) == face);
    CuAssertTrue(testCase, cap_getTopFace(cap_getReverse(rootCap)) == face);
    cactusCapTestTeardown(testCase);
}

void testCap_getBottomAndTopFaceEnd(CuTest* testCase) {
    cactusCapTestSetup(testCase);

    End *end1 = end_construct(0, flower);
    Cap *cap1T = cap_construct(end1, rootEvent);
    Cap *cap1L = cap_construct(end1, leafEvent);
    cap_makeParentAndChild(cap1T, cap1L);
    end_setRootInstance(end1, cap1T);

    End *end2 = end_construct(0, flower);
    Cap *cap2T = cap_construct(end2, rootEvent);
    Cap *cap2L = cap_construct(end2, leafEvent);
    cap_makeParentAndChild(cap2T, cap2L);
    end_setRootInstance(end2, cap2T);

    CuAssertTrue(testCase, cap_getBottomFaceEnd(cap1T) == NULL);
    CuAssertTrue(testCase, cap_getBottomFaceEnd(cap1L) == NULL);
    CuAssertTrue(testCase, cap_getTopFaceEnd(cap1T) == NULL);
    CuAssertTrue(testCase, cap_getTopFaceEnd(cap1L) == NULL);

    cap_makeAdjacent(cap1L, cap2L);
    cap_makeAdjacent(cap1T, cap2T);

    //Now make the face
    Face *face = face_construct(flower);
    face_allocateSpace(face, 2);
    face_setTopNode(face, 0, cap1T);
    face_setTopNode(face, 1, cap2T);
    face_setBottomNodeNumber(face, 0, 1);
    face_setBottomNodeNumber(face, 1, 1);
    face_addBottomNode(face, 0, cap1L);
    face_addBottomNode(face, 1, cap2L);

    CuAssertTrue(testCase, cap_getTopFaceEnd(cap1L) == NULL);
    CuAssertTrue(testCase, cap_getTopFaceEnd(cap2L) == NULL);
    FaceEnd *faceEnd1 = cap_getTopFaceEnd(cap1T);
    FaceEnd *faceEnd2 = cap_getTopFaceEnd(cap2T);
    CuAssertTrue(testCase, faceEnd_getTopNode(faceEnd1) == cap1T);
    CuAssertTrue(testCase, faceEnd_getTopNode(faceEnd2) == cap2T);

    CuAssertTrue(testCase, cap_getBottomFaceEnd(cap1T) == NULL);
    CuAssertTrue(testCase, cap_getBottomFaceEnd(cap2T) == NULL);
    CuAssertTrue(testCase, cap_getBottomFaceEnd(cap1L) == faceEnd1);
    CuAssertTrue(testCase, cap_getBottomFaceEnd(cap2L) == faceEnd2);

    cactusCapTestTeardown(testCase);
}

void testCap_getParent(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getParent(rootCap) == NULL);
    CuAssertTrue(testCase, cap_getParent(leaf1Cap) == rootCap);
    CuAssertTrue(testCase, cap_getParent(leaf2Cap) == cap_getReverse(rootCap));
    cactusCapTestTeardown(testCase);
}

void testCap_getChildNumber(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertIntEquals(testCase, 3, cap_getChildNumber(rootCap));
    CuAssertIntEquals(testCase, 0, cap_getChildNumber(leaf1Cap));
    CuAssertIntEquals(testCase, 0, cap_getChildNumber(leaf2Cap));
    CuAssertIntEquals(testCase, 0, cap_getChildNumber(leaf3Cap));
    cactusCapTestTeardown(testCase);
}

void testCap_getChild(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_getChild(rootCap, 0) == leaf1Cap);
    if (!nestedTest) {
        CuAssertTrue(testCase, cap_getChild(rootCap, 1) == cap_getReverse(leaf2Cap));
    } else {
        // leaf2Cap is at the end of the child list when it's been
        // serialized, deleted, and unserialized.
        CuAssertTrue(testCase, cap_getChild(rootCap, 2) == cap_getReverse(leaf2Cap));
    }
    cactusCapTestTeardown(testCase);
}

void testCap_isInternal(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    CuAssertTrue(testCase, cap_isInternal(rootCap));
    CuAssertTrue(testCase, !cap_isInternal(leaf1Cap));
    CuAssertTrue(testCase, !cap_isInternal(leaf2Cap));
    cactusCapTestTeardown(testCase);
}

void testCap_serialisation(CuTest* testCase) {
    cactusCapTestSetup(testCase);
    int64_t i;
    void
            *vA =
                    binaryRepresentation_makeBinaryRepresentation(leaf2Cap,
                            (void(*)(void *, void(*)(const void *, size_t,
                                    size_t))) cap_writeBinaryRepresentation, &i);
    CuAssertTrue(testCase, i > 0);
    cap_destruct(leaf2Cap);
    void *vA2 = vA;
    leaf2Cap = cap_loadFromBinaryRepresentation(&vA2, end);
    free(vA);
    nestedTest = 1;
    testCap_getName(testCase);
    testCap_getOrientation(testCase);
    testCap_getReverse(testCase);
    testCap_getEvent(testCase);
    testCap_getEnd(testCase);
    testCap_getSegment(testCase);
    testCap_getOtherSegmentCap(testCase);
    testCap_segmentCoordinates(testCase);
    testCap_segmentCoordinatesReverseStrand(testCase);
    testCap_getCoordinate(testCase);
    testCap_setCoordinate(testCase);
    testCap_getStrand(testCase);
    testCap_getSide(testCase);
    testCap_getSequence(testCase);
    testCap_adjacent(testCase);
    testCap_getTopCap(testCase);
    testCap_getTopFace(testCase);
    testCap_getBottomAndTopFaceEnd(testCase);
    testCap_getParent(testCase);
    testCap_getChildNumber(testCase);
    testCap_getChild(testCase);
    testCap_isInternal(testCase);
    nestedTest = 0;
    cactusCapTestTeardown(testCase);
}

CuSuite* cactusCapTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testCap_getName);
    SUITE_ADD_TEST(suite, testCap_getOrientation);
    SUITE_ADD_TEST(suite, testCap_getReverse);
    SUITE_ADD_TEST(suite, testCap_getEvent);
    SUITE_ADD_TEST(suite, testCap_getEnd);
    SUITE_ADD_TEST(suite, testCap_getSegment);
    SUITE_ADD_TEST(suite, testCap_getOtherSegmentCap);
    SUITE_ADD_TEST(suite, testCap_segmentCoordinates);
    SUITE_ADD_TEST(suite, testCap_segmentCoordinatesReverseStrand);
    SUITE_ADD_TEST(suite, testCap_getCoordinate);
    SUITE_ADD_TEST(suite, testCap_setCoordinate);
    SUITE_ADD_TEST(suite, testCap_getStrand);
    SUITE_ADD_TEST(suite, testCap_getSide);
    SUITE_ADD_TEST(suite, testCap_getSequence);
    SUITE_ADD_TEST(suite, testCap_adjacent);
    SUITE_ADD_TEST(suite, testCap_getTopCap);
    SUITE_ADD_TEST(suite, testCap_getTopFace);
    SUITE_ADD_TEST(suite, testCap_getBottomAndTopFaceEnd);
    SUITE_ADD_TEST(suite, testCap_getParent);
    SUITE_ADD_TEST(suite, testCap_getChildNumber);
    SUITE_ADD_TEST(suite, testCap_getChild);
    SUITE_ADD_TEST(suite, testCap_isInternal);
    SUITE_ADD_TEST(suite, testCap_serialisation);
    SUITE_ADD_TEST(suite, testCap_construct);
    return suite;
}
