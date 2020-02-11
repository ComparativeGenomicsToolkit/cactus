/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusFacesTestShared.h"

static bool nestedTest = 0;

void cactusFaceTestSetup(CuTest* testCase) {
    if (!nestedTest) {
        cactusFacesTestSharedSetup(testCase);
    }
}

void cactusFaceTestTeardown(CuTest* testCase) {
    if (!nestedTest) {
        cactusFacesTestSharedTeardown(testCase);
    }
}

void testFace_construct(CuTest* testCase) {
    nestedTest = 0;
    cactusFaceTestSetup(testCase);
    CuAssertTrue(testCase, face != NULL);
    cactusFaceTestTeardown(testCase);
}

void testFace_getCardinal(CuTest* testCase) {
    cactusFaceTestSetup(testCase);
    CuAssertTrue(testCase, face_getCardinal(face) == 4);
    cactusFaceTestTeardown(testCase);
}

void testFace_faceEndIterator(CuTest* testCase) {
    cactusFaceTestSetup(testCase);
    Face_FaceEndIterator *iterator = face_getFaceEndIterator(face);
    FaceEnd *faceEnd1 = face_getNextFaceEnd(iterator);
    CuAssertTrue(testCase, faceEnd1 != NULL);
    CuAssertTrue(testCase, faceEnd_getTopNode(faceEnd1) == cap_getPositiveOrientation(topCap1));
    FaceEnd *faceEnd2 = face_getNextFaceEnd(iterator);
    CuAssertTrue(testCase, faceEnd2 != NULL);
    CuAssertTrue(testCase, faceEnd_getTopNode(faceEnd2) == cap_getPositiveOrientation(topCap2));
    FaceEnd *faceEnd3 = face_getNextFaceEnd(iterator);
    CuAssertTrue(testCase, faceEnd3 != NULL);
    CuAssertTrue(testCase, faceEnd_getTopNode(faceEnd3) == cap_getPositiveOrientation(topCap3));
    FaceEnd *faceEnd4 = face_getNextFaceEnd(iterator);
    CuAssertTrue(testCase, faceEnd4 != NULL);
    CuAssertTrue(testCase, faceEnd_getTopNode(faceEnd4) == cap_getPositiveOrientation(topCap4));
    CuAssertTrue(testCase, face_getNextFaceEnd(iterator) == NULL);
    CuAssertTrue(testCase, face_getNextFaceEnd(iterator) == NULL);
    CuAssertTrue(testCase, face_getNextFaceEnd(iterator) == NULL);
    Face_FaceEndIterator *iterator2 = face_copyFaceEndIterator(iterator);
    CuAssertTrue(testCase, face_getPreviousFaceEnd(iterator) == faceEnd4);
    CuAssertTrue(testCase, face_getPreviousFaceEnd(iterator2) == faceEnd4); //test copied iterator
    face_destructFaceEndIterator(iterator2);
    CuAssertTrue(testCase, face_getPreviousFaceEnd(iterator) == faceEnd3);
    CuAssertTrue(testCase, face_getPreviousFaceEnd(iterator) == faceEnd2);
    CuAssertTrue(testCase, face_getPreviousFaceEnd(iterator) == faceEnd1);
    CuAssertTrue(testCase, face_getPreviousFaceEnd(iterator) == NULL);
    CuAssertTrue(testCase, face_getPreviousFaceEnd(iterator) == NULL);
    CuAssertTrue(testCase, face_getNextFaceEnd(iterator) == faceEnd1);
    face_destructFaceEndIterator(iterator);
    cactusFaceTestTeardown(testCase);
}

void testFace_getTopNode(CuTest* testCase) {
    cactusFaceTestSetup(testCase);
    CuAssertTrue(testCase, face_getTopNode(face, 0) == cap_getPositiveOrientation(topCap1));
    CuAssertTrue(testCase, face_getTopNode(face, 1) == cap_getPositiveOrientation(topCap2));
    CuAssertTrue(testCase, face_getTopNode(face, 2) == cap_getPositiveOrientation(topCap3));
    CuAssertTrue(testCase, face_getTopNode(face, 3) == cap_getPositiveOrientation(topCap4));
    cactusFaceTestTeardown(testCase);
}

void testFace_getDerivedDestination(CuTest* testCase) {
    cactusFaceTestSetup(testCase);
    CuAssertTrue(testCase, face_getDerivedDestinationAtIndex(face, 0, 0) == cap_getPositiveOrientation(topCap4));
    CuAssertTrue(testCase, face_getDerivedDestinationAtIndex(face, 1, 0) == cap_getPositiveOrientation(topCap3));
    CuAssertTrue(testCase, face_getDerivedDestinationAtIndex(face, 2, 0) == cap_getPositiveOrientation(topCap2));
    CuAssertTrue(testCase, face_getDerivedDestinationAtIndex(face, 3, 0) == cap_getPositiveOrientation(topCap1));
    CuAssertTrue(testCase, face_getDerivedDestination(face, 0) == cap_getPositiveOrientation(topCap4));
    CuAssertTrue(testCase, face_getDerivedDestination(face, 1) == cap_getPositiveOrientation(topCap3));
    CuAssertTrue(testCase, face_getDerivedDestination(face, 2) == cap_getPositiveOrientation(topCap2));
    CuAssertTrue(testCase, face_getDerivedDestination(face, 3) == cap_getPositiveOrientation(topCap1));
    cactusFaceTestTeardown(testCase);
}

void testFace_getBottomNodeNumber(CuTest* testCase) {
    cactusFaceTestSetup(testCase);
    CuAssertTrue(testCase, face_getBottomNodeNumber(face, 0) == 1);
    CuAssertTrue(testCase, face_getBottomNodeNumber(face, 1) == 1);
    CuAssertTrue(testCase, face_getBottomNodeNumber(face, 2) == 1);
    CuAssertTrue(testCase, face_getBottomNodeNumber(face, 3) == 1);
    cactusFaceTestTeardown(testCase);
}

void testFace_getBottomNode(CuTest* testCase) {
    cactusFaceTestSetup(testCase);
    CuAssertTrue(testCase, face_getBottomNode(face, 0, 0) == cap_getPositiveOrientation(bottomCap1));
    CuAssertTrue(testCase, face_getBottomNode(face, 1, 0) == cap_getPositiveOrientation(bottomCap2));
    CuAssertTrue(testCase, face_getBottomNode(face, 2, 0) == cap_getPositiveOrientation(bottomCap3));
    CuAssertTrue(testCase, face_getBottomNode(face, 3, 0) == cap_getPositiveOrientation(bottomCap4));
    cactusFaceTestTeardown(testCase);
}

CuSuite *cactusFaceTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testFace_getTopNode);
    SUITE_ADD_TEST(suite, testFace_getCardinal);
    SUITE_ADD_TEST(suite, testFace_faceEndIterator);
    SUITE_ADD_TEST(suite, testFace_getDerivedDestination);
    SUITE_ADD_TEST(suite, testFace_getBottomNodeNumber);
    SUITE_ADD_TEST(suite, testFace_getBottomNode);
    SUITE_ADD_TEST(suite, testFace_construct);
    return suite;
}
