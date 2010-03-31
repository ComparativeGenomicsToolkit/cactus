#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <stdlib.h>
#include <sstream>
#include <iostream>

#include "xmlParser.h"
#include "Argument_helper.h"
#include "reconstructionTree.h"

extern "C" {
	#include "pinchGraph.h"
	#include "cactusGraph.h"
	#include "commonC.h"
	#include "fastCMaths.h"
	#include "bioioC.h"
	#include "hashTableC.h"
	#include "net.h"
	#include "pairwiseAlignment.h"
};

int main(int argc, char *argv[]) {
	/*
	 * The script checks a reconstruction tree.
	 *
	 */

	///////////////////////////////////////////////////////////////////////////
	// (0) Parse the inputs handed by genomeCactus.py / setup stuff.
	///////////////////////////////////////////////////////////////////////////

	std::string logLevelString = "None";
	std::string absolutePathPrefix = "None";
	std::string relativeReconstructionProblemFile = "None";
	bool checkRecursive = false;
	bool checkAdjacencies = true;

	dsr::Argument_helper ah;

	ah.new_named_string('a', "logLevel", "", "Set the log level", logLevelString);
	ah.new_named_string('c', "absolutePathPrefix", "", "The absolute file path to the reconstruction tree hierarchy", absolutePathPrefix);
	ah.new_named_string('d', "reconstructionProblem", "", "The relative path to the file containing the reconstruction problem to solve", relativeReconstructionProblemFile);

	ah.new_flag('e', "recursive", "Check reconstruction tree recursively (default is false)", checkRecursive);
	ah.new_flag('f', "dontCheckAdjacencies", "Don't check the adjacencies/operations of the reconstruction tree", checkAdjacencies);

	ah.set_description("Checks a reconstruction file");
	ah.set_author("Benedict Paten, benedict@soe.ucsc.edu");
	ah.set_version("0.1");

	ah.process(argc, argv);

	std::string absoluteReconstructionProblemFile = (absolutePathPrefix + "/" + relativeReconstructionProblemFile).c_str();

	//////////////////////////////////////////////
	//Set up logging/log inputs
	//////////////////////////////////////////////

	if(strcmp(logLevelString.c_str(), "INFO") == 0) {
		setLogLevel(LOGGING_INFO);
	}
	if(strcmp(logLevelString.c_str(), "DEBUG") == 0) {
		setLogLevel(LOGGING_DEBUG);
	}

	logInfo("Recursion problem file : %s\n", absoluteReconstructionProblemFile.c_str());

	///////////////////////////////////////////////////////////////////////////
	// Check the reconstruction problem
	///////////////////////////////////////////////////////////////////////////

	XMLNode xMainNode=XMLNode::openFileHelper(absoluteReconstructionProblemFile.c_str(), "reconstruction_problem");
	logInfo("Parsed the reconstruction problem from the XML file\n");

	checkReconstructionTree(absolutePathPrefix.c_str(), xMainNode, checkRecursive ? TRUE : FALSE, checkAdjacencies ? TRUE : FALSE);
	logInfo("Checked the input reconstruction problem file\n");

	return 0;
}
