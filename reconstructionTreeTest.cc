#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <stdlib.h>
#include <sstream>
#include <iostream>

#include "XMLTools.h"
#include "xmlParser.h"
#include "reconstructionTree.h"

extern "C" {
	#include "pinchGraph.h"
	#include "commonC.h"
	#include "fastCMaths.h"
	#include "bioioC.h"
	#include "hashTableC.h"
	#include "net.h"
	#include "pairwiseAlignment.h"
};

/*
 * I've split out the tests of the reconstruction tree from the functions
 * manipulating it because it was just getting to be too much code in one file.
 */

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods to check a reconstruction tree is okay.
//
//I am not claiming this is fully comprehensive.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * String conventions for names and instances of caps and atoms:
 *
 * Basic names are alpha-numeric only
 *
 * Stubs start with $ sign
 *
 * Instances of an element (cap/atom) are denoted by n.m where n is the name of the element
 * (string) and m is the instance (string) of the element. Thus n.m is an instance of n.
 *
 * Caps are the ends of atoms or, at the top level, telomeres. [n represents the right cap of n,
 * similarly n] represents the left side of n.
 *
 * That's it!
 */

void checkName(const char *name) {
	/*
	 * Checks names is valid.
	 */

	int32_t i;
	for(i=0; i<(int32_t)strlen(name); i++) {
		assert(isalpha(name[i]) || isdigit(name[i])); //names must be alphanumeric
	}

}

void checkCapName(const char *capName) {
	/*
	 * Checks that the cap name is valid.
	 */

	int32_t i;
	char *cA;
	assert(strlen(capName) > 1);
	if(capName[0] == '$') { //is a stub
		checkName(capName+1);
	}
	else {
		i = strlen(capName)-1;
		assert(capName[0] == '[' || capName[i] == ']'); //non stub cap names must begin with '[' or end with ']'
		if(capName[i] == ']') {
			cA = stringCopy(capName);
			cA[i] = '\0';
			checkName(cA);
			free(cA);
		}
		else {
			checkName(capName+1);
		}
	}

}

int32_t checkInstanceName(const char *instanceName) {
	/*
	 * Checks that the instance name is valid.
	 */

	int32_t i;
	i = strlen(instanceName)-1;
	assert(strlen(instanceName) >= 3); //an instance name of the form n.m must be of length greater than or equal to 3, because both n and m must be greater than zero length.
	assert(instanceName[i] != '.'); //the instance component m of an instance name n.m must have non zero length.
	while(instanceName[i] != '.') {
		assert(i > 0); //for name of an instance n.m the name n must be of greater than zero length.
		assert(isalpha(instanceName[i]) || isdigit(instanceName[i])); //instance names must be of the form n.m, where m is an alpha numeric string.
		i--;
	}
	return i;

}

void checkCapInstanceName(const char *capInstanceName) {
	/*
	 * Checks that the cap instance name is valid.
	 */

	int32_t i;
	char *cA;

	cA = stringCopy(capInstanceName);
	i = checkInstanceName(capInstanceName);
	cA[i] = '\0';
	checkCapName(cA);
	free(cA);

}

void checkAtomName(const char *atomName) {
	/*
	 * Checks the atom name is valid.
	 */

	checkName(atomName);
}

void checkAtomInstanceName(const char *atomInstanceName) {
	/*
	 * Checks that the atom instance name is valid.
	 */

	int32_t i;
	char *cA;

	cA = stringCopy(atomInstanceName);
	i = checkInstanceName(atomInstanceName);
	cA[i] = '\0';
	checkAtomName(cA);
	free(cA);

}

void checkInstanceNameIsConsistentWithName(const char *name, const char *instanceName) {
	/*
	 * Checks that the instance name is consistent with the name of the elements.
	 */

	char *cA;
	cA = removeInstance(instanceName);
	assert(strcmp(cA, name) == 0); //for an instance name of the form n.m, n must equal the name of the element.
	free(cA);

}

char *checkReconstructionTree_removeCap(const char *name) {
	/*
	 * Returns the name of the element from an element instance string.
	 */

	char *cA;
	assert(name[0] == '[' || name[strlen(name)-1] == ']'); //a non stub cap name must start with [ or end with ].
	if(name[0] == '[') {
		return stringCopy(name+1);
	}
	cA = stringCopy(name);
	cA[strlen(cA)-1] = '\0';
	return cA;

}

void checkReconstructionTree_checkStringIsUnique(struct hashtable *hash, const char *string) {
	/*
	 * Checks that the string has not been seen before (in the hash).
	 */

	char *cA;

	cA = stringCopy(string);
	assert(hashtable_search(hash, cA) == NULL); //a duplicate name has been detected, names must be unique
	hashtable_insert(hash, cA, cA);

}

void getAtomInstanceEndNames(const char *atomInstanceName, char **leftCapInstanceName, char **rightCapInstanceName) {
	/*
	 * Gets the atom instance names for an atom instance, putting them in the left and right cap instance name arguments.
	 */
	char **cA;
	char *cA2;
	char *cA3;

	if(atomInstanceName[0] == '-') { //deal with inverted atom instance
		atomInstanceName += 1;
		cA = leftCapInstanceName;
		leftCapInstanceName = rightCapInstanceName;
		rightCapInstanceName = cA;
	}

	//split instance and name
	cA2 = removeInstance(atomInstanceName);
	cA3 = getInstance(atomInstanceName);

	//make left hand string
	*leftCapInstanceName = (char *)mallocLocal(sizeof(char)*(2 + strlen(atomInstanceName)));
	sprintf(*leftCapInstanceName, "%s].%s", cA2, cA3);

	*rightCapInstanceName = (char *)mallocLocal(sizeof(char)*(2 + strlen(atomInstanceName)));
	sprintf(*rightCapInstanceName, "[%s.%s", cA2, cA3);

	free(cA2);
	free(cA3);
}


////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods to check the trees of a reconstruction problem.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void checkEventTree(XMLNode eventTreeNode, struct hashtable *hash) {

	/*
	 * Checks an event tree.
	 */
	char *cA;
	int_least32_t i;

	cA = stringCopy(eventTreeNode.getAttribute("event"));

	checkName(cA); //check event name is valid
	assert(hashtable_search(hash, cA) == NULL); //event names must be unique.
	hashtable_insert(hash, cA, cA); //put in the hash

	for(i=0; i<eventTreeNode.nChildNode("clade"); i++) {
		checkEventTree(eventTreeNode.getChildNode("clade", i), hash);
	}

}

int32_t checkCapP(XMLNode capTreeNode, struct hashtable *validEvents,
		struct hashtable *instancesSeen,
		const char *name) {

	/*
	 * Checks a cap tree.
	 */
	const char *event;
	const char *instance;
	char *cA;
	int32_t i, j;

	//check event is okay
	event = capTreeNode.getAttribute("event");
	checkName(event); //check event name is valid
	assert(hashtable_search(validEvents, (void *)event) != NULL); //check event name is in the event tree.

	//check instance is known and unique.
	instance = (char *)capTreeNode.getAttribute("instance");
	assert(instance != NULL); //cap tree node must have an instance attribute
	checkCapInstanceName(instance);
	assert(hashtable_search(instancesSeen, (void *)instance) == NULL); //cap tree node instance names must be unique within the its tree.
	hashtable_insert(instancesSeen, stringCopy(instance), stringCopy(event));
	cA = removeInstance(instance);
	assert(strcmp(name, cA) == 0); //for an instance name of the form n.m the name n must the same as the name of the element
	free(cA);

	//now do recursion
	j = 0;
	for(i=0; i<capTreeNode.nChildNode("clade"); i++) {
		j += checkCapP(capTreeNode.getChildNode("clade", i),
					   validEvents, instancesSeen, name);
	}
	if(capTreeNode.getAttribute("augmented") != NULL &&
	   strcmp(capTreeNode.getAttribute("augmented"), "true") == 0) {
		return j;
	}
	else {
		return j+1;
	}

}

int32_t checkCap(XMLNode endNode, struct hashtable *validEvents,
		struct hashtable *instancesSeen, const char *rootEventName) {
	/*
	 * Checks that the cap is okay, looking at the events, instances and adjacencies.
	 */

	int32_t i;

	//check the name of the cap
	checkCapName(endNode.getText());

	assert(endNode.nChildNode("phylogeny") == 1); //check we have just one end tree in a cap.
	XMLNode endTreeNode = endNode.getChildNode("phylogeny");
	assert(endTreeNode.nChildNode("clade") == 1); //check we have a single root clade in a end tree.
	//JM assert(strcmp(endTreeNode.getChildNode("clade").getAttribute("event"), rootEventName) == 0); //check that the root event in an end tree is the 'root' event of the event tree for the reconstruction.

	//check the tree
	i = checkCapP(endTreeNode.getChildNode("clade"), validEvents, instancesSeen,
			endNode.getText());
	return i;

}

void checkEventTreeAndChildEventTree(XMLNode eventTreeNode, XMLNode childEventTreeNode) {

	/*
	 * Checks child event tree against a parent.
	 */
	//for each binary node in the child event tree check that it is the same as the event tree node of the parent present (does not assume the same node ordering, hence the messy code)
	int32_t i, j;

	while(eventTreeNode.nChildNode("clade") == 1) {
		eventTreeNode = eventTreeNode.getChildNode("clade");
	}
	while(childEventTreeNode.nChildNode("clade") == 1) {
		childEventTreeNode = childEventTreeNode.getChildNode("clade");
	}

	assert(strcmp(eventTreeNode.getAttribute("event"), childEventTreeNode.getAttribute("event")) == 0); //check that the root event name of the event tree in a parent reconstruction is the same as the root event of a child reconstruction problem.
	assert(eventTreeNode.nChildNode("clade") == childEventTreeNode.nChildNode("clade")); //checks that the multiplicity of the root event of the event tree in a parent reconstruction is the same as the root event of a child reconstruction problem.
	for(i=0; i<eventTreeNode.nChildNode("clade"); i++) {
		XMLNode eventTreeNode2 = eventTreeNode.getChildNode("clade", i);
		while(eventTreeNode2.nChildNode("clade") == 1) {
			eventTreeNode2 = eventTreeNode2.getChildNode("clade");
		}
		XMLNode childEventTreeNode2;
		for(j=0; j<childEventTreeNode.nChildNode("clade"); j++) {
			childEventTreeNode2 = childEventTreeNode.getChildNode("clade", j);
			while(childEventTreeNode2.nChildNode("clade") == 1) {
				childEventTreeNode2 = childEventTreeNode2.getChildNode("clade");
			}
			if(strcmp(eventTreeNode2.getAttribute("event"), childEventTreeNode2.getAttribute("event")) == 0) {
				break;
			}
		}
		checkEventTreeAndChildEventTree(eventTreeNode2, childEventTreeNode2);
	}

}

///Methods for checking adjacencies and operations.

void constructParentMap(XMLNode capInstanceNode, struct hashtable *parentMap) {
	int32_t i;

	assert(capInstanceNode.isAttributeSet("instance")); //
	for(i=0; i<capInstanceNode.nChildNode("clade"); i++) {
		XMLNode childCapInstanceNode = capInstanceNode.getChildNode("clade", i);
		assert(childCapInstanceNode.isAttributeSet("instance"));
		assert(hashtable_search(parentMap, (void *)childCapInstanceNode.getAttribute("instance")) == NULL);
		hashtable_insert(parentMap,
				stringCopy(childCapInstanceNode.getAttribute("instance")),
				stringCopy(capInstanceNode.getAttribute("instance")));
		constructParentMap(childCapInstanceNode, parentMap);
	}
}

void constructAdjacencyMap(XMLNode capInstanceNode, struct hashtable *adjacencyMap,
		int32_t checkInternalAdjacencies) {
	int32_t i;

	for(i=0; i<capInstanceNode.nChildNode("clade"); i++) {
		XMLNode childCapInstanceNode = capInstanceNode.getChildNode("clade", i);
		constructAdjacencyMap(childCapInstanceNode, adjacencyMap, checkInternalAdjacencies);
	}

	if(capInstanceNode.nChildNode("clade") == 0 || checkInternalAdjacencies == TRUE) {
		assert(capInstanceNode.isAttributeSet("instance")); //check cap instance node has an instance attribute
		assert(capInstanceNode.isAttributeSet("adjacency")); //check cap instance node has an adjacency attribute
		assert(hashtable_search(adjacencyMap, (void *)capInstanceNode.getAttribute("instance")) == NULL); //check that the adjacency is unique.

		hashtable_insert(adjacencyMap,
				stringCopy(capInstanceNode.getAttribute("instance")),
				stringCopy(capInstanceNode.getAttribute("adjacency")));
	}
}

void checkCapTreeAdjacenciesP(XMLNode parentCapInstanceNode, XMLNode capInstanceNode,
		struct hashtable *parentMap, struct hashtable *adjacencyMap, XMLNode operationsNode) {

	/*
	 * Checks the adjacencies of a cap tree node, see the main function for overview of checks.
	 */
	int32_t i, j, k, l;
	const char *adjacency;
	const char *parentAdjacency;
	const char *cA;
	const char *cA2;
	struct List *list;
	struct List *list2;

	//Check children recursively
	for(i=0; i<capInstanceNode.nChildNode("clade"); i++) {
		checkCapTreeAdjacenciesP(capInstanceNode, capInstanceNode.getChildNode("clade", i), parentMap, adjacencyMap, operationsNode);
	}

	//check had adjacency
	assert(capInstanceNode.isAttributeSet("adjacency")); //checks that the cap instance node has an adjacency attribute.
	adjacency = capInstanceNode.getAttribute("adjacency");

	//check adjacency is reciprocated
	assert(hashtable_search(adjacencyMap, (void *)adjacency) != NULL); //checks that the adjacency the cap instance has is linked back to it.
	//JM assert(strcmp((char *)hashtable_search(adjacencyMap, (void *)adjacency), capInstanceNode.getAttribute("instance")) == 0); //checks that the adjacency the cap instance has is linked back to it.

	//check adjacency map is okay
	//JM ?? assert(hashtable_search(parentMap, (void *)adjacency) != NULL);

	//check parent has an an adjacency
	assert(parentCapInstanceNode.isAttributeSet("adjacency")); //checks that the cap instance node has an adjacency attribute.
	parentAdjacency = parentCapInstanceNode.getAttribute("adjacency");

	//check parent adjacency is reciprocated
	assert(hashtable_search(adjacencyMap, (void *)parentAdjacency) != NULL); //checks that the adjacency the cap instance has is linked back to it.
	//JM assert(strcmp((char *)hashtable_search(adjacencyMap, (void *)parentAdjacency), parentCapInstanceNode.getAttribute("instance")) == 0); //checks that the adjacency the cap instance has is linked back to it.

	if(!capInstanceNode.isAttributeSet("operation_index")) { //no operation
		//check does not have second adjacency
		assert(!capInstanceNode.isAttributeSet("adjaceny_2")); //check does not have second adjacency, given that it is not part of an operation.

		//check parent does not have second adjacency
		assert(!parentCapInstanceNode.isAttributeSet("adjaceny_2")); //check does not have second adjacency

		//check parent adjacency is equivalent
		//JM assert(strcmp((char *)hashtable_search(parentMap, (void *)adjacency), parentAdjacency) == 0); //check that parent and child end instances have an equivalent adjacency, given that there is no operation.
	}
	else { //has an op
		//check their are exactly two sets of non-equivalent adjacencies
		list = constructEmptyList(0, NULL);
		listAppend(list, hashtable_search(parentMap, (void *)adjacency));
		if(capInstanceNode.isAttributeSet("adjacency_2")) {
			assert(strcmp(adjacency, capInstanceNode.getAttribute("adjacency_2")) != 0); //check that second adjacency is not the same as the first adjacency.
			listAppend(list, hashtable_search(parentMap, (void *)capInstanceNode.getAttribute("adjacency_2")));
		}

		list2 = constructEmptyList(0, NULL);
		listAppend(list2, (void *)parentAdjacency);
		if(parentCapInstanceNode.isAttributeSet("adjacency_2")) {
			assert(strcmp(parentAdjacency, parentCapInstanceNode.getAttribute("adjacency_2")) != 0); //check that second adjacency is not the same as the first adjacency.
			listAppend(list2, (void *)parentCapInstanceNode.getAttribute("adjacency_2"));
		}

		k = 0;
		l = 0;
		for(i=0; i<list->length; i++) {
			cA = (char *)list->list[i];
			for(j=0; j<list2->length; j++) {
				cA2 = (char *)list2->list[j];
				if(strcmp(cA, cA2) == 0) {
					k++;
				}
				else {
					l++;
				}
			}
		}
		assert(k == l || (k == 0 && l == 1));
		destructList(list);
		destructList(list2);

		//check the operation exists
		assert(sscanf(capInstanceNode.getAttribute("operation_index"), "%i", &i) == 1); //get the operation index.
		assert(i >= 0); //check the operation index is greater than or equal to zero.
		assert(i < operationsNode.nChildNode("operation")); //check operation index is less than total number of operations in the reconstruction problem.
		XMLNode operationNode = operationsNode.getChildNode("operation", i);
		assert(strcmp(operationNode.getAttribute("operation_index"), capInstanceNode.getAttribute("operation_index")) == 0); //check that the ith operation in the operation list has an operation index equal to i.

		//and that the branch lies on the switch
		//point range, or on a branch off the switch point range.

		//if the cap instance is part of the configuration of an operation check
		//it has one adjacency, and that the adjacency is covered by the operation.
		assert(operationNode.nChildNode("configuration") == 2); //check operation has two configurations.
		for(i=0; i<operationNode.nChildNode("configuration"); i++) {
			XMLNode configurationNode = operationNode.getChildNode("configuration", i);
			if(strcmp(capInstanceNode.getAttribute("event"), configurationNode.getAttribute("event")) == 0) {
				assert(!capInstanceNode.isAttributeSet("adjacency_2")); //check that the adjacencies in the end configurations of an operation are unambiguous.
				assert(configurationNode.nChildNode("adjacency_pairs") == 1); //check operation configuration node has only one list of adjacency pairs.
				XMLNode adjacencyPairsNode = configurationNode.getChildNode("adjacency_pairs");
				k = FALSE;
				assert(adjacencyPairsNode.nChildNode("adjacency_pair") > 0); //check that the list of adjacency pairs in a operation configuration is non empty.
				for(j=0; j<adjacencyPairsNode.nChildNode("adjacency_pair"); j++) {
					XMLNode adjacencyPairNode = adjacencyPairsNode.getChildNode("adjacency_pair", j);
					if(strcmp(adjacencyPairNode.getAttribute("left"), capInstanceNode.getAttribute("instance")) == 0) {
						assert(k == FALSE); //checks we have seen only one adjacency pair in an operation configuration that includes a given cap instance end.
						k = TRUE;
						assert(strcmp(adjacency, adjacencyPairNode.getAttribute("right")) == 0); //checks that the adjacency pair in an operation configuration matches that observed in the adjacencies of the end instances.
					}
					if(strcmp(adjacencyPairNode.getAttribute("right"), capInstanceNode.getAttribute("instance")) == 0) {
						assert(k == FALSE); //checks we have seen only one adjacency pair in an operation configuration that includes a given cap instance end.
						k = TRUE;
						assert(strcmp(adjacency, adjacencyPairNode.getAttribute("left")) == 0);  //checks that the adjacency pair in an operation configuration matches that observed in the adjacencies of the end instances.
					}
				}
				assert(k == TRUE); //checks that the an adjacency observed between two ends that are part of an operation is represented in a configuration of the operation.
			}
		}
	}

}

void checkCapTreeAdjacencies(XMLNode capNode, struct hashtable *parentMap, struct hashtable *adjacencyMap, XMLNode operationsNode) {

	int32_t i;

	XMLNode capInstanceNode = capNode.getChildNode("phylogeny").getChildNode("clade");
	for(i=0; i<capInstanceNode.nChildNode("clade"); i++) {
		checkCapTreeAdjacenciesP(capInstanceNode, capInstanceNode.getChildNode("clade", i), parentMap, adjacencyMap, operationsNode);
	}

}

void checkReconstructionTree(const char *absoluteFilePrefix,
		XMLNode xMainNode,
		int32_t checkChildrenRecursively,
		int32_t checkInternalAdjacencies) {

	/*
	 * Checks if the reconstruction tree is okay.
	 */
	int32_t i, j, k, l, m, n, o;
	char *cA;
	char *cA2;
	struct hashtable *capTreesHash;
	struct hashtable *leafEventsHash;
	struct hashtable *validEvents;
	struct hashtable *atomInstanceCoordinates;
	struct hashtable *instancesSeen;
	struct hashtable *instancesSeen_InChildren;
	struct hashtable *adjacencyMap;
	struct hashtable *parentMap;
	char *leftCapInstanceName;
	char *rightCapInstanceName;
	int32_t *iA;
	char *cA3;
	char *cA4;
	char *cA5;
	const char *atomName;
	const char *rootEventName;

	//Do basics
	assert(xMainNode.nChildNode("event_tree") == 1); //check reconstruction problem has just one
	XMLNode eventTreeNode = xMainNode.getChildNode("event_tree");

	assert(xMainNode.nChildNode("caps") == 1); //check reconstruction problem has just one
	XMLNode capsNode = xMainNode.getChildNode("caps");

	assert(xMainNode.nChildNode("atoms") == 1); //check reconstruction problem has just one
	XMLNode atomsNode = xMainNode.getChildNode("atoms");

	assert(xMainNode.nChildNode("strings") == 1); //check reconstruction problem has just one
	XMLNode stringsNode = xMainNode.getChildNode("strings");

	assert(xMainNode.nChildNode("sequences") == 1); //check reconstruction problem has just one
	XMLNode sequencesNode = xMainNode.getChildNode("sequences");

	assert(xMainNode.nChildNode("adjacency_components") == 1); //check reconstruction problem has just one
	XMLNode adjacencyComponentsNode = xMainNode.getChildNode("adjacency_components");

	assert(xMainNode.nChildNode("chains") == 1); //check reconstruction problem has just one
	XMLNode chainsNode = xMainNode.getChildNode("chains");

	assert(xMainNode.nChildNode("operations") == 1); //check reconstruction problem has just one
	XMLNode operationsNode = xMainNode.getChildNode("operations");

	logDebug("Checked the basics\n");

	/*
	 * Event tree
	 */

	//check that there is an event tree and that each name is unique.
	//also use the hash for testing the cap trees have valid events with valid event names.
	validEvents = create_hashtable(LARGE_CHUNK_SIZE,
							hashtable_stringHashKey, hashtable_stringEqualKey,
							free, NULL);

	assert(eventTreeNode.nChildNode("phylogeny") == 1); //check event tree has one phylogeny tag.
	XMLNode eventTreeRootNode = eventTreeNode.getChildNode("phylogeny");
	assert(eventTreeRootNode.nChildNode("clade") == 1); //check that there is only one tree in the event tree.
	checkEventTree(eventTreeRootNode.getChildNode("clade"), validEvents); //checks each event tree event name is unique
	rootEventName = eventTreeRootNode.getChildNode("clade").getAttribute("event");
	//and valid and puts them in a hash

	logDebug("Checked the event tree\n");

	/*
	 * Caps
	 */

	capTreesHash = create_hashtable(LARGE_CHUNK_SIZE,
									hashtable_stringHashKey, hashtable_stringEqualKey,
											   free, free);

	instancesSeen = create_hashtable(LARGE_CHUNK_SIZE,
												hashtable_stringHashKey, hashtable_stringEqualKey,
												free, free);

	//check each cap and atom tree has a tree
	//and check that each node in a cap and atom tree has an event in the event tree.
	//check each cap tree leaf name is a cap instance. check the internal nodes have valid names.

	for(i=0; i<capsNode.nChildNode("cap"); i++) {
		XMLNode capNode = capsNode.getChildNode("cap", i);
		checkCap(capNode, validEvents, instancesSeen, rootEventName);
		hashtable_insert(capTreesHash, stringCopy(capNode.getText()), capNode.getChildNode("phylogeny").createXMLString());
	}

	logDebug("Checked the caps\n");

	/*
	 * Atoms
	 */

	atomInstanceCoordinates = create_hashtable(LARGE_CHUNK_SIZE,
			hashtable_stringHashKey, hashtable_stringEqualKey,
			free, (void (*)(void *))destructIntPair);

	leafEventsHash = createLeafEvents(xMainNode);

	for(i=0; i<atomsNode.nChildNode("atom"); i++) {
		XMLNode atomNode = atomsNode.getChildNode("atom", i);
		atomName = atomNode.getText();

		logDebug("Checking the atom %s\n", atomName);
		checkAtomName(atomName);

		assert(atomNode.nChildNode("cap") == 2); //checks that the atom has two caps.

		XMLNode capNode = atomNode.getChildNode("cap", 0);
		cA = checkReconstructionTree_removeCap(capNode.getText());
		assert(strcmp(cA, atomName) == 0); //check that the cap name of an atom named n is either [n or n]
		free(cA);
		j = checkCap(capNode, validEvents, instancesSeen, rootEventName);
		hashtable_insert(capTreesHash, stringCopy(capNode.getText()), capNode.getChildNode("phylogeny").createXMLString());

		capNode = atomNode.getChildNode("cap", 1);
		cA = checkReconstructionTree_removeCap(capNode.getText());
		assert(strcmp(cA, atomName) == 0); //check that the cap name of an atom named n is either [n or n]
		free(cA);
		assert(j == checkCap(capNode, validEvents, instancesSeen, rootEventName)); //check that both cap trees of an atom have the same number of non augmented instances.
		hashtable_insert(capTreesHash, stringCopy(capNode.getText()), capNode.getChildNode("phylogeny").createXMLString());

		assert(atomNode.nChildNode("atom_instances") == 1); //check atom has one list of atom instances.
		XMLNode atomInstancesNode = atomNode.getChildNode("atom_instances");
		assert(atomInstancesNode.nChildNode("atom_instance") > 0); //check atom has one or more atom instances.
		assert(atomInstancesNode.nChildNode("atom_instance") == j); //check an atom has equal number of atom instances in its two cap trees as in its list of atom instances.
		m = 0;
		for(j=0; j<atomInstancesNode.nChildNode("atom_instance"); j++) {
			XMLNode atomInstanceNode = atomInstancesNode.getChildNode("atom_instance", j);

			checkAtomInstanceName((char *)atomInstanceNode.getText());
			checkInstanceNameIsConsistentWithName((char *)atomNode.getText(), (char *)atomInstanceNode.getText());

			//no duplicates
			assert(hashtable_search(instancesSeen, (char *)atomInstanceNode.getText()) == NULL); //check atom atom instance name is unique.
			hashtable_insert(instancesSeen, stringCopy(atomInstanceNode.getText()), stringCopy(atomInstanceNode.getText()));

			//check that the caps have been seen
			getAtomInstanceEndNames(atomInstanceNode.getText(), &leftCapInstanceName, &rightCapInstanceName);
			//check that caps have been seen
			assert(hashtable_search(instancesSeen, leftCapInstanceName) != NULL); //check that an atom instance in list of atom instances has caps in the atom's cap trees.
			assert(hashtable_search(instancesSeen, rightCapInstanceName) != NULL); //check that an atom instance in list of atom instances has caps in the atom's cap trees.
			//check events are consistent.
			assert(strcmp((char *)hashtable_search(instancesSeen, leftCapInstanceName), atomInstanceNode.getAttribute("event")) == 0); //check that the events of an atom instance and its caps are all equal.
			assert(strcmp((char *)hashtable_search(instancesSeen, rightCapInstanceName), atomInstanceNode.getAttribute("event")) == 0); //check that the events of an atom instance and its caps are all equal.
			free(hashtable_remove(instancesSeen, leftCapInstanceName, TRUE));
			free(hashtable_remove(instancesSeen, rightCapInstanceName, TRUE));

			free(leftCapInstanceName);
			free(rightCapInstanceName);

			if(hashtable_search(leafEventsHash, (void *)atomInstanceNode.getAttribute("event")) != NULL) {
				m++;
				//put the atom instances in a hash for checking against the strings later on.
				assert(sscanf(atomInstanceNode.getChildNode("coordinates").getAttribute("start"), "%i", &k) == 1);
				assert(sscanf(atomNode.getAttribute("length"), "%i", &l) == 1);
				assert(k >= 0); //check the start coordinate of an atom instance is greater than zero.
				assert(l > 0); //check the length of an atom is greater than zero.
				hashtable_insert(atomInstanceCoordinates, stringCopy(atomInstanceNode.getText()), constructIntPair(k, l));

				//we check that each instance is within exactly one string
				l = FALSE;
				for(k=0; k<stringsNode.nChildNode("string"); k++) {
					XMLNode stringNode = stringsNode.getChildNode("string", k);
					if(stringNode.getText() != NULL) {
						cA2 = stringCopy(stringNode.getText());
						cA = strtok(cA2, " ");
						while(cA != NULL) {
							if(strcmp(cA[0] == '-' ? cA+1 : cA, atomInstanceNode.getText()) == 0) {
								assert(l == FALSE); //check that atom instance is contained in only one string of atoms.
								l = TRUE;

								//check strands are matched
								if(atomInstanceNode.getChildNode("coordinates").getAttribute("strand")[0] == '+') {
									assert(cA[0] != '-'); //check that strand of atom instance is consistent between the atom strings and the actual atom instance.
								}
								else {
									assert(cA[0] == '-'); //check that strand of atom instance is consistent between the atom strings and the actual atom instance.
								}
							}
							cA = strtok(NULL, " ");
						}
						free(cA2);
					}
				}
				assert(l == TRUE); //check that each atom instance is within exactly one string
			}
		}
	assert(m*2 == atomInstancesNode.nChildNode("atom_instance")); //check that the number of atom instances with coordinates tags is the same as the number of atom instances which are leaves.
	}

	logDebug("Checked the atoms\n");

	/*
	 * Check the adjacencies and operations.
	 */

	///////////
	//A cap instance branch has a derived and an ancestral cap instance. (the edge plus the two nodes)
	//An adjacency for the derived cap instance of a cap instance branch is equivalent to the ancestral cap instance's adjacency
	//if the two adjacent cap instances are themselves together a cap instance branch.
	//If either node of a cap branch has two adjacencies (the maximum possible), or if the ancestral and derived cap instance's adjacencies are not
	//equivalent then an operation must be given to explain the possible adjacency change.
	//Furthermore only one possible adjacency exchange is possible per branch thus only two sets of non-equivalent adjacencies are allowed per branch.
	////////////

	//build a map of cap instances and their parents
	parentMap = create_hashtable(LARGE_CHUNK_SIZE,
			 hashtable_stringHashKey, hashtable_stringEqualKey,
			 free, free);
	adjacencyMap = create_hashtable(LARGE_CHUNK_SIZE,
				 hashtable_stringHashKey, hashtable_stringEqualKey,
				 free, free);

	for(i=0; i<capsNode.nChildNode("cap"); i++) {
		XMLNode capInstanceNode = capsNode.getChildNode("cap", i).getChildNode("phylogeny").getChildNode("clade");
		constructParentMap(capInstanceNode, parentMap);
		constructAdjacencyMap(capInstanceNode, adjacencyMap, checkInternalAdjacencies);
	}
	for(i=0; i<atomsNode.nChildNode("atom"); i++) {
		XMLNode atomNode = atomsNode.getChildNode("atom", i);
		for(j=0; j<2; j++) {
			XMLNode capInstanceNode = atomNode.getChildNode("cap", j).getChildNode("phylogeny").getChildNode("clade");
			constructParentMap(capInstanceNode, parentMap);
			constructAdjacencyMap(capInstanceNode, adjacencyMap, checkInternalAdjacencies);
		}
	}
	logDebug("Constructed the internal adjacencies and parent map\n");

	if(checkInternalAdjacencies == TRUE) {
		//now do the actual adjacency checking
		for(i=0; i<capsNode.nChildNode("cap"); i++) {
			checkCapTreeAdjacencies(capsNode.getChildNode("cap", i), parentMap, adjacencyMap, operationsNode);
		}
		for(i=0; i<atomsNode.nChildNode("atom"); i++) {
			XMLNode atomNode = atomsNode.getChildNode("atom", i);
			checkCapTreeAdjacencies(atomNode.getChildNode("cap", 0), parentMap, adjacencyMap, operationsNode);
			checkCapTreeAdjacencies(atomNode.getChildNode("cap", 1), parentMap, adjacencyMap, operationsNode);
		}
	}

	////////////
	//Check the operations:
	//For each operation trace from one end of the switch point range checking the cap instance branches are labelled
	//with the operation and that at the beginning and end the correct adjacencies are created.
	//Check branches off the switch point range contain both adjacencies all the way to the leaves.
	//////////////

	for(i=0; i<operationsNode.nChildNode("operation"); i++) {
		XMLNode operationNode = operationsNode.getChildNode("operation", i);
		assert(operationNode.nChildNode("configuration") == 2); //check each operation has two configuration tags.
		for(j=0; j<operationNode.nChildNode("configuration"); j++) {
			XMLNode configurationNode = operationNode.getChildNode("configuration", j);
			XMLNode adjacencyPairsNode = configurationNode.getChildNode("adjacency_pairs");
			assert(adjacencyPairsNode.nChildNode("adjacency_pair") > 0); //check each operation configuration has at least one pair adjacencies in it.
			for(k=0; k<adjacencyPairsNode.nChildNode("adjacency_pair"); k++) {
				XMLNode adjacencyPairNode = adjacencyPairsNode.getChildNode("adjacency_pair", k);
				assert(strcmp(adjacencyPairNode.getAttribute("left"), (char *)hashtable_search(adjacencyMap, (void *)adjacencyPairNode.getAttribute("right"))) == 0); //check each adjacency pair's adjacency is also represented in the adjacencies encoded in the end trees.
				assert(strcmp(adjacencyPairNode.getAttribute("right"), (char *)hashtable_search(adjacencyMap, (void *)adjacencyPairNode.getAttribute("left"))) == 0); //check each adjacency pair's adjacency is also represented in the adjacencies encoded in the end trees.
			}
		}
	}


	logDebug("Checked the adjacencies and operations\n");

	/*
	 * Strings
	 */

	//check strings
	//we check each string has an encompassing sequence and that coordinates are okay.
	//we also check that is does not overlap with any other string
	for(i=0; i<stringsNode.nChildNode("string"); i++) {
		XMLNode stringNode = stringsNode.getChildNode("string", i);

		logDebug("Checking the string: %s %s\n", stringNode.getAttribute("contig"), stringNode.getAttribute("start"));

		//check we've seen the caps
		assert(hashtable_search(instancesSeen, (void *)stringNode.getAttribute("left")) != NULL);
		assert(hashtable_search(instancesSeen, (void *)stringNode.getAttribute("right")) != NULL);

		//we check each string has an encompassing sequence and that the coordinates are okay.
		l = FALSE;
		for(k=0; k<sequencesNode.nChildNode("sequence"); k++) {
			XMLNode sequenceNode = sequencesNode.getChildNode("sequence", k);
			if(strcmp(stringNode.getAttribute("contig"), sequenceNode.getAttribute("contig")) == 0) {
				assert(l == FALSE);
				l = TRUE;

				assert(sscanf(sequenceNode.getAttribute("length"), "%i", &m) == 1);
				assert(sscanf(stringNode.getAttribute("start"), "%i", &n) == 1);
				assert(sscanf(stringNode.getAttribute("length"), "%i", &o) == 1);

				assert(n >= 0); //start of string is greater than zero
				assert(o >= 0); //length is greater than equal to zero
				assert(n + o <= m); //end is before end of sequence
			}
		}
		assert(l == TRUE);

		logDebug("Checked the string has okay coordinates\n");

		//check does not overlap with any other string
		for(k=0; k<stringsNode.nChildNode("string"); k++) {
			if(k != i) {
				XMLNode stringNode2 = stringsNode.getChildNode("string", k);
				if(strcmp(stringNode.getAttribute("contig"), stringNode2.getAttribute("contig")) == 0) {
					//check do not overlap
					assert(sscanf(stringNode.getAttribute("start"), "%i", &l) == 1);
					assert(sscanf(stringNode.getAttribute("length"), "%i", &m) == 1);
					assert(sscanf(stringNode2.getAttribute("start"), "%i", &n) == 1);
					assert(sscanf(stringNode2.getAttribute("length"), "%i", &o) == 1);

					assert(l != n);
					if(l < n) { //check does not overlap
						assert(l + m <= n);
					}
					else {
						assert(n + o <= l);
					}
				}
			}
		}

		logDebug("Checked the string does not overlap with any other\n");

		if(stringNode.getText() != NULL) {
			cA2 = stringCopy(stringNode.getText());
			cA = strtok(cA2, " ");
			assert(sscanf(stringNode.getAttribute("start"), "%i", &l) == 1);
			assert(sscanf(stringNode.getAttribute("length"), "%i", &m) == 1);
			assert(l >= 0);
			assert(m >= 0);
			n = l;

			logDebug("Checking the string of atom instances is proper\n");

			while(cA != NULL) { //check string of atoms is okay
				iA = (int32_t *)hashtable_search(atomInstanceCoordinates, cA[0] == '-' ? cA+1 : cA);
				assert(iA != NULL);
				assert(iA[0] >= l);
				assert(iA[0] + iA[1] <= l + m);
				assert(iA[0] >= n);
				n = iA[0] + iA[1];
				cA = strtok(NULL, " ");
			}
			free(cA2);
		}

		logDebug("Checked the string of atom instances is proper\n");

		//now check the adjacencies are okay
		cA3 = stringCopy(stringNode.getAttribute("left"));
		assert(hashtable_search(adjacencyMap, cA3) != NULL);

		if(stringNode.getText() != NULL) {
			cA2 = stringCopy(stringNode.getText());
			cA = strtok(cA2, " ");
			do {
				getAtomInstanceEndNames(cA, &cA4, &cA5);
				assert(hashtable_search(adjacencyMap, cA4) != NULL);
				assert(hashtable_search(adjacencyMap, cA5) != NULL);
		//JM		assert(strcmp((char *)hashtable_search(adjacencyMap, cA3), cA4) == 0);
		//JM		assert(strcmp((char *)hashtable_search(adjacencyMap, cA4), cA3) == 0);
				free(cA3);
				free(cA4);
				cA3 = cA5;
				cA = strtok(NULL, " ");
			} while(cA != NULL);
			free(cA2);
		}

		cA4 = (char *)stringNode.getAttribute("right");
		assert(hashtable_search(adjacencyMap, cA4) != NULL);
		//JM assert(strcmp((char *)hashtable_search(adjacencyMap, cA3), cA4) == 0);
		//JM assert(strcmp((char *)hashtable_search(adjacencyMap, cA4), cA3) == 0);
		free(cA3);

		logDebug("Checked the string of atom instances has the same adjacency structure as that given by the atom instances\n");
	}

	logDebug("Checked the strings\n");

	/*
	 * Sequences
	 */

	//check sequences
	//we check that each sequence has one or more referring strings and is unique.
	for(i=0; i<sequencesNode.nChildNode("sequence"); i++) {
		XMLNode sequenceNode = sequencesNode.getChildNode("sequence", i);

		logDebug("Checking the sequence: %s\n", sequenceNode.getAttribute("contig"));

		k = FALSE;
		for(j=0; j<stringsNode.nChildNode("string"); j++) {
			XMLNode stringNode = stringsNode.getChildNode("string", j);
			if(strcmp(sequenceNode.getAttribute("contig"), stringNode.getAttribute("contig")) == 0) {
				k = TRUE;
			}
		}
		assert(k == TRUE);

		for(j=0; j<sequencesNode.nChildNode("sequence"); j++) {
			XMLNode sequenceNode2 = sequencesNode.getChildNode("sequence", j);
			if(j != i) {
				assert(strcmp(sequenceNode.getAttribute("contig"), sequenceNode2.getAttribute("contig")) != 0);
			}
		}

		logDebug("Checked the sequence\n");
	}

	logDebug("Checked the sequences\n");

	/*
	 * Adjacency components
	 */

	//check adjacency components
	//we check that each cap is part of the caps of the reconstruction problem.
	//we check no duplicates
	//we parse the child file and check that its caps matches the parent
	instancesSeen_InChildren = create_hashtable(LARGE_CHUNK_SIZE,
							 hashtable_stringHashKey, hashtable_stringEqualKey,
							 free, NULL);
	int32_t adjacencyComponentsTotalCapNumber = 0;
	for(i=0; i<adjacencyComponentsNode.nChildNode(); i++) {
		//must be in correct order
		XMLNode adjacencyComponentNode = adjacencyComponentsNode.getChildNode("adjacency_component", i);
		assert(sscanf(adjacencyComponentNode.getAttribute("adjacency_component_index"), "%i", &j) == 1);
		assert(j == i);

		logDebug("Checking the adjacency component: %s\n", adjacencyComponentNode.getText());
		assert(adjacencyComponentNode.getText() != NULL);

		//check adjacency component is correctly linked to a child file that contains a history for all the caps.
		cA = getAbsoluteFilePath(absoluteFilePrefix, adjacencyComponentNode.getAttribute("child_file"));
		XMLNode childNode=XMLNode::openFileHelper(cA, "reconstruction_problem");
		free(cA);

		cA2 = stringCopy(adjacencyComponentNode.getText());
		cA = strtok(cA2, " ");
		assert(cA != NULL);

		l = 0;
		while(cA != NULL) {
			l++;
			adjacencyComponentsTotalCapNumber++;
			XMLNode childCapsNode = childNode.getChildNode("caps");

			//checks that the cap is in the child
			k = FALSE;
			for(j=0; j<childCapsNode.nChildNode("cap"); j++) {
				XMLNode childCap = childCapsNode.getChildNode("cap", j);
				if(strcmp(childCap.getText(), cA) == 0) {
					//not currently checking for equivalence with the parent cap.
					assert(k == FALSE);
					k = TRUE;
				}
			}
			assert(k == TRUE);

			//check that the caps are unique (no duplicates)
			checkReconstructionTree_checkStringIsUnique(instancesSeen_InChildren, cA);
			cA = strtok(NULL, " ");
		}
		//clean up
		free(cA2);

		//checks adjacency component is non empty and has equal or greater number of caps as child problem (maybe greater because of added stubs)
		assert(l > 0);
		assert(l <= childNode.getChildNode("caps").nChildNode("cap"));

		//check the child event tree.
		assert(childNode.nChildNode("event_tree") == 1);
		XMLNode childEventTree = childNode.getChildNode("event_tree");
		checkEventTreeAndChildEventTree(eventTreeRootNode.getChildNode("clade"), childEventTree.getChildNode("phylogeny").getChildNode("clade"));

		logDebug("Checked the adjacency component\n");

	}

	//check that non terminal adjacency components have all there caps contained in child problems and
	//that terminal problems have no adjacency components. -- this test is wrong if we allow the creation of stubs as we go
	/*if(capsNode.nChildNode("cap") + atomsNode.nChildNode("atom")*2 != adjacencyComponentsTotalCapNumber) {
		assert(atomsNode.nChildNode("atom") == 0);
		assert(adjacencyComponentsNode.nChildNode() == 0);
	}*/

	logDebug("Checked the adjacency components\n");

	/*
	 * Chains
	 */

	//check chains
	//we check each chain refers to correct adjacency components, and that ends are okay and form a chain of atoms (internally, not at the ends)
	for(i=0; i<chainsNode.nChildNode("chain"); i++) {
		XMLNode chainNode = chainsNode.getChildNode("chain", i);
		//must be in correct order

		logDebug("Checking the chain: %s\n", chainNode.getAttribute("chain_index"));

		assert(sscanf(chainNode.getAttribute("chain_index"), "%i", &j) == 1);
		assert(j == i);

		//length (number of links) of the chain must be 1 or greater
		//assert(sscanf(chainNode.getAttribute("length"), "%i", &j) == 1);
		//assert(j > 0);

		//must have links and number must be equal to length
		XMLNode adjacencyPairsNode = chainNode.getChildNode("adjacency_pairs");
		//assert(j == linksNode.nChildNode("link"));

		for(j=0; j<adjacencyPairsNode.nChildNode(); j++) {
			XMLNode adjacencyPairNode = adjacencyPairsNode.getChildNode("adjacency_pair", j);
			//must be in correct order
			//assert(sscanf(linkNode.getAttribute("link_index"), "%i", &k) == 1);
			//assert(k == j);

			assert(sscanf(adjacencyPairNode.getAttribute("adjacency_component_index"), "%i", &k) == 1);
			//check left and right caps are in the adjacency component
			XMLNode adjacencyComponentNode = adjacencyComponentsNode.getChildNode("adjacency_component", k);
			cA2 = stringCopy(adjacencyComponentNode.getText());
			cA = strtok(cA2, " ");
			k = FALSE;
			l = FALSE;
			while(cA != NULL) {
				if(strcmp(cA, adjacencyPairNode.getAttribute("left")) == 0) {
					assert(k == FALSE);
					k = TRUE;
				}
				if(strcmp(cA, adjacencyPairNode.getAttribute("right")) == 0) {
					assert(l == FALSE);
					l = TRUE;
				}
				cA = strtok(NULL, " ");
			}
			assert(k == TRUE);
			assert(l == TRUE);
			free(cA2);

			if(j > 0) { //check previous right cap is same atom as next left cap
				XMLNode previousAdjacencyPairNode = adjacencyPairsNode.getChildNode("adjacency_pair", j-1);
				cA3 = checkReconstructionTree_removeCap(previousAdjacencyPairNode.getAttribute("right"));
				cA4 = checkReconstructionTree_removeCap(adjacencyPairNode.getAttribute("left"));
				assert(strcmp(previousAdjacencyPairNode.getAttribute("right"), adjacencyPairNode.getAttribute("left")) != 0); //not the same ends
				assert(strcmp(cA3, cA4) == 0);
				free(cA3);
				free(cA4);
			}
		}
	}

	logDebug("Checked the chains\n");

	/*
	 * Cleaning up.
	 */

	//cleanup
	hashtable_destroy(capTreesHash, TRUE, TRUE);
	hashtable_destroy(leafEventsHash, FALSE, TRUE);
	hashtable_destroy(validEvents, FALSE, TRUE);
	hashtable_destroy(instancesSeen, TRUE, TRUE);
	hashtable_destroy(instancesSeen_InChildren, FALSE, TRUE);
	hashtable_destroy(atomInstanceCoordinates, TRUE, TRUE);
	hashtable_destroy(parentMap, TRUE, TRUE);
	hashtable_destroy(adjacencyMap, TRUE, TRUE);

	logDebug("Cleaned up, now running the recursive jobs\n");

	if(checkChildrenRecursively == TRUE) {
		for(i=0; i<adjacencyComponentsNode.nChildNode("adjacency_component"); i++) {
			XMLNode adjacencyComponentNode = adjacencyComponentsNode.getChildNode("adjacency_component", i);
			
			cA = getAbsoluteFilePath(absoluteFilePrefix, adjacencyComponentNode.getAttribute("child_file"));
			XMLNode childNode=XMLNode::openFileHelper(cA, "reconstruction_problem");
			free(cA);
			
			checkReconstructionTree(absoluteFilePrefix, childNode, checkChildrenRecursively, checkInternalAdjacencies);
		}
	}


	logDebug("Finished checking reconstruction problem\n");


}
