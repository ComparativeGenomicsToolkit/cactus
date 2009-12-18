
void usage() {
	fprintf(stderr, "cactus_colinearAligner [net-names, ordered by order they should be processed], version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-b --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
}

int32_t getOrientedEndInstancesP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return netMisc_nameCompare(endInstance_getName((EndInstance *)o1), endInstance_getName((EndInstance *)o2));
}

struct avl_table *getOrientedEndInstances(Net *net, End *leftEnd, End *rightEnd) {
	/*
	 * Gets a list of end instances, oriented such that each following sequence
	 * can be aligned together to create a co-linear alignment.
	 */
	struct avl_table *endInstances = avl_create(getOrientedEndInstancesP, NULL, NULL);

	End_InstanceIterator *iterator;
	EndInstance *endInstance;

	//check they are both caps, inherited from the parent chain.
	assert(end_isCap(leftEnd));
	assert(end_isCap(rightEnd));

	//Get sequences which stretch from the left end to either
	//(1) the right end
	//(2) or a stub.
	//(3) back to itself.
	iterator = end_getInstanceIterator(leftEnd);
	while((endInstance = end_getNext(iterator)) != NULL) {
		assert(avl_find(endInstances, endInstance) == NULL);

		avl_insert(endInstances, endInstance);
	}
	end_destructInstanceIterator(iterator);

	//Get sequences which stretch from the right end to either
	//(1) or a stub.
	//(2) back to itself.
	iterator = end_getInstanceIterator(rightEnd);
	while((endInstance = end_getNext(iterator)) != NULL) {
		assert(avl_find(endInstances, endInstance) == NULL);

		endInstance = endInstance_getReverse(endInstance); //reverse, so in the correct left-to-right orientation
		if(avl_find(endInstances, endInstance) == NULL) {
			avl_insert(endInstances, endInstance);
		}
		else { //must be a stub.
			assert(end_isStub(endInstance_getEnd(endInstance)));
		}
	}
	end_destructInstanceIterator(iterator);

	//Get end instances which are not linked to either the left or right ends, i.e.
	//those which have stubs at both ends.
	Net_EndIterator *endIterator = net_getEndIterator(net);
	End *end;
	while((end = net_getNextEnd(endIterator)) != NULL) {
		assert(end_getOrientation(end)); //check we only get stuff back positively oriented.
		if(end != end_getPositiveOrientation(leftEnd) && end != end_getPositiveOrientation(rightEnd)) {
			assert(end != leftEnd && end != rightEnd);
			assert(end_isStub(end));

			iterator = end_getInstanceIterator(end);
			while((endInstance = end_getNext(iterator)) != NULL) {
				if(avl_find(endInstances, endInstance) == NULL &&
				   avl_find(endInstances, endInstance) == NULL) {
					assert(end_isStub(endInstance_getEnd(endInstance)));
					assert(end_isStub(endInstance_getEnd(endInstance_getAdjacency(endInstance))));

					//at this point we need to decide which direction the stubs should be oriented.
					avl_insert(endInstances, endInstance);
				}
			}
			end_destructInstanceIterator(iterator);
		}
	}
	net_destructEndIterator(endIterator);

	return endInstances;
}

char **getStrings(struct avl_table *endInstances) {
	/*
	 * Gets the strings associated with the adjacency interval for each end instance.
	 */
	char **strings = malloc(sizeof(void *) * avl_count(endInstances));
	struct avl_traverser *iterator;
	iterator = mallocLocal(sizeof(struct avl_traverser));
	avl_t_init(iterator, items);
	EndInstance *endInstance;
	int32_t i = 0;
	while((endInstance = avl_t_next(iterator)) != NULL) {
		Sequence *sequence = endInstance_getSequence(endInstance);
		assert(sequence != NULL);
		EndInstance *endInstance2 = endInstance_getAdjacency(endInstance);
		int32_t strand = 1;
		if(!endInstance_getSide(endInstance)) {
			assert(endInstance_getSide(endInstance2));
		}
		else {
			assert(!endInstance_getSide(endInstance2));
			strand = 0;
			EndInstance *endInstance3 = endInstance;
			endInstance = endInstance2;
			endInstance2 = endInstance3;
		}

		int32_t length = endInstance_getCoordinate(endInstance2)-endInstance_getCoordinate(endInstance)-1;
		assert(length >= 0);
		strings[i] = sequence_getString(sequence, endInstance_getCoordinate(endInstance)+1, length, strand);

	}
	free(iterator);
	return strings;
}

char **runPecan(const char **sequences, int32_t sequenceNumber, const char *tempDir,
		int32_t *alignmentLength) {
	/*
	 *	Runs Pecan on the set of sequences and returns a multiple alignment of the results, in row, column format.
	 */
	FILE *fileHandle;

	//Put the sequences in temporary files
	int32_t i;
	struct List *tempFiles = constructEmptyList(0, free);
	for(i=0; i<sequenceNumber; i++) {
		char *name = printString("%i.fa", i);
		char *tempFile = pathJoin(tempDir, cA);
		free(name);
		fileHandle = fopen(tempFile, "w");
		fastaWrite(sequences[i], cA, fileHandle);
		fclose(fileHandle);
		listAppend(tempFiles, tempFile);
	}

	//Run pecan.
	outputFile = pathJoin(tempDir, "output.mfa");
	const char *tempFilesString = joinStrings(" ", tempFiles->list, tempFiles->length);
	char *command = printString("java bp.pecan.Pecan -F %s -G %s", tempFilesString, outputFile);
	exitOnFailure(system(command), "Something went wrong running Pecan");
	free(command);

	//Read the alignment
	fileHandle = fopen(outputFile, "r");
	struct List *alignedSequences = constructEmptyList(0, free);
	struct List *alignedSequencesLength = constructEmptyList(0, (void (*)(void *))destructInt);
	struct List *fastaNames = constructEmptyList(0, free);
	fastaRead(fileHandle, alignedSequences, alignedSequencesLength, fastaNames);
	fclose(fileHandle);

	//Convert into char array.
	char **alignment = malloc(sizeof(void *) * sequenceNumber);
	assert(alignedSequences->length == sequenceNumber);
	if(sequenceNumber > 0) {
		*alignmentLength = *((int32_t *)alignedSequences->list[0]);
	}
	else {
		*alignmentLength = 0;
	}
	for(i=0; i<sequenceNumber; i++) {
		assert(*alignmentLength == *((int32_t *)alignedSequences->list[i]));
		alignment[i] = alignedSequences->list[i];
	}

	//Cleanup.
	free(outputFile);
	destructList(tempFiles);
	destructList(fastaNames);
	alignedSequences->destructElement = NULL; //don't clean up the strings!
	destructList(alignedSequences);
	destructList(alignedSequencesLength);

	return alignment;
}

void convertAlignmentToAtoms(struct List *endInstances, const char **sequences, int32_t alignmentLength, Net *net) {
	/*
	 * Convert the alignment to atoms.
	 */

	//Setup coordinates for each sequence from the set of end instances.



	//Set the atom length to 1.

	//For 2nd column of the alignment if previous column has same members, extend
	//atom length by 1 else create atom for previous columns, and update
	//coordinates, and set atom length to 1.

	//For last column create atom if alignment of greater than zero length.
}

int main(int argc, char *argv[]) {

	char * logLevelString = NULL;
	char * netDiskName = NULL;

	/*
	 * Parse the options.
	 */
	while(1) {
		static struct option long_options[] = {
				{ "logLevel", required_argument, 0, 'a' },
				{ "netDisk", required_argument, 0, 'b' },
				{ "netName", required_argument, 0, 'c' },
				{ "help", no_argument, 0, 'h' },
				{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		int key = getopt_long(argc, argv, "a:b:c:h", long_options, &option_index);

		if(key == -1) {
			break;
		}

		switch(key) {
			case 'a':
				logLevelString = stringCopy(optarg);
				break;
			case 'b':
				netDiskName = stringCopy(optarg);
				break;
			case 'c':
				netName = stringCopy(optarg);
				break;
			case 'h':
				usage();
				return 0;
			default:
				usage();
				return 1;
		}
	}

	if(logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
		setLogLevel(LOGGING_INFO);
	}
	if(logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
		setLogLevel(LOGGING_DEBUG);
	}

	/*
	 * Load the netdisk
	 */

	netDisk = netDisk_construct(netDiskName);
	logInfo("Set up the net disk\n");

	/*
	 * For each net.
	 */
	for (j = optind; j < argc; j++) {
		/*
		 * Read the net.
		 */
		const char *netName = argv[j];
		logInfo("Processing the net named: %s", netName);
		net = netDisk_getNet(netDisk, netMisc_stringToName(netName));
		logInfo("Parsed the net to be aligned\n");

		/*
		 * Get out the sequences.
		 */
		struct List *sequences = getOrientedSequences(net);
		logInfo("Got the sequences\n");

		/*
		 * Call Pecan.
		 */
		int32_t alignmentLength;
		struct List *alignment = runPecan(sequences, &alignmentLength);
		logInfo("Got the alignment\n");

		/*
		 * Convert the alignment into atoms.
		 */
		convertAlignmentToAtoms(alignment, net, alignmentLength);
		logInfo("Converted the alignment to atoms\n");

		/*
		 * Cleanup
		 */
		destructList(sequences);
		destructList(alignment);
	}

	/*
	 * Write the net back to disk.
	 */
	startTime = time(NULL);
	netDisk_write(netDisk);
	logInfo("Updated the net on disk in: %i seconds\n", time(NULL) - startTime);
}

