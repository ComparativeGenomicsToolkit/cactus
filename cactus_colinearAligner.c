
void usage() {
	fprintf(stderr, "cactus_colinearAligner [net-names, ordered by order they should be processed], version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-b --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
}

double net_getTotalBaseLength(Net *net) {
	Net_EndIterator *endIterator = net_getEndIterator(net);
	End *end;
	double totalLength = 0.0;
	while((end = net_getNextEnd(endIterator)) != NULL) {
		if(!end_isAtomEnd(end)) {
			End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
			EndInstance *endInstance;
			while((endInstance = end_getNext(instanceIterator)) != NULL) {
				endInstance = endInstance_getStrand(endInstance) ? endInstance : endInstance_getReverse(endInstance);
				if(!endInstance_getSide(endInstance)) {
					EndInstance *endInstance2 = endInstance_getAdjacency(endInstance);
					while(end_isAtomEnd(endInstance_getEnd(endInstance2))) {
						AtomInstance *atomInstance = endInstance_getAtomInstance(endInstance2);
						assert(atomInstance != NULL);
						assert(atomInstance_get5End(atomInstance) == endInstance2);
						endInstance2 = endInstance_getAdjacency(atomInstance_get3End(atomInstance));
						assert(endInstance_getStrand(endInstance2));
						assert(endInstance_getSide(endInstance2));
					}
					assert(endInstance_getStrand(endInstance2));
					assert(endInstance_getSide(endInstance2));
					int32_t length = endInstance_getCoordinate(endInstance2) - endInstance_getCoordinate(endInstance) - 1;
					assert(length >= 0);
					totalLength += length;
				}
			}
			end_destructInstanceIterator(instanceIterator);
		}
	}
	net_destructEndIterator(endIterator);
	return totalLength;
}

struct List *getOrientedSequences(Net *net) {

}

struct List *runPecan(struct List *sequences) {

}

void convertAlignmentToAtoms(struct List *alignment, int32_t alignmentLength, Net *net) {

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

