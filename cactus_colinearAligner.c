
void usage() {
	fprintf(stderr, "cactus_colinearAligner [net-names, ordered by order they should be processed], version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-b --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
	/*
	 * Parse the options.
	 */
	while(1) {
		static struct option long_options[] = {
				{ "logLevel", required_argument, 0, 'a' },
				{ "netDisk", required_argument, 0, 'c' },
				{ "help", no_argument, 0, 'h' },
				{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		int key = getopt_long(argc, argv, "a:b:h", long_options, &option_index);

		if(key == -1) {
			break;
		}

		switch(key) {
		case 'a':
			logLevelString = stringCopy(optarg);
			break;
		case 'c':
			netDiskName = stringCopy(optarg);
			break;
		case 'h':
			usage();
			return 0;
		default:
			usage();
			return 1;
		}
	}

	/*
	 * Load the netdisk
	 */


	/*
	 * Get out the sequences.
	 */

	/*
	 * Call Pecan.
	 */

	/*
	 * Read in the alignment.
	 */

	/*
	 * Convert the alignment into atoms.
	 */

	/*
	 * Write the net back to disk.
	 */
}

