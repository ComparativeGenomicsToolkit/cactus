README for cactus


Installation:

Run the the doc/README.txt to get the main installation of the package done.
You will in addition need to install lastz (a pairwise alignment program), and put it on your path.

lastz (http://www.bx.psu.edu/miller_lab/)

Running the pipeline:

This is all subject to change.

Currently you run 'cactus_workflow.py'

...

cactus_workflow.py --help
Usage: cactus_workflow.py [options] contigFilexN

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  --logInfo             Turn on logging at INFO level
  --logDebug            Turn on logging at DEBUG level
  --logLevel=LOGLEVEL   Log at level (may be either INFO/DEBUG/CRITICAL)
  --tempDirRoot=TEMPDIRROOT
                        Path to where temporary directory containing all temp
                        files are created, by default uses the current working
                        directory as the base
  --logFile=LOGFILE     File to log in
  --noRotatingLogging   Turn off rotating logging, which prevents log files
                        getting too big
  --job=JOBFILE         Job file containing command to run
  --speciesTree=SPECIESTREE
                        The species tree relating the input sequences
  --reconstructionTree=RECONSTRUCTIONTREE
                        Top level directory that will be created to write the
                        reconstruction tree structure in
  --aligner=ALIGNER     The program to build alignments from (is used as
                        prefix string, and so may contain additional arguments
                        to the program (such as a configuration file)
  --alignmentIterations=ALIGNMENTITERATIONS
                        The number of recursive alignment iterations to emply
  --treeBuilder=TREEBUILDER
                        The program to build trees from
  --adjacencyBuilder=ADJACENCYBUILDER
                        The program to build adjacencies from
                        
The workflow is designed to be run using a cluster, though it can run on a single machine 
fine. However, it must currently be run from within jobTree, the batch management 
system we use for parallelising the code.

This you run the workflow as a command to jobTree.py. For example:

jobTree.py --command "cactus_workflow.py --speciesTree '(human:0.01, chimp:0.01);' human chimp  
--reconstructionTree FOO/reconProblem --job JOB_FILE" --jobTree FOO/jobTree 
--logDebug

For notes on running commands with jobTree see the jobTree documentation.

The remaining details the arguments to the actual workflow command (which forms the command argument to jobTree).
The essential arguments are the list of contig files/ directories containing contig files,
a species tree relating the sequences, and a place to put the reconstruction tree.

For example, say you have two files, chimp1.fa and human1.fa and the tree '(h:0.01, c:0.01);', you 
would run the command:

cactus_workflow.py --speciesTree '(h:0.01, c:0.01);' human1.fa chimp1.fa  --reconstructionTree ./reconstructionTree --job JOB_FILE

This would create a reconstruction tree directory in the dir ./reconstructionTree, using the given
species tree and sequence files.

Notes:

Sequences are associated with the leaves of the species tree by their left to right order with respect to the
the order of leaves in the species tree.

Thus running the above command associates the sequences in human1.fa with the leaf 'h' in the species tree.
Similarly it associates the sequences in chimp1.fa with the leaf 'c' in the species tree.

However, running:

cactus_workflow.py --speciesTree.py '(h:0.01, c:0.01);' chimp1.fa human1.fa  --reconstructionTree ./reconstructionTree --job JOB_FILE

Would associate the sequences in human1.fa with the leaf 'c' and similarly the sequences in chimp1.fa with 'h'.

If you have more than one sequence file to associate with a leaf of the tree you can pass a directory containing these files
as the positional argument instead. Thus chimp1.fa could be substituted with a chimp1/ directory containing chimp1/chimp1.fa and any
other sequences that are chimp.



