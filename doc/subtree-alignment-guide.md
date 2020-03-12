# Guide to running independent sub-tree alignments

Running Cactus on traditional HPC clusters can be problematic due to limited
support in Toil and the idiosyncrasies of clusters.  This guide describes ways
to run a Cactus alignment as independent sub-tree alignments on multiple
machines and then combine the result.

This approach requires two largish memory systems, each with a largish number
of cores.  The exact system configuration is dependent on the genomes being
aligned.  This uses Toil single-machine mode to do the sub-tree alignments.
This mode uses a lot of threads, around 30 per core.  You may need to increase
the max user processes limit (ulimit -u) to avoid run the alignments.

The goal is to split the alignments into two subtrees, each with around half
of the genomes.  The sequences of the ancestors of the sub-trees are then aligned g



                          /-hg38
               /----------|
               |          \-mm10
               |                                /-bosTau8
               |                     /----------|
    -----------|                     |          \-susScr11
               |          /----------|                     /-felCat8
               |          |          |          /----------|
               |          |          |          |          \-canFam3
               |          |          \----------|
               \----------|                     |          /-rAeg
                          |                     \----------|
                          |                                \-mMyo
                          |          /-sorAra2
                          \----------|
                                     \-conCri1




Assuming the former, the typical way to do this is to *try* to find whole
subtrees representing ~half the genomes each and split into three
alignments, two of which can be done in parallel and one representing the
join. The join can be done just by taking the fasta of the ancestors from
the 2 subtree alignments. (It's okay if the join actually aligns a few leaf
genomes as well as the 2 ancestors. For example the join would probably
have at least hg38 in this case). This can be done many times recursively
but if you value your time/sanity it's not worth doing this many times.

It's a much better idea to use the --root option to do the 2 split
alignments rather than creating separate alignment inputs. This way you can
be sure to get decent outgroups for each. Naively creating 2 totally
separate subtree alignments will give terrible results.




- large number of threads
