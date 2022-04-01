# Updating Cactus alignments
This document provides examples of how to use `cactus-update-prepare`, a wrapper tool of `cactus-prepare` yielding a list of command-lines to update cactus alignments in a step-by-step fashion. 
## Adding Genomes Examples

### Adding to a node

This example takes the `simChimp` genome from the evolverPrimates [example](https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverPrimates.txt) and adds it as a **child** of the `mr` node in the evolverMammals alignment [example](https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverMammals.txt). It updates the evolverMammals alignment  at a cost of **one** Cactus alignment only. The update procedure is a [two-fold](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/updating-alignments.md#adding-to-a-node) step. First, it  recreates the `mr` alignment by re-aligning `simMouse_chr6` and `simRat_chr6` (current `mr`'s children) with `simChimp`. Then, it uses the `halReplaceGenome` (taken from the [HAL API](https://github.com/ComparativeGenomicsToolkit/hal)) to replace the old `mr` genome with the brand-new one in the given HAL-format file whilst keeping the `mr`'s ancestor (parent) genome sequence intact (without any changes).

For this example, run:
```
cactus-update-prepare node evolverMammals.hal input.txt --genome mr --outDir steps --jobStore jobstore
```
where 
- `node` indicate the `adding-to-a-node` action
- `evolverMammals.hal` is the HAL-format file containing the current `evolverMammals` alignment
- `input.txt` is a text tab-spaced file providing a list of genomes to be added to a node (in this example, the `mr` node)
- `mr` is the target genome node receiving new genomes
- `--outDir` and `--jobStore` are options inherited from [`cactus-prepare`](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/src/cactus/progressive/cactus_prepare.py)
    - `steps` is the name of the directory where assemblies and all cactus-call dependencies will be placed
    - `jobstore` is the name of the directory used by [Toil](https://toil.readthedocs.io/en/latest/) (Cactus internal scheduler) to store temporary data of job executions

This input file is a tab-spaced text file containing the locations of the input sequences as well as the branch lengths with the target genome. It is formatted as follows:
```
name1 path1 length1
name2 path2 length1
...
nameN pathN lengthN
```


For this example, the `input.txt` file mentioned above looks like:
```
simGorilla https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/primates/loci1/simChimp.chr6 0.075
```

Before the execution of this example, the original tree of the `evolverMammals` alignment looks like:
```
((simHuman_chr6:0.144018,(simMouse_chr6:0.084509,simRat_chr6:0.091589)mr:0.271974)Anc1:0.020593,(simCow_chr6:0.18908,simDog_chr6:0.16303)Anc2:0.032898)Anc0;
```


The tree of`evolverMammals` alignment  after the addition of `simChimp` to the `mr` node will look like:
```
((simHuman_chr6:0.144018,(simMouse_chr6:0.084509,simRat_chr6:0.091589,simGorilla:0.075)mr:0.271974)Anc1:0.020593,(simCow_chr6:0.18908,simDog_chr6:0.16303)Anc2:0.032898)Anc0;
```


### Adding to a branch

This example takes the `simChimp` genome from the evolverPrimates [example](https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverPrimates.txt) and adds it to the **branch** connecting `Anc1` with `simHuman_chr6` in the evolverMammals alignment [example](https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverMammals.txt). It updates the evolverMammals alignment  at a cost of **two** Cactus alignments.  The update procedure is a [four-fold](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/updating-alignments.md#adding-to-a-branch) step. First, it infers a new ancestor (labelled in this example as `ch`) between `Anc1` and `simHuman_chr6`, where `Anc1` becomes `ch`'s parent and `simHuman_chr6` becomes `ch`'s child. Next, the `ch` genome sequence is  created by aligning the `ch`'s children (i.e., `simHuman_chr6` and `simChimp`) to produce the HAL-format file `ch.hal`. After that, another Cactus call is triggered to make the alignment among `ch` and its new siblings (`Anc1`'s children) to yield `Anc1.hal` file. Finally, the  given `evolverMammals.hal` is updated using `halAddToBranch` (taken from the [HAL API](https://github.com/ComparativeGenomicsToolkit/hal)). This tool reads the `ch.hal` and `Anc1.hal` files containing updated alignment information about `Anc1` and `mr` to patch the given `evolverMammals.hal` file.

For this example, run:
```
cactus-update-prepare branch evolverMammals.hal input.txt --parentGenome Anc1 --childGenome simHuman_chr6 --ancestorName cg --topBranchLength 0.10 --outDir steps --jobStore jobstore
```
where:
- `branch` indicate the `adding-to-a-branch` action
- `evolverMammals.hal` is the HAL-format file containing the current `evolverMammals` alignment
- `input.txt` is a text tab-spaced file providing a list of genomes to be added to a node (in this example, the `mr` node)
- `--parentGenome` is the up-level genome name of the target branch
- `--childGenome` is the low-level genome of the target branch
- `--ancestorName` is a desired name for the new ancestor inferred between `parentGenome` and `childGenome`
- `--topBranchLength` is the branch length between the new node being inferred and `parentGenome`
    - The branch length between the new node and `childGenome` is  equal to the difference between the existing branch length  connecting `parentGenome` with `childGenome` and  the given `topBranchLength` if  `topBranchLength` is less than or equal to the existing  `parentGenome`-`childGenome` branch length
    - Otherwise the parameter `--forceBottomBranchLength` must be employed to force a specific branch length value
- `--outDir` and `--jobStore` are options inherited from [`cactus-prepare`](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/src/cactus/progressive/cactus_prepare.py)
    - `steps` is the name of the directory where assemblies and all cactus-call dependencies will be placed
    - `jobstore` is the name of the directory used by [Toil](https://toil.readthedocs.io/en/latest/) (Cactus internal scheduler) to store temporary data of job executions

For this example, the `input.txt` file mentioned above looks like:
```
simGorilla https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/primates/loci1/simChimp.chr6 0.075
```

Before the execution of this example, the original tree of the `evolverMammals` alignment looks like:
```
((simHuman_chr6:0.144018,(simMouse_chr6:0.084509,simRat_chr6:0.091589)mr:0.271974)Anc1:0.020593,(simCow_chr6:0.18908,simDog_chr6:0.16303)Anc2:0.032898)Anc0;
```

The tree of the `evolverMammals` alignment after the addition of `simChimp` to the branch connecting `Anc1` to `simHuman_chr6` will look like:
```
(((simMouse_chr6:0.084509,simRat_chr6:0.091589)mr:0.271974,(simHuman_chr6:0.044018,simGorilla:0.075)cg:0.1)Anc1:0.020593,(simCow_chr6:0.18908,simDog_chr6:0.16303)Anc2:0.032898)Anc0;
```