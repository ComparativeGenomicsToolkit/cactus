# Updating Cactus alignments
This document provides examples of how to use `cactus-update-prepare`, a wrapper tool of `cactus-prepare` yielding a list of command-lines to update cactus alignments in a step-by-step fashion. 

## EvolverPrimates Example
Before show examples of updating an alignment, let's create an alignment. Let's take evolverPrimates [example](https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverPrimates.txt) to create `evolverMammals.hal`, 

```
cactus jobStore ./evolverMammals.txt ./evolverMammals.hal
```
After running the example above, the original `evolverMammals` alignment tree looks like:
```
((simHuman_chr6:0.144018,(simMouse_chr6:0.084509,simRat_chr6:0.091589)mr:0.271974)Anc1:0.020593,(simCow_chr6:0.18908,simDog_chr6:0.16303)Anc2:0.032898)Anc0;
```

## Adding Genomes Examples


Both examples described in this section attempts to add the `simGorilla` genome into the `evolverMammals` example, but in different places. For demonstration purposes, the `input.txt` file mentioned in both cases is the same. This file is a tab-spaced text file containing the input sequences physical locations as well as the branch lengths within the target genome (an existing node in the  `adding-to-a-node` case or a new node inferred in the `adding-to-a-branch` case). 

This file is formatted as follows:
```
name1 path1 length1
name2 path2 length1
...
nameN pathN lengthN
```

The `input.txt` file mentioned below looks like:
```
simGorilla https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/primates/loci1/simGorilla.chr6 0.075
```

### Adding to a node

This example takes the `simGorilla` genome from the evolverPrimates [example](https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverPrimates.txt) and adds it as a **child** of the `mr` node in the evolverMammals alignment [example](https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverMammals.txt). It updates the evolverMammals alignment  at a cost of **one** Cactus alignment only. The update procedure is a [two-fold](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/updating-alignments.md#adding-to-a-node) step. First, it  recreates the `mr` alignment by re-aligning `simMouse_chr6` and `simRat_chr6` (current `mr`'s children) with `simGorilla`. Then, it uses the `halReplaceGenome` (taken from the [HAL API](https://github.com/ComparativeGenomicsToolkit/hal)) to replace the old `mr` genome with the brand-new one in the given HAL-format file whilst keeping the `mr`'s ancestor (parent) genome sequence intact (without any changes).

For this example, run:
```
cactus-update-prepare node ./evolverMammals.hal ./input.txt --genome mr --outDir ./steps --jobStore ./jobstore
```
where 
- `node` indicate the `adding-to-a-node` action
- `evolverMammals.hal` is the HAL-format file containing the current `evolverMammals` alignment
- `input.txt` is a text tab-spaced file providing a list of genomes to be added to a node (in this example, the `mr` node)
- `mr` is the target genome node receiving new genomes
- `--outDir` and `--jobStore` are options inherited from [`cactus-prepare`](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/src/cactus/progressive/cactus_prepare.py)
    - `steps` is the name of the directory where assemblies and all cactus-call dependencies will be placed
    - `jobstore` is the name of the directory used by [Toil](https://toil.readthedocs.io/en/latest/) (Cactus internal scheduler) to store temporary data of job executions



The output of `cactus-update-prepare` is expected to be as follows:
```
## Preprocessor
cactus-preprocess ./jobstore/2 ./steps/seq_file.in ./steps/seq_file.out --inputNames simGorilla --realTimeLogging --logInfo --retryCount 0

## Alignment

### Round 0
cactus-blast ./jobstore/3 ./steps/seq_file.out ./steps/mr.cigar --root mr --realTimeLogging --logInfo --retryCount 0 --includeRoot
cactus-align ./jobstore/4 ./steps/seq_file.out ./steps/mr.cigar ./steps/mr.hal --root mr --realTimeLogging --logInfo --retryCount 0 --maxCores 2 --includeRoot

## Alignment update
halReplaceGenome --bottomAlignmentFile ./steps/mr.hal --topAlignmentFile ./evolverMammals.hal ./evolverMammals.hal mr --hdf5InMemory

## Aligment validation
halValidate --genome mr ./evolverMammals.hal --hdf5InMemory
```


Afterwards, the `evolverMammals` alignment tree  will contain  `simGorilla` as a child of the `mr` node as follows:
```
((simHuman_chr6:0.144018,(simMouse_chr6:0.084509,simRat_chr6:0.091589,simGorilla:0.075)mr:0.271974)Anc1:0.020593,(simCow_chr6:0.18908,simDog_chr6:0.16303)Anc2:0.032898)Anc0;
```

### Adding to a branch

This example now takes the `simGorilla` genome from the evolverPrimates [example](https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverPrimates.txt) and adds it to the **branch** connecting `Anc1` with `mr` in the evolverMammals alignment. It updates the evolverMammals alignment  at a cost of **two** Cactus alignments.  The update procedure is a [four-fold](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/updating-alignments.md#adding-to-a-branch) step. First, it infers a new ancestor (labelled in this example as `AncGorilla`) between `Anc1` and `mr`, where `Anc1` becomes `AncGorilla`'s parent and `mr` becomes `AncGorilla`'s child. Next, the `AncGorilla` genome sequence is  created by aligning the `AncGorilla`'s children (i.e., `mr` and `simGorilla`) to produce the HAL-format file `AncGorilla.hal`. After that, another Cactus call is triggered to make the alignment among `AncGorilla` and its new siblings (`Anc1`'s children) to yield the `Anc1.hal` file. In the final step, the original file `evolverMammals.hal` is updated using `halAddToBranch` (taken from the [HAL API](https://github.com/ComparativeGenomicsToolkit/hal)). It reads the `AncGorilla.hal` and `Anc1.hal` containing updates to patch `evolverMammals.hal`.

For this example, run:
```
cactus-update-prepare branch --parentGenome Anc1 --childGenome mr ./evolverMammals.hal ./input.txt --cactus-prepare-options '--alignCores 4' --outDir ./steps --jobStore ./jobstore --ancestorName AncGorilla --topBranchLength 0.10
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


The output of `cactus-prepare-update` above will looks like as follows:
```
## Preprocessor
cactus-preprocess ./jobstore/0 ./steps/seq_file.in ./steps/seq_file.out --inputNames simGorilla --realTimeLogging --logInfo --retryCount 0

## Alignment

### Round 0
cactus-blast ./jobstore/1 ./steps/seq_file.out ./steps/AncGorilla.cigar --root AncGorilla --realTimeLogging --logInfo --retryCount 0
cactus-align ./jobstore/2 ./steps/seq_file.out ./steps/AncGorilla.cigar ./steps/AncGorilla.hal --root AncGorilla --realTimeLogging --logInfo --retryCount 0 --maxCores 4
hal2fasta ./steps/AncGorilla.hal AncGorilla --hdf5InMemory > ./steps/AncGorilla.fa

### Round 1
cactus-blast ./jobstore/3 ./steps/seq_file.out ./steps/Anc1.cigar --root Anc1 --realTimeLogging --logInfo --retryCount 0 --includeRoot
cactus-align ./jobstore/4 ./steps/seq_file.out ./steps/Anc1.cigar ./steps/Anc1.hal --root Anc1 --realTimeLogging --logInfo --retryCount 0 --maxCores 4 --includeRoot

## Alignment update
halAddToBranch ./evolverMammals.hal ./steps/AncGorilla.hal ./steps/Anc1.hal Anc1 AncGorilla mr simGorilla 0.1 0.075 --hdf5InMemory

## Alignment validation
halValidate --genome Anc1 ./evolverMammals.hal --hdf5InMemory
halValidate --genome AncGorilla ./evolverMammals.hal --hdf5InMemory
halValidate --genome simHuman_chr6 ./evolverMammals.hal --hdf5InMemory
halValidate --genome simGorilla ./evolverMammals.hal --hdf5InMemory

```


Afterwards, the `evolverMammals` alignment tree will contain `simGorilla` as a child of `AncGorilla` (new node) and as a sibling of `simHuman_chr6`, and it will look like as follows:
```
(((simMouse_chr6:0.084509,simRat_chr6:0.091589)mr:0.271974,(simHuman_chr6:0.044018,simGorilla:0.075)AncGorilla:0.1)Anc1:0.020593,(simCow_chr6:0.18908,simDog_chr6:0.16303)Anc2:0.032898)Anc0;
```