# Updating Cactus alignments
This document provides examples of how to use `cactus-update-prepare`, a wrapper tool of `cactus-prepare` yielding a list of command-lines to update cactus alignments in a step-by-step fashion. 

## EvolverPrimates Example
Before show examples of updating an alignment, let's create an alignment first using the evolverPrimates [example](https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverPrimates.txt) as follows:

```
cactus jobStore ./evolverMammals.txt ./evolverMammals.hal
```
After running the command line above, the `evolverMammals.hal` file is created containing an alignment which the tree looks like as follows:
```
((simHuman_chr6:0.144018,(simMouse_chr6:0.084509,simRat_chr6:0.091589)mr:0.271974)Anc1:0.020593,(simCow_chr6:0.18908,simDog_chr6:0.16303)Anc2:0.032898)Anc0;
```

## Adding Genomes Examples


The examples described on this page add the `simGorilla` genome into the `evolverMammals` alignment created above. For this purpose, a tab-spaced text file (named  `input.txt`) is created containing the `simGorilla` location and the branch length connecting the `simGorilla` genome with its parent (an existing node in the tree for the  `adding-to-a-node` case or a new node being inferred for the `adding-to-a-branch` case). 

This file is formatted as follows:
```
name1 path1 length1
name2 path2 length1
...
nameN pathN lengthN
```

The `input.txt` file mentioned looks like as follows:
```
simGorilla https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/primates/loci1/simGorilla.chr6 0.075
```

### Adding-to-a-node Case

This example takes the `simGorilla` genome from the evolverPrimates [example](https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverPrimates.txt) and adds it as a **child** of the `mr` node in the evolverMammals alignment [example](https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverMammals.txt). It updates the given alignment  at a cost of **one** Cactus alignment only. The update procedure is a [two-fold](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/updating-alignments.md#adding-to-a-node) step. First, it update the `mr` node's alignment blocks by re-aligning its current children (`simMouse_chr6` and `simRat_chr6`) with `simGorilla`. Then, it uses the `halReplaceGenome` (taken from the [HAL API](https://github.com/ComparativeGenomicsToolkit/hal)) to update the `mr` node's alignment blocks. Bearing in mind that this update procedure does *NOT* change the `mr` genome. The genome represented by the `mr` node remains the same.

For this example, run:
```
cactus-update-prepare add node ./evolverMammals.hal ./input.txt --genome mr --outDir ./steps --jobStore ./jobstore
```
where 
- `add` indicates the action of adding a new genome to the given alignment
- `node` indicates the approach to update the given alignment
- `evolverMammals.hal` is the HAL-format file containing the current `evolverMammals` alignment
- `input.txt` is the text tab-spaced file providing a list of genomes to be added to an existing node
- `mr` is the target genome node receiving new genome(s)
- `--outDir` and `--jobStore` are options inherited from [`cactus-prepare`](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/src/cactus/progressive/cactus_prepare.py)
    - `steps` is the name of the directory where assemblies and all cactus-call dependencies will be placed
    - `jobstore` is the name of the directory used by [Toil](https://toil.readthedocs.io/en/latest/) (Cactus internal scheduler) to store temporary data of job executions



After running the above, the output of `cactus-update-prepare` is expected to be as follows:
```
## Preprocessor
cactus-preprocess ./jobstore/2 ./steps/seq_file.in ./steps/seq_file.out --inputNames simGorilla --realTimeLogging --logInfo --retryCount 0

## Alignment

### Round 0
cactus-blast ./jobstore/3 ./steps/seq_file.out ./steps/mr.cigar --root mr --realTimeLogging --logInfo --retryCount 0 --includeRoot 
cactus-align ./jobstore/4 ./steps/seq_file.out ./steps/mr.cigar ./steps/mr.hal --root mr --realTimeLogging --logInfo --retryCount 0 --maxCores 2 --includeRoot 

## Alignment update
halReplaceGenome --bottomAlignmentFile ./steps/mr.hal --topAlignmentFile ./evolverMammals.hal ./evolverMammals.hal mr --hdf5InMemory 

## Alignment validation
halValidate --genome mr ./evolverMammals.hal --hdf5InMemory
```


The updated `evolverMammals` alignment tree, now containing the `simGorilla` node as a child of the `mr` node, looks like as follows: 
```
((simHuman_chr6:0.144018,(simMouse_chr6:0.084509,simRat_chr6:0.091589,simGorilla:0.075)mr:0.271974)Anc1:0.020593,(simCow_chr6:0.18908,simDog_chr6:0.16303)Anc2:0.032898)Anc0;
```

### Adding-to-a-branch Case

This example takes the `simGorilla` genome from the evolverPrimates [example](https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverPrimates.txt) and adds it to the **branch** connecting  the `Anc1` and `mr`  nodes in the evolverMammals alignment [example](https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/master/examples/evolverMammals.txt). This update procedure is a [four-fold](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/updating-alignments.md#adding-to-a-branch) step and updates the given alignment  at a cost of **two** Cactus runs . First, a new ancestor (labelled in this example as `AncGorilla`) between `Anc1` and `mr` is designed. By doing so, `Anc1` becomes the parent of `AncGorilla` and `mr` becomes the child of `AncGorilla`. Secondly, Cactus is called to infer a genome for the  `AncGorilla` node. This step creates a HAL-format file `AncGorilla.hal` containing the alignment information between `simGorilla` and `mr`. Next, Cactus is called once again to align the recently generated `AncGorilla` with `simHuman_chr6` (a children of `Anc1`). This step produces `Anc1.hal` containing updated alignment information for the `Anc1` node. Lastly, the `halAddToBranch` tool (taken from the [HAL API](https://github.com/ComparativeGenomicsToolkit/hal)) comes into play to update the  alignment blocks of `Anc1` in the initial alignment HAL file. Bearing in mind that this update procedure does *NOT* change the `Anc1` and `mr` genomes. The genomes represented by these nodes remains the same.

For this example, run:
```
cactus-update-prepare add branch --parentGenome Anc1 --childGenome mr ./evolverMammals.hal ./input.txt --cactus-prepare-options '--alignCores 4' --outDir ./steps --jobStore ./jobstore --ancestorName AncGorilla --topBranchLength 0.10
```
where:
- `add` indicates the action of adding a new genome to the given alignment
- `branch` indicates the approach to update the given alignment
- `evolverMammals.hal` is the HAL-format file containing the current `evolverMammals` alignment
- `input.txt` is the text tab-spaced file providing a list of genomes to be added to a branch
- `--parentGenome` is the top-level genome name of the target branch
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
halValidate --genome mr ./evolverMammals.hal --hdf5InMemory
halValidate --genome simGorilla ./evolverMammals.hal --hdf5InMemory

```


The updated `evolverMammals` alignment tree, now containing the `simGorilla` within the branch connecting  `Anc1` and `mr` nodes, looks like as follows:
```
(((simMouse_chr6:0.084509,simRat_chr6:0.091589)mr:0.271974,(simHuman_chr6:0.044018,simGorilla:0.075)AncGorilla:0.1)Anc1:0.020593,(simCow_chr6:0.18908,simDog_chr6:0.16303)Anc2:0.032898)Anc0;
```

## Replacing a genome

The replacement of a genome (due to an assembly version update) is carried out by removing a genome and then re-adding the new version by using the `add-to-a-node` procedure as a child of its parent.

The example below shows how to update the `simHuman_chr6` genome (which parent is `Anc1`) in the `evolverMammals` example:
```
cactus-update-prepare replace ./evolverMammals.hal ./input.txt --genome simHuman_chr6 --outDir ./steps --jobStore ./jobstore
```
where 
- `replace` indicates the action of replacing a genome to the given alignment
- `evolverMammals.hal` is the HAL-format file containing the current `evolverMammals` alignment
- `input.txt` is the text tab-spaced file providing the genome to be replaced
- `simHuman_chr6` is the target genome being replaced
- `--outDir` and `--jobStore` are options inherited from [`cactus-prepare`](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/src/cactus/progressive/cactus_prepare.py)
    - `steps` is the name of the directory where assemblies and all cactus-call dependencies will be placed
    - `jobstore` is the name of the directory used by [Toil](https://toil.readthedocs.io/en/latest/) (Cactus internal scheduler) to store temporary data of job executions

The output of `cactus-prepare-update` above will looks like as follows:

```
## Preprocessor
cactus-preprocess ./jobstore/1 ./steps/seq_file.in ./steps/seq_file.out --inputNames simHuman_chr6 --realTimeLogging --logInfo --retryCount 0

## Alignment

### Round 0
cactus-blast ./jobstore/2 ./steps/seq_file.out ./steps/Anc1.cigar --root Anc1 --realTimeLogging --logInfo --retryCount 0 --includeRoot 
cactus-align ./jobstore/3 ./steps/seq_file.out ./steps/Anc1.cigar ./steps/Anc1.hal --root Anc1 --realTimeLogging --logInfo --retryCount 0 --maxCores 2 --includeRoot 

## Alignment update
halReplaceGenome --bottomAlignmentFile ./steps/Anc1.hal --topAlignmentFile ./evolverMammals.hal ./evolverMammals.hal Anc1 --hdf5InMemory 

## Alignment validation
halValidate --genome Anc1 ./evolverMammals.hal --hdf5InMemory
```