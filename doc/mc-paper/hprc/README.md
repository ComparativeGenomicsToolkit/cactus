# HPRC Pangenome Methods

## Construction

The graphs were constructed as described [here](../pangenome.md#hprc-version-10-graphs).  They are available to download [here](https://github.com/human-pangenomics/hpp_pangenome_resources/).

The versions of the graph used to create Supplementary Figure 3 in order to show sequence removed when not using dna-brnn were created with the commands listed [here](../pangenome.md#hprc-graph).

The non-reference sequence was computed using [count-vg-hap-cov](https://github.com/ComparativeGenomicsToolkit/hal2vg/blob/f3d9a1838d1fb5582b6e1cd509792daee51fd2a9/count-vg-hap-cov.cpp) on the [GRCh38](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=pangenomes/scratch/2021_08_11_minigraph_cactus/GRCh38-chrom-graphs/) and [CHM13](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=pangenomes/scratch/2021_08_11_minigraph_cactus/CHM13-chrom-graphs/) chromosome vg graphs.

## Snakemake pipeline and scripts

[Snakemake](https://snakemake.github.io/) is used for multiple analysis, to help run the same analysis across many samples or approaches.
Our pipeline, implemented in the Snakefiles (`Snakefile_*`), is set up to work with our data located in our private S3 buckets.
Hence, to reproduce the analysis, one would have to replace the parts specifying those S3 paths (`S3.remote(...)`), with local paths where the raw data was downloaded.

The pipeline also uses custom scripts or additional files (`resources/..'` in the Snakefile). 
Those files have been deposited in the [`resources`](resources) folder.

## Non-reference sequence

[count-vg-hap-cov](https://github.com/ComparativeGenomicsToolkit/hal2vg/releases/download/v1.1.2/count-vg-hap-cov) was used to compute the amount of sequence (reference and non-reference) in the graph according to how many haplotypes cover it. It was run on the using `-r GRCh38` (`CHM13`) to specify the reference of the GRCh38 (CHM13) based graphs. The chromosome vg output of `cactus-graphmap-join` was used as input.

## Snarl summary statistics

The distance index contains the snarl structure and records the length of the maximum/minimum path through each.
We extract this information using `vg view -B`, counting the number of snarls for each profile (size x depth x minimum length x maximum length).

### Minigraph pangenome

The GFA is converted to a VG graph which is then indexed.

```sh
## Minigraph pangenome
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph/hprc-v1.0-minigraph-grch38.gfa.gz
gunzip hprc-v1.0-minigraph-grch38.gfa.gz
vg view -vF hprc-v1.0-minigraph-grch38.gfa > hprc-v1.0-minigraph-grch38.vg
vg snarls -T hprc-v1.0-minigraph-grch38.vg > hprc-v1.0-minigraph-grch38.snarls
vg index -s hprc-v1.0-minigraph-grch38.snarls -j hprc-v1.0-minigraph-grch38.dist hprc-v1.0-minigraph-grch38.vg
```

The snarl statistics are then compiled:

```sh
echo -e "n\tnode_count\tdepth\tmin_length\tmax_length" > hprc-v1.0-minigraph-grch38.dist-stats.tsv
vg view -B hprc-v1.0-minigraph-grch38.dist | jq -r 'select(.type=="snarl") | select(.node_count>2) | [.node_count, .depth, .minimum_length, .maximum_length] | @tsv' | sort | uniq -c >> hprc-v1.0-minigraph-grch38.dist-stats.tsv
gzip hprc-v1.0-minigraph-grch38.dist-stats.tsv
```

### Minigraph-Cactus graph

The distance is already available for this pangenome:

```sh
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/filtered/hprc-v1.0-mc-grch38-maxdel.10mb.dist
echo -e "n\tnode_count\tdepth\tmin_length\tmax_length" > hprc-v1.0-mc-grch38-maxdel.10mb.dist-stats.tsv
vg view -B hprc-v1.0-mc-grch38-maxdel.10mb.dist | jq -r 'select(.type=="snarl") | select(.node_count>2) | [.node_count, .depth, .minimum_length, .maximum_length] | @tsv' | sort | uniq -c >> hprc-v1.0-mc-grch38-maxdel.10mb.dist-stats.tsv
gzip hprc-v1.0-mc-grch38-maxdel.10mb.dist-stats.tsv
```

### Figures

The two `*.dist-stats.tsv` files created above were read by the [`snarls-stats.R`](snarls-stats.R) to make the figure shown in the manuscript.

## Mapping statistics

We mapped short reads (30x PCR-free novaseq) for HG001/2/5 to the reference genomes and the different pangenomes.
Those alignments were then parsed to compute the mapping statistics.

- The input pangenomes can be found at: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=pangenomes/freeze/freeze1/minigraph-cactus/filtered/
- The input FASTQs were downloaded from the [DeepVariant benchmarking dataset bucket](https://console.cloud.google.com/storage/browser/deepvariant/benchmarking/fastq/wgs_pcr_free/30x).

For both mapping and compiling the statistics we used the [`Snakefile_mapstats`](Snakefile_mapstats):

```sh
snakemake -s Snakefile_mapstats --cores 16 -p mapstats_bwa
snakemake -s Snakefile_mapstats --cores 16 -p mapstats_giraffe
snakemake -s Snakefile_mapstats --cores 16 --config graph=hprc_v1_0_mc_chm13_minaf_0_1 -p mapstats_giraffe
snakemake -s Snakefile_mapstats --cores 16 --config graph=hprc_v1_0_mc_grch38_minaf_0_1_noclip -p mapstats_giraffe
```

Internally, the workflow will use the [`resources/compute_mapping_stats.py`](resources/compute_mapping_stats.py) script to compute and tally reads in each profile in the BAM file.
For example, it will parse both the CIGAR field and the MD tag to figure out if the read is aligning perfectly to the reference genome.

The `<SAMPLE>.<MAPPER>.<REF>.mapstats.txt` files created by these runs contains the number of reads for each profile (*mapping quality* x *perfect alignment*).
They are read by the [`mapstats-analysis.R`](mapstats-analysis.R) script to make the graphs.
A TAR file containing all of those `*mapstats.txt` files for HG001,2,5 and the different reference (pan)genomes is shared at [https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/hprc-human/mapstats.hg001_2_5.giraffe_hprc.bwa_grch38.tar.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/hprc-human/mapstats.hg001_2_5.giraffe_hprc.bwa_grch38.tar.gz).

## Small variant calling with DeepVariant

### Calling the small variants on Terra with the WDL workflow

The small variants were called on [Terra](https://app.terra.bio/) using the public [GiraffeDeepVariant (tag:giraffedv-v1.3-mc) WDL workflow](https://dockstore.org/workflows/github.com/vgteam/vg_wdl/GiraffeDeepVariant:giraffedv-v1.3-mc?tab=info) deposited on Dockstore.

Examples of the inputs are provided as JSON files in the [`terra-files`](terra-files) folder:
- Aligning to the GRCh38-based pangenome and projecting the reads to GRCh38: [`terra-files/giraffe-deepvariant-grch38_pangenome-grch38_projection.json`](terra-files/giraffe-deepvariant-grch38_pangenome-grch38_projection.json)
- Aligning to the CHM13-based pangenome and projecting the reads to GRCh38: [`terra-files/giraffe-deepvariant-chm13_pangenome-grch38_projection.json`](terra-files/giraffe-deepvariant-chm13_pangenome-grch38_projection.json)
- Aligning to the CHM13-based pangenome and projecting the reads to CHM13: [`terra-files/`](terra-files/)
- Aligning to GRCh38 with BWA-MEM: [`terra-files/bwa_deepvariant.wdl`](terra-files/bwa_deepvariant.wdl) and [`terra-files/bwa-deepvariant-chm13.json`](terra-files/bwa-deepvariant-chm13.json)

Other files used for the runs were also place in [`terra-files`](terra-files):
- [`terra-files/CHM13.path_list.txt`](terra-files/CHM13.path_list.txt)
- [`terra-files/GRCh38.path_list.txt`](terra-files/GRCh38.path_list.txt)

`chm13v2.0.plus_hs38d1_analysis_set.compact_decoys.fa` was created by adding the [hs38d1 decoy sequences](https://www.ncbi.nlm.nih.gov/assembly/GCA_000786075.2/) to [CHM12 v2.0](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz).

The DeepVariant models, trained on the HPRC mappings, are available at [https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=publications/PANGENOME_2022/DeepVariant/models/DEEPVARIANT_MC_Y1/](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=publications/PANGENOME_2022/DeepVariant/models/DEEPVARIANT_MC_Y1/).
As for the mapping statistics, we used FASTQs for the HG001/2/5 samples (30x PCR-free novaseq) from the [DeepVariant benchmarking dataset bucket](https://console.cloud.google.com/storage/browser/deepvariant/benchmarking/fastq/wgs_pcr_free/30x).

The VCF files produced are available at:

1. [`https://public.gi.ucsc.edu/~jmonlong/hprc/minigraph-cactus-manuscript/hprc-human/hprc-calls-on-grch38/`](https://public.gi.ucsc.edu/~jmonlong/hprc/minigraph-cactus-manuscript/hprc-human/hprc-calls-on-grch38/) for the calls on the GRCh38 reference.
1. [`https://public.gi.ucsc.edu/~jmonlong/hprc/minigraph-cactus-manuscript/hprc-human/hprc-calls-on-chm13/`](https://public.gi.ucsc.edu/~jmonlong/hprc/minigraph-cactus-manuscript/hprc-human/hprc-calls-on-chm13/) for the calls on the CHM13 reference

### Evaluating the calls against the GIAB/CMRG truthsets

The VCFs with the calls for each sample and methods was downloaded to our S3 bucket at `s3://vg-k8s/vgamb/wg/giraffedv/giab-calls/{sample}.30x_novaseq_pcrfree.{method}.vcf.gz`.
This is our private S3 bucket. 
We mention this because that's where the Snakemake pipeline will look for files.
Hence, that is the path that should be changed in order to reproduce this analysis.

From there, we ran the evaluation (using [hap.py](https://github.com/Illumina/hap.py)) with Snakemake:

```sh
snakemake -s Snakefile_eval_grch38 --configfile snakemake_config_hg002.yaml --cores 16 -p main
snakemake -s Snakefile_eval_grch38 --configfile snakemake_config_hg002.yaml --cores 1 -p main_roc
snakemake -s Snakefile_eval_grch38 --configfile snakemake_config_giab1-2-5.yaml --cores 16 -p main
snakemake -s Snakefile_eval_grch38 --configfile snakemake_config_giab1-2-5.yaml --cores 1 -p main_roc
```

The files from the GIAB v4.2.1 and CMRG v1.0 truthsets are downloaded by the workflow, see the `dwl_*` rules in the `Snakefile_eval*` files.

#### Lifting calls from CHM13 to GRCh38

The calls from both BWA-DeepVariant and Giraffe-DeepVariant were first lifted from CHM13 v2.0 to GRCh38.

We download the chain file and prepare the genome for Picard:

```sh
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-grch38.chain
java -jar ~/build/picard.jar CreateSequenceDictionary R=hg38.fa
```

Then, for each method and sample, we download split the variants into biallelic records in the VCF, lift it over, and merge back the records into multiallelic ones.
The output VCF is uploaded to our private S3 bucket to be used later by the evaluation workflow.

```sh
for SAMP in HG001 HG002 HG005
do
    ## Liftover BWA calls GRch38 -> CHM13 v2
    bcftools norm -Oz -o ${SAMP}.bwadv.chm13.bial.vcf.gz -m -any ${SAMP}.bwadv.chm13.vcf.gz
    java -jar ~/build/picard.jar LiftoverVcf I=${SAMP}.bwadv.chm13.bial.vcf.gz \
         O=${SAMP}.bwadv.chm13.bial.liftedGRCh38.vcf.gz CHAIN=chm13v2-grch38.chain \
         REJECT=${SAMP}.bwadv.chm13.notliftedGRCh38.vcf.gz R=hg38.fa RECOVER_SWAPPED_REF_ALT=true
    bcftools norm -Oz -o ${SAMP}.bwadv.chm13.liftedGRCh38.vcf.gz -m +any ${SAMP}.bwadv.chm13.bial.liftedGRCh38.vcf.gz
    aws s3 cp ${SAMP}.bwadv.chm13.liftedGRCh38.vcf.gz s3://vg-k8s/vgamb/wg/giraffedv/giab-calls/${SAMP}.30x_novaseq_pcrfree.bwadv_chm13_lifted.vcf.gz
    aws s3 cp ${SAMP}.bwadv.chm13.vcf.gz s3://vg-k8s/vgamb/wg/giraffedv/giab-calls-chm13/${SAMP}.30x_novaseq_pcrfree.bwadv.vcf.gz
    ## liftover Giraffe-DV calls GRch38 -> CHM13 v2
    bcftools norm -Oz -o ${SAMP}.giraffedv.chm13.bial.vcf.gz -m -any ${SAMP}.giraffedv.chm13.vcf.gz
    java -jar ~/build/picard.jar LiftoverVcf I=${SAMP}.giraffedv.chm13.bial.vcf.gz \
         O=${SAMP}.giraffedv.chm13.bial.liftedGRCh38.vcf.gz CHAIN=chm13v2-grch38.chain \
         REJECT=${SAMP}.giraffedv.chm13.notliftedGRCh38.vcf.gz R=hg38.fa RECOVER_SWAPPED_REF_ALT=true
    bcftools norm -Oz -o ${SAMP}.giraffedv.chm13.liftedGRCh38.vcf.gz -m +any ${SAMP}.giraffedv.chm13.bial.liftedGRCh38.vcf.gz
    aws s3 cp ${SAMP}.giraffedv.chm13.liftedGRCh38.vcf.gz s3://vg-k8s/vgamb/wg/giraffedv/giab-calls/${SAMP}.30x_novaseq_pcrfree.giraffedv_chm13_lifted.vcf.gz
    aws s3 cp ${SAMP}.giraffedv.chm13.vcf.gz s3://vg-k8s/vgamb/wg/giraffedv/giab-calls-chm13/${SAMP}.30x_novaseq_pcrfree.giraffedv.vcf.gz
done
```

Some variants, that would appear when comparing to GRCh38, won't be present in the CHM13-based callset because they are homozygous for the CHM13 allele.
To remove those "invisible" variants from the GRCh38 truthset, we lift them over, from GRCh38 to CHM13, and identify which ones become homozygous ref.

```sh
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/grch38-chm13v2.chain
java -jar ~/build/picard.jar CreateSequenceDictionary R=chm13v2.0.plus_hs38d1_analysis_set.compact_decoys.fa

## liftover GIAB truthset (to filter out variants that can't be lifted or result in hom ref on CHM13)
for SAMP in HG001 HG002 HG005
do
    aws s3 cp s3://vg-k8s/vgamb/wg/giraffedv/giab-evaluation/truth/GIAB_4_2_1.${SAMP}.vcf.gz .
    bcftools norm -m -any GIAB_4_2_1.${SAMP}.vcf.gz | \
        bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o GIAB_4_2_1.${SAMP}.bial.ids.vcf.gz 
    java -jar ~/build/picard.jar LiftoverVcf I=GIAB_4_2_1.${SAMP}.bial.ids.vcf.gz  \
         O=GIAB_4_2_1.${SAMP}.bial.ids.lifted.vcf.gz CHAIN=grch38-chm13v2.chain \
         REJECT=GIAB_4_2_1.${SAMP}.bial.ids.notlifted.vcf.gz \
         R=chm13v2.0.plus_hs38d1_analysis_set.compact_decoys.fa RECOVER_SWAPPED_REF_ALT=true 2> LiftOverVcf.GIAB_4_2_1.${SAMP}.log
    zcat GIAB_4_2_1.${SAMP}.bial.ids.lifted.vcf.gz | \
        grep SwappedAlleles | grep -e "0/0" -e "0|0" | cut -f3 > discard_ids.txt
    zcat GIAB_4_2_1.${SAMP}.bial.ids.notlifted.vcf.gz | cut -f3 >> discard_ids.txt
    bcftools filter -e 'ID=@discard_ids.txt' GIAB_4_2_1.${SAMP}.bial.ids.vcf.gz | \
        bcftools norm -Oz -o GIAB_4_2_1_chm13visible.${SAMP}.vcf.gz -m +any 
    aws s3 cp GIAB_4_2_1_chm13visible.${SAMP}.vcf.gz s3://vg-k8s/vgamb/wg/giraffedv/giab-evaluation/truth/GIAB_4_2_1_chm13visible.${SAMP}.vcf.gz
    aws s3 cp s3://vg-k8s/vgamb/wg/giraffedv/giab-evaluation/truth/GIAB_4_2_1.${SAMP}.wg_noinconsistent.bed.gz s3://vg-k8s/vgamb/wg/giraffedv/giab-evaluation/truth/GIAB_4_2_1_chm13visible.${SAMP}.wg_noinconsistent.bed.gz
done

## liftover CMRG truthset (to filter out variants that can't be lifted or result in hom ref on CHM13)
aws s3 cp s3://vg-k8s/vgamb/wg/giraffedv/giab-evaluation/truth/CMRG_1_0.HG002.vcf.gz .
bcftools norm -m -any CMRG_1_0.HG002.vcf.gz | \
    bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o CMRG_1_0.HG002.bial.ids.vcf.gz 
java -jar ~/build/picard.jar LiftoverVcf I=CMRG_1_0.HG002.bial.ids.vcf.gz  \
     O=CMRG_1_0.HG002.bial.ids.lifted.vcf.gz CHAIN=grch38-chm13v2.chain \
     REJECT=CMRG_1_0.HG002.bial.ids.notlifted.vcf.gz \
     R=chm13v2.0.plus_hs38d1_analysis_set.compact_decoys.fa RECOVER_SWAPPED_REF_ALT=true 2> LiftOverVcf.CMRG_1_0.HG002.log
bcftools view -C 0 CMRG_1_0.HG002.bial.ids.lifted.vcf.gz | \
    grep SwappedAlleles | cut -f3 > discard_ids.txt
zcat CMRG_1_0.HG002.bial.ids.notlifted.vcf.gz | grep -v "#" | cut -f3 >> discard_ids.txt
bcftools filter -e 'ID=@discard_ids.txt' CMRG_1_0.HG002.bial.ids.vcf.gz | \
    bcftools norm -Oz -o CMRG_1_0_chm13visible.HG002.vcf.gz -m +any 
aws s3 cp CMRG_1_0_chm13visible.HG002.vcf.gz s3://vg-k8s/vgamb/wg/giraffedv/giab-evaluation/truth/CMRG_1_0_chm13visible.HG002.vcf.gz
aws s3 cp s3://vg-k8s/vgamb/wg/giraffedv/giab-evaluation/truth/CMRG_1_0.HG002.wg_noinconsistent.bed.gz s3://vg-k8s/vgamb/wg/giraffedv/giab-evaluation/truth/CMRG_1_0_chm13visible.HG002.wg_noinconsistent.bed.gz
```

Once the VCFs and BED files lifted to GRCh38, the evaluation was run with:

```sh
snakemake -s Snakefile_eval_grch38 --configfile snakemake_config_lifted.yaml --cores 8 -p main
snakemake -s Snakefile_eval_grch38 --configfile snakemake_config_lifted_cmrg.yaml --cores 8 -p main
```

#### Evaluation on CHM13

To evaluate the calls on CHM13, we used the CMRG truthset from CHM13 v1.0 which we lifted to CHM13 v1.1. 

```sh
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/changes/v1.0_to_v1.1/v1.0_to_v1.1.chain

## Liftover CMRG VCF v1.0 to v1.1
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/CHM13v1.0/SmallVariant/HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.vcf.gz
bcftools norm -Oz -o HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.bial.vcf.gz -m -any HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.vcf.gz
java -jar ~/build/picard.jar LiftoverVcf I=HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.bial.vcf.gz \
     O=HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.bial.liftedGRCh38.vcf.gz CHAIN=v1.0_to_v1.1.chain \
     REJECT=HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.notliftedGRCh38.vcf.gz \
     R=chm13v2.0.plus_hs38d1_analysis_set.compact_decoys.fa RECOVER_SWAPPED_REF_ALT=true
zcat HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.bial.liftedGRCh38.vcf.gz | awk '{if($5!="*"){print $0}}' | \
    bcftools norm -Oz -o HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.liftedGRCh38.vcf.gz -m +any 
aws s3 cp HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.liftedGRCh38.vcf.gz s3://vg-k8s/vgamb/wg/giraffedv/giab-evaluation-chm13/truth/CMRG_1_0.HG002.vcf.gz

## Liftover CMRG confident regions v1.0 to v1.1
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/CHM13v1.0/SmallVariant/HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.bed
liftOver HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.bed v1.0_to_v1.1.chain \
         HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.lifted.bed \
         HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.notlifted.bed
gzip HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.lifted.bed
aws s3 cp HG002_CHM13v1.0_CMRG_smallvar_v1.00_draft.lifted.bed.gz s3://vg-k8s/vgamb/wg/giraffedv/giab-evaluation-chm13/truth/CMRG_1_0.HG002.wg_noinconsistent.bed.gz

## Liftover CMRG "dipcall" regions v1.0 to v1.1
wget --quiet https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/CHM13v1.0/SupplementaryFiles/HG002v11-align2-CHM13v1.0/HG002v11-align2-CHM13v1.0.dip.bed
liftOver HG002v11-align2-CHM13v1.0.dip.bed v1.0_to_v1.1.chain \
         HG002v11-align2-CHM13v1.0.dip.lifted.bed \
         HG002v11-align2-CHM13v1.0.dip.notlifted.bed
gzip HG002v11-align2-CHM13v1.0.dip.lifted.bed
aws s3 cp HG002v11-align2-CHM13v1.0.dip.lifted.bed.gz s3://vg-k8s/vgamb/wg/giraffedv/giab-evaluation-chm13/truth/GIAB_4_2_1_lifted.HG002.dipcall.bed.gz
```

We also lifted over the confident regions from GIAB v4.2.1 from GRCh38 to CHM13 v2.0.

```sh
aws s3 cp s3://vg-k8s/vgamb/wg/giraffedv/giab-evaluation/truth/GIAB_4_2_1.HG002.wg_noinconsistent.bed.gz .
gunzip GIAB_4_2_1.HG002.wg_noinconsistent.bed.gz
liftOver GIAB_4_2_1.HG002.wg_noinconsistent.bed grch38-chm13v2.chain \
         GIAB_4_2_1.HG002.wg_noinconsistent.lifted.bed \
         GIAB_4_2_1.HG002.wg_noinconsistent.notlifted.bed
gzip GIAB_4_2_1.HG002.wg_noinconsistent.lifted.bed
aws s3 cp GIAB_4_2_1.HG002.wg_noinconsistent.lifted.bed.gz s3://vg-k8s/vgamb/wg/giraffedv/giab-evaluation-chm13/truth/GIAB_4_2_1_lifted.HG002.lifted.bed.gz
```

Once the lifting over complete, we evaluate the calls with:

```sh
snakemake -s Snakefile_eval_chm13 --configfile snakemake_config_chm13.yaml --cores 8 -p main
```

#### Figures

The [`calls-evaluation-exploration.R`](calls-evaluation-exploration.R) script makes figures from the combined evaluation summary files (`eval-summary-*.tsv` and `eval-roc-summary-*.tsv`).

The combined summary files are available at [https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/hprc-human/eval-summary-hprc.tar.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/mc_2022/hprc-human/eval-summary-hprc.tar.gz).

#### Stratification of the GRCh38-based vs CHM13-based pangenomes performance

The stratifications from GIAB were first downloaded following the links in [https://github.com/genome-in-a-bottle/genome-stratifications/blob/master/GRCh38/v3.1-GRCh38-stratifications-all-except-genome-specific-stratifications.tsv](https://github.com/genome-in-a-bottle/genome-stratifications/blob/master/GRCh38/v3.1-GRCh38-stratifications-all-except-genome-specific-stratifications.tsv).
The [`resources/prepare_stratification_regions.R`](resources/prepare_stratification_regions.R) script downloads them all and compile them in an R object, then saved in `regions.all.giab.stratifications.RDS`.

The [`resources/get_stratification_on_sample.R`](resources/get_stratification_on_sample.R) then reads a tabulated summary of the VCF and overlaps it with the stratification regions to estimate the performance in each region set.

```sh
for SAMP in HG001 HG002 HG005
do
    ## GRCh38-based pangenome
    aws s3 cp s3://vg-k8s/vgamb/wg/giraffedv/giab-evaluation/happy_runs/GIAB_4_2_1.${SAMP}.30x_novaseq_pcrfree.giraffedv_chm13.nofalsedup_in_chm13/GIAB_4_2_1.${SAMP}.30x_novaseq_pcrfree.giraffedv_chm13.nofalsedup_in_chm13.vcf.gz happy_runs/
    bcftools query -f '%CHROM\t%POS\t%END[\t%BD\t%BVT]\n' \
             happy_runs/GIAB_4_2_1.${SAMP}.30x_novaseq_pcrfree.giraffedv_chm13.nofalsedup_in_chm13.vcf.gz | \
        awk '{if(($4!="." && $4 != "N") || ($6!="." && $6 != "N")){print $0}}' | \
        gzip > happy_runs/GIAB_4_2_1.${SAMP}.30x_novaseq_pcrfree.giraffedv_chm13.nofalsedup_in_chm13.locs.eval.tsv.gz
    Rscript get_stratification_on_sample.R \
            happy_runs/GIAB_4_2_1.${SAMP}.30x_novaseq_pcrfree.giraffedv.nofalsedup_in_chm13.locs.eval.tsv.gz \
            happy_runs/GIAB_4_2_1.${SAMP}.30x_novaseq_pcrfree.giraffedv.nofalsedup_in_chm13.stratification.rds
    ## CHM13-based pangenome
    aws s3 cp s3://vg-k8s/vgamb/wg/giraffedv/giab-evaluation/happy_runs/GIAB_4_2_1.${SAMP}.30x_novaseq_pcrfree.giraffedv.nofalsedup_in_chm13/GIAB_4_2_1.${SAMP}.30x_novaseq_pcrfree.giraffedv.nofalsedup_in_chm13.vcf.gz happy_runs/
    bcftools query -f '%CHROM\t%POS\t%END[\t%BD\t%BVT]\n' \
             happy_runs/GIAB_4_2_1.${SAMP}.30x_novaseq_pcrfree.giraffedv.nofalsedup_in_chm13.vcf.gz | \
        awk '{if(($4!="." && $4 != "N") || ($6!="." && $6 != "N")){print $0}}' | \
        gzip > happy_runs/GIAB_4_2_1.${SAMP}.30x_novaseq_pcrfree.giraffedv.nofalsedup_in_chm13.locs.eval.tsv.gz
    Rscript get_stratification_on_sample.R \
            happy_runs/GIAB_4_2_1.${SAMP}.30x_novaseq_pcrfree.giraffedv_chm13.nofalsedup_in_chm13.locs.eval.tsv.gz \
            happy_runs/GIAB_4_2_1.${SAMP}.30x_novaseq_pcrfree.giraffedv_chm13.nofalsedup_in_chm13.stratification.rds
done
```

Finally, those files are read by [`eval-stratification-grch38-chm13.R`](eval-stratification-grch38-chm13.R) to make the supplementary figure shown in the manuscript.
