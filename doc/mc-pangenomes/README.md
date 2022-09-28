# Minigraph-Cactus Pangenomes

These are some pangenomes created using [Minigraph-Cactus](../pangenome.md), along with all steps used to produce them.

The commits used can be found in the logs in the Output links. Data is provided in HAL, vg, GFA and VCF formats. Please read about [VCF output](../pangenome.md#vcf-output) before using it.  

Please cite the [Minigraph-Cactus paper](in prep) when using these pangenomes. For the HPRC pangenomes, please also cite the [HPRC Paper](https://doi.org/10.1101/2022.07.09.499321 ).

|<sub>**Name**</sub>| <sub>**Species**</sub> | <sub> **Date** </sub> |<sub>**Haplotypes**</sub> | <sub>**Reference**</sub> | <sub>**SeqFile**</sub> | <sub>**Commands**</sub> | <sub>**Output**</sub>|
| :-------- | :-------- | :-------- | :------ | :------ |:------ | :------ | :------ |
| <sub> Cow </sub> | <sub> *Bos taurus* </sub> | <sub> 2022-09-22 </sub> | <sub> 5 </sub> | <sub> bosTau9 </sub> | <sub> [5-cow-pg-2022-09-22-seqfile.txt](5-cow-pg-2022-09-22-seqfile.txt) </sub> | <sub> [5-cow-pg-2022-09-22-commands.md](5-cow-pg-2022-09-22-commands.md) </sub> | <sub> s3://vg-k8s/users/hickey/cow-pangenome// </sub> | 
| <sub> Chicken </sub> | <sub> *Gallus gallus* </sub> | <sub> 2022-09-23 </sub> | <sub> 10 </sub> | <sub> galGal6 </sub> | <sub> [10-chicken-pg-2022-09-23-seqfile.txt](10-chicken-pg-2022-09-23-seqfile.txt) </sub> | <sub> [10-chicken-pg-2022-09-23-commands.md](10-chicken-pg-2022-09-23-commands.md) </sub> | <sub> s3://vg-k8s/users/hickey/chicken-pangenome// </sub> | 
| <sub> Dog </sub> | <sub> *Canis lupus familiaris* </sub> | <sub> 2022-09-23 </sub> | <sub> 9 </sub> | <sub> canFam4 </sub> | <sub> [9-dog-pg-2022-09-23-seqfile.txt](9-dog-pg-2022-09-23-seqfile.txt) </sub> | <sub> [9-dog-pg-2022-09-23-commands.md](9-dog-pg-2022-09-23-commands.md) </sub> | <sub> s3://vg-k8s/users/hickey/dog-pangenome// </sub> | 
| <sub> Fruit fly </sub> | <sub> *Drosophila melanogaster* </sub> | <sub> 2022-05-26 </sub> | <sub> 16 </sub> | <sub> dm6 </sub> | <sub> seqfile </sub> | <sub> steps </sub> | | <sub> data </sub> |
| <sub> HPRC v1.0 </sub> | <sub> *Homo sapiens* </sub> | <sub> 2021-08-11 </sub> | <sub> 90 </sub> | <sub> GRCh38 </sub> | <sub> [see here](https://github.com/human-pangenomics/hpp_pangenome_resources/) </sub> | <sub> [see here](../pangenome.md#hprc-version-1.0-graphs) </sub> | <sub> [see here](https://github.com/human-pangenomics/hpp_pangenome_resources/) </sub> |
| <sub> HPRC v1.0 </sub> | <sub> *Homo sapiens* </sub> | <sub> 2021-08-11 </sub> | <sub> 90 </sub> | <sub> CHM13 </sub> | <sub> [see here](https://github.com/human-pangenomics/hpp_pangenome_resources/) </sub> | <sub> [see here](../pangenome.md#hprc-version-1.0-graphs) </sub> | <sub> [see here](https://github.com/human-pangenomics/hpp_pangenome_resources/) </sub> |
| <sub> Mouse </sub> | <sub> *Mus musculus* </sub> | <sub> 2022-09-23 </sub> | <sub> 30 </sub> | <sub> mm39 </sub> | <sub> [30-mouse-pg-2022-09-23-seqfile.txt](30-mouse-pg-2022-09-23-seqfile.txt) </sub> | <sub> [30-mouse-pg-2022-09-23-commands.md](30-mouse-pg-2022-09-23-commands.md) </sub> | <sub> s3://vg-k8s/users/hickey/mouse-pangenome// </sub> | 
| <sub> Soybean </sub> | <sub> *Glycine max* </sub> | <sub> 2022-09-26 </sub> | <sub> 17 </sub> | <sub> Glycine_max_v4.0 </sub> | <sub> [17-soybean-pg-2022-09-26-seqfile.txt](17-soybean-pg-2022-09-26-seqfile.txt) </sub> | <sub> [17-soybean-pg-2022-09-26-commands.md](17-soybean-pg-2022-09-26-commands.md) </sub> | <sub> s3://vg-k8s/users/hickey/soybean-pangenome// </sub> |

Creating an HPRC seqfile and running the current pipeline on it is described [here](../pangenome.md#hprc-graph).
