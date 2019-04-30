# MerlinsManuscript
This repository houses the complete materials for reproducing microbiota-related results published in: <br/>

**Host and gut microbial metabolism of *Bifidobacterium longum*-fermented rice bran and rice bran in healthy mice** <br/>
Nealon et al. 2019 <br/>
Journal: *Beneficial Microbes* <br/>
doi: (INSERT FINAL INFO) <br/>
citation: (INSERT FINAL INFO) <br/>

**directory/file descriptions are as follows:** <br/>
`class/` = contains reference sequences used to train taxonomy classifiers (see: `Code_S1.R`) <br/>
`dada2/` = contains outputs for dada2 feature table construction (see: `Code_S1.R`) <br/>
`main/` = pdfs for main text figures <br/>
`raw_fastq/` = demultiplexed (individual) fastq sequence files for each sample <br/>
`sra/` = items relevant for submission to the NCBI Sequence Read Archive ([BioProject Accession PRJNA516457](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA516457)) <br/>
`submission/` = items produced for original submission (pre-revisions) to *Beneficial Microbes* <br/>
`supplement/` = includes figures/tables/files for supplementary materials <br/>
`taxa/` = contains Greengenes/SILVA taxonomic classifications for each FeatureID (see: Code_S1.R) <br/>
`vault/` = central storage for items produced by `Code_S1.R` and `Code_S2.R` <br/>
`Code_S1.R` = manifest creation and QIIME 2 processing <br/>
`Code_S2.R` = post QIIME 2 analysis conducted in R <br/>
`MerlinsManuscript.Rproj` = base of operations for running `Code_S1.R` and `Code_S2.R` in R Studio <br/>
`MetadataFile_S1.txt` = sample metadata file <br/>
`manifest_R1.csv` = manifest file for importing demultplexed fastq files into QIIME 2 (see: `Code_S1.R`) <br/>
`seqs_R1_.qzv` = sequence read visualization artifact <br/>
`seqs_R1_trim.qzv` = trimmed sequence read visualization artifact <br/>

A few files were unable to be uploaded int their respective locations due to file size limitations. To overcome this restriction, a master directory with complete materials was created, compressed, and uploaded to Google Drive. This directory can be downloaded via clicking "Download anyway" at the following: [MerlinsManuscript.tar.gz](https://drive.google.com/uc?export=download=w4JS&id=1w3rJhchSeyjtWiUCoUCXveMOkqpX6bZ5) <br/>
