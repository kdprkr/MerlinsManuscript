## MerlinsManuscript
This repository houses the complete materials for reproducing microbiota-related results published in: <br/>

## Health implications for dietary supplementation with *Bifidobacterium longum*-fermented rice bran and rice bran as examined through  gut microbiomes and metabolomics 
Nealon et al. 2019 <br/>
Journal: *Beneficial Microbes* <br/>
[pubmed link ~ pending](https://github.com/kdprkr/MerlinsManuscript) <br/>
[journal link](https://www.wageningenacademic.com/doi/abs/10.3920/BM2019.0017) <br/>
citation:
Nealon NJ, KD Parker, P Lahaie, H Ibrahim, AK Maurya, K Raina, EP Ryan. 2019. Bifidobacterium longum-fermented rice bran and rice bran supplementation affects the gut microbiome and metabolome. Beneficial Microbes 10(8) 823-839. doi: 10.3920/BM2019.0017. <br/>

**Abstract:** This study investigated gut microbiota composition along with food, host, and microbial derived metabolites in the colon and systemic circulation of healthy mice following dietary rice bran and fermented rice bran intake. Adult male BALB/c mice were fed a control diet or one of two experimental diets containing 10% w/w rice bran fermented by *Bifidobacterium longum* or 10% w/w non-fermented rice bran for 15 weeks. Metabolomics was performed on the study diets (food), the murine colon and whole blood. These were analysed in concert with 16S rRNA amplicon sequencing of faeces, caecum, and colon microbiomes. Principal components analysis of murine microbiota composition displayed marked separation between control and experimental diets, and between faecal and tissue (caecum and colon) microbiomes. Colon and caecal microbiomes in both experimental diet groups showed enrichment of *Roseburia*, *Lachnospiraceae*, and Clostridiales related amplicon sequence variants compared to control. Bacterial composition was largely similar between experimental diets. Metabolite profiling revealed 530 small molecules comprising of 39% amino acids and 21% lipids that had differential abundances across food, colon, and blood matrices, and statistically significant between the control, rice bran, and fermented rice bran groups. The amino acid metabolite, N-delta-acetylornithine, was notably increased by *B. longum* rice bran fermentation when compared to non-fermented rice bran in food, colon, and blood. These findings support that dietary intake of rice bran fermented with *B. longum* modulates multiple metabolic pathways important to the gut and overall health.  <br/>

## Overview
`class/` = contains reference sequences used to train taxonomy classifiers (see: `Code_S1.R`) <br/>
`dada2/` = contains outputs for dada2 feature table construction (see: `Code_S1.R`) <br/>
`main/` = pdfs for main text figures <br/>
`raw_fastq/` = demultiplexed (individual) fastq sequence files for each sample <br/>
`revision_01/` = items produced for revised submission (1st round) to *Beneficial Microbes* <br/>
`revision_02/` = items produced for revised submission (2nd round) to *Beneficial Microbes* <br/>
`sra/` = items relevant for submission to the NCBI Sequence Read Archive ([BioProject Accession PRJNA516457](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA516457)) <br/>
`submission/` = items produced for original submission (pre-revisions) to *Beneficial Microbes* <br/>
`supplement/` = includes figures/tables/files for supplementary materials <br/>
`taxa/` = contains Greengenes/SILVA taxonomic classifications for each FeatureID (see: Code_S1.R) <br/>
`vault/` = central storage for items produced by `Code_S1.R` and `Code_S2.R` <br/>
`Code_S1.R` = manifest creation and QIIME 2 processing <br/>
`Code_S2.R` = post QIIME 2 analysis conducted in R <br/>
`MerlinsManuscript.Rproj` = base of operations for running `Code_S1.R` and `Code_S2.R` in R Studio <br/>
`MetadataFile_S1.txt` = sample data file <br/>
`manifest_R1.csv` = manifest file for importing demultplexed fastq files into QIIME 2 (see: `Code_S1.R`) <br/>
`seqs_R1.qzv` = sequence read visualization artifact <br/>
`seqs_R1_trim.qzv` = trimmed sequence read visualization artifact <br/>

Due to file size restrictions, a few items were unable to be uploaded into their respective locations above. To overcome this limitation, the master directory with complete materials was compressed and uploaded to Google Drive. This directory can be downloaded via clicking "Download anyway" after visiting the following [link](https://drive.google.com/uc?id=1vfk0-mTub7uidIflR8YNGLPOgoFeqGGf&export=download). To setup proper file paths, download, uncompress, and move `MerlinsManuscript/` onto the `Desktop/`. <br/>

If any files cannot be located or if any links are broken, **please open an issue**. <br/>

If you have any questions regarding any items here, **please contact me via the email listed on my profile**. <br/>
