### ************************************
### PREFACE ----
### ************************************

# this file houses the entirety of code for analysis performed in QIIME2*
# *after creating a manifest using R

# plan for paths accordingly:
# path from home directory on my machine = ~/Desktop/MerlinsManuscript/
# the working directory for QIIME2 will ALWAYS be: ~/Desktop/MerlinsManuscript/
# therefore, paths will ALWAYS be relative to MerlinsManuscript/
# raw fastq files are located in MerlinsManuscript/raw_fastq/
# database files for the tax classifier are located in the home directory (~/)
# several items are stored in the following directories:
# MerlinsManuscript/supplement/
# MerlinsManuscript/vault/

# *gasp*
setwd("~/Desktop/")

# relative paths from wd for inputs
ifp_met <- "MerlinsManuscript/MetadataFile_S1.txt"
M_ifp_fsq <- "MerlinsManuscript/raw_fastq/"

# relative paths from wd for outputs headed for storage in the "vault/"
M_ofv_mft <- "MerlinsManuscript/vault/manifest_R1.csv"
M_ofv_WS <- "MerlinsManuscript/vault/WS_ManifestR1.RData"

# relative paths from wd for outputs headed for storage in "MerlinsManuscript/"
ofp_mft <- "MerlinsManuscript/manifest_R1.csv"

### ************************************
### SECTION M - create manifest for QIIME2 ----
### ************************************

# purpose:
# this section creates a comma-separated "Fastq manifest" file (.csv)
# needed to import individual (demultiplexed) fastq files into QIIME2

# description of steps:
# STEP 1 = creates three columns (sample-id, absolute-filepath, direction)
# STEP 2 = creates a manifest and checks the manifest for errors

### ************************************
### M - STEP 1 - create columns  ----
### ************************************

# for column sample-id:
# obtain extension-less file names and truncate to retain SampleIDs only
M_fname <- gsub(pattern = ".fastq.gz", replacement = "",
                list.files(path = M_ifp_fsq, full.names = FALSE))

M_sID_0 <- sapply(M_fname, function(x) {strsplit(x, split = "_S0")[[1]][1]})
M_sID_1 <- sapply(M_sID_0, function(x) {strsplit(x, split = "_S1")[[1]][1]})
M_sID_2 <- sapply(M_sID_1, function(x) {strsplit(x, split = "_S2")[[1]][1]})
M_sID_3 <- sapply(M_sID_2, function(x) {strsplit(x, split = "_S3")[[1]][1]})
M_sID_4 <- sapply(M_sID_3, function(x) {strsplit(x, split = "_S4")[[1]][1]})
M_sID_5 <- sapply(M_sID_4, function(x) {strsplit(x, split = "_S5")[[1]][1]})
M_sID_6 <- sapply(M_sID_5, function(x) {strsplit(x, split = "_S6")[[1]][1]})
M_sID_7 <- sapply(M_sID_6, function(x) {strsplit(x, split = "_S7")[[1]][1]})
M_sID_8 <- sapply(M_sID_7, function(x) {strsplit(x, split = "_S8")[[1]][1]})
M_sID_9 <- sapply(M_sID_8, function(x) {strsplit(x, split = "_S9")[[1]][1]})

# for column absolute-filepath:
# obtain file names with extensions and prepend names with filepath
M_afp <- gsub(pattern = M_ifp_fsq, replacement = "$PWD/raw_fastq",
              list.files(path = M_ifp_fsq, full.names = TRUE))

# for column direction:
# match pattern in file names to specify read direction
M_drc <- sapply(M_afp, function(x) {
  ifelse(grepl(pattern = "R1", x, ignore.case = FALSE),
         yes = "forward", no = "reverse")})

### ************************************
### M - STEP 2 - create manifest ----
### ************************************

# create manifest df containing info from the vetors created above
M_mft_0 <- data.frame("sample_id" = M_sID_9, 
                      "absolute_filepath" = M_afp,
                      "direction" = M_drc)

# logical check to confirm manifest SampleIDs match meta SampleIDs:
# (1) remove duplicate rows from manifest and create a sorted vector for check
M_mft_unq_0 <- data.frame("SampleID" = unique(M_mft_0$sample_id))
M_mft_unq_1 <- as.character(sort(M_mft_unq_0$SampleID))

# (2) read in metadata file and create a sorted vector for check
M_met_0 <- read.table(ifp_met, header = TRUE, sep = "\t", as.is = TRUE)
M_met_1 <- as.character(sort(M_met_0$SampleID))

# (3) store logical check and print to confirm check = TRUE
M_sID_LOGICAL <- all(M_mft_unq_1 == M_met_1)

# format manifest prior to output
M_mft_1 <- M_mft_0
names(M_mft_1) <- gsub(pattern = "\\_", replacement = "-", x = names(M_mft_1))

### ************************************
### M - WRITE OUTPUTS ----
### ************************************

# print any LOGICAL
print(M_sID_LOGICAL)

write.csv(row.names = FALSE, quote = FALSE, file = M_ofv_mft, x = M_mft_1)
write.csv(row.names = FALSE, quote = FALSE, file = ofp_mft, x = M_mft_1)

save.image(file = M_ofv_WS)

### ************************************
### SECTION Q - QIIME2 analysis ----
### ************************************

# purpose:
# this section stores all QIIME 2 commands
# the R Studio terminal is capable of running these commands;
# however, I did not set the script up in that way 
# rather, I am partial to the R Studio interface
# and to scripts written in the .R format
# so, *gasp* copy and paste is the approach here
# if viewing this script in R Studio, I would recommend unchecking the box here:
# Tools --> Global Options --> Code --> Diagnostics --> Show diagnostics for R

# computational details
# Model: MacBook Pro (15-inch, 2017)
# Processor: 2.8 GHz Intel Core i7 (quad core)
# Memory: 16 GB 2133 MHz LPDDR3
# OS: Mojave 10.14.2

# QIIME 2 version = 2018.11
# run with raw fastq single-end forward (R1) reads

### ************************************
### Q - STEP 1 - import raw fastq seqs ----
### ************************************

# import FASTQ sequences into QIIME 2
# runtime: ~2 minutes
qiime tools import \
--input-path manifest_R1.csv \
--type 'SampleData[SequencesWithQuality]' \
--input-format SingleEndFastqManifestPhred33 \
--output-path seqs_R1.qza

# visualize summary statistics
# runtime: ~2 minutes
qiime demux summarize \
--i-data seqs_R1.qza \
--o-visualization seqs_R1.qzv

# view
qiime tools view \
seqs_R1.qzv

# quit viewer
q

### ************************************
### Q - STEP 2 - remove adapter sequences ----
### ************************************

# remove 515F (Parada) primer and any preceding bases
# runtime: ~4 minutes
qiime cutadapt trim-single \
--i-demultiplexed-sequences seqs_R1.qza \
--p-cores 3 \
--p-front GTGYCAGCMGCCGCGGTAA \
--p-error-rate 0.1 \
--p-match-adapter-wildcards \
--o-trimmed-sequences seqs_R1_trim.qza

# visualize summary statistics
# runtime: ~2 minutes
qiime demux summarize \
--i-data seqs_R1_trim.qza \
--o-visualization seqs_R1_trim.qzv

# view
qiime tools view \
seqs_R1_trim.qzv

# quit viewer
q

### ************************************
### Q - STEP 3 - dada2 ----
### ************************************

# make new directory to store outputs
mkdir dada2/

# runtime: ~11 hours
# change from default: --p-n-reads-learn 1000000 --> 10000000
qiime dada2 denoise-single \
--p-n-threads 0 \
--i-demultiplexed-seqs seqs_R1_trim.qza \
--p-trunc-len 248 \
--p-trim-left 0 \
--p-max-ee 2 \
--p-trunc-q 2 \
--p-chimera-method consensus \
--p-min-fold-parent-over-abundance 1 \
--p-n-reads-learn 10000000 \
--o-table dada2/fts_tab.qza \
--o-representative-sequences dada2/rep_seq.qza \
--o-denoising-stats dada2/stats.qza

# visualize summary statistics
qiime feature-table summarize \
--i-table dada2/fts_tab.qza \
--o-visualization dada2/fts_tab.qzv

# view
qiime tools view \
dada2/fts_tab.qzv

# quit viewer
q

### ************************************
### Q - STEP 4 - tabulate rep seqs ----
### ************************************

# tabulate seqs for easy 'BLASTing'
qiime feature-table tabulate-seqs \
--i-data dada2/rep_seq.qza \
--o-visualization dada2/rep_seq.qzv

# view
qiime tools view \
dada2/rep_seq.qzv

# quit viewer
q

### ************************************
### Q - STEP 5 - train 'in-house' taxonomy classifiers ----
### ************************************

# NOTE: these classifiers use the 515F/806R (Parada/Apprill) primer pair

# make new directory to store outputs
mkdir class/
  
# import relevant files
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path ../../gg_13_8_otus/rep_set/99_otus.fasta \
--output-path class/ggs_otu_99.qza

qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path ../../SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna \
--output-path class/slv_otu_99.qza

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path ../../gg_13_8_otus/taxonomy/99_otu_taxonomy.txt \
--output-path class/ggs_tax_99.qza

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path ../../SILVA_132_QIIME_release/taxonomy/16S_only/99/taxonomy_7_levels.txt \
--output-path class/slv_tax_99.qza

# extract reference reads
# runtime: ~11 minutes
qiime feature-classifier extract-reads \
--i-sequences class/ggs_otu_99.qza \
--p-f-primer GTGYCAGCMGCCGCGGTAA \
--p-r-primer GGACTACNVGGGTWTCTAAT \
--p-trunc-len 248 \
--p-min-length 50 \
--p-max-length 0 \
--o-reads class/ggs_13_8_99_515Par_806App_248bp_ref_seqs.qza

# runtime: ~21 minutes
qiime feature-classifier extract-reads \
--i-sequences class/slv_otu_99.qza \
--p-f-primer GTGYCAGCMGCCGCGGTAA \
--p-r-primer GGACTACNVGGGTWTCTAAT \
--p-trunc-len 248 \
--p-min-length 50 \
--p-max-length 0 \
--o-reads class/slv_13_8_99_515Par_806App_248bp_ref_seqs.qza

# train the classifiers
# runtime: ~4 minutes
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads class/ggs_13_8_99_515Par_806App_248bp_ref_seqs.qza \
--i-reference-taxonomy class/ggs_tax_99.qza \
--o-classifier class/ggs_13_8_99_515Par_806App_248bp_nb_classifier.qza

# runtime: ~2 hours
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads class/slv_13_8_99_515Par_806App_248bp_ref_seqs.qza \
--i-reference-taxonomy class/slv_tax_99.qza \
--o-classifier class/slv_13_8_99_515Par_806App_248bp_nb_classifier.qza

### ************************************
### Q - STEP 6 - assign taxonomy ----
### ************************************

# make new directory to store outputs
mkdir taxa/

# assign taxonomy with in-house (int) classifiers
# runtime: ~1 minute
qiime feature-classifier classify-sklearn \
--p-n-jobs -1 \
--p-confidence 0.75 \
--i-reads dada2/rep_seq.qza \
--i-classifier class/ggs_13_8_99_515Par_806App_248bp_nb_classifier.qza \
--o-classification taxa/int_ggs.qza

# runtime: ~21 minutes
qiime feature-classifier classify-sklearn \
--p-n-jobs -1 \
--p-confidence 0.75 \
--i-reads dada2/rep_seq.qza \
--i-classifier class/slv_13_8_99_515Par_806App_248bp_nb_classifier.qza \
--o-classification taxa/int_slv.qza

### ************************************
### Q - EXPORT ARTIFACTS ----
### ************************************

# NOTE: all of the exported files will be copied to the "vault/"

# feature table
qiime tools export \
--input-path dada2/fts_tab.qza \
--output-path dada2/
  
# convert from biom format to classic format (tab-delimited) for use with R
biom convert -i dada2/feature-table.biom -o dada2/fts_tab.tsv --to-tsv

# use preferred text editor to:
# remove the first line '# Constructed from biom file'
# and change '#OTU ID' --> 'FeatureID' (one word, no hash)
vim dada2/fts_tab.tsv

# copy feature table to the "vault/"
cp dada2/fts_tab.tsv vault/q2_fts_tab.tsv

# rep sequences
qiime tools export \
--input-path dada2/rep_seq.qza \
--output-path dada2/
  
# copy representative sequences to the "vault/" (and rename)
cp dada2/dna-sequences.fasta vault/q2_rep_seq.fasta

# in-house trained greengenes taxonomy table
qiime tools export \
--input-path taxa/int_ggs.qza \
--output-path taxa/
  
# use preferred text editor to:
# change 'Feature ID' to 'FeatureID' (one word, no space)
vim taxa/taxonomy.tsv

# copy to rename, then copy greengenes taxonomy table to the "vault/"
cp taxa/taxonomy.tsv taxa/int_ggs.tsv
cp taxa/int_ggs.tsv vault/q2_int_ggs.tsv

# in-house trained silva taxonomy table
qiime tools export \
--input-path taxa/int_slv.qza \
--output-path taxa/
  
# use preferred text editor to:
# change 'Feature ID' to 'FeatureID' (one word, no space)
vim taxa/taxonomy.tsv

# copy to rename, then copy silva taxonomy table to the "vault/"
cp taxa/taxonomy.tsv taxa/int_slv.tsv
cp taxa/int_slv.tsv vault/q2_int_slv.tsv

# proceed to Code_S2.R
