### ************************************
### PREFACE ----
### ************************************

# this file houses the entirety of code for analysis performed in R*
# *after analysis performed in QIIME 2 (see Code_S1.R)

# I would not recommend hitting the source button for this file...
# ... although if everything is setup per the instructions below...
# ... it will run all the way through... took about 10 minutes on my computer

# path from home directory on my machine
# ~/Desktop/MerlinsManuscript/

# *gasp*
setwd("~/Desktop/")

# relative paths from wd for inputs
ifp_met <- "MerlinsManuscript/MetadataFile_S1.txt"
A_ifp_fts_tab <- "MerlinsManuscript/vault/q2_fts_tab.tsv"
A_ifp_rep_seq <- "MerlinsManuscript/vault/q2_rep_seq.fasta"
A_ifp_int_ggs <- "MerlinsManuscript/vault/q2_int_ggs.tsv"
A_ifp_int_slv <- "MerlinsManuscript/vault/q2_int_slv.tsv"

# relative paths from wd for outputs headed for storage in "main/"
ofm_hth_EDPf_pca_plot <- "MerlinsManuscript/main/Figure_3.pdf"
ofm_hth_EDP_RvF_fts_plot <- "MerlinsManuscript/main/Figure_4.pdf"
ofm_study_timeline <- "MerlinsManuscript/main/Figure_1.pdf"

# relative paths from wd for outputs headed for storage in "supplement/"
ofs_ebtks_abs <- "MerlinsManuscript/supplement/DataFile_S1.txt"
ofs_hth_abs_nEDPf <- "MerlinsManuscript/supplement/DataFile_S2.txt"
ofs_hth_abs_EDPf <- "MerlinsManuscript/supplement/DataFile_S3.txt"
ofs_hth_abs_EDP <- "MerlinsManuscript/supplement/DataFile_S4.txt"
ofs_hth_bif_plot <- "MerlinsManuscript/supplement/Figure_S1.pdf"
ofs_hth_nEDPf_pca_plot <- "MerlinsManuscript/supplement/Figure_S2.pdf"
ofs_D_EDP_q1_sigs_tab <- "MerlinsManuscript/supplement/Table_S3.txt"
ofs_hth_EDP_CvR_CvF_fts_plot <- "MerlinsManuscript/supplement/Figure_S3.pdf"

# relative paths for outputs headed for storage in the "vault/"
A_ofv_ebtks_abs <- "MerlinsManuscript/vault/table_ebtks_abs.txt"
A_ofv_hth_abs_nEDPf <- "MerlinsManuscript/vault/table_hth_abs_nEDPf.txt"
A_ofv_hth_abs_EDPf <- "MerlinsManuscript/vault/table_hth_abs_EDPf.txt"
A_ofv_hth_abs_EDP <- "MerlinsManuscript/vault/table_hth_abs_EDP.txt"
B_ofv_hth_bif_plot <- "MerlinsManuscript/vault/plot_bifido.pdf"
C_ofv_hth_EDPf_pca_plot <- "MerlinsManuscript/vault/plot_EDPf_pca.pdf"
C_ofv_hth_nEDPf_pca_plot <- "MerlinsManuscript/vault/plot_nEDPf_pca.pdf"
D_ofv_D_EDP_q1_sigs <- "MerlinsManuscript/vault/table_EDP_q1_sigs.txt"
E_ofv_hth_EDP_RvF_fts <- "MerlinsManuscript/vault/plot_EDP_RvF_fts.pdf"
E_ofv_hth_EDP_CvR_CvF_fts <- "MerlinsManuscript/vault/plot_EDP_CvR_CvF_fts.pdf"
F_ofv_timeline <- "MerlinsManuscript/vault/plot_timeline.pdf"
A_ofv_WS <- "MerlinsManuscript/vault/WS_Section_A.RData"
B_ofv_WS <- "MerlinsManuscript/vault/WS_Section_B.RData"
C_ofv_WS <- "MerlinsManuscript/vault/WS_Section_C.RData"
D_ofv_WS <- "MerlinsManuscript/vault/WS_Section_D.RData"
E_ofv_WS <- "MerlinsManuscript/vault/WS_Section_E.RData"
F_ofv_WS <- "MerlinsManuscript/vault/WS_Section_F.RData"
ofv_WS_all <- "MerlinsManuscript/vault/WS_Final_Code_S2.Rdata"

### ************************************
### PACKAGE INFO ----
### ************************************

# dependencies accessed via pkg::name
# dplyr
# ggbiplot
# ggpubr
# reshape2
# zCompositions

# dependencies access via require
require(ALDEx2)
require(ggplot2)
require(grid)

# version info
# R 3.5.2 'Eggshell Igloo' (2018-12-20)
# R studio 1.1.463
# ALDEx2 1.14.0
# dplyr 0.7.8
# ggbiplot 0.55
# ggplot2 3.1.0
# ggpubr 0.2
# grid 3.5.2
# reshape2 1.4.3
# zCompositions 1.1.2

### ************************************
### SECTION A - create master table ----
### ************************************

# purpose:
# take outputs from QIIME 2 and combine them into a master table...
# ... to serve as the entry point for all downstream processing

# description for each step found under each subsection

# description:
# found under each subsection

# main/supplementary items generated:
# Data File S1
# Data File S2
# Data File S3
# Data File S4

### ************************************
### A - INPUTS - provenance ----
### ************************************

# inputs accessed via: "$PWD/MerlinsManuscript/vault/"
# q2_fts_tab.tsv
# q2_rep_seq.fasta
# q2_int_ggs.tsv
# q2_int_slv.tsv

# provenance for steps prior to DADA2 feature table construction
# seqs.qza -> seqs_trim.qza (listed as X)
# qiime tools import -> qiime cutadapt-trim single (listed as Y)
# import raw FASTQ reads -> remove adapters on 5' end (listed as Z)

# provenance for q2_fts_tab.tsv (DADA2 feature table)
# X -> fts_tab.qza -> feature_table.biom -> feature_table.tsv -> q2_fts_tab.tsv
# Y -> denoise with DADA2 -> export -> convert to tsv -> copy to rename
# Z -> qiime dada2 denoise-single -> qiime tools export -> biom convert -> cp

# provenance for q2_rep_seq.fasta (DADA2 representative sequences)
# X -> rep_seq.qza -> dna-sequences.fasta -> q2_rep_seq.fasta
# Y -> denoise with DADA2 -> export -> copy to rename
# Z -> qiime dada2 denoise-single -> qiime tools export -> cp

# provenance for q2_int_ggs.tsv (Greengenes taxonomic assignments)
# int_ggs.qza -> taxonomy.tsv -> q2_int_ggs.tsv
# assign taxonomy to rep_seqs.qza -> export -> copy to rename
# qiime feature-classifier classify-sklearn -> qiime tools export -> cp

# provenance for q2_int_slv.tsv (SILVA taxonomic assignments)
# int_slv.qza -> taxonomy.tsv -> q2_int_slv.tsv
# assign taxonomy to rep_seqs.qza -> export -> copy to rename
# qiime feature-classifier classify-sklearn -> qiime tools export -> cp

### ************************************
### A - FUNCTIONS ----
### ************************************

# function takes an input .fasta file in two-line format which looks like:
# line 1: >header
# line 2: DNA sequence
# ... and creates a two column dataframe
fasta_to_df <- function(fasta_file, hdr_colname = "", seq_colname = "") {
  # read in the .fasta file line by line
  fasta <- readLines(fasta_file)
  # identify header lines by finding line numbers with the '>' character
  hlines <- grep(pattern = ">", x = fasta)
  # create a three column df containing:
  # header line numbers (hdr)
  # sequence line number beginning (beg)
  # sequence line number end (end)
  slines <- data.frame(hdr = hlines,
                       beg = hlines+1,
                       end = c((hlines-1)[-1], length(fasta)))
  # create a vector of identical length to hdr_lines to be used for
  # storing values obtained in the loop below
  storage_vec <- rep(NA, length(hlines))
  # loop to obtain sequences
  for(i in 1:length(hlines)) {
    storage_vec[i] <- paste(fasta[slines$beg[i]:slines$end[i]], collapse = "")
  }
  # create a two column df containing:
  # the header with '>' character replaced (hdr)
  # representative sequence associated with the header (seq)
  new_df_v0 <- data.frame("hdr" = gsub(pattern = ">", replacement = "",
                                       x = fasta[hlines]),
                          "seq" = storage_vec)
  # replace generic column names with names specified in the function input
  new_df_v1 <- new_df_v0
  hdr_colnum <- which(names(new_df_v1) == "hdr")
  seq_colnum <- which(names(new_df_v1) == "seq")
  names(new_df_v1)[hdr_colnum] <- hdr_colname
  names(new_df_v1)[seq_colnum] <- seq_colname
  return(new_df_v1)
}
#
# example usage:
#df_fasta <- fasta_to_df(fasta_file = ifp_fasta,
#                        hdr_colname = "FeatureID",
#                        seq_colname = "RepSeq")

## *********************************** ## *********************************** ##
## *********************************** ## *********************************** ##

# trunc_tax() function requires an input df with a column containing either:
# greengenes formatted taxonomic strings, e.g.:
# k__; p__ ; c__; o__; f__; g__; s__"
# or silva taxonomic strings, e.g.:
# "D_0__;D_1__;D_2__;D_3__;D_4__;D_5__;D_6__"
# the output df will contain a new column with:
# the desired taxonomic level for the specified database naming convention
# NOTE: if the above strings are not in consecutive order...
# ... the output may not be correct i.e.:
# if string input is missing D_2__ = "D_0__;D_1__Example;D_3__Example;D_4__"
# output for ;D_1__ = "Example;D_3__Example;D_4__" rather than just "Example"
# output for ;D_2__ = "Unassigned"
# output for ;D_3__ = "Example"
trunc_tax <- function(data = data.frame, tax_type = c("Greengenes", "SILVA"),
                      tax_lvl = c("Kingdom", "Phylum", "Class", "Order",
                                  "Family", "Genus", "Species"),
                      input_col = "", str1 = "", str2 = "", classifier = "") {
  
  # internal checks to ensure correct input classes
  if (!inherits(data, "data.frame")) {
    stop("input data must be class 'data.frame'")
  }
  if (!tax_type == "Greengenes" && !tax_type == "SILVA") {
    stop("tax_type should be either 'Greengenes' or 'SILVA'")
  }
  if (!tax_lvl == "Kingdom" && !tax_lvl == "Phylum" && !tax_lvl == "Class" &&
      !tax_lvl == "Order" && !tax_lvl == "Family" && !tax_lvl == "Genus" &&
      !tax_lvl == "Species") {
    stop("tax_lvl should be one of: 'Kingdom', 'Phylum', 'Class', 'Order',
         'Family', 'Genus', 'Species'")
  }
  if (!inherits(input_col, "character")) {
    stop("input_col must be character")
  }
  if (!inherits(str1, "character")) {
    stop("str1 must be character")
  }
  if (!inherits(str2, "character")) {
    stop("str2 must be character")
  }
  if (!inherits(classifier, "character")) {
    stop("classifier must be character")
  }
  
  # if tax_lvl is Species, pasting the string split is handled differently
  num <- ifelse(tax_lvl == "Species", yes = 2, no = 1)
  
  # store original input df and format the new df
  new_df1 <- data
  new_df1$new_col <- 1
  
  # loop through data and truncate the taxonomic strings
  for(row in 1:nrow(new_df1)) {
    # split the string in the specified input column with input of str1
    splt1 <- strsplit(new_df1[row, input_col], str1)[[1]][2]
    # goal of first ifelse test:
    # yes = tax level was not assigned
    # no = tax level was likely assigned, split the string with input of str2
    splt2 <- ifelse(is.na(splt1), yes = "Unassigned",
                    no = paste("", strsplit(splt1, str2)[[1]][num], sep = ""))
    # goal of second ifelse test:
    # yes = tax level was not assigned
    # no = tax level was actually assigned
    new_df1[row, "new_col"] <- ifelse(nchar(splt2) == 0 | "NA" %in% splt2,
                                      yes = "Unassigned", no = splt2)
    # if tax_type is Greengenes, tax_lvl is Species and row is not Unassigned,
    # append Genus in front of Species assignment
    if (tax_type == "Greengenes" & tax_lvl == "Species" &
        !new_df1[row, "new_col"] == "Unassigned") {
      new_df1[row, "new_col"] <- paste("", strsplit(splt1, str2)[[1]][1],
                                       strsplit(splt1, str2)[[1]][2], sep = " ")
    }
  } # close loop
  
  # internal warning in the event that all rows in the new column are Unassigned
  # this outcome could be correct if the specified tax lvl had no assignemnts;
  # however, if assignments are expected for the specified tax lvl, this outcome
  # is incorrect and the function inputs need to be double checked
  if (isTRUE(all(grepl(pattern = "Unassigned", x = new_df1[, "new_col"],
                       ignore.case = FALSE, fixed = TRUE)))) {
    warning("in new_col all row values in output column are 'Unassigned'")
  }
  
  # format the new df
  new_df2 <- new_df1
  new_col_num <- which(names(new_df2) == 'new_col')
  if (tax_type == "Greengenes") {
    names(new_df2)[new_col_num] <- paste(classifier, "ggs", tax_lvl, sep = ".")
  }
  if (tax_type == "SILVA") {
    names(new_df2)[new_col_num] <- paste(classifier, "slv", tax_lvl, sep = ".")
  }
  return(new_df2)
  }

# inputs are incorrect, prints warning
#new_df <- trunc_tax(data = a.df, tax_type = "silva",
#                   tax_lvl = "Species", input_col = "ggTaxon",
#                   str1 = "; g__", str2 = "; s__")

# corrected inputs, no warning
#new_df <- trunc_tax(data = a.df, tax_type = "greengenes",
#                   tax_lvl = "Species", input_col = "ggTaxon",
#                   str1 = "; g__", str2 = "; s__")

### ************************************
### A - STEP 1 - format taxa tables ----
### ************************************

# this step truncates Greengenes & SILVA taxonomic lineages, placing them...
# ...into new cols: .Kingdom; .Phylum; .Class; .Order; .Family; .Genus; .Species
# this step also determines the lowest assigned level, placing that info...
# ... into new cols: .lws.txn and .lws.lvl

# read in taxonomy tables
A_int_ggs <- read.table(A_ifp_int_ggs, header = TRUE, sep = "\t", as.is = TRUE,
                        col.names = c("FeatureID", "int.ggs.tax",
                                      "int.ggs.cnf"))
A_int_slv <- read.table(A_ifp_int_slv, header = TRUE, sep = "\t", as.is = TRUE,
                        col.names = c("FeatureID", "int.slv.tax",
                                      "int.slv.cnf"))

# truncate taxonomic lineages for each level, and reorder columns
A_int_ggs_L1 <- trunc_tax(data = A_int_ggs, tax_type = "Greengenes",
                          tax_lvl = "Kingdom", input_col = "int.ggs.tax",
                          str1 = "k__", str2 = "; p__", classifier = "int")
A_int_ggs_L2 <- trunc_tax(data = A_int_ggs_L1, tax_type = "Greengenes",
                          tax_lvl = "Phylum", input_col = "int.ggs.tax",
                          str1 = "; p__", str2 = "; c__", classifier = "int")
A_int_ggs_L3 <- trunc_tax(data = A_int_ggs_L2, tax_type = "Greengenes",
                          tax_lvl = "Class", input_col = "int.ggs.tax",
                          str1 = "; c__", str2 = "; o__", classifier = "int")
A_int_ggs_L4 <- trunc_tax(data = A_int_ggs_L3, tax_type = "Greengenes",
                          tax_lvl = "Order", input_col = "int.ggs.tax",
                          str1 = "; o__", str2 = "; f__", classifier = "int")
A_int_ggs_L5 <- trunc_tax(data = A_int_ggs_L4, tax_type = "Greengenes",
                          tax_lvl = "Family", input_col = "int.ggs.tax",
                          str1 = "; f__", str2 = "; g__", classifier = "int")
A_int_ggs_L6 <- trunc_tax(data = A_int_ggs_L5, tax_type = "Greengenes",
                          tax_lvl = "Genus", input_col = "int.ggs.tax",
                          str1 = "; g__", str2 = "; s__", classifier = "int")
A_int_ggs_L7 <- trunc_tax(data = A_int_ggs_L6, tax_type = "Greengenes",
                          tax_lvl = "Species", input_col = "int.ggs.tax",
                          str1 = "; g__", str2 = "; s__", classifier = "int")
A_int_ggs_trunc <- dplyr::select(A_int_ggs_L7, FeatureID, int.ggs.tax,
                                 int.ggs.Kingdom, int.ggs.Phylum, int.ggs.Class,
                                 int.ggs.Order, int.ggs.Family, int.ggs.Genus,
                                 int.ggs.Species, int.ggs.cnf)

A_int_slv_L1 <- trunc_tax(data = A_int_slv, tax_type = "SILVA",
                          tax_lvl = "Kingdom", input_col = "int.slv.tax",
                          str1 = "D_0__", str2 = ";D_1__", classifier = "int")
A_int_slv_L2 <- trunc_tax(data = A_int_slv_L1, tax_type = "SILVA",
                          tax_lvl = "Phylum", input_col = "int.slv.tax",
                          str1 = ";D_1__", str2 = ";D_2__", classifier = "int")
A_int_slv_L3 <- trunc_tax(data = A_int_slv_L2, tax_type = "SILVA",
                          tax_lvl = "Class", input_col = "int.slv.tax",
                          str1 = ";D_2__", str2 = ";D_3__", classifier = "int")
A_int_slv_L4 <- trunc_tax(data = A_int_slv_L3, tax_type = "SILVA",
                          tax_lvl = "Order", input_col = "int.slv.tax",
                          str1 = ";D_3__", str2 = ";D_4__", classifier = "int")
A_int_slv_L5 <- trunc_tax(data = A_int_slv_L4, tax_type = "SILVA",
                          tax_lvl = "Family", input_col = "int.slv.tax",
                          str1 = ";D_4__", str2 = ";D_5__", classifier = "int")
A_int_slv_L6 <- trunc_tax(data = A_int_slv_L5, tax_type = "SILVA",
                          tax_lvl = "Genus", input_col = "int.slv.tax",
                          str1 = ";D_5__", str2 = ";D_6__", classifier = "int")
A_int_slv_L7 <- trunc_tax(data = A_int_slv_L6, tax_type = "SILVA",
                          tax_lvl = "Species", input_col = "int.slv.tax",
                          str1 = ";D_5__", str2 = ";D_6__", classifier = "int")
A_int_slv_trunc <- dplyr::select(A_int_slv_L7, FeatureID, int.slv.tax,
                                 int.slv.Kingdom, int.slv.Phylum, int.slv.Class,
                                 int.slv.Order, int.slv.Family, int.slv.Genus,
                                 int.slv.Species, int.slv.cnf)

# create vectors useful to determine lowest assigned level for each feature
A_int_ggs_col <- c("int.ggs.Kingdom", "int.ggs.Phylum", "int.ggs.Class",
                   "int.ggs.Order", "int.ggs.Family","int.ggs.Genus",
                   "int.ggs.Species")

A_int_slv_col <- c("int.slv.Kingdom", "int.slv.Phylum", "int.slv.Class",
                   "int.slv.Order", "int.slv.Family","int.slv.Genus",
                   "int.slv.Species")

A_int_ggs_str <- "int.ggs."
A_int_slv_str <- "int.slv."

A_val_una <- "Unassigned"

A_flt_col <- "FeatureID"

# SILVA has assignments that are not of interest
# i.e. Family = uncultured Rubrobacteria bacterium, among others
# these tend to begin at the Family level and proceed into lower levels
# however, we will play it safe and assume that this could happen at any level
# our goal here is to create a vector of unwanted values
# this is akin to the vector A_val_una we defined above for Greengenes

# create vectors of unique values from rows in all of the truncated level cols
A_int_slv_K_unq <- unique(A_int_slv_trunc[, A_int_slv_col[1]])
A_int_slv_P_unq <- unique(A_int_slv_trunc[, A_int_slv_col[2]])
A_int_slv_C_unq <- unique(A_int_slv_trunc[, A_int_slv_col[3]])
A_int_slv_O_unq <- unique(A_int_slv_trunc[, A_int_slv_col[4]])
A_int_slv_F_unq <- unique(A_int_slv_trunc[, A_int_slv_col[5]])
A_int_slv_G_unq <- unique(A_int_slv_trunc[, A_int_slv_col[6]])
A_int_slv_S_unq <- unique(A_int_slv_trunc[, A_int_slv_col[7]])

# define a vector of values we do not want with partial strings of:
# unassigned - uncultured - unidentified - gut metagenome/human gut metagenome
# a few of the uncultured organisms were found to give issue so...
# ... they are included with full strings
A_int_slv_pat <- paste(c("unass", "uncul", "unide", "gut",
                         "Clostridium sp. Clone-49",
                         "Clostridium sp. Culture-54",
                         "Clostridium sp. Culture-41",
                         "Clostridium sp. Culture-27",
                         "Clostridium phoceensis",
                         "Azospirillum sp. 47_25"), collapse = "|")

# grep and return vectors with the unwanted character values
A_int_slv_K_unq_unw <- grep(pattern = A_int_slv_pat, x = A_int_slv_K_unq,
                            ignore.case = TRUE, value = TRUE)
A_int_slv_K_unq_unw <- grep(pattern = A_int_slv_pat, x = A_int_slv_K_unq,
                            ignore.case = TRUE, value = TRUE)
A_int_slv_P_unq_unw <- grep(pattern = A_int_slv_pat, x = A_int_slv_P_unq,
                            ignore.case = TRUE, value = TRUE)
A_int_slv_C_unq_unw <- grep(pattern = A_int_slv_pat, x = A_int_slv_C_unq,
                            ignore.case = TRUE, value = TRUE)
A_int_slv_O_unq_unw <- grep(pattern = A_int_slv_pat, x = A_int_slv_O_unq,
                            ignore.case = TRUE, value = TRUE)
A_int_slv_F_unq_unw <- grep(pattern = A_int_slv_pat, x = A_int_slv_F_unq,
                            ignore.case = TRUE, value = TRUE)
A_int_slv_G_unq_unw <- grep(pattern = A_int_slv_pat, x = A_int_slv_G_unq,
                            ignore.case = TRUE, value = TRUE)
A_int_slv_S_unq_unw <- grep(pattern = A_int_slv_pat, x = A_int_slv_S_unq,
                            ignore.case = TRUE, value = TRUE)

# combine the vectors, passing unique again, and also collapsing everything ...
# .... into a vector for use in grepl pattern matching
A_val_slv <- paste(unique(c(A_int_slv_K_unq_unw, A_int_slv_P_unq_unw,
                            A_int_slv_C_unq_unw, A_int_slv_O_unq_unw,
                            A_int_slv_F_unq_unw, A_int_slv_G_unq_unw,
                            A_int_slv_S_unq_unw)), collapse = "|")

# now that we have our vectors, process each of the truncated level cols...
# ... to find the lowest assignment level

# Unassigned Kingdom:
A_int_ggs_unK_1 <- dplyr::filter_at(A_int_ggs_trunc,
                                    dplyr::vars(A_int_ggs_col),
                                    dplyr::all_vars(. == A_val_una))
A_int_slv_unK_1 <- dplyr::filter_at(A_int_slv_trunc,
                                    dplyr::vars(A_int_slv_col),
                                    dplyr::all_vars(
                                      grepl(pattern = A_val_slv, x = .)))

# Assigned Kingdom; Unassigned Phylum:
A_int_ggs_K_1 <- dplyr::inner_join(
  dplyr::filter_at(A_int_ggs_trunc, dplyr::vars(A_int_ggs_col[1]),
                   dplyr::all_vars(!. == A_val_una)),
  dplyr::filter_at(A_int_ggs_trunc, dplyr::vars(A_int_ggs_col[2:7]),
                   dplyr::all_vars(. == A_val_una)))

A_int_slv_K_1 <- dplyr::inner_join(
  dplyr::filter_at(A_int_slv_trunc, dplyr::vars(A_int_slv_col[1]),
                   dplyr::all_vars(
                     !grepl(pattern = A_val_slv, x = .))),
  dplyr::filter_at(A_int_slv_trunc, dplyr::vars(A_int_slv_col[2:7]),
                   dplyr::all_vars(
                     grepl(pattern = A_val_slv, x = .))))

# Assigned Kingdom-Phylum; Unassigned Class:
A_int_ggs_P_1 <- dplyr::inner_join(
  dplyr::filter_at(A_int_ggs_trunc, dplyr::vars(A_int_ggs_col[1:2]),
                   dplyr::all_vars(!. == A_val_una)),
  dplyr::filter_at(A_int_ggs_trunc, dplyr::vars(A_int_ggs_col[3:7]),
                   dplyr::all_vars(. == A_val_una)))

A_int_slv_P_1 <- dplyr::inner_join(
  dplyr::filter_at(A_int_slv_trunc, dplyr::vars(A_int_slv_col[1:2]),
                   dplyr::all_vars(
                     !grepl(pattern = A_val_slv, x = .))),
  dplyr::filter_at(A_int_slv_trunc, dplyr::vars(A_int_slv_col[3:7]),
                   dplyr::all_vars(
                     grepl(pattern = A_val_slv, x = .))))

# Assigned Kingdom-Phylum-Class; Unassigned Order:
A_int_ggs_C_1 <- dplyr::inner_join(
  dplyr::filter_at(A_int_ggs_trunc, dplyr::vars(A_int_ggs_col[1:3]),
                   dplyr::all_vars(!. == A_val_una)),
  dplyr::filter_at(A_int_ggs_trunc, dplyr::vars(A_int_ggs_col[4:7]),
                   dplyr::all_vars(. == A_val_una)))

A_int_slv_C_1 <- dplyr::inner_join(
  dplyr::filter_at(A_int_slv_trunc, dplyr::vars(A_int_slv_col[1:3]),
                   dplyr::all_vars(
                     !grepl(pattern = A_val_slv, x = .))),
  dplyr::filter_at(A_int_slv_trunc, dplyr::vars(A_int_slv_col[4:7]),
                   dplyr::all_vars(
                     grepl(pattern = A_val_slv, x = .))))

# Assigned Kingdom-Phylum-Class-Order; Unassigned Family:
A_int_ggs_O_1 <- dplyr::inner_join(
  dplyr::filter_at(A_int_ggs_trunc, dplyr::vars(A_int_ggs_col[1:4]),
                   dplyr::all_vars(!. == A_val_una)),
  dplyr::filter_at(A_int_ggs_trunc, dplyr::vars(A_int_ggs_col[5:7]),
                   dplyr::all_vars(. == A_val_una)))

A_int_slv_O_1 <- dplyr::inner_join(
  dplyr::filter_at(A_int_slv_trunc, dplyr::vars(A_int_slv_col[1:4]),
                   dplyr::all_vars(
                     !grepl(pattern = A_val_slv, x = .))),
  dplyr::filter_at(A_int_slv_trunc, dplyr::vars(A_int_slv_col[5:7]),
                   dplyr::all_vars(
                     grepl(pattern = A_val_slv, x = .))))

# Assigned Kingdom-Phylum-Class-Order-Family; Unassigned Genus:
A_int_ggs_F_1 <- dplyr::inner_join(
  dplyr::filter_at(A_int_ggs_trunc, dplyr::vars(A_int_ggs_col[1:5]),
                   dplyr::all_vars(!. == A_val_una)),
  dplyr::filter_at(A_int_ggs_trunc, dplyr::vars(A_int_ggs_col[6:7]),
                   dplyr::all_vars(. == A_val_una)))

A_int_slv_F_1 <- dplyr::inner_join(
  dplyr::filter_at(A_int_slv_trunc, dplyr::vars(A_int_slv_col[1:5]),
                   dplyr::all_vars(
                     !grepl(pattern = A_val_slv, x = .))),
  dplyr::filter_at(A_int_slv_trunc, dplyr::vars(A_int_slv_col[6:7]),
                   dplyr::all_vars(
                     grepl(pattern = A_val_slv, x = .))))

# Assigned Kingdom-Phylum-Class-Order-Family-Genus; Unassigned Species:
A_int_ggs_G_1 <- dplyr::inner_join(
  dplyr::filter_at(A_int_ggs_trunc, dplyr::vars(A_int_ggs_col[1:6]),
                   dplyr::all_vars(!. == A_val_una)),
  dplyr::filter_at(A_int_ggs_trunc, dplyr::vars(A_int_ggs_col[7]),
                   dplyr::all_vars(. == A_val_una)))

A_int_slv_G_1 <- dplyr::inner_join(
  dplyr::filter_at(A_int_slv_trunc, dplyr::vars(A_int_slv_col[1:6]),
                   dplyr::all_vars(
                     !grepl(pattern = A_val_slv, x = .))),
  dplyr::filter_at(A_int_slv_trunc, dplyr::vars(A_int_slv_col[7]),
                   dplyr::all_vars(
                     grepl(pattern = A_val_slv, x = .))))

# Assigned Kingdom-Phylum-Class-Order-Family-Genus-Species:
A_int_ggs_S_1 <- dplyr::filter_at(A_int_ggs_trunc,
                                  dplyr::vars(A_int_ggs_col[1:7]),
                                  dplyr::all_vars(!. == A_val_una))

A_int_slv_S_1 <- dplyr::filter_at(A_int_slv_trunc,
                                  dplyr::vars(A_int_slv_col[1:7]),
                                  dplyr::all_vars(
                                    !grepl(pattern = A_val_slv, x = .)))

# retain relevant columns
A_int_ggs_unK_2 <- dplyr::select(A_int_ggs_unK_1,
                                 dplyr::one_of(A_flt_col, A_int_ggs_col[1]))
A_int_slv_unK_2 <- dplyr::select(A_int_slv_unK_1,
                                 dplyr::one_of(A_flt_col, A_int_slv_col[1]))

A_int_ggs_K_2 <- dplyr::select(A_int_ggs_K_1,
                               dplyr::one_of(A_flt_col, A_int_ggs_col[1]))
A_int_slv_K_2 <- dplyr::select(A_int_slv_K_1,
                               dplyr::one_of(A_flt_col, A_int_slv_col[1]))

A_int_ggs_P_2 <- dplyr::select(A_int_ggs_P_1,
                               dplyr::one_of(A_flt_col, A_int_ggs_col[2]))
A_int_slv_P_2 <- dplyr::select(A_int_slv_P_1,
                               dplyr::one_of(A_flt_col, A_int_slv_col[2]))

A_int_ggs_C_2 <- dplyr::select(A_int_ggs_C_1,
                               dplyr::one_of(A_flt_col, A_int_ggs_col[3]))
A_int_slv_C_2 <- dplyr::select(A_int_slv_C_1,
                               dplyr::one_of(A_flt_col, A_int_slv_col[3]))

A_int_ggs_O_2 <- dplyr::select(A_int_ggs_O_1,
                               dplyr::one_of(A_flt_col, A_int_ggs_col[4]))
A_int_slv_O_2 <- dplyr::select(A_int_slv_O_1,
                               dplyr::one_of(A_flt_col, A_int_slv_col[4]))

A_int_ggs_F_2 <- dplyr::select(A_int_ggs_F_1,
                               dplyr::one_of(A_flt_col, A_int_ggs_col[5]))
A_int_slv_F_2 <- dplyr::select(A_int_slv_F_1,
                               dplyr::one_of(A_flt_col, A_int_slv_col[5]))

A_int_ggs_G_2 <- dplyr::select(A_int_ggs_G_1,
                               dplyr::one_of(A_flt_col, A_int_ggs_col[6]))
A_int_slv_G_2 <- dplyr::select(A_int_slv_G_1,
                               dplyr::one_of(A_flt_col, A_int_slv_col[6]))

A_int_ggs_S_2 <- dplyr::select(A_int_ggs_S_1,
                               dplyr::one_of(A_flt_col, A_int_ggs_col[7]))
A_int_slv_S_2 <- dplyr::select(A_int_slv_S_1,
                               dplyr::one_of(A_flt_col, A_int_slv_col[7]))

# rename column int.ggs.L to int.ggs.lws.txn
A_int_ggs_unK_3 <- dplyr::rename(A_int_ggs_unK_2,
                                 int.ggs.lws.txn = A_int_ggs_col[1])
A_int_slv_unK_3 <- dplyr::rename(A_int_slv_unK_2,
                                 int.slv.lws.txn = A_int_slv_col[1])

A_int_ggs_K_3 <- dplyr::rename(A_int_ggs_K_2,
                               int.ggs.lws.txn = A_int_ggs_col[1])
A_int_slv_K_3 <- dplyr::rename(A_int_slv_K_2,
                               int.slv.lws.txn = A_int_slv_col[1])

A_int_ggs_P_3 <- dplyr::rename(A_int_ggs_P_2,
                               int.ggs.lws.txn = A_int_ggs_col[2])
A_int_slv_P_3 <- dplyr::rename(A_int_slv_P_2,
                               int.slv.lws.txn = A_int_slv_col[2])

A_int_ggs_C_3 <- dplyr::rename(A_int_ggs_C_2,
                               int.ggs.lws.txn = A_int_ggs_col[3])
A_int_slv_C_3 <- dplyr::rename(A_int_slv_C_2,
                               int.slv.lws.txn = A_int_slv_col[3])

A_int_ggs_O_3 <- dplyr::rename(A_int_ggs_O_2,
                               int.ggs.lws.txn = A_int_ggs_col[4])
A_int_slv_O_3 <- dplyr::rename(A_int_slv_O_2,
                               int.slv.lws.txn = A_int_slv_col[4])

A_int_ggs_F_3 <- dplyr::rename(A_int_ggs_F_2,
                               int.ggs.lws.txn = A_int_ggs_col[5])
A_int_slv_F_3 <- dplyr::rename(A_int_slv_F_2,
                               int.slv.lws.txn = A_int_slv_col[5])

A_int_ggs_G_3 <- dplyr::rename(A_int_ggs_G_2,
                               int.ggs.lws.txn = A_int_ggs_col[6])
A_int_slv_G_3 <- dplyr::rename(A_int_slv_G_2,
                               int.slv.lws.txn = A_int_slv_col[6])

A_int_ggs_S_3 <- dplyr::rename(A_int_ggs_S_2,
                               int.ggs.lws.txn = A_int_ggs_col[7])
A_int_slv_S_3 <- dplyr::rename(A_int_slv_S_2,
                               int.slv.lws.txn = A_int_slv_col[7])

# add new column indicating the assignment level for the lowest assigned taxon
A_int_ggs_unK_4 <- A_int_ggs_unK_3
A_int_slv_unK_4 <- A_int_slv_unK_3

A_int_ggs_K_4 <- A_int_ggs_K_3
A_int_slv_K_4 <- A_int_slv_K_3

A_int_ggs_P_4 <- A_int_ggs_P_3
A_int_slv_P_4 <- A_int_slv_P_3

A_int_ggs_C_4 <- A_int_ggs_C_3
A_int_slv_C_4 <- A_int_slv_C_3

A_int_ggs_O_4 <- A_int_ggs_O_3
A_int_slv_O_4 <- A_int_slv_O_3

A_int_ggs_F_4 <- A_int_ggs_F_3
A_int_slv_F_4 <- A_int_slv_F_3

A_int_ggs_G_4 <- A_int_ggs_G_3
A_int_slv_G_4 <- A_int_slv_G_3

A_int_ggs_S_4 <- A_int_ggs_S_3
A_int_slv_S_4 <- A_int_slv_S_3

A_int_ggs_unK_4$int.ggs.lws.lvl <- A_val_una
A_int_slv_unK_4$int.slv.lws.lvl <- A_val_una

A_int_ggs_K_4$int.ggs.lws.lvl <- strsplit(A_int_ggs_col[1],
                                          A_int_ggs_str)[[1]][2]
A_int_slv_K_4$int.slv.lws.lvl <- strsplit(A_int_slv_col[1],
                                          A_int_slv_str)[[1]][2]

A_int_ggs_P_4$int.ggs.lws.lvl <- strsplit(A_int_ggs_col[2],
                                          A_int_ggs_str)[[1]][2]
A_int_slv_P_4$int.slv.lws.lvl <- strsplit(A_int_slv_col[2],
                                          A_int_slv_str)[[1]][2]

A_int_ggs_C_4$int.ggs.lws.lvl <- strsplit(A_int_ggs_col[3],
                                          A_int_ggs_str)[[1]][2]
A_int_slv_C_4$int.slv.lws.lvl <- strsplit(A_int_slv_col[3],
                                          A_int_slv_str)[[1]][2]

A_int_ggs_O_4$int.ggs.lws.lvl <- strsplit(A_int_ggs_col[4],
                                          A_int_ggs_str)[[1]][2]
A_int_slv_O_4$int.slv.lws.lvl <- strsplit(A_int_slv_col[4],
                                          A_int_slv_str)[[1]][2]

A_int_ggs_F_4$int.ggs.lws.lvl <- strsplit(A_int_ggs_col[5],
                                          A_int_ggs_str)[[1]][2]
A_int_slv_F_4$int.slv.lws.lvl <- strsplit(A_int_slv_col[5],
                                          A_int_slv_str)[[1]][2]

A_int_ggs_G_4$int.ggs.lws.lvl <- strsplit(A_int_ggs_col[6],
                                          A_int_ggs_str)[[1]][2]
A_int_slv_G_4$int.slv.lws.lvl <- strsplit(A_int_slv_col[6],
                                          A_int_slv_str)[[1]][2]

A_int_ggs_S_4$int.ggs.lws.lvl <- strsplit(A_int_ggs_col[7],
                                          A_int_ggs_str)[[1]][2]
A_int_slv_S_4$int.slv.lws.lvl <- strsplit(A_int_slv_col[7],
                                          A_int_slv_str)[[1]][2]

# combine the above dfs
A_int_ggs_lws <- rbind(A_int_ggs_unK_4, A_int_ggs_K_4, A_int_ggs_P_4,
                       A_int_ggs_C_4, A_int_ggs_O_4, A_int_ggs_F_4,
                       A_int_ggs_G_4, A_int_ggs_S_4)
A_int_slv_lws <- rbind(A_int_slv_unK_4, A_int_slv_K_4, A_int_slv_P_4,
                       A_int_slv_C_4, A_int_slv_O_4, A_int_slv_F_4,
                       A_int_slv_G_4, A_int_slv_S_4)

# check to make sure all Features are accounted for in the lowest level dfs
A_int_ggs_test <- dplyr::anti_join(x = A_int_ggs_trunc, y = A_int_ggs_lws,
                                   by = "FeatureID")
A_int_slv_test <- dplyr::anti_join(x = A_int_slv_trunc, y = A_int_slv_lws,
                                   by = "FeatureID")

A_int_ggs_LOGICAL <- ifelse(nrow(A_int_ggs_test) == 0, yes = TRUE, no = FALSE)
A_int_slv_LOGICAL <- ifelse(nrow(A_int_slv_test) == 0, yes = TRUE, no = FALSE)

# merge the lws dfs with the trunc dfs
A_int_ggs_frm <- merge(A_int_ggs_trunc, A_int_ggs_lws, by = "FeatureID",
                       sort = FALSE)
A_int_slv_frm <- merge(A_int_slv_trunc, A_int_slv_lws, by = "FeatureID",
                       sort = FALSE)

# merge the formatted Greengenes and SILVA dfs
A_tax_mrg <- merge(A_int_ggs_frm, A_int_slv_frm, by = "FeatureID",
                   sort = FALSE)

### ************************************
### A - STEP 2 - create master table (ebtks) ----
### ************************************

## this step creates a master dataframe (ebtks_) containing:
# 'FeatureID' = MD5 sums generated via QIIME 2
# 'RepSeq' = representative DNA sequence for respective FeatureID
# 'int.ggs.tax' = taxonomic lineage with Greengenes naming convention
# 'int.ggs.cnf'= confidence in Greengenes taxonomic classification
# 'int.ggs.L' = taxon for respective level (L)
# 'int.ggs.lws.txn' = lowest assignment by Greengenes
# 'int.ggs.lws.lvl' = level for lowest assignment by Greengenes
# 'int.slv.tax' = taxonomic lineage with SILVA naming convention
# 'int.slv.cnf' = confidence in SILVA taxonomic classification
# 'int.slv.L' = taxon for respective level (L)
# 'int.ggs.lws.txn' = lowest assignment by SILVA
# 'int.ggs.lws.lvl' = level for lowest assignment by SILVA
# remaining columns = SampleID; row values = absolute counts of features

# read in dada2 feature table with absolute counts
A_fts_abs <- read.table(A_ifp_fts_tab, header = TRUE, sep = "\t", as.is = TRUE,
                        check.names = FALSE)

# read in representative sequences as df
A_rep_seq <- fasta_to_df(fasta_file = A_ifp_rep_seq, hdr_colname = "FeatureID",
                         seq_colname = "RepSeq")

# convert absolute feature counts to relative abundances
# indexing with the numeric output of the which() function preserves column...
# ... FeatureID no matter where it is located in the df
A_fts_rel <- A_fts_abs
A_fts_num <- which(names(A_fts_rel) == "FeatureID")
A_fts_rel[, -A_fts_num] <- lapply(A_fts_rel[, -A_fts_num],
                                  function(x) {x/sum(x)})

# optional, check conversion by running: colSums(A_fts_rel[, -A_fts_num])

# merge taxa, rep seqs, and fts dfs
A_seq_tax <- merge(A_rep_seq, A_tax_mrg, by = "FeatureID", sort = FALSE)
A_ebtks_abs <- merge(A_seq_tax, A_fts_abs, by = "FeatureID", sort = FALSE)
A_ebtks_rel <- merge(A_seq_tax, A_fts_rel, by = "FeatureID", sort = FALSE)

# ensure that column FeatureID is a character and not a factor
A_ebtks_abs$FeatureID <- as.character(A_ebtks_abs$FeatureID)
A_ebtks_rel$FeatureID <- as.character(A_ebtks_rel$FeatureID)

# A_ebtks_abs == Data File S1

### ************************************
### A - STEP 3 - remove contaminant features ----
### ************************************

## this step uses taxonomic assignments in A_tax_mrg ...
# ... to identify and remove features likely to be contamination
# e.g. chloroplast, mitohondria, eukaryotic, and Unassigned Kingdom

A_clp_int_ggs <- dplyr::filter(A_tax_mrg, grepl(pattern = "chloropl",
                                                x = A_tax_mrg$int.ggs.tax,
                                                ignore.case = TRUE))
A_clp_int_slv <- dplyr::filter(A_tax_mrg, grepl(pattern = "chloropl",
                                                x = A_tax_mrg$int.slv.tax,
                                                ignore.case = TRUE))
A_mtc_int_ggs <- dplyr::filter(A_tax_mrg, grepl(pattern = "mitochon",
                                                x = A_tax_mrg$int.ggs.tax,
                                                ignore.case = TRUE))
A_mtc_int_slv <- dplyr::filter(A_tax_mrg, grepl(pattern = "mitochon",
                                                x = A_tax_mrg$int.slv.tax,
                                                ignore.case = TRUE))
A_euk_int_ggs <- dplyr::filter(A_tax_mrg, grepl(pattern = "eukaryot",
                                                x = A_tax_mrg$int.ggs.tax,
                                                ignore.case = TRUE))
A_euk_int_slv <- dplyr::filter(A_tax_mrg, grepl(pattern = "eukaryot",
                                                x = A_tax_mrg$int.slv.tax,
                                                ignore.case = TRUE))

# careful here - cannot be partial match
A_unk_int_ggs <- dplyr::filter(A_tax_mrg, grepl(pattern = "Unassigned",
                                                x = A_tax_mrg$int.ggs.tax,
                                                fixed = TRUE))
A_unk_int_slv <- dplyr::filter(A_tax_mrg, grepl(pattern = "Unassigned",
                                                x = A_tax_mrg$int.slv.tax,
                                                fixed = TRUE))

# succesively merge Greengenes and SILVA contaminant dfs
# unequal df lengths are handled by not specifying a 'by' variable
A_clp_int_mrg <- merge(A_clp_int_ggs, A_clp_int_slv, all = TRUE, sort = FALSE)
A_mtc_int_mrg <- merge(A_mtc_int_ggs, A_mtc_int_slv, all = TRUE, sort = FALSE)
A_euk_int_mrg <- merge(A_euk_int_ggs, A_euk_int_slv, all = TRUE, sort = FALSE)
A_unk_int_mrg <- merge(A_unk_int_ggs, A_unk_int_slv, all = TRUE, sort = FALSE)

# combine contaminant dfs
A_ctm <- rbind(A_clp_int_mrg, A_mtc_int_mrg, A_euk_int_mrg, A_unk_int_mrg)

# remove contaminants from ebtks_abs df
A_ebtks_abs_ctm_rmv <- dplyr::anti_join(x = A_ebtks_abs, y = A_ctm,
                                        by = "FeatureID")

# create ctm dfs with absolute or relative abundances
A_ctm_abs <- dplyr::semi_join(x = A_ebtks_abs, y = A_ctm, by = "FeatureID")
A_ctm_rel <- dplyr::semi_join(x = A_ebtks_rel, y = A_ctm, by = "FeatureID")

### ************************************
### A - STEP 4 - remove samples ----
### ************************************

# this step removes samples with:
# low total feature count for features not considered contaminants
# high total relative abundance for features considered contaminants

# first, identify samples with low total feature count...
# remove unneeded columns
A_ebtks_abs_ctm_rmv_frm <- dplyr::select(A_ebtks_abs_ctm_rmv,
                                         -FeatureID, -RepSeq,
                                         -dplyr::contains("ggs"),
                                         -dplyr::contains("slv"))

# sum absolute feature counts for each sample
A_ebtks_abs_ctm_rmv_fts_sum <- data.frame(
  "SampleID" = names(A_ebtks_abs_ctm_rmv_frm),
  "Total" =  colSums(A_ebtks_abs_ctm_rmv_frm),
  row.names = NULL)

# create a vector of SampleIDs with low total feature count (to be removed)
A_low_fts_smp <- dplyr::filter(A_ebtks_abs_ctm_rmv_fts_sum, Total < 400)
A_low_fts_sID <- as.character(A_low_fts_smp$SampleID)

# second, identify samples with high contaminant relative abundance
# remove unneeded columns
A_ctm_rel_frm <- dplyr::select(A_ctm_rel,
                               -FeatureID, -RepSeq,
                               -dplyr::contains("ggs"),
                               -dplyr::contains("slv"))

# sum relative abundances for each sample
A_ctm_rel_frm_sum <- data.frame("SampleID" = names(A_ctm_rel_frm),
                                "Total" = colSums(A_ctm_rel_frm),
                                row.names = NULL)

# create a vector of SampleIDs with high contam rel abundance (to be removed)
A_hgh_ctm_smp <- dplyr::filter(A_ctm_rel_frm_sum, Total > 0.05)
A_hgh_ctm_sID <- as.character(A_hgh_ctm_smp$SampleID)

# remove the samples identified in the vectors above
A_ebtks_abs_smp_rmv <- dplyr::select(A_ebtks_abs_ctm_rmv,
                                     -dplyr::one_of(A_low_fts_sID),
                                     -dplyr::one_of(A_hgh_ctm_sID))

### ************************************
### A - STEP 5 - split by cohort ----
### ************************************

# this step splits A_ebtks_abs_smp_rmv in order to isolate:
# negative controls
# healthy mouse cohort; with and without negative controls
# NOTE: the aom/dss cohort is not considered; the focus here healthy mice

A_neg_abs <- dplyr::select(A_ebtks_abs_smp_rmv, FeatureID, RepSeq,
                           dplyr::contains("ggs"), 
                           dplyr::contains("slv"),
                           dplyr::contains("neg"))

A_hth_abs <- dplyr::select(A_ebtks_abs_smp_rmv, FeatureID, RepSeq,
                           dplyr::contains("ggs"), 
                           dplyr::contains("slv"),
                           dplyr::contains("HTH"))

A_hth_neg_abs <- dplyr::select(A_ebtks_abs_smp_rmv, FeatureID, RepSeq,
                               dplyr::contains("ggs"),
                               dplyr::contains("slv"),
                               dplyr::contains("HTH"), 
                               dplyr::contains("neg"))

# check that the correct number of samples were retained
A_neg_abs_chk <- dplyr::select(A_neg_abs, -FeatureID, -RepSeq,
                               -dplyr::contains("ggs"),
                               -dplyr::contains("slv"))

A_hth_abs_chk <- dplyr::select(A_hth_abs, -FeatureID, -RepSeq,
                               -dplyr::contains("ggs"),
                               -dplyr::contains("slv"))

A_hth_neg_abs_chk <- dplyr::select(A_hth_neg_abs, -FeatureID, -RepSeq,
                                   -dplyr::contains("ggs"),
                                   -dplyr::contains("slv"))

A_neg_abs_CHECK <- ncol(A_neg_abs_chk)
A_hth_abs_CHECK <- ncol(A_hth_abs_chk)
A_hth_neg_abs_CHECK <- ncol(A_hth_neg_abs_chk)

### ************************************
### A - STEP 6 - split by sample type ----
### ************************************

# this step splits relevant dfs from A - STEP 5 to isolate:
# nEDPf = neg control, caecum, distal and proximal colon, and faecal samples
# EDPf = caecum, distal and proximal colon, and faecal samples
# EDP = caecum and distal and proximal colon samples

A_hth_abs_nEDPf_0 <- dplyr::select(A_hth_neg_abs, FeatureID, RepSeq,
                                   dplyr::contains("ggs"),
                                   dplyr::contains("slv"),
                                   dplyr::contains("CEM"),
                                   dplyr::contains("CO"),
                                   dplyr::contains("FEC"),
                                   dplyr::contains("neg"))

A_hth_abs_EDPf_0 <- dplyr::select(A_hth_abs, FeatureID, RepSeq,
                                  dplyr::contains("ggs"),
                                  dplyr::contains("slv"),
                                  dplyr::contains("CEM"),
                                  dplyr::contains("CO"),
                                  dplyr::contains("FEC"))

A_hth_abs_EDP_0 <- dplyr::select(A_hth_abs, FeatureID, RepSeq,
                                 dplyr::contains("ggs"),
                                 dplyr::contains("slv"),
                                 dplyr::contains("CEM"),
                                 dplyr::contains("CO"))

# check that the correct number of samples were retained
A_hth_abs_nEDPf_0_chk <- dplyr::select(A_hth_abs_nEDPf_0, 
                                       -FeatureID, -RepSeq,
                                       -dplyr::contains("ggs"),
                                       -dplyr::contains("slv"))

A_hth_abs_EDPf_0_chk <- dplyr::select(A_hth_abs_EDPf_0, 
                                      -FeatureID, -RepSeq,
                                      -dplyr::contains("ggs"),
                                      -dplyr::contains("slv"))

A_hth_abs_EDP_0_chk <- dplyr::select(A_hth_abs_EDP_0, 
                                     -FeatureID, -RepSeq,
                                     -dplyr::contains("ggs"),
                                     -dplyr::contains("slv"))

A_hth_abs_nEDPf_0_CHECK <- ncol(A_hth_abs_nEDPf_0_chk)
A_hth_abs_EDPf_0_CHECK <- ncol(A_hth_abs_EDPf_0_chk)
A_hth_abs_EDP_0_CHECK <- ncol(A_hth_abs_EDP_0_chk)

### ************************************
### A - STEP 7 - remove features ----
### ************************************

# following the removal of samples...
# ... features with total count less than 2 are removed
# this reflects the threshold of 2 seen in the original dada2 table

# create new versions of dfs so the originals are not overwritten
A_hth_abs_nEDPf_1 <- A_hth_abs_nEDPf_0
A_hth_abs_EDPf_1 <- A_hth_abs_EDPf_0
A_hth_abs_EDP_1 <- A_hth_abs_EDP_0

# convert column FeatureID into row names and then remove unneeded columns
row.names(A_hth_abs_nEDPf_1) <- A_hth_abs_nEDPf_1$FeatureID
row.names(A_hth_abs_EDPf_1) <- A_hth_abs_EDPf_1$FeatureID
row.names(A_hth_abs_EDP_1) <- A_hth_abs_EDP_1$FeatureID

A_hth_abs_nEDPf_2 <- dplyr::select(A_hth_abs_nEDPf_1,
                                   -FeatureID, -RepSeq,
                                   -dplyr::contains("ggs"),
                                   -dplyr::contains("slv"))

A_hth_abs_EDPf_2 <- dplyr::select(A_hth_abs_EDPf_1,
                                  -FeatureID, -RepSeq,
                                  -dplyr::contains("ggs"),
                                  -dplyr::contains("slv"))

A_hth_abs_EDP_2 <- dplyr::select(A_hth_abs_EDP_1,
                                 -FeatureID, -RepSeq,
                                 -dplyr::contains("ggs"),
                                 -dplyr::contains("slv"))

# sum absolute counts for each feature
A_hth_abs_nEDPf_sum <- data.frame("FeatureID" = row.names(A_hth_abs_nEDPf_2),
                                  "FeatureTotal" = rowSums(A_hth_abs_nEDPf_2),
                                  row.names = NULL)

A_hth_abs_EDPf_sum <- data.frame("FeatureID" = row.names(A_hth_abs_EDPf_2),
                                 "FeatureTotal" = rowSums(A_hth_abs_EDPf_2),
                                 row.names = NULL)

A_hth_abs_EDP_sum <- data.frame("FeatureID" = row.names(A_hth_abs_EDP_2),
                                "FeatureTotal" = rowSums(A_hth_abs_EDP_2),
                                row.names = NULL)

# retain features from the _sum dfs with total absolute count >= 2
A_hth_abs_nEDPf_rtn <- dplyr::filter(A_hth_abs_nEDPf_sum, FeatureTotal >= 2)
A_hth_abs_nEDPf_rtn$FeatureID <- as.character(A_hth_abs_nEDPf_rtn$FeatureID)

A_hth_abs_EDPf_rtn <- dplyr::filter(A_hth_abs_EDPf_sum, FeatureTotal >= 2)
A_hth_abs_EDPf_rtn$FeatureID <- as.character(A_hth_abs_EDPf_rtn$FeatureID)

A_hth_abs_EDP_rtn <- dplyr::filter(A_hth_abs_EDP_sum, FeatureTotal >= 2)
A_hth_abs_EDP_rtn$FeatureID <- as.character(A_hth_abs_EDP_rtn$FeatureID)

# join the original subset dfs with the retain dfs
A_hth_abs_nEDPf <- dplyr::semi_join(x = A_hth_abs_nEDPf_0,
                                    y = A_hth_abs_nEDPf_rtn,
                                    by = "FeatureID")

A_hth_abs_EDPf <- dplyr::semi_join(x = A_hth_abs_EDPf_0,
                                   y = A_hth_abs_EDPf_rtn,
                                   by = "FeatureID")
A_hth_abs_EDP <- dplyr::semi_join(x = A_hth_abs_EDP_0,
                                  y = A_hth_abs_EDP_rtn,
                                  by = "FeatureID")

# check the number of features retained
A_hth_abs_nEDPf_CHECK <- nrow(A_hth_abs_nEDPf)
A_hth_abs_EDPf_CHECK <- nrow(A_hth_abs_EDPf)
A_hth_abs_EDP_CHECK <- nrow(A_hth_abs_EDP)

# A_hth_abs_nEDPf == Data File S2
# A_hth_abs_EDPf == Data File S3
# A_hth_abs_EDP == Data File S4

### ************************************
### A - WRITE OUTPUTS ----
### ************************************

# print any LOGICALs or CHECKs
print(A_int_ggs_LOGICAL) # TRUE
print(A_int_slv_LOGICAL) # TRUE
print(A_neg_abs_CHECK) # 16
print(A_hth_abs_CHECK) # 49
print(A_hth_neg_abs_CHECK) # 65
print(A_hth_abs_nEDPf_0_CHECK) # 57
print(A_hth_abs_EDPf_0_CHECK) # 41
print(A_hth_abs_EDP_0_CHECK) # 26
print(A_hth_abs_nEDPf_CHECK) # 1339
print(A_hth_abs_EDPf_CHECK) # 1118
print(A_hth_abs_EDP_CHECK) # 790

# relative paths from wd for outputs headed for storage in "supplement/"
write.table(sep = "\t", row.names = FALSE, x = A_ebtks_abs,
            file = ofs_ebtks_abs)
write.table(sep = "\t", row.names = FALSE, x = A_hth_abs_nEDPf,
            file = ofs_hth_abs_nEDPf)
write.table(sep = "\t", row.names = FALSE, x = A_hth_abs_EDPf,
            file = ofs_hth_abs_EDPf)
write.table(sep = "\t", row.names = FALSE, x = A_hth_abs_EDP,
            file = ofs_hth_abs_EDP)

# relative paths for outputs headed for storage in the "vault/"
write.table(sep = "\t", row.names = FALSE, x = A_ebtks_abs,
            file = A_ofv_ebtks_abs)
write.table(sep = "\t", row.names = FALSE, x = A_hth_abs_nEDPf,
            file = A_ofv_hth_abs_nEDPf)
write.table(sep = "\t", row.names = FALSE, x = A_hth_abs_EDPf,
            file = A_ofv_hth_abs_EDPf)
write.table(sep = "\t", row.names = FALSE, x = A_hth_abs_EDP,
            file = A_ofv_hth_abs_EDP)

# create a list of objects from section and save them as a workspace file
A_obj <- c(ls(pattern = "A_"), "fasta_to_df", "trunc_tax")
save(list = A_obj, file = A_ofv_WS)

### ************************************
### SECTION B - create Bifidobacterium plot ----
### ************************************

# purpose:
# visualize sample relative abundance for Bifidobacterium spp.

# main/supplementary items generated:
# Figure S1

### ************************************
### B - STEP 1 - format input dfs ----
### ************************************

# this step converts absolute to relative abundances...
# ... and filters to isolate taxa assigned as Bifidobacterium

# create a new version of the input df from Section A
B_hth_abs_nEDPf_0 <- A_hth_abs_nEDPf

# store taxonomic information (added  ack in later)
B_tax <- dplyr::select(B_hth_abs_nEDPf_0,
                       FeatureID,
                       dplyr::contains("ggs"),
                       dplyr::contains("slv"))

# convert column FeatureID into row names and then remove unneeded columns
B_hth_abs_nEDPf_1 <- B_hth_abs_nEDPf_0
row.names(B_hth_abs_nEDPf_1) <- B_hth_abs_nEDPf_1$FeatureID

B_hth_abs_nEDPf_2 <- dplyr::select(B_hth_abs_nEDPf_1,
                                   -FeatureID, -RepSeq,
                                   -dplyr::contains("ggs"),
                                   -dplyr::contains("slv"))

# convert absolute feature counts to relative abundances (% out of 100)
B_hth_rel_nEDPf_0 <- as.data.frame(apply(B_hth_abs_nEDPf_2, 2,
                                         function(x) {(x/sum(x) * 100)}))

# optional, check conversion by running: colSums(B_hth_rel_nEDPf_0) # == 100

# convert FeatureID back into a col and add the taxonomic info back
B_hth_rel_nEDPf_1 <- B_hth_rel_nEDPf_0
B_hth_rel_nEDPf_1$FeatureID <- row.names(B_hth_rel_nEDPf_1)
B_hth_rel_nEDPf_2 <- B_hth_rel_nEDPf_1
row.names(B_hth_rel_nEDPf_2) <- 1:nrow(B_hth_rel_nEDPf_2)

B_hth_rel_nEDPf_tax <- merge(x = B_tax, y = B_hth_rel_nEDPf_2,
                             by = "FeatureID", sort = FALSE)

# filter to isolate all Bifidobacterium assignments
B_hth_rel_bif_ggs <- dplyr::filter(B_hth_rel_nEDPf_tax,
                                   grepl(pattern = "bifido",
                                         x = B_hth_rel_nEDPf_tax$int.ggs.Genus,
                                         ignore.case = TRUE))

B_hth_rel_bif_slv <- dplyr::filter(B_hth_rel_nEDPf_tax,
                                   grepl(pattern = "bifido",
                                         x = B_hth_rel_nEDPf_tax$int.slv.Genus,
                                         ignore.case = TRUE))

# check the number of features isolated
B_hth_rel_bif_ggs_CHECK <- nrow(B_hth_rel_bif_ggs) # 9
B_hth_rel_bif_slv_CHECK <- nrow(B_hth_rel_bif_slv) # 12

# merge Greengenes and SILVA Bifidobacterium dfs
# unequal df lengths are handled by not specifying a 'by' variable
B_hth_rel_bif_mrg <- merge(B_hth_rel_bif_ggs, B_hth_rel_bif_slv,
                           all = TRUE, sort = FALSE)

# check the number of features in merged df
B_hth_rel_bif_mrg_CHECK <- nrow(B_hth_rel_bif_mrg) # 12

# two features were assigned as Bifidobacterium by SILVA; however...
# ... Greengenes assigned them as Peptostreptococcaceae & Clostridiales...
# ... remove these features before plotting
B_hth_rel_bif_flt_0 <- dplyr::filter(B_hth_rel_bif_mrg,
                                     int.ggs.Genus == "Bifidobacterium",
                                     int.slv.Genus == "Bifidobacterium")

# check the number of features retained
B_hth_rel_bif_flt_0_CHECK <- nrow(B_hth_rel_bif_flt_0) # 9

# create a vector of unwanted columns and then remove them
B_rmv_col <- c("int.ggs.tax", "int.ggs.Kingdom", "int.ggs.Phylum",
               "int.ggs.Class", "int.ggs.Order", "int.ggs.Family",
               "int.ggs.Genus", "int.ggs.Species", "int.ggs.cnf",
               "int.ggs.lws.lvl",
               "int.slv.tax", "int.slv.Kingdom", "int.slv.Phylum",
               "int.slv.Class", "int.slv.Order", "int.slv.Family",
               "int.slv.Genus", "int.slv.Species", "int.slv.cnf",
               "int.slv.lws.lvl")

B_hth_rel_bif_flt_1 <- dplyr::select(B_hth_rel_bif_flt_0,
                                     -dplyr::one_of(B_rmv_col))

# sum rels of identical taxa using Greengenes, (i.e. collapse by taxon)
# this step chains together several dplyr functions as follows:
# select to remove unwanted columns
# group_by to group data by unique variables (row values) in the Greengenes col
# summarise_all to collapse the table by adding all row values for groups
B_hth_rel_bif_ggs_sum <- dplyr::summarise_all(
  dplyr::group_by(
    dplyr::select(
      B_hth_rel_bif_flt_1, -FeatureID, -int.slv.lws.txn),
    int.ggs.lws.txn),
  .fun = sum)

# check the number of taxa after collapse
B_hth_rel_bif_ggs_sum_CHECK <- nrow(B_hth_rel_bif_ggs_sum) # 3

# reshape the df to convert SampleID (currently cols) into rows...
# ... with correseponding rels for each taxon
B_hth_rel_bif_ggs_sum_rsh <- reshape2::melt(B_hth_rel_bif_ggs_sum,
                                            id.vars = "int.ggs.lws.txn",
                                            variable.name = "SampleID",
                                            value.name = "rel")

# read in metadata file and merge with reshaped df
B_met <- read.table(ifp_met, header = TRUE, sep = "\t", as.is = TRUE)

B_hth_rel_bif_met_0 <- merge(x = B_hth_rel_bif_ggs_sum_rsh, y = B_met,
                             all = FALSE, sort = FALSE, by = "SampleID")

# create a vector of wanted columns and then retain them
B_kep_col <- c("DietID", "int.ggs.lws.txn", "rel", "Diet",
               "DietOneLetter", "Type", "TypeOneLetter")

B_hth_rel_bif_met_1 <- dplyr::select(B_hth_rel_bif_met_0,
                                     dplyr::one_of(B_kep_col))

# split by diet group
B_hth_rel_bif_n_0 <- dplyr::filter(B_hth_rel_bif_met_1,
                                   Diet == "Negative.control")

B_hth_rel_bif_C_0 <- dplyr::filter(B_hth_rel_bif_met_1,
                                   Diet == "Control")

B_hth_rel_bif_R_0 <- dplyr::filter(B_hth_rel_bif_met_1,
                                   Diet == "Rice.bran")

B_hth_rel_bif_F_0 <- dplyr::filter(B_hth_rel_bif_met_1,
                                   Diet == "Fermented.rice.bran")

# check that the correct number of samples were retained
B_hth_rel_bif_n_0_CHECK <- length(unique(B_hth_rel_bif_n_0$DietID)) # 16 samples
B_hth_rel_bif_C_0_CHECK <- length(unique(B_hth_rel_bif_C_0$DietID)) # 15 samples
B_hth_rel_bif_R_0_CHECK <- length(unique(B_hth_rel_bif_R_0$DietID)) # 13 samples
B_hth_rel_bif_F_0_CHECK <- length(unique(B_hth_rel_bif_F_0$DietID)) # 13 samples

### ************************************
### B - STEP 2 - format for ggplots ----
### ************************************

# add in scale for x-axis for negative control samples:
# create sorted df from unique values in DietID col, add new col, add Scale
B_n_unq <- data.frame("DietID" = sort(unique(B_hth_rel_bif_n_0$DietID)))
B_n_unq$Scale <- 1:length(B_n_unq$DietID)
B_hth_rel_bif_n_1 <- merge(x = B_hth_rel_bif_n_0, y = B_n_unq,
                           by = "DietID")

# add in scale for x-axis for mouse samples:
# split by sample type, then proceed with same steps as used above
# NOTE: faeces will be handled differently
# NOTE: scale is orderd by sample type E - DP - f with 1 space in between each
B_C_E <- dplyr::filter(B_hth_rel_bif_C_0, TypeOneLetter == "E")
B_R_E <- dplyr::filter(B_hth_rel_bif_R_0, TypeOneLetter == "E")
B_F_E <- dplyr::filter(B_hth_rel_bif_F_0, TypeOneLetter == "E")

B_C_D <- dplyr::filter(B_hth_rel_bif_C_0, TypeOneLetter == "D")
B_R_D <- dplyr::filter(B_hth_rel_bif_R_0, TypeOneLetter == "D")
B_F_D <- dplyr::filter(B_hth_rel_bif_F_0, TypeOneLetter == "D")

B_C_P <- dplyr::filter(B_hth_rel_bif_C_0, TypeOneLetter == "P")
B_R_P <- dplyr::filter(B_hth_rel_bif_R_0, TypeOneLetter == "P")
B_F_P <- dplyr::filter(B_hth_rel_bif_F_0, TypeOneLetter == "P")

B_C_E_unq <- data.frame("DietID" = sort(unique(B_C_E$DietID)))
B_R_E_unq <- data.frame("DietID" = sort(unique(B_R_E$DietID)))
B_F_E_unq <- data.frame("DietID" = sort(unique(B_F_E$DietID)))

B_C_D_unq <- data.frame("DietID" = sort(unique(B_C_D$DietID)))
B_R_D_unq <- data.frame("DietID" = sort(unique(B_R_D$DietID)))
B_F_D_unq <- data.frame("DietID" = sort(unique(B_F_D$DietID)))

B_C_P_unq <- data.frame("DietID" = sort(unique(B_C_P$DietID)))
B_R_P_unq <- data.frame("DietID" = sort(unique(B_R_P$DietID)))
B_F_P_unq <- data.frame("DietID" = sort(unique(B_F_P$DietID)))

# successively number the scales for each diet, and then bind each diet together
B_C_E_unq$Scale <- c(1:5)
B_C_D_unq$Scale <- c(7:8)
B_C_P_unq$Scale <- c(9:11)
B_C_EDP_unq <- rbind(B_C_E_unq, B_C_D_unq, B_C_P_unq)

B_R_E_unq$Scale <- c(1:4)
B_R_D_unq$Scale <- c(6:7)
B_R_P_unq$Scale <- c(8:9)
B_R_EDP_unq <- rbind(B_R_E_unq, B_R_D_unq, B_R_P_unq)

B_F_E_unq$Scale <- c(1:4)
B_F_D_unq$Scale <- c(6:7)
B_F_P_unq$Scale <- c(8:9)
B_F_EDP_unq <- rbind(B_F_E_unq, B_F_D_unq, B_F_P_unq)

B_hth_rel_bif_C_1 <- merge(x = B_hth_rel_bif_C_0, y = B_C_EDP_unq,
                           all = TRUE, by = "DietID")

B_hth_rel_bif_R_1 <- merge(x = B_hth_rel_bif_R_0, y = B_R_EDP_unq,
                           all = TRUE, by = "DietID")

B_hth_rel_bif_F_1 <- merge(x = B_hth_rel_bif_F_0, y = B_F_EDP_unq,
                           all = TRUE, by = "DietID")

# now for faeces which begins at 13 for C diet, and 11 for R and F diets
B_hth_rel_bif_C_1$Scale[B_hth_rel_bif_C_1$TypeOneLetter == "f1"] <- 13
B_hth_rel_bif_C_1$Scale[B_hth_rel_bif_C_1$TypeOneLetter == "f2"] <- 14
B_hth_rel_bif_C_1$Scale[B_hth_rel_bif_C_1$TypeOneLetter == "f6"] <- 15
B_hth_rel_bif_C_1$Scale[B_hth_rel_bif_C_1$TypeOneLetter == "f10"] <- 16
B_hth_rel_bif_C_1$Scale[B_hth_rel_bif_C_1$TypeOneLetter == "f14"] <- 17

B_hth_rel_bif_R_1$Scale[B_hth_rel_bif_R_1$TypeOneLetter == "f1"] <- 11
B_hth_rel_bif_R_1$Scale[B_hth_rel_bif_R_1$TypeOneLetter == "f2"] <- 12
B_hth_rel_bif_R_1$Scale[B_hth_rel_bif_R_1$TypeOneLetter == "f6"] <- 13
B_hth_rel_bif_R_1$Scale[B_hth_rel_bif_R_1$TypeOneLetter == "f10"] <- 14
B_hth_rel_bif_R_1$Scale[B_hth_rel_bif_R_1$TypeOneLetter == "f14"] <- 15

B_hth_rel_bif_F_1$Scale[B_hth_rel_bif_F_1$TypeOneLetter == "f1"] <- 11
B_hth_rel_bif_F_1$Scale[B_hth_rel_bif_F_1$TypeOneLetter == "f2"] <- 12
B_hth_rel_bif_F_1$Scale[B_hth_rel_bif_F_1$TypeOneLetter == "f6"] <- 13
B_hth_rel_bif_F_1$Scale[B_hth_rel_bif_F_1$TypeOneLetter == "f10"] <- 14
B_hth_rel_bif_F_1$Scale[B_hth_rel_bif_F_1$TypeOneLetter == "f14"] <- 15

# sort the above dfs by Scale
B_hth_rel_bif_C_2 <- B_hth_rel_bif_C_1[order(B_hth_rel_bif_C_1$Scale),]
B_hth_rel_bif_R_2 <- B_hth_rel_bif_R_1[order(B_hth_rel_bif_R_1$Scale),]
B_hth_rel_bif_F_2 <- B_hth_rel_bif_F_1[order(B_hth_rel_bif_F_1$Scale),]

# add new col specifying shapes for sample type
B_hth_rel_bif_C_3 <- B_hth_rel_bif_C_2
B_hth_rel_bif_R_3 <- B_hth_rel_bif_R_2
B_hth_rel_bif_F_3 <- B_hth_rel_bif_F_2

B_hth_rel_bif_C_3$shape[B_hth_rel_bif_C_3$TypeOneLetter == "E"] <- 17
B_hth_rel_bif_R_3$shape[B_hth_rel_bif_R_3$TypeOneLetter == "E"] <- 17
B_hth_rel_bif_F_3$shape[B_hth_rel_bif_F_3$TypeOneLetter == "E"] <- 17

B_hth_rel_bif_C_3$shape[B_hth_rel_bif_C_3$TypeOneLetter == "D"] <- 19
B_hth_rel_bif_R_3$shape[B_hth_rel_bif_R_3$TypeOneLetter == "D"] <- 19
B_hth_rel_bif_F_3$shape[B_hth_rel_bif_F_3$TypeOneLetter == "D"] <- 19

B_hth_rel_bif_C_3$shape[B_hth_rel_bif_C_3$TypeOneLetter == "P"] <- 1
B_hth_rel_bif_R_3$shape[B_hth_rel_bif_R_3$TypeOneLetter == "P"] <- 1
B_hth_rel_bif_F_3$shape[B_hth_rel_bif_F_3$TypeOneLetter == "P"] <- 1

B_hth_rel_bif_C_3$shape[B_hth_rel_bif_C_3$Type == "feces"] <- 4
B_hth_rel_bif_R_3$shape[B_hth_rel_bif_R_3$Type == "feces"] <- 4
B_hth_rel_bif_F_3$shape[B_hth_rel_bif_F_3$Type == "feces"] <- 4

# define stack locations for the three taxa (top, middle, bottom):
# reorder the df with all rels most to least % create a vector of unique taxa
# NOTE: add rev() to place most abuntant on bottom
B_odr_rel <- B_hth_rel_bif_met_1[order(B_hth_rel_bif_met_1$rel,
                                       decreasing = TRUE),]
B_odr_tax <- rev(unique(B_odr_rel$int.ggs.lws.txn))

# using the vector above, define stack location for each taxon
B_hth_rel_bif_n_1$int.ggs.lws.txn <- factor(B_hth_rel_bif_n_1$int.ggs.lws.txn,
                                            levels = B_odr_tax)
B_hth_rel_bif_C_3$int.ggs.lws.txn <- factor(B_hth_rel_bif_C_3$int.ggs.lws.txn,
                                            levels = B_odr_tax)
B_hth_rel_bif_R_3$int.ggs.lws.txn <- factor(B_hth_rel_bif_R_3$int.ggs.lws.txn,
                                            levels = B_odr_tax)
B_hth_rel_bif_F_3$int.ggs.lws.txn <- factor(B_hth_rel_bif_F_3$int.ggs.lws.txn,
                                            levels = B_odr_tax)

# hex code for black and grey
B_blk <- "#000000"
B_gry <- "#bbbbbb"

# font family
B_fnt_fam <- "Courier"

# axis parameters
B_x_lim_n <- c(0, 17)
B_x_lim_C <- c(0, 18)
B_x_lim_RF <- c(0, 16)
B_y_lim <- c(-12, 112)
B_y_lim_CRF <- c(-21, 112)
B_y_brk <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
B_y_lab <- c("0", "", "20", "", "40", "", "60", "", "80", "", "100")

# various text based labels
B_top_txt_n <- data.frame(x = B_x_lim_n[2] / 2, y = B_y_lim[2] - 6,
                          label = "Negative control")
B_top_lne_n <- data.frame(x = B_x_lim_n[1], xend = B_x_lim_n[2],
                          y = 100, yend = 100)

B_top_txt_C <- data.frame(x = B_x_lim_C[1] + 9, y = B_y_lim_CRF[2] - 6,
                          label = "Control")
B_top_lne_C <- data.frame(x = B_x_lim_C[1], xend = B_x_lim_C[2],
                          y = 100, yend = 100)

B_top_txt_R <- data.frame(x = B_x_lim_RF[1] + 8, y = B_y_lim_CRF[2] - 6,
                          label = "Rice bran")

B_top_txt_F <- data.frame(x = B_x_lim_RF[1] + 8, y = B_y_lim_CRF[2] - 6)

B_top_lne_RF <- data.frame(x = B_x_lim_RF[1], xend = B_x_lim_RF[2],
                           y = 100, yend = 100)

B_btm_txt_C_f <- data.frame(x = c(8, 10, 13:17), y = B_y_lim_CRF[1] * 0.5,
                            label = c("distal", "proximal",
                                      "1", "2", "6", "10", "14"))

B_btm_txt_RF_f <- data.frame(x = c(6.5, 8.5, 11:15), y = B_y_lim_CRF[1] * 0.5,
                             label = c("distal", "proximal",
                                       "1", "2", "6", "10", "14"))

B_btm_typ_txt_C <- data.frame(x = c(3, 9, 15),
                              y = B_y_lim_CRF[1] * 0.75,
                              label = c("caecum", "colon", "faeces"))

B_btm_typ_txt_RF <- data.frame(x = c(2.5, 7.5, 13),
                               y = B_y_lim_CRF[1] * 0.75,
                               label = c("caecum", "colon", "faeces"))

# plot and save parameters
B_lne_size <- 0.5
B_ggp_top_txt_size <- 3.3
B_ggp_txt_size <- 2.0
B_ggp_shp_size <- 1.3
B_ggp_shp_strk <- 0.2
B_ggp_axs_lab_size <- 8
B_ggp_axs_txt_size <- 7
B_lab_fnt <- list(size = 9, family = "Courier", face = "bold")

B_hex_n <- c("#a6a6a6", "#595959", "#262626") # 3 shades of gray (#969696)
B_hex_C <- c("#acd5a9", "#509b4b", "#2c562a") # 3 shades of green
B_hex_R <- c("#9ec8e0", "#3783ae", "#1f4961") # 3 shades of blue
B_hex_F <- c("#feb580", "#e56001", "#7f3601") # 3 shades of orange

# create within plot legend
B_tax_txt_lgn <- data.frame(x = 1.1, y = c(76, 82, 88),
                            label = c(" = Bifidobacterium pseudolongum",
                                      " = Bifidobacterium bifidum",
                                      " = Bifidobacterium"))
B_tax_shp_lgn <- data.frame(x = 1, y = B_tax_txt_lgn$y)

### ************************************
### B - STEP 3 - create ggplots ----
### ************************************

# negative controls (B_tax_txt_lgn plotted twice to darken)
B_ggp_hth_rel_bif_n <- ggplot(data = B_hth_rel_bif_n_1,
                              aes(x = Scale, y = rel, fill = int.ggs.lws.txn)) +
  geom_text(size = 1.6,
            inherit.aes = FALSE, data = B_tax_txt_lgn,
            aes(x = x, y = y, label = label), hjust = 0,
            color = B_hex_n, family = B_fnt_fam, fontface = "italic") +
  geom_text(size = 1.6,
            inherit.aes = FALSE, data = B_tax_txt_lgn,
            aes(x = x, y = y, label = label), hjust = 0,
            color = B_hex_n, family = B_fnt_fam, fontface = "italic") +
  geom_point(size = 2.0, stroke = 0,
             inherit.aes = FALSE, data = B_tax_shp_lgn, shape = 15,
             aes(x = x, y = y), color = B_hex_n) +
  geom_segment(size = B_lne_size,
               inherit.aes = FALSE, data = B_top_lne_n,
               aes(x = x, xend = xend, y = y, yend = yend),
               color = "#969696", lineend = "square", linetype = "dotted") +
  geom_text(size = B_ggp_top_txt_size,
            inherit.aes = FALSE, data = B_top_txt_n,
            aes(x = x, y = y, label = label),
            color = "#969696", family = B_fnt_fam) +
  geom_text(size = B_ggp_txt_size,
            aes(x = Scale, y = B_y_lim[1] * 0.5, label = TypeOneLetter),
            color = B_blk, family = B_fnt_fam) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = B_hex_n) +
  theme(panel.border = element_rect(fill = NA, size = B_lne_size,
                                    color = B_blk),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = B_ggp_axs_txt_size, color = B_blk),
        axis.ticks.y = element_line(size = B_lne_size, color = B_blk),
        text = element_text(family = B_fnt_fam)) +
  scale_y_continuous(limits = B_y_lim, breaks = B_y_brk, labels = B_y_lab,
                     expand = c(0,0)) +
  scale_x_continuous(limits = B_x_lim_n, expand = c(0,0))

# control diet (B_tax_txt_lgn plotted twice to darken)
B_ggp_hth_rel_bif_C <- ggplot(data = B_hth_rel_bif_C_3,
                              aes(x = Scale, y = rel, fill = int.ggs.lws.txn)) +
  geom_text(size = 1.6,
            inherit.aes = FALSE, data = B_tax_txt_lgn,
            aes(x = x, y = y, label = label), hjust = 0,
            color = B_hex_C, family = B_fnt_fam, fontface = "italic") +
  geom_text(size = 1.6,
            inherit.aes = FALSE, data = B_tax_txt_lgn,
            aes(x = x, y = y, label = label), hjust = 0,
            color = B_hex_C, family = B_fnt_fam, fontface = "italic") +
  geom_point(size = 2.0, stroke = 0,
             inherit.aes = FALSE, data = B_tax_shp_lgn, shape = 15,
             aes(x = x, y = y), color = B_hex_C) +
  geom_segment(size = B_lne_size,
               inherit.aes = FALSE, data = B_top_lne_C,
               aes(x = x, xend = xend, y = y, yend = yend),
               color = "#33a02c", lineend = "square", linetype = "dotdash") +
  geom_text(size = B_ggp_top_txt_size,
            inherit.aes = FALSE, data = B_top_txt_C,
            aes(x = x, y = y, label = label),
            color = "#33a02c", family = B_fnt_fam) +
  geom_text(size = B_ggp_txt_size,
            inherit.aes = FALSE, data = B_btm_typ_txt_C,
            aes(x = x, y = y, label = label),
            color = B_blk, family = B_fnt_fam) +
  geom_text(size = B_ggp_txt_size,
            inherit.aes = FALSE, data = B_btm_txt_C_f,
            aes(x = x, y = y, label = label),
            color = B_blk, family = B_fnt_fam) +
  geom_point(size = B_ggp_shp_size, stroke = B_ggp_shp_strk,
             inherit.aes = FALSE, shape = B_hth_rel_bif_C_3$shape,
             x = B_hth_rel_bif_C_3$Scale, y = B_y_lim[1] * 0.25,
             color = B_blk) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = B_hex_C) +
  theme(panel.border = element_rect(fill = NA, size = B_lne_size,
                                    color = B_blk),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = B_ggp_axs_txt_size, color = B_blk),
        axis.ticks.y = element_line(size = B_lne_size, color = B_blk),
        text = element_text(family = B_fnt_fam)) +
  scale_y_continuous(limits = B_y_lim_CRF, breaks = B_y_brk, labels = B_y_lab,
                     expand = c(0,0)) +
  scale_x_continuous(limits = B_x_lim_C, expand = c(0,0))

# rice bran diet (B_tax_txt_lgn plotted twice to darken)
B_ggp_hth_rel_bif_R <- ggplot(data = B_hth_rel_bif_R_3,
                              aes(x = Scale, y = rel, fill = int.ggs.lws.txn)) +
  geom_text(size = 1.6,
            inherit.aes = FALSE, data = B_tax_txt_lgn,
            aes(x = x, y = y, label = label), hjust = 0,
            color = B_hex_R, family = B_fnt_fam, fontface = "italic") +
  geom_text(size = 1.6,
            inherit.aes = FALSE, data = B_tax_txt_lgn,
            aes(x = x, y = y, label = label), hjust = 0,
            color = B_hex_R, family = B_fnt_fam, fontface = "italic") +
  geom_point(size = 2.0, stroke = 0,
             inherit.aes = FALSE, data = B_tax_shp_lgn, shape = 15,
             aes(x = x, y = y), color = B_hex_R) +
  geom_segment(size = B_lne_size,
               inherit.aes = FALSE, data = B_top_lne_RF,
               aes(x = x, xend = xend, y = y, yend = yend),
               color = "#1f78b4", lineend = "square", linetype = "dashed") +
  geom_text(size = B_ggp_top_txt_size,
            inherit.aes = FALSE, data = B_top_txt_R,
            aes(x = x, y = y, label = label),
            color = "#1f78b4", family = B_fnt_fam) +
  geom_text(size = B_ggp_txt_size,
            inherit.aes = FALSE, data = B_btm_txt_RF_f,
            aes(x = x, y = y, label = label),
            color = B_blk, family = B_fnt_fam) +
  geom_text(size = B_ggp_txt_size,
            inherit.aes = FALSE, data = B_btm_typ_txt_RF,
            aes(x = x, y = y, label = label),
            color = B_blk, family = B_fnt_fam) +
  geom_point(size = B_ggp_shp_size, stroke = B_ggp_shp_strk,
             inherit.aes = FALSE, shape = B_hth_rel_bif_R_3$shape,
             x = B_hth_rel_bif_R_3$Scale, y = B_y_lim[1] * 0.25,
             color = B_blk) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = B_hex_R) +
  theme(panel.border = element_rect(fill = NA, size = B_lne_size,
                                    color = B_blk),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = B_ggp_axs_txt_size, color = B_blk),
        axis.ticks.y = element_line(size = B_lne_size, color = B_blk),
        text = element_text(family = B_fnt_fam)) +
  scale_y_continuous(limits = B_y_lim_CRF, breaks = B_y_brk, labels = B_y_lab,
                     expand = c(0,0)) +
  scale_x_continuous(limits = B_x_lim_RF, expand = c(0,0))

# fermented rice bran diet (B_tax_txt_lgn plotted twice to darken)
B_ggp_hth_rel_bif_F <- ggplot(data = B_hth_rel_bif_F_3,
                              aes(x = Scale, y = rel, fill = int.ggs.lws.txn)) +
  geom_text(size = 1.6,
            inherit.aes = FALSE, data = B_tax_txt_lgn,
            aes(x = x, y = y, label = label), hjust = 0,
            color = B_hex_F, family = B_fnt_fam, fontface = "italic") +
  geom_text(size = 1.6,
            inherit.aes = FALSE, data = B_tax_txt_lgn,
            aes(x = x, y = y, label = label), hjust = 0,
            color = B_hex_F, family = B_fnt_fam, fontface = "italic") +
  geom_point(size = 2.0, stroke = 0,
             inherit.aes = FALSE, data = B_tax_shp_lgn, shape = 15,
             aes(x = x, y = y), color = B_hex_F) +
  geom_segment(size = B_lne_size,
               inherit.aes = FALSE, data = B_top_lne_RF,
               aes(x = x, xend = xend, y = y, yend = yend),
               color = "#e66101", lineend = "square", linetype = "solid") +
  geom_text(size = B_ggp_top_txt_size,
            inherit.aes = FALSE, data = B_top_txt_F,
            aes(x = x, y = y),
            color = "#e66101", family = B_fnt_fam,
            label = expression(
              paste(italic("B. longum"), "-fermented rice bran"))) +
  geom_text(size = B_ggp_txt_size,
            inherit.aes = FALSE, data = B_btm_txt_RF_f,
            aes(x = x, y = y, label = label),
            color = B_blk, family = B_fnt_fam) +
  geom_text(size = B_ggp_txt_size,
            inherit.aes = FALSE, data = B_btm_typ_txt_RF,
            aes(x = x, y = y, label = label),
            color = B_blk, family = B_fnt_fam) +
  geom_point(size = B_ggp_shp_size, stroke = B_ggp_shp_strk,
             inherit.aes = FALSE, shape = B_hth_rel_bif_F_3$shape,
             x = B_hth_rel_bif_F_3$Scale, y = B_y_lim[1] * 0.25,
             color = B_blk) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = B_hex_F) +
  theme(panel.border = element_rect(fill = NA, size = B_lne_size,
                                    color = B_blk),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = B_ggp_axs_txt_size, color = B_blk),
        axis.ticks.y = element_line(size = B_lne_size, color = B_blk),
        text = element_text(family = B_fnt_fam)) +
  scale_y_continuous(limits = B_y_lim_CRF, breaks = B_y_brk,
                     labels = B_y_lab, expand = c(0,0)) +
  scale_x_continuous(limits = B_x_lim_RF, expand = c(0,0))

### ************************************
### B - STEP 4 - arrange ggplots into panels ----
### ************************************

# arrange plots
B_gga_plt_0 <- ggpubr::ggarrange(B_ggp_hth_rel_bif_C,
                                 B_ggp_hth_rel_bif_R,
                                 B_ggp_hth_rel_bif_F,
                                 B_ggp_hth_rel_bif_n, 
                                 labels = "AUTO", font.label = B_lab_fnt,
                                 ncol = 1, nrow = 4, heights = c(1, 1, 1, 1),
                                 align = "hv")

# add universal y axis
B_gga_ytt_lab <- ggpubr::text_grob("Relative abundance (%)", color = B_blk,
                                   rot = 90, family = B_fnt_fam,
                                   size = 8)

B_gga_bif_plot <- ggpubr::annotate_figure(B_gga_plt_0, left = B_gga_ytt_lab)

### ************************************
### B - WRITE OUTPUTS ----
### ************************************

# print any LOGICALs or CHECKs
print(B_hth_rel_bif_ggs_CHECK) # 9
print(B_hth_rel_bif_slv_CHECK) # 12
print(B_hth_rel_bif_mrg_CHECK) # 12
print(B_hth_rel_bif_flt_0_CHECK) # 9
print(B_hth_rel_bif_ggs_sum_CHECK) # 3
print(B_hth_rel_bif_n_0_CHECK) # 16
print(B_hth_rel_bif_C_0_CHECK) # 15
print(B_hth_rel_bif_R_0_CHECK) # 13
print(B_hth_rel_bif_F_0_CHECK) # 13

# relative paths from wd for outputs headed for storage in "supplement/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = 165, height = 225,
       filename = ofs_hth_bif_plot, plot = B_gga_bif_plot)

# relative paths for outputs headed for storage in the "vault/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = 165, height = 225,
       filename = B_ofv_hth_bif_plot, plot = B_gga_bif_plot)

# create a list of objects from section and save them as a workspace file
B_obj <- c(ls(pattern = "B_"), "A_hth_abs_nEDPf")
save(list = B_obj, file = B_ofv_WS)

### ************************************
### SECTION C - principal components analysis ----
### ************************************

# purpose:
# this section creates a dataaset with count zero multiplicative replacement...
# ... of rounded zeros, followed by transformation using the ...
# ... centered log-ratio (clr) transformation over absolute abundances...
# ... and then creates a 'form' biplot from principal components analysis
# NOTE: form biplot is chosen as rays (representing features) are turned off

# main/supplementary items generated:
# Figure 3
# Figure S2

### ************************************
### C - FUNCTION ----
### ************************************

# function takes an input generated by the prcomp() and isolates...
# the value for the amount of variation explained by the PC specified
# NOTE: function isolates 1 value at a time
iso_var_val <- function(data = prcomp, digits, which.pc) {
  # internal check to ensure correct input class
  if (!inherits(data, "prcomp")) {
    stop("data must be class 'prcomp'")
  }
  # isolate standard deviations for all pcs
  pcs_sdv <- data$sdev
  # square each sdev, sum the squares, and divide each sdev by the sum
  pcs_sdv_sqr <- pcs_sdv^2
  pcs_sdv_sqr_sum <- sum(pcs_sdv_sqr)
  vals_raw <- pcs_sdv_sqr/pcs_sdv_sqr_sum
  # convert to percentage and round to specified number of digits
  vals_rnd <- (round(vals_raw, digits = digits)) * 100
  # isolate the specified PC
  var_val <- vals_rnd[which.pc]
  return(var_val)
}

# example usage: pc1 <- iso_var_val(data = prcomp, digits = 3, which.pc = 1)
# or place inside the paste function to create a text label useful for plotting
# lab_pc1 <- paste("PC 1 (", iso_var_val(data = prcomp, digits = 3,
#                                       which.pc = 1), "%)", sep ="")

### ************************************
### C - STEP 1 - CZM/clr/PCA ----
### ************************************

# create a new version of the input df from Section A (do this twice)
C_hth_abs_EDPf_0 <- A_hth_abs_EDPf
C_hth_abs_nEDPf_0 <- A_hth_abs_nEDPf

C_hth_abs_EDPf_1 <- C_hth_abs_EDPf_0
C_hth_abs_nEDPf_1 <- C_hth_abs_nEDPf_0

# convert column FeatureID into row names and then remove unneeded columns
row.names(C_hth_abs_EDPf_1) <- C_hth_abs_EDPf_1$FeatureID
row.names(C_hth_abs_nEDPf_1) <- C_hth_abs_nEDPf_1$FeatureID

C_hth_abs_EDPf_2 <- dplyr::select(C_hth_abs_EDPf_1, -FeatureID, -RepSeq,
                                  -dplyr::contains("ggs"),
                                  -dplyr::contains("slv"))
C_hth_abs_nEDPf_2 <- dplyr::select(C_hth_abs_nEDPf_1, -FeatureID, -RepSeq,
                                   -dplyr::contains("ggs"),
                                   -dplyr::contains("slv"))

# zero replacement
# input dfs (from above) are formatted with samples as cols and features as rows
# t() the input to cmultRepl() to convert features to cols and samples as rows
# t() for the output to maintain samples as cols and features as rows
C_hth_abs_EDPf_czm <- as.data.frame(
  t(zCompositions::cmultRepl(t(C_hth_abs_EDPf_2),
                             method = "CZM", output = "counts")))

C_hth_abs_nEDPf_czm <- as.data.frame(
  t(zCompositions::cmultRepl(t(C_hth_abs_nEDPf_2),
                             method = "CZM", output = "counts")))

# No. corrected values: 26 (C_hth_abs_EDPf_czm)
# No. corrected values: 5894 (C_hth_abs_nEDPf_czm)

# the input dfs are formatted with Samples as cols and Features as rows
C_hth_EDPf_clr_0 <- as.data.frame(apply(C_hth_abs_EDPf_czm, 2,
                                        function(x) {log2(x) - mean(log2(x))}))
C_hth_nEDPf_clr_0 <- as.data.frame(apply(C_hth_abs_nEDPf_czm, 2,
                                         function(x) {log2(x) - mean(log2(x))}))


# transpose dfs to convert features to cols and samples to rows and compute pcs
C_hth_EDPf_clr_1 <- as.data.frame(t(C_hth_EDPf_clr_0))
C_hth_nEDPf_clr_1 <- as.data.frame(t(C_hth_nEDPf_clr_0))

C_hth_EDPf_clr_pca <- prcomp(C_hth_EDPf_clr_1)
C_hth_nEDPf_clr_pca <- prcomp(C_hth_nEDPf_clr_1)

### ************************************
### C - STEP 2 - format for ggplots ----
### ************************************

# create axis labels for the amount of variation explained by pca1 and pca2
C_hth_EDPf_lab_pc1 <- paste("PC 1 (", iso_var_val(data = C_hth_EDPf_clr_pca,
                                                  digits = 3, which.pc = 1),
                            "%)", sep = "")
C_hth_EDPf_lab_pc2 <- paste("PC 2 (", iso_var_val(data = C_hth_EDPf_clr_pca,
                                                  digits = 3, which.pc = 2),
                            "%)", sep = "")
C_hth_nEDPf_lab_pc1 <- paste("PC 1 (", iso_var_val(data = C_hth_nEDPf_clr_pca,
                                                   digits = 3, which.pc = 1),
                             "%)", sep = "")
C_hth_nEDPf_lab_pc2 <- paste("PC 2 (", iso_var_val(data = C_hth_nEDPf_clr_pca,
                                                   digits = 3, which.pc = 2),
                             "%)", sep = "")


C_hth_EDPf_gbp <- ggbiplot::ggbiplot(C_hth_EDPf_clr_pca, 
                                     scale = 0, var.axes = FALSE,
                                     labels = row.names(C_hth_EDPf_clr_1))
C_hth_nEDPf_gbp <- ggbiplot::ggbiplot(C_hth_nEDPf_clr_pca,, 
                                      scale = 0, var.axes = FALSE,
                                      labels = row.names(C_hth_nEDPf_clr_1))

#print(C_hth_EDPf_gbp)
#print(C_hth_nEDPf_gbp)

# the above plots need a spot of work to create a clear and meaningful visual
# if we isolate the underyling plot data and combine it with our metadata...
# ... we have much greater control over customizing the aesthetics

# create a data.frame from the default ggbiplot objects that contains:
# xvar-yvar-label (i.e. the x and y coordinates for each sample)

C_hth_EDPf_crd <- C_hth_EDPf_gbp[["data"]]
C_hth_nEDPf_crd <- C_hth_nEDPf_gbp[["data"]]

# read in metadata file and merge with _crd dfs
C_met <- read.table(ifp_met, header = TRUE, sep = "\t", as.is = TRUE)

C_hth_EDPf_crd_met <- merge(x = C_hth_EDPf_crd, y = C_met, all = FALSE,
                            sort = FALSE, by.x = "labels", by.y = "SampleID")
C_hth_nEDPf_crd_met <- merge(x = C_hth_nEDPf_crd, y = C_met, all = FALSE,
                             sort = FALSE, by.x = "labels", by.y = "SampleID")

# isolate negative control samples from C_hth_nEDPf_crd_met
C_hth_n_crd_met <- dplyr::filter(C_hth_nEDPf_crd_met,
                                 Cohort == "negative.control")

# define color codes and shapes for samples
# plot order is: E D P f
C_hex_type <- c("#2d004b", "#8073ac", "#542788", "#8c510a")
C_shp_type <- c(17, 19, 1, 4)

# plot order is: E D P f n
C_hex_type_n <- c("#2d004b", "#8073ac", "#542788", "#8c510a", "#636363")
C_shp_type_n <- c(17, 19, 1, 4, 32)

# plot order is: C F R
C_hex_diet <- c("#33a02c", "#e66101", "#1f78b4")

# plot order is: C F n R
C_hex_diet_n <- c("#33a02c", "#e66101", "#636363", "#1f78b4")

# hex code for black and grey
C_blk <- "#000000"
C_gry <- "#bbbbbb"

# font family
C_fnt_fam <- "Courier"

# axis parameters
C_x_lim_EDPf <- c(-58, 58)
C_x_brk_EDPf <- c(-40, -20, 0, 20, 40)
C_x_lab_EDPf <- c("-40", "-20", "0", "20", "40")

C_y_lim_EDPf <- c(-55.5, 55.5)
C_y_brk_EDPf <- c(-40, -20, 0, 20, 40)
C_y_lab_EDPf <- c("-40", "-20", "0", "20", "40")

C_x_lim_nEDPf <- c(-83, 55)
C_x_brk_nEDPf <- c(-80, -60, -40, -20, 0, 20, 40)
C_x_lab_nEDPf <- c("-80", "-60", "-40", "-20", "0", "20", "40")

C_y_lim_nEDPf <- c(-51, 66)
C_y_brk_nEDPf <- c(-40, -20, 0, 20, 40)
C_y_lab_nEDPf <- c("-40", "-20", "0", "20", "40")

# create data.frame to plot dotted lines for x and y axes
C_dot_seg_EDPf <- data.frame(x = c(C_x_lim_EDPf[1], C_x_lim_EDPf[2], 0, 0),
                             y = c(0, 0, C_y_lim_EDPf[1], C_y_lim_EDPf[2]), 
                             xend = 0, yend = 0)

C_dot_seg_nEDPf <- data.frame(x = c(C_x_lim_nEDPf[1], C_x_lim_nEDPf[2], 0, 0),
                              y = c(0, 0, C_y_lim_nEDPf[1], C_y_lim_nEDPf[2]), 
                              xend = 0, yend = 0)

# plot and save parameters
C_ggp_txt_size <- 2.4
C_ggp_shp_size <- 1.6
C_ggp_shp_strk <- 0.4
C_ggp_pcs_lab_size <- 8
C_ggp_axs_txt_size <- 6
C_lab_fnt <- list(size = 13, family = "Courier", face = "bold")

### ************************************
### C - STEP 3 - create ggplots ----
### ************************************

# EDPf colored/shaped by sample type
C_ggp_hth_EDPf_type <- ggplot(data = C_hth_EDPf_crd_met) +
  geom_segment(inherit.aes = FALSE, data = C_dot_seg_EDPf, size = 0.3,
               aes(x = x, xend = xend, y = y, yend = yend),
               color = C_gry, linetype = "dotted") +
  geom_point(size = C_ggp_shp_size, stroke = C_ggp_shp_strk,
             (aes(x = xvar, y = yvar, shape = Type, color = Type))) +
  labs(x = C_hth_EDPf_lab_pc1, y = C_hth_EDPf_lab_pc2) +
  scale_color_manual(values = C_hex_type) +
  scale_shape_manual(values = C_shp_type) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, color = C_blk),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = C_ggp_pcs_lab_size, color = C_blk),
        axis.text = element_text(size = C_ggp_axs_txt_size, color = C_blk),
        axis.ticks = element_line(size = 0.25, color = C_blk),
        text = element_text(family = C_fnt_fam)) +
  scale_x_continuous(limits = C_x_lim_EDPf, breaks = C_x_brk_EDPf, 
                     labels = C_x_lab_EDPf, expand = c(0,0)) +
  scale_y_continuous(limits = C_y_lim_EDPf, breaks = C_y_brk_EDPf, 
                     labels = C_y_lab_EDPf, expand = c(0,0))

# EDPf colored by study diet (text drawn thrice to avoid using bold fontface)
C_ggp_hth_EDPf_diet <- ggplot(data = C_hth_EDPf_crd_met) +
  geom_segment(inherit.aes = FALSE, data = C_dot_seg_EDPf, size = 0.3,
               aes(x = x, xend = xend, y = y, yend = yend),
               color = C_gry, linetype = "dotted") +
  geom_text(size = C_ggp_txt_size, family = C_fnt_fam,
            (aes(x = xvar, y = yvar, label = DietOneLetter, color = Diet))) +
  geom_text(size = C_ggp_txt_size, family = C_fnt_fam,
            (aes(x = xvar, y = yvar, label = DietOneLetter, color = Diet))) +
  geom_text(size = C_ggp_txt_size, family = C_fnt_fam,
            (aes(x = xvar, y = yvar, label = DietOneLetter, color = Diet))) +
  labs(x = C_hth_EDPf_lab_pc1, y = C_hth_EDPf_lab_pc2) +
  scale_color_manual(values = C_hex_diet) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, color = C_blk),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = C_ggp_pcs_lab_size, color = C_blk),
        axis.text = element_text(size = C_ggp_axs_txt_size, color = C_blk),
        axis.ticks = element_line(size = 0.25, color = C_blk),
        text = element_text(family = C_fnt_fam)) +
  scale_x_continuous(limits = C_x_lim_EDPf, breaks = C_x_brk_EDPf, 
                     labels = C_x_lab_EDPf, expand = c(0,0)) +
  scale_y_continuous(limits = C_y_lim_EDPf, breaks = C_y_brk_EDPf, 
                     labels = C_y_lab_EDPf, expand = c(0,0))

# nEDPf colored/shaped by sample type (text drawn thrice)
C_ggp_hth_nEDPf_type <- ggplot(data = C_hth_nEDPf_crd_met) +
  geom_segment(inherit.aes = FALSE, data = C_dot_seg_nEDPf, size = 0.3,
               aes(x = x, xend = xend, y = y, yend = yend),
               color = C_gry, linetype = "dotted") +
  geom_point(size = C_ggp_shp_size, stroke = C_ggp_shp_strk,
             (aes(x = xvar, y = yvar, shape = Type, color = Type))) +
  geom_text(inherit.aes = FALSE, data = C_hth_n_crd_met, family = C_fnt_fam,
            size = C_ggp_txt_size, color = "#636363",
            (aes(x = xvar, y = yvar, label = TypeOneLetter))) +
  geom_text(inherit.aes = FALSE, data = C_hth_n_crd_met, family = C_fnt_fam,
            size = C_ggp_txt_size, color = "#636363",
            (aes(x = xvar, y = yvar, label = TypeOneLetter))) +
  geom_text(inherit.aes = FALSE, data = C_hth_n_crd_met, family = C_fnt_fam,
            size = C_ggp_txt_size, color = "#636363",
            (aes(x = xvar, y = yvar, label = TypeOneLetter))) +
  labs(x = C_hth_nEDPf_lab_pc1, y = C_hth_nEDPf_lab_pc2) +
  scale_color_manual(values = C_hex_type_n) +
  scale_shape_manual(values = C_shp_type_n) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, color = C_blk),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = C_ggp_pcs_lab_size, color = C_blk),
        axis.text = element_text(size = C_ggp_axs_txt_size, color = C_blk),
        axis.ticks = element_line(size = 0.25, color = C_blk),
        text = element_text(family = C_fnt_fam)) +
  scale_x_continuous(limits = C_x_lim_nEDPf, breaks = C_x_brk_nEDPf, 
                     labels = C_x_lab_nEDPf, expand = c(0,0)) +
  scale_y_continuous(limits = C_y_lim_nEDPf, breaks = C_y_brk_nEDPf, 
                     labels = C_y_lab_nEDPf, expand = c(0,0))

# nEDPf colored by study diet (text drawn thrice to avoid using bold fontface)
C_ggp_hth_nEDPf_diet <- ggplot(data = C_hth_nEDPf_crd_met) +
  geom_segment(inherit.aes = FALSE, data = C_dot_seg_nEDPf, size = 0.3,
               aes(x = x, xend = xend, y = y, yend = yend),
               color = C_gry, linetype = "dotted") +
  geom_text(size = C_ggp_txt_size, family = C_fnt_fam,
            (aes(x = xvar, y = yvar, label = DietOneLetter, color = Diet))) +
  geom_text(size = C_ggp_txt_size, family = C_fnt_fam,
            (aes(x = xvar, y = yvar, label = DietOneLetter, color = Diet))) +
  geom_text(size = C_ggp_txt_size, family = C_fnt_fam,
            (aes(x = xvar, y = yvar, label = DietOneLetter, color = Diet))) +
  labs(x = C_hth_nEDPf_lab_pc1, y = C_hth_nEDPf_lab_pc2) +
  scale_color_manual(values = C_hex_diet_n) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, color = C_blk),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = C_ggp_pcs_lab_size, color = C_blk),
        axis.text = element_text(size = C_ggp_axs_txt_size, color = C_blk),
        axis.ticks = element_line(size = 0.25, color = C_blk),
        text = element_text(family = C_fnt_fam)) +
  scale_x_continuous(limits = C_x_lim_nEDPf, breaks = C_x_brk_nEDPf, 
                     labels = C_x_lab_nEDPf, expand = c(0,0)) +
  scale_y_continuous(limits = C_y_lim_nEDPf, breaks = C_y_brk_nEDPf, 
                     labels = C_y_lab_nEDPf, expand = c(0,0))

### ************************************
### C - STEP 4 - create ggplot legends ----
### ************************************

# NOTE: legends with geom_text require alterations to ggplot2's default...
# ... legend paramaters, so store the original text key parameters

# store default (original) key
C_ori_key <- GeomText$draw_key

# create legend for sample type w/o negative controls
# create df
C_lgn_type <- data.frame(x = 1, y = 1, "label" = c("E", "D", "P", "f"))

# create ggplot
C_ggp_lgn_type <- ggplot(data = C_lgn_type, aes(x = x, y = y)) +
  geom_point(size = 1.2, stroke = 0.3, aes(shape = label, color = label)) +
  scale_color_manual(values = C_hex_type, name = NULL,
                     labels = c("caecum",
                                "colon (distal)",
                                "colon (proximal)",
                                "faeces")) +
  scale_shape_manual(values = C_shp_type, name = NULL,
                     labels = c("caecum",
                                "colon (distal)",
                                "colon (proximal)",
                                "faeces")) +
  theme_void(base_family = C_fnt_fam) +
  theme(legend.position = c(0.5, 0.5),
        legend.direction = "vertical",
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = NA, color = NA),
        legend.key.size = unit(0, "mm"),
        legend.text = element_text(size = 5.5, color = C_blk, hjust = 0,
                                   margin = margin(-0.1, 0, -0.1, -1.4, 
                                                   "mm"))) +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0))

# nest the above legend within respective plot
C_ggp_hth_EDPf_type_lgn_0 <- C_ggp_hth_EDPf_type + 
  annotation_custom(ggplotGrob(C_ggp_lgn_type), 
                    xmin = C_x_lim_EDPf[1] * 0.645, 
                    xmax = C_x_lim_EDPf[1] * 0.645, 
                    ymin = C_y_lim_EDPf[2] * 0.865, 
                    ymax = C_y_lim_EDPf[2] * 0.865)

# create legend for study diet w/o negative controls
# define vectors for a new key
C_txt_var_diet <- c("C", "F", "R")
C_txt_col_diet <- c("#33a02c", "#e66101", "#1f78b4")

# thank you to stackoverlow user20650 for the function below
# https://stackoverflow.com/questions/49965758/
# change-geom-texts-default-a-legend-to-label-string-itself
GeomText$draw_key <- function(data, params, size,
                              var = C_txt_var_diet,
                              cols = C_txt_col_diet) {
  # sort as ggplot sorts these alphanumerically / or levels of factor
  txt <- if(is.factor(var)) levels(var) else sort(var)
  txt <- txt[match(data$colour, cols)]
  textGrob(txt, 0.5, 0.5, just = "center", gp = gpar(
    col = alpha(data$colour, data$alpha),
    fontfamily = data$family,
    fontface = data$fontface,
    fontsize = data$size * .pt))
}

# create df
C_lgn_diet <- data.frame(x = 1, y = 1, "label" = c("C", "R", "F"))

# create ggplot (text drawn thrice to avoid using bold fontface)
C_ggp_lgn_diet <- ggplot(data = C_lgn_diet, aes(x = x, y = y)) +
  geom_text(size = 2.0, aes(label = label, color = label), family = C_fnt_fam) +
  geom_text(size = 2.0, aes(label = label, color = label), family = C_fnt_fam) +
  geom_text(size = 2.0, aes(label = label, color = label), family = C_fnt_fam) +
  scale_color_manual(values = c("#33a02c", "#1f78b4", "#e66101"), name = NULL,
                     labels = c("Control",
                                "Rice bran",
                                expression(paste("Fermented (", 
                                                 italic("B.longum"),")")))) +
  theme_void(base_family = C_fnt_fam) +
  theme(legend.position = c(0.5, 0.5),
        legend.direction = "vertical",
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = NA, color = NA),
        legend.key.size = unit(0, "mm"),
        legend.text = element_text(size = 5.5, color = C_blk, hjust = 0,
                                   margin = margin(-0.6, 0, -0.6, -1.4, 
                                                   "mm"))) +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0))

# nest the above legend within respective plot
C_ggp_hth_EDPf_diet_lgn_0 <- C_ggp_hth_EDPf_diet + 
  annotation_custom(ggplotGrob(C_ggp_lgn_diet), 
                    xmin = C_x_lim_EDPf[1] * 0.573, 
                    xmax = C_x_lim_EDPf[1] * 0.573, 
                    ymin = C_y_lim_EDPf[2] * 0.899, 
                    ymax = C_y_lim_EDPf[2] * 0.899)

# create legend for sample type w/ negative controls
# define vectors for a new key
C_txt_var_type_n <- c("", "", "", "", "n")
C_txt_col_type_n <- c("#2d004b", "#8073ac", "#542788", "#8c510a", "#636363")

# thank you to stackoverlow user20650 for the function below
# https://stackoverflow.com/questions/49965758/
# change-geom-texts-default-a-legend-to-label-string-itself
GeomText$draw_key <- function(data, params, size,
                              var = C_txt_var_type_n,
                              cols = C_txt_col_type_n) {
  # sort as ggplot sorts these alphanumerically / or levels of factor
  txt <- if(is.factor(var)) levels(var) else sort(var)
  txt <- txt[match(data$colour, cols)]
  textGrob(txt, 0.5, 0.5, just = "center", gp = gpar(
    col = alpha(data$colour, data$alpha),
    fontfamily = data$family,
    fontface = data$fontface,
    fontsize = data$size * .pt))
}

# create df
C_lgn_type_n <- data.frame(x = 1, y = 1, "label" = c("E", "D", "P", "f", "n"))

# create ggplot (text drawn thrice to avoid using bold fontface)
C_ggp_lgn_type_n <- ggplot(data = C_lgn_type_n, aes(x = x, y = y)) +
  geom_point(size = 1.2, stroke = 0.3, aes(shape = label, color = label)) +
  geom_text(size = 2.0, aes(label = label, color = label), family = C_fnt_fam) +
  geom_text(size = 2.0, aes(label = label, color = label), family = C_fnt_fam) +
  geom_text(size = 2.0, aes(label = label, color = label), family = C_fnt_fam) +
  scale_color_manual(values = C_hex_type_n, name = NULL,
                     labels = c("caecum",
                                "colon (distal)",
                                "colon (proximal)",
                                "faeces",
                                "negative control")) +
  scale_shape_manual(values = c(C_shp_type, 32), 
                     name = NULL,
                     labels = c("caecum",
                                "colon (distal)",
                                "colon (proximal)",
                                "faeces",
                                "negative control")) +
  theme_void(base_family = C_fnt_fam) +
  theme(legend.position = c(0.5, 0.5),
        legend.direction = "vertical",
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = NA, color = NA),
        legend.key.size = unit(0, "mm"),
        legend.text = element_text(size = 5.5, color = C_blk, hjust = 0,
                                   margin = margin(-0.1, 0, -0.1, -1.4, 
                                                   "mm"))) +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0))

# nest the above legend within respective plot
C_ggp_hth_nEDPf_type_lgn_0 <- C_ggp_hth_nEDPf_type + 
  annotation_custom(ggplotGrob(C_ggp_lgn_type_n), 
                    xmin = C_x_lim_nEDPf[1] * 0.702, 
                    xmax = C_x_lim_nEDPf[1] * 0.702, 
                    ymin = C_y_lim_nEDPf[2] * 0.84, 
                    ymax = C_y_lim_nEDPf[2] * 0.84)

# create legend for study diet w/ negative controls
# define vectors for a new key
C_txt_var_diet_n <- c("C", "F", "n", "R")
C_txt_col_diet_n <- c("#33a02c", "#e66101", "#636363", "#1f78b4")

GeomText$draw_key <- function(data, params, size,
                              var = C_txt_var_diet_n,
                              cols = C_txt_col_diet_n) {
  # sort as ggplot sorts these alphanumerically / or levels of factor
  txt <- if(is.factor(var)) levels(var) else sort(var)
  txt <- txt[match(data$colour, cols)]
  textGrob(txt, 0.5, 0.5, just = "center", gp = gpar(
    col = alpha(data$colour, data$alpha),
    fontfamily = data$family,
    fontface = data$fontface,
    fontsize = data$size * .pt))
}

# create df
C_lgn_diet_n <- data.frame(x = 1, y = 1, "label" = c("a", "b", "c", "d"))

# create ggplot (text drawn thrice to avoid using bold fontface)
C_ggp_lgn_diet_n <- ggplot(data = C_lgn_diet_n, aes(x = x, y = y)) +
  geom_text(size = 2.0, aes(label = label, color = label), family = C_fnt_fam) +
  geom_text(size = 2.0, aes(label = label, color = label), family = C_fnt_fam) +
  geom_text(size = 2.0, aes(label = label, color = label), family = C_fnt_fam) +
  scale_color_manual(values = c("#33a02c", "#1f78b4", "#e66101", "#636363"),
                     name = NULL,
                     labels = c("Control",
                                "Rice bran",
                                expression(paste("Fermented (", 
                                                 italic("B.longum"),")")),
                                "negative control")) +
  theme_void(base_family = C_fnt_fam) +
  theme(legend.position = c(0.5, 0.5),
        legend.direction = "vertical",
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = NA, color = NA),
        legend.key.size = unit(0, "mm"),
        legend.text = element_text(size = 5.5, color = C_blk, hjust = 0,
                                   margin = margin(-0.6, 0, -0.6, -1.4, 
                                                   "mm"))) +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0))

# nest the above legend within respective plot
C_ggp_hth_nEDPf_diet_lgn_0 <- C_ggp_hth_nEDPf_diet +
  annotation_custom(ggplotGrob(C_ggp_lgn_diet_n), 
                    xmin = C_x_lim_nEDPf[1] * 0.64, 
                    xmax = C_x_lim_nEDPf[1] * 0.64, 
                    ymin = C_y_lim_nEDPf[2] * 0.87, 
                    ymax = C_y_lim_nEDPf[2] * 0.87)

### ************************************
### C - STEP 5 - format/arrange ggplots into panels ----
### ************************************

# create dfs with relevant coordinates to make a box around nested legends
C_lgn_seg_type <- data.frame(x = c(C_x_lim_EDPf[1], 
                                   C_x_lim_EDPf[1] + 43.5), 
                             xend = c(C_x_lim_EDPf[1] + 43.5,
                                      C_x_lim_EDPf[1] + 43.5), 
                             y = c(C_y_lim_EDPf[2] - 19,
                                   C_y_lim_EDPf[2]),
                             yend = c(C_y_lim_EDPf[2] - 19,
                                      C_y_lim_EDPf[2] - 19))

C_lgn_seg_diet <- data.frame(x = c(C_x_lim_EDPf[1], 
                                   C_x_lim_EDPf[1] + 53.5), 
                             xend = c(C_x_lim_EDPf[1] + 53.5,
                                      C_x_lim_EDPf[1] + 53.5), 
                             y = c(C_y_lim_EDPf[2] - 15.5,
                                   C_y_lim_EDPf[2]),
                             yend = c(C_y_lim_EDPf[2] - 15.5,
                                      C_y_lim_EDPf[2] - 15.5))

C_lgn_seg_type_n <- data.frame(x = c(C_x_lim_nEDPf[1], 
                                     C_x_lim_nEDPf[1] + 54), 
                               xend = c(C_x_lim_nEDPf[1] + 54,
                                        C_x_lim_nEDPf[1] + 54), 
                               y = c(C_y_lim_nEDPf[2] - 26,
                                     C_y_lim_nEDPf[2]),
                               yend = c(C_y_lim_nEDPf[2] - 26,
                                        C_y_lim_nEDPf[2] - 26))

C_lgn_seg_diet_n <- data.frame(x = c(C_x_lim_nEDPf[1], 
                                     C_x_lim_nEDPf[1] + 64.5), 
                               xend = c(C_x_lim_nEDPf[1] + 64.5,
                                        C_x_lim_nEDPf[1] + 64.5), 
                               y = c(C_y_lim_nEDPf[2] - 21.5,
                                     C_y_lim_nEDPf[2]),
                               yend = c(C_y_lim_nEDPf[2] - 21.5,
                                        C_y_lim_nEDPf[2] - 21.5))

# add boxes to ggplots
C_ggp_hth_EDPf_type_lgn_1 <- C_ggp_hth_EDPf_type_lgn_0 +
  geom_segment(data = C_lgn_seg_type, size = 0.275, color = C_blk,
               aes(x = x, y = y, xend = xend, yend = yend),
               lineend = "square", linejoin = "mitre")

C_ggp_hth_EDPf_diet_lgn_1 <- C_ggp_hth_EDPf_diet_lgn_0 +
  geom_segment(data = C_lgn_seg_diet, size = 0.275, color = C_blk,
               aes(x = x, y = y, xend = xend, yend = yend),
               lineend = "square", linejoin = "mitre")

C_ggp_hth_nEDPf_type_lgn_1 <- C_ggp_hth_nEDPf_type_lgn_0 +
  geom_segment(data = C_lgn_seg_type_n, size = 0.275, color = C_blk,
               aes(x = x, y = y, xend = xend, yend = yend),
               lineend = "square", linejoin = "mitre")

C_ggp_hth_nEDPf_diet_lgn_1 <- C_ggp_hth_nEDPf_diet_lgn_0 +
  geom_segment(data = C_lgn_seg_diet_n, size = 0.275, color = C_blk,
               aes(x = x, y = y, xend = xend, yend = yend),
               lineend = "square", linejoin = "mitre")

# arrange plots
C_gga_EDPf_pca <- ggpubr::ggarrange(C_ggp_hth_EDPf_type_lgn_1,
                                    C_ggp_hth_EDPf_diet_lgn_1,
                                    labels = "AUTO", font.label = C_lab_fnt,
                                    ncol = 1, nrow = 2, widths = c(1,1),
                                    align = "hv")

C_gga_nEDPf_pca <- ggpubr::ggarrange(C_ggp_hth_nEDPf_type_lgn_1,
                                     C_ggp_hth_nEDPf_diet_lgn_1,
                                     labels = "AUTO", font.label = C_lab_fnt,
                                     ncol = 2, nrow = 1, widths = c(1,1),
                                     align = "hv")

### ************************************
### C - WRITE OUTPUTS ----
### ************************************

# reset the default ggplot key
GeomText$draw_key <- C_ori_key

# relative paths from wd for outputs headed for storage in "main/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = 75, height = 120,
       filename = ofm_hth_EDPf_pca_plot, plot = C_gga_EDPf_pca)

# relative paths from wd for outputs headed for storage in "supplement/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = 150, height = 60,
       filename = ofs_hth_nEDPf_pca_plot, plot = C_gga_nEDPf_pca)

# relative paths for outputs headed for storage in the "vault/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = 75, height = 120,
       filename = C_ofv_hth_EDPf_pca_plot, plot = C_gga_EDPf_pca)

ggsave(device = "pdf", dpi = 600, units = "mm", width = 150, height = 60,
       filename = C_ofv_hth_nEDPf_pca_plot, plot = C_gga_nEDPf_pca)

# create a list of objects from section and save them as a workspace file
C_obj <- c(ls(pattern = "C_"), "A_hth_abs_EDPf", "A_hth_abs_nEDPf")
save(list = C_obj, file = C_ofv_WS)

### ************************************
### SECTION D - differential abundance testing ----
### ************************************

# purpose:
# this section tests for differentially abundant features (SVs) between...
# ... diet groups (e.g. Control vs Rice bran; Control vs Fermented rice bran...
# ... and Rice bran vs Fermented rice bran)

# main/supplementary items generated:
# Table S3

### ************************************
### D - STEP 1 - format for ALDEx2 ----
### ************************************

# create a new version of the input df from Section A
D_hth_abs_EDP_0 <- A_hth_abs_EDP

# create df with taxonomic classifications for fts (useful downstream)
D_hth_EDP_tax <- dplyr::select(D_hth_abs_EDP_0, FeatureID,
                               dplyr::contains("ggs"),
                               dplyr::contains("slv"))

# remove unneeded columns
D_hth_abs_EDP_1 <- dplyr::select(D_hth_abs_EDP_0, -RepSeq,
                                 -dplyr::contains("ggs"),
                                 -dplyr::contains("slv"))

# convert FeatureID to row.names and remove FeatureID column
D_hth_abs_EDP_2 <- as.data.frame(D_hth_abs_EDP_1)
row.names(D_hth_abs_EDP_2) <- D_hth_abs_EDP_2$FeatureID
D_hth_abs_EDP_3 <- dplyr::select(D_hth_abs_EDP_2, -FeatureID)

# transpose to convert samples to rows and fts to cols
D_hth_abs_EDP_4 <- as.data.frame(t(D_hth_abs_EDP_3))

# create column from SampleIDs
D_hth_abs_EDP_5 <- D_hth_abs_EDP_4
D_hth_abs_EDP_5$SampleID <- row.names(D_hth_abs_EDP_5)

# read in metadata file and merge with the above df
D_met <- read.table(ifp_met, header = TRUE, sep = "\t", as.is = TRUE)

D_hth_abs_EDP_met_0 <- merge(x = D_hth_abs_EDP_5, y = D_met,
                             all = FALSE, sort = FALSE, by = "SampleID")

# create dfs for pairwise comparisons of diet groups
# convert column DietID into row.names
D_hth_abs_EDP_met_1 <- D_hth_abs_EDP_met_0
row.names(D_hth_abs_EDP_met_1) <- D_hth_abs_EDP_met_1$DietID

# remove uneeded columns
D_col_rmv <- c("SampleID", "Cohort", "CohortOneLetter", "CohortID",
               "Diet", "DietOneLetter", "DietID",
               "Type", "TypeOneLetter", "TypeID", "ID")

D_hth_abs_EDP_met_2 <- dplyr::select(D_hth_abs_EDP_met_1,
                                     -dplyr::one_of(D_col_rmv))

# transpose to convert samples to rowns and fts to cols
D_hth_abs_EDP_met_3 <- as.data.frame(t(D_hth_abs_EDP_met_2))

# filter to isolate pairwise conditions
D_hth_EDP_CvR <- dplyr::select(D_hth_abs_EDP_met_3,
                               dplyr::starts_with("C"),
                               dplyr::starts_with("R"))
D_hth_EDP_CvF <- dplyr::select(D_hth_abs_EDP_met_3,
                               dplyr::starts_with("C"),
                               dplyr::starts_with("F"))
D_hth_EDP_RvF <- dplyr::select(D_hth_abs_EDP_met_3,
                               dplyr::starts_with("R"),
                               dplyr::starts_with("F"))

# create vectors defining the pairwise conditions
# NOTE: must be in exact order as columns in the above dfs
D_hth_EDP_CvR_cnd <- c(rep("C", length(grep("C", names(D_hth_EDP_CvR)))),
                       rep("R", length(grep("R", names(D_hth_EDP_CvR)))))

D_hth_EDP_CvF_cnd <- c(rep("C", length(grep("C", names(D_hth_EDP_CvF)))),
                       rep("F", length(grep("F", names(D_hth_EDP_CvF)))))

D_hth_EDP_RvF_cnd <- c(rep("R", length(grep("R", names(D_hth_EDP_RvF)))),
                       rep("F", length(grep("F", names(D_hth_EDP_RvF)))))

### ************************************
### D - STEP 2 - perform tests ----
### ************************************

# aldex clr with monte carlo simulations
D_hth_EDP_CvR_clr <- aldex.clr(D_hth_EDP_CvR, conds = D_hth_EDP_CvR_cnd,
                               mc.samples = 1000)
D_hth_EDP_CvF_clr <- aldex.clr(D_hth_EDP_CvF, conds = D_hth_EDP_CvF_cnd,
                               mc.samples = 1000)
D_hth_EDP_RvF_clr <- aldex.clr(D_hth_EDP_RvF, conds = D_hth_EDP_RvF_cnd,
                               mc.samples = 1000)

# Welch's t-test and Wilcoxon rank sum tests with Benjamani Hochberg correction
D_hth_EDP_CvR_sts <- aldex.ttest(D_hth_EDP_CvR_clr,
                                 conditions = D_hth_EDP_CvR_cnd)
D_hth_EDP_CvF_sts <- aldex.ttest(D_hth_EDP_CvF_clr,
                                 conditions = D_hth_EDP_CvF_cnd)
D_hth_EDP_RvF_sts <- aldex.ttest(D_hth_EDP_RvF_clr,
                                 conditions = D_hth_EDP_RvF_cnd)

# calculate effect sizes for each mc instance, report the expected value
D_hth_EDP_CvR_eff <- aldex.effect(D_hth_EDP_CvR_clr,
                                  conditions = D_hth_EDP_CvR_cnd,
                                  include.sample.summary = TRUE)
D_hth_EDP_CvF_eff <- aldex.effect(D_hth_EDP_CvF_clr,
                                  conditions = D_hth_EDP_CvF_cnd,
                                  include.sample.summary = TRUE)
D_hth_EDP_RvF_eff <- aldex.effect(D_hth_EDP_RvF_clr,
                                  conditions = D_hth_EDP_RvF_cnd,
                                  include.sample.summary = TRUE)

### ************************************
### D - STEP 3 - format results ----
### ************************************

# combine stats and effects dfs
D_hth_EDP_CvR_rst_0 <- data.frame(D_hth_EDP_CvR_sts, D_hth_EDP_CvR_eff)
D_hth_EDP_CvF_rst_0 <- data.frame(D_hth_EDP_CvF_sts, D_hth_EDP_CvF_eff)
D_hth_EDP_RvF_rst_0 <- data.frame(D_hth_EDP_RvF_sts, D_hth_EDP_RvF_eff)

# add back in taxonomic info
D_hth_EDP_CvR_rst_1 <- D_hth_EDP_CvR_rst_0
D_hth_EDP_CvF_rst_1 <- D_hth_EDP_CvF_rst_0
D_hth_EDP_RvF_rst_1 <- D_hth_EDP_RvF_rst_0

D_hth_EDP_CvR_rst_1$FeatureID <- row.names(D_hth_EDP_CvR_rst_1)
D_hth_EDP_CvF_rst_1$FeatureID <- row.names(D_hth_EDP_CvF_rst_1)
D_hth_EDP_RvF_rst_1$FeatureID <- row.names(D_hth_EDP_RvF_rst_1)

D_hth_EDP_CvR_rst_2 <- dplyr::select(D_hth_EDP_CvR_rst_1, FeatureID,
                                     dplyr::everything())
D_hth_EDP_CvF_rst_2 <- dplyr::select(D_hth_EDP_CvF_rst_1, FeatureID,
                                     dplyr::everything())
D_hth_EDP_RvF_rst_2 <- dplyr::select(D_hth_EDP_RvF_rst_1, FeatureID,
                                     dplyr::everything())

D_hth_EDP_CvR_rst_tax <- merge(x = D_hth_EDP_CvR_rst_2, y = D_hth_EDP_tax,
                               all.y = FALSE, sort = FALSE, by = "FeatureID")
D_hth_EDP_CvF_rst_tax <- merge(x = D_hth_EDP_CvF_rst_2, y = D_hth_EDP_tax,
                               all.y = FALSE, sort = FALSE, by = "FeatureID")
D_hth_EDP_RvF_rst_tax <- merge(x = D_hth_EDP_RvF_rst_2, y = D_hth_EDP_tax,
                               all.y = FALSE, sort = FALSE, by = "FeatureID")

# filter to retain Features with wi.eBH values less than or equal to 0.1
# wi.eBH henceforth referred to as q value
D_hth_EDP_CvR_q1_sigs_0 <- dplyr::filter(D_hth_EDP_CvR_rst_tax, wi.eBH <= 0.1)
D_hth_EDP_CvF_q1_sigs_0 <- dplyr::filter(D_hth_EDP_CvF_rst_tax, wi.eBH <= 0.1)
D_hth_EDP_RvF_q1_sigs_0 <- dplyr::filter(D_hth_EDP_RvF_rst_tax, wi.eBH <= 0.1)

# count number of differentially abundant features at the specified threshold
# NOTE: numbers may change between runs (see vignette for ALDEx2)
D_hth_EDP_CvR_q1_sigs_num <- nrow(D_hth_EDP_CvR_q1_sigs_0)
D_hth_EDP_CvF_q1_sigs_num <- nrow(D_hth_EDP_CvF_q1_sigs_0)
D_hth_EDP_RvF_q1_sigs_num <- nrow(D_hth_EDP_RvF_q1_sigs_0)

# format the above tables to obtain a final table with columns:
# "Comparison", "Median log2 fold difference", "Direction of increase",	
# "P", "q", "Greengenes", "SILVA", "FeatureID"

# retain columns of interest
D_hth_EDP_CvR_q1_sigs_1 <- dplyr::select(D_hth_EDP_CvR_q1_sigs_0, 
                                         rab.win.C, rab.win.R, diff.btw, wi.ep, 
                                         wi.eBH, int.ggs.lws.txn, 
                                         int.slv.lws.txn, FeatureID)

D_hth_EDP_CvF_q1_sigs_1 <- dplyr::select(D_hth_EDP_CvF_q1_sigs_0, 
                                         rab.win.C, rab.win.F, diff.btw, wi.ep, 
                                         wi.eBH, int.ggs.lws.txn, 
                                         int.slv.lws.txn, FeatureID)

D_hth_EDP_RvF_q1_sigs_1 <- dplyr::select(D_hth_EDP_RvF_q1_sigs_0, 
                                         rab.win.R, rab.win.F, diff.btw, wi.ep, 
                                         wi.eBH, int.ggs.lws.txn, 
                                         int.slv.lws.txn, FeatureID)

# determine direction of increase (i.e. which diet group was SV increased in?)
# this determination will use values in column 'diff.btw' as follows:
# comparison C vs R: + value = increased in R; - value = increased in C
# comparison C vs F: + value = increased in F; - value = increased in C
# comparison R vs F: + value = increased in R; - value = increased in F
# NOTE: this can be doublechecked/confirmed by looking at rows in 'rab.win.'
D_hth_EDP_CvR_q1_sigs_2 <- D_hth_EDP_CvR_q1_sigs_1
D_hth_EDP_CvF_q1_sigs_2 <- D_hth_EDP_CvF_q1_sigs_1
D_hth_EDP_RvF_q1_sigs_2 <- D_hth_EDP_RvF_q1_sigs_1

D_hth_EDP_CvR_q1_sigs_2$Direction.of.increase <- ifelse(
  D_hth_EDP_CvR_q1_sigs_2[, "diff.btw"] > 0, 
  yes = "Rice bran", no = "Control")

D_hth_EDP_CvF_q1_sigs_2$Direction.of.increase <- ifelse(
  D_hth_EDP_CvF_q1_sigs_2[, "diff.btw"] > 0, 
  yes = "B. longum-fermented rice bran", no = "Control")

D_hth_EDP_RvF_q1_sigs_2$Direction.of.increase <- ifelse(
  D_hth_EDP_RvF_q1_sigs_2[, "diff.btw"] > 0, 
  yes = "Rice bran", no = "B. longum-fermented rice bran")

# create new column specifying the comparison made
# NOTE: this is written to logically follow the direction of increase
# e.g. CvR is written as 'R vs C' to logically follow R - C = + or - value
# where + values are increased in R and negative values are increased in C
D_hth_EDP_CvR_q1_sigs_3 <- D_hth_EDP_CvR_q1_sigs_2
D_hth_EDP_CvF_q1_sigs_3 <- D_hth_EDP_CvF_q1_sigs_2
D_hth_EDP_RvF_q1_sigs_3 <- D_hth_EDP_RvF_q1_sigs_2

D_hth_EDP_CvR_q1_sigs_3$Comparison <- "R vs C"
D_hth_EDP_CvF_q1_sigs_3$Comparison <- "F vs C"
D_hth_EDP_RvF_q1_sigs_3$Comparison <- "R vs F"

# reorder, rename, and drop unneeded columns
D_hth_EDP_CvR_q1_sigs_4 <- dplyr::select(D_hth_EDP_CvR_q1_sigs_3, 
                                         Comparison,
                                         Median.log2.fold.difference = diff.btw,
                                         Direction.of.increase, 
                                         P = wi.ep, q = wi.eBH, 
                                         Greengenes = int.ggs.lws.txn, 
                                         SILVA = int.slv.lws.txn, 
                                         FeatureID)

D_hth_EDP_CvF_q1_sigs_4 <- dplyr::select(D_hth_EDP_CvF_q1_sigs_3, 
                                         Comparison,
                                         Median.log2.fold.difference = diff.btw,
                                         Direction.of.increase, 
                                         P = wi.ep, q = wi.eBH, 
                                         Greengenes = int.ggs.lws.txn, 
                                         SILVA = int.slv.lws.txn, 
                                         FeatureID)

D_hth_EDP_RvF_q1_sigs_4 <- dplyr::select(D_hth_EDP_RvF_q1_sigs_3, 
                                         Comparison,
                                         Median.log2.fold.difference = diff.btw,
                                         Direction.of.increase, 
                                         P = wi.ep, q = wi.eBH, 
                                         Greengenes = int.ggs.lws.txn, 
                                         SILVA = int.slv.lws.txn, 
                                         FeatureID)

# sort the above dfs by direction of increase
D_hth_EDP_CvR_q1_sigs_5 <- D_hth_EDP_CvR_q1_sigs_4[
  order(D_hth_EDP_CvR_q1_sigs_4$Direction.of.increase, decreasing = TRUE), ]

D_hth_EDP_CvF_q1_sigs_5 <- D_hth_EDP_CvF_q1_sigs_4[
  order(D_hth_EDP_CvF_q1_sigs_4$Direction.of.increase), ]

D_hth_EDP_RvF_q1_sigs_5 <- D_hth_EDP_RvF_q1_sigs_4[
  order(D_hth_EDP_RvF_q1_sigs_4$Direction.of.increase, decreasing = TRUE), ]

# combine above formatted dfs
D_EDP_q1_sigs_0 <- rbind(D_hth_EDP_CvR_q1_sigs_5, 
                         D_hth_EDP_CvF_q1_sigs_5,
                         D_hth_EDP_RvF_q1_sigs_5)

# round all numeric values the three decimal places
D_EDP_q1_sigs_1 <- dplyr::mutate_if(D_EDP_q1_sigs_0, is.numeric, 
                                    round, digits = 3)

# format columns prior to output
D_EDP_q1_sigs <- D_EDP_q1_sigs_1

names(D_EDP_q1_sigs) <- gsub(pattern = "\\.", replacement = " ", 
                             x = names(D_EDP_q1_sigs))

### ************************************
### D - WRITE OUTPUTS ----
### ************************************

# print the number of differentially abundundant SVs from above
print(D_hth_EDP_CvR_q1_sigs_num) # 29 (sometimes 30)
print(D_hth_EDP_CvF_q1_sigs_num) # 58 (sometimes 56 or 57)
print(D_hth_EDP_RvF_q1_sigs_num) # 2 (was 2 for every instance I ran)

# relative paths from wd for outputs headed for storage in "supplement/"
write.table(sep = "\t", row.names = FALSE, x = D_EDP_q1_sigs,
            file = ofs_D_EDP_q1_sigs_tab)

# relative paths for outputs headed for storage in the "vault/"
write.table(sep = "\t", row.names = FALSE, x = D_EDP_q1_sigs,
            file = D_ofv_D_EDP_q1_sigs)

# create a list of objects from section and save them as a workspace file
D_obj <- c(ls(pattern = "D_"), "A_hth_abs_EDP")
save(list = D_obj, file = D_ofv_WS)

### ************************************
### SECTION E - visualize features of interest ----
### ************************************

# purpose:
# this section creates ggplots for differentially abundant SVs... 
# of potential interest

# main/supplementary items generated:
# Figure 4
# Figure S3

# NOTE: for a dplyr filtering step, FeatureIDs are explicilty referenced
# if starting from Code_S1.R which means you ran through QIIME
# and you will not have the same FeatureIDs
# as a result, the code in this section WILL NOT work
# it is written in my notes to alter this in later versions; however,
# if you are reading this note then this file has not yet been updated

### ************************************
### E - STEP 1 - format ALDEx2 results for ggplot ----
### ************************************

# create new versions of the input dfs from Section D
E_hth_EDP_RvF_q1_sigs_0 <- D_hth_EDP_RvF_q1_sigs_0
E_hth_EDP_CvR_q1_sigs_0 <- D_hth_EDP_CvR_q1_sigs_0
E_hth_EDP_CvF_q1_sigs_0 <- D_hth_EDP_CvF_q1_sigs_0

# reorder, rename, and drop unneeded columns
E_hth_EDP_RvF_q1_sigs_1 <- dplyr::select(E_hth_EDP_RvF_q1_sigs_0, 
                                         FeatureID, q = wi.eBH, 
                                         rab.win.R, rab.win.F, 
                                         dplyr::contains("sample"), 
                                         ggs.txn = int.ggs.lws.txn, 
                                         slv.txn = int.slv.lws.txn)

E_hth_EDP_CvR_q1_sigs_1 <- dplyr::select(E_hth_EDP_CvR_q1_sigs_0, 
                                         FeatureID, q = wi.eBH, 
                                         rab.win.C, rab.win.R, 
                                         dplyr::contains("sample"), 
                                         ggs.txn = int.ggs.lws.txn, 
                                         slv.txn = int.slv.lws.txn)

E_hth_EDP_CvF_q1_sigs_1 <- dplyr::select(E_hth_EDP_CvF_q1_sigs_0, 
                                         FeatureID, q = wi.eBH, 
                                         rab.win.C, rab.win.F, 
                                         dplyr::contains("sample"), 
                                         ggs.txn = int.ggs.lws.txn, 
                                         slv.txn = int.slv.lws.txn)

# for CvR and CvF, filter to retain the six Features increased in R and F
E_hth_EDP_CvR_q1_sigs_2 <- dplyr::filter(
  E_hth_EDP_CvR_q1_sigs_1,
  FeatureID == "2d0b61d3962846ad597598bd8657bdc3" |
    FeatureID == "73f14d03d440b253e6e671db77fb578a" |
    FeatureID == "f03b21b0d3b38bc8b7ab248a02d1073e" |
    FeatureID == "243a725661674fde1c0c2afd5f24393d" |
    FeatureID == "15ac133af11f957ad7bf0943b14a4993" |
    FeatureID == "2d0559e9488391eb0918ff0a97e14037")

E_hth_EDP_CvF_q1_sigs_2 <- dplyr::filter(
  E_hth_EDP_CvF_q1_sigs_1,
  FeatureID == "2d0b61d3962846ad597598bd8657bdc3" |
    FeatureID == "73f14d03d440b253e6e671db77fb578a" |
    FeatureID == "f03b21b0d3b38bc8b7ab248a02d1073e" |
    FeatureID == "243a725661674fde1c0c2afd5f24393d" |
    FeatureID == "15ac133af11f957ad7bf0943b14a4993" |
    FeatureID == "2d0559e9488391eb0918ff0a97e14037")


# create two new dfs, one for plotting group data, one for plotting sample data
E_RvF_grp_0 <- dplyr::select(E_hth_EDP_RvF_q1_sigs_1, 
                             -dplyr::contains("sample"))
E_CvR_grp_0 <- dplyr::select(E_hth_EDP_CvR_q1_sigs_2, 
                             -dplyr::contains("sample"))
E_CvF_grp_0 <- dplyr::select(E_hth_EDP_CvF_q1_sigs_2, 
                             -dplyr::contains("sample"))

E_RvF_smp_0 <- dplyr::select(E_hth_EDP_RvF_q1_sigs_1, 
                             -dplyr::contains("win"), -q)
E_CvR_smp_0 <- dplyr::select(E_hth_EDP_CvR_q1_sigs_2, 
                             -dplyr::contains("win"), -q)
E_CvF_smp_0 <- dplyr::select(E_hth_EDP_CvF_q1_sigs_2, 
                             -dplyr::contains("win"), -q)

# format dfs to remove 'rab.win.' and 'rab.sample.'
E_RvF_grp_1 <- E_RvF_grp_0
E_CvR_grp_1 <- E_CvR_grp_0
E_CvF_grp_1 <- E_CvF_grp_0

E_RvF_smp_1 <- E_RvF_smp_0
E_CvR_smp_1 <- E_CvR_smp_0
E_CvF_smp_1 <- E_CvF_smp_0

names(E_RvF_grp_1) <- gsub(pattern = "rab.win.", replacement = "", 
                           x = names(E_RvF_grp_1))
names(E_CvR_grp_1) <- gsub(pattern = "rab.win.", replacement = "", 
                           x = names(E_CvR_grp_1))
names(E_CvF_grp_1) <- gsub(pattern = "rab.win.", replacement = "", 
                           x = names(E_CvF_grp_1))

names(E_RvF_smp_1) <- gsub(pattern = "rab.sample.", replacement = "", 
                           x = names(E_RvF_smp_1))
names(E_CvR_smp_1) <- gsub(pattern = "rab.sample.", replacement = "", 
                           x = names(E_CvR_smp_1))
names(E_CvF_smp_1) <- gsub(pattern = "rab.sample.", replacement = "", 
                           x = names(E_CvF_smp_1))

# reshape dfs to convert columns with sample information into rows...
# ... with correseponding column data for each sample
E_id_vars <- c("FeatureID", "ggs.txn", "slv.txn", "q")

E_RvF_grp_rsh_0 <- reshape2::melt(E_RvF_grp_1, id.vars = E_id_vars, 
                                  variable.name = "DietOneLetter", 
                                  value.name = "med.clr")
E_CvF_grp_rsh_0 <- reshape2::melt(E_CvF_grp_1, id.vars = E_id_vars, 
                                  variable.name = "DietOneLetter", 
                                  value.name = "med.clr")
E_CvR_grp_rsh_0 <- reshape2::melt(E_CvR_grp_1, id.vars = E_id_vars, 
                                  variable.name = "DietOneLetter", 
                                  value.name = "med.clr")

E_RvF_smp_rsh <- reshape2::melt(E_RvF_smp_1, id.vars = E_id_vars[1:3], 
                                variable.name = "DietID", 
                                value.name = "med.clr")
E_CvR_smp_rsh <- reshape2::melt(E_CvR_smp_1, id.vars = E_id_vars[1:3], 
                                variable.name = "DietID", 
                                value.name = "med.clr")
E_CvF_smp_rsh <- reshape2::melt(E_CvF_smp_1, id.vars = E_id_vars[1:3], 
                                variable.name = "DietID", 
                                value.name = "med.clr")

# read in metadata file, retain needed columns, and merge w/ reshaped sample dfs
E_met_0 <- read.table(ifp_met, header = TRUE, sep = "\t", as.is = TRUE)

E_met_1 <- dplyr::select(E_met_0, DietID, DietOneLetter, TypeOneLetter)

E_RvF_smp_met_0 <- merge(x = E_met_1, y = E_RvF_smp_rsh, by = "DietID",
                         all = FALSE, sort = TRUE)
E_CvR_smp_met_0 <- merge(x = E_met_1, y = E_CvR_smp_rsh, by = "DietID",
                         all = FALSE, sort = TRUE)
E_CvF_smp_met_0 <- merge(x = E_met_1, y = E_CvF_smp_rsh, by = "DietID",
                         all = FALSE, sort = TRUE)

# add in x-axis scale for plotting
E_RvF_grp_rsh_1 <- E_RvF_grp_rsh_0
E_RvF_smp_met_1 <- E_RvF_smp_met_0
E_CvR_grp_rsh_1 <- E_CvR_grp_rsh_0
E_CvF_grp_rsh_1 <- E_CvF_grp_rsh_0
E_CvR_smp_met_1 <- E_CvR_smp_met_0
E_CvF_smp_met_1 <- E_CvF_smp_met_0

E_RvF_grp_rsh_1$Scale[E_RvF_grp_rsh_1$DietOneLetter == "R"] <- 1
E_RvF_smp_met_1$Scale[E_RvF_smp_met_1$DietOneLetter == "R"] <- 1
E_RvF_grp_rsh_1$Scale[E_RvF_grp_rsh_1$DietOneLetter == "F"] <- 1
E_RvF_smp_met_1$Scale[E_RvF_smp_met_1$DietOneLetter == "F"] <- 1

E_CvR_grp_rsh_1$Scale[E_CvR_grp_rsh_1$DietOneLetter == "C"] <- 1
E_CvR_smp_met_1$Scale[E_CvR_smp_met_1$DietOneLetter == "C"] <- 1
E_CvR_grp_rsh_1$Scale[E_CvR_grp_rsh_1$DietOneLetter == "R"] <- 1
E_CvR_smp_met_1$Scale[E_CvR_smp_met_1$DietOneLetter == "R"] <- 1

E_CvF_grp_rsh_1$Scale[E_CvF_grp_rsh_1$DietOneLetter == "C"] <- 2
E_CvF_smp_met_1$Scale[E_CvF_smp_met_1$DietOneLetter == "C"] <- 2
E_CvF_grp_rsh_1$Scale[E_CvF_grp_rsh_1$DietOneLetter == "F"] <- 2
E_CvF_smp_met_1$Scale[E_CvF_smp_met_1$DietOneLetter == "F"] <- 2

# create new dfs with individual features
E_RvF_grp_ft1 <- dplyr::filter(E_RvF_grp_rsh_1, 
                               FeatureID == "6bc1c99d965f0989dfcfcfeacf57ccc7")
E_RvF_smp_ft1 <- dplyr::filter(E_RvF_smp_met_1, 
                               FeatureID == "6bc1c99d965f0989dfcfcfeacf57ccc7")
E_RvF_grp_ft2 <- dplyr::filter(E_RvF_grp_rsh_1, 
                               FeatureID == "083e6b5d9d344047f0553b3052e0ccf3")
E_RvF_smp_ft2 <- dplyr::filter(E_RvF_smp_met_1, 
                               FeatureID == "083e6b5d9d344047f0553b3052e0ccf3")

E_CvR_grp_ft1 <- dplyr::filter(E_CvR_grp_rsh_1, 
                               FeatureID == "2d0b61d3962846ad597598bd8657bdc3")
E_CvR_smp_ft1 <- dplyr::filter(E_CvR_smp_met_1, 
                               FeatureID == "2d0b61d3962846ad597598bd8657bdc3")
E_CvR_grp_ft2 <- dplyr::filter(E_CvR_grp_rsh_1, 
                               FeatureID == "73f14d03d440b253e6e671db77fb578a") 
E_CvR_smp_ft2 <- dplyr::filter(E_CvR_smp_met_1, 
                               FeatureID == "73f14d03d440b253e6e671db77fb578a") 
E_CvR_grp_ft3 <- dplyr::filter(E_CvR_grp_rsh_1, 
                               FeatureID == "f03b21b0d3b38bc8b7ab248a02d1073e")
E_CvR_smp_ft3 <- dplyr::filter(E_CvR_smp_met_1, 
                               FeatureID == "f03b21b0d3b38bc8b7ab248a02d1073e")
E_CvR_grp_ft4 <- dplyr::filter(E_CvR_grp_rsh_1, 
                               FeatureID == "243a725661674fde1c0c2afd5f24393d")
E_CvR_smp_ft4 <- dplyr::filter(E_CvR_smp_met_1, 
                               FeatureID == "243a725661674fde1c0c2afd5f24393d")
E_CvR_grp_ft5 <- dplyr::filter(E_CvR_grp_rsh_1, 
                               FeatureID == "15ac133af11f957ad7bf0943b14a4993")
E_CvR_smp_ft5 <- dplyr::filter(E_CvR_smp_met_1, 
                               FeatureID == "15ac133af11f957ad7bf0943b14a4993")
E_CvR_grp_ft6 <- dplyr::filter(E_CvR_grp_rsh_1, 
                               FeatureID == "2d0559e9488391eb0918ff0a97e14037")
E_CvR_smp_ft6 <- dplyr::filter(E_CvR_smp_met_1, 
                               FeatureID == "2d0559e9488391eb0918ff0a97e14037")


E_CvF_grp_ft1 <- dplyr::filter(E_CvF_grp_rsh_1, 
                               FeatureID == "2d0b61d3962846ad597598bd8657bdc3")
E_CvF_smp_ft1 <- dplyr::filter(E_CvF_smp_met_1, 
                               FeatureID == "2d0b61d3962846ad597598bd8657bdc3")
E_CvF_grp_ft2 <- dplyr::filter(E_CvF_grp_rsh_1, 
                               FeatureID == "73f14d03d440b253e6e671db77fb578a") 
E_CvF_smp_ft2 <- dplyr::filter(E_CvF_smp_met_1, 
                               FeatureID == "73f14d03d440b253e6e671db77fb578a") 
E_CvF_grp_ft3 <- dplyr::filter(E_CvF_grp_rsh_1, 
                               FeatureID == "f03b21b0d3b38bc8b7ab248a02d1073e")
E_CvF_smp_ft3 <- dplyr::filter(E_CvF_smp_met_1, 
                               FeatureID == "f03b21b0d3b38bc8b7ab248a02d1073e")
E_CvF_grp_ft4 <- dplyr::filter(E_CvF_grp_rsh_1, 
                               FeatureID == "243a725661674fde1c0c2afd5f24393d")
E_CvF_smp_ft4 <- dplyr::filter(E_CvF_smp_met_1, 
                               FeatureID == "243a725661674fde1c0c2afd5f24393d")
E_CvF_grp_ft5 <- dplyr::filter(E_CvF_grp_rsh_1, 
                               FeatureID == "15ac133af11f957ad7bf0943b14a4993")
E_CvF_smp_ft5 <- dplyr::filter(E_CvF_smp_met_1, 
                               FeatureID == "15ac133af11f957ad7bf0943b14a4993")
E_CvF_grp_ft6 <- dplyr::filter(E_CvF_grp_rsh_1, 
                               FeatureID == "2d0559e9488391eb0918ff0a97e14037")
E_CvF_smp_ft6 <- dplyr::filter(E_CvF_smp_met_1, 
                               FeatureID == "2d0559e9488391eb0918ff0a97e14037")

# for CvR and CvF, combine dfs for the same features
E_CvR_CvF_grp_ft1 <- rbind(E_CvR_grp_ft1, E_CvF_grp_ft1)
E_CvR_CvF_smp_ft1 <- rbind(E_CvR_smp_ft1, E_CvF_smp_ft1)
E_CvR_CvF_grp_ft2 <- rbind(E_CvR_grp_ft2, E_CvF_grp_ft2)
E_CvR_CvF_smp_ft2 <- rbind(E_CvR_smp_ft2, E_CvF_smp_ft2)
E_CvR_CvF_grp_ft3 <- rbind(E_CvR_grp_ft3, E_CvF_grp_ft3)
E_CvR_CvF_smp_ft3 <- rbind(E_CvR_smp_ft3, E_CvF_smp_ft3)
E_CvR_CvF_grp_ft4 <- rbind(E_CvR_grp_ft4, E_CvF_grp_ft4)
E_CvR_CvF_smp_ft4 <- rbind(E_CvR_smp_ft4, E_CvF_smp_ft4)
E_CvR_CvF_grp_ft5 <- rbind(E_CvR_grp_ft5, E_CvF_grp_ft5)
E_CvR_CvF_smp_ft5 <- rbind(E_CvR_smp_ft5, E_CvF_smp_ft5)
E_CvR_CvF_grp_ft6 <- rbind(E_CvR_grp_ft6, E_CvF_grp_ft6)
E_CvR_CvF_smp_ft6 <- rbind(E_CvR_smp_ft6, E_CvF_smp_ft6)

### ************************************
### E - STEP 2 - format for ggplots ----
### ************************************

# define shape, color, and linetype codes
# plot order is: D E P
E_shp_type <- c(19, 17, 1)

E_hex_FR <- c("#e66101", "#1f78b4")
E_hex_CFR <- c("#33a02c", "#e66101", "#1f78b4")
E_hex_RC_FC <- c("#1f78b4", "#33a02c", "#e66101", "#33a02c")

E_lne_FR <- c(1, 2)
E_lne_RC_FC <- c(2, 4, 1, 4)

# hex code for black and grey
E_blk <- "#000000"
E_gry <- "#bbbbbb"

# font family
E_fnt_fam <- "Courier"

# axis parameters
E_x_lim_RvF <- c(0.82, 1.18)
E_y_lim_RvF <- c(-10.5, 10.5)
E_y_brk_RvF <- c(-9:9)
E_y_lab_RvF <- c("", "-8", "", "-6", "", "-4", "", "-2", "", "0", 
                 "", "2", "", "4", "", "6", "", "8", "")

E_x_lim_CvR_CvF <- c(0.5, 2.5)
E_y_lim_CvR_CvF <- c(-13.5, 13.5)
E_y_brk_CvR_CvF <- c(-12:12)
E_y_lab_CvR_CvF <- c("-12", "", "-10", "", "-8", "", "-6", "", 
                     "-4", "", "-2", "", "0", "", "2", "", 
                     "4", "", "6", "", "8", "", "10", "", "12")

# create dfs to plot lines (horizontal at 0, q value, line above taxon)
E_seg_x_RvF <- data.frame(xb = E_x_lim_RvF[1], yb = 0,
                          xe = E_x_lim_RvF[2], ye = 0)
E_seg_q_RvF <- data.frame(xb = 0.9, yb = E_y_lim_RvF[2] - 0.9, 
                          xe = 1.1, ye = E_y_lim_RvF[2] - 0.9)
E_seg_t_RvF <- data.frame(xb = E_x_lim_RvF[1], yb = E_y_lim_RvF[1] + 1.2, 
                          xe = E_x_lim_RvF[2], ye = E_y_lim_RvF[1] + 1.2)

E_seg_x_CvR_CvF <- data.frame(xb = E_x_lim_CvR_CvF[1], yb = 0,
                              xe = E_x_lim_CvR_CvF[2], ye = 0)
E_seg_q_CvR_CvF <- data.frame(xb = c(0.8, 1.8), xe = c(1.2, 2.2),
                              yb = E_y_lim_CvR_CvF[2] - 0.9, 
                              ye = E_y_lim_CvR_CvF[2] - 0.9)
E_seg_t_CvR_CvF <- data.frame(xb = E_x_lim_CvR_CvF[1], 
                              xe = E_x_lim_CvR_CvF[2],
                              yb = E_y_lim_CvR_CvF[1] + 1.3, 
                              ye = E_y_lim_CvR_CvF[1] + 1.3)

# round and paste q-values into vectors useful to check the q-value labels
E_RvF_ft1_q_CHECK <- as.character(round(unique(E_RvF_grp_ft1[, "q"]), 
                                        digits = 3))
E_RvF_ft2_q_CHECK <- as.character(round(unique(E_RvF_grp_ft2[, "q"]), 
                                        digits = 3))

E_CvR_CvF_ft1_q_CHECK <- as.character(round(unique(E_CvR_CvF_grp_ft1[, "q"]), 
                                            digits = 3))
E_CvR_CvF_ft2_q_CHECK <- as.character(round(unique(E_CvR_CvF_grp_ft2[, "q"]), 
                                            digits = 3))
E_CvR_CvF_ft3_q_CHECK <- as.character(round(unique(E_CvR_CvF_grp_ft3[, "q"]), 
                                            digits = 3))
E_CvR_CvF_ft4_q_CHECK <- as.character(round(unique(E_CvR_CvF_grp_ft4[, "q"]), 
                                            digits = 3))
E_CvR_CvF_ft5_q_CHECK <- as.character(round(unique(E_CvR_CvF_grp_ft5[, "q"]), 
                                            digits = 3))
E_CvR_CvF_ft6_q_CHECK <- as.character(round(unique(E_CvR_CvF_grp_ft6[, "q"]), 
                                            digits = 3))

# plot and save parameters
E_ggp_txt_size <- 1.9
E_ggp_shp_size <- 1.3
E_ggp_shp_strk <- 0.3
E_shp_jitr_RvF <- 0.11
E_shp_doge_RvF <- 0.3
E_shp_jitr_CvR_CvF <- 0.16
E_shp_doge_CvR_CvF <- 0.84
E_ggp_ebr_size <- 0.6
E_ebr_widt_RvF <- 0.3
E_ebr_widt_CvR_CvF <- 0.8
E_ggp_axs_txt_size <- 6
E_lab_fnt <- list(size = 13, family = "Courier", face = "bold")
E_ggp_marg_RvF <- theme(plot.margin = unit(c(1, 1, 0.5, 1), "mm"))
E_ggp_marg_CvR_CvF <- theme(plot.margin = unit(c(2, 2, 0.75, 2), "mm"))

### ************************************
### E - STEP 3 - create ggplots ----
### ************************************

# RvF panel A
E_ggp_RvF_ft1 <- ggplot(data = E_RvF_smp_ft1, aes(x = Scale, y = med.clr, 
                                                  shape = TypeOneLetter, 
                                                  color = DietOneLetter,
                                                  group = rev(DietOneLetter))) +
  geom_segment(inherit.aes = FALSE, data = E_seg_x_RvF, size = 0.3,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_gry, linetype = 3) +
  geom_errorbar(inherit.aes = FALSE, data = E_RvF_grp_ft1,
                aes(x = Scale, ymin = med.clr, ymax = med.clr,
                    group = DietOneLetter),
                position = position_dodge(width = E_shp_doge_RvF),
                size = E_ggp_ebr_size, width = E_ebr_widt_RvF, alpha = 0.4,
                color = E_hex_FR, linetype = E_lne_FR) +
  geom_jitter(size = E_ggp_shp_size, stroke = E_ggp_shp_strk, 
              position = position_jitterdodge(jitter.width = E_shp_jitr_RvF, 
                                              dodge.width = E_shp_doge_RvF,
                                              seed = 123456789)) +
  geom_segment(inherit.aes = FALSE, data = E_seg_q_RvF, size = 0.25,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_blk, linetype = 1) +
  annotate("text", x = 1, y = E_y_lim_RvF[2] - 0.5, size = E_ggp_txt_size, 
           color = E_blk, family = E_fnt_fam, 
           label = expression(paste(italic("q"), " = 0.076"))) +
  geom_segment(inherit.aes = FALSE, data = E_seg_t_RvF, size = 0.25,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_blk, linetype = 1) +
  annotate("text", x = 1, y = E_y_lim_RvF[1] + 0.88, size = E_ggp_txt_size,
           color = E_blk, family = E_fnt_fam, label = "Lachnospiraceae/") +
  annotate("text", x = 1, y = E_y_lim_RvF[1] + 0.38, size = E_ggp_txt_size,
           color = E_blk, family = E_fnt_fam, 
           label = expression(italic("Roseburia"))) +
  scale_shape_manual(values = E_shp_type) +
  scale_color_manual(values = E_hex_FR) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, color = E_blk),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = E_ggp_axs_txt_size, color = E_blk),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.25, color = E_blk),
        text = element_text(family = E_fnt_fam)) +
  scale_x_continuous(limits = E_x_lim_RvF, expand = c(0, 0)) +
  scale_y_continuous(limits = E_y_lim_RvF, breaks = E_y_brk_RvF, 
                     labels = E_y_lab_RvF, expand = c(0, 0))

# RvF panel B
E_ggp_RvF_ft2 <- ggplot(data = E_RvF_smp_ft2, aes(x = Scale, y = med.clr, 
                                                  shape = TypeOneLetter, 
                                                  color = DietOneLetter,
                                                  group = rev(DietOneLetter))) +
  geom_segment(inherit.aes = FALSE, data = E_seg_x_RvF, size = 0.3,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_gry, linetype = 3) +
  geom_errorbar(inherit.aes = FALSE, data = E_RvF_grp_ft2,
                aes(x = Scale, ymin = med.clr, ymax = med.clr,
                    group = DietOneLetter),
                position = position_dodge(width = E_shp_doge_RvF),
                size = E_ggp_ebr_size, width = E_ebr_widt_RvF, alpha = 0.4,
                color = E_hex_FR, linetype = E_lne_FR) +
  geom_jitter(size = E_ggp_shp_size, stroke = E_ggp_shp_strk, 
              position = position_jitterdodge(jitter.width = E_shp_jitr_RvF, 
                                              dodge.width = E_shp_doge_RvF,
                                              seed = 123456789))  +
  geom_segment(inherit.aes = FALSE, data = E_seg_q_RvF, size = 0.25,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_blk, linetype = 1) +
  annotate("text", x = 1, y = E_y_lim_RvF[2] - 0.5, size = E_ggp_txt_size, 
           color = E_blk, family = E_fnt_fam, 
           label = expression(paste(italic("q"), " = 0.069"))) +
  geom_segment(inherit.aes = FALSE, data = E_seg_t_RvF, size = 0.25,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_blk, linetype = 1) +
  annotate("text", x = 1, y = E_y_lim_RvF[1] + 0.88, size = E_ggp_txt_size,
           color = E_blk, family = E_fnt_fam, 
           label = "Clostridiales/") +
  annotate("text", x = 1, y = E_y_lim_RvF[1] + 0.38, size = E_ggp_txt_size,
           color = E_blk, family = E_fnt_fam, 
           label = "Clostridiales (vadinBB60)") +
  scale_shape_manual(values = E_shp_type) +
  scale_color_manual(values = E_hex_FR) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, color = E_blk),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = E_ggp_axs_txt_size, color = E_blk),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.25, color = E_blk),
        text = element_text(family = E_fnt_fam)) +
  scale_x_continuous(limits = E_x_lim_RvF, expand = c(0, 0)) +
  scale_y_continuous(limits = E_y_lim_RvF, breaks = E_y_brk_RvF, 
                     labels = E_y_lab_RvF, expand = c(0, 0))

# CvR_CvF panel A
E_ggp_CvR_CvF_ft1 <- ggplot(data = E_CvR_CvF_smp_ft1, 
                            aes(x = Scale, y = med.clr, 
                                shape = TypeOneLetter,
                                color = DietOneLetter, 
                                group = DietOneLetter))  +
  geom_segment(inherit.aes = FALSE, data = E_seg_x_CvR_CvF, size = 0.3,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_gry, linetype = 3) +
  geom_errorbar(inherit.aes = FALSE, data = E_CvR_CvF_grp_ft1,
                aes(x = Scale, ymin = med.clr, ymax = med.clr,
                    group = DietOneLetter),
                position = position_dodge(width = E_shp_doge_CvR_CvF),
                size = E_ggp_ebr_size, width = E_ebr_widt_CvR_CvF, alpha = 0.4,
                color = E_hex_RC_FC, linetype = E_lne_RC_FC) +
  geom_jitter(size = E_ggp_shp_size, stroke = E_ggp_shp_strk, 
              position = position_jitterdodge(jitter.width = E_shp_jitr_CvR_CvF, 
                                              dodge.width = E_shp_doge_CvR_CvF,
                                              seed = 123456789)) +
  geom_segment(inherit.aes = FALSE, data = E_seg_q_CvR_CvF, size = 0.25,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_blk, linetype = 1) +
  annotate("text", x = c(1, 2), y = E_y_lim_CvR_CvF[2] - 0.5, 
           size = E_ggp_txt_size, color = E_blk, family = E_fnt_fam, 
           label = c(expression(paste(italic("q"), " = 0.042")), 
                     expression(paste(italic("q"), " = 0.038")))) +
  geom_segment(inherit.aes = FALSE, data = E_seg_t_CvR_CvF, size = 0.25,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_blk, linetype = 1) +
  annotate("text", x = 1.5, y = E_y_lim_CvR_CvF[1] + 0.93, 
           size = E_ggp_txt_size,
           color = E_blk, family = E_fnt_fam, 
           label = "Lachnospiraceae/") +
  annotate("text", x = 1.5, y = E_y_lim_CvR_CvF[1] + 0.45,
           size = E_ggp_txt_size, color = E_blk, family = E_fnt_fam, 
           label = expression(italic("Roseburia"))) +
  scale_shape_manual(values = E_shp_type) +
  scale_color_manual(values = E_hex_CFR) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, color = E_blk),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = E_ggp_axs_txt_size, color = E_blk),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.25, color = E_blk),
        text = element_text(family = E_fnt_fam)) +
  scale_x_continuous(limits = E_x_lim_CvR_CvF, expand = c(0, 0)) +
  scale_y_continuous(limits = E_y_lim_CvR_CvF, breaks = E_y_brk_CvR_CvF, 
                     labels = E_y_lab_CvR_CvF, expand = c(0, 0))

# CvR_CvF panel B
E_ggp_CvR_CvF_ft2 <- ggplot(data = E_CvR_CvF_smp_ft2, 
                            aes(x = Scale, y = med.clr, 
                                shape = TypeOneLetter,
                                color = DietOneLetter, 
                                group = DietOneLetter))  +
  geom_segment(inherit.aes = FALSE, data = E_seg_x_CvR_CvF, size = 0.3,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_gry, linetype = 3) +
  geom_errorbar(inherit.aes = FALSE, data = E_CvR_CvF_grp_ft2,
                aes(x = Scale, ymin = med.clr, ymax = med.clr,
                    group = DietOneLetter),
                position = position_dodge(width = E_shp_doge_CvR_CvF),
                size = E_ggp_ebr_size, width = E_ebr_widt_CvR_CvF, alpha = 0.4,
                color = E_hex_RC_FC, linetype = E_lne_RC_FC) +
  geom_jitter(size = E_ggp_shp_size, stroke = E_ggp_shp_strk, 
              position = position_jitterdodge(jitter.width = E_shp_jitr_CvR_CvF, 
                                              dodge.width = E_shp_doge_CvR_CvF,
                                              seed = 123456789)) +
  geom_segment(inherit.aes = FALSE, data = E_seg_q_CvR_CvF, size = 0.25,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_blk, linetype = 1) +
  annotate("text", x = c(1, 2), y = E_y_lim_CvR_CvF[2] - 0.5, 
           size = E_ggp_txt_size, color = E_blk, family = E_fnt_fam, 
           label = c(expression(paste(italic("q"), " = 0.017")), 
                     expression(paste(italic("q"), " = 0.016")))) +
  geom_segment(inherit.aes = FALSE, data = E_seg_t_CvR_CvF, size = 0.25,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_blk, linetype = 1) +
  annotate("text", x = 1.5, y = E_y_lim_CvR_CvF[1] + 0.93, 
           size = E_ggp_txt_size, color = E_blk, family = E_fnt_fam, 
           label = "Lachnospiraceae/") +
  annotate("text", x = 1.5, y = E_y_lim_CvR_CvF[1] + 0.45, 
           size = E_ggp_txt_size, color = E_blk, family = E_fnt_fam, 
           label = "Lachnospiraceae") +
  scale_shape_manual(values = E_shp_type) +
  scale_color_manual(values = E_hex_CFR) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, color = E_blk),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = E_ggp_axs_txt_size, color = E_blk),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.25, color = E_blk),
        text = element_text(family = E_fnt_fam)) +
  scale_x_continuous(limits = E_x_lim_CvR_CvF, expand = c(0, 0)) +
  scale_y_continuous(limits = E_y_lim_CvR_CvF, breaks = E_y_brk_CvR_CvF, 
                     labels = E_y_lab_CvR_CvF, expand = c(0, 0))

# CvR_CvF panel C
E_ggp_CvR_CvF_ft3 <- ggplot(data = E_CvR_CvF_smp_ft3, 
                            aes(x = Scale, y = med.clr, 
                                shape = TypeOneLetter,
                                color = DietOneLetter, 
                                group = DietOneLetter))  +
  geom_segment(inherit.aes = FALSE, data = E_seg_x_CvR_CvF, size = 0.3,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_gry, linetype = 3) +
  geom_errorbar(inherit.aes = FALSE, data = E_CvR_CvF_grp_ft3,
                aes(x = Scale, ymin = med.clr, ymax = med.clr,
                    group = DietOneLetter),
                position = position_dodge(width = E_shp_doge_CvR_CvF),
                size = E_ggp_ebr_size, width = E_ebr_widt_CvR_CvF, alpha = 0.4,
                color = E_hex_RC_FC, linetype = E_lne_RC_FC) +
  geom_jitter(size = E_ggp_shp_size, stroke = E_ggp_shp_strk, 
              position = position_jitterdodge(jitter.width = E_shp_jitr_CvR_CvF, 
                                              dodge.width = E_shp_doge_CvR_CvF,
                                              seed = 123456789)) +
  geom_segment(inherit.aes = FALSE, data = E_seg_q_CvR_CvF, size = 0.25,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_blk, linetype = 1) +
  annotate("text", x = c(1, 2), y = E_y_lim_CvR_CvF[2] - 0.5, 
           size = E_ggp_txt_size, color = E_blk, family = E_fnt_fam, 
           label = c(expression(paste(italic("q"), " = 0.014")), 
                     expression(paste(italic("q"), " = 0.050")))) +
  geom_segment(inherit.aes = FALSE, data = E_seg_t_CvR_CvF, size = 0.25,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_blk, linetype = 1) +
  annotate("text", x = 1.5, y = E_y_lim_CvR_CvF[1] + 0.93, 
           size = E_ggp_txt_size, color = E_blk, family = E_fnt_fam, 
           label = "Clostridiales/") +
  annotate("text", x = 1.5, y = E_y_lim_CvR_CvF[1] + 0.45, 
           size = E_ggp_txt_size, color = E_blk, family = E_fnt_fam, 
           label = "Lachnospiraceae (A2)") +
  scale_shape_manual(values = E_shp_type) +
  scale_color_manual(values = E_hex_CFR) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, color = E_blk),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = E_ggp_axs_txt_size, color = E_blk),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.25, color = E_blk),
        text = element_text(family = E_fnt_fam)) +
  scale_x_continuous(limits = E_x_lim_CvR_CvF, expand = c(0, 0)) +
  scale_y_continuous(limits = E_y_lim_CvR_CvF, breaks = E_y_brk_CvR_CvF, 
                     labels = E_y_lab_CvR_CvF, expand = c(0, 0))

# CvR_CvF panel D
E_ggp_CvR_CvF_ft4 <- ggplot(data = E_CvR_CvF_smp_ft4, 
                            aes(x = Scale, y = med.clr, 
                                shape = TypeOneLetter,
                                color = DietOneLetter, 
                                group = DietOneLetter))  +
  geom_segment(inherit.aes = FALSE, data = E_seg_x_CvR_CvF, size = 0.3,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_gry, linetype = 3) +
  geom_errorbar(inherit.aes = FALSE, data = E_CvR_CvF_grp_ft4,
                aes(x = Scale, ymin = med.clr, ymax = med.clr,
                    group = DietOneLetter),
                position = position_dodge(width = E_shp_doge_CvR_CvF),
                size = E_ggp_ebr_size, width = E_ebr_widt_CvR_CvF, alpha = 0.4,
                color = E_hex_RC_FC, linetype = E_lne_RC_FC) +
  geom_jitter(size = E_ggp_shp_size, stroke = E_ggp_shp_strk, 
              position = position_jitterdodge(jitter.width = E_shp_jitr_CvR_CvF, 
                                              dodge.width = E_shp_doge_CvR_CvF,
                                              seed = 123456789)) +
  geom_segment(inherit.aes = FALSE, data = E_seg_q_CvR_CvF, size = 0.25,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_blk, linetype = 1) +
  annotate("text", x = c(1, 2), y = E_y_lim_CvR_CvF[2] - 0.5, 
           size = E_ggp_txt_size, color = E_blk, family = E_fnt_fam, 
           label = c(expression(paste(italic("q"), " = 0.052")), 
                     expression(paste(italic("q"), " = 0.048")))) +
  geom_segment(inherit.aes = FALSE, data = E_seg_t_CvR_CvF, size = 0.25,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_blk, linetype = 1) +
  annotate("text", x = 1.5, y = E_y_lim_CvR_CvF[1] + 0.93, 
           size = E_ggp_txt_size, color = E_blk, family = E_fnt_fam, 
           label = "Clostridiales/") +
  annotate("text", x = 1.5, y = E_y_lim_CvR_CvF[1] + 0.45, 
           size = E_ggp_txt_size, color = E_blk, family = E_fnt_fam, 
           label = "Lachnospiraceae (NK4B4)") +
  scale_shape_manual(values = E_shp_type) +
  scale_color_manual(values = E_hex_CFR) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, color = E_blk),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = E_ggp_axs_txt_size, color = E_blk),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.25, color = E_blk),
        text = element_text(family = E_fnt_fam)) +
  scale_x_continuous(limits = E_x_lim_CvR_CvF, expand = c(0, 0)) +
  scale_y_continuous(limits = E_y_lim_CvR_CvF, breaks = E_y_brk_CvR_CvF, 
                     labels = E_y_lab_CvR_CvF, expand = c(0, 0))

# CvR_CvF panel E
E_ggp_CvR_CvF_ft5 <- ggplot(data = E_CvR_CvF_smp_ft5, 
                            aes(x = Scale, y = med.clr, 
                                shape = TypeOneLetter,
                                color = DietOneLetter, 
                                group = DietOneLetter))  +
  geom_segment(inherit.aes = FALSE, data = E_seg_x_CvR_CvF, size = 0.3,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_gry, linetype = 3) +
  geom_errorbar(inherit.aes = FALSE, data = E_CvR_CvF_grp_ft5,
                aes(x = Scale, ymin = med.clr, ymax = med.clr,
                    group = DietOneLetter),
                position = position_dodge(width = E_shp_doge_CvR_CvF),
                size = E_ggp_ebr_size, width = E_ebr_widt_CvR_CvF, alpha = 0.4,
                color = E_hex_RC_FC, linetype = E_lne_RC_FC) +
  geom_jitter(size = E_ggp_shp_size, stroke = E_ggp_shp_strk, 
              position = position_jitterdodge(jitter.width = E_shp_jitr_CvR_CvF, 
                                              dodge.width = E_shp_doge_CvR_CvF,
                                              seed = 123456789)) +
  geom_segment(inherit.aes = FALSE, data = E_seg_q_CvR_CvF, size = 0.25,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_blk, linetype = 1) +
  annotate("text", x = c(1, 2), y = E_y_lim_CvR_CvF[2] - 0.5, 
           size = E_ggp_txt_size, color = E_blk, family = E_fnt_fam, 
           label = c(expression(paste(italic("q"), " = 0.014")), 
                     expression(paste(italic("q"), " = 0.058")))) +
  geom_segment(inherit.aes = FALSE, data = E_seg_t_CvR_CvF, size = 0.25,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_blk, linetype = 1) +
  annotate("text", x = 1.5, y = E_y_lim_CvR_CvF[1] + 0.93, 
           size = E_ggp_txt_size, color = E_blk, family = E_fnt_fam, 
           label = "Clostridiales/") +
  annotate("text", x = 1.5, y = E_y_lim_CvR_CvF[1] + 0.45, 
           size = E_ggp_txt_size, color = E_blk, family = E_fnt_fam, 
           label = "Lachnospiraceae (NK4A136)") +
  scale_shape_manual(values = E_shp_type) +
  scale_color_manual(values = E_hex_CFR) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, color = E_blk),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = E_ggp_axs_txt_size, color = E_blk),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.25, color = E_blk),
        text = element_text(family = E_fnt_fam)) +
  scale_x_continuous(limits = E_x_lim_CvR_CvF, expand = c(0, 0)) +
  scale_y_continuous(limits = E_y_lim_CvR_CvF, breaks = E_y_brk_CvR_CvF, 
                     labels = E_y_lab_CvR_CvF, expand = c(0, 0))

# CvR_CvF panel F
E_ggp_CvR_CvF_ft6 <- ggplot(data = E_CvR_CvF_smp_ft6, 
                            aes(x = Scale, y = med.clr, 
                                shape = TypeOneLetter,
                                color = DietOneLetter, 
                                group = DietOneLetter))  +
  geom_segment(inherit.aes = FALSE, data = E_seg_x_CvR_CvF, size = 0.3,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_gry, linetype = 3) +
  geom_errorbar(inherit.aes = FALSE, data = E_CvR_CvF_grp_ft6,
                aes(x = Scale, ymin = med.clr, ymax = med.clr,
                    group = DietOneLetter),
                position = position_dodge(width = E_shp_doge_CvR_CvF),
                size = E_ggp_ebr_size, width = E_ebr_widt_CvR_CvF, alpha = 0.4,
                color = E_hex_RC_FC, linetype = E_lne_RC_FC) +
  geom_jitter(size = E_ggp_shp_size, stroke = E_ggp_shp_strk, 
              position = position_jitterdodge(jitter.width = E_shp_jitr_CvR_CvF, 
                                              dodge.width = E_shp_doge_CvR_CvF,
                                              seed = 123456789)) +
  geom_segment(inherit.aes = FALSE, data = E_seg_q_CvR_CvF, size = 0.25,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_blk, linetype = 1) +
  annotate("text", x = c(1, 2), y = E_y_lim_CvR_CvF[2] - 0.5, 
           size = E_ggp_txt_size, color = E_blk, family = E_fnt_fam, 
           label = c(expression(paste(italic("q"), " = 0.031")), 
                     expression(paste(italic("q"), " = 0.010")))) +
  geom_segment(inherit.aes = FALSE, data = E_seg_t_CvR_CvF, size = 0.25,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = E_blk, linetype = 1) +
  annotate("text", x = 1.5, y = E_y_lim_CvR_CvF[1] + 0.93, 
           size = E_ggp_txt_size, color = E_blk, family = E_fnt_fam, 
           label = "Clostridiales/") +
  annotate("text", x = 1.5, y = E_y_lim_CvR_CvF[1] + 0.45, 
           size = E_ggp_txt_size, color = E_blk, family = E_fnt_fam, 
           label = "Clostridiales (CIEAF 020)") +
  scale_shape_manual(values = E_shp_type) +
  scale_color_manual(values = E_hex_CFR) +
  theme(panel.border = element_rect(fill = NA, size = 0.5, color = E_blk),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = E_ggp_axs_txt_size, color = E_blk),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 0.25, color = E_blk),
        text = element_text(family = E_fnt_fam)) +
  scale_x_continuous(limits = E_x_lim_CvR_CvF, expand = c(0, 0)) +
  scale_y_continuous(limits = E_y_lim_CvR_CvF, breaks = E_y_brk_CvR_CvF, 
                     labels = E_y_lab_CvR_CvF, expand = c(0, 0))

### ************************************
### E - STEP 4 - create ggplot legends ----
### ************************************

# create legend for sample type
# create df
E_lgn_type <- data.frame(x = 1, y = 1, "label" = c("E", "D", "P"))

# create ggplot
E_ggp_lgn_type <- ggplot(data = E_lgn_type, aes(x = x, y = y)) +
  geom_point(size = 1, stroke = 0.3, aes(shape = label), color = E_blk) +
  scale_shape_manual(values = c(17, 19, 1), name = NULL,
                     labels = c("caecum",
                                "colon (distal)",
                                "colon (proximal)")) +
  theme_void(base_family = E_fnt_fam) +
  theme(legend.position = c(0.5, 0.5),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = NA, color = NA),
        legend.key.size = unit(0, "mm"),
        legend.text = element_text(size = 5, color = E_blk, hjust = 0,
                                   margin = margin(-0.6, 0, -0.6, -1.4, 
                                                   "mm"))) +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0))

# create legend for line types (study diet)
# create df
E_lgn_diet_RF <- data.frame(x = c(0.9, 0.95, 0.9, 0.95), y = 1,
                            "color" = c("a", "b"))

# create ggplot
E_ggp_lgn_diet_RF <- ggplot(data = E_lgn_diet_RF, 
                            aes(x = x, y = y, color = color)) +
  geom_line(size = 0.3, aes(linetype = color)) + 
  scale_linetype_manual(values = c(2, 1), name = NULL,
                        labels = c("Rice bran",
                                   expression(paste(italic("B. longum"), 
                                                    "-fermented rice bran")))) +
  scale_color_manual(values = c("#1f78b4", "#e66101"), name = NULL, 
                     labels = c("Rice bran",
                                expression(paste(italic("B. longum"), 
                                                 "-fermented rice bran")))) +
  theme_void(base_family = E_fnt_fam) +
  theme(legend.position = c(0.5, 0.5),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = NA, color = NA),
        legend.key.height = unit(0, "mm"),
        legend.key.width = unit(6, "mm"),
        legend.text = element_text(size = 5, color = E_blk, hjust = 0,
                                   margin = margin(-0.6, 0, -0.6, -1.4, 
                                                   "mm"))) +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0))

# create legend for line types (study diet)
# create df
E_lgn_diet_CFR <- data.frame(x = c(0.9, 0.95, 0.9, 0.95, 0.9, 0.95), y = 1,
                             "color" = c("a", "b", "c"))

# create ggplot
E_ggp_lgn_diet_CFR <- ggplot(data = E_lgn_diet_CFR, 
                             aes(x = x, y = y, color = color)) +
  geom_line(size = 0.3, aes(linetype = color)) + 
  scale_linetype_manual(values = c(4, 2, 1), name = NULL,
                        labels = c("Control",
                                   "Rice bran",
                                   expression(paste(italic("B. longum"), 
                                                    "-fermented rice bran")))) +
  scale_color_manual(values = c("#33a02c", "#1f78b4", "#e66101"), name = NULL, 
                     labels = c("Control",
                                "Rice bran",
                                expression(paste(italic("B. longum"), 
                                                 "-fermented rice bran")))) +
  theme_void(base_family = E_fnt_fam) +
  theme(legend.position = c(0.5, 0.5),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = NA, color = NA),
        legend.key.height = unit(0, "mm"),
        legend.key.width = unit(6, "mm"),
        legend.text = element_text(size = 5, color = E_blk, hjust = 0,
                                   margin = margin(-0.6, 0, -0.6, -1.4, 
                                                   "mm"))) +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2), expand = c(0, 0))

### ************************************
### E - STEP 5 - arrange ggplots into panels ----
### ************************************

# arrange legends
E_gga_RvF_lgn <- ggpubr::ggarrange(E_ggp_lgn_diet_RF, E_ggp_lgn_type,
                                   ncol = 1, nrow = 2, align = "h")
E_gga_CvR_CvF_lgn <- ggpubr::ggarrange(E_ggp_lgn_diet_CFR, E_ggp_lgn_type,
                                       ncol = 1, nrow = 2, align = "h")

# arrange plots and add universal y axis label
# for CvR_CvF arrange three at a time, add universal y, then combine
E_gga_RvF_fts <- ggpubr::ggarrange((E_ggp_RvF_ft1 + E_ggp_marg_RvF),
                                   (E_ggp_RvF_ft2 + E_ggp_marg_RvF),
                                   labels = "AUTO", font.label = E_lab_fnt,
                                   ncol = 2, nrow = 1, widths = c(1, 1),
                                   align = "hv")

E_gga_CvR_CvF_abc <- ggpubr::ggarrange((E_ggp_CvR_CvF_ft1 + E_ggp_marg_CvR_CvF),
                                       (E_ggp_CvR_CvF_ft2 + E_ggp_marg_CvR_CvF),
                                       (E_ggp_CvR_CvF_ft3 + E_ggp_marg_CvR_CvF),
                                       labels = c("A", "B", "C"), 
                                       font.label = E_lab_fnt,
                                       ncol = 3, nrow = 1, widths = c(1, 1, 1),
                                       align = "hv")

E_gga_CvR_CvF_def <- ggpubr::ggarrange((E_ggp_CvR_CvF_ft4 + E_ggp_marg_CvR_CvF),
                                       (E_ggp_CvR_CvF_ft5 + E_ggp_marg_CvR_CvF),
                                       (E_ggp_CvR_CvF_ft6 + E_ggp_marg_CvR_CvF),
                                       labels = c("D", "E", "F"), 
                                       font.label = E_lab_fnt,
                                       ncol = 3, nrow = 1, widths = c(1, 1, 1),
                                       align = "hv")

E_gga_ytt_lab <- ggpubr::text_grob("Median clr value", color = E_blk,
                                   rot = 90, family = E_fnt_fam, size = 7)

E_gga_RvF_fts_plot <- ggpubr::annotate_figure(E_gga_RvF_fts, 
                                              left = E_gga_ytt_lab)

E_gga_CvR_CvF_abc_plot <- ggpubr::annotate_figure(E_gga_CvR_CvF_abc, 
                                                  left = E_gga_ytt_lab)
E_gga_CvR_CvF_def_plot <- ggpubr::annotate_figure(E_gga_CvR_CvF_def, 
                                                  left = E_gga_ytt_lab)

E_gga_CvR_CvF_fts_plot <- ggpubr::ggarrange(E_gga_CvR_CvF_abc_plot, 
                                            E_gga_CvR_CvF_def_plot,
                                            ncol = 1, nrow = 2, 
                                            widths = c(1, 1, 1),
                                            align = "hv")

# arrange plots with legends
E_gga_RvF_fts_plot_lgn <- ggpubr::ggarrange(E_gga_RvF_fts_plot, E_gga_RvF_lgn,
                                            ncol = 1, nrow = 2, 
                                            heights = c(14, 1),
                                            align = "h")

E_gga_CvR_CvF_fts_plot_lgn <- ggpubr::ggarrange(E_gga_CvR_CvF_fts_plot, 
                                                E_gga_CvR_CvF_lgn,
                                                ncol = 1, nrow = 2, 
                                                heights = c(28, 1),
                                                align = "h")

### ************************************
### E - WRITE OUTPUTS ----
### ************************************

# print LOGICALs to check that q values in plots match those in table
print(E_RvF_ft1_q_CHECK)
print(E_RvF_ft2_q_CHECK)
print(E_CvR_CvF_ft1_q_CHECK)
print(E_CvR_CvF_ft2_q_CHECK)
print(E_CvR_CvF_ft3_q_CHECK)
print(E_CvR_CvF_ft4_q_CHECK)
print(E_CvR_CvF_ft5_q_CHECK)
print(E_CvR_CvF_ft6_q_CHECK)

# relative paths from wd for outputs headed for storage in "main/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = 78, height = 78,
       filename = ofm_hth_EDP_RvF_fts_plot, 
       plot = E_gga_RvF_fts_plot_lgn)

# relative paths from wd for outputs headed for storage in "supplement/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = 121, height = 156,
       filename = ofs_hth_EDP_CvR_CvF_fts_plot, 
       plot = E_gga_CvR_CvF_fts_plot_lgn)

# relative paths for outputs headed for storage in the "vault/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = 78, height = 78,
       filename = E_ofv_hth_EDP_RvF_fts, 
       plot = E_gga_RvF_fts_plot_lgn)
ggsave(device = "pdf", dpi = 600, units = "mm", width = 121, height = 156,
       filename = E_ofv_hth_EDP_CvR_CvF_fts, 
       plot = E_gga_CvR_CvF_fts_plot_lgn)

# create a list of objects from section and save them as a workspace file
E_obj <- c(ls(pattern = "E_"),
           "D_hth_EDP_RvF_q1_sigs_0", 
           "D_hth_EDP_CvR_q1_sigs_0",
           "D_hth_EDP_CvF_q1_sigs_0")
save(list = E_obj, file = E_ofv_WS)

### ************************************
### SECTION F - create study design timeline ----
### ************************************

# purpose:
# this section creates a visually meaningful representation of the study design

# main/supplementary items generated:
# Figure 1

## ************************************
### F - STEP 1 - create vectors and dfs ----
### ************************************

# hex code for various colors
F_blk <- "#000000"
F_gry <- "#bdbdbd"

# order is F - R - C
F_hex_diet <- c("#e66101", "#1f78b4", "#33a02c")

# global font family
F_fnt_fam <- "Courier"

F_out_box_size <- 0.5
F_lab_out_box_size <- 2.9
F_inn_box_size <- 0.36
F_inn_box_shd_size <- 0.4
F_shp_size <- 1.3
F_shp_stroke <- 0.36
F_txt_size <- 2.2
F_arw_typ <- grid::arrow(length = grid::unit(0.16, "picas"), type = "closed")

# x and y axis limits
F_x_lim <- c(-2.5, 19)
F_y_lim <- c(0, 3.75)

# x positions for sampling timepoints
F_smp_x <- c(0, 2/7, 2, 6, 10, 14, 15)

# y positions for study diets
F_diet_y <- c(1, 1.25, 1.5)

# outer box vertical lines
F_out_box_vrt <- data.frame(xb = c(F_x_lim[1], F_x_lim[2]),
                            xe = c(F_x_lim[1], F_x_lim[2]),
                            yb = F_y_lim[1], ye = F_y_lim[2])

# outer box horizontal lines and segments
F_out_box_hrz <- data.frame(xb = F_x_lim[1], yb = c(F_y_lim[1], F_y_lim[2]),
                            xe = F_x_lim[2], ye = c(F_y_lim[1], F_y_lim[2]))

# inner box shading for first 48 hours post diet intervention
F_inn_box_shd <- data.frame(xmin = F_smp_x[1], ymin = F_diet_y[1] - 0.25,
                            xmax = F_smp_x[2], ymax = F_diet_y[3] + 0.25)

# inner box horizontal cap lines
F_inn_box_hrz <- data.frame(xb = F_smp_x[1], xe = F_smp_x[7], 
                            yb = c(F_diet_y[1] - 0.25, F_diet_y[3] + 0.25),
                            ye = c(F_diet_y[1] - 0.25, F_diet_y[3] + 0.25))

#inner box vertical lines (excluding diet intervention)
F_inn_box_vrt <- data.frame(xb = F_smp_x[2:7], yb = F_diet_y[1] - 0.25, 
                            xe = F_smp_x[2:7], ye = F_diet_y[3] + 0.45)

# study diet horizontal lines (pre diet intervention)
F_diet_pre <- data.frame(xb = c(F_smp_x[1] - 2.25,  
                                rep(F_smp_x[1] - 1.125, times = 3)),
                         xe = c(F_smp_x[1] - 1.125, 
                                rep(F_smp_x[1], times = 3)),
                         yb = F_diet_y[2],
                         ye = c(F_diet_y[2], F_diet_y))

# study diet horizontal lines (post diet intervention)
F_diet_pst <- data.frame(xb = F_smp_x[1], yb = F_diet_y, 
                         xe = F_smp_x[7], ye = F_diet_y)

# inner box vertical cap line at diet intervention
F_inn_box_cap <- data.frame(xb = F_smp_x[1], 
                            xe = F_smp_x[1],
                            yb = F_diet_y[1] - 0.25,
                            ye = F_diet_y[3] + 0.25)

# labels for week number (time)
F_lab_week_num <- data.frame(x = F_smp_x[2:7], y = F_diet_y[3] + 0.45,
                             label = c("1", "2", "6", "10", "14", "15"))

# label denoting week
F_lab_week <- data.frame(x = F_smp_x[1] - 1.125, y = F_diet_y[3] + 0.45,
                         label = "week:")

# single letters for study diet
F_ltr_diet <- data.frame(x = F_smp_x[7], y = F_diet_y, 
                         label = c("C", "R", "F"))

# sample type collecton times
F_shp_type <- data.frame(x = c(F_smp_x[2:7], rep(F_smp_x[7], times = 3)),
                         y = c(rep(F_diet_y[3] + 0.75, times = 5),
                               F_diet_y[3] + 1.0, F_diet_y[3] + 1.25,
                               F_diet_y[3] + 1.5, F_diet_y[3] + 1.75),
                         shape = c(rep(4, times = 5), 1, 19, 17, 9),
                         color = c(rep("#8c510a", times = 5), "#542788",
                                   "#8073ac", "#2d004b", "#e41a1c"))

# labels for shape legend
F_lab_type <- data.frame(x = F_smp_x[7] + 0.36,
                         y = c(F_diet_y[3] + 0.75, F_diet_y[3] + 1.0, 
                               F_diet_y[3] + 1.25, F_diet_y[3] + 1.5, 
                               F_diet_y[3] + 1.75),
                         label = c("faeces", "colon (proximal)", 
                                   "colon (distal)", "caecum", "blood"),
                         color = c("#8c510a", "#542788", 
                                   "#8073ac", "#2d004b", "#e41a1c"))

# arrows for diet intervention and end of study
F_arw <- data.frame(xb = c(F_smp_x[1], F_smp_x[7]), yb = F_diet_y[1] - 0.60,
                    xe = c(F_smp_x[1], F_smp_x[7]), ye = F_diet_y[1] - 0.33)

# labels for diet intervention and end of study arrows
F_lab_arw <- data.frame(x = c(F_smp_x[1], F_smp_x[7]), y = F_diet_y[1] - 0.75,
                        label = c("diet intervention", "end of study"))

# strain and mouse numbers
F_mouse <- data.frame(x = F_x_lim[1] + 0.25, y = F_y_lim[2] - 0.25,
                      label = "Male BALB/c mice:")

# labels for study diet
F_lab_diet <- data.frame(x = F_x_lim[1] + 0.50, 
                         y = c(F_y_lim[2] - 1.00, F_y_lim[2] - 0.75, 
                               F_y_lim[2] - 0.50))

## ************************************
### F - STEP 2 - format and create ggplot ----
### ************************************

# create ggplot
F_ggp_timeline <- ggplot(data = F_out_box_vrt) +
  # outer box vertical lines
  geom_segment(aes(x = xb, xend = xe, y = yb, yend = ye), size = F_out_box_size,
               color = F_blk, linetype = 1, lineend = "square", 
               linejoin = "mitre") +
  # outer box horizontal lines and segments
  geom_segment(inherit.aes = F, data = F_out_box_hrz,
               aes(x = xb, xend = xe, y = yb, yend = ye), size = F_out_box_size,
               color = F_blk, linetype = 1, lineend = "square", 
               linejoin = "mitre") +
  # inner box shading for first 48 hours post diet intervention
  geom_rect(inherit.aes = F, data = F_inn_box_shd, fill = F_gry,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  # inner box horizontal cap lines
  geom_segment(size = F_inn_box_size,
               inherit.aes = F, data = F_inn_box_hrz,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = F_blk, linetype = 1, lineend = "square", 
               linejoin = "mitre") +
  # inner box vertical lines (excluding diet intervention)
  geom_segment(size = F_inn_box_size,
               inherit.aes = F, data = F_inn_box_vrt,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = F_blk, linetype = 1, lineend = "square", 
               linejoin = "mitre") +
  # study diet horizontal lines (pre diet intervention)
  geom_segment(inherit.aes = F, data = F_diet_pre,
               aes(x = xb, xend = xe, y = yb, yend = ye), size = F_inn_box_size,
               color = F_hex_diet[3], linetype = 4, lineend = "square", 
               linejoin = "mitre") +
  # study diet horizontal lines (post diet intervention)
  geom_segment(inherit.aes = F, data = F_diet_pst,
               aes(x = xb, xend = xe, y = yb, yend = ye), size = F_inn_box_size,
               color = F_hex_diet, linetype = c(1, 2, 4), lineend = "square", 
               linejoin = "mitre") +
  # inner box vertical cap line at diet intervention
  geom_segment(size = F_inn_box_size,
               inherit.aes = F, data = F_inn_box_cap,
               aes(x = xb, xend = xe, y = yb, yend = ye),
               color = F_blk, linetype = 1, lineend = "square", 
               linejoin = "mitre") +
  # rectangle border for week number labels - first plot black
  geom_point(inherit.aes = F, data = F_lab_week_num,
             aes(x = x, y = y), size = 3.25, stroke = 0,
             shape = 15, color = F_blk) +
  # rectangle border for week number labels - second plot blank
  geom_point(inherit.aes = F, data = F_lab_week_num,
             aes(x = x, y = y), size = 2.5, stroke = 0,
             shape = 15, color = "white") +
  # labels for week number (time)
  geom_text(inherit.aes = F, data = F_lab_week_num,
            aes(x = x, y = y, label = label), size = 1.75,
            color = F_blk, family = F_fnt_fam) +
  # label denoting week
  geom_text(inherit.aes = F, data = F_lab_week,
            aes(x = x, y = y, label = label), size = 1.75,
            color = F_blk, family = F_fnt_fam) +
  # single letters for study diet - first plot blank
  geom_label(inherit.aes = F, data = F_ltr_diet,
             aes(x = x, y = y, label = rev(label)), size = F_txt_size,
             label.padding = unit(0.1, "mm"), label.r = unit(0, "mm"),
             label.size = NA, vjust = "middle", hjust = "middle",
             fill = "white", color = NA, family = F_fnt_fam) +
  # single letters for study diet (drawn thrice to avoid using bold fontface)
  geom_text(inherit.aes = F, data = F_ltr_diet,
            aes(x = x, y = y, label = rev(label)), size = F_txt_size,
            color = F_hex_diet, family = F_fnt_fam) +
  geom_text(inherit.aes = F, data = F_ltr_diet,
            aes(x = x, y = y, label = rev(label)), size = F_txt_size,
            color = F_hex_diet, family = F_fnt_fam) +
  geom_text(inherit.aes = F, data = F_ltr_diet,
            aes(x = x, y = y, label = rev(label)), size = F_txt_size,
            color = F_hex_diet, family = F_fnt_fam) +
  # labels for study diet (drawn thrice to avoid using bold fontface)
  geom_text(inherit.aes = F, data = F_lab_diet,
            aes(x = x, y = y), size = 1.7,
            color = F_hex_diet, hjust = 0, family = F_fnt_fam,
            label = rev(c("[n=5] C: Control", "[n=4] R: Rice bran",
                          expression(paste("[n=4] F: ", italic("B. longum"), 
                                           "-fermented rice bran"))))) +
  geom_text(inherit.aes = F, data = F_lab_diet,
            aes(x = x, y = y), size = 1.7,
            color = F_hex_diet, hjust = 0, family = F_fnt_fam,
            label = rev(c("[n=5] C: Control", "[n=4] R: Rice bran",
                          expression(paste("[n=4] F: ", italic("B. longum"), 
                                           "-fermented rice bran"))))) +
  geom_text(inherit.aes = F, data = F_lab_diet,
            aes(x = x, y = y), size = 1.7,
            color = F_hex_diet, hjust = 0, family = F_fnt_fam,
            label = rev(c("[n=5] C: Control", "[n=4] R: Rice bran",
                          expression(paste("[n=4] F: ", italic("B. longum"), 
                                           "-fermented rice bran"))))) +
  # sample type collecton times
  geom_point(inherit.aes = F, data = F_shp_type,
             aes(x = x, y = y), size = F_shp_size, stroke = F_shp_stroke,
             shape = F_shp_type$shape, color = F_shp_type$color) +
  # labels for shapes (drawn thrice to avoid using bold fontface)
  geom_text(inherit.aes = F, data = F_lab_type,
            aes(x = x, y = y, label = label), size = F_txt_size,
            color = F_lab_type$color, hjust = 0, family = F_fnt_fam) +
  geom_text(inherit.aes = F, data = F_lab_type,
            aes(x = x, y = y, label = label), size = F_txt_size,
            color = F_lab_type$color, hjust = 0, family = F_fnt_fam) +
  geom_text(inherit.aes = F, data = F_lab_type,
            aes(x = x, y = y, label = label), size = F_txt_size,
            color = F_lab_type$color, hjust = 0, family = F_fnt_fam) +
  # arrows for diet intervention and end of study
  geom_segment(inherit.aes = F, data = F_arw,
               aes(x = xb, xend = xe, y = yb, yend = ye), size = F_inn_box_size,
               color = F_blk, linejoin = "mitre", arrow = F_arw_typ) +
  # labels for diet intervention and end of study arrows
  geom_text(inherit.aes = F, data = F_lab_arw,
            aes(x = x, y = y, label = label), size = 1.75,
            color = F_blk, family = F_fnt_fam) +
  # mouse info
  geom_text(inherit.aes = F, data = F_mouse,
            aes(x = x, y = y, label = label), size = F_txt_size,
            color = F_blk, hjust = 0, family = F_fnt_fam) +
  scale_x_continuous(limits = c(F_x_lim), expand = c(0, 0.2)) +
  scale_y_continuous(limits = c(F_y_lim), expand = c(0, 0.2)) +
  theme_void()

### ************************************
### SECTION F - WRITE OUTPUTS ----
### ************************************

# relative paths from wd for outputs headed for storage in "main/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = 125, height = 37.5,
       filename = ofm_study_timeline, plot = F_ggp_timeline)

# relative paths for outputs headed for storage in the "vault/"
ggsave(device = "pdf", dpi = 600, units = "mm", width = 125, height = 37.5,
       filename = F_ofv_timeline, plot = F_ggp_timeline)

# create a list of objects from section and save them as a workspace file
F_obj <- c(ls(pattern = "F_"))
save(list = F_obj, file = F_ofv_WS)

### ************************************
### EPILOGUE ----
### ************************************

# ... Epilogue... Nineteen Years Later ...
# yes it is part of the book's canon and no one can debate that
# but it does not have to exist as part of my canon as a reader
# it was the single most unsatisfying chapter I have ever read
# it did not add anything of meaning to the stories that I hold so dear

# anyway, the final workspace is saved here, it may take a little while
save.image(file = ofv_WS_all)

# if for some reason something fails or something is missing or you have an idea
# or you just want to chat about how terrible Nineteen Years Later was...
# open an issue on the project's page at: github.com/kdprkr/MerlinsManuscript/
# or contact me via email at: kristopher.parker@colostate.edu

# oh and don't worry about those pesky warning messages that read:
# In is.na(x) : is.na() applied to non-(list or vector) of type 'expression'
# they come from the italic text expressions in some of the ggplots
# ... do be worried if there are other warnings, however

# All was well [without that epilogue, Ms. Rowling].