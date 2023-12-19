# mane_seq_logo.R

library(data.table)
library(ggplot2)
library(ggseqlogo)

# Read in MANE Select transcript data once transformed in to a data.table
# i.e. tsv file contain |transcript_id|cnda_seq|five_prime_utr_length|
b <- fread("...")

# Create the 11-bp context
b[, context:= toupper(substr(seq, five_prime_utr_length-5, five_prime_utr_length+5))]

# Filter to context 
b %<>% .[nchar(context) == 11]

# Remove ATG for plotting 
b[, non_atg_context := paste0(substr(context, 1,6), "---", substr(context, 10, 11))]

# Create a matrix
seqs <- b$non_atg_context
seq_matrix <- matrix(unlist(strsplit(seqs, "")), ncol = 11, byrow=T)
modified_seqs <- apply(seq_matrix, 1, paste, collapse = "")

# Plot using ggseqlogo
ggplot() + geom_logo(modified_seqs, method="bits")+theme_logo()