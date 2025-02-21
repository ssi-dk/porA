#!/usr/bin/env Rscript

# Load libraries
suppressWarnings(suppressPackageStartupMessages(library(future.apply)))
suppressWarnings(suppressPackageStartupMessages(library(optparse)))
suppressWarnings(suppressPackageStartupMessages(library(dada2)))
suppressWarnings(suppressPackageStartupMessages(library(phyloseq)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(ape)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))
suppressWarnings(suppressPackageStartupMessages(library(stringr)))
suppressWarnings(suppressPackageStartupMessages(library(DECIPHER)))
suppressWarnings(suppressPackageStartupMessages(library(phangorn)))
suppressWarnings(suppressPackageStartupMessages(library(ggtree)))
suppressWarnings(suppressPackageStartupMessages(library(igraph)))
suppressWarnings(suppressPackageStartupMessages(library(ggraph)))
suppressWarnings(suppressPackageStartupMessages(library(dendextend)))
suppressWarnings(suppressPackageStartupMessages(library(forcats))) # For fct_reorder
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))
suppressWarnings(suppressPackageStartupMessages(library(scales)))
suppressWarnings(suppressPackageStartupMessages(library(seqinr)))
suppressWarnings(suppressPackageStartupMessages(library(msa)))


# Create an OptionParser object for running automatic with Rscript 
option_list = list(
  make_option("--seqtab", type="character", default=NULL, 
              help="Comma seperated list of saved seqtab.rds paths", metavar="character"),
  make_option("--out", type="character", default=NULL, 
              help="Output folder name, could be runid to add to the same as the other plots, or a combination of two runids", metavar="character"),
  make_option("--barplot_height", type="integer", default=7, metavar = "interger", 
              help="Output barplot height, default 7"),
  make_option("--barplot_width", type="integer", default=NULL, metavar = "interger",
              help="Output barplot width, default number of samples times 0.5"),
  make_option("--mds_height", type="integer", default=8, metavar = "interger", 
              help="Output plot height, default 8"),
  make_option("--mds_width", type="integer", default=8, metavar = "interger",
              help="Output plot width, default 8"),
  make_option("--tree_height", type="integer", default=NULL, metavar = "interger", 
              help="Output tree plot height, default number of types found times 0.1"),
  make_option("--tree_width", type="integer", default=10, metavar = "interger",
              help="Output tree plot width, default 10"),
  make_option("--prefix", type="character", default=NULL, 
              help="Path to result folder", metavar="character"),
  make_option("--threads", type="integer", default=6, metavar = "interger", 
              help="Number of threads, default 6")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#make sure all arguments are there as input 
if (is.null(opt$seqtab) || is.null(opt$prefix) || is.null(opt$out)) {
  stop("Not all input required have been entered. Run with comma seperated list of seqtab.rds to include, runid or other name for the folder to save the results in and prefix, the base folder where all results should be saved")
}

#check if output folders already exist otherwise create them 
if (!dir.exists(paste(opt$prefix, "tables", opt$out, sep="/"))) {
  dir.create(paste(opt$prefix, "tables", opt$out, sep="/"), recursive = TRUE)
}

if (!dir.exists(paste(opt$prefix, "plots", opt$out, sep="/"))) {
  dir.create(paste(opt$prefix, "plots", opt$out, sep="/"), recursive = TRUE)
}

# Set up parallel for future apply 
plan(multisession, workers = opt$threads)


#read seqtabs, split if multiple 
if(str_count(opt$seqtab, ',')>=1){
  seqtabs_path<-unlist(strsplit(opt$seqtab, ","))
  seqtabs<-future_lapply(seqtabs_path, readRDS)
  seqtab_merge <- mergeSequenceTables(tables=seqtabs)
}  else if(str_count(opt$seqtab, ',')==0){
  seqtab_merge=readRDS(opt$seqtab) 
} else{
  stop("Enter at least on seqtab.rds file")
}


##remove Flok if in samplename 
rownames(seqtab_merge)<-gsub('Flok', '', rownames(seqtab_merge))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab_merge, method="pooled", multithread=TRUE, verbose = TRUE)

#write ASV table
asv <- seqtab.nochim %>% t() %>% data.frame(check.names=FALSE) %>% rownames_to_column("ASV")
write.table(asv, file=paste(opt$prefix, "tables", opt$out, "ASV.tsv", sep="/"), sep="\t", quote=FALSE, col.names = NA, row.names = TRUE)
#write.table(asv, file="Scratch/ASV.tsv", sep="\t", quote=FALSE)

##phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE))

##move sequence to refseq and rename to number 
dna <- Biostrings::DNAStringSet(taxa_names(ps))##get DNA seq
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- seq(ntaxa(ps))

##add types 
types<-matrix(names(refseq(ps)))
rownames(types)<- names(refseq(ps))
colnames(types)<- "Type"

##add samples
sampleinfo<-data.frame('Sampleid'=rownames(otu_table(ps)), extra=rep(1, length(rownames(otu_table(ps)))))##minimum 2 cols for it to work 
rownames(sampleinfo)<- rownames(otu_table(ps))

##merge all into one obj
ps<-merge_phyloseq(ps, tax_table(types), sample_data(sampleinfo))

#save amplicons to fasta file 
writeXStringSet(refseq(ps), filepath = paste(opt$prefix, "tables", opt$out, "amplicons.fa", sep="/"))  

##clean otu - add all vars only 1nt from the most abundant 
##align amplicons
system(paste("mafft --thread", opt$threads, "--auto",
             paste(opt$prefix, "tables", opt$out, "amplicons.fa", sep="/"), 
             ">", 
             paste(opt$prefix, "tables", opt$out, "amplicons_mafft.fa", sep="/")))
ali <- readDNAStringSet(paste(opt$prefix, "tables", opt$out, "amplicons_mafft.fa", sep="/"))


##calcualte pairwise distances 
pd<-dist.dna(as.DNAbin(ali), model="N", as.matrix=T, pairwise.deletion=TRUE)
write.table(pd, file=paste(opt$prefix, "tables", opt$out, "ASV_dist.tsv", sep="/"), sep="\t", quote=FALSE, col.names = NA, row.names = TRUE)

# Ensure OTU data is valid
otu_data <- as.data.frame(otu_table(ps))  # Extract OTU table from phyloseq object
if (nrow(otu_data) == 0 || ncol(otu_data) == 0) {
  stop("OTU table is empty or invalid. Check your phyloseq object.")
}

# Precompute similar ASVs
similar_asvs <- future_lapply(rownames(pd), function(otu) {
  names(which(pd[otu, ] <= 1 & pd[otu, ] > 0))  # OTUs within 1 nucleotide difference
})
names(similar_asvs) <- rownames(pd)

# function to process each sample - add all vars with 1 diff to the most abundant var
process_sample <- function(sample, otu_data, similar_asvs) {
  # Ensure sample_reads is a named numeric vector
  sample_reads <- as.numeric(otu_data[sample, ])
  names(sample_reads) <- colnames(otu_data)

  while (any(sample_reads > 0)) {
    # Identify the most abundant variable
    max_var <- names(sample_reads)[which.max(sample_reads)]

    # Get all similar variables with non-zero reads
    similar_vars <- similar_asvs[[max_var]]
    similar_vars <- similar_vars[sample_reads[similar_vars] > 0]  # Filter non-zero vars

    # Merge reads from similar variables into max_var
    for (var in similar_vars) {
      if (var != max_var) {
        sample_reads[max_var] <- sample_reads[max_var] + sample_reads[var]
        sample_reads[var] <- 0  # Zero out the similar variable
      }
    }
    
    # Mark max_var as processed (but retain its final count)
    sample_reads[max_var] <- -sample_reads[max_var]  # Temporarily mark as processed
  }
  
  # Convert negative values (processed variables) back to positive
  sample_reads <- abs(sample_reads)
  return(sample_reads)
}


# Apply the function to all samples in parallel
test_result <- future_lapply(rownames(otu_data), function(sample) {
  process_sample(sample, otu_data, similar_asvs)
})

# Combine the results into a single data frame
updated_abundance_data <- do.call(rbind, test_result)
rownames(updated_abundance_data) <- rownames(otu_data)
colnames(updated_abundance_data) <- colnames(otu_data)


# Update the phyloseq object
otu_table(ps) <- otu_table(as.matrix(updated_abundance_data), taxa_are_rows = FALSE)
write.table(otu_table(ps), file=paste(opt$prefix, "tables", opt$out, "OTU_table_processed.tsv", sep="/"), sep="\t", quote=FALSE, col.names = NA, row.names = TRUE)

##only keep values above 20 for plots sake 
otu_table(ps) <- otu_table(as.matrix(updated_abundance_data) * (as.matrix(updated_abundance_data) >= 20), taxa_are_rows = FALSE)

##split by if it starts with a digit or what letter it starts with, thereby splitting by flok 
samples_d <- grep("^[A-Za-z]{1}-", rownames(otu_table(ps)), value = TRUE, invert = TRUE)##would be all starting with digits (our first mock and other samples with multiple letters, not just flok)
samples_d<-sort(samples_d)
samples_l <- grep("^[A-Za-z]{1}-", rownames(otu_table(ps)), value = TRUE) ##just our floks, example N-S1 
samples_l<-sort(samples_l)
unique_letters <- unique(toupper(substring(samples_l, 1, 1)))

sorted_samples <- c(samples_d, samples_l)

# Reorder sample data in the phyloseq object and split by flok
sample_data(ps)$Sampleid <- factor(rownames(sample_data(ps)), levels = sorted_samples)

split_list <- list()
if (length(samples_d)>0){
  split_list[["mock_other"]] <- samples_d
}

if (length(samples_l)>0){
  for (letter in unique_letters) {
    split_list[[letter]] <- grep(paste0("^", letter, "-"), rownames(otu_table(ps)), value = TRUE, ignore.case = TRUE)
  }
}

phyloseq_list <- future_lapply(split_list, function(samples) prune_samples(samples, ps))


##colors for barplots 
full_df<-psmelt(ps)

type_abundance <- aggregate(Abundance ~ Type, full_df, sum)

type_abundance$Group <- with(type_abundance, case_when(
  Abundance >= 10000 ~ "High",
  Abundance > 100 & Abundance < 10000 ~ "Mid",
  Abundance <= 100 ~ "Low"
))

# Split Types into groups
high_types <- type_abundance$Type[type_abundance$Group == "High"]
mid_types <- type_abundance$Type[type_abundance$Group == "Mid"]
low_types <- type_abundance$Type[type_abundance$Group == "Low"]

# Assign different color palettes for each group
high_colors <- colorRampPalette(brewer.pal(n = 9, name = "Set1"))(length(high_types))
mid_colors  <- colorRampPalette(brewer.pal(n = 9, name = "Paired"))(length(mid_types))
low_colors  <- colorRampPalette(brewer.pal(n = 9, name = "Set3"))(length(low_types)) 

# Create a named vector mapping Type -> Color
color_mapping <- setNames(c(high_colors, mid_colors, low_colors),
                          c(high_types, mid_types, low_types))


##make plots for each group 
# Parallelize the generation of barplots
future_lapply(names(phyloseq_list), function(name) {
  sub_ps <- phyloseq_list[[name]]

  df <- psmelt(sub_ps)  # Convert phyloseq object to a data frame
  

  # Reorder by abundance
  df$Type <- fct_reorder(df$Type, df$Abundance)  
  
  
  # Generate the barplot for counts
  if (is.null(opt$barplot_width)){
    opt$barplot_width<-nrow(sample_data(sub_ps))*0.5
  }


  pdf(paste(opt$prefix, "/plots/", opt$out, "/", name, "_type_counts_barplot.pdf", sep=""), height=opt$barplot_height, width=opt$barplot_width)
    p <- ggplot(df, aes(x=Sample, y=Abundance, fill=Type)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=color_mapping) +
    geom_text(aes(label=ifelse(Abundance > 200, as.character(Type), "")),
              vjust=0.5,
              position=position_stack(vjust=0.5)) +
    labs(y="Reads") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle=90, hjust=1))  # Rotate sample names
  
  print(p)
  dev.off()
  
  
    ##proportion
    sub_ps <- transform_sample_counts(sub_ps, function(OTU) OTU/sum(OTU))  # Normalize to %

    df_prop <- psmelt(sub_ps)  # Convert to dataframe
    df_prop$Type <- fct_reorder(df_prop$Type, df_prop$Abundance)  # Reorder by abundance
    
    pdf(paste(opt$prefix, "/plots/", opt$out, "/", name, "_type_proportion_barplot.pdf", sep=""), height=opt$barplot_height, width=opt$barplot_width)
      p <- ggplot(df_prop, aes(x=Sample, y=Abundance, fill=Type)) +
      geom_bar(stat="identity") +
      scale_fill_manual(values=color_mapping) +
      geom_text(aes(label=ifelse(Abundance > 0.01, as.character(Type), "")),
                vjust=0.5,
                position=position_stack(vjust=0.5)) +
      labs(y="Abundance (%)") +
      theme(legend.position = "none",
            axis.text.x = element_text(angle=90, hjust=1))  # Rotate sample names

    print(p)
    dev.off()
})



##per sample type dendrograms

# Convert `pd` to a full numeric matrix
pd_mat <- if (inherits(pd, "dist")) as.matrix(pd) else pd

# Function to plot type dendrogram for one or all samples
plot_dendrogram_for_sample <- function(ps, pd_mat, sample_name = NULL) {

  otu_counts <- as(otu_table(ps), "matrix")  # Extract OTU table

  # Ensure we only process the selected sample
  if (!is.null(sample_name)) {
    if (!(sample_name %in% rownames(otu_counts))) {
      stop("Sample not found in dataset.")
    }
    sample_list <- sample_name
  } else {
    sample_list <- rownames(otu_counts)
  }

  for (sample in sample_list) {

    # Get OTUs with non-zero counts
    sample_otus <- which(otu_counts[sample, ] > 0)

    if (length(sample_otus) == 0) {
      cat("No OTUs with counts > 0 in sample:", sample, "\n")
      next
    }

    if (length(sample_otus) < 2) {
      cat("Skipping sample:", sample, "- Not enough OTUs for clustering\n")
      next
    }

    # Subset `pd_mat` to only the OTUs present in the sample
    pd_sample <- pd_mat[sample_otus, sample_otus]

    # Convert `pd_sample` into a phylogenetic tree (dendrogram format)
    hc <- hclust(as.dist(pd_sample), method="average")  # Hierarchical clustering
    tree <- phyloseq::as.phylo(hc)  # Convert to phylo format

    # Scale edge lengths (normalize to number of segregating sites)
    if ("seg.sites" %in% ls()) {
      tree$edge.length <- tree$edge.length / sum(tree$edge.length) * length(seg.sites(ali))
    }

    # Extract OTU types if available
    if (!is.null(tax_table(ps, errorIfNULL=FALSE))) {
      otu_types <- data.frame(label = taxa_names(ps), Type = as.character(tax_table(ps)[, 1]))  # Assuming first column is the type
    } else {
      otu_types <- data.frame(label = tree$tip.label, Type = rep("Unknown", length(tree$tip.label)))  # Default to "Unknown"
    }

    # Create dataframe for node sizes based on abundance **for this sample**
    otu_sizes <- otu_counts[sample, sample_otus]  # Get only for this sample
    otu_sizes_df <- data.frame(label = names(otu_sizes), abundance = otu_sizes)

    # Combine data for visualization
    tip_data <- merge(otu_types, otu_sizes_df, by="label", all.x=TRUE)

    # Create edge dataframe for branch labels
    edges <- data.frame(tree$edge, edgeLab = tree$edge.length)
    colnames(edges) <- c("parent", "node", "edgeLab")
    edges$edgeLab <- round(edges$edgeLab, 2)  # Round distances
    edges$edgeLab[edges$edgeLab < 0.01] <- NA  # Remove very small distances

    # Plot dendrogram with ggtree
    p <- ggtree(tree) +
      theme_tree2() +  # Ensure proper alignment
      ggtitle(paste("Sample", sample))

    # Add tip labels with colors (e.g., types)
    p <- p %<+% tip_data +
      geom_tiplab(aes(label=Type), hjust=-0.00, vjust=-0.5, size=4, color="black") +  # Type labels
      geom_tippoint(aes(size=abundance), color="plum4") +  # Node size by abundance
      scale_size(range = c(2, 6))  # Scale node sizes

    # Add branch labels (distance values)
    p <- p %<+% edges +
      geom_text(aes(x=branch, label=edgeLab), size=3, color="gray40", vjust=-0.5, na.rm=TRUE)

    # Show plot
    pdf(paste(opt$prefix, "/plots/", opt$out, "/", sample, "_type_dendrogram.pdf", sep=""))
    print(p)
    dev.off()
  }
}


# Run for a single sample:
#plot_dendrogram_for_sample(ps, pd_mat, "MC-3_seq92")
# Run for all samples:
plot_dendrogram_for_sample(ps, pd_mat)


