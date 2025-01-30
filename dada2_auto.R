#!/usr/bin/env Rscript

# Load libraries:
suppressWarnings(suppressPackageStartupMessages(library(optparse)))
suppressWarnings(suppressPackageStartupMessages(library(dada2)))
suppressWarnings(suppressPackageStartupMessages(library(phyloseq)))
suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(tibble)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(ape)))



# Create an OptionParser object - run directly from terminal 
option_list = list(
  make_option("--r1", type="character", default=NULL, 
              help="Comma seperated list of forward read(1) paths", metavar="character"),
  make_option("--r2", type="character", default=NULL, 
              help="Comma seperated list of reverse read(2) paths", metavar="character"),
  make_option("--prefix", type="character", default=NULL, 
              help="Path to result folder", metavar="character"),
  make_option("--trunc1", type="integer", default=280, 
              help="Truncation length for read 1 (fwd)", metavar="integer"),
  make_option("--trunc2", type="integer", default=260, 
              help="Truncation length for read 2 (rev)", metavar="integer"),
  make_option("--run", type="character", default=NULL, 
              help="Name of runid", metavar="character")

); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


##get reads
str1<-opt$r1
r1<-unlist(strsplit(str1, ","))

str2<-opt$r2
r2<-unlist(strsplit(str2, ","))

##height and width for plot
len<-length(r1)
h<-round((2.5*len))
w<-round((4*len))

##raw read1 quality plot 
pdf(paste(opt$prefix, "plots", opt$run, "Read1_quality.pdf", sep="/"), height=h, width=w)
plotQualityProfile(r1) + theme(text = element_text(size = 10))
dev.off()

##one plot per sample
# for (i in r1){
#   fname<-basename(i)
#   subname<-sub("_rmprimer.*", "", fname)
#   pdf(paste(opt$prefix, "plots", opt$run, sprintf("%s_read1_quality.pdf", subname), sep="/"), height=2.5, width=4)
#   print(plotQualityProfile(i) + theme(text = element_text(size = 3)))
#   dev.off()
# }

##raw read2 quality plot 
pdf(paste(opt$prefix, "plots", opt$run, "Read2_quality.pdf", sep="/"), height=h, width=w)
plotQualityProfile(r2) + theme(text = element_text(size = 10))
dev.off()

##read qc filtering 
filts <- file.path(opt$prefix, "data/filtered", opt$run)

# Filter and trim
out <- filterAndTrim(fwd=r1, filt=filts, rev=r2, filt.rev=filts,
                     maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE, truncLen=c(opt$trunc1,opt$trunc2),
                     compress=TRUE, multithread=TRUE, matchIDs = TRUE) # On Windows set multithread=FALSE

row.names(out)<-sapply(strsplit(basename(row.names(out)), "_rmprimer"), `[`, 1)

#sort by rowname
out<-out[ order(row.names(out)), ]

filtFs<-list.files(filts, pattern='_R1', full.names=TRUE)
filtRs<-list.files(filts, pattern='_R2', full.names=TRUE)
names(filtFs) <- sapply(strsplit(basename(filtFs), "_rmprimer"), `[`, 1)
names(filtRs) <- sapply(strsplit(basename(filtRs), "_rmprimer"), `[`, 1)

##filtered read1 quality plot
pdf(paste(opt$prefix, "plots", opt$run, "filter_read1_quality.pdf", sep="/"), height=h, width=w)
plotQualityProfile(filtFs) + theme(text = element_text(size = 10))
dev.off()

##filtered read2 quality plot
pdf(paste(opt$prefix, "plots", opt$run, "filter_read2_quality.pdf", sep="/"), height=h, width=w)
plotQualityProfile(filtRs) + theme(text = element_text(size = 10))
dev.off()

##dereplicate - not necessary in new version - done on the go with learnError and dada
#derepF <- derepFastq(filtFs, verbose=TRUE)
#derepR <- derepFastq(filtRs, verbose=TRUE)


# Learn the error rates - parametric model of the errors introduced by PCR amplification and sequencing.
set.seed(100)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

##error rate plot read1
pdf(paste(opt$prefix, "plots", opt$run, "Read1_errorRate.pdf", sep="/"))
plotErrors(errF, nominalQ=TRUE)
dev.off()

##error rate plot read2
pdf(paste(opt$prefix, "plots", opt$run, "Read2_errorRate.pdf", sep="/"))
plotErrors(errR, nominalQ=TRUE)
dev.off()

# Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool=TRUE)


# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)

getN <- function(x) sum(getUniques(x))

# Ensure row names are consistent and unique for all datasets
sample_names <- unique(c(rownames(out), names(dadaFs), names(dadaRs), names(mergers), rownames(seqtab)))

# Create a master data frame with all sample names
track <- data.frame(Sample = sample_names)

# Add columns for each step, filling missing samples with NA
track$out <- out[match(sample_names, rownames(out)), ]
track$dadaFs <- sapply(sample_names, function(s) ifelse(s %in% names(dadaFs), getN(dadaFs[[s]]), NA))
track$dadaRs <- sapply(sample_names, function(s) ifelse(s %in% names(dadaRs), getN(dadaRs[[s]]), NA))
track$mergers <- sapply(sample_names, function(s) ifelse(s %in% names(mergers), getN(mergers[[s]]), NA))
track$seqtab <- sapply(sample_names, function(s) ifelse(s %in% rownames(seqtab), rowSums(seqtab[s, , drop = FALSE]), NA))

write.table(track, file=paste(opt$prefix, "tables", opt$run, "reads_overview.tsv", sep="/"), sep="\t", quote=FALSE, row.names = FALSE)

rownames(seqtab)<-gsub("_rmprimer.*", "", rownames(seqtab), perl=TRUE)
saveRDS(seqtab, paste(opt$prefix, "data/seqtabs", opt$run, "seqtab.rds", sep="/"))


