#!/usr/bin/env Rscript

# This is a script creating plots showing the distribution of some QC values, 
# and where the current normal and tumor sample lies within those distributions.
# QC metrics shown are read count, duplication rate, coverage, fold enrichment, on-bait rate, insert size & contamination


library(ggplot2)
library(optparse)
library(data.table)

# read in command line options
option_list <- list(
  make_option(c("-s", "--samples"), action = "store", type = "character",
              default = NULL, help = "Samples of interest to match with qc file names. If more than one they must be separated with colon, e.g. tumor:normal"),
  make_option(c("-o", "--outfile"), action = "store", type = "character",
              default = "QC_overview.pdf", help = "Path to output pdf file"),
  make_option(c("-m", "--mainpath"), action = "store", type = "character",
              default = "/nfs/PROBIO/autoseq-output", help = "Path to output pdf file")
)

opt <- parse_args(OptionParser(option_list = option_list))
sample_string = opt$samples
outfile = opt$outfile
main_path = opt$mainpath


# find the qc files for all samples
cat("Find all available qc files...\n")
hsmetrics_files = dir(Sys.glob(paste0(main_path, "/*/*/qc/picard")), 
                pattern = "picard-hsmetrics.txt$", recursive = TRUE, full.names = TRUE)
markduplicates_files = dir(Sys.glob(paste0(main_path, "/*/*/qc/picard")), 
                           pattern = "markdups-metrics.txt$", recursive = TRUE, full.names = TRUE)
insertsize_files = dir(Sys.glob(paste0(main_path, "/*/*/qc/picard")), 
                           pattern = "picard-insertsize.txt$", recursive = TRUE, full.names = TRUE)
contest_files = dir(Sys.glob(paste0(main_path, "/*/*/contamination")), 
                    pattern = "contest.txt$", recursive = TRUE, full.names = TRUE)

# remove files from old pipeline (modified before 2019-04-30)
hsmetrics_files = hsmetrics_files[file.mtime(hsmetrics_files) > as.POSIXct("2019-04-30")]
markduplicates_files = markduplicates_files[file.mtime(markduplicates_files) > as.POSIXct("2019-04-30")]
insertsize_files = insertsize_files[file.mtime(insertsize_files) > as.POSIXct("2019-04-30")]
contest_files = contest_files[file.mtime(contest_files) > as.POSIXct("2019-04-30")]


# read in the files
cat("Read in the files...\n")
HsMetrics = data.frame()
for (f in hsmetrics_files) {
  SAMP = strsplit(basename(f), split = "\\.")[[1]][1]
  HsMetrics = rbind(HsMetrics, cbind(SAMP, read.table(f, skip = 6, nrow = 1, sep = "\t", 
                                                        header = TRUE, stringsAsFactors = FALSE), 
                                     stringsAsFactors = FALSE))
}
MarkDuplicates = data.frame()
for (f in markduplicates_files) {
  SAMP = sub("-markdups-metrics.txt", "", basename(f))
  MarkDuplicates = rbind(MarkDuplicates, cbind(SAMP, read.table(f, skip = 6, nrow = 1, sep = "\t", 
                                                                  header = TRUE, stringsAsFactors = FALSE), 
                                               stringsAsFactors = FALSE))
}
InsertSize = data.frame()
for (f in insertsize_files) {
  SAMP = strsplit(basename(f), split = "\\.")[[1]][1]
  InsertSize = rbind(InsertSize, cbind(SAMP, read.table(f, skip = 6, nrow = 1, sep = "\t", 
                                                        header = TRUE, stringsAsFactors = FALSE), 
                                     stringsAsFactors = FALSE))
}
ContEst = data.frame()
for (f in contest_files) {
  SAMP = strsplit(basename(f), split = "\\.")[[1]][1]
  ContEst = rbind(ContEst, cbind(SAMP, read.table(f, nrow = 1, sep = "\t", header = TRUE, stringsAsFactors = FALSE), 
                                       stringsAsFactors = FALSE))
}

# merge the QC tables and add sample type and capture kit
qc_merge = merge(merge(merge(HsMetrics, MarkDuplicates, by = "SAMP"), InsertSize, by = "SAMP"), ContEst, by = "SAMP")
qc_merge$sample_type = sapply(strsplit(qc_merge$SAMP, split = "-"), "[", 4)
qc_merge$capture = sapply(strsplit(qc_merge$SAMP, split = "-"), "[", 7)

# label the samples of interest (soi)
samples = strsplit(sample_string, split = ":")[[1]]
qc_merge$soi = qc_merge$SAMP %in% samples


# create an ouput table for the samples of interest
soi_table = data.table(qc_merge)[i = soi==TRUE, 
                                 j =list(SAMP, MEAN_TARGET_COVERAGE, FOLD_ENRICHMENT, dedupped_on_bait_rate=ON_BAIT_BASES/PF_BASES_ALIGNED, 
                                         READ_PAIRS_EXAMINED, PERCENT_DUPLICATION, "contamination_%"=contamination, MEDIAN_INSERT_SIZE)]
table_outfile = sub("pdf$", "txt", outfile)
write.table(x = soi_table, file = table_outfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




##################################################################################
#   ____   __       ____   _______
#  | o  \ | |     /     \ |__  __|
#  | __/  | |    |  / \ |   | |
#  ||     | |___ | \ /  |   | |
#  ||     |____| \____ /    |_|
#################################################################################

# Center plot titles (must be specified as of ggplot v 2.2.0)
theme_update(plot.title = element_text(hjust = 0.5))
theme_update(plot.subtitle = element_text(hjust = 0.5))

# plotting function to create desired histograms 
my_histogram = function(x, binwidth, sample_binwidth, ybreaks, xbreaks, xbreaks_minor, x_string, title_string) {
  p = ggplot(qc_merge, aes_string(x = x)) +
    geom_histogram(aes(group = capture, fill = capture), binwidth = binwidth, alpha = 0.7, color = "black") +
    # geom_freqpoly(binwidth = binwidth, alpha = 0.5) +
    geom_histogram(data = subset(qc_merge, soi), binwidth = sample_binwidth, fill = "red") + # the sample of interest
    scale_fill_manual(values = c("antiquewhite", "aliceblue")) +
    scale_y_continuous(breaks = ybreaks) +
    scale_x_continuous(name = x_string, breaks = xbreaks, minor_breaks = xbreaks_minor) +
    facet_wrap(~sample_type, ncol = 1) +
    ggtitle(paste0(title_string, " (binwidth ", binwidth, ")"), 
            subtitle = paste0("Red piles show present samples (binwidth ", sample_binwidth, ")"))
  print(p)
}

# plotting function to create desired scatter plots
my_scatter = function(x, y, xbreaks, ybreaks, x_string, y_string, title_string) {
  ggplot(qc_merge, aes(shape = capture)) +
    geom_point(aes_string(x = x, y = y), size = 3) +
    geom_point(data = subset(qc_merge, soi), aes_string(x = x, y = y), fill = "red", size = 3) +
    scale_alpha_manual(values = c(0.7, 1)) +
    # scale_color_manual(values = c("antiquewhite", "aliceblue")) +
    # scale_fill_manual(values = c("antiquewhite", "aliceblue")) +
    scale_shape_manual(values = c(24, 25), guide = guide_legend(override.aes = list(fill = NA))) +
    # scale_size_manual(values = c(3, 7), labels = NULL, guide = FALSE) +
    scale_x_continuous(name = x_string, breaks = xbreaks) +
    scale_y_continuous(name = y_string, breaks = ybreaks) +
    facet_wrap(~sample_type, ncol = 1) +
    ggtitle(title_string, subtitle = "Red symbols show present samples")
}

# Plot histograms and scatter plots
cat(paste0("Create plots (saved in ", outfile, ") ...\n"))
pdf(file = outfile, width=14)

# read count histogram
my_histogram(x = "READ_PAIRS_EXAMINED", binwidth = 5e6, sample_binwidth = 1e5, ybreaks = seq(0,nrow(qc_merge), 1),
             xbreaks = seq(0, 1e12, 1e7), xbreaks_minor = seq(0, 1e12, 5e6), x_string = "number of read pairs", 
             title_string = "Read count")

# duplication histogram
my_histogram(x = "PERCENT_DUPLICATION", binwidth = 0.05, sample_binwidth = 0.01, ybreaks = seq(0,nrow(qc_merge), 1),
             xbreaks = seq(0, 1, 0.1), xbreaks_minor = seq(0, 1, 0.05), x_string = "duplication rate", 
             title_string = "Duplication rate")

# duplication vs read count scatter plot
my_scatter(x = "READ_PAIRS_EXAMINED", y = "PERCENT_DUPLICATION", xbreaks = seq(0, 1e12, 1e7), ybreaks = seq(0, 1, 0.1),
           x_string = "number of read pairs", y_string = "duplication rate", title_string = "Duplication rate vs Read count")

# coverage histogram
my_histogram(x = "MEAN_TARGET_COVERAGE", binwidth = 100, sample_binwidth = 1, ybreaks = seq(0,nrow(qc_merge), 1),
             xbreaks = seq(0, 5000, 500), xbreaks_minor = seq(0, 5000, 100), x_string = "mean target coverage", 
             title_string = "Coverage")

# coverage vs read count scatter plot
my_scatter(x = "READ_PAIRS_EXAMINED", y = "MEAN_TARGET_COVERAGE", xbreaks = seq(0, 1e12, 1e7), ybreaks = seq(0, 5000, 500),
           x_string = "number of read pairs", y_string = "mean target coverage", title_string = "Coverage rate vs Read count")

# fold enrichment histogram
my_histogram(x = "FOLD_ENRICHMENT", binwidth = 50, sample_binwidth = 1, ybreaks = seq(0,nrow(qc_merge), 1),
             xbreaks = seq(0, 5000, 100), xbreaks_minor = seq(0, 5000, 50), x_string = "fold enrichment", 
             title_string = "Fold enrichment")

# fold enrichment vs read count scatter plot
my_scatter(x = "READ_PAIRS_EXAMINED", y = "FOLD_ENRICHMENT", xbreaks = seq(0, 1e12, 1e7), ybreaks = seq(0, 5000, 100),
           x_string = "number of read pairs", y_string = "fold enrichment", title_string = "Fold enrichment rate vs Read count")

# on-bait rate histogram
my_histogram(x = "ON_BAIT_BASES/PF_BASES_ALIGNED", binwidth = 0.05, sample_binwidth = 0.01, ybreaks = seq(0,nrow(qc_merge), 1),
             xbreaks = seq(0, 1, 0.05), xbreaks_minor = seq(0, 1, 0.01), x_string = "dedupped on-bait rate", 
             title_string = "Dedupped on-bait rate")

# on-bait rate vs read count scatter plot
my_scatter(x = "READ_PAIRS_EXAMINED", y = "ON_BAIT_BASES/PF_BASES_ALIGNED", xbreaks = seq(0, 1e12, 1e7), ybreaks = seq(0, 1, 0.05),
           x_string = "number of read pairs", y_string = "dedupped on-bait rate", title_string = "Dedupped on-bait rate rate vs Read count")

# insert size histogram
my_histogram(x = "MEDIAN_INSERT_SIZE", binwidth = 2, sample_binwidth = 1, ybreaks = seq(0,nrow(qc_merge), 1),
             xbreaks = seq(0, 1000, 10), xbreaks_minor = seq(0, 1000, 2), x_string = "median insert size", 
             title_string = "Insert size")

# contamination histogram
my_histogram(x = "contamination", binwidth = 0.1, sample_binwidth = 0.1, ybreaks = seq(0,nrow(qc_merge), 2),
             xbreaks = seq(0, 100, 0.5), xbreaks_minor = seq(0, 100, 0.1), x_string = "contamination, %", 
             title_string = "Contamination")


dev.off()

cat("Completed.\n")


