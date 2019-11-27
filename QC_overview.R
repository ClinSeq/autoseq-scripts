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
  make_option(c("-d", "--analysisdir"), action = "store", type = "character",
              default = NULL, help = "Analysis directory for samples of interest, in order to separate between the same samples that possibly are being run in multiple analyses"),
  make_option(c("-o", "--outfile"), action = "store", type = "character",
              default = "QC_overview.pdf", help = "Path to output pdf file"),
  make_option(c("-m", "--mainpath"), action = "store", type = "character",
              default = "/nfs/PROBIO/autoseq-output", help = "Path to output pdf file")
)

opt <- parse_args(OptionParser(option_list = option_list))
sample_string = opt$samples
analysis_dir = opt$analysisdir
outfile = opt$outfile
main_path = opt$mainpath

# for testing:
# system("umount ~/anchorageMountpointR2")
# system("sshfs rebber@anchorage:/nfs/CLINSEQ/PSFF/autoseq-output/ ~/anchorageMountpointR2")
# sample_string = "PSFF-P-00368719-N-03627591-KH-PN:PSFF-P-00368719-CFDNA-03627592-KH-PN"
# analysis_dir = "~/anchorageMountPointR2/P-00368719/PSFF-P-00368719-CFDNA-03627592-KH20191120-PN20191121_PSFF-P-00368719-N-03627591-KH20191120-PN20191121"
# outfile = "~/scratch/QC_overview.pdf"
# main_path = "~/anchorageMountPointR2"


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
msings_files = dir(path = Sys.glob(paste0(main_path, "/*/*/msings*")),
                   pattern = paste0("MSI_Analysis.txt$"), full.names = TRUE, recursive = TRUE)

# remove files from old pipeline (modified before 2019-04-30)
hsmetrics_files = hsmetrics_files[file.mtime(hsmetrics_files) > as.POSIXct("2019-04-30")]
markduplicates_files = markduplicates_files[file.mtime(markduplicates_files) > as.POSIXct("2019-04-30")]
insertsize_files = insertsize_files[file.mtime(insertsize_files) > as.POSIXct("2019-04-30")]
contest_files = contest_files[file.mtime(contest_files) > as.POSIXct("2019-04-30")]
msings_files = msings_files[file.mtime(msings_files) > as.POSIXct("2019-04-30")]


# read in the files
cat("Read in the files...\n")
HsMetrics = data.frame()
for (f in hsmetrics_files) {
  SAMP = strsplit(basename(f), split = "\\.")[[1]][1]
  DIR = dirname(dirname(dirname(dirname(f))))
  HsMetrics = rbind(HsMetrics, cbind(SAMP, DIR, read.table(f, skip = 6, nrow = 1, sep = "\t", 
                                                           header = TRUE, stringsAsFactors = FALSE), 
                                     stringsAsFactors = FALSE))
}
MarkDuplicates = data.frame()
for (f in markduplicates_files) {
  SAMP = sub("-markdups-metrics.txt", "", basename(f))
  DIR = dirname(dirname(dirname(dirname(f))))
  MarkDuplicates = rbind(MarkDuplicates, cbind(SAMP, DIR, read.table(f, skip = 6, nrow = 1, sep = "\t", 
                                                                     header = TRUE, stringsAsFactors = FALSE), 
                                               stringsAsFactors = FALSE))
}
InsertSize = data.frame()
InsertSize_histogram = data.frame(stringsAsFactors = FALSE)
for (f in insertsize_files) {
  SAMP = strsplit(basename(f), split = "\\.")[[1]][1]
  DIR = dirname(dirname(dirname(dirname(f))))
  InsertSize = rbind(InsertSize, cbind(SAMP, DIR, read.table(f, skip = 6, nrow = 1, sep = "\t", 
                                                             header = TRUE, stringsAsFactors = FALSE), 
                                       stringsAsFactors = FALSE))
  InsertSize_histogram = rbind(InsertSize_histogram, 
                               cbind(SAMP, DIR, read.table(f, skip = 10, sep = "\t",
                                                           header = TRUE, stringsAsFactors = FALSE),
                                     stringsAsFactors = FALSE))
}
ContEst = data.frame()
for (f in contest_files) {
  SAMP = strsplit(basename(f), split = "\\.")[[1]][1]
  DIR = dirname(dirname(f))
  ContEst = rbind(ContEst, cbind(SAMP, DIR, read.table(f, nrow = 1, sep = "\t", header = TRUE, stringsAsFactors = FALSE), 
                                 stringsAsFactors = FALSE))
}
msings = data.frame()
for (f in msings_files) {
  SAMP = sub("-nodups.MSI_Analysis.txt", "", basename(f))
  DIR = dirname(dirname(dirname(f)))
  msings = rbind(msings, cbind(SAMP, DIR, t(read.table(f, header = FALSE, nrows = 5, sep = "\t", stringsAsFactors = FALSE)))[2,], stringsAsFactors=FALSE)
}
colnames(msings) = c("SAMP", "DIR", read.table(f, header = FALSE, nrows = 5, sep = "\t", stringsAsFactors = FALSE)[,1])
msings$`msi status`[which(msings$msing_score<0.2)] = "NEG"  # use cut-off 0.2 for MSI-H (default 0.1 is too low)
msings$msing_score = as.numeric(msings$msing_score)

# add missing values in the insert size histogram table
for (d in unique(InsertSize_histogram$DIR)) {
  d1 = subset(InsertSize_histogram, DIR == d)
  for (s in unique(d1$SAMP)) {
    d2 = subset(d1, SAMP == s)
    missing = which(!seq(2,max(d2$insert_size)) %in% d2$insert_size) + 1
    if (length(missing)>0) {
      InsertSize_histogram = rbind(InsertSize_histogram,
                                   data.frame(SAMP = s, DIR = d, insert_size = missing, All_Reads.fr_count = 0))
    }
    # print(paste0("Sample ", s, " missing sizes: ", paste0(missing, collapse = ", ")))
  }
}



# merge the QC tables 
qc_merge = merge(merge(merge(merge(HsMetrics, MarkDuplicates, by = c("SAMP", "DIR")), InsertSize, by = c("SAMP", "DIR")), ContEst, by = c("SAMP", "DIR")),
                 msings, by = c("SAMP", "DIR"), all.x = TRUE)

# add sample type and capture kit
qc_merge$sample_type = sapply(strsplit(qc_merge$SAMP, split = "-"), "[", 4)
qc_merge$capture = sapply(strsplit(qc_merge$SAMP, split = "-"), "[", 7)
InsertSize_histogram$sample_type = sapply(strsplit(InsertSize_histogram$SAMP, split = "-"), "[", 4)
InsertSize_histogram$capture = sapply(strsplit(InsertSize_histogram$SAMP, split = "-"), "[", 7)

# label the samples of interest (soi) and analysis directory of interest (doi)
samples = strsplit(sample_string, split = ":")[[1]]
qc_merge$soi = qc_merge$SAMP %in% samples
qc_merge$doi = qc_merge$DIR == Sys.glob(analysis_dir)
InsertSize_histogram$soi = InsertSize_histogram$SAMP %in% samples
InsertSize_histogram$doi = InsertSize_histogram$DIR == Sys.glob(analysis_dir)


# create an ouput table for the samples of interest
soi_table = data.table(qc_merge)[i = soi&doi, 
                                 j =list(SAMP, MEAN_TARGET_COVERAGE, FOLD_ENRICHMENT, dedupped_on_bait_rate=ON_BAIT_BASES/PF_BASES_ALIGNED, 
                                         READ_PAIRS_EXAMINED, PERCENT_DUPLICATION, "contamination_%"=contamination, MEDIAN_INSERT_SIZE, msing_score)]
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
    geom_histogram(data = subset(qc_merge, soi), binwidth = sample_binwidth, fill = "blue", show.legend = FALSE) + # the sample of interest
    geom_histogram(data = subset(qc_merge, soi&doi), binwidth = sample_binwidth, fill = "red", show.legend = FALSE) + # the sample of interest
    scale_fill_manual(values = c("antiquewhite", "aliceblue", "lightpink", "palegreen")) +
    scale_y_continuous(breaks = ybreaks) +
    scale_x_continuous(name = x_string, breaks = xbreaks, minor_breaks = xbreaks_minor) +
    facet_wrap(~sample_type, ncol = 1) +
    ggtitle(paste0(title_string, " (binwidth ", binwidth, ")"), 
            subtitle = paste0("Red piles show present samples, blue piles show these samples if run earlier (binwidth ", sample_binwidth, ")"))
  print(p)
}

# plotting function to create desired scatter plots
my_scatter = function(x, y, xbreaks, ybreaks, x_string, y_string, title_string) {
  ggplot(qc_merge, aes(shape = capture)) +
    geom_point(aes_string(x = x, y = y), size = 3) +
    geom_point(data = subset(qc_merge, soi), aes_string(x = x, y = y), fill = "blue", size = 3, show.legend = FALSE) +
    geom_point(data = subset(qc_merge, soi&doi), aes_string(x = x, y = y), fill = "red", size = 3, show.legend = FALSE) +
    scale_alpha_manual(values = c(0.7, 1)) +
    scale_shape_manual(values = c(24, 25, 21, 22), guide = guide_legend(override.aes = list(fill = NA))) +
    scale_x_continuous(name = x_string, breaks = xbreaks) +
    scale_y_continuous(name = y_string, breaks = ybreaks) +
    facet_wrap(~sample_type, ncol = 1) +
    ggtitle(title_string, subtitle = "Red symbols show present samples, blue symbols show these samples if run earlier")
}

# Plot histograms and scatter plots
cat(paste0("Create plots (saved in ", outfile, ") ...\n"))
pdf(file = outfile, width=14)

# duplication vs read count scatter plot
my_scatter(x = "READ_PAIRS_EXAMINED", y = "PERCENT_DUPLICATION", xbreaks = seq(0, 1e12, 1e7), ybreaks = seq(0, 1, 0.1),
           x_string = "number of read pairs", y_string = "duplication rate", title_string = "Duplication rate vs Read count")

# coverage vs read count scatter plot
my_scatter(x = "READ_PAIRS_EXAMINED", y = "MEAN_TARGET_COVERAGE", xbreaks = seq(0, 1e12, 1e7), ybreaks = seq(0, 5000, 500),
           x_string = "number of read pairs", y_string = "mean target coverage", title_string = "Coverage vs Read count")

# duplication vs fold enrichment scatter plot
my_scatter(x = "FOLD_ENRICHMENT", y = "PERCENT_DUPLICATION", xbreaks = seq(0, 5000, 50), ybreaks = waiver(),
           x_string = "fold enrichment", y_string = "duplication rate", title_string = "Duplication rate vs Fold enrichment")

# on-bait rate vs duplication scatter plot
my_scatter(x = "PERCENT_DUPLICATION", y = "ON_BAIT_BASES/PF_BASES_ALIGNED", xbreaks = waiver(), ybreaks = waiver(),
           x_string = "duplication rate", y_string = "dedupped on-bait rate", title_string = "Dedupped on-bait rate rate vs Duplication rate")

# insert size histogram
ggplot(InsertSize_histogram, aes(x = insert_size, y = All_Reads.fr_count, group = interaction (SAMP, DIR),
                                 linetype = capture)) +
  geom_line(color = "black", alpha = 0.5) +
  geom_line(data = subset(InsertSize_histogram, soi&!doi), color = "blue", show.legend = FALSE) +
  geom_line(data = subset(InsertSize_histogram, soi&doi), color = "red", show.legend = FALSE) +
  geom_vline(xintercept = 167, alpha = 0.8, linetype = 3) +
  geom_label(data=data.frame(sample_type="CFDNA"), aes(x = 167, y = Inf), label = "167 bases, length of one chromatosome", 
             inherit.aes = FALSE, vjust = "top", hjust = 0, nudge_x = 2) +
  # scale_linetype(guide = guide_legend(override.aes = list(color = "black"))) +
  scale_x_continuous(name = "insert size", breaks = seq(0,2000,100)) +
  scale_y_continuous(name = "count") +
  facet_wrap(~sample_type, ncol = 1) +
  ggtitle("Insert size", subtitle = "Red lines show present samples, blue lines show these samples if run earlier.")

# contamination histogram
if (max(qc_merge$contamination) < 5) {
  xbreaks = seq(0.1, 6, 0.2)
  xbreaks_minor = seq(0, 6, 0.1)
} else {
  xbreaks = waiver()
  xbreaks_minor = waiver()
}
my_histogram(x = "contamination", binwidth = 0.1, sample_binwidth = 0.1, ybreaks = seq(0,nrow(qc_merge), 2),
             xbreaks = xbreaks, xbreaks_minor = xbreaks_minor, x_string = "contamination, %", 
             title_string = "Contamination")

# msings score vs read count scatter plot
my_scatter(x = "READ_PAIRS_EXAMINED", y = "msing_score", xbreaks = seq(0, 1e12, 1e7), ybreaks = waiver(),
             x_string = "number of read pairs", y_string = "mSINGS score", title_string = "mSINGS score vs Read count")



dev.off()

cat("Completed.\n")


