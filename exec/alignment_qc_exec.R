# Executable script called by alignment QC shell scripts 
# 1st CL argument should be project name (i.e. name of dir in /ag_sanders/work/data)
# 2nd CL argument shiuld be no. of threads to use (e.g. 32)

# load libraries
library(breakpointR)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(doParallel)


# define functions
load2object <- function (filename){
    if (file.exists(filename)) 
        return(eval(parse(text = load(filename))))
    cat(paste("error! the file ", filename, 
              " was not found :("))
    NULL
}

plot_metric = function(plot_input = combined_qc_stats, colname2plot=NA, plotting_colours = metric_cols,
                       plot_type = "violin"){
    if(is.na(colname2plot)) stop("plot_metric() define the colname you want to plot!")
    if(!plot_type %in% c("boxplot","violin")) stop("plot_metric() plot_type must = boxplot or violin")
    if(plot_type=="violin"){
        ggplot(plot_input) + 
            geom_jitter(aes(x = "", y = !!as.name(colname2plot)), size = 3.5, width = .25) + 
            geom_violin(aes(x = "", y = !!as.name(colname2plot), fill = colname2plot), alpha = 0.5) + 
            labs(x = "", y = colname2plot) + 
            scale_fill_manual(values = plotting_colours) + 
            theme_bw() +
            theme(legend.position = "none", axis.text.y = element_text(size=20), axis.title.y = element_text(size=20),
                  panel.border = element_blank(), axis.line = element_line(size = 1),
                  axis.ticks.y = element_line(size = 1), axis.ticks.x = element_blank(), panel.grid = element_blank())
    }
    if(plot_type=="boxplot"){
        ggplot(plot_input) + 
            geom_jitter(aes(x = "", y = !!as.name(colname2plot)), size = 3.5, width = .25) + 
            geom_boxplot(aes(x = "", y = !!as.name(colname2plot), fill = colname2plot), alpha = 0.5, outlier.shape = NA) + 
            labs(x = "", y = colname2plot) + 
            scale_fill_manual(values = plotting_colours) + 
            theme_bw() +
            theme(legend.position = "none", axis.text.y = element_text(size=25), axis.title.y = element_text(size=25),
                  panel.border = element_blank(), axis.line = element_line(size = 1),
                  axis.ticks.y = element_line(size = 1), axis.ticks.x = element_blank(), panel.grid = element_blank())
    }
}

bpR_calcs = function (bamfile, ID = basename(bamfile), pairedEndReads = TRUE, 
                      chromosomes = NULL, windowsize = 2e+06, binMethod = "size", 
                      multi.sizes = NULL, trim = 10, peakTh = 0.33, zlim = 3.291, 
                      background = 0.05, min.mapq = 10, pair2frgm = FALSE, filtAlt = FALSE, 
                      genoT = "fisher", minReads = 10, maskRegions = NULL, conf = 0.99) {
    requireNamespace(c("breakpointR","GenomicRanges","cowplot","breakpointRdata","methods",
                                            "utils","grDevices","stats","S4Vectors","GenomeInfoDb",
                                            "IRanges","Rsamtools","GenomicAlignments","ggplot2",
                                            "BiocGenerics","gtools","doParallel","foreach"))
    suppressWarnings(fragments <- readBamFileAsGRanges(bamfile, pairedEndReads = T, chromosomes = mychroms,
                                      min.mapq = 10, pair2frgm = F, filtAlt = F))
    
    chroms.in.data <- seqlevels(fragments)
    if (is.null(mychroms)) {
        chromosomes <- chroms.in.data
    }
    chroms2use <- intersect(mychroms, chroms.in.data)
    fragments <- keepSeqlevels(fragments, value = mychroms, 
                               pruning.mode = "coarse")
    # }
    filename <- basename(bamfile)
    message("Working on ", filename)
    maskRegions = NULL
    if (!is.null(maskRegions)) {
        mask <- IRanges::findOverlaps(maskRegions, fragments)
        fragments <- fragments[-S4Vectors::subjectHits(mask)]
    }
    reads.all.chroms <- fragments
    deltas.all.chroms <- GenomicRanges::GRangesList()
    GenomeInfoDb::seqlevels(deltas.all.chroms) <- GenomeInfoDb::seqlevels(fragments)
    GenomeInfoDb::seqlengths(deltas.all.chroms) <- GenomeInfoDb::seqlengths(fragments)
    breaks.all.chroms <- GenomicRanges::GRangesList()
    GenomeInfoDb::seqlevels(breaks.all.chroms) <- GenomeInfoDb::seqlevels(fragments)
    GenomeInfoDb::seqlengths(breaks.all.chroms) <- GenomeInfoDb::seqlengths(fragments)
    confint.all.chroms <- GenomicRanges::GRangesList()
    GenomeInfoDb::seqlevels(confint.all.chroms) <- GenomeInfoDb::seqlevels(fragments)
    GenomeInfoDb::seqlengths(confint.all.chroms) <- GenomeInfoDb::seqlengths(fragments)
    counts.all.chroms <- GenomicRanges::GRangesList()
    GenomeInfoDb::seqlevels(counts.all.chroms) <- GenomeInfoDb::seqlevels(fragments)
    GenomeInfoDb::seqlengths(counts.all.chroms) <- GenomeInfoDb::seqlengths(fragments)
    
    binMethod = "size"
    windowsize = 2e+6
    multi.sizes = NULL
    peakTh = 0.33
    
    for (chr in unique(seqnames(fragments))) {
        # message("  Working on chromosome ", chr)
        fragments.chr <- fragments[seqnames(fragments) == chr]
        fragments.chr <- GenomeInfoDb::keepSeqlevels(fragments.chr, 
                                                     chr)
        # print("    calculating deltaWs ...")
        if (binMethod == "size") {
            if (windowsize > GenomeInfoDb::seqlengths(fragments)[chr]) {
                # stopTimedMessage(ptm)
                next
            }
            tiles <- unlist(GenomicRanges::tileGenome(seqlengths(fragments)[chr], 
                                                      tilewidth = windowsize))
            counts <- GenomicRanges::countOverlaps(tiles, fragments.chr)
            reads.per.window <- max(10, round(mean(counts[counts > 
                                                              0], trim = 0.05)))
            if (!is.null(multi.sizes) & length(multi.sizes) > 
                1) {
                dw <- deltaWCalculatorVariousWindows(frags = fragments.chr, 
                                                     reads.per.window = reads.per.window, multi.sizes = multi.sizes)
            } else {
                dw <- deltaWCalculator(frags = fragments.chr, 
                                       reads.per.window = reads.per.window)
            }
        }
        
        deltaWs <- dw[seqnames(dw) == chr]
        breaks <- suppressMessages(breakSeekr(deltaWs, trim = 10, 
                                              peakTh = peakTh, zlim = 3.291))
        if (length(breaks) > 0) {
            maxiter <- 10
            iter <- 1
            utils::flush.console()
            newBreaks <- GenotypeBreaks(breaks = breaks, fragments = fragments, 
                                        background = 0.05, minReads = 10, 
                                        genoT = "fisher")
            prev.breaks <- breaks
            breaks <- newBreaks
            while (length(prev.breaks) > length(newBreaks) && 
                   !is.null(breaks)) {
                utils::flush.console()
                iter <- iter + 1
                newBreaks <- GenotypeBreaks(breaks = breaks, 
                                            fragments = fragments, background = 0.05, 
                                            minReads = 10, genoT = "fisher")
                prev.breaks <- breaks
                breaks <- newBreaks
                if (iter == maxiter) {
                    break
                }
            }
        }
	if(length(breaks)==0) newBreaks = NULL
        if (is.null(newBreaks)) {
            Ws <- length(which(as.logical(strand(fragments.chr) == 
                                              "-")))
            Cs <- length(which(as.logical(strand(fragments.chr) == 
                                              "+")))
            chrRange <- GenomicRanges::GRanges(seqnames = chr, 
                                               ranges = IRanges(start = 1, end = GenomeInfoDb::seqlengths(fragments.chr)[chr]))
            counts <- cbind(Ws, Cs)
            mcols(chrRange) <- counts
            seqlengths(chrRange) <- GenomeInfoDb::seqlengths(fragments.chr)[chr]
            WC.ratio <- (chrRange$Ws - chrRange$Cs)/sum(c(chrRange$Ws, 
                                                          chrRange$Cs))
            if (WC.ratio > 0.8) {
                state <- "ww"
            }
            else if (WC.ratio < -0.8) {
                state <- "cc"
            }
            else if (WC.ratio < 0.2 & WC.ratio > -0.2) {
                state <- "wc"
            }
            else {
                state <- "?"
            }
            chrRange$states <- state
            suppressWarnings(counts.all.chroms[[chr]] <- chrRange)
            confint <- NULL
        }else {
            breaks.strand <- newBreaks
            strand(breaks.strand) <- "*"
            breakrange <- GenomicRanges::gaps(breaks.strand)
            breakrange <- breakrange[strand(breakrange) == "*"]
            strand(breakrange) <- "-"
            Ws <- GenomicRanges::countOverlaps(breakrange, fragments)
            strand(breakrange) <- "+"
            Cs <- GenomicRanges::countOverlaps(breakrange, fragments)
            strand(breakrange) <- "*"
            concat.states <- paste(newBreaks$genoT, collapse = "-")
            split.states <- unlist(strsplit(concat.states, "-"))
            states.idx <- seq(from = 1, to = length(split.states), 
                              by = 2)
            states.idx <- c(states.idx, length(split.states))
            states <- split.states[states.idx]
            counts <- cbind(Ws, Cs)
            mcols(breakrange) <- counts
            breakrange$states <- states
            suppressWarnings(counts.all.chroms[[chr]] <- breakrange)
        }
        if (length(deltaWs) > 0) {
            deltas.all.chroms[[chr]] <- deltaWs[, "deltaW"]
        }
        if (length(newBreaks) > 0) {
            breaks.all.chroms[[chr]] <- newBreaks
        }
    }
    deltas.all.chroms <- unlist(deltas.all.chroms, use.names = FALSE)
    breaks.all.chroms <- unlist(breaks.all.chroms, use.names = FALSE)
    confint.all.chroms <- unlist(confint.all.chroms, use.names = FALSE)
    counts.all.chroms <- unlist(counts.all.chroms, use.names = FALSE)
    ww <- counts.all.chroms[counts.all.chroms$states == "ww"]
    cc <- counts.all.chroms[counts.all.chroms$states == "cc"]
    if (length(ww) > 0) {
        ww$Cs <- ww$Cs + 1
        ww$Ws <- ww$Ws + 1
        bg.estim.ww <- sum(ww$Cs)/sum(ww$Ws)
    }else {
        bg.estim.ww <- 0
    }
    if (length(cc) > 0) {
        cc$Ws <- cc$Ws + 1
        cc$Cs <- cc$Cs + 1
        bg.estim.cc <- sum(cc$Ws)/sum(cc$Cs)
    } else {
        bg.estim.cc <- 0
    }
    bg.estim <- mean(bg.estim.ww, bg.estim.cc)
    chr.lengths <- seqlengths(fragments)[!is.na(seqlengths(fragments))]
    if (length(chr.lengths) > 0) {
        tiles <- unlist(GenomicRanges::tileGenome(chr.lengths, 
                                                  tilewidth = 1e+06))
        counts <- GenomicRanges::countOverlaps(tiles, fragments)
        reads.MB <- round(median(counts))
    } else {
        reads.MB <- NA
    }
    red.frags <- GenomicRanges::reduce(fragments)
    perc.cov <- (sum(as.numeric(width(red.frags)))/sum(as.numeric(seqlengths(red.frags)))) * 
        100
    perc.cov <- round(perc.cov, digits = 2)
    libname = gsub(".sorted.mdup.bam.RData|.sort.mdup.bam.RData|.sorted.mdup.bam|.sort.mdup.bam",
                   "",filename)
    library.metrics <- c(library = libname, background.estimate = as.numeric(bg.estim), med.reads.per.MB = as.numeric(reads.MB), 
                         perc.coverage = as.numeric(perc.cov))
    return(library.metrics)
}

# load command line args
project_name = NA ; n_threads = NA
args = commandArgs(trailingOnly=T)
project_name = args[1] 
n_threads= as.numeric(args[2])

if(is.na(project_name)) stop('set project_name in command line argument')
if(is.na(n_threads)) stop('set n_threads in command line argument')

print("terminal arguments loaded", quote = F)
print(paste("project_name =",project_name), quote = F)
print(paste("n_threads =",n_threads), quote = F)

# run breakpointR function to calc background
bpr_indir = file.path("//fast/groups/ag_sanders/work/data",project_name,"bam")
mychroms = paste0(rep("chr",1,24), c(1:22,"X","Y"))

print(paste("Running breakpointR on bam files in", bpr_indir), quote = F)

cl <- parallel::makeCluster(n_threads)
doParallel::registerDoParallel(cl)
message("cluster set up")

bpr_inputfiles = list.files(bpr_indir)[!grepl(".bai",list.files(bpr_indir))]
bpR_stats_par = foreach(mybam = bpr_inputfiles, .combine=rbind, .packages = "breakpointR") %dopar% {
    bpR_calcs(file.path(bpr_indir, mybam))
}
parallel::stopCluster(cl)
message("cluster closed")

bpR_stats = bpR_stats_par %>% 
    as.data.frame() %>% 
    mutate(across(contains("."), .fns = as.numeric))
rownames(bpR_stats) = NULL


# load mosaicatcher output and calculate entropy and spikines
mosaiccatcher_load = read.table(file.path("//fast/groups/ag_sanders/work/data",project_name,"qc/mosaicatcher/counts.txt.gz"),
                               header = T) 
samplecolname = ifelse("cell" %in% names(mosaiccatcher_load),"cell","sample")

mosaiccatcher_out = mosaiccatcher_load %>% 
    mutate(bin_dp = w + c) %>% 
    group_by(!!as.name(samplecolname)) %>% 
    mutate(total_counts = sum(bin_dp),
           Spikiness = sum(abs(diff(bin_dp))) / total_counts, # spikiness equation from AneuFinder (not exported)
           n = bin_dp / total_counts,
           Entropy = -sum(n*log(n), na.rm = T)) %>% # entropy also from AneuFinder (Bakker et al. 2016 Genome Biology)
    select(library=!!as.name(samplecolname), Spikiness, Entropy) %>% 
    distinct()

print(paste("Mosaicatcher data loaded from", file.path("//fast/groups/ag_sanders/work/data",project_name,"qc/mosaicatcher/counts.txt.gz")), quote = F)

# load coverage
coveragedir=file.path("//fast/groups/ag_sanders/work/data",project_name,"qc/alignment_stats/meandepthbychrom")
for(myfile in list.files(coveragedir)){
    currfile = read.table(file.path(coveragedir,myfile), header = F)
    currdf = data.frame(library = gsub("_mean_depth_bychrom.txt","",myfile),
                        mean_depth = mean(currfile$V4))
    if(myfile == list.files(coveragedir)[1]){
        coverage_df = currdf
    }else{
        coverage_df = rbind(coverage_df, currdf)
    }
}
print(paste("Depth of coverage data loaded from", coveragedir), quote = F)


# combine all qc stats
combined_qc_stats = read.table(file.path("//fast/groups/ag_sanders/work/data",project_name,"qc/alignment_stats/all_samples_qc_metrics.txt"),
           header = T, sep = "\t") %>% 
    left_join(coverage_df) %>% 
    left_join(bpR_stats) %>% 
    left_join(mosaiccatcher_out) %>% 
    arrange(library)

resdir = file.path("//fast/groups/ag_sanders/work/data",project_name,"qc/alignment_stats")
write.table(combined_qc_stats, quote = F, sep = "\t",row.names = F,
            file = file.path(resdir,"all_samples_all_qc_metrics.txt"))

print(paste("All QC metrics table saved to", file.path(resdir,"all_samples_all_qc_metrics.txt")), quote = F)



# plot metrics 
print(paste("Plotting QC metrics to", resdir), quote = F)

# change colnames
combined_qc_stats = combined_qc_stats %>% 
    rename(`GC content (%)` = gc_content, `No. of reads` = n_reads, `No. of reads mapped` = n_reads_mapped,
           `Mean insert size (bp)` = mean_insert_size, `Mean depth` = mean_depth,
           Background = background.estimate, `Median reads per Mb` = med.reads.per.MB,
           `Coverage (%)` = perc.coverage)

# define colours for plotting
metric_names = names(combined_qc_stats)[!names(combined_qc_stats) %in% c("library")] # as not paired end data 
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
metric_cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
metric_cols = metric_cols[sample.int(length(metric_cols),length(metric_names))]; names(metric_cols) = metric_names

pdf(file.path(resdir,"qc_metrics_boxplot.pdf"), width = 5, height = 6)
for(metric in metric_names){
    myplot = plot_metric(colname2plot = metric, plot_type = "boxplot")
    print(myplot)
}
dev.off()

# pdf(file.path(resdir,"qc_metrics_violins.pdf"), width = 5, height = 6)
# for(metric in metric_names){
#     myplot = plot_metric(colname2plot = metric, plot_type = "violin")
#     print(myplot)
# }
# dev.off()
