mute <- suppressPackageStartupMessages
mute(library(ggplot2))
mute(library(stringr))
mute(library(magrittr))
mute(library(WriteXLS))
mute(library(tidyr))
mute(library(dplyr))
mute(library(plotly))
mute(library(Signac))
mute(library(Seurat))
mute(library(cluster))
mute(library(clustree))
mute(library(mclust))
mute(library(cowplot))
mute(library(gridExtra))
mute(library(ggrastr))
mute(library(viridis))
mute(library(GenomicRanges))
mute(library(GenomeInfoDb))
mute(library(BSgenome.Mmusculus.UCSC.mm10))
mute(library(EnsDb.Hsapiens.v86))
mute(library(data.table))
mute(library(patchwork))
mute(library(foreach))
mute(library(doParallel))
mute(library(Matrix))
mute(library(optparse))


##!!!! please change the species parameter before running the script


peak_dir = "../data/peak"

## define helper function

readSummits <- function(file){
    df <- data.frame(readr::read_tsv(file,
                                     col_names = c("chr","start","end","name","score")))
    df <- df[,c(1,2,3,5)] #do not keep name column it can make the size really large
    return(makeGRangesFromDataFrame(df = df,
                                    keep.extra.columns = TRUE,
                                    starts.in.df.are.0based = TRUE))
}

clusterGRanges <- function(gr, filter = TRUE, by = "score", decreasing = TRUE){
      gr <- sort(sortSeqlevels(gr))
      r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
      o <- findOverlaps(gr,r)
      mcols(gr)$cluster <- subjectHits(o)
      gr <- gr[order(mcols(gr)[,by], decreasing = decreasing),]
      gr <- gr[!duplicated(mcols(gr)$cluster),]
      gr <- sort(sortSeqlevels(gr))
      mcols(gr)$cluster <- NULL
      return(gr)
}

nonOverlappingGRanges <- function(gr, by = "score", decreasing = TRUE, verbose = FALSE){
    stopifnot(by %in% colnames(mcols(gr)))
    i <-  0
    gr_converge <- gr
    while(length(gr_converge) > 0){
      if(verbose){
        message(".", appendLF = FALSE)
      }
      i <-  i + 1
      gr_selected <- clusterGRanges(gr = gr_converge, filter = TRUE, by = by, decreasing = decreasing)
      gr_converge <- subsetByOverlaps(gr_converge ,gr_selected, invert=TRUE) #blacklist selected gr
      if(i == 1){ #if i=1 then set gr_all to clustered
        gr_all <- gr_selected
      }else{
        gr_all <- c(gr_all, gr_selected)
      } 
    }
    if(verbose){
      message("\nSelected ", length(gr_all), " from ", length(gr))
    }
    gr_all <- sort(sortSeqlevels(gr_all))
    return(gr_all)

}


countInsertions <- function(query, fragments, by = "RG"){
    #Count By Fragments Insertions
    inserts <- c(
        GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), RG = mcols(fragments)[,by]),
        GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), RG = mcols(fragments)[,by])
    )
    by <- "RG"
    overlapDF <- DataFrame(findOverlaps(query, inserts, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
    overlapDF$name <- mcols(inserts)[overlapDF[, 2], by]
    overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
    #Calculate Overlap Stats
    inPeaks <- table(overlapDF$name)
    total <- table(mcols(inserts)[, by])
    total <- total[names(inPeaks)]
    frip <- inPeaks / total
    #Summarize
    sparseM <- Matrix::sparseMatrix(
        i = overlapTDF[, 1],
        j = overlapTDF[, 4],
        x = rep(1, nrow(overlapTDF)),
        dims = c(length(query), length(unique(overlapDF$name))))
    colnames(sparseM) <- unique(overlapDF$name)
    return(sparseM)
}




## Make Non-Overlapping Peak Set
blacklist_file = "../../Blacklists/mm10-blacklist.v2.bed.gz" ## https://github.com/Boyle-Lab/Blacklist/tree/master/lists
blacklist <- rtracklayer::import.bed(blacklist_file)

chromSizes <- GRanges(names(seqlengths(BSgenome.Mmusculus.UCSC.mm10)), 
                      IRanges(1, seqlengths(BSgenome.Mmusculus.UCSC.mm10)))
chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes,
                                                    pruning.mode = "coarse")
peaks_files <- list.files(peaks_dir,
                          pattern = "\\_summits.bed", 
                          full.names = TRUE)

gr_list <- GenomicRangesList(lapply(peaks_files, function(x){
        extended_summits <- readSummits(x) %>%
            resize(., width = 2 * 250 + 1, fix = "center") %>%
            subsetByOverlaps(.,chromSizes,type="within") %>%
              subsetByOverlaps(.,blacklist,invert=TRUE) %>%
              nonOverlappingGRanges(., by="score", decreasing=TRUE)
        extended_summits <- extended_summits[order(extended_summits$score, 
                                                   decreasing=TRUE)]
        extended_summits <- head(extended_summits, 100000)
        mcols(extended_summits)$scoreQuantile <-trunc(rank(mcols(extended_summits)$score)) / length(mcols(extended_summits)$score)
        extended_summits
        }))
unionPeaks <- nonOverlappingGRanges(unlist(gr_list),
                                    by = "scoreQuantile", 
                                    decreasing = TRUE)
unionPeaks <- sort(sortSeqlevels(unionPeaks))

unionPeaks <- unionPeaks[seqnames(unionPeaks) %in% paste0("chr",c(1:19,"X"))]
unionPeaks <- keepSeqlevels(unionPeaks, paste0("chr",c(1:19,"X")))

df <- data.frame(seqnames=seqnames(unionPeaks),
                 starts=start(unionPeaks)-1,
                 ends=end(unionPeaks),
                 name=paste0("peaks_", 1:length(unionPeaks)),
                 score=score(unionPeaks))

write.table(df, file = paste0(dir.out, "/unionPeaks.bed"), 
            sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE)


## Count Insertions to create peak by cell matrix
fragment.path <- file.path(dir_1round,"Fragments/fragments.tsv")
cell_df <- read.csv(file.path(dir_1round, "singlecell.csv"), row.names=1)
cells <- rownames(cell_df)

message("Reading in fragment files...", date())
fragments <- data.frame(readr::read_tsv(fragment.path, col_names=FALSE))
fragments <- GRanges(
    seqnames = fragments[,1],
    IRanges(fragments[,2]+1, fragments[,3]),
    RG = fragments[,4],
    N = fragments[,5]
)

fragments <- fragments[fragments$RG %in% cells]

#Create Counts matirx
message("creating countmatrix...", sample, "  ", date())
counts <- countInsertions(unionPeaks, fragments, by = "RG")
rownames(counts) <- paste(seqnames(unionPeaks),
                          start(unionPeaks),
                          end(unionPeaks),
                          sep="_")
counts <- counts[rowSums(counts) > 0, ]
Ypeaks <- which(grepl("chrY", rownames(counts)))
if(length(Ypeaks) > 0){
  counts <- counts[-(Ypeaks),]
}

message("done countmatrix...", sample, "  ", date())
saveRDS(counts, file = file.path(dir.out, "unionPeaks_matrix.Rds"))



mtx_dir <- file.path(dir.out, "filtered_peak_bc_matrix")
if(!dir.exists(mtx_dir)){
  dir.create(mtx_dir)
}


writeMM(counts, file = file.path(mtx_dir, "matrix.mtx"))
barcodes <- as.data.frame(colnames(counts))
peaks <- as.data.frame(stringr::str_split_fixed(rownames(counts), ":|-|_", 3))

write.table(barcodes, file =  file.path(mtx_dir,"barcodes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(peaks, file =  file.path(mtx_dir,"peaks.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## Session information

sessionInfo()

