#plot.overlaps.circos <- function(filename, circos.dir="results/circos.overlap/"){
  # Generate a circos plot showing the 10 samples with WGS/meth/gxn data
  # This file creates a circos.conf file for use with circos software available at http://circos.ca
  
  circos.cmd <- "/Users/adam/Software/circos-0.69-6/bin/circos"
  if(!file.exists(circos.cmd)){
    stop("ERROR: circos not installed. Please see http://circos.ca/ and configure plot_functions/plot.genomic.circos.R")
  }
  dir.create(circos.dir, recursive = T, showWarnings = F)
  conf <- paste(circos.dir, "circos.conf", sep="")
  
  # Load overlapping sample data
  samples <- read.xls("resources/overlap.samples.xlsx", stringsAsFactors = F)
  samples <- samples[which(samples$Whole.Genome.Sequencing == "YES" & samples$Gene.expression == "YES" & samples$Methylation == "YES"),]
  
  samples.gdata <- gdata[,samples$Sample.Number..GXN.]
  samples.mdata <- mdata[,make.names(samples$Sample.Number..Meth.)]
  samples.cnas  <- cnas.segmented[,c("chr", "start", "end",samples$Sample.Number..WGS.)]
  samples.muts <- muts.all[which(muts.all$patient %in% samples$Sample.Number..WGS.),]
  
  # Create three data track files:
  # TCGA mutation frequencies by gene
  # CIS mutation frequencies by gene
  # Driver mutation locations
  
  library(GenomicRanges)
  library(Homo.sapiens)
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  # Get genes
  refgenome <- BSgenome.Hsapiens.UCSC.hg19
  gene_coord <- genes(Homo.sapiens, columns=c("GENEID", "SYMBOL"))
  genes <- data.frame(gene_coord)
  
  # Only look at chr1-23
  genes$seqnames <- as.character(genes$seqnames)
  genes <- genes[grep("^chr([0-9]+|[X-Y])$", genes$seqnames),]
  # Manipulate to circos format
  genes$seqnames <- gsub("chr", "hs", genes$seqnames)
  
  # TCGA mutation rate data frame
  track <- data.frame(
    chr=genes$seqnames,
    start=genes$start,
    end=genes$end,
    value=0,
    gene=as.character(genes$SYMBOL)
  )
  
  # Make a per-patient mutation track
  mut.cols <- list(
    'missense'='black',
    'start_lost'='black',
    'stop_lost'='black',
    'frameshift'='vlblue',
    'nonsense'='dgreen',
    'silent'='dgrey',
    'splice_region'="lgreen",
    "ess_splice"="lgreen"
  )
  for(i in 1:dim(samples)[1]){
    name <- samples$Sample.Number..WGS.[i]
    x <- muts.all[which(muts.all$patient == name),]
    x$type <- as.character(x$type)
    # Filter out non-coding mutations
    x <- x[which(x$type %in% names(mut.cols)),]
    if(dim(x)[1] != 0){
      x <- data.frame(
        chr=paste("hs", x$chr, sep=""),
        start=x$start,
        end=x$end,
        value=sample(1:100, dim(x)[1], replace=T), # This value is irrelevant as we set the colour manually
        col=paste("color=", mut.cols[x$type], sep="")
      )
      write.table(x, sep="\t", quote=F, col.names=F, row.names = F, file=paste(circos.dir, name, "_muts.txt", sep=""))
    }
    
    # Make a CNA track
    x <- cnas.segmented[,c("chr", "start", "end", name)]
    x[,name] <- x[,name] / wgs.pheno$ploidy[match(name, wgs.pheno$name)]
    x$chr <- paste("hs", x$chr, sep="")
    x <- x[which(!(x$chr %in% c("hsX", "hsY"))),]
    x$col <- 'color=vvlgrey'
    x$col[which(x[,name] < 1)] <- "color=lblue"
    x$col[which(x[,name] == 0)] <- "color=vdblue"
    x$col[which(x[,name] > 1)] <- "color=lred"
    x$col[which(x[,name] >= 2)] <- "color=vdred"
    write.table(x, sep="\t", quote=F, col.names=F, row.names = F, file=paste(circos.dir, name, "_cnas.txt", sep=""))
    
    # Make a GXN track
    x <- track
    x <- x[which(x$gene %in% rownames(gdata)),]
    x$value <- gdata[x$gene,samples$Sample.Number..GXN.[i]]
    x <- x[which(!is.na(x$value)),]
    x$gene <- NULL
    write.table(x, sep="\t", quote=F, col.names=F, row.names = F, file=paste(circos.dir, name, "_gxn.txt", sep=""))
  }
  
  
  # GXN track per patient:
  
  
  # Create the circos config file
  text <- paste("
                karyotype = data/karyotype/karyotype.human.txt
                
                <ideogram>
                
                <spacing>
                default = 0.005r
                </spacing>
                
                radius    = 0.95r
                thickness = 2p
                fill      = yes
                color     = dgrey
                
                show_label       = yes
                # see etc/fonts.conf for list of font names
                label_font       = default 
                label_radius     = 1r + 5p
                label_size       = 30
                label_parallel   = yes
                
                </ideogram>
                
                <plots>
                
                ", sep="")
  
  # Add mutation tracks per-sample
  # Placed at 0.64-0.65; 0.62-0.63 etc.
  for(i in 1:dim(samples)[1]){
    f <- paste(circos.dir, samples$Sample.Number..WGS.[i], "_muts.txt", sep="")
    if(!file.exists(f)){next}
    text <- paste(text, "
                  <plot>
                  type=heatmap
                  file=",f,"
                  r1=",0.95-0.08*i-0.0025,"r
                  r0=",0.95-0.08*i-0.02,"r
                  </plot>
                  ", sep="")
  }
  # Add CNA tracks per-sample, interlaced with above
  for(i in 1:dim(samples)[1]){
    f <- paste(circos.dir, samples$Sample.Number..WGS.[i], "_cnas.txt", sep="")
    if(!file.exists(f)){next}
    text <- paste(text, "
                  <plot>
                  type=heatmap
                  file=",f,"
                  r1=",0.95-0.08*i-0.02,"r
                  r0=",0.95-0.08*i-0.04,"r
                  </plot>
                  ", sep="")
  }
  # Add GXN tracks per-sample, interlaced with above
  for(i in 1:dim(samples)[1]){
    f <- paste(circos.dir, samples$Sample.Number..WGS.[i], "_gxn.txt", sep="")
    if(!file.exists(f)){next}
    text <- paste(text, "
                  <plot>
                  type=histogram
                  file=",f,"
                  r1=",0.95-0.08*i-0.04,"r
                  r0=",0.95-0.08*i-0.06,"r
                  </plot>
                  ", sep="")
  }
  
  text <- paste(text, "
                </plots>
                ################################################################
                # The remaining content is standard and required. It is imported 
                # from default files in the Circos distribution.
                #
                # These should be present in every Circos configuration file and
                # overridden as required. To see the content of these files, 
                # look in etc/ in the Circos distribution.
                
                <image>
                # Included from Circos distribution.
                <<include etc/image.conf>>
                </image>
                
                # RGB/HSV color definitions, color lists, location of fonts, fill patterns.
                # Included from Circos distribution.
                <<include etc/colors_fonts_patterns.conf>>
                
                # Debugging, I/O an dother system parameters
                # Included from Circos distribution.
                # This file need to be edited to set max_points_per_track to 250000 (for the TCGA CNA track)
                <<include etc/housekeeping.conf>>
                ",sep="")
  
  write(text, file=conf)
  
  wd <- getwd() 
  setwd(circos.dir)
  system(circos.cmd) # conf file should be in this directory already
  # file.copy("circos.png", paste(wd, filename, sep="/"))
  setwd(wd)
#}