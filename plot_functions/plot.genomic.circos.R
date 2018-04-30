source('utility_functions/removeOutliers.R')
plot.genomic.circos <- function(filename, circos.dir=paste(getwd(), "results/circos/", sep="/")){
  # Generate a circos plot comparing TCGA and CIS data
  # This file creates a circos.conf file for use with circos software available at http://circos.ca
  
  circos.cmd <- "/Users/adam/Software/circos-0.69-6/bin/circos"
  if(!file.exists(circos.cmd)){
    stop("ERROR: circos not installed. Please see http://circos.ca/ and configure plot_functions/plot.genomic.circos.R")
  }
  
  dir.create(circos.dir, recursive = T, showWarnings = F)
  
  conf <- paste(circos.dir, "circos.conf", sep="")
  circos.tcga <- paste(circos.dir, "circos.tcga.txt", sep="")
  circos.cis <- paste(circos.dir, "circos.cis.txt", sep="")
  circos.drivers <- paste(circos.dir, "circos.drivers.txt", sep="")
  circos.driver.labels <- paste(circos.dir, "circos.driver.labels.txt", sep="")
  circos.cnas <- paste(circos.dir, "circos.cnas.txt", sep="")
  circos.gxn <- paste(circos.dir, "circos.gxn.txt", sep="")
  circos.meth <- paste(circos.dir, "circos.meth.txt", sep="")
  circos.cnas.tcga <- paste(circos.dir, "circos.cnas.tcga.txt", sep="")
  circos.gxn.tcga <- paste(circos.dir, "circos.gxn.tcga.txt", sep="")
  circos.meth.tcga <- paste(circos.dir, "circos.meth.tcga.txt", sep="")
  
  hmcol <- colorRampPalette(c("Green","Black","Red"))(256)
  hmcol2 <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(256)
  hmcol.list <- paste( apply(col2rgb(hmcol), 2, function(x){paste(x, collapse=",")}) , collapse="),(")
  hmcol.list <- paste0("(", hmcol.list, ")")
  
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
  tcga.freq <- data.frame(
    chr=genes$seqnames,
    start=genes$start,
    end=genes$end,
    value=0,
    gene=as.character(genes$SYMBOL)
  )
  cis.freq <- tcga.freq
  drivers.plot <- tcga.freq
  driver.labels <- tcga.freq
  cis.gxn <- tcga.freq
  cis.meth <- tcga.freq
  samples.muts <- list()
  for(name in wgs.pheno$name){
    samples.muts[[name]] <- tcga.freq
  }
  
  # Add data for TCGA mutation rate
  sel <- which(names(tcga.snvs.rates) %in% tcga.freq$gene)
  tcga.freq$value[match(names(tcga.snvs.rates)[sel], tcga.freq$gene)] <- as.numeric(tcga.snvs.rates)[sel]
  tcga.freq$gene <- NULL
  write.table(tcga.freq, sep="\t", quote=F, col.names=F, row.names = F, file=circos.tcga)
  
  # Add data for CIS mutation rates
  muts.freq <- apply(muts, 1, function(x){length(which(x > 0))}) / dim(muts)[2]
  sel <- which(names(muts.freq) %in% cis.freq$gene)
  cis.freq$value[match(names(muts.freq)[sel], cis.freq$gene)] <- as.numeric(muts.freq)[sel]
  cis.freq$gene <- NULL
  write.table(cis.freq, sep="\t", quote=F, col.names=F, row.names = F, file=circos.cis)
  
  # Drivers are simply 1 for drivers, 0 otherwise
  driver.genes <- read.csv('resources/driver_genes.csv', stringsAsFactors = F)
  driver.genes <- unique(driver.genes$Gene[which(driver.genes$Cancer %in% c("LUSC", "PANCAN"))])
  drivers.plot$value[which(drivers.plot$gene %in% driver.genes)] <- 1
  drivers.plot$gene <- NULL
  write.table(drivers.plot, sep="\t", quote=F, col.names=F, row.names = F, file=circos.drivers)
  
  # Label driver mutations
  driver.labels <- driver.labels[which(driver.labels$gene %in% driver.genes),]
  driver.labels$value <- driver.labels$gene
  driver.labels$gene <- NULL
  write.table(driver.labels, sep="\t", quote=F, col.names=F, row.names = F, file=circos.driver.labels)
  
  # Make a copy number heatmap from cnas.segmented
  cnas.track <- cnas.segmented[,1:3]
  cnas.track$value <- apply(cnas.segmented[,wgs.pheno$name] / wgs.pheno$ploidy, 1, mean)
  cnas.track$chr <- paste("hs", cnas.track$chr, sep="")
  # Remove sex chromosomes
  cnas.track <- cnas.track[which(!(cnas.track$chr %in% c('hsX', "hsY"))),]
  # Add colours manually
  cnas.track$col <- "color=vvlgrey"
  cnas.track$col[which(cnas.track$value < 0.9)] <- "color=lblue"
  cnas.track$col[which(cnas.track$value < 0.4)] <- "color=vdblue"
  cnas.track$col[which(cnas.track$value > 1.1)] <- "color=lred"
  cnas.track$col[which(cnas.track$value > 2)] <- "color=vdred"
  write.table(cnas.track, sep="\t", quote=F, col.names=F, row.names = F, file=circos.cnas)
  
  # Repeat for TCGA copy number
  tcga.cnas.track <- tcga.cnas.segmented.mean
  colnames(tcga.cnas.track) <- c("chr", "start", "end", "value")
  tcga.cnas.track$chr <- paste("hs", tcga.cnas.track$chr, sep="")
  # Remove sex chromosomes
  tcga.cnas.track <- tcga.cnas.track[which(!(tcga.cnas.track$chr %in% c('hsX', "hsY"))),]
  # Add colours manually
  tcga.cnas.track$col <- "color=vvlgrey"
  tcga.cnas.track$col[which(tcga.cnas.track$value < 0.9)] <- "color=lblue"
  tcga.cnas.track$col[which(tcga.cnas.track$value < 0.4)] <- "color=vdblue"
  tcga.cnas.track$col[which(tcga.cnas.track$value > 1.1)] <- "color=lred"
  tcga.cnas.track$col[which(tcga.cnas.track$value > 2)] <- "color=vdred"
  write.table(tcga.cnas.track, sep="\t", quote=F, col.names=F, row.names = F, file=circos.cnas.tcga)
  
  # These TCGA values are relative copy numbers so choose appropriate cutoffs for plot similarity based on:
  # par(mfrow=c(1,2))
  # plot(cnas.track$value, ylim=c(-1.5,3.5))
  # abline(h=0.4)
  # abline(h=0.9)
  # abline(h=1.1)
  # abline(h=2)
  # plot(tcga.cnas.track$value, ylim=c(-0.5,0.5))
  # abline(h=-0.3)
  # abline(h=-0.4)
  # abline(h=-0.3)
  # abline(h=0.3)
  # abline(h=0.4)
  # par(mfrow=c(1,1))
  #
  
  # Make a gxn track
  cis.gxn$gene <- as.character(cis.gxn$gene)
  cis.gxn <- cis.gxn[which(cis.gxn$gene %in% rownames(gdata)),]
  cis.gxn$value <- as.numeric(apply(gdata[cis.gxn$gene,], 1, function(x){ mean(as.numeric(x)) }))
  cis.gxn$gene <- NULL
  # Cap outliers for a better plot
  cis.gxn$value <- removeOutliers(cis.gxn$value)
  write.table(cis.gxn, sep="\t", quote=F, col.names=F, row.names = F, file=circos.gxn)
  
  
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
  for(i in 1:dim(wgs.pheno)[1]){
    name <- wgs.pheno$name[i]
    x <- muts.all[which(muts.all$patient == name),]
    x$type <- as.character(x$type)
    # Filter out non-coding mutations
    x <- x[which(x$type %in% names(mut.cols)),]
    if(dim(x)[1] == 0){next}
    x <- data.frame(
      chr=paste("hs", x$chr, sep=""),
      start=x$start,
      end=x$end,
      value=sample(1:100, dim(x)[1], replace=T), # This value is irrelevant as we set the colour manually
      col=paste("color=", mut.cols[x$type], sep="")
    )
    write.table(x, sep="\t", quote=F, col.names=F, row.names = F, file=paste(circos.dir, name, "_muts.txt", sep=""))
  }
  
  for(i in 1:dim(wgs.pheno)[1]){
    name <- wgs.pheno$name[i]
    # Make a CNA track
    x <- cnas.segmented[,c("chr", "start", "end", name)]
    x[,name] <- x[,name] / wgs.pheno$ploidy[i]
    x$chr <- paste("hs", x$chr, sep="")
    x <- x[which(!(x$chr %in% c("hsX", "hsY"))),]
    x$col <- 'color=vvlgrey'
    x$col[which(x[,name] < 1)] <- "color=lblue"
    x$col[which(x[,name] == 0)] <- "color=vdblue"
    x$col[which(x[,name] > 1)] <- "color=lred"
    x$col[which(x[,name] >= 2)] <- "color=vdred"
    write.table(x, sep="\t", quote=F, col.names=F, row.names = F, file=paste(circos.dir, name, "_cnas.txt", sep=""))
  }
  
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
                
                # TCGA mutation frequencies
                <plot>
                type=histogram
                file=",circos.tcga,"
                r1   = 0.88r
                r0   = 0.81r
                min = 0
                max = 1
                extend_bin=no
                <axes>
                show=yes
                <axis>
                color=vlgrey
                spacing=0.2r
                </axis>
                </axes>
                </plot>
                
                # CIS mutation frequencies
                <plot>
                type=histogram
                file=",circos.cis,"
                r1   = 0.76r
                r0   = 0.69r
                orientation = in
                min = 0
                max = 1
                extend_bin=no
                <axes>
                show=yes
                <axis>
                color=vlgrey
                spacing=0.2r
                </axis>
                </axes>
                </plot>
                
                #<plot>
                #type=histogram
                #file=",circos.drivers,"
                #r1   = 0.79r
                #r0   = 0.78r
                #color = red
                #fill_color=red
                #</plot>
                
                # Labels for drivers
                <plot>
                type = text
                color            = black
                file             = ",circos.driver.labels,"
                r0 = 0.88r
                r1 = 0.95r
                
                show_links     = yes
                link_dims      = 4p,4p,8p,4p,4p
                link_thickness = 2p
                link_color     = red
                
                label_size   = 20p
                label_font   = condensed
                
                padding  = 0p
                rpadding = 0p
                
                label_snuggle = yes
                
                </plot>
                
                # Heatmap for CIS CNAs
                <plot>
                type  = heatmap
                file  = ",circos.cnas,"
                r1    = 0.78r
                r0    = 0.76r
                </plot>
                
                # Heatmap for TCGA CNAs - AWAITING SANGER DATA
                #<plot>
                #type  = heatmap
                #file  = ",circos.cnas.tcga,"
                #r1    = 0.81r
                #r0    = 0.79r
                #</plot>
                
                # Heatmap for CIS GXNs
                #<plot>
                #type  = heatmap
                #file  = ",circos.gxn,"
                #r1    = 0.69r
                #r0    = 0.67r
                #color = ",hmcol.list,"
                #</plot>
                
                # Heatmap for CIS meth data
                #<plot>
                #type  = heatmap
                #file  = ",circos.meth,"
                #r1    = 0.67r
                #r0    = 0.65r
                #color = blue,lblue,yellow,lred,red
                #</plot>
                
                ", sep="")
  
  # Add mutation tracks per-sample
  # Placed at 0.64-0.65; 0.62-0.63 etc.
  for(i in 1:dim(wgs.pheno)[1]){
    f <- paste(circos.dir, wgs.pheno$name[i], "_muts.txt", sep="")
    if(!file.exists(f)){next}
    text <- paste(text, "
                  <plot>
                  type=heatmap
                  file=",f,"
                  r1=",0.67-0.02*i-0.0025,"r
                  r0=",0.67-0.02*i-0.01,"r
                  </plot>
                  ", sep="")
  }
  # Add CNA tracks per-sample, interlaced with above
  for(i in 1:dim(wgs.pheno)[1]){
    f <- paste(circos.dir, wgs.pheno$name[i], "_cnas.txt", sep="")
    if(!file.exists(f)){next}
    text <- paste(text, "
                  <plot>
                  type=heatmap
                  file=",f,"
                  r1=",0.67-0.02*i-0.01,"r
                  r0=",0.67-0.02*i-0.02,"r
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
  # We need different behaviour for relative and absolute file paths
  if(substr(filename, 1, 1) == "/"){
    file.copy("circos.png", filename, overwrite = T)
  }else{
    file.copy("circos.png", paste(wd, filename, sep="/"), overwrite = T)
  }
  setwd(wd)
}