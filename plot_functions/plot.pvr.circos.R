#plot.pvr.circos <- function(filename, circos.dir="results/circos.pvr/"){
# Generate a circos plot showing the 10 samples with WGS/meth/gxn data
# This file creates a circos.conf file for use with circos software available at http://circos.ca
source('utility_functions/removeOutliers.R')
plot.pvr.circos <- function(filename, circos.dir=paste(getwd(), "results/circos.pvr/", sep="/")){
  
  circos.cmd <- "/Users/adam/Software/circos-0.69-6/bin/circos"
  if(!file.exists(circos.cmd)){
    stop("ERROR: circos not installed. Please see http://circos.ca/ and configure plot_functions/plot.genomic.circos.R")
  }
  dir.create(circos.dir, recursive = T, showWarnings = F)
  conf <- paste(circos.dir, "circos.conf", sep="")
  wd <- getwd() 
  
  gxn.circos <- paste(circos.dir, "gxn.circosdata.txt", sep="")
  gxn.circos.tcga <- paste(circos.dir, "gxn.circosdata.tcga.txt", sep="")
  meth.circos <- paste(circos.dir, "meth.circosdata.txt", sep="")
  meth.circos.tcga <- paste(circos.dir, "meth.circosdata.tcga.txt", sep="")
  cnas.prog.circos <- paste(circos.dir, "cna.prog.circosdata.txt", sep="")
  cnas.reg.circos <- paste(circos.dir, "cna.reg.circosdata.txt", sep="")
  
  # Here we plot PvR data around the genome:
  #   Differentially expressed genes (heatmap)
  #   DMRs (heatmap)
  #   CNAs?
  #   Mutation profiles - P outward, R inward
  
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
  track$gene <- as.character(track$gene)
  
  # GXN - only DE genes highlighted
  gdiff <- limmaCompare(gdata.d, gpheno.d, fdr_limit = 1)
  gxn.track <- track[which(as.character(track$gene) %in% rownames(gdiff)[1:100]),]
  gxn.track$value <- gdiff[as.character(gxn.track$gene),]$fc
  gxn.track$gene <- NULL
  gxn.track <- gxn.track[which(!is.na(gxn.track$value)),]
  write.table(gxn.track, sep="\t", quote=F, col.names=F, row.names = F, file=gxn.circos)
  
  # GXN for TCGA Ca vs Control - only DE genes highlighted
  source('utility_functions/limmaCompare.R')
  gdiff.tcga <- limmaCompare(tcga.gdata, tcga.gpheno, fdr_limit = 0.01)
  gxn.track.tcga <- track[which(track$gene %in% rownames(gdiff.tcga)[1:100]),]
  gxn.track.tcga$value <- gdiff.tcga[gxn.track.tcga$gene,]$fc
  gxn.track.tcga$gene <- NULL
  gxn.track.tcga <- gxn.track.tcga[which(!is.na(gxn.track.tcga$value)),]
  write.table(gxn.track.tcga, sep="\t", quote=F, col.names=F, row.names = F, file=gxn.circos.tcga)
  
  if(!exists("dmrs")){
    dmrs <- champ.DMR(beta=as.matrix(mdata.d), pheno=mpheno.d$Sample_Group, compare.group = c("Progressive", "Regressive"), method="ProbeLasso")
  }
  
  meth.track <- data.frame(
    chr=gsub("chr", "hs", dmrs$ProbeLassoDMR$seqnames),
    start=dmrs$ProbeLassoDMR$start,
    end=dmrs$ProbeLassoDMR$end,
    value=dmrs$ProbeLassoDMR$betaAv_Progressive - dmrs$ProbeLassoDMR$betaAv_Regressive
  )
  write.table(meth.track, sep="\t", quote=F, col.names=F, row.names = F, file=meth.circos)
  
  # Generate DMRs from tcga data for comparison:
  if(!exists("dmrs.tcga")){
    tcga.mpheno.tmp <- tcga.mpheno
    tcga.mpheno.tmp$Sample_Group <- make.names(tcga.mpheno.tmp$Sample_Group)
    # Impute to remove NAs
    tcga.mdata.imputed <- champ.impute(beta=as.matrix(tcga.mdata), SampleCutoff = 0.5, ProbeCutoff = 0.5, pd=tcga.mpheno.tmp)
    tcga.mpheno.tmp <- tcga.mdata.imputed$pd
    tcga.mdata.imputed <- tcga.mdata.imputed$beta
    # Strange bug - only use probes in package data(illumina450Gr) for DM (removes very few probes)
    data(illumina450Gr)
    tcga.mdata.imputed <- tcga.mdata.imputed[which(rownames(tcga.mdata.imputed) %in% names(illumina450Gr)),]
    # Find DMRs
    dmrs.tcga <- champ.DMR(beta=tcga.mdata.imputed, pheno=tcga.mpheno.tmp$Sample_Group, compare.group=c("TCGA.SqCC", "TCGA.Control"), method="ProbeLasso")
  }
  meth.track.tcga <- data.frame(
    chr=gsub("chr", "hs", dmrs.tcga$ProbeLassoDMR$seqnames),
    start=dmrs.tcga$ProbeLassoDMR$start,
    end=dmrs.tcga$ProbeLassoDMR$end,
    value=dmrs.tcga$ProbeLassoDMR$betaAv_TCGA.SqCC - dmrs.tcga$ProbeLassoDMR$betaAv_TCGA.Control
  )
  write.table(meth.track.tcga, sep="\t", quote=F, col.names=F, row.names = F, file=meth.circos.tcga)
  
  # Create a comparative CNA track, showing only regions which differ significantly
  # Due to low n numbers we can't find significant regions using a method like this:
  # cnadiff.data <- cnas.segmented[,wgs.pheno$name]
  # rownames(cnadiff.data) <- paste(cnas.segmented$chr, cnas.segmented$start, cnas.segmented$end, sep="-")
  # cnadiff <- limmaCompare(cnadiff.data, wgs.pheno, fdr_limit = 1)
  # df <- data.frame(
  #   row=rownames(cnadiff.data),
  #   p=unlist(lapply(rownames(cnadiff.data), function(x){
  #     wilcox.test(
  #       as.numeric(cnadiff.data[x, which(wgs.pheno$progression == 1)]),
  #       as.numeric(cnadiff.data[x, which(wgs.pheno$progression == 0)])
  #     )$p.value
  #   }))
  # )
  # df$p.adj <- p.adjust(df$p)
  
  # We will therefore plot prog and reg CNA tracks adjacently
  cnas.segmented.data <- cnas.segmented[,wgs.pheno$name] / wgs.pheno$ploidy
  sel <- which(cnas.segmented$chr %in% 1:22)
  cna.prog.track <- data.frame(
    chr=paste0("hs",as.character(cnas.segmented$chr[sel])),
    start=cnas.segmented$start[sel],
    end=cnas.segmented$end[sel],
    value=apply(cnas.segmented.data[sel,which(wgs.pheno$progression == 1)], 1, mean)
  )
  # Add colours manually
  cna.prog.track$col <- "color=vvlgrey"
  cna.prog.track$col[which(cna.prog.track$value < 0.9)] <- "color=lblue"
  cna.prog.track$col[which(cna.prog.track$value < 0.4)] <- "color=vdblue"
  cna.prog.track$col[which(cna.prog.track$value > 1.1)] <- "color=lred"
  cna.prog.track$col[which(cna.prog.track$value > 2)] <- "color=vdred"
  write.table(cna.prog.track, sep="\t", quote=F, col.names=F, row.names = F, file=cnas.prog.circos)
  # Repeat for reg
  cna.reg.track <- data.frame(
    chr=paste0("hs",as.character(cnas.segmented$chr[sel])),
    start=cnas.segmented$start[sel],
    end=cnas.segmented$end[sel],
    value=apply(cnas.segmented.data[sel,which(wgs.pheno$progression == 0)], 1, mean)
  )
  cna.reg.track$col <- "color=vvlgrey"
  cna.reg.track$col[which(cna.reg.track$value < 0.9)] <- "color=lblue"
  cna.reg.track$col[which(cna.reg.track$value < 0.4)] <- "color=vdblue"
  cna.reg.track$col[which(cna.reg.track$value > 1.1)] <- "color=lred"
  cna.reg.track$col[which(cna.reg.track$value > 2)] <- "color=vdred"
  write.table(cna.reg.track, sep="\t", quote=F, col.names=F, row.names = F, file=cnas.reg.circos)
  
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
              
# CIS DE genes
              <plot>
              type = histogram
              file = ",gxn.circos,"
              r1   = 0.975r
                r0   = 0.885r
              min = -3
              max = 3
              fill_color = blue
              extend_bin=no
              thickness=5
              <axes>
                show=yes
              <axis>
              color=vdgrey
              position=0.5r
              </axis>
              </axes>
              <backgrounds>
              <background>
              color=vvlred
              y0=0.5r
              </background>
              <background>
              color=vvlgreen
              y1=0.5r
              </background>
              </backgrounds>
              </plot>

# CIS methylation DMRs
              <plot>
              type = histogram
              file = ",meth.circos,"
              r1   = 0.875r
              r0   = 0.785r
              min = -0.5
              max = 0.5
              fill_color = yellow
              thickness = 5
              extend_bin=no
              <axes>
              show=yes
              <axis>
              color=vdgrey
              position=0.5r
              </axis>
              </axes>
              <backgrounds>
              <background>
              color=vvlyellow
              y0=0.5r
              </background>
              <background>
              color=vvlblue
              y1=0.5r
              </background>
              </backgrounds>
              </plot>

# CIS Prog CNAs
              <plot>
              type = heatmap
              file = ",cnas.prog.circos,"
              r1   = 0.775r
              r0   = 0.75r
              </plot>
# CIS Reg CNAs
              <plot>
              type = heatmap
              file = ",cnas.reg.circos,"
              r1   = 0.75r
              r0   = 0.725r
              </plot>

#################
# TCGA data
#################
# TCGA GXN DEs
              <plot>
              type = histogram
              file = ",gxn.circos.tcga,"
              r1   = 0.625r
              r0   = 0.535r
              min = -3
              max = 3
              fill_color = blue
              extend_bin=no
              thickness=5
              <axes>
              show=yes
              <axis>
              color=vdgrey
              position=0.5r
              </axis>
              </axes>
              <backgrounds>
              <background>
              color=vvlred
              y0=0.5r
              </background>
              <background>
              color=vvlgreen
              y1=0.5r
              </background>
              </backgrounds>
              </plot>

              
# TCGA methylation DMRs
              <plot>
              type = histogram
              file = ",meth.circos.tcga,"
              r1   = 0.525r
              r0   = 0.435r
              min = -0.5
              max = 0.5
              fill_color = yellow
              thickness = 5
              extend_bin=no
              <axes>
              show=yes
              <axis>
              color=vdgrey
              position=0.5r
              </axis>
              </axes>
              <backgrounds>
              <background>
              color=vvlyellow
              y0=0.5r
              </background>
              <background>
              color=vvlblue
              y1=0.5r
              </background>
              </backgrounds>
              </plot>

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
