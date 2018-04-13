# Generate a circos plot comparing TCGA and CIS data
# This file creates a circos.conf file for use with circos software available at http://circos.ca

conf <- "results/circos/circos.conf"
circos.tcga <- "results/circos/circos.tcga.txt"
circos.cis <- "results/circos/circos.cis.txt"
circos.drivers <- "results/circos/circos.drivers.txt"
circos.driver.labels <- "results/circos/circos.driver.labels.txt"
circos.cnas <- "results/circos/circos.cnas.txt"

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

# TCGA data frame
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
samples.muts <- list()
for(name in wgs.pheno$name){
  samples.muts[[name]] <- tcga.freq
}

# Add data for TCGA
sel <- which(names(tcga.snvs.rates) %in% tcga.freq$gene)
tcga.freq$value[match(names(tcga.snvs.rates)[sel], tcga.freq$gene)] <- as.numeric(tcga.snvs.rates)[sel]
tcga.freq$gene <- NULL
write.table(tcga.freq, sep="\t", quote=F, col.names=F, row.names = F, file=circos.tcga)

# Add data for CIS
sel <- which(names(muts.freq) %in% cis.freq$gene)
cis.freq$value[match(names(muts.freq)[sel], cis.freq$gene)] <- as.numeric(muts.freq)[sel]
cis.freq$gene <- NULL
write.table(cis.freq, sep="\t", quote=F, col.names=F, row.names = F, file=circos.cis)

# Drivers are simply 1 for drivers, 0 otherwise
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
cnas.track$col[which(cnas.track$value < 0.7)] <- "color=lblue"
cnas.track$col[which(cnas.track$value < 0.4)] <- "color=vdblue"
cnas.track$col[which(cnas.track$value > 1.3)] <- "color=lred"
cnas.track$col[which(cnas.track$value > 2)] <- "color=vdred"
write.table(cnas.track, sep="\t", quote=F, col.names=F, row.names = F, file=circos.cnas)


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
  write.table(x, sep="\t", quote=F, col.names=F, row.names = F, file=paste("results/circos/", name, "_muts.txt", sep=""))
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
  write.table(x, sep="\t", quote=F, col.names=F, row.names = F, file=paste("results/circos/", name, "_cnas.txt", sep=""))
}

# Create the circos config file
text <- paste("
karyotype = data/karyotype/karyotype.human.txt

<ideogram>
  
  <spacing>
    default = 0.005r
  </spacing>
  
  radius    = 0.95r
  thickness = 20p
  fill      = yes

  show_label       = yes
  # see etc/fonts.conf for list of font names
  label_font       = default 
  label_radius     = 1r + 75p
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
</plot>

# CIS mutation frequencies
<plot>
type=histogram
file=",circos.cis,"
r1   = 0.78r
r0   = 0.71r
orientation = in
</plot>

<plot>
type=histogram
file=",circos.drivers,"
r1   = 0.81r
r0   = 0.78r
color = red
fill_color=red
</plot>

# Labels for drivers
<plot>
type = text
color            = black
file             = ",circos.driver.labels,"
r0 = 0.88r
r1 = 0.95r
              
show_links     = no
#link_dims      = 4p,4p,8p,4p,4p
#link_thickness = 2p
#link_color     = red
              
label_size   = 18p
label_font   = condensed
              
padding  = 0p
rpadding = 0p

label_snuggle = yes
              
</plot>

# Heatmap for CIS CNAs
<plot>
type  = heatmap
file  = ",circos.cnas,"
r1    = 0.71r
r0    = 0.68r

</plot>
", sep="")

# Add mutation tracks per-sample
# Placed at 0.64-0.65; 0.62-0.63 etc.
for(i in 1:dim(wgs.pheno)[1]){
  f <- paste("results/circos/", wgs.pheno$name[i], "_muts.txt", sep="")
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
  f <- paste("results/circos/", wgs.pheno$name[i], "_cnas.txt", sep="")
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
<<include etc/housekeeping.conf>>
",sep="")

write(text, file=conf)

system(paste("/Users/adam/Software/circos-0.69-6/bin/circos -conf", conf))