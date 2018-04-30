plot.genomic.signatures <- function(filename){
  
  if(!exists("mut_mat")){
    stop("ERROR: mut_mat not found. Have you run signature analysis code from full_analysis.R?")
  }
  
  select <- which(rownames(fit_res$contribution) %in% sigs.to.analyse)
  # Order by number of mutations
  o <- order(colSums(mut_mat))
  
  
  # Repeat plots split into prog/reg
  #dev.new(width=5, height=5*length(which(wgs.pheno$progression == 0)) / dim(wgs.pheno)[1])
  png(paste0(filename, "_reg_rel.png"), width=960, height=960*length(which(wgs.pheno$progression == 0)) / dim(wgs.pheno)[1])
  plot_contribution(fit_res$contribution[select,o],
                    cancer_signatures[,select],
                    coord_flip = TRUE,
                    mode = "relative",
                    palette = col_vector[1:dim(cancer_signatures)[2]],
                    index = which(wgs.pheno[o,]$progression == 0)
  )
  dev.off()
  png(paste0(filename, "_reg_abs.png"), width=960, height=960*length(which(wgs.pheno$progression == 0)) / dim(wgs.pheno)[1])
  plot_contribution(fit_res$contribution[select,o],
                    cancer_signatures[,select],
                    coord_flip = TRUE,
                    mode = "absolute",
                    palette = col_vector[1:dim(cancer_signatures)[2]],
                    index = which(wgs.pheno[o,]$progression == 0)
  )
  dev.off()
  png(paste0(filename, "_prog_rel.png"), width=960, height=960*length(which(wgs.pheno$progression == 1)) / dim(wgs.pheno)[1])
  plot_contribution(fit_res$contribution[select,o],
                    cancer_signatures[,select],
                    coord_flip = TRUE,
                    mode = "relative",
                    palette = col_vector[1:dim(cancer_signatures)[2]],
                    index = which(wgs.pheno[o,]$progression == 1)
  )
  dev.off()
  png(paste0(filename, "_prog_abs.png"), width=960, height=960*length(which(wgs.pheno$progression == 1)) / dim(wgs.pheno)[1])
  plot_contribution(fit_res$contribution[select,o],
                    cancer_signatures[,select],
                    coord_flip = TRUE,
                    mode = "absolute",
                    palette = col_vector[1:dim(cancer_signatures)[2]],
                    index = which(wgs.pheno[o,]$progression == 1)
  )
  dev.off()
  
  
  
  # Plot a mean value
  # df <- data.frame(
  #   row.names = rownames(fit_res$contribution),
  #   cis.mean = as.numeric(apply(fit_res$contribution, 1, mean))
  # )
  # plot_contribution(df,
  #                   cancer_signatures,
  #                   coord_flip = FALSE,
  #                   mode = "absolute",
  #                   palette = col_vector[1:dim(cancer_signatures)[2]])
  # 
  # # Try TCGA data
  # tcga.pts <- unique(tcga.snvs$case_id)
  # tcga.vcfs <- list()
  # for(pt in tcga.pts){
  #   pt.snvs <- tcga.snvs[tcga.snvs$case_id == pt,]
  #   range <- GRanges(
  #     seqnames = gsub("chr", "", pt.snvs$Chromosome),
  #     ranges=IRanges(
  #       start=pt.snvs$Start_Position,
  #       end=pt.snvs$End_Position
  #     ),
  #     paramRangeID=NA,
  #     REF=DNAStringSet(
  #       pt.snvs$Reference_Allele
  #     ),
  #     ALT=DNAStringSetList(lapply(pt.snvs$Allele, function(x){ DNAStringSet(x)} )),
  #     QUAL=NA,
  #     FILTER="PASS"
  #   )
  #   tcga.vcfs[[length(tcga.vcfs) + 1]] <- range
  # }
  # tcga.vcfs <- GRangesList(tcga.vcfs)
  # 
  # tcga.ref_genome <- "BSgenome.Hsapiens.NCBI.GRCh38"
  # library(tcga.ref_genome, character.only = T)
  # 
  # # Plot the spectrum across all samples
  # tcga.type_occurrences <- mut_type_occurrences(tcga.vcfs, tcga.ref_genome)
  # plot_spectrum(tcga.type_occurrences)
  # 
  # tcga.mut_mat <- mut_matrix(vcf_list = tcga.vcfs[1:10], ref_genome = tcga.ref_genome)
  # tcga.fit_res <- fit_to_signatures(tcga.mut_mat, cancer_signatures)
  # 
  # # Plot contribution barplot
  # plot_contribution(tcga.fit_res$contribution[select,o],
  #                   cancer_signatures[,select],
  #                   coord_flip = FALSE,
  #                   mode = "relative",
  #                   palette = col_vector[1:dim(cancer_signatures)[2]]
  # )
  # # Repeat as absolute plot
  # plot_contribution(tcga.fit_res$contribution[select,o],
  #                   cancer_signatures[,select],
  #                   coord_flip = FALSE,
  #                   mode = "absolute",
  #                   palette = col_vector[1:dim(cancer_signatures)[2]]
  # )
  
  pdf(paste0(filename, "_all.pdf"))
  # Plot contribution barplot
  plot_contribution(fit_res$contribution[select,o],
                    cancer_signatures[,select],
                    coord_flip = TRUE,
                    mode = "relative",
                    palette = col_vector[1:dim(cancer_signatures)[2]]
  )
  # Repeat as absolute plot
  plot_contribution(fit_res$contribution[select,o],
                    cancer_signatures[,select],
                    coord_flip = TRUE,
                    mode = "absolute",
                    palette = col_vector[1:dim(cancer_signatures)[2]]
  )
  # Plot as a signature heatmap
  plot_cosine_heatmap(cos_sim_samples_signatures,
                      col_order = cosmic_order,
                      cluster_rows = TRUE)
  
  # Plot prog vs reg for individual signatures
  # fit_res$contribution has approximate number of mutations per signature. colSums gives roughly the correct mutation count.
  plotdata <- wgs.pheno
  for(sig in sigs.to.analyse){
    # Compare absolute mutations in this signature
    plotdata[,sig] <- fit_res$contribution[sig,wgs.pheno$name]
    
    compare.fn(dependent_variable = sig, compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years"),
               random_effects = "Patient", modelinfo=plotdata, title=paste0(sig, " (abs)"))
    
    # Compare relative proportion of mutations in this signature
    sig.rel <- paste0(sig, ".rel")
    plotdata[,sig.rel] <- 100*fit_res$contribution[sig, wgs.pheno$name] / colSums(fit_res$contribution[,wgs.pheno$name])
    
    compare.fn(dependent_variable = sig.rel, compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years"),
               random_effects = "Patient", modelinfo=plotdata, title=paste0(sig, " (rel)"))
  }
  
  
  dev.off() 
  
  
}