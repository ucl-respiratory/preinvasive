plot.genomic.signatures <- function(filename){
  
  if(!exists("mut_mat")){
    stop("ERROR: mut_mat not found. Have you run signature analysis code from full_analysis.R?")
  }
  
  select <- which(rownames(fit_res$contribution) %in% sigs.to.analyse)
  # Order by number of mutations
  o <- order(colSums(mut_mat))
  
  
  # Repeat plots split into prog/reg
  #dev.new(width=5, height=5*length(which(wgs.pheno$progression == 0)) / dim(wgs.pheno)[1])
  pdf(paste0(filename, "_reg_rel.pdf"), width=9.60, height=9.60*length(which(wgs.pheno$progression == 0)) / dim(wgs.pheno)[1])
  plot_contribution(fit_res$contribution[select,o],
                    cancer_signatures[,select],
                    coord_flip = TRUE,
                    mode = "relative",
                    palette = col_vector[1:dim(cancer_signatures)[2]],
                    index = which(wgs.pheno[o,]$progression == 0)
  )
  dev.off()
  pdf(paste0(filename, "_reg_abs.pdf"), width=9.60, height=9.60*length(which(wgs.pheno$progression == 0)) / dim(wgs.pheno)[1])
  plot_contribution(fit_res$contribution[select,o],
                    cancer_signatures[,select],
                    coord_flip = TRUE,
                    mode = "absolute",
                    palette = col_vector[1:dim(cancer_signatures)[2]],
                    index = which(wgs.pheno[o,]$progression == 0)
  )
  dev.off()
  pdf(paste0(filename, "_prog_rel.pdf"), width=9.60, height=9.60*length(which(wgs.pheno$progression == 1)) / dim(wgs.pheno)[1])
  plot_contribution(fit_res$contribution[select,o],
                    cancer_signatures[,select],
                    coord_flip = TRUE,
                    mode = "relative",
                    palette = col_vector[1:dim(cancer_signatures)[2]],
                    index = which(wgs.pheno[o,]$progression == 1)
  )
  dev.off()
  pdf(paste0(filename, "_prog_abs.pdf"), width=9.60, height=9.60*length(which(wgs.pheno$progression == 1)) / dim(wgs.pheno)[1])
  plot_contribution(fit_res$contribution[select,o],
                    cancer_signatures[,select],
                    coord_flip = TRUE,
                    mode = "absolute",
                    palette = col_vector[1:dim(cancer_signatures)[2]],
                    index = which(wgs.pheno[o,]$progression == 1)
  )
  dev.off()
  
  
  
  
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
  plotdata <- plotdata[-which(plotdata$query.reg == 1),]
  for(sig in sigs.to.analyse){
    # Compare absolute mutations in this signature
    plotdata[,sig] <- fit_res$contribution[sig,plotdata$name]
    
    compare.fn(dependent_variable = sig, compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years", "purity"),
               random_effects = "Patient", modelinfo=plotdata, title=paste0(sig, " (abs)"))
    
    # Compare relative proportion of mutations in this signature
    sig.rel <- paste0(sig, ".rel")
    plotdata[,sig.rel] <- 100*fit_res$contribution[sig, plotdata$name] / colSums(fit_res$contribution[,plotdata$name])
    
    compare.fn(dependent_variable = sig.rel, compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years", "purity"),
               random_effects = "Patient", modelinfo=plotdata, title=paste0(sig, " (rel)"))
  }
  
  
  dev.off() 
  
  
  # TCGA comparison:
  pdf(paste0(filename, "_with_tcga.pdf"), width=960, height=960*3 / dim(wgs.pheno)[1])
  
  # Plot mean signature contributions for CIS and TCGA
  cont.combined <- data.frame(row.names = row.names(fit_res_exonic$contribution))
  cont.combined$cis  <- apply(fit_res_exonic$contribution, 1, mean)
  cont.combined$tcga <- apply(fit_res_tcga$contribution, 1, mean)
  plot_contribution(as.matrix(cont.combined[select,]), cancer_signatures[,select], coord_flip = T, mode="relative", palette=col_vector[1:dim(cancer_signatures)[2]])

  dev.off()
  
  # Plot cosine similarity matrix
  pdf(paste0(filename, "_tcga_similarity.pdf"))
  cs <- cos_sim_matrix(tcga_mut_mat, mut_mat_exonic)
  rownames(cs) <- 1:dim(cs)[1]
  colnames(cs) <- paste0(colnames(cs), " (", c("R", "P")[wgs.pheno$progression+1], ")")
  plot_cosine_heatmap(cs, cluster_rows = T)
  dev.off()
  
}