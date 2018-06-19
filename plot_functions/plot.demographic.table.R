# Table 1 - demographics
plot.demographic.table <- function(filename){
  df <- data.frame(matrix(nrow=15, ncol=12))
  rownames(df) <- c(
    "Patients", "Lesions profiled", 
    "Age at specimen profiled: mean", "Age at specimen profiled: median", "Age at specimen profiled: range",
    "Smoking History: mean", "Smoking History: median", "Smoking History: range",
    "Gender: F", "Gender: M",
    "COPD: N", "COPD: Y", "COPD: ?",
    "Previous history of lung cancer: Y", "Previous history of lung cancer: N"
  )
  colnames(df) <- c(
    "gxn.d.prog","gxn.d.reg","gxn.v.prog","gxn.v.reg",
    "methyl.d.prog", "methyl.d.reg",  "methyl.d.cont", "methyl.v.prog", "methyl.v.reg",  "methyl.v.cont", 
    "seq.prog", "seq.reg"
  )
  
  
  
  ################################################################################################
  # Gene Expression
  ################################################################################################
  gpheno.d.p   <- gpheno.d[which(gpheno.d$progression == 1),]
  sel <- which(duplicated(gpheno.d.p$Patient))
  if(length(sel) > 0){ gpheno.d.p.u <- gpheno.d.p[-sel,] } else { gpheno.d.p.u <- gpheno.d.p }
  gpheno.d.r   <- gpheno.d[which(gpheno.d$progression == 0),]
  sel <- which(duplicated(gpheno.d.r$Patient))
  if(length(sel) > 0){ gpheno.d.r.u <- gpheno.d.r[-sel,] } else { gpheno.d.r.u <- gpheno.d.r }
  
  gpheno.v.p   <- gpheno.v[which(gpheno.v$progression == 1),]
  sel <- which(duplicated(gpheno.v.p$Patient))
  if(length(sel) > 0){ gpheno.v.p.u <- gpheno.v.p[-sel,] } else { gpheno.v.p.u <- gpheno.v.p }
  gpheno.v.r   <- gpheno.v[which(gpheno.v$progression == 0),]
  sel <- which(duplicated(gpheno.v.r$Patient))
  if(length(sel) > 0){ gpheno.v.r.u <- gpheno.v.r[-sel,] } else { gpheno.v.r.u <- gpheno.v.r }
  df["Patients", "gxn.d.prog"] <- dim(gpheno.d.p.u)[1]
  df["Patients", "gxn.d.reg"]  <- dim(gpheno.d.r.u)[1]
  df["Patients", "gxn.v.prog"] <- dim(gpheno.v.p.u)[1]
  df["Patients", "gxn.v.reg"]  <- dim(gpheno.v.r.u)[1]
  
  df["Lesions profiled", "gxn.d.prog"] <- length(which(gpheno.d$progression == 1))
  df["Lesions profiled", "gxn.d.reg"]  <- length(which(gpheno.d$progression == 0))
  df["Lesions profiled", "gxn.v.prog"] <- length(which(gpheno.v$progression == 1))
  df["Lesions profiled", "gxn.v.reg"]  <- length(which(gpheno.v$progression == 0))
  
  df["Age at specimen profiled: mean", "gxn.d.prog"] <- mean(gpheno.d.p.u$Age.at.specimen.profiled, na.rm=T)
  df["Age at specimen profiled: mean", "gxn.d.reg"]  <- mean(gpheno.d.r.u$Age.at.specimen.profiled, na.rm=T)
  df["Age at specimen profiled: mean", "gxn.v.prog"] <- mean(gpheno.v.p.u$Age.at.specimen.profiled, na.rm=T)
  df["Age at specimen profiled: mean", "gxn.v.reg"]  <- mean(gpheno.v.r.u$Age.at.specimen.profiled, na.rm=T)
  
  df["Age at specimen profiled: median", "gxn.d.prog"] <- median(gpheno.d.p.u$Age.at.specimen.profiled, na.rm=T)
  df["Age at specimen profiled: median", "gxn.d.reg"]  <- median(gpheno.d.r.u$Age.at.specimen.profiled, na.rm=T)
  df["Age at specimen profiled: median", "gxn.v.prog"] <- median(gpheno.v.p.u$Age.at.specimen.profiled, na.rm=T)
  df["Age at specimen profiled: median", "gxn.v.reg"]  <- median(gpheno.v.r.u$Age.at.specimen.profiled, na.rm=T)
  
  df["Age at specimen profiled: range", "gxn.d.prog"] <- paste(range(gpheno.d.p.u$Age.at.specimen.profiled, na.rm=T),collapse="-")
  df["Age at specimen profiled: range", "gxn.d.reg"]  <- paste(range(gpheno.d.r.u$Age.at.specimen.profiled, na.rm=T),collapse="-")
  df["Age at specimen profiled: range", "gxn.v.prog"] <- paste(range(gpheno.v.p.u$Age.at.specimen.profiled, na.rm=T),collapse="-")
  df["Age at specimen profiled: range", "gxn.v.reg"]  <- paste(range(gpheno.v.r.u$Age.at.specimen.profiled, na.rm=T),collapse="-")
  
  df["Smoking History: mean", "gxn.d.prog"] <- mean(gpheno.d.p.u$Pack.years, na.rm=T)
  df["Smoking History: mean", "gxn.d.reg"]  <- mean(gpheno.d.r.u$Pack.years, na.rm=T)
  df["Smoking History: mean", "gxn.v.prog"] <- mean(gpheno.v.p.u$Pack.years, na.rm=T)
  df["Smoking History: mean", "gxn.v.reg"]  <- mean(gpheno.v.r.u$Pack.years, na.rm=T)
  
  df["Smoking History: median", "gxn.d.prog"] <- median(gpheno.d.p.u$Pack.years, na.rm=T)
  df["Smoking History: median", "gxn.d.reg"]  <- median(gpheno.d.r.u$Pack.years, na.rm=T)
  df["Smoking History: median", "gxn.v.prog"] <- median(gpheno.v.p.u$Pack.years, na.rm=T)
  df["Smoking History: median", "gxn.v.reg"]  <- median(gpheno.v.r.u$Pack.years, na.rm=T)
  
  df["Smoking History: range", "gxn.d.prog"] <- paste(range(gpheno.d.p.u$Pack.years, na.rm=T),collapse="-")
  df["Smoking History: range", "gxn.d.reg"]  <- paste(range(gpheno.d.r.u$Pack.years, na.rm=T),collapse="-")
  df["Smoking History: range", "gxn.v.prog"] <- paste(range(gpheno.v.p.u$Pack.years, na.rm=T),collapse="-")
  df["Smoking History: range", "gxn.v.reg"]  <- paste(range(gpheno.v.r.u$Pack.years, na.rm=T),collapse="-")
  
  df["Gender: M", "gxn.d.prog"] <- length(which(gpheno.d.p.u$Gender == "M"))
  df["Gender: M", "gxn.d.reg"]  <- length(which(gpheno.d.r.u$Gender == "M"))
  df["Gender: M", "gxn.v.prog"] <- length(which(gpheno.v.p.u$Gender == "M"))
  df["Gender: M", "gxn.v.reg"]  <- length(which(gpheno.v.r.u$Gender == "M"))
  df["Gender: F", "gxn.d.prog"] <- length(which(gpheno.d.p.u$Gender == "F"))
  df["Gender: F", "gxn.d.reg"]  <- length(which(gpheno.d.r.u$Gender == "F"))
  df["Gender: F", "gxn.v.prog"] <- length(which(gpheno.v.p.u$Gender == "F"))
  df["Gender: F", "gxn.v.reg"]  <- length(which(gpheno.v.r.u$Gender == "F"))
  
  df["COPD: Y", "gxn.d.prog"] <- length(which(gpheno.d.p.u$COPD == "Y"))
  df["COPD: Y", "gxn.d.reg"]  <- length(which(gpheno.d.r.u$COPD == "Y"))
  df["COPD: Y", "gxn.v.prog"] <- length(which(gpheno.v.p.u$COPD == "Y"))
  df["COPD: Y", "gxn.v.reg"]  <- length(which(gpheno.v.r.u$COPD == "Y"))
  df["COPD: N", "gxn.d.prog"] <- length(which(gpheno.d.p.u$COPD == "N"))
  df["COPD: N", "gxn.d.reg"]  <- length(which(gpheno.d.r.u$COPD == "N"))
  df["COPD: N", "gxn.v.prog"] <- length(which(gpheno.v.p.u$COPD == "N"))
  df["COPD: N", "gxn.v.reg"]  <- length(which(gpheno.v.r.u$COPD == "N"))
  df["COPD: ?", "gxn.d.prog"] <- length(which(gpheno.d.p.u$COPD == ""))
  df["COPD: ?", "gxn.d.reg"]  <- length(which(gpheno.d.r.u$COPD == ""))
  df["COPD: ?", "gxn.v.prog"] <- length(which(gpheno.v.p.u$COPD == ""))
  df["COPD: ?", "gxn.v.reg"]  <- length(which(gpheno.v.r.u$COPD == ""))
  
  df["Previous history of lung cancer: Y", "gxn.d.prog"] <- length(which(gpheno.d.p.u$Prev.History.of.LC == 1))
  df["Previous history of lung cancer: Y", "gxn.d.reg"]  <- length(which(gpheno.d.r.u$Prev.History.of.LC == 1))
  df["Previous history of lung cancer: Y", "gxn.v.prog"] <- length(which(gpheno.v.p.u$Prev.History.of.LC == 1))
  df["Previous history of lung cancer: Y", "gxn.v.reg"]  <- length(which(gpheno.v.r.u$Prev.History.of.LC == 1))
  df["Previous history of lung cancer: N", "gxn.d.prog"] <- length(which(gpheno.d.p.u$Prev.History.of.LC == 0))
  df["Previous history of lung cancer: N", "gxn.d.reg"]  <- length(which(gpheno.d.r.u$Prev.History.of.LC == 0))
  df["Previous history of lung cancer: N", "gxn.v.prog"] <- length(which(gpheno.v.p.u$Prev.History.of.LC == 0))
  df["Previous history of lung cancer: N", "gxn.v.reg"]  <- length(which(gpheno.v.r.u$Prev.History.of.LC == 0))
  
  
  ################################################################################################
  # Methylation
  ################################################################################################
  mpheno.d.p   <- mpheno.d[which(mpheno.d$Sample_Group == "Progressive"),]
  sel <- which(duplicated(mpheno.d.p$Patient))
  if(length(sel) > 0){ mpheno.d.p.u <- mpheno.d.p[-sel,] } else { mpheno.d.p.u <- mpheno.d.p }
  mpheno.d.r   <- mpheno.d[which(mpheno.d$Sample_Group == "Regressive"),]
  sel <- which(duplicated(mpheno.d.r$Patient))
  if(length(sel) > 0){ mpheno.d.r.u <- mpheno.d.r[-sel,] } else { mpheno.d.r.u <- mpheno.d.r }
  mpheno.d.c   <- mpheno.d[which(mpheno.d$Sample_Group == "Control"),]
  sel <- which(duplicated(mpheno.d.c$Patient))
  if(length(sel) > 0){ mpheno.d.c.u <- mpheno.d.c[-sel,] } else { mpheno.d.c.u <- mpheno.d.c }
  
  mpheno.v.p   <- mpheno.v[which(mpheno.v$Sample_Group == "Progressive"),]
  sel <- which(duplicated(mpheno.v.p$Patient))
  if(length(sel) > 0){ mpheno.v.p.u <- mpheno.v.p[-sel,] } else { mpheno.v.p.u <- mpheno.v.p }
  mpheno.v.r   <- mpheno.v[which(mpheno.v$Sample_Group == "Regressive"),]
  sel <- which(duplicated(mpheno.v.r$Patient))
  if(length(sel) > 0){ mpheno.v.r.u <- mpheno.v.r[-sel,] } else { mpheno.v.r.u <- mpheno.v.r }
  mpheno.v.c   <- mpheno.v[which(mpheno.v$Sample_Group == "Control"),]
  sel <- which(duplicated(mpheno.v.c$Patient))
  if(length(sel) > 0){ mpheno.v.c.u <- mpheno.v.c[-sel,] } else { mpheno.v.c.u <- mpheno.v.c }
  
  df["Patients", "methyl.d.prog"] <- dim(mpheno.d.p.u)[1]
  df["Patients", "methyl.d.reg"]  <- dim(mpheno.d.r.u)[1]
  df["Patients", "methyl.d.cont"]  <- dim(mpheno.d.c.u)[1]
  df["Patients", "methyl.v.prog"] <- dim(mpheno.v.p.u)[1]
  df["Patients", "methyl.v.reg"]  <- dim(mpheno.v.r.u)[1]
  df["Patients", "methyl.v.cont"]  <- dim(mpheno.v.c.u)[1]
  
  df["Lesions profiled", "methyl.d.prog"] <- length(which(mpheno.d$Sample_Group == "Progressive"))
  df["Lesions profiled", "methyl.d.reg"]  <- length(which(mpheno.d$Sample_Group == "Regressive"))
  df["Lesions profiled", "methyl.d.cont"]  <- length(which(mpheno.d$Sample_Group == "Control"))
  df["Lesions profiled", "methyl.v.prog"] <- length(which(mpheno.v$Sample_Group == "Progressive"))
  df["Lesions profiled", "methyl.v.reg"]  <- length(which(mpheno.v$Sample_Group == "Regressive"))
  df["Lesions profiled", "methyl.v.cont"]  <- length(which(mpheno.v$Sample_Group == "Control"))
  
  df["Age at specimen profiled: mean", "methyl.d.prog"] <- mean(mpheno.d.p.u$Age.at.specimen.collected, na.rm=T)
  df["Age at specimen profiled: mean", "methyl.d.reg"]  <- mean(mpheno.d.r.u$Age.at.specimen.collected, na.rm=T)
  df["Age at specimen profiled: mean", "methyl.d.cont"]  <- mean(mpheno.d.c.u$Age.at.specimen.collected, na.rm=T)
  df["Age at specimen profiled: mean", "methyl.v.prog"] <- mean(mpheno.v.p.u$Age.at.specimen.collected, na.rm=T)
  df["Age at specimen profiled: mean", "methyl.v.reg"]  <- mean(mpheno.v.r.u$Age.at.specimen.collected, na.rm=T)
  df["Age at specimen profiled: mean", "methyl.v.cont"]  <- mean(mpheno.v.c.u$Age.at.specimen.collected, na.rm=T)
  
  df["Age at specimen profiled: median", "methyl.d.prog"] <- median(mpheno.d.p.u$Age.at.specimen.collected, na.rm=T)
  df["Age at specimen profiled: median", "methyl.d.reg"]  <- median(mpheno.d.r.u$Age.at.specimen.collected, na.rm=T)
  df["Age at specimen profiled: median", "methyl.d.cont"]  <- median(mpheno.d.c.u$Age.at.specimen.collected, na.rm=T)
  df["Age at specimen profiled: median", "methyl.v.prog"] <- median(mpheno.v.p.u$Age.at.specimen.collected, na.rm=T)
  df["Age at specimen profiled: median", "methyl.v.reg"]  <- median(mpheno.v.r.u$Age.at.specimen.collected, na.rm=T)
  df["Age at specimen profiled: median", "methyl.v.cont"]  <- median(mpheno.v.c.u$Age.at.specimen.collected, na.rm=T)
  
  df["Age at specimen profiled: range", "methyl.d.prog"] <- paste(range(mpheno.d.p.u$Age.at.specimen.collected, na.rm=T),collapse="-")
  df["Age at specimen profiled: range", "methyl.d.reg"]  <- paste(range(mpheno.d.r.u$Age.at.specimen.collected, na.rm=T),collapse="-")
  df["Age at specimen profiled: range", "methyl.d.cont"]  <- paste(range(mpheno.d.c.u$Age.at.specimen.collected, na.rm=T),collapse="-")
  df["Age at specimen profiled: range", "methyl.v.prog"] <- paste(range(mpheno.v.p.u$Age.at.specimen.collected, na.rm=T),collapse="-")
  df["Age at specimen profiled: range", "methyl.v.reg"]  <- paste(range(mpheno.v.r.u$Age.at.specimen.collected, na.rm=T),collapse="-")
  df["Age at specimen profiled: range", "methyl.v.cont"]  <- paste(range(mpheno.v.c.u$Age.at.specimen.collected, na.rm=T),collapse="-")
  
  df["Smoking History: mean", "methyl.d.prog"] <- mean(mpheno.d.p.u$Pack.years, na.rm=T)
  df["Smoking History: mean", "methyl.d.reg"]  <- mean(mpheno.d.r.u$Pack.years, na.rm=T)
  df["Smoking History: mean", "methyl.d.cont"]  <- mean(mpheno.d.c.u$Pack.years, na.rm=T)
  df["Smoking History: mean", "methyl.v.prog"] <- mean(mpheno.v.p.u$Pack.years, na.rm=T)
  df["Smoking History: mean", "methyl.v.reg"]  <- mean(mpheno.v.r.u$Pack.years, na.rm=T)
  df["Smoking History: mean", "methyl.v.cont"]  <- mean(mpheno.v.c.u$Pack.years, na.rm=T)
  
  df["Smoking History: median", "methyl.d.prog"] <- median(mpheno.d.p.u$Pack.years, na.rm=T)
  df["Smoking History: median", "methyl.d.reg"]  <- median(mpheno.d.r.u$Pack.years, na.rm=T)
  df["Smoking History: median", "methyl.d.cont"]  <- median(mpheno.d.c.u$Pack.years, na.rm=T)
  df["Smoking History: median", "methyl.v.prog"] <- median(mpheno.v.p.u$Pack.years, na.rm=T)
  df["Smoking History: median", "methyl.v.reg"]  <- median(mpheno.v.r.u$Pack.years, na.rm=T)
  df["Smoking History: median", "methyl.v.cont"]  <- median(mpheno.v.c.u$Pack.years, na.rm=T)
  
  df["Smoking History: range", "methyl.d.prog"] <- paste(range(mpheno.d.p.u$Pack.years, na.rm=T),collapse="-")
  df["Smoking History: range", "methyl.d.reg"]  <- paste(range(mpheno.d.r.u$Pack.years, na.rm=T),collapse="-")
  df["Smoking History: range", "methyl.d.cont"]  <- paste(range(mpheno.d.c.u$Pack.years, na.rm=T),collapse="-")
  df["Smoking History: range", "methyl.v.prog"] <- paste(range(mpheno.v.p.u$Pack.years, na.rm=T),collapse="-")
  df["Smoking History: range", "methyl.v.reg"]  <- paste(range(mpheno.v.r.u$Pack.years, na.rm=T),collapse="-")
  df["Smoking History: range", "methyl.v.cont"]  <- paste(range(mpheno.v.c.u$Pack.years, na.rm=T),collapse="-")
  
  df["Gender: M", "methyl.d.prog"] <- length(which(mpheno.d.p.u$Gender == "M"))
  df["Gender: M", "methyl.d.reg"]  <- length(which(mpheno.d.r.u$Gender == "M"))
  df["Gender: M", "methyl.d.cont"]  <- length(which(mpheno.d.c.u$Gender == "M"))
  df["Gender: M", "methyl.v.prog"] <- length(which(mpheno.v.p.u$Gender == "M"))
  df["Gender: M", "methyl.v.reg"]  <- length(which(mpheno.v.r.u$Gender == "M"))
  df["Gender: M", "methyl.v.cont"]  <- length(which(mpheno.v.c.u$Gender == "M"))
  df["Gender: F", "methyl.d.prog"] <- length(which(mpheno.d.p.u$Gender == "F"))
  df["Gender: F", "methyl.d.reg"]  <- length(which(mpheno.d.r.u$Gender == "F"))
  df["Gender: F", "methyl.d.cont"]  <- length(which(mpheno.d.c.u$Gender == "F"))
  df["Gender: F", "methyl.v.prog"] <- length(which(mpheno.v.p.u$Gender == "F"))
  df["Gender: F", "methyl.v.reg"]  <- length(which(mpheno.v.r.u$Gender == "F"))
  df["Gender: F", "methyl.v.cont"]  <- length(which(mpheno.v.c.u$Gender == "F"))
  
  df["COPD: Y", "methyl.d.prog"] <- length(which(mpheno.d.p.u$COPD == "YES"))
  df["COPD: Y", "methyl.d.reg"]  <- length(which(mpheno.d.r.u$COPD == "YES"))
  df["COPD: Y", "methyl.d.cont"]  <- length(which(mpheno.d.c.u$COPD == "YES"))
  df["COPD: Y", "methyl.v.prog"] <- length(which(mpheno.v.p.u$COPD == "YES"))
  df["COPD: Y", "methyl.v.reg"]  <- length(which(mpheno.v.r.u$COPD == "YES"))
  df["COPD: Y", "methyl.v.cont"]  <- length(which(mpheno.v.c.u$COPD == "YES"))
  df["COPD: N", "methyl.d.prog"] <- length(which(mpheno.d.p.u$COPD == "NO"))
  df["COPD: N", "methyl.d.reg"]  <- length(which(mpheno.d.r.u$COPD == "NO"))
  df["COPD: N", "methyl.d.cont"]  <- length(which(mpheno.d.c.u$COPD == "NO"))
  df["COPD: N", "methyl.v.prog"] <- length(which(mpheno.v.p.u$COPD == "NO"))
  df["COPD: N", "methyl.v.reg"]  <- length(which(mpheno.v.r.u$COPD == "NO"))
  df["COPD: N", "methyl.v.cont"]  <- length(which(mpheno.v.c.u$COPD == "NO"))
  df["COPD: ?", "methyl.d.prog"] <- length(which(mpheno.d.p.u$COPD == ""))
  df["COPD: ?", "methyl.d.reg"]  <- length(which(mpheno.d.r.u$COPD == ""))
  df["COPD: ?", "methyl.d.cont"]  <- length(which(mpheno.d.c.u$COPD == ""))
  df["COPD: ?", "methyl.v.prog"] <- length(which(mpheno.v.p.u$COPD == ""))
  df["COPD: ?", "methyl.v.reg"]  <- length(which(mpheno.v.r.u$COPD == ""))
  df["COPD: ?", "methyl.v.cont"]  <- length(which(mpheno.v.c.u$COPD == ""))
  
  df["Previous history of lung cancer: Y", "methyl.d.prog"] <- length(which(mpheno.d.p.u$Previous.Lung.CA == 1))
  df["Previous history of lung cancer: Y", "methyl.d.reg"]  <- length(which(mpheno.d.r.u$Previous.Lung.CA == 1))
  df["Previous history of lung cancer: Y", "methyl.d.cont"]  <- length(which(mpheno.d.c.u$Previous.Lung.CA == 1))
  df["Previous history of lung cancer: Y", "methyl.v.prog"] <- length(which(mpheno.v.p.u$Previous.Lung.CA == 1))
  df["Previous history of lung cancer: Y", "methyl.v.reg"]  <- length(which(mpheno.v.r.u$Previous.Lung.CA == 1))
  df["Previous history of lung cancer: Y", "methyl.v.cont"]  <- length(which(mpheno.v.c.u$Previous.Lung.CA == 1))
  df["Previous history of lung cancer: N", "methyl.d.prog"] <- length(which(mpheno.d.p.u$Previous.Lung.CA == 0))
  df["Previous history of lung cancer: N", "methyl.d.reg"]  <- length(which(mpheno.d.r.u$Previous.Lung.CA == 0))
  df["Previous history of lung cancer: N", "methyl.d.cont"]  <- length(which(mpheno.d.c.u$Previous.Lung.CA == 0))
  df["Previous history of lung cancer: N", "methyl.v.prog"] <- length(which(mpheno.v.p.u$Previous.Lung.CA == 0))
  df["Previous history of lung cancer: N", "methyl.v.reg"]  <- length(which(mpheno.v.r.u$Previous.Lung.CA == 0))
  df["Previous history of lung cancer: N", "methyl.v.cont"]  <- length(which(mpheno.v.c.u$Previous.Lung.CA == 0))
  
  
  ################################################################################################
  # WGS
  ################################################################################################
  wgs.pheno.p   <- wgs.pheno[which(wgs.pheno$progression == 1),]
  sel <- which(duplicated(wgs.pheno.p$Patient))
  if(length(sel) > 0){ wgs.pheno.p.u <- wgs.pheno.p[-sel,] } else { wgs.pheno.p.u <- wgs.pheno.p }
  wgs.pheno.r   <- wgs.pheno[which(wgs.pheno$progression == 0),]
  sel <- which(duplicated(wgs.pheno.r$Patient))
  if(length(sel) > 0){ wgs.pheno.r.u <- wgs.pheno.r[-sel,] } else { wgs.pheno.r.u <- wgs.pheno.r }
  
  df["Patients", "seq.prog"] <- dim(wgs.pheno.p.u)[1]
  df["Patients", "seq.reg"]  <- dim(wgs.pheno.r.u)[1]
  
  df["Lesions profiled", "seq.prog"] <- length(which(wgs.pheno.p$progression == 1))
  df["Lesions profiled", "seq.reg"]  <- length(which(wgs.pheno.r$progression == 0))
  
  df["Age at specimen profiled: mean", "seq.prog"] <- mean(wgs.pheno.p.u$Age.at.specimen.profiled, na.rm=T)
  df["Age at specimen profiled: mean", "seq.reg"]  <- mean(wgs.pheno.r.u$Age.at.specimen.profiled, na.rm=T)
  
  df["Age at specimen profiled: median", "seq.prog"] <- median(wgs.pheno.p.u$Age.at.specimen.profiled, na.rm=T)
  df["Age at specimen profiled: median", "seq.reg"]  <- median(wgs.pheno.r.u$Age.at.specimen.profiled, na.rm=T)
  
  df["Age at specimen profiled: range", "seq.prog"] <- paste(range(wgs.pheno.p.u$Age.at.specimen.profiled, na.rm=T),collapse="-")
  df["Age at specimen profiled: range", "seq.reg"]  <- paste(range(wgs.pheno.r.u$Age.at.specimen.profiled, na.rm=T),collapse="-")
  
  df["Smoking History: mean", "seq.prog"] <- mean(wgs.pheno.p.u$Pack.years, na.rm=T)
  df["Smoking History: mean", "seq.reg"]  <- mean(wgs.pheno.r.u$Pack.years, na.rm=T)
  
  df["Smoking History: median", "seq.prog"] <- median(wgs.pheno.p.u$Pack.years, na.rm=T)
  df["Smoking History: median", "seq.reg"]  <- median(wgs.pheno.r.u$Pack.years, na.rm=T)
  
  df["Smoking History: range", "seq.prog"] <- paste(range(wgs.pheno.p.u$Pack.years, na.rm=T),collapse="-")
  df["Smoking History: range", "seq.reg"]  <- paste(range(wgs.pheno.r.u$Pack.years, na.rm=T),collapse="-")
  
  df["Gender: M", "seq.prog"] <- length(which(wgs.pheno.p.u$Gender == "M"))
  df["Gender: M", "seq.reg"]  <- length(which(wgs.pheno.r.u$Gender == "M"))
  df["Gender: F", "seq.prog"] <- length(which(wgs.pheno.p.u$Gender == "F"))
  df["Gender: F", "seq.reg"]  <- length(which(wgs.pheno.r.u$Gender == "F"))
  
  df["COPD: Y", "seq.prog"] <- length(which(wgs.pheno.p.u$COPD == "YES"))
  df["COPD: Y", "seq.reg"]  <- length(which(wgs.pheno.r.u$COPD == "YES"))
  df["COPD: N", "seq.prog"] <- length(which(wgs.pheno.p.u$COPD == "NO"))
  df["COPD: N", "seq.reg"]  <- length(which(wgs.pheno.r.u$COPD == "NO"))
  df["COPD: ?", "seq.prog"] <- length(which(wgs.pheno.p.u$COPD == ""))
  df["COPD: ?", "seq.reg"]  <- length(which(wgs.pheno.r.u$COPD == ""))
  
  df["Previous history of lung cancer: Y", "seq.prog"] <- length(which(wgs.pheno.p.u$Prev.History.of.LC == 1))
  df["Previous history of lung cancer: Y", "seq.reg"]  <- length(which(wgs.pheno.r.u$Prev.History.of.LC == 1))
  df["Previous history of lung cancer: N", "seq.prog"] <- length(which(wgs.pheno.p.u$Prev.History.of.LC == 0))
  df["Previous history of lung cancer: N", "seq.reg"]  <- length(which(wgs.pheno.r.u$Prev.History.of.LC == 0))
  
  library(WriteXLS)
  WriteXLS('df', ExcelFileName=filename, row.names=T)
  
}

