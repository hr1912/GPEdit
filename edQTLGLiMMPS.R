
library('magrittr')
library('lme4')
library('parallel')
#library('microbenchmark')
library('compiler')

options(warn=-1)

numCores <- detectCores()
# using all availiable cores
# numCores

###############
#  GLiMMPS method orginated from 
# generalized linear mixed model regression  with LRT test for fixed effect 

##################

glmm.sQTL <- function(data) {
  #######  glmm 
  # library(lme4)
  #options(warn=2) #this turns all warnings into errors, KZ change
  
  snp.pval <- 1
  betas <-  c(0,0)
  
  if ( length(unique(data$SNP[!is.na(data$y) & data$n>0 & !is.na(data$SNP)] )) >1) {
    nomissing <- (!is.na(data$y) & data$n>0 & !is.na(data$SNP))
    response <- cbind(data$y[nomissing], data$n[nomissing]-data$y[nomissing])
    SNP <- data$SNP[nomissing]
    obs <- seq(1,length(SNP))
    
    testglm = try( suppressMessages( glmer(response ~ SNP +(1|obs), family=binomial, nAGQ = 0)),silent=TRUE)  #invidual-level random effect for overdispersion. 
    testglm0 = try(suppressMessages( glmer(response ~ 1 +(1|obs), family=binomial, nAGQ = 0)),silent=TRUE) # 
    
    #  if (!( inherits(testglm,"try-error") || inherits(testglm0,"try-error") )) { # only when it converges ## KZchange
    betas <- try(fixef(testglm), silent=TRUE)
    if (! inherits(betas,"try-error") & !is.na(betas[2])) {
      testanova <- try(suppressMessages(anova(testglm, testglm0)),silent=TRUE)
      if (! inherits(testanova,"try-error") & length(testanova$"Pr(>Chisq)"[2]) > 0) {
        snp.pval <- testanova$"Pr(>Chisq)"[2]
      }
      #snp.pval <-  anova(testglm, testglm0)$"Pr(>Chisq)"[2]
    } else betas=c(0, 0)
    
    #  } ## KZchange
    #psi.fitted <- results$fitted.value
  } 
  #  options(error = NULL, warn = 0) ## KZchange
  
  return ( list(betas =betas,  pval = snp.pval ) )  
}


glmm.sQTL.cmp <- cmpfun(glmm.sQTL)

args <- commandArgs(TRUE)
ct <- args[1]

# cancer type
#ct <- "CHOL"

rootdir <- '/extraspace/Projects/edQTL'

CTDIR <- file.path(rootdir, "05.edQTL", ct)
# 05.edQTL is the work directory for edQTL calling 

if (!file.exists( CTDIR )) {system(paste ("mkdir",CTDIR)) }

GENODIR <- file.path(CTDIR, 'geno')
PHIDIR<- file.path(CTDIR, 'phi')

if (!file.exists( GENODIR )) {system(paste ("mkdir",GENODIR)) }
if (!file.exists( PHIDIR )) {system(paste ("mkdir",PHIDIR)) }

setwd(PHIDIR)
system(paste0("ln -s ", file.path(rootdir, "04.phi", ct), "/*.rds ./"))
# 04.phi is the directory where all editing frequencies are located 

setwd(GENODIR)
system(paste0("ln -s ", file.path(rootdir, "01.imputed", paste0(ct, "_impute_parsed.txt")),
              " ./"))
system(paste0("ln -s ", file.path(rootdir, "00.data", paste0("dbsnp.hg38.",ct,".txt")),
              " ./"))

# 00.data is the directory where all the dbsnp locataion are located
# 01.imputed is the directory where all the imputed genotype data are located

setwd(CTDIR)


#editing_idx <- 1:6215

cis_distance <- 200000 #  within 200kb defined as cis regulatory 
pval_cutoff <- .05 # a rough pvalue cutoff 
valid_chr <- seq(1:22) # autosomes

editing_info <- readr::read_rds(file.path(PHIDIR,paste0(ct,"_edInfo.rds"))) %>% 
  dplyr::filter(chr %in% paste0("chr", valid_chr)) %>% 
  dplyr::select(chr, pos, mut, edID) %>% dplyr::distinct_all()

editing_ids <- editing_info$edID ### editing information with in each cancer type

###

dbsnp <- readr::read_tsv(file.path(GENODIR, paste0("dbsnp.hg38.",ct,".txt")))

geno_info <- readr::read_tsv(file.path(GENODIR, paste0(ct, "_impute_parsed.txt")))
colnames(geno_info) <- paste0(colnames(geno_info), "A")

geno_info %>% dplyr::rename(ID = IDA) -> geno_info
geno_info %>% dplyr::select(-ID) %>% as.matrix() -> geno_info_mat

mclapply(1:nrow(geno_info_mat) ,
         function(xx) min(table(geno_info_mat[xx,])) > 3, 
         mc.cores = numCores) -> geno_info_valid_idx 

geno_info[unlist(geno_info_valid_idx), ] -> geno_info

dbsnp %>% dplyr::right_join(geno_info) %>% 
  dplyr::rename(SNP = ID, CHROM = `#CHROM`)  %>% 
  dplyr::arrange(CHROM, POS) -> combined_tbl

combined_tbl %>% dplyr::group_by(CHROM) -> combined_tbl_gp
dplyr::group_split(combined_tbl_gp) -> combined_tbl_split


all_reads_mat <- readr::read_rds(file.path(PHIDIR,paste0(ct,"_allreadsCount.rds"))) ### coverage reads number
ed_reads_mat <- readr::read_rds(file.path(PHIDIR,paste0(ct,"_edreadsCount.rds"))) ### edited reads number 

n_samples <- dim(all_reads_mat)[2]
n_ed <- dim(all_reads_mat)[1]

print(n_samples)
print(n_ed)

all_reads_mat %>% colnames() -> IDs.pheno


TMPASSODIR <- file.path(CTDIR,"tmpasso")
if (!file.exists( TMPASSODIR )) {system(paste ("mkdir",TMPASSODIR)) }

for (c in 1:22) {
  cmmd <- paste("mkdir ",TMPASSODIR,"/chr",c,sep="")
  if (!file.exists( paste(TMPASSODIR,"/chr",c,sep="")) ) {system(cmmd)}
}


# test 
#editing_idx <- 3

### function for calling one edQTL

one_ed_run <- function (editing_idx) { 
  
  targetedID <- editing_ids[editing_idx]
  
  editing_info[editing_idx,] -> targeted_info
  
  chr <- targeted_info$chr %>% sub("chr", "", .)
  
  SNP_startpos <- targeted_info$pos - cis_distance
  SNP_endpos <- targeted_info$pos + cis_distance
  
  
  #print(targeted_info)
  
  # get one tibble
  
  combined_tbl_split[[as.numeric(chr)]] %>% 
    dplyr::filter(POS >= SNP_startpos & POS <= SNP_endpos) -> one_tmp_tbl
  
  nsnps <- dim(one_tmp_tbl)[1]
  
  if (nsnps == 0) {return(0)}
  
  IDs.geno <- colnames(one_tmp_tbl)[-c(1:5)] 
  IDs.common <- intersect(IDs.geno, IDs.pheno)
  
  sub <-  match(IDs.common, IDs.pheno)
  
  n <- all_reads_mat[targetedID,sub]
  y <- ed_reads_mat[targetedID,sub]
  

  geno <- one_tmp_tbl %>% dplyr::select(-c(1:5)) %>% as.matrix() 
  rownames(geno) <- one_tmp_tbl$SNP

  geno <- geno[, IDs.common]
  #colnames(geno) <- names(n)
  
  pvals.glmm <- rep(NA,nsnps)
  betas.glmm <-  rep(NA,nsnps) 

  #microbenchmark(apply(geno, 1, glmm.test), times = 10L)
  # 
  
  
   # glmm.test <- function(SNP) {
   #   results.glmm<-glmm.sQTL(list(n=n,y=y, SNP=SNP))
   #   return(c(results.glmm$pval, results.glmm$betas[2]))
   # }
   # 
   # glmm.test.cmp <- cmpfun(glmm.test)
   # microbenchmark(apply(geno, 1, glmm.test.cmp), times = 10L)

  if (nsnps == 1) {

    results.glmm <- glmm.sQTL.cmp(list(n=n, y=y, SNP= geno))
    pvals.glmm <- results.glmm$pval
    betas.glmm <- results.glmm$betas[2]

  } else {
  
    for (gi in 1:nsnps) {
      
     # onedata <- list(n=n, y=y, SNP= geno[gi,])
      results.glmm <- glmm.sQTL.cmp(list(n=n, y=y, SNP= geno[gi,]))
      pvals.glmm[gi] <- results.glmm$pval
      betas.glmm[gi] <- results.glmm$betas[2] 
      
    }
  }  
  
  if (any(pvals.glmm <= pval_cutoff)) {
  
    if (nsnps > 1) {
      
      tmpout <- cbind(one_tmp_tbl$CHROM, one_tmp_tbl$SNP, one_tmp_tbl$POS,
                      formatC(pvals.glmm,format="e",digits=3),  
                      round(betas.glmm,3), rep(targetedID, nsnps))
      
      colnames(tmpout) <- c("Chr","SNPID","Pos",
                            "pvals.glmm",
                            "Beta.glmm",
                            "edID")
      
      write.table(tmpout  , paste(TMPASSODIR,"/chr",chr,"/",targetedID,".asso",sep=""),
                  row.names=F,col.names=T, quote=F,sep="\t" )
    
    } else {
      
      tmpout <- cbind(one_tmp_tbl$CHROM, one_tmp_tbl$SNP, one_tmp_tbl$POS,
                      formatC(pvals.glmm,format="e",digits=3),  
                      round(betas.glmm,3), rep(targetedID, nsnps))
      
      names(tmpout) <- c("Chr","SNPID","Pos",
                         "pvals.glmm",
                         "Beta.glmm", 
                         "edID")
      
      write.table(t(as.data.frame(tmpout)), paste(TMPASSODIR,"/chr",chr,"/",targetedID,".asso",sep=""),
                  row.names=F,col.names=T, quote=F,sep="\t")
      
    }
  }
  
  #return(1)
}
#TMPGENODIR <- file.path(CTDIR, "tmpgeno")

#system.time(one_ed_run(3))

mclapply(1:length(editing_ids), one_ed_run, mc.cores = numCores)
