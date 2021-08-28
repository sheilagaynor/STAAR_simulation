library(STAAR); library(rje); library(matrixStats) 
library(gdsfmt); library(SeqArray); library(STAAR)
library(SeqVarTools); library(dplyr); library(data.table)

genofile <- seqOpen('../EU_20K_gds/sim_geno_EU_20K.gds')
all_snp_pos <- as.numeric(seqGetData(genofile,"position"))
all_var_id <- seqGetData(genofile, "variant.id")

full_agg_list <- fread(file="../EU_20K_gds/sim_agg_units.csv.gz")
set.seed(1001)
selected_set_name <- unique(full_agg_list$group_id)[sample(1:length(unique(full_agg_list$group_id)),1)]
agg_select <- full_agg_list[full_agg_list$group_id==selected_set_name,]
variant.id.gene <- all_var_id[all_snp_pos %in% agg_select$pos]
seqSetFilter(genofile,variant.id=variant.id.gene)

# simulate causal variants using annotations
delta0 <- rje::logit(0.18)
delta15 <- rep(log(5),10)
num_snp <- length(all_snp_pos)
set.seed(2021)
annot_all <- matrix( rnorm(num_snp*10), num_snp, 10)
annot <- annot_all[which(all_snp_pos %in% agg_select$pos),]
causalprob <- apply(annot,1,function(annot){
      set.seed(43)
      ind <- sample(1:10,5)
      rje::expit(delta0 + delta15[ind] %*% annot[ind])})
causal_ind <- rbinom(dim(annot)[1],1,causalprob)

# Simulate phenotype data
num_samp <- 20000
set.seed(28)
X1 <- rnorm(num_samp)
X2 <- rbinom(num_samp,1,0.5)
eps <- rnorm(num_samp)
maf <- 1-seqAlleleFreq(genofile)
beta <- -0.1 * log10(maf[which(causal_ind==1)])
beta <- beta * (2 * (runif(length(beta)) < 0.8) - 1)
genotype <- seqGetData(genofile, "$dosage")
Y <- 0.5 * X1 + 0.5 * X2 + genotype[,which(causal_ind==1)] %*% beta + eps
ID <- seqGetData(genofile,"sample.id")
pheno <- data.frame(ID,Y,X1,X2)
write.csv(pheno, '../EU_20K_gds/sim_pheno.csv', row.names = F)
