# In 27 parallelized jobs, merges genotype ouput in R
# Saves output as tped/tfam file for use in Plink

# Load libraries
arrayid <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) 
chunks <- split(1:243, ceiling(seq_along(1:243)/9))
library(data.table)
library(gdsfmt); library(SeqArray); library(STAAR)
library(SeqVarTools); library(dplyr); library(matrixStats)

iter_kk <- chunks[[arrayid]]
print(iter_kk)
geno_in <- load(paste0("../EU_20K/geno_EU_20K_",iter_kk[1],".Rdata"))
pos <- snploc$CHROM_POS + seq(0,1000000*243,1000000)[iter_kk[1]]
pos_table <- data.table(cbind(2, paste0('rs',pos), 0, pos))
all_geno <- cbind(pos_table, t(genotype_all))
all_maf <- maf

for ( kk in iter_kk[2:length(iter_kk)] ){
  load(paste0("../EU_20K/geno_EU_20K_",kk,".Rdata"))
  # Update position to get continuous chrom 
  pos <- snploc$CHROM_POS + seq(0,1000000*243,1000000)[kk]
  pos_table <- data.table(cbind(2, paste0('rs',pos), 0, pos))
  # Add to output
  all_geno <- rbind(all_geno, cbind(pos_table, t(genotype_all)))
  all_maf <- c(all_maf, maf)
}
rm(genotype_all); rm(snploc); rm(maf); gc()

# Generate tfam and tped
tfam <- data.frame( FID = seq(20000), IID = seq(20000), FatherID = 0,
                        MotherID = 0, Sex = 1, Phenotype = 1)

# Output tped, tfam
fwrite(all_geno, paste0("../EU_20K_gds/geno_EU_20K_",arrayid,".tped"), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
write.table(tfam, paste0("../EU_20K_gds/geno_EU_20K_",arrayid,".tfam"), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

# Get in bim bed fam from tped/fam
command_text <- paste0("module load plink/1.90-fasrc01; plink --tfile ../EU_20K_gds/geno_EU_20K_", arrayid, 
                       " --make-bed --out ../EU_20K_gds/geno_EU_20K_",arrayid)
system(command_text)
