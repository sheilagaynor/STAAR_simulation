# Get in bim bed fam from tped/fam
for (arrayid in 1:27){
  command_text <- paste0("module load plink/1.90-fasrc01; plink --tfile ../EU_20K_gds/geno_EU_20K_", arrayid, 
                         " --make-bed --out ../EU_20K_gds/geno_EU_20K_",arrayid)
  system(command_text)  
}


# Merge to a single plink file
command_mergelist <- paste0('seq 2 27 | xargs -I _num echo ../EU_20K_gds/geno_EU_20K__num > ../EU_20K_gds/mergelist.txt')
system(command_mergelist)
command_merge <- paste0('module load plink/1.90-fasrc01; plink --bfile ../EU_20K_gds/geno_EU_20K_1 --merge-list ../EU_20K_gds/mergelist.txt --make-bed --out ../EU_20K_gds/sim_geno_EU_20K')
system(command_merge)


# Get into gds
library(SeqArray)
seqBED2GDS('../EU_20K_gds/sim_geno_EU_20K.bed',
           '../EU_20K_gds/sim_geno_EU_20K.fam',
           '../EU_20K_gds/sim_geno_EU_20K.bim',
           '../EU_20K_gds/sim_geno_EU_20K.gds')


# Fix scientific notation messed up positions
library(gdsfmt); library(SeqArray); library(STAAR)
library(SeqVarTools); library(dplyr); library(matrixStats)
genofile <- seqOpen('../EU_20K_gds/sim_geno_EU_20K.gds', readonly = FALSE)
all_snp_pos <- as.numeric(seqGetData(genofile,"position"))
all_rs_id <- seqGetData(genofile, "annotation/id")
select_pos <- which(grepl( 'e', all_rs_id, fixed = TRUE))
all_snp_pos[select_pos] <- as.integer(as.numeric(substr(all_rs_id[select_pos], 3, nchar(all_rs_id[select_pos]))) + 
                                        all_snp_pos[select_pos])
seqAddValue(genofile, "position", all_snp_pos, replace=TRUE)
seqClose(genofile)


# Get the aggregation units
genofile <- seqOpen('../EU_20K_gds/sim_geno_EU_20K.gds', readonly = FALSE)
all_var_id <- seqGetData(genofile, "variant.id")
all_snp_pos <- as.numeric(seqGetData(genofile,"position"))
sigLength <- 5000; agg_list <- c()
for ( i in 1:1200 ){
  # Randomly select a signal region
  set.seed(i)
  startloc <- sample(1:(max(all_snp_pos)-sigLength+1),1)
  endloc <- startloc + sigLength - 1
  # Extract the RVs in signal region 
  var_region <- all_var_id[all_snp_pos >= startloc & all_snp_pos <= endloc]
  seqSetFilter(genofile,variant.id=var_region)
  set_var_pos <- as.numeric(seqGetData(genofile,"position"))
  rare_var_region <- (seqAlleleFreq(genofile) <= 0.01) | (seqAlleleFreq(genofile) >= 0.99)
  # Add to agg list
  agg_list <- rbind(agg_list, cbind(set_var_pos[rare_var_region],paste0('simgene_',i)))
  seqResetFilter(genofile)
}
full_agg_list <- data.frame(cbind(2,agg_list[,1],'C','G',agg_list[,2]), stringsAsFactors = F)
names(full_agg_list) <- c('chr', 'pos', 'ref', 'alt', 'group_id')
full_agg_list$pos <- as.integer(as.character(full_agg_list$pos))
write.csv(full_agg_list, file=gzfile("../EU_20K_gds/sim_agg_units.csv.gz"), row.names = F)


#Add annotation and additional channels
num_snp <- length(all_snp_pos)
#seqAddValue(genofile, "annotation/filter", rep(as.character('PASS'),num_snp), replace=TRUE)
seqAddValue(genofile, "annotation/info/AVGDP", rep(25,num_snp), replace=TRUE)
set.seed(2021)
annot <- matrix( rnorm(num_snp*10), num_snp, 10)
annot_rank <- colRanks(annot,preserveShape = TRUE)
PHRED <- -10*log10(1-annot_rank/dim(annot_rank)[1])
addfolder.gdsn(index.gdsn(genofile, "annotation/info"), "sim_annotation")
for (ii in 1:10){
  seqAddValue(genofile, paste0("annotation/info/sim_annotation/annot_",ii), PHRED[,ii], replace=TRUE)
}

seqClose(genofile)