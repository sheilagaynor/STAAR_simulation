arrayid <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) 
chunks <- split(1:243, ceiling(seq_along(1:243)/12))
chunks[[20]] <- c(chunks[[20]],chunks[[21]]); chunks <- chunks[1:20]
library(data.table)

recode_to_alleles <- function( geno_needs_recoding ){
  geno_needs_recoding <- geno_needs_recoding[, names(geno_needs_recoding) := lapply(.SD, as.character)]
  geno_needs_recoding <- geno_needs_recoding[, replace(.SD, .SD == '0', 'C C')]
  geno_needs_recoding <- geno_needs_recoding[, replace(.SD, .SD == '1', 'C G')]
  geno_needs_recoding <- geno_needs_recoding[, replace(.SD, .SD == '2', 'G G')]
  return(geno_needs_recoding)
}

for(kk in chunks[[arrayid]]){
	set.seed(kk)
	command_text <- paste0("cd ../cosi_1.2/examples/bestfit; perl run_EU_20K_",arrayid,".pl")
	system(command_text)

	##Read in output
	snploc <- read.table(paste0("../cosi_1.2/examples/bestfit/geno_EU_20K_",arrayid,".pos-1"),header=T)
	genotype_input <- read.table(paste0("../cosi_1.2/examples/bestfit/geno_EU_20K_",arrayid,".hap-1"),header=F)
	
	##Convert format
	genotype_input <- genotype_input[,3:dim(genotype_input)[2]] #drop first 2 cols
	min_allele <- ifelse(snploc$FREQ1 > snploc$FREQ2, 2, 1) #getting minor allele
	genotype_input <- rbind(min_allele,genotype_input)
	genotype_TF <- apply(genotype_input,2, function(z) z[2:length(z)]==z[1]) #get counts of minor alleles
	genotype_dosage <- data.table(apply(genotype_TF,2, function(z) z[1:(length(z)/2)]+z[(length(z)/2+1):length(z)])) #get dosage
	rm(genotype_input); rm(genotype_TF)
	maf <- apply(genotype_dosage,2,mean)/2
	gc()

	##Save genotypes with MAF>0
	genotype_dosage <- genotype_dosage[,maf>0,with=F]
	genotype_all <- recode_to_alleles(genotype_dosage)
	snploc <- snploc[maf>0,];	maf <- maf[maf>0]
	rm(genotype_dosage);
	save(genotype_all, snploc, maf, file=paste0("../EU_20K/geno_EU_20K_",kk,".Rdata"))
}