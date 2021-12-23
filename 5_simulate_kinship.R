# Generates kinship matrix
# Ref: https://github.com/xihaoli/STAAR/blob/master/docs/STAAR_vignette.html

library(kinship2)
grid <- 1
Npercell <- 10000
ndiv <- 1
vfam <- 0.5
N <- round(grid*grid*Npercell/ndiv)

unitmat <- matrix(0.5, 4, 4)
diag(unitmat) <- 1
unitmat[1,2] <- unitmat[2,1] <- 0
ped <- data.frame(famid = rep(as.integer(1:2500), each=4), id = as.integer(1:20000L), fa = rep(0, 20000), mo = rep(0, 20000))
for(i in 1:2500) {
  ped$fa[4*i-(0:1)] <- ped$id[4*i-3]
  ped$mo[4*i-(0:1)] <- ped$id[4*i-2]
}
ped$famid[10001:20000] <- as.integer(2501:12500)
kins <- makekinship(ped$famid, ped$id, ped$fa, ped$mo)
saveRDS(kins,'../EU_20K_gds/sim_kinship_sparse.Rds')
saveRDS(as.matrix(kins),'../EU_20K_gds/sim_kinship_dense.Rds')