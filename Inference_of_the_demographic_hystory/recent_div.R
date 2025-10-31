# Script for generate reference table of recent_div scenario
###module load statistics/R/4.3.0
###module load devel/python/Python-3.7.9

library(reticulate)
library(poolfstat)
library(KScorrect)
library(Rcpp)
#library(progress)
#library(doParallel)
use_python("/tools/devel/python/Python-3.7.9/bin/python")
source_python(file = "simu.py")
source_python(file = "compute_stat.py")

args <- commandArgs(trailingOnly = TRUE)

k=as.integer(args[1])
j=k+as.numeric(args[2])-1
for (seed in k:j) {
	nsample=1 ; 
	nchr=1000L ; chr_length=1e4 #1,000 chromosome of 10,000 pb
	lambda.cov=50 ; eps=0.001 #parameter for simulate pool data
	maf.thr=0.00 ; min.rc=2
	exp.eps=0.01 ; nthreads=1
	set.seed(seed)

	run.simu_recent_div<-function(Ne_ancpp,Ne_fu_ancfound,Ne_wp_ancfound,Ne_sp_found,Ne_wp_found,Ne_fu_anc,Ne_wp_anc,
               			 Ne_bot_fu,Ne_bot_sp,Ne_bot_wp,Ne_sp,Ne_wp,Ne_fu,tsplit_PP,tsplit_WP,t_bot,m_ancpp,
                		 m_spwp_anc,m_spwp_rec,m_spfu_anc,m_spfu_rec,m_wpfu_anc,m_wpfu_rec,
						 n_sp,n_wp,n_fu,nchr,chr_length,min.rc,maf.thr,seed,n_cpu,seed_dr){    

    # First run simulaton for have count data
	G = simu_recent_div(Ne_ancpp=Ne_ancpp,Ne_fu=Ne_fu,Ne_wp=Ne_wp,Ne_sp=Ne_sp,Ne_sp_found=Ne_sp_found,Ne_wp_found=Ne_wp_found,Ne_fu_ancfound=Ne_fu_ancfound,Ne_fu_anc=Ne_fu_anc,Ne_wp_anc=Ne_wp_anc,
	               Ne_bot_sp=Ne_bot_sp,Ne_bot_wp=Ne_bot_wp,Ne_bot_fu=Ne_bot_fu,Ne_wp_ancfound=Ne_wp_ancfound,m_spfu_rec=m_spfu_rec,m_wpfu_rec=m_wpfu_rec,m_ancpp=m_ancpp,
	               tsplit_PP=tsplit_PP,tsplit_WP=tsplit_WP,t_bot=t_bot,m_spwp_anc=m_spwp_anc,m_spwp_rec=m_spwp_rec,m_spfu_anc=m_spfu_anc,m_wpfu_anc=m_wpfu_anc,
		       n_sp=n_sp,n_wp=n_wp,n_fu=n_fu,nchr=nchr,chr_length=chr_length,seed=seed,n_cpu=nthreads)
	  
    # Remove non bialleliq sites
    G <- subset(G, !rowSums(G == 2 | G == 3) > 0)

    # Index populations colums
	  sp.idx=1:(2*n_sp)
	  wp.idx=(2*n_sp+1):(2*(n_sp+n_wp))
	  fu.idx=(2*(n_sp+n_wp)+1):(2*(n_sp+n_wp+n_fu))

    # Create a countdata object for poolfstats
	counts=new("countdata")
    counts@nsnp=nrow(G) ;  counts@npops=3
	counts@popnames=c("LSP","LWP","FU")
	counts@snp.info=data.frame(Chromosome=sample(1:nchr,nrow(G),replace = TRUE),
								 Position=sample(1:chr_length,nrow(G),replace = TRUE),
								 RefAllele=rep("poly",nrow(G)),
								 AltAllele=rep(NA,nrow(G)))
	counts@refallele.count=cbind(rowSums(G[,sp.idx]==0),rowSums(G[,wp.idx]==0),rowSums(G[,fu.idx]==0))
	counts@total.count=matrix(rep(c(2*n_sp,2*n_wp,2*n_fu),each=counts@nsnp),ncol=3)

    # Simulate pool data
	pool=sim.readcounts(counts,lambda.cov = rep(lambda.cov,counts@npops),overdisp=1,seq.eps = eps,exp.eps = exp.eps,
                                                   maf.thr=maf.thr,min.rc=min.rc,genome.size=(chr_length*nchr))
	pooldata=pooldata.subset(pool, min.cov.per.pool=10)

    # Summary statistic on simulate pool data
	pool.fstats=compute.fstats(pooldata,verbose=FALSE)
	fst.id=pool.fstats@fst.values$Estimate
	heteros=pool.fstats@heterozygosities$Estimate
	f2=pool.fstats@f2.values$Estimate
	f3=pool.fstats@f3.values$Estimate
	fst.wc=computeFST(pooldata)$Fst[[1]]
	  
	maf.SP=0.5-abs(0.5-(pooldata@refallele.readcount[,1]/pooldata@readcoverage[,1]))
	maf.WP=0.5-abs(0.5-(pooldata@refallele.readcount[,2]/pooldata@readcoverage[,2]))
	maf.FU=0.5-abs(0.5-(pooldata@refallele.readcount[,3]/pooldata@readcoverage[,3]))
	pfix=c(mean(maf.SP==0),mean(maf.WP==0),mean(maf.FU==0))
	  
	alt.SP=0.5-(0.5-(pooldata@refallele.readcount[,1]/pooldata@readcoverage[,1]))
    alt.WP=0.5-(0.5-(pooldata@refallele.readcount[,2]/pooldata@readcoverage[,2]))
    alt.FU=0.5-(0.5-(pooldata@refallele.readcount[,3]/pooldata@readcoverage[,3]))
    pfix.alt=c(mean(alt.SP==0),mean(alt.WP==0),mean(alt.FU==0))
    pfix.ref=c(mean(alt.SP==1),mean(alt.WP==1),mean(alt.FU==1))	  

    # Simulate individual data for summary statistics on individual data
	simu_recent_div_ind(Ne_ancpp=Ne_ancpp,Ne_fu=Ne_fu,Ne_wp=Ne_wp,Ne_sp=Ne_sp,Ne_sp_found=Ne_sp_found,Ne_wp_found=Ne_wp_found,Ne_fu_ancfound=Ne_fu_ancfound,Ne_fu_anc=Ne_fu_anc,Ne_wp_anc=Ne_wp_anc,
                        Ne_bot_sp=Ne_bot_sp,Ne_bot_wp=Ne_bot_wp,Ne_bot_fu=Ne_bot_fu,Ne_wp_ancfound=Ne_wp_ancfound,m_spfu_rec=m_spfu_rec,m_wpfu_rec=m_wpfu_rec,m_ancpp=m_ancpp,
                        tsplit_PP=tsplit_PP,tsplit_WP=tsplit_WP,t_bot=t_bot,m_spwp_anc=m_spwp_anc,m_spwp_rec=m_spwp_rec,m_spfu_anc=m_spfu_anc,m_wpfu_anc=m_wpfu_anc,
                        n_sp=25L,n_wp=18L,n_fu=10L,nchr=nchr,chr_length=chr_length,seed=seed,seed_dr=seed_dr)
    
    # Simulate haplotypes data for summary statistics describe in Navascués et al. (2014)  
    H <- simu_recent_div_Uhl(Ne_ancpp=Ne_ancpp,Ne_fu=Ne_fu,Ne_wp=Ne_wp,Ne_sp=Ne_sp,Ne_sp_found=Ne_sp_found,Ne_wp_found=Ne_wp_found,Ne_fu_ancfound=Ne_fu_ancfound,Ne_fu_anc=Ne_fu_anc,Ne_wp_anc=Ne_wp_anc,
                       Ne_bot_sp=Ne_bot_sp,Ne_bot_wp=Ne_bot_wp,Ne_bot_fu=Ne_bot_fu,Ne_wp_ancfound=Ne_wp_ancfound,m_spfu_rec=m_spfu_rec,m_wpfu_rec=m_wpfu_rec,m_ancpp=m_ancpp,
                       tsplit_PP=tsplit_PP,tsplit_WP=tsplit_WP,t_bot=t_bot,m_spwp_anc=m_spwp_anc,m_spwp_rec=m_spwp_rec,m_spfu_anc=m_spfu_anc,m_wpfu_anc=m_wpfu_anc,
                       n_sp=25L,n_wp=18L,n_fu=10L,nchr=nchr,chr_length=chr_length,seed=seed)

	if ("SP" %in% names(H)) {
  		sp_haplotype_matrix = H[["SP"]]
	} else {
  		sp_haplotype_matrix = NULL
	}

	if ("WP" %in% names(H)) {
  		wp_haplotype_matrix = H[["WP"]]
	} else {
  		wp_haplotype_matrix = NULL
	}

	if ("FU" %in% names(H)) {
  		fu_haplotype_matrix = H[["FU"]]
	} else {
  		fu_haplotype_matrix = NULL
	}

	is_non_empty_matrix <- function(mat) {
  		return(!is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0)
	}

    # Remove non biallelic sites
	if (is_non_empty_matrix(sp_haplotype_matrix)) {
  		valid_sp_rows = rowSums(sp_haplotype_matrix == 2 | sp_haplotype_matrix == 3) == 0
	} else {
  		valid_sp_rows = logical(0)  
	}

	if (is_non_empty_matrix(wp_haplotype_matrix)) {
  		valid_wp_rows = rowSums(wp_haplotype_matrix == 2 | wp_haplotype_matrix == 3) == 0
	} else {
  		valid_wp_rows = logical(0) 
	}

	if (is_non_empty_matrix(fu_haplotype_matrix)) {
  		valid_fu_rows = rowSums(fu_haplotype_matrix == 2 | fu_haplotype_matrix == 3) == 0
	} else {
  		valid_fu_rows = logical(0)  
	}

	valid_rows = valid_sp_rows & valid_wp_rows & valid_fu_rows

	if (is_non_empty_matrix(sp_haplotype_matrix)) {
  		sp_haplotype_matrix <- sp_haplotype_matrix[valid_rows, , drop = FALSE]
	}

	if (is_non_empty_matrix(wp_haplotype_matrix)) {
  		wp_haplotype_matrix <- wp_haplotype_matrix[valid_rows, , drop = FALSE]
	}

	if (is_non_empty_matrix(fu_haplotype_matrix)) {
  		fu_haplotype_matrix <- fu_haplotype_matrix[valid_rows, , drop = FALSE]
	}

	cat("SP Genotype Matrix (Filtered):", dim(sp_haplotype_matrix), "\n")
    cat("WP Genotype Matrix (Filtered):", dim(wp_haplotype_matrix), "\n")
    cat("FU Genotype Matrix (Filtered):", dim(fu_haplotype_matrix), "\n")

	is_non_empty_matrix <- function(mat) {
   	 	return(nrow(mat) > 0 && ncol(mat) > 0)
	}

	if (is_non_empty_matrix(sp_haplotype_matrix) && is_non_empty_matrix(wp_haplotype_matrix) && nrow(sp_haplotype_matrix)>100) {
  		stat_uhl_spwp <- get_stat_pop_from_haplo(sp_haplotype_matrix, wp_haplotype_matrix)
		stat_uhl_spwp_vec <- unlist(stat_uhl_spwp)
		stat_uhl_spwp_vec <- c(stat_uhl_spwp_vec[1], stat_uhl_spwp_vec[2], tail(stat_uhl_spwp_vec, 8))
	} else {
  		stat_uhl_spwp_vec <- rep("NA",10)
	}
	if (is_non_empty_matrix(sp_haplotype_matrix) && is_non_empty_matrix(fu_haplotype_matrix) && nrow(sp_haplotype_matrix)>100) {
  		stat_uhl_spfu <- get_stat_pop_from_haplo(sp_haplotype_matrix, fu_haplotype_matrix)
	  	stat_uhl_spfu_vec <- unlist(stat_uhl_spfu)
	    stat_uhl_spfu_vec <- c(stat_uhl_spfu_vec[1], stat_uhl_spfu_vec[2], tail(stat_uhl_spfu_vec, 8))
	} else {
  	 	stat_uhl_spfu_vec <- rep("NA",10)
	}

	if (is_non_empty_matrix(wp_haplotype_matrix) && is_non_empty_matrix(fu_haplotype_matrix) && nrow(sp_haplotype_matrix)>100) {
  		stat_uhl_wpfu <- get_stat_pop_from_haplo(wp_haplotype_matrix, fu_haplotype_matrix)
	  	stat_uhl_wpfu_vec <- unlist(stat_uhl_wpfu)
	    stat_uhl_wpfu_vec <- c(stat_uhl_wpfu_vec[1], stat_uhl_wpfu_vec[2], tail(stat_uhl_wpfu_vec, 8))
	} else {
  		stat_uhl_wpfu_vec <- rep("NA",10)
	}

	print(stat_uhl_spwp_vec)
	print(stat_uhl_spfu_vec)
	print(stat_uhl_wpfu_vec)

    # Estimate summary statistics from simulated individual data
    system(paste0("bash ind.sh ",seed))

    # Heterozygosity
    file_path_SP=paste0("dr_",k,"/SP_",seed,".txt")
    if (file.exists(file_path_SP)) {
           file_info <-file.info(file_path_SP)
     if (file_info$size > 0) {
           hsp <- read.table(file_path_SP)
           hetsp <- mean(hsp[, 1])
     } else {
           cat("File is empty:", file_path_SP, "\n")
           hetsp <- NA
     }
    } else {
           cat("File not found:", file_path_SP, "\n")
           hetsp <- NA
    }

    file_path_WP=paste0("dr_",k,"/WP_",seed,".txt")
    if (file.exists(file_path_WP)) {
           file_info<-file.info(file_path_WP)
     if (file_info$size > 0) {
            hwp <- read.table(file_path_WP)
            hetwp <- mean(hwp[, 1])
     } else {
            cat("File is empty:", file_path_WP, "\n")
            hetwp <- NA
     }
    } else {
            cat("File not found:", file_path_WP, "\n")
            hetwp <- NA
    }

    file_path_FU=paste0("dr_",k,"/FU_",seed,".txt")
    if (file.exists(file_path_FU)) {
            file_info<-file.info(file_path_FU)
    if (file_info$size > 0) {
            hfu <- read.table(file_path_FU)
            hetfu <- mean(hfu[, 1])
     } else {
            cat("File is empty:", file_path_FU, "\n")
            hetfu <- NA
      }
     } else {
            cat("File not found:", file_path_FU, "\n")
            hetfu <- NA
     }

    # TajimaD
	file_path_SP=paste0("dr_",k,"/SP_",seed,".TajimaD.txt")
          if (file.exists(file_path_SP)) {
                file_info <-file.info(file_path_SP)
            if (file_info$size > 0) {
              TDsp <- read.table(file_path_SP)
              TajDsp <- mean(TDsp[, 1], na.rm=TRUE)
            } else {
              cat("File is empty:", file_path_SP, "\n")
              TajDsp <- NA
            }
          } else {
            cat("File not found:", file_path_SP, "\n")
            TajDsp <- NA
          }

    file_path_WP=paste0("dr_",k,"/WP_",seed,".TajimaD.txt")
          if (file.exists(file_path_WP)) {
                file_info<-file.info(file_path_WP)
            if (file_info$size > 0) {
                TDwp <- read.table(file_path_WP)
                TajDwp <- mean(TDwp[, 1], na.rm=TRUE)
            } else {
              cat("File is empty:", file_path_WP, "\n")
              TajDwp <- NA
            }
          } else {
            cat("File not found:", file_path_WP, "\n")
            TajDwp <- NA
          }

    file_path_FU=paste0("dr_",k,"/FU_",seed,".TajimaD.txt")
          if (file.exists(file_path_FU)) {
                  file_info<-file.info(file_path_FU)
            if (file_info$size > 0) {
              TDfu <- read.table(file_path_FU)
              TajDfu <- mean(TDfu[, 1],na.rm = TRUE)
            } else {
              cat("File is empty:", file_path_FU, "\n")
              TajDfu <- NA
            }
          } else {
            cat("File not found:", file_path_FU, "\n")
            TajDfu <- NA
          }
          
    # SFS	  
    file_path_SP=paste0("dr_",k,"/SP_",seed,".SFS.txt")
          if (file.exists(file_path_SP)) {
                file_info <-file.info(file_path_SP)
            if (file_info$size > 0) {
              sp <- read.table(file_path_SP)
              total_sum_sp <- sum(sp$V1)
              sfssp <- sp$V1/total_sum_sp
            } else {
              cat("File is empty:", file_path_SP, "\n")
              sfssp <- rep("NA",51)
            }
          } else {
            cat("File not found:", file_path_SP, "\n")
            sfssp <- rep("NA",51)
          }

    file_path_WP=paste0("dr_",k,"/WP_",seed,".SFS.txt")
          if (file.exists(file_path_WP)) {
                file_info<-file.info(file_path_WP)
            if (file_info$size > 0) {
              wp <- read.table(file_path_WP)
              total_sum_wp <- sum(wp$V1)
              sfswp <- wp$V1/total_sum_wp
            } else {
              cat("File is empty:", file_path_WP, "\n")
              sfswp <- rep("NA",37)
            }
          } else {
            cat("File not found:", file_path_WP, "\n")
            sfswp <- rep("NA",37)
          }

    file_path_FU=paste0("dr_",k,"/FU_",seed,".SFS.txt")
          if (file.exists(file_path_FU)) {
                  file_info<-file.info(file_path_FU)
            if (file_info$size > 0) {
              fu <- read.table(file_path_FU)
              total_sum_fu <- sum(fu$V1)
              sfsfu <- fu$V1/total_sum_fu
            } else {
              cat("File is empty:", file_path_FU, "\n")
              sfsfu <- rep("NA",21)
            }
          } else {
            cat("File not found:", file_path_FU, "\n")
            sfsfu <- rep("NA",21)
          }
	  
	system(paste0("rm SP_",seed,".txt WP_",seed,".txt FU_",seed,".txt SP_",seed,".TajimaD.txt WP_",seed,".TajimaD.txt FU_",seed,".TajimaD.txt SP_",
	seed,".SFS.txt WP_",seed,".SFS.txt FU_",seed,".SFS.txt"))
	tmp.out=c(fst.wc,fst.id,heteros,hetsp,hetwp,hetfu,f2,f3,pfix,pfix.alt,pfix.ref,TajDsp,TajDwp,TajDfu,stat_uhl_spwp_vec,stat_uhl_spfu_vec,stat_uhl_wpfu_vec,sfssp,sfswp,sfsfu)
	  
	SP_classes <- paste("SPclasse", 0:50, sep = "")
	WP_classes <- paste("WPclasse", 0:36, sep = "")
	FU_classes <- paste("FUclasse", 0:20, sep = "")
	names(tmp.out)=c("FSTwc","FSTidSPWP","FSTidSPFU","FSTidWPFU","HetSP","HetWP","HetFU","het.indSP","het.indWP","het.indFU","f2SPWP","f2SPFU",
	"f2WPFU","f3SP","f3WP","f3FU","PfixSP","PfixWP","PfixFU","Pfix.altSP","Pfix.altWP","Pfix.altFU","Pfix.refSP","Pfix.refWP","Pfix.refFU",
	"DTajimaSP","DTajimaWP","DTajimaFU",'Rf_stat_pop_spwp','Rs_stat_pop_spwp','Wx2s1_stat_pop_spwp','Wx1s2_stat_pop_spwp','Wx1F_stat_pop_spwp','Wx2F_stat_pop_spwp','Wx1_new_stat_pop_spwp',
	'Wx2_new_stat_pop_spwp','Wx1F_new_stat_pop_spwp','Wx2F_new_stat_pop_spwp','Rf_stat_pop_spfu','Rs_stat_pop_spfu','Wx2s1_stat_pop_spfu','Wx1s2_stat_pop_spfu','Wx1F_stat_pop_spfu',
	'Wx2F_stat_pop_spfu','Wx1_new_stat_pop_spfu','Wx2_new_stat_pop_spfu','Wx1F_new_stat_pop_spfu','Wx2F_new_stat_pop_spfu','Rf_stat_pop_wpfu','Rs_stat_pop_wpfu','Wx2s1_stat_pop_wpfu',
	'Wx1s2_stat_pop_wpfu','Wx1F_stat_pop_wpfu','Wx2F_stat_pop_wpfu','Wx1_new_stat_pop_wpfu','Wx2_new_stat_pop_wpfu','Wx1F_new_stat_pop_wpfu','Wx2F_new_stat_pop_wpfu',SP_classes,WP_classes,FU_classes)
	tmp.out
	}

  	# Set prior distribution
	set.seed(seed)

  	# Effective population size
	Ne_ancpp=rlunif(nsample,100000,1000000, base=10)
	Ne_sp=rlunif(nsample,100,100000, base=10)
	Ne_fu=rlunif(nsample,100,200000, base=10)
	Ne_wp=rlunif(nsample,100,10000, base=10)
	Ne_fu_anc=rlunif(nsample,100,200000, base=10)
	Ne_wp_anc=rlunif(nsample,100,100000, base=10)
	Ne_sp_found=rep(0,nsample)
	Ne_wp_found=rep(0,nsample)
	Ne_fu_ancfound=rep(0,nsample)
	Ne_wp_ancfound=rep(0,nsample)
	Ne_bot_sp=rep(0,nsample)
	Ne_bot_wp=rep(0,nsample)
	Ne_bot_fu=rep(0,nsample)
	for(i in 1:nsample){
	  Ne_sp_found[i]=rlunif(1,1,Ne_wp_anc[i], base=10)
	  Ne_wp_found[i]=rlunif(1,1,Ne_wp_anc[i], base=10)
	  Ne_fu_ancfound[i]=rlunif(1,100,Ne_ancpp[i], base=10)
	  Ne_wp_ancfound[i]=rlunif(1,100,Ne_ancpp[i], base=10)
	  Ne_bot_sp[i]=rlunif(1,1,Ne_sp[i], base=10)
	  Ne_bot_wp[i]=rlunif(1,1,Ne_wp[i], base=10)
	  Ne_bot_fu[i]=rlunif(1,1,Ne_fu_anc[i], base=10)
	}

  	# Timing of divergence
  	tsplit_PP=rlunif(nsample,500,20000, base=10)
  	if(tsplit_PP[1]<2000){
	   	tsplit_WP=rlunif(nsample,10,tsplit_PP[1]-200, base=10)
  	}else{
		tsplit_WP=rlunif(nsample,10,tsplit_PP[1]-500, base=10)
  	}
  	t_bot=rlunif(nsample,10,tsplit_WP[1], base=10)

  	# Migration rate
	Nmancpp=(1/rlunif(nsample,0.001,0.5, base=10)-1)/4
	Nmspwp_A=(1/rlunif(nsample,0.001,0.5, base=10)-1)/4
	Nmspfu_A=(1/rlunif(nsample,0.001,0.5, base=10)-1)/4
	Nmfuwp_A=(1/rlunif(nsample,0.001,0.5, base=10)-1)/4
	Nmspwp_R=(1/rlunif(nsample,0.001,0.5, base=10)-1)/4
	Nmspfu_R=(1/rlunif(nsample,0.001,0.5, base=10)-1)/4
	Nmfuwp_R=(1/rlunif(nsample,0.001,0.5, base=10)-1)/4
	m_ancpp=Nmancpp/Ne_fu_anc
	m_spwp_rec=Nmspwp_R/Ne_wp
  	m_spwp_anc=Nmspwp_A/Ne_wp
  	m_spfu_rec=Nmspfu_R/Ne_fu
  	m_spfu_anc=Nmspfu_A/Ne_fu
	m_wpfu_rec=Nmfuwp_R/Ne_fu
  	m_wpfu_anc=Nmfuwp_A/Ne_fu

  	# Creation of reference table for this simulation
	file.name=paste0("reftable.recent_div.",seed)
	SP_classes <- paste("SPclasse", 0:50, sep = "")
	WP_classes <- paste("WPclasse", 0:36, sep = "")
	FU_classes <- paste("FUclasse", 0:20, sep = "")
	cat(file=file.name,c(
	  "sim.Ne_ancpp","sim.Ne_sp","sim.Ne_wp","sim.Ne_fu","sim.Ne_sp_found","sim.Ne_wp_found","sim.Ne_fu_ancfound","sim.Ne_wp_ancfound","sim.Ne_fu_anc","sim.Ne_wp_anc","sim.Ne_bot_sp",
	  "sim.Ne_bot_wp","sim.Ne_bot_fu","sim.tsplit_PP","sim.tsplit_WP","sim.t_bot","sim.m_ancpp","sim.m_spwp_anc","sim.m_spwp_rec","sim.m_spfu_anc","sim.m_spfu_rec","sim.m_wpfu_anc","sim.m_wpfu_rec",
	  "FSTwc","FSTidSPWP","FSTidSPFU","FSTidWPFU","HetSP","HetWP","HetFU","het.indSP","het.indWP","het.indFU","f2SPWP","f2SPFU","f2WPFU","f3SP","f3WP","f3FU","PfixSP","PfixWP","PfixFU",
	  "Pfix.altSP","Pfix.altWP","Pfix.altFU","Pfix.refSP","Pfix.refWP","Pfix.refFU","DTajimaSP","DTajimaWP","DTajimaFU",'Rf_stat_pop_spwp','Rs_stat_pop_spwp','Wx2s1_stat_pop_spwp','Wx1s2_stat_pop_spwp','Wx1F_stat_pop_spwp','Wx2F_stat_pop_spwp','Wx1_new_stat_pop_spwp',
	  'Wx2_new_stat_pop_spwp','Wx1F_new_stat_pop_spwp','Wx2F_new_stat_pop_spwp','Rf_stat_pop_spfu','Rs_stat_pop_spfu','Wx2s1_stat_pop_spfu','Wx1s2_stat_pop_spfu','Wx1F_stat_pop_spfu',
	  'Wx2F_stat_pop_spfu','Wx1_new_stat_pop_spfu','Wx2_new_stat_pop_spfu','Wx1F_new_stat_pop_spfu','Wx2F_new_stat_pop_spfu','Rf_stat_pop_wpfu','Rs_stat_pop_wpfu','Wx2s1_stat_pop_wpfu',
	  'Wx1s2_stat_pop_wpfu','Wx1F_stat_pop_wpfu','Wx2F_stat_pop_wpfu','Wx1_new_stat_pop_wpfu','Wx2_new_stat_pop_wpfu','Wx1F_new_stat_pop_wpfu','Wx2F_new_stat_pop_wpfu',SP_classes,WP_classes,FU_classes),"\n")

  	# Run the simulation
	for(i in 1:nsample){
	  simu=run.simu_recent_div(Ne_ancpp=Ne_ancpp[i],Ne_fu=Ne_fu[i],Ne_wp=Ne_wp[i],Ne_sp=Ne_sp[i],Ne_sp_found=Ne_sp_found[i],Ne_wp_found=Ne_wp_found[i],Ne_fu_ancfound=Ne_fu_ancfound[i],
      			Ne_wp_ancfound=Ne_wp_ancfound[i],Ne_fu_anc=Ne_fu_anc[i],Ne_wp_anc=Ne_wp_anc[i],tsplit_PP=tsplit_PP[i],tsplit_WP=tsplit_WP[i],t_bot=t_bot[i],m_ancpp=m_ancpp[i],
			      m_spwp_anc=m_spwp_anc[i],m_spwp_rec=m_spwp_rec[i],m_spfu_anc=m_spfu_anc[i],m_wpfu_anc=m_wpfu_anc[i],Ne_bot_sp=Ne_bot_sp[i],Ne_bot_wp=Ne_bot_wp[i],Ne_bot_fu=Ne_bot_fu[i],
			      m_spfu_rec=m_spfu_rec[i],m_wpfu_rec=m_wpfu_rec[i],n_sp=50L,n_wp=50L,n_fu=40L,nchr=nchr,chr_length=chr_length,min.rc=min.rc,maf.thr=maf.thr,
			      seed=seed,seed_dr=k, n_cpu=nthreads)
	  simu.par=c(Ne_ancpp[i],Ne_sp[i],Ne_wp[i],Ne_fu[i],Ne_sp_found[i],Ne_wp_found[i],Ne_fu_ancfound[i],Ne_wp_ancfound[i],Ne_fu_anc[i],Ne_wp_anc[i],Ne_bot_sp[i],
		     Ne_bot_wp[i],Ne_bot_fu[i],tsplit_PP[i],tsplit_WP[i],t_bot[i],m_ancpp[i],m_spwp_anc[i],m_spwp_rec[i],m_spfu_anc[i],m_spfu_rec[i],m_wpfu_anc[i],m_wpfu_rec[i])
	  cat(file=file.name,c(simu.par,simu),"\n",append=TRUE)
	}
  rm(H,G)
  gc()
}
