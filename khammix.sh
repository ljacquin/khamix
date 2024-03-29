#!/bin/sh
#============================================================#
#  kernelized qtl-haplotype mapping by mixed model (khammix) #										
#============================================================#
# source modules for R
# source module_load_r.sh

#--------------------------------------------------------------------------#
#1. setting haplotype size, in number of markers, and trait to be analyzed #
#--------------------------------------------------------------------------#
trait_name="LLGTH"
nb_snp_hap=6			
nb_chromosomes=12		# 12
kernel_index=1			# 1 for VanRaden linear kernel and 2 for Gaussian kernel (i.e. RBF)		        
local_or_cluster=1		# 1 for local computation and 2 for parallelize computation on a cluster 
	
#echo "$trait_name" > trait_name
echo "$trait_name" > trait_name
echo "$nb_snp_hap" > nb_snp_hap
echo "$nb_chromosomes" > nb_chromosomes
echo "$kernel_index" > kernel_index

mv trait_name		data_parameters/
mv nb_snp_hap		data_parameters/
mv nb_chromosomes	data_parameters/
mv kernel_index		data_parameters/

#-------------------------------------------------------------------------#
#3. genome scan, for each chromosome, associated with the analyzed trait  #
#-------------------------------------------------------------------------#

for chromo_num_k in $(seq 1 1 $nb_chromosomes)
 do
  echo "$chromo_num_k">chromo_num_k

  #-------------------------------------------------------------------------#
  #3.1 creation of the directory for chromosome chromo_num_k to be analyzed #
  #-------------------------------------------------------------------------#

	mkdir genome_scan_chromo_num_$chromo_num_k
	
	# copy data into the directory for chromosome chromo_num_k
	mv chromo_num_k					genome_scan_chromo_num_$chromo_num_k
	cp data_parameters/trait_name			genome_scan_chromo_num_$chromo_num_k
	cp data_parameters/phenotypes.txt		genome_scan_chromo_num_$chromo_num_k
	cp data_parameters/incidence_fixed_effects.txt	genome_scan_chromo_num_$chromo_num_k
	cp data_parameters/nb_snp_hap			genome_scan_chromo_num_$chromo_num_k
	cp data_parameters/nb_chromosomes		genome_scan_chromo_num_$chromo_num_k
	cp data_parameters/genotypes.txt		genome_scan_chromo_num_$chromo_num_k
	cp data_parameters/phased_genotypes.txt		genome_scan_chromo_num_$chromo_num_k
	cp data_parameters/physical_map.txt		genome_scan_chromo_num_$chromo_num_k 
	cp data_parameters/kernel_index			genome_scan_chromo_num_$chromo_num_k

	# copy programs into the directory for chromosome chromo_num_k
	cp programs/get_data_trait_name.R		genome_scan_chromo_num_$chromo_num_k
	cp programs/get_data_fixed_effects.R		genome_scan_chromo_num_$chromo_num_k
	cp programs/compute_k_matrix.R			genome_scan_chromo_num_$chromo_num_k
	cp programs/get_data_chromo_num_k.R		genome_scan_chromo_num_$chromo_num_k
	cp programs/ibs_haplotypes_window		genome_scan_chromo_num_$chromo_num_k
	cp programs/compute_z_h_matrix			genome_scan_chromo_num_$chromo_num_k
	cp programs/emmreml_rlrt.R			genome_scan_chromo_num_$chromo_num_k
  	cp programs/scan_chromo_num_k.sh		genome_scan_chromo_num_$chromo_num_k  	
  	
	# extract data associated to trait, fixed effects and compute Gram matrix (i.e. kernel matrix)
	cd genome_scan_chromo_num_$chromo_num_k
	R -q --vanilla<get_data_trait_name.R
	R -q --vanilla<get_data_fixed_effects.R
	R -q --vanilla<compute_k_matrix.R	

  	# perform genome scan by sliding window for chromosome $chromo_num_k 
	if [ "$kernel_index" -gt 1 ] ; then
		echo 	
		#qsub -q normal.q -v scan_chromo_num_k.sh		# a cc2 queue
		qsub -q workq scan_chromo_num_k.sh			# a cc2 queue
	else
		echo
		./scan_chromo_num_k.sh
	fi
	   
	# moving to the next directory "genome_scan_chromo_num_$chromo_num_k",
	# i.e. moving up the directory tree
	cd ../

done  
