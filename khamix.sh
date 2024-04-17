#!/bin/sh
#=======================================================================#
#  kernelized haplotype-based mixed model association mapping (KHAMIX)  #
#=======================================================================#
# source modules for R
# source module_load_r.sh

#--------------------------------------------------------------------------#
#1. setting haplotype size, in number of markers, and trait to be analyzed #
#--------------------------------------------------------------------------#
trait_name="LLGTH"
nb_snp_hap=6
nb_chromosomes=12		# total number of chromosomes for the analyzed data set
kernel_index=1			# 1 for VanRaden linear kernel and 2 for Gaussian kernel (i.e. RBF)
signif_level=0.01		# significance level for the restricted likelihood ratio test (RLRT)
local_or_cluster=2		# 1 for local computation and 2 for parallelize computation on a cluster

echo "$trait_name" > trait_name.txt
echo "$nb_snp_hap" > nb_snp_hap.txt
echo "$nb_chromosomes" > nb_chromosomes.txt
echo "$kernel_index" > kernel_index.txt
echo "$signif_level" > signif_level.txt

mv trait_name.txt		data_parameters/
mv nb_snp_hap.txt		data_parameters/
mv nb_chromosomes.txt		data_parameters/
mv kernel_index.txt		data_parameters/
mv signif_level.txt		data_parameters/

#------------------------------------------------#
#2. compute estimates under null hypothesis (h0) #
#------------------------------------------------#
mkdir estimates_h0/
cp data_parameters/trait_name.txt				estimates_h0/
cp data_parameters/phenotypes.txt				estimates_h0/
if [ -f data_parameters/incidence_fixed_effects.txt ] ; then
    cp data_parameters/incidence_fixed_effects.txt		estimates_h0/
fi
if [ -f data_parameters/incidence_polygenic_effects.txt ] ; then
    cp data_parameters/incidence_polygenic_effects.txt	estimates_h0/
fi
cp data_parameters/kernel_index.txt				estimates_h0/
cp data_parameters/genotypes.txt				estimates_h0/
cp programs/get_data_trait_name.R				estimates_h0/
cp programs/get_data_incidence.R				estimates_h0/
cp programs/compute_k_matrix.R					estimates_h0/
cp programs/compute_estimates_h0.R				estimates_h0/
cd estimates_h0/
	R -q --vanilla<get_data_trait_name.R
	R -q --vanilla<get_data_incidence.R
	R -q --vanilla<compute_k_matrix.R
	R -q --vanilla<compute_estimates_h0.R
cd ../

#-------------------------------------------------------------------------#
#3. genome scan, for each chromosome, associated with the analyzed trait  #
#-------------------------------------------------------------------------#

for chromo_num_k in $(seq 1 1 $nb_chromosomes)
do
	echo "$chromo_num_k">chromo_num_k.txt

	#-------------------------------------------------------------------------#
	#3.1 creation of the directory for chromosome chromo_num_k to be analyzed #
	#-------------------------------------------------------------------------#
	mkdir genome_scan_chromo_num_$chromo_num_k

	# copy data into the directory for chromosome chromo_num_k
	mv chromo_num_k.txt					genome_scan_chromo_num_$chromo_num_k
	cp data_parameters/nb_snp_hap.txt			genome_scan_chromo_num_$chromo_num_k
	cp data_parameters/nb_chromosomes.txt			genome_scan_chromo_num_$chromo_num_k
	cp data_parameters/genotypes.txt			genome_scan_chromo_num_$chromo_num_k
	cp data_parameters/phased_genotypes.txt			genome_scan_chromo_num_$chromo_num_k
	cp data_parameters/physical_map.txt			genome_scan_chromo_num_$chromo_num_k

	# copy data and estimates for model under null hypothesis (h0) into directory for chromosome chromo_num_k
	cp estimates_h0/phenotypes_trait_name.txt		genome_scan_chromo_num_$chromo_num_k
	cp estimates_h0/x_matrix.txt				genome_scan_chromo_num_$chromo_num_k
	cp estimates_h0/z_u_matrix.txt				genome_scan_chromo_num_$chromo_num_k
	cp estimates_h0/k_matrix.txt				genome_scan_chromo_num_$chromo_num_k
	cp estimates_h0/emmreml_h0				genome_scan_chromo_num_$chromo_num_k

	# copy programs into the directory for chromosome chromo_num_k
	cp programs/get_data_chromo_num_k.R			genome_scan_chromo_num_$chromo_num_k
	cp programs/ibs_haplotypes_window			genome_scan_chromo_num_$chromo_num_k
	cp programs/compute_z_h_matrix				genome_scan_chromo_num_$chromo_num_k
	cp programs/emmreml_rlrt.R				genome_scan_chromo_num_$chromo_num_k
  	cp programs/scan_chromo_num_k.sh			genome_scan_chromo_num_$chromo_num_k

	# extract data associated to trait, fixed effects and compute Gram matrix (i.e. kernel matrix)
	cd genome_scan_chromo_num_$chromo_num_k

  	# perform genome scan by sliding window for chromosome $chromo_num_k
	if [ "$local_or_cluster" -gt 1 ] ; then
		# qsub -q normal.q -v scan_chromo_num_k.sh  # a cc2 queue
		sbatch scan_chromo_num_k.sh
		# qsub -q workq scan_chromo_num_k.sh
	else
		./scan_chromo_num_k.sh
	fi

	# moving to the next directory "genome_scan_chromo_num_$chromo_num_k",
	# i.e. moving up the directory tree
	cd ../
done

