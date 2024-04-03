#!/bin/sh
#======================================#
# get results for complete genome scan #
#======================================#
# source R modules
# source module_load_r.sh

# start program for reformatting results
mkdir results_genome_scan
cd data_parameters/

nb_chromosomes=$(cat "nb_chromosomes")    
nb_snp_hap=$(cat "nb_snp_hap")    
kernel_index=$(cat "kernel_index")

cp trait_name                    ../results_genome_scan
cp nb_snp_hap                    ../results_genome_scan
cp nb_chromosomes                ../results_genome_scan
cp physical_map.txt              ../results_genome_scan
cp phased_genotypes.txt          ../results_genome_scan
cp kernel_index                  ../results_genome_scan
cp signif_level                  ../results_genome_scan
 
cd ../programs/
cp plot_results_scans.R          ../results_genome_scan
cp plot_nb_hap_scans.R           ../results_genome_scan
cp get_significant_haplotypes.R  ../results_genome_scan
cd ../

#--------------------------------#
# moving through all directories #
#--------------------------------#

for chromo_num_k in $(seq 1 1 $nb_chromosomes)
do
	cd genome_scan_chromo_num_$chromo_num_k
	cp vect_rlrt_value_chromo_num_$chromo_num_k	../results_genome_scan
	cp vect_nb_hap_window_chromo_num_$chromo_num_k	../results_genome_scan
	cd ../
done

cd results_genome_scan
R -q --vanilla < plot_results_scans.R
R -q --vanilla < plot_nb_hap_scans.R
R -q --vanilla < get_significant_haplotypes.R


if [ "$kernel_index" -gt 1 ] ; then

	if [ "$nb_snp_hap" -gt 1 ] ; then

		for chromo_num_k in $(seq $nb_chromosomes -1 1)
		do
			mkdir results_chromo_num_$chromo_num_k
			mv flanking_markers_of_tested_positions_in_kb_with_rlrt_value_on_chromosome_$chromo_num_k*		results_chromo_num_$chromo_num_k
			mv flanking_markers_of_tested_positions_in_kb_with_significant_rlrt_value_on_chromosome_$chromo_num_k*	results_chromo_num_$chromo_num_k	
			mv significant_haplotypes_chromo_num_$chromo_num_k*							results_chromo_num_$chromo_num_k
			mv vect_rlrt_value_chromo_num_$chromo_num_k*								results_chromo_num_$chromo_num_k 
			mv vect_nb_hap_window_chromo_num_$chromo_num_k*								results_chromo_num_$chromo_num_k 
			mv number_of_haplotypes_per_window_for_chromosome_$chromo_num_k*					results_chromo_num_$chromo_num_k
			mv kernelized_haplotype_based_genome_scan_of_chromosome_$chromo_num_k*					results_chromo_num_$chromo_num_k
		done

		mkdir results_all_chromosomes
		mv number_of_haplotypes_per_window_for_complete_genome_scan*							results_all_chromosomes
		mv kernelized_haplotype_based_genome_scan_for_*									results_all_chromosomes

	else

		for chromo_num_k in $(seq $nb_chromosomes -1 1)
		do
			mkdir results_chromo_num_$chromo_num_k
			mv markers_in_kb_with_rlrt_value_on_chromosome_$chromo_num_k*						results_chromo_num_$chromo_num_k
			mv markers_in_kb_with_significant_rlrt_value_on_chromosome_$chromo_num_k*				results_chromo_num_$chromo_num_k
			mv significant_snps_chromo_num_$chromo_num_k*								results_chromo_num_$chromo_num_k
			mv vect_rlrt_value_chromo_num_$chromo_num_k*								results_chromo_num_$chromo_num_k 
			mv kernelized_gwas_of_chromosome_$chromo_num_k*								results_chromo_num_$chromo_num_k
		done

		mkdir results_all_chromosomes                                              
		mv kernelized_gwas_for_*											results_all_chromosomes

	fi

 else
 
	if [ "$nb_snp_hap" -gt 1 ] ; then

		for chromo_num_k in $(seq $nb_chromosomes -1 1)
		do
			mkdir results_chromo_num_$chromo_num_k
			mv flanking_markers_of_tested_positions_in_kb_with_rlrt_value_on_chromosome_$chromo_num_k*		results_chromo_num_$chromo_num_k
			mv flanking_markers_of_tested_positions_in_kb_with_significant_rlrt_value_on_chromosome_$chromo_num_k*	results_chromo_num_$chromo_num_k	
			mv significant_haplotypes_chromo_num_$chromo_num_k*							results_chromo_num_$chromo_num_k
			mv vect_rlrt_value_chromo_num_$chromo_num_k*								results_chromo_num_$chromo_num_k 
			mv vect_nb_hap_window_chromo_num_$chromo_num_k*								results_chromo_num_$chromo_num_k 
			mv number_of_haplotypes_per_window_for_chromosome_$chromo_num_k*					results_chromo_num_$chromo_num_k
			mv haplotype_based_genome_scan_of_chromosome_$chromo_num_k*						results_chromo_num_$chromo_num_k
		done

		mkdir results_all_chromosomes
		mv number_of_haplotypes_per_window_for_complete_genome_scan*							results_all_chromosomes
		mv haplotype_based_genome_scan_for_*										results_all_chromosomes

	else

		for chromo_num_k in $(seq $nb_chromosomes -1 1)
		do
			mkdir results_chromo_num_$chromo_num_k
			mv markers_in_kb_with_rlrt_value_on_chromosome_$chromo_num_k*						results_chromo_num_$chromo_num_k
			mv markers_in_kb_with_significant_rlrt_value_on_chromosome_$chromo_num_k*				results_chromo_num_$chromo_num_k
			mv significant_snps_chromo_num_$chromo_num_k*								results_chromo_num_$chromo_num_k
			mv vect_rlrt_value_chromo_num_$chromo_num_k*								results_chromo_num_$chromo_num_k 
			mv gwas_of_chromosome_$chromo_num_k*									results_chromo_num_$chromo_num_k
		done

		mkdir results_all_chromosomes                                              
		mv gwas_for_*													results_all_chromosomes

	fi

    

fi    
clear

# End of program for reformatting results

