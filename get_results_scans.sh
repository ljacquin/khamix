#!/bin/sh
#======================================#
# get results for complete genome scan #
#======================================#
# source R modules
# source module_load_r.sh

# start program for reformatting results
if [ ! -d results_genome_scan ]; then
 	mkdir results_genome_scan
fi
cd data_parameters/

nb_chromosomes=$(cat "nb_chromosomes.txt")
nb_snp_hap=$(cat "nb_snp_hap.txt")
kernel_index=$(cat "kernel_index.txt")

modify_signif_level=true	# true or false (note lower case)
new_signif_level=0.01

if [ "$modify_signif_level"=true ] ; then
    echo "$new_signif_level" > signif_level.txt
	if [ -d results_genome_scan ]; then
		rm -rf results_genome_scan/*
	fi
fi


cp trait_name.txt		../results_genome_scan
cp nb_snp_hap.txt		../results_genome_scan
cp nb_chromosomes.txt		../results_genome_scan
cp physical_map.txt		../results_genome_scan
cp phased_genotypes.txt		../results_genome_scan
cp genotypes.txt		../results_genome_scan
cp kernel_index.txt		../results_genome_scan
cp signif_level.txt		../results_genome_scan

cd ../programs/
cp mult_test_correction.R	../results_genome_scan
cp get_results_scans.R		../results_genome_scan
cp plot_nb_hap_scans.R		../results_genome_scan
cp plot_manhattan_scan.R	../results_genome_scan
cd ../

#--------------------------------#
# moving through all directories #
#--------------------------------#

for chromo_num_k in $(seq 1 1 $nb_chromosomes)
do
	cd genome_scan_chromo_num_$chromo_num_k
	cp vect_rlrt_value_chromo_num_$chromo_num_k.txt		../results_genome_scan
	cp vect_nb_hap_window_chromo_num_$chromo_num_k.txt	../results_genome_scan
	cd ../
done

cd results_genome_scan

R -q --vanilla < mult_test_correction.R
R -q --vanilla < get_results_scans.R
R -q --vanilla < plot_nb_hap_scans.R
R -q --vanilla < plot_manhattan_scan.R

if [ "$kernel_index" -gt 1 ] ; then

	if [ "$nb_snp_hap" -gt 1 ] ; then

		for chromo_num_k in $(seq $nb_chromosomes -1 1)
		do
			if [ -d results_chromo_num_$chromo_num_k ]; then
				rm -rf results_chromo_num_$chromo_num_k/*
			else
				mkdir results_chromo_num_$chromo_num_k
			fi
			mv flanking_markers_of_tested_positions_in_kb_with_rlrt_value_on_chromosome_$chromo_num_k*		results_chromo_num_$chromo_num_k
			mv flanking_markers_of_tested_positions_in_kb_with_significant_rlrt_value_on_chromosome_$chromo_num_k*	results_chromo_num_$chromo_num_k
			if ls significant_haplotypes_chromo_num_$chromo_num_k* 1> /dev/null 2>&1; then
				mv significant_haplotypes_chromo_num_$chromo_num_k*				results_chromo_num_$chromo_num_k
			fi
			mv vect_rlrt_value_chromo_num_$chromo_num_k*						results_chromo_num_$chromo_num_k
			mv vect_nb_hap_window_chromo_num_$chromo_num_k*						results_chromo_num_$chromo_num_k
			mv number_of_haplotypes_per_window_for_chromosome_$chromo_num_k*			results_chromo_num_$chromo_num_k
			mv kernelized_haplotype_based_genome_scan_of_chromosome_$chromo_num_k*			results_chromo_num_$chromo_num_k
		done

		if [ -d results_all_chromosomes ]; then
			rm -rf results_all_chromosomes/*
		else
			mkdir results_all_chromosomes
		fi
		mv number_of_haplotypes_per_window_for_complete_genome_scan*					results_all_chromosomes
		mv kernelized_haplotype_based_genome_scan_for_*							results_all_chromosomes
		mv qq_plot_kernelized_haplotype_*								results_all_chromosomes
		mv flanking_markers_of_tested_positions_with_statistics*					results_all_chromosomes

	else

		for chromo_num_k in $(seq $nb_chromosomes -1 1)
		do
			if [ -d results_chromo_num_$chromo_num_k ]; then
				rm -rf results_chromo_num_$chromo_num_k/*
			else
				mkdir results_chromo_num_$chromo_num_k
			fi
			mv markers_in_kb_with_rlrt_value_on_chromosome_$chromo_num_k*				results_chromo_num_$chromo_num_k
			mv markers_in_kb_with_significant_rlrt_value_on_chromosome_$chromo_num_k*		results_chromo_num_$chromo_num_k
			if ls significant_snps_chromo_num_$chromo_num_k* 1> /dev/null 2>&1; then
				mv significant_snps_chromo_num_$chromo_num_k*					results_chromo_num_$chromo_num_k
			fi
			mv vect_rlrt_value_chromo_num_$chromo_num_k*						results_chromo_num_$chromo_num_k
			mv kernelized_gwas_of_chromosome_$chromo_num_k*						results_chromo_num_$chromo_num_k
		done

		if [ -d results_all_chromosomes ]; then
			rm -rf results_all_chromosomes/*
		else
			mkdir results_all_chromosomes
		fi
		mv kernelized_gwas_for_*									results_all_chromosomes
		mv qq_plot_kernelized_gwas_*									results_all_chromosomes
		mv markers_of_tested_positions_with_statistics.txt*						results_all_chromosomes
		rm vect_nb_hap_window_chromo_num_*

	fi

 else

	if [ "$nb_snp_hap" -gt 1 ] ; then

		for chromo_num_k in $(seq $nb_chromosomes -1 1)
		do
			if [ -d results_chromo_num_$chromo_num_k ]; then
				rm -rf results_chromo_num_$chromo_num_k/*
			else
				mkdir results_chromo_num_$chromo_num_k
			fi
			mv flanking_markers_of_tested_positions_in_kb_with_rlrt_value_on_chromosome_$chromo_num_k*		results_chromo_num_$chromo_num_k
			mv flanking_markers_of_tested_positions_in_kb_with_significant_rlrt_value_on_chromosome_$chromo_num_k*	results_chromo_num_$chromo_num_k
			if ls significant_haplotypes_chromo_num_$chromo_num_k* 1> /dev/null 2>&1; then
				mv significant_haplotypes_chromo_num_$chromo_num_k*				results_chromo_num_$chromo_num_k
			fi
			mv vect_rlrt_value_chromo_num_$chromo_num_k*						results_chromo_num_$chromo_num_k
			mv vect_nb_hap_window_chromo_num_$chromo_num_k*						results_chromo_num_$chromo_num_k
			mv number_of_haplotypes_per_window_for_chromosome_$chromo_num_k*			results_chromo_num_$chromo_num_k
			mv haplotype_based_genome_scan_of_chromosome_$chromo_num_k*				results_chromo_num_$chromo_num_k
		done

		if [ -d results_all_chromosomes ]; then
			rm -rf results_all_chromosomes/*
		else
			mkdir results_all_chromosomes
		fi
		mv number_of_haplotypes_per_window_for_complete_genome_scan*					results_all_chromosomes
		mv haplotype_based_genome_scan_for_*								results_all_chromosomes
		mv qq_plot_haplotype_*										results_all_chromosomes
		mv flanking_markers_of_tested_positions_with_statistics*					results_all_chromosomes

	else

		for chromo_num_k in $(seq $nb_chromosomes -1 1)
		do
			if [ -d results_chromo_num_$chromo_num_k ]; then
				rm -rf results_chromo_num_$chromo_num_k/*
			else
				mkdir results_chromo_num_$chromo_num_k
			fi
			mv markers_in_kb_with_rlrt_value_on_chromosome_$chromo_num_k*				results_chromo_num_$chromo_num_k
			mv markers_in_kb_with_significant_rlrt_value_on_chromosome_$chromo_num_k*		results_chromo_num_$chromo_num_k
			if ls significant_snps_chromo_num_$chromo_num_k* 1> /dev/null 2>&1; then
				mv significant_snps_chromo_num_$chromo_num_k*					results_chromo_num_$chromo_num_k
			fi
			mv vect_rlrt_value_chromo_num_$chromo_num_k*						results_chromo_num_$chromo_num_k
			mv gwas_of_chromosome_$chromo_num_k*							results_chromo_num_$chromo_num_k
		done

		if [ -d results_all_chromosomes ]; then
			rm -rf results_all_chromosomes/*
		else
			mkdir results_all_chromosomes
		fi
		mv gwas_for_*											results_all_chromosomes
		mv qq_plot_gwas_*										results_all_chromosomes
		mv markers_of_tested_positions_with_statistics.txt*						results_all_chromosomes
		rm vect_nb_hap_window_chromo_num_*

	fi

fi
clear

rm *.R
rm *.txt

# End of program for reformatting results

