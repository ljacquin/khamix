 # extract haplotype data for analyzed chromosome
R -q --vanilla<get_data_chromo_num_k.R

chromo_num_k=$(cat "chromo_num_k")
nb_snp_chromo_num_k=$(cat "nb_snp_chromo_num_k")
nb_snp_hap=$(cat "nb_snp_hap")

nb_scanned_positions=$(( nb_snp_chromo_num_k-	$((nb_snp_hap-1)) )) 

# genome scan for chromosome $chromo_num_k
for scanned_position in $( seq 1 1 $nb_scanned_positions )  
do        	   
	echo "$scanned_position">scanned_position 
	
	#compute ibs status of haplotypes
	./ibs_haplotypes_window
	./compute_z_h_matrix
	  
	# compute restricted likelihood ratio from emmreml estimated parameters
	R -q --vanilla<emmreml_rlrt.R		  
	cat rlrt_value>>vect_rlrt_value_chromo_num_$chromo_num_k
	cat nb_col_h_matrix>>vect_nb_hap_window_chromo_num_$chromo_num_k
done 
 # end of genome scan for chromosome $chromo_num_k
