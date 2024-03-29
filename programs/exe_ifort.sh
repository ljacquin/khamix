 # create executables
 # module load compiler/intel-2013.3.174

 echo 

  /opt/intel/oneapi/compiler/latest/bin/ifort -traceback -O3  ibs_haplotypes_window.f90
  mv a.out ibs_haplotypes_window
 
  /opt/intel/oneapi/compiler/latest/bin/ifort -traceback -O3  compute_z_h_matrix.f90
  mv a.out compute_z_h_matrix

 

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
