program ibs_haplotypes_window

use kinds

implicit none

!========================!
!== declared variables ==!	
!========================!
integer::i, j, k, m, s, nb_hap, io, rankH, max_haplo, index_col, deb_fen, nb_snp_hap, lg, bool, nb_elem, p, indic, nb_snp_chromo_num_k	
integer,allocatable::haplotypes(:,:),haplotypes_temp(:,:),W(:,:)
integer,allocatable::num_ind(:,:)
integer,allocatable::ibs_window(:),tab(:,:),classe_dip(:)

!=====================!
!== read parameters ==!	
!=====================!
open(1,file='scanned_position')
read(1,*)index_col
close(1)

open(28,file='nb_snp_chromo_num_k')
read(28,*) nb_snp_chromo_num_k
close(28)


open(50,file='nb_snp_hap')
read(50,*)nb_snp_hap
close(50)

!============================!
!== read phased haplotypes ==!	
!============================!
open(2, file='haplotypes_chromo_num_k')
nb_hap=0
	    do
		read(2,*,iostat=io)
		  if(io/=0)exit
		  nb_hap=nb_hap+1
	    end do	
	       
rewind(2)
allocate( haplotypes_temp(nb_hap, nb_snp_chromo_num_k) )			
   do i=1,nb_hap		
	read(2,*) haplotypes_temp(i,:)			
   end do        
 close(2)

!=========================================================================================================!
!-search for segregating haplotypes that are identical by state (IBS) over a window of nb_snp_hap markers !
!=========================================================================================================!
allocate(haplotypes(nb_hap,nb_snp_hap))
haplotypes=0

haplotypes(1:nb_hap,1:nb_snp_hap) = haplotypes_temp( 1:nb_hap, index_col:(index_col+nb_snp_hap-1) )

allocate(ibs_window(nb_hap))	
ibs_window=0
rankH=1   
s=1       

do i=1,nb_hap 
             if (ibs_window(i)==0) then
          		  ibs_window(i)=s			
		           if (s>rankH) rankH=s		    
		    	   s=s+1
              end if 	      	       	      
         do k=1,nb_hap    
	   lg=0	     		     
		 do j=1,nb_snp_hap
		    if (haplotypes(i,j)==haplotypes(k,j)) lg=lg+1
		  end do  
		        	              
                 if ((ibs_window(k)==0).and.(lg==nb_snp_hap)) then		 
	    		ibs_window(k)=ibs_window(i)			 
	         end if            
	 end do  	 
end do       

!============================================!
!- write the IBS status of these haplotypes -!
!============================================!
open(2, file='ibs_status_haplotypes_window')
   do i=1,nb_hap   		
	write(2,'(1(i5,1x))')(ibs_window(i))
   end do         
 close(2)

open(25, file='nb_col_h_matrix')
rankH=maxval(ibs_window)
write(25,'(i3)')(rankH)
 close(25)

!****************
do j=1,nb_hap
  do i=1,nb_hap
  	  if ((ibs_window(j)==ibs_window(i)).and.(j/=i)) then
	                      ibs_window(i)=0
	   end if
   end do
 end do
!****************

open(8, file='haplotype_classes')
   do i=1,rankH      	 
	   do j=1,nb_hap
	        
	        if    (i==ibs_window(j)) then
			write(8,'(1000i2)')(haplotypes(j,k),k=1,nb_snp_hap)
		end if  
			  	      
	  end do	     
    end do      	          
 close(8)

deallocate(ibs_window, haplotypes, haplotypes_temp)

end program ibs_haplotypes_window
