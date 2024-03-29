program compute_z_h_matrix

use kinds

implicit none

!========================!
!== declared variables ==!	
!========================!
integer::i,j,nbhaps,io,indice,rank_h_matrix,indic,k,nb_elem,p 	  !deb_fen= starting index of window 
integer,allocatable::ibs_window(:),tab(:,:),W(:,:)
real(r8),allocatable::h_matrix(:,:)

!====================================================!
!- read the IBS status of haplotypes for the window -!
!====================================================!
open(1, file='ibs_status_haplotypes_window')
nbhaps=0
	    do
		read(1,*,iostat=io)
		  if(io/=0)exit
		  nbhaps=nbhaps+1
	    end do		 
	! print*,nbhaps	
rewind(1)
allocate(ibs_window(nbhaps))

   do i=1,nbhaps   		
	read(1,*) ibs_window(i)
	!print*, ibs_window(i)
   end do          
 close(1)

rank_h_matrix=maxval(ibs_window)

!**************** 
nb_elem=nbhaps/2
allocate(tab(nb_elem,2))
p=1
k=2   
   do i=1,nbhaps          
      if(k==2) then
         tab(p,1)=ibs_window(i)
	 tab(p,2)=ibs_window(i+1)
	  p=p+1
          k=0
      end if
    k=k+1   
   end do 
!****************

!====================================================!
!- write the incidence matrix for haplotype effects -!
!====================================================!
open(24,file='z_h_matrix')
allocate(W(nb_elem,rank_h_matrix))
W=0 
do p=1,nb_elem    
   do i=1,rank_h_matrix           
            if (i==tab(p,1)) then
               W(p,i)=1
	       indic=i
             end if     
            if (i==tab(p,2)) then
               W(p,i)=1
	      indic=i
             end if     
   end do     
 if ((sum(W(p,:)))/=2) W(p,indic)=2  
end do   
 do i=1,nb_elem
     write(24,'(1000000(i2,1x))')(W(i,j),j=1,rank_h_matrix)
 end do       
 close(24) 

!=======================================================!
!- write the covariance h matrix for haplotype effects -!
!=======================================================!
allocate(h_matrix(rank_h_matrix,rank_h_matrix))
open(4,file='h_matrix')
do i=1,rank_h_matrix
     do j=1,rank_h_matrix
	
	  if (i==j) then
	         h_matrix(i,j)=1.0
	    else
	         h_matrix(i,j)=0.0
	  end if
     
       end do
end do
   do i=1,rank_h_matrix
	write(4,'(100(30F8.4,1x))')(h_matrix(i,j),j=1,rank_h_matrix)
	
   end do
 close(4)

deallocate(ibs_window)
deallocate(h_matrix)
deallocate(tab,W)

end program compute_z_h_matrix
