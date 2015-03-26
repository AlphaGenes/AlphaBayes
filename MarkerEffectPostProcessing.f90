!###########################################################################################################################################################

subroutine MarkerEffectPostProcessing
use Global
implicit none

real(4) :: sdot,myone,myzero,TmpVal
real(4) :: Correlation
integer :: i,j,One


myone=1.d0
myzero=0.d0
One=1

! Rescale back to phenotype scale
G(:,1)=G(:,1)*sqrt(VarY)

!Output section
open (unit=1002,file="SnpSolutions.txt",status="unknown")

allocate(SnpOut(nSnpExternal))

SnpOut=0.d0
j=0
do i=1,nSnpExternal
	if (FixedSnp(i)==1) then
		j=j+1
		SnpOut(i)=G(j,1)
	endif
	write (1002,*) i,SnpOut(i)
enddo

close(1002)
flush(1002)

open (unit=1001,file="TbvEbv.txt",status="unknown")

call sgemm('n','n',nAnisTe,One,nSnp,myone,GenosTe,nAnisTe,G,nSnp,myzero,Ebv,nAnisTe)
do i=1,nAnisTe
	write(1001,*) i,Tbv(i,1),Ebv(i,1)
enddo
close(1001)

call PearsnR4 (Tbv(:,1),Ebv(:,1),nAnisTe,Correlation)
print*, Correlation
open (unit=1001,file="TbvEbvCorrelation.txt",status="unknown")
  write(1001,*) Correlation
flush(1002)
close(1001)

!call system("sleep 10")

end subroutine MarkerEffectPostProcessing

!###########################################################################################################################################################
