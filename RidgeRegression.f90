!###########################################################################################################################################################

subroutine RidgeRegression
use Global
implicit none

real(4) :: sdot,eps,InvLhs,Rhs,Lhs,SolOld,myone,myzero,TmpVal
real(4) :: Correlation
integer :: i,h,j,snpid,RandomOrdering(nSnp),One

allocate(XpX(nSnp,1))
allocate(Xg(nAnisTr,1))
allocate(Lambda(nSnp))

myone=1.d0
myzero=0.d0
One=1

!Construct XpX
do j=1,nSnp
	XpX(j,1)=sdot(nAnisTr, GenosTr(:,j), 1, GenosTr(:,j), 1) + 0.000000000000001
	Lambda(j)=VarE/(VarA/nSnp)
enddo  

Mu=0.00
E(:,1)=Phen(:,1)-Mu !GG

do h=1,nRound
	eps=0.0	

	!Intercept
	E(:,1)=E(:,1)+Mu
	Rhs=sum(E(:,1))
	InvLhs=1.0/float(nAnisTr)
	Mu=Rhs*InvLhs
	E(:,1)=E(:,1)-Mu

	!Snp effects
	call random_order(RandomOrdering,nSnp,idum)
	do j=1,nSnp
		snpid=RandomOrdering(j)

		E(:,1)=E(:,1)+(GenosTr(:,snpid)*G(snpid,1)) 

		Lhs=XpX(snpid,1)+Lambda(snpid)
		Rhs=sdot(nAnisTr, GenosTr(:,snpid), 1, E(:,1), 1)

		SolOld=G(snpid,1)
		G(snpid,1)=Rhs/Lhs

		E(:,1)=E(:,1)-(GenosTr(:,snpid)*G(snpid,1))
		Eps=Eps+(G(snpid,1)-SolOld)**2			

	enddo

	if (mod(h,200)==0) then
		call sgemm('n','n',nAnisTr,One,nSnp,myone,GenosTr,nAnisTr,G,nSnp,myzero,Xg,nAnisTr)
		E(:,1)=Phen(:,1)-Xg(:,1)-Mu
	endif

	if (Eps.lt.1e-16) then
		ConvergedRounds=h
		exit
	endif

enddo

print *, "Converged in ", ConvergedRounds," rounds"

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
close(1001)

end subroutine RidgeRegression

!###########################################################################################################################################################
