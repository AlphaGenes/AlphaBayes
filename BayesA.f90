!###########################################################################################################################################################

subroutine BayesA
use Global
implicit none

real(8) :: sdot,InvLhs,Rhs,Lhs,SolOld,myone,myzero,vdf,ScaleCoefficient,LambdaSnp,EpE
real(8),allocatable,dimension(:,:) :: GAccum,Gvar
real(8) :: random_gamma,gasdev,ShapeCoefficient,sc,gampe1,gamp2
integer :: i,h,j,snpid,RandomOrdering(nSnp),One,TotalRoundsContributing

allocate(XpX(nSnp,1))
allocate(Xg(nAnisTr,1))
allocate(Lambda(nSnp))
allocate(Gvar(nSnp,1))
allocate(GAccum(nSnp,1))

myone=1.d0
myzero=0.d0
One=1
TotalRoundsContributing=0
GAccum=0
G=0.1

vdf=4.012
ShapeCoefficient=(vdf/2) !Ben
ScaleCoefficient=((vdf-2)*VarA)/Sum2pq
sc=2.0
gampe1=float((nAnisTr)/2)-2.0
gamp2=2.0

!Construct XpX
do j=1,nSnp
	XpX(j,1)=sdot(nAnisTr, GenosTr(:,j), 1, GenosTr(:,j), 1) + 0.000000000000001
	Lambda(j)=VarE/(VarA/nSnp)
enddo

Mu=0.00
E(:,1)=Phen(:,1)-Mu !GG


open (unit=6,form='formatted',CARRIAGECONTROL='FORTRAN') 

do h=1,nRound

	EpE=sdot(nAnisTr, E(:,1), 1, E(:,1), 1) + 0.000000000000001
	VarE=EpE/random_gamma(idum,gampe1,gamp2,.TRUE.)

	write(6, 100) "   BayesA round   ",h
	100 format ('+', a17,i10)    

	do j=1,nSnp
		GVar(j,1)=(ScaleCoefficient+(G(j,1)**2))/random_gamma(idum,ShapeCoefficient,sc,.TRUE.)	
	enddo	

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
		
		LambdaSnp=VarE/gvar(snpid,1)
		Lhs=XpX(snpid,1)+LambdaSnp
		Rhs=sdot(nAnisTr, GenosTr(:,snpid), 1, E(:,1), 1)
		G(snpid,1)=(Rhs/Lhs)+(gasdev(idum)*sqrt(1/Lhs))

		E(:,1)=E(:,1)-(GenosTr(:,snpid)*G(snpid,1))
		

	enddo
	if (h>nBurn) then
		GAccum=GAccum+G
		TotalRoundsContributing=TotalRoundsContributing+1
	endif	


	if (mod(h,200)==0) then
		call sgemm('n','n',nAnisTr,One,nSnp,myone,GenosTr,nAnisTr,G,nSnp,myzero,Xg,nAnisTr)
		E(:,1)=Phen(:,1)-Xg(:,1)-Mu
	endif

enddo

G=GAccum/TotalRoundsContributing

end subroutine BayesA

!###########################################################################################################################################################
