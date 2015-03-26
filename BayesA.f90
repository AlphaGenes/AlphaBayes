!###########################################################################################################################################################

subroutine BayesA
use Global
implicit none

real(4) :: sdot,InvLhs,Rhs,Lhs,SolOld,myone,myzero,LambdaSnp,EpE
real(4),allocatable,dimension(:,:) :: GAccum,Gvar
real(8) :: random_gamma,gasdev,vdf,ShapeCoefficient,ScaleCoefficient,sc,gampe1,gamp2
integer :: i,h,j,snpid,RandomOrdering(nSnp),One,nSamp

allocate(XpX(nSnp,1))
allocate(Xg(nAnisTr,1))
allocate(Lambda(nSnp))
allocate(Gvar(nSnp,1))
allocate(GAccum(nSnp,1))

myone=1.0
myzero=0.0
One=1
Mu=0.0
G=0.1
GAccum=0

vdf=4.012
ShapeCoefficient=(vdf/2) !Ben
ScaleCoefficient=((vdf-2)*VarA)/Sum2pq
sc=2.0
gampe1=float((nAnisTr)/2)-2.0
gamp2=2.0

!Construct XpX
do j=1,nSnp
	XpX(j,1)=sdot(nAnisTr, GenosTr(:,j), 1, GenosTr(:,j), 1) + 0.000000000000001
enddo

E(:,1)=Phen(:,1)

open (unit=6,form='formatted',CARRIAGECONTROL='FORTRAN') 

nSamp=nRound-nBurn

do h=1,nRound

	EpE=sdot(nAnisTr, E(:,1), 1, E(:,1), 1) + 0.000000000000001
	VarE=EpE/real(random_gamma(idum,gampe1,gamp2,.TRUE.))

	write(6, 100) "   BayesA round   ",h
	100 format ('+', a17,i10)    

	do j=1,nSnp
		GVar(j,1)=(ScaleCoefficient+(G(j,1)**2))/real(random_gamma(idum,ShapeCoefficient,sc,.TRUE.))
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
		GAccum=GAccum+G/nSamp
	endif	


	if (mod(h,200)==0) then
		call sgemm('n','n',nAnisTr,One,nSnp,myone,GenosTr,nAnisTr,G,nSnp,myzero,Xg,nAnisTr)
		E(:,1)=Phen(:,1)-Xg(:,1)-Mu
	endif

enddo

G=GAccum

end subroutine BayesA

!###########################################################################################################################################################
