!###########################################################################################################################################################

subroutine RidgeRegression
use Global
implicit none

real(4) :: sdot,Eps,InvLhs,Rhs,Lhs,SolOld,myone,myzero
integer :: i,h,j,snpid,RandomOrdering(nSnp),One

allocate(XpX(nSnp,1))
allocate(Xg(nAnisTr,1))
allocate(Lambda(nSnp))

myone=1.0
myzero=0.0
One=1
Mu=0.0
G=0.000001

!Construct XpX
do j=1,nSnp
	XpX(j,1)=sdot(nAnisTr, GenosTr(:,j), 1, GenosTr(:,j), 1) + 0.000000000000001
	Lambda(j)=VarE/(VarA/nSnp)
enddo

E(:,1)=Phen(:,1)

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

	if (eps.lt.1e-8) then
		exit
	endif
enddo

print *, "Stopped after ", h," rounds out of ", nRound

end subroutine RidgeRegression

!###########################################################################################################################################################
