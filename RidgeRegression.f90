!###############################################################################

subroutine RidgeRegression
	use Global
	implicit none

	integer :: i,h,j,snpid,RandomOrdering(nSnp)

	real(4) :: sdot,Eps,InvLhs,Rhs,Lhs,SolOld,OneR,ZeroR,nAnisTrR,LambdaSnp
	real(4),allocatable :: XpX(:,:),Xg(:,:),E(:,:)

	allocate(XpX(nSnp,1))
	allocate(Xg(nAnisTr,1))
	allocate(E(nAnisTr,1))

	OneR=1.0
	ZeroR=0.0
	Mu=0.0
	G=0.000001
	nAnisTrR=float(nAnisTr)

	! Construct XpX and Lambda
	do j=1,nSnp
		XpX(j,1)=sdot(nAnisTr,GenosTr(:,j),1,GenosTr(:,j),1) + 0.000000000000001
	enddo

	LambdaSnp=VarE/(VarA/float(nSnp))

	! Working phenotypes
	E(:,1)=PhenTr(:,1)

	open(unit=6,form="formatted",CARRIAGECONTROL="FORTRAN")
	do h=1,nRound
		!             123456789 123456789 123456789 123456789 123456789 12
		write(6,100) " Ridge regression (fixed variance components) round ",h
		100 format("+",a52,i10)

		!Intercept
		E(:,1)=E(:,1)+Mu
		Rhs=sum(E(:,1))
		InvLhs=1.0/nAnisTrR
		Mu=Rhs*InvLhs
		E(:,1)=E(:,1)-Mu

		!Snp effects
		eps=0.0
		call random_order(RandomOrdering,nSnp,idum)
		do j=1,nSnp
			snpid=RandomOrdering(j)
			E(:,1)=E(:,1)+(GenosTr(:,snpid)*G(snpid,1))
			Lhs=XpX(snpid,1)+LambdaSnp
			Rhs=sdot(nAnisTr,GenosTr(:,snpid),1,E(:,1),1)
			SolOld=G(snpid,1)
			G(snpid,1)=Rhs/Lhs
			E(:,1)=E(:,1)-(GenosTr(:,snpid)*G(snpid,1))
			Eps=Eps+(G(snpid,1)-SolOld)**2
		enddo

		if (mod(h,200)==0) then
			call sgemm("n","n",nAnisTr,1,nSnp,OneR,GenosTr,nAnisTr,G,nSnp,ZeroR,Xg,nAnisTr)
			E(:,1)=PhenTr(:,1)-Xg(:,1)-Mu
		endif

		if (eps.lt.1e-8) then
			exit
		endif
	enddo
	close(6)

	deallocate(XpX)
	deallocate(Xg)
	deallocate(E)
end subroutine RidgeRegression

!###############################################################################