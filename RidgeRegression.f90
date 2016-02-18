!###########################################################################################################################################################

subroutine RidgeRegression
	use Global
	implicit none

	real(4) :: sdot,Eps,InvLhs,Rhs,Lhs,SolOld,myone,myzero,nAnisTrR
	integer :: i,h,j,snpid,RandomOrdering(nSnp),One

	allocate(XpX(nSnp,1))
	allocate(Xg(nAnisTr,1))
	allocate(Lambda(nSnp))

	myone=1.0
	myzero=0.0
	One=1
	Mu=0.0
	G=0.000001
	nAnisTrR=float(nAnisTr)

	! Construct XpX and Lambda
	do j=1,nSnp
		XpX(j,1)=sdot(nAnisTr,GenosTr(:,j),One,GenosTr(:,j),One) + 0.000000000000001
		Lambda(j)=VarE/(VarA/float(nSnp))
	enddo

	! Working phenotypes
	E(:,1)=Phen(:,1)

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
			Lhs=XpX(snpid,1)+Lambda(snpid)
			Rhs=sdot(nAnisTr,GenosTr(:,snpid),One,E(:,1),One)
			SolOld=G(snpid,1)
			G(snpid,1)=Rhs/Lhs
			E(:,1)=E(:,1)-(GenosTr(:,snpid)*G(snpid,1))
			Eps=Eps+(G(snpid,1)-SolOld)**2
		enddo

		if (mod(h,200)==0) then
			call sgemm("n","n",nAnisTr,One,nSnp,myone,GenosTr,nAnisTr,G,nSnp,myzero,Xg,nAnisTr)
			E(:,1)=Phen(:,1)-Xg(:,1)-Mu
		endif

		if (eps.lt.1e-8) then
			exit
		endif
	enddo
	close(6)

end subroutine RidgeRegression

!###########################################################################################################################################################
