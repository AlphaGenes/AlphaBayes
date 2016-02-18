!###########################################################################################################################################################

subroutine RidgeRegressionMCMC
	use Global
	implicit none

	real(4) :: sdot,InvLhs,Rhs,Lhs,SolOld,myone,myzero,nAnisTrR,LambdaSnp,EpE,VarG,GpG
	real(4) :: VarEAccum,VarGAccum,nSampR
	real(4),allocatable :: GAccum(:,:),SX2(:),MX2(:)
	real(8) :: random_gamma,gasdev,R2,EDF0,EDF2,GDF0,GDF2,ES0,GS0,MSX
	integer :: i,h,j,snpid,RandomOrdering(nSnp),One

	allocate(XpX(nSnp,1))
	allocate(Xg(nAnisTr,1))
	allocate(GAccum(nSnp,1))
	allocate(SX2(nSnp))
	allocate(MX2(nSnp))

	myone=1.0
	myzero=0.0
	One=1
	Mu=0.0
	G=0.1
	GAccum=0.0

	nAnisTrR=float(nAnisTr)
	R2=0.5

	EDF0=5.0
	EDF2=(float(nAnisTr)+EDF0)/2.0
	ES0=(VarY*(1.0-R2))*(EDF0+2.0)

	GDF0=5.0
	GDF2=(float(nSnp)+GDF0)/2.0
	SX2(:)=0.0
	MX2(:)=0.0
	do i=1,nSnp
		SX2(i)=sum(GenosTr(:,i)*GenosTr(:,i))
		MX2(i)=sum(GenosTr(:,i))/nAnisTrR
		MX2(i)=MX2(i)*MX2(i)
	enddo
	MSX=sum(SX2)/nAnisTrR-sum(MX2)
	GS0=((VarY*R2)/MSX)*(GDF0+2.0)

	! Construct XpX
	do j=1,nSnp
		XpX(j,1)=sdot(nAnisTr,GenosTr(:,j),One,GenosTr(:,j),One) + 0.000000000000001
	enddo

	! Working phenotypes
	E(:,1)=Phen(:,1)

	open(unit=6,form="formatted",CARRIAGECONTROL="FORTRAN")
	nSampR=float(nRound-nBurn)
	do h=1,nRound
		!             123456789 123456789 123456789 123456789 123456789 123456
		write(6,100) " Ridge regression (MCMC for all model parameters) round ",h
		100 format("+", a56,i10)

		! Residual variance
		EpE=sdot(nAnisTr,E(:,1),1,E(:,1),1) + ES0 + 0.000000000000001
		VarE=EpE/real(random_gamma(idum,EDF2,2.0d0,.true.))

		! Snp variance
		GpG=sdot(nSnp,G(:,1),1,G(:,1),1) + GS0 + 0.000000000000001
		VarG=GpG/real(random_gamma(idum,GDF2,2.0d0,.true.))

		! Ratio of variances
		LambdaSnp=VarE/VarG

		!Intercept
		E(:,1)=E(:,1)+Mu
		Rhs=sum(E(:,1))
		InvLhs=1.0/nAnisTrR
		Mu=Rhs*InvLhs
		E(:,1)=E(:,1)-Mu

		!Snp effects
		call random_order(RandomOrdering,nSnp,idum)
		do j=1,nSnp
			snpid=RandomOrdering(j)
			E(:,1)=E(:,1)+(GenosTr(:,snpid)*G(snpid,1))
			Lhs=XpX(snpid,1)+LambdaSnp
			Rhs=sdot(nAnisTr,GenosTr(:,snpid),One,E(:,1),One)
			G(snpid,1)=(Rhs/Lhs)+(gasdev(idum)*sqrt(1.0/Lhs))
			E(:,1)=E(:,1)-(GenosTr(:,snpid)*G(snpid,1))
		enddo
		if (h>nBurn) then
			GAccum=GAccum+G/nSampR
			VarGAccum=VarGAccum+VarE/nSampR
			VarEAccum=VarEAccum+VarG/nSampR
		endif

		! Recompute residuals to avoid rounding errors
		if (mod(h,200)==0) then
			call sgemm("n","n",nAnisTr,One,nSnp,myone,GenosTr,nAnisTr,G,nSnp,myzero,Xg,nAnisTr)
			E(:,1)=Phen(:,1)-Xg(:,1)-Mu
		endif

	enddo
	close(6)

	G=GAccum
	VarG=VarGAccum
	VarE=VarEAccum

end subroutine RidgeRegressionMCMC

!###########################################################################################################################################################
