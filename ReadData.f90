!###############################################################################

subroutine ReadData
	use Global
	implicit none

	integer :: i,j,nNotMissing,UnitGenoTr,UnitPhenoTr,UnitAlleleFreq

	real(4) :: ave,adev,sdev,var,skew,curt

	character(len=100) :: DumC

	open(newunit=UnitGenoTr,file=trim(GenoTrFile),status="old")
	open(newunit=UnitPhenoTr,file=trim(PhenoTrFile),status="old")
	open(newunit=UnitAlleleFreq,file="AlleleFreq.txt",status="unknown")

	allocate(SnpTmp(nSnpExternal))
	allocate(GenosTr(nAnisTr,nSnp))
	allocate(PhenTr(nAnisTr,1))
	allocate(G(nSnp,1))
	allocate(AlleleFreq(nSnp))
	allocate(SnpPosition(nSnp))

	j=0
	do i=1,nSnpExternal
		if (FixedSnp(i)==1) then
			j=j+1
			SnpPosition(j)=i
		endif
	enddo

	do i=1,nAnisTr
		read(UnitGenoTr,*) DumC,SnpTmp(:)
		GenosTr(i,:)=SnpTmp(SnpPosition(:))
		read(UnitPhenoTr,*) DumC,PhenTr(i,1)
	enddo

	call momentR4(PhenTr(:,1),nAnisTr,ave,adev,sdev,var,skew,curt)

	PhenTr(:,1)=(PhenTr(:,1)-ave)/sdev
	VarY=var

	SumExpVarX=0.0
	SumObsVarX=0.0
	Sum2pq=0.0
	do j=1,nSnp

		! Compute allele freqs
		AlleleFreq(j)=0.0
	  	nNotMissing=0
	  	do i=1,nAnisTr
			if ((GenosTr(i,j)>-0.1).and.(GenosTr(i,j)<2.1)) then
		  		AlleleFreq(j)=AlleleFreq(j)+GenosTr(i,j)
		  		nNotMissing=nNotMissing+1
			endif
		enddo
		if (nNotMissing/=0) then
			AlleleFreq(j)=AlleleFreq(j)/float(2*nNotMissing)
		else
			AlleleFreq(j)=0.0
		endif
		write(UnitAlleleFreq,"(i8,f11.8)") j,AlleleFreq(j)

		! Fix any odd data
	  	do i=1,nAnisTr
			if ((GenosTr(i,j)<-0.1).or.(GenosTr(i,j)>2.1)) then
				GenosTr(i,j)=2.0*AlleleFreq(j)
			endif
		enddo

		! Standardize
		ExpVarX=2.0*(1.0-AlleleFreq(j))*AlleleFreq(j)+0.00001
		Sum2pq=Sum2pq+ExpVarX
		SumExpVarX=SumExpVarX+ExpVarX
	  	ObsVarX=var+0.00001
	  	SumObsVarX=SumObsVarX+ObsVarX

	  	! ... center
		GenosTr(:,j)=GenosTr(:,j)-(2.0*AlleleFreq(j))

		! ... scale
		if (ScalingOpt==2) then
			! Scale by marker specific variance - expected
		  	ExpVarX=sqrt(ExpVarX)
		  	GenosTr(:,j)=GenosTr(:,j)/ExpVarX
		endif

		if (ScalingOpt==3) then
			! Scale by marker specific variance - observed
		  	ObsVarX=sqrt(ObsVarX)
		  	GenosTr(:,j)=GenosTr(:,j)/ObsVarX
		endif

	enddo

	if (ScalingOpt==4) then
	  	! Scale by average marker variance - expected
	  	ExpVarX=sqrt(SumExpVarX/float(nSnp))
	  	GenosTr(:,:)=GenosTr(:,:)/ExpVarX
	endif

	if (ScalingOpt==5) then
	  	! Scale by average marker variance - observed
	  	ObsVarX=sqrt(SumObsVarX/float(nSnp))
	  	GenosTr(:,:)=GenosTr(:,:)/ObsVarX
	endif

	close(UnitGenoTr)
	close(UnitPhenoTr)
	close(UnitAlleleFreq)
end subroutine ReadData

!###############################################################################