!##################################################################################################################################

subroutine ReadData
	use Global
	implicit none

	integer :: i,j,nMissing,nNotMissing
	integer :: UnitGenoTr,UnitGenoTe,UnitPhenoTr,UnitPhenoTe,UnitAlleleFreq
	character(len=100) :: DumC
	real(4) :: ave,adev,sdev,var,skew,curt
	real(4) :: TmpAlleleFreq,TmpExpVarX,SumExpVarX,TmpObsVarX,SumObsVarX

	open(newunit=UnitGenoTr,file=trim(GenoTrFile),status="old")
	open(newunit=UnitPhenoTr,file=trim(PhenoTrFile),status="old")
	open(newunit=UnitGenoTe,file=trim(GenoTeFile),status="old")
	open(newunit=UnitPhenoTe,file=trim(TbvFile),status="old")
	open(newunit=UnitAlleleFreq,file="AlleleFreq.txt",status="unknown")

	allocate(SnpTmp(nSnpExternal))
	allocate(GenosTr(nAnisTr,nSnp))
	allocate(GenosTe(nAnisTe,nSnp))
	allocate(GenoTeId(nAnisTe))
	allocate(Phen(nAnisTr,1))
	allocate(E(nAnisTr,1))
	allocate(G(nSnp,1))
	allocate(Tbv(nAnisTe,1))
	allocate(Ebv(nAnisTe,1))
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
		read(UnitPhenoTr,*) DumC,Phen(i,1)
	enddo

	call momentR4(Phen(:,1),nAnisTr,ave,adev,sdev,var,skew,curt)

	Phen(:,1)=(Phen(:,1)-ave)/sdev
	VarY=var

	do i=1,nAnisTe
		read(UnitGenoTe,*) GenoTeId(i),SnpTmp(:)
		GenosTe(i,:)=SnpTmp(SnpPosition(:))
		read(UnitPhenoTe,*) DumC,Tbv(i,1)
	enddo

	SumExpVarX=0.0
	SumObsVarX=0.0
	Sum2pq=0.0
	do j=1,nSnp

		AlleleFreq(j)=0.0
	  	nMissing=0
	  	nNotMissing=0
	  	do i=1,nAnisTr
			if ((GenosTr(i,j)>-0.1).and.(GenosTr(i,j)<2.1)) then
		  		AlleleFreq(j)=AlleleFreq(j)+GenosTr(i,j)
		  		nNotMissing=nNotMissing+1
			else
		  		nMissing=nMissing+1
			endif
		enddo
		if (nNotMissing/=0) then
			AlleleFreq(j)=AlleleFreq(j)/nNotMissing
			AlleleFreq(j)=AlleleFreq(j)/2
		else
			AlleleFreq(j)=0.0
		endif

		TmpAlleleFreq=AlleleFreq(j)
		write(UnitAlleleFreq,"(i8,f11.8)") j,AlleleFreq(j)

	  	do i=1,nAnisTr
			if ((GenosTr(i,j)<-0.1).or.(GenosTr(i,j)>2.1)) then
				GenosTr(i,j)=2.0*AlleleFreq(j)
			endif
		enddo

	  	do i=1,nAnisTe
			if ((GenosTe(i,j)<-0.1).or.(GenosTe(i,j)>2.1)) then
				GenosTe(i,j)=2.0*AlleleFreq(j)
			endif
		enddo

		TmpExpVarX=2.0*(1.0-TmpAlleleFreq)*TmpAlleleFreq+0.00001
		Sum2pq=Sum2pq+TmpExpVarX
		SumExpVarX=SumExpVarX+TmpExpVarX
	  	TmpObsVarX=var+0.00001
	  	SumObsVarX=SumObsVarX+TmpObsVarX

	  	! Center
		GenosTr(:,j)=GenosTr(:,j)-(2.0*AlleleFreq(j))
		GenosTe(:,j)=GenosTe(:,j)-(2.0*AlleleFreq(j))

		if(ScalingOpt==2) then
			! Scale by marker specific variance - expected
		  	TmpExpVarX=sqrt(TmpExpVarX)
		  	GenosTr(:,j)=GenosTr(:,j)/TmpExpVarX
		  	GenosTe(:,j)=GenosTe(:,j)/TmpExpVarX
		endif

		if(ScalingOpt==3) then
			! Scale by marker specific variance - observed
		  	TmpObsVarX=sqrt(TmpObsVarX)
		  	GenosTr(:,j)=GenosTr(:,j)/TmpObsVarX
		  	GenosTe(:,j)=GenosTe(:,j)/TmpObsVarX
		endif

	enddo

	if(ScalingOpt==4) then
	  	! Scale by average marker variance - expected
	  	TmpExpVarX=sqrt(SumExpVarX/nSnp)
	  	GenosTr(:,:)=GenosTr(:,:)/TmpExpVarX
	  	GenosTe(:,:)=GenosTe(:,:)/TmpExpVarX
	endif
	if(ScalingOpt==5) then
	  	! Scale by average marker variance - observed
	  	TmpObsVarX=sqrt(SumObsVarX/nSnp)
	  	GenosTr(:,:)=GenosTr(:,:)/TmpObsVarX
	  	GenosTe(:,:)=GenosTe(:,:)/TmpObsVarX
	endif

	close(UnitGenoTr)
	close(UnitPhenoTr)
	close(UnitGenoTe)
	close(UnitPhenoTe)
	close(UnitAlleleFreq)
end subroutine ReadData

!###########################################################################################################################################################
