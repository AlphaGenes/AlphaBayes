!###############################################################################

subroutine Prediction
	use Global
	implicit none

	integer :: i,j,Pop,UnitCor,UnitGenoTe,UnitPhenoTe,UnitEbv
	integer,allocatable :: IdTe(:)

	real(4) :: sdot,OneR,ZeroR,Correlation
	real(4),allocatable :: Ebv(:,:),PhenoTe(:,:),GenosTe(:,:)

	character(len=100) :: DumC

	OneR=1.0
	ZeroR=0.0

	open(newunit=UnitCor,file="TbvEbvCorrelation.txt",status="unknown")
	do Pop=1,nTePop
		allocate(IdTe(nAnisTe(Pop)))
		allocate(GenosTe(nAnisTe(Pop),nSnp))
		allocate(PhenoTe(nAnisTe(Pop),1))
		allocate(Ebv(nAnisTe(Pop),1))

		open(newunit=UnitGenoTe,file=trim(GenoTeFile(Pop)),status="old")
		open(newunit=UnitPhenoTe,file=trim(PhenoTeFile(Pop)),status="old")

		do i=1,nAnisTe(Pop)
			read(UnitGenoTe,*) IdTe(i),SnpTmp(:)
			GenosTe(i,:)=SnpTmp(SnpPosition(:))
			read(UnitPhenoTe,*) DumC,PhenoTe(i,1)
		enddo

		do j=1,nSnp
			! Fix any odd data
		  	do i=1,nAnisTe(Pop)
				if ((GenosTe(i,j)<-0.1).or.(GenosTe(i,j)>2.1)) then
					GenosTe(i,j)=2.0*AlleleFreq(j)
				endif
			enddo

			! Standardize

			! ... center
			GenosTe(:,j)=GenosTe(:,j)-(2.0*AlleleFreq(j))

			! ... scale
			if (ScalingOpt==2) then
				! Scale by marker specific variance - expected
			  	GenosTe(:,j)=GenosTe(:,j)/ExpVarX
			endif

			if (ScalingOpt==3) then
				! Scale by marker specific variance - observed
			  	GenosTe(:,j)=GenosTe(:,j)/ObsVarX
			endif
		enddo

		if (ScalingOpt==4) then
		  	! Scale by average marker variance - expected
		  	GenosTe(:,:)=GenosTe(:,:)/ExpVarX
		endif

		if (ScalingOpt==5) then
		  	! Scale by average marker variance - observed
		  	GenosTe(:,:)=GenosTe(:,:)/ObsVarX
		endif

		close(UnitGenoTe)
		close(UnitPhenoTe)

		open(newunit=UnitEbv,file="Ebv.txt",status="unknown")
		call sgemm("n","n",nAnisTe(Pop),1,nSnp,OneR,GenosTe,nAnisTe(Pop),G,nSnp,ZeroR,Ebv,nAnisTe(Pop))
		do i=1,nAnisTe(Pop)
			write(UnitEbv,"(i20,2f20.10)") IdTe(i),Ebv(i,1)
		enddo
		close(UnitEbv)
		flush(UnitEbv)

		call PearsnR4(PhenoTe(:,1),Ebv(:,1),nAnisTe(Pop),Correlation)
	    write(UnitCor,"(f13.10)") Pop,Correlation

		close(UnitGenoTe)
		close(UnitPhenoTe)
		deallocate(IdTe)
		deallocate(GenosTe)
		deallocate(PhenoTe)
		deallocate(Ebv)
	enddo
	close(UnitCor)
	flush(UnitCor)
end subroutine Prediction

!###############################################################################