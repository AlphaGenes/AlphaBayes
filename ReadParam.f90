!###############################################################################

subroutine ReadParam
	use Global
	implicit none

	integer :: i,UnitSpec,UnitFixedSnp

	character(len=100) :: DumC

	open(newunit=UnitSpec,file="AlphaBayesSpec.txt",status="old")

	read(UnitSpec,*) DumC,GenoTrFile
	read(UnitSpec,*) DumC,PhenoTrFile

	read(UnitSpec,*) DumC,nTePop
	allocate(PhenoTeFile(nTePop))
	allocate(GenoTeFile(nTePop))
	allocate(nAnisTe(nTePop))
	read(UnitSpec,*) DumC,GenoTeFile(:)
	read(UnitSpec,*) DumC,PhenoTeFile(:)

	read(UnitSpec,*) DumC,FileFixedSnp
	read(UnitSpec,*) DumC,nSnpExternal

	if (trim(FileFixedSnp)/="None") then
		open(newunit=UnitFixedSnp,file=trim(FileFixedSnp),status="old")
		allocate(FixedSnp(nSnpExternal))
		do i=1,nSnpExternal
			read(UnitFixedSnp,*) FixedSnp(i)
		enddo
		nSnp=sum(FixedSnp(:))
	else
		nSnp=nSnpExternal
		allocate(FixedSnp(nSnpExternal))
		FixedSnp=1
	endif
	close(UnitFixedSnp)

	read(UnitSpec,*) DumC,nAnisTr
	read(UnitSpec,*) DumC,nAnisTe(:)

	read(UnitSpec,*) DumC,nRound
	read(UnitSpec,*) DumC,nBurn

	read(UnitSpec,*) DumC,VarA
	read(UnitSpec,*) DumC,VarE

	read(UnitSpec,*) DumC,nProcessors
	call OMP_SET_NUM_THREADS(nProcessors)

	read(UnitSpec,*) DumC,ScalingOpt

	read(UnitSpec,*) DumC,MarkerSolver
	close(UnitSpec)
end subroutine ReadParam

!###############################################################################