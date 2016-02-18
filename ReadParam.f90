!##################################################################################################################################

subroutine ReadParam
	use Global
	implicit none
	!include 'omp_lib.h' ! for omp_get_time

	character(len=100) :: DumC
	integer :: i,UnitSpec,UnitFixedSnp

	open(newunit=UnitSpec,file="AlphaBayesSpec.txt",status="old")

	read(UnitSpec,*) DumC,GenoTrFile
	read(UnitSpec,*) DumC,GenoTeFile
	read(UnitSpec,*) DumC,PhenoTrFile
	read(UnitSpec,*) DumC,PhenoTeFile
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

	read(UnitSpec,*) DumC,nAnisTr
	read(UnitSpec,*) DumC,nAnisTe
	read(UnitSpec,*) DumC,nRound
	read(UnitSpec,*) DumC,nBurn
	read(UnitSpec,*) DumC,VarA
	read(UnitSpec,*) DumC,VarE
	read(UnitSpec,*) DumC,nProcessors
	read(UnitSpec,*) DumC,ScalingOpt
	read(UnitSpec,*) DumC,MarkerSolver

	call OMP_SET_NUM_THREADS(nProcessors)
end subroutine ReadParam

!##################################################################################################################################
