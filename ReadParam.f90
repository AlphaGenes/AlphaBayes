!##################################################################################################################################

subroutine ReadParam
use Global
implicit none
include 'omp_lib.h' ! for omp_get_time

character(len=100) :: dumC
integer :: i


open (unit=100,file="AlphaBayesSpec.txt",status="old")

read (100,*) dumC,GenoTrFile
read (100,*) dumC,GenoTeFile
read (100,*) dumC,PhenoTrFile
read (100,*) dumC,TbvFile
read (100,*) dumC,FileFixedSnp
read (100,*) dumC,nSnpExternal

if (trim(FileFixedSnp)/="None") then
	open (unit=10,file=trim(FileFixedSnp),status="old")
	allocate(FixedSnp(nSnpExternal))
	do i=1,nSnpExternal
		read (10,*) FixedSnp(i)
	enddo
	nSnp=sum(FixedSnp(:))
else
	nSnp=nSnpExternal
	allocate(FixedSnp(nSnpExternal))
	FixedSnp=1
endif

read (100,*) dumC,nAnisTr
read (100,*) dumC,nAnisTe
read (100,*) dumC,nRound
read (100,*) dumC,nBurn
read (100,*) dumC,VarA
read (100,*) dumC,VarE
read (100,*) dumC,nProcessors
read (100,*) dumC,ScalingOpt
read (100,*) dumC,MissingGenoCode

call OMP_SET_NUM_THREADS(nProcessors)

end subroutine ReadParam

!##################################################################################################################################
