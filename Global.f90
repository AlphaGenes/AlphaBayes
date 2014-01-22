!##################################################################################################################################

module Global
implicit none

integer :: idum,nSnp,nSnpExternal,nAnisTr,nAnisTe,nRound,nBurn,nProcessors,ConvergedRounds,ScalingOpt,MissingGenoCode
real(4) :: VarA,VarE,Mu
real(4),allocatable,dimension(:) :: SnpTmp,Lambda,SnpOut
real(4),allocatable,dimension(:,:) :: GenosTr,GenosTe,Phen,E,G,Ebv,Tbv,XpX,Xg
integer,allocatable,dimension(:) :: FixedSnp,SnpPosition


character(len=1000) :: GenoTrFile,GenoTeFile,PhenoTrFile,TbvFile,FileFixedSnp

end module Global

!##################################################################################################################################
