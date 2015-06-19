!##################################################################################################################################

module Global
implicit none

integer :: idum,nSnp,nSnpExternal,nAnisTr,nAnisTe,nRound,nBurn,nProcessors,ScalingOpt,MissingGenoCode
real(4) :: VarY,VarA,VarE,Mu,Sum2pq
real(4),allocatable,dimension(:) :: SnpTmp,Lambda,SnpOut,AlleleFreq,GenoTeId
real(4),allocatable,dimension(:,:) :: GenosTr,GenosTe,Phen,E,G,Ebv,Tbv,XpX,Xg
integer,allocatable,dimension(:) :: FixedSnp,SnpPosition

character(len=1000) :: GenoTrFile,GenoTeFile,PhenoTrFile,TbvFile,FileFixedSnp,MarkerSolver

end module Global

!##################################################################################################################################
