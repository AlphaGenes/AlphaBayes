!###############################################################################

module Global
	implicit none

	integer :: idum,nSnp,nSnpExternal,nAnisTr,nRound,nBurn,nProcessors,ScalingOpt,nTePop
	integer,allocatable :: nAnisTe(:),FixedSnp(:),SnpPosition(:)

	real(4) :: VarY,VarA,VarE,Mu,Sum2pq
	real(4) :: ExpVarX,SumExpVarX,ObsVarX,SumObsVarX
	real(4),allocatable :: SnpTmp(:),AlleleFreq(:),GenosTr(:,:),PhenTr(:,:),G(:,:)

	character(len=1000) :: GenoTrFile,PhenoTrFile,FileFixedSnp,MarkerSolver
	character(len=1000),allocatable :: GenoTeFile(:),PhenoTeFile(:)
end module Global

!###############################################################################