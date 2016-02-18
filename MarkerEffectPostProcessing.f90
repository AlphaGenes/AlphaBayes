!###########################################################################################################################################################

subroutine MarkerEffectPostProcessing
	use Global
	implicit none

	real(4) :: sdot,myone,myzero,TmpVal
	real(4) :: Correlation
	integer :: i,j,One,UnitSnpSol,UnitEbv,UnitCor

	myone=1.0
	myzero=0.0
	One=1

	! Rescale back to phenotype scale
	G(:,1)=G(:,1)*sqrt(VarY)

	! Output section
	open(unit=UnitSnpSol,file="SnpSolutions.txt",status="unknown")
	allocate(SnpOut(nSnpExternal))
	SnpOut=0.0
	j=0
	do i=1,nSnpExternal
		if (FixedSnp(i)==1) then
			j=j+1
			SnpOut(i)=G(j,1)
		endif
		write(UnitSnpSol,"(i20,f20.10)") i,SnpOut(i)
	enddo
	close(UnitSnpSol)
	flush(UnitSnpSol)

	open(unit=UnitEbv,file="Ebv.txt",status="unknown")				! Before was "TbvEbv.txt"
	call sgemm("n","n",nAnisTe,One,nSnp,myone,GenosTe,nAnisTe,G,nSnp,myzero,Ebv,nAnisTe)
	do i=1,nAnisTe
		write(UnitEbv,"(i20,2f20.10)") GenoTeId(i),Ebv(i,1)			! To be directly used in AlphaDrop, the TBV are no more printed in this file. The TBV can be found
																	! in the files "SelectionPedTbvTesting.txt" and "SelectionPedTbvTestingRecomb.txt", which are used by AlphaBayes.
	enddo
	close(UnitEbv)
	flush(UnitEbv)

	open(unit=1001,file="TbvEbvCorrelation.txt",status="unknown")
	call PearsnR4 (Tbv(:,1),Ebv(:,1),nAnisTe,Correlation)
    write(UnitCor,"(f13.10)") Correlation
	close(UnitCor)
	flush(UnitCor)

	!call system("sleep 10")

end subroutine MarkerEffectPostProcessing

!###########################################################################################################################################################
