!###############################################################################

subroutine MarkerEffectPostProcessing
	use Global
	implicit none
	integer :: i,j,UnitSnpSol

	! Rescale back to phenotype scale
	G(:,1)=G(:,1)*sqrt(VarY)

	! Output
	open(newunit=UnitSnpSol,file="SnpSolutions.txt",status="unknown")
	j=0
	do i=1,nSnpExternal
		if (FixedSnp(i)==1) then
			j=j+1
			write(UnitSnpSol,"(i20,f20.10)") i,G(j,1)
		else
			write(UnitSnpSol,"(i20,f20.10)") i,0.0
		endif
	enddo
	close(UnitSnpSol)
	flush(UnitSnpSol)
end subroutine MarkerEffectPostProcessing

!###############################################################################