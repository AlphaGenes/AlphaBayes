!##################################################################################################################################

program AlphaBayes
	use Global
	implicit none

	call InitiateSeed
	call ReadParam
	call ReadData
	if (trim(MarkerSolver)=="RidgeMCMC") call RidgeRegressionMCMC
	if (trim(MarkerSolver)=="Ridge") call RidgeRegression
	if (trim(MarkerSolver)=="BayesA") call BayesA
	call MarkerEffectPostProcessing
end program AlphaBayes

!##################################################################################################################################
