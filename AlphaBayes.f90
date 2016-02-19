#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

!###############################################################################

program AlphaBayes
	use Global
	implicit none

    write(stdout,"(a)") ""
    write(stdout,"(a30,a,a30)") " ","**********************"," "
    write(stdout,"(a30,a,a30)") " ","*                    *"," "
    write(stdout,"(a30,a,a30)") " ","*    AlphaBayes      *"," "
    write(stdout,"(a30,a,a30)") " ","*                    *"," "
    write(stdout,"(a30,a,a30)") " ","**********************"
    write(stdout,"(a30,a,a30)") " ","VERSION:"//TOSTRING(VERS)," "
    write(stdout,"(a)") ""
    write(stdout,"(a35,a)")     " ","No Liability"
    write(stdout,"(a25,a)")     " ","Bugs to John.Hickey@roslin.ed.ac.uk"
    write(stdout,"(a)") ""

	call InitiateSeed
	call ReadParam
	call ReadData
	if (trim(MarkerSolver)=="RidgeMCMC") call RidgeRegressionMCMC
	if (trim(MarkerSolver)=="Ridge") call RidgeRegression
	if (trim(MarkerSolver)=="BayesA") call BayesA
	call MarkerEffectPostProcessing
	call Prediction
end program AlphaBayes

!###############################################################################