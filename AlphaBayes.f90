!##################################################################################################################################

program AlphaBayes
use Global
implicit none


call InitiateSeed
call ReadParam
call ReadData
call RidgeRegression

end program AlphaBayes

!##################################################################################################################################
