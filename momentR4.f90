!###########################################################################################################################################################

SUBROUTINE momentR4(DATA,n,ave,adev,sdev,var,skew,curt)
IMPLICIT NONE
INTEGER n
real(4) :: adev,ave,curt,sdev,skew,var,DATA(n)
INTEGER j
real(4) :: p,s,ep
IF (n.le.1) PAUSE 'n must be at least 2 in moment'
s=0
DO j= 1,n
        s=s+DATA(j)
END DO

ave=s/n
adev=0
var=0
skew=0
curt=0
ep=0

DO j=1,n
        s=DATA(j)-ave
        ep=ep+s
        adev=adev+ABS(s)
        p=s*s
        var=var+p
        p=p*s
        skew=skew+p
        p=p*s
        curt=curt+p
END DO

adev=adev/n
var=(var-ep**2/n)/(n-1)
sdev=SQRT(var)
IF(var.ne.0)then
        skew=skew/(n*sdev**3)
        curt=curt/(n*var**2)-3
ELSE
        !PRINT*, 'no skew or kurtosis when zero variance in moment'
        !PAUSE 'no skew or kurtosis when zero variance in moment'
END IF
RETURN
END SUBROUTINE momentR4

!#############################################################################################################################################################################################################################
