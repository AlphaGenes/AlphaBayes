!###########################################################################################################################################################

SUBROUTINE momentR4(DATA,n,ave,adev,sdev,var,skew,curt)
    IMPLICIT NONE
    INTEGER n,j
    real(4) :: adev,ave,curt,sdev,skew,var,DATA(n)
    real(4) :: p,s,ep
    IF (n.le.1) PAUSE 'n must be at least 2 in moment'
    s=0.0
    DO j=1,n
        s=s+DATA(j)
    END DO

    ave=s/float(n)
    adev=0.0
    var=0.0
    skew=0.0
    curt=0.0
    ep=0.0

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

    adev=adev/float(n)
    var=(var-ep**2.0/float(n))/(float(n)-1.0)
    sdev=SQRT(var)
    IF (var.ne.0) then
        skew=skew/(float(n)*sdev**3.0)
        curt=curt/(float(n)*var**2.0)-3.0
    ELSE
        !PRINT*, 'no skew or kurtosis when zero variance in moment'
        !PAUSE 'no skew or kurtosis when zero variance in moment'
    END IF
    RETURN
END SUBROUTINE momentR4

!#############################################################################################################################################################################################################################
