!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!  (C) Copr. 1986-92 Numerical Recipes Software 6

 SUBROUTINE random_order(order,n,idum)
 IMPLICIT NONE

!     Generate a random ordering of the integers 1 ... n.

INTEGER, INTENT(IN)  :: n
INTEGER, INTENT(OUT) :: order(n)
INTEGER :: idum
DOUBLE PRECISION ran1

!     Local variables

INTEGER :: i, j, k
double precision    :: wk

DO i = 1, n
  order(i) = i
END DO

!     Starting at the end, swap the current last indicator with one
!     randomly chosen from those preceeding it.

DO i = n, 2, -1
  wk=ran1(idum)
  j = 1 + i * wk
  IF (j < i) THEN
    k = order(i)
    order(i) = order(j)
    order(j) = k
  END IF
END DO

RETURN
END SUBROUTINE random_order

!#############################################################################################################################################################################################################################
