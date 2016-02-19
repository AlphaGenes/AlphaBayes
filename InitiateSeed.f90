!###############################################################################

subroutine InitiateSeed
use Global
implicit none
integer :: edum
DOUBLE PRECISION :: W(1),GASDEV

open (unit=3,file="Seed.txt",status="old")

!READ AND WRITE SEED BACK TO FILE
READ (3,*) idum
W(1)=GASDEV(idum)
!Code to write new seed to file
IF (idum>=0) THEN
	edum=(-1*idum)
ELSE
	edum=idum
END IF
REWIND (3)
WRITE (3,*) edum
idum=edum

end subroutine InitiateSeed

!###############################################################################