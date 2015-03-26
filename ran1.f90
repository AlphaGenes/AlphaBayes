!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! This Function returns a uniform random deviate between 0.0 and 1.0.
! Set IDUM to any negative value to initialize or reinitialize the sequence.
!MODIFIED FOR REAL

FUNCTION ran1(idum)
IMPLICIT NONE
 INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
 DOUBLE PRECISION ran1,AM,EPS,RNMX
 PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
 INTEGER j,k,iv(NTAB),iy
 SAVE iv,iy
 DATA iv /NTAB*0/, iy /0/
  IF (idum.le.0.or.iy.eq.0) then
      idum=max(-idum,1)
  DO 11 j=NTAB+8,1,-1
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
  IF (idum.lt.0) idum=idum+IM
  IF (j.le.NTAB) iv(j)=idum

11 CONTINUE
     iy=iv(1)
  END IF
     k=idum/IQ
     idum=IA*(idum-k*IQ)-IR*k
  IF (idum.lt.0) idum=idum+IM
     j=1+iy/NDIV
     iy=iv(j)
     iv(j)=idum
     ran1=min(AM*iy,RNMX)
  RETURN
END

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
