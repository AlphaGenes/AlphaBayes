!##################################################################################################################################

subroutine ReadParam
use Global
implicit none
include 'omp_lib.h' ! for omp_get_time

character(len=100) :: dumC
integer :: i


open (unit=100,file="AlphaBayesSpec.txt",status="old")

read (100,*) dumC,GenoTrFile
read (100,*) dumC,GenoTeFile
read (100,*) dumC,PhenoTrFile
read (100,*) dumC,TbvFile
read (100,*) dumC,FileFixedSnp
read (100,*) dumC,nSnpExternal

if (trim(FileFixedSnp)/="None") then
	open (unit=10,file=trim(FileFixedSnp),status="old")
	allocate(FixedSnp(nSnpExternal))
	do i=1,nSnpExternal
		read (10,*) FixedSnp(i)
	enddo
	nSnp=sum(FixedSnp(:))
else
	nSnp=nSnpExternal
	allocate(FixedSnp(nSnpExternal))
	FixedSnp=1
endif

read (100,*) dumC,nAnisTr
read (100,*) dumC,nAnisTe
read (100,*) dumC,nRound
read (100,*) dumC,nBurn	
read (100,*) dumC,VarA
read (100,*) dumC,VarE		
read (100,*) dumC,nProcessors

call OMP_SET_NUM_THREADS(nProcessors)

end subroutine ReadParam

!##################################################################################################################################

subroutine ReadData
use Global
implicit none

integer :: i,j
character(len=100) :: dumC
real (4) :: ave,adev,sdev,var,skew,curt

open (unit=101,file=trim(GenoTrFile),status="old")
open (unit=102,file=trim(PhenoTrFile),status="old")
open (unit=103,file=trim(GenoTeFile),status="old")
open (unit=104,file=trim(TbvFile),status="old")
open (unit=105,file="AlleleFreq.txt",status="unknown")

allocate(SnpTmp(nSnpExternal))
allocate(GenosTr(nAnisTr,nSnp))
allocate(GenosTe(nAnisTe,nSnp))
allocate(Phen(nAnisTr,1))
allocate(E(nAnisTr,1))
allocate(G(nSnp,1))
allocate(Tbv(nAnisTe,1))
allocate(Ebv(nAnisTe,1))

allocate(SnpPosition(nSnp))
j=0
do i=1,nSnpExternal
	if (FixedSnp(i)==1) then
		j=j+1
		SnpPosition(j)=i
	endif
enddo	

do i=1,nAnisTr
	read (101,*) dumC,SnpTmp(:)
	GenosTr(i,:)=SnpTmp(SnpPosition(:))
	read (102,*) dumC,Phen(i,1)
enddo	

call momentR4(Phen(:,1),nAnisTr,ave,adev,sdev,var,skew,curt)

Phen(:,1)=(Phen(:,1)-ave)/sdev

do i=1,nAnisTe
	read (103,*) dumC,SnpTmp(:)
	GenosTe(i,:)=SnpTmp(SnpPosition(:))
	read (104,*) dumC,Tbv(i,1)
enddo	

do j=1,nSnp
	call momentR4(GenosTr(:,j),nAnisTr,ave,adev,sdev,var,skew,curt)
	GenosTr(:,j)=(GenosTr(:,j)-ave)/sdev
	GenosTe(:,j)=(GenosTe(:,j)-ave)/sdev
enddo

end subroutine ReadData

!###########################################################################################################################################################

subroutine RidgeRegression
use Global
implicit none

real(4) :: sdot,eps,InvLhs,Rhs,Lhs,SolOld,myone,myzero, TmpVal,Correlation
integer :: i,h,j,snpid,RandomOrdering(nSnp),One

allocate(XpX(nSnp,1))
allocate(Xg(nAnisTr,1))
allocate(Lambda(nSnp))

myone=1.d0
myzero=0.d0
One=1

!Construct XpX
do j=1,nSnp
	XpX(j,1)=sdot(nAnisTr, GenosTr(:,j), 1, GenosTr(:,j), 1) + 0.000000000000001
	Lambda(j)=VarE/(VarA/nSnp)
enddo  

Mu=0.00
E(:,1)=Phen(:,1)-Mu !GG

do h=1,nRound
	eps=0.0	

	!Intercept
	E(:,1)=E(:,1)+Mu
	Rhs=sum(E(:,1))
	InvLhs=1.0/float(nAnisTr)
	Mu=Rhs*InvLhs
	E(:,1)=E(:,1)-Mu

	!Snp effects
	call random_order(RandomOrdering,nSnp,idum)
	do j=1,nSnp
		snpid=RandomOrdering(j)

		E(:,1)=E(:,1)+(GenosTr(:,snpid)*G(snpid,1)) 

		Lhs=XpX(snpid,1)+Lambda(snpid)
		Rhs=sdot(nAnisTr, GenosTr(:,snpid), 1, E(:,1), 1)

		SolOld=G(snpid,1)
		G(snpid,1)=Rhs/Lhs

		E(:,1)=E(:,1)-(GenosTr(:,snpid)*G(snpid,1))
		Eps=Eps+(G(snpid,1)-SolOld)**2			

	enddo

	if (mod(h,200)==0) then
		call sgemm('n','n',nAnisTr,One,nSnp,myone,GenosTr,nAnisTr,G,nSnp,myzero,Xg,nAnisTr)
		E(:,1)=Phen(:,1)-Xg(:,1)-Mu
	endif

	if (Eps.lt. 1e-10) then
		ConvergedRounds=h
		exit
	endif	

enddo

print *, "Converged in ", ConvergedRounds," rounds"

!Output section
open (unit=1002,file="SnpSolutions.txt",status="unknown")

allocate(SnpOut(nSnpExternal))

SnpOut=0.d0
j=0
do i=1,nSnpExternal
	if (FixedSnp(i)==1) then
		j=j+1
		SnpOut(i)=G(j,1) 
	endif
	write (1002,*) i,SnpOut(i)
enddo	


close(1002)

open (unit=1001,file="TbvEbv.txt",status="unknown")

call sgemm('n','n',nAnisTe,One,nSnp,myone,GenosTe,nAnisTe,G,nSnp,myzero,Ebv,nAnisTe)
do i=1,nAnisTe
	write(1001,*) i,Tbv(i,1),Ebv(i,1)
enddo
close(1001)

call PearsnR4 (Tbv(:,1),Ebv(:,1),nAnisTe,Correlation)
print*, Correlation
open (unit=1001,file="TbvEbvCorrelation.txt",status="unknown")
  write(1001,*) Correlation 
close(1001)




end subroutine RidgeRegression

!###########################################################################################################################################################

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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUNCTION gasdev(idum)
IMPLICIT NONE
!C USES ran1
!Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
!as the source of uniform deviates.

INTEGER idum
DOUBLE PRECISION :: gasdev
INTEGER iset
DOUBLE PRECISION fac,gset,rsq,v1,v2,ran1
SAVE iset,gset
DATA iset/0/
if (idum.lt.0) iset=0
if (iset.eq.0) then
1 	v1=2.*ran1(idum)-1.
	v2=2.*ran1(idum)-1.
	rsq=v1**2+v2**2
	if(rsq.ge.1..or.rsq.eq.0.)goto 1
	fac=sqrt(-2.*log(rsq)/rsq)
	gset=v1*fac
	gasdev=v2*fac
	iset=1
else
	gasdev=gset
	iset=0
endif
return
END

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

subroutine PearsnR4 (x,y,n,r)

implicit none
integer n
real(4) prob,r,z,x(n),y(n),TINY
parameter (tiny=1.e-20)
integer j
real(4) ax,ay,df,sxx,sxy,syy,t,xt,betai,yt

ax=0.0
ay=0.0
DO j=1,n
        ax=ax+x(j)
        ay=ay+y(j)
END DO
ax=ax/n
ay=ay/n
sxx=0.
syy=0.
sxy=0.
DO j=1,n
        xt=x(j)-ax
        yt=y(j)-ay
        sxx=sxx+xt**2
        syy=syy+yt**2
        sxy=sxy+xt*yt
END DO
r=sxy/(SQRT(sxx*syy)+TINY)
z=0.5*LOG(((1.+r)+TINY)/((1.-r)+TINY))
df=n-2
t=r*SQRT(df/(((1.-r)+TINY)*((1.+r)+TINY)))
!prob=betai(0.5*df,0.5,df/(df+t**2))
!prob=erfcc(ABS(z*SQRT(n-1.))/1.4142136)
prob=0
return

end subroutine PearsnR4

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
