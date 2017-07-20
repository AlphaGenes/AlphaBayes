
!###############################################################################

module AlphaStatMod

  ! Basic statistical stuff

  use ISO_Fortran_Env, STDIN=>input_unit,STDOUT=>output_unit,STDERR=>error_unit

  implicit none

  private
  ! Types
  public :: DescStatS,DescStatD
  public :: DescStatMatrixS,DescStatMatrixD
  public :: CorrelationS,CorrelationD
  ! Methods
  public :: CalcMean,CalcVar,CalcSD
  public :: CalcDescStat
  public :: CalcDescStatMatrix,CalcDescStatSymMatrix,CalcDescStatLowTriMatrix
  public :: CalcCorrelation

  interface CalcMean
    module procedure CalcMeanS, CalcMeanD
  end interface

  interface CalcVar
    module procedure CalcVarS, CalcVarD
  end interface

  interface CalcSD
    module procedure CalcSDS, CalcSDD
  end interface

  interface CalcDescStat
    module procedure CalcDescStatS, CalcDescStatD
  end interface

  type :: DescStatS
    integer(int32) :: n
    real(real32)   :: Mean
    real(real32)   :: Var
    real(real32)   :: SD
    real(real32)   :: Skew ! TODO: do we ever need skewness?
    real(real32)   :: Curt ! TODO: do we ever need curtosis?
    real(real32)   :: Min
    real(real32)   :: Max
  end type

  type :: DescStatD
    integer(int32) :: n
    real(real64)   :: Mean
    real(real64)   :: Var
    real(real64)   :: SD
    real(real64)   :: Skew
    real(real64)   :: Curt
    real(real64)   :: Min
    real(real64)   :: Max
  end type

  interface CalcDescStatSymMatrix
    module procedure CalcDescStatSymMatrixS, CalcDescStatSymMatrixD
  end interface

  interface CalcDescStatLowTriMatrix
    module procedure CalcDescStatSymMatrixS, CalcDescStatSymMatrixD
  end interface

  interface CalcDescStatMatrix
    module procedure CalcDescStatMatrixS, CalcDescStatMatrixD
  end interface

  type :: DescStatMatrixS
    type(DescStatS) :: All
    type(DescStatS) :: Diag
    type(DescStatS) :: OffDiag
  end type

  type :: DescStatMatrixD
    type(DescStatD) :: All
    type(DescStatD) :: Diag
    type(DescStatD) :: OffDiag
  end type

  interface CalcCorrelation
    module procedure CalcCorrelationS, CalcCorrelationD
  end interface

  type :: CorrelationS
    real(real32)   :: Cor
    real(real32)   :: Cov
    real(real32)   :: Var1
    real(real32)   :: Var2
  end type

  type :: CorrelationD
    real(real64)   :: Cor
    real(real64)   :: Cov
    real(real64)   :: Var1
    real(real64)   :: Var2
  end type

  contains

    !###########################################################################

    function CalcMeanS(x) result(Res)
      implicit none

      ! Arguments
      real(real32),intent(in)  :: x(:)
      real(real32)             :: Res

      Res=sum(x(:))/real(size(x))
      return
    end function

    !###########################################################################

    function CalcMeanD(x) result(Res)
      implicit none

      ! Arguments
      real(real64),intent(in)  :: x(:)
      real(real64)             :: Res

      Res=sum(x(:))/dble(size(x))
      return
    end function

    !###########################################################################

    function CalcVarS(x,Mean) result(Res)
      implicit none

      ! Arguments
      real(real32),intent(in)          :: x(:)
      real(real32),intent(in),optional :: Mean
      real(real32)                     :: Res
      ! Other
      integer(int32) i,n

      real(real32) :: nR,MeanIn,Dev

      n=size(x)
      if (n < 1) then
        write(STDERR,"(a)") "ERROR: number of records must be at least 2 for CalcVar"
        write(STDERR,"(a)") " "
        stop 1
      end if

      nR=real(n)
      if (.not.present(Mean)) then
        MeanIn=sum(x(:))/nR
      else
        MeanIn=Mean
      end if
      Res=0.0

      do i=1,n
        Dev=x(i)-MeanIn
        Res=Res+Dev*Dev
      end do

      Res=Res/(nR-1.0)
      return
    end function

    !###########################################################################

    function CalcVarD(x,Mean) result(Res)
      implicit none

      ! Arguments
      real(real64),intent(in)          :: x(:)
      real(real64),intent(in),optional :: Mean
      real(real64)                     :: Res

      ! Other
      integer(int32) i,n

      real(real64) :: nR,MeanIn,Dev

      n=size(x)
      if (n < 1) then
        write(STDERR,"(a)") "ERROR: number of records must be at least 2 for CalcVar"
        write(STDERR,"(a)") " "
        stop 1
      end if

      nR=dble(n)
      if (.not.present(Mean)) then
        MeanIn=sum(x(:))/nR
      else
        MeanIn=Mean
      end if
      Res=0.0d0

      do i=1,n
        Dev=x(i)-MeanIn
        Res=Res+Dev*Dev
      end do

      Res=Res/(nR-1.0d0)
      return
    end function

    !###########################################################################

    function CalcSDS(x,Mean) result(Res)
      implicit none

      ! Arguments
      real(real32),intent(in)          :: x(:)
      real(real32),intent(in),optional :: Mean
      real(real32)                     :: Res

      Res=sqrt(CalcVarS(x,Mean))
      return
    end function

    !###########################################################################

    function CalcSDD(x,Mean) result(Res)
      implicit none

      ! Arguments
      real(real64),intent(in)          :: x(:)
      real(real64),intent(in),optional :: Mean
      real(real64)                     :: Res

      Res=sqrt(CalcVarD(x,Mean))
      return
    end function

    !###########################################################################

    function CalcDescStatS(x) result(Res)
      implicit none

      ! Arguments
      real(real32),intent(in) :: x(:)
      type(DescStatS)         :: Res

      ! Other
      integer(int32) i,n
      real(real32) :: nR,Dev,Mul

      n=size(x)
      if (n < 1) then
        write(STDERR,"(a)") "ERROR: number of records must be at least 2 for CalcDescStat"
        write(STDERR,"(a)") " "
        stop 1
      end if

      nR=real(n)
      Res%Mean=sum(x(:))/nR
      Res%Var=0.0
      Res%Skew=0.0
      Res%Curt=0.0

      do i=1,n
        Dev=x(i)-Res%Mean
        Mul=Dev*Dev
        Res%Var=Res%Var+Mul
        Mul=Mul*Dev
        Res%Skew=Res%Skew+Mul
        Mul=Mul*Dev
        Res%Curt=Res%Curt+Mul
      end do

      Res%Var=Res%Var/(nR-1.0)
      Res%SD=sqrt(Res%Var)
      if (Res%Var > 0.0) then
        Res%Skew=Res%Skew/(nR*Res%SD**3)
        Res%Curt=Res%Curt/(nR*Res%Var**2)-3.0
      end if

      Res%Min=minval(x(:))
      Res%Max=maxval(x(:))
      return
    end function

    !###########################################################################

    function CalcDescStatD(x) result(Res)
      implicit none

      ! Arguments
      real(real64),intent(in) :: x(:)
      type(DescStatD)         :: Res

      ! Other
      integer(int32) i,n

      real(real64) :: nR,Dev,Mul

      n=size(x)
      if (n < 1) then
        write(STDERR,"(a)") "ERROR: number of records must be at least 2 for CalcDescStat"
        write(STDERR,"(a)") " "
        stop 1
      end if

      nR=dble(n)
      Res%Mean=sum(x(:))/nR
      Res%Var=0.0d0
      Res%Skew=0.0d0
      Res%Curt=0.0d0

      do i=1,n
        Dev=x(i)-Res%Mean
        Mul=Dev*Dev
        Res%Var=Res%Var+Mul
        Mul=Mul*Dev
        Res%Skew=Res%Skew+Mul
        Mul=Mul*Dev
        Res%Curt=Res%Curt+Mul
      end do

      Res%Var=Res%Var/(nR-1.0d0)
      Res%SD=sqrt(Res%Var)
      if (Res%Var > 0.0d0) then
        Res%Skew=Res%Skew/(nR*Res%SD**3)
        Res%Curt=Res%Curt/(nR*Res%Var**2)-3.0d0
      end if

      Res%Min=minval(x(:))
      Res%Max=maxval(x(:))
      return
    end function

    !###########################################################################

    function CalcDescStatSymMatrixS(x,Diag) result(Res)
      implicit none

      ! Arguments
      real(real32),intent(in)     :: x(:,:)
      logical,intent(in),optional :: Diag
      type(DescStatMatrixS)       :: Res

      ! Other
      integer(int32) i,j,k,n,p
      real(real32),allocatable :: DiagVal(:)
      real(real32),allocatable :: OffDiagVal(:)
      logical                  :: DiagInternal

      n=size(x,1)
      p=size(x,2)
      if (n /= p) then
        write(STDERR,"(a)") "ERROR: CalcDescStatSymMatrix work only with symmetric matrices!"
        write(STDERR,"(a)") " "
        stop 1
      end if

      if (present(Diag)) then
        DiagInternal=Diag
      else
        DiagInternal=.true.
      end if

      ! Diagonal
      if (DiagInternal) then
        allocate(DiagVal(n))
        do i=1,n
          DiagVal(i)=x(i,i)
        end do
        Res%Diag=CalcDescStatS(DiagVal)
      end if

      ! Off-diagonal (lower-triangle only!!!)
      allocate(OffDiagVal(nint(real(n*n)/2.0-real(n)/2.0))) ! n*n/2 is half of a matrix, n/2 removes half of diagonal
      k=0
      do j=1,(p-1)
        do i=(j+1),n
          k=k+1
          OffDiagVal(k)=x(i,j)
        end do
      end do
      Res%OffDiag=CalcDescStatS(OffDiagVal)

      if (DiagInternal) then
        Res%All=CalcDescStatS([DiagVal,OffDiagVal])
      end if
      return
    end function

    !###########################################################################

    function CalcDescStatSymMatrixD(x,Diag) result(Res)
      implicit none

      ! Arguments
      real(real64),intent(in)     :: x(:,:)
      logical,intent(in),optional :: Diag
      type(DescStatMatrixD)       :: Res

      ! Other
      integer(int32) i,j,k,n,p
      real(real64),allocatable :: DiagVal(:)
      real(real64),allocatable :: OffDiagVal(:)
      logical                  :: DiagInternal

      n=size(x,1)
      p=size(x,2)
      if (n /= p) then
        write(STDERR,"(a)") "ERROR: CalcDescStatSymMatrix work only with symmetric matrices!"
        write(STDERR,"(a)") " "
        stop 1
      end if

      if (present(Diag)) then
        DiagInternal=Diag
      else
        DiagInternal=.true.
      end if

      ! Diagonal
      if (DiagInternal) then
        allocate(DiagVal(n))
        do i=1,n
          DiagVal(i)=x(i,i)
        end do
        Res%Diag=CalcDescStatD(DiagVal)
      end if

      ! Off-diagonal (lower-triangle only!!!)
      allocate(OffDiagVal(nint(real(n*n)/2.0-real(n)/2.0))) ! n*n/2 is half of a matrix, n/2 removes half of diagonal
      k=0
      do j=1,(p-1)
        do i=(j+1),n
          k=k+1
          OffDiagVal(k)=x(i,j)
        end do
      end do
      Res%OffDiag=CalcDescStatD(OffDiagVal)

      if (DiagInternal) then
        Res%All=CalcDescStatD([DiagVal,OffDiagVal])
      end if
      return
    end function

    !###########################################################################

    function CalcDescStatMatrixS(x) result(Res)
      implicit none

      ! Arguments
      real(real32),intent(in) :: x(:,:)
      type(DescStatMatrixS)   :: Res

      ! Other
      integer(int32) i,j,k,l,n,p,MinNP
      real(real32),allocatable :: Diag(:)
      real(real32),allocatable :: OffDiag(:)

      n=size(x,1)
      p=size(x,2)

      MinNP=minval([n,p])
      allocate(Diag(MinNP))
      allocate(OffDiag(n*p-MinNP))

      k=0
      l=0
      do j=1,p
        do i=1,n
          if (i == j) then
            k=k+1
            Diag(k)=x(i,j)
          else
            l=l+1
            OffDiag(l)=x(i,j)
          end if
        end do
      end do

      Res%Diag=CalcDescStatS(Diag)
      Res%OffDiag=CalcDescStatS(OffDiag)
      Res%All=CalcDescStatS([Diag,OffDiag])
      return
    end function

    !###########################################################################

    function CalcDescStatMatrixD(x) result(Res)
      implicit none

      ! Arguments
      real(real64),intent(in) :: x(:,:)
      type(DescStatMatrixD)   :: Res

      ! Other
      integer(int32) i,j,k,l,n,p,MinNP
      real(real64),allocatable :: Diag(:)
      real(real64),allocatable :: OffDiag(:)

      n=size(x,1)
      p=size(x,2)

      MinNP=minval([n,p])
      allocate(Diag(MinNP))
      allocate(OffDiag(n*p-MinNP))

      k=0
      l=0
      do j=1,p
        do i=1,n
          if (i == j) then
            k=k+1
            Diag(k)=x(i,j)
          else
            l=l+1
            OffDiag(l)=x(i,j)
          end if
        end do
      end do

      Res%Diag=CalcDescStatD(Diag)
      Res%OffDiag=CalcDescStatD(OffDiag)
      Res%All=CalcDescStatD([Diag,OffDiag])
      return
    end function

    !###########################################################################

    function CalcCorrelationS(x,y) result(Res)
      implicit none
      ! Arguments
      real(real32),intent(in) :: x(:)
      real(real32),intent(in) :: y(:)
      type(CorrelationS)      :: Res

      ! Other
      integer(int32) :: i,n
      real(real32) :: nR,MeanX,MeanY,DevX,DevY,SumXX,SumXY,SumYY

      n=size(x)
      nR=real(n)
      MeanX=sum(x(:))/nR
      MeanY=sum(y(:))/nR
      SumXX=0.0
      SumXY=0.0
      SumYY=0.0
      do i=1,n
        DevX=x(i)-MeanX
        DevY=y(i)-MeanY
        SumXX=SumXX+DevX*DevX
        SumXY=SumXY+DevX*DevY
        SumYY=SumYY+DevY*DevY
      end do
      Res%Cor=SumXY/(sqrt(SumXX*SumYY)+tiny(x))
      Res%Cov=SumXY/(nR-1.0)
      Res%Var1=SumXX/(nR-1.0)
      Res%Var2=SumYY/(nR-1.0)
      return
    end function

    !###########################################################################

    function CalcCorrelationD(x,y) result(Res)
      implicit none
      ! Arguments
      real(real64),intent(in) :: x(:)
      real(real64),intent(in) :: y(:)
      type(CorrelationD)      :: Res

      ! Other
      integer(int32) :: i,n
      real(real64) :: nR,MeanX,MeanY,DevX,DevY,SumXX,SumXY,SumYY

      n=size(x)
      nR=dble(n)
      MeanX=sum(x(:))/nR
      MeanY=sum(y(:))/nR
      SumXX=0.0d0
      SumXY=0.0d0
      SumYY=0.0d0
      do i=1,n
        DevX=x(i)-MeanX
        DevY=y(i)-MeanY
        SumXX=SumXX+DevX*DevX
        SumXY=SumXY+DevX*DevY
        SumYY=SumYY+DevY*DevY
      end do
      Res%Cor=SumXY/(sqrt(SumXX*SumYY)+tiny(x))
      Res%Cov=SumXY/(nR-1.0d0)
      Res%Var1=SumXX/(nR-1.0d0)
      Res%Var2=SumYY/(nR-1.0d0)
      return
    end function

    !###########################################################################
end module

!###############################################################################
