
!###############################################################################

function SampleIntelMultinomialI(n,p) result(Res)

  ! Sample n values from a Multinomial(p) distribution
  ! n input (integer), number of samples to generate (default 1)
  ! p input (real), probabilities for the different categories

  ! TODO: This is not from Intel so might consider some other name

  implicit none

  ! Arguments
  integer(int32),intent(in),optional :: n
  real(real64),intent(in)            :: p(:)
  integer(int32),allocatable         :: Res(:)

  ! Other
  integer(int32) :: nOpt,i,j,k
  integer(int32) :: b(1)

  real(real64) :: pi,psum,psumtmp
  real(real64),allocatable :: pInternal(:)

  if (present(n)) then
    nOpt=n
  else
    nOpt=1
  end if

  k=size(p)
  allocate(pInternal(k))
  pInternal(:)=p(:)

  psum=sum(pInternal)
  if (abs(psum - 1.0d0) > 1e-7) then
    pInternal=pInternal/psum ! rescale
  end if
  psum=1.0d0

  allocate(Res(nOpt))

  ! Over samples
  do j=1,nOpt

    ! Over categories of a sample
    psumtmp=psum
    do i=1,k

      if (pInternal(i) > 0.0d0) then ! sample only if needed (border case)

        if (pInternal(i) < 1.0d0) then ! likewise
          pi=pInternal(i)/psumtmp

          if (pi < 1.0d0) then ! likewise
            b=SampleIntelBernoulliI(p=pi)

            if (b(1) > 0) then
              Res(j)=i
              exit
            end if

            psumtmp=psumtmp-pInternal(i)

          else
            Res(j)=i
            exit
          end if

        else
          Res(j)=i
          exit
        end if

      end if

    end do

  end do

  return

end function

!###############################################################################
