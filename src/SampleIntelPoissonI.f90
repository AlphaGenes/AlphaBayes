
!###############################################################################

function SampleIntelPoissonI(n,lambda) result(Res)

  ! Sample n values from a Poisson(lambda) distribution
  ! n input (integer), number of samples to generate (default 1)
  ! lambda input (real), mean and variance of distribution (default 1.0)

  ! https://software.intel.com/en-us/node/470694 (2014-12-01)

  implicit none

  ! Arguments
  integer(int32),intent(in),optional :: n
  real(real64),intent(in),optional   :: lambda
  integer(int32),allocatable         :: Res(:)

  ! Other
  integer(int32) :: nOpt

  real(real64) :: lambdaOpt

  if (present(n)) then
    nOpt=n
  else
    nOpt=1
  end if

  if (present(lambda)) then
    lambdaOpt=lambda
  else
    lambdaOpt=1.0d0
  end if

  allocate(Res(nOpt))

  if (lambdaOpt < 27.0d0) then ! based on Intel documentation https://software.intel.com/en-us/node/470598 (2014-12-01)
    RNGMethod=VSL_RNG_METHOD_POISSON_POISNORM
  else
    RNGMethod=VSL_RNG_METHOD_POISSON_PTPE
  end if
  RNGErrCode=virngpoisson(RNGMethod,RNGStream,nOpt,Res,lambdaOpt)
  if (RNGErrCode /= vsl_status_ok) then
    write(STDERR,"(a)") "ERROR: SampleIntelPoissonI failed"
    write(STDERR,"(a)") " "
    stop 1
  end if

  return

end function

!###############################################################################
