
!###############################################################################

function SampleIntelGaussS(n,mu,sigma2) result(Res)

  ! Sample n values from a Gauss(mu,sigma2) distribution; single precision
  ! n input (integer), number of samples to generate (default 1)
  ! mu input (real), mean (default 0.0)
  ! sigma2 input (real), variance (default 1.0)

  ! https://software.intel.com/en-us/node/470654 (2016-03-07)

  implicit none

  ! Arguments
  integer(int32),intent(in),optional :: n
  real(real32),intent(in),optional   :: mu
  real(real32),intent(in),optional   :: sigma2
  real(real32),allocatable           :: Res(:)

  ! Other
  integer(int32) :: nOpt

  real(real32) :: muOpt,sigma

  if (present(n)) then
    nOpt=n
  else
    nOpt=1
  end if

  if (present(mu)) then
    muOpt=mu
  else
    muOpt=0.0
  end if

  if (present(sigma2)) then
    sigma=sqrt(sigma2)
  else
    sigma=1.0
  end if

  allocate(Res(nOpt))

  RNGMethod=VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
  RNGErrCode=vsrnggaussian(RNGMethod,RNGStream,nOpt,Res,muOpt,sigma)
  if (RNGErrCode /= vsl_status_ok) then
    write(STDERR,"(a)") "ERROR: SampleIntelGaussS failed"
    write(STDERR,"(a)") " "
    stop 1
  end if

  return

end function

!###############################################################################
