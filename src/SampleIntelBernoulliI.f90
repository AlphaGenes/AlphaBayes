
!###############################################################################

function SampleIntelBernoulliI(n,p) result(Res)

  ! Sample n values from a Bernoulli(p) distribution
  ! n input (integer), number of samples to generate (default 1)
  ! p input (real), probability of of a success (default 0.5)

  ! https://software.intel.com/en-us/node/470686 (2014-11-25)

  implicit none

  ! Arguments
  integer(int32),intent(in),optional :: n
  real(real64),intent(in),optional   :: p
  integer(int32),allocatable         :: Res(:)

  ! Other
  integer(int32) :: nOpt

  real(real64) :: pOpt

  if (present(n)) then
    nOpt=n
  else
    nOpt=1
  end if

  if (present(p)) then
    pOpt=p
  else
    pOpt=0.5d0
  end if

  allocate(Res(nOpt))

  RNGMethod=VSL_RNG_METHOD_BERNOULLI_ICDF
  RNGErrCode=virngbernoulli(RNGMethod,RNGStream,nOpt,Res,pOpt)
  if (RNGErrCode /= vsl_status_ok) then
    write(STDERR,"(a)") "ERROR: SampleIntelBernoulliI failed"
    write(STDERR,"(a)") " "
    stop 1
  end if

  return

end function

!###############################################################################
