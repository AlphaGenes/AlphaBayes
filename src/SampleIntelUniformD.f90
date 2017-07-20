
!###############################################################################

function SampleIntelUniformD(n,a,b) result(Res)

  ! Sample n values from a Uniform(a,b) distribution; double precision
  ! n input (integer), number of samples to generate (default 1)
  ! a input (real), minimal value (default 0.0)
  ! b input (real), maximal value (default 1.0)

  ! https://software.intel.com/en-us/node/470652 (2014-11-25)

  implicit none

  ! Arguments
  integer(int32),intent(in),optional :: n
  real(real64),intent(in),optional   :: a
  real(real64),intent(in),optional   :: b
  real(real64),allocatable           :: Res(:)

  ! Other
  integer(int32) :: nOpt

  real(real64) :: aOpt,bOpt

  if (present(n)) then
    nOpt=n
  else
    nOpt=1
  end if

  if (present(a)) then
    aOpt=a
  else
    aOpt=0.d0
  end if

  if (present(b)) then
    bOpt=b
  else
    bOpt=1.d0
  end if

  allocate(Res(nOpt))

  RNGMethod=VSL_RNG_METHOD_UNIFORM_STD_ACCURATE
  RNGErrCode=vdrnguniform(RNGMethod,RNGStream,nOpt,Res,aOpt,bOpt)
  if (RNGErrCode /= vsl_status_ok) then
    write(STDERR,"(a)") "ERROR: SampleIntelUniformD failed"
    write(STDERR,"(a)") " "
    stop 1
  end if

  return

end function

!###############################################################################
