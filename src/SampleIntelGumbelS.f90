
!###############################################################################

function SampleIntelGumbelS(n,a,b) result(Res)

  ! Sample n values from a Gumbel(a,b) (=type-I extreme value) distribution; single precision
  ! n input (integer), number of samples to generate (default 1)
  ! a input (real), location (default 0.0)
  ! b input (real), scale (default 1.0)

  ! https://software.intel.com/en-us/node/470670 (2014-11-25)

  implicit none

  ! Arguments
  integer(int32),intent(in),optional :: n
  real(real32),intent(in),optional   :: a
  real(real32),intent(in),optional   :: b
  real(real32),allocatable           :: Res(:)

  ! Other
  integer(int32) :: nOpt

  real(real32) :: aOpt,bOpt

  if (present(n)) then
    nOpt=n
  else
    nOpt=1
  end if

  if (present(a)) then
    aOpt=a
  else
    aOpt=0.0
  end if

  if (present(b)) then
    bOpt=b
  else
    bOpt=1.0
  end if

  allocate(Res(nOpt))

  RNGMethod=VSL_RNG_METHOD_GUMBEL_ICDF
  RNGErrCode=vsrnggumbel(RNGMethod,RNGStream,nOpt,Res,aOpt,bOpt)
  if (RNGErrCode /= vsl_status_ok) then
    write(STDERR,"(a)") "ERROR: SampleIntelGumbelS failed"
    write(STDERR,"(a)") " "
    stop 1
  end if

  return

end function

!###############################################################################
