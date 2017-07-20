
!###############################################################################

function SampleIntelGammaS(n,alpha,beta,shape,scale,rate,shift) result(Res)

  ! Sample n values from a Gamma(alpha,a,beta) distribution; double precision
  ! n input (integer), number of samples to generate (default 1)
  ! alpha input (real), shape
  ! beta input (real), scale
  ! shape input (real), the same as alpha, but to make life easier
  ! scale input (real), the same as beta, but to make life easier
  ! rate input (real),  the same as 1/beta, but to make life easier
  ! a input (real), shift or displacement (default 0.0)

  ! https://software.intel.com/en-us/node/470672 (2016-03-07)

  implicit none

  ! Arguments
  integer(int32),intent(in),optional :: n
  real(real32),intent(in),optional   :: alpha
  real(real32),intent(in),optional   :: beta
  real(real32),intent(in),optional   :: shape
  real(real32),intent(in),optional   :: scale
  real(real32),intent(in),optional   :: rate
  real(real32),intent(in),optional   :: shift
  real(real32),allocatable           :: Res(:)

  ! Other
  integer(int32) :: nOpt

  real(real32) :: shapeOpt,scaleOpt,shiftOpt

  if (present(n)) then
    nOpt=n
  else
    nOpt=1
  end if

  if (.not.present(alpha) .and. .not.present(shape)) then
      write(STDERR,"(a)") "ERROR: SampleIntelGammaS requires either alpha (shape) or shape argument"
      write(STDERR,"(a)") " "
  end if

  if (.not.present(beta) .and. .not.present(scale) .and. .not.present(rate)) then
      write(STDERR,"(a)") "ERROR: SampleIntelGammaS requires either beta (scale), scale, or rate argument"
      write(STDERR,"(a)") " "
  end if

  if (present(alpha)) then
    if (.not. (alpha > 0.0)) then
      write(STDERR,"(a)") "ERROR: SampleIntelGammaS requires alpha (shape) parameter to be greater than zero"
      write(STDERR,"(a)") " "
      stop 1
    end if
    shapeOpt=alpha
  end if

  if (present(shape)) then
    if (.not. (shape > 0.0)) then
      write(STDERR,"(a)") "ERROR: SampleIntelGammaS requires shape parameter to be greater than zero"
      write(STDERR,"(a)") " "
      stop 1
    end if
    shapeOpt=shape
  end if

  if (present(beta)) then
    if (.not. (beta > 0.0)) then
      write(STDERR,"(a)") "ERROR: SampleIntelGammaS requires beta (scale) parameter to be greater than zero"
      write(STDERR,"(a)") " "
      stop 1
    end if
    scaleOpt=beta
  end if

  if (present(scale)) then
    if (.not. (scale > 0.0)) then
      write(STDERR,"(a)") "ERROR: SampleIntelGammaS requires scale parameter to be greater than zero"
      write(STDERR,"(a)") " "
      stop 1
    end if
    scaleOpt=scale
  end if

  if (present(rate)) then
    if (.not. (rate > 0.0)) then
      write(STDERR,"(a)") "ERROR: SampleIntelGammaS requires rate parameter to be greater than zero"
      write(STDERR,"(a)") " "
      stop 1
    end if
    scaleOpt=1.0/rate
    if (.not. (scaleOpt > 0.0)) then
      write(STDERR,"(a)") "ERROR: SampleIntelGammaD requires rate parameter that gives non-zero scale (=1/rate) parameter"
      write(STDERR,"(a,e)") "ERROR: rate: ",rate
      write(STDERR,"(a,e)") "ERROR: scale: ",scaleOpt
      write(STDERR,"(a)") " "
      stop 1
    end if
  end if

  if (present(shift)) then
    shiftOpt=shift
  else
    shiftOpt=0.0
  end if

  allocate(Res(nOpt))

  RNGMethod=VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE
  RNGErrCode=vsrnggamma(RNGMethod,RNGStream,nOpt,Res,shapeOpt,shiftOpt,scaleOpt)
  if (RNGErrCode /= vsl_status_ok) then
    write(STDERR,"(a)") "ERROR: SampleIntelGammaS failed"
    write(STDERR,"(a)") " "
    stop 1
  end if

  return

end function

!###############################################################################
