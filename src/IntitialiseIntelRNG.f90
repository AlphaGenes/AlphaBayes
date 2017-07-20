
!###############################################################################

subroutine IntitialiseIntelRNG(Seed,SeedFile,Out)
  use iso_fortran_env
  ! Start an RNG stream

  ! https://software.intel.com/en-us/node/470610 (2014-11-25)
  ! https://software.intel.com/en-us/node/470612 (2014-11-25)

  implicit none

  ! Arguments
  integer(int32),intent(in),optional  :: Seed     ! A number to initialize RNG
  character(len=*),optional           :: SeedFile ! File to save the seed in
  integer(int32),intent(out),optional :: Out      ! Make the seed value available outside

  ! Other
  integer(int32) :: Size,Unit,BRNG
  integer(int32),allocatable :: SeedList(:)

  ! Get the size of seed array
  call random_seed(size=Size)
  allocate(SeedList(Size))

  ! Set seed
  if (present(Seed)) then ! using the given value
    SeedList(1)=Seed
  else                    ! using system/compiler value
    call random_seed
    call random_seed(get=SeedList)
  end if

  ! Save to a file
  if (present(SeedFile)) then
    open(newunit=Unit,file=trim(SeedFile),status="unknown")
    write(Unit,*) SeedList(1)
    close(Unit)
  end if

  ! Output
  if (present(Out)) then
      Out=SeedList(1)
  end if

  ! Start a RNG stream
  BRNG=VSL_BRNG_MT19937
  RNGErrCode=vslnewstream(RNGStream,brng,SeedList(1))
  if (RNGErrCode /= VSL_STATUS_OK) then
    write(STDERR,"(a)") "ERROR: IntitialiseIntelRNG failed"
    write(STDERR,"(a)") " "
    stop 1
  end if

end subroutine

!###############################################################################
