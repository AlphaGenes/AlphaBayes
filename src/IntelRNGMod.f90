
!###############################################################################

include "mkl_vsl.f90"

module IntelRNGMod

  ! Random Number Generation (RNG) module
  !   "interface" to the Intel MKL Vector Statistical Library (VSL) RNG capabilities

  ! https://software.intel.com/en-us/node/470592 (2014-11-25)

  use mkl_vsl_type
  use mkl_vsl
  use ISO_Fortran_Env, STDIN=>input_unit,STDOUT=>output_unit,STDERR=>error_unit

  implicit none

  integer(int32) RNGErrCode,RNGMethod

  type(vsl_stream_state) :: RNGStream

  private

  ! RNG stream management
  public :: IntitialiseIntelRNG,UnintitialiseIntelRNG

  ! Discrete
  public :: SampleIntelUniformI
  public :: SampleIntelBernoulliI
  public :: SampleIntelMultinomialI
  public :: SampleIntelPoissonI

  ! Continuous
  ! TODO: should we make an interface and have generic for either single or double precision
  !       but do we determine single or double based on inputs or???
  public :: SampleIntelUniformS,SampleIntelUniformD
  public :: SampleIntelGaussS,SampleIntelGaussD
  public :: SampleIntelGammaS,SampleIntelGammaD
  public :: SampleIntelGumbelS,SampleIntelGumbelD

  contains

    ! RNG stream management
    include "IntitialiseIntelRNG.f90"
    include "UnintitialiseIntelRNG.f90"

    ! Discrete
    include "SampleIntelUniformI.f90"
    include "SampleIntelBernoulliI.f90"
    include "SampleIntelMultinomialI.f90"
    include "SampleIntelPoissonI.f90"

    ! Continuous
    include "SampleIntelUniformS.f90"
    include "SampleIntelUniformD.f90"
    include "SampleIntelGaussS.f90"
    include "SampleIntelGaussD.f90"
    include "SampleIntelGammaS.f90"
    include "SampleIntelGammaD.f90"
    include "SampleIntelGumbelS.f90"
    include "SampleIntelGumbelD.f90"

end module

!###############################################################################
