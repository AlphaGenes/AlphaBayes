#ifdef _WIN32

#define STRINGIFY(x)#x
#define TOSTRING(x) STRINGIFY(x)

#define DASH "\"
#define COPY "copy"
#define MD "md"
#define RMDIR "RMDIR /S /Q"
#define RM "del"
#define RENAME "MOVE /Y"
#define SH "BAT"
#define EXE ".exe"
#define NULL " >NUL"

#else

#define STRINGIFY(x)#x
#define TOSTRING(x) STRINGIFY(x)

#define DASH "/"
#define COPY "cp"
#define MD "mkdir"
#define RMDIR "rm -r"
#define RM "rm"
#define RENAME "mv"
#define SH "sh"
#define EXE ""
#define NULL ""

#endif

!###############################################################################

!-------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-------------------------------------------------------------------------------
!
!> @file     AlphaBayesModule.f90
!
! DESCRIPTION:
!> @brief    Genome-wide marker regression
!
!> @details  Genome-wide marker regression
!
!> @author   Gregor Gorjanc, gregor.gorjanc@roslin.ed.ac.uk
!
!> @date     2018-04-05
!
!> @version  0.1.0 (alpha)
!
!-------------------------------------------------------------------------------
module AlphaBayesModule

  use ISO_Fortran_Env, STDIN=>input_unit,STDOUT=>output_unit,STDERR=>error_unit
  use IntelRNGMod
  use ConstantModule, only : FILELENGTH, SPECOPTIONLENGTH, IDLENGTH
  use AlphaHouseMod
  use AlphaStatMod
  use Blas95, only : dot, & ! https://software.intel.com/en-us/mkl-developer-reference-fortran-dot  dot(x,y) does x'y
                     gemv   ! https://software.intel.com/en-us/mkl-developer-reference-fortran-gemv gemv(A,x,y) does y=Ax
  use F95_precision

  implicit none

  private
  ! Functions
  public :: AlphaBayesTitle,ReadParam,ReadData,Analysis,Prediction

  ! Module global variables :(((
  integer(int32) :: nRecTr,nMar,nCovar,nFixed,nRandom,nEffMax
  integer(int32) :: nIter,nBurn,nProcessor,GenoScaleMethod,nGenoPart,nTePop
  integer(int32),allocatable :: nRecTe(:),nMarPerGenoPart(:)
  integer(int32),allocatable :: MarPerGenoPart(:,:)

  real(real64) :: PhenMean,PhenVar,PhenSd,VarBv,VarRandom,VarMar,VarErr,MuEst,nMarR,nRecTrR
  real(real64) :: ExpVarX,SumExpVarX,ObsVarX,SumObsVarX,EpsTolerance
  real(real64),allocatable :: AlleleFreq(:),MarSd(:),PhenoTr(:),Err(:)
  real(real64),allocatable :: CovarEst(:),FixedEst(:),RandomEst(:),MarEst(:),BvEst(:)
  real(real64),allocatable :: CovarMean(:),CovarVar(:),CovarSd(:)
  real(real64),allocatable :: XpXCovar(:),XpXFixed(:),XpXRandom(:),XpXMar(:)
  real(real64),allocatable :: FixedTr(:,:),RandomTr(:,:),CovarTr(:,:),GenosTr(:,:)

  character(len=FILELENGTH) :: GenoTrFile,PhenoTrFile,CovarTrFile,FixedTrFile,RandomTrFile,EstimationMethod
  character(len=FILELENGTH),allocatable :: GenoTeFile(:),PhenoTeFile(:),GenoPartFile(:)
  character(len=IDLENGTH),allocatable :: IdTr(:)

  logical :: EstimateVariances

  ! Module parameters
  REAL(REAL64),PARAMETER :: ONER=1.0d0,ZEROR=0.0d0

  CHARACTER(LEN=100),PARAMETER :: SPECFILE="AlphaBayesSpec.txt"

  CHARACTER(LEN=100),PARAMETER :: INTERCEPTESTIMATEFILE="InterceptEstimate.txt"
  CHARACTER(LEN=100),PARAMETER :: INTERCEPTSAMPLESFILE="InterceptSamples.txt"

  CHARACTER(LEN=100),PARAMETER :: COVARIATEESTIMATEFILE="CovariateEstimate.txt"
  CHARACTER(LEN=100),PARAMETER :: COVARIATESAMPLESFILE="CovariateSamples.txt"

  CHARACTER(LEN=100),PARAMETER :: FIXEDEFFECTESTIMATEFILE="FixedEffectEstimate.txt"
  CHARACTER(LEN=100),PARAMETER :: FIXEDEFFECTSAMPLESFILE="FixedEffectSamples.txt"

  CHARACTER(LEN=100),PARAMETER :: RANDOMEFFECTESTIMATEFILE="RandomEffectEstimate.txt"
  CHARACTER(LEN=100),PARAMETER :: RANDOMEFFECTSAMPLESFILE="RandomEffectSamples.txt"

  CHARACTER(LEN=100),PARAMETER :: RANDOMEFFECTVARIANCEESTIMATEFILE="RandomEffectVarianceEstimate.txt"
  CHARACTER(LEN=100),PARAMETER :: RANDOMEFFECTVARIANCESAMPLESFILE="RandomEffectVarianceSamples.txt"

  CHARACTER(LEN=100),PARAMETER :: MARKERESTIMATEFILE="MarkerEstimate.txt"
  CHARACTER(LEN=100),PARAMETER :: MARKERSAMPLESFILE="MarkerSamples.txt"

  CHARACTER(LEN=100),PARAMETER :: MARKERVARIANCEESTIMATEFILE="MarkerVarianceEstimate.txt"
  CHARACTER(LEN=100),PARAMETER :: MARKERVARIANCESAMPLESFILE="MarkerVarianceSamples.txt"

  CHARACTER(LEN=100),PARAMETER :: RESIDUALVARIANCEESTIMATEFILE="ResidualVarianceEstimate.txt"
  CHARACTER(LEN=100),PARAMETER :: RESIDUALVARIANCESAMPLESFILE="ResidualVarianceSamples.txt"

  CHARACTER(LEN=100),PARAMETER :: GENOPARTFILESTART="GenomePartition"
  CHARACTER(LEN=100),PARAMETER :: GENOPARTESTIMATEFILEEND="Estimate.txt"
  CHARACTER(LEN=100),PARAMETER :: GENOPARTSAMPLESFILEEND="Samples.txt"

  contains

    !###########################################################################

    subroutine AlphaBayesTitle
      implicit none
      write(STDOUT, "(a)") ""
      write(STDOUT, "(a)") "                            ***********************                           "
      write(STDOUT, "(a)") "                            *                     *                           "
      write(STDOUT, "(a)") "                            *      AlphaBayes     *                           "
      write(STDOUT, "(a)") "                            *                     *                           "
      write(STDOUT, "(a)") "                            ***********************                           "
      write(STDOUT, "(a)") "                                                                              "
      write(STDOUT, "(a)") "          Software for estimating marker effects and genomic prediction       "
      write(STDOUT, "(a)") "                       http://AlphaGenes.Roslin.ed.ac.uk                      "
      write(STDOUT, "(a)") "                                 No liability                                 "
      write(STDOUT, "(a)") ""
      write(STDOUT, "(a)") "                       Commit:   "//TOSTRING(COMMIT)//"                       "
      write(STDOUT, "(a)") "                       Compiled: "//__DATE__//", "//__TIME__
      write(STDOUT, "(a)") ""
    end subroutine

    !###########################################################################

    subroutine ReadParam
      implicit none

      integer(int32) :: SpecUnit, Stat, GTePop, PTePop, nRecTeI, GenoPart

      real(real64) :: TmpR

      logical :: LogStdoutInternal, PhenoTrFileGiven, CovarTrFileGiven, FixedTrFileGiven
      logical :: RandomTrFileGiven, GenoTrFileGiven

      character(len=:), allocatable :: DumString
      character(len=SPECOPTIONLENGTH) :: Line
      character(len=SPECOPTIONLENGTH) :: First
      character(len=SPECOPTIONLENGTH), allocatable, dimension(:) :: Second

      open(newunit=SpecUnit, file=SPECFILE, status="old")

      ! Defaults
      LogStdoutInternal = .true.
      PhenoTrFile = ""
      nRecTr = 0
      PhenoTrFileGiven = .false.
      CovarTrFile = ""
      nCovar = 0
      CovarTrFileGiven = .false.
      FixedTrFile = ""
      nFixed = 0
      FixedTrFileGiven = .false.
      ! RandomTrFile allocatable
      nRandom = 0
      RandomTrFileGiven = .false.
      GenoTrFile = ""
      nMar = 0
      GenoTrFileGiven = .false.
      ! nRecTe allocatable
      nTePop = 0
      ! PhenoTeFile allocatable
      ! GenoTeFile allocatable
      GTePop = 0
      PTePop = 0
      nRecTeI = 0
      nIter = 10000
      nBurn = 1000
      VarRandom = -1.0d0
      VarBv = -1.0d0
      VarErr = -1.0d0
      nProcessor = 1
      GenoScaleMethod = 4
      EstimationMethod = "RidgeSolve"
      EstimateVariances = .false.
      nGenoPart = 0
      GenoPart = 0
      EpsTolerance=1.0E-8

      ! Process spec file
      Stat = 0
      ReadSpec: do while (Stat .eq. 0)
        read(SpecUnit, "(a)", iostat=Stat) Line
        if (len_trim(Line) .eq. 0) then
          cycle
        end if
        call SplitLineIntoTwoParts(trim(adjustl(Line)), First, Second)
        DumString = ParseToFirstWhitespace(First)
        ! @todo why (len_trim(Line) .eq. 0)? if we use (len_trim(Line) .eq. 0) above
        if (First(1:1) .eq. "=" .or. len_trim(Line) .eq. 0) then
          cycle
        else
          select case (ToLower(trim(DumString)))
            ! Inputs
            case ("phenotypetrainfile")
              if (allocated(Second)) then
                PhenoTrFileGiven = .true.
                write(PhenoTrFile, *) trim(adjustl(Second(1)))
                PhenoTrFile = adjustl(PhenoTrFile)
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Phenotype train file: "//trim(PhenoTrFile)
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for PhenotypeTrainFile, i.e., PhenotypeTrainFile, PhenotypesTrain.txt"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberoftrainrecords")
              if (allocated(Second)) then
                nRecTr = Char2Int(trim(adjustl(Second(1))))
                nRecTrR = dble(nRecTr)
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Number of train records: "//trim(Int2Char(nRecTr))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfTrainRecords, i.e., NumberOfTrainRecords, 100"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("covariatetrainfile")
              if (allocated(Second)) then
                CovarTrFileGiven = .true.
                write(CovarTrFile, *) trim(adjustl(Second(1)))
                CovarTrFile = adjustl(CovarTrFile)
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Covariate train file: "//trim(CovarTrFile)
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for CovariateTrainFile, i.e., CovariateTrainFile, CovariatesTrain.txt"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberofcovariates")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .ne. "0") then
                  nCovar = Char2Int(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Number of covariates: "//trim(Int2Char(nTePop))
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfCovariates, i.e., NumberOfCovariates, 2"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("fixedeffecttrainfile")
              if (allocated(Second)) then
                FixedTrFileGiven = .true.
                write(FixedTrFile, *) trim(adjustl(Second(1)))
                FixedTrFile = adjustl(FixedTrFile)
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Fixed effect train file: "//trim(FixedTrFile)
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for FixedEffectTrainFile, i.e., FixedEffectTrainFile, FixedEffectsTrain.txt"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberoffixedeffectlevels")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .ne. "0") then
                  nFixed = Char2Int(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Number of fixed effect levels: "//trim(Int2Char(nFixed))
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfFixedEffectLevels, i.e., NumberOfFixedEffectLevels, 2"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("randomeffecttrainfile")
              if (allocated(Second)) then
                RandomTrFileGiven = .true.
                write(RandomTrFile, *) trim(adjustl(Second(1)))
                RandomTrFile = adjustl(RandomTrFile)
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Random effect train file: "//trim(RandomTrFile)
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for RandomTrFileEffectTrainFile, i.e., RandomEffectTrainFile, RandomEffectsTrain.txt"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberofrandomeffectlevels")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .ne. "0") then
                  nRandom = Char2Int(trim(adjustl(Second(1))))
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Number of random effect levels: "//trim(Int2Char(nRandom))
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfRandomEffectLevels, i.e., NumberOfRandomEffectLevels, 20"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("genotypetrainfile")
              if (allocated(Second)) then
                GenoTrFileGiven = .true.
                write(GenoTrFile, *) trim(adjustl(Second(1)))
                GenoTrFile = adjustl(GenoTrFile)
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Genotype train file: "//trim(GenoTrFile)
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for GenotypeTrainFile, i.e., GenotypeTrainFile, GenotypesTrain.txt"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberofmarkers")
              if (allocated(Second)) then
                nMar = Char2Int(trim(adjustl(Second(1))))
                nMarR = dble(nMar)
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Number of markers: "//trim(Int2Char(nMar))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfMarkers, i.e., NumberOfMarkers, 100"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberofpredictsets")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .ne. "0") then
                  nTePop = Char2Int(trim(adjustl(Second(1))))
                  allocate(PhenoTeFile(nTePop))
                  allocate(GenoTeFile(nTePop))
                  allocate(nRecTe(nTePop))
                  nRecTe = 0
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Number of prediction sets: "//trim(Int2Char(nTePop))
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfPredictSets, i.e., NumberOfPredictSets, 2"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("phenotypepredictfile")
              if (allocated(Second)) then
                PTePop = PTePop + 1
                if (PTePop .gt. nTePop) then
                  write(STDERR, "(a)") " ERROR: Too many PhenotypePredictFile specifications vs. NumberOfPredictSets"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
                write(PhenoTeFile(PTePop), *) trim(adjustl(Second(1)))
                PhenoTeFile(PTePop) = adjustl(PhenoTeFile(PTePop))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Phenotype predict file: "//trim(PhenoTeFile(PTePop))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for PhenotypePredictFile, i.e., PhenotypePredictFile, PhenotypesPredict.txt"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberofpredictrecords")
              if (allocated(Second)) then
                nRecTeI = nRecTeI + 1
                if (nRecTeI .gt. nTePop) then
                  write(STDERR, "(a)") " ERROR: Too many NumberOfPredictRecords specifications vs. NumberOfPredictSets"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
                nRecTe(nRecTeI) = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Number of prediction records: "//trim(Int2Char(nRecTe(nRecTeI)))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfPredictRecords, i.e., NumberOfPredictRecords, 100"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("genotypepredictfile")
              if (allocated(Second)) then
                GTePop = GTePop + 1
                if (GTePop .gt. nTePop) then
                  write(STDERR, "(a)") " ERROR: Too many GenotypePredictFile specifications vs. NumberOfPredictSets"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
                write(GenoTeFile(GTePop), *) trim(adjustl(Second(1)))
                GenoTeFile(GTePop) = adjustl(GenoTeFile(GTePop))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Genotype predict file: "//trim(GenoTeFile(GTePop))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for GenotypePredictFile, i.e., GenotypePredictFile, GenotypesPredict.txt"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberofiterations")
              if (allocated(Second)) then
                nIter = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Number of iterations: "//trim(Int2Char(nIter))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfIterations, i.e., NumberOfIterations, 10000"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberofburniniterations")
              if (allocated(Second)) then
                nBurn = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Number of burn-in iterations: "//trim(Int2Char(nBurn))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfBurnInIterations, i.e., NumberOfBurnInIterations, 1000"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("geneticvariance")
              if (allocated(Second)) then
                VarBv = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Genetic variance: "//trim(Real2Char(VarBv))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for GeneticVariance, i.e., GeneticVariance, 0.5"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("residualvariance")
              if (allocated(Second)) then
                VarErr = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Residual variance: "//trim(Real2Char(VarErr))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for ResidualVariance, i.e., ResidualVariance, 0.5"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberofprocessors")
              if (allocated(Second)) then
                nProcessor = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Number of burn-in iterations: "//trim(Int2Char(nProcessor))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfProcessors, i.e., NumberOfProcessors, 4"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("genotypescalingmethod")
              if (allocated(Second)) then
                GenoScaleMethod = Char2Int(trim(adjustl(Second(1))))
                if (GenoScaleMethod .lt. 1 .or. GenoScaleMethod .gt. 5) then
                  write(STDERR,"(a)") " ERROR: GenotypeScalingMethod must be between 1 and 5"
                  write(STDERR,"(a)") " "
                  stop 1
                endif
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Genotype scaling method: "//trim(Int2Char(GenoScaleMethod))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for GenotypeScalingMethod, i.e., GenotypeScalingMethod, 4"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("estimationmethod")
              if (allocated(Second)) then
                write(EstimationMethod, *) trim(adjustl(Second(1)))
                EstimationMethod = adjustl(EstimationMethod)
                if (.not.(trim(EstimationMethod) .eq. "RidgeSolve" .or. trim(EstimationMethod) .eq. "RidgeSample")) then
                  write(STDERR,"(a)") " ERROR: EstimationMethod must be either RidgeSolve or RidgeSample"
                  write(STDERR,"(a)") " "
                  stop 1
                endif
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Estimation method: "//trim(EstimationMethod)
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a string for EstimationMethod, i.e., EstimationMethod, RidgeSolve"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("estimatevariances")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .eq. "yes") then
                  EstimateVariances = .true.
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Estimate variances: Yes"
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a string for EstimateVariances, i.e., EstimateVariances, Yes"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("convergencetolerance")
              if (allocated(Second)) then
                EpsTolerance = Char2Double(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Convergence tolerance: "//trim(Real2Char(EpsTolerance))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a value for ConvergenceTolerance, i.e., ConvergenceTolerance, 1.0E-8"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberofgenomepartitions")
              if (allocated(Second)) then
                nGenoPart = Char2Int(trim(adjustl(Second(1))))
                allocate(GenoPartFile(nGenoPart))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Number of genome partitions: "//trim(Int2Char(nGenoPart))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfGenomePartitions, i.e., NumberOfGenomePartitions, 4"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("genomepartitionfile")
              if (allocated(Second)) then
                GenoPart = GenoPart + 1
                if (GenoPart .gt. nGenoPart) then
                  write(STDERR, "(a)") " ERROR: Too many GenomePartitionFile specifications vs. NumberOfGenomePartitions"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
                write(GenoPartFile(GenoPart), *) trim(adjustl(Second(1)))
                GenoPartFile(GenoPart) = adjustl(GenoPartFile(GenoPart))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Genome partition file: "//trim(GenoPartFile(GenoPart))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a file for GenomePartitionFile, i.e., GenomePartitionFile, Partition1.txt"
                write(STDERR, "(a)") " "
                stop 1
              end if

              case ("stop")
              if (LogStdoutInternal) then
                write(STDOUT, "(3a)") " NOTE: Encountered Stop specification - the rest of specifications will be ignored"
              end if
              exit

            case default
              if (LogStdoutInternal) then
                write(STDOUT, "(a)") " NOTE: Specification '"//trim(Line)//"' was ignored"
                write(STDOUT, "(a)") " "
              end if
          end select
        end if
      end do ReadSpec
      close(SpecUnit)

      if (.not. PhenoTrFileGiven) then
        write(STDERR, "(a)") " ERROR: Must specify a phenotype train file, i.e., PhenotypeTrainFile, PhenotypesTrain.txt"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (CovarTrFileGiven .and. nCovar .eq. 0) then
        write(STDERR, "(a)") " ERROR: Must specify number of covariates in the covariate file, i.e., NumberOfCovariates, 2"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (FixedTrFileGiven .and. nFixed .eq. 0) then
        write(STDERR, "(a)") " ERROR: Must specify number of levels in the fixed effects file, i.e., NumberOfFixedEffectLevels, 2"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (RandomTrFileGiven .and. nRandom .eq. 0) then
        write(STDERR, "(a)") " ERROR: Must specify number of levels in the random effects file, i.e., NumberOfRandomEffectLevels, 20"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (GenoTrFileGiven .and. nMar .eq. 0) then
        write(STDERR, "(a)") " ERROR: Must specify number of markers in the genotype file, i.e., NumberOfMarkers, 100"
        write(STDERR, "(a)") " "
        stop 1
      end if

      TmpR = 1.0d0
      if (RandomTrFileGiven .and. VarRandom .lt. 0) then
        TmpR = TmpR + 1.0d0
      end if
      if (GenoTrFileGiven .and. VarBv .lt. 0) then
        TmpR = TmpR + 1.0d0
      end if

      if (RandomTrFileGiven .and. VarRandom .lt. 0) then
        VarRandom = 1.0d0 / TmpR
        write(STDOUT, "(a)") " NOTE: Variance of the random effect set to: "//trim(Real2Char(VarRandom))
        write(STDERR, "(a)") " "
      end if
      if (GenoTrFileGiven .and. VarBv .lt. 0) then
        VarBv = 1.0d0 / TmpR
        write(STDOUT, "(a)") " NOTE: Variance of the genetic effect set to: "//trim(Real2Char(VarBv))
        write(STDERR, "(a)") " "
      end if
      if (VarErr .lt. 0) then
        VarErr = 1.0d0 / TmpR
        write(STDOUT, "(a)") " NOTE: Variance of the residual set to: "//trim(Real2Char(VarErr))
        write(STDERR, "(a)") " "
      end if

      nEffMax = maxval([nCovar,nFixed,nRandom,nMar])
      call OMP_SET_NUM_THREADS(nProcessor)
    end subroutine

    !###########################################################################

    subroutine ReadData
      implicit none

      integer(int32) :: i,j,nNotMissing,PhenoTrUnit,CovarTrUnit,FixedTrUnit
      integer(int32) :: RandomTrUnit,GenoTrUnit,Unit!,AlleleFreqUnit

      character(len=100) :: DumC

      real(real64) :: TmpStat

      if (nRecTr.eq.0) then
        nRecTr = CountLines(PhenoTrFile)
        nRecTrR = dble(nRecTr)
      end if

      open(newunit=PhenoTrUnit,file=trim(PhenoTrFile),status="old")

      allocate(PhenoTr(nRecTr))
      allocate(IdTr(nRecTr))
      if (nCovar.gt.0) then
        allocate(CovarTr(nRecTr,nCovar))
        allocate(CovarEst(nCovar))
        allocate(CovarMean(nCovar))
        allocate(CovarVar(nCovar))
        allocate(CovarSd(nCovar))
        open(newunit=CovarTrUnit,file=trim(CovarTrFile),status="old")
      end if
      if (nFixed.gt.0) then
        open(newunit=FixedTrUnit,file=trim(FixedTrFile),status="old")
        allocate(FixedTr(nRecTr,nFixed))
        allocate(FixedEst(nFixed))
      end if
      if (nRandom.gt.0) then
        open(newunit=RandomTrUnit,file=trim(RandomTrFile),status="old")
        allocate(RandomTr(nRecTr,nRandom))
        allocate(RandomEst(nRandom))
      end if
      if (nMar.gt.0) then
        allocate(GenosTr(nRecTr,nMar))
        allocate(MarEst(nMar))
        allocate(AlleleFreq(nMar))
        allocate(MarSd(nMar))
        open(newunit=GenoTrUnit,file=trim(GenoTrFile),status="old")
        !open(newunit=AlleleFreqUnit,file="AlleleFreq.txt",status="unknown")
        !write(AlleleFreqUnit,"(a11,1x,a11)") "Mar","AlleleFreq"
      end if

      do i=1,nRecTr
        read(PhenoTrUnit,*) IdTr(i),PhenoTr(i)
        if (nCovar.gt.0) then
          read(CovarTrUnit,*) DumC,CovarTr(i,:)
          if (trim(IdTr(i)).ne.trim(DumC)) then
            write(STDERR,"(a)") " ERROR: Individual identifications in the phenotype and covariate files do not match"
            write(STDERR,"(a,i)") " ERROR: Line: ",i
            write(STDERR,"(a,a)") " ERROR: Phenotype file identification: ",trim(IdTr(i))
            write(STDERR,"(a,a)") " ERROR: Covariate file identification: ",trim(DumC)
            write(STDERR,"(a)") " "
            stop 1
          end if
        end if
        if (nFixed.gt.0) then
          read(FixedTrUnit,*) DumC,FixedTr(i,:)
          if (trim(IdTr(i)).ne.trim(DumC)) then
            write(STDERR,"(a)") " ERROR: Individual identifications in the phenotype and fixed effect files do not match"
            write(STDERR,"(a,i)") " ERROR: Line: ",i
            write(STDERR,"(a,a)") " ERROR: Phenotype file identification: ",trim(IdTr(i))
            write(STDERR,"(a,a)") " ERROR: Fixed effect file identification: ",trim(DumC)
            write(STDERR,"(a)") " "
            stop 1
          end if
        end if
        if (nRandom.gt.0) then
          read(RandomTrUnit,*) DumC,RandomTr(i,:)
          if (trim(IdTr(i)).ne.trim(DumC)) then
            write(STDERR,"(a)") " ERROR: Individual identifications in the phenotype and random effect files do not match"
            write(STDERR,"(a,i)") " ERROR: Line: ",i
            write(STDERR,"(a,a)") " ERROR: Phenotype file identification: ",trim(IdTr(i))
            write(STDERR,"(a,a)") " ERROR: Random effect file identification: ",trim(DumC)
            write(STDERR,"(a)") " "
            stop 1
          end if
        end if
        if (nMar.gt.0) then
          read(GenoTrUnit,*) DumC,GenosTr(i,:)
          if (trim(IdTr(i)).ne.trim(DumC)) then
            write(STDERR,"(a)") " ERROR: Individual identifications in the phenotype and genotype files do not match in the training set"
            write(STDERR,"(a,i)") " ERROR: Line: ",i
            write(STDERR,"(a,a)") " ERROR: Phenotype file identification: ",trim(IdTr(i))
            write(STDERR,"(a,a)") " ERROR: Genotype file identification: ",trim(DumC)
            write(STDERR,"(a)") " "
            stop 1
          end if
        end if
      end do
      close(GenoTrUnit)
      if (nCovar.gt.0) then
        close(CovarTrUnit)
      end if
      if (nFixed.gt.0) then
        close(FixedTrUnit)
      end if
      if (nRandom.gt.0) then
        close(RandomTrUnit)
      end if
      if (nMar.gt.0) then
        close(PhenoTrUnit)
        ! flush(AlleleFreqUnit)
        ! close(AlleleFreqUnit)
      end if

      PhenMean=Mean(PhenoTr)
      PhenVar=Var(PhenoTr,PhenMean)
      PhenSd=sqrt(PhenVar)
      PhenoTr=(PhenoTr-PhenMean)/PhenSd
      if (nCovar.gt.0) then
        do i=1,nCovar
          CovarMean(i)=Mean(CovarTr(:,i))
          CovarVar(i)=Var(CovarTr(:,i),CovarMean(i))
          CovarSd(i)=sqrt(CovarVar(i))
          CovarTr(:,i)=(CovarTr(:,i)-CovarMean(i))/CovarSd(i)
        end do
      end if
      if (nMar.gt.0) then
        SumExpVarX=0.0d0
        SumObsVarX=0.0d0
        do j=1,nMar

          ! Compute allele freqs
          AlleleFreq(j)=0.0d0
          nNotMissing=0
          do i=1,nRecTr
            if ((GenosTr(i,j).ge.0.0).and.(GenosTr(i,j).le.2.0)) then
              AlleleFreq(j)=AlleleFreq(j)+GenosTr(i,j)
              nNotMissing=nNotMissing+1
            end if
          end do
          if (nNotMissing.ne.0) then
            AlleleFreq(j)=AlleleFreq(j)/dble(2*nNotMissing)
          else
            AlleleFreq(j)=0.0d0
          end if
          !write(AlleleFreqUnit,"(i11,1x,f11.8)") j,AlleleFreq(j)

          ! Fix any odd data
          do i=1,nRecTr
            if ((GenosTr(i,j).lt.0.0).or.(GenosTr(i,j).gt.2.0)) then
              GenosTr(i,j)=2.0d0*AlleleFreq(j)
            end if
          end do

          ! Standardize
          TmpStat=Var(GenosTr(:,j))
          ExpVarX=2.0d0*(1.0d0-AlleleFreq(j))*AlleleFreq(j)+0.00001d0 ! if p=0.00001, then 2*p*q=0.00001
          SumExpVarX=SumExpVarX+ExpVarX
          ObsVarX=TmpStat+0.00001d0
          SumObsVarX=SumObsVarX+ObsVarX

          ! ... center
          GenosTr(:,j)=GenosTr(:,j)-(2.0d0*AlleleFreq(j))

          ! ... scale
          MarSd(j)=1.0d0

          if (GenoScaleMethod.eq.2) then
            ! Scale by marker specific variance - expected
            ExpVarX=sqrt(ExpVarX)
            MarSd(j)=ExpVarX
            GenosTr(:,j)=GenosTr(:,j)/MarSd(j)
          end if

          if (GenoScaleMethod.eq.3) then
            ! Scale by marker specific variance - observed
            ObsVarX=sqrt(ObsVarX)
            MarSd(j)=ObsVarX
            GenosTr(:,j)=GenosTr(:,j)/MarSd(j)
          end if

        end do

        if (GenoScaleMethod.eq.4) then
          ! Scale by average marker variance - expected
          ExpVarX=sqrt(SumExpVarX/nMarR)
          MarSd=ExpVarX
          GenosTr=GenosTr/MarSd(1)
        end if

        if (GenoScaleMethod.eq.5) then
          ! Scale by average marker variance - observed
          ObsVarX=sqrt(SumObsVarX/nMarR)
          MarSd=ObsVarX
          GenosTr=GenosTr/MarSd(1)
        end if
      end if

      ! A note about scaling and its influence on marker estimates - if we
      ! want to get their estimate for the observed phenotype and genotype scale
      ! y     = mu     + b*x     + e
      ! y     = mu     + b*x/SDX + e
      ! y     = mu     + b*z     + e
      ! y/SDY = mu/SDY + b*z/SDY + e/SDY
      ! y'    = mu'    + b'*z    + e'
      ! b'*z=b*z/SDY
      ! b*z =b'*z*SDY
      ! b*x/SDX=b'*x/SDX*SDY
      ! b=b'/SDX*SDY

      ! open(newunit=GenoTrUnit,file="GenoTrainProcessed.txt",status="unknown")
      ! open(newunit=PhenoTrUnit,file="PhenoTrainProcessed.txt",status="unknown")
      ! do i=1,nRecTr
      !   write(GenoTrUnit,"("//Int2Char(nMar)//"(1x,f20.16))") GenosTr(i,:)
      !   write(PhenoTrUnit,*) PhenoTr(i)
      ! end do
      ! flush(GenoTrUnit)
      ! close(GenoTrUnit)
      ! flush(PhenoTrUnit)
      ! close(PhenoTrUnit)

      if (nMar.gt.0.and.nGenoPart.gt.0) then
        allocate(nMarPerGenoPart(nGenoPart))
        do i=1,nGenoPart
          nMarPerGenoPart(i)=CountLines(GenoPartFile(i))
          ! write(STDOUT,*) i,nMarPerGenoPart(i)
        end do
        allocate(MarPerGenoPart(maxval(nMarPerGenoPart),nGenoPart))
        do i=1,nGenoPart
          open(newunit=Unit,file=trim(GenoPartFile(i)),status="unknown")
          do j=1,nMarPerGenoPart(i)
            read(Unit,*) MarPerGenoPart(j,i)
          end do
          close(Unit)
          ! write(STDOUT,*) i,nMarPerGenoPart(i),MarPerGenoPart(1:nMarPerGenoPart(i),i)
        end do
      end if

      ! Here as both RidgeRegressionSolve and RidgeRegressionSample use them
      allocate(Err(nRecTr))
      if (nCovar.gt.0) then
        allocate(XpXCovar(nCovar))
        do j=1,nCovar
          XpXCovar(j)=dot(x=CovarTr(:,j),y=CovarTr(:,j)) + tiny(CovarTr(1,1))
        end do
      end if
      if (nFixed.gt.0) then
        allocate(XpXFixed(nFixed))
        do j=1,nFixed
          XpXFixed(j)=dot(x=FixedTr(:,j),y=FixedTr(:,j)) + tiny(FixedTr(1,1))
        end do
      end if
      if (nRandom.gt.0) then
        allocate(XpXRandom(nRandom))
        do j=1,nRandom
          XpXRandom(j)=dot(x=RandomTr(:,j),y=RandomTr(:,j)) + tiny(RandomTr(1,1))
        end do
      end if
      if (nMar.gt.0) then
        allocate(XpXMar(nMar))
        do j=1,nMar
          XpXMar(j)=dot(x=GenosTr(:,j),y=GenosTr(:,j)) + tiny(GenosTr(1,1))
        end do
      end if
    end subroutine

    !###########################################################################

    subroutine Analysis
      implicit none
      if (trim(EstimationMethod).eq."RidgeSolve") then
        call RidgeRegressionSolve
      end if
      if (trim(EstimationMethod).eq."RidgeSample") then
        call RidgeRegressionSolve
        call RidgeRegressionSample
      end if
    end subroutine

    !###########################################################################

    subroutine RidgeRegressionSolve
      implicit none

      integer(int32) :: Iter,i,j,k,RandomOrdering(nEffMax),Unit

      real(real64) :: PreMar,PreRandom,Rhs,Lhs,Sol,Diff,Eps
      real(real64),allocatable :: BvEstGenoPart(:)

      write(STDOUT, "(a)") ""
      write(STDOUT, "(a)") " Running estimation of marker effects with provided variance components"
      write(STDOUT, "(a)") ""

      allocate(BvEst(nRecTr))
      if (nGenoPart.gt.0) then
        allocate(BvEstGenoPart(nRecTr))
      end if

      ! Working phenotypes (centered) --> residuals
      Err=PhenoTr

      ! Initialise
      MuEst=0.0d0

      ! Initialise and construct XpXTauE
      if (nCovar.gt.0) then
        CovarEst=0.0d0
        XpXCovar=XpXCovar/VarErr ! can do it only once for all rounds!!!
      end if
      if (nFixed.gt.0) then
        FixedEst=0.0d0
        XpXFixed=XpXFixed/VarErr ! can do it only once for all rounds!!!
      end if
      if (nRandom.gt.0) then
        RandomEst=0.0d0
        XpXRandom=XpXRandom/VarErr ! can do it only once for all rounds!!!
        PreRandom=1.0d0/VarRandom
      end if
      if (nMar.gt.0) then
        MarEst=0.0d0
        XpXMar=XpXMar/VarErr ! can do it only once for all rounds!!!
        VarMar=VarBv/nMarR ! approximate variance of allele substitution effects
        PreMar=1.0d0/VarMar
      end if

      ! Iterate
      do Iter=1,nIter
        Eps=0.0d0

        ! Intercept
        Lhs=nRecTrR/VarErr
        Rhs=sum(Err)/VarErr + nRecTrR*MuEst/VarErr
        Sol=Rhs/Lhs
        Diff=Sol-MuEst
        Err=Err-Diff
        MuEst=Sol
        Eps=Eps+Diff*Diff

        ! Covariates
        if (nCovar.gt.0) then
          RandomOrdering(1:nCovar)=RandomOrder(nCovar)
          do j=1,nCovar
            k=RandomOrdering(j)
            Lhs=XpXCovar(k)
            Rhs=dot(x=CovarTr(:,k),y=Err)/VarErr + XpXCovar(k)*CovarEst(k)
            Sol=Rhs/Lhs
            Diff=Sol-CovarEst(k)
            Err=Err-CovarTr(:,k)*Diff
            CovarEst(k)=Sol
            Eps=Eps+Diff*Diff
          end do
        end if

        ! Fixed effects
        if (nFixed.gt.0) then
          RandomOrdering(1:nFixed)=RandomOrder(nFixed)
          do j=1,nFixed
            k=RandomOrdering(j)
            Lhs=XpXFixed(k)
            Rhs=dot(x=FixedTr(:,k),y=Err)/VarErr + XpXFixed(k)*FixedEst(k)
            Sol=Rhs/Lhs
            Diff=Sol-FixedEst(k)
            Err=Err-FixedTr(:,k)*Diff
            FixedEst(k)=Sol
            Eps=Eps+Diff*Diff
          end do
        end if

        ! Random effects
        if (nRandom.gt.0) then
          RandomOrdering(1:nRandom)=RandomOrder(nRandom)
          do j=1,nRandom
            k=RandomOrdering(j)
            Lhs=XpXRandom(k) + PreRandom
            Rhs=dot(x=RandomTr(:,k),y=Err)/VarErr + XpXRandom(k)*RandomEst(k)
            Sol=Rhs/Lhs
            Diff=Sol-RandomEst(k)
            Err=Err-RandomTr(:,k)*Diff
            RandomEst(k)=Sol
            Eps=Eps+Diff*Diff
          end do
        end if

        ! Markers
        if (nMar.gt.0) then
          RandomOrdering(1:nMar)=RandomOrder(nMar)
          do j=1,nMar
            k=RandomOrdering(j)
            Lhs=XpXMar(k) + PreMar
            Rhs=dot(x=GenosTr(:,k),y=Err)/VarErr + XpXMar(k)*MarEst(k)
            Sol=Rhs/Lhs
            Diff=Sol-MarEst(k)
            Err=Err-GenosTr(:,k)*Diff
            MarEst(k)=Sol
            Eps=Eps+Diff*Diff
          end do
        end if

        ! Recompute working residuals to avoid rounding errors
        if (mod(Iter,100).eq.0) then
          Err=PhenoTr-MuEst
          if (nCovar.gt.0) then
            call gemv(A=CovarTr,x=CovarEst,y=BvEst)
            Err=Err-BvEst
          end if
          if (nFixed.gt.0) then
            call gemv(A=FixedTr,x=FixedEst,y=BvEst)
            Err=Err-BvEst
          end if
          if (nRandom.gt.0) then
            call gemv(A=RandomTr,x=RandomEst,y=BvEst)
            Err=Err-BvEst
          end if
          if (nMar.gt.0) then
            call gemv(A=GenosTr,x=MarEst,y=BvEst)
            Err=Err-BvEst
          end if
        end if

        ! Stopping criteria
        if (Eps.lt.EpsTolerance) then
          exit
        end if
      end do

      open(newunit=Unit,file=INTERCEPTESTIMATEFILE,status="unknown")
      write(Unit,*) MuEst*PhenSd + PhenMean
      flush(Unit)
      close(Unit)

      if (nCovar.gt.0) then
        open(newunit=Unit,file=COVARIATEESTIMATEFILE,status="unknown")
        do i=1,nCovar
          write(Unit,*) CovarEst(i)/CovarSd(i)*PhenSd
        end do
        flush(Unit)
        close(Unit)
      end if

      if (nFixed.gt.0) then
        open(newunit=Unit,file=FIXEDEFFECTESTIMATEFILE,status="unknown")
        do i=1,nFixed
          write(Unit,*) FixedEst(i)*PhenSd
        end do
        flush(Unit)
        close(Unit)
      end if

      if (nRandom.gt.0) then
        open(newunit=Unit,file=RANDOMEFFECTESTIMATEFILE,status="unknown")
        do i=1,nRandom
          write(Unit,*) RandomEst(i)*PhenSd
        end do
        flush(Unit)
        close(Unit)
      end if

      if (nMar.gt.0) then
        open(newunit=Unit,file=MARKERESTIMATEFILE,status="unknown")
        do i=1,nMar
          write(Unit,*) MarEst(i)/MarSd(i)*PhenSd
        end do
        flush(Unit)
        close(Unit)

        if (nGenoPart.gt.0) then
          call gemv(A=GenosTr,x=MarEst/MarSd*PhenSd,y=BvEst)
          do i=1,nGenoPart
            BvEstGenoPart=0.0d0
            ! Might have used gemv here, but for MCMC this would mean calling gemv nIter*nGenoPart times!!!
            do j=1,nMarPerGenoPart(i)
              k=MarPerGenoPart(j,i)
              BvEstGenoPart=BvEstGenoPart+(GenosTr(:,k)*MarEst(k)/MarSd(k)*PhenSd)
            end do
            BvEst=BvEst-BvEstGenoPart
            open(newunit=Unit,file=trim(GENOPARTFILESTART)//trim(Int2Char(i))//trim(GENOPARTESTIMATEFILEEND),status="unknown")
            do j=1,nRecTr
              write(Unit,*) BvEstGenoPart(j)
            end do
            flush(Unit)
            close(Unit)
          end do
          open(newunit=Unit,file=trim(GENOPARTFILESTART)//trim("Remainder")//trim(GENOPARTESTIMATEFILEEND),status="unknown")
          do j=1,nRecTr
            write(Unit,*) BvEst(j)
          end do
          flush(Unit)
          close(Unit)
        end if
        ! Convert XpXTauE=XpX/VarErr back to XpX
        if (nCovar.gt.0) then
          XpXCovar=XpXCovar*VarErr
        end if
        if (nFixed.gt.0) then
          XpXFixed=XpXFixed*VarErr
        end if
        if (nRandom.gt.0) then
          XpXRandom=XpXRandom*VarErr
        end if
        if (nMar.gt.0) then
          XpXMar=XpXMar*VarErr
        end if
        if (nGenoPart.gt.0) then
          deallocate(BvEstGenoPart)
        end if
      end if
    end subroutine

    !###########################################################################

    subroutine RidgeRegressionSample
      implicit none

      integer(int32) :: Iter,i,j,k,RandomOrdering(nEffMax)
      integer(int32) :: Unit,MuUnit,CovarUnit,FixedUnit,RandomUnit,MarUnit
      integer(int32) :: VarMarUnit,VarRandomUnit,VarErrUnit
      integer(int32),allocatable :: GenoPartUnit(:)

      real(real64) :: TmpR,Rhs,Lhs,Sol,Diff,nSampR
      real(real64) :: MuSamp,VarRandomSamp,VarRandomEst,VarMarSamp,VarMarEst,VarErrSamp,VarErrEst
      real(real64) :: R2,RandomDF0,RandomDF,RandomS0,MarDF0,MarDF,MarS0,ErrDF0,ErrDF,ErrS0,MSX
      real(real64),allocatable :: CovarSamp(:),FixedSamp(:),RandomSamp(:),MarSamp(:)
      real(real64),allocatable :: SX2(:),MX2(:),BvSamp(:),BvSampGenoPart(:)
      real(real64),allocatable :: GaussDevMu(:),GaussDevCovar(:),GaussDevFixed(:)
      real(real64),allocatable :: GaussDevRandom(:),GaussDevMar(:)
      real(real64),allocatable :: GammaDevRandom(:),GammaDevMar(:),GammaDevErr(:)

      character(len=100) :: CovarSampFmt,FixedSampFmt,RandomSampFmt,MarSampFmt,BvSampFmt

      if (EstimateVariances) then
        write(STDOUT, "(a)") ""
        write(STDOUT, "(a)") " Running estimation of marker effects and variance components using MCMC"
        write(STDOUT, "(a)") ""
      else
        write(STDOUT, "(a)") ""
        write(STDOUT, "(a)") " Running estimation of marker effects with provided variance components using MCMC"
        write(STDOUT, "(a)") ""
      end if

      BvSampFmt="("//trim(Int2Char(nRecTr))//trim("f)")
      if (nCovar.gt.0) then
        CovarSampFmt="("//trim(Int2Char(nCovar))//trim("f)")
      end if
      if (nFixed.gt.0) then
        FixedSampFmt="("//trim(Int2Char(nFixed))//trim("f)")
      end if
      if (nRandom.gt.0) then
        RandomSampFmt="("//trim(Int2Char(nRandom))//trim("f)")
      end if
      if (nMar.gt.0) then
        MarSampFmt="("//trim(Int2Char(nMar))//trim("f)")
      end if
      allocate(BvSamp(nRecTr))
      allocate(GaussDevMu(nIter))
      if (nCovar.gt.0) then
        allocate(CovarSamp(nCovar))
      end if
      if (nFixed.gt.0) then
        allocate(FixedSamp(nFixed))
      end if
      if (nRandom.gt.0) then
        allocate(RandomSamp(nRandom))
      end if
      if (nMar.gt.0) then
        allocate(MarSamp(nMar))
        allocate(GaussDevMar(nMar))
      end if
      if (EstimateVariances) then
        allocate(SX2(nMar))
        allocate(MX2(nMar))
        if (nMar.gt.0) then
          allocate(GammaDevMar(nIter))
        end if
        allocate(GammaDevErr(nIter))
      end if
      if (nMar.gt.0.and.nGenoPart.gt.0) then
        allocate(GenoPartUnit(nGenoPart+1))
        allocate(BvSampGenoPart(nRecTr))
      end if

      ! Initialise
      MuEst=0.0d0
      MuSamp=0.0d0
      if (nCovar.gt.0) then
        CovarEst=0.0d0
        CovarSamp=0.0d0
      end if
      if (nFixed.gt.0) then
        FixedEst=0.0d0
        FixedSamp=0.0d0
      end if
      if (nRandom.gt.0) then
        RandomEst=0.0d0
        RandomSamp=0.0d0
        VarRandomEst=0.0d0
      end if
      if (nMar.gt.0) then
        MarEst=0.0d0
        MarSamp=0.0d0
        VarMarEst=0.0d0
      end if
      VarErrEst=0.0d0

      if (EstimateVariances) then
        ! These prior parameters are modelled as in BGLR
        ! R2=0.5d0 with error and markers
        TmpR = 1.0d0
        if (nRandom.gt.0) then
          TmpR = TmpR + 1.0d0
        end if
        if (nMar.gt.0) then
          TmpR = TmpR + 1.0d0
        end if
        R2 = 1/TmpR

        TmpR=1 ! =Var(Err) ! should be 1 when Phen is standardized

        ErrDF0=5.0d0
        ErrDF=nRecTrR+ErrDF0
        VarErrSamp=TmpR*R2
        ErrS0=VarErrSamp*(ErrDF0+2.0d0)

        if (nRandom.gt.0) then
          RandomDF0=5.0d0
          RandomDF=nRecTrR+ErrDF0
          VarRandomSamp=TmpR*R2
          RandomS0=VarRandomSamp*(RandomDF0+2.0d0)
        end if

        if (nMar.gt.0) then
          MarDF0=5.0d0
          MarDF=nMarR+MarDF0
          SX2=0.0d0
          MX2=0.0d0
          do j=1,nMar
            SX2(j)=sum(GenosTr(:,j)*GenosTr(:,j))
            MX2(j)=Mean(GenosTr(:,j))
            MX2(j)=MX2(j)*MX2(j)
          end do
          MSX=sum(SX2)/nRecTrR-sum(MX2)
          VarMarSamp=TmpR*R2/MSX
          MarS0=VarMarSamp*(MarDF0+2.0d0)
        end if
      else
        if (nRandom.gt.0) then
          VarRandomSamp=VarRandom/PhenVar
        end if
        if (nMar.gt.0) then
          VarMarSamp=(VarBv/nMarR)/PhenVar
        end if
        VarErrSamp=VarErr/PhenVar
      end if

      ! Gauss and Gamma deviates
      GaussDevMu=SampleIntelGaussD(n=nIter)
      if (EstimateVariances) then
        GammaDevErr=SampleIntelGammaD(n=nIter,shape=ErrDF/2.0d0,scale=2.0d0)
        if (nRandom.gt.0) then
          GammaDevRandom=SampleIntelGammaD(n=nIter,shape=RandomDF/2.0d0,scale=2.0d0)
        end if
        if (nMar.gt.0) then
          GammaDevMar=SampleIntelGammaD(n=nIter,shape=MarDF/2.0d0,scale=2.0d0)
        end if
      end if

      ! Samples files
      open(newunit=MuUnit,file=INTERCEPTSAMPLESFILE,status="unknown")
      if (nCovar.gt.0) then
        open(newunit=CovarUnit,file=COVARIATESAMPLESFILE,status="unknown")
      end if
      if (nFixed.gt.0) then
        open(newunit=FixedUnit,file=FIXEDEFFECTSAMPLESFILE,status="unknown")
      end if
      if (nRandom.gt.0) then
        open(newunit=RandomUnit,file=RANDOMEFFECTSAMPLESFILE,status="unknown")
      end if
      if (nMar.gt.0) then
        open(newunit=MarUnit,file=MARKERSAMPLESFILE,status="unknown")
      end if
      if (EstimateVariances) then
        open(newunit=VarErrUnit,file=RESIDUALVARIANCESAMPLESFILE,status="unknown")
        if (nRandom.gt.0) then
          open(newunit=VarRandomUnit,file=RANDOMEFFECTVARIANCESAMPLESFILE,status="unknown")
        end if
        if (nMar.gt.0) then
          open(newunit=VarMarUnit,file=MARKERVARIANCESAMPLESFILE,status="unknown")
        end if
      end if
      if (nMar.gt.0.and.nGenoPart.gt.0) then
        do i=1,nGenoPart
          open(newunit=GenoPartUnit(i),file=trim(GENOPARTFILESTART)//trim(Int2Char(i))//trim(GENOPARTSAMPLESFILEEND),status="unknown")
        end do
        open(newunit=GenoPartUnit(i),file=trim(GENOPARTFILESTART)//trim("Remainder")//trim(GENOPARTSAMPLESFILEEND),status="unknown")
      end if

      ! Iterate
      nSampR=dble(nIter-nBurn)
      do Iter=1,nIter
        ! Intercept
        Lhs=nRecTrR/VarErrSamp
        Rhs=sum(Err)/VarErrSamp + nRecTrR*MuSamp/VarErrSamp
        Sol=Rhs/Lhs + GaussDevMu(Iter)/sqrt(Lhs)
        Diff=Sol-MuSamp
        Err=Err-Diff
        MuSamp=Sol
        if (Iter.gt.nBurn) then
          write(MuUnit,*) MuSamp*PhenSd + PhenMean
          MuEst=MuEst+MuSamp/nSampR
        end if

        ! Covariates
        if (nCovar.gt.0) then
          RandomOrdering(1:nCovar)=RandomOrder(nCovar)
          GaussDevCovar=SampleIntelGaussD(n=nCovar)
          do j=1,nCovar
            k=RandomOrdering(j)
            Lhs=XpXCovar(k)/VarErrSamp
            Rhs=dot(x=CovarTr(:,k),y=Err)/VarErrSamp + XpXCovar(k)*CovarSamp(k)/VarErrSamp
            Sol=Rhs/Lhs + GaussDevCovar(j)/sqrt(Lhs)
            Diff=Sol-CovarSamp(k)
            Err=Err-CovarTr(:,k)*Diff
            CovarSamp(k)=Sol
          end do
          if (Iter.gt.nBurn) then
            write(CovarUnit,CovarSampFmt) CovarSamp/CovarSd*PhenSd
            CovarEst=CovarEst+CovarSamp/nSampR
          end if
        end if

        ! Fixed effects
        if (nFixed.gt.0) then
          RandomOrdering(1:nFixed)=RandomOrder(nFixed)
          GaussDevFixed=SampleIntelGaussD(n=nFixed)
          do j=1,nFixed
            k=RandomOrdering(j)
            Lhs=XpXFixed(k)/VarErrSamp
            Rhs=dot(x=FixedTr(:,k),y=Err)/VarErrSamp + XpXFixed(k)*FixedSamp(k)/VarErrSamp
            Sol=Rhs/Lhs + GaussDevFixed(j)/sqrt(Lhs)
            Diff=Sol-FixedSamp(k)
            Err=Err-FixedTr(:,k)*Diff
            FixedSamp(k)=Sol
          end do
          if (Iter.gt.nBurn) then
            write(FixedUnit,FixedSampFmt) FixedSamp*PhenSd
            FixedEst=FixedEst+FixedSamp/nSampR
          end if
        end if

        ! Random effects
        if (nRandom.gt.0) then
          RandomOrdering(1:nRandom)=RandomOrder(nRandom)
          GaussDevRandom=SampleIntelGaussD(n=nRandom)
          do j=1,nRandom
            k=RandomOrdering(j)
            Lhs=XpXRandom(k)/VarErrSamp + 1.0d0/VarRandomSamp
            Rhs=dot(x=RandomTr(:,k),y=Err)/VarErrSamp + XpXRandom(k)*RandomSamp(k)/VarErrSamp
            Sol=Rhs/Lhs + GaussDevRandom(j)/sqrt(Lhs)
            Diff=Sol-RandomSamp(k)
            Err=Err-RandomTr(:,k)*Diff
            RandomSamp(k)=Sol
          end do
          if (Iter.gt.nBurn) then
            write(RandomUnit,RandomSampFmt) RandomSamp*PhenSd
            RandomEst=RandomEst+RandomSamp/nSampR
          end if
          ! Random effect variance
          if (EstimateVariances) then
            VarRandomSamp=(dot(x=RandomSamp,y=RandomSamp)+RandomS0)/GammaDevRandom(Iter)
            if (Iter.gt.nBurn) then
              write(VarRandomUnit,*) VarRandomSamp*PhenVar
              VarRandomEst=VarRandomEst+VarRandomSamp/nSampR
            end if
          end if
        end if

        ! Markers
        if (nMar.gt.0) then
          RandomOrdering(1:nMar)=RandomOrder(nMar)
          GaussDevMar=SampleIntelGaussD(n=nMar)
          do j=1,nMar
            k=RandomOrdering(j)
            Lhs=XpXMar(k)/VarErrSamp + 1.0d0/VarMarSamp
            Rhs=dot(x=GenosTr(:,k),y=Err)/VarErrSamp + XpXMar(k)*MarSamp(k)/VarErrSamp
            Sol=Rhs/Lhs + GaussDevMar(j)/sqrt(Lhs)
            Diff=Sol-MarSamp(k)
            Err=Err-GenosTr(:,k)*Diff
            MarSamp(k)=Sol
          end do
          if (Iter.gt.nBurn) then
            write(MarUnit,MarSampFmt) MarSamp/MarSd*PhenSd
            MarEst=MarEst+MarSamp/nSampR
            ! Genome partitions
            if (nGenoPart.gt.0) then
              call gemv(A=GenosTr,x=MarSamp/MarSd*PhenSd,y=BvSamp)
              do i=1,nGenoPart
                BvSampGenoPart=0.0d0
                ! Might have used gemv here, but for MCMC this would mean calling gemv nIter*nGenoPart times!!!
                do j=1,nMarPerGenoPart(i)
                  k=MarPerGenoPart(j,i)
                  BvSampGenoPart=BvSampGenoPart+(GenosTr(:,k)*MarSamp(k)/MarSd(k)*PhenSd)
                end do
                BvSamp=BvSamp-BvSampGenoPart
                write(GenoPartUnit(i),BvSampFmt) BvSampGenoPart
              end do
              write(GenoPartUnit(i),BvSampFmt) BvSamp
            end if
          end if
          ! Marker variance
          if (EstimateVariances) then
            VarMarSamp=(dot(x=MarSamp,y=MarSamp)+MarS0)/GammaDevMar(Iter)
            if (Iter.gt.nBurn) then
              write(VarMarUnit,*) VarMarSamp*PhenVar
              VarMarEst=VarMarEst+VarMarSamp/nSampR
            end if
          end if
        end if

        ! Recompute working residuals to avoid rounding errors
        if (mod(Iter,100).eq.0) then
          Err=PhenoTr-MuSamp
          if (nCovar.gt.0) then
            call gemv(A=CovarTr,x=CovarSamp,y=BvSamp)
            Err=Err-BvSamp
          end if
          if (nFixed.gt.0) then
            call gemv(A=FixedTr,x=FixedSamp,y=BvSamp)
            Err=Err-BvSamp
          end if
          if (nRandom.gt.0) then
            call gemv(A=RandomTr,x=RandomSamp,y=BvSamp)
            Err=Err-BvSamp
          end if
          if (nMar.gt.0) then
            call gemv(A=GenosTr,x=MarSamp,y=BvSamp)
            Err=Err-BvSamp
          end if
        end if
        ! Residual variance
        if (EstimateVariances) then
          VarErrSamp=(dot(x=Err,y=Err)+ErrS0)/GammaDevErr(Iter)
          if (Iter.gt.nBurn) then
            write(VarErrUnit,*) VarErrSamp*PhenVar
            VarErrEst=VarErrEst+VarErrSamp/nSampR
          end if
        end if
      end do

      flush(MuUnit)
      close(MuUnit)
      if (nCovar.gt.0) then
        flush(CovarUnit)
        close(CovarUnit)
      end if
      if (nFixed.gt.0) then
        flush(FixedUnit)
        close(FixedUnit)
      end if
      if (nRandom.gt.0) then
        flush(RandomUnit)
        close(RandomUnit)
      end if
      if (nMar.gt.0) then
        flush(MarUnit)
        close(MarUnit)
      end if
      if (EstimateVariances) then
        if (nRandom.gt.0) then
          flush(VarRandomUnit)
          close(VarRandomUnit)
        end if
        if (nMar.gt.0) then
          flush(VarMarUnit)
          close(VarMarUnit)
        end if
        flush(VarErrUnit)
        close(VarErrUnit)
      end if
      if (nMar.gt.0.and.nGenoPart.gt.0) then
        do i=1,nGenoPart+1
          flush(GenoPartUnit(i))
          close(GenoPartUnit(i))
        end do
      end if

      ! Output posterior means
      open(newunit=Unit,file=INTERCEPTESTIMATEFILE,status="unknown")
      write(Unit,*) MuEst*PhenSd + PhenMean
      flush(Unit)
      close(Unit)

      if (nCovar.gt.0) then
        open(newunit=Unit,file=COVARIATEESTIMATEFILE,status="unknown")
        do i=1,nCovar
          write(Unit,*) CovarEst(i)/CovarSd(i)*PhenSd
        end do
        flush(Unit)
        close(Unit)
      end if

      if (nFixed.gt.0) then
        open(newunit=Unit,file=FIXEDEFFECTESTIMATEFILE,status="unknown")
        do i=1,nFixed
          write(Unit,*) FixedEst(i)*PhenSd
        end do
        flush(Unit)
        close(Unit)
      end if

      if (nRandom.gt.0) then
        open(newunit=Unit,file=RANDOMEFFECTESTIMATEFILE,status="unknown")
        do i=1,nRandom
          write(Unit,*) RandomEst(i)*PhenSd
        end do
        flush(Unit)
        close(Unit)
      end if

      if (nMar.gt.0) then
        open(newunit=Unit,file=MARKERESTIMATEFILE,status="unknown")
        do i=1,nMar
          write(Unit,*) MarEst(i)/MarSd(i)*PhenSd
        end do
        flush(Unit)
        close(Unit)
      end if

      if (EstimateVariances) then
        open(newunit=Unit,file=RESIDUALVARIANCEESTIMATEFILE,status="unknown")
        write(Unit,*) VarErrEst*PhenVar
        flush(Unit)
        close(Unit)

        if (nRandom.gt.0) then
          open(newunit=Unit,file=RANDOMEFFECTVARIANCEESTIMATEFILE,status="unknown")
          write(Unit,*) VarRandomEst*PhenVar
          flush(Unit)
          close(Unit)
        end if

        if (nMar.gt.0) then
          open(newunit=Unit,file=MARKERVARIANCEESTIMATEFILE,status="unknown")
          write(Unit,*) VarMarEst*PhenVar
          flush(Unit)
          close(Unit)
        end if
      end if

      if (nMar.gt.0.and.nGenoPart.gt.0) then
        ! These are not BvEst, not Samp, but we simply reuse the already allocated variable
        call gemv(A=GenosTr,x=MarEst/MarSd*PhenSd,y=BvSamp)
        do i=1,nGenoPart
          BvSampGenoPart=0.0d0
          ! Might have used gemv here, but for MCMC this would mean calling gemv nIter*nGenoPart times!!!
          do j=1,nMarPerGenoPart(i)
            k=MarPerGenoPart(j,i)
            BvSampGenoPart=BvSampGenoPart+(GenosTr(:,k)*MarEst(k)/MarSd(k)*PhenSd)
          end do
          BvSamp=BvSamp-BvSampGenoPart
          open(newunit=Unit,file=trim(GENOPARTFILESTART)//trim(Int2Char(i))//trim(GENOPARTESTIMATEFILEEND),status="unknown")
          do j=1,nRecTr
            write(Unit,*) BvSampGenoPart(j)
          end do
          flush(Unit)
          close(Unit)
        end do
        open(newunit=Unit,file=trim(GENOPARTFILESTART)//trim("Remainder")//trim(GENOPARTESTIMATEFILEEND),status="unknown")
        do j=1,nRecTr
          write(Unit,*) BvSamp(j)
        end do
        flush(Unit)
        close(Unit)
      end if

      deallocate(BvSamp)
      deallocate(GaussDevMu)
      if (nCovar.gt.0) then
        deallocate(CovarSamp)
        deallocate(GaussDevCovar)
      end if
      if (nFixed.gt.0) then
        deallocate(FixedSamp)
        deallocate(GaussDevFixed)
      end if
      if (nRandom.gt.0) then
        deallocate(RandomSamp)
        deallocate(GaussDevRandom)
      end if
      if (nMar.gt.0) then
        deallocate(MarSamp)
        deallocate(GaussDevMar)
      end if
      if (EstimateVariances) then
        deallocate(SX2)
        deallocate(MX2)
        deallocate(GammaDevErr)
        if (nRandom.gt.0) then
          deallocate(GammaDevRandom)
        end if
        if (nMar.gt.0) then
          deallocate(GammaDevMar)
        end if
      end if
      if (nMar.gt.0.and.nGenoPart.gt.0) then
        deallocate(GenoPartUnit)
        deallocate(BvSampGenoPart)
      end if
    end subroutine

    !###########################################################################

    subroutine Prediction
      implicit none

      integer(int32) :: i,j,Pop,Unit,SummaryUnit,GenoTeUnit,PhenoTeUnit

      real(real64),allocatable :: PhenoTe(:),GenosTe(:,:)

      character(len=100) :: DumC,File
      character(len=20),allocatable :: IdTe(:)

      type(CorrelationReal64) :: Cors

      if (nMar.gt.0) then

        open(newunit=SummaryUnit,file="PredictionsSummary.txt",status="unknown")
        write(SummaryUnit,"(a14,3(1x,a12),3(1x,a20))") "Set","CorsObsEst","SlopeObsEst","SlopeEstObs","CovObsEst","VarObs","VarEst"

        open(newunit=Unit,file="PredictionsForSetTrain.txt",status="unknown")
        call gemv(A=GenosTr,x=MarEst,y=BvEst)
        PhenoTr=PhenMean+PhenoTr*PhenSd
        write(Unit,"(a20,2(1x,a20))") "Id","Observed","Estimate"
        do i=1,nRecTr
          write(Unit,"(a20,2(1x,f20.10))") IdTr(i),PhenoTr(i),BvEst(i)
        end do
        flush(Unit)
        close(Unit)

        Cors=Cor(PhenoTr,BvEst)
        DumC="Train"
        write(SummaryUnit,"(a14,3(1x,f12.4),3(1x,f20.10))") adjustl(DumC),Cors%Cor,Cors%Cov/Cors%Var2,Cors%Cov/Cors%Var1,Cors%Cov,Cors%Var1,Cors%Var2

        deallocate(BvEst)

        MarEst=MarEst/MarSd*PhenSd

        do Pop=1,nTePop
          if (nRecTe(Pop).eq.0) then
            nRecTe(Pop) = CountLines(PhenoTeFile(Pop))
          end if
        end do

        i=maxval(nRecTe(:))
        allocate(IdTe(i))
        allocate(GenosTe(i,nMar))
        allocate(PhenoTe(i))
        allocate(BvEst(i))

        do Pop=1,nTePop

          open(newunit=GenoTeUnit,file=trim(GenoTeFile(Pop)),status="old")
          open(newunit=PhenoTeUnit,file=trim(PhenoTeFile(Pop)),status="old")

          do i=1,nRecTe(Pop)
            read(GenoTeUnit,*) IdTe(i),GenosTe(i,:)
            read(PhenoTeUnit,*) DumC,PhenoTe(i)
            if (trim(DumC).ne.trim(IdTe(i))) then
              write(STDERR,"(a,i)") " ERROR: Individual identifications in the genotype and phenotype files do not match in prediction set ",Pop
              write(STDERR,"(a,i)") " ERROR: Line: ",i
              write(STDERR,"(a,a)") " ERROR: Genotype file identification: ",trim(DumC)
              write(STDERR,"(a,a)") " ERROR: Phenotype file identification: ",trim(IdTr(i))
              write(STDERR,"(a)") " "
              stop 1
            end if
          end do

          do j=1,nMar
            ! Fix any odd data
            do i=1,nRecTe(Pop)
              if ((GenosTe(i,j).lt.0.0).or.(GenosTe(i,j).gt.2.0)) then
                GenosTe(i,j)=2.0d0*AlleleFreq(j)
              end if
            end do

            ! Standardize

            ! ... center
            GenosTe(:,j)=GenosTe(:,j)-(2.0d0*AlleleFreq(j))

            ! ... scale
            ! no need because we scaled the marker solutions
          end do

          close(GenoTeUnit)
          close(PhenoTeUnit)

          File="PredictionsForSetPredict"//Int2Char(Pop)//".txt"
          open(newunit=Unit,file=trim(File),status="unknown")
          call gemv(A=GenosTe(1:nRecTe(Pop),:),x=MarEst,y=BvEst(1:nRecTe(Pop)))
          write(Unit,"(a20,2(1x,a20))") "Id","Observed","Estimate"
          do i=1,nRecTe(Pop)
            write(Unit,"(a20,2(1x,f20.10))") IdTe(i),PhenoTe(i),BvEst(i)
          end do
          flush(Unit)
          close(Unit)

          Cors=cor(PhenoTe(1:nRecTe(Pop)),BvEst(1:nRecTe(Pop)))
          DumC="Predict"//Int2Char(Pop)
          write(SummaryUnit,"(a14,3(1x,f12.4),3(1x,f20.10))") adjustl(DumC),Cors%Cor,Cors%Cov/Cors%Var2,Cors%Cov/Cors%Var1,Cors%Cov,Cors%Var1,Cors%Var2

          ! open(newunit=GenoTeUnit,file="GenoTest"//Int2Char(Pop)//"Processed.txt",status="unknown")
          ! open(newunit=PhenoTeUnit,file="PhenoTest"//Int2Char(Pop)//"Processed.txt",status="unknown")
          ! do i=1,nRecTe(Pop)
          !   write(GenoTeUnit,"("//Int2Char(nMar)//"(1x,f20.16))") GenosTe(i,:)
          !   write(PhenoTeUnit,*) PhenoTe(i)
          ! end do
          ! flush(GenoTeUnit)
          ! close(GenoTeUnit)
          ! flush(PhenoTeUnit)
          ! close(PhenoTeUnit)
        end do

        deallocate(IdTe)
        deallocate(GenosTe)
        deallocate(PhenoTe)
        deallocate(BvEst)

        flush(SummaryUnit)
        close(SummaryUnit)

      end if
    end subroutine

    !###########################################################################
end module

!###############################################################################

program AlphaBayes

  use ISO_Fortran_Env, STDIN=>input_unit,STDOUT=>output_unit,STDERR=>error_unit
  use AlphaHouseMod, only : PrintElapsedTime
  use AlphaBayesModule
  use IntelRNGMod
  implicit none

  real(real32) :: StartTime, EndTime

  call cpu_time(StartTime)
  call AlphaBayesTitle
  call IntitialiseIntelRNG
  call ReadParam
  call ReadData
  call Analysis
  call Prediction
  call UnintitialiseIntelRNG
  call cpu_time(EndTime)
  !call AlphaBayesTitle
  call PrintElapsedTime(StartTime, EndTime)

end program

!###############################################################################
