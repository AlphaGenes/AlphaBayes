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
  integer(int32) :: nMar,nAnisTr,nIter,nBurn,nProcessor,GenoScaleMethod,nGenoPart,nTePop
  integer(int32),allocatable :: nAnisTe(:),nMarPerGenoPart(:),MarPerGenoPart(:,:)

  real(real64) :: MeanY,VarY,SdY,VarBv,VarMar,VarErr,MuEst,nMarR,nAnisTrR,ExpVarX,SumExpVarX,ObsVarX,SumObsVarX,EpsTolerance
  real(real64),allocatable :: AlleleFreq(:),ScaleCoef(:),PhenoTr(:),MarEst(:),GenosTr(:,:),Err(:),XpX(:)

  character(len=FILELENGTH) :: GenoTrFile,PhenoTrFile,EstimationMethod
  character(len=FILELENGTH),allocatable :: GenoTeFile(:),PhenoTeFile(:),GenoPartFile(:)
  character(len=IDLENGTH),allocatable :: IdTr(:)

  logical :: EstimateVariances

  ! Module parameters
  REAL(REAL64),PARAMETER :: ONER=1.0d0,ZEROR=0.0d0

  CHARACTER(LEN=100),PARAMETER :: SPECFILE="AlphaBayesSpec.txt"

  CHARACTER(LEN=100),PARAMETER :: INTERCEPTESTIMATEFILE="InterceptEstimate.txt"
  CHARACTER(LEN=100),PARAMETER :: INTERCEPTSAMPLESFILE="InterceptSamples.txt"

  CHARACTER(LEN=100),PARAMETER :: MARKERESTIMATEFILE="MarkerEstimate.txt"
  CHARACTER(LEN=100),PARAMETER :: MARKERSAMPLESFILE="MarkerSamples.txt"

  CHARACTER(LEN=100),PARAMETER :: RESIDUALVARIANCEESTIMATEFILE="ResidualVarianceEstimate.txt"
  CHARACTER(LEN=100),PARAMETER :: RESIDUALVARIANCESAMPLESFILE="ResidualVarianceSamples.txt"
  CHARACTER(LEN=100),PARAMETER :: MARKERVARIANCEESTIMATEFILE="MarkerVarianceEstimate.txt"
  CHARACTER(LEN=100),PARAMETER :: MARKERVARIANCESAMPLESFILE="MarkerVarianceSamples.txt"

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

      integer(int32) :: SpecUnit, Stat, GTePop, PTePop, nAnisTeI, GenoPart
      logical :: LogStdoutInternal
      character(len=:), allocatable :: DumString
      character(len=SPECOPTIONLENGTH) :: Line
      character(len=SPECOPTIONLENGTH) :: First
      character(len=SPECOPTIONLENGTH), allocatable, dimension(:) :: Second

      open(newunit=SpecUnit, file=SPECFILE, status="old")

      ! Defaults
      LogStdoutInternal = .true.
      GenoTrFile = ""
      PhenoTrFile = ""
      nTePop = 0
      GTePop = 0
      PTePop = 0
      nAnisTr = 0
      nAnisTeI = 0
      nIter = 10000
      nBurn = 1000
      VarBv = 0.5
      VarErr = 0.5
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
            case ("genotypetrainfile")
              if (allocated(Second)) then
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

            case ("phenotypetrainfile")
              if (allocated(Second)) then
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

            case ("numberofpredictsets")
              if (allocated(Second)) then
                if (ToLower(trim(adjustl(Second(1)))) .ne. "0") then
                  nTePop = Char2Int(trim(adjustl(Second(1))))
                  allocate(PhenoTeFile(nTePop))
                  allocate(GenoTeFile(nTePop))
                  allocate(nAnisTe(nTePop))
                  nAnisTe = 0
                  if (LogStdoutInternal) then
                    write(STDOUT, "(a)") " Number of prediction sets: "//trim(Int2Char(nTePop))
                  end if
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfPredictSets, i.e., NumberOfPredictSets, 2"
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

            case ("numberoftrainrecords")
              if (allocated(Second)) then
                nAnisTr = Char2Int(trim(adjustl(Second(1))))
                nAnisTrR = dble(nAnisTr)
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Number of train records: "//trim(Int2Char(nAnisTr))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfTrainRecords, i.e., NumberOfTrainRecords, 100"
                write(STDERR, "(a)") " "
                stop 1
              end if

            case ("numberofpredictrecords")
              if (allocated(Second)) then
                nAnisTeI = nAnisTeI + 1
                if (nAnisTeI .gt. nTePop) then
                  write(STDERR, "(a)") " ERROR: Too many NumberOfPredictRecords specifications vs. NumberOfPredictSets"
                  write(STDERR, "(a)") " "
                  stop 1
                end if
                nAnisTe(nAnisTeI) = Char2Int(trim(adjustl(Second(1))))
                if (LogStdoutInternal) then
                  write(STDOUT, "(a)") " Number of prediction records: "//trim(Int2Char(nAnisTe(nAnisTeI)))
                end if
              else
                write(STDERR, "(a)") " ERROR: Must specify a number for NumberOfPredictRecords, i.e., NumberOfPredictRecords, 100"
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

      if (trim(GenoTrFile) .eq. "") then
        write(STDERR, "(a)") " ERROR: Must specify a genotype train file, i.e., GenotypeTrainFile, GenotypesTrain.txt"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (trim(PhenoTrFile) .eq. "") then
        write(STDERR, "(a)") " ERROR: Must specify a phenotype train file, i.e., PhenotypeTrainFile, PhenotypesTrain.txt"
        write(STDERR, "(a)") " "
        stop 1
      end if

      if (nMar .eq. 0) then
        write(STDERR, "(a)") " ERROR: Must specify the number of markers, i.e., NumberOfMarkers, 1000"
        write(STDERR, "(a)") " "
        stop 1
      end if

      call OMP_SET_NUM_THREADS(nProcessor)
    end subroutine

    !###########################################################################

    subroutine ReadData
      implicit none

      integer(int32) :: i,j,nNotMissing,GenoTrUnit,PhenoTrUnit,Unit!,AlleleFreqUnit

      character(len=100) :: DumC

      real(real64) :: TmpStat

      if (nAnisTr.eq.0) then
        nAnisTr = CountLines(PhenoTrFile)
        nAnisTrR = dble(nAnisTr)
      end if

      open(newunit=GenoTrUnit,file=trim(GenoTrFile),status="old")
      open(newunit=PhenoTrUnit,file=trim(PhenoTrFile),status="old")
      !open(newunit=AlleleFreqUnit,file="AlleleFreq.txt",status="unknown")
      !write(AlleleFreqUnit,"(a11,1x,a11)") "Mar","AlleleFreq"

      allocate(GenosTr(nAnisTr,nMar))
      allocate(PhenoTr(nAnisTr))
      allocate(MarEst(nMar))
      allocate(AlleleFreq(nMar))
      allocate(ScaleCoef(nMar))
      allocate(IdTr(nAnisTr))

      do i=1,nAnisTr
        read(GenoTrUnit,*) DumC,GenosTr(i,:)
        read(PhenoTrUnit,*) IdTr(i),PhenoTr(i)
        if (trim(DumC).ne.trim(IdTr(i))) then
          write(STDERR,"(a)") " ERROR: Individual identifications in the genotype and phenotype files do not match in the training set"
          write(STDERR,"(a,i)") " ERROR: Line: ",i
          write(STDERR,"(a,a)") " ERROR: Genotype file identification: ",trim(DumC)
          write(STDERR,"(a,a)") " ERROR: Phenotype file identification: ",trim(IdTr(i))
          write(STDERR,"(a)") " "
          stop 1
        end if
      end do

      MeanY=Mean(PhenoTr)
      VarY=Var(PhenoTr,MeanY)
      SdY=sqrt(VarY)
      PhenoTr=(PhenoTr-MeanY)/SdY

      SumExpVarX=0.0d0
      SumObsVarX=0.0d0
      do j=1,nMar

        ! Compute allele freqs
        AlleleFreq(j)=0.0d0
        nNotMissing=0
        do i=1,nAnisTr
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
        do i=1,nAnisTr
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
        ScaleCoef(j)=1.0d0

        if (GenoScaleMethod.eq.2) then
          ! Scale by marker specific variance - expected
          ExpVarX=sqrt(ExpVarX)
          ScaleCoef(j)=ExpVarX
          GenosTr(:,j)=GenosTr(:,j)/ScaleCoef(j)
        end if

        if (GenoScaleMethod.eq.3) then
          ! Scale by marker specific variance - observed
          ObsVarX=sqrt(ObsVarX)
          ScaleCoef(j)=ObsVarX
          GenosTr(:,j)=GenosTr(:,j)/ScaleCoef(j)
        end if

      end do

      if (GenoScaleMethod.eq.4) then
        ! Scale by average marker variance - expected
        ExpVarX=sqrt(SumExpVarX/nMarR)
        ScaleCoef(:)=ExpVarX
        GenosTr(:,:)=GenosTr(:,:)/ScaleCoef(1)
      end if

      if (GenoScaleMethod.eq.5) then
        ! Scale by average marker variance - observed
        ObsVarX=sqrt(SumObsVarX/nMarR)
        ScaleCoef(:)=ObsVarX
        GenosTr(:,:)=GenosTr(:,:)/ScaleCoef(1)
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

      close(GenoTrUnit)
      close(PhenoTrUnit)
      ! flush(AlleleFreqUnit)
      !close(AlleleFreqUnit)

      ! open(newunit=GenoTrUnit,file="GenoTrainProcessed.txt",status="unknown")
      ! open(newunit=PhenoTrUnit,file="PhenoTrainProcessed.txt",status="unknown")
      ! do i=1,nAnisTr
      !   write(GenoTrUnit,"("//Int2Char(nMar)//"(1x,f20.16))") GenosTr(i,:)
      !   write(PhenoTrUnit,*) PhenoTr(i)
      ! end do
      ! flush(GenoTrUnit)
      ! close(GenoTrUnit)
      ! flush(PhenoTrUnit)
      ! close(PhenoTrUnit)

      if (nGenoPart.gt.0) then
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
      allocate(Err(nAnisTr))
      allocate(XpX(nMar))

      ! Construct XpX
      do j=1,nMar
        XpX(j)=dot(x=GenosTr(:,j),y=GenosTr(:,j)) + tiny(GenosTr(1,1))
      end do
    end subroutine

    !###########################################################################

    subroutine Analysis
      implicit none
      if (trim(EstimationMethod).eq."RidgeSolve") then
        write(STDOUT, "(a)") ""
        write(STDOUT, "(a)") " Running estimation of marker effects with provided variance components"
        write(STDOUT, "(a)") ""
        call RidgeRegressionSolve
      end if
      if (trim(EstimationMethod).eq."RidgeSample") then
        write(STDOUT, "(a)") ""
        write(STDOUT, "(a)") " Running estimation of marker effects and variance components with MCMC"
        write(STDOUT, "(a)") ""
        call RidgeRegressionSolve
        call RidgeRegressionSample
      end if
    end subroutine

    !###########################################################################

    subroutine RidgeRegressionSolve
      implicit none

      integer(int32) :: Iter,i,j,Mar,RandomOrdering(nMar),Unit

      real(real64) :: PreMar,Rhs,Lhs,Sol,Diff,Eps
      real(real64),allocatable :: BvEst(:),BvEstGenoPart(:)

      allocate(BvEst(nAnisTr))
      if (nGenoPart.gt.0) then
        allocate(BvEstGenoPart(nAnisTr))
      end if

      ! Working phenotypes
      Err=PhenoTr

      ! Initialise
      MuEst=0.0d0
      MarEst=0.0d0

      ! Construct XpXTauE
      XpX=XpX/VarErr ! can do it only once for all rounds!!!
      VarMar=VarBv/nMarR ! approximate variance of allele substitution effects
      PreMar=1.0d0/VarMar

      ! Iterate
      do Iter=1,nIter
        Eps=0.0d0

        ! Intercept
        Lhs=nAnisTrR/VarErr
        Rhs=sum(Err)/VarErr + nAnisTrR*MuEst/VarErr
        Sol=Rhs/Lhs
        Diff=Sol-MuEst
        Err=Err-Diff
        MuEst=Sol
        Eps=Eps+Diff*Diff

        ! Marker effects
        RandomOrdering=RandomOrder(nMar)
        do j=1,nMar
          Mar=RandomOrdering(j)
          Lhs=XpX(Mar) + PreMar
          Rhs=dot(x=GenosTr(:,Mar),y=Err)/VarErr + XpX(Mar)*MarEst(Mar)
          Sol=Rhs/Lhs
          Diff=Sol-MarEst(Mar)
          Err=Err-GenosTr(:,Mar)*Diff
          MarEst(Mar)=Sol
          Eps=Eps+Diff*Diff
        end do

        ! Recompute residuals to avoid rounding errors
        if (mod(Iter,100).eq.0) then
          call gemv(A=GenosTr,x=MarEst,y=BvEst)
          Err=PhenoTr-MuEst-BvEst
        end if

        ! Stopping criteria
        if (Eps.lt.EpsTolerance) then
          exit
        end if
      end do

      open(newunit=Unit,file=INTERCEPTESTIMATEFILE,status="unknown")
      write(Unit,*) MuEst*SdY + MeanY
      flush(Unit)
      close(Unit)

      open(newunit=Unit,file=MARKERESTIMATEFILE,status="unknown")
      do i=1,nMar
        write(Unit,*) MarEst(i)/ScaleCoef(i)*SdY
      end do
      flush(Unit)
      close(Unit)

      ! Recalculate breeding values and residuals for later use
      if (nGenoPart.gt.0.or.EstimationMethod.eq."RidgeSample") then
        call gemv(A=GenosTr,x=MarEst/ScaleCoef*SdY,y=BvEst)
        Err=PhenoTr-MuEst-BvEst
      end if

      if (nGenoPart.gt.0) then
        do i=1,nGenoPart
          BvEstGenoPart=0.0d0
          ! Might have used gemv here, but for MCMC this would mean calling gemv nIter*nGenoPart times!!!
          do j=1,nMarPerGenoPart(i)
            Mar=MarPerGenoPart(j,i)
            BvEstGenoPart=BvEstGenoPart+(GenosTr(:,Mar)*MarEst(Mar)/ScaleCoef(Mar)*SdY)
          end do
          BvEst=BvEst-BvEstGenoPart
          open(newunit=Unit,file=trim(GENOPARTFILESTART)//trim(Int2Char(i))//trim(GENOPARTESTIMATEFILEEND),status="unknown")
          do j=1,nAnisTr
            write(Unit,*) BvEstGenoPart(j)
          end do
          flush(Unit)
          close(Unit)
        end do
        open(newunit=Unit,file=trim(GENOPARTFILESTART)//trim("Remainder")//trim(GENOPARTESTIMATEFILEEND),status="unknown")
        do j=1,nAnisTr
          write(Unit,*) BvEst(j)
        end do
        flush(Unit)
        close(Unit)
      end if

      ! Convert XpXTauE=XpX/VarErr back to XpX
      XpX=XpX*VarErr
      deallocate(BvEst)
      if (nGenoPart.gt.0) then
        deallocate(BvEstGenoPart)
      end if
    end subroutine

    !###########################################################################

    subroutine RidgeRegressionSample
      implicit none

      integer(int32) :: Iter,i,j,Mar,RandomOrdering(nMar),Unit,MuUnit,MarUnit,VarMarUnit,VarErrUnit
      integer(int32),allocatable :: GenoPartUnit(:)

      real(real64) :: TmpR,Rhs,Lhs,Sol,Diff,nSampR
      real(real64) :: MuSamp,VarErrSamp,VarMarSamp,VarErrEst,VarMarEst
      real(real64) :: R2,ErrDF0,ErrDF,MarDF0,MarDF,ErrS0,MarS0,MSX
      real(real64),allocatable :: MarSamp(:)
      real(real64),allocatable :: SX2(:),MX2(:),BvSamp(:),BvSampGenoPart(:)
      real(real64),allocatable :: GaussDevMu(:),GaussDevMar(:),GammaDevErr(:),GammaDevMar(:)

      character(len=100) :: MarSampFmt, BvSampFmt

      MarSampFmt="("//trim(Int2Char(nMar))//trim("f)")
      BvSampFmt="("//trim(Int2Char(nAnisTr))//trim("f)")

      allocate(BvSamp(nAnisTr))
      allocate(MarSamp(nMar))
      allocate(GaussDevMu(nIter))
      allocate(GaussDevMar(nMar))
      if (EstimateVariances) then
        allocate(SX2(nMar))
        allocate(MX2(nMar))
        allocate(GammaDevErr(nIter))
        allocate(GammaDevMar(nIter))
      end if
      if (nGenoPart.gt.0) then
        allocate(GenoPartUnit(nGenoPart+1))
        allocate(BvSampGenoPart(nAnisTr))
      end if

      ! Initialise
      MuEst=0.0d0
      MuSamp=0.0d0
      MarEst=0.0d0
      MarSamp=0.0d0
      VarMarEst=0.0d0
      VarErrEst=0.0d0

      if (EstimateVariances) then
        ! These prior parameters are modelled as in BGLR
        R2=0.5d0

        TmpR=1 ! =Var(Err) ! should be 1 when Phen is standardized
        VarErrSamp=TmpR*(1.0d0-R2)
        ErrDF0=5.0d0
        ErrDF=nAnisTrR+ErrDF0
        ErrS0=VarErrSamp*(ErrDF0+2.0d0)

        MarDF0=5.0d0
        MarDF=nMarR+MarDF0
        SX2(:)=0.0d0
        MX2(:)=0.0d0
        do j=1,nMar
          SX2(j)=sum(GenosTr(:,j)*GenosTr(:,j))
          MX2(j)=Mean(GenosTr(:,j))
          MX2(j)=MX2(j)*MX2(j)
        end do
        MSX=sum(SX2)/nAnisTrR-sum(MX2)
        VarMarSamp=TmpR*R2/MSX
        MarS0=VarMarSamp*(MarDF0+2.0d0)
      else
        VarMarSamp=(VarBv/nMarR)/VarY
        VarErrSamp=VarErr/VarY
      end if

      ! Gauss and Gamma deviates
      GaussDevMu(:)=SampleIntelGaussD(n=nIter)
      if (EstimateVariances) then
        GammaDevErr(:)=SampleIntelGammaD(n=nIter,shape=ErrDF/2.0d0,scale=2.0d0)
        GammaDevMar(:)=SampleIntelGammaD(n=nIter,shape=MarDF/2.0d0,scale=2.0d0)
      end if

      ! Samples files
      open(newunit=MuUnit,file=INTERCEPTSAMPLESFILE,status="unknown")
      open(newunit=MarUnit,file=MARKERSAMPLESFILE,status="unknown")
      if (EstimateVariances) then
        open(newunit=VarErrUnit,file=RESIDUALVARIANCESAMPLESFILE,status="unknown")
        open(newunit=VarMarUnit,file=MARKERVARIANCESAMPLESFILE,status="unknown")
      end if
      if (nGenoPart.gt.0) then
        do i=1,nGenoPart
          open(newunit=GenoPartUnit(i),file=trim(GENOPARTFILESTART)//trim(Int2Char(i))//trim(GENOPARTSAMPLESFILEEND),status="unknown")
        end do
        open(newunit=GenoPartUnit(i),file=trim(GENOPARTFILESTART)//trim("Remainder")//trim(GENOPARTSAMPLESFILEEND),status="unknown")
      end if

      ! Iterate
      nSampR=dble(nIter-nBurn)
      do Iter=1,nIter
        ! Intercept
        Lhs=nAnisTrR/VarErrSamp
        Rhs=sum(Err)/VarErrSamp + nAnisTrR*MuSamp/VarErrSamp
        Sol=Rhs/Lhs + GaussDevMu(Iter)/sqrt(Lhs)
        Diff=Sol-MuSamp
        Err=Err-Diff
        MuSamp=Sol
        if (Iter.gt.nBurn) then
          write(MuUnit,*) MuSamp*SdY + MeanY
          MuEst=MuEst+MuSamp/nSampR
        end if

        ! Marker effects
        RandomOrdering=RandomOrder(nMar)
        GaussDevMar(:)=SampleIntelGaussD(n=nMar)
        do j=1,nMar
          Mar=RandomOrdering(j)
          Lhs=XpX(Mar)/VarErrSamp + 1.0d0/VarMarSamp
          Rhs=dot(x=GenosTr(:,Mar),y=Err)/VarErrSamp + XpX(Mar)*MarSamp(Mar)/VarErrSamp
          Sol=Rhs/Lhs + GaussDevMar(j)/sqrt(Lhs)
          Diff=Sol-MarSamp(Mar)
          Err=Err-GenosTr(:,Mar)*Diff
          MarSamp(Mar)=Sol
        end do
        if (Iter.gt.nBurn) then
          write(MarUnit,MarSampFmt) MarSamp/ScaleCoef*SdY
          MarEst=MarEst+MarSamp/nSampR
        end if

        ! Marker variance
        if (EstimateVariances) then
          VarMarSamp=(dot(x=MarSamp,y=MarSamp)+MarS0)/GammaDevMar(Iter)
          if (Iter.gt.nBurn) then
            write(VarMarUnit,*) VarMarSamp*VarY
            VarMarEst=VarMarEst+VarMarSamp/nSampR
          end if
        end if

        ! Recompute residuals to avoid rounding errors
        if (mod(Iter,100).eq.0) then
          call gemv(A=GenosTr,x=MarSamp,y=BvSamp)
          Err=PhenoTr-MuSamp-BvSamp
        end if

        ! Residual variance
        if (EstimateVariances) then
          VarErrSamp=(dot(x=Err,y=Err)+ErrS0)/GammaDevErr(Iter)
          if (Iter.gt.nBurn) then
            write(VarErrUnit,*) VarErrSamp*VarY
            VarErrEst=VarErrEst+VarErrSamp/nSampR
          end if
        end if

        ! Genome partitions
        if (Iter.gt.nBurn) then
          if (nGenoPart.gt.0) then
            call gemv(A=GenosTr,x=MarSamp/ScaleCoef*SdY,y=BvSamp)
            do i=1,nGenoPart
              BvSampGenoPart=0.0d0
              ! Might have used gemv here, but for MCMC this would mean calling gemv nIter*nGenoPart times!!!
              do j=1,nMarPerGenoPart(i)
                Mar=MarPerGenoPart(j,i)
                BvSampGenoPart=BvSampGenoPart+(GenosTr(:,Mar)*MarSamp(Mar)/ScaleCoef(Mar)*SdY)
              end do
              BvSamp=BvSamp-BvSampGenoPart
              write(GenoPartUnit(i),BvSampFmt) BvSampGenoPart
            end do
            write(GenoPartUnit(i),BvSampFmt) BvSamp
          end if
        end if

      end do

      flush(MuUnit)
      close(MuUnit)
      flush(MarUnit)
      close(MarUnit)
      if (EstimateVariances) then
        flush(VarMarUnit)
        close(VarMarUnit)
        flush(VarErrUnit)
        close(VarErrUnit)
      end if
      if (nGenoPart.gt.0) then
        do i=1,nGenoPart+1
          flush(GenoPartUnit(i))
          close(GenoPartUnit(i))
        end do
      end if

      ! Output posterior means
      open(newunit=Unit,file=INTERCEPTESTIMATEFILE,status="unknown")
      write(Unit,*) MuEst*SdY + MeanY
      flush(Unit)
      close(Unit)

      open(newunit=Unit,file=MARKERESTIMATEFILE,status="unknown")
      do i=1,nMar
        write(Unit,*) MarEst(i)/ScaleCoef(i)*SdY
      end do
      flush(Unit)
      close(Unit)

      if (EstimateVariances) then
        open(newunit=Unit,file=RESIDUALVARIANCEESTIMATEFILE,status="unknown")
        write(Unit,*) VarErrEst*VarY
        flush(Unit)
        close(Unit)

        open(newunit=Unit,file=MARKERVARIANCEESTIMATEFILE,status="unknown")
        write(Unit,*) VarMarEst*VarY
        flush(Unit)
        close(Unit)
      end if

      if (nGenoPart.gt.0) then
        ! These are not BvEst, not Samp, but we simply reuse the already allocated variable
        call gemv(A=GenosTr,x=MarEst/ScaleCoef*SdY,y=BvSamp)
        do i=1,nGenoPart
          BvSampGenoPart=0.0d0
          ! Might have used gemv here, but for MCMC this would mean calling gemv nIter*nGenoPart times!!!
          do j=1,nMarPerGenoPart(i)
            Mar=MarPerGenoPart(j,i)
            BvSampGenoPart=BvSampGenoPart+(GenosTr(:,Mar)*MarEst(Mar)/ScaleCoef(Mar)*SdY)
          end do
          BvSamp=BvSamp-BvSampGenoPart
          open(newunit=Unit,file=trim(GENOPARTFILESTART)//trim(Int2Char(i))//trim(GENOPARTESTIMATEFILEEND),status="unknown")
          do j=1,nAnisTr
            write(Unit,*) BvSampGenoPart(j)
          end do
          flush(Unit)
          close(Unit)
        end do
        open(newunit=Unit,file=trim(GENOPARTFILESTART)//trim("Remainder")//trim(GENOPARTESTIMATEFILEEND),status="unknown")
        do j=1,nAnisTr
          write(Unit,*) BvSamp(j)
        end do
        flush(Unit)
        close(Unit)
      end if

      deallocate(XpX)
      deallocate(BvSamp)
      deallocate(Err)
      deallocate(MarSamp)
      deallocate(GaussDevMu)
      deallocate(GaussDevMar)
      if (EstimateVariances) then
        deallocate(SX2)
        deallocate(MX2)
        deallocate(GammaDevErr)
        deallocate(GammaDevMar)
      end if
      if (nGenoPart.gt.0) then
        deallocate(GenoPartUnit)
        deallocate(BvSampGenoPart)
      end if
    end subroutine

    !###########################################################################

    subroutine Prediction
      implicit none

      integer(int32) :: i,j,Pop,Unit,SummaryUnit,GenoTeUnit,PhenoTeUnit

      real(real64),allocatable :: BvEst(:),PhenoTe(:),GenosTe(:,:)

      character(len=100) :: DumC,File
      character(len=20),allocatable :: IdTe(:)

      type(CorrelationReal64) :: Cors

      open(newunit=SummaryUnit,file="PredictionsSummary.txt",status="unknown")
      write(SummaryUnit,"(a14,3(1x,a12),3(1x,a20))") "Set","CorsObsEst","SlopeObsEst","SlopeEstObs","CovObsEst","VarObs","VarEst"

      allocate(BvEst(nAnisTr))

      open(newunit=Unit,file="PredictionsForSetTrain.txt",status="unknown")
      call gemv(A=GenosTr,x=MarEst,y=BvEst)
      PhenoTr=MeanY+PhenoTr*SdY
      write(Unit,"(a20,2(1x,a20))") "Id","Observed","Estimate"
      do i=1,nAnisTr
        write(Unit,"(a20,2(1x,f20.10))") IdTr(i),PhenoTr(i),BvEst(i)
      end do
      flush(Unit)
      close(Unit)

      Cors=Cor(PhenoTr,BvEst)
      DumC="Train"
      write(SummaryUnit,"(a14,3(1x,f12.4),3(1x,f20.10))") adjustl(DumC),Cors%Cor,Cors%Cov/Cors%Var2,Cors%Cov/Cors%Var1,Cors%Cov,Cors%Var1,Cors%Var2

      deallocate(BvEst)

      MarEst=MarEst/ScaleCoef(:)*SdY

      do Pop=1,nTePop
        if (nAnisTe(Pop).eq.0) then
          nAnisTe(Pop) = CountLines(PhenoTeFile(Pop))
        end if
      end do

      i=maxval(nAnisTe(:))
      allocate(IdTe(i))
      allocate(GenosTe(i,nMar))
      allocate(PhenoTe(i))
      allocate(BvEst(i))

      do Pop=1,nTePop

        open(newunit=GenoTeUnit,file=trim(GenoTeFile(Pop)),status="old")
        open(newunit=PhenoTeUnit,file=trim(PhenoTeFile(Pop)),status="old")

        do i=1,nAnisTe(Pop)
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
          do i=1,nAnisTe(Pop)
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
        call gemv(A=GenosTe(1:nAnisTe(Pop),:),x=MarEst,y=BvEst(1:nAnisTe(Pop)))
        write(Unit,"(a20,2(1x,a20))") "Id","Observed","Estimate"
        do i=1,nAnisTe(Pop)
          write(Unit,"(a20,2(1x,f20.10))") IdTe(i),PhenoTe(i),BvEst(i)
        end do
        flush(Unit)
        close(Unit)

        Cors=cor(PhenoTe(1:nAnisTe(Pop)),BvEst(1:nAnisTe(Pop)))
        DumC="Predict"//Int2Char(Pop)
        write(SummaryUnit,"(a14,3(1x,f12.4),3(1x,f20.10))") adjustl(DumC),Cors%Cor,Cors%Cov/Cors%Var2,Cors%Cov/Cors%Var1,Cors%Cov,Cors%Var1,Cors%Var2

        ! open(newunit=GenoTeUnit,file="GenoTest"//Int2Char(Pop)//"Processed.txt",status="unknown")
        ! open(newunit=PhenoTeUnit,file="PhenoTest"//Int2Char(Pop)//"Processed.txt",status="unknown")
        ! do i=1,nAnisTe(Pop)
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
