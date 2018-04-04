#ifdef BINARY
#define BINFILE ,form="unformatted"
#else
#define BINFILE
#endif

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#ifdef OS_UNIX
#define DASH "/"
#define COPY "cp"
#define MKDIR "mkdir -p"
#define RMDIR "rm -r"
#define RM "rm"
#define RENAME "mv"
#else
#define DASH "\"
#define COPY "copy"
#define MKDIR "md"
#define RMDIR "rmdir /S"
#define RM "del"
#define RENAME "move /Y"
#endif

!###############################################################################

module AlphaBayesModule

  use ISO_Fortran_Env, STDIN=>input_unit,STDOUT=>output_unit,STDERR=>error_unit
  use IntelRNGMod
  use AlphaHouseMod
  use AlphaStatMod
  use Blas95, only : dot, & ! https://software.intel.com/en-us/mkl-developer-reference-fortran-dot  dot(x,y) does x'y
                     gemv   ! https://software.intel.com/en-us/mkl-developer-reference-fortran-gemv gemv(A,x,y) does y=Ax
  use F95_precision

  implicit none

  integer(int32) :: nSnp,nAnisTr,nIter,nBurn,nProcessor,ScalingOpt,nPartition,nTePop
  integer(int32),allocatable :: nAnisTe(:),nSnpPerPartition(:),SnpPerPartition(:,:)

  real(real64) :: MeanY,VarY,SdY,VarG,VarE,Mu,nSnpR,nAnisTrR,ExpVarX,SumExpVarX,ObsVarX,SumObsVarX
  real(real64),allocatable :: AlleleFreq(:),ScaleCoef(:),PhenoTr(:),G(:),GenosTr(:,:)

  character(len=1000) :: GenoTrFile,PhenoTrFile,Method
  character(len=1000),allocatable :: GenoTeFile(:),PhenoTeFile(:)
  character(len=20),allocatable :: IdTr(:)

  REAL(REAL64),PARAMETER :: ONER=1.0d0,ZEROR=0.0d0

  CHARACTER(LEN=100),PARAMETER :: SPECFILE="AlphaBayesSpec.txt"

  CHARACTER(LEN=100),PARAMETER :: INTERCEPTESTIMATEFILE="InterceptEstimate.txt"
  CHARACTER(LEN=100),PARAMETER :: INTERCEPTSAMPLESFILE="InterceptSamples.txt"

  CHARACTER(LEN=100),PARAMETER :: SNPESTIMATEFILE="SnpEstimate.txt"
  CHARACTER(LEN=100),PARAMETER :: SNPSAMPLESFILE="SnpSamples.txt"

  CHARACTER(LEN=100),PARAMETER :: RESIDUALVARIANCEESTIMATEFILE="ResidualVarianceEstimate.txt"
  CHARACTER(LEN=100),PARAMETER :: RESIDUALVARIANCESAMPLESFILE="ResidualVarianceSamples.txt"
  CHARACTER(LEN=100),PARAMETER :: SNPVARIANCEESTIMATEFILE="SnpVarianceEstimate.txt"
  CHARACTER(LEN=100),PARAMETER :: SNPVARIANCESAMPLESFILE="SnpVarianceSamples.txt"

  private
  public :: ReadParam,ReadData,Method,Prediction
  public :: RidgeRegression,RidgeRegressionMCMC

  contains

    !###########################################################################

    subroutine ReadParam
      implicit none

      integer(int32) :: i,SpecUnit

      character(len=100) :: DumC

      open(newunit=SpecUnit,file=SPECFILE,status="old")

      read(SpecUnit,*) DumC,GenoTrFile
      read(SpecUnit,*) DumC,PhenoTrFile

      read(SpecUnit,*) DumC,nTePop
      if (nTePop.gt.0) then
        allocate(PhenoTeFile(nTePop))
        allocate(GenoTeFile(nTePop))
        allocate(nAnisTe(nTePop))
        read(SpecUnit,*) DumC,GenoTeFile(:)
        read(SpecUnit,*) DumC,PhenoTeFile(:)
      end if

      read(SpecUnit,*) DumC,nSnp
      nSnpR=dble(nSnp)

      read(SpecUnit,*) DumC,nAnisTr
      nAnisTrR=dble(nAnisTr)
      read(SpecUnit,*) DumC,nAnisTe(:)

      read(SpecUnit,*) DumC,nIter
      read(SpecUnit,*) DumC,nBurn

      read(SpecUnit,*) DumC,VarG ! here it is meant as genetic variance!!!
      read(SpecUnit,*) DumC,VarE

      read(SpecUnit,*) DumC,nProcessor
      call OMP_SET_NUM_THREADS(nProcessor)

      read(SpecUnit,*) DumC,ScalingOpt
      if (ScalingOpt.lt.1.or.ScalingOpt.gt.5) then
        write(STDERR,"(a)") "ERROR: ScalingOption must be between 1 and 5"
        write(STDERR,"(a)") " "
        stop 1
      endif

      read(SpecUnit,*) DumC,Method

      close(SpecUnit)
    end subroutine

    !###########################################################################

    subroutine ReadData
      implicit none

      integer(int32) :: i,j,nNotMissing,GenoTrUnit,PhenoTrUnit!,AlleleFreqUnit

      character(len=100) :: DumC

      real(real64) :: TmpStat

      open(newunit=GenoTrUnit,file=trim(GenoTrFile),status="old")
      open(newunit=PhenoTrUnit,file=trim(PhenoTrFile),status="old")
      !open(newunit=AlleleFreqUnit,file="AlleleFreq.txt",status="unknown")
      !write(AlleleFreqUnit,"(a11,1x,a11)") "Snp","AlleleFreq"

      allocate(GenosTr(nAnisTr,nSnp))
      allocate(PhenoTr(nAnisTr))
      allocate(G(nSnp))
      allocate(AlleleFreq(nSnp))
      allocate(ScaleCoef(nSnp))
      allocate(IdTr(nAnisTr))

      do i=1,nAnisTr
        read(GenoTrUnit,*) DumC,GenosTr(i,:)
        read(PhenoTrUnit,*) IdTr(i),PhenoTr(i)
        if (trim(DumC).ne.trim(IdTr(i))) then
          write(STDERR,"(a)") "ERROR: Individual identifications in the genotype and phenotype files do not match in the training set"
          write(STDERR,"(a,i)") "ERROR: Line: ",i
          write(STDERR,"(a,a)") "ERROR: Genotype file identification: ",trim(DumC)
          write(STDERR,"(a,a)") "ERROR: Phenotype file identification: ",trim(IdTr(i))
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
      do j=1,nSnp

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

        if (ScalingOpt.eq.2) then
          ! Scale by marker specific variance - expected
          ExpVarX=sqrt(ExpVarX)
          ScaleCoef(j)=ExpVarX
          GenosTr(:,j)=GenosTr(:,j)/ScaleCoef(j)
        end if

        if (ScalingOpt.eq.3) then
          ! Scale by marker specific variance - observed
          ObsVarX=sqrt(ObsVarX)
          ScaleCoef(j)=ObsVarX
          GenosTr(:,j)=GenosTr(:,j)/ScaleCoef(j)
        end if

      end do

      if (ScalingOpt.eq.4) then
        ! Scale by average marker variance - expected
        ExpVarX=sqrt(SumExpVarX/nSnpR)
        ScaleCoef(:)=ExpVarX
        GenosTr(:,:)=GenosTr(:,:)/ScaleCoef(1)
      end if

      if (ScalingOpt.eq.5) then
        ! Scale by average marker variance - observed
        ObsVarX=sqrt(SumObsVarX/nSnpR)
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
      !   write(GenoTrUnit,"("//Int2Char(nSnp)//"(1x,f20.16))") GenosTr(i,:)
      !   write(PhenoTrUnit,*) PhenoTr(i)
      ! end do
      ! flush(GenoTrUnit)
      ! close(GenoTrUnit)
      ! flush(PhenoTrUnit)
      ! close(PhenoTrUnit)
    end subroutine

    !###########################################################################

    subroutine RidgeRegression
      implicit none

      integer(int32) :: Iter,j,Snp,RandomOrdering(nSnp),Unit

      real(real64) :: VarS,Rhs,Lhs,Sol,Diff,Eps
      real(real64),allocatable :: XpXTauE(:),Xg(:),E(:),XgPartition(:)

      allocate(XpXTauE(nSnp))
      allocate(Xg(nAnisTr))
      allocate(E(nAnisTr))

      ! Working phenotypes
      E=PhenoTr

      ! Initialise
      Mu=0.0d0
      G=0.0d0

      ! Construct XpXTauE
      do j=1,nSnp
        XpXTauE(j)=dot(x=GenosTr(:,j),y=GenosTr(:,j)) + tiny(GenosTr(1,1))
      end do
      XpXTauE=XpXTauE/VarE ! can do it only once for all rounds!!!
      VarS=VarG/nSnpR ! approximate variance of allele substitution effects

      ! Iterate
      do Iter=1,nIter
        Eps=0.0d0

        ! Intercept
        Lhs=nAnisTrR/VarE
        Rhs=sum(E)/VarE + nAnisTrR*Mu/VarE
        Sol=Rhs/Lhs
        Diff=Sol-Mu
        E=E-Diff
        Mu=Sol
        Eps=Eps+Diff*Diff

        ! Snp effects
        RandomOrdering=RandomOrder(nSnp)
        do j=1,nSnp
          Snp=RandomOrdering(j)
          Lhs=XpXTauE(Snp) + 1.0d0/VarS
          Rhs=dot(x=GenosTr(:,Snp),y=E)/VarE + XpXTauE(Snp)*G(Snp)
          Sol=Rhs/Lhs
          Diff=Sol-G(Snp)
          E=E-GenosTr(:,Snp)*Diff
          G(Snp)=Sol
          Eps=Eps+Diff*Diff
        end do

        ! Recompute residuals to avoid rounding errors
        if (mod(Iter,100).eq.0) then
          call gemv(A=GenosTr,x=G,y=Xg)
          E=PhenoTr-Mu-Xg
        end if

        ! Stopping criteria
        if (eps.lt.1.0E-8) then
          exit
        end if
      end do

      open(newunit=Unit,file=INTERCEPTESTIMATEFILE,status="unknown")
      write(Unit,*) Mu*SdY
      flush(Unit)
      close(Unit)

      open(newunit=Unit,file=SNPESTIMATEFILE,status="unknown")
      do i=1,nSnp
        write(Unit,*) G(i)/ScaleCoef(i)*SdY
      end do
      flush(Unit)
      close(Unit)

      deallocate(XpXTauE)
      deallocate(Xg)
      deallocate(E)
    end subroutine

    !###########################################################################

    subroutine RidgeRegressionMCMC
      implicit none

      integer(int32) :: Iter,i,j,Snp,RandomOrdering(nSnp),Unit

      real(real64) :: TmpR,Rhs,Lhs,Sol,Diff,nSampR
      real(real64) :: MuSamp,VarESamp,VarGSamp
      real(real64) :: R2,EDF0,EDF,GDF0,GDF,ES0,GS0,MSX,EpE,GpG
      real(real64),allocatable :: GSamp(:),MuAll(:),GAll(:,:),VarEAll(:),VarGAll(:)
      real(real64),allocatable :: SX2(:),MX2(:),XpX(:),Xg(:),E(:)
      real(real64),allocatable :: GaussDevMu(:),GaussDevSnp(:),GammaDevE(:),GammaDevG(:)

      allocate(XpX(nSnp))
      allocate(Xg(nAnisTr))
      allocate(E(nAnisTr))
      allocate(GSamp(nSnp))
      allocate(MuAll(nIter))
      allocate(GAll(nSnp,nIter))
      allocate(VarEAll(nIter))
      allocate(VarGAll(nIter))
      allocate(SX2(nSnp))
      allocate(MX2(nSnp))
      allocate(GaussDevMu(nIter))
      allocate(GaussDevSnp(nSnp))
      allocate(GammaDevE(nIter))
      allocate(GammaDevG(nIter))

      ! Working phenotypes
      E=PhenoTr

      ! Initialise
      Mu=0.0d0
      MuSamp=0.0d0
      G=0.0d0
      GSamp=0.0d0
      VarG=0.0d0! VarG here is variance of allele substitution effects!!!
      VarGSamp=0.0d0
      VarE=0.0d0
      VarESamp=0.0d0

      ! These prior parameters are modelled as in BGLR
      R2=0.5d0

      TmpR=Var(E) ! should be 1 if Phen is standardized
      VarESamp=TmpR*(1.0d0-R2)
      EDF0=5.0d0
      EDF=nAnisTrR+EDF0
      ES0=VarESamp*(EDF0+2.0d0)

      GDF0=5.0d0
      GDF=nSnpR+GDF0
      SX2(:)=0.0d0
      MX2(:)=0.0d0
      do j=1,nSnp
        SX2(j)=sum(GenosTr(:,j)*GenosTr(:,j))
        MX2(j)=Mean(GenosTr(:,j))
        MX2(j)=MX2(j)*MX2(j)
      end do
      MSX=sum(SX2)/nAnisTrR-sum(MX2)
      VarGSamp=TmpR*R2/MSX
      GS0=VarGSamp*(GDF0+2.0d0)

      deallocate(SX2)
      deallocate(MX2)

      ! Construct XpX
      do j=1,nSnp
        XpX(j)=dot(x=GenosTr(:,j),y=GenosTr(:,j)) + tiny(GenosTr(1,1))
      end do

      ! Gauss and Gamma deviates
      GaussDevMu(:)=SampleIntelGaussD(n=nIter)
      GammaDevE(:)=SampleIntelGammaD(n=nIter,shape=EDF/2.0d0,scale=2.0d0)
      GammaDevG(:)=SampleIntelGammaD(n=nIter,shape=GDF/2.0d0,scale=2.0d0)

      ! Iterate
      nSampR=dble(nIter-nBurn)
      do Iter=1,nIter
        ! Intercept
        Lhs=nAnisTrR/VarESamp
        Rhs=sum(E)/VarESamp + nAnisTrR*MuSamp/VarESamp
        Sol=Rhs/Lhs + GaussDevMu(Iter)/sqrt(Lhs)
        Diff=Sol-MuSamp
        E=E-Diff
        MuSamp=Sol

        ! Snp effects
        RandomOrdering=RandomOrder(nSnp)
        GaussDevSnp(:)=SampleIntelGaussD(n=nSnp)
        do j=1,nSnp
          Snp=RandomOrdering(j)
          Lhs=XpX(Snp)/VarESamp + 1.0d0/VarGSamp
          Rhs=dot(x=GenosTr(:,Snp),y=E)/VarESamp + XpX(Snp)*GSamp(Snp)/VarESamp
          Sol=Rhs/Lhs + GaussDevSnp(j)/sqrt(Lhs)
          Diff=Sol-GSamp(Snp)
          E=E-GenosTr(:,Snp)*Diff
          GSamp(Snp)=Sol
        end do

        ! Snp variance
        GpG=dot(x=GSamp,y=GSamp)+GS0
        VarGSamp=GpG/GammaDevG(Iter)

        ! Recompute residuals to avoid rounding errors
        if (mod(Iter,100).eq.0) then
          call gemv(A=GenosTr,x=GSamp,y=Xg)
          E=PhenoTr-MuSamp-Xg
        end if

        ! Residual variance
        EpE=dot(x=E,y=E)+ES0
        VarESamp=EpE/GammaDevE(Iter)

        MuAll(Iter)=MuSamp
        GAll(:,Iter)=GSamp
        VarEAll(Iter)=VarESamp
        VarGAll(Iter)=VarGSamp

      end do

      ! Posterior means
      do j=nBurn+1,nIter
        Mu=Mu+MuAll(j)/nSampR
        G(:)=G(:)+GAll(:,j)/nSampR
        VarE=VarE+VarEAll(j)/nSampR
        VarG=VarG+VarGAll(j)/nSampR
      end do

      ! Outputs
      open(newunit=Unit,file=INTERCEPTSAMPLESFILE,status="unknown")
      write(Unit,*) MuAll(:)*SdY
      flush(Unit)
      close(Unit)
      open(newunit=Unit,file=INTERCEPTESTIMATEFILE,status="unknown")
      write(Unit,*) Mu*SdY
      flush(Unit)
      close(Unit)

      open(newunit=Unit,file=SNPSAMPLESFILE,status="unknown")
      do i=1,nSnp
        write(Unit,*) GAll(i,:)/ScaleCoef(i)*SdY
      end do
      flush(Unit)
      close(Unit)
      open(newunit=Unit,file=SNPESTIMATEFILE,status="unknown")
      do i=1,nSnp
        write(Unit,*) G(i)/ScaleCoef(i)*SdY
      end do
      flush(Unit)
      close(Unit)

      open(newunit=Unit,file=RESIDUALVARIANCESAMPLESFILE,status="unknown")
      write(Unit,*) VarEAll(:)*VarY
      flush(Unit)
      close(Unit)
      open(newunit=Unit,file=RESIDUALVARIANCEESTIMATEFILE,status="unknown")
      write(Unit,*) VarE*VarY
      flush(Unit)
      close(Unit)

      open(newunit=Unit,file=SNPVARIANCESAMPLESFILE,status="unknown")
      write(Unit,*) VarGAll(:)*VarY
      flush(Unit)
      close(Unit)
      open(newunit=Unit,file=SNPVARIANCEESTIMATEFILE,status="unknown")
      write(Unit,*) VarG*VarY
      flush(Unit)
      close(Unit)

      deallocate(XpX)
      deallocate(Xg)
      deallocate(E)
      deallocate(MuAll)
      deallocate(GAll)
      deallocate(VarEAll)
      deallocate(VarGAll)
    end subroutine

    !###########################################################################

    subroutine Prediction
      implicit none

      integer(int32) :: i,j,Pop,Unit,SummaryUnit,GenoTeUnit,PhenoTeUnit

      real(real64),allocatable :: Ebv(:),PhenoTe(:),GenosTe(:,:)

      character(len=100) :: DumC,File
      character(len=20),allocatable :: IdTe(:)

      type(CorrelationReal64) :: Cors

      open(newunit=SummaryUnit,file="PredictionsSummary.txt",status="unknown")
      write(SummaryUnit,"(a14,3(1x,a12),3(1x,a20))") "Set","CorsObsEst","SlopeObsEst","SlopeEstObs","CovObsEst","VarObs","VarEst"

      allocate(Ebv(nAnisTr))

      open(newunit=Unit,file="PredictionsForSetTrain.txt",status="unknown")
      call gemv(A=GenosTr,x=G,y=Ebv)
      PhenoTr=MeanY+PhenoTr*SdY
      write(Unit,"(a20,2(1x,a20))") "Id","Observed","Estimate"
      do i=1,nAnisTr
        write(Unit,"(a20,2(1x,f20.10))") IdTr(i),PhenoTr(i),Ebv(i)
      end do
      flush(Unit)
      close(Unit)

      Cors=Cor(PhenoTr,Ebv)
      DumC="Train"
      write(SummaryUnit,"(a14,3(1x,f12.4),3(1x,f20.10))") adjustl(DumC),Cors%Cor,Cors%Cov/Cors%Var2,Cors%Cov/Cors%Var1,Cors%Cov,Cors%Var1,Cors%Var2

      deallocate(Ebv)

      G(:)=G(:)/ScaleCoef(:)*SdY

      i=maxval(nAnisTe(:))
      allocate(IdTe(i))
      allocate(GenosTe(i,nSnp))
      allocate(PhenoTe(i))
      allocate(Ebv(i))

      do Pop=1,nTePop

        open(newunit=GenoTeUnit,file=trim(GenoTeFile(Pop)),status="old")
        open(newunit=PhenoTeUnit,file=trim(PhenoTeFile(Pop)),status="old")

        do i=1,nAnisTe(Pop)
          read(GenoTeUnit,*) IdTe(i),GenosTe(i,:)
          read(PhenoTeUnit,*) DumC,PhenoTe(i)
          if (trim(DumC).ne.trim(IdTe(i))) then
            write(STDERR,"(a,i)") "ERROR: Individual identifications in the genotype and phenotype files do not match in prediction set ",Pop
            write(STDERR,"(a,i)") "ERROR: Line: ",i
            write(STDERR,"(a,a)") "ERROR: Genotype file identification: ",trim(DumC)
            write(STDERR,"(a,a)") "ERROR: Phenotype file identification: ",trim(IdTr(i))
            write(STDERR,"(a)") " "
            stop 1
          end if
        end do

        do j=1,nSnp
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
        call gemv(A=GenosTe(1:nAnisTe(Pop),:),x=G,y=Ebv(1:nAnisTe(Pop)))
        write(Unit,"(a20,2(1x,a20))") "Id","Observed","Estimate"
        do i=1,nAnisTe(Pop)
          write(Unit,"(a20,2(1x,f20.10))") IdTe(i),PhenoTe(i),Ebv(i)
        end do
        flush(Unit)
        close(Unit)

        Cors=cor(PhenoTe(1:nAnisTe(Pop)),Ebv(1:nAnisTe(Pop)))
        DumC="Predict"//Int2Char(Pop)
        write(SummaryUnit,"(a14,3(1x,f12.4),3(1x,f20.10))") adjustl(DumC),Cors%Cor,Cors%Cov/Cors%Var2,Cors%Cov/Cors%Var1,Cors%Cov,Cors%Var1,Cors%Var2

        ! open(newunit=GenoTeUnit,file="GenoTest"//Int2Char(Pop)//"Processed.txt",status="unknown")
        ! open(newunit=PhenoTeUnit,file="PhenoTest"//Int2Char(Pop)//"Processed.txt",status="unknown")
        ! do i=1,nAnisTe(Pop)
        !   write(GenoTeUnit,"("//Int2Char(nSnp)//"(1x,f20.16))") GenosTe(i,:)
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
      deallocate(Ebv)

      flush(SummaryUnit)
      close(SummaryUnit)
    end subroutine

    !###########################################################################
end module

!###############################################################################

program AlphaBayes

  use ISO_Fortran_Env, STDIN=>input_unit,STDOUT=>output_unit,STDERR=>error_unit
  use AlphaBayesModule
  use IntelRNGMod
  implicit none

  real(real32) :: Start,Finish

  call cpu_time(Start)

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

  call IntitialiseIntelRNG

  call ReadParam
  call ReadData
  if (trim(Method).eq."Ridge") then
    call RidgeRegression
  end if
  if (trim(Method).eq."RidgeMCMC") then
    call RidgeRegressionMCMC
  end if
  call Prediction

  call UnintitialiseIntelRNG

  call cpu_time(Finish)

  write(STDOUT,"(a,f20.4,a)") "Time duration of AlphaBayes: ",Finish-Start," seconds"
  write(STDOUT,"(a)") " "
end program

!###############################################################################
