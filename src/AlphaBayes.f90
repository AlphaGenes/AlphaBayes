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

module AlphaBayesMod

  use ISO_Fortran_Env, STDIN=>input_unit,STDOUT=>output_unit,STDERR=>error_unit
  use IntelRNGMod
  use AlphaHouseMod
  use AlphaStatMod

  implicit none

  integer(int32) :: nSnp,nSnpExternal,nAnisTr,nIter,nBurn,nProcessors,ScalingOpt,nTePop
  integer(int32),allocatable :: nAnisTe(:),FixedSnp(:),SnpPosition(:)

  real(real64) :: MeanY,VarY,SdY,VarG,VarE,Mu,nSnpR,nAnisTrR,ExpVarX,SumExpVarX,ObsVarX,SumObsVarX
  real(real64),allocatable :: SnpTmp(:),AlleleFreq(:),ScaleCoef(:),GenosTr(:,:),PhenoTr(:,:),G(:,:)

  character(len=1000) :: GenoTrFile,PhenoTrFile,FileFixedSnp,Method
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

      integer(int32) :: i,UnitSpec,UnitFixedSnp

      character(len=100) :: DumC

      open(newunit=UnitSpec,file=SPECFILE,status="old")

      read(UnitSpec,*) DumC,GenoTrFile
      read(UnitSpec,*) DumC,PhenoTrFile

      read(UnitSpec,*) DumC,nTePop
      allocate(PhenoTeFile(nTePop))
      allocate(GenoTeFile(nTePop))
      allocate(nAnisTe(nTePop))
      read(UnitSpec,*) DumC,GenoTeFile(:)
      read(UnitSpec,*) DumC,PhenoTeFile(:)

      read(UnitSpec,*) DumC,FileFixedSnp
      read(UnitSpec,*) DumC,nSnpExternal

      if (trim(FileFixedSnp)/="None") then
        open(newunit=UnitFixedSnp,file=trim(FileFixedSnp),status="old")
        allocate(FixedSnp(nSnpExternal))
        do i=1,nSnpExternal
          read(UnitFixedSnp,*) FixedSnp(i)
        end do
        nSnp=sum(FixedSnp(:))
        close(UnitFixedSnp)
      else
        nSnp=nSnpExternal
        allocate(FixedSnp(nSnpExternal))
        FixedSnp(:)=1
      end if
      nSnpR=dble(nSnp)

      read(UnitSpec,*) DumC,nAnisTr
      nAnisTrR=dble(nAnisTr)
      read(UnitSpec,*) DumC,nAnisTe(:)

      read(UnitSpec,*) DumC,nIter
      read(UnitSpec,*) DumC,nBurn

      read(UnitSpec,*) DumC,VarG ! here it is meant as genetic variance!!!
      read(UnitSpec,*) DumC,VarE

      read(UnitSpec,*) DumC,nProcessors
      call OMP_SET_NUM_THREADS(nProcessors)

      read(UnitSpec,*) DumC,ScalingOpt
      if (ScalingOpt<1 .or. ScalingOpt>5) then
        write(STDERR,"(a)") "ERROR: ScalingOption must be between 1 and 5"
        write(STDERR,"(a)") " "
        stop 1
      endif

      read(UnitSpec,*) DumC,Method
      close(UnitSpec)
    end subroutine

    !###########################################################################

    subroutine ReadData
      implicit none

      integer(int32) :: i,j,nNotMissing,UnitGenoTr,UnitPhenoTr!,UnitAlleleFreq

      character(len=100) :: DumC

      real(real64) :: TmpStat

      open(newunit=UnitGenoTr,file=trim(GenoTrFile),status="old")
      open(newunit=UnitPhenoTr,file=trim(PhenoTrFile),status="old")
      !open(newunit=UnitAlleleFreq,file="AlleleFreq.txt",status="unknown")
      !write(UnitAlleleFreq,"(a11,1x,a11)") "Snp","AlleleFreq"

      allocate(SnpTmp(nSnpExternal))
      allocate(GenosTr(nAnisTr,nSnp))
      allocate(PhenoTr(nAnisTr,1))
      allocate(G(nSnp,1))
      allocate(AlleleFreq(nSnp))
      allocate(ScaleCoef(nSnp))
      allocate(SnpPosition(nSnp))
      allocate(IdTr(nAnisTr))

      j=0
      do i=1,nSnpExternal
        if (FixedSnp(i)==1) then
          j=j+1
          SnpPosition(j)=i
        end if
      end do

      do i=1,nAnisTr
        read(UnitGenoTr,*) DumC,SnpTmp(:)
        GenosTr(i,:)=SnpTmp(SnpPosition(:))
        read(UnitPhenoTr,*) IdTr(i),PhenoTr(i,1)
        if (trim(DumC) /= trim(IdTr(i))) then
          write(STDERR,"(a)") "ERROR: Individual identifications in the genotype and phenotype files do not match in the training set"
          write(STDERR,"(a,i)") "ERROR: Line: ",i
          write(STDERR,"(a,a)") "ERROR: Genotype file identification: ",trim(DumC)
          write(STDERR,"(a,a)") "ERROR: Phenotype file identification: ",trim(IdTr(i))
          write(STDERR,"(a)") " "
          stop 1
        end if
      end do

      MeanY=Mean(PhenoTr(:,1))
      VarY=Var(PhenoTr(:,1),MeanY)
      SdY=sqrt(VarY)
      PhenoTr(:,1)=(PhenoTr(:,1)-MeanY)/SdY

      SumExpVarX=0.0d0
      SumObsVarX=0.0d0
      do j=1,nSnp

        ! Compute allele freqs
        AlleleFreq(j)=0.0d0
        nNotMissing=0
        do i=1,nAnisTr
          if ((GenosTr(i,j)>-0.1).and.(GenosTr(i,j)<2.1)) then
            AlleleFreq(j)=AlleleFreq(j)+GenosTr(i,j)
            nNotMissing=nNotMissing+1
          end if
        end do
        if (nNotMissing/=0) then
          AlleleFreq(j)=AlleleFreq(j)/dble(2*nNotMissing)
        else
          AlleleFreq(j)=0.0d0
        end if
        !write(UnitAlleleFreq,"(i11,1x,f11.8)") j,AlleleFreq(j)

        ! Fix any odd data
        do i=1,nAnisTr
          if ((GenosTr(i,j)<-0.1).or.(GenosTr(i,j)>2.1)) then
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

        if (ScalingOpt==2) then
          ! Scale by marker specific variance - expected
          ExpVarX=sqrt(ExpVarX)
          ScaleCoef(j)=ExpVarX
          GenosTr(:,j)=GenosTr(:,j)/ScaleCoef(j)
        end if

        if (ScalingOpt==3) then
          ! Scale by marker specific variance - observed
          ObsVarX=sqrt(ObsVarX)
          ScaleCoef(j)=ObsVarX
          GenosTr(:,j)=GenosTr(:,j)/ScaleCoef(j)
        end if

      end do

      if (ScalingOpt==4) then
        ! Scale by average marker variance - expected
        ExpVarX=sqrt(SumExpVarX/nSnpR)
        ScaleCoef(:)=ExpVarX
        GenosTr(:,:)=GenosTr(:,:)/ScaleCoef(1)
      end if

      if (ScalingOpt==5) then
        ! Scale by average marker variance - observed
        ObsVarX=sqrt(SumObsVarX/nSnpR)
        ScaleCoef(:)=ObsVarX
        GenosTr(:,:)=GenosTr(:,:)/ScaleCoef(1)
      end if

      close(UnitGenoTr)
      close(UnitPhenoTr)
      !close(UnitAlleleFreq)

      ! open(newunit=UnitGenoTr,file="GenoTrainProcessed.txt",status="unknown")
      ! open(newunit=UnitPhenoTr,file="PhenoTrainProcessed.txt",status="unknown")
      ! do i=1,nAnisTr
      !   write(UnitGenoTr,"("//Int2Char(nSnp)//"(1x,f20.16))") GenosTr(i,:)
      !   write(UnitPhenoTr,*) PhenoTr(i,1)
      ! end do
      ! close(UnitGenoTr)
      ! close(UnitPhenoTr)
    end subroutine

    !###########################################################################

    subroutine RidgeRegression
      implicit none

      integer(int32) :: Iter,j,Snp,RandomOrdering(nSnp),UnitVar

      real(real64) :: ddot,VarS,Rhs,Lhs,Sol,Diff,Eps
      real(real64),allocatable :: XpXTauE(:,:),Xg(:,:),E(:,:)

      allocate(XpXTauE(nSnp,1))
      allocate(Xg(nAnisTr,1))
      allocate(E(nAnisTr,1))

      ! Working phenotypes
      E(:,1)=PhenoTr(:,1)

      ! Initialise
      Mu=0.0d0
      G=0.0d0

      ! Construct XpXTauE
      do j=1,nSnp
        XpXTauE(j,1)=ddot(nAnisTr,GenosTr(:,j),1,GenosTr(:,j),1) + tiny(GenosTr(1,1))
      end do
      XpXTauE(:,1)=XpXTauE(:,1)/VarE ! can do it only once for all rounds!!!
      VarS=VarG/nSnpR ! approximate variance of allele substitution effects

      ! Iterate
      do Iter=1,nIter
        Eps=0.0d0

        ! Intercept
        Lhs=nAnisTrR/VarE
        Rhs=sum(E(:,1))/VarE + nAnisTrR*Mu/VarE
        Sol=Rhs/Lhs
        Diff=Sol-Mu
        E(:,1)=E(:,1)-Diff
        Mu=Sol
        Eps=Eps + Diff*Diff

        ! Snp effects
        RandomOrdering=RandomOrder(nSnp)
        do j=1,nSnp
          Snp=RandomOrdering(j)
          Lhs=XpXTauE(Snp,1) + 1.0d0/VarS
          Rhs=ddot(nAnisTr,GenosTr(:,Snp),1,E(:,1),1)/VarE + XpXTauE(Snp,1)*G(Snp,1)
          Sol=Rhs/Lhs
          Diff=Sol-G(Snp,1)
          E(:,1)=E(:,1) - GenosTr(:,Snp)*Diff
          G(Snp,1)=Sol
          Eps=Eps + Diff*Diff
        end do

        ! Recompute residuals to avoid rounding errors
        if (mod(Iter,200)==0) then
          call dgemm("n","n",nAnisTr,1,nSnp,ONER,GenosTr,nAnisTr,G,nSnp,ZEROR,Xg,nAnisTr)
          E(:,1)=PhenoTr(:,1)-Mu-Xg(:,1)
        end if

        ! Stopping criteria
        if (eps < 1.0E-8) then
          exit
        end if
      end do

      open(newunit=UnitVar,file=INTERCEPTESTIMATEFILE,status="unknown")
      write(UnitVar,"(f20.10)") Mu*SdY
      close(UnitVar)
      open(newunit=UnitVar,file=SNPESTIMATEFILE,status="unknown")
      do j=1,nSnp
        write(UnitVar,"(f20.10)") G(j,1)*SdY
      end do
      close(UnitVar)

      deallocate(XpXTauE)
      deallocate(Xg)
      deallocate(E)
    end subroutine

    !###########################################################################

    subroutine RidgeRegressionMCMC
      implicit none

      integer(int32) :: Iter,i,j,Snp,RandomOrdering(nSnp),UnitVar

      real(real64) :: TmpR,ddot,Rhs,Lhs,Sol,Diff,nSampR
      real(real64) :: MuSamp,VarESamp,VarGSamp
      real(real64) :: R2,EDF0,EDF,GDF0,GDF,ES0,GS0,MSX,EpE,GpG
      real(real64),allocatable :: GSamp(:,:),MuAll(:),GAll(:,:),VarEAll(:),VarGAll(:)
      real(real64),allocatable :: SX2(:),MX2(:),XpX(:,:),Xg(:,:),E(:,:)
      real(real64),allocatable :: GaussDevMu(:),GaussDevSnp(:),GammaDevE(:),GammaDevG(:)

      allocate(XpX(nSnp,1))
      allocate(Xg(nAnisTr,1))
      allocate(E(nAnisTr,1))
      allocate(GSamp(nSnp,1))
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
      E(:,1)=PhenoTr(:,1)

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

      TmpR=Var(E(:,1)) ! should be 1 if Phen is standardized
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
        XpX(j,1)=ddot(nAnisTr,GenosTr(:,j),1,GenosTr(:,j),1) + tiny(GenosTr(1,1))
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
        Rhs=sum(E(:,1))/VarESamp + nAnisTrR*MuSamp/VarESamp
        Sol=Rhs/Lhs + GaussDevMu(Iter)/sqrt(Lhs)
        Diff=Sol-MuSamp
        E(:,1)=E(:,1)-Diff
        MuSamp=Sol

        ! Snp effects
        RandomOrdering=RandomOrder(nSnp)
        GaussDevSnp(:)=SampleIntelGaussD(n=nSnp)
        do j=1,nSnp
          Snp=RandomOrdering(j)
          ! This does not work (got too big variances!), and I can not see why
          !Lhs=XpX(Snp,1) + VarESamp/VarGSamp
          !Rhs=ddot(nAnisTr,GenosTr(:,Snp),1,E(:,1),1) + XpX(Snp,1)*GSamp(Snp,1)
          Lhs=XpX(Snp,1)/VarESamp + 1.0d0/VarGSamp
          Rhs=ddot(nAnisTr,GenosTr(:,Snp),1,E(:,1),1)/VarESamp + XpX(Snp,1)*GSamp(Snp,1)/VarESamp
          Sol=Rhs/Lhs + GaussDevSnp(j)/sqrt(Lhs)
          Diff=Sol-GSamp(Snp,1)
          E(:,1)=E(:,1) - GenosTr(:,Snp)*Diff
          GSamp(Snp,1)=Sol
        end do

        ! Snp variance
        GpG=ddot(nSnp,GSamp(:,1),1,GSamp(:,1),1)+GS0
        VarGSamp=GpG/GammaDevG(Iter)

        ! Recompute residuals to avoid rounding errors
        if (mod(Iter,200)==0) then
          call dgemm("n","n",nAnisTr,1,nSnp,ONER,GenosTr,nAnisTr,GSamp,nSnp,ZEROR,Xg,nAnisTr)
          E(:,1)=PhenoTr(:,1)-MuSamp-Xg(:,1)
        end if

        ! Residual variance
        EpE=ddot(nAnisTr,E(:,1),1,E(:,1),1)+ES0
        VarESamp=EpE/GammaDevE(Iter)

        MuAll(Iter)=MuSamp
        GAll(:,Iter)=GSamp(:,1)
        VarEAll(Iter)=VarESamp
        VarGAll(Iter)=VarGSamp

      end do

      ! Posterior means
      do j=nBurn+1,nIter
        Mu=Mu+MuAll(j)/nSampR
        G(:,1)=G(:,1)+GAll(:,j)/nSampR
        VarE=VarE+VarEAll(j)/nSampR
        VarG=VarG+VarGAll(j)/nSampR
      end do

      ! Outputs
      open(newunit=UnitVar,file=INTERCEPTSAMPLESFILE,status="unknown")
      write(UnitVar,*) MuAll(:)*SdY
      close(UnitVar)

      open(newunit=UnitVar,file=SNPSAMPLESFILE,status="unknown")
      GaussDevSnp(:)=0.0d0 ! just a placeholder for zero output of fixed markers
      j=0
      do i=1,nSnpExternal
          if (FixedSnp(i)==1) then
            j=j+1
            write(UnitVar,*) GAll(j,:)/ScaleCoef(j)*SdY
          else
            write(UnitVar,*) GaussDevSnp(:)
          end if
      end do
      close(UnitVar)

      open(newunit=UnitVar,file=RESIDUALVARIANCEESTIMATEFILE,status="unknown")
      write(UnitVar,"(f20.10)") VarE*VarY
      close(UnitVar)
      open(newunit=UnitVar,file=RESIDUALVARIANCESAMPLESFILE,status="unknown")
      write(UnitVar,*) VarEAll(:)*VarY
      close(UnitVar)

      open(newunit=UnitVar,file=SNPVARIANCEESTIMATEFILE,status="unknown")
      write(UnitVar,"(f20.10)") VarG*VarY
      close(UnitVar)
      open(newunit=UnitVar,file=SNPVARIANCESAMPLESFILE,status="unknown")
      write(UnitVar,*) VarGAll(:)*VarY
      close(UnitVar)

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

      integer(int32) :: i,j,Pop,UnitSum,UnitSnpSol,UnitGenoTe,UnitPhenoTe,UnitEbv

      real(real64),allocatable :: Ebv(:,:),PhenoTe(:,:),GenosTe(:,:)

      character(len=100) :: DumC,File
      character(len=20),allocatable :: IdTe(:)

      type(CorrelationReal64) :: Cors

      open(newunit=UnitSum,file="PredictionsSummary.txt",status="unknown")
      write(UnitSum,"(a14,3(1x,a12),3(1x,a20))") "Set","CorsObsEst","SlopeObsEst","SlopeEstObs","CovObsEst","VarObs","VarEst"

      allocate(Ebv(nAnisTr,1))

      open(newunit=UnitEbv,file="PredictionsForSetTrain.txt",status="unknown")
      call dgemm("n","n",nAnisTr,1,nSnp,ONER,GenosTr,nAnisTr,G,nSnp,ZEROR,Ebv,nAnisTr)
      PhenoTr(:,1)=MeanY+PhenoTr(:,1)*SdY
      write(UnitEbv,"(a20,2(1x,a20))") "Id","Observed","Estimate"
      do i=1,nAnisTr
        write(UnitEbv,"(a20,2(1x,f20.10))") IdTr(i),PhenoTr(i,1),Ebv(i,1)
      end do
      flush(UnitEbv)
      close(UnitEbv)

      Cors=Cor(PhenoTr(:,1),Ebv(:,1))
      DumC="Train"
      write(UnitSum,"(a14,3(1x,f12.4),3(1x,f20.10))") adjustl(DumC),Cors%Cor,Cors%Cov/Cors%Var2,Cors%Cov/Cors%Var1,Cors%Cov,Cors%Var1,Cors%Var2

      deallocate(Ebv)

      ! Rescale marker estimates back to observed phenotype and genotype scale

      ! y     = mu     + b*x     + e
      ! y     = mu     + b*x/SDX + e
      ! y     = mu     + b*z     + e
      ! y/SDY = mu/SDY + b*z/SDY + e/SDY
      ! y'    = mu'    + b'*z    + e'
      ! b'*z=b*z/SDY
      ! b*z =b'*z*SDY
      ! b*x/SDX=b'*x/SDX*SDY
      ! b=b'/SDX*SDY

      G(:,1)=G(:,1)/ScaleCoef(:)*SdY

      ! Output
      open(newunit=UnitSnpSol,file=INTERCEPTESTIMATEFILE,status="unknown")
      write(UnitSnpSol,"(f20.10)") Mu
      flush(UnitSnpSol)
      close(UnitSnpSol)

      open(newunit=UnitSnpSol,file=SNPESTIMATEFILE,status="unknown")
      write(UnitSnpSol,"(a20,1x,a20)") "Snp","Estimate"
      j=0
      do i=1,nSnpExternal
          if (FixedSnp(i)==1) then
            j=j+1
            write(UnitSnpSol,"(i20,1x,f20.10)") i,G(j,1)
          else
            write(UnitSnpSol,"(i20,1x,f20.10)") i,0.0d0
          end if
      end do
      flush(UnitSnpSol)
      close(UnitSnpSol)

      i=maxval(nAnisTe(:))
      allocate(IdTe(i))
      allocate(GenosTe(i,nSnp))
      allocate(PhenoTe(i,1))
      allocate(Ebv(i,1))

      do Pop=1,nTePop

        open(newunit=UnitGenoTe,file=trim(GenoTeFile(Pop)),status="old")
        open(newunit=UnitPhenoTe,file=trim(PhenoTeFile(Pop)),status="old")

        do i=1,nAnisTe(Pop)
          read(UnitGenoTe,*) IdTe(i),SnpTmp(:)
          GenosTe(i,:)=SnpTmp(SnpPosition(:))
          read(UnitPhenoTe,*) DumC,PhenoTe(i,1)
          if (trim(DumC) /= trim(IdTe(i))) then
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
            if ((GenosTe(i,j)<-0.1).or.(GenosTe(i,j)>2.1)) then
              GenosTe(i,j)=2.0d0*AlleleFreq(j)
            end if
          end do

          ! Standardize

          ! ... center
          GenosTe(:,j)=GenosTe(:,j)-(2.0d0*AlleleFreq(j))

          ! ... scale
          ! no need because we scaled the marker solutions
        end do

        close(UnitGenoTe)
        close(UnitPhenoTe)

        File="PredictionsForSetPredict"//Int2Char(Pop)//".txt"
        open(newunit=UnitEbv,file=trim(File),status="unknown")
        call dgemm("n","n",nAnisTe(Pop),1,nSnp,ONER,GenosTe(1:nAnisTe(Pop),:),nAnisTe(Pop),G,nSnp,ZEROR,Ebv(1:nAnisTe(Pop),1),nAnisTe(Pop))
        write(UnitEbv,"(a20,2(1x,a20))") "Id","Observed","Estimate"
        do i=1,nAnisTe(Pop)
          write(UnitEbv,"(a20,2(1x,f20.10))") IdTe(i),PhenoTe(i,1),Ebv(i,1)
        end do
        flush(UnitEbv)
        close(UnitEbv)

        Cors=cor(PhenoTe(1:nAnisTe(Pop),1),Ebv(1:nAnisTe(Pop),1))
        DumC="Predict"//Int2Char(Pop)
        write(UnitSum,"(a14,3(1x,f12.4),3(1x,f20.10))") adjustl(DumC),Cors%Cor,Cors%Cov/Cors%Var2,Cors%Cov/Cors%Var1,Cors%Cov,Cors%Var1,Cors%Var2

        ! open(newunit=UnitGenoTe,file="GenoTest"//Int2Char(Pop)//"Processed.txt",status="unknown")
        ! open(newunit=UnitPhenoTe,file="PhenoTest"//Int2Char(Pop)//"Processed.txt",status="unknown")
        ! do i=1,nAnisTe(Pop)
        !   write(UnitGenoTe,"("//Int2Char(nSnp)//"(1x,f20.16))") GenosTe(i,:)
        !   write(UnitPhenoTe,*) PhenoTe(i,1)
        ! end do
        ! close(UnitGenoTe)
        ! close(UnitPhenoTe)
      end do

      deallocate(IdTe)
      deallocate(GenosTe)
      deallocate(PhenoTe)
      deallocate(Ebv)

      flush(UnitSum)
      close(UnitSum)
    end subroutine

    !###########################################################################
end module

!###############################################################################

program AlphaBayes

  use ISO_Fortran_Env, STDIN=>input_unit,STDOUT=>output_unit,STDERR=>error_unit
  use AlphaBayesMod
  use IntelRNGMod
  implicit none

  real(real32) :: Start,Finish

  include "mkl_vml.f90"

  call cpu_time(Start)

  write(STDOUT,"(a)") ""
  write(STDOUT,"(a30,a,a30)") " ","**********************"," "
  write(STDOUT,"(a30,a,a30)") " ","*                    *"," "
  write(STDOUT,"(a30,a,a30)") " ","*     AlphaBayes     *"," "
  write(STDOUT,"(a30,a,a30)") " ","*                    *"," "
  write(STDOUT,"(a30,a,a30)") " ","**********************"
  write(STDOUT,"(a30,a,a30)") " ","VERSION:"//TOSTRING(VERS)," "
  write(STDOUT,"(a)") ""
  write(STDOUT,"(a35,a)")     " ","No Liability"
  write(STDOUT,"(a25,a)")     " ","Bugs to John.Hickey@roslin.ed.ac.uk"
  write(STDOUT,"(a)") ""

  call IntitialiseIntelRNG

  call ReadParam
  call ReadData
  if (trim(Method)=="Ridge") then
    call RidgeRegression
  end if
  if (trim(Method)=="RidgeMCMC") then
    call RidgeRegressionMCMC
  end if
  call Prediction

  call UnintitialiseIntelRNG

  call cpu_time(Finish)

  write(STDOUT,"(a,f20.4,a)") "Time duration of AlphaBayes: ",Finish-Start," seconds"
  write(STDOUT,"(a)") " "
end program

!###############################################################################
