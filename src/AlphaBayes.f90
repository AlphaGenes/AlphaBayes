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
  use AlphaHouseMod, only : RandomOrder,Int2Char
  use AlphaStatMod, only : CalcMean,CalcVar,CalcCorrelation,CorrelationD

  implicit none

  integer(int32) :: nSnp,nSnpExternal,nAnisTr,nIter,nBurn,nProcessors,ScalingOpt,nTePop
  integer(int32),allocatable :: nAnisTe(:),FixedSnp(:),SnpPosition(:)

  real(real64) :: MeanY,VarY,VarG,VarE,Mu,nAnisTrR,ExpVarX,SumExpVarX,ObsVarX,SumObsVarX
  real(real64),allocatable :: SnpTmp(:),AlleleFreq(:),ScaleCoef(:),GenosTr(:,:),PhenoTr(:,:),G(:,:)

  character(len=1000) :: GenoTrFile,PhenoTrFile,FileFixedSnp,Method
  character(len=1000),allocatable :: GenoTeFile(:),PhenoTeFile(:)
  character(len=20),allocatable :: IdTr(:)

  REAL(REAL64),PARAMETER :: ONER=1.0d0,ZEROR=0.0d0

  private
  public :: ReadParam,ReadData,Method,Prediction
  public :: RidgeRegression,RidgeRegressionMCMC!,BayesA

  contains

    !###########################################################################

    subroutine ReadParam
      implicit none

      integer(int32) :: i,UnitSpec,UnitFixedSnp

      character(len=100) :: DumC

      open(newunit=UnitSpec,file="AlphaBayesSpec.txt",status="old")

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

      read(UnitSpec,*) DumC,nAnisTr
      nAnisTrR=dble(nAnisTr)
      read(UnitSpec,*) DumC,nAnisTe(:)

      read(UnitSpec,*) DumC,nIter
      read(UnitSpec,*) DumC,nBurn

      read(UnitSpec,*) DumC,VarG
      read(UnitSpec,*) DumC,VarE

      read(UnitSpec,*) DumC,nProcessors
      call OMP_SET_NUM_THREADS(nProcessors)

      read(UnitSpec,*) DumC,ScalingOpt

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

      call CalcMean(PhenoTr(:,1),MeanY)
      call CalcVar(PhenoTr(:,1),VarY,MeanY)
      PhenoTr(:,1)=(PhenoTr(:,1)-MeanY)/sqrt(VarY)

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
        call CalcVar(GenosTr(:,j),TmpStat)
        ExpVarX=2.0d0*(1.0d0-AlleleFreq(j))*AlleleFreq(j)+0.00001d0 ! if p=0.00001, then 2*p*q=0.00001
        SumExpVarX=SumExpVarX+ExpVarX
        ObsVarX=TmpStat+0.00001d0
        SumObsVarX=SumObsVarX+ObsVarX

        ! ... center
        GenosTr(:,j)=GenosTr(:,j)-(2.0d0*AlleleFreq(j))

        ! ... scale
        if (ScalingOpt==1) then
          ScaleCoef(j)=1.0d0
        end if

        if (ScalingOpt==2) then
          ! Scale by marker specific variance - expected
          ExpVarX=sqrt(ExpVarX)
          GenosTr(:,j)=GenosTr(:,j)/ExpVarX
          ScaleCoef(j)=ExpVarX
        end if

        if (ScalingOpt==3) then
          ! Scale by marker specific variance - observed
          ObsVarX=sqrt(ObsVarX)
          GenosTr(:,j)=GenosTr(:,j)/ObsVarX
          ScaleCoef(j)=ObsVarX
        end if

      end do

      if (ScalingOpt==4) then
        ! Scale by average marker variance - expected
        ExpVarX=sqrt(SumExpVarX/dble(nSnp))
        GenosTr(:,:)=GenosTr(:,:)/ExpVarX
        ScaleCoef(:)=ExpVarX
      end if

      if (ScalingOpt==5) then
        ! Scale by average marker variance - observed
        ObsVarX=sqrt(SumObsVarX/dble(nSnp))
        GenosTr(:,:)=GenosTr(:,:)/ObsVarX
        ScaleCoef(:)=ObsVarX
      end if

      close(UnitGenoTr)
      close(UnitPhenoTr)
      !close(UnitAlleleFreq)
    end subroutine

    !###########################################################################

    subroutine RidgeRegression
      implicit none

      integer(int32) :: Iter,j,SnpId,RandomOrdering(nSnp)

      real(real64) :: ddot,Eps,Rhs,Lhs,SolOld,LambdaSnp
      real(real64),allocatable :: XpX(:,:),Xg(:,:),E(:,:)

      allocate(XpX(nSnp,1))
      allocate(Xg(nAnisTr,1))
      allocate(E(nAnisTr,1))

      Mu=0.0d0
      G=0.000001d0

      ! Construct XpX and Lambda
      do j=1,nSnp
        XpX(j,1)=ddot(nAnisTr,GenosTr(:,j),1,GenosTr(:,j),1) + 0.000000000000001d0
      end do

      LambdaSnp=VarE/(VarG/dble(nSnp))

      ! Working phenotypes
      E(:,1)=PhenoTr(:,1)

      do Iter=1,nIter
        ! Intercept
        E(:,1)=E(:,1)+Mu
        Mu=sum(E(:,1))/nAnisTrR
        E(:,1)=E(:,1)-Mu

        ! Snp effects
        eps=0.0d0
        call RandomOrder(RandomOrdering,nSnp)
        do j=1,nSnp
          SnpId=RandomOrdering(j)
          E(:,1)=E(:,1)+(GenosTr(:,SnpId)*G(SnpId,1))
          Lhs=XpX(SnpId,1)+LambdaSnp
          Rhs=ddot(nAnisTr,GenosTr(:,SnpId),1,E(:,1),1)
          SolOld=G(SnpId,1)
          G(SnpId,1)=Rhs/Lhs
          E(:,1)=E(:,1)-(GenosTr(:,SnpId)*G(SnpId,1))
          Eps=Eps+(G(SnpId,1)-SolOld)**2
        end do

        if (mod(Iter,200)==0) then
          call dgemm("n","n",nAnisTr,1,nSnp,ONER,GenosTr,nAnisTr,G,nSnp,ZEROR,Xg,nAnisTr)
          E(:,1)=PhenoTr(:,1)-Xg(:,1)-Mu
        end if

        if (eps < 1e-8) then
          exit
        end if
      end do

      deallocate(XpX)
      deallocate(Xg)
      deallocate(E)
    end subroutine

    !###########################################################################

    subroutine RidgeRegressionMCMC
      implicit none

      integer(int32) :: i,Iter,j,SnpId,RandomOrdering(nSnp),UnitVar

      real(real64) :: ddot,Rhs,Lhs,LambdaSnp,nSampR,VarESamp,VarGSamp,VarGAccum,VarEAccum
      real(real64) :: R2,EDF0,EDF2,GDF0,GDF2,ES0,GS0,MSX,EpE,GpG
      real(real64),allocatable :: GAccum(:,:),SX2(:),MX2(:),XpX(:,:),Xg(:,:),E(:,:)
      real(real64),allocatable :: GaussDev(:),GammaDevE(:),GammaDevG(:)

      allocate(XpX(nSnp,1))
      allocate(Xg(nAnisTr,1))
      allocate(E(nAnisTr,1))
      allocate(GAccum(nSnp,1))
      allocate(SX2(nSnp))
      allocate(MX2(nSnp))
      allocate(GaussDev(nSnp))
      allocate(GammaDevE(nIter))
      allocate(GammaDevG(nIter))

      Mu=0.0d0
      G=0.1d0
      GAccum=0.0d0
      VarGAccum=0.0d0
      VarEAccum=0.0d0

      ! These prior parameters are modelled as in BGLR
      R2=0.5d0

      EDF0=5.0d0
      EDF2=(dble(nAnisTr)+EDF0)/2.0d0
      ES0=(VarY*(1.0d0-R2))*(EDF0+2.0d0)

      GDF0=5.0d0
      GDF2=(dble(nSnp)+GDF0)/2.0d0
      SX2(:)=0.0d0
      MX2(:)=0.0d0
      do i=1,nSnp
        SX2(i)=sum(GenosTr(:,i)*GenosTr(:,i))
        MX2(i)=sum(GenosTr(:,i))/nAnisTrR
        MX2(i)=MX2(i)*MX2(i)
      end do
      MSX=sum(SX2)/nAnisTrR-sum(MX2)
      GS0=((VarY*R2)/MSX)*(GDF0+2.0d0)

      deallocate(SX2)
      deallocate(MX2)

      ! Construct XpX
      do j=1,nSnp
        XpX(j,1)=ddot(nAnisTr,GenosTr(:,j),1,GenosTr(:,j),1) + 0.000000000000001d0
      end do

      ! Working phenotypes
      E(:,1)=PhenoTr(:,1)

      ! Gamma deviates
      GammaDevE(:)=SampleIntelGammaD(n=nIter,alpha=EDF2,beta=2.0d0)
      GammaDevG(:)=SampleIntelGammaD(n=nIter,alpha=GDF2,beta=2.0d0)

      nSampR=dble(nIter-nBurn)
      do Iter=1,nIter
        ! Residual variance
        EpE=ddot(nAnisTr,E(:,1),1,E(:,1),1) + ES0 + 0.000000000000001d0
        VarESamp=EpE/GammaDevE(Iter)

        ! Snp variance
        GpG=ddot(nSnp,G(:,1),1,G(:,1),1) + GS0 + 0.000000000000001d0
        VarGSamp=GpG/GammaDevG(Iter)

        ! Ratio of variances
        LambdaSnp=VarESamp/VarGSamp

        ! Intercept
        E(:,1)=E(:,1)+Mu
        Mu=sum(E(:,1))/nAnisTrR
        E(:,1)=E(:,1)-Mu

        ! Snp effects
        call RandomOrder(RandomOrdering,nSnp)
        GaussDev(:)=SampleIntelGaussD(n=nSnp)
        do j=1,nSnp
          SnpId=RandomOrdering(j)
          E(:,1)=E(:,1)+(GenosTr(:,SnpId)*G(SnpId,1))
          Lhs=XpX(SnpId,1)+LambdaSnp
          Rhs=ddot(nAnisTr,GenosTr(:,SnpId),1,E(:,1),1)
          G(SnpId,1)=(Rhs/Lhs)+(GaussDev(SnpId)/sqrt(Lhs))
          E(:,1)=E(:,1)-(GenosTr(:,SnpId)*G(SnpId,1))
        end do
        if (Iter>nBurn) then
          GAccum=GAccum+G/nSampR
          VarGAccum=VarGAccum+VarGSamp/nSampR
          VarEAccum=VarEAccum+VarESamp/nSampR
        end if

        ! Recompute residuals to avoid rounding errors
        if (mod(Iter,200)==0) then
          call dgemm("n","n",nAnisTr,1,nSnp,ONER,GenosTr,nAnisTr,G,nSnp,ZEROR,Xg,nAnisTr)
          E(:,1)=PhenoTr(:,1)-Xg(:,1)-Mu
        end if

      end do

      G=GAccum
      open(newunit=UnitVar,file="ResidualVarianceEstimate.txt",status="unknown")
      write(UnitVar,"(f20.10)") VarEAccum
      close(UnitVar)
      open(newunit=UnitVar,file="SnpVarianceEstimate.txt",status="unknown")
      write(UnitVar,"(f20.10)") VarGAccum
      close(UnitVar)

      deallocate(XpX)
      deallocate(Xg)
      deallocate(E)
      deallocate(GAccum)
    end subroutine

    !###########################################################################

    subroutine BayesA
      implicit none

      integer(int32) :: Iter,j,SnpId,RandomOrdering(nSnp)

      real(real64) :: ddot,Rhs,Lhs,LambdaSnp,EpE,nSampR,VarESamp
      real(real64),allocatable :: GAccum(:,:),VarGSamp(:,:),XpX(:,:),Xg(:,:),E(:,:)
      real(real64) :: vdf,ShapeCoefficient,ScaleCoefficient,sc,gampe1,gampe2
      real(real64),allocatable :: GaussDev(:),GammaDevE(:),GammaDevG(:)

      allocate(XpX(nSnp,1))
      allocate(Xg(nAnisTr,1))
      allocate(E(nAnisTr,1))
      allocate(VarGSamp(nSnp,1))
      allocate(GAccum(nSnp,1))
      allocate(GaussDev(nSnp))
      allocate(GammaDevE(nIter))
      allocate(GammaDevG(nSnp))

      Mu=0.0d0
      G=0.1d0
      GAccum=0.0d0

      ! TODO: model the prior parameters as in BGLR
      vdf=4.012d0
      ShapeCoefficient=(vdf/2.0d0) !Ben
      ScaleCoefficient=((vdf-2.0d0)*VarG)/SumExpVarX
      sc=2.0d0
      gampe1=(dble(nAnisTr)/2.0d0)-2.0d0
      gampe2=2.0d0

      ! Construct XpX
      do j=1,nSnp
          XpX(j,1)=ddot(nAnisTr,GenosTr(:,j),1,GenosTr(:,j),1) + 0.000000000000001d0
      end do

      ! Working phenotypes
      E(:,1)=PhenoTr(:,1)

      nSampR=dble(nIter-nBurn)
      GammaDevE(:)=SampleIntelGammaD(n=nIter,alpha=gampe1,beta=gampe2)
      do Iter=1,nIter
        ! Residual variance
        EpE=ddot(nAnisTr,E(:,1),1,E(:,1),1) + 0.000000000000001d0
        VarESamp=EpE/GammaDevE(Iter)

        ! Snp variances
        GammaDevG(:)=SampleIntelGammaD(n=nSnp,alpha=ShapeCoefficient,beta=sc)
        VarGSamp(:,1)=(ScaleCoefficient+(G(:,1)*G(:,1)))/GammaDevG(:)

        !Intercept
        E(:,1)=E(:,1)+Mu
        Mu=sum(E(:,1))/nAnisTrR
        E(:,1)=E(:,1)-Mu

        !Snp effects
        call RandomOrder(RandomOrdering,nSnp)
        GaussDev(:)=SampleIntelGaussD(n=nSnp)
        do j=1,nSnp
          SnpId=RandomOrdering(j)
          E(:,1)=E(:,1)+(GenosTr(:,SnpId)*G(SnpId,1))
          LambdaSnp=VarESamp/VarGSamp(SnpId,1)
          Lhs=XpX(SnpId,1)+LambdaSnp
          Rhs=ddot(nAnisTr,GenosTr(:,SnpId),1,E(:,1),1)
          G(SnpId,1)=(Rhs/Lhs)+(GaussDev(SnpId)/sqrt(Lhs))
          E(:,1)=E(:,1)-(GenosTr(:,SnpId)*G(SnpId,1))
        end do
        if (Iter>nBurn) then
          GAccum=GAccum+G/nSampR
        end if

        ! Recompute residuals to avoid rounding errors
        if (mod(Iter,200)==0) then
          call dgemm('n','n',nAnisTr,1,nSnp,ONER,GenosTr,nAnisTr,G,nSnp,ZEROR,Xg,nAnisTr)
          E(:,1)=PhenoTr(:,1)-Xg(:,1)-Mu
        end if

      end do

      G=GAccum

      deallocate(XpX)
      deallocate(Xg)
      deallocate(E)
      deallocate(VarGSamp)
      deallocate(GAccum)
    end subroutine

    !###########################################################################

    subroutine Prediction
      implicit none

      integer(int32) :: i,j,Pop,UnitSum,UnitSnpSol,UnitGenoTe,UnitPhenoTe,UnitEbv

      real(real64),allocatable :: Ebv(:,:),PhenoTe(:,:),GenosTe(:,:)

      character(len=100) :: DumC,File
      character(len=20),allocatable :: IdTe(:)

      type(CorrelationD) :: Cor

      open(newunit=UnitSum,file="PredictionsSummary.txt",status="unknown")
      write(UnitSum,"(a14,3(1x,a12),3(1x,a20))") "Set","CorObsEst","SlopeObsEst","SlopeEstObs","CovObsEst","VarObs","VarEst"

      allocate(Ebv(nAnisTr,1))

      open(newunit=UnitEbv,file="PredictionsForSetTrain.txt",status="unknown")
      call dgemm("n","n",nAnisTr,1,nSnp,ONER,GenosTr,nAnisTr,G,nSnp,ZEROR,Ebv,nAnisTr)
      PhenoTr(:,1)=PhenoTr(:,1)*sqrt(VarY)+MeanY
      write(UnitEbv,"(a20,2(1x,a20))") "Id","Observed","Estimate"
      do i=1,nAnisTr
        write(UnitEbv,"(a20,2(1x,f20.10))") IdTr(i),PhenoTr(i,1),Ebv(i,1)
      end do
      flush(UnitEbv)
      close(UnitEbv)

      call CalcCorrelation(PhenoTr(:,1),Ebv(:,1),Cor)
      DumC="Train"
      write(UnitSum,"(a14,3(1x,f12.4),3(1x,f20.10))") adjustl(DumC),Cor%Cor,Cor%Cov/Cor%Var2,Cor%Cov/Cor%Var1,Cor%Cov,Cor%Var1,Cor%Var2

      deallocate(Ebv)

      ! Rescale marker estimates back to observed phenotype and genotype scale

      ! y     = mu     + b*x     + e
      ! y     = mu     + b*x/SDX + e
      ! y     = mu     + b*z     + e
      ! y/SDY = mu/SDY + b*z/SDY + e/SDY
      ! y'    = mu'    + b'*z    + e'
      ! b'*z=b /SDY*z
      ! b*z =b'*SDY*z
      ! b*x/SDX =b'*SDY*x/SDX
      ! b=b'*SDY/SDX

      G(:,1)=G(:,1)*sqrt(VarY)/ScaleCoef(:)

      ! Output
      open(newunit=UnitSnpSol,file="SnpEstimates.txt",status="unknown")
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

        call CalcCorrelation(PhenoTe(1:nAnisTe(Pop),1),Ebv(1:nAnisTe(Pop),1),Cor)
        DumC="Predict"//Int2Char(Pop)
        write(UnitSum,"(a14,3(1x,f12.4),3(1x,f20.10))") adjustl(DumC),Cor%Cor,Cor%Cov/Cor%Var2,Cor%Cov/Cor%Var1,Cor%Cov,Cor%Var1,Cor%Var2

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
  if (trim(Method)=="RidgeMCMC") then
    call RidgeRegressionMCMC
  end if
  if (trim(Method)=="Ridge") then
    call RidgeRegression
  end if
  ! if (trim(Method)=="BayesA") then
  !   ! call BayesA
  ! end if
  call Prediction

  call UnintitialiseIntelRNG

  call cpu_time(Finish)

  write(STDOUT,"(a,f20.4,a)") "Time duration of AlphaBayes: ",Finish-Start," seconds"
  write(STDOUT,"(a)") " "
end program

!###############################################################################
