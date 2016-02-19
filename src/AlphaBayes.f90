#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

!###############################################################################

module Global
    implicit none

    integer :: idum,nSnp,nSnpExternal,nAnisTr,nRound,nBurn,nProcessors,ScalingOpt,nTePop
    integer,allocatable :: nAnisTe(:),FixedSnp(:),SnpPosition(:)

    real(4) :: VarY,VarA,VarE,Mu,Sum2pq
    real(4) :: ExpVarX,SumExpVarX,ObsVarX,SumObsVarX
    real(4),allocatable :: SnpTmp(:),AlleleFreq(:),GenosTr(:,:),PhenTr(:,:),G(:,:)

    character(len=1000) :: GenoTrFile,PhenoTrFile,FileFixedSnp,MarkerSolver
    character(len=1000),allocatable :: GenoTeFile(:),PhenoTeFile(:)

    contains

        subroutine PearsnR4(x,y,n,r,sxx,syy,sxy)
            implicit none
            ! Arguments
            integer,intent(in) :: n
            real(4),intent(in) :: x(n),y(n)
            real(4),intent(out) :: r
            real(4),intent(out) :: sxx
            real(4),intent(out) :: syy
            real(4),intent(out) :: sxy
            ! Other
            integer :: j
            real(4) :: TINY,ax,ay,xt,yt!,betai,prob,z,df,t

            tiny=1.0e-20
            ax=0.0
            ay=0.0
            DO j=1,n
                ax=ax+x(j)
                ay=ay+y(j)
            END DO
            ax=ax/float(n)
            ay=ay/float(n)
            sxx=0.
            syy=0.
            sxy=0.
            DO j=1,n
                xt=x(j)-ax
                yt=y(j)-ay
                sxx=sxx+xt*xt
                syy=syy+yt*yt
                sxy=sxy+xt*yt
            END DO
            r=sxy/(SQRT(sxx*syy)+TINY)
            !z=0.5*LOG(((1.+r)+TINY)/((1.-r)+TINY))
            !df=n-2
            !t=r*SQRT(df/(((1.-r)+TINY)*((1.+r)+TINY)))
            !prob=betai(0.5*df,0.5,df/(df+t**2))
            !prob=erfcc(ABS(z*SQRT(n-1.))/1.4142136)
            !prob=0
            return
        end subroutine PearsnR4

    !###########################################################################
end module Global

!###############################################################################

subroutine ReadParam
    use Global
    implicit none

    integer :: i,UnitSpec,UnitFixedSnp

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
        enddo
        nSnp=sum(FixedSnp(:))
        close(UnitFixedSnp)
    else
        nSnp=nSnpExternal
        allocate(FixedSnp(nSnpExternal))
        FixedSnp=1
    endif

    read(UnitSpec,*) DumC,nAnisTr
    read(UnitSpec,*) DumC,nAnisTe(:)

    read(UnitSpec,*) DumC,nRound
    read(UnitSpec,*) DumC,nBurn

    read(UnitSpec,*) DumC,VarA
    read(UnitSpec,*) DumC,VarE

    read(UnitSpec,*) DumC,nProcessors
    call OMP_SET_NUM_THREADS(nProcessors)

    read(UnitSpec,*) DumC,ScalingOpt

    read(UnitSpec,*) DumC,MarkerSolver
    close(UnitSpec)
end subroutine ReadParam

!###############################################################################

subroutine ReadData
    use Global
    implicit none

    integer :: i,j,nNotMissing,UnitGenoTr,UnitPhenoTr,UnitAlleleFreq

    real(4) :: ave,adev,sdev,var,skew,curt

    character(len=100) :: DumC

    open(newunit=UnitGenoTr,file=trim(GenoTrFile),status="old")
    open(newunit=UnitPhenoTr,file=trim(PhenoTrFile),status="old")
    open(newunit=UnitAlleleFreq,file="AlleleFreq.txt",status="unknown")

    allocate(SnpTmp(nSnpExternal))
    allocate(GenosTr(nAnisTr,nSnp))
    allocate(PhenTr(nAnisTr,1))
    allocate(G(nSnp,1))
    allocate(AlleleFreq(nSnp))
    allocate(SnpPosition(nSnp))

    j=0
    do i=1,nSnpExternal
        if (FixedSnp(i)==1) then
            j=j+1
            SnpPosition(j)=i
        endif
    enddo

    do i=1,nAnisTr
        read(UnitGenoTr,*) DumC,SnpTmp(:)
        GenosTr(i,:)=SnpTmp(SnpPosition(:))
        read(UnitPhenoTr,*) DumC,PhenTr(i,1)
    enddo

    call momentR4(PhenTr(:,1),nAnisTr,ave,adev,sdev,var,skew,curt)

    PhenTr(:,1)=(PhenTr(:,1)-ave)/sdev
    VarY=var

    SumExpVarX=0.0
    SumObsVarX=0.0
    Sum2pq=0.0
    do j=1,nSnp

        ! Compute allele freqs
        AlleleFreq(j)=0.0
        nNotMissing=0
        do i=1,nAnisTr
            if ((GenosTr(i,j)>-0.1).and.(GenosTr(i,j)<2.1)) then
                AlleleFreq(j)=AlleleFreq(j)+GenosTr(i,j)
                nNotMissing=nNotMissing+1
            endif
        enddo
        if (nNotMissing/=0) then
            AlleleFreq(j)=AlleleFreq(j)/float(2*nNotMissing)
        else
            AlleleFreq(j)=0.0
        endif
        write(UnitAlleleFreq,"(i8,f11.8)") j,AlleleFreq(j)

        ! Fix any odd data
        do i=1,nAnisTr
            if ((GenosTr(i,j)<-0.1).or.(GenosTr(i,j)>2.1)) then
                GenosTr(i,j)=2.0*AlleleFreq(j)
            endif
        enddo

        ! Standardize
        ExpVarX=2.0*(1.0-AlleleFreq(j))*AlleleFreq(j)+0.00001
        Sum2pq=Sum2pq+ExpVarX
        SumExpVarX=SumExpVarX+ExpVarX
        ObsVarX=var+0.00001
        SumObsVarX=SumObsVarX+ObsVarX

        ! ... center
        GenosTr(:,j)=GenosTr(:,j)-(2.0*AlleleFreq(j))

        ! ... scale
        if (ScalingOpt==2) then
            ! Scale by marker specific variance - expected
            ExpVarX=sqrt(ExpVarX)
            GenosTr(:,j)=GenosTr(:,j)/ExpVarX
        endif

        if (ScalingOpt==3) then
            ! Scale by marker specific variance - observed
            ObsVarX=sqrt(ObsVarX)
            GenosTr(:,j)=GenosTr(:,j)/ObsVarX
        endif

    enddo

    if (ScalingOpt==4) then
        ! Scale by average marker variance - expected
        ExpVarX=sqrt(SumExpVarX/float(nSnp))
        GenosTr(:,:)=GenosTr(:,:)/ExpVarX
    endif

    if (ScalingOpt==5) then
        ! Scale by average marker variance - observed
        ObsVarX=sqrt(SumObsVarX/float(nSnp))
        GenosTr(:,:)=GenosTr(:,:)/ObsVarX
    endif

    close(UnitGenoTr)
    close(UnitPhenoTr)
    close(UnitAlleleFreq)
end subroutine ReadData

!###############################################################################

subroutine RidgeRegression
    use Global
    implicit none

    integer :: h,j,snpid,RandomOrdering(nSnp)

    real(4) :: sdot,Eps,InvLhs,Rhs,Lhs,SolOld,OneR,ZeroR,nAnisTrR,LambdaSnp
    real(4),allocatable :: XpX(:,:),Xg(:,:),E(:,:)

    allocate(XpX(nSnp,1))
    allocate(Xg(nAnisTr,1))
    allocate(E(nAnisTr,1))

    OneR=1.0
    ZeroR=0.0
    Mu=0.0
    G=0.000001
    nAnisTrR=float(nAnisTr)

    ! Construct XpX and Lambda
    do j=1,nSnp
        XpX(j,1)=sdot(nAnisTr,GenosTr(:,j),1,GenosTr(:,j),1) + 0.000000000000001
    enddo

    LambdaSnp=VarE/(VarA/float(nSnp))

    ! Working phenotypes
    E(:,1)=PhenTr(:,1)

    !open(unit=6,form="formatted",CARRIAGECONTROL="FORTRAN")
    do h=1,nRound
        !             123456789 123456789 123456789 123456789 123456789 12
        !write(6,100) " Ridge regression (fixed variance components) round ",h
        !100 format("+",a52,i10)

        !Intercept
        E(:,1)=E(:,1)+Mu
        Rhs=sum(E(:,1))
        InvLhs=1.0/nAnisTrR
        Mu=Rhs*InvLhs
        E(:,1)=E(:,1)-Mu

        !Snp effects
        eps=0.0
        call random_order(RandomOrdering,nSnp,idum)
        do j=1,nSnp
            snpid=RandomOrdering(j)
            E(:,1)=E(:,1)+(GenosTr(:,snpid)*G(snpid,1))
            Lhs=XpX(snpid,1)+LambdaSnp
            Rhs=sdot(nAnisTr,GenosTr(:,snpid),1,E(:,1),1)
            SolOld=G(snpid,1)
            G(snpid,1)=Rhs/Lhs
            E(:,1)=E(:,1)-(GenosTr(:,snpid)*G(snpid,1))
            Eps=Eps+(G(snpid,1)-SolOld)**2
        enddo

        if (mod(h,200)==0) then
            call sgemm("n","n",nAnisTr,1,nSnp,OneR,GenosTr,nAnisTr,G,nSnp,ZeroR,Xg,nAnisTr)
            E(:,1)=PhenTr(:,1)-Xg(:,1)-Mu
        endif

        if (eps.lt.1e-8) then
            exit
        endif
    enddo
    !close(6)

    deallocate(XpX)
    deallocate(Xg)
    deallocate(E)
end subroutine RidgeRegression

!###############################################################################

subroutine RidgeRegressionMCMC
    use Global
    implicit none

    integer :: i,h,j,snpid,RandomOrdering(nSnp)

    real(4) :: sdot,InvLhs,Rhs,Lhs,OneR,ZeroR,nAnisTrR,LambdaSnp,EpE,VarG,GpG
    real(4) :: nSampR!,VarEAccum,VarGAccum
    real(4),allocatable :: GAccum(:,:),SX2(:),MX2(:),XpX(:,:),Xg(:,:),E(:,:)
    real(8) :: random_gamma,gasdev,R2,EDF0,EDF2,GDF0,GDF2,ES0,GS0,MSX

    allocate(XpX(nSnp,1))
    allocate(Xg(nAnisTr,1))
    allocate(GAccum(nSnp,1))
    allocate(SX2(nSnp))
    allocate(MX2(nSnp))
    allocate(E(nAnisTr,1))

    OneR=1.0
    ZeroR=0.0
    Mu=0.0
    G=0.1
    GAccum=0.0

    nAnisTrR=float(nAnisTr)
    R2=0.5

    EDF0=5.0
    EDF2=(float(nAnisTr)+EDF0)/2.0
    ES0=(VarY*(1.0-R2))*(EDF0+2.0)

    GDF0=5.0
    GDF2=(float(nSnp)+GDF0)/2.0
    SX2(:)=0.0
    MX2(:)=0.0
    do i=1,nSnp
        SX2(i)=sum(GenosTr(:,i)*GenosTr(:,i))
        MX2(i)=sum(GenosTr(:,i))/nAnisTrR
        MX2(i)=MX2(i)*MX2(i)
    enddo
    MSX=sum(SX2)/nAnisTrR-sum(MX2)
    GS0=((VarY*R2)/MSX)*(GDF0+2.0)

    ! Construct XpX
    do j=1,nSnp
        XpX(j,1)=sdot(nAnisTr,GenosTr(:,j),1,GenosTr(:,j),1) + 0.000000000000001
    enddo

    ! Working phenotypes
    E(:,1)=PhenTr(:,1)

    !open(unit=6,form="formatted",CARRIAGECONTROL="FORTRAN")
    nSampR=float(nRound-nBurn)
    do h=1,nRound
        !             123456789 123456789 123456789 123456789 123456789 123456
        !write(6,100) " Ridge regression (MCMC for all model parameters) round ",h
        !100 format("+",a56,i10)

        ! Residual variance
        EpE=sdot(nAnisTr,E(:,1),1,E(:,1),1) + ES0 + 0.000000000000001
        VarE=EpE/real(random_gamma(idum,EDF2,2.0d0,.true.))

        ! Snp variance
        GpG=sdot(nSnp,G(:,1),1,G(:,1),1) + GS0 + 0.000000000000001
        VarG=GpG/real(random_gamma(idum,GDF2,2.0d0,.true.))

        ! Ratio of variances
        LambdaSnp=VarE/VarG

        !Intercept
        E(:,1)=E(:,1)+Mu
        Rhs=sum(E(:,1))
        InvLhs=1.0/nAnisTrR
        Mu=Rhs*InvLhs
        E(:,1)=E(:,1)-Mu

        !Snp effects
        call random_order(RandomOrdering,nSnp,idum)
        do j=1,nSnp
            snpid=RandomOrdering(j)
            E(:,1)=E(:,1)+(GenosTr(:,snpid)*G(snpid,1))
            Lhs=XpX(snpid,1)+LambdaSnp
            Rhs=sdot(nAnisTr,GenosTr(:,snpid),1,E(:,1),1)
            G(snpid,1)=(Rhs/Lhs)+(gasdev(idum)*sqrt(1.0/Lhs))
            E(:,1)=E(:,1)-(GenosTr(:,snpid)*G(snpid,1))
        enddo
        if (h>nBurn) then
            GAccum=GAccum+G/nSampR
            !VarGAccum=VarGAccum+VarE/nSampR
            !VarEAccum=VarEAccum+VarG/nSampR
        endif

        ! Recompute residuals to avoid rounding errors
        if (mod(h,200)==0) then
            call sgemm("n","n",nAnisTr,1,nSnp,OneR,GenosTr,nAnisTr,G,nSnp,ZeroR,Xg,nAnisTr)
            E(:,1)=PhenTr(:,1)-Xg(:,1)-Mu
        endif

    enddo
    close(6)

    G=GAccum
    !VarG=VarGAccum
    !VarE=VarEAccum

    deallocate(XpX)
    deallocate(Xg)
    deallocate(GAccum)
    deallocate(SX2)
    deallocate(MX2)
    deallocate(E)
end subroutine RidgeRegressionMCMC

!###############################################################################

subroutine BayesA
    use Global
    implicit none

    integer :: h,j,snpid,RandomOrdering(nSnp)

    real(4) :: sdot,InvLhs,Rhs,Lhs,OneR,ZeroR,LambdaSnp,EpE,nSampR,nAnisTrR
    real(4),allocatable :: GAccum(:,:),GVar(:,:),XpX(:,:),Xg(:,:),E(:,:)
    real(8) :: random_gamma,gasdev,vdf,ShapeCoefficient,ScaleCoefficient,sc,gampe1,gampe2

    allocate(XpX(nSnp,1))
    allocate(Xg(nAnisTr,1))
    allocate(GVar(nSnp,1))
    allocate(GAccum(nSnp,1))
    allocate(E(nAnisTr,1))

    OneR=1.0
    ZeroR=0.0
    Mu=0.0
    G=0.1
    GAccum=0
    nAnisTrR=float(nAnisTr)

    vdf=4.012
    ShapeCoefficient=(vdf/2.0) !Ben
    ScaleCoefficient=((vdf-2.0)*VarA)/Sum2pq
    sc=2.0
    gampe1=(float(nAnisTr)/2.0)-2.0
    gampe2=2.0

    ! Construct XpX
    do j=1,nSnp
        XpX(j,1)=sdot(nAnisTr,GenosTr(:,j),1,GenosTr(:,j),1) + 0.000000000000001
    enddo

    ! Working phenotypes
    E(:,1)=PhenTr(:,1)

    !open(unit=6,form='formatted',CARRIAGECONTROL='FORTRAN')
    nSampR=float(nRound-nBurn)
    do h=1,nRound
        !             123456789 1234
        !write(6,100) " BayesA round ",h
        !100 format('+',a14,i10)

        ! Residual variance
        EpE=sdot(nAnisTr, E(:,1), 1, E(:,1), 1) + 0.000000000000001
        VarE=EpE/real(random_gamma(idum,gampe1,gampe2,.true.))

        ! Snp variances
        do j=1,nSnp
            GVar(j,1)=(ScaleCoefficient+(G(j,1)**2))/real(random_gamma(idum,ShapeCoefficient,sc,.true.))
        enddo

        !Intercept
        E(:,1)=E(:,1)+Mu
        Rhs=sum(E(:,1))
        InvLhs=1.0/nAnisTrR
        Mu=Rhs*InvLhs
        E(:,1)=E(:,1)-Mu

        !Snp effects
        call random_order(RandomOrdering,nSnp,idum)
        do j=1,nSnp
            snpid=RandomOrdering(j)
            E(:,1)=E(:,1)+(GenosTr(:,snpid)*G(snpid,1))
            LambdaSnp=VarE/gvar(snpid,1)
            Lhs=XpX(snpid,1)+LambdaSnp
            Rhs=sdot(nAnisTr,GenosTr(:,snpid),1,E(:,1),1)
            G(snpid,1)=(Rhs/Lhs)+(gasdev(idum)*sqrt(1.0/Lhs))
            E(:,1)=E(:,1)-(GenosTr(:,snpid)*G(snpid,1))
        enddo
        if (h>nBurn) then
            GAccum=GAccum+G/nSampR
        endif

        ! Recompute residuals to avoid rounding errors
        if (mod(h,200)==0) then
            call sgemm('n','n',nAnisTr,1,nSnp,OneR,GenosTr,nAnisTr,G,nSnp,ZeroR,Xg,nAnisTr)
            E(:,1)=PhenTr(:,1)-Xg(:,1)-Mu
        endif

    enddo

    G=GAccum

    deallocate(XpX)
    deallocate(Xg)
    deallocate(GVar)
    deallocate(GAccum)
    deallocate(E)
end subroutine BayesA

!###############################################################################

subroutine MarkerEffectPostProcessing
    use Global
    implicit none
    integer :: i,j,UnitSnpSol

    ! Rescale back to phenotype scale
    G(:,1)=G(:,1)*sqrt(VarY)

    ! Output
    open(newunit=UnitSnpSol,file="SnpSolutions.txt",status="unknown")
    j=0
    do i=1,nSnpExternal
        if (FixedSnp(i)==1) then
            j=j+1
            write(UnitSnpSol,"(i20,f20.10)") i,G(j,1)
        else
            write(UnitSnpSol,"(i20,f20.10)") i,0.0
        endif
    enddo
    flush(UnitSnpSol)
    close(UnitSnpSol)
end subroutine MarkerEffectPostProcessing

!###############################################################################

subroutine Prediction
    use Global
    implicit none

    integer :: i,j,Pop,UnitCor,UnitGenoTe,UnitPhenoTe,UnitEbv
    integer,allocatable :: IdTe(:)

    real(4) :: OneR,ZeroR,Cor,Var1,Var2,Cov12
    real(4),allocatable :: Ebv(:,:),PhenoTe(:,:),GenosTe(:,:)

    character(len=100) :: DumC

    OneR=1.0
    ZeroR=0.0

    open(newunit=UnitCor,file="TbvEbvCorrelation.txt",status="unknown")
    do Pop=1,nTePop
        allocate(IdTe(nAnisTe(Pop)))
        allocate(GenosTe(nAnisTe(Pop),nSnp))
        allocate(PhenoTe(nAnisTe(Pop),1))
        allocate(Ebv(nAnisTe(Pop),1))

        open(newunit=UnitGenoTe,file=trim(GenoTeFile(Pop)),status="old")
        open(newunit=UnitPhenoTe,file=trim(PhenoTeFile(Pop)),status="old")

        do i=1,nAnisTe(Pop)
            read(UnitGenoTe,*) IdTe(i),SnpTmp(:)
            GenosTe(i,:)=SnpTmp(SnpPosition(:))
            read(UnitPhenoTe,*) DumC,PhenoTe(i,1)
        enddo

        do j=1,nSnp
            ! Fix any odd data
            do i=1,nAnisTe(Pop)
                if ((GenosTe(i,j)<-0.1).or.(GenosTe(i,j)>2.1)) then
                    GenosTe(i,j)=2.0*AlleleFreq(j)
                endif
            enddo

            ! Standardize

            ! ... center
            GenosTe(:,j)=GenosTe(:,j)-(2.0*AlleleFreq(j))

            ! ... scale
            if (ScalingOpt==2) then
                ! Scale by marker specific variance - expected
                GenosTe(:,j)=GenosTe(:,j)/ExpVarX
            endif

            if (ScalingOpt==3) then
                ! Scale by marker specific variance - observed
                GenosTe(:,j)=GenosTe(:,j)/ObsVarX
            endif
        enddo

        if (ScalingOpt==4) then
            ! Scale by average marker variance - expected
            GenosTe(:,:)=GenosTe(:,:)/ExpVarX
        endif

        if (ScalingOpt==5) then
            ! Scale by average marker variance - observed
            GenosTe(:,:)=GenosTe(:,:)/ObsVarX
        endif

        close(UnitGenoTe)
        close(UnitPhenoTe)

        open(newunit=UnitEbv,file="Ebv.txt",status="unknown")
        call sgemm("n","n",nAnisTe(Pop),1,nSnp,OneR,GenosTe,nAnisTe(Pop),G,nSnp,ZeroR,Ebv,nAnisTe(Pop))
        do i=1,nAnisTe(Pop)
            write(UnitEbv,"(i20,2f20.10)") IdTe(i),Ebv(i,1)
        enddo
        flush(UnitEbv)
        close(UnitEbv)

        call PearsnR4(PhenoTe(:,1),Ebv(:,1),nAnisTe(Pop),Cor,Var1,Var2,Cov12)
        write(UnitCor,"(i4,3f7.3,3f20.10)") Pop,Cor,Cov12/Var2,Cov12/Var1,Var1/float(nAnisTe(Pop)),Var2/float(nAnisTe(Pop)),Cov12/float(nAnisTe(Pop))

        deallocate(IdTe)
        deallocate(GenosTe)
        deallocate(PhenoTe)
        deallocate(Ebv)
    enddo
    flush(UnitCor)
    close(UnitCor)
end subroutine Prediction

!###############################################################################

SUBROUTINE momentR4(DATA,n,ave,adev,sdev,var,skew,curt)
    IMPLICIT NONE
    INTEGER n,j
    real(4) :: adev,ave,curt,sdev,skew,var,DATA(n)
    real(4) :: p,s,ep
    IF (n.le.1) PAUSE 'n must be at least 2 in moment'
    s=0.0
    DO j=1,n
        s=s+DATA(j)
    END DO

    ave=s/float(n)
    adev=0.0
    var=0.0
    skew=0.0
    curt=0.0
    ep=0.0

    DO j=1,n
        s=DATA(j)-ave
        ep=ep+s
        adev=adev+ABS(s)
        p=s*s
        var=var+p
        p=p*s
        skew=skew+p
        p=p*s
        curt=curt+p
    END DO

    adev=adev/float(n)
    var=(var-ep**2.0/float(n))/(float(n)-1.0)
    sdev=SQRT(var)
    IF (var.ne.0) then
        skew=skew/(float(n)*sdev**3.0)
        curt=curt/(float(n)*var**2.0)-3.0
    ELSE
        !PRINT*, 'no skew or kurtosis when zero variance in moment'
        !PAUSE 'no skew or kurtosis when zero variance in moment'
    END IF
    RETURN
END SUBROUTINE momentR4

!###############################################################################

! This Function returns a uniform random deviate between 0.0 and 1.0.
! Set IDUM to any negative value to initialize or reinitialize the sequence.
!MODIFIED FOR REAL

FUNCTION ran1(idum)
IMPLICIT NONE
 INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
 DOUBLE PRECISION ran1,AM,EPS,RNMX
 PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
 INTEGER j,k,iv(NTAB),iy
 SAVE iv,iy
 DATA iv /NTAB*0/, iy /0/
  IF (idum.le.0.or.iy.eq.0) then
      idum=max(-idum,1)
  DO 11 j=NTAB+8,1,-1
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
  IF (idum.lt.0) idum=idum+IM
  IF (j.le.NTAB) iv(j)=idum

11 CONTINUE
     iy=iv(1)
  END IF
     k=idum/IQ
     idum=IA*(idum-k*IQ)-IR*k
  IF (idum.lt.0) idum=idum+IM
     j=1+iy/NDIV
     iy=iv(j)
     iv(j)=idum
     ran1=min(AM*iy,RNMX)
  RETURN
END

!###############################################################################

FUNCTION random_gamma(idum,s, b, first) RESULT(fn_val)
    IMPLICIT NONE

    ! Adapted from Fortran 77 code from the book:
    !     Dagpunar, J. 'Principles of random variate generation'
    !     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

    !     N.B. This version is in `double precision' and includes scaling

    !     FUNCTION GENERATES A RANDOM GAMMA VARIATE.
    !     CALLS EITHER random_gamma1 (S > 1.0)
    !     OR random_exponential (S = 1.0)
    !     OR random_gamma2 (S < 1.0).

    !     S = SHAPE PARAMETER OF DISTRIBUTION (0 < REAL).
    !     B = Scale parameter

    !IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER :: idum

    DOUBLE PRECISION, INTENT(IN)  :: s, b
    LOGICAL, INTENT(IN)    :: first
    DOUBLE PRECISION              :: fn_val

    ! Local parameters
    DOUBLE PRECISION, PARAMETER  :: one = 1.0_dp, zero = 0.0_dp

    IF (s <= zero) THEN
      WRITE(*, *) 'SHAPE PARAMETER VALUE MUST BE POSITIVE'
      STOP
    END IF

    IF (s >= one) THEN
      fn_val = random_gamma1(s, first)
    ELSE IF (s < one) THEN
      fn_val = random_gamma2(s, first)
    END IF

    ! Now scale the random variable
    fn_val = b * fn_val
    RETURN

    CONTAINS

    FUNCTION random_gamma1(s, first) RESULT(fn_val)
        IMPLICIT NONE

        ! Adapted from Fortran 77 code from the book:
        !     Dagpunar, J. 'Principles of random variate generation'
        !     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

        ! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
        ! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO GAMMA**(S-1)*EXP(-GAMMA),
        ! BASED UPON BEST'S T DISTRIBUTION METHOD

        !     S = SHAPE PARAMETER OF DISTRIBUTION
        !          (1.0 < REAL)

        DOUBLE PRECISION, INTENT(IN)  :: s
        LOGICAL, INTENT(IN)    :: first
        DOUBLE PRECISION              :: fn_val

        !     Local variables
        DOUBLE PRECISION             :: d, r, g, f, x
        DOUBLE PRECISION, SAVE       :: b, h
        DOUBLE PRECISION, PARAMETER  :: sixty4 = 64.0_dp, three = 3.0_dp, pt75 = 0.75_dp,  &
                                 two = 2.0_dp, half = 0.5_dp
        DOUBLE PRECISION :: ran1

        IF (s <= one) THEN
          WRITE(*, *) 'IMPERMISSIBLE SHAPE PARAMETER VALUE'
          STOP
        END IF

        IF (first) THEN                        ! Initialization, if necessary
          b = s - one
          h = SQRT(three*s - pt75)
        END IF

        DO
          r=ran1(idum)
          g = r - r*r
          IF (g <= zero) CYCLE
          f = (r - half)*h/SQRT(g)
          x = b + f
          IF (x <= zero) CYCLE
          r=ran1(idum)
          d = sixty4*g*(r*g)**2
          IF (d <= zero) EXIT
          IF (d*x < x - two*f*f) EXIT
          IF (LOG(d) < two*(b*LOG(x/b) - f)) EXIT
        END DO
        fn_val = x

        RETURN
    END FUNCTION random_gamma1

    FUNCTION random_gamma2(s, first) RESULT(fn_val)
        IMPLICIT NONE

        ! Adapted from Fortran 77 code from the book:
        !     Dagpunar, J. 'Principles of random variate generation'
        !     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

        ! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
        ! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO
        ! GAMMA2**(S-1) * EXP(-GAMMA2),
        ! USING A SWITCHING METHOD.

        !    S = SHAPE PARAMETER OF DISTRIBUTION
        !          (REAL < 1.0)

        DOUBLE PRECISION, INTENT(IN)  :: s
        LOGICAL, INTENT(IN)    :: first
        DOUBLE PRECISION              :: fn_val

        !     Local variables
        DOUBLE PRECISION            :: r, x, w
        DOUBLE PRECISION, SAVE       :: a, p, c, uf, vr, d
        DOUBLE PRECISION, PARAMETER  :: vsmall = EPSILON(one)
        DOUBLE PRECISION :: ran1

        IF (s <= zero .OR. s >= one) THEN
          WRITE(*, *) 'SHAPE PARAMETER VALUE OUTSIDE PERMITTED RANGE'
          STOP
        END IF

        IF (first) THEN                        ! Initialization, if necessary
          a = one - s
          p = a/(a + s*EXP(-a))
          IF (s < vsmall) THEN
            WRITE(*, *) 'SHAPE PARAMETER VALUE TOO SMALL'
            PRINT*, s
            STOP
          END IF
          c = one/s
          uf = p*(vsmall/a)**s
          vr = one - vsmall
          d = a*LOG(a)
        END IF

        DO
          r=ran1(idum)
          IF (r >= vr) THEN
            CYCLE
          ELSE IF (r > p) THEN
            x = a - LOG((one - r)/(one - p))
            w = a*LOG(x)-d
          ELSE IF (r > uf) THEN
            x = a*(r/p)**c
            w = x
          ELSE
            fn_val = zero
            RETURN
          END IF

          r=ran1(idum)
          IF (one-r <= w .AND. r > zero) THEN
            IF (r*(w + one) >= one) CYCLE
            IF (-LOG(r) <= w) CYCLE
          END IF
          EXIT
        END DO

        fn_val = x
        RETURN
    END FUNCTION random_gamma2
END FUNCTION random_gamma

!###############################################################################

SUBROUTINE random_order(order,n,idum)
    IMPLICIT NONE

    !  (C) Copr. 1986-92 Numerical Recipes Software 6
    !     Generate a random ordering of the integers 1 ... n.

    INTEGER, INTENT(IN)  :: n
    INTEGER, INTENT(OUT) :: order(n)
    INTEGER :: idum
    DOUBLE PRECISION ran1

    !     Local variables

    INTEGER :: i, j, k
    double precision    :: wk

    DO i = 1, n
      order(i) = i
    END DO

    !     Starting at the end, swap the current last indicator with one
    !     randomly chosen from those preceeding it.

    DO i = n, 2, -1
      wk=ran1(idum)
      j = 1 + i * wk
      IF (j < i) THEN
        k = order(i)
        order(i) = order(j)
        order(j) = k
      END IF
    END DO

    RETURN
END SUBROUTINE random_order

!###############################################################################

FUNCTION gasdev(idum)
    IMPLICIT NONE
    !C USES ran1
    !Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
    !as the source of uniform deviates.

    INTEGER idum
    DOUBLE PRECISION :: gasdev
    INTEGER iset
    DOUBLE PRECISION fac,gset,rsq,v1,v2,ran1
    SAVE iset,gset
    DATA iset/0/
    rsq=0.
    if (idum.lt.0) iset=0
    if (iset.eq.0) then
        do while (rsq.ge.1..or.rsq.eq.0.)
            v1=2.*ran1(idum)-1.
            v2=2.*ran1(idum)-1.
            rsq=v1**2+v2**2
        enddo
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
    else
        gasdev=gset
        iset=0
    endif
    return
END

!###############################################################################

subroutine InitiateSeed
    use Global
    implicit none
    integer :: edum
    DOUBLE PRECISION :: W(1),GASDEV

    open (unit=3,file="Seed.txt",status="old")

    !READ AND WRITE SEED BACK TO FILE
    READ (3,*) idum
    W(1)=GASDEV(idum)
    !Code to write new seed to file
    IF (idum>=0) THEN
        edum=(-1*idum)
    ELSE
        edum=idum
    END IF
    REWIND (3)
    WRITE (3,*) edum
    idum=edum
end subroutine InitiateSeed

!###############################################################################

program AlphaBayes
	use Global
	implicit none

    write(6,"(a)") ""
    write(6,"(a30,a,a30)") " ","**********************"," "
    write(6,"(a30,a,a30)") " ","*                    *"," "
    write(6,"(a30,a,a30)") " ","*     AlphaBayes     *"," "
    write(6,"(a30,a,a30)") " ","*                    *"," "
    write(6,"(a30,a,a30)") " ","**********************"
    write(6,"(a30,a,a30)") " ","VERSION:"//TOSTRING(VERS)," "
    write(6,"(a)") ""
    write(6,"(a35,a)")     " ","No Liability"
    write(6,"(a25,a)")     " ","Bugs to John.Hickey@roslin.ed.ac.uk"
    write(6,"(a)") ""

	call InitiateSeed
	call ReadParam
	call ReadData
	if (trim(MarkerSolver)=="RidgeMCMC") call RidgeRegressionMCMC
	if (trim(MarkerSolver)=="Ridge") call RidgeRegression
	if (trim(MarkerSolver)=="BayesA") call BayesA
	call MarkerEffectPostProcessing
	call Prediction
end program AlphaBayes

!###############################################################################