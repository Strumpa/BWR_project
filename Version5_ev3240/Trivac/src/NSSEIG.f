*DECK NSSEIG
      SUBROUTINE NSSEIG(NMAX,NMAY,NMAZ,LL4F,NDIM,NEL,NMIX,NG,MAT,IDL,
     > VOL,MUX,MUY,MUZ,IMAX,IMAY,IMAZ,IPY,IPZ,CHI,SIGF,SCAT,A11X,A11Y,
     > A11Z,EPSTHR,MAXTHR,NADI,EPSOUT,MAXOUT,ICL1,ICL2,ITER,EVECT,
     > FKEFF,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of a multigroup eigenvalue system for the calculation of the
* direct neutron flux in Trivac. Use the preconditioned power method
* with a two-parameter SVAT acceleration technique. CMFD solution.
*
*Copyright:
* Copyright (C) 2023 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NMAX    first dimension of array A11X.
* NMAY    first dimension of array A11Y.
* NMAZ    first dimension of array A11Z.
* LL4F    number of unknowns per energy group.
* NDIM    number of dimensions (1, 2 or 3).
* NEL     number of nodes.
* NMIX    number of mixtures in the nodal calculation.
* NG      number of energy groups.
* MAT     material mixtures.
* IDL     position of averaged fluxes in unknown vector.
* VOL     node volumes.
* MUX     X-oriented compressed storage mode indices.
* MUY     Y-oriented compressed storage mode indices.
* MUZ     Z-oriented compressed storage mode indices.
* IMAX    X-oriented position of each first non-zero column element.
* IMAY    Y-oriented position of each first non-zero column element.
* IMAZ    Z-oriented position of each first non-zero column element.
* IPY     Y-oriented permutation matrices.
* IPZ     Z-oriented permutation matrices.
* CHI     fission spectra.
* SIGF    nu times fission cross section.
* SCAT    scattering cross section.
* A11X    X-oriented sparse coefficient matrix.
* A11Y    Y-oriented sparse coefficient matrix.
* A11Z    Z-oriented sparse coefficient matrix.
* EPSTHR  thermal iteration epsilon.
* MAXTHR  maximum number of thermal iterations.
* NADI    number of inner ADI iterations.
* EPSOUT  convergence epsilon for the power method.
* MAXOUT  maximum number of iterations for the power method.
* ICL1    number of free iretations in one cycle of the up-scattering
*         iterations.
* ICL2    number of accelerated up-scattering iterations in one cycle.
* EVECT   initial estimate of fundamental eigenvalue.
* IMPX    print parameter.
* FKEFF   initial estimate of fundamental eigenvalue.
*
*Parameters: output
* ITER    number of iterations.
* EVECT   corresponding eigenvector.
* FKEFF   fundamental eigenvalue.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER, INTENT(IN) :: NMAX,NMAY,NMAZ,LL4F,NDIM,NEL,NMIX,NG,
     > MAT(NEL),IDL(NEL),MUX(LL4F),MUY(LL4F),MUZ(LL4F),IMAX(LL4F),
     > IMAY(LL4F),IMAZ(LL4F),IPY(LL4F),IPZ(LL4F),MAXTHR,NADI,MAXOUT,
     > ICL1,ICL2,IMPX
      REAL, INTENT(IN) :: VOL(NEL),CHI(NMIX,NG),SIGF(NMIX,NG),
     > SCAT(NMIX,NG,NG),A11X(NMAX,NG),A11Y(NMAY,NG),A11Z(NMAZ,NG)
      INTEGER, INTENT(OUT) :: ITER
      REAL, INTENT(IN) :: EPSTHR,EPSOUT
      REAL, INTENT(INOUT) :: EVECT(LL4F,NG),FKEFF
*----
*  LOCAL VARIABLES
*----
      REAL, PARAMETER :: EPS1=1.0E-5
      REAL(KIND=8), PARAMETER :: ALP_TAB(24) = (/ 0.2, 0.4, 0.6,
     1  0.8, 1.0, 1.2, 1.5, 2.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0,
     2  40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0 /)
      REAL(KIND=8), PARAMETER :: BET_TAB(11) = (/ -1.0, -0.8, -0.6,
     1 -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 /)
      REAL(KIND=8) :: AEAE,AEAG,AEAH,AGAG,AGAH,AHAH,BEBE,BEBG,BEBH,
     1 BGBG,BGBH,BHBH,AEBE,AEBG,AEBH,AGBE,AGBG,AGBH,AHBE,AHBG,AHBH,
     2 X,DXDA,DXDB,Y,DYDA,DYDB,Z,DZDA,DZDB,F,D2F(2,3),EVAL,ALP,BET,
     3 FMIN,VVV
      LOGICAL LOGTES
      CHARACTER(LEN=3) :: TEXT3
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: S2,F1,GARM1,GARM2
      REAL, ALLOCATABLE, DIMENSION(:,:) :: S,GRAD1,GRAD2,GAR1,GAR2,
     1 GAR3,GAF1,GAF2,GAF3
      REAL, ALLOCATABLE, DIMENSION(:,:) :: IA11X,IA11Y,IA11Z
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: WORK
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(S2(LL4F),S(LL4F,NG),GRAD1(LL4F,NG),GRAD2(LL4F,NG),
     > GAR1(LL4F,NG),GAR2(LL4F,NG),GAR3(LL4F,NG),GAF1(LL4F,NG),
     > GAF2(LL4F,NG),GAF3(LL4F,NG),WORK(LL4F,NG,3),IA11X(NMAX,NG),
     > IA11Y(NMAY,NG),IA11Z(NMAZ,NG))
*----
*  LU MATRIX FACTORIZATION
*----
      IA11X(:NMAX,:NG)=A11X(:NMAX,:NG)
      DO IG=1,NG
        CALL ALLUF(LL4F,IA11X(1,IG),MUX,IMAX)
      ENDDO
      IF(NDIM.GT.1) THEN
        IA11Y(:NMAY,:NG)=A11Y(:NMAY,:NG)
        DO IG=1,NG
          CALL ALLUF(LL4F,IA11Y(1,IG),MUY,IMAY)
        ENDDO
      ENDIF
      IF(NDIM.EQ.3) THEN
        IA11Z(:NMAZ,:NG)=A11Z(:NMAZ,:NG)
        DO IG=1,NG
          CALL ALLUF(LL4F,IA11Z(1,IG),MUZ,IMAZ)
        ENDDO
      ENDIF
*----
*  POWER METHOD
*----
      NCTOT=ICL1+ICL2
      IF(ICL2.EQ.0) THEN
         NCPTM=NCTOT+1
      ELSE
         NCPTM=ICL1
      ENDIF
      EVAL=1.0D0/FKEFF
      VVV=0.0D0
      ISTART=1
      NNADI=NADI
      TEST=0.0
      IF(IMPX.GE.2) WRITE (6,600) NADI
      ITER=0
      DO
        ITER=ITER+1
        IF(ITER > MAXOUT) CALL XABORT('NSSEIG: OUTER ITER. FAILURE.')
*----
*  EIGENVALUE EVALUATION
*----
        CALL NSSMPA(NMAX,NMAY,NMAZ,LL4F,NDIM,NEL,NMIX,NG,MAT,IDL,
     >  VOL,MUX,MUY,MUZ,IMAX,IMAY,IMAZ,IPY,IPZ,SCAT,A11X,A11Y,A11Z,
     >  EVECT,WORK(1,1,1))
        CALL NSSMPB(LL4F,NEL,NMIX,NG,MAT,IDL,VOL,CHI,SIGF,EVECT,
     >  WORK(1,1,2))
        AEBE=0.0D0
        BEBE=0.0D0
        DO IG=1,NG
          DO I=1,LL4F
            AEBE=AEBE+WORK(I,IG,1)*WORK(I,IG,2)
            BEBE=BEBE+WORK(I,IG,2)**2
          ENDDO
        ENDDO
        EVAL=AEBE/BEBE
        S(:LL4F,:NG)=REAL(EVAL)*WORK(:LL4F,:NG,2)-WORK(:LL4F,:NG,1)
*----
*  PERFORM THERMAL (UP-SCATTERING) ITERATIONS
*----
        WORK(:LL4F,:NG,:3)=0.0D0
        IGDEB=1
        TEXT3='NO '
        JTER=1
        ALLOCATE(F1(LL4F),GARM1(LL4F),GARM2(LL4F))
        DO
          WORK(:LL4F,:NG,1)=WORK(:LL4F,:NG,2)
          WORK(:LL4F,:NG,2)=WORK(:LL4F,:NG,3)
          WORK(:LL4F,:NG,3)=0.0D0
          GRAD1(:LL4F,:NG)=0.0D0
          DO IG=IGDEB,NG
            S2(:LL4F)=S(:LL4F,IG)
            DO JG=1,NG
              IF(JG.EQ.IG) CYCLE
              DO IEL=1,NEL
                IBM=MAT(IEL)
                IF(IBM.LE.0) CYCLE
                IND=IDL(IEL)
                IF(IND.EQ.0) CYCLE
                S2(IND)=S2(IND)+VOL(IEL)*SCAT(IBM,IG,JG)*GRAD1(IND,JG)
              ENDDO
            ENDDO
*
            WORK(:LL4F,IG,3)=0.0
            DO IADI=1,NNADI
              IF(IADI.EQ.1) THEN
                F1(:LL4F)=S2(:LL4F)
              ELSE
*               scalar multiplication for a x-oriented matrix.
                CALL ALLUM(LL4F,A11X(1,IG),WORK(1,IG,3),F1,MUX,IMAX,1)
                IF(NDIM.GE.2) THEN
*                  scalar multiplication for a y-oriented matrix.
                   GARM1(IPY(:LL4F))=WORK(:LL4F,IG,3)
                   GARM2(IPY(:LL4F))=F1(:LL4F)
                   CALL ALLUM(LL4F,A11Y(1,IG),GARM1,GARM2,MUY,IMAY,2)
                   F1(:LL4F)=GARM2(IPY(:LL4F))
                ENDIF
                IF(NDIM.EQ.3) THEN
*                  scalar multiplication for a z-oriented matrix.
                   GARM1(IPZ(:LL4F))=WORK(:LL4F,IG,3)
                   GARM2(IPZ(:LL4F))=F1(:LL4F)
                   CALL ALLUM(LL4F,A11Z(1,IG),GARM1,GARM2,MUZ,IMAZ,2)
                   F1(:LL4F)=GARM2(IPZ(:LL4F))
                ENDIF
                F1(:LL4F)=S2(:LL4F)-F1(:LL4F)
              ENDIF
*             scalar solution for a x-oriented linear system.
              CALL ALLUS(LL4F,MUX,IMAX,IA11X(1,IG),F1)
              IF(NDIM.GE.2) THEN
*               scalar solution for a y-oriented linear system.
                DO I=1,LL4F
                  II=IPY(I)
                  GARM1(II)=F1(I)*A11Y(MUY(II),IG)
                ENDDO
                CALL ALLUS(LL4F,MUY,IMAY,IA11Y(1,IG),GARM1)
                F1(:LL4F)=GARM1(IPY(:LL4F))
              ENDIF
              IF(NDIM.EQ.3) THEN
*               scalar solution for a z-oriented linear system.
                DO I=1,LL4F
                  II=IPZ(I)
                  GARM1(II)=F1(I)*A11Z(MUZ(II),IG)
                ENDDO
                CALL ALLUS(LL4F,MUZ,IMAZ,IA11Z(1,IG),GARM1)
                F1(:LL4F)=GARM1(IPZ(:LL4F))
              ENDIF
              WORK(:LL4F,IG,3)=WORK(:LL4F,IG,3)+F1(:LL4F)
              GRAD1(:LL4F,IG)=WORK(:LL4F,IG,3)
            ENDDO
          
          ENDDO
          IF(MAXTHR.EQ.0) EXIT
          IF(MOD(JTER-1,NCTOT).GE.NCPTM) THEN
            CALL NSS2AC(NG,LL4F,IGDEB,WORK,ZMU)
          ELSE
            ZMU=1.0D0
          ENDIF
          IGDEBO=IGDEB
          DO IG=IGDEBO,NG
            GINN=0.0D0
            FINN=0.0D0
            DO I=1,LL4F
              GINN=MAX(GINN,ABS(WORK(I,IG,2)-WORK(I,IG,3)))
              FINN=MAX(FINN,ABS(WORK(I,IG,3)))
            ENDDO
            GINN=GINN/FINN
            IF((GINN.LT.EPSTHR).AND.(IGDEB.EQ.IG)) IGDEB=IGDEB+1
          ENDDO
          IF(GINN.LT.EPSTHR) TEXT3='YES'
          IF(IMPX.GT.2) WRITE(6,610) JTER,GINN,EPSTHR,IGDEB,ZMU,TEXT3
          IF((GINN.LT.EPSTHR).OR.(JTER.EQ.MAXTHR)) EXIT
          JTER=JTER+1
        ENDDO
        DEALLOCATE(GARM2,GARM1,F1)
*----
*  DISPLACEMENT EVALUATION
*----
        F=0.0D0
        DELS=ABS(REAL((EVAL-VVV)/EVAL))
        VVV=EVAL
*----
*  EVALUATION OF THE TWO ACCELERATION PARAMETERS ALP AND BET
*----
        ALP=1.0D0
        BET=0.0D0
        N=0
        AEAE=0.0D0
        AEAG=0.0D0
        AEAH=0.0D0
        AGAG=0.0D0
        AGAH=0.0D0
        AHAH=0.0D0
        BEBG=0.0D0
        BEBH=0.0D0
        BGBG=0.0D0
        BGBH=0.0D0
        BHBH=0.0D0
        AEBG=0.0D0
        AEBH=0.0D0
        AGBE=0.0D0
        AGBG=0.0D0
        AGBH=0.0D0
        AHBE=0.0D0
        AHBG=0.0D0
        AHBH=0.0D0
        CALL NSSMPA(NMAX,NMAY,NMAZ,LL4F,NDIM,NEL,NMIX,NG,MAT,IDL,VOL,
     >  MUX,MUY,MUZ,IMAX,IMAY,IMAZ,IPY,IPZ,SCAT,A11X,A11Y,A11Z,EVECT,
     >  GAR1)
        CALL NSSMPA(NMAX,NMAY,NMAZ,LL4F,NDIM,NEL,NMIX,NG,MAT,IDL,VOL,
     >  MUX,MUY,MUZ,IMAX,IMAY,IMAZ,IPY,IPZ,SCAT,A11X,A11Y,A11Z,GRAD1,
     >  GAR2)
        IF(1+MOD(ITER-ISTART,ICL1+ICL2).GT.ICL1) THEN
          CALL NSSMPB(LL4F,NEL,NMIX,NG,MAT,IDL,VOL,CHI,SIGF,EVECT,GAF1)
          CALL NSSMPB(LL4F,NEL,NMIX,NG,MAT,IDL,VOL,CHI,SIGF,GRAD1,GAF2)
          CALL NSSMPB(LL4F,NEL,NMIX,NG,MAT,IDL,VOL,CHI,SIGF,GRAD2,GAF3)
          DO IG=1,NG
            DO I=1,LL4F
*             COMPUTE (A ,A )
              AEAE=AEAE+GAR1(I,IG)**2
              AEAG=AEAG+GAR1(I,IG)*GAR2(I,IG)
              AEAH=AEAH+GAR1(I,IG)*GAR3(I,IG)
              AGAG=AGAG+GAR2(I,IG)**2
              AGAH=AGAH+GAR2(I,IG)*GAR3(I,IG)
              AHAH=AHAH+GAR3(I,IG)**2
*             COMPUTE (B ,B )
              BEBG=BEBG+GAF1(I,IG)*GAF2(I,IG)
              BEBH=BEBH+GAF1(I,IG)*GAF3(I,IG)
              BGBG=BGBG+GAF2(I,IG)**2
              BGBH=BGBH+GAF2(I,IG)*GAF3(I,IG)
              BHBH=BHBH+GAF3(I,IG)**2
*             COMPUTE (A ,B )
              AEBG=AEBG+GAR1(I,IG)*GAF2(I,IG)
              AEBH=AEBH+GAR1(I,IG)*GAF3(I,IG)
              AGBE=AGBE+GAR2(I,IG)*GAF1(I,IG)
              AGBG=AGBG+GAR2(I,IG)*GAF2(I,IG)
              AGBH=AGBH+GAR2(I,IG)*GAF3(I,IG)
              AHBE=AHBE+GAR3(I,IG)*GAF1(I,IG)
              AHBG=AHBG+GAR3(I,IG)*GAF2(I,IG)
              AHBH=AHBH+GAR3(I,IG)*GAF3(I,IG)
            ENDDO
          ENDDO
*
  210     N=N+1
          IF(N.GT.10) GO TO 215
*         COMPUTE X(ITER+1)
          X=BEBE+ALP*ALP*BGBG+BET*BET*BHBH+2.0D0*(ALP*BEBG+BET*BEBH
     >    +ALP*BET*BGBH)
          DXDA=2.0D0*(BEBG+ALP*BGBG+BET*BGBH)
          DXDB=2.0D0*(BEBH+ALP*BGBH+BET*BHBH)
*         COMPUTE Y(ITER+1)
          Y=AEAE+ALP*ALP*AGAG+BET*BET*AHAH+2.0D0*(ALP*AEAG+BET*AEAH
     >    +ALP*BET*AGAH)
          DYDA=2.0D0*(AEAG+ALP*AGAG+BET*AGAH)
          DYDB=2.0D0*(AEAH+ALP*AGAH+BET*AHAH)
*         COMPUTE Z(ITER+1)
          Z=AEBE+ALP*ALP*AGBG+BET*BET*AHBH+ALP*(AEBG+AGBE)
     >    +BET*(AEBH+AHBE)+ALP*BET*(AGBH+AHBG)
          DZDA=AEBG+AGBE+2.0D0*ALP*AGBG+BET*(AGBH+AHBG)
          DZDB=AEBH+AHBE+ALP*(AGBH+AHBG)+2.0D0*BET*AHBH
*         COMPUTE F(ITER+1)
          F=X*Y-Z*Z
          D2F(1,1)=2.0D0*(BGBG*Y+DXDA*DYDA+X*AGAG-DZDA**2-2.0D0*Z*AGBG)
          D2F(1,2)=2.0D0*BGBH*Y+DXDA*DYDB+DXDB*DYDA+2.0D0*X*AGAH
     >    -2.0D0*DZDA*DZDB-2.0D0*Z*(AGBH+AHBG)
          D2F(2,2)=2.0D0*(BHBH*Y+DXDB*DYDB+X*AHAH-DZDB**2-2.0D0*Z*AHBH)
          D2F(2,1)=D2F(1,2)
          D2F(1,3)=DXDA*Y+X*DYDA-2.0D0*Z*DZDA
          D2F(2,3)=DXDB*Y+X*DYDB-2.0D0*Z*DZDB
*         SOLUTION OF A LINEAR SYSTEM.
          CALL ALSBD(2,1,D2F,IER,2)
          IF(IER.NE.0) GO TO 215
          ALP=ALP-D2F(1,3)
          BET=BET-D2F(2,3)
          IF(ALP.GT.100.0D0) GO TO 215
          IF((ABS(D2F(1,3)).LE.1.0D-4).AND.(ABS(D2F(2,3)).LE.1.0D-4))
     >    GO TO 220
          GO TO 210
*
*         alternative algorithm in case of Newton-Raphton failure
  215     IF(IMPX.GT.0) WRITE(6,'(/30H NSSEIG: FAILURE OF THE NEWTON,
     >    55H-RAPHTON ALGORIHTHM FOR COMPUTING THE OVERRELAXATION PA,
     >    9HRAMETERS.)')
          IAMIN=999
          IBMIN=999
          FMIN=HUGE(FMIN)
          DO IA=1,SIZE(ALP_TAB)
            ALP=ALP_TAB(IA)
            DO IB=1,SIZE(BET_TAB)
              BET=BET_TAB(IB)
*             COMPUTE X
              X=BEBE+ALP*ALP*BGBG+BET*BET*BHBH+2.0D0*(ALP*BEBG+BET*BEBH
     >        +ALP*BET*BGBH)
*             COMPUTE Y
              Y=AEAE+ALP*ALP*AGAG+BET*BET*AHAH+2.0D0*(ALP*AEAG+BET*AEAH
     >        +ALP*BET*AGAH)
*             COMPUTE Z
              Z=AEBE+ALP*ALP*AGBG+BET*BET*AHBH+ALP*(AEBG+AGBE)
     >        +BET*(AEBH+AHBE)+ALP*BET*(AGBH+AHBG)
*             COMPUTE F
              F=X*Y-Z*Z
              IF(F.LT.FMIN) THEN
                IAMIN=IA
                IBMIN=IB
                FMIN=F
              ENDIF
            ENDDO
          ENDDO
          ALP=ALP_TAB(IAMIN)
          BET=BET_TAB(IBMIN)
  220     BET=BET/ALP
          IF((ALP.LT.1.0D0).AND.(ALP.GT.0.0D0)) THEN
            ALP=1.0D0
            BET=0.0D0
          ELSE IF(ALP.LE.0.0D0) THEN
            ISTART=ITER+1
            ALP=1.0D0
            BET=0.0D0
          ENDIF
          DO IG=1,NG
            DO I=1,LL4F
              GRAD1(I,IG)=REAL(ALP)*(GRAD1(I,IG)+REAL(BET)*GRAD2(I,IG))
              GAR2(I,IG)=REAL(ALP)*(GAR2(I,IG)+REAL(BET)*GAR3(I,IG))
            ENDDO
          ENDDO
        ENDIF
*
        LOGTES=(ITER.LT.ICL1).OR.(MOD(ITER-ISTART,ICL1+ICL2).EQ.ICL1-1)
        DELT=0.0D0
        IF(LOGTES.AND.(DELS.LE.EPS1)) THEN
          CALL NSSMPB(LL4F,NEL,NMIX,NG,MAT,IDL,VOL,CHI,SIGF,EVECT,GAF1)
          CALL NSSMPB(LL4F,NEL,NMIX,NG,MAT,IDL,VOL,CHI,SIGF,GRAD1,GAF2)
          DO IG=1,NG
            DELN=0.0D0
            DELD=0.0D0
            DO I=1,LL4F
              EVECT(I,IG)=EVECT(I,IG)+GRAD1(I,IG)
              GAR1(I,IG)=GAR1(I,IG)+GAR2(I,IG)
              GRAD2(I,IG)=GRAD1(I,IG)
              GAR3(I,IG)=GAR2(I,IG)
              DELN=MAX(DELN,ABS(GAF2(I,IG)))
              DELD=MAX(DELD,ABS(GAF1(I,IG)))
            ENDDO
            IF(DELD.NE.0.0D0) DELT=MAX(DELT,DELN/DELD)
          ENDDO
          IF(IMPX.GE.2) WRITE (6,620) ITER,AEAE,AEAG,AEAH,AGAG,AGAH,
     >    AHAH,BEBE,ALP,BET,EVAL,F,DELS,DELT,N,BEBG,BEBH,BGBG,BGBH,
     >    BHBH,AEBE,AEBG,AEBH,AGBE,AGBG,AGBH,AHBE,AHBG,AHBH
          IF(DELT.LE.EPSOUT) EXIT
        ELSE
          DO IG=1,NG
            DO I=1,LL4F
              EVECT(I,IG)=EVECT(I,IG)+GRAD1(I,IG)
              GAR1(I,IG)=GAR1(I,IG)+GAR2(I,IG)
              GRAD2(I,IG)=GRAD1(I,IG)
              GAR3(I,IG)=GAR2(I,IG)
            ENDDO
          ENDDO
          IF(IMPX.GE.2) WRITE (6,620) ITER,AEAE,AEAG,AEAH,AGAG,AGAH,
     >    AHAH,BEBE,ALP,BET,EVAL,F,DELS,DELT,N,BEBG,BEBH,BGBG,BGBH,
     >    BHBH,AEBE,AEBG,AEBH,AGBE,AGBG,AGBH,AHBE,AHBG,AHBH
        ENDIF
*
        IF(ITER.EQ.1) TEST=DELS
        IF((ITER.GT.5).AND.(DELS.GT.TEST)) CALL XABORT('NSSEIG: CONVER'
     >  //'GENCE FAILURE.')
        IF(ITER.GE.MAXOUT) THEN
          WRITE (6,630)
          EXIT
        ENDIF
        IF(MOD(ITER,36).EQ.0) THEN
           ISTART=ITER+1
           NNADI=NNADI+1
           IF(IMPX.GE.1) WRITE (6,650) NNADI
        ENDIF
      ENDDO
*----
*  FLUX NORMALIZATION
*----
      FMAX=MAXVAL(EVECT(:LL4F,:NG))
      EVECT(:LL4F,:NG)=EVECT(:LL4F,:NG)/FMAX
*----
*  SOLUTION EDITION
*----
      FKEFF=REAL(1.0D0/EVAL)
      IF(IMPX.GE.1) WRITE (6,640) ITER,FKEFF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IA11Z,IA11Y,IA11X)
      DEALLOCATE(WORK,GAF3,GAF2,GAF1,GAR3,GAR2,GAR1,GRAD2,GRAD1,S,S2)
      RETURN
*----
*  FORMATS
*----
  600 FORMAT(1H1/50H NSSEIG: ITERATIVE PROCEDURE BASED ON PRECONDITION,
     > 17HED POWER METHOD (,I2,37H ADI ITERATIONS PER OUTER ITERATION)./
     > 9X,16HDIRECT EQUATION.)
  610 FORMAT (10X,3HIN(,I3,6H) FLX:,5H PRC=,1P,E9.2,5H TAR=,E9.2,
     > 7H IGDEB=, I13,6H ACCE=,0P,F12.5,12H  CONVERGED=,A3)
  620 FORMAT(1X,I3,1P,7E9.1,0P,2F8.3,E14.6,3E10.2,I4/(4X,1P,7E9.1))
  630 FORMAT(/53H NSSEIG: ***WARNING*** THE MAXIMUM NUMBER OF OUTER IT,
     > 20HERATIONS IS REACHED.)
  640 FORMAT(/23H NSSEIG: CONVERGENCE IN,I4,12H ITERATIONS.//
     > 42H NSSEIG: EFFECTIVE MULTIPLICATION FACTOR =,1P,E17.10/)
  650 FORMAT(/53H NSSEIG: INCREASING THE NUMBER OF INNER ITERATIONS TO,
     1 I3,36H ADI ITERATIONS PER OUTER ITERATION./)
      !
      CONTAINS
        SUBROUTINE NSSMPA(NMAX,NMAY,NMAZ,LL4F,NDIM,NEL,NMIX,NG,MAT,IDL,
     >  VOL,MUX,MUY,MUZ,IMAX,IMAY,IMAZ,IPY,IPZ,SCAT,A11X,A11Y,A11Z,
     >  EVECT,S2)
        !
        ! A*EVECT MULTIPLICATION
        !
        INTEGER, INTENT(IN) :: NMAX,NMAY,NMAZ,LL4F,NDIM,NEL,NMIX,NG,
     >  MAT(NEL),IDL(NEL),MUX(LL4F),MUY(LL4F),MUZ(LL4F),IMAX(LL4F),
     >  IMAY(LL4F),IMAZ(LL4F),IPY(LL4F),IPZ(LL4F)
        REAL, INTENT(IN) :: VOL(NEL),SCAT(NMIX,NG,NG)
        REAL, INTENT(IN) :: EVECT(LL4F,NG),A11X(NMAX,NG),A11Y(NMAY,NG),
     >  A11Z(NMAZ,NG)
        REAL, INTENT(OUT) :: S2(LL4F,NG)
        REAL, ALLOCATABLE, DIMENSION(:) :: GAR1,GAR2
        !
        ALLOCATE(GAR1(LL4F),GAR2(LL4F))
        DO IG=1,NG
*         scalar multiplication for a x-oriented matrix.
          CALL ALLUM(LL4F,A11X(1,IG),EVECT(1,IG),S2(1,IG),MUX,IMAX,1)
          IF(NDIM.GE.2) THEN
*           scalar multiplication for a y-oriented matrix.
            GAR1(IPY(:LL4F))=EVECT(:LL4F,IG)
            GAR2(IPY(:LL4F))=S2(:LL4F,IG)
            CALL ALLUM(LL4F,A11Y(1,IG),GAR1,GAR2,MUY,IMAY,2)
            S2(:LL4F,IG)=GAR2(IPY(:LL4F))
          ENDIF
          IF(NDIM.EQ.3) THEN
*           scalar multiplication for a z-oriented matrix.
            GAR1(IPZ(:LL4F))=EVECT(:LL4F,IG)
            GAR2(IPZ(:LL4F))=S2(:LL4F,IG)
            CALL ALLUM(LL4F,A11Z(1,IG),GAR1,GAR2,MUZ,IMAZ,2)
            S2(:LL4F,IG)=GAR2(IPZ(:LL4F))
          ENDIF
          DO JG=1,NG
            IF(JG.EQ.IG) CYCLE
            DO IEL=1,NEL
              IBM=MAT(IEL)
              IF(IBM.LE.0) CYCLE
              IND=IDL(IEL)
              IF(IND.EQ.0) CYCLE
              S2(IND,IG)=S2(IND,IG)-VOL(IEL)*SCAT(IBM,IG,JG)*
     >        EVECT(IND,JG)
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(GAR2,GAR1)
        END SUBROUTINE NSSMPA
        !
        SUBROUTINE NSSMPB(LL4F,NEL,NMIX,NG,MAT,IDL,VOL,CHI,SIGF,EVECT,
     >  S2)
        !
        ! B*EVECT MULTIPLICATION
        !
        INTEGER, INTENT(IN) :: LL4F,NEL,NMIX,NG,MAT(NEL),IDL(NEL)
        REAL, INTENT(IN) :: VOL(NEL),CHI(NMIX,NG),SIGF(NMIX,NG)
        REAL, INTENT(IN) :: EVECT(LL4F,NG)
        REAL, INTENT(OUT) :: S2(LL4F,NG)
        !
        S2(:LL4F,:NG)=0.0D0
        DO IG=1,NG
          DO JG=1,NG ! IG <-- JG
            DO IEL=1,NEL
              IBM=MAT(IEL)
              IF(IBM.LE.0) CYCLE
              IND=IDL(IEL)
              IF(IND.EQ.0) CYCLE
              S2(IND,IG)=S2(IND,IG)+VOL(IEL)*CHI(IBM,IG)*SIGF(IBM,JG)*
     >        EVECT(IND,JG)
            ENDDO
          ENDDO
        ENDDO
        END SUBROUTINE NSSMPB
      END SUBROUTINE NSSEIG
