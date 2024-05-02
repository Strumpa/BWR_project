*DECK BREERA
      SUBROUTINE BREERA(IPMAC1,NC,NG,LX1,NMIX1,IMIX,ICODE,ISPH,ZKEFF,
     1 B2,ENER,XXX1,VOL1,FLX1,DC1,TOT1,CHI1,SIGF1,SCAT1,JXM,JXP,FHETXM,
     2 FHETXP,ADF1,NGET,ADFREF,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Implement the 1D ERM-ANM reflector model.
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
* IPMAC1  nodal macrolib.
* NC      number of sn macrolibs.
* NG      number of energy groups.
* LX1     number of nodes in the reflector model.
* NMIX1   number of mixtures in the nodal calculation.
* IMIX    mix index of each node.
* ICODE   physical albedo index on each side of the domain.
* ISPH    SPH flag (=0: use discontinuity factors; =1: use SPH factors).
* ZKEFF   effective multiplication factor.
* B2      buckling.
* ENER    energy limits.
* XXX1    spatial mesh.
* VOL1    volumes.
* FLX1    averaged fluxes
* DC1     diffusion coefficients.
* TOT1    total cross sections.
* CHI1    fission spectra.
* SIGF1   nu*fission cross sections.
* SCAT1   scattering P0 cross sections.
* JXM     left boundary currents.
* JXP     right boundary currents.
* FHETXM  left boundary fluxes.
* FHETXP  right boundary fluxes.
* ADF1    assembly discontinuity factors from macrolib.
* NGET    type of NGET normalization if discontinuity factors
*         (=0: simple; =1: imposed ADF on fuel assembly; =2: recover
*         fuel assembly ADF from input macrolib).
* ADFREF  imposed ADF values on fuel assembly side.
* IPRINT  edition flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC1
      INTEGER NC,NG,LX1,NMIX1,IMIX(LX1),ICODE(2),ISPH,NGET,IPRINT
      REAL ZKEFF(NC),B2(NC),ENER(NG+1),XXX1(LX1+1),VOL1(NMIX1,NC),
     1 FLX1(NMIX1,NG,NC),DC1(NMIX1,NG,NC),TOT1(NMIX1,NG,NC),
     2 CHI1(NMIX1,NG,NC),SIGF1(NMIX1,NG,NC),SCAT1(NMIX1,NG,NG,NC),
     3 JXM(NMIX1,NG,NC),JXP(NMIX1,NG,NC),FHETXM(NMIX1,NG,NC),
     4 FHETXP(NMIX1,NG,NC),ADF1(NMIX1,NG,NC),ADFREF(NG)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE)
      CHARACTER(LEN=8) HADF(2)
      TYPE(C_PTR) JPMAC1,KPMAC1
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK1D,WORK1,WORK2,WORK4,WORK5,
     1 VOLTOT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FLX,DC,TOT,CHI,SIGF,
     1 ADF,AFACTOR,BETA,WORK3
      REAL, ALLOCATABLE, DIMENSION(:,:,:) ::FDXM,FDXP,SCAT
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: TAU,B,X
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: WORK2D
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: FHOMM,FHOMP,L,R
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WORK1(NG),WORK2(NG),WORK4(NG),WORK5(NG),VOLTOT(NMIX1),
     1 FLX(NMIX1,NG),DC(NMIX1,NG),TOT(NMIX1,NG),CHI(NMIX1,NG),
     1 SIGF(NMIX1,NG),ADF(NMIX1,NG),AFACTOR(NG,NG),BETA(NG,NG))
      ALLOCATE(FDXM(NMIX1,NG,NG),FDXP(NMIX1,NG,NG),SCAT(NMIX1,NG,NG),
     1 WORK3(NG,NG))
      ALLOCATE(FHOMM(NC,NG,NMIX1),FHOMP(NC,NG,NMIX1),L(NG,2*NG,LX1),
     1 R(NG,2*NG,LX1))
*----
*  AVERAGE THE OUTPUT NODAL MACROLIB
*----
      VOLTOT(:)=0.0
      FLX(:,:)=0.0
      DC(:,:)=0.0
      TOT(:,:)=0.0
      CHI(:,:)=0.0
      SIGF(:,:)=0.0
      SCAT(:,:,:)=0.0
      ADF(:,:)=0.0
      FHOMM(:NC,:NG,:NMIX1)=0.0D0
      FHOMP(:NC,:NG,:NMIX1)=0.0D0
      DO IC=1,NC
        DO IBM=1,NMIX1
          VOLTOT(IBM)=VOLTOT(IBM)+VOL1(IBM,IC)
          DO IGR=1,NG
            FLX(IBM,IGR)=FLX(IBM,IGR)+FLX1(IBM,IGR,IC)
            DC(IBM,IGR)=DC(IBM,IGR)+DC1(IBM,IGR,IC)
            TOT(IBM,IGR)=TOT(IBM,IGR)+TOT1(IBM,IGR,IC)
            CHI(IBM,IGR)=CHI(IBM,IGR)+CHI1(IBM,IGR,IC)
            SIGF(IBM,IGR)=SIGF(IBM,IGR)+SIGF1(IBM,IGR,IC)
            DO JGR=1,NG
              SCAT(IBM,IGR,JGR)=SCAT(IBM,IGR,JGR)+SCAT1(IBM,IGR,JGR,IC)
            ENDDO
            ADF(IBM,IGR)=ADF(IBM,IGR)+ADF1(IBM,IGR,IC)
          ENDDO
        ENDDO
      ENDDO
      VOLTOT(:)=VOLTOT(:)/REAL(NC)
      FLX(:,:)=FLX(:,:)/REAL(NC)
      DC(:,:)=DC(:,:)/REAL(NC)
      TOT(:,:)=TOT(:,:)/REAL(NC)
      CHI(:,:)=CHI(:,:)/REAL(NC)
      SIGF(:,:)=SIGF(:,:)/REAL(NC)
      SCAT(:,:,:)=SCAT(:,:,:)/REAL(NC)
      ADF(:,:)=ADF(:,:)/REAL(NC)
*----
*  LOOP OVER CASES
*----
      IF(ISPH.EQ.1) CALL XABORT('BREERA: SPH OPTION NOT IMPLEMENTED.')
      J_FUEL=0
      DO IC=1,NC
*----
*  SET AND SOLVE ANALYTIC NODAL SYSTEM
*----
        DO I=1,LX1
          IBM=IMIX(I)
          IF(IBM.EQ.0) CYCLE
          WORK1(:NG)=DC1(IBM,:NG,IC)
          WORK3(:NG,:NG)=SCAT1(IBM,:NG,:NG,IC)
          WORK4(:NG)=CHI1(IBM,:NG,IC)
          WORK5(:NG)=SIGF1(IBM,:NG,IC)
          DO IGR=1,NG
            IF(SIGF1(IBM,IGR,IC).GT.0.0) J_FUEL=I
            DIFF=DC1(IBM,IGR,IC)
            WORK2(IGR)=TOT1(IBM,IGR,IC)+B2(IC)*DC1(IBM,IGR,IC)-
     1      SCAT1(IBM,IGR,IGR,IC)
          ENDDO
          VOL=XXX1(I+1)-XXX1(I)
          CALL NSSLR1(ZKEFF(IC),NG,VOL,WORK1,WORK2,WORK3,WORK4,WORK5,
     1    L(1,1,I),R(1,1,I))
          !
          FHOMM(IC,:NG,IBM)=FHOMM(IC,:NG,IBM)+REAL(MATMUL(L(:NG,:NG,I),
     1    FLX1(IBM,:NG,IC))+
     2    MATMUL(L(:NG,NG+1:2*NG,I),JXM(IBM,:NG,IC)),4)*VOL
          FHOMP(IC,:NG,IBM)=FHOMP(IC,:NG,IBM)+REAL(MATMUL(R(:NG,:NG,I),
     1    FLX1(IBM,:NG,IC))+
     2    MATMUL(R(:NG,NG+1:2*NG,I),JXP(IBM,:NG,IC)),4)*VOL
        ENDDO
        DO IBM=1,NMIX1
          FHOMM(IC,:NG,IBM)=FHOMM(IC,:NG,IBM)/VOLTOT(IBM)
          FHOMP(IC,:NG,IBM)=FHOMP(IC,:NG,IBM)/VOLTOT(IBM)
        ENDDO
        IF(IPRINT.GT.0) THEN
          WRITE(6,'(/39H BREERA: NODAL SURFACE FLUXES FOR CASE=,I5)') IC
          DO IBM=1,NMIX1
            WRITE(6,'(/9H MIXTURE=,I5)') IBM
            WRITE(6,20) 'FHOMM',FHOMM(IC,:NG,IBM)
            WRITE(6,20) 'FHOMP',FHOMP(IC,:NG,IBM)
          ENDDO
        ENDIF
*----
*  END OF LOOP OVER CASES
*----
      ENDDO
*----
*  COMPUTE DISCONTINUITY AND ALBEDO FACTORS
*----
      AFACTOR(:,:)=0.0
      DO IBM=1,NMIX1
        IF(NC.EQ.1) THEN
          ! DF-NEM approach
          FDXM(IBM,:,:)=0.0
          FDXP(IBM,:,:)=0.0
          DO IGR=1,NG
            FDXM(IBM,IGR,IGR)=FHETXM(IBM,IGR,1)/REAL(FHOMM(1,IGR,IBM))
            FDXP(IBM,IGR,IGR)=FHETXP(IBM,IGR,1)/REAL(FHOMP(1,IGR,IBM))
          ENDDO
          IF(IBM.EQ.NMIX1) THEN
            DO IGR=1,NG
              AFACTOR(IGR,IGR)=JXP(IBM,IGR,1)/REAL(FHOMP(1,IGR,IBM))
            ENDDO
          ENDIF
        ELSE IF(NC.LT.NG) THEN
          CALL XABORT('BREERA: DEGENERATE SYSTEM')
        ELSE IF(NC.EQ.NG) THEN
          ! ERM-ANM approach: linear system resolution
          ALLOCATE(WORK2D(NC,2*NG))
          DO IGR=1,NG
            DO IC=1,NC
              WORK2D(IC,IGR)=FHOMM(IC,IGR,IBM)
              WORK2D(IC,NG+IGR)=FHETXM(IBM,IGR,IC)
            ENDDO
          ENDDO
          CALL ALSBD(NC,NG,WORK2D,IER,NC)
          IF(IER.NE.0) CALL XABORT('BREERA: SINGULAR MATRIX(1).')
          DO IGR=1,NG
            DO IC=1,NC
              FDXM(IBM,IGR,IC)=REAL(WORK2D(IC,NG+IGR))
            ENDDO
          ENDDO
          DO IGR=1,NG
            DO IC=1,NC
              WORK2D(IC,IGR)=FHOMP(IC,IGR,IBM)
              WORK2D(IC,NG+IGR)=FHETXP(IBM,IGR,IC)
            ENDDO
          ENDDO
          CALL ALSBD(NC,NG,WORK2D,IER,NC)
          IF(IER.NE.0) CALL XABORT('BREERA: SINGULAR MATRIX(2).')
          DO IGR=1,NG
            DO IC=1,NC
              FDXP(IBM,IGR,IC)=REAL(WORK2D(IC,NG+IGR))
            ENDDO
          ENDDO
          IF(IBM.EQ.NMIX1) THEN
            DO IGR=1,NG
              DO IC=1,NC
                WORK2D(IC,IGR)=FHOMP(IC,IGR,IBM)
                WORK2D(IC,NG+IGR)=JXP(IBM,IGR,IC)
              ENDDO
            ENDDO
            CALL ALSBD(NC,NG,WORK2D,IER,NC)
            IF(IER.NE.0) CALL XABORT('BREERA: SINGULAR MATRIX(3).')
            DO IGR=1,NG
              DO JGR=1,NG
                AFACTOR(IGR,JGR)=REAL(WORK2D(JGR,NG+IGR))
              ENDDO
            ENDDO
          ENDIF
          DEALLOCATE(WORK2D)
        ELSE IF(NC.GE.NG) THEN
          ! ERM-ANM approach: pseudo inversion
          ALLOCATE(TAU(NG),B(NC),X(NG))
          CALL ALST2F(NC,NC,NG,FHOMM(1,1,IBM),TAU)
          DO IGR=1,NG
            B(:)=FHETXM(IBM,IGR,:)
            CALL ALST2S(NC,NC,NG,FHOMM(1,1,IBM),TAU,B,X)
            FDXM(IBM,IGR,:)=REAL(X(:))
          ENDDO
          CALL ALST2F(NC,NC,NG,FHOMP(1,1,IBM),TAU)
          DO IGR=1,NG
            B(:)=FHETXP(IBM,IGR,:)
            CALL ALST2S(NC,NC,NG,FHOMP(1,1,IBM),TAU,B,X)
            FDXP(IBM,IGR,:)=REAL(X(:))
          ENDDO
          IF(IBM.EQ.NMIX1) THEN
            DO IGR=1,NG
              B(:)=JXP(IBM,IGR,:)
              CALL ALST2S(NC,NC,NG,FHOMP(1,1,IBM),TAU,B,X)
              AFACTOR(IGR,:)=REAL(X(:))
            ENDDO
          ENDIF
          DEALLOCATE(X,B,TAU)
        ENDIF
      ENDDO
      IF(IPRINT.GT.0) THEN
        WRITE(6,'(/48H BREERA: DISCONTINUITY FACTORS BEFORE NORMALIZAT,
     1  3HION)')
        DO IBM=1,NMIX1
          WRITE(6,'(/9H MIXTURE=,I5)') IBM
          WRITE(6,20) 'FDXM',FDXM(IBM,:NG,:NG)
          WRITE(6,20) 'FDXP',FDXP(IBM,:NG,:NG)
        ENDDO
      ENDIF
*----
*  COMPUTE ALBEDOS
*----
      IF(ICODE(2).NE.0) THEN
        BETA(:,:)=0.0
        DO IGR=1,NG
          DO JGR=1,NG
            BETA(IGR,JGR)=(1.0-2.0*AFACTOR(IGR,JGR))/(1.0+2.0*
     1      AFACTOR(IGR,JGR))
          ENDDO
        ENDDO
        IF(IPRINT.GT.0) THEN
          WRITE(6,'(/16H BREERA: ALBEDOS)')
          WRITE(6,20) 'BETA',BETA(:NG,:NG)
        ENDIF
      ENDIF
*----
*  NGET NORMALIZATION OF THE DISCONTINUITY FACTORS
*----
      ALLOCATE(WORK2D(NG,2*NG))
      DO J=1,LX1-1
        IBM=IMIX(J)
        IBMP=IMIX(J+1)
        IF((IBM.EQ.0).OR.(IBMP.EQ.0)) CYCLE
        DO IGR=1,NG
          DO JGR=1,NG
            WORK2D(IGR,JGR)=FDXP(IBM,IGR,JGR)
            WORK2D(IGR,NG+JGR)=FDXM(IBMP,IGR,JGR)
          ENDDO
        ENDDO
        CALL ALSBD(NG,NG,WORK2D,IER,NG)
        IF(IER.NE.0) CALL XABORT('BREERA: SINGULAR MATRIX(3).')
        DO IGR=1,NG
          ! impose the adf on the fuel assembly side
          IF((J.EQ.J_FUEL).AND.(NGET.EQ.1)) THEN
            FNORM=ADFREF(IGR)
          ELSE IF((J.EQ.J_FUEL).AND.(NGET.EQ.2)) THEN
            FNORM=ADF(IBM,IGR)
          ELSE
            FNORM=FDXP(IBM,IGR,IGR)
          ENDIF
          FDXP(IBM,IGR,:)=0.0
          FDXP(IBM,IGR,IGR)=FNORM
          DO JGR=1,NG
            FDXM(IBMP,IGR,JGR)=REAL(WORK2D(IGR,NG+JGR))*FNORM
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(WORK2D)
      IF(J_FUEL.GT.0) THEN
        DO J=J_FUEL,1,-1
          IBM=IMIX(J)
          IF(IBM.EQ.0) CYCLE
          DO IGR=1,NG
            FNORM=FDXP(IBM,IGR,IGR)/FDXM(IBM,IGR,IGR)
            DO JGR=1,NG
              IF(J>1) THEN
               IBMM=IMIX(J-1)
               IF(IBMM.GT.0) FDXP(IBMM,IGR,JGR)=FDXP(IBMM,IGR,JGR)*FNORM
              ENDIF
              FDXM(IBM,IGR,JGR)=FDXM(IBM,IGR,JGR)*FNORM
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      DO J=J_FUEL+1,LX1
        IBM=IMIX(J)
        IF(IBM.EQ.0) CYCLE
        DO IGR=1,NG
          FNORM=FDXM(IBM,IGR,IGR)/FDXP(IBM,IGR,IGR)
          DO JGR=1,NG
            IF(J<LX1) THEN
              IBMP=IMIX(J+1)
              IF(IBMP.GT.0) FDXM(IBMP,IGR,JGR)=FDXM(IBMP,IGR,JGR)*FNORM
            ENDIF
            FDXP(IBM,IGR,JGR)=FDXP(IBM,IGR,JGR)*FNORM
          ENDDO
        ENDDO
      ENDDO
      IF(IPRINT.GT.0) THEN
        WRITE(6,'(/48H BREERA: DISCONTINUITY FACTORS AFTER NGET NORMAL,
     1  7HIZATION)')
        DO IBM=1,NMIX1
          WRITE(6,'(/9H MIXTURE=,I5)') IBM
          WRITE(6,20) 'FDXM',FDXM(IBM,:NG,:NG)
          WRITE(6,20) 'FDXP',FDXP(IBM,:NG,:NG)
        ENDDO
      ENDIF
*----
*  SAVE THE OUTPUT NODAL MACROLIB
*----
      ALLOCATE(IJJ(NMIX1),NJJ(NMIX1),IPOS(NMIX1),WORK1D(NMIX1*NG))
      ISTATE(:)=0
      ISTATE(1)=NG
      ISTATE(2)=NMIX1
      ISTATE(3)=1
      IF(J_FUEL.GT.0) ISTATE(4)=1
      IF(ICODE(2).NE.0) ISTATE(8)=1  ! physical matrix albedo info
      ISTATE(9)=1  ! diffusion coefficient information
      IF(ISPH.EQ.0) ISTATE(12)=4 ! discontinuity factor information
      CALL LCMPUT(IPMAC1,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPMAC1,'ENERGY',NG+1,2,ENER)
      CALL LCMPUT(IPMAC1,'VOLUME',NMIX1,2,VOLTOT)
      IF(ICODE(2).NE.0) CALL LCMPUT(IPMAC1,'ALBEDO',NG*NG,2,BETA)
      IF(ISPH.EQ.0) THEN
        CALL LCMSIX(IPMAC1,'ADF',1)
          NTYPE=2
          HADF(1)='ERM_M'
          HADF(2)='ERM_P'
          CALL LCMPUT(IPMAC1,'NTYPE',1,1,NTYPE)
          CALL LCMPTC(IPMAC1,'HADF',8,NTYPE,HADF)
          CALL LCMPUT(IPMAC1,HADF(1),NMIX1*NG*NG,2,FDXM)
          CALL LCMPUT(IPMAC1,HADF(2),NMIX1*NG*NG,2,FDXP)
        CALL LCMSIX(IPMAC1,' ',2)
      ENDIF
      JPMAC1=LCMLID(IPMAC1,'GROUP',NG)
      DO IGR=1,NG
        KPMAC1=LCMDIL(JPMAC1,IGR)
        DO IBM=1,NMIX1
          WORK1D(IBM)=VOLTOT(IBM)*FLX(IBM,IGR)
        ENDDO
        CALL LCMPUT(KPMAC1,'FLUX-INTG',NMIX1,2,WORK1D)
        CALL LCMPUT(KPMAC1,'NTOT0',NMIX1,2,TOT(:,IGR))
        CALL LCMPUT(KPMAC1,'DIFF',NMIX1,2,DC(:,IGR))
        DO IBM=1,NMIX1
          WORK1D(IBM)=SCAT(IBM,IGR,IGR)
        ENDDO
        CALL LCMPUT(KPMAC1,'SIGW00',NMIX1,2,WORK1D)
        CALL LCMPUT(KPMAC1,'CHI',NMIX1,2,CHI(:,IGR))
        CALL LCMPUT(KPMAC1,'NUSIGF',NMIX1,2,SIGF(:,IGR))
        IPOSDE=0
        DO IBM=1,NMIX1
          J2=IGR
          J1=IGR
          DO JGR=1,NG
            IF(SCAT(IBM,IGR,JGR).NE.0.0) THEN
              J2=MAX(J2,JGR)
              J1=MIN(J1,JGR)
            ENDIF
          ENDDO
          NJJ(IBM)=J2-J1+1
          IJJ(IBM)=J2
          IPOS(IBM)=IPOSDE+1
          DO JGR=J2,J1,-1
            IPOSDE=IPOSDE+1
            IF(IPOSDE.GT.NG*NMIX1) CALL XABORT('BREERA: SCAT OVERFLOW.')
            WORK1D(IPOSDE)=SCAT(IBM,IGR,JGR)
          ENDDO
        ENDDO
        CALL LCMPUT(KPMAC1,'SCAT00',IPOSDE,2,WORK1D)
        CALL LCMPUT(KPMAC1,'NJJS00',NMIX1,1,NJJ)
        CALL LCMPUT(KPMAC1,'IJJS00',NMIX1,1,IJJ)
        CALL LCMPUT(KPMAC1,'IPOS00',NMIX1,1,IPOS)
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(R,L,FHOMP,FHOMM)
      DEALLOCATE(SCAT,FDXP,FDXM)
      DEALLOCATE(WORK3,BETA,AFACTOR,ADF,SIGF,CHI,TOT,DC,FLX,VOLTOT,
     1 WORK5,WORK4,WORK2,WORK1)
      RETURN
   20 FORMAT(1X,A9,1P,10E12.4,/(10X,10E12.4))
      END
