*DECK BREANM
      SUBROUTINE BREANM(IPMAC1,NG,LX1,NMIX1,IMIX,ICODE,ISPH,ZKEFF,B2,
     1 ENER,XXX1,VOL1,FLX1,DC1,TOT1,CHI1,SIGF1,SCAT1,JXM,JXP,FHETXM,
     2 FHETXP,ADF1,NGET,ADFREF,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Implement the 1D FD-ANM reflector model.
*
*Copyright:
* Copyright (C) 2022 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPMAC1  nodal macrolib.
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
      INTEGER NG,LX1,NMIX1,IMIX(LX1),ICODE(2),ISPH,NGET,IPRINT
      REAL ZKEFF,B2,ENER(NG+1),XXX1(LX1+1),VOL1(NMIX1),FLX1(NMIX1,NG),
     1 DC1(NMIX1,NG),TOT1(NMIX1,NG),CHI1(NMIX1,NG),SIGF1(NMIX1,NG),
     2 SCAT1(NMIX1,NG,NG),JXM(NMIX1,NG),JXP(NMIX1,NG),FHETXM(NMIX1,NG),
     3 FHETXP(NMIX1,NG),ADF1(NMIX1,NG),ADFREF(NG)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE)
      CHARACTER HADF*8
      TYPE(C_PTR) JPMAC1,KPMAC1
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK,AFACTOR,BETA,WORK1,WORK2,
     1 WORK4,WORK5,VOLTOT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FDXM,FDXP,FHOMM,FHOMP,WORK3,
     1 ZCODE2
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: L,R
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(L(NG,2*NG,LX1),R(NG,2*NG,LX1))
      ALLOCATE(FHOMM(NMIX1,NG),FHOMP(NMIX1,NG),FDXM(NMIX1,NG),
     1 FDXP(NMIX1,NG),AFACTOR(NG),BETA(NG),WORK1(NG),WORK2(NG),
     2 WORK3(NG,NG),WORK4(NG),WORK5(NG),VOLTOT(NMIX1))
      ALLOCATE(ZCODE2(6,NG))
*----
*  COMPUTE BOUNDARY FLUXES
*----
      FDXM(:NMIX1,:NG)=0.0
      FDXP(:NMIX1,:NG)=0.0
      FHOMM(:NMIX1,:NG)=0.0
      FHOMP(:NMIX1,:NG)=0.0
      VOLTOT(:NMIX1)=0.0
      J_FUEL=0
      DO I=1,LX1
        IBM=IMIX(I)
        IF(IBM.EQ.0) CYCLE
        WORK1(:NG)=DC1(IBM,:NG)
        WORK3(:NG,:NG)=SCAT1(IBM,:NG,:NG)
        WORK4(:NG)=CHI1(IBM,:NG)
        WORK5(:NG)=SIGF1(IBM,:NG)
        DO IGR=1,NG
          IF(SIGF1(IBM,IGR).GT.0.0) J_FUEL=I
          DIFF=DC1(IBM,IGR)
          WORK2(IGR)=TOT1(IBM,IGR)+B2*DC1(IBM,IGR)-SCAT1(IBM,IGR,IGR)
        ENDDO
        VOL=XXX1(I+1)-XXX1(I)
        CALL NSSLR1(ZKEFF,NG,VOL,WORK1,WORK2,WORK3,WORK4,WORK5,
     1  L(1,1,I),R(1,1,I))
        !
        VOLTOT(IBM)=VOLTOT(IBM)+VOL
        FHOMM(IBM,:NG)=FHOMM(IBM,:NG)+REAL(MATMUL(L(:NG,:NG,I),
     1  FLX1(IBM,:NG))+MATMUL(L(:NG,NG+1:2*NG,I),JXM(IBM,:NG)),4)*VOL
        FHOMP(IBM,:NG)=FHOMP(IBM,:NG)+REAL(MATMUL(R(:NG,:NG,I),
     1  FLX1(IBM,:NG))+MATMUL(R(:NG,NG+1:2*NG,I),JXP(IBM,:NG)),4)*VOL
      ENDDO
      IF(IPRINT.GT.0) THEN
        WRITE(6,'(/23H BREANM: SURFACE FLUXES)')
        DO I=1,LX1
          IBM=IMIX(I)
          IF(IBM.EQ.0) CYCLE
          WRITE(6,'(/8H REGION=,I5)') I
          WRITE(6,20) 'fluxm',REAL(MATMUL(L(:NG,:NG,I),
     1    FLX1(IBM,:NG))+MATMUL(L(:NG,NG+1:2*NG,I),JXM(IBM,:NG)),4)
          WRITE(6,20) 'fluxp',REAL(MATMUL(R(:NG,:NG,I),
     1    FLX1(IBM,:NG))+MATMUL(R(:NG,NG+1:2*NG,I),JXP(IBM,:NG)),4)
        ENDDO
      ENDIF
      DO IBM=1,NMIX1
        DO IGR=1,NG
          FDXM(IBM,IGR)=VOLTOT(IBM)*FHETXM(IBM,IGR)/FHOMM(IBM,IGR)
          FDXP(IBM,IGR)=VOLTOT(IBM)*FHETXP(IBM,IGR)/FHOMP(IBM,IGR)
        ENDDO
      ENDDO
      IF(IPRINT.GT.0) THEN
        WRITE(6,'(/48H BREANM: DISCONTINUITY FACTORS BEFORE NORMALIZAT,
     1  3HION)')
        DO IBM=1,NMIX1
          WRITE(6,'(/9H MIXTURE=,I5)') IBM
          WRITE(6,20) 'FDXM',FDXM(IBM,:NG)
          WRITE(6,20) 'FDXP',FDXP(IBM,:NG)
        ENDDO
      ENDIF
*----
*  COMPUTE ALBEDOS
*----
      IF(ICODE(2).NE.0) THEN
        BETA(:)=0.0
        IBM=IMIX(LX1)
        DO IGR=1,NG
          IF(IBM.EQ.0) CYCLE
          AFACTOR(IGR)=FDXP(IBM,IGR)*JXP(IBM,IGR)/FHETXP(IBM,IGR)
          BETA(IGR)=(1.0-2.0*AFACTOR(IGR))/(1.0+2.0*AFACTOR(IGR))
        ENDDO
        IF(IPRINT.GT.0) THEN
          WRITE(6,'(/16H BREANM: ALBEDOS)')
          WRITE(6,20) 'BETA',BETA(:NG)
        ENDIF
      ENDIF
*----
*  NGET NORMALIZATION OF THE DISCONTINUITY FACTORS
*----
      IF(J_FUEL.GT.0) THEN
        IF(NGET.GT.0) THEN
          IBM=IMIX(J_FUEL)
          DO IGR=1,NG
            ! impose the adf on the fuel assembly side
            IF(IBM.EQ.0) CYCLE
            IF(NGET.EQ.1) THEN
              FNORM=ADFREF(IGR)/FDXP(IBM,IGR)
            ELSE
              FNORM=ADF1(IBM,IGR)/FDXP(IBM,IGR)
            ENDIF
            FDXP(IBM,IGR)=FDXP(IBM,IGR)*FNORM
            IF(J_FUEL<LX1) THEN
              IBMP=IMIX(J_FUEL+1)
              IF(IBMP.GT.0) FDXM(IBMP,IGR)=FDXM(IBMP,IGR)*FNORM
            ENDIF
          ENDDO
        ENDIF
        DO J=J_FUEL,1,-1
          IBM=IMIX(J)
          IF(IBM.EQ.0) CYCLE
          DO IGR=1,NG
            IF(J>1) THEN
              IBMM=IMIX(J-1)
              IF(IBMM.GT.0) FDXP(IBMM,IGR)=FDXP(IBMM,IGR)*FDXP(IBM,IGR)/
     1        FDXM(IBM,IGR)
            ENDIF
            FDXM(IBM,IGR)=FDXP(IBM,IGR)
          ENDDO
        ENDDO
      ENDIF
      DO J=J_FUEL+1,LX1
        IBM=IMIX(J)
        IF(IBM.EQ.0) CYCLE
        DO IGR=1,NG
          IF(J<LX1) THEN
            IBMP=IMIX(J+1)
            IF(IBMP.GT.0) FDXM(IBMP,IGR)=FDXM(IBMP,IGR)*FDXM(IBM,IGR)/
     1      FDXP(IBM,IGR)
          ENDIF
          FDXP(IBM,IGR)=FDXM(IBM,IGR)
        ENDDO
      ENDDO
      IF(IPRINT.GT.0) THEN
        WRITE(6,'(/48H BREANM: DISCONTINUITY FACTORS AFTER NGET NORMAL,
     1  7HIZATION)')
        DO IBM=1,NMIX1
          WRITE(6,'(/9H MIXTURE=,I5)') IBM
          WRITE(6,20) 'FDX',FDXM(IBM,:NG)
        ENDDO
      ENDIF
*----
*  APPLY SPH FACTORS
*----
      IF(ISPH.EQ.1) THEN
        DO IGR=1,NG
          DO J=1,LX1
            IBM=IMIX(J)
            IF(IBM.EQ.0) CYCLE
            TOT1(IBM,IGR)=TOT1(IBM,IGR)/FDXM(IBM,IGR)
            DC1(IBM,IGR)=DC1(IBM,IGR)/FDXM(IBM,IGR)
            SIGF1(IBM,IGR)=SIGF1(IBM,IGR)/FDXM(IBM,IGR)
            DO JGR=1,NG
              SCAT1(IBM,IGR,JGR)=SCAT1(IBM,IGR,JGR)/FDXM(IBM,JGR)
            ENDDO
          ENDDO
        ENDDO
        IF(ICODE(2).NE.0) THEN
          BETA(:)=0.0
          IF(ICODE(2).NE.0) THEN
            IBM=IMIX(LX1)
            DO IGR=1,NG
              IF(IBM.EQ.0) CYCLE
              AFACTOR(IGR)=AFACTOR(IGR)/FDXM(IBM,IGR)
              BETA(IGR)=(1.0-2.0*AFACTOR(IGR))/(1.0+2.0*AFACTOR(IGR))
            ENDDO
          ENDIF
          IF(IPRINT.GT.0) THEN
            WRITE(6,'(/30H BREANM: SPH CORRECTED ALBEDOS)')
            WRITE(6,20) 'BETA',BETA(:NG)
          ENDIF
        ENDIF
      ENDIF
      IF(IPRINT.GT.0) THEN
        WRITE(6,'(/31H BREANM: DIFFUSION COEFFICIENTS)')
        DO IBM=1,NMIX1
          WRITE(6,'(/9H MIXTURE=,I5)') IBM
          WRITE(6,20) 'DIFF',DC1(IBM,:NG)
        ENDDO
      ENDIF
*----
*  SAVE THE OUTPUT NODAL MACROLIB
*----
      ALLOCATE(IJJ(NMIX1),NJJ(NMIX1),IPOS(NMIX1),WORK(NMIX1*NG))
      ISTATE(:)=0
      ISTATE(1)=NG
      ISTATE(2)=NMIX1
      ISTATE(3)=1
      IF(J_FUEL.GT.0) ISTATE(4)=1
      IF(ICODE(2).NE.0) ISTATE(8)=1  ! physical albedo information
      ISTATE(9)=1  ! diffusion coefficient information
      IF(ISPH.EQ.0) ISTATE(12)=3 ! discontinuity factor information
      CALL LCMPUT(IPMAC1,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPMAC1,'ENERGY',NG+1,2,ENER)
      CALL LCMPUT(IPMAC1,'VOLUME',NMIX1,2,VOL1)
      IF(ICODE(2).NE.0) CALL LCMPUT(IPMAC1,'ALBEDO',NG,2,BETA)
      IF(ISPH.EQ.0) THEN
        CALL LCMSIX(IPMAC1,'ADF',1)
          NTYPE=1
          HADF='FD_B'
          CALL LCMPUT(IPMAC1,'NTYPE',1,1,NTYPE)
          CALL LCMPTC(IPMAC1,'HADF',8,1,HADF)
          CALL LCMPUT(IPMAC1,HADF,NMIX1*NG,2,FDXM)
        CALL LCMSIX(IPMAC1,' ',2)
      ENDIF
      JPMAC1=LCMLID(IPMAC1,'GROUP',NG)
      DO IGR=1,NG
        KPMAC1=LCMDIL(JPMAC1,IGR)
        DO IBM=1,NMIX1
          WORK(IBM)=VOL1(IBM)*FLX1(IBM,IGR)
        ENDDO
        CALL LCMPUT(KPMAC1,'FLUX-INTG',NMIX1,2,WORK)
        CALL LCMPUT(KPMAC1,'NTOT0',NMIX1,2,TOT1(:,IGR))
        CALL LCMPUT(KPMAC1,'DIFF',NMIX1,2,DC1(:,IGR))
        DO IBM=1,NMIX1
          WORK(IBM)=SCAT1(IBM,IGR,IGR)
        ENDDO
        CALL LCMPUT(KPMAC1,'SIGW00',NMIX1,2,WORK)
        CALL LCMPUT(KPMAC1,'CHI',NMIX1,2,CHI1(:,IGR))
        CALL LCMPUT(KPMAC1,'NUSIGF',NMIX1,2,SIGF1(:,IGR))
        IF(ISPH.EQ.1) THEN
          DO IBM=1,NMIX1
            WORK(IBM)=1.0/FDXM(IBM,IGR)
          ENDDO
          CALL LCMPUT(KPMAC1,'NSPH',NMIX1,2,WORK)
        ENDIF
        IPOSDE=0
        DO IBM=1,NMIX1
          J2=IGR
          J1=IGR
          DO JGR=1,NG
            IF(SCAT1(IBM,IGR,JGR).NE.0.0) THEN
              J2=MAX(J2,JGR)
              J1=MIN(J1,JGR)
            ENDIF
          ENDDO
          NJJ(IBM)=J2-J1+1
          IJJ(IBM)=J2
          IPOS(IBM)=IPOSDE+1
          DO JGR=J2,J1,-1
            IPOSDE=IPOSDE+1
            IF(IPOSDE.GT.NG*NMIX1) CALL XABORT('BREANM: SCAT OVERFLOW.')
            WORK(IPOSDE)=SCAT1(IBM,IGR,JGR)
          ENDDO
        ENDDO
        CALL LCMPUT(KPMAC1,'SCAT00',IPOSDE,2,WORK)
        CALL LCMPUT(KPMAC1,'NJJS00',NMIX1,1,NJJ)
        CALL LCMPUT(KPMAC1,'IJJS00',NMIX1,1,IJJ)
        CALL LCMPUT(KPMAC1,'IPOS00',NMIX1,1,IPOS)
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ZCODE2)
      DEALLOCATE(WORK,IPOS,NJJ,IJJ)
      DEALLOCATE(VOLTOT,WORK5,WORK4,WORK3,WORK2,WORK1,BETA,AFACTOR,
     1 FDXP,FDXM,FHOMP,FHOMM)
      DEALLOCATE(R,L)
      RETURN
   20 FORMAT(1X,A9,1P,10E12.4,/(10X,10E12.4))
      END
