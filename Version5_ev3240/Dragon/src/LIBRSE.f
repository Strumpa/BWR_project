*DECK LIBRSE
      SUBROUTINE LIBRSE(IPLIB,IPTMP,MAXTRA,HNAMIS,LBIN,NGRP,NL,NED,
     1 NDEL,HVECT,NFS,IMPX,DELI,AWR,IALTER,SVDEPS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Process snapshots for the resonance spectrum expansion method.
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
* IPLIB   pointer to the isotopic directory in microlib.
* IPTMP   pointer to the multi-dilution internal library.
* MAXTRA  maximum number of energy bins of size DELI.
* HNAMIS  local name of the isotope:
*         HNAMIS(1:8)  is the local isotope name;
*         HNAMIS(9:12) is a suffix function of the mix number.
* LBIN    number of fine energy groups.
* NGRP    number of coarse energy groups.
* NL      number of Legendre orders required in the calculation
*         (NL=1 or higher).
* NED     number of extra vector edits.
* NDEL    number of delayed neutron precursor groups.
* HVECT   names of the extra vector edits.
* NFS     number of fine energy groups in each coarse energy group.
* IMPX    print flag (equal to zero for no print).
* DELI    elementary lethargy width.
* AWR     mass ratio for current isotope.
* IALTER  type of approximation (=0: use exponentials; =1: use Taylor
*         expansions).
* SVDEPS  rank accuracy of the singular value decomposition.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPTMP
      INTEGER MAXTRA,LBIN,NGRP,NL,NED,NDEL,NFS(NGRP),IMPX,IALTER
      REAL DELI,AWR,SVDEPS
      CHARACTER HNAMIS*12,HVECT(NED)*8
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPLIB1,JPLIB2,KPLIB,JPTMP,KPTMP
      PARAMETER (MAXITER=100,MAXNOR=12)
      CHARACTER TEXT12*12
      LOGICAL LGOLD
      DOUBLE PRECISION DQQ,GAR0,GAR1(3),FFF
      CHARACTER(LEN=4),DIMENSION(4),PARAMETER ::
     1   HPART=(/'NWT0','NTOT','SIGF','SIGS'/)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NJJ,MRANK
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISMIN,ISMAX,ISM
      REAL, ALLOCATABLE, DIMENSION(:) :: EBIN,UUU,DEL,STR,SIGS,SIGT,
     1 PRI,STIS,KSN,DILUT,DELTG,GOLD,SIGP
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FLUX,TOTAL,SIGF
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SIGSD,SADD,ZDEL
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: SCAT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: W,PHIGAR,DGAR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: V,XSDIL,TTT,DDD
      DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: U,BBB,SSIGS,SSIGT
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LSCAT,LADD
      COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: CT,CD
      TYPE VECTOR_ARRAY
        DOUBLE PRECISION, POINTER, DIMENSION(:) :: VECTOR
      END TYPE VECTOR_ARRAY
      TYPE MATRIX_ARRAY
        DOUBLE PRECISION, POINTER, DIMENSION(:,:) :: MATRIX
      END TYPE MATRIX_ARRAY
      TYPE(VECTOR_ARRAY), ALLOCATABLE, DIMENSION(:) :: UU_V,U_V,SIGT_V,
     1 XI_V,GAMMA_V,PSI_V
      TYPE(MATRIX_ARRAY), ALLOCATABLE, DIMENSION(:) :: PHI_M,BBB_M,U_M,
     1 V_M,SIGT_M,T_M,SIGP_M,EFF_M
      TYPE(MATRIX_ARRAY), ALLOCATABLE, DIMENSION(:,:) :: SCAT_M,TSCAT_M
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NJJ(NGRP),MRANK(NGRP))
      ALLOCATE(EBIN(LBIN+1),UUU(LBIN+1),DEL(LBIN),STR(LBIN),SIGS(LBIN),
     1 SIGT(LBIN),PRI(MAXTRA),STIS(LBIN),EFF_M(NGRP),PHI_M(NGRP),
     2 UU_V(NGRP),U_V(NGRP),U_M(NGRP),V_M(NGRP),SIGT_M(NGRP),
     3 SCAT_M(NGRP,NGRP),T_M(NGRP),SIGT_V(NGRP),XI_V(NGRP),
     4 GAMMA_V(NGRP),PSI_V(NGRP),TSCAT_M(NGRP,NGRP),SIGP_M(NGRP),
     5 ISM(2,NL))
      CALL LCMLEN(IPTMP,'ISOTOPESLIST',NDIL,ITYLCM)
      ALLOCATE(FLUX(NGRP,NDIL),TOTAL(NGRP,NDIL),SIGF(NGRP,NDIL),
     1 SIGSD(NGRP,NL,NDIL),SCAT(NGRP,NGRP,NL,NDIL),SADD(NGRP,NED,NDIL),
     2 ZDEL(NGRP,NDEL,NDIL),DELTG(NGRP),LSCAT(NL),LADD(NED),GOLD(NGRP),
     3 ISMIN(NL,NGRP),ISMAX(NL,NGRP))
*----
*  ALLOCATE DILUTION-DEPENDENT ARRAYS
*----
      ALLOCATE(DILUT(NDIL),KSN(NGRP*NDIL))
      CALL LCMGET(IPTMP,'ISOTOPESDSN',KSN)
      DO IDIL=1,NDIL
        DILUT(IDIL)=KSN((IDIL-1)*NGRP+1)
      ENDDO
      DEALLOCATE(KSN)
      IF(IMPX.GT.2) THEN
        WRITE(6,'(/32H LIBRSE: DILUTIONS FOR ISOTOPE '',A12,2H'':)')
     1  HNAMIS
        WRITE(6,'(1X,1P,12E12.4)') DILUT(:NDIL)
      ENDIF
*----
*  RECOVER INFORMATION FROM *TEMPORARY* LCM OBJECT.
*----
      NDIL=NDIL-1
      CALL LIBEXT(IPTMP,NGRP,NL,NDIL,NED,HVECT,NDEL,.TRUE.,IMPX,DILUT,
     1 MDIL,LSCAT,LSIGF,LADD,LGOLD,FLUX,TOTAL,SIGF,SIGSD,SCAT,SADD,
     2 ZDEL,DELTG,GOLD,ISMIN,ISMAX)
      NDIL=NDIL+1
*----
*  COPY INFINITE DILUTION DATA FROM IPTMP TO IPLIB.
*----
      JPTMP=LCMGID(IPTMP,'ISOTOPESLIST')
      CALL LCMLEL(JPTMP,NDIL,ILENG,ITYLCM)
      IF(ILENG.EQ.0) THEN
        TEXT12=HNAMIS(1:8)
        WRITE(TEXT12(9:12),'(I4.4)') NDIL
        CALL XABORT('LIBRSE: MISSING LIST ITEM FOR '//TEXT12)
      ENDIF
      KPTMP=LCMGIL(JPTMP,NDIL) ! set NDIL-th isotope
      CALL LCMLEN(KPTMP,'LAMBDA-D',NDEL,ITYLCM)
      CALL LCMEQU(KPTMP,IPLIB)
*
*     DESTROY THE MULTI-DILUTION INTERNAL LIBRARY.
      CALL LCMCL(IPTMP,2)
*----
*  RECOVER AUTOLIB DATA.
*----
      CALL LCMGET(IPLIB,'BIN-ENERGY',EBIN)
      CALL LCMGET(IPLIB,'BIN-NTOT0',SIGT)
      CALL LCMGET(IPLIB,'BIN-SIGS00',SIGS)
*----
*  NULLIFY POINTERS
*----
      DO IG=1,NGRP
        NULLIFY(SIGT_M(IG)%MATRIX)
        DO JG=1,NGRP
          NULLIFY(SCAT_M(IG,JG)%MATRIX)
          NULLIFY(TSCAT_M(IG,JG)%MATRIX)
        ENDDO
      ENDDO
*----
*  ELASTIC SCATTERING INFORMATION USED TO REBUILD THE SCAT MATRIX.
*----
      IBIN=0
      DELMIN=1.0E10
      UUU(1)=0.0
      DO IG=1,NGRP
        DO LI=1,NFS(IG)
          DELM=LOG(EBIN(IBIN+LI)/EBIN(IBIN+LI+1))
          UUU(IBIN+LI+1)=LOG(EBIN(1)/EBIN(IBIN+LI+1))
          DEL(IBIN+LI)=DELM
          DELMIN=MIN(DELMIN,DELM)
        ENDDO
        IBIN=IBIN+NFS(IG)
      ENDDO
      CALL LCMLEN(IPLIB,'BIN-DELI',LENGT,ITYLCM)
      IF((LENGT.EQ.1).AND.(ITYLCM.EQ.2)) THEN
        CALL LCMGET(IPLIB,'BIN-DELI',DELI)
      ELSE
        DELI=1.0/REAL(INT(1.00001/DELMIN))
      ENDIF
      PRI(:MAXTRA)=0.0
      CALL LIBPRI(MAXTRA,DELI,AWR,IALTER,0,NEXT,PRI)
*----
*  SOLVE FLUX CALCULATOR PROBLEMS FOR MANY DILUTIONS
*----
      LLL=0
      DO IG=1,NGRP
        LGBIN=NFS(IG)
        IF(LGBIN.EQ.0) CYCLE
        ALLOCATE(EFF_M(IG)%MATRIX(NDIL,3))
        ALLOCATE(PHI_M(IG)%MATRIX(NFS(IG),NDIL))
        LLL=LLL+LGBIN
      ENDDO
      ALLOCATE(PHIGAR(LBIN))
      DO IDIL=1,NDIL
        PHIGAR(:)=0.0D0
        LLL=0
        DO IG=1,NGRP
          IF(IMPX.GE.9) WRITE(6,'(29H LIBRSE: coarse energy group=,I8,
     1    10H dilution=,1P,E12.4)') IG,DILUT(IDIL)
          LGBIN=NFS(IG)
          IF(LGBIN.EQ.0) CYCLE
          LLL1=LLL
          GAR0=0.0D0
          GAR1(:3)=0.0D0
          DO LI=1,LGBIN
            LLL=LLL+1
            III=1
            STR(:LBIN)=0.0
            CALL LIBECT(MAXTRA,LLL,PRI,UUU(2),DELI,DEL,NEXT,III,MML,
     1      STIS)
            DO MM=1,MML
              STR(MM)=STR(MM)+STIS(MM)*SIGS(LLL-MM+1)
            ENDDO
            SIGMA=SIGT(LLL)-STR(1)
            DQQ=DILUT(IDIL)*DEL(LLL)
            DO MM=2,MIN(LLL,MML)
              DQQ=DQQ+STR(MM)*PHIGAR(LLL-MM+1)
            ENDDO
            PHIGAR(LLL)=DQQ/(DILUT(IDIL)+SIGMA)
            GAR0=GAR0+(UUU(LLL+1)-UUU(LLL))
            GAR1(1)=GAR1(1)+PHIGAR(LLL)
            GAR1(2)=GAR1(2)+SIGT(LLL)*PHIGAR(LLL)
            GAR1(3)=GAR1(3)+SIGS(LLL)*PHIGAR(LLL)
          ENDDO
          FFF=FLUX(IG,IDIL)*GAR0/GAR1(1)
          PHI_M(IG)%MATRIX(:LGBIN,IDIL)=PHIGAR(LLL1+1:LLL1+LGBIN)*FFF
          EFF_M(IG)%MATRIX(IDIL,:3)=GAR1(:3)*FFF/GAR0
        ENDDO
      ENDDO
      DEALLOCATE(PHIGAR)
*----
* LOOP OVER RESONANT GROUPS
*----
      MRANK(:NGRP)=0
      ALLOCATE(W(NDIL),V(NDIL,NDIL))
      LLL=0
      DO IG=1,NGRP
        NJJ(IG)=0
        LGBIN=NFS(IG)
        IF(LGBIN.EQ.0) CYCLE
        ALLOCATE(U(LGBIN,NDIL))
        U(:LGBIN,:NDIL)=PHI_M(IG)%MATRIX(:LGBIN,:NDIL)
        CALL ALSVDF(U,LGBIN,NDIL,LGBIN,NDIL,W,V)
        U_M(IG)%MATRIX=>U
        DO IDIL=1,NDIL
          IF(W(IDIL).LE.SVDEPS*DELI) THEN
            EXIT
          ELSE
            MRANK(IG)=MRANK(IG)+1
            IF(MRANK(IG).GT.20) CALL XABORT('LIBRSE: MRANK OVERFLOW.')
          ENDIF
        ENDDO
        MI=MRANK(IG)
        IF(IMPX.GE.2) THEN
          WRITE(6,'(/15H LIBRSE: RANK('',A12,2H'',,I4,2H)=,I3)') HNAMIS,
     1    IG,MI
        ENDIF
        IF(MI.EQ.0) CYCLE
        ALLOCATE(V_M(IG)%MATRIX(NDIL,MI),BBB_M(NGRP))
        DO JG=1,NGRP
          NULLIFY(BBB_M(JG)%MATRIX)
        ENDDO
        DO IM=1,MI
          V_M(IG)%MATRIX(:NDIL,IM)=V(:NDIL,IM)/W(IM)
        ENDDO
        !
        ! compute UU_V and U_V
        ALLOCATE(UU_V(IG)%VECTOR(MI),U_V(IG)%VECTOR(MI))
        DO IMR1=1,MI
          UU_V(IG)%VECTOR(IMR1)=SUM(U_M(IG)%MATRIX(:LGBIN,IMR1))
          U_V(IG)%VECTOR(IMR1)=0.0D0
          DO LI=1,LGBIN
            U_V(IG)%VECTOR(IMR1)=U_V(IG)%VECTOR(IMR1)+DEL(LLL+LI)*
     1      U_M(IG)%MATRIX(LI,IMR1)
          ENDDO
        ENDDO
        !
        ! compute SIGT_M
        ALLOCATE(SSIGT(MI,MI))
        SSIGT(:,:)=0.0D0
        DO LI=1,LGBIN
          SIGMA=SIGT(LLL+LI)
          DO IMR1=1,MI
            DO IMR2=1,MI
              SSIGT(IMR1,IMR2)=SSIGT(IMR1,IMR2)+U(LI,IMR1)*SIGMA*
     1        U(LI,IMR2)
            ENDDO
          ENDDO
        ENDDO
        SIGT_M(IG)%MATRIX=>SSIGT
        NULLIFY(SSIGT)
        !
        ! compute SCAT_M
        DO LI=1,LGBIN
          III=1
          STR(:LBIN)=0.0
          CALL LIBECT(MAXTRA,LLL+LI,PRI,UUU(2),DELI,DEL,NEXT,III,MML,
     1    STIS)
          JG=IG
          LJ=LI+1
          DO MM=1,MML
            LJ=LJ-1
            STR(MM)=STIS(MM)*SIGS(LLL+LI-MM+1)
   10       IF(.NOT.ASSOCIATED(BBB_M(JG)%MATRIX)) THEN
              NJJ(IG)=NJJ(IG)+1
              ALLOCATE(BBB(NFS(IG),NFS(JG)))
              BBB(:NFS(IG),:NFS(JG))=0.0D0
              BBB_M(JG)%MATRIX => BBB
            ENDIF
            IF(LJ.LE.0) THEN
              JG=JG-1
              IF(JG.LE.0) CALL XABORT('LIBRSE: SIGS MATRIX FAILURE(1).')
              LJ=LJ+NFS(JG)
              IF(LJ.LE.0) CALL XABORT('LIBRSE: SIGS MATRIX FAILURE(2).')
              GO TO 10
            ENDIF
            IF(.NOT.ASSOCIATED(BBB_M(JG)%MATRIX)) THEN
              CALL XABORT('LIBRSE: BBB_M%MATRIX NOT ASSOCIATED.')
            ENDIF
            BBB_M(JG)%MATRIX(LI,LJ)=STR(MM)
          ENDDO
        ENDDO
        DO JG=IG-NJJ(IG)+1,IG
          IF(.NOT.ASSOCIATED(U_M(JG)%MATRIX)) THEN
            CALL XABORT('LIBRSE: U_M%MATRIX NOT ASSOCIATED.')
          ENDIF
          BBB => BBB_M(JG)%MATRIX
          MJ=MRANK(JG)
          ALLOCATE(SSIGS(MI,MJ),DDD(LGBIN,MJ))
          DDD=MATMUL(BBB,U_M(JG)%MATRIX(:NFS(JG),:MJ))
          SSIGS=MATMUL(TRANSPOSE(U(:LGBIN,:MI)),DDD)
          SCAT_M(IG,JG)%MATRIX => SSIGS
          NULLIFY(SSIGS)
          DEALLOCATE(DDD)
        ENDDO
        DO JG=1,IG
          IF(ASSOCIATED(BBB_M(JG)%MATRIX)) THEN
            DEALLOCATE(BBB_M(JG)%MATRIX)
          ENDIF
        ENDDO
        DEALLOCATE(BBB_M)
*----
*  LINEAR TRANSFORMATION
*----
        ALLOCATE(T_M(IG)%MATRIX(MI,MI),SIGT_V(IG)%VECTOR(MI))
        ALLOCATE(XI_V(IG)%VECTOR(MI),GAMMA_V(IG)%VECTOR(MI),
     1  PSI_V(IG)%VECTOR(MI))
        ALLOCATE(CT(MI,MI),CD(MI,MI))
        CALL ALHQR(MI,MI,SIGT_M(IG)%MATRIX,MAXITER,ITER,CT,CD)
        DO LI=1,MI
          IF(AIMAG(CD(LI,LI)) /= 0.0D0) THEN
            CALL XABORT('LIBRSE: COMPLEX EIGENVALUE FOUND.')
          ENDIF
          SIGT_V(IG)%VECTOR(LI)=REAL(CD(LI,LI),8)
          DO LJ=1,MI
            T_M(IG)%MATRIX(LI,LJ)=REAL(CT(LI,LJ),8)
          ENDDO
        ENDDO
        DEALLOCATE(CD,CT)
        ALLOCATE(TTT(MI,MI))
        TTT(:MI,:MI)=TRANSPOSE(T_M(IG)%MATRIX(:MI,:MI))
        XI_V(IG)%VECTOR=MATMUL(TTT,UU_V(IG)%VECTOR)
        GAMMA_V(IG)%VECTOR=MATMUL(TTT,U_V(IG)%VECTOR)
        DO JG=1,IG
          IF(ASSOCIATED(SCAT_M(IG,JG)%MATRIX)) THEN
            MJ=MRANK(JG)
            ALLOCATE(TSCAT_M(IG,JG)%MATRIX(MI,MJ))
            ALLOCATE(SSIGS(MI,MJ),DDD(MI,MJ))
            DDD=MATMUL(SCAT_M(IG,JG)%MATRIX,T_M(JG)%MATRIX)
            SSIGS=MATMUL(TTT,DDD)
            TSCAT_M(IG,JG)%MATRIX => SSIGS
            NULLIFY(SSIGS)
            DEALLOCATE(DDD)
          ENDIF
        ENDDO
        DEALLOCATE(TTT)
        LLL=LLL+LGBIN
        IF(IMPX.GE.3) THEN
          WRITE(6,'(/34H LIBRSE: NWT0 BASE POINTS IN GROUP,I4,1H:)') IG
          WRITE(6,'(1X,1P,12E12.4)') XI_V(IG)%VECTOR(:)/DELTG(IG)
          WRITE(6,'(/34H LIBRSE: SIGT BASE POINTS IN GROUP,I4,1H:)') IG
          WRITE(6,'(1X,1P,12E12.4)') SIGT_V(IG)%VECTOR(:)
        ENDIF
      ENDDO
*----
*  COMPUTE RESONANCE SPECTRUM EXPANSION TABLES
*  XSDIL dilution dependent self-shielded cross sections:
*        XSDIL(1,:NDIL) self-shielded fluxes;
*        XSDIL(2,:NDIL) total self-shielded cross sections;
*        XSDIL(3,:NDIL) nu*fission self-shielded cross sections;
*        XSDIL(4,:NDIL) P0 scattering cross sections;
*        etc.
*        XSDIL(j,NDIL) are the infinite dilution values.
*----
      ALLOCATE(DGAR(NDIL))
      DO IG=1,NGRP
        MI=MRANK(IG)
        IF(MI.EQ.0) CYCLE
        NPART=3+NL+NED+NDEL
        DO IL=1,NL
          NPART=NPART+MAX(ISMAX(IL,IG)-ISMIN(IL,IG)+1,0)
        ENDDO
        ALLOCATE(SIGP_M(IG)%MATRIX(NPART,MI),XSDIL(NPART,NDIL))
        XSDIL(1,:NDIL)=EFF_M(IG)%MATRIX(:NDIL,1)
        XSDIL(2,:NDIL)=EFF_M(IG)%MATRIX(:NDIL,2)
        XSDIL(3,:NDIL)=SIGF(IG,:NDIL)*XSDIL(1,:NDIL)
        XSDIL(4,:NDIL)=EFF_M(IG)%MATRIX(:NDIL,3)
        DGAR(:NDIL)=XSDIL(4,:NDIL)/(SIGSD(IG,1,:NDIL)*XSDIL(1,:NDIL))
        DO IL=2,NL
         XSDIL(3+IL,:NDIL)=DGAR(:NDIL)*SIGSD(IG,IL,:NDIL)*XSDIL(1,:NDIL)
        ENDDO
        IF(NPART.EQ.3+NL) EXIT
        IOF=3+NL
        DO IL=1,NL
          IF(LSCAT(IL)) THEN
            DO JG=ISMIN(IL,IG),ISMAX(IL,IG)
              IOF=IOF+1
              XSDIL(IOF,:NDIL)=DGAR(:NDIL)*SCAT(JG,IG,IL,:NDIL)*
     1        XSDIL(1,:NDIL)
            ENDDO
          ENDIF
        ENDDO
        DO IED=1,NED
          IOF=IOF+1
          IF((HVECT(IED).EQ.'NINEL').OR.(HVECT(IED).EQ.'NELAS').OR.
     1       (HVECT(IED).EQ.'N2N').OR.(HVECT(IED).EQ.'N3N').OR.
     2       (HVECT(IED).EQ.'N4N').OR.(HVECT(IED).EQ.'NX').OR.
     3       (HVECT(IED).EQ.'STRD')) THEN
            XSDIL(IOF,:NDIL)=SADD(IG,IED,:NDIL)*XSDIL(1,:NDIL)
          ELSE
            XSDIL(IOF,:NDIL)=SADD(IG,IED,:NDIL)*XSDIL(1,:NDIL)
          ENDIF
        ENDDO
        DO IDEL=1,NDEL
          IOF=IOF+1
          XSDIL(IOF,:NDIL)=ZDEL(IG,IDEL,:NDIL)*XSDIL(1,:NDIL)
        ENDDO
        IF(IOF.NE.NPART) CALL XABORT('LIBRSE: INVALID NPART.')
        ALLOCATE(DDD(NPART,MI))
        DDD(:NPART,:MI)=MATMUL(XSDIL(:NPART,:NDIL),
     1           V_M(IG)%MATRIX(:NDIL,:MI))
        SIGP_M(IG)%MATRIX(:NPART,:MI)=MATMUL(DDD(:NPART,:MI),
     1           T_M(IG)%MATRIX(:MI,:MI))
        DEALLOCATE(DDD,XSDIL)
        IF(IMPX.GE.2) THEN
          DO IPART=1,4
            WRITE(6,'(/17H LIBRSE: REACTION,A5,21H BASE POINTS IN GROUP,
     1      I4,1H:)') HPART(IPART),IG
            WRITE(6,'(1X,1P,12E12.4)') SIGP_M(IG)%MATRIX(IPART,:MI)
          ENDDO
        ENDIF
      ENDDO
      DEALLOCATE(DGAR)
*----
*  SAVE INFORMATION IN IPLIB
*----
      NPOS=0
      DO IG=1,NGRP
        DO JG=1,IG
          IF(ASSOCIATED(SCAT_M(IG,JG)%MATRIX)) NPOS=NPOS+1
        ENDDO
      ENDDO
      CALL LCMSIX(IPLIB,'PT-TABLE',1)
      CALL LCMPUT(IPLIB,'NDEL',1,1,NDEL)
      JPLIB1=LCMLID(IPLIB,'GROUP-RSE',NGRP)
      JPLIB2=LCMLID(IPLIB,HNAMIS,NPOS) ! holds TSCAT_M information
      IPOS=0
      DO IG=1,NGRP
        MI=MRANK(IG)
        IF(MI.EQ.0) CYCLE
        NPART=SIZE(SIGP_M(IG)%MATRIX,1)
        NJJ(IG)=0
        DO JG=IG,1,-1
          IF(ASSOCIATED(SCAT_M(IG,JG)%MATRIX)) THEN
            MJ=MRANK(JG)
            IPOS=IPOS+1
            IF(IPOS.GT.NPOS) CALL XABORT('LIBRSE: NPOS OVERFLOW.')
            NJJ(IG)=NJJ(IG)+1
            CALL LCMPDL(JPLIB2,IPOS,MI*MJ,4,TSCAT_M(IG,JG)%MATRIX)
            DEALLOCATE(TSCAT_M(IG,JG)%MATRIX)
          ENDIF
        ENDDO
        KPLIB=LCMDIL(JPLIB1,IG)
        LGBIN=NFS(IG)
        CALL LCMPUT(KPLIB,'SIGT_V',MI,4,SIGT_V(IG)%VECTOR)
        CALL LCMPUT(KPLIB,'XI_V',MI,4,XI_V(IG)%VECTOR)
        CALL LCMPUT(KPLIB,'GAMMA_V',MI,4,GAMMA_V(IG)%VECTOR)
        CALL LCMPUT(KPLIB,'T_M',MI*MI,4,T_M(IG)%MATRIX)
        CALL LCMPUT(KPLIB,'U_M',LGBIN*MI,4,U_M(IG)%MATRIX)
        CALL LCMPUT(KPLIB,'RSE-TABLE',NPART*MI,4,SIGP_M(IG)%MATRIX)
        DEALLOCATE(SIGT_V(IG)%VECTOR,GAMMA_V(IG)%VECTOR,XI_V(IG)%VECTOR,
     1  SIGP_M(IG)%MATRIX)
        DO IL=1,NL
          ISM(1,IL)=ISMIN(IL,IG)
          ISM(2,IL)=ISMAX(IL,IG)
        ENDDO
        CALL LCMPUT(KPLIB,'ISM-LIMITS',2*NL,1,ISM)
      ENDDO
      CALL LCMPUT(IPLIB,'SVDEPS',1,2,SVDEPS)
      CALL LCMPUT(IPLIB,'NJJS00',NGRP,1,NJJ)
*----
*  PROCESS UNRESOLVED ENERGY DOMAIN.
*----
      JPLIB1=LCMLID(IPLIB,'GROUP-PT',NGRP)
      DO IG=1,NGRP
        IF(NFS(IG).GT.0) CYCLE
        NPART=3+NL+NED+NDEL
        DO IL=1,NL
          NPART=NPART+MAX(ISMAX(IL,IG)-ISMIN(IL,IG)+1,0)
        ENDDO
        ALLOCATE(SIGP(MAXNOR*NPART))
        CALL XDRSET(SIGP,MAXNOR*NPART,0.0)
        NDIL=NDIL-1
        CALL LIBTAB(IG,NGRP,NL,NDIL,NPART,NED,NDEL,HNAMIS,IMPX,LSCAT,
     1  LSIGF,LADD,DILUT,TOTAL,SIGF,SIGSD,SCAT,SADD,ZDEL,1.0,ISMIN,
     2  ISMAX,MRANK(IG),SIGP)
        NDIL=NDIL+1
*
        IF(MRANK(IG).GT.1) THEN
*          SAVE THE PROBABILITY TABLE INTO IPLIB.
           KPLIB=LCMDIL(JPLIB1,IG)
           CALL LCMPUT(KPLIB,'PROB-TABLE',MAXNOR*NPART,2,SIGP)
           DO IL=1,NL
             ISM(1,IL)=ISMIN(IL,IG)
             ISM(2,IL)=ISMAX(IL,IG)
           ENDDO
           CALL LCMPUT(KPLIB,'ISM-LIMITS',2*NL,1,ISM)
        ENDIF
        DEALLOCATE(SIGP)
      ENDDO
      CALL LCMPUT(IPLIB,'NOR',NGRP,1,MRANK)
      CALL LCMSIX(IPLIB,' ',2)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DO IG=1,NGRP
        IF(NFS(IG).EQ.0) CYCLE
        DEALLOCATE(U_M(IG)%MATRIX,U_V(IG)%VECTOR,UU_V(IG)%VECTOR)
        DEALLOCATE(PSI_V(IG)%VECTOR)
        DEALLOCATE(PHI_M(IG)%MATRIX,EFF_M(IG)%MATRIX)
        DEALLOCATE(T_M(IG)%MATRIX,V_M(IG)%MATRIX)
      ENDDO
      DEALLOCATE(ISMAX,ISMIN)
      DEALLOCATE(GOLD,LADD,LSCAT,DELTG,ZDEL,SADD,SCAT,SIGSD,SIGF,TOTAL,
     1 FLUX)
      DEALLOCATE(DILUT,ISM,SIGP_M,TSCAT_M,PSI_V,XI_V,GAMMA_V,SIGT_V,T_M,
     1 SCAT_M,SIGT_M,V_M,U_M,UU_V,U_V,PHI_M,EFF_M,V,W,STIS,PRI,SIGT,
     2 SIGS,STR,DEL,UUU,EBIN)
      DEALLOCATE(MRANK,NJJ)
      RETURN
      END
