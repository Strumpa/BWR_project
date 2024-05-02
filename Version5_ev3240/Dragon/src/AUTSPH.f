*DECK AUTSPH
      SUBROUTINE AUTSPH(IPLI0,IPTRK,IFTRAK,NREG,NUN,NBMIX,NBISO,NIRES,
     1 NL,NED,NDEL,HCAL,MAT,VOL,KEYFLX,CDOOR,LEAKSW,IMPX,DEN,MIX,IAPT,
     2 ITRANC,IPHASE,NGRP,MASKG,IREX,TITR,SIGGAR,UNGAR,PHGAR,STGAR,
     3 SFGAR,SSGAR,S0GAR,SAGAR,SDGAR,DELTAU,SPH)
*
*-----------------------------------------------------------------------
*
*Purpose:
* SPH equivalence procedure over the self-shielded cross sections. Use
* all the standard solution doors of Dragon. Autosecol specific version.
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
* IPLI0   pointer to the LCM object containing subgroup-related
*         information.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  file unit number used to store the tracks.
* NREG    number of regions.
* NUN     number of unknowns per energy group.
* NBMIX   number of mixtures in the internal library.
* NBISO   number of isotopes.
* NIRES   number of correlated resonant isotopes.
* NL      number of Legendre orders required in the calculation
*         (NL=1 or higher).
* NED     number of extra vector edits.
* NDEL    number of delayed neutron precursor groups.
* HCAL    name of the self-shielding calculation.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  pointers of fluxes in unknown vector.
* CDOOR   name of the geometry/solution operator.
* LEAKSW  leakage flag (LEAKSW=.TRUE. if neutron leakage through
*         external boundary is present).
* IMPX    print flag (equal to zero for no print).
* DEN     density of each isotope.
* MIX     mix number of each isotope (can be zero).
* IAPT    resonant isotope index associated with isotope I. Mixed
*         moderator if IAPT(I)=NIRES+1. Out-of-fuel isotope if
*         IAPT(I)=0.
* ITRANC  type of transport correction.
* IPHASE  type of flux solution (=1 use a native flux solution door;
*         =2 use collision probabilities).
* NGRP    number of energy groups.
* MASKG   energy group mask pointing on self-shielded groups.
* IREX    fuel region index assigned to each mixture. Equal to zero
*         in non-resonant mixtures or in mixtures not used.
* TITR    title.
* SIGGAR  macroscopic x-s of the non-resonant isotopes in each mixture:
*         (*,*,*,1) total; (*,*,*,2) transport correction; 
*         (*,*,*,3) P0 scattering.
* UNGAR   averaged fluxes per volume.
* STGAR   microscopic self-shielded total x-s.
*
*Parameters: input/output
* PHGAR   uncorrected and SPH-corrected averaged fluxes.
* SFGAR   uncorrected and SPH-corrected microscopic self-shielded 
*         fission x-s.
* SSGAR   uncorrected and SPH-corrected microscopic 
*         self-shielded scattering x-s.
* S0GAR   uncorrected and SPH-corrected microscopic 
*         transfer scattering x-s
*         (isotope,secondary,primary).
* SAGAR   uncorrected and SPH-corrected microscopic 
*         additional x-s.
* SDGAR   uncorrected and SPH-corrected microscopic 
*         self-shielded delayed nu-sigf x-s.
* DELTAU  lethargy width of each energy group.
*
*Parameters: output
* SPH     SPH factors.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLI0,IPTRK
      INTEGER IFTRAK,NREG,NUN,NBMIX,NBISO,NIRES,NL,NED,NDEL,MAT(NREG),
     1 KEYFLX(NREG),IMPX,MIX(NBISO),IAPT(NBISO),ITRANC,IPHASE,NGRP,
     2 IREX(NBMIX)
      REAL VOL(NREG),DEN(NBISO),SIGGAR(NBMIX,0:NIRES,NGRP,3),
     1 UNGAR(NREG,NGRP),PHGAR(NIRES,NGRP),STGAR(NIRES,NGRP),
     2 SFGAR(NIRES,NGRP),SSGAR(NIRES,NL,NGRP),S0GAR(NIRES,NL,NGRP,NGRP),
     3 SAGAR(NIRES,NED,NGRP),SDGAR(NIRES,NDEL,NGRP),DELTAU(NGRP),
     4 SPH(NIRES,NGRP)
      LOGICAL LEAKSW,MASKG(NGRP)
      CHARACTER CDOOR*12,HCAL*12,TITR*72
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPLI0,KPLI0,IPMACR,IPSOU
      LOGICAL LHOMOG,LPROB,LTIT,LEXAC,REBFLG
      INTEGER NALBP
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NPSYS
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGTXS,SIGS0X
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SUNKNO,FUNKNO,SIGTI
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LVOL
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NPSYS(NGRP))
      ALLOCATE(SIGTI(NBMIX,5),SIGTXS(0:NBMIX),SIGS0X(0:NBMIX),
     1 SUNKNO(NUN,NGRP),FUNKNO(NUN,NGRP),LVOL(NREG))
*----
*  SET LHOMOG.
*----
      NALBP=0
      LHOMOG=.TRUE.
      DO 10 I=1,NREG
      IBM=MAT(I)
      IF(IBM.EQ.0) GO TO 10
      IF(IREX(IBM).EQ.0) LHOMOG=.FALSE.
   10 CONTINUE
      IF(LHOMOG) GO TO 260
*----
* SET MACRO CALCULATION
*----
      ICPIJ=0
      CALL KDRCPU(TK1)
      CALL LCMSIX(IPLI0,'SHIBA_SG',1)
      CALL LCMSIX(IPLI0,HCAL,1)
      LTIT=.TRUE.
      JPLI0=LCMLID(IPLI0,'GROUP',NGRP)
*----
*  LOOP OVER SELF-SHIELDED ENERGY GROUPS.
*----
      CALL XDRSET(FUNKNO,NUN*NGRP,0.0)
      CALL XDRSET(SUNKNO,NUN*NGRP,0.0)
      CALL XDISET(NPSYS,NGRP,0)
      DO 80 IGRP=1,NGRP
      IF(.NOT.MASKG(IGRP)) GO TO 80
      NPSYS(IGRP)=IGRP
*----
*  SET THE MIXTURE-DEPENDENT MACROSCOPIC XS.
*----
      CALL XDRSET(SIGTI,NBMIX*5,0.0)
      DO 50 IBM=1,NBMIX
      DO 40 IRES=0,NIRES
      IF(IRES.EQ.0) THEN
         SIGTI(IBM,1)=SIGTI(IBM,1)+SIGGAR(IBM,0,IGRP,1)
         SIGTI(IBM,3)=SIGTI(IBM,3)+SIGGAR(IBM,0,IGRP,3)
         IF(ITRANC.NE.0) SIGTI(IBM,2)=SIGTI(IBM,2)+
     1   SIGGAR(IBM,0,IGRP,2)
      ELSE IF((IRES.GT.0).AND.(IREX(IBM).NE.0)) THEN
         DENN=0.0
         DO 20 ISO=1,NBISO
         IF((IAPT(ISO).EQ.IRES).AND.(MIX(ISO).EQ.IBM)) DENN=DEN(ISO)
   20    CONTINUE
         SIGTI(IBM,5)=SIGTI(IBM,5)+DENN*STGAR(IRES,IGRP)
         DO 30 JGRP=1,NGRP
           SIGTI(IBM,4)=SIGTI(IBM,4)+PHGAR(IRES,JGRP)*DENN*
     1     S0GAR(IRES,1,IGRP,JGRP)*DELTAU(JGRP)/DELTAU(IGRP)
   30    CONTINUE
      ENDIF
   40 CONTINUE
   50 CONTINUE
*----
*  COMPUTE THE SOURCES.
*----
      DO 60 I=1,NREG
      IBM=MAT(I)
      IF(IBM.EQ.0) GO TO 60
      SUNKNO(KEYFLX(I),IGRP)=SIGTI(IBM,3)
      IF(IREX(IBM).GT.0) THEN
         SUNKNO(KEYFLX(I),IGRP)=SUNKNO(KEYFLX(I),IGRP)+SIGTI(IBM,4)-
     1   UNGAR(I,IGRP)*SIGTI(IBM,5)
         IF(.NOT.LHOMOG) SUNKNO(KEYFLX(I),IGRP)=SUNKNO(KEYFLX(I),IGRP)-
     1   UNGAR(I,IGRP)*SIGTI(IBM,1)
      ENDIF
   60 CONTINUE
*
      IF(NPSYS(IGRP).NE.0) THEN
        ICPIJ=ICPIJ+1
        SIGTXS(0:)=0.0
        SIGS0X(0:)=0.0
        DO 70 IBM=1,NBMIX
        SIGTXS(IBM)=SIGTI(IBM,1)-SIGTI(IBM,2)
        IND=IREX(IBM)
        IF(IND.EQ.0) THEN
*          REMOVE TRANSPORT CORRECTION.
           SIGS0X(IBM)=-SIGTI(IBM,2)
        ELSE IF(IND.GT.0) THEN
*          BELL ACCELERATION.
           SIGTXS(IBM)=SIGTXS(IBM)+SIGTI(IBM,5)
           SIGS0X(IBM)=SIGTXS(IBM)
           IF(LHOMOG) SIGS0X(IBM)=SIGS0X(IBM)-SIGTI(IBM,1)
        ENDIF
   70   CONTINUE
        KPLI0=LCMDIL(JPLI0,IGRP)
        CALL LCMPUT(KPLI0,'DRAGON-TXSC',NBMIX+1,2,SIGTXS)
        CALL LCMPUT(KPLI0,'DRAGON-S0XSC',NBMIX+1,2,SIGS0X)
      ENDIF
   80 CONTINUE
*----
*  SOLVE FOR THE FLUX USING DIRECT SELF-SHIELDED CROSS SECTIONS
*----
      ISTRM=1
      NANI=1
      NW=0
      KNORM=1
      IMPY=MAX(0,IMPX-3)
      IF(IPHASE.EQ.1) THEN
*        USE A NATIVE DOOR.
         CALL DOORAV(CDOOR,JPLI0,NPSYS,IPTRK,IFTRAK,IMPY,NGRP,NREG,
     1   NBMIX,NANI,NW,MAT,VOL,KNORM,LEAKSW,TITR,NALBP,ISTRM)
      ELSE IF(IPHASE.EQ.2) THEN
*        USE A COLLISION PROBABILITY DOOR.
         IPIJK=1
         CALL DOORPV(CDOOR,JPLI0,NPSYS,IPTRK,IFTRAK,IMPY,NGRP,NREG,
     1   NBMIX,NANI,MAT,VOL,KNORM,IPIJK,LEAKSW,.FALSE.,TITR,NALBP)
      ENDIF
      IDIR=0
      LEXAC=.FALSE.
      IPMACR=C_NULL_PTR
      IPSOU=C_NULL_PTR
      REBFLG=.FALSE.
      CALL DOORFV(CDOOR,JPLI0,NPSYS,IPTRK,IFTRAK,IMPX,NGRP,NBMIX,IDIR,
     1 NREG,NUN,IPHASE,LEXAC,MAT,VOL,KEYFLX,TITR,SUNKNO,FUNKNO,IPMACR,
     2 IPSOU,REBFLG)
*----
*  LOOP OVER THE RESONANT ISOTOPES.
*----
      CALL XDLSET(LVOL,NREG,.FALSE.)
      CALL XDRSET(SPH,NGRP*NIRES,1.0)
      DO 240 IRES=1,NIRES
*----
*  HOMOGENIZE THE FLUX
*----
      VOLMER=0.0
      DO 100 I=1,NREG
      IBM=MAT(I)
      IF(IBM.EQ.0) GO TO 100
      DO 90 ISO=1,NBISO
      IF((IAPT(ISO).EQ.IRES).AND.(MIX(ISO).EQ.IBM)) LVOL(I)=.TRUE.
   90 CONTINUE
      IF(LVOL(I)) VOLMER=VOLMER+VOL(I)
  100 CONTINUE
      DO 230 IGRP=1,NGRP
      IF(NPSYS(IGRP).NE.0) THEN
        FLNEW=0.0
        DO 110 I=1,NREG
        IF(LVOL(I)) FLNEW=FLNEW+FUNKNO(KEYFLX(I),IGRP)*VOL(I)
  110   CONTINUE
        FLNEW=FLNEW/VOLMER
*----
*  SPH FACTOR CONTROL.
*----
        SPHNEW=PHGAR(IRES,IGRP)/FLNEW
        LPROB=(SPHNEW.LE.0.0).OR.(SPHNEW.GT.1.0).OR.(FLNEW.LT.0.05)
        IF(LPROB) SPHNEW=1.0
        SPH(IRES,IGRP)=SPHNEW
      ENDIF
      IF(MASKG(IGRP)) THEN
        SPHNEW=SPH(IRES,IGRP)
        PHGAR(IRES,IGRP)=PHGAR(IRES,IGRP)/SPHNEW
        SFGAR(IRES,IGRP)=SFGAR(IRES,IGRP)*SPHNEW
        DO 170 IL=1,NL
        IF(MOD(IL-1,2).EQ.0) THEN
          SSGAR(IRES,IL,IGRP)=SSGAR(IRES,IL,IGRP)*SPHNEW+
     1    STGAR(IRES,IGRP)*(1.0-SPHNEW)
        ELSE
          SSGAR(IRES,IL,IGRP)=0.0
        ENDIF
        DO 160 JGRP=1,NGRP
        IF(MOD(IL-1,2).EQ.0) THEN
          IF(IGRP.EQ.JGRP) THEN
            S0GAR(IRES,IL,IGRP,IGRP)=S0GAR(IRES,IL,IGRP,IGRP)*
     1      SPHNEW+STGAR(IRES,IGRP)*(1.0-SPHNEW)
          ELSE
            S0GAR(IRES,IL,JGRP,IGRP)=S0GAR(IRES,IL,JGRP,IGRP)*SPHNEW
          ENDIF
        ELSE
          IF(IGRP.EQ.JGRP) THEN
            S0GAR(IRES,IL,IGRP,IGRP)=S0GAR(IRES,IL,IGRP,IGRP)/
     1      SPHNEW+STGAR(IRES,IGRP)*(1.0-1.0/SPHNEW)
          ELSE
            S0GAR(IRES,IL,JGRP,IGRP)=S0GAR(IRES,IL,JGRP,IGRP)/
     1      SPH(IRES,JGRP)
          ENDIF
        ENDIF
        IF(MOD(IL-1,2).EQ.1) THEN
          SSGAR(IRES,IL,IGRP)=SSGAR(IRES,IL,IGRP)+
     1    S0GAR(IRES,IL,JGRP,IGRP)
        ENDIF
  160   CONTINUE
  170   CONTINUE
        DO 180 IED=1,NED
        SAGAR(IRES,IED,IGRP)=SAGAR(IRES,IED,IGRP)*SPHNEW
  180   CONTINUE
        DO 190 IDEL=1,NDEL
        SDGAR(IRES,IDEL,IGRP)=SDGAR(IRES,IDEL,IGRP)*SPHNEW
  190   CONTINUE
      ENDIF
  230 CONTINUE
  240 CONTINUE
*     ***************************************************************
      CALL LCMSIX(IPLI0,' ',2)
      CALL LCMSIX(IPLI0,' ',2)
      CALL KDRCPU(TK2)
      IF(IMPX.GT.1) WRITE(6,'(/34H AUTSPH: CPU TIME SPENT TO COMPUTE,
     1 18H THE SPH FACTORS =,F8.1,8H SECOND./9X,17HNUMBER OF ASSEMBL,
     2 15HY DOORS CALLS =,I5,1H.)') TK2-TK1,ICPIJ
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  260 DEALLOCATE(LVOL,FUNKNO,SUNKNO,SIGS0X,SIGTXS,SIGTI)
      DEALLOCATE(NPSYS)
      RETURN
      END
