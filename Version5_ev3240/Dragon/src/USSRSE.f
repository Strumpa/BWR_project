*DECK USSRSE
      SUBROUTINE USSRSE(IPTRK,IPLIB,IPLI0,IFTRAK,NREG,NUN,NBMIX,NBISO,
     1 NIRES,NL,NED,NDEL,ISONAM,ISOBIS,HCAL,MAT,VOL,KEYFLX,CDOOR,
     2 LEAKSW,IMPX,DEN,MIX,IAPT,IPHASE,NGRP,IGRMIN,IGRMAX,NBNRS,IREX,
     3 TITR,ICORR,MAXST,GOLD,UNGAR,PHGAR,STGAR,SFGAR,SSGAR,S0GAR,SAGAR,
     4 SDGAR,MASKG,SIGGAR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the self-shielded cross sections in each energy group using
* the resonance spectrum expansion method.
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
* IPTRK   pointer to the tracking (L_TRACK signature).
* IPLIB   pointer to the internal microscopic cross section library
*         with subgroups (L_LIBRARY signature).
* IPLI0   pointer to the internal microscopic cross section library
*         builded by the self-shielding module.
* IFTRAK  file unit number used to store the tracks.
* NREG    number of regions.
* NUN     number of unknowns per energy group and band.
* NBMIX   number of mixtures in the internal library.
* NBISO   number of isotopes.
* NIRES   number of correlated resonant isotopes.
* NL      number of legendre orders required in the calculation
*         (NL=1 or higher).
* NED     number of extra vector edits.
* NDEL    number of delayed neutron precursor groups.
* ISONAM  alias name of isotopes in IPLIB.
* ISOBIS  alias name of isotopes in IPLI0.
* HCAL    name of the self-shielding calculation.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  pointers of fluxes in unknown vector.
* CDOOR   name of the geometry/solution operator.
* LEAKSW  leakage flag (LEAKSW=.true. if neutron leakage through
*         external boundary is present).
* IMPX    print flag (equal to zero for no print).
* DEN     density of each isotope.
* MIX     mix number of each isotope (can be zero).
* IAPT    resonant isotope index associated with isotope I. Mixed
*         moderator if IAPT(I)=NIRES+1. Out-of-fuel isotope if
*         IAPT(I)=0.
* IPHASE  type of flux solution (=1 use a native flux solution door;
*         =2 use collision probabilities).
* NGRP    number of energy groups.
* IGRMIN  first group where the self-shielding is applied.
* IGRMAX  most thermal group where the self-shielding is applied.
* NBNRS   number of correlated fuel regions. Note that NBNRS=max(IREX).
* IREX    fuel region index assigned to each mixture. Equal to zero
*         in non-resonant mixtures or in mixtures not used.
* TITR    title.
* ICORR   mutual resonance shielding flag (=1 to suppress the model
*         in cases it is required in LIB operator).
* MAXST   maximum number of fixed point iterations for the ST scattering
*         source.
*
*Parameters: output
* GOLD    Goldstein-Cohen parameters.
* UNGAR   averaged fluxes per volume.
* PHGAR   averaged fluxes.
* STGAR   microscopic self-shielded total x-s.
* SFGAR   microscopic self-shielded fission x-s.
* SSGAR   microscopic self-shielded scattering x-s.
* S0GAR   microscopic transfer scattering xs (isotope,secondary,
*         primary).
* SAGAR   microscopic self-shielded additional xs.
* SDGAR   microscopic self-shielded delayed nu-sigf xs.
* MASKG   energy group mask pointing on self-shielded groups.
* SIGGAR  macroscopic x-s of the non-resonant isotopes in each mixture:
*         (*,*,*,1) total; (*,*,*,2) transport correction; 
*         (*,*,*,3) P0 scattering; (*,*,*,4) flux times P0 scattering.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPLIB,IPLI0
      INTEGER IFTRAK,NREG,NUN,NBMIX,NBISO,NIRES,NL,NED,NDEL,
     1 ISONAM(3,NBISO),ISOBIS(3,NBISO),MAT(NREG),KEYFLX(NREG),IMPX,
     2 MIX(NBISO),IAPT(NBISO),IPHASE,NGRP,IGRMIN,IGRMAX,NBNRS,
     3 IREX(NBMIX),ICORR,MAXST
      REAL VOL(NREG),DEN(NBISO),GOLD(NIRES,NGRP),UNGAR(NREG,NIRES,NGRP),
     1 PHGAR(NBNRS,NIRES,NGRP),STGAR(NBNRS,NIRES,NGRP),
     2 SFGAR(NBNRS,NIRES,NGRP),SSGAR(NBNRS,NIRES,NL,NGRP),
     3 S0GAR(NBNRS,NIRES,NL,NGRP,NGRP),SAGAR(NBNRS,NIRES,NED,NGRP),
     4 SDGAR(NBNRS,NIRES,NDEL,NGRP),SIGGAR(NBMIX,0:NIRES,NGRP,4)
      LOGICAL LEAKSW,MASKG(NGRP,NIRES)
      CHARACTER HCAL*12,CDOOR*12,TITR*72
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPP,KPLIB,LPLIB,MPLIB,JPLI0,KPLI0,IOFSET
      LOGICAL LLIB
      PARAMETER (MAXED=50,MAXNOR=20)
      CHARACTER TEXT12*12,HVECT(MAXED)*8,CBDPNM*12,HSMG*131
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPPT1,IPISO1,IPISO2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IWRK
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NOR,IPPT2,ISM
      REAL, ALLOCATABLE, DIMENSION(:) :: GAS,GA1,VOLMER,DELTAU,GOLD2
      REAL, ALLOCATABLE, DIMENSION(:,:) :: GA2,CONR,XFLUX
      TYPE VECTOR_ARRAY
        DOUBLE PRECISION, POINTER, DIMENSION(:) :: VECTOR
      END TYPE VECTOR_ARRAY
      TYPE(VECTOR_ARRAY), ALLOCATABLE, DIMENSION(:,:) :: GAMMA_V
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IPPT1(NIRES))
      ALLOCATE(NOR(NIRES,NGRP),IPPT2(NIRES,5),IWRK(NGRP),ISM(2,NL))
      ALLOCATE(GAS(NGRP),GA1(NGRP),GA2(NGRP,NGRP),CONR(NBNRS,NIRES),
     1 VOLMER(0:NBNRS),DELTAU(NGRP),GAMMA_V(NGRP,NIRES))
      ALLOCATE(IPISO1(NBISO),IPISO2(NBISO))
*
      CALL KDRCPU(TK1)
      CALL XDRSET(PHGAR,NBNRS*NIRES*NGRP,1.0)
      NOR(:NIRES,1)=-1
*
      IF(NED.GT.0) THEN
        IF(NED.GT.MAXED) CALL XABORT('USSRSE: INVALID VALUE OF MAXED.')
        CALL LCMGTC(IPLIB,'ADDXSNAME-P0',8,NED,HVECT)
      ENDIF
*
      CALL LIBIPS(IPLIB,NBISO,IPISO1)
      CALL LIBIPS(IPLI0,NBISO,IPISO2)
      CALL XDRSET(SIGGAR,NBMIX*(NIRES+1)*NGRP*4,0.0)
      DO 190 ISO=1,NBISO
      IBM=MIX(ISO)
      DO 30 I=1,NREG
      IF(MAT(I).EQ.IBM) GO TO 35
   30 CONTINUE
      GO TO 190
   35 IRES=IAPT(ISO)
      DENN=DEN(ISO)
      JRES=IRES
      IF(IRES.EQ.NIRES+1) JRES=0
*----
*  RECOVER INFINITE DILUTION OR SELF-SHIELDED CROSS SECTIONS AND
*  COMPUTE OUT-OF-FUEL MACROSCOPIC CROSS SECTIONS.
*----
      KPLI0=IPISO2(ISO) ! set ISO-th isotope
      IF(C_ASSOCIATED(KPLI0)) THEN
         CALL LCMLEN(KPLI0,'NTOT0',ILEN0,ITYLCM)
         IF(ILEN0.NE.0) THEN
            LLIB=.FALSE.
            IPP=KPLI0
         ELSE
            LLIB=.TRUE.
            IPP=IPISO1(ISO) ! set ISO-th isotope
         ENDIF
      ELSE
         LLIB=.TRUE.
         IPP=IPISO1(ISO) ! set ISO-th isotope
      ENDIF
      IF(LLIB.AND.(.NOT.C_ASSOCIATED(IPP))) THEN
         WRITE(HSMG,'(18H USSRSE: ISOTOPE '',3A4,7H'' (ISO=,I8,5H) IS ,
     1   39HNOT AVAILABLE IN THE ORIGINAL MICROLIB.)') (ISONAM(I0,ISO),
     2   I0=1,3),ISO
         CALL XABORT(HSMG)
      ENDIF
      IF((.NOT.LLIB).AND.(IMPX.GT.2)) WRITE(6,'(/18H USSRSE: RECOVER I,
     1 8HSOTOPE '',3A4,23H'' FROM THE NEW LIBRARY.)') (ISOBIS(I0,ISO),
     2 I0=1,3)
      IF((DENN.NE.0.0).AND.(IBM.NE.0)) THEN
         CALL LCMLEN(IPP,'NTOT0',ILENGT,ITYLCM)
         IF(ILENGT.NE.NGRP) THEN
           CALL LCMLIB(IPP)
           CALL XABORT('USSRSE: INVALID X-SECTIONS.')
         ENDIF
         CALL LCMGET(IPP,'NTOT0',GA1)
         CALL XDRLGS(IPP,-1,IMPX,0,0,1,NGRP,GAS,GA2,ITYPRO)
         DO 40 IGRP=1,NGRP
         SIGGAR(IBM,JRES,IGRP,1)=SIGGAR(IBM,JRES,IGRP,1)+DENN*GA1(IGRP)
   40    CONTINUE
         CALL LCMGET(IPP,'SIGS00',GA1)
         CALL LCMLEN(IPP,'NWT0',ILENGT,ITYLCM)
         IF(ILENGT.GT.0) THEN
            CALL LCMGET(IPP,'NWT0',GAS)
         ELSE
            CALL XDRSET(GAS,NGRP,1.0)
         ENDIF
         DO 45 IGRP=1,NGRP
         SIGGAR(IBM,JRES,IGRP,3)=SIGGAR(IBM,JRES,IGRP,3)+DENN*GA1(IGRP)
         DO 44 JGRP=1,IGRP
         SIGGAR(IBM,JRES,IGRP,4)=SIGGAR(IBM,JRES,IGRP,4)+DENN*
     1   GA2(IGRP,JGRP)*GAS(JGRP)
   44    CONTINUE
   45    CONTINUE
         CALL LCMLEN(IPP,'TRANC',ILENGT,ITYLCM)
         IF(ILENGT.GT.0) THEN
            CALL LCMGET(IPP,'TRANC',GA1)
         ELSE
            CALL XDRSET(GA1,NGRP,0.0)
         ENDIF
         DO 50 IGRP=1,NGRP
         SIGGAR(IBM,JRES,IGRP,2)=SIGGAR(IBM,JRES,IGRP,2)+DENN*GA1(IGRP)
   50    CONTINUE
      ENDIF
      CALL LCMGET(IPLI0,'DELTAU',DELTAU)
*----
*  RECOVER PROBABILITY TABLE INFORMATION.
*----
      IF((IRES.GT.0).AND.(IRES.LE.NIRES)) THEN
         IF(NOR(IRES,1).EQ.-1) THEN
            KPLIB=IPISO1(ISO) ! set ISO-th isotope
*
*           RECOVER INFINITE DILUTION VALUES.
            CALL LCMGET(KPLIB,'NTOT0',GAS)
            DO 55 IGRP=1,NGRP
            STGAR(:NBNRS,IRES,IGRP)=0.0
            STGAR(:NBNRS,IRES,IGRP)=GAS(IGRP)
            SFGAR(:NBNRS,IRES,IGRP)=0.0
            SAGAR(:NBNRS,IRES,:NED,IGRP)=0.0
            SDGAR(:NBNRS,IRES,:NDEL,IGRP)=0.0
   55       CONTINUE
            CALL LCMLEN(KPLIB,'NUSIGF',ILENGT,ITYLCM)
            IF(ILENGT.GT.0) THEN
               CALL LCMGET(KPLIB,'NUSIGF',GAS)
               DO 60 IGRP=1,NGRP
               SFGAR(:NBNRS,IRES,IGRP)=GAS(IGRP)
   60          CONTINUE
            ENDIF
            DO 80 IL=1,NL
            CALL XDRLGS(KPLIB,-1,IMPX,IL-1,IL-1,1,NGRP,GAS,GA2,ITYPRO)
*           JG IS THE SECONDARY GROUP.
            DO 72 IGRP=1,NGRP
            SSGAR(:NBNRS,IRES,IL,IGRP)=GAS(IGRP)
            DO 70 JGRP=1,NGRP
            S0GAR(:NBNRS,IRES,IL,JGRP,IGRP)=GA2(JGRP,IGRP)
   70       CONTINUE
   72       CONTINUE
   80       CONTINUE
            DO 110 IED=1,NED
            CALL LCMLEN(KPLIB,HVECT(IED),ILENGT,ITYLCM)
            IF(ILENGT.GT.0) THEN
               CALL LCMGET(KPLIB,HVECT(IED),GAS)
               DO 105 IGRP=1,NGRP
               SAGAR(:NBNRS,IRES,IED,IGRP)=GAS(IGRP)
  105          CONTINUE
            ENDIF
  110       CONTINUE
            DO 130 IDEL=1,NDEL
            WRITE(TEXT12,'(6HNUSIGF,I2.2)') IDEL
            CALL LCMLEN(KPLIB,TEXT12,ILENGT,ITYLCM)
            IF(ILENGT.GT.0) THEN
               CALL LCMGET(KPLIB,TEXT12,GAS)
               DO 125 IGRP=1,NGRP
               SDGAR(:NBNRS,IRES,IDEL,IGRP)=GAS(IGRP)
  125          CONTINUE
            ENDIF
  130       CONTINUE
*
            GOLD(IRES,:NGRP)=1.0
            NOR(IRES,:NGRP)=0
            CALL LCMLEN(KPLIB,'NGOLD',ILENGT,ITYLCM)
            IF(ILENGT.GT.0) THEN
              ALLOCATE(GOLD2(NGRP))
              CALL LCMGET(KPLIB,'NGOLD',GOLD2)
              GOLD(IRES,IGRMIN:IGRMAX)=GOLD2(IGRMIN:IGRMAX)
              DEALLOCATE(GOLD2)
            ENDIF
            CALL LCMLEN(KPLIB,'PT-TABLE',ILENGT,ITYLCM)
            IF(ILENGT.EQ.-1) THEN
              CALL LCMSIX(KPLIB,'PT-TABLE',1)
                CALL LCMGET(KPLIB,'NOR',IWRK)
                CALL LCMLEN(KPLIB,'GROUP-RSE',ILENGT,ITYLCM)
                IF(ILENGT.GT.0) THEN
                  LPLIB=LCMGID(KPLIB,'GROUP-RSE')
                  DO IGRP=IGRMIN,IGRMAX
                    MI=NOR(IRES,IGRP)
                    CALL LCMLEL(LPLIB,IGRP,ILENGT,ITYLCM)
                    IF(ILENGT.EQ.0) GOLD(IRES,IGRP)=1.0
                    NOR(IRES,IGRP)=IWRK(IGRP)
                    CALL LCMLEL(LPLIB,IGRP,ILENGT,ITYLCM)
                    IF(GOLD(IRES,IGRP).EQ.-1001.0) THEN
                      MPLIB=LCMGIL(LPLIB,IGRP)
                      CALL LCMGPD(MPLIB,'GAMMA_V',IOFSET)
                      CALL C_F_POINTER(IOFSET,GAMMA_V(IGRP,IRES)%VECTOR,
     1                (/MI/))
                    ENDIF
                  ENDDO
                ENDIF
              CALL LCMSIX(KPLIB,' ',2)
            ENDIF
         ENDIF
      ENDIF
  190 CONTINUE
      CALL KDRCPU(TK2)
      IF(IMPX.GT.1) WRITE(6,'(/34H USSRSE: CPU TIME SPENT TO RECOVER,
     1 23H INFINITE-DILUTION XS =,F8.1,8H SECOND./)') TK2-TK1
*
      CALL KDRCPU(TK1)
      TK4=0.0
      TK5=0.0
      ICPIJ=0
*----
*  COMPUTE THE MERGED VOLUMES AND NUMBER DENSITIES.
*----
      VOLMER(0:NBNRS)=0.0
      DO 210 I=1,NREG
      IBM=MAT(I)
      IF(IBM.GT.0) VOLMER(IREX(IBM))=VOLMER(IREX(IBM))+VOL(I)
  210 CONTINUE
      CONR(:NBNRS,:NIRES)=0.0
      DO 240 ISO=1,NBISO
      JRES=IAPT(ISO)
      IF((JRES.GT.0).AND.(JRES.LE.NIRES)) THEN
         DENN=DEN(ISO)
         DO 230 IREG=1,NREG
         IBM=MAT(IREG)
         IF(MIX(ISO).EQ.IBM) THEN
            IND=IREX(IBM)
            IF(IND.EQ.0) CALL XABORT('USSRSE: IREX FAILURE.')
            CONR(IND,JRES)=CONR(IND,JRES)+DENN*VOL(IREG)/VOLMER(IND)
         ENDIF
  230    CONTINUE
      ENDIF
  240 CONTINUE
*----
*  RECOVER POSITION OF PROBABILITY TABLES AND NAME OF RESONANT ISOTOPE.
*----
      DO 270 IRES=1,NIRES
      ISOT=0
      DO 250 JSOT=1,NBISO
      IF(IAPT(JSOT).EQ.IRES) THEN
         ISOT=JSOT
         GO TO 260
      ENDIF
  250 CONTINUE
      CALL XABORT('USSRSE: UNABLE TO FIND A RESONANT ISOTOPE.')
  260 KPLIB=IPISO1(ISOT) ! set ISOT-th isotope
      CALL LCMLEN(KPLIB,'PT-TABLE',ILONG,ITYLCM)
      IF(ILONG.EQ.0) CALL XABORT('USSRSE: BUG1.')
      CALL LCMSIX(KPLIB,'PT-TABLE',1)
      CALL LCMGET(KPLIB,'NDEL',NDEL0)
      IF(NDEL0.GT.NDEL) CALL XABORT('USSRSE: NDEL OVERFLOW.')
      IPPT1(IRES)=KPLIB
      CALL LCMSIX(KPLIB,' ',2)
      IPPT2(IRES,1)=IREX(MIX(ISOT))
      IPPT2(IRES,2)=ISONAM(1,ISOT)
      IPPT2(IRES,3)=ISONAM(2,ISOT)
      IPPT2(IRES,4)=ISONAM(3,ISOT)
      IPPT2(IRES,5)=NDEL0
      IF(IPPT2(IRES,1).LE.0) CALL XABORT('USSRSE: BUG3.')
  270 CONTINUE
*----
*  DETERMINE WHICH GROUPS ARE SELF-SHIELDED.
*----
      DO 290 IGRP=1,NGRP
      DO 280 IRES=1,NIRES
      MASKG(IGRP,IRES)=((IGRP.GE.IGRMIN).AND.(IGRP.LE.IGRMAX).AND.
     1 (NOR(IRES,IGRP).GT.0))
  280 CONTINUE
  290 CONTINUE
*----
*  INITIALIZATION OF THE MULTIBAND FLUXES AND SOURCES.
*----
      CALL LCMSIX(IPLI0,'SHIBA_SG',1)
      CALL LCMSIX(IPLI0,HCAL,1)
      DO 310 IRES=1,NIRES
      WRITE(CBDPNM,'(3HCOR,I4.4,1H/,I4.4)') IRES,NIRES
      CALL LCMSIX(IPLI0,CBDPNM,1)
      JPLI0=LCMLID(IPLI0,'NWT0-PT',NGRP)
      DO 300 IGRP=1,NGRP
      IF(MASKG(IGRP,IRES)) THEN
         CALL LCMLEL(JPLI0,IGRP,ILENGT,ITYLCM)
         IF(ILENGT.EQ.0) THEN
            MI=NOR(IRES,IGRP)
            ALLOCATE(XFLUX(NBNRS,MI))
            IF(GOLD(IRES,IGRP).EQ.1.0) THEN
              XFLUX(:NBNRS,:MI)=1.0
            ELSE IF(GOLD(IRES,IGRP).EQ.-1001.0) THEN
              DO IM=1,MI
                XFLUX(:NBNRS,IM)=REAL(GAMMA_V(IGRP,IRES)%VECTOR(IM))
              ENDDO
            ENDIF
            CALL LCMPDL(JPLI0,IGRP,NBNRS*MI,2,XFLUX)
            DEALLOCATE(XFLUX)
         ENDIF
      ENDIF
  300 CONTINUE
      CALL LCMSIX(IPLI0,' ',2)
  310 CONTINUE
*
      DO 340 IRES=1,NIRES
      DO 330 IGRP=1,NGRP
      IF(MASKG(IGRP,IRES)) ICPIJ=ICPIJ+NOR(IRES,IGRP)
  330 CONTINUE
      CALL KDRCPU(TKA)
*----
*  ITERATIVE APPROACH FOR THE HELIOS/WIMS-7 METHOD.
*----
      MAX_R=12
      CALL USSIT1(MAX_R,NGRP,MASKG(1,IRES),IRES,IPLI0,IPTRK,IFTRAK,
     1 CDOOR,IMPX,NBMIX,NREG,NUN,NL,IPHASE,MAXST,MAT,VOL,KEYFLX,LEAKSW,
     2 IREX,SIGGAR,TITR,NIRES,NBNRS,NOR,CONR,GOLD,IPPT1,IPPT2,STGAR,
     3 SSGAR,VOLMER,UNGAR)
*----
*  ITERATIVE APPROACH FOR THE RESONANCE SPECTRUM EXPANSION METHOD.
*----
      CALL USSIT3(MAXNOR,NGRP,MASKG(1,IRES),IRES,IPLI0,IPTRK,IFTRAK,
     1 CDOOR,IMPX,NBMIX,NREG,NUN,IPHASE,MAXST,MAT,VOL,KEYFLX,LEAKSW,
     2 IREX,SIGGAR,TITR,ICORR,NIRES,NBNRS,CONR,GOLD,IPPT1,IPPT2,
     3 VOLMER,DELTAU,UNGAR)
      CALL KDRCPU(TKB)
      TK4=TK4+(TKB-TKA)
  340 CONTINUE
*----
*  COMPUTE THE SELF-SHIELDED REACTION RATES.
*----
      CALL USSIT4(MAXNOR,IPLI0,IPPT1,IPPT2,NGRP,NIRES,NBNRS,NL,NED,
     1 NDEL,PHGAR,STGAR,SFGAR,SSGAR,S0GAR,SAGAR,SDGAR)
      CALL LCMSIX(IPLI0,' ',2)
      CALL LCMSIX(IPLI0,' ',2)
      CALL LCMVAL(IPLI0,' ')
*----
*  RESET MASKG FOR SPH CALCULATION IN SMALL LETHARGY WIDTH GROUPS.
*----
      IF(NIRES.GT.1) THEN
        DO 360 IGRP=1,NGRP
        DO 350 IRES=1,NIRES
        IF(MASKG(IGRP,IRES)) THEN
          MASKG(IGRP,IRES)=.NOT.(GOLD(IRES,IGRP).EQ.-1001.)
          IF(DELTAU(IGRP).GT.0.1) MASKG(IGRP,IRES)=.TRUE.
        ENDIF
  350   CONTINUE
  360   CONTINUE
      ENDIF
      CALL KDRCPU(TK2)
      IF(IMPX.GT.1) WRITE(6,'(/34H USSRSE: CPU TIME SPENT TO COMPUTE,
     1 31H SELF-SHIELDED REACTION RATES =,F8.1,19H SECOND, INCLUDING:
     2 /9X,F8.1,46H SECOND TO BUILD/SOLVE SUBGROUP MATRIX SYSTEM;/9X,
     4 F8.1,38H SECOND TO COMPUTE THE REACTION RATES./9X,9HNUMBER OF,
     5 23H ASSEMBLY DOORS CALLS =,I5,1H.)') TK2-TK1,TK4,TK5,ICPIJ
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IPISO2,IPISO1)
      DEALLOCATE(GAMMA_V,DELTAU,VOLMER,CONR,GA2,GA1,GAS)
      DEALLOCATE(ISM,IWRK,IPPT2,NOR)
      DEALLOCATE(IPPT1)
      RETURN
      END
