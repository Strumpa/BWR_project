*DECK AUTFLU
      SUBROUTINE AUTFLU(IPTRK,IPLIB,IPLI0,IFTRAK,NREG,NUN,NBMIX,NBISO,
     1 NIRES,MAT,VOL,KEYFLX,CDOOR,LEAKSW,IMPX,DEN,MIX,IAPT,IPHASE,NGRP,
     2 IGRMIN,IGRRES,IGRMAX,DIL,TITR,IALTER,DELI,LBIN,NBIN,EBIN,MAXTRA,
     3 ISEED,ITRANC,UUU,PHI,SIGT,SIGS,SIGS1,SIGF,SIGGAR,MASKG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the flux in the Autolib fine groups using the Autosecol
* method.
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
* IGRRES  first resolved energy group.
* IGRMAX  most thermal group where the self-shielding is applied.
* DIL     microscopic dilution cross section of each isotope.
* TITR    title.
* IALTER  type of elastic slowing-down kernel (=0: use exact kernel;
*         =1: use an approximate kernel for the resonant isotopes).
* DELI    elementary lethargy width used by the elastic kernel.
* LBIN    total number of fine energy groups in the Autolib.
* NBIN    number of fine energy groups in each coarse energy group.
* EBIN    energy limits of the Autolib fine groups.
* MAXTRA  maximum number of down-scattering terms.
* ISEED   the seed for the generation of random numbers in the
*         unresolved energy domain.
* ITRANC  type of transport correction.
*
*Parameters: output
* UUU     lethargy limits of the Autolib fine groups.
* PHI     flux in the Autolib fine groups.
* SIGT    total microscopic x-s.
* SIGS    P0 scattering microscopic x-s.
* SIGS1   P1 scattering microscopic x-s.
* SIGF  nu*fission microscopic x-s.
* SIGGAR  macroscopic x-s of the non-resonant isotopes in each mixture:
*         (*,*,*,1) total; (*,*,*,2) transport correction; 
*         (*,*,*,3) P0 scattering.
* MASKG   energy group mask pointing on self-shielded groups.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPLIB,IPLI0
      INTEGER IFTRAK,NREG,NUN,NBMIX,NBISO,NIRES,MAT(NREG),KEYFLX(NREG),
     1 IMPX,MIX(NBISO),IAPT(NBISO),IPHASE,NGRP,IGRMIN,IGRRES,IGRMAX,
     2 IALTER,LBIN,NBIN(NGRP),MAXTRA,ISEED,ITRANC
      REAL VOL(NREG),DEN(NBISO),DIL(NBISO),DELI,EBIN(LBIN+1),
     1 UUU(LBIN+1),PHI(LBIN,NREG),SIGT(LBIN,NBISO),SIGS(LBIN,NBISO),
     2 SIGS1(LBIN,NBISO),SIGF(LBIN,NBISO),SIGGAR(NBMIX,0:NIRES,NGRP,3)
      LOGICAL LEAKSW,MASKG(NGRP)
      CHARACTER CDOOR*12,TITR*72
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPP,KPLIB,KPLI0,IPSYS
      LOGICAL LLIB
      DOUBLE PRECISION DUUU
      CHARACTER HNAMIS*12,HNABIS*12,HSMG*131
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO1,IPISO2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NJJ,IJJ,NEXT,III
      REAL, ALLOCATABLE, DIMENSION(:) :: GA1,DELTAU,SGAR,PRI
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIGINF,SCAT,GA2,CONC
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASKH
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(GA1(NGRP),GA2(NGRP,NGRP),CONC(NBMIX,NBISO),DELTAU(NGRP),
     1 MASKH(NGRP),SIGINF(NGRP,3))
      ALLOCATE(IPISO1(NBISO),IPISO2(NBISO))
      ALLOCATE(NJJ(NGRP),IJJ(NGRP),SGAR(NGRP**2),SCAT(NGRP,NGRP))
      ALLOCATE(PRI(MAXTRA),NEXT(NBISO),III(NBISO+1))
      NEXT(:NBISO)=0
*
      CALL KDRCPU(TK1)
      SIGGAR(:NBMIX,0:NIRES,:NGRP,:3)=0.0
      CALL LIBIPS(IPLIB,NBISO,IPISO1)
      CALL LIBIPS(IPLI0,NBISO,IPISO2)
      CALL LCMGET(IPLI0,'DELTAU',DELTAU)
      III(1)=1
      DO 120 ISO=1,NBISO
      SIGT(:LBIN,ISO)=0.0
      SIGS(:LBIN,ISO)=0.0
      SIGS1(:LBIN,ISO)=0.0
      SIGF(:LBIN,ISO)=0.0
      IF(IMPX.GT.1) WRITE(6,'(/32H AUTFLU: RECOVER XS FOR ISOTOPE=,I5)')
     1 ISO
      IBM=MIX(ISO)
      DO 10 NRE=1,NREG
      IF(MAT(NRE).EQ.IBM) GO TO 20
   10 CONTINUE
      III(ISO+1)=III(ISO)
      GO TO 120
   20 KPLIB=IPISO1(ISO) ! infinite dilution isotope
      KPLI0=IPISO2(ISO) ! self-shielded isotope
      CALL LCMGTC(KPLIB,'ALIAS',12,1,HNAMIS)
      IRES=IAPT(ISO)
      JRES=IRES
      IF(IRES.EQ.NIRES+1) JRES=0
      DENN=DEN(ISO)
      CALL LCMGET(KPLIB,'AWR',AWR)
      CALL LCMLEN(KPLIB,'BIN-NFS',LENGT,ITYLCM)
      IF((IRES.GT.0).AND.(IRES.LE.NIRES).AND.(LENGT.GT.0)) THEN
        ! resonant isotope
        IF(IMPX.GT.2) WRITE(6,'(26H AUTFLU: PROCESS AUTOLIB '',A12,
     1  1H'')') HNAMIS
        IF(LBIN.EQ.0) CALL XABORT('AUTFLU: MISSING AUTOLIB DATA.')
        DUUU=0.0D0
        UUU(1)=0.0
        DO 40 IGR=1,LBIN
        DUUU=DUUU+LOG(EBIN(IGR)/EBIN(IGR+1))
        UUU(IGR+1)=REAL(DUUU)
   40   CONTINUE
        ! recover unresolved xs values
        LLL=0
        IF(IGRRES.GT.IGRMIN) THEN
          CALL LCMLEN(KPLIB,'NUSIGF',N10,ITYLCM)
          CALL LCMLEN(KPLIB,'SIGS00',N12,ITYLCM)
          SIGINF(:NGRP,:3)=0.0
          CALL LCMGET(KPLIB,'NTOT0',SIGINF(1,1))
          IF(N10.GT.0) CALL LCMGET(KPLIB,'NUSIGF',SIGINF(1,2))
          IF(N12.GT.0) CALL LCMGET(KPLIB,'SIGS00',SIGINF(1,3))
          CALL AUTTAB(KPLIB,HNAMIS,IGRMIN,IGRRES,NGRP,LBIN,NBIN,UUU,
     1    ISEED,SIGINF,LLL,SIGT(1,ISO),SIGS(1,ISO),SIGF(1,ISO))
        ENDIF
        ! recover resolved xs values
        IF(LLL.LT.LBIN) THEN
          CALL LCMLEN(KPLIB,'BIN-NTOT0',LENG,ITYLCM)
          IF(LENG.GT.LBIN) CALL XABORT('AUTFLU: LBIN OVERFLOW.')
          CALL LCMGET(KPLIB,'BIN-NTOT0',SIGT(LLL+1,ISO))
          CALL LCMGET(KPLIB,'BIN-SIGS00',SIGS(LLL+1,ISO))
          CALL LCMLEN(KPLIB,'BIN-NUSIGF',ILEN,ITYLCM)
          IF(ILEN.GT.0) CALL LCMGET(KPLIB,'BIN-NUSIGF',SIGF(LLL+1,ISO))
        ENDIF
*----
*  ELASTIC SCATTERING INFORMATION.
*----
        MAXIII=MAXTRA-III(ISO)+1
        CALL LIBPRI(MAXIII,DELI,AWR,IALTER,0,NEXT(ISO),PRI(III(ISO)))
      ELSE
        ! infinite dilution isotope
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
          WRITE(HSMG,'(18H AUTFLU: ISOTOPE '',A12,7H'' (ISO=,I8,5H) IS ,
     1    39HNOT AVAILABLE IN THE ORIGINAL MICROLIB.)') HNAMIS,ISO
          CALL XABORT(HSMG)
        ENDIF
        IF((.NOT.LLIB).AND.(IMPX.GT.2)) THEN
          CALL LCMGTC(KPLI0,'ALIAS',12,1,HNABIS)
          WRITE(6,'(26H AUTFLU: RECOVER ISOTOPE '',A12,11H'' FROM THE ,
     1    12HNEW LIBRARY.)') HNABIS
        ELSE IF(LLIB.AND.(IMPX.GT.2)) THEN
          WRITE(6,'(26H AUTFLU: RECOVER ISOTOPE '',A12,11H'' FROM THE ,
     1    12HOLD LIBRARY.)') HNAMIS
        ENDIF
        NSCAR=0
        SCAT(:NGRP,:NGRP)=0.0
        CALL LCMLEN(IPP,'NTOT0',N13,ITYLCM)
        IF(N13.EQ.0) THEN
           CALL LCMLIB(IPP)
           CALL XABORT('AUTFLU: NO INFINITE-DILUTION TOTAL XS.')
        ELSE IF(N13.GT.LBIN) THEN
           CALL XABORT('AUTFLU: LBIN OVERFLOW.')
        ELSE IF(N13.NE.NGRP) THEN
           CALL XABORT('AUTFLU: INVALID X-SECTIONS.')
        ENDIF
        CALL LCMGET(IPP,'NTOT0',SIGT(1,ISO))
        CALL LCMLEN(IPP,'NUSIGF',N10,ITYLCM)
        IF(N10.GT.0) CALL LCMGET(IPP,'NUSIGF',SIGF(1,ISO))
        CALL LCMLEN(IPP,'SIGS00',N12,ITYLCM)
        IF(N12.GT.0) THEN
          CALL LCMGET(IPP,'SIGS00',SIGS(1,ISO))
          CALL LCMLEN(IPP,'SIGS01',N14,ITYLCM)
          IF(N14.GT.0) THEN
            CALL LCMGET(IPP,'SIGS01',SIGS1(1,ISO))
          ELSE
            SIGS1(:NGRP,ISO)=0.0
          ENDIF
        ELSE
          CALL LCMGET(IPP,'SCAT00',SGAR)
          CALL LCMGET(IPP,'NJJS00',NJJ)
          CALL LCMGET(IPP,'IJJS00',IJJ)
          IGAR1=0
          DO IG2=1,NGRP
            DO IG1=IJJ(IG2),IJJ(IG2)-NJJ(IG2)+1,-1
              IGAR1=IGAR1+1
              SCAT(IG1,IG2)=SGAR(IGAR1)
            ENDDO
          ENDDO
          DO IG1=1,NGRP
            SUMSC=0.0D0
            DO IG2=1,NGRP
              SUMSC=SUMSC+SCAT(IG1,IG2)
            ENDDO
            SIGS(IG1,ISO)=REAL(SUMSC)
          ENDDO
          CALL LCMLEN(IPP,'SCAT01',N87,ITYLCM)
          IF(N87.GT.0) THEN
             CALL LCMGET(IPP,'SCAT01',SGAR)
             CALL LCMGET(IPP,'NJJS01',NJJ)
             CALL LCMGET(IPP,'IJJS01',IJJ)
             IGAR1=0
             DO IG2=1,NGRP
               DO IG1=IJJ(IG2),IJJ(IG2)-NJJ(IG2)+1,-1
                 IGAR1=IGAR1+1
                 SCAT(IG1,IG2)=SGAR(IGAR1)
               ENDDO
             ENDDO
             DO IG1=1,NGRP
               SUMSC=0.0D0
               DO IG2=1,NGRP
                 SUMSC=SUMSC+SCAT(IG1,IG2)
               ENDDO
               SIGS1(IG1,ISO)=REAL(SUMSC)
             ENDDO
          ENDIF
        ENDIF
        ! compute SIGGAR used by SPH equivalence
        IF((DENN.NE.0.0).AND.(IBM.NE.0).AND.(JRES.EQ.0)) THEN
          DO 70 IG1=1,NGRP
          SIGGAR(IBM,JRES,IG1,1)=SIGGAR(IBM,JRES,IG1,1)+DENN*
     1    SIGT(IG1,ISO)
          SIGGAR(IBM,JRES,IG1,3)=SIGGAR(IBM,JRES,IG1,3)+DENN*
     1    SIGS(IG1,ISO)
   70     CONTINUE
          CALL LCMLEN(IPP,'TRANC',LENGT,ITYLCM)
          IF(LENGT.GT.0) THEN
            CALL LCMGET(IPP,'TRANC',GA1)
            DO 80 IG1=1,NGRP
            SIGGAR(IBM,JRES,IG1,2)=SIGGAR(IBM,JRES,IG1,2)+DENN*GA1(IG1)
   80       CONTINUE
          ENDIF
        ENDIF
        ! expand xs of non-resonant isotopes
        CALL AUTPRD(NGRP,LBIN,NBIN,SIGS(1,ISO))
        CALL AUTPRD(NGRP,LBIN,NBIN,SIGT(1,ISO))
        CALL AUTPRD(NGRP,LBIN,NBIN,SIGF(1,ISO))
        CALL AUTPRD(NGRP,LBIN,NBIN,SIGS1(1,ISO))
      ENDIF
      III(ISO+1)=III(ISO)+NEXT(ISO)
*----
*  CROSS SECTION EDITION.
*----
      IF(IMPX.GT.7) THEN
         CALL LCMGTC(KPLIB,'ALIAS',12,1,HNAMIS)
         WRITE(6,540) HNAMIS,(K,SIGS(K,ISO),SIGT(K,ISO),SIGF(K,ISO),
     1   SIGS1(K,ISO),K=1,LBIN)
         I2=III(ISO+1)-1
         WRITE(6,550) III(ISO),I2,(PRI(K),K=III(ISO),I2)
      ENDIF
  120 CONTINUE
      CALL KDRCPU(TK2)
      IF(IMPX.GT.1) WRITE(6,'(/36H AUTFLU: CPU TIME SPENT TO RECOVER A,
     1 33HUTOLIB AND INFINITE-DILUTION XS =,F8.1,8H SECOND./)') TK2-TK1
*
      CALL KDRCPU(TK1)
      TK4=0.0
      TK5=0.0
      ICPIJ=0
*----
*  SET THE NUMBER DENSITIES.
*----
      CONC(:NBMIX,:NBISO)=0.0
      DO 130 ISO=1,NBISO
      IBM=MIX(ISO)
      IF(IBM.LE.0) CYCLE
      CONC(IBM,ISO)=DEN(ISO)
  130 CONTINUE
*----
*  DETERMINE WHICH GROUPS ARE SELF-SHIELDED.
*----
       MASKG(:NGRP)=.FALSE.
       MASKG(IGRMIN:IGRMAX)=.TRUE.
*----
*  COMPUTE THE NEUTRON FLUX.
*----
      CALL KDRCPU(TKA)
      CALL KDRCPU(TKA)
      DO 330 IG1=1,NGRP
      ICPIJ=ICPIJ+NBIN(IG1)
  330 CONTINUE
      CALL LCMOP(IPSYS,'**SYSTEM**',0,1,0)
      KNORM=1
      IMPY=MAX(0,IMPX-3)
      IF(IPHASE.EQ.1) THEN
        CALL AUTIT2(IPTRK,IFTRAK,IPSYS,MAXTRA,KNORM,NUN,LBIN,NREG,
     1  NBMIX,NBISO,MAT,VOL,KEYFLX,NIRES,IAPT,CDOOR,LEAKSW,TITR,IMPY,
     2  CONC,SIGS,SIGT,SIGS1,DIL,PRI,UUU,DELI,ITRANC,NEXT,III,PHI)
      ELSE IF(IPHASE.EQ.2) THEN
        CALL AUTIT1(IPTRK,IFTRAK,IPSYS,MAXTRA,KNORM,LBIN,NREG,
     1  NBMIX,NBISO,MAT,VOL,NIRES,IAPT,CDOOR,LEAKSW,TITR,IMPY,CONC,
     2  SIGS,SIGT,SIGS1,DIL,PRI,UUU,DELI,ITRANC,NEXT,III,PHI)
      ENDIF
      CALL LCMCL(IPSYS,2)
      CALL KDRCPU(TKB)
      TK4=TK4+(TKB-TKA)
*     ***************************************************************
      CALL LCMVAL(IPLI0,' ')
*----
*  RESET MASKG FOR SPH CALCULATION IN SMALL LETHARGY WIDTH GROUPS.
*----
      DO 380 IG1=1,NGRP
      IF(MASKG(IG1)) THEN
        IF(DELTAU(IG1).GT.0.1) MASKG(IG1)=.TRUE.
      ENDIF
  380 CONTINUE
      CALL KDRCPU(TK2)
      IF(IMPX.GT.1) WRITE(6,'(/34H AUTFLU: CPU TIME SPENT TO COMPUTE,
     1 31H SELF-SHIELDED REACTION RATES =,F8.1,19H SECOND, INCLUDING:
     2 /9X,F8.1,36H SECOND TO SOLVE FOR AUTOSECOL FLUX;/9X,7HNUMBER ,
     3 25HOF ASSEMBLY DOORS CALLS =,I5,1H.)') TK2-TK1,TK4,ICPIJ
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(III,NEXT,PRI)
      DEALLOCATE(SCAT,SGAR,IJJ,NJJ)
      DEALLOCATE(IPISO2,IPISO1)
      DEALLOCATE(SIGINF,MASKH,DELTAU,CONC,GA2,GA1)
      RETURN
  540 FORMAT(/18H AUTFLU: ISOTOPE ',A12,1H'/12X,4HSIGS,16X,4HSIGT,
     1 16X,4HSIGF,16X,5HSIGS1/(I5,1P,4E20.7))
  550 FORMAT(/29H AUTFLU: SCATTERING ELEMENTS:,2I10/(1P,10E13.5))
      END
