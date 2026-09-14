*DECK LIBDRA
      SUBROUTINE LIBDRA (IPLIB,IPDRL,NAMFIL,NGRO,NBISO,NL,ISONAM,
     1 ISONRF,IPISO,TN,SN,SB,MASKI,NED,HVECT,IMPX,NGF,NGFR,NDEL,NBESP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transcription of the useful interpolated microscopic cross section
* data from a microscopic x-section library (draglib format) to LCM
* data structures.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* IPDRL   pointer to the draglib (L_DRAGLIB signature).
* NAMFIL  name of the Dragon library file.
* NGRO    number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* NL      number of Legendre orders required in the calculation
*         NL=1 or higher.
* ISONAM  alias name of isotopes.
* ISONRF  library name of isotopes.
* IPISO   pointer array towards microlib isotopes.
* TN      temperature of each isotope.
* SN      dilution cross section in each energy group of each
*         isotope. a value of 1.0e10 is used for infinite dilution.
* SB      dilution cross section as used in livolant and jeanpierre
*         normalization.
* MASKI   isotopic mask. Isotope with index I is processed if
*         MASKI(I)=.true.
* NED     number of extra vector edits.
* HVECT   names of the extra vector edits.
* IMPX    print flag.
*
*Parameters: output
* NGF     number of fast groups without self-shielding.
* NGFR    number of fast and resonance groups.
* NDEL    number of precursor groups for delayed neutrons.
* NBESP   number of energy-dependent fission spectra.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      PARAMETER(MAXDEL=10,MAXESP=4)
      CHARACTER*(*) HVECT(NED),NAMFIL
      TYPE(C_PTR) IPLIB,IPDRL,IPISO(NBISO)
      INTEGER NGRO,NBISO,NL,ISONAM(3,NBISO),ISONRF(3,NBISO),NED,IMPX,
     1 NGF,NGFR,NDEL,NBESP
      REAL TN(NBISO),SN(NGRO,NBISO),SB(NGRO,NBISO)
      LOGICAL MASKI(NBISO)
*----
*  LOCAL VARIABLES
*----
      CHARACTER CD*4,HSMG*131,HVERS*12,HNISOR*12,HNAMIS*12,HNUSIG*12,
     1 HCHI*12
      PARAMETER (IOUT=6,MAXTMP=50,NOTX=3)
      TYPE(C_PTR) KPLIB
      LOGICAL LSIGF,LGOLD,LOGT,LNZERO
      INTEGER IESP(MAXESP+1)
      DOUBLE PRECISION FACTOR,TERP(MAXTMP),DDELI
      REAL TEMP(MAXTMP),ZLAMB(MAXDEL),EESP(MAXESP+1)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NFS,ITYPRO
      REAL, ALLOCATABLE, DIMENSION(:) :: AWR,DELTA,TOTAL,GOLD,ENER,
     1 BIN,EBIN,SIGS2,SCAT2,TOTAL2,SIGF2,CHI2,SADD2,GOLD2,BIN2,CHI4G2
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIGS,SIGF,CHI,SADD,CHI4G
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SCAT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ZNPHI,ZNPHI2
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LSCAT,LADD
      CHARACTER(LEN=72), ALLOCATABLE, DIMENSION(:) :: HTITLE
      REAL NFTOT(NGRO),HFACT(NGRO),KAPF
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NFS(NGRO),ITYPRO(NL))
      ALLOCATE(AWR(NBISO),DELTA(NGRO),SIGS(NGRO,NL),SCAT(NGRO,NGRO,NL),
     1 TOTAL(NGRO),SIGF(NGRO,0:MAXDEL),CHI(NGRO,0:MAXDEL),
     2 SADD(NGRO,NED),ENER(NGRO+1),GOLD(NGRO),ZNPHI(NGRO))
      ALLOCATE(LSCAT(NL),LADD(NED))
*----
*  RECOVER THE GROUP STRUCTURE.
*----
      NGF=NGRO+1
      NGFR=0
      NDEL=0
      HVERS='**UNKNOWN**'
      CALL LCMLEN(IPDRL,'VERSION',LENGT,ITYLCM)
      IF(LENGT.NE.0) CALL LCMGTC(IPDRL,'VERSION',12,HVERS)
      IF(IMPX.GT.0) WRITE (IOUT,900) TRIM(NAMFIL),TRIM(HVERS)
      IF(HVERS.EQ.'RELEASE_2003') THEN
        WRITE(IOUT,'(46H LIBDRA: ***WARNING*** RELEASE_2003 DRAGLIBS A,
     1  15HRE DEPRECIATED.)')
      ENDIF
      CALL LCMLEN(IPDRL,'README',LENGT,ITYLCM)
      IF((IMPX.GT.0).AND.(LENGT.GT.0)) THEN
         LENGT=(LENGT-1)/18+1
         ALLOCATE(HTITLE(LENGT))
         CALL LCMGTC(IPDRL,'README',72,LENGT,HTITLE)
         WRITE (IOUT,940)
         DO 10 J=1,LENGT
         WRITE (IOUT,'(1X,A)') TRIM(HTITLE(J))
   10    CONTINUE
         DEALLOCATE(HTITLE)
         WRITE (IOUT,'(40H LIBDRA: NUMBER OF ISOTOPES IN MICROLIB=,I6)')
     1   NBISO
      ENDIF
      CALL LCMLEN(IPDRL,'ENERGY',LENGT,ITYLCM)
      LENGT=LENGT-1
      IF(LENGT.NE.NGRO) CALL XABORT('LIBDRA: INVALID GROUP STRUCTURE.')
      CALL LCMGET(IPDRL,'ENERGY',ENER)
      CALL LCMLEN(IPDRL,'DELTAU',LENGT,ITYLCM)
      IF(LENGT.EQ.NGRO) THEN
         CALL LCMGET(IPDRL,'DELTAU',DELTA)
      ELSE IF(LENGT.EQ.0) THEN
         IF(ENER(NGRO+1).EQ.0.0) ENER(NGRO+1)=1.0E-5
         DO 15 J=1,NGRO
         DELTA(J)=LOG(ENER(J)/ENER(J+1))
   15    CONTINUE
      ENDIF
      CALL LCMPUT(IPLIB,'ENERGY',NGRO+1,2,ENER)
      CALL LCMPUT(IPLIB,'DELTAU',NGRO,2,DELTA)
      CALL LCMLEN(IPDRL,'CHI-LIMITS',NBESP,ITYLCM)
      IF(NBESP.GT.0) THEN
         NBESP=NBESP-1
         IF(NBESP.GT.MAXESP) CALL XABORT('LIBDRA: MAXESP OVERFLOW.')
         CALL LCMGET(IPDRL,'CHI-LIMITS',IESP)
         CALL LCMPUT(IPLIB,'CHI-LIMITS',NBESP+1,1,IESP)
         CALL LCMGET(IPDRL,'CHI-ENERGY',EESP)
         CALL LCMPUT(IPLIB,'CHI-ENERGY',NBESP+1,2,EESP)
      ENDIF
      ALLOCATE(CHI4G(NGRO,NBESP))
*----
*  READ THROUGH DRAGON FILE AND ACCUMULATE CROSS SECTIONS FOR THIS RANGE
*  OF MATS, LEGENDRE ORDERS, AND GROUPS.
*----
      DO 400 IMX=1,NBISO
      IF(MASKI(IMX)) THEN
         WRITE(HNAMIS,'(3A4)') (ISONAM(I0,IMX),I0=1,3)
         WRITE(HNISOR,'(3A4)') (ISONRF(I0,IMX),I0=1,3)
         CALL LCMLEN(IPDRL,HNISOR,LENGT,ITYLCM)
         IF(LENGT.EQ.0) THEN
            CALL LCMLIB(IPDRL)
            WRITE (HSMG,910) HNAMIS,HNISOR,NAMFIL,IMX
            CALL XABORT(HSMG)
         ENDIF
         IF(IMPX.GT.0) WRITE (IOUT,920) HNAMIS,HNISOR
         CALL LCMSIX(IPDRL,HNISOR,1)
*
         CALL LCMGET(IPDRL,'AWR',AWR(IMX))
         CALL LCMLEN(IPDRL,'README',LTITLE,ITYLCM)
         IF(LTITLE.GT.0) THEN
           LTITLE=(LTITLE-1)/18+1
           ALLOCATE(HTITLE(LTITLE))
           CALL LCMGTC(IPDRL,'README',72,LTITLE,HTITLE)
           IF(IMPX.GT.0) THEN
             WRITE (IOUT,930)
             DO 20 J=1,LTITLE
             WRITE (IOUT,'(1X,A)') TRIM(HTITLE(J))
   20        CONTINUE
           ENDIF
         ENDIF
*----
*  RECOVER BIN TYPE INFORMATION (IF AVAILABLE).
*----
         LBIN=0
         CALL LCMLEN (IPDRL,'BIN-NFS',LENGT,ITYXSM)
         IF(LENGT.GT.0) THEN
            CALL LCMGET (IPDRL,'BIN-NFS',NFS)
            DO 30 I=1,NGRO
            LBIN=LBIN+NFS(I)
   30       CONTINUE
            ALLOCATE(BIN(3*LBIN),EBIN(LBIN+1))
            CALL LCMGET (IPDRL,'BIN-ENERGY',EBIN)
            CALL LCMLEN(IPDRL,'BIN-DELI',LENDEL,ITYLCM)
            IF((LENDEL.EQ.1).AND.(ITYLCM.EQ.2)) THEN
              CALL LCMGET (IPDRL,'BIN-DELI',RDELI)
            ELSE IF((LENDEL.EQ.1).AND.(ITYLCM.EQ.4)) THEN
              CALL LCMGET (IPDRL,'BIN-DELI',DDELI)
              RDELI=REAL(DDELI)
            ENDIF
         ENDIF
*----
*  RECOVER ECCOLIB INFORMATION (IF AVAILABLE).
*----
         DELECC=0.0
         IGECCO=0
         CALL LCMLEN(IPDRL,'ONFLIGHTDEL',LENGT,ITYLCM)
         IF(LENGT.GT.0) THEN
            CALL LCMGET(IPDRL,'ONFLIGHTDEL',DELECC)
            CALL LCMGET(IPDRL,'ONFLIGHTIGR',IGECCO)
         ENDIF
*
         CALL LCMLEN (IPDRL,'TEMPERATURE',NTMP,ITYLCM)
         IF(NTMP.GT.MAXTMP) CALL XABORT('LIBDRA: MAXTMP OVERFLOW.')
         IF(NTMP.EQ.0) THEN
            CALL LCMLEN (IPDRL,'LAMBDA-D',NDEL0,ITYLCM)
            NDEL=MAX(NDEL,NDEL0)
            IF(NDEL0.GT.MAXDEL) CALL XABORT('LIBDRA: MAXDEL OVERFLOW.')
            IF(NDEL0.GT.0) CALL LCMGET (IPDRL,'LAMBDA-D',ZLAMB)
            CALL LIBDRB (IPDRL,HNAMIS,NGRO,NL,NDEL0,NBESP,ENER,
     1      SN(1,IMX),SB(1,IMX),NED,HVECT,DELTA,LBIN,NFS,EBIN,AWR(IMX),
     2      DELECC,IGECCO,IMPX,NGF,NGFR,LSCAT,LSIGF,LADD,LGOLD,
     3      SIGS(1,1),SCAT(1,1,1),TOTAL,ZNPHI,SIGF(1,1),CHI(1,1),
     4      CHI4G(1,1),SADD(1,1),GOLD(1),BIN(1))
         ELSE
*----
*  PERFORM TEMPERATURE LAGRANGIAN INTERPOLATION (ORDER ABS(NOTX)).
*----
            CALL LCMSIX (IPDRL,'SUBTMP0001',1)
            CALL LCMLEN (IPDRL,'LAMBDA-D',NDEL0,ITYLCM)
            NDEL=MAX(NDEL,NDEL0)
            IF(NDEL0.GT.MAXDEL) CALL XABORT('LIBDRA: MAXDEL OVERFLOW.')
            IF(NDEL0.GT.0) CALL LCMGET (IPDRL,'LAMBDA-D',ZLAMB)
            CALL LCMSIX (IPDRL,' ',2)
            CALL LCMGET (IPDRL,'TEMPERATURE',TEMP)
            CALL LIBLEX(NTMP,TN(IMX),TEMP,NOTX,TERP)
            TOTAL(:NGRO)=0.0
            ZNPHI(:NGRO)=0.0D0
            SIGF(:NGRO,:NDEL0)=0.0
            CHI(:NGRO,:NDEL0)=0.0
            CHI4G(:NGRO,:NBESP)=0.0
            GOLD(:NGRO)=0.0
            SIGS(:NGRO,:NL)=0.0
            SCAT(:NGRO,:NGRO,:NL)=0.0
            SADD(:NGRO,:NED)=0.0
            IF(LBIN.GT.0) BIN(:3*LBIN)=0.0
            ALLOCATE(SIGS2(NGRO*NL),SCAT2(NGRO*NGRO*NL),TOTAL2(NGRO),
     1      SIGF2(NGRO*(NDEL0+1)),CHI2(NGRO*(NDEL0+1)),SADD2(NGRO*NED),
     2      GOLD2(NGRO),BIN2(3*LBIN),ZNPHI2(NGRO),CHI4G2(NGRO*NBESP))
            FACTOR=1.0D0
            DO 210 ITM=1,NTMP
            TERPM=REAL(TERP(ITM))
            FACTOR=FACTOR-TERP(ITM)
            IF(TERPM.EQ.0.0) GO TO 210
            IF(IMPX.GT.4) WRITE(6,'(/30H DRAGLIB ACCESS AT TEMPERATURE,
     >      1P,E12.4,18H KELVIN. FACTOR = ,E12.4)') TEMP(ITM),TERPM
            WRITE (CD,'(I4.4)') ITM
            CALL LCMSIX (IPDRL,'SUBTMP'//CD,1)
            CALL LIBDRB (IPDRL,HNAMIS,NGRO,NL,NDEL0,NBESP,ENER,
     1      SN(1,IMX),SB(1,IMX),NED,HVECT,DELTA,LBIN,NFS,EBIN,AWR(IMX),
     2      DELECC,IGECCO,IMPX,NGF,NGFR,LSCAT,LSIGF,LADD,LGOLD,SIGS2(1),
     3      SCAT2(1),TOTAL2,ZNPHI2,SIGF2(1),CHI2(1),CHI4G2(1),SADD2(1),
     4      GOLD2(1),BIN2(1))
            CALL LCMSIX (IPDRL,' ',2)
            DO 130 IG=1,NGRO
              TOTAL(IG)=TOTAL(IG)+TERPM*TOTAL2(IG)
              ZNPHI(IG)=ZNPHI(IG)+TERPM*ZNPHI2(IG)
  130       CONTINUE
            IF(LSIGF) THEN
              DO 141 IDEL=0,NDEL0
              DO 140 IG=1,NGRO
                IOFSET=IDEL*NGRO+IG-1
                SIGF(IG,IDEL)=SIGF(IG,IDEL)+TERPM*SIGF2(IOFSET+1)
                CHI(IG,IDEL)=CHI(IG,IDEL)+TERPM*CHI2(IOFSET+1)
  140         CONTINUE
  141         CONTINUE
              DO 146 ISP=1,NBESP
              DO 145 IG=1,NGRO
                IOFSET=(ISP-1)*NGRO+IG-1
                CHI4G(IG,ISP)=CHI4G(IG,ISP)+TERPM*CHI4G2(IOFSET+1)
  145         CONTINUE
  146         CONTINUE
            ENDIF
            DO 160 IL=1,NL
              IF(LSCAT(IL)) THEN
                DO 150 IG2=1,NGRO
                  SIGS(IG2,IL)=SIGS(IG2,IL)+TERPM*SIGS2((IL-1)*NGRO+IG2)
                  IOF=(IL-1)*NGRO*NGRO+(IG2-1)*NGRO
                  DO 151 IG1=1,NGRO
                    SCAT(IG1,IG2,IL)=SCAT(IG1,IG2,IL)+TERPM*
     >              SCAT2(IOF+IG1)
  151             CONTINUE
  150           CONTINUE
              ENDIF
  160       CONTINUE
            DO 180 IED=1,NED
            IF(LADD(IED)) THEN
               DO 170 IG=1,NGRO
                 SADD(IG,IED)=SADD(IG,IED)+TERPM*SADD2((IED-1)*NGRO+IG)
  170          CONTINUE
            ENDIF
  180       CONTINUE
            IF(LGOLD) THEN
               DO 190 IG=1,NGRO
               GOLD(IG)=GOLD(IG)+TERPM*GOLD2(IG)
  190          CONTINUE
            ENDIF
            DO 200 IG=1,3*LBIN
            BIN(IG)=BIN(IG)+TERPM*BIN2(IG)
  200       CONTINUE
  210       CONTINUE
            DEALLOCATE(CHI4G2,ZNPHI2,BIN2,GOLD2,SADD2,CHI2,SIGF2,TOTAL2,
     >      SCAT2,SIGS2)
            IF(ABS(FACTOR).GT.1.0D-4) CALL XABORT('LIBDRA: TERP ERROR')
         ENDIF
         CALL LCMSIX(IPDRL,' ',2)

*----
*        SAVE CROSS SECTION DATA ON LCM.
*----
         KPLIB=IPISO(IMX) ! set IMX-th isotope
         CALL LCMPTC(KPLIB,'ALIAS',12,HNAMIS)
         CALL LCMPUT(KPLIB,'AWR',1,2,AWR(IMX))
         IF(LTITLE.GT.0) THEN
            CALL LCMPTC(KPLIB,'README',72,LTITLE,HTITLE)
            DEALLOCATE(HTITLE)
         ENDIF
         DO 220 IG=1,NGRO
         IF(TOTAL(IG).LT.0.0) THEN
            WRITE(HSMG,'(42HLIBDRA: NEGATIVE TOTAL CROSS SECTION IN GR,
     1      3HOUP,I4,14H FOR ISOTOPE '',A12,2H''.)') IG,HNAMIS
            CALL XABORT(HSMG)
         ELSE IF(ZNPHI(IG).LT.0.0) THEN
            WRITE(HSMG,'(41HLIBDRA: NEGATIVE INTEGRATED FLUX IN GROUP,
     1      I4,14H FOR ISOTOPE '',A12,2H''.)') IG,HNAMIS
            CALL XABORT(HSMG)
         ENDIF
  220    CONTINUE
         CALL LCMPUT(KPLIB,'NTOT0',NGRO,2,TOTAL)
         TOTAL(:NGRO)=REAL(ZNPHI(:NGRO))
         CALL LCMPUT(KPLIB,'NWT0',NGRO,2,TOTAL)
         IF(NDEL0.GT.0) CALL LCMPUT (KPLIB,'LAMBDA-D',NDEL0,2,ZLAMB)
         IF(LSIGF) THEN
            DO 250 IDEL=0,NDEL0
            IF(IDEL.EQ.0) THEN
               HNUSIG='NUSIGF'
            ELSE
               WRITE(HNUSIG,'(6HNUSIGF,I2.2)') IDEL
            ENDIF
            CALL LCMPUT(KPLIB,HNUSIG,NGRO,2,SIGF(1,IDEL))
            IF(IDEL.EQ.0) THEN
               IF(NBESP.GT.0) GO TO 250
               HCHI='CHI'
            ELSE
               WRITE(HCHI,'(3HCHI,I2.2)') IDEL
            ENDIF
            CALL LCMPUT(KPLIB,HCHI,NGRO,2,CHI(1,IDEL))
  250       CONTINUE
            DO 260 ISP=1,NBESP
               LNZERO=.FALSE.
               DO 255 IG=1,NGRO
               LNZERO=LNZERO.OR.(CHI4G(IG,ISP).NE.0.0)
  255          CONTINUE
               IF(LNZERO) THEN
                 WRITE(HCHI,'(5HCHI--,I2.2)') ISP
                 CALL LCMPUT(KPLIB,HCHI,NGRO,2,CHI4G(1,ISP))
               ENDIF
  260       CONTINUE
         ENDIF
         CALL XDRLGS(KPLIB,1,0,0,NL-1,1,NGRO,SIGS,SCAT,ITYPRO)
         CALL LCMLEN(KPLIB,'H-FACTOR',LENGT,ITYLCM)
         IF(LENGT.NE.0) CALL LCMDEL(KPLIB,'H-FACTOR')
         DO 340 IED=1,NED
         IF(LADD(IED).AND.(HVECT(IED)(:3).NE.'CHI')
     1               .AND.(HVECT(IED)(:2).NE.'NU')
     2               .AND.(HVECT(IED).NE.'NTOT0')
     3               .AND.(HVECT(IED)(:3).NE.'NWT')) THEN
            CALL LCMPUT(KPLIB,HVECT(IED),NGRO,2,SADD(1,IED))
         ENDIF
  340    CONTINUE
         IF(LGOLD) CALL LCMPUT(KPLIB,'NGOLD',NGRO,2,GOLD)
         IF(LBIN.GT.0) THEN
            CALL LCMPUT(KPLIB,'BIN-NFS',NGRO,1,NFS)
            CALL LCMPUT(KPLIB,'BIN-ENERGY',LBIN+1,2,EBIN)
            CALL LCMPUT(KPLIB,'BIN-NTOT0',LBIN,2,BIN)
            CALL LCMPUT(KPLIB,'BIN-SIGS00',LBIN,2,BIN(LBIN+1))
            LOGT=.FALSE.
            DO 350 I=1,LBIN
            LOGT=LOGT.OR.(BIN(2*LBIN+I).NE.0.0)
  350       CONTINUE
            IF(LOGT) THEN
               CALL LCMPUT(KPLIB,'BIN-NUSIGF',LBIN,2,BIN(2*LBIN+1))
            ENDIF
            DEALLOCATE(EBIN,BIN)
            IF(LENDEL.EQ.1) CALL LCMPUT(KPLIB,'BIN-DELI',1,2,RDELI)
         ENDIF
         IF(IMPX.GT.9) CALL LCMLIB(KPLIB)
      ENDIF
  400 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(LADD,LSCAT)
      DEALLOCATE(CHI4G,ZNPHI,GOLD,ENER,SADD,CHI,SIGF,TOTAL,SCAT,SIGS,
     1 DELTA,AWR)
      DEALLOCATE(ITYPRO,NFS)
      RETURN
*
  900 FORMAT(/33H PROCESSING DRAGON LIBRARY NAMED ,A,9H VERSION ,A,1H.)
  910 FORMAT(26HLIBDRA: MATERIAL/ISOTOPE ',A12,5H' = ',A12,9H' IS MISS,
     1 25HING ON DRAGON FILE NAMED ,A12,10H (ISOTOPE=,I10,2H).)
  920 FORMAT(/30H PROCESSING ISOTOPE/MATERIAL ',A12,11H' (HNISOR=',A12,
     1 3H').)
  930 FORMAT(/23H ISOTOPE/MATERIAL INFO:)
  940 FORMAT(/24H X-SECTION LIBRARY INFO:)
      END
