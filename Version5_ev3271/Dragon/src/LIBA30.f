*DECK LIBA30
      SUBROUTINE LIBA30 (IPLIB,NAMFIL,NGRO,NBISO,NL,ISONAM,ISONRF,
     1 IPISO,MASKI,TN,LSHI,SN,SB,IMPX,NGF,NGFR,NDEL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transcription of the useful interpolated microscopic cross section
* data from APOLIB-3 to LCM data structures.
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
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* NAMFIL  name of the APOLIB-3 file in HDF5 format.
* NGRO    number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* NL      number of Legendre orders required in the calculation
*         NL=1 or higher.
* ISONAM  alias name of isotopes.
* ISONRF  library reference name of isotopes.
* IPISO   pointer array towards microlib isotopes.
* MASKI   isotopic mask. Isotope with index I is processed if
*         MASKI(I)=.true.
* TN      temperature of each isotope.
* LSHI    resonant region number associated with each isotope.
*         Infinite dilution will be assumed if LSHI(i)=0.
* SN      dilution cross section in each energy group of each
*         isotope. a value of 1.0E10 is used for infinite dilution.
* SB      dilution cross section as used by Livolant and Jeanpierre
*         normalization.
* IMPX    print flag.
*
*Parameters: output
* NGF     number of fast groups without self-shielding.
* NGFR    number of fast and resonance groups.
* NDEL    number of precursor groups for delayed neutrons.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPISO(NBISO)
      INTEGER NGRO,NBISO,NL,ISONAM(3,NBISO),ISONRF(3,NBISO),LSHI(NBISO),
     1 IMPX,NGF,NGFR,NDEL
      REAL TN(NBISO),SN(NGRO,NBISO),SB(NGRO,NBISO)
      CHARACTER NAMFIL*(*)
      LOGICAL MASKI(NBISO)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPAP1,IPAP2
      PARAMETER (IOUT=6)
      TYPE(C_PTR) KPLIB
      CHARACTER RECNAM*80,RECNA2*80,TEXT80*80,HNAMIS*12,HNISOR*12,
     1 HSMG*131,TEXT12*12,CFILNA1*64,CFILNA2*64
      LOGICAL L104,LSIGS,LABSO,LFISS,LDIF
      INTEGER RANK,TYPE,NBYTE,DIMSR(5)
      DOUBLE PRECISION XDRCST,DSUM
      REAL TKT(5)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITYPRO,ORANIS,ENRANG,
     1 FSTTMP,TMPMON,ADDTMP,ITEMPA,ISPAOF,IAFAG,IFAGR,FLXADD
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IPR
      REAL, ALLOCATABLE, DIMENSION(:) :: ENERG,DELTA,SECT,XSTOT,TAUX,
     1 AMASS,TEMP,TEMPM,XS,WGTFLX,BGXS,ABSOXS,DIFFXS,FISSXS,DK104
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIGS
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SCAT
      CHARACTER(LEN=24), ALLOCATABLE, DIMENSION(:) :: NOM,NOMS,HREANM
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IPR(2,NBISO),ITYPRO(NL))
      ALLOCATE(SECT(NGRO),SIGS(NGRO,NL),SCAT(NGRO,NGRO,NL),XSTOT(NGRO))
*
      ANEUT=REAL(XDRCST('Neutron mass','amu'))
      NGF=NGRO+1
      NGFR=0
      NDEL=0
      IF(IMPX.GT.0) WRITE (IOUT,800) NAMFIL
*----
*  OPEN THE APOLIB-3 FILE.
*----
      IND = INDEX(NAMFIL, ":")
      IF(IND.EQ.0) THEN
        CFILNA1=NAMFIL
        CFILNA2=" "
      ELSE
        CFILNA1=NAMFIL(:IND-1)
        CFILNA2=NAMFIL(IND+1:)
      ENDIF
      CALL hdf5_open_file(CFILNA1, IPAP1, .TRUE.)
      IF(IMPX.GT.0) THEN
        CALL hdf5_read_data(IPAP1,"Head/LibraryInfo",TEXT80)
        WRITE (IOUT,810) TEXT80
        WRITE (IOUT,'(40H LIBA30: NUMBER OF ISOTOPES IN MICROLIB=,I6)')
     1  NBISO
      ENDIF
      CALL hdf5_read_data(IPAP1,"Head/nbIs",NISOT)
      CALL hdf5_read_data(IPAP1,"Head/IsNames",NOM)
      IF(IMPX.GE.10) THEN
        DO ISO=1,NISOT
          WRITE(IOUT,'(8H -----> ,A)') TRIM(NOM(ISO))
        ENDDO
      ENDIF
      NISOTS=0
      IF(CFILNA2.NE.' ') THEN
        CALL hdf5_open_file(CFILNA2, IPAP2, .TRUE.)
        CALL hdf5_read_data(IPAP2,"Isotopes/NIsotope",NISOTS)
        CALL hdf5_read_data(IPAP2,"Isotopes/IsoNames",NOMS)
        IF(IMPX.GE.10) THEN
          DO ISO=1,NISOTS
            WRITE(IOUT,'(8H SS---> ,A)') TRIM(NOMS(ISO))
          ENDDO
        ENDIF
      ENDIF
*----
*  RECOVER INFORMATION FROM EnergyMesh GROUP
*----
      CALL hdf5_read_data(IPAP1, "EnergyMesh/nbGr", NGRI)
      CALL hdf5_read_data(IPAP1, "EnergyMesh/EnMshInMeV", ENERG)
      CALL hdf5_read_data(IPAP1, "EnergyMesh/EnGrInLtg", DELTA)
      ENERG(:NGRO+1)=ENERG(:NGRO+1)*1.E6
      IF(NGRI.NE.NGRO) CALL XABORT('LIBA30: INVALIB NGRO.')
      CALL LCMPUT(IPLIB,'ENERGY',NGRO+1,2,ENERG)
      CALL LCMPUT(IPLIB,'DELTAU',NGRO,2,DELTA)
*----
*  RECOVER INFORMATION FROM PhysicalData GROUP
*----
      CALL hdf5_read_data(IPAP1, "PhysicalData/AtomicMass", AMASS)
      IF(SIZE(AMASS).NE.NISOT) CALL XABORT('LIBA30: INVALIB NISOT.')
      DO IA=1,NISOT
        AMASS(IA)=AMASS(IA)/ANEUT
      ENDDO
*----
*  RECOVER INFORMATION FROM WeightFlux GROUP
*----
      CALL hdf5_read_data(IPAP1, "WeightFlux/nbFluxTypes", NBFLX)
      IF(NBFLX.GT.0) THEN
        CALL hdf5_read_data(IPAP1, "WeightFlux/FlxAdd", FLXADD)
        CALL hdf5_read_data(IPAP1, "WeightFlux/WgtFlx", WGTFLX)
      ENDIF
*----
*  SET THE CORRESPONDANCE BETWEEN THE APOLIB AND THE LIST OF ISOTOPES.
*----
      IF(IMPX.GT.1) WRITE(IOUT,820) NISOT,NISOTS
      CALL XDISET(IPR,2*NBISO,0)
      CALL KDRCPU(TK1)
      DO 50 IMX=1,NBISO
      IF(MASKI(IMX)) THEN
         WRITE(HNISOR,'(3A4)') (ISONRF(I0,IMX),I0=1,3)
         KISO=0
         DO 10 ISO=1,NISOT
         IF(HNISOR.EQ.NOM(ISO)) THEN
           KISO=ISO
           GO TO 20
         ENDIF
   10    CONTINUE
         WRITE (HSMG,780) HNISOR,TRIM(CFILNA1)
         CALL XABORT(HSMG)
   20    IPR(1,IMX)=KISO
*
         IF((NISOTS.GT.0).AND.(LSHI(IMX).GT.0)) THEN
           KISO=0
           DO 30 ISO=1,NISOTS
           IF(HNISOR.EQ.NOMS(ISO)) THEN
             KISO=ISO
             GO TO 40
           ENDIF
   30      CONTINUE
           WRITE (HSMG,790) HNISOR,TRIM(CFILNA2)
           CALL XABORT(HSMG)
   40      IPR(2,IMX)=KISO
         ENDIF
      ENDIF
   50 CONTINUE
      DEALLOCATE(NOM)
      IF(NISOTS.GT.0) DEALLOCATE(NOMS)
      CALL KDRCPU(TK2)
      TKT(1)=TK2-TK1
*----
*  RECOVER INFORMATION FROM TemperatureM GROUP
*----
      CALL hdf5_read_data(IPAP1, "TemperatureM/TempMshAdd", ITEMPA)
      CALL hdf5_read_data(IPAP1, "TemperatureM/TempMesh", TEMPM)
*----
*  PROCESS INFINITE DILUTION INFORMATION.
*----
      CALL KDRCPU(TK1)
      DO 560 IMX=1,NBISO
      KISEG=IPR(1,IMX)
      IF(KISEG.GT.0) THEN
        WRITE(HNAMIS,'(3A4)') (ISONAM(I0,IMX),I0=1,3)
        WRITE(HNISOR,'(3A4)') (ISONRF(I0,IMX),I0=1,3)
        IF(IMPX.GT.0) WRITE (IOUT,830) HNAMIS,HNISOR
        IF(IMPX.GT.1) WRITE(IOUT,'(/29H LIBA30: PROCESSING ISOTOPE '',
     1  A,2H''.)') TRIM(HNISOR)
        WRITE(TEXT80,'(18HAPOLIB-3 ISOTOPE: ,A)') TRIM(HNISOR)
*----
*  RECOVER INFORMATION FROM Dimensions GROUP
*----
        WRITE(RECNAM,'(10HIsotopeXS/,A,12H/Dimensions/)') TRIM(HNISOR)
        CALL hdf5_read_data(IPAP1, TRIM(RECNAM)//"nbRea", NBREA)
        CALL hdf5_read_data(IPAP1, TRIM(RECNAM)//"orAnis", ORANIS)
        CALL hdf5_read_data(IPAP1, TRIM(RECNAM)//"nbTemp", NBTMP)
        CALL hdf5_read_data(IPAP1, TRIM(RECNAM)//"nbTypEner", NBTENR)
*----
*  RECOVER INFORMATION FROM Info GROUP
*----
        WRITE(RECNAM,'(10HIsotopeXS/,A,6H/Info/)') TRIM(HNISOR)
        CALL hdf5_read_data(IPAP1, TRIM(RECNAM)//"WgtFlxON", IWFLON)
        CALL hdf5_read_data(IPAP1, TRIM(RECNAM)//"isFissile", ISFIS)
        CALL hdf5_read_data(IPAP1, TRIM(RECNAM)//"isTranProb", ITPROB)
        CALL hdf5_read_data(IPAP1, TRIM(RECNAM)//"FstTmpDepGr", FSTTMP)
        CALL hdf5_read_data(IPAP1, TRIM(RECNAM)//"EnergyRange", ENRANG)
        CALL hdf5_read_data(IPAP1, TRIM(RECNAM)//"ReaNames", HREANM)
        CALL hdf5_read_data(IPAP1, TRIM(RECNAM)//"ChiErgMshInd", ICHIEG)
        CALL hdf5_read_data(IPAP1, TRIM(RECNAM)//"TempMshON", TMPMON)
        CALL hdf5_read_data(IPAP1, TRIM(RECNAM)//"addrTempIntp", ADDTMP)
        CALL hdf5_read_data(IPAP1, TRIM(RECNAM)//"isPartOf", ISPAOF)
        WRITE(RECNAM,'(10HIsotopeXS/,A,12H/ReactionXS/)') TRIM(HNISOR)
        WRITE(RECNA2,'(10HIsotopeXS/,A,19H/Profile/SCATTProf/)')
     1  TRIM(HNISOR)
        DO JMX=IMX,NBISO
          IF(IPR(1,JMX).EQ.KISEG) THEN
            KPLIB=IPISO(JMX) ! set JMX-th isotope
            CALL LCMLEN(KPLIB,'ALIAS',ILENG,ITYLCM)
            IF(ILENG.EQ.0) THEN
              WRITE(HNAMIS,'(3A4)') (ISONAM(I0,JMX),I0=1,3)
              CALL LCMPTC(KPLIB,'ALIAS',12,1,HNAMIS)
              IF(IPR(1,JMX).LE.0) CALL XABORT('LIBA30: BAD AWR.')
              CALL LCMPUT(KPLIB,'AWR',1,2,AMASS(IPR(1,JMX)))
              CALL LCMPTC(KPLIB,'README',80,1,TEXT80)
            ENDIF
            IF(NBFLX.GT.0) THEN
              IOF=FLXADD(IWFLON+1)+1
              SECT(:NGRO)=WGTFLX(IOF:IOF+NGRO-1)
              CALL LCMPUT(KPLIB,'NWT0',NGRO,2,SECT)
            ENDIF
            LSIGS=.FALSE.
            LABSO=.FALSE.
            DO I=1,NBREA
              IGR0=ENRANG(2*I-1)+1
              NBGR=ENRANG(2*I)
              IFGTD=FSTTMP(I)
              IF(IFGTD.GE.1) THEN
                NTDG=NBGR-IFGTD+1 ! number of temp-dependent groups
                NBTMP2=NBTMP
                MSHIND=TMPMON(I)+1
                IADD=ITEMPA(MSHIND)
                IF(ITEMPA(MSHIND+1)-ITEMPA(MSHIND).NE.NBTMP) THEN
                  CALL XABORT('LIBA30: INVALID NBTMP.')
                ENDIF
              ELSE
                NTDG=0
                NBTMP2=1
                IADD=0
              ENDIF
              NGDG=NBGR-IGR0+1 ! number of groups in energy range
              IF(IMPX.GT.2) THEN
                WRITE(IOUT,860) TRIM(HREANM(I)),NGDG,NTDG
                IF(ISPAOF(I).GE.0) WRITE(IOUT,870) HREANM(ISPAOF(I)+1)
                IF(IFGTD.GE.1) WRITE(IOUT,880) TEMPM(IADD+1:IADD+NBTMP)
              ENDIF
              IND=LEN(TRIM(HREANM(I)))
              CALL hdf5_read_data(IPAP1, TRIM(RECNAM)//HREANM(I), XS)
              NSECT0=SIZE(XS)
              IF(HREANM(I)(IND-3:IND).EQ.'TOTA') THEN
                IF(NSECT0.NE.NGDG+(NBTMP2-1)*NTDG) THEN
                  WRITE(HSMG,'(33HLIBA30: INVALID SIZE FOR ISOTOPE ,A,
     1            14H AND REACTION ,A,7H. SIZE=,I6,11H SHOULD BE=,I6,
     2            7H. NGDF=,I6,6H NTDG=,I6,7H NBTMP=,I6,1H.)') 
     3            TRIM(HNISOR),TRIM(HREANM(I)),NSECT0,
     4            NGDG+(NBTMP2-1)*NTDG,NGDG,NTDG,NBTMP2
                  WRITE(IOUT,'(/1X,A)') HSMG
                  GO TO 550
                ENDIF
                SECT(:NGRO)=0.0
                IF(IFGTD.GE.1) THEN
                  CALL LIBA22(NGDG,TN(JMX),NBTMP,NSECT0,IFGTD,
     1            TEMPM(IADD+1),XS,SECT(IGR0))
                ELSE
                  IF(NSECT0.NE.NGDG) CALL XABORT('LIBA30: INVALID NSEC'
     1            //'T0(1).')
                  SECT(IGR0:IGR0+NGDG-1)=XS(:NSECT0)
                ENDIF
                IF(HREANM(I).EQ.'ABSO-TOTA') LABSO=.TRUE.
                TEXT12=HREANM(I)(:12)
                IF(TEXT12.EQ.'MT16-TOTA') TEXT12='N2N'
                IF(TEXT12.EQ.'MT17-TOTA') TEXT12='N3N'
                IF(TEXT12.EQ.'MT28-TOTA') TEXT12='NNP'
                IF(TEXT12.EQ.'MT37-TOTA') TEXT12='N4N'
                IF(TEXT12.EQ.'MT102-TOTA') TEXT12='NG'
                IF(TEXT12.EQ.'MT103-TOTA') TEXT12='NP'
                IF(TEXT12.EQ.'MT104-TOTA') TEXT12='ND'
                IF(TEXT12.EQ.'MT105-TOTA') TEXT12='NT'
                IF(TEXT12.EQ.'MT107-TOTA') TEXT12='NA'
                IF(TEXT12.EQ.'MT108-TOTA') TEXT12='N2A'
                IF(TEXT12.EQ.'FISS-TOTA') TEXT12='NFTOT'
                IF(TEXT12.EQ.'NUFISS-TOTA') TEXT12='NUSIGF'
                IF(TEXT12.EQ.'CHI-TOTA') TEXT12='CHI'
                CALL LCMPUT(KPLIB,TEXT12,NGRO,2,SECT)
              ELSE IF(HREANM(I)(IND-3:IND).EQ.'JUMP') THEN
                IF(ORANIS(I).LE.0) THEN
                  CALL XABORT('LIBA30: INVALID JUMP ANISOTROPY.')
                ELSE IF(NSECT0.NE.(NGDG+(NBTMP2-1)*NTDG)*ORANIS(I)) THEN
                  CALL XABORT('LIBA30: INVALID JUMP SIZE.')
                ENDIF
                IF(HREANM(I)(:4).EQ.'SCAT') THEN
                  SIGS(:NGRO,:NL)=0.0
                  DO IL=1,MIN(ORANIS(I),NL)
                    IOF1=(IL-1)*(NGDG+(NBTMP2-1)*NTDG)+1
                    IOF2=IL*(NGDG+(NBTMP2-1)*NTDG)
                    IF(IFGTD.GE.1) THEN
                      CALL LIBA22(NGDG,TN(JMX),NBTMP,NSECT0,IFGTD,
     1                TEMPM(IADD+1),XS(IOF1),SIGS(IGR0,IL))
                    ELSE
                      IF(NSECT0.NE.NGDG*ORANIS(I)) CALL XABORT('LIBA30'
     1                //': INVALID NSECT0(2).')
                      SIGS(IGR0:IGR0+NGDG-1,IL)=XS(IOF1:IOF2)
                    ENDIF
                  ENDDO
                  LSIGS=.TRUE.
                  IF(.NOT.LABSO) CALL XABORT('LIBA30: NO ABSO-TOTA.')
                  SECT(:NGRO)=0.0
                  CALL LCMGET(KPLIB,'ABSO-TOTA',SECT)
                  DO IG=1,NGRO
                    SECT(IG)=SECT(IG)+SIGS(IG,1)
                  ENDDO
                  CALL LCMPUT(KPLIB,'NTOT0',NGRO,2,SECT)
                  CALL LCMLEN(KPLIB,'NXN-TOTA',ILENG,ITYLCM)
                  IF(ILENG.GT.0) THEN
                    CALL LCMGET(KPLIB,'NXN-TOTA',SECT)
                    DO IG=1,NGRO
                      SIGS(IG,1)=SIGS(IG,1)+SECT(IG)
                    ENDDO
                  ENDIF
                ENDIF
              ELSE IF(HREANM(I)(IND-3:IND).EQ.'PROF') THEN
                IF(HREANM(I)(:4).EQ.'SCAT') THEN
                 IF(.NOT.LSIGS) CALL XABORT('LIBA30: SIGS NOT SET.')
                 CALL hdf5_read_data(IPAP1,TRIM(RECNA2)//"AddressFAG",
     1           IAFAG)
                 IF(SIZE(IAFAG).NE.(NBGR+1)*ORANIS(I)) CALL XABORT('LI'
     1           //'BA30: INVALID AddressFAG SIZE.')
                 CALL hdf5_read_data(IPAP1,TRIM(RECNA2)//"FstArrGroup",
     1           IFAGR)
                 NV=0
                 DO IL=1,ORANIS(I)
                   DO IG=1,NGRO ! departure group
                     NV=NV+(IAFAG((IL-1)*(NGRO+1)+IG+1)-IAFAG((IL-1)*
     1               (NGRO+1)+IG))
                   ENDDO
                   IF(IFGTD.GE.1) THEN
                     DO IG=IFGTD,NGRO ! departure group
                       NV=NV+(NBTMP2-1)*(IAFAG((IL-1)*(NGRO+1)+IG+1)-
     1                 IAFAG((IL-1)*(NGRO+1)+IG))
                     ENDDO
                   ENDIF
                 ENDDO
                 IF(NSECT0.NE.NV) CALL XABORT('LIBA30: INVALID NSECTO('
     1           //'3).')
                 ILMIN=MIN(ORANIS(I),NL)
                 CALL LIBA33(NBGR,ILMIN,TN(JMX),NBTMP,NSECT0,IFGTD,
     1           TEMPM(IADD+1),IAFAG,IFAGR,XS,SCAT)
                 IF(ITPROB.NE.0) THEN
                   DO IL=1,ILMIN
                     DO IG=1,NBGR
                       SCAT(:NBGR,IG,IL)=SCAT(:NBGR,IG,IL)*SIGS(IG,IL)
                     ENDDO
                   ENDDO
                 ELSE
                   DO IL=1,ILMIN
                     DO IG=1,NBGR
                       DSUM=SUM(SCAT(:NGRO,IG,IL))
                       SCAT(:NBGR,IG,IL)=SCAT(:NBGR,IG,IL)*SIGS(IG,IL)/
     1                 REAL(DSUM)
                     ENDDO
                   ENDDO
                 ENDIF
                 DEALLOCATE(IFAGR,IAFAG)
                 CALL XDRLGS(KPLIB,1,IMPX,0,ILMIN-1,1,NBGR,SIGS,SCAT,
     1           ITYPRO)
                ENDIF
              ELSE
                CALL XABORT('LIBA30: TOTA/JUMP/PROF SUFFIX EXPECTED.')
              ENDIF
  550         DEALLOCATE(XS)
            ENDDO
            IF(IMPX.GT.1) CALL LCMLIB(KPLIB)
          ENDIF
        ENDDO
        DO JMX=IMX,NBISO
          IF(IPR(1,JMX).EQ.KISEG) IPR(1,JMX)=0
        ENDDO
        DEALLOCATE(ISPAOF,ADDTMP,TMPMON,HREANM,ENRANG,FSTTMP,ORANIS)
      ENDIF
  560 CONTINUE
      DEALLOCATE(TEMPM,ITEMPA,AMASS)
      CALL KDRCPU(TK2)
      TKT(2)=TK2-TK1
*----
*  PROCESS SELF-SHIELDING DATA.
*----
      L104=.FALSE.
      LABSO=.TRUE.
      LDIF=.TRUE.
      CALL KDRCPU(TK1)
      DO 600 IMX=1,NBISO
      KISEG=IPR(2,IMX)
      IF(KISEG.GT.0) THEN
        WRITE(HNISOR,'(3A4)') (ISONRF(I0,IMX),I0=1,3)
        IF(IMPX.GT.1) WRITE(IOUT,'(/31H LIBA30: PROCESSING SELF-SHIELD,
     1  12HED ISOTOPE '',A,2H''.)') TRIM(HNISOR)
*----
*  RECOVER INFORMATION FROM Dimensions GROUP
*----
        WRITE(RECNAM,'(9HIsotopes/,A,11H/HomoRates/)') TRIM(HNISOR)
        IF(.NOT.hdf5_group_exists(IPAP2,TRIM(RECNAM))) THEN
          WRITE(HSMG,'(38HLIBA30: missing HomoRates in group ,A,1H.)')
     1    TRIM(RECNAM)
          CALL XABORT(HSMG)
        ENDIF
        CALL hdf5_read_data(IPAP2, TRIM(RECNAM)//"FirstGrp", IGR0)
        CALL hdf5_read_data(IPAP2, TRIM(RECNAM)//"LastGrp", JGR0)
        CALL hdf5_read_data(IPAP2, TRIM(RECNAM)//"NbOfGrp", NBGR)
        CALL hdf5_read_data(IPAP2, TRIM(RECNAM)//"Temp", TEMP)
        CALL hdf5_read_data(IPAP2, TRIM(RECNAM)//"BgXS", BGXS)
        CALL hdf5_read_data(IPAP2, TRIM(RECNAM)//"AbsoRate", ABSOXS)
        CALL hdf5_read_data(IPAP2, TRIM(RECNAM)//"DiffRate", DIFFXS)
        CALL hdf5_info(IPAP2,TRIM(RECNAM)//"FissRate",RANK,TYPE,NBYTE,
     1   DIMSR)
        NGF=MIN(NGF,IGR0)
        NGFR=MAX(NGFR,JGR0)
        LFISS=(TYPE.NE.99)
        NBTMP=SIZE(TEMP)
        NBDIL=SIZE(BGXS)
        IF(IMPX.GT.1) THEN
          WRITE(IOUT,910) (BGXS(I),I=1,NBDIL)
          WRITE(IOUT,920) (TEMP(I),I=1,NBTMP)
          WRITE(IOUT,930) IGR0,JGR0,NBGR,NBDIL,NBTMP
        ENDIF
        IF(LFISS) THEN
          CALL hdf5_read_data(IPAP2, TRIM(RECNAM)//"FissRate",FISSXS)
        ELSE
          ALLOCATE(FISSXS(NBDIL*NBGR*NBTMP))
          FISSXS(:NBDIL*NBGR*NBTMP)=0.0
        ENDIF
        ALLOCATE(TAUX(7*NBGR),DK104(NBDIL*NBGR*NBTMP))
        DK104(:NBDIL*NBGR*NBTMP)=0.0
        DO JMX=IMX,NBISO
          IF(IPR(2,JMX).EQ.KISEG) THEN
            WRITE(HNAMIS,'(3A4)') (ISONAM(I0,JMX),I0=1,3)
            KPLIB=IPISO(JMX) ! set JMX-th isotope
            IF(IMPX.GT.3) WRITE(IOUT,'(/17H LIBA30: PROCESS ,A12,1H:)')
     1      HNAMIS
            CALL LIBA34(HNAMIS,NGRO,IGR0,NBGR,NBDIL,NBTMP,LFISS,L104,
     1      BGXS,TEMP,TN(JMX),SN(1,JMX),ABSOXS,DIFFXS,FISSXS,DK104,
     2      IMPX,TAUX)
*
*           COMPUTE THE SELF-SHIELDED FLUX AND CROSS SECTIONS.
            CALL LIBA25(KPLIB,LABSO,LDIF,LFISS,L104,NGRO,IGR0,NBGR,
     1      NBDIL,NL,BGXS,SN(1,JMX),SB(1,JMX),DELTA,ISONAM(1,JMX),
     2      TAUX,IMPX)
          ENDIF
        ENDDO
        DO JMX=IMX,NBISO
          IF(IPR(2,JMX).EQ.KISEG) IPR(2,JMX)=0
        ENDDO
        DEALLOCATE(TAUX,DK104,FISSXS,DIFFXS,ABSOXS,BGXS,TEMP)
      ENDIF
  600 CONTINUE
      CALL KDRCPU(TK2)
      TKT(3)=TK2-TK1
*----
*  CHECK IF ALL REACTIONS HAVE BEEN PROCESSED.
*----
      DO 575 IMX=1,NBISO
      DO 570 I=1,2
      IF(IPR(I,IMX).NE.0) THEN
         WRITE(HSMG,950) I,(ISONAM(I0,IMX),I0=1,3)
         CALL XABORT(HSMG)
      ENDIF
  570 CONTINUE
  575 CONTINUE
      IF(IMPX.GT.2) WRITE(IOUT,940) (TKT(I),I=1,3)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      IF(NBFLX.GT.0) DEALLOCATE(WGTFLX,FLXADD)
      DEALLOCATE(DELTA,ENERG)
      DEALLOCATE(XSTOT,SCAT,SIGS,SECT)
      DEALLOCATE(ITYPRO,IPR)
      RETURN
*
  780 FORMAT(26HLIBA30: MATERIAL/ISOTOPE ',A,22H' IS MISSING ON APOLIB,
     1 13H-3 FILE NAME ,A,1H.)
  790 FORMAT(49HLIBA30: SELF-SHIELDING DATA OF MATERIAL/ISOTOPE ',A12,
     1 35H' IS MISSING ON APOLIB-3 FILE NAME ,A,1H.)
  800 FORMAT(/43H LIBA30: PROCESSING APOLIB-3 LIBRARY NAME: ,A,1H.)
  810 FORMAT(/32H LIBA30: X-SECTION LIBRARY INFO:/9X,A80/)
  820 FORMAT(/35H LIBA30: PROBING THE APOLIB-3 FILE./9X,11HNUMBER OF I,
     1 29HSOTOPES AT INFINITE DILUTION=,I8/9X,21HNUMBER OF SELF-SHIELD,
     2 12HED ISOTOPES=,I8)
  830 FORMAT(/30H PROCESSING ISOTOPE/MATERIAL ',A12,11H' (HNISOR=',A12,
     1 3H').)
  860 FORMAT(/9X,5H---- ,A,5H ----/9X,29HNUMBER OF GROUPS IN ENERGY RA,
     1 4HNGE=,I5/10X,32HNUMBER OF TEMP-DEPENDENT GROUPS=,I5)
  870 FORMAT(9X,21HGLOBAL REACTION NAME=,A)
  880 FORMAT(9X,13HTEMPERATURES=,1P,9E12.4/(22X,9E12.4))
  910 FORMAT(/9X,10HDILUTIONS=,1P,9E12.4/(19X,9E12.4))
  920 FORMAT(/9X,28HSELF-SHIELDING TEMPERATURES=,1P,7E12.4/(37X,7E12.4))
  930 FORMAT(/9X,5HIGR0=,I4,6H JGR0=,I4,6H NBGR=,I4,7H NBDIL=,I4,
     1 7H NBTMP=,I4)
  940 FORMAT(/26H LIBA30: CPU TIME USAGE --,F10.2,9H INDEXING/26X,
     1 F10.2,24H INFINITE DILUTION P0 XS/26X,F10.2,16H DILUTION-DEPEND,
     2 11HENT XS DATA)
  950 FORMAT(26HLIBA30: REMAINING REACTION,I3,14H FOR ISOTOPE ',3A4,
     1 2H'.)
      END
