*DECK MCRLIB
      SUBROUTINE MCRLIB(MAXNIS,MAXISO,IPLIB,IPMPO,IACCS,NMIX,NGRP,IMPX,
     1 HEQUI,HMASL,NCAL,HEDIT,ITER,MY1,MY2,NBISO,TERP,NISO,LISO,HISO,
     2 CONC,ITODO,MIXC,LRES,LPURE,ILUPS,B2,VTOT,YLDS,DECAYC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Build the Microlib by scanning the NCAL elementary calculations in
* a MPO file and weighting them with TERP factors.
*
*Copyright:
* Copyright (C) 2022 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* MAXNIS  maximum value of NISO(I) in user data.
* MAXISO  maximum allocated space for output Microlib TOC information.
* IPLIB   address of the output Microlib LCM object.
* IPMPO   pointer to the MPO file.
* IACCS   =0 Microlib is created; =1 ... is updated.
* NMIX    maximum number of material mixtures in the Microlib.
* NGRP    number of energy groups.
* IMPX    print parameter (equal to zero for no print).
* HEQUI   keyword of SPH-factor set to be recovered.
* HMASL   keyword of MASL data set to be recovered.
* NCAL    number of elementary calculations in the MPO file.
* HEDIT   name of output group for a (multigroup mesh, output geometry)
*         couple (generally equal to 'output_0').
* ITER    completion flag (=0: compute the macrolib).
* MY1     number of fissile isotopes including macroscopic sets.
* MY2     number of fission fragment.
* NBISO   number of particularized isotopes in the MPO file.
* TERP    interpolation factors.
* NISO    number of user-selected isotopes.
* LISO    type of treatment (=.true.: ALL; =.false.: ONLY).
* HISO    name of the user-selected isotopes.
* CONC    user-defined number density of the user-selected isotopes. A
*         value of -99.99 is set to indicate that the MPO file value is
*         used.
* ITODO   non-depletion mask (=1 to force a user-selected isotope to be
*         non-depleting)
* MIXC    mixture index in the MPO file corresponding to each Microlib
*         mixture. Equal to zero if a Microlib mixture is not updated.
* LRES    =.true. if the interpolation is done without updating isotopic 
*         densities
* LPURE   =.true. if the interpolation is a pure linear interpolation 
*         with TERP factors.
* ILUPS   up-scattering removing flag (=1 to remove up-scattering from
*         output cross-sections).
* B2      buckling
* VTOT    volume of updated core.
* YLDS    fission yields.
* DECAYC  radioactive decay constants.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPMPO
      INTEGER MAXNIS,MAXISO,IACCS,NMIX,NGRP,IMPX,NCAL,ITER,MY1,MY2,
     > NBISO,NISO(NMIX),ITODO(NMIX,MAXNIS),MIXC(NMIX),ILUPS
      REAL TERP(NCAL,NMIX),CONC(NMIX,MAXNIS),B2
      DOUBLE PRECISION VTOT,YLDS(MY1,MY2),DECAYC(NBISO)
      LOGICAL LISO(NMIX),LRES,LPURE
      CHARACTER(LEN=80) HEQUI,HMASL
      CHARACTER(LEN=12) HEDIT
      CHARACTER(LEN=8) HISO(NMIX,NBISO)
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      INTEGER, PARAMETER::MAXREA=50
      INTEGER, PARAMETER::NSTATE=40
      INTEGER, PARAMETER::MAXFRD=4
      TYPE(C_PTR) JPLIB,KPLIB
      REAL FACT0, WEIGHT, DEN
      INTEGER I, J, I0, IBM, IBMOLD, ICAL, IED2, IFISS, IGR, ILONG, IDF,
     > IPRC, IREA, IREAB, IREAF, ISO, ITRANC, ITSTMP, ITYLCM, IY1, IY2,
     > JSO, KSO, KSO1, LMY1, LSO, MAXMIX, NBISO2, NBISO2I, NBS1, NCALS,
     > NED2, NL, NMIL, NPAR, NPRC, NREA, NSURFD, NISOF, NISOP, NISOS,
     > RANK, NBYTE, TYPE, DIMSR(5), ILOC, NADDRXS, NLOC, NMGF, ID, ID_E,
     > ID_G, NENERG, NGEOME, ADDRZI, ISOM, NISOM, IGRC, NALBP, NALBP2
      CHARACTER RECNAM*80,RECNA2*80,TEXT8*8,TEXT12*12,HSMG*131,
     > HVECT2(MAXREA)*8,HRESID*8
      INTEGER ISTATE(NSTATE),INAME(2),IHRES(2)
      REAL TMPDAY(3)
      LOGICAL LUSER,LSTRD,LSPH,LMASL,LALBG
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IMIX2,ISONA,ISOMI,ITOD2,
     > ISTY1,ISTY2,IPIFI,IMICR,ITOD1,JJSO,IPYMIX,LOCAD,REACTION,ISOTOPE,
     > ADDRISO,IGYELD,IADRY,DIMS_MPO
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: HUSE2,HNAM2,OUPUTID
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: ADDRXS
      REAL, ALLOCATABLE, DIMENSION(:) :: DENS2,DENS3,VOL2,VOLMI2,SPH,
     > ENER,XVOLM,CONCE,TAUXFI,NWT0,FLUXS,DENIS,GAR1,GAR2,LAMB,BETAR,
     > INVELS,BETARB,INVELSB,DECAY2,RVALO
      REAL, ALLOCATABLE, DIMENSION(:,:) :: DENS1,FACT,CHIRS,CHIRSB,
     > TAUXGF
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: XS,SIGS,DENS0,FLUX,YLDS2
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: SS2D
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: YLDSM
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LXS,MASK,MASKL
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HPYNAM
      CHARACTER(LEN=24), ALLOCATABLE, DIMENSION(:) :: TEXT24,NOMREA,
     > NOMISO
      CHARACTER(LEN=80), ALLOCATABLE, DIMENSION(:) :: LOCTYP,LOCKEY
*----
*  RECOVER MPO FILE CHARACTERISTICS
*----
      I=0
      CALL MPOTOC(IPMPO,HEDIT,0,NREA,I0,NMIL,NPAR,NLOC,NISOF,NISOP,
     > NISOS,NCALS,I,NSURFD,NALBP,NPRC)
      IF(NBISO.NE.I0) CALL XABORT('MCRLIB: INVALID VALUE OF NBISO.')
      IF(NGRP.NE.I) CALL XABORT('MCRLIB: INVALID VALUE OF NGRP.')
      IF(NREA+2.GT.MAXREA) CALL XABORT('MCRLIB: MAXREA OVERFLOW')
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IMIX2(MAXISO),ITOD2(MAXISO),ISTY1(MAXISO),ISTY2(MAXISO),
     > HUSE2(3,MAXISO),HNAM2(3,MAXISO))
      ALLOCATE(DENS2(MAXISO),DENS3(MAXISO),VOL2(MAXISO),VOLMI2(NMIX),
     > FLUX(NMIX,NGRP,2),SPH(NGRP))
*----
*  MICROLIB INITIALIZATION
*----
      VOLMI2(:NMIX)=0.0
      DENS2(:MAXISO)=0.0
      VOL2(:MAXISO)=0.0
      IMIX2(:MAXISO)=0
      ITOD2(:MAXISO)=0
      ISTY2(:MAXISO)=0
      IF(IACCS.EQ.0) THEN
        IF(LRES) CALL XABORT('MCRLIB: RES OPTION IS INVALID.')
        NBISO2=0
        NED2=0
        TEXT12='L_LIBRARY'
        CALL LCMPTC(IPLIB,'SIGNATURE',12,1,TEXT12)
      ELSE
        CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
        IF(ISTATE(1).NE.NMIX) CALL XABORT('MCRLIB: INVALID NUMBER OF '
     1  //'MATERIAL MIXTURES IN THE MICROLIB.')
        IF(ISTATE(3).NE.NGRP) CALL XABORT('MCRLIB: INVALID NUMBER OF '
     1  //'ENERGY GROUPS IN THE MICROLIB.')
        NBISO2=ISTATE(2)
        IF(NBISO2.GT.MAXISO) CALL XABORT('MCRLIB: MAXISO OVERFLOW(1).')
        NED2=ISTATE(13)
        IF(NED2.GT.MAXREA) CALL XABORT('MCRLIB: MAXREA OVERFLOW.')
        CALL LCMLEN(IPLIB,'MIXTURESVOL',ILONG,ITYLCM)
        IF(ILONG.GT.0) THEN
          CALL LCMGET(IPLIB,'MIXTURESVOL',VOLMI2)
        ELSE
          VOLMI2(:NMIX)=0.0
        ENDIF
        CALL LCMGET(IPLIB,'ISOTOPESUSED',HUSE2)
        CALL LCMGET(IPLIB,'ISOTOPERNAME',HNAM2)
        CALL LCMGET(IPLIB,'ISOTOPESDENS',DENS2)
        CALL LCMGET(IPLIB,'ISOTOPESVOL',VOL2)
        CALL LCMGET(IPLIB,'ISOTOPESMIX',IMIX2)
        CALL LCMGET(IPLIB,'ISOTOPESTODO',ITOD2)
        CALL LCMGET(IPLIB,'ISOTOPESTYPE',ISTY2)
        IF(NED2.GT.0) CALL LCMGTC(IPLIB,'ADDXSNAME-P0',8,NED2,HVECT2)
      ENDIF
*----
*  SET EQUIVALENCE AND HEAVY DENSITY FLAGS.
*----
      LSPH=.FALSE.
      LMASL=.FALSE.
      NLOC=0
      IF(hdf5_group_exists(IPMPO,"/local_values/")) THEN
        CALL hdf5_read_data(IPMPO,"/local_values/LOCVALTYPE",LOCTYP)
        CALL hdf5_read_data(IPMPO,"/local_values/LOCVALNAME",LOCKEY)
        NLOC=SIZE(LOCTYP,1)
        DO ILOC=1,NLOC
          LSPH=LSPH.OR.((LOCTYP(ILOC).EQ.'EQUI').AND.
     >                  (LOCKEY(ILOC).EQ.HEQUI))
          LMASL=LMASL.OR.((LOCTYP(ILOC).EQ.'HEAVY_METAL_DENSITY').AND.
     >                    (LOCKEY(ILOC).EQ.HMASL))
        ENDDO
      ENDIF
*----
*  FIND SCATTERING ANISOTROPY.
*----
      WRITE(RECNAM,'(8H/output/,A,6H/info/)') TRIM(HEDIT)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"NADDRXS",NADDRXS)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ADDRXS",ADDRXS)
      NL=0
      DO I=1,NADDRXS
        DO ISO=1,NBISO
          NL=MAX(NL,ADDRXS(NREA+1,ISO,I))
          NL=MAX(NL,ADDRXS(NREA+2,ISO,I))
        ENDDO
      ENDDO
      IF(IMPX.GT.1) THEN
        WRITE(6,'(37H MCRLIB: number of legendre orders  =,I4)') NL
      ENDIF
*----
*  SET ENERGY MESH AND ZONE VOLUMES
*----
      CALL hdf5_read_data(IPMPO,"/energymesh/NENERGYMESH",NENERG)
      CALL hdf5_read_data(IPMPO,"/geometry/NGEOMETRY",NGEOME)
      CALL hdf5_read_data(IPMPO,"/output/OUPUTID",OUPUTID)
      READ(HEDIT,'(7X,I2)') ID
      ID_G=0
      ID_E=0
      DO I=1,NGEOME
        DO J=1,NENERG
          IF(OUPUTID(J,I).EQ.ID) THEN
            ID_G=I-1
            ID_E=J-1
            GO TO 10
          ENDIF
        ENDDO
      ENDDO
      CALL XABORT('MCRLIB: no ID found in /output/OUPUTID.')
   10 WRITE(RECNAM,'(23H/energymesh/energymesh_,I0)') ID_E
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"/ENERGY",ENER)
      IF(SIZE(ENER,1)-1.NE.NGRP) CALL XABORT('MCRLIB: INVALID NGRP VAL'
     > //'UE.')
      DO IGR=1,NGRP+1
        ENER(IGR)=ENER(IGR)/1.0E-6
      ENDDO
      WRITE(RECNAM,'(19H/geometry/geometry_,I0,1H/)') ID_G
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ZONEVOLUME",XVOLM)
      CALL LCMPUT(IPLIB,'ENERGY',NGRP+1,2,ENER)
      DO IGR=1,NGRP
        ENER(IGR)=LOG(ENER(IGR)/ENER(IGR+1))
      ENDDO
      CALL LCMPUT(IPLIB,'DELTAU',NGRP,2,ENER)
      DEALLOCATE(ENER)
*----
*  RECOVER INFORMATION ON REACTIONS AND ISOTOPE NAMES
*----
      IREAB=0
      IREAF=0
      WRITE(RECNAM,'(8H/output/,A,6H/info/)') TRIM(HEDIT)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"REACTION",REACTION)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ISOTOPE",ISOTOPE)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ADDRISO",ADDRISO)
      NBISO=ADDRISO(SIZE(ADDRISO,1))
      IF(NBISO.EQ.0) CALL XABORT('MCRLIB: NO CROSS SECTIONS.')
      ALLOCATE(NOMREA(NREA+2),NOMISO(NBISO))
      IF(NREA.GT.0) THEN
        CALL hdf5_read_data(IPMPO,"/contents/reactions/REACTIONAME",
     >  TEXT24)
        DO I=1,NREA
          NOMREA(I)=TEXT24(REACTION(I)+1)
        ENDDO
        DEALLOCATE(TEXT24,REACTION)
        DO IREA=1,NREA
          IF(NOMREA(IREA).EQ.'Absorption') THEN
            IREAB=IREA
            EXIT
          ENDIF
        ENDDO
        DO IREA=1,NREA
          IF(NOMREA(IREA).EQ.'NuFission') THEN
            IREAF=IREA
            EXIT
          ENDIF
        ENDDO
      ENDIF
      NOMREA(NREA+1)='Total'
      NOMREA(NREA+2)='Leakage'
      NREA=NREA+2
      CALL hdf5_read_data(IPMPO,"/contents/isotopes/ISOTOPENAME",TEXT24)
      DO I=1,NBISO
        NOMISO(I)=TEXT24(ISOTOPE(I)+1)
      ENDDO
      DEALLOCATE(TEXT24,ADDRISO,ISOTOPE)
      IF(IMPX.GT.1) THEN
        WRITE(6,'(/24H MCRLIB: reaction names:)')
        DO I=1,NREA
          WRITE(6,'(5X,7HNOMREA(,I3,2H)=,A)') I,TRIM(NOMREA(I))
        ENDDO
        WRITE(6,'(/23H MCRLIB: isotope names:)')
        DO I=1,NBISO
          WRITE(6,'(5X,7HNOMISO(,I3,2H)=,A)') I,TRIM(NOMISO(I))
        ENDDO
      ENDIF
*----
*  LOOP OVER MPO MIXTURES TO COMPUTE DENS0(NMIL,NCAL,NBISO)
*----
      ALLOCATE(DENS0(NMIL,NCAL,NBISO))
      DENS0(:NMIL,:NCAL,:NBISO)=0.0
      DO 30 IBMOLD=1,NMIL
      DO ICAL=1,NCAL
        DO IBM=1,NMIX
          IF((TERP(ICAL,IBM).NE.0.0).AND.(MIXC(IBM).EQ.IBMOLD)) GO TO 15
        ENDDO
        CYCLE
   15   WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,6H/zone_,I0,1H/)')
     >  TRIM(HEDIT),ICAL-1,IBMOLD-1
        IF(NBISO.GT.0) THEN
          CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"CONCENTRATION",CONCE)
          DO 20 ISO=1,NBISO
            DENS0(IBMOLD,ICAL,ISO)=CONCE(ISO)
   20     CONTINUE
          DEALLOCATE(CONCE)
        ENDIF
      ENDDO
   30 CONTINUE
*----
*  LOOP OVER MICROLIB MIXTURES
*----
      YLDS(:MY1,:MY2)=0.0D0
      DECAYC(:NBISO)=0.0D0
      VTOT=0.0D0
      DO 40 IBM=1,NMIX
      IBMOLD=MIXC(IBM)
      IF(IBMOLD.NE.0) VTOT=VTOT+XVOLM(IBMOLD)
   40 CONTINUE
      ALLOCATE(JJSO(NBISO),YLDSM(MY1,MY2),ITOD1(NBISO))
      ALLOCATE(TAUXFI(NBISO),NWT0(NGRP),SIGS(NGRP,NL,NBISO),
     > SS2D(NGRP,NGRP,NL,NBISO),XS(NGRP,NREA,NBISO))
      ALLOCATE(LXS(NREA))
      ALLOCATE(CHIRS(NGRP,NPRC),BETAR(NPRC),INVELS(NGRP))
      CHIRS(:NGRP,:NPRC)=0.0
      BETAR(:NPRC)=0.0
      INVELS(:NGRP)=0.0
      ALLOCATE(BETARB(NPRC),INVELSB(NGRP))
      ALLOCATE(DENS1(NBISO,NCAL),FACT(NBISO,NCAL))
      JPLIB=LCMLID(IPLIB,'ISOTOPESLIST',NBISO*NMIX)
*
      DO 180 IBM=1,NMIX
      IBMOLD=MIXC(IBM)
      IF(IBMOLD.EQ.0) GO TO 180
      IF(NISO(IBM).GT.MAXNIS) CALL XABORT('MCRLIB: MAXNIS OVERFLOW.')
      VOLMI2(IBM)=XVOLM(IBMOLD)
*----
*  RECOVER ITOD1(NBISO) INDICES.
*----
      ITOD1(:NBISO)=0
      DO 50 ISO=1,NBISO ! MPO file isotope
        DO KSO=1,NISO(IBM) ! user-selected isotope
          IF(NOMISO(ISO).EQ.HISO(IBM,KSO)) THEN
            ITOD1(ISO)=ITODO(IBM,KSO)
            GO TO 50
          ENDIF
        ENDDO
   50 CONTINUE
*----
*  COMPUTE THE NUMBER DENSITIES OF EACH ELEMENTARY CALCULATION.
*----
      DENS1(:NBISO,:NCAL)=0.0
      DENS3(:NBISO)=0.0
      DO ICAL=1,NCAL
        WEIGHT=TERP(ICAL,IBM)
        IF(WEIGHT.EQ.0.0) CYCLE
        DO ISO=1,NBISO
          LUSER=.FALSE.
          KSO1=0
          DO KSO=1,NISO(IBM) ! user-selected isotope
            IF(NOMISO(ISO).EQ.HISO(IBM,KSO)) THEN
              KSO1=KSO
              LUSER=(CONC(IBM,KSO1).NE.-99.99)
              GO TO 60
            ENDIF
          ENDDO
   60     IF(LUSER) THEN
            DENS1(ISO,ICAL)=CONC(IBM,KSO1)
            CYCLE
          ENDIF
          IF(.NOT.LISO(IBM)) CYCLE
          DENS1(ISO,ICAL)=DENS0(IBMOLD,ICAL,ISO)
        ENDDO
        DO ISO=1,NBISO
          DENS3(ISO)=DENS3(ISO)+WEIGHT*DENS1(ISO,ICAL)
        ENDDO
      ENDDO
      FACT(:NBISO,:NCAL)=1.0
      IF(.NOT.LPURE) THEN
        DO ICAL=1,NCAL
          IF(TERP(ICAL,IBM).EQ.0.0) CYCLE
          DO ISO=1,NBISO
            IF(DENS3(ISO).GT.DENS1(ISO,ICAL)*1.0E-9) THEN
              FACT(ISO,ICAL)=DENS1(ISO,ICAL)/DENS3(ISO)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
*----
*  INITIALIZE WORKING ARRAYS.
*----
      TAUXFI(:NBISO)=0.0
      NWT0(:NGRP)=0.0
      SIGS(:NGRP,:NL,:NBISO)=0.0
      SS2D(:NGRP,:NGRP,:NL,:NBISO)=0.0
      XS(:NGRP,:NREA,:NBISO)=0.0
      LXS(:NREA)=.FALSE.
      YLDSM(:MY1,:MY2)=0.0D0
*----
*  MAIN LOOP OVER ELEMENTARY CALCULATIONS
*----
      TEXT12='*MAC*RES'
      READ(TEXT12,'(2A4)') IHRES(1),IHRES(2)
      LSTRD=.FALSE.
      DO 80 ICAL=1,NCAL
      WEIGHT=TERP(ICAL,IBM)
      IF(WEIGHT.EQ.0.0) GO TO 80
*----
*  SELECT THE HDF5 GROUP CORRESPONDING TO ICAL AND IBMOLD
*----
      WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,6H/zone_,I0,1H/)')
     > TRIM(HEDIT),ICAL-1,IBMOLD-1
      NMGF=0
      IF(hdf5_group_exists(IPMPO,TRIM(RECNAM)//"yields")) THEN
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"yields/NMGF",NMGF)
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"yields/YIELDGROUP",
     >  IGYELD)
      ENDIF
      ALLOCATE(TAUXGF(NMGF,NBISO))
*----
*  RECOVER INFORMATION FROM caldir GROUP.
*----
      WRITE(RECNA2,'(A,9Hkinetics/)') TRIM(RECNAM)
      CALL hdf5_info(IPMPO,TRIM(RECNA2)//"LAMBDA",RANK,TYPE,NBYTE,DIMSR)
      NPRC=0
      IF(TYPE.NE.99) THEN
        NPRC=DIMSR(1)
        CALL hdf5_read_data(IPMPO,TRIM(RECNA2)//"LAMBDA",LAMB)
        CALL hdf5_read_data(IPMPO,TRIM(RECNA2)//"CHIDA",CHIRSB)
        CALL hdf5_read_data(IPMPO,TRIM(RECNA2)//"BETADA",BETARB)
        CALL hdf5_read_data(IPMPO,TRIM(RECNA2)//"INVELA",INVELSB)
      ENDIF
*----
*  RECOVER SPH FACTORS.
*----
      SPH(:NGRP)=1.0
      IF(HEQUI.NE.' ') THEN
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"LOCALVALUE",RVALO)
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"LOCALVALADDR",LOCAD)
        DO ILOC=1,NLOC
          IF((LOCTYP(ILOC).EQ.'EQUI').AND.(LOCKEY(ILOC).EQ.HEQUI))
     >    THEN
            IF(LOCAD(ILOC+1)-LOCAD(ILOC).NE.NGRP) THEN
              CALL XABORT('MCRLIB: INVALID NUMBER OF COMPONENTS FOR '
     >        //'SPH FACTORS')
            ENDIF
            DO IGR=1,NGRP
              SPH(IGR)=RVALO(LOCAD(ILOC)+IGR-1)
            ENDDO
          ENDIF
        ENDDO
        DEALLOCATE(LOCAD,RVALO)
      ENDIF
*----
*  RECOVER FLUXES.
*----
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ZONEFLUX",FLUXS)
      DO I=1,NGRP
        NWT0(I)=NWT0(I)+WEIGHT*FLUXS(I)/SPH(I)
      ENDDO
*----
*  RECOVER MICROSCOPIC CROSS SECTIONS.
*----
      DO ISO=1,NBISO
        FACT0=FACT(ISO,ICAL)
        DEN=DENS0(IBMOLD,ICAL,ISO)
        CALL MCRSX2(IPMPO,HEDIT,RECNAM,NREA,NGRP,NMGF,NL,ISO,NOMREA,
     1  NOMISO(ISO),DEN,FACT0,WEIGHT,SPH,FLUXS,IREAB,IREAF,LPURE,IGYELD,
     2  LXS,XS(1,1,ISO),SIGS(1,1,ISO),SS2D(1,1,1,ISO),TAUXFI(ISO),
     3  TAUXGF(1,ISO))
      ENDDO
      IF(NMGF.GT.0) DEALLOCATE(IGYELD)
      DEALLOCATE(FLUXS)
*
      IF(NPRC.GT.0) THEN
        DO IGR=1,NGRP
          INVELS(IGR)=INVELS(IGR)+SPH(IGR)*WEIGHT*INVELSB(IGR)
          DO IPRC=1,NPRC
            CHIRS(IGR,IPRC)=CHIRS(IGR,IPRC)+WEIGHT*CHIRSB(IGR,IPRC)
          ENDDO
        ENDDO
        DO IPRC=1,NPRC
          BETAR(IPRC)=BETAR(IPRC)+WEIGHT*BETARB(IPRC)
        ENDDO
      ENDIF
*----
*  COMPUTE DEPLETION CHAIN DATA
*----
      IF(NISOF*NISOP.GT.0) THEN
        IF(NMGF.EQ.0) CALL XABORT('MCRLIB: INVALID NMGF.')
        WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,6H/zone_,I0,1H/)')
     >  TRIM(HEDIT),ICAL-1,IBMOLD-1
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"yields/NISF",I0)
        IF(I0.NE.NISOF) CALL XABORT('MCRLIB: INVALID NISOF.')
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"yields/NISP",I0)
        IF(I0.NE.NISOP) CALL XABORT('MCRLIB: INVALID NISOP.')
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"yields/YIELD",YLDS2)
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"yields/ADDRY",
     >  DIMS_MPO)
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ADDRZI",ADDRZI)
        WRITE(RECNAM,'(8H/output/,A,6H/info/)') TRIM(HEDIT)
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ADDRISO",ADDRISO)
        NISOM=ADDRISO(ADDRZI+2)-ADDRISO(ADDRZI+1)
        DO IY1=1,NISOF
          ISO=0
          DO ISOM=1,NISOM
            IF(DIMS_MPO(ISOM).EQ.IY1) THEN
              ISO=ADDRISO(ADDRZI+1)+ISOM
              EXIT
            ENDIF
          ENDDO
          IF(ISO.EQ.0) CALL XABORT('MCRLIB: UNABLE TO FIND ISO.')
          DEN=0.0
          DO IGRC=1,NMGF
            DEN=DEN+TAUXGF(IGRC,ISO)
            DO IY2=1,NISOP
              YLDSM(IY1,IY2)=YLDSM(IY1,IY2)+WEIGHT*TAUXGF(IGRC,ISO)*
     >        YLDS2(IY1,IY2,IGRC)
              YLDS(IY1,IY2)=YLDS(IY1,IY2)+WEIGHT*TAUXGF(IGRC,ISO)*
     >        YLDS2(IY1,IY2,IGRC)*VOLMI2(IBM)/VTOT
            ENDDO
          ENDDO
          IF(DEN.EQ.0.0) CYCLE
          DO IY2=1,NISOP
            YLDSM(IY1,IY2)=YLDSM(IY1,IY2)/DEN
            YLDS(IY1,IY2)=YLDS(IY1,IY2)/DEN
          ENDDO
        ENDDO
        DEALLOCATE(ADDRISO,DIMS_MPO,YLDS2)
      ENDIF
      CALL hdf5_info(IPMPO,"/contents/isotopes/DECAYCONST",RANK,TYPE,
     1 NBYTE,DIMSR)
      IF(TYPE.NE.99) THEN
        CALL hdf5_read_data(IPMPO,"/contents/isotopes/DECAYCONST",
     >  DECAY2)
        DO ISO=1,NBISO
          DECAYC(ISO)=DECAYC(ISO)+WEIGHT*DECAY2(ISO)*VOLMI2(IBM)/VTOT
        ENDDO
        DEALLOCATE(DECAY2)
      ENDIF
      DEALLOCATE(TAUXGF)
   80 CONTINUE ! end of loop over elementary calculations.
*----
*  IDENTIFY SPECIAL FLUX EDITS
*----
      DO 100 IREA=1,NREA
        IF((NOMREA(IREA).EQ.'Total').or.
     &     (NOMREA(IREA).EQ.'Absorption').or.
     &     (NOMREA(IREA).EQ.'CaptureEnergyCapture').or.
     &     (NOMREA(IREA).EQ.'Diffusion').or.
     &     (NOMREA(IREA).EQ.'FissionEnergyFission').or.
     &     (NOMREA(IREA).EQ.'FissionSpectrum').or.
     &     (NOMREA(IREA).EQ.'NuFission').or.
     &     (NOMREA(IREA).EQ.'Scattering')) CYCLE
        DO 90 IED2=1,NED2
        IF(HVECT2(IED2).EQ.NOMREA(IREA)(:8)) GO TO 100
        IF(HVECT2(IED2).EQ.'NFTOT') GO TO 100
   90   CONTINUE
        NED2=NED2+1
        IF(NED2.GT.MAXREA) CALL XABORT('MCRLIB: MAXREA OVERFLOW.')
        IF(NOMREA(IREA).EQ.'Fission') THEN
          HVECT2(NED2)='NFTOT'
        ELSE
          HVECT2(NED2)=NOMREA(IREA)(:8)
        ENDIF
  100 CONTINUE
*----
*  SET FLAG LSTRD
*----
      LSTRD=.TRUE.
      DO IREA=1,NREA
        IF(NOMREA(IREA).EQ.'Leakage') THEN
          IF(LXS(IREA)) LSTRD=.FALSE.
          EXIT
        ENDIF
      ENDDO
*----
*  SET IADRY FOR MIXTURE IBMOLD
*----
      ALLOCATE(IADRY(NBISO))
      IADRY(:NBISO)=0
      DO ICAL=NCAL,1,-1
      IF(TERP(ICAL,IBM).EQ.0.0) CYCLE
        WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,6H/zone_,I0,1H/)')
     1  TRIM(HEDIT),ICAL-1,IBMOLD-1
        IF((hdf5_group_exists(IPMPO,TRIM(RECNAM)//"yields")).AND.
     1    (NISOP.GT.0)) THEN
          CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"yields/ADDRY",
     1    DIMS_MPO)
          CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ADDRZI",ADDRZI)
          WRITE(RECNAM,'(8H/output/,A,6H/info/)') TRIM(HEDIT)
          CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ADDRISO",ADDRISO)
          NISOM=ADDRISO(ADDRZI+2)-ADDRISO(ADDRZI+1)
          DO ISOM=1,NISOM
            ISO=ADDRISO(ADDRZI+1)+ISOM
            IADRY(ISO)=DIMS_MPO(ISOM)
          ENDDO
          DEALLOCATE(ADDRISO,DIMS_MPO)
        ENDIF
        EXIT
      ENDDO
*----
*  SAVE CROSS SECTIONS IN MICROLIB FOR MIXTURE IBM
*----
      ISTY1(:NBISO)=0
      JJSO(:NBISO)=0
      NBISO2I=NBISO2
      HRESID=' '
      DO ISO=1,NBISO
        READ(NOMISO(ISO),'(2A4)') INAME(:2)
        CALL SCRFND(MAXISO,NBISO2I,NBISO2,INAME,IBM,HRESID,HUSE2,
     1  HNAM2,IMIX2,JJSO(ISO))
        KPLIB=LCMDIL(JPLIB,JJSO(ISO)) ! step up isot JJSO(ISO)
        CALL MCRISO(KPLIB,NREA,NGRP,NL,NPRC,NOMREA,NWT0,XS(1,1,ISO),
     1  SIGS(1,1,ISO),SS2D(1,1,1,ISO),TAUXFI(ISO),LXS,LAMB,CHIRS,
     2  BETAR,INVELS,INAME,LSTRD,LPURE,ILUPS,ITRANC,IFISS)
        IF(MY1*MY2.GT.0) CALL MCRNDF(IMPX,NBISO,ISO,IBM,NOMISO,KPLIB,
     1  MY1,MY2,YLDSM,IADRY,ISTY1(ISO))
      ENDDO
      DEALLOCATE(IADRY)
*----
*  SET NUMBER DENSITIES AND VOLUMES IN OUTPUT MICROLIB
*----
      IF(LRES) THEN
*       -- Number densities are left unchanged except if they are
*       -- listed in HISO array.
        DO 110 KSO=1,NISO(IBM) ! user-selected isotope
          DO JSO=1,NBISO2 ! microlib isotope
            IF(IMIX2(JSO).NE.IBM) CYCLE
            WRITE(TEXT8,'(2A4)') HUSE2(1,JSO),HUSE2(2,JSO)
            IF(HISO(IBM,KSO).EQ.TEXT8) THEN
              ITOD2(JSO)=ITODO(IBM,KSO)
              IF(CONC(IBM,KSO).EQ.-99.99) THEN
*               -- Only number densities of isotopes set with "MICR" and
*               -- "*" keywords are interpolated
                DENS2(JSO)=0.0
                DO ISO=1,NBISO ! MPO file isotope
                  IF(JJSO(ISO).EQ.JSO) DENS2(JSO)=DENS2(JSO)+DENS3(ISO)
                ENDDO
              ELSE IF(CONC(IBM,KSO).NE.-99.99) THEN
*               -- Number densities of isotopes set with "MICR" and
*               -- fixed value are forced to this value
                DENS2(JSO)=CONC(IBM,KSO)
              ENDIF
              GO TO 110
            ENDIF
          ENDDO
          WRITE(HSMG,'(31HMCRLIB: UNABLE TO FIND ISOTOPE ,A8,6H IN MI,
     1    5HXTURE,I8,1H.)') HISO(IBM,KSO),IBM
          CALL XABORT(HSMG)
  110   CONTINUE
      ELSE
*       -- Number densities are interpolated or not according to
*       -- ALL/ONLY option
        DO JSO=1,NBISO2 ! microlib isotope
          WRITE(TEXT8,'(2A4)') HUSE2(1,JSO),HUSE2(2,JSO)
          IF(IBM.EQ.IMIX2(JSO)) THEN
            DO ISO=1,NBISO ! MPO file isotope
              IF(NOMISO(ISO).EQ.TEXT8) THEN
                DENS2(JSO)=0.0
                VOL2(JSO)=0.0
                CYCLE
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        DO 130 ISO=1,NBISO ! MPO file isotope
        IF(.NOT.LISO(IBM)) THEN
*         --ONLY option
          DO KSO=1,NISO(IBM) ! user-selected isotope
            IF(NOMISO(ISO).EQ.HISO(IBM,KSO)) GO TO 120
          ENDDO
          GO TO 130
        ENDIF
  120   JSO=JJSO(ISO)
        IF(JSO.GT.0) THEN
          ITOD2(JSO)=ITOD1(ISO)
          ISTY2(JSO)=ISTY1(ISO)
          DENS2(JSO)=DENS2(JSO)+DENS3(ISO)
          VOL2(JSO)=VOL2(JSO)+XVOLM(IBMOLD)
        ENDIF
  130   CONTINUE
      ENDIF
*----
*  SET PIFI INFORMATION
*----
      ALLOCATE(IMICR(NBISO))
      IMICR(:NBISO)=0
      NBS1=0
      DO 140 JSO=1,NBISO2 ! microlib isotope
      IF(IMIX2(JSO).EQ.IBM) THEN
        NBS1=NBS1+1
        IF(NBS1.GT.NBISO) CALL XABORT('MCRLIB: NBISO OVERFLOW.')
        IMICR(NBS1)=JSO
      ENDIF
  140 CONTINUE
      DO 170 ISO=1,NBS1 ! MPO file isotope
      JSO=IMICR(ISO)
      KPLIB=LCMDIL(JPLIB,JSO) ! step up isot JSO
      CALL LCMLEN(KPLIB,'PYIELD',LMY1,ITYLCM)
      IF(LMY1.GT.0) THEN
        ALLOCATE(HPYNAM(LMY1),IPYMIX(LMY1),IPIFI(LMY1))
        IPIFI(:LMY1)=0
        CALL LCMGTC(KPLIB,'PYNAM',8,LMY1,HPYNAM)
        CALL LCMGET(KPLIB,'PYMIX',IPYMIX)
        DO 160 IY1=1,LMY1
        IF(HPYNAM(IY1).NE.' ') THEN
          DO 150 KSO=1,NBS1
          LSO=IMICR(KSO)
          WRITE(TEXT8,'(2A4)') HUSE2(:2,LSO)
          IF((HPYNAM(IY1).EQ.TEXT8).AND.(IPYMIX(IY1).EQ.IMIX2(LSO)))
     >    THEN
            IPIFI(IY1)=LSO
            GO TO 160
          ENDIF
  150     CONTINUE
          IF(IPIFI(IY1).EQ.0) THEN
            WRITE(HSMG,'(40HMCRLIB: FAILURE TO FIND FISSILE ISOTOPE ,
     1      A12,25H AMONG MICROLIB ISOTOPES.)') HPYNAM(IY1)
            CALL XABORT(HSMG)
          ENDIF
        ENDIF
  160   CONTINUE
        CALL LCMPUT(KPLIB,'PIFI',LMY1,1,IPIFI)
        DEALLOCATE(IPIFI,IPYMIX,HPYNAM)
      ENDIF
  170 CONTINUE
      DEALLOCATE(IMICR)
  180 CONTINUE ! end of loop over microlib mixtures.
*----
*  RELEASE MEMORY
*----
      DEALLOCATE(FACT,DENS1)
      IF(NPRC.GT.0) DEALLOCATE(INVELSB,BETARB,CHIRSB,INVELS,BETAR,
     > CHIRS,LAMB)
      DEALLOCATE(LXS,XS,SS2D,SIGS,NWT0,TAUXFI)
      DEALLOCATE(ITOD1,YLDSM)
      DEALLOCATE(JJSO,DENS0,XVOLM,NOMISO)
*----
*  MICROLIB FINALIZATION
*----
      IF(.NOT.LRES) THEN
        ISTATE(:NSTATE)=0
        ISTATE(1)=NMIX
        ISTATE(2)=NBISO2
        ISTATE(3)=NGRP
        ISTATE(4)=NL
        ISTATE(5)=ITRANC
        ISTATE(7)=1
        IF(ITER.EQ.3) ISTATE(12)=NMIX
        ISTATE(13)=NED2
        ISTATE(14)=NMIX
        ISTATE(18)=1
        ISTATE(19)=NPRC
        ISTATE(20)=MY1
        ISTATE(22)=MAXISO/NMIX
        IF(NBISO2.EQ.0) CALL XABORT('MCRLIB: NBISO2=0.')
        CALL LCMPUT(IPLIB,'STATE-VECTOR',NSTATE,1,ISTATE)
        CALL LCMPUT(IPLIB,'MIXTURESVOL',NMIX,2,VOLMI2)
        CALL LCMPUT(IPLIB,'ISOTOPESUSED',3*NBISO2,3,HUSE2)
        CALL LCMPUT(IPLIB,'ISOTOPERNAME',3*NBISO2,3,HNAM2)
        CALL LCMPUT(IPLIB,'ISOTOPESDENS',NBISO2,2,DENS2)
        CALL LCMPUT(IPLIB,'ISOTOPESMIX',NBISO2,1,IMIX2)
        CALL LCMPUT(IPLIB,'ISOTOPESVOL',NBISO2,2,VOL2)
        IF(NED2.GT.0) CALL LCMPTC(IPLIB,'ADDXSNAME-P0',8,NED2,HVECT2)
        CALL LCMPUT(IPLIB,'ISOTOPESTODO',NBISO2,1,ITOD2)
        CALL LCMPUT(IPLIB,'ISOTOPESTYPE',NBISO2,1,ISTY2)
      ELSE IF(LRES.AND.(NBISO.GT.0)) THEN
        CALL LCMPUT(IPLIB,'ISOTOPESDENS',NBISO2,2,DENS2)
        CALL LCMPUT(IPLIB,'ISOTOPESVOL',NBISO2,2,VOL2)
      ENDIF
      IF(IMPX.GT.5) CALL LCMLIB(IPLIB)
*----
*  COMPUTE THE MACROSCOPIC X-SECTIONS
*----
      IF((ITER.NE.0).AND.(ITER.NE.3)) GO TO 280
      CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
      MAXMIX=ISTATE(1)
      IF(MAXMIX.NE.NMIX) CALL XABORT('MCRLIB: INVALID NMIX.')
      NBISO=ISTATE(2)
      ALLOCATE(MASK(MAXMIX),MASKL(NGRP))
      ALLOCATE(ISONA(3*NBISO),ISOMI(NBISO),DENIS(NBISO))
      CALL LCMGET(IPLIB,'ISOTOPESUSED',ISONA)
      CALL LCMGET(IPLIB,'ISOTOPESMIX',ISOMI)
      CALL LCMGET(IPLIB,'ISOTOPESDENS',DENIS)
      MASK(:MAXMIX)=.TRUE.
      MASKL(:NGRP)=.TRUE.
      ITSTMP=0
      TMPDAY(1)=0.0
      TMPDAY(2)=0.0
      TMPDAY(3)=0.0
      CALL LCMLEN(IPLIB,'MACROLIB',ILONG,ITYLCM)
      IF(ILONG.NE.0) CALL LCMDEL(IPLIB,'MACROLIB')
      CALL LIBMIX(IPLIB,MAXMIX,NGRP,NBISO,ISONA,ISOMI,DENIS,MASK,MASKL,
     > ITSTMP,TMPDAY)
      DEALLOCATE(MASKL,MASK)
      DEALLOCATE(DENIS,ISOMI,ISONA)
*----
*  INCLUDE LEAKAGE IN THE MACROLIB (USED ONLY FOR NON-REGRESSION TESTS)
*----
      IF(B2.NE.0.0) THEN
        IF(IMPX.GT.0) WRITE(IOUT,'(/31H MCRLIB: INCLUDE LEAKAGE IN THE,
     >  14H MACROLIB (B2=,1P,E12.5,2H).)') B2
        CALL LCMSIX(IPLIB,'MACROLIB',1)
        JPLIB=LCMGID(IPLIB,'GROUP')
        ALLOCATE(GAR1(NMIX),GAR2(NMIX))
        DO 270 IGR=1,NGRP
          KPLIB=LCMGIL(JPLIB,IGR)
          CALL LCMGET(KPLIB,'NTOT0',GAR1)
          CALL LCMGET(KPLIB,'DIFF',GAR2)
          DO 260 IBM=1,NMIX
            IF(MIXC(IBM).NE.0) GAR1(IBM)=GAR1(IBM)+B2*GAR2(IBM)
  260     CONTINUE
          CALL LCMPUT(KPLIB,'NTOT0',NMIX,2,GAR1)
  270   CONTINUE
        DEALLOCATE(GAR2,GAR1)
        CALL LCMSIX(IPLIB,' ',2)
      ENDIF
*----
*  PROCESS ADF and physical albedos (if required)
*----
  280 LALBG=.TRUE.
      IDF=0
      IF(NALBP.GT.0) THEN
        WRITE(RECNAM,'(8H/output/,A,16H/statept_0/flux/)') TRIM(HEDIT)
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"NALBP",NALBP2)
        IF(NALBP2.NE.NALBP) CALL XABORT('MCRLIB: INVALID NALBP.')
        CALL hdf5_info(IPMPO,TRIM(RECNAM)//"ALBEDOGxG",RANK,TYPE,NBYTE,
     &  DIMSR)
        IF(TYPE.NE.99) LALBG=.FALSE.
      ENDIF
      CALL LCMSIX(IPLIB,'MACROLIB',1)
      CALL MCRAGF(IPLIB,IPMPO,IACCS,NMIL,NMIX,NGRP,NALBP,LALBG,IMPX,
     1 NCAL,TERP,MIXC,NSURFD,HEDIT,VOLMI2,IDF)
      IF(NSURFD.GT.0) THEN
        CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
        ISTATE(12)=IDF ! ADF information
        CALL LCMPUT(IPLIB,'STATE-VECTOR',NSTATE,1,ISTATE)
      ENDIF
      CALL LCMSIX(IPLIB,' ',2)
      IF(NSURFD.GT.0) THEN
        CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
        ISTATE(24)=IDF ! ADF information
        CALL LCMPUT(IPLIB,'STATE-VECTOR',NSTATE,1,ISTATE)
      ENDIF
      IACCS=1
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SPH,FLUX,VOLMI2,VOL2,DENS3,DENS2)
      DEALLOCATE(HNAM2,HUSE2,ISTY2,ISTY1,ITOD2,IMIX2)
      RETURN
      END
