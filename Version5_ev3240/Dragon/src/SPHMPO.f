*DECK SPHMPO
      SUBROUTINE SPHMPO(IPMPO,IPMAC,ICAL,IMPX,HEQUI,HMASL,NMIL,NALBP,
     1 NGROUP,HEDIT,VOSAP,ILUPS,SPH,B2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Extract a Macrolib corresponding to an elementary calculation in a
* MPO file
*
*Copyright:
* Copyright (C) 2022 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPMPO   pointer to the MPO file.
* ICAL    index of the elementary calculation being considered.
* IMPX    print parameter (equal to zero for no print).
* HEQUI   keyword of SPH-factor set to be recovered.
* HMASL   keyword of MASL data set to be recovered.
* NMIL    number of mixtures in the elementary calculation.
* NALBP   number of physical albedos per energy group.
* NGROUP  number of energy groups in the elementary calculation.
* HEDIT   name of output group for a (multigroup mesh, output geometry)
*         couple (generally equal to 'output_0').
* VOSAP   mixture volumes in the MPO file.
* ILUPS   up-scattering removing flag (=1 to remove up-scattering from
*         output cross-sections).
* B2      imposed buckling.
*
*Parameters: output
* IPMAC   pointer to the Macrolib (L_MACROLIB signature).
* SPH     SPH-factor set extracted from the MPO file.
* B2      buckling recovered from MPO file.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMPO,IPMAC
      INTEGER ICAL,IMPX,NMIL,NALBP,NGROUP,ILUPS
      REAL VOSAP(NMIL),SPH(NMIL+NALBP,NGROUP),B2
      CHARACTER(LEN=80) HEQUI,HMASL
      CHARACTER(LEN=12) HEDIT
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::NSTATE=40
      INTEGER ISTATE(NSTATE)
      INTEGER I,NADDRXS,NREA,NISO,NISOM,ADDRZX,ADDRZI,IBM,ISO,ISOM,
     & NBISO,NBYTE,RANK,TYPE,DIMSR(5),NL,ILOC,NLOC,IPROF,JOFS,NL1,
     & NL2,IPRC,NPRC,IOF,IOF1,IOF2,IGR,JGR,ITRANC,NED,IFISS,IGMAX,
     & IGMIN,IPOSDE,IL,IREA,NSURFD,IDF,NALBP2,TYPE2,TYPE4
      REAL FLOTT,DEN,ZIL,FF,CSCAT
      LOGICAL LSPH,LMASL,LSTRD,LDIFF,LHFACT
      CHARACTER RECNAM*80,RECNA2*80,TEXT12*12,CM*2,HSMG*131
      TYPE(C_PTR) JPMAC,KPMAC
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: REACTION,ISOTOPE,IDATAP,
     & ADDRISO,LOCAD,FAG,ADR
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: ADDRXS
      REAL, ALLOCATABLE, DIMENSION(:) :: FLUXS,CONCEN,RDATAX,RVALO,
     & FMASL,LAMBDAD,BETADF,INVEL,LAMB,VREAL,VFLUX,SURF,LG
      REAL, ALLOCATABLE, DIMENSION(:,:) :: NWT0,EFACT,CHID,SIGS0,TOTAL,
     & DIFF,BETAR,INVELS,SURFLX,ALBP
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: XS,SIGS,CHIRS,ALBP_ERM
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: SS2D
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HEDI,HADF
      CHARACTER(LEN=24), ALLOCATABLE, DIMENSION(:) :: TEXT24,NOMREA,
     & NOMISO
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LXS
      CHARACTER(LEN=80), ALLOCATABLE, DIMENSION(:) :: LOCTYP,LOCKEY
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPOS,NJJM,IJJM
      REAL, ALLOCATABLE, DIMENSION(:) :: SCAT,GAR
*----
*  SCRATCH STORAGE ALLOCATION
*   SIGS0    P0 scattering cross sections.
*   TOTAL    total cross sections.
*   DIFF     diffusion coefficients.
*   FMASL    heavy element mass.
*----
      ALLOCATE(IPOS(NMIL),NJJM(NMIL),IJJM(NMIL))
      ALLOCATE(SIGS0(NMIL,NGROUP),TOTAL(NMIL,NGROUP),DIFF(NMIL,NGROUP),
     & FMASL(NMIL))
      SIGS0(:NMIL,:NGROUP)=0.0
      TOTAL(:NMIL,:NGROUP)=0.0
      DIFF(:NMIL,:NGROUP)=0.0
      FMASL(:NMIL)=0.0
*----
*  RECOVER MPO FILE CHARACTERISTICS
*----
      WRITE(RECNAM,'(8H/output/,A,6H/info/)') TRIM(HEDIT)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"NADDRXS",NADDRXS)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"NREA",NREA)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"NISO",NISO)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ADDRXS",ADDRXS)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"REACTION",REACTION)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ISOTOPE",ISOTOPE)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"TRANSPROFILE",IDATAP)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ADDRISO",ADDRISO)
      NBISO=ADDRISO(SIZE(ADDRISO,1))
      IF(NBISO.EQ.0) CALL XABORT('SPHMPO: NO CROSS SECTIONS.')
      ALLOCATE(NOMREA(NREA),NOMISO(NBISO))
      CALL hdf5_read_data(IPMPO,"contents/reactions/REACTIONAME",TEXT24)
      DO I=1,NREA
        NOMREA(I)=TEXT24(REACTION(I)+1)
      ENDDO
      DEALLOCATE(TEXT24)
      CALL hdf5_read_data(IPMPO,"contents/isotopes/ISOTOPENAME",TEXT24)
      DO I=1,NBISO
        NOMISO(I)=TEXT24(ISOTOPE(I)+1)
      ENDDO
      DEALLOCATE(TEXT24)
      IF(IMPX.GT.1) THEN
        WRITE(6,'(/24H SPHMPO: reaction names:)')
        DO I=1,NREA
          WRITE(6,'(5X,7HNOMREA(,I3,2H)=,A)') I,TRIM(NOMREA(I))
        ENDDO
        WRITE(6,'(/23H SPHMPO: isotope names:)')
        DO I=1,NBISO
          WRITE(6,'(5X,7HNOMISO(,I3,2H)=,A)') I,TRIM(NOMISO(I))
        ENDDO
        WRITE(6,'(/34H SPHMPO: number of energy groups =,I4)') NGROUP
        WRITE(6,'(30H SPHMPO: number of mixtures  =,I4)') NMIL
        WRITE(6,'(30H SPHMPO: number of reactions =,I4)') NREA
        WRITE(6,'(30H SPHMPO: number of isotopes  =,I4)') NBISO
      ENDIF
      WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0)') TRIM(HEDIT),ICAL-1
      IF(.not.hdf5_group_exists(IPMPO,TRIM(RECNAM))) THEN
        WRITE(HSMG,'(38HSPHMPO: missing elementary calculation,I5,
     &  10H in group ,A,1H.)') ICAL,TRIM(HEDIT)
        CALL XABORT(HSMG)
      ENDIF
*----
*  SET EQUIVALENCE AND HEAVY DENSITY FLAGS.
*----
      LSPH=.FALSE.
      LMASL=.FALSE.
      NLOC=0
      IF(hdf5_group_exists(IPMPO,"/local_values/")) THEN
        CALL hdf5_read_data(IPMPO,"local_values/LOCVALTYPE",LOCTYP)
        CALL hdf5_read_data(IPMPO,"local_values/LOCVALNAME",LOCKEY)
        NLOC=SIZE(LOCTYP,1)
        DO ILOC=1,NLOC
          LSPH=LSPH.OR.((LOCTYP(ILOC).EQ.'EQUI').AND.
     &                  (LOCKEY(ILOC).EQ.HEQUI))
          LMASL=LMASL.OR.((LOCTYP(ILOC).EQ.'HEAVY_METAL_DENSITY').AND.
     &                    (LOCKEY(ILOC).EQ.HMASL))
        ENDDO
      ENDIF
*----
*  FIND SCATTERING ANISOTROPY.
*----
      NL=0
      DO I=1,NADDRXS
        DO ISO=1,NISO
          NL=MAX(NL,ADDRXS(NREA+1,ISO,I))
          NL=MAX(NL,ADDRXS(NREA+2,ISO,I))
        ENDDO
      ENDDO
      IF(IMPX.GT.1) THEN
        WRITE(6,'(37H SPHMPO: number of legendre orders  =,I4)') NL
      ENDIF
*----
*  RECOVER GENERAL INFORMATION
*----
      LSTRD=(B2.EQ.0.0)
      WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,8H/addons/)')
     & TRIM(HEDIT),ICAL-1
      CALL hdf5_info(IPMPO,TRIM(RECNAM)//"KEFF",RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) THEN
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"KEFF",FLOTT)
        CALL LCMPUT(IPMAC,'K-EFFECTIVE',1,2,FLOTT)
        IF(IMPX.GT.1) THEN
          WRITE(6,'(22H SPHMPO: K-EFFECTIVE =,1P,E13.6)') FLOTT
        ENDIF
      ENDIF
      CALL hdf5_info(IPMPO,TRIM(RECNAM)//"KINF",RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) THEN
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"KINF",FLOTT)
        CALL LCMPUT(IPMAC,'K-INFINITY',1,2,FLOTT)
        IF(IMPX.GT.1) THEN
          WRITE(6,'(21H SPHMPO: K-INFINITY =,1P,E13.6)') FLOTT
        ENDIF
      ENDIF
      CALL hdf5_info(IPMPO,TRIM(RECNAM)//"B2",RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) THEN
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"B2",B2)
        LSTRD=(B2.EQ.0.0)
        CALL LCMPUT(IPMAC,'B2  B1HOM',1,2,B2)
        IF(IMPX.GT.1) THEN
          WRITE(6,'(13H SPHMPO: B2 =,1P,E14.6)') B2
        ENDIF
      ENDIF
*----
*  RECOVER DISCONTINUITY FACTOR INFORMATION
*----
      IDF=0
      WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,6H/flux/)')
     & TRIM(HEDIT),ICAL-1
      CALL hdf5_info(IPMPO,TRIM(RECNAM)//"NSURF",RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) THEN
        CALL hdf5_info(IPMPO,TRIM(RECNAM)//"SURFFLUX",RANK,TYPE2,NBYTE,
     &  DIMSR)
        CALL hdf5_info(IPMPO,TRIM(RECNAM)//"SURFFLUXGxG",RANK,TYPE4,
     &  NBYTE,DIMSR)
        IF(TYPE2.NE.99) THEN
          IDF=2 ! boundary flux information
          CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"SURFFLUX",SURFLX)
        ELSE IF(TYPE4.NE.99) THEN
          IDF=4 ! matrix discontinuity factor information
          CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"SURFFLUXGxG",SURFLX)
        ELSE
          CALL hdf5_list(IPMPO,TRIM(RECNAM))
          CALL XABORT('SPHMPO: UNABLE TO SET TYPE OF DF.')
        ENDIF
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"NSURF",NSURFD)
        CALL hdf5_info(IPMPO,TRIM(RECNAM)//"SURF",RANK,TYPE,NBYTE,DIMSR)
        IF(TYPE.NE.99) THEN
          CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"SURF",SURF)
          IF(DIMSR(1)/NSURFD.NE.NMIL) THEN
            WRITE(HSMG,'(24HSPHMPO: INVALID LENGTH (,I5,11H) FOR SURF ,
     &      14HGROUP. LENGTH=,I5,10H EXPECTED.)') DIMSR(1),NMIL*NSURFD
            CALL XABORT(HSMG)
          ENDIF
        ELSE
*         temporary.....
          CALL hdf5_read_data(IPMPO,"/geometry/geometry_0/COORDINATE",
     &    LG)
          ALLOCATE(SURF(NMIL*NSURFD))
          SURF(:NMIL*NSURFD)=LG(2)
          DEALLOCATE(LG)
        ENDIF
        CALL LCMSIX(IPMAC,'ADF',1)
        CALL LCMPUT(IPMAC,'NTYPE',1,1,NSURFD)
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"TOTALFLUX",VFLUX)
        DO IGR=1,NGROUP
          DO IBM=1,NMIL
            IOF=(IGR-1)*NMIL+IBM
            VFLUX(IOF)=VFLUX(IOF)/VOSAP(IBM)
          ENDDO
        ENDDO
        CALL LCMPUT(IPMAC,'AVG_FLUX',NMIL*NGROUP,2,VFLUX)
        ALLOCATE(HADF(NSURFD))
        IF(IDF.EQ.2) THEN
          ALLOCATE(VREAL(NMIL*NGROUP))
          DO I=1,NSURFD
            WRITE(HADF(I),'(3HFD_,I5.5)') I
            DO IGR=1,NGROUP
              DO IBM=1,NMIL
                IOF=(IGR-1)*NMIL+IBM
                VREAL(IOF)=SURFLX(I,IOF)/SURF((I-1)*NMIL+IBM)
              ENDDO
            ENDDO
            CALL LCMPUT(IPMAC,HADF(I),NMIL*NGROUP,2,VREAL)
          ENDDO
          DEALLOCATE(VREAL)
        ELSE IF(IDF.EQ.4) THEN
          ALLOCATE(VREAL(NMIL*NGROUP*NGROUP))
          DO I=1,NSURFD
            WRITE(HADF(I),'(3HFD_,I5.5)') I
            DO JGR=1,NGROUP
              DO IGR=1,NGROUP
                DO IBM=1,NMIL
                  IOF1=(IGR-1)*NMIL+IBM
                  IOF2=((JGR-1)*NGROUP+IGR-1)*NMIL+IBM
                  VREAL(IOF2)=SURFLX(I,IOF2)/SURF((I-1)*NMIL+IBM)/
     &            VFLUX(IOF1)
                ENDDO
              ENDDO
            ENDDO
            CALL LCMPUT(IPMAC,HADF(I),NMIL*NGROUP*NGROUP,2,VREAL)
          ENDDO
          DEALLOCATE(VREAL)
        ENDIF
        CALL LCMPTC(IPMAC,'HADF',8,NSURFD,HADF)
        DEALLOCATE(VFLUX,HADF,SURF,SURFLX)
        CALL LCMSIX(IPMAC,' ',2)
      ENDIF
*----
*  RECOVER ALBEDO INFORMATION
*----
      CALL hdf5_info(IPMPO,TRIM(RECNAM)//"ALBEDOG",RANK,TYPE,NBYTE,
     & DIMSR)
      IF(TYPE.NE.99) THEN
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"NALBP",NALBP2)
        IF(NALBP2.NE.NALBP) CALL XABORT('SPHMPO: INVALID NALBP(1).')
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ALBEDOG",ALBP)
        CALL LCMPUT(IPMAC,'ALBEDO',NALBP*NGROUP,2,ALBP)
        DEALLOCATE(ALBP)        
      ENDIF
      CALL hdf5_info(IPMPO,TRIM(RECNAM)//"ALBEDOGxG",RANK,TYPE,NBYTE,
     & DIMSR)
      IF(TYPE.NE.99) THEN
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"NALBP",NALBP2)
        IF(NALBP2.NE.NALBP) CALL XABORT('SPHMPO: INVALID NALBP(2).')
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ALBEDOGxG",ALBP_ERM)
        CALL LCMPUT(IPMAC,'ALBEDO',NALBP*NGROUP*NGROUP,2,ALBP_ERM)
        DEALLOCATE(ALBP_ERM)        
      ENDIF
*----
*  ALLOCATE MACROLIB WORKING ARRAYS.
*----
      ALLOCATE(LXS(NREA),NWT0(NMIL,NGROUP),EFACT(NMIL,NGROUP),
     1 SIGS(NMIL,NGROUP,NL),SS2D(NMIL,NGROUP,NGROUP,NL),
     2 XS(NMIL,NGROUP,NREA),FAG(NGROUP),ADR(NGROUP))
      NWT0(:NMIL,:NGROUP)=0.0
      EFACT(:NMIL,:NGROUP)=0.0
      SIGS(:NMIL,:NGROUP,:NL)=0.0
      SS2D(:NMIL,:NGROUP,:NGROUP,:NL)=0.0
      XS(:NMIL,:NGROUP,:NREA)=0.0
      LXS(:NREA)=.FALSE.
      LDIFF=.FALSE.
*----
*  ALLOCATE DELAYED NEUTRON WORKING ARRAYS.
*----
      WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,
     & 24H/zone_0/kinetics/LAMBDAD)') TRIM(HEDIT),ICAL-1
      NPRC=0
      CALL hdf5_info(IPMPO,TRIM(RECNAM),RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) NPRC=DIMSR(1)
      IF(IMPX.GT.1) THEN
        WRITE(6,'(37H SPHMPO: number of precursor groups =,I4)') NPRC
      ENDIF
      ALLOCATE(LAMB(NPRC),CHIRS(NGROUP,NPRC,NMIL),BETAR(NPRC,NMIL),
     & INVELS(NGROUP,NMIL))
      LAMB(:NPRC)=0.0
      CHIRS(:NGROUP,:NPRC,:NMIL)=0.0
      BETAR(:NPRC,:NMIL)=0.0
      INVELS(:NGROUP,:NMIL)=0.0
*----
*  LOOP OVER MPO MIXTURES
*----
      DO IBM=1,NMIL
        WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,6H/zone_,I0,1H/)')
     &  TRIM(HEDIT),ICAL-1,IBM-1
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ZONEFLUX",FLUXS)
        DO I=1,NGROUP
          NWT0(IBM,I)=NWT0(IBM,I)+FLUXS(I)
        ENDDO
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ADDRZI",ADDRZI)
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ADDRZX",ADDRZX)
        NISOM=ADDRISO(ADDRZI+2)-ADDRISO(ADDRZI+1)
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"CONCENTRATION",CONCEN)
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"CROSSECTION",RDATAX)
        DO ISOM=1,NISOM
          ! loop over the isotopes present in the mix
          NL1=ADDRXS(NREA+1,ISOM,ADDRZX+1)
          NL2=ADDRXS(NREA+2,ISOM,ADDRZX+1)
          IF((NL1.GT.NL).OR.(NL2.GT.NL)) THEN
            CALL XABORT('SPHMPO: NL OVERFLOW.')
          ENDIF
          DEN=CONCEN(ISOM)
          IF(DEN.NE.0.0) THEN
            IF(IMPX.GT.3) THEN
              ISO=ADDRISO(ADDRZI+1)+ISOM
              WRITE(6,'(10H  mixture=,I5,15H concentration(,A,2H)=,1P,
     &        E12.4)') IBM,TRIM(NOMISO(ISO)),DEN
            ENDIF
            DO IREA=1,NREA
              IOF=ADDRXS(IREA,ISOM,ADDRZX+1)
              IF(IOF.LT.0) CYCLE
              IF(NOMREA(IREA).EQ.'Diffusion') THEN
                DO IL=1,NL1
                  DO IGR=1,NGROUP
                    FLOTT=DEN*RDATAX(IOF+(IL-1)*NGROUP+IGR)
                    SIGS(IBM,IGR,IL)=SIGS(IBM,IGR,IL)+FLOTT
                    LXS(IREA)=LXS(IREA).OR.(FLOTT.NE.0.0)
                  ENDDO
                ENDDO
              ELSE IF(NOMREA(IREA).EQ.'Scattering') THEN
                IPROF=ADDRXS(NREA+3,ISOM,ADDRZX+1)
                DO IGR=1,NGROUP
                  FAG(IGR)=IDATAP(IPROF+IGR)+1
                  ADR(IGR)=IDATAP(IPROF+NGROUP+IGR)
                ENDDO
                ADR(NGROUP+1)=IDATAP(IPROF+1+2*NGROUP)
                JOFS=0
                DO IL=1,NL2
                  ZIL=REAL(2*IL-1)
                  DO IGR=1,NGROUP
                    DO JGR=FAG(IGR),FAG(IGR)+(ADR(IGR+1)-ADR(IGR))-1
                      IF(JGR.GT.NGROUP) CALL XABORT('SPHMPO: SS2D OVER'
     &                //'FLOW.')
                      FLOTT=DEN*RDATAX(IOF+JOFS+1)/ZIL
                      SS2D(IBM,JGR,IGR,IL)=SS2D(IBM,JGR,IGR,IL)+FLOTT ! JGR <-- IGR
                      JOFS=JOFS+1
                      LXS(IREA)=LXS(IREA).OR.(FLOTT.NE.0.0)
                    ENDDO
                  ENDDO
                ENDDO
              ELSE IF(NOMREA(IREA).EQ.'FissionSpectrum') THEN
                DO IGR=1,NGROUP
                  XS(IBM,IGR,IREA)=RDATAX(IOF+IGR)
                  LXS(IREA)=LXS(IREA).OR.(RDATAX(IOF+IGR).NE.0.0)
                ENDDO
              ELSE
                DO IGR=1,NGROUP
                  XS(IBM,IGR,IREA)=XS(IBM,IGR,IREA)+DEN*RDATAX(IOF+IGR)
                  LXS(IREA)=LXS(IREA).OR.(DEN*RDATAX(IOF+IGR).NE.0.0)
                ENDDO
              ENDIF
            ENDDO ! end of loop over reactions
          ENDIF
        ENDDO ! end of loop over isotopes
        DEALLOCATE(RDATAX,CONCEN,FLUXS)
*
        IF(hdf5_group_exists(IPMPO,TRIM(RECNAM)//"kinetics")) THEN
          WRITE(RECNA2,'(A,9Hkinetics/)') TRIM(RECNAM)
          CALL hdf5_list(IPMPO,TRIM(RECNA2))
          CALL hdf5_read_data(IPMPO,TRIM(RECNA2)//"LAMBDAD",LAMBDAD)
          IF(SIZE(LAMBDAD,1).NE.NPRC) CALL XABORT('SPHMPO: WRONG NPRC.')
          CALL hdf5_read_data(IPMPO,TRIM(RECNA2)//"CHID",CHID)
          CALL hdf5_read_data(IPMPO,TRIM(RECNA2)//"BETADF",BETADF)
          CALL hdf5_read_data(IPMPO,TRIM(RECNA2)//"INVERSESPEED",INVEL)
          CHIRS(:NGROUP,:NPRC,IBM)=CHID(:NGROUP,:NPRC)
          BETAR(:NPRC,IBM)=BETADF(:NPRC)
          INVELS(:NGROUP,IBM)=INVEL(:NGROUP)
          DEALLOCATE(INVEL,BETADF,CHID,LAMBDAD)
        ENDIF
*
*       UP-SCATTERING CORRECTION OF THE MACROLIB.
        IF(ILUPS.EQ.1) THEN
          DO JGR=2,NGROUP
            DO IGR=1,JGR-1 ! IGR < JGR
              FF=NWT0(IBM,JGR)/NWT0(IBM,IGR)
              CSCAT=SS2D(IBM,IGR,JGR,1) ! IGR < JGR
              DO IL=1,NL
                CSCAT=SS2D(IBM,IGR,JGR,IL)
                SIGS(IBM,IGR,IL)=SIGS(IBM,IGR,IL)-CSCAT*FF
                SIGS(IBM,JGR,IL)=SIGS(IBM,JGR,IL)-CSCAT
                SS2D(IBM,JGR,IGR,IL)=SS2D(IBM,JGR,IGR,IL)-CSCAT*FF
                SS2D(IBM,IGR,JGR,IL)=0.0
              ENDDO
            ENDDO
          ENDDO
        ENDIF
*
        IF(LSPH) THEN
          CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"LOCALVALUE",RVALO)
          CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"LOCALVALADDR",LOCAD)
          DO ILOC=1,NLOC
            IF((LOCTYP(ILOC).EQ.'EQUI').AND.(LOCKEY(ILOC).EQ.HEQUI))
     &      THEN
              IF(LOCAD(ILOC+1)-LOCAD(ILOC).NE.NGROUP) THEN
                CALL XABORT('SPHMPO: INVALID NUMBER OF COMPONENTS FOR '
     &          //'SPH FACTORS')
              ENDIF
              DO IGR=1,NGROUP
                SPH(IBM,IGR)=RVALO(LOCAD(ILOC)+IGR-1)
              ENDDO
            ENDIF
          ENDDO
          DEALLOCATE(LOCAD,RVALO)
        ELSE
          SPH(IBM,:NGROUP)=1.0
        ENDIF
        IF(LMASL) THEN
          CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"LOCALVALUE",RVALO)
          CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"LOCALVALADDR",LOCAD)
          DO ILOC=1,NLOC
            IF((LOCTYP(ILOC).EQ.'HEAVY_METAL_DENSITY').AND.
     &      (LOCKEY(ILOC).EQ.HMASL)) THEN
              IF(LOCAD(ILOC+1)-LOCAD(ILOC).NE.1) THEN
                CALL XABORT('SPHMPO: INVALID NUMBER OF COMPONENTS FOR '
     &          //'HEAVY_METAL_DENSITY')
              ENDIF
              FMASL(IBM)=RVALO(LOCAD(ILOC))
            ENDIF
          ENDDO
          DEALLOCATE(RVALO)
        ENDIF
*----
*  RECOVER DIFFUSION COEFFICIENT INFORMATION
*----
        IF(hdf5_group_exists(IPMPO,TRIM(RECNAM)//"leakage")) THEN
          CALL hdf5_info(IPMPO,TRIM(RECNAM)//"leakage/DIFFCOEF",RANK,
     &    TYPE,NBYTE,DIMSR)
          IF(TYPE.NE.99) THEN
            CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"leakage/DIFFCOEF",
     &      VREAL)
            DO IGR=1,NGROUP
              DIFF(IBM,IGR)=VREAL(IGR)
            ENDDO
            DEALLOCATE(VREAL)
            LDIFF=.TRUE.
            LSTRD=.FALSE.
            GO TO 10
          ENDIF
          CALL hdf5_info(IPMPO,TRIM(RECNAM)//"leakage/DB2",RANK,TYPE,
     &    NBYTE,DIMSR)
          IF(TYPE.NE.99) THEN
            CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"leakage/DB2",VREAL)
            DO IGR=1,NGROUP
              DIFF(IBM,IGR)=VREAL(IGR)/B2
            ENDDO
            DEALLOCATE(VREAL)
            LDIFF=.TRUE.
            LSTRD=.FALSE.
          ENDIF
        ENDIF
   10   CONTINUE
      ENDDO ! end of loop over mixtures
      IF(NALBP.GT.0) THEN
        SPH(NMIL+1:NMIL+NALBP,:NGROUP)=1.0 ! assigned to albedo function
      ENDIF
      DEALLOCATE(ADDRISO,IDATAP,ISOTOPE,REACTION,ADDRXS)
*----
*  IDENTIFY SPECIAL FLUX EDITS
*----
      ALLOCATE(HEDI(NREA))
      NED=0
      DO IREA=1,NREA
        IF((NOMREA(IREA).EQ.'Total').or.
     &     (NOMREA(IREA).EQ.'Absorption').or.
     &     (NOMREA(IREA).EQ.'CaptureEnergyCapture').or.
     &     (NOMREA(IREA).EQ.'Diffusion').or.
     &     (NOMREA(IREA).EQ.'FissionEnergyFission').or.
     &     (NOMREA(IREA).EQ.'FissionSpectrum').or.
     &     (NOMREA(IREA).EQ.'NuFission').or.
     &     (NOMREA(IREA).EQ.'Scattering')) CYCLE
        NED=NED+1
        IF(NOMREA(IREA).EQ.'Fission') THEN
          HEDI(NED)='NFTOT'
        ELSE
          HEDI(NED)=NOMREA(IREA)(:8)
        ENDIF
      ENDDO
*----
*  STORE MACROLIB.
*----
      CALL LCMPUT(IPMAC,'VOLUME',NMIL,2,VOSAP)
      IF(LMASL) CALL LCMPUT(IPMAC,'MASL',NMIL,2,FMASL)
      IFISS=0
      ITRANC=0
      LHFACT=.FALSE.
      ALLOCATE(VREAL(NMIL))
      JPMAC=LCMLID(IPMAC,'GROUP',NGROUP)
      DO IGR=1,NGROUP
        KPMAC=LCMDIL(JPMAC,IGR)
        CALL LCMPUT(KPMAC,'FLUX-INTG',NMIL,2,NWT0(1,IGR))
        IF(NPRC.GT.0) THEN
          DO IBM=1,NMIL
            VREAL(IBM)=INVELS(IGR,IBM)
          ENDDO
          CALL LCMPUT(KPMAC,'OVERV',NMIL,2,VREAL)
        ENDIF
        DO IREA=1,NREA
          IF(.NOT.LXS(IREA)) CYCLE
          IF(NOMREA(IREA).EQ.'Absorption') THEN
            TOTAL(:,IGR)=TOTAL(:,IGR)+XS(:,IGR,IREA)
          ELSE IF(NOMREA(IREA).EQ.'Nexcess') THEN
*           correct scattering XS with excess XS
            SIGS0(:,IGR)=SIGS0(:,IGR)+XS(:,IGR,IREA)
            CALL LCMPUT(KPMAC,'N2N',NMIL,2,XS(1,IGR,IREA))
          ELSE IF(NOMREA(IREA).EQ.'Fission') THEN
            CALL LCMPUT(KPMAC,'NFTOT',NMIL,2,XS(1,IGR,IREA))
          ELSE IF(NOMREA(IREA).EQ.'FissionSpectrum') THEN
            CALL LCMPUT(KPMAC,'CHI',NMIL,2,XS(1,IGR,IREA))
            DO IPRC=1,NPRC
              DO IBM=1,NMIL
                VREAL(IBM)=CHIRS(IGR,IPRC,IBM)
              ENDDO
              WRITE(TEXT12,'(A3,I2.2)') 'CHI',IPRC
              CALL LCMPUT(KPMAC,TEXT12,NMIL,2,VREAL)
            ENDDO
          ELSE IF(NOMREA(IREA).EQ.'NuFission') THEN
            IFISS=1
            CALL LCMPUT(KPMAC,'NUSIGF',NMIL,2,XS(1,IGR,IREA))
            DO IPRC=1,NPRC
              DO IBM=1,NMIL
                VREAL(IBM)=XS(IBM,IGR,IREA)*BETAR(IPRC,IBM)
              ENDDO
              WRITE(TEXT12,'(A6,I2.2)') 'NUSIGF',IPRC
              CALL LCMPUT(KPMAC,TEXT12,NMIL,2,VREAL)
            ENDDO
          ELSE IF(NOMREA(IREA).EQ.'CaptureEnergyCapture') THEN
            LHFACT=.TRUE.
            EFACT(:,IGR)=EFACT(:,IGR)+XS(:,IGR,IREA)
          ELSE IF(NOMREA(IREA).EQ.'FissionEnergyFission') THEN
            LHFACT=.TRUE.
            EFACT(:,IGR)=EFACT(:,IGR)+XS(:,IGR,IREA)
          ELSE IF(NOMREA(IREA).EQ.'TransportCorrection') THEN
            ITRANC=2
            CALL LCMPUT(KPMAC,'TRANC',NMIL,2,XS(1,IGR,IREA))
          ELSE IF(NOMREA(IREA).EQ.'Diffusion') THEN
            DO IL=1,NL
              WRITE(CM,'(I2.2)') IL-1
              IF(IL.EQ.1) THEN
                DO IBM=1,NMIL
                  SIGS0(IBM,IGR)=SIGS0(IBM,IGR)+SIGS(IBM,IGR,IL)
                  TOTAL(IBM,IGR)=TOTAL(IBM,IGR)+SIGS(IBM,IGR,IL)
                ENDDO
              ELSE
                CALL LCMPUT(KPMAC,'SIGS'//CM,NMIL,2,SIGS(1,IGR,IL))
              ENDIF
            ENDDO
          ELSE IF(NOMREA(IREA).EQ.'Scattering') THEN
            ALLOCATE(SCAT(NGROUP*NMIL),GAR(NMIL))
            DO IL=1,NL
              WRITE(CM,'(I2.2)') IL-1
              IPOSDE=0
              DO IBM=1,NMIL
                IPOS(IBM)=IPOSDE+1
                IGMIN=IGR
                IGMAX=IGR
                DO JGR=NGROUP,1,-1
                  IF(SS2D(IBM,IGR,JGR,IL).NE.0.0) THEN
                    IGMIN=MIN(IGMIN,JGR)
                    IGMAX=MAX(IGMAX,JGR)
                  ENDIF
                ENDDO
                IJJM(IBM)=IGMAX
                NJJM(IBM)=IGMAX-IGMIN+1
                DO JGR=IGMAX,IGMIN,-1
                  IPOSDE=IPOSDE+1
                  SCAT(IPOSDE)=SS2D(IBM,IGR,JGR,IL) ! IGR <-- JGR
                ENDDO
                GAR(IBM)=SCAT(IPOS(IBM)+IJJM(IBM)-IGR)
              ENDDO
              CALL LCMPUT(KPMAC,'SCAT'//CM,IPOSDE,2,SCAT)
              CALL LCMPUT(KPMAC,'NJJS'//CM,NMIL,1,NJJM)
              CALL LCMPUT(KPMAC,'IJJS'//CM,NMIL,1,IJJM)
              CALL LCMPUT(KPMAC,'IPOS'//CM,NMIL,1,IPOS)
              CALL LCMPUT(KPMAC,'SIGW'//CM,NMIL,2,GAR)
            ENDDO
            DEALLOCATE(GAR,SCAT)
          ELSE
            CALL LCMPUT(KPMAC,NOMREA(IREA),NMIL,2,XS(1,IGR,IREA))
          ENDIF
        ENDDO ! end of loop over reactions
        IF(LSTRD) THEN
          IF((ITRANC.EQ.0).AND.(NL.GT.1)) THEN
*           Apollo-type transport correction
            VREAL(:)=TOTAL(:,IGR)-SIGS(:,IGR,2)
          ELSE
            VREAL(:)=TOTAL(:,IGR)
          ENDIF
          DO IBM=1,NMIL
            DIFF(IBM,IGR)=1.0/(3.0*VREAL(IBM))
          ENDDO
          LDIFF=.TRUE.
        ENDIF
        IF((ITRANC.EQ.0).AND.(NL.GT.1)) THEN
*         Apollo-type transport correction
          IF(IGR.EQ.NGROUP) ITRANC=2
          CALL LCMPUT(KPMAC,'TRANC',NMIL,2,SIGS(1,IGR,2))
        ENDIF
        IF(LDIFF) CALL LCMPUT(KPMAC,'DIFF',NMIL,2,DIFF(1,IGR))
        IF(LHFACT) CALL LCMPUT(KPMAC,'H-FACTOR',NMIL,2,EFACT(1,IGR))
      ENDDO
      DEALLOCATE(VREAL)
*----
*  RELEASE MEMORY
*----
      DEALLOCATE(INVELS,BETAR,CHIRS,LAMB)
      DEALLOCATE(ADR,FAG,XS,SS2D,SIGS,EFACT,NWT0,LXS)
      DEALLOCATE(NOMISO,NOMREA)
*----
*  SAVE SCATTERING P0 AND TOTAL CROSS SECTION INFO
*----
      DO IGR=1,NGROUP
        KPMAC=LCMDIL(JPMAC,IGR)
        CALL LCMPUT(KPMAC,'SIGS00',NMIL,2,SIGS0(1,IGR))
        CALL LCMPUT(KPMAC,'NTOT0',NMIL,2,TOTAL(1,IGR))
      ENDDO
*----
*  WRITE STATE VECTOR
*----
      IF(IMPX.GT.1) THEN
        WRITE(6,'(32H SPHMPO: fissile isotope index =,I4)') IFISS
        WRITE(6,'(37H SPHMPO: transport correction index =,I4)') ITRANC
      ENDIF
      TEXT12='L_MACROLIB'
      CALL LCMPTC(IPMAC,'SIGNATURE',12,1,TEXT12)
      ISTATE(:NSTATE)=0
      ISTATE(1)=NGROUP
      ISTATE(2)=NMIL
      ISTATE(3)=NL ! 1+scattering anisotropy
      ISTATE(4)=IFISS
      ISTATE(5)=NED
      ISTATE(6)=ITRANC
      ISTATE(8)=NALBP
      IF(LDIFF) ISTATE(9)=1
      ISTATE(12)=IDF ! ADF information
      CALL LCMPUT(IPMAC,'STATE-VECTOR',NSTATE,1,ISTATE)
      IF(NED.GT.0) CALL LCMPTC(IPMAC,'ADDXSNAME-P0',8,NED,HEDI)
      DEALLOCATE(HEDI)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FMASL,DIFF,TOTAL,SIGS0)
      DEALLOCATE(IJJM,NJJM,IPOS)
      RETURN
      END
