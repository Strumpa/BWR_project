*DECK APXCA2
      SUBROUTINE APXCA2(IPAPX,IPMICR,NREA,NISO,NMAC,NED,NPRC,NG,NL,
     1 ITRANC,NALBP,IMC,NMIL,NBISO,ICAL,IMPX,FNORM,NMILNR,NISFS,NISPS,
     2 VOLMIL,FLXMIL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the cross sections of an elementary calculation in the Apex
* file.
*
*Copyright:
* Copyright (C) 2025 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPAPX   pointer to the Apex file.
* IPMICR  pointer to the microlib to include (L_LIBRARY signature).
* NREA    number of requested reactions.
* NISO    number of particularized isotopes.
* NMAC    number of macros.
* NED     number of additional edition cross sections.
* NPRC    number of delayed neutron precursors.
* NG      number of condensed energy groups.
* NL      number of Legendre orders.
* ITRANC  type of transport correction.
* NALBP   number of physical albedos per energy group.
* IMC     type of macro-calculation (1 for diffusion or SPN;
*         2 other method).
* NMIL    number of mixtures in the Apex file.
* NBISO   number of isotopes in the condensed microlib of the edition
*         object. A given isotope may appear in many mixtures.
* ICAL    index of the current elementary calculation.
* FNORM   flux normalization factor.
* IMPX    print parameter.
*
*Parameters: output
* NMILNR  number of mixtures with delayed neutron data.
* NISFS   number of particularized fissile isotopes.
* NISPS   number of particularized fission products.
* VOLMIL  mixture volumes.
* FLXMIL  averaged flux of mixtures.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPAPX,IPMICR
      INTEGER NREA,NISO,NMAC,NED,NPRC,NG,NL,ITRANC,NALBP,IMC,NMIL,NBISO,
     1 ICAL,IMPX,NMILNR,NISFS,NISPS
      REAL FNORM,VOLMIL(NMIL),FLXMIL(NMIL,NG)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPEDIT,KPEDIT,IPTEMP
      CHARACTER RECNAM*80,RECNAM2*80,TEXT8*8,TEXT12*12,HSMG*131
      LOGICAL EXIST,LSPH
      DOUBLE PRECISION CONV
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MIX,ITYPE,IDEPL
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISONAM
      REAL, ALLOCATABLE, DIMENSION(:) :: OVERV,WORKD,WORK1,WORK2,DEN,
     1 DENISO,ENRGS,VOLMIX,WORK1D
      REAL, ALLOCATABLE, DIMENSION(:,:) :: DNUSIG,DCHI,SPH,CONCES,DECAYC
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:,:) :: IPERM
      CHARACTER(LEN=4), ALLOCATABLE, DIMENSION(:) :: TYPISO
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: NOMISO,NOMMAC
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: NOMREA
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ISONAM(3,NBISO),MIX(NBISO),ITYPE(NBISO),IDEPL(NBISO))
      ALLOCATE(OVERV(NG),DNUSIG(NG,NPRC+1),DCHI(NG,NPRC),WORKD(NPRC),
     1 WORK1(NG*NMIL+1),WORK2(NG),DEN(NBISO),DENISO(NISO),
     2 CONCES(NISO,NMIL),IPERM(NISO,NMIL),VOLMIX(NMIL))
*
      CONV=1.0D6 ! convert MeV to eV in H-FACTOR
*----
*  RECOVER INFORMATION FROM THE 'explicit' GROUP.
*----
      IF(NREA.GT.0) CALL hdf5_read_data(IPAPX,"/explicit/REANAME",
     1 NOMREA)
      IF(NMAC.GT.0) CALL hdf5_read_data(IPAPX,"/explicit/MACNAME",
     1 NOMMAC)
      IF(NISO.GT.0) THEN
        CALL hdf5_read_data(IPAPX,"/physconst/ISOTA",NOMISO)
        CALL hdf5_read_data(IPAPX,"/physconst/ISOTYP",TYPISO)
      ENDIF
*----
*  SAVE INFORMATION TO THE 'physconst' GROUP.
*----
      IF(ICAL.EQ.1) THEN
         ALLOCATE(ENRGS(NG+1))
         CALL LCMLEN(IPMICR,'ENERGY',ILONG,ITYLCM)
         IF(ILONG.EQ.0) THEN
            CALL LCMSIX(IPMICR,'MACROLIB',1)
            CALL LCMLEN(IPMICR,'ENERGY',ILONG,ITYLCM)
            IF(ILONG.NE.NG+1) CALL XABORT('APXCA2: BAD VALUE OF NG(1).')
            CALL LCMGET(IPMICR,'ENERGY',ENRGS)
            CALL LCMSIX(IPMICR,' ',2)
         ELSE
            IF(ILONG.NE.NG+1) CALL XABORT('APXCA2: BAD VALUE OF NG(2).')
            CALL LCMGET(IPMICR,'ENERGY',ENRGS)
         ENDIF
         ENRGS(:NG+1)=ENRGS(:NG+1)*1.0E-6
         CALL hdf5_write_data(IPAPX,"/physconst/ENRGS",ENRGS)
         DEALLOCATE(ENRGS)
         CALL LCMSIX(IPMICR,'MACROLIB',1)
         CALL LCMLEN(IPMICR,'VOLUME',ILONG,ITYLCM)
         IF(ILONG.NE.NMIL) CALL XABORT('APXCA2: INCORRECT VOLUME.')
         CALL LCMGET(IPMICR,'VOLUME',VOLMIL)
         CALL LCMSIX(IPMICR,' ',2)
      ENDIF
*----
*  RECOVER INVERSE OF SPH EQUIVALENCE FACTORS.
*----
      CALL LCMSIX(IPMICR,'MACROLIB',1)
      CALL LCMGET(IPMICR,'VOLUME',VOLMIX)
      JPEDIT=LCMGID(IPMICR,'GROUP')
      LSPH=.FALSE.
      ALLOCATE(SPH(NMIL+NALBP,NG))
      DO 80 IGR=1,NG
      KPEDIT=LCMGIL(JPEDIT,IGR)
      CALL LCMLEN(KPEDIT,'NSPH',ILONG,ITYLCM)
      IF(ILONG.GT.0) THEN
         LSPH=.TRUE.
         CALL LCMGET(KPEDIT,'NSPH',WORK1)
         DO 70 IMIL=1,NMIL
         SPH(IMIL,IGR)=1.0/WORK1(IMIL)
   70    CONTINUE
         DO 75 IALB=1,NALBP
         SPH(NMIL+IALB,IGR)=1.0
   75    CONTINUE
      ELSE
         SPH(:NMIL+NALBP,IGR)=1.0
      ENDIF
   80 CONTINUE
      CALL LCMSIX(IPMICR,' ',2)
*----
*  CREATE A SPH-UNCORRECTED MICROLIB.
*----
      CALL LCMOP(IPTEMP,'*TEMPORARY*',0,1,0)
      CALL LCMEQU(IPMICR,IPTEMP)
      IF(LSPH) THEN
        IF(IMC.EQ.0) CALL XABORT('APXCA2: UNDEFINED TYPE OF SPH.')
        NW=1 ! NTOT1 cross section present
        CALL SPHCMI(IPTEMP,0,IMC,NMIL,NBISO,NG,NL,NW,NED,NPRC,NALBP,SPH)
      ENDIF
      DEALLOCATE(SPH)
*----
*  FIND THE NUMBER AND NAMES OF THE ISOTOPES IN THE OUTPUT TABLES.
*----
      IF(NISO.GT.0) THEN
        IPERM(:NISO,:NMIL)=C_NULL_PTR
        CONCES(:NISO,:NMIL)=0.0
        IF(NBISO.GT.0) THEN
          ALLOCATE(IPISO(NBISO))
          CALL LCMGET(IPTEMP,'ISOTOPESUSED',ISONAM)
          CALL LCMGET(IPTEMP,'ISOTOPESMIX',MIX)
          CALL LCMGET(IPTEMP,'ISOTOPESDENS',DEN)
          CALL LCMGET(IPMICR,'ISOTOPESTYPE',ITYPE)
          CALL LCMGET(IPMICR,'ISOTOPESTODO',IDEPL)
          CALL LIBIPS(IPTEMP,NBISO,IPISO)
          DO IBISO=1,NBISO
            IF(IDEPL(IBISO).EQ.1) CYCLE
            IMIL=MIX(IBISO)
            IF(IMIL.EQ.0) CYCLE
            WRITE(TEXT12,'(3A4)') (ISONAM(I0,IBISO),I0=1,3)
            DO ISO=1,NISO
              IF(NOMISO(ISO).EQ.TEXT12(:8)) THEN
                IPERM(ISO,IMIL)=IPISO(IBISO)
                CONCES(ISO,IMIL)=DEN(IBISO)
                CYCLE
              ENDIF
            ENDDO
          ENDDO
          DEALLOCATE(IPISO)
        ENDIF
        DO ISO=1,NISO
          DO IMIL=1,NMIL
            IF(CONCES(ISO,IMIL).EQ.0.0) GO TO 10
            IF(C_ASSOCIATED(IPERM(ISO,IMIL))) GO TO 10
          ENDDO
          WRITE(HSMG,'(17HAPXCA2: ISOTOPE '',A8,7H'' (ISO=,I8,3H) I,
     1    32HS NOT AVAILABLE IN THE MICROLIB.)') NOMISO(ISO),ISO
          CALL XABORT(HSMG)
   10     CONTINUE
        ENDDO
*----
*  RECOVER RADIOACTIVE DECAY CONSTANTS.
*----
        IF(ICAL.EQ.1) THEN
          ALLOCATE(DECAYC(1,NISO),IPISO(NBISO))
          CALL LIBIPS(IPTEMP,NBISO,IPISO)
          DECAYC(1,:NISO)=0.0
          DO 40 ISO=1,NISO
          IISOTS=0
          DO 20 IBISO=1,NBISO
          IF(IDEPL(IBISO).EQ.1) GO TO 20
          IF(MIX(IBISO).EQ.0) GO TO 20
          IISOTS=ISO
          WRITE(TEXT12,'(3A4)') (ISONAM(I0,IBISO),I0=1,3)
          IF(TEXT12(:8).EQ.NOMISO(ISO)) GO TO 30
   20     CONTINUE
          WRITE(6,'(48H APXCA2: **WARNING** UNABLE TO FIND THE RADIOACT,
     1    31HIVE DECAY CONSTANT FOR ISOTOPE ,A,1H.)') TRIM(NOMISO(ISO))
          GO TO 40
   30     JPEDIT=IPISO(IISOTS)
          IF(.NOT.C_ASSOCIATED(JPEDIT)) GO TO 40
          CALL LCMLEN(JPEDIT,'DECAY',ILONG,ITYLCM)
          IF(ILONG.EQ.1) CALL LCMGET(JPEDIT,'DECAY',DECAYC(1,ISO))
   40     CONTINUE
          DECAYC(1,:NISO)=DECAYC(1,:NISO)*1.0E-8
          CALL hdf5_write_data(IPAPX,"/physconst/DECAYC",DECAYC)
          DEALLOCATE(IPISO,DECAYC)
        ENDIF
      ENDIF
*----
*  FILL miscellaneous GROUP
*----
      WRITE(RECNAM,'(4Hcalc,I8,15H/miscellaneous/)') ICAL
      CALL LCMSIX(IPTEMP,'MACROLIB',1)
      NVDIV=0
      CALL LCMLEN(IPTEMP,'K-EFFECTIVE',ILONG,ITYLCM)
      IF(ILONG.EQ.1) THEN
         CALL LCMGET(IPTEMP,'K-EFFECTIVE',FLOTT)
         CALL hdf5_write_data(IPAPX,TRIM(RECNAM)//"KEFF",FLOTT)
      ENDIF
      CALL LCMLEN(IPTEMP,'K-INFINITY',ILONG,ITYLCM)
      IF(ILONG.EQ.1) THEN
         CALL LCMGET(IPTEMP,'K-INFINITY',FLOTT)
         CALL hdf5_write_data(IPAPX,TRIM(RECNAM)//"KINF",FLOTT)
      ENDIF
      CALL LCMLEN(IPTEMP,'B2  B1HOM',ILONG,ITYLCM)
      IF(ILONG.EQ.1) THEN
         CALL LCMGET(IPTEMP,'B2  B1HOM',B2)
      ELSE
         B2=0.0
      ENDIF
      IF(B2.EQ.0.0) B2=1.0E-10
      CALL hdf5_write_data(IPAPX,TRIM(RECNAM)//"B2",B2)
      CALL LCMSIX(IPTEMP,' ',2)
*----
*  LOOP OVER APEX MIXTURES.
*----
      NMILNR=0
      DO 500 IMIL=1,NMIL
      IF(NMIL.EQ.1) THEN
        WRITE(RECNAM,'(4Hcalc,I8,4H/xs/)') ICAL
      ELSE
        WRITE(RECNAM,'(4Hcalc,I8,3H/xs,I8,1H/)') ICAL,IMIL
      ENDIF
      CALL hdf5_create_group(IPAPX,RECNAM)
*----
*  RECOVER APEX VOLUMES AND INTEGRATED FLUXES.
*----
      CALL hdf5_write_data(IPAPX,TRIM(RECNAM)//"MEDIA_VOLUME",
     1 VOLMIX(IMIL))
      WORK2(:NG)=0.0
      CALL LCMSIX(IPTEMP,'MACROLIB',1)
      JPEDIT=LCMGID(IPTEMP,'GROUP')
      DO IGR=1,NG
        KPEDIT=LCMGIL(JPEDIT,IGR)
        CALL LCMLEN(KPEDIT,'FLUX-INTG',ILONG,ITYLCM)
        IF(ILONG.EQ.0) CYCLE
        ALLOCATE(WORK1D(NMIL))
        CALL LCMGET(KPEDIT,'FLUX-INTG',WORK1D)
        WORK2(IGR)=WORK1D(IMIL)
        DEALLOCATE(WORK1D)
      ENDDO
      CALL hdf5_write_data(IPAPX,TRIM(RECNAM)//"FLUX",WORK2)
      CALL LCMSIX(IPTEMP,' ',2)
*----
*  RECOVER APEX CROSS SECTIONS
*----
      IF(NISO.GT.0) THEN
        CALL hdf5_create_group(IPAPX,TRIM(RECNAM)//"mic")
        RECNAM2=TRIM(RECNAM)//"mic/CONC"
        CALL hdf5_write_data(IPAPX,TRIM(RECNAM2),CONCES(:NISO,IMIL))
      ENDIF
      CALL hdf5_create_group(IPAPX,TRIM(RECNAM)//"mac")
      DO IREA=1,NREA
        CALL APXSX2(IPAPX,IPTEMP,NG,NL,NMAC,NISO,NMIL,IMIL,ITRANC,
     1  RECNAM,NOMMAC,TYPISO,NOMREA(IREA),IPERM(1,IMIL),CONCES(1,IMIL),
     2  B2)
      ENDDO
      IF(IMPX.GT.0) THEN
        CALL hdf5_list(IPAPX,TRIM(RECNAM))
        IF(NISO.GT.0) CALL hdf5_list(IPAPX,TRIM(RECNAM)//"mic")
        CALL hdf5_list(IPAPX,TRIM(RECNAM)//"mac")
      ENDIF
      IOR=0
      IOI=0
      IIS=0
      NISMAX=NMAC
*
      CALL LCMSIX(IPTEMP,'MACROLIB',1)
      JPEDIT=LCMGID(IPTEMP,'GROUP')
      DO 150 IGR=1,NG
      KPEDIT=LCMGIL(JPEDIT,IGR)
*----
*  RECOVER THE NEUTRON FLUX.
*----
      CALL LCMGET(KPEDIT,'FLUX-INTG',WORK1)
      IF(FNORM.NE.1.0) THEN
        FLXMIL(IMIL,IGR)=WORK1(IMIL)*FNORM*1.0E13
      ELSE
        FLXMIL(IMIL,IGR)=WORK1(IMIL)
      ENDIF
*----
*  RECOVER DELAYED NEUTRON INFORMATION.
*----
      CALL LCMLEN(KPEDIT,'NUSIGF',ILONG,ITYLCM)
      IF((NPRC.GT.0).AND.(ILONG.NE.0)) THEN
         CALL LCMGET(KPEDIT,'NUSIGF',WORK1)
         DNUSIG(IGR,NPRC+1)=WORK1(IMIL)
         CALL LCMGET(KPEDIT,'OVERV',WORK1)
         OVERV(IGR)=WORK1(IMIL)
         DO 90 IPRC=1,NPRC
         WRITE(TEXT12,'(6HNUSIGF,I2.2)') IPRC
         CALL LCMGET(KPEDIT,TEXT12,WORK1)
         DNUSIG(IGR,IPRC)=WORK1(IMIL)
         WRITE(TEXT12,'(3HCHI,I2.2)') IPRC
         CALL LCMGET(KPEDIT,TEXT12,WORK1)
         DCHI(IGR,IPRC)=WORK1(IMIL)
   90    CONTINUE
      ELSE
         DNUSIG(IGR,:NPRC+1)=0.0
      ENDIF
  150 CONTINUE
      CALL LCMSIX(IPTEMP,' ',2)
*----
*  STORE INFORMATION IN THE calc_id/kinetics GROUP.
*----
      IF(NPRC.GT.0) THEN
         EXIST=.FALSE.
         DO 455 IPRC=1,NPRC
         DO 450 IGR=1,NG
         EXIST=EXIST.OR.(DNUSIG(IGR,IPRC).NE.0.0)
  450    CONTINUE
  455    CONTINUE
         IF(EXIST) THEN
            NMILNR=NMILNR+1
            RECNAM2=TRIM(RECNAM)//"kinetics/"
            CALL hdf5_create_group(IPAPX,TRIM(RECNAM2))
            CALL hdf5_write_data(IPAPX,TRIM(RECNAM2)//"NBGRD",NPRC)
            CALL hdf5_write_data(IPAPX,TRIM(RECNAM2)//"CHIDA",DCHI)
            CALL hdf5_write_data(IPAPX,TRIM(RECNAM2)//"INVELA",OVERV)
            CALL LCMSIX(IPTEMP,'MACROLIB',1)
            CALL LCMLEN(IPTEMP,'LAMBDA-D',ILONG,ITYLCM)
            IF(ILONG.GT.NPRC) CALL XABORT('APXCA2: LAMBDA-D OVERFLOW.')
            CALL LCMGET(IPTEMP,'LAMBDA-D',WORKD)
            CALL LCMSIX(IPTEMP,' ',2)
            CALL hdf5_write_data(IPAPX,TRIM(RECNAM2)//"LAMBDA",WORKD)
            TGENRS=0.0
            DENOM=0.0
            DO 460 IGR=1,NG
            TGENRS=TGENRS+OVERV(IGR)*FLXMIL(IMIL,IGR)
            DENOM=DENOM+DNUSIG(IGR,NPRC+1)*FLXMIL(IMIL,IGR)
  460       CONTINUE
            TGENRS=TGENRS/DENOM
            DO 480 IPRC=1,NPRC
            WORKD(IPRC)=0.0
            DO 470 IGR=1,NG
            WORKD(IPRC)=WORKD(IPRC)+DNUSIG(IGR,IPRC)*FLXMIL(IMIL,IGR)
  470       CONTINUE
            WORKD(IPRC)=WORKD(IPRC)/DENOM
  480       CONTINUE
            CALL hdf5_write_data(IPAPX,TRIM(RECNAM2)//"BETADA",WORKD)
            CALL hdf5_write_data(IPAPX,TRIM(RECNAM2)//"NGENT",TGENRS)
            IF(IMPX.GT.0) CALL hdf5_list(IPAPX,TRIM(RECNAM2))
         ENDIF
      ENDIF
  500 CONTINUE
*----
*  COMPUTE NISFS AND NISPS
*----
      NISFS=0
      NISPS=0
      DO 530 ISO=1,NISO
      DO 510 IBISO=1,NBISO
      WRITE(TEXT8,'(2A4)') (ISONAM(I0,IBISO),I0=1,2)
      IF(NOMISO(ISO).EQ.TEXT8) THEN
         ITY=ITYPE(IBISO)
         GO TO 520
      ENDIF
  510 CONTINUE
      GO TO 530
  520 IF(ITY.EQ.2) THEN
         NISFS=NISFS+1
      ELSE IF(ITY.EQ.3) THEN
         NISPS=NISPS+1
      ENDIF
  530 CONTINUE
      CALL LCMCL(IPTEMP,2)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      IF(NISO.GT.0) DEALLOCATE(TYPISO,NOMISO)
      IF(NMAC.GT.0) DEALLOCATE(NOMMAC)
      IF(NREA.GT.0) DEALLOCATE(NOMREA)
      DEALLOCATE(VOLMIX,IPERM,CONCES,DENISO,DEN,WORK2,WORK1,WORKD,DCHI,
     1 DNUSIG,OVERV,IDEPL,ITYPE,MIX,ISONAM)
      RETURN
      END
