*DECK MCRAGF
      SUBROUTINE MCRAGF(IPMAC,IPMPO,IACCS,NMIL,NMIX,NGRP,NALBP,LALBG,
     1 IMPX,NCAL,TERP,MIXC,NSURFD,HEDIT,VOLMI2,IDF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Build the macrolib by scanning the NCAL elementary calculations and
* weighting them with TERP factors. ADF and physical albedos part.
*
*Copyright:
* Copyright (C) 2022 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPMAC   address of the output macrolib LCM object.
* IPMPO   address of the MPO file.
* IACCS   =0 macrolib is created; =1 ... is updated.
* NMIL    number of material mixtures in the MPO file.
* NMIX    maximum number of material mixtures in the macrolib.
* NGRP    number of energy groups.
* NALBP   number of physical albedos per energy group.
* LALBG   type of physical albedos (.true.: diagonal; .false.: GxG).
* IMPX    print parameter (equal to zero for no print).
* NCAL    number of elementary calculations in the MPO file.
* TERP    interpolation factors.
* MIXC    mixture index in the MPO file corresponding to each macrolib
*         mixture. Equal to zero if a macrolib mixture is not updated.
* NSURFD  number of ADF.
* HEDIT   name of output group for a (multigroup mesh, output geometry)
*         couple (generally equal to 'output_0').
* VOLMI2  mixture volumes in the macrolib.
*
*Parameters: output
* IDF     type of discontinuity factors (DF) in the macrolib (=0: not
*         used; =3: DF used; =4 matrix DF used). 
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC,IPMPO
      INTEGER IACCS,NMIL,NMIX,NGRP,NALBP,IMPX,NCAL,MIXC(NMIX),NSURFD,IDF
      REAL TERP(NCAL,NMIX),VOLMI2(NMIX)
      LOGICAL LALBG
      CHARACTER(LEN=12) HEDIT
*----
*  LOCAL VARIABLES
*----
      REAL WEIGHT,FACTOR
      CHARACTER RECNAM*80
      INTEGER IKEFF,IKINF,I,IBM,IBMOLD,ICAL,IGR,JGR,ILONG,ITYLCM,IAL,
     1 RANK,NBYTE,TYPE,TYPE2,TYPE4,DIMSR(5),IOF,NSURFD_OLD
      DOUBLE PRECISION DGAR1,DGAR2
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR4,ZKINF,ZKEFF,VREAL,SURF,LG
      REAL, ALLOCATABLE, DIMENSION(:,:) :: GAR6,ALBP,AVGFL2,SURFLX
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ADF2,ALBP2,ALBP_ERM
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: ADF2M,ALBP2_E
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HADF
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(GAR4(NGRP*NGRP),GAR6(NGRP,2),ALBP2(NMIX,NALBP,NGRP),
     1 ALBP2_E(NMIX,NALBP,NGRP,NGRP),ZKINF(NMIX),ZKEFF(NMIX),
     2 HADF(NSURFD),ADF2(NMIX,NGRP,NSURFD),ADF2M(NMIX,NGRP,NGRP,NSURFD),
     3 AVGFL2(NMIX,NGRP))
*----
*  MICROLIB INITIALIZATION
*----
      IKINF=0
      IKEFF=0
      IDF=0
      IF(NSURFD.GT.0) THEN
        WRITE(RECNAM,'(8H/output/,A,16H/statept_0/flux/)') TRIM(HEDIT)
        CALL hdf5_info(IPMPO,TRIM(RECNAM)//"SURFFLUX",RANK,TYPE2,NBYTE,
     &  DIMSR)
        CALL hdf5_info(IPMPO,TRIM(RECNAM)//"SURFFLUXGxG",RANK,TYPE4,
     &  NBYTE,DIMSR)
        IF(TYPE2.NE.99) THEN
          IDF=2 ! boundary flux information
        ELSE IF(TYPE4.NE.99) THEN
          IDF=4 ! matrix discontinuity factor information
        ELSE
          CALL hdf5_list(IPMPO,TRIM(RECNAM))
          CALL XABORT('MCRAGF: UNABLE TO SET TYPE OF DF.')
        ENDIF
        ADF2(:NMIX,:NGRP,:NSURFD)=0.0
        ADF2M(:NMIX,:NGRP,:NGRP,:NSURFD)=0.0
      ENDIF
      AVGFL2(:NMIX,:NGRP)=0.0
      IF(NALBP.NE.0) ALBP2(:NMIX,:NALBP,:NGRP)=0.0
      ZKINF(:NMIX)=0.0
      ZKEFF(:NMIX)=0.0
      NSURFD_OLD=NSURFD
      IF(IACCS.NE.0) THEN
        ! Recover ADF, GFF and physical albedos previously computed
        CALL LCMLEN(IPMAC,'VOLUME',ILONG,ITYLCM)
        IF(ILONG.EQ.0) CALL XABORT('MCRAGF: NO VOLUMES IN MACROLIB.')
        CALL LCMGET(IPMAC,'VOLUME',VOLMI2)
        IF((NALBP.GT.0).AND.LALBG) THEN
*         diagonal albedo matrix
          ALLOCATE(ALBP(NALBP,NGRP))
          CALL LCMLEN(IPMAC,'ALBEDO',ILONG,ITYLCM)
          IF(ILONG.EQ.NALBP*NGRP) THEN
            CALL LCMGET(IPMAC,'ALBEDO',ALBP)
            DO IBM=1,NMIX ! mixtures in Macrolib
              ALBP2(IBM,:NALBP,:NGRP)=ALBP(:NALBP,:NGRP)
            ENDDO
          ENDIF
          DEALLOCATE(ALBP)
        ELSE IF(NALBP.GT.0) THEN
*         GxG albedo matrix
          ALLOCATE(ALBP_ERM(NALBP,NGRP,NGRP))
          CALL LCMLEN(IPMAC,'ALBEDO',ILONG,ITYLCM)
          IF(ILONG.EQ.NALBP*NGRP*NGRP) THEN
            CALL LCMGET(IPMAC,'ALBEDO',ALBP_ERM)
            DO IBM=1,NMIX ! mixtures in Macrolib
              ALBP2_E(IBM,:NALBP,:NGRP,:NGRP)=
     &                  ALBP_ERM(:NALBP,:NGRP,:NGRP)
            ENDDO
            DEALLOCATE(ALBP_ERM)
          ENDIF
        ENDIF
        IF(NSURFD_OLD.GT.0) THEN
          CALL LCMLEN(IPMAC,'ADF',ILONG,ITYLCM)
          IF(ILONG.EQ.0) THEN
            CALL LCMLIB(IPMAC)
            CALL XABORT('MCRAGF: UNABLE TO FIND DIRECTORY ADF.')
          ENDIF
          CALL LCMSIX(IPMAC,'ADF',1)
          CALL LCMGET(IPMAC,'NTYPE',NSURFD)
          IF(NSURFD.GT.NSURFD_OLD) THEN
            DEALLOCATE(ADF2M,ADF2,HADF)
            ALLOCATE(HADF(NSURFD),ADF2(NMIX,NGRP,NSURFD),
     &      ADF2M(NMIX,NGRP,NGRP,NSURFD))
            ADF2(:NMIX,:NGRP,:NSURFD)=0.0
            ADF2M(:NMIX,:NGRP,:NGRP,:NSURFD)=0.0
          ENDIF
          CALL LCMGTC(IPMAC,'HADF',8,NSURFD,HADF)
          DO I=1,NSURFD
            IF(IDF.EQ.2) THEN
              CALL LCMGET(IPMAC,HADF(I),ADF2(1,1,I))
            ELSE IF(IDF.EQ.4) THEN
              CALL LCMGET(IPMAC,HADF(I),ADF2M(1,1,1,I))
            ENDIF
          ENDDO
          CALL LCMGET(IPMAC,'AVG_FLUX',AVGFL2)
          CALL LCMSIX(IPMAC,' ',2)
        ENDIF
        DO IBM=1,NMIX
          IF(MIXC(IBM).NE.0) THEN
            IF(NALBP.NE.0) THEN
              IF(LALBG) THEN
                ALBP2(IBM,:NALBP,:NGRP)=0.0
              ELSE
                ALBP2_E(IBM,:NALBP,:NGRP,:NGRP)=0.0
              ENDIF
            ENDIF
            IF((NSURFD.GT.0).AND.(IDF.EQ.2)) THEN
              ADF2(IBM,:NGRP,:NSURFD)=0.0
            ELSE IF((NSURFD.GT.0).AND.(IDF.EQ.4)) THEN
              ADF2M(IBM,:NGRP,:NGRP,:NSURFD)=0.0
            ENDIF
            AVGFL2(IBM,:NGRP)=0.0
            ZKINF(IBM)=0.0
            ZKEFF(IBM)=0.0
          ENDIF
        ENDDO
      ENDIF
*----
*  OVERALL ELEMENTARY CALCULATION LOOP
*----
      DO 40 ICAL=1,NCAL
      WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,6H/flux/)')
     & TRIM(HEDIT),ICAL-1
      IF(NSURFD_OLD.GT.0) THEN
        CALL hdf5_info(IPMPO,TRIM(RECNAM)//"SURF",RANK,TYPE,NBYTE,DIMSR)
        IF(TYPE.NE.99) THEN
          CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"SURF",SURF)
        ELSE
*         temporary.....
          CALL hdf5_read_data(IPMPO,"/geometry/geometry_0/COORDINATE",
     &    LG)
          ALLOCATE(SURF(NSURFD_OLD))
          SURF(:NSURFD_OLD)=LG(2)
          DEALLOCATE(LG)
        ENDIF
        IF(IDF.EQ.2) THEN
          CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"SURFFLUX",SURFLX)
        ELSE IF(IDF.EQ.4) THEN
          CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"SURFFLUXGxG",SURFLX)
        ENDIF
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"TOTALFLUX",VREAL)
        DO 10 IBM=1,NMIX ! mixtures in Macrolib
        IBMOLD=MIXC(IBM)
        IF(IBMOLD.EQ.0) GO TO 10
        WEIGHT=TERP(ICAL,IBM)
        IF(WEIGHT.EQ.0.0) GO TO 10
        IF(IDF.EQ.2) THEN
          DO I=1,NSURFD_OLD
            WRITE(HADF(I),'(3HFD_,I5.5)') I
            DO IGR=1,NGRP
              IOF=(IGR-1)*NMIL+IBMOLD
              ADF2(IBM,IGR,I)=ADF2(IBM,IGR,I)+WEIGHT*SURFLX(I,IOF)/
     &        SURF(I)
            ENDDO
          ENDDO
        ELSE IF(IDF.EQ.4) THEN
          DO I=1,NSURFD_OLD
            WRITE(HADF(I),'(3HFD_,I5.5)') I
            DO JGR=1,NGRP
              DO IGR=1,NGRP
                IOF=((JGR-1)*NGRP+IGR-1)*NMIL+IBMOLD
                ADF2M(IBM,IGR,JGR,I)=ADF2M(IBM,IGR,JGR,I)+WEIGHT*
     &          SURFLX(I,IOF)/SURF(I)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        DO IGR=1,NGRP
          IOF=(IGR-1)*NMIL+IBMOLD
          AVGFL2(IBM,IGR)=AVGFL2(IBM,IGR)+WEIGHT*VREAL(IOF)/VOLMI2(IBM)
        ENDDO
   10   CONTINUE
        DEALLOCATE(SURF,VREAL,SURFLX)
      ENDIF
*----
*  PROCESS PHYSICAL ALBEDO INFORMATION
*----
      IF((NALBP.GT.0).AND.LALBG) THEN
*       diagonal albedo matrix
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ALBEDOG",ALBP)
        DO 20 IBM=1,NMIX ! mixtures in Macrolib
        IBMOLD=MIXC(IBM)
        IF(IBMOLD.EQ.0) GO TO 20
        WEIGHT=TERP(ICAL,IBM)
        IF(WEIGHT.EQ.0.0) GO TO 20
        DO IGR=1,NGRP
          DO IAL=1,NALBP
            FACTOR=(1.0-ALBP(IAL,IGR))/(1.0+ALBP(IAL,IGR))
            ALBP2(IBM,IAL,IGR)=ALBP2(IBM,IAL,IGR)+WEIGHT*FACTOR
          ENDDO
        ENDDO
   20     CONTINUE
        DEALLOCATE(ALBP)        
      ELSE IF(NALBP.GT.0) THEN
*       GxG albedo matrix
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ALBEDOGxG",ALBP_ERM)
        DO 25 IBM=1,NMIX ! mixtures in Macrolib
        IBMOLD=MIXC(IBM)
        IF(IBMOLD.EQ.0) GO TO 25
        WEIGHT=TERP(ICAL,IBM)
        IF(WEIGHT.EQ.0.0) GO TO 25
        DO IGR=1,NGRP
          DO JGR=1,NGRP
            DO IAL=1,NALBP
              FACTOR=(1.0-ALBP_ERM(IAL,IGR,JGR))/(1.0+
     1        ALBP_ERM(IAL,IGR,JGR))
              ALBP2_E(IBM,IAL,IGR,JGR)=ALBP2_E(IBM,IAL,IGR,JGR)+WEIGHT*
     1        FACTOR
            ENDDO
          ENDDO
        ENDDO
   25   CONTINUE
        DEALLOCATE(ALBP_ERM)        
      ENDIF
*----
*  PROCESS KINF AND KEFF
*----
      DO 30 IBM=1,NMIX ! mixtures in Macrolib
      IBMOLD=MIXC(IBM)
      IF(IBMOLD.EQ.0) GO TO 30
      WEIGHT=TERP(ICAL,IBM)
      IF(WEIGHT.EQ.0.0) GO TO 30
      WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,8H/addons/)')
     & TRIM(HEDIT),ICAL-1
      CALL hdf5_info(IPMPO,TRIM(RECNAM)//"KINF",RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) THEN
        IKINF=1
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"KINF",FACTOR)
        ZKINF(IBM)=ZKINF(IBM)+WEIGHT*FACTOR
      ENDIF
      CALL hdf5_info(IPMPO,TRIM(RECNAM)//"KEFF",RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) THEN
        IKEFF=1
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"KEFF",FACTOR)
        ZKEFF(IBM)=ZKEFF(IBM)+WEIGHT*FACTOR
      ENDIF
   30 CONTINUE
   40 CONTINUE
*----
*  SAVE ADF INFORMATION
*----
      IF((NSURFD.GT.0).AND.(IDF.EQ.2)) THEN
        CALL LCMSIX(IPMAC,'ADF',1)
        CALL LCMPUT(IPMAC,'NTYPE',1,1,NSURFD)
        CALL LCMPTC(IPMAC,'HADF',8,NSURFD,HADF)
        DO IBM=1,NMIX
          IF(MIXC(IBM).EQ.0) CYCLE
          IF(NSURFD.GT.NSURFD_OLD) THEN
            IF(NSURFD_OLD.NE.1) CALL XABORT('MCRAGF: INVALID NSURFD.')
*           assign the same ADF to all sides.
            DO I=2,NSURFD
              ADF2(IBM,:NGRP,I)=ADF2(IBM,:NGRP,1)
            ENDDO
          ENDIF
        ENDDO
        DO I=1,NSURFD
          CALL LCMPUT(IPMAC,HADF(I),NMIX*NGRP,2,ADF2(1,1,I))
        ENDDO
        CALL LCMPUT(IPMAC,'AVG_FLUX',NMIX*NGRP,2,AVGFL2)
        CALL LCMSIX(IPMAC,' ',2)
        IF(IMPX.GT.1) THEN
          DO IBM=1,NMIX
            IF(MIXC(IBM).EQ.0) CYCLE
            WRITE(6,'(/40H MCRAGF: DISCONTINUITY FACTORS - MIXTURE,I5)')
     1      IBM
            DO I=1,NSURFD
              WRITE(6,'(1X,A,1H:,1P,(5X,10E12.4))') TRIM(HADF(I)),
     1        (ADF2(IBM,IGR,I)/AVGFL2(IBM,IGR),IGR=1,NGRP)
            ENDDO
          ENDDO
        ENDIF
      ELSE IF((NSURFD.GT.0).AND.(IDF.EQ.4)) THEN
        CALL LCMSIX(IPMAC,'ADF',1)
        CALL LCMPUT(IPMAC,'NTYPE',1,1,NSURFD)
        CALL LCMPTC(IPMAC,'HADF',8,NSURFD,HADF)
        DO IBM=1,NMIX
          IF(MIXC(IBM).EQ.0) CYCLE
          DO JGR=1,NGRP
            DO IGR=1,NGRP
              ADF2M(IBM,IGR,JGR,:NSURFD)=ADF2M(IBM,IGR,JGR,:NSURFD)/
     1        AVGFL2(IBM,IGR)
            ENDDO
          ENDDO
          IF(NSURFD.GT.NSURFD_OLD) THEN
            IF(NSURFD_OLD.NE.1) CALL XABORT('MCRAGF: INVALID NSURFD.')
*           assign the same matrix ADF to all sides.
            DO I=2,NSURFD
              ADF2M(IBM,:NGRP,:NGRP,I)=ADF2M(IBM,:NGRP,:NGRP,1)
            ENDDO
          ENDIF
        ENDDO
        DO I=1,NSURFD
          CALL LCMPUT(IPMAC,HADF(I),NMIX*NGRP*NGRP,2,ADF2M(1,1,1,I))
        ENDDO
        CALL LCMPUT(IPMAC,'AVG_FLUX',NMIX*NGRP,2,AVGFL2)
        CALL LCMSIX(IPMAC,' ',2)
        IF(IMPX.GT.1) THEN
          DO IBM=1,NMIX
            IF(MIXC(IBM).EQ.0) CYCLE
            WRITE(6,'(/40H MCRAGF: DISCONTINUITY FACTORS - MIXTURE,I5)')
     1      IBM
            DO I=1,NSURFD
              WRITE(6,'(1X,A,1H:,1P,(5X,10E12.4))') TRIM(HADF(I)),
     1        ((ADF2M(IBM,IGR,JGR,I),IGR=1,NGRP),JGR=1,NGRP)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
*----
*  AVERAGE PHYSICAL ALBEDO INFORMATION
*----
      IF((NALBP.GT.0).AND.LALBG) THEN
*       diagonal albedo matrix
        ALLOCATE(ALBP(NALBP,NGRP))
        DO IGR=1,NGRP
          DO IAL=1,NALBP
            DGAR1=0.0D0
            DGAR2=0.0D0
            DO IBM=1,NMIX
              DGAR1=DGAR1+ALBP2(IBM,IAL,IGR)*VOLMI2(IBM)
              DGAR2=DGAR2+VOLMI2(IBM)
            ENDDO
            ALBP(IAL,IGR)=REAL((1.D0-DGAR1/DGAR2)/(1.D0+DGAR1/DGAR2))
          ENDDO
        ENDDO
        CALL LCMPUT(IPMAC,'ALBEDO',NALBP*NGRP,2,ALBP)
        DEALLOCATE(ALBP)
      ELSE IF(NALBP.GT.0) THEN
*       GxG albedo matrix
        ALLOCATE(ALBP_ERM(NALBP,NGRP,NGRP))
        DO IGR=1,NGRP
          DO JGR=1,NGRP
            DO IAL=1,NALBP
              DGAR1=0.0D0
              DGAR2=0.0D0
              DO IBM=1,NMIX
                DGAR1=DGAR1+ALBP2_E(IBM,IAL,IGR,JGR)*VOLMI2(IBM)
                DGAR2=DGAR2+VOLMI2(IBM)
              ENDDO
              ALBP_ERM(IAL,IGR,JGR)=REAL((1.D0-DGAR1/DGAR2)/(1.D0+
     1        DGAR1/DGAR2))
            ENDDO
          ENDDO
        ENDDO
        CALL LCMPUT(IPMAC,'ALBEDO',NALBP*NGRP*NGRP,2,ALBP_ERM)
        DEALLOCATE(ALBP_ERM)
      ENDIF
*----
*  AVERAGE KINF
*----
      IF(IKINF.EQ.1) THEN
        DGAR1=0.0D0
        DGAR2=0.0D0
        DO IBM=1,NMIX
          DGAR1=DGAR1+ZKINF(IBM)*VOLMI2(IBM)
          DGAR2=DGAR2+VOLMI2(IBM)
        ENDDO
        FACTOR=REAL(DGAR1/DGAR2)
        CALL LCMPUT(IPMAC,'K-INFINITY',1,2,FACTOR)
      ENDIF
*----
*  AVERAGE KEFF
*----
      IF(IKEFF.EQ.1) THEN
        DGAR1=0.0D0
        DGAR2=0.0D0
        DO IBM=1,NMIX
          DGAR1=DGAR1+ZKEFF(IBM)*VOLMI2(IBM)
          DGAR2=DGAR2+VOLMI2(IBM)
        ENDDO
        FACTOR=REAL(DGAR1/DGAR2)
        CALL LCMPUT(IPMAC,'K-EFFECTIVE',1,2,FACTOR)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(AVGFL2,ADF2M,ADF2,HADF,ZKEFF,ZKINF,ALBP2_E,ALBP2,
     1 GAR6,GAR4)
      RETURN
      END
