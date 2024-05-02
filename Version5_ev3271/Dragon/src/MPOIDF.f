*DECK MPOIDF
      SUBROUTINE MPOIDF(IPMPO,IPEDIT,HEDIT,NG,NMIL,ICAL,IDF,NALBP,
     1 FNORM,VOLMIL,FLXMIL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To store discontinuity factor and albedo information in the MPO file.
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
* IPMPO   pointer to the MPO file.
* IPEDIT  pointer to the edition object (L_EDIT signature).
* HEDIT   name of output group for a (multigroup mesh, output geometry)
*         couple (generally equal to 'output_0').
* NG      number of condensed energy groups.
* NMIL    number of mixtures.
* ICAL    index of the current elementary calculation.
* IDF     type of surfacic information (2/3: boundary flux/DF).
* NALBP   number of physical albedos per energy group.
* FNORM   flux normalization factor.
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
      TYPE(C_PTR) IPMPO,IPEDIT
      INTEGER NG,NMIL,ICAL,IDF,NALBP
      REAL FNORM,VOLMIL(NMIL),FLXMIL(NMIL,NG)
      CHARACTER(LEN=12) HEDIT
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*131,RECNAM*80
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: VREAL,SURF
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SURFLX,ALBP
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ALBP2
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HADF
*----
*  CREATE GROUP flux.
*----
      WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,6H/flux/)')
     1 TRIM(HEDIT),ICAL-1
      CALL hdf5_create_group(IPMPO,TRIM(RECNAM))
*----
*  RECOVER DISCONTINUITY FACTOR INFORMATION FROM MACROLIB
*----
      CALL LCMSIX(IPEDIT,'MACROLIB',1)
      CALL LCMLEN(IPEDIT,'ADF',ILONG,ITYLCM)
      IF(ILONG.NE.0) THEN
        CALL LCMSIX(IPEDIT,'ADF',1)
        CALL LCMGET(IPEDIT,'NTYPE',NSURFD)
        IF((IDF.EQ.2).OR.(IDF.EQ.3)) THEN
          ALLOCATE(SURFLX(NSURFD,NMIL*NG),SURF(NMIL*NG))
        ELSE IF(IDF.EQ.4) THEN
          ALLOCATE(SURFLX(NSURFD,NMIL*NG*NG),SURF(NMIL*NG*NG))
        ELSE
          CALL XABORT('MPOIDF: INVALID ADF OPTION.')
        ENDIF
        ALLOCATE(HADF(NSURFD))
        CALL LCMGTC(IPEDIT,'HADF',8,NSURFD,HADF)
        DO I=1,NSURFD
          CALL LCMLEN(IPEDIT,HADF(I),ILONG,ITYLCM)
          IF(IDF.EQ.2) THEN
*           boundary flux information
            IF(ILONG.NE.NMIL*NG) THEN
              WRITE(HSMG,'(16HMPOIDF: INVALID ,A,8H LENGTH=,I5,
     1        10H EXPECTED=,I5,4H.(1))') HADF(I),ILONG,NMIL*NG
              CALL XABORT(HSMG)
            ENDIF
            CALL LCMGET(IPEDIT,HADF(I),SURF)
            DO IMIL=1,NMIL
              DO IGR=1,NG
                IF(FNORM.NE.1.0) THEN
                  SURFLX(I,(IGR-1)*NMIL+IMIL)=SURF((IGR-1)*NMIL+IMIL)*
     1            FNORM*1.0E13
                ELSE
                  SURFLX(I,(IGR-1)*NMIL+IMIL)=SURF((IGR-1)*NMIL+IMIL)
                ENDIF
              ENDDO
            ENDDO
          ELSE IF(IDF.EQ.3) THEN
*           discontinuity factor information
            IF(ILONG.NE.NMIL*NG) THEN
              WRITE(HSMG,'(16HMPOIDF: INVALID ,A,8H LENGTH=,I5,
     1        10H EXPECTED=,I5,4H.(2))') HADF(I),ILONG,NMIL*NG
              CALL XABORT(HSMG)
            ENDIF
            CALL LCMGET(IPEDIT,HADF(I),SURF)
            DO IMIL=1,NMIL
              DO IGR=1,NG
                IOF=(IGR-1)*NMIL+IMIL
                SURFLX(I,IOF)=SURF(IOF)*FLXMIL(IMIL,IGR)
              ENDDO
            ENDDO
          ELSE IF(IDF.EQ.4) THEN
*           matrix discontinuity factor information
            IF(ILONG.NE.NMIL*NG*NG) THEN
              WRITE(HSMG,'(16HMPOIDF: INVALID ,A,8H LENGTH=,I5,
     1        10H EXPECTED=,I5,4H.(3))') HADF(I),ILONG,NMIL*NG*NG
              CALL XABORT(HSMG)
            ENDIF
            CALL LCMGET(IPEDIT,HADF(I),SURF)
            DO IMIL=1,NMIL
              DO IGR=1,NG
                DO JGR=1,NG
                  IOF=((JGR-1)*NG+IGR-1)*NMIL+IMIL
                  SURFLX(I,IOF)=SURF(IOF)*FLXMIL(IMIL,IGR)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        DEALLOCATE(HADF,SURF)
        CALL LCMSIX(IPEDIT,' ',2)
*----
*  MOVE TO THE /statept_id/flux GROUP.
*----
        ALLOCATE(SURF(NMIL*NSURFD))
        SURF(:)=1.0
        CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"NSURF",NSURFD)
        CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"SURF",SURF)
        IF((IDF.EQ.2).OR.(IDF.EQ.3)) THEN
          CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"SURFFLUX",SURFLX)
        ELSE IF(IDF.EQ.4) THEN
          CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"SURFFLUXGxG",SURFLX)
        ENDIF
        ALLOCATE(VREAL(NMIL*NG))
        DO IGR=1,NG
          DO IMIL=1,NMIL
            VREAL((IGR-1)*NMIL+IMIL)=FLXMIL(IMIL,IGR)*VOLMIL(IMIL)
          ENDDO
        ENDDO
        CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"TOTALFLUX",VREAL)
        DEALLOCATE(VREAL,SURF,SURFLX)
      ENDIF
*----
*  RECOVER AND SAVE ALBEDO INFORMATION
*----
      IF(NALBP.NE.0) THEN
        CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"NALBP",NALBP)
        CALL LCMLEN(IPEDIT,'ALBEDO',ILONG,ITYLCM)
        IF(ILONG.EQ.NALBP*NG) THEN
*         diagonal physical albedos
          ALLOCATE(ALBP(NALBP,NG))
          CALL LCMGET(IPEDIT,'ALBEDO',ALBP)
          CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"ALBEDOG",ALBP)
          DEALLOCATE(ALBP)
        ELSE IF(ILONG.EQ.NALBP*NG*NG) THEN
*         matrix physical albedos
          ALLOCATE(ALBP2(NALBP,NG,NG))
          CALL LCMGET(IPEDIT,'ALBEDO',ALBP2)
          CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"ALBEDOGxG",ALBP2)
          DEALLOCATE(ALBP2)
        ELSE
          CALL XABORT('MPOIDF: INCONSISTENT ALBEDO INFORMATION.')
        ENDIF
      ENDIF
      CALL LCMSIX(IPEDIT,' ',2)
      RETURN
      END
