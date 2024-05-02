*DECK SPHSTM
      SUBROUTINE SPHSTM(IPAPX,ICAL,IMPX,LNEW,HEQUI,HEDIT,NMIL,NGROUP,
     1 SPH)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Store a new set of SPH factors for an elementary calculation in a
* MPO file.
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
* IPAPX   pointer to the MPO file.
* ICAL    index of the elementary calculation being considered.
* IMPX    print parameter (equal to zero for no print).
* LNEW    flag set to .TRUE. to allow the overwriting of the existing
*         set of SPH factors named HEQUI.
* HEQUI   LOCKEY name of SPH-factor set to be stored.
* HEDIT   name of output group for a (multigroup mesh, output geometry)
*         couple (generally equal to 'output_0').
* NMIL    number of mixtures in the elementary calculation.
* NGROUP  number of energy groups in the elementary calculation.
* SPH     SPH-factor set to be stored the MPO file.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPAPX
      INTEGER ICAL,IMPX,NMIL,NGROUP
      REAL SPH(NMIL,NGROUP)
      LOGICAL LNEW
      CHARACTER(LEN=80) HEQUI
      CHARACTER(LEN=12) HEDIT
*----
*  LOCAL VARIABLES
*----
      CHARACTER RECNAM*80
*----
*  SLLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: LOCAD,LOCA_OLD
      REAL, ALLOCATABLE, DIMENSION(:) :: RVALO,RVALO_OLD
      CHARACTER(LEN=80), ALLOCATABLE, DIMENSION(:) :: LOCTYP,LOCKEY,
     & LOCTYP_OLD,LOCKEY_OLD
*----
*  RECOVER MPO FILE CHARACTERISTICS
*----
      NLOC_OLD=0
      IF(hdf5_group_exists(IPAPX,"/local_values/")) THEN
        CALL hdf5_read_data(IPAPX,"local_values/LOCVALTYPE",LOCTYP_OLD)
        CALL hdf5_read_data(IPAPX,"local_values/LOCVALNAME",LOCKEY_OLD)
        NLOC_OLD=SIZE(LOCTYP_OLD,1)
        NLOC=NLOC_OLD
        DO ILOC=1,NLOC_OLD
          IF((LOCTYP_OLD(ILOC).EQ.'EQUI').AND.
     &                  (LOCKEY_OLD(ILOC).EQ.HEQUI)) THEN
*           SET HEQUI EXISTS.
            IF(LNEW) THEN
              IF(IMPX.GT.0) WRITE(6,'(28H SPHSTM: OVERWRITE SPH-FACTO,
     &        12HR SET NAMED ,A)') HEQUI
              JLOC=ILOC
              GO TO 10
            ELSE
              CALL XABORT('SPHSTM: THIS SPH FACTOR SET EXISTS: '//HEQUI)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
      NLOC=NLOC_OLD+1
      JLOC=NLOC
   10 ALLOCATE(LOCTYP(NLOC),LOCKEY(NLOC))
      IF(NLOC_OLD.GT.0) THEN
        LOCTYP(:NLOC_OLD)=LOCTYP_OLD(:NLOC_OLD)
        LOCKEY(:NLOC_OLD)=LOCKEY_OLD(:NLOC_OLD)
        DEALLOCATE(LOCKEY_OLD,LOCTYP_OLD)
      ENDIF
      LOCTYP(JLOC)='EQUI'
      LOCKEY(JLOC)=HEQUI
      CALL hdf5_delete(IPAPX,"local_values/LOCVALTYPE")
      CALL hdf5_delete(IPAPX,"local_values/LOCVALNAME")
      CALL hdf5_write_data(IPAPX,"local_values/LOCVALTYPE",LOCTYP)
      CALL hdf5_write_data(IPAPX,"local_values/LOCVALNAME",LOCKEY)
*----
*  LOOP OVER MIXTURES.
*----
      DO IBM=1,NMIL
        WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,6H/zone_,I0,1H/)')
     &  TRIM(HEDIT),ICAL-1,IBM-1
        CALL hdf5_read_data(IPAPX,TRIM(RECNAM)//"LOCVALADDR",LOCA_OLD)
        CALL hdf5_read_data(IPAPX,TRIM(RECNAM)//"LOCALVALUE",RVALO_OLD)
        ALLOCATE(LOCAD(NLOC+1))
        LOCAD(:NLOC_OLD+1)=LOCA_OLD(:NLOC_OLD+1)
        LOCAD(JLOC+1)=LOCAD(JLOC)+NGROUP
        ALLOCATE(RVALO(LOCAD(NLOC+1)))
        RVALO(:LOCA_OLD(NLOC_OLD+1))=RVALO_OLD(:LOCA_OLD(NLOC_OLD+1))
        DEALLOCATE(LOCA_OLD,RVALO_OLD)
        DO IGR=1,NGROUP
          RVALO(LOCAD(JLOC)+IGR)=SPH(IBM,IGR)
        ENDDO
        CALL hdf5_delete(IPAPX,TRIM(RECNAM)//"LOCVALADDR")
        CALL hdf5_delete(IPAPX,TRIM(RECNAM)//"LOCALVALUE")
        CALL hdf5_write_data(IPAPX,TRIM(RECNAM)//"LOCVALADDR",LOCAD)
        CALL hdf5_write_data(IPAPX,TRIM(RECNAM)//"LOCALVALUE",RVALO)
        DEALLOCATE(RVALO,LOCAD)
      ENDDO
      RETURN
      END
