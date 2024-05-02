*DECK MPOPAV
      SUBROUTINE MPOPAV(IPMPO,HEDIT,IPAR,NPAR,TTYPE,RVAL,IVAL,CVAL,IV,
     1 LGNEW)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To return the index of a global parameter value. Reorganize the
* parameters group if required.
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
* HEDIT   name of output group for a (multigroup mesh, output geometry)
*         couple (generally equal to 'output_0').
* IPAR    index of the global parameter.
* NPAR    total number of global parameters.
* TTYPE   type of the global parameter value.
* RVAL    global parameter value if TTYPE='FLOAT'.
* IVAL    global parameter value if TTYPE='INTEGER'.
* CVAL    global parameter value if TTYPE='STRING'.
*
*Parameters: output
* IV      index of the global parameter value (IV >= 0).
* LGNEW   new parameter flag (=.true. if the parameter value is new).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMPO
      INTEGER IPAR,NPAR,IV,IVAL
      REAL RVAL
      LOGICAL LGNEW,LSHIFT
      CHARACTER HEDIT*12,TTYPE*8,CVAL*(*)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (REPS=1.0E-5)
      INTEGER RANK,TYPE,NBYTE,DIMSR(5)
      CHARACTER RECNAM*72
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NVALUE,VINTE,VINTE_OLD,
     1 DIMS_MPO,MUPLET
      REAL, ALLOCATABLE, DIMENSION(:) :: VREAL,VREAL_OLD
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: VCHAR,VCHAR_OLD
*
      IF(IPAR.GT.NPAR) CALL XABORT('MPOPAV: NPAR OVERFLOW.')
      CALL hdf5_read_data(IPMPO,"/parameters/info/NVALUE",NVALUE)
      WRITE(RECNAM,'(25H/parameters/values/PARAM_,I0)') IPAR-1
*
      LGNEW=.TRUE.
      LSHIFT=.FALSE.
      IF(TTYPE.EQ.'FLOAT') THEN
         IF(NVALUE(IPAR).EQ.0) THEN
           ALLOCATE(VREAL(1))
           IV=0
           VREAL(IV+1)=RVAL
           NVALUE(IPAR)=1
         ELSE
           CALL hdf5_get_shape(IPMPO,TRIM(RECNAM),DIMS_MPO)
           ILONG=DIMS_MPO(1)
           DEALLOCATE(DIMS_MPO)
           IF(ILONG.GT.NVALUE(IPAR)) CALL XABORT('MPOPAV: NVALUE OVER'
     1     //'FLOW(1).')
           CALL hdf5_read_data(IPMPO,TRIM(RECNAM),VREAL_OLD)
           DO 10 I=1,NVALUE(IPAR)
             IF(RVAL.LE.VREAL_OLD(I)*(1.+REPS))THEN
               IV=I-1
               LGNEW=RVAL.LT.VREAL_OLD(IV+1)*(1.-REPS)
               GO TO 20
             ENDIF
   10      CONTINUE
           IV=NVALUE(IPAR)
   20      ALLOCATE(VREAL(NVALUE(IPAR)+1))
           VREAL(:NVALUE(IPAR))=VREAL_OLD(:NVALUE(IPAR))
           IF(LGNEW) THEN
             LSHIFT=IV.LT.NVALUE(IPAR)
             NVALUE(IPAR)=NVALUE(IPAR)+1
             DO 30 J=NVALUE(IPAR)-1,IV+1,-1
                VREAL(J+1)=VREAL_OLD(J)
   30        CONTINUE
             VREAL(IV+1)=RVAL
           ENDIF
           DEALLOCATE(VREAL_OLD)
         ENDIF
         IF(LGNEW) THEN
           CALL hdf5_info(IPMPO,TRIM(RECNAM),RANK,TYPE,NBYTE,DIMSR)
           IF(TYPE.NE.99) CALL hdf5_delete(IPMPO,TRIM(RECNAM))
           CALL hdf5_write_data(IPMPO,TRIM(RECNAM),VREAL)
         ENDIF
         DEALLOCATE(VREAL)
      ELSE IF(TTYPE.EQ.'INTEGER') THEN
         IF(NVALUE(IPAR).EQ.0) THEN
           ALLOCATE(VINTE(1))
           IV=0
           VINTE(IV+1)=IVAL
           NVALUE(IPAR)=1
         ELSE
           CALL hdf5_get_shape(IPMPO,TRIM(RECNAM),DIMS_MPO)
           ILONG=DIMS_MPO(1)
           DEALLOCATE(DIMS_MPO)
           IF(ILONG.GT.NVALUE(IPAR)) CALL XABORT('MPOPAV: NVALUE OVER'
     1     //'FLOW(2).')
           CALL hdf5_read_data(IPMPO,TRIM(RECNAM),VINTE_OLD)
           DO 40 I=1,NVALUE(IPAR)
             IF(IVAL.LE.VINTE_OLD(I))THEN
               IV=I-1
               LGNEW=IVAL.LT.VINTE_OLD(IV+1)
               GO TO 50
             ENDIF
   40      CONTINUE
           IV=NVALUE(IPAR)
   50      ALLOCATE(VINTE(NVALUE(IPAR)+1))
           VINTE(:NVALUE(IPAR))=VINTE_OLD(:NVALUE(IPAR))
           IF(LGNEW) THEN
             NVALUE(IPAR)=NVALUE(IPAR)+1
             DO 60 J=NVALUE(IPAR)-1,IV+1,-1
               VINTE(J+1)=VINTE_OLD(J)
   60        CONTINUE
             VINTE(IV+1)=IVAL
           ENDIF
           DEALLOCATE(VINTE_OLD)
         ENDIF
         IF(LGNEW) THEN
           CALL hdf5_info(IPMPO,TRIM(RECNAM),RANK,TYPE,NBYTE,DIMSR)
           IF(TYPE.NE.99) CALL hdf5_delete(IPMPO,TRIM(RECNAM))
           CALL hdf5_write_data(IPMPO,TRIM(RECNAM),VINTE)
         ENDIF
         DEALLOCATE(VINTE)
      ELSE IF(TTYPE.EQ.'STRING') THEN
         IF(NVALUE(IPAR).EQ.0) THEN
           ALLOCATE(VCHAR(1))
           IV=0
           VCHAR(IV+1)=CVAL
           NVALUE(IPAR)=1
         ELSE
           CALL hdf5_get_shape(IPMPO,TRIM(RECNAM),DIMS_MPO)
           ILONG=DIMS_MPO(1)
           DEALLOCATE(DIMS_MPO)
           IF(ILONG.GT.NVALUE(IPAR)) CALL XABORT('MPOPAV: NVALUE OVER'
     1     //'FLOW(3).')
           CALL hdf5_read_data(IPMPO,TRIM(RECNAM),VCHAR_OLD)
           DO 70 I=1,NVALUE(IPAR)
             IF(CVAL.EQ.VCHAR_OLD(I))THEN
               IV=I-1
               LGNEW=.FALSE.
               GO TO 80
             ENDIF
   70      CONTINUE
           IV=NVALUE(IPAR)
   80      ALLOCATE(VCHAR(NVALUE(IPAR)+1))
           VCHAR(:NVALUE(IPAR))=VCHAR_OLD(:NVALUE(IPAR))
           IF(LGNEW) THEN
             NVALUE(IPAR)=NVALUE(IPAR)+1
             VCHAR(NVALUE(IPAR))=CVAL
           ENDIF
           DEALLOCATE(VCHAR_OLD)
         ENDIF
         IF(LGNEW) THEN
           CALL hdf5_info(IPMPO,TRIM(RECNAM),RANK,TYPE,NBYTE,DIMSR)
           IF(TYPE.NE.99) CALL hdf5_delete(IPMPO,TRIM(RECNAM))
           CALL hdf5_write_data(IPMPO,TRIM(RECNAM),VCHAR)
         ENDIF
         DEALLOCATE(VCHAR)
      ELSE
         CALL XABORT('MPOPAV: UNKNOWN TYPE='//TTYPE//'.')
      ENDIF
*
      IF(LGNEW) THEN
        CALL hdf5_delete(IPMPO,"/parameters/info/NVALUE")
        CALL hdf5_write_data(IPMPO,"/parameters/info/NVALUE",NVALUE)
      ENDIF
      IF(LSHIFT) THEN
        CALL hdf5_read_data(IPMPO,"/parameters/tree/NSTATEPOINT",NCALAR)
        DO 90 ICAL=1,NCALAR
        WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,14H/PARAMVALUEORD)')
     1  TRIM(HEDIT),ICAL-1
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM),MUPLET)
        IF(MUPLET(IPAR).GE.IV) THEN
          MUPLET(IPAR)=MUPLET(IPAR)+1
          CALL hdf5_delete(IPMPO,TRIM(RECNAM))
          CALL hdf5_write_data(IPMPO,TRIM(RECNAM),MUPLET)
        ENDIF
        DEALLOCATE(MUPLET)
  90    CONTINUE
      ENDIF
      DEALLOCATE(NVALUE)
      RETURN
      END
