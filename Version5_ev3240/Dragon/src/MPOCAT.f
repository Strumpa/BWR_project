*DECK MPOCAT
      SUBROUTINE MPOCAT(IPMPO,IPRHS,NPAR,MUPCPO,LGNCPO,LWARN,HEDIT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To catenate a RHS MPO file into the output MPO file.
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
* IPMPO   pointer to the output MPO file.
* IPRHS   pointer to the rhs MPO file (contains the new calculations).
* NPAR    number of global parameters in the output MPO file.
* MUPCPO  tuple of the new global parameters in the output MPO file.
* LGNCPO  LGNEW value of the new global parameters in the output MPO
*         file.
* LWARN   logical used in case if an elementary calculation in the RHS
*         is already present in MPO file. If LWARN=.true. a warning is
*         send and the MPO file values are kept otherwise XABORT is
*         called (default).
* HEDIT   name of output group for a (multigroup mesh, output geometry)
*         couple (generally equal to 'output_0').
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMPO,IPRHS
      INTEGER NPAR,MUPCPO(NPAR)
      LOGICAL LGNCPO(NPAR),LWARN
      CHARACTER(LEN=12) HEDIT
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*131,RECNAM*72,RECNA2*72,TEXT24*24
      LOGICAL EQUAL
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: VINTE,NVALUE,MUPLET,MUPRHS,
     1 MUBASE
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: OUPUTID
      REAL, ALLOCATABLE, DIMENSION(:) :: VREAL
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LGNEW
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: PARFMT
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: VCHAR
      CHARACTER(LEN=24), ALLOCATABLE, DIMENSION(:) :: PARCPO,PARKEY
*----
*  CHECK THE COMPATIBILITY OF THE NEW RHS MPO
*     NCAL is the number of calculations in the output MPO file
*     NCALR is the number of calculations in RHS MPO file
*----
      NGR=0
      CALL MPOTOC(IPRHS,HEDIT,IMPX,NREA,NBISO,NMILR,NPARR,NLOCR,NISOF,
     1 NISOP,NISOS,NCALR,NGR,NSURFD,NALBP,NPRC)
      IF(NCALR.EQ.0) CALL XABORT('MPOCAT: NO CALCULATION IN RHS MPO.')
      NCAL=0
      IF(hdf5_group_exists(IPMPO,"/parameters/tree")) THEN
        CALL hdf5_read_data(IPMPO,"/parameters/tree/NSTATEPOINT",NCAL)
      ENDIF
      IF(NCAL.GT.0) THEN
        NG=0
        CALL MPOTOC(IPMPO,HEDIT,IMPX,NREA,NBISO,NMIL,NPAR1,NLOC,NISOF,
     1  NISOP,NISOS,NCAL,NG,NSURFD,NALBP,NPRC)
        IF(NGR.NE.NG) THEN
          WRITE(HSMG,'(42HMPOCAT: ELEMENTARY CALCULATION WITH AN INV,
     1    20HALIB NB. OF GROUPS =,I7,3H NE,I7,1H.)') NGR,NG
          CALL XABORT(HSMG)
        ELSE IF(NMILR.NE.NMIL) THEN
          WRITE(HSMG,'(42HMPOCAT: ELEMENTARY CALCULATION WITH AN INV,
     1    22HALIB NB. OF MIXTURES =,I7,3H NE,I7,1H.)') NMILR,NMIL
          CALL XABORT(HSMG)
        ELSE IF(NPAR1.NE.NPAR) THEN
          WRITE(HSMG,'(42HMPOCAT: ELEMENTARY CALCULATION WITH AN INV,
     1    31HALIB NB. OF GLOBAL PARAMETERS =,I7,3H NE,I7,1H.)') NPAR1,
     2    NPAR
          CALL XABORT(HSMG)
        ELSE IF(NPARR.GT.NPAR) THEN
          WRITE(HSMG,'(42HMPOCAT: ELEMENTARY CALCULATION WITH AN INV,
     1    31HALIB NB. OF GLOBAL PARAMETERS =,I7,3H GT,I7,1H.)') NPARR,
     2    NPAR
          CALL XABORT(HSMG)
        ELSE IF(NLOCR.NE.NLOC) THEN
          WRITE(HSMG,'(42HMPOCAT: ELEMENTARY CALCULATION WITH AN INV,
     1    30HALIB NB. OF LOCAL PARAMETERS =,I7,3H NE,I7,1H.)') NLOCR,
     2    NLOC
          CALL XABORT(HSMG)
        ENDIF
      ENDIF
*----
*  MAIN LOOP OVER THE NCALR ELEMENTARY CALCULATIONS OF THE RHS MPO
*----
      IDEM=0
      NCALS=NCAL
      DO 140 ICALR=1,NCALR
*----
*  COMPUTE THE MUPLET VECTOR FROM THE RHS MPO
*----
      ALLOCATE(MUPRHS(NPARR),MUPLET(NPAR),LGNEW(NPAR))
      WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0)') TRIM(HEDIT),ICALR-1
      CALL hdf5_read_data(IPRHS,TRIM(RECNAM)//"/PARAMVALUEORD",VINTE)
      IF(SIZE(VINTE).NE.NPARR) THEN
        WRITE(HSMG,'(43HMPOCAT: INCONSISTENT PARAMVALUEORD LENGTH (,
     1  I5,3H VS,I5,2H).)') SIZE(VINTE),NPARR
        CALL XABORT(HSMG)
      ENDIF
      DO 20 IPAR=1,NPARR
      MUPRHS(IPAR)=VINTE(IPAR)
   20 CONTINUE
      DEALLOCATE(VINTE)
*----
*  RECOVER THE GLOBAL PARAMETERS
*----
      DO 30 I=1,NPAR
      MUPLET(I)=MUPCPO(I)
      LGNEW(I)=LGNCPO(I)
   30 CONTINUE
      IF(NPAR.GT.0) THEN
        CALL hdf5_read_data(IPMPO,"/parameters/info/PARAMNAME",PARCPO)
      ENDIF
      CALL hdf5_read_data(IPRHS,"/parameters/info/PARAMNAME",PARKEY)
      CALL hdf5_read_data(IPRHS,"/parameters/info/PARAMFORM",PARFMT)
      CALL hdf5_read_data(IPRHS,"/parameters/info/NVALUE",NVALUE)
      DO 100 IPAR=1,NPARR
        DO 80 I0=1,NPAR
        IF(PARKEY(IPAR).EQ.PARCPO(I0)) THEN
          IPARN=I0
          GO TO 90
        ENDIF
   80   CONTINUE
        CALL XABORT('MPOCAT: UNABLE TO FIND '//PARKEY(IPAR)//'.')
   90   WRITE(RECNAM,'(25H/parameters/values/PARAM_,I0)') IPAR-1
        IVAL=MUPRHS(IPAR)+1
        IF(PARFMT(IPAR).EQ.'FLOAT') THEN
          CALL hdf5_read_data(IPRHS,RECNAM,VREAL)
          FLOTT=VREAL(IVAL)
          DEALLOCATE(VREAL)
        ELSE IF(PARFMT(IPAR).EQ.'INTEGER') THEN
          CALL hdf5_read_data(IPRHS,RECNAM,VINTE)
          NITMA=VINTE(IVAL)
          DEALLOCATE(VINTE)
        ELSE IF(PARFMT(IPAR).EQ.'STRING') THEN
          CALL hdf5_read_data(IPRHS,RECNAM,VCHAR)
          TEXT24=VCHAR(IVAL)
          DEALLOCATE(VCHAR)
        ENDIF
        CALL MPOPAV(IPMPO,HEDIT,IPARN,NPAR,PARFMT(IPAR),FLOTT,NITMA,
     1  TEXT24,MUPLET(IPARN),LGNEW(IPARN))
  100 CONTINUE
      IF(NPAR.GT.0) DEALLOCATE(PARCPO)
*----
*  CHECK IF THE NEW MUPLET ALREADY EXISTS IN THE MPO FILE
*----
      DO 120 ICAL=1,NCALS
        WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,14H/PARAMVALUEORD)')
     1  TRIM(HEDIT),ICAL-1
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM),MUBASE)
        EQUAL=.FALSE.
        DO 110 I=1,NPAR
          EQUAL=(MUPLET(I).EQ.MUBASE(I))
          IF(.NOT.EQUAL) GO TO 115
  110   CONTINUE
        WRITE(6,'(/25H MPOCAT: ADD CALCULATION=,I8/8H MUPLET=,20I7:/
     1  (8X,20I7))') ICALR-1,MUPLET(:)
        WRITE(6,'(/48H MPOCAT: ELEMENTARY CALCULATION HAS THE SAME PAR,
     1  36HAMETERS AS ELEMENTARY CALCULATION NB,I8,1H.)') ICAL-1
        IF(LWARN) THEN
          IDEM=IDEM+1
          DEALLOCATE(MUBASE)
          GOTO 130
        ELSE
          CALL XABORT('MPOCAT: WARNING-ONLY FLAG NOT SET.')
        ENDIF
  115   DEALLOCATE(MUBASE)
  120 CONTINUE
      NCALS=NCALS+1
      IF(NCALS.NE.NCAL+ICALR-IDEM) CALL XABORT('MPOCAT: INVALID NCALS.')
*----
*  RECOVER THE ELEMENTARY CALCULATION
*----
      IF(NCALS.EQ.1) THEN
        CALL hdf5_create_group(IPMPO,"/output")
        CALL hdf5_read_data(IPRHS,"/output/NOUTPUT",NOUTPUT)
        CALL hdf5_write_data(IPMPO,"/output/NOUTPUT",NOUTPUT)
        CALL hdf5_read_data(IPRHS,"/output/OUPUTID",OUPUTID)
        CALL hdf5_write_data(IPMPO,"/output/OUPUTID",OUPUTID)
        DEALLOCATE(OUPUTID)
        WRITE(RECNAM,'(8H/output/,A)') TRIM(HEDIT)
        CALL hdf5_create_group(IPMPO,TRIM(RECNAM))
        CALL hdf5_delete(IPMPO,"/energymesh")
        call hdf5_copy(IPRHS,"/energymesh",IPMPO,"/energymesh") ! IPRHS -> IPMPO
        CALL hdf5_delete(IPMPO,"/geometry")
        call hdf5_copy(IPRHS,"/geometry",IPMPO,"/geometry") ! IPRHS -> IPMPO
        WRITE(RECNAM,'(8H/output/,A,5H/info)') TRIM(HEDIT)
        call hdf5_copy(IPRHS,TRIM(RECNAM),IPMPO,TRIM(RECNAM)) ! IPRHS -> IPMPO
      ENDIF
      WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0)') TRIM(HEDIT),NCALS-1
      WRITE(RECNA2,'(8H/output/,A,9H/statept_,I0)') TRIM(HEDIT),ICALR-1
      call hdf5_copy(IPRHS,RECNA2,IPMPO,RECNAM) ! IPRHS -> IPMPO
      CALL hdf5_delete(IPMPO,TRIM(RECNAM)//"/PARAMVALUEORD")
      CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"/PARAMVALUEORD",MUPLET)
      CALL hdf5_delete(IPMPO,"/parameters/tree/NSTATEPOINT")
      CALL hdf5_write_data(IPMPO,"/parameters/tree/NSTATEPOINT",NCALS)
  130 DEALLOCATE(LGNEW,MUPLET,MUPRHS,NVALUE,PARFMT,PARKEY)
  140 CONTINUE
* END OF LOOP ON RHS ELEMENTARY CALCULATIONS. ********************
      RETURN
      END
