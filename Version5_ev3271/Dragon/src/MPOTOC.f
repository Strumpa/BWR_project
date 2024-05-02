*DECK MPOTOC
      SUBROUTINE MPOTOC(IPMPO,HEDIT,IMPX,NREA,NBISO,NMIL,NPAR,NLOC,
     1 NISOF,NISOP,NISOS,NCAL,NGRP,NSURFD,NALBP,NPRC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the table of content of an MPO file.
*
*Copyright:
* Copyright (C) 2022 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPMPO   address of the MPO file.
* HEDIT   name of output group for a (multigroup mesh, output geometry)
*         couple (generally equal to 'output_0').
* IMPX    print parameter (equal to zero for no print).
*
*Parameters: output
* NREA    number of neutron-induced reaction
* NBISO   number of particularized isotopes
* NMIL    number of mixtures in the MPO file
* NPAR    number of global parameters
* NLOC    number of local parameters
* NISOF   number of particularized fissile isotopes
* NISOP   number of particularized fission products
* NISOS   number of particularized stable isotopes
* NCAL    number of elementary calculations
* NGRP    number of energy groups
* NSURFD  number of discontinuity factors values in the MPO file
* NALBP   number of physical albedos per energy group
* NPRC    number of precursors
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMPO
      INTEGER IMPX,NREA,NBISO,NMIL,NPAR,NLOC,NISOF,NISOP,NISOS,NCAL,
     1 NGRP,NSURFD,NALBP,NPRC
      CHARACTER(LEN=12) HEDIT
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      INTEGER I,J,NENERG,NGEOME,ID_G,ID_E,ID,IBM,NGRP2,RANK,TYPE,NBYTE,
     1 DIMSR(5)
      CHARACTER HSMG*131,RECNAM*80,HFORMAT*132
      CHARACTER(LEN=100), ALLOCATABLE, DIMENSION(:) :: LIST
      INTEGER, ALLOCATABLE, DIMENSION(:) :: DIMS_MPO,ADDRISO
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: OUPUTID
*----
*  LIST GROUPS AND DATASETS ON THE ROOT FILE
*----
      IF(IMPX.GT.0) THEN
        CALL hdf5_list_groups(IPMPO, '/', LIST)
        WRITE(*,*)
        WRITE(*,*) 'MPOTOC: GROUP TABLE OF CONTENTS'
        DO I=1,SIZE(LIST)
          WRITE(*,*) TRIM(LIST(I))
        ENDDO
        DEALLOCATE(LIST)
      ENDIF
*----
*  RECOVER MPO PARAMETERS
*----
      ID_G=-1
      ID_E=-1
      CALL hdf5_read_data(IPMPO,"/parameters/tree/NSTATEPOINT",NCAL)
      CALL hdf5_read_data(IPMPO,"/energymesh/NENERGYMESH",NENERG)
      CALL hdf5_read_data(IPMPO,"/geometry/NGEOMETRY",NGEOME)
      IF((NENERG.GT.0).AND.(NGEOME.GT.0)) THEN
        CALL hdf5_read_data(IPMPO,"/output/OUPUTID",OUPUTID)
        READ(HEDIT,'(7X,I2)') ID
        DO I=1,NGEOME
          DO J=1,NENERG
            IF(OUPUTID(J,I).EQ.ID) THEN
              ID_G=I-1
              ID_E=J-1
              GO TO 10
            ENDIF
          ENDDO
        ENDDO
        CALL XABORT('MPOTOC: no ID found in /output/OUPUTID.')
   10   WRITE(RECNAM,'(23H/energymesh/energymesh_,I0,1H/)') ID_E
        IF(IMPX.GT.1) THEN
          HFORMAT='(/42H MPOTOC: Process MPO multiparameter file o,'//
     >    '9Hn output=,A)'
          WRITE(IOUT,HFORMAT) TRIM(HEDIT)
          WRITE(IOUT,'(24H MPOTOC:   energy group=,A)') TRIM(RECNAM)
        ENDIF
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"NG",NGRP2)
        IF(NGRP.EQ.0) THEN
          NGRP=NGRP2
        ELSE IF(NGRP2.NE.NGRP) THEN
          WRITE(HSMG,'(44H MPOTOC: THE MPO FILE HAS AN INVALID NUMBER ,
     1    18HOF ENERGY GROUPS (,I4,3H VS,I5,2H).)') NGRP2,NGRP
          CALL XABORT(HSMG)
        ENDIF
        DEALLOCATE(OUPUTID)
        WRITE(RECNAM,'(19H/geometry/geometry_,I0,1H/)') ID_G
        IF(IMPX.GT.1) THEN
          WRITE(IOUT,'(24H MPOTOC: geometry group=,A)') TRIM(RECNAM)
        ENDIF
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"NZONE",NMIL)
      ENDIF
      WRITE(RECNAM,'(8H/output/,A,6H/info/)') TRIM(HEDIT)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"NREA",NREA)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ADDRISO",ADDRISO)
      NBISO=ADDRISO(SIZE(ADDRISO,1))
*----
*  SET NPAR
*----
      NPAR=0
      CALL hdf5_info(IPMPO,"/parameters/info/NVALUE",RANK,TYPE,NBYTE,
     1 DIMSR)
      IF(RANK.GT.0) NPAR=DIMSR(1)
*----
*  SET NLOC
*----
      IF(hdf5_group_exists(IPMPO,"/local_values")) THEN
        CALL hdf5_get_shape(IPMPO,"/local_values/LOCVALNAME",DIMS_MPO)
        NLOC=DIMS_MPO(1)
        DEALLOCATE(DIMS_MPO)
      ELSE
        NLOC=0
      ENDIF
*----
*  SET NISOF AND NISOP
*----
      NISOF=0
      NISOP=0
      IF(NBISO.GT.0) THEN
        DO IBM=1,NMIL
          WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,6H/zone_,I0,1H/)')
     1    TRIM(HEDIT),0,IBM-1
          IF(hdf5_group_exists(IPMPO,TRIM(RECNAM)//"yields")) THEN
            CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"yields/NISF",NISOF)
            CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"yields/NISP",NISOP)
            EXIT
          ENDIF
        ENDDO
      ENDIF
      NISOS=NBISO-(NISOF+NISOP)
      DEALLOCATE(ADDRISO)
*----
*  SET NSURFD AND NALBP
*----
      WRITE(RECNAM,'(8H/output/,A,15H/statept_0/flux)') TRIM(HEDIT)
      NSURFD=0
      NALBP=0
      CALL hdf5_info(IPMPO,TRIM(RECNAM)//"/NSURF",RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"/NSURF",
     1 NSURFD)
      CALL hdf5_info(IPMPO,TRIM(RECNAM)//"/NALBP",RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"/NALBP",
     1 NALBP)
*----
*  SET NPRC
*----
      NPRC=0
      WRITE(RECNAM,'(8H/output/,A,27H/statept_0/zone_0/kinetics/)')
     1 TRIM(HEDIT)
      IF(hdf5_group_exists(IPMPO,RECNAM)) THEN
        CALL hdf5_get_shape(IPMPO,TRIM(RECNAM)//"LAMBDAD",DIMS_MPO)
        NPRC=DIMS_MPO(1)
        DEALLOCATE(DIMS_MPO)
      ENDIF
*----
*  PRINT MPO PARAMETERS
*----
      IF(IMPX.GT.0) THEN
        WRITE(IOUT,'(/38H MPOTOC: table of content information:)')
        WRITE(IOUT,'(36H   nb of neutron-induced reactions =,I3)') NREA
        WRITE(IOUT,'(34H   nb of particularized isotopes =,I4)') NBISO
        WRITE(IOUT,'(19H   nb of mixtures =,I5)') NMIL
        WRITE(IOUT,'(28H   nb of global parameters =,I4)') NPAR
        WRITE(IOUT,'(27H   nb of local parameters =,I4)') NLOC
        WRITE(IOUT,'(42H   nb of particularized fissile isotopes =,I4)')
     1  NISOF
        WRITE(IOUT,'(42H   nb of particularized fission products =,I4)')
     1  NISOP
        WRITE(IOUT,'(41H   nb of particularized stable isotopes =,I4)')
     1  NISOS
        WRITE(IOUT,'(23H   nb of calculations =,I9)') NCAL
        WRITE(IOUT,'(24H   nb of energy groups =,I4)') NGRP
        WRITE(IOUT,'(38H   nb of discontinuity factor values =,I4)')
     1  NSURFD
        WRITE(IOUT,'(44H   nb of physical albedos per energy group =,
     1  I4)') NALBP
        WRITE(IOUT,'(21H   nb of precursors =,I4/)') NPRC
      ENDIF
      RETURN
      END
