*DECK LIBE3R
      SUBROUTINE LIBE3R(CFILNA1,CFILNA2,MAXR,NEL,NBESP,IMPX,ITNAM,
     1 ITZEA,KPAX,BPAX,ENER)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read depletion data on an APOLIB-3 formatted library.
*
*Copyright:
* Copyright (C) 2022 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert
*
*Parameters: input
* CFILNA1 APOLIB-3 cross section file name.
* CFILNA2 APOLIB-3 depletion file name.
* MAXR    number of reaction types.
* NEL     number of isotopes on library.
* NBESP   number of energy-dependent fission yield matrices.
* IMPX    print flag.
*
*Parameters: output
* ITNAM   reactive isotope names in chain.
* ITZEA   6-digit nuclide identifier:
*         atomic number z*10000 (digits) + mass number a*10 +
*         energy state (0 = ground state, 1 = first state, etc.).
* KPAX    complete reaction type matrix.
* BPAX    complete branching ratio matrix.
* ENER    output energy mesh corresponding to a yield macrogroup.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER CFILNA1*(*),CFILNA2*(*)
      INTEGER MAXR,NEL,NBESP,IMPX,ITNAM(3,NEL),ITZEA(NEL),
     1 KPAX(NEL+MAXR,NEL)
      REAL BPAX(NBESP,NEL+MAXR,NEL),ENER(NBESP+1)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPAP1,IPAP2
      DOUBLE PRECISION SUM
      PARAMETER (IOUT=6,MAXR2=12)
      PARAMETER (KDECAY=1,KFISSP=2,KCAPTU=3,KN2N=4,KN3N=5,KN4N=6)
      CHARACTER RECNAM*80,RECNAM2*80,HSMG*131,NMDEPA(MAXR2)*6
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IA,IZ
      REAL, ALLOCATABLE, DIMENSION(:) :: BRANCH,FSYIELDS,YIELDEN
      CHARACTER(LEN=100), ALLOCATABLE, DIMENSION(:) :: LIST
      CHARACTER(LEN=24), ALLOCATABLE, DIMENSION(:) :: NAMES,DEPLNAMES,
     > FSISNAMES,FSPRNAMES,PNAMES,REACID
*----
*  DATA STATEMENTS
*----
      SAVE      NMDEPA
      DATA      NMDEPA/'DECAY ','NUFISS','MT102 ','MT16  ',
     >                 'MT17  ','MT37  ','MT107 ','MT103 ',
     >                 'MT108 ','MT28  ','MT104 ','MT105 '/
*----
*  OPEN APOLIB FILES
*----
      CALL hdf5_open_file(CFILNA1,IPAP1,.true.)
      IF(IMPX.GT.1) THEN
        CALL hdf5_list_groups(IPAP1,"/",LIST)
        WRITE(*,*)
        WRITE(*,*) 'LIBE3R: GROUP TABLE OF CONTENTS FOR FILE ',CFILNA1
        DO I=1,SIZE(LIST)
          WRITE(*,*) TRIM(LIST(I))
        ENDDO
        DEALLOCATE(LIST)
      ENDIF
      IPAP2=C_NULL_PTR
      IF(CFILNA2.NE.' ') THEN
        CALL hdf5_open_file(CFILNA2,IPAP2,.true.)
        IF(IMPX.GT.1) THEN
          CALL hdf5_list_groups(IPAP2,"/",LIST)
          WRITE(*,*)
          WRITE(*,*) 'LIBE3R: GROUP TABLE OF CONTENTS FOR FILE ',CFILNA2
          DO I=1,SIZE(LIST)
            WRITE(*,*) TRIM(LIST(I))
          ENDDO
          DEALLOCATE(LIST)
        ENDIF
      ENDIF
*----
*  RECOVER INFORMATION FROM Head AND PhysicalData GROUPS IN IPAP1
*----
      CALL hdf5_read_data(IPAP1,"Head/nbIs",NISOT)
      IF(NISOT.NE.NEL) CALL XABORT('LIBE3R: INVALID VALUE OF NEL.')
      CALL hdf5_read_data(IPAP1,"Head/IsNames",NAMES)
      DO ISO=1,NISOT
        READ(NAMES(ISO),'(3A4)') (ITNAM(II,ISO),II=1,3)
      ENDDO
      CALL hdf5_read_data(IPAP1,"PhysicalData/MassNb",IA)
      CALL hdf5_read_data(IPAP1,"PhysicalData/AtomicNb",IZ)
*----
*  RECOVER INFORMATION FROM "Yields GROUP IN IPAP2
*----
      IF(C_ASSOCIATED(IPAP2)) THEN
        CALL hdf5_list_groups(IPAP2,"Chain",DEPLNAMES)
        CALL hdf5_read_data(IPAP2,"Yields/FsIsNames",FSISNAMES)
        CALL hdf5_read_data(IPAP2,"Yields/FsPrNames",FSPRNAMES)
        CALL hdf5_read_data(IPAP2,"Yields/FsYields",FSYIELDS)
        CALL hdf5_read_data(IPAP2,"Yields/YieldEnMshInMeV",YIELDEN)
        NDEP=SIZE(DEPLNAMES)
        NDFI=SIZE(FSISNAMES)
        NDFP=SIZE(FSPRNAMES)
*----
*  FISSION YIELD NORMALIZATION
*----
        DO IGF=1,NBESP
          DO IFI=1,NDFI
            SUM=0.0D0
            DO IFP=1,NDFP
              IOF=((IFP-1)*NDFI+IFI-1)*NBESP+IGF
              SUM=SUM+FSYIELDS(IOF)
            ENDDO
            DO IFP=1,NDFP
              IOF=((IFP-1)*NDFI+IFI-1)*NBESP+IGF
              FSYIELDS(IOF)=2.0*FSYIELDS(IOF)/REAL(SUM)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        NDEP=0
        NDFI=0
        NDFP=0
      ENDIF
*----
*  MAIN LOOP OVER ISOTOPES
*----
      NDFP2=0
      DO ISO=1,NISOT
        ITZEA(ISO)=IZ(ISO)*10000+IA(ISO)*10
        II=LEN(TRIM(NAMES(ISO)))
        IF(NAMES(ISO)(II:II).EQ."M") ITZEA(ISO)=ITZEA(ISO)+1
*----
*  DECAY AND BURNOUT OF ISOTOPE ISO
*----
        WRITE(RECNAM,'(10HIsotopeXS/,A,10H/DecayData)') TRIM(NAMES(ISO))
        IF(hdf5_group_exists(IPAP1,TRIM(RECNAM))) THEN
          CALL hdf5_read_data(IPAP1,TRIM(RECNAM)//"/Lambda",DECAY)
          IF(DECAY.NE.0.0) THEN
            KPAX(NEL+KDECAY,ISO)=1
            BPAX(:,NEL+KDECAY,ISO)=DECAY*1.E8
          ENDIF
        ENDIF
        WRITE(RECNAM,'(10HIsotopeXS/,A,12H/ReactionXS/)')
     1  TRIM(NAMES(ISO))
        IF(hdf5_group_exists(IPAP1,TRIM(RECNAM))) THEN
          CALL hdf5_list_datasets(IPAP1,TRIM(RECNAM),LIST)
          DO IREAC=2,MAXR
            II=LEN(TRIM(NMDEPA(IREAC)))
            DO I=1,SIZE(LIST)
              IF(LIST(I)(:II).EQ.NMDEPA(IREAC)) THEN
                KPAX(NEL+IREAC,ISO)=1
                EXIT
              ENDIF
            ENDDO
          ENDDO
          DEALLOCATE(LIST)
        ENDIF
        WRITE(RECNAM,'(10HIsotopeXS/,A,8H/Energy/)') TRIM(NAMES(ISO))
        IF(hdf5_group_exists(IPAP1,TRIM(RECNAM))) THEN
          IF(KPAX(NEL+KFISSP,ISO).EQ.1) THEN
            WRITE(RECNAM2,'(A,16HFISS/EnergyValue)') TRIM(RECNAM)
            CALL hdf5_read_data(IPAP1,TRIM(RECNAM2),VALUE)
            BPAX(:,NEL+KFISSP,ISO)=VALUE
          ENDIF
          IF(KPAX(NEL+KCAPTU,ISO).EQ.1) THEN
            WRITE(RECNAM2,'(A,18HMT-102/EnergyValue)') TRIM(RECNAM)
            CALL hdf5_read_data(IPAP1,TRIM(RECNAM2),VALUE)
            BPAX(:,NEL+KCAPTU,ISO)=VALUE
          ENDIF
        ENDIF
        IF(IMPX.GT.2) THEN
          WRITE(IOUT,100) NAMES(ISO),BPAX(1,NEL+KDECAY,ISO),
     1    BPAX(1,NEL+KFISSP,ISO),BPAX(1,NEL+KCAPTU,ISO)
          WRITE(IOUT,110) (NMDEPA(I),KPAX(NEL+I,ISO),I=1,MAXR)
        ENDIF
*----
*  PARENT REACTIONS OF ISOTOPE ISO
*----
        IF(.NOT.C_ASSOCIATED(IPAP2)) CYCLE
        DO I=1,NDEP
          IF(NAMES(ISO).EQ.DEPLNAMES(I)) GO TO 10
        ENDDO
        GO TO 25
   10   WRITE(RECNAM,'(6HChain/,A,9H/NBPARENT)') TRIM(NAMES(ISO))
        CALL hdf5_read_data(IPAP2, RECNAM, NBPAR)
        IF(NBPAR.EQ.0) GO TO 25
        WRITE(RECNAM,'(6HChain/,A,1H/)') TRIM(NAMES(ISO))
        CALL hdf5_read_data(IPAP2,TRIM(RECNAM)//"BRANCHRATIO",BRANCH)
        CALL hdf5_read_data(IPAP2,TRIM(RECNAM)//"PARENTNAME",PNAMES)
        CALL hdf5_read_data(IPAP2,TRIM(RECNAM)//"REACTIONID",REACID)
        IF(IMPX.GT.2) WRITE(IOUT,120) (PNAMES(IPAR),IPAR=1,NBPAR)
        DO IPAR=1,NBPAR
          JSO=0
          DO I=1,NISOT
            IF(PNAMES(IPAR).EQ.NAMES(I)) THEN
              JSO=I
              GO TO 20
            ENDIF
          ENDDO
          WRITE(HSMG,'(38HLIBE3R: UNABLE TO FIND PARENT ISOTOPE ,A,
     1    8H OF SON ,A,19H IN DEPLETION LIST.)') TRIM(PNAMES(IPAR)),
     2    TRIM(NAMES(ISO))
          CALL XABORT(HSMG)
   20     IF(REACID(IPAR)(:5).EQ.'DRTYP') THEN
            KPAX(ISO,JSO)=KDECAY
          ELSE IF(REACID(IPAR).EQ.'REAMT102') THEN
            KPAX(ISO,JSO)=KCAPTU
          ELSE IF(REACID(IPAR).EQ.'REAMT16') THEN
            KPAX(ISO,JSO)=KN2N
          ELSE IF(REACID(IPAR).EQ.'REAMT17') THEN
            KPAX(ISO,JSO)=KN3N
          ELSE IF(REACID(IPAR).EQ.'REAMT37') THEN
            KPAX(ISO,JSO)=KN4N
          ELSE
            WRITE(HSMG,'(36HLIBE3R: UNKNOWN PRODUCTION REACTION ,A)')
     1      TRIM(REACID(IPAR))
            CALL XABORT(HSMG)
          ENDIF
          BPAX(:,ISO,JSO)=BRANCH(IPAR)
        ENDDO
        DEALLOCATE(REACID,PNAMES,BRANCH)
*----
*  FISSION YIELD OF ISOTOPE ISO
*----
   25   IFP=0
        DO I=1,NDFP
          IF(NAMES(ISO).EQ.FSPRNAMES(I)) THEN
            IFP=I
            GO TO 30
          ENDIF
        ENDDO
        GO TO 50
   30   DO IPAR=1,NDFI
          JSO=0
          DO I=1,NISOT
            IF(FSISNAMES(IPAR).EQ.NAMES(I)) THEN
              JSO=I ! fissile isotope
              GO TO 40
            ENDIF
          ENDDO
          WRITE(HSMG,'(39HLIBE3R: UNABLE TO FIND FISSILE ISOTOPE ,A,
     1    8H OF SON ,A,19H IN DEPLETION LIST.)') TRIM(PNAMES(IPAR)),
     2    TRIM(NAMES(ISO))
          CALL XABORT(HSMG)
   40     KPAX(ISO,JSO)=KFISSP
          DO I=1,NBESP
            IOF=((IFP-1)*NDFI+IPAR-1)*NBESP+I
            BPAX(I,ISO,JSO)=FSYIELDS(IOF)
            ENER(I)=YIELDEN(I)*1.E6
          ENDDO
          ENER(NBESP+1)=YIELDEN(NBESP+1)*1.E6
        ENDDO
   50   DO IPAR=1,NDFP
          IF(FSPRNAMES(IPAR).EQ.NAMES(ISO)) THEN
            NDFP2=NDFP2+1
            GO TO 60
          ENDIF
        ENDDO
   60   CONTINUE
      ENDDO
      IF(NDFP2.NE.NDFP) CALL XABORT('LIBE3R: MISSING FISSION PRODUCT.')
      IF(C_ASSOCIATED(IPAP2)) DEALLOCATE(YIELDEN,FSYIELDS,FSPRNAMES,
     1 FSISNAMES,DEPLNAMES)
      DEALLOCATE(IZ,IA,NAMES)
      CALL hdf5_close_file(IPAP1)
      CALL hdf5_close_file(IPAP2)
*----
*  FIND FISSION PRODUCTS
*----
      DO ISO=1,NISOT
        DO JSO=1,NISOT
          IF(KPAX(JSO,ISO).EQ.KFISSP) KPAX(NEL+KFISSP,JSO)=-1
        ENDDO
      ENDDO
      RETURN
*
  100 FORMAT(/44H LIBE3R: DECAY AND BURNOUT DATA FOR ISOTOPE=,A/
     1 5X,6HDECAY=,1P,E12.5,7H E-8 /S,16H FISSION ENERGY=,E12.5,4H MEV,
     1 16H CAPTURE ENERGY=,E12.5,4H MEV)
  110 FORMAT(5X,12(A6,2H= ,I1,2X))
  120 FORMAT(5X,14HPARENT NAMES: ,12A8/(19X,12A8))
      END
