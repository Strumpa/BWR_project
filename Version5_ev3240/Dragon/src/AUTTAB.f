*DECK AUTFLU
      SUBROUTINE AUTTAB(KPLIB,HNAMIS,IGRMIN,IGRRES,NGRP,LBIN,NBIN,UUU,
     1 ISEED,SIGINF,LLL,SIGT,SIGS,SIGF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover resonant Autolib data in the unresolved energy domain.
*
*Copyright:
* Copyright (C) 2023 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* KPLIB   isotope subdirectory in the internal microscopic cross-section
*         library with subgroups.
* HNAMIS  character*12 name of the resonant isotope.
* IGRMIN  first group where the self-shielding is applied.
* IGRRES  first resolved group where the self-shielding is applied.
* NGRP    number of energy groups.
* LBIN    total number of fine energy groups in the Autolib.
* NBIN    number of fine energy groups in each coarse energy group.
* UUU     lethargy limits of the groups.
* ISEED   the seed for the generation of random numbers in the
*         unresolved energy domain.
* SIGINF  infinite dilution x-s values.
*
*Parameters: output
* LLL     number of fine energy groups in the unresolved domain.
* SIGT    total microscopic x-s.
* SIGS    P0 scattering microscopic x-s.
* SIGF    nu*fission microscopic x-s.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) KPLIB
      CHARACTER HNAMIS*12
      INTEGER IGRMIN,IGRRES,NGRP,LBIN,NBIN(NGRP),ISEED,LLL
      REAL UUU(LBIN+1),SIGINF(NGRP,3),SIGT(LBIN),SIGS(LBIN),SIGF(LBIN)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) LPLIB,MPLIB
      PARAMETER(MAXNOR=12)
      CHARACTER HSMG*131
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NOR
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGP
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NOR(NGRP))
*----
*  SET THE RANDOM NUMBER GENERATOR
*----
      IFIRST=1
      IF(ISEED.EQ.0) THEN
         CALL CLETIM(IT1,IT2)
         ISEED=IT1+IT2
         DO 10 JJ=0,MOD(ISEED,10)
         CALL RANDF(ISEED,IFIRST,RAND)
   10    CONTINUE
      ENDIF
*
      CALL LCMLEN(KPLIB,'PT-TABLE',LENG,ITYLCM)
      IF(LENG.EQ.0) THEN
        WRITE(HSMG,'(38HAUTTAB: NO PT-TABLE DATA FOR ISOTOPE '',A12,
     1  23H'' FOR UNRESOLVED GROUPS,2I5,1H.)') HNAMIS,IGRMIN,IGRRES-1
        CALL XABORT(HSMG)
      ENDIF
      CALL LCMSIX(KPLIB,'PT-TABLE',1)
      CALL LCMGET(KPLIB,'NOR',NOR)
      LLL=0
      DO 20 IGRP=1,IGRMIN-1
      LLL=LLL+NBIN(IGRP)
   20 CONTINUE
      LPLIB=LCMGID(KPLIB,'GROUP-PT')
      DO 80 IGRP=IGRMIN,IGRRES-1
      IF(NOR(IGRP).LE.0) THEN
        WRITE(HSMG,'(42HAUTTAB: NO PROBABILITY TABLE DATA IN GROUP,I5,
     1  13H OF ISOTOPE '',A12,2H''.)') IGRP,HNAMIS
        CALL XABORT(HSMG)
      ELSE IF(NBIN(IGRP).LE.0) THEN
        WRITE(HSMG,'(32HAUTTAB: NO AUTOLIB MESH IN GROUP,I5,1H.)') IGRP
        CALL XABORT(HSMG)
      ENDIF
      IF(NOR(IGRP).EQ.1) THEN
        DO 30 IBIN=LLL+1,LLL+NBIN(IGRP)
        SIGT(IBIN)=SIGINF(IGRP,1)
        SIGF(IBIN)=SIGINF(IGRP,2)
        SIGS(IBIN)=SIGINF(IGRP,3)
   30   CONTINUE
      ELSE
        MPLIB=LCMGIL(LPLIB,IGRP)
        CALL LCMLEN(MPLIB,'PROB-TABLE',LENG,ITYLCM)
        NPART=LENG/MAXNOR
        IF(NPART.LT.2) THEN
          CALL LCMLIB(MPLIB)
          CALL XABORT('AUTTAB: SCATTERING INFO MISSING.')
        ENDIF
        DELG=UUU(LLL+NBIN(IGRP)+1)-UUU(LLL+1)
        ALLOCATE(SIGP(MAXNOR*NPART))
        CALL LCMGET(MPLIB,'PROB-TABLE',SIGP)
        ADSIGT=0.0
        ADSIGF=0.0
        ADSIGS=0.0
        DO 60 IBIN=LLL+1,LLL+NBIN(IGRP)
        CALL RANDF(ISEED,IFIRST,RAND)
        WW=0.0
        DO 40 INOR=1,NOR(IGRP)
        WW=WW+SIGP(INOR)
        IF(RAND.LE.WW) THEN
          SIGT(IBIN)=SIGP(MAXNOR+INOR)
          SIGF(IBIN)=SIGP(2*MAXNOR+INOR)
          SIGS(IBIN)=SIGP(3*MAXNOR+INOR)
          GO TO 50
        ENDIF
   40   CONTINUE
        WRITE(HSMG,'(43HAUTTAB: WEIGHT NORMALIZATION ISSUE IN GROUP,I5,
     1  1H.)') IGRP
        CALL XABORT(HSMG)
   50   ADSIGT=ADSIGT+SIGT(IBIN)*(UUU(IBIN+1)-UUU(IBIN))/DELG
        ADSIGF=ADSIGF+SIGF(IBIN)*(UUU(IBIN+1)-UUU(IBIN))/DELG
        ADSIGS=ADSIGS+SIGS(IBIN)*(UUU(IBIN+1)-UUU(IBIN))/DELG
   60   CONTINUE
        FACTT=SIGINF(IGRP,1)/ADSIGT
        IF(ADSIGF.NE.0.0) THEN
          FACTF=SIGINF(IGRP,2)/ADSIGF
        ELSE
          FACTF=0.0
        ENDIF
        FACTS=SIGINF(IGRP,3)/ADSIGS
        DO 70 IBIN=LLL+1,LLL+NBIN(IGRP)
        SIGT(IBIN)=SIGT(IBIN)*FACTT
        SIGF(IBIN)=SIGF(IBIN)*FACTF
        SIGS(IBIN)=SIGS(IBIN)*FACTS
   70   CONTINUE
        DEALLOCATE(SIGP)
      ENDIF
      LLL=LLL+NBIN(IGRP)
   80 CONTINUE
      CALL LCMSIX(KPLIB,' ',2)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(NOR)
      RETURN
      END
