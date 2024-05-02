*DECK USSIT4
      SUBROUTINE USSIT4(MAXNOR,IPLI0,IPPT1,IPPT2,NGRP,NIRES,NBNRS,NL,
     1 NED,NDEL,PHGAR,STGAR,SFGAR,SSGAR,S0GAR,SAGAR,SDGAR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute self-shielded microscopic cross sections for the RSE method.
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
* MAXNOR  maximum number of base points.
* IPLI0   pointer to the internal microscopic cross section library
*         builded by the self-shielding module.
* IPPT1   pointer to LCM directory of each resonant isotope.
* IPPT2   information related to each resonant isotope:
*         IPPT2(:,1) index of a resonant region (used with infinite
*         dilution case);
*         IPPT2(:,2:4) alias name of resonant isotope;
*         IPPT2(:,5) number of delayed neutron groups.
* NGRP    number of energy groups.
* NIRES   exact number of resonant isotopes.
* NBNRS   number of correlated fuel regions.
* NL      number of Legendre orders required in the calculation
*         (NL=1 or higher).
* NED     number of extra vector edits.
* NDEL    number of delayed neutron precursor groups.
*
*Parameters: output
* PHGAR   averaged flux.
* STGAR   averaged microscopic total xs in resonant region.
* SFGAR   averaged nu*microscopic fission xs in resonant region.
* SSGAR   averaged microscopic scattering xs in resonant region.
* S0GAR   averaged microscopic transfer scattering xs in resonant
*         region for primary neutrons in current group.
* SAGAR   averaged microscopic self-shielded additional xs.
* SDGAR   microscopic self-shielded delayed nu-sigf xs.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLI0,IPPT1(NIRES)
      INTEGER MAXNOR,IPPT2(NIRES,5),NGRP,NIRES,NBNRS,NL,NED,NDEL
      REAL PHGAR(NBNRS,NIRES,NGRP),STGAR(NBNRS,NIRES,NGRP),
     1 SFGAR(NBNRS,NIRES,NGRP),SSGAR(NBNRS,NIRES,NL,NGRP),
     2 S0GAR(NBNRS,NIRES,NL,NGRP,NGRP),SAGAR(NBNRS,NIRES,NED,NGRP),
     3 SDGAR(NBNRS,NIRES,NDEL,NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(MAX_R=12)
      TYPE(C_PTR) IPLIB,JPLIB1,JPLIB2,KPLIB,JPLI0
      CHARACTER HSMG*131,TEXT12*12,CBDPNM*12
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MRANK
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISM
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIGP_R
      REAL, ALLOCATABLE, DIMENSION(:,:) :: XFLUX
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: CGAR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: SIGP
*----
*  SCRATCH STORAGE ALLOCATION.
*----
      ALLOCATE(XFLUX(NBNRS,MAXNOR),ISM(2,NL),MRANK(NGRP))
*----
*  RECOVER INFORMATION FROM THE INTERNAL MICROSCOPIC LIBRARY.
*----
      DO IRES=1,NIRES
        IPLIB=IPPT1(IRES)
        CALL LCMLEN(IPLIB,'NOR',ILONG,ITYLCM)
        IF(ILONG.NE.NGRP) THEN
          CALL LCMLIB(IPLIB)
          CALL XABORT('USSIT4: RANK ARRAY MISSING.')
        ENDIF
        JPLIB1=LCMGID(IPLIB,'GROUP-RSE')
        JPLIB2=LCMGID(IPLIB,'GROUP-PT')
*----
*  WEIGHT DILUTION-DEPENDENT DATA.
*----
        CALL LCMGET(IPLIB,'NOR',MRANK)
        WRITE(CBDPNM,'(3HCOR,I4.4,1H/,I4.4)') IRES,NIRES
        CALL LCMSIX(IPLI0,CBDPNM,1)
        CALL LCMLEN(IPLI0,'NWT0-PT',ILONG,ITYLCM)
        JPLI0=LCMGID(IPLI0,'NWT0-PT')
        DO IGRP=1,NGRP
          MI=MRANK(IGRP)
          IF(MI.LE.1) CYCLE
          CALL LCMLEL(JPLI0,IGRP,ILONG,ITYLCM)
          IF(ILONG.EQ.0) CYCLE
          CALL LCMGDL(JPLI0,IGRP,XFLUX)
          CALL LCMLEL(JPLIB1,IGRP,ILONG1,ITYLCM)
          CALL LCMLEL(JPLIB2,IGRP,ILONG2,ITYLCM)
          IF((ILONG1.EQ.-1).AND.(ILONG2.EQ.0)) THEN
            ! recover a RSE table
            KPLIB=LCMGIL(JPLIB1,IGRP)
            CALL LCMLEN(KPLIB,'RSE-TABLE',ILONG,ITYLCM)
            IF(ILONG.EQ.0) CALL XABORT('USSIT4: MISSING SIGP INFO(1).')
            NPART=ILONG/MI
            CALL LCMGET(KPLIB,'ISM-LIMITS',ISM)
            ALLOCATE(SIGP(NPART,MI),CGAR(NPART))
            CALL LCMGET(KPLIB,'RSE-TABLE',SIGP)
          ELSE IF((ILONG1.EQ.0).AND.(ILONG2.EQ.-1)) THEN
            ! recover a physical probability table
            IF(MI.GT.MAX_R) CALL XABORT('USSIT4: MAX_R OVERFLOW.')
            KPLIB=LCMGIL(JPLIB2,IGRP)
            CALL LCMLEN(KPLIB,'PROB-TABLE',ILONG,ITYLCM)
            IF(ILONG.EQ.0) CALL XABORT('USSIT4: MISSING SIGP INFO(2).')
            NPART=ILONG/MAX_R
            CALL LCMGET(KPLIB,'ISM-LIMITS',ISM)
            ALLOCATE(SIGP(NPART,MI),CGAR(NPART))
            SIGP(:NPART,:MI)=0.0D0
            ALLOCATE(SIGP_R(MAX_R,NPART))
            CALL LCMGET(KPLIB,'PROB-TABLE',SIGP_R)
            DO IPP=1,NPART
              SIGP(IPP,:MI)=SIGP_R(:MI,IPP)
            ENDDO
            DO I=1,NBNRS
              XFLUX(I,:MI)=XFLUX(I,:MI)*SIGP_R(:MI,1)
            ENDDO
            SIGP(1,:MI)=1.0D0
            DEALLOCATE(SIGP_R)
          ELSE
            CALL XABORT('USSIT4: TWO TYPES OF PROBABILITY TABLES.')
          ENDIF
          ! perform weighting of PT or RSE table
          NDEL0=IPPT2(IRES,5)
          DO I=1,NBNRS
            CGAR(1)=0.0D0
            DO IM=1,MI
              CGAR(1)=CGAR(1)+XFLUX(I,IM)*SIGP(1,IM)
            ENDDO
            CGAR(2:NPART)=MATMUL(SIGP(2:NPART,:MI),XFLUX(I,:MI))
            PHGAR(I,IRES,IGRP)=REAL(CGAR(1),4)
            STGAR(I,IRES,IGRP)=REAL(CGAR(2)/CGAR(1),4)
            SFGAR(I,IRES,IGRP)=REAL(CGAR(3)/CGAR(1),4)
            IPP=3
            DO IL=1,NL
              IPP=IPP+1
              SSGAR(I,IRES,IL,IGRP)=REAL(CGAR(IPP)/CGAR(1),4)
            ENDDO
            DO IL=1,NL
              S0GAR(I,IRES,IL,:NGRP,IGRP)=0.0
              DO JGRP=ISM(1,IL),ISM(2,IL)
                IPP=IPP+1
                S0GAR(I,IRES,IL,JGRP,IGRP)=REAL(CGAR(IPP)/CGAR(1),4)
              ENDDO
            ENDDO
            DO IED=1,NED
              IPP=IPP+1
              SAGAR(I,IRES,IED,IGRP)=REAL(CGAR(IPP)/CGAR(1),4)
            ENDDO
            DO IDEL=1,NDEL0
              IPP=IPP+1
              SDGAR(I,IRES,IDEL,IGRP)=REAL(CGAR(IPP)/CGAR(1),4)
            ENDDO
            IF(NDEL0.NE.0) IPP=IPP+NDEL-NDEL0
            IF(IPP.NE.NPART) THEN
              WRITE(TEXT12,'(3A4)') (IPPT2(IRES,J0),J0=2,4)
              WRITE(HSMG,'(26HUSSIT4: FAILURE. ISOTOPE='',A12,
     1        7H'' (IPP=,I6,7H NPART=,I6,6H IGRP=,I6,2H).)') TEXT12,
     2        IPP,NPART,IGRP
              CALL XABORT(HSMG)
            ENDIF
          ENDDO
          DEALLOCATE(CGAR,SIGP)
        ENDDO
        CALL LCMSIX(IPLI0,' ',2)
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION.
*----
      DEALLOCATE(MRANK,ISM,XFLUX)
      RETURN
      END
