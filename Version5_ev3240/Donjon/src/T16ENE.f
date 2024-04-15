*DECK T16ENE
      SUBROUTINE T16ENE(IPRINT,MXGRP ,NG    ,NGCOND,NGMTR ,NGREAC,
     >                  NGCCPO,IFGCND,IFGMTR,IFGEDI,ENECPO,ENET16,
     >                  VELMTR)
*
*----
*
*Purpose:
*  Generate and analyse energy structure.
*
*Author(s): 
* G. Marleau
*
*Parameters: input
* IPRINT  print level.
* MXGRP   maximum number of groups.
* NG      number of groups in library.
* NGCOND  number of condensed groups.
* NGMTR   numbre of main transport groups.
* NGREAC  numbre of edit groups.
*
*Parameters: input/output
* NGCCPO  numbre of edit groups.
* IFGCND  reference/exit condensation few groups.
* IFGMTR  reference/exit main transport few groups.
* IFGEDI  reference/exit edit few groups.
* ENECPO  final energy group structure for CPO.
*
*Parameters: output
* ENET16  energy group structure for tape16.
* VELMTR  velocity for main transport.
*
*----
*
      IMPLICIT         NONE
      INTEGER          IPRINT,MXGRP,NG,NGCOND,NGMTR,NGREAC,NGCCPO
      INTEGER          IFGCND(MXGRP),IFGMTR(MXGRP),
     >                 IFGEDI(MXGRP)
      REAL             ENECPO(MXGRP+1),ENET16(MXGRP+1),
     >                 VELMTR(MXGRP)
*----
*  LOCAL VARIABLES
*  FOR AVERAGED NEUTRON VELOCITY
*  V=SQRT(2*ENER/M)=SQRT(2/M)*SQRT(ENER)
*  SQFMAS=SQRT(2/M) IN CM/S/SQRT(EV) FOR V IN CM/S AND E IN EV
*        =SQRT(2*1.602189E-19(J/EV)* 1.0E4(CM2/M2) /1.67495E-27 (KG))
*        =1383155.30602 CM/S/SQRT(EV)
*----
      INTEGER          IOUT,MGELIB,MGWLIB
      CHARACTER        NAMSBR*6
      REAL             SQFMAS,PRECIS
      PARAMETER       (IOUT=6,MGELIB=89,MGWLIB=69,
     >                 NAMSBR='T16ENE',SQFMAS=1383155.30602,
     >                 PRECIS=1.0E-5)
      INTEGER          IGR,IGC,IGD,IGF
      REAL             EAVG
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IFGCPO
      REAL, ALLOCATABLE, DIMENSION(:) :: VELEDI,VELT16
*----
*  DATA
*----
      REAL             ENEELB(MGELIB+1),ENEWLB(MGWLIB+1)
      SAVE             ENEELB,ENEWLB
      DATA             ENEELB
     >/1.0000E+7,7.7880E+6,6.0653E+6,4.7237E+6,3.6788E+6,2.8650E+6,
     > 2.2313E+6,1.7377E+6,1.3534E+6,1.0540E+6,8.2085E+5,6.3928E+5,
     > 4.9787E+5,3.8774E+5,3.0197E+5,2.3518E+5,1.8316E+5,1.4264E+5,
     > 1.1109E+5,8.6517E+4,6.7379E+4,4.0868E+4,2.4788E+4,1.5034E+4,
     > 9.1188E+3,5.5308E+3,3.3546E+3,2.0347E+3,1.2341E+3,7.4852E+2,
     > 4.5400E+2,2.7536E+2,1.6702E+2,1.3007E+2,1.0130E+2,7.8893E+1,
     > 6.1442E+1,4.7851E+1,3.7267E+1,2.9023E+1,2.2603E+1,1.7603E+1,
     > 1.3710E+1,1.0677E+1,8.3153E+0,6.4760E+0,5.0435E+0,4.0000E+0,
     > 3.3000E+0,2.6000E+0,2.1000E+0,1.5000E+0,1.3000E+0,1.1500E+0,
     > 1.1230E+0,1.0970E+0,1.0710E+0,1.0450E+0,1.0200E+0,9.9600E-1,
     > 9.7200E-1,9.5000E-1,9.1000E-1,8.5000E-1,7.8000E-1,6.2500E-1,
     > 5.0000E-1,4.0000E-1,3.5000E-1,3.2000E-1,3.0000E-1,2.8000E-1,
     > 2.5000E-1,2.2000E-1,1.8000E-1,1.4000E-1,1.0000E-1,8.0000E-2,
     > 6.7000E-2,5.8000E-2,5.0000E-2,4.2000E-2,3.5000E-2,3.0000E-2,
     > 2.5000E-2,2.0000E-2,1.5000E-2,1.0000E-2,5.0000E-3,2.0000E-4/
      DATA             ENEWLB
     >/1.0000E+7,6.0655E+6,3.6790E+6,2.2310E+6,1.3530E+6,8.2100E+5,
     > 5.0000E+5,3.0250E+5,1.8300E+5,1.1100E+5,6.7340E+4,4.0850E+4,
     > 2.4780E+4,1.5030E+4,9.1180E+3,5.5300E+3,3.5191E+3,2.2394E+3,
     > 1.4251E+3,9.0690E+2,3.6726E+2,1.4873E+2,7.5501E+1,4.8052E+1,
     > 2.7700E+1,1.5968E+1,9.8770E+0,4.0000E+0,3.3000E+0,2.6000E+0,
     > 2.1000E+0,1.5000E+0,1.3000E+0,1.1500E+0,1.1230E+0,1.0970E+0,
     > 1.0710E+0,1.0450E+0,1.0200E+0,9.9600E-1,9.7200E-1,9.5000E-1,
     > 9.1000E-1,8.5000E-1,7.8000E-1,6.2500E-1,5.0000E-1,4.0000E-1,
     > 3.5000E-1,3.2000E-1,3.0000E-1,2.8000E-1,2.5000E-1,2.2000E-1,
     > 1.8000E-1,1.4000E-1,1.0000E-1,8.0000E-2,6.7000E-2,5.8000E-2,
     > 5.0000E-2,4.2000E-2,3.5000E-2,3.0000E-2,2.5000E-2,2.0000E-2,
     > 1.5000E-2,1.0000E-2,5.0000E-3,1.0000E-5/
*----
*  STORE ORIGINAL GROUP STRUCTURE IN ENET16
*----
      ALLOCATE(IFGCPO(MXGRP),VELEDI(MXGRP),VELT16(MXGRP))
      IF(NG .EQ. MGELIB) THEN
        ENET16(1)=ENEELB(1)
        DO IGR=2,MGELIB+1
          ENET16(IGR)=ENEELB(IGR)
          EAVG=SQRT(ENET16(IGR)*ENET16(IGR-1))
          VELT16(IGR-1)=SQFMAS*SQRT(EAVG)
        ENDDO
      ELSE IF(NG .EQ. MGWLIB) THEN
        ENET16(1)=ENEWLB(1)
        DO IGR=2,MGWLIB+1
          ENET16(IGR)=ENEWLB(IGR)
          EAVG=SQRT(ENET16(IGR)*ENET16(IGR-1))
          VELT16(IGR-1)=SQFMAS*SQRT(EAVG)
        ENDDO
      ELSE
        CALL XABORT(NAMSBR//
     >  ': INVALID TAPE16 GROUP STRUCTURE')
      ENDIF
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6010) NG,NGMTR,NGREAC,NGCOND,NGCCPO
        WRITE(IOUT,6030) (ENET16(IGC),IGC=1,NG+1)
        WRITE(IOUT,6040) (VELT16(IGC),IGC=1,NG)
      ENDIF
*----
*  COMPUTE AVERAGED NEUTRON GROUP VELOCITY
*  AVERAGED NEUTRON ENERGY ENER=SQRT(E(G+1)*E(G))
*  V=SQRT(2*ENER/M)=SQRT(2/M)*SQRT(ENER)
*   =SQFMAS*SQRT(ENER)
*----
      IF(NGMTR .GT. 0) THEN
        IGF=1
        DO IGR=1,NGMTR
          IGD=IGF
          IGF=IFGMTR(IGR)+1
          EAVG=SQRT(ENET16(IGD)*ENET16(IGF))
          VELMTR(IGR)=SQFMAS*SQRT(EAVG)
        ENDDO
      ENDIF
      IF(IPRINT .GE. 10 .AND. NGMTR .GT. 0) THEN
        WRITE(IOUT,6020) (IFGMTR(IGC),IGC=1,NGMTR)
        WRITE(IOUT,6031) ENET16(1),
     >  (ENET16(IFGMTR(IGC)+1),IGC=1,NGMTR)
        WRITE(IOUT,6041) (VELMTR(IGC),IGC=1,NGMTR)
      ENDIF
      IF(NGREAC .GT. 0) THEN
        IGF=1
        DO IGR=1,NGREAC
          IGD=IGF
          IGF=IFGEDI(IGR)+1
          EAVG=SQRT(ENET16(IGD)*ENET16(IGF))
          VELEDI(IGR)=SQFMAS*SQRT(EAVG)
        ENDDO
*----
*    TEST IF CONDENSATION STRUCTURE PROVIDED BY IFGEDI
*    COMPATIBLE WITH IFGMTR
*----
        IF(NGMTR .GT. 0) THEN
          DO IGC=1,NGREAC
            DO IGR=IGC,NGMTR
              IF(IFGEDI(IGC) .EQ. IFGMTR(IGR)) THEN
                GO TO 105
              ENDIF
            ENDDO
            CALL XABORT(NAMSBR//
     >      ': IFGEDI AND IFGMTR NOT COMPATIBLE')
 105        CONTINUE
          ENDDO
        ENDIF
      ENDIF
      IF(IPRINT .GE. 10 .AND. NGREAC .GT. 0) THEN
        WRITE(IOUT,6021) (IFGEDI(IGC),IGC=1,NGREAC)
        WRITE(IOUT,6032) ENET16(1),
     >  (ENET16(IFGEDI(IGC)+1),IGC=1,NGREAC)
        WRITE(IOUT,6042) (VELEDI(IGC),IGC=1,NGREAC)
      ENDIF
*----
*  IF NGCCPO > 0 FIND IFGCPO FROM ENECPO
*----
      IF(NGCCPO .GT. 0) THEN
        IF(ABS(ENECPO(1)-ENET16(1)) .GT. PRECIS ) CALL XABORT(NAMSBR//
     >  ': ENECPO(1) SHOULD BE IDENTICAL TO ENET16(1)')
        DO IGC=2,NGCCPO+1
          DO IGR=IGC,NG+1
            IF(ABS(ENECPO(IGC)-ENET16(IGR)) .LT. PRECIS ) THEN
              IFGCPO(IGC-1)=IGR-1
              GO TO 115
            ENDIF
          ENDDO
 115      CONTINUE
        ENDDO
        IF(NGCOND .GT. 0) THEN
*----
*  IF NGCOND > 0
*  IFGCPO AND IFGCND NUST BE IDENTICAL
*----
          IF(NGCCPO .NE. NGCOND) CALL XABORT(NAMSBR//
     >    ': NGCCPO AND NGCOND MUST BE IDENTICAL')
          DO IGC=1,NGCCPO
            IF(IFGCPO(IGC) .NE. IFGCND(IGC))
     >      CALL XABORT(NAMSBR//
     >      ': IFGCPO AND IFGCND MUST BE IDENTICAL')
          ENDDO
        ENDIF
      ELSE
*----
*  IF NGCCPO =0
*----
        IF(NGCOND .GT. 0) THEN
*----
*  IF NGCOND > 0
*  IFGCPO = IFGCND
*----
          NGCCPO=NGCOND
          ENECPO(1)=ENET16(1)
          DO IGC=1,NGCCPO
            IFGCPO(IGC)=IFGCND(IGC)
            ENECPO(IGC+1)=ENET16(IFGCPO(IGC)+1)
          ENDDO
        ELSE IF(NGREAC .GT. 0) THEN
*----
*  IF NGCOND = 0
*  AND NGREAC > 0
*  IFGCPO = IFGEDI
*----
          NGCCPO=NGREAC
          ENECPO(1)=ENET16(1)
          DO IGC=1,NGCCPO
            IFGCPO(IGC)=IFGEDI(IGC)
            ENECPO(IGC+1)=ENET16(IFGCPO(IGC)+1)
          ENDDO
        ELSE
*----
*  IF NGCOND = 0
*  AND NGREAC = 0
*  IFGCPO = IFGMTR
*----
          NGCCPO=NGMTR
          ENECPO(1)=ENET16(1)
          DO IGC=1,NGCCPO
            IFGCPO(IGC)=IFGMTR(IGC)
            ENECPO(IGC+1)=ENET16(IFGCPO(IGC)+1)   
          ENDDO
        ENDIF
      ENDIF
      IF(NGREAC .GT. 0) THEN
*----
*  IF NGREAC > 0
*    TEST IF CONDENSATION STRUCTURE PROVIDED BY IFGEDI
*    COMPATIBLE WITH IFGCPO AND IFGMTR
*  ENDIF
*----
        
        DO IGC=1,NGCCPO
          DO IGR=IGC,NGREAC
            IF(IFGCPO(IGC) .EQ. IFGEDI(IGR)) THEN
              IFGEDI(IGC)=IGR
              GO TO 135
            ENDIF
          ENDDO
          CALL XABORT(NAMSBR//
     >    ': IFGCPO AND IFGEDI NOT COMPATIBLE')
 135      CONTINUE
        ENDDO
      ENDIF
*----
*  NGMTR > 0
*    TEST IF CONDENSATION STRUCTURE PROVIDED BY IFGMTR
*    COMPATIBLE WITH IFGCPO
*  ENDIF
*----
      DO IGC=1,NGCCPO
        DO IGR=IGC,NGMTR
          IF(IFGCPO(IGC) .EQ. IFGMTR(IGR)) THEN
            IFGMTR(IGC)=IGR
            GO TO 155
          ENDIF
        ENDDO
        CALL XABORT(NAMSBR//
     >  ': IFGCPO AND IFGMTR NOT COMPATIBLE')
 155    CONTINUE
      ENDDO
      IF(IPRINT .GE.10) THEN
        IF(NGCOND .GT. 0)
     >    WRITE(IOUT,6022) (IFGCND(IGC),IGC=1,NGCOND)
        IF(NGCCPO .GT. 0) THEN
          WRITE(IOUT,6023) (IFGCPO(IGC),IGC=1,NGCCPO)
          WRITE(IOUT,6033) (ENECPO(IGC),IGC=1,NGCCPO+1)
        ENDIF
        WRITE(IOUT,6001)
      ENDIF
      DEALLOCATE(VELT16,VELEDI,IFGCPO)
      RETURN
*----
*  PRINT FORMAT
*----
 6000 FORMAT(1X,5('*'),' OUTPUT FROM ',A6,1X,5('*'))
 6001 FORMAT(1X,30('*'))
 6010 FORMAT(6X,'NUMBER OF LIBRARY GROUPS          = ',I10/
     >       6X,'NUMBER OF MAIN TRANSPORT GROUPS   = ',I10/
     >       6X,'NUMBER OF EDITING GROUPS          = ',I10/
     >       6X,'NUMBER OF CONDENSATION GROUPS     = ',I10/
     >       6X,'NUMBER OF CPO GROUPS              = ',I10)
 6020 FORMAT(6X,'MAIN TRANSPORT FEW GROUPS IDENTIFIER '/
     >10(2X,I6))
 6021 FORMAT(6X,'EDIT FEW GROUPS IDENTIFIER '/
     >10(2X,I6))
 6022 FORMAT(6X,'CONDENSATION FEW GROUPS IDENTIFIER '/
     >10(2X,I6))
 6023 FORMAT(6X,'CPO FEW GROUPS IDENTIFIER '/
     >10(2X,I6))
 6030 FORMAT(6X,'INITIAL ENERGY STRUCTURE (EV)'/
     >1P,10(2X,E10.3))
 6031 FORMAT(6X,'ENERGY STRUCTURE IN MAIN GROUPS (EV)'/
     >1P,10(2X,E10.3))
 6032 FORMAT(6X,'ENERGY STRUCTURE IN EDIT GROUPS (EV)'/
     >1P,10(2X,E10.3))
 6033 FORMAT(6X,'FINAL ENERGY STRUCTURE (EV)'/
     >1P,10(2X,E10.3))
 6040 FORMAT(6X,'AVERAGE VELOCITY IN INITIAL GROUPS (CM/S)'/
     >1P,10(2X,E10.3))
 6041 FORMAT(6X,'AVERAGE VELOCITY IN MAIN GROUPS (CM/S)'/
     >1P,10(2X,E10.3))
 6042 FORMAT(6X,'AVERAGE VELOCITY IN EDIT GROUPS (CM/S)'/
     >1P,10(2X,E10.3))
      END