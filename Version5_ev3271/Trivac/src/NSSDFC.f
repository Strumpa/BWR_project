*DECK NSSDFC
      SUBROUTINE NSSDFC(IMPX,IDIM,NX,NY,NZ,NCODE,ICODE,ZCODE,MAT,XXX,
     1 YYY,ZZZ,LL4F,LL4X,LL4Y,LL4Z,VOL,XX,YY,ZZ,IDL,KN,QFR,IQFR,MUX,
     2 MUY,MUZ,IMAX,IMAY,IMAZ,IPY,IPZ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a coarse mesh finite difference (NEM
* type) in a 3-D geometry.
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
* IMPX    print parameter.
* IDIM    number of Cartesian dimensions.
* NX      number of elements along the X axis.
* NY      number of elements along the Y axis.
* NZ      number of elements along the Z axis.
* NCODE   type of boundary condition applied on each side:
*         I=1: X-; I=2: X+; I=3: Y-; I=4: Y+; I=5: Z-; I=6: Z+;
*         NCODE(I)=1: VOID;  NCODE(I)=2: REFL;  NCODE(I)=4: TRAN;
*         NCODE(I)=7: ZERO.
* ICODE   physical albedo index on each side of the domain.
* ZCODE   albedo corresponding to boundary condition 'VOID' on each
*         side (ZCODE(i)=0.0 by default).
* MAT     mixture index assigned to each element.
* XXX     Cartesian coordinates along the X axis.
* YYY     Cartesian coordinates along the Y axis.
* ZZZ     Cartesian coordinates along the Z axis.
* LL4F    total number of averaged flux unknown per energy group.
*
*Parameters: output
* LL4X    total number of X-direccted interface net currents.
* LL4Y    total number of Y-direccted interface net currents.
* LL4Z    total number of Z-direccted interface net currents.
* VOL     volume of each element.
* XX      X-directed mesh spacings.
* YY      Y-directed mesh spacings.
* ZZ      Z-directed mesh spacings.
* IDL     position of averaged fluxes in unknown vector.
* KN      element-ordered interface net current unknown list.
* QFR     element-ordered boundary conditions.
* IQFR    element-ordered physical albedo indices.
* MUX     X-oriented compressed storage mode indices.
* MUY     Y-oriented compressed storage mode indices.
* MUZ     Z-oriented compressed storage mode indices.
* IMAX    X-oriented position of each first non-zero column element.
* IMAY    Y-oriented position of each first non-zero column element.
* IMAZ    Z-oriented position of each first non-zero column element.
* IPY     Y-oriented permutation matrices.
* IPZ     Z-oriented permutation matrices.
*
*-----------------------------------------------------------------------
*
      INTEGER IMPX,IDIM,NX,NY,NZ,NCODE(6),ICODE(6),MAT(NX,NY,NZ),LL4F,
     1 LL4X,LL4Y,LL4Z,IDL(NX,NY,NZ),KN(6,NX,NY,NZ),IQFR(6,NX,NY,NZ),
     2 MUX(LL4F),MUY(LL4F),MUZ(LL4F),IMAX(LL4F),IMAY(LL4F),IMAZ(LL4F),
     3 IPY(LL4F),IPZ(LL4F)
      REAL ZCODE(6),XXX(NX+1),YYY(NY+1),ZZZ(NZ+1),VOL(NX,NY,NZ),
     1 XX(NX,NY,NZ),YY(NX,NY,NZ),ZZ(NX,NY,NZ),QFR(6,NX,NY,NZ)
*----
*  LOCAL VARIABLES
*----
      LOGICAL LL1,LALB
      INTEGER, ALLOCATABLE, DIMENSION(:) :: JPX,JPY,JPZ
*
      ALB(X)=0.5*(1.0-X)/(1.0+X)
*----
*  IDENTIFICATION OF THE NON VIRTUAL NODES
*----
      IF(IMPX.GT.0) WRITE(6,700) NX,NY,NZ
      ALLOCATE(JPX((NX+1)*NY*NZ),JPY((NY+1)*NX*NZ),JPZ((NZ+1)*NX*NY))
      JPX(:)=0
      JPY(:)=0
      JPZ(:)=0
      IND=0
      DO K0=1,NZ
        DO K1=1,NY
          DO K2=1,NX
            IDL(K2,K1,K0)=0
            KN(:6,K2,K1,K0)=0
            IF(MAT(K2,K1,K0).EQ.0) CYCLE
            IND=IND+1
            IDL(K2,K1,K0)=IND
            KN(1,K2,K1,K0)=K2    +(NX+1)*(K1-1)+(NX+1)*NY*(K0-1)
            KN(2,K2,K1,K0)=(K2+1)+(NX+1)*(K1-1)+(NX+1)*NY*(K0-1)
            KN(3,K2,K1,K0)=K1    +(NY+1)*(K0-1)+(NY+1)*NZ*(K2-1)
            KN(4,K2,K1,K0)=(K1+1)+(NY+1)*(K0-1)+(NY+1)*NZ*(K2-1)
            KN(5,K2,K1,K0)=K0    +(NZ+1)*(K2-1)+(NZ+1)*NX*(K1-1)
            KN(6,K2,K1,K0)=(K0+1)+(NZ+1)*(K2-1)+(NZ+1)*NX*(K1-1)
            JPX(KN(1:2,K2,K1,K0))=1
            JPY(KN(3:4,K2,K1,K0))=1
            JPZ(KN(5:6,K2,K1,K0))=1
          ENDDO
        ENDDO
      ENDDO
      IF(IND.NE.LL4F) CALL XABORT('NSSDFC: WRONG VALUE OF LL4F.')
      LL4X=0
      DO I=1,(NX+1)*NY*NZ
        IF(JPX(I).EQ.1) THEN
          LL4X=LL4X+1
          JPX(I)=LL4X
        ENDIF
      ENDDO
      LL4Y=0
      DO I=1,(NY+1)*NX*NZ
        IF(JPY(I).EQ.1) THEN
          LL4Y=LL4Y+1
          JPY(I)=LL4Y
        ENDIF
      ENDDO
      LL4Z=0
      DO I=1,(NZ+1)*NX*NY
        IF(JPZ(I).EQ.1) THEN
          LL4Z=LL4Z+1
          JPZ(I)=LL4Z
        ENDIF
      ENDDO
      DO K0=1,NZ
        DO K1=1,NY
          DO K2=1,NX
            IF(MAT(K2,K1,K0).EQ.0) CYCLE
            KN(1:2,K2,K1,K0)=JPX(KN(1:2,K2,K1,K0))
            KN(3:4,K2,K1,K0)=JPY(KN(3:4,K2,K1,K0))
            KN(5:6,K2,K1,K0)=JPZ(KN(5:6,K2,K1,K0))
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(JPZ,JPY,JPX)
*----
*  IDENTIFICATION OF THE GEOMETRY. MAIN LOOP OVER THE NODES
*----
      QFR(:6,:NX,:NY,:NZ)=0.0
      IQFR(:6,:NX,:NY,:NZ)=-99
      DO K0=1,NZ
        DO K1=1,NY
          DO K2=1,NX
            XX(K2,K1,K0)=0.0
            YY(K2,K1,K0)=0.0
            ZZ(K2,K1,K0)=0.0
            VOL(K2,K1,K0)=0.0
            IF(MAT(K2,K1,K0).LE.0) CYCLE
            XX(K2,K1,K0)=XXX(K2+1)-XXX(K2)
            YY(K2,K1,K0)=YYY(K1+1)-YYY(K1)
            ZZ(K2,K1,K0)=ZZZ(K0+1)-ZZZ(K0)
*----
*  VOID, REFL OR ZERO BOUNDARY CONTITION
*----
            IQFR(:2,K2,K1,K0)=0
            IF(K2.EQ.1) THEN
              LL1=.TRUE.
            ELSE
              LL1=(MAT(K2-1,K1,K0).EQ.0)
            ENDIF
            IF(LL1) THEN
              LALB=(NCODE(1).EQ.1).OR.(NCODE(1).EQ.6)
              IF(LALB.AND.(ICODE(1).EQ.0)) THEN
                QFR(1,K2,K1,K0)=ALB(ZCODE(1))
                IQFR(1,K2,K1,K0)=-1
              ELSE IF(LALB) THEN
                QFR(1,K2,K1,K0)=1.0
                IQFR(1,K2,K1,K0)=ICODE(1)
              ELSE IF(NCODE(1).EQ.2) THEN
                IQFR(1,K2,K1,K0)=-2
              ELSE IF(NCODE(1).EQ.7) THEN
                IQFR(1,K2,K1,K0)=-3
              ELSE IF(NCODE(1).EQ.5) THEN
                CALL XABORT('NSSDFC: SYME NOT IMPLEMENTED(1).')
              ENDIF
            ENDIF
*
            IF(K2.EQ.NX) THEN
              LL1=.TRUE.
            ELSE
              LL1=(MAT(K2+1,K1,K0).EQ.0)
            ENDIF
            IF(LL1) THEN
              LALB=(NCODE(2).EQ.1).OR.(NCODE(2).EQ.6)
              IF(LALB.AND.(ICODE(2).EQ.0)) THEN
                QFR(2,K2,K1,K0)=ALB(ZCODE(2))
                IQFR(2,K2,K1,K0)=-1
              ELSE IF(LALB) THEN
                QFR(2,K2,K1,K0)=1.0
                IQFR(2,K2,K1,K0)=ICODE(2)
              ELSE IF(NCODE(2).EQ.2) THEN
                IQFR(2,K2,K1,K0)=-2
              ELSE IF(NCODE(2).EQ.7) THEN
                IQFR(2,K2,K1,K0)=-3
              ELSE IF(NCODE(1).EQ.5) THEN
                CALL XABORT('NSSDFC: SYME NOT IMPLEMENTED(2).')
              ENDIF
            ENDIF
*
            IF(IDIM == 1) GO TO 100
            IQFR(3:4,K2,K1,K0)=0
            IF(K1.EQ.1) THEN
              LL1=.TRUE.
            ELSE
              LL1=(MAT(K2,K1-1,K0).EQ.0)
            ENDIF
            IF(LL1) THEN
              LALB=(NCODE(3).EQ.1).OR.(NCODE(3).EQ.6)
              IF(LALB.AND.(ICODE(3).EQ.0)) THEN
                QFR(3,K2,K1,K0)=ALB(ZCODE(3))
                IQFR(3,K2,K1,K0)=-1
              ELSE IF(LALB) THEN
                QFR(3,K2,K1,K0)=1.0
                IQFR(3,K2,K1,K0)=ICODE(3)
              ELSE IF(NCODE(3).EQ.2) THEN
                IQFR(3,K2,K1,K0)=-2
              ELSE IF(NCODE(3).EQ.7) THEN
                IQFR(3,K2,K1,K0)=-3
              ELSE IF(NCODE(1).EQ.5) THEN
                CALL XABORT('NSSDFC: SYME NOT IMPLEMENTED(3).')
              ENDIF
            ENDIF
*
            IF(K1.EQ.NY) THEN
              LL1=.TRUE.
            ELSE
              LL1=(MAT(K2,K1+1,K0).EQ.0)
            ENDIF
            IF(LL1) THEN
              LALB=(NCODE(4).EQ.1).OR.(NCODE(4).EQ.6)
              IF(LALB.AND.(ICODE(4).EQ.0)) THEN
                QFR(4,K2,K1,K0)=ALB(ZCODE(4))
                IQFR(4,K2,K1,K0)=-1
              ELSE IF(LALB) THEN
                QFR(4,K2,K1,K0)=1.0
                IQFR(4,K2,K1,K0)=ICODE(4)
              ELSE IF(NCODE(4).EQ.2) THEN
                IQFR(4,K2,K1,K0)=-2
              ELSE IF(NCODE(4).EQ.7) THEN
                IQFR(4,K2,K1,K0)=-3
              ELSE IF(NCODE(1).EQ.5) THEN
                CALL XABORT('NSSDFC: SYME NOT IMPLEMENTED(4).')
              ENDIF
            ENDIF
*
            IF(IDIM == 2) GO TO 100
            IQFR(5:6,K2,K1,K0)=0
            IF(K0.EQ.1) THEN
              LL1=.TRUE.
            ELSE
              LL1=(MAT(K2,K1,K0-1).EQ.0)
            ENDIF
            IF(LL1) THEN
              LALB=(NCODE(5).EQ.1).OR.(NCODE(5).EQ.6)
              IF(LALB.AND.(ICODE(5).EQ.0)) THEN
                QFR(5,K2,K1,K0)=ALB(ZCODE(5))
                IQFR(5,K2,K1,K0)=-1
              ELSE IF(LALB) THEN
                QFR(5,K2,K1,K0)=1.0
                IQFR(5,K2,K1,K0)=ICODE(5)
              ELSE IF(NCODE(5).EQ.2) THEN
                IQFR(5,K2,K1,K0)=-2
              ELSE IF(NCODE(5).EQ.7) THEN
                IQFR(5,K2,K1,K0)=-3
              ELSE IF(NCODE(1).EQ.5) THEN
                CALL XABORT('NSSDFC: SYME NOT IMPLEMENTED(5).')
              ENDIF
            ENDIF
*
            IF(K0.EQ.NZ) THEN
              LL1=.TRUE.
            ELSE
              LL1=(MAT(K2,K1,K0+1).EQ.0)
            ENDIF
            IF(LL1) THEN
              LALB=(NCODE(6).EQ.1).OR.(NCODE(6).EQ.6)
              IF(LALB.AND.(ICODE(6).EQ.0)) THEN
                QFR(6,K2,K1,K0)=ALB(ZCODE(6))
                IQFR(6,K2,K1,K0)=-1
              ELSE IF(LALB) THEN
                QFR(6,K2,K1,K0)=1.0
                IQFR(6,K2,K1,K0)=ICODE(6)
              ELSE IF(NCODE(6).EQ.2) THEN
                IQFR(6,K2,K1,K0)=-2
              ELSE IF(NCODE(6).EQ.7) THEN
                IQFR(6,K2,K1,K0)=-3
              ELSE IF(NCODE(1).EQ.5) THEN
                CALL XABORT('NSSDFC: SYME NOT IMPLEMENTED(6).')
              ENDIF
            ENDIF
*----
*  TRAN BOUNDARY CONDITION
*----
  100       IF((K2.EQ.1).AND.(NCODE(1).EQ.4)) THEN
              KN(1,K2,K1,K0)=KN(2,NX,K1,K0)
            ENDIF
            IF((K2.EQ.NX).AND.(NCODE(2).EQ.4)) THEN
              KN(2,K2,K1,K0)=KN(1,1,K1,K0)
            ENDIF
            IF((K1.EQ.1).AND.(NCODE(3).EQ.4)) THEN
              KN(3,K2,K1,K0)=KN(2,K2,NY,K0)
            ENDIF
            IF((K1.EQ.NY).AND.(NCODE(4).EQ.4)) THEN
              KN(4,K2,K1,K0)=KN(1,K2,1,K0)
            ENDIF
            IF((K0.EQ.1).AND.(NCODE(5).EQ.4)) THEN
              KN(5,K2,K1,K0)=KN(6,K2,K1,NZ)
            ENDIF
            IF((K0.EQ.NZ).AND.(NCODE(6).EQ.4)) THEN
              KN(6,K2,K1,K0)=KN(5,K2,K1,1)
            ENDIF
*
            VOL(K2,K1,K0)=XX(K2,K1,K0)*YY(K2,K1,K0)*ZZ(K2,K1,K0)
          ENDDO
        ENDDO
      ENDDO
* END OF THE MAIN LOOP OVER NODES.
*
      IF(IMPX.GE.2) THEN
         WRITE(6,720) VOL(:NX,:NY,:NZ)
         WRITE(6,750)
         DO K0=1,NZ
           DO K1=1,NY
             DO K2=1,NX
               IF(MAT(K2,K1,K0).LE.0) CYCLE
               KEL=(K0-1)*NX*NY+(K1-1)*NX+K2
               WRITE (6,760) KEL,(KN(I,K2,K1,K0),I=1,6),
     1         (QFR(I,K2,K1,K0),I=1,6),(IQFR(I,K2,K1,K0),I=1,6)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
*----
*  COMPUTE THE PERMUTATION VECTORS IPY AND IPZ
*----
      IF(IDIM.GE.2) THEN
         INX1=0
         DO K2=1,NX
           DO K0=1,NZ
             DO K1=1,NY
               INX2=IDL(K2,K1,K0)
               IF(INX2.LE.0) CYCLE
               INX1=INX1+1
               IPY(INX2)=INX1
             ENDDO
           ENDDO
         ENDDO
         IF(INX1.NE.IND) CALL XABORT('NSSDFC: FAILURE OF THE RENUMBERI'
     1   //'NG ALGORITHM(1)')
         IF(IDIM.EQ.3) THEN
            INX1=0
            DO K1=1,NY
              DO K2=1,NX
                DO K0=1,NZ
                  INX2=IDL(K2,K1,K0)
                  IF(INX2.LE.0) CYCLE
                  INX1=INX1+1
                  IPZ(INX2)=INX1
                ENDDO
              ENDDO
            ENDDO
            IF(INX1.NE.IND) CALL XABORT('NSSDFC: FAILURE OF THE RENUMB'
     1      //'ERING ALGORITHM(2)')
         ENDIF
      ENDIF
*----
*  COMPUTE VECTOR MUX
*----
      MUX(:LL4F)=1
      DO K0=1,NZ
        DO K1=1,NY
*         X- SIDE:
          DO K2=2,NX
            KEL=IDL(K2,K1,K0)
            IF(KEL.EQ.0) CYCLE
            KK1=IDL(K2-1,K1,K0)
            IF(KK1.GT.0) MUX(KEL)=MAX0(MUX(KEL),KEL-KK1+1)
          ENDDO
*         X+ SIDE:
          DO K2=1,NX-1
            KEL=IDL(K2,K1,K0)
            IF(KEL.EQ.0) CYCLE
            KK2=IDL(K2+1,K1,K0)
            IF(KK2.GT.0) MUX(KEL)=MAX0(MUX(KEL),KEL-KK2+1)
          ENDDO
        ENDDO
      ENDDO
*----
*  COMPUTE VECTOR MUY
*----
      IF(IDIM.GE.2) THEN
        MUY(:LL4F)=1
        DO K2=1,NX
          DO K0=1,NZ
*           Y- SIDE:
            DO K1=2,NY
              KEL=IDL(K2,K1,K0)
              IF(KEL.EQ.0) CYCLE
              INY1=IPY(KEL)
              KK3=IDL(K2,K1-1,K0)
              IF(KK3.GT.0) MUY(INY1)=MAX0(MUY(INY1),INY1-IPY(KK3)+1)
            ENDDO
*           Y- SIDE:
            DO K1=1,NY-1
              KEL=IDL(K2,K1,K0)
              IF(KEL.EQ.0) CYCLE
              INY1=IPY(KEL)
              KK4=IDL(K2,K1+1,K0)
              IF(KK4.GT.0) MUY(INY1)=MAX0(MUY(INY1),INY1-IPY(KK4)+1)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        MUY(:LL4F)=0
      ENDIF
*----
*  COMPUTE VECTOR MUZ
*----
      IF(IDIM.EQ.3) THEN
        MUZ(:LL4F)=1
        DO K1=1,NY
          DO K2=1,NX
*           Z- SIDE:
            DO K0=2,NZ
              KEL=IDL(K2,K1,K0)
              IF(KEL.EQ.0) CYCLE
              INZ1=IPZ(KEL)
              KK5=IDL(K2,K1,K0-1)
              IF(KK5.GT.0) MUZ(INZ1)=MAX0(MUZ(INZ1),INZ1-IPZ(KK5)+1)
            ENDDO
*           Z+ SIDE:
            DO K0=1,NZ-1
              KEL=IDL(K2,K1,K0)
              IF(KEL.EQ.0) CYCLE
              INZ1=IPZ(KEL)
              KK6=IDL(K2,K1,K0+1)
              IF(KK6.GT.0) MUZ(INZ1)=MAX0(MUZ(INZ1),INZ1-IPZ(KK6)+1)
            ENDDO
          ENDDO
        ENDDO
      ELSE
        MUZ(:LL4F)=0
      ENDIF
*
      MUXMAX=0
      MUYMAX=0
      MUZMAX=0
      IIMAXX=0
      IIMAXY=0
      IIMAXZ=0
      DO I=1,LL4F
        MUXMAX=MAX(MUXMAX,MUX(I))
        MUYMAX=MAX(MUYMAX,MUY(I))
        MUZMAX=MAX(MUZMAX,MUZ(I))
        IBAND=MUX(I)
        IIMAXX=IIMAXX+IBAND
        MUX(I)=IIMAXX
        IIMAXX=IIMAXX+IBAND-1
        IMAX(I)=IIMAXX
        IBAND=MUY(I)
        IIMAXY=IIMAXY+IBAND
        MUY(I)=IIMAXY
        IIMAXY=IIMAXY+IBAND-1
        IMAY(I)=IIMAXY
        IBAND=MUZ(I)
        IIMAXZ=IIMAXZ+IBAND
        MUZ(I)=IIMAXZ
        IIMAXZ=IIMAXZ+IBAND-1
        IMAZ(I)=IIMAXZ
      ENDDO
      IF(IMPX.GT.0) WRITE (6,770) MUXMAX,MUYMAX,MUZMAX
      RETURN
*
  700 FORMAT(/46H NSSDFC: COARSE MESH FINITE DIFFERENCE METHOD.//3H NU,
     1 28HMBER OF NODES ALONG X AXIS =,I3/17X,14HALONG Y AXIS =,I3/
     2 17X,14HALONG Z AXIS =,I3)
  720 FORMAT(/17H VOLUMES PER NODE/(1X,1P,10E13.4))
  750 FORMAT(/22H NUMBERING OF UNKNOWNS/1X,21(1H-)//4X,4HNODE,5X,3HINT,
     1 26HERFACE NET CURRENT INDICES,28X,23HVOID BOUNDARY CONDITION)
  760 FORMAT(1X,I6,7X,6I8,6X,6F9.2/68X,6I9)
  770 FORMAT(/41H NSSDFC: MAXIMUM BANDWIDTH ALONG X AXIS =,I5/
     1 27X,14HALONG Y AXIS =,I5/27X,14HALONG Z AXIS =,I5)
      END
