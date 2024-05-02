*DECK NSSCO
      SUBROUTINE NSSCO(NX,NY,NZ,NMIX,I,J,K,MAT,XX,YY,ZZ,DIFF,IQFR,QFR,
     1 COEF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the mesh centered finite difference coefficients.
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
* NX      number of X-directed nodes.
* NY      number of Y-directed nodes.
* NZ      number of Z-directed nodes.
* NMIX    number of material mixtures.
* I       X-index of node under consideration.
* J       Y-index of node under consideration.
* K       Z-index of node under consideration.
* MAT     mixture index assigned to each node.
* DIF     diffusion coefficients.
* XX      X-directed mesh spacings.
* YY      Y-directed mesh spacings.
* ZZ      Z-directed mesh spacings.
* IQFR    node-ordered unknown list:
*         =0: neighbour exists;
*         =-1: void/albedo boundary condition;
*         =-2: reflection boundary condition;
*         =-3: ZERO flux boundary condition;
*         =-4: SYME boundary condition (axial symmetry).
* QFR     element-ordered boundary conditions.
*
*Parameters: output
* COEF    mesh centered finite difference coefficients.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NX,NY,NZ,I,J,K,MAT(NX,NY,NZ),IQFR(6)
      REAL DIFF(NMIX),XX(NX,NY,NZ),YY(NX,NY,NZ),ZZ(NX,NY,NZ),QFR(6),
     1 COEF(6)
*----
*  LOCAL VARIABLES
*----
      DHARM(X1,X2,DIF1,DIF2)=2.0*DIF1*DIF2/(X1*DIF2+X2*DIF1)
*
      IBM=MAT(I,J,K)
      DX=XX(I,J,K)
      DY=YY(I,J,K)
      DZ=ZZ(I,J,K)
      KK1=IQFR(1)
      KK2=IQFR(2)
      KK3=IQFR(3)
      KK4=IQFR(4)
      KK5=IQFR(5)
      KK6=IQFR(6)
      COEF(:6)=0.
      ! x- side:
      IF(KK1 == 0) THEN
        COEF(1)=DHARM(DX,XX(I-1,J,K),DIFF(IBM),DIFF(MAT(I-1,J,K)))
      ELSE IF((KK1 > 0).OR.(KK1 == -1)) THEN
        COEF(1)=DHARM(DX,DX,DIFF(IBM),DX*QFR(1)/2.0)
      ELSE IF(KK1 == -2) THEN
        COEF(1)=0.0
      ELSE IF(KK1 == -3) THEN
        COEF(1)=2.0*DHARM(DX,DX,DIFF(IBM),DIFF(IBM))
      ENDIF
      ! x+ side:
      IF(KK2 == 0) THEN
        COEF(2)=DHARM(DX,XX(I+1,J,K),DIFF(IBM),DIFF(MAT(I+1,J,K)))
      ELSE IF((KK2 > 0).OR.(KK2 == -1)) THEN
        COEF(2)=DHARM(DX,DX,DIFF(IBM),DX*QFR(2)/2.0)
      ELSE IF(KK2 == -2) THEN
        COEF(2)=0.0
      ELSE IF(KK2 == -3) THEN
        COEF(2)=2.0*DHARM(DX,DX,DIFF(IBM),DIFF(IBM))
      ELSE IF(KK2 == -4) THEN
        IF(KK1 == -4) CALL XABORT('NSSCO: INCONSISTENT SYME (1).')
        COEF(2)=COEF(1)
      ENDIF
      IF(KK1 == -4) THEN
        IF(KK2 == -4) CALL XABORT('NSSCO: INCONSISTENT SYME (2).')
        COEF(1)=COEF(2)
      ENDIF
      ! y- side:
      IF(KK3 == 0) THEN
        COEF(3)=DHARM(DY,YY(I,J-1,K),DIFF(IBM),DIFF(MAT(I,J-1,K)))
      ELSE IF((KK3 > 0).OR.(KK3 == -1)) THEN
        COEF(3)=DHARM(DY,DY,DIFF(IBM),DY*QFR(3)/2.0)
      ELSE IF(KK3 == -2) THEN
        COEF(3)=0.0
      ELSE IF(KK3 == -3) THEN
        COEF(3)=2.0*DHARM(DY,DY,DIFF(IBM),DIFF(IBM))
      ENDIF
      ! y+ side:
      IF(KK4 == 0) THEN
        COEF(4)=DHARM(DY,YY(I,J+1,K),DIFF(IBM),DIFF(MAT(I,J+1,K)))
      ELSE IF((KK4 > 0).OR.(KK4 == -1)) THEN
        COEF(4)=DHARM(DY,DY,DIFF(IBM),DY*QFR(4)/2.0)
      ELSE IF(KK4 == -2) THEN
        COEF(4)=0.0
      ELSE IF(KK4 == -3) THEN
        COEF(4)=2.0*DHARM(DY,DY,DIFF(IBM),DIFF(IBM))
      ELSE IF(KK4 == -4) THEN
        IF(KK3 == -4) CALL XABORT('NSSCO: INCONSISTENT SYME (3).')
        COEF(4)=COEF(3)
      ENDIF
      IF(KK3 == -4) THEN
        IF(KK4 == -4) CALL XABORT('NSSCO: INCONSISTENT SYME (4).')
        COEF(3)=COEF(4)
      ENDIF
      ! z- side:
      IF(KK5 == 0) THEN
        COEF(5)=DHARM(DZ,ZZ(I,J,K-1),DIFF(IBM),DIFF(MAT(I,J,K-1)))
      ELSE IF((KK5 > 0).OR.(KK5 == -1)) THEN
        COEF(5)=DHARM(DZ,DZ,DIFF(IBM),DZ*QFR(5)/2.0)
      ELSE IF(KK5 == -2) THEN
        COEF(5)=0.0
      ELSE IF(KK5 == -3) THEN
        COEF(5)=2.0*DHARM(DZ,DZ,DIFF(IBM),DIFF(IBM))
      ENDIF
      ! z+ side:
      IF(KK6 == 0) THEN
        COEF(6)=DHARM(DZ,ZZ(I,J,K+1),DIFF(IBM),DIFF(MAT(I,J,K+1)))
      ELSE IF((KK6 > 0).OR.(KK6 == -1)) THEN
        COEF(6)=DHARM(DZ,DZ,DIFF(IBM),DZ*QFR(6)/2.0)
      ELSE IF(KK6 == -2) THEN
        COEF(6)=0.0
      ELSE IF(KK6 == -3) THEN
        COEF(6)=2.0*DHARM(DZ,DZ,DIFF(IBM),DIFF(IBM))
      ELSE IF(KK6 == -4) THEN
        IF(KK5 == -4) CALL XABORT('NSSCO: INCONSISTENT SYME (5).')
        COEF(6)=COEF(5)
      ENDIF
      IF(KK5 == -4) THEN
        IF(KK6 == -4) CALL XABORT('NSSCO: INCONSISTENT SYME (6).')
        COEF(5)=COEF(6)
      ENDIF
      RETURN
      END
