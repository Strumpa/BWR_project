*DECK VALU5
      SUBROUTINE VALU5 (KPMAC,NX,NY,NZ,LL4F,LL4X,LL4Y,NUN,NMIX,X,Y,Z,
     1 XXX,YYY,ZZZ,EVT,ISS,KFLX,KN,IXLG,IYLG,IZLG,ICORN,AXYZ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Interpolation of the flux distribution for nodal method in 3D.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* KPMAC   group directory in the macrolib.
* NX      number of elements along the X axis.
* NY      number of elements along the Y axis.
* NY      number of elements along the Z axis.
* LL4F    number of averaged flux unknowns.
* LL4X    number of X-directed net currents.
* LL4Y    number of Y-directed net currents.
* NUN     dimension of unknown array EVT.
* NMIX    number of mixtures.
* X       Cartesian coordinates along the X axis where the flux is
*         interpolated.
* Y       Cartesian coordinates along the Y axis where the flux is
*         interpolated.
* Z       Cartesian coordinates along the Z axis where the flux is
*         interpolated.
* XXX     Cartesian coordinates along the X axis.
* YYY     Cartesian coordinates along the Y axis.
* ZZZ     Cartesian coordinates along the Z axis.
* EVT     reconstruction coefficients of the flux.
* ISS     mixture index assigned to each element.
* KFLX    correspondence between local and global numbering.
* KN      element-ordered interface net current unknown list.
* IXLG    number of interpolated points according to X.
* IYLG    number of interpolated points according to Y.
* IZLG    number of interpolated points according to Z.
* ICORN   flag to activate corner flux correction (0/1: ON/OFF).
*                                                                      
*Parameters: output
* AXYZ    interpolated fluxes.
*                                                                      
*----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) KPMAC
      INTEGER NX,NY,NZ,LL4F,LL4X,LL4Y,NUN,NMIX,ISS(NX*NY*NZ),
     1 KFLX(NX*NY*NZ),KN(6,NX,NY,NZ),IXLG,IYLG,IZLG,ICORN
      REAL X(IXLG),Y(IYLG),Z(IZLG),XXX(NX+1),YYY(NY+1),ZZZ(NZ+1),
     1 EVT(NUN),AXYZ(IXLG,IYLG,IZLG)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION  WORK1(4,5),FC2(8)
      DOUBLE PRECISION GAR,COEFX,COEFY,COEFZ,U,V,W,P2U,P2V,P2W
      LOGICAL LOGC1,LOGC2,LOGC3
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: DIFF
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: FCORN,DELC
*----
*  RECOVER DIFFUSION COEFFICIENTS
*----
      ALLOCATE(DIFF(NMIX))
      CALL LCMGET(KPMAC,'DIFF',DIFF)
*----
*  COMPUTE CORNER FLUXES
*----
      ALLOCATE(DELC(8,NX,NY,NZ))
      DELC(:8,:NX,:NY,:NZ)=0.D0
      IOFY=7*LL4F+LL4X
      IOFZ=7*LL4F+LL4X+LL4Y
      IF(ICORN==1) THEN
        ALLOCATE(FCORN(8,NX,NY,NZ))
        FCORN(:8,:NX,:NY,:NZ)=0.D0
        DO KS=1,NZ
          DO JS=1,NY
            DO IS=1,NX
              IEL=(KS-1)*NX*NY+(JS-1)*NX+IS
              IND1=KFLX(IEL)
              IF(IND1.EQ.0) CYCLE
              IBM=ISS(IEL)
              IF(IBM.LE.0) CYCLE
              JXM=KN(1,IS,JS,KS) ; JXP=KN(2,IS,JS,KS)
              JYM=KN(3,IS,JS,KS) ; JYP=KN(4,IS,JS,KS)
              JZM=KN(5,IS,JS,KS) ; JZP=KN(6,IS,JS,KS)
              COEFX=DIFF(IBM)/(XXX(IS+1)-XXX(IS))
              COEFY=DIFF(IBM)/(YYY(JS+1)-YYY(JS))
              COEFZ=DIFF(IBM)/(ZZZ(KS+1)-ZZZ(KS))
*
              WORK1(:,:)=0.0
              WORK1(1,1)=-0.5
              WORK1(1,2)=0.5
              WORK1(1,5)=EVT(LL4F+IND1)-EVT(IND1)
              WORK1(2,1)=0.5
              WORK1(2,2)=0.5
              WORK1(2,5)=EVT(2*LL4F+IND1)-EVT(IND1)
              WORK1(3,1)=-COEFX
              WORK1(3,2)=3.0*COEFX
              IF(JXM.NE.0) WORK1(3,5)=EVT(7*LL4F+JXM)
              WORK1(4,1)=-COEFX
              WORK1(4,2)=-3.0*COEFX
              IF(JXP.NE.0) WORK1(4,5)=EVT(7*LL4F+JXP)
              WORK1(3,3)=-0.5*COEFX
              WORK1(3,4)=0.2*COEFX
              WORK1(4,3)=-0.5*COEFX
              WORK1(4,4)=-0.2*COEFX
              CALL ALSBD(4,1,WORK1,IER,4)
              IF(IER.NE.0) CALL XABORT('VALU5: SINGULAR MATRIX(1).')
              DO IC=1,8
                SELECT CASE(IC)
                CASE(1,3,5,7)
                  U=-0.5
                CASE DEFAULT
                  U=0.5
                END SELECT
                GAR=EVT(IND1)+WORK1(1,5)*U+WORK1(2,5)*(3.0*U**2-0.25)
                GAR=GAR+WORK1(3,5)*(U**2-0.25)*U+WORK1(4,5)*(U**2-0.25)*
     1          (U**2-0.05)
                FCORN(IC,IS,JS,KS)=GAR
              ENDDO
*
              WORK1(:,:)=0.0
              WORK1(1,1)=-0.5
              WORK1(1,2)=0.5
              WORK1(1,5)=EVT(3*LL4F+IND1)-EVT(IND1)
              WORK1(2,1)=0.5
              WORK1(2,2)=0.5
              WORK1(2,5)=EVT(4*LL4F+IND1)-EVT(IND1)
              WORK1(3,1)=-COEFY
              WORK1(3,2)=3.0*COEFY
              IF(JYM.NE.0) WORK1(3,5)=EVT(IOFY+JYM)
              WORK1(4,1)=-COEFY
              WORK1(4,2)=-3.0*COEFY
              IF(JYP.NE.0) WORK1(4,5)=EVT(IOFY+JYP)
              WORK1(3,3)=-0.5*COEFY
              WORK1(3,4)=0.2*COEFY
              WORK1(4,3)=-0.5*COEFY
              WORK1(4,4)=-0.2*COEFY
              CALL ALSBD(4,1,WORK1,IER,4)
              IF(IER.NE.0) CALL XABORT('VALU5: SINGULAR MATRIX(2).')
              DO IC=1,8
                SELECT CASE(IC)
                CASE(1,2,5,6)
                  V=-0.5
                CASE DEFAULT
                  V=0.5
                END SELECT
                GAR=FCORN(IC,IS,JS,KS)+WORK1(1,5)*V+WORK1(2,5)*
     1          (3.0*V**2-0.25)
                GAR=GAR+WORK1(3,5)*(V**2-0.25)*V+WORK1(4,5)*(V**2-0.25)*
     1          (V**2-0.05)
                FCORN(IC,IS,JS,KS)=GAR
              ENDDO
*
              WORK1(:,:)=0.0
              WORK1(1,1)=-0.5
              WORK1(1,2)=0.5
              WORK1(1,5)=EVT(7*LL4F+IND1)-EVT(IND1)
              WORK1(2,1)=0.5
              WORK1(2,2)=0.5
              WORK1(2,5)=EVT(6*LL4F+IND1)-EVT(IND1)
              WORK1(3,1)=-COEFZ
              WORK1(3,2)=3.0*COEFZ
              IF(JZM.NE.0) WORK1(3,5)=EVT(IOFZ+JZM)
              WORK1(4,1)=-COEFZ
              WORK1(4,2)=-3.0*COEFZ
              IF(JZP.NE.0) WORK1(4,5)=EVT(IOFZ+JZP)
              WORK1(3,3)=-0.5*COEFZ
              WORK1(3,4)=0.2*COEFZ
              WORK1(4,3)=-0.5*COEFZ
              WORK1(4,4)=-0.2*COEFZ
              CALL ALSBD(4,1,WORK1,IER,4)
              IF(IER.NE.0) CALL XABORT('VALU5: SINGULAR MATRIX(3).')
              DO IC=1,8
                SELECT CASE(IC)
                CASE(1,2,3,4)
                  W=-0.5
                CASE DEFAULT
                  W=0.5
                END SELECT
                GAR=FCORN(IC,IS,JS,KS)+WORK1(1,5)*W+WORK1(2,5)*
     1          (3.0*W**2-0.25)
                GAR=GAR+WORK1(3,5)*(W**2-0.25)*W+WORK1(4,5)*(W**2-0.25)*
     1          (W**2-0.05)
                FCORN(IC,IS,JS,KS)=GAR
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO KS=1,NZ
          DO JS=1,NY
            DO IS=1,NX
              IEL=(KS-1)*NX*NY+(JS-1)*NX+IS
              IND1=KFLX(IEL)
              IF(IND1.EQ.0) CYCLE
              ! corner 1
              NB=1 ; GAR=FCORN(1,IS,JS,KS)
              LOGC1=(IS>1) ; LOGC2=(JS>1) ; LOGC3=(KS>1)
              IF(LOGC1) THEN
                IF(KFLX((KS-1)*NX*NY+(JS-1)*NX+IS-1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(2,IS-1,JS,KS)
                ENDIF
              ENDIF
              IF(LOGC2) THEN
                IF(KFLX((KS-1)*NX*NY+(JS-2)*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(3,IS,JS-1,KS)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC2) THEN
                IF(KFLX((KS-1)*NX*NY+(JS-2)*NX+IS-1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(4,IS-1,JS-1,KS)
                ENDIF
              ENDIF
              IF(LOGC3) THEN
                IF(KFLX((KS-2)*NX*NY+(JS-1)*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(5,IS,JS,KS-1)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC3) THEN
                IF(KFLX((KS-2)*NX*NY+(JS-1)*NX+IS-1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(6,IS-1,JS,KS-1)
                ENDIF
              ENDIF
              IF(LOGC2.AND.LOGC3) THEN
                IF(KFLX((KS-2)*NX*NY+(JS-2)*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(7,IS,JS-1,KS-1)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC2.AND.LOGC3) THEN
                IF(KFLX((KS-2)*NX*NY+(JS-2)*NX+IS-1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(8,IS-1,JS-1,KS-1)
                ENDIF
              ENDIF
              FC2(1)=GAR/REAL(NB)-FCORN(1,IS,JS,KS)
              ! corner 2
              NB=1 ; GAR=FCORN(2,IS,JS,KS)
              LOGC1=(IS<NX); LOGC2=(JS>1) ; LOGC3=(KS>1)
              IF(LOGC1) THEN
                IF(KFLX((KS-1)*NX*NY+(JS-1)*NX+IS+1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(1,IS+1,JS,KS)
                ENDIF
              ENDIF
              IF(LOGC2) THEN
                IF(KFLX((KS-1)*NX*NY+(JS-2)*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(4,IS,JS-1,KS)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC2) THEN
                IF(KFLX((KS-1)*NX*NY+(JS-2)*NX+IS+1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(3,IS+1,JS-1,KS)
                ENDIF
              ENDIF
              IF(LOGC3) THEN
                IF(KFLX((KS-2)*NX*NY+(JS-1)*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(6,IS,JS,KS-1)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC3) THEN
                IF(KFLX((KS-2)*NX*NY+(JS-1)*NX+IS+1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(5,IS+1,JS,KS-1)
                ENDIF
              ENDIF
              IF(LOGC2.AND.LOGC3) THEN
                IF(KFLX((KS-2)*NX*NY+(JS-2)*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(8,IS,JS-1,KS-1)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC2.AND.LOGC3) THEN
                IF(KFLX((KS-2)*NX*NY+(JS-2)*NX+IS+1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(7,IS+1,JS-1,KS-1)
                ENDIF
              ENDIF
              FC2(2)=GAR/REAL(NB)-FCORN(2,IS,JS,KS)
              ! corner 3
              NB=1 ; GAR=FCORN(3,IS,JS,KS)
              LOGC1=(IS>1) ;  LOGC2=(JS<NY) ; LOGC3=(KS>1)
              IF(LOGC1) THEN
                IF(KFLX((KS-1)*NX*NY+(JS-1)*NX+IS-1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(4,IS-1,JS,KS)
                ENDIF
              ENDIF
              IF(LOGC2) THEN
                IF(KFLX((KS-1)*NX*NY+JS*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(1,IS,JS+1,KS)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC2) THEN
                IF(KFLX((KS-1)*NX*NY+JS*NX+IS-1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(2,IS-1,JS+1,KS)
                ENDIF
              ENDIF
              IF(LOGC3) THEN
                IF(KFLX((KS-2)*NX*NY+(JS-1)*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(7,IS,JS,KS-1)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC3) THEN
                IF(KFLX((KS-2)*NX*NY+(JS-1)*NX+IS-1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(8,IS-1,JS,KS-1)
                ENDIF
              ENDIF
              IF(LOGC2.AND.LOGC3) THEN
                IF(KFLX((KS-2)*NX*NY+JS*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(5,IS,JS+1,KS-1)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC2.AND.LOGC3) THEN
                IF(KFLX((KS-2)*NX*NY+JS*NX+IS-1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(6,IS-1,JS+1,KS-1)
                ENDIF
              ENDIF
              FC2(3)=GAR/REAL(NB)-FCORN(3,IS,JS,KS)
              ! corner 4
              NB=1 ; GAR=FCORN(4,IS,JS,KS)
              LOGC1=(IS<NX) ; LOGC2=(JS<NY) ; LOGC3=(KS>1)
              IF(LOGC1) THEN
                IF(KFLX((KS-1)*NX*NY+(JS-1)*NX+IS+1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(3,IS+1,JS,KS)
                ENDIF
              ENDIF
              IF(LOGC2) THEN
                IF(KFLX((KS-1)*NX*NY+JS*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(2,IS,JS+1,KS)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC2) THEN
                IF(KFLX((KS-1)*NX*NY+JS*NX+IS+1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(1,IS+1,JS+1,KS)
                ENDIF
              ENDIF
              IF(LOGC3) THEN
                IF(KFLX((KS-2)*NX*NY+(JS-1)*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(8,IS,JS,KS-1)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC3) THEN
                IF(KFLX((KS-2)*NX*NY+(JS-1)*NX+IS+1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(7,IS+1,JS,KS-1)
                ENDIF
              ENDIF
              IF(LOGC2.AND.LOGC3) THEN
                IF(KFLX((KS-2)*NX*NY+JS*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(6,IS,JS+1,KS-1)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC2.AND.LOGC3) THEN
                IF(KFLX((KS-2)*NX*NY+JS*NX+IS+1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(5,IS+1,JS+1,KS-1)
                ENDIF
              ENDIF
              FC2(4)=GAR/REAL(NB)-FCORN(4,IS,JS,KS)
              ! corner 5
              NB=1 ; GAR=FCORN(5,IS,JS,KS)
              LOGC1=(IS>1) ; LOGC2=(JS>1) ; LOGC3=(KS<NZ)
              IF(LOGC1) THEN
                IF(KFLX((KS-1)*NX*NY+(JS-1)*NX+IS-1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(6,IS-1,JS,KS)
                ENDIF
              ENDIF
              IF(LOGC2) THEN
                IF(KFLX((KS-1)*NX*NY+(JS-2)*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(7,IS,JS-1,KS)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC2) THEN
                IF(KFLX((KS-1)*NX*NY+(JS-2)*NX+IS-1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(8,IS-1,JS-1,KS)
                ENDIF
              ENDIF
              IF(LOGC3) THEN
                IF(KFLX(KS*NX*NY+(JS-1)*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(1,IS,JS,KS+1)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC3) THEN
                IF(KFLX(KS*NX*NY+(JS-1)*NX+IS-1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(2,IS-1,JS,KS+1)
                ENDIF
              ENDIF
              IF(LOGC2.AND.LOGC3) THEN
                IF(KFLX(KS*NX*NY+(JS-2)*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(3,IS,JS-1,KS+1)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC2.AND.LOGC3) THEN
                IF(KFLX(KS*NX*NY+(JS-2)*NX+IS-1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(4,IS-1,JS-1,KS+1)
                ENDIF
              ENDIF
              FC2(5)=GAR/REAL(NB)-FCORN(5,IS,JS,KS)
              ! corner 6
              NB=1 ; GAR=FCORN(6,IS,JS,KS)
              LOGC1=(IS<NX); LOGC2=(JS>1) ; LOGC3=(KS<NZ)
              IF(LOGC1) THEN
                IF(KFLX((KS-1)*NX*NY+(JS-1)*NX+IS+1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(5,IS+1,JS,KS)
                ENDIF
              ENDIF
              IF(LOGC2) THEN
                IF(KFLX((KS-1)*NX*NY+(JS-2)*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(8,IS,JS-1,KS)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC2) THEN
                IF(KFLX((KS-1)*NX*NY+(JS-2)*NX+IS+1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(7,IS+1,JS-1,KS)
                ENDIF
              ENDIF
              IF(LOGC3) THEN
                IF(KFLX(KS*NX*NY+(JS-1)*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(2,IS,JS,KS+1)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC3) THEN
                IF(KFLX(KS*NX*NY+(JS-1)*NX+IS+1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(1,IS+1,JS,KS+1)
                ENDIF
              ENDIF
              IF(LOGC2.AND.LOGC3) THEN
                IF(KFLX(KS*NX*NY+(JS-2)*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(4,IS,JS-1,KS+1)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC2.AND.LOGC3) THEN
                IF(KFLX(KS*NX*NY+(JS-2)*NX+IS+1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(3,IS+1,JS-1,KS+1)
                ENDIF
              ENDIF
              FC2(6)=GAR/REAL(NB)-FCORN(6,IS,JS,KS)
              ! corner 7
              NB=1 ; GAR=FCORN(7,IS,JS,KS)
              LOGC1=(IS>1) ;  LOGC2=(JS<NY) ; LOGC3=(KS<NZ)
              IF(LOGC1) THEN
                IF(KFLX((KS-1)*NX*NY+(JS-1)*NX+IS-1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(8,IS-1,JS,KS)
                ENDIF
              ENDIF
              IF(LOGC2) THEN
                IF(KFLX((KS-1)*NX*NY+JS*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(5,IS,JS+1,KS)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC2) THEN
                IF(KFLX((KS-1)*NX*NY+JS*NX+IS-1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(6,IS-1,JS+1,KS)
                ENDIF
              ENDIF
              IF(LOGC3) THEN
                IF(KFLX(KS*NX*NY+(JS-1)*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(3,IS,JS,KS+1)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC3) THEN
                IF(KFLX(KS*NX*NY+(JS-1)*NX+IS-1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(4,IS-1,JS,KS+1)
                ENDIF
              ENDIF
              IF(LOGC2.AND.LOGC3) THEN
                IF(KFLX(KS*NX*NY+JS*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(1,IS,JS+1,KS+1)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC2.AND.LOGC3) THEN
                IF(KFLX(KS*NX*NY+JS*NX+IS-1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(2,IS-1,JS+1,KS+1)
                ENDIF
              ENDIF
              FC2(7)=GAR/REAL(NB)-FCORN(7,IS,JS,KS)
              ! corner 8
              NB=1 ; GAR=FCORN(8,IS,JS,KS)
              LOGC1=(IS<NX) ; LOGC2=(JS<NY) ; LOGC3=(KS<NZ)
              IF(LOGC1) THEN
                IF(KFLX((KS-1)*NX*NY+(JS-1)*NX+IS+1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(7,IS+1,JS,KS)
                ENDIF
              ENDIF
              IF(LOGC2) THEN
                IF(KFLX((KS-1)*NX*NY+JS*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(6,IS,JS+1,KS)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC2) THEN
                IF(KFLX((KS-1)*NX*NY+JS*NX+IS+1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(5,IS+1,JS+1,KS)
                ENDIF
              ENDIF
              IF(LOGC3) THEN
                IF(KFLX(KS*NX*NY+(JS-1)*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(4,IS,JS,KS+1)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC3) THEN
                IF(KFLX(KS*NX*NY+(JS-1)*NX+IS+1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(3,IS+1,JS,KS+1)
                ENDIF
              ENDIF
              IF(LOGC2.AND.LOGC3) THEN
                IF(KFLX(KS*NX*NY+JS*NX+IS)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(2,IS,JS+1,KS+1)
                ENDIF
              ENDIF
              IF(LOGC1.AND.LOGC2.AND.LOGC3) THEN
                IF(KFLX(KS*NX*NY+JS*NX+IS+1)>0) THEN
                  NB=NB+1 ; GAR=GAR+FCORN(1,IS+1,JS+1,KS+1)
                ENDIF
              ENDIF
              FC2(8)=GAR/REAL(NB)-FCORN(8,IS,JS,KS)
              ! polynomial coefficients of correction terms
              DELC(1,IS,JS,KS)=-FC2(1)+FC2(2)+FC2(3)-FC2(4)+FC2(5)-
     1        FC2(6)-FC2(7)+FC2(8)
              DELC(2,IS,JS,KS)= FC2(1)+FC2(2)-FC2(3)-FC2(4)-FC2(5)-
     1        FC2(6)+FC2(7)+FC2(8)
              DELC(3,IS,JS,KS)= FC2(1)-FC2(2)+FC2(3)-FC2(4)-FC2(5)+
     1        FC2(6)-FC2(7)+FC2(8)
              DELC(4,IS,JS,KS)=-FC2(1)-FC2(2)-FC2(3)-FC2(4)+FC2(5)+
     1        FC2(6)+FC2(7)+FC2(8)
              DELC(5,IS,JS,KS)= FC2(1)-FC2(2)-FC2(3)+FC2(4)+FC2(5)-
     1        FC2(6)-FC2(7)+FC2(8)
              DELC(6,IS,JS,KS)=-FC2(1)-FC2(2)+FC2(3)+FC2(4)-FC2(5)-
     1        FC2(6)+FC2(7)+FC2(8)
              DELC(7,IS,JS,KS)=-FC2(1)+FC2(2)-FC2(3)+FC2(4)-FC2(5)+
     1        FC2(6)-FC2(7)+FC2(8)
              DELC(8,IS,JS,KS)= FC2(1)+FC2(2)+FC2(3)+FC2(4)+FC2(5)+
     1        FC2(6)+FC2(7)+FC2(8)
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(FCORN)
      ENDIF
*----
*  PERFORM INTERPOLATION
*----
      DO K=1,IZLG
        COTE=Z(K)
        DO J=1,IYLG
          ORDO=Y(J)
          DO I=1,IXLG
            ABSC=X(I)
            GAR=0.0D0
               AXYZ(I,J,K)=REAL(GAR)
*                                                          
*           Find the node index containing the interpolation point
            IS=0
            JS=0
            KS=0
            DO L=1,NX
              IS=L
              IF((ABSC.GE.XXX(L)).AND.(ABSC.LE.XXX(L+1))) GO TO 10
            ENDDO
            CALL XABORT('VALU5: WRONG INTERPOLATION(1).')
   10       DO L=1,NY
              JS=L
              IF((ORDO.GE.YYY(L)).AND.(ORDO.LE.YYY(L+1))) GO TO 20
            ENDDO
            CALL XABORT('VALU5: WRONG INTERPOLATION(2).')
   20       DO L=1,NZ
              KS=L
              IF((COTE.GE.ZZZ(L)).AND.(COTE.LE.ZZZ(L+1))) GO TO 30
            ENDDO
            CALL XABORT('VALU5: WRONG INTERPOLATION(3).')
   30       IEL=(KS-1)*NX*NY+(JS-1)*NX+IS
            IND1=KFLX(IEL)
            IF(IND1.EQ.0) GO TO 40
            IBM=ISS(IEL)
            IF(IBM.LE.0) GO TO 40
            JXM=KN(1,IS,JS,KS) ; JXP=KN(2,IS,JS,KS)
            JYM=KN(3,IS,JS,KS) ; JYP=KN(4,IS,JS,KS)
            JZM=KN(5,IS,JS,KS) ; JZP=KN(6,IS,JS,KS)
            COEFX=DIFF(IBM)/(XXX(IS+1)-XXX(IS))
            COEFY=DIFF(IBM)/(YYY(JS+1)-YYY(JS))
            COEFZ=DIFF(IBM)/(ZZZ(KS+1)-ZZZ(KS))
            U=(ABSC-XXX(IS))/(XXX(IS+1)-XXX(IS))-0.5
            V=(ORDO-YYY(JS))/(YYY(JS+1)-YYY(JS))-0.5
            W=(COTE-ZZZ(KS))/(ZZZ(KS+1)-ZZZ(KS))-0.5
            GAR=EVT(IND1)
*
            WORK1(:,:)=0.0
            WORK1(1,1)=-0.5
            WORK1(1,2)=0.5
            WORK1(1,5)=EVT(LL4F+IND1)-EVT(IND1)
            WORK1(2,1)=0.5
            WORK1(2,2)=0.5
            WORK1(2,5)=EVT(2*LL4F+IND1)-EVT(IND1)
            WORK1(3,1)=-COEFX
            WORK1(3,2)=3.0*COEFX
            IF(JXM.NE.0) WORK1(3,5)=EVT(7*LL4F+JXM)
            WORK1(4,1)=-COEFX
            WORK1(4,2)=-3.0*COEFX
            IF(JXP.NE.0) WORK1(4,5)=EVT(7*LL4F+JXP)
            WORK1(3,3)=-0.5*COEFX
            WORK1(3,4)=0.2*COEFX
            WORK1(4,3)=-0.5*COEFX
            WORK1(4,4)=-0.2*COEFX
            CALL ALSBD(4,1,WORK1,IER,4)
            IF(IER.NE.0) CALL XABORT('VALU5: SINGULAR MATRIX(4).')
            GAR=GAR+WORK1(1,5)*U+WORK1(2,5)*(3.0*U**2-0.25)
            GAR=GAR+WORK1(3,5)*(U**2-0.25)*U+WORK1(4,5)*(U**2-0.25)*
     1     (U**2-0.05)
*
            WORK1(:,:)=0.0
            WORK1(1,1)=-0.5
            WORK1(1,2)=0.5
            WORK1(1,5)=EVT(3*LL4F+IND1)-EVT(IND1)
            WORK1(2,1)=0.5
            WORK1(2,2)=0.5
            WORK1(2,5)=EVT(4*LL4F+IND1)-EVT(IND1)
            WORK1(3,1)=-COEFY
            WORK1(3,2)=3.0*COEFY
            IF(JYM.NE.0) WORK1(3,5)=EVT(IOFY+JYM)
            WORK1(4,1)=-COEFY
            WORK1(4,2)=-3.0*COEFY
            IF(JYP.NE.0) WORK1(4,5)=EVT(IOFY+JYP)
            WORK1(3,3)=-0.5*COEFY
            WORK1(3,4)=0.2*COEFY
            WORK1(4,3)=-0.5*COEFY
            WORK1(4,4)=-0.2*COEFY
            CALL ALSBD(4,1,WORK1,IER,4)
            IF(IER.NE.0) CALL XABORT('VALU5: SINGULAR MATRIX(5).')
            GAR=GAR+WORK1(1,5)*V+WORK1(2,5)*(3.0*V**2-0.25)
            GAR=GAR+WORK1(3,5)*(V**2-0.25)*V+WORK1(4,5)*(V**2-0.25)*
     1      (V**2-0.05)
*
            WORK1(:,:)=0.0
            WORK1(1,1)=-0.5
            WORK1(1,2)=0.5
            WORK1(1,5)=EVT(5*LL4F+IND1)-EVT(IND1)
            WORK1(2,1)=0.5
            WORK1(2,2)=0.5
            WORK1(2,5)=EVT(6*LL4F+IND1)-EVT(IND1)
            WORK1(3,1)=-COEFZ
            WORK1(3,2)=3.0*COEFZ
            IF(JZM.NE.0) WORK1(3,5)=EVT(IOFZ+JZM)
            WORK1(4,1)=-COEFZ
            WORK1(4,2)=-3.0*COEFZ
            IF(JZP.NE.0) WORK1(4,5)=EVT(IOFZ+JZP)
            WORK1(3,3)=-0.5*COEFZ
            WORK1(3,4)=0.2*COEFZ
            WORK1(4,3)=-0.5*COEFZ
            WORK1(4,4)=-0.2*COEFZ
            CALL ALSBD(4,1,WORK1,IER,4)
            IF(IER.NE.0) CALL XABORT('VALU5: SINGULAR MATRIX(6).')
            GAR=GAR+WORK1(1,5)*W+WORK1(2,5)*(3.0*W**2-0.25)
            GAR=GAR+WORK1(3,5)*(W**2-0.25)*W+WORK1(4,5)*(W**2-0.25)*
     1      (W**2-0.05)
*
            IF(ICORN==1) THEN
              ! perform interpolation of corner flux correction
              P2U=3.0*U**2-0.25 ; P2V=3.0*V**2-0.25 ; P2W=3.0*W**2-0.25
              GAR=GAR+DELC(1,IS,JS,KS)*U*V*W + DELC(2,IS,JS,KS)*P2U*V*W+
     1        DELC(3,IS,JS,KS)*U*P2V*W + DELC(4,IS,JS,KS)*P2U*P2V*W+
     2        DELC(5,IS,JS,KS)*U*V*P2W + DELC(6,IS,JS,KS)*P2U*V*P2W+
     3        DELC(7,IS,JS,KS)*U*P2V*P2W + DELC(8,IS,JS,KS)*P2U*P2V*P2W
            ENDIF
   40       AXYZ(I,J,K)=REAL(GAR)
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(DELC,DIFF)
      RETURN
      END
