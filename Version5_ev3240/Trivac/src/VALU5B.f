*DECK VALU5B
      SUBROUTINE VALU5B (KPMAC,NX,NY,LL4F,LL4X,NUN,NMIX,X,Y,XXX,YYY,
     1 EVT,ISS,KFLX,KN,IXLG,IYLG,ICORN,AXY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Interpolation of the flux distribution for nodal method in 2D.
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
* LL4F    number of averaged flux unknowns.
* LL4X    number of X-directed net currents.
* NUN     dimension of unknown array EVT.
* NMIX    number of mixtures.
* X       Cartesian coordinates along the X axis where the flux is
*         interpolated.
* Y       Cartesian coordinates along the Y axis where the flux is
*         interpolated.
* XXX     Cartesian coordinates along the X axis.
* YYY     Cartesian coordinates along the Y axis.
* EVT     reconstruction coefficients of the flux.
* ISS     mixture index assigned to each element.
* KFLX    correspondence between local and global numbering.
* KN      element-ordered interface net current unknown list.
* IXLG    number of interpolated points according to X.
* IYLG    number of interpolated points according to Y.
* ICORN   flag to activate corner flux correction (0/1: OFF/ON).
*                                                                      
*Parameters: output
* AXY     interpolated fluxes.
*                                                                      
*----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) KPMAC
      INTEGER NX,NY,LL4F,LL4X,NUN,NMIX,ISS(NX*NY),KFLX(NX*NY),
     1 KN(6,NX,NY),IXLG,IYLG,ICORN
      REAL X(IXLG),Y(IYLG),XXX(NX+1),YYY(NY+1),EVT(NUN),AXY(IXLG,IYLG)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION  WORK1(4,5),FC2(4)
      DOUBLE PRECISION GAR,COEFX,COEFY,U,V,P2U,P2V
      LOGICAL LOGC1,LOGC2
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: DIFF
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: FCORN,DELC
*----
*  RECOVER DIFFUSION COEFFICIENTS
*----
      ALLOCATE(DIFF(NMIX))
      CALL LCMGET(KPMAC,'DIFF',DIFF)
*----
*  COMPUTE CORNER FLUXES
*----
      ALLOCATE(DELC(4,NX,NY))
      DELC(:4,:NX,:NY)=0.D0
      IF(ICORN==1) THEN
        ALLOCATE(FCORN(4,NX,NY))
        FCORN(:4,:NX,:NY)=0.D0
        DO JS=1,NY
          DO IS=1,NX
            IEL=(JS-1)*NX+IS
            IND1=KFLX(IEL)
            IF(IND1.EQ.0) CYCLE
            IBM=ISS(IEL)
            IF(IBM.LE.0) CYCLE
            JXM=KN(1,IS,JS) ; JXP=KN(2,IS,JS)
            JYM=KN(3,IS,JS) ; JYP=KN(4,IS,JS)
            COEFX=DIFF(IBM)/(XXX(IS+1)-XXX(IS))
            COEFY=DIFF(IBM)/(YYY(JS+1)-YYY(JS))
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
            IF(JXM.NE.0) WORK1(3,5)=EVT(5*LL4F+JXM)
            WORK1(4,1)=-COEFX
            WORK1(4,2)=-3.0*COEFX
            IF(JXP.NE.0) WORK1(4,5)=EVT(5*LL4F+JXP)
            WORK1(3,3)=-0.5*COEFX
            WORK1(3,4)=0.2*COEFX
            WORK1(4,3)=-0.5*COEFX
            WORK1(4,4)=-0.2*COEFX
            CALL ALSBD(4,1,WORK1,IER,4)
            IF(IER.NE.0) CALL XABORT('VALU5B: SINGULAR MATRIX(1).')
            DO IC=1,4
              SELECT CASE(IC)
              CASE(1,3)
                U=-0.5
              CASE DEFAULT
                U=0.5
              END SELECT
              GAR=EVT(IND1)+WORK1(1,5)*U+WORK1(2,5)*(3.0*U**2-0.25)
              GAR=GAR+WORK1(3,5)*(U**2-0.25)*U+WORK1(4,5)*(U**2-0.25)*
     1        (U**2-0.05)
              FCORN(IC,IS,JS)=GAR
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
            IF(JYM.NE.0) WORK1(3,5)=EVT(5*LL4F+LL4X+JYM)
            WORK1(4,1)=-COEFY
            WORK1(4,2)=-3.0*COEFY
            IF(JYP.NE.0) WORK1(4,5)=EVT(5*LL4F+LL4X+JYP)
            WORK1(3,3)=-0.5*COEFY
            WORK1(3,4)=0.2*COEFY
            WORK1(4,3)=-0.5*COEFY
            WORK1(4,4)=-0.2*COEFY
            CALL ALSBD(4,1,WORK1,IER,4)
            IF(IER.NE.0) CALL XABORT('VALU5B: SINGULAR MATRIX(2).')
            DO IC=1,4
              SELECT CASE(IC)
              CASE(1,2)
                V=-0.5
              CASE DEFAULT
                V=0.5
              END SELECT
              GAR=FCORN(IC,IS,JS)+WORK1(1,5)*V+WORK1(2,5)*
     1        (3.0*V**2-0.25)
              GAR=GAR+WORK1(3,5)*(V**2-0.25)*V+WORK1(4,5)*(V**2-0.25)*
     1        (V**2-0.05)
              FCORN(IC,IS,JS)=GAR
            ENDDO
          ENDDO
        ENDDO
        DO JS=1,NY
          DO IS=1,NX
            IEL=(JS-1)*NX+IS
            IND1=KFLX(IEL)
            IF(IND1.EQ.0) CYCLE
            ! corner 1
            NB=1; GAR=FCORN(1,IS,JS)
            LOGC1=(IS>1) ; LOGC2=(JS>1)
            IF(LOGC2) LOGC2=(KFLX((JS-2)*NX+IS)>0)
            IF(LOGC1) THEN
              IF(KFLX((JS-1)*NX+IS-1)>0) THEN
                NB=NB+1 ;GAR=GAR+FCORN(2,IS-1,JS)
              ENDIF
            ENDIF
            IF(LOGC2) THEN
              IF(KFLX((JS-2)*NX+IS)>0) THEN
                NB=NB+1 ;GAR=GAR+FCORN(3,IS,JS-1)
              ENDIF
            ENDIF
            IF(LOGC1.AND.LOGC2) THEN
              IF(KFLX((JS-2)*NX+IS-1)>0) THEN
                NB=NB+1 ;GAR=GAR+FCORN(4,IS-1,JS-1)
              ENDIF
            ENDIF
            FC2(1)=GAR/REAL(NB)-FCORN(1,IS,JS)
            ! corner 2
            NB=1 ;GAR=FCORN(2,IS,JS)
            LOGC1=(IS<NX) ; LOGC2=(JS>1)
            IF(LOGC1) THEN
              IF(KFLX((JS-1)*NX+IS+1)>0) THEN
                NB=NB+1 ;GAR=GAR+FCORN(1,IS+1,JS)
              ENDIF
            ENDIF
            IF(LOGC2) THEN
              IF(KFLX((JS-2)*NX+IS)>0) THEN
                NB=NB+1 ;GAR=GAR+FCORN(4,IS,JS-1)
              ENDIF
            ENDIF
            IF(LOGC1.AND.LOGC2) THEN
              IF(KFLX((JS-2)*NX+IS+1)>0) THEN
                NB=NB+1 ;GAR=GAR+FCORN(3,IS+1,JS-1)
              ENDIF
            ENDIF
            FC2(2)=GAR/REAL(NB)-FCORN(2,IS,JS)
            ! corner 3
            NB=1 ; GAR=FCORN(3,IS,JS)
            LOGC1=(IS>1) ; LOGC2=(JS<NY)
            IF(LOGC1) THEN
              IF(KFLX((JS-1)*NX+IS-1)>0) THEN
                NB=NB+1 ; GAR=GAR+FCORN(4,IS-1,JS)
              ENDIF
            ENDIF
            IF(LOGC2) THEN
              IF(KFLX(JS*NX+IS)>0) THEN
                NB=NB+1 ; GAR=GAR+FCORN(1,IS,JS+1)
              ENDIF
            ENDIF
            IF(LOGC1.AND.LOGC2) THEN
              IF(KFLX(JS*NX+IS-1)>0) THEN
                NB=NB+1 ; GAR=GAR+FCORN(2,IS-1,JS+1)
              ENDIF
            ENDIF
            FC2(3)=GAR/REAL(NB)-FCORN(3,IS,JS)
            ! corner 4
            NB=1
            GAR=FCORN(4,IS,JS)
            LOGC1=(IS<NX)
            IF(LOGC1) LOGC1=(KFLX((JS-1)*NX+IS+1)>0)
            LOGC2=(JS<NY)
            IF(LOGC2) LOGC2=(KFLX(JS*NX+IS)>0)
            IF(LOGC1) THEN
              IF(KFLX((JS-1)*NX+IS+1)>0) THEN
                NB=NB+1 ; GAR=GAR+FCORN(3,IS+1,JS)
              ENDIF
            ENDIF
            IF(LOGC2) THEN
              IF(KFLX(JS*NX+IS)>0) THEN
                NB=NB+1 ; GAR=GAR+FCORN(2,IS,JS+1)
              ENDIF
            ENDIF
            IF(LOGC1.AND.LOGC2) THEN
              IF(KFLX(JS*NX+IS+1)>0) THEN
                NB=NB+1 ; GAR=GAR+FCORN(1,IS+1,JS+1)
              ENDIF
            ENDIF
            FC2(4)=GAR/REAL(NB)-FCORN(4,IS,JS)
            ! polynomial coefficients of correction terms
            DELC(1,IS,JS)= FC2(1)-FC2(2)-FC2(3)+FC2(4)
            DELC(2,IS,JS)=-FC2(1)-FC2(2)+FC2(3)+FC2(4)
            DELC(3,IS,JS)=-FC2(1)+FC2(2)-FC2(3)+FC2(4)
            DELC(4,IS,JS)= FC2(1)+FC2(2)+FC2(3)+FC2(4)
          ENDDO
        ENDDO
        DEALLOCATE(FCORN)
      ENDIF
*----
*  PERFORM INTERPOLATION
*----
      DO J=1,IYLG
        ORDO=Y(J)
        DO I=1,IXLG
          ABSC=X(I)
          GAR=0.0D0
*                                                          
*         Find the node index containing the interpolation point
          IS=0; JS=0
          DO L=1,NX
            IS=L
            IF((ABSC.GE.XXX(L)).AND.(ABSC.LE.XXX(L+1))) GO TO 10
          ENDDO
          CALL XABORT('VALU5B: WRONG INTERPOLATION(1).')
   10     DO L=1,NY
            JS=L
            IF((ORDO.GE.YYY(L)).AND.(ORDO.LE.YYY(L+1))) GO TO 20
          ENDDO
          CALL XABORT('VALU5B: WRONG INTERPOLATION(2).')
   20     IEL=(JS-1)*NX+IS
          IND1=KFLX(IEL)
          IF(IND1.EQ.0) GO TO 30
          IBM=ISS(IEL)
          IF(IBM.LE.0) GO TO 30
          JXM=KN(1,IS,JS) ; JXP=KN(2,IS,JS)
          JYM=KN(3,IS,JS) ; JYP=KN(4,IS,JS)
          COEFX=DIFF(IBM)/(XXX(IS+1)-XXX(IS))
          COEFY=DIFF(IBM)/(YYY(JS+1)-YYY(JS))
          U=(ABSC-XXX(IS))/(XXX(IS+1)-XXX(IS))-0.5
          V=(ORDO-YYY(JS))/(YYY(JS+1)-YYY(JS))-0.5
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
          IF(JXM.NE.0) WORK1(3,5)=EVT(5*LL4F+JXM)
          WORK1(4,1)=-COEFX
          WORK1(4,2)=-3.0*COEFX
          IF(JXP.NE.0) WORK1(4,5)=EVT(5*LL4F+JXP)
          WORK1(3,3)=-0.5*COEFX
          WORK1(3,4)=0.2*COEFX
          WORK1(4,3)=-0.5*COEFX
          WORK1(4,4)=-0.2*COEFX
          CALL ALSBD(4,1,WORK1,IER,4)
          IF(IER.NE.0) CALL XABORT('VALU5B: SINGULAR MATRIX(3).')
          GAR=GAR+WORK1(1,5)*U+WORK1(2,5)*(3.0*U**2-0.25)
          GAR=GAR+WORK1(3,5)*(U**2-0.25)*U+WORK1(4,5)*(U**2-0.25)*
     1    (U**2-0.05)
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
          IF(JYM.NE.0) WORK1(3,5)=EVT(5*LL4F+LL4X+JYM)
          WORK1(4,1)=-COEFY
          WORK1(4,2)=-3.0*COEFY
          IF(JYP.NE.0) WORK1(4,5)=EVT(5*LL4F+LL4X+JYP)
          WORK1(3,3)=-0.5*COEFY
          WORK1(3,4)=0.2*COEFY
          WORK1(4,3)=-0.5*COEFY
          WORK1(4,4)=-0.2*COEFY
          CALL ALSBD(4,1,WORK1,IER,4)
          IF(IER.NE.0) CALL XABORT('VALU5B: SINGULAR MATRIX(4).')
          GAR=GAR+WORK1(1,5)*V+WORK1(2,5)*(3.0*V**2-0.25)
          GAR=GAR+WORK1(3,5)*(V**2-0.25)*V+WORK1(4,5)*(V**2-0.25)*
     1    (V**2-0.05)
*
          IF(ICORN==1) THEN
            ! perform interpolation of corner flux correction
            P2U=3.0*U**2-0.25 ; P2V=3.0*V**2-0.25
            GAR=GAR+DELC(1,IS,JS)*U*V + DELC(2,IS,JS)*P2U*V+
     1      DELC(3,IS,JS)*U*P2V + DELC(4,IS,JS)*P2U*P2V
          ENDIF
   30     AXY(I,J)=REAL(GAR)
        ENDDO
      ENDDO
      DEALLOCATE(DELC,DIFF)
      RETURN
      END
