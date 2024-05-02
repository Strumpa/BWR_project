*DECK NSSFL5
      SUBROUTINE NSSFL5(IPFLX,NUN,NG,NX,NY,NZ,LL4F,LL4X,LL4Y,LL4Z,NMIX,
     > NALB,ICL1,ICL2,NADI,EPSNOD,MAXNOD,EPSTHR,MAXTHR,EPSOUT,MAXOUT,
     > MAT,XX,YY,ZZ,XXX,YYY,ZZZ,IDL,VOL,KN,IQFR,QFR,DIFF,SIGR,CHI,SIGF,
     > SCAT,BETA,FD,BNDTL,NPASS,MUX,MUY,MUZ,IMAX,IMAY,IMAZ,IPY,IPZ,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Flux calculation for the analytic nodal method in Cartesian 3D
* geometry using the nodal correction iteration strategy.
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
* IPFLX   nodal flux.
* NUN     number of unknowns per energy group.
* NG      number of energy groups.
* NX      number of nodes in the X direction.
* NY      number of nodes in the Y direction.
* NZ      number of nodes in the Z direction.
* LL4F    number of nodal flux unknowns.
* LL4X    number of nodal X-directed net currents unknowns.
* LL4Y    number of nodal Y-directed net currents unknowns.
* LL4Z    number of nodal Z-directed net currents unknowns.
* NMIX    number of mixtures in the nodal calculation.
* NALB    number of physical albedos.
* ICL1    number of free iterations in one cycle of the inverse power
*         method (used for thermal iterations).
* ICL2    number of accelerated iterations in one cycle.
* NADI    number of inner ADI iterations.
* EPSNOD  nodal correction epsilon.
* MAXNOD  maximum number of nodal correction iterations.
* EPSTHR  thermal iteration epsilon.
* MAXTHR  maximum number of thermal iterations.
* EPSOUT  convergence epsilon for the power method.
* MAXOUT  maximum number of iterations for the power method.
* MAT     material mixtures.
* XX      mesh spacings in the X direction.
* YY      mesh spacings in the Y direction.
* ZZ      mesh spacings in the Z direction.
* XXX     Cartesian coordinates along the X axis.
* YYY     Cartesian coordinates along the Y axis.
* ZZZ     Cartesian coordinates along the Z axis.
* IDL     position of averaged fluxes in unknown vector.
* VOL     node volumes.
* KN      node-ordered interface net current unknown list.
* IQFR    boundary condition information.
* QFR     albedo function information.
* DIFF    diffusion coefficients
* SIGR    removal cross sections.
* CHI     fission spectra.
* SIGF    nu times fission cross section.
* SCAT    scattering cross section.
* BETA    albedos.
* FD      discontinuity factors.
* BNDTL   set to 'flat' or 'quadratic'.
* NPASS   number of transverse current iterations.
* MUX     X-oriented compressed storage mode indices.
* MUY     Y-oriented compressed storage mode indices.
* MUZ     Z-oriented compressed storage mode indices.
* IMAX    X-oriented position of each first non-zero column element.
* IMAY    Y-oriented position of each first non-zero column element.
* IMAZ    Z-oriented position of each first non-zero column element.
* IPY     Y-oriented permutation matrices.
* IPZ     Z-oriented permutation matrices.
* IMPX    edition flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPFLX
      INTEGER NUN,NG,NX,NY,LL4F,LL4X,LL4Y,NMIX,NALB,ICL1,ICL2,
     1 NADI,MAXNOD,MAXTHR,MAXOUT,IMPX,MAT(NX*NY*NZ),IDL(NX*NY*NZ),
     2 KN(6,NX,NY,NZ),IQFR(6,NX,NY,NZ),NPASS,MUX(LL4F),MUY(LL4F),
     3 MUZ(LL4F),IMAX(LL4F),IMAY(LL4F),IMAZ(LL4F),IPY(LL4F),IPZ(LL4F)
      REAL EPSNOD,EPSTHR,EPSOUT,XX(NX*NY*NZ),YY(NX*NY*NZ),ZZ(NX*NY*NZ),
     1 XXX(NX+1),YYY(NY+1),ZZZ(NZ+1),VOL(NX*NY*NZ),QFR(6,NX*NY*NZ),
     2 DIFF(NMIX,NG),SIGR(NMIX,NG),CHI(NMIX,NG),SIGF(NMIX,NG),
     3 SCAT(NMIX,NG,NG),BETA(NALB,NG,NG),FD(NMIX,6,NG,NG)
      CHARACTER(LEN=12) :: BNDTL
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPFLX
      INTEGER, PARAMETER :: NDIM=3
      REAL :: COEF(6),KEFF,KEFF_OLD
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: EVECT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: A11X,A11Y,A11Z
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: QFR2,DRIFT
      REAL, POINTER, DIMENSION(:,:) :: SAVG
*
      ALB(X)=0.5*(1.0-X)/(1.0+X)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      NEL=NX*NY*NZ
      N=LL4F*NG
      ALLOCATE(QFR2(6,NEL,NG),EVECT(N))
      ALLOCATE(DRIFT(6,NEL,NG),SAVG(NUN,NG))
*----
*  ALBEDO PROCESSING
*----
      QFR2(:6,:NEL,:NG)=0.0
      DO IG=1,NG
        DO IQW=1,6
          DO I=1,NX
            DO J=1,NY
              DO K=1,NZ
                IEL=(K-1)*NX*NY+(J-1)*NX+I
                IALB=IQFR(IQW,I,J,K)
                IF(IALB > 0) THEN
                  IF(IALB.GT.NALB) CALL XABORT('NSSFL5: BETA OVERFLOW.')
                  QFR2(IQW,IEL,IG)=QFR(IQW,IEL)*ALB(BETA(IALB,IG,IG))
                ELSE
                  QFR2(IQW,IEL,IG)=QFR(IQW,IEL)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
*----
*  INITIALIZATIONS
*----
      KEFF_OLD=0.0
      KEFF=1.0
      CALL LCMLEN(IPFLX,'FLUX',ILONG,ITYLCM)
      IF(ILONG == 0) THEN
        JPFLX=LCMLID(IPFLX,'FLUX',NG)
        SAVG(:NUN,:NG)=1.0
      ELSE
        JPFLX=LCMGID(IPFLX,'FLUX')
        DO IG=1,NG
          CALL LCMLEL(JPFLX,IG,ILONG,ITYLCM)
          IF(ILONG /= NUN) CALL XABORT('NSSFL5: INVALID FLUX.')
          CALL LCMGDL(JPFLX,IG,SAVG(1,IG))
        ENDDO
        CALL LCMGET(IPFLX,'K-EFFECTIVE',KEFF)
      ENDIF
      CALL LCMLEN(IPFLX,'DRIFT',ILONG,ITYLCM)
      IF(ILONG == 0) THEN
        JPFLX=LCMLID(IPFLX,'DRIFT',6*NEL)
        DRIFT(:6,:NEL,:NG)=0.0
      ELSE
        JPFLX=LCMGID(IPFLX,'DRIFT')
        DO IG=1,NG
          CALL LCMLEL(JPFLX,IG,ILONG,ITYLCM)
          IF(ILONG /= 6*NEL) CALL XABORT('NSSFL5: INVALID DRIFT.')
          CALL LCMGDL(JPFLX,IG,DRIFT(1,1,IG))
        ENDDO
      ENDIF
      DO IND1=1,LL4F
        DO IG=1,NG
          EVECT((IG-1)*LL4F+IND1)=SAVG(IND1,IG)
        ENDDO
      ENDDO
*----
*  NODAL CORRECTION LOOP
*----
      NMAX=IMAX(LL4F)
      NMAY=IMAY(LL4F)
      NMAZ=IMAZ(LL4F)
      ALLOCATE(A11X(NMAX,NG),A11Y(NMAY,NG),A11Z(NMAZ,NG))
      JTER=0
      SAVG(:NUN,:NG)=0.0
      IOFY=7*LL4F+LL4X
      IOFZ=7*LL4F+LL4X+LL4Y
      DO WHILE((ABS(KEFF_OLD-KEFF) >= EPSNOD).OR.(JTER==0))
        JTER=JTER+1
        IF(IMPX > 0) THEN
          WRITE(6,'(36H NSSFL5: Nodal correction iteration=,I5)')
     >    JTER
        ENDIF
        IF(JTER > MAXNOD) THEN
          WRITE(6,'(/22H ACCURACY AT ITERATION,I4,2H =,1P,E12.5)')
     >    JTER,ABS(KEFF_OLD-KEFF)
          CALL XABORT('NSSFL5: NODAL ITERATION FAILURE')
        ENDIF
        !
        ! set CMFD matrices
        IOF=0
        DO IG=1,NG
          CALL NSSMXYZ(LL4F,NDIM,NX,NY,NZ,NMIX,MAT,XX,YY,ZZ,IDL,VOL,
     >    IQFR,QFR2(1,1,IG),DIFF(1,IG),DRIFT(1,1,IG),SIGR(1,IG),
     >    MUX,MUY,MUZ,IMAX,IMAY,IMAZ,IPY,IPZ,A11X(1,IG),A11Y(1,IG),
     >    A11Z(1,IG))
        ENDDO
        !
        ! CMFD power iteration
        DELTA=ABS(KEFF_OLD-KEFF)
        KEFF_OLD=KEFF
        CALL NSSEIG(NMAX,NMAY,NMAZ,LL4F,NDIM,NEL,NMIX,NG,MAT,IDL,VOL,
     >  MUX,MUY,MUZ,IMAX,IMAY,IMAZ,IPY,IPZ,CHI,SIGF,SCAT,A11X,A11Y,A11Z,
     >  EPSTHR,MAXTHR,NADI,EPSOUT,MAXOUT,ICL1,ICL2,ITER,EVECT,KEFF,IMPX)
        IF(IMPX > 0) WRITE(6,10) JTER,KEFF,ITER,DELTA
        IF(IMPX > 2) THEN
          WRITE(6,'(1X,A)') 'NSSFL5: EVECT='
          IOF=0
          DO IG=1,NG
            WRITE(6,'(1X,1P,14E12.4)') EVECT(IOF+1:IOF+LL4F)
            IOF=IOF+LL4F
          ENDDO
        ENDIF
        !
        ! begin construct SAVG
        IF(NUN /= IOFZ+LL4Z) CALL XABORT('NSSFL5: INVALID NUN.')
        DO IEL=1,LL4F
          DO IG=1,NG
            SAVG(IEL,IG)=EVECT((IG-1)*LL4F+IEL)
          ENDDO
        ENDDO
        !
        ! one- and two-node anm relations
        CALL NSSANM3(NUN,NX,NY,NZ,LL4F,LL4X,LL4Y,NG,BNDTL,NPASS,NMIX,
     >  IDL,KN,IQFR,QFR2,MAT,XXX,YYY,ZZZ,KEFF,DIFF,SIGR,CHI,SIGF,SCAT,
     >  FD,SAVG)
        !
        ! compute new drift coefficients
        DRIFT(:6,:NEL,:NG)=0.0
        DO IG=1,NG
          DO K=1,NZ
            DO J=1,NY
              DO I=1,NX
                IEL=(K-1)*NX*NY+(J-1)*NX+I
                IND1=IDL(IEL)
                IF(IND1 == 0) CYCLE
                KK1=IQFR(1,I,J,K)
                KK2=IQFR(2,I,J,K)
                JXM=KN(1,I,J,K) ; JXP=KN(2,I,J,K)
                JYM=KN(3,I,J,K) ; JYP=KN(4,I,J,K)
                JZM=KN(5,I,J,K) ; JZP=KN(6,I,J,K)
                CALL NSSCO(NX,NY,NZ,NMIX,I,J,K,MAT,XX,YY,ZZ,DIFF(1,IG),
     >          IQFR(1,I,J,K),QFR2(1,IEL,IG),COEF)
                IF((KK1 < 0) .AND. (KK2 < 0)) THEN
                  DRIFT(1,IEL,IG)=-(SAVG(7*LL4F+JXM,IG)+COEF(1)*
     >            SAVG(IND1,IG))/SAVG(IND1,IG)
                  DRIFT(2,IEL,IG)=-(SAVG(7*LL4F+JXP,IG)-COEF(2)*
     >            SAVG(IND1,IG))/SAVG(IND1,IG)
                ELSE IF(KK1 < 0) THEN
                  DRIFT(1,IEL,IG)=-(SAVG(7*LL4F+JXM,IG)+COEF(1)*
     >            SAVG(IND1,IG))/SAVG(IND1,IG)
                  IND3=IDL((K-1)*NX*NY+(J-1)*NX+I+1)
                  IF(IND3 /= 0) DRIFT(2,IEL,IG)=-(SAVG(7*LL4F+JXP,IG)+
     >            COEF(2)*(SAVG(IND3,IG)-SAVG(IND1,IG)))/(SAVG(IND3,IG)+
     >            SAVG(IND1,IG))
                ELSE IF(KK2 < 0) THEN
                  IND2=IDL((K-1)*NX*NY+(J-1)*NX+I-1)
                  IF(IND2 /= 0) DRIFT(1,IEL,IG)=-(SAVG(7*LL4F+JXM,IG)+
     >            COEF(1)*(SAVG(IND1,IG)-SAVG(IND2,IG)))/(SAVG(IND1,IG)+
     >            SAVG(IND2,IG))
                  DRIFT(2,IEL,IG)=-(SAVG(7*LL4F+JXP,IG)-COEF(2)*
     >            SAVG(IND1,IG))/SAVG(IND1,IG)
                ELSE
                  IND2=IDL((K-1)*NX*NY+(J-1)*NX+I-1)
                  IND3=IDL((K-1)*NX*NY+(J-1)*NX+I+1)
                  IF(IND2 /= 0) DRIFT(1,IEL,IG)=-(SAVG(7*LL4F+JXM,IG)+
     >            COEF(1)*(SAVG(IND1,IG)-SAVG(IND2,IG)))/(SAVG(IND1,IG)+
     >            SAVG(IND2,IG))
                  IF(IND3 /= 0) DRIFT(2,IEL,IG)=-(SAVG(7*LL4F+JXP,IG)+
     >            COEF(2)*(SAVG(IND3,IG)-SAVG(IND1,IG)))/(SAVG(IND3,IG)+
     >            SAVG(IND1,IG))
                ENDIF
                KK3=IQFR(3,I,J,K)
                KK4=IQFR(4,I,J,K)
                IF((KK3 < 0).AND.(KK4 < 0)) THEN
                  DRIFT(3,IEL,IG)=-(SAVG(IOFY+JYM,IG)+COEF(3)*
     >            SAVG(IND1,IG))/SAVG(IND1,IG)
                  DRIFT(4,IEL,IG)=-(SAVG(IOFY+JYP,IG)-COEF(4)*
     >            SAVG(IND1,IG))/SAVG(IND1,IG)
                ELSE IF(KK3 < 0) THEN
                  DRIFT(3,IEL,IG)=-(SAVG(IOFY+JYM,IG)+COEF(3)*
     >            SAVG(IND1,IG))/SAVG(IND1,IG)
                  IND3=IDL((K-1)*NX*NY+J*NX+I)
                  IF(IND3 /= 0) DRIFT(4,IEL,IG)=-(SAVG(IOFY+JYP,IG)+
     >            COEF(4)*(SAVG(IND3,IG)-SAVG(IND1,IG)))/(SAVG(IND3,IG)+
     >            SAVG(IND1,IG))
                ELSE IF(KK4 < 0) THEN
                  IND2=IDL((K-1)*NX*NY+(J-2)*NX+I)
                  IF(IND2 /= 0) DRIFT(3,IEL,IG)=-(SAVG(IOFY+JYM,IG)+
     >            COEF(3)*(SAVG(IND1,IG)-SAVG(IND2,IG)))/(SAVG(IND1,IG)+
     >            SAVG(IND2,IG))
                  DRIFT(4,IEL,IG)=-(SAVG(IOFY+JYP,IG)-COEF(4)*
     >            SAVG(IND1,IG))/SAVG(IND1,IG)
                ELSE
                  IND2=IDL((K-1)*NX*NY+(J-2)*NX+I)
                  IND3=IDL((K-1)*NX*NY+J*NX+I)
                  IF(IND2 /= 0) DRIFT(3,IEL,IG)=-(SAVG(IOFY+JYM,IG)+
     >            COEF(3)*(SAVG(IND1,IG)-SAVG(IND2,IG)))/(SAVG(IND1,IG)+
     >            SAVG(IND2,IG))
                  IF(IND3 /= 0) DRIFT(4,IEL,IG)=-(SAVG(IOFY+JYP,IG)+
     >            COEF(4)*(SAVG(IND3,IG)-SAVG(IND1,IG)))/(SAVG(IND3,IG)+
     >            SAVG(IND1,IG))
                ENDIF
                KK5=IQFR(5,I,J,K)
                KK6=IQFR(6,I,J,K)
                IF((KK5 < 0).AND.(KK6 < 0)) THEN
                  DRIFT(5,IEL,IG)=-(SAVG(IOFZ+JZM,IG)+COEF(5)*
     >            SAVG(IND1,IG))/SAVG(IND1,IG)
                  DRIFT(6,IEL,IG)=-(SAVG(IOFZ+JZP,IG)-COEF(6)*
     >            SAVG(IND1,IG))/SAVG(IND1,IG)
                ELSE IF(KK5 < 0) THEN
                  DRIFT(5,IEL,IG)=-(SAVG(IOFZ+JZM,IG)+COEF(5)*
     >            SAVG(IND1,IG))/SAVG(IND1,IG)
                  IND3=IDL(K*NX*NY+(J-1)*NX+I)
                  IF(IND3 /= 0) DRIFT(6,IEL,IG)=-(SAVG(IOFZ+JZP,IG)+
     >            COEF(6)*(SAVG(IND3,IG)-SAVG(IND1,IG)))/(SAVG(IND3,IG)+
     >            SAVG(IND1,IG))
                ELSE IF(KK6 < 0) THEN
                  IND2=IDL((K-2)*NX*NY+(J-1)*NX+I)
                  IF(IND2 /= 0) DRIFT(5,IEL,IG)=-(SAVG(IOFZ+JZM,IG)+
     >            COEF(5)*(SAVG(IND1,IG)-SAVG(IND2,IG)))/(SAVG(IND1,IG)+
     >            SAVG(IND2,IG))
                  DRIFT(6,IEL,IG)=-(SAVG(IOFZ+JZP,IG)-COEF(6)*
     >            SAVG(IND1,IG))/SAVG(IND1,IG)
                ELSE
                  IND2=IDL((K-2)*NX*NY+(J-1)*NX+I)
                  IND3=IDL(K*NX*NY+(J-1)*NX+I)
                  IF(IND2 /= 0) DRIFT(5,IEL,IG)=-(SAVG(IOFZ+JZM,IG)+
     >            COEF(5)*(SAVG(IND1,IG)-SAVG(IND2,IG)))/(SAVG(IND1,IG)+
     >            SAVG(IND2,IG))
                  IF(IND3 /= 0) DRIFT(6,IEL,IG)=-(SAVG(IOFZ+JZP,IG)+
     >            COEF(6)*(SAVG(IND3,IG)-SAVG(IND1,IG)))/(SAVG(IND3,IG)+
     >            SAVG(IND1,IG))
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        IF(IMPX > 5) THEN
          DO IG=1,NG
            WRITE(6,'(28H NSSFL5: DRIFT COEFFICIENTS(,I5,2H):)') IG
            DO IEL=1,NX*NY*NZ
              WRITE(6,'(1P,I7,6E12.4)') IEL,DRIFT(:6,IEL,IG)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      DEALLOCATE(A11Z,A11Y,A11X)
*----
*  END OF NODAL CORRECTION LOOP
*----
      IF(IMPX.GT.0) WRITE(6,20) KEFF,JTER
      IF(IMPX > 2) THEN
        WRITE(6,'(/21H NSSFL5: UNKNOWNS----)')
        DO IG=1,NG
          WRITE(6,'(14H NSSFL5: SAVG(,I4,2H)=)') IG
          WRITE(6,'(1P,12E12.4)') SAVG(:LL4F,IG)
          WRITE(6,'(19H X-BOUNDARY FLUXES:)')
          WRITE(6,'(1P,12E12.4)') SAVG(LL4F+1:2*LL4F,IG)
          WRITE(6,'(1P,12E12.4)') SAVG(2*LL4F+1:3*LL4F,IG)
          WRITE(6,'(19H Y-BOUNDARY FLUXES:)')
          WRITE(6,'(1P,12E12.4)') SAVG(3*LL4F+1:4*LL4F,IG)
          WRITE(6,'(1P,12E12.4)') SAVG(4*LL4F+1:5*LL4F,IG)
          WRITE(6,'(19H Z-BOUNDARY FLUXES:)')
          WRITE(6,'(1P,12E12.4)') SAVG(5*LL4F+1:6*LL4F,IG)
          WRITE(6,'(1P,12E12.4)') SAVG(6*LL4F+1:7*LL4F,IG)
          WRITE(6,'(12H X-CURRENTS:)')
          WRITE(6,'(1P,12E12.4)') SAVG(7*LL4F+1:IOFY,IG)
          WRITE(6,'(12H Y-CURRENTS:)')
          WRITE(6,'(1P,12E12.4)') SAVG(IOFY+1:IOFZ,IG)
          WRITE(6,'(12H Z-CURRENTS:)')
          WRITE(6,'(1P,12E12.4)') SAVG(IOFZ+1:NUN,IG)
          WRITE(6,'(5H ----)')
        ENDDO
      ENDIF
*----
*  SAVE SOLUTION
*----
      JPFLX=LCMGID(IPFLX,'FLUX')
      DO IG=1,NG
        CALL LCMPDL(JPFLX,IG,NUN,2,SAVG(1,IG))
      ENDDO
      JPFLX=LCMGID(IPFLX,'DRIFT')
      DO IG=1,NG
        CALL LCMPDL(JPFLX,IG,6*NEL,2,DRIFT(1,1,IG))
      ENDDO
      CALL LCMPUT(IPFLX,'K-EFFECTIVE',1,2,KEFF)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(EVECT,QFR2)
      DEALLOCATE(SAVG,DRIFT)
      RETURN
*
   10 FORMAT(14H NSSFL5: JTER=,I4,11H CMFD KEFF=,1P E13.6,
     1 12H OBTAINED IN,I4,28H CMFD ITERATIONS WITH ERROR=,
     2 1P,E11.4,1H.)
   20 FORMAT(18H NSSFL5: ANM KEFF=,F11.8,12H OBTAINED IN,I5,
     1 29H NODAL CORRECTION ITERATIONS.)
      END
