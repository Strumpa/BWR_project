*DECK NSSFL4
      SUBROUTINE NSSFL4(IPFLX,NUN,NG,NX,NY,LL4F,LL4X,LL4Y,NMIX,NALB,
     > ICL1,ICL2,NADI,EPSNOD,MAXNOD,EPSTHR,MAXTHR,EPSOUT,MAXOUT,MAT,
     > XX,YY,XXX,YYY,IDL,VOL,KN,IQFR,QFR,DIFF,SIGR,CHI,SIGF,SCAT,BETA,
     > FD,BNDTL,NPASS,MUX,MUY,IMAX,IMAY,IPY,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Flux calculation for the analytic nodal method in Cartesian 2D
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
* LL4F    number of nodal flux unknowns.
* LL4X    number of nodal X-directed net currents unknowns.
* LL4Y    number of nodal Y-directed net currents unknowns.
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
* XXX     Cartesian coordinates along the X axis.
* YYY     Cartesian coordinates along the Y axis.
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
* IMAX    X-oriented position of each first non-zero column element.
* IMAY    Y-oriented position of each first non-zero column element.
* IPY     Y-oriented permutation matrices.
* IMPX    edition flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPFLX
      INTEGER NUN,NG,NX,NY,LL4F,LL4X,LL4Y,NMIX,NALB,ICL1,ICL2,NADI,
     1 MAXNOD,MAXTHR,MAXOUT,IMPX,MAT(NX*NY),IDL(NX*NY),KN(6,NX,NY),
     2 NPASS,IQFR(6,NX,NY),MUX(LL4F),MUY(LL4F),IMAX(LL4F),IMAY(LL4F),
     3 IPY(LL4F)
      REAL EPSNOD,EPSTHR,EPSOUT,XX(NX*NY),YY(NX*NY),XXX(NX+1),YYY(NY+1),
     1 VOL(NX*NY),QFR(6,NX*NY),DIFF(NMIX,NG),SIGR(NMIX,NG),CHI(NMIX,NG),
     2 SIGF(NMIX,NG),SCAT(NMIX,NG,NG),BETA(NALB,NG,NG),FD(NMIX,4,NG,NG)
      CHARACTER(LEN=12) :: BNDTL
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPFLX
      INTEGER, PARAMETER :: NZ=1,NDIM=2
      INTEGER :: MUZ(1),IMAZ(1),IPZ(1)
      REAL :: COEF(6),KEFF,KEFF_OLD
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: ZZ,EVECT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: A11X,A11Y,A11Z,SAVG
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: QFR2,DRIFT
*
      ALB(X)=0.5*(1.0-X)/(1.0+X)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      NEL=NX*NY
      N=LL4F*NG
      ALLOCATE(QFR2(6,NEL,NG),ZZ(NEL),EVECT(N))
      ALLOCATE(DRIFT(6,NEL,NG),SAVG(NUN,NG))
*----
*  ALBEDO PROCESSING
*----
      QFR2(:6,:NEL,:NG)=0.0
      DO IG=1,NG
        DO IQW=1,4
          DO I=1,NX
            DO J=1,NY
              IEL=(J-1)*NX+I
              IALB=IQFR(IQW,I,J)
              IF(IALB > 0) THEN
                IF(IALB.GT.NALB) CALL XABORT('NSSFL4: BETA OVERFLOW.')
                QFR2(IQW,IEL,IG)=QFR(IQW,IEL)*ALB(BETA(IALB,IG,IG))
              ELSE
                QFR2(IQW,IEL,IG)=QFR(IQW,IEL)
              ENDIF
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
          IF(ILONG /= NUN) CALL XABORT('NSSFL4: INVALID FLUX.')
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
          IF(ILONG /= 6*NEL) CALL XABORT('NSSFL4: INVALID DRIFT.')
          CALL LCMGDL(JPFLX,IG,DRIFT(1,1,IG))
        ENDDO
      ENDIF
      DO IEL=1,LL4F
        DO IG=1,NG
          EVECT((IG-1)*LL4F+IEL)=SAVG(IEL,IG)
        ENDDO
      ENDDO
*----
*  NODAL CORRECTION LOOP
*----
      NMAX=IMAX(LL4F)
      NMAY=IMAY(LL4F)
      NMAZ=1
      ALLOCATE(A11X(NMAX,NG),A11Y(NMAY,NG),A11Z(NMAZ,NG))
      ZZ(:NEL)=1.0
      MUZ(1)=0
      IMAZ(1)=0
      IPZ(1)=0
      JTER=0
      SAVG(:NUN,:NG)=0.0
      IOFY=5*LL4F+LL4X
      DO WHILE((ABS(KEFF_OLD-KEFF) >= EPSNOD).OR.(JTER==0))
        JTER=JTER+1
        IF(IMPX.GT.0) THEN
          WRITE(6,'(36H NSSFL4: Nodal correction iteration=,I5)')
     >    JTER
        ENDIF
        IF(JTER > MAXNOD) THEN
          WRITE(6,'(/22H ACCURACY AT ITERATION,I4,2H =,1P,E12.5)')
     >    JTER,ABS(KEFF_OLD-KEFF)
          CALL XABORT('NSSFL4: NODAL ITERATION FAILURE')
        ENDIF
        !
        ! set coarse mesh finite difference matrix
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
          WRITE(6,'(1X,A)') 'NSSFL4: EVECT='
          IOF=0
          DO IG=1,NG
            WRITE(6,'(1X,1P,14E12.4)') EVECT(IOF+1:IOF+LL4F)
            IOF=IOF+LL4F
          ENDDO
        ENDIF
        !
        ! begin construct SAVG
        IF(NUN /= IOFY+LL4Y) CALL XABORT('NSSFL4: INVALID NUN.')
        DO IND1=1,LL4F
          DO IG=1,NG
            SAVG(IND1,IG)=EVECT((IG-1)*LL4F+IND1)
          ENDDO
        ENDDO
        !
        ! one- and two-node anm relations
        CALL NSSANM2(NUN,NX,NY,LL4F,LL4X,NG,BNDTL,NPASS,NMIX,IDL,KN,
     >  IQFR,QFR2,MAT,XXX,YYY,KEFF,DIFF,SIGR,CHI,SIGF,SCAT,FD,SAVG)
        !
        ! compute new drift coefficients
        DRIFT(:6,:NEL,:NG)=0.0
        DO IG=1,NG
          DO J=1,NY
            DO I=1,NX
              IEL=(J-1)*NX+I
              IND1=IDL(IEL)
              IF(IND1 == 0) CYCLE
              KK1=IQFR(1,I,J)
              KK2=IQFR(2,I,J)
              JXM=KN(1,I,J) ; JXP=KN(2,I,J)
              JYM=KN(3,I,J) ; JYP=KN(4,I,J)
              CALL NSSCO(NX,NY,NZ,NMIX,I,J,1,MAT,XX,YY,ZZ,DIFF(1,IG),
     >        IQFR(1,I,J),QFR2(1,IEL,IG),COEF)
              IF((KK1 < 0).AND.(KK2 < 0)) THEN
                DRIFT(1,IEL,IG)=-(SAVG(5*LL4F+JXM,IG)+COEF(1)*
     >          SAVG(IND1,IG))/SAVG(IND1,IG)
                DRIFT(2,IEL,IG)=-(SAVG(5*LL4F+JXP,IG)-COEF(2)*
     >          SAVG(IND1,IG))/SAVG(IND1,IG)
              ELSE IF(KK1 < 0) THEN
                DRIFT(1,IEL,IG)=-(SAVG(5*LL4F+JXM,IG)+COEF(1)*
     >          SAVG(IND1,IG))/SAVG(IND1,IG)
                IND3=IDL((J-1)*NX+I+1)
                IF(IND3 /= 0) DRIFT(2,IEL,IG)=-(SAVG(5*LL4F+JXP,IG)+
     >          COEF(2)*(SAVG(IND3,IG)-SAVG(IND1,IG)))/(SAVG(IND3,IG)+
     >          SAVG(IND1,IG))
              ELSE IF(KK2 < 0) THEN
                IND2=IDL((J-1)*NX+I-1)
                IF(IND2 /= 0) DRIFT(1,IEL,IG)=-(SAVG(5*LL4F+JXM,IG)+
     >          COEF(1)*(SAVG(IND1,IG)-SAVG(IND2,IG)))/(SAVG(IND1,IG)+
     >          SAVG(IND2,IG))
                DRIFT(2,IEL,IG)=-(SAVG(5*LL4F+JXP,IG)-COEF(2)*
     >          SAVG(IND1,IG))/SAVG(IND1,IG)
              ELSE
                IND2=IDL((J-1)*NX+I-1)
                IND3=IDL((J-1)*NX+I+1)
                IF(IND2 /= 0) DRIFT(1,IEL,IG)=-(SAVG(5*LL4F+JXM,IG)+
     >          COEF(1)*(SAVG(IND1,IG)-SAVG(IND2,IG)))/(SAVG(IND1,IG)+
     >          SAVG(IND2,IG))
                IF(IND3 /= 0) DRIFT(2,IEL,IG)=-(SAVG(5*LL4F+JXP,IG)+
     >          COEF(2)*(SAVG(IND3,IG)-SAVG(IND1,IG)))/(SAVG(IND3,IG)+
     >          SAVG(IND1,IG))
              ENDIF
              KK3=IQFR(3,I,J)
              KK4=IQFR(4,I,J)
              IF((KK3 < 0).AND.(KK4 < 0)) THEN
                DRIFT(3,IEL,IG)=-(SAVG(IOFY+JYM,IG)+COEF(3)*
     >          SAVG(IND1,IG))/SAVG(IND1,IG)
                DRIFT(4,IEL,IG)=-(SAVG(IOFY+JYP,IG)-COEF(4)*
     >          SAVG(IND1,IG))/SAVG(IND1,IG)
              ELSE IF(KK3 < 0) THEN
                DRIFT(3,IEL,IG)=-(SAVG(IOFY+JYM,IG)+COEF(3)*
     >          SAVG(IND1,IG))/SAVG(IND1,IG)
                IND3=IDL(J*NX+I)
                IF(IND3 /= 0) DRIFT(4,IEL,IG)=-(SAVG(IOFY+JYP,IG)+
     >          COEF(4)*(SAVG(IND3,IG)-SAVG(IND1,IG)))/(SAVG(IND3,IG)+
     >          SAVG(IND1,IG))
              ELSE IF(KK4 < 0) THEN
                IND2=IDL((J-2)*NX+I)
                IF(IND2 /= 0) DRIFT(3,IEL,IG)=-(SAVG(IOFY+JYM,IG)+
     >          COEF(3)*(SAVG(IND1,IG)-SAVG(IND2,IG)))/(SAVG(IND1,IG)+
     >          SAVG(IND2,IG))
                DRIFT(4,IEL,IG)=-(SAVG(IOFY+JYP,IG)-COEF(4)*
     >          SAVG(IND1,IG))/SAVG(IND1,IG)
              ELSE
                IND2=IDL((J-2)*NX+I)
                IND3=IDL(J*NX+I)
                IF(IND2 /= 0) DRIFT(3,IEL,IG)=-(SAVG(IOFY+JYM,IG)+
     >          COEF(3)*(SAVG(IND1,IG)-SAVG(IND2,IG)))/(SAVG(IND1,IG)+
     >          SAVG(IND2,IG))
                IF(IND3 /= 0) DRIFT(4,IEL,IG)=-(SAVG(IOFY+JYP,IG)+
     >          COEF(4)*(SAVG(IND3,IG)-SAVG(IND1,IG)))/(SAVG(IND3,IG)+
     >          SAVG(IND1,IG))
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        IF(IMPX > 5) THEN
          DO IG=1,NG
            WRITE(6,'(28H NSSFL4: DRIFT COEFFICIENTS(,I5,2H):)') IG
            DO IEL=1,NX*NY
              WRITE(6,'(1P,I7,4E12.4)') IEL,DRIFT(:4,IEL,IG)
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
        WRITE(6,'(/21H NSSFL4: UNKNOWNS----)')
        DO IG=1,NG
          WRITE(6,'(14H NSSFL4: SAVG(,I4,2H)=)') IG
          WRITE(6,'(1P,12E12.4)') SAVG(:LL4F,IG)
          WRITE(6,'(19H X-BOUNDARY FLUXES:)')
          WRITE(6,'(1P,12E12.4)') SAVG(LL4F+1:2*LL4F,IG)
          WRITE(6,'(1P,12E12.4)') SAVG(2*LL4F+1:3*LL4F,IG)
          WRITE(6,'(19H Y-BOUNDARY FLUXES:)')
          WRITE(6,'(1P,12E12.4)') SAVG(3*LL4F+1:4*LL4F,IG)
          WRITE(6,'(1P,12E12.4)') SAVG(4*LL4F+1:5*LL4F,IG)
          WRITE(6,'(12H X-CURRENTS:)')
          WRITE(6,'(1P,12E12.4)') SAVG(5*LL4F+1:IOFY,IG)
          WRITE(6,'(12H Y-CURRENTS:)')
          WRITE(6,'(1P,12E12.4)') SAVG(IOFY+1:NUN,IG)
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
      DEALLOCATE(SAVG,DRIFT)
      DEALLOCATE(EVECT,ZZ,QFR2)
      RETURN
*
   10 FORMAT(14H NSSFL4: JTER=,I4,11H CMFD KEFF=,1P E13.6,
     1 12H OBTAINED IN,I4,28H CMFD ITERATIONS WITH ERROR=,
     2 1P,E11.4,1H.)
   20 FORMAT(18H NSSFL4: ANM KEFF=,F11.8,12H OBTAINED IN,I5,
     1 11H ITERATIONS)
      END
