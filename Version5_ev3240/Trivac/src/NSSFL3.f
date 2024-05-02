*DECK NSSFL3
      SUBROUTINE NSSFL3(IPFLX,NUN,NG,NEL,NMIX,NALB,EPSNOD,MAXNOD,
     1 EPSOUT,MAXOUT,MAT,XX,XXX,IDL,IQFR,QFR,DIFF,SIGR,CHI,SIGF,SCAT,
     2 BETA,FD,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Flux calculation for the analytic nodal method in Cartesian 1D
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
* NUN     number of unknowns per energy group (=4*NEL+1).
* NG      number of energy groups.
* NEL     number of nodes in the nodal calculation.
* NMIX    number of mixtures in the nodal calculation.
* NALB    number of physical albedos.
* EPSNOD  nodal correction epsilon.
* MAXNOD  maximum number of nodal correction iterations.
* EPSOUT  convergence epsilon for the power method.
* MAXOUT  maximum number of iterations for the power method.
* MAT     material mixtures.
* XX      mesh spacings.
* XXX     Cartesian coordinates along the X axis.
* IDL     position of averaged fluxes in unknown vector.
* IQFR    boundary condition information.
* QFR     albedo function information.
* DIFF    diffusion coefficients
* SIGR    removal cross sections.
* CHI     fission spectra.
* SIGF    nu times fission cross section.
* SCAT    scattering cross section.
* BETA    albedos.
* FD      discontinuity factors.
* IMPX    edition flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPFLX
      INTEGER NUN,NG,NEL,NMIX,NALB,MAXNOD,MAXOUT,IMPX,MAT(NEL),IDL(NEL),
     1 IQFR(6,NEL)
      REAL EPSNOD,EPSOUT,XX(NEL),XXX(NEL+1),QFR(6,NEL),DIFF(NMIX,NG),
     1 SIGR(NMIX,NG),CHI(NMIX,NG),SIGF(NMIX,NG),SCAT(NMIX,NG,NG),
     2 BETA(NALB,NG,NG),FD(NMIX,2,NG,NG)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPFLX
      INTEGER, PARAMETER :: NY=1,NZ=1,NDIM=1
      REAL :: COEF(6),CODR(6),KEFF,KEFF_OLD
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: YY,ZZ,EVECT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: A,SAVG
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: QFR2,DRIFT
*
      ALB(X)=0.5*(1.0-X)/(1.0+X)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      N=NEL*NG
      ALLOCATE(QFR2(6,NEL,NG),YY(NEL),ZZ(NEL),A(N,2*N),EVECT(N))
      ALLOCATE(DRIFT(6,NEL,NG),SAVG(NUN,NG))
*----
*  ALBEDO PROCESSING
*----
      QFR2(:6,:NEL,:NG)=0.0
      DO IG=1,NG
        DO IQW=1,2
          DO IEL=1,NEL
            IALB=IQFR(IQW,IEL)
            IF(IALB > 0) THEN
              IF(IALB.GT.NALB) CALL XABORT('NSSFL3: BETA OVERFLOW.')
              QFR2(IQW,IEL,IG)=QFR(IQW,IEL)*ALB(BETA(IALB,IG,IG))
            ELSE
              QFR2(IQW,IEL,IG)=QFR(IQW,IEL)
            ENDIF
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
          IF(ILONG /= NUN) CALL XABORT('NSSFL3: INVALID FLUX.')
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
          IF(ILONG /= 6*NEL) CALL XABORT('NSSFL3: INVALID DRIFT.')
          CALL LCMGDL(JPFLX,IG,DRIFT(1,1,IG))
        ENDDO
      ENDIF
      DO IEL=1,NEL
        DO IG=1,NG
          EVECT((IG-1)*NEL+IEL)=SAVG(IEL,IG)
        ENDDO
      ENDDO
*----
*  NODAL CORRECTION LOOP
*----
      YY(:NEL)=1.0
      ZZ(:NEL)=1.0
      JTER=0
      DO WHILE((ABS(KEFF_OLD-KEFF) >= EPSNOD).OR.(JTER==0))
        JTER=JTER+1
        IF(IMPX > 0) THEN
          WRITE(6,'(36H NSSFL3: Nodal correction iteration=,I5)')
     1    JTER
        ENDIF
        IF(JTER > MAXNOD) THEN
          WRITE(6,'(/22H ACCURACY AT ITERATION,I4,2H =,1P,E12.5)')
     1    JTER,ABS(KEFF_OLD-KEFF)
          CALL XABORT('NSSFL3: NODAL ITERATION FAILURE')
        ENDIF
        !
        ! set CMFD matrix for x-directed couplings
        A(:N,:2*N)=0.D0
        IOF=0
        DO IG=1,NG
          DO IEL=1,NEL
            IBM=MAT(IEL)
            IF(IBM <= 0) CYCLE
            KEL=IDL(IEL)
            IF(KEL == 0) CYCLE
            VOL0=XX(IEL)
            CALL NSSCO(NX,NY,NZ,NMIX,IEL,1,1,MAT,XX,YY,ZZ,DIFF(1,IG),
     >      IQFR(1,IEL),QFR2(1,IEL,IG),COEF)
            COEF(1:2)=COEF(1:2)*VOL0/XX(IEL)
            CODR(1:2)=DRIFT(1:2,IEL,IG)*VOL0/XX(IEL)
            KEL2=0
            KK1=IQFR(1,IEL)
            IF(KK1 == -4) THEN
              KEL2=IDL(NX)
            ELSE IF(KK1 == 0) THEN
              KEL2=IDL(IEL-1)
            ENDIF
            IF(KEL2 /= 0) THEN
              A(IOF+KEL,IOF+KEL2)=A(IOF+KEL,IOF+KEL2)-COEF(1)+CODR(1)
            ENDIF
            KEL2=0
            KK2=IQFR(2,IEL)
            IF(KK2 == -4) THEN
              KEL2=IDL(1)
            ELSE IF(KK2 == 0) THEN
              KEL2=IDL(IEL+1)
            ENDIF
            IF(KEL2 /= 0) THEN
              A(IOF+KEL,IOF+KEL2)=A(IOF+KEL,IOF+KEL2)-COEF(2)-CODR(2)
            ENDIF
            A(IOF+KEL,IOF+KEL)=A(IOF+KEL,IOF+KEL)+COEF(1)+CODR(1)+
     >      COEF(2)-CODR(2)
            A(IOF+KEL,IOF+KEL)=A(IOF+KEL,IOF+KEL)+SIGR(IBM,IG)*VOL0
          ENDDO
          JOF=0
          DO JG=1,NG ! IG <-- JG
            DO IEL=1,NEL
              IBM=MAT(IEL)
              IF(IBM <= 0) CYCLE
              KEL=IDL(IEL)
              IF(KEL == 0) CYCLE
              IF(IG /= JG) A(IOF+KEL,JOF+KEL)=-XX(IEL)*SCAT(IBM,IG,JG)
              A(IOF+KEL,N+JOF+KEL)=XX(IEL)*CHI(IBM,IG)*SIGF(IBM,JG)
            ENDDO
            JOF=JOF+NEL
          ENDDO
          IOF=IOF+NEL
        ENDDO
        CALL ALSB(N,N,A,IER,N)
        IF(IER /= 0) CALL XABORT('NSSFL3: SINGULAR MATRIX.')
        !
        ! CMFD power iteration (use double precision)
        DELTA=ABS(KEFF_OLD-KEFF)
        KEFF_OLD=KEFF
        CALL AL1EIG(N,A(1,N+1),EPSOUT,MAXOUT,ITER,EVECT,KEFF,IMPX)
*----
*  FLUX NORMALIZATION
*----
        FMAX=MAXVAL(EVECT(:N))
        EVECT(:N)=EVECT(:N)/FMAX
        IF(IMPX > 0) WRITE(6,10) JTER,KEFF,ITER,DELTA
        IF(IMPX > 2) THEN
          WRITE(6,'(1X,A)') 'NSSFL3: EVECT='
          IOF=0
          DO IG=1,NG
            WRITE(6,'(1X,1P,14E12.4)') EVECT(IOF+1:IOF+NEL)
            IOF=IOF+NEL
          ENDDO
        ENDIF
        !
        ! begin construct SAVG
        IF(NUN /= 4*NEL+1) CALL XABORT('NSSFL3: INVALID NUN.')
        SAVG(:NUN,:NG)=0.0
        DO IEL=1,NEL
          DO IG=1,NG
            SAVG(IEL,IG)=EVECT((IG-1)*NEL+IEL)
          ENDDO
        ENDDO
        !
        ! one- and two-node anm relations
        CALL NSSANM1(NEL,NG,NMIX,IQFR,QFR2,MAT,XXX,KEFF,DIFF,SIGR,CHI,
     1  SIGF,SCAT,FD,SAVG)
        !
        ! compute new drift coefficients
        DO IG=1,NG
          DO IEL=1,NEL
            IBM=MAT(IEL)
            IF(IBM == 0) CYCLE
            CALL NSSCO(NX,NY,NZ,NMIX,IEL,1,1,MAT,XX,YY,ZZ,DIFF(1,IG),
     1      IQFR(1,IEL),QFR2(1,IEL,IG),COEF)
            IF(IEL == 1) THEN
              DRIFT(1,IEL,IG)=-(SAVG(3*NEL+IEL,IG)+COEF(1)*SAVG(IEL,IG))
     1        /SAVG(IEL,IG)
              DRIFT(2,IEL,IG)=-(SAVG(3*NEL+IEL+1,IG)+COEF(2)*
     1        (SAVG(IEL+1,IG)-SAVG(IEL,IG)))/(SAVG(IEL+1,IG)+
     2        SAVG(IEL,IG))
            ELSE IF(IEL < NEL) THEN
              DRIFT(1,IEL,IG)=-(SAVG(3*NEL+IEL,IG)+COEF(1)*(SAVG(IEL,IG)
     1        -SAVG(IEL-1,IG)))/(SAVG(IEL,IG)+SAVG(IEL-1,IG))
              DRIFT(2,IEL,IG)=-(SAVG(3*NEL+IEL+1,IG)+COEF(2)*
     1        (SAVG(IEL+1,IG)-SAVG(IEL,IG)))/(SAVG(IEL+1,IG)+
     2        SAVG(IEL,IG))
            ELSE
              DRIFT(1,IEL,IG)=-(SAVG(3*NEL+IEL,IG)+COEF(1)*(SAVG(IEL,IG)
     1        -SAVG(IEL-1,IG)))/(SAVG(IEL,IG)+SAVG(IEL-1,IG))
              DRIFT(2,IEL,IG)=-(SAVG(3*NEL+IEL+1,IG)-COEF(2)*
     1        SAVG(IEL,IG))/SAVG(IEL,IG)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
*----
*  END OF NODAL CORRECTION LOOP
*----
      IF(IMPX.GT.0) WRITE(6,20) KEFF,JTER
      IF(IMPX > 2) THEN
        WRITE(6,'(/21H NSSFL3: UNKNOWNS----)')
        DO IG=1,NG
          WRITE(6,'(14H NSSFL3: SAVG(,I4,2H)=)') IG
          WRITE(6,'(1P,12E12.4)') SAVG(:NEL,IG)
          WRITE(6,'(19H X-BOUNDARY FLUXES:)')
          WRITE(6,'(1P,12E12.4)') SAVG(NEL+1:2*NEL,IG)
          WRITE(6,'(1P,12E12.4)') SAVG(2*NEL+1:3*NEL,IG)
          WRITE(6,'(12H X-CURRENTS:)')
          WRITE(6,'(1P,12E12.4)') SAVG(3*NEL+1:,IG)
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
      DEALLOCATE(EVECT,A,ZZ,YY,QFR2)
      RETURN
*
   10 FORMAT(14H NSSFL3: JTER=,I4,11H CMFD KEFF=,1P E13.6,
     1 12H OBTAINED IN,I4,28H CMFD ITERATIONS WITH ERROR=,
     2 1P,E11.4,1H.)
   20 FORMAT(18H NSSFL3: ANM KEFF=,F11.8,12H OBTAINED IN,I5,
     1 28H NODAL CORRECTION ITERATIONS)
      END
