*DECK NSSDRV
      SUBROUTINE NSSDRV(IPTRK,IPMAC,IPFLX,ICHX,IDIM,NUN,NG,NEL,NMIX,
     1 ITRIAL,ICL1,ICL2,NADI,EPSNOD,MAXNOD,EPSTHR,MAXTHR,EPSOUT,MAXOUT,
     2 LNODF,BNDTL,NPASS,BB2,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for the flux calculation with the nodal expansion method.
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
* IPTRK   nodal tracking.
* IPMAC   nodal macrolib.
* IPFLX   nodal flux.
* ICHX    solution flag (=4.:CMFD; =5: NEM; =6: ANM).
* IDIM    number of dimensions (1, 2 or 3).
* NUN     number of unknowns per energy group.
* NG      number of energy groups.
* NEL     number of nodes in the nodal calculation.
* NMIX    number of mixtures in the nodal calculation.
* ITRIAL  type of expansion functions in the nodal calculation
*         (=1: polynomial; =2: hyperbolic).
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
* LNODF   flag set to .true. to force discontinuity factors to one.
* BNDTL   set to 'flat', 'linear' or 'quadratic' in 2D cases.
* BB2     imposed leakage used in non-regression tests.
* NPASS   number of transverse current iterations.
* IPRINT  edition flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPMAC,IPFLX
      INTEGER ICHX,IDIM,NUN,NG,NEL,NMIX,ITRIAL(NMIX,NG),ICL1,ICL2,
     > MAXNOD,NADI,MAXTHR,MAXOUT,NPASS,IPRINT
      REAL EPSNOD,EPSTHR,EPSOUT,BB2
      LOGICAL LNODF
      CHARACTER*12 BNDTL
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE),ICODE(6)
      TYPE(C_PTR) JPMAC,KPMAC
      CHARACTER(LEN=8) HADF(6)
      CHARACTER(LEN=72) TITLE
      CHARACTER(LEN=132) HSMG
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAT,IJJ,NJJ,IPOS,IDL,MUX,
     1 MUY,MUZ,IMAX,IMAY,IMAZ,IPY,IPZ
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: KN,IQFR
      REAL, ALLOCATABLE, DIMENSION(:) :: XX,YY,ZZ,XXX,YYY,ZZZ,WORK,VOL
      REAL, ALLOCATABLE, DIMENSION(:,:) :: DIFF,SIGR,CHI,SIGF,QFR,ALBP,
     1 GAR2,GAR3
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: BETA,SCAT,FDXM,FDXP
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: FD
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(MAT(NEL),IDL(NEL),KN(6,NEL),IQFR(6,NEL))
      ALLOCATE(XX(NEL),YY(NEL),ZZ(NEL),VOL(NEL),DIFF(NMIX,NG),
     1 SIGR(NMIX,NG),CHI(NMIX,NG),SIGF(NMIX,NG),SCAT(NMIX,NG,NG),
     2 QFR(6,NEL),FD(NMIX,2*IDIM,NG,NG))
*----
*  RECOVER TRACKING INFORMATION
*----
      TITLE=' '
      CALL LCMLEN(IPTRK,'TITLE',LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
        CALL LCMGTC(IPTRK,'TITLE',72,1,TITLE)
        IF(IPRINT.GT.0) WRITE(6,'(/9H NSSDRV: ,A72)') TITLE
      ENDIF
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NX=ISTATE(14)
      NY=ISTATE(15)
      NZ=ISTATE(16)
      LL4F=ISTATE(25)
      LL4X=ISTATE(27)
      LL4Y=ISTATE(28)
      LL4Z=ISTATE(29)
      ALLOCATE(MUX(LL4F),MUY(LL4F),MUZ(LL4F),IMAX(LL4F),IMAY(LL4F),
     1 IMAZ(LL4F),IPY(LL4F),IPZ(LL4F))
      CALL LCMGET(IPTRK,'ICODE',ICODE)
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'KEYFLX',IDL)
      CALL LCMGET(IPTRK,'VOLUME',VOL)
      CALL LCMGET(IPTRK,'XX',XX)
      CALL LCMGET(IPTRK,'KN',KN)
      IF(IDIM.GE.2) CALL LCMGET(IPTRK,'YY',YY)
      IF(IDIM.EQ.3) CALL LCMGET(IPTRK,'ZZ',ZZ)
      CALL LCMGET(IPTRK,'QFR',QFR)
      CALL LCMGET(IPTRK,'IQFR',IQFR)
      ALLOCATE(XXX(NX+1),YYY(NY+1),ZZZ(NZ+1))
      CALL LCMGET(IPTRK,'XXX',XXX)
      CALL LCMGET(IPTRK,'MUX',MUX)
      CALL LCMGET(IPTRK,'IMAX',IMAX)
      IF(IDIM.GE.2) THEN
        CALL LCMGET(IPTRK,'YYY',YYY)
        CALL LCMGET(IPTRK,'MUY',MUY)
        CALL LCMGET(IPTRK,'IMAY',IMAY)
        CALL LCMGET(IPTRK,'IPY',IPY)
      ENDIF
      IF(IDIM.EQ.3) THEN
        CALL LCMGET(IPTRK,'ZZZ',ZZZ)
        CALL LCMGET(IPTRK,'MUZ',MUZ)
        CALL LCMGET(IPTRK,'IMAZ',IMAZ)
        CALL LCMGET(IPTRK,'IPZ',IPZ)
      ENDIF
*----
*  RECOVER MACROLIB INFORMATION
*----
      IF(BB2.NE.0.0) THEN
        IF(IPRINT.GT.0) WRITE(6,'(/32H NSSDRV: INCLUDE LEAKAGE IN THE ,
     >  13HMACROLIB (B2=,1P,E12.5,2H).)') BB2
      ENDIF
      CALL LCMGET(IPMAC,'STATE-VECTOR',ISTATE)
      NALB=ISTATE(8) ! number of physical albedos
      JPMAC=LCMGID(IPMAC,'GROUP')
      ALLOCATE(WORK(NMIX*NG),IJJ(NMIX),NJJ(NMIX),IPOS(NMIX))
      DO IGR=1,NG
        KPMAC=LCMGIL(JPMAC,IGR)
        CALL LCMGET(KPMAC,'NTOT0',SIGR(1,IGR))
        CALL LCMGET(KPMAC,'DIFF',DIFF(1,IGR))
        CALL LCMGET(KPMAC,'CHI',CHI(1,IGR))
        CALL LCMGET(KPMAC,'NUSIGF',SIGF(1,IGR))
        CALL LCMGET(KPMAC,'IJJS00',IJJ)
        CALL LCMGET(KPMAC,'NJJS00',NJJ)
        CALL LCMGET(KPMAC,'IPOS00',IPOS)
        CALL LCMGET(KPMAC,'SCAT00',WORK)
        DO IBM=1,NMIX
          SCAT(IBM,IGR,:)=0.0
          IPOSDE=IPOS(IBM)-1
          DO JGR=IJJ(IBM),IJJ(IBM)-NJJ(IBM)+1,-1
            IPOSDE=IPOSDE+1
            IF(IPOSDE.GT.NMIX*NG) CALL XABORT('NSSDRV: SCAT OVERFLOW.')
            SCAT(IBM,IGR,JGR)=WORK(IPOSDE) ! IGR <-- JGR
          ENDDO
          SIGR(IBM,IGR)=SIGR(IBM,IGR)-SCAT(IBM,IGR,IGR)
        ENDDO
        IF(BB2.NE.0.0) THEN
          DO IBM=1,NMIX
            SIGR(IBM,IGR)=SIGR(IBM,IGR)+BB2*DIFF(IBM,IGR)
          ENDDO
        ENDIF
        DO IBM=1,NMIX
          IF(SIGR(IBM,IGR).LE.0.0) CALL XABORT('NSSDRV: SIGR<=0.')
        ENDDO
      ENDDO
      DEALLOCATE(IPOS,NJJ,IJJ,WORK)
      ALLOCATE(FDXM(NMIX,NG,NG),FDXP(NMIX,NG,NG),BETA(NALB,NG,NG),
     > GAR2(NMIX,NG),GAR3(NMIX,NG))
      IF(NALB.GT.0) THEN
        CALL LCMLEN(IPMAC,'ALBEDO',ILONG,ITYLCM)
        IF(ILONG.EQ.NALB*NG) THEN
          ALLOCATE(ALBP(NALB,NG))
          CALL LCMGET(IPMAC,'ALBEDO',ALBP)
          BETA(:,:,:)=1.0
          DO IGR=1,NG
            BETA(:NALB,IGR,IGR)=ALBP(:NALB,IGR)
          ENDDO
          DEALLOCATE(ALBP)
        ELSE IF(ILONG.EQ.NALB*NG*NG) THEN
          CALL LCMGET(IPMAC,'ALBEDO',BETA)
        ELSE
          CALL XABORT('NSSDRV: INVALID ALBEDO LENGTH.')
        ENDIF
        IF(IPRINT.GT.1) THEN
          DO IALB=1,NALB
            WRITE(6,'(/35H NSSDRV: PHYSICAL ALBEDO MATRIX ID=,I4)') IALB
            DO IGR=1,NG
              WRITE(6,'(5X,1P,10E12.4)') BETA(IALB,IGR,:)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
      FD(:,:,:,:)=0.0
      IF(LNODF.OR.ISTATE(12).EQ.0) THEN
        DO IBM=1,NMIX
          DO IGR=1,NG
            FD(IBM,1:2,IGR,IGR)=1.0
          ENDDO
          IF(IDIM.GE.2) THEN
            DO IGR=1,NG
              FD(IBM,3:4,IGR,IGR)=1.0
            ENDDO
          ENDIF
          IF(IDIM.EQ.3) THEN
            DO IGR=1,NG
              FD(IBM,5:6,IGR,IGR)=1.0
            ENDDO
          ENDIF
        ENDDO
      ELSE IF(ISTATE(12).EQ.2) then
        ! APOLLO3 case with 4 equisurf values
        CALL LCMSIX(IPMAC,'ADF',1)
        CALL LCMGET(IPMAC,'NTYPE',NSURFD)
        CALL LCMGET(IPMAC,'AVG_FLUX',GAR3)
        CALL LCMGTC(IPMAC,'HADF',8,NSURFD,HADF)
        DO I=1,NSURFD
          IF(HADF(I)(1:3).NE.'FD_') CYCLE
          CALL LCMLEN(IPMAC,HADF(I),ILONG,ITYLCM)
          IF(ILONG.NE.NG*NMIX) THEN
            WRITE(HSMG,'(27HNSSDRV: INVALID LENGTH FOR ,A)') HADF(I)
            CALL XABORT(HSMG)
          ENDIF
          CALL LCMGET(IPMAC,HADF(I),GAR2)
          DO IBM=1,NMIX
            DO IGR=1,NG
              FD(IBM,I,IGR,IGR)=1.0
              IF((GAR2(IBM,IGR).NE.0.0).AND.(HADF(I)(:3).EQ.'FD_')) THEN
                IF(GAR3(IBM,IGR).NE.0.0) THEN
                  FD(IBM,I,IGR,IGR)=GAR2(IBM,IGR)/GAR3(IBM,IGR)
                ELSE
                  CALL XABORT('NSSDRV: ZERO AVERAGED FLUX.')
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        CALL LCMSIX(IPMAC,' ',2)
      ELSE IF(ISTATE(12).EQ.3) THEN
        CALL LCMSIX(IPMAC,'ADF',1)
        CALL LCMGET(IPMAC,'NTYPE',NSURFD)
        CALL LCMGTC(IPMAC,'HADF',8,1,HADF(1))
        CALL LCMGET(IPMAC,HADF(1),GAR2)
        CALL LCMSIX(IPMAC,' ',2)
        DO IBM=1,NMIX
          DO IGR=1,NG
            FD(IBM,1:2,IGR,IGR)=GAR2(IBM,IGR)
          ENDDO
          IF(IDIM.GE.2) THEN
            DO IGR=1,NG
              FD(IBM,3:4,IGR,IGR)=GAR2(IBM,IGR)
            ENDDO
          ENDIF
          IF(IDIM.EQ.3) THEN
            DO IGR=1,NG
              FD(IBM,5:6,IGR,IGR)=GAR2(IBM,IGR)
            ENDDO
          ENDIF
        ENDDO
      ELSE IF(ISTATE(12).EQ.4) THEN
        CALL LCMSIX(IPMAC,'ADF',1)
        CALL LCMGET(IPMAC,'NTYPE',NSURFD)
        IF(IDIM.EQ.1) THEN
          CALL LCMGTC(IPMAC,'HADF',8,2,HADF)
        ELSE IF(IDIM.EQ.2) THEN
          CALL LCMGTC(IPMAC,'HADF',8,4,HADF)
        ELSE IF(IDIM.EQ.3) THEN
          CALL LCMGTC(IPMAC,'HADF',8,6,HADF)
        ENDIF
        CALL LCMGET(IPMAC,HADF(1),FDXM)
        CALL LCMGET(IPMAC,HADF(2),FDXP)
        DO JGR=1,NG
          DO IGR=1,NG
            FD(:NMIX,1,IGR,JGR)=FDXM(:NMIX,IGR,JGR)
            FD(:NMIX,2,IGR,JGR)=FDXP(:NMIX,IGR,JGR)
          ENDDO
        ENDDO
        IF(IDIM.GE.2) THEN
          CALL LCMGET(IPMAC,HADF(3),FDXM)
          CALL LCMGET(IPMAC,HADF(4),FDXP)
          DO JGR=1,NG
            DO IGR=1,NG
              FD(:NMIX,3,IGR,JGR)=FDXM(:NMIX,IGR,JGR)
              FD(:NMIX,4,IGR,JGR)=FDXP(:NMIX,IGR,JGR)
            ENDDO
          ENDDO
        ENDIF
        IF(IDIM.EQ.3) THEN
          CALL LCMGET(IPMAC,HADF(5),FDXM)
          CALL LCMGET(IPMAC,HADF(6),FDXP)
          DO JGR=1,NG
            DO IGR=1,NG
              FD(:NMIX,5,IGR,JGR)=FDXM(:NMIX,IGR,JGR)
              FD(:NMIX,6,IGR,JGR)=FDXP(:NMIX,IGR,JGR)
            ENDDO
          ENDDO
        ENDIF
        CALL LCMSIX(IPMAC,' ',2)
      ELSE
        WRITE(6,'(13H NSSDRV: IDF=,I3)') ISTATE(12)
        CALL XABORT('NSSDRV: FLUX/CURRENT INFORMATION NOT SUPPORTED.')
      ENDIF
      IF(IPRINT.GT.3) THEN
        DO I=1,NSURFD
          WRITE(6,'(/31H NSSDRV: discontinuity factors ,A8)') HADF(I)
          DO IBM=1,NMIX
            DO IGR=1,NG
              WRITE(6,'(4H FD(,2I4,2H)=,1p,12E12.4/(8X,12E12.4))')
     1        IBM,IGR,FD(IBM,:,IGR,IGR)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      DEALLOCATE(GAR3,GAR2,FDXP,FDXM)
*----
*  COMPUTE THE FLUX AND STORE NODAL SOLUTION IN IPFLX
*----
      IF(ICHX.EQ.5) THEN ! NEM
        CALL NSSFL1(IPFLX,NUN,NG,NEL,NMIX,NALB,ITRIAL,EPSOUT,MAXOUT,
     1  MAT,XX,IQFR,QFR,DIFF,SIGR,CHI,SIGF,SCAT,BETA,FD,IPRINT)
      ELSE IF(ICHX.EQ.4) THEN ! CMFD
        CALL NSSFL2(IPFLX,NUN,NG,NEL,NMIX,NALB,EPSOUT,MAXOUT,MAT,XX,
     1  IQFR,QFR,DIFF,SIGR,CHI,SIGF,SCAT,BETA,FD,IPRINT)
      ELSE IF((ICHX.EQ.6).AND.(IDIM.EQ.1)) THEN ! ANM-1D
        CALL NSSFL3(IPFLX,NUN,NG,NEL,NMIX,NALB,EPSNOD,MAXNOD,EPSOUT,
     1  MAXOUT,MAT,XX,XXX,IDL,IQFR,QFR,DIFF,SIGR,CHI,SIGF,SCAT,BETA,
     2  FD,IPRINT)
      ELSE IF((ICHX.EQ.6).AND.(IDIM.EQ.2)) THEN ! ANM-2D
        CALL NSSFL4(IPFLX,NUN,NG,NX,NY,LL4F,LL4X,LL4Y,NMIX,NALB,ICL1,
     1  ICL2,NADI,EPSNOD,MAXNOD,EPSTHR,MAXTHR,EPSOUT,MAXOUT,MAT,XX,YY,
     2  XXX,YYY,IDL,VOL,KN,IQFR,QFR,DIFF,SIGR,CHI,SIGF,SCAT,BETA,FD,
     3  BNDTL,NPASS,MUX,MUY,IMAX,IMAY,IPY,IPRINT)
      ELSE IF((ICHX.EQ.6).AND.(IDIM.EQ.3)) THEN ! ANM-3D
        CALL NSSFL5(IPFLX,NUN,NG,NX,NY,NZ,LL4F,LL4X,LL4Y,LL4Z,NMIX,NALB,
     1  ICL1,ICL2,NADI,EPSNOD,MAXNOD,EPSTHR,MAXTHR,EPSOUT,MAXOUT,MAT,XX,
     2  YY,ZZ,XXX,YYY,ZZZ,IDL,VOL,KN,IQFR,QFR,DIFF,SIGR,CHI,SIGF,SCAT,
     3  BETA,FD,BNDTL,NPASS,MUX,MUY,MUZ,IMAX,IMAY,IMAZ,IPY,IPZ,IPRINT)
      ELSE
        CALL XABORT('NSSDRV: OPTION NOT AVAILABLE.')
      ENDIF
      ISTATE(:)=0
      ISTATE(1)=NG
      ISTATE(2)=NUN
      ISTATE(6)=2
      CALL LCMPUT(IPFLX,'STATE-VECTOR',NSTATE,1,ISTATE)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IPZ,IPY,IMAZ,IMAY,IMAX,MUZ,MUY,MUX)
      DEALLOCATE(BETA)
      DEALLOCATE(ZZZ,YYY,XXX)
      DEALLOCATE(FD,QFR,SCAT,SIGF,CHI,SIGR,DIFF,VOL,ZZ,YY,XX)
      DEALLOCATE(IQFR,KN,IDL,MAT)
      RETURN
      END
