*DECK SNFE2D
      SUBROUTINE SNFE2D(NUN,NGEFF,IMPX,INCONV,NGIND,LX,LY,IELEM,
     1 EELEM,NM,NME,NMX,NMY,NMAT,NPQ,NSCT,MAT,VOL,TOTAL,ESTOPW,
     2 NCODE,ZCODE,DELTAE,QEXT,LFIXUP,DU,DE,W,MRM,MRMY,DB,DA,PL,FUNKNO,
     3 ISLG,FLUXC,ISBS,NBS,ISBSM,BS,MAXL,WX,WY,WE,CST,ISADPT,IBFP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 2D Cartesian
* geometry for the HODD method. Energy-angle multithreading. Albedo
* boundary conditions. Boltzmann-Fokker-Planck (BFP) discretization.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert, A. A. Calloo and C. Bienvenue
*
*Parameters: input
* NUN     total number of unknowns in vector FUNKNO.
* NGEFF   number of energy groups processed in parallel.
* IMPX    print flag (equal to zero for no print).
* INCONV  energy group convergence flag (set to .FALSE. if converged).
* NGIND   energy group indices assign to the NGEFF set.
* LX      number of meshes along X axis.
* LY      number of meshes along Y axis.
* IELEM   measure of order of the spatial approximation polynomial:
*         =1 constant - default for HODD;
*         =2 linear - default for DG;
*         >3 higher orders.
* EELEM   measure of order of the energy approximation polynomial:
*         =1 constant - default for HODD;
*         =2 linear - default for DG;
*         >3 higher orders.
* NM      number of moments in space and energy for flux components
* NME     number of moments for energy boundaries components
* NMX     number of moments for X axis boundaries components
* NMY     number of moments for Y axis boundaries components
* NMAT    number of material mixtures.
* NPQ     number of SN directions in four octants (including zero-weight
*         directions).
* NSCT    maximum number of spherical harmonics moments of the flux.
* MAT     material mixture index in each region.
* VOL     volumes of each region.
* TOTAL   macroscopic total cross sections.
* ESTOPW  stopping power.
* NCODE   boundary condition indices.
* ZCODE   albedos.
* DELTAE  energy group width in MeV.
* QEXT    Legendre components of the fixed source.
* LFIXUP  flag to enable negative flux fixup.
* DU      first direction cosines ($\\mu$).
* DE      second direction cosines ($\\eta$).
* W       weights.
* MRM     quadrature index.
* MRMY    quadrature index.
* DB      diamond-scheme parameter.
* DA      diamond-scheme parameter.
* PL      discrete values of the spherical harmonics corresponding
*         to the 2D SN quadrature.
* ISBS    flag to indicate the presence or not of boundary fixed
*         sources.
* NBS     number of boundary fixed sources.
* ISBSM   flag array to indicate the presence or not of boundary fixed
*         source in each unit surface.
* BS      boundary source array with their intensities.
* MAXL    maximum size of boundary source array.
* WX      spatial X axis closure relation weighting factors.
* WY      spatial Y axis closure relation weighting factors.
* WE      energy closure relation weighting factors.
* CST     constants for the polynomial approximations.
* ISADPTX flag to enable/disable adaptive X axis flux calculations.
* ISADPTY flag to enable/disable adaptive Y axis flux calculations.
* ISADPTE flag to enable/disable adaptive energy flux calculations.
* IBFP    type of energy proparation relation:
*         =1 Galerkin type;
*         =2 heuristic Przybylski and Ligou type.
*
*Parameters: input/output
* FUNKNO  Legendre components of the flux and boundary fluxes.
* FLUXC   flux at the cutoff energy.
*
*-----------------------------------------------------------------------
#if defined(_OPENMP)
      USE omp_lib
#endif
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NUN,NGEFF,IMPX,NGIND(NGEFF),LX,LY,IELEM,EELEM,
     1 NM,NME,NMX,NMY,NMAT,NPQ,
     2 NSCT,MAT(LX,LY),NCODE(4),MRM(NPQ),MRMY(NPQ),ISLG(NGEFF),ISBS,
     3 NBS,ISBSM(4*ISBS,NPQ*ISBS,NGEFF*ISBS),MAXL
      LOGICAL INCONV(NGEFF)
      REAL VOL(LX,LY),TOTAL(0:NMAT,NGEFF),ESTOPW(0:NMAT,2,NGEFF),
     1 ZCODE(4),DELTAE(NGEFF),QEXT(NUN,NGEFF),DU(NPQ),DE(NPQ),W(NPQ),
     2 DB(LX,NPQ),DA(LX,LY,NPQ),PL(NSCT,NPQ),FUNKNO(NUN,NGEFF),
     3 FLUXC(LX,LY),BS(MAXL*ISBS,NBS*ISBS),WX(IELEM+1),WY(IELEM+1),
     4 WE(EELEM+1),CST(MAX(IELEM,EELEM))
      LOGICAL LFIXUP,ISADPT(3)
*----
*  LOCAL VARIABLES
*----
      INTEGER NPQD(4),IIND(4)
      REAL BM,BP,TB,WX0(IELEM+1),WY0(IELEM+1),WE0(EELEM+1)
      DOUBLE PRECISION Q(NM),Q2(NM,NM+1),FEP(NME),
     1 XNJ(NMY),V
      PARAMETER(IUNOUT=6,RLOG=1.0E-8,PI=3.141592654)
      LOGICAL ISFIX(3)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDANG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: FLUX,FLUX0
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: FLUX_G,
     1 FLUX0_G
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: XNI
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDANG(NPQ,4))
      ALLOCATE(XNI(NMX,LY),FLUX(NM,NSCT,LX,LY),
     1 FLUX0(NME,NPQ,LX,LY))
      ALLOCATE(FLUX_G(NM,NSCT,LX,LY,NGEFF),
     1 FLUX0_G(NME,NPQ,LX,LY,NGEFF))
*----
*  LENGTH OF FUNKNO COMPONENTS (IN ORDER)
*----
      LFLX=NM*LX*LY*NSCT
      LXNI=NMX*LY*NPQ
      LXNJ=NMY*LX*NPQ
      LFEP=NME*LX*LY*NPQ
*----
*  SET OCTANT SWAPPING ORDER.
*----
      NPQD(:4)=0
      INDANG(:NPQ,:4)=0
      DO M=1,NPQ
        VU=DU(M)
        VE=DE(M)
        IF((VU.GE.0.0).AND.(VE.GE.0.0)) THEN
          IND=1
          JND=4
        ELSE IF((VU.LE.0.0).AND.(VE.GE.0.0)) THEN
          IND=2
          JND=3
        ELSE IF((VU.LE.0.0).AND.(VE.LE.0.0)) THEN
          IND=3
          JND=1
        ELSE
          IND=4
          JND=2
        ENDIF
        IIND(JND)=IND
        NPQD(IND)=NPQD(IND)+1
        INDANG(NPQD(IND),IND)=M
      ENDDO
*----
*  MAIN LOOP OVER OCTANTS.
*----

      FLUX_G(:NM,:NSCT,:LX,:LY,:NGEFF)=0.0D0
      FLUX0_G(:NME,:NPQ,:LX,:LY,:NGEFF)=0.0D0
      WE0=WE
      WX0=WX
      WY0=WY
      
      DO 190 JND=1,4
      IND=IIND(JND)
*----
*  PRELIMINARY LOOPS FOR SETTING BOUNDARY CONDITIONS.
*----

*$OMP  PARALLEL DO
*$OMP+ PRIVATE(M,IG,VU,VE,M1,IOF,JOF,IEL,I,J,IPQD)
*$OMP+ SHARED(FUNKNO) COLLAPSE(2)

      DO 70 IG=1,NGEFF
      DO 60 IPQD=1,NPQD(IND)
      IF(.NOT.INCONV(IG)) GO TO 60
      M=INDANG(IPQD,IND)
      VU=DU(M)
      VE=DE(M)
      ! X-BOUNDARY
      IF(VU.GT.0.0)THEN
         M1=MRM(M)
         IF((NCODE(1).NE.4))THEN
            DO IEL=1,NMX
               DO J=1,LY
                  IOF=((M-1)*LY+(J-1))*NMX+IEL
                  JOF=((M1-1)*LY+(J-1))*NMX+IEL
                  FUNKNO(LFLX+IOF,IG)=FUNKNO(LFLX+JOF,IG)
               ENDDO
            ENDDO
         ENDIF
      ELSEIF(VU.LT.0.0)THEN
         M1=MRM(M)
         IF((NCODE(2).NE.4))THEN
            DO IEL=1,NMX
               DO J=1,LY
                  IOF=((M-1)*LY+(J-1))*NMX+IEL
                  JOF=((M1-1)*LY+(J-1))*NMX+IEL
                  FUNKNO(LFLX+IOF,IG)=FUNKNO(LFLX+JOF,IG)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
      ! Y-BOUNDARY
      IF(VE.GT.0.0)THEN
         M1=MRMY(M)
         IF((NCODE(3).NE.4))THEN
            DO IEL=1,NMY
               DO I=1,LX
                  IOF=((M-1)*LX+(I-1))*NMY+IEL
                  JOF=((M1-1)*LX+(I-1))*NMY+IEL
                  FUNKNO(LFLX+LXNI+IOF,IG)=
     >               FUNKNO(LFLX+LXNI+JOF,IG)
               ENDDO
            ENDDO
         ENDIF
      ELSEIF(VE.LT.0.0)THEN
         M1=MRMY(M)
         IF((NCODE(4).NE.4))THEN
            DO IEL=1,NMY
               DO I=1,LX
                  IOF=((M-1)*LX+(I-1))*NMY+IEL
                  JOF=((M1-1)*LX+(I-1))*NMY+IEL
                  FUNKNO(LFLX+LXNI+IOF,IG)=
     >               FUNKNO(LFLX+LXNI+JOF,IG)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
   60 CONTINUE
   70 CONTINUE

*$OMP END PARALLEL DO

*----
*  MAIN SWAPPING LOOPS FOR SN FLUX CALCULATION
*----

*$OMP  PARALLEL DO
*$OMP+ PRIVATE(ITID,FLUX,M,IG,XNI,XNJ,Q,Q2,IOF,IER,II,JJ,IEL,I,J,L)
*$OMP+ PRIVATE(FEP,FLUX0,BM,BP,IIE,IIX,IIY,IE,IX,IY)
*$OMP+ PRIVATE(ISFIX,JX,JY,JE,TB,V,SIGMA,IBM,J0,I0,IPQD)
*$OMP+ FIRSTPRIVATE(WE,WX,WY,WE0,WX0,WY0) SHARED(FUNKNO) 
*$OMP+ REDUCTION(+:FLUX_G,FLUX0_G,FLUXC) COLLAPSE(2)

      ! LOOP FOR GROUPS TO EXECUTE IN PARALLEL
      DO 180 IG=1,NGEFF

      ! LOOP OVER ALL DIRECTIONS
      DO 170 IPQD=1,NPQD(IND)
      IF(.NOT.INCONV(IG)) GO TO 170
      M=INDANG(IPQD,IND)
      IF(W(M).EQ.0.0) GO TO 170

      ! GET AND PRINT THREAD NUMBER
#if defined(_OPENMP)
      ITID=omp_get_thread_num()
#else
      ITID=0
#endif
      IF(IMPX.GT.5) WRITE(IUNOUT,400) ITID,NGIND(IG),IPQD

      ! INITIALIZE FLUXES
      FLUX(:NM,:NSCT,:LX,:LY)=0.0D0
      FLUX0(:NME,:NPQ,:LX,:LY)=0.0D0

*----
*  LOOP OVER X- AND Y-DIRECTED AXES.
*----

      ! X-AXIS LOOP
      DO 155 I0=1,LX
      I=I0
      IF((IND.EQ.2).OR.(IND.EQ.3)) I=LX+1-I
 
      ! Y-BOUNDARIES CONDITIONS
      XNJ=0.0
      DO IEL=1,NMY
      IOF=(M-1)*NMY*LX+(I-1)*NMY+IEL
      IF((IND.EQ.1).OR.(IND.EQ.2)) THEN
         XNJ(IEL)=FUNKNO(LFLX+LXNI+IOF,IG)*ZCODE(3)
      ELSE
         XNJ(IEL)=FUNKNO(LFLX+LXNI+IOF,IG)*ZCODE(4)
      ENDIF
      ENDDO

      ! Y-BOUNDARIES FIXED SOURCES
      IF(ISBS.EQ.1) THEN
        IF((IND.EQ.3.OR.IND.EQ.4).AND.ISBSM(4,M,IG).NE.0) THEN
          XNJ(1)=XNJ(1)+BS(I,ISBSM(4,M,IG))
        ELSE IF((IND.EQ.1.OR.IND.EQ.2).AND.ISBSM(3,M,IG).NE.0) THEN
          XNJ(1)=XNJ(1)+BS(I,ISBSM(3,M,IG))
        ENDIF
      ENDIF

      ! Y-AXIS LOOP
      DO 140 J0=1,LY
      J=J0
      IF((IND.EQ.3).OR.(IND.EQ.4)) J=LY+1-J

      ! X-BOUNDARIES CONDITIONS
      IF(I0.EQ.1) THEN
      XNI(:NMX,J)=0.0
      DO IEL=1,NMX
         IOF=(M-1)*NMX*LY+(J-1)*NMX+IEL
         IF((IND.EQ.1).OR.(IND.EQ.4)) THEN
            XNI(IEL,J)=FUNKNO(LFLX+IOF,IG)*ZCODE(1)
         ELSE
            XNI(IEL,J)=FUNKNO(LFLX+IOF,IG)*ZCODE(2)
         ENDIF
      ENDDO
      ENDIF

      ! X-BOUNDARIES FIXED SOURCES
      IF(ISBS.EQ.1.AND.I0.EQ.1) THEN
        IF((IND.EQ.2.OR.IND.EQ.3).AND.ISBSM(2,M,IG).NE.0) THEN
          XNI(1,J)=XNI(1,J)+BS(J,ISBSM(2,M,IG))
        ELSE IF((IND.EQ.1.OR.IND.EQ.4).AND.ISBSM(1,M,IG).NE.0) THEN
          XNI(1,J)=XNI(1,J)+BS(J,ISBSM(1,M,IG))
        ENDIF
      ENDIF 

      ! DATA
      IBM=MAT(I,J)
      IF(IBM.EQ.0) GO TO 140
      SIGMA=TOTAL(IBM,IG)
      BM=ESTOPW(IBM,1,IG)/DELTAE(IG)
      BP=ESTOPW(IBM,2,IG)/DELTAE(IG)
      V=VOL(I,J)

      ! TYPE OF ENERGY PROPAGATION FACTOR
      IF(IBFP.EQ.1) THEN ! GALERKIN TYPE
        TB=BM/BP
        WE(1)=WE(1)*TB
        WE(2:EELEM+1)=(WE(2:EELEM+1)-1)*TB+1
      ELSE ! PRZYBYLSKI AND LIGOU TYPE
        TB=1.0
      ENDIF

      ! SOURCE DENSITY TERM
      DO IEL=1,NM
      Q(IEL)=0.0D0
      DO L=1,NSCT
      IOF=((J-1)*LX*NSCT+(I-1)*NSCT+(L-1))*NM+IEL
      Q(IEL)=Q(IEL)+QEXT(IOF,IG)*PL(L,M)/(4.0D0*PI)
      ENDDO
      ENDDO

      ! ENERGY GROUP UPPER BOUNDARY INCIDENT FLUX
      DO IEL=1,NME
      IOF=((J-1)*LX*NPQ+(I-1)*NPQ+(M-1))*NME+IEL
      FEP(IEL)=QEXT(LFLX+LXNI+LXNJ+IOF,IG)
      ENDDO

      ISFIX=.FALSE.
      DO WHILE (.NOT.ALL(ISFIX)) ! LOOP FOR ADAPTIVE CALCULATION

      ! FLUX MOMENT COEFFICIENTS MATRIX
      CALL XDDSET(Q2,NM*(NM+1),0.0D0)
      
      DO IY=1,IELEM
      DO JY=1,IELEM
      DO IX=1,IELEM
      DO JX=1,IELEM
      DO IE=1,EELEM
      DO JE=1,EELEM
        II=IELEM*EELEM*(IY-1)+EELEM*(IX-1)+IE
        JJ=IELEM*EELEM*(JY-1)+EELEM*(JX-1)+JE

        ! DIAGONAL TERMS
        IF(II.EQ.JJ) THEN
          Q2(II,JJ)=(SIGMA+CST(IE)**2*WE(JE+1)*BP+(IE-1)*(BM-BP))*V
     1              +CST(IX)**2*WX(JX+1)*ABS(DA(I,J,M))
     2              +CST(IY)**2*WY(JY+1)*ABS(DB(I,M))

        ! UPPER DIAGONAL TERMS
        ELSEIF(II.LT.JJ) THEN
          IF(IY.EQ.JY) THEN
          ! ENERGY TERMS
          IF(IX.EQ.JX) THEN
          IF(MOD(IE+JE,2).EQ.1) THEN
            Q2(II,JJ)=-CST(IE)*CST(JE)*WE(JE+1)*BP*V
          ELSE
            Q2(II,JJ)=CST(IE)*CST(JE)*WE(JE+1)*BP*V
          ENDIF
          ! X-SPACE TERMS
          ELSEIF(IE.EQ.JE) THEN
          IF(MOD(IX+JX,2).EQ.1) THEN
            Q2(II,JJ)=CST(IX)*CST(JX)*WX(JX+1)*DA(I,J,M)
          ELSE
            Q2(II,JJ)=CST(IX)*CST(JX)*WX(JX+1)*ABS(DA(I,J,M))
          ENDIF
          ENDIF
          ! Y-SPACE TERMS
          ELSEIF(IX.EQ.JX.AND.IE.EQ.JE) THEN
          IF(MOD(IY+JY,2).EQ.1) THEN
            Q2(II,JJ)=CST(IY)*CST(JY)*WY(JY+1)*DB(I,M)
          ELSE
            Q2(II,JJ)=CST(IY)*CST(JY)*WY(JY+1)*ABS(DB(I,M))
          ENDIF
          ENDIF 

        ! UNDER DIAGONAL TERMS
        ELSE

          IF(IY.EQ.JY) THEN
          ! ENERGY TERMS
          IF(IX.EQ.JX) THEN
          IF(MOD(IE+JE,2).EQ.1) THEN
            Q2(II,JJ)=-CST(IE)*CST(JE)*(WE(JE+1)*BP-BM-BP)*V
          ELSE
            Q2(II,JJ)=CST(IE)*CST(JE)*(WE(JE+1)*BP+BM-BP)*V
          ENDIF
          ! X-SPACE TERMS
          ELSEIF(IE.EQ.JE) THEN
          IF(MOD(IX+JX,2).EQ.1) THEN
            Q2(II,JJ)=CST(IX)*CST(JX)*(WX(JX+1)-2.0D0)*DA(I,J,M)
          ELSE
            Q2(II,JJ)=CST(IX)*CST(JX)*WX(JX+1)*ABS(DA(I,J,M))
          ENDIF
          ENDIF
          ! Y-SPACE TERMS
          ELSEIF(IX.EQ.JX.AND.IE.EQ.JE) THEN
          IF(MOD(IY+JY,2).EQ.1) THEN
            Q2(II,JJ)=CST(IY)*CST(JY)*(WY(JY+1)-2.0D0)*DB(I,M)
          ELSE
            Q2(II,JJ)=CST(IY)*CST(JY)*WY(JY+1)*ABS(DB(I,M))
          ENDIF
          ENDIF
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      ! FLUX SOURCE VECTOR
      DO IY=1,IELEM
      DO IX=1,IELEM
      DO IE=1,EELEM
        II=IELEM*EELEM*(IY-1)+EELEM*(IX-1)+IE
        IIE=IELEM*(IY-1)+IX
        IIX=EELEM*(IY-1)+IE
        IIY=EELEM*(IX-1)+IE
        Q2(II,NM+1)=Q(II)*V
        ! ENERGY TERMS
        IF(MOD(IE,2).EQ.1) THEN
          Q2(II,NM+1)=Q2(II,NM+1)+CST(IE)*(BM-WE(1)*BP)*FEP(IIE)*V
        ELSE
          Q2(II,NM+1)=Q2(II,NM+1)+CST(IE)*(BM+WE(1)*BP)*FEP(IIE)*V
        ENDIF
        ! X-SPACE TERMS
        IF(MOD(IX,2).EQ.1) THEN
          Q2(II,NM+1)=Q2(II,NM+1)+CST(IX)*(1-WX(1))
     1                *XNI(IIX,J)*ABS(DA(I,J,M))
        ELSE
          Q2(II,NM+1)=Q2(II,NM+1)-CST(IX)*(1+WX(1))
     1                *XNI(IIX,J)*DA(I,J,M)
        ENDIF
        ! Y-SPACE TERMS
        IF(MOD(IY,2).EQ.1) THEN
          Q2(II,NM+1)=Q2(II,NM+1)+CST(IY)*(1-WY(1))
     1                *XNJ(IIY)*ABS(DB(I,M))
        ELSE
          Q2(II,NM+1)=Q2(II,NM+1)-CST(IY)*(1+WY(1))
     1                *XNJ(IIY)*DB(I,M)
        ENDIF
      ENDDO
      ENDDO
      ENDDO

      CALL ALSBD(NM,1,Q2,IER,NM)
      IF(IER.NE.0) CALL XABORT('SNFE2D: SINGULAR MATRIX.')

      ! ADAPTIVE CORRECTION OF WEIGHTING PARAMETERS
      IF(ANY(ISADPT)) THEN
        IF(ISADPT(1)) THEN
          CALL SNADPT(EELEM,NM,IELEM**2,Q2(1:EELEM:1,NM+1),
     1    FEP,TB,WE,ISFIX(1))
        ELSE
          ISFIX(1)=.TRUE.
        ENDIF
        IF(ISADPT(2)) THEN
          CALL SNADPT(IELEM,NM,EELEM*IELEM,Q2(1:IELEM*EELEM:IELEM,NM+1),
     1    XNI(:NMX,J),1.0,WX,ISFIX(2))
        ELSE
          ISFIX(2)=.TRUE.
        ENDIF
        IF(ISADPT(3)) THEN
          CALL SNADPT(IELEM,NM,EELEM*IELEM,Q2(1:NM:IELEM*EELEM,NM+1),
     1    XNJ,1.0,WY,ISFIX(3))
        ELSE
          ISFIX(3)=.TRUE.
        ENDIF
      ELSE
        ISFIX=.TRUE.
      ENDIF

      END DO ! END OF ADAPTIVE LOOP

      ! CLOSURE RELATIONS
      IF(IELEM.EQ.1.AND.LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
      XNI(:NMX,J)=WX(1)*XNI(:NMX,J)
      XNJ(:NMY)=WY(1)*XNJ(:NMY)
      FEP(:NME)=WE(1)*FEP(:NME)
      DO IY=1,IELEM
      DO IX=1,IELEM
      DO IE=1,EELEM
        II=IELEM*EELEM*(IY-1)+EELEM*(IX-1)+IE
        IIE=IELEM*(IY-1)+IX
        IIX=EELEM*(IY-1)+IE
        IIY=EELEM*(IX-1)+IE
        ! ENERGY 
        IF(MOD(IE,2).EQ.1) THEN
          FEP(IIE)=FEP(IIE)+CST(IE)*WE(IE+1)*Q2(II,NM+1)
        ELSE
          FEP(IIE)=FEP(IIE)-CST(IE)*WE(IE+1)*Q2(II,NM+1)
        ENDIF
        ! X-SPACE
        IF(MOD(IX,2).EQ.1) THEN
          XNI(IIX,J)=XNI(IIX,J)+CST(IX)*WX(IX+1)
     1                 *Q2(II,NM+1)
        ELSE
          XNI(IIX,J)=XNI(IIX,J)+CST(IX)*WX(IX+1)
     1                 *Q2(II,NM+1)*SIGN(1.0,DA(I,J,M))
        ENDIF
        ! Y-SPACE
        IF(MOD(IY,2).EQ.1) THEN
          XNJ(IIY)=XNJ(IIY)+CST(IY)*WY(IY+1)
     1                 *Q2(II,NM+1)
        ELSE
          XNJ(IIY)=XNJ(IIY)+CST(IY)*WY(IY+1)
     1                 *Q2(II,NM+1)*SIGN(1.0,DB(I,M))
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      IF(IELEM.EQ.1.AND.LFIXUP.AND.(XNI(1,J).LE.RLOG)) XNI(1,J)=0.0
      IF(IELEM.EQ.1.AND.LFIXUP.AND.(XNJ(1).LE.RLOG)) XNJ(1)=0.0
      WE=WE0
      WX=WX0
      WY=WY0

      ! SAVE ENERGY GROUP LOWER BOUNDARY OUTGOING FLUX 
      FLUX0(:NME,M,I,J)=REAL(FEP(:NME))/DELTAE(IG)

      ! SAVE LAST GROUP LOWER BOUNDARY FLUX
      IF(ISLG(IG).EQ.1) THEN
      FLUXC(I,J)=FLUXC(I,J)+2.0*W(M)*REAL(FLUX0(1,M,I,J))*PL(1,M)
      ENDIF

      ! SAVE LEGENDRE MOMENT OF THE FLUX
      DO L=1,NSCT
      DO IEL=1,NM
      FLUX(IEL,L,I,J)=FLUX(IEL,L,I,J)+2.0D0*W(M)
     1 *Q2(IEL,NM+1)*PL(L,M)
      ENDDO
      ENDDO

  140 CONTINUE ! END OF Y-LOOP

      ! SAVE Y-BOUNDARY CONDITIONS
      DO IEL=1,NMY
      IOF=(M-1)*NMY*LX+(I-1)*NMY+IEL
      FUNKNO(LFLX+LXNI+IOF,IG)=REAL(XNJ(IEL))
      ENDDO

  155 CONTINUE ! END OF X-LOOP

      ! SAVE BOUNDARY CONDITIONS
      DO J=1,LY
      DO IEL=1,NMX
      IOF=(M-1)*NMX*LY+(J-1)*NMX+IEL
      FUNKNO(LFLX+IOF,IG)=REAL(XNI(IEL,J))
      ENDDO
      ENDDO

      ! SAVE FLUX INFORMATION
      FLUX_G(:,:,:,:,IG)=FLUX_G(:,:,:,:,IG)+FLUX(:,:,:,:)
      FLUX0_G(:,:,:,:,IG)=FLUX0_G(:,:,:,:,IG)+FLUX0(:,:,:,:)

  170 CONTINUE ! END OF DIRECTION LOOP
  180 CONTINUE ! END OF ENERGY LOOP
*$OMP END PARALLEL DO
  190 CONTINUE ! END OF OCTANT LOOP
      ! SAVE FLUX INFORMATION
      DO 200 IG=1,NGEFF
        IF(.NOT.INCONV(IG)) GO TO 200
        FUNKNO(:LFLX,IG)=
     1  RESHAPE(REAL(FLUX_G(:NM,:NSCT,:LX,:LY,IG)),
     2   (/ LFLX /) )
        FUNKNO(LFLX+LXNI+LXNJ+1:LFLX+LXNI+LXNJ+LFEP,IG)=
     1  RESHAPE(REAL(FLUX0_G(:NME,:NPQ,:LX,:LY,IG)), (/ LFEP /) )
  200 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XNI,FLUX0_G,FLUX_G,FLUX0,FLUX,INDANG)
      RETURN
  400 FORMAT(16H SNFP12: thread=,I8,12H --->(group=,I4,7H angle=,I4,1H))
      END
