*DECK SNFE1P
      SUBROUTINE SNFE1D(LX,NMAT,IELEM,EELEM,NM,NLF,NSCT,U,W,PL,MAT,
     1 VOL,TOTAL,ESTOPW,NCODE,ZCODE,DELTAE,QEXT,LFIXUP,LSHOOT,
     2 FUNKNO,ISBS,NBS,ISBSM,BS,WX,WE,CST,ISADPT,IBFP,NUN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 1D slab
* geometry. Albedo boundary conditions. Boltzmann-Fokker-Planck (BFP)
* discretization.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert, A. A. Calloo and C. Bienvenue
*
*Parameters: input
* LX    number of regions.
* NMAT    number of material mixtures.
* IELEM   measure of order of the spatial approximation polynomial:
*         =1 constant - default for HODD;
*         =2 linear - default for DG;
*         >3 higher orders.
* EELEM   measure of order of the energy approximation polynomial:
*         =1 constant - default for HODD;
*         =2 linear - default for DG;
*         >3 higher orders.
* NM      total number of moments of the flux in space and
*         energy.
* NLF     number of $\\mu$ levels.
* NSCT    number of Legendre components in the flux:
*         =1: isotropic sources;
*         =2: linearly anisotropic sources.
* U       base points in $\\mu$ of the SN quadrature.
* W       weights of the SN quadrature.
* PL      discrete values of the Legendre polynomials corresponding
*         to the SN quadrature.
* MAT     material mixture index in each region.
* VOL     volumes of each region.
* TOTAL   macroscopic total cross sections.
* ESTOPW  stopping power.
* NCODE   boundary condition indices.
* ZCODE   albedos.
* DELTAE  energy group width in MeV.
* QEXT    Legendre components of the fixed source.
* QEXT0   initial slowing-down angular fluxes.
* LFIXUP  flag to enable negative flux fixup.
* LSHOOT  flag to enable/disable shooting method.
* ISBS    flag to indicate the presence or not of boundary fixed
*         sources.
* NBS     number of boundary fixed sources.
* ISBSM   flag array to indicate the presence or not of boundary fixed
*         source in each unit surface.
* BS      boundary source array with their intensities.
* WX      spatial closure relation weighting factors.
* WE      energy closure relation weighting factors.
* CST     constants for the polynomial approximations.
* ISADPTX flag to enable/disable spatial adaptive flux calculations.
* ISADPTE flag to enable/disable energy adaptive flux calculations.
* IBFP    type of energy proparation relation:
*         =1 Galerkin type;
*         =2 heuristic Przybylski and Ligou type.
* NUN     total number of unknowns in vector FUNKNO   
*
*Parameters: input/output
* FUNKNO  Legendre components of the flux and boundary fluxes.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER LX,NMAT,IELEM,EELEM,NLF,NSCT,MAT(LX),
     1 NCODE(2),ISBS,NBS,ISBSM(2*ISBS,NLF*ISBS),NM,IBFP,NUN
      REAL U(NLF),W(NLF),PL(NSCT,NLF),VOL(LX),TOTAL(0:NMAT),
     1 ESTOPW(0:NMAT,2),ZCODE(2),DELTAE,QEXT(NUN),
     2 FUNKNO(NUN),
     3 BS(NBS*ISBS),WX(IELEM+1),WE(EELEM+1),
     4 CST(MAX(IELEM,EELEM))
      LOGICAL LFIXUP,LSHOOT,ISADPT(2)
*----
*  LOCAL VARIABLES
*----
      REAL BM,BP,TB,WX0(IELEM+1),WE0(EELEM+1)
      DOUBLE PRECISION XNI(EELEM),FEP(IELEM),XNI1(EELEM),XNI2(EELEM),
     1 XNIA(EELEM),XNIB(EELEM),XNIA1(EELEM),XNIA2(EELEM),XNIB1(EELEM),
     2 XNIB2(EELEM)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: Q
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Q2
      PARAMETER(RLOG=1.0E-8)
      LOGICAL ISSHOOT,ISFIX(2)
*----
*  ALLOCATABLE ARRAYS
*----
      ALLOCATE(Q(NM),Q2(NM,NM+1))
*----
*  LENGTH OF FUNKNO COMPONENTS (IN ORDER)
*----
      LFLX=NM*LX*NSCT
      LXNI=EELEM*NLF
      IF(LSHOOT) LXNI=0
      LFEP=IELEM*NLF*LX
*----
*  INNER ITERATION.
*----

      CALL XDRSET(FUNKNO(:LFLX),LFLX,0.0)
      XNI=0.0D0
      WX0=WX
      WE0=WE

      ! SHOOTING METHOD WHEN THERE IS A NON-VACUUM RIGHT
      ! BOUNDARY CONDITION.
      ISSHOOT=(ZCODE(2).NE.0.0).AND.LSHOOT
      IF(ISSHOOT) THEN
        NS=6
      ELSE
        NS=2
      ENDIF

      ! LOOP OVER ALL DIRECTIONS
      DO 200 M0=1,NLF/2

      ! LOOP FOR SHOOTING METHOD
      DO 500 IS=1,NS

      ! CHOOSE DIRECTION
      IF(MOD(IS,2).EQ.0) THEN
        M=NLF-M0+1 ! FORWARD
      ELSE
        M=M0       ! BACKWARD
      ENDIF

      ! SHOOTING METHOD BOUNDARY CONDITIONS.
      IF(ISSHOOT) THEN
        ! 1ST BACKWARD SWEEP
        IF(IS.EQ.1) THEN
          XNI(:EELEM)=0.0D0
          XNI1(:EELEM)=0.0D0
          XNI2(:EELEM)=0.0D0
        ! 1ST FORWARD SWEEP
        ELSEIF(IS.EQ.2) THEN
          XNIA1=0.0D0
          IF(NCODE(1).EQ.4) THEN
          XNIA1(:EELEM)=REAL(XNI(:EELEM))
          XNI(:EELEM)=0.0D0
          ELSE
          XNI(:EELEM)=ZCODE(1)*REAL(XNI(:EELEM))
          ENDIF
        ! 2ND BACKWARD SWEEP
        ELSEIF(IS.EQ.3) THEN
          XNIA2(:EELEM)=0.0D0
          XNIA(:EELEM)=0.0D0
          IF(NCODE(1).EQ.4) THEN
          XNIA2(:EELEM)=REAL(XNI(:EELEM))
          ELSE
          XNIA(:EELEM)=REAL(XNI(:EELEM))
          ENDIF
          XNI(:EELEM)=1.0D0
        ! 2ND FORWARD SWEEP
        ELSEIF(IS.EQ.4) THEN
          IF(NCODE(1).EQ.4) THEN
          XNIB1(:EELEM)=REAL(XNI(:EELEM))
          XNI1(:EELEM)=XNIA1(:EELEM)/(1.0D0+XNIA1(:EELEM)-XNIB1(:EELEM))
          XNI(:EELEM)=1.0D0
          ELSE
          XNI(:EELEM)=ZCODE(1)*REAL(XNI(:EELEM))
          ENDIF
        ! 3RD BACKWARD SWEEP
        ELSEIF(IS.EQ.5) THEN
          IF(NCODE(1).EQ.4) THEN
          XNIB2(:EELEM)=REAL(XNI(:EELEM))
          XNI2(:EELEM)=XNIA2(:EELEM)/(1.0D0+XNIA2(:EELEM)-XNIB2(:EELEM))
          XNI(:EELEM)=XNI1(:EELEM)
          ELSE
          XNIB(:EELEM)=REAL(XNI(:EELEM))
          XNI(:EELEM)=ZCODE(2)*XNIA(:EELEM)/(1.0D0+ZCODE(2)
     1                *(XNIA(:EELEM)-XNIB(:EELEM)))
          ENDIF
        ! 3RD FORWARD SWEEP
        ELSEIF(IS.EQ.6) THEN
          XNI(:EELEM)=ZCODE(1)*XNI(:EELEM)
          IF(NCODE(1).EQ.4) XNI(:EELEM)=XNI2
        ENDIF
      ! NO SHOOTING METHOD BOUNDARY CONDITIONS
      ELSE
        IF(.NOT.LSHOOT) THEN
        IF(U(M).GT.0.0) THEN
          IF(NCODE(1).NE.4) THEN
            DO IEL=1,EELEM
            IOF=(M-1)*EELEM+IEL
            FUNKNO(LFLX+IOF)=FUNKNO(LFLX+LXNI-IOF+1)
            ENDDO
          ENDIF
        ELSE
          IF(NCODE(2).NE.4) THEN
            DO IEL=1,EELEM
            IOF=(M-1)*EELEM+IEL
            FUNKNO(LFLX+IOF)=FUNKNO(LFLX+LXNI-IOF+1)
            ENDDO
          ENDIF
        ENDIF
        XNI(:EELEM)=0.0D0
        ELSE
          IF(IS.EQ.1) THEN
            XNI(:EELEM)=0.0D0
          ELSE
            XNI(:EELEM)=ZCODE(1)*XNI(:EELEM)
          ENDIF
        ENDIF
      ENDIF

      ! BOUNDARY FIXED SOURCES
      IF(U(M).GT.0) THEN
        IF(ISBS.EQ.1.AND.ISBSM(1,M).NE.0) THEN
          XNI(1)=XNI(1)+BS(ISBSM(1,M))
        ENDIF
      ELSE
        IF(ISBS.EQ.1.AND.ISBSM(2,M).NE.0) THEN
          XNI(1)=XNI(1)+BS(ISBSM(2,M))
        ENDIF
      ENDIF
      
      ! SWEEPING OVER ALL VOXELS
      DO 30 I0=1,LX
      I=I0
      IF(U(M).LT.0) I=LX+1-I

      ! X-BOUNDARIES CONDITIONS (NO SHOOTING)
      IF(.NOT.LSHOOT.AND.I0.EQ.1) THEN
        DO IEL=1,EELEM
        IOF=(M-1)*EELEM+IEL
        IF(U(M).GT.0.0) THEN
          XNI=FUNKNO(LFLX+IOF)*ZCODE(1)
        ELSE
          XNI=FUNKNO(LFLX+IOF)*ZCODE(2)
        ENDIF
        ENDDO
      ENDIF
      
      ! DATA
      IBM=MAT(I)
      SIGMA=TOTAL(IBM)
      BM=ESTOPW(IBM,1)/DELTAE
      BP=ESTOPW(IBM,2)/DELTAE

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
      Q(IEL)=0.0
      DO L=1,NSCT
      IOF=(I-1)*NSCT*NM+(L-1)*NM+IEL
      Q(IEL)=Q(IEL)+QEXT(IOF)*PL(L,M)/2.0
      ENDDO
      ENDDO

      ! ENERGY GROUP UPPER BOUNDARY INCIDENT FLUX
      DO IEL=1,IELEM
        IOF=(I-1)*NLF*IELEM+(M-1)*IELEM+IEL
        FEP(IEL)=QEXT(LFLX+LXNI+IOF)
      ENDDO

      ISFIX=.FALSE.
      DO WHILE (.NOT.ALL(ISFIX)) ! LOOP FOR ADAPTIVE CALCULATION

      !FLUX MOMENT COEFFICIENTS MATRIX
      Q2=0.0D0
      DO IX=1,IELEM
      DO JX=1,IELEM
      DO IE=1,EELEM
      DO JE=1,EELEM
        II=EELEM*(IX-1)+IE
        JJ=EELEM*(JX-1)+JE

        ! DIAGONAL TERMS
        IF(II.EQ.JJ) THEN
          Q2(II,JJ)=(SIGMA+CST(IE)**2*WE(JE+1)*BP
     1              +(IE-1)*(BM-BP))*VOL(I)
     2              +CST(IX)**2*WX(JX+1)*ABS(U(M))

        ! UPPER DIAGONAL TERMS
        ELSEIF(II.LT.JJ) THEN
          ! ENERGY TERMS
          IF(IX.EQ.JX) THEN
          IF(MOD(IE+JE,2).EQ.1) THEN 
            Q2(II,JJ)=-CST(IE)*CST(JE)*WE(JE+1)*BP*VOL(I)
          ELSE
            Q2(II,JJ)=CST(IE)*CST(JE)*WE(JE+1)*BP*VOL(I)
          ENDIF
          ! X-SPACE TERMS
          ELSEIF(IE.EQ.JE) THEN
          IF(MOD(IX+JX,2).EQ.1) THEN 
            Q2(II,JJ)=CST(IX)*CST(JX)*WX(JX+1)*U(M)
          ELSE
            Q2(II,JJ)=CST(IX)*CST(JX)*WX(JX+1)*ABS(U(M))
          ENDIF
          ENDIF

        ! UNDER DIAGONAL TERMS
        ELSE
          ! ENERGY TERMS
          IF(IX.EQ.JX) THEN
          IF(MOD(IE+JE,2).EQ.1) THEN 
            Q2(II,JJ)=-CST(IE)*CST(JE)*(WE(JE+1)*BP-BM-BP)*VOL(I)
          ELSE
            Q2(II,JJ)=CST(IE)*CST(JE)*(WE(JE+1)*BP+BM-BP)*VOL(I)
          ENDIF
          ! X-SPACE TERMS
          ELSEIF(IE.EQ.JE) THEN
          IF(MOD(IX+JX,2).EQ.1) THEN 
            Q2(II,JJ)=CST(IX)*CST(JX)*(WX(JX+1)-2.0D0)*U(M)
          ELSE
            Q2(II,JJ)=CST(IX)*CST(JX)*WX(JX+1)*ABS(U(M))
          ENDIF
          ENDIF
        ENDIF 
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      ! FLUX SOURCE VECTOR
      DO IX=1,IELEM
      DO IE=1,EELEM
        II=EELEM*(IX-1)+IE
        Q2(II,NM+1)=Q(II)*VOL(I)
        ! ENERGY TERMS
        IF(MOD(IE,2).EQ.1) THEN
          Q2(II,NM+1)=Q2(II,NM+1)+CST(IE)*(BM-WE(1)*BP)*FEP(IX)*VOL(I)
        ELSE
          Q2(II,NM+1)=Q2(II,NM+1)+CST(IE)*(BM+WE(1)*BP)*FEP(IX)*VOL(I)
        ENDIF
        ! X-SPACE TERMS
        IF(MOD(IX,2).EQ.1) THEN
          Q2(II,NM+1)=Q2(II,NM+1)+CST(IX)*(1-WX(1))*XNI(IE)*ABS(U(M))
        ELSE
          Q2(II,NM+1)=Q2(II,NM+1)-CST(IX)*(1+WX(1))*XNI(IE)*U(M)
        ENDIF
      ENDDO
      ENDDO

      CALL ALSBD(NM,1,Q2,IER,NM)
      IF(IER.NE.0) CALL XABORT('SNFE1D: SINGULAR MATRIX.')

      ! ADAPTIVE CORRECTION OF WEIGHTING PARAMETERS
      IF(ANY(ISADPT)) THEN
        IF(ISADPT(1)) THEN
          CALL SNADPT(EELEM,NM,IELEM,Q2(1:NM:1,NM+1),FEP,
     1    TB,WE,ISFIX(1))
        ELSE
          ISFIX(1)=.TRUE.
        ENDIF
        IF(ISADPT(2)) THEN
          CALL SNADPT(IELEM,NM,EELEM,Q2(1:NM:EELEM,NM+1),XNI,
     1    1.0,WX,ISFIX(2))
        ELSE
          ISFIX(2)=.TRUE.
        ENDIF
      ELSE
        ISFIX=.TRUE.
      ENDIF

      END DO ! END OF ADAPTIVE LOOP
        
      ! CLOSURE RELATIONS
      IF(IELEM.EQ.1.AND.LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
      XNI(:EELEM)=WX(1)*XNI(:EELEM)
      FEP(:IELEM)=WE(1)*FEP(:IELEM) 
      DO IX=1,IELEM
      DO IE=1,EELEM
        II=EELEM*(IX-1)+IE
        ! ENERGY TERMS
        IF(MOD(IE,2).EQ.1) THEN
          FEP(IX)=FEP(IX)+CST(IE)*WE(IE+1)*Q2(II,NM+1)
        ELSE
          FEP(IX)=FEP(IX)-CST(IE)*WE(IE+1)*Q2(II,NM+1)
        ENDIF
        ! X-SPACE TERMS
        IF(MOD(IX,2).EQ.1) THEN
          XNI(IE)=XNI(IE)+CST(IX)*WX(IX+1)*Q2(II,NM+1)
        ELSE
          XNI(IE)=XNI(IE)+CST(IX)*WX(IX+1)*Q2(II,NM+1)*SIGN(1.0,U(M))
        ENDIF
      ENDDO
      ENDDO
      IF(IELEM.EQ.1.AND.LFIXUP.AND.(XNI(1).LE.RLOG)) XNI(1)=0.0
      WX=WX0
      WE=WE0

      IF(ISSHOOT.AND.IS.LT.5) GO TO 30

      ! SAVE ENERGY GROUP LOWER BOUNDARY OUTGOING FLUX
      DO IEL=1,IELEM
        IOF=(I-1)*NLF*IELEM+(M-1)*IELEM+IEL
        FUNKNO(LFLX+LXNI+IOF)=REAL(FEP(IEL))/DELTAE
      ENDDO

      ! SAVE LEGENDRE MOMENT OF THE FLUX
      DO L=1,NSCT
      DO IEL=1,NM
      IOF=(I-1)*NSCT*NM+(L-1)*NM+IEL
      FUNKNO(IOF)=FUNKNO(IOF)+W(M)*REAL(Q2(IEL,NM+1))*PL(L,M)
      ENDDO
      ENDDO

   30 CONTINUE ! END OF X-LOOP

      ! SAVE BOUNDARIES FLUX
      IF(.NOT.LSHOOT) THEN
        DO IEL=1,EELEM
        IOF=(M-1)*EELEM+IEL
        IF(.NOT.LSHOOT) FUNKNO(LFLX+IOF)=REAL(XNI(IEL))
        ENDDO
      ENDIF

  500 CONTINUE ! END OF SHOOTING METHOD LOOP 
  200 CONTINUE ! END OF DIRECTION LOOP

      DEALLOCATE(Q,Q2)
      RETURN
      END
