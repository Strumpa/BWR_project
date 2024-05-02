*DECK MCCGF
      SUBROUTINE MCCGF(KPSYS,IPTRK,IFTRAK,IPMACR,IMPX,NGRP,NGEFF,NGIND,
     1                 IDIR,NBREG,NBMIX,NUNKNO,LEXAC,MAT,VOL,KEYFLX,
     2                 FUNKNO2,SUNKNO,TITR,REBFLG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solve N-group transport equation for fluxes using the method of
* characteristics (vectorial version).
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert and R. Le Tellier
*
*Parameters: input/output
* KPSYS   pointer to the assembly LCM object (L_PIJ signature). KPSYS is
*         an array of directories.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IPMACR  pointer to the macrolib LCM object.
* IFTRAK  tracking file unit number.
* IMPX    print flag (equal to zero for no print).
* NGRP    number of energy groups.
* NGEFF   number of energy groups processed in parallel.
* NGIND   energy group indices assign to the NGEFF set.
* IDIR    direction of fundamental current for TIBERE with MoC 
*         (=0,1,2,3). 
* NBREG   total number of volumes for which specific values of the
*         neutron flux and reactions rates are required.
* NBMIX   number of mixtures (NBMIX=max(MAT(i))).
* NUNKNO  total number of unknowns in vectors SUNKNO and FUNKNO.
* LEXAC   type of exponential function calculation (=.false. to compute
*         exponential functions using tables).
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  position of flux elements in FUNKNO vector.
* FUNKNO  unknown vector.
* SUNKNO  input source vector.
* TITR    title.
* REBFLG  ACA or SCR rebalancing flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) KPSYS(NGEFF),IPTRK,IPMACR
      INTEGER NGRP,NGEFF,NGIND(NGEFF),IFTRAK,IMPX,IDIR,NBREG,NBMIX,
     1 NUNKNO,MAT(NBREG),KEYFLX(NBREG)
      REAL VOL(NBREG),FUNKNO2(NUNKNO,NGEFF),SUNKNO(NUNKNO,NGRP)
      CHARACTER TITR*72
      LOGICAL LEXAC,REBFLG
*----
*  GENERIC INTERFACES
*----
      INTERFACE
        SUBROUTINE MOCFFI_TEMPLATE(SUBSCH,NR,NS,NUN,MT,LINE,SEGLEN,
     1             NRSEG,NE,MATALB,SIGANG,KEYFLX,YG,FLUX,EXPT,EXP2,
     2             FLM,FLP,CYM,CYP,IDIR,OMG2)
          INTEGER NR,NS,NUN,MT,LINE,NRSEG(LINE),NE,MATALB(-NS:NR),
     1    KEYFLX(NR),IDIR
          REAL SIGANG(-6:MT),YG(NE)
          DOUBLE PRECISION SEGLEN(LINE),FLUX(NUN),EXPT(NE,LINE),
     1     EXP2(NE,LINE),FLM(NE,LINE),FLP(NE,LINE),CYM(NE,LINE),
     2     CYP(NE,LINE),OMG2(NE,3)
          EXTERNAL SUBSCH
        END SUBROUTINE MOCFFI_TEMPLATE
        !
        SUBROUTINE MOCFFA_TEMPLATE(SUBSCH,NR,NS,NUN,MT,LINE,SEGLEN,
     1             NRSEG,NE,NF,MATALB,SIGANG,KEYFLX,YG,FLUX,EXPT,EXP2,
     2             FLM,FLP,CYM,CYP,NPHI,NSUB,KANGL,TRHAR)
          INTEGER NR,NS,NUN,MT,LINE,NRSEG(LINE),NE,NF,MATALB(-NS:NR),
     1    KEYFLX(NR,NF),NPHI,NSUB,KANGL(NSUB)
          REAL SIGANG(-6:MT),YG(NE),TRHAR(NE,NF,NPHI,2)
          DOUBLE PRECISION SEGLEN(LINE),FLUX(NUN),EXPT(NE,LINE),
     1    EXP2(NE,LINE),FLM(NE,LINE),FLP(NE,LINE),CYM(NE,LINE),
     2    CYP(NE,LINE)
          EXTERNAL SUBSCH
        END SUBROUTINE MOCFFA_TEMPLATE
        !
        SUBROUTINE MOCSCH_TEMPLATE(N,NREG,NSOUT,M,NOM,NZON,H,SIGANG,
     1             EXPT,EXP2,NMU,ZMU)
          INTEGER N,NREG,NSOUT,M,NOM(N),NZON(-NSOUT:NREG),NMU
          REAL SIGANG(-6:M),ZMU(NMU)
          DOUBLE PRECISION H(N),EXPT(NMU,N),EXP2(2,NMU,N)
        END SUBROUTINE MOCSCH_TEMPLATE
        !
        SUBROUTINE MCGFFI_TEMPLATE(SUBSCH,K,KPN,M,N,H,NOM,NZON,XST,S,
     1             NREG,KEYFLX,KEYCUR,F,B,W,OMEGA2,IDIR,NSOUT,XSI)
          INTEGER K,KPN,M,N,NOM(N),NZON(K),NREG,KEYFLX(NREG,1),
     1     KEYCUR(K-NREG),IDIR,NSOUT
          REAL XST(0:M)
          DOUBLE PRECISION W,H(N),S(KPN),F(KPN),B(N),OMEGA2(3),
     1    XSI(NSOUT)
          EXTERNAL SUBSCH
        END SUBROUTINE MCGFFI_TEMPLATE
        !
        SUBROUTINE MCGFFA_TEMPLATE(SUBSCH,K,KPN,M,N,H,NOM,NZON,XST,SP,
     1             SM,NREG,NMU,NANI,NFUNL,NMOD,TRHAR,KEYFLX,KEYCUR,IMU,
     2             F,B,MODP,MODM)
          INTEGER K,KPN,M,N,NOM(N),NZON(K),NMU,NFUNL,NMOD,NREG,
     1     KEYFLX(NREG,NFUNL),KEYCUR(K-NREG),IMU,NANI,MODP,MODM
          REAL XST(0:M),TRHAR(NMU,NFUNL,NMOD)
          DOUBLE PRECISION H(N),SP(N),SM(N),F(KPN),B(2,N)
          EXTERNAL SUBSCH
        END SUBROUTINE MCGFFA_TEMPLATE
        !
        SUBROUTINE MCGSCH_TEMPLATE(N,K,M,NOM,NZON,H,XST,B)
          INTEGER N,K,M,NOM(N),NZON(K)
          REAL XST(0:M)
          DOUBLE PRECISION H(N),B(N)
        END SUBROUTINE MCGSCH_TEMPLATE
      END INTERFACE
      PROCEDURE(MOCFFI_TEMPLATE), POINTER :: MOCFFI
      PROCEDURE(MOCFFA_TEMPLATE), POINTER :: MOCFFA
      PROCEDURE(MOCSCH_TEMPLATE), POINTER :: MOCSCH
      PROCEDURE(MCGFFI_TEMPLATE), POINTER :: MCGFFI
      PROCEDURE(MCGFFA_TEMPLATE), POINTER :: MCGFFA
      PROCEDURE(MCGSCH_TEMPLATE), POINTER :: MCGSCH
      PROCEDURE(MOCFFI_TEMPLATE) :: MOCFFIS,MOCFFIR,MOCFFIT
      PROCEDURE(MOCFFA_TEMPLATE) :: MOCFFAS,MOCFFAR,MOCFFAT
      PROCEDURE(MOCSCH_TEMPLATE) :: MOCSCAS,MOCDDFS,MOCSCES,MOCSCA,
     1 MOCDDF,MOCSCE,MOCSCAT,MOCDDFT,MOCSCET,MOCSCEL,MOCDDFL,MOCSCAL
      PROCEDURE(MCGFFI_TEMPLATE) :: MCGFFIS,MCGFFIR,MCGFFIT
      PROCEDURE(MCGFFA_TEMPLATE) :: MCGFFAS,MCGFFAR,MCGFFAT
      PROCEDURE(MCGSCH_TEMPLATE) :: MCGSCAS,MCGDDFS,MCGSCES,MCGSCA,
     1 MCGDDF,MCGSCE,MCGSCAT,MCGDDFT,MCGSCET,MCGSCEL,MCGDDFL,MCGSCAL
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6,NSTATE=40,MXNMU=64)
      CHARACTER   TEXT4*4,TOPT*72
      INTEGER     JPAR(NSTATE),PACA,STIS,TRTY,IGB(8)
      REAL        ZREAL(4),HDD,DELU,FACSYM
      LOGICAL     CYCLIC,LVOID,LEXF,LFORW,LPRISM
      EXTERNAL    MOCFFAL,MCGFFAL
      INTEGER, TARGET, SAVE, DIMENSION(1) :: IDUMMY
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR) WZMU_PTR,ZMU_PTR,V_PTR,NZON_PTR,KEY_PTR,KEYCUR_PTR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITST,MATALB
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGAL,REPS,EPS,CPO
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: CAZ0,CAZ1,CAZ2
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: INCONV
      INTEGER, POINTER, DIMENSION(:) :: NZON,KEY,KEYCUR
      REAL, POINTER, DIMENSION(:) :: WZMU,ZMU,V
*
      IF(MAT(1).LT.0) CALL XABORT('MCCGF: EXPECTING MAT(1)>=0')
      IF(VOL(1).LT.0.0) CALL XABORT('MCCGF: EXPECTING VOL(1)>=0')
      IF(IMPX.GT.3) WRITE(IUNOUT,'(//8H MCCGF: ,A72/)') TITR
*----
*  RECOVER MCCG SPECIFIC PARAMETERS
*----
*     check for cross-sections in SYS object
      CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
      CALL LCMLEN(KPSYS(1),'DRAGON-TXSC',ILENG,ITYLCM)
      IF(ILENG.NE.NBMIX+1) CALL XABORT('MCCGF: INVALID VALUE OF NBMIX.')
      IF(JPAR(4).GT.NBMIX) CALL XABORT('MCCGF: MIXTURE OVERFLOW.')
*     check for a tracking binary file
      IF(IFTRAK.LE.0) CALL XABORT('MCCGF: INVALID TRACKING FILE.')
*     recover state-vector information
      IF(JPAR(40).EQ.1) THEN
         CALL LCMSIX(IPTRK,'BIHET',1)
         CALL LCMGET(IPTRK,'PARAM',IGB)
         NREG=IGB(3)
         CALL LCMSIX(IPTRK,' ',2)
      ELSE
         NREG=JPAR(1)
      ENDIF
      NSOU=JPAR(5)
      NFI=NREG+NSOU
      IF(JPAR(2).GT.NUNKNO) 
     1 CALL XABORT('MCCGF: UNKNOWN VECTOR OVERFLOW.')
      NANI=JPAR(6)
      TRTY=JPAR(9)
      IF(TRTY.EQ.1) THEN
         CYCLIC=.TRUE.
         NLONG=NREG
      ELSE
         CYCLIC=.FALSE.
         NLONG=NFI
      ENDIF
      NZP=JPAR(39)
      LPRISM=(NZP.NE.0)
      CALL LCMGET(IPTRK,'MCCG-STATE',JPAR)
      NMU=JPAR(2)
      IF(NMU.GT.MXNMU)
     1 CALL XABORT('MCCGF: POLAR ANGLE QUADRATURE OVERFLOW') 
      NMAX=JPAR(5)
      MAXI=JPAR(13)
      STIS=JPAR(15)
      LC=JPAR(6)
      IAAC=JPAR(7)
      KRYL=JPAR(3)
      IDIFC=JPAR(4)
      ISCR=JPAR(8)
      LPS=JPAR(9)
      PACA=JPAR(10)
      LEXF=(JPAR(12).EQ.1)
      LFORW=(JPAR(18).EQ.0)
      NFUNL=JPAR(19)
      NLIN=JPAR(20)
*     to be coherent with the exponential function used for the Pjj calculation
      IF((LEXAC).AND.(.NOT.LEXF).AND.(STIS.EQ.1)) STIS=0 
      NPJJM=JPAR(16)
*     recover real parameters
      CALL LCMGET(IPTRK,'REAL-PARAM',ZREAL)
      EPSI=ZREAL(1)
      DELU=ZREAL(3)
      FACSYM=ZREAL(4)
!!! temporary
      HDD=ZREAL(2)
      IF(HDD.GT.0.0) THEN
         ISCH=0
      ELSEIF(LEXF) THEN
         ISCH=-1
      ELSE
         ISCH=1
      ENDIF
*----
* RECOVER TRACKING FILE INFORMATION
*----
      REWIND IFTRAK
      READ(IFTRAK) TEXT4,NCOMNT,NBTR,IFMT
      DO ICOM=1,NCOMNT
         READ(IFTRAK)
      ENDDO
      READ(IFTRAK) NDIM,ISPEC,N2REG,N2SOU,NALBG,NCOR,NANGL,MXSUB,MXSEG
      IF(NCOR.NE.1) 
     1 CALL XABORT('MCCGF: INVALID TRACKING FILE: NCOR.NE.1')
      ALLOCATE(MATALB(N2REG+N2SOU+1))
      READ(IFTRAK)
      READ(IFTRAK) (MATALB(JJ),JJ=1,N2REG+N2SOU+1)
      READ(IFTRAK)
      READ(IFTRAK)
      ALLOCATE(CAZ0(NANGL),CAZ1(NANGL),CAZ2(NANGL),CPO(NMU))
      IF(NDIM.EQ.2) THEN
         CALL LCMGET(IPTRK,'XMU$MCCG',CPO)
         READ(IFTRAK) (CAZ1(JJ),CAZ2(JJ),JJ=1,NANGL)
      ELSE ! NDIM.EQ.3
**        correction Sylvie Musongela, december 2019
         READ(IFTRAK) (CAZ1(JJ),CAZ2(JJ),CAZ0(JJ),JJ=1,NANGL)
         DO JJ=1,NANGL
            CAZ1(JJ)=CAZ1(JJ)/SQRT(1.0D0-CAZ0(JJ)*CAZ0(JJ))
            CAZ2(JJ)=CAZ2(JJ)/SQRT(1.0D0-CAZ0(JJ)*CAZ0(JJ))
         ENDDO
      ENDIF
*----
* RECOVER TRACKING TABLE INFORMATION
*----
*     recover polar quadrature
      CALL LCMGPD(IPTRK,'WZMU$MCCG',WZMU_PTR)
      CALL LCMGPD(IPTRK,'ZMU$MCCG',ZMU_PTR)
*     recover modified MATALB, VOLSUR and KEYFLX
      CALL LCMGPD(IPTRK,'V$MCCG',V_PTR)
      CALL LCMGPD(IPTRK,'NZON$MCCG',NZON_PTR)
      CALL LCMGPD(IPTRK,'KEYFLX$ANIS',KEY_PTR)
*     recover index for the currents in FUNKNO (non-cyclic case)
      IF(.NOT.CYCLIC) CALL LCMGPD(IPTRK,'KEYCUR$MCCG',KEYCUR_PTR)
*
      CALL C_F_POINTER(WZMU_PTR,WZMU,(/ NMU /))
      CALL C_F_POINTER(ZMU_PTR,ZMU,(/ NMU /))
      CALL C_F_POINTER(V_PTR,V,(/ NLONG /))
      CALL C_F_POINTER(NZON_PTR,NZON,(/ NLONG /))
      CALL C_F_POINTER(KEY_PTR,KEY,(/ NREG*NLIN*NFUNL /))
      IF(.NOT.CYCLIC) THEN
         CALL C_F_POINTER(KEYCUR_PTR,KEYCUR,(/ NLONG-NBREG /))
      ELSE
         KEYCUR=>IDUMMY
      ENDIF
*----
*  CONSTRUCT TOTAL CROSS SECTIONS ARRAY AND CHECK FOR ZERO CROSS SECTION
*----
      CALL LCMLEN(KPSYS(1),'ALBEDO',NALBP,ITYLCM)
      ALLOCATE(SIGAL((NBMIX+7)*NGEFF))
      CALL MCGSIG(IPTRK,NBMIX,NGEFF,NALBP,KPSYS,SIGAL,LVOID)
      IF((LVOID).AND.(STIS.EQ.-1)) THEN
         IF(IMPX.GT.0) 
     1     WRITE(IUNOUT,*) 'VOID EXISTS -> STIS SET TO 1 INSTEAD OF -1'
         STIS=1
      ENDIF
      ISCH=ISCH+10*STIS+100*(NLIN-1)
*----
*  ASSIGN GENERIC INTERFACES
*----
      NULLIFY(MOCFFI)
      NULLIFY(MOCFFA)
      NULLIFY(MOCSCH)
      NULLIFY(MCGFFI)
      NULLIFY(MCGFFA)
      NULLIFY(MCGSCH)
      IF(CYCLIC) THEN
*     --------------------------------
*     Method of Cyclic Characteristics
*     --------------------------------
*********'Source Term Isolation' Strategy turned off 
         IF(ISCH.EQ.1) THEN
*          Step-Characteristics Scheme with Tabulated Exponentials
           TOPT='CYCLIC - STIS 0 - SC SCHEME - TABULATED EXP'
           MOCFFI => MOCFFIS
           MOCFFA => MOCFFAS
           MOCSCH => MOCSCAS
         ELSEIF(ISCH.EQ.0) THEN
*          Diamond-Differencing Scheme
           TOPT='CYCLIC - STIS 0 - DD0 SCHEME'
           MOCFFI => MOCFFIS
           MOCFFA => MOCFFAS
           MOCSCH => MOCDDFS
         ELSEIF(ISCH.EQ.-1) THEN
*          Step-Characteristics Scheme with Exact Exponentials
           TOPT='CYCLIC - STIS 0 - SC SCHEME - EXACT EXP'
           MOCFFI => MOCFFIS
           MOCFFA => MOCFFAS
           MOCSCH => MOCSCES
*********'Source Term Isolation' Strategy turned on
         ELSEIF(ISCH.EQ.11) THEN
*          Step-Characteristics Scheme with Tabulated Exponentials
           TOPT='CYCLIC - STIS 1 - SC SCHEME - TABULATED EXP'
           MOCFFI => MOCFFIR
           MOCFFA => MOCFFAR
           MOCSCH => MOCSCA
         ELSEIF(ISCH.EQ.10) THEN
*          Diamond-Differencing Scheme
           TOPT='CYCLIC - STIS 1 - DD0 SCHEME'
           MOCFFI => MOCFFIR
           MOCFFA => MOCFFAR
           MOCSCH => MOCDDF
         ELSEIF(ISCH.EQ.9) THEN
*          Step-Characteristics Scheme with Exact Exponentials
           TOPT='CYCLIC - STIS 1 - SC SCHEME - EXACT EXP'
           MOCFFI => MOCFFIR
           MOCFFA => MOCFFAR
           MOCSCH => MOCSCE
*********'MOCC/MCI' Iterative Strategy
         ELSEIF(ISCH.EQ.-9) THEN
*          Step-Characteristics Scheme with Tabulated Exponentials
           TOPT='CYCLIC - STIS -1 - SC SCHEME - TABULATED EXP'
           MOCFFI => MOCFFIT
           MOCFFA => MOCFFAT
           MOCSCH => MOCSCAT
         ELSEIF(ISCH.EQ.-10) THEN
*          Diamond-Differencing Scheme
           TOPT='CYCLIC - STIS -1 - DD0 SCHEME'
           MOCFFI => MOCFFIT
           MOCFFA => MOCFFAT
           MOCSCH => MOCDDFT
         ELSEIF(ISCH.EQ.-11) THEN
*          Step-Characteristics Scheme with Exact Exponentials
           TOPT='CYCLIC - STIS -1 - SC SCHEME - EXACT EXP'
           MOCFFI => MOCFFIT
           MOCFFA => MOCFFAT
           MOCSCH => MOCSCET
         ELSEIF(ISCH.EQ.199) THEN
*          Lin.-Disc.-Characteristics Scheme with Exact Exponentials
           TOPT='CYCLIC - STIS 0 - LDC SCHEME - EXACT EXP'
           MOCFFI => MOCFFIT
           MOCFFA => MOCFFAT
           MOCSCH => MOCSCEL
         ELSEIF(ISCH.EQ.200) THEN
*          Lin.-Disc.-Characteristics Scheme with Exact Exponentials
           TOPT='CYCLIC - STIS 0 - LDC SCHEME - DD1 SCHEME'
           MOCFFI => MOCFFIT
           MOCFFA => MOCFFAT
           MOCSCH => MOCDDFL
         ELSEIF(ISCH.EQ.201) THEN
*          Lin.-Disc.-Characteristics Scheme with Exact Exponentials
           TOPT='CYCLIC - STIS 0 - LDC SCHEME - TABULATED EXP'
           MOCFFI => MOCFFIT
           MOCFFA => MOCFFAT
           MOCSCH => MOCSCAL
         ELSE
           CALL XABORT('MCCGF: CYCLIC SCHEME NOT IMPLEMENTED')
         ENDIF   
      ELSE
*     ------------------------------------
*     Method of Non-Cyclic Characteristics
*     ------------------------------------
         IF(ISCH.EQ.1) THEN
*          Step-Characteristics Scheme with Tabulated Exponentials
           TOPT='NON CYCLIC - STIS 0 - SC SCHEME - TABULATED EXP'
           MCGFFI => MCGFFIS
           MCGFFA => MCGFFAS
           MCGSCH => MCGSCAS
         ELSEIF(ISCH.EQ.0) THEN
*          Diamond-Differencing Scheme
           TOPT='NON CYCLIC - STIS 0 - DD0 SCHEME'
           MCGFFI => MCGFFIS
           MCGFFA => MCGFFAS
           MCGSCH => MCGDDFS
         ELSEIF(ISCH.EQ.-1) THEN
*          Step-Characteristics Scheme with Exact Exponentials
           TOPT='NON CYCLIC - STIS 0 - SC SCHEME - EXACT EXP'
           MCGFFI => MCGFFIS
           MCGFFA => MCGFFAS
           MCGSCH => MCGSCES
         ELSEIF(ISCH.EQ.11) THEN
*          Step-Characteristics Scheme with Tabulated Exponentials
           TOPT='NON CYCLIC - STIS 1 - SC SCHEME - TABULATED EXP'
           MCGFFI => MCGFFIR
           MCGFFA => MCGFFAR
           MCGSCH => MCGSCA
         ELSEIF(ISCH.EQ.10) THEN
*          Diamond-Differencing Scheme
           TOPT='NON CYCLIC - STIS 1 - DD0 SCHEME'
           MCGFFI => MCGFFIR
           MCGFFA => MCGFFAR
           MCGSCH => MCGDDF
         ELSEIF(ISCH.EQ.9) THEN
*          Step-Characteristics Scheme with Exact Exponentials
           TOPT='NON CYCLIC - STIS 1 - SC SCHEME - EXACT EXP'
           MCGFFI => MCGFFIR
           MCGFFA => MCGFFAR
           MCGSCH => MCGSCE
         ELSEIF(ISCH.EQ.-9) THEN
*          Step-Characteristics Scheme with Tabulated Exponentials
           TOPT='NON CYCLIC - STIS -1 - SC SCHEME - TABULATED EXP'
           MCGFFI => MCGFFIT
           MCGFFA => MCGFFAT
           MCGSCH => MCGSCAT
         ELSEIF(ISCH.EQ.-10) THEN
*          Diamond-Differencing Scheme
           TOPT='NON CYCLIC - STIS -1 - DD0 SCHEME'
           MCGFFI => MCGFFIT
           MCGFFA => MCGFFAT
           MCGSCH => MCGDDFT
         ELSEIF(ISCH.EQ.-11) THEN
*          Step-Characteristics Scheme with Exact Exponentials
           TOPT='NON CYCLIC - STIS -1 - SC SCHEME - EXACT EXP'
           MCGFFI => MCGFFIT
           MCGFFA => MCGFFAT
           MCGSCH => MCGSCET
         ELSEIF(ISCH.EQ.199) THEN
*          Lin.-Disc.-Characteristics Scheme with Exact Exponentials
           TOPT='NON CYCLIC - STIS 0 - LDC SCHEME - EXACT EXP'
           MCGFFI => MCGFFIT
           MCGFFA => MCGFFAT
           MCGSCH => MCGSCEL
         ELSEIF(ISCH.EQ.200) THEN
*          Diamond-Differencing Scheme
           TOPT='NON CYCLIC - STIS 0 - DD1 SCHEME'
           MCGFFI => MCGFFIT
           MCGFFA => MCGFFAT
           MCGSCH => MCGDDFL
         ELSEIF(ISCH.EQ.201) THEN
*          Lin.-Disc.-Characteristics Scheme with Tabulated Exponentials
           TOPT='NON CYCLIC - STIS 0 - LDC SCHEME - TABULATED EXP'
           MCGFFI => MCGFFIT
           MCGFFA => MCGFFAT
           MCGSCH => MCGSCAL
         ELSE
           CALL XABORT('MCCGF: NON-CYCLIC SCHEME NOT IMPLEMENTED')
         ENDIF
      ENDIF
*----
*  PERFORM INNER ITERATIONS TO COMPUTE THE NEUTRON FLUX IN THE DIFFERENT
*  GROUPS
*----
      ALLOCATE(REPS(MAXI*NGEFF),EPS(NGEFF),ITST(NGEFF),INCONV(NGEFF))
      CALL XDLSET(INCONV,NGEFF,.TRUE.)
      LNCONV=NGEFF
*
      IF(IDIFC.EQ.1) THEN
*     ------------------------------------
*     ACA-Simplified Transport Calculation
*     ------------------------------------
         TOPT='ACA-SIMPLIFIED TRANSPORT OPERATOR'
         CALL MCGFLS(IMPX,IPTRK,IPMACR,NUNKNO,NFI,NBREG,NLONG,NBMIX,
     1        NGRP,NGEFF,LC,LFORW,PACA,NZON,KEY,KEYCUR,NGIND,KPSYS,
     2        INCONV,EPSI,MAXI,FUNKNO2,SUNKNO)
*     ------------------------------------
      ELSE
        IF(CYCLIC) THEN
*     --------------------------------
*     Method of Cyclic Characteristics
*     --------------------------------
          CALL MCGFLX(MOCFFI,MOCFFA,MOCSCH,MOCFFAL,CYCLIC,KPSYS,
     1       IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2       NSOU,NGRP,NGEFF,NGIND,NZON,MATALB,V,FUNKNO2,SUNKNO,
     3       NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4       CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5       KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6       STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
        ELSE
*     ------------------------------------
*     Method of Non-Cyclic Characteristics
*     ------------------------------------
          CALL MCGFLX(MCGFFI,MCGFFA,MCGSCH,MCGFFAL,CYCLIC,KPSYS,
     1        IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2        NSOU,NGRP,NGEFF,NGIND,NZON,MATALB,V,FUNKNO2,SUNKNO,
     3        NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4        CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5        KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6        STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
        ENDIF
      ENDIF
      DEALLOCATE(INCONV,SIGAL,CPO,CAZ2,CAZ1,CAZ0,MATALB)
*---
* PRINT RESULTS
*---
      IF(IMPX.GT.0) WRITE(IUNOUT,50) TOPT
      IF((IDIFC.EQ.0).AND.(MAXI.GT.1).AND.(IMPX.GT.1)) THEN     
      DO II=1,NGEFF
         IG=NGIND(II)
         WRITE(IUNOUT,100) IG
         IF(IMPX.GT.3) 
     1      WRITE(IUNOUT,200) (SUNKNO(KEYFLX(I),IG),I=1,NBREG)
         ITEMP=ITST(II)
         TEMP=EPS(II)
         IF((ITEMP.EQ.MAXI).AND.(TEMP.GT.EPSI)) WRITE(IUNOUT,60) 
         WRITE(IUNOUT,70) ITEMP,TEMP
         IF(IMPX.GT.2) THEN
           WRITE(IUNOUT,150) (REPS(MAXI*(II-1)+I),I=1,MIN(ITEMP,100))
         ENDIF
      ENDDO
      ENDIF
      DEALLOCATE(ITST,EPS,REPS)
      RETURN
*
 50   FORMAT(9X,18H M O C PARAMETERS:,2X,A72)
 60   FORMAT(49H *** WARNING *** MAXIMUM NUMBER OF MCCG ITERATION,
     1       10HS REACHED.)
 70   FORMAT(34H MCCGF: NUMBER OF MCCG ITERATIONS=,I4,
     1       11H  ACCURACY=,1P,E11.4,1H.)
 100  FORMAT(9X,8H GROUP (,I4,4H) : )
 150  FORMAT('----  EPS  ----'/(1P,6E16.6))
 200  FORMAT(/33H N E U T R O N    S O U R C E S :/(1P,6(5X,E15.7)))
      END
