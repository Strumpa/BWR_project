*DECK THMSAL
      SUBROUTINE THMSAL(IMPX,ITIME,I,J,K,K0,MFLOW,HMAVG,ENT,HD,SNAME,
     > SCOMP,IHCONV,KHCONV,ISUBM,RADCL,ZF,PHI,XFL,EPS,SLIP,DZ,TCALO,
     > RHO,RHOLAV,TSCLAD,KWA)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Adaptation of THMH2O.f for convection of Molten Salts using Gnielinski
* correlation.
*
*Copyright:
* Copyright (C) 2023 Ecole Polytechnique de Montreal.
*
*Author(s): 
* Cristian Garrido Tamm (cristian.garrido@idom.com)
*
*Parameters: input
* IMPX    printing index (=0 for no print).
* ITIME   type of calculation  (0=steady-state; 1=transient).
* I       position of channel alon X-axis
* J       position of channel alon Y-axis
* K       position of channel alon Z-axis
* K0      onser of nuclear boiling point
* MFLOW   massic coolant flow rate in Kg/m^2/s
* HMAVG   averaged enthalpy
* ENT     four values of enthalpy in J/Kg to be used in Gaussian
*         integration
* HD      hydraulic diameter in m
* SNAME   Name of the molten salt (e.g. "LiF-BeF2")
* SCOMP   Composition of the molten salt (e.g. "0.66-0.34")
* IHCONV  flag indicating HCONV chosen (0=default/1=user-provided).
* KHCONV  fixed user-provided HCONV value in W/m^2/K.
* ISUBM   subcooling model (0: one-phase; 1: Jens-Lottes model;
*         2: Saha-Zuber model).
* RADCL   outer clad radius in m
* ZF      parameters used to compute heat flux on clad surface in
*         transient cases.
* PHI     heat flow exchanged between clad and fluid in W/m^2.
*         Given in steady-state cases.
* XFL     input coolant flow quality
* EPS     input coolant void fraction
* SLIP    input slip ratio of vapor phase speed to liquid phase speed.
* DZ      axial mesh width in m.
*
*Parameters: output
* PHI     heat flow exchanged between clad and fluid in W/m^2.
*         Computed in transient cases.
* XFL     output coolant flow quality
* EPS     output coolant void fraction
* SLIP    output slip ratio of vapor phase speed to liquid phase speed.
* TCALO   coolant temperature in K
* RHO     coolant density in Kg/m^3
* RHOLAV  liquid density in kg/m^3
* TSCLAD  clad temperature in K
* KWA     flow regime (=0: single-phase; =1: subcooled; =2: nucleate
*         boiling; =3 superheated steam)
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER I,J,K,K0,IHCONV,ISUBM,KWA
      REAL MFLOW,HMAVG,ENT(4),HD,KHCONV,RADCL,ZF(2),PHI,TCALO,RHO,
     > RHOLAV,TSCLAD,XFL,EPS,SLIP,DZ
*----
*  LOCAL VARIABLES
*----
      REAL W(4),HL(4)
      CHARACTER HSMG*131
*----
*  SAVE VARIABLES
*----
      SAVE W
      DATA W /0.347855,0.652145,0.652145,0.347855/
*----
*  COMPUTE THE DENSITY AND TEMPERATURE OF THE LIQUID
*----
      IF(HMAVG.LT.0.0) CALL XABORT('THMSAL: NEGATIVE INPUT ENTHALPY.')
      IF(ISUBM.NE.0) CALL XABORT('THMSAL: NOT A ONE PHASE FLOW.')
      CALL THMSST(SNAME,SCOMP,TSAT,IMPX)
      HL(1)=ENT(1)
      HL(2)=ENT(2)
      HL(3)=ENT(3)
      HL(4)=ENT(4)
      CALL THMSH(SNAME,SCOMP,HL(1),R11,TL1,IMPX)
      CALL THMSH(SNAME,SCOMP,HL(2),R11,TL2,IMPX)
      CALL THMSH(SNAME,SCOMP,HL(3),R11,TL3,IMPX)
      CALL THMSH(SNAME,SCOMP,HL(4),R11,TL4,IMPX)
      CALL THMSPT(SNAME,SCOMP,TL1,RHO1,R2,R3,R4,CP1,IMPX)
      CALL THMSPT(SNAME,SCOMP,TL2,RHO2,R2,R3,R4,CP2,IMPX)
      CALL THMSPT(SNAME,SCOMP,TL3,RHO3,R2,R3,R4,CP3,IMPX)
      CALL THMSPT(SNAME,SCOMP,TL4,RHO4,R2,R3,R4,CP4,IMPX)
      TL=0.5*(W(1)*TL1+W(2)*TL2+W(3)*TL3+W(4)*TL4)
      RHOLAV=0.5*(W(1)*RHO1+W(2)*RHO2+W(3)*RHO3+W(4)*RHO4)
      CPLAV=0.5*(W(1)*CP1+W(2)*CP2+W(3)*CP3+W(4)*CP4)
*----
*  COMPUTE THE FLUID PROPERTIES
*  RHO: fluid density
*  REL: Reynolds number of liquid phase
*  PRL: Prandtl number of liquid phase
*----
      IF(XFL.NE.0.0) THEN
        CALL XABORT('THMSAL: INVALID VALUE OF FLOW QUALITY')
      ENDIF
*     One phase liquid
      TB=TSAT
      IF(TL.LT.TB) THEN
        TCALO=TL
      ELSE
        TCALO=TB
      ENDIF
      CALL THMSPT(SNAME,SCOMP,TCALO,R1,R2,ZKONE,ZMUONE,CPONE,IMPX)
      RHO=RHOLAV
      REL=MFLOW*HD/ZMUONE
      PRL=ZMUONE*CPONE/ZKONE
      ZKL=ZKONE
      XFL0=XFL
      EPS0=EPS
      SLIP0=SLIP
*----
*  THERMAL EXCHANGE BETWEEN CLAD AND FLUID USING THE DITTUS AND BOELTER
*  CORRELATION (SINGLE PHASE) OR CHEN CORRELATION (SATURATED BOILING)
*----
      IF(IHCONV.EQ.0) THEN
        ITER=0
        KWA=99
*CGT CHECK IF REYNOLDS AND PRANDTL ARE IN RANGE OF VALIDITY OF
*    GNIELINSKI CORRELATION
        TSCLAD=TCALO
        IF((REL.LT.2300).OR.(REL.GT.1E6)) THEN
          WRITE(6,*) "  THMSAL: ***WARNING*** REYNOLDS OUT RANGE."
        ENDIF
        IF((PRL.LT.0.6).OR.(PRL.GT.1E5)) THEN
          WRITE(6,*) "  THMSAL: ***WARNING*** PRANDTL OUT RANGE."
        ENDIF
        DO
          ITER=ITER+1
          IF(ITER.GT.50) THEN
            WRITE(HSMG,'(30HTHMSAL: HCONV FAILURE IN SLICE,I5,1H.)') K
            CALL XABORT(HSMG)
          ENDIF
*CGT Changed Dittus-Boelter by Gnielinski correlation
*CGT PRW: Prandtl number of liquid at wall temperature
          CALL THMSPT(SNAME,SCOMP,TSCLAD,R1,R2,ZKONE,ZMUONE,CPONE,IMPX)
          PRW=ZMUONE*CPONE/ZKONE
          HA=(ZKL/HD)*0.012*(REL**0.87-280)*PRL**0.8*(1+(HD/DZ)
     >    **(2.0/3.0))*(PRL/PRW)**0.11
          IF(IMPX.GT.4) THEN
            WRITE(6,*) 'THMSAL: REL,PRL,PRW,HA=',REL,PRL,PRW,HA
          ENDIF
          F=1.0
          S=1.0
          IF((XFL.EQ.XFL0).OR.(TSCLAD.LE.TSAT).OR.(KWA.EQ.0)) THEN
*           Single-phase convection. Use Gnielinski correlation
            KWA=0
            HB=0.0
            K0=0
            XFL=XFL0
            EPS=EPS0
            SLIP=SLIP0
          ELSE 
            CALL XABORT('THMSAL: INVALID HEAT TRANSFER REGIME')
          ENDIF
*         Chen correlation
          HCONV=F*HA+S*HB
          IF(HCONV.LE.0.0) THEN
            WRITE(HSMG,'(34HTHMSAL: DRY OUT REACHED IN CHANNEL,3I5)')
     >      I,J,K
            CALL XABORT(HSMG)
          ENDIF
          IF(ITIME.EQ.0) THEN
            TWAL=(PHI+S*HB*TSAT+F*HA*TCALO)/(S*HB+F*HA)
          ELSE
            ZNUM=ZF(1)+RADCL*S*HB*TSAT+RADCL*F*HA*TCALO
            ZDEN=ZF(2)+RADCL*S*HB+RADCL*F*HA
            TWAL=MAX(273.15,ZNUM/ZDEN)
            PHI=MAX(0.0,(ZF(1)-TWAL*ZF(2))/RADCL)
          ENDIF
          IF(ABS(TSCLAD-TWAL).GT.1.0E-5*TSCLAD) THEN
            TSCLAD=TWAL
          ELSE
            EXIT
          ENDIF
        ENDDO
      ELSE IF(IHCONV.EQ.1) THEN
        IF(ITIME.EQ.0) THEN
          TSCLAD=TCALO+PHI/KHCONV
        ELSE
          RCHC=RADCL*KHCONV
          TSCLAD=MAX(273.15,(ZF(1)+RCHC*TCALO)/(ZF(2)+RCHC))
          PHI=(ZF(1)-TSCLAD*ZF(2))/RADCL
        ENDIF
      ENDIF
      RETURN
      END
