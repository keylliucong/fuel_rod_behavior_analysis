C
C     ---------------------------------------------------------------
C
C     FUNCTION CPSF(PP)
C
C     FUNCTION CPSF IS FOR SPECIFIC HEAT OF SATURATED WATER
C
C     PP     : PRESSURE (MPA)                             (INPUT)
C     CPSF   : SPECIFIC HEAT OF SATURATED WATER(KJ/KG.K)  (OUTPUT)
C              LAST REVISION : APRIL 28 , 1992
C              MAXIMUN RELATIVE ERROR < 0.83%
C
C     ---------------------------------------------------------------
C
      FUNCTION CPSF(PP)
      IMPLICIT REAL(8) (A-H,O-Z)
      P=PP
      IF(P-2.3198D0) 10,20,30
  10  P=2.3198D0
      GOTO 20
  30  IF(P-18.675D0) 20,20,40
  40  P=18.675D0
  20  CP= .3538033720D+01
     &   +.8042237160D+00*P   -.2098582680D+00*P*P
     &   +.3018394110D-01*P**3-.1948064890D-02*P**4
     &   +.4982855530D-04*P**5-.8537371340D-07*P**6
      CPSF=CP
      RETURN
      END
C
C     ---------------------------------------------------------------
C
C     FUNCTION SHPT(PP,TT)
C
C     FUNCTION SHPT CALCULATS LIQUID SPECIFIC HEAT
C
C     PP     : PRESSURE (MPA)                             (INPUT)
C     TT     : TEMPRATURE (K)                             (INPUT)
C     SHPT   : SPECIFIC HEAT  (KJ/KG K)                   (OUTPUT)
C              THE MAXIMUM ERROR < 0.7%
C
C     ---------------------------------------------------------------
C
      FUNCTION SHPT(PP,TT)
      IMPLICIT REAL(8) (A-H,O-Z)
      TS=TSATP(PP)
      IF (TT.GE.TS-0.1D0) THEN
      SHPT=CPSF(PP)
      RETURN
      ELSE
      H1=HLPT(PP,TT)
      H2=HLPT(PP,TT-1.0D0)
      SHPT=H1-H2
      ENDIF
      RETURN
      END
C
C     -----------------------------------------------------------
C
C     FUNCTION  TSATP(P)
C
C     FUNCTION  TSATP CALCULATE STATURATION TEMPERATURE AT THE
C     PRESSURE
C
C     TSATP   : STATURATION TEMPERATURE(K)            (OUTPUT)
C     P       : PRESSURE(MPA)                          (INPUT)
C
C     REFERENCE : THE BWR PLANT ANALYZER , NUREG/CR - 3943
C     LAST REVISION : JULY 21 , 1991
C     REMARK  :   MAXIMUN RELATIVE ERROR LESS THAN 0.3%
C
C     -----------------------------------------------------------
C
      FUNCTION TSATP(P)
      IMPLICIT REAL(8) (A-H,O-Z)
      PP=P*1.D6
      FZ=0.0D0
      FM=0.0D0
      FZ= -3.84576467D+02       +0.57071646D+00*PP
     &    +1.04091792D-04*PP**2 +1.02949324D-09*PP**3
     &    +8.52096126D-16*PP**4 +5.56170847D-23*PP**5
      FM=  1.99967470D+01       +2.09487630D-02*PP
     &    +1.30393669D-06*PP**2 +6.75716812D-12*PP**3
     &    +3.17424682D-18*PP**4 +1.00000000D-25*PP**5
      TSAT=FZ/FM
      TSATP=TSAT+273.15D0
      RETURN
      END
C
C     ---------------------------------------------------
C
C     FUNCTION PSATT(TT)
C
C     FUNCTION PSATT CALCULATES THE SATURATED PRESSURE AT THE
C     TEMPERATURE
C
C     PSATT    : SATURATED PRESSURE(MPA)       (OUTPUT)
C     TT       : TEMPERATURE(K)                 (INPUT)
C
C     LAST REVISION : APRIL 28 ,1992
C     MAX. RELATIVE ERROR : 0.7%
C
C     ---------------------------------------------------
C
      FUNCTION PSATT(TT)
      IMPLICIT REAL(8) (A-H,O-Z)
      T=TT
      IF(T-453.15D0) 10,20,30
  10  T=453.15D0
      GOTO 20
  30  IF(T-543.15D0) 20,20,40
  20  PSATT =  .3285669320D+02
     &        -.1067699190D+00*T    -.3145188380D-04*T*T
     &        +.2387639020D-06*T**3 +.1894991280D-10*T**4
      GOTO 70
  40  IF(T-638.15D0) 50,50,60
  60  T=638.15D0
  50  PSATT =  .8335653680D+02
     &        -.2142781610D+00*T    -.1482528160D-03*T*T
     &        +.5137015930D-06*T**3
  70  CONTINUE
      RETURN
      END
C     ---------------------------------------------------
C
C     FUNCTION PSATT1(TT)
C
C     FUNCTION PSATT1 CALCULATES THE SATURATED PRESSURE AT THE
C     TEMPERATURE
C
C     PSATT1    : SATURATED PRESSURE(MPA)       (OUTPUT)
C     TT       : TEMPERATURE(K)                 (INPUT)
C
C     LAST REVISION : APRIL 28 ,1992
C     MAX. RELATIVE ERROR : 0.7%
C
C     ---------------------------------------------------
C
      FUNCTION PSATT1(TT)
      IMPLICIT REAL(8) (A-H,O-Z)
      T=TT
	IF(T.LE.273.15D0)T=273.15D0
	IF(T.GE.643.15D0)T=643.15D0
      IF(T.GT.543.15D0)THEN
        PSATT1 =  .8335653680D+02
     &        -.2142781610D+00*T    -.1482528160D-03*T*T
     &        +.5137015930D-06*T**3
	ELSE IF(T.GT.453.15D0)THEN
        PSATT1 =  .3285669320D+02
     &        -.1067699190D+00*T    -.3145188380D-04*T*T
     &        +.2387639020D-06*T**3 +.1894991280D-10*T**4
	ELSE
        PSATT1=2.728729496634253D-002+9.833874653893427D-003*T
	1        -1.466436282479518D-004*T**2+8.089164216185229D-007*T**3
     2        -1.983881976407186D-009*T**4+1.832465471847988D-012*T**5
	ENDIF
C
      RETURN
      END
C
C
C     --------------------------------------------------------------
C
C     FUNCTION TPH (PP,HH)
C
C     FUNCTION TPH IS FOR TEMPERATURE AS A FUNCTION OF PRESSURE AND
C     ENTHALPY
C
C     TPH    : TEMPERATURE AS A FUNCTION OF PRESSURE AND   (OUTPUT)
C              ENTHALPY(K)
C     PP   -- PRESSURE(MPA)                                  (INPUT)
C     HH   -- ENTHALPY(KJ/KG)                                (INPUT)
C
C     LAST REVISION   JULY 4 , 1991
C     REFERANCE  RETRAN--02--A PROGRAM FOR TRANSIENT THERMALHYDRAULIC
C     ANALYSIS  OF  COMPLEX  FLUID FLOW SYSTEMS .  VOLUME 1 :
C     EQUATIONS AND NUMERICS
C
C     --------------------------------------------------------------
C
      FUNCTION TPH(PP,HH)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION CT1(2,4),CT2(5,5),CT3(5,5),CT4(5,5)
      COMMON /CCC1/CT1,CT2,CT3,CT4
      PCRIT=3208.2D0
      HCRIT=906.0D0
      HFS=SFHP(PP)*0.429922613D0
      HGS=SGHP(PP)*0.429922613D0
      P=PP*145.05D0
      H=HH*0.429922613D0
      TPH=0.0D0
      IF(P.GT.PCRIT) GOTO 20
      IF(H.GT.HFS) GOTO 11
      DO 10 I=1,2
      DO 10 J=1,4
   10 TPH=TPH+CT1(I,J)*P**(I-1)*H**(J-1)
      GOTO 30
   11 IF(H.GT.HGS) GOTO 15
      DO 12 I1=1,2
      DO 12 J1=1,4
   12 TPH=TPH+CT1(I1,J1)*P**(I1-1)*HFS**(J1-1)
      GOTO 30
   15 DO 16 I2=1,5
      DO 16 J2=1,5
   16 TPH=TPH+CT3(I2,J2)*P**(I2-1)*H**(J2-1)
      GOTO 30
   20 HFCRIT=SFHP(PCRIT/145.05D0)*0.429922613D0
      HGCRIT=SGHP(PCRIT/145.05D0)*0.429922613D0
      IF(H.GT.HFCRIT) GOTO 25
      DO 22 I3=1,5
      DO 22 J3=1,5
   22 TPH=TPH+CT2(I3,J3)*P**(I3-1)*H**(J3-1)
      GOTO 30
  25  IF(H.LT.HGCRIT) GOTO 28
         DO 26 I4=1,5
      DO 26 J4=1,5
   26 TPH=TPH+CT4(I4,J4)*P**(I4-1)*H**(J4-1)
      GOTO 30
   28 WRITE(*,*) 'P (=',P,') > PCRIT (=3208.2PSIA) BUT '
      WRITE(*,*) 'HFCRIT (=',HFCRIT,') < H (=',H,') < HGCRIT (=',HGCRIT,
     &')'
      GOTO 31
   30 TPH=(TPH+459.67D0)/1.8D0
   31 RETURN
      END
C
C     -----------------------------------------------------------
C
C     FUNCTION HLPT(PP,TT)
C
C     FUNCTION HLPT CALCULATES LIQUID ENTHALPY
C
C     HLPT  : LIQUID ENTHALPY(KJ/KG)                (OUTPUT)
C     P     : PRESSURE(MPA)                          (INPUT)
C     T     : TEMPERATURE(K)                         (INPUT)
C
C     REMARK : RELATIVE ERROR < 1.D-2 IF P < 20 MPA
C     LAST REVISION : MAY 23 , 1992
C
C     -----------------------------------------------------------
C
      FUNCTION HLPT(PP,TT)
      IMPLICIT REAL(8) (A-H,O-Z)
      P=145.038D0*PP
      T=TT
      TS=TSATP(PP)
      IF(TT.GT.TS) T=TS
      T=1.8D0*T-459.67D0
      HLPT=0.001618D0*P-47.62D0
     &    +(0.965566D0-4.982D-6*P)*T
     &    +(8.2445D0*P-8843.4D0)/(802.0D0-T)
     &    +(6.6612D0*P-20392.0D0)/(T-750.0D0)
      HLPT=HLPT*2.32596004D0
      RETURN
      END
C
C     -------------------------------------------------------------
C
C     FUNCTION SFSV (PP)
C
C     FUNCTION SFSV CALCULATES SATURATED WATER SPECIFIC VOLUME WHEN
C     PRESSURE IS GIVEN
C
C     SFSV   : SATURATED WATER SPECIFIC VOLUME(M**3/KG)  (OUTPUT)
C     PP     : PRESSURE(MPA)                              (INPUT)
C
C     REMARK     MAXIMUN RELATIVE ERROR < 2.0%
C     REVISION   JUNE 11 , 1990
C
C     -------------------------------------------------------------
C
      FUNCTION SFSV(PP)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION CN1(3,5)
      COMMON /CCC2/CN1
      HFS=SFHP(PP)*0.429922613D0
      P=PP*145.05D0
      SFSV=0.0D0
      DO 5 I=1,3
      DO 5 J=1,5
    5 SFSV=SFSV+CN1(I,J)*P**(I-1)*HFS**(J-1)
      SFSV=DEXP(SFSV)
      SFSV=SFSV*0.0624280D0
      RETURN
      END
C
C     -------------------------------------------------------------
C
C     FUNCTION   SFSVDP (PP)
C
C     FUNCTION   SFSVDP MEANS D(SFSV)/D(PP) DERIVATIVE TO PRESSURE
C     FOR SATURATED WATER SPECIFIC VOLUME
C
C     SFSVDP  : DERIVATIVE TO PRESSURE FOR SATURATED    (OUTPUT)
C               WATER SPECIFIC VOLUME (M**3/KG.MPA)      (INPUT)
C     PP      : PRESSURE        (MPA)
C
C     LAST REVISION  :   JULY  5 , 1991
C
C     -------------------------------------------------------------
C
      FUNCTION SFSVDP(PP)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION CN1(3,5)
      COMMON /CCC2/CN1
      HFS=SFHP(PP)*0.429922613D0
      HFSDP=SFHPDP(PP)*0.429922613D0/145.05D0
      P=PP*145.05D0
      SFSV=0.0D0
      FVE=0.0D0
      DO 5 I=1,3
      DO 5 J=1,5
      SFSV=SFSV+CN1(I,J)*P**(I-1)*HFS**(J-1)
    5 FVE=FVE+CN1(I,J)*((I-1)*P**(I-2)*HFS**(J-1)
     &       +P**(I-1)*(J-1)*HFS**(J-2)*HFSDP)
      SFSVDP=DEXP(SFSV)*FVE
      SFSVDP=SFSVDP*0.0624280D0*145.0D0
      RETURN
      END
C
C     -------------------------------------------------------------
C
C     FUNCTION  SGSV (PP)
C
C     FUNCTION  SGSVCALCULATES THE  SATURATED STEAM SPECIFIC VOLUME
C     WHEN THE PRESSURE IS KNOWN
C
C     SGSV  : SATURATED STEAM SPECIFIC VOLUME(M**3/KG)  (OUTPUT)
C     PP    : PRESSURE(MPA)                              (INPUT)
C
C     REMARK     MAXIMUN RELATIVE ERROR < +/-2.5%
C     REVISION   JULY 4 , 1991
C
C     -------------------------------------------------------------
C
      FUNCTION SGSV(PP)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION CN2(4,3)
      COMMON /CCC3/CN2
      HGS=SGHP(PP)*0.429922613D0
      P=PP*145.05D0
      SGSV=0.0D0
      DO 5 I=1,4
      DO 5 J=1,3
    5 SGSV=SGSV+CN2(I,J)*P**(I-2)*HGS**(J-1)
      SGSV=SGSV*0.0624280D0
      RETURN
      END
C
C     -------------------------------------------------------------
C
C     FUNCTION   SGSVDP (PP)
C
C     FUNCTION SGSVDP CALCULATES  DERIVATIVE FOR SATURATED STEAM
C     SPECIFIC VOLUME
C
C     SGSVDP  : DERIVATIVE FOR SATURATED STEAM SPECIFIC  (OUTPUT)
C               VOLUME(M**3/KG.MPA)
C     PP      : PRESSURE(MPA)                             (INPUT)
C
C     LAST REVISION :  JULY 5 , 1991
C
C     ------------------------------------------------------------
C
      FUNCTION SGSVDP(PP)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION CN2(4,3)
      COMMON /CCC3/CN2
      HGS=SGHP(PP)*0.429922613D0
      HGSDP=SGHPDP(PP)*0.429922613D0/145.0D0
      P=PP*145.05D0
      SGSVDP=0.0D0
      DO 5 I=1,4
      DO 5 J=1,3
    5 SGSVDP=SGSVDP+CN2(I,J)*((I-2)*P**(I-3)*HGS**(J-1)
     &                        +P**(I-2)*(J-1)*HGS**(J-2)*HGSDP)
      SGSVDP=SGSVDP*0.0624280D0*145.0D0
      RETURN
      END
C
C     -------------------------------------------------------------
C
C     FUNCTION  SVPH (PP,HH)
C
C     FUNCTION SVPH IS FOR SPECIFIC VOLUME AS A FUNCTION OF PRESSURE
C     AND ENTHALPY
C
C     SVPH   : SPECIFIC VOLUME AS A FUNCTION OF         (OUTPUT)
C              PRESSURE AND ENTHALPY(M**3/KG)
C     PP     : PRESSURE(MPA)                             (INPUT)
C     HH     : ENTHALPY(KJ/KG)                           (INPUT)
C
C     LAST REVISION   JULY 4 , 1991
C     REFERANCE  RETRAN--02--A PROGRAM FOR TRANSIENT THERMAL-
C     HYDRAULIC  ANALYSIS  OF  COMPLEX  FLUID FLOW SYSTEMS .
C     VOLUME 1 : EQUATIONS AND NUMERICS
C
C     -------------------------------------------------------------
C
      FUNCTION SVPH(PP,HH)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION CG1(3,5),CG2(4,3)
      COMMON /CCC4/CG1,CG2
      PCRIT=3208.2D0
      P=PP*145.05D0
      HFS=SFHP(PP)*0.429922613D0
      HGS=SGHP(PP)*0.429922613D0
      H=HH*0.429922613D0
  334 CONTINUE
      SVPH=0.0D0
      IF(P.GE.PCRIT) GOTO 20
      IF(H.GT.HFS) GOTO 10
    2 DO 5 I=1,3
      DO 5 J=1,5
    5 SVPH=SVPH+CG1(I,J)*P**(I-1)*H**(J-1)
      SVPH=DEXP(SVPH)
      GOTO 30
   10 IF(H.LT.HGS) GOTO 25
   12 DO 15 I1=1,4
      DO 15 J1=1,3
   15 SVPH=SVPH+CG2(I1,J1)*P**(I1-2)*H**(J1-1)
      GOTO 30
   20 IF(H.GT.1000.0D0) GOTO 12
      GOTO 2
   25 XPH=(H-HFS)/(HGS-HFS)
      SFSV=0.0D0
      DO 26 I2=1,3
      DO 26 J2=1,5
   26 SFSV=SFSV+CG1(I2,J2)*P**(I2-1)*HFS**(J2-1)
      SFSV=DEXP(SFSV)
      SGSV=0.0D0
      DO 27 I3=1,4
      DO 27 J3=1,3
   27 SGSV=SGSV+CG2(I3,J3)*P**(I3-2)*HGS**(J3-1)
      SVPH=SFSV+XPH*(SGSV-SFSV)
   30 SVPH=SVPH*0.0624280D0
      RETURN
      END
C
C     -------------------------------------------------------------
C
C     FUNCTION  SVPHDP (PP,HH)
C
C     FUNCTION SVPHDP MEANS D(SVPH)/D(PP) DERIVATIVE TO PRESSURE
C     FOR SPECIFIC VOLUME AS A FUNCTION OF PRESSURE AND ENTHALPY
C
C     SVPHDP : D(SVPH)/D(PP) DERIVATIVE TO PRESSURE   (OUTPUT)
C              FOR SPECIFIC VOLUME AS A FUNCTION OF
C              PRESSURE AND ENTHALPY(M**3/KG.MPA)
C     PP     : PRESSURE(MPA)                           (INPUT)
C     HH     : ENTHALPY(KJ/KG)                         (INPUT)
C
C     REVISION        JULY 5 , 1991
C
C     -------------------------------------------------------------
C
      FUNCTION SVPHDP(PP,HH)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION CG1(3,5),CG2(4,3)
      COMMON /CCC4/CG1,CG2
      PCRIT=3208.2D0
      P=PP*145.05D0
      H=HH*0.429922613D0
      HFS=SFHP(PP)*0.429922613D0
      HGS=SGHP(PP)*0.429922613D0
      HHSF=SFHP(PP)*0.429922613D0
      HHSG=SGHP(PP)*0.429922613D0
      H=HH*0.429922613D0
      IF(H.LT.0.0D0) H=HHSF*0.75D0
      IF(H.GT.1.5D0*HHSG) H=1.5D0*HHSG
      SVPH=0.0D0
      SVPHDP=0.0D0
      IF(P.GE.PCRIT) GOTO 20
      IF(H.GT.HFS) GOTO 10
    2 DO 5 I=1,3
      DO 5 J=1,5
      SVPH=SVPH+CG1(I,J)*P**(I-1)*H**(J-1)
    5 SVPHDP=SVPHDP+(I-1)*CG1(I,J)*P**(I-2)*H**(J-1)
      SVPHDP=DEXP(SVPH)*SVPHDP
      GOTO 30
   10 IF(H.LT.HGS) GOTO 25
   12 DO 15 I1=1,4
      DO 15 J1=1,3
   15 SVPHDP=SVPHDP+(I1-2)*CG2(I1,J1)*P**(I1-3)*H**(J1-1)
      GOTO 30
   20 IF(H.GT.1000.0D0) GOTO 12
      GOTO 2
   25 CONTINUE
C      HFSDP=SFHPDP(PP)*2.9639615
C      HGSDP=SGHPDP(PP)*2.9639615
      HFSDP=SFHPDP(PP)*0.429922613D0/145.05D0
      HGSDP=SGHPDP(PP)*0.429922613D0/145.05D0
      XPH=(H-HFS)/(HGS-HFS)
      XPHDP=(-HFSDP*(HGS-HFS)-(H-HFS)*(HGSDP-HFSDP))/((HGS-HFS)**2)
      SFSV1=0.0D0
      SFSVDP=0.0D0
      DO 26 I2=1,3
      DO 26 J2=1,5
      SFSV1=SFSV1+CG1(I2,J2)*P**(I2-1)*HFS**(J2-1)
   26 SFSVDP=SFSVDP+CG1(I2,J2)*((I2-1)*P**(I2-2)*HFS**(J2-1)
     &                         +(J2-1)*P**(I2-1)*HFS**(J2-2)*HFSDP)
      SFSV=DEXP(SFSV1)
      SFSVDP=SFSV*SFSVDP
      SGSV=0.0D0
      SGSVDP=0.0D0
      DO 27 I3=1,4
      DO 27 J3=1,3
      SGSV=SGSV+CG2(I3,J3)*P**(I3-2)*HGS**(J3-1)
   27 SGSVDP=SGSVDP+CG2(I3,J3)*((I3-2)*P**(I3-3)*HGS**(J3-1)
     &                         +(J3-1)*P**(I3-2)*HGS**(J3-2)*HGSDP)
	SVPHDP=SFSVDP+XPHDP*(SGSV-SFSV)+XPH*(SGSVDP-SFSVDP)
   30 SVPHDP=SVPHDP*0.0624280D0*145.05D0
      SVPHDP=-DABS(SVPHDP)
      RETURN
      END
C
C     -------------------------------------------------------------
C
      FUNCTION SVPHDPT(PP,HH)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION CG1(3,5),CG2(4,3)
      COMMON /CCC4/CG1,CG2
      HFS=SFHP(PP)
      HGS=SGHP(PP)
      HHSF=SFHP(PP)
      HHSG=SGHP(PP)
	SVF=SFSV(PP)
	SVG=SGSV(PP)
	SVFDP=SFSVDP(PP)
	SVGDP=SGSVDP(PP)
      HFSDP=SFHPDP(PP)
      HGSDP=SGHPDP(PP)
	X=(HH-HHSF)/(HHSG-HHSF)
	DXDP=-1.D0/(HHSG-HHSF)*((1.D0-X)*HFSDP+X*HGSDP)
	SVPHDPT=DXDP*(SVG-SVF)+(1.D0-X)*SVFDP+X*SVGDP
      RETURN
      END
C
C     -------------------------------------------------------------
C
C     FUNCTION  SVPHDH (PP,HH)
C
C     FUNCTION SVPHDH MEANS D(SVPH)/D(HH) DERIVATIVE TO ENTHALPY
C     FOR SPECIFIC VOLUME AS A FUNCTION OF PRESSURE AND ENTHALPY
C
C     SVPHDH  : D(SVPH)/D(HH) DERIVATIVE TO ENTHALPY     (OUTPUT)
C                FOR SPECIFIC VOLUME AS A FUNCTION OF
C                PRESSURE AND ENTHALPY  (M**3/KJ)
C     PP      : PRESSURE(MPA)                             (INPUT)
C     HH      : ENTHALPY(KJ/KG)                           (INPUT)
C
C     REVISION        JULY 5 , 1991
C
C     -------------------------------------------------------------
C
      FUNCTION SVPHDH(PP,HH)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION CG1(3,5),CG2(4,3)
      COMMON /CCC4/CG1,CG2
      PCRIT=3208.2D0
      P=PP*145.05D0
      HFS=SFHP(PP)*0.429922613D0
      HGS=SGHP(PP)*0.429922613D0
      H=HH*0.429922613D0
      IF(H.LT.0.0D0) H=HFS*0.75D0
      IF(H.GT.1.5D0*HGS) H=1.5D0*HGS
      SVPH=0.0D0
      SVPHDH=0.0D0
      IF(P.GE.PCRIT) GOTO 20
      IF(H.GT.HFS) GOTO 10
    2 DO 5 I=1,3
      DO 5 J=1,5
      SVPH=SVPH+CG1(I,J)*P**(I-1)*H**(J-1)
    5 SVPHDH=SVPHDH+(J-1)*CG1(I,J)*P**(I-1)*H**(J-2)
      SVPHDH=DEXP(SVPH)*SVPHDH
      GOTO 30
   10 IF(H.LT.HGS) GOTO 25
   12 DO 15 I1=1,4
      DO 15 J1=1,3
   15 SVPHDH=SVPHDH+(J1-1)*CG2(I1,J1)*P**(I1-2)*H**(J1-2)
      GOTO 30
   20 IF(H.GT.1000.0D0) GOTO 12
      GOTO 2
   25 XPHDH=1.D0/(HGS-HFS)
      SFSV=0.0D0
      DO 26 I2=1,3
      DO 26 J2=1,5
   26 SFSV=SFSV+CG1(I2,J2)*P**(I2-1)*HFS**(J2-1)
      SFSV=DEXP(SFSV)
      SGSV=0.0D0
      DO 27 I3=1,4
      DO 27 J3=1,3
   27 SGSV=SGSV+CG2(I3,J3)*P**(I3-2)*HGS**(J3-1)
      SVPHDH=XPHDH*(SGSV-SFSV)
   30 SVPHDH=SVPHDH*0.0624280D0*0.429922613D0
      SVPHDH=DABS(SVPHDH)
      RETURN
      END
C
C     -------------------------------------------------------------
C
      FUNCTION SVPHDHT(PP)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION CG1(3,5),CG2(4,3)
      COMMON /CCC4/CG1,CG2
      HFS=SFHP(PP)
      HGS=SGHP(PP)
	SVF=SFSV(PP)
	SVG=SGSV(PP)
      SVPHDHT=(SVG-SVF)/(HGS-HFS)
      RETURN
      END
C
C     -------------------------------------------------------------
C
C     FUNCTION SFHP (PP)
C
C     FUNCTION SFHP CALCULATES SATURATED LIQUID ENTHALPY
C
C     SFHP  : SATURATED LIQUID ENTHALPY(KJ/KG)     (OUTPUT)
C     PP    : PRESSURE(MPA)                         (INPUT)
C
C     LAST REVISION   JULY 4 , 1991
C     REMARK     MAXIMUN RELATIVE ERROR LESS THAN 2.6%
C
C     -------------------------------------------------------------
C
      FUNCTION SFHP(PP)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION CSF1(9),CSF2(9),CSF3(9)
      COMMON /CCC5/CSF1,CSF2,CSF3
      IF(PP.GE.10.0D0.AND.PP.LE.20.0D0) THEN
      SFHP=SFHP1(PP)
      ELSE
      P=PP*145.05D0
      PCRIT=3208.506D0
      IF(P.LT..1D0) GOTO 6
      IF(P.GT.3208.5D0) GOTO 7
      SFHP=0.0D0
      IF(P.LT.950.0D0) GOTO 2
      IF(P.LT.2350.0D0) GOTO 4
      DO 1 I1=1,9
    1 SFHP=SFHP+CSF3(I1)*((PCRIT-P)**0.41D0)**(I1-1)
      GOTO 8
    2 DO 3 I2=1,9
    3 SFHP=SFHP+CSF1(I2)*(DLOG(P))**(I2-1)
      GOTO 8
    4 DO 5 I3=1,9
    5 SFHP=SFHP+CSF2(I3)*(DLOG(P))**(I3-1)
      GOTO 8
    6 WRITE(*,*) 'P < .1D0 PSIA OR 0.000689 MPA'
      GOTO 9
    7 WRITE(*,*) 'P > 3208.5 PSIA OR 22.11789 MPA'
      GOTO 9
    8 SFHP=SFHP*2.326D0
    9 RETURN
      ENDIF
      RETURN
      END

      FUNCTION SFHP1(PP)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION P(51),H(51)
      DATA P/10.D0,10.2D0,10.4D0,10.6D0,10.8D0,
     &       11.D0,11.2D0,11.4D0,11.6D0,11.8D0,
     &       12.D0,12.2D0,12.4D0,12.6D0,12.8D0,
     &       13.D0,13.2D0,13.4D0,13.6D0,13.8D0,
     &       14.D0,14.2D0,14.4D0,14.6D0,14.8D0,
     &       15.D0,15.2D0,15.4D0,15.6D0,15.8D0,
     &       16.D0,16.2D0,16.4D0,16.6D0,16.8D0,
     &       17.D0,17.2D0,17.4D0,17.6D0,17.8D0,
     &       18.D0,18.2D0,18.4D0,18.6D0,18.8D0,
     &       19.D0,19.2D0,19.4D0,19.6D0,19.8D0,
     &       20.D0/
      DATA H/1408.5D0,1416.8D0,1425.3D0,1433.7D0,1442.2D0,
     &       1450.8D0,1458.9D0,1467.2D0,1475.5D0,1483.9D0,
     &       1492.2D0,1500.2D0,1508.2D0,1516.1D0,1524.3D0,
     &       1532.5D0,1540.6D0,1548.4D0,1556.0D0,1564.2D0,
     &       1571.7D0,1580.0D0,1587.6D0,1595.7D0,1603.1D0,
     &       1611.0D0,1618.7D0,1627.2D0,1634.8D0,1642.7D0,
     &       1650.9D0,1658.7D0,1666.7D0,1674.5D0,1683.0D0,
     &       1691.7D0,1700.4D0,1709.0D0,1717.6D0,1726.2D0,
     &       1734.8D0,1743.4D0,1752.1D0,1760.9D0,1769.7D0,
     &       1778.7D0,1787.8D0,1797.0D0,1806.5D0,1816.3D0,
     &       1826.5D0/
      SFHP1=1826.5D0
      DO 5 I=1,50
      IF(PP.GE.P(I).AND.PP.LT.P(I+1)) THEN
      SFHP1=H(I)+(H(I+1)-H(I))/(P(I+1)-P(I))*(PP-P(I))
      ENDIF
  5   CONTINUE
      RETURN
      END

      FUNCTION SGHP1(PP)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION P(51),H(51)
      DATA P/10.D0,10.2D0,10.4D0,10.6D0,10.8D0,
     &       11.D0,11.2D0,11.4D0,11.6D0,11.8D0,
     &       12.D0,12.2D0,12.4D0,12.6D0,12.8D0,
     &       13.D0,13.2D0,13.4D0,13.6D0,13.8D0,
     &       14.D0,14.2D0,14.4D0,14.6D0,14.8D0,
     &       15.D0,15.2D0,15.4D0,15.6D0,15.8D0,
     &       16.D0,16.2D0,16.4D0,16.6D0,16.8D0,
     &       17.D0,17.2D0,17.4D0,17.6D0,17.8D0,
     &       18.D0,18.2D0,18.4D0,18.6D0,18.8D0,
     &       19.D0,19.2D0,19.4D0,19.6D0,19.8D0,
     &       20.D0/
      DATA H/2727.7D0,2724.2D0,2720.6D0,2716.9D0,2713.1D0,
     &       2709.3D0,2705.4D0,2701.5D0,2697.4D0,2693.3D0,
     &       2689.1D0,2684.9D0,2680.6D0,2676.1D0,2671.6D0,
     &       2667.0D0,2662.3D0,2657.5D0,2652.5D0,2647.6D0,
     &       2642.4D0,2637.2D0,2631.7D0,2626.3D0,2620.6D0,
     &       2614.9D0,2609.3D0,2603.2D0,2597.3D0,2591.2D0,
     &       2584.8D0,2578.5D0,2571.9D0,2565.0D0,2558.2D0,
     &       2551.3D0,2544.2D0,2536.9D0,2529.4D0,2521.7D0,
     &       2513.8D0,2505.7D0,2497.3D0,2488.7D0,2479.8D0,
     &       2470.6D0,2461.1D0,2451.1D0,2440.7D0,2429.8D0,
     &       2418.4D0/
      SGHP1=2418.4D0
      DO 5 I=1,50
      IF(PP.GE.P(I).AND.PP.LT.P(I+1)) THEN
      SGHP1=H(I)+(H(I+1)-H(I))/(P(I+1)-P(I))*(PP-P(I))
      ENDIF
  5   CONTINUE
      RETURN
      END
C
C     -------------------------------------------------------------
C
C              FUNCTION  SFHPDP (PP)
C
C     FUNCTION SFHPDP MEANS DSFHP/DP, DERIVATIVE FOR SATURATED
C     LIQUID ENTHALPY
C
C     SFHPDP  : DSFHP/DP DERIVATIVE FOR SATURATED   (OUTPUT)
C               LIQUID ENTHALPY  (KJ/KG.MPA)
C     PP      : PRESSURE(MPA)                        (INPUT)
C
C     LAST REVISION  :  JULY 5 , 1991
C
C     -------------------------------------------------------------
C
      FUNCTION SFHPDP(PP)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION CSF1(9),CSF2(9),CSF3(9)
      COMMON /CCC5/CSF1,CSF2,CSF3
      P=PP*145.05D0
      PCRIT=3208.506D0
      IF(P.LT..1D0) GOTO 6
      IF(P.GT.3208.5D0) GOTO 7
      SFHPDP=0.0D0
      IF(P.LT.950.0D0) GOTO 2
      IF(P.LT.2350.0D0) GOTO 4
      DO 1 I1=1,9
    1 SFHPDP=SFHPDP+(I1-1)*CSF3(I1)*((PCRIT-P)**0.41D0)**(I1-2)
     &     *0.41D0*((PCRIT-P)**(-0.59D0))*(-1.D0)
      GOTO 8
    2 DO 3 I2=1,9
    3 SFHPDP=SFHPDP+(I2-1)*CSF1(I2)*(DLOG(P))**(I2-2)/P
      GOTO 8
    4 DO 5 I3=1,9
    5 SFHPDP=SFHPDP+(I3-1)*CSF2(I3)*(DLOG(P))**(I3-2)/P
      GOTO 8
    6 WRITE(*,*) 'P < .1D0 PSIA OR 0.000689 MPA'
      GOTO 9
    7 WRITE(*,*) 'P > 3208.5 PSIA OR 22.11789 MPA'
      GOTO 9
    8 SFHPDP=SFHPDP*2.326D0*145.0D0
    9 RETURN
      END
C
C     -------------------------------------------------------------
C
C     FUNCTION  SGHP (PP)
C
C     FUNCTION SGHP CALCULATES SATURATED GAS ENTHALPY
C
C     SGHP   : SATURATED GAS ENTHALPY  (KJ/KG)      (OUTPUT)
C     PP     : PRESSURE(MPA)                         (INPUT)
C
C     LAST REVISION   JULY 4 , 1991
C     REMARK     MAXIMUN RELATIVE ERROR LESS THAN 1.0%
C
C     -------------------------------------------------------------
C
      FUNCTION SGHP(PP)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION CSG1(12),CSG2(9),CSG3(7)
      COMMON /CCC6/CSG1,CSG2,CSG3
      IF(PP.GE.10.0D0.AND.PP.LE.20.0D0) THEN
      SGHP=SGHP1(PP)
      ELSE
      P=PP*145.05D0
      PCRIT=3208.506D0
      IF(P.LT..1D0) GOTO 6
      IF(P.GT.3208.5D0) GOTO 7
      SGHP=0.0D0
      IF(P.LT.950.0D0) GOTO 2
      IF(P.LT.2650.0D0) GOTO 4
      DO 1 I1=1,7
    1 SGHP=SGHP+CSG3(I1)*((PCRIT-P)**0.41D0)**(I1-1)
      GOTO 8
    2 DO 3 I2=1,12
    3 SGHP=SGHP+CSG1(I2)*(DLOG(P))**(I2-1)
      GOTO 8
    4 DO 5 I3=1,9
    5 SGHP=SGHP+CSG2(I3)*(DLOG(P))**(I3-1)
      GOTO 8
    6 WRITE(*,*) 'P < .1D0 PSIA OR 0.000689 MPA'
      GOTO 9
    7 WRITE(*,*) 'P > 3208.5 PSIA OR 22.11789 MPA'
      GOTO 9
    8 SGHP=SGHP*2.326D0
    9 RETURN
      ENDIF
      RETURN
      END
C
C     -------------------------------------------------------------
C
C     FUNCTION SGHPDP (PP)
C
C     FUNCTION SGHP CALCULATES DERRIVATIVE FOR SATURATED STEAM
C     ENTHALPY
C
C     SGHP  : DERRIVATIVE FOR SATURATED STEAM      (OUTPUT)
C             ENTHALPY  (KJ/KG.MPA)
C     PP    : PRESSURE(MPA)                         (INPUT)
C
C     LAST REVISION  :  JULY 5 , 1991
C
C     -------------------------------------------------------------
C
      FUNCTION SGHPDP(PP)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION CSG1(12),CSG2(9),CSG3(7)
      COMMON /CCC6/CSG1,CSG2,CSG3
      P=PP*145.05D0
      PCRIT=3208.506D0
      IF(P.LT..1D0) GOTO 6
      IF(P.GT.3208.5D0) GOTO 7
      SGHPDP=0.0D0
      IF(P.LT.950.0D0) GOTO 2
      IF(P.LT.2650.0D0) GOTO 4
      DO 1 I1=1,7
    1 SGHPDP=SGHPDP+(I1-1)*CSG3(I1)*((PCRIT-P)**0.41D0)**(I1-2)
     &              *0.41D0*((PCRIT-P)**(-0.59D0))*(-1.D0)
      GOTO 8
    2 DO 3 I2=1,12
    3 SGHPDP=SGHPDP+(I2-1)*CSG1(I2)*(DLOG(P))**(I2-2)/P
      GOTO 8
    4 DO 5 I3=1,9
    5 SGHPDP=SGHPDP+(I3-1)*CSG2(I3)*(DLOG(P))**(I3-2)/P
      GOTO 8
    6 WRITE(*,*) 'P < .1D0 PSIA OR 0.000689 MPA'
      GOTO 9
    7 WRITE(*,*) 'P > 3208.5 PSIA OR 22.11789 MPA'
      GOTO 9
    8 SGHPDP=SGHPDP*2.326D0*145.0D0
    9 RETURN
      END
C
C     -------------------------------------------------------------
C
C     FUNCTION   THCOND(T,RHO)
C
C     FUNCTION THCOND CALCULATES THERMAL CONDUCTIVITY COEFFICIENT
C
C     THCOND  : THERMAL CONDUCTIVITY COEFFICIENT(W/K.M)  (OUTPUT)
C     T       : TEMPERATURE(K) (<1773.15D0 K)               (INPUT)
C     DENS    : DENSITY(KG/M**3)                          (INPUT)
C
C     REFRANCE   "PROPORTIES OF WATER AND STEAM IN SI-UNITS"
C                SECOND , REVISED AND , UPDATED . <<CHINESE
C                TRANSLATION>> , 1983 . PAGE 189 .
C     LAST REVISION :  JUNE 26 , 1992
C
C     -------------------------------------------------------------
C
      FUNCTION THCOND(T,RHO)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION A(4),B(3),BB(2),C(6),D(4)
      DATA FACT/647.3D0/,FACE/317.7D0/
      DATA A/1.02811D-2,2.99621D-2,1.56146D-2,-4.22464D-3/
      DATA B/-3.9707D-1,4.00302D-1,1.06D0/
      DATA BB/-1.71587D-1,2.39219D0/
      DATA C/6.42857D-1,-4.11717D0,-6.17937D0,3.08976D-3,
     &    8.22994D-2,1.00932D1/
      DATA D/7.01309D-2,1.1852D-2,1.69937D-3,-1.02D0/
      TH=T/FACT
      RH=RHO/FACE
      DE=0.
      DO 10 K=1,4
10    DE=DE+A(K)*TH**(K-1)
      VLO=TH**0.5D0*DE
      VLA=B(1)
      VLA=VLA+B(2)*RH
      X=(RH+BB(2))**2
      VLA=VLA+B(3)*DEXP(BB(1)*X)
      TE=DABS(TH-1.D0)+C(4)
      Q=2.+C(5)*TE**(-0.6D0)
      R=Q+1.D0
      IF(TH-1.D0) 1,2,2
1     S=C(6)*TE**(-0.6D0)
      GOTO 3
2     S=TE**(-1.D0)
3     FE=(D(1)*(1.D0/TH)**10+D(2))
     1   *RH**1.8D0*DEXP(C(1)*(1.D0-RH**2.8D0))
      FG=D(3)*S*RH**Q*DEXP(Q*(1.D0-RH**R)/R)
      Y=C(2)*TH**1.5D0+C(3)*(1.D0/RH)**5
      FH=D(4)*DEXP(Y)
      VLB=FE+FG+FH
      THCOND=VLO+VLA+VLB
      RETURN
      END
C
C     ---------------------------------------------------
C
C     FUNCTION PRF(PP)
C
C     FUNCTION PRF CALCULATES THE PRANDEL NUMBER OF
C     SATURATED WATER
C
C     PRF  : PRANDEL NUMBER OF SATURATED WATER    (OUTPUT)
C     PP   : PRESSURE(MPA)                (INPUT)
C
C     REMARK : MAX. RELATIVE ERROR < 6.0%
C     LAST REVISOIN : JUNE 16 , 1992
C
C     ---------------------------------------------------
C
      FUNCTION PRF(PP)
      IMPLICIT REAL(8) (A-H,O-Z)
      P=PP
      IF(P-1.555D0) 10,20,30
  10  P=1.555D0
      GOTO 20
  30  IF(P-18.67D0) 20,20,40
  40  P=18.67D0
  20  PRF= .8337492340D+00
     &   +.4784865670D-01*P-.1042182550D-01*P*P
     &   +.6044807370D-03*P**3
      RETURN
      END
C
C     ----------------------------------------------------------
C
C     FUNCTION   THCONL(PT,N)
C
C     FUNCTION THCONL IS FOR THERMAL CONDUCTIVITY OF STATURATION
C     WATER AS A FUNCTION OF PRESSURE OR STATURATION TEMPERATURE
C
C     THCONL  : THERMAL CONDUCTIVITY OF STATURATION WATER
C               AS A FUNCTION OF PRESSURE OR STATURATION
C               TEMPERATURE(W/M.K)
C     PT      : PRESSURE (MPA) , WHEN N=1 ;
C               STATURATION TEMPERATURE (K) , WHEN N=2 .
C
C     REFERENCE : THE BWR PLANT ANALYZER . NUREG/CR-3943
C     LAST REVISION : 22 JULY ,1991
C     REMARK : MAXIMUN RELATIVE ERROR LESS THAN 5.0% IF PT
C              LESS THAN 645.15D0 K
C
C     ----------------------------------------------------------
C
      FUNCTION THCONL(PT,N)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION A(5),B(4)
      DATA A/ 6.09937000D-01,-2.05611490D-03,-9.67565800D-06,
     &        3.94689000D-08,-3.17009600D-11/
      DATA B/ 1.00000000D000,-5.38878000D-03, 3.64155290D-06,
     &        5.99807000D-09/
	COMMON /SG9/TIME
      IF(N.NE.1) GOTO 5
      TSAT=TSATP(PT)-273.15D0
      GOTO 10
    5 TSAT=PT-273.15D0
   10 IF(TSAT.GE.350.0D0) GOTO 25
      FZ=0.0D0
      FM=0.0D0
      DO 15 I=1,5
   15 FZ=FZ+A(I)*TSAT**(I-1)
      DO 20 J=1,4
   20 FM=FM+B(J)*TSAT**(J-1)
      THCONL=DABS(FZ/FM)
      GOTO 30
   25 THCONL=0.435D0
      IF(TSAT.GT.372.0D0) WRITE(*,*)'IN THCONL: WARNNING! TSAT>645.15K'
   30 CONTINUE
      THCONL=THCONL
C      IF(TIME.GE.8.2D-2)THEN
C	  WRITE(*,*)'PT=',PT,' TSAT=',TSAT
C	  WRITE(*,*)'FZ=',FZ,' FM=',FM
C	ENDIF
	 
      RETURN
      END
C
C     ----------------------------------------------------------
C
      FUNCTION THCONL1(PT,N)
      IMPLICIT REAL(8) (A-H,O-Z)
C      DIMENSION A(5),B(4)
	COMMON /SG9/TIME
C
	IF(N.EQ.1)THEN
        T=TSATP(PT)
	ELSE
	  T=PT
	ENDIF
C
	IF(T.GE.643.15D0)T=643.15D0
	IF(T.LE.273.15D0)T=273.15D0
C     
      IF(T.GE.593.15D0)THEN
        THCONL1=-3.71954720928594D0+1.658526029220823D-002*T
	1          -1.595013313454172D-005*T**2
	ELSEIF(T.GE.553.15D0)THEN
        THCONL1=-1.17282457409001D0+7.987731314700212D-003*T
	1          -8.715410441243242D-006*T**2
	ELSE
        THCONL1=-3.47249539042321D0+3.446136984474435D-002*T 
	1          -1.073574155569477D-004*T**2
     2          +1.509429063610818D-007*T**3
     3          -8.234595285143100D-011*T**4
	ENDIF

      RETURN
      END
C     ----------------------------------------------------------
C
C     FUNCTION   VISSG(PT)
C
C     FUNCTION VISSG IS FOR DYNAMIC VISCOSITY OF SATURATION STEAM
C     AS A FUNCTION OF PRESSURE
C
C     VISSG   : DYNAMIC VISCOSITY OF SATURATION STEAM   (OUTPUT)
C               AS A FUNCTION OF PREAAURE(10-6 -KG/S.M)
C     PT      : PRESSURE (MPA)                           (INPUT)
C
C     REFERENCE : THE BWR PLANT ANALYZER . NUREG/CR-3943
C     LAST REVISION : 23 JULY ,1991
C     REMARK : MAXIMUN RELATIVE ERROR LESS THAN 3.0%
C
C     ----------------------------------------------------------
C
      FUNCTION VISSG(PT)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION A(3),G(4),F(4)
      DATA A/ 3.53D-08, 6.765D-11, 1.021D-14/
      DATA G/126.0,-1.60,4.8D-03,-4.7407D-06/
      DATA F/-2.885D-06,2.427D-08,-6.7893D-11,6.317037D-14/
      TSAT=TSATP(PT)-273.15D0
      VISSG0=0.407D-07*TSAT+8.04D-06
      IF(TSAT.GT.400.0D0) GOTO 10
      DENG=1.D0/SGSV(PT)
      VISSG=VISSG0-DENG*(1.858D-07-5.90D-10*TSAT)
      GOTO 30
   10 F1=0.D0
      F2=0.0D0
      F3=0.0D0
      PP=PT*1.0D+06
      DO 15 I=1,3
   15 F1=F1+A(I)*PP**(I-1)
      DO 20 J=1,4
   20 F2=F2+G(J)*TSAT**(J-1)
      DO 25 K=1,4
   25 F3=F3+F(K)*TSAT**(K-1)
      DENG=1.D0/SGSV(PT)
      VISSG=VISSG0+(F3+F1*F2)
   30 CONTINUE
         VISSG=VISSG*1.D+6
      RETURN
      END
C
C     ---------------------------------------------------
C
C     FUNCTION PRG(PP)
C
C     FUNCTION PRG RETURNS PRANDEL NUMBER OF SATURATED STEAM
C
C     PRG   : PRANDEL NUMBER OF SATURATED STEAM   (OUTPUT)
C     PP    : PRESSURE(MPA)                        (INPUT)
C
C     REMARK : MAXIMUN RELATIVE ERROR < 0.5%
C     LAST REVISION : JUNE 3 , 1992
C
C     ---------------------------------------------------
C
      FUNCTION PRG(PP)
      IMPLICIT REAL(8) (A-H,O-Z)
      P=PP
      IF(P-0.476D0) 10,20,30
  10  P=0.476D0
      GOTO 20
  30  IF(P-18.67D0) 20,20,40
  40  P=18.67D0
  20  PRG= .9985297910D+00
     &    +.1137216310D+00*P   -.1287209610D-01*P*P
     &    +.9052224340D-03*P**3-.6741919600D-05*P**4
      RETURN
      END
C
C     ----------------------------------------------------------
C     FUNCTION   THCONG(PT,N)
C
C     FUNCTION THCONG RETURNS THERMAL CONDUCTIVITY OF SATURATION
C     STEAM AS A FUNCTION OF PRESSURE
C
C     THCONG  : THERMAL CONDUCTIVITY OF SATURATION STEAM (OUTPUT)
C               AS A FUNCTION OF PRESSURE  (W/M.K)
C     PT      : PRESSURE (MPA)  , WHEN N=1 ;              (INPUT)
C               SATURATION TEMPERATURE   (K) , WHEN N=2 .
C
C     REFERENCE : THE BWR PLANT ANALYZER . NUREG/CR-3943
C     LAST REVISION : 23 JULY ,1991
C     REMARK : MAXIMUN RELATIVE ERROR LESS THAN 4.0% IF
C              PT < 633.15D0 K OR PT < 18.67 MPA
C
C     ----------------------------------------------------------
C
      FUNCTION THCONG(PT,N)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION A(4),B(3)
      DATA A/ 1.76D-02, 5.87D-05, 1.04D-07,-4.51D-11/
      DATA B/ 1.0351D-04, 4.198D-07,-2.771D-11/
      DATA C/ 2.1482D+05/
      IF(N.NE.1) GOTO 5
      TSAT=TSATP(PT)-273.15D0
      GOTO 10
    5 TSAT=PT-273.15D0
         PT=PSATT(PT)
   10 X1=0.0D0
      DO 15 I=1,4
   15 X1=X1+A(I)*TSAT**(I-1)
      X2=0.0D0
      DO 20 J=1,3
   20 X2=X2+B(J)*TSAT**(J-1)
      DG=1.D0/SGSV(PT)
      THCONG=X1+DG*(X2+C*DG*TSAT**(-4.2D0))
      THCONG=THCONG
      RETURN
      END
C
C     ----------------------------------------------------------
C
C     FUNCTION   VISCL(PT,N)
C
C     FUNCTION VISCL RETURNS DYNAMIC VISCOSITY OF SATURATION
C     WATER AS A FUNCTION OF PRESSURE OR SATURATION
C
C     VISCL   : DYNAMIC VISCOSITY OF SATURATION WATER AS (OUTPUT)
C               A FUNCTION OF PRESSURE OR SATURATION
C               TEMPRATURE      (10-6-KG/S.M)
C     PT      : PRESSURE (MPA)  , WHEN N=1 ;              (INPUT)
C               SATURATION TEMPERATURE   (K) , WHEN N=2 .
C
C     LAST REVISION : 2 JUNE ,1992
C     REMARK : MAXIMUN RELATIVE ERROR LESS THAN  2.4%
C
C     ----------------------------------------------------------
C
      FUNCTION VISCL(PT,N)
      IMPLICIT REAL(8) (A-H,O-Z)
      IF(N.NE.1) GOTO 5
      TSAT=TSATP(PT)-273.15D0
      GOTO 10
    5 TSAT=PT-273.15D0
   10 CONTINUE
      VISCL =  .7805527340D+03
     &        -.8066720960D+01*TSAT**1 +.3920077900D-01*TSAT**2
     &        -.9558879530D-04*TSAT**3 +.1125612520D-06*TSAT**4
     &        -.5448403880D-10*TSAT**5
      VISCL=VISCL
      RETURN
      END
C
C     ---------------------------------------------------
C
C     FUNCTION VISC(TT,RHO)
C
C     FUNCTION VISC COMPUTE THE VISCOSITY OF SUBCOOLED AND
C     SATURATED LIQUID AND SATURATED AND SUPERHEATED STEAM
C
C     TT   :  FLUID TEMPERATURE (K)             (INPUT)
C     RO   :  DENSITY           (KG/M3)         (INPUT)
C     VISC :  FLUID VISCOSITY (10-6 PA-SEC)     (OUTPUT)
C
C     THE EQUATION IN THIS FUNCTION IS BASED ON THE EQUATION
C     GIVEN BY UNKNOWN, 'NEW VALUES FOR THE VISCOSITY OF
C     WATER SUBSTANCE',
C
C     MECHANICAL ENGINEERING, JULY, 1976
C     (COPY FROM RETRAN-02)
C
C     ---------------------------------------------------
C
      FUNCTION VISC(T,RHO)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION B(5,6),A(4)
      DATA FACT1/647.27D0/,FACR1/0.003147D0/
      DATA A/0.0181583D0,0.0177624D0,0.0105287D0,-0.0036744D0/
      DATA B/0.501938D0,0.235622D0,-0.274637D0,0.145831D0,-0.0270448D0,
     &      0.162888D0,0.789393D0,-0.743539D0,0.263129D0,-0.0253093D0,
     &     -0.130356D0,0.673665D0,-0.959456D0,0.347247D0,-0.0267758D0,
     &      0.907919D0,1.207552D0,-0.687343D0,0.213486D0,-0.0822904D0,
     &     -0.551119D0,0.0670665D0,-0.497089D0,0.100754D0,0.0602253D0,
     &     0.146543D0,-0.084337D0,0.195286D0,-0.032932D0,-0.0202595D0/
      TH=T/FACT1
      C=TH**0.5D0
      D=0.D0
      DO 10 K=1,4
10    D=D+A(K)*(1.D0/TH)**(K-1)
      ETAO=C/D
      SUM=0.D0
      DO 20 I=1,6
      DO 20 J=1,5
20    SUM=SUM+B(J,I)*(1.D0/TH-1.D0)**(I-1)*(RHO*FACR1-1.D0)
     &   **(J-1)
      VISC=ETAO*DEXP(RHO*FACR1*SUM)
      RETURN
      END
C
C     ----------------------------------------------------------------
C
C       FUNCTION ZPK(TT)
C
C       FUNCTION ZPK RETURNS SATURATED PRESSURE
C
C       ZPK   : SATURATED PRESSURE (MPa)      (OUTPUT)
C       TT    : TEMPRATURE (K)                 (INPUT)
C
C     -----------------------------------------------------------------
C
      FUNCTION ZPK(TT)
      IMPLICIT REAL(8) (A-H,O-Z)
      REAL(8) KK
      DIMENSION F(10)
      DATA F/0.0D0,-7.691234564D0,-2.608023696D1,-1.681706546D2,
     1 6.423285504D1,-1.189646225D2,4.16711732D0,
     1 2.09750676D1,1.0D9,6.0D0/
      ZT=TT/647.3D0
      ZTT=1.D0-ZT
      KK=0.0D0
      DO 10 II=1,6
      I=7-II
10    KK=KK*ZTT+F(I)
      KK=KK/ZT
      KK=KK/(1.D0+F(7)*ZTT+F(8)*ZTT**2)
      KK=KK-ZTT/(F(9)*ZTT**2+F(10))
      ZPK=DEXP(KK)*221.2D0
      RETURN
      END
C
C     ----------------------------------------------------------------
C
C       FUNCTION ZTK(P1)
C
C       FUNCTION ZTK RETURNS SATURATED TEMPRATURE
C
C       ZTK   : SATURATED TEMPRATURE (K)         (OUTPUT)
C       P1    : PRESSURE (MPa)                    (INPUT)
C
C     -----------------------------------------------------------------
C
      FUNCTION ZTK(P1)
      IMPLICIT REAL(8) (A-H,O-Z)
      PP=P1*10.0D0
      TA=100.0D0*PP**0.25D0
11    PB=ZPK(TA+273.15D0)
      IF(DABS((PP-PB)/PP).LT.0.000005D0) GOTO 12
      TA=TA+25.0D0*(PP-PB)/PB**7.5D-1
      GOTO 11
12    ZTK=TA+273.15D0
      RETURN
      END
C
C     -------------------------------------------------------------
C       FUNCTION  HFPT(P1,T1)
C
C       FUNCTION HFPT SUBCOOLED AND SATURATED WATER ENTHALPY
C
C       HGPT  : SUBCOOLED AND SATURATED WATER ENTHALPY(KJ/KG) (OUTPUT)
C       P1    : PRESSURE (MPa)                                 (INPUT)
C       T1    : TEMPRATURE (K)                                 (INPUT)
C
C     -------------------------------------------------------------
C
        FUNCTION HFPT(P1,T1)
      IMPLICIT REAL(8) (A-H,O-Z)
        DIMENSION A(23),E(12)
        DATA A/6.824687741D3,-5.422063673D2,-2.096666205D4,
     *3.941286787D4,-6.733277739D4,9.902381028D4,-1.093911774D5,
     *8.590841667D4,-4.511168742D4,1.418138926D4,-2.017271113D3,
     *7.982692717D0,-2.616571843D-2,1.52241179D-3,2.284279054D-2,
     *2.421647003D2,1.269716088D-10,2.074838328D-7,2.174020350D-8,
     *1.105710498D-9,1.293441934D1,1.308119072D-5,6.047626338D-14/
        DATA E/8.438375405D-1,5.362162162D-4,1.72000D0,7.342278489D-2,
     *4.975858870D-2,6.53715430D-1,1.15D-6,1.5108D-5,1.4188D-1,
     *7.002753165D0,2.995284926D-4,2.04D-1/
      PP=P1*10.0D0
      TT=T1-273.15D0
        ZP=PP/221.2D0
        ZT=(273.15D0+TT)/647.3D0
        Y=1.D0-E(1)*ZT**2-E(2)/ZT**6
        Z=Y+DSQRT(E(3)*Y**2-2.0D0*E(4)*ZT+2.0D0*E(5)*ZP)
        YP=6.0D0*E(2)/ZT**7-2.0D0*E(1)*ZT
        ZH=0.0D0
        DO 5 II=1,10
        I=11-II
5       ZH=ZH*ZT+(I-2)*A(I+1)
        ZH=A(1)*ZT-ZH+A(12)*(Z*(17.0D0*(Z/29.0D0-Y/12.0D0)
	1  +5.0D0*ZT*YP/12.0D0)
     *+E(4)*ZT-(E(3)-1.D0)*ZT*Y*YP)/Z**0.2941176D0
        ZH=ZH+(A(13)-A(15)*ZT**2+A(16)*(9.0D0*ZT+E(6))*(E(6)
     *-ZT)**9+A(17)*(20.0D0*ZT**19+E(7))/(E(7)+ZT**19)**2)*ZP
        ZH=ZH-(12.0D0*ZT**11+E(8))*(A(18)*ZP+A(19)*ZP**2
     *+A(20)*ZP**3)/(E(8)+ZT**11)**2
        ZH=ZH+A(21)*ZT**18*(17.0D0*E(9)+19.0D0*ZT**2)*(1.D0/(E(10)+ZP)
     ***3+E(11)*ZP)
        ZH=ZH+A(22)*E(12)*ZP**3+21.D0*A(23)*ZP**4/ZT**20
        HFPT=ZH*70.1204D0
        RETURN
        END
C
C     -------------------------------------------------------------
C       FUNCTION  HGPT(P1,T1)
C
C       FUNCTION HGPT RETURNS SATURATED AND SUPER-HEATED
C       STEAM ENTHALPY
C
C       HGPT   : SATURATED AND SUPER-HEATED STEAM ENTHALPY (OUTPUT)
C                (KJ/KG)
C       P1     : PRESSURE (MPa)                             (INPUT)
C       T1     : TEMPRATURE (K)                             (INPUT)
C
C     -------------------------------------------------------------
C
        FUNCTION HGPT(P1,T1)
      IMPLICIT REAL(8) (A-H,O-Z)
        REAL(8) L0,L1,L2
        DIMENSION BUV(8,3),ZUV(8,3),B0(5),BUL(3,3),XUL(3,3),B9(7)
        DATA B0/2.856067796D1,-5.438923329D1,4.330662834D-1,
     *-6.547711697D-1,8.565182058D-2/
        DATA B9/1.936587558D2,-1.388522425D3,4.126607219D3,
     *-6.508211677D3,5.745984054D3,-2.693088365D3,5.235718623D2/
        DATA BUV/6.670375918D-2,8.390104328D-2,4.520918904D-1,
     *-5.975336707D-1,5.958051609D-1,1.190610271D-1,1.683998803D-1,
     *6.552390126D-3,1.3889883801D0,2.614670893D-2,1.069036614D-1,
     *-8.847535804D-2,-5.159303373D-1,-9.867174132D-2,
     *-5.809438001D-2,5.710218649D-4,0.0D0,-3.373439453D-2,0.0D0,
     *0.0D0,2.075021122D-1,0.0D0,0.0D0,0.0D0/
        DATA ZUV/13.0D0,1.8D1,1.8D1,2.5D1,3.2D1,12.0D0,2.4D1,2.4D1,
     *3.0D0,2.0D0,1.0D1,1.4D1,2.8D1,1.1D1,1.8D1,1.4D1,
     *0.0D0,1.0D0,0.0D0,0.0D0,2.4D1,0.0D0,0.0D0,0.0D0/
      DATA XUL/1.4D1,1.9D1,5.4D1,0.0D0,0.0D0,2.7D1,0.0D0,0.0D0,0.0D0/
        DATA BUL/4.006073948D-1,8.636081627D-2,-8.532322921D-1,
     *0.0D0,0.0D0,3.460208861D-1,0.0D0,0.0D0,0.0D0/
        DATA L0,L1,L2,B,B1,C1/1.574373327D1,-34.17061978D0,
     *1.931380707D1,7.633333333D-1,1.683599274D1,4.260321148D0/
        PP=P1*10.0D0
        TT=T1-273.15D0
        ZP=PP/221.2D0
        ZT=(273.15D0+TT)/647.3D0
        BL=L0+L1*ZT+L2*ZT**2
        BLP=L1+2.0D0*L2*ZT
        X=DEXP(B*(1.D0-ZT))
        Z0=0.0D0
        DO 70 II=1,5
        I=6-II
70      Z0=Z0*ZT+(I-2)*B0(I)
        Z1=0.0D0
        DO 80 I=1,5
        Z2=0.0D0
        DO 90 J=1,3
90      Z2=Z2+BUV(I,J)*(1.D0+B*ZUV(I,J)*ZT)*X**ZUV(I,J)
80      Z1=Z1+Z2*ZP**I
        Z5=0.0D0
        IF(ZP.LT.0.005) GOTO 100
        DO 1030 I=1,3
        Z4=0.0D0
        DO 1020 J=1,3
        Z2=0.0D0
        Z3=0.0D0
        DO 1025 K=1,3
        Z2=Z2+XUL(I,K)*BUL(I,K)*X**XUL(I,K)
1025    Z3=Z3+BUL(I,K)*X**XUL(I,K)
        Z4=Z4+BUV(I+5,J)*X**ZUV(I+5,J)*((1.D0+ZUV(I+5,J)*B*ZT)
     *-B*ZT*Z2/(Z3+1.D0/ZP**(I+3)))
1020    CONTINUE
1030    Z5=Z5+Z4/(Z3+1.D0/ZP**(I+3))
100     Z6=0.0D0
        IF(ZP.LT..1D0) GOTO 1100
        DO 1110 II=1,7
        I=8-II
1110    Z6=Z6*X+B9(I)*(1.D0+ZT*(10.0D0*BLP/BL+B*(I-1)))
        Z6=Z6*ZP*(ZP/BL)**10
1100    ZH=B1*ZT-Z0-Z1-Z5+Z6
        HGPT=ZH*70.1204D0
        RETURN
        END
C
C     -------------------------------------------------------------
C       FUNCTION  SWHP(PP)
C
C       FUNCTION SWHP IS FOR SATURATED WATER ENTHALPY
C
C       SWHP   : SATURATED WATER ENTHALPY (KJ/KG)   (OUTPUT)
C       PP     : PRESSURE (MPa)                      (INPUT)
C
C     -------------------------------------------------------------
C
        FUNCTION SWHP(PP)
        IMPLICIT REAL(8) (A-H,O-Z)
        TS=ZTK(PP)
        SWHP=HFPT(PP,TS)
        RETURN
        END
C
C     -------------------------------------------------------------
C       FUNCTION  SSHP(PP)
C
C       FUNCTION SSHP IS FOR SATURATED STEAM ENTHALPY
C
C       SSHP   : SATURATED STEAM ENTHALPY (KJ/KG)  (OUTPUT)
C       PP     : PRESSURE (MPa)                     (INPUT)
C
C     -------------------------------------------------------------
C
      FUNCTION SSHP(PP)
      IMPLICIT REAL(8) (A-H,O-Z)
      TS=ZTK(PP)
      SSHP=HGPT(PP,TS)
      RETURN
      END
C
C     -------------------------------------------------------------
C       FUNCTION PRC(P,T)
C
C       FUNCTION PRC RETURNS LIQUID PRANDTL NUMBER
C
C       PRC   : LIQUID PRANDTL NUMBER          (OUTPUT)
C       P     : PRESSURE (MPa)                  (INPUT)
C       T     : TEMPRATURE(K)                   (INPUT)
C
C     -------------------------------------------------------------
C
      FUNCTION PRC(P,T)
      IMPLICIT REAL(8) (A-H,O-Z)
      CP=(HLPT(P,T)-HLPT(P,T-2.0D0))/2.0D0
      RO=1.D0/SVPH(P,HLPT(P,T))
      PRC=VISC(T,RO)*1.0D-6*CP*1000.0D0/THCOND(T,RO)
      RETURN
      END
C
C     -------------------------------------------------------------
C
C     FUNCTION  SIGMA(T)
C
C     FUNCTION SIGMA SURFACE TENSION BETWEEN WATER AND VAPOR
C
C     SIGMA  : SURFACE TENSION BETWEEN WATER AND VAPOR(N/M) (OUTPUT)
C     T      : TEMPERATURE(K)                                (INPUT)
C
C     LAST REVISION   JULY 9 , 1991
C
C     -------------------------------------------------------------
C
      FUNCTION SIGMA(T)
      IMPLICIT REAL(8) (A-H,O-Z)
      DIMENSION A(5)
      DATA A/1.160936807D-1,1.121404688D-3,-5.752805180D-6,
     &                      1.286274650D-8,-1.149719290D-11/
      DATA TK,BETA/647.3D0,0.83D0/
      SUM=0.D0
      DO 1 I=2,5
    1 SUM=SUM+A(I)*(TK-T)**I
      SIGMA=SUM+A(1)*(TK-T)**2/(1.D0+BETA*(TK-T))
      SIGMA=SIGMA*1.D-3
      RETURN
      END
C
C     ---------------------------------------------------
C
C     FUNCTION TCUO2(TT)
C
C     FUNCTION TCUO2 RETURNS THERMAL CONDUCTIVITY OF UO2
C
C     TCUO2   : THERMAL CONDUCTIVITY OF UO2(W/K.M)  (OUTPUT)
C     TT      :  TEMPERATURE (K)                     (INPUT)
C
C     LAST REVISION : APRIL 28 , 1992
C     MAXIMUN RELATIVE ERROR < 3.8%
C
C     ---------------------------------------------------
C
      FUNCTION TCUO2(TT)
      IMPLICIT REAL(8) (A-H,O-Z)
      T=TT
      IF(T-273.15D0) 10,20,30
  10  T=273.15D0
      GOTO 20
  30  IF(T-3088.71D0) 20,20,40
  40  T=3088.71D0
  20  TC= .1288788790D+02
     &   -.2113481610D-01*T+.1830378640D-04*T*T
     &   -.8458677310D-08*T**3+.2008269200D-11*T**4
     &   -.1877097800D-15*T**5
      TCUO2=TC
      RETURN
      END
C
C     ---------------------------------------------------
C
C     FUNCTION TCNUO2(TT)
C
C     FUNCTION TCNUO2 RETURNS THERMAL CONDUCTIVITY OF UO2
C
C     TCNUO2   : THERMAL CONDUCTIVITY OF UO2( W/K.M)  (OUTPUT)
C     TT       : TEMPERATURE    (K)                    (INPUT)
C
C     LAST REVISION : APRIL 28 , 1992
C     MAXIMUN RELATIVE ERROR < 3.8%
C
C     ---------------------------------------------------
C
      FUNCTION TCNUO2(TT)
      IMPLICIT REAL(8) (A-H,O-Z)
      T=TT-273.15D0
      TC=38.24D0/(T+402.55D0)+4.788D-13*TT*TT*TT
      TCNUO2=TC*100.0D0
      RETURN
      END
C
C     ---------------------------------------------------
C
C     FUNCTION TCCLA(TT)
C
C     FUNCTION TCCLA IS FOR THERMAL CONDUCTIVITY OF CLAD
C
C     TCCLA    : THERMAL CONDUCTIVITY OF CLAD( W/K.M)  (OUTPUT)
C     TT       : TEMPERATURE    (K)                      (INPUT)
C
C     LAST REVISION : APRIL 28 , 1992
C     MAXIMUN RELATIVE ERROR < 4.4%
C
C     ---------------------------------------------------
C
      FUNCTION TCCLA(TT)
      IMPLICIT REAL(8) (A-H,O-Z)
      T=TT
      IF(T-273.15D0) 10,20,30
  10  T=273.15D0
      GOTO 20
  30  IF(T-2122.04D0) 20,20,40
  40  T=2122.04D0
  20  TC= .1163511080D+02
     &   +.4135427060D-02*T+.5157591650D-05*T*T
     &   -.9851657260D-09*T**3+.1153472280D-11*T**4
      TCCLA=TC
      RETURN
      END
C
C     ---------------------------------------------------
C     FUNCTION TCINC(TT)
C
C     FUNCTION TCINC IS FOR THERMAL CONDUCTIVITY OF INCONEL
C
C     TCINC : THERMAL CONDUCTIVITY OF INCONEL(W/K.M)  (OUTPUT)
C     TT    : TEMPERATURE    (K)                       (INPUT)
C
C     LAST REVISION : APRIL 28 , 1992
C     MAXIMUN RELATIVE ERROR < 2.1%
C
C     ---------------------------------------------------
C
      FUNCTION TCINC(TT)
      IMPLICIT REAL(8) (A-H,O-Z)
      T=TT
      IF(T-373.15D0) 10,20,30
  10  T=373.15D0
      GOTO 20
  30  IF(T-1273.15D0) 20,20,40
  40  T=1273.15D0
  20  TC= .1866690820D+02
     &   -.4085551940D-02*T   -.1036019960D-04*T*T
     &   +.3841511640D-07*T**3-.1706647670D-10*T**4
      TCINC=TC
      RETURN
      END
C
C     ------------------------------------------------------
C     FUNCTION TCINC8(TT)
C
C     FUNCTION TCINC8 THERMAL CONDUCTIVITY OF INCOLOY800
C
C     TCINC8 : THERMAL CONDUCTIVITY OF INCOLOY800(W/K.M) (OUTPUT)
C     TT     : TEMPERATURE    (K)  (360-1173K)            (INPUT)
C
C     LAST REVISION : APRIL 30 , 1993
C     MAXIMUN RELATIVE ERROR < 4.4%
C
C     ------------------------------------------------------
C
         FUNCTION TCINC8(TT)
      IMPLICIT REAL(8) (A-H,O-Z)
         T=TT-273.15D0
         TC=11.5D0+0.016D0*T
         TCINC8=TC
         RETURN
         END
C
C     ------------------------------------------------------
C
C     FUNCTION TCINC6(TT)
C
C     FUNCTION TCINC6 THERMAL CONDUCTIVITY OF INCONEL600
C
C     TCINC6 : THERMAL CONDUCTIVITY OF INCONEL600(W/K.M) (OUTPUT)
C     TT     : TEMPERATURE(K)  (360-1173K)                (INPUT)
C
C     LAST REVISION : APRIL 30 , 1993
C     MAXIMUN RELATIVE ERROR < 4.4%
C
C     ------------------------------------------------------
C
         FUNCTION TCINC6(TT)
      IMPLICIT REAL(8) (A-H,O-Z)
         T=TT-273.15D0
         TC=13.89D0+0.017D0*T
         TCINC6=TC
         RETURN
         END
C
C     ------------------------------------------------------
C     FUNCTION TCCLA4(TT)
C
C     FUNCTION TCCLA4 RETURNS THERMAL CONDUCTIVITY OF CLAD(Kr-4)
C
C     TCCLA4 : THERMAL CONDUCTIVITY OF CLAD(Kr-4)(W/K.M)  (OUTPUT)
C     TT     : TEMPERATURE(K)  (360-1123K)                 (INPUT)
C
C     LAST REVISION : APRIL 30 , 1993
C     MAXIMUN RELATIVE ERROR < 4.4%
C
C     ------------------------------------------------------
C
         FUNCTION TCCLA4(TT)
      IMPLICIT REAL(8) (A-H,O-Z)
         T=TT-273.15D0
         TC=7.73D0+3.15D-2*T-2.87D-5*T*T+1.552D-8*T*T*T
         TCCLA4=TC
         RETURN
         END
C
C     ------------------------------------------------------
C
C     FUNCTION TCCLA2(TT)
C
C     FUNCTION TCCLA2 RETURNS THERMAL CONDUCTIVITY OF CLAD(Kr-2)
C
C     TCCLA2  : THERMAL CONDUCTIVITY OF CLAD(Kr-2)(W/K.M)  (OUTPUT)
C     TT      : TEMPERATURE    (K)  (360-1123K)             (INPUT)
C
C     LAST REVISION : APRIL 30 , 1993
C     MAXIMUN RELATIVE ERROR < 4.4%
C
C     ------------------------------------------------------
C
      FUNCTION TCCLA2(TT)
      IMPLICIT REAL(8) (A-H,O-Z)
      T=TT-273.15D0
      TC=95.D0*(0.17D0+1.04D-4*T+1.08D-7*T*T)
      TCCLA2=TC
      RETURN
      END
C
C     ---------------------------------------------------
C
C     FUNCTION VHCUO2(TT)
C
C     FUNCTION VHCUO2 RETURNS VOLU. HEAT CAPACITY OF UO2
C
C     VHCUO2  : VOLU. HEAT CAPACITY OF UO2(KJ/K.M^3)    (OUTPUT)
C      TT     : TEMPERATURE   (K)                        (INPUT)
C
C     LAST REVISION : APRIL 28 , 1992
C     MAXIMUN RELATIVE ERROR < 2.0%
C
C     ---------------------------------------------------
C
      FUNCTION VHCUO2(TT)
      IMPLICIT REAL(8) (A-H,O-Z)
      T=TT
      IF(T-373.15D0) 10,20,30
  10  T=373.15D0
      GOTO 20
  30  IF(T-2373.15D0) 20,20,40
  40  IF(T-3113.15D0) 50,50,60
  60  VHC=6802.5527D0
      GOTO 70
  20  VHC= .1472484980D+04
     &    +.4770160190D+01*T-.4261192860D-02*T*T
     &    +.1457196500D-05*T**3-.1090418860D-09*T**4
      GOTO 70
  50  VHC=-.8826350000D+05
     &    +.9059268180D+02*T-.2918222540D-01*T*T
     &    +.3178810400D-05*T**3
  70  VHCUO2=VHC
      RETURN
      END
C
C     ---------------------------------------------------
C
C     FUNCTION CPUO2(TT)
C
C     FUNCTION CPUO2 IS FOR PRES. HEAT CAPACITY OF UO2
C
C     CPUO2 : PRES. HEAT CAPACITY OF UO2(KJ/K.KG )  (OUTPUT)
C     TT    : TEMPERATURE   (K) (300-3073K)          (INPUT)
C
C     LAST REVISION : APRIL 30 , 1993
C     MAXIMUN RELATIVE ERROR < 2.0%
C
C     ---------------------------------------------------
C
      FUNCTION CPUO2(TT)
      IMPLICIT REAL(8) (A-H,O-Z)
      T=TT-273.15D0
      IF(T-25.0D0) 10,20,30
  10  T=25.0D0
      GOTO 20
  30  IF(T-1226.0D0) 20,20,40
  40  IF(T-2800.0D0) 50,50,60
  60  T=2800.0D0
      GOTO 50
  20  VHC=304.38D0+2.51D-2*T-6.0D+6/((T+273.15D0)*(T+273.15D0))
      GOTO 70
  50  VHC=-712.25D0+2.789D0*T-2.71D-3*T*T+1.12D-6*T*T*T
     &    -1.59D-10*T*T*T*T
  70  CPUO2=VHC*1.D-3
      RETURN
      END
C
C     ---------------------------------------------------
C
C     FUNCTION VHCCLA(TT)
C
C     FUNCTION VHCCLA IS FOR VOLU. HEAT CAPACITY ODF CLAD
C
C     VHCCLA  : VOLU. HEAT CAPACITY ODF CLAD(KJ/K.M^3)  (OUTPUT)
C     TT      : TEMPERATURE   (K)                        (INPUT)
C
C     LAST REVISION : APRIL 28 , 1992
C     MAXIMUN RELATIVE ERROR < 1.2%
C
C     ---------------------------------------------------
C
      FUNCTION VHCCLA(TT)
      IMPLICIT REAL(8) (A-H,O-Z)
      T=TT
      IF(T-255.37D0) 10,20,30
  10  T=255.37D0
      GOTO 20
  30  IF(T-1248.43D0) 20,20,40
  40  VHC=2312.183D0
      GOTO 50
  20  VHC= .1640218010D+04
     &    +.1156234380D+01*T-.4841190520D-03*T*T
  50  VHCCLA=VHC
      RETURN
      END
C
C     ---------------------------------------------------
C
C     FUNCTION CPCLA2(TT)
C
C     FUNCTION CPCLA2 IS FOR PRES. HEAT CAPACITY ODF CLAD(Kr2)
C
C     CPCLA2  : PRES. HEAT CAPACITY ODF CLAD(Kr2)(KJ/K.KG) (OUTPUT)
C     TT      : TEMPERATURE   (K) (273-1320K)               (INPUT)
C
C     LAST REVISION : APRIL 30 , 1993
C     MAXIMUN RELATIVE ERROR < 1.2%
C
C     ---------------------------------------------------
C
      FUNCTION CPCLA2(TT)
      IMPLICIT REAL(8) (A-H,O-Z)
      T=TT-273.15D0
      IF(T-0.0) 10,10,30
  10  T=0.0D0
      GOTO 20
  30  IF(T-633.0D0) 20,20,40
  40  IF(T-900.0D0) 60,60,70
  70  IF(T-1050.0D0) 80,80,90
  90  T=1050.0D0
         GOTO 80
  20  VHC=285.0D0+9994.7D-5*(1.8D0*T+32.0D0)
         GOTO 50
  60  VHC=359.6D0+9994.7D-5*(1.8D0*T+32.0D0)
         GOTO 50
  80  VHC=357.9D0+9994.7D-5*(1.8D0*T+32.0D0)
  50  CPCLA2=VHC*1.D-3
      RETURN
      END
C
C     ---------------------------------------------------
C
C     FUNCTION CPCLA4(TT)
C
C     FUNCTION CPCLA4 IS FOR PRES. HEAT CAPACITY OF CLAD (Kr4)
C
C     CPCLA4 : PRES. HEAT CAPACITY OF CLAD (Kr4)(KJ/KG.K)  (OUTPUT)
C     TT     : TEMPERATURE   (K)                            (INPUT)
C
C     LAST REVISION : APRIL 30 , 1993
C     MAXIMUN RELATIVE ERROR < 1.2%
C
C     ---------------------------------------------------
C
      FUNCTION CPCLA4(TT)
      IMPLICIT REAL(8) (A-H,O-Z)
      T=TT-273.15D0
      IF(T-0.0D0) 10,10,30
  10  T=0.0D0
      GOTO 20
  30  IF(T-750.0D0) 20,20,40
  40  VHC=360.0D0
      GOTO 50
  20  VHC=286.5D0+.1D0*T
  50  CPCLA4=VHC*1.D-3
      RETURN
      END
C
C     ---------------------------------------------------
C
C     FUNCTION VHCINC(TT)
C
C     FUNCTION VHCINC RETURNS VOLU. HEAT CAPACITY OF INCONEL
C
C     VHCINC  : VOLU. HEAT CAPACITY OF INCONEL(KJ/K.M^3) (OUTPUT)
C     TT      : TEMPERATURE   (K)                         (INPUT)
C
C     LAST REVISION : APRIL 28 , 1992
C     MAXIMUN RELATIVE ERROR < 0.72%
C
C     ---------------------------------------------------
C
      FUNCTION VHCINC(TT)
      IMPLICIT REAL(8) (A-H,O-Z)
      T=TT
      IF(T-366.48D0) 10,20,30
  10  T=366.48D0
      GOTO 20
  30  IF(T-810.92D0) 20,20,40
  40  T=810.92D0
  20  VHC= .3352891840D+04
     &    +.1342795490D+01*T+.8341889770D-04*T*T
      VHCINC=VHC
      RETURN
      END
C
      FUNCTION COEXPAN(T)
      IMPLICIT REAL(8) (A-H,O-Z)
C      DIMENSION A(5),B(4)
C
	IF(T.GE.643.15D0)T=643.15D0
	IF(T.LE.273.15D0)T=273.15D0
C     
      IF(T.LE.383.15D0)THEN
	  COEXPAN=-262.116011700509200D0+2.182947516426936D0*T
	1        -6.052288896101953D-003*T**2
	1        +5.730030497013489D-006*T**3
	ELSEIF(T.LE.493.15D0)THEN
        COEXPAN=9.861355244750287D0-5.838389569095244D-002*T
	1         +1.399628494641906D-004*T**2
	ELSE
        COEXPAN=-240.3834905199119D0+2.094546524020144D0*T
	1         -5.434267699611387D-003*T**2
     1         +4.540807585537368D-006*T**3
	ENDIF
	COEXPAN=COEXPAN*1.D-4
C
      RETURN
      END
C     ----------------------------------------------------------


      BLOCK DATA
      IMPLICIT REAL(8) (A-H,O-Z)
	DIMENSION CT1(2,4),CT2(5,5),CT3(5,5),CT4(5,5)
      DIMENSION CN1(3,5),CN2(4,3),CG1(3,5),CG2(4,3)
      DIMENSION CSF1(9),CSF2(9),CSF3(9),CSG1(12),CSG2(9),CSG3(7)
      COMMON /CCC1/CT1,CT2,CT3,CT4
      COMMON /CCC2/CN1
      COMMON /CCC3/CN2
      COMMON /CCC4/CG1,CG2
      COMMON /CCC5/CSF1,CSF2,CSF3
      COMMON /CCC6/CSG1,CSG2,CSG3
      DATA CT1(1,1),CT1(1,2),CT1(1,3),CT1(1,4)/
     &          0.3276275552D002,0.9763617000D000,
     &          0.1857226027D-03,-.4682674330D-06/
      DATA CT1(2,1),CT1(2,2),CT1(2,3),CT1(2,4)/
     &          0.3360880214D-02,-.5595281760D-04,
     &          0.1618595991D-06,-.1180204381D-09/
      DATA CT2(1,1),CT2(1,2),CT2(1,3),CT2(1,4),CT2(1,5)/
     &          0.6390801208D003,-.3055217235D001,0.8713231868D-02,
     &          -.6269403683D-05,-.9844700000D-11/
      DATA CT2(2,1),CT2(2,2),CT2(2,3),CT2(2,4),CT2(2,5)/
     &          -.4302857237D000,0.2673303422D-02,-.5198380474D-05,
     &          0.3037825558D-08,0.3309704045D-12/
      DATA CT2(3,1),CT2(3,2),CT2(3,3),CT2(3,4),CT2(3,5)/
     &          0.1174524584D-03,-.6839200986D-06,0.1168011772D-08,
     &          -.4260074181D-12,-.2732087163D-15/
      DATA CT2(4,1),CT2(4,2),CT2(4,3),CT2(4,4),CT2(4,5)/
     &          -.1473977290D-07,0.8018858166D-10,-.1164901355D-12,
     &          0.4559400289D-17,0.5825511142D-19/
      DATA CT2(5,1),CT2(5,2),CT2(5,3),CT2(5,4),CT2(5,5)/
     &          0.7104327342D-12,-.3649500626D-14,0.4457387575D-17,
     &          0.1678398723D-20,-.3756804091D-23/
      DATA CT3(1,1),CT3(1,2),CT3(1,3),CT3(1,4),CT3(1,5)/
     &          -.1179100862D005,0.2829274345D002,-.2678181564D-01,
     &          0.1218742752D-04,-.2092033147D-08/
      DATA CT3(2,1),CT3(2,2),CT3(2,3),CT3(2,4),CT3(2,5)/
     &          0.1256160907D003,-.3333448495D000,0.3326901268D-03,
     &          -.1477890326D-06,0.2463258371D-10/
      DATA CT3(3,1),CT3(3,2),CT3(3,3),CT3(3,4),CT3(3,5)/
     &          -.1083713369D000,0.2928177730D-03,-.2972436458D-06,
     &          0.1342639113D-09,-.2275585718D-13/
      DATA CT3(4,1),CT3(4,2),CT3(4,3),CT3(4,4),CT3(4,5)/
     &          0.3278071846D-04,-.8970959364D-07,0.9246248312D-10,
     &          -.4249155515D-13,0.7338316751D-17/
      DATA CT3(5,1),CT3(5,2),CT3(5,3),CT3(5,4),CT3(5,5)/
     &          -.3425564927D-08,0.9527692453D-11,-.1001409043D-13,
     &          0.4703914404D-17,-.8315044742D-21/
      DATA CT4(1,1),CT4(1,2),CT4(1,3),CT4(1,4),CT4(1,5)/
     &          0.3795256853D004,-.6347031007D001,0.2867228326D-02,
     &          0.5953599813D-08,0.4798207438D-10/
      DATA CT4(2,1),CT4(2,2),CT4(2,3),CT4(2,4),CT4(2,5)/
     &          -.3910086240D001,0.1222747819D-01,-.1404664699D-04,
     &          0.7505679464D-08,-.1608693653D-11/
      DATA CT4(3,1),CT4(3,2),CT4(3,3),CT4(3,4),CT4(3,5)/
     &          0.3410500159D-04,0.7010900113D-09,-.1030201866D-09,
     &          0.5731099333D-14,0.3720795449D-16/
      DATA CT4(4,1),CT4(4,2),CT4(4,3),CT4(4,4),CT4(4,5)/
     &          0.1527377542D-06,-.5356866315D-09,0.6823225984D-12,
     &          -.3668096142D-15,0.6946004624D-19/
      DATA CT4(5,1),CT4(5,2),CT4(5,3),CT4(5,4),CT4(5,5)/
     &          -.1437179752D-10,0.5006731336D-13,-.6365519546D-16,
     &          0.3473711350D-19,-.6842306083D-23/
      DATA CN1(1,1),CN1(1,2),CN1(1,3),CN1(1,4),CN1(1,5)/
     &          -.4117961750D001,-.3811924543D-03,0.4308265942D-05,
     &          -.9160120130D-08,0.8017924673D-11/
      DATA CN1(2,1),CN1(2,2),CN1(2,3),CN1(2,4),CN1(2,5)/
     &          -.4816067020D-05,0.7744786733D-07,-.6988467605D-09,
     &          0.1916720525D-11,-.1760288590D-14/
      DATA CN1(3,1),CN1(3,2),CN1(3,3),CN1(3,4),CN1(3,5)/
     &          -.1820625039D-08,0.1440785930D-10,-.2082170753D-13,
     &          -.3603625114D-16,0.7407124321D-19/
      DATA CN2(1,1),CN2(1,2),CN2(1,3)/
     &     -.1403086182D004,0.1802594763D001,-.2097279215D-03/
      DATA CN2(2,1),CN2(2,2),CN2(2,3)/
     &     0.3817195017D000,-.5394444747D-03,0.1855203702D-06/
      DATA CN2(3,1),CN2(3,2),CN2(3,3)/
     &     -.6449501159D-04,0.8437637660D-07,-.2713755001D-10/
      DATA CN2(4,1),CN2(4,2),CN2(4,3)/
     &     0.7823817858D-08,-.1053834646D-10,0.3629590764D-14/
      DATA CG1(1,1),CG1(1,2),CG1(1,3),CG1(1,4),CG1(1,5)/
     &          -.4117961750D001,-.3811924543D-03,0.4308265942D-05,
     &          -.9160120130D-08,0.8017924673D-11/
      DATA CG1(2,1),CG1(2,2),CG1(2,3),CG1(2,4),CG1(2,5)/
     &          -.4816067020D-05,0.7744786733D-07,-.6988467605D-09,
     &          0.1916720525D-11,-.1760288590D-14/
      DATA CG1(3,1),CG1(3,2),CG1(3,3),CG1(3,4),CG1(3,5)/
     &          -.1820625039D-08,0.1440785930D-10,-.2082170753D-13,
     &          -.3603625114D-16,0.7407124321D-19/
      DATA CG2(1,1),CG2(1,2),CG2(1,3)/
     &     -.1403086182D004,0.1802594763D001,-.2097279215D-03/
      DATA CG2(2,1),CG2(2,2),CG2(2,3)/
     &     0.3817195017D000,-.5394444747D-03,0.1855203702D-06/
      DATA CG2(3,1),CG2(3,2),CG2(3,3)/
     &     -.6449501159D-04,0.8437637660D-07,-.2713755001D-10/
      DATA CG2(4,1),CG2(4,2),CG2(4,3)/
     &     0.7823817858D-08,-.1053834646D-10,0.3629590764D-14/
      DATA CSF1/0.6970887859D2,0.3337529994D2,0.2318240735D1,
     &         0.1840599513D0,-.5245502284D-2,0.2878007027D-2,
     &         0.1753652324D-2,-.4334859629D-3,0.3325699282D-4/
      DATA CSF2/0.8408618802D6,0.3637413208D6,-.4634506669D6,
     &         0.1130306339D6,-.4350217298D3,-.3898988188D4,
     &         0.6697399434D3,-.4730726377D2,0.1265125057D1/
      DATA CSF3/0.9060030436D3,-.1426813520D2,0.1522233257D1,
     &         -.6973992961D0,0.1743091663D0,-.2319717696D-1,
     &         0.1694019149D-2,-.6454771710D-4,0.1003003098D-5/
      DATA CSG1/0.1105836875D4,0.1436943768D2,0.8018288621D0,
     &         0.1617232913D-1, -0.1501147505D-2,0.0D0,0.0D0,0.0D0,
     &         0.0D0,-.1237675562D-4,0.3004773304D-5,-.2062390734D-6/
      DATA CSG2/-0.2234264997D7,0.1231247634D7,-0.1978847871D6,
     &          0.1859988044D2,-0.2765701318D1,0.1036033878D4,
     &         -0.2143423131D3,0.1690507762D2,-0.4864322134D0/
      DATA CSG3/0.9059978254D3,0.5561957539D1,0.3434189609D1,
     &        -0.6406390628D0,0.5918579484D-1,-0.2725378570D-2,
     &         0.5006336938D-4/
      END
