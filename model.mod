; H. Koshimichi, Analyzing the worldwide progression of COVID-19 cases and deaths using nonlinear mixed-effects model
; NONMEM code, 2024

$PROBLEM    Covid-19
$INPUT      ID WEEK=TIME PD DV1 LDV=DV MDV EXDMDV DC0
$INPUT      BLQFLG NLLOQ DEAVG DEUPPER DELOWER NWCDEATH LOCKNFLG LOCKCFLG PNVAC VAC_TRANS
$INPUT      V_ALPHA V_BETA V_GAMMA V_DELTA V_OMICRON V_OTHERS
$INPUT      LPOP LCPOP LGDP LATITUDE AGED65 CARDIODR 
$DATA       nmdata.csv  IGNORE=@
            IGNORE=(MDV.NE.0) IGNORE=(WEEK.LT.0) IGNORE=(LPOP.LT.6) IGNORE=(WEEK.GE.96)
            IGNORE=(CARDIODR.LE.0) IGNORE=(V_OMICRON.GE.25)
$SUBROUTINE ADVAN13 TOL=5
$ABBREVIATED NOCHECKMU
$MODEL      NCOMP=10
            COMP=(SUSC)      ; Susceptible
            COMP=(EXINC)     ; Exposed (incubation period, without infectiousness)
            COMP=(ASY)       ; Asymptomatic
            COMP=(ASYQ)      ; Asymptomatic Quarantined
            COMP=(PRESY)     ; Presymptomatic
            COMP=(PRESYQ)    ; Presymptomatic Quarantined
            COMP=(SY1)       ; Symptomatic to recovery
            COMP=(SY1Q)      ; Symptomatic to recovery Quarantined
            COMP=(SY2)      ; Symptomatic to death
            COMP=(SY2Q)     ; Symptomatic to death Quarantined

$PK
      MU_1 = THETA(4)
      MU_2 = THETA(6) + THETA(27)*(LCPOP-3.24)/3.24 + THETA(28)*(CARDIODR-263)/263
      MU_3 = THETA(8)
      MU_4 = THETA(9)
      MU_5 = THETA(11)
      MU_6 = THETA(13) + THETA(24)*(LGDP-3.98)/3.98 + THETA(25)*(AGED65-8.74)/8.74 + THETA(26)*(LATITUDE-19.5)/19.5

      IF(NEWIND.NE.2)THEN
        LF = 0
        LFD = 0
        LTIME = 0
        LTIMED = 0
      ENDIF

      RHO0 = 0.80            ; Proportion of symptomatic
      ZASY = 0.35           ; Infectiousness asymptomatic vs symptomatic
      ZPSY = 0.63           ; Infectiousness presymptomatic vs symptomatic
      TONSET = 4.8/7          ; Time from infection to symptom onset (incubation period)
      TSINF = 2.5/7           ; Time of infection start before symptom onset
      TINC = TONSET - TSINF ; Time from infection to start time of having infectiousness
      TFINF = 8/7             ; Time of infection finish after symptom onset

      RHO0_LOGIT = LOG(RHO0/(1-RHO0))
      RHO_LOGIT = RHO0_LOGIT - (1-PNVAC/100)*LOG(1.70)
      RHO = 1/(1+EXP(-RHO_LOGIT))

      MTCDT  = EXP(MU_1 + ETA(1)) ; relate to alpha1
      MTEXDT = MTCDT

      SUSC0   = 1.0E8
      EXINC0  = DC0*EXP(THETA(5))

      GAMMA1 = LOG(2)/TINC
      GAMMA2 = LOG(2)/TSINF
      DELTA  = LOG(2)/TFINF
      ALPHA1 = LOG(2)/MTCDT
      ALPHA2 = LOG(2)/MTEXDT

      DR0_LOGIT = MU_2 + ETA(2)
      DR0 = 1/RHO/(1+EXP(-DR0_LOGIT))

      PRESY0 = RHO*GAMMA1*EXINC0
      ASY0 = (1-RHO)/RHO*PRESY0

      ASYQ0 = 0
      PRESYQ0 = 0
      SY1Q0_LOGIT = THETA(7)
      SY2Q0_LOGIT = SY1Q0_LOGIT

      EPSI1_MAX = EXP(MU_3 + ETA(3))
      EPSI1_INTE = MU_4 + ETA(4)
      EPSI1_SLP = EXP(THETA(10))
      EPSI1 = EPSI1_MAX*EXP(-EXP(EPSI1_INTE - EPSI1_SLP*TIME))
      EPSI2 = EPSI1_MAX
      THET = 1.0E-6
      
      SY0 = DC0/(EPSI1*(1-DR0)*(1-1/(1+EXP(-SY1Q0_LOGIT))) + EPSI2*DR0*(1-1/(1+EXP(-SY2Q0_LOGIT))))/(1+THET*GAMMA1*EXINC0)
      SY10 = (1-DR0)*SY0
      SY20 = DR0*SY0
      SY1Q0 = SY10/(1+EXP(-SY1Q0_LOGIT))
      SY2Q0 = SY20/(1+EXP(-SY2Q0_LOGIT))

      VAC2D_OTHER = LOG(THETA(14))
      VAC2D_DELTA = LOG(THETA(15))
      VACD_ODDS = (VAC2D_OTHER*(1-V_DELTA/100) + VAC2D_DELTA*V_DELTA/100)*(1-PNVAC/100)

      V_OTHERS2 = V_OTHERS + V_OMICRON
      VAR_ALPHA = LOG(THETA(16))
      VAR_BETA = LOG(THETA(17))
      VAR_GAMMA = LOG(THETA(18))
      VAR_DELTA = LOG(THETA(19))
      LVAR = VAR_ALPHA*(V_ALPHA/100) + VAR_BETA*(V_BETA/100) + VAR_GAMMA*(V_GAMMA/100) + VAR_DELTA*(V_DELTA/100)
      VAR = EXP(LVAR)

      VARD_ALPHA = LOG(THETA(20))
      VARD_BETA = LOG(THETA(21))
      VARD_GAMMA = LOG(THETA(22))
      VARD_DELTA = LOG(THETA(23))
      VARD_ODDS = (VARD_ALPHA*V_ALPHA + VARD_BETA*V_BETA + VARD_GAMMA*V_GAMMA + VARD_DELTA*V_DELTA)/100

      DR_LOGIT = DR0_LOGIT + VACD_ODDS + VARD_ODDS
      DR_ALL = 1/RHO/(1+EXP(-DR_LOGIT))

      LOCK1_LOGIT = MU_5 + ETA(5)
      LOCKDOWN1 = 1/(1+EXP(-LOCK1_LOGIT))
      LOCKDOWN2 = 0
      TLAST = 78
      LOCKDOWN = (LOCKDOWN2 - LOCKDOWN1)*TIME/TLAST + LOCKDOWN1
      IF(TIME.GT.TLAST) LOCKDOWN = LOCKDOWN2
      LOLMD = EXP(THETA(12))

      VAC = VAC_TRANS

      OC1 = 0
      IF(WEEK.LT.4) OC1 = 1
      OC2 = 0
      IF(WEEK.GE.4.AND.WEEK.LT.8) OC2 = 1
      OC3 = 0
      IF(WEEK.GE.8.AND.WEEK.LT.12) OC3 = 1
      OC4 = 0
      IF(WEEK.GE.12.AND.WEEK.LT.16) OC4 = 1
      OC5 = 0
      IF(WEEK.GE.16.AND.WEEK.LT.20) OC5 = 1
      OC6 = 0
      IF(WEEK.GE.20.AND.WEEK.LT.24) OC6 = 1
      OC7 = 0
      IF(WEEK.GE.24.AND.WEEK.LT.28) OC7 = 1
      OC8 = 0
      IF(WEEK.GE.28.AND.WEEK.LT.32) OC8 = 1
      OC9 = 0
      IF(WEEK.GE.32.AND.WEEK.LT.36) OC9 = 1
      OC10 = 0
      IF(WEEK.GE.36.AND.WEEK.LT.40) OC10 = 1
      OC11 = 0
      IF(WEEK.GE.40.AND.WEEK.LT.44) OC11 = 1
      OC12 = 0
      IF(WEEK.GE.44.AND.WEEK.LT.48) OC12 = 1
      OC13 = 0
      IF(WEEK.GE.48.AND.WEEK.LT.52) OC13 = 1
      OC14 = 0
      IF(WEEK.GE.52.AND.WEEK.LT.56) OC14 = 1
      OC15 = 0
      IF(WEEK.GE.56.AND.WEEK.LT.60) OC15 = 1
      OC16 = 0
      IF(WEEK.GE.60.AND.WEEK.LT.64) OC16 = 1
      OC17 = 0
      IF(WEEK.GE.64.AND.WEEK.LT.68) OC17 = 1
      OC18 = 0
      IF(WEEK.GE.68.AND.WEEK.LT.72) OC18 = 1
      OC19 = 0
      IF(WEEK.GE.72.AND.WEEK.LT.76) OC19 = 1
      OC20 = 0
      IF(WEEK.GE.76.AND.WEEK.LT.80) OC20 = 1
      OC21 = 0
      IF(WEEK.GE.80.AND.WEEK.LT.84) OC21 = 1
      OC22 = 0
      IF(WEEK.GE.84.AND.WEEK.LT.88) OC22 = 1
      OC23 = 0
      IF(WEEK.GE.88.AND.WEEK.LT.92) OC23 = 1
      OC24 = 0
      IF(WEEK.GE.92.AND.WEEK.LT.96) OC24 = 1

      BETA0 =  EXP(MU_6 + ETA(6))
      KPBETA = OC1*ETA(7) + OC2*ETA(8) + OC3*ETA(9) + OC4*ETA(10) + OC5*ETA(11) + OC6*ETA(12)
      KPBETA = KPBETA + OC7*ETA(13) + OC8*ETA(14) + OC9*ETA(15) + OC10*ETA(16) + OC11*ETA(17) + OC12*ETA(18)
      KPBETA = KPBETA + OC13*ETA(19) + OC14*ETA(20) + OC15*ETA(21) + OC16*ETA(22) + OC17*ETA(23) + OC18*ETA(24)
      KPBETA = KPBETA + OC19*ETA(25) + OC20*ETA(26) + OC21*ETA(27) + OC22*ETA(28) + OC23*ETA(29) + OC24*ETA(30)
      BETA = BETA0*EXP(KPBETA)

      IF(WEEK.LT.4) BETA1 = BETA
      IF(WEEK.GE.4.AND.WEEK.LT.8) BETA2 = BETA
      IF(WEEK.GE.8.AND.WEEK.LT.12) BETA3 = BETA
      IF(WEEK.GE.12.AND.WEEK.LT.16) BETA4 = BETA
      IF(WEEK.GE.16.AND.WEEK.LT.20) BETA5 = BETA
      IF(WEEK.GE.20.AND.WEEK.LT.24) BETA6 = BETA
      IF(WEEK.GE.24.AND.WEEK.LT.28) BETA7 = BETA
      IF(WEEK.GE.28.AND.WEEK.LT.32) BETA8 = BETA
      IF(WEEK.GE.32.AND.WEEK.LT.36) BETA9 = BETA
      IF(WEEK.GE.36.AND.WEEK.LT.40) BETA10 = BETA
      IF(WEEK.GE.40.AND.WEEK.LT.44) BETA11 = BETA
      IF(WEEK.GE.44.AND.WEEK.LT.48) BETA12 = BETA
      IF(WEEK.GE.48.AND.WEEK.LT.52) BETA13 = BETA
      IF(WEEK.GE.52.AND.WEEK.LT.56) BETA14 = BETA
      IF(WEEK.GE.56.AND.WEEK.LT.60) BETA15 = BETA
      IF(WEEK.GE.60.AND.WEEK.LT.64) BETA16 = BETA
      IF(WEEK.GE.64.AND.WEEK.LT.68) BETA17 = BETA
      IF(WEEK.GE.68.AND.WEEK.LT.72) BETA18 = BETA
      IF(WEEK.GE.72.AND.WEEK.LT.76) BETA19 = BETA
      IF(WEEK.GE.76.AND.WEEK.LT.80) BETA20 = BETA
      IF(WEEK.GE.80.AND.WEEK.LT.84) BETA21 = BETA
      IF(WEEK.GE.84.AND.WEEK.LT.88) BETA22 = BETA
      IF(WEEK.GE.88.AND.WEEK.LT.92) BETA23 = BETA
      IF(WEEK.GE.92.AND.WEEK.LT.96) BETA24 = BETA

;---- parameter initiation ---
      A_0(1) = SUSC0
      A_0(2) = EXINC0
      A_0(3) = ASY0
      A_0(4) = ASYQ0
      A_0(5) = PRESY0
      A_0(6) = PRESYQ0
      A_0(7) = SY10
      A_0(8) = SY1Q0
      A_0(9) = SY20
      A_0(10) = SY2Q0

;----- Compartment initialization -----

$DES
      IF(LOCKNFLG.EQ.1.OR.LOCKCFLG.EQ.1)THEN
        LFD = 1
        LTIMED = TIME
      ENDIF
      LOCKDOWND = 1 - LOCKDOWN*LFD
      IF(LOCKNFLG.EQ.0.AND.LOCKCFLG.EQ.0) LOCKDOWND = 1 - LOCKDOWN*EXP(-LOLMD*(TIME-LTIMED))*LFD

      DADT(1) = -BETA*LOCKDOWND*VAR*VAC*A(1)*(ZASY*(A(3) - A(4)) + ZPSY*(A(5) - A(6)) + (A(7) - A(8)) + (A(9) - A(10)))
      DADT(2) = BETA*LOCKDOWND*VAR*VAC*A(1)*(ZASY*(A(3) - A(4)) + ZPSY*(A(5) - A(6)) + (A(7) - A(8)) + (A(9) - A(10))) - GAMMA1*A(2)
      DADT(3) = (1-RHO)*GAMMA1*A(2) - DELTA*A(3)
      DADT(4) = - DELTA*A(4) + THET*(EPSI1*(A(7) - A(8)) + EPSI2*(A(9) - A(10)))*(A(3) - A(4))
      DADT(5) = RHO*GAMMA1*A(2) - GAMMA2*A(5)
      DADT(6) = - GAMMA2*A(6) + THET*(EPSI1*(A(7) - A(8)) + EPSI2*(A(9) - A(10)))*(A(5) - A(6))
      DADT(7) = (1-DR_ALL)*GAMMA2*A(5) - DELTA*A(7)
      DADT(8) = (1-DR_ALL)*GAMMA2*A(6) - DELTA*A(8) + EPSI1*(A(7) - A(8))
      DADT(9) = DR_ALL*GAMMA2*A(5) - ALPHA1*A(10) - ALPHA2*(A(9) - A(10))
      DADT(10) = DR_ALL*GAMMA2*A(6) - ALPHA1*A(10) + EPSI2*(A(9) - A(10))

$ERROR
      WCS1 = EPSI1*(A(7) - A(8))
      WCS2 = EPSI2*(A(9) - A(10))
      WCAS = THET*(WCS1 + WCS2)*(A(3) - A(4))
      WCPS = THET*(WCS1 + WCS2)*(A(5) - A(6))
      WKI = WCAS + WCPS + WCS1 + WCS2

      IF(LOCKNFLG.EQ.1.OR.LOCKCFLG.EQ.1)THEN
        LF = 1
        LTIME = TIME
      ENDIF
      LOCKDOWNT = 1 - LOCKDOWN*LF
      IF(LOCKNFLG.EQ.0.AND.LOCKCFLG.EQ.0) LOCKDOWNT = 1 - LOCKDOWN*EXP(-LOLMD*(TIME-LTIME))*LF

      WKCD = ALPHA1*A(10)
      WKEDT = ALPHA2*(A(9) - A(10))
      WKDT = DEAVG + WKCD + WKEDT

      A1 = A(1)
      A2 = A(2)
      A3 = A(3)
      A4 = A(4)
      A5 = A(5)
      A6 = A(6)
      A7 = A(7)
      A8 = A(8)
      A9 = A(9)
      A10 = A(10)

      W1 = THETA(1)
      W2 = THETA(2)
      W3 = THETA(3)

      ;-- M3 Method ---
      DEL = 1.0E-10
      ; PD = 1: weekly confirmed cases
        LLOQ = LOG10(ABS(NLLOQ) + DEL)
        IPRED = LOG10(ABS(WKI) + DEL)
        WTMP = W1
        EPSTMP = EPS(1)
      ; PD = 2: weekly confirmed deaths
      IF(PD.EQ.2) THEN
        LLOQ = LOG10(ABS(NLLOQ) + DEL)
        IPRED = LOG10(ABS(WKCD) + DEL)
        WTMP = W2
        EPSTMP = EPS(2)
      ENDIF
      ; PD = 3: weekly total deaths, including weekly confirmed deaths, weekly excess deaths, and estimated weekly covid-19 unrelated deaths
      IF(PD.EQ.3) THEN
        LLOQ = LOG10(ABS(DEUPPER) + DEL)
        IPRED = LOG10(ABS(WKDT) + DEL)
        WTMP = W3
        EPSTMP = EPS(3)
      ENDIF
      IF (COMACT==1) PREDV = IPRED
      DUM = (LLOQ - IPRED)/WTMP
      ;-- PHI(X) is the integral from .INF to X of N(0,1) .
      DUM2 = PHI(DUM) + DEL
      TYPE = 1
      IF (BLQFLG.EQ.1) TYPE = 2
      IF (MDV==1) TYPE = 0
      IF (TYPE.EQ.2) DV_LOQ = LLOQ

      Y = 0
      IRES = 0
      IWRES = 0

      ;-- Prediction DV>=LLOQ ---
      IF(TYPE.NE.2.OR.NPDE_MODE==1) THEN
        F_FLAG = 0
        Y = IPRED + WTMP*EPSTMP
        IRES = DV - IPRED
        IWRES = (DV - IPRED)/WTMP
      ENDIF

      ;-- Likelihood DV<LLOQ ---
      IF(TYPE.EQ.2.AND.NPDE_MODE==0) THEN
        F_FLAG = 1
        Y = DUM2
        IWRES = 0
        MDVRES = 1
      ENDIF

$THETA  (0.02, 0.2, 2) (0.03, 0.26, 3) (0.003, 0.026, 0.3) ; Proportional error PD1,2,3
$THETA  (-10, 0.15, 10) ; MTCDT
$THETA  (0.01, 3, 15) ; EXINC0
$THETA  (-10, -4.5, 10) ; DR_LOGIT
$THETA  (-10, -3.3, 10) ; SY1Q0_LOGIT
$THETA  (-15, -2.2, 10) (-10, -0.2, 10) (-10, -2.5, 10) ; EPSI1_MAX,_INT,_SLP
$THETA  (-10, -0.5, 10) ; LOCK1_LOGIT
$THETA  (-10, -2, 10) ; LOLMD
$THETA  (-25, -18.65, -10) ; BETA0
$THETA  (0.01, 0.1, 10) (0.01, 0.35, 10) ; VAC2D_OTHER,DELTA
$THETA  (0.01, 1, 10) (0.01, 1, 10) (0.01, 1, 10) (0.01, 1, 10); VAR_ALPHA,BETA,GAMMA,DELTA
$THETA  (0.01, 1, 10) (0.01, 1, 10) (0.01, 1, 10) (0.01, 1, 10); VARD_ALPHA,BETA,GAMMA,DELTA
$THETA  (-10, 0.1, 10) (-10, 0.1, 10) (-10, 0.1, 10); Cov on BETA0, LGDP,AGED65,LATITUDE
$THETA  (-10, 0.1, 10) (-10, 0.1, 10); Cov on DR0_LOGIT, LCPOP,CARDIODR

$OMEGA  BLOCK(6)
  0.01
  0 0.01
  0 0 0.01
  0 0 0 0.01
  0 0 0 0 0.01
  0 0 0 0 0 0.01
$OMEGA  BLOCK(1)
  0.01
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME
$OMEGA  BLOCK(1) SAME

$SIGMA BLOCK(3)
    1  FIXED
    0  1
    0  0  1


$PRIOR NWPRI
$THETAP
  0.2 FIXED 0.26 FIXED 0.026 FIXED
  0.15 FIXED
  3 FIXED
  -4.5 FIXED
  -3.3 FIXED
  -2.2 FIXED -0.2 FIXED -2.5 FIXED
  -0.5 FIXED
  -2 FIXED
  -18.65 FIXED
  0.1 FIXED 0.35 FIXED
  1 FIXED 1 FIXED 1 FIXED 1 FIXED
  1 FIXED 1 FIXED 1 FIXED 1 FIXED
  0.1 FIXED 0.1 FIXED 0.1 FIXED
  0.1 FIXED 0.1 FIXED

$THETAPV  BLOCK(28)
  10000 FIXED
  0 10000
  0 0 10000
  0 0 0 10000
  0 0 0 0 10000
  0 0 0 0 0 10000
  0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 0 0 10000 
  0 0 0 0 0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10000  
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10000
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10000

$OMEGAP  BLOCK(7)
  0.2 FIXED
  0 0.2
  0 0 0.2
  0 0 0 0.2
  0 0 0 0 0.2
  0 0 0 0 0 0.2
  0 0 0 0 0 0 0.2
$OMEGAPD (7 FIXED)

$SIGMAP BLOCK(3)
  0.06 FIXED
  0 0.06
  0 0 0.06

$SIGMAPD (3 FIXED)

$EST METHOD=BAYES INTERACTION LAPLACE NUMERICAL SLOW NBURN=100000 NITER=300000 THIN=10
FILE=output.ext SEED=12345 PRINT=100
CTYPE=0 CINTERVAL=100 CITER=10 CALPHA=0.05 NOABORT
$COV PRINT=E MATRIX=S UNCONDITIONAL
$TABLE NOAPPEND ONEHEADER NOPRINT FILE=output.csv
ID TIME Y DV1 DV MDV EVID PRED IPRED CWRES CWRESI PD
RHO0 RHO ZASY ZPSY TONSET TSINF TINC TFINF
SUSC0 EXINC0 ASY0 ASYQ0 PRESY0 PRESYQ0 SY10 SY1Q0 SY20 SY2Q0 SY0
SY1Q0_LOGIT SY2Q0_LOGIT GAMMA1 GAMMA2 DELTA ALPHA1 ALPHA2
EPSI1 EPSI1_MAX EPSI1_INTE EPSI1_SLP EPSI2 THET MTCDT MTEXDT
DR0_LOGIT DR0 DR_LOGIT DR_ALL 
BETA BETA0 BETA1 BETA2 BETA3 BETA4 BETA5 BETA6 BETA7 BETA8 BETA9 BETA10
BETA11 BETA12 BETA13 BETA14 BETA15 BETA16 BETA17 BETA18 BETA19 BETA20
BETA21 BETA22 BETA23 BETA24
LOCK1_LOGIT LOCKDOWN1 LOCKDOWN2 LOLMD LOCKDOWN LOCKDOWNT LOCKNFLG LOCKCFLG LF LTIME
VAC2D_OTHER VAC2D_DELTA VACD_ODDS VAC VAR_ALPHA VAR_BETA VAR_GAMMA VAR_DELTA VAR
VARD_ALPHA VARD_BETA VARD_GAMMA VARD_DELTA VARD_ODDS
WKI WKCD NWCDEATH WKEDT WKDT DEAVG DEUPPER DELOWER NLLOQ DC0
ETAS(1:LAST)
