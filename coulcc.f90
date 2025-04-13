MODULE COULCC_M
!-----------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only: spi=>int32, dpf=>real64, qpf=>real128 &
                                        &, stdin=>input_unit,stdout=>output_unit
!-----------------------------------------------------------------------
  CONTAINS
!-----------------------------------------------------------------------
!  COMPLEX COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD            
!-----------------------------------------------------------------------
!  A CONTINUED-FRACTION ALGORITHM FOR COULOMB                           
!  FUNCTIONS OF COMPLEX ORDER WITH COMPLEX ARGUMENTS                    
!-----------------------------------------------------------------------
!  REF. IN COMP. PHYS. COMMUN. 36 (1985) 363                            
!-----------------------------------------------------------------------
!  A. R. Barnett           Manchester  March   1981                     
!  modified I.J. Thompson  Daresbury, Sept. 1983 for Complex Functions  
!  modified A.R. Flores    WashU, April 2025 for Modern Compilers       
!-----------------------------------------------------------------------
!  original program  RCWFN       in    CPC  8 (1974) 377-395            
!                 +  RCWFF       in    CPC 11 (1976) 141-142            
!                 +  COULFG      in    CPC 27 (1982) 147-166            
!  description of real algorithm in    CPC 21 (1981) 297-314            
!  description of complex algorithm    JCP XX (1985) YYY-ZZZ            
!  this version written up       in    CPC XX (1985) YYY-ZZZ            
!-----------------------------------------------------------------------
!  COULCC returns F,G,G',G',SIG for complex XX, ETA1, and ZLMIN,        
!   for NL integer-spaced lambda values ZLMIN to ZLMIN+NL-1 inclusive,  
!   thus giving  complex-energy solutions to the Coulomb Schrodinger    
!   equation,to the Klein-Gordon equation and to suitable forms of      
!   the Dirac equation ,also spherical & cylindrical Bessel equations   
!                                                                       
!  if /MODE1/= 1  get F,G,F',G'   for integer-spaced lambda values      
!            = 2      F,G      unused arrays must be dimensioned in     
!            = 3      F,  F'          call to at least length (1)       
!            = 4      F                                                 
!            = 11 get F,H+,F',H+' ) if KFN=0, H+ = G + i.F        )     
!            = 12     F,H+        )       >0, H+ = J + i.Y = H(1) ) in  
!            = 21 get F,H-,F',H-' ) if KFN=0, H- = G - i.F        ) GC  
!            = 22     F,H-        )       >0, H- = J - i.Y = H(2) )     
!                                                                       
!     if MODE1<0 then the values returned are scaled by an exponential  
!                factor (dependent only on XX) to bring nearer unity    
!                the functions for large /XX/, small ETA & /ZL/ < /XX/  
!        Define SCALE = (  0        if MODE1 > 0                        
!                       (  IMAG(XX) if MODE1 < 0  &  KFN < 3            
!                       (  REAL(XX) if MODE1 < 0  &  KFN = 3            
!        then FC = EXP(-ABS(SCALE)) * ( F, j, J, or I)                  
!         and GC = EXP(-ABS(SCALE)) * ( G, y, or Y )                    
!               or EXP(SCALE)       * ( H+, H(1), or K)                 
!               or EXP(-SCALE)      * ( H- or H(2) )                    
!                                                                       
!  if  KFN  =  0,-1  complex Coulomb functions are returned   F & G     
!           =  1   spherical Bessel      "      "     "       j & y     
!           =  2 cylindrical Bessel      "      "     "       J & Y     
!           =  3 modified cyl. Bessel    "      "     "       I & K     
!                                                                       
!          and where Coulomb phase shifts put in SIG if KFN=0 (not -1)  
!                                                                       
!  The use of MODE and KFN is independent                               
!    (except that for KFN=3,  H(1) & H(2) are not given)                
!                                                                       
!  With negative orders lambda, COULCC can still be used but with       
!    reduced accuracy as CF1 becomes unstable. The user is thus         
!    strongly advised to use reflection formulae based on               
!    H+-(ZL,,) = H+-(-ZL-1,,) * exp +-i(sig(ZL)-sig(-ZL-1)-(ZL+1/2)pi)  
!                                                                       
!  Precision:  results to within 2-3 decimals of 'machine accuracy',    
!               but if CF1A fails because X too small or ETA too large  
!               the F solution  is less accurate if it decreases with   
!               decreasing lambda (e.g. for lambda.LE.-1 & ETA.NE.0)    
!              RERR in COMMON/STEED/ traces the main roundoff errors.   
!                                                                       
!   COULCC is coded for real*8 on IBM or equivalent  ACCUR >= 10**-14   
!          with a section of doubled REAL*16 for less roundoff errors.  
!          (If no doubled precision available, increase JMAX to eg 100) 
!   Use IMPLICIT COMPLEX*32 & REAL*16 on VS compiler ACCUR >= 10**-32   
!   For single precision CDC (48 bits) reassign REAL*8=REAL etc.        
!                                                                       
!   IFAIL  on input   = 0 : no printing of error messages               
!                    ne 0 : print error messages on file 6              
!   IFAIL  in output = -2 : argument out of range                       
!                    = -1 : one of the continued fractions failed,      
!                           or arithmetic check before final recursion  
!                    =  0 : All Calculations satisfactory               
!                    ge 0 : results available for orders up to & at     
!                             position NL-IFAIL in the output arrays.   
!                    = -3 : values at ZLMIN not found as over/underflow 
!                    = -4 : roundoff errors make results meaningless    
!-----------------------------------------------------------------------
!     Machine dependent constants :                                     
!                                                                       
!     ACCUR    target bound on relative error (except near 0 crossings) 
!               (ACCUR should be at least 100 * ACC8)                   
!     ACC8     smallest number with 1+ACC8 .ne.1 in REAL*8  arithmetic  
!     ACC16    smallest number with 1+ACC16.ne.1 in REAL*16 arithmetic  
!     FPMAX    magnitude of largest floating point number * ACC8        
!     FPMIN    magnitude of smallest floating point number / ACC8       
!     FPLMAX   LOG(FPMAX)                                               
!     FPLMIN   LOG(FPMIN)                                               
!                                                                       
!     ROUTINES CALLED :       LOGAM/CLOGAM/CDIGAM,                      
!                             F20, CF1A, RCF, CF1C, CF2, F11, CF1R      
!     Intrinsic functions :   MIN, MAX, SQRT, REAL, IMAG, ABS, LOG, EXP
!      (Generic names)        NINT, MOD, ATAN, ATAN2, COS, SIN, DCMPLX, 
!                             SIGN, CONJG, INT, TANH                    
!     Note: Statement fntn.   NINTC = integer nearest to a complex no.  
!                                                                       
!     Parameters determining region of calculations :                   
!                                                                       
!        R20      estimate of (2F0 iterations)/(CF2 iterations)         
!        ASYM     minimum X/(ETA**2+L) for CF1A to converge easily      
!        XNEAR    minimum ABS(X) for CF2 to converge accurately         
!        LIMIT    maximum no. iterations for CF1, CF2, and 1F1 series   
!        JMAX     size of work arrays for Pade accelerations            
!        NDROP    number of successive decrements to define instability 
!-----------------------------------------------------------------------
    SUBROUTINE COULCC(XX,ETA1,ZLMIN,NL,FC,GC,FCP,GCP,SIG,MODE1,KFN,IFAIL)
!-----------------------------------------------------------------------
      IMPLICIT COMPLEX(dpf) (A-H,O-Z)                                   
      PARAMETER(JMAX=50)                                                
      DIMENSION FC(NL),GC(NL),FCP(NL),GCP(NL),SIG(NL),XRCF(JMAX,4)      
      LOGICAL PR,ETANE0,IFCP,RLEL,DONEM,UNSTAB,ZLNEG,AXIAL,NOCF2,NPINT  
      REAL(dpf) ERR,RERR,ABSC,ACCUR,ACCT,ACC8,ACCH,ACC16,ACCB, XNEAR    
      REAL(dpf) ZERO,ONE,TWO,HALF,HPI,TLOG,FPMAX,FPMIN,FPLMIN,FPLMAX    
      REAL(dpf) PACCQ,EPS,OFF,SCALE,SF,SFSH,TA,RK,OMEGA,R20,ASYM,ABSX   
!-----------------------------------------------------------------------
      COMMON /STEED/ RERR,NFP,N11,NPQ(2),N20,KAS(2)                     
!***  common blocks are for information & storage only.                 
!     (they are not essential to working of the code)                   
      COMMON /RCFCM1/ PK,EK,CLGAA,CLGAB,CLGBB,DSIG,TPK1,W,RL,FCL1,Q &   
      & ,GAM,HCL,HPL,FCM,HCL1,ALPHA,BETA,PL                             
      EQUIVALENCE            (PK,XRCF(1,1))                             
!-----------------------------------------------------------------------
      DATA ZERO,ONE,TWO,LIMIT /0.0D+0, 1.0D+0, 2.0D+0, 20000 /          
      DATA HALF, CI / 0.5D+0, (0D+0, 1D+0) /                            
      DATA FPMAX,FPMIN,FPLMAX,FPLMIN / 1D+60,1D-60 ,140D+0, -140D+0 /   
      DATA R20,ASYM,XNEAR,NDROP / 3., 3., .5, 5 /                       
      DATA ACCUR, ACC8, ACC16 / 1D-14, 2D-16, 3D-33 /                   
!-----------------------------------------------------------------------
      NINTC(W) = NINT(REAL(REAL(W)))                                    
      ABSC(W) = ABS(REAL(W)) + ABS(IMAG(W))                             
      NPINT(W,ACCB) = ABSC(NINTC(W)-W).LT.ACCB .AND. REAL(W).LT.HALF    
!-----------------------------------------------------------------------
      MODE = MOD(ABS(MODE1),10)                                         
      IFCP = MOD(MODE,2).EQ.1                                           
      PR = IFAIL.NE.0                                                   
      IFAIL = -2                                                        
      N11   = 0                                                         
      NFP   = 0                                                         
      KAS(1)= 0                                                      
      KAS(2)= 0                                                      
      NPQ(1)= 0                                                      
      NPQ(2)= 0                                                      
      N20 = 0                                                           
      HPI = TWO*ATAN(ONE)                                               
      TLOG = LOG(TWO)                                                   
      ACCUR = MAX(ACCUR, 50*ACC8)                                       
      ACCT = ACCUR * .5                                                 
!-----------------------------------------------------------------------
!     initialise the log-gamma function
!-----------------------------------------------------------------------
!     DUMMY = LOGAM(ACC8)                                              
      ACCH  = SQRT(ACCUR)                                               
      ACCB  = SQRT(ACCH)                                                
      RERR = ACCT                                                       
!-----------------------------------------------------------------------
      CIK = ONE                                                         
      IF (KFN.GE.3) CIK = CI * SIGN(ONE,ACC8-IMAG(XX))                
      X     = XX * CIK                                                  
      ETA   = ETA1                                                      
      IF (KFN .GT. 0) ETA = ZERO                                         
      ETANE0  = ABSC(ETA).GT.ACC8                                    
      ETAI = ETA*CI                                                  
      DELL  = ZERO                                                      
      IF (KFN .GE. 2)  DELL = HALF                                       
      ZM1   = ZLMIN - DELL                                              
      SCALE = ZERO                                                      
      IF (MODE1.LT.0) SCALE = IMAG(X)                                    
!-----------------------------------------------------------------------
      M1 = 1                                                            
      L1  = M1 + NL - 1                                                 
      RLEL = ABS(IMAG(ETA)) + ABS(IMAG(ZM1)) .LT. ACC8                  
      ABSX = ABS(X)                                                     
      AXIAL = RLEL .AND. ABS(IMAG(X)) .LT. ACC8 * ABSX                  
      IF (MODE.LE.2 .AND. ABSX.LT.FPMIN) GO TO 310                       
      XI  = ONE/X                                                       
      XLOG = LOG(X)                                                     
!-----------------------------------------------------------------------
! log with cut along the negative real axis see also OMEGA
!-----------------------------------------------------------------------
      ID = 1                                                            
      DONEM = .FALSE.                                                   
      UNSTAB = .FALSE.                                               
      LF = M1                                                           
      IFAIL = -1                                                        
  10  ZLM = ZM1 + LF - M1                                            
      ZLL = ZM1 + L1 - M1                                            
!-----------------------------------------------------------------------
! *** ZLL  is final lambda value, or 0.5 smaller for J,Y Bessels
!-----------------------------------------------------------------------
      Z11 = ZLL                                                 
      IF (ID.LT.0) Z11 = ZLM                                     
      P11 = CI*SIGN(ONE,ACC8-IMAG(ETA))                         
      LAST = L1                                                         
!-----------------------------------------------------------------------
! *** Find phase shifts and Gamow factor at lambda = ZLL
!-----------------------------------------------------------------------
      PK = ZLL + ONE                                                    
      AA = PK - ETAI                                                    
      AB = PK + ETAI                                                    
      BB = TWO*PK                                                       
      ZLNEG = NPINT(BB,ACCB)                                         
      CLGAA = CLOGAM(AA)                                 
      CLGAB = CLGAA                                      
      IF (ETANE0.AND..NOT.RLEL)  CLGAB = CLOGAM(AB)                   
      IF (ETANE0.AND.     RLEL)  CLGAB = CONJG(CLGAA)                 
      SIGMA = (CLGAA - CLGAB) * CI*HALF                              
      IF (KFN.EQ.0) SIG(L1) = SIGMA                                   
      IF(.NOT.ZLNEG) CLL = ZLL*TLOG- HPI*ETA - CLOGAM(BB) + (CLGAA+CLGAB)*HALF
      THETA  = X - ETA*(XLOG+TLOG) - ZLL*HPI + SIGMA            
!-----------------------------------------------------------------------
      TA = (IMAG(AA)**2+IMAG(AB)**2+ABS(REAL(AA))+ABS(REAL(AB)))*HALF 
      IF (ID.GT.0 .AND. ABSX .LT. TA*ASYM .AND. .NOT.ZLNEG) GO TO 20     
!-----------------------------------------------------------------------
! ***  use CF1 instead of CF1A, if predicted to converge faster, 
!          (otherwise using CF1A as it treats negative lambda &  
!           recurrence-unstable cases properly)                  
!-----------------------------------------------------------------------
      RK = SIGN(ONE, REAL(X) + ACC8)                               
      P =  THETA                                                   
      IF (RK.LT.0) P = -X + ETA*(LOG(-X)+TLOG)-ZLL*HPI-SIGMA        
      F = RK * CF1A(X*RK,ETA*RK,ZLL,P,ACCT,JMAX,NFP,FEST,ERR,FPMAX,XRCF,XRCF(1,3), XRCF(1,4))
      FESL = LOG(FEST) + ABS(IMAG(X))                                   
      NFP = - NFP                                                    
      IF (NFP.LT.0   .OR.(UNSTAB.AND.ERR.LT.ACCB)) GO TO 40              
      IF(.NOT.ZLNEG .OR. UNSTAB.AND.ERR.GT.ACCB)  GO TO 20              
      IF (PR) WRITE(STDOUT,1060) '-L',ERR                                  
      IF (ERR.GT.ACCB) GO TO 280                                      
      GO TO 40                                                       
!-----------------------------------------------------------------------
! ***  evaluate CF1  =  f   =  F'(ZLL,ETA,X)/F(ZLL,ETA,X)             
!-----------------------------------------------------------------------
  20  IF (AXIAL) THEN                                                    
!  REAL VERSION   
        F = CF1R(REAL(X),REAL(ETA),REAL(ZLL),ACC8,SF,RK,ETANE0,LIMIT,ERR,NFP,ACCH,FPMIN,FPMAX,PR,'COULCC')
        FCL = SF                                                      
        TPK1= RK                                                      
      ELSE                                                           
!  COMPLEX VERSION
        F = CF1C(X,ETA,ZLL,ACC8,FCL,TPK1,ETANE0,LIMIT,ERR,NFP,ACCH,FPMIN,FPMAX,PR,'COULCC')
      END IF                                                          
!-----------------------------------------------------------------------
      IF (ERR.GT.ONE) GO TO 390                                          
!-----------------------------------------------------------------------
! ***  Make a simple check for CF1 being badly unstable:                
!-----------------------------------------------------------------------
      IF (ID.LT.0) GO TO 30                                              
      UNSTAB = REAL((ONE-ETA*XI)*CI*IMAG(THETA)/F).GT.ZERO     &        
        & .AND..NOT.AXIAL .AND. ABS(IMAG(THETA)).GT.-LOG(ACC8)*.5 &        
        & .AND. ABSC(ETA)+ABSC(ZLL).LT.ABSC(X)                             
      IF (UNSTAB) GO TO 60                                               
!-----------------------------------------------------------------------
! *** compare accumulated phase FCL with asymptotic phase for G(k+1) :  
!     to determine estimate of F(ZLL) (with correct sign) to start recur
!-----------------------------------------------------------------------
  30  W    =  X*X  *(HALF/TPK1 + ONE/TPK1**2) + ETA*(ETA-TWO*X)/TPK1     
      FESL = (ZLL+ONE) * XLOG + CLL - W - LOG(FCL)                    
  40  FESL = FESL - ABS(SCALE)                                          
      RK   =        MAX(REAL(FESL), FPLMIN*HALF)                    
      FESL = DCMPLX(MIN(RK,   FPLMAX*HALF ) , IMAG(FESL))           
      FEST = EXP(FESL)                                                   
!-----------------------------------------------------------------------
      RERR = MAX(RERR, ERR, ACC8 * ABS(REAL(THETA)) )              
!-----------------------------------------------------------------------
      FCL = FEST                                                        
      FPL = FCL*F                                                       
      IF (IFCP) FCP(L1) = FPL                                            
      FC (L1) = FCL                                                    
!-----------------------------------------------------------------------
! *** downward recurrence to lambda = ZLM. array GC,if present,stores RL
!-----------------------------------------------------------------------
      I  = MAX(-ID, 0)                                                  
      ZL  = ZLL + I                                                     
      MONO = 0                                                       
      OFF = ABS(FCL)                                                  
      TA = ABSC(SIGMA)                                               
      DO 70  L  = L1-ID,LF,-ID                                          
        IF (ETANE0) THEN                                                
          IF (RLEL) THEN                                            
            DSIG = ATAN2(REAL(ETA),REAL(ZL))                    
            RL = SQRT(REAL(ZL)**2 + REAL(ETA)**2)               
          ELSE                                                  
            AA = ZL - ETAI                                      
            BB = ZL + ETAI                                      
            IF (ABSC(AA).LT.ACCH.OR.ABSC(BB).LT.ACCH) GOTO 50    
            DSIG = (LOG(AA) - LOG(BB)) * CI*HALF                
            RL = AA * EXP(CI*DSIG)                              
          END IF                                                  
          IF (ABSC(SIGMA).LT.TA*HALF) THEN                            
!-----------------------------------------------------------------------
! re-calculate SIGMA because of accumulating roundoffs:                 
!-----------------------------------------------------------------------
            SL =(CLOGAM(ZL+I-ETAI)-CLOGAM(ZL+I+ETAI))*CI*HALF       
            RL = (ZL - ETAI) * EXP(CI*ID*(SIGMA - SL))              
            SIGMA = SL                                              
            TA = ZERO                                               
          ELSE                                                      
            SIGMA = SIGMA - DSIG * ID                               
          END IF                                                     
          TA = MAX(TA, ABSC(SIGMA))                               
          SL =  ETA  + ZL*ZL*XI                                   
          PL = ZERO                                               
          IF (ABSC(ZL).GT.ACCH) PL = (SL*SL - RL*RL)/ZL            
          FCL1  = (FCL *SL + ID*ZL*FPL)/RL                           
          SF = ABS(FCL1)                                            
          IF (SF.GT.FPMAX) GO TO 350                        
          FPL   = (FPL *SL + ID*PL*FCL)/RL                           
          IF (MODE .LE. 1) GCP(L+ID)= PL * ID                         
        ELSE                                                            
!-----------------------------------------------------------------------
!  ETA = 0, including Bessels.  NB RL==SL                               
!-----------------------------------------------------------------------
          RL = ZL* XI                                                  
          FCL1 = FCL * RL + FPL*ID                                     
          SF = ABS(FCL1)                                            
          IF (SF.GT.FPMAX) GO TO 350                         
          FPL  =(FCL1* RL - FCL) * ID                                  
        END IF                                                           
        IF (SF.LT.OFF) THEN                                        
          MONO = MONO + 1                                        
        ELSE                                                    
          MONO = 0                                               
        END IF                                                   
        FCL =  FCL1                                                  
        OFF = SF                                                     
        FC(L) =  FCL                                                   
        IF (IFCP) FCP(L)  = FPL                                         
        IF (KFN.EQ.0) SIG(L) = SIGMA                                  
        IF (MODE .LE. 2) GC(L+ID) = RL                                
        ZL = ZL - ID                                                      
        IF (MONO.LT.NDROP) GO TO 70                                        
        IF (AXIAL .OR. REAL(ZLM)*ID.GT.-NDROP.AND..NOT.ETANE0) GO TO 70    
        UNSTAB = .TRUE.                                                
!-----------------------------------------------------------------------
! ***    take action if cannot or should not recur below this ZL:       
!-----------------------------------------------------------------------
  50    ZLM = ZL                                                       
        LF = L                                                         
        IF (ID.LT.0) GO TO 380                                       
        IF(.NOT.UNSTAB) LF = L + 1                                     
        IF (L+MONO.LT.L1-2 .OR. ID.LT.0 .OR. .NOT.UNSTAB) GO TO 80      
!-----------------------------------------------------------------------
! otherwise, all L values (for stability) should be done    
!            in the reverse direction:                      
!-----------------------------------------------------------------------
        GO TO 60                                                       
  70  CONTINUE                                                          
      GO TO 80                                                          
  60  ID = -1                                                     
      LF = L1                                                     
      L1 = M1                                                     
      RERR = ACCT                                                 
      GO TO 10                                                    
  80  IF (FCL .EQ. ZERO) FCL = + ACC8                                    
      F = FPL/FCL                                                      
!-----------------------------------------------------------------------
! *** Check, if second time around, that the 'f' values agree]          
!-----------------------------------------------------------------------
      IF (ID.GT.0) FIRST = F                                             
      IF (DONEM) RERR = MAX(RERR, ABSC(F-FIRST)/ABSC(F))                 
      IF (DONEM) GO TO 90                                                
!-----------------------------------------------------------------------
      NOCF2 = .FALSE.                                                  
      THETAM  = X - ETA*(XLOG+TLOG) - ZLM*HPI + SIGMA                   
!-----------------------------------------------------------------------
! *** on left x-plane, determine OMEGA by requiring cut on -x axis      
!     on right x-plane, choose OMEGA (using estimate based on THETAM)   
!       so H(omega) is smaller and recurs upwards accurately.           
!     (x-plane boundary is shifted to give CF2(LH) a chance to converge)
!-----------------------------------------------------------------------
      OMEGA = SIGN(ONE,IMAG(X)+ACC8)               
      IF (REAL(X).GE.XNEAR) OMEGA = SIGN(ONE,IMAG(THETAM)+ACC8)          
!-----------------------------------------------------------------------
      SFSH = EXP(OMEGA*SCALE - ABS(SCALE))                           
      OFF = EXP(MIN(TWO*MAX(ABS(IMAG(X)),ABS(IMAG(THETAM)),ABS(IMAG(ZLM))*3),FPLMAX))
      EPS = MAX(ACC8 , ACCT * HALF / OFF)                           
!-----------------------------------------------------------------------
! ***    Try first estimated omega, then its opposite,                  
!        to find the H(omega) linearly independent of F                 
!        i.e. maximise  CF1-CF2 = 1/(F H(omega)) , to minimise H(omega) 
!-----------------------------------------------------------------------
  90  DO 100 L=1,2                                                      
        LH = 1                                                          
        IF (OMEGA.LT.ZERO) LH = 2                                        
        PM = CI*OMEGA                                                   
        ETAP = ETA * PM                                                 
        IF (DONEM) GO TO 130                                             
        PQ1 = ZERO                                                      
        PACCQ = ONE                                                     
        KASE = 0                                                        
!-----------------------------------------------------------------------
! ***   Check for small X, i.e. whether to avoid CF2 :                  
!-----------------------------------------------------------------------
        IF (MODE.GE.3 .AND. ABSX.LT.ONE ) GO TO 190                      
        IF (MODE.LT.3 .AND. (NOCF2 .OR. ABSX.LT.XNEAR .AND.  &           
          &  ABSC(ETA)*ABSX .LT. 5 .AND. ABSC(ZLM).LT.4)) THEN          
          KASE = 5                                                      
          GO TO 120                                                     
        END IF                                                          
!-----------------------------------------------------------------------
! ***  Evaluate   CF2 : PQ1 = p + i.omega.q  at lambda = ZLM            
!-----------------------------------------------------------------------
        PQ1 = CF2(X,ETA,ZLM,PM,EPS,LIMIT,ERR,NPQ(LH),ACC8,ACCH,PR,ACCUR,DELL,'COULCC')
!-----------------------------------------------------------------------
        ERR = ERR * MAX(ONE,ABSC(PQ1)/MAX(ABSC(F-PQ1),ACC8) )            
        IF (ERR.LT.ACCH) GO TO 110                                  
!-----------------------------------------------------------------------
! *** check if impossible to get F-PQ accurately because of cancellation
!-----------------------------------------------------------------------
        NOCF2 = REAL(X).LT.XNEAR .AND. ABS(IMAG(X)).LT.-LOG(ACC8)
!-----------------------------------------------------------------------
!     original guess for OMEGA (based on THETAM) was wrong              
!     Use KASE 5 or 6 if necessary if Re(X) < XNEAR                     
!-----------------------------------------------------------------------
  100 OMEGA = - OMEGA                                        
!-----------------------------------------------------------------------
      IF (UNSTAB) GO TO 360                                    
      IF (REAL(X).LT.-XNEAR .AND. PR) WRITE(STDOUT,1060) '-X',ERR   
  110 RERR = MAX(RERR,ERR)                                          
!-----------------------------------------------------------------------
! ***  establish case of calculation required for irregular solution    
!-----------------------------------------------------------------------
  120 IF (KASE.GE.5) GO TO 130                                         
      IF (REAL(X) .GT. XNEAR) THEN                                     
!-----------------------------------------------------------------------
!  estimate errors if KASE 2 or 3 were to be used:                      
!-----------------------------------------------------------------------
        PACCQ = EPS * OFF * ABSC(PQ1) / MAX(ABS(IMAG(PQ1)),ACC8)      
      END IF                                                          
      IF (PACCQ .LT. ACCUR) THEN                                       
        KASE = 2                                                      
        IF (AXIAL) KASE = 3                                            
      ELSE                                                            
        KASE = 1                                                      
        IF (NPQ(1) * R20 .LT. JMAX) KASE = 4                       
!-----------------------------------------------------------------------
!  i.e. change to kase=4 if the 2F0 predicted to converge               
!-----------------------------------------------------------------------
      END IF                                                          
  130 GO TO (190,140,150,170,190,190), ABS(KASE)                      
  140 IF(.NOT.DONEM) PQ2 = CF2(X,ETA,ZLM,-PM,EPS,LIMIT,ERR,NPQ(3-LH),ACC8,ACCH,PR,ACCUR,DELL,'COULCC')
!-----------------------------------------------------------------------
! ***  Evaluate   CF2 : PQ2 = p - i.omega.q  at lambda = ZLM   (Kase 2) 
!-----------------------------------------------------------------------
      P = (PQ2 + PQ1) * HALF                                      
      Q = (PQ2 - PQ1) * HALF*PM                                   
!-----------------------------------------------------------------------
      GO TO 160                                                       
!-----------------------------------------------------------------------
  150 P = REAL(PQ1)                                               
      Q = IMAG(PQ1)                                               
!-----------------------------------------------------------------------
! ***   With Kase = 3 on the real axes, P and Q are real & PQ2 = PQ1*   
!-----------------------------------------------------------------------
      PQ2 = CONJG(PQ1)                                                
!-----------------------------------------------------------------------
! *** solve for FCM = F at lambda = ZLM,then find norm factor W=FCM/FCL 
!-----------------------------------------------------------------------
  160 W   = (PQ1 - F) * (PQ2 - F)                                       
      SF  = EXP(-ABS(SCALE))                                          
      FCM = SQRT(Q / W) * SF                                            
!-----------------------------------------------------------------------
!  any SQRT given here is corrected by                                  
!  using sign for FCM nearest to phase of FCL                           
!-----------------------------------------------------------------------
      IF (REAL(FCM/FCL).LT.ZERO) FCM  = - FCM                            
      GAM = (F - P)/Q                                                   
      TA = ABSC(GAM + PM)                                            
      PACCQ= EPS * MAX(TA,ONE/TA)                                    
      HCL = FCM * (GAM + PM) * (SFSH/(SF*SF))                           
!-----------------------------------------------------------------------
      IF (PACCQ.GT.ACCUR .AND. KASE.GT.0) THEN                           
!-----------------------------------------------------------------------
!  Consider a KASE = 1 Calculation                                      
!-----------------------------------------------------------------------
        F11V= F11(X,ETA,Z11,P11,ACCT,LIMIT,0,ERR,N11,FPMAX,ACC8,ACC16)
        IF (ERR.LT.PACCQ) GO TO 200                                    
      END IF                                                         
      RERR=MAX(RERR,PACCQ)                                              
!-----------------------------------------------------------------------
      GO TO 230                                                         
!-----------------------------------------------------------------------
! *** Arrive here if KASE = 4                                           
!     to evaluate the exponentially decreasing H(LH) directly.          
!-----------------------------------------------------------------------
  170 IF (DONEM) GO TO 180                                              
      AA = ETAP - ZLM                                                   
      BB = ETAP + ZLM + ONE                                             
      F20V = F20(AA,BB,-HALF*PM*XI, ACCT,JMAX,ERR,FPMAX,N20,XRCF)       
      IF (N20.LE.0) GO TO 190                                          
      RERR = MAX(RERR,ERR)                                            
      HCL = FPMIN                                                    
      IF (ABS(REAL(PM*THETAM)+OMEGA*SCALE).GT.FPLMAX) GO TO 330       
  180 HCL = F20V * EXP(PM * THETAM + OMEGA*SCALE)                     
      FCM = SFSH / ((F - PQ1) * HCL )                                 
!-----------------------------------------------------------------------
      GO TO 230                                                       
!-----------------------------------------------------------------------
! *** Arrive here if KASE=1   (or if 2F0 tried mistakenly & failed)     
!-----------------------------------------------------------------------
!  for small values of X, calculate F(X,SL) directly from 1F1           
!      using REAL*16 arithmetic if possible.                            
!  where Z11 = ZLL if ID>0, or = ZLM if ID<0                            
!-----------------------------------------------------------------------
  190 F11V = F11(X,ETA,Z11,P11,ACCT,LIMIT,0,ERR,N11,FPMAX,ACC8,ACC16) 
!-----------------------------------------------------------------------
  200 IF (N11.LT.0) THEN                                           
!-----------------------------------------------------------------------
!    F11 failed from BB = negative integer                              
!-----------------------------------------------------------------------
        WRITE(STDOUT,1060) '-L',ONE                                   
        GO TO 390                                                
      END IF                                                    
!-----------------------------------------------------------------------
      IF (ERR.GT.PACCQ .AND. PACCQ.LT.ACCB) THEN                   
!-----------------------------------------------------------------------
!  Consider a KASE 2 or 3 calculation                                   
!-----------------------------------------------------------------------
        KASE = -2                                               
        IF (AXIAL) KASE = -3                                     
        GO TO 130                                               
      END IF                                                   
      RERR = MAX(RERR, ERR)                                          
      IF (ERR.GT.FPMAX) GO TO 370                                     
      IF (ID.LT.0) CLL = Z11*TLOG-HPI*ETA-CLOGAM(BB)+CLOGAM(Z11+ONE+P11*ETA)-P11*SIGMA
      EK = (Z11+ONE)*XLOG - P11*X + CLL  - ABS(SCALE)                 
      IF (ID.GT.0) EK = EK - FESL + LOG(FCL)                             
      IF (REAL(EK).GT.FPLMAX) GO TO 350                               
      IF (REAL(EK).LT.FPLMIN) GO TO 340                               
      FCM = F11V * EXP(EK)                                              
!-----------------------------------------------------------------------
      IF (KASE.GE.5) THEN                                                
        IF (ABSC(ZLM+ZLM-NINTC(ZLM+ZLM)).LT.ACCH) KASE = 6               
!-----------------------------------------------------------------------
! ***  For abs(X) < XNEAR, then CF2 may not converge accurately, so     
! ***      use an expansion for irregular soln from origin :            
!-----------------------------------------------------------------------
        SL = ZLM                                                       
        ZLNEG = REAL(ZLM) .LT. -ONE + ACCB                          
        IF (KASE.EQ.5 .OR. ZLNEG) SL = - ZLM - ONE                      
        PK = SL + ONE                                                  
        AA = PK - ETAP                                              
        AB = PK + ETAP                                              
        BB = TWO*PK                                                 
        CLGAA = CLOGAM(AA)                                 
        CLGAB = CLGAA                                      
        IF (ETANE0) CLGAB = CLOGAM(AB)                                 
        CLGBB = CLOGAM(BB)                                 
        IF (KASE.EQ.6 .AND. .NOT.ZLNEG) THEN                          
          IF (NPINT(AA,ACCUR)) CLGAA = CLGAB - TWO*PM*SIGMA          
          IF (NPINT(AB,ACCUR)) CLGAB = CLGAA + TWO*PM*SIGMA          
        END IF                                                      
        CLL = SL*TLOG- HPI*ETA - CLGBB + (CLGAA + CLGAB) * HALF       
        DSIG = (CLGAA - CLGAB) * PM*HALF                              
        IF (KASE.EQ.6) P11 = - PM                                   
        EK = PK * XLOG - P11*X + CLL  - ABS(SCALE)                   
        SF = EXP(-ABS(SCALE))                              
        CHI = ZERO                                         
        IF(.NOT.( KASE.EQ.5 .OR. ZLNEG ) ) GO TO 210                  
!-----------------------------------------------------------------------
! *** Use  G(l)  =  (cos(CHI) * F(l) - F(-l-1)) /  sin(CHI)             
!-----------------------------------------------------------------------
!      where CHI = sig(l) - sig(-l-1) - (2l+1)*pi/2                     
!-----------------------------------------------------------------------
        CHI = SIGMA - DSIG - (ZLM-SL) * HPI                            
        F11V=F11(X,ETA,SL,P11,ACCT,LIMIT,0,ERR,NPQ(1),FPMAX,ACC8,ACC16)
        RERR = MAX(RERR,ERR)                                
        IF (KASE.EQ.6) GO TO 210                                     
        FESL = F11V * EXP( EK )                                        
        FCL1 = EXP(PM*CHI) * FCM                                       
        HCL = FCL1 - FESL                                              
        RERR=MAX(RERR,ACCT*MAX(ABSC(FCL1),ABSC(FESL))/ABSC(HCL)) 
        HCL = HCL / SIN(CHI) * (SFSH/(SF*SF))                          
!-----------------------------------------------------------------------
        GO TO 220                                                        
!-----------------------------------------------------------------------
! *** Use the logarithmic expansion for the irregular solution (KASE 6) 
!        for the case that BB is integral so sin(CHI) would be zero.    
!-----------------------------------------------------------------------
  210   RL = BB - ONE                                                  
        N  = NINTC(RL)                                                 
        ZLOG = XLOG + TLOG - PM*HPI                                    
        CHI = CHI + PM * THETAM + OMEGA * SCALE + AB * ZLOG            
        AA  = ONE - AA                                              
        IF (NPINT(AA,ACCUR)) THEN                                       
          HCL = ZERO                                                  
        ELSE                                                           
          IF (ID.GT.0 .AND. .NOT.ZLNEG) F11V = FCM * EXP(-EK)       
          HCL = EXP(CHI-CLGBB-CLOGAM(AA))*(-1)**(N+1)*(F11V*ZLOG &   
            & +F11(X,ETA,SL,-PM,ACCT,LIMIT,2,ERR,NPQ(2),FPMAX,ACC8,ACC16))
          RERR = MAX(RERR,ERR)                                    
        END IF                                                       
        IF (N.GT.0) THEN                                                
           EK  = CHI + CLOGAM(RL) - CLGAB - RL*ZLOG                    
           DF  = F11(X,ETA,-SL-ONE,-PM,ZERO,N,0,ERR,L,FPMAX,ACC8,ACC16) 
           HCL = HCL + EXP(EK) * DF                                   
        END IF                                                       
!-----------------------------------------------------------------------
  220   PQ1 = F - SFSH/(FCM * HCL)                                    
!-----------------------------------------------------------------------
      ELSE !KASE<=4                                                   
        IF (MODE.LE.2) HCL = SFSH/((F - PQ1) * FCM)                   
        KASE = 1                                                     
      END IF                                                          
!-----------------------------------------------------------------------
! ***  Now have absolute normalisations for Coulomb Functions           
!          FCM & HCL  at lambda = ZLM                                   
!      so determine linear transformations for Functions required :     
!-----------------------------------------------------------------------
  230 IH = ABS(MODE1) / 10                                              
      IF (KFN.EQ.3) IH = (3-IMAG(CIK))/2  + HALF                       
      P11 = ONE                                                         
      IF (IH.EQ.1) P11 = CI                                              
      IF (IH.EQ.2) P11 = -CI                                             
      DF = - PM                                             
      IF (IH.GE.1) DF = - PM + P11                                       
      IF (ABSC(DF).LT.ACCH) DF = ZERO                                
!-----------------------------------------------------------------------
! *** Normalisations for spherical or cylindrical Bessel functions      
!-----------------------------------------------------------------------
      ALPHA = ZERO                                  
      IF (KFN .EQ. 1) ALPHA = XI                                    
      IF (KFN .GE. 2) ALPHA = XI*HALF                               
      BETA  = ONE                                   
      IF (KFN .EQ. 1) BETA  = XI                                    
      IF (KFN .GE. 2) BETA  = SQRT(XI/HPI)                          
      IF (KFN .GE. 2 .AND. REAL(BETA).LT.ZERO) BETA  = - BETA       
!-----------------------------------------------------------------------
      AA = ONE                                                          
      IF (KFN.GT.0) AA = -P11 * BETA                                     
      IF (KFN.GE.3) THEN                                                 
!-----------------------------------------------------------------------
! Calculate rescaling factors for I & K output                          
!-----------------------------------------------------------------------
        P = EXP((ZLM+DELL) * HPI * CIK)                                
        AA= BETA * HPI * P                                             
        BETA = BETA / P                                                
        Q = CIK * ID                                                   
      END IF                                                           
!-----------------------------------------------------------------------
! Calculate rescaling factors for GC output      
!-----------------------------------------------------------------------
      IF (IH.EQ.0) THEN                                                  
        TA = ABS(SCALE) + IMAG(PM)*SCALE                               
        RK = ZERO                                                      
        IF (TA.LT.FPLMAX) RK = EXP(-TA)                                 
      ELSE                                                             
        TA = ABS(SCALE) + IMAG(P11)*SCALE                              
!-----------------------------------------------------------------------
        IF (ABSC(DF).GT.ACCH .AND. TA.GT.FPLMAX) GO TO 320              
        IF (ABSC(DF).GT.ACCH) DF = DF * EXP(TA)                         
        SF = TWO * (LH-IH) * SCALE                                     
        RK = ZERO                                                      
        IF (SF.GT.FPLMAX) GO TO 320                                     
        IF (SF.GT.FPLMIN) RK = EXP(SF)                                  
      END IF                                                             
!-----------------------------------------------------------------------
      KAS((3-ID)/2) = KASE                                           
      W = FCM / FCL                                                     
      IF (LOG(ABSC(W))+LOG(ABSC(FC(LF))) .LT. FPLMIN) GO TO 340       
      IF (MODE.GE.3) GO TO 240                                        
      IF (ABSC(F-PQ1) .LT. ACCH*ABSC(F) .AND. PR) WRITE(STDOUT,1020) LH,ZLM+DELL
      HPL = HCL * PQ1                                                   
      IF (ABSC(HPL).LT.FPMIN.OR.ABSC(HCL).LT.FPMIN) GO TO 330         
!-----------------------------------------------------------------------
! *** IDward recurrence from HCL,HPL(LF) (stored GC(L) is RL if reqd)   
! *** renormalise FC,FCP at each lambda                                 
! ***    ZL   = ZLM - MIN(ID,0) here                                    
!-----------------------------------------------------------------------
  240 DO 270 L = LF,L1,ID                                               
!-----------------------------------------------------------------------
        FCL = W* FC(L)                                     
        IF (ABSC(FCL).LT.FPMIN) GO TO 340                  
        IF (IFCP) FPL = W*FCP(L)                                     
        FC(L)  = BETA * FCL                                
        IF (IFCP) FCP(L) = BETA * (FPL - ALPHA * FCL) * CIK          
        FC(L)  = TIDY(FC(L),ACCUR)                         
        IF (IFCP) FCP(L) = TIDY(FCP(L),ACCUR)                        
        IF (MODE .GE. 3) GO TO 260                                        
        IF (L.EQ.LF)  GO TO 250                                           
        ZL = ZL + ID                                      
        ZID= ZL * ID                                      
        RL = GC(L)                                        
        IF (ETANE0) THEN                                              
          SL = ETA + ZL*ZL*XI                               
          IF (MODE.EQ.1) THEN                                          
            PL = GCP(L)                                               
          ELSE                                                        
            PL = ZERO                                                 
            IF (ABSC(ZL).GT.ACCH) PL = (SL*SL - RL*RL)/ZID             
          END IF                                                       
          HCL1= (SL*HCL - ZID*HPL) / RL                           
          HPL = (SL*HPL - PL *HCL) / RL                           
        ELSE                                                           
          HCL1= RL * HCL - HPL * ID                                   
          HPL = (HCL - RL * HCL1) * ID                                
        END IF                                                          
        HCL = HCL1                                                
        IF (ABSC(HCL).GT.FPMAX) GO TO 320                               
  250   GC(L) = AA * (RK * HCL + DF * FCL)                             
        IF (MODE.EQ.1) GCP(L) = (AA*(RK*HPL+DF*FPL)-ALPHA*GC(L))*CIK
        GC(L) = TIDY(GC(L),ACCUR)                                      
        IF (MODE.EQ.1) GCP(L) = TIDY(GCP(L),ACCUR)                         
        IF (KFN.GE.3) AA = AA * Q                                       
  260   IF (KFN.GE.3) BETA = - BETA * Q                                 
!-----------------------------------------------------------------------
  270 LAST = MIN(LAST,(L1-L)*ID)                                     
!-----------------------------------------------------------------------
! *** Come here after all soft errors to determine how many L values ok 
!-----------------------------------------------------------------------
  280 IF (ID.GT.0 .OR.  LAST.EQ.0) IFAIL = LAST                         
      IF (ID.LT.0 .AND. LAST.NE.0) IFAIL = -3                           
!-----------------------------------------------------------------------
! *** Come here after ALL errors for this L range (ZLM,ZLL)             
!-----------------------------------------------------------------------
  290 IF (ID.GT.0 .AND. LF.NE.M1) GO TO 300                              
      IF (IFAIL.LT.0) RETURN                                          
      IF (RERR.GT.ACCB) WRITE(STDOUT,1070) RERR                            
      IF (RERR.GT.0.1) IFAIL = -4                                     
!-----------------------------------------------------------------------
      RETURN                                                         
!-----------------------------------------------------------------------
! *** so on first block, 'F' started decreasing monotonically,          
!                        or hit bound states for low ZL.                
!     thus redo M1 to LF-1 in reverse direction                         
!      i.e. do CF1A at ZLMIN & CF2 at ZLM (midway between ZLMIN & ZLMAX)
!-----------------------------------------------------------------------
  300 ID = -1                                                           
      IF(.NOT.UNSTAB) LF = LF - 1                                       
      DONEM = UNSTAB                                                    
      LF = MIN(LF,L1)                                                   
      L1 = M1                                                           
!-----------------------------------------------------------------------
      GO TO 10                                                          
!-----------------------------------------------------------------------
! ***    error messages                                                 
!-----------------------------------------------------------------------
 1000 FORMAT(/' COULCC: CANNOT CALCULATE IRREGULAR SOLUTIONS FOR X =',1P,2D10.2,', AS ABS(X) IS TOO SMALL'/)
 1010 FORMAT(' COULCC: AT ZL =',2F8.3,' ',A2,'REGULAR SOLUTION (',1P,2E10.1,') WILL BE ',A4,' THAN',E10.1)
 1020 FORMAT('0COULCC WARNING: LINEAR INDEPENDENCE BETWEEN ''F'' AND ''H(',I1,')'' &
             & IS LOST AT ZL =',2F7.2,' (EG. COULOMB EIGENSTATE, OR CF1 UNSTABLE)'/)
 1030 FORMAT(' COULCC: (ETA&L)/X TOO LARGE FOR CF1A, AND CF1 UNSTABLE ATL =',2F8.2)
 1040 FORMAT(' COULCC: OVERFLOW IN 1F1 SERIES AT ZL =',2F8.3,' AT TERM',I5)
 1050 FORMAT(' COULCC: BOTH BOUND-STATE POLES AND F-INSTABILITIES OCCUR'&
        &  ,', OR MULTIPLE INSTABILITIES PRESENT.'                         &
        & ,/,' TRY CALLING TWICE,  FIRST FOR ZL FROM',2F8.3,' TO',2F8.3,   &
        &  ' (INCL.)',/,20X,     'SECOND FOR ZL FROM',2F8.3,' TO',2F8.3)   
 1060 FORMAT('0COULCC WARNING: AS ''',A2,''' REFLECTION RULES NOT USED,ERRORS CAN BE UP TO',1P,D12.2/)
 1070 FORMAT('0COULCC WARNING: OVERALL ROUNDOFF ERROR APPROX.',1P,E11.1)
!-----------------------------------------------------------------------
  310 IF (PR) WRITE (6,1000) XX                                          
      RETURN                                                            
  320 IF (PR) WRITE(STDOUT,1010) ZL+DELL,'IR',HCL,'MORE',FPMAX                
      GO TO 280                                                         
  330 IF (PR) WRITE(STDOUT,1010) ZL+DELL,'IR',HCL,'LESS',FPMIN                
      GO TO 280                                                         
  340 IF (PR) WRITE(STDOUT,1010) ZL+DELL,'  ',FCL,'LESS',FPMIN                
      GO TO 280                                                         
  350 IF (PR) WRITE(STDOUT,1010) ZL+DELL,'  ',FCL,'MORE',FPMAX                
      GO TO 280                                                         
  360 IF (PR) WRITE(STDOUT,1030) ZLL+DELL                                     
      GO TO 280                                                         
  370 IF (PR) WRITE(STDOUT,1040) Z11,I                                        
      GO TO 390                                                         
  380 IF (PR) WRITE(STDOUT,1050) ZLMIN,ZLM,ZLM+ONE,ZLMIN+NL-ONE               
  390 IFAIL = -1                                                        
      GO TO 290                                                         
!-----------------------------------------------------------------------
    END SUBROUTINE COULCC                                             
!-----------------------------------------------------------------------
    COMPLEX(dpf) FUNCTION CF1C(X,ETA,ZL,EPS,FCL,TPK1,ETANE0,LIMIT,ERR,NFP &
                            &, ACCH,FPMIN,FPMAX,PR,CALLER)              
!-----------------------------------------------------------------------
      IMPLICIT COMPLEX(dpf)(A-H,O-Z)                                      
      LOGICAL PR,ETANE0                                                 
      REAL(dpf) ONE,TWO,EPS,ERR,ACCH,FPMIN,FPMAX,ABSC,SMALL,RK,PX          
      CHARACTER(len=6) CALLER                                                
      DATA ONE,TWO / 1D+0, 2D+0 /                                       
      ABSC(W) = ABS(REAL(W)) + ABS(IMAG(W))                             
!-----------------------------------------------------------------------
 1000 FORMAT(/' ',A6,': CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,ETA,X = ',/1X,1P,13D9.2/)
 1010 FORMAT(' ',A6,': CF1 HAS FAILED TO CONVERGE AFTER ',I10  ,' ITERATIONS AS ABS(X) =',F15.0)
!-----------------------------------------------------------------------
! *** Evaluate CF1  =  F   =  F'(ZL,ETA,X)/F(ZL,ETA,X)               
!     using complex arithmetic                                       
!-----------------------------------------------------------------------
      FCL = ONE                                                         
      XI  = ONE/X                                                        
      PK  = ZL + ONE                                                    
      PX  = PK  + LIMIT                                                 
  10  EK  = ETA / PK                                                    
      RK2 = ONE + EK*EK                                      
      F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI                           
      PK1 =  PK + ONE                                                   
      TPK1 = PK + PK1                                                
      TK  = TPK1*(XI + EK/PK1)                                          
      IF (ETANE0) THEN                                                   
!-----------------------------------------------------------------------
! ***   test ensures b1 .ne. zero for negative ETA etc.; fixup is exact.
!-----------------------------------------------------------------------
        IF (ABSC(TK) .GT. ACCH) GO TO 20                           
        FCL  = RK2/(ONE + (ETA/PK1)**2)                            
        SL   = TPK1*XI * (TPK1+TWO)*XI                             
        PK   =  TWO + PK                                           
        GO TO 10                                                   
      END IF                                                          
  20  D  =  ONE/TK                                                     
      DF = -FCL*RK2*D                                                  
      IF (REAL(PK).GT.REAL(ZL)+TWO) FCL = - RK2 * SL               
      FCL = FCL * D * TPK1 * XI                                   
      F   =  F  + DF                                                    
!-----------------------------------------------------------------------
! *** begin CF1 loop on PK = k = lambda + 1                           
!-----------------------------------------------------------------------
      RK    = ONE                                                       
      SMALL = SQRT(FPMIN)                                            
  30  PK    = PK1                                                       
      PK1 = PK1 + ONE                                                 
      TPK1 = PK + PK1                                                
      IF (ETANE0) THEN                                                
        EK  = ETA / PK                                               
        RK2 = ONE + EK*EK                                   
      END IF                                                         
      TK = TPK1*(XI + EK/PK1)                                        
      D  =  TK - D*RK2                                               
      IF (ABSC(D) .GT. ACCH) GO TO 40                
      IF (PR) WRITE (6,1000) CALLER,D,DF,ACCH,PK,EK,ETA,X        
      RK= RK +   ONE                                            
      IF( RK .GT. TWO) GO TO 50               
  40  D = ONE/D                                                     
      FCL = FCL * D * TPK1*XI                                     
      IF (ABSC(FCL).LT.SMALL) FCL = FCL / SMALL                    
      IF (ABSC(FCL).GT.FPMAX) FCL = FCL / FPMAX                    
      DF = DF*(D*TK - ONE)                                           
      F  = F  + DF                                                   
      IF( REAL(PK) .GT. PX ) GO TO 50                           
      IF (ABSC(DF) .GE. ABSC(F)*EPS) GO TO 30                
      NFP = PK - ZL - 1                                       
      ERR = EPS * SQRT(REAL(NFP))                           
      CF1C = F                                                          
!-----------------------------------------------------------------------
      RETURN                                                            
!-----------------------------------------------------------------------
  50  IF (PR) WRITE (6,1010) CALLER,LIMIT,ABS(X)                         
      ERR = TWO                                                         
!-----------------------------------------------------------------------
    END FUNCTION CF1C                                                              
!-----------------------------------------------------------------------
    COMPLEX(dpf) FUNCTION CF2(X,ETA,ZL,PM,EPS,LIMIT,ERR,NPQ,ACC8,ACCH  & 
                            &, PR,ACCUR,DELL,CALLER)                                
!-----------------------------------------------------------------------
      IMPLICIT COMPLEX(dpf)(A-H,O-Z)                                    
      LOGICAL PR                                                        
      REAL(dpf) EPS,ERR,ACC8,ACCH,ACCUR,TA,RK                           
      REAL(dpf) ABSC,ZERO,HALF,ONE,TWO                                  
      CHARACTER(len=6) CALLER                                           
      DATA ZERO,HALF,ONE,TWO / 0D+0, .5D+0, 1D+0, 2D+0 /                
      ABSC(W) = ABS(REAL(W)) + ABS(IMAG(W))                             
!-----------------------------------------------------------------------
!                                    (omega)        (omega)             
! *** Evaluate  CF2  = p + PM.q  =  H   (ETA,X)' / H   (ETA,X)          
!                                    ZL             ZL                  
!     where PM = omega.i                                                
!-----------------------------------------------------------------------
      TA = TWO*LIMIT                                                    
      E2MM1 = ETA*ETA + ZL*ZL + ZL                                      
      ETAP = ETA * PM                                                   
      XI = ONE/X                                                        
      WI = TWO*ETAP                                                     
      RK = ZERO                                                         
      PQ = (ONE - ETA*XI) * PM                                          
      AA = -E2MM1 + ETAP                                                
      BB = TWO*(X - ETA + PM)                                           
      RL = XI * PM                                                   
      IF (ABSC(BB).LT.ACCH) THEN                                         
        RL = RL * AA / (AA + RK + WI)                                  
        PQ = PQ + RL * (BB + TWO*PM)                                   
        AA = AA + TWO*(RK+ONE+WI)                                   
        BB = BB + (TWO+TWO)*PM                                      
        RK = RK + (TWO+TWO)                                         
      END IF                                                          
      DD = ONE/BB                                                       
      DL = AA*DD* RL                                                    
  10  PQ = PQ + DL                                                   
      RK = RK + TWO                                                  
      AA = AA + RK + WI                                              
      BB = BB + TWO*PM                                               
      DD = ONE/(AA*DD + BB)                                          
      DL = DL*(BB*DD - ONE)                                          
      ERR = ABSC(DL)/ABSC(PQ)                                     
      IF (ERR.GE.MAX(EPS,ACC8*RK*HALF) .AND. RK.LE.TA) GO TO 10       
!-----------------------------------------------------------------------
      NPQ = RK/TWO                                                 
      PQ  = PQ + DL                                                
!-----------------------------------------------------------------------
 1000 FORMAT(' ',A6,': CF2(',I2,') NOT CONVERGED FULLY IN ',I7,         &
     & ' ITERATIONS, SO ERROR IN IRREGULAR SOLUTION =',1P,D11.2,' AT ZL &
     & =', 0P,2F8.3)  
!-----------------------------------------------------------------------
      IF (PR.AND.NPQ.GE.LIMIT-1 .AND. ERR.GT.ACCUR) WRITE(STDOUT,1000) CALLER,INT(IMAG(PM)),NPQ,ERR,ZL+DELL
      CF2 = PQ                                                          
!-----------------------------------------------------------------------
    END FUNCTION CF2                                                  
!-----------------------------------------------------------------------
    COMPLEX(dpf) FUNCTION F11(X,ETA,ZL,P,EPS,LIMIT,KIND,ERR,NITS,FPMAX,ACC8,ACC16)
!-----------------------------------------------------------------------
      IMPLICIT REAL(dpf)(A-H,O-Z)                                       
      COMPLEX(dpf) X,ETA,ZL,P,AA,BB,Z,CI                       
      COMPLEX(dpf) DD,G,F,AI,BI,T                                       
      LOGICAL ZLLIN                                                     
      REAL(qpf) AR,BR,GR,GI,DR,DI,TR,TI,UR,UI,FI,FI1,DEN                
      DATA ZERO,ONE,TWO / 0D+0, 1D+0, 2D+0 /, CI / (0D+0, 1D+0) /       
      ABSC(AA) = ABS(REAL(AA)) + ABS(IMAG(AA))                          
      NINTC(AA) = NINT(REAL(REAL(AA)))                                  
!-----------------------------------------------------------------------
! *** evaluate the HYPERGEOMETRIC FUNCTION 1F1                          
!                                        i                              
!            F (AA;BB; Z) = SUM  (AA)   Z / ( (BB)  i] )                
!           1 1              i       i            i                     
!                                                                       
!     to accuracy EPS with at most LIMIT terms.                         
!  If KIND = 0 : using extended precision but real arithmetic only,     
!            1 : using normal precision in complex arithmetic,          
!   or       2 : using normal complex arithmetic, but with CDIGAM factor
!-----------------------------------------------------------------------
!  where                                                                
      AA = ZL+ONE - ETA*P                                            
      BB = TWO*(ZL+ONE)                                              
!  and                                                                  
      Z  = TWO*P*X                                                   
!-----------------------------------------------------------------------
      ZLLIN = REAL(BB).LE.ZERO .AND. ABS(BB-NINTC(BB)).LT.ACC8**0.25 
      IF(.NOT.ZLLIN.OR.REAL(BB)+LIMIT.LT.1.5) GO TO 10           
      NITS = -1                                               
      RETURN                                                  
  10  IF (LIMIT.LE.0) THEN                                               
        F11 = ZERO                                                     
        ERR = ZERO                                                     
        NITS= 1                                                        
        RETURN                                                         
      END IF                                                          
      TA = ONE                                                          
      RK = ONE                                                          
      IF (KIND.LE.0.AND.ABSC(Z)*ABSC(AA).GT.ABSC(BB)*1.0) THEN         
        DR = ONE                                                       
        DI = ZERO                                                      
        GR = ONE                                                       
        GI = ZERO                                                      
        AR = REAL(AA)                                                  
        BR = REAL(BB)                                                  
        FI = ZERO                                                      
        DO 20 I=2,LIMIT                                                   
          FI1 = FI + ONE                                                 
          TR = BR * FI1                                                  
          TI = IMAG(BB) * FI1                                            
          DEN= ONE / (TR*TR + TI*TI)                                     
          UR = (AR*TR + IMAG(AA)*TI) * DEN                               
          UI = (IMAG(AA)*TR - AR*TI) * DEN                               
          TR = UR*GR - UI*GI                                             
          TI = UR*GI + UI*GR                                             
          GR = REAL(Z) * TR - IMAG(Z)*TI                                 
          GI = REAL(Z) * TI + IMAG(Z)*TR                                 
          DR = DR + GR                                                   
          DI = DI + GI                                                   
          ERR = ABS(GR) + ABS(GI)                                     
          IF (ERR.GT.FPMAX) GO TO 60                                
          RK  = ABS(DR) + ABS(DI)                                     
          TA = MAX(TA,RK)                                             
          IF (ERR.LT.RK*EPS .OR. I.GE.4.AND.ERR.LT.ACC16) GO TO 30        
          FI = FI1                                                       
          AR = AR + ONE                                                  
  20    BR = BR + ONE                                                  
!-----------------------------------------------------------------------
  30    F11 = DR + CI * DI                                             
        ERR = ACC16 * TA / RK                                          
!-----------------------------------------------------------------------
      ELSE                                                              
!-----------------------------------------------------------------------
!* -------------------alternative code----------------------------------
!-----------------------------------------------------------------------
!*    If REAL*16 arithmetic is not available, (or already using it]),   
!*    then use KIND > 0                                                 
!-----------------------------------------------------------------------
        G = ONE                                                        
        F = ONE                                                       
        IF (KIND.GE.2) F = CDIGAM(AA) - CDIGAM(BB) - CDIGAM(G)         
        DD = F                                                         
        DO 40 I=2,LIMIT                                                
          AI = AA + (I-2)                                             
          BI = BB + (I-2)                                             
          R  = I-ONE                                                  
          G = G * Z * AI / (BI * R)                                      
!-----------------------------------------------------------------------
!  multiply by (psi(a+r)-psi(b+r)-psi(1+r))                             
!-----------------------------------------------------------------------
          IF (KIND.GE.2) F = F + ONE/AI - ONE/BI - ONE/R                 
          T  = G * F                                                     
          DD = DD + T                                                    
          ERR = ABSC(T)                                               
          IF (ERR.GT.FPMAX) GO TO 60                                
          RK = ABSC(DD)                                               
          TA = MAX(TA,RK)                                                
          IF (ERR.LT.RK*EPS.OR.ERR.LT.ACC8.AND.I.GE.4) GO TO 50           
  40    CONTINUE                                                       
!-----------------------------------------------------------------------
  50    ERR = ACC8 * TA / RK                                           
        F11 = DD                                                       
!* --------------------------------------------- end of alternative code
      END IF                                                             
  60  NITS = I                                                       
!-----------------------------------------------------------------------
    END FUNCTION F11                                                              
!-----------------------------------------------------------------------
    REAL(dpf) FUNCTION CF1R(X,ETA,ZL,EPS,FCL,TPK1,ETANE0,LIMIT,ERR,NFP &       
                         &,  ACCH,FPMIN,FPMAX,PR,CALLER)                         
!-----------------------------------------------------------------------
      IMPLICIT REAL(dpf)(A-H,O-Z)                                          
      LOGICAL PR,ETANE0                                                 
      CHARACTER(len=6) CALLER                                                
      DATA ONE,TWO / 1D+0, 2D+0 /                                       
!-----------------------------------------------------------------------
 1000 FORMAT(/' ',A6,': CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,ETA,X = ',/1X,1P,7D9.2/)
 1010 FORMAT(' ',A6,': CF1 HAS FAILED TO CONVERGE AFTER ',I10  ,' ITERATIONS AS ABS(X) =',F15.0)
!-----------------------------------------------------------------------
! ***    Evaluate CF1  =  F   =  F'(ZL,ETA,X)/F(ZL,ETA,X)               
!        using real arithmetic                                          
!-----------------------------------------------------------------------
      FCL = ONE                                                         
      XI = ONE/X                                                        
      PK  = ZL + ONE                                                    
      PX  = PK  + LIMIT                                                 
  10  EK  = ETA / PK                                                    
      RK2 = ONE + EK*EK                                      
      F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI                           
      PK1 =  PK + ONE                                                   
      TPK1 = PK + PK1                                                
      TK  = TPK1*(XI + EK/PK1)                                          
      IF (ETANE0) THEN                                                   
!-----------------------------------------------------------------------
! ***   test ensures b1 .ne. zero for negative ETA etc.; fixup is exact.
!-----------------------------------------------------------------------
        IF (ABS(TK) .GT. ACCH)  GO TO 20                            
        FCL  = RK2/(ONE + (ETA/PK1)**2)                            
        SL   = TPK1*XI * (TPK1+TWO)*XI                             
        PK   =  TWO + PK                                           
        GO TO 10                                                   
      END IF                                                          
  20  D =  ONE/TK                                                     
      DF= -FCL*RK2*D                                                  
      IF (PK.GT.ZL+TWO) FCL = - RK2 * SL                           
      FCL=FCL * D * TPK1 * XI                                   
      F  = F  + DF                                                    
!-----------------------------------------------------------------------
! ***   begin CF1 loop on PK = k = lambda + 1                           
!-----------------------------------------------------------------------
      RK    = ONE                                                       
      SMALL = SQRT(FPMIN)                                            
  30  PK    = PK1                                                       
      PK1   = PK1 + ONE                                                 
      TPK1 = PK + PK1                                                
      IF (ETANE0) THEN                                                
        EK  = ETA / PK                                               
        RK2 = ONE + EK*EK                                   
      END IF                                                         
      TK=TPK1*(XI + EK/PK1)                                        
      D = TK - D*RK2                                               
      IF (ABS(D) .GT. ACCH) GO TO 40                 
      IF (PR) WRITE (6,1000) CALLER,D,DF,ACCH,PK,EK,ETA,X        
      RK= RK + ONE                                            
      IF( RK .GT. TWO ) GO TO 50               
  40  D = ONE/D                                                     
      FCL = FCL * D * TPK1*XI                                     
      IF (ABS(FCL).LT.SMALL) FCL = FCL / SMALL                     
      IF (ABS(FCL).GT.FPMAX) FCL = FCL / FPMAX                     
      DF  = DF*(D*TK - ONE)                                           
      F   = F  + DF                                                   
      IF( PK .GT. PX ) GO TO 50                                 
      IF (ABS(DF) .GE. ABS(F)*EPS) GO TO 30                  
      NFP = PK - ZL - 1                                       
      ERR = EPS * SQRT(REAL(NFP))                           
      CF1R = F                                                          
!-----------------------------------------------------------------------
      RETURN                                                            
!-----------------------------------------------------------------------
  50  IF (PR) WRITE (6,1010) CALLER,LIMIT,ABS(X)                         
      ERR = TWO                                                         
!-----------------------------------------------------------------------
    END FUNCTION CF1R                                                   
!-----------------------------------------------------------------------
    COMPLEX(dpf) FUNCTION F20(AA,BB,Z,EPS,JMAX,RE,FPMAX,N,X)            
!-----------------------------------------------------------------------
!     evaluate the HYPERGEOMETRIC FUNCTION 2F0                          
!-----------------------------------------------------------------------
!            F (AA,BB;;Z) = SUM  (AA)  (BB)  Z / i]                     
!           2 0              i       i     i                            
!-----------------------------------------------------------------------
!     to accuracy EPS with at most JMAX terms.                          
!-----------------------------------------------------------------------
!     if the terms start diverging,                                     
!     the corresponding continued fraction is found by RCF              
!     & evaluated progressively by Steed's method to obtain convergence.
!-----------------------------------------------------------------------
!      useful number also input:  FPMAX = near-largest f.p. number      
!-----------------------------------------------------------------------
      IMPLICIT COMPLEX(dpf)(A-H,O-Z)                                      
      DIMENSION X(JMAX,4)                                               
      LOGICAL FINITE                                                    
      REAL(dpf) EP,EPS,AT,ATL,ABSC,RE,FPMAX                                
      DATA ONE,ZERO / (1D+0,0D+0), (0D+0,0D+0) /                        
!-----------------------------------------------------------------------
      ABSC(W) = ABS(REAL(W)) + ABS(IMAG(W))                             
      NINTC(W) = NINT(REAL(REAL(W)))                                    
!-----------------------------------------------------------------------
      RE = 0.0                                                          
      X(1,1) = ONE                                                      
      SUM = X(1,1)                                                      
      ATL = ABSC(X(1,1))                                                
      F = SUM                                                     
      D = ONE                                                        
      DF= SUM                                                     
      J = 0                                                             
      EP = EPS * JMAX *10.                                              
      MA = - NINTC(AA)                                                  
      MB = - NINTC(BB)                                                  
      FINITE = ABS(ABS(REAL(AA))-MA).LT.EP .AND. ABS(IMAG(AA)).LT.EP &  
     &    .OR. ABS(ABS(REAL(BB))-MB).LT.EP .AND. ABS(IMAG(BB)).LT.EP    
      IMAX = JMAX                                                       
      IF (FINITE.AND.MA.GE.0) IMAX = MIN(MA+1,IMAX)                      
      IF (FINITE.AND.MB.GE.0) IMAX = MIN(MB+1,IMAX)                      
      DO 10 I=2,IMAX                                                    
        X(I,1) = X(I-1,1) * Z * (AA+I-2) * (BB+I-2) / (I-1)               
        IF (ABSC(X(I,1)).GT.FPMAX) GO TO 40                             
        AT = ABSC(X(I,1))                                                 
        IF (J.EQ.0) THEN                                                
          SUM = SUM + X(I,1)                                     
          IF (AT .LT. ABSC(SUM)*EPS) GO TO 20                     
        END IF                                                    
        IF (FINITE) GO TO 10                                               
        IF (J.GT.0 .OR. AT.GT.ATL .OR. I.GE.JMAX-2) J = J + 1              
        IF (J.EQ.0) GO TO 10                                            
        CALL RCF(X(1,1),X(1,2),J,I,X(1,3),EPS)                         
        IF (I.LT.0) GO TO 40                                       
        DO 50 K=MAX(J,2),I                                          
          D = ONE/(D*X(K,2) + ONE)                                    
          DF = DF*(D - ONE)                                           
          F = F + DF                                                  
          IF (ABSC(DF) .LT. ABSC(F)*EPS) GO TO 30                      
          IF (DF.EQ.ZERO.AND.F.EQ.ZERO.AND.I.GE.4) GO TO 30            
    50  CONTINUE                                                    
!-----------------------------------------------------------------------
        J = I                                                          
  10  ATL = AT                                                          
!-----------------------------------------------------------------------
      IF(.NOT.FINITE) I = -JMAX                                         
  20  N = I                                                             
      F20 = SUM                                                        
      IF(.NOT.FINITE) RE  = AT / ABSC(SUM)                             
      RETURN                                                           
!-----------------------------------------------------------------------
  30  F20 = F                                                           
      RE = ABSC(DF) / ABSC(F)                                           
      N = K                                                             
      RETURN                                                            
!-----------------------------------------------------------------------
  40  I = 0                                                             
      GO TO 20                                                          
!-----------------------------------------------------------------------
    END FUNCTION F20                                                  
!-----------------------------------------------------------------------
    COMPLEX(dpf) FUNCTION CF1A(RHO,ETA,XL,PSI,EPS,NMAX,NUSED,FCL,RE,FPMAX,XX,G,C)
!-----------------------------------------------------------------------
!     evaluate the ASYMPTOTIC EXPANSION for the                         
!            LOGARITHMIC DERIVATIVE OF THE REGULAR SOLUTION             
!-----------------------------------------------------------------------
! ***        CF1A  =  f   =  F'(XL,ETA,RHO)/F(XL,ETA,RHO)               
!-----------------------------------------------------------------------
!      that is valid for REAL(RHO)>0, and best for RHO >> ETA**2, XL,   
!      and is derived from the 2F0 expansions for H+ and H-             
!      e.g. by Froeberg (Rev. Mod. Physics Vol 27, p399 , 1955)         
!      Some lines of this subprogram are for convenience copied from    
!           Takemasa, Tamura & Wolter CPC 17 (1979) 351.                
!-----------------------------------------------------------------------
!     Evaluate to accuracy EPS with at most NMAX terms.                 
!-----------------------------------------------------------------------
!     If the terms start diverging,                                     
!     the corresponding continued fraction is found by RCF              
!     & evaluated progressively by Steed's method to obtain convergence.
!-----------------------------------------------------------------------
!      useful number also input:  FPMAX = near-largest f.p. number      
!-----------------------------------------------------------------------
      IMPLICIT COMPLEX(dpf)(A-H,O-Z)                                      
      DIMENSION XX(2,NMAX),G(NMAX),C(NMAX)                              
      REAL(dpf) RE,EPS,T1,T2,T3,ZERO,ONE,TWO,AT,ATL,ABSC,FPMAX             
      DATA ZERO,ONE,TWO,CI / 0D+0, 1D+0, 2D+0, (0D+0,1D+0) /            
!-----------------------------------------------------------------------
      ABSC(W) = ABS(REAL(W)) + ABS(IMAG(W))                             
!-----------------------------------------------------------------------
      HPI = TWO*ATAN(ONE)                                               
      T1 = SIN(REAL(PSI))                                               
      T2 = COS(REAL(PSI))                                               
      ATL= TANH(IMAG(PSI))                                              
!-----------------------------------------------------------------------
! GIVE COS(PSI)/COSH(IM(PSI)), WHICH ALWAYS HAS CORRECT SIGN
!-----------------------------------------------------------------------
      COSL = DCMPLX( T2 , -T1 * ATL )                               
      TANL = DCMPLX(T1,T2*ATL) / COSL                                   
      RE = ZERO                                                         
      XLL1= XL*(XL+ONE)                                                 
      ETASQ = ETA*ETA                                                   
      SL1=ONE                                                           
      SL=SL1                                                            
      SC1=ZERO                                                          
      SC=SC1                                                            
      TL1=SC                                                            
      TL=TL1                                                            
      TC1=ONE-ETA/RHO                                                   
      TC=TC1                                                            
      FCL  = TL + SL*TANL                                               
      G(1) = (TC + SC*TANL) / FCL                                       
      GLAST = G(1)                                                      
      ATL = ABSC(GLAST)                                                 
      F    = GLAST                                                   
      D = ONE                                                        
      DF   = GLAST                                                   
      J = 0                                                             
      DO 10 N=2,NMAX                                                    
        T1=N-1                                                            
        T2=TWO*T1-ONE                                                     
        T3=T1*(T1-ONE)                                                    
        DENOM=TWO*RHO*T1                                                  
        C1=(ETA*T2)/DENOM                                                 
        C2=(ETASQ+XLL1-T3)/DENOM                                          
        SL2=C1*SL1-C2*TL1                                                 
        TL2=C1*TL1+C2*SL1                                                 
        SC2=C1*SC1-C2*TC1-SL2/RHO                                         
        TC2=C1*TC1+C2*SC1-TL2/RHO                                         
        SL=SL+SL2                                                         
        TL=TL+TL2                                                         
        SC=SC+SC2                                                         
        TC=TC+TC2                                                         
        SL1=SL2                                                           
        TL1=TL2                                                           
        SC1=SC2                                                           
        TC1=TC2                                                           
        FCL  =  TL + SL*TANL                                              
        IF (ABSC(FCL).GT.FPMAX .OR. ABSC(FCL).LT.1./FPMAX) GO TO 40     
        GSUM = (TC + SC*TANL) / FCL                                       
        G(N) = GSUM - GLAST                                               
        GLAST = GSUM                                                      
        AT = ABSC(G(N))                                                
        IF (AT.LT.ABSC(GSUM)*EPS) GO TO 20                              
        IF (J.GT.0 .OR. AT.GT.ATL .OR. N.GE.NMAX-2) J = J + 1              
        IF (J.EQ.0) GO TO 10                                            
        CALL RCF(G,C,J,N,XX,EPS)                                    
        IF (N.LT.0) GO TO 40                                       
        DO 60 K=MAX(J,2),N                                          
          D = ONE/(D*C(K) + ONE)                                   
          DF = DF*(D - ONE)                                        
          F = F + DF                                               
          IF (ABSC(DF) .LT. ABSC(F)*EPS) GO TO 30                         
          IF (DF.EQ.ZERO.AND.F.EQ.ZERO.AND.N.GE.4) GO TO 30               
  60    CONTINUE                                                  
        J = N                                                          
  10  ATL = AT                                                       
      K = -NMAX                                                         
      GO TO 30                                                          
!-----------------------------------------------------------------------
  20  FCL = FCL * COSL                                                  
      CF1A = GSUM                                                    
      RE = AT / ABSC(GSUM)                                           
      NUSED = N                                                      
      RETURN                                                         
!-----------------------------------------------------------------------
  30  CF1A = F                                                          
      FCL = FCL * COSL                                                  
      RE = ABSC(DF) / ABSC(F)                                        
      NUSED = K                                                      
      RETURN                                                            
!-----------------------------------------------------------------------
  40  CF1A = G(1)                                                       
      FCL = 1.0                                                         
      RE = 1.0                                                          
      NUSED = 0                                                         
!-----------------------------------------------------------------------
    END FUNCTION CF1A                                                 
!-----------------------------------------------------------------------
    SUBROUTINE RCF(A,B,IBEG,INUM,XX,EPS)                              
!-----------------------------------------------------------------------
!  RCF converts polynomial A to the corresponding continued             
!         fraction, in 'normal'  form with coefficients B               
!         by the 'P algorithmn' of Patry & Gupta                        
!-----------------------------------------------------------------------
!   A(z) = A1/z + A2/z**3 + A3/z**5 + ... + An/z**(2n-1)                
!-----------------------------------------------------------------------
!   B(z) = B1/z+ B2/z+ B3/z+ .../(z+ Bn/z)                              
!-----------------------------------------------------------------------
!  data:                                                                
!   A     vector A(k), k=1,INUM         input                           
!   B     vector B(k), k=IBEG,INUM      output                          
!   IBEG  order of first coef. calc.    input                           
!   INUM  order of A, even or odd       input                           
!   XX    auxiliary vector of length .ge. length of vector B            
!         caller provides space for A,B,XX                              
!     Note that neither of the first two terms A(1) A(2) should be zero 
!             & the user can start the calculation with any value of    
!                IBEG provided the c.f. coefs have been already         
!                calculated up to INUM = IBEG-1                         
!             & the method breaks down as soon as the absolute value    
!                of a c.f. coef. is less than EPS.    At the time of the
!                break up XX(1) has been replaced by 1E-50, and INUM has
!                been replaced by minus times the number of this coef.  
!   algorithm: J.Patry & S.Gupta,                                       
!              EIR-bericht nr. 247,                                     
!              Eidg. Institut fur Reaktorforschung Wuerenlingen         
!              Wueringlingen, Schweiz.                                  
!              November 1973                                            
!   see also:  Haenggi,Roesel & Trautmann,                              
!              Jnl. Computational Physics, vol 137, pp242-258 (1980)    
!   note:      restart procedure modified by I.J.Thompson               
!-----------------------------------------------------------------------
      IMPLICIT COMPLEX(dpf)(A-H,O-Z)                                      
      DIMENSION A(100),B(100),XX(2,100)                                 
      LOGICAL EVEN                                                      
      REAL(dpf) EPS                                                        
!-----------------------------------------------------------------------
      COMMON /RCFCM2/ X1,M2M1,MP12,EVEN,M                               
!-----------------------------------------------------------------------
!     ibn = ibeg + inum - 1                                             
!-----------------------------------------------------------------------
      IBN = INUM                                                        
!-----------------------------------------------------------------------
!     B(IBN) is last value set on this call                             
!-----------------------------------------------------------------------
      IF (IBEG.GT.4 .AND. M .NE. IBEG-1) GO TO 90                        
!-----------------------------------------------------------------------
!     B(M) is last value set in previous call                           
!-----------------------------------------------------------------------
      IF (IBEG.GT.4) GO TO 50                                            
      IF (IBEG.EQ.4) GO TO 20                                            
      B(1) = A(1)                                                       
      IF (IBN.GE.2) B(2) = - A(2)/A(1)                                   
      IF (IBN.LT.3) GO TO 10                                             
      X0 = A(3) / A(2)                                                  
      XX(2,1) = B(2)                                                    
      XX(1,1) = - X0                                                    
      XX(1,2) = 0.                                                      
      B(3) = -X0 - B(2)                                                 
      X0 = -B(3) * A(2)                                                 
      M = 3                                                             
      MP12 = 2                                                          
      EVEN = .TRUE.                                                     
      IF (IBN.GT.3) GO TO 20                                             
!-----------------------------------------------------------------------
  10  RETURN                                                            
!-----------------------------------------------------------------------
  20  IF (ABS(B(3)) .LT. EPS*ABS(X0)) GO TO 80                            
      M = 4                                                             
  30  X1 = A(M)                                                         
      M2M1 = MP12                                                       
      MP12 = M2M1 + 1                                                   
      IF (EVEN) MP12 = M2M1                                              
      DO 40 K=2,MP12                                                    
  40  X1 = X1 + A(M-K+1) * XX(1,K-1)                                    
      B(M) = - X1/X0                                                    
!-----------------------------------------------------------------------
      IF (M.GE.IBN) RETURN                                               
!-----------------------------------------------------------------------
  50  IF (ABS(B(M)).LT.EPS*ABS(X0)) GO TO 80                             
      K = M2M1                                                          
  60  XX(2,K) = XX(1,K) + B(M) * XX(2,K-1)                              
      K = K-1                                                           
      IF (K.GT.1) GO TO 60                                               
      XX(2,1) = XX(1,1) + B(M)                                          
      DO 70 K=1,M2M1                                                    
        X0 = XX(2,K)                                                      
        XX(2,K) = XX(1,K)                                                 
  70  XX(1,K) = X0                                                      
      X0 = X1                                                           
      XX(1,M2M1+1) = 0.                                                 
      M = M+1                                                           
      EVEN = .NOT.EVEN                                                  
      GO TO 30                                                          
  80  INUM = -M                                                         
!-----------------------------------------------------------------------
      RETURN                                                            
!-----------------------------------------------------------------------
 1000 FORMAT('0RCF: LAST CALL SET M =',I4,', BUT RESTART REQUIRES',I4)  
   90 WRITE(STDOUT,1000) M,IBEG-1                                       
      STOP ("RCF HAS FAILED")                                           
!-----------------------------------------------------------------------
    END SUBROUTINE RCF                                                  
!-----------------------------------------------------------------------
!   this routine computes the logarithm of the gamma function gamma(z)
!   for any complex argument 'Z' to any accuracy preset by CALL LOGAM 
!-----------------------------------------------------------------------
    COMPLEX(dpf) FUNCTION CLOGAM(Z)                         
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(dpf),   PARAMETER :: ACC8=(2._dpf)*(10._dpf)**(-16_dpf)  !2D-16
      REAL(dpf),   PARAMETER :: ACC16=(3._dpf)*(10._dpf)**(-33_dpf) !3D-33_dpf
      REAL(dpf),   PARAMETER :: FPLMIN=-140._dpf 
      REAL(dpf),   PARAMETER :: ZERO=0._dpf
      REAL(dpf),   PARAMETER :: HALF=0.5_dpf
      REAL(dpf),   PARAMETER :: QUART=0.25_dpf
      REAL(dpf),   PARAMETER :: ONE=1._dpf
      REAL(dpf),   PARAMETER :: TWO=2._dpf
      REAL(dpf),   PARAMETER :: FOUR=4._dpf
      REAL(dpf),   PARAMETER :: PI=FOUR*ATAN(ONE)
      REAL(dpf),   PARAMETER :: ALPI = LOG(PI)
      REAL(dpf),   PARAMETER :: HL2P = LOG(TWO*PI)*HALF
      INTEGER(spi),PARAMETER :: NB=15
      REAL(dpf),DIMENSION(NB) :: B
      COMPLEX(dpf) :: Z,V,H,R,SER
      INTEGER(spi) :: NX0,N,J,K,NT,I,MX
      REAL(dpf) :: ACCUR,ACC,X0,X,A,C,D,E,F,T,F21,ERR
      LOGICAL :: ACCUR_REACHED
!-----------------------------------------------------------------------
      REAL(dpf),DIMENSION(NB),PARAMETER :: BN= [             +1._dpf &
                                               &,            -1._dpf &
                                               &,            +1._dpf &
                                               &,            -1._dpf &
                                               &,            +5._dpf &
                                               &,          -691._dpf &
                                               &,          +  7._dpf &
                                               &,         -3617._dpf &
                                               &,         43867._dpf &
                                               &,       -174611._dpf &
                                               &,        854513._dpf &
                                               &,    -236364091._dpf &
                                               &,     + 8553103._dpf &
                                               &,  -23749461029._dpf &
                                               &, 8615841276005._dpf] 
      REAL(dpf),DIMENSION(NB),PARAMETER :: BD= [              6._dpf &
                                               &,            30._dpf &
                                               &,            42._dpf &
                                               &,            30._dpf &
                                               &,            66._dpf &
                                               &,          2730._dpf &
                                               &,             6._dpf &
                                               &,           510._dpf &
                                               &,           798._dpf &
                                               &,           330._dpf &
                                               &,           138._dpf &
                                               &,          2730._dpf &
                                               &,             6._dpf &
                                               &,           870._dpf &
                                               &,         14322._dpf]
!-----------------------------------------------------------------------
 1000 FORMAT(1X,A6,' ... ARGUMENT IS NON POSITIVE INTEGER = ',F20.2)    
!-----------------------------------------------------------------------
! INLINE INIT FOR NOW (logam entry)
!-----------------------------------------------------------------------
      ACC=ACC8
!-----------------------------------------------------------------------
      NX0 = 6                                                           
      X0  = NX0 + ONE                                                   
      ACCUR = ACC                                                       
      ACCUR_REACHED=.FALSE.
      LOOP_120: DO K=1,NB                                                     
        F21 = K*2 - ONE                                                  
        B(K) = BN(K) / (BD(K) * K*TWO * F21)                             
        ERR = ABS(B(K)) * K*TWO / X0**F21                                
        IF (ERR.LT.ACC) THEN
          ACCUR_REACHED=.TRUE.
          EXIT LOOP_120
        END IF                                          
      END DO LOOP_120                                         
      IF (.NOT.ACCUR_REACHED) THEN
        NX0 = INT((ERR/ACC)**(ONE/F21) * X0)                             
        K = NB                                                           
      END IF
      NT = K                                                            
!-----------------------------------------------------------------------
      X=Z%RE                                                         
      T=Z%IM                                                         
      MX=INT(REAL(ACCUR*100._dpf-X,KIND=dpf),KIND=spi)                                     
      IF (ABS(ABS(X)-MX)+ABS(T).LT.ACCUR*50) THEN
        WRITE(STDOUT,1000) 'CLOGAM',X                                       
        CLOGAM = ZERO                                                     
      ELSE                  
        F=ABS(T)                                                          
        V=CMPLX(X,F,KIND=dpf)                                                     
        IF (X.LT.ZERO) V=ONE-V                                           
        H=ZERO                                                            
        C=V%RE                                                         
        N=NX0-INT(C,KIND=spi)                                                      
        IF (N.GE.0_spi) THEN                                             
          H=V                                                               
          D=V%IM                                                         
          A=ATAN2(D,C)                                                      
          IF (N.NE.0) THEN                                             
            DO I = 1,N                                                     
              C=C+ONE                                                           
              V=CMPLX(C,D,KIND=dpf)                                                     
              H=H*V                                                             
              A=A+ATAN2(D,C)
            END DO
          END IF                                                    
          H=CMPLX(HALF*LOG(H%RE**2+H%IM**2),A,KIND=dpf)                       
          V=V+ONE                                                           
        END IF
        R=ONE/V**2                                                        
        SER = B(NT)                                                       
        DO J=2,NT                                                      
          K = NT+1 - J                                                    
          SER = B(K) + R*SER
        END DO                                                
!-----------------------------------------------------------------------
        CLOGAM = HL2P+(V-HALF)*LOG(V)-V + SER/V - H                       
!-----------------------------------------------------------------------
        IF (X.LT.ZERO) THEN                                          
!-----------------------------------------------------------------------
          A= INT(X,KIND=spi)-ONE                                                     
          C=PI*(X-A)                                                        
          D=PI*F                                                            
          E = ZERO                                                        
          F = -TWO*D                                                      
          IF (F.GT.FPLMIN) E = EXP(F)                                      
          F=SIN(C)                                                          
          E= D + HALF*LOG(E*F**2+QUART*(ONE-E)**2)                          
          F=ATAN2(COS(C)*TANH(D),F)-A*PI                                    
          CLOGAM=ALPI-CMPLX(E,F,KIND=dpf)-CLOGAM
        END IF                                  
!-----------------------------------------------------------------------
        IF (SIGN(ONE,T).LT.-HALF) CLOGAM=CONJG(CLOGAM)                 
      END IF
!-----------------------------------------------------------------------
    END FUNCTION CLOGAM                                                 
!-----------------------------------------------------------------------
!   this routine computes the logarithmic derivative of the gamma     
!   function  psi(Z) = digamma(Z) = d (ln gamma(Z))/dZ  for any       
!   complex argument Z, to any accuracy preset by CALL LOGAM(ACC)     
!-----------------------------------------------------------------------
    COMPLEX(dpf) FUNCTION CDIGAM(Z)                                     
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(dpf),   PARAMETER :: ACC8=(2._dpf)*(10._dpf)**(-16_dpf)  !2D-16
      REAL(dpf),   PARAMETER :: ACC16=(3._dpf)*(10._dpf)**(-33_dpf) !3D-33_dpf
      REAL(dpf),   PARAMETER :: ZERO=0._dpf
      REAL(dpf),   PARAMETER :: HALF=0.5_dpf
      REAL(dpf),   PARAMETER :: QUART=0.25_dpf
      REAL(dpf),   PARAMETER :: ONE=1._dpf
      REAL(dpf),   PARAMETER :: TWO=2._dpf
      REAL(dpf),   PARAMETER :: FOUR=4._dpf
      REAL(dpf),   PARAMETER :: PI=FOUR*ATAN(ONE)
      INTEGER(spi),PARAMETER :: NB=15
      REAL(dpf),DIMENSION(NB) :: B
      COMPLEX(dpf) :: Z,U,V,H,R,SER
      INTEGER(spi) :: NX0,N,J,K,NT,I
      REAL(dpf) :: ACCUR,ACC,X0,X,A,F21,ERR
      LOGICAL :: ACCUR_REACHED
!-----------------------------------------------------------------------                           
      REAL(dpf),DIMENSION(NB),PARAMETER :: BN= [             +1._dpf &
                                               &,            -1._dpf &
                                               &,            +1._dpf &
                                               &,            -1._dpf &
                                               &,            +5._dpf &
                                               &,          -691._dpf &
                                               &,          +  7._dpf &
                                               &,         -3617._dpf &
                                               &,         43867._dpf &
                                               &,       -174611._dpf &
                                               &,        854513._dpf &
                                               &,    -236364091._dpf &
                                               &,     + 8553103._dpf &
                                               &,  -23749461029._dpf &
                                               &, 8615841276005._dpf] 
      REAL(dpf),DIMENSION(NB),PARAMETER :: BD= [              6._dpf &
                                               &,            30._dpf &
                                               &,            42._dpf &
                                               &,            30._dpf &
                                               &,            66._dpf &
                                               &,          2730._dpf &
                                               &,             6._dpf &
                                               &,           510._dpf &
                                               &,           798._dpf &
                                               &,           330._dpf &
                                               &,           138._dpf &
                                               &,          2730._dpf &
                                               &,             6._dpf &
                                               &,           870._dpf &
                                               &,         14322._dpf]
!-----------------------------------------------------------------------
 1000 FORMAT(1X,A6,' ... ARGUMENT IS NON POSITIVE INTEGER = ',F20.2)    
!-----------------------------------------------------------------------
! INLINE INIT FOR NOW (logam entry)
!-----------------------------------------------------------------------
      ACC=ACC8
!-----------------------------------------------------------------------
      NX0 = 6                                                           
      X0  = NX0 + ONE                                                   
      ACCUR = ACC                                                       
      ACCUR_REACHED=.FALSE.
      LOOP_120: DO K=1,NB                                                     
        F21 = K*2 - ONE                                                  
        B(K) = BN(K) / (BD(K) * K*TWO * F21)                             
        ERR = ABS(B(K)) * K*TWO / X0**F21                                
        IF (ERR.LT.ACC) THEN
          ACCUR_REACHED=.TRUE.
          EXIT LOOP_120
        END IF                                          
      END DO LOOP_120                                         
      IF (.NOT.ACCUR_REACHED) THEN
        NX0 = INT((ERR/ACC)**(ONE/F21) * X0)                             
        K = NB                                                           
      END IF
      NT = K                                                            
!-----------------------------------------------------------------------
      U=Z                                                               
      X=U%RE                                                         
      A=ABS(X)                                                          
      IF (ABS(U%IM) + ABS(A+INT(X,KIND=spi)).LT.ACCUR) THEN              
!-----------------------------------------------------------------------
        WRITE(STDOUT,1000) 'CDIGAM',X                                     
        CDIGAM=ZERO
!-----------------------------------------------------------------------
      ELSE
        IF (X.LT.ZERO) U=-U                                              
        V=U                                                               
        H=ZERO                                                            
        N=NX0-INT(A,KIND=spi)                                             
        IF (N.GE.0_spi) THEN                                             
          H=ONE/V                                                           
          IF (N.NE.0) THEN                                             
            DO I = 1,N                                                     
              V=V+ONE                                                           
              H=H+ONE/V                                                         
            END DO
          END IF                                                        
          V=V+ONE                                                           
        END IF
        R=ONE/V**2                                                        
        SER = B(NT) * (2*NT-1)                                            
        DO J=2,NT                                                     
          K = NT+1 - J                                                    
          SER = B(K)*(2*K-1) + R*SER
        END DO                                        
!-----------------------------------------------------------------------
        CDIGAM = LOG(V) - HALF/V - R*SER - H                              
!-----------------------------------------------------------------------
        IF (X.LT.ZERO) THEN                                            
          H=PI*U                                                            
          CDIGAM = CDIGAM + ONE/U + PI*COS(H)/SIN(H) 
        END IF                     
!-----------------------------------------------------------------------
      END IF                                                            
!-----------------------------------------------------------------------
    END FUNCTION CDIGAM
!-----------------------------------------------------------------------
!   TIDY A COMPLEX NUMBER                                             
!-----------------------------------------------------------------------
    COMPLEX(dpf) FUNCTION TIDY(Z,ACC)                                   
      IMPLICIT NONE
!-----------------------------------------------------------------------
      REAL(dpf) :: X,Y,ACC,AZ                                              
      COMPLEX(dpf) :: Z                                                 
      REAL(dpf),PARAMETER :: ZERO=0._dpf
!-----------------------------------------------------------------------
      X = Z%RE                                                       
      Y = Z%IM                                                       
      AZ= (ABS(X) + ABS(Y)) * ACC * 5                                   
      IF (ABS(X) .LT. AZ) X = ZERO                                       
      IF (ABS(Y) .LT. AZ) Y = ZERO                                       
      TIDY = CMPLX(X,Y,KIND=dpf)                                                
!-----------------------------------------------------------------------
    END FUNCTION TIDY                                                   
!-----------------------------------------------------------------------
END MODULE COULCC_M
!-----------------------------------------------------------------------
PROGRAM CCTEST
!-----------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only: spi=>int32, dpf=>real64 &
                              &,stdin=>input_unit,stdout=>output_unit
!-----------------------------------------------------------------------
  USE COULCC_M                                               
!-----------------------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------------------
      LOGICAL :: WHIT                                                      
      INTEGER(spi) :: I,L,NL,MODE,IFAIL,MD,IH,KFN,KFIN
      COMPLEX(dpf) :: X,ETA,ZLMIN,ZL,WS
      COMPLEX(dpf),DIMENSION(201) :: FC,GC,FCP,GCP,SIG               
!-----------------------------------------------------------------------
      REAL(dpf) :: RERR
      INTEGER(spi),DIMENSION(2) :: NFP,KASE
      INTEGER(spi),DIMENSION(3) :: NPQ
      COMMON /STEED/ RERR,NFP,NPQ,KASE                                  
!-----------------------------------------------------------------------
      CHARACTER(len=20) :: NOTE                                         
      CHARACTER(len=4)  :: IRREG,REG                                  
!-----------------------------------------------------------------------
      INTEGER(spi),PARAMETER,DIMENSION(8) :: LDISP=[1,2,4,11,31,101,301,1001]
      REAL(dpf),   PARAMETER :: ZERO=0._dpf
      REAL(dpf),   PARAMETER :: HALF=0.5_dpf
      REAL(dpf),   PARAMETER :: ONE=1._dpf
      REAL(dpf),   PARAMETER :: FOUR=4._dpf
      COMPLEX(dpf),PARAMETER :: CI=CMPLX(0._dpf,1._dpf,KIND=dpf)
      REAL(dpf),   PARAMETER :: PI=FOUR*ATAN(ONE)
!-----------------------------------------------------------------------
      CHARACTER(len=4),DIMENSION(2,4,3) :: WHO                          
      DATA WHO / 'F','G','j','y','J','Y','I','K' ,         &           
       &         'F','H+','j','h(1)','J','H(1)','?','?' ,  &           
       &         'F','H-','j','h(2)','J','H(2)','?','?' /               
!-----------------------------------------------------------------------
 1000 FORMAT('1TEST OF THE CONTINUED-FRACTION COULOMB & BESSEL ROUTINES'/)
 1010 FORMAT(/'0X =',F12.4,F9.4,', ETA =',2F8.3,', ZLMIN =',2F8.3,'  NL=',I4,  '  MODE = ',I3,'  KFN =',I2,8X,A20)
 1020 FORMAT(' COULCC  :: IFAIL =',I4,'  RERR =',1P,E10.1,'  ITS =',I6,I4,2I6,I4,2I2/)
 1030 FORMAT(' ZL =',2F8.3,' :: FC = ',A4,' =',1P,2D20.12,',  GC = ',A4,' =',2D20.12)
 1040 FORMAT(24X,' FC''=       ',1P,2D20.12,',  GC'' =',6X ,2D20.12)    
 1050 FORMAT(24X,' SIG=   ',   2F20.12)                                 
!-----------------------------------------------------------------------
      WRITE(STDOUT,1000)                                                     
!-----------------------------------------------------------------------
  10  READ(STDIN,*,END=40) X,ETA,ZLMIN,NL,MODE,KFN,WHIT,NOTE                
      IF (NL.LE.0) GO TO 40                                              
      IFAIL = 1                                                         
      MD = MOD(ABS(MODE),10)                                            
      IH =     ABS(MODE)/10                                             
      KFIN = KFN                                                        
      IF (WHIT) KFN=MAX(KFN,0)                                           
      WRITE(STDOUT,1010) X,ETA,ZLMIN,NL,MODE,KFN,NOTE                        
!-----------------------------------------------------------------------
      CALL COULCC(X,ETA,ZLMIN,NL,FC,GC,FCP,GCP,SIG,MODE,KFN,IFAIL)      
!-----------------------------------------------------------------------
      WRITE(STDOUT,1020) IFAIL,RERR,NFP,NPQ,KASE                             
      IF (IFAIL.LT.0) GO TO 30                                           
      DO 20 I=1,8                                                       
        L = LDISP(I)                                                      
        IF (L.GT.NL-IFAIL) GO TO 20                                        
        ZL = ZLMIN + L - 1                                             
        IF (KFN.NE.0) SIG(L) = ZERO                                     
        IRREG = WHO(2,MAX(KFN+1,1),IH+1)                               
        REG = WHO(1,MAX(KFN+1,1),1)                                  
        IF (WHIT) THEN                                                  
          IRREG = 'WHIT'                                              
          WS = EXP(-HALF*PI*(ETA - CI*ZL)  - CI*SIG(L))               
          GC(L)  = WS * GC(L)                                         
          IF (MD.EQ.1) GCP(L) = CI*WS * GCP(L)                         
          FC(L)  = CI/WS * FC(L)                                      
          IF (MOD(MD,2).EQ.1) FCP(L)  = FCP(L) / WS                    
          REG = 'WH-F'                                               
        END IF                                                        
        WRITE(STDOUT,1030) ZL,REG,FC(L),IRREG,GC(L)                            
        IF (MD.EQ.1)WRITE(STDOUT,1040) FCP(L),GCP(L)                 
        IF (SIG(L).NE.ZERO.AND.KFIN.EQ.0) WRITE(STDOUT,1050) SIG(L)      
  20  CONTINUE                                                          
  30  CONTINUE                                                          
      GO TO 10                                                          
  40  STOP                                                              
!-----------------------------------------------------------------------
END PROGRAM CCTEST