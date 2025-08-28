ACDPCOULCC.  COULCC, A CONTINUED-FRACTION ALGORITHM FOR COULOMB         ACDP0000
1   FUNCTIONS OF COMPLEX ORDER WITH COMPLEX ARGUMENTS.  I.J. THOMPSON,  ACDP0000
2   A.R. BARNETT.                                                       ACDP0000
REF. IN COMP. PHYS. COMMUN. 36 (1985) 363                               ACDP0000
//TPCPC27 JOB (45800,TP),IAN,MSGLEVEL=(2,0),NOTIFY=TP,MSGCLASS=T        ACDP0001
// EXEC FVCLG,REGION.C=1000K,TIME.C=(0,15),PARM.C='OPTIMIZE(0)',        ACDP0002
//   TIME.L=(0,2),REGION.G=400K,TIME.G=(0,12)                           ACDP0003
//C.SYSIN DD *                                                          ACDP0004
      SUBROUTINE COULCC(XX,ETA1,ZLMIN,NL, FC,GC,FCP,GCP, SIG,           ACDP0005
     X                  MODE1,KFN,IFAIL)                                ACDP0006
C                                                                       ACDP0007
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACDP0008
C                                                                      CACDP0009
C  COMPLEX COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD           CACDP0010
C                                                                      CACDP0011
C  A. R. Barnett           Manchester  March   1981                    CACDP0012
C  modified I.J. Thompson  Daresbury, Sept. 1983 for Complex Functions CACDP0013
C                                                                      CACDP0014
C  original program  RCWFN       in    CPC  8 (1974) 377-395           CACDP0015
C                 +  RCWFF       in    CPC 11 (1976) 141-142           CACDP0016
C                 +  COULFG      in    CPC 27 (1982) 147-166           CACDP0017
C  description of real algorithm in    CPC 21 (1981) 297-314           CACDP0018
C  description of complex algorithm    JCP XX (1985) YYY-ZZZ           CACDP0019
C  this version written up       in    CPC XX (1985) YYY-ZZZ           CACDP0020
C                                                                      CACDP0021
C  COULCC returns F,G,G',G',SIG for complex XX, ETA1, and ZLMIN,       CACDP0022
C   for NL integer-spaced lambda values ZLMIN to ZLMIN+NL-1 inclusive, CACDP0023
C   thus giving  complex-energy solutions to the Coulomb Schrodinger   CACDP0024
C   equation,to the Klein-Gordon equation and to suitable forms of     CACDP0025
C   the Dirac equation ,also spherical & cylindrical Bessel equations  CACDP0026
C                                                                      CACDP0027
C  if /MODE1/= 1  get F,G,F',G'   for integer-spaced lambda values     CACDP0028
C            = 2      F,G      unused arrays must be dimensioned in    CACDP0029
C            = 3      F,  F'          call to at least length (1)      CACDP0030
C            = 4      F                                                CACDP0031
C            = 11 get F,H+,F',H+' ) if KFN=0, H+ = G + i.F        )    CACDP0032
C            = 12     F,H+        )       >0, H+ = J + i.Y = H(1) ) in CACDP0033
C            = 21 get F,H-,F',H-' ) if KFN=0, H- = G - i.F        ) GC CACDP0034
C            = 22     F,H-        )       >0, H- = J - i.Y = H(2) )    CACDP0035
C                                                                      CACDP0036
C     if MODE1<0 then the values returned are scaled by an exponential CACDP0037
C                factor (dependent only on XX) to bring nearer unity   CACDP0038
C                the functions for large /XX/, small ETA & /ZL/ < /XX/ CACDP0039
C        Define SCALE = (  0        if MODE1 > 0                       CACDP0040
C                       (  IMAG(XX) if MODE1 < 0  &  KFN < 3           CACDP0041
C                       (  REAL(XX) if MODE1 < 0  &  KFN = 3           CACDP0042
C        then FC = EXP(-ABS(SCALE)) * ( F, j, J, or I)                 CACDP0043
C         and GC = EXP(-ABS(SCALE)) * ( G, y, or Y )                   CACDP0044
C               or EXP(SCALE)       * ( H+, H(1), or K)                CACDP0045
C               or EXP(-SCALE)      * ( H- or H(2) )                   CACDP0046
C                                                                      CACDP0047
C  if  KFN  =  0,-1  complex Coulomb functions are returned   F & G    CACDP0048
C           =  1   spherical Bessel      "      "     "       j & y    CACDP0049
C           =  2 cylindrical Bessel      "      "     "       J & Y    CACDP0050
C           =  3 modified cyl. Bessel    "      "     "       I & K    CACDP0051
C                                                                      CACDP0052
C          and where Coulomb phase shifts put in SIG if KFN=0 (not -1) CACDP0053
C                                                                      CACDP0054
C  The use of MODE and KFN is independent                              CACDP0055
C    (except that for KFN=3,  H(1) & H(2) are not given)               CACDP0056
C                                                                      CACDP0057
C  With negative orders lambda, COULCC can still be used but with      CACDP0058
C    reduced accuracy as CF1 becomes unstable. The user is thus        CACDP0059
C    strongly advised to use reflection formulae based on              CACDP0060
C    H+-(ZL,,) = H+-(-ZL-1,,) * exp +-i(sig(ZL)-sig(-ZL-1)-(ZL+1/2)pi) CACDP0061
C                                                                      CACDP0062
C  Precision:  results to within 2-3 decimals of 'machine accuracy',   CACDP0063
C               but if CF1A fails because X too small or ETA too large CACDP0064
C               the F solution  is less accurate if it decreases with  CACDP0065
C               decreasing lambda (e.g. for lambda.LE.-1 & ETA.NE.0)   CACDP0066
C              RERR in COMMON/STEED/ traces the main roundoff errors.  CACDP0067
C                                                                      CACDP0068
C   COULCC is coded for real*8 on IBM or equivalent  ACCUR >= 10**-14  CACDP0069
C          with a section of doubled REAL*16 for less roundoff errors. CACDP0070
C          (If no doubled precision available, increase JMAX to eg 100)CACDP0071
C   Use IMPLICIT COMPLEX*32 & REAL*16 on VS compiler ACCUR >= 10**-32  CACDP0072
C   For single precision CDC (48 bits) reassign REAL*8=REAL etc.       CACDP0073
C                                                                      CACDP0074
C   IFAIL  on input   = 0 : no printing of error messages              CACDP0075
C                    ne 0 : print error messages on file 6             CACDP0076
C   IFAIL  in output = -2 : argument out of range                      CACDP0077
C                    = -1 : one of the continued fractions failed,     CACDP0078
C                           or arithmetic check before final recursion CACDP0079
C                    =  0 : All Calculations satisfactory              CACDP0080
C                    ge 0 : results available for orders up to & at    CACDP0081
C                             position NL-IFAIL in the output arrays.  CACDP0082
C                    = -3 : values at ZLMIN not found as over/underflowCACDP0083
C                    = -4 : roundoff errors make results meaningless   CACDP0084
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACDP0085
C                                                                      CACDP0086
C     Machine dependent constants :                                    CACDP0087
C                                                                      CACDP0088
C     ACCUR    target bound on relative error (except near 0 crossings)CACDP0089
C               (ACCUR should be at least 100 * ACC8)                  CACDP0090
C     ACC8     smallest number with 1+ACC8 .ne.1 in REAL*8  arithmetic CACDP0091
C     ACC16    smallest number with 1+ACC16.ne.1 in REAL*16 arithmetic CACDP0092
C     FPMAX    magnitude of largest floating point number * ACC8       CACDP0093
C     FPMIN    magnitude of smallest floating point number / ACC8      CACDP0094
C     FPLMAX   LOG(FPMAX)                                              CACDP0095
C     FPLMIN   LOG(FPMIN)                                              CACDP0096
C                                                                      CACDP0097
C     ROUTINES CALLED :       LOGAM/CLOGAM/CDIGAM,                     CACDP0098
C                             F20, CF1A, RCF, CF1C, CF2, F11, CF1R     CACDP0099
C     Intrinsic functions :   MIN, MAX, SQRT, REAL, IMAG, ABS, LOG, EXP,ACDP0100
C      (Generic names)        NINT, MOD, ATAN, ATAN2, COS, SIN, DCMPLX, ACDP0101
C                             SIGN, CONJG, INT, TANH                   CACDP0102
C     Note: Statement fntn.   NINTC = integer nearest to a complex no. CACDP0103
C                                                                      CACDP0104
C     Parameters determining region of calculations :                  CACDP0105
C                                                                      CACDP0106
C        R20      estimate of (2F0 iterations)/(CF2 iterations)        CACDP0107
C        ASYM     minimum X/(ETA**2+L) for CF1A to converge easily     CACDP0108
C        XNEAR    minimum ABS(X) for CF2 to converge accurately        CACDP0109
C        LIMIT    maximum no. iterations for CF1, CF2, and 1F1 series  CACDP0110
C        JMAX     size of work arrays for Pade accelerations           CACDP0111
C        NDROP    number of successive decrements to define instabilityCACDP0112
C                                                                      CACDP0113
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACDP0114
C                                                                       ACDP0115
      IMPLICIT COMPLEX*16 (A-H,O-Z)                                     ACDP0116
      PARAMETER(JMAX=50)                                                ACDP0117
      DIMENSION FC(NL),GC(NL),FCP(NL),GCP(NL),SIG(NL),XRCF(JMAX,4)      ACDP0118
      LOGICAL PR,ETANE0,IFCP,RLEL,DONEM,UNSTAB,ZLNEG,AXIAL,NOCF2,NPINT  ACDP0119
      REAL*8 ERR,RERR,ABSC,ACCUR,ACCT,ACC8,ACCH,ACC16,ACCB, XNEAR,CF1R, ACDP0120
     X       ZERO,ONE,TWO,HALF,HPI,TLOG,FPMAX,FPMIN,FPLMIN,FPLMAX,      ACDP0121
     X       PACCQ,EPS,OFF,SCALE,SF,SFSH,TA,RK,OMEGA,R20,ASYM,ABSX      ACDP0122
C                                                                       ACDP0123
      COMMON       /STEED/ RERR,NFP,N11,NPQ(2),N20,KAS(2)               ACDP0124
C***  common blocks are for information & storage only.                 ACDP0125
C     (they are not essential to working of the code)                   ACDP0126
      COMMON /RCFCM1/ PK,EK,CLGAA,CLGAB,CLGBB,DSIG,TPK1,W,RL,FCL1,Q,GAM,ACDP0127
     X                HCL,HPL,FCM,HCL1,ALPHA,BETA,PL                    ACDP0128
      EQUIVALENCE            (PK,XRCF(1,1))                             ACDP0129
C                                                                       ACDP0130
      DATA ZERO,ONE,TWO,LIMIT /0.0D+0, 1.0D+0, 2.0D+0, 20000 /,         ACDP0131
     X     HALF, CI / 0.5D+0, (0D+0, 1D+0) /,                           ACDP0132
     X     FPMAX,FPMIN,FPLMAX,FPLMIN / 1D+60,1D-60 ,140D+0, -140D+0 /,  ACDP0133
     X     R20,ASYM,XNEAR,NDROP / 3., 3., .5, 5 /,                      ACDP0134
     X     ACCUR, ACC8, ACC16 / 1D-14, 2D-16, 3D-33 /                   ACDP0135
      NINTC(W) = NINT(REAL(REAL(W)))                                    ACDP0136
      ABSC(W) = ABS(REAL(W)) + ABS(IMAG(W))                             ACDP0137
      NPINT(W,ACCB) = ABSC(NINTC(W)-W).LT.ACCB .AND. REAL(W).LT.HALF    ACDP0138
C                                                                       ACDP0139
      MODE = MOD(ABS(MODE1),10)                                         ACDP0140
      IFCP = MOD(MODE,2).EQ.1                                           ACDP0141
      PR = IFAIL.NE.0                                                   ACDP0142
      IFAIL = -2                                                        ACDP0143
      N11   = 0                                                         ACDP0144
      NFP   = 0                                                         ACDP0145
      KAS(1)   = 0                                                      ACDP0146
      KAS(2)   = 0                                                      ACDP0147
      NPQ(1)   = 0                                                      ACDP0148
      NPQ(2)   = 0                                                      ACDP0149
      N20 = 0                                                           ACDP0150
      HPI = TWO*ATAN(ONE)                                               ACDP0151
      TLOG = LOG(TWO)                                                   ACDP0152
      ACCUR = MAX(ACCUR, 50*ACC8)                                       ACDP0153
      ACCT = ACCUR * .5                                                 ACDP0154
C                       initialise the log-gamma function :             ACDP0155
      CALL LOGAM(ACC8)                                                  ACDP0156
      ACCH  = SQRT(ACCUR)                                               ACDP0157
      ACCB  = SQRT(ACCH)                                                ACDP0158
      RERR = ACCT                                                       ACDP0159
C                                                                       ACDP0160
      CIK = ONE                                                         ACDP0161
         IF(KFN.GE.3) CIK = CI * SIGN(ONE,ACC8-IMAG(XX))                ACDP0162
      X     = XX * CIK                                                  ACDP0163
      ETA   = ETA1                                                      ACDP0164
      IF(KFN .GT. 0) ETA = ZERO                                         ACDP0165
         ETANE0  = ABSC(ETA).GT.ACC8                                    ACDP0166
         ETAI = ETA*CI                                                  ACDP0167
      DELL  = ZERO                                                      ACDP0168
      IF(KFN .GE. 2)  DELL = HALF                                       ACDP0169
      ZM1   = ZLMIN - DELL                                              ACDP0170
      SCALE = ZERO                                                      ACDP0171
      IF(MODE1.LT.0) SCALE = IMAG(X)                                    ACDP0172
C                                                                       ACDP0173
      M1 = 1                                                            ACDP0174
      L1  = M1 + NL - 1                                                 ACDP0175
      RLEL = ABS(IMAG(ETA)) + ABS(IMAG(ZM1)) .LT. ACC8                  ACDP0176
      ABSX = ABS(X)                                                     ACDP0177
      AXIAL = RLEL .AND. ABS(IMAG(X)) .LT. ACC8 * ABSX                  ACDP0178
      IF(MODE.LE.2 .AND. ABSX.LT.FPMIN) GO TO 310                       ACDP0179
      XI  = ONE/X                                                       ACDP0180
      XLOG = LOG(X)                                                     ACDP0181
C            log with cut along the negative real axis] see also OMEGA  ACDP0182
      ID = 1                                                            ACDP0183
      DONEM = .FALSE.                                                   ACDP0184
         UNSTAB = .FALSE.                                               ACDP0185
      LF = M1                                                           ACDP0186
      IFAIL = -1                                                        ACDP0187
   10    ZLM = ZM1 + LF - M1                                            ACDP0188
         ZLL = ZM1 + L1 - M1                                            ACDP0189
C                                                                       ACDP0190
C ***       ZLL  is final lambda value, or 0.5 smaller for J,Y Bessels  ACDP0191
C                                                                       ACDP0192
              Z11 = ZLL                                                 ACDP0193
              IF(ID.LT.0) Z11 = ZLM                                     ACDP0194
              P11 = CI*SIGN(ONE,ACC8-IMAG(ETA))                         ACDP0195
      LAST = L1                                                         ACDP0196
C                                                                       ACDP0197
C ***       Find phase shifts and Gamow factor at lambda = ZLL          ACDP0198
C                                                                       ACDP0199
      PK = ZLL + ONE                                                    ACDP0200
      AA = PK - ETAI                                                    ACDP0201
      AB = PK + ETAI                                                    ACDP0202
      BB = TWO*PK                                                       ACDP0203
         ZLNEG = NPINT(BB,ACCB)                                         ACDP0204
                     CLGAA = CLOGAM(AA)                                 ACDP0205
                     CLGAB = CLGAA                                      ACDP0206
         IF(ETANE0.AND..NOT.RLEL)  CLGAB = CLOGAM(AB)                   ACDP0207
         IF(ETANE0.AND.     RLEL)  CLGAB = CONJG(CLGAA)                 ACDP0208
         SIGMA = (CLGAA - CLGAB) * CI*HALF                              ACDP0209
         IF(KFN.EQ.0) SIG(L1) = SIGMA                                   ACDP0210
         IF(.NOT.ZLNEG) CLL = ZLL*TLOG- HPI*ETA - CLOGAM(BB)            ACDP0211
     X                                          + (CLGAA+CLGAB)*HALF    ACDP0212
              THETA  = X - ETA*(XLOG+TLOG) - ZLL*HPI + SIGMA            ACDP0213
C                                                                       ACDP0214
        TA = (IMAG(AA)**2+IMAG(AB)**2+ABS(REAL(AA))+ABS(REAL(AB)))*HALF ACDP0215
      IF(ID.GT.0 .AND. ABSX .LT. TA*ASYM .AND. .NOT.ZLNEG) GO TO 20     ACDP0216
C                                                                       ACDP0217
C ***         use CF1 instead of CF1A, if predicted to converge faster, ACDP0218
C                 (otherwise using CF1A as it treats negative lambda &  ACDP0219
C                  recurrence-unstable cases properly)                  ACDP0220
C                                                                       ACDP0221
           RK = SIGN(ONE, REAL(X) + ACC8)                               ACDP0222
           P =  THETA                                                   ACDP0223
           IF(RK.LT.0) P = -X + ETA*(LOG(-X)+TLOG)-ZLL*HPI-SIGMA        ACDP0224
      F = RK * CF1A(X*RK,ETA*RK,ZLL,P,ACCT,JMAX,NFP,FEST,ERR,FPMAX,XRCF,ACDP0225
     X                                      XRCF(1,3), XRCF(1,4))       ACDP0226
      FESL = LOG(FEST) + ABS(IMAG(X))                                   ACDP0227
         NFP = - NFP                                                    ACDP0228
      IF(NFP.LT.0   .OR.(UNSTAB.AND.ERR.LT.ACCB)) GO TO 40              ACDP0229
      IF(.NOT.ZLNEG .OR. UNSTAB.AND.ERR.GT.ACCB)  GO TO 20              ACDP0230
         IF(PR) WRITE(6,1060) '-L',ERR                                  ACDP0231
         IF(ERR.GT.ACCB) GO TO 280                                      ACDP0232
         GO TO 40                                                       ACDP0233
C                                                                       ACDP0234
C ***    evaluate CF1  =  f   =  F'(ZLL,ETA,X)/F(ZLL,ETA,X)             ACDP0235
C                                                                       ACDP0236
   20 IF(AXIAL) THEN                                                    ACDP0237
C                                                        REAL VERSION   ACDP0238
      F = CF1R(X,ETA,ZLL,ACC8,SF ,RK,  ETANE0,LIMIT,ERR,NFP,            ACDP0239
     X         ACCH,FPMIN,FPMAX,PR,'COULCC')                            ACDP0240
          FCL = SF                                                      ACDP0241
          TPK1= RK                                                      ACDP0242
         ELSE                                                           ACDP0243
C                                                        COMPLEX VERSIONACDP0244
      F = CF1C(X,ETA,ZLL,ACC8,FCL,TPK1,ETANE0,LIMIT,ERR,NFP,            ACDP0245
     X         ACCH,FPMIN,FPMAX,PR,'COULCC')                            ACDP0246
         ENDIF                                                          ACDP0247
      IF(ERR.GT.ONE) GO TO 390                                          ACDP0248
C                                                                       ACDP0249
C ***  Make a simple check for CF1 being badly unstable:                ACDP0250
C                                                                       ACDP0251
      IF(ID.LT.0) GO TO 30                                              ACDP0252
      UNSTAB = REAL((ONE-ETA*XI)*CI*IMAG(THETA)/F).GT.ZERO              ACDP0253
     X .AND..NOT.AXIAL .AND. ABS(IMAG(THETA)).GT.-LOG(ACC8)*.5          ACDP0254
     X .AND. ABSC(ETA)+ABSC(ZLL).LT.ABSC(X)                             ACDP0255
      IF(UNSTAB) GO TO 60                                               ACDP0256
C                                                                       ACDP0257
C *** compare accumulated phase FCL with asymptotic phase for G(k+1) :  ACDP0258
C     to determine estimate of F(ZLL) (with correct sign) to start recurACDP0259
C                                                                       ACDP0260
   30 W   =  X*X  *(HALF/TPK1 + ONE/TPK1**2) + ETA*(ETA-TWO*X)/TPK1     ACDP0261
      FESL   = (ZLL+ONE) * XLOG + CLL - W - LOG(FCL)                    ACDP0262
   40 FESL = FESL - ABS(SCALE)                                          ACDP0263
          RK   =        MAX(REAL(FESL), FPLMIN*HALF)                    ACDP0264
          FESL = DCMPLX(MIN(RK,   FPLMAX*HALF ) , IMAG(FESL))           ACDP0265
      FEST= EXP(FESL)                                                   ACDP0266
C                                                                       ACDP0267
           RERR = MAX(RERR, ERR, ACC8 * ABS(REAL(THETA)) )              ACDP0268
C                                                                       ACDP0269
      FCL = FEST                                                        ACDP0270
      FPL = FCL*F                                                       ACDP0271
      IF(IFCP) FCP(L1) = FPL                                            ACDP0272
               FC (L1) = FCL                                            ACDP0273
C                                                                       ACDP0274
C *** downward recurrence to lambda = ZLM. array GC,if present,stores RLACDP0275
C                                                                       ACDP0276
      I  = MAX(-ID, 0)                                                  ACDP0277
      ZL  = ZLL + I                                                     ACDP0278
         MONO = 0                                                       ACDP0279
        OFF = ABS(FCL)                                                  ACDP0280
         TA = ABSC(SIGMA)                                               ACDP0281
      DO 70  L  = L1-ID,LF,-ID                                          ACDP0282
         IF(ETANE0) THEN                                                ACDP0283
               IF(RLEL) THEN                                            ACDP0284
                    DSIG = ATAN2(REAL(ETA),REAL(ZL))                    ACDP0285
                    RL = SQRT(REAL(ZL)**2 + REAL(ETA)**2)               ACDP0286
                  ELSE                                                  ACDP0287
                    AA = ZL - ETAI                                      ACDP0288
                    BB = ZL + ETAI                                      ACDP0289
                    IF(ABSC(AA).LT.ACCH.OR.ABSC(BB).LT.ACCH) GOTO 50    ACDP0290
                    DSIG = (LOG(AA) - LOG(BB)) * CI*HALF                ACDP0291
                    RL = AA * EXP(CI*DSIG)                              ACDP0292
                 ENDIF                                                  ACDP0293
             IF(ABSC(SIGMA).LT.TA*HALF) THEN                            ACDP0294
C               re-calculate SIGMA because of accumulating roundoffs:   ACDP0295
                SL =(CLOGAM(ZL+I-ETAI)-CLOGAM(ZL+I+ETAI))*CI*HALF       ACDP0296
                RL = (ZL - ETAI) * EXP(CI*ID*(SIGMA - SL))              ACDP0297
                SIGMA = SL                                              ACDP0298
                TA = ZERO                                               ACDP0299
              ELSE                                                      ACDP0300
                SIGMA = SIGMA - DSIG * ID                               ACDP0301
              ENDIF                                                     ACDP0302
                TA = MAX(TA, ABSC(SIGMA))                               ACDP0303
             SL    =  ETA  + ZL*ZL*XI                                   ACDP0304
                PL = ZERO                                               ACDP0305
                IF(ABSC(ZL).GT.ACCH) PL = (SL*SL - RL*RL)/ZL            ACDP0306
             FCL1  = (FCL *SL + ID*ZL*FPL)/RL                           ACDP0307
              SF = ABS(FCL1)                                            ACDP0308
                       IF(SF.GT.FPMAX) GO TO 350                        ACDP0309
             FPL   = (FPL *SL + ID*PL*FCL)/RL                           ACDP0310
             IF(MODE .LE. 1) GCP(L+ID)= PL * ID                         ACDP0311
        ELSE                                                            ACDP0312
C                               ETA = 0, including Bessels.  NB RL==SL  ACDP0313
           RL = ZL* XI                                                  ACDP0314
           FCL1 = FCL * RL + FPL*ID                                     ACDP0315
              SF = ABS(FCL1)                                            ACDP0316
                      IF(SF.GT.FPMAX) GO TO 350                         ACDP0317
           FPL  =(FCL1* RL - FCL) * ID                                  ACDP0318
        ENDIF                                                           ACDP0319
C             IF(ABSC(FCL1).LT.ABSC(FCL)) THEN                          ACDP0320
              IF(SF.LT.OFF) THEN                                        ACDP0321
                 MONO = MONO + 1                                        ACDP0322
                ELSE                                                    ACDP0323
                 MONO = 0                                               ACDP0324
                ENDIF                                                   ACDP0325
         FCL   =  FCL1                                                  ACDP0326
           OFF = SF                                                     ACDP0327
         FC(L) =  FCL                                                   ACDP0328
         IF(IFCP) FCP(L)  = FPL                                         ACDP0329
           IF(KFN.EQ.0) SIG(L) = SIGMA                                  ACDP0330
           IF(MODE .LE. 2) GC(L+ID) = RL                                ACDP0331
      ZL = ZL - ID                                                      ACDP0332
      IF(MONO.LT.NDROP) GO TO 70                                        ACDP0333
      IF(AXIAL .OR. REAL(ZLM)*ID.GT.-NDROP.AND..NOT.ETANE0) GO TO 70    ACDP0334
         UNSTAB = .TRUE.                                                ACDP0335
C                                                                       ACDP0336
C ***    take action if cannot or should not recur below this ZL:       ACDP0337
   50    ZLM = ZL                                                       ACDP0338
         LF = L                                                         ACDP0339
            IF(ID.LT.0) GO TO 380                                       ACDP0340
         IF(.NOT.UNSTAB) LF = L + 1                                     ACDP0341
         IF(L+MONO.LT.L1-2 .OR. ID.LT.0 .OR. .NOT.UNSTAB) GO TO 80      ACDP0342
C             otherwise, all L values (for stability) should be done    ACDP0343
C                        in the reverse direction:                      ACDP0344
         GO TO 60                                                       ACDP0345
   70 CONTINUE                                                          ACDP0346
      GO TO 80                                                          ACDP0347
   60       ID = -1                                                     ACDP0348
            LF = L1                                                     ACDP0349
            L1 = M1                                                     ACDP0350
            RERR = ACCT                                                 ACDP0351
            GO TO 10                                                    ACDP0352
   80 IF(FCL .EQ. ZERO) FCL = + ACC8                                    ACDP0353
      F  = FPL/FCL                                                      ACDP0354
C                                                                       ACDP0355
C *** Check, if second time around, that the 'f' values agree]          ACDP0356
C                                                                       ACDP0357
      IF(ID.GT.0) FIRST = F                                             ACDP0358
      IF(DONEM) RERR = MAX(RERR, ABSC(F-FIRST)/ABSC(F))                 ACDP0359
      IF(DONEM) GO TO 90                                                ACDP0360
C                                                                       ACDP0361
       NOCF2 = .FALSE.                                                  ACDP0362
      THETAM  = X - ETA*(XLOG+TLOG) - ZLM*HPI + SIGMA                   ACDP0363
C                                                                       ACDP0364
C *** on left x-plane, determine OMEGA by requiring cut on -x axis      ACDP0365
C     on right x-plane, choose OMEGA (using estimate based on THETAM)   ACDP0366
C       so H(omega) is smaller and recurs upwards accurately.           ACDP0367
C     (x-plane boundary is shifted to give CF2(LH) a chance to converge)ACDP0368
C                                                                       ACDP0369
                           OMEGA = SIGN(ONE,IMAG(X)+ACC8)               ACDP0370
      IF(REAL(X).GE.XNEAR) OMEGA = SIGN(ONE,IMAG(THETAM)+ACC8)          ACDP0371
C                                                                       ACDP0372
         SFSH = EXP(OMEGA*SCALE - ABS(SCALE))                           ACDP0373
         OFF=EXP(MIN(TWO * MAX(ABS(IMAG(X)),ABS(IMAG(THETAM)),          ACDP0374
     X                         ABS(IMAG(ZLM))*3 ) , FPLMAX) )           ACDP0375
          EPS = MAX(ACC8 , ACCT * HALF / OFF)                           ACDP0376
C                                                                       ACDP0377
C ***    Try first estimated omega, then its opposite,                  ACDP0378
C        to find the H(omega) linearly independent of F                 ACDP0379
C        i.e. maximise  CF1-CF2 = 1/(F H(omega)) , to minimise H(omega) ACDP0380
C                                                                       ACDP0381
   90 DO 100 L=1,2                                                      ACDP0382
         LH = 1                                                         ACDP0383
         IF(OMEGA.LT.ZERO) LH = 2                                       ACDP0384
      PM = CI*OMEGA                                                     ACDP0385
      ETAP = ETA * PM                                                   ACDP0386
         IF(DONEM) GO TO 130                                            ACDP0387
         PQ1 = ZERO                                                     ACDP0388
         PACCQ = ONE                                                    ACDP0389
         KASE = 0                                                       ACDP0390
C                                                                       ACDP0391
C ***            Check for small X, i.e. whether to avoid CF2 :         ACDP0392
C                                                                       ACDP0393
      IF(MODE.GE.3 .AND. ABSX.LT.ONE ) GO TO 190                        ACDP0394
      IF(MODE.LT.3 .AND. (NOCF2 .OR. ABSX.LT.XNEAR .AND.                ACDP0395
     X   ABSC(ETA)*ABSX .LT. 5 .AND. ABSC(ZLM).LT.4)) THEN              ACDP0396
        KASE = 5                                                        ACDP0397
        GO TO 120                                                       ACDP0398
        ENDIF                                                           ACDP0399
C                                                                       ACDP0400
C ***  Evaluate   CF2 : PQ1 = p + i.omega.q  at lambda = ZLM            ACDP0401
C                                                                       ACDP0402
         PQ1 = CF2(X,ETA,ZLM,PM,EPS,LIMIT,ERR,NPQ(LH),ACC8,ACCH,        ACDP0403
     X             PR,ACCUR,DELL,'COULCC')                              ACDP0404
C                                                                       ACDP0405
       ERR = ERR * MAX(ONE,ABSC(PQ1)/MAX(ABSC(F-PQ1),ACC8) )            ACDP0406
       IF(ERR.LT.ACCH)       GO TO 110                                  ACDP0407
C                                                                       ACDP0408
C *** check if impossible to get F-PQ accurately because of cancellationACDP0409
               NOCF2 = REAL(X).LT.XNEAR .AND. ABS(IMAG(X)).LT.-LOG(ACC8)ACDP0410
C                original guess for OMEGA (based on THETAM) was wrong   ACDP0411
C                Use KASE 5 or 6 if necessary if Re(X) < XNEAR          ACDP0412
  100            OMEGA = - OMEGA                                        ACDP0413
                IF(UNSTAB) GO TO 360                                    ACDP0414
                IF(REAL(X).LT.-XNEAR .AND. PR) WRITE(6,1060) '-X',ERR   ACDP0415
  110     RERR = MAX(RERR,ERR)                                          ACDP0416
C                                                                       ACDP0417
C ***  establish case of calculation required for irregular solution    ACDP0418
C                                                                       ACDP0419
  120 IF(KASE.GE.5) GO TO 130                                           ACDP0420
      IF(REAL(X) .GT. XNEAR) THEN                                       ACDP0421
C          estimate errors if KASE 2 or 3 were to be used:              ACDP0422
         PACCQ = EPS * OFF * ABSC(PQ1) / MAX(ABS(IMAG(PQ1)),ACC8)       ACDP0423
        ENDIF                                                           ACDP0424
      IF(PACCQ .LT. ACCUR) THEN                                         ACDP0425
          KASE = 2                                                      ACDP0426
          IF(AXIAL) KASE = 3                                            ACDP0427
      ELSE                                                              ACDP0428
          KASE = 1                                                      ACDP0429
          IF(NPQ(1) * R20 .LT. JMAX)     KASE = 4                       ACDP0430
C             i.e. change to kase=4 if the 2F0 predicted to converge    ACDP0431
      ENDIF                                                             ACDP0432
  130 GO TO (190,140,150,170,190,190),  ABS(KASE)                       ACDP0433
  140    IF(.NOT.DONEM)                                                 ACDP0434
C                                                                       ACDP0435
C ***  Evaluate   CF2 : PQ2 = p - i.omega.q  at lambda = ZLM   (Kase 2) ACDP0436
C                                                                       ACDP0437
     X  PQ2 = CF2(X,ETA,ZLM,-PM,EPS,LIMIT,ERR,NPQ(3-LH),ACC8,ACCH,      ACDP0438
     X             PR,ACCUR,DELL,'COULCC')                              ACDP0439
C                                                                       ACDP0440
        P     = (PQ2 + PQ1) * HALF                                      ACDP0441
        Q     = (PQ2 - PQ1) * HALF*PM                                   ACDP0442
      GO TO 160                                                         ACDP0443
  150   P     = REAL(PQ1)                                               ACDP0444
        Q     = IMAG(PQ1)                                               ACDP0445
C                                                                       ACDP0446
C ***   With Kase = 3 on the real axes, P and Q are real & PQ2 = PQ1*   ACDP0447
C                                                                       ACDP0448
        PQ2 = CONJG(PQ1)                                                ACDP0449
C                                                                       ACDP0450
C *** solve for FCM = F at lambda = ZLM,then find norm factor W=FCM/FCL ACDP0451
C                                                                       ACDP0452
  160 W   = (PQ1 - F) * (PQ2 - F)                                       ACDP0453
         SF = EXP(-ABS(SCALE))                                          ACDP0454
      FCM = SQRT(Q / W) * SF                                            ACDP0455
C                  any SQRT given here is corrected by                  ACDP0456
C                  using sign for FCM nearest to phase of FCL           ACDP0457
      IF(REAL(FCM/FCL).LT.ZERO) FCM  = - FCM                            ACDP0458
      GAM = (F - P)/Q                                                   ACDP0459
         TA = ABSC(GAM + PM)                                            ACDP0460
         PACCQ= EPS * MAX(TA,ONE/TA)                                    ACDP0461
      HCL = FCM * (GAM + PM) * (SFSH/(SF*SF))                           ACDP0462
C                                                                       ACDP0463
      IF(PACCQ.GT.ACCUR .AND. KASE.GT.0) THEN                           ACDP0464
C                                    Consider a KASE = 1 Calculation    ACDP0465
          F11V= F11(X,ETA,Z11,P11,ACCT,LIMIT,0,ERR,N11,FPMAX,ACC8,ACC16)ACDP0466
          IF(ERR.LT.PACCQ) GO TO 200                                    ACDP0467
          ENDIF                                                         ACDP0468
      RERR=MAX(RERR,PACCQ)                                              ACDP0469
      GO TO 230                                                         ACDP0470
C                                                                       ACDP0471
C *** Arrive here if KASE = 4                                           ACDP0472
C     to evaluate the exponentially decreasing H(LH) directly.          ACDP0473
C                                                                       ACDP0474
  170  IF(DONEM) GO TO 180                                              ACDP0475
      AA = ETAP - ZLM                                                   ACDP0476
      BB = ETAP + ZLM + ONE                                             ACDP0477
      F20V = F20(AA,BB,-HALF*PM*XI, ACCT,JMAX,ERR,FPMAX,N20,XRCF)       ACDP0478
        IF(N20.LE.0) GO TO 190                                          ACDP0479
        RERR = MAX(RERR,ERR)                                            ACDP0480
         HCL = FPMIN                                                    ACDP0481
         IF(ABS(REAL(PM*THETAM)+OMEGA*SCALE).GT.FPLMAX) GO TO 330       ACDP0482
  180 HCL = F20V * EXP(PM * THETAM + OMEGA*SCALE)                       ACDP0483
      FCM = SFSH / ((F - PQ1) * HCL )                                   ACDP0484
      GO TO 230                                                         ACDP0485
C                                                                       ACDP0486
C *** Arrive here if KASE=1   (or if 2F0 tried mistakenly & failed)     ACDP0487
C                                                                       ACDP0488
C           for small values of X, calculate F(X,SL) directly from 1F1  ACDP0489
C               using REAL*16 arithmetic if possible.                   ACDP0490
C           where Z11 = ZLL if ID>0, or = ZLM if ID<0                   ACDP0491
C                                                                       ACDP0492
  190 F11V = F11(X,ETA,Z11,P11,ACCT,LIMIT,0,ERR,N11,FPMAX,ACC8,ACC16)   ACDP0493
C                                                                       ACDP0494
  200       IF(N11.LT.0) THEN                                           ACDP0495
C                               F11 failed from BB = negative integer   ACDP0496
               WRITE(6,1060) '-L',ONE                                   ACDP0497
               GO TO 390                                                ACDP0498
               ENDIF                                                    ACDP0499
            IF(ERR.GT.PACCQ .AND. PACCQ.LT.ACCB) THEN                   ACDP0500
C                               Consider a KASE 2 or 3 calculation :    ACDP0501
                KASE = -2                                               ACDP0502
                IF(AXIAL) KASE = -3                                     ACDP0503
                GO TO 130                                               ACDP0504
                ENDIF                                                   ACDP0505
         RERR = MAX(RERR, ERR)                                          ACDP0506
         IF(ERR.GT.FPMAX) GO TO 370                                     ACDP0507
         IF(ID.LT.0) CLL = Z11*TLOG- HPI*ETA - CLOGAM(BB)               ACDP0508
     X                       + CLOGAM(Z11 + ONE + P11*ETA) - P11*SIGMA  ACDP0509
      EK   = (Z11+ONE)*XLOG - P11*X + CLL  - ABS(SCALE)                 ACDP0510
      IF(ID.GT.0) EK = EK - FESL + LOG(FCL)                             ACDP0511
         IF(REAL(EK).GT.FPLMAX) GO TO 350                               ACDP0512
         IF(REAL(EK).LT.FPLMIN) GO TO 340                               ACDP0513
      FCM = F11V * EXP(EK)                                              ACDP0514
C                                                                       ACDP0515
      IF(KASE.GE.5) THEN                                                ACDP0516
        IF(ABSC(ZLM+ZLM-NINTC(ZLM+ZLM)).LT.ACCH) KASE = 6               ACDP0517
C                                                                       ACDP0518
C ***  For abs(X) < XNEAR, then CF2 may not converge accurately, so     ACDP0519
C ***      use an expansion for irregular soln from origin :            ACDP0520
C                                                                       ACDP0521
         SL = ZLM                                                       ACDP0522
            ZLNEG = REAL(ZLM) .LT. -ONE + ACCB                          ACDP0523
         IF(KASE.EQ.5 .OR. ZLNEG) SL = - ZLM - ONE                      ACDP0524
         PK = SL + ONE                                                  ACDP0525
            AA = PK - ETAP                                              ACDP0526
            AB = PK + ETAP                                              ACDP0527
            BB = TWO*PK                                                 ACDP0528
                     CLGAA = CLOGAM(AA)                                 ACDP0529
                     CLGAB = CLGAA                                      ACDP0530
         IF(ETANE0)  CLGAB = CLOGAM(AB)                                 ACDP0531
                     CLGBB = CLOGAM(BB)                                 ACDP0532
           IF(KASE.EQ.6 .AND. .NOT.ZLNEG) THEN                          ACDP0533
              IF(NPINT(AA,ACCUR)) CLGAA = CLGAB - TWO*PM*SIGMA          ACDP0534
              IF(NPINT(AB,ACCUR)) CLGAB = CLGAA + TWO*PM*SIGMA          ACDP0535
             ENDIF                                                      ACDP0536
          CLL = SL*TLOG- HPI*ETA - CLGBB + (CLGAA + CLGAB) * HALF       ACDP0537
          DSIG = (CLGAA - CLGAB) * PM*HALF                              ACDP0538
             IF(KASE.EQ.6) P11 = - PM                                   ACDP0539
          EK  = PK * XLOG - P11*X + CLL  - ABS(SCALE)                   ACDP0540
                     SF = EXP(-ABS(SCALE))                              ACDP0541
                     CHI = ZERO                                         ACDP0542
       IF(.NOT.( KASE.EQ.5 .OR. ZLNEG ) ) GO TO 210                     ACDP0543
C                                                                       ACDP0544
C *** Use  G(l)  =  (cos(CHI) * F(l) - F(-l-1)) /  sin(CHI)             ACDP0545
C                                                                       ACDP0546
C      where CHI = sig(l) - sig(-l-1) - (2l+1)*pi/2                     ACDP0547
C                                                                       ACDP0548
         CHI = SIGMA - DSIG - (ZLM-SL) * HPI                            ACDP0549
         F11V=F11(X,ETA,SL,P11,ACCT,LIMIT,0,ERR,NPQ(1),FPMAX,ACC8,ACC16)ACDP0550
                    RERR = MAX(RERR,ERR)                                ACDP0551
            IF(KASE.EQ.6) GO TO 210                                     ACDP0552
         FESL = F11V * EXP( EK )                                        ACDP0553
         FCL1 = EXP(PM*CHI) * FCM                                       ACDP0554
         HCL = FCL1 - FESL                                              ACDP0555
               RERR=MAX(RERR,ACCT*MAX(ABSC(FCL1),ABSC(FESL))/ABSC(HCL)) ACDP0556
         HCL = HCL / SIN(CHI) * (SFSH/(SF*SF))                          ACDP0557
       GO TO 220                                                        ACDP0558
C                                                                       ACDP0559
C *** Use the logarithmic expansion for the irregular solution (KASE 6) ACDP0560
C        for the case that BB is integral so sin(CHI) would be zero.    ACDP0561
C                                                                       ACDP0562
  210    RL = BB - ONE                                                  ACDP0563
         N  = NINTC(RL)                                                 ACDP0564
         ZLOG = XLOG + TLOG - PM*HPI                                    ACDP0565
         CHI = CHI + PM * THETAM + OMEGA * SCALE + AB * ZLOG            ACDP0566
            AA  = ONE - AA                                              ACDP0567
         IF(NPINT(AA,ACCUR)) THEN                                       ACDP0568
            HCL = ZERO                                                  ACDP0569
         ELSE                                                           ACDP0570
               IF(ID.GT.0 .AND. .NOT.ZLNEG) F11V = FCM * EXP(-EK)       ACDP0571
            HCL = EXP(CHI - CLGBB - CLOGAM(AA)) * (-1)**(N+1)           ACDP0572
     X              * ( F11V * ZLOG +                                   ACDP0573
     X      F11(X,ETA,SL,-PM,ACCT,LIMIT,2,ERR,NPQ(2),FPMAX,ACC8,ACC16)) ACDP0574
                RERR = MAX(RERR,ERR)                                    ACDP0575
            ENDIF                                                       ACDP0576
         IF(N.GT.0) THEN                                                ACDP0577
             EK = CHI + CLOGAM(RL) - CLGAB - RL*ZLOG                    ACDP0578
             DF =F11(X,ETA,-SL-ONE,-PM,ZERO,N,0,ERR,L,FPMAX,ACC8,ACC16) ACDP0579
             HCL = HCL + EXP(EK) * DF                                   ACDP0580
            ENDIF                                                       ACDP0581
C                                                                       ACDP0582
  220    PQ1 = F - SFSH/(FCM * HCL)                                     ACDP0583
      ELSE                                                              ACDP0584
           IF(MODE.LE.2) HCL = SFSH/((F - PQ1) * FCM)                   ACDP0585
           KASE = 1                                                     ACDP0586
      ENDIF                                                             ACDP0587
C                                                                       ACDP0588
C ***  Now have absolute normalisations for Coulomb Functions           ACDP0589
C          FCM & HCL  at lambda = ZLM                                   ACDP0590
C      so determine linear transformations for Functions required :     ACDP0591
C                                                                       ACDP0592
  230 IH = ABS(MODE1) / 10                                              ACDP0593
        IF(KFN.EQ.3) IH = (3-IMAG(CIK))/2  + HALF                       ACDP0594
      P11 = ONE                                                         ACDP0595
      IF(IH.EQ.1) P11 = CI                                              ACDP0596
      IF(IH.EQ.2) P11 = -CI                                             ACDP0597
                  DF = - PM                                             ACDP0598
      IF(IH.GE.1) DF = - PM + P11                                       ACDP0599
          IF(ABSC(DF).LT.ACCH) DF = ZERO                                ACDP0600
C                                                                       ACDP0601
C *** Normalisations for spherical or cylindrical Bessel functions      ACDP0602
C                                                                       ACDP0603
                          ALPHA = ZERO                                  ACDP0604
          IF(KFN  .EQ. 1) ALPHA = XI                                    ACDP0605
          IF(KFN  .GE. 2) ALPHA = XI*HALF                               ACDP0606
                          BETA  = ONE                                   ACDP0607
          IF(KFN  .EQ. 1) BETA  = XI                                    ACDP0608
          IF(KFN  .GE. 2) BETA  = SQRT(XI/HPI)                          ACDP0609
          IF(KFN  .GE. 2 .AND. REAL(BETA).LT.ZERO) BETA  = - BETA       ACDP0610
C                                                                       ACDP0611
      AA = ONE                                                          ACDP0612
      IF(KFN.GT.0) AA = -P11 * BETA                                     ACDP0613
      IF(KFN.GE.3) THEN                                                 ACDP0614
C                        Calculate rescaling factors for I & K output   ACDP0615
         P = EXP((ZLM+DELL) * HPI * CIK)                                ACDP0616
         AA= BETA * HPI * P                                             ACDP0617
         BETA = BETA / P                                                ACDP0618
         Q = CIK * ID                                                   ACDP0619
        ENDIF                                                           ACDP0620
C                        Calculate rescaling factors for GC output      ACDP0621
      IF(IH.EQ.0) THEN                                                  ACDP0622
         TA = ABS(SCALE) + IMAG(PM)*SCALE                               ACDP0623
         RK = ZERO                                                      ACDP0624
         IF(TA.LT.FPLMAX) RK = EXP(-TA)                                 ACDP0625
       ELSE                                                             ACDP0626
         TA = ABS(SCALE) + IMAG(P11)*SCALE                              ACDP0627
C                                                                       ACDP0628
         IF(ABSC(DF).GT.ACCH .AND. TA.GT.FPLMAX) GO TO 320              ACDP0629
         IF(ABSC(DF).GT.ACCH) DF = DF * EXP(TA)                         ACDP0630
         SF = TWO * (LH-IH) * SCALE                                     ACDP0631
         RK = ZERO                                                      ACDP0632
         IF(SF.GT.FPLMAX) GO TO 320                                     ACDP0633
         IF(SF.GT.FPLMIN) RK = EXP(SF)                                  ACDP0634
      ENDIF                                                             ACDP0635
C                                                                       ACDP0636
         KAS((3-ID)/2) = KASE                                           ACDP0637
      W = FCM / FCL                                                     ACDP0638
         IF(LOG(ABSC(W))+LOG(ABSC(FC(LF))) .LT. FPLMIN) GO TO 340       ACDP0639
         IF(MODE.GE.3) GO TO 240                                        ACDP0640
            IF(ABSC(F-PQ1) .LT. ACCH*ABSC(F) .AND. PR)                  ACDP0641
     X                             WRITE(6,1020) LH,ZLM+DELL            ACDP0642
      HPL = HCL * PQ1                                                   ACDP0643
         IF(ABSC(HPL).LT.FPMIN.OR.ABSC(HCL).LT.FPMIN) GO TO 330         ACDP0644
C                                                                       ACDP0645
C *** IDward recurrence from HCL,HPL(LF) (stored GC(L) is RL if reqd)   ACDP0646
C *** renormalise FC,FCP at each lambda                                 ACDP0647
C ***    ZL   = ZLM - MIN(ID,0) here                                    ACDP0648
C                                                                       ACDP0649
  240 DO 270 L = LF,L1,ID                                               ACDP0650
                     FCL = W* FC(L)                                     ACDP0651
                      IF(ABSC(FCL).LT.FPMIN) GO TO 340                  ACDP0652
            IF(IFCP) FPL = W*FCP(L)                                     ACDP0653
                     FC(L)  = BETA * FCL                                ACDP0654
            IF(IFCP) FCP(L) = BETA * (FPL - ALPHA * FCL) * CIK          ACDP0655
                     FC(L)  = TIDY(FC(L),ACCUR)                         ACDP0656
            IF(IFCP) FCP(L) = TIDY(FCP(L),ACCUR)                        ACDP0657
       IF(MODE .GE. 3) GO TO 260                                        ACDP0658
       IF(L.EQ.LF)  GO TO 250                                           ACDP0659
                      ZL = ZL + ID                                      ACDP0660
                      ZID= ZL * ID                                      ACDP0661
                      RL = GC(L)                                        ACDP0662
         IF(ETANE0)   THEN                                              ACDP0663
                      SL = ETA + ZL*ZL*XI                               ACDP0664
            IF(MODE.EQ.1) THEN                                          ACDP0665
              PL = GCP(L)                                               ACDP0666
            ELSE                                                        ACDP0667
              PL = ZERO                                                 ACDP0668
              IF(ABSC(ZL).GT.ACCH) PL = (SL*SL - RL*RL)/ZID             ACDP0669
            ENDIF                                                       ACDP0670
           HCL1     = (SL*HCL - ZID*HPL) / RL                           ACDP0671
           HPL      = (SL*HPL - PL *HCL) / RL                           ACDP0672
         ELSE                                                           ACDP0673
           HCL1 = RL * HCL - HPL * ID                                   ACDP0674
           HPL  = (HCL - RL * HCL1) * ID                                ACDP0675
         ENDIF                                                          ACDP0676
         HCL      = HCL1                                                ACDP0677
         IF(ABSC(HCL).GT.FPMAX) GO TO 320                               ACDP0678
  250    GC(L) = AA * (RK * HCL + DF * FCL)                             ACDP0679
      IF(MODE.EQ.1) GCP(L) = (AA *(RK*HPL +DF*FPL) - ALPHA * GC(L)) *CIKACDP0680
         GC(L) = TIDY(GC(L),ACCUR)                                      ACDP0681
      IF(MODE.EQ.1) GCP(L) = TIDY(GCP(L),ACCUR)                         ACDP0682
         IF(KFN.GE.3) AA = AA * Q                                       ACDP0683
  260    IF(KFN.GE.3) BETA = - BETA * Q                                 ACDP0684
  270  LAST = MIN(LAST,(L1 - L)*ID)                                     ACDP0685
C                                                                       ACDP0686
C *** Come here after all soft errors to determine how many L values ok ACDP0687
C                                                                       ACDP0688
  280  IF(ID.GT.0 .OR.  LAST.EQ.0) IFAIL = LAST                         ACDP0689
       IF(ID.LT.0 .AND. LAST.NE.0) IFAIL = -3                           ACDP0690
C                                                                       ACDP0691
C *** Come here after ALL errors for this L range (ZLM,ZLL)             ACDP0692
C                                                                       ACDP0693
  290 IF(ID.GT.0 .AND. LF.NE.M1) GO TO 300                              ACDP0694
         IF(IFAIL.LT.0) RETURN                                          ACDP0695
         IF(RERR.GT.ACCB) WRITE(6,1070) RERR                            ACDP0696
         IF(RERR.GT.0.1) IFAIL = -4                                     ACDP0697
         RETURN                                                         ACDP0698
C                                                                       ACDP0699
C *** so on first block, 'F' started decreasing monotonically,          ACDP0700
C                        or hit bound states for low ZL.                ACDP0701
C     thus redo M1 to LF-1 in reverse direction                         ACDP0702
C      i.e. do CF1A at ZLMIN & CF2 at ZLM (midway between ZLMIN & ZLMAX)ACDP0703
C                                                                       ACDP0704
  300 ID = -1                                                           ACDP0705
      IF(.NOT.UNSTAB) LF = LF - 1                                       ACDP0706
      DONEM = UNSTAB                                                    ACDP0707
      LF = MIN(LF,L1)                                                   ACDP0708
      L1 = M1                                                           ACDP0709
      GO TO 10                                                          ACDP0710
C                                                                       ACDP0711
C ***    error messages                                                 ACDP0712
C                                                                       ACDP0713
  310 IF(PR) WRITE (6,1000) XX                                          ACDP0714
 1000 FORMAT(/' COULCC: CANNOT CALCULATE IRREGULAR SOLUTIONS FOR X =',  ACDP0715
     X 1P,2D10.2,', AS ABS(X) IS TOO SMALL'/)                           ACDP0716
      RETURN                                                            ACDP0717
  320 IF(PR) WRITE(6,1010) ZL+DELL,'IR',HCL,'MORE',FPMAX                ACDP0718
 1010 FORMAT(' COULCC: AT ZL =',2F8.3,' ',A2,'REGULAR SOLUTION (',1P,   ACDP0719
     X 2E10.1,') WILL BE ',A4,' THAN',E10.1)                            ACDP0720
      GO TO 280                                                         ACDP0721
  330 IF(PR) WRITE(6,1010) ZL+DELL,'IR',HCL,'LESS',FPMIN                ACDP0722
      GO TO 280                                                         ACDP0723
  340 IF(PR) WRITE(6,1010) ZL+DELL,'  ',FCL,'LESS',FPMIN                ACDP0724
      GO TO 280                                                         ACDP0725
  350 IF(PR) WRITE(6,1010) ZL+DELL,'  ',FCL,'MORE',FPMAX                ACDP0726
      GO TO 280                                                         ACDP0727
 1020 FORMAT('0COULCC WARNING: LINEAR INDEPENDENCE BETWEEN ''F'' AND ''HACDP0728
     X(',I1,')'' IS LOST AT ZL =',2F7.2,' (EG. COULOMB EIGENSTATE, OR CFACDP0729
     X1 UNSTABLE)'/)                                                    ACDP0730
  360 IF(PR) WRITE(6,1030) ZLL+DELL                                     ACDP0731
 1030 FORMAT(' COULCC: (ETA&L)/X TOO LARGE FOR CF1A, AND CF1 UNSTABLE ATACDP0732
     X L =',2F8.2)                                                      ACDP0733
      GO TO 280                                                         ACDP0734
  370 IF(PR) WRITE(6,1040) Z11,I                                        ACDP0735
 1040 FORMAT(' COULCC: OVERFLOW IN 1F1 SERIES AT ZL =',2F8.3,' AT TERM',ACDP0736
     X I5)                                                              ACDP0737
      GO TO 390                                                         ACDP0738
  380 IF(PR) WRITE(6,1050) ZLMIN,ZLM,ZLM+ONE,ZLMIN+NL-ONE               ACDP0739
 1050 FORMAT(' COULCC: BOTH BOUND-STATE POLES AND F-INSTABILITIES OCCUR'ACDP0740
     X ,', OR MULTIPLE INSTABILITIES PRESENT.'                          ACDP0741
     X,/,' TRY CALLING TWICE,  FIRST FOR ZL FROM',2F8.3,' TO',2F8.3,    ACDP0742
     X ' (INCL.)',/,20X,     'SECOND FOR ZL FROM',2F8.3,' TO',2F8.3)    ACDP0743
C     GO TO 390                                                         ACDP0744
  390 IFAIL = -1                                                        ACDP0745
      GO TO 290                                                         ACDP0746
 1060 FORMAT('0COULCC WARNING: AS ''',A2,''' REFLECTION RULES NOT USED, ACDP0747
     #ERRORS CAN BE UP TO',1P,D12.2/)                                   ACDP0748
 1070 FORMAT('0COULCC WARNING: OVERALL ROUNDOFF ERROR APPROX.',1P,E11.1)ACDP0749
      END                                                               ACDP0750
      FUNCTION CF1C(X,ETA,ZL,EPS,FCL,TPK1,ETANE0,LIMIT,ERR,NFP,         ACDP0751
     X              ACCH,FPMIN,FPMAX,PR,CALLER)                         ACDP0752
      IMPLICIT COMPLEX*16(A-H,O-Z)                                      ACDP0753
      LOGICAL PR,ETANE0                                                 ACDP0754
      REAL*8 ONE,TWO,EPS,ERR,ACCH,FPMIN,FPMAX,ABSC,SMALL,RK,PX          ACDP0755
      CHARACTER*6 CALLER                                                ACDP0756
      DATA ONE,TWO / 1D+0, 2D+0 /                                       ACDP0757
      ABSC(W) = ABS(REAL(W)) + ABS(IMAG(W))                             ACDP0758
C                                                                       ACDP0759
C                                                                       ACDP0760
C ***    Evaluate CF1  =  F   =  F'(ZL,ETA,X)/F(ZL,ETA,X)               ACDP0761
C                                                                       ACDP0762
C        using complex arithmetic                                       ACDP0763
C                                                                       ACDP0764
      FCL = ONE                                                         ACDP0765
      XI = ONE/X                                                        ACDP0766
      PK  = ZL + ONE                                                    ACDP0767
      PX  = PK  + LIMIT                                                 ACDP0768
   10 EK  = ETA / PK                                                    ACDP0769
        RK2 =          ONE + EK*EK                                      ACDP0770
      F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI                           ACDP0771
      PK1 =  PK + ONE                                                   ACDP0772
         TPK1 = PK + PK1                                                ACDP0773
      TK  = TPK1*(XI + EK/PK1)                                          ACDP0774
      IF(ETANE0) THEN                                                   ACDP0775
C ***   test ensures b1 .ne. zero for negative ETA etc.; fixup is exact.ACDP0776
             IF(ABSC(TK) .GT. ACCH)  GO TO 20                           ACDP0777
             FCL  = RK2/(ONE + (ETA/PK1)**2)                            ACDP0778
             SL   = TPK1*XI * (TPK1+TWO)*XI                             ACDP0779
             PK   =  TWO + PK                                           ACDP0780
             GO TO 10                                                   ACDP0781
         ENDIF                                                          ACDP0782
   20 D   =  ONE/TK                                                     ACDP0783
      DF  = -FCL*RK2*D                                                  ACDP0784
            IF(REAL(PK).GT.REAL(ZL)+TWO) FCL = - RK2 * SL               ACDP0785
            FCL = FCL * D * TPK1 * XI                                   ACDP0786
      F   =  F  + DF                                                    ACDP0787
C                                                                       ACDP0788
C ***   begin CF1 loop on PK = k = lambda + 1                           ACDP0789
C                                                                       ACDP0790
      RK    = ONE                                                       ACDP0791
      SMALL    = SQRT(FPMIN)                                            ACDP0792
   30 PK    = PK1                                                       ACDP0793
        PK1 = PK1 + ONE                                                 ACDP0794
         TPK1 = PK + PK1                                                ACDP0795
         IF(ETANE0) THEN                                                ACDP0796
           EK  = ETA / PK                                               ACDP0797
           RK2 =          ONE + EK*EK                                   ACDP0798
          ENDIF                                                         ACDP0799
        TK  = TPK1*(XI + EK/PK1)                                        ACDP0800
        D   =  TK - D*RK2                                               ACDP0801
              IF(ABSC(D) .GT. ACCH)             GO TO 40                ACDP0802
              IF(PR) WRITE (6,1000) CALLER,D,DF,ACCH,PK,EK,ETA,X        ACDP0803
              RK= RK +   ONE                                            ACDP0804
              IF( RK .GT. TWO )                  GO TO 50               ACDP0805
   40 D     = ONE/D                                                     ACDP0806
            FCL = FCL * D * TPK1*XI                                     ACDP0807
            IF(ABSC(FCL).LT.SMALL) FCL = FCL / SMALL                    ACDP0808
            IF(ABSC(FCL).GT.FPMAX) FCL = FCL / FPMAX                    ACDP0809
        DF  = DF*(D*TK - ONE)                                           ACDP0810
        F   = F  + DF                                                   ACDP0811
              IF( REAL(PK) .GT. PX ) GO TO 50                           ACDP0812
      IF(ABSC(DF) .GE. ABSC(F)*EPS)             GO TO 30                ACDP0813
                NFP = PK - ZL - 1                                       ACDP0814
                  ERR = EPS * SQRT(REAL(NFP))                           ACDP0815
      CF1C = F                                                          ACDP0816
      RETURN                                                            ACDP0817
 1000 FORMAT(/' ',A6,': CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,ETA,X = ', ACDP0818
     X    /1X,1P,13D9.2/)                                               ACDP0819
   50 IF(PR) WRITE (6,1010) CALLER,LIMIT,ABS(X)                         ACDP0820
 1010 FORMAT(' ',A6,': CF1 HAS FAILED TO CONVERGE AFTER ',I10  ,' ITERATACDP0821
     XIONS AS ABS(X) =',F15.0)                                          ACDP0822
      ERR = TWO                                                         ACDP0823
      RETURN                                                            ACDP0824
      END                                                               ACDP0825
      FUNCTION CF2(X,ETA,ZL,PM,EPS,LIMIT,ERR,NPQ,ACC8,ACCH,             ACDP0826
     X             PR,ACCUR,DELL,CALLER)                                ACDP0827
      IMPLICIT COMPLEX*16(A-H,O-Z)                                      ACDP0828
      LOGICAL PR                                                        ACDP0829
      REAL*8 EPS,ERR,ACC8,ACCH,ACCUR,TA,RK,                             ACDP0830
     X       ABSC,ZERO,HALF,ONE,TWO                                     ACDP0831
      CHARACTER*6 CALLER                                                ACDP0832
      DATA ZERO,HALF,ONE,TWO / 0D+0, .5D+0, 1D+0, 2D+0 /                ACDP0833
      ABSC(W) = ABS(REAL(W)) + ABS(IMAG(W))                             ACDP0834
C                                                                       ACDP0835
C                                    (omega)        (omega)             ACDP0836
C *** Evaluate  CF2  = p + PM.q  =  H   (ETA,X)' / H   (ETA,X)          ACDP0837
C                                    ZL             ZL                  ACDP0838
C     where PM = omega.i                                                ACDP0839
C                                                                       ACDP0840
      TA = TWO*LIMIT                                                    ACDP0841
      E2MM1 = ETA*ETA + ZL*ZL + ZL                                      ACDP0842
      ETAP = ETA * PM                                                   ACDP0843
      XI = ONE/X                                                        ACDP0844
      WI = TWO*ETAP                                                     ACDP0845
      RK = ZERO                                                         ACDP0846
      PQ = (ONE - ETA*XI) * PM                                          ACDP0847
      AA = -E2MM1 + ETAP                                                ACDP0848
      BB = TWO*(X - ETA + PM)                                           ACDP0849
         RL = XI * PM                                                   ACDP0850
      IF(ABSC(BB).LT.ACCH) THEN                                         ACDP0851
         RL = RL * AA / (AA + RK + WI)                                  ACDP0852
         PQ = PQ + RL * (BB + TWO*PM)                                   ACDP0853
            AA = AA + TWO*(RK+ONE+WI)                                   ACDP0854
            BB = BB + (TWO+TWO)*PM                                      ACDP0855
            RK = RK + (TWO+TWO)                                         ACDP0856
         ENDIF                                                          ACDP0857
      DD = ONE/BB                                                       ACDP0858
      DL = AA*DD* RL                                                    ACDP0859
   10 PQ    = PQ + DL                                                   ACDP0860
         RK = RK + TWO                                                  ACDP0861
         AA = AA + RK + WI                                              ACDP0862
         BB = BB + TWO*PM                                               ACDP0863
         DD = ONE/(AA*DD + BB)                                          ACDP0864
         DL = DL*(BB*DD - ONE)                                          ACDP0865
            ERR = ABSC(DL)/ABSC(PQ)                                     ACDP0866
         IF(ERR.GE.MAX(EPS,ACC8*RK*HALF) .AND. RK.LE.TA) GO TO 10       ACDP0867
C                                                                       ACDP0868
         NPQ   = RK/TWO                                                 ACDP0869
         PQ    = PQ + DL                                                ACDP0870
           IF(PR.AND.NPQ.GE.LIMIT-1 .AND. ERR.GT.ACCUR)                 ACDP0871
     X             WRITE(6,1000) CALLER,INT(IMAG(PM)),NPQ,ERR,ZL+DELL   ACDP0872
 1000 FORMAT(' ',A6,': CF2(',I2,') NOT CONVERGED FULLY IN ',I7,         ACDP0873
     X' ITERATIONS, SO ERROR IN IRREGULAR SOLUTION =',1P,D11.2,' AT ZL  ACDP0874
     X=', 0P,2F8.3)                                                     ACDP0875
      CF2 = PQ                                                          ACDP0876
      RETURN                                                            ACDP0877
      END                                                               ACDP0878
      FUNCTION F11(X,ETA,ZL,P,EPS,LIMIT,KIND,ERR,NITS,FPMAX,ACC8,ACC16) ACDP0879
      IMPLICIT REAL*8(A-H,O-Z)                                          ACDP0880
      COMPLEX*16 X,ETA,ZL,P,AA,BB,Z,F11,CDIGAM,CI                       ACDP0881
       COMPLEX*16 DD,G,F,AI,BI,T                                        ACDP0882
      LOGICAL ZLLIN                                                     ACDP0883
      REAL*16 AR,BR,GR,GI,DR,DI,TR,TI,UR,UI,FI,FI1,DEN                  ACDP0884
      DATA ZERO,ONE,TWO / 0D+0, 1D+0, 2D+0 /, CI / (0D+0, 1D+0) /       ACDP0885
      ABSC(AA) = ABS(REAL(AA)) + ABS(IMAG(AA))                          ACDP0886
      NINTC(AA) = NINT(REAL(REAL(AA)))                                  ACDP0887
C                                                                       ACDP0888
C *** evaluate the HYPERGEOMETRIC FUNCTION 1F1                          ACDP0889
C                                        i                              ACDP0890
C            F (AA;BB; Z) = SUM  (AA)   Z / ( (BB)  i] )                ACDP0891
C           1 1              i       i            i                     ACDP0892
C                                                                       ACDP0893
C     to accuracy EPS with at most LIMIT terms.                         ACDP0894
C  If KIND = 0 : using extended precision but real arithmetic only,     ACDP0895
C            1 : using normal precision in complex arithmetic,          ACDP0896
C   or       2 : using normal complex arithmetic, but with CDIGAM factorACDP0897
C                                                                       ACDP0898
C  where                                                                ACDP0899
         AA = ZL+ONE - ETA*P                                            ACDP0900
         BB = TWO*(ZL+ONE)                                              ACDP0901
C  and                                                                  ACDP0902
         Z  = TWO*P*X                                                   ACDP0903
C                                                                       ACDP0904
         ZLLIN = REAL(BB).LE.ZERO .AND. ABS(BB-NINTC(BB)).LT.ACC8**0.25 ACDP0905
             IF(.NOT.ZLLIN.OR.REAL(BB)+LIMIT.LT.1.5) GO TO 10           ACDP0906
                NITS = -1                                               ACDP0907
                RETURN                                                  ACDP0908
   10 IF(LIMIT.LE.0) THEN                                               ACDP0909
         F11 = ZERO                                                     ACDP0910
         ERR = ZERO                                                     ACDP0911
         NITS= 1                                                        ACDP0912
         RETURN                                                         ACDP0913
         ENDIF                                                          ACDP0914
      TA = ONE                                                          ACDP0915
      RK = ONE                                                          ACDP0916
      IF(KIND.LE.0.AND.ABSC(Z)*ABSC(AA).GT.ABSC(BB) * 1.0) THEN         ACDP0917
         DR = ONE                                                       ACDP0918
         DI = ZERO                                                      ACDP0919
         GR = ONE                                                       ACDP0920
         GI = ZERO                                                      ACDP0921
         AR = REAL(AA)                                                  ACDP0922
         BR = REAL(BB)                                                  ACDP0923
         FI = ZERO                                                      ACDP0924
      DO 20 I=2,LIMIT                                                   ACDP0925
         FI1 = FI + ONE                                                 ACDP0926
         TR = BR * FI1                                                  ACDP0927
         TI = IMAG(BB) * FI1                                            ACDP0928
         DEN= ONE / (TR*TR + TI*TI)                                     ACDP0929
         UR = (AR*TR + IMAG(AA)*TI) * DEN                               ACDP0930
         UI = (IMAG(AA)*TR - AR*TI) * DEN                               ACDP0931
         TR = UR*GR - UI*GI                                             ACDP0932
         TI = UR*GI + UI*GR                                             ACDP0933
         GR = REAL(Z) * TR - IMAG(Z)*TI                                 ACDP0934
         GI = REAL(Z) * TI + IMAG(Z)*TR                                 ACDP0935
         DR = DR + GR                                                   ACDP0936
         DI = DI + GI                                                   ACDP0937
            ERR = ABS(GR) + ABS(GI)                                     ACDP0938
               IF(ERR.GT.FPMAX) GO TO 60                                ACDP0939
            RK  = ABS(DR) + ABS(DI)                                     ACDP0940
            TA = MAX(TA,RK)                                             ACDP0941
         IF(ERR.LT.RK*EPS .OR. I.GE.4.AND.ERR.LT.ACC16) GO TO 30        ACDP0942
         FI = FI1                                                       ACDP0943
         AR = AR + ONE                                                  ACDP0944
   20    BR = BR + ONE                                                  ACDP0945
C                                                                       ACDP0946
   30    F11 = DR + CI * DI                                             ACDP0947
         ERR = ACC16 * TA / RK                                          ACDP0948
C                                                                       ACDP0949
      ELSE                                                              ACDP0950
C* ---------------------------------- alternative code                  ACDP0951
C*    If REAL*16 arithmetic is not available, (or already using it]),   ACDP0952
C*    then use KIND > 0                                                 ACDP0953
         G = ONE                                                        ACDP0954
          F = ONE                                                       ACDP0955
          IF(KIND.GE.2) F = CDIGAM(AA) - CDIGAM(BB) - CDIGAM(G)         ACDP0956
         DD = F                                                         ACDP0957
         DO 40 I=2,LIMIT                                                ACDP0958
            AI = AA + (I-2)                                             ACDP0959
            BI = BB + (I-2)                                             ACDP0960
            R  = I-ONE                                                  ACDP0961
         G = G * Z * AI / (BI * R)                                      ACDP0962
         IF(KIND.GE.2)                                                  ACDP0963
C                              multiply by (psi(a+r)-psi(b+r)-psi(1+r)) ACDP0964
     X        F = F + ONE/AI - ONE/BI - ONE/R                           ACDP0965
         T  = G * F                                                     ACDP0966
         DD = DD + T                                                    ACDP0967
            ERR = ABSC(T)                                               ACDP0968
               IF(ERR.GT.FPMAX) GO TO 60                                ACDP0969
            RK = ABSC(DD)                                               ACDP0970
         TA = MAX(TA,RK)                                                ACDP0971
         IF(ERR.LT.RK*EPS.OR.ERR.LT.ACC8.AND.I.GE.4) GO TO 50           ACDP0972
   40    CONTINUE                                                       ACDP0973
                                                                        ACDP0974
   50    ERR = ACC8 * TA / RK                                           ACDP0975
         F11 = DD                                                       ACDP0976
C* ------------------------------------------- end of alternative code  ACDP0977
      ENDIF                                                             ACDP0978
   60    NITS = I                                                       ACDP0979
      RETURN                                                            ACDP0980
      END                                                               ACDP0981
      FUNCTION CF1R(X,ETA,ZL,EPS,FCL,TPK1,ETANE0,LIMIT,ERR,NFP,         ACDP0982
     X              ACCH,FPMIN,FPMAX,PR,CALLER)                         ACDP0983
      IMPLICIT REAL*8(A-H,O-Z)                                          ACDP0984
      LOGICAL PR,ETANE0                                                 ACDP0985
      CHARACTER*6 CALLER                                                ACDP0986
      DATA ONE,TWO / 1D+0, 2D+0 /                                       ACDP0987
C                                                                       ACDP0988
C                                                                       ACDP0989
C ***    Evaluate CF1  =  F   =  F'(ZL,ETA,X)/F(ZL,ETA,X)               ACDP0990
C                                                                       ACDP0991
C        using real arithmetic                                          ACDP0992
C                                                                       ACDP0993
      FCL = ONE                                                         ACDP0994
      XI = ONE/X                                                        ACDP0995
      PK  = ZL + ONE                                                    ACDP0996
      PX  = PK  + LIMIT                                                 ACDP0997
   10 EK  = ETA / PK                                                    ACDP0998
        RK2 =          ONE + EK*EK                                      ACDP0999
      F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI                           ACDP1000
      PK1 =  PK + ONE                                                   ACDP1001
         TPK1 = PK + PK1                                                ACDP1002
      TK  = TPK1*(XI + EK/PK1)                                          ACDP1003
      IF(ETANE0) THEN                                                   ACDP1004
C ***   test ensures b1 .ne. zero for negative ETA etc.; fixup is exact.ACDP1005
             IF(ABS(TK) .GT. ACCH)  GO TO 20                            ACDP1006
             FCL  = RK2/(ONE + (ETA/PK1)**2)                            ACDP1007
             SL   = TPK1*XI * (TPK1+TWO)*XI                             ACDP1008
             PK   =  TWO + PK                                           ACDP1009
             GO TO 10                                                   ACDP1010
         ENDIF                                                          ACDP1011
   20 D   =  ONE/TK                                                     ACDP1012
      DF  = -FCL*RK2*D                                                  ACDP1013
            IF(PK.GT.ZL+TWO) FCL = - RK2 * SL                           ACDP1014
            FCL = FCL * D * TPK1 * XI                                   ACDP1015
      F   =  F  + DF                                                    ACDP1016
C                                                                       ACDP1017
C ***   begin CF1 loop on PK = k = lambda + 1                           ACDP1018
C                                                                       ACDP1019
      RK    = ONE                                                       ACDP1020
      SMALL    = SQRT(FPMIN)                                            ACDP1021
   30 PK    = PK1                                                       ACDP1022
        PK1 = PK1 + ONE                                                 ACDP1023
         TPK1 = PK + PK1                                                ACDP1024
         IF(ETANE0) THEN                                                ACDP1025
           EK  = ETA / PK                                               ACDP1026
           RK2 =          ONE + EK*EK                                   ACDP1027
          ENDIF                                                         ACDP1028
        TK  = TPK1*(XI + EK/PK1)                                        ACDP1029
        D   =  TK - D*RK2                                               ACDP1030
              IF(ABS(D) .GT. ACCH)             GO TO 40                 ACDP1031
              IF(PR) WRITE (6,1000) CALLER,D,DF,ACCH,PK,EK,ETA,X        ACDP1032
              RK= RK +   ONE                                            ACDP1033
              IF( RK .GT. TWO )                  GO TO 50               ACDP1034
   40 D     = ONE/D                                                     ACDP1035
            FCL = FCL * D * TPK1*XI                                     ACDP1036
            IF(ABS(FCL).LT.SMALL) FCL = FCL / SMALL                     ACDP1037
            IF(ABS(FCL).GT.FPMAX) FCL = FCL / FPMAX                     ACDP1038
        DF  = DF*(D*TK - ONE)                                           ACDP1039
        F   = F  + DF                                                   ACDP1040
              IF( PK .GT. PX ) GO TO 50                                 ACDP1041
      IF(ABS(DF) .GE. ABS(F)*EPS)             GO TO 30                  ACDP1042
                NFP = PK - ZL - 1                                       ACDP1043
                  ERR = EPS * SQRT(REAL(NFP))                           ACDP1044
      CF1R = F                                                          ACDP1045
      RETURN                                                            ACDP1046
 1000 FORMAT(/' ',A6,': CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,ETA,X = ', ACDP1047
     X    /1X,1P,7D9.2/)                                                ACDP1048
   50 IF(PR) WRITE (6,1010) CALLER,LIMIT,ABS(X)                         ACDP1049
 1010 FORMAT(' ',A6,': CF1 HAS FAILED TO CONVERGE AFTER ',I10  ,' ITERATACDP1050
     XIONS AS ABS(X) =',F15.0)                                          ACDP1051
      ERR = TWO                                                         ACDP1052
      RETURN                                                            ACDP1053
      END                                                               ACDP1054
      FUNCTION F20(AA,BB,Z,EPS,JMAX,RE,FPMAX,N,X)                       ACDP1055
C                                                                       ACDP1056
C     evaluate the HYPERGEOMETRIC FUNCTION 2F0                          ACDP1057
C                                             i                         ACDP1058
C            F (AA,BB;;Z) = SUM  (AA)  (BB)  Z / i]                     ACDP1059
C           2 0              i       i     i                            ACDP1060
C                                                                       ACDP1061
C     to accuracy EPS with at most JMAX terms.                          ACDP1062
C                                                                       ACDP1063
C     if the terms start diverging,                                     ACDP1064
C     the corresponding continued fraction is found by RCF              ACDP1065
C     & evaluated progressively by Steed's method to obtain convergence.ACDP1066
C                                                                       ACDP1067
C      useful number also input:  FPMAX = near-largest f.p. number      ACDP1068
C                                                                       ACDP1069
      IMPLICIT COMPLEX*16(A-H,O-Z)                                      ACDP1070
      DIMENSION X(JMAX,4)                                               ACDP1071
      LOGICAL FINITE                                                    ACDP1072
      REAL*8 EP,EPS,AT,ATL,ABSC,RE,FPMAX                                ACDP1073
      DATA ONE,ZERO / (1D+0,0D+0), (0D+0,0D+0) /                        ACDP1074
      ABSC(W) = ABS(REAL(W)) + ABS(IMAG(W))                             ACDP1075
      NINTC(W) = NINT(REAL(REAL(W)))                                    ACDP1076
C                                                                       ACDP1077
      RE = 0.0                                                          ACDP1078
      X(1,1) = ONE                                                      ACDP1079
      SUM = X(1,1)                                                      ACDP1080
      ATL = ABSC(X(1,1))                                                ACDP1081
         F    = SUM                                                     ACDP1082
         D = ONE                                                        ACDP1083
         DF   = SUM                                                     ACDP1084
      J = 0                                                             ACDP1085
      EP = EPS * JMAX *10.                                              ACDP1086
      MA = - NINTC(AA)                                                  ACDP1087
      MB = - NINTC(BB)                                                  ACDP1088
      FINITE = ABS(ABS(REAL(AA))-MA).LT.EP .AND. ABS(IMAG(AA)).LT.EP    ACDP1089
     X    .OR. ABS(ABS(REAL(BB))-MB).LT.EP .AND. ABS(IMAG(BB)).LT.EP    ACDP1090
      IMAX = JMAX                                                       ACDP1091
      IF(FINITE.AND.MA.GE.0) IMAX = MIN(MA+1,IMAX)                      ACDP1092
      IF(FINITE.AND.MB.GE.0) IMAX = MIN(MB+1,IMAX)                      ACDP1093
      DO 10 I=2,IMAX                                                    ACDP1094
      X(I,1) = X(I-1,1) * Z * (AA+I-2) * (BB+I-2) / (I-1)               ACDP1095
         IF(ABSC(X(I,1)).GT.FPMAX) GO TO 40                             ACDP1096
      AT = ABSC(X(I,1))                                                 ACDP1097
         IF(J.EQ.0) THEN                                                ACDP1098
                 SUM = SUM + X(I,1)                                     ACDP1099
                 IF(AT .LT. ABSC(SUM)*EPS) GO TO 20                     ACDP1100
               ENDIF                                                    ACDP1101
      IF(FINITE) GO TO 10                                               ACDP1102
      IF(J.GT.0 .OR. AT.GT.ATL .OR. I.GE.JMAX-2) J = J + 1              ACDP1103
         IF(J.EQ.0) GO TO 10                                            ACDP1104
         CALL RCF(X(1,1),X(1,2),J,I,X(1,3),EPS)                         ACDP1105
              IF(I.LT.0) GO TO 40                                       ACDP1106
            DO 50 K=MAX(J,2),I                                          ACDP1107
            D = ONE/(D*X(K,2) + ONE)                                    ACDP1108
            DF = DF*(D - ONE)                                           ACDP1109
            F = F + DF                                                  ACDP1110
            IF(ABSC(DF) .LT. ABSC(F)*EPS) GO TO 30                      ACDP1111
            IF(DF.EQ.ZERO.AND.F.EQ.ZERO.AND.I.GE.4) GO TO 30            ACDP1112
   50       CONTINUE                                                    ACDP1113
         J = I                                                          ACDP1114
   10 ATL = AT                                                          ACDP1115
      IF(.NOT.FINITE) I = -JMAX                                         ACDP1116
   20 N = I                                                             ACDP1117
       F20 = SUM                                                        ACDP1118
       IF(.NOT.FINITE) RE  = AT / ABSC(SUM)                             ACDP1119
       RETURN                                                           ACDP1120
   30 F20 = F                                                           ACDP1121
      RE = ABSC(DF) / ABSC(F)                                           ACDP1122
      N = K                                                             ACDP1123
      RETURN                                                            ACDP1124
   40 I = 0                                                             ACDP1125
      GO TO 20                                                          ACDP1126
      END                                                               ACDP1127
      FUNCTION CF1A(RHO,ETA,XL,PSI,EPS,NMAX,NUSED,FCL,RE,FPMAX,XX,G,C)  ACDP1128
C                                                                       ACDP1129
C     evaluate the ASYMPTOTIC EXPANSION for the                         ACDP1130
C            LOGARITHMIC DERIVATIVE OF THE REGULAR SOLUTION             ACDP1131
C                                                                       ACDP1132
C ***        CF1A  =  f   =  F'(XL,ETA,RHO)/F(XL,ETA,RHO)               ACDP1133
C                                                                       ACDP1134
C      that is valid for REAL(RHO)>0, and best for RHO >> ETA**2, XL,   ACDP1135
C      and is derived from the 2F0 expansions for H+ and H-             ACDP1136
C      e.g. by Froeberg (Rev. Mod. Physics Vol 27, p399 , 1955)         ACDP1137
C      Some lines of this subprogram are for convenience copied from    ACDP1138
C           Takemasa, Tamura & Wolter CPC 17 (1979) 351.                ACDP1139
C                                                                       ACDP1140
C     Evaluate to accuracy EPS with at most NMAX terms.                 ACDP1141
C                                                                       ACDP1142
C     If the terms start diverging,                                     ACDP1143
C     the corresponding continued fraction is found by RCF              ACDP1144
C     & evaluated progressively by Steed's method to obtain convergence.ACDP1145
C                                                                       ACDP1146
C      useful number also input:  FPMAX = near-largest f.p. number      ACDP1147
C                                                                       ACDP1148
      IMPLICIT COMPLEX*16(A-H,O-Z)                                      ACDP1149
      DIMENSION XX(2,NMAX),G(NMAX),C(NMAX)                              ACDP1150
      REAL*8 RE,EPS,T1,T2,T3,ZERO,ONE,TWO,AT,ATL,ABSC,FPMAX             ACDP1151
      DATA ZERO,ONE,TWO,CI / 0D+0, 1D+0, 2D+0, (0D+0,1D+0) /            ACDP1152
      ABSC(W) = ABS(REAL(W)) + ABS(IMAG(W))                             ACDP1153
C                                                                       ACDP1154
      HPI = TWO*ATAN(ONE)                                               ACDP1155
      T1 = SIN(REAL(PSI))                                               ACDP1156
      T2 = COS(REAL(PSI))                                               ACDP1157
      ATL= TANH(IMAG(PSI))                                              ACDP1158
C             GIVE COS(PSI)/COSH(IM(PSI)), WHICH ALWAYS HAS CORRECT SIGNACDP1159
          COSL = DCMPLX( T2 , -T1 * ATL )                               ACDP1160
      TANL = DCMPLX(T1,T2*ATL) / COSL                                   ACDP1161
      RE = ZERO                                                         ACDP1162
      XLL1= XL*(XL+ONE)                                                 ACDP1163
      ETASQ = ETA*ETA                                                   ACDP1164
      SL1=ONE                                                           ACDP1165
      SL=SL1                                                            ACDP1166
      SC1=ZERO                                                          ACDP1167
      SC=SC1                                                            ACDP1168
      TL1=SC                                                            ACDP1169
      TL=TL1                                                            ACDP1170
      TC1=ONE-ETA/RHO                                                   ACDP1171
      TC=TC1                                                            ACDP1172
      FCL  = TL + SL*TANL                                               ACDP1173
      G(1) = (TC + SC*TANL) / FCL                                       ACDP1174
      GLAST = G(1)                                                      ACDP1175
      ATL = ABSC(GLAST)                                                 ACDP1176
         F    = GLAST                                                   ACDP1177
         D = ONE                                                        ACDP1178
         DF   = GLAST                                                   ACDP1179
      J = 0                                                             ACDP1180
      DO 10 N=2,NMAX                                                    ACDP1181
      T1=N-1                                                            ACDP1182
      T2=TWO*T1-ONE                                                     ACDP1183
      T3=T1*(T1-ONE)                                                    ACDP1184
      DENOM=TWO*RHO*T1                                                  ACDP1185
      C1=(ETA*T2)/DENOM                                                 ACDP1186
      C2=(ETASQ+XLL1-T3)/DENOM                                          ACDP1187
      SL2=C1*SL1-C2*TL1                                                 ACDP1188
      TL2=C1*TL1+C2*SL1                                                 ACDP1189
      SC2=C1*SC1-C2*TC1-SL2/RHO                                         ACDP1190
      TC2=C1*TC1+C2*SC1-TL2/RHO                                         ACDP1191
      SL=SL+SL2                                                         ACDP1192
      TL=TL+TL2                                                         ACDP1193
      SC=SC+SC2                                                         ACDP1194
      TC=TC+TC2                                                         ACDP1195
      SL1=SL2                                                           ACDP1196
      TL1=TL2                                                           ACDP1197
      SC1=SC2                                                           ACDP1198
      TC1=TC2                                                           ACDP1199
      FCL  =  TL + SL*TANL                                              ACDP1200
         IF(ABSC(FCL).GT.FPMAX .OR. ABSC(FCL).LT.1./FPMAX) GO TO 40     ACDP1201
      GSUM = (TC + SC*TANL) / FCL                                       ACDP1202
      G(N) = GSUM - GLAST                                               ACDP1203
      GLAST = GSUM                                                      ACDP1204
         AT = ABSC(G(N))                                                ACDP1205
         IF(AT.LT.ABSC(GSUM)*EPS) GO TO 20                              ACDP1206
      IF(J.GT.0 .OR. AT.GT.ATL .OR. N.GE.NMAX-2) J = J + 1              ACDP1207
         IF(J.EQ.0) GO TO 10                                            ACDP1208
            CALL RCF(G,C,J,N,XX,EPS)                                    ACDP1209
              IF(N.LT.0) GO TO 40                                       ACDP1210
            DO 60 K=MAX(J,2),N                                          ACDP1211
               D = ONE/(D*C(K) + ONE)                                   ACDP1212
               DF = DF*(D - ONE)                                        ACDP1213
               F = F + DF                                               ACDP1214
         IF(ABSC(DF) .LT. ABSC(F)*EPS) GO TO 30                         ACDP1215
         IF(DF.EQ.ZERO.AND.F.EQ.ZERO.AND.N.GE.4) GO TO 30               ACDP1216
   60         CONTINUE                                                  ACDP1217
         J = N                                                          ACDP1218
   10    ATL = AT                                                       ACDP1219
      K = -NMAX                                                         ACDP1220
      GO TO 30                                                          ACDP1221
   20 FCL = FCL * COSL                                                  ACDP1222
         CF1A = GSUM                                                    ACDP1223
         RE = AT / ABSC(GSUM)                                           ACDP1224
         NUSED = N                                                      ACDP1225
         RETURN                                                         ACDP1226
   30 CF1A = F                                                          ACDP1227
      FCL = FCL * COSL                                                  ACDP1228
         RE = ABSC(DF) / ABSC(F)                                        ACDP1229
         NUSED = K                                                      ACDP1230
      RETURN                                                            ACDP1231
   40 CF1A = G(1)                                                       ACDP1232
      FCL = 1.0                                                         ACDP1233
      RE = 1.0                                                          ACDP1234
      NUSED = 0                                                         ACDP1235
      RETURN                                                            ACDP1236
      END                                                               ACDP1237
      SUBROUTINE RCF(A,B,IBEG,INUM,XX,EPS)                              ACDP1238
C                                                                       ACDP1239
C*******************************************************************    ACDP1240
C                                                                       ACDP1241
C  RCF converts polynomial A to the corresponding continued             ACDP1242
C         fraction, in 'normal'  form with coefficients B               ACDP1243
C         by the 'P algorithmn' of Patry & Gupta                        ACDP1244
C                                                                       ACDP1245
C   A(z) = A1/z + A2/z**3 + A3/z**5 + ... + An/z**(2n-1)                ACDP1246
C                                                                       ACDP1247
C   B(z) = B1/z+ B2/z+ B3/z+ .../(z+ Bn/z)                              ACDP1248
C                                                                       ACDP1249
C  data:                                                                ACDP1250
C   A     vector A(k), k=1,INUM         input                           ACDP1251
C   B     vector B(k), k=IBEG,INUM      output                          ACDP1252
C   IBEG  order of first coef. calc.    input                           ACDP1253
C   INUM  order of A, even or odd       input                           ACDP1254
C   XX    auxiliary vector of length .ge. length of vector B            ACDP1255
C         caller provides space for A,B,XX                              ACDP1256
C     Note that neither of the first two terms A(1) A(2) should be zero ACDP1257
C             & the user can start the calculation with any value of    ACDP1258
C                IBEG provided the c.f. coefs have been already         ACDP1259
C                calculated up to INUM = IBEG-1                         ACDP1260
C             & the method breaks down as soon as the absolute value    ACDP1261
C                of a c.f. coef. is less than EPS.    At the time of theACDP1262
C                break up XX(1) has been replaced by 1E-50, and INUM hasACDP1263
C                been replaced by minus times the number of this coef.  ACDP1264
C   algorithm: J.Patry & S.Gupta,                                       ACDP1265
C              EIR-bericht nr. 247,                                     ACDP1266
C              Eidg. Institut fur Reaktorforschung Wuerenlingen         ACDP1267
C              Wueringlingen, Schweiz.                                  ACDP1268
C              November 1973                                            ACDP1269
C   see also:  Haenggi,Roesel & Trautmann,                              ACDP1270
C              Jnl. Computational Physics, vol 137, pp242-258 (1980)    ACDP1271
C   note:      restart procedure modified by I.J.Thompson               ACDP1272
C                                                                       ACDP1273
C*******************************************************************    ACDP1274
C                                                                       ACDP1275
      IMPLICIT COMPLEX*16(A-H,O-Z)                                      ACDP1276
      DIMENSION A(100),B(100),XX(2,100)                                 ACDP1277
      LOGICAL EVEN                                                      ACDP1278
      REAL*8 EPS                                                        ACDP1279
      COMMON /RCFCM2/ X1,M2M1,MP12,EVEN,M                               ACDP1280
C     ibn = ibeg + inum - 1                                             ACDP1281
      IBN = INUM                                                        ACDP1282
C                             B(IBN) is last value set on this call     ACDP1283
      IF(IBEG.GT.4 .AND. M .NE. IBEG-1) GO TO 90                        ACDP1284
C                             B(M) is last value set in previous call   ACDP1285
      IF(IBEG.GT.4) GO TO 50                                            ACDP1286
      IF(IBEG.EQ.4) GO TO 20                                            ACDP1287
      B(1) = A(1)                                                       ACDP1288
      IF(IBN.GE.2) B(2) = - A(2)/A(1)                                   ACDP1289
      IF(IBN.LT.3) GO TO 10                                             ACDP1290
      X0 = A(3) / A(2)                                                  ACDP1291
      XX(2,1) = B(2)                                                    ACDP1292
      XX(1,1) = - X0                                                    ACDP1293
      XX(1,2) = 0.                                                      ACDP1294
      B(3) = -X0 - B(2)                                                 ACDP1295
      X0 = -B(3) * A(2)                                                 ACDP1296
      M = 3                                                             ACDP1297
      MP12 = 2                                                          ACDP1298
      EVEN = .TRUE.                                                     ACDP1299
      IF(IBN.GT.3) GO TO 20                                             ACDP1300
   10 RETURN                                                            ACDP1301
   20 IF(ABS(B(3)) .LT. EPS*ABS(X0)) GOTO 80                            ACDP1302
      M = 4                                                             ACDP1303
   30 X1 = A(M)                                                         ACDP1304
      M2M1 = MP12                                                       ACDP1305
      MP12 = M2M1 + 1                                                   ACDP1306
      IF(EVEN) MP12 = M2M1                                              ACDP1307
      DO 40 K=2,MP12                                                    ACDP1308
   40 X1 = X1 + A(M-K+1) * XX(1,K-1)                                    ACDP1309
      B(M) = - X1/X0                                                    ACDP1310
      IF(M.GE.IBN) RETURN                                               ACDP1311
   50 IF(ABS(B(M)).LT.EPS*ABS(X0)) GO TO 80                             ACDP1312
      K = M2M1                                                          ACDP1313
   60 XX(2,K) = XX(1,K) + B(M) * XX(2,K-1)                              ACDP1314
      K = K-1                                                           ACDP1315
      IF(K.GT.1) GO TO 60                                               ACDP1316
      XX(2,1) = XX(1,1) + B(M)                                          ACDP1317
      DO 70 K=1,M2M1                                                    ACDP1318
      X0 = XX(2,K)                                                      ACDP1319
      XX(2,K) = XX(1,K)                                                 ACDP1320
   70 XX(1,K) = X0                                                      ACDP1321
      X0 = X1                                                           ACDP1322
      XX(1,M2M1+1) = 0.                                                 ACDP1323
      M = M+1                                                           ACDP1324
      EVEN = .NOT.EVEN                                                  ACDP1325
      GO TO 30                                                          ACDP1326
   80 INUM = -M                                                         ACDP1327
C     XX(1,1) = 1.E-50                                                  ACDP1328
C     PRINT 1000,M                                                      ACDP1329
C1000 FORMAT('0RCF: ZERO CF COEFFICIENT AT POSITION ',I4/)              ACDP1330
      RETURN                                                            ACDP1331
   90 PRINT 1000,M,IBEG-1                                               ACDP1332
 1000 FORMAT('0RCF: LAST CALL SET M =',I4,', BUT RESTART REQUIRES',I4)  ACDP1333
      STOP                                                              ACDP1334
      END                                                               ACDP1335
      FUNCTION CLOGAM(Z)                                                ACDP1336
C                                                                       ACDP1337
C     this routine computes the logarithm of the gamma function gamma(z)ACDP1338
C     for any complex argument 'Z' to any accuracy preset by CALL LOGAM ACDP1339
C                                                                       ACDP1340
      IMPLICIT REAL*8(A-H,O-Z)                                          ACDP1341
      COMPLEX*16 Z,U,V,H,R,CLOGAM,CDIGAM,SER                            ACDP1342
      DIMENSION B(15),BN(15),BD(15)                                     ACDP1343
C                                                                       ACDP1344
      DATA LERR /6/, NX0 /6/, NB /15/,                                  ACDP1345
     X  ZERO,ONE,TWO,FOUR,HALF,QUART /0D+0,1D+0,2D+0,4D+0,.5D+0,.25D+0/ ACDP1346
      DATA BN(1),BD(1)    / +1D+0,   6D+0 /,                            ACDP1347
     X     BN(2),BD(2)    / -1D+0,  30D+0 /,                            ACDP1348
     X     BN(3),BD(3)    / +1D+0,  42D+0 /,                            ACDP1349
     X     BN(4),BD(4)    / -1D+0,  30D+0 /,                            ACDP1350
     X     BN(5),BD(5)    / +5D+0,  66D+0 /,                            ACDP1351
     X     BN(6),BD(6)    /          -691D+0,  2730D+0/,                ACDP1352
     X     BN(7),BD(7)    /          +  7D+0,     6D+0/,                ACDP1353
     X     BN(8),BD(8)    /         -3617D+0,   510D+0/,                ACDP1354
     X     BN(9),BD(9)    /         43867D+0,   798D+0/,                ACDP1355
     X     BN(10),BD(10)  /       -174611D+0,   330D+0/,                ACDP1356
     X     BN(11),BD(11)  /        854513D+0,   138D+0/,                ACDP1357
     X     BN(12),BD(12)  /    -236364091D+0,  2730D+0/,                ACDP1358
     X     BN(13),BD(13)  /     + 8553103D+0,     6D+0/,                ACDP1359
     X     BN(14),BD(14)  /  -23749461029D+0,   870D+0/,                ACDP1360
     X     BN(15),BD(15)  / 8615841276005D+0, 14322D+0/                 ACDP1361
      DATA FPLMIN / -140D+0 /                                           ACDP1362
C                                                                       ACDP1363
      X=REAL(Z)                                                         ACDP1364
      T=IMAG(Z)                                                         ACDP1365
      MX = INT(REAL(ACCUR*100 - X))                                     ACDP1366
      IF(ABS(ABS(X)-MX) + ABS(T).LT.ACCUR*50) GO TO 60                  ACDP1367
      F=ABS(T)                                                          ACDP1368
      V=DCMPLX(X,F)                                                     ACDP1369
      IF(X .LT. ZERO) V=ONE-V                                           ACDP1370
      H=ZERO                                                            ACDP1371
      C=REAL(V)                                                         ACDP1372
      N=NX0-INT(C)                                                      ACDP1373
      IF(N .LT. 0) GO TO 30                                             ACDP1374
      H=V                                                               ACDP1375
      D=IMAG(V)                                                         ACDP1376
      A=ATAN2(D,C)                                                      ACDP1377
      IF(N .EQ. 0) GO TO 20                                             ACDP1378
      DO 10 I = 1,N                                                     ACDP1379
      C=C+ONE                                                           ACDP1380
      V=DCMPLX(C,D)                                                     ACDP1381
      H=H*V                                                             ACDP1382
   10 A=A+ATAN2(D,C)                                                    ACDP1383
   20 H=DCMPLX(HALF*LOG(REAL(H)**2+IMAG(H)**2),A)                       ACDP1384
      V=V+ONE                                                           ACDP1385
   30 R=ONE/V**2                                                        ACDP1386
      SER = B(NT)                                                       ACDP1387
      DO 40 J=2,NT                                                      ACDP1388
        K = NT+1 - J                                                    ACDP1389
   40 SER = B(K) + R*SER                                                ACDP1390
      CLOGAM = HL2P+(V-HALF)*LOG(V)-V + SER/V - H                       ACDP1391
      IF(X .GE. ZERO) GO TO 50                                          ACDP1392
C                                                                       ACDP1393
      A= INT(X)-ONE                                                     ACDP1394
      C=PI*(X-A)                                                        ACDP1395
      D=PI*F                                                            ACDP1396
C     E=EXP(-TWO*D)                                                     ACDP1397
        E = ZERO                                                        ACDP1398
        F = -TWO*D                                                      ACDP1399
        IF(F.GT.FPLMIN) E = EXP(F)                                      ACDP1400
      F=SIN(C)                                                          ACDP1401
      E= D + HALF*LOG(E*F**2+QUART*(ONE-E)**2)                          ACDP1402
      F=ATAN2(COS(C)*TANH(D),F)-A*PI                                    ACDP1403
      CLOGAM=ALPI-DCMPLX(E,F)-CLOGAM                                    ACDP1404
C                                                                       ACDP1405
   50 IF(SIGN(ONE,T) .LT. -HALF) CLOGAM=CONJG(CLOGAM)                   ACDP1406
      RETURN                                                            ACDP1407
C                                                                       ACDP1408
   60 WRITE(LERR,1000) 'CLOGAM',X                                       ACDP1409
 1000 FORMAT(1X,A6,' ... ARGUMENT IS NON POSITIVE INTEGER = ',F20.2)    ACDP1410
      CLOGAM = ZERO                                                     ACDP1411
      RETURN                                                            ACDP1412
C                                                                       ACDP1413
      ENTRY CDIGAM(Z)                                                   ACDP1414
C                                                                       ACDP1415
C     this routine computes the logarithmic derivative of the gamma     ACDP1416
C     function  psi(Z) = digamma(Z) = d (ln gamma(Z))/dZ  for any       ACDP1417
C     complex argument Z, to any accuracy preset by CALL LOGAM(ACC)     ACDP1418
C                                                                       ACDP1419
      U=Z                                                               ACDP1420
      X=REAL(U)                                                         ACDP1421
      A=ABS(X)                                                          ACDP1422
      IF(ABS(IMAG(U)) + ABS(A + INT(X)) .LT. ACCUR) GO TO 110           ACDP1423
      IF(X .LT. ZERO) U=-U                                              ACDP1424
      V=U                                                               ACDP1425
      H=ZERO                                                            ACDP1426
      N=NX0-INT(A)                                                      ACDP1427
      IF(N .LT. 0) GO TO 90                                             ACDP1428
      H=ONE/V                                                           ACDP1429
      IF(N .EQ. 0) GO TO 80                                             ACDP1430
      DO 70 I = 1,N                                                     ACDP1431
      V=V+ONE                                                           ACDP1432
   70 H=H+ONE/V                                                         ACDP1433
   80 V=V+ONE                                                           ACDP1434
   90 R=ONE/V**2                                                        ACDP1435
      SER = B(NT) * (2*NT-1)                                            ACDP1436
      DO 100 J=2,NT                                                     ACDP1437
        K = NT+1 - J                                                    ACDP1438
  100 SER = B(K)*(2*K-1) + R*SER                                        ACDP1439
      CDIGAM = LOG(V) - HALF/V - R*SER - H                              ACDP1440
      IF(X .GE. ZERO) RETURN                                            ACDP1441
      H=PI*U                                                            ACDP1442
      CDIGAM = CDIGAM + ONE/U + PI*COS(H)/SIN(H)                        ACDP1443
      RETURN                                                            ACDP1444
C                                                                       ACDP1445
  110 WRITE(LERR,1000) 'CDIGAM',X                                       ACDP1446
      CDIGAM=ZERO                                                       ACDP1447
      RETURN                                                            ACDP1448
C                                                                       ACDP1449
      ENTRY LOGAM(ACC)                                                  ACDP1450
C                                                                       ACDP1451
C      initialisation call for calculations to accuracy 'ACC'           ACDP1452
C                                                                       ACDP1453
      NX0 = 6                                                           ACDP1454
      X0  = NX0 + ONE                                                   ACDP1455
      PI = FOUR*ATAN(ONE)                                               ACDP1456
      ALPI = LOG(PI)                                                    ACDP1457
      HL2P = LOG(TWO*PI) * HALF                                         ACDP1458
      ACCUR = ACC                                                       ACDP1459
      DO 120 K=1,NB                                                     ACDP1460
       F21 = K*2 - ONE                                                  ACDP1461
       B(K) = BN(K) / (BD(K) * K*TWO * F21)                             ACDP1462
       ERR = ABS(B(K)) * K*TWO / X0**F21                                ACDP1463
  120 IF(ERR.LT.ACC) GO TO 130                                          ACDP1464
       NX0 = INT((ERR/ACC)**(ONE/F21) * X0)                             ACDP1465
       K = NB                                                           ACDP1466
  130 NT = K                                                            ACDP1467
C     print *,' logam requires k = ',k ,' with cutoff at x =',nx0+1     ACDP1468
      RETURN                                                            ACDP1469
      END                                                               ACDP1470
      FUNCTION TIDY(Z,ACC)                                              ACDP1471
C                     TIDY A COMPLEX NUMBER                             ACDP1472
      REAL*8 X,Y,ACC,AZ                                                 ACDP1473
      COMPLEX*16 Z,TIDY                                                 ACDP1474
C                                                                       ACDP1475
      X = REAL(Z)                                                       ACDP1476
      Y = IMAG(Z)                                                       ACDP1477
      AZ= (ABS(X) + ABS(Y)) * ACC * 5                                   ACDP1478
      IF(ABS(X) .LT. AZ) X = 0D+0                                       ACDP1479
      IF(ABS(Y) .LT. AZ) Y = 0D+0                                       ACDP1480
      TIDY = DCMPLX(X,Y)                                                ACDP1481
      RETURN                                                            ACDP1482
      END                                                               ACDP1483
CCCTEST OF COULCC27                                                     ACDP1484
      IMPLICIT REAL*8 (A-H,O-Z)                                         ACDP1485
      COMPLEX*16 X,ETA,ZLMIN,FC(201),GC(201),FCP(201),GCP(201),         ACDP1486
     X           SIG(201),ZL,WS,CI                                      ACDP1487
      INTEGER LDISP(8)                                                  ACDP1488
      LOGICAL WHIT                                                      ACDP1489
      COMMON       /STEED/ RERR,NFP(2),NPQ(3),KASE(2)                   ACDP1490
      CHARACTER*20 NOTE                                                 ACDP1491
      CHARACTER*4 WHO(2,4,3),IRREG,REG                                  ACDP1492
      DATA ZERO,HALF,ONE,FOUR,CI / 0D+0,0.5D+0,1D+0,4D+0,(0D+0,1D+0) /  ACDP1493
      DATA WHO / 'F','G','j','y','J','Y','I','K' ,                      ACDP1494
     X           'F','H+','j','h(1)','J','H(1)','?','?' ,               ACDP1495
     X           'F','H-','j','h(2)','J','H(2)','?','?' /               ACDP1496
      DATA LDISP / 1,2,4,11,31,101,301,1001 /                           ACDP1497
C                                                                       ACDP1498
      PI = FOUR * ATAN(ONE)                                             ACDP1499
      WRITE(6,1000)                                                     ACDP1500
C                                                                       ACDP1501
   10 READ(5,*,END=40) X,ETA,ZLMIN,NL,MODE,KFN,WHIT,NOTE                ACDP1502
      IF(NL.LE.0) GO TO 40                                              ACDP1503
      IFAIL = 1                                                         ACDP1504
      MD = MOD(ABS(MODE),10)                                            ACDP1505
      IH =     ABS(MODE)/10                                             ACDP1506
      KFIN = KFN                                                        ACDP1507
      IF(WHIT) KFN=MAX(KFN,0)                                           ACDP1508
      WRITE(6,1010) X,ETA,ZLMIN,NL,MODE,KFN,NOTE                        ACDP1509
C                                                                       ACDP1510
      CALL COULCC(X,ETA,ZLMIN,NL,FC,GC,FCP,GCP,SIG,MODE,KFN,IFAIL)      ACDP1511
C                                                                       ACDP1512
      WRITE(6,1020) IFAIL,RERR,NFP,NPQ,KASE                             ACDP1513
      IF(IFAIL.LT.0) GO TO 30                                           ACDP1514
      DO 20 I=1,8                                                       ACDP1515
      L = LDISP(I)                                                      ACDP1516
      IF(L.GT.NL-IFAIL) GO TO 20                                        ACDP1517
         ZL = ZLMIN + L - 1                                             ACDP1518
         IF(KFN.NE.0) SIG(L) = ZERO                                     ACDP1519
         IRREG = WHO(2,MAX(KFN+1,1),IH+1)                               ACDP1520
           REG = WHO(1,MAX(KFN+1,1),1)                                  ACDP1521
         IF(WHIT) THEN                                                  ACDP1522
            IRREG = 'WHIT'                                              ACDP1523
            WS = EXP(-HALF*PI*(ETA - CI*ZL)  - CI*SIG(L))               ACDP1524
            GC(L)  = WS * GC(L)                                         ACDP1525
            IF(MD.EQ.1) GCP(L) = CI*WS * GCP(L)                         ACDP1526
            FC(L)  = CI/WS * FC(L)                                      ACDP1527
            IF(MOD(MD,2).EQ.1) FCP(L)  = FCP(L) / WS                    ACDP1528
             REG = 'WH-F'                                               ACDP1529
           ENDIF                                                        ACDP1530
      WRITE(6,1030) ZL,REG,FC(L),IRREG,GC(L)                            ACDP1531
      IF(MD.EQ.1)WRITE(6,1040)            FCP(L),GCP(L)                 ACDP1532
      IF(SIG(L).NE.ZERO.AND.KFIN.EQ.0) WRITE(6,1050) SIG(L)             ACDP1533
   20 CONTINUE                                                          ACDP1534
   30 CONTINUE                                                          ACDP1535
      GO TO 10                                                          ACDP1536
   40 STOP                                                              ACDP1537
 1000 FORMAT('1TEST OF THE CONTINUED-FRACTION COULOMB & BESSEL ROUTINES'ACDP1538
     X /)                                                               ACDP1539
 1010 FORMAT(/'0X =',F12.4,F9.4,', ETA =',2F8.3,', ZLMIN =',2F8.3,'  NL ACDP1540
     X=',I4,  '  MODE = ',I3,'  KFN =',I2,8X,A20)                       ACDP1541
 1020 FORMAT(' COULCC  :: IFAIL =',I4,                                  ACDP1542
     X   '  RERR =',1P,E10.1,'  ITS =',I6,I4,2I6,I4,2I2/)               ACDP1543
 1030 FORMAT(' ZL =',2F8.3,' :: FC = ',A4,' =',1P,2D20.12,',  GC = ',   ACDP1544
     X   A4,' =',2D20.12)                                               ACDP1545
 1040 FORMAT(24X,' FC''=       ',1P,2D20.12,',  GC'' =',6X ,2D20.12)    ACDP1546
 1050 FORMAT(24X,' SIG=   ',   2F20.12)                                 ACDP1547
      END                                                               ACDP1548
//L.SYSPRINT DD DUMMY                                                   ACDP1549
//G.SYSIN DD *                                                          ACDP1550
(0,10),(0,0),(0,0),2,12,2,F,   'ARDILL: J,H(1)'                         ACDP1551
(0,15),(0,0),(0,0),2,12,2,F,   'ARDILL: J,H(1),?'                       ACDP1552
(.01,0),(0,0),(0.2,0),11,-1,3,F,   'CAMPBELL: I,K'                      ACDP1553
(12.2,13.3),(0,0),(0.1,0),31,-2,3,F,   'CAMPBELL: I,K'                  ACDP1554
(0,19.2),(0,0),(0.728,0),11,-2,3,F,   'CAMPBELL: I'                     ACDP1555
(.01,0),(0,0),(0.10,0),11,2,2,F,   'CAMPBELL: J,Y'                      ACDP1556
(20,-.001),(0,0),(0,0),11,2,1,F,   'LENTZ: J,Y SPH.BESSEL'              ACDP1557
(20,-100),(0,0),(0,0),11,2,1,F,   'LENTZ: J,Y SPH.BESSEL'               ACDP1558
(1000,-10),(0,0),(0,0),11,2,1,F,   'LENTZ: J,Y SPH.BESSEL'              ACDP1559
(20,0),(15,0),(0,0),1,2,0,F,   'KLEIN: WKB#1'                           ACDP1560
(80,0),(50,0),(0,0),1,2,0,F,   'KLEIN: WKB#2'                           ACDP1561
(20,0),(10,0),(1,0),1,2,0,F,   'KLEIN: WKB#3'                           ACDP1562
(.001,0),(50,0),(0,0),1,2,0,F,   'KLEIN: WKB'                           ACDP1563
(0.238,0),(50,0),(1,0),1,2,0,F,   'KLEIN: WKB'                          ACDP1564
(1999,0),(1000,0),(0,0),1,2,0,F,   'KLEIN: WKB'                         ACDP1565
(10,0),(10,0),(10,0),1,1,0,F,   'KLEIN'                                 ACDP1566
(8.4751695,-.28219245),(3.3281266,.11081456),(5,0),1,1,0,F,'TAMURA #1,3'ACDP1567
(4.914581,-.48663942), (5.6899111,.56341223),(10,0),1,1,0,F, 'TAMURA #6'ACDP1568
(45,0),(9,0),(18,0),1,1,0,F, 'TAKEMASA 1.1'                             ACDP1569
(45,0),(9,0),(17.3083637,3.2538931),1,1,0,F, 'TAKEMASA 1.2'             ACDP1570
(45,0),(9,0),(6.000,6.92820323),  1,1,0,F, 'TAKEMASA 1.6'               ACDP1571
(25,0),(9,0),(17.3083637,3.2538931),1,1,0,F, 'TAKEMASA 2.2'             ACDP1572
(25,0),(9,0),(6.000,6.92820323),  1,1,0,F, 'TAKEMASA 2.6'               ACDP1573
(.5,0),(.25,0),(0,0),1,1,-1,F, 'BARNETT: COULFG H.P.'                   ACDP1574
(1.0,0),(1.,0),(0,0),1,1,-1,F, 'BARNETT: COULFG H.P.'                   ACDP1575
(68.,0),(32,0),(0,0),1,1,-1,F, 'BARNETT: COULFG H.P.'                   ACDP1576
(0,37.5),(0,6.66666666666),(0,0),31,12,-1,T,'BELL&SCOTT: WHIT#3'        ACDP1577
(0,12.8),(0,5),(0,0),11,12,0,T,'B&S#2 BOUND STATES'                     ACDP1578
(0,31.62277660),(0,.6324555321),(0,0),11,12,-1,T, 'NOBLE & IJT, #1'     ACDP1579
(0,26.47640459),(0,.7553895746),(0,0),11,11,-1,T, 'NOBLE & IJT, #2'     ACDP1580
(0,3.984184201314),(0,13.8045826249),(0,0),11,12,-1,T, 'NOBLE & IJT, #6'ACDP1581
(0,0.1),(0,-1),(0,0),3,12,-1,T,  'HEBBARD&ROBSON: WHIT'                 ACDP1582
(0,0.1),(0,-10),(0,0),2,12,-1,T,  'HEBBARD&ROBSON: WHIT'                ACDP1583
(0,3.0),(0,-1),(0,0),3,12,-1,T,  'HEBBARD&ROBSON: WHIT'                 ACDP1584
(0,3.0),(0,-10),(0,0),3,12,-1,T,  'HEBBARD&ROBSON: WHIT'                ACDP1585
(0,18.),(0,-10),(0,0),3,12,-1,T,  'HEBBARD&ROBSON: WHIT'                ACDP1586
(-.7,.7),   (0,0),(0,-20),11,-2, 2,F,   'SHOULD USE -X RULES'           ACDP1587
(-1,1),   (0,0),(-30,0),61,-2, 2,F,   'SHOULD USE -L RULES'             ACDP1588
(-1,1),   (0,0),(-30.1,0),61,-2, 2,F,   'SHOULD USE -L RULES'           ACDP1589
(60,0),(50,0),(0,0),200,2,0,F,   'EXP RANGE (LARGE L)'                  ACDP1590
(0,100),(0,10),(0,0),51,-12,0,T,'POLES & INSTABITIES'                   ACDP1591
(920,0),(500,0),(1,0),1,2,0,F,   '1F1 OVERFLOW'                         ACDP1592
(30000,0),(150,0),(1,0),1,2,0,F,   'NO CF1 CONVERGENCE'                 ACDP1593
(0,0),(0,0),(0,0),0,0,0,F,   ' '                                        ACDP1594
(1920,0),(1000,0),(0,0),1,2,0,F,   'KLEIN: WKB'                         ACDP1595
// EXEC NOTIFY,OPT=COND                                                 ACDP1596
//                                                                      ACDP1597
                                                                        ACDP****
