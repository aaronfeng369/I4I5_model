C     Dec-26-2013===YF===Add updated DWDE, D2WDE2
C 	  Mar-17-2015===YF,CL===change mu to dmu, ka to dka
C							use DWDE and D2W2DE2 symmetric properties
C     Mar-18-2015===YF,CL===modify the component for simplification
C     Mar-31-2015===YF===Modify the Gaussian point = (sample element + indenter element)*8
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TIME(2)
C      
      COMMON RTOTAL(8,9)
      COMMON DROTINC(8,9)
      
      DIMENSION R0(3,3), DROT(3,3), R(3,3)

C 
     1 
     2 
C INITIALIZE COMMON ROTATION VARIABLES
C 
     1 
     2 
C At the start of the analysis.
      IF(LOP.EQ.0) THEN
C For all integration points.    
         DO K1=1,8
            RTOTAL(K1,1) = 1.D0
            RTOTAL(K1,2) = 0.D0
            RTOTAL(K1,3) = 0.D0

            RTOTAL(K1,4) = 0.D0
            RTOTAL(K1,5) = 1.D0
            RTOTAL(K1,6) = 0.D0
   
            RTOTAL(K1,7) = 0.D0
            RTOTAL(K1,8) = 0.D0
            RTOTAL(K1,9) = 1.D0

            DROTINC(K1,1) = 1.D0
            DROTINC(K1,2) = 0.D0
            DROTINC(K1,3) = 0.D0

            DROTINC(K1,4) = 0.D0
            DROTINC(K1,5) = 1.D0
            DROTINC(K1,6) = 0.D0

            DROTINC(K1,7) = 0.D0
            DROTINC(K1,8) = 0.D0
            DROTINC(K1,9) = 1.D0
         END DO
      END IF
C           
C 
     1 
     2 
C UPDATE TOTAL ROTATION
C 
     1 
     2 
C At the end of the current increment
      IF(LOP.EQ.2) THEN
C For all integration points    
         DO K1=1,8
            R0(1,1) = RTOTAL(K1,1)
            R0(2,1) = RTOTAL(K1,2)
            R0(3,1) = RTOTAL(K1,3)

            R0(1,2) = RTOTAL(K1,4)
            R0(2,2) = RTOTAL(K1,5)
            R0(3,2) = RTOTAL(K1,6)

            R0(1,3) = RTOTAL(K1,7)
            R0(2,3) = RTOTAL(K1,8)
            R0(3,3) = RTOTAL(K1,9)

            DROT(1,1) = DROTINC(K1,1)
            DROT(2,1) = DROTINC(K1,2)
            DROT(3,1) = DROTINC(K1,3)

            DROT(1,2) = DROTINC(K1,4)
            DROT(2,2) = DROTINC(K1,5)
            DROT(3,2) = DROTINC(K1,6)

            DROT(1,3) = DROTINC(K1,7)
            DROT(2,3) = DROTINC(K1,8)
            DROT(3,3) = DROTINC(K1,9)

C Update the total rotation  
            CALL MULMAT(DROT,R0,R)
    
            RTOTAL(K1,1) = R(1,1)
            RTOTAL(K1,2) = R(2,1)
            RTOTAL(K1,3) = R(3,1)

            RTOTAL(K1,4) = R(1,2)
            RTOTAL(K1,5) = R(2,2)
            RTOTAL(K1,6) = R(3,2)

            RTOTAL(K1,7) = R(1,3)
            RTOTAL(K1,8) = R(2,3)
            RTOTAL(K1,9) = R(3,3)
         END DO
      END IF
C
      RETURN
      END
C
C ------------------------------------------------------------------------------
       SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
       INCLUDE 'ABA_PARAM.INC'
C
       CHARACTER*8 CMNAME
       DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
       COMMON RTOTAL(8,9)
       COMMON DROTINC(8,9)
       DIMENSION CBAR1(3,3), CL1(3,3,3,3), CL2(3,3,3,3), C(3,3),
     2 DEVS(3,3), DI(3,3), DWDE(3,3), D2WDE2(3,3,3,3),
     3 FBAR1(3,3), FBAR1T(3,3), F1(3,3), F2(3,3), FE(3,3), FT(3,3),
     4 E1(3,3), E2(3,3), E3(3,3), FWC(3,3), FWCF1(3,3),
     5 GI(3,3),
     6 R0(3,3), R(3,3), RT(3,3), RI(3,3), RIT(3,3),
     7 S1(3,3), S2(3,3), S3(3,3), S4(3,3), S4A(3,3), S5(3,3), S6(3,3),
     8 SOUT1(3,3), SOUT2(3,3), SYM(3,3,3,3),
     9 WC(3,3)
C ------------------------------------------------------------------------------
C MATERIAL PROPERTIES
C ------------------------------------------------------------------------------
        dmu = PROPS(1)
        ze = PROPS(2)
        phi = PROPS(3)
        dka = PROPS(4)
        A1 = PROPS(5)
        A2 = PROPS(6)
        A3 = PROPS(7)
        
C        C1 = PROPS(22)
C        D1 = PROPS(23)
        T0 = 0.D0
        TMAX = 1.D0
        AG1T = 1.D0
        AG2T = 1.D0
        AG3T = 1.D0
        A11 = 1.D0
        A21 = 0.D0
        A31 = 0.D0
        A12 = 0.D0
        A22 = 1.D0
        A32 = 0.D0
        A13 = 0.D0
        A23 = 0.D0
        A33 = 1.D0
        E = 2.718281828459046
        NLINE = 8*(NOEL - 1) + NPT
        CT = TIME(2) + DTIME
C ------------------------------------------------------------------------------
C INITIALIZE MATRICES
C ------------------------------------------------------------------------------
C DI = Identity matrix
       DO K1=1,3
          DO K2=1,3
             IF(K1.EQ.K2) THEN
                DI(K1,K2) = 1.D0
             ELSE
                DI(K1,K2) = 0.D0
             END IF
                DWDE(K1,K2) = 0.D0
                WC(K1,K2) = 0.D0
          END DO
        END DO
C Fourth order symmetric matrix
       DO K1=1,3
          DO K2=1,3
             DO K3=1,3
                DO K4=1,3
                    SYM(K1,K2,K3,K4) = (DI(K1,K4)*DI(K2,K3) +
     1 DI(K1,K3)*DI(K2,K4))/2.D0
                END DO
             END DO
          END DO
       END DO
C ------------------------------------------------------------------------------
C MATERIAL ORIENTATION
C ------------------------------------------------------------------------------
       RI(1,1) = A11
       RI(2,1) = A21
       RI(3,1) = A31
       RI(1,2) = A12
       RI(2,2) = A22
       RI(3,2) = A32
       RI(1,3) = A13
       RI(2,3) = A23
       RI(3,3) = A33
       CALL TRMAT(RI,RIT)
       CALL MULMAT(RI,DFGRD1,F1)
       CALL MULMAT(F1,RIT,F2)
C ------------------------------------------------------------------------------
C GROWTH
C ------------------------------------------------------------------------------
       IF(CT.LT.T0) THEN
       AG1 = 1.D0
       AG2 = 1.D0
       AG3 = 1.D0
       END IF
       IF(CT.GT.TMAX) THEN
       AG1 = AG1T
       AG2 = AG2T
       AG3 = AG2T
       END IF
       IF(CT.GE.T0.AND.CT.LE.TMAX) THEN
       AG1 = 1.D0 + (AG1T - 1.D0)*(CT - T0)/(TMAX - T0)
       AG2 = 1.D0 + (AG2T - 1.D0)*(CT - T0)/(TMAX - T0)
       AG3 = 1.D0 + (AG3T - 1.D0)*(CT - T0)/(TMAX - T0)
       END IF
       GI(1,1) = 1.D0/AG1
       GI(2,2) = 1.D0/AG2
       GI(3,3) = 1.D0/AG3
       GI(1,2) = 0.D0
       GI(1,3) = 0.D0
       GI(2,3) = 0.D0
       GI(2,1) = 0.D0
       GI(3,1) = 0.D0
       GI(3,2) = 0.D0       
C Remove growth part
       CALL MULMAT(F2, GI, FE)
C ------------------------------------------------------------------------------
C VOLUMETRIC CONSTANTS
C ------------------------------------------------------------------------------
C DUDJ = dU/dJ, D2UDJ2 = (d/dJ)(dU/dJ)
       CALL DETMAT(FE,AJ1)
       DUDJ = dka*(AJ1-1.D0)
       D2UDJ2 = dka
C ------------------------------------------------------------------------------
C DISTORTIONAL DEFORMATION GRADIENT
C ------------------------------------------------------------------------------
       DO K1=1,3
        DO K2=1,3
            FBAR1(K1,K2) = (AJ1**(-1.D0/3.D0))*FE(K1,K2)
        END DO
       END DO
C ------------------------------------------------------------------------------
C DISTORTIONAL STRAIN COMPONENTS
C ------------------------------------------------------------------------------
       CALL TRMAT(FBAR1, FBAR1T)
       CALL MULMAT(FBAR1T, FBAR1, CBAR1)
       
       Ebar1 = (CBAR1(1,1) - 1.D0)*0.5d0
       Ebar2 = (CBAR1(2,2) - 1.D0)*0.5d0
       Ebar3 = (CBAR1(3,3) - 1.D0)*0.5d0
       Ebar4 = CBAR1(1,2)*0.5d0
       Ebar5 = CBAR1(1,3)*0.5d0
       Ebar6 = CBAR1(2,3)*0.5d0
       
       !Ebar1 = (CBAR1(1,1) - 1.D0)/2.D0
       !Ebar2 = (CBAR1(2,2) - 1.D0)/2.D0
       !Ebar3 = (CBAR1(3,3) - 1.D0)/2.D0
       !Ebar4 = CBAR1(1,2)/2.D0
       !Ebar5 = CBAR1(1,3)/2.D0
       !Ebar6 = CBAR1(2,3)/2.D0
C Anisotropy Constant
C       AQ = b1111*Ebar1*Ebar1 + b1122*Ebar1*Ebar2 + b1133*Ebar1*Ebar3 +
C     1 b1112*Ebar1*Ebar4 + b1113*Ebar1*Ebar5 + b1123*Ebar1*Ebar6 +
C     2 b2222*Ebar2*Ebar2 + b2233*Ebar2*Ebar3 + b2212*Ebar2*Ebar4 +
C     3 b2213*Ebar2*Ebar5 + b2223*Ebar2*Ebar6 + b3333*Ebar3*Ebar3 +
C     4 b3312*Ebar3*Ebar4 + b3313*Ebar3*Ebar5 + b3323*Ebar3*Ebar6 +
C     5 b1212*Ebar4*Ebar4 + b1213*Ebar4*Ebar5 + b1223*Ebar4*Ebar6 +
C     6 b1313*Ebar5*Ebar5 + b1323*Ebar5*Ebar6 + b2323*Ebar6*Ebar6
C ------------------------------------------------------------------------------
C CAUCHY STRESS
C ------------------------------------------------------------------------------
C DWDE = dW/dE
      
      Ebar12 = Ebar1 * 2.d0
      Ebar14 = Ebar1 * 4.d0
      Ebar18 = Ebar1 * 8.d0
      
      Ebar22 = Ebar2 * 2.d0
      Ebar24 = Ebar2 * 4.d0
      Ebar28 = Ebar2 * 8.d0
      
      Ebar32 = Ebar3 * 2.d0
      Ebar34 = Ebar3 * 4.d0
      Ebar38 = Ebar3 * 8.d0
      
      Ebar42 = Ebar4 * 2.d0
      Ebar44 = Ebar4 * 4.d0
      Ebar48 = Ebar4 * 8.d0
      
      Ebar52 = Ebar5 * 2.d0
      Ebar54 = Ebar5 * 4.d0
      Ebar58 = Ebar5 * 8.d0
      
      Ebar62 = Ebar6 * 2.d0
      Ebar64 = Ebar6 * 4.d0
      Ebar68 = Ebar6 * 8.d0      
      
      A12 = A1**2
      A13 = A1**3
      A14 = A1**4
      
      A22 = A2**2
      A23 = A2**3
      A24 = A2**4
      
      A32 = A3**2
      A33 = A3**3
      A34 = A3**4
      
      A1A2 = A1*A2
      A1A3 = A1*A3
      A2A3 = A2*A3
      
      dfac1 = 4.d0 * dmu * phi
      dfac2 = 4.d0 * dmu * ze
      
      dfac3 = dmu * phi
      dfac4 = dmu * ze
      
      dfac5 = dfac2 - dfac1
      
      
      DWDE(1,1) = dmu+dfac3*(A1*(A1*(Ebar18+4.D0)+A2*Ebar44+A3
     +*Ebar54)-A12*(A1*(A1*(Ebar12+1.D0)+A2*Ebar42+A3
     +*Ebar52)+A2*(A2*(Ebar22+1.D0)+A1*Ebar42+A3*Ebar62
     +)+A3*(A3*(Ebar32+1.D0)+A1*Ebar52+A2*Ebar62))*4.d0
     ++A1A2*Ebar44+A1A3*Ebar54)*0.5d0+A12*dfac4
     +*(A1*(A1*(Ebar12+1.D0)+A2*Ebar42+A3*Ebar52)+A2*(A2
     +*(Ebar22+1.D0)+A1*Ebar42+A3*Ebar62)+A3*(A3*(Ebar32
     ++1.D0)+A1*Ebar52+A2*Ebar62)-1.D0)*2.D0
     
!      DWDE(1,1) = dmu+dmu*phi*(A1*(A1*(Ebar1*8.D0+4.D0)+A2*Ebar4*4.D0+A
!     +3*Ebar5*4.D0)-A1**2*(A1*(A1*(Ebar1*2.D0+1.D0)+A2*Ebar4*2.D0+A3
!     +*Ebar5*2.D0)+A2*(A2*(Ebar2*2.D0+1.D0)+A1*Ebar4*2.D0+A3*Ebar6*2
!     +.D0)+A3*(A3*(Ebar3*2.D0+1.D0)+A1*Ebar5*2.D0+A2*Ebar6*2.D0))*4
!     +.D0+A1*A2*Ebar4*4.D0+A1*A3*Ebar5*4.D0)*(1.D0/2.D0)+A1**2*dmu*z
!     +e*(A1*(A1*(Ebar1*2.D0+1.D0)+A2*Ebar4*2.D0+A3*Ebar5*2.D0)+A2*(A
!     +2*(Ebar2*2.D0+1.D0)+A1*Ebar4*2.D0+A3*Ebar6*2.D0)+A3*(A3*(Ebar3
!     +*2.D0+1.D0)+A1*Ebar5*2.D0+A2*Ebar6*2.D0)-1.D0)*2.D0

      DWDE(1,2) = dfac3*(A1*(A2*(Ebar14+Ebar24+4.D0)+A1*Ebar44
     ++A3*Ebar64)+A22*Ebar44+A2A3*Ebar54-A1A2*
     +(A1*(A1*(Ebar12+1.D0)+A2*Ebar42+A3*Ebar52)+A2*(A2*
     +(Ebar22+1.D0)+A1*Ebar42+A3*Ebar62)+A3*(A3*(Ebar32
     ++1.D0)+A1*Ebar52+A2*Ebar62))*4.D0)*0.5d0+A1
     +*A2*dfac4*(A1*(A1*(Ebar12+1.D0)+A2*Ebar42+A3*Ebar52
     +)+A2*(A2*(Ebar22+1.D0)+A1*Ebar42+A3*Ebar62)+A3*(
     +A3*(Ebar32+1.D0)+A1*Ebar52+A2*Ebar62)-1.D0)*2.D0
     
!      DWDE(1,2) = dmu*phi*(A1*(A2*(Ebar1*4.D0+Ebar2*4.D0+4.D0)+A1*Ebar
!     +4*4.D0+A3*Ebar6*4.D0)+A2**2*Ebar4*4.D0+A2*A3*Ebar5*4.D0-A1*A2*
!     +(A1*(A1*(Ebar1*2.D0+1.D0)+A2*Ebar4*2.D0+A3*Ebar5*2.D0)+A2*(A2*
!     +(Ebar2*2.D0+1.D0)+A1*Ebar4*2.D0+A3*Ebar6*2.D0)+A3*(A3*(Ebar3*2
!     +.D0+1.D0)+A1*Ebar5*2.D0+A2*Ebar6*2.D0))*4.D0)*(1.D0/2.D0)+A
!     +1*A2*dmu*ze*(A1*(A1*(Ebar1*2.D0+1.D0)+A2*Ebar4*2.D0+A3*Ebar5*2.0
!     +D0)+A2*(A2*(Ebar2*2.D0+1.D0)+A1*Ebar4*2.D0+A3*Ebar6*2.D0)+A3*(
!     +A3*(Ebar3*2.D0+1.D0)+A1*Ebar5*2.D0+A2*Ebar6*2.D0)-1.D0)*2.D0

      DWDE(1,3) = dfac3*(A1*(A3*(Ebar14+Ebar34+4.D0)+A1*Ebar54
     ++A2*Ebar64)+A32*Ebar54+A2A3*Ebar44-A1A3*
     +(A1*(A1*(Ebar12+1.D0)+A2*Ebar42+A3*Ebar52)+A2*(A2*
     +(Ebar22+1.D0)+A1*Ebar42+A3*Ebar62)+A3*(A3*(Ebar32
     ++1.D0)+A1*Ebar52+A2*Ebar62))*4.D0)*0.5d0+A1
     +*A3*dfac4*(A1*(A1*(Ebar12+1.D0)+A2*Ebar42+A3*Ebar52
     +)+A2*(A2*(Ebar22+1.D0)+A1*Ebar42+A3*Ebar62)+A3*(
     +A3*(Ebar32+1.D0)+A1*Ebar52+A2*Ebar62)-1.D0)*2.D0
     
!      DWDE(1,3) = dmu*phi*(A1*(A3*(Ebar1*4.D0+Ebar3*4.D0+4.D0)+A1*Ebar
!     +5*4.D0+A2*Ebar6*4.D0)+A3**2*Ebar5*4.D0+A2*A3*Ebar4*4.D0-A1*A3*
!     +(A1*(A1*(Ebar1*2.D0+1.D0)+A2*Ebar4*2.D0+A3*Ebar5*2.D0)+A2*(A2*
!     +(Ebar2*2.D0+1.D0)+A1*Ebar4*2.D0+A3*Ebar6*2.D0)+A3*(A3*(Ebar3*2
!     +.D0+1.D0)+A1*Ebar5*2.D0+A2*Ebar6*2.D0))*4.D0)*(1.D0/2.D0)+A
!     +1*A3*dmu*ze*(A1*(A1*(Ebar1*2.D0+1.D0)+A2*Ebar4*2.D0+A3*Ebar5*2.0
!     +D0)+A2*(A2*(Ebar2*2.D0+1.D0)+A1*Ebar4*2.D0+A3*Ebar6*2.D0)+A3*(
!     +A3*(Ebar3*2.D0+1.D0)+A1*Ebar5*2.D0+A2*Ebar6*2.D0)-1.D0)*2.D0
      
      DWDE(2,1) = DWDE(1,2)
!      DWDE(2,1) = dmu*phi*(A2*(A1*(Ebar1*4.D0+Ebar2*4.D0+4.D0)+A2*Ebar
!     +4*4.D0+A3*Ebar5*4.D0)+A1**2*Ebar4*4.D0+A1*A3*Ebar6*4.D0-A1*A2*
!     +(A1*(A1*(Ebar1*2.D0+1.D0)+A2*Ebar4*2.D0+A3*Ebar5*2.D0)+A2*(A2*
!     +(Ebar2*2.D0+1.D0)+A1*Ebar4*2.D0+A3*Ebar6*2.D0)+A3*(A3*(Ebar3*2
!     +.D0+1.D0)+A1*Ebar5*2.D0+A2*Ebar6*2.D0))*4.D0)*(1.D0/2.D0)+A
!     +1*A2*dmu*ze*(A1*(A1*(Ebar1*2.D0+1.D0)+A2*Ebar4*2.D0+A3*Ebar5*2.0
!     +D0)+A2*(A2*(Ebar2*2.D0+1.D0)+A1*Ebar4*2.D0+A3*Ebar6*2.D0)+A3*(
!     +A3*(Ebar3*2.D0+1.D0)+A1*Ebar5*2.D0+A2*Ebar6*2.D0)-1.D0)*2.D0

      DWDE(2,2) = dmu+dfac3*(A2*(A2*(Ebar28+4.D0)+A1*Ebar44+A3
     +*Ebar64)-A22*(A1*(A1*(Ebar12+1.D0)+A2*Ebar42+A3
     +*Ebar52)+A2*(A2*(Ebar22+1.D0)+A1*Ebar42+A3*Ebar62
     +)+A3*(A3*(Ebar32+1.D0)+A1*Ebar52+A2*Ebar62))*4.d0
     ++A1A2*Ebar44+A2A3*Ebar64)*0.5d0+A22*dfac4
     +*(A1*(A1*(Ebar12+1.D0)+A2*Ebar42+A3*Ebar52)+A2*(A2
     +*(Ebar22+1.D0)+A1*Ebar42+A3*Ebar62)+A3*(A3*(Ebar32
     ++1.D0)+A1*Ebar52+A2*Ebar62)-1.D0)*2.D0
     
!      DWDE(2,2) = dmu+dmu*phi*(A2*(A2*(Ebar2*8.D0+4.D0)+A1*Ebar4*4.D0+A
!     +3*Ebar6*4.D0)-A2**2*(A1*(A1*(Ebar1*2.D0+1.D0)+A2*Ebar4*2.D0+A3
!     +*Ebar5*2.D0)+A2*(A2*(Ebar2*2.D0+1.D0)+A1*Ebar4*2.D0+A3*Ebar6*2
!     +.D0)+A3*(A3*(Ebar3*2.D0+1.D0)+A1*Ebar5*2.D0+A2*Ebar6*2.D0))*4
!     +.D0+A1*A2*Ebar4*4.D0+A2*A3*Ebar6*4.D0)*(1.D0/2.D0)+A2**2*dmu*z
!     +e*(A1*(A1*(Ebar1*2.D0+1.D0)+A2*Ebar4*2.D0+A3*Ebar5*2.D0)+A2*(A
!     +2*(Ebar2*2.D0+1.D0)+A1*Ebar4*2.D0+A3*Ebar6*2.D0)+A3*(A3*(Ebar3
!     +*2.D0+1.D0)+A1*Ebar5*2.D0+A2*Ebar6*2.D0)-1.D0)*2.D0

      DWDE(2,3) = dfac3*(A2*(A3*(Ebar24+Ebar34+4.D0)+A1*Ebar54
     ++A2*Ebar64)+A32*Ebar64+A1A3*Ebar44-A2A3*
     +(A1*(A1*(Ebar12+1.D0)+A2*Ebar42+A3*Ebar52)+A2*(A2*
     +(Ebar22+1.D0)+A1*Ebar42+A3*Ebar62)+A3*(A3*(Ebar32
     ++1.D0)+A1*Ebar52+A2*Ebar62))*4.D0)*0.5d0+A2
     +*A3*dfac4*(A1*(A1*(Ebar12+1.D0)+A2*Ebar42+A3*Ebar52
     +)+A2*(A2*(Ebar22+1.D0)+A1*Ebar42+A3*Ebar62)+A3*(
     +A3*(Ebar32+1.D0)+A1*Ebar52+A2*Ebar62)-1.D0)*2.D0
     
!      DWDE(2,3) = dmu*phi*(A2*(A3*(Ebar2*4.D0+Ebar3*4.D0+4.D0)+A1*Ebar
!     +5*4.D0+A2*Ebar6*4.D0)+A3**2*Ebar6*4.D0+A1*A3*Ebar4*4.D0-A2*A3*
!     +(A1*(A1*(Ebar1*2.D0+1.D0)+A2*Ebar4*2.D0+A3*Ebar5*2.D0)+A2*(A2*
!     +(Ebar2*2.D0+1.D0)+A1*Ebar4*2.D0+A3*Ebar6*2.D0)+A3*(A3*(Ebar3*2
!     +.D0+1.D0)+A1*Ebar5*2.D0+A2*Ebar6*2.D0))*4.D0)*(1.D0/2.D0)+A
!     +2*A3*dmu*ze*(A1*(A1*(Ebar1*2.D0+1.D0)+A2*Ebar4*2.D0+A3*Ebar5*2.0
!     +D0)+A2*(A2*(Ebar2*2.D0+1.D0)+A1*Ebar4*2.D0+A3*Ebar6*2.D0)+A3*(
!     +A3*(Ebar3*2.D0+1.D0)+A1*Ebar5*2.D0+A2*Ebar6*2.D0)-1.D0)*2.D0
      
      DWDE(3,1) = DWDE(1,3)
!      DWDE(3,1) = dmu*phi*(A3*(A1*(Ebar1*4.D0+Ebar3*4.D0+4.D0)+A2*Ebar
!     +4*4.D0+A3*Ebar5*4.D0)+A1**2*Ebar5*4.D0+A1*A2*Ebar6*4.D0-A1*A3*
!     +(A1*(A1*(Ebar1*2.D0+1.D0)+A2*Ebar4*2.D0+A3*Ebar5*2.D0)+A2*(A2*
!     +(Ebar2*2.D0+1.D0)+A1*Ebar4*2.D0+A3*Ebar6*2.D0)+A3*(A3*(Ebar3*2
!     +.D0+1.D0)+A1*Ebar5*2.D0+A2*Ebar6*2.D0))*4.D0)*(1.D0/2.D0)+A
!     +1*A3*dmu*ze*(A1*(A1*(Ebar1*2.D0+1.D0)+A2*Ebar4*2.D0+A3*Ebar5*2.0
!     +D0)+A2*(A2*(Ebar2*2.D0+1.D0)+A1*Ebar4*2.D0+A3*Ebar6*2.D0)+A3*(
!     +A3*(Ebar3*2.D0+1.D0)+A1*Ebar5*2.D0+A2*Ebar6*2.D0)-1.D0)*2.D0
      
      DWDE(3,2) = DWDE(2,3)
!      DWDE(3,2) = dmu*phi*(A3*(A2*(Ebar2*4.D0+Ebar3*4.D0+4.D0)+A1*Ebar
!     +4*4.D0+A3*Ebar6*4.D0)+A2**2*Ebar6*4.D0+A1*A2*Ebar5*4.D0-A2*A3*
!     +(A1*(A1*(Ebar1*2.D0+1.D0)+A2*Ebar4*2.D0+A3*Ebar5*2.D0)+A2*(A2*
!     +(Ebar2*2.D0+1.D0)+A1*Ebar4*2.D0+A3*Ebar6*2.D0)+A3*(A3*(Ebar3*2
!     +.D0+1.D0)+A1*Ebar5*2.D0+A2*Ebar6*2.D0))*4.D0)*(1.D0/2.D0)+A
!     +2*A3*dmu*ze*(A1*(A1*(Ebar1*2.D0+1.D0)+A2*Ebar4*2.D0+A3*Ebar5*2.0
!     +D0)+A2*(A2*(Ebar2*2.D0+1.D0)+A1*Ebar4*2.D0+A3*Ebar6*2.D0)+A3*(
!     +A3*(Ebar3*2.D0+1.D0)+A1*Ebar5*2.D0+A2*Ebar6*2.D0)-1.D0)*2.D0

      DWDE(3,3) = dmu+dfac3*(A3*(A3*(Ebar38+4.D0)+A1*Ebar54+A2
     +*Ebar64)-A32*(A1*(A1*(Ebar12+1.D0)+A2*Ebar42+A3
     +*Ebar52)+A2*(A2*(Ebar22+1.D0)+A1*Ebar42+A3*Ebar62
     +)+A3*(A3*(Ebar32+1.D0)+A1*Ebar52+A2*Ebar62))*4.d0
     ++A1A3*Ebar54+A2A3*Ebar64)*0.5d0+A32*dfac4
     +*(A1*(A1*(Ebar12+1.D0)+A2*Ebar42+A3*Ebar52)+A2*(A2
     +*(Ebar22+1.D0)+A1*Ebar42+A3*Ebar62)+A3*(A3*(Ebar32
     ++1.D0)+A1*Ebar52+A2*Ebar62)-1.D0)*2.D0
     
!      DWDE(3,3) = dmu+dmu*phi*(A3*(A3*(Ebar3*8.D0+4.D0)+A1*Ebar5*4.D0+A
!     +2*Ebar6*4.D0)-A3**2*(A1*(A1*(Ebar1*2.D0+1.D0)+A2*Ebar4*2.D0+A3
!     +*Ebar5*2.D0)+A2*(A2*(Ebar2*2.D0+1.D0)+A1*Ebar4*2.D0+A3*Ebar6*2
!     +.D0)+A3*(A3*(Ebar3*2.D0+1.D0)+A1*Ebar5*2.D0+A2*Ebar6*2.D0))*4
!     +.D0+A1*A3*Ebar5*4.D0+A2*A3*Ebar6*4.D0)*(1.D0/2.D0)+A3**2*dmu*z
!     +e*(A1*(A1*(Ebar1*2.D0+1.D0)+A2*Ebar4*2.D0+A3*Ebar5*2.D0)+A2*(A
!     +2*(Ebar2*2.D0+1.D0)+A1*Ebar4*2.D0+A3*Ebar6*2.D0)+A3*(A3*(Ebar3
!     +*2.D0+1.D0)+A1*Ebar5*2.D0+A2*Ebar6*2.D0)-1.D0)*2.D0

C
       CALL MULMAT(FBAR1, DWDE, S1)
       CALL MULMAT(S1, FBAR1T, S2)
       CALL DEV(S2, S3)
       DO K1=1,3
         DO K2=1,3
            S4(K1,K2) = S3(K1,K2)/AJ1 + DUDJ*DI(K1,K2)
            S4A(K1,K2) = S3(K1,K2)/AJ1
         END DO
       END DO
C Rotate stresses back to global system
       CALL MULMAT(RIT,S4,S5)
       CALL MULMAT(S5,RI,S6)
       STRESS(1) = S6(1,1)
       STRESS(2) = S6(2,2)
       STRESS(3) = S6(3,3)
       STRESS(4) = S6(1,2)
       STRESS(5) = S6(1,3)
       STRESS(6) = S6(2,3)
C END IF
C ------------------------------------------------------------------------------
C DDSDDE
C ------------------------------------------------------------------------------
C D2WDE2 = d2Wbar/(dEbar)-dyad-(dEbar)
C            
      
      D2WDE2(1,1,1,1) = dfac1*(A12-A14) + A14*dfac2
      !(dmu*phi*(8*A1^2 - 8*A1^4))/2 + 4*A1^4*dmu*ze
      
      D2WDE2(1,1,2,2) = dfac5*A12*A22
      !D2WDE2(1,1,2,2) = dfac2*A12*A22 - dfac1*A12*A22
      !4*A1^2*A2^2*dmu*ze - 4*A1^2*A2^2*dmu*phi
      
      D2WDE2(1,1,3,3) = dfac5*A12*A32
      !D2WDE2(1,1,3,3) = dfac2*A12*A32 - dfac1*A12*A32     
      !4*A1^2*A3^2*dmu*ze - 4*A1^2*A3^2*dmu*phi
      
      D2WDE2(1,1,1,2) = dfac2*A13*A2 - dfac3*A2*2.d0*(2.d0*A13-A1)
      !D2WDE2(1,1,1,2) = dfac2*A13*A2 - dfac3*(4.d0*A13*A2-2.d0*A1*A2)
      !4*A1^3*A2*dmu*ze - (dmu*phi*(8*A1^3*A2 - 4*A1*A2))/2
      
      D2WDE2(1,1,1,3) = dfac2*A13*A3 - dfac3*A3*2.d0*(2.d0*A13-A1)
      !D2WDE2(1,1,1,3) = dfac2*A13*A3 - dfac3*(4.d0*A13*A3-2.d0*A1*A3)
      !4*A1^3*A3*dmu*ze - (dmu*phi*(8*A1^3*A3 - 4*A1*A3))/2
      
      D2WDE2(1,1,2,3) = dfac5*A12*A2*A3
      !D2WDE2(1,1,2,3) = dfac2*A12*A2*A3 - dfac1*A12*A2*A3
      !4*A1^2*A2*A3*dmu*ze - 4*A1^2*A2*A3*dmu*phi
      
      D2WDE2(2,2,2,2) = dfac1*(A22-A24) + A24*dfac2
      !(dmu*phi*(8*A2^2 - 8*A2^4))/2 + 4*A2^4*dmu*ze
      
      D2WDE2(2,2,3,3) = dfac5*A22*A32
      !D2WDE2(2,2,3,3) = dfac2*A22*A32 - dfac1*A22*A32
      !4*A2^2*A3^2*dmu*ze - 4*A2^2*A3^2*dmu*phi
      
      D2WDE2(2,2,1,2) = dfac2*A1*A23 - dfac3*2.d0*A1*(2.d0*A23-A2)
      !D2WDE2(2,2,1,2) = dfac2*A1*A23 - dfac3*(4.d0*A1*A23-2.d0*A1*A2)
      !4*A1*A2^3*dmu*ze - (dmu*phi*(8*A1*A2^3 - 4*A1*A2))/2
      
      D2WDE2(2,2,1,3) = dfac5*A1*A22*A3
      !D2WDE2(2,2,1,3) = dfac2*A1*A22*A3 - dfac1*A1*A22*A3
      !4*A1*A2^2*A3*dmu*ze - 4*A1*A2^2*A3*dmu*phi
      
      D2WDE2(2,2,2,3) = dfac2*A23*A3 - dfac3*2.d0*A3*(2.d0*A23-A2)
      !D2WDE2(2,2,2,3) = dfac2*A23*A3 - dfac3*(4.d0*A23*A3-2.d0*A2*A3)
      !4*A2^3*A3*dmu*ze - (dmu*phi*(8*A2^3*A3 - 4*A2*A3))/2
      
      D2WDE2(3,3,3,3) = dfac1*(A32-A34) + A34*dfac2
      !(dmu*phi*(8*A3^2 - 8*A3^4))/2 + 4*A3^4*dmu*ze
                           
      D2WDE2(3,3,1,2) = dfac5*A1*A2*A32
      !D2WDE2(3,3,1,2) = dfac2*A1*A2*A32 - dfac2*A1*A2*A32   
      !4*A1*A2*A3^2*dmu*ze - 4*A1*A2*A3^2*dmu*phi

      D2WDE2(3,3,1,3) = dfac2*A1*A33 - dfac3*2.d0*A1*(2.d0*A33-A3)
      !D2WDE2(3,3,1,3) = dfac2*A1*A33 - dfac3*(4.d0*A1*A33-2.d0*A1*A3)
      !4*A1*A3^3*dmu*ze - (dmu*phi*(8*A1*A3^3 - 4*A1*A3))/2

      D2WDE2(3,3,2,3) = dfac2*A2*A33 - dfac3*2.d0*A2*(2.d0*A33-A3)
      !D2WDE2(3,3,2,3) = dfac2*A2*A33 - dfac3*(4.d0*A2*A33-2.d0*A2*A3)
      !4*A2*A3^3*dmu*ze - (dmu*phi*(8*A2*A3^3 - 4*A2*A3))/2
      
      D2WDE2(1,2,1,2) = dfac5*A12*A22
      !D2WDE2(1,2,1,2) = dfac2*A12*A22 - dfac1*A12*A22
      !4*A1^2*A2^2*dmu*ze - 4*A1^2*A2^2*dmu*phi
      
      D2WDE2(1,2,1,3) = dfac5*A12*A2*A3
      !D2WDE2(1,2,1,3) = dfac2*A12*A2*A3 - dfac1*A12*A2*A3 
      !4*A1^2*A2*A3*dmu*ze - 4*A1^2*A2*A3*dmu*phi
      
      D2WDE2(1,2,2,3) = A1*A3*(dfac3*(2.d0-4.d0*A22) + dfac2*A22)
!      D2WDE2(1,2,2,3) = dfac3*(2.d0*A1*A3-4.d0*A1*A22*A3) + 
!     &                  dfac2*A1*A22*A3
      !(dmu*phi*(4*A1*A3 - 8*A1*A2^2*A3))/2 + 4*A1*A2^2*A3*dmu*ze
      
      D2WDE2(1,3,1,3) = dfac5*A12*A32
      !D2WDE2(1,3,1,3) = dfac2*A12*A32 - dfac1*A12*A32
      !4*A1^2*A3^2*dmu*ze - 4*A1^2*A3^2*dmu*phi
      
      D2WDE2(1,3,2,3) = dfac5*A1*A2*A32
      !D2WDE2(1,3,2,3) = dfac2*A1*A2*A32 - dfac1*A1*A2*A32
      !4*A1*A2*A3^2*dmu*ze - 4*A1*A2*A3^2*dmu*phi
      
      D2WDE2(2,3,2,3) = dfac5*A22*A32
      !D2WDE2(2,3,2,3) = dfac2*A22*A32 - dfac1*A22*A32
      !4*A2^2*A3^2*dmu*ze - 4*A2^2*A3^2*dmu*phi
      
      
      D2WDE2(2,2,1,1) = D2WDE2(1,1,2,2)
      
      D2WDE2(3,3,1,1) = D2WDE2(1,1,3,3)
      
      D2WDE2(1,1,2,1) = D2WDE2(1,1,1,2)
      D2WDE2(1,2,1,1) = D2WDE2(1,1,1,2)
      D2WDE2(2,1,1,1) = D2WDE2(1,1,1,2)
      
      D2WDE2(1,1,3,1) = D2WDE2(1,1,1,3)
      D2WDE2(1,3,1,1) = D2WDE2(1,1,1,3)
      D2WDE2(3,1,1,1) = D2WDE2(1,1,1,3)
      
      D2WDE2(1,1,3,2) = D2WDE2(1,1,2,3)
      D2WDE2(2,3,1,1) = D2WDE2(1,1,2,3)
      D2WDE2(3,2,1,1) = D2WDE2(1,1,2,3)
      
      D2WDE2(3,3,2,2) = D2WDE2(2,2,3,3)
      
      D2WDE2(2,2,2,1) = D2WDE2(2,2,1,2)
      D2WDE2(1,2,2,2) = D2WDE2(2,2,1,2)
      D2WDE2(2,1,2,2) = D2WDE2(2,2,1,2)
      
      D2WDE2(2,2,3,1) = D2WDE2(2,2,1,3)
      D2WDE2(1,3,2,2) = D2WDE2(2,2,1,3)
      D2WDE2(3,1,2,2) = D2WDE2(2,2,1,3)
      
      D2WDE2(2,2,3,2) = D2WDE2(2,2,2,3)
      D2WDE2(2,3,2,2) = D2WDE2(2,2,2,3)
      D2WDE2(3,2,2,2) = D2WDE2(2,2,2,3)
            
      D2WDE2(3,3,2,1) = D2WDE2(3,3,1,2)
      D2WDE2(1,2,3,3) = D2WDE2(3,3,1,2)
      D2WDE2(2,1,3,3) = D2WDE2(3,3,1,2)
      
      D2WDE2(3,3,3,1) = D2WDE2(3,3,1,3)
      D2WDE2(1,3,3,3) = D2WDE2(3,3,1,3)
      D2WDE2(3,1,3,3) = D2WDE2(3,3,1,3)
      
      D2WDE2(3,3,3,2) = D2WDE2(3,3,2,3)
      D2WDE2(2,3,3,3) = D2WDE2(3,3,2,3)
      D2WDE2(3,2,3,3) = D2WDE2(3,3,2,3)
      
      D2WDE2(1,2,2,1) = D2WDE2(1,2,1,2)
      D2WDE2(2,1,1,2) = D2WDE2(1,2,1,2)
      D2WDE2(2,1,2,1) = D2WDE2(1,2,1,2)
      
      D2WDE2(1,2,3,1) = D2WDE2(1,2,1,3)
      D2WDE2(2,1,1,3) = D2WDE2(1,2,1,3)
      D2WDE2(2,1,3,1) = D2WDE2(1,2,1,3)
      D2WDE2(1,3,1,2) = D2WDE2(1,2,1,3)
      D2WDE2(1,3,2,1) = D2WDE2(1,2,1,3)
      D2WDE2(3,1,1,2) = D2WDE2(1,2,1,3)
      D2WDE2(3,1,2,1) = D2WDE2(1,2,1,3)
      
      D2WDE2(1,2,3,2) = D2WDE2(1,2,2,3)
      D2WDE2(2,1,2,3) = D2WDE2(1,2,2,3)
      D2WDE2(2,1,3,2) = D2WDE2(1,2,2,3)
      D2WDE2(2,3,1,2) = D2WDE2(1,2,2,3)
      D2WDE2(2,3,2,1) = D2WDE2(1,2,2,3)
      D2WDE2(3,2,1,2) = D2WDE2(1,2,2,3)
      D2WDE2(3,2,2,1) = D2WDE2(1,2,2,3)
      
      D2WDE2(1,3,3,1) = D2WDE2(1,3,1,3)
      D2WDE2(3,1,1,3) = D2WDE2(1,3,1,3)
      D2WDE2(3,1,3,1) = D2WDE2(1,3,1,3)
      
      D2WDE2(1,3,3,2) = D2WDE2(1,3,2,3)
      D2WDE2(3,1,2,3) = D2WDE2(1,3,2,3)
      D2WDE2(3,1,3,2) = D2WDE2(1,3,2,3)
      D2WDE2(2,3,1,3) = D2WDE2(1,3,2,3)
      D2WDE2(2,3,3,1) = D2WDE2(1,3,2,3)
      D2WDE2(3,2,1,3) = D2WDE2(1,3,2,3)
      D2WDE2(3,2,3,1) = D2WDE2(1,3,2,3)
      
      D2WDE2(2,3,3,2) = D2WDE2(2,3,2,3)
      D2WDE2(3,2,2,3) = D2WDE2(2,3,2,3)
      D2WDE2(3,2,3,2) = D2WDE2(2,3,2,3)
      
C ------------------------------------------------------------------------------
C TENSOR OF ELASTICITY LOOP
C ------------------------------------------------------------------------------
C DEVIATORIC STRESS DEVS = DEV(S4)
       CALL DEV(S4,DEVS)
C WC = D2WDE2:C     
       DO K1=1,3
        DO K2=1,3
            DO K3=1,3
                DO K4=1,3
                WC(K1,K2) = WC(K1,K2) + D2WDE2(K1,K2,K3,K4)*CBAR1(K3,K4)
                END DO
            END DO
        END DO
       END DO
C FWCF1 = F.(D2WDE2:C).TRANSPOSE(F)
       CALL MULMAT(FBAR1,WC,FWC)
       CALL MULMAT(FWC,FBAR1T,FWCF1)
C FULL COMPONENT LOOP
       DO K1=1,3
        DO K2=1,3
          DO K3=1,3
            DO K4=1,3
                SUMCSP1 = 0.D0
                SUMCWC1 = 0.D0
                DO L1=1,3
                    SUMCSP2 = 0.D0
                    DO L2=1,3
                        SUMCSP3 = 0.D0
                        SUMCWC2 = 0.D0
                        DO L3=1,3
                            SUMCSP4 = 0.D0
                                DO L4=1,3
                                    SUMCSP4 = SUMCSP4 +
     1 D2WDE2(L1,L2,L3,L4)*FBAR1(K4,L4)
                                    SUMCWC2 = SUMCWC2 +
     1 D2WDE2(L1,L2,L3,L4)*CBAR1(L3,L4)
                                END DO
                            SUMCSP3 = SUMCSP3 + SUMCSP4*FBAR1(K3,L3)
                        END DO
                        SUMCSP2 = SUMCSP2 + SUMCSP3*FBAR1(K2,L2)
                        SUMWCW1 = SUMWCW1 + SUMWCW2*CBAR1(L1,L2)
                    END DO
                    SUMCSP1 = SUMCSP1 + SUMCSP2*FBAR1(K1,L1)
                END DO
                CSP = SUMCSP1/AJ1
                CWC = DI(K1,K2)*DI(K3,K4)*SUMCWC1/(9.D0*AJ1)
                FWCF = DI(K1,K2)*FWCF1(K3,K4)/(3.D0*AJ1)
     1 + FWCF1(K1,K2)*DI(K3,K4)/(3.D0*AJ1)
                CS1 = (2.D0/3.D0)*(DEVS(K1,K2)*DI(K3,K4) +
     1 DI(K1,K2)*DEVS(K3,K4))
                CS2 = (2.D0/3.D0)*(S2(1,1) +
     1 S2(2,2) + S2(3,3))*(SYM(K1,K2,K3,K4) -
     2 (1.D0/3.D0)*DI(K1,K2)*DI(K3,K4))/AJ1
                U1 = DUDJ*(DI(K1,K2)*DI(K3,K4) -
     1 2.D0*SYM(K1,K2,K3,K4))
                U2 = AJ1*D2UDJ2*DI(K1,K2)*DI(K3,K4)
                DSIG = (S4(K1,K4)*DI(K2,K3) + S4(K2,K3)*DI(K1,K4) +
     1 S4(K1,K3)*DI(K2,K4) + S4(K2,K4)*DI(K1,K3))/2.D0
                CL1(K1,K2,K3,K4) = CSP + CWC - FWCF - CS1 +
     1 CS2 + U1 + U2 + DSIG
            END DO
          END DO
        END DO
       END DO
C Rotate CL back to global system
       DO K1=1,3
        DO K2=1,3
         DO K3=1,3
            DO K4=1,3
                SUM1 = 0.D0
                DO L1=1,3
                    SUM2 = 0.D0
                    DO L2=1,3
                        SUM3 = 0.D0
                        DO L3=1,3
                            SUM4 = 0.D0
                            DO L4=1,3
                                SUM4 = SUM4 + CL1(L1,L2,L3,L4)*RI(L4,K4)
                            END DO
                            SUM3 = SUM3 + SUM4*RI(L3,K3)
                        END DO
                        SUM2 = SUM2 + SUM3*RI(L2,K2)
                    END DO
                    SUM1 = SUM1 + SUM2*RI(L1,K1)
                END DO
                CL2(K1,K2,K3,K4) = SUM1
            END DO
          END DO
        END DO
       END DO
C Components of DDSDDE
        DDSDDE(1,1) = CL2(1,1,1,1)
        DDSDDE(1,2) = CL2(1,1,2,2)
        DDSDDE(1,3) = CL2(1,1,3,3)
        DDSDDE(1,4) = CL2(1,1,1,2)
        DDSDDE(1,5) = CL2(1,1,1,3)
        DDSDDE(1,6) = CL2(1,1,2,3)
        DDSDDE(2,2) = CL2(2,2,2,2)
        DDSDDE(2,3) = CL2(2,2,3,3)
        DDSDDE(2,4) = CL2(2,2,1,2)
        DDSDDE(2,5) = CL2(2,2,1,3)
        DDSDDE(2,6) = CL2(2,2,2,3)
        DDSDDE(3,3) = CL2(3,3,3,3)
        DDSDDE(3,4) = CL2(3,3,1,2)
        DDSDDE(3,5) = CL2(3,3,1,3)
        DDSDDE(3,6) = CL2(3,3,2,3)
        DDSDDE(4,4) = CL2(1,2,1,2)
        DDSDDE(4,5) = CL2(1,2,1,3)
        DDSDDE(4,6) = CL2(1,2,2,3)
        DDSDDE(5,5) = CL2(1,3,1,3)
        DDSDDE(5,6) = CL2(1,3,2,3)
        DDSDDE(6,6) = CL2(2,3,2,3)
C Fill symmetric parts of DDSDDE
       DO K1=2,6
            K3 = K1-1
            DO K2=1,K3
                DDSDDE(K1,K2) = DDSDDE(K2,K1)
            END DO
       END DO
C ------------------------------------------------------------------------------
C SPECIFIC ELASTIC STRAIN ENERGY: SSE
C ------------------------------------------------------------------------------
        SSE = (C1*(-1.D0 + EXP(AQ)))/2.D0
C ------------------------------------------------------------------------------
C OUTPUT ROTATION
C ------------------------------------------------------------------------------
C Use the A_ij bases in the initial increment.
         IF(KINC.EQ.1) THEN
         R0(1,1) = A11
         R0(2,1) = A21
         R0(3,1) = A31
         R0(1,2) = A12
         R0(2,2) = A22
         R0(3,2) = A32
         R0(1,3) = A13
         R0(2,3) = A23
         R0(3,3) = A33
C Use the updated bases in later increments.
            ELSE
          R0(1,1) = RTOTAL(NLINE,1)
          R0(2,1) = RTOTAL(NLINE,2)
          R0(3,1) = RTOTAL(NLINE,3)
          R0(1,2) = RTOTAL(NLINE,4)
          R0(2,2) = RTOTAL(NLINE,5)
          R0(3,2) = RTOTAL(NLINE,6)
          R0(1,3) = RTOTAL(NLINE,7)
          R0(2,3) = RTOTAL(NLINE,8)
          R0(3,3) = RTOTAL(NLINE,9)
         END IF
       CALL MULMAT(DROT,R0,R)
C Save initial total rotation
        RTOTAL(NLINE,1) = R0(1,1)
        RTOTAL(NLINE,2) = R0(2,1)
        RTOTAL(NLINE,3) = R0(3,1)
        RTOTAL(NLINE,4) = R0(1,2)
        RTOTAL(NLINE,5) = R0(2,2)
        RTOTAL(NLINE,6) = R0(3,2)
        RTOTAL(NLINE,7) = R0(1,3)
        RTOTAL(NLINE,8) = R0(2,3)
        RTOTAL(NLINE,9) = R0(3,3)
C Save incremental rotation
        DROTINC(NLINE,1) = DROT(1,1)
        DROTINC(NLINE,2) = DROT(2,1)
        DROTINC(NLINE,3) = DROT(3,1)
        DROTINC(NLINE,4) = DROT(1,2)
        DROTINC(NLINE,5) = DROT(2,2)
        DROTINC(NLINE,6) = DROT(3,2)
        DROTINC(NLINE,7) = DROT(1,3)
        DROTINC(NLINE,8) = DROT(2,3)
        DROTINC(NLINE,9) = DROT(3,3)
C Calculate current increment total rotation transposed.
        CALL TRMAT(R,RT)
C ------------------------------------------------------------------------------
C LOCAL STRESSES (STATEV1-STATEV6) AND LOCAL GL STRAINS (STATEV7-STATEV12)
C ------------------------------------------------------------------------------
       CALL MULMAT(RT,S4A,SOUT1)
       CALL MULMAT(SOUT1,R,SOUT2)
       STATEV(1) = SOUT2(1,1)
       STATEV(2) = SOUT2(2,2)
       STATEV(3) = SOUT2(3,3)
       STATEV(4) = SOUT2(1,2)
       STATEV(5) = SOUT2(1,3)
       STATEV(6) = SOUT2(2,3)
       CALL TRMAT(F2,FT)
       CALL MULMAT(FT,F2,C)
       DO K1=1,3
          DO K2=1,3
            IF(K1.EQ.K2) THEN
              E1(K1,K2) = (C(K1,K2) - 1.D0)/2.D0
            ELSE
              E1(K1,K2) = C(K1,K2)/2.D0
           END IF
        END DO
       END DO
       STATEV(7) = E1(1,1)
       STATEV(8) = E1(2,2)
       STATEV(9) = E1(3,3)
       STATEV(10) = E1(1,2)
       STATEV(11) = E1(1,3)
       STATEV(12) = E1(2,3)
C -----------------------------------------------------------------------------
C END UMAT
C ------------------------------------------------------------------------------
       RETURN
       END
C      
C ------------------------------------------------------------------------------
C AUXILIARY FUNCTIONS
C
C TRMAT - Gives B, the 3x3 matrix transpose of A.
C MULMAT - Gives C the 3x3 matrix multiplication of A.B.
C MATINV - Gives B, the 3x3 matrix inverse of A.
C DETMAT - Gives the determinant of a 3x3 matrix A.
C DEV - Gives B, the 3x3 devaitoric part of 3x3 matrix A.
C ------------------------------------------------------------------------------
C TRMAT
C Subroutine to transpose a 3x3 matrix A into a 3x3 matrix AT.
C INPUTS: A - 3x3 matrix.
C AT - The 3x3 matrix transpose of A.
C
       SUBROUTINE TRMAT(A,AT)      
       INCLUDE 'ABA_PARAM.INC'     
       DIMENSION A(3,3), AT(3,3)
       COMMON RTOTAL(8,9)
       COMMON DROTINC(8,9)
C
       DO K1=1,3
          DO K2=1,3
             AT(K1,K2) = A(K2,K1)
          END DO
       END DO
       
       RETURN
       END
C      
C      
C ------------------------------------------------------------------------------
C MULMAT
C Subroutine which computes the matrix C = A*B, where C, A, and B are
C 3x3 matrices.
C INPUTS: A - 3x3 matrix.
C B - 3x3 matrix.
C C - 3x3 matrix, such that C = A*B.
C
       SUBROUTINE MULMAT(A,B,C)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION A(3,3), B(3,3), C(3,3)
       COMMON RTOTAL(8,9)
       COMMON DROTINC(8,9)
C
       C(1,1) = A(1,1)*B(1,1) + A(1,2)*B(2,1) + A(1,3)*B(3,1)
       C(1,2) = A(1,1)*B(1,2) + A(1,2)*B(2,2) + A(1,3)*B(3,2)
       C(1,3) = A(1,1)*B(1,3) + A(1,2)*B(2,3) + A(1,3)*B(3,3)
       C(2,1) = A(2,1)*B(1,1) + A(2,2)*B(2,1) + A(2,3)*B(3,1)
       C(2,2) = A(2,1)*B(1,2) + A(2,2)*B(2,2) + A(2,3)*B(3,2)
       C(2,3) = A(2,1)*B(1,3) + A(2,2)*B(2,3) + A(2,3)*B(3,3)
       C(3,1) = A(3,1)*B(1,1) + A(3,2)*B(2,1) + A(3,3)*B(3,1)
       C(3,2) = A(3,1)*B(1,2) + A(3,2)*B(2,2) + A(3,3)*B(3,2)
       C(3,3) = A(3,1)*B(1,3) + A(3,2)*B(2,3) + A(3,3)*B(3,3)
C      
       RETURN
       END
C      
C      
C ------------------------------------------------------------------------------
C MATINV
C Subroutine which calculates the inverse of a 3x3 matrix.
C INPUTS: A - 3x3 initial matrix.
C AINV - 3x3 matrix, which is the inverse of A.
C
       SUBROUTINE MATINV(A,AINV)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION A(3,3), AINV(3,3), COFA(3,3), ADJA(3,3)
       COMMON RTOTAL(8,9)
       COMMON DROTINC(8,9)
C
C Compute the cofactor of A:
C
       COFA(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)
       COFA(1,2) = (-1.D0)*(A(2,1)*A(3,3) - A(3,1)*A(2,3))
       COFA(1,3) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
       COFA(2,1) = (-1.D0)*(A(1,2)*A(3,3) - A(1,3)*A(3,2))
       COFA(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
       COFA(2,3) = (-1.D0)*(A(1,1)*A(3,2) - A(1,2)*A(3,1))
       COFA(3,1) = A(1,2)*A(2,3) - A(2,2)*A(1,3)
       COFA(3,2) = (-1.D0)*(A(1,1)*A(2,3) - A(2,1)*A(1,3))
       COFA(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)
C Compute the adjoint of matrix A:
       DO K1=1,3
          DO K2=1,3
             ADJA(K1,K1) = COFA(K2,K1)
          END DO
       END DO
C Compute the determinant of A:
       DETA = A(1,1)*COFA(1,1) + A(1,2)*COFA(1,2) + A(1,3)*COFA(1,3)
C Compute the inverse of A:
       DO K1=1,3
          DO K2=1,3
             AINV(K1,K2) = (1.D0/DETA)*ADJA(K1,K2)
          END DO
       END DO
C      
       RETURN
       END
C ------------------------------------------------------------------------------
C DETMAT
C Subroutine to calculate the determinant of a 3x3 matrix
C INPUTS: A - 3x3 matrix.
C DETA - determinant of A.
C
       SUBROUTINE DETMAT(A, DETA)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION A(3,3)
       COMMON RTOTAL(8,9)
       COMMON DROTINC(8,9)
       DETA = A(1,1)*( A(2,2)*A(3,3) - A(2,3)*A(3,2) ) -
     1 A(1,2)*( A(2,1)*A(3,3) - A(2,3)*A(3,1) ) +
     2 A(1,3)*( A(2,1)*A(3,2) - A(2,2)*A(3,1) )
C     
       RETURN
       END
C      
C      
C ------------------------------------------------------------------------------
C DEV
C Subroutine to calculate the deviatoric components
C of a 3x3 matrix
C INPUTS: A - 3x3 matrix.
C DEVA - deviatoric part of A.
C
       SUBROUTINE DEV(A, DEVA)
       INCLUDE 'ABA_PARAM.INC'
       DIMENSION A(3,3), DEVA(3,3)
       COMMON RTOTAL(8,9)
       COMMON DROTINC(8,9)
C
       TRA = A(1,1) + A(2,2) + A(3,3)
       DO K1=1,3
          DO K2=1,3
             IF(K1.EQ.K2) THEN
             DEVA(K1,K2) = A(K1,K2) - (1.D0/3.D0)*TRA
             ELSE
             DEVA(K1,K2) = A(K1,K2)
             END IF
          END DO
       END DO
C     
       RETURN
       END
C ------------------------------------------------------------------------------
C END SUBROUTINES
C ------------------------------------------------------------------------------
       
       
       