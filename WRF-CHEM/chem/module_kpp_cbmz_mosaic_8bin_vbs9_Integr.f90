MODULE cbmz_mosaic_8bin_vbs9_Integrator
 USE cbmz_mosaic_8bin_vbs9_Parameters
 USE cbmz_mosaic_8bin_vbs9_Precision
 USE cbmz_mosaic_8bin_vbs9_JacobianSP
  IMPLICIT NONE
  INTEGER, PARAMETER :: ifun=1, ijac=2, istp=3, iacc=4, &
    irej=5, idec=6, isol=7, isng=8, itexit=1, ihexit=2
  CHARACTER(LEN=50), PARAMETER, DIMENSION(-8:1) :: IERR_NAMES = (/ &
    'Matrix is repeatedly singular                     ', &
    'Step size too small                               ', &
    'No of steps exceeds maximum bound                 ', &
    'Improper tolerance values                         ', &
    'FacMin/FacMax/FacRej must be positive             ', &
    'Hmin/Hmax/Hstart must be positive                 ', &
    'Selected Rosenbrock method not implemented        ', &
    'Improper value for maximal no of steps            ', &
    '                                                  ', &
    'Success                                           ' /)
CONTAINS
SUBROUTINE cbmz_mosaic_8bin_vbs9_INTEGRATE( TIN, TOUT, &
  FIX, VAR, RCONST, ATOL, RTOL, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )
   USE cbmz_mosaic_8bin_vbs9_Parameters
   IMPLICIT NONE
   REAL(kind=dp), INTENT(INOUT), DIMENSION(NFIX) :: FIX
   REAL(kind=dp), INTENT(INOUT), DIMENSION(NVAR) :: VAR
   REAL(kind=dp), INTENT(IN), DIMENSION(NSPEC) :: ATOL, RTOL
   REAL(kind=dp), INTENT(IN), DIMENSION(NREACT) :: RCONST
   REAL(kind=dp), INTENT(IN) :: TIN
   REAL(kind=dp), INTENT(IN) :: TOUT
   INTEGER, INTENT(IN), OPTIONAL :: ICNTRL_U(20)
   REAL(kind=dp), INTENT(IN), OPTIONAL :: RCNTRL_U(20)
   INTEGER, INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   REAL(kind=dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER, INTENT(OUT), OPTIONAL :: IERR_U
   REAL(kind=dp) :: STEPMIN
   INTEGER :: N_stp, N_acc, N_rej, N_sng
   SAVE N_stp, N_acc, N_rej, N_sng
   INTEGER :: i, IERR
   REAL(kind=dp) :: RCNTRL(20), RSTATUS(20)
   INTEGER :: ICNTRL(20), ISTATUS(20)
   ICNTRL(:) = 0
   RCNTRL(:) = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF
   CALL cbmz_mosaic_8bin_vbs9_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT, &
         ATOL,RTOL, &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
   STEPMIN = RCNTRL(ihexit)
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U)) IERR_U = IERR
END SUBROUTINE cbmz_mosaic_8bin_vbs9_INTEGRATE
SUBROUTINE cbmz_mosaic_8bin_vbs9_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol, &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
  USE cbmz_mosaic_8bin_vbs9_Parameters
  IMPLICIT NONE
   REAL(kind=dp), INTENT(INOUT) :: Y(NVAR)
   REAL(kind=dp), INTENT(IN), DIMENSION(NFIX) :: FIX
   REAL(kind=dp), INTENT(IN), DIMENSION(NREACT) :: RCONST
   REAL(kind=dp), INTENT(IN) :: Tstart,Tend
   REAL(kind=dp), INTENT(IN) :: AbsTol(NVAR),RelTol(NVAR)
   INTEGER, INTENT(IN) :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN) :: RCNTRL(20)
   INTEGER, INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT) :: IERR
   INTEGER, PARAMETER :: Smax = 6
   INTEGER :: Method, ros_S
   REAL(kind=dp), DIMENSION(Smax) :: ros_M, ros_E, ros_Alpha, ros_Gamma
   REAL(kind=dp), DIMENSION(Smax*(Smax-1)/2) :: ros_A, ros_C
   REAL(kind=dp) :: ros_ELO
   LOGICAL, DIMENSION(Smax) :: ros_NewF
   CHARACTER(LEN=12) :: ros_Name
  INTEGER :: Nfun,Njac,Nstp,Nacc,Nrej,Ndec,Nsol,Nsng
   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart, Hexit
   REAL(kind=dp) :: Texit
   INTEGER :: i, UplimTol, Max_no_steps
   LOGICAL :: Autonomous, VectorTol
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp
   Nfun = ISTATUS(ifun)
   Njac = ISTATUS(ijac)
   Nstp = ISTATUS(istp)
   Nacc = ISTATUS(iacc)
   Nrej = ISTATUS(irej)
   Ndec = ISTATUS(idec)
   Nsol = ISTATUS(isol)
   Nsng = ISTATUS(isng)
   Autonomous = .NOT.(ICNTRL(1) == 0)
   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
         UplimTol = NVAR
   ELSE
      VectorTol = .FALSE.
         UplimTol = 1
   END IF
   IF (ICNTRL(3) == 0) THEN
      Method = 4
   ELSEIF ( (ICNTRL(3) >= 1).AND.(ICNTRL(3) <= 5) ) THEN
      Method = ICNTRL(3)
   ELSE
      PRINT * , 'User-selected Rosenbrock method: ICNTRL(3)=', Method
      CALL cbmz_mosaic_8bin_vbs9_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL cbmz_mosaic_8bin_vbs9_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF
   Roundoff = cbmz_mosaic_8bin_vbs9_WLAMCH('E')
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL cbmz_mosaic_8bin_vbs9_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL cbmz_mosaic_8bin_vbs9_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL cbmz_mosaic_8bin_vbs9_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL cbmz_mosaic_8bin_vbs9_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL cbmz_mosaic_8bin_vbs9_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL cbmz_mosaic_8bin_vbs9_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL cbmz_mosaic_8bin_vbs9_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL cbmz_mosaic_8bin_vbs9_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO
   SELECT CASE (Method)
     CASE (1)
       CALL cbmz_mosaic_8bin_vbs9_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL cbmz_mosaic_8bin_vbs9_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL cbmz_mosaic_8bin_vbs9_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL cbmz_mosaic_8bin_vbs9_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL cbmz_mosaic_8bin_vbs9_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL cbmz_mosaic_8bin_vbs9_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT
   CALL cbmz_mosaic_8bin_vbs9_ros_Integrator(Y,Tstart,Tend,Texit, &
        AbsTol, RelTol, &
        ros_S, ros_M, ros_E, ros_A, ros_C, &
        ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, &
        Autonomous, VectorTol, Max_no_steps, &
        Roundoff, Hmin, Hmax, Hstart, Hexit, &
        FacMin, FacMax, FacRej, FacSafe, &
        IERR, &
         Nfun,Njac,Nstp,Nacc,Nrej,Ndec,Nsol,Nsng,&
         RCONST, FIX &
)
   ISTATUS(ifun) = Nfun
   ISTATUS(ijac) = Njac
   ISTATUS(istp) = Nstp
   ISTATUS(iacc) = Nacc
   ISTATUS(irej) = Nrej
   ISTATUS(idec) = Ndec
   ISTATUS(isol) = Nsol
   ISTATUS(isng) = Nsng
   RSTATUS(itexit) = Texit
   RSTATUS(ihexit) = Hexit
CONTAINS
 SUBROUTINE cbmz_mosaic_8bin_vbs9_ros_ErrorMsg(Code,T,H,IERR)
   USE cbmz_mosaic_8bin_vbs9_Precision
   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN) :: Code
   INTEGER, INTENT(OUT) :: IERR
   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:'
   IF ((Code>=-8).AND.(Code<=-1)) THEN
     PRINT *, IERR_NAMES(Code)
   ELSE
     PRINT *, 'Unknown Error code: ', Code
   ENDIF
   PRINT *, "T=", T, "and H=", H
 END SUBROUTINE cbmz_mosaic_8bin_vbs9_ros_ErrorMsg
 SUBROUTINE cbmz_mosaic_8bin_vbs9_ros_Integrator (Y, Tstart, Tend, T, &
        AbsTol, RelTol, &
        ros_S, ros_M, ros_E, ros_A, ros_C, &
        ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, &
        Autonomous, VectorTol, Max_no_steps, &
        Roundoff, Hmin, Hmax, Hstart, Hexit, &
        FacMin, FacMax, FacRej, FacSafe, &
        IERR, &
        Nfun,Njac,Nstp,Nacc,Nrej,Ndec,Nsol,Nsng, &
        RCONST, FIX &
 )
  IMPLICIT NONE
   REAL(kind=dp), INTENT(INOUT) :: Y(NVAR)
   REAL(kind=dp), INTENT(IN) :: Tstart,Tend
   REAL(kind=dp), INTENT(OUT) :: T
   REAL(kind=dp), INTENT(IN) :: AbsTol(NVAR), RelTol(NVAR)
   INTEGER, INTENT(IN) :: ros_S
   REAL(kind=dp), INTENT(IN) :: ros_M(ros_S), ros_E(ros_S), &
       ros_Alpha(ros_S), ros_A(ros_S*(ros_S-1)/2), &
       ros_Gamma(ros_S), ros_C(ros_S*(ros_S-1)/2), ros_ELO
   LOGICAL, INTENT(IN) :: ros_NewF(ros_S)
   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp), INTENT(OUT) :: Hexit
   INTEGER, INTENT(OUT) :: IERR
   REAL(kind=dp), INTENT(IN), DIMENSION(NFIX) :: FIX
   REAL(kind=dp), INTENT(IN), DIMENSION(NREACT) :: RCONST
  INTEGER, INTENT(INOUT) :: Nfun,Njac,Nstp,Nacc,Nrej,Ndec,Nsol,Nsng
   REAL(kind=dp) :: Ynew(NVAR), Fcn0(NVAR), Fcn(NVAR)
   REAL(kind=dp) :: K(NVAR*ros_S), dFdT(NVAR)
   REAL(kind=dp) :: Jac0(LU_NONZERO), Ghimj(LU_NONZERO)
   REAL(kind=dp) :: H, Hnew, HC, HG, Fac, Tau
   REAL(kind=dp) :: Err, Yerr(NVAR)
   INTEGER :: Pivot(NVAR), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp
   T = Tstart
   Hexit = 0.0_dp
   H = MIN(Hstart,Hmax)
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin
   IF (Tend >= Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF
   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.
TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )
   IF ( Nstp > Max_no_steps ) THEN
      CALL cbmz_mosaic_8bin_vbs9_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN
      CALL cbmz_mosaic_8bin_vbs9_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF
   Hexit = H
   H = MIN(H,ABS(Tend-T))
   CALL cbmz_mosaic_8bin_vbs9_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF (.NOT.Autonomous) THEN
      CALL cbmz_mosaic_8bin_vbs9_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF
   CALL cbmz_mosaic_8bin_vbs9_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)
UntilAccepted: DO
   CALL cbmz_mosaic_8bin_vbs9_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec, Nsng )
   IF (Singular) THEN
       CALL cbmz_mosaic_8bin_vbs9_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF
Stage: DO istage = 1, ros_S
       ioffset = NVAR*(istage-1)
       IF ( istage == 1 ) THEN
         CALL cbmz_mosaic_8bin_vbs9_WCOPY(NVAR,Fcn0,1,Fcn,1)
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL cbmz_mosaic_8bin_vbs9_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL cbmz_mosaic_8bin_vbs9_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL cbmz_mosaic_8bin_vbs9_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF
       CALL cbmz_mosaic_8bin_vbs9_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL cbmz_mosaic_8bin_vbs9_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL cbmz_mosaic_8bin_vbs9_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL cbmz_mosaic_8bin_vbs9_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)
   END DO Stage
   CALL cbmz_mosaic_8bin_vbs9_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL cbmz_mosaic_8bin_vbs9_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO
   CALL cbmz_mosaic_8bin_vbs9_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL cbmz_mosaic_8bin_vbs9_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = cbmz_mosaic_8bin_vbs9_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )
   Fac = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac
   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN
      Nacc = Nacc+1
      CALL cbmz_mosaic_8bin_vbs9_WCOPY(NVAR,Ynew,1,Y,1)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN
         Hnew = MIN(Hnew,H)
      END IF
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted
   ELSE
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (Nacc >= 1) THEN
         Nrej = Nrej+1
      END IF
   END IF
   END DO UntilAccepted
   END DO TimeLoop
   IERR = 1
  END SUBROUTINE cbmz_mosaic_8bin_vbs9_ros_Integrator
  REAL(kind=dp) FUNCTION cbmz_mosaic_8bin_vbs9_ros_ErrorNorm ( Y, Ynew, Yerr, &
               AbsTol, RelTol, VectorTol )
   IMPLICIT NONE
   REAL(kind=dp), INTENT(IN) :: Y(NVAR), Ynew(NVAR), &
          Yerr(NVAR), AbsTol(NVAR), RelTol(NVAR)
   LOGICAL, INTENT(IN) :: VectorTol
   REAL(kind=dp) :: Err, Scale, Ymax
   INTEGER :: i
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp
   Err = ZERO
   DO i=1,NVAR
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err = SQRT(Err/NVAR)
    cbmz_mosaic_8bin_vbs9_ros_ErrorNorm = Err
  END FUNCTION cbmz_mosaic_8bin_vbs9_ros_ErrorNorm
  SUBROUTINE cbmz_mosaic_8bin_vbs9_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   IMPLICIT NONE
   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)
   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)
   INTEGER, INTENT(INOUT) ::Nfun
   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp
   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL cbmz_mosaic_8bin_vbs9_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL cbmz_mosaic_8bin_vbs9_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL cbmz_mosaic_8bin_vbs9_WSCAL(NVAR,(ONE/Delta),dFdT,1)
  END SUBROUTINE cbmz_mosaic_8bin_vbs9_ros_FunTimeDeriv
  SUBROUTINE cbmz_mosaic_8bin_vbs9_ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular, Ndec, Nsng )
   IMPLICIT NONE
   REAL(kind=dp), INTENT(IN) :: Jac0(LU_NONZERO)
   REAL(kind=dp), INTENT(IN) :: gam
   INTEGER, INTENT(IN) :: Direction
   REAL(kind=dp), INTENT(OUT) :: Ghimj(LU_NONZERO)
   LOGICAL, INTENT(OUT) :: Singular
   INTEGER, INTENT(OUT) :: Pivot(NVAR)
   REAL(kind=dp), INTENT(INOUT) :: H
   INTEGER, INTENT(INOUT) :: Ndec, Nsng
   INTEGER :: i, ising, Nconsecutive
   REAL(kind=dp) :: ghinv
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, HALF = 0.5_dp
   Nconsecutive = 0
   Singular = .TRUE.
   DO WHILE (Singular)
     CALL cbmz_mosaic_8bin_vbs9_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL cbmz_mosaic_8bin_vbs9_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO
     CALL cbmz_mosaic_8bin_vbs9_ros_Decomp( Ghimj, Pivot, ising, Ndec )
     IF (ising == 0) THEN
        Singular = .FALSE.
     ELSE
        Nsng = Nsng+1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ising = ',ising
        IF (Nconsecutive <= 5) THEN
           H = H*HALF
        ELSE
           RETURN
        END IF
      END IF
   END DO
  END SUBROUTINE cbmz_mosaic_8bin_vbs9_ros_PrepareMatrix
  SUBROUTINE cbmz_mosaic_8bin_vbs9_ros_Decomp( A, Pivot, ising, Ndec )
   IMPLICIT NONE
   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)
   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec
CALL decomp_cbmz_mosaic_8bin_vbs9 ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1
  END SUBROUTINE cbmz_mosaic_8bin_vbs9_ros_Decomp
  SUBROUTINE cbmz_mosaic_8bin_vbs9_ros_Solve( A, Pivot, b, Nsol )
   IMPLICIT NONE
   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)
   INTEGER, INTENT(INOUT) :: nsol
   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)
   CALL cbmz_mosaic_8bin_vbs9_KppSolve( A, b )
   Nsol = Nsol+1
  END SUBROUTINE cbmz_mosaic_8bin_vbs9_ros_Solve
  SUBROUTINE cbmz_mosaic_8bin_vbs9_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
            ros_Gamma,ros_NewF,ros_ELO,ros_Name)
  IMPLICIT NONE
   INTEGER, PARAMETER :: S=2
   INTEGER, INTENT(OUT) :: ros_S
   REAL(kind=dp), DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name
    REAL(kind=dp) :: g
    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)
    ros_Name = 'ROS-2'
    ros_S = S
    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)
    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)
    ros_ELO = 2.0_dp
    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp
    ros_Gamma(1) = g
    ros_Gamma(2) =-g
 END SUBROUTINE cbmz_mosaic_8bin_vbs9_Ros2
  SUBROUTINE cbmz_mosaic_8bin_vbs9_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
           ros_Gamma,ros_NewF,ros_ELO,ros_Name)
  IMPLICIT NONE
   INTEGER, PARAMETER :: S=3
   INTEGER, INTENT(OUT) :: ros_S
   REAL(kind=dp), DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name
   ros_Name = 'ROS-3'
   ros_S = S
   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp
   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) = 0.40759956452537699824805835358067E+01_dp
   ros_C(3) = 0.92076794298330791242156818474003E+01_dp
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.
   ros_M(1) = 0.1E+01_dp
   ros_M(2) = 0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514E+00_dp
   ros_E(1) = 0.5E+00_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) = 0.22354069897811569627360909276199E+00_dp
   ros_ELO = 3.0_dp
   ros_Alpha(1)= 0.0E+00_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356E+00_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356E+00_dp
   ros_Gamma(1)= 0.43586652150845899941601945119356E+00_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314E+00_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp
  END SUBROUTINE cbmz_mosaic_8bin_vbs9_Ros3
  SUBROUTINE cbmz_mosaic_8bin_vbs9_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
           ros_Gamma,ros_NewF,ros_ELO,ros_Name)
  IMPLICIT NONE
   INTEGER, PARAMETER :: S=4
   INTEGER, INTENT(OUT) :: ros_S
   REAL(kind=dp), DIMENSION(4), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(6), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(4), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name
   REAL(kind=dp) :: g
   ros_Name = 'ROS-4'
   ros_S = S
   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156E+00_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp
   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975E+00_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626E+00_dp
   ros_C(6) =-0.6949742501781779E+00_dp
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .TRUE.
   ros_NewF(4) = .FALSE.
   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792E+00_dp
   ros_M(3) = 0.4353179431840180E+00_dp
   ros_M(4) = 0.1093502252409163E+01_dp
   ros_E(1) =-0.2815431932141155E+00_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311E+00_dp
   ros_E(4) =-0.1093502252409163E+01_dp
   ros_ELO = 4.0_dp
   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900E+00_dp
   ros_Alpha(4) = ros_Alpha(3)
   ros_Gamma(1) = 0.5728200000000000E+00_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482E+00_dp
   ros_Gamma(4) =-0.1049021087100450E+00_dp
  END SUBROUTINE cbmz_mosaic_8bin_vbs9_Ros4
  SUBROUTINE cbmz_mosaic_8bin_vbs9_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
            ros_Gamma,ros_NewF,ros_ELO,ros_Name)
  IMPLICIT NONE
   INTEGER, PARAMETER :: S=4
   INTEGER, INTENT(OUT) :: ros_S
   REAL(kind=dp), DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name
   REAL(kind=dp) :: g
   ros_Name = 'RODAS-3'
   ros_S = S
   ros_A(1) = 0.0E+00_dp
   ros_A(2) = 2.0E+00_dp
   ros_A(3) = 0.0E+00_dp
   ros_A(4) = 2.0E+00_dp
   ros_A(5) = 0.0E+00_dp
   ros_A(6) = 1.0E+00_dp
   ros_C(1) = 4.0E+00_dp
   ros_C(2) = 1.0E+00_dp
   ros_C(3) =-1.0E+00_dp
   ros_C(4) = 1.0E+00_dp
   ros_C(5) =-1.0E+00_dp
   ros_C(6) =-(8.0E+00_dp/3.0E+00_dp)
   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .FALSE.
   ros_NewF(3) = .TRUE.
   ros_NewF(4) = .TRUE.
   ros_M(1) = 2.0E+00_dp
   ros_M(2) = 0.0E+00_dp
   ros_M(3) = 1.0E+00_dp
   ros_M(4) = 1.0E+00_dp
   ros_E(1) = 0.0E+00_dp
   ros_E(2) = 0.0E+00_dp
   ros_E(3) = 0.0E+00_dp
   ros_E(4) = 1.0E+00_dp
   ros_ELO = 3.0E+00_dp
   ros_Alpha(1) = 0.0E+00_dp
   ros_Alpha(2) = 0.0E+00_dp
   ros_Alpha(3) = 1.0E+00_dp
   ros_Alpha(4) = 1.0E+00_dp
   ros_Gamma(1) = 0.5E+00_dp
   ros_Gamma(2) = 1.5E+00_dp
   ros_Gamma(3) = 0.0E+00_dp
   ros_Gamma(4) = 0.0E+00_dp
  END SUBROUTINE cbmz_mosaic_8bin_vbs9_Rodas3
  SUBROUTINE cbmz_mosaic_8bin_vbs9_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
             ros_Gamma,ros_NewF,ros_ELO,ros_Name)
  IMPLICIT NONE
   INTEGER, PARAMETER :: S=6
   INTEGER, INTENT(OUT) :: ros_S
   REAL(kind=dp), DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name
    REAL(kind=dp) :: g
    ros_Name = 'RODAS-4'
    ros_S = 6
    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp
    ros_Gamma(1) = 0.2500000000000000E+00_dp
    ros_Gamma(2) =-0.1043000000000000E+00_dp
    ros_Gamma(3) = 0.1035000000000000E+00_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp
    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826E+00_dp
    ros_A(3) = 0.2557011698983284E+00_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817E+00_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950E+00_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0E+00_dp
    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915E+00_dp
    ros_C(4) =-0.1073529058151375E+00_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp
    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0E+00_dp
    ros_M(6) = 1.0E+00_dp
    ros_E(1) = 0.0E+00_dp
    ros_E(2) = 0.0E+00_dp
    ros_E(3) = 0.0E+00_dp
    ros_E(4) = 0.0E+00_dp
    ros_E(5) = 0.0E+00_dp
    ros_E(6) = 1.0E+00_dp
    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.
    ros_ELO = 4.0_dp
  END SUBROUTINE cbmz_mosaic_8bin_vbs9_Rodas4
END SUBROUTINE cbmz_mosaic_8bin_vbs9_Rosenbrock
SUBROUTINE cbmz_mosaic_8bin_vbs9_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )
   USE cbmz_mosaic_8bin_vbs9_Parameters
   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)
   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun
   CALL cbmz_mosaic_8bin_vbs9_Fun( Y, FIX, RCONST, Ydot )
   Nfun = Nfun+1
END SUBROUTINE cbmz_mosaic_8bin_vbs9_FunTemplate
SUBROUTINE cbmz_mosaic_8bin_vbs9_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )
 USE cbmz_mosaic_8bin_vbs9_Parameters
 USE cbmz_mosaic_8bin_vbs9_Jacobian
    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)
    INTEGER :: Njac
    REAL(kind=dp) :: Jcb(LU_NONZERO)
    REAL(kind=dp) :: Told
    CALL cbmz_mosaic_8bin_vbs9_Jac_SP( Y, FIX, RCONST, Jcb )
    Njac = Njac+1
END SUBROUTINE cbmz_mosaic_8bin_vbs9_JacTemplate
SUBROUTINE cbmz_mosaic_8bin_vbs9_Fun ( V, F, RCT, Vdot )
  REAL(kind=dp) :: V(NVAR)
  REAL(kind=dp) :: F(NFIX)
  REAL(kind=dp) :: RCT(NREACT)
  REAL(kind=dp) :: Vdot(NVAR)
  REAL(kind=dp) :: A(NREACT)
  A(1) = RCT(1)*V(129)
  A(2) = RCT(2)*V(136)
  A(3) = RCT(3)*V(102)
  A(4) = RCT(4)*V(107)
  A(5) = RCT(5)*V(96)
  A(6) = RCT(6)*V(90)
  A(7) = RCT(7)*V(127)
  A(8) = RCT(8)*V(127)
  A(9) = RCT(9)*V(88)
  A(10) = RCT(10)*V(116)
  A(11) = RCT(11)*V(116)
  A(12) = RCT(12)*V(100)
  A(13) = RCT(13)*V(101)
  A(14) = RCT(14)*V(123)
  A(15) = RCT(15)*V(118)
  A(16) = RCT(16)*V(121)
  A(17) = RCT(17)*V(110)
  A(18) = RCT(18)*V(135)
  A(19) = RCT(19)*V(126)
  A(20) = RCT(20)*V(125)
  A(21) = RCT(21)*V(36)*F(2)
  A(22) = RCT(22)*V(36)*F(1)
  A(23) = RCT(23)*V(106)*F(2)
  A(24) = RCT(24)*V(106)*V(127)
  A(25) = RCT(25)*V(106)*V(129)
  A(26) = RCT(26)*V(106)*V(129)
  A(27) = RCT(27)*V(106)*V(133)
  A(28) = RCT(28)*V(127)*V(133)
  A(29) = RCT(29)*V(127)*V(129)
  A(30) = RCT(30)*V(127)*V(131)
  A(31) = RCT(31)*V(127)*V(132)
  A(32) = RCT(32)*V(131)*F(2)
  A(33) = RCT(33)*V(131)*V(133)
  A(34) = RCT(34)*V(129)*V(131)
  A(35) = RCT(35)*V(131)*V(136)
  A(36) = RCT(36)*V(102)*V(131)
  A(37) = RCT(37)*V(107)*V(131)
  A(38) = RCT(38)*V(96)*V(131)
  A(39) = RCT(39)*V(131)*V(132)
  A(40) = RCT(40)*V(88)*V(131)
  A(41) = RCT(41)*V(132)*V(132)
  A(42) = RCT(42)*V(132)*V(132)*F(1)
  A(43) = RCT(43)*V(132)*V(133)
  A(44) = RCT(44)*V(129)*V(132)
  A(45) = RCT(45)*V(129)*V(132)
  A(46) = RCT(46)*V(96)
  A(47) = RCT(47)*V(133)*V(136)
  A(48) = RCT(48)*V(129)*V(136)
  A(49) = RCT(49)*V(129)*V(136)
  A(50) = RCT(50)*V(136)*V(136)
  A(51) = RCT(51)*V(132)*V(136)
  A(52) = RCT(52)*V(90)*F(1)
  A(53) = RCT(53)*V(90)
  A(54) = RCT(54)*V(108)*V(131)
  A(55) = RCT(55)*V(37)*V(131)
  A(56) = RCT(56)*V(93)*V(131)
  A(57) = RCT(57)*V(98)*V(131)
  A(58) = RCT(58)*V(109)*V(131)
  A(59) = RCT(59)*V(116)*V(131)
  A(60) = RCT(60)*V(116)*V(136)
  A(61) = RCT(61)*V(123)*V(131)
  A(62) = RCT(62)*V(123)*V(136)
  A(63) = RCT(63)*V(118)*V(131)
  A(64) = RCT(64)*V(121)*V(131)
  A(65) = RCT(65)*V(121)*V(136)
  A(66) = RCT(66)*V(103)*V(127)
  A(67) = RCT(67)*V(103)*V(131)
  A(68) = RCT(68)*V(114)*V(127)
  A(69) = RCT(69)*V(119)*V(127)
  A(70) = RCT(70)*V(114)*V(131)
  A(71) = RCT(71)*V(119)*V(131)
  A(72) = RCT(72)*V(114)*V(136)
  A(73) = RCT(73)*V(119)*V(136)
  A(74) = RCT(74)*V(89)*V(131)
  A(75) = RCT(75)*V(91)*V(131)
  A(76) = RCT(76)*V(97)*V(133)
  A(77) = RCT(77)*V(105)*V(131)
  A(78) = RCT(78)*V(105)*V(136)
  A(79) = RCT(79)*V(92)*V(129)
  A(80) = RCT(80)*V(110)*V(131)
  A(81) = RCT(81)*V(110)*V(127)
  A(82) = RCT(82)*V(115)*V(131)
  A(83) = RCT(83)*V(115)*V(127)
  A(84) = RCT(84)*V(115)*V(136)
  A(85) = RCT(85)*V(125)*V(131)
  A(86) = RCT(86)*V(125)*V(127)
  A(87) = RCT(87)*V(125)*V(136)
  A(88) = RCT(88)*V(113)*V(133)
  A(89) = RCT(89)*V(111)*V(133)
  A(90) = RCT(90)*V(112)*V(133)
  A(91) = RCT(91)*V(113)*V(132)
  A(92) = RCT(92)*V(111)*V(132)
  A(93) = RCT(93)*V(112)*V(132)
  A(94) = RCT(94)*V(100)*V(131)
  A(95) = RCT(95)*V(101)*V(131)
  A(96) = RCT(96)*V(131)*V(135)
  A(97) = RCT(97)*V(126)*V(131)
  A(98) = RCT(98)*V(129)*V(130)
  A(99) = RCT(99)*V(87)
  A(100) = RCT(100)*V(104)*V(131)
  A(101) = RCT(101)*V(124)*V(133)
  A(102) = RCT(102)*V(122)*V(133)
  A(103) = RCT(103)*V(128)*V(133)
  A(104) = RCT(104)*V(130)*V(133)
  A(105) = RCT(105)*V(133)*V(134)
  A(106) = RCT(106)*V(120)*V(133)
  A(107) = RCT(107)*V(117)*V(133)
  A(108) = RCT(108)*V(124)*V(136)
  A(109) = RCT(109)*V(122)*V(136)
  A(110) = RCT(110)*V(128)*V(136)
  A(111) = RCT(111)*V(130)*V(136)
  A(112) = RCT(112)*V(134)*V(136)
  A(113) = RCT(113)*V(120)*V(136)
  A(114) = RCT(114)*V(117)*V(136)
  A(115) = RCT(115)*V(124)*V(132)
  A(116) = RCT(116)*V(122)*V(132)
  A(117) = RCT(117)*V(128)*V(132)
  A(118) = RCT(118)*V(130)*V(132)
  A(119) = RCT(119)*V(132)*V(134)
  A(120) = RCT(120)*V(120)*V(132)
  A(121) = RCT(121)*V(117)*V(132)
  A(122) = RCT(122)*V(48)*V(131)
  A(123) = RCT(123)*V(99)*V(109)
  A(127) = RCT(127)*V(124)
  A(128) = RCT(128)*V(122)
  A(129) = RCT(129)*V(128)
  A(130) = RCT(130)*V(130)
  A(131) = RCT(131)*V(134)
  A(132) = RCT(132)*V(120)
  A(133) = RCT(133)*V(113)
  A(134) = RCT(134)*V(111)
  A(135) = RCT(135)*V(112)
  A(136) = RCT(136)*V(117)
  A(137) = RCT(137)*V(94)*V(131)
  A(138) = RCT(138)*V(94)*V(136)
  A(139) = RCT(139)*V(94)*V(127)
  A(140) = RCT(140)*V(95)*V(131)
  A(141) = RCT(141)*V(95)*V(136)
  A(142) = RCT(142)*V(95)*V(127)
  A(143) = RCT(143)*V(4)*V(131)
  A(144) = RCT(144)*V(5)*V(131)
  A(145) = RCT(145)*V(39)*V(131)
  A(146) = RCT(146)*V(61)*V(131)
  A(147) = RCT(147)*V(6)*V(131)
  A(148) = RCT(148)*V(7)*V(131)
  A(149) = RCT(149)*V(64)*V(131)
  A(150) = RCT(150)*V(85)*V(131)
  A(151) = RCT(151)*V(38)*V(131)
  A(152) = RCT(152)*V(20)*V(131)
  A(153) = RCT(153)*V(60)*V(131)
  A(154) = RCT(154)*V(62)*V(131)
  A(155) = RCT(155)*V(40)*V(131)
  A(156) = RCT(156)*V(21)*V(131)
  A(157) = RCT(157)*V(59)*V(131)
  A(158) = RCT(158)*V(58)*V(131)
  A(159) = RCT(159)*V(41)*V(131)
  A(160) = RCT(160)*V(22)*V(131)
  A(161) = RCT(161)*V(57)*V(131)
  A(162) = RCT(162)*V(56)*V(131)
  A(163) = RCT(163)*V(42)*V(131)
  A(164) = RCT(164)*V(23)*V(131)
  A(165) = RCT(165)*V(55)*V(131)
  A(166) = RCT(166)*V(54)*V(131)
  A(167) = RCT(167)*V(43)*V(131)
  A(168) = RCT(168)*V(24)*V(131)
  A(169) = RCT(169)*V(53)*V(131)
  A(170) = RCT(170)*V(52)*V(131)
  A(171) = RCT(171)*V(44)*V(131)
  A(172) = RCT(172)*V(25)*V(131)
  A(173) = RCT(173)*V(51)*V(131)
  A(174) = RCT(174)*V(50)*V(131)
  A(175) = RCT(175)*V(45)*V(131)
  A(176) = RCT(176)*V(26)*V(131)
  A(177) = RCT(177)*V(49)*V(131)
  A(178) = RCT(178)*V(47)*V(131)
  A(179) = RCT(179)*V(46)*V(131)
  A(180) = RCT(180)*V(27)*V(131)
  A(181) = RCT(181)*V(63)*V(131)
  A(182) = RCT(182)*V(28)*V(131)
  A(183) = RCT(183)*V(84)*V(131)
  A(184) = RCT(184)*V(86)*V(131)
  A(185) = RCT(185)*V(65)*V(131)
  A(186) = RCT(186)*V(29)*V(131)
  A(187) = RCT(187)*V(83)*V(131)
  A(188) = RCT(188)*V(82)*V(131)
  A(189) = RCT(189)*V(66)*V(131)
  A(190) = RCT(190)*V(30)*V(131)
  A(191) = RCT(191)*V(81)*V(131)
  A(192) = RCT(192)*V(80)*V(131)
  A(193) = RCT(193)*V(67)*V(131)
  A(194) = RCT(194)*V(31)*V(131)
  A(195) = RCT(195)*V(79)*V(131)
  A(196) = RCT(196)*V(78)*V(131)
  A(197) = RCT(197)*V(68)*V(131)
  A(198) = RCT(198)*V(32)*V(131)
  A(199) = RCT(199)*V(77)*V(131)
  A(200) = RCT(200)*V(76)*V(131)
  A(201) = RCT(201)*V(69)*V(131)
  A(202) = RCT(202)*V(33)*V(131)
  A(203) = RCT(203)*V(75)*V(131)
  A(204) = RCT(204)*V(74)*V(131)
  A(205) = RCT(205)*V(70)*V(131)
  A(206) = RCT(206)*V(34)*V(131)
  A(207) = RCT(207)*V(73)*V(131)
  A(208) = RCT(208)*V(72)*V(131)
  A(209) = RCT(209)*V(71)*V(131)
  A(210) = RCT(210)*V(35)*V(131)
  Vdot(1) = A(55)
  Vdot(2) = 0
  Vdot(3) = 0
  Vdot(4) = 0
  Vdot(5) = 0
  Vdot(6) = 0
  Vdot(7) = 0
  Vdot(8) = 0.52*A(66)+0.22*A(68)
  Vdot(9) = 0.09*A(68)+0.16*A(69)+0.39*A(83)+0.46*A(86)+0.4*A(118)
  Vdot(10) = 0.039*A(74)+0.039*A(75)+0.039*A(77)+0.039*A(78)
  Vdot(11) = 0.108*A(74)+0.108*A(75)+0.108*A(77)+0.108*A(78)
  Vdot(12) = 0.012*A(58)
  Vdot(13) = 0.008*A(68)+0.008*A(69)+0.008*A(70)+0.008*A(71)+0.008*A(72)+0.008*A(73)
  Vdot(14) = 0.0064*A(137)+0.022*A(139)
  Vdot(15) = 0.055*A(137)+0.19*A(139)
  Vdot(16) = 0.037*A(140)+0.13*A(142)
  Vdot(17) = 0.056*A(140)+0.19*A(142)
  Vdot(18) = A(30)+A(32)+A(33)+A(34)+A(35)+A(36)+A(37)+A(38)+A(39)+A(40)+A(54)+A(55)+A(56)+A(57)+A(58)+A(59)+A(61)+A(63)&
               &+A(64)+A(67)+A(70)+A(71)+A(74)+A(75)+A(77)+A(80)+A(82)+A(85)+A(94)+A(95)+A(96)+A(97)+A(100)+A(122)
  Vdot(19) = A(137)+A(140)+A(143)+A(144)+A(145)+A(146)+A(147)+A(148)+A(149)+A(150)+A(151)+A(152)+A(153)+A(154)+A(155)&
               &+A(156)+A(157)+A(158)+A(159)+A(160)+A(161)+A(162)+A(163)+A(164)+A(165)+A(166)+A(167)+A(168)+A(169)+A(170)&
               &+A(171)+A(172)+A(173)+A(174)+A(175)+A(176)+A(177)+A(178)+A(179)+A(180)+A(181)+A(182)+A(183)+A(184)+A(185)&
               &+A(186)+A(187)+A(188)+A(189)+A(190)+A(191)+A(192)+A(193)+A(194)+A(195)+A(196)+A(197)+A(198)+A(199)+A(200)&
               &+A(201)+A(202)+A(203)+A(204)+A(205)+A(206)+A(207)+A(208)+A(209)+A(210)
  Vdot(20) = -A(152)
  Vdot(21) = -A(156)
  Vdot(22) = -A(160)
  Vdot(23) = -A(164)
  Vdot(24) = -A(168)
  Vdot(25) = -A(172)
  Vdot(26) = -A(176)
  Vdot(27) = -A(180)
  Vdot(28) = -A(182)
  Vdot(29) = -A(186)
  Vdot(30) = -A(190)
  Vdot(31) = -A(194)
  Vdot(32) = -A(198)
  Vdot(33) = -A(202)
  Vdot(34) = -A(206)
  Vdot(35) = -A(210)
  Vdot(36) = A(8)-A(21)-A(22)
  Vdot(37) = -A(55)
  Vdot(38) = -A(151)
  Vdot(39) = A(151)+A(153)
  Vdot(40) = -A(155)
  Vdot(41) = -A(159)
  Vdot(42) = -A(163)
  Vdot(43) = -A(167)
  Vdot(44) = -A(171)
  Vdot(45) = -A(175)
  Vdot(46) = -A(179)
  Vdot(47) = -A(178)+0.075*A(179)+A(180)
  Vdot(48) = -A(122)
  Vdot(49) = -A(177)+A(179)
  Vdot(50) = -A(174)+0.075*A(175)+A(176)+0.075*A(177)+A(178)
  Vdot(51) = -A(173)+A(175)+A(177)
  Vdot(52) = -A(170)+0.075*A(171)+A(172)+0.075*A(173)+A(174)
  Vdot(53) = -A(169)+A(171)+A(173)
  Vdot(54) = -A(166)+0.075*A(167)+A(168)+0.075*A(169)+A(170)
  Vdot(55) = -A(165)+A(167)+A(169)
  Vdot(56) = -A(162)+0.075*A(163)+A(164)+0.075*A(165)+A(166)
  Vdot(57) = -A(161)+A(163)+A(165)
  Vdot(58) = -A(158)+0.075*A(159)+A(160)+0.075*A(161)+A(162)
  Vdot(59) = -A(157)+A(159)+A(161)
  Vdot(60) = -A(153)+A(155)+A(157)
  Vdot(61) = 0.075*A(151)+A(152)+0.075*A(153)+A(154)
  Vdot(62) = -A(154)+0.075*A(155)+A(156)+0.075*A(157)+A(158)
  Vdot(63) = -A(181)
  Vdot(64) = A(181)+A(183)
  Vdot(65) = -A(185)
  Vdot(66) = -A(189)
  Vdot(67) = -A(193)
  Vdot(68) = -A(197)
  Vdot(69) = -A(201)
  Vdot(70) = -A(205)
  Vdot(71) = -A(209)
  Vdot(72) = -A(208)+0.075*A(209)+A(210)
  Vdot(73) = -A(207)+A(209)
  Vdot(74) = -A(204)+0.075*A(205)+A(206)+0.075*A(207)+A(208)
  Vdot(75) = -A(203)+A(205)+A(207)
  Vdot(76) = -A(200)+0.075*A(201)+A(202)+0.075*A(203)+A(204)
  Vdot(77) = -A(199)+A(201)+A(203)
  Vdot(78) = -A(196)+0.075*A(197)+A(198)+0.075*A(199)+A(200)
  Vdot(79) = -A(195)+A(197)+A(199)
  Vdot(80) = -A(192)+0.075*A(193)+A(194)+0.075*A(195)+A(196)
  Vdot(81) = -A(191)+A(193)+A(195)
  Vdot(82) = -A(188)+0.075*A(189)+A(190)+0.075*A(191)+A(192)
  Vdot(83) = -A(187)+A(189)+A(191)
  Vdot(84) = -A(183)+A(185)+A(187)
  Vdot(85) = 0.075*A(181)+A(182)+0.075*A(183)+A(184)
  Vdot(86) = -A(184)+0.075*A(185)+A(186)+0.075*A(187)+A(188)
  Vdot(87) = A(98)-A(99)
  Vdot(88) = -A(9)-A(40)+A(41)+A(42)
  Vdot(89) = -A(74)
  Vdot(90) = -A(6)+A(49)-A(52)-A(53)
  Vdot(91) = -A(75)
  Vdot(92) = 0.4*A(77)+A(78)-A(79)
  Vdot(93) = -A(56)+0.06*A(68)+0.08*A(69)
  Vdot(94) = -A(137)-A(138)-A(139)
  Vdot(95) = -A(140)-A(141)-A(142)
  Vdot(96) = -A(5)-A(38)+A(44)-A(46)
  Vdot(97) = 0.8*A(74)+0.45*A(75)-A(76)
  Vdot(98) = -A(57)+0.01*A(68)+0.01*A(69)+0.2*A(128)
  Vdot(99) = 1.98*A(18)+1.98*A(19)+1.06*A(68)+2.26*A(69)+A(70)+2.23*A(71)+0.42*A(96)+1.68*A(103)+A(106)+1.98*A(110)&
               &+A(113)-A(123)+1.25*A(129)+A(132)
  Vdot(100) = -A(12)-A(94)+A(115)
  Vdot(101) = -A(13)-A(95)+A(116)
  Vdot(102) = -A(3)+A(33)-A(36)+A(45)
  Vdot(103) = -A(66)-A(67)
  Vdot(104) = 0.03*A(68)+0.04*A(69)-A(100)+0.34*A(127)
  Vdot(105) = 0.12*A(74)+0.05*A(75)-A(77)-A(78)
  Vdot(106) = A(1)+0.89*A(2)+A(7)+A(21)-A(23)-A(24)-A(25)-A(26)-A(27)
  Vdot(107) = -A(4)+A(34)-A(37)+0.3*A(51)+2*A(52)+A(60)+A(62)+A(65)+A(78)+0.07*A(87)
  Vdot(108) = A(10)+A(11)+A(14)+A(16)+A(17)+0.33*A(20)-A(54)+A(59)+A(60)+A(65)+0.24*A(66)+0.31*A(68)+0.3*A(69)+2*A(80)&
                &+0.69*A(81)+0.07*A(83)+0.16*A(86)+0.64*A(87)+0.59*A(90)
  Vdot(109) = -A(58)+1.1*A(75)+1.86*A(87)+0.18*A(88)+1.6*A(89)+2*A(92)-A(123)+2*A(134)
  Vdot(110) = -A(17)+0.95*A(76)+0.3*A(77)-A(80)-A(81)
  Vdot(111) = A(84)-A(89)-A(92)-A(134)
  Vdot(112) = 0.5*A(85)-A(90)-A(93)-A(135)
  Vdot(113) = A(82)-A(88)-A(91)-A(133)
  Vdot(114) = -A(68)-A(70)-A(72)
  Vdot(115) = -A(82)-A(83)-A(84)
  Vdot(116) = -A(10)-A(11)+A(12)+0.2*A(20)-A(59)-A(60)+A(66)+1.56*A(67)+0.57*A(68)+A(70)+A(80)+0.7*A(81)+0.6*A(83)+0.15&
                &*A(86)+0.28*A(87)+0.63*A(88)+0.25*A(90)+0.3*A(94)+A(100)+A(101)+A(105)+0.5*A(106)+A(108)+A(112)+0.5*A(113)&
                &+0.66*A(127)+0.7*A(131)+0.5*A(132)
  Vdot(117) = 0.4*A(18)+0.41*A(19)+A(64)+A(67)+A(70)+A(71)+0.08*A(74)+0.5*A(75)+0.6*A(77)+A(80)+0.03*A(81)+0.08*A(82)&
                &+0.2*A(83)+0.2*A(85)+0.07*A(86)+0.93*A(87)+0.34*A(103)-A(107)+0.4*A(110)-A(114)-A(121)+0.24*A(129)-A(136)
  Vdot(118) = -A(15)+0.74*A(18)+0.74*A(19)+0.03*A(20)-A(63)+0.07*A(69)+0.23*A(71)+0.09*A(86)+0.63*A(90)+0.62*A(103)+0.74&
                &*A(110)+0.57*A(129)+0.15*A(131)+0.5*A(135)
  Vdot(119) = -A(69)-A(71)-A(73)
  Vdot(120) = A(72)+A(73)+A(97)-A(106)-A(113)-A(120)-A(132)
  Vdot(121) = -A(16)-A(64)-A(65)+0.04*A(68)+0.07*A(69)+0.8*A(75)+0.2*A(81)+0.85*A(86)+0.34*A(90)+0.19*A(96)+0.15*A(131)
  Vdot(122) = 0.1*A(18)+0.1*A(19)+A(57)+0.06*A(68)+0.05*A(69)+0.7*A(95)-A(102)+0.08*A(103)-A(109)+0.1*A(110)-A(116)&
                &-A(128)+0.06*A(129)
  Vdot(123) = A(13)-A(14)+0.3*A(18)+0.3*A(19)+0.07*A(20)-A(61)-A(62)+0.22*A(67)+0.47*A(68)+1.03*A(69)+A(70)+1.77*A(71)&
                &+0.03*A(81)+0.15*A(83)+0.02*A(86)+0.28*A(87)+0.8*A(89)+0.55*A(90)+0.3*A(95)+0.04*A(96)+A(102)+0.25*A(103)&
                &+0.5*A(106)+A(109)+0.3*A(110)+0.5*A(113)+A(122)+0.8*A(128)+0.21*A(129)+0.5*A(132)+A(134)+0.5*A(135)
  Vdot(124) = A(14)+A(15)+0.7*A(20)+A(56)+0.07*A(68)+0.1*A(69)+0.05*A(86)+0.7*A(94)-A(101)+A(104)-A(108)+A(111)-A(115)&
                &-A(127)+A(130)
  Vdot(125) = -A(20)+0.65*A(83)-A(85)-A(86)-A(87)+0.91*A(88)+0.2*A(89)+A(133)
  Vdot(126) = -A(19)+0.05*A(76)+A(79)+0.93*A(87)+0.09*A(88)+0.8*A(89)+A(92)-A(97)+0.16*A(103)+0.5*A(106)+0.5*A(113)&
                &+A(120)+0.5*A(132)+A(134)
  Vdot(127) = -A(7)-A(8)+A(23)-A(24)-A(28)-A(29)-A(30)-A(31)-A(66)-A(68)-A(69)-A(81)-A(83)-A(86)+0.4*A(118)-A(139)&
                &-A(142)
  Vdot(128) = A(58)+0.03*A(68)+0.09*A(69)+0.77*A(96)-A(103)-A(110)-A(117)-A(129)
  Vdot(129) = -A(1)+0.89*A(2)+A(4)+A(5)+A(6)+A(19)-A(25)-A(26)+A(27)+A(28)-A(29)-A(34)+A(35)+A(36)+A(38)+A(43)-A(44)&
                &-A(45)+A(46)+2*A(47)-A(49)+2*A(50)+0.7*A(51)+A(53)+0.95*A(76)-A(79)+0.91*A(88)+1.2*A(89)+A(90)-A(98)+A(99)&
                &+A(101)+A(102)+0.84*A(103)+A(104)+A(105)+1.5*A(106)+A(107)+A(108)+A(109)+A(110)+A(111)+A(112)+1.5*A(113)&
                &+A(114)+0.5*A(132)
  Vdot(130) = A(15)+A(16)+A(17)+0.97*A(20)+A(61)+A(62)+A(64)+A(65)+0.13*A(68)+0.19*A(69)+A(80)+0.62*A(81)+0.2*A(83)+0.5&
                &*A(85)+0.11*A(86)+0.07*A(87)-A(98)+A(99)-A(104)+A(105)-A(111)+A(112)-A(118)-A(130)+0.7*A(131)
  Vdot(131) = A(3)+A(4)+2*A(9)+A(12)+A(13)+A(18)+2*A(22)-A(30)+A(31)-A(32)-A(33)-A(34)-A(35)-A(36)-A(37)-A(38)-A(39)&
                &-A(40)+A(43)+0.7*A(51)-A(54)-A(55)-A(56)-A(57)-A(58)-A(59)-A(61)-A(63)-A(64)+0.12*A(66)-A(67)+0.33*A(68)&
                &+0.6*A(69)-A(70)-A(71)-A(74)-A(75)-A(77)-A(80)+0.08*A(81)-A(82)+0.27*A(83)-A(85)+0.27*A(86)-0.7*A(94)-0.7&
                &*A(95)-0.77*A(96)-A(97)-A(100)-A(122)-A(137)-A(140)-A(143)-A(144)-A(145)-A(146)-A(147)-A(148)-A(149)-A(150)&
                &-A(151)-A(152)-A(153)-A(154)-A(155)-A(156)-A(157)-A(158)-A(159)-A(160)-A(161)-A(162)-A(163)-A(164)-A(165)&
                &-A(166)-A(167)-A(168)-A(169)-A(170)-A(171)-A(172)-A(173)-A(174)-A(175)-A(176)-A(177)-A(178)-A(179)-A(180)&
                &-A(181)-A(182)-A(183)-A(184)-A(185)-A(186)-A(187)-A(188)-A(189)-A(190)-A(191)-A(192)-A(193)-A(194)-A(195)&
                &-A(196)-A(197)-A(198)-A(199)-A(200)-A(201)-A(202)-A(203)-A(204)-A(205)-A(206)-A(207)-A(208)-A(209)-A(210)
  Vdot(132) = A(5)+2*A(10)+A(12)+A(13)+A(14)+A(16)+A(17)+0.9*A(18)+0.9*A(19)+0.33*A(20)+A(30)-A(31)+A(32)+A(35)-A(39)&
                &+A(40)-2*A(41)-2*A(42)-A(43)-A(44)-A(45)+A(46)-A(51)+A(54)+A(55)+A(59)+A(60)+0.22*A(66)+A(67)+0.26*A(68)&
                &+0.22*A(69)+A(70)+A(71)+0.2*A(74)+0.55*A(75)+0.95*A(76)+0.6*A(77)+2*A(80)+0.76*A(81)+0.07*A(83)+0.1*A(86)&
                &+0.93*A(87)+0.91*A(88)+0.8*A(89)+A(90)-A(91)-A(92)-A(93)+A(100)+A(101)+A(102)+0.76*A(103)+0.5*A(106)+A(108)&
                &+A(109)+0.9*A(110)+0.5*A(113)-A(115)-A(116)-A(117)-A(118)-A(119)-A(120)-A(121)+A(122)+0.32*A(127)+0.6&
                &*A(128)+0.54*A(129)
  Vdot(133) = A(1)+0.11*A(2)+A(3)+A(25)-A(27)-A(28)-A(33)-A(43)-A(47)+A(48)-A(76)-A(88)-A(89)-A(90)-A(101)-A(102)-A(103)&
                &-A(104)-A(105)-A(106)-A(107)
  Vdot(134) = A(63)+0.11*A(69)-A(105)-A(112)-A(119)-A(131)
  Vdot(135) = -A(18)+A(91)+A(93)-A(96)+A(117)+A(119)
  Vdot(136) = -A(2)+A(6)+A(26)+A(29)-A(35)+A(37)-A(47)-A(48)-A(49)-2*A(50)-A(51)+A(53)-A(60)-A(62)-A(65)-A(72)-A(73)&
                &-A(78)-A(84)-A(87)-A(108)-A(109)-A(110)-A(111)-A(112)-A(113)-A(114)-A(138)-A(141)
END SUBROUTINE cbmz_mosaic_8bin_vbs9_Fun
SUBROUTINE cbmz_mosaic_8bin_vbs9_Jac_SP ( V, F, RCT, JVS )
  REAL(kind=dp) :: V(NVAR)
  REAL(kind=dp) :: F(NFIX)
  REAL(kind=dp) :: RCT(NREACT)
  REAL(kind=dp) :: JVS(LU_NONZERO)
  REAL(kind=dp) :: B(382)
  B(1) = RCT(1)
  B(2) = RCT(2)
  B(3) = RCT(3)
  B(4) = RCT(4)
  B(5) = RCT(5)
  B(6) = RCT(6)
  B(7) = RCT(7)
  B(8) = RCT(8)
  B(9) = RCT(9)
  B(10) = RCT(10)
  B(11) = RCT(11)
  B(12) = RCT(12)
  B(13) = RCT(13)
  B(14) = RCT(14)
  B(15) = RCT(15)
  B(16) = RCT(16)
  B(17) = RCT(17)
  B(18) = RCT(18)
  B(19) = RCT(19)
  B(20) = RCT(20)
  B(21) = RCT(21)*F(2)
  B(23) = RCT(22)*F(1)
  B(25) = RCT(23)*F(2)
  B(27) = RCT(24)*V(127)
  B(28) = RCT(24)*V(106)
  B(29) = RCT(25)*V(129)
  B(30) = RCT(25)*V(106)
  B(31) = RCT(26)*V(129)
  B(32) = RCT(26)*V(106)
  B(33) = RCT(27)*V(133)
  B(34) = RCT(27)*V(106)
  B(35) = RCT(28)*V(133)
  B(36) = RCT(28)*V(127)
  B(37) = RCT(29)*V(129)
  B(38) = RCT(29)*V(127)
  B(39) = RCT(30)*V(131)
  B(40) = RCT(30)*V(127)
  B(41) = RCT(31)*V(132)
  B(42) = RCT(31)*V(127)
  B(43) = RCT(32)*F(2)
  B(45) = RCT(33)*V(133)
  B(46) = RCT(33)*V(131)
  B(47) = RCT(34)*V(131)
  B(48) = RCT(34)*V(129)
  B(49) = RCT(35)*V(136)
  B(50) = RCT(35)*V(131)
  B(51) = RCT(36)*V(131)
  B(52) = RCT(36)*V(102)
  B(53) = RCT(37)*V(131)
  B(54) = RCT(37)*V(107)
  B(55) = RCT(38)*V(131)
  B(56) = RCT(38)*V(96)
  B(57) = RCT(39)*V(132)
  B(58) = RCT(39)*V(131)
  B(59) = RCT(40)*V(131)
  B(60) = RCT(40)*V(88)
  B(61) = RCT(41)*2*V(132)
  B(62) = RCT(42)*2*V(132)*F(1)
  B(64) = RCT(43)*V(133)
  B(65) = RCT(43)*V(132)
  B(66) = RCT(44)*V(132)
  B(67) = RCT(44)*V(129)
  B(68) = RCT(45)*V(132)
  B(69) = RCT(45)*V(129)
  B(70) = RCT(46)
  B(71) = RCT(47)*V(136)
  B(72) = RCT(47)*V(133)
  B(73) = RCT(48)*V(136)
  B(74) = RCT(48)*V(129)
  B(75) = RCT(49)*V(136)
  B(76) = RCT(49)*V(129)
  B(77) = RCT(50)*2*V(136)
  B(78) = RCT(51)*V(136)
  B(79) = RCT(51)*V(132)
  B(80) = RCT(52)*F(1)
  B(82) = RCT(53)
  B(83) = RCT(54)*V(131)
  B(84) = RCT(54)*V(108)
  B(85) = RCT(55)*V(131)
  B(86) = RCT(55)*V(37)
  B(87) = RCT(56)*V(131)
  B(88) = RCT(56)*V(93)
  B(89) = RCT(57)*V(131)
  B(90) = RCT(57)*V(98)
  B(91) = RCT(58)*V(131)
  B(92) = RCT(58)*V(109)
  B(93) = RCT(59)*V(131)
  B(94) = RCT(59)*V(116)
  B(95) = RCT(60)*V(136)
  B(96) = RCT(60)*V(116)
  B(97) = RCT(61)*V(131)
  B(98) = RCT(61)*V(123)
  B(99) = RCT(62)*V(136)
  B(100) = RCT(62)*V(123)
  B(101) = RCT(63)*V(131)
  B(102) = RCT(63)*V(118)
  B(103) = RCT(64)*V(131)
  B(104) = RCT(64)*V(121)
  B(105) = RCT(65)*V(136)
  B(106) = RCT(65)*V(121)
  B(107) = RCT(66)*V(127)
  B(108) = RCT(66)*V(103)
  B(109) = RCT(67)*V(131)
  B(110) = RCT(67)*V(103)
  B(111) = RCT(68)*V(127)
  B(112) = RCT(68)*V(114)
  B(113) = RCT(69)*V(127)
  B(114) = RCT(69)*V(119)
  B(115) = RCT(70)*V(131)
  B(116) = RCT(70)*V(114)
  B(117) = RCT(71)*V(131)
  B(118) = RCT(71)*V(119)
  B(119) = RCT(72)*V(136)
  B(120) = RCT(72)*V(114)
  B(121) = RCT(73)*V(136)
  B(122) = RCT(73)*V(119)
  B(123) = RCT(74)*V(131)
  B(124) = RCT(74)*V(89)
  B(125) = RCT(75)*V(131)
  B(126) = RCT(75)*V(91)
  B(127) = RCT(76)*V(133)
  B(128) = RCT(76)*V(97)
  B(129) = RCT(77)*V(131)
  B(130) = RCT(77)*V(105)
  B(131) = RCT(78)*V(136)
  B(132) = RCT(78)*V(105)
  B(133) = RCT(79)*V(129)
  B(134) = RCT(79)*V(92)
  B(135) = RCT(80)*V(131)
  B(136) = RCT(80)*V(110)
  B(137) = RCT(81)*V(127)
  B(138) = RCT(81)*V(110)
  B(139) = RCT(82)*V(131)
  B(140) = RCT(82)*V(115)
  B(141) = RCT(83)*V(127)
  B(142) = RCT(83)*V(115)
  B(143) = RCT(84)*V(136)
  B(144) = RCT(84)*V(115)
  B(145) = RCT(85)*V(131)
  B(146) = RCT(85)*V(125)
  B(147) = RCT(86)*V(127)
  B(148) = RCT(86)*V(125)
  B(149) = RCT(87)*V(136)
  B(150) = RCT(87)*V(125)
  B(151) = RCT(88)*V(133)
  B(152) = RCT(88)*V(113)
  B(153) = RCT(89)*V(133)
  B(154) = RCT(89)*V(111)
  B(155) = RCT(90)*V(133)
  B(156) = RCT(90)*V(112)
  B(157) = RCT(91)*V(132)
  B(158) = RCT(91)*V(113)
  B(159) = RCT(92)*V(132)
  B(160) = RCT(92)*V(111)
  B(161) = RCT(93)*V(132)
  B(162) = RCT(93)*V(112)
  B(163) = RCT(94)*V(131)
  B(164) = RCT(94)*V(100)
  B(165) = RCT(95)*V(131)
  B(166) = RCT(95)*V(101)
  B(167) = RCT(96)*V(135)
  B(168) = RCT(96)*V(131)
  B(169) = RCT(97)*V(131)
  B(170) = RCT(97)*V(126)
  B(171) = RCT(98)*V(130)
  B(172) = RCT(98)*V(129)
  B(173) = RCT(99)
  B(174) = RCT(100)*V(131)
  B(175) = RCT(100)*V(104)
  B(176) = RCT(101)*V(133)
  B(177) = RCT(101)*V(124)
  B(178) = RCT(102)*V(133)
  B(179) = RCT(102)*V(122)
  B(180) = RCT(103)*V(133)
  B(181) = RCT(103)*V(128)
  B(182) = RCT(104)*V(133)
  B(183) = RCT(104)*V(130)
  B(184) = RCT(105)*V(134)
  B(185) = RCT(105)*V(133)
  B(186) = RCT(106)*V(133)
  B(187) = RCT(106)*V(120)
  B(188) = RCT(107)*V(133)
  B(189) = RCT(107)*V(117)
  B(190) = RCT(108)*V(136)
  B(191) = RCT(108)*V(124)
  B(192) = RCT(109)*V(136)
  B(193) = RCT(109)*V(122)
  B(194) = RCT(110)*V(136)
  B(195) = RCT(110)*V(128)
  B(196) = RCT(111)*V(136)
  B(197) = RCT(111)*V(130)
  B(198) = RCT(112)*V(136)
  B(199) = RCT(112)*V(134)
  B(200) = RCT(113)*V(136)
  B(201) = RCT(113)*V(120)
  B(202) = RCT(114)*V(136)
  B(203) = RCT(114)*V(117)
  B(204) = RCT(115)*V(132)
  B(205) = RCT(115)*V(124)
  B(206) = RCT(116)*V(132)
  B(207) = RCT(116)*V(122)
  B(208) = RCT(117)*V(132)
  B(209) = RCT(117)*V(128)
  B(210) = RCT(118)*V(132)
  B(211) = RCT(118)*V(130)
  B(212) = RCT(119)*V(134)
  B(213) = RCT(119)*V(132)
  B(214) = RCT(120)*V(132)
  B(215) = RCT(120)*V(120)
  B(216) = RCT(121)*V(132)
  B(217) = RCT(121)*V(117)
  B(218) = RCT(122)*V(131)
  B(219) = RCT(122)*V(48)
  B(220) = RCT(123)*V(109)
  B(221) = RCT(123)*V(99)
  B(222) = 1
  B(223) = 1
  B(224) = 1
  B(225) = RCT(127)
  B(226) = RCT(128)
  B(227) = RCT(129)
  B(228) = RCT(130)
  B(229) = RCT(131)
  B(230) = RCT(132)
  B(231) = RCT(133)
  B(232) = RCT(134)
  B(233) = RCT(135)
  B(234) = RCT(136)
  B(235) = RCT(137)*V(131)
  B(236) = RCT(137)*V(94)
  B(237) = RCT(138)*V(136)
  B(238) = RCT(138)*V(94)
  B(239) = RCT(139)*V(127)
  B(240) = RCT(139)*V(94)
  B(241) = RCT(140)*V(131)
  B(242) = RCT(140)*V(95)
  B(243) = RCT(141)*V(136)
  B(244) = RCT(141)*V(95)
  B(245) = RCT(142)*V(127)
  B(246) = RCT(142)*V(95)
  B(247) = RCT(143)*V(131)
  B(248) = RCT(143)*V(4)
  B(249) = RCT(144)*V(131)
  B(250) = RCT(144)*V(5)
  B(251) = RCT(145)*V(131)
  B(252) = RCT(145)*V(39)
  B(253) = RCT(146)*V(131)
  B(254) = RCT(146)*V(61)
  B(255) = RCT(147)*V(131)
  B(256) = RCT(147)*V(6)
  B(257) = RCT(148)*V(131)
  B(258) = RCT(148)*V(7)
  B(259) = RCT(149)*V(131)
  B(260) = RCT(149)*V(64)
  B(261) = RCT(150)*V(131)
  B(262) = RCT(150)*V(85)
  B(263) = RCT(151)*V(131)
  B(264) = RCT(151)*V(38)
  B(265) = RCT(152)*V(131)
  B(266) = RCT(152)*V(20)
  B(267) = RCT(153)*V(131)
  B(268) = RCT(153)*V(60)
  B(269) = RCT(154)*V(131)
  B(270) = RCT(154)*V(62)
  B(271) = RCT(155)*V(131)
  B(272) = RCT(155)*V(40)
  B(273) = RCT(156)*V(131)
  B(274) = RCT(156)*V(21)
  B(275) = RCT(157)*V(131)
  B(276) = RCT(157)*V(59)
  B(277) = RCT(158)*V(131)
  B(278) = RCT(158)*V(58)
  B(279) = RCT(159)*V(131)
  B(280) = RCT(159)*V(41)
  B(281) = RCT(160)*V(131)
  B(282) = RCT(160)*V(22)
  B(283) = RCT(161)*V(131)
  B(284) = RCT(161)*V(57)
  B(285) = RCT(162)*V(131)
  B(286) = RCT(162)*V(56)
  B(287) = RCT(163)*V(131)
  B(288) = RCT(163)*V(42)
  B(289) = RCT(164)*V(131)
  B(290) = RCT(164)*V(23)
  B(291) = RCT(165)*V(131)
  B(292) = RCT(165)*V(55)
  B(293) = RCT(166)*V(131)
  B(294) = RCT(166)*V(54)
  B(295) = RCT(167)*V(131)
  B(296) = RCT(167)*V(43)
  B(297) = RCT(168)*V(131)
  B(298) = RCT(168)*V(24)
  B(299) = RCT(169)*V(131)
  B(300) = RCT(169)*V(53)
  B(301) = RCT(170)*V(131)
  B(302) = RCT(170)*V(52)
  B(303) = RCT(171)*V(131)
  B(304) = RCT(171)*V(44)
  B(305) = RCT(172)*V(131)
  B(306) = RCT(172)*V(25)
  B(307) = RCT(173)*V(131)
  B(308) = RCT(173)*V(51)
  B(309) = RCT(174)*V(131)
  B(310) = RCT(174)*V(50)
  B(311) = RCT(175)*V(131)
  B(312) = RCT(175)*V(45)
  B(313) = RCT(176)*V(131)
  B(314) = RCT(176)*V(26)
  B(315) = RCT(177)*V(131)
  B(316) = RCT(177)*V(49)
  B(317) = RCT(178)*V(131)
  B(318) = RCT(178)*V(47)
  B(319) = RCT(179)*V(131)
  B(320) = RCT(179)*V(46)
  B(321) = RCT(180)*V(131)
  B(322) = RCT(180)*V(27)
  B(323) = RCT(181)*V(131)
  B(324) = RCT(181)*V(63)
  B(325) = RCT(182)*V(131)
  B(326) = RCT(182)*V(28)
  B(327) = RCT(183)*V(131)
  B(328) = RCT(183)*V(84)
  B(329) = RCT(184)*V(131)
  B(330) = RCT(184)*V(86)
  B(331) = RCT(185)*V(131)
  B(332) = RCT(185)*V(65)
  B(333) = RCT(186)*V(131)
  B(334) = RCT(186)*V(29)
  B(335) = RCT(187)*V(131)
  B(336) = RCT(187)*V(83)
  B(337) = RCT(188)*V(131)
  B(338) = RCT(188)*V(82)
  B(339) = RCT(189)*V(131)
  B(340) = RCT(189)*V(66)
  B(341) = RCT(190)*V(131)
  B(342) = RCT(190)*V(30)
  B(343) = RCT(191)*V(131)
  B(344) = RCT(191)*V(81)
  B(345) = RCT(192)*V(131)
  B(346) = RCT(192)*V(80)
  B(347) = RCT(193)*V(131)
  B(348) = RCT(193)*V(67)
  B(349) = RCT(194)*V(131)
  B(350) = RCT(194)*V(31)
  B(351) = RCT(195)*V(131)
  B(352) = RCT(195)*V(79)
  B(353) = RCT(196)*V(131)
  B(354) = RCT(196)*V(78)
  B(355) = RCT(197)*V(131)
  B(356) = RCT(197)*V(68)
  B(357) = RCT(198)*V(131)
  B(358) = RCT(198)*V(32)
  B(359) = RCT(199)*V(131)
  B(360) = RCT(199)*V(77)
  B(361) = RCT(200)*V(131)
  B(362) = RCT(200)*V(76)
  B(363) = RCT(201)*V(131)
  B(364) = RCT(201)*V(69)
  B(365) = RCT(202)*V(131)
  B(366) = RCT(202)*V(33)
  B(367) = RCT(203)*V(131)
  B(368) = RCT(203)*V(75)
  B(369) = RCT(204)*V(131)
  B(370) = RCT(204)*V(74)
  B(371) = RCT(205)*V(131)
  B(372) = RCT(205)*V(70)
  B(373) = RCT(206)*V(131)
  B(374) = RCT(206)*V(34)
  B(375) = RCT(207)*V(131)
  B(376) = RCT(207)*V(73)
  B(377) = RCT(208)*V(131)
  B(378) = RCT(208)*V(72)
  B(379) = RCT(209)*V(131)
  B(380) = RCT(209)*V(71)
  B(381) = RCT(210)*V(131)
  B(382) = RCT(210)*V(35)
  JVS(1) = 0
  JVS(2) = B(85)
  JVS(3) = B(86)
  JVS(4) = 0
  JVS(5) = 0
  JVS(6) = 0
  JVS(7) = 0
  JVS(8) = 0
  JVS(9) = 0
  JVS(10) = 0
  JVS(11) = 0.52*B(107)
  JVS(12) = 0.22*B(111)
  JVS(13) = 0.52*B(108)+0.22*B(112)
  JVS(14) = 0
  JVS(15) = 0.09*B(111)
  JVS(16) = 0.39*B(141)
  JVS(17) = 0.16*B(113)
  JVS(18) = 0.46*B(147)
  JVS(19) = 0.09*B(112)+0.16*B(114)+0.39*B(142)+0.46*B(148)
  JVS(20) = 0.4*B(210)
  JVS(21) = 0.4*B(211)
  JVS(22) = 0
  JVS(23) = 0.039*B(123)
  JVS(24) = 0.039*B(125)
  JVS(25) = 0.039*B(129)+0.039*B(131)
  JVS(26) = 0.039*B(124)+0.039*B(126)+0.039*B(130)
  JVS(27) = 0.039*B(132)
  JVS(28) = 0
  JVS(29) = 0.108*B(123)
  JVS(30) = 0.108*B(125)
  JVS(31) = 0.108*B(129)+0.108*B(131)
  JVS(32) = 0.108*B(124)+0.108*B(126)+0.108*B(130)
  JVS(33) = 0.108*B(132)
  JVS(34) = 0
  JVS(35) = 0.012*B(91)
  JVS(36) = 0.012*B(92)
  JVS(37) = 0
  JVS(38) = 0.008*B(111)+0.008*B(115)+0.008*B(119)
  JVS(39) = 0.008*B(113)+0.008*B(117)+0.008*B(121)
  JVS(40) = 0.008*B(112)+0.008*B(114)
  JVS(41) = 0.008*B(116)+0.008*B(118)
  JVS(42) = 0.008*B(120)+0.008*B(122)
  JVS(43) = 0
  JVS(44) = 0.0064*B(235)+0.022*B(239)
  JVS(45) = 0.022*B(240)
  JVS(46) = 0.0064*B(236)
  JVS(47) = 0
  JVS(48) = 0.055*B(235)+0.19*B(239)
  JVS(49) = 0.19*B(240)
  JVS(50) = 0.055*B(236)
  JVS(51) = 0
  JVS(52) = 0.037*B(241)+0.13*B(245)
  JVS(53) = 0.13*B(246)
  JVS(54) = 0.037*B(242)
  JVS(55) = 0
  JVS(56) = 0.056*B(241)+0.19*B(245)
  JVS(57) = 0.19*B(246)
  JVS(58) = 0.056*B(242)
  JVS(59) = 0
  JVS(60) = B(85)
  JVS(61) = B(218)
  JVS(62) = B(59)
  JVS(63) = B(123)
  JVS(64) = B(125)
  JVS(65) = B(87)
  JVS(66) = B(55)
  JVS(67) = B(89)
  JVS(68) = B(163)
  JVS(69) = B(165)
  JVS(70) = B(51)
  JVS(71) = B(109)
  JVS(72) = B(174)
  JVS(73) = B(129)
  JVS(74) = B(53)
  JVS(75) = B(83)
  JVS(76) = B(91)
  JVS(77) = B(135)
  JVS(78) = B(115)
  JVS(79) = B(139)
  JVS(80) = B(93)
  JVS(81) = B(101)
  JVS(82) = B(117)
  JVS(83) = B(103)
  JVS(84) = B(97)
  JVS(85) = B(145)
  JVS(86) = B(169)
  JVS(87) = B(39)
  JVS(88) = B(47)
  JVS(89) = B(40)+B(43)+B(45)+B(48)+B(49)+B(52)+B(54)+B(56)+B(57)+B(60)+B(84)+B(86)+B(88)+B(90)+B(92)+B(94)+B(98)+B(102)&
              &+B(104)+B(110)+B(116)+B(118)+B(124)+B(126)+B(130)+B(136)+B(140)+B(146)+B(164)+B(166)+B(167)+B(170)+B(175)&
              &+B(219)
  JVS(90) = B(58)
  JVS(91) = B(46)
  JVS(92) = B(168)
  JVS(93) = B(50)
  JVS(94) = B(247)
  JVS(95) = B(249)
  JVS(96) = B(255)
  JVS(97) = B(257)
  JVS(98) = 0
  JVS(99) = B(265)
  JVS(100) = B(273)
  JVS(101) = B(281)
  JVS(102) = B(289)
  JVS(103) = B(297)
  JVS(104) = B(305)
  JVS(105) = B(313)
  JVS(106) = B(321)
  JVS(107) = B(325)
  JVS(108) = B(333)
  JVS(109) = B(341)
  JVS(110) = B(349)
  JVS(111) = B(357)
  JVS(112) = B(365)
  JVS(113) = B(373)
  JVS(114) = B(381)
  JVS(115) = B(263)
  JVS(116) = B(251)
  JVS(117) = B(271)
  JVS(118) = B(279)
  JVS(119) = B(287)
  JVS(120) = B(295)
  JVS(121) = B(303)
  JVS(122) = B(311)
  JVS(123) = B(319)
  JVS(124) = B(317)
  JVS(125) = B(315)
  JVS(126) = B(309)
  JVS(127) = B(307)
  JVS(128) = B(301)
  JVS(129) = B(299)
  JVS(130) = B(293)
  JVS(131) = B(291)
  JVS(132) = B(285)
  JVS(133) = B(283)
  JVS(134) = B(277)
  JVS(135) = B(275)
  JVS(136) = B(267)
  JVS(137) = B(253)
  JVS(138) = B(269)
  JVS(139) = B(323)
  JVS(140) = B(259)
  JVS(141) = B(331)
  JVS(142) = B(339)
  JVS(143) = B(347)
  JVS(144) = B(355)
  JVS(145) = B(363)
  JVS(146) = B(371)
  JVS(147) = B(379)
  JVS(148) = B(377)
  JVS(149) = B(375)
  JVS(150) = B(369)
  JVS(151) = B(367)
  JVS(152) = B(361)
  JVS(153) = B(359)
  JVS(154) = B(353)
  JVS(155) = B(351)
  JVS(156) = B(345)
  JVS(157) = B(343)
  JVS(158) = B(337)
  JVS(159) = B(335)
  JVS(160) = B(327)
  JVS(161) = B(261)
  JVS(162) = B(329)
  JVS(163) = B(235)
  JVS(164) = B(241)
  JVS(165) = B(236)+B(242)+B(248)+B(250)+B(252)+B(254)+B(256)+B(258)+B(260)+B(262)+B(264)+B(266)+B(268)+B(270)+B(272)&
               &+B(274)+B(276)+B(278)+B(280)+B(282)+B(284)+B(286)+B(288)+B(290)+B(292)+B(294)+B(296)+B(298)+B(300)+B(302)&
               &+B(304)+B(306)+B(308)+B(310)+B(312)+B(314)+B(316)+B(318)+B(320)+B(322)+B(324)+B(326)+B(328)+B(330)+B(332)&
               &+B(334)+B(336)+B(338)+B(340)+B(342)+B(344)+B(346)+B(348)+B(350)+B(352)+B(354)+B(356)+B(358)+B(360)+B(362)&
               &+B(364)+B(366)+B(368)+B(370)+B(372)+B(374)+B(376)+B(378)+B(380)+B(382)
  JVS(166) = -B(265)
  JVS(167) = -B(266)
  JVS(168) = -B(273)
  JVS(169) = -B(274)
  JVS(170) = -B(281)
  JVS(171) = -B(282)
  JVS(172) = -B(289)
  JVS(173) = -B(290)
  JVS(174) = -B(297)
  JVS(175) = -B(298)
  JVS(176) = -B(305)
  JVS(177) = -B(306)
  JVS(178) = -B(313)
  JVS(179) = -B(314)
  JVS(180) = -B(321)
  JVS(181) = -B(322)
  JVS(182) = -B(325)
  JVS(183) = -B(326)
  JVS(184) = -B(333)
  JVS(185) = -B(334)
  JVS(186) = -B(341)
  JVS(187) = -B(342)
  JVS(188) = -B(349)
  JVS(189) = -B(350)
  JVS(190) = -B(357)
  JVS(191) = -B(358)
  JVS(192) = -B(365)
  JVS(193) = -B(366)
  JVS(194) = -B(373)
  JVS(195) = -B(374)
  JVS(196) = -B(381)
  JVS(197) = -B(382)
  JVS(198) = -B(21)-B(23)
  JVS(199) = B(8)
  JVS(200) = -B(85)
  JVS(201) = -B(86)
  JVS(202) = -B(263)
  JVS(203) = -B(264)
  JVS(204) = B(263)
  JVS(205) = 0
  JVS(206) = B(267)
  JVS(207) = B(264)+B(268)
  JVS(208) = -B(271)
  JVS(209) = -B(272)
  JVS(210) = -B(279)
  JVS(211) = -B(280)
  JVS(212) = -B(287)
  JVS(213) = -B(288)
  JVS(214) = -B(295)
  JVS(215) = -B(296)
  JVS(216) = -B(303)
  JVS(217) = -B(304)
  JVS(218) = -B(311)
  JVS(219) = -B(312)
  JVS(220) = -B(319)
  JVS(221) = -B(320)
  JVS(222) = B(321)
  JVS(223) = 0.075*B(319)
  JVS(224) = -B(317)
  JVS(225) = -B(318)+0.075*B(320)+B(322)
  JVS(226) = -B(218)
  JVS(227) = -B(219)
  JVS(228) = B(319)
  JVS(229) = -B(315)
  JVS(230) = -B(316)+B(320)
  JVS(231) = B(313)
  JVS(232) = 0.075*B(311)
  JVS(233) = B(317)
  JVS(234) = 0.075*B(315)
  JVS(235) = -B(309)
  JVS(236) = -B(310)+0.075*B(312)+B(314)+0.075*B(316)+B(318)
  JVS(237) = B(311)
  JVS(238) = B(315)
  JVS(239) = -B(307)
  JVS(240) = -B(308)+B(312)+B(316)
  JVS(241) = B(305)
  JVS(242) = 0.075*B(303)
  JVS(243) = B(309)
  JVS(244) = 0.075*B(307)
  JVS(245) = -B(301)
  JVS(246) = -B(302)+0.075*B(304)+B(306)+0.075*B(308)+B(310)
  JVS(247) = B(303)
  JVS(248) = B(307)
  JVS(249) = -B(299)
  JVS(250) = -B(300)+B(304)+B(308)
  JVS(251) = B(297)
  JVS(252) = 0.075*B(295)
  JVS(253) = B(301)
  JVS(254) = 0.075*B(299)
  JVS(255) = -B(293)
  JVS(256) = -B(294)+0.075*B(296)+B(298)+0.075*B(300)+B(302)
  JVS(257) = B(295)
  JVS(258) = B(299)
  JVS(259) = -B(291)
  JVS(260) = -B(292)+B(296)+B(300)
  JVS(261) = B(289)
  JVS(262) = 0.075*B(287)
  JVS(263) = B(293)
  JVS(264) = 0.075*B(291)
  JVS(265) = -B(285)
  JVS(266) = -B(286)+0.075*B(288)+B(290)+0.075*B(292)+B(294)
  JVS(267) = B(287)
  JVS(268) = B(291)
  JVS(269) = -B(283)
  JVS(270) = -B(284)+B(288)+B(292)
  JVS(271) = B(281)
  JVS(272) = 0.075*B(279)
  JVS(273) = B(285)
  JVS(274) = 0.075*B(283)
  JVS(275) = -B(277)
  JVS(276) = -B(278)+0.075*B(280)+B(282)+0.075*B(284)+B(286)
  JVS(277) = B(279)
  JVS(278) = B(283)
  JVS(279) = -B(275)
  JVS(280) = -B(276)+B(280)+B(284)
  JVS(281) = B(271)
  JVS(282) = B(275)
  JVS(283) = -B(267)
  JVS(284) = -B(268)+B(272)+B(276)
  JVS(285) = B(265)
  JVS(286) = 0.075*B(263)
  JVS(287) = 0.075*B(267)
  JVS(288) = 0
  JVS(289) = B(269)
  JVS(290) = 0.075*B(264)+B(266)+0.075*B(268)+B(270)
  JVS(291) = B(273)
  JVS(292) = 0.075*B(271)
  JVS(293) = B(277)
  JVS(294) = 0.075*B(275)
  JVS(295) = -B(269)
  JVS(296) = -B(270)+0.075*B(272)+B(274)+0.075*B(276)+B(278)
  JVS(297) = -B(323)
  JVS(298) = -B(324)
  JVS(299) = B(323)
  JVS(300) = 0
  JVS(301) = B(327)
  JVS(302) = B(324)+B(328)
  JVS(303) = -B(331)
  JVS(304) = -B(332)
  JVS(305) = -B(339)
  JVS(306) = -B(340)
  JVS(307) = -B(347)
  JVS(308) = -B(348)
  JVS(309) = -B(355)
  JVS(310) = -B(356)
  JVS(311) = -B(363)
  JVS(312) = -B(364)
  JVS(313) = -B(371)
  JVS(314) = -B(372)
  JVS(315) = -B(379)
  JVS(316) = -B(380)
  JVS(317) = B(381)
  JVS(318) = 0.075*B(379)
  JVS(319) = -B(377)
  JVS(320) = -B(378)+0.075*B(380)+B(382)
  JVS(321) = B(379)
  JVS(322) = -B(375)
  JVS(323) = -B(376)+B(380)
  JVS(324) = B(373)
  JVS(325) = 0.075*B(371)
  JVS(326) = B(377)
  JVS(327) = 0.075*B(375)
  JVS(328) = -B(369)
  JVS(329) = -B(370)+0.075*B(372)+B(374)+0.075*B(376)+B(378)
  JVS(330) = B(371)
  JVS(331) = B(375)
  JVS(332) = -B(367)
  JVS(333) = -B(368)+B(372)+B(376)
  JVS(334) = B(365)
  JVS(335) = 0.075*B(363)
  JVS(336) = B(369)
  JVS(337) = 0.075*B(367)
  JVS(338) = -B(361)
  JVS(339) = -B(362)+0.075*B(364)+B(366)+0.075*B(368)+B(370)
  JVS(340) = B(363)
  JVS(341) = B(367)
  JVS(342) = -B(359)
  JVS(343) = -B(360)+B(364)+B(368)
  JVS(344) = B(357)
  JVS(345) = 0.075*B(355)
  JVS(346) = B(361)
  JVS(347) = 0.075*B(359)
  JVS(348) = -B(353)
  JVS(349) = -B(354)+0.075*B(356)+B(358)+0.075*B(360)+B(362)
  JVS(350) = B(355)
  JVS(351) = B(359)
  JVS(352) = -B(351)
  JVS(353) = -B(352)+B(356)+B(360)
  JVS(354) = B(349)
  JVS(355) = 0.075*B(347)
  JVS(356) = B(353)
  JVS(357) = 0.075*B(351)
  JVS(358) = -B(345)
  JVS(359) = -B(346)+0.075*B(348)+B(350)+0.075*B(352)+B(354)
  JVS(360) = B(347)
  JVS(361) = B(351)
  JVS(362) = -B(343)
  JVS(363) = -B(344)+B(348)+B(352)
  JVS(364) = B(341)
  JVS(365) = 0.075*B(339)
  JVS(366) = B(345)
  JVS(367) = 0.075*B(343)
  JVS(368) = -B(337)
  JVS(369) = -B(338)+0.075*B(340)+B(342)+0.075*B(344)+B(346)
  JVS(370) = B(339)
  JVS(371) = B(343)
  JVS(372) = -B(335)
  JVS(373) = -B(336)+B(340)+B(344)
  JVS(374) = B(331)
  JVS(375) = B(335)
  JVS(376) = -B(327)
  JVS(377) = -B(328)+B(332)+B(336)
  JVS(378) = B(325)
  JVS(379) = 0.075*B(323)
  JVS(380) = 0.075*B(327)
  JVS(381) = 0
  JVS(382) = B(329)
  JVS(383) = 0.075*B(324)+B(326)+0.075*B(328)+B(330)
  JVS(384) = B(333)
  JVS(385) = 0.075*B(331)
  JVS(386) = B(337)
  JVS(387) = 0.075*B(335)
  JVS(388) = -B(329)
  JVS(389) = -B(330)+0.075*B(332)+B(334)+0.075*B(336)+B(338)
  JVS(390) = -B(173)
  JVS(391) = B(171)
  JVS(392) = B(172)
  JVS(393) = -B(9)-B(59)
  JVS(394) = -B(60)
  JVS(395) = B(61)+B(62)
  JVS(396) = -B(123)
  JVS(397) = -B(124)
  JVS(398) = -B(6)-B(80)-B(82)
  JVS(399) = B(75)
  JVS(400) = B(76)
  JVS(401) = -B(125)
  JVS(402) = -B(126)
  JVS(403) = -B(133)
  JVS(404) = 0.4*B(129)+B(131)
  JVS(405) = -B(134)
  JVS(406) = 0.4*B(130)
  JVS(407) = B(132)
  JVS(408) = -B(87)
  JVS(409) = 0.06*B(111)
  JVS(410) = 0.08*B(113)
  JVS(411) = 0.06*B(112)+0.08*B(114)
  JVS(412) = -B(88)
  JVS(413) = -B(235)-B(237)-B(239)
  JVS(414) = -B(240)
  JVS(415) = -B(236)
  JVS(416) = -B(238)
  JVS(417) = -B(241)-B(243)-B(245)
  JVS(418) = -B(246)
  JVS(419) = -B(242)
  JVS(420) = -B(244)
  JVS(421) = -B(5)-B(55)-B(70)
  JVS(422) = B(66)
  JVS(423) = -B(56)
  JVS(424) = B(67)
  JVS(425) = 0.8*B(123)
  JVS(426) = 0.45*B(125)
  JVS(427) = -B(127)
  JVS(428) = 0.8*B(124)+0.45*B(126)
  JVS(429) = -B(128)
  JVS(430) = -B(89)
  JVS(431) = 0.01*B(111)
  JVS(432) = 0.01*B(113)
  JVS(433) = 0.2*B(226)
  JVS(434) = 0.01*B(112)+0.01*B(114)
  JVS(435) = -B(90)
  JVS(436) = -B(220)
  JVS(437) = -B(221)
  JVS(438) = 1.06*B(111)+B(115)
  JVS(439) = 2.26*B(113)+2.23*B(117)
  JVS(440) = B(186)+B(200)+B(230)
  JVS(441) = 1.98*B(19)
  JVS(442) = 1.06*B(112)+2.26*B(114)
  JVS(443) = 1.68*B(180)+1.98*B(194)+1.25*B(227)
  JVS(444) = B(116)+2.23*B(118)+0.42*B(167)
  JVS(445) = 1.68*B(181)+B(187)
  JVS(446) = 1.98*B(18)+0.42*B(168)
  JVS(447) = 1.98*B(195)+B(201)
  JVS(448) = -B(12)-B(163)
  JVS(449) = B(204)
  JVS(450) = -B(164)
  JVS(451) = B(205)
  JVS(452) = -B(13)-B(165)
  JVS(453) = B(206)
  JVS(454) = -B(166)
  JVS(455) = B(207)
  JVS(456) = -B(3)-B(51)
  JVS(457) = B(68)
  JVS(458) = B(45)-B(52)
  JVS(459) = B(69)
  JVS(460) = B(46)
  JVS(461) = -B(107)-B(109)
  JVS(462) = -B(108)
  JVS(463) = -B(110)
  JVS(464) = -B(174)
  JVS(465) = 0.03*B(111)
  JVS(466) = 0.04*B(113)
  JVS(467) = 0.34*B(225)
  JVS(468) = 0.03*B(112)+0.04*B(114)
  JVS(469) = -B(175)
  JVS(470) = 0.12*B(123)
  JVS(471) = 0.05*B(125)
  JVS(472) = -B(129)-B(131)
  JVS(473) = 0.12*B(124)+0.05*B(126)-B(130)
  JVS(474) = -B(132)
  JVS(475) = B(21)
  JVS(476) = -B(25)-B(27)-B(29)-B(31)-B(33)
  JVS(477) = B(7)-B(28)
  JVS(478) = B(1)-B(30)-B(32)
  JVS(479) = -B(34)
  JVS(480) = 0.89*B(2)
  JVS(481) = 2*B(80)
  JVS(482) = B(131)
  JVS(483) = -B(4)-B(53)
  JVS(484) = B(95)
  JVS(485) = B(105)
  JVS(486) = B(99)
  JVS(487) = 0.07*B(149)
  JVS(488) = B(47)
  JVS(489) = B(48)-B(54)
  JVS(490) = 0.3*B(78)
  JVS(491) = 0.3*B(79)+B(96)+B(100)+B(106)+B(132)+0.07*B(150)
  JVS(492) = 0.24*B(107)
  JVS(493) = -B(83)
  JVS(494) = B(17)+2*B(135)+0.69*B(137)
  JVS(495) = 0.59*B(155)
  JVS(496) = 0.31*B(111)
  JVS(497) = 0.07*B(141)
  JVS(498) = B(10)+B(11)+B(93)+B(95)
  JVS(499) = 0.3*B(113)
  JVS(500) = B(16)+B(105)
  JVS(501) = B(14)
  JVS(502) = 0.33*B(20)+0.16*B(147)+0.64*B(149)
  JVS(503) = 0.24*B(108)+0.31*B(112)+0.3*B(114)+0.69*B(138)+0.07*B(142)+0.16*B(148)
  JVS(504) = -B(84)+B(94)+2*B(136)
  JVS(505) = 0.59*B(156)
  JVS(506) = B(96)+B(106)+0.64*B(150)
  JVS(507) = 1.1*B(125)
  JVS(508) = -B(220)
  JVS(509) = -B(91)-B(221)
  JVS(510) = 1.6*B(153)+2*B(159)+2*B(232)
  JVS(511) = 0.18*B(151)
  JVS(512) = 0
  JVS(513) = 0
  JVS(514) = 0
  JVS(515) = 1.86*B(149)
  JVS(516) = 0
  JVS(517) = 0
  JVS(518) = 0
  JVS(519) = -B(92)+1.1*B(126)
  JVS(520) = 2*B(160)
  JVS(521) = 0.18*B(152)+1.6*B(154)
  JVS(522) = 0
  JVS(523) = 1.86*B(150)
  JVS(524) = 0.95*B(127)
  JVS(525) = 0.3*B(129)
  JVS(526) = -B(17)-B(135)-B(137)
  JVS(527) = -B(138)
  JVS(528) = 0.3*B(130)-B(136)
  JVS(529) = 0.95*B(128)
  JVS(530) = 0
  JVS(531) = -B(153)-B(159)-B(232)
  JVS(532) = B(143)
  JVS(533) = -B(160)
  JVS(534) = -B(154)
  JVS(535) = B(144)
  JVS(536) = -B(155)-B(161)-B(233)
  JVS(537) = 0.5*B(145)
  JVS(538) = 0.5*B(146)
  JVS(539) = -B(162)
  JVS(540) = -B(156)
  JVS(541) = -B(151)-B(157)-B(231)
  JVS(542) = B(139)
  JVS(543) = B(140)
  JVS(544) = -B(158)
  JVS(545) = -B(152)
  JVS(546) = -B(111)-B(115)-B(119)
  JVS(547) = -B(112)
  JVS(548) = -B(116)
  JVS(549) = -B(120)
  JVS(550) = -B(139)-B(141)-B(143)
  JVS(551) = -B(142)
  JVS(552) = -B(140)
  JVS(553) = -B(144)
  JVS(554) = B(12)+0.3*B(163)
  JVS(555) = B(107)+1.56*B(109)
  JVS(556) = B(174)
  JVS(557) = B(135)+0.7*B(137)
  JVS(558) = 0.25*B(155)
  JVS(559) = 0.63*B(151)
  JVS(560) = 0.57*B(111)+B(115)
  JVS(561) = 0.6*B(141)
  JVS(562) = -B(10)-B(11)-B(93)-B(95)
  JVS(563) = 0
  JVS(564) = 0.5*B(186)+0.5*B(200)+0.5*B(230)
  JVS(565) = B(176)+B(190)+0.66*B(225)
  JVS(566) = 0.2*B(20)+0.15*B(147)+0.28*B(149)
  JVS(567) = B(108)+0.57*B(112)+0.7*B(138)+0.6*B(142)+0.15*B(148)
  JVS(568) = -B(94)+1.56*B(110)+B(116)+B(136)+0.3*B(164)+B(175)
  JVS(569) = 0
  JVS(570) = 0.63*B(152)+0.25*B(156)+B(177)+B(184)+0.5*B(187)
  JVS(571) = B(185)+B(198)+0.7*B(229)
  JVS(572) = -B(96)+0.28*B(150)+B(191)+B(199)+0.5*B(201)
  JVS(573) = 0.08*B(123)
  JVS(574) = 0.5*B(125)
  JVS(575) = B(109)
  JVS(576) = 0.6*B(129)
  JVS(577) = B(135)+0.03*B(137)
  JVS(578) = B(115)
  JVS(579) = 0.08*B(139)+0.2*B(141)
  JVS(580) = -B(188)-B(202)-B(216)-B(234)
  JVS(581) = B(117)
  JVS(582) = B(103)
  JVS(583) = 0.2*B(145)+0.07*B(147)+0.93*B(149)
  JVS(584) = 0.41*B(19)
  JVS(585) = 0.03*B(138)+0.2*B(142)+0.07*B(148)
  JVS(586) = 0.34*B(180)+0.4*B(194)+0.24*B(227)
  JVS(587) = B(104)+B(110)+B(116)+B(118)+0.08*B(124)+0.5*B(126)+0.6*B(130)+B(136)+0.08*B(140)+0.2*B(146)
  JVS(588) = -B(217)
  JVS(589) = 0.34*B(181)-B(189)
  JVS(590) = 0.4*B(18)
  JVS(591) = 0.93*B(150)+0.4*B(195)-B(203)
  JVS(592) = 0.63*B(155)+0.5*B(233)
  JVS(593) = -B(15)-B(101)
  JVS(594) = 0.07*B(113)+0.23*B(117)
  JVS(595) = 0.03*B(20)+0.09*B(147)
  JVS(596) = 0.74*B(19)
  JVS(597) = 0.07*B(114)+0.09*B(148)
  JVS(598) = 0.62*B(180)+0.74*B(194)+0.57*B(227)
  JVS(599) = -B(102)+0.23*B(118)
  JVS(600) = 0
  JVS(601) = 0.63*B(156)+0.62*B(181)
  JVS(602) = 0.15*B(229)
  JVS(603) = 0.74*B(18)
  JVS(604) = 0.74*B(195)
  JVS(605) = -B(113)-B(117)-B(121)
  JVS(606) = -B(114)
  JVS(607) = -B(118)
  JVS(608) = -B(122)
  JVS(609) = B(119)
  JVS(610) = B(121)
  JVS(611) = -B(186)-B(200)-B(214)-B(230)
  JVS(612) = B(169)
  JVS(613) = 0
  JVS(614) = B(170)
  JVS(615) = -B(215)
  JVS(616) = -B(187)
  JVS(617) = B(120)+B(122)-B(201)
  JVS(618) = 0.8*B(125)
  JVS(619) = 0.2*B(137)
  JVS(620) = 0.34*B(155)
  JVS(621) = 0.04*B(111)
  JVS(622) = 0.07*B(113)
  JVS(623) = -B(16)-B(103)-B(105)
  JVS(624) = 0.85*B(147)
  JVS(625) = 0.04*B(112)+0.07*B(114)+0.2*B(138)+0.85*B(148)
  JVS(626) = -B(104)+0.8*B(126)+0.19*B(167)
  JVS(627) = 0
  JVS(628) = 0.34*B(156)
  JVS(629) = 0.15*B(229)
  JVS(630) = 0.19*B(168)
  JVS(631) = -B(106)
  JVS(632) = B(89)
  JVS(633) = 0.7*B(165)
  JVS(634) = 0.06*B(111)
  JVS(635) = 0.05*B(113)
  JVS(636) = -B(178)-B(192)-B(206)-B(226)
  JVS(637) = 0.1*B(19)
  JVS(638) = 0.06*B(112)+0.05*B(114)
  JVS(639) = 0.08*B(180)+0.1*B(194)+0.06*B(227)
  JVS(640) = B(90)+0.7*B(166)
  JVS(641) = -B(207)
  JVS(642) = -B(179)+0.08*B(181)
  JVS(643) = 0.1*B(18)
  JVS(644) = -B(193)+0.1*B(195)
  JVS(645) = B(218)
  JVS(646) = B(13)+0.3*B(165)
  JVS(647) = 0.22*B(109)
  JVS(648) = 0.03*B(137)
  JVS(649) = 0.8*B(153)+B(232)
  JVS(650) = 0.55*B(155)+0.5*B(233)
  JVS(651) = 0.47*B(111)+B(115)
  JVS(652) = 0.15*B(141)
  JVS(653) = 1.03*B(113)+1.77*B(117)
  JVS(654) = 0.5*B(186)+0.5*B(200)+0.5*B(230)
  JVS(655) = B(178)+B(192)+0.8*B(226)
  JVS(656) = -B(14)-B(97)-B(99)
  JVS(657) = 0.07*B(20)+0.02*B(147)+0.28*B(149)
  JVS(658) = 0.3*B(19)
  JVS(659) = 0.47*B(112)+1.03*B(114)+0.03*B(138)+0.15*B(142)+0.02*B(148)
  JVS(660) = 0.25*B(180)+0.3*B(194)+0.21*B(227)
  JVS(661) = -B(98)+0.22*B(110)+B(116)+1.77*B(118)+0.3*B(166)+0.04*B(167)+B(219)
  JVS(662) = 0
  JVS(663) = 0.8*B(154)+0.55*B(156)+B(179)+0.25*B(181)+0.5*B(187)
  JVS(664) = 0.3*B(18)+0.04*B(168)
  JVS(665) = -B(100)+0.28*B(150)+B(193)+0.3*B(195)+0.5*B(201)
  JVS(666) = B(87)
  JVS(667) = 0.7*B(163)
  JVS(668) = 0.07*B(111)
  JVS(669) = B(15)
  JVS(670) = 0.1*B(113)
  JVS(671) = B(14)
  JVS(672) = -B(176)-B(190)-B(204)-B(225)
  JVS(673) = 0.7*B(20)+0.05*B(147)
  JVS(674) = 0
  JVS(675) = 0.07*B(112)+0.1*B(114)+0.05*B(148)
  JVS(676) = 0
  JVS(677) = B(182)+B(196)+B(228)
  JVS(678) = B(88)+0.7*B(164)
  JVS(679) = -B(205)
  JVS(680) = -B(177)+B(183)
  JVS(681) = 0
  JVS(682) = 0
  JVS(683) = -B(191)+B(197)
  JVS(684) = 0.2*B(153)
  JVS(685) = 0.91*B(151)+B(231)
  JVS(686) = 0.65*B(141)
  JVS(687) = -B(20)-B(145)-B(147)-B(149)
  JVS(688) = 0.65*B(142)-B(148)
  JVS(689) = -B(146)
  JVS(690) = 0
  JVS(691) = 0.91*B(152)+0.2*B(154)
  JVS(692) = -B(150)
  JVS(693) = B(133)
  JVS(694) = 0.05*B(127)
  JVS(695) = 0
  JVS(696) = 0.8*B(153)+B(159)+B(232)
  JVS(697) = 0.09*B(151)
  JVS(698) = 0
  JVS(699) = 0.5*B(186)+0.5*B(200)+B(214)+0.5*B(230)
  JVS(700) = 0.93*B(149)
  JVS(701) = -B(19)-B(169)
  JVS(702) = 0
  JVS(703) = 0.16*B(180)
  JVS(704) = B(134)
  JVS(705) = -B(170)
  JVS(706) = B(160)+B(215)
  JVS(707) = 0.05*B(128)+0.09*B(152)+0.8*B(154)+0.16*B(181)+0.5*B(187)
  JVS(708) = 0.93*B(150)+0.5*B(201)
  JVS(709) = -B(239)
  JVS(710) = -B(245)
  JVS(711) = -B(107)
  JVS(712) = B(25)-B(27)
  JVS(713) = -B(137)
  JVS(714) = -B(111)
  JVS(715) = -B(141)
  JVS(716) = -B(113)
  JVS(717) = -B(147)
  JVS(718) = -B(7)-B(8)-B(28)-B(35)-B(37)-B(39)-B(41)-B(108)-B(112)-B(114)-B(138)-B(142)-B(148)-B(240)-B(246)
  JVS(719) = -B(38)
  JVS(720) = 0.4*B(210)
  JVS(721) = -B(40)
  JVS(722) = -B(42)+0.4*B(211)
  JVS(723) = -B(36)
  JVS(724) = 0
  JVS(725) = B(91)
  JVS(726) = 0
  JVS(727) = 0
  JVS(728) = 0.03*B(111)
  JVS(729) = 0
  JVS(730) = 0.09*B(113)
  JVS(731) = 0
  JVS(732) = 0
  JVS(733) = 0
  JVS(734) = 0.03*B(112)+0.09*B(114)
  JVS(735) = -B(180)-B(194)-B(208)-B(227)
  JVS(736) = 0
  JVS(737) = 0
  JVS(738) = B(92)+0.77*B(167)
  JVS(739) = -B(209)
  JVS(740) = -B(181)
  JVS(741) = 0.77*B(168)
  JVS(742) = -B(195)
  JVS(743) = B(173)
  JVS(744) = B(6)+B(82)
  JVS(745) = -B(133)
  JVS(746) = B(5)+B(55)+B(70)
  JVS(747) = 0.95*B(127)
  JVS(748) = B(51)
  JVS(749) = 0
  JVS(750) = -B(29)-B(31)+B(33)
  JVS(751) = B(4)
  JVS(752) = 1.2*B(153)
  JVS(753) = B(155)
  JVS(754) = 0.91*B(151)
  JVS(755) = 0
  JVS(756) = 0
  JVS(757) = B(188)+B(202)
  JVS(758) = 0
  JVS(759) = 1.5*B(186)+1.5*B(200)+0.5*B(230)
  JVS(760) = 0
  JVS(761) = B(178)+B(192)
  JVS(762) = 0
  JVS(763) = B(176)+B(190)
  JVS(764) = 0
  JVS(765) = B(19)
  JVS(766) = B(35)-B(37)
  JVS(767) = 0.84*B(180)+B(194)
  JVS(768) = -B(1)-B(30)-B(32)-B(38)-B(47)-B(66)-B(68)-B(75)-B(134)-B(171)
  JVS(769) = -B(172)+B(182)+B(196)
  JVS(770) = -B(48)+B(49)+B(52)+B(56)
  JVS(771) = B(64)-B(67)-B(69)+0.7*B(78)
  JVS(772) = B(34)+B(36)+B(65)+2*B(71)+0.95*B(128)+0.91*B(152)+1.2*B(154)+B(156)+B(177)+B(179)+0.84*B(181)+B(183)+B(184)&
               &+1.5*B(187)+B(189)
  JVS(773) = B(185)+B(198)
  JVS(774) = 0
  JVS(775) = 0.89*B(2)+B(50)+2*B(72)-B(76)+2*B(77)+0.7*B(79)+B(191)+B(193)+B(195)+B(197)+B(199)+1.5*B(201)+B(203)
  JVS(776) = B(173)
  JVS(777) = B(17)+B(135)+0.62*B(137)
  JVS(778) = 0.13*B(111)
  JVS(779) = 0.2*B(141)
  JVS(780) = B(15)
  JVS(781) = 0.19*B(113)
  JVS(782) = B(16)+B(103)+B(105)
  JVS(783) = B(97)+B(99)
  JVS(784) = 0.97*B(20)+0.5*B(145)+0.11*B(147)+0.07*B(149)
  JVS(785) = 0
  JVS(786) = 0.13*B(112)+0.19*B(114)+0.62*B(138)+0.2*B(142)+0.11*B(148)
  JVS(787) = 0
  JVS(788) = -B(171)
  JVS(789) = -B(172)-B(182)-B(196)-B(210)-B(228)
  JVS(790) = B(98)+B(104)+B(136)+0.5*B(146)
  JVS(791) = -B(211)
  JVS(792) = -B(183)+B(184)
  JVS(793) = B(185)+B(198)+0.7*B(229)
  JVS(794) = 0
  JVS(795) = B(100)+B(106)+0.07*B(150)-B(197)+B(199)
  JVS(796) = -B(247)
  JVS(797) = -B(249)
  JVS(798) = -B(255)
  JVS(799) = -B(257)
  JVS(800) = -B(265)
  JVS(801) = -B(273)
  JVS(802) = -B(281)
  JVS(803) = -B(289)
  JVS(804) = -B(297)
  JVS(805) = -B(305)
  JVS(806) = -B(313)
  JVS(807) = -B(321)
  JVS(808) = -B(325)
  JVS(809) = -B(333)
  JVS(810) = -B(341)
  JVS(811) = -B(349)
  JVS(812) = -B(357)
  JVS(813) = -B(365)
  JVS(814) = -B(373)
  JVS(815) = -B(381)
  JVS(816) = 2*B(23)
  JVS(817) = -B(85)
  JVS(818) = -B(263)
  JVS(819) = -B(251)
  JVS(820) = -B(271)
  JVS(821) = -B(279)
  JVS(822) = -B(287)
  JVS(823) = -B(295)
  JVS(824) = -B(303)
  JVS(825) = -B(311)
  JVS(826) = -B(319)
  JVS(827) = -B(317)
  JVS(828) = -B(218)
  JVS(829) = -B(315)
  JVS(830) = -B(309)
  JVS(831) = -B(307)
  JVS(832) = -B(301)
  JVS(833) = -B(299)
  JVS(834) = -B(293)
  JVS(835) = -B(291)
  JVS(836) = -B(285)
  JVS(837) = -B(283)
  JVS(838) = -B(277)
  JVS(839) = -B(275)
  JVS(840) = -B(267)
  JVS(841) = -B(253)
  JVS(842) = -B(269)
  JVS(843) = -B(323)
  JVS(844) = -B(259)
  JVS(845) = -B(331)
  JVS(846) = -B(339)
  JVS(847) = -B(347)
  JVS(848) = -B(355)
  JVS(849) = -B(363)
  JVS(850) = -B(371)
  JVS(851) = -B(379)
  JVS(852) = -B(377)
  JVS(853) = -B(375)
  JVS(854) = -B(369)
  JVS(855) = -B(367)
  JVS(856) = -B(361)
  JVS(857) = -B(359)
  JVS(858) = -B(353)
  JVS(859) = -B(351)
  JVS(860) = -B(345)
  JVS(861) = -B(343)
  JVS(862) = -B(337)
  JVS(863) = -B(335)
  JVS(864) = -B(327)
  JVS(865) = -B(261)
  JVS(866) = -B(329)
  JVS(867) = 2*B(9)-B(59)
  JVS(868) = -B(123)
  JVS(869) = -B(125)
  JVS(870) = -B(87)
  JVS(871) = -B(235)
  JVS(872) = -B(241)
  JVS(873) = -B(55)
  JVS(874) = -B(89)
  JVS(875) = B(12)-0.7*B(163)
  JVS(876) = B(13)-0.7*B(165)
  JVS(877) = B(3)-B(51)
  JVS(878) = 0.12*B(107)-B(109)
  JVS(879) = -B(174)
  JVS(880) = -B(129)
  JVS(881) = B(4)-B(53)
  JVS(882) = -B(83)
  JVS(883) = -B(91)
  JVS(884) = -B(135)+0.08*B(137)
  JVS(885) = 0
  JVS(886) = 0
  JVS(887) = 0
  JVS(888) = 0.33*B(111)-B(115)
  JVS(889) = -B(139)+0.27*B(141)
  JVS(890) = -B(93)
  JVS(891) = -B(101)
  JVS(892) = 0.6*B(113)-B(117)
  JVS(893) = 0
  JVS(894) = -B(103)
  JVS(895) = 0
  JVS(896) = -B(97)
  JVS(897) = 0
  JVS(898) = -B(145)+0.27*B(147)
  JVS(899) = -B(169)
  JVS(900) = -B(39)+B(41)+0.12*B(108)+0.33*B(112)+0.6*B(114)+0.08*B(138)+0.27*B(142)+0.27*B(148)
  JVS(901) = 0
  JVS(902) = -B(47)
  JVS(903) = 0
  JVS(904) = -B(40)-B(43)-B(45)-B(48)-B(49)-B(52)-B(54)-B(56)-B(57)-B(60)-B(84)-B(86)-B(88)-B(90)-B(92)-B(94)-B(98)&
               &-B(102)-B(104)-B(110)-B(116)-B(118)-B(124)-B(126)-B(130)-B(136)-B(140)-B(146)-0.7*B(164)-0.7*B(166)-0.77&
               &*B(167)-B(170)-B(175)-B(219)-B(236)-B(242)-B(248)-B(250)-B(252)-B(254)-B(256)-B(258)-B(260)-B(262)-B(264)&
               &-B(266)-B(268)-B(270)-B(272)-B(274)-B(276)-B(278)-B(280)-B(282)-B(284)-B(286)-B(288)-B(290)-B(292)-B(294)&
               &-B(296)-B(298)-B(300)-B(302)-B(304)-B(306)-B(308)-B(310)-B(312)-B(314)-B(316)-B(318)-B(320)-B(322)-B(324)&
               &-B(326)-B(328)-B(330)-B(332)-B(334)-B(336)-B(338)-B(340)-B(342)-B(344)-B(346)-B(348)-B(350)-B(352)-B(354)&
               &-B(356)-B(358)-B(360)-B(362)-B(364)-B(366)-B(368)-B(370)-B(372)-B(374)-B(376)-B(378)-B(380)-B(382)
  JVS(905) = B(42)-B(58)+B(64)+0.7*B(78)
  JVS(906) = -B(46)+B(65)
  JVS(907) = 0
  JVS(908) = B(18)-0.77*B(168)
  JVS(909) = -B(50)+0.7*B(79)
  JVS(910) = B(85)
  JVS(911) = B(218)
  JVS(912) = B(59)
  JVS(913) = 0.2*B(123)
  JVS(914) = 0.55*B(125)
  JVS(915) = B(5)+B(70)
  JVS(916) = 0.95*B(127)
  JVS(917) = B(12)
  JVS(918) = B(13)
  JVS(919) = 0.22*B(107)+B(109)
  JVS(920) = B(174)
  JVS(921) = 0.6*B(129)
  JVS(922) = B(83)
  JVS(923) = B(17)+2*B(135)+0.76*B(137)
  JVS(924) = 0.8*B(153)-B(159)
  JVS(925) = B(155)-B(161)
  JVS(926) = 0.91*B(151)-B(157)
  JVS(927) = 0.26*B(111)+B(115)
  JVS(928) = 0.07*B(141)
  JVS(929) = 2*B(10)+B(93)+B(95)
  JVS(930) = -B(216)
  JVS(931) = 0.22*B(113)+B(117)
  JVS(932) = 0.5*B(186)+0.5*B(200)-B(214)
  JVS(933) = B(16)
  JVS(934) = B(178)+B(192)-B(206)+0.6*B(226)
  JVS(935) = B(14)
  JVS(936) = B(176)+B(190)-B(204)+0.32*B(225)
  JVS(937) = 0.33*B(20)+0.1*B(147)+0.93*B(149)
  JVS(938) = 0.9*B(19)
  JVS(939) = B(39)-B(41)+0.22*B(108)+0.26*B(112)+0.22*B(114)+0.76*B(138)+0.07*B(142)+0.1*B(148)
  JVS(940) = 0.76*B(180)+0.9*B(194)-B(208)+0.54*B(227)
  JVS(941) = -B(66)-B(68)
  JVS(942) = -B(210)
  JVS(943) = B(40)+B(43)+B(49)-B(57)+B(60)+B(84)+B(86)+B(94)+B(110)+B(116)+B(118)+0.2*B(124)+0.55*B(126)+0.6*B(130)+2&
               &*B(136)+B(175)+B(219)
  JVS(944) = -B(42)-B(58)-2*B(61)-2*B(62)-B(64)-B(67)-B(69)-B(78)-B(158)-B(160)-B(162)-B(205)-B(207)-B(209)-B(211)&
               &-B(212)-B(215)-B(217)
  JVS(945) = -B(65)+0.95*B(128)+0.91*B(152)+0.8*B(154)+B(156)+B(177)+B(179)+0.76*B(181)+0.5*B(187)
  JVS(946) = -B(213)
  JVS(947) = 0.9*B(18)
  JVS(948) = B(50)-B(79)+B(96)+0.93*B(150)+B(191)+B(193)+0.9*B(195)+0.5*B(201)
  JVS(949) = -B(127)
  JVS(950) = B(3)
  JVS(951) = B(29)-B(33)
  JVS(952) = -B(153)
  JVS(953) = -B(155)
  JVS(954) = -B(151)
  JVS(955) = 0
  JVS(956) = -B(188)
  JVS(957) = 0
  JVS(958) = -B(186)
  JVS(959) = 0
  JVS(960) = -B(178)
  JVS(961) = -B(176)
  JVS(962) = 0
  JVS(963) = 0
  JVS(964) = -B(35)
  JVS(965) = -B(180)
  JVS(966) = B(1)+B(30)+B(73)
  JVS(967) = -B(182)
  JVS(968) = -B(45)
  JVS(969) = -B(64)
  JVS(970) = -B(34)-B(36)-B(46)-B(65)-B(71)-B(128)-B(152)-B(154)-B(156)-B(177)-B(179)-B(181)-B(183)-B(184)-B(187)-B(189)
  JVS(971) = -B(185)
  JVS(972) = 0
  JVS(973) = 0.11*B(2)-B(72)+B(74)
  JVS(974) = B(101)
  JVS(975) = 0.11*B(113)
  JVS(976) = 0
  JVS(977) = 0
  JVS(978) = 0.11*B(114)
  JVS(979) = 0
  JVS(980) = 0
  JVS(981) = 0
  JVS(982) = B(102)
  JVS(983) = -B(212)
  JVS(984) = -B(184)
  JVS(985) = -B(185)-B(198)-B(213)-B(229)
  JVS(986) = 0
  JVS(987) = -B(199)
  JVS(988) = B(161)
  JVS(989) = B(157)
  JVS(990) = 0
  JVS(991) = 0
  JVS(992) = 0
  JVS(993) = B(208)
  JVS(994) = 0
  JVS(995) = 0
  JVS(996) = -B(167)
  JVS(997) = B(158)+B(162)+B(209)+B(212)
  JVS(998) = 0
  JVS(999) = B(213)
  JVS(1000) = -B(18)-B(168)
  JVS(1001) = 0
  JVS(1002) = B(6)+B(82)
  JVS(1003) = -B(237)
  JVS(1004) = -B(243)
  JVS(1005) = -B(131)
  JVS(1006) = B(31)
  JVS(1007) = B(53)
  JVS(1008) = -B(119)
  JVS(1009) = -B(143)
  JVS(1010) = -B(95)
  JVS(1011) = -B(202)
  JVS(1012) = -B(121)
  JVS(1013) = -B(200)
  JVS(1014) = -B(105)
  JVS(1015) = -B(192)
  JVS(1016) = -B(99)
  JVS(1017) = -B(190)
  JVS(1018) = -B(149)
  JVS(1019) = 0
  JVS(1020) = B(37)
  JVS(1021) = -B(194)
  JVS(1022) = B(32)+B(38)-B(73)-B(75)
  JVS(1023) = -B(196)
  JVS(1024) = -B(49)+B(54)
  JVS(1025) = -B(78)
  JVS(1026) = -B(71)
  JVS(1027) = -B(198)
  JVS(1028) = 0
  JVS(1029) = -B(2)-B(50)-B(72)-B(74)-B(76)-2*B(77)-B(79)-B(96)-B(100)-B(106)-B(120)-B(122)-B(132)-B(144)-B(150)-B(191)&
                &-B(193)-B(195)-B(197)-B(199)-B(201)-B(203)-B(238)-B(244)
END SUBROUTINE cbmz_mosaic_8bin_vbs9_Jac_SP
SUBROUTINE cbmz_mosaic_8bin_vbs9_KppDecomp( JVS, IER )
      INTEGER :: IER
      REAL(kind=dp) :: JVS(1029), W(136), a
      INTEGER :: k, kk, j, jj
      a = 0.
      IER = 0
      DO k=1,NVAR
        IF ( ABS(JVS(LU_DIAG(k))) < TINY(a) ) THEN
            IER = k
            RETURN
        END IF
        DO kk = LU_CROW(k), LU_CROW(k+1)-1
              W( LU_ICOL(kk) ) = JVS(kk)
        END DO
        DO kk = LU_CROW(k), LU_DIAG(k)-1
            j = LU_ICOL(kk)
            a = -W(j) / JVS( LU_DIAG(j) )
            W(j) = -a
            DO jj = LU_DIAG(j)+1, LU_CROW(j+1)-1
               W( LU_ICOL(jj) ) = W( LU_ICOL(jj) ) + a*JVS(jj)
            END DO
         END DO
         DO kk = LU_CROW(k), LU_CROW(k+1)-1
            JVS(kk) = W( LU_ICOL(kk) )
         END DO
      END DO
END SUBROUTINE cbmz_mosaic_8bin_vbs9_KppDecomp
SUBROUTINE cbmz_mosaic_8bin_vbs9_KppDecompCmplx( JVS, IER )
      INTEGER :: IER
      DOUBLE COMPLEX :: JVS(1029), W(136), a
      INTEGER :: k, kk, j, jj
      IER = 0
      DO k=1,NVAR
        IF ( JVS( LU_DIAG(k) ) .EQ. 0. ) THEN
            IER = k
            RETURN
        END IF
        DO kk = LU_CROW(k), LU_CROW(k+1)-1
              W( LU_ICOL(kk) ) = JVS(kk)
        END DO
        DO kk = LU_CROW(k), LU_DIAG(k)-1
            j = LU_ICOL(kk)
            a = -W(j) / JVS( LU_DIAG(j) )
            W(j) = -a
            DO jj = LU_DIAG(j)+1, LU_CROW(j+1)-1
               W( LU_ICOL(jj) ) = W( LU_ICOL(jj) ) + a*JVS(jj)
            END DO
         END DO
         DO kk = LU_CROW(k), LU_CROW(k+1)-1
            JVS(kk) = W( LU_ICOL(kk) )
         END DO
      END DO
END SUBROUTINE cbmz_mosaic_8bin_vbs9_KppDecompCmplx
SUBROUTINE cbmz_mosaic_8bin_vbs9_KppSolveIndirect( JVS, X )
      INTEGER i, j
      REAL(kind=dp) JVS(1029), X(136), sum
      DO i=1,NVAR
         DO j = LU_CROW(i), LU_DIAG(i)-1
             X(i) = X(i) - JVS(j)*X(LU_ICOL(j));
         END DO
      END DO
      DO i=NVAR,1,-1
        sum = X(i);
        DO j = LU_DIAG(i)+1, LU_CROW(i+1)-1
          sum = sum - JVS(j)*X(LU_ICOL(j));
        END DO
        X(i) = sum/JVS(LU_DIAG(i));
      END DO
END SUBROUTINE cbmz_mosaic_8bin_vbs9_KppSolveIndirect
SUBROUTINE cbmz_mosaic_8bin_vbs9_KppSolveCmplx( JVS, X )
      INTEGER i, j
      DOUBLE COMPLEX JVS(1029), X(136), sum
      DO i=1,NVAR
         DO j = LU_CROW(i), LU_DIAG(i)-1
             X(i) = X(i) - JVS(j)*X(LU_ICOL(j));
         END DO
      END DO
      DO i=NVAR,1,-1
        sum = X(i);
        DO j = LU_DIAG(i)+1, LU_CROW(i+1)-1
          sum = sum - JVS(j)*X(LU_ICOL(j));
        END DO
        X(i) = sum/JVS(LU_DIAG(i));
      END DO
END SUBROUTINE cbmz_mosaic_8bin_vbs9_KppSolveCmplx
SUBROUTINE cbmz_mosaic_8bin_vbs9_KppSolve ( JVS, X )
  REAL(kind=dp) :: JVS(LU_NONZERO)
  REAL(kind=dp) :: X(NVAR)
  X(19) = X(19)-JVS(94)*X(4)-JVS(95)*X(5)-JVS(96)*X(6)-JVS(97)*X(7)
  X(39) = X(39)-JVS(204)*X(38)
  X(47) = X(47)-JVS(222)*X(27)-JVS(223)*X(46)
  X(49) = X(49)-JVS(228)*X(46)
  X(50) = X(50)-JVS(231)*X(26)-JVS(232)*X(45)-JVS(233)*X(47)-JVS(234)*X(49)
  X(51) = X(51)-JVS(237)*X(45)-JVS(238)*X(49)
  X(52) = X(52)-JVS(241)*X(25)-JVS(242)*X(44)-JVS(243)*X(50)-JVS(244)*X(51)
  X(53) = X(53)-JVS(247)*X(44)-JVS(248)*X(51)
  X(54) = X(54)-JVS(251)*X(24)-JVS(252)*X(43)-JVS(253)*X(52)-JVS(254)*X(53)
  X(55) = X(55)-JVS(257)*X(43)-JVS(258)*X(53)
  X(56) = X(56)-JVS(261)*X(23)-JVS(262)*X(42)-JVS(263)*X(54)-JVS(264)*X(55)
  X(57) = X(57)-JVS(267)*X(42)-JVS(268)*X(55)
  X(58) = X(58)-JVS(271)*X(22)-JVS(272)*X(41)-JVS(273)*X(56)-JVS(274)*X(57)
  X(59) = X(59)-JVS(277)*X(41)-JVS(278)*X(57)
  X(60) = X(60)-JVS(281)*X(40)-JVS(282)*X(59)
  X(61) = X(61)-JVS(285)*X(20)-JVS(286)*X(38)-JVS(287)*X(60)
  X(62) = X(62)-JVS(291)*X(21)-JVS(292)*X(40)-JVS(293)*X(58)-JVS(294)*X(59)
  X(64) = X(64)-JVS(299)*X(63)
  X(72) = X(72)-JVS(317)*X(35)-JVS(318)*X(71)
  X(73) = X(73)-JVS(321)*X(71)
  X(74) = X(74)-JVS(324)*X(34)-JVS(325)*X(70)-JVS(326)*X(72)-JVS(327)*X(73)
  X(75) = X(75)-JVS(330)*X(70)-JVS(331)*X(73)
  X(76) = X(76)-JVS(334)*X(33)-JVS(335)*X(69)-JVS(336)*X(74)-JVS(337)*X(75)
  X(77) = X(77)-JVS(340)*X(69)-JVS(341)*X(75)
  X(78) = X(78)-JVS(344)*X(32)-JVS(345)*X(68)-JVS(346)*X(76)-JVS(347)*X(77)
  X(79) = X(79)-JVS(350)*X(68)-JVS(351)*X(77)
  X(80) = X(80)-JVS(354)*X(31)-JVS(355)*X(67)-JVS(356)*X(78)-JVS(357)*X(79)
  X(81) = X(81)-JVS(360)*X(67)-JVS(361)*X(79)
  X(82) = X(82)-JVS(364)*X(30)-JVS(365)*X(66)-JVS(366)*X(80)-JVS(367)*X(81)
  X(83) = X(83)-JVS(370)*X(66)-JVS(371)*X(81)
  X(84) = X(84)-JVS(374)*X(65)-JVS(375)*X(83)
  X(85) = X(85)-JVS(378)*X(28)-JVS(379)*X(63)-JVS(380)*X(84)
  X(86) = X(86)-JVS(384)*X(29)-JVS(385)*X(65)-JVS(386)*X(82)-JVS(387)*X(83)
  X(97) = X(97)-JVS(425)*X(89)-JVS(426)*X(91)
  X(105) = X(105)-JVS(470)*X(89)-JVS(471)*X(91)
  X(106) = X(106)-JVS(475)*X(36)
  X(107) = X(107)-JVS(481)*X(90)-JVS(482)*X(105)
  X(108) = X(108)-JVS(492)*X(103)
  X(109) = X(109)-JVS(507)*X(91)-JVS(508)*X(99)
  X(110) = X(110)-JVS(524)*X(97)-JVS(525)*X(105)
  X(116) = X(116)-JVS(554)*X(100)-JVS(555)*X(103)-JVS(556)*X(104)-JVS(557)*X(110)-JVS(558)*X(112)-JVS(559)*X(113)&
             &-JVS(560)*X(114)-JVS(561)*X(115)
  X(117) = X(117)-JVS(573)*X(89)-JVS(574)*X(91)-JVS(575)*X(103)-JVS(576)*X(105)-JVS(577)*X(110)-JVS(578)*X(114)-JVS(579)&
             &*X(115)
  X(118) = X(118)-JVS(592)*X(112)
  X(120) = X(120)-JVS(609)*X(114)-JVS(610)*X(119)
  X(121) = X(121)-JVS(618)*X(91)-JVS(619)*X(110)-JVS(620)*X(112)-JVS(621)*X(114)-JVS(622)*X(119)
  X(122) = X(122)-JVS(632)*X(98)-JVS(633)*X(101)-JVS(634)*X(114)-JVS(635)*X(119)
  X(123) = X(123)-JVS(645)*X(48)-JVS(646)*X(101)-JVS(647)*X(103)-JVS(648)*X(110)-JVS(649)*X(111)-JVS(650)*X(112)&
             &-JVS(651)*X(114)-JVS(652)*X(115)-JVS(653)*X(119)-JVS(654)*X(120)-JVS(655)*X(122)
  X(124) = X(124)-JVS(666)*X(93)-JVS(667)*X(100)-JVS(668)*X(114)-JVS(669)*X(118)-JVS(670)*X(119)-JVS(671)*X(123)
  X(125) = X(125)-JVS(684)*X(111)-JVS(685)*X(113)-JVS(686)*X(115)
  X(126) = X(126)-JVS(693)*X(92)-JVS(694)*X(97)-JVS(695)*X(105)-JVS(696)*X(111)-JVS(697)*X(113)-JVS(698)*X(115)-JVS(699)&
             &*X(120)-JVS(700)*X(125)
  X(127) = X(127)-JVS(709)*X(94)-JVS(710)*X(95)-JVS(711)*X(103)-JVS(712)*X(106)-JVS(713)*X(110)-JVS(714)*X(114)-JVS(715)&
             &*X(115)-JVS(716)*X(119)-JVS(717)*X(125)
  X(128) = X(128)-JVS(725)*X(109)-JVS(726)*X(111)-JVS(727)*X(113)-JVS(728)*X(114)-JVS(729)*X(115)-JVS(730)*X(119)&
             &-JVS(731)*X(120)-JVS(732)*X(125)-JVS(733)*X(126)-JVS(734)*X(127)
  X(129) = X(129)-JVS(743)*X(87)-JVS(744)*X(90)-JVS(745)*X(92)-JVS(746)*X(96)-JVS(747)*X(97)-JVS(748)*X(102)-JVS(749)&
             &*X(105)-JVS(750)*X(106)-JVS(751)*X(107)-JVS(752)*X(111)-JVS(753)*X(112)-JVS(754)*X(113)-JVS(755)*X(115)&
             &-JVS(756)*X(116)-JVS(757)*X(117)-JVS(758)*X(119)-JVS(759)*X(120)-JVS(760)*X(121)-JVS(761)*X(122)-JVS(762)&
             &*X(123)-JVS(763)*X(124)-JVS(764)*X(125)-JVS(765)*X(126)-JVS(766)*X(127)-JVS(767)*X(128)
  X(130) = X(130)-JVS(776)*X(87)-JVS(777)*X(110)-JVS(778)*X(114)-JVS(779)*X(115)-JVS(780)*X(118)-JVS(781)*X(119)&
             &-JVS(782)*X(121)-JVS(783)*X(123)-JVS(784)*X(125)-JVS(785)*X(126)-JVS(786)*X(127)-JVS(787)*X(128)-JVS(788)&
             &*X(129)
  X(131) = X(131)-JVS(796)*X(4)-JVS(797)*X(5)-JVS(798)*X(6)-JVS(799)*X(7)-JVS(800)*X(20)-JVS(801)*X(21)-JVS(802)*X(22)&
             &-JVS(803)*X(23)-JVS(804)*X(24)-JVS(805)*X(25)-JVS(806)*X(26)-JVS(807)*X(27)-JVS(808)*X(28)-JVS(809)*X(29)&
             &-JVS(810)*X(30)-JVS(811)*X(31)-JVS(812)*X(32)-JVS(813)*X(33)-JVS(814)*X(34)-JVS(815)*X(35)-JVS(816)*X(36)&
             &-JVS(817)*X(37)-JVS(818)*X(38)-JVS(819)*X(39)-JVS(820)*X(40)-JVS(821)*X(41)-JVS(822)*X(42)-JVS(823)*X(43)&
             &-JVS(824)*X(44)-JVS(825)*X(45)-JVS(826)*X(46)-JVS(827)*X(47)-JVS(828)*X(48)-JVS(829)*X(49)-JVS(830)*X(50)&
             &-JVS(831)*X(51)-JVS(832)*X(52)-JVS(833)*X(53)-JVS(834)*X(54)-JVS(835)*X(55)-JVS(836)*X(56)-JVS(837)*X(57)&
             &-JVS(838)*X(58)-JVS(839)*X(59)-JVS(840)*X(60)-JVS(841)*X(61)-JVS(842)*X(62)-JVS(843)*X(63)-JVS(844)*X(64)&
             &-JVS(845)*X(65)-JVS(846)*X(66)-JVS(847)*X(67)-JVS(848)*X(68)-JVS(849)*X(69)-JVS(850)*X(70)-JVS(851)*X(71)&
             &-JVS(852)*X(72)-JVS(853)*X(73)-JVS(854)*X(74)-JVS(855)*X(75)-JVS(856)*X(76)-JVS(857)*X(77)-JVS(858)*X(78)&
             &-JVS(859)*X(79)-JVS(860)*X(80)-JVS(861)*X(81)-JVS(862)*X(82)-JVS(863)*X(83)-JVS(864)*X(84)-JVS(865)*X(85)&
             &-JVS(866)*X(86)-JVS(867)*X(88)-JVS(868)*X(89)-JVS(869)*X(91)-JVS(870)*X(93)-JVS(871)*X(94)-JVS(872)*X(95)&
             &-JVS(873)*X(96)-JVS(874)*X(98)-JVS(875)*X(100)-JVS(876)*X(101)-JVS(877)*X(102)-JVS(878)*X(103)-JVS(879)*X(104)&
             &-JVS(880)*X(105)-JVS(881)*X(107)-JVS(882)*X(108)-JVS(883)*X(109)-JVS(884)*X(110)-JVS(885)*X(111)-JVS(886)&
             &*X(112)-JVS(887)*X(113)-JVS(888)*X(114)-JVS(889)*X(115)-JVS(890)*X(116)-JVS(891)*X(118)-JVS(892)*X(119)&
             &-JVS(893)*X(120)-JVS(894)*X(121)-JVS(895)*X(122)-JVS(896)*X(123)-JVS(897)*X(124)-JVS(898)*X(125)-JVS(899)&
             &*X(126)-JVS(900)*X(127)-JVS(901)*X(128)-JVS(902)*X(129)-JVS(903)*X(130)
  X(132) = X(132)-JVS(910)*X(37)-JVS(911)*X(48)-JVS(912)*X(88)-JVS(913)*X(89)-JVS(914)*X(91)-JVS(915)*X(96)-JVS(916)&
             &*X(97)-JVS(917)*X(100)-JVS(918)*X(101)-JVS(919)*X(103)-JVS(920)*X(104)-JVS(921)*X(105)-JVS(922)*X(108)&
             &-JVS(923)*X(110)-JVS(924)*X(111)-JVS(925)*X(112)-JVS(926)*X(113)-JVS(927)*X(114)-JVS(928)*X(115)-JVS(929)&
             &*X(116)-JVS(930)*X(117)-JVS(931)*X(119)-JVS(932)*X(120)-JVS(933)*X(121)-JVS(934)*X(122)-JVS(935)*X(123)&
             &-JVS(936)*X(124)-JVS(937)*X(125)-JVS(938)*X(126)-JVS(939)*X(127)-JVS(940)*X(128)-JVS(941)*X(129)-JVS(942)&
             &*X(130)-JVS(943)*X(131)
  X(133) = X(133)-JVS(949)*X(97)-JVS(950)*X(102)-JVS(951)*X(106)-JVS(952)*X(111)-JVS(953)*X(112)-JVS(954)*X(113)&
             &-JVS(955)*X(115)-JVS(956)*X(117)-JVS(957)*X(119)-JVS(958)*X(120)-JVS(959)*X(121)-JVS(960)*X(122)-JVS(961)&
             &*X(124)-JVS(962)*X(125)-JVS(963)*X(126)-JVS(964)*X(127)-JVS(965)*X(128)-JVS(966)*X(129)-JVS(967)*X(130)&
             &-JVS(968)*X(131)-JVS(969)*X(132)
  X(134) = X(134)-JVS(974)*X(118)-JVS(975)*X(119)-JVS(976)*X(125)-JVS(977)*X(126)-JVS(978)*X(127)-JVS(979)*X(128)&
             &-JVS(980)*X(129)-JVS(981)*X(130)-JVS(982)*X(131)-JVS(983)*X(132)-JVS(984)*X(133)
  X(135) = X(135)-JVS(988)*X(112)-JVS(989)*X(113)-JVS(990)*X(115)-JVS(991)*X(125)-JVS(992)*X(127)-JVS(993)*X(128)&
             &-JVS(994)*X(129)-JVS(995)*X(130)-JVS(996)*X(131)-JVS(997)*X(132)-JVS(998)*X(133)-JVS(999)*X(134)
  X(136) = X(136)-JVS(1002)*X(90)-JVS(1003)*X(94)-JVS(1004)*X(95)-JVS(1005)*X(105)-JVS(1006)*X(106)-JVS(1007)*X(107)&
             &-JVS(1008)*X(114)-JVS(1009)*X(115)-JVS(1010)*X(116)-JVS(1011)*X(117)-JVS(1012)*X(119)-JVS(1013)*X(120)&
             &-JVS(1014)*X(121)-JVS(1015)*X(122)-JVS(1016)*X(123)-JVS(1017)*X(124)-JVS(1018)*X(125)-JVS(1019)*X(126)&
             &-JVS(1020)*X(127)-JVS(1021)*X(128)-JVS(1022)*X(129)-JVS(1023)*X(130)-JVS(1024)*X(131)-JVS(1025)*X(132)&
             &-JVS(1026)*X(133)-JVS(1027)*X(134)-JVS(1028)*X(135)
  X(136) = X(136)/JVS(1029)
  X(135) = (X(135)-JVS(1001)*X(136))/(JVS(1000))
  X(134) = (X(134)-JVS(986)*X(135)-JVS(987)*X(136))/(JVS(985))
  X(133) = (X(133)-JVS(971)*X(134)-JVS(972)*X(135)-JVS(973)*X(136))/(JVS(970))
  X(132) = (X(132)-JVS(945)*X(133)-JVS(946)*X(134)-JVS(947)*X(135)-JVS(948)*X(136))/(JVS(944))
  X(131) = (X(131)-JVS(905)*X(132)-JVS(906)*X(133)-JVS(907)*X(134)-JVS(908)*X(135)-JVS(909)*X(136))/(JVS(904))
  X(130) = (X(130)-JVS(790)*X(131)-JVS(791)*X(132)-JVS(792)*X(133)-JVS(793)*X(134)-JVS(794)*X(135)-JVS(795)*X(136))&
             &/(JVS(789))
  X(129) = (X(129)-JVS(769)*X(130)-JVS(770)*X(131)-JVS(771)*X(132)-JVS(772)*X(133)-JVS(773)*X(134)-JVS(774)*X(135)&
             &-JVS(775)*X(136))/(JVS(768))
  X(128) = (X(128)-JVS(736)*X(129)-JVS(737)*X(130)-JVS(738)*X(131)-JVS(739)*X(132)-JVS(740)*X(133)-JVS(741)*X(135)&
             &-JVS(742)*X(136))/(JVS(735))
  X(127) = (X(127)-JVS(719)*X(129)-JVS(720)*X(130)-JVS(721)*X(131)-JVS(722)*X(132)-JVS(723)*X(133)-JVS(724)*X(136))&
             &/(JVS(718))
  X(126) = (X(126)-JVS(702)*X(127)-JVS(703)*X(128)-JVS(704)*X(129)-JVS(705)*X(131)-JVS(706)*X(132)-JVS(707)*X(133)&
             &-JVS(708)*X(136))/(JVS(701))
  X(125) = (X(125)-JVS(688)*X(127)-JVS(689)*X(131)-JVS(690)*X(132)-JVS(691)*X(133)-JVS(692)*X(136))/(JVS(687))
  X(124) = (X(124)-JVS(673)*X(125)-JVS(674)*X(126)-JVS(675)*X(127)-JVS(676)*X(128)-JVS(677)*X(130)-JVS(678)*X(131)&
             &-JVS(679)*X(132)-JVS(680)*X(133)-JVS(681)*X(134)-JVS(682)*X(135)-JVS(683)*X(136))/(JVS(672))
  X(123) = (X(123)-JVS(657)*X(125)-JVS(658)*X(126)-JVS(659)*X(127)-JVS(660)*X(128)-JVS(661)*X(131)-JVS(662)*X(132)&
             &-JVS(663)*X(133)-JVS(664)*X(135)-JVS(665)*X(136))/(JVS(656))
  X(122) = (X(122)-JVS(637)*X(126)-JVS(638)*X(127)-JVS(639)*X(128)-JVS(640)*X(131)-JVS(641)*X(132)-JVS(642)*X(133)&
             &-JVS(643)*X(135)-JVS(644)*X(136))/(JVS(636))
  X(121) = (X(121)-JVS(624)*X(125)-JVS(625)*X(127)-JVS(626)*X(131)-JVS(627)*X(132)-JVS(628)*X(133)-JVS(629)*X(134)&
             &-JVS(630)*X(135)-JVS(631)*X(136))/(JVS(623))
  X(120) = (X(120)-JVS(612)*X(126)-JVS(613)*X(127)-JVS(614)*X(131)-JVS(615)*X(132)-JVS(616)*X(133)-JVS(617)*X(136))&
             &/(JVS(611))
  X(119) = (X(119)-JVS(606)*X(127)-JVS(607)*X(131)-JVS(608)*X(136))/(JVS(605))
  X(118) = (X(118)-JVS(594)*X(119)-JVS(595)*X(125)-JVS(596)*X(126)-JVS(597)*X(127)-JVS(598)*X(128)-JVS(599)*X(131)&
             &-JVS(600)*X(132)-JVS(601)*X(133)-JVS(602)*X(134)-JVS(603)*X(135)-JVS(604)*X(136))/(JVS(593))
  X(117) = (X(117)-JVS(581)*X(119)-JVS(582)*X(121)-JVS(583)*X(125)-JVS(584)*X(126)-JVS(585)*X(127)-JVS(586)*X(128)&
             &-JVS(587)*X(131)-JVS(588)*X(132)-JVS(589)*X(133)-JVS(590)*X(135)-JVS(591)*X(136))/(JVS(580))
  X(116) = (X(116)-JVS(563)*X(119)-JVS(564)*X(120)-JVS(565)*X(124)-JVS(566)*X(125)-JVS(567)*X(127)-JVS(568)*X(131)&
             &-JVS(569)*X(132)-JVS(570)*X(133)-JVS(571)*X(134)-JVS(572)*X(136))/(JVS(562))
  X(115) = (X(115)-JVS(551)*X(127)-JVS(552)*X(131)-JVS(553)*X(136))/(JVS(550))
  X(114) = (X(114)-JVS(547)*X(127)-JVS(548)*X(131)-JVS(549)*X(136))/(JVS(546))
  X(113) = (X(113)-JVS(542)*X(115)-JVS(543)*X(131)-JVS(544)*X(132)-JVS(545)*X(133))/(JVS(541))
  X(112) = (X(112)-JVS(537)*X(125)-JVS(538)*X(131)-JVS(539)*X(132)-JVS(540)*X(133))/(JVS(536))
  X(111) = (X(111)-JVS(532)*X(115)-JVS(533)*X(132)-JVS(534)*X(133)-JVS(535)*X(136))/(JVS(531))
  X(110) = (X(110)-JVS(527)*X(127)-JVS(528)*X(131)-JVS(529)*X(133)-JVS(530)*X(136))/(JVS(526))
  X(109) = (X(109)-JVS(510)*X(111)-JVS(511)*X(113)-JVS(512)*X(114)-JVS(513)*X(119)-JVS(514)*X(120)-JVS(515)*X(125)&
             &-JVS(516)*X(126)-JVS(517)*X(127)-JVS(518)*X(128)-JVS(519)*X(131)-JVS(520)*X(132)-JVS(521)*X(133)-JVS(522)&
             &*X(135)-JVS(523)*X(136))/(JVS(509))
  X(108) = (X(108)-JVS(494)*X(110)-JVS(495)*X(112)-JVS(496)*X(114)-JVS(497)*X(115)-JVS(498)*X(116)-JVS(499)*X(119)&
             &-JVS(500)*X(121)-JVS(501)*X(123)-JVS(502)*X(125)-JVS(503)*X(127)-JVS(504)*X(131)-JVS(505)*X(133)-JVS(506)&
             &*X(136))/(JVS(493))
  X(107) = (X(107)-JVS(484)*X(116)-JVS(485)*X(121)-JVS(486)*X(123)-JVS(487)*X(125)-JVS(488)*X(129)-JVS(489)*X(131)&
             &-JVS(490)*X(132)-JVS(491)*X(136))/(JVS(483))
  X(106) = (X(106)-JVS(477)*X(127)-JVS(478)*X(129)-JVS(479)*X(133)-JVS(480)*X(136))/(JVS(476))
  X(105) = (X(105)-JVS(473)*X(131)-JVS(474)*X(136))/(JVS(472))
  X(104) = (X(104)-JVS(465)*X(114)-JVS(466)*X(119)-JVS(467)*X(124)-JVS(468)*X(127)-JVS(469)*X(131))/(JVS(464))
  X(103) = (X(103)-JVS(462)*X(127)-JVS(463)*X(131))/(JVS(461))
  X(102) = (X(102)-JVS(457)*X(129)-JVS(458)*X(131)-JVS(459)*X(132)-JVS(460)*X(133))/(JVS(456))
  X(101) = (X(101)-JVS(453)*X(122)-JVS(454)*X(131)-JVS(455)*X(132))/(JVS(452))
  X(100) = (X(100)-JVS(449)*X(124)-JVS(450)*X(131)-JVS(451)*X(132))/(JVS(448))
  X(99) = (X(99)-JVS(437)*X(109)-JVS(438)*X(114)-JVS(439)*X(119)-JVS(440)*X(120)-JVS(441)*X(126)-JVS(442)*X(127)&
            &-JVS(443)*X(128)-JVS(444)*X(131)-JVS(445)*X(133)-JVS(446)*X(135)-JVS(447)*X(136))/(JVS(436))
  X(98) = (X(98)-JVS(431)*X(114)-JVS(432)*X(119)-JVS(433)*X(122)-JVS(434)*X(127)-JVS(435)*X(131))/(JVS(430))
  X(97) = (X(97)-JVS(428)*X(131)-JVS(429)*X(133))/(JVS(427))
  X(96) = (X(96)-JVS(422)*X(129)-JVS(423)*X(131)-JVS(424)*X(132))/(JVS(421))
  X(95) = (X(95)-JVS(418)*X(127)-JVS(419)*X(131)-JVS(420)*X(136))/(JVS(417))
  X(94) = (X(94)-JVS(414)*X(127)-JVS(415)*X(131)-JVS(416)*X(136))/(JVS(413))
  X(93) = (X(93)-JVS(409)*X(114)-JVS(410)*X(119)-JVS(411)*X(127)-JVS(412)*X(131))/(JVS(408))
  X(92) = (X(92)-JVS(404)*X(105)-JVS(405)*X(129)-JVS(406)*X(131)-JVS(407)*X(136))/(JVS(403))
  X(91) = (X(91)-JVS(402)*X(131))/(JVS(401))
  X(90) = (X(90)-JVS(399)*X(129)-JVS(400)*X(136))/(JVS(398))
  X(89) = (X(89)-JVS(397)*X(131))/(JVS(396))
  X(88) = (X(88)-JVS(394)*X(131)-JVS(395)*X(132))/(JVS(393))
  X(87) = (X(87)-JVS(391)*X(129)-JVS(392)*X(130))/(JVS(390))
  X(86) = (X(86)-JVS(389)*X(131))/(JVS(388))
  X(85) = (X(85)-JVS(382)*X(86)-JVS(383)*X(131))/(JVS(381))
  X(84) = (X(84)-JVS(377)*X(131))/(JVS(376))
  X(83) = (X(83)-JVS(373)*X(131))/(JVS(372))
  X(82) = (X(82)-JVS(369)*X(131))/(JVS(368))
  X(81) = (X(81)-JVS(363)*X(131))/(JVS(362))
  X(80) = (X(80)-JVS(359)*X(131))/(JVS(358))
  X(79) = (X(79)-JVS(353)*X(131))/(JVS(352))
  X(78) = (X(78)-JVS(349)*X(131))/(JVS(348))
  X(77) = (X(77)-JVS(343)*X(131))/(JVS(342))
  X(76) = (X(76)-JVS(339)*X(131))/(JVS(338))
  X(75) = (X(75)-JVS(333)*X(131))/(JVS(332))
  X(74) = (X(74)-JVS(329)*X(131))/(JVS(328))
  X(73) = (X(73)-JVS(323)*X(131))/(JVS(322))
  X(72) = (X(72)-JVS(320)*X(131))/(JVS(319))
  X(71) = (X(71)-JVS(316)*X(131))/(JVS(315))
  X(70) = (X(70)-JVS(314)*X(131))/(JVS(313))
  X(69) = (X(69)-JVS(312)*X(131))/(JVS(311))
  X(68) = (X(68)-JVS(310)*X(131))/(JVS(309))
  X(67) = (X(67)-JVS(308)*X(131))/(JVS(307))
  X(66) = (X(66)-JVS(306)*X(131))/(JVS(305))
  X(65) = (X(65)-JVS(304)*X(131))/(JVS(303))
  X(64) = (X(64)-JVS(301)*X(84)-JVS(302)*X(131))/(JVS(300))
  X(63) = (X(63)-JVS(298)*X(131))/(JVS(297))
  X(62) = (X(62)-JVS(296)*X(131))/(JVS(295))
  X(61) = (X(61)-JVS(289)*X(62)-JVS(290)*X(131))/(JVS(288))
  X(60) = (X(60)-JVS(284)*X(131))/(JVS(283))
  X(59) = (X(59)-JVS(280)*X(131))/(JVS(279))
  X(58) = (X(58)-JVS(276)*X(131))/(JVS(275))
  X(57) = (X(57)-JVS(270)*X(131))/(JVS(269))
  X(56) = (X(56)-JVS(266)*X(131))/(JVS(265))
  X(55) = (X(55)-JVS(260)*X(131))/(JVS(259))
  X(54) = (X(54)-JVS(256)*X(131))/(JVS(255))
  X(53) = (X(53)-JVS(250)*X(131))/(JVS(249))
  X(52) = (X(52)-JVS(246)*X(131))/(JVS(245))
  X(51) = (X(51)-JVS(240)*X(131))/(JVS(239))
  X(50) = (X(50)-JVS(236)*X(131))/(JVS(235))
  X(49) = (X(49)-JVS(230)*X(131))/(JVS(229))
  X(48) = (X(48)-JVS(227)*X(131))/(JVS(226))
  X(47) = (X(47)-JVS(225)*X(131))/(JVS(224))
  X(46) = (X(46)-JVS(221)*X(131))/(JVS(220))
  X(45) = (X(45)-JVS(219)*X(131))/(JVS(218))
  X(44) = (X(44)-JVS(217)*X(131))/(JVS(216))
  X(43) = (X(43)-JVS(215)*X(131))/(JVS(214))
  X(42) = (X(42)-JVS(213)*X(131))/(JVS(212))
  X(41) = (X(41)-JVS(211)*X(131))/(JVS(210))
  X(40) = (X(40)-JVS(209)*X(131))/(JVS(208))
  X(39) = (X(39)-JVS(206)*X(60)-JVS(207)*X(131))/(JVS(205))
  X(38) = (X(38)-JVS(203)*X(131))/(JVS(202))
  X(37) = (X(37)-JVS(201)*X(131))/(JVS(200))
  X(36) = (X(36)-JVS(199)*X(127))/(JVS(198))
  X(35) = (X(35)-JVS(197)*X(131))/(JVS(196))
  X(34) = (X(34)-JVS(195)*X(131))/(JVS(194))
  X(33) = (X(33)-JVS(193)*X(131))/(JVS(192))
  X(32) = (X(32)-JVS(191)*X(131))/(JVS(190))
  X(31) = (X(31)-JVS(189)*X(131))/(JVS(188))
  X(30) = (X(30)-JVS(187)*X(131))/(JVS(186))
  X(29) = (X(29)-JVS(185)*X(131))/(JVS(184))
  X(28) = (X(28)-JVS(183)*X(131))/(JVS(182))
  X(27) = (X(27)-JVS(181)*X(131))/(JVS(180))
  X(26) = (X(26)-JVS(179)*X(131))/(JVS(178))
  X(25) = (X(25)-JVS(177)*X(131))/(JVS(176))
  X(24) = (X(24)-JVS(175)*X(131))/(JVS(174))
  X(23) = (X(23)-JVS(173)*X(131))/(JVS(172))
  X(22) = (X(22)-JVS(171)*X(131))/(JVS(170))
  X(21) = (X(21)-JVS(169)*X(131))/(JVS(168))
  X(20) = (X(20)-JVS(167)*X(131))/(JVS(166))
  X(19) = (X(19)-JVS(99)*X(20)-JVS(100)*X(21)-JVS(101)*X(22)-JVS(102)*X(23)-JVS(103)*X(24)-JVS(104)*X(25)-JVS(105)*X(26)&
            &-JVS(106)*X(27)-JVS(107)*X(28)-JVS(108)*X(29)-JVS(109)*X(30)-JVS(110)*X(31)-JVS(111)*X(32)-JVS(112)*X(33)&
            &-JVS(113)*X(34)-JVS(114)*X(35)-JVS(115)*X(38)-JVS(116)*X(39)-JVS(117)*X(40)-JVS(118)*X(41)-JVS(119)*X(42)&
            &-JVS(120)*X(43)-JVS(121)*X(44)-JVS(122)*X(45)-JVS(123)*X(46)-JVS(124)*X(47)-JVS(125)*X(49)-JVS(126)*X(50)&
            &-JVS(127)*X(51)-JVS(128)*X(52)-JVS(129)*X(53)-JVS(130)*X(54)-JVS(131)*X(55)-JVS(132)*X(56)-JVS(133)*X(57)&
            &-JVS(134)*X(58)-JVS(135)*X(59)-JVS(136)*X(60)-JVS(137)*X(61)-JVS(138)*X(62)-JVS(139)*X(63)-JVS(140)*X(64)&
            &-JVS(141)*X(65)-JVS(142)*X(66)-JVS(143)*X(67)-JVS(144)*X(68)-JVS(145)*X(69)-JVS(146)*X(70)-JVS(147)*X(71)&
            &-JVS(148)*X(72)-JVS(149)*X(73)-JVS(150)*X(74)-JVS(151)*X(75)-JVS(152)*X(76)-JVS(153)*X(77)-JVS(154)*X(78)&
            &-JVS(155)*X(79)-JVS(156)*X(80)-JVS(157)*X(81)-JVS(158)*X(82)-JVS(159)*X(83)-JVS(160)*X(84)-JVS(161)*X(85)&
            &-JVS(162)*X(86)-JVS(163)*X(94)-JVS(164)*X(95)-JVS(165)*X(131))/(JVS(98))
  X(18) = (X(18)-JVS(60)*X(37)-JVS(61)*X(48)-JVS(62)*X(88)-JVS(63)*X(89)-JVS(64)*X(91)-JVS(65)*X(93)-JVS(66)*X(96)&
            &-JVS(67)*X(98)-JVS(68)*X(100)-JVS(69)*X(101)-JVS(70)*X(102)-JVS(71)*X(103)-JVS(72)*X(104)-JVS(73)*X(105)&
            &-JVS(74)*X(107)-JVS(75)*X(108)-JVS(76)*X(109)-JVS(77)*X(110)-JVS(78)*X(114)-JVS(79)*X(115)-JVS(80)*X(116)&
            &-JVS(81)*X(118)-JVS(82)*X(119)-JVS(83)*X(121)-JVS(84)*X(123)-JVS(85)*X(125)-JVS(86)*X(126)-JVS(87)*X(127)&
            &-JVS(88)*X(129)-JVS(89)*X(131)-JVS(90)*X(132)-JVS(91)*X(133)-JVS(92)*X(135)-JVS(93)*X(136))/(JVS(59))
  X(17) = (X(17)-JVS(56)*X(95)-JVS(57)*X(127)-JVS(58)*X(131))/(JVS(55))
  X(16) = (X(16)-JVS(52)*X(95)-JVS(53)*X(127)-JVS(54)*X(131))/(JVS(51))
  X(15) = (X(15)-JVS(48)*X(94)-JVS(49)*X(127)-JVS(50)*X(131))/(JVS(47))
  X(14) = (X(14)-JVS(44)*X(94)-JVS(45)*X(127)-JVS(46)*X(131))/(JVS(43))
  X(13) = (X(13)-JVS(38)*X(114)-JVS(39)*X(119)-JVS(40)*X(127)-JVS(41)*X(131)-JVS(42)*X(136))/(JVS(37))
  X(12) = (X(12)-JVS(35)*X(109)-JVS(36)*X(131))/(JVS(34))
  X(11) = (X(11)-JVS(29)*X(89)-JVS(30)*X(91)-JVS(31)*X(105)-JVS(32)*X(131)-JVS(33)*X(136))/(JVS(28))
  X(10) = (X(10)-JVS(23)*X(89)-JVS(24)*X(91)-JVS(25)*X(105)-JVS(26)*X(131)-JVS(27)*X(136))/(JVS(22))
  X(9) = (X(9)-JVS(15)*X(114)-JVS(16)*X(115)-JVS(17)*X(119)-JVS(18)*X(125)-JVS(19)*X(127)-JVS(20)*X(130)-JVS(21)*X(132))&
           &/(JVS(14))
  X(8) = (X(8)-JVS(11)*X(103)-JVS(12)*X(114)-JVS(13)*X(127))/(JVS(10))
  X(7) = X(7)/JVS(9)
  X(6) = X(6)/JVS(8)
  X(5) = X(5)/JVS(7)
  X(4) = X(4)/JVS(6)
  X(3) = X(3)/JVS(5)
  X(2) = X(2)/JVS(4)
  X(1) = (X(1)-JVS(2)*X(37)-JVS(3)*X(131))/(JVS(1))
END SUBROUTINE cbmz_mosaic_8bin_vbs9_KppSolve
      SUBROUTINE cbmz_mosaic_8bin_vbs9_WCOPY(N,X,incX,Y,incY)
      INTEGER i,incX,incY,M,MP1,N
      REAL(kind=dp) X(N),Y(N)
      IF (N.LE.0) RETURN
      M = MOD(N,8)
      IF( M .NE. 0 ) THEN
        DO i = 1,M
          Y(i) = X(i)
        END DO
        IF( N .LT. 8 ) RETURN
      END IF
      MP1 = M+1
      DO i = MP1,N,8
        Y(i) = X(i)
        Y(i + 1) = X(i + 1)
        Y(i + 2) = X(i + 2)
        Y(i + 3) = X(i + 3)
        Y(i + 4) = X(i + 4)
        Y(i + 5) = X(i + 5)
        Y(i + 6) = X(i + 6)
        Y(i + 7) = X(i + 7)
      END DO
      END SUBROUTINE cbmz_mosaic_8bin_vbs9_WCOPY
      SUBROUTINE cbmz_mosaic_8bin_vbs9_WAXPY(N,Alpha,X,incX,Y,incY)
      INTEGER i,incX,incY,M,MP1,N
      REAL(kind=dp) X(N),Y(N),Alpha
      REAL(kind=dp) ZERO
      PARAMETER( ZERO = 0.0_dp )
      IF (Alpha .EQ. ZERO) RETURN
      IF (N .LE. 0) RETURN
      M = MOD(N,4)
      IF( M .NE. 0 ) THEN
        DO i = 1,M
          Y(i) = Y(i) + Alpha*X(i)
        END DO
        IF( N .LT. 4 ) RETURN
      END IF
      MP1 = M + 1
      DO i = MP1,N,4
        Y(i) = Y(i) + Alpha*X(i)
        Y(i + 1) = Y(i + 1) + Alpha*X(i + 1)
        Y(i + 2) = Y(i + 2) + Alpha*X(i + 2)
        Y(i + 3) = Y(i + 3) + Alpha*X(i + 3)
      END DO
      END SUBROUTINE cbmz_mosaic_8bin_vbs9_WAXPY
      SUBROUTINE cbmz_mosaic_8bin_vbs9_WSCAL(N,Alpha,X,incX)
      INTEGER i,incX,M,MP1,N
      REAL(kind=dp) X(N),Alpha
      REAL(kind=dp) ZERO, ONE
      PARAMETER( ZERO = 0.0_dp )
      PARAMETER( ONE = 1.0_dp )
      IF (Alpha .EQ. ONE) RETURN
      IF (N .LE. 0) RETURN
      M = MOD(N,5)
      IF( M .NE. 0 ) THEN
        IF (Alpha .EQ. (-ONE)) THEN
          DO i = 1,M
            X(i) = -X(i)
          END DO
        ELSEIF (Alpha .EQ. ZERO) THEN
          DO i = 1,M
            X(i) = ZERO
          END DO
        ELSE
          DO i = 1,M
            X(i) = Alpha*X(i)
          END DO
        END IF
        IF( N .LT. 5 ) RETURN
      END IF
      MP1 = M + 1
      IF (Alpha .EQ. (-ONE)) THEN
        DO i = MP1,N,5
          X(i) = -X(i)
          X(i + 1) = -X(i + 1)
          X(i + 2) = -X(i + 2)
          X(i + 3) = -X(i + 3)
          X(i + 4) = -X(i + 4)
        END DO
      ELSEIF (Alpha .EQ. ZERO) THEN
        DO i = MP1,N,5
          X(i) = ZERO
          X(i + 1) = ZERO
          X(i + 2) = ZERO
          X(i + 3) = ZERO
          X(i + 4) = ZERO
        END DO
      ELSE
        DO i = MP1,N,5
          X(i) = Alpha*X(i)
          X(i + 1) = Alpha*X(i + 1)
          X(i + 2) = Alpha*X(i + 2)
          X(i + 3) = Alpha*X(i + 3)
          X(i + 4) = Alpha*X(i + 4)
        END DO
      END IF
      END SUBROUTINE cbmz_mosaic_8bin_vbs9_WSCAL
      REAL(kind=dp) FUNCTION cbmz_mosaic_8bin_vbs9_WLAMCH( C )
      CHARACTER C
      INTEGER i
      REAL(kind=dp) ONE, HALF, Eps, Sum
      PARAMETER (ONE = 1.0_dp)
      PARAMETER (HALF = 0.5_dp)
      LOGICAL First
      SAVE First, Eps
      DATA First /.TRUE./
      IF (First) THEN
        First = .FALSE.
        Eps = HALF**(16)
        DO i = 17, 80
          Eps = Eps*HALF
          CALL cbmz_mosaic_8bin_vbs9_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10 Eps = Eps*2
        i = i-1
      END IF
      cbmz_mosaic_8bin_vbs9_WLAMCH = Eps
      END FUNCTION cbmz_mosaic_8bin_vbs9_WLAMCH
      SUBROUTINE cbmz_mosaic_8bin_vbs9_WLAMCH_ADD( A, B, Sum )
      REAL(kind=dp) A, B, Sum
      Sum = A + B
      END SUBROUTINE cbmz_mosaic_8bin_vbs9_WLAMCH_ADD
      SUBROUTINE cbmz_mosaic_8bin_vbs9_SET2ZERO(N,Y)
      INTEGER :: i,M,MP1,N
      REAL(kind=dp) :: Y(N)
      REAL(kind=dp), PARAMETER :: ZERO = 0.0d0
      IF (N.LE.0) RETURN
      M = MOD(N,8)
      IF( M .NE. 0 ) THEN
        DO i = 1,M
          Y(i) = ZERO
        END DO
        IF( N .LT. 8 ) RETURN
      END IF
      MP1 = M+1
      DO i = MP1,N,8
        Y(i) = ZERO
        Y(i + 1) = ZERO
        Y(i + 2) = ZERO
        Y(i + 3) = ZERO
        Y(i + 4) = ZERO
        Y(i + 5) = ZERO
        Y(i + 6) = ZERO
        Y(i + 7) = ZERO
      END DO
      END SUBROUTINE cbmz_mosaic_8bin_vbs9_SET2ZERO
      REAL(kind=dp) FUNCTION cbmz_mosaic_8bin_vbs9_WDOT (N, DX, incX, DY, incY)
      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N)
      INTEGER :: i, IX, IY, M, MP1, NS
      cbmz_mosaic_8bin_vbs9_WDOT = 0.0D0
      IF (N .LE. 0) RETURN
      IF (incX .EQ. incY) IF (incX-1) 5,20,60
    5 IX = 1
      IY = 1
      IF (incX .LT. 0) IX = (-N+1)*incX + 1
      IF (incY .LT. 0) IY = (-N+1)*incY + 1
      DO i = 1,N
        cbmz_mosaic_8bin_vbs9_WDOT = cbmz_mosaic_8bin_vbs9_WDOT + DX(IX)*DY(IY)
        IX = IX + incX
        IY = IY + incY
      END DO
      RETURN
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO i = 1,M
         cbmz_mosaic_8bin_vbs9_WDOT = cbmz_mosaic_8bin_vbs9_WDOT + DX(i)*DY(i)
      END DO
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO i = MP1,N,5
          cbmz_mosaic_8bin_vbs9_WDOT = cbmz_mosaic_8bin_vbs9_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) + &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)
      END DO
      RETURN
   60 NS = N*incX
      DO i = 1,NS,incX
        cbmz_mosaic_8bin_vbs9_WDOT = cbmz_mosaic_8bin_vbs9_WDOT + DX(i)*DY(i)
      END DO
      END FUNCTION cbmz_mosaic_8bin_vbs9_WDOT
   SUBROUTINE decomp_cbmz_mosaic_8bin_vbs9( JVS, IER )
     IMPLICIT NONE
      INTEGER :: IER
      REAL(kind=dp) :: JVS(LU_NONZERO), W(NVAR), a
  a = 0._dp
  ier = 0
  IF ( ABS( JVS( 1 )) < TINY(a) ) THEN
         IER = 1
         RETURN
  END IF
   W( 1 ) = JVS( 1 )
   W( 37 ) = JVS( 2 )
   W( 131 ) = JVS( 3 )
  JVS( 1) = W( 1 )
  JVS( 2) = W( 37 )
  JVS( 3) = W( 131 )
  IF ( ABS( JVS( 4 )) < TINY(a) ) THEN
         IER = 2
         RETURN
  END IF
   W( 2 ) = JVS( 4 )
  JVS( 4) = W( 2 )
  IF ( ABS( JVS( 5 )) < TINY(a) ) THEN
         IER = 3
         RETURN
  END IF
   W( 3 ) = JVS( 5 )
  JVS( 5) = W( 3 )
  IF ( ABS( JVS( 6 )) < TINY(a) ) THEN
         IER = 4
         RETURN
  END IF
   W( 4 ) = JVS( 6 )
  JVS( 6) = W( 4 )
  IF ( ABS( JVS( 7 )) < TINY(a) ) THEN
         IER = 5
         RETURN
  END IF
   W( 5 ) = JVS( 7 )
  JVS( 7) = W( 5 )
  IF ( ABS( JVS( 8 )) < TINY(a) ) THEN
         IER = 6
         RETURN
  END IF
   W( 6 ) = JVS( 8 )
  JVS( 8) = W( 6 )
  IF ( ABS( JVS( 9 )) < TINY(a) ) THEN
         IER = 7
         RETURN
  END IF
   W( 7 ) = JVS( 9 )
  JVS( 9) = W( 7 )
  IF ( ABS( JVS( 10 )) < TINY(a) ) THEN
         IER = 8
         RETURN
  END IF
   W( 8 ) = JVS( 10 )
   W( 103 ) = JVS( 11 )
   W( 114 ) = JVS( 12 )
   W( 127 ) = JVS( 13 )
  JVS( 10) = W( 8 )
  JVS( 11) = W( 103 )
  JVS( 12) = W( 114 )
  JVS( 13) = W( 127 )
  IF ( ABS( JVS( 14 )) < TINY(a) ) THEN
         IER = 9
         RETURN
  END IF
   W( 9 ) = JVS( 14 )
   W( 114 ) = JVS( 15 )
   W( 115 ) = JVS( 16 )
   W( 119 ) = JVS( 17 )
   W( 125 ) = JVS( 18 )
   W( 127 ) = JVS( 19 )
   W( 130 ) = JVS( 20 )
   W( 132 ) = JVS( 21 )
  JVS( 14) = W( 9 )
  JVS( 15) = W( 114 )
  JVS( 16) = W( 115 )
  JVS( 17) = W( 119 )
  JVS( 18) = W( 125 )
  JVS( 19) = W( 127 )
  JVS( 20) = W( 130 )
  JVS( 21) = W( 132 )
  IF ( ABS( JVS( 22 )) < TINY(a) ) THEN
         IER = 10
         RETURN
  END IF
   W( 10 ) = JVS( 22 )
   W( 89 ) = JVS( 23 )
   W( 91 ) = JVS( 24 )
   W( 105 ) = JVS( 25 )
   W( 131 ) = JVS( 26 )
   W( 136 ) = JVS( 27 )
  JVS( 22) = W( 10 )
  JVS( 23) = W( 89 )
  JVS( 24) = W( 91 )
  JVS( 25) = W( 105 )
  JVS( 26) = W( 131 )
  JVS( 27) = W( 136 )
  IF ( ABS( JVS( 28 )) < TINY(a) ) THEN
         IER = 11
         RETURN
  END IF
   W( 11 ) = JVS( 28 )
   W( 89 ) = JVS( 29 )
   W( 91 ) = JVS( 30 )
   W( 105 ) = JVS( 31 )
   W( 131 ) = JVS( 32 )
   W( 136 ) = JVS( 33 )
  JVS( 28) = W( 11 )
  JVS( 29) = W( 89 )
  JVS( 30) = W( 91 )
  JVS( 31) = W( 105 )
  JVS( 32) = W( 131 )
  JVS( 33) = W( 136 )
  IF ( ABS( JVS( 34 )) < TINY(a) ) THEN
         IER = 12
         RETURN
  END IF
   W( 12 ) = JVS( 34 )
   W( 109 ) = JVS( 35 )
   W( 131 ) = JVS( 36 )
  JVS( 34) = W( 12 )
  JVS( 35) = W( 109 )
  JVS( 36) = W( 131 )
  IF ( ABS( JVS( 37 )) < TINY(a) ) THEN
         IER = 13
         RETURN
  END IF
   W( 13 ) = JVS( 37 )
   W( 114 ) = JVS( 38 )
   W( 119 ) = JVS( 39 )
   W( 127 ) = JVS( 40 )
   W( 131 ) = JVS( 41 )
   W( 136 ) = JVS( 42 )
  JVS( 37) = W( 13 )
  JVS( 38) = W( 114 )
  JVS( 39) = W( 119 )
  JVS( 40) = W( 127 )
  JVS( 41) = W( 131 )
  JVS( 42) = W( 136 )
  IF ( ABS( JVS( 43 )) < TINY(a) ) THEN
         IER = 14
         RETURN
  END IF
   W( 14 ) = JVS( 43 )
   W( 94 ) = JVS( 44 )
   W( 127 ) = JVS( 45 )
   W( 131 ) = JVS( 46 )
  JVS( 43) = W( 14 )
  JVS( 44) = W( 94 )
  JVS( 45) = W( 127 )
  JVS( 46) = W( 131 )
  IF ( ABS( JVS( 47 )) < TINY(a) ) THEN
         IER = 15
         RETURN
  END IF
   W( 15 ) = JVS( 47 )
   W( 94 ) = JVS( 48 )
   W( 127 ) = JVS( 49 )
   W( 131 ) = JVS( 50 )
  JVS( 47) = W( 15 )
  JVS( 48) = W( 94 )
  JVS( 49) = W( 127 )
  JVS( 50) = W( 131 )
  IF ( ABS( JVS( 51 )) < TINY(a) ) THEN
         IER = 16
         RETURN
  END IF
   W( 16 ) = JVS( 51 )
   W( 95 ) = JVS( 52 )
   W( 127 ) = JVS( 53 )
   W( 131 ) = JVS( 54 )
  JVS( 51) = W( 16 )
  JVS( 52) = W( 95 )
  JVS( 53) = W( 127 )
  JVS( 54) = W( 131 )
  IF ( ABS( JVS( 55 )) < TINY(a) ) THEN
         IER = 17
         RETURN
  END IF
   W( 17 ) = JVS( 55 )
   W( 95 ) = JVS( 56 )
   W( 127 ) = JVS( 57 )
   W( 131 ) = JVS( 58 )
  JVS( 55) = W( 17 )
  JVS( 56) = W( 95 )
  JVS( 57) = W( 127 )
  JVS( 58) = W( 131 )
  IF ( ABS( JVS( 59 )) < TINY(a) ) THEN
         IER = 18
         RETURN
  END IF
   W( 18 ) = JVS( 59 )
   W( 37 ) = JVS( 60 )
   W( 48 ) = JVS( 61 )
   W( 88 ) = JVS( 62 )
   W( 89 ) = JVS( 63 )
   W( 91 ) = JVS( 64 )
   W( 93 ) = JVS( 65 )
   W( 96 ) = JVS( 66 )
   W( 98 ) = JVS( 67 )
   W( 100 ) = JVS( 68 )
   W( 101 ) = JVS( 69 )
   W( 102 ) = JVS( 70 )
   W( 103 ) = JVS( 71 )
   W( 104 ) = JVS( 72 )
   W( 105 ) = JVS( 73 )
   W( 107 ) = JVS( 74 )
   W( 108 ) = JVS( 75 )
   W( 109 ) = JVS( 76 )
   W( 110 ) = JVS( 77 )
   W( 114 ) = JVS( 78 )
   W( 115 ) = JVS( 79 )
   W( 116 ) = JVS( 80 )
   W( 118 ) = JVS( 81 )
   W( 119 ) = JVS( 82 )
   W( 121 ) = JVS( 83 )
   W( 123 ) = JVS( 84 )
   W( 125 ) = JVS( 85 )
   W( 126 ) = JVS( 86 )
   W( 127 ) = JVS( 87 )
   W( 129 ) = JVS( 88 )
   W( 131 ) = JVS( 89 )
   W( 132 ) = JVS( 90 )
   W( 133 ) = JVS( 91 )
   W( 135 ) = JVS( 92 )
   W( 136 ) = JVS( 93 )
  JVS( 59) = W( 18 )
  JVS( 60) = W( 37 )
  JVS( 61) = W( 48 )
  JVS( 62) = W( 88 )
  JVS( 63) = W( 89 )
  JVS( 64) = W( 91 )
  JVS( 65) = W( 93 )
  JVS( 66) = W( 96 )
  JVS( 67) = W( 98 )
  JVS( 68) = W( 100 )
  JVS( 69) = W( 101 )
  JVS( 70) = W( 102 )
  JVS( 71) = W( 103 )
  JVS( 72) = W( 104 )
  JVS( 73) = W( 105 )
  JVS( 74) = W( 107 )
  JVS( 75) = W( 108 )
  JVS( 76) = W( 109 )
  JVS( 77) = W( 110 )
  JVS( 78) = W( 114 )
  JVS( 79) = W( 115 )
  JVS( 80) = W( 116 )
  JVS( 81) = W( 118 )
  JVS( 82) = W( 119 )
  JVS( 83) = W( 121 )
  JVS( 84) = W( 123 )
  JVS( 85) = W( 125 )
  JVS( 86) = W( 126 )
  JVS( 87) = W( 127 )
  JVS( 88) = W( 129 )
  JVS( 89) = W( 131 )
  JVS( 90) = W( 132 )
  JVS( 91) = W( 133 )
  JVS( 92) = W( 135 )
  JVS( 93) = W( 136 )
  IF ( ABS( JVS( 98 )) < TINY(a) ) THEN
         IER = 19
         RETURN
  END IF
   W( 4 ) = JVS( 94 )
   W( 5 ) = JVS( 95 )
   W( 6 ) = JVS( 96 )
   W( 7 ) = JVS( 97 )
   W( 19 ) = JVS( 98 )
   W( 20 ) = JVS( 99 )
   W( 21 ) = JVS( 100 )
   W( 22 ) = JVS( 101 )
   W( 23 ) = JVS( 102 )
   W( 24 ) = JVS( 103 )
   W( 25 ) = JVS( 104 )
   W( 26 ) = JVS( 105 )
   W( 27 ) = JVS( 106 )
   W( 28 ) = JVS( 107 )
   W( 29 ) = JVS( 108 )
   W( 30 ) = JVS( 109 )
   W( 31 ) = JVS( 110 )
   W( 32 ) = JVS( 111 )
   W( 33 ) = JVS( 112 )
   W( 34 ) = JVS( 113 )
   W( 35 ) = JVS( 114 )
   W( 38 ) = JVS( 115 )
   W( 39 ) = JVS( 116 )
   W( 40 ) = JVS( 117 )
   W( 41 ) = JVS( 118 )
   W( 42 ) = JVS( 119 )
   W( 43 ) = JVS( 120 )
   W( 44 ) = JVS( 121 )
   W( 45 ) = JVS( 122 )
   W( 46 ) = JVS( 123 )
   W( 47 ) = JVS( 124 )
   W( 49 ) = JVS( 125 )
   W( 50 ) = JVS( 126 )
   W( 51 ) = JVS( 127 )
   W( 52 ) = JVS( 128 )
   W( 53 ) = JVS( 129 )
   W( 54 ) = JVS( 130 )
   W( 55 ) = JVS( 131 )
   W( 56 ) = JVS( 132 )
   W( 57 ) = JVS( 133 )
   W( 58 ) = JVS( 134 )
   W( 59 ) = JVS( 135 )
   W( 60 ) = JVS( 136 )
   W( 61 ) = JVS( 137 )
   W( 62 ) = JVS( 138 )
   W( 63 ) = JVS( 139 )
   W( 64 ) = JVS( 140 )
   W( 65 ) = JVS( 141 )
   W( 66 ) = JVS( 142 )
   W( 67 ) = JVS( 143 )
   W( 68 ) = JVS( 144 )
   W( 69 ) = JVS( 145 )
   W( 70 ) = JVS( 146 )
   W( 71 ) = JVS( 147 )
   W( 72 ) = JVS( 148 )
   W( 73 ) = JVS( 149 )
   W( 74 ) = JVS( 150 )
   W( 75 ) = JVS( 151 )
   W( 76 ) = JVS( 152 )
   W( 77 ) = JVS( 153 )
   W( 78 ) = JVS( 154 )
   W( 79 ) = JVS( 155 )
   W( 80 ) = JVS( 156 )
   W( 81 ) = JVS( 157 )
   W( 82 ) = JVS( 158 )
   W( 83 ) = JVS( 159 )
   W( 84 ) = JVS( 160 )
   W( 85 ) = JVS( 161 )
   W( 86 ) = JVS( 162 )
   W( 94 ) = JVS( 163 )
   W( 95 ) = JVS( 164 )
   W( 131 ) = JVS( 165 )
  a = -W( 4 ) / JVS( 6 )
  W( 4 ) = -a
  a = -W( 5 ) / JVS( 7 )
  W( 5 ) = -a
  a = -W( 6 ) / JVS( 8 )
  W( 6 ) = -a
  a = -W( 7 ) / JVS( 9 )
  W( 7 ) = -a
  JVS( 94) = W( 4 )
  JVS( 95) = W( 5 )
  JVS( 96) = W( 6 )
  JVS( 97) = W( 7 )
  JVS( 98) = W( 19 )
  JVS( 99) = W( 20 )
  JVS( 100) = W( 21 )
  JVS( 101) = W( 22 )
  JVS( 102) = W( 23 )
  JVS( 103) = W( 24 )
  JVS( 104) = W( 25 )
  JVS( 105) = W( 26 )
  JVS( 106) = W( 27 )
  JVS( 107) = W( 28 )
  JVS( 108) = W( 29 )
  JVS( 109) = W( 30 )
  JVS( 110) = W( 31 )
  JVS( 111) = W( 32 )
  JVS( 112) = W( 33 )
  JVS( 113) = W( 34 )
  JVS( 114) = W( 35 )
  JVS( 115) = W( 38 )
  JVS( 116) = W( 39 )
  JVS( 117) = W( 40 )
  JVS( 118) = W( 41 )
  JVS( 119) = W( 42 )
  JVS( 120) = W( 43 )
  JVS( 121) = W( 44 )
  JVS( 122) = W( 45 )
  JVS( 123) = W( 46 )
  JVS( 124) = W( 47 )
  JVS( 125) = W( 49 )
  JVS( 126) = W( 50 )
  JVS( 127) = W( 51 )
  JVS( 128) = W( 52 )
  JVS( 129) = W( 53 )
  JVS( 130) = W( 54 )
  JVS( 131) = W( 55 )
  JVS( 132) = W( 56 )
  JVS( 133) = W( 57 )
  JVS( 134) = W( 58 )
  JVS( 135) = W( 59 )
  JVS( 136) = W( 60 )
  JVS( 137) = W( 61 )
  JVS( 138) = W( 62 )
  JVS( 139) = W( 63 )
  JVS( 140) = W( 64 )
  JVS( 141) = W( 65 )
  JVS( 142) = W( 66 )
  JVS( 143) = W( 67 )
  JVS( 144) = W( 68 )
  JVS( 145) = W( 69 )
  JVS( 146) = W( 70 )
  JVS( 147) = W( 71 )
  JVS( 148) = W( 72 )
  JVS( 149) = W( 73 )
  JVS( 150) = W( 74 )
  JVS( 151) = W( 75 )
  JVS( 152) = W( 76 )
  JVS( 153) = W( 77 )
  JVS( 154) = W( 78 )
  JVS( 155) = W( 79 )
  JVS( 156) = W( 80 )
  JVS( 157) = W( 81 )
  JVS( 158) = W( 82 )
  JVS( 159) = W( 83 )
  JVS( 160) = W( 84 )
  JVS( 161) = W( 85 )
  JVS( 162) = W( 86 )
  JVS( 163) = W( 94 )
  JVS( 164) = W( 95 )
  JVS( 165) = W( 131 )
  IF ( ABS( JVS( 166 )) < TINY(a) ) THEN
         IER = 20
         RETURN
  END IF
   W( 20 ) = JVS( 166 )
   W( 131 ) = JVS( 167 )
  JVS( 166) = W( 20 )
  JVS( 167) = W( 131 )
  IF ( ABS( JVS( 168 )) < TINY(a) ) THEN
         IER = 21
         RETURN
  END IF
   W( 21 ) = JVS( 168 )
   W( 131 ) = JVS( 169 )
  JVS( 168) = W( 21 )
  JVS( 169) = W( 131 )
  IF ( ABS( JVS( 170 )) < TINY(a) ) THEN
         IER = 22
         RETURN
  END IF
   W( 22 ) = JVS( 170 )
   W( 131 ) = JVS( 171 )
  JVS( 170) = W( 22 )
  JVS( 171) = W( 131 )
  IF ( ABS( JVS( 172 )) < TINY(a) ) THEN
         IER = 23
         RETURN
  END IF
   W( 23 ) = JVS( 172 )
   W( 131 ) = JVS( 173 )
  JVS( 172) = W( 23 )
  JVS( 173) = W( 131 )
  IF ( ABS( JVS( 174 )) < TINY(a) ) THEN
         IER = 24
         RETURN
  END IF
   W( 24 ) = JVS( 174 )
   W( 131 ) = JVS( 175 )
  JVS( 174) = W( 24 )
  JVS( 175) = W( 131 )
  IF ( ABS( JVS( 176 )) < TINY(a) ) THEN
         IER = 25
         RETURN
  END IF
   W( 25 ) = JVS( 176 )
   W( 131 ) = JVS( 177 )
  JVS( 176) = W( 25 )
  JVS( 177) = W( 131 )
  IF ( ABS( JVS( 178 )) < TINY(a) ) THEN
         IER = 26
         RETURN
  END IF
   W( 26 ) = JVS( 178 )
   W( 131 ) = JVS( 179 )
  JVS( 178) = W( 26 )
  JVS( 179) = W( 131 )
  IF ( ABS( JVS( 180 )) < TINY(a) ) THEN
         IER = 27
         RETURN
  END IF
   W( 27 ) = JVS( 180 )
   W( 131 ) = JVS( 181 )
  JVS( 180) = W( 27 )
  JVS( 181) = W( 131 )
  IF ( ABS( JVS( 182 )) < TINY(a) ) THEN
         IER = 28
         RETURN
  END IF
   W( 28 ) = JVS( 182 )
   W( 131 ) = JVS( 183 )
  JVS( 182) = W( 28 )
  JVS( 183) = W( 131 )
  IF ( ABS( JVS( 184 )) < TINY(a) ) THEN
         IER = 29
         RETURN
  END IF
   W( 29 ) = JVS( 184 )
   W( 131 ) = JVS( 185 )
  JVS( 184) = W( 29 )
  JVS( 185) = W( 131 )
  IF ( ABS( JVS( 186 )) < TINY(a) ) THEN
         IER = 30
         RETURN
  END IF
   W( 30 ) = JVS( 186 )
   W( 131 ) = JVS( 187 )
  JVS( 186) = W( 30 )
  JVS( 187) = W( 131 )
  IF ( ABS( JVS( 188 )) < TINY(a) ) THEN
         IER = 31
         RETURN
  END IF
   W( 31 ) = JVS( 188 )
   W( 131 ) = JVS( 189 )
  JVS( 188) = W( 31 )
  JVS( 189) = W( 131 )
  IF ( ABS( JVS( 190 )) < TINY(a) ) THEN
         IER = 32
         RETURN
  END IF
   W( 32 ) = JVS( 190 )
   W( 131 ) = JVS( 191 )
  JVS( 190) = W( 32 )
  JVS( 191) = W( 131 )
  IF ( ABS( JVS( 192 )) < TINY(a) ) THEN
         IER = 33
         RETURN
  END IF
   W( 33 ) = JVS( 192 )
   W( 131 ) = JVS( 193 )
  JVS( 192) = W( 33 )
  JVS( 193) = W( 131 )
  IF ( ABS( JVS( 194 )) < TINY(a) ) THEN
         IER = 34
         RETURN
  END IF
   W( 34 ) = JVS( 194 )
   W( 131 ) = JVS( 195 )
  JVS( 194) = W( 34 )
  JVS( 195) = W( 131 )
  IF ( ABS( JVS( 196 )) < TINY(a) ) THEN
         IER = 35
         RETURN
  END IF
   W( 35 ) = JVS( 196 )
   W( 131 ) = JVS( 197 )
  JVS( 196) = W( 35 )
  JVS( 197) = W( 131 )
  IF ( ABS( JVS( 198 )) < TINY(a) ) THEN
         IER = 36
         RETURN
  END IF
   W( 36 ) = JVS( 198 )
   W( 127 ) = JVS( 199 )
  JVS( 198) = W( 36 )
  JVS( 199) = W( 127 )
  IF ( ABS( JVS( 200 )) < TINY(a) ) THEN
         IER = 37
         RETURN
  END IF
   W( 37 ) = JVS( 200 )
   W( 131 ) = JVS( 201 )
  JVS( 200) = W( 37 )
  JVS( 201) = W( 131 )
  IF ( ABS( JVS( 202 )) < TINY(a) ) THEN
         IER = 38
         RETURN
  END IF
   W( 38 ) = JVS( 202 )
   W( 131 ) = JVS( 203 )
  JVS( 202) = W( 38 )
  JVS( 203) = W( 131 )
  IF ( ABS( JVS( 205 )) < TINY(a) ) THEN
         IER = 39
         RETURN
  END IF
   W( 38 ) = JVS( 204 )
   W( 39 ) = JVS( 205 )
   W( 60 ) = JVS( 206 )
   W( 131 ) = JVS( 207 )
  a = -W( 38 ) / JVS( 202 )
  W( 38 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 203 )
  JVS( 204) = W( 38 )
  JVS( 205) = W( 39 )
  JVS( 206) = W( 60 )
  JVS( 207) = W( 131 )
  IF ( ABS( JVS( 208 )) < TINY(a) ) THEN
         IER = 40
         RETURN
  END IF
   W( 40 ) = JVS( 208 )
   W( 131 ) = JVS( 209 )
  JVS( 208) = W( 40 )
  JVS( 209) = W( 131 )
  IF ( ABS( JVS( 210 )) < TINY(a) ) THEN
         IER = 41
         RETURN
  END IF
   W( 41 ) = JVS( 210 )
   W( 131 ) = JVS( 211 )
  JVS( 210) = W( 41 )
  JVS( 211) = W( 131 )
  IF ( ABS( JVS( 212 )) < TINY(a) ) THEN
         IER = 42
         RETURN
  END IF
   W( 42 ) = JVS( 212 )
   W( 131 ) = JVS( 213 )
  JVS( 212) = W( 42 )
  JVS( 213) = W( 131 )
  IF ( ABS( JVS( 214 )) < TINY(a) ) THEN
         IER = 43
         RETURN
  END IF
   W( 43 ) = JVS( 214 )
   W( 131 ) = JVS( 215 )
  JVS( 214) = W( 43 )
  JVS( 215) = W( 131 )
  IF ( ABS( JVS( 216 )) < TINY(a) ) THEN
         IER = 44
         RETURN
  END IF
   W( 44 ) = JVS( 216 )
   W( 131 ) = JVS( 217 )
  JVS( 216) = W( 44 )
  JVS( 217) = W( 131 )
  IF ( ABS( JVS( 218 )) < TINY(a) ) THEN
         IER = 45
         RETURN
  END IF
   W( 45 ) = JVS( 218 )
   W( 131 ) = JVS( 219 )
  JVS( 218) = W( 45 )
  JVS( 219) = W( 131 )
  IF ( ABS( JVS( 220 )) < TINY(a) ) THEN
         IER = 46
         RETURN
  END IF
   W( 46 ) = JVS( 220 )
   W( 131 ) = JVS( 221 )
  JVS( 220) = W( 46 )
  JVS( 221) = W( 131 )
  IF ( ABS( JVS( 224 )) < TINY(a) ) THEN
         IER = 47
         RETURN
  END IF
   W( 27 ) = JVS( 222 )
   W( 46 ) = JVS( 223 )
   W( 47 ) = JVS( 224 )
   W( 131 ) = JVS( 225 )
  a = -W( 27 ) / JVS( 180 )
  W( 27 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 181 )
  a = -W( 46 ) / JVS( 220 )
  W( 46 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 221 )
  JVS( 222) = W( 27 )
  JVS( 223) = W( 46 )
  JVS( 224) = W( 47 )
  JVS( 225) = W( 131 )
  IF ( ABS( JVS( 226 )) < TINY(a) ) THEN
         IER = 48
         RETURN
  END IF
   W( 48 ) = JVS( 226 )
   W( 131 ) = JVS( 227 )
  JVS( 226) = W( 48 )
  JVS( 227) = W( 131 )
  IF ( ABS( JVS( 229 )) < TINY(a) ) THEN
         IER = 49
         RETURN
  END IF
   W( 46 ) = JVS( 228 )
   W( 49 ) = JVS( 229 )
   W( 131 ) = JVS( 230 )
  a = -W( 46 ) / JVS( 220 )
  W( 46 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 221 )
  JVS( 228) = W( 46 )
  JVS( 229) = W( 49 )
  JVS( 230) = W( 131 )
  IF ( ABS( JVS( 235 )) < TINY(a) ) THEN
         IER = 50
         RETURN
  END IF
   W( 26 ) = JVS( 231 )
   W( 45 ) = JVS( 232 )
   W( 47 ) = JVS( 233 )
   W( 49 ) = JVS( 234 )
   W( 50 ) = JVS( 235 )
   W( 131 ) = JVS( 236 )
  a = -W( 26 ) / JVS( 178 )
  W( 26 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 179 )
  a = -W( 45 ) / JVS( 218 )
  W( 45 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 219 )
  a = -W( 47 ) / JVS( 224 )
  W( 47 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 225 )
  a = -W( 49 ) / JVS( 229 )
  W( 49 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 230 )
  JVS( 231) = W( 26 )
  JVS( 232) = W( 45 )
  JVS( 233) = W( 47 )
  JVS( 234) = W( 49 )
  JVS( 235) = W( 50 )
  JVS( 236) = W( 131 )
  IF ( ABS( JVS( 239 )) < TINY(a) ) THEN
         IER = 51
         RETURN
  END IF
   W( 45 ) = JVS( 237 )
   W( 49 ) = JVS( 238 )
   W( 51 ) = JVS( 239 )
   W( 131 ) = JVS( 240 )
  a = -W( 45 ) / JVS( 218 )
  W( 45 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 219 )
  a = -W( 49 ) / JVS( 229 )
  W( 49 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 230 )
  JVS( 237) = W( 45 )
  JVS( 238) = W( 49 )
  JVS( 239) = W( 51 )
  JVS( 240) = W( 131 )
  IF ( ABS( JVS( 245 )) < TINY(a) ) THEN
         IER = 52
         RETURN
  END IF
   W( 25 ) = JVS( 241 )
   W( 44 ) = JVS( 242 )
   W( 50 ) = JVS( 243 )
   W( 51 ) = JVS( 244 )
   W( 52 ) = JVS( 245 )
   W( 131 ) = JVS( 246 )
  a = -W( 25 ) / JVS( 176 )
  W( 25 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 177 )
  a = -W( 44 ) / JVS( 216 )
  W( 44 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 217 )
  a = -W( 50 ) / JVS( 235 )
  W( 50 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 236 )
  a = -W( 51 ) / JVS( 239 )
  W( 51 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 240 )
  JVS( 241) = W( 25 )
  JVS( 242) = W( 44 )
  JVS( 243) = W( 50 )
  JVS( 244) = W( 51 )
  JVS( 245) = W( 52 )
  JVS( 246) = W( 131 )
  IF ( ABS( JVS( 249 )) < TINY(a) ) THEN
         IER = 53
         RETURN
  END IF
   W( 44 ) = JVS( 247 )
   W( 51 ) = JVS( 248 )
   W( 53 ) = JVS( 249 )
   W( 131 ) = JVS( 250 )
  a = -W( 44 ) / JVS( 216 )
  W( 44 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 217 )
  a = -W( 51 ) / JVS( 239 )
  W( 51 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 240 )
  JVS( 247) = W( 44 )
  JVS( 248) = W( 51 )
  JVS( 249) = W( 53 )
  JVS( 250) = W( 131 )
  IF ( ABS( JVS( 255 )) < TINY(a) ) THEN
         IER = 54
         RETURN
  END IF
   W( 24 ) = JVS( 251 )
   W( 43 ) = JVS( 252 )
   W( 52 ) = JVS( 253 )
   W( 53 ) = JVS( 254 )
   W( 54 ) = JVS( 255 )
   W( 131 ) = JVS( 256 )
  a = -W( 24 ) / JVS( 174 )
  W( 24 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 175 )
  a = -W( 43 ) / JVS( 214 )
  W( 43 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 215 )
  a = -W( 52 ) / JVS( 245 )
  W( 52 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 246 )
  a = -W( 53 ) / JVS( 249 )
  W( 53 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 250 )
  JVS( 251) = W( 24 )
  JVS( 252) = W( 43 )
  JVS( 253) = W( 52 )
  JVS( 254) = W( 53 )
  JVS( 255) = W( 54 )
  JVS( 256) = W( 131 )
  IF ( ABS( JVS( 259 )) < TINY(a) ) THEN
         IER = 55
         RETURN
  END IF
   W( 43 ) = JVS( 257 )
   W( 53 ) = JVS( 258 )
   W( 55 ) = JVS( 259 )
   W( 131 ) = JVS( 260 )
  a = -W( 43 ) / JVS( 214 )
  W( 43 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 215 )
  a = -W( 53 ) / JVS( 249 )
  W( 53 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 250 )
  JVS( 257) = W( 43 )
  JVS( 258) = W( 53 )
  JVS( 259) = W( 55 )
  JVS( 260) = W( 131 )
  IF ( ABS( JVS( 265 )) < TINY(a) ) THEN
         IER = 56
         RETURN
  END IF
   W( 23 ) = JVS( 261 )
   W( 42 ) = JVS( 262 )
   W( 54 ) = JVS( 263 )
   W( 55 ) = JVS( 264 )
   W( 56 ) = JVS( 265 )
   W( 131 ) = JVS( 266 )
  a = -W( 23 ) / JVS( 172 )
  W( 23 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 173 )
  a = -W( 42 ) / JVS( 212 )
  W( 42 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 213 )
  a = -W( 54 ) / JVS( 255 )
  W( 54 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 256 )
  a = -W( 55 ) / JVS( 259 )
  W( 55 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 260 )
  JVS( 261) = W( 23 )
  JVS( 262) = W( 42 )
  JVS( 263) = W( 54 )
  JVS( 264) = W( 55 )
  JVS( 265) = W( 56 )
  JVS( 266) = W( 131 )
  IF ( ABS( JVS( 269 )) < TINY(a) ) THEN
         IER = 57
         RETURN
  END IF
   W( 42 ) = JVS( 267 )
   W( 55 ) = JVS( 268 )
   W( 57 ) = JVS( 269 )
   W( 131 ) = JVS( 270 )
  a = -W( 42 ) / JVS( 212 )
  W( 42 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 213 )
  a = -W( 55 ) / JVS( 259 )
  W( 55 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 260 )
  JVS( 267) = W( 42 )
  JVS( 268) = W( 55 )
  JVS( 269) = W( 57 )
  JVS( 270) = W( 131 )
  IF ( ABS( JVS( 275 )) < TINY(a) ) THEN
         IER = 58
         RETURN
  END IF
   W( 22 ) = JVS( 271 )
   W( 41 ) = JVS( 272 )
   W( 56 ) = JVS( 273 )
   W( 57 ) = JVS( 274 )
   W( 58 ) = JVS( 275 )
   W( 131 ) = JVS( 276 )
  a = -W( 22 ) / JVS( 170 )
  W( 22 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 171 )
  a = -W( 41 ) / JVS( 210 )
  W( 41 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 211 )
  a = -W( 56 ) / JVS( 265 )
  W( 56 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 266 )
  a = -W( 57 ) / JVS( 269 )
  W( 57 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 270 )
  JVS( 271) = W( 22 )
  JVS( 272) = W( 41 )
  JVS( 273) = W( 56 )
  JVS( 274) = W( 57 )
  JVS( 275) = W( 58 )
  JVS( 276) = W( 131 )
  IF ( ABS( JVS( 279 )) < TINY(a) ) THEN
         IER = 59
         RETURN
  END IF
   W( 41 ) = JVS( 277 )
   W( 57 ) = JVS( 278 )
   W( 59 ) = JVS( 279 )
   W( 131 ) = JVS( 280 )
  a = -W( 41 ) / JVS( 210 )
  W( 41 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 211 )
  a = -W( 57 ) / JVS( 269 )
  W( 57 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 270 )
  JVS( 277) = W( 41 )
  JVS( 278) = W( 57 )
  JVS( 279) = W( 59 )
  JVS( 280) = W( 131 )
  IF ( ABS( JVS( 283 )) < TINY(a) ) THEN
         IER = 60
         RETURN
  END IF
   W( 40 ) = JVS( 281 )
   W( 59 ) = JVS( 282 )
   W( 60 ) = JVS( 283 )
   W( 131 ) = JVS( 284 )
  a = -W( 40 ) / JVS( 208 )
  W( 40 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 209 )
  a = -W( 59 ) / JVS( 279 )
  W( 59 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 280 )
  JVS( 281) = W( 40 )
  JVS( 282) = W( 59 )
  JVS( 283) = W( 60 )
  JVS( 284) = W( 131 )
  IF ( ABS( JVS( 288 )) < TINY(a) ) THEN
         IER = 61
         RETURN
  END IF
   W( 20 ) = JVS( 285 )
   W( 38 ) = JVS( 286 )
   W( 60 ) = JVS( 287 )
   W( 61 ) = JVS( 288 )
   W( 62 ) = JVS( 289 )
   W( 131 ) = JVS( 290 )
  a = -W( 20 ) / JVS( 166 )
  W( 20 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 167 )
  a = -W( 38 ) / JVS( 202 )
  W( 38 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 203 )
  a = -W( 60 ) / JVS( 283 )
  W( 60 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 284 )
  JVS( 285) = W( 20 )
  JVS( 286) = W( 38 )
  JVS( 287) = W( 60 )
  JVS( 288) = W( 61 )
  JVS( 289) = W( 62 )
  JVS( 290) = W( 131 )
  IF ( ABS( JVS( 295 )) < TINY(a) ) THEN
         IER = 62
         RETURN
  END IF
   W( 21 ) = JVS( 291 )
   W( 40 ) = JVS( 292 )
   W( 58 ) = JVS( 293 )
   W( 59 ) = JVS( 294 )
   W( 62 ) = JVS( 295 )
   W( 131 ) = JVS( 296 )
  a = -W( 21 ) / JVS( 168 )
  W( 21 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 169 )
  a = -W( 40 ) / JVS( 208 )
  W( 40 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 209 )
  a = -W( 58 ) / JVS( 275 )
  W( 58 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 276 )
  a = -W( 59 ) / JVS( 279 )
  W( 59 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 280 )
  JVS( 291) = W( 21 )
  JVS( 292) = W( 40 )
  JVS( 293) = W( 58 )
  JVS( 294) = W( 59 )
  JVS( 295) = W( 62 )
  JVS( 296) = W( 131 )
  IF ( ABS( JVS( 297 )) < TINY(a) ) THEN
         IER = 63
         RETURN
  END IF
   W( 63 ) = JVS( 297 )
   W( 131 ) = JVS( 298 )
  JVS( 297) = W( 63 )
  JVS( 298) = W( 131 )
  IF ( ABS( JVS( 300 )) < TINY(a) ) THEN
         IER = 64
         RETURN
  END IF
   W( 63 ) = JVS( 299 )
   W( 64 ) = JVS( 300 )
   W( 84 ) = JVS( 301 )
   W( 131 ) = JVS( 302 )
  a = -W( 63 ) / JVS( 297 )
  W( 63 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 298 )
  JVS( 299) = W( 63 )
  JVS( 300) = W( 64 )
  JVS( 301) = W( 84 )
  JVS( 302) = W( 131 )
  IF ( ABS( JVS( 303 )) < TINY(a) ) THEN
         IER = 65
         RETURN
  END IF
   W( 65 ) = JVS( 303 )
   W( 131 ) = JVS( 304 )
  JVS( 303) = W( 65 )
  JVS( 304) = W( 131 )
  IF ( ABS( JVS( 305 )) < TINY(a) ) THEN
         IER = 66
         RETURN
  END IF
   W( 66 ) = JVS( 305 )
   W( 131 ) = JVS( 306 )
  JVS( 305) = W( 66 )
  JVS( 306) = W( 131 )
  IF ( ABS( JVS( 307 )) < TINY(a) ) THEN
         IER = 67
         RETURN
  END IF
   W( 67 ) = JVS( 307 )
   W( 131 ) = JVS( 308 )
  JVS( 307) = W( 67 )
  JVS( 308) = W( 131 )
  IF ( ABS( JVS( 309 )) < TINY(a) ) THEN
         IER = 68
         RETURN
  END IF
   W( 68 ) = JVS( 309 )
   W( 131 ) = JVS( 310 )
  JVS( 309) = W( 68 )
  JVS( 310) = W( 131 )
  IF ( ABS( JVS( 311 )) < TINY(a) ) THEN
         IER = 69
         RETURN
  END IF
   W( 69 ) = JVS( 311 )
   W( 131 ) = JVS( 312 )
  JVS( 311) = W( 69 )
  JVS( 312) = W( 131 )
  IF ( ABS( JVS( 313 )) < TINY(a) ) THEN
         IER = 70
         RETURN
  END IF
   W( 70 ) = JVS( 313 )
   W( 131 ) = JVS( 314 )
  JVS( 313) = W( 70 )
  JVS( 314) = W( 131 )
  IF ( ABS( JVS( 315 )) < TINY(a) ) THEN
         IER = 71
         RETURN
  END IF
   W( 71 ) = JVS( 315 )
   W( 131 ) = JVS( 316 )
  JVS( 315) = W( 71 )
  JVS( 316) = W( 131 )
  IF ( ABS( JVS( 319 )) < TINY(a) ) THEN
         IER = 72
         RETURN
  END IF
   W( 35 ) = JVS( 317 )
   W( 71 ) = JVS( 318 )
   W( 72 ) = JVS( 319 )
   W( 131 ) = JVS( 320 )
  a = -W( 35 ) / JVS( 196 )
  W( 35 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 197 )
  a = -W( 71 ) / JVS( 315 )
  W( 71 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 316 )
  JVS( 317) = W( 35 )
  JVS( 318) = W( 71 )
  JVS( 319) = W( 72 )
  JVS( 320) = W( 131 )
  IF ( ABS( JVS( 322 )) < TINY(a) ) THEN
         IER = 73
         RETURN
  END IF
   W( 71 ) = JVS( 321 )
   W( 73 ) = JVS( 322 )
   W( 131 ) = JVS( 323 )
  a = -W( 71 ) / JVS( 315 )
  W( 71 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 316 )
  JVS( 321) = W( 71 )
  JVS( 322) = W( 73 )
  JVS( 323) = W( 131 )
  IF ( ABS( JVS( 328 )) < TINY(a) ) THEN
         IER = 74
         RETURN
  END IF
   W( 34 ) = JVS( 324 )
   W( 70 ) = JVS( 325 )
   W( 72 ) = JVS( 326 )
   W( 73 ) = JVS( 327 )
   W( 74 ) = JVS( 328 )
   W( 131 ) = JVS( 329 )
  a = -W( 34 ) / JVS( 194 )
  W( 34 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 195 )
  a = -W( 70 ) / JVS( 313 )
  W( 70 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 314 )
  a = -W( 72 ) / JVS( 319 )
  W( 72 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 320 )
  a = -W( 73 ) / JVS( 322 )
  W( 73 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 323 )
  JVS( 324) = W( 34 )
  JVS( 325) = W( 70 )
  JVS( 326) = W( 72 )
  JVS( 327) = W( 73 )
  JVS( 328) = W( 74 )
  JVS( 329) = W( 131 )
  IF ( ABS( JVS( 332 )) < TINY(a) ) THEN
         IER = 75
         RETURN
  END IF
   W( 70 ) = JVS( 330 )
   W( 73 ) = JVS( 331 )
   W( 75 ) = JVS( 332 )
   W( 131 ) = JVS( 333 )
  a = -W( 70 ) / JVS( 313 )
  W( 70 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 314 )
  a = -W( 73 ) / JVS( 322 )
  W( 73 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 323 )
  JVS( 330) = W( 70 )
  JVS( 331) = W( 73 )
  JVS( 332) = W( 75 )
  JVS( 333) = W( 131 )
  IF ( ABS( JVS( 338 )) < TINY(a) ) THEN
         IER = 76
         RETURN
  END IF
   W( 33 ) = JVS( 334 )
   W( 69 ) = JVS( 335 )
   W( 74 ) = JVS( 336 )
   W( 75 ) = JVS( 337 )
   W( 76 ) = JVS( 338 )
   W( 131 ) = JVS( 339 )
  a = -W( 33 ) / JVS( 192 )
  W( 33 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 193 )
  a = -W( 69 ) / JVS( 311 )
  W( 69 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 312 )
  a = -W( 74 ) / JVS( 328 )
  W( 74 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 329 )
  a = -W( 75 ) / JVS( 332 )
  W( 75 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 333 )
  JVS( 334) = W( 33 )
  JVS( 335) = W( 69 )
  JVS( 336) = W( 74 )
  JVS( 337) = W( 75 )
  JVS( 338) = W( 76 )
  JVS( 339) = W( 131 )
  IF ( ABS( JVS( 342 )) < TINY(a) ) THEN
         IER = 77
         RETURN
  END IF
   W( 69 ) = JVS( 340 )
   W( 75 ) = JVS( 341 )
   W( 77 ) = JVS( 342 )
   W( 131 ) = JVS( 343 )
  a = -W( 69 ) / JVS( 311 )
  W( 69 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 312 )
  a = -W( 75 ) / JVS( 332 )
  W( 75 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 333 )
  JVS( 340) = W( 69 )
  JVS( 341) = W( 75 )
  JVS( 342) = W( 77 )
  JVS( 343) = W( 131 )
  IF ( ABS( JVS( 348 )) < TINY(a) ) THEN
         IER = 78
         RETURN
  END IF
   W( 32 ) = JVS( 344 )
   W( 68 ) = JVS( 345 )
   W( 76 ) = JVS( 346 )
   W( 77 ) = JVS( 347 )
   W( 78 ) = JVS( 348 )
   W( 131 ) = JVS( 349 )
  a = -W( 32 ) / JVS( 190 )
  W( 32 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 191 )
  a = -W( 68 ) / JVS( 309 )
  W( 68 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 310 )
  a = -W( 76 ) / JVS( 338 )
  W( 76 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 339 )
  a = -W( 77 ) / JVS( 342 )
  W( 77 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 343 )
  JVS( 344) = W( 32 )
  JVS( 345) = W( 68 )
  JVS( 346) = W( 76 )
  JVS( 347) = W( 77 )
  JVS( 348) = W( 78 )
  JVS( 349) = W( 131 )
  IF ( ABS( JVS( 352 )) < TINY(a) ) THEN
         IER = 79
         RETURN
  END IF
   W( 68 ) = JVS( 350 )
   W( 77 ) = JVS( 351 )
   W( 79 ) = JVS( 352 )
   W( 131 ) = JVS( 353 )
  a = -W( 68 ) / JVS( 309 )
  W( 68 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 310 )
  a = -W( 77 ) / JVS( 342 )
  W( 77 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 343 )
  JVS( 350) = W( 68 )
  JVS( 351) = W( 77 )
  JVS( 352) = W( 79 )
  JVS( 353) = W( 131 )
  IF ( ABS( JVS( 358 )) < TINY(a) ) THEN
         IER = 80
         RETURN
  END IF
   W( 31 ) = JVS( 354 )
   W( 67 ) = JVS( 355 )
   W( 78 ) = JVS( 356 )
   W( 79 ) = JVS( 357 )
   W( 80 ) = JVS( 358 )
   W( 131 ) = JVS( 359 )
  a = -W( 31 ) / JVS( 188 )
  W( 31 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 189 )
  a = -W( 67 ) / JVS( 307 )
  W( 67 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 308 )
  a = -W( 78 ) / JVS( 348 )
  W( 78 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 349 )
  a = -W( 79 ) / JVS( 352 )
  W( 79 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 353 )
  JVS( 354) = W( 31 )
  JVS( 355) = W( 67 )
  JVS( 356) = W( 78 )
  JVS( 357) = W( 79 )
  JVS( 358) = W( 80 )
  JVS( 359) = W( 131 )
  IF ( ABS( JVS( 362 )) < TINY(a) ) THEN
         IER = 81
         RETURN
  END IF
   W( 67 ) = JVS( 360 )
   W( 79 ) = JVS( 361 )
   W( 81 ) = JVS( 362 )
   W( 131 ) = JVS( 363 )
  a = -W( 67 ) / JVS( 307 )
  W( 67 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 308 )
  a = -W( 79 ) / JVS( 352 )
  W( 79 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 353 )
  JVS( 360) = W( 67 )
  JVS( 361) = W( 79 )
  JVS( 362) = W( 81 )
  JVS( 363) = W( 131 )
  IF ( ABS( JVS( 368 )) < TINY(a) ) THEN
         IER = 82
         RETURN
  END IF
   W( 30 ) = JVS( 364 )
   W( 66 ) = JVS( 365 )
   W( 80 ) = JVS( 366 )
   W( 81 ) = JVS( 367 )
   W( 82 ) = JVS( 368 )
   W( 131 ) = JVS( 369 )
  a = -W( 30 ) / JVS( 186 )
  W( 30 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 187 )
  a = -W( 66 ) / JVS( 305 )
  W( 66 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 306 )
  a = -W( 80 ) / JVS( 358 )
  W( 80 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 359 )
  a = -W( 81 ) / JVS( 362 )
  W( 81 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 363 )
  JVS( 364) = W( 30 )
  JVS( 365) = W( 66 )
  JVS( 366) = W( 80 )
  JVS( 367) = W( 81 )
  JVS( 368) = W( 82 )
  JVS( 369) = W( 131 )
  IF ( ABS( JVS( 372 )) < TINY(a) ) THEN
         IER = 83
         RETURN
  END IF
   W( 66 ) = JVS( 370 )
   W( 81 ) = JVS( 371 )
   W( 83 ) = JVS( 372 )
   W( 131 ) = JVS( 373 )
  a = -W( 66 ) / JVS( 305 )
  W( 66 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 306 )
  a = -W( 81 ) / JVS( 362 )
  W( 81 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 363 )
  JVS( 370) = W( 66 )
  JVS( 371) = W( 81 )
  JVS( 372) = W( 83 )
  JVS( 373) = W( 131 )
  IF ( ABS( JVS( 376 )) < TINY(a) ) THEN
         IER = 84
         RETURN
  END IF
   W( 65 ) = JVS( 374 )
   W( 83 ) = JVS( 375 )
   W( 84 ) = JVS( 376 )
   W( 131 ) = JVS( 377 )
  a = -W( 65 ) / JVS( 303 )
  W( 65 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 304 )
  a = -W( 83 ) / JVS( 372 )
  W( 83 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 373 )
  JVS( 374) = W( 65 )
  JVS( 375) = W( 83 )
  JVS( 376) = W( 84 )
  JVS( 377) = W( 131 )
  IF ( ABS( JVS( 381 )) < TINY(a) ) THEN
         IER = 85
         RETURN
  END IF
   W( 28 ) = JVS( 378 )
   W( 63 ) = JVS( 379 )
   W( 84 ) = JVS( 380 )
   W( 85 ) = JVS( 381 )
   W( 86 ) = JVS( 382 )
   W( 131 ) = JVS( 383 )
  a = -W( 28 ) / JVS( 182 )
  W( 28 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 183 )
  a = -W( 63 ) / JVS( 297 )
  W( 63 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 298 )
  a = -W( 84 ) / JVS( 376 )
  W( 84 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 377 )
  JVS( 378) = W( 28 )
  JVS( 379) = W( 63 )
  JVS( 380) = W( 84 )
  JVS( 381) = W( 85 )
  JVS( 382) = W( 86 )
  JVS( 383) = W( 131 )
  IF ( ABS( JVS( 388 )) < TINY(a) ) THEN
         IER = 86
         RETURN
  END IF
   W( 29 ) = JVS( 384 )
   W( 65 ) = JVS( 385 )
   W( 82 ) = JVS( 386 )
   W( 83 ) = JVS( 387 )
   W( 86 ) = JVS( 388 )
   W( 131 ) = JVS( 389 )
  a = -W( 29 ) / JVS( 184 )
  W( 29 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 185 )
  a = -W( 65 ) / JVS( 303 )
  W( 65 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 304 )
  a = -W( 82 ) / JVS( 368 )
  W( 82 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 369 )
  a = -W( 83 ) / JVS( 372 )
  W( 83 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 373 )
  JVS( 384) = W( 29 )
  JVS( 385) = W( 65 )
  JVS( 386) = W( 82 )
  JVS( 387) = W( 83 )
  JVS( 388) = W( 86 )
  JVS( 389) = W( 131 )
  IF ( ABS( JVS( 390 )) < TINY(a) ) THEN
         IER = 87
         RETURN
  END IF
   W( 87 ) = JVS( 390 )
   W( 129 ) = JVS( 391 )
   W( 130 ) = JVS( 392 )
  JVS( 390) = W( 87 )
  JVS( 391) = W( 129 )
  JVS( 392) = W( 130 )
  IF ( ABS( JVS( 393 )) < TINY(a) ) THEN
         IER = 88
         RETURN
  END IF
   W( 88 ) = JVS( 393 )
   W( 131 ) = JVS( 394 )
   W( 132 ) = JVS( 395 )
  JVS( 393) = W( 88 )
  JVS( 394) = W( 131 )
  JVS( 395) = W( 132 )
  IF ( ABS( JVS( 396 )) < TINY(a) ) THEN
         IER = 89
         RETURN
  END IF
   W( 89 ) = JVS( 396 )
   W( 131 ) = JVS( 397 )
  JVS( 396) = W( 89 )
  JVS( 397) = W( 131 )
  IF ( ABS( JVS( 398 )) < TINY(a) ) THEN
         IER = 90
         RETURN
  END IF
   W( 90 ) = JVS( 398 )
   W( 129 ) = JVS( 399 )
   W( 136 ) = JVS( 400 )
  JVS( 398) = W( 90 )
  JVS( 399) = W( 129 )
  JVS( 400) = W( 136 )
  IF ( ABS( JVS( 401 )) < TINY(a) ) THEN
         IER = 91
         RETURN
  END IF
   W( 91 ) = JVS( 401 )
   W( 131 ) = JVS( 402 )
  JVS( 401) = W( 91 )
  JVS( 402) = W( 131 )
  IF ( ABS( JVS( 403 )) < TINY(a) ) THEN
         IER = 92
         RETURN
  END IF
   W( 92 ) = JVS( 403 )
   W( 105 ) = JVS( 404 )
   W( 129 ) = JVS( 405 )
   W( 131 ) = JVS( 406 )
   W( 136 ) = JVS( 407 )
  JVS( 403) = W( 92 )
  JVS( 404) = W( 105 )
  JVS( 405) = W( 129 )
  JVS( 406) = W( 131 )
  JVS( 407) = W( 136 )
  IF ( ABS( JVS( 408 )) < TINY(a) ) THEN
         IER = 93
         RETURN
  END IF
   W( 93 ) = JVS( 408 )
   W( 114 ) = JVS( 409 )
   W( 119 ) = JVS( 410 )
   W( 127 ) = JVS( 411 )
   W( 131 ) = JVS( 412 )
  JVS( 408) = W( 93 )
  JVS( 409) = W( 114 )
  JVS( 410) = W( 119 )
  JVS( 411) = W( 127 )
  JVS( 412) = W( 131 )
  IF ( ABS( JVS( 413 )) < TINY(a) ) THEN
         IER = 94
         RETURN
  END IF
   W( 94 ) = JVS( 413 )
   W( 127 ) = JVS( 414 )
   W( 131 ) = JVS( 415 )
   W( 136 ) = JVS( 416 )
  JVS( 413) = W( 94 )
  JVS( 414) = W( 127 )
  JVS( 415) = W( 131 )
  JVS( 416) = W( 136 )
  IF ( ABS( JVS( 417 )) < TINY(a) ) THEN
         IER = 95
         RETURN
  END IF
   W( 95 ) = JVS( 417 )
   W( 127 ) = JVS( 418 )
   W( 131 ) = JVS( 419 )
   W( 136 ) = JVS( 420 )
  JVS( 417) = W( 95 )
  JVS( 418) = W( 127 )
  JVS( 419) = W( 131 )
  JVS( 420) = W( 136 )
  IF ( ABS( JVS( 421 )) < TINY(a) ) THEN
         IER = 96
         RETURN
  END IF
   W( 96 ) = JVS( 421 )
   W( 129 ) = JVS( 422 )
   W( 131 ) = JVS( 423 )
   W( 132 ) = JVS( 424 )
  JVS( 421) = W( 96 )
  JVS( 422) = W( 129 )
  JVS( 423) = W( 131 )
  JVS( 424) = W( 132 )
  IF ( ABS( JVS( 427 )) < TINY(a) ) THEN
         IER = 97
         RETURN
  END IF
   W( 89 ) = JVS( 425 )
   W( 91 ) = JVS( 426 )
   W( 97 ) = JVS( 427 )
   W( 131 ) = JVS( 428 )
   W( 133 ) = JVS( 429 )
  a = -W( 89 ) / JVS( 396 )
  W( 89 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 397 )
  a = -W( 91 ) / JVS( 401 )
  W( 91 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 402 )
  JVS( 425) = W( 89 )
  JVS( 426) = W( 91 )
  JVS( 427) = W( 97 )
  JVS( 428) = W( 131 )
  JVS( 429) = W( 133 )
  IF ( ABS( JVS( 430 )) < TINY(a) ) THEN
         IER = 98
         RETURN
  END IF
   W( 98 ) = JVS( 430 )
   W( 114 ) = JVS( 431 )
   W( 119 ) = JVS( 432 )
   W( 122 ) = JVS( 433 )
   W( 127 ) = JVS( 434 )
   W( 131 ) = JVS( 435 )
  JVS( 430) = W( 98 )
  JVS( 431) = W( 114 )
  JVS( 432) = W( 119 )
  JVS( 433) = W( 122 )
  JVS( 434) = W( 127 )
  JVS( 435) = W( 131 )
  IF ( ABS( JVS( 436 )) < TINY(a) ) THEN
         IER = 99
         RETURN
  END IF
   W( 99 ) = JVS( 436 )
   W( 109 ) = JVS( 437 )
   W( 114 ) = JVS( 438 )
   W( 119 ) = JVS( 439 )
   W( 120 ) = JVS( 440 )
   W( 126 ) = JVS( 441 )
   W( 127 ) = JVS( 442 )
   W( 128 ) = JVS( 443 )
   W( 131 ) = JVS( 444 )
   W( 133 ) = JVS( 445 )
   W( 135 ) = JVS( 446 )
   W( 136 ) = JVS( 447 )
  JVS( 436) = W( 99 )
  JVS( 437) = W( 109 )
  JVS( 438) = W( 114 )
  JVS( 439) = W( 119 )
  JVS( 440) = W( 120 )
  JVS( 441) = W( 126 )
  JVS( 442) = W( 127 )
  JVS( 443) = W( 128 )
  JVS( 444) = W( 131 )
  JVS( 445) = W( 133 )
  JVS( 446) = W( 135 )
  JVS( 447) = W( 136 )
  IF ( ABS( JVS( 448 )) < TINY(a) ) THEN
         IER = 100
         RETURN
  END IF
   W( 100 ) = JVS( 448 )
   W( 124 ) = JVS( 449 )
   W( 131 ) = JVS( 450 )
   W( 132 ) = JVS( 451 )
  JVS( 448) = W( 100 )
  JVS( 449) = W( 124 )
  JVS( 450) = W( 131 )
  JVS( 451) = W( 132 )
  IF ( ABS( JVS( 452 )) < TINY(a) ) THEN
         IER = 101
         RETURN
  END IF
   W( 101 ) = JVS( 452 )
   W( 122 ) = JVS( 453 )
   W( 131 ) = JVS( 454 )
   W( 132 ) = JVS( 455 )
  JVS( 452) = W( 101 )
  JVS( 453) = W( 122 )
  JVS( 454) = W( 131 )
  JVS( 455) = W( 132 )
  IF ( ABS( JVS( 456 )) < TINY(a) ) THEN
         IER = 102
         RETURN
  END IF
   W( 102 ) = JVS( 456 )
   W( 129 ) = JVS( 457 )
   W( 131 ) = JVS( 458 )
   W( 132 ) = JVS( 459 )
   W( 133 ) = JVS( 460 )
  JVS( 456) = W( 102 )
  JVS( 457) = W( 129 )
  JVS( 458) = W( 131 )
  JVS( 459) = W( 132 )
  JVS( 460) = W( 133 )
  IF ( ABS( JVS( 461 )) < TINY(a) ) THEN
         IER = 103
         RETURN
  END IF
   W( 103 ) = JVS( 461 )
   W( 127 ) = JVS( 462 )
   W( 131 ) = JVS( 463 )
  JVS( 461) = W( 103 )
  JVS( 462) = W( 127 )
  JVS( 463) = W( 131 )
  IF ( ABS( JVS( 464 )) < TINY(a) ) THEN
         IER = 104
         RETURN
  END IF
   W( 104 ) = JVS( 464 )
   W( 114 ) = JVS( 465 )
   W( 119 ) = JVS( 466 )
   W( 124 ) = JVS( 467 )
   W( 127 ) = JVS( 468 )
   W( 131 ) = JVS( 469 )
  JVS( 464) = W( 104 )
  JVS( 465) = W( 114 )
  JVS( 466) = W( 119 )
  JVS( 467) = W( 124 )
  JVS( 468) = W( 127 )
  JVS( 469) = W( 131 )
  IF ( ABS( JVS( 472 )) < TINY(a) ) THEN
         IER = 105
         RETURN
  END IF
   W( 89 ) = JVS( 470 )
   W( 91 ) = JVS( 471 )
   W( 105 ) = JVS( 472 )
   W( 131 ) = JVS( 473 )
   W( 136 ) = JVS( 474 )
  a = -W( 89 ) / JVS( 396 )
  W( 89 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 397 )
  a = -W( 91 ) / JVS( 401 )
  W( 91 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 402 )
  JVS( 470) = W( 89 )
  JVS( 471) = W( 91 )
  JVS( 472) = W( 105 )
  JVS( 473) = W( 131 )
  JVS( 474) = W( 136 )
  IF ( ABS( JVS( 476 )) < TINY(a) ) THEN
         IER = 106
         RETURN
  END IF
   W( 36 ) = JVS( 475 )
   W( 106 ) = JVS( 476 )
   W( 127 ) = JVS( 477 )
   W( 129 ) = JVS( 478 )
   W( 133 ) = JVS( 479 )
   W( 136 ) = JVS( 480 )
  a = -W( 36 ) / JVS( 198 )
  W( 36 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 199 )
  JVS( 475) = W( 36 )
  JVS( 476) = W( 106 )
  JVS( 477) = W( 127 )
  JVS( 478) = W( 129 )
  JVS( 479) = W( 133 )
  JVS( 480) = W( 136 )
  IF ( ABS( JVS( 483 )) < TINY(a) ) THEN
         IER = 107
         RETURN
  END IF
   W( 90 ) = JVS( 481 )
   W( 105 ) = JVS( 482 )
   W( 107 ) = JVS( 483 )
   W( 116 ) = JVS( 484 )
   W( 121 ) = JVS( 485 )
   W( 123 ) = JVS( 486 )
   W( 125 ) = JVS( 487 )
   W( 129 ) = JVS( 488 )
   W( 131 ) = JVS( 489 )
   W( 132 ) = JVS( 490 )
   W( 136 ) = JVS( 491 )
  a = -W( 90 ) / JVS( 398 )
  W( 90 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 399 )
  W( 136 ) = W( 136 ) + a*JVS( 400 )
  a = -W( 105 ) / JVS( 472 )
  W( 105 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 473 )
  W( 136 ) = W( 136 ) + a*JVS( 474 )
  JVS( 481) = W( 90 )
  JVS( 482) = W( 105 )
  JVS( 483) = W( 107 )
  JVS( 484) = W( 116 )
  JVS( 485) = W( 121 )
  JVS( 486) = W( 123 )
  JVS( 487) = W( 125 )
  JVS( 488) = W( 129 )
  JVS( 489) = W( 131 )
  JVS( 490) = W( 132 )
  JVS( 491) = W( 136 )
  IF ( ABS( JVS( 493 )) < TINY(a) ) THEN
         IER = 108
         RETURN
  END IF
   W( 103 ) = JVS( 492 )
   W( 108 ) = JVS( 493 )
   W( 110 ) = JVS( 494 )
   W( 112 ) = JVS( 495 )
   W( 114 ) = JVS( 496 )
   W( 115 ) = JVS( 497 )
   W( 116 ) = JVS( 498 )
   W( 119 ) = JVS( 499 )
   W( 121 ) = JVS( 500 )
   W( 123 ) = JVS( 501 )
   W( 125 ) = JVS( 502 )
   W( 127 ) = JVS( 503 )
   W( 131 ) = JVS( 504 )
   W( 133 ) = JVS( 505 )
   W( 136 ) = JVS( 506 )
  a = -W( 103 ) / JVS( 461 )
  W( 103 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 462 )
  W( 131 ) = W( 131 ) + a*JVS( 463 )
  JVS( 492) = W( 103 )
  JVS( 493) = W( 108 )
  JVS( 494) = W( 110 )
  JVS( 495) = W( 112 )
  JVS( 496) = W( 114 )
  JVS( 497) = W( 115 )
  JVS( 498) = W( 116 )
  JVS( 499) = W( 119 )
  JVS( 500) = W( 121 )
  JVS( 501) = W( 123 )
  JVS( 502) = W( 125 )
  JVS( 503) = W( 127 )
  JVS( 504) = W( 131 )
  JVS( 505) = W( 133 )
  JVS( 506) = W( 136 )
  IF ( ABS( JVS( 509 )) < TINY(a) ) THEN
         IER = 109
         RETURN
  END IF
   W( 91 ) = JVS( 507 )
   W( 99 ) = JVS( 508 )
   W( 109 ) = JVS( 509 )
   W( 111 ) = JVS( 510 )
   W( 113 ) = JVS( 511 )
   W( 114 ) = JVS( 512 )
   W( 119 ) = JVS( 513 )
   W( 120 ) = JVS( 514 )
   W( 125 ) = JVS( 515 )
   W( 126 ) = JVS( 516 )
   W( 127 ) = JVS( 517 )
   W( 128 ) = JVS( 518 )
   W( 131 ) = JVS( 519 )
   W( 132 ) = JVS( 520 )
   W( 133 ) = JVS( 521 )
   W( 135 ) = JVS( 522 )
   W( 136 ) = JVS( 523 )
  a = -W( 91 ) / JVS( 401 )
  W( 91 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 402 )
  a = -W( 99 ) / JVS( 436 )
  W( 99 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 437 )
  W( 114 ) = W( 114 ) + a*JVS( 438 )
  W( 119 ) = W( 119 ) + a*JVS( 439 )
  W( 120 ) = W( 120 ) + a*JVS( 440 )
  W( 126 ) = W( 126 ) + a*JVS( 441 )
  W( 127 ) = W( 127 ) + a*JVS( 442 )
  W( 128 ) = W( 128 ) + a*JVS( 443 )
  W( 131 ) = W( 131 ) + a*JVS( 444 )
  W( 133 ) = W( 133 ) + a*JVS( 445 )
  W( 135 ) = W( 135 ) + a*JVS( 446 )
  W( 136 ) = W( 136 ) + a*JVS( 447 )
  JVS( 507) = W( 91 )
  JVS( 508) = W( 99 )
  JVS( 509) = W( 109 )
  JVS( 510) = W( 111 )
  JVS( 511) = W( 113 )
  JVS( 512) = W( 114 )
  JVS( 513) = W( 119 )
  JVS( 514) = W( 120 )
  JVS( 515) = W( 125 )
  JVS( 516) = W( 126 )
  JVS( 517) = W( 127 )
  JVS( 518) = W( 128 )
  JVS( 519) = W( 131 )
  JVS( 520) = W( 132 )
  JVS( 521) = W( 133 )
  JVS( 522) = W( 135 )
  JVS( 523) = W( 136 )
  IF ( ABS( JVS( 526 )) < TINY(a) ) THEN
         IER = 110
         RETURN
  END IF
   W( 97 ) = JVS( 524 )
   W( 105 ) = JVS( 525 )
   W( 110 ) = JVS( 526 )
   W( 127 ) = JVS( 527 )
   W( 131 ) = JVS( 528 )
   W( 133 ) = JVS( 529 )
   W( 136 ) = JVS( 530 )
  a = -W( 97 ) / JVS( 427 )
  W( 97 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 428 )
  W( 133 ) = W( 133 ) + a*JVS( 429 )
  a = -W( 105 ) / JVS( 472 )
  W( 105 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 473 )
  W( 136 ) = W( 136 ) + a*JVS( 474 )
  JVS( 524) = W( 97 )
  JVS( 525) = W( 105 )
  JVS( 526) = W( 110 )
  JVS( 527) = W( 127 )
  JVS( 528) = W( 131 )
  JVS( 529) = W( 133 )
  JVS( 530) = W( 136 )
  IF ( ABS( JVS( 531 )) < TINY(a) ) THEN
         IER = 111
         RETURN
  END IF
   W( 111 ) = JVS( 531 )
   W( 115 ) = JVS( 532 )
   W( 132 ) = JVS( 533 )
   W( 133 ) = JVS( 534 )
   W( 136 ) = JVS( 535 )
  JVS( 531) = W( 111 )
  JVS( 532) = W( 115 )
  JVS( 533) = W( 132 )
  JVS( 534) = W( 133 )
  JVS( 535) = W( 136 )
  IF ( ABS( JVS( 536 )) < TINY(a) ) THEN
         IER = 112
         RETURN
  END IF
   W( 112 ) = JVS( 536 )
   W( 125 ) = JVS( 537 )
   W( 131 ) = JVS( 538 )
   W( 132 ) = JVS( 539 )
   W( 133 ) = JVS( 540 )
  JVS( 536) = W( 112 )
  JVS( 537) = W( 125 )
  JVS( 538) = W( 131 )
  JVS( 539) = W( 132 )
  JVS( 540) = W( 133 )
  IF ( ABS( JVS( 541 )) < TINY(a) ) THEN
         IER = 113
         RETURN
  END IF
   W( 113 ) = JVS( 541 )
   W( 115 ) = JVS( 542 )
   W( 131 ) = JVS( 543 )
   W( 132 ) = JVS( 544 )
   W( 133 ) = JVS( 545 )
  JVS( 541) = W( 113 )
  JVS( 542) = W( 115 )
  JVS( 543) = W( 131 )
  JVS( 544) = W( 132 )
  JVS( 545) = W( 133 )
  IF ( ABS( JVS( 546 )) < TINY(a) ) THEN
         IER = 114
         RETURN
  END IF
   W( 114 ) = JVS( 546 )
   W( 127 ) = JVS( 547 )
   W( 131 ) = JVS( 548 )
   W( 136 ) = JVS( 549 )
  JVS( 546) = W( 114 )
  JVS( 547) = W( 127 )
  JVS( 548) = W( 131 )
  JVS( 549) = W( 136 )
  IF ( ABS( JVS( 550 )) < TINY(a) ) THEN
         IER = 115
         RETURN
  END IF
   W( 115 ) = JVS( 550 )
   W( 127 ) = JVS( 551 )
   W( 131 ) = JVS( 552 )
   W( 136 ) = JVS( 553 )
  JVS( 550) = W( 115 )
  JVS( 551) = W( 127 )
  JVS( 552) = W( 131 )
  JVS( 553) = W( 136 )
  IF ( ABS( JVS( 562 )) < TINY(a) ) THEN
         IER = 116
         RETURN
  END IF
   W( 100 ) = JVS( 554 )
   W( 103 ) = JVS( 555 )
   W( 104 ) = JVS( 556 )
   W( 110 ) = JVS( 557 )
   W( 112 ) = JVS( 558 )
   W( 113 ) = JVS( 559 )
   W( 114 ) = JVS( 560 )
   W( 115 ) = JVS( 561 )
   W( 116 ) = JVS( 562 )
   W( 119 ) = JVS( 563 )
   W( 120 ) = JVS( 564 )
   W( 124 ) = JVS( 565 )
   W( 125 ) = JVS( 566 )
   W( 127 ) = JVS( 567 )
   W( 131 ) = JVS( 568 )
   W( 132 ) = JVS( 569 )
   W( 133 ) = JVS( 570 )
   W( 134 ) = JVS( 571 )
   W( 136 ) = JVS( 572 )
  a = -W( 100 ) / JVS( 448 )
  W( 100 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 449 )
  W( 131 ) = W( 131 ) + a*JVS( 450 )
  W( 132 ) = W( 132 ) + a*JVS( 451 )
  a = -W( 103 ) / JVS( 461 )
  W( 103 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 462 )
  W( 131 ) = W( 131 ) + a*JVS( 463 )
  a = -W( 104 ) / JVS( 464 )
  W( 104 ) = -a
  W( 114 ) = W( 114 ) + a*JVS( 465 )
  W( 119 ) = W( 119 ) + a*JVS( 466 )
  W( 124 ) = W( 124 ) + a*JVS( 467 )
  W( 127 ) = W( 127 ) + a*JVS( 468 )
  W( 131 ) = W( 131 ) + a*JVS( 469 )
  a = -W( 110 ) / JVS( 526 )
  W( 110 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 527 )
  W( 131 ) = W( 131 ) + a*JVS( 528 )
  W( 133 ) = W( 133 ) + a*JVS( 529 )
  W( 136 ) = W( 136 ) + a*JVS( 530 )
  a = -W( 112 ) / JVS( 536 )
  W( 112 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 537 )
  W( 131 ) = W( 131 ) + a*JVS( 538 )
  W( 132 ) = W( 132 ) + a*JVS( 539 )
  W( 133 ) = W( 133 ) + a*JVS( 540 )
  a = -W( 113 ) / JVS( 541 )
  W( 113 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 542 )
  W( 131 ) = W( 131 ) + a*JVS( 543 )
  W( 132 ) = W( 132 ) + a*JVS( 544 )
  W( 133 ) = W( 133 ) + a*JVS( 545 )
  a = -W( 114 ) / JVS( 546 )
  W( 114 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 547 )
  W( 131 ) = W( 131 ) + a*JVS( 548 )
  W( 136 ) = W( 136 ) + a*JVS( 549 )
  a = -W( 115 ) / JVS( 550 )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 136 ) = W( 136 ) + a*JVS( 553 )
  JVS( 554) = W( 100 )
  JVS( 555) = W( 103 )
  JVS( 556) = W( 104 )
  JVS( 557) = W( 110 )
  JVS( 558) = W( 112 )
  JVS( 559) = W( 113 )
  JVS( 560) = W( 114 )
  JVS( 561) = W( 115 )
  JVS( 562) = W( 116 )
  JVS( 563) = W( 119 )
  JVS( 564) = W( 120 )
  JVS( 565) = W( 124 )
  JVS( 566) = W( 125 )
  JVS( 567) = W( 127 )
  JVS( 568) = W( 131 )
  JVS( 569) = W( 132 )
  JVS( 570) = W( 133 )
  JVS( 571) = W( 134 )
  JVS( 572) = W( 136 )
  IF ( ABS( JVS( 580 )) < TINY(a) ) THEN
         IER = 117
         RETURN
  END IF
   W( 89 ) = JVS( 573 )
   W( 91 ) = JVS( 574 )
   W( 103 ) = JVS( 575 )
   W( 105 ) = JVS( 576 )
   W( 110 ) = JVS( 577 )
   W( 114 ) = JVS( 578 )
   W( 115 ) = JVS( 579 )
   W( 117 ) = JVS( 580 )
   W( 119 ) = JVS( 581 )
   W( 121 ) = JVS( 582 )
   W( 125 ) = JVS( 583 )
   W( 126 ) = JVS( 584 )
   W( 127 ) = JVS( 585 )
   W( 128 ) = JVS( 586 )
   W( 131 ) = JVS( 587 )
   W( 132 ) = JVS( 588 )
   W( 133 ) = JVS( 589 )
   W( 135 ) = JVS( 590 )
   W( 136 ) = JVS( 591 )
  a = -W( 89 ) / JVS( 396 )
  W( 89 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 397 )
  a = -W( 91 ) / JVS( 401 )
  W( 91 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 402 )
  a = -W( 103 ) / JVS( 461 )
  W( 103 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 462 )
  W( 131 ) = W( 131 ) + a*JVS( 463 )
  a = -W( 105 ) / JVS( 472 )
  W( 105 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 473 )
  W( 136 ) = W( 136 ) + a*JVS( 474 )
  a = -W( 110 ) / JVS( 526 )
  W( 110 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 527 )
  W( 131 ) = W( 131 ) + a*JVS( 528 )
  W( 133 ) = W( 133 ) + a*JVS( 529 )
  W( 136 ) = W( 136 ) + a*JVS( 530 )
  a = -W( 114 ) / JVS( 546 )
  W( 114 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 547 )
  W( 131 ) = W( 131 ) + a*JVS( 548 )
  W( 136 ) = W( 136 ) + a*JVS( 549 )
  a = -W( 115 ) / JVS( 550 )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 136 ) = W( 136 ) + a*JVS( 553 )
  JVS( 573) = W( 89 )
  JVS( 574) = W( 91 )
  JVS( 575) = W( 103 )
  JVS( 576) = W( 105 )
  JVS( 577) = W( 110 )
  JVS( 578) = W( 114 )
  JVS( 579) = W( 115 )
  JVS( 580) = W( 117 )
  JVS( 581) = W( 119 )
  JVS( 582) = W( 121 )
  JVS( 583) = W( 125 )
  JVS( 584) = W( 126 )
  JVS( 585) = W( 127 )
  JVS( 586) = W( 128 )
  JVS( 587) = W( 131 )
  JVS( 588) = W( 132 )
  JVS( 589) = W( 133 )
  JVS( 590) = W( 135 )
  JVS( 591) = W( 136 )
  IF ( ABS( JVS( 593 )) < TINY(a) ) THEN
         IER = 118
         RETURN
  END IF
   W( 112 ) = JVS( 592 )
   W( 118 ) = JVS( 593 )
   W( 119 ) = JVS( 594 )
   W( 125 ) = JVS( 595 )
   W( 126 ) = JVS( 596 )
   W( 127 ) = JVS( 597 )
   W( 128 ) = JVS( 598 )
   W( 131 ) = JVS( 599 )
   W( 132 ) = JVS( 600 )
   W( 133 ) = JVS( 601 )
   W( 134 ) = JVS( 602 )
   W( 135 ) = JVS( 603 )
   W( 136 ) = JVS( 604 )
  a = -W( 112 ) / JVS( 536 )
  W( 112 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 537 )
  W( 131 ) = W( 131 ) + a*JVS( 538 )
  W( 132 ) = W( 132 ) + a*JVS( 539 )
  W( 133 ) = W( 133 ) + a*JVS( 540 )
  JVS( 592) = W( 112 )
  JVS( 593) = W( 118 )
  JVS( 594) = W( 119 )
  JVS( 595) = W( 125 )
  JVS( 596) = W( 126 )
  JVS( 597) = W( 127 )
  JVS( 598) = W( 128 )
  JVS( 599) = W( 131 )
  JVS( 600) = W( 132 )
  JVS( 601) = W( 133 )
  JVS( 602) = W( 134 )
  JVS( 603) = W( 135 )
  JVS( 604) = W( 136 )
  IF ( ABS( JVS( 605 )) < TINY(a) ) THEN
         IER = 119
         RETURN
  END IF
   W( 119 ) = JVS( 605 )
   W( 127 ) = JVS( 606 )
   W( 131 ) = JVS( 607 )
   W( 136 ) = JVS( 608 )
  JVS( 605) = W( 119 )
  JVS( 606) = W( 127 )
  JVS( 607) = W( 131 )
  JVS( 608) = W( 136 )
  IF ( ABS( JVS( 611 )) < TINY(a) ) THEN
         IER = 120
         RETURN
  END IF
   W( 114 ) = JVS( 609 )
   W( 119 ) = JVS( 610 )
   W( 120 ) = JVS( 611 )
   W( 126 ) = JVS( 612 )
   W( 127 ) = JVS( 613 )
   W( 131 ) = JVS( 614 )
   W( 132 ) = JVS( 615 )
   W( 133 ) = JVS( 616 )
   W( 136 ) = JVS( 617 )
  a = -W( 114 ) / JVS( 546 )
  W( 114 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 547 )
  W( 131 ) = W( 131 ) + a*JVS( 548 )
  W( 136 ) = W( 136 ) + a*JVS( 549 )
  a = -W( 119 ) / JVS( 605 )
  W( 119 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 606 )
  W( 131 ) = W( 131 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  JVS( 609) = W( 114 )
  JVS( 610) = W( 119 )
  JVS( 611) = W( 120 )
  JVS( 612) = W( 126 )
  JVS( 613) = W( 127 )
  JVS( 614) = W( 131 )
  JVS( 615) = W( 132 )
  JVS( 616) = W( 133 )
  JVS( 617) = W( 136 )
  IF ( ABS( JVS( 623 )) < TINY(a) ) THEN
         IER = 121
         RETURN
  END IF
   W( 91 ) = JVS( 618 )
   W( 110 ) = JVS( 619 )
   W( 112 ) = JVS( 620 )
   W( 114 ) = JVS( 621 )
   W( 119 ) = JVS( 622 )
   W( 121 ) = JVS( 623 )
   W( 125 ) = JVS( 624 )
   W( 127 ) = JVS( 625 )
   W( 131 ) = JVS( 626 )
   W( 132 ) = JVS( 627 )
   W( 133 ) = JVS( 628 )
   W( 134 ) = JVS( 629 )
   W( 135 ) = JVS( 630 )
   W( 136 ) = JVS( 631 )
  a = -W( 91 ) / JVS( 401 )
  W( 91 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 402 )
  a = -W( 110 ) / JVS( 526 )
  W( 110 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 527 )
  W( 131 ) = W( 131 ) + a*JVS( 528 )
  W( 133 ) = W( 133 ) + a*JVS( 529 )
  W( 136 ) = W( 136 ) + a*JVS( 530 )
  a = -W( 112 ) / JVS( 536 )
  W( 112 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 537 )
  W( 131 ) = W( 131 ) + a*JVS( 538 )
  W( 132 ) = W( 132 ) + a*JVS( 539 )
  W( 133 ) = W( 133 ) + a*JVS( 540 )
  a = -W( 114 ) / JVS( 546 )
  W( 114 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 547 )
  W( 131 ) = W( 131 ) + a*JVS( 548 )
  W( 136 ) = W( 136 ) + a*JVS( 549 )
  a = -W( 119 ) / JVS( 605 )
  W( 119 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 606 )
  W( 131 ) = W( 131 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  JVS( 618) = W( 91 )
  JVS( 619) = W( 110 )
  JVS( 620) = W( 112 )
  JVS( 621) = W( 114 )
  JVS( 622) = W( 119 )
  JVS( 623) = W( 121 )
  JVS( 624) = W( 125 )
  JVS( 625) = W( 127 )
  JVS( 626) = W( 131 )
  JVS( 627) = W( 132 )
  JVS( 628) = W( 133 )
  JVS( 629) = W( 134 )
  JVS( 630) = W( 135 )
  JVS( 631) = W( 136 )
  IF ( ABS( JVS( 636 )) < TINY(a) ) THEN
         IER = 122
         RETURN
  END IF
   W( 98 ) = JVS( 632 )
   W( 101 ) = JVS( 633 )
   W( 114 ) = JVS( 634 )
   W( 119 ) = JVS( 635 )
   W( 122 ) = JVS( 636 )
   W( 126 ) = JVS( 637 )
   W( 127 ) = JVS( 638 )
   W( 128 ) = JVS( 639 )
   W( 131 ) = JVS( 640 )
   W( 132 ) = JVS( 641 )
   W( 133 ) = JVS( 642 )
   W( 135 ) = JVS( 643 )
   W( 136 ) = JVS( 644 )
  a = -W( 98 ) / JVS( 430 )
  W( 98 ) = -a
  W( 114 ) = W( 114 ) + a*JVS( 431 )
  W( 119 ) = W( 119 ) + a*JVS( 432 )
  W( 122 ) = W( 122 ) + a*JVS( 433 )
  W( 127 ) = W( 127 ) + a*JVS( 434 )
  W( 131 ) = W( 131 ) + a*JVS( 435 )
  a = -W( 101 ) / JVS( 452 )
  W( 101 ) = -a
  W( 122 ) = W( 122 ) + a*JVS( 453 )
  W( 131 ) = W( 131 ) + a*JVS( 454 )
  W( 132 ) = W( 132 ) + a*JVS( 455 )
  a = -W( 114 ) / JVS( 546 )
  W( 114 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 547 )
  W( 131 ) = W( 131 ) + a*JVS( 548 )
  W( 136 ) = W( 136 ) + a*JVS( 549 )
  a = -W( 119 ) / JVS( 605 )
  W( 119 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 606 )
  W( 131 ) = W( 131 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  JVS( 632) = W( 98 )
  JVS( 633) = W( 101 )
  JVS( 634) = W( 114 )
  JVS( 635) = W( 119 )
  JVS( 636) = W( 122 )
  JVS( 637) = W( 126 )
  JVS( 638) = W( 127 )
  JVS( 639) = W( 128 )
  JVS( 640) = W( 131 )
  JVS( 641) = W( 132 )
  JVS( 642) = W( 133 )
  JVS( 643) = W( 135 )
  JVS( 644) = W( 136 )
  IF ( ABS( JVS( 656 )) < TINY(a) ) THEN
         IER = 123
         RETURN
  END IF
   W( 48 ) = JVS( 645 )
   W( 101 ) = JVS( 646 )
   W( 103 ) = JVS( 647 )
   W( 110 ) = JVS( 648 )
   W( 111 ) = JVS( 649 )
   W( 112 ) = JVS( 650 )
   W( 114 ) = JVS( 651 )
   W( 115 ) = JVS( 652 )
   W( 119 ) = JVS( 653 )
   W( 120 ) = JVS( 654 )
   W( 122 ) = JVS( 655 )
   W( 123 ) = JVS( 656 )
   W( 125 ) = JVS( 657 )
   W( 126 ) = JVS( 658 )
   W( 127 ) = JVS( 659 )
   W( 128 ) = JVS( 660 )
   W( 131 ) = JVS( 661 )
   W( 132 ) = JVS( 662 )
   W( 133 ) = JVS( 663 )
   W( 135 ) = JVS( 664 )
   W( 136 ) = JVS( 665 )
  a = -W( 48 ) / JVS( 226 )
  W( 48 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 227 )
  a = -W( 101 ) / JVS( 452 )
  W( 101 ) = -a
  W( 122 ) = W( 122 ) + a*JVS( 453 )
  W( 131 ) = W( 131 ) + a*JVS( 454 )
  W( 132 ) = W( 132 ) + a*JVS( 455 )
  a = -W( 103 ) / JVS( 461 )
  W( 103 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 462 )
  W( 131 ) = W( 131 ) + a*JVS( 463 )
  a = -W( 110 ) / JVS( 526 )
  W( 110 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 527 )
  W( 131 ) = W( 131 ) + a*JVS( 528 )
  W( 133 ) = W( 133 ) + a*JVS( 529 )
  W( 136 ) = W( 136 ) + a*JVS( 530 )
  a = -W( 111 ) / JVS( 531 )
  W( 111 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 532 )
  W( 132 ) = W( 132 ) + a*JVS( 533 )
  W( 133 ) = W( 133 ) + a*JVS( 534 )
  W( 136 ) = W( 136 ) + a*JVS( 535 )
  a = -W( 112 ) / JVS( 536 )
  W( 112 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 537 )
  W( 131 ) = W( 131 ) + a*JVS( 538 )
  W( 132 ) = W( 132 ) + a*JVS( 539 )
  W( 133 ) = W( 133 ) + a*JVS( 540 )
  a = -W( 114 ) / JVS( 546 )
  W( 114 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 547 )
  W( 131 ) = W( 131 ) + a*JVS( 548 )
  W( 136 ) = W( 136 ) + a*JVS( 549 )
  a = -W( 115 ) / JVS( 550 )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 136 ) = W( 136 ) + a*JVS( 553 )
  a = -W( 119 ) / JVS( 605 )
  W( 119 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 606 )
  W( 131 ) = W( 131 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 120 ) / JVS( 611 )
  W( 120 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 612 )
  W( 127 ) = W( 127 ) + a*JVS( 613 )
  W( 131 ) = W( 131 ) + a*JVS( 614 )
  W( 132 ) = W( 132 ) + a*JVS( 615 )
  W( 133 ) = W( 133 ) + a*JVS( 616 )
  W( 136 ) = W( 136 ) + a*JVS( 617 )
  a = -W( 122 ) / JVS( 636 )
  W( 122 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 637 )
  W( 127 ) = W( 127 ) + a*JVS( 638 )
  W( 128 ) = W( 128 ) + a*JVS( 639 )
  W( 131 ) = W( 131 ) + a*JVS( 640 )
  W( 132 ) = W( 132 ) + a*JVS( 641 )
  W( 133 ) = W( 133 ) + a*JVS( 642 )
  W( 135 ) = W( 135 ) + a*JVS( 643 )
  W( 136 ) = W( 136 ) + a*JVS( 644 )
  JVS( 645) = W( 48 )
  JVS( 646) = W( 101 )
  JVS( 647) = W( 103 )
  JVS( 648) = W( 110 )
  JVS( 649) = W( 111 )
  JVS( 650) = W( 112 )
  JVS( 651) = W( 114 )
  JVS( 652) = W( 115 )
  JVS( 653) = W( 119 )
  JVS( 654) = W( 120 )
  JVS( 655) = W( 122 )
  JVS( 656) = W( 123 )
  JVS( 657) = W( 125 )
  JVS( 658) = W( 126 )
  JVS( 659) = W( 127 )
  JVS( 660) = W( 128 )
  JVS( 661) = W( 131 )
  JVS( 662) = W( 132 )
  JVS( 663) = W( 133 )
  JVS( 664) = W( 135 )
  JVS( 665) = W( 136 )
  IF ( ABS( JVS( 672 )) < TINY(a) ) THEN
         IER = 124
         RETURN
  END IF
   W( 93 ) = JVS( 666 )
   W( 100 ) = JVS( 667 )
   W( 114 ) = JVS( 668 )
   W( 118 ) = JVS( 669 )
   W( 119 ) = JVS( 670 )
   W( 123 ) = JVS( 671 )
   W( 124 ) = JVS( 672 )
   W( 125 ) = JVS( 673 )
   W( 126 ) = JVS( 674 )
   W( 127 ) = JVS( 675 )
   W( 128 ) = JVS( 676 )
   W( 130 ) = JVS( 677 )
   W( 131 ) = JVS( 678 )
   W( 132 ) = JVS( 679 )
   W( 133 ) = JVS( 680 )
   W( 134 ) = JVS( 681 )
   W( 135 ) = JVS( 682 )
   W( 136 ) = JVS( 683 )
  a = -W( 93 ) / JVS( 408 )
  W( 93 ) = -a
  W( 114 ) = W( 114 ) + a*JVS( 409 )
  W( 119 ) = W( 119 ) + a*JVS( 410 )
  W( 127 ) = W( 127 ) + a*JVS( 411 )
  W( 131 ) = W( 131 ) + a*JVS( 412 )
  a = -W( 100 ) / JVS( 448 )
  W( 100 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 449 )
  W( 131 ) = W( 131 ) + a*JVS( 450 )
  W( 132 ) = W( 132 ) + a*JVS( 451 )
  a = -W( 114 ) / JVS( 546 )
  W( 114 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 547 )
  W( 131 ) = W( 131 ) + a*JVS( 548 )
  W( 136 ) = W( 136 ) + a*JVS( 549 )
  a = -W( 118 ) / JVS( 593 )
  W( 118 ) = -a
  W( 119 ) = W( 119 ) + a*JVS( 594 )
  W( 125 ) = W( 125 ) + a*JVS( 595 )
  W( 126 ) = W( 126 ) + a*JVS( 596 )
  W( 127 ) = W( 127 ) + a*JVS( 597 )
  W( 128 ) = W( 128 ) + a*JVS( 598 )
  W( 131 ) = W( 131 ) + a*JVS( 599 )
  W( 132 ) = W( 132 ) + a*JVS( 600 )
  W( 133 ) = W( 133 ) + a*JVS( 601 )
  W( 134 ) = W( 134 ) + a*JVS( 602 )
  W( 135 ) = W( 135 ) + a*JVS( 603 )
  W( 136 ) = W( 136 ) + a*JVS( 604 )
  a = -W( 119 ) / JVS( 605 )
  W( 119 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 606 )
  W( 131 ) = W( 131 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 123 ) / JVS( 656 )
  W( 123 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 657 )
  W( 126 ) = W( 126 ) + a*JVS( 658 )
  W( 127 ) = W( 127 ) + a*JVS( 659 )
  W( 128 ) = W( 128 ) + a*JVS( 660 )
  W( 131 ) = W( 131 ) + a*JVS( 661 )
  W( 132 ) = W( 132 ) + a*JVS( 662 )
  W( 133 ) = W( 133 ) + a*JVS( 663 )
  W( 135 ) = W( 135 ) + a*JVS( 664 )
  W( 136 ) = W( 136 ) + a*JVS( 665 )
  JVS( 666) = W( 93 )
  JVS( 667) = W( 100 )
  JVS( 668) = W( 114 )
  JVS( 669) = W( 118 )
  JVS( 670) = W( 119 )
  JVS( 671) = W( 123 )
  JVS( 672) = W( 124 )
  JVS( 673) = W( 125 )
  JVS( 674) = W( 126 )
  JVS( 675) = W( 127 )
  JVS( 676) = W( 128 )
  JVS( 677) = W( 130 )
  JVS( 678) = W( 131 )
  JVS( 679) = W( 132 )
  JVS( 680) = W( 133 )
  JVS( 681) = W( 134 )
  JVS( 682) = W( 135 )
  JVS( 683) = W( 136 )
  IF ( ABS( JVS( 687 )) < TINY(a) ) THEN
         IER = 125
         RETURN
  END IF
   W( 111 ) = JVS( 684 )
   W( 113 ) = JVS( 685 )
   W( 115 ) = JVS( 686 )
   W( 125 ) = JVS( 687 )
   W( 127 ) = JVS( 688 )
   W( 131 ) = JVS( 689 )
   W( 132 ) = JVS( 690 )
   W( 133 ) = JVS( 691 )
   W( 136 ) = JVS( 692 )
  a = -W( 111 ) / JVS( 531 )
  W( 111 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 532 )
  W( 132 ) = W( 132 ) + a*JVS( 533 )
  W( 133 ) = W( 133 ) + a*JVS( 534 )
  W( 136 ) = W( 136 ) + a*JVS( 535 )
  a = -W( 113 ) / JVS( 541 )
  W( 113 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 542 )
  W( 131 ) = W( 131 ) + a*JVS( 543 )
  W( 132 ) = W( 132 ) + a*JVS( 544 )
  W( 133 ) = W( 133 ) + a*JVS( 545 )
  a = -W( 115 ) / JVS( 550 )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 136 ) = W( 136 ) + a*JVS( 553 )
  JVS( 684) = W( 111 )
  JVS( 685) = W( 113 )
  JVS( 686) = W( 115 )
  JVS( 687) = W( 125 )
  JVS( 688) = W( 127 )
  JVS( 689) = W( 131 )
  JVS( 690) = W( 132 )
  JVS( 691) = W( 133 )
  JVS( 692) = W( 136 )
  IF ( ABS( JVS( 701 )) < TINY(a) ) THEN
         IER = 126
         RETURN
  END IF
   W( 92 ) = JVS( 693 )
   W( 97 ) = JVS( 694 )
   W( 105 ) = JVS( 695 )
   W( 111 ) = JVS( 696 )
   W( 113 ) = JVS( 697 )
   W( 115 ) = JVS( 698 )
   W( 120 ) = JVS( 699 )
   W( 125 ) = JVS( 700 )
   W( 126 ) = JVS( 701 )
   W( 127 ) = JVS( 702 )
   W( 128 ) = JVS( 703 )
   W( 129 ) = JVS( 704 )
   W( 131 ) = JVS( 705 )
   W( 132 ) = JVS( 706 )
   W( 133 ) = JVS( 707 )
   W( 136 ) = JVS( 708 )
  a = -W( 92 ) / JVS( 403 )
  W( 92 ) = -a
  W( 105 ) = W( 105 ) + a*JVS( 404 )
  W( 129 ) = W( 129 ) + a*JVS( 405 )
  W( 131 ) = W( 131 ) + a*JVS( 406 )
  W( 136 ) = W( 136 ) + a*JVS( 407 )
  a = -W( 97 ) / JVS( 427 )
  W( 97 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 428 )
  W( 133 ) = W( 133 ) + a*JVS( 429 )
  a = -W( 105 ) / JVS( 472 )
  W( 105 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 473 )
  W( 136 ) = W( 136 ) + a*JVS( 474 )
  a = -W( 111 ) / JVS( 531 )
  W( 111 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 532 )
  W( 132 ) = W( 132 ) + a*JVS( 533 )
  W( 133 ) = W( 133 ) + a*JVS( 534 )
  W( 136 ) = W( 136 ) + a*JVS( 535 )
  a = -W( 113 ) / JVS( 541 )
  W( 113 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 542 )
  W( 131 ) = W( 131 ) + a*JVS( 543 )
  W( 132 ) = W( 132 ) + a*JVS( 544 )
  W( 133 ) = W( 133 ) + a*JVS( 545 )
  a = -W( 115 ) / JVS( 550 )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 136 ) = W( 136 ) + a*JVS( 553 )
  a = -W( 120 ) / JVS( 611 )
  W( 120 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 612 )
  W( 127 ) = W( 127 ) + a*JVS( 613 )
  W( 131 ) = W( 131 ) + a*JVS( 614 )
  W( 132 ) = W( 132 ) + a*JVS( 615 )
  W( 133 ) = W( 133 ) + a*JVS( 616 )
  W( 136 ) = W( 136 ) + a*JVS( 617 )
  a = -W( 125 ) / JVS( 687 )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 688 )
  W( 131 ) = W( 131 ) + a*JVS( 689 )
  W( 132 ) = W( 132 ) + a*JVS( 690 )
  W( 133 ) = W( 133 ) + a*JVS( 691 )
  W( 136 ) = W( 136 ) + a*JVS( 692 )
  JVS( 693) = W( 92 )
  JVS( 694) = W( 97 )
  JVS( 695) = W( 105 )
  JVS( 696) = W( 111 )
  JVS( 697) = W( 113 )
  JVS( 698) = W( 115 )
  JVS( 699) = W( 120 )
  JVS( 700) = W( 125 )
  JVS( 701) = W( 126 )
  JVS( 702) = W( 127 )
  JVS( 703) = W( 128 )
  JVS( 704) = W( 129 )
  JVS( 705) = W( 131 )
  JVS( 706) = W( 132 )
  JVS( 707) = W( 133 )
  JVS( 708) = W( 136 )
  IF ( ABS( JVS( 718 )) < TINY(a) ) THEN
         IER = 127
         RETURN
  END IF
   W( 94 ) = JVS( 709 )
   W( 95 ) = JVS( 710 )
   W( 103 ) = JVS( 711 )
   W( 106 ) = JVS( 712 )
   W( 110 ) = JVS( 713 )
   W( 114 ) = JVS( 714 )
   W( 115 ) = JVS( 715 )
   W( 119 ) = JVS( 716 )
   W( 125 ) = JVS( 717 )
   W( 127 ) = JVS( 718 )
   W( 129 ) = JVS( 719 )
   W( 130 ) = JVS( 720 )
   W( 131 ) = JVS( 721 )
   W( 132 ) = JVS( 722 )
   W( 133 ) = JVS( 723 )
   W( 136 ) = JVS( 724 )
  a = -W( 94 ) / JVS( 413 )
  W( 94 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 414 )
  W( 131 ) = W( 131 ) + a*JVS( 415 )
  W( 136 ) = W( 136 ) + a*JVS( 416 )
  a = -W( 95 ) / JVS( 417 )
  W( 95 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 418 )
  W( 131 ) = W( 131 ) + a*JVS( 419 )
  W( 136 ) = W( 136 ) + a*JVS( 420 )
  a = -W( 103 ) / JVS( 461 )
  W( 103 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 462 )
  W( 131 ) = W( 131 ) + a*JVS( 463 )
  a = -W( 106 ) / JVS( 476 )
  W( 106 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 477 )
  W( 129 ) = W( 129 ) + a*JVS( 478 )
  W( 133 ) = W( 133 ) + a*JVS( 479 )
  W( 136 ) = W( 136 ) + a*JVS( 480 )
  a = -W( 110 ) / JVS( 526 )
  W( 110 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 527 )
  W( 131 ) = W( 131 ) + a*JVS( 528 )
  W( 133 ) = W( 133 ) + a*JVS( 529 )
  W( 136 ) = W( 136 ) + a*JVS( 530 )
  a = -W( 114 ) / JVS( 546 )
  W( 114 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 547 )
  W( 131 ) = W( 131 ) + a*JVS( 548 )
  W( 136 ) = W( 136 ) + a*JVS( 549 )
  a = -W( 115 ) / JVS( 550 )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 136 ) = W( 136 ) + a*JVS( 553 )
  a = -W( 119 ) / JVS( 605 )
  W( 119 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 606 )
  W( 131 ) = W( 131 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 125 ) / JVS( 687 )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 688 )
  W( 131 ) = W( 131 ) + a*JVS( 689 )
  W( 132 ) = W( 132 ) + a*JVS( 690 )
  W( 133 ) = W( 133 ) + a*JVS( 691 )
  W( 136 ) = W( 136 ) + a*JVS( 692 )
  JVS( 709) = W( 94 )
  JVS( 710) = W( 95 )
  JVS( 711) = W( 103 )
  JVS( 712) = W( 106 )
  JVS( 713) = W( 110 )
  JVS( 714) = W( 114 )
  JVS( 715) = W( 115 )
  JVS( 716) = W( 119 )
  JVS( 717) = W( 125 )
  JVS( 718) = W( 127 )
  JVS( 719) = W( 129 )
  JVS( 720) = W( 130 )
  JVS( 721) = W( 131 )
  JVS( 722) = W( 132 )
  JVS( 723) = W( 133 )
  JVS( 724) = W( 136 )
  IF ( ABS( JVS( 735 )) < TINY(a) ) THEN
         IER = 128
         RETURN
  END IF
   W( 109 ) = JVS( 725 )
   W( 111 ) = JVS( 726 )
   W( 113 ) = JVS( 727 )
   W( 114 ) = JVS( 728 )
   W( 115 ) = JVS( 729 )
   W( 119 ) = JVS( 730 )
   W( 120 ) = JVS( 731 )
   W( 125 ) = JVS( 732 )
   W( 126 ) = JVS( 733 )
   W( 127 ) = JVS( 734 )
   W( 128 ) = JVS( 735 )
   W( 129 ) = JVS( 736 )
   W( 130 ) = JVS( 737 )
   W( 131 ) = JVS( 738 )
   W( 132 ) = JVS( 739 )
   W( 133 ) = JVS( 740 )
   W( 135 ) = JVS( 741 )
   W( 136 ) = JVS( 742 )
  a = -W( 109 ) / JVS( 509 )
  W( 109 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 510 )
  W( 113 ) = W( 113 ) + a*JVS( 511 )
  W( 114 ) = W( 114 ) + a*JVS( 512 )
  W( 119 ) = W( 119 ) + a*JVS( 513 )
  W( 120 ) = W( 120 ) + a*JVS( 514 )
  W( 125 ) = W( 125 ) + a*JVS( 515 )
  W( 126 ) = W( 126 ) + a*JVS( 516 )
  W( 127 ) = W( 127 ) + a*JVS( 517 )
  W( 128 ) = W( 128 ) + a*JVS( 518 )
  W( 131 ) = W( 131 ) + a*JVS( 519 )
  W( 132 ) = W( 132 ) + a*JVS( 520 )
  W( 133 ) = W( 133 ) + a*JVS( 521 )
  W( 135 ) = W( 135 ) + a*JVS( 522 )
  W( 136 ) = W( 136 ) + a*JVS( 523 )
  a = -W( 111 ) / JVS( 531 )
  W( 111 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 532 )
  W( 132 ) = W( 132 ) + a*JVS( 533 )
  W( 133 ) = W( 133 ) + a*JVS( 534 )
  W( 136 ) = W( 136 ) + a*JVS( 535 )
  a = -W( 113 ) / JVS( 541 )
  W( 113 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 542 )
  W( 131 ) = W( 131 ) + a*JVS( 543 )
  W( 132 ) = W( 132 ) + a*JVS( 544 )
  W( 133 ) = W( 133 ) + a*JVS( 545 )
  a = -W( 114 ) / JVS( 546 )
  W( 114 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 547 )
  W( 131 ) = W( 131 ) + a*JVS( 548 )
  W( 136 ) = W( 136 ) + a*JVS( 549 )
  a = -W( 115 ) / JVS( 550 )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 136 ) = W( 136 ) + a*JVS( 553 )
  a = -W( 119 ) / JVS( 605 )
  W( 119 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 606 )
  W( 131 ) = W( 131 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 120 ) / JVS( 611 )
  W( 120 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 612 )
  W( 127 ) = W( 127 ) + a*JVS( 613 )
  W( 131 ) = W( 131 ) + a*JVS( 614 )
  W( 132 ) = W( 132 ) + a*JVS( 615 )
  W( 133 ) = W( 133 ) + a*JVS( 616 )
  W( 136 ) = W( 136 ) + a*JVS( 617 )
  a = -W( 125 ) / JVS( 687 )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 688 )
  W( 131 ) = W( 131 ) + a*JVS( 689 )
  W( 132 ) = W( 132 ) + a*JVS( 690 )
  W( 133 ) = W( 133 ) + a*JVS( 691 )
  W( 136 ) = W( 136 ) + a*JVS( 692 )
  a = -W( 126 ) / JVS( 701 )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 702 )
  W( 128 ) = W( 128 ) + a*JVS( 703 )
  W( 129 ) = W( 129 ) + a*JVS( 704 )
  W( 131 ) = W( 131 ) + a*JVS( 705 )
  W( 132 ) = W( 132 ) + a*JVS( 706 )
  W( 133 ) = W( 133 ) + a*JVS( 707 )
  W( 136 ) = W( 136 ) + a*JVS( 708 )
  a = -W( 127 ) / JVS( 718 )
  W( 127 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 719 )
  W( 130 ) = W( 130 ) + a*JVS( 720 )
  W( 131 ) = W( 131 ) + a*JVS( 721 )
  W( 132 ) = W( 132 ) + a*JVS( 722 )
  W( 133 ) = W( 133 ) + a*JVS( 723 )
  W( 136 ) = W( 136 ) + a*JVS( 724 )
  JVS( 725) = W( 109 )
  JVS( 726) = W( 111 )
  JVS( 727) = W( 113 )
  JVS( 728) = W( 114 )
  JVS( 729) = W( 115 )
  JVS( 730) = W( 119 )
  JVS( 731) = W( 120 )
  JVS( 732) = W( 125 )
  JVS( 733) = W( 126 )
  JVS( 734) = W( 127 )
  JVS( 735) = W( 128 )
  JVS( 736) = W( 129 )
  JVS( 737) = W( 130 )
  JVS( 738) = W( 131 )
  JVS( 739) = W( 132 )
  JVS( 740) = W( 133 )
  JVS( 741) = W( 135 )
  JVS( 742) = W( 136 )
  IF ( ABS( JVS( 768 )) < TINY(a) ) THEN
         IER = 129
         RETURN
  END IF
   W( 87 ) = JVS( 743 )
   W( 90 ) = JVS( 744 )
   W( 92 ) = JVS( 745 )
   W( 96 ) = JVS( 746 )
   W( 97 ) = JVS( 747 )
   W( 102 ) = JVS( 748 )
   W( 105 ) = JVS( 749 )
   W( 106 ) = JVS( 750 )
   W( 107 ) = JVS( 751 )
   W( 111 ) = JVS( 752 )
   W( 112 ) = JVS( 753 )
   W( 113 ) = JVS( 754 )
   W( 115 ) = JVS( 755 )
   W( 116 ) = JVS( 756 )
   W( 117 ) = JVS( 757 )
   W( 119 ) = JVS( 758 )
   W( 120 ) = JVS( 759 )
   W( 121 ) = JVS( 760 )
   W( 122 ) = JVS( 761 )
   W( 123 ) = JVS( 762 )
   W( 124 ) = JVS( 763 )
   W( 125 ) = JVS( 764 )
   W( 126 ) = JVS( 765 )
   W( 127 ) = JVS( 766 )
   W( 128 ) = JVS( 767 )
   W( 129 ) = JVS( 768 )
   W( 130 ) = JVS( 769 )
   W( 131 ) = JVS( 770 )
   W( 132 ) = JVS( 771 )
   W( 133 ) = JVS( 772 )
   W( 134 ) = JVS( 773 )
   W( 135 ) = JVS( 774 )
   W( 136 ) = JVS( 775 )
  a = -W( 87 ) / JVS( 390 )
  W( 87 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 391 )
  W( 130 ) = W( 130 ) + a*JVS( 392 )
  a = -W( 90 ) / JVS( 398 )
  W( 90 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 399 )
  W( 136 ) = W( 136 ) + a*JVS( 400 )
  a = -W( 92 ) / JVS( 403 )
  W( 92 ) = -a
  W( 105 ) = W( 105 ) + a*JVS( 404 )
  W( 129 ) = W( 129 ) + a*JVS( 405 )
  W( 131 ) = W( 131 ) + a*JVS( 406 )
  W( 136 ) = W( 136 ) + a*JVS( 407 )
  a = -W( 96 ) / JVS( 421 )
  W( 96 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 422 )
  W( 131 ) = W( 131 ) + a*JVS( 423 )
  W( 132 ) = W( 132 ) + a*JVS( 424 )
  a = -W( 97 ) / JVS( 427 )
  W( 97 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 428 )
  W( 133 ) = W( 133 ) + a*JVS( 429 )
  a = -W( 102 ) / JVS( 456 )
  W( 102 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 457 )
  W( 131 ) = W( 131 ) + a*JVS( 458 )
  W( 132 ) = W( 132 ) + a*JVS( 459 )
  W( 133 ) = W( 133 ) + a*JVS( 460 )
  a = -W( 105 ) / JVS( 472 )
  W( 105 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 473 )
  W( 136 ) = W( 136 ) + a*JVS( 474 )
  a = -W( 106 ) / JVS( 476 )
  W( 106 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 477 )
  W( 129 ) = W( 129 ) + a*JVS( 478 )
  W( 133 ) = W( 133 ) + a*JVS( 479 )
  W( 136 ) = W( 136 ) + a*JVS( 480 )
  a = -W( 107 ) / JVS( 483 )
  W( 107 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 484 )
  W( 121 ) = W( 121 ) + a*JVS( 485 )
  W( 123 ) = W( 123 ) + a*JVS( 486 )
  W( 125 ) = W( 125 ) + a*JVS( 487 )
  W( 129 ) = W( 129 ) + a*JVS( 488 )
  W( 131 ) = W( 131 ) + a*JVS( 489 )
  W( 132 ) = W( 132 ) + a*JVS( 490 )
  W( 136 ) = W( 136 ) + a*JVS( 491 )
  a = -W( 111 ) / JVS( 531 )
  W( 111 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 532 )
  W( 132 ) = W( 132 ) + a*JVS( 533 )
  W( 133 ) = W( 133 ) + a*JVS( 534 )
  W( 136 ) = W( 136 ) + a*JVS( 535 )
  a = -W( 112 ) / JVS( 536 )
  W( 112 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 537 )
  W( 131 ) = W( 131 ) + a*JVS( 538 )
  W( 132 ) = W( 132 ) + a*JVS( 539 )
  W( 133 ) = W( 133 ) + a*JVS( 540 )
  a = -W( 113 ) / JVS( 541 )
  W( 113 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 542 )
  W( 131 ) = W( 131 ) + a*JVS( 543 )
  W( 132 ) = W( 132 ) + a*JVS( 544 )
  W( 133 ) = W( 133 ) + a*JVS( 545 )
  a = -W( 115 ) / JVS( 550 )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 136 ) = W( 136 ) + a*JVS( 553 )
  a = -W( 116 ) / JVS( 562 )
  W( 116 ) = -a
  W( 119 ) = W( 119 ) + a*JVS( 563 )
  W( 120 ) = W( 120 ) + a*JVS( 564 )
  W( 124 ) = W( 124 ) + a*JVS( 565 )
  W( 125 ) = W( 125 ) + a*JVS( 566 )
  W( 127 ) = W( 127 ) + a*JVS( 567 )
  W( 131 ) = W( 131 ) + a*JVS( 568 )
  W( 132 ) = W( 132 ) + a*JVS( 569 )
  W( 133 ) = W( 133 ) + a*JVS( 570 )
  W( 134 ) = W( 134 ) + a*JVS( 571 )
  W( 136 ) = W( 136 ) + a*JVS( 572 )
  a = -W( 117 ) / JVS( 580 )
  W( 117 ) = -a
  W( 119 ) = W( 119 ) + a*JVS( 581 )
  W( 121 ) = W( 121 ) + a*JVS( 582 )
  W( 125 ) = W( 125 ) + a*JVS( 583 )
  W( 126 ) = W( 126 ) + a*JVS( 584 )
  W( 127 ) = W( 127 ) + a*JVS( 585 )
  W( 128 ) = W( 128 ) + a*JVS( 586 )
  W( 131 ) = W( 131 ) + a*JVS( 587 )
  W( 132 ) = W( 132 ) + a*JVS( 588 )
  W( 133 ) = W( 133 ) + a*JVS( 589 )
  W( 135 ) = W( 135 ) + a*JVS( 590 )
  W( 136 ) = W( 136 ) + a*JVS( 591 )
  a = -W( 119 ) / JVS( 605 )
  W( 119 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 606 )
  W( 131 ) = W( 131 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 120 ) / JVS( 611 )
  W( 120 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 612 )
  W( 127 ) = W( 127 ) + a*JVS( 613 )
  W( 131 ) = W( 131 ) + a*JVS( 614 )
  W( 132 ) = W( 132 ) + a*JVS( 615 )
  W( 133 ) = W( 133 ) + a*JVS( 616 )
  W( 136 ) = W( 136 ) + a*JVS( 617 )
  a = -W( 121 ) / JVS( 623 )
  W( 121 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 624 )
  W( 127 ) = W( 127 ) + a*JVS( 625 )
  W( 131 ) = W( 131 ) + a*JVS( 626 )
  W( 132 ) = W( 132 ) + a*JVS( 627 )
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 135 ) = W( 135 ) + a*JVS( 630 )
  W( 136 ) = W( 136 ) + a*JVS( 631 )
  a = -W( 122 ) / JVS( 636 )
  W( 122 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 637 )
  W( 127 ) = W( 127 ) + a*JVS( 638 )
  W( 128 ) = W( 128 ) + a*JVS( 639 )
  W( 131 ) = W( 131 ) + a*JVS( 640 )
  W( 132 ) = W( 132 ) + a*JVS( 641 )
  W( 133 ) = W( 133 ) + a*JVS( 642 )
  W( 135 ) = W( 135 ) + a*JVS( 643 )
  W( 136 ) = W( 136 ) + a*JVS( 644 )
  a = -W( 123 ) / JVS( 656 )
  W( 123 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 657 )
  W( 126 ) = W( 126 ) + a*JVS( 658 )
  W( 127 ) = W( 127 ) + a*JVS( 659 )
  W( 128 ) = W( 128 ) + a*JVS( 660 )
  W( 131 ) = W( 131 ) + a*JVS( 661 )
  W( 132 ) = W( 132 ) + a*JVS( 662 )
  W( 133 ) = W( 133 ) + a*JVS( 663 )
  W( 135 ) = W( 135 ) + a*JVS( 664 )
  W( 136 ) = W( 136 ) + a*JVS( 665 )
  a = -W( 124 ) / JVS( 672 )
  W( 124 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 673 )
  W( 126 ) = W( 126 ) + a*JVS( 674 )
  W( 127 ) = W( 127 ) + a*JVS( 675 )
  W( 128 ) = W( 128 ) + a*JVS( 676 )
  W( 130 ) = W( 130 ) + a*JVS( 677 )
  W( 131 ) = W( 131 ) + a*JVS( 678 )
  W( 132 ) = W( 132 ) + a*JVS( 679 )
  W( 133 ) = W( 133 ) + a*JVS( 680 )
  W( 134 ) = W( 134 ) + a*JVS( 681 )
  W( 135 ) = W( 135 ) + a*JVS( 682 )
  W( 136 ) = W( 136 ) + a*JVS( 683 )
  a = -W( 125 ) / JVS( 687 )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 688 )
  W( 131 ) = W( 131 ) + a*JVS( 689 )
  W( 132 ) = W( 132 ) + a*JVS( 690 )
  W( 133 ) = W( 133 ) + a*JVS( 691 )
  W( 136 ) = W( 136 ) + a*JVS( 692 )
  a = -W( 126 ) / JVS( 701 )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 702 )
  W( 128 ) = W( 128 ) + a*JVS( 703 )
  W( 129 ) = W( 129 ) + a*JVS( 704 )
  W( 131 ) = W( 131 ) + a*JVS( 705 )
  W( 132 ) = W( 132 ) + a*JVS( 706 )
  W( 133 ) = W( 133 ) + a*JVS( 707 )
  W( 136 ) = W( 136 ) + a*JVS( 708 )
  a = -W( 127 ) / JVS( 718 )
  W( 127 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 719 )
  W( 130 ) = W( 130 ) + a*JVS( 720 )
  W( 131 ) = W( 131 ) + a*JVS( 721 )
  W( 132 ) = W( 132 ) + a*JVS( 722 )
  W( 133 ) = W( 133 ) + a*JVS( 723 )
  W( 136 ) = W( 136 ) + a*JVS( 724 )
  a = -W( 128 ) / JVS( 735 )
  W( 128 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 736 )
  W( 130 ) = W( 130 ) + a*JVS( 737 )
  W( 131 ) = W( 131 ) + a*JVS( 738 )
  W( 132 ) = W( 132 ) + a*JVS( 739 )
  W( 133 ) = W( 133 ) + a*JVS( 740 )
  W( 135 ) = W( 135 ) + a*JVS( 741 )
  W( 136 ) = W( 136 ) + a*JVS( 742 )
  JVS( 743) = W( 87 )
  JVS( 744) = W( 90 )
  JVS( 745) = W( 92 )
  JVS( 746) = W( 96 )
  JVS( 747) = W( 97 )
  JVS( 748) = W( 102 )
  JVS( 749) = W( 105 )
  JVS( 750) = W( 106 )
  JVS( 751) = W( 107 )
  JVS( 752) = W( 111 )
  JVS( 753) = W( 112 )
  JVS( 754) = W( 113 )
  JVS( 755) = W( 115 )
  JVS( 756) = W( 116 )
  JVS( 757) = W( 117 )
  JVS( 758) = W( 119 )
  JVS( 759) = W( 120 )
  JVS( 760) = W( 121 )
  JVS( 761) = W( 122 )
  JVS( 762) = W( 123 )
  JVS( 763) = W( 124 )
  JVS( 764) = W( 125 )
  JVS( 765) = W( 126 )
  JVS( 766) = W( 127 )
  JVS( 767) = W( 128 )
  JVS( 768) = W( 129 )
  JVS( 769) = W( 130 )
  JVS( 770) = W( 131 )
  JVS( 771) = W( 132 )
  JVS( 772) = W( 133 )
  JVS( 773) = W( 134 )
  JVS( 774) = W( 135 )
  JVS( 775) = W( 136 )
  IF ( ABS( JVS( 789 )) < TINY(a) ) THEN
         IER = 130
         RETURN
  END IF
   W( 87 ) = JVS( 776 )
   W( 110 ) = JVS( 777 )
   W( 114 ) = JVS( 778 )
   W( 115 ) = JVS( 779 )
   W( 118 ) = JVS( 780 )
   W( 119 ) = JVS( 781 )
   W( 121 ) = JVS( 782 )
   W( 123 ) = JVS( 783 )
   W( 125 ) = JVS( 784 )
   W( 126 ) = JVS( 785 )
   W( 127 ) = JVS( 786 )
   W( 128 ) = JVS( 787 )
   W( 129 ) = JVS( 788 )
   W( 130 ) = JVS( 789 )
   W( 131 ) = JVS( 790 )
   W( 132 ) = JVS( 791 )
   W( 133 ) = JVS( 792 )
   W( 134 ) = JVS( 793 )
   W( 135 ) = JVS( 794 )
   W( 136 ) = JVS( 795 )
  a = -W( 87 ) / JVS( 390 )
  W( 87 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 391 )
  W( 130 ) = W( 130 ) + a*JVS( 392 )
  a = -W( 110 ) / JVS( 526 )
  W( 110 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 527 )
  W( 131 ) = W( 131 ) + a*JVS( 528 )
  W( 133 ) = W( 133 ) + a*JVS( 529 )
  W( 136 ) = W( 136 ) + a*JVS( 530 )
  a = -W( 114 ) / JVS( 546 )
  W( 114 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 547 )
  W( 131 ) = W( 131 ) + a*JVS( 548 )
  W( 136 ) = W( 136 ) + a*JVS( 549 )
  a = -W( 115 ) / JVS( 550 )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 136 ) = W( 136 ) + a*JVS( 553 )
  a = -W( 118 ) / JVS( 593 )
  W( 118 ) = -a
  W( 119 ) = W( 119 ) + a*JVS( 594 )
  W( 125 ) = W( 125 ) + a*JVS( 595 )
  W( 126 ) = W( 126 ) + a*JVS( 596 )
  W( 127 ) = W( 127 ) + a*JVS( 597 )
  W( 128 ) = W( 128 ) + a*JVS( 598 )
  W( 131 ) = W( 131 ) + a*JVS( 599 )
  W( 132 ) = W( 132 ) + a*JVS( 600 )
  W( 133 ) = W( 133 ) + a*JVS( 601 )
  W( 134 ) = W( 134 ) + a*JVS( 602 )
  W( 135 ) = W( 135 ) + a*JVS( 603 )
  W( 136 ) = W( 136 ) + a*JVS( 604 )
  a = -W( 119 ) / JVS( 605 )
  W( 119 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 606 )
  W( 131 ) = W( 131 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 121 ) / JVS( 623 )
  W( 121 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 624 )
  W( 127 ) = W( 127 ) + a*JVS( 625 )
  W( 131 ) = W( 131 ) + a*JVS( 626 )
  W( 132 ) = W( 132 ) + a*JVS( 627 )
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 135 ) = W( 135 ) + a*JVS( 630 )
  W( 136 ) = W( 136 ) + a*JVS( 631 )
  a = -W( 123 ) / JVS( 656 )
  W( 123 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 657 )
  W( 126 ) = W( 126 ) + a*JVS( 658 )
  W( 127 ) = W( 127 ) + a*JVS( 659 )
  W( 128 ) = W( 128 ) + a*JVS( 660 )
  W( 131 ) = W( 131 ) + a*JVS( 661 )
  W( 132 ) = W( 132 ) + a*JVS( 662 )
  W( 133 ) = W( 133 ) + a*JVS( 663 )
  W( 135 ) = W( 135 ) + a*JVS( 664 )
  W( 136 ) = W( 136 ) + a*JVS( 665 )
  a = -W( 125 ) / JVS( 687 )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 688 )
  W( 131 ) = W( 131 ) + a*JVS( 689 )
  W( 132 ) = W( 132 ) + a*JVS( 690 )
  W( 133 ) = W( 133 ) + a*JVS( 691 )
  W( 136 ) = W( 136 ) + a*JVS( 692 )
  a = -W( 126 ) / JVS( 701 )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 702 )
  W( 128 ) = W( 128 ) + a*JVS( 703 )
  W( 129 ) = W( 129 ) + a*JVS( 704 )
  W( 131 ) = W( 131 ) + a*JVS( 705 )
  W( 132 ) = W( 132 ) + a*JVS( 706 )
  W( 133 ) = W( 133 ) + a*JVS( 707 )
  W( 136 ) = W( 136 ) + a*JVS( 708 )
  a = -W( 127 ) / JVS( 718 )
  W( 127 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 719 )
  W( 130 ) = W( 130 ) + a*JVS( 720 )
  W( 131 ) = W( 131 ) + a*JVS( 721 )
  W( 132 ) = W( 132 ) + a*JVS( 722 )
  W( 133 ) = W( 133 ) + a*JVS( 723 )
  W( 136 ) = W( 136 ) + a*JVS( 724 )
  a = -W( 128 ) / JVS( 735 )
  W( 128 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 736 )
  W( 130 ) = W( 130 ) + a*JVS( 737 )
  W( 131 ) = W( 131 ) + a*JVS( 738 )
  W( 132 ) = W( 132 ) + a*JVS( 739 )
  W( 133 ) = W( 133 ) + a*JVS( 740 )
  W( 135 ) = W( 135 ) + a*JVS( 741 )
  W( 136 ) = W( 136 ) + a*JVS( 742 )
  a = -W( 129 ) / JVS( 768 )
  W( 129 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 769 )
  W( 131 ) = W( 131 ) + a*JVS( 770 )
  W( 132 ) = W( 132 ) + a*JVS( 771 )
  W( 133 ) = W( 133 ) + a*JVS( 772 )
  W( 134 ) = W( 134 ) + a*JVS( 773 )
  W( 135 ) = W( 135 ) + a*JVS( 774 )
  W( 136 ) = W( 136 ) + a*JVS( 775 )
  JVS( 776) = W( 87 )
  JVS( 777) = W( 110 )
  JVS( 778) = W( 114 )
  JVS( 779) = W( 115 )
  JVS( 780) = W( 118 )
  JVS( 781) = W( 119 )
  JVS( 782) = W( 121 )
  JVS( 783) = W( 123 )
  JVS( 784) = W( 125 )
  JVS( 785) = W( 126 )
  JVS( 786) = W( 127 )
  JVS( 787) = W( 128 )
  JVS( 788) = W( 129 )
  JVS( 789) = W( 130 )
  JVS( 790) = W( 131 )
  JVS( 791) = W( 132 )
  JVS( 792) = W( 133 )
  JVS( 793) = W( 134 )
  JVS( 794) = W( 135 )
  JVS( 795) = W( 136 )
  IF ( ABS( JVS( 904 )) < TINY(a) ) THEN
         IER = 131
         RETURN
  END IF
   W( 4 ) = JVS( 796 )
   W( 5 ) = JVS( 797 )
   W( 6 ) = JVS( 798 )
   W( 7 ) = JVS( 799 )
   W( 20 ) = JVS( 800 )
   W( 21 ) = JVS( 801 )
   W( 22 ) = JVS( 802 )
   W( 23 ) = JVS( 803 )
   W( 24 ) = JVS( 804 )
   W( 25 ) = JVS( 805 )
   W( 26 ) = JVS( 806 )
   W( 27 ) = JVS( 807 )
   W( 28 ) = JVS( 808 )
   W( 29 ) = JVS( 809 )
   W( 30 ) = JVS( 810 )
   W( 31 ) = JVS( 811 )
   W( 32 ) = JVS( 812 )
   W( 33 ) = JVS( 813 )
   W( 34 ) = JVS( 814 )
   W( 35 ) = JVS( 815 )
   W( 36 ) = JVS( 816 )
   W( 37 ) = JVS( 817 )
   W( 38 ) = JVS( 818 )
   W( 39 ) = JVS( 819 )
   W( 40 ) = JVS( 820 )
   W( 41 ) = JVS( 821 )
   W( 42 ) = JVS( 822 )
   W( 43 ) = JVS( 823 )
   W( 44 ) = JVS( 824 )
   W( 45 ) = JVS( 825 )
   W( 46 ) = JVS( 826 )
   W( 47 ) = JVS( 827 )
   W( 48 ) = JVS( 828 )
   W( 49 ) = JVS( 829 )
   W( 50 ) = JVS( 830 )
   W( 51 ) = JVS( 831 )
   W( 52 ) = JVS( 832 )
   W( 53 ) = JVS( 833 )
   W( 54 ) = JVS( 834 )
   W( 55 ) = JVS( 835 )
   W( 56 ) = JVS( 836 )
   W( 57 ) = JVS( 837 )
   W( 58 ) = JVS( 838 )
   W( 59 ) = JVS( 839 )
   W( 60 ) = JVS( 840 )
   W( 61 ) = JVS( 841 )
   W( 62 ) = JVS( 842 )
   W( 63 ) = JVS( 843 )
   W( 64 ) = JVS( 844 )
   W( 65 ) = JVS( 845 )
   W( 66 ) = JVS( 846 )
   W( 67 ) = JVS( 847 )
   W( 68 ) = JVS( 848 )
   W( 69 ) = JVS( 849 )
   W( 70 ) = JVS( 850 )
   W( 71 ) = JVS( 851 )
   W( 72 ) = JVS( 852 )
   W( 73 ) = JVS( 853 )
   W( 74 ) = JVS( 854 )
   W( 75 ) = JVS( 855 )
   W( 76 ) = JVS( 856 )
   W( 77 ) = JVS( 857 )
   W( 78 ) = JVS( 858 )
   W( 79 ) = JVS( 859 )
   W( 80 ) = JVS( 860 )
   W( 81 ) = JVS( 861 )
   W( 82 ) = JVS( 862 )
   W( 83 ) = JVS( 863 )
   W( 84 ) = JVS( 864 )
   W( 85 ) = JVS( 865 )
   W( 86 ) = JVS( 866 )
   W( 88 ) = JVS( 867 )
   W( 89 ) = JVS( 868 )
   W( 91 ) = JVS( 869 )
   W( 93 ) = JVS( 870 )
   W( 94 ) = JVS( 871 )
   W( 95 ) = JVS( 872 )
   W( 96 ) = JVS( 873 )
   W( 98 ) = JVS( 874 )
   W( 100 ) = JVS( 875 )
   W( 101 ) = JVS( 876 )
   W( 102 ) = JVS( 877 )
   W( 103 ) = JVS( 878 )
   W( 104 ) = JVS( 879 )
   W( 105 ) = JVS( 880 )
   W( 107 ) = JVS( 881 )
   W( 108 ) = JVS( 882 )
   W( 109 ) = JVS( 883 )
   W( 110 ) = JVS( 884 )
   W( 111 ) = JVS( 885 )
   W( 112 ) = JVS( 886 )
   W( 113 ) = JVS( 887 )
   W( 114 ) = JVS( 888 )
   W( 115 ) = JVS( 889 )
   W( 116 ) = JVS( 890 )
   W( 118 ) = JVS( 891 )
   W( 119 ) = JVS( 892 )
   W( 120 ) = JVS( 893 )
   W( 121 ) = JVS( 894 )
   W( 122 ) = JVS( 895 )
   W( 123 ) = JVS( 896 )
   W( 124 ) = JVS( 897 )
   W( 125 ) = JVS( 898 )
   W( 126 ) = JVS( 899 )
   W( 127 ) = JVS( 900 )
   W( 128 ) = JVS( 901 )
   W( 129 ) = JVS( 902 )
   W( 130 ) = JVS( 903 )
   W( 131 ) = JVS( 904 )
   W( 132 ) = JVS( 905 )
   W( 133 ) = JVS( 906 )
   W( 134 ) = JVS( 907 )
   W( 135 ) = JVS( 908 )
   W( 136 ) = JVS( 909 )
  a = -W( 4 ) / JVS( 6 )
  W( 4 ) = -a
  a = -W( 5 ) / JVS( 7 )
  W( 5 ) = -a
  a = -W( 6 ) / JVS( 8 )
  W( 6 ) = -a
  a = -W( 7 ) / JVS( 9 )
  W( 7 ) = -a
  a = -W( 20 ) / JVS( 166 )
  W( 20 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 167 )
  a = -W( 21 ) / JVS( 168 )
  W( 21 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 169 )
  a = -W( 22 ) / JVS( 170 )
  W( 22 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 171 )
  a = -W( 23 ) / JVS( 172 )
  W( 23 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 173 )
  a = -W( 24 ) / JVS( 174 )
  W( 24 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 175 )
  a = -W( 25 ) / JVS( 176 )
  W( 25 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 177 )
  a = -W( 26 ) / JVS( 178 )
  W( 26 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 179 )
  a = -W( 27 ) / JVS( 180 )
  W( 27 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 181 )
  a = -W( 28 ) / JVS( 182 )
  W( 28 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 183 )
  a = -W( 29 ) / JVS( 184 )
  W( 29 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 185 )
  a = -W( 30 ) / JVS( 186 )
  W( 30 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 187 )
  a = -W( 31 ) / JVS( 188 )
  W( 31 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 189 )
  a = -W( 32 ) / JVS( 190 )
  W( 32 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 191 )
  a = -W( 33 ) / JVS( 192 )
  W( 33 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 193 )
  a = -W( 34 ) / JVS( 194 )
  W( 34 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 195 )
  a = -W( 35 ) / JVS( 196 )
  W( 35 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 197 )
  a = -W( 36 ) / JVS( 198 )
  W( 36 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 199 )
  a = -W( 37 ) / JVS( 200 )
  W( 37 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 201 )
  a = -W( 38 ) / JVS( 202 )
  W( 38 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 203 )
  a = -W( 39 ) / JVS( 205 )
  W( 39 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 206 )
  W( 131 ) = W( 131 ) + a*JVS( 207 )
  a = -W( 40 ) / JVS( 208 )
  W( 40 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 209 )
  a = -W( 41 ) / JVS( 210 )
  W( 41 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 211 )
  a = -W( 42 ) / JVS( 212 )
  W( 42 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 213 )
  a = -W( 43 ) / JVS( 214 )
  W( 43 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 215 )
  a = -W( 44 ) / JVS( 216 )
  W( 44 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 217 )
  a = -W( 45 ) / JVS( 218 )
  W( 45 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 219 )
  a = -W( 46 ) / JVS( 220 )
  W( 46 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 221 )
  a = -W( 47 ) / JVS( 224 )
  W( 47 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 225 )
  a = -W( 48 ) / JVS( 226 )
  W( 48 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 227 )
  a = -W( 49 ) / JVS( 229 )
  W( 49 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 230 )
  a = -W( 50 ) / JVS( 235 )
  W( 50 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 236 )
  a = -W( 51 ) / JVS( 239 )
  W( 51 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 240 )
  a = -W( 52 ) / JVS( 245 )
  W( 52 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 246 )
  a = -W( 53 ) / JVS( 249 )
  W( 53 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 250 )
  a = -W( 54 ) / JVS( 255 )
  W( 54 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 256 )
  a = -W( 55 ) / JVS( 259 )
  W( 55 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 260 )
  a = -W( 56 ) / JVS( 265 )
  W( 56 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 266 )
  a = -W( 57 ) / JVS( 269 )
  W( 57 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 270 )
  a = -W( 58 ) / JVS( 275 )
  W( 58 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 276 )
  a = -W( 59 ) / JVS( 279 )
  W( 59 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 280 )
  a = -W( 60 ) / JVS( 283 )
  W( 60 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 284 )
  a = -W( 61 ) / JVS( 288 )
  W( 61 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 289 )
  W( 131 ) = W( 131 ) + a*JVS( 290 )
  a = -W( 62 ) / JVS( 295 )
  W( 62 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 296 )
  a = -W( 63 ) / JVS( 297 )
  W( 63 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 298 )
  a = -W( 64 ) / JVS( 300 )
  W( 64 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 301 )
  W( 131 ) = W( 131 ) + a*JVS( 302 )
  a = -W( 65 ) / JVS( 303 )
  W( 65 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 304 )
  a = -W( 66 ) / JVS( 305 )
  W( 66 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 306 )
  a = -W( 67 ) / JVS( 307 )
  W( 67 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 308 )
  a = -W( 68 ) / JVS( 309 )
  W( 68 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 310 )
  a = -W( 69 ) / JVS( 311 )
  W( 69 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 312 )
  a = -W( 70 ) / JVS( 313 )
  W( 70 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 314 )
  a = -W( 71 ) / JVS( 315 )
  W( 71 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 316 )
  a = -W( 72 ) / JVS( 319 )
  W( 72 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 320 )
  a = -W( 73 ) / JVS( 322 )
  W( 73 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 323 )
  a = -W( 74 ) / JVS( 328 )
  W( 74 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 329 )
  a = -W( 75 ) / JVS( 332 )
  W( 75 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 333 )
  a = -W( 76 ) / JVS( 338 )
  W( 76 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 339 )
  a = -W( 77 ) / JVS( 342 )
  W( 77 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 343 )
  a = -W( 78 ) / JVS( 348 )
  W( 78 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 349 )
  a = -W( 79 ) / JVS( 352 )
  W( 79 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 353 )
  a = -W( 80 ) / JVS( 358 )
  W( 80 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 359 )
  a = -W( 81 ) / JVS( 362 )
  W( 81 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 363 )
  a = -W( 82 ) / JVS( 368 )
  W( 82 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 369 )
  a = -W( 83 ) / JVS( 372 )
  W( 83 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 373 )
  a = -W( 84 ) / JVS( 376 )
  W( 84 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 377 )
  a = -W( 85 ) / JVS( 381 )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 382 )
  W( 131 ) = W( 131 ) + a*JVS( 383 )
  a = -W( 86 ) / JVS( 388 )
  W( 86 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 389 )
  a = -W( 88 ) / JVS( 393 )
  W( 88 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 394 )
  W( 132 ) = W( 132 ) + a*JVS( 395 )
  a = -W( 89 ) / JVS( 396 )
  W( 89 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 397 )
  a = -W( 91 ) / JVS( 401 )
  W( 91 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 402 )
  a = -W( 93 ) / JVS( 408 )
  W( 93 ) = -a
  W( 114 ) = W( 114 ) + a*JVS( 409 )
  W( 119 ) = W( 119 ) + a*JVS( 410 )
  W( 127 ) = W( 127 ) + a*JVS( 411 )
  W( 131 ) = W( 131 ) + a*JVS( 412 )
  a = -W( 94 ) / JVS( 413 )
  W( 94 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 414 )
  W( 131 ) = W( 131 ) + a*JVS( 415 )
  W( 136 ) = W( 136 ) + a*JVS( 416 )
  a = -W( 95 ) / JVS( 417 )
  W( 95 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 418 )
  W( 131 ) = W( 131 ) + a*JVS( 419 )
  W( 136 ) = W( 136 ) + a*JVS( 420 )
  a = -W( 96 ) / JVS( 421 )
  W( 96 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 422 )
  W( 131 ) = W( 131 ) + a*JVS( 423 )
  W( 132 ) = W( 132 ) + a*JVS( 424 )
  a = -W( 98 ) / JVS( 430 )
  W( 98 ) = -a
  W( 114 ) = W( 114 ) + a*JVS( 431 )
  W( 119 ) = W( 119 ) + a*JVS( 432 )
  W( 122 ) = W( 122 ) + a*JVS( 433 )
  W( 127 ) = W( 127 ) + a*JVS( 434 )
  W( 131 ) = W( 131 ) + a*JVS( 435 )
  a = -W( 100 ) / JVS( 448 )
  W( 100 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 449 )
  W( 131 ) = W( 131 ) + a*JVS( 450 )
  W( 132 ) = W( 132 ) + a*JVS( 451 )
  a = -W( 101 ) / JVS( 452 )
  W( 101 ) = -a
  W( 122 ) = W( 122 ) + a*JVS( 453 )
  W( 131 ) = W( 131 ) + a*JVS( 454 )
  W( 132 ) = W( 132 ) + a*JVS( 455 )
  a = -W( 102 ) / JVS( 456 )
  W( 102 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 457 )
  W( 131 ) = W( 131 ) + a*JVS( 458 )
  W( 132 ) = W( 132 ) + a*JVS( 459 )
  W( 133 ) = W( 133 ) + a*JVS( 460 )
  a = -W( 103 ) / JVS( 461 )
  W( 103 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 462 )
  W( 131 ) = W( 131 ) + a*JVS( 463 )
  a = -W( 104 ) / JVS( 464 )
  W( 104 ) = -a
  W( 114 ) = W( 114 ) + a*JVS( 465 )
  W( 119 ) = W( 119 ) + a*JVS( 466 )
  W( 124 ) = W( 124 ) + a*JVS( 467 )
  W( 127 ) = W( 127 ) + a*JVS( 468 )
  W( 131 ) = W( 131 ) + a*JVS( 469 )
  a = -W( 105 ) / JVS( 472 )
  W( 105 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 473 )
  W( 136 ) = W( 136 ) + a*JVS( 474 )
  a = -W( 107 ) / JVS( 483 )
  W( 107 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 484 )
  W( 121 ) = W( 121 ) + a*JVS( 485 )
  W( 123 ) = W( 123 ) + a*JVS( 486 )
  W( 125 ) = W( 125 ) + a*JVS( 487 )
  W( 129 ) = W( 129 ) + a*JVS( 488 )
  W( 131 ) = W( 131 ) + a*JVS( 489 )
  W( 132 ) = W( 132 ) + a*JVS( 490 )
  W( 136 ) = W( 136 ) + a*JVS( 491 )
  a = -W( 108 ) / JVS( 493 )
  W( 108 ) = -a
  W( 110 ) = W( 110 ) + a*JVS( 494 )
  W( 112 ) = W( 112 ) + a*JVS( 495 )
  W( 114 ) = W( 114 ) + a*JVS( 496 )
  W( 115 ) = W( 115 ) + a*JVS( 497 )
  W( 116 ) = W( 116 ) + a*JVS( 498 )
  W( 119 ) = W( 119 ) + a*JVS( 499 )
  W( 121 ) = W( 121 ) + a*JVS( 500 )
  W( 123 ) = W( 123 ) + a*JVS( 501 )
  W( 125 ) = W( 125 ) + a*JVS( 502 )
  W( 127 ) = W( 127 ) + a*JVS( 503 )
  W( 131 ) = W( 131 ) + a*JVS( 504 )
  W( 133 ) = W( 133 ) + a*JVS( 505 )
  W( 136 ) = W( 136 ) + a*JVS( 506 )
  a = -W( 109 ) / JVS( 509 )
  W( 109 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 510 )
  W( 113 ) = W( 113 ) + a*JVS( 511 )
  W( 114 ) = W( 114 ) + a*JVS( 512 )
  W( 119 ) = W( 119 ) + a*JVS( 513 )
  W( 120 ) = W( 120 ) + a*JVS( 514 )
  W( 125 ) = W( 125 ) + a*JVS( 515 )
  W( 126 ) = W( 126 ) + a*JVS( 516 )
  W( 127 ) = W( 127 ) + a*JVS( 517 )
  W( 128 ) = W( 128 ) + a*JVS( 518 )
  W( 131 ) = W( 131 ) + a*JVS( 519 )
  W( 132 ) = W( 132 ) + a*JVS( 520 )
  W( 133 ) = W( 133 ) + a*JVS( 521 )
  W( 135 ) = W( 135 ) + a*JVS( 522 )
  W( 136 ) = W( 136 ) + a*JVS( 523 )
  a = -W( 110 ) / JVS( 526 )
  W( 110 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 527 )
  W( 131 ) = W( 131 ) + a*JVS( 528 )
  W( 133 ) = W( 133 ) + a*JVS( 529 )
  W( 136 ) = W( 136 ) + a*JVS( 530 )
  a = -W( 111 ) / JVS( 531 )
  W( 111 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 532 )
  W( 132 ) = W( 132 ) + a*JVS( 533 )
  W( 133 ) = W( 133 ) + a*JVS( 534 )
  W( 136 ) = W( 136 ) + a*JVS( 535 )
  a = -W( 112 ) / JVS( 536 )
  W( 112 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 537 )
  W( 131 ) = W( 131 ) + a*JVS( 538 )
  W( 132 ) = W( 132 ) + a*JVS( 539 )
  W( 133 ) = W( 133 ) + a*JVS( 540 )
  a = -W( 113 ) / JVS( 541 )
  W( 113 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 542 )
  W( 131 ) = W( 131 ) + a*JVS( 543 )
  W( 132 ) = W( 132 ) + a*JVS( 544 )
  W( 133 ) = W( 133 ) + a*JVS( 545 )
  a = -W( 114 ) / JVS( 546 )
  W( 114 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 547 )
  W( 131 ) = W( 131 ) + a*JVS( 548 )
  W( 136 ) = W( 136 ) + a*JVS( 549 )
  a = -W( 115 ) / JVS( 550 )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 136 ) = W( 136 ) + a*JVS( 553 )
  a = -W( 116 ) / JVS( 562 )
  W( 116 ) = -a
  W( 119 ) = W( 119 ) + a*JVS( 563 )
  W( 120 ) = W( 120 ) + a*JVS( 564 )
  W( 124 ) = W( 124 ) + a*JVS( 565 )
  W( 125 ) = W( 125 ) + a*JVS( 566 )
  W( 127 ) = W( 127 ) + a*JVS( 567 )
  W( 131 ) = W( 131 ) + a*JVS( 568 )
  W( 132 ) = W( 132 ) + a*JVS( 569 )
  W( 133 ) = W( 133 ) + a*JVS( 570 )
  W( 134 ) = W( 134 ) + a*JVS( 571 )
  W( 136 ) = W( 136 ) + a*JVS( 572 )
  a = -W( 118 ) / JVS( 593 )
  W( 118 ) = -a
  W( 119 ) = W( 119 ) + a*JVS( 594 )
  W( 125 ) = W( 125 ) + a*JVS( 595 )
  W( 126 ) = W( 126 ) + a*JVS( 596 )
  W( 127 ) = W( 127 ) + a*JVS( 597 )
  W( 128 ) = W( 128 ) + a*JVS( 598 )
  W( 131 ) = W( 131 ) + a*JVS( 599 )
  W( 132 ) = W( 132 ) + a*JVS( 600 )
  W( 133 ) = W( 133 ) + a*JVS( 601 )
  W( 134 ) = W( 134 ) + a*JVS( 602 )
  W( 135 ) = W( 135 ) + a*JVS( 603 )
  W( 136 ) = W( 136 ) + a*JVS( 604 )
  a = -W( 119 ) / JVS( 605 )
  W( 119 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 606 )
  W( 131 ) = W( 131 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 120 ) / JVS( 611 )
  W( 120 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 612 )
  W( 127 ) = W( 127 ) + a*JVS( 613 )
  W( 131 ) = W( 131 ) + a*JVS( 614 )
  W( 132 ) = W( 132 ) + a*JVS( 615 )
  W( 133 ) = W( 133 ) + a*JVS( 616 )
  W( 136 ) = W( 136 ) + a*JVS( 617 )
  a = -W( 121 ) / JVS( 623 )
  W( 121 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 624 )
  W( 127 ) = W( 127 ) + a*JVS( 625 )
  W( 131 ) = W( 131 ) + a*JVS( 626 )
  W( 132 ) = W( 132 ) + a*JVS( 627 )
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 135 ) = W( 135 ) + a*JVS( 630 )
  W( 136 ) = W( 136 ) + a*JVS( 631 )
  a = -W( 122 ) / JVS( 636 )
  W( 122 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 637 )
  W( 127 ) = W( 127 ) + a*JVS( 638 )
  W( 128 ) = W( 128 ) + a*JVS( 639 )
  W( 131 ) = W( 131 ) + a*JVS( 640 )
  W( 132 ) = W( 132 ) + a*JVS( 641 )
  W( 133 ) = W( 133 ) + a*JVS( 642 )
  W( 135 ) = W( 135 ) + a*JVS( 643 )
  W( 136 ) = W( 136 ) + a*JVS( 644 )
  a = -W( 123 ) / JVS( 656 )
  W( 123 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 657 )
  W( 126 ) = W( 126 ) + a*JVS( 658 )
  W( 127 ) = W( 127 ) + a*JVS( 659 )
  W( 128 ) = W( 128 ) + a*JVS( 660 )
  W( 131 ) = W( 131 ) + a*JVS( 661 )
  W( 132 ) = W( 132 ) + a*JVS( 662 )
  W( 133 ) = W( 133 ) + a*JVS( 663 )
  W( 135 ) = W( 135 ) + a*JVS( 664 )
  W( 136 ) = W( 136 ) + a*JVS( 665 )
  a = -W( 124 ) / JVS( 672 )
  W( 124 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 673 )
  W( 126 ) = W( 126 ) + a*JVS( 674 )
  W( 127 ) = W( 127 ) + a*JVS( 675 )
  W( 128 ) = W( 128 ) + a*JVS( 676 )
  W( 130 ) = W( 130 ) + a*JVS( 677 )
  W( 131 ) = W( 131 ) + a*JVS( 678 )
  W( 132 ) = W( 132 ) + a*JVS( 679 )
  W( 133 ) = W( 133 ) + a*JVS( 680 )
  W( 134 ) = W( 134 ) + a*JVS( 681 )
  W( 135 ) = W( 135 ) + a*JVS( 682 )
  W( 136 ) = W( 136 ) + a*JVS( 683 )
  a = -W( 125 ) / JVS( 687 )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 688 )
  W( 131 ) = W( 131 ) + a*JVS( 689 )
  W( 132 ) = W( 132 ) + a*JVS( 690 )
  W( 133 ) = W( 133 ) + a*JVS( 691 )
  W( 136 ) = W( 136 ) + a*JVS( 692 )
  a = -W( 126 ) / JVS( 701 )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 702 )
  W( 128 ) = W( 128 ) + a*JVS( 703 )
  W( 129 ) = W( 129 ) + a*JVS( 704 )
  W( 131 ) = W( 131 ) + a*JVS( 705 )
  W( 132 ) = W( 132 ) + a*JVS( 706 )
  W( 133 ) = W( 133 ) + a*JVS( 707 )
  W( 136 ) = W( 136 ) + a*JVS( 708 )
  a = -W( 127 ) / JVS( 718 )
  W( 127 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 719 )
  W( 130 ) = W( 130 ) + a*JVS( 720 )
  W( 131 ) = W( 131 ) + a*JVS( 721 )
  W( 132 ) = W( 132 ) + a*JVS( 722 )
  W( 133 ) = W( 133 ) + a*JVS( 723 )
  W( 136 ) = W( 136 ) + a*JVS( 724 )
  a = -W( 128 ) / JVS( 735 )
  W( 128 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 736 )
  W( 130 ) = W( 130 ) + a*JVS( 737 )
  W( 131 ) = W( 131 ) + a*JVS( 738 )
  W( 132 ) = W( 132 ) + a*JVS( 739 )
  W( 133 ) = W( 133 ) + a*JVS( 740 )
  W( 135 ) = W( 135 ) + a*JVS( 741 )
  W( 136 ) = W( 136 ) + a*JVS( 742 )
  a = -W( 129 ) / JVS( 768 )
  W( 129 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 769 )
  W( 131 ) = W( 131 ) + a*JVS( 770 )
  W( 132 ) = W( 132 ) + a*JVS( 771 )
  W( 133 ) = W( 133 ) + a*JVS( 772 )
  W( 134 ) = W( 134 ) + a*JVS( 773 )
  W( 135 ) = W( 135 ) + a*JVS( 774 )
  W( 136 ) = W( 136 ) + a*JVS( 775 )
  a = -W( 130 ) / JVS( 789 )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 790 )
  W( 132 ) = W( 132 ) + a*JVS( 791 )
  W( 133 ) = W( 133 ) + a*JVS( 792 )
  W( 134 ) = W( 134 ) + a*JVS( 793 )
  W( 135 ) = W( 135 ) + a*JVS( 794 )
  W( 136 ) = W( 136 ) + a*JVS( 795 )
  JVS( 796) = W( 4 )
  JVS( 797) = W( 5 )
  JVS( 798) = W( 6 )
  JVS( 799) = W( 7 )
  JVS( 800) = W( 20 )
  JVS( 801) = W( 21 )
  JVS( 802) = W( 22 )
  JVS( 803) = W( 23 )
  JVS( 804) = W( 24 )
  JVS( 805) = W( 25 )
  JVS( 806) = W( 26 )
  JVS( 807) = W( 27 )
  JVS( 808) = W( 28 )
  JVS( 809) = W( 29 )
  JVS( 810) = W( 30 )
  JVS( 811) = W( 31 )
  JVS( 812) = W( 32 )
  JVS( 813) = W( 33 )
  JVS( 814) = W( 34 )
  JVS( 815) = W( 35 )
  JVS( 816) = W( 36 )
  JVS( 817) = W( 37 )
  JVS( 818) = W( 38 )
  JVS( 819) = W( 39 )
  JVS( 820) = W( 40 )
  JVS( 821) = W( 41 )
  JVS( 822) = W( 42 )
  JVS( 823) = W( 43 )
  JVS( 824) = W( 44 )
  JVS( 825) = W( 45 )
  JVS( 826) = W( 46 )
  JVS( 827) = W( 47 )
  JVS( 828) = W( 48 )
  JVS( 829) = W( 49 )
  JVS( 830) = W( 50 )
  JVS( 831) = W( 51 )
  JVS( 832) = W( 52 )
  JVS( 833) = W( 53 )
  JVS( 834) = W( 54 )
  JVS( 835) = W( 55 )
  JVS( 836) = W( 56 )
  JVS( 837) = W( 57 )
  JVS( 838) = W( 58 )
  JVS( 839) = W( 59 )
  JVS( 840) = W( 60 )
  JVS( 841) = W( 61 )
  JVS( 842) = W( 62 )
  JVS( 843) = W( 63 )
  JVS( 844) = W( 64 )
  JVS( 845) = W( 65 )
  JVS( 846) = W( 66 )
  JVS( 847) = W( 67 )
  JVS( 848) = W( 68 )
  JVS( 849) = W( 69 )
  JVS( 850) = W( 70 )
  JVS( 851) = W( 71 )
  JVS( 852) = W( 72 )
  JVS( 853) = W( 73 )
  JVS( 854) = W( 74 )
  JVS( 855) = W( 75 )
  JVS( 856) = W( 76 )
  JVS( 857) = W( 77 )
  JVS( 858) = W( 78 )
  JVS( 859) = W( 79 )
  JVS( 860) = W( 80 )
  JVS( 861) = W( 81 )
  JVS( 862) = W( 82 )
  JVS( 863) = W( 83 )
  JVS( 864) = W( 84 )
  JVS( 865) = W( 85 )
  JVS( 866) = W( 86 )
  JVS( 867) = W( 88 )
  JVS( 868) = W( 89 )
  JVS( 869) = W( 91 )
  JVS( 870) = W( 93 )
  JVS( 871) = W( 94 )
  JVS( 872) = W( 95 )
  JVS( 873) = W( 96 )
  JVS( 874) = W( 98 )
  JVS( 875) = W( 100 )
  JVS( 876) = W( 101 )
  JVS( 877) = W( 102 )
  JVS( 878) = W( 103 )
  JVS( 879) = W( 104 )
  JVS( 880) = W( 105 )
  JVS( 881) = W( 107 )
  JVS( 882) = W( 108 )
  JVS( 883) = W( 109 )
  JVS( 884) = W( 110 )
  JVS( 885) = W( 111 )
  JVS( 886) = W( 112 )
  JVS( 887) = W( 113 )
  JVS( 888) = W( 114 )
  JVS( 889) = W( 115 )
  JVS( 890) = W( 116 )
  JVS( 891) = W( 118 )
  JVS( 892) = W( 119 )
  JVS( 893) = W( 120 )
  JVS( 894) = W( 121 )
  JVS( 895) = W( 122 )
  JVS( 896) = W( 123 )
  JVS( 897) = W( 124 )
  JVS( 898) = W( 125 )
  JVS( 899) = W( 126 )
  JVS( 900) = W( 127 )
  JVS( 901) = W( 128 )
  JVS( 902) = W( 129 )
  JVS( 903) = W( 130 )
  JVS( 904) = W( 131 )
  JVS( 905) = W( 132 )
  JVS( 906) = W( 133 )
  JVS( 907) = W( 134 )
  JVS( 908) = W( 135 )
  JVS( 909) = W( 136 )
  IF ( ABS( JVS( 944 )) < TINY(a) ) THEN
         IER = 132
         RETURN
  END IF
   W( 37 ) = JVS( 910 )
   W( 48 ) = JVS( 911 )
   W( 88 ) = JVS( 912 )
   W( 89 ) = JVS( 913 )
   W( 91 ) = JVS( 914 )
   W( 96 ) = JVS( 915 )
   W( 97 ) = JVS( 916 )
   W( 100 ) = JVS( 917 )
   W( 101 ) = JVS( 918 )
   W( 103 ) = JVS( 919 )
   W( 104 ) = JVS( 920 )
   W( 105 ) = JVS( 921 )
   W( 108 ) = JVS( 922 )
   W( 110 ) = JVS( 923 )
   W( 111 ) = JVS( 924 )
   W( 112 ) = JVS( 925 )
   W( 113 ) = JVS( 926 )
   W( 114 ) = JVS( 927 )
   W( 115 ) = JVS( 928 )
   W( 116 ) = JVS( 929 )
   W( 117 ) = JVS( 930 )
   W( 119 ) = JVS( 931 )
   W( 120 ) = JVS( 932 )
   W( 121 ) = JVS( 933 )
   W( 122 ) = JVS( 934 )
   W( 123 ) = JVS( 935 )
   W( 124 ) = JVS( 936 )
   W( 125 ) = JVS( 937 )
   W( 126 ) = JVS( 938 )
   W( 127 ) = JVS( 939 )
   W( 128 ) = JVS( 940 )
   W( 129 ) = JVS( 941 )
   W( 130 ) = JVS( 942 )
   W( 131 ) = JVS( 943 )
   W( 132 ) = JVS( 944 )
   W( 133 ) = JVS( 945 )
   W( 134 ) = JVS( 946 )
   W( 135 ) = JVS( 947 )
   W( 136 ) = JVS( 948 )
  a = -W( 37 ) / JVS( 200 )
  W( 37 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 201 )
  a = -W( 48 ) / JVS( 226 )
  W( 48 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 227 )
  a = -W( 88 ) / JVS( 393 )
  W( 88 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 394 )
  W( 132 ) = W( 132 ) + a*JVS( 395 )
  a = -W( 89 ) / JVS( 396 )
  W( 89 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 397 )
  a = -W( 91 ) / JVS( 401 )
  W( 91 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 402 )
  a = -W( 96 ) / JVS( 421 )
  W( 96 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 422 )
  W( 131 ) = W( 131 ) + a*JVS( 423 )
  W( 132 ) = W( 132 ) + a*JVS( 424 )
  a = -W( 97 ) / JVS( 427 )
  W( 97 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 428 )
  W( 133 ) = W( 133 ) + a*JVS( 429 )
  a = -W( 100 ) / JVS( 448 )
  W( 100 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 449 )
  W( 131 ) = W( 131 ) + a*JVS( 450 )
  W( 132 ) = W( 132 ) + a*JVS( 451 )
  a = -W( 101 ) / JVS( 452 )
  W( 101 ) = -a
  W( 122 ) = W( 122 ) + a*JVS( 453 )
  W( 131 ) = W( 131 ) + a*JVS( 454 )
  W( 132 ) = W( 132 ) + a*JVS( 455 )
  a = -W( 103 ) / JVS( 461 )
  W( 103 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 462 )
  W( 131 ) = W( 131 ) + a*JVS( 463 )
  a = -W( 104 ) / JVS( 464 )
  W( 104 ) = -a
  W( 114 ) = W( 114 ) + a*JVS( 465 )
  W( 119 ) = W( 119 ) + a*JVS( 466 )
  W( 124 ) = W( 124 ) + a*JVS( 467 )
  W( 127 ) = W( 127 ) + a*JVS( 468 )
  W( 131 ) = W( 131 ) + a*JVS( 469 )
  a = -W( 105 ) / JVS( 472 )
  W( 105 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 473 )
  W( 136 ) = W( 136 ) + a*JVS( 474 )
  a = -W( 108 ) / JVS( 493 )
  W( 108 ) = -a
  W( 110 ) = W( 110 ) + a*JVS( 494 )
  W( 112 ) = W( 112 ) + a*JVS( 495 )
  W( 114 ) = W( 114 ) + a*JVS( 496 )
  W( 115 ) = W( 115 ) + a*JVS( 497 )
  W( 116 ) = W( 116 ) + a*JVS( 498 )
  W( 119 ) = W( 119 ) + a*JVS( 499 )
  W( 121 ) = W( 121 ) + a*JVS( 500 )
  W( 123 ) = W( 123 ) + a*JVS( 501 )
  W( 125 ) = W( 125 ) + a*JVS( 502 )
  W( 127 ) = W( 127 ) + a*JVS( 503 )
  W( 131 ) = W( 131 ) + a*JVS( 504 )
  W( 133 ) = W( 133 ) + a*JVS( 505 )
  W( 136 ) = W( 136 ) + a*JVS( 506 )
  a = -W( 110 ) / JVS( 526 )
  W( 110 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 527 )
  W( 131 ) = W( 131 ) + a*JVS( 528 )
  W( 133 ) = W( 133 ) + a*JVS( 529 )
  W( 136 ) = W( 136 ) + a*JVS( 530 )
  a = -W( 111 ) / JVS( 531 )
  W( 111 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 532 )
  W( 132 ) = W( 132 ) + a*JVS( 533 )
  W( 133 ) = W( 133 ) + a*JVS( 534 )
  W( 136 ) = W( 136 ) + a*JVS( 535 )
  a = -W( 112 ) / JVS( 536 )
  W( 112 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 537 )
  W( 131 ) = W( 131 ) + a*JVS( 538 )
  W( 132 ) = W( 132 ) + a*JVS( 539 )
  W( 133 ) = W( 133 ) + a*JVS( 540 )
  a = -W( 113 ) / JVS( 541 )
  W( 113 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 542 )
  W( 131 ) = W( 131 ) + a*JVS( 543 )
  W( 132 ) = W( 132 ) + a*JVS( 544 )
  W( 133 ) = W( 133 ) + a*JVS( 545 )
  a = -W( 114 ) / JVS( 546 )
  W( 114 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 547 )
  W( 131 ) = W( 131 ) + a*JVS( 548 )
  W( 136 ) = W( 136 ) + a*JVS( 549 )
  a = -W( 115 ) / JVS( 550 )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 136 ) = W( 136 ) + a*JVS( 553 )
  a = -W( 116 ) / JVS( 562 )
  W( 116 ) = -a
  W( 119 ) = W( 119 ) + a*JVS( 563 )
  W( 120 ) = W( 120 ) + a*JVS( 564 )
  W( 124 ) = W( 124 ) + a*JVS( 565 )
  W( 125 ) = W( 125 ) + a*JVS( 566 )
  W( 127 ) = W( 127 ) + a*JVS( 567 )
  W( 131 ) = W( 131 ) + a*JVS( 568 )
  W( 132 ) = W( 132 ) + a*JVS( 569 )
  W( 133 ) = W( 133 ) + a*JVS( 570 )
  W( 134 ) = W( 134 ) + a*JVS( 571 )
  W( 136 ) = W( 136 ) + a*JVS( 572 )
  a = -W( 117 ) / JVS( 580 )
  W( 117 ) = -a
  W( 119 ) = W( 119 ) + a*JVS( 581 )
  W( 121 ) = W( 121 ) + a*JVS( 582 )
  W( 125 ) = W( 125 ) + a*JVS( 583 )
  W( 126 ) = W( 126 ) + a*JVS( 584 )
  W( 127 ) = W( 127 ) + a*JVS( 585 )
  W( 128 ) = W( 128 ) + a*JVS( 586 )
  W( 131 ) = W( 131 ) + a*JVS( 587 )
  W( 132 ) = W( 132 ) + a*JVS( 588 )
  W( 133 ) = W( 133 ) + a*JVS( 589 )
  W( 135 ) = W( 135 ) + a*JVS( 590 )
  W( 136 ) = W( 136 ) + a*JVS( 591 )
  a = -W( 119 ) / JVS( 605 )
  W( 119 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 606 )
  W( 131 ) = W( 131 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 120 ) / JVS( 611 )
  W( 120 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 612 )
  W( 127 ) = W( 127 ) + a*JVS( 613 )
  W( 131 ) = W( 131 ) + a*JVS( 614 )
  W( 132 ) = W( 132 ) + a*JVS( 615 )
  W( 133 ) = W( 133 ) + a*JVS( 616 )
  W( 136 ) = W( 136 ) + a*JVS( 617 )
  a = -W( 121 ) / JVS( 623 )
  W( 121 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 624 )
  W( 127 ) = W( 127 ) + a*JVS( 625 )
  W( 131 ) = W( 131 ) + a*JVS( 626 )
  W( 132 ) = W( 132 ) + a*JVS( 627 )
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 135 ) = W( 135 ) + a*JVS( 630 )
  W( 136 ) = W( 136 ) + a*JVS( 631 )
  a = -W( 122 ) / JVS( 636 )
  W( 122 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 637 )
  W( 127 ) = W( 127 ) + a*JVS( 638 )
  W( 128 ) = W( 128 ) + a*JVS( 639 )
  W( 131 ) = W( 131 ) + a*JVS( 640 )
  W( 132 ) = W( 132 ) + a*JVS( 641 )
  W( 133 ) = W( 133 ) + a*JVS( 642 )
  W( 135 ) = W( 135 ) + a*JVS( 643 )
  W( 136 ) = W( 136 ) + a*JVS( 644 )
  a = -W( 123 ) / JVS( 656 )
  W( 123 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 657 )
  W( 126 ) = W( 126 ) + a*JVS( 658 )
  W( 127 ) = W( 127 ) + a*JVS( 659 )
  W( 128 ) = W( 128 ) + a*JVS( 660 )
  W( 131 ) = W( 131 ) + a*JVS( 661 )
  W( 132 ) = W( 132 ) + a*JVS( 662 )
  W( 133 ) = W( 133 ) + a*JVS( 663 )
  W( 135 ) = W( 135 ) + a*JVS( 664 )
  W( 136 ) = W( 136 ) + a*JVS( 665 )
  a = -W( 124 ) / JVS( 672 )
  W( 124 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 673 )
  W( 126 ) = W( 126 ) + a*JVS( 674 )
  W( 127 ) = W( 127 ) + a*JVS( 675 )
  W( 128 ) = W( 128 ) + a*JVS( 676 )
  W( 130 ) = W( 130 ) + a*JVS( 677 )
  W( 131 ) = W( 131 ) + a*JVS( 678 )
  W( 132 ) = W( 132 ) + a*JVS( 679 )
  W( 133 ) = W( 133 ) + a*JVS( 680 )
  W( 134 ) = W( 134 ) + a*JVS( 681 )
  W( 135 ) = W( 135 ) + a*JVS( 682 )
  W( 136 ) = W( 136 ) + a*JVS( 683 )
  a = -W( 125 ) / JVS( 687 )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 688 )
  W( 131 ) = W( 131 ) + a*JVS( 689 )
  W( 132 ) = W( 132 ) + a*JVS( 690 )
  W( 133 ) = W( 133 ) + a*JVS( 691 )
  W( 136 ) = W( 136 ) + a*JVS( 692 )
  a = -W( 126 ) / JVS( 701 )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 702 )
  W( 128 ) = W( 128 ) + a*JVS( 703 )
  W( 129 ) = W( 129 ) + a*JVS( 704 )
  W( 131 ) = W( 131 ) + a*JVS( 705 )
  W( 132 ) = W( 132 ) + a*JVS( 706 )
  W( 133 ) = W( 133 ) + a*JVS( 707 )
  W( 136 ) = W( 136 ) + a*JVS( 708 )
  a = -W( 127 ) / JVS( 718 )
  W( 127 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 719 )
  W( 130 ) = W( 130 ) + a*JVS( 720 )
  W( 131 ) = W( 131 ) + a*JVS( 721 )
  W( 132 ) = W( 132 ) + a*JVS( 722 )
  W( 133 ) = W( 133 ) + a*JVS( 723 )
  W( 136 ) = W( 136 ) + a*JVS( 724 )
  a = -W( 128 ) / JVS( 735 )
  W( 128 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 736 )
  W( 130 ) = W( 130 ) + a*JVS( 737 )
  W( 131 ) = W( 131 ) + a*JVS( 738 )
  W( 132 ) = W( 132 ) + a*JVS( 739 )
  W( 133 ) = W( 133 ) + a*JVS( 740 )
  W( 135 ) = W( 135 ) + a*JVS( 741 )
  W( 136 ) = W( 136 ) + a*JVS( 742 )
  a = -W( 129 ) / JVS( 768 )
  W( 129 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 769 )
  W( 131 ) = W( 131 ) + a*JVS( 770 )
  W( 132 ) = W( 132 ) + a*JVS( 771 )
  W( 133 ) = W( 133 ) + a*JVS( 772 )
  W( 134 ) = W( 134 ) + a*JVS( 773 )
  W( 135 ) = W( 135 ) + a*JVS( 774 )
  W( 136 ) = W( 136 ) + a*JVS( 775 )
  a = -W( 130 ) / JVS( 789 )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 790 )
  W( 132 ) = W( 132 ) + a*JVS( 791 )
  W( 133 ) = W( 133 ) + a*JVS( 792 )
  W( 134 ) = W( 134 ) + a*JVS( 793 )
  W( 135 ) = W( 135 ) + a*JVS( 794 )
  W( 136 ) = W( 136 ) + a*JVS( 795 )
  a = -W( 131 ) / JVS( 904 )
  W( 131 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 905 )
  W( 133 ) = W( 133 ) + a*JVS( 906 )
  W( 134 ) = W( 134 ) + a*JVS( 907 )
  W( 135 ) = W( 135 ) + a*JVS( 908 )
  W( 136 ) = W( 136 ) + a*JVS( 909 )
  JVS( 910) = W( 37 )
  JVS( 911) = W( 48 )
  JVS( 912) = W( 88 )
  JVS( 913) = W( 89 )
  JVS( 914) = W( 91 )
  JVS( 915) = W( 96 )
  JVS( 916) = W( 97 )
  JVS( 917) = W( 100 )
  JVS( 918) = W( 101 )
  JVS( 919) = W( 103 )
  JVS( 920) = W( 104 )
  JVS( 921) = W( 105 )
  JVS( 922) = W( 108 )
  JVS( 923) = W( 110 )
  JVS( 924) = W( 111 )
  JVS( 925) = W( 112 )
  JVS( 926) = W( 113 )
  JVS( 927) = W( 114 )
  JVS( 928) = W( 115 )
  JVS( 929) = W( 116 )
  JVS( 930) = W( 117 )
  JVS( 931) = W( 119 )
  JVS( 932) = W( 120 )
  JVS( 933) = W( 121 )
  JVS( 934) = W( 122 )
  JVS( 935) = W( 123 )
  JVS( 936) = W( 124 )
  JVS( 937) = W( 125 )
  JVS( 938) = W( 126 )
  JVS( 939) = W( 127 )
  JVS( 940) = W( 128 )
  JVS( 941) = W( 129 )
  JVS( 942) = W( 130 )
  JVS( 943) = W( 131 )
  JVS( 944) = W( 132 )
  JVS( 945) = W( 133 )
  JVS( 946) = W( 134 )
  JVS( 947) = W( 135 )
  JVS( 948) = W( 136 )
  IF ( ABS( JVS( 970 )) < TINY(a) ) THEN
         IER = 133
         RETURN
  END IF
   W( 97 ) = JVS( 949 )
   W( 102 ) = JVS( 950 )
   W( 106 ) = JVS( 951 )
   W( 111 ) = JVS( 952 )
   W( 112 ) = JVS( 953 )
   W( 113 ) = JVS( 954 )
   W( 115 ) = JVS( 955 )
   W( 117 ) = JVS( 956 )
   W( 119 ) = JVS( 957 )
   W( 120 ) = JVS( 958 )
   W( 121 ) = JVS( 959 )
   W( 122 ) = JVS( 960 )
   W( 124 ) = JVS( 961 )
   W( 125 ) = JVS( 962 )
   W( 126 ) = JVS( 963 )
   W( 127 ) = JVS( 964 )
   W( 128 ) = JVS( 965 )
   W( 129 ) = JVS( 966 )
   W( 130 ) = JVS( 967 )
   W( 131 ) = JVS( 968 )
   W( 132 ) = JVS( 969 )
   W( 133 ) = JVS( 970 )
   W( 134 ) = JVS( 971 )
   W( 135 ) = JVS( 972 )
   W( 136 ) = JVS( 973 )
  a = -W( 97 ) / JVS( 427 )
  W( 97 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 428 )
  W( 133 ) = W( 133 ) + a*JVS( 429 )
  a = -W( 102 ) / JVS( 456 )
  W( 102 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 457 )
  W( 131 ) = W( 131 ) + a*JVS( 458 )
  W( 132 ) = W( 132 ) + a*JVS( 459 )
  W( 133 ) = W( 133 ) + a*JVS( 460 )
  a = -W( 106 ) / JVS( 476 )
  W( 106 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 477 )
  W( 129 ) = W( 129 ) + a*JVS( 478 )
  W( 133 ) = W( 133 ) + a*JVS( 479 )
  W( 136 ) = W( 136 ) + a*JVS( 480 )
  a = -W( 111 ) / JVS( 531 )
  W( 111 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 532 )
  W( 132 ) = W( 132 ) + a*JVS( 533 )
  W( 133 ) = W( 133 ) + a*JVS( 534 )
  W( 136 ) = W( 136 ) + a*JVS( 535 )
  a = -W( 112 ) / JVS( 536 )
  W( 112 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 537 )
  W( 131 ) = W( 131 ) + a*JVS( 538 )
  W( 132 ) = W( 132 ) + a*JVS( 539 )
  W( 133 ) = W( 133 ) + a*JVS( 540 )
  a = -W( 113 ) / JVS( 541 )
  W( 113 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 542 )
  W( 131 ) = W( 131 ) + a*JVS( 543 )
  W( 132 ) = W( 132 ) + a*JVS( 544 )
  W( 133 ) = W( 133 ) + a*JVS( 545 )
  a = -W( 115 ) / JVS( 550 )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 136 ) = W( 136 ) + a*JVS( 553 )
  a = -W( 117 ) / JVS( 580 )
  W( 117 ) = -a
  W( 119 ) = W( 119 ) + a*JVS( 581 )
  W( 121 ) = W( 121 ) + a*JVS( 582 )
  W( 125 ) = W( 125 ) + a*JVS( 583 )
  W( 126 ) = W( 126 ) + a*JVS( 584 )
  W( 127 ) = W( 127 ) + a*JVS( 585 )
  W( 128 ) = W( 128 ) + a*JVS( 586 )
  W( 131 ) = W( 131 ) + a*JVS( 587 )
  W( 132 ) = W( 132 ) + a*JVS( 588 )
  W( 133 ) = W( 133 ) + a*JVS( 589 )
  W( 135 ) = W( 135 ) + a*JVS( 590 )
  W( 136 ) = W( 136 ) + a*JVS( 591 )
  a = -W( 119 ) / JVS( 605 )
  W( 119 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 606 )
  W( 131 ) = W( 131 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 120 ) / JVS( 611 )
  W( 120 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 612 )
  W( 127 ) = W( 127 ) + a*JVS( 613 )
  W( 131 ) = W( 131 ) + a*JVS( 614 )
  W( 132 ) = W( 132 ) + a*JVS( 615 )
  W( 133 ) = W( 133 ) + a*JVS( 616 )
  W( 136 ) = W( 136 ) + a*JVS( 617 )
  a = -W( 121 ) / JVS( 623 )
  W( 121 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 624 )
  W( 127 ) = W( 127 ) + a*JVS( 625 )
  W( 131 ) = W( 131 ) + a*JVS( 626 )
  W( 132 ) = W( 132 ) + a*JVS( 627 )
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 135 ) = W( 135 ) + a*JVS( 630 )
  W( 136 ) = W( 136 ) + a*JVS( 631 )
  a = -W( 122 ) / JVS( 636 )
  W( 122 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 637 )
  W( 127 ) = W( 127 ) + a*JVS( 638 )
  W( 128 ) = W( 128 ) + a*JVS( 639 )
  W( 131 ) = W( 131 ) + a*JVS( 640 )
  W( 132 ) = W( 132 ) + a*JVS( 641 )
  W( 133 ) = W( 133 ) + a*JVS( 642 )
  W( 135 ) = W( 135 ) + a*JVS( 643 )
  W( 136 ) = W( 136 ) + a*JVS( 644 )
  a = -W( 124 ) / JVS( 672 )
  W( 124 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 673 )
  W( 126 ) = W( 126 ) + a*JVS( 674 )
  W( 127 ) = W( 127 ) + a*JVS( 675 )
  W( 128 ) = W( 128 ) + a*JVS( 676 )
  W( 130 ) = W( 130 ) + a*JVS( 677 )
  W( 131 ) = W( 131 ) + a*JVS( 678 )
  W( 132 ) = W( 132 ) + a*JVS( 679 )
  W( 133 ) = W( 133 ) + a*JVS( 680 )
  W( 134 ) = W( 134 ) + a*JVS( 681 )
  W( 135 ) = W( 135 ) + a*JVS( 682 )
  W( 136 ) = W( 136 ) + a*JVS( 683 )
  a = -W( 125 ) / JVS( 687 )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 688 )
  W( 131 ) = W( 131 ) + a*JVS( 689 )
  W( 132 ) = W( 132 ) + a*JVS( 690 )
  W( 133 ) = W( 133 ) + a*JVS( 691 )
  W( 136 ) = W( 136 ) + a*JVS( 692 )
  a = -W( 126 ) / JVS( 701 )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 702 )
  W( 128 ) = W( 128 ) + a*JVS( 703 )
  W( 129 ) = W( 129 ) + a*JVS( 704 )
  W( 131 ) = W( 131 ) + a*JVS( 705 )
  W( 132 ) = W( 132 ) + a*JVS( 706 )
  W( 133 ) = W( 133 ) + a*JVS( 707 )
  W( 136 ) = W( 136 ) + a*JVS( 708 )
  a = -W( 127 ) / JVS( 718 )
  W( 127 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 719 )
  W( 130 ) = W( 130 ) + a*JVS( 720 )
  W( 131 ) = W( 131 ) + a*JVS( 721 )
  W( 132 ) = W( 132 ) + a*JVS( 722 )
  W( 133 ) = W( 133 ) + a*JVS( 723 )
  W( 136 ) = W( 136 ) + a*JVS( 724 )
  a = -W( 128 ) / JVS( 735 )
  W( 128 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 736 )
  W( 130 ) = W( 130 ) + a*JVS( 737 )
  W( 131 ) = W( 131 ) + a*JVS( 738 )
  W( 132 ) = W( 132 ) + a*JVS( 739 )
  W( 133 ) = W( 133 ) + a*JVS( 740 )
  W( 135 ) = W( 135 ) + a*JVS( 741 )
  W( 136 ) = W( 136 ) + a*JVS( 742 )
  a = -W( 129 ) / JVS( 768 )
  W( 129 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 769 )
  W( 131 ) = W( 131 ) + a*JVS( 770 )
  W( 132 ) = W( 132 ) + a*JVS( 771 )
  W( 133 ) = W( 133 ) + a*JVS( 772 )
  W( 134 ) = W( 134 ) + a*JVS( 773 )
  W( 135 ) = W( 135 ) + a*JVS( 774 )
  W( 136 ) = W( 136 ) + a*JVS( 775 )
  a = -W( 130 ) / JVS( 789 )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 790 )
  W( 132 ) = W( 132 ) + a*JVS( 791 )
  W( 133 ) = W( 133 ) + a*JVS( 792 )
  W( 134 ) = W( 134 ) + a*JVS( 793 )
  W( 135 ) = W( 135 ) + a*JVS( 794 )
  W( 136 ) = W( 136 ) + a*JVS( 795 )
  a = -W( 131 ) / JVS( 904 )
  W( 131 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 905 )
  W( 133 ) = W( 133 ) + a*JVS( 906 )
  W( 134 ) = W( 134 ) + a*JVS( 907 )
  W( 135 ) = W( 135 ) + a*JVS( 908 )
  W( 136 ) = W( 136 ) + a*JVS( 909 )
  a = -W( 132 ) / JVS( 944 )
  W( 132 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 945 )
  W( 134 ) = W( 134 ) + a*JVS( 946 )
  W( 135 ) = W( 135 ) + a*JVS( 947 )
  W( 136 ) = W( 136 ) + a*JVS( 948 )
  JVS( 949) = W( 97 )
  JVS( 950) = W( 102 )
  JVS( 951) = W( 106 )
  JVS( 952) = W( 111 )
  JVS( 953) = W( 112 )
  JVS( 954) = W( 113 )
  JVS( 955) = W( 115 )
  JVS( 956) = W( 117 )
  JVS( 957) = W( 119 )
  JVS( 958) = W( 120 )
  JVS( 959) = W( 121 )
  JVS( 960) = W( 122 )
  JVS( 961) = W( 124 )
  JVS( 962) = W( 125 )
  JVS( 963) = W( 126 )
  JVS( 964) = W( 127 )
  JVS( 965) = W( 128 )
  JVS( 966) = W( 129 )
  JVS( 967) = W( 130 )
  JVS( 968) = W( 131 )
  JVS( 969) = W( 132 )
  JVS( 970) = W( 133 )
  JVS( 971) = W( 134 )
  JVS( 972) = W( 135 )
  JVS( 973) = W( 136 )
  IF ( ABS( JVS( 985 )) < TINY(a) ) THEN
         IER = 134
         RETURN
  END IF
   W( 118 ) = JVS( 974 )
   W( 119 ) = JVS( 975 )
   W( 125 ) = JVS( 976 )
   W( 126 ) = JVS( 977 )
   W( 127 ) = JVS( 978 )
   W( 128 ) = JVS( 979 )
   W( 129 ) = JVS( 980 )
   W( 130 ) = JVS( 981 )
   W( 131 ) = JVS( 982 )
   W( 132 ) = JVS( 983 )
   W( 133 ) = JVS( 984 )
   W( 134 ) = JVS( 985 )
   W( 135 ) = JVS( 986 )
   W( 136 ) = JVS( 987 )
  a = -W( 118 ) / JVS( 593 )
  W( 118 ) = -a
  W( 119 ) = W( 119 ) + a*JVS( 594 )
  W( 125 ) = W( 125 ) + a*JVS( 595 )
  W( 126 ) = W( 126 ) + a*JVS( 596 )
  W( 127 ) = W( 127 ) + a*JVS( 597 )
  W( 128 ) = W( 128 ) + a*JVS( 598 )
  W( 131 ) = W( 131 ) + a*JVS( 599 )
  W( 132 ) = W( 132 ) + a*JVS( 600 )
  W( 133 ) = W( 133 ) + a*JVS( 601 )
  W( 134 ) = W( 134 ) + a*JVS( 602 )
  W( 135 ) = W( 135 ) + a*JVS( 603 )
  W( 136 ) = W( 136 ) + a*JVS( 604 )
  a = -W( 119 ) / JVS( 605 )
  W( 119 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 606 )
  W( 131 ) = W( 131 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 125 ) / JVS( 687 )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 688 )
  W( 131 ) = W( 131 ) + a*JVS( 689 )
  W( 132 ) = W( 132 ) + a*JVS( 690 )
  W( 133 ) = W( 133 ) + a*JVS( 691 )
  W( 136 ) = W( 136 ) + a*JVS( 692 )
  a = -W( 126 ) / JVS( 701 )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 702 )
  W( 128 ) = W( 128 ) + a*JVS( 703 )
  W( 129 ) = W( 129 ) + a*JVS( 704 )
  W( 131 ) = W( 131 ) + a*JVS( 705 )
  W( 132 ) = W( 132 ) + a*JVS( 706 )
  W( 133 ) = W( 133 ) + a*JVS( 707 )
  W( 136 ) = W( 136 ) + a*JVS( 708 )
  a = -W( 127 ) / JVS( 718 )
  W( 127 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 719 )
  W( 130 ) = W( 130 ) + a*JVS( 720 )
  W( 131 ) = W( 131 ) + a*JVS( 721 )
  W( 132 ) = W( 132 ) + a*JVS( 722 )
  W( 133 ) = W( 133 ) + a*JVS( 723 )
  W( 136 ) = W( 136 ) + a*JVS( 724 )
  a = -W( 128 ) / JVS( 735 )
  W( 128 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 736 )
  W( 130 ) = W( 130 ) + a*JVS( 737 )
  W( 131 ) = W( 131 ) + a*JVS( 738 )
  W( 132 ) = W( 132 ) + a*JVS( 739 )
  W( 133 ) = W( 133 ) + a*JVS( 740 )
  W( 135 ) = W( 135 ) + a*JVS( 741 )
  W( 136 ) = W( 136 ) + a*JVS( 742 )
  a = -W( 129 ) / JVS( 768 )
  W( 129 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 769 )
  W( 131 ) = W( 131 ) + a*JVS( 770 )
  W( 132 ) = W( 132 ) + a*JVS( 771 )
  W( 133 ) = W( 133 ) + a*JVS( 772 )
  W( 134 ) = W( 134 ) + a*JVS( 773 )
  W( 135 ) = W( 135 ) + a*JVS( 774 )
  W( 136 ) = W( 136 ) + a*JVS( 775 )
  a = -W( 130 ) / JVS( 789 )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 790 )
  W( 132 ) = W( 132 ) + a*JVS( 791 )
  W( 133 ) = W( 133 ) + a*JVS( 792 )
  W( 134 ) = W( 134 ) + a*JVS( 793 )
  W( 135 ) = W( 135 ) + a*JVS( 794 )
  W( 136 ) = W( 136 ) + a*JVS( 795 )
  a = -W( 131 ) / JVS( 904 )
  W( 131 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 905 )
  W( 133 ) = W( 133 ) + a*JVS( 906 )
  W( 134 ) = W( 134 ) + a*JVS( 907 )
  W( 135 ) = W( 135 ) + a*JVS( 908 )
  W( 136 ) = W( 136 ) + a*JVS( 909 )
  a = -W( 132 ) / JVS( 944 )
  W( 132 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 945 )
  W( 134 ) = W( 134 ) + a*JVS( 946 )
  W( 135 ) = W( 135 ) + a*JVS( 947 )
  W( 136 ) = W( 136 ) + a*JVS( 948 )
  a = -W( 133 ) / JVS( 970 )
  W( 133 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 971 )
  W( 135 ) = W( 135 ) + a*JVS( 972 )
  W( 136 ) = W( 136 ) + a*JVS( 973 )
  JVS( 974) = W( 118 )
  JVS( 975) = W( 119 )
  JVS( 976) = W( 125 )
  JVS( 977) = W( 126 )
  JVS( 978) = W( 127 )
  JVS( 979) = W( 128 )
  JVS( 980) = W( 129 )
  JVS( 981) = W( 130 )
  JVS( 982) = W( 131 )
  JVS( 983) = W( 132 )
  JVS( 984) = W( 133 )
  JVS( 985) = W( 134 )
  JVS( 986) = W( 135 )
  JVS( 987) = W( 136 )
  IF ( ABS( JVS( 1000 )) < TINY(a) ) THEN
         IER = 135
         RETURN
  END IF
   W( 112 ) = JVS( 988 )
   W( 113 ) = JVS( 989 )
   W( 115 ) = JVS( 990 )
   W( 125 ) = JVS( 991 )
   W( 127 ) = JVS( 992 )
   W( 128 ) = JVS( 993 )
   W( 129 ) = JVS( 994 )
   W( 130 ) = JVS( 995 )
   W( 131 ) = JVS( 996 )
   W( 132 ) = JVS( 997 )
   W( 133 ) = JVS( 998 )
   W( 134 ) = JVS( 999 )
   W( 135 ) = JVS( 1000 )
   W( 136 ) = JVS( 1001 )
  a = -W( 112 ) / JVS( 536 )
  W( 112 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 537 )
  W( 131 ) = W( 131 ) + a*JVS( 538 )
  W( 132 ) = W( 132 ) + a*JVS( 539 )
  W( 133 ) = W( 133 ) + a*JVS( 540 )
  a = -W( 113 ) / JVS( 541 )
  W( 113 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 542 )
  W( 131 ) = W( 131 ) + a*JVS( 543 )
  W( 132 ) = W( 132 ) + a*JVS( 544 )
  W( 133 ) = W( 133 ) + a*JVS( 545 )
  a = -W( 115 ) / JVS( 550 )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 136 ) = W( 136 ) + a*JVS( 553 )
  a = -W( 125 ) / JVS( 687 )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 688 )
  W( 131 ) = W( 131 ) + a*JVS( 689 )
  W( 132 ) = W( 132 ) + a*JVS( 690 )
  W( 133 ) = W( 133 ) + a*JVS( 691 )
  W( 136 ) = W( 136 ) + a*JVS( 692 )
  a = -W( 127 ) / JVS( 718 )
  W( 127 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 719 )
  W( 130 ) = W( 130 ) + a*JVS( 720 )
  W( 131 ) = W( 131 ) + a*JVS( 721 )
  W( 132 ) = W( 132 ) + a*JVS( 722 )
  W( 133 ) = W( 133 ) + a*JVS( 723 )
  W( 136 ) = W( 136 ) + a*JVS( 724 )
  a = -W( 128 ) / JVS( 735 )
  W( 128 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 736 )
  W( 130 ) = W( 130 ) + a*JVS( 737 )
  W( 131 ) = W( 131 ) + a*JVS( 738 )
  W( 132 ) = W( 132 ) + a*JVS( 739 )
  W( 133 ) = W( 133 ) + a*JVS( 740 )
  W( 135 ) = W( 135 ) + a*JVS( 741 )
  W( 136 ) = W( 136 ) + a*JVS( 742 )
  a = -W( 129 ) / JVS( 768 )
  W( 129 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 769 )
  W( 131 ) = W( 131 ) + a*JVS( 770 )
  W( 132 ) = W( 132 ) + a*JVS( 771 )
  W( 133 ) = W( 133 ) + a*JVS( 772 )
  W( 134 ) = W( 134 ) + a*JVS( 773 )
  W( 135 ) = W( 135 ) + a*JVS( 774 )
  W( 136 ) = W( 136 ) + a*JVS( 775 )
  a = -W( 130 ) / JVS( 789 )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 790 )
  W( 132 ) = W( 132 ) + a*JVS( 791 )
  W( 133 ) = W( 133 ) + a*JVS( 792 )
  W( 134 ) = W( 134 ) + a*JVS( 793 )
  W( 135 ) = W( 135 ) + a*JVS( 794 )
  W( 136 ) = W( 136 ) + a*JVS( 795 )
  a = -W( 131 ) / JVS( 904 )
  W( 131 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 905 )
  W( 133 ) = W( 133 ) + a*JVS( 906 )
  W( 134 ) = W( 134 ) + a*JVS( 907 )
  W( 135 ) = W( 135 ) + a*JVS( 908 )
  W( 136 ) = W( 136 ) + a*JVS( 909 )
  a = -W( 132 ) / JVS( 944 )
  W( 132 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 945 )
  W( 134 ) = W( 134 ) + a*JVS( 946 )
  W( 135 ) = W( 135 ) + a*JVS( 947 )
  W( 136 ) = W( 136 ) + a*JVS( 948 )
  a = -W( 133 ) / JVS( 970 )
  W( 133 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 971 )
  W( 135 ) = W( 135 ) + a*JVS( 972 )
  W( 136 ) = W( 136 ) + a*JVS( 973 )
  a = -W( 134 ) / JVS( 985 )
  W( 134 ) = -a
  W( 135 ) = W( 135 ) + a*JVS( 986 )
  W( 136 ) = W( 136 ) + a*JVS( 987 )
  JVS( 988) = W( 112 )
  JVS( 989) = W( 113 )
  JVS( 990) = W( 115 )
  JVS( 991) = W( 125 )
  JVS( 992) = W( 127 )
  JVS( 993) = W( 128 )
  JVS( 994) = W( 129 )
  JVS( 995) = W( 130 )
  JVS( 996) = W( 131 )
  JVS( 997) = W( 132 )
  JVS( 998) = W( 133 )
  JVS( 999) = W( 134 )
  JVS( 1000) = W( 135 )
  JVS( 1001) = W( 136 )
  IF ( ABS( JVS( 1029 )) < TINY(a) ) THEN
         IER = 136
         RETURN
  END IF
   W( 90 ) = JVS( 1002 )
   W( 94 ) = JVS( 1003 )
   W( 95 ) = JVS( 1004 )
   W( 105 ) = JVS( 1005 )
   W( 106 ) = JVS( 1006 )
   W( 107 ) = JVS( 1007 )
   W( 114 ) = JVS( 1008 )
   W( 115 ) = JVS( 1009 )
   W( 116 ) = JVS( 1010 )
   W( 117 ) = JVS( 1011 )
   W( 119 ) = JVS( 1012 )
   W( 120 ) = JVS( 1013 )
   W( 121 ) = JVS( 1014 )
   W( 122 ) = JVS( 1015 )
   W( 123 ) = JVS( 1016 )
   W( 124 ) = JVS( 1017 )
   W( 125 ) = JVS( 1018 )
   W( 126 ) = JVS( 1019 )
   W( 127 ) = JVS( 1020 )
   W( 128 ) = JVS( 1021 )
   W( 129 ) = JVS( 1022 )
   W( 130 ) = JVS( 1023 )
   W( 131 ) = JVS( 1024 )
   W( 132 ) = JVS( 1025 )
   W( 133 ) = JVS( 1026 )
   W( 134 ) = JVS( 1027 )
   W( 135 ) = JVS( 1028 )
   W( 136 ) = JVS( 1029 )
  a = -W( 90 ) / JVS( 398 )
  W( 90 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 399 )
  W( 136 ) = W( 136 ) + a*JVS( 400 )
  a = -W( 94 ) / JVS( 413 )
  W( 94 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 414 )
  W( 131 ) = W( 131 ) + a*JVS( 415 )
  W( 136 ) = W( 136 ) + a*JVS( 416 )
  a = -W( 95 ) / JVS( 417 )
  W( 95 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 418 )
  W( 131 ) = W( 131 ) + a*JVS( 419 )
  W( 136 ) = W( 136 ) + a*JVS( 420 )
  a = -W( 105 ) / JVS( 472 )
  W( 105 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 473 )
  W( 136 ) = W( 136 ) + a*JVS( 474 )
  a = -W( 106 ) / JVS( 476 )
  W( 106 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 477 )
  W( 129 ) = W( 129 ) + a*JVS( 478 )
  W( 133 ) = W( 133 ) + a*JVS( 479 )
  W( 136 ) = W( 136 ) + a*JVS( 480 )
  a = -W( 107 ) / JVS( 483 )
  W( 107 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 484 )
  W( 121 ) = W( 121 ) + a*JVS( 485 )
  W( 123 ) = W( 123 ) + a*JVS( 486 )
  W( 125 ) = W( 125 ) + a*JVS( 487 )
  W( 129 ) = W( 129 ) + a*JVS( 488 )
  W( 131 ) = W( 131 ) + a*JVS( 489 )
  W( 132 ) = W( 132 ) + a*JVS( 490 )
  W( 136 ) = W( 136 ) + a*JVS( 491 )
  a = -W( 114 ) / JVS( 546 )
  W( 114 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 547 )
  W( 131 ) = W( 131 ) + a*JVS( 548 )
  W( 136 ) = W( 136 ) + a*JVS( 549 )
  a = -W( 115 ) / JVS( 550 )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 136 ) = W( 136 ) + a*JVS( 553 )
  a = -W( 116 ) / JVS( 562 )
  W( 116 ) = -a
  W( 119 ) = W( 119 ) + a*JVS( 563 )
  W( 120 ) = W( 120 ) + a*JVS( 564 )
  W( 124 ) = W( 124 ) + a*JVS( 565 )
  W( 125 ) = W( 125 ) + a*JVS( 566 )
  W( 127 ) = W( 127 ) + a*JVS( 567 )
  W( 131 ) = W( 131 ) + a*JVS( 568 )
  W( 132 ) = W( 132 ) + a*JVS( 569 )
  W( 133 ) = W( 133 ) + a*JVS( 570 )
  W( 134 ) = W( 134 ) + a*JVS( 571 )
  W( 136 ) = W( 136 ) + a*JVS( 572 )
  a = -W( 117 ) / JVS( 580 )
  W( 117 ) = -a
  W( 119 ) = W( 119 ) + a*JVS( 581 )
  W( 121 ) = W( 121 ) + a*JVS( 582 )
  W( 125 ) = W( 125 ) + a*JVS( 583 )
  W( 126 ) = W( 126 ) + a*JVS( 584 )
  W( 127 ) = W( 127 ) + a*JVS( 585 )
  W( 128 ) = W( 128 ) + a*JVS( 586 )
  W( 131 ) = W( 131 ) + a*JVS( 587 )
  W( 132 ) = W( 132 ) + a*JVS( 588 )
  W( 133 ) = W( 133 ) + a*JVS( 589 )
  W( 135 ) = W( 135 ) + a*JVS( 590 )
  W( 136 ) = W( 136 ) + a*JVS( 591 )
  a = -W( 119 ) / JVS( 605 )
  W( 119 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 606 )
  W( 131 ) = W( 131 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 120 ) / JVS( 611 )
  W( 120 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 612 )
  W( 127 ) = W( 127 ) + a*JVS( 613 )
  W( 131 ) = W( 131 ) + a*JVS( 614 )
  W( 132 ) = W( 132 ) + a*JVS( 615 )
  W( 133 ) = W( 133 ) + a*JVS( 616 )
  W( 136 ) = W( 136 ) + a*JVS( 617 )
  a = -W( 121 ) / JVS( 623 )
  W( 121 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 624 )
  W( 127 ) = W( 127 ) + a*JVS( 625 )
  W( 131 ) = W( 131 ) + a*JVS( 626 )
  W( 132 ) = W( 132 ) + a*JVS( 627 )
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 135 ) = W( 135 ) + a*JVS( 630 )
  W( 136 ) = W( 136 ) + a*JVS( 631 )
  a = -W( 122 ) / JVS( 636 )
  W( 122 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 637 )
  W( 127 ) = W( 127 ) + a*JVS( 638 )
  W( 128 ) = W( 128 ) + a*JVS( 639 )
  W( 131 ) = W( 131 ) + a*JVS( 640 )
  W( 132 ) = W( 132 ) + a*JVS( 641 )
  W( 133 ) = W( 133 ) + a*JVS( 642 )
  W( 135 ) = W( 135 ) + a*JVS( 643 )
  W( 136 ) = W( 136 ) + a*JVS( 644 )
  a = -W( 123 ) / JVS( 656 )
  W( 123 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 657 )
  W( 126 ) = W( 126 ) + a*JVS( 658 )
  W( 127 ) = W( 127 ) + a*JVS( 659 )
  W( 128 ) = W( 128 ) + a*JVS( 660 )
  W( 131 ) = W( 131 ) + a*JVS( 661 )
  W( 132 ) = W( 132 ) + a*JVS( 662 )
  W( 133 ) = W( 133 ) + a*JVS( 663 )
  W( 135 ) = W( 135 ) + a*JVS( 664 )
  W( 136 ) = W( 136 ) + a*JVS( 665 )
  a = -W( 124 ) / JVS( 672 )
  W( 124 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 673 )
  W( 126 ) = W( 126 ) + a*JVS( 674 )
  W( 127 ) = W( 127 ) + a*JVS( 675 )
  W( 128 ) = W( 128 ) + a*JVS( 676 )
  W( 130 ) = W( 130 ) + a*JVS( 677 )
  W( 131 ) = W( 131 ) + a*JVS( 678 )
  W( 132 ) = W( 132 ) + a*JVS( 679 )
  W( 133 ) = W( 133 ) + a*JVS( 680 )
  W( 134 ) = W( 134 ) + a*JVS( 681 )
  W( 135 ) = W( 135 ) + a*JVS( 682 )
  W( 136 ) = W( 136 ) + a*JVS( 683 )
  a = -W( 125 ) / JVS( 687 )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 688 )
  W( 131 ) = W( 131 ) + a*JVS( 689 )
  W( 132 ) = W( 132 ) + a*JVS( 690 )
  W( 133 ) = W( 133 ) + a*JVS( 691 )
  W( 136 ) = W( 136 ) + a*JVS( 692 )
  a = -W( 126 ) / JVS( 701 )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 702 )
  W( 128 ) = W( 128 ) + a*JVS( 703 )
  W( 129 ) = W( 129 ) + a*JVS( 704 )
  W( 131 ) = W( 131 ) + a*JVS( 705 )
  W( 132 ) = W( 132 ) + a*JVS( 706 )
  W( 133 ) = W( 133 ) + a*JVS( 707 )
  W( 136 ) = W( 136 ) + a*JVS( 708 )
  a = -W( 127 ) / JVS( 718 )
  W( 127 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 719 )
  W( 130 ) = W( 130 ) + a*JVS( 720 )
  W( 131 ) = W( 131 ) + a*JVS( 721 )
  W( 132 ) = W( 132 ) + a*JVS( 722 )
  W( 133 ) = W( 133 ) + a*JVS( 723 )
  W( 136 ) = W( 136 ) + a*JVS( 724 )
  a = -W( 128 ) / JVS( 735 )
  W( 128 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 736 )
  W( 130 ) = W( 130 ) + a*JVS( 737 )
  W( 131 ) = W( 131 ) + a*JVS( 738 )
  W( 132 ) = W( 132 ) + a*JVS( 739 )
  W( 133 ) = W( 133 ) + a*JVS( 740 )
  W( 135 ) = W( 135 ) + a*JVS( 741 )
  W( 136 ) = W( 136 ) + a*JVS( 742 )
  a = -W( 129 ) / JVS( 768 )
  W( 129 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 769 )
  W( 131 ) = W( 131 ) + a*JVS( 770 )
  W( 132 ) = W( 132 ) + a*JVS( 771 )
  W( 133 ) = W( 133 ) + a*JVS( 772 )
  W( 134 ) = W( 134 ) + a*JVS( 773 )
  W( 135 ) = W( 135 ) + a*JVS( 774 )
  W( 136 ) = W( 136 ) + a*JVS( 775 )
  a = -W( 130 ) / JVS( 789 )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 790 )
  W( 132 ) = W( 132 ) + a*JVS( 791 )
  W( 133 ) = W( 133 ) + a*JVS( 792 )
  W( 134 ) = W( 134 ) + a*JVS( 793 )
  W( 135 ) = W( 135 ) + a*JVS( 794 )
  W( 136 ) = W( 136 ) + a*JVS( 795 )
  a = -W( 131 ) / JVS( 904 )
  W( 131 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 905 )
  W( 133 ) = W( 133 ) + a*JVS( 906 )
  W( 134 ) = W( 134 ) + a*JVS( 907 )
  W( 135 ) = W( 135 ) + a*JVS( 908 )
  W( 136 ) = W( 136 ) + a*JVS( 909 )
  a = -W( 132 ) / JVS( 944 )
  W( 132 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 945 )
  W( 134 ) = W( 134 ) + a*JVS( 946 )
  W( 135 ) = W( 135 ) + a*JVS( 947 )
  W( 136 ) = W( 136 ) + a*JVS( 948 )
  a = -W( 133 ) / JVS( 970 )
  W( 133 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 971 )
  W( 135 ) = W( 135 ) + a*JVS( 972 )
  W( 136 ) = W( 136 ) + a*JVS( 973 )
  a = -W( 134 ) / JVS( 985 )
  W( 134 ) = -a
  W( 135 ) = W( 135 ) + a*JVS( 986 )
  W( 136 ) = W( 136 ) + a*JVS( 987 )
  a = -W( 135 ) / JVS( 1000 )
  W( 135 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 1001 )
  JVS( 1002) = W( 90 )
  JVS( 1003) = W( 94 )
  JVS( 1004) = W( 95 )
  JVS( 1005) = W( 105 )
  JVS( 1006) = W( 106 )
  JVS( 1007) = W( 107 )
  JVS( 1008) = W( 114 )
  JVS( 1009) = W( 115 )
  JVS( 1010) = W( 116 )
  JVS( 1011) = W( 117 )
  JVS( 1012) = W( 119 )
  JVS( 1013) = W( 120 )
  JVS( 1014) = W( 121 )
  JVS( 1015) = W( 122 )
  JVS( 1016) = W( 123 )
  JVS( 1017) = W( 124 )
  JVS( 1018) = W( 125 )
  JVS( 1019) = W( 126 )
  JVS( 1020) = W( 127 )
  JVS( 1021) = W( 128 )
  JVS( 1022) = W( 129 )
  JVS( 1023) = W( 130 )
  JVS( 1024) = W( 131 )
  JVS( 1025) = W( 132 )
  JVS( 1026) = W( 133 )
  JVS( 1027) = W( 134 )
  JVS( 1028) = W( 135 )
  JVS( 1029) = W( 136 )
   END SUBROUTINE decomp_cbmz_mosaic_8bin_vbs9
END MODULE cbmz_mosaic_8bin_vbs9_Integrator
