MODULE saprc99_mosaic_4bin_vbs9_Integrator
 USE saprc99_mosaic_4bin_vbs9_Parameters
 USE saprc99_mosaic_4bin_vbs9_Precision
 USE saprc99_mosaic_4bin_vbs9_JacobianSP
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
SUBROUTINE saprc99_mosaic_4bin_vbs9_INTEGRATE( TIN, TOUT, &
  FIX, VAR, RCONST, ATOL, RTOL, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )
   USE saprc99_mosaic_4bin_vbs9_Parameters
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
   CALL saprc99_mosaic_4bin_vbs9_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT, &
         ATOL,RTOL, &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
   STEPMIN = RCNTRL(ihexit)
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U)) IERR_U = IERR
END SUBROUTINE saprc99_mosaic_4bin_vbs9_INTEGRATE
SUBROUTINE saprc99_mosaic_4bin_vbs9_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol, &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
  USE saprc99_mosaic_4bin_vbs9_Parameters
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
      CALL saprc99_mosaic_4bin_vbs9_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL saprc99_mosaic_4bin_vbs9_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF
   Roundoff = saprc99_mosaic_4bin_vbs9_WLAMCH('E')
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL saprc99_mosaic_4bin_vbs9_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL saprc99_mosaic_4bin_vbs9_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL saprc99_mosaic_4bin_vbs9_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL saprc99_mosaic_4bin_vbs9_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL saprc99_mosaic_4bin_vbs9_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL saprc99_mosaic_4bin_vbs9_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL saprc99_mosaic_4bin_vbs9_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL saprc99_mosaic_4bin_vbs9_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO
   SELECT CASE (Method)
     CASE (1)
       CALL saprc99_mosaic_4bin_vbs9_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL saprc99_mosaic_4bin_vbs9_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL saprc99_mosaic_4bin_vbs9_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL saprc99_mosaic_4bin_vbs9_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL saprc99_mosaic_4bin_vbs9_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL saprc99_mosaic_4bin_vbs9_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT
   CALL saprc99_mosaic_4bin_vbs9_ros_Integrator(Y,Tstart,Tend,Texit, &
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
 SUBROUTINE saprc99_mosaic_4bin_vbs9_ros_ErrorMsg(Code,T,H,IERR)
   USE saprc99_mosaic_4bin_vbs9_Precision
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
 END SUBROUTINE saprc99_mosaic_4bin_vbs9_ros_ErrorMsg
 SUBROUTINE saprc99_mosaic_4bin_vbs9_ros_Integrator (Y, Tstart, Tend, T, &
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
      CALL saprc99_mosaic_4bin_vbs9_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN
      CALL saprc99_mosaic_4bin_vbs9_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF
   Hexit = H
   H = MIN(H,ABS(Tend-T))
   CALL saprc99_mosaic_4bin_vbs9_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF (.NOT.Autonomous) THEN
      CALL saprc99_mosaic_4bin_vbs9_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF
   CALL saprc99_mosaic_4bin_vbs9_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)
UntilAccepted: DO
   CALL saprc99_mosaic_4bin_vbs9_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec, Nsng )
   IF (Singular) THEN
       CALL saprc99_mosaic_4bin_vbs9_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF
Stage: DO istage = 1, ros_S
       ioffset = NVAR*(istage-1)
       IF ( istage == 1 ) THEN
         CALL saprc99_mosaic_4bin_vbs9_WCOPY(NVAR,Fcn0,1,Fcn,1)
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL saprc99_mosaic_4bin_vbs9_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL saprc99_mosaic_4bin_vbs9_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL saprc99_mosaic_4bin_vbs9_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF
       CALL saprc99_mosaic_4bin_vbs9_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL saprc99_mosaic_4bin_vbs9_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL saprc99_mosaic_4bin_vbs9_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL saprc99_mosaic_4bin_vbs9_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)
   END DO Stage
   CALL saprc99_mosaic_4bin_vbs9_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL saprc99_mosaic_4bin_vbs9_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO
   CALL saprc99_mosaic_4bin_vbs9_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL saprc99_mosaic_4bin_vbs9_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = saprc99_mosaic_4bin_vbs9_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )
   Fac = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac
   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN
      Nacc = Nacc+1
      CALL saprc99_mosaic_4bin_vbs9_WCOPY(NVAR,Ynew,1,Y,1)
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
  END SUBROUTINE saprc99_mosaic_4bin_vbs9_ros_Integrator
  REAL(kind=dp) FUNCTION saprc99_mosaic_4bin_vbs9_ros_ErrorNorm ( Y, Ynew, Yerr, &
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
    saprc99_mosaic_4bin_vbs9_ros_ErrorNorm = Err
  END FUNCTION saprc99_mosaic_4bin_vbs9_ros_ErrorNorm
  SUBROUTINE saprc99_mosaic_4bin_vbs9_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   IMPLICIT NONE
   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)
   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)
   INTEGER, INTENT(INOUT) ::Nfun
   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp
   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL saprc99_mosaic_4bin_vbs9_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL saprc99_mosaic_4bin_vbs9_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL saprc99_mosaic_4bin_vbs9_WSCAL(NVAR,(ONE/Delta),dFdT,1)
  END SUBROUTINE saprc99_mosaic_4bin_vbs9_ros_FunTimeDeriv
  SUBROUTINE saprc99_mosaic_4bin_vbs9_ros_PrepareMatrix ( H, Direction, gam, &
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
     CALL saprc99_mosaic_4bin_vbs9_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL saprc99_mosaic_4bin_vbs9_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO
     CALL saprc99_mosaic_4bin_vbs9_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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
  END SUBROUTINE saprc99_mosaic_4bin_vbs9_ros_PrepareMatrix
  SUBROUTINE saprc99_mosaic_4bin_vbs9_ros_Decomp( A, Pivot, ising, Ndec )
   IMPLICIT NONE
   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)
   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec
CALL decomp_saprc99_mosaic_4bin_vbs9 ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1
  END SUBROUTINE saprc99_mosaic_4bin_vbs9_ros_Decomp
  SUBROUTINE saprc99_mosaic_4bin_vbs9_ros_Solve( A, Pivot, b, Nsol )
   IMPLICIT NONE
   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)
   INTEGER, INTENT(INOUT) :: nsol
   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)
   CALL saprc99_mosaic_4bin_vbs9_KppSolve( A, b )
   Nsol = Nsol+1
  END SUBROUTINE saprc99_mosaic_4bin_vbs9_ros_Solve
  SUBROUTINE saprc99_mosaic_4bin_vbs9_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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
 END SUBROUTINE saprc99_mosaic_4bin_vbs9_Ros2
  SUBROUTINE saprc99_mosaic_4bin_vbs9_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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
  END SUBROUTINE saprc99_mosaic_4bin_vbs9_Ros3
  SUBROUTINE saprc99_mosaic_4bin_vbs9_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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
  END SUBROUTINE saprc99_mosaic_4bin_vbs9_Ros4
  SUBROUTINE saprc99_mosaic_4bin_vbs9_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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
  END SUBROUTINE saprc99_mosaic_4bin_vbs9_Rodas3
  SUBROUTINE saprc99_mosaic_4bin_vbs9_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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
  END SUBROUTINE saprc99_mosaic_4bin_vbs9_Rodas4
END SUBROUTINE saprc99_mosaic_4bin_vbs9_Rosenbrock
SUBROUTINE saprc99_mosaic_4bin_vbs9_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )
   USE saprc99_mosaic_4bin_vbs9_Parameters
   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)
   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun
   CALL saprc99_mosaic_4bin_vbs9_Fun( Y, FIX, RCONST, Ydot )
   Nfun = Nfun+1
END SUBROUTINE saprc99_mosaic_4bin_vbs9_FunTemplate
SUBROUTINE saprc99_mosaic_4bin_vbs9_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )
 USE saprc99_mosaic_4bin_vbs9_Parameters
 USE saprc99_mosaic_4bin_vbs9_Jacobian
    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)
    INTEGER :: Njac
    REAL(kind=dp) :: Jcb(LU_NONZERO)
    REAL(kind=dp) :: Told
    CALL saprc99_mosaic_4bin_vbs9_Jac_SP( Y, FIX, RCONST, Jcb )
    Njac = Njac+1
END SUBROUTINE saprc99_mosaic_4bin_vbs9_JacTemplate
SUBROUTINE saprc99_mosaic_4bin_vbs9_Fun ( V, F, RCT, Vdot )
  REAL(kind=dp) :: V(NVAR)
  REAL(kind=dp) :: F(NFIX)
  REAL(kind=dp) :: RCT(NREACT)
  REAL(kind=dp) :: Vdot(NVAR)
  REAL(kind=dp) :: A(NREACT)
  A(1) = RCT(1)*V(167)
  A(2) = RCT(2)*V(152)*F(2)
  A(3) = RCT(3)*V(152)*V(156)
  A(4) = RCT(4)*V(152)*V(160)*F(2)
  A(5) = RCT(5)*V(152)*V(167)
  A(6) = RCT(6)*V(152)*V(167)
  A(7) = RCT(7)*V(156)*V(160)
  A(8) = RCT(8)*V(156)*V(167)
  A(9) = RCT(9)*V(157)*V(160)
  A(10) = RCT(10)*V(160)*V(160)*F(2)
  A(11) = RCT(11)*V(157)*V(167)
  A(12) = RCT(12)*V(110)
  A(13) = RCT(13)*V(110)*F(1)
  A(14) = RCT(14)*V(157)*V(167)
  A(15) = RCT(15)*V(157)
  A(16) = RCT(16)*V(157)
  A(17) = RCT(17)*V(156)
  A(18) = RCT(18)*V(156)
  A(19) = RCT(19)*V(90)*F(1)
  A(20) = RCT(20)*V(90)*F(2)
  A(21) = RCT(21)*V(160)*V(165)
  A(22) = RCT(22)*V(111)
  A(23) = RCT(23)*V(111)
  A(24) = RCT(24)*V(111)*V(165)
  A(25) = RCT(25)*V(165)*V(167)
  A(26) = RCT(26)*V(157)*V(165)
  A(27) = RCT(27)*V(134)*V(165)
  A(28) = RCT(28)*V(134)
  A(29) = RCT(29)*V(133)*V(165)
  A(30) = RCT(30)*V(156)*V(165)
  A(31) = RCT(31)*V(160)*V(161)
  A(32) = RCT(32)*V(161)*V(167)
  A(33) = RCT(33)*V(119)
  A(34) = RCT(34)*V(119)
  A(35) = RCT(35)*V(119)*V(165)
  A(36) = RCT(36)*V(156)*V(161)
  A(37) = RCT(37)*V(161)*V(161)
  A(38) = RCT(38)*V(161)*V(161)*F(1)
  A(39) = RCT(39)*V(157)*V(161)
  A(40) = RCT(40)*V(157)*V(157)
  A(41) = RCT(41)*V(106)
  A(42) = RCT(42)*V(106)*V(165)
  A(43) = RCT(43)*V(161)*V(165)
  A(44) = RCT(44)*V(100)*V(165)
  A(45) = RCT(45)*V(165)*F(2)
  A(46) = RCT(46)*V(160)*V(168)
  A(47) = RCT(47)*V(161)*V(168)
  A(48) = RCT(48)*V(157)*V(168)
  A(49) = RCT(49)*V(168)*V(168)
  A(50) = RCT(50)*V(168)*V(168)
  A(51) = RCT(51)*V(160)*V(166)
  A(52) = RCT(52)*V(161)*V(166)
  A(53) = RCT(53)*V(157)*V(166)
  A(54) = RCT(54)*V(166)*V(168)
  A(55) = RCT(55)*V(166)*V(166)
  A(56) = RCT(56)*V(141)*V(160)
  A(57) = RCT(57)*V(141)*V(161)
  A(58) = RCT(58)*V(141)*V(157)
  A(59) = RCT(59)*V(141)*V(168)
  A(60) = RCT(60)*V(141)*V(166)
  A(61) = RCT(61)*V(141)*V(141)
  A(62) = RCT(62)*V(158)*V(160)
  A(63) = RCT(63)*V(158)*V(161)
  A(64) = RCT(64)*V(158)*V(168)
  A(65) = RCT(65)*V(157)*V(158)
  A(66) = RCT(66)*V(158)*V(166)
  A(67) = RCT(67)*V(141)*V(158)
  A(68) = RCT(68)*V(158)*V(158)
  A(69) = RCT(69)*V(163)*V(167)
  A(70) = RCT(70)*V(104)
  A(71) = RCT(71)*V(160)*V(163)
  A(72) = RCT(72)*V(161)*V(163)
  A(73) = RCT(73)*V(157)*V(163)
  A(74) = RCT(74)*V(163)*V(168)
  A(75) = RCT(75)*V(163)*V(166)
  A(76) = RCT(76)*V(141)*V(163)
  A(77) = RCT(77)*V(158)*V(163)
  A(78) = RCT(78)*V(163)*V(163)
  A(79) = RCT(79)*V(162)*V(167)
  A(80) = RCT(80)*V(102)
  A(81) = RCT(81)*V(160)*V(162)
  A(82) = RCT(82)*V(161)*V(162)
  A(83) = RCT(83)*V(157)*V(162)
  A(84) = RCT(84)*V(162)*V(168)
  A(85) = RCT(85)*V(162)*V(166)
  A(86) = RCT(86)*V(141)*V(162)
  A(87) = RCT(87)*V(158)*V(162)
  A(88) = RCT(88)*V(162)*V(163)
  A(89) = RCT(89)*V(162)*V(162)
  A(90) = RCT(90)*V(164)*V(167)
  A(91) = RCT(91)*V(105)
  A(92) = RCT(92)*V(160)*V(164)
  A(93) = RCT(93)*V(161)*V(164)
  A(94) = RCT(94)*V(157)*V(164)
  A(95) = RCT(95)*V(164)*V(168)
  A(96) = RCT(96)*V(164)*V(166)
  A(97) = RCT(97)*V(141)*V(164)
  A(98) = RCT(98)*V(158)*V(164)
  A(99) = RCT(99)*V(163)*V(164)
  A(100) = RCT(100)*V(162)*V(164)
  A(101) = RCT(101)*V(164)*V(164)
  A(102) = RCT(102)*V(159)*V(167)
  A(103) = RCT(103)*V(103)
  A(104) = RCT(104)*V(159)*V(160)
  A(105) = RCT(105)*V(159)*V(161)
  A(106) = RCT(106)*V(157)*V(159)
  A(107) = RCT(107)*V(159)*V(168)
  A(108) = RCT(108)*V(159)*V(166)
  A(109) = RCT(109)*V(141)*V(159)
  A(110) = RCT(110)*V(158)*V(159)
  A(111) = RCT(111)*V(159)*V(163)
  A(112) = RCT(112)*V(159)*V(162)
  A(113) = RCT(113)*V(159)*V(164)
  A(114) = RCT(114)*V(159)*V(159)
  A(115) = RCT(115)*V(114)*V(167)
  A(116) = RCT(116)*V(114)
  A(117) = RCT(117)*V(138)*V(167)
  A(118) = RCT(118)*V(138)*V(161)
  A(119) = RCT(119)*V(138)
  A(120) = RCT(120)*V(118)*V(167)
  A(121) = RCT(121)*V(118)*V(161)
  A(122) = RCT(122)*V(118)
  A(123) = RCT(123)*V(150)
  A(124) = RCT(124)*V(150)
  A(125) = RCT(125)*V(150)*V(165)
  A(126) = RCT(126)*V(150)*V(161)
  A(127) = RCT(127)*V(117)
  A(128) = RCT(128)*V(117)*V(160)
  A(129) = RCT(129)*V(150)*V(157)
  A(130) = RCT(130)*V(149)*V(165)
  A(131) = RCT(131)*V(149)
  A(132) = RCT(132)*V(149)*V(157)
  A(133) = RCT(133)*V(153)*V(165)
  A(134) = RCT(134)*V(153)
  A(135) = RCT(135)*V(153)*V(157)
  A(136) = RCT(136)*V(136)*V(165)
  A(137) = RCT(137)*V(136)
  A(138) = RCT(138)*V(154)*V(165)
  A(139) = RCT(139)*V(154)
  A(140) = RCT(140)*V(120)*V(165)
  A(141) = RCT(141)*V(108)*V(165)
  A(142) = RCT(142)*V(116)*V(165)
  A(143) = RCT(143)*V(116)
  A(144) = RCT(144)*V(129)*V(165)
  A(145) = RCT(145)*V(129)
  A(146) = RCT(146)*V(145)
  A(147) = RCT(147)*V(145)
  A(148) = RCT(148)*V(145)*V(165)
  A(149) = RCT(149)*V(145)*V(157)
  A(150) = RCT(150)*V(132)
  A(151) = RCT(151)*V(132)*V(165)
  A(152) = RCT(152)*V(132)*V(157)
  A(153) = RCT(153)*V(109)
  A(154) = RCT(154)*V(131)*V(165)
  A(155) = RCT(155)*V(131)*V(157)
  A(156) = RCT(156)*V(124)*V(165)
  A(157) = RCT(157)*V(124)*V(157)
  A(158) = RCT(158)*V(128)*V(157)
  A(159) = RCT(159)*V(130)*V(165)
  A(160) = RCT(160)*V(130)
  A(161) = RCT(161)*V(130)*V(157)
  A(162) = RCT(162)*V(142)*V(165)
  A(163) = RCT(163)*V(142)*V(156)
  A(164) = RCT(164)*V(142)*V(157)
  A(165) = RCT(165)*V(142)*V(152)
  A(166) = RCT(166)*V(142)
  A(167) = RCT(167)*V(148)*V(165)
  A(168) = RCT(168)*V(148)*V(156)
  A(169) = RCT(169)*V(148)*V(152)
  A(170) = RCT(170)*V(148)
  A(171) = RCT(171)*V(146)*V(165)
  A(172) = RCT(172)*V(146)*V(156)
  A(173) = RCT(173)*V(146)*V(157)
  A(174) = RCT(174)*V(146)
  A(175) = RCT(175)*V(155)*V(165)
  A(176) = RCT(176)*V(155)
  A(177) = RCT(177)*V(151)*V(165)
  A(178) = RCT(178)*V(151)
  A(179) = RCT(179)*V(127)*V(165)
  A(180) = RCT(180)*V(127)*V(156)
  A(181) = RCT(181)*V(122)*V(165)
  A(182) = RCT(182)*V(122)
  A(183) = RCT(183)*V(123)*V(165)
  A(184) = RCT(184)*V(123)
  A(185) = RCT(185)*V(99)*V(165)
  A(186) = RCT(186)*V(135)*V(165)
  A(187) = RCT(187)*V(135)*V(156)
  A(188) = RCT(188)*V(135)*V(157)
  A(189) = RCT(189)*V(135)*V(152)
  A(190) = RCT(190)*V(140)*V(165)
  A(191) = RCT(191)*V(140)*V(156)
  A(192) = RCT(192)*V(140)*V(157)
  A(193) = RCT(193)*V(140)*V(152)
  A(194) = RCT(194)*V(144)*V(165)
  A(195) = RCT(195)*V(144)*V(156)
  A(196) = RCT(196)*V(144)*V(157)
  A(197) = RCT(197)*V(144)*V(152)
  A(198) = RCT(198)*V(143)*V(165)
  A(199) = RCT(199)*V(143)*V(156)
  A(200) = RCT(200)*V(143)*V(157)
  A(201) = RCT(201)*V(143)*V(152)
  A(202) = RCT(202)*V(101)*V(165)
  A(203) = RCT(203)*V(107)*V(165)
  A(204) = RCT(204)*V(126)*V(165)
  A(205) = RCT(205)*V(113)*V(165)
  A(206) = RCT(206)*V(125)*V(165)
  A(207) = RCT(207)*V(112)*V(165)
  A(208) = RCT(208)*V(121)*V(165)
  A(209) = RCT(209)*V(115)*V(165)
  A(210) = RCT(210)*V(139)*V(165)
  A(211) = RCT(211)*V(139)*V(156)
  A(212) = RCT(212)*V(139)*V(157)
  A(213) = RCT(213)*V(139)*V(152)
  A(214) = RCT(214)*V(147)*V(165)
  A(215) = RCT(215)*V(147)*V(156)
  A(216) = RCT(216)*V(147)*V(157)
  A(217) = RCT(217)*V(147)*V(152)
  A(218) = RCT(218)*V(126)*V(156)
  A(219) = RCT(219)*V(137)*V(165)
  A(220) = RCT(220)*V(137)*V(156)
  A(221) = RCT(221)*V(137)*V(157)
  A(222) = RCT(222)*V(137)*V(152)
  A(223) = RCT(223)*V(100)
  A(224) = RCT(224)*V(161)
  A(225) = RCT(225)*V(100)
  A(226) = RCT(226)*V(1)
  A(227) = RCT(227)*V(134)
  A(228) = RCT(228)*V(106)
  A(229) = RCT(229)*V(2)
  A(230) = RCT(230)*V(125)*V(165)
  A(231) = RCT(231)*V(125)*V(165)
  A(232) = RCT(232)*V(112)*V(165)
  A(233) = RCT(233)*V(112)*V(165)
  A(234) = RCT(234)*V(139)*V(165)
  A(235) = RCT(235)*V(139)*V(165)
  A(236) = RCT(236)*V(139)*V(165)
  A(237) = RCT(237)*V(139)*V(165)
  A(238) = RCT(238)*V(139)*V(165)
  A(239) = RCT(239)*V(139)*V(165)
  A(240) = RCT(240)*V(139)*V(165)
  A(241) = RCT(241)*V(139)*V(165)
  A(242) = RCT(242)*V(147)*V(165)
  A(243) = RCT(243)*V(147)*V(165)
  A(244) = RCT(244)*V(147)*V(165)
  A(245) = RCT(245)*V(147)*V(165)
  A(246) = RCT(246)*V(147)*V(165)
  A(247) = RCT(247)*V(147)*V(165)
  A(248) = RCT(248)*V(147)*V(165)
  A(249) = RCT(249)*V(147)*V(165)
  A(250) = RCT(250)*V(121)*V(165)
  A(251) = RCT(251)*V(121)*V(165)
  A(252) = RCT(252)*V(121)*V(165)
  A(253) = RCT(253)*V(121)*V(165)
  A(254) = RCT(254)*V(121)*V(165)
  A(255) = RCT(255)*V(121)*V(165)
  A(256) = RCT(256)*V(121)*V(165)
  A(257) = RCT(257)*V(121)*V(165)
  A(258) = RCT(258)*V(115)*V(165)
  A(259) = RCT(259)*V(115)*V(165)
  A(260) = RCT(260)*V(115)*V(165)
  A(261) = RCT(261)*V(115)*V(165)
  A(262) = RCT(262)*V(115)*V(165)
  A(263) = RCT(263)*V(115)*V(165)
  A(264) = RCT(264)*V(115)*V(165)
  A(265) = RCT(265)*V(115)*V(165)
  A(266) = RCT(266)*V(140)*V(165)
  A(267) = RCT(267)*V(140)*V(165)
  A(268) = RCT(268)*V(140)*V(165)
  A(269) = RCT(269)*V(140)*V(165)
  A(270) = RCT(270)*V(140)*V(165)
  A(271) = RCT(271)*V(140)*V(165)
  A(272) = RCT(272)*V(144)*V(165)
  A(273) = RCT(273)*V(144)*V(165)
  A(274) = RCT(274)*V(144)*V(165)
  A(275) = RCT(275)*V(144)*V(165)
  A(276) = RCT(276)*V(144)*V(165)
  A(277) = RCT(277)*V(144)*V(165)
  A(278) = RCT(278)*V(144)*V(165)
  A(279) = RCT(279)*V(144)*V(165)
  A(280) = RCT(280)*V(143)*V(165)
  A(281) = RCT(281)*V(143)*V(165)
  A(282) = RCT(282)*V(143)*V(165)
  A(283) = RCT(283)*V(143)*V(165)
  A(284) = RCT(284)*V(3)*V(165)
  A(285) = RCT(285)*V(4)*V(165)
  A(286) = RCT(286)*V(74)*V(165)
  A(287) = RCT(287)*V(34)*V(165)
  A(288) = RCT(288)*V(5)*V(165)
  A(289) = RCT(289)*V(6)*V(165)
  A(290) = RCT(290)*V(91)*V(165)
  A(291) = RCT(291)*V(50)*V(165)
  A(292) = RCT(292)*V(66)*V(165)
  A(293) = RCT(293)*V(35)*V(165)
  A(294) = RCT(294)*V(75)*V(165)
  A(295) = RCT(295)*V(36)*V(165)
  A(296) = RCT(296)*V(67)*V(165)
  A(297) = RCT(297)*V(37)*V(165)
  A(298) = RCT(298)*V(76)*V(165)
  A(299) = RCT(299)*V(38)*V(165)
  A(300) = RCT(300)*V(68)*V(165)
  A(301) = RCT(301)*V(39)*V(165)
  A(302) = RCT(302)*V(77)*V(165)
  A(303) = RCT(303)*V(40)*V(165)
  A(304) = RCT(304)*V(69)*V(165)
  A(305) = RCT(305)*V(41)*V(165)
  A(306) = RCT(306)*V(78)*V(165)
  A(307) = RCT(307)*V(42)*V(165)
  A(308) = RCT(308)*V(70)*V(165)
  A(309) = RCT(309)*V(43)*V(165)
  A(310) = RCT(310)*V(79)*V(165)
  A(311) = RCT(311)*V(44)*V(165)
  A(312) = RCT(312)*V(71)*V(165)
  A(313) = RCT(313)*V(45)*V(165)
  A(314) = RCT(314)*V(80)*V(165)
  A(315) = RCT(315)*V(46)*V(165)
  A(316) = RCT(316)*V(72)*V(165)
  A(317) = RCT(317)*V(47)*V(165)
  A(318) = RCT(318)*V(81)*V(165)
  A(319) = RCT(319)*V(48)*V(165)
  A(320) = RCT(320)*V(73)*V(165)
  A(321) = RCT(321)*V(49)*V(165)
  A(322) = RCT(322)*V(82)*V(165)
  A(323) = RCT(323)*V(51)*V(165)
  A(324) = RCT(324)*V(92)*V(165)
  A(325) = RCT(325)*V(52)*V(165)
  A(326) = RCT(326)*V(83)*V(165)
  A(327) = RCT(327)*V(53)*V(165)
  A(328) = RCT(328)*V(93)*V(165)
  A(329) = RCT(329)*V(54)*V(165)
  A(330) = RCT(330)*V(84)*V(165)
  A(331) = RCT(331)*V(55)*V(165)
  A(332) = RCT(332)*V(94)*V(165)
  A(333) = RCT(333)*V(56)*V(165)
  A(334) = RCT(334)*V(85)*V(165)
  A(335) = RCT(335)*V(57)*V(165)
  A(336) = RCT(336)*V(95)*V(165)
  A(337) = RCT(337)*V(58)*V(165)
  A(338) = RCT(338)*V(86)*V(165)
  A(339) = RCT(339)*V(59)*V(165)
  A(340) = RCT(340)*V(96)*V(165)
  A(341) = RCT(341)*V(60)*V(165)
  A(342) = RCT(342)*V(87)*V(165)
  A(343) = RCT(343)*V(61)*V(165)
  A(344) = RCT(344)*V(97)*V(165)
  A(345) = RCT(345)*V(62)*V(165)
  A(346) = RCT(346)*V(88)*V(165)
  A(347) = RCT(347)*V(63)*V(165)
  A(348) = RCT(348)*V(98)*V(165)
  A(349) = RCT(349)*V(64)*V(165)
  A(350) = RCT(350)*V(89)*V(165)
  A(351) = RCT(351)*V(65)*V(165)
  Vdot(1) = A(44)+A(223)-A(226)
  Vdot(2) = 0.5*A(218)+0.135*A(220)-A(229)
  Vdot(3) = 0
  Vdot(4) = 0
  Vdot(5) = 0
  Vdot(6) = 0
  Vdot(7) = A(128)+0.333*A(163)+0.351*A(168)+0.1*A(172)+0.37*A(187)+0.204*A(191)+0.103*A(195)+0.103*A(199)+0.297*A(204)&
              &+0.185*A(211)+0.073*A(215)+0.185*A(220)
  Vdot(8) = 0.25*A(72)+A(74)+A(75)+A(77)+0.05*A(211)+0.129*A(215)+0.17*A(220)
  Vdot(9) = 0.25*A(82)+A(84)+A(85)+A(87)+0.25*A(93)+A(95)+A(96)+A(98)+0.25*A(105)+A(107)+A(108)+2*A(110)+0.372*A(172)&
              &+0.15*A(191)+0.189*A(195)+0.189*A(199)+0.119*A(211)+0.247*A(215)
  Vdot(10) = 0.75*A(72)
  Vdot(11) = 0.75*A(82)+0.75*A(93)+0.75*A(105)
  Vdot(12) = 2*A(120)+A(221)
  Vdot(13) = 6*A(120)+7*A(160)+0.048*A(219)+0.07*A(220)+2.693*A(221)+0.55*A(222)
  Vdot(14) = A(46)+A(48)+A(51)+A(53)+A(56)+A(58)+A(62)+A(65)
  Vdot(15) = A(47)+A(49)+A(50)+A(52)+A(54)+A(55)+A(57)+A(59)+A(60)+A(61)+A(63)+A(64)+A(66)+A(67)+A(68)
  Vdot(16) = A(234)+A(242)+A(254)+A(262)
  Vdot(17) = A(230)+A(232)+A(235)+A(243)+A(255)+A(263)
  Vdot(18) = A(236)+A(244)+A(256)+A(264)
  Vdot(19) = A(237)+A(245)+A(257)+A(265)
  Vdot(20) = A(238)+A(246)+A(250)+A(258)
  Vdot(21) = A(231)+A(233)+A(239)+A(247)+A(251)+A(259)
  Vdot(22) = A(240)+A(248)+A(252)+A(260)
  Vdot(23) = A(241)+A(249)+A(253)+A(261)
  Vdot(24) = A(266)+A(276)
  Vdot(25) = A(267)+A(277)
  Vdot(26) = A(268)+A(278)
  Vdot(27) = A(279)
  Vdot(28) = A(269)+A(272)+A(280)
  Vdot(29) = A(270)+A(273)+A(281)
  Vdot(30) = A(271)+A(274)+A(282)
  Vdot(31) = A(275)+A(283)
  Vdot(32) = A(21)+A(24)+A(25)+A(26)+A(27)+A(29)+A(30)+A(43)+A(44)+A(45)+A(125)+A(130)+A(133)+A(136)+A(138)+A(140)&
               &+A(141)+A(142)+A(144)+A(148)+A(151)+A(154)+A(156)+A(159)+A(162)+A(167)+A(171)+A(175)+A(177)+A(179)+A(181)&
               &+A(183)+A(185)+A(186)+A(190)+A(194)+A(198)+A(202)+A(203)+A(204)+A(205)+A(206)+A(207)+A(208)+A(209)+A(210)&
               &+A(214)+A(219)
  Vdot(33) = A(284)+A(285)+A(286)+A(287)+A(288)+A(289)+A(290)+A(291)+A(292)+A(293)+A(294)+A(295)+A(296)+A(297)+A(298)&
               &+A(299)+A(300)+A(301)+A(302)+A(303)+A(304)+A(305)+A(306)+A(307)+A(308)+A(309)+A(310)+A(311)+A(312)+A(313)&
               &+A(314)+A(315)+A(316)+A(317)+A(318)+A(319)+A(320)+A(321)+A(322)+A(323)+A(324)+A(325)+A(326)+A(327)+A(328)&
               &+A(329)+A(330)+A(331)+A(332)+A(333)+A(334)+A(335)+A(336)+A(337)+A(338)+A(339)+A(340)+A(341)+A(342)+A(343)&
               &+A(344)+A(345)+A(346)+A(347)+A(348)+A(349)+A(350)+A(351)
  Vdot(34) = 0.15*A(292)+A(293)+0.15*A(294)+A(295)
  Vdot(35) = -A(293)
  Vdot(36) = -A(295)+0.15*A(296)+A(297)+0.15*A(298)+A(299)
  Vdot(37) = -A(297)
  Vdot(38) = -A(299)+0.15*A(300)+A(301)+0.15*A(302)+A(303)
  Vdot(39) = -A(301)
  Vdot(40) = -A(303)+0.15*A(304)+A(305)+0.15*A(306)+A(307)
  Vdot(41) = -A(305)
  Vdot(42) = -A(307)+0.15*A(308)+A(309)+0.15*A(310)+A(311)
  Vdot(43) = -A(309)
  Vdot(44) = -A(311)+0.15*A(312)+A(313)+0.15*A(314)+A(315)
  Vdot(45) = -A(313)
  Vdot(46) = -A(315)+0.15*A(316)+A(317)+0.15*A(318)+A(319)
  Vdot(47) = -A(317)
  Vdot(48) = -A(319)+0.15*A(320)+A(321)
  Vdot(49) = -A(321)
  Vdot(50) = 0.15*A(322)+A(323)+0.15*A(324)+A(325)
  Vdot(51) = -A(323)
  Vdot(52) = -A(325)+0.15*A(326)+A(327)+0.15*A(328)+A(329)
  Vdot(53) = -A(327)
  Vdot(54) = -A(329)+0.15*A(330)+A(331)+0.15*A(332)+A(333)
  Vdot(55) = -A(331)
  Vdot(56) = -A(333)+0.15*A(334)+A(335)+0.15*A(336)+A(337)
  Vdot(57) = -A(335)
  Vdot(58) = -A(337)+0.15*A(338)+A(339)+0.15*A(340)+A(341)
  Vdot(59) = -A(339)
  Vdot(60) = -A(341)+0.15*A(342)+A(343)+0.15*A(344)+A(345)
  Vdot(61) = -A(343)
  Vdot(62) = -A(345)+0.15*A(346)+A(347)+0.15*A(348)+A(349)
  Vdot(63) = -A(347)
  Vdot(64) = -A(349)+0.15*A(350)+A(351)
  Vdot(65) = -A(351)
  Vdot(66) = -A(292)
  Vdot(67) = -A(296)
  Vdot(68) = -A(300)
  Vdot(69) = -A(304)
  Vdot(70) = -A(308)
  Vdot(71) = -A(312)
  Vdot(72) = -A(316)
  Vdot(73) = -A(320)
  Vdot(74) = A(292)+A(294)
  Vdot(75) = -A(294)+A(296)+A(298)
  Vdot(76) = -A(298)+A(300)+A(302)
  Vdot(77) = -A(302)+A(304)+A(306)
  Vdot(78) = -A(306)+A(308)+A(310)
  Vdot(79) = -A(310)+A(312)+A(314)
  Vdot(80) = -A(314)+A(316)+A(318)
  Vdot(81) = -A(318)+A(320)
  Vdot(82) = -A(322)
  Vdot(83) = -A(326)
  Vdot(84) = -A(330)
  Vdot(85) = -A(334)
  Vdot(86) = -A(338)
  Vdot(87) = -A(342)
  Vdot(88) = -A(346)
  Vdot(89) = -A(350)
  Vdot(90) = A(18)-A(19)-A(20)
  Vdot(91) = A(322)+A(324)
  Vdot(92) = -A(324)+A(326)+A(328)
  Vdot(93) = -A(328)+A(330)+A(332)
  Vdot(94) = -A(332)+A(334)+A(336)
  Vdot(95) = -A(336)+A(338)+A(340)
  Vdot(96) = -A(340)+A(342)+A(344)
  Vdot(97) = -A(344)+A(346)+A(348)
  Vdot(98) = -A(348)+A(350)
  Vdot(99) = -A(185)
  Vdot(100) = -A(44)-A(223)-A(225)
  Vdot(101) = -A(202)
  Vdot(102) = A(79)-A(80)
  Vdot(103) = A(102)-A(103)
  Vdot(104) = A(69)-A(70)
  Vdot(105) = A(90)-A(91)
  Vdot(106) = A(37)+A(38)-A(41)-A(42)-A(228)
  Vdot(107) = -A(203)
  Vdot(108) = -A(141)
  Vdot(109) = -A(153)+0.031*A(195)+0.031*A(199)+0.087*A(209)
  Vdot(110) = A(11)-A(12)-A(13)
  Vdot(111) = A(21)-A(22)-A(23)-A(24)
  Vdot(112) = -A(207)
  Vdot(113) = -A(205)
  Vdot(114) = -A(115)-A(116)+0.236*A(205)
  Vdot(115) = -A(209)
  Vdot(116) = A(47)-A(142)-A(143)
  Vdot(117) = A(126)-A(127)-A(128)
  Vdot(118) = -A(120)-A(121)-A(122)+A(158)
  Vdot(119) = A(32)-A(33)-A(34)-A(35)
  Vdot(120) = A(49)+0.25*A(54)+0.25*A(64)-A(140)
  Vdot(121) = -A(208)
  Vdot(122) = -A(181)-A(182)+0.108*A(208)+0.099*A(209)
  Vdot(123) = -A(183)-A(184)+0.051*A(208)+0.093*A(209)
  Vdot(124) = -A(156)-A(157)+0.207*A(208)+0.187*A(209)
  Vdot(125) = -A(206)
  Vdot(126) = -A(204)-A(218)
  Vdot(127) = -A(179)-A(180)+0.491*A(208)+0.561*A(209)
  Vdot(128) = A(117)+A(121)+A(122)-A(158)
  Vdot(129) = A(52)+A(63)-A(144)-A(145)
  Vdot(130) = -A(159)-A(160)-A(161)+0.059*A(208)+0.05*A(209)+0.061*A(214)+0.042*A(215)+0.015*A(216)
  Vdot(131) = A(118)+A(119)-A(154)-A(155)+0.017*A(208)
  Vdot(132) = -A(150)-A(151)-A(152)+0.23*A(156)+0.084*A(162)+0.9*A(163)+0.3*A(167)+0.95*A(168)+0.174*A(171)+0.742*A(172)&
                &+0.008*A(173)+0.5*A(182)+0.5*A(184)+0.119*A(208)+0.287*A(209)
  Vdot(133) = -A(29)+A(123)+A(124)+A(125)+A(129)+A(131)+0.034*A(133)+A(134)+2*A(146)+A(147)+1.26*A(148)+1.26*A(149)&
                &+A(150)+A(151)+A(152)+0.416*A(162)+0.45*A(163)+0.5*A(164)+0.67*A(166)+0.475*A(168)+0.7*A(170)+0.336*A(171)&
                &+0.498*A(172)+0.572*A(173)+1.233*A(174)+A(179)+1.5*A(180)+A(182)+A(184)+0.5*A(187)+0.491*A(189)+0.275&
                &*A(191)+0.157*A(195)+0.157*A(199)+0.393*A(204)+0.002*A(206)+0.345*A(211)+0.265*A(215)+0.012*A(217)+1.5&
                &*A(218)+0.51*A(220)
  Vdot(134) = 2*A(13)+A(25)-A(27)-A(28)+0.2*A(39)+A(129)+A(132)+A(135)+A(149)+A(152)+A(155)+A(157)+A(158)+A(161)+0.5&
                &*A(164)+0.15*A(173)-A(227)
  Vdot(135) = -A(186)-A(187)-A(188)-A(189)
  Vdot(136) = A(116)-A(136)-A(137)+0.006*A(177)+0.02*A(178)+0.13*A(195)+0.13*A(199)+0.704*A(203)+0.024*A(205)+0.452&
                &*A(206)+0.072*A(207)+0.005*A(210)+0.001*A(211)+0.024*A(212)+0.127*A(214)+0.045*A(215)+0.102*A(216)
  Vdot(137) = -A(219)-A(220)-A(221)-A(222)
  Vdot(138) = A(92)+A(94)+A(99)+A(100)+2*A(101)+A(113)-A(117)-A(118)-A(119)+0.24*A(154)+A(155)+0.24*A(156)+A(157)
  Vdot(139) = -A(210)-A(211)-A(212)-A(213)
  Vdot(140) = -A(190)-A(191)-A(192)-A(193)
  Vdot(141) = -A(56)-A(57)-A(58)-A(59)-A(60)-A(67)-A(76)-A(86)+A(92)+A(94)-A(97)+A(99)+A(100)+2*A(101)-A(109)+A(113)&
                &+A(136)+0.616*A(138)+0.675*A(167)+0.515*A(176)+0.596*A(177)+0.152*A(178)+A(181)+A(182)+A(183)+A(184)+0.079&
                &*A(190)+0.126*A(191)+0.187*A(192)+0.24*A(193)+0.5*A(194)+0.729*A(195)+0.75*A(196)+0.5*A(198)+0.729*A(199)&
                &+0.75*A(200)+0.559*A(205)+0.936*A(206)+0.948*A(207)+0.205*A(210)+0.488*A(212)+0.001*A(214)+0.137*A(215)&
                &+0.711*A(216)
  Vdot(142) = -A(162)-A(163)-A(164)-A(165)-A(166)+0.23*A(190)+0.39*A(191)+0.025*A(214)+0.026*A(215)+0.012*A(217)
  Vdot(143) = -A(198)-A(199)-A(200)-A(201)
  Vdot(144) = -A(194)-A(195)-A(196)-A(197)
  Vdot(145) = -A(146)-A(147)-A(148)-A(149)+0.23*A(154)+0.15*A(171)+0.023*A(172)+A(180)+0.5*A(182)+0.5*A(184)+0.009&
                &*A(189)+0.001*A(195)+0.001*A(199)+0.607*A(204)+0.118*A(208)+0.097*A(209)
  Vdot(146) = -A(171)-A(172)-A(173)-A(174)+0.357*A(190)+0.936*A(192)+0.025*A(214)
  Vdot(147) = -A(214)-A(215)-A(216)-A(217)
  Vdot(148) = -A(167)-A(168)-A(169)-A(170)+0.32*A(190)+0.16*A(191)+0.019*A(215)+0.048*A(216)
  Vdot(149) = A(81)+A(83)+A(88)+2*A(89)+A(100)+A(112)-A(130)-A(131)-A(132)+0.034*A(133)+A(134)+0.482*A(138)+A(139)+0.96&
                &*A(141)+0.129*A(171)+0.047*A(172)+0.467*A(174)+0.084*A(175)+0.246*A(176)+0.439*A(177)+0.431*A(178)+0.195&
                &*A(186)+0.25*A(189)+A(202)+0.445*A(205)+0.455*A(206)+0.099*A(207)+0.294*A(210)+0.154*A(211)+0.009*A(212)&
                &+0.732*A(214)+0.456*A(215)+0.507*A(216)+0.984*A(219)+0.5*A(220)
  Vdot(150) = A(46)+A(48)+A(49)+2*A(50)+0.75*A(54)+0.75*A(64)+A(74)+A(84)+A(95)+A(104)+A(106)+A(107)+A(111)+A(112)&
                &+A(113)+2*A(114)-A(123)-A(124)-A(125)-A(126)+A(127)-A(129)+A(136)+0.115*A(138)+A(140)+0.081*A(141)+0.35&
                &*A(142)+A(143)+A(147)+0.084*A(162)+0.2*A(163)+0.67*A(166)+0.3*A(167)+0.1*A(168)+0.055*A(171)+0.125*A(172)&
                &+0.227*A(173)+0.3*A(174)+0.213*A(175)+0.506*A(176)+0.01*A(177)+0.134*A(178)+1.61*A(186)+A(187)+0.191*A(189)&
                &+0.624*A(190)+0.592*A(191)+0.24*A(193)+0.276*A(194)+0.235*A(195)+0.276*A(198)+0.235*A(199)+0.096*A(204)&
                &+0.026*A(205)+0.024*A(206)+0.026*A(207)+0.732*A(210)+0.5*A(211)+0.244*A(214)+0.269*A(215)+0.079*A(216)&
                &+0.984*A(219)+0.5*A(220)
  Vdot(151) = A(62)+A(115)+0.572*A(173)-0.69*A(177)-A(178)+0.276*A(196)+0.276*A(200)+0.511*A(212)+0.321*A(216)
  Vdot(152) = A(1)-A(2)-A(3)-A(4)-A(5)-A(6)+A(16)+A(17)+A(20)-A(165)-A(169)-A(189)-A(193)-A(197)-A(201)-A(213)-A(217)&
                &-A(222)
  Vdot(153) = -A(133)-A(134)-A(135)+0.37*A(138)+A(144)+A(145)+A(165)+0.675*A(167)+0.45*A(169)+0.013*A(171)+0.218*A(173)&
                &+0.558*A(175)+0.71*A(176)+0.213*A(177)+0.147*A(178)+A(179)+A(181)+A(183)+A(188)+0.474*A(194)+0.205*A(195)&
                &+0.474*A(196)+0.147*A(197)+0.474*A(198)+0.205*A(199)+0.474*A(200)+0.147*A(201)+0.261*A(203)+0.122*A(205)&
                &+0.244*A(206)+0.204*A(207)+0.497*A(210)+0.363*A(211)+0.037*A(212)+0.45*A(213)+0.511*A(214)+0.305*A(215)&
                &+0.151*A(216)+0.069*A(217)+0.45*A(222)
  Vdot(154) = 0.5*A(64)+A(65)+0.5*A(66)+A(68)-A(138)-A(139)+0.416*A(162)+0.55*A(169)+0.15*A(171)+0.21*A(172)+0.233&
                &*A(174)+0.115*A(175)+0.177*A(177)+0.243*A(178)+0.332*A(205)+0.11*A(206)+0.089*A(207)+0.437*A(213)+0.072&
                &*A(214)+0.026*A(215)+0.001*A(216)+0.659*A(217)+0.55*A(222)
  Vdot(155) = 0.5*A(64)+0.5*A(66)+A(68)+A(77)+A(87)+A(98)+0.7*A(170)+0.332*A(171)-0.671*A(175)-A(176)+0.048*A(177)+0.435&
                &*A(178)+0.1*A(191)+0.75*A(193)+0.276*A(194)+0.276*A(195)+0.853*A(197)+0.276*A(198)+0.276*A(199)+0.853&
                &*A(201)+0.125*A(206)+0.417*A(207)+0.055*A(208)+0.119*A(210)+0.215*A(211)+0.113*A(213)+0.043*A(215)+0.259&
                &*A(217)
  Vdot(156) = A(2)-A(3)-A(7)-A(8)-A(17)-A(18)-A(30)-A(36)+0.25*A(72)+0.25*A(82)+0.25*A(93)+0.25*A(105)-A(163)-A(168)&
                &-A(172)-A(180)-A(187)-A(191)-A(195)-A(199)-A(211)-A(215)-A(218)-A(220)
  Vdot(157) = A(6)+A(8)-A(9)-A(11)+A(12)-A(14)-A(15)-A(16)-A(26)+A(27)+0.39*A(34)-A(39)-2*A(40)-A(48)-A(53)-A(58)-A(65)&
                &-A(73)-A(83)-A(94)-A(106)-A(129)-A(132)-A(135)-A(149)-A(152)-A(155)-A(157)-A(158)-A(161)-A(164)-A(173)&
                &-A(188)-A(192)-A(196)-A(200)-A(212)-A(216)-A(221)
  Vdot(158) = -A(62)-A(63)-A(64)-A(65)-A(66)-2*A(68)-A(77)-A(87)-A(98)-A(110)+0.001*A(133)+0.042*A(138)+0.025*A(167)&
                &+0.041*A(171)+0.051*A(173)+0.07*A(175)+0.04*A(176)+0.173*A(177)+0.095*A(178)+0.093*A(190)+0.008*A(191)&
                &+0.064*A(192)+0.01*A(193)+0.25*A(194)+0.18*A(195)+0.25*A(196)+0.25*A(198)+0.18*A(199)+0.25*A(200)+0.035&
                &*A(203)+0.07*A(205)+0.143*A(206)+0.347*A(207)+0.011*A(208)+0.009*A(209)+0.09*A(210)+0.001*A(211)+0.176&
                &*A(212)+0.082*A(214)+0.002*A(215)+0.136*A(216)+0.001*A(217)+0.016*A(219)+0.051*A(221)
  Vdot(159) = -A(102)+A(103)-A(104)-A(105)-A(106)-A(107)-A(108)-A(110)-A(111)-A(112)-A(113)-2*A(114)+0.5*A(162)+0.5&
                &*A(164)+0.33*A(166)+0.3*A(170)+0.289*A(171)+0.15*A(173)+0.192*A(191)+0.24*A(193)
  Vdot(160) = A(1)-A(4)+A(5)-A(7)-A(9)-2*A(10)+A(14)+A(15)-A(21)+A(22)-A(31)-A(46)-A(51)-A(56)-A(62)-A(71)-A(81)-A(92)&
                &-A(104)-A(128)
  Vdot(161) = A(23)+A(26)+A(29)+A(30)-A(31)-A(32)+A(33)+0.61*A(34)-A(36)-2*A(37)-2*A(38)-A(39)+A(42)-A(43)+A(44)+A(45)&
                &+A(46)-A(47)+A(48)+2*A(50)+A(51)-A(52)+A(53)+A(54)+A(55)-A(63)+A(64)+A(65)+A(66)+A(68)-A(72)-A(82)-A(93)&
                &-A(105)-A(118)-A(121)+2*A(123)+A(125)-A(126)+A(127)+A(128)+A(129)+A(131)+A(134)+A(140)+0.95*A(141)+A(143)&
                &+A(145)+2*A(146)+0.63*A(148)+0.63*A(149)+A(150)+0.008*A(163)+0.34*A(166)+0.064*A(168)+0.4*A(172)+1.233&
                &*A(174)+0.379*A(175)+0.113*A(177)+0.341*A(178)+1.5*A(180)+0.5*A(182)+0.5*A(184)+0.12*A(187)+0.5*A(189)&
                &+0.033*A(195)+0.033*A(199)+0.297*A(204)+0.224*A(208)+0.187*A(209)+0.056*A(211)+0.003*A(215)+0.013*A(217)&
                &+1.5*A(218)+0.06*A(220)-A(224)
  Vdot(162) = -A(79)+A(80)-A(81)-A(82)-A(83)-A(84)-A(85)-A(87)-A(88)-2*A(89)-A(100)-A(112)+0.965*A(133)+A(135)+0.096&
                &*A(138)+0.37*A(148)+0.37*A(149)+0.1*A(163)+0.05*A(168)+0.048*A(172)+0.3*A(174)+0.049*A(175)+0.333*A(176)&
                &+0.201*A(195)+0.201*A(199)+0.006*A(215)
  Vdot(163) = -A(69)+A(70)-A(71)-A(72)-A(73)-A(74)-A(75)-A(77)-2*A(78)-A(88)-A(99)+A(104)+A(106)+A(112)+A(113)+2*A(114)&
                &+A(130)+A(132)+A(136)+A(137)+0.492*A(138)+A(139)+A(150)+A(151)+A(152)+2*A(153)+0.67*A(166)+0.675*A(167)&
                &+0.467*A(174)+0.029*A(175)+0.667*A(176)+A(181)+0.5*A(182)+A(183)+0.5*A(184)+0.123*A(195)+0.123*A(199)+0.011&
                &*A(206)+0.137*A(215)
  Vdot(164) = -A(90)+A(91)-A(92)-A(93)-A(94)-A(95)-A(96)-A(98)-A(99)-A(100)-2*A(101)-A(113)+A(159)+A(161)
  Vdot(165) = 2*A(19)-A(21)+A(22)-A(24)-A(25)-A(26)-A(27)+A(28)-A(29)-A(30)+A(31)+0.39*A(34)-A(35)+A(36)+0.8*A(39)+2&
                &*A(41)-A(42)-A(43)-A(44)-A(45)-A(125)-A(130)-A(133)-A(136)-A(138)-A(140)-A(141)-0.65*A(142)+A(143)-0.34&
                &*A(144)+A(145)-A(148)-A(151)-A(154)-A(156)-A(159)-A(162)+0.208*A(163)+0.33*A(166)-A(167)+0.164*A(168)&
                &-A(171)+0.285*A(172)-A(175)-A(177)-A(179)+0.5*A(180)-A(181)-A(183)-A(185)-A(186)+0.12*A(187)-A(190)+0.266&
                &*A(191)-A(194)+0.567*A(195)-A(198)+0.567*A(199)-A(202)-A(203)-0.397*A(204)-A(205)-A(206)-A(207)-A(208)&
                &-A(209)-A(210)+0.155*A(211)-A(214)+0.378*A(215)+0.5*A(218)-A(219)+0.32*A(220)-A(284)-A(286)-A(288)-A(290)&
                &-A(292)-A(294)-A(296)-A(298)-A(300)-A(302)-A(304)-A(306)-A(308)-A(310)-A(312)-A(314)-A(316)-A(318)-A(320)&
                &-A(322)-A(324)-A(326)-A(328)-A(330)-A(332)-A(334)-A(336)-A(338)-A(340)-A(342)-A(344)-A(346)-A(348)-A(350)
  Vdot(166) = -A(51)-A(52)-A(53)-A(54)-2*A(55)-A(66)-A(75)+A(81)+A(83)-A(85)+A(88)+2*A(89)-A(96)+A(100)-A(108)+A(112)&
                &+0.034*A(133)+A(134)+0.37*A(138)+A(139)+0.05*A(141)+0.34*A(144)+0.76*A(154)+0.76*A(156)+0.5*A(162)+0.1&
                &*A(163)+0.5*A(164)+0.33*A(166)+0.3*A(167)+0.05*A(168)+0.67*A(171)+0.048*A(172)+0.799*A(173)+0.473*A(175)&
                &+0.96*A(176)+0.376*A(177)+0.564*A(178)+A(179)+A(182)+A(184)+A(186)+A(188)+0.2*A(189)+0.907*A(190)+0.066&
                &*A(191)+0.749*A(192)+0.75*A(194)+0.031*A(195)+0.276*A(196)+0.75*A(198)+0.031*A(199)+0.276*A(200)+A(202)&
                &+0.965*A(203)+0.1*A(204)+0.695*A(205)+0.835*A(206)+0.653*A(207)+0.765*A(208)+0.804*A(209)+0.91*A(210)+0.022&
                &*A(211)+0.824*A(212)+0.918*A(214)+0.033*A(215)+0.442*A(216)+0.012*A(217)+0.984*A(219)+0.949*A(221)
  Vdot(167) = -A(1)+A(4)-A(5)-A(6)+A(7)-A(8)+2*A(9)+2*A(10)-A(11)+A(12)+A(16)+A(23)+A(24)-A(25)+A(26)+A(28)+A(31)-A(32)&
                &+A(33)+0.61*A(34)+A(35)+0.8*A(39)+2*A(40)+A(46)+A(48)+A(51)+A(53)+A(56)+A(58)+A(65)-A(69)+A(70)+A(71)+A(73)&
                &-A(79)+A(80)+A(81)+A(83)-A(90)+A(91)+A(92)+A(94)-A(102)+A(103)+A(104)+A(106)-A(115)-A(117)-A(120)+A(128)&
                &+0.338*A(177)+A(178)+0.187*A(192)+0.474*A(196)+0.474*A(200)+0.391*A(216)
  Vdot(168) = -A(46)-A(47)-A(48)-2*A(49)-2*A(50)-A(54)-A(64)+A(71)+A(73)-A(74)+2*A(78)-A(84)+A(88)-A(95)+A(99)-A(107)&
                &+A(111)+A(116)+A(131)+A(137)+0.65*A(142)+0.3*A(170)+A(185)+0.3*A(189)+0.25*A(193)+0.011*A(206)+0.076*A(211)&
                &+0.197*A(215)+0.03*A(216)+0.26*A(220)
END SUBROUTINE saprc99_mosaic_4bin_vbs9_Fun
SUBROUTINE saprc99_mosaic_4bin_vbs9_Jac_SP ( V, F, RCT, JVS )
  REAL(kind=dp) :: V(NVAR)
  REAL(kind=dp) :: F(NFIX)
  REAL(kind=dp) :: RCT(NREACT)
  REAL(kind=dp) :: JVS(LU_NONZERO)
  REAL(kind=dp) :: B(645)
  B(1) = RCT(1)
  B(2) = RCT(2)*F(2)
  B(4) = RCT(3)*V(156)
  B(5) = RCT(3)*V(152)
  B(6) = RCT(4)*V(160)*F(2)
  B(7) = RCT(4)*V(152)*F(2)
  B(9) = RCT(5)*V(167)
  B(10) = RCT(5)*V(152)
  B(11) = RCT(6)*V(167)
  B(12) = RCT(6)*V(152)
  B(13) = RCT(7)*V(160)
  B(14) = RCT(7)*V(156)
  B(15) = RCT(8)*V(167)
  B(16) = RCT(8)*V(156)
  B(17) = RCT(9)*V(160)
  B(18) = RCT(9)*V(157)
  B(19) = RCT(10)*2*V(160)*F(2)
  B(21) = RCT(11)*V(167)
  B(22) = RCT(11)*V(157)
  B(23) = RCT(12)
  B(24) = RCT(13)*F(1)
  B(26) = RCT(14)*V(167)
  B(27) = RCT(14)*V(157)
  B(28) = RCT(15)
  B(29) = RCT(16)
  B(30) = RCT(17)
  B(31) = RCT(18)
  B(32) = RCT(19)*F(1)
  B(34) = RCT(20)*F(2)
  B(36) = RCT(21)*V(165)
  B(37) = RCT(21)*V(160)
  B(38) = RCT(22)
  B(39) = RCT(23)
  B(40) = RCT(24)*V(165)
  B(41) = RCT(24)*V(111)
  B(42) = RCT(25)*V(167)
  B(43) = RCT(25)*V(165)
  B(44) = RCT(26)*V(165)
  B(45) = RCT(26)*V(157)
  B(46) = RCT(27)*V(165)
  B(47) = RCT(27)*V(134)
  B(48) = RCT(28)
  B(49) = RCT(29)*V(165)
  B(50) = RCT(29)*V(133)
  B(51) = RCT(30)*V(165)
  B(52) = RCT(30)*V(156)
  B(53) = RCT(31)*V(161)
  B(54) = RCT(31)*V(160)
  B(55) = RCT(32)*V(167)
  B(56) = RCT(32)*V(161)
  B(57) = RCT(33)
  B(58) = RCT(34)
  B(59) = RCT(35)*V(165)
  B(60) = RCT(35)*V(119)
  B(61) = RCT(36)*V(161)
  B(62) = RCT(36)*V(156)
  B(63) = RCT(37)*2*V(161)
  B(64) = RCT(38)*2*V(161)*F(1)
  B(66) = RCT(39)*V(161)
  B(67) = RCT(39)*V(157)
  B(68) = RCT(40)*2*V(157)
  B(69) = RCT(41)
  B(70) = RCT(42)*V(165)
  B(71) = RCT(42)*V(106)
  B(72) = RCT(43)*V(165)
  B(73) = RCT(43)*V(161)
  B(74) = RCT(44)*V(165)
  B(75) = RCT(44)*V(100)
  B(76) = RCT(45)*F(2)
  B(78) = RCT(46)*V(168)
  B(79) = RCT(46)*V(160)
  B(80) = RCT(47)*V(168)
  B(81) = RCT(47)*V(161)
  B(82) = RCT(48)*V(168)
  B(83) = RCT(48)*V(157)
  B(84) = RCT(49)*2*V(168)
  B(85) = RCT(50)*2*V(168)
  B(86) = RCT(51)*V(166)
  B(87) = RCT(51)*V(160)
  B(88) = RCT(52)*V(166)
  B(89) = RCT(52)*V(161)
  B(90) = RCT(53)*V(166)
  B(91) = RCT(53)*V(157)
  B(92) = RCT(54)*V(168)
  B(93) = RCT(54)*V(166)
  B(94) = RCT(55)*2*V(166)
  B(95) = RCT(56)*V(160)
  B(96) = RCT(56)*V(141)
  B(97) = RCT(57)*V(161)
  B(98) = RCT(57)*V(141)
  B(99) = RCT(58)*V(157)
  B(100) = RCT(58)*V(141)
  B(101) = RCT(59)*V(168)
  B(102) = RCT(59)*V(141)
  B(103) = RCT(60)*V(166)
  B(104) = RCT(60)*V(141)
  B(105) = RCT(61)*2*V(141)
  B(106) = RCT(62)*V(160)
  B(107) = RCT(62)*V(158)
  B(108) = RCT(63)*V(161)
  B(109) = RCT(63)*V(158)
  B(110) = RCT(64)*V(168)
  B(111) = RCT(64)*V(158)
  B(112) = RCT(65)*V(158)
  B(113) = RCT(65)*V(157)
  B(114) = RCT(66)*V(166)
  B(115) = RCT(66)*V(158)
  B(116) = RCT(67)*V(158)
  B(117) = RCT(67)*V(141)
  B(118) = RCT(68)*2*V(158)
  B(119) = RCT(69)*V(167)
  B(120) = RCT(69)*V(163)
  B(121) = RCT(70)
  B(122) = RCT(71)*V(163)
  B(123) = RCT(71)*V(160)
  B(124) = RCT(72)*V(163)
  B(125) = RCT(72)*V(161)
  B(126) = RCT(73)*V(163)
  B(127) = RCT(73)*V(157)
  B(128) = RCT(74)*V(168)
  B(129) = RCT(74)*V(163)
  B(130) = RCT(75)*V(166)
  B(131) = RCT(75)*V(163)
  B(132) = RCT(76)*V(163)
  B(133) = RCT(76)*V(141)
  B(134) = RCT(77)*V(163)
  B(135) = RCT(77)*V(158)
  B(136) = RCT(78)*2*V(163)
  B(137) = RCT(79)*V(167)
  B(138) = RCT(79)*V(162)
  B(139) = RCT(80)
  B(140) = RCT(81)*V(162)
  B(141) = RCT(81)*V(160)
  B(142) = RCT(82)*V(162)
  B(143) = RCT(82)*V(161)
  B(144) = RCT(83)*V(162)
  B(145) = RCT(83)*V(157)
  B(146) = RCT(84)*V(168)
  B(147) = RCT(84)*V(162)
  B(148) = RCT(85)*V(166)
  B(149) = RCT(85)*V(162)
  B(150) = RCT(86)*V(162)
  B(151) = RCT(86)*V(141)
  B(152) = RCT(87)*V(162)
  B(153) = RCT(87)*V(158)
  B(154) = RCT(88)*V(163)
  B(155) = RCT(88)*V(162)
  B(156) = RCT(89)*2*V(162)
  B(157) = RCT(90)*V(167)
  B(158) = RCT(90)*V(164)
  B(159) = RCT(91)
  B(160) = RCT(92)*V(164)
  B(161) = RCT(92)*V(160)
  B(162) = RCT(93)*V(164)
  B(163) = RCT(93)*V(161)
  B(164) = RCT(94)*V(164)
  B(165) = RCT(94)*V(157)
  B(166) = RCT(95)*V(168)
  B(167) = RCT(95)*V(164)
  B(168) = RCT(96)*V(166)
  B(169) = RCT(96)*V(164)
  B(170) = RCT(97)*V(164)
  B(171) = RCT(97)*V(141)
  B(172) = RCT(98)*V(164)
  B(173) = RCT(98)*V(158)
  B(174) = RCT(99)*V(164)
  B(175) = RCT(99)*V(163)
  B(176) = RCT(100)*V(164)
  B(177) = RCT(100)*V(162)
  B(178) = RCT(101)*2*V(164)
  B(179) = RCT(102)*V(167)
  B(180) = RCT(102)*V(159)
  B(181) = RCT(103)
  B(182) = RCT(104)*V(160)
  B(183) = RCT(104)*V(159)
  B(184) = RCT(105)*V(161)
  B(185) = RCT(105)*V(159)
  B(186) = RCT(106)*V(159)
  B(187) = RCT(106)*V(157)
  B(188) = RCT(107)*V(168)
  B(189) = RCT(107)*V(159)
  B(190) = RCT(108)*V(166)
  B(191) = RCT(108)*V(159)
  B(192) = RCT(109)*V(159)
  B(193) = RCT(109)*V(141)
  B(194) = RCT(110)*V(159)
  B(195) = RCT(110)*V(158)
  B(196) = RCT(111)*V(163)
  B(197) = RCT(111)*V(159)
  B(198) = RCT(112)*V(162)
  B(199) = RCT(112)*V(159)
  B(200) = RCT(113)*V(164)
  B(201) = RCT(113)*V(159)
  B(202) = RCT(114)*2*V(159)
  B(203) = RCT(115)*V(167)
  B(204) = RCT(115)*V(114)
  B(205) = RCT(116)
  B(206) = RCT(117)*V(167)
  B(207) = RCT(117)*V(138)
  B(208) = RCT(118)*V(161)
  B(209) = RCT(118)*V(138)
  B(210) = RCT(119)
  B(211) = RCT(120)*V(167)
  B(212) = RCT(120)*V(118)
  B(213) = RCT(121)*V(161)
  B(214) = RCT(121)*V(118)
  B(215) = RCT(122)
  B(216) = RCT(123)
  B(217) = RCT(124)
  B(218) = RCT(125)*V(165)
  B(219) = RCT(125)*V(150)
  B(220) = RCT(126)*V(161)
  B(221) = RCT(126)*V(150)
  B(222) = RCT(127)
  B(223) = RCT(128)*V(160)
  B(224) = RCT(128)*V(117)
  B(225) = RCT(129)*V(157)
  B(226) = RCT(129)*V(150)
  B(227) = RCT(130)*V(165)
  B(228) = RCT(130)*V(149)
  B(229) = RCT(131)
  B(230) = RCT(132)*V(157)
  B(231) = RCT(132)*V(149)
  B(232) = RCT(133)*V(165)
  B(233) = RCT(133)*V(153)
  B(234) = RCT(134)
  B(235) = RCT(135)*V(157)
  B(236) = RCT(135)*V(153)
  B(237) = RCT(136)*V(165)
  B(238) = RCT(136)*V(136)
  B(239) = RCT(137)
  B(240) = RCT(138)*V(165)
  B(241) = RCT(138)*V(154)
  B(242) = RCT(139)
  B(243) = RCT(140)*V(165)
  B(244) = RCT(140)*V(120)
  B(245) = RCT(141)*V(165)
  B(246) = RCT(141)*V(108)
  B(247) = RCT(142)*V(165)
  B(248) = RCT(142)*V(116)
  B(249) = RCT(143)
  B(250) = RCT(144)*V(165)
  B(251) = RCT(144)*V(129)
  B(252) = RCT(145)
  B(253) = RCT(146)
  B(254) = RCT(147)
  B(255) = RCT(148)*V(165)
  B(256) = RCT(148)*V(145)
  B(257) = RCT(149)*V(157)
  B(258) = RCT(149)*V(145)
  B(259) = RCT(150)
  B(260) = RCT(151)*V(165)
  B(261) = RCT(151)*V(132)
  B(262) = RCT(152)*V(157)
  B(263) = RCT(152)*V(132)
  B(264) = RCT(153)
  B(265) = RCT(154)*V(165)
  B(266) = RCT(154)*V(131)
  B(267) = RCT(155)*V(157)
  B(268) = RCT(155)*V(131)
  B(269) = RCT(156)*V(165)
  B(270) = RCT(156)*V(124)
  B(271) = RCT(157)*V(157)
  B(272) = RCT(157)*V(124)
  B(273) = RCT(158)*V(157)
  B(274) = RCT(158)*V(128)
  B(275) = RCT(159)*V(165)
  B(276) = RCT(159)*V(130)
  B(277) = RCT(160)
  B(278) = RCT(161)*V(157)
  B(279) = RCT(161)*V(130)
  B(280) = RCT(162)*V(165)
  B(281) = RCT(162)*V(142)
  B(282) = RCT(163)*V(156)
  B(283) = RCT(163)*V(142)
  B(284) = RCT(164)*V(157)
  B(285) = RCT(164)*V(142)
  B(286) = RCT(165)*V(152)
  B(287) = RCT(165)*V(142)
  B(288) = RCT(166)
  B(289) = RCT(167)*V(165)
  B(290) = RCT(167)*V(148)
  B(291) = RCT(168)*V(156)
  B(292) = RCT(168)*V(148)
  B(293) = RCT(169)*V(152)
  B(294) = RCT(169)*V(148)
  B(295) = RCT(170)
  B(296) = RCT(171)*V(165)
  B(297) = RCT(171)*V(146)
  B(298) = RCT(172)*V(156)
  B(299) = RCT(172)*V(146)
  B(300) = RCT(173)*V(157)
  B(301) = RCT(173)*V(146)
  B(302) = RCT(174)
  B(303) = RCT(175)*V(165)
  B(304) = RCT(175)*V(155)
  B(305) = RCT(176)
  B(306) = RCT(177)*V(165)
  B(307) = RCT(177)*V(151)
  B(308) = RCT(178)
  B(309) = RCT(179)*V(165)
  B(310) = RCT(179)*V(127)
  B(311) = RCT(180)*V(156)
  B(312) = RCT(180)*V(127)
  B(313) = RCT(181)*V(165)
  B(314) = RCT(181)*V(122)
  B(315) = RCT(182)
  B(316) = RCT(183)*V(165)
  B(317) = RCT(183)*V(123)
  B(318) = RCT(184)
  B(319) = RCT(185)*V(165)
  B(320) = RCT(185)*V(99)
  B(321) = RCT(186)*V(165)
  B(322) = RCT(186)*V(135)
  B(323) = RCT(187)*V(156)
  B(324) = RCT(187)*V(135)
  B(325) = RCT(188)*V(157)
  B(326) = RCT(188)*V(135)
  B(327) = RCT(189)*V(152)
  B(328) = RCT(189)*V(135)
  B(329) = RCT(190)*V(165)
  B(330) = RCT(190)*V(140)
  B(331) = RCT(191)*V(156)
  B(332) = RCT(191)*V(140)
  B(333) = RCT(192)*V(157)
  B(334) = RCT(192)*V(140)
  B(335) = RCT(193)*V(152)
  B(336) = RCT(193)*V(140)
  B(337) = RCT(194)*V(165)
  B(338) = RCT(194)*V(144)
  B(339) = RCT(195)*V(156)
  B(340) = RCT(195)*V(144)
  B(341) = RCT(196)*V(157)
  B(342) = RCT(196)*V(144)
  B(343) = RCT(197)*V(152)
  B(344) = RCT(197)*V(144)
  B(345) = RCT(198)*V(165)
  B(346) = RCT(198)*V(143)
  B(347) = RCT(199)*V(156)
  B(348) = RCT(199)*V(143)
  B(349) = RCT(200)*V(157)
  B(350) = RCT(200)*V(143)
  B(351) = RCT(201)*V(152)
  B(352) = RCT(201)*V(143)
  B(353) = RCT(202)*V(165)
  B(354) = RCT(202)*V(101)
  B(355) = RCT(203)*V(165)
  B(356) = RCT(203)*V(107)
  B(357) = RCT(204)*V(165)
  B(358) = RCT(204)*V(126)
  B(359) = RCT(205)*V(165)
  B(360) = RCT(205)*V(113)
  B(361) = RCT(206)*V(165)
  B(362) = RCT(206)*V(125)
  B(363) = RCT(207)*V(165)
  B(364) = RCT(207)*V(112)
  B(365) = RCT(208)*V(165)
  B(366) = RCT(208)*V(121)
  B(367) = RCT(209)*V(165)
  B(368) = RCT(209)*V(115)
  B(369) = RCT(210)*V(165)
  B(370) = RCT(210)*V(139)
  B(371) = RCT(211)*V(156)
  B(372) = RCT(211)*V(139)
  B(373) = RCT(212)*V(157)
  B(374) = RCT(212)*V(139)
  B(375) = RCT(213)*V(152)
  B(376) = RCT(213)*V(139)
  B(377) = RCT(214)*V(165)
  B(378) = RCT(214)*V(147)
  B(379) = RCT(215)*V(156)
  B(380) = RCT(215)*V(147)
  B(381) = RCT(216)*V(157)
  B(382) = RCT(216)*V(147)
  B(383) = RCT(217)*V(152)
  B(384) = RCT(217)*V(147)
  B(385) = RCT(218)*V(156)
  B(386) = RCT(218)*V(126)
  B(387) = RCT(219)*V(165)
  B(388) = RCT(219)*V(137)
  B(389) = RCT(220)*V(156)
  B(390) = RCT(220)*V(137)
  B(391) = RCT(221)*V(157)
  B(392) = RCT(221)*V(137)
  B(393) = RCT(222)*V(152)
  B(394) = RCT(222)*V(137)
  B(395) = RCT(223)
  B(396) = RCT(224)
  B(397) = RCT(225)
  B(398) = RCT(226)
  B(399) = RCT(227)
  B(400) = RCT(228)
  B(401) = RCT(229)
  B(402) = RCT(230)*V(165)
  B(403) = RCT(230)*V(125)
  B(404) = RCT(231)*V(165)
  B(405) = RCT(231)*V(125)
  B(406) = RCT(232)*V(165)
  B(407) = RCT(232)*V(112)
  B(408) = RCT(233)*V(165)
  B(409) = RCT(233)*V(112)
  B(410) = RCT(234)*V(165)
  B(411) = RCT(234)*V(139)
  B(412) = RCT(235)*V(165)
  B(413) = RCT(235)*V(139)
  B(414) = RCT(236)*V(165)
  B(415) = RCT(236)*V(139)
  B(416) = RCT(237)*V(165)
  B(417) = RCT(237)*V(139)
  B(418) = RCT(238)*V(165)
  B(419) = RCT(238)*V(139)
  B(420) = RCT(239)*V(165)
  B(421) = RCT(239)*V(139)
  B(422) = RCT(240)*V(165)
  B(423) = RCT(240)*V(139)
  B(424) = RCT(241)*V(165)
  B(425) = RCT(241)*V(139)
  B(426) = RCT(242)*V(165)
  B(427) = RCT(242)*V(147)
  B(428) = RCT(243)*V(165)
  B(429) = RCT(243)*V(147)
  B(430) = RCT(244)*V(165)
  B(431) = RCT(244)*V(147)
  B(432) = RCT(245)*V(165)
  B(433) = RCT(245)*V(147)
  B(434) = RCT(246)*V(165)
  B(435) = RCT(246)*V(147)
  B(436) = RCT(247)*V(165)
  B(437) = RCT(247)*V(147)
  B(438) = RCT(248)*V(165)
  B(439) = RCT(248)*V(147)
  B(440) = RCT(249)*V(165)
  B(441) = RCT(249)*V(147)
  B(442) = RCT(250)*V(165)
  B(443) = RCT(250)*V(121)
  B(444) = RCT(251)*V(165)
  B(445) = RCT(251)*V(121)
  B(446) = RCT(252)*V(165)
  B(447) = RCT(252)*V(121)
  B(448) = RCT(253)*V(165)
  B(449) = RCT(253)*V(121)
  B(450) = RCT(254)*V(165)
  B(451) = RCT(254)*V(121)
  B(452) = RCT(255)*V(165)
  B(453) = RCT(255)*V(121)
  B(454) = RCT(256)*V(165)
  B(455) = RCT(256)*V(121)
  B(456) = RCT(257)*V(165)
  B(457) = RCT(257)*V(121)
  B(458) = RCT(258)*V(165)
  B(459) = RCT(258)*V(115)
  B(460) = RCT(259)*V(165)
  B(461) = RCT(259)*V(115)
  B(462) = RCT(260)*V(165)
  B(463) = RCT(260)*V(115)
  B(464) = RCT(261)*V(165)
  B(465) = RCT(261)*V(115)
  B(466) = RCT(262)*V(165)
  B(467) = RCT(262)*V(115)
  B(468) = RCT(263)*V(165)
  B(469) = RCT(263)*V(115)
  B(470) = RCT(264)*V(165)
  B(471) = RCT(264)*V(115)
  B(472) = RCT(265)*V(165)
  B(473) = RCT(265)*V(115)
  B(474) = RCT(266)*V(165)
  B(475) = RCT(266)*V(140)
  B(476) = RCT(267)*V(165)
  B(477) = RCT(267)*V(140)
  B(478) = RCT(268)*V(165)
  B(479) = RCT(268)*V(140)
  B(480) = RCT(269)*V(165)
  B(481) = RCT(269)*V(140)
  B(482) = RCT(270)*V(165)
  B(483) = RCT(270)*V(140)
  B(484) = RCT(271)*V(165)
  B(485) = RCT(271)*V(140)
  B(486) = RCT(272)*V(165)
  B(487) = RCT(272)*V(144)
  B(488) = RCT(273)*V(165)
  B(489) = RCT(273)*V(144)
  B(490) = RCT(274)*V(165)
  B(491) = RCT(274)*V(144)
  B(492) = RCT(275)*V(165)
  B(493) = RCT(275)*V(144)
  B(494) = RCT(276)*V(165)
  B(495) = RCT(276)*V(144)
  B(496) = RCT(277)*V(165)
  B(497) = RCT(277)*V(144)
  B(498) = RCT(278)*V(165)
  B(499) = RCT(278)*V(144)
  B(500) = RCT(279)*V(165)
  B(501) = RCT(279)*V(144)
  B(502) = RCT(280)*V(165)
  B(503) = RCT(280)*V(143)
  B(504) = RCT(281)*V(165)
  B(505) = RCT(281)*V(143)
  B(506) = RCT(282)*V(165)
  B(507) = RCT(282)*V(143)
  B(508) = RCT(283)*V(165)
  B(509) = RCT(283)*V(143)
  B(510) = RCT(284)*V(165)
  B(511) = RCT(284)*V(3)
  B(512) = RCT(285)*V(165)
  B(513) = RCT(285)*V(4)
  B(514) = RCT(286)*V(165)
  B(515) = RCT(286)*V(74)
  B(516) = RCT(287)*V(165)
  B(517) = RCT(287)*V(34)
  B(518) = RCT(288)*V(165)
  B(519) = RCT(288)*V(5)
  B(520) = RCT(289)*V(165)
  B(521) = RCT(289)*V(6)
  B(522) = RCT(290)*V(165)
  B(523) = RCT(290)*V(91)
  B(524) = RCT(291)*V(165)
  B(525) = RCT(291)*V(50)
  B(526) = RCT(292)*V(165)
  B(527) = RCT(292)*V(66)
  B(528) = RCT(293)*V(165)
  B(529) = RCT(293)*V(35)
  B(530) = RCT(294)*V(165)
  B(531) = RCT(294)*V(75)
  B(532) = RCT(295)*V(165)
  B(533) = RCT(295)*V(36)
  B(534) = RCT(296)*V(165)
  B(535) = RCT(296)*V(67)
  B(536) = RCT(297)*V(165)
  B(537) = RCT(297)*V(37)
  B(538) = RCT(298)*V(165)
  B(539) = RCT(298)*V(76)
  B(540) = RCT(299)*V(165)
  B(541) = RCT(299)*V(38)
  B(542) = RCT(300)*V(165)
  B(543) = RCT(300)*V(68)
  B(544) = RCT(301)*V(165)
  B(545) = RCT(301)*V(39)
  B(546) = RCT(302)*V(165)
  B(547) = RCT(302)*V(77)
  B(548) = RCT(303)*V(165)
  B(549) = RCT(303)*V(40)
  B(550) = RCT(304)*V(165)
  B(551) = RCT(304)*V(69)
  B(552) = RCT(305)*V(165)
  B(553) = RCT(305)*V(41)
  B(554) = RCT(306)*V(165)
  B(555) = RCT(306)*V(78)
  B(556) = RCT(307)*V(165)
  B(557) = RCT(307)*V(42)
  B(558) = RCT(308)*V(165)
  B(559) = RCT(308)*V(70)
  B(560) = RCT(309)*V(165)
  B(561) = RCT(309)*V(43)
  B(562) = RCT(310)*V(165)
  B(563) = RCT(310)*V(79)
  B(564) = RCT(311)*V(165)
  B(565) = RCT(311)*V(44)
  B(566) = RCT(312)*V(165)
  B(567) = RCT(312)*V(71)
  B(568) = RCT(313)*V(165)
  B(569) = RCT(313)*V(45)
  B(570) = RCT(314)*V(165)
  B(571) = RCT(314)*V(80)
  B(572) = RCT(315)*V(165)
  B(573) = RCT(315)*V(46)
  B(574) = RCT(316)*V(165)
  B(575) = RCT(316)*V(72)
  B(576) = RCT(317)*V(165)
  B(577) = RCT(317)*V(47)
  B(578) = RCT(318)*V(165)
  B(579) = RCT(318)*V(81)
  B(580) = RCT(319)*V(165)
  B(581) = RCT(319)*V(48)
  B(582) = RCT(320)*V(165)
  B(583) = RCT(320)*V(73)
  B(584) = RCT(321)*V(165)
  B(585) = RCT(321)*V(49)
  B(586) = RCT(322)*V(165)
  B(587) = RCT(322)*V(82)
  B(588) = RCT(323)*V(165)
  B(589) = RCT(323)*V(51)
  B(590) = RCT(324)*V(165)
  B(591) = RCT(324)*V(92)
  B(592) = RCT(325)*V(165)
  B(593) = RCT(325)*V(52)
  B(594) = RCT(326)*V(165)
  B(595) = RCT(326)*V(83)
  B(596) = RCT(327)*V(165)
  B(597) = RCT(327)*V(53)
  B(598) = RCT(328)*V(165)
  B(599) = RCT(328)*V(93)
  B(600) = RCT(329)*V(165)
  B(601) = RCT(329)*V(54)
  B(602) = RCT(330)*V(165)
  B(603) = RCT(330)*V(84)
  B(604) = RCT(331)*V(165)
  B(605) = RCT(331)*V(55)
  B(606) = RCT(332)*V(165)
  B(607) = RCT(332)*V(94)
  B(608) = RCT(333)*V(165)
  B(609) = RCT(333)*V(56)
  B(610) = RCT(334)*V(165)
  B(611) = RCT(334)*V(85)
  B(612) = RCT(335)*V(165)
  B(613) = RCT(335)*V(57)
  B(614) = RCT(336)*V(165)
  B(615) = RCT(336)*V(95)
  B(616) = RCT(337)*V(165)
  B(617) = RCT(337)*V(58)
  B(618) = RCT(338)*V(165)
  B(619) = RCT(338)*V(86)
  B(620) = RCT(339)*V(165)
  B(621) = RCT(339)*V(59)
  B(622) = RCT(340)*V(165)
  B(623) = RCT(340)*V(96)
  B(624) = RCT(341)*V(165)
  B(625) = RCT(341)*V(60)
  B(626) = RCT(342)*V(165)
  B(627) = RCT(342)*V(87)
  B(628) = RCT(343)*V(165)
  B(629) = RCT(343)*V(61)
  B(630) = RCT(344)*V(165)
  B(631) = RCT(344)*V(97)
  B(632) = RCT(345)*V(165)
  B(633) = RCT(345)*V(62)
  B(634) = RCT(346)*V(165)
  B(635) = RCT(346)*V(88)
  B(636) = RCT(347)*V(165)
  B(637) = RCT(347)*V(63)
  B(638) = RCT(348)*V(165)
  B(639) = RCT(348)*V(98)
  B(640) = RCT(349)*V(165)
  B(641) = RCT(349)*V(64)
  B(642) = RCT(350)*V(165)
  B(643) = RCT(350)*V(89)
  B(644) = RCT(351)*V(165)
  B(645) = RCT(351)*V(65)
  JVS(1) = -B(398)
  JVS(2) = B(74)+B(395)
  JVS(3) = B(75)
  JVS(4) = -B(401)
  JVS(5) = 0.5*B(385)
  JVS(6) = 0.135*B(389)
  JVS(7) = 0.5*B(386)+0.135*B(390)
  JVS(8) = 0
  JVS(9) = 0
  JVS(10) = 0
  JVS(11) = 0
  JVS(12) = 0
  JVS(13) = B(223)
  JVS(14) = 0.297*B(357)
  JVS(15) = 0.37*B(323)
  JVS(16) = 0.185*B(389)
  JVS(17) = 0.185*B(371)
  JVS(18) = 0.204*B(331)
  JVS(19) = 0.333*B(282)
  JVS(20) = 0.103*B(347)
  JVS(21) = 0.103*B(339)
  JVS(22) = 0.1*B(298)
  JVS(23) = 0.073*B(379)
  JVS(24) = 0.351*B(291)
  JVS(25) = 0.333*B(283)+0.351*B(292)+0.1*B(299)+0.37*B(324)+0.204*B(332)+0.103*B(340)+0.103*B(348)+0.185*B(372)+0.073&
              &*B(380)+0.185*B(390)
  JVS(26) = B(224)
  JVS(27) = 0.297*B(358)
  JVS(28) = 0
  JVS(29) = 0.17*B(389)
  JVS(30) = 0.05*B(371)
  JVS(31) = 0.129*B(379)
  JVS(32) = 0.05*B(372)+0.129*B(380)+0.17*B(390)
  JVS(33) = B(134)
  JVS(34) = 0.25*B(124)
  JVS(35) = 0.25*B(125)+B(128)+B(130)+B(135)
  JVS(36) = B(131)
  JVS(37) = B(129)
  JVS(38) = 0
  JVS(39) = 0.119*B(371)
  JVS(40) = 0.15*B(331)
  JVS(41) = 0.189*B(347)
  JVS(42) = 0.189*B(339)
  JVS(43) = 0.372*B(298)
  JVS(44) = 0.247*B(379)
  JVS(45) = 0.372*B(299)+0.15*B(332)+0.189*B(340)+0.189*B(348)+0.119*B(372)+0.247*B(380)
  JVS(46) = B(152)+B(172)+2*B(194)
  JVS(47) = 0.25*B(184)+B(188)+B(190)+2*B(195)
  JVS(48) = 0.25*B(142)+0.25*B(162)+0.25*B(185)
  JVS(49) = 0.25*B(143)+B(146)+B(148)+B(153)
  JVS(50) = 0.25*B(163)+B(166)+B(168)+B(173)
  JVS(51) = B(149)+B(169)+B(191)
  JVS(52) = B(147)+B(167)+B(189)
  JVS(53) = 0
  JVS(54) = 0.75*B(124)
  JVS(55) = 0.75*B(125)
  JVS(56) = 0
  JVS(57) = 0.75*B(184)
  JVS(58) = 0.75*B(142)+0.75*B(162)+0.75*B(185)
  JVS(59) = 0.75*B(143)
  JVS(60) = 0.75*B(163)
  JVS(61) = 0
  JVS(62) = 2*B(211)
  JVS(63) = B(391)
  JVS(64) = B(392)
  JVS(65) = 2*B(212)
  JVS(66) = 0
  JVS(67) = 6*B(211)
  JVS(68) = 7*B(277)
  JVS(69) = 0.048*B(387)+0.07*B(389)+2.693*B(391)+0.55*B(393)
  JVS(70) = 0.55*B(394)
  JVS(71) = 0.07*B(390)
  JVS(72) = 2.693*B(392)
  JVS(73) = 0.048*B(388)
  JVS(74) = 6*B(212)
  JVS(75) = 0
  JVS(76) = B(95)+B(99)
  JVS(77) = B(82)+B(90)+B(100)+B(112)
  JVS(78) = B(106)+B(113)
  JVS(79) = B(78)+B(86)+B(96)+B(107)
  JVS(80) = B(87)+B(91)
  JVS(81) = B(79)+B(83)
  JVS(82) = 0
  JVS(83) = B(97)+B(101)+B(103)+B(105)+B(116)
  JVS(84) = B(108)+B(110)+B(114)+B(117)+B(118)
  JVS(85) = B(80)+B(88)+B(98)+B(109)
  JVS(86) = B(89)+B(92)+B(94)+B(104)+B(115)
  JVS(87) = B(81)+B(84)+B(85)+B(93)+B(102)+B(111)
  JVS(88) = 0
  JVS(89) = B(466)
  JVS(90) = B(450)
  JVS(91) = B(410)
  JVS(92) = B(426)
  JVS(93) = B(411)+B(427)+B(451)+B(467)
  JVS(94) = 0
  JVS(95) = B(406)
  JVS(96) = B(468)
  JVS(97) = B(452)
  JVS(98) = B(402)
  JVS(99) = B(412)
  JVS(100) = B(428)
  JVS(101) = B(403)+B(407)+B(413)+B(429)+B(453)+B(469)
  JVS(102) = 0
  JVS(103) = B(470)
  JVS(104) = B(454)
  JVS(105) = B(414)
  JVS(106) = B(430)
  JVS(107) = B(415)+B(431)+B(455)+B(471)
  JVS(108) = 0
  JVS(109) = B(472)
  JVS(110) = B(456)
  JVS(111) = B(416)
  JVS(112) = B(432)
  JVS(113) = B(417)+B(433)+B(457)+B(473)
  JVS(114) = 0
  JVS(115) = B(458)
  JVS(116) = B(442)
  JVS(117) = B(418)
  JVS(118) = B(434)
  JVS(119) = B(419)+B(435)+B(443)+B(459)
  JVS(120) = 0
  JVS(121) = B(408)
  JVS(122) = B(460)
  JVS(123) = B(444)
  JVS(124) = B(404)
  JVS(125) = B(420)
  JVS(126) = B(436)
  JVS(127) = B(405)+B(409)+B(421)+B(437)+B(445)+B(461)
  JVS(128) = 0
  JVS(129) = B(462)
  JVS(130) = B(446)
  JVS(131) = B(422)
  JVS(132) = B(438)
  JVS(133) = B(423)+B(439)+B(447)+B(463)
  JVS(134) = 0
  JVS(135) = B(464)
  JVS(136) = B(448)
  JVS(137) = B(424)
  JVS(138) = B(440)
  JVS(139) = B(425)+B(441)+B(449)+B(465)
  JVS(140) = 0
  JVS(141) = B(474)
  JVS(142) = B(494)
  JVS(143) = B(475)+B(495)
  JVS(144) = 0
  JVS(145) = B(476)
  JVS(146) = B(496)
  JVS(147) = B(477)+B(497)
  JVS(148) = 0
  JVS(149) = B(478)
  JVS(150) = B(498)
  JVS(151) = B(479)+B(499)
  JVS(152) = 0
  JVS(153) = B(500)
  JVS(154) = B(501)
  JVS(155) = 0
  JVS(156) = B(480)
  JVS(157) = B(502)
  JVS(158) = B(486)
  JVS(159) = B(481)+B(487)+B(503)
  JVS(160) = 0
  JVS(161) = B(482)
  JVS(162) = B(504)
  JVS(163) = B(488)
  JVS(164) = B(483)+B(489)+B(505)
  JVS(165) = 0
  JVS(166) = B(484)
  JVS(167) = B(506)
  JVS(168) = B(490)
  JVS(169) = B(485)+B(491)+B(507)
  JVS(170) = 0
  JVS(171) = B(508)
  JVS(172) = B(492)
  JVS(173) = B(493)+B(509)
  JVS(174) = 0
  JVS(175) = B(319)
  JVS(176) = B(74)
  JVS(177) = B(353)
  JVS(178) = B(355)
  JVS(179) = B(245)
  JVS(180) = B(40)
  JVS(181) = B(363)
  JVS(182) = B(359)
  JVS(183) = B(367)
  JVS(184) = B(247)
  JVS(185) = B(243)
  JVS(186) = B(365)
  JVS(187) = B(313)
  JVS(188) = B(316)
  JVS(189) = B(269)
  JVS(190) = B(361)
  JVS(191) = B(357)
  JVS(192) = B(309)
  JVS(193) = B(250)
  JVS(194) = B(275)
  JVS(195) = B(265)
  JVS(196) = B(260)
  JVS(197) = B(49)
  JVS(198) = B(46)
  JVS(199) = B(321)
  JVS(200) = B(237)
  JVS(201) = B(387)
  JVS(202) = B(369)
  JVS(203) = B(329)
  JVS(204) = B(280)
  JVS(205) = B(345)
  JVS(206) = B(337)
  JVS(207) = B(255)
  JVS(208) = B(296)
  JVS(209) = B(377)
  JVS(210) = B(289)
  JVS(211) = B(227)
  JVS(212) = B(218)
  JVS(213) = B(306)
  JVS(214) = B(232)
  JVS(215) = B(240)
  JVS(216) = B(303)
  JVS(217) = B(51)
  JVS(218) = B(44)
  JVS(219) = B(36)
  JVS(220) = B(72)
  JVS(221) = B(37)+B(41)+B(42)+B(45)+B(47)+B(50)+B(52)+B(73)+B(75)+B(76)+B(219)+B(228)+B(233)+B(238)+B(241)+B(244)&
               &+B(246)+B(248)+B(251)+B(256)+B(261)+B(266)+B(270)+B(276)+B(281)+B(290)+B(297)+B(304)+B(307)+B(310)+B(314)&
               &+B(317)+B(320)+B(322)+B(330)+B(338)+B(346)+B(354)+B(356)+B(358)+B(360)+B(362)+B(364)+B(366)+B(368)+B(370)&
               &+B(378)+B(388)
  JVS(222) = B(43)
  JVS(223) = B(510)
  JVS(224) = B(512)
  JVS(225) = B(518)
  JVS(226) = B(520)
  JVS(227) = 0
  JVS(228) = B(516)
  JVS(229) = B(528)
  JVS(230) = B(532)
  JVS(231) = B(536)
  JVS(232) = B(540)
  JVS(233) = B(544)
  JVS(234) = B(548)
  JVS(235) = B(552)
  JVS(236) = B(556)
  JVS(237) = B(560)
  JVS(238) = B(564)
  JVS(239) = B(568)
  JVS(240) = B(572)
  JVS(241) = B(576)
  JVS(242) = B(580)
  JVS(243) = B(584)
  JVS(244) = B(524)
  JVS(245) = B(588)
  JVS(246) = B(592)
  JVS(247) = B(596)
  JVS(248) = B(600)
  JVS(249) = B(604)
  JVS(250) = B(608)
  JVS(251) = B(612)
  JVS(252) = B(616)
  JVS(253) = B(620)
  JVS(254) = B(624)
  JVS(255) = B(628)
  JVS(256) = B(632)
  JVS(257) = B(636)
  JVS(258) = B(640)
  JVS(259) = B(644)
  JVS(260) = B(526)
  JVS(261) = B(534)
  JVS(262) = B(542)
  JVS(263) = B(550)
  JVS(264) = B(558)
  JVS(265) = B(566)
  JVS(266) = B(574)
  JVS(267) = B(582)
  JVS(268) = B(514)
  JVS(269) = B(530)
  JVS(270) = B(538)
  JVS(271) = B(546)
  JVS(272) = B(554)
  JVS(273) = B(562)
  JVS(274) = B(570)
  JVS(275) = B(578)
  JVS(276) = B(586)
  JVS(277) = B(594)
  JVS(278) = B(602)
  JVS(279) = B(610)
  JVS(280) = B(618)
  JVS(281) = B(626)
  JVS(282) = B(634)
  JVS(283) = B(642)
  JVS(284) = B(522)
  JVS(285) = B(590)
  JVS(286) = B(598)
  JVS(287) = B(606)
  JVS(288) = B(614)
  JVS(289) = B(622)
  JVS(290) = B(630)
  JVS(291) = B(638)
  JVS(292) = B(511)+B(513)+B(515)+B(517)+B(519)+B(521)+B(523)+B(525)+B(527)+B(529)+B(531)+B(533)+B(535)+B(537)+B(539)&
               &+B(541)+B(543)+B(545)+B(547)+B(549)+B(551)+B(553)+B(555)+B(557)+B(559)+B(561)+B(563)+B(565)+B(567)+B(569)&
               &+B(571)+B(573)+B(575)+B(577)+B(579)+B(581)+B(583)+B(585)+B(587)+B(589)+B(591)+B(593)+B(595)+B(597)+B(599)&
               &+B(601)+B(603)+B(605)+B(607)+B(609)+B(611)+B(613)+B(615)+B(617)+B(619)+B(621)+B(623)+B(625)+B(627)+B(629)&
               &+B(631)+B(633)+B(635)+B(637)+B(639)+B(641)+B(643)+B(645)
  JVS(293) = 0
  JVS(294) = B(528)
  JVS(295) = B(532)
  JVS(296) = 0.15*B(526)
  JVS(297) = 0.15*B(530)
  JVS(298) = 0.15*B(527)+B(529)+0.15*B(531)+B(533)
  JVS(299) = -B(528)
  JVS(300) = -B(529)
  JVS(301) = -B(532)
  JVS(302) = B(536)
  JVS(303) = B(540)
  JVS(304) = 0.15*B(534)
  JVS(305) = 0.15*B(538)
  JVS(306) = -B(533)+0.15*B(535)+B(537)+0.15*B(539)+B(541)
  JVS(307) = -B(536)
  JVS(308) = -B(537)
  JVS(309) = -B(540)
  JVS(310) = B(544)
  JVS(311) = B(548)
  JVS(312) = 0.15*B(542)
  JVS(313) = 0.15*B(546)
  JVS(314) = -B(541)+0.15*B(543)+B(545)+0.15*B(547)+B(549)
  JVS(315) = -B(544)
  JVS(316) = -B(545)
  JVS(317) = -B(548)
  JVS(318) = B(552)
  JVS(319) = B(556)
  JVS(320) = 0.15*B(550)
  JVS(321) = 0.15*B(554)
  JVS(322) = -B(549)+0.15*B(551)+B(553)+0.15*B(555)+B(557)
  JVS(323) = -B(552)
  JVS(324) = -B(553)
  JVS(325) = -B(556)
  JVS(326) = B(560)
  JVS(327) = B(564)
  JVS(328) = 0.15*B(558)
  JVS(329) = 0.15*B(562)
  JVS(330) = -B(557)+0.15*B(559)+B(561)+0.15*B(563)+B(565)
  JVS(331) = -B(560)
  JVS(332) = -B(561)
  JVS(333) = -B(564)
  JVS(334) = B(568)
  JVS(335) = B(572)
  JVS(336) = 0.15*B(566)
  JVS(337) = 0.15*B(570)
  JVS(338) = -B(565)+0.15*B(567)+B(569)+0.15*B(571)+B(573)
  JVS(339) = -B(568)
  JVS(340) = -B(569)
  JVS(341) = -B(572)
  JVS(342) = B(576)
  JVS(343) = B(580)
  JVS(344) = 0.15*B(574)
  JVS(345) = 0.15*B(578)
  JVS(346) = -B(573)+0.15*B(575)+B(577)+0.15*B(579)+B(581)
  JVS(347) = -B(576)
  JVS(348) = -B(577)
  JVS(349) = -B(580)
  JVS(350) = B(584)
  JVS(351) = 0.15*B(582)
  JVS(352) = -B(581)+0.15*B(583)+B(585)
  JVS(353) = -B(584)
  JVS(354) = -B(585)
  JVS(355) = 0
  JVS(356) = B(588)
  JVS(357) = B(592)
  JVS(358) = 0.15*B(586)
  JVS(359) = 0.15*B(590)
  JVS(360) = 0.15*B(587)+B(589)+0.15*B(591)+B(593)
  JVS(361) = -B(588)
  JVS(362) = -B(589)
  JVS(363) = -B(592)
  JVS(364) = B(596)
  JVS(365) = B(600)
  JVS(366) = 0.15*B(594)
  JVS(367) = 0.15*B(598)
  JVS(368) = -B(593)+0.15*B(595)+B(597)+0.15*B(599)+B(601)
  JVS(369) = -B(596)
  JVS(370) = -B(597)
  JVS(371) = -B(600)
  JVS(372) = B(604)
  JVS(373) = B(608)
  JVS(374) = 0.15*B(602)
  JVS(375) = 0.15*B(606)
  JVS(376) = -B(601)+0.15*B(603)+B(605)+0.15*B(607)+B(609)
  JVS(377) = -B(604)
  JVS(378) = -B(605)
  JVS(379) = -B(608)
  JVS(380) = B(612)
  JVS(381) = B(616)
  JVS(382) = 0.15*B(610)
  JVS(383) = 0.15*B(614)
  JVS(384) = -B(609)+0.15*B(611)+B(613)+0.15*B(615)+B(617)
  JVS(385) = -B(612)
  JVS(386) = -B(613)
  JVS(387) = -B(616)
  JVS(388) = B(620)
  JVS(389) = B(624)
  JVS(390) = 0.15*B(618)
  JVS(391) = 0.15*B(622)
  JVS(392) = -B(617)+0.15*B(619)+B(621)+0.15*B(623)+B(625)
  JVS(393) = -B(620)
  JVS(394) = -B(621)
  JVS(395) = -B(624)
  JVS(396) = B(628)
  JVS(397) = B(632)
  JVS(398) = 0.15*B(626)
  JVS(399) = 0.15*B(630)
  JVS(400) = -B(625)+0.15*B(627)+B(629)+0.15*B(631)+B(633)
  JVS(401) = -B(628)
  JVS(402) = -B(629)
  JVS(403) = -B(632)
  JVS(404) = B(636)
  JVS(405) = B(640)
  JVS(406) = 0.15*B(634)
  JVS(407) = 0.15*B(638)
  JVS(408) = -B(633)+0.15*B(635)+B(637)+0.15*B(639)+B(641)
  JVS(409) = -B(636)
  JVS(410) = -B(637)
  JVS(411) = -B(640)
  JVS(412) = B(644)
  JVS(413) = 0.15*B(642)
  JVS(414) = -B(641)+0.15*B(643)+B(645)
  JVS(415) = -B(644)
  JVS(416) = -B(645)
  JVS(417) = -B(526)
  JVS(418) = -B(527)
  JVS(419) = -B(534)
  JVS(420) = -B(535)
  JVS(421) = -B(542)
  JVS(422) = -B(543)
  JVS(423) = -B(550)
  JVS(424) = -B(551)
  JVS(425) = -B(558)
  JVS(426) = -B(559)
  JVS(427) = -B(566)
  JVS(428) = -B(567)
  JVS(429) = -B(574)
  JVS(430) = -B(575)
  JVS(431) = -B(582)
  JVS(432) = -B(583)
  JVS(433) = B(526)
  JVS(434) = 0
  JVS(435) = B(530)
  JVS(436) = B(527)+B(531)
  JVS(437) = B(534)
  JVS(438) = -B(530)
  JVS(439) = B(538)
  JVS(440) = -B(531)+B(535)+B(539)
  JVS(441) = B(542)
  JVS(442) = -B(538)
  JVS(443) = B(546)
  JVS(444) = -B(539)+B(543)+B(547)
  JVS(445) = B(550)
  JVS(446) = -B(546)
  JVS(447) = B(554)
  JVS(448) = -B(547)+B(551)+B(555)
  JVS(449) = B(558)
  JVS(450) = -B(554)
  JVS(451) = B(562)
  JVS(452) = -B(555)+B(559)+B(563)
  JVS(453) = B(566)
  JVS(454) = -B(562)
  JVS(455) = B(570)
  JVS(456) = -B(563)+B(567)+B(571)
  JVS(457) = B(574)
  JVS(458) = -B(570)
  JVS(459) = B(578)
  JVS(460) = -B(571)+B(575)+B(579)
  JVS(461) = B(582)
  JVS(462) = -B(578)
  JVS(463) = -B(579)+B(583)
  JVS(464) = -B(586)
  JVS(465) = -B(587)
  JVS(466) = -B(594)
  JVS(467) = -B(595)
  JVS(468) = -B(602)
  JVS(469) = -B(603)
  JVS(470) = -B(610)
  JVS(471) = -B(611)
  JVS(472) = -B(618)
  JVS(473) = -B(619)
  JVS(474) = -B(626)
  JVS(475) = -B(627)
  JVS(476) = -B(634)
  JVS(477) = -B(635)
  JVS(478) = -B(642)
  JVS(479) = -B(643)
  JVS(480) = -B(32)-B(34)
  JVS(481) = B(31)
  JVS(482) = B(586)
  JVS(483) = 0
  JVS(484) = B(590)
  JVS(485) = B(587)+B(591)
  JVS(486) = B(594)
  JVS(487) = -B(590)
  JVS(488) = B(598)
  JVS(489) = -B(591)+B(595)+B(599)
  JVS(490) = B(602)
  JVS(491) = -B(598)
  JVS(492) = B(606)
  JVS(493) = -B(599)+B(603)+B(607)
  JVS(494) = B(610)
  JVS(495) = -B(606)
  JVS(496) = B(614)
  JVS(497) = -B(607)+B(611)+B(615)
  JVS(498) = B(618)
  JVS(499) = -B(614)
  JVS(500) = B(622)
  JVS(501) = -B(615)+B(619)+B(623)
  JVS(502) = B(626)
  JVS(503) = -B(622)
  JVS(504) = B(630)
  JVS(505) = -B(623)+B(627)+B(631)
  JVS(506) = B(634)
  JVS(507) = -B(630)
  JVS(508) = B(638)
  JVS(509) = -B(631)+B(635)+B(639)
  JVS(510) = B(642)
  JVS(511) = -B(638)
  JVS(512) = -B(639)+B(643)
  JVS(513) = -B(319)
  JVS(514) = -B(320)
  JVS(515) = -B(74)-B(395)-B(397)
  JVS(516) = -B(75)
  JVS(517) = -B(353)
  JVS(518) = -B(354)
  JVS(519) = -B(139)
  JVS(520) = B(137)
  JVS(521) = B(138)
  JVS(522) = -B(181)
  JVS(523) = B(179)
  JVS(524) = B(180)
  JVS(525) = -B(121)
  JVS(526) = B(119)
  JVS(527) = B(120)
  JVS(528) = -B(159)
  JVS(529) = B(157)
  JVS(530) = B(158)
  JVS(531) = -B(69)-B(70)-B(400)
  JVS(532) = B(63)+B(64)
  JVS(533) = -B(71)
  JVS(534) = -B(355)
  JVS(535) = -B(356)
  JVS(536) = -B(245)
  JVS(537) = -B(246)
  JVS(538) = -B(264)
  JVS(539) = 0.087*B(367)
  JVS(540) = 0.031*B(347)
  JVS(541) = 0.031*B(339)
  JVS(542) = 0.031*B(340)+0.031*B(348)
  JVS(543) = 0.087*B(368)
  JVS(544) = -B(23)-B(24)
  JVS(545) = B(21)
  JVS(546) = B(22)
  JVS(547) = -B(38)-B(39)-B(40)
  JVS(548) = B(36)
  JVS(549) = B(37)-B(41)
  JVS(550) = -B(363)
  JVS(551) = -B(364)
  JVS(552) = -B(359)
  JVS(553) = -B(360)
  JVS(554) = 0.236*B(359)
  JVS(555) = -B(203)-B(205)
  JVS(556) = 0.236*B(360)
  JVS(557) = -B(204)
  JVS(558) = -B(367)
  JVS(559) = -B(368)
  JVS(560) = -B(247)-B(249)
  JVS(561) = B(80)
  JVS(562) = -B(248)
  JVS(563) = B(81)
  JVS(564) = -B(222)-B(223)
  JVS(565) = B(220)
  JVS(566) = -B(224)
  JVS(567) = B(221)
  JVS(568) = -B(211)-B(213)-B(215)
  JVS(569) = B(273)
  JVS(570) = B(274)
  JVS(571) = -B(214)
  JVS(572) = -B(212)
  JVS(573) = -B(57)-B(58)-B(59)
  JVS(574) = B(55)
  JVS(575) = -B(60)
  JVS(576) = B(56)
  JVS(577) = -B(243)
  JVS(578) = 0.25*B(110)
  JVS(579) = -B(244)
  JVS(580) = 0.25*B(92)
  JVS(581) = B(84)+0.25*B(93)+0.25*B(111)
  JVS(582) = -B(365)
  JVS(583) = -B(366)
  JVS(584) = 0.099*B(367)
  JVS(585) = 0.108*B(365)
  JVS(586) = -B(313)-B(315)
  JVS(587) = -B(314)+0.108*B(366)+0.099*B(368)
  JVS(588) = 0.093*B(367)
  JVS(589) = 0.051*B(365)
  JVS(590) = -B(316)-B(318)
  JVS(591) = -B(317)+0.051*B(366)+0.093*B(368)
  JVS(592) = 0.187*B(367)
  JVS(593) = 0.207*B(365)
  JVS(594) = -B(269)-B(271)
  JVS(595) = -B(272)
  JVS(596) = -B(270)+0.207*B(366)+0.187*B(368)
  JVS(597) = -B(361)
  JVS(598) = -B(362)
  JVS(599) = -B(357)-B(385)
  JVS(600) = -B(386)
  JVS(601) = -B(358)
  JVS(602) = 0.561*B(367)
  JVS(603) = 0.491*B(365)
  JVS(604) = -B(309)-B(311)
  JVS(605) = -B(312)
  JVS(606) = -B(310)+0.491*B(366)+0.561*B(368)
  JVS(607) = B(213)+B(215)
  JVS(608) = -B(273)
  JVS(609) = B(206)
  JVS(610) = -B(274)
  JVS(611) = B(214)
  JVS(612) = B(207)
  JVS(613) = -B(250)-B(252)
  JVS(614) = B(108)
  JVS(615) = B(88)+B(109)
  JVS(616) = -B(251)
  JVS(617) = B(89)
  JVS(618) = 0.05*B(367)
  JVS(619) = 0.059*B(365)
  JVS(620) = -B(275)-B(277)-B(278)
  JVS(621) = 0.061*B(377)+0.042*B(379)+0.015*B(381)
  JVS(622) = 0.042*B(380)
  JVS(623) = -B(279)+0.015*B(382)
  JVS(624) = -B(276)+0.059*B(366)+0.05*B(368)+0.061*B(378)
  JVS(625) = 0.017*B(365)
  JVS(626) = -B(265)-B(267)
  JVS(627) = B(208)+B(210)
  JVS(628) = -B(268)
  JVS(629) = B(209)
  JVS(630) = -B(266)+0.017*B(366)
  JVS(631) = 0.287*B(367)
  JVS(632) = 0.119*B(365)
  JVS(633) = 0.5*B(315)
  JVS(634) = 0.5*B(318)
  JVS(635) = 0.23*B(269)
  JVS(636) = -B(259)-B(260)-B(262)
  JVS(637) = 0.084*B(280)+0.9*B(282)
  JVS(638) = 0.174*B(296)+0.742*B(298)+0.008*B(300)
  JVS(639) = 0.3*B(289)+0.95*B(291)
  JVS(640) = 0.9*B(283)+0.95*B(292)+0.742*B(299)
  JVS(641) = -B(263)+0.008*B(301)
  JVS(642) = -B(261)+0.23*B(270)+0.084*B(281)+0.3*B(290)+0.174*B(297)+0.119*B(366)+0.287*B(368)
  JVS(643) = B(315)
  JVS(644) = B(318)
  JVS(645) = 0.002*B(361)
  JVS(646) = 0.393*B(357)+1.5*B(385)
  JVS(647) = B(309)+1.5*B(311)
  JVS(648) = B(259)+B(260)+B(262)
  JVS(649) = -B(49)
  JVS(650) = 0.5*B(323)+0.491*B(327)
  JVS(651) = 0.51*B(389)
  JVS(652) = 0.345*B(371)
  JVS(653) = 0.275*B(331)
  JVS(654) = 0.416*B(280)+0.45*B(282)+0.5*B(284)+0.67*B(288)
  JVS(655) = 0.157*B(347)
  JVS(656) = 0.157*B(339)
  JVS(657) = 2*B(253)+B(254)+1.26*B(255)+1.26*B(257)
  JVS(658) = 0.336*B(296)+0.498*B(298)+0.572*B(300)+1.233*B(302)
  JVS(659) = 0.265*B(379)+0.012*B(383)
  JVS(660) = 0.475*B(291)+0.7*B(295)
  JVS(661) = B(229)
  JVS(662) = B(216)+B(217)+B(218)+B(225)
  JVS(663) = 0.491*B(328)+0.012*B(384)
  JVS(664) = 0.034*B(232)+B(234)
  JVS(665) = 0.45*B(283)+0.475*B(292)+0.498*B(299)+1.5*B(312)+0.5*B(324)+0.275*B(332)+0.157*B(340)+0.157*B(348)+0.345&
               &*B(372)+0.265*B(380)+1.5*B(386)+0.51*B(390)
  JVS(666) = B(226)+1.26*B(258)+B(263)+0.5*B(285)+0.572*B(301)
  JVS(667) = -B(50)+B(219)+0.034*B(233)+1.26*B(256)+B(261)+0.416*B(281)+0.336*B(297)+B(310)+0.393*B(358)+0.002*B(362)
  JVS(668) = 2*B(24)
  JVS(669) = B(271)
  JVS(670) = B(273)
  JVS(671) = B(278)
  JVS(672) = B(267)
  JVS(673) = B(262)
  JVS(674) = -B(46)-B(48)-B(399)
  JVS(675) = 0
  JVS(676) = 0.5*B(284)
  JVS(677) = B(257)
  JVS(678) = 0.15*B(300)
  JVS(679) = 0
  JVS(680) = 0
  JVS(681) = B(230)
  JVS(682) = B(225)
  JVS(683) = B(235)
  JVS(684) = 0
  JVS(685) = 0.2*B(66)+B(226)+B(231)+B(236)+B(258)+B(263)+B(268)+B(272)+B(274)+B(279)+0.5*B(285)+0.15*B(301)
  JVS(686) = 0.2*B(67)
  JVS(687) = B(42)-B(47)
  JVS(688) = B(43)
  JVS(689) = -B(321)-B(323)-B(325)-B(327)
  JVS(690) = -B(328)
  JVS(691) = -B(324)
  JVS(692) = -B(326)
  JVS(693) = -B(322)
  JVS(694) = 0.704*B(355)
  JVS(695) = 0.072*B(363)
  JVS(696) = 0.024*B(359)
  JVS(697) = B(205)
  JVS(698) = 0.452*B(361)
  JVS(699) = -B(237)-B(239)
  JVS(700) = 0.005*B(369)+0.001*B(371)+0.024*B(373)
  JVS(701) = 0.13*B(347)
  JVS(702) = 0.13*B(339)
  JVS(703) = 0.127*B(377)+0.045*B(379)+0.102*B(381)
  JVS(704) = 0.006*B(306)+0.02*B(308)
  JVS(705) = 0.13*B(340)+0.13*B(348)+0.001*B(372)+0.045*B(380)
  JVS(706) = 0.024*B(374)+0.102*B(382)
  JVS(707) = -B(238)+0.006*B(307)+0.704*B(356)+0.024*B(360)+0.452*B(362)+0.072*B(364)+0.005*B(370)+0.127*B(378)
  JVS(708) = 0
  JVS(709) = -B(387)-B(389)-B(391)-B(393)
  JVS(710) = -B(394)
  JVS(711) = -B(390)
  JVS(712) = -B(392)
  JVS(713) = -B(388)
  JVS(714) = 0.24*B(269)+B(271)
  JVS(715) = 0.24*B(265)+B(267)
  JVS(716) = -B(206)-B(208)-B(210)
  JVS(717) = B(164)+B(268)+B(272)
  JVS(718) = B(200)
  JVS(719) = B(160)
  JVS(720) = -B(209)
  JVS(721) = B(176)
  JVS(722) = B(174)
  JVS(723) = B(161)+B(165)+B(175)+B(177)+2*B(178)+B(201)
  JVS(724) = 0.24*B(266)+0.24*B(270)
  JVS(725) = -B(207)
  JVS(726) = -B(369)-B(371)-B(373)-B(375)
  JVS(727) = -B(376)
  JVS(728) = -B(372)
  JVS(729) = -B(374)
  JVS(730) = -B(370)
  JVS(731) = -B(329)-B(331)-B(333)-B(335)
  JVS(732) = -B(336)
  JVS(733) = -B(332)
  JVS(734) = -B(334)
  JVS(735) = -B(330)
  JVS(736) = 0.948*B(363)
  JVS(737) = 0.559*B(359)
  JVS(738) = B(313)+B(315)
  JVS(739) = B(316)+B(318)
  JVS(740) = 0.936*B(361)
  JVS(741) = B(237)
  JVS(742) = 0.205*B(369)+0.488*B(373)
  JVS(743) = 0.079*B(329)+0.126*B(331)+0.187*B(333)+0.24*B(335)
  JVS(744) = -B(95)-B(97)-B(99)-B(101)-B(103)-B(116)-B(132)-B(150)-B(170)-B(192)
  JVS(745) = 0.5*B(345)+0.729*B(347)+0.75*B(349)
  JVS(746) = 0.5*B(337)+0.729*B(339)+0.75*B(341)
  JVS(747) = 0.001*B(377)+0.137*B(379)+0.711*B(381)
  JVS(748) = 0.675*B(289)
  JVS(749) = 0.596*B(306)+0.152*B(308)
  JVS(750) = 0.24*B(336)
  JVS(751) = 0.616*B(240)
  JVS(752) = 0.515*B(305)
  JVS(753) = 0.126*B(332)+0.729*B(340)+0.729*B(348)+0.137*B(380)
  JVS(754) = -B(100)+B(164)+0.187*B(334)+0.75*B(342)+0.75*B(350)+0.488*B(374)+0.711*B(382)
  JVS(755) = -B(117)
  JVS(756) = -B(193)+B(200)
  JVS(757) = -B(96)+B(160)
  JVS(758) = -B(98)
  JVS(759) = -B(151)+B(176)
  JVS(760) = -B(133)+B(174)
  JVS(761) = B(161)+B(165)-B(171)+B(175)+B(177)+2*B(178)+B(201)
  JVS(762) = B(238)+0.616*B(241)+0.675*B(290)+0.596*B(307)+B(314)+B(317)+0.079*B(330)+0.5*B(338)+0.5*B(346)+0.559*B(360)&
               &+0.936*B(362)+0.948*B(364)+0.205*B(370)+0.001*B(378)
  JVS(763) = -B(104)
  JVS(764) = 0
  JVS(765) = -B(102)
  JVS(766) = 0.23*B(329)+0.39*B(331)
  JVS(767) = -B(280)-B(282)-B(284)-B(286)-B(288)
  JVS(768) = 0.025*B(377)+0.026*B(379)+0.012*B(383)
  JVS(769) = -B(287)+0.012*B(384)
  JVS(770) = -B(283)+0.39*B(332)+0.026*B(380)
  JVS(771) = -B(285)
  JVS(772) = -B(281)+0.23*B(330)+0.025*B(378)
  JVS(773) = -B(345)-B(347)-B(349)-B(351)
  JVS(774) = -B(352)
  JVS(775) = -B(348)
  JVS(776) = -B(350)
  JVS(777) = -B(346)
  JVS(778) = -B(337)-B(339)-B(341)-B(343)
  JVS(779) = -B(344)
  JVS(780) = -B(340)
  JVS(781) = -B(342)
  JVS(782) = -B(338)
  JVS(783) = 0.097*B(367)
  JVS(784) = 0.118*B(365)
  JVS(785) = 0.5*B(315)
  JVS(786) = 0.5*B(318)
  JVS(787) = 0.607*B(357)
  JVS(788) = B(311)
  JVS(789) = 0.23*B(265)
  JVS(790) = 0.009*B(327)
  JVS(791) = 0
  JVS(792) = 0.001*B(347)
  JVS(793) = 0.001*B(339)
  JVS(794) = -B(253)-B(254)-B(255)-B(257)
  JVS(795) = 0.15*B(296)+0.023*B(298)
  JVS(796) = 0.009*B(328)
  JVS(797) = 0.023*B(299)+B(312)+0.001*B(340)+0.001*B(348)
  JVS(798) = -B(258)
  JVS(799) = 0
  JVS(800) = 0
  JVS(801) = 0
  JVS(802) = 0
  JVS(803) = 0
  JVS(804) = 0
  JVS(805) = -B(256)+0.23*B(266)+0.15*B(297)+0.607*B(358)+0.118*B(366)+0.097*B(368)
  JVS(806) = 0
  JVS(807) = 0.357*B(329)+0.936*B(333)
  JVS(808) = -B(296)-B(298)-B(300)-B(302)
  JVS(809) = 0.025*B(377)
  JVS(810) = 0
  JVS(811) = -B(299)
  JVS(812) = -B(301)+0.936*B(334)
  JVS(813) = -B(297)+0.357*B(330)+0.025*B(378)
  JVS(814) = -B(377)-B(379)-B(381)-B(383)
  JVS(815) = -B(384)
  JVS(816) = -B(380)
  JVS(817) = -B(382)
  JVS(818) = -B(378)
  JVS(819) = 0.32*B(329)+0.16*B(331)
  JVS(820) = 0.019*B(379)+0.048*B(381)
  JVS(821) = -B(289)-B(291)-B(293)-B(295)
  JVS(822) = -B(294)
  JVS(823) = -B(292)+0.16*B(332)+0.019*B(380)
  JVS(824) = 0.048*B(382)
  JVS(825) = -B(290)+0.32*B(330)
  JVS(826) = B(353)
  JVS(827) = 0.96*B(245)
  JVS(828) = 0.099*B(363)
  JVS(829) = 0.445*B(359)
  JVS(830) = 0.455*B(361)
  JVS(831) = 0.195*B(321)+0.25*B(327)
  JVS(832) = 0.984*B(387)+0.5*B(389)
  JVS(833) = 0.294*B(369)+0.154*B(371)+0.009*B(373)
  JVS(834) = 0.129*B(296)+0.047*B(298)+0.467*B(302)
  JVS(835) = 0.732*B(377)+0.456*B(379)+0.507*B(381)
  JVS(836) = -B(227)-B(229)-B(230)
  JVS(837) = 0.439*B(306)+0.431*B(308)
  JVS(838) = 0.25*B(328)
  JVS(839) = 0.034*B(232)+B(234)
  JVS(840) = 0.482*B(240)+B(242)
  JVS(841) = 0.084*B(303)+0.246*B(305)
  JVS(842) = 0.047*B(299)+0.154*B(372)+0.456*B(380)+0.5*B(390)
  JVS(843) = B(144)-B(231)+0.009*B(374)+0.507*B(382)
  JVS(844) = B(198)
  JVS(845) = B(140)
  JVS(846) = B(141)+B(145)+B(154)+2*B(156)+B(176)+B(199)
  JVS(847) = B(155)
  JVS(848) = B(177)
  JVS(849) = -B(228)+0.034*B(233)+0.482*B(241)+0.96*B(246)+0.129*B(297)+0.084*B(304)+0.439*B(307)+0.195*B(322)+B(354)&
               &+0.445*B(360)+0.455*B(362)+0.099*B(364)+0.294*B(370)+0.732*B(378)+0.984*B(388)
  JVS(850) = 0.081*B(245)
  JVS(851) = 0.026*B(363)
  JVS(852) = 0.026*B(359)
  JVS(853) = 0.35*B(247)+B(249)
  JVS(854) = B(222)
  JVS(855) = B(243)
  JVS(856) = 0.024*B(361)
  JVS(857) = 0.096*B(357)
  JVS(858) = 1.61*B(321)+B(323)+0.191*B(327)
  JVS(859) = B(237)
  JVS(860) = 0.984*B(387)+0.5*B(389)
  JVS(861) = 0.732*B(369)+0.5*B(371)
  JVS(862) = 0.624*B(329)+0.592*B(331)+0.24*B(335)
  JVS(863) = 0.084*B(280)+0.2*B(282)+0.67*B(288)
  JVS(864) = 0.276*B(345)+0.235*B(347)
  JVS(865) = 0.276*B(337)+0.235*B(339)
  JVS(866) = B(254)
  JVS(867) = 0.055*B(296)+0.125*B(298)+0.227*B(300)+0.3*B(302)
  JVS(868) = 0.244*B(377)+0.269*B(379)+0.079*B(381)
  JVS(869) = 0.3*B(289)+0.1*B(291)
  JVS(870) = -B(216)-B(217)-B(218)-B(220)-B(225)
  JVS(871) = 0.01*B(306)+0.134*B(308)
  JVS(872) = 0.191*B(328)+0.24*B(336)
  JVS(873) = 0.115*B(240)
  JVS(874) = 0.213*B(303)+0.506*B(305)
  JVS(875) = 0.2*B(283)+0.1*B(292)+0.125*B(299)+B(324)+0.592*B(332)+0.235*B(340)+0.235*B(348)+0.5*B(372)+0.269*B(380)&
               &+0.5*B(390)
  JVS(876) = B(82)+B(186)-B(226)+0.227*B(301)+0.079*B(382)
  JVS(877) = 0.75*B(110)
  JVS(878) = B(182)+B(187)+B(188)+B(196)+B(198)+B(200)+2*B(202)
  JVS(879) = B(78)+B(183)
  JVS(880) = -B(221)
  JVS(881) = B(146)+B(199)
  JVS(882) = B(128)+B(197)
  JVS(883) = B(166)+B(201)
  JVS(884) = -B(219)+B(238)+0.115*B(241)+B(244)+0.081*B(246)+0.35*B(248)+0.084*B(281)+0.3*B(290)+0.055*B(297)+0.213&
               &*B(304)+0.01*B(307)+1.61*B(322)+0.624*B(330)+0.276*B(338)+0.276*B(346)+0.096*B(358)+0.026*B(360)+0.024&
               &*B(362)+0.026*B(364)+0.732*B(370)+0.244*B(378)+0.984*B(388)
  JVS(885) = 0.75*B(92)
  JVS(886) = 0
  JVS(887) = B(79)+B(83)+B(84)+2*B(85)+0.75*B(93)+0.75*B(111)+B(129)+B(147)+B(167)+B(189)
  JVS(888) = B(203)
  JVS(889) = 0.511*B(373)
  JVS(890) = 0.276*B(349)
  JVS(891) = 0.276*B(341)
  JVS(892) = 0.572*B(300)
  JVS(893) = 0.321*B(381)
  JVS(894) = -0.69*B(306)-B(308)
  JVS(895) = 0
  JVS(896) = 0
  JVS(897) = 0.572*B(301)+0.276*B(342)+0.276*B(350)+0.511*B(374)+0.321*B(382)
  JVS(898) = B(106)
  JVS(899) = B(107)
  JVS(900) = -0.69*B(307)
  JVS(901) = B(204)
  JVS(902) = B(34)
  JVS(903) = -B(327)
  JVS(904) = -B(393)
  JVS(905) = -B(375)
  JVS(906) = -B(335)
  JVS(907) = -B(286)
  JVS(908) = -B(351)
  JVS(909) = -B(343)
  JVS(910) = -B(383)
  JVS(911) = -B(293)
  JVS(912) = -B(2)-B(4)-B(6)-B(9)-B(11)-B(287)-B(294)-B(328)-B(336)-B(344)-B(352)-B(376)-B(384)-B(394)
  JVS(913) = -B(5)+B(30)
  JVS(914) = B(29)
  JVS(915) = -B(7)
  JVS(916) = 0
  JVS(917) = B(1)-B(10)-B(12)
  JVS(918) = 0.261*B(355)
  JVS(919) = 0.204*B(363)
  JVS(920) = 0.122*B(359)
  JVS(921) = B(313)
  JVS(922) = B(316)
  JVS(923) = 0.244*B(361)
  JVS(924) = B(309)
  JVS(925) = B(250)+B(252)
  JVS(926) = B(325)
  JVS(927) = 0.45*B(393)
  JVS(928) = 0.497*B(369)+0.363*B(371)+0.037*B(373)+0.45*B(375)
  JVS(929) = B(286)
  JVS(930) = 0.474*B(345)+0.205*B(347)+0.474*B(349)+0.147*B(351)
  JVS(931) = 0.474*B(337)+0.205*B(339)+0.474*B(341)+0.147*B(343)
  JVS(932) = 0.013*B(296)+0.218*B(300)
  JVS(933) = 0.511*B(377)+0.305*B(379)+0.151*B(381)+0.069*B(383)
  JVS(934) = 0.675*B(289)+0.45*B(293)
  JVS(935) = 0.213*B(306)+0.147*B(308)
  JVS(936) = B(287)+0.45*B(294)+0.147*B(344)+0.147*B(352)+0.45*B(376)+0.069*B(384)+0.45*B(394)
  JVS(937) = -B(232)-B(234)-B(235)
  JVS(938) = 0.37*B(240)
  JVS(939) = 0.558*B(303)+0.71*B(305)
  JVS(940) = 0.205*B(340)+0.205*B(348)+0.363*B(372)+0.305*B(380)
  JVS(941) = -B(236)+0.218*B(301)+B(326)+0.474*B(342)+0.474*B(350)+0.037*B(374)+0.151*B(382)
  JVS(942) = 0
  JVS(943) = 0
  JVS(944) = 0
  JVS(945) = -B(233)+0.37*B(241)+B(251)+0.675*B(290)+0.013*B(297)+0.558*B(304)+0.213*B(307)+B(310)+B(314)+B(317)+0.474&
               &*B(338)+0.474*B(346)+0.261*B(356)+0.122*B(360)+0.244*B(362)+0.204*B(364)+0.497*B(370)+0.511*B(378)
  JVS(946) = 0
  JVS(947) = 0
  JVS(948) = 0.089*B(363)
  JVS(949) = 0.332*B(359)
  JVS(950) = 0.11*B(361)
  JVS(951) = 0.55*B(393)
  JVS(952) = 0.437*B(375)
  JVS(953) = 0.416*B(280)
  JVS(954) = 0.15*B(296)+0.21*B(298)+0.233*B(302)
  JVS(955) = 0.072*B(377)+0.026*B(379)+0.001*B(381)+0.659*B(383)
  JVS(956) = 0.55*B(293)
  JVS(957) = 0.177*B(306)+0.243*B(308)
  JVS(958) = 0.55*B(294)+0.437*B(376)+0.659*B(384)+0.55*B(394)
  JVS(959) = -B(240)-B(242)
  JVS(960) = 0.115*B(303)
  JVS(961) = 0.21*B(299)+0.026*B(380)
  JVS(962) = B(112)+0.001*B(382)
  JVS(963) = 0.5*B(110)+B(113)+0.5*B(114)+B(118)
  JVS(964) = 0
  JVS(965) = -B(241)+0.416*B(281)+0.15*B(297)+0.115*B(304)+0.177*B(307)+0.332*B(360)+0.11*B(362)+0.089*B(364)+0.072&
               &*B(378)
  JVS(966) = 0.5*B(115)
  JVS(967) = 0
  JVS(968) = 0.5*B(111)
  JVS(969) = 0.417*B(363)
  JVS(970) = 0.055*B(365)
  JVS(971) = 0.125*B(361)
  JVS(972) = 0.119*B(369)+0.215*B(371)+0.113*B(375)
  JVS(973) = 0.1*B(331)+0.75*B(335)
  JVS(974) = 0.276*B(345)+0.276*B(347)+0.853*B(351)
  JVS(975) = 0.276*B(337)+0.276*B(339)+0.853*B(343)
  JVS(976) = 0.332*B(296)
  JVS(977) = 0.043*B(379)+0.259*B(383)
  JVS(978) = 0.7*B(295)
  JVS(979) = 0.048*B(306)+0.435*B(308)
  JVS(980) = 0.75*B(336)+0.853*B(344)+0.853*B(352)+0.113*B(376)+0.259*B(384)
  JVS(981) = -0.671*B(303)-B(305)
  JVS(982) = 0.1*B(332)+0.276*B(340)+0.276*B(348)+0.215*B(372)+0.043*B(380)
  JVS(983) = 0
  JVS(984) = 0.5*B(110)+0.5*B(114)+B(118)+B(134)+B(152)+B(172)
  JVS(985) = 0
  JVS(986) = B(153)
  JVS(987) = B(135)
  JVS(988) = B(173)
  JVS(989) = 0.332*B(297)-0.671*B(304)+0.048*B(307)+0.276*B(338)+0.276*B(346)+0.125*B(362)+0.417*B(364)+0.055*B(366)&
               &+0.119*B(370)
  JVS(990) = 0.5*B(115)
  JVS(991) = 0
  JVS(992) = 0.5*B(111)
  JVS(993) = -B(385)
  JVS(994) = -B(311)
  JVS(995) = -B(323)
  JVS(996) = -B(389)
  JVS(997) = -B(371)
  JVS(998) = -B(331)
  JVS(999) = -B(282)
  JVS(1000) = -B(347)
  JVS(1001) = -B(339)
  JVS(1002) = -B(298)
  JVS(1003) = -B(379)
  JVS(1004) = -B(291)
  JVS(1005) = B(2)-B(4)
  JVS(1006) = -B(5)-B(13)-B(15)-B(30)-B(31)-B(51)-B(61)-B(283)-B(292)-B(299)-B(312)-B(324)-B(332)-B(340)-B(348)-B(372)&
                &-B(380)-B(386)-B(390)
  JVS(1007) = 0
  JVS(1008) = 0.25*B(184)
  JVS(1009) = -B(14)
  JVS(1010) = -B(62)+0.25*B(124)+0.25*B(142)+0.25*B(162)+0.25*B(185)
  JVS(1011) = 0.25*B(143)
  JVS(1012) = 0.25*B(125)
  JVS(1013) = 0.25*B(163)
  JVS(1014) = -B(52)
  JVS(1015) = -B(16)
  JVS(1016) = B(23)
  JVS(1017) = 0.39*B(58)
  JVS(1018) = -B(271)
  JVS(1019) = -B(273)
  JVS(1020) = -B(278)
  JVS(1021) = -B(267)
  JVS(1022) = -B(262)
  JVS(1023) = B(46)
  JVS(1024) = -B(325)
  JVS(1025) = -B(391)
  JVS(1026) = 0
  JVS(1027) = -B(373)
  JVS(1028) = -B(333)
  JVS(1029) = -B(99)
  JVS(1030) = -B(284)
  JVS(1031) = -B(349)
  JVS(1032) = -B(341)
  JVS(1033) = -B(257)
  JVS(1034) = -B(300)
  JVS(1035) = -B(381)
  JVS(1036) = 0
  JVS(1037) = -B(230)
  JVS(1038) = -B(225)
  JVS(1039) = 0
  JVS(1040) = B(11)
  JVS(1041) = -B(235)
  JVS(1042) = 0
  JVS(1043) = 0
  JVS(1044) = B(15)
  JVS(1045) = -B(17)-B(21)-B(26)-B(28)-B(29)-B(44)-B(66)-2*B(68)-B(82)-B(90)-B(100)-B(112)-B(126)-B(144)-B(164)-B(186)&
                &-B(226)-B(231)-B(236)-B(258)-B(263)-B(268)-B(272)-B(274)-B(279)-B(285)-B(301)-B(326)-B(334)-B(342)-B(350)&
                &-B(374)-B(382)-B(392)
  JVS(1046) = -B(113)
  JVS(1047) = -B(187)
  JVS(1048) = -B(18)
  JVS(1049) = -B(67)
  JVS(1050) = -B(145)
  JVS(1051) = -B(127)
  JVS(1052) = -B(165)
  JVS(1053) = -B(45)+B(47)
  JVS(1054) = -B(91)
  JVS(1055) = B(12)+B(16)-B(22)-B(27)
  JVS(1056) = -B(83)
  JVS(1057) = 0.035*B(355)
  JVS(1058) = 0.347*B(363)
  JVS(1059) = 0.07*B(359)
  JVS(1060) = 0.009*B(367)
  JVS(1061) = 0.011*B(365)
  JVS(1062) = 0.143*B(361)
  JVS(1063) = 0.016*B(387)+0.051*B(391)
  JVS(1064) = 0.09*B(369)+0.001*B(371)+0.176*B(373)
  JVS(1065) = 0.093*B(329)+0.008*B(331)+0.064*B(333)+0.01*B(335)
  JVS(1066) = 0.25*B(345)+0.18*B(347)+0.25*B(349)
  JVS(1067) = 0.25*B(337)+0.18*B(339)+0.25*B(341)
  JVS(1068) = 0.041*B(296)+0.051*B(300)
  JVS(1069) = 0.082*B(377)+0.002*B(379)+0.136*B(381)+0.001*B(383)
  JVS(1070) = 0.025*B(289)
  JVS(1071) = 0.173*B(306)+0.095*B(308)
  JVS(1072) = 0.01*B(336)+0.001*B(384)
  JVS(1073) = 0.001*B(232)
  JVS(1074) = 0.042*B(240)
  JVS(1075) = 0.07*B(303)+0.04*B(305)
  JVS(1076) = 0.008*B(332)+0.18*B(340)+0.18*B(348)+0.001*B(372)+0.002*B(380)
  JVS(1077) = -B(112)+0.051*B(301)+0.064*B(334)+0.25*B(342)+0.25*B(350)+0.176*B(374)+0.136*B(382)+0.051*B(392)
  JVS(1078) = -B(106)-B(108)-B(110)-B(113)-B(114)-2*B(118)-B(134)-B(152)-B(172)-B(194)
  JVS(1079) = -B(195)
  JVS(1080) = -B(107)
  JVS(1081) = -B(109)
  JVS(1082) = -B(153)
  JVS(1083) = -B(135)
  JVS(1084) = -B(173)
  JVS(1085) = 0.001*B(233)+0.042*B(241)+0.025*B(290)+0.041*B(297)+0.07*B(304)+0.173*B(307)+0.093*B(330)+0.25*B(338)+0.25&
                &*B(346)+0.035*B(356)+0.07*B(360)+0.143*B(362)+0.347*B(364)+0.011*B(366)+0.009*B(368)+0.09*B(370)+0.082&
                &*B(378)+0.016*B(388)
  JVS(1086) = -B(115)
  JVS(1087) = 0
  JVS(1088) = -B(111)
  JVS(1089) = B(181)
  JVS(1090) = 0.192*B(331)+0.24*B(335)
  JVS(1091) = 0.5*B(280)+0.5*B(284)+0.33*B(288)
  JVS(1092) = 0.289*B(296)+0.15*B(300)
  JVS(1093) = 0
  JVS(1094) = 0.3*B(295)
  JVS(1095) = 0.24*B(336)
  JVS(1096) = 0.192*B(332)
  JVS(1097) = -B(186)+0.5*B(285)+0.15*B(301)
  JVS(1098) = -B(194)
  JVS(1099) = -B(179)-B(182)-B(184)-B(187)-B(188)-B(190)-B(195)-B(196)-B(198)-B(200)-2*B(202)
  JVS(1100) = -B(183)
  JVS(1101) = -B(185)
  JVS(1102) = -B(199)
  JVS(1103) = -B(197)
  JVS(1104) = -B(201)
  JVS(1105) = 0.5*B(281)+0.289*B(297)
  JVS(1106) = -B(191)
  JVS(1107) = -B(180)
  JVS(1108) = -B(189)
  JVS(1109) = B(38)
  JVS(1110) = -B(223)
  JVS(1111) = -B(95)
  JVS(1112) = 0
  JVS(1113) = 0
  JVS(1114) = 0
  JVS(1115) = 0
  JVS(1116) = 0
  JVS(1117) = 0
  JVS(1118) = -B(6)+B(9)
  JVS(1119) = 0
  JVS(1120) = 0
  JVS(1121) = -B(13)
  JVS(1122) = -B(17)+B(26)+B(28)
  JVS(1123) = -B(106)
  JVS(1124) = -B(182)
  JVS(1125) = -B(7)-B(14)-B(18)-2*B(19)-B(36)-B(53)-B(78)-B(86)-B(96)-B(107)-B(122)-B(140)-B(160)-B(183)-B(224)
  JVS(1126) = -B(54)
  JVS(1127) = -B(141)
  JVS(1128) = -B(123)
  JVS(1129) = -B(161)
  JVS(1130) = -B(37)
  JVS(1131) = -B(87)
  JVS(1132) = B(1)+B(10)+B(27)
  JVS(1133) = -B(79)
  JVS(1134) = B(74)
  JVS(1135) = B(70)
  JVS(1136) = 0.95*B(245)
  JVS(1137) = B(39)
  JVS(1138) = 0.187*B(367)
  JVS(1139) = B(249)
  JVS(1140) = B(222)+B(223)
  JVS(1141) = -B(213)
  JVS(1142) = B(57)+0.61*B(58)
  JVS(1143) = B(243)
  JVS(1144) = 0.224*B(365)
  JVS(1145) = 0.5*B(315)
  JVS(1146) = 0.5*B(318)
  JVS(1147) = 0.297*B(357)+1.5*B(385)
  JVS(1148) = 1.5*B(311)
  JVS(1149) = 0
  JVS(1150) = B(252)
  JVS(1151) = B(259)
  JVS(1152) = B(49)
  JVS(1153) = 0.12*B(323)+0.5*B(327)
  JVS(1154) = 0.06*B(389)
  JVS(1155) = -B(208)
  JVS(1156) = 0.056*B(371)
  JVS(1157) = 0
  JVS(1158) = 0.008*B(282)+0.34*B(288)
  JVS(1159) = 0.033*B(347)
  JVS(1160) = 0.033*B(339)
  JVS(1161) = 2*B(253)+0.63*B(255)+0.63*B(257)
  JVS(1162) = 0.4*B(298)+1.233*B(302)
  JVS(1163) = 0.003*B(379)+0.013*B(383)
  JVS(1164) = 0.064*B(291)
  JVS(1165) = B(229)
  JVS(1166) = 2*B(216)+B(218)-B(220)+B(225)
  JVS(1167) = 0.113*B(306)+0.341*B(308)
  JVS(1168) = 0.5*B(328)+0.013*B(384)
  JVS(1169) = B(234)
  JVS(1170) = 0
  JVS(1171) = 0.379*B(303)
  JVS(1172) = B(51)-B(61)+0.008*B(283)+0.064*B(292)+0.4*B(299)+1.5*B(312)+0.12*B(324)+0.033*B(340)+0.033*B(348)+0.056&
                &*B(372)+0.003*B(380)+1.5*B(386)+0.06*B(390)
  JVS(1173) = B(44)-B(66)+B(82)+B(90)+B(112)+B(226)+0.63*B(258)
  JVS(1174) = -B(108)+B(110)+B(113)+B(114)+B(118)
  JVS(1175) = -B(184)
  JVS(1176) = -B(53)+B(78)+B(86)+B(224)
  JVS(1177) = -B(54)-B(55)-B(62)-2*B(63)-2*B(64)-B(67)-B(72)-B(80)-B(88)-B(109)-B(124)-B(142)-B(162)-B(185)-B(209)&
                &-B(214)-B(221)-B(396)
  JVS(1178) = -B(143)
  JVS(1179) = -B(125)
  JVS(1180) = -B(163)
  JVS(1181) = B(45)+B(50)+B(52)+B(71)-B(73)+B(75)+B(76)+B(219)+B(244)+0.95*B(246)+0.63*B(256)+0.379*B(304)+0.113*B(307)&
                &+0.297*B(358)+0.224*B(366)+0.187*B(368)
  JVS(1182) = B(87)-B(89)+B(91)+B(92)+B(94)+B(115)
  JVS(1183) = -B(56)
  JVS(1184) = B(79)-B(81)+B(83)+2*B(85)+B(93)+B(111)
  JVS(1185) = B(139)
  JVS(1186) = 0.1*B(282)
  JVS(1187) = 0.201*B(347)
  JVS(1188) = 0.201*B(339)
  JVS(1189) = 0.37*B(255)+0.37*B(257)
  JVS(1190) = 0.048*B(298)+0.3*B(302)
  JVS(1191) = 0.006*B(379)
  JVS(1192) = 0.05*B(291)
  JVS(1193) = 0
  JVS(1194) = 0.965*B(232)+B(235)
  JVS(1195) = 0.096*B(240)
  JVS(1196) = 0.049*B(303)+0.333*B(305)
  JVS(1197) = 0.1*B(283)+0.05*B(292)+0.048*B(299)+0.201*B(340)+0.201*B(348)+0.006*B(380)
  JVS(1198) = -B(144)+B(236)+0.37*B(258)
  JVS(1199) = -B(152)
  JVS(1200) = -B(198)
  JVS(1201) = -B(140)
  JVS(1202) = -B(142)
  JVS(1203) = -B(137)-B(141)-B(143)-B(145)-B(146)-B(148)-B(153)-B(154)-2*B(156)-B(176)-B(199)
  JVS(1204) = -B(155)
  JVS(1205) = -B(177)
  JVS(1206) = 0.965*B(233)+0.096*B(241)+0.37*B(256)+0.049*B(304)
  JVS(1207) = -B(149)
  JVS(1208) = -B(138)
  JVS(1209) = -B(147)
  JVS(1210) = B(121)
  JVS(1211) = 2*B(264)
  JVS(1212) = 0
  JVS(1213) = B(313)+0.5*B(315)
  JVS(1214) = B(316)+0.5*B(318)
  JVS(1215) = 0.011*B(361)
  JVS(1216) = B(259)+B(260)+B(262)
  JVS(1217) = B(237)+B(239)
  JVS(1218) = 0
  JVS(1219) = 0.67*B(288)
  JVS(1220) = 0.123*B(347)
  JVS(1221) = 0.123*B(339)
  JVS(1222) = 0.467*B(302)
  JVS(1223) = 0.137*B(379)
  JVS(1224) = 0.675*B(289)
  JVS(1225) = B(227)+B(230)
  JVS(1226) = 0
  JVS(1227) = 0
  JVS(1228) = 0
  JVS(1229) = 0.492*B(240)+B(242)
  JVS(1230) = 0.029*B(303)+0.667*B(305)
  JVS(1231) = 0.123*B(340)+0.123*B(348)+0.137*B(380)
  JVS(1232) = -B(126)+B(186)+B(231)+B(263)
  JVS(1233) = -B(134)
  JVS(1234) = B(182)+B(187)+B(198)+B(200)+2*B(202)
  JVS(1235) = -B(122)+B(183)
  JVS(1236) = -B(124)
  JVS(1237) = -B(154)+B(199)
  JVS(1238) = -B(119)-B(123)-B(125)-B(127)-B(128)-B(130)-B(135)-2*B(136)-B(155)-B(174)
  JVS(1239) = -B(175)+B(201)
  JVS(1240) = B(228)+B(238)+0.492*B(241)+B(261)+0.675*B(290)+0.029*B(304)+B(314)+B(317)+0.011*B(362)
  JVS(1241) = -B(131)
  JVS(1242) = -B(120)
  JVS(1243) = -B(129)
  JVS(1244) = B(159)
  JVS(1245) = B(275)+B(278)
  JVS(1246) = 0
  JVS(1247) = 0
  JVS(1248) = 0
  JVS(1249) = -B(164)+B(279)
  JVS(1250) = -B(172)
  JVS(1251) = -B(200)
  JVS(1252) = -B(160)
  JVS(1253) = -B(162)
  JVS(1254) = -B(176)
  JVS(1255) = -B(174)
  JVS(1256) = -B(157)-B(161)-B(163)-B(165)-B(166)-B(168)-B(173)-B(175)-B(177)-2*B(178)-B(201)
  JVS(1257) = B(276)
  JVS(1258) = -B(169)
  JVS(1259) = -B(158)
  JVS(1260) = -B(167)
  JVS(1261) = -B(510)
  JVS(1262) = -B(518)
  JVS(1263) = -B(526)
  JVS(1264) = -B(534)
  JVS(1265) = -B(542)
  JVS(1266) = -B(550)
  JVS(1267) = -B(558)
  JVS(1268) = -B(566)
  JVS(1269) = -B(574)
  JVS(1270) = -B(582)
  JVS(1271) = -B(514)
  JVS(1272) = -B(530)
  JVS(1273) = -B(538)
  JVS(1274) = -B(546)
  JVS(1275) = -B(554)
  JVS(1276) = -B(562)
  JVS(1277) = -B(570)
  JVS(1278) = -B(578)
  JVS(1279) = -B(586)
  JVS(1280) = -B(594)
  JVS(1281) = -B(602)
  JVS(1282) = -B(610)
  JVS(1283) = -B(618)
  JVS(1284) = -B(626)
  JVS(1285) = -B(634)
  JVS(1286) = -B(642)
  JVS(1287) = 2*B(32)
  JVS(1288) = -B(522)
  JVS(1289) = -B(590)
  JVS(1290) = -B(598)
  JVS(1291) = -B(606)
  JVS(1292) = -B(614)
  JVS(1293) = -B(622)
  JVS(1294) = -B(630)
  JVS(1295) = -B(638)
  JVS(1296) = -B(319)
  JVS(1297) = -B(74)
  JVS(1298) = -B(353)
  JVS(1299) = 2*B(69)-B(70)
  JVS(1300) = -B(355)
  JVS(1301) = -B(245)
  JVS(1302) = B(38)-B(40)
  JVS(1303) = -B(363)
  JVS(1304) = -B(359)
  JVS(1305) = -B(367)
  JVS(1306) = -0.65*B(247)+B(249)
  JVS(1307) = 0.39*B(58)-B(59)
  JVS(1308) = -B(243)
  JVS(1309) = -B(365)
  JVS(1310) = -B(313)
  JVS(1311) = -B(316)
  JVS(1312) = -B(269)
  JVS(1313) = -B(361)
  JVS(1314) = -0.397*B(357)+0.5*B(385)
  JVS(1315) = -B(309)+0.5*B(311)
  JVS(1316) = -0.34*B(250)+B(252)
  JVS(1317) = -B(275)
  JVS(1318) = -B(265)
  JVS(1319) = -B(260)
  JVS(1320) = -B(49)
  JVS(1321) = -B(46)+B(48)
  JVS(1322) = -B(321)+0.12*B(323)
  JVS(1323) = -B(237)
  JVS(1324) = -B(387)+0.32*B(389)
  JVS(1325) = 0
  JVS(1326) = -B(369)+0.155*B(371)
  JVS(1327) = -B(329)+0.266*B(331)
  JVS(1328) = -B(280)+0.208*B(282)+0.33*B(288)
  JVS(1329) = -B(345)+0.567*B(347)
  JVS(1330) = -B(337)+0.567*B(339)
  JVS(1331) = -B(255)
  JVS(1332) = -B(296)+0.285*B(298)
  JVS(1333) = -B(377)+0.378*B(379)
  JVS(1334) = -B(289)+0.164*B(291)
  JVS(1335) = -B(227)
  JVS(1336) = -B(218)
  JVS(1337) = -B(306)
  JVS(1338) = 0
  JVS(1339) = -B(232)
  JVS(1340) = -B(240)
  JVS(1341) = -B(303)
  JVS(1342) = -B(51)+B(61)+0.208*B(283)+0.164*B(292)+0.285*B(299)+0.5*B(312)+0.12*B(324)+0.266*B(332)+0.567*B(340)+0.567&
                &*B(348)+0.155*B(372)+0.378*B(380)+0.5*B(386)+0.32*B(390)
  JVS(1343) = -B(44)+0.8*B(66)
  JVS(1344) = 0
  JVS(1345) = 0
  JVS(1346) = -B(36)+B(53)
  JVS(1347) = B(54)+B(62)+0.8*B(67)-B(72)
  JVS(1348) = 0
  JVS(1349) = 0
  JVS(1350) = 0
  JVS(1351) = -B(37)-B(41)-B(42)-B(45)-B(47)-B(50)-B(52)-B(60)-B(71)-B(73)-B(75)-B(76)-B(219)-B(228)-B(233)-B(238)&
                &-B(241)-B(244)-B(246)-0.65*B(248)-0.34*B(251)-B(256)-B(261)-B(266)-B(270)-B(276)-B(281)-B(290)-B(297)&
                &-B(304)-B(307)-B(310)-B(314)-B(317)-B(320)-B(322)-B(330)-B(338)-B(346)-B(354)-B(356)-0.397*B(358)-B(360)&
                &-B(362)-B(364)-B(366)-B(368)-B(370)-B(378)-B(388)-B(511)-B(515)-B(519)-B(523)-B(527)-B(531)-B(535)-B(539)&
                &-B(543)-B(547)-B(551)-B(555)-B(559)-B(563)-B(567)-B(571)-B(575)-B(579)-B(583)-B(587)-B(591)-B(595)-B(599)&
                &-B(603)-B(607)-B(611)-B(615)-B(619)-B(623)-B(627)-B(631)-B(635)-B(639)-B(643)
  JVS(1352) = 0
  JVS(1353) = -B(43)
  JVS(1354) = 0
  JVS(1355) = B(353)
  JVS(1356) = 0.965*B(355)
  JVS(1357) = 0.05*B(245)
  JVS(1358) = 0.653*B(363)
  JVS(1359) = 0.695*B(359)
  JVS(1360) = 0.804*B(367)
  JVS(1361) = 0.765*B(365)
  JVS(1362) = B(315)
  JVS(1363) = B(318)
  JVS(1364) = 0.76*B(269)
  JVS(1365) = 0.835*B(361)
  JVS(1366) = 0.1*B(357)
  JVS(1367) = B(309)
  JVS(1368) = 0.34*B(250)
  JVS(1369) = 0.76*B(265)
  JVS(1370) = B(321)+B(325)+0.2*B(327)
  JVS(1371) = 0.984*B(387)+0.949*B(391)
  JVS(1372) = 0
  JVS(1373) = 0.91*B(369)+0.022*B(371)+0.824*B(373)
  JVS(1374) = 0.907*B(329)+0.066*B(331)+0.749*B(333)
  JVS(1375) = 0.5*B(280)+0.1*B(282)+0.5*B(284)+0.33*B(288)
  JVS(1376) = 0.75*B(345)+0.031*B(347)+0.276*B(349)
  JVS(1377) = 0.75*B(337)+0.031*B(339)+0.276*B(341)
  JVS(1378) = 0.67*B(296)+0.048*B(298)+0.799*B(300)
  JVS(1379) = 0.918*B(377)+0.033*B(379)+0.442*B(381)+0.012*B(383)
  JVS(1380) = 0.3*B(289)+0.05*B(291)
  JVS(1381) = 0.376*B(306)+0.564*B(308)
  JVS(1382) = 0.2*B(328)+0.012*B(384)
  JVS(1383) = 0.034*B(232)+B(234)
  JVS(1384) = 0.37*B(240)+B(242)
  JVS(1385) = 0.473*B(303)+0.96*B(305)
  JVS(1386) = 0.1*B(283)+0.05*B(292)+0.048*B(299)+0.066*B(332)+0.031*B(340)+0.031*B(348)+0.022*B(372)+0.033*B(380)
  JVS(1387) = -B(90)+B(144)+0.5*B(285)+0.799*B(301)+B(326)+0.749*B(334)+0.276*B(342)+0.276*B(350)+0.824*B(374)+0.442&
                &*B(382)+0.949*B(392)
  JVS(1388) = -B(114)
  JVS(1389) = -B(190)+B(198)
  JVS(1390) = -B(86)+B(140)
  JVS(1391) = -B(88)
  JVS(1392) = B(141)+B(145)-B(148)+B(154)+2*B(156)+B(176)+B(199)
  JVS(1393) = -B(130)+B(155)
  JVS(1394) = -B(168)+B(177)
  JVS(1395) = 0.034*B(233)+0.37*B(241)+0.05*B(246)+0.34*B(251)+0.76*B(266)+0.76*B(270)+0.5*B(281)+0.3*B(290)+0.67*B(297)&
                &+0.473*B(304)+0.376*B(307)+B(310)+B(322)+0.907*B(330)+0.75*B(338)+0.75*B(346)+B(354)+0.965*B(356)+0.1&
                &*B(358)+0.695*B(360)+0.835*B(362)+0.653*B(364)+0.765*B(366)+0.804*B(368)+0.91*B(370)+0.918*B(378)+0.984&
                &*B(388)
  JVS(1396) = -B(87)-B(89)-B(91)-B(92)-2*B(94)-B(115)-B(131)-B(149)-B(169)-B(191)
  JVS(1397) = 0
  JVS(1398) = -B(93)
  JVS(1399) = B(139)
  JVS(1400) = B(181)
  JVS(1401) = B(121)
  JVS(1402) = B(159)
  JVS(1403) = B(23)
  JVS(1404) = B(39)+B(40)
  JVS(1405) = -B(203)
  JVS(1406) = B(223)
  JVS(1407) = -B(211)
  JVS(1408) = B(57)+0.61*B(58)+B(59)
  JVS(1409) = 0
  JVS(1410) = B(48)
  JVS(1411) = -B(206)
  JVS(1412) = 0.187*B(333)
  JVS(1413) = B(95)+B(99)
  JVS(1414) = 0
  JVS(1415) = 0.474*B(349)
  JVS(1416) = 0.474*B(341)
  JVS(1417) = 0
  JVS(1418) = 0
  JVS(1419) = 0.391*B(381)
  JVS(1420) = 0
  JVS(1421) = 0
  JVS(1422) = 0
  JVS(1423) = 0.338*B(306)+B(308)
  JVS(1424) = B(6)-B(9)-B(11)
  JVS(1425) = 0
  JVS(1426) = 0
  JVS(1427) = 0
  JVS(1428) = B(13)-B(15)
  JVS(1429) = 2*B(17)-B(21)+B(29)+B(44)+0.8*B(66)+2*B(68)+B(82)+B(90)+B(100)+B(112)+B(126)+B(144)+B(164)+B(186)+0.187&
                &*B(334)+0.474*B(342)+0.474*B(350)+0.391*B(382)
  JVS(1430) = B(113)
  JVS(1431) = -B(179)+B(182)+B(187)
  JVS(1432) = B(7)+B(14)+2*B(18)+2*B(19)+B(53)+B(78)+B(86)+B(96)+B(122)+B(140)+B(160)+B(183)+B(224)
  JVS(1433) = B(54)-B(55)+0.8*B(67)
  JVS(1434) = -B(137)+B(141)+B(145)
  JVS(1435) = -B(119)+B(123)+B(127)
  JVS(1436) = -B(157)+B(161)+B(165)
  JVS(1437) = B(41)-B(42)+B(45)+B(60)+0.338*B(307)
  JVS(1438) = B(87)+B(91)
  JVS(1439) = -B(1)-B(10)-B(12)-B(16)-B(22)-B(43)-B(56)-B(120)-B(138)-B(158)-B(180)-B(204)-B(207)-B(212)
  JVS(1440) = B(79)+B(83)
  JVS(1441) = B(319)
  JVS(1442) = B(205)
  JVS(1443) = 0.65*B(247)
  JVS(1444) = 0.011*B(361)
  JVS(1445) = 0.3*B(327)
  JVS(1446) = B(239)
  JVS(1447) = 0.26*B(389)
  JVS(1448) = 0.076*B(371)
  JVS(1449) = 0.25*B(335)
  JVS(1450) = 0
  JVS(1451) = 0
  JVS(1452) = 0.197*B(379)+0.03*B(381)
  JVS(1453) = 0.3*B(295)
  JVS(1454) = B(229)
  JVS(1455) = 0
  JVS(1456) = 0.3*B(328)+0.25*B(336)
  JVS(1457) = 0
  JVS(1458) = 0
  JVS(1459) = 0
  JVS(1460) = 0.076*B(372)+0.197*B(380)+0.26*B(390)
  JVS(1461) = -B(82)+B(126)+0.03*B(382)
  JVS(1462) = -B(110)
  JVS(1463) = -B(188)+B(196)
  JVS(1464) = -B(78)+B(122)
  JVS(1465) = -B(80)
  JVS(1466) = -B(146)+B(154)
  JVS(1467) = B(123)+B(127)-B(128)+2*B(136)+B(155)+B(174)+B(197)
  JVS(1468) = -B(166)+B(175)
  JVS(1469) = 0.65*B(248)+B(320)+0.011*B(362)
  JVS(1470) = -B(92)
  JVS(1471) = 0
  JVS(1472) = -B(79)-B(81)-B(83)-2*B(84)-2*B(85)-B(93)-B(111)-B(129)-B(147)-B(167)-B(189)
END SUBROUTINE saprc99_mosaic_4bin_vbs9_Jac_SP
SUBROUTINE saprc99_mosaic_4bin_vbs9_KppDecomp( JVS, IER )
      INTEGER :: IER
      REAL(kind=dp) :: JVS(1472), W(168), a
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
END SUBROUTINE saprc99_mosaic_4bin_vbs9_KppDecomp
SUBROUTINE saprc99_mosaic_4bin_vbs9_KppDecompCmplx( JVS, IER )
      INTEGER :: IER
      DOUBLE COMPLEX :: JVS(1472), W(168), a
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
END SUBROUTINE saprc99_mosaic_4bin_vbs9_KppDecompCmplx
SUBROUTINE saprc99_mosaic_4bin_vbs9_KppSolveIndirect( JVS, X )
      INTEGER i, j
      REAL(kind=dp) JVS(1472), X(168), sum
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
END SUBROUTINE saprc99_mosaic_4bin_vbs9_KppSolveIndirect
SUBROUTINE saprc99_mosaic_4bin_vbs9_KppSolveCmplx( JVS, X )
      INTEGER i, j
      DOUBLE COMPLEX JVS(1472), X(168), sum
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
END SUBROUTINE saprc99_mosaic_4bin_vbs9_KppSolveCmplx
SUBROUTINE saprc99_mosaic_4bin_vbs9_KppSolve ( JVS, X )
  REAL(kind=dp) :: JVS(LU_NONZERO)
  REAL(kind=dp) :: X(NVAR)
  X(33) = X(33)-JVS(223)*X(3)-JVS(224)*X(4)-JVS(225)*X(5)-JVS(226)*X(6)
  X(74) = X(74)-JVS(433)*X(66)
  X(75) = X(75)-JVS(437)*X(67)
  X(76) = X(76)-JVS(441)*X(68)
  X(77) = X(77)-JVS(445)*X(69)
  X(78) = X(78)-JVS(449)*X(70)
  X(79) = X(79)-JVS(453)*X(71)
  X(80) = X(80)-JVS(457)*X(72)
  X(81) = X(81)-JVS(461)*X(73)
  X(91) = X(91)-JVS(482)*X(82)
  X(92) = X(92)-JVS(486)*X(83)
  X(93) = X(93)-JVS(490)*X(84)
  X(94) = X(94)-JVS(494)*X(85)
  X(95) = X(95)-JVS(498)*X(86)
  X(96) = X(96)-JVS(502)*X(87)
  X(97) = X(97)-JVS(506)*X(88)
  X(98) = X(98)-JVS(510)*X(89)
  X(114) = X(114)-JVS(554)*X(113)
  X(122) = X(122)-JVS(584)*X(115)-JVS(585)*X(121)
  X(123) = X(123)-JVS(588)*X(115)-JVS(589)*X(121)
  X(124) = X(124)-JVS(592)*X(115)-JVS(593)*X(121)
  X(127) = X(127)-JVS(602)*X(115)-JVS(603)*X(121)
  X(128) = X(128)-JVS(607)*X(118)
  X(130) = X(130)-JVS(618)*X(115)-JVS(619)*X(121)
  X(131) = X(131)-JVS(625)*X(121)
  X(132) = X(132)-JVS(631)*X(115)-JVS(632)*X(121)-JVS(633)*X(122)-JVS(634)*X(123)-JVS(635)*X(124)
  X(133) = X(133)-JVS(643)*X(122)-JVS(644)*X(123)-JVS(645)*X(125)-JVS(646)*X(126)-JVS(647)*X(127)-JVS(648)*X(132)
  X(134) = X(134)-JVS(668)*X(110)-JVS(669)*X(124)-JVS(670)*X(128)-JVS(671)*X(130)-JVS(672)*X(131)-JVS(673)*X(132)
  X(136) = X(136)-JVS(694)*X(107)-JVS(695)*X(112)-JVS(696)*X(113)-JVS(697)*X(114)-JVS(698)*X(125)
  X(138) = X(138)-JVS(714)*X(124)-JVS(715)*X(131)
  X(141) = X(141)-JVS(736)*X(112)-JVS(737)*X(113)-JVS(738)*X(122)-JVS(739)*X(123)-JVS(740)*X(125)-JVS(741)*X(136)&
             &-JVS(742)*X(139)-JVS(743)*X(140)
  X(142) = X(142)-JVS(766)*X(140)
  X(145) = X(145)-JVS(783)*X(115)-JVS(784)*X(121)-JVS(785)*X(122)-JVS(786)*X(123)-JVS(787)*X(126)-JVS(788)*X(127)&
             &-JVS(789)*X(131)-JVS(790)*X(135)-JVS(791)*X(138)-JVS(792)*X(143)-JVS(793)*X(144)
  X(146) = X(146)-JVS(807)*X(140)
  X(148) = X(148)-JVS(819)*X(140)-JVS(820)*X(147)
  X(149) = X(149)-JVS(826)*X(101)-JVS(827)*X(108)-JVS(828)*X(112)-JVS(829)*X(113)-JVS(830)*X(125)-JVS(831)*X(135)&
             &-JVS(832)*X(137)-JVS(833)*X(139)-JVS(834)*X(146)-JVS(835)*X(147)
  X(150) = X(150)-JVS(850)*X(108)-JVS(851)*X(112)-JVS(852)*X(113)-JVS(853)*X(116)-JVS(854)*X(117)-JVS(855)*X(120)&
             &-JVS(856)*X(125)-JVS(857)*X(126)-JVS(858)*X(135)-JVS(859)*X(136)-JVS(860)*X(137)-JVS(861)*X(139)-JVS(862)&
             &*X(140)-JVS(863)*X(142)-JVS(864)*X(143)-JVS(865)*X(144)-JVS(866)*X(145)-JVS(867)*X(146)-JVS(868)*X(147)&
             &-JVS(869)*X(148)
  X(151) = X(151)-JVS(888)*X(114)-JVS(889)*X(139)-JVS(890)*X(143)-JVS(891)*X(144)-JVS(892)*X(146)-JVS(893)*X(147)
  X(152) = X(152)-JVS(902)*X(90)-JVS(903)*X(135)-JVS(904)*X(137)-JVS(905)*X(139)-JVS(906)*X(140)-JVS(907)*X(142)&
             &-JVS(908)*X(143)-JVS(909)*X(144)-JVS(910)*X(147)-JVS(911)*X(148)
  X(153) = X(153)-JVS(918)*X(107)-JVS(919)*X(112)-JVS(920)*X(113)-JVS(921)*X(122)-JVS(922)*X(123)-JVS(923)*X(125)&
             &-JVS(924)*X(127)-JVS(925)*X(129)-JVS(926)*X(135)-JVS(927)*X(137)-JVS(928)*X(139)-JVS(929)*X(142)-JVS(930)&
             &*X(143)-JVS(931)*X(144)-JVS(932)*X(146)-JVS(933)*X(147)-JVS(934)*X(148)-JVS(935)*X(151)-JVS(936)*X(152)
  X(154) = X(154)-JVS(948)*X(112)-JVS(949)*X(113)-JVS(950)*X(125)-JVS(951)*X(137)-JVS(952)*X(139)-JVS(953)*X(142)&
             &-JVS(954)*X(146)-JVS(955)*X(147)-JVS(956)*X(148)-JVS(957)*X(151)-JVS(958)*X(152)
  X(155) = X(155)-JVS(969)*X(112)-JVS(970)*X(121)-JVS(971)*X(125)-JVS(972)*X(139)-JVS(973)*X(140)-JVS(974)*X(143)&
             &-JVS(975)*X(144)-JVS(976)*X(146)-JVS(977)*X(147)-JVS(978)*X(148)-JVS(979)*X(151)-JVS(980)*X(152)
  X(156) = X(156)-JVS(993)*X(126)-JVS(994)*X(127)-JVS(995)*X(135)-JVS(996)*X(137)-JVS(997)*X(139)-JVS(998)*X(140)&
             &-JVS(999)*X(142)-JVS(1000)*X(143)-JVS(1001)*X(144)-JVS(1002)*X(146)-JVS(1003)*X(147)-JVS(1004)*X(148)&
             &-JVS(1005)*X(152)
  X(157) = X(157)-JVS(1016)*X(110)-JVS(1017)*X(119)-JVS(1018)*X(124)-JVS(1019)*X(128)-JVS(1020)*X(130)-JVS(1021)*X(131)&
             &-JVS(1022)*X(132)-JVS(1023)*X(134)-JVS(1024)*X(135)-JVS(1025)*X(137)-JVS(1026)*X(138)-JVS(1027)*X(139)&
             &-JVS(1028)*X(140)-JVS(1029)*X(141)-JVS(1030)*X(142)-JVS(1031)*X(143)-JVS(1032)*X(144)-JVS(1033)*X(145)&
             &-JVS(1034)*X(146)-JVS(1035)*X(147)-JVS(1036)*X(148)-JVS(1037)*X(149)-JVS(1038)*X(150)-JVS(1039)*X(151)&
             &-JVS(1040)*X(152)-JVS(1041)*X(153)-JVS(1042)*X(154)-JVS(1043)*X(155)-JVS(1044)*X(156)
  X(158) = X(158)-JVS(1057)*X(107)-JVS(1058)*X(112)-JVS(1059)*X(113)-JVS(1060)*X(115)-JVS(1061)*X(121)-JVS(1062)*X(125)&
             &-JVS(1063)*X(137)-JVS(1064)*X(139)-JVS(1065)*X(140)-JVS(1066)*X(143)-JVS(1067)*X(144)-JVS(1068)*X(146)&
             &-JVS(1069)*X(147)-JVS(1070)*X(148)-JVS(1071)*X(151)-JVS(1072)*X(152)-JVS(1073)*X(153)-JVS(1074)*X(154)&
             &-JVS(1075)*X(155)-JVS(1076)*X(156)-JVS(1077)*X(157)
  X(159) = X(159)-JVS(1089)*X(103)-JVS(1090)*X(140)-JVS(1091)*X(142)-JVS(1092)*X(146)-JVS(1093)*X(147)-JVS(1094)*X(148)&
             &-JVS(1095)*X(152)-JVS(1096)*X(156)-JVS(1097)*X(157)-JVS(1098)*X(158)
  X(160) = X(160)-JVS(1109)*X(111)-JVS(1110)*X(117)-JVS(1111)*X(141)-JVS(1112)*X(143)-JVS(1113)*X(144)-JVS(1114)*X(147)&
             &-JVS(1115)*X(148)-JVS(1116)*X(150)-JVS(1117)*X(151)-JVS(1118)*X(152)-JVS(1119)*X(154)-JVS(1120)*X(155)&
             &-JVS(1121)*X(156)-JVS(1122)*X(157)-JVS(1123)*X(158)-JVS(1124)*X(159)
  X(161) = X(161)-JVS(1134)*X(100)-JVS(1135)*X(106)-JVS(1136)*X(108)-JVS(1137)*X(111)-JVS(1138)*X(115)-JVS(1139)*X(116)&
             &-JVS(1140)*X(117)-JVS(1141)*X(118)-JVS(1142)*X(119)-JVS(1143)*X(120)-JVS(1144)*X(121)-JVS(1145)*X(122)&
             &-JVS(1146)*X(123)-JVS(1147)*X(126)-JVS(1148)*X(127)-JVS(1149)*X(128)-JVS(1150)*X(129)-JVS(1151)*X(132)&
             &-JVS(1152)*X(133)-JVS(1153)*X(135)-JVS(1154)*X(137)-JVS(1155)*X(138)-JVS(1156)*X(139)-JVS(1157)*X(140)&
             &-JVS(1158)*X(142)-JVS(1159)*X(143)-JVS(1160)*X(144)-JVS(1161)*X(145)-JVS(1162)*X(146)-JVS(1163)*X(147)&
             &-JVS(1164)*X(148)-JVS(1165)*X(149)-JVS(1166)*X(150)-JVS(1167)*X(151)-JVS(1168)*X(152)-JVS(1169)*X(153)&
             &-JVS(1170)*X(154)-JVS(1171)*X(155)-JVS(1172)*X(156)-JVS(1173)*X(157)-JVS(1174)*X(158)-JVS(1175)*X(159)&
             &-JVS(1176)*X(160)
  X(162) = X(162)-JVS(1185)*X(102)-JVS(1186)*X(142)-JVS(1187)*X(143)-JVS(1188)*X(144)-JVS(1189)*X(145)-JVS(1190)*X(146)&
             &-JVS(1191)*X(147)-JVS(1192)*X(148)-JVS(1193)*X(152)-JVS(1194)*X(153)-JVS(1195)*X(154)-JVS(1196)*X(155)&
             &-JVS(1197)*X(156)-JVS(1198)*X(157)-JVS(1199)*X(158)-JVS(1200)*X(159)-JVS(1201)*X(160)-JVS(1202)*X(161)
  X(163) = X(163)-JVS(1210)*X(104)-JVS(1211)*X(109)-JVS(1212)*X(115)-JVS(1213)*X(122)-JVS(1214)*X(123)-JVS(1215)*X(125)&
             &-JVS(1216)*X(132)-JVS(1217)*X(136)-JVS(1218)*X(139)-JVS(1219)*X(142)-JVS(1220)*X(143)-JVS(1221)*X(144)&
             &-JVS(1222)*X(146)-JVS(1223)*X(147)-JVS(1224)*X(148)-JVS(1225)*X(149)-JVS(1226)*X(151)-JVS(1227)*X(152)&
             &-JVS(1228)*X(153)-JVS(1229)*X(154)-JVS(1230)*X(155)-JVS(1231)*X(156)-JVS(1232)*X(157)-JVS(1233)*X(158)&
             &-JVS(1234)*X(159)-JVS(1235)*X(160)-JVS(1236)*X(161)-JVS(1237)*X(162)
  X(164) = X(164)-JVS(1244)*X(105)-JVS(1245)*X(130)-JVS(1246)*X(147)-JVS(1247)*X(152)-JVS(1248)*X(156)-JVS(1249)*X(157)&
             &-JVS(1250)*X(158)-JVS(1251)*X(159)-JVS(1252)*X(160)-JVS(1253)*X(161)-JVS(1254)*X(162)-JVS(1255)*X(163)
  X(165) = X(165)-JVS(1261)*X(3)-JVS(1262)*X(5)-JVS(1263)*X(66)-JVS(1264)*X(67)-JVS(1265)*X(68)-JVS(1266)*X(69)&
             &-JVS(1267)*X(70)-JVS(1268)*X(71)-JVS(1269)*X(72)-JVS(1270)*X(73)-JVS(1271)*X(74)-JVS(1272)*X(75)-JVS(1273)&
             &*X(76)-JVS(1274)*X(77)-JVS(1275)*X(78)-JVS(1276)*X(79)-JVS(1277)*X(80)-JVS(1278)*X(81)-JVS(1279)*X(82)&
             &-JVS(1280)*X(83)-JVS(1281)*X(84)-JVS(1282)*X(85)-JVS(1283)*X(86)-JVS(1284)*X(87)-JVS(1285)*X(88)-JVS(1286)&
             &*X(89)-JVS(1287)*X(90)-JVS(1288)*X(91)-JVS(1289)*X(92)-JVS(1290)*X(93)-JVS(1291)*X(94)-JVS(1292)*X(95)&
             &-JVS(1293)*X(96)-JVS(1294)*X(97)-JVS(1295)*X(98)-JVS(1296)*X(99)-JVS(1297)*X(100)-JVS(1298)*X(101)-JVS(1299)&
             &*X(106)-JVS(1300)*X(107)-JVS(1301)*X(108)-JVS(1302)*X(111)-JVS(1303)*X(112)-JVS(1304)*X(113)-JVS(1305)*X(115)&
             &-JVS(1306)*X(116)-JVS(1307)*X(119)-JVS(1308)*X(120)-JVS(1309)*X(121)-JVS(1310)*X(122)-JVS(1311)*X(123)&
             &-JVS(1312)*X(124)-JVS(1313)*X(125)-JVS(1314)*X(126)-JVS(1315)*X(127)-JVS(1316)*X(129)-JVS(1317)*X(130)&
             &-JVS(1318)*X(131)-JVS(1319)*X(132)-JVS(1320)*X(133)-JVS(1321)*X(134)-JVS(1322)*X(135)-JVS(1323)*X(136)&
             &-JVS(1324)*X(137)-JVS(1325)*X(138)-JVS(1326)*X(139)-JVS(1327)*X(140)-JVS(1328)*X(142)-JVS(1329)*X(143)&
             &-JVS(1330)*X(144)-JVS(1331)*X(145)-JVS(1332)*X(146)-JVS(1333)*X(147)-JVS(1334)*X(148)-JVS(1335)*X(149)&
             &-JVS(1336)*X(150)-JVS(1337)*X(151)-JVS(1338)*X(152)-JVS(1339)*X(153)-JVS(1340)*X(154)-JVS(1341)*X(155)&
             &-JVS(1342)*X(156)-JVS(1343)*X(157)-JVS(1344)*X(158)-JVS(1345)*X(159)-JVS(1346)*X(160)-JVS(1347)*X(161)&
             &-JVS(1348)*X(162)-JVS(1349)*X(163)-JVS(1350)*X(164)
  X(166) = X(166)-JVS(1355)*X(101)-JVS(1356)*X(107)-JVS(1357)*X(108)-JVS(1358)*X(112)-JVS(1359)*X(113)-JVS(1360)*X(115)&
             &-JVS(1361)*X(121)-JVS(1362)*X(122)-JVS(1363)*X(123)-JVS(1364)*X(124)-JVS(1365)*X(125)-JVS(1366)*X(126)&
             &-JVS(1367)*X(127)-JVS(1368)*X(129)-JVS(1369)*X(131)-JVS(1370)*X(135)-JVS(1371)*X(137)-JVS(1372)*X(138)&
             &-JVS(1373)*X(139)-JVS(1374)*X(140)-JVS(1375)*X(142)-JVS(1376)*X(143)-JVS(1377)*X(144)-JVS(1378)*X(146)&
             &-JVS(1379)*X(147)-JVS(1380)*X(148)-JVS(1381)*X(151)-JVS(1382)*X(152)-JVS(1383)*X(153)-JVS(1384)*X(154)&
             &-JVS(1385)*X(155)-JVS(1386)*X(156)-JVS(1387)*X(157)-JVS(1388)*X(158)-JVS(1389)*X(159)-JVS(1390)*X(160)&
             &-JVS(1391)*X(161)-JVS(1392)*X(162)-JVS(1393)*X(163)-JVS(1394)*X(164)-JVS(1395)*X(165)
  X(167) = X(167)-JVS(1399)*X(102)-JVS(1400)*X(103)-JVS(1401)*X(104)-JVS(1402)*X(105)-JVS(1403)*X(110)-JVS(1404)*X(111)&
             &-JVS(1405)*X(114)-JVS(1406)*X(117)-JVS(1407)*X(118)-JVS(1408)*X(119)-JVS(1409)*X(128)-JVS(1410)*X(134)&
             &-JVS(1411)*X(138)-JVS(1412)*X(140)-JVS(1413)*X(141)-JVS(1414)*X(142)-JVS(1415)*X(143)-JVS(1416)*X(144)&
             &-JVS(1417)*X(145)-JVS(1418)*X(146)-JVS(1419)*X(147)-JVS(1420)*X(148)-JVS(1421)*X(149)-JVS(1422)*X(150)&
             &-JVS(1423)*X(151)-JVS(1424)*X(152)-JVS(1425)*X(153)-JVS(1426)*X(154)-JVS(1427)*X(155)-JVS(1428)*X(156)&
             &-JVS(1429)*X(157)-JVS(1430)*X(158)-JVS(1431)*X(159)-JVS(1432)*X(160)-JVS(1433)*X(161)-JVS(1434)*X(162)&
             &-JVS(1435)*X(163)-JVS(1436)*X(164)-JVS(1437)*X(165)-JVS(1438)*X(166)
  X(168) = X(168)-JVS(1441)*X(99)-JVS(1442)*X(114)-JVS(1443)*X(116)-JVS(1444)*X(125)-JVS(1445)*X(135)-JVS(1446)*X(136)&
             &-JVS(1447)*X(137)-JVS(1448)*X(139)-JVS(1449)*X(140)-JVS(1450)*X(143)-JVS(1451)*X(144)-JVS(1452)*X(147)&
             &-JVS(1453)*X(148)-JVS(1454)*X(149)-JVS(1455)*X(151)-JVS(1456)*X(152)-JVS(1457)*X(153)-JVS(1458)*X(154)&
             &-JVS(1459)*X(155)-JVS(1460)*X(156)-JVS(1461)*X(157)-JVS(1462)*X(158)-JVS(1463)*X(159)-JVS(1464)*X(160)&
             &-JVS(1465)*X(161)-JVS(1466)*X(162)-JVS(1467)*X(163)-JVS(1468)*X(164)-JVS(1469)*X(165)-JVS(1470)*X(166)&
             &-JVS(1471)*X(167)
  X(168) = X(168)/JVS(1472)
  X(167) = (X(167)-JVS(1440)*X(168))/(JVS(1439))
  X(166) = (X(166)-JVS(1397)*X(167)-JVS(1398)*X(168))/(JVS(1396))
  X(165) = (X(165)-JVS(1352)*X(166)-JVS(1353)*X(167)-JVS(1354)*X(168))/(JVS(1351))
  X(164) = (X(164)-JVS(1257)*X(165)-JVS(1258)*X(166)-JVS(1259)*X(167)-JVS(1260)*X(168))/(JVS(1256))
  X(163) = (X(163)-JVS(1239)*X(164)-JVS(1240)*X(165)-JVS(1241)*X(166)-JVS(1242)*X(167)-JVS(1243)*X(168))/(JVS(1238))
  X(162) = (X(162)-JVS(1204)*X(163)-JVS(1205)*X(164)-JVS(1206)*X(165)-JVS(1207)*X(166)-JVS(1208)*X(167)-JVS(1209)&
             &*X(168))/(JVS(1203))
  X(161) = (X(161)-JVS(1178)*X(162)-JVS(1179)*X(163)-JVS(1180)*X(164)-JVS(1181)*X(165)-JVS(1182)*X(166)-JVS(1183)*X(167)&
             &-JVS(1184)*X(168))/(JVS(1177))
  X(160) = (X(160)-JVS(1126)*X(161)-JVS(1127)*X(162)-JVS(1128)*X(163)-JVS(1129)*X(164)-JVS(1130)*X(165)-JVS(1131)*X(166)&
             &-JVS(1132)*X(167)-JVS(1133)*X(168))/(JVS(1125))
  X(159) = (X(159)-JVS(1100)*X(160)-JVS(1101)*X(161)-JVS(1102)*X(162)-JVS(1103)*X(163)-JVS(1104)*X(164)-JVS(1105)*X(165)&
             &-JVS(1106)*X(166)-JVS(1107)*X(167)-JVS(1108)*X(168))/(JVS(1099))
  X(158) = (X(158)-JVS(1079)*X(159)-JVS(1080)*X(160)-JVS(1081)*X(161)-JVS(1082)*X(162)-JVS(1083)*X(163)-JVS(1084)*X(164)&
             &-JVS(1085)*X(165)-JVS(1086)*X(166)-JVS(1087)*X(167)-JVS(1088)*X(168))/(JVS(1078))
  X(157) = (X(157)-JVS(1046)*X(158)-JVS(1047)*X(159)-JVS(1048)*X(160)-JVS(1049)*X(161)-JVS(1050)*X(162)-JVS(1051)*X(163)&
             &-JVS(1052)*X(164)-JVS(1053)*X(165)-JVS(1054)*X(166)-JVS(1055)*X(167)-JVS(1056)*X(168))/(JVS(1045))
  X(156) = (X(156)-JVS(1007)*X(157)-JVS(1008)*X(159)-JVS(1009)*X(160)-JVS(1010)*X(161)-JVS(1011)*X(162)-JVS(1012)*X(163)&
             &-JVS(1013)*X(164)-JVS(1014)*X(165)-JVS(1015)*X(167))/(JVS(1006))
  X(155) = (X(155)-JVS(982)*X(156)-JVS(983)*X(157)-JVS(984)*X(158)-JVS(985)*X(160)-JVS(986)*X(162)-JVS(987)*X(163)&
             &-JVS(988)*X(164)-JVS(989)*X(165)-JVS(990)*X(166)-JVS(991)*X(167)-JVS(992)*X(168))/(JVS(981))
  X(154) = (X(154)-JVS(960)*X(155)-JVS(961)*X(156)-JVS(962)*X(157)-JVS(963)*X(158)-JVS(964)*X(160)-JVS(965)*X(165)&
             &-JVS(966)*X(166)-JVS(967)*X(167)-JVS(968)*X(168))/(JVS(959))
  X(153) = (X(153)-JVS(938)*X(154)-JVS(939)*X(155)-JVS(940)*X(156)-JVS(941)*X(157)-JVS(942)*X(158)-JVS(943)*X(160)&
             &-JVS(944)*X(161)-JVS(945)*X(165)-JVS(946)*X(166)-JVS(947)*X(167))/(JVS(937))
  X(152) = (X(152)-JVS(913)*X(156)-JVS(914)*X(157)-JVS(915)*X(160)-JVS(916)*X(165)-JVS(917)*X(167))/(JVS(912))
  X(151) = (X(151)-JVS(895)*X(152)-JVS(896)*X(156)-JVS(897)*X(157)-JVS(898)*X(158)-JVS(899)*X(160)-JVS(900)*X(165)&
             &-JVS(901)*X(167))/(JVS(894))
  X(150) = (X(150)-JVS(871)*X(151)-JVS(872)*X(152)-JVS(873)*X(154)-JVS(874)*X(155)-JVS(875)*X(156)-JVS(876)*X(157)&
             &-JVS(877)*X(158)-JVS(878)*X(159)-JVS(879)*X(160)-JVS(880)*X(161)-JVS(881)*X(162)-JVS(882)*X(163)-JVS(883)&
             &*X(164)-JVS(884)*X(165)-JVS(885)*X(166)-JVS(886)*X(167)-JVS(887)*X(168))/(JVS(870))
  X(149) = (X(149)-JVS(837)*X(151)-JVS(838)*X(152)-JVS(839)*X(153)-JVS(840)*X(154)-JVS(841)*X(155)-JVS(842)*X(156)&
             &-JVS(843)*X(157)-JVS(844)*X(159)-JVS(845)*X(160)-JVS(846)*X(162)-JVS(847)*X(163)-JVS(848)*X(164)-JVS(849)&
             &*X(165))/(JVS(836))
  X(148) = (X(148)-JVS(822)*X(152)-JVS(823)*X(156)-JVS(824)*X(157)-JVS(825)*X(165))/(JVS(821))
  X(147) = (X(147)-JVS(815)*X(152)-JVS(816)*X(156)-JVS(817)*X(157)-JVS(818)*X(165))/(JVS(814))
  X(146) = (X(146)-JVS(809)*X(147)-JVS(810)*X(152)-JVS(811)*X(156)-JVS(812)*X(157)-JVS(813)*X(165))/(JVS(808))
  X(145) = (X(145)-JVS(795)*X(146)-JVS(796)*X(152)-JVS(797)*X(156)-JVS(798)*X(157)-JVS(799)*X(159)-JVS(800)*X(160)&
             &-JVS(801)*X(161)-JVS(802)*X(162)-JVS(803)*X(163)-JVS(804)*X(164)-JVS(805)*X(165)-JVS(806)*X(167))/(JVS(794))
  X(144) = (X(144)-JVS(779)*X(152)-JVS(780)*X(156)-JVS(781)*X(157)-JVS(782)*X(165))/(JVS(778))
  X(143) = (X(143)-JVS(774)*X(152)-JVS(775)*X(156)-JVS(776)*X(157)-JVS(777)*X(165))/(JVS(773))
  X(142) = (X(142)-JVS(768)*X(147)-JVS(769)*X(152)-JVS(770)*X(156)-JVS(771)*X(157)-JVS(772)*X(165))/(JVS(767))
  X(141) = (X(141)-JVS(745)*X(143)-JVS(746)*X(144)-JVS(747)*X(147)-JVS(748)*X(148)-JVS(749)*X(151)-JVS(750)*X(152)&
             &-JVS(751)*X(154)-JVS(752)*X(155)-JVS(753)*X(156)-JVS(754)*X(157)-JVS(755)*X(158)-JVS(756)*X(159)-JVS(757)&
             &*X(160)-JVS(758)*X(161)-JVS(759)*X(162)-JVS(760)*X(163)-JVS(761)*X(164)-JVS(762)*X(165)-JVS(763)*X(166)&
             &-JVS(764)*X(167)-JVS(765)*X(168))/(JVS(744))
  X(140) = (X(140)-JVS(732)*X(152)-JVS(733)*X(156)-JVS(734)*X(157)-JVS(735)*X(165))/(JVS(731))
  X(139) = (X(139)-JVS(727)*X(152)-JVS(728)*X(156)-JVS(729)*X(157)-JVS(730)*X(165))/(JVS(726))
  X(138) = (X(138)-JVS(717)*X(157)-JVS(718)*X(159)-JVS(719)*X(160)-JVS(720)*X(161)-JVS(721)*X(162)-JVS(722)*X(163)&
             &-JVS(723)*X(164)-JVS(724)*X(165)-JVS(725)*X(167))/(JVS(716))
  X(137) = (X(137)-JVS(710)*X(152)-JVS(711)*X(156)-JVS(712)*X(157)-JVS(713)*X(165))/(JVS(709))
  X(136) = (X(136)-JVS(700)*X(139)-JVS(701)*X(143)-JVS(702)*X(144)-JVS(703)*X(147)-JVS(704)*X(151)-JVS(705)*X(156)&
             &-JVS(706)*X(157)-JVS(707)*X(165)-JVS(708)*X(167))/(JVS(699))
  X(135) = (X(135)-JVS(690)*X(152)-JVS(691)*X(156)-JVS(692)*X(157)-JVS(693)*X(165))/(JVS(689))
  X(134) = (X(134)-JVS(675)*X(138)-JVS(676)*X(142)-JVS(677)*X(145)-JVS(678)*X(146)-JVS(679)*X(147)-JVS(680)*X(148)&
             &-JVS(681)*X(149)-JVS(682)*X(150)-JVS(683)*X(153)-JVS(684)*X(156)-JVS(685)*X(157)-JVS(686)*X(161)-JVS(687)&
             &*X(165)-JVS(688)*X(167))/(JVS(674))
  X(133) = (X(133)-JVS(650)*X(135)-JVS(651)*X(137)-JVS(652)*X(139)-JVS(653)*X(140)-JVS(654)*X(142)-JVS(655)*X(143)&
             &-JVS(656)*X(144)-JVS(657)*X(145)-JVS(658)*X(146)-JVS(659)*X(147)-JVS(660)*X(148)-JVS(661)*X(149)-JVS(662)&
             &*X(150)-JVS(663)*X(152)-JVS(664)*X(153)-JVS(665)*X(156)-JVS(666)*X(157)-JVS(667)*X(165))/(JVS(649))
  X(132) = (X(132)-JVS(637)*X(142)-JVS(638)*X(146)-JVS(639)*X(148)-JVS(640)*X(156)-JVS(641)*X(157)-JVS(642)*X(165))&
             &/(JVS(636))
  X(131) = (X(131)-JVS(627)*X(138)-JVS(628)*X(157)-JVS(629)*X(161)-JVS(630)*X(165))/(JVS(626))
  X(130) = (X(130)-JVS(621)*X(147)-JVS(622)*X(156)-JVS(623)*X(157)-JVS(624)*X(165))/(JVS(620))
  X(129) = (X(129)-JVS(614)*X(158)-JVS(615)*X(161)-JVS(616)*X(165)-JVS(617)*X(166))/(JVS(613))
  X(128) = (X(128)-JVS(609)*X(138)-JVS(610)*X(157)-JVS(611)*X(161)-JVS(612)*X(167))/(JVS(608))
  X(127) = (X(127)-JVS(605)*X(156)-JVS(606)*X(165))/(JVS(604))
  X(126) = (X(126)-JVS(600)*X(156)-JVS(601)*X(165))/(JVS(599))
  X(125) = (X(125)-JVS(598)*X(165))/(JVS(597))
  X(124) = (X(124)-JVS(595)*X(157)-JVS(596)*X(165))/(JVS(594))
  X(123) = (X(123)-JVS(591)*X(165))/(JVS(590))
  X(122) = (X(122)-JVS(587)*X(165))/(JVS(586))
  X(121) = (X(121)-JVS(583)*X(165))/(JVS(582))
  X(120) = (X(120)-JVS(578)*X(158)-JVS(579)*X(165)-JVS(580)*X(166)-JVS(581)*X(168))/(JVS(577))
  X(119) = (X(119)-JVS(574)*X(161)-JVS(575)*X(165)-JVS(576)*X(167))/(JVS(573))
  X(118) = (X(118)-JVS(569)*X(128)-JVS(570)*X(157)-JVS(571)*X(161)-JVS(572)*X(167))/(JVS(568))
  X(117) = (X(117)-JVS(565)*X(150)-JVS(566)*X(160)-JVS(567)*X(161))/(JVS(564))
  X(116) = (X(116)-JVS(561)*X(161)-JVS(562)*X(165)-JVS(563)*X(168))/(JVS(560))
  X(115) = (X(115)-JVS(559)*X(165))/(JVS(558))
  X(114) = (X(114)-JVS(556)*X(165)-JVS(557)*X(167))/(JVS(555))
  X(113) = (X(113)-JVS(553)*X(165))/(JVS(552))
  X(112) = (X(112)-JVS(551)*X(165))/(JVS(550))
  X(111) = (X(111)-JVS(548)*X(160)-JVS(549)*X(165))/(JVS(547))
  X(110) = (X(110)-JVS(545)*X(157)-JVS(546)*X(167))/(JVS(544))
  X(109) = (X(109)-JVS(539)*X(115)-JVS(540)*X(143)-JVS(541)*X(144)-JVS(542)*X(156)-JVS(543)*X(165))/(JVS(538))
  X(108) = (X(108)-JVS(537)*X(165))/(JVS(536))
  X(107) = (X(107)-JVS(535)*X(165))/(JVS(534))
  X(106) = (X(106)-JVS(532)*X(161)-JVS(533)*X(165))/(JVS(531))
  X(105) = (X(105)-JVS(529)*X(164)-JVS(530)*X(167))/(JVS(528))
  X(104) = (X(104)-JVS(526)*X(163)-JVS(527)*X(167))/(JVS(525))
  X(103) = (X(103)-JVS(523)*X(159)-JVS(524)*X(167))/(JVS(522))
  X(102) = (X(102)-JVS(520)*X(162)-JVS(521)*X(167))/(JVS(519))
  X(101) = (X(101)-JVS(518)*X(165))/(JVS(517))
  X(100) = (X(100)-JVS(516)*X(165))/(JVS(515))
  X(99) = (X(99)-JVS(514)*X(165))/(JVS(513))
  X(98) = (X(98)-JVS(512)*X(165))/(JVS(511))
  X(97) = (X(97)-JVS(508)*X(98)-JVS(509)*X(165))/(JVS(507))
  X(96) = (X(96)-JVS(504)*X(97)-JVS(505)*X(165))/(JVS(503))
  X(95) = (X(95)-JVS(500)*X(96)-JVS(501)*X(165))/(JVS(499))
  X(94) = (X(94)-JVS(496)*X(95)-JVS(497)*X(165))/(JVS(495))
  X(93) = (X(93)-JVS(492)*X(94)-JVS(493)*X(165))/(JVS(491))
  X(92) = (X(92)-JVS(488)*X(93)-JVS(489)*X(165))/(JVS(487))
  X(91) = (X(91)-JVS(484)*X(92)-JVS(485)*X(165))/(JVS(483))
  X(90) = (X(90)-JVS(481)*X(156))/(JVS(480))
  X(89) = (X(89)-JVS(479)*X(165))/(JVS(478))
  X(88) = (X(88)-JVS(477)*X(165))/(JVS(476))
  X(87) = (X(87)-JVS(475)*X(165))/(JVS(474))
  X(86) = (X(86)-JVS(473)*X(165))/(JVS(472))
  X(85) = (X(85)-JVS(471)*X(165))/(JVS(470))
  X(84) = (X(84)-JVS(469)*X(165))/(JVS(468))
  X(83) = (X(83)-JVS(467)*X(165))/(JVS(466))
  X(82) = (X(82)-JVS(465)*X(165))/(JVS(464))
  X(81) = (X(81)-JVS(463)*X(165))/(JVS(462))
  X(80) = (X(80)-JVS(459)*X(81)-JVS(460)*X(165))/(JVS(458))
  X(79) = (X(79)-JVS(455)*X(80)-JVS(456)*X(165))/(JVS(454))
  X(78) = (X(78)-JVS(451)*X(79)-JVS(452)*X(165))/(JVS(450))
  X(77) = (X(77)-JVS(447)*X(78)-JVS(448)*X(165))/(JVS(446))
  X(76) = (X(76)-JVS(443)*X(77)-JVS(444)*X(165))/(JVS(442))
  X(75) = (X(75)-JVS(439)*X(76)-JVS(440)*X(165))/(JVS(438))
  X(74) = (X(74)-JVS(435)*X(75)-JVS(436)*X(165))/(JVS(434))
  X(73) = (X(73)-JVS(432)*X(165))/(JVS(431))
  X(72) = (X(72)-JVS(430)*X(165))/(JVS(429))
  X(71) = (X(71)-JVS(428)*X(165))/(JVS(427))
  X(70) = (X(70)-JVS(426)*X(165))/(JVS(425))
  X(69) = (X(69)-JVS(424)*X(165))/(JVS(423))
  X(68) = (X(68)-JVS(422)*X(165))/(JVS(421))
  X(67) = (X(67)-JVS(420)*X(165))/(JVS(419))
  X(66) = (X(66)-JVS(418)*X(165))/(JVS(417))
  X(65) = (X(65)-JVS(416)*X(165))/(JVS(415))
  X(64) = (X(64)-JVS(412)*X(65)-JVS(413)*X(89)-JVS(414)*X(165))/(JVS(411))
  X(63) = (X(63)-JVS(410)*X(165))/(JVS(409))
  X(62) = (X(62)-JVS(404)*X(63)-JVS(405)*X(64)-JVS(406)*X(88)-JVS(407)*X(98)-JVS(408)*X(165))/(JVS(403))
  X(61) = (X(61)-JVS(402)*X(165))/(JVS(401))
  X(60) = (X(60)-JVS(396)*X(61)-JVS(397)*X(62)-JVS(398)*X(87)-JVS(399)*X(97)-JVS(400)*X(165))/(JVS(395))
  X(59) = (X(59)-JVS(394)*X(165))/(JVS(393))
  X(58) = (X(58)-JVS(388)*X(59)-JVS(389)*X(60)-JVS(390)*X(86)-JVS(391)*X(96)-JVS(392)*X(165))/(JVS(387))
  X(57) = (X(57)-JVS(386)*X(165))/(JVS(385))
  X(56) = (X(56)-JVS(380)*X(57)-JVS(381)*X(58)-JVS(382)*X(85)-JVS(383)*X(95)-JVS(384)*X(165))/(JVS(379))
  X(55) = (X(55)-JVS(378)*X(165))/(JVS(377))
  X(54) = (X(54)-JVS(372)*X(55)-JVS(373)*X(56)-JVS(374)*X(84)-JVS(375)*X(94)-JVS(376)*X(165))/(JVS(371))
  X(53) = (X(53)-JVS(370)*X(165))/(JVS(369))
  X(52) = (X(52)-JVS(364)*X(53)-JVS(365)*X(54)-JVS(366)*X(83)-JVS(367)*X(93)-JVS(368)*X(165))/(JVS(363))
  X(51) = (X(51)-JVS(362)*X(165))/(JVS(361))
  X(50) = (X(50)-JVS(356)*X(51)-JVS(357)*X(52)-JVS(358)*X(82)-JVS(359)*X(92)-JVS(360)*X(165))/(JVS(355))
  X(49) = (X(49)-JVS(354)*X(165))/(JVS(353))
  X(48) = (X(48)-JVS(350)*X(49)-JVS(351)*X(73)-JVS(352)*X(165))/(JVS(349))
  X(47) = (X(47)-JVS(348)*X(165))/(JVS(347))
  X(46) = (X(46)-JVS(342)*X(47)-JVS(343)*X(48)-JVS(344)*X(72)-JVS(345)*X(81)-JVS(346)*X(165))/(JVS(341))
  X(45) = (X(45)-JVS(340)*X(165))/(JVS(339))
  X(44) = (X(44)-JVS(334)*X(45)-JVS(335)*X(46)-JVS(336)*X(71)-JVS(337)*X(80)-JVS(338)*X(165))/(JVS(333))
  X(43) = (X(43)-JVS(332)*X(165))/(JVS(331))
  X(42) = (X(42)-JVS(326)*X(43)-JVS(327)*X(44)-JVS(328)*X(70)-JVS(329)*X(79)-JVS(330)*X(165))/(JVS(325))
  X(41) = (X(41)-JVS(324)*X(165))/(JVS(323))
  X(40) = (X(40)-JVS(318)*X(41)-JVS(319)*X(42)-JVS(320)*X(69)-JVS(321)*X(78)-JVS(322)*X(165))/(JVS(317))
  X(39) = (X(39)-JVS(316)*X(165))/(JVS(315))
  X(38) = (X(38)-JVS(310)*X(39)-JVS(311)*X(40)-JVS(312)*X(68)-JVS(313)*X(77)-JVS(314)*X(165))/(JVS(309))
  X(37) = (X(37)-JVS(308)*X(165))/(JVS(307))
  X(36) = (X(36)-JVS(302)*X(37)-JVS(303)*X(38)-JVS(304)*X(67)-JVS(305)*X(76)-JVS(306)*X(165))/(JVS(301))
  X(35) = (X(35)-JVS(300)*X(165))/(JVS(299))
  X(34) = (X(34)-JVS(294)*X(35)-JVS(295)*X(36)-JVS(296)*X(66)-JVS(297)*X(75)-JVS(298)*X(165))/(JVS(293))
  X(33) = (X(33)-JVS(228)*X(34)-JVS(229)*X(35)-JVS(230)*X(36)-JVS(231)*X(37)-JVS(232)*X(38)-JVS(233)*X(39)-JVS(234)&
            &*X(40)-JVS(235)*X(41)-JVS(236)*X(42)-JVS(237)*X(43)-JVS(238)*X(44)-JVS(239)*X(45)-JVS(240)*X(46)-JVS(241)*X(47)&
            &-JVS(242)*X(48)-JVS(243)*X(49)-JVS(244)*X(50)-JVS(245)*X(51)-JVS(246)*X(52)-JVS(247)*X(53)-JVS(248)*X(54)&
            &-JVS(249)*X(55)-JVS(250)*X(56)-JVS(251)*X(57)-JVS(252)*X(58)-JVS(253)*X(59)-JVS(254)*X(60)-JVS(255)*X(61)&
            &-JVS(256)*X(62)-JVS(257)*X(63)-JVS(258)*X(64)-JVS(259)*X(65)-JVS(260)*X(66)-JVS(261)*X(67)-JVS(262)*X(68)&
            &-JVS(263)*X(69)-JVS(264)*X(70)-JVS(265)*X(71)-JVS(266)*X(72)-JVS(267)*X(73)-JVS(268)*X(74)-JVS(269)*X(75)&
            &-JVS(270)*X(76)-JVS(271)*X(77)-JVS(272)*X(78)-JVS(273)*X(79)-JVS(274)*X(80)-JVS(275)*X(81)-JVS(276)*X(82)&
            &-JVS(277)*X(83)-JVS(278)*X(84)-JVS(279)*X(85)-JVS(280)*X(86)-JVS(281)*X(87)-JVS(282)*X(88)-JVS(283)*X(89)&
            &-JVS(284)*X(91)-JVS(285)*X(92)-JVS(286)*X(93)-JVS(287)*X(94)-JVS(288)*X(95)-JVS(289)*X(96)-JVS(290)*X(97)&
            &-JVS(291)*X(98)-JVS(292)*X(165))/(JVS(227))
  X(32) = (X(32)-JVS(175)*X(99)-JVS(176)*X(100)-JVS(177)*X(101)-JVS(178)*X(107)-JVS(179)*X(108)-JVS(180)*X(111)-JVS(181)&
            &*X(112)-JVS(182)*X(113)-JVS(183)*X(115)-JVS(184)*X(116)-JVS(185)*X(120)-JVS(186)*X(121)-JVS(187)*X(122)&
            &-JVS(188)*X(123)-JVS(189)*X(124)-JVS(190)*X(125)-JVS(191)*X(126)-JVS(192)*X(127)-JVS(193)*X(129)-JVS(194)&
            &*X(130)-JVS(195)*X(131)-JVS(196)*X(132)-JVS(197)*X(133)-JVS(198)*X(134)-JVS(199)*X(135)-JVS(200)*X(136)&
            &-JVS(201)*X(137)-JVS(202)*X(139)-JVS(203)*X(140)-JVS(204)*X(142)-JVS(205)*X(143)-JVS(206)*X(144)-JVS(207)&
            &*X(145)-JVS(208)*X(146)-JVS(209)*X(147)-JVS(210)*X(148)-JVS(211)*X(149)-JVS(212)*X(150)-JVS(213)*X(151)&
            &-JVS(214)*X(153)-JVS(215)*X(154)-JVS(216)*X(155)-JVS(217)*X(156)-JVS(218)*X(157)-JVS(219)*X(160)-JVS(220)&
            &*X(161)-JVS(221)*X(165)-JVS(222)*X(167))/(JVS(174))
  X(31) = (X(31)-JVS(171)*X(143)-JVS(172)*X(144)-JVS(173)*X(165))/(JVS(170))
  X(30) = (X(30)-JVS(166)*X(140)-JVS(167)*X(143)-JVS(168)*X(144)-JVS(169)*X(165))/(JVS(165))
  X(29) = (X(29)-JVS(161)*X(140)-JVS(162)*X(143)-JVS(163)*X(144)-JVS(164)*X(165))/(JVS(160))
  X(28) = (X(28)-JVS(156)*X(140)-JVS(157)*X(143)-JVS(158)*X(144)-JVS(159)*X(165))/(JVS(155))
  X(27) = (X(27)-JVS(153)*X(144)-JVS(154)*X(165))/(JVS(152))
  X(26) = (X(26)-JVS(149)*X(140)-JVS(150)*X(144)-JVS(151)*X(165))/(JVS(148))
  X(25) = (X(25)-JVS(145)*X(140)-JVS(146)*X(144)-JVS(147)*X(165))/(JVS(144))
  X(24) = (X(24)-JVS(141)*X(140)-JVS(142)*X(144)-JVS(143)*X(165))/(JVS(140))
  X(23) = (X(23)-JVS(135)*X(115)-JVS(136)*X(121)-JVS(137)*X(139)-JVS(138)*X(147)-JVS(139)*X(165))/(JVS(134))
  X(22) = (X(22)-JVS(129)*X(115)-JVS(130)*X(121)-JVS(131)*X(139)-JVS(132)*X(147)-JVS(133)*X(165))/(JVS(128))
  X(21) = (X(21)-JVS(121)*X(112)-JVS(122)*X(115)-JVS(123)*X(121)-JVS(124)*X(125)-JVS(125)*X(139)-JVS(126)*X(147)&
            &-JVS(127)*X(165))/(JVS(120))
  X(20) = (X(20)-JVS(115)*X(115)-JVS(116)*X(121)-JVS(117)*X(139)-JVS(118)*X(147)-JVS(119)*X(165))/(JVS(114))
  X(19) = (X(19)-JVS(109)*X(115)-JVS(110)*X(121)-JVS(111)*X(139)-JVS(112)*X(147)-JVS(113)*X(165))/(JVS(108))
  X(18) = (X(18)-JVS(103)*X(115)-JVS(104)*X(121)-JVS(105)*X(139)-JVS(106)*X(147)-JVS(107)*X(165))/(JVS(102))
  X(17) = (X(17)-JVS(95)*X(112)-JVS(96)*X(115)-JVS(97)*X(121)-JVS(98)*X(125)-JVS(99)*X(139)-JVS(100)*X(147)-JVS(101)&
            &*X(165))/(JVS(94))
  X(16) = (X(16)-JVS(89)*X(115)-JVS(90)*X(121)-JVS(91)*X(139)-JVS(92)*X(147)-JVS(93)*X(165))/(JVS(88))
  X(15) = (X(15)-JVS(83)*X(141)-JVS(84)*X(158)-JVS(85)*X(161)-JVS(86)*X(166)-JVS(87)*X(168))/(JVS(82))
  X(14) = (X(14)-JVS(76)*X(141)-JVS(77)*X(157)-JVS(78)*X(158)-JVS(79)*X(160)-JVS(80)*X(166)-JVS(81)*X(168))/(JVS(75))
  X(13) = (X(13)-JVS(67)*X(118)-JVS(68)*X(130)-JVS(69)*X(137)-JVS(70)*X(152)-JVS(71)*X(156)-JVS(72)*X(157)-JVS(73)&
            &*X(165)-JVS(74)*X(167))/(JVS(66))
  X(12) = (X(12)-JVS(62)*X(118)-JVS(63)*X(137)-JVS(64)*X(157)-JVS(65)*X(167))/(JVS(61))
  X(11) = (X(11)-JVS(57)*X(159)-JVS(58)*X(161)-JVS(59)*X(162)-JVS(60)*X(164))/(JVS(56))
  X(10) = (X(10)-JVS(54)*X(161)-JVS(55)*X(163))/(JVS(53))
  X(9) = (X(9)-JVS(39)*X(139)-JVS(40)*X(140)-JVS(41)*X(143)-JVS(42)*X(144)-JVS(43)*X(146)-JVS(44)*X(147)-JVS(45)*X(156)&
           &-JVS(46)*X(158)-JVS(47)*X(159)-JVS(48)*X(161)-JVS(49)*X(162)-JVS(50)*X(164)-JVS(51)*X(166)-JVS(52)*X(168))&
           &/(JVS(38))
  X(8) = (X(8)-JVS(29)*X(137)-JVS(30)*X(139)-JVS(31)*X(147)-JVS(32)*X(156)-JVS(33)*X(158)-JVS(34)*X(161)-JVS(35)*X(163)&
           &-JVS(36)*X(166)-JVS(37)*X(168))/(JVS(28))
  X(7) = (X(7)-JVS(13)*X(117)-JVS(14)*X(126)-JVS(15)*X(135)-JVS(16)*X(137)-JVS(17)*X(139)-JVS(18)*X(140)-JVS(19)*X(142)&
           &-JVS(20)*X(143)-JVS(21)*X(144)-JVS(22)*X(146)-JVS(23)*X(147)-JVS(24)*X(148)-JVS(25)*X(156)-JVS(26)*X(160)&
           &-JVS(27)*X(165))/(JVS(12))
  X(6) = X(6)/JVS(11)
  X(5) = X(5)/JVS(10)
  X(4) = X(4)/JVS(9)
  X(3) = X(3)/JVS(8)
  X(2) = (X(2)-JVS(5)*X(126)-JVS(6)*X(137)-JVS(7)*X(156))/(JVS(4))
  X(1) = (X(1)-JVS(2)*X(100)-JVS(3)*X(165))/(JVS(1))
END SUBROUTINE saprc99_mosaic_4bin_vbs9_KppSolve
      SUBROUTINE saprc99_mosaic_4bin_vbs9_WCOPY(N,X,incX,Y,incY)
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
      END SUBROUTINE saprc99_mosaic_4bin_vbs9_WCOPY
      SUBROUTINE saprc99_mosaic_4bin_vbs9_WAXPY(N,Alpha,X,incX,Y,incY)
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
      END SUBROUTINE saprc99_mosaic_4bin_vbs9_WAXPY
      SUBROUTINE saprc99_mosaic_4bin_vbs9_WSCAL(N,Alpha,X,incX)
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
      END SUBROUTINE saprc99_mosaic_4bin_vbs9_WSCAL
      REAL(kind=dp) FUNCTION saprc99_mosaic_4bin_vbs9_WLAMCH( C )
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
          CALL saprc99_mosaic_4bin_vbs9_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10 Eps = Eps*2
        i = i-1
      END IF
      saprc99_mosaic_4bin_vbs9_WLAMCH = Eps
      END FUNCTION saprc99_mosaic_4bin_vbs9_WLAMCH
      SUBROUTINE saprc99_mosaic_4bin_vbs9_WLAMCH_ADD( A, B, Sum )
      REAL(kind=dp) A, B, Sum
      Sum = A + B
      END SUBROUTINE saprc99_mosaic_4bin_vbs9_WLAMCH_ADD
      SUBROUTINE saprc99_mosaic_4bin_vbs9_SET2ZERO(N,Y)
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
      END SUBROUTINE saprc99_mosaic_4bin_vbs9_SET2ZERO
      REAL(kind=dp) FUNCTION saprc99_mosaic_4bin_vbs9_WDOT (N, DX, incX, DY, incY)
      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N)
      INTEGER :: i, IX, IY, M, MP1, NS
      saprc99_mosaic_4bin_vbs9_WDOT = 0.0D0
      IF (N .LE. 0) RETURN
      IF (incX .EQ. incY) IF (incX-1) 5,20,60
    5 IX = 1
      IY = 1
      IF (incX .LT. 0) IX = (-N+1)*incX + 1
      IF (incY .LT. 0) IY = (-N+1)*incY + 1
      DO i = 1,N
        saprc99_mosaic_4bin_vbs9_WDOT = saprc99_mosaic_4bin_vbs9_WDOT + DX(IX)*DY(IY)
        IX = IX + incX
        IY = IY + incY
      END DO
      RETURN
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO i = 1,M
         saprc99_mosaic_4bin_vbs9_WDOT = saprc99_mosaic_4bin_vbs9_WDOT + DX(i)*DY(i)
      END DO
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO i = MP1,N,5
          saprc99_mosaic_4bin_vbs9_WDOT = saprc99_mosaic_4bin_vbs9_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) + &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)
      END DO
      RETURN
   60 NS = N*incX
      DO i = 1,NS,incX
        saprc99_mosaic_4bin_vbs9_WDOT = saprc99_mosaic_4bin_vbs9_WDOT + DX(i)*DY(i)
      END DO
      END FUNCTION saprc99_mosaic_4bin_vbs9_WDOT
   SUBROUTINE decomp_saprc99_mosaic_4bin_vbs9( JVS, IER )
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
   W( 100 ) = JVS( 2 )
   W( 165 ) = JVS( 3 )
  JVS( 1) = W( 1 )
  JVS( 2) = W( 100 )
  JVS( 3) = W( 165 )
  IF ( ABS( JVS( 4 )) < TINY(a) ) THEN
         IER = 2
         RETURN
  END IF
   W( 2 ) = JVS( 4 )
   W( 126 ) = JVS( 5 )
   W( 137 ) = JVS( 6 )
   W( 156 ) = JVS( 7 )
  JVS( 4) = W( 2 )
  JVS( 5) = W( 126 )
  JVS( 6) = W( 137 )
  JVS( 7) = W( 156 )
  IF ( ABS( JVS( 8 )) < TINY(a) ) THEN
         IER = 3
         RETURN
  END IF
   W( 3 ) = JVS( 8 )
  JVS( 8) = W( 3 )
  IF ( ABS( JVS( 9 )) < TINY(a) ) THEN
         IER = 4
         RETURN
  END IF
   W( 4 ) = JVS( 9 )
  JVS( 9) = W( 4 )
  IF ( ABS( JVS( 10 )) < TINY(a) ) THEN
         IER = 5
         RETURN
  END IF
   W( 5 ) = JVS( 10 )
  JVS( 10) = W( 5 )
  IF ( ABS( JVS( 11 )) < TINY(a) ) THEN
         IER = 6
         RETURN
  END IF
   W( 6 ) = JVS( 11 )
  JVS( 11) = W( 6 )
  IF ( ABS( JVS( 12 )) < TINY(a) ) THEN
         IER = 7
         RETURN
  END IF
   W( 7 ) = JVS( 12 )
   W( 117 ) = JVS( 13 )
   W( 126 ) = JVS( 14 )
   W( 135 ) = JVS( 15 )
   W( 137 ) = JVS( 16 )
   W( 139 ) = JVS( 17 )
   W( 140 ) = JVS( 18 )
   W( 142 ) = JVS( 19 )
   W( 143 ) = JVS( 20 )
   W( 144 ) = JVS( 21 )
   W( 146 ) = JVS( 22 )
   W( 147 ) = JVS( 23 )
   W( 148 ) = JVS( 24 )
   W( 156 ) = JVS( 25 )
   W( 160 ) = JVS( 26 )
   W( 165 ) = JVS( 27 )
  JVS( 12) = W( 7 )
  JVS( 13) = W( 117 )
  JVS( 14) = W( 126 )
  JVS( 15) = W( 135 )
  JVS( 16) = W( 137 )
  JVS( 17) = W( 139 )
  JVS( 18) = W( 140 )
  JVS( 19) = W( 142 )
  JVS( 20) = W( 143 )
  JVS( 21) = W( 144 )
  JVS( 22) = W( 146 )
  JVS( 23) = W( 147 )
  JVS( 24) = W( 148 )
  JVS( 25) = W( 156 )
  JVS( 26) = W( 160 )
  JVS( 27) = W( 165 )
  IF ( ABS( JVS( 28 )) < TINY(a) ) THEN
         IER = 8
         RETURN
  END IF
   W( 8 ) = JVS( 28 )
   W( 137 ) = JVS( 29 )
   W( 139 ) = JVS( 30 )
   W( 147 ) = JVS( 31 )
   W( 156 ) = JVS( 32 )
   W( 158 ) = JVS( 33 )
   W( 161 ) = JVS( 34 )
   W( 163 ) = JVS( 35 )
   W( 166 ) = JVS( 36 )
   W( 168 ) = JVS( 37 )
  JVS( 28) = W( 8 )
  JVS( 29) = W( 137 )
  JVS( 30) = W( 139 )
  JVS( 31) = W( 147 )
  JVS( 32) = W( 156 )
  JVS( 33) = W( 158 )
  JVS( 34) = W( 161 )
  JVS( 35) = W( 163 )
  JVS( 36) = W( 166 )
  JVS( 37) = W( 168 )
  IF ( ABS( JVS( 38 )) < TINY(a) ) THEN
         IER = 9
         RETURN
  END IF
   W( 9 ) = JVS( 38 )
   W( 139 ) = JVS( 39 )
   W( 140 ) = JVS( 40 )
   W( 143 ) = JVS( 41 )
   W( 144 ) = JVS( 42 )
   W( 146 ) = JVS( 43 )
   W( 147 ) = JVS( 44 )
   W( 156 ) = JVS( 45 )
   W( 158 ) = JVS( 46 )
   W( 159 ) = JVS( 47 )
   W( 161 ) = JVS( 48 )
   W( 162 ) = JVS( 49 )
   W( 164 ) = JVS( 50 )
   W( 166 ) = JVS( 51 )
   W( 168 ) = JVS( 52 )
  JVS( 38) = W( 9 )
  JVS( 39) = W( 139 )
  JVS( 40) = W( 140 )
  JVS( 41) = W( 143 )
  JVS( 42) = W( 144 )
  JVS( 43) = W( 146 )
  JVS( 44) = W( 147 )
  JVS( 45) = W( 156 )
  JVS( 46) = W( 158 )
  JVS( 47) = W( 159 )
  JVS( 48) = W( 161 )
  JVS( 49) = W( 162 )
  JVS( 50) = W( 164 )
  JVS( 51) = W( 166 )
  JVS( 52) = W( 168 )
  IF ( ABS( JVS( 53 )) < TINY(a) ) THEN
         IER = 10
         RETURN
  END IF
   W( 10 ) = JVS( 53 )
   W( 161 ) = JVS( 54 )
   W( 163 ) = JVS( 55 )
  JVS( 53) = W( 10 )
  JVS( 54) = W( 161 )
  JVS( 55) = W( 163 )
  IF ( ABS( JVS( 56 )) < TINY(a) ) THEN
         IER = 11
         RETURN
  END IF
   W( 11 ) = JVS( 56 )
   W( 159 ) = JVS( 57 )
   W( 161 ) = JVS( 58 )
   W( 162 ) = JVS( 59 )
   W( 164 ) = JVS( 60 )
  JVS( 56) = W( 11 )
  JVS( 57) = W( 159 )
  JVS( 58) = W( 161 )
  JVS( 59) = W( 162 )
  JVS( 60) = W( 164 )
  IF ( ABS( JVS( 61 )) < TINY(a) ) THEN
         IER = 12
         RETURN
  END IF
   W( 12 ) = JVS( 61 )
   W( 118 ) = JVS( 62 )
   W( 137 ) = JVS( 63 )
   W( 157 ) = JVS( 64 )
   W( 167 ) = JVS( 65 )
  JVS( 61) = W( 12 )
  JVS( 62) = W( 118 )
  JVS( 63) = W( 137 )
  JVS( 64) = W( 157 )
  JVS( 65) = W( 167 )
  IF ( ABS( JVS( 66 )) < TINY(a) ) THEN
         IER = 13
         RETURN
  END IF
   W( 13 ) = JVS( 66 )
   W( 118 ) = JVS( 67 )
   W( 130 ) = JVS( 68 )
   W( 137 ) = JVS( 69 )
   W( 152 ) = JVS( 70 )
   W( 156 ) = JVS( 71 )
   W( 157 ) = JVS( 72 )
   W( 165 ) = JVS( 73 )
   W( 167 ) = JVS( 74 )
  JVS( 66) = W( 13 )
  JVS( 67) = W( 118 )
  JVS( 68) = W( 130 )
  JVS( 69) = W( 137 )
  JVS( 70) = W( 152 )
  JVS( 71) = W( 156 )
  JVS( 72) = W( 157 )
  JVS( 73) = W( 165 )
  JVS( 74) = W( 167 )
  IF ( ABS( JVS( 75 )) < TINY(a) ) THEN
         IER = 14
         RETURN
  END IF
   W( 14 ) = JVS( 75 )
   W( 141 ) = JVS( 76 )
   W( 157 ) = JVS( 77 )
   W( 158 ) = JVS( 78 )
   W( 160 ) = JVS( 79 )
   W( 166 ) = JVS( 80 )
   W( 168 ) = JVS( 81 )
  JVS( 75) = W( 14 )
  JVS( 76) = W( 141 )
  JVS( 77) = W( 157 )
  JVS( 78) = W( 158 )
  JVS( 79) = W( 160 )
  JVS( 80) = W( 166 )
  JVS( 81) = W( 168 )
  IF ( ABS( JVS( 82 )) < TINY(a) ) THEN
         IER = 15
         RETURN
  END IF
   W( 15 ) = JVS( 82 )
   W( 141 ) = JVS( 83 )
   W( 158 ) = JVS( 84 )
   W( 161 ) = JVS( 85 )
   W( 166 ) = JVS( 86 )
   W( 168 ) = JVS( 87 )
  JVS( 82) = W( 15 )
  JVS( 83) = W( 141 )
  JVS( 84) = W( 158 )
  JVS( 85) = W( 161 )
  JVS( 86) = W( 166 )
  JVS( 87) = W( 168 )
  IF ( ABS( JVS( 88 )) < TINY(a) ) THEN
         IER = 16
         RETURN
  END IF
   W( 16 ) = JVS( 88 )
   W( 115 ) = JVS( 89 )
   W( 121 ) = JVS( 90 )
   W( 139 ) = JVS( 91 )
   W( 147 ) = JVS( 92 )
   W( 165 ) = JVS( 93 )
  JVS( 88) = W( 16 )
  JVS( 89) = W( 115 )
  JVS( 90) = W( 121 )
  JVS( 91) = W( 139 )
  JVS( 92) = W( 147 )
  JVS( 93) = W( 165 )
  IF ( ABS( JVS( 94 )) < TINY(a) ) THEN
         IER = 17
         RETURN
  END IF
   W( 17 ) = JVS( 94 )
   W( 112 ) = JVS( 95 )
   W( 115 ) = JVS( 96 )
   W( 121 ) = JVS( 97 )
   W( 125 ) = JVS( 98 )
   W( 139 ) = JVS( 99 )
   W( 147 ) = JVS( 100 )
   W( 165 ) = JVS( 101 )
  JVS( 94) = W( 17 )
  JVS( 95) = W( 112 )
  JVS( 96) = W( 115 )
  JVS( 97) = W( 121 )
  JVS( 98) = W( 125 )
  JVS( 99) = W( 139 )
  JVS( 100) = W( 147 )
  JVS( 101) = W( 165 )
  IF ( ABS( JVS( 102 )) < TINY(a) ) THEN
         IER = 18
         RETURN
  END IF
   W( 18 ) = JVS( 102 )
   W( 115 ) = JVS( 103 )
   W( 121 ) = JVS( 104 )
   W( 139 ) = JVS( 105 )
   W( 147 ) = JVS( 106 )
   W( 165 ) = JVS( 107 )
  JVS( 102) = W( 18 )
  JVS( 103) = W( 115 )
  JVS( 104) = W( 121 )
  JVS( 105) = W( 139 )
  JVS( 106) = W( 147 )
  JVS( 107) = W( 165 )
  IF ( ABS( JVS( 108 )) < TINY(a) ) THEN
         IER = 19
         RETURN
  END IF
   W( 19 ) = JVS( 108 )
   W( 115 ) = JVS( 109 )
   W( 121 ) = JVS( 110 )
   W( 139 ) = JVS( 111 )
   W( 147 ) = JVS( 112 )
   W( 165 ) = JVS( 113 )
  JVS( 108) = W( 19 )
  JVS( 109) = W( 115 )
  JVS( 110) = W( 121 )
  JVS( 111) = W( 139 )
  JVS( 112) = W( 147 )
  JVS( 113) = W( 165 )
  IF ( ABS( JVS( 114 )) < TINY(a) ) THEN
         IER = 20
         RETURN
  END IF
   W( 20 ) = JVS( 114 )
   W( 115 ) = JVS( 115 )
   W( 121 ) = JVS( 116 )
   W( 139 ) = JVS( 117 )
   W( 147 ) = JVS( 118 )
   W( 165 ) = JVS( 119 )
  JVS( 114) = W( 20 )
  JVS( 115) = W( 115 )
  JVS( 116) = W( 121 )
  JVS( 117) = W( 139 )
  JVS( 118) = W( 147 )
  JVS( 119) = W( 165 )
  IF ( ABS( JVS( 120 )) < TINY(a) ) THEN
         IER = 21
         RETURN
  END IF
   W( 21 ) = JVS( 120 )
   W( 112 ) = JVS( 121 )
   W( 115 ) = JVS( 122 )
   W( 121 ) = JVS( 123 )
   W( 125 ) = JVS( 124 )
   W( 139 ) = JVS( 125 )
   W( 147 ) = JVS( 126 )
   W( 165 ) = JVS( 127 )
  JVS( 120) = W( 21 )
  JVS( 121) = W( 112 )
  JVS( 122) = W( 115 )
  JVS( 123) = W( 121 )
  JVS( 124) = W( 125 )
  JVS( 125) = W( 139 )
  JVS( 126) = W( 147 )
  JVS( 127) = W( 165 )
  IF ( ABS( JVS( 128 )) < TINY(a) ) THEN
         IER = 22
         RETURN
  END IF
   W( 22 ) = JVS( 128 )
   W( 115 ) = JVS( 129 )
   W( 121 ) = JVS( 130 )
   W( 139 ) = JVS( 131 )
   W( 147 ) = JVS( 132 )
   W( 165 ) = JVS( 133 )
  JVS( 128) = W( 22 )
  JVS( 129) = W( 115 )
  JVS( 130) = W( 121 )
  JVS( 131) = W( 139 )
  JVS( 132) = W( 147 )
  JVS( 133) = W( 165 )
  IF ( ABS( JVS( 134 )) < TINY(a) ) THEN
         IER = 23
         RETURN
  END IF
   W( 23 ) = JVS( 134 )
   W( 115 ) = JVS( 135 )
   W( 121 ) = JVS( 136 )
   W( 139 ) = JVS( 137 )
   W( 147 ) = JVS( 138 )
   W( 165 ) = JVS( 139 )
  JVS( 134) = W( 23 )
  JVS( 135) = W( 115 )
  JVS( 136) = W( 121 )
  JVS( 137) = W( 139 )
  JVS( 138) = W( 147 )
  JVS( 139) = W( 165 )
  IF ( ABS( JVS( 140 )) < TINY(a) ) THEN
         IER = 24
         RETURN
  END IF
   W( 24 ) = JVS( 140 )
   W( 140 ) = JVS( 141 )
   W( 144 ) = JVS( 142 )
   W( 165 ) = JVS( 143 )
  JVS( 140) = W( 24 )
  JVS( 141) = W( 140 )
  JVS( 142) = W( 144 )
  JVS( 143) = W( 165 )
  IF ( ABS( JVS( 144 )) < TINY(a) ) THEN
         IER = 25
         RETURN
  END IF
   W( 25 ) = JVS( 144 )
   W( 140 ) = JVS( 145 )
   W( 144 ) = JVS( 146 )
   W( 165 ) = JVS( 147 )
  JVS( 144) = W( 25 )
  JVS( 145) = W( 140 )
  JVS( 146) = W( 144 )
  JVS( 147) = W( 165 )
  IF ( ABS( JVS( 148 )) < TINY(a) ) THEN
         IER = 26
         RETURN
  END IF
   W( 26 ) = JVS( 148 )
   W( 140 ) = JVS( 149 )
   W( 144 ) = JVS( 150 )
   W( 165 ) = JVS( 151 )
  JVS( 148) = W( 26 )
  JVS( 149) = W( 140 )
  JVS( 150) = W( 144 )
  JVS( 151) = W( 165 )
  IF ( ABS( JVS( 152 )) < TINY(a) ) THEN
         IER = 27
         RETURN
  END IF
   W( 27 ) = JVS( 152 )
   W( 144 ) = JVS( 153 )
   W( 165 ) = JVS( 154 )
  JVS( 152) = W( 27 )
  JVS( 153) = W( 144 )
  JVS( 154) = W( 165 )
  IF ( ABS( JVS( 155 )) < TINY(a) ) THEN
         IER = 28
         RETURN
  END IF
   W( 28 ) = JVS( 155 )
   W( 140 ) = JVS( 156 )
   W( 143 ) = JVS( 157 )
   W( 144 ) = JVS( 158 )
   W( 165 ) = JVS( 159 )
  JVS( 155) = W( 28 )
  JVS( 156) = W( 140 )
  JVS( 157) = W( 143 )
  JVS( 158) = W( 144 )
  JVS( 159) = W( 165 )
  IF ( ABS( JVS( 160 )) < TINY(a) ) THEN
         IER = 29
         RETURN
  END IF
   W( 29 ) = JVS( 160 )
   W( 140 ) = JVS( 161 )
   W( 143 ) = JVS( 162 )
   W( 144 ) = JVS( 163 )
   W( 165 ) = JVS( 164 )
  JVS( 160) = W( 29 )
  JVS( 161) = W( 140 )
  JVS( 162) = W( 143 )
  JVS( 163) = W( 144 )
  JVS( 164) = W( 165 )
  IF ( ABS( JVS( 165 )) < TINY(a) ) THEN
         IER = 30
         RETURN
  END IF
   W( 30 ) = JVS( 165 )
   W( 140 ) = JVS( 166 )
   W( 143 ) = JVS( 167 )
   W( 144 ) = JVS( 168 )
   W( 165 ) = JVS( 169 )
  JVS( 165) = W( 30 )
  JVS( 166) = W( 140 )
  JVS( 167) = W( 143 )
  JVS( 168) = W( 144 )
  JVS( 169) = W( 165 )
  IF ( ABS( JVS( 170 )) < TINY(a) ) THEN
         IER = 31
         RETURN
  END IF
   W( 31 ) = JVS( 170 )
   W( 143 ) = JVS( 171 )
   W( 144 ) = JVS( 172 )
   W( 165 ) = JVS( 173 )
  JVS( 170) = W( 31 )
  JVS( 171) = W( 143 )
  JVS( 172) = W( 144 )
  JVS( 173) = W( 165 )
  IF ( ABS( JVS( 174 )) < TINY(a) ) THEN
         IER = 32
         RETURN
  END IF
   W( 32 ) = JVS( 174 )
   W( 99 ) = JVS( 175 )
   W( 100 ) = JVS( 176 )
   W( 101 ) = JVS( 177 )
   W( 107 ) = JVS( 178 )
   W( 108 ) = JVS( 179 )
   W( 111 ) = JVS( 180 )
   W( 112 ) = JVS( 181 )
   W( 113 ) = JVS( 182 )
   W( 115 ) = JVS( 183 )
   W( 116 ) = JVS( 184 )
   W( 120 ) = JVS( 185 )
   W( 121 ) = JVS( 186 )
   W( 122 ) = JVS( 187 )
   W( 123 ) = JVS( 188 )
   W( 124 ) = JVS( 189 )
   W( 125 ) = JVS( 190 )
   W( 126 ) = JVS( 191 )
   W( 127 ) = JVS( 192 )
   W( 129 ) = JVS( 193 )
   W( 130 ) = JVS( 194 )
   W( 131 ) = JVS( 195 )
   W( 132 ) = JVS( 196 )
   W( 133 ) = JVS( 197 )
   W( 134 ) = JVS( 198 )
   W( 135 ) = JVS( 199 )
   W( 136 ) = JVS( 200 )
   W( 137 ) = JVS( 201 )
   W( 139 ) = JVS( 202 )
   W( 140 ) = JVS( 203 )
   W( 142 ) = JVS( 204 )
   W( 143 ) = JVS( 205 )
   W( 144 ) = JVS( 206 )
   W( 145 ) = JVS( 207 )
   W( 146 ) = JVS( 208 )
   W( 147 ) = JVS( 209 )
   W( 148 ) = JVS( 210 )
   W( 149 ) = JVS( 211 )
   W( 150 ) = JVS( 212 )
   W( 151 ) = JVS( 213 )
   W( 153 ) = JVS( 214 )
   W( 154 ) = JVS( 215 )
   W( 155 ) = JVS( 216 )
   W( 156 ) = JVS( 217 )
   W( 157 ) = JVS( 218 )
   W( 160 ) = JVS( 219 )
   W( 161 ) = JVS( 220 )
   W( 165 ) = JVS( 221 )
   W( 167 ) = JVS( 222 )
  JVS( 174) = W( 32 )
  JVS( 175) = W( 99 )
  JVS( 176) = W( 100 )
  JVS( 177) = W( 101 )
  JVS( 178) = W( 107 )
  JVS( 179) = W( 108 )
  JVS( 180) = W( 111 )
  JVS( 181) = W( 112 )
  JVS( 182) = W( 113 )
  JVS( 183) = W( 115 )
  JVS( 184) = W( 116 )
  JVS( 185) = W( 120 )
  JVS( 186) = W( 121 )
  JVS( 187) = W( 122 )
  JVS( 188) = W( 123 )
  JVS( 189) = W( 124 )
  JVS( 190) = W( 125 )
  JVS( 191) = W( 126 )
  JVS( 192) = W( 127 )
  JVS( 193) = W( 129 )
  JVS( 194) = W( 130 )
  JVS( 195) = W( 131 )
  JVS( 196) = W( 132 )
  JVS( 197) = W( 133 )
  JVS( 198) = W( 134 )
  JVS( 199) = W( 135 )
  JVS( 200) = W( 136 )
  JVS( 201) = W( 137 )
  JVS( 202) = W( 139 )
  JVS( 203) = W( 140 )
  JVS( 204) = W( 142 )
  JVS( 205) = W( 143 )
  JVS( 206) = W( 144 )
  JVS( 207) = W( 145 )
  JVS( 208) = W( 146 )
  JVS( 209) = W( 147 )
  JVS( 210) = W( 148 )
  JVS( 211) = W( 149 )
  JVS( 212) = W( 150 )
  JVS( 213) = W( 151 )
  JVS( 214) = W( 153 )
  JVS( 215) = W( 154 )
  JVS( 216) = W( 155 )
  JVS( 217) = W( 156 )
  JVS( 218) = W( 157 )
  JVS( 219) = W( 160 )
  JVS( 220) = W( 161 )
  JVS( 221) = W( 165 )
  JVS( 222) = W( 167 )
  IF ( ABS( JVS( 227 )) < TINY(a) ) THEN
         IER = 33
         RETURN
  END IF
   W( 3 ) = JVS( 223 )
   W( 4 ) = JVS( 224 )
   W( 5 ) = JVS( 225 )
   W( 6 ) = JVS( 226 )
   W( 33 ) = JVS( 227 )
   W( 34 ) = JVS( 228 )
   W( 35 ) = JVS( 229 )
   W( 36 ) = JVS( 230 )
   W( 37 ) = JVS( 231 )
   W( 38 ) = JVS( 232 )
   W( 39 ) = JVS( 233 )
   W( 40 ) = JVS( 234 )
   W( 41 ) = JVS( 235 )
   W( 42 ) = JVS( 236 )
   W( 43 ) = JVS( 237 )
   W( 44 ) = JVS( 238 )
   W( 45 ) = JVS( 239 )
   W( 46 ) = JVS( 240 )
   W( 47 ) = JVS( 241 )
   W( 48 ) = JVS( 242 )
   W( 49 ) = JVS( 243 )
   W( 50 ) = JVS( 244 )
   W( 51 ) = JVS( 245 )
   W( 52 ) = JVS( 246 )
   W( 53 ) = JVS( 247 )
   W( 54 ) = JVS( 248 )
   W( 55 ) = JVS( 249 )
   W( 56 ) = JVS( 250 )
   W( 57 ) = JVS( 251 )
   W( 58 ) = JVS( 252 )
   W( 59 ) = JVS( 253 )
   W( 60 ) = JVS( 254 )
   W( 61 ) = JVS( 255 )
   W( 62 ) = JVS( 256 )
   W( 63 ) = JVS( 257 )
   W( 64 ) = JVS( 258 )
   W( 65 ) = JVS( 259 )
   W( 66 ) = JVS( 260 )
   W( 67 ) = JVS( 261 )
   W( 68 ) = JVS( 262 )
   W( 69 ) = JVS( 263 )
   W( 70 ) = JVS( 264 )
   W( 71 ) = JVS( 265 )
   W( 72 ) = JVS( 266 )
   W( 73 ) = JVS( 267 )
   W( 74 ) = JVS( 268 )
   W( 75 ) = JVS( 269 )
   W( 76 ) = JVS( 270 )
   W( 77 ) = JVS( 271 )
   W( 78 ) = JVS( 272 )
   W( 79 ) = JVS( 273 )
   W( 80 ) = JVS( 274 )
   W( 81 ) = JVS( 275 )
   W( 82 ) = JVS( 276 )
   W( 83 ) = JVS( 277 )
   W( 84 ) = JVS( 278 )
   W( 85 ) = JVS( 279 )
   W( 86 ) = JVS( 280 )
   W( 87 ) = JVS( 281 )
   W( 88 ) = JVS( 282 )
   W( 89 ) = JVS( 283 )
   W( 91 ) = JVS( 284 )
   W( 92 ) = JVS( 285 )
   W( 93 ) = JVS( 286 )
   W( 94 ) = JVS( 287 )
   W( 95 ) = JVS( 288 )
   W( 96 ) = JVS( 289 )
   W( 97 ) = JVS( 290 )
   W( 98 ) = JVS( 291 )
   W( 165 ) = JVS( 292 )
  a = -W( 3 ) / JVS( 8 )
  W( 3 ) = -a
  a = -W( 4 ) / JVS( 9 )
  W( 4 ) = -a
  a = -W( 5 ) / JVS( 10 )
  W( 5 ) = -a
  a = -W( 6 ) / JVS( 11 )
  W( 6 ) = -a
  JVS( 223) = W( 3 )
  JVS( 224) = W( 4 )
  JVS( 225) = W( 5 )
  JVS( 226) = W( 6 )
  JVS( 227) = W( 33 )
  JVS( 228) = W( 34 )
  JVS( 229) = W( 35 )
  JVS( 230) = W( 36 )
  JVS( 231) = W( 37 )
  JVS( 232) = W( 38 )
  JVS( 233) = W( 39 )
  JVS( 234) = W( 40 )
  JVS( 235) = W( 41 )
  JVS( 236) = W( 42 )
  JVS( 237) = W( 43 )
  JVS( 238) = W( 44 )
  JVS( 239) = W( 45 )
  JVS( 240) = W( 46 )
  JVS( 241) = W( 47 )
  JVS( 242) = W( 48 )
  JVS( 243) = W( 49 )
  JVS( 244) = W( 50 )
  JVS( 245) = W( 51 )
  JVS( 246) = W( 52 )
  JVS( 247) = W( 53 )
  JVS( 248) = W( 54 )
  JVS( 249) = W( 55 )
  JVS( 250) = W( 56 )
  JVS( 251) = W( 57 )
  JVS( 252) = W( 58 )
  JVS( 253) = W( 59 )
  JVS( 254) = W( 60 )
  JVS( 255) = W( 61 )
  JVS( 256) = W( 62 )
  JVS( 257) = W( 63 )
  JVS( 258) = W( 64 )
  JVS( 259) = W( 65 )
  JVS( 260) = W( 66 )
  JVS( 261) = W( 67 )
  JVS( 262) = W( 68 )
  JVS( 263) = W( 69 )
  JVS( 264) = W( 70 )
  JVS( 265) = W( 71 )
  JVS( 266) = W( 72 )
  JVS( 267) = W( 73 )
  JVS( 268) = W( 74 )
  JVS( 269) = W( 75 )
  JVS( 270) = W( 76 )
  JVS( 271) = W( 77 )
  JVS( 272) = W( 78 )
  JVS( 273) = W( 79 )
  JVS( 274) = W( 80 )
  JVS( 275) = W( 81 )
  JVS( 276) = W( 82 )
  JVS( 277) = W( 83 )
  JVS( 278) = W( 84 )
  JVS( 279) = W( 85 )
  JVS( 280) = W( 86 )
  JVS( 281) = W( 87 )
  JVS( 282) = W( 88 )
  JVS( 283) = W( 89 )
  JVS( 284) = W( 91 )
  JVS( 285) = W( 92 )
  JVS( 286) = W( 93 )
  JVS( 287) = W( 94 )
  JVS( 288) = W( 95 )
  JVS( 289) = W( 96 )
  JVS( 290) = W( 97 )
  JVS( 291) = W( 98 )
  JVS( 292) = W( 165 )
  IF ( ABS( JVS( 293 )) < TINY(a) ) THEN
         IER = 34
         RETURN
  END IF
   W( 34 ) = JVS( 293 )
   W( 35 ) = JVS( 294 )
   W( 36 ) = JVS( 295 )
   W( 66 ) = JVS( 296 )
   W( 75 ) = JVS( 297 )
   W( 165 ) = JVS( 298 )
  JVS( 293) = W( 34 )
  JVS( 294) = W( 35 )
  JVS( 295) = W( 36 )
  JVS( 296) = W( 66 )
  JVS( 297) = W( 75 )
  JVS( 298) = W( 165 )
  IF ( ABS( JVS( 299 )) < TINY(a) ) THEN
         IER = 35
         RETURN
  END IF
   W( 35 ) = JVS( 299 )
   W( 165 ) = JVS( 300 )
  JVS( 299) = W( 35 )
  JVS( 300) = W( 165 )
  IF ( ABS( JVS( 301 )) < TINY(a) ) THEN
         IER = 36
         RETURN
  END IF
   W( 36 ) = JVS( 301 )
   W( 37 ) = JVS( 302 )
   W( 38 ) = JVS( 303 )
   W( 67 ) = JVS( 304 )
   W( 76 ) = JVS( 305 )
   W( 165 ) = JVS( 306 )
  JVS( 301) = W( 36 )
  JVS( 302) = W( 37 )
  JVS( 303) = W( 38 )
  JVS( 304) = W( 67 )
  JVS( 305) = W( 76 )
  JVS( 306) = W( 165 )
  IF ( ABS( JVS( 307 )) < TINY(a) ) THEN
         IER = 37
         RETURN
  END IF
   W( 37 ) = JVS( 307 )
   W( 165 ) = JVS( 308 )
  JVS( 307) = W( 37 )
  JVS( 308) = W( 165 )
  IF ( ABS( JVS( 309 )) < TINY(a) ) THEN
         IER = 38
         RETURN
  END IF
   W( 38 ) = JVS( 309 )
   W( 39 ) = JVS( 310 )
   W( 40 ) = JVS( 311 )
   W( 68 ) = JVS( 312 )
   W( 77 ) = JVS( 313 )
   W( 165 ) = JVS( 314 )
  JVS( 309) = W( 38 )
  JVS( 310) = W( 39 )
  JVS( 311) = W( 40 )
  JVS( 312) = W( 68 )
  JVS( 313) = W( 77 )
  JVS( 314) = W( 165 )
  IF ( ABS( JVS( 315 )) < TINY(a) ) THEN
         IER = 39
         RETURN
  END IF
   W( 39 ) = JVS( 315 )
   W( 165 ) = JVS( 316 )
  JVS( 315) = W( 39 )
  JVS( 316) = W( 165 )
  IF ( ABS( JVS( 317 )) < TINY(a) ) THEN
         IER = 40
         RETURN
  END IF
   W( 40 ) = JVS( 317 )
   W( 41 ) = JVS( 318 )
   W( 42 ) = JVS( 319 )
   W( 69 ) = JVS( 320 )
   W( 78 ) = JVS( 321 )
   W( 165 ) = JVS( 322 )
  JVS( 317) = W( 40 )
  JVS( 318) = W( 41 )
  JVS( 319) = W( 42 )
  JVS( 320) = W( 69 )
  JVS( 321) = W( 78 )
  JVS( 322) = W( 165 )
  IF ( ABS( JVS( 323 )) < TINY(a) ) THEN
         IER = 41
         RETURN
  END IF
   W( 41 ) = JVS( 323 )
   W( 165 ) = JVS( 324 )
  JVS( 323) = W( 41 )
  JVS( 324) = W( 165 )
  IF ( ABS( JVS( 325 )) < TINY(a) ) THEN
         IER = 42
         RETURN
  END IF
   W( 42 ) = JVS( 325 )
   W( 43 ) = JVS( 326 )
   W( 44 ) = JVS( 327 )
   W( 70 ) = JVS( 328 )
   W( 79 ) = JVS( 329 )
   W( 165 ) = JVS( 330 )
  JVS( 325) = W( 42 )
  JVS( 326) = W( 43 )
  JVS( 327) = W( 44 )
  JVS( 328) = W( 70 )
  JVS( 329) = W( 79 )
  JVS( 330) = W( 165 )
  IF ( ABS( JVS( 331 )) < TINY(a) ) THEN
         IER = 43
         RETURN
  END IF
   W( 43 ) = JVS( 331 )
   W( 165 ) = JVS( 332 )
  JVS( 331) = W( 43 )
  JVS( 332) = W( 165 )
  IF ( ABS( JVS( 333 )) < TINY(a) ) THEN
         IER = 44
         RETURN
  END IF
   W( 44 ) = JVS( 333 )
   W( 45 ) = JVS( 334 )
   W( 46 ) = JVS( 335 )
   W( 71 ) = JVS( 336 )
   W( 80 ) = JVS( 337 )
   W( 165 ) = JVS( 338 )
  JVS( 333) = W( 44 )
  JVS( 334) = W( 45 )
  JVS( 335) = W( 46 )
  JVS( 336) = W( 71 )
  JVS( 337) = W( 80 )
  JVS( 338) = W( 165 )
  IF ( ABS( JVS( 339 )) < TINY(a) ) THEN
         IER = 45
         RETURN
  END IF
   W( 45 ) = JVS( 339 )
   W( 165 ) = JVS( 340 )
  JVS( 339) = W( 45 )
  JVS( 340) = W( 165 )
  IF ( ABS( JVS( 341 )) < TINY(a) ) THEN
         IER = 46
         RETURN
  END IF
   W( 46 ) = JVS( 341 )
   W( 47 ) = JVS( 342 )
   W( 48 ) = JVS( 343 )
   W( 72 ) = JVS( 344 )
   W( 81 ) = JVS( 345 )
   W( 165 ) = JVS( 346 )
  JVS( 341) = W( 46 )
  JVS( 342) = W( 47 )
  JVS( 343) = W( 48 )
  JVS( 344) = W( 72 )
  JVS( 345) = W( 81 )
  JVS( 346) = W( 165 )
  IF ( ABS( JVS( 347 )) < TINY(a) ) THEN
         IER = 47
         RETURN
  END IF
   W( 47 ) = JVS( 347 )
   W( 165 ) = JVS( 348 )
  JVS( 347) = W( 47 )
  JVS( 348) = W( 165 )
  IF ( ABS( JVS( 349 )) < TINY(a) ) THEN
         IER = 48
         RETURN
  END IF
   W( 48 ) = JVS( 349 )
   W( 49 ) = JVS( 350 )
   W( 73 ) = JVS( 351 )
   W( 165 ) = JVS( 352 )
  JVS( 349) = W( 48 )
  JVS( 350) = W( 49 )
  JVS( 351) = W( 73 )
  JVS( 352) = W( 165 )
  IF ( ABS( JVS( 353 )) < TINY(a) ) THEN
         IER = 49
         RETURN
  END IF
   W( 49 ) = JVS( 353 )
   W( 165 ) = JVS( 354 )
  JVS( 353) = W( 49 )
  JVS( 354) = W( 165 )
  IF ( ABS( JVS( 355 )) < TINY(a) ) THEN
         IER = 50
         RETURN
  END IF
   W( 50 ) = JVS( 355 )
   W( 51 ) = JVS( 356 )
   W( 52 ) = JVS( 357 )
   W( 82 ) = JVS( 358 )
   W( 92 ) = JVS( 359 )
   W( 165 ) = JVS( 360 )
  JVS( 355) = W( 50 )
  JVS( 356) = W( 51 )
  JVS( 357) = W( 52 )
  JVS( 358) = W( 82 )
  JVS( 359) = W( 92 )
  JVS( 360) = W( 165 )
  IF ( ABS( JVS( 361 )) < TINY(a) ) THEN
         IER = 51
         RETURN
  END IF
   W( 51 ) = JVS( 361 )
   W( 165 ) = JVS( 362 )
  JVS( 361) = W( 51 )
  JVS( 362) = W( 165 )
  IF ( ABS( JVS( 363 )) < TINY(a) ) THEN
         IER = 52
         RETURN
  END IF
   W( 52 ) = JVS( 363 )
   W( 53 ) = JVS( 364 )
   W( 54 ) = JVS( 365 )
   W( 83 ) = JVS( 366 )
   W( 93 ) = JVS( 367 )
   W( 165 ) = JVS( 368 )
  JVS( 363) = W( 52 )
  JVS( 364) = W( 53 )
  JVS( 365) = W( 54 )
  JVS( 366) = W( 83 )
  JVS( 367) = W( 93 )
  JVS( 368) = W( 165 )
  IF ( ABS( JVS( 369 )) < TINY(a) ) THEN
         IER = 53
         RETURN
  END IF
   W( 53 ) = JVS( 369 )
   W( 165 ) = JVS( 370 )
  JVS( 369) = W( 53 )
  JVS( 370) = W( 165 )
  IF ( ABS( JVS( 371 )) < TINY(a) ) THEN
         IER = 54
         RETURN
  END IF
   W( 54 ) = JVS( 371 )
   W( 55 ) = JVS( 372 )
   W( 56 ) = JVS( 373 )
   W( 84 ) = JVS( 374 )
   W( 94 ) = JVS( 375 )
   W( 165 ) = JVS( 376 )
  JVS( 371) = W( 54 )
  JVS( 372) = W( 55 )
  JVS( 373) = W( 56 )
  JVS( 374) = W( 84 )
  JVS( 375) = W( 94 )
  JVS( 376) = W( 165 )
  IF ( ABS( JVS( 377 )) < TINY(a) ) THEN
         IER = 55
         RETURN
  END IF
   W( 55 ) = JVS( 377 )
   W( 165 ) = JVS( 378 )
  JVS( 377) = W( 55 )
  JVS( 378) = W( 165 )
  IF ( ABS( JVS( 379 )) < TINY(a) ) THEN
         IER = 56
         RETURN
  END IF
   W( 56 ) = JVS( 379 )
   W( 57 ) = JVS( 380 )
   W( 58 ) = JVS( 381 )
   W( 85 ) = JVS( 382 )
   W( 95 ) = JVS( 383 )
   W( 165 ) = JVS( 384 )
  JVS( 379) = W( 56 )
  JVS( 380) = W( 57 )
  JVS( 381) = W( 58 )
  JVS( 382) = W( 85 )
  JVS( 383) = W( 95 )
  JVS( 384) = W( 165 )
  IF ( ABS( JVS( 385 )) < TINY(a) ) THEN
         IER = 57
         RETURN
  END IF
   W( 57 ) = JVS( 385 )
   W( 165 ) = JVS( 386 )
  JVS( 385) = W( 57 )
  JVS( 386) = W( 165 )
  IF ( ABS( JVS( 387 )) < TINY(a) ) THEN
         IER = 58
         RETURN
  END IF
   W( 58 ) = JVS( 387 )
   W( 59 ) = JVS( 388 )
   W( 60 ) = JVS( 389 )
   W( 86 ) = JVS( 390 )
   W( 96 ) = JVS( 391 )
   W( 165 ) = JVS( 392 )
  JVS( 387) = W( 58 )
  JVS( 388) = W( 59 )
  JVS( 389) = W( 60 )
  JVS( 390) = W( 86 )
  JVS( 391) = W( 96 )
  JVS( 392) = W( 165 )
  IF ( ABS( JVS( 393 )) < TINY(a) ) THEN
         IER = 59
         RETURN
  END IF
   W( 59 ) = JVS( 393 )
   W( 165 ) = JVS( 394 )
  JVS( 393) = W( 59 )
  JVS( 394) = W( 165 )
  IF ( ABS( JVS( 395 )) < TINY(a) ) THEN
         IER = 60
         RETURN
  END IF
   W( 60 ) = JVS( 395 )
   W( 61 ) = JVS( 396 )
   W( 62 ) = JVS( 397 )
   W( 87 ) = JVS( 398 )
   W( 97 ) = JVS( 399 )
   W( 165 ) = JVS( 400 )
  JVS( 395) = W( 60 )
  JVS( 396) = W( 61 )
  JVS( 397) = W( 62 )
  JVS( 398) = W( 87 )
  JVS( 399) = W( 97 )
  JVS( 400) = W( 165 )
  IF ( ABS( JVS( 401 )) < TINY(a) ) THEN
         IER = 61
         RETURN
  END IF
   W( 61 ) = JVS( 401 )
   W( 165 ) = JVS( 402 )
  JVS( 401) = W( 61 )
  JVS( 402) = W( 165 )
  IF ( ABS( JVS( 403 )) < TINY(a) ) THEN
         IER = 62
         RETURN
  END IF
   W( 62 ) = JVS( 403 )
   W( 63 ) = JVS( 404 )
   W( 64 ) = JVS( 405 )
   W( 88 ) = JVS( 406 )
   W( 98 ) = JVS( 407 )
   W( 165 ) = JVS( 408 )
  JVS( 403) = W( 62 )
  JVS( 404) = W( 63 )
  JVS( 405) = W( 64 )
  JVS( 406) = W( 88 )
  JVS( 407) = W( 98 )
  JVS( 408) = W( 165 )
  IF ( ABS( JVS( 409 )) < TINY(a) ) THEN
         IER = 63
         RETURN
  END IF
   W( 63 ) = JVS( 409 )
   W( 165 ) = JVS( 410 )
  JVS( 409) = W( 63 )
  JVS( 410) = W( 165 )
  IF ( ABS( JVS( 411 )) < TINY(a) ) THEN
         IER = 64
         RETURN
  END IF
   W( 64 ) = JVS( 411 )
   W( 65 ) = JVS( 412 )
   W( 89 ) = JVS( 413 )
   W( 165 ) = JVS( 414 )
  JVS( 411) = W( 64 )
  JVS( 412) = W( 65 )
  JVS( 413) = W( 89 )
  JVS( 414) = W( 165 )
  IF ( ABS( JVS( 415 )) < TINY(a) ) THEN
         IER = 65
         RETURN
  END IF
   W( 65 ) = JVS( 415 )
   W( 165 ) = JVS( 416 )
  JVS( 415) = W( 65 )
  JVS( 416) = W( 165 )
  IF ( ABS( JVS( 417 )) < TINY(a) ) THEN
         IER = 66
         RETURN
  END IF
   W( 66 ) = JVS( 417 )
   W( 165 ) = JVS( 418 )
  JVS( 417) = W( 66 )
  JVS( 418) = W( 165 )
  IF ( ABS( JVS( 419 )) < TINY(a) ) THEN
         IER = 67
         RETURN
  END IF
   W( 67 ) = JVS( 419 )
   W( 165 ) = JVS( 420 )
  JVS( 419) = W( 67 )
  JVS( 420) = W( 165 )
  IF ( ABS( JVS( 421 )) < TINY(a) ) THEN
         IER = 68
         RETURN
  END IF
   W( 68 ) = JVS( 421 )
   W( 165 ) = JVS( 422 )
  JVS( 421) = W( 68 )
  JVS( 422) = W( 165 )
  IF ( ABS( JVS( 423 )) < TINY(a) ) THEN
         IER = 69
         RETURN
  END IF
   W( 69 ) = JVS( 423 )
   W( 165 ) = JVS( 424 )
  JVS( 423) = W( 69 )
  JVS( 424) = W( 165 )
  IF ( ABS( JVS( 425 )) < TINY(a) ) THEN
         IER = 70
         RETURN
  END IF
   W( 70 ) = JVS( 425 )
   W( 165 ) = JVS( 426 )
  JVS( 425) = W( 70 )
  JVS( 426) = W( 165 )
  IF ( ABS( JVS( 427 )) < TINY(a) ) THEN
         IER = 71
         RETURN
  END IF
   W( 71 ) = JVS( 427 )
   W( 165 ) = JVS( 428 )
  JVS( 427) = W( 71 )
  JVS( 428) = W( 165 )
  IF ( ABS( JVS( 429 )) < TINY(a) ) THEN
         IER = 72
         RETURN
  END IF
   W( 72 ) = JVS( 429 )
   W( 165 ) = JVS( 430 )
  JVS( 429) = W( 72 )
  JVS( 430) = W( 165 )
  IF ( ABS( JVS( 431 )) < TINY(a) ) THEN
         IER = 73
         RETURN
  END IF
   W( 73 ) = JVS( 431 )
   W( 165 ) = JVS( 432 )
  JVS( 431) = W( 73 )
  JVS( 432) = W( 165 )
  IF ( ABS( JVS( 434 )) < TINY(a) ) THEN
         IER = 74
         RETURN
  END IF
   W( 66 ) = JVS( 433 )
   W( 74 ) = JVS( 434 )
   W( 75 ) = JVS( 435 )
   W( 165 ) = JVS( 436 )
  a = -W( 66 ) / JVS( 417 )
  W( 66 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 418 )
  JVS( 433) = W( 66 )
  JVS( 434) = W( 74 )
  JVS( 435) = W( 75 )
  JVS( 436) = W( 165 )
  IF ( ABS( JVS( 438 )) < TINY(a) ) THEN
         IER = 75
         RETURN
  END IF
   W( 67 ) = JVS( 437 )
   W( 75 ) = JVS( 438 )
   W( 76 ) = JVS( 439 )
   W( 165 ) = JVS( 440 )
  a = -W( 67 ) / JVS( 419 )
  W( 67 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 420 )
  JVS( 437) = W( 67 )
  JVS( 438) = W( 75 )
  JVS( 439) = W( 76 )
  JVS( 440) = W( 165 )
  IF ( ABS( JVS( 442 )) < TINY(a) ) THEN
         IER = 76
         RETURN
  END IF
   W( 68 ) = JVS( 441 )
   W( 76 ) = JVS( 442 )
   W( 77 ) = JVS( 443 )
   W( 165 ) = JVS( 444 )
  a = -W( 68 ) / JVS( 421 )
  W( 68 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 422 )
  JVS( 441) = W( 68 )
  JVS( 442) = W( 76 )
  JVS( 443) = W( 77 )
  JVS( 444) = W( 165 )
  IF ( ABS( JVS( 446 )) < TINY(a) ) THEN
         IER = 77
         RETURN
  END IF
   W( 69 ) = JVS( 445 )
   W( 77 ) = JVS( 446 )
   W( 78 ) = JVS( 447 )
   W( 165 ) = JVS( 448 )
  a = -W( 69 ) / JVS( 423 )
  W( 69 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 424 )
  JVS( 445) = W( 69 )
  JVS( 446) = W( 77 )
  JVS( 447) = W( 78 )
  JVS( 448) = W( 165 )
  IF ( ABS( JVS( 450 )) < TINY(a) ) THEN
         IER = 78
         RETURN
  END IF
   W( 70 ) = JVS( 449 )
   W( 78 ) = JVS( 450 )
   W( 79 ) = JVS( 451 )
   W( 165 ) = JVS( 452 )
  a = -W( 70 ) / JVS( 425 )
  W( 70 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 426 )
  JVS( 449) = W( 70 )
  JVS( 450) = W( 78 )
  JVS( 451) = W( 79 )
  JVS( 452) = W( 165 )
  IF ( ABS( JVS( 454 )) < TINY(a) ) THEN
         IER = 79
         RETURN
  END IF
   W( 71 ) = JVS( 453 )
   W( 79 ) = JVS( 454 )
   W( 80 ) = JVS( 455 )
   W( 165 ) = JVS( 456 )
  a = -W( 71 ) / JVS( 427 )
  W( 71 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 428 )
  JVS( 453) = W( 71 )
  JVS( 454) = W( 79 )
  JVS( 455) = W( 80 )
  JVS( 456) = W( 165 )
  IF ( ABS( JVS( 458 )) < TINY(a) ) THEN
         IER = 80
         RETURN
  END IF
   W( 72 ) = JVS( 457 )
   W( 80 ) = JVS( 458 )
   W( 81 ) = JVS( 459 )
   W( 165 ) = JVS( 460 )
  a = -W( 72 ) / JVS( 429 )
  W( 72 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 430 )
  JVS( 457) = W( 72 )
  JVS( 458) = W( 80 )
  JVS( 459) = W( 81 )
  JVS( 460) = W( 165 )
  IF ( ABS( JVS( 462 )) < TINY(a) ) THEN
         IER = 81
         RETURN
  END IF
   W( 73 ) = JVS( 461 )
   W( 81 ) = JVS( 462 )
   W( 165 ) = JVS( 463 )
  a = -W( 73 ) / JVS( 431 )
  W( 73 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 432 )
  JVS( 461) = W( 73 )
  JVS( 462) = W( 81 )
  JVS( 463) = W( 165 )
  IF ( ABS( JVS( 464 )) < TINY(a) ) THEN
         IER = 82
         RETURN
  END IF
   W( 82 ) = JVS( 464 )
   W( 165 ) = JVS( 465 )
  JVS( 464) = W( 82 )
  JVS( 465) = W( 165 )
  IF ( ABS( JVS( 466 )) < TINY(a) ) THEN
         IER = 83
         RETURN
  END IF
   W( 83 ) = JVS( 466 )
   W( 165 ) = JVS( 467 )
  JVS( 466) = W( 83 )
  JVS( 467) = W( 165 )
  IF ( ABS( JVS( 468 )) < TINY(a) ) THEN
         IER = 84
         RETURN
  END IF
   W( 84 ) = JVS( 468 )
   W( 165 ) = JVS( 469 )
  JVS( 468) = W( 84 )
  JVS( 469) = W( 165 )
  IF ( ABS( JVS( 470 )) < TINY(a) ) THEN
         IER = 85
         RETURN
  END IF
   W( 85 ) = JVS( 470 )
   W( 165 ) = JVS( 471 )
  JVS( 470) = W( 85 )
  JVS( 471) = W( 165 )
  IF ( ABS( JVS( 472 )) < TINY(a) ) THEN
         IER = 86
         RETURN
  END IF
   W( 86 ) = JVS( 472 )
   W( 165 ) = JVS( 473 )
  JVS( 472) = W( 86 )
  JVS( 473) = W( 165 )
  IF ( ABS( JVS( 474 )) < TINY(a) ) THEN
         IER = 87
         RETURN
  END IF
   W( 87 ) = JVS( 474 )
   W( 165 ) = JVS( 475 )
  JVS( 474) = W( 87 )
  JVS( 475) = W( 165 )
  IF ( ABS( JVS( 476 )) < TINY(a) ) THEN
         IER = 88
         RETURN
  END IF
   W( 88 ) = JVS( 476 )
   W( 165 ) = JVS( 477 )
  JVS( 476) = W( 88 )
  JVS( 477) = W( 165 )
  IF ( ABS( JVS( 478 )) < TINY(a) ) THEN
         IER = 89
         RETURN
  END IF
   W( 89 ) = JVS( 478 )
   W( 165 ) = JVS( 479 )
  JVS( 478) = W( 89 )
  JVS( 479) = W( 165 )
  IF ( ABS( JVS( 480 )) < TINY(a) ) THEN
         IER = 90
         RETURN
  END IF
   W( 90 ) = JVS( 480 )
   W( 156 ) = JVS( 481 )
  JVS( 480) = W( 90 )
  JVS( 481) = W( 156 )
  IF ( ABS( JVS( 483 )) < TINY(a) ) THEN
         IER = 91
         RETURN
  END IF
   W( 82 ) = JVS( 482 )
   W( 91 ) = JVS( 483 )
   W( 92 ) = JVS( 484 )
   W( 165 ) = JVS( 485 )
  a = -W( 82 ) / JVS( 464 )
  W( 82 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 465 )
  JVS( 482) = W( 82 )
  JVS( 483) = W( 91 )
  JVS( 484) = W( 92 )
  JVS( 485) = W( 165 )
  IF ( ABS( JVS( 487 )) < TINY(a) ) THEN
         IER = 92
         RETURN
  END IF
   W( 83 ) = JVS( 486 )
   W( 92 ) = JVS( 487 )
   W( 93 ) = JVS( 488 )
   W( 165 ) = JVS( 489 )
  a = -W( 83 ) / JVS( 466 )
  W( 83 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 467 )
  JVS( 486) = W( 83 )
  JVS( 487) = W( 92 )
  JVS( 488) = W( 93 )
  JVS( 489) = W( 165 )
  IF ( ABS( JVS( 491 )) < TINY(a) ) THEN
         IER = 93
         RETURN
  END IF
   W( 84 ) = JVS( 490 )
   W( 93 ) = JVS( 491 )
   W( 94 ) = JVS( 492 )
   W( 165 ) = JVS( 493 )
  a = -W( 84 ) / JVS( 468 )
  W( 84 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 469 )
  JVS( 490) = W( 84 )
  JVS( 491) = W( 93 )
  JVS( 492) = W( 94 )
  JVS( 493) = W( 165 )
  IF ( ABS( JVS( 495 )) < TINY(a) ) THEN
         IER = 94
         RETURN
  END IF
   W( 85 ) = JVS( 494 )
   W( 94 ) = JVS( 495 )
   W( 95 ) = JVS( 496 )
   W( 165 ) = JVS( 497 )
  a = -W( 85 ) / JVS( 470 )
  W( 85 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 471 )
  JVS( 494) = W( 85 )
  JVS( 495) = W( 94 )
  JVS( 496) = W( 95 )
  JVS( 497) = W( 165 )
  IF ( ABS( JVS( 499 )) < TINY(a) ) THEN
         IER = 95
         RETURN
  END IF
   W( 86 ) = JVS( 498 )
   W( 95 ) = JVS( 499 )
   W( 96 ) = JVS( 500 )
   W( 165 ) = JVS( 501 )
  a = -W( 86 ) / JVS( 472 )
  W( 86 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 473 )
  JVS( 498) = W( 86 )
  JVS( 499) = W( 95 )
  JVS( 500) = W( 96 )
  JVS( 501) = W( 165 )
  IF ( ABS( JVS( 503 )) < TINY(a) ) THEN
         IER = 96
         RETURN
  END IF
   W( 87 ) = JVS( 502 )
   W( 96 ) = JVS( 503 )
   W( 97 ) = JVS( 504 )
   W( 165 ) = JVS( 505 )
  a = -W( 87 ) / JVS( 474 )
  W( 87 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 475 )
  JVS( 502) = W( 87 )
  JVS( 503) = W( 96 )
  JVS( 504) = W( 97 )
  JVS( 505) = W( 165 )
  IF ( ABS( JVS( 507 )) < TINY(a) ) THEN
         IER = 97
         RETURN
  END IF
   W( 88 ) = JVS( 506 )
   W( 97 ) = JVS( 507 )
   W( 98 ) = JVS( 508 )
   W( 165 ) = JVS( 509 )
  a = -W( 88 ) / JVS( 476 )
  W( 88 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 477 )
  JVS( 506) = W( 88 )
  JVS( 507) = W( 97 )
  JVS( 508) = W( 98 )
  JVS( 509) = W( 165 )
  IF ( ABS( JVS( 511 )) < TINY(a) ) THEN
         IER = 98
         RETURN
  END IF
   W( 89 ) = JVS( 510 )
   W( 98 ) = JVS( 511 )
   W( 165 ) = JVS( 512 )
  a = -W( 89 ) / JVS( 478 )
  W( 89 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 479 )
  JVS( 510) = W( 89 )
  JVS( 511) = W( 98 )
  JVS( 512) = W( 165 )
  IF ( ABS( JVS( 513 )) < TINY(a) ) THEN
         IER = 99
         RETURN
  END IF
   W( 99 ) = JVS( 513 )
   W( 165 ) = JVS( 514 )
  JVS( 513) = W( 99 )
  JVS( 514) = W( 165 )
  IF ( ABS( JVS( 515 )) < TINY(a) ) THEN
         IER = 100
         RETURN
  END IF
   W( 100 ) = JVS( 515 )
   W( 165 ) = JVS( 516 )
  JVS( 515) = W( 100 )
  JVS( 516) = W( 165 )
  IF ( ABS( JVS( 517 )) < TINY(a) ) THEN
         IER = 101
         RETURN
  END IF
   W( 101 ) = JVS( 517 )
   W( 165 ) = JVS( 518 )
  JVS( 517) = W( 101 )
  JVS( 518) = W( 165 )
  IF ( ABS( JVS( 519 )) < TINY(a) ) THEN
         IER = 102
         RETURN
  END IF
   W( 102 ) = JVS( 519 )
   W( 162 ) = JVS( 520 )
   W( 167 ) = JVS( 521 )
  JVS( 519) = W( 102 )
  JVS( 520) = W( 162 )
  JVS( 521) = W( 167 )
  IF ( ABS( JVS( 522 )) < TINY(a) ) THEN
         IER = 103
         RETURN
  END IF
   W( 103 ) = JVS( 522 )
   W( 159 ) = JVS( 523 )
   W( 167 ) = JVS( 524 )
  JVS( 522) = W( 103 )
  JVS( 523) = W( 159 )
  JVS( 524) = W( 167 )
  IF ( ABS( JVS( 525 )) < TINY(a) ) THEN
         IER = 104
         RETURN
  END IF
   W( 104 ) = JVS( 525 )
   W( 163 ) = JVS( 526 )
   W( 167 ) = JVS( 527 )
  JVS( 525) = W( 104 )
  JVS( 526) = W( 163 )
  JVS( 527) = W( 167 )
  IF ( ABS( JVS( 528 )) < TINY(a) ) THEN
         IER = 105
         RETURN
  END IF
   W( 105 ) = JVS( 528 )
   W( 164 ) = JVS( 529 )
   W( 167 ) = JVS( 530 )
  JVS( 528) = W( 105 )
  JVS( 529) = W( 164 )
  JVS( 530) = W( 167 )
  IF ( ABS( JVS( 531 )) < TINY(a) ) THEN
         IER = 106
         RETURN
  END IF
   W( 106 ) = JVS( 531 )
   W( 161 ) = JVS( 532 )
   W( 165 ) = JVS( 533 )
  JVS( 531) = W( 106 )
  JVS( 532) = W( 161 )
  JVS( 533) = W( 165 )
  IF ( ABS( JVS( 534 )) < TINY(a) ) THEN
         IER = 107
         RETURN
  END IF
   W( 107 ) = JVS( 534 )
   W( 165 ) = JVS( 535 )
  JVS( 534) = W( 107 )
  JVS( 535) = W( 165 )
  IF ( ABS( JVS( 536 )) < TINY(a) ) THEN
         IER = 108
         RETURN
  END IF
   W( 108 ) = JVS( 536 )
   W( 165 ) = JVS( 537 )
  JVS( 536) = W( 108 )
  JVS( 537) = W( 165 )
  IF ( ABS( JVS( 538 )) < TINY(a) ) THEN
         IER = 109
         RETURN
  END IF
   W( 109 ) = JVS( 538 )
   W( 115 ) = JVS( 539 )
   W( 143 ) = JVS( 540 )
   W( 144 ) = JVS( 541 )
   W( 156 ) = JVS( 542 )
   W( 165 ) = JVS( 543 )
  JVS( 538) = W( 109 )
  JVS( 539) = W( 115 )
  JVS( 540) = W( 143 )
  JVS( 541) = W( 144 )
  JVS( 542) = W( 156 )
  JVS( 543) = W( 165 )
  IF ( ABS( JVS( 544 )) < TINY(a) ) THEN
         IER = 110
         RETURN
  END IF
   W( 110 ) = JVS( 544 )
   W( 157 ) = JVS( 545 )
   W( 167 ) = JVS( 546 )
  JVS( 544) = W( 110 )
  JVS( 545) = W( 157 )
  JVS( 546) = W( 167 )
  IF ( ABS( JVS( 547 )) < TINY(a) ) THEN
         IER = 111
         RETURN
  END IF
   W( 111 ) = JVS( 547 )
   W( 160 ) = JVS( 548 )
   W( 165 ) = JVS( 549 )
  JVS( 547) = W( 111 )
  JVS( 548) = W( 160 )
  JVS( 549) = W( 165 )
  IF ( ABS( JVS( 550 )) < TINY(a) ) THEN
         IER = 112
         RETURN
  END IF
   W( 112 ) = JVS( 550 )
   W( 165 ) = JVS( 551 )
  JVS( 550) = W( 112 )
  JVS( 551) = W( 165 )
  IF ( ABS( JVS( 552 )) < TINY(a) ) THEN
         IER = 113
         RETURN
  END IF
   W( 113 ) = JVS( 552 )
   W( 165 ) = JVS( 553 )
  JVS( 552) = W( 113 )
  JVS( 553) = W( 165 )
  IF ( ABS( JVS( 555 )) < TINY(a) ) THEN
         IER = 114
         RETURN
  END IF
   W( 113 ) = JVS( 554 )
   W( 114 ) = JVS( 555 )
   W( 165 ) = JVS( 556 )
   W( 167 ) = JVS( 557 )
  a = -W( 113 ) / JVS( 552 )
  W( 113 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 553 )
  JVS( 554) = W( 113 )
  JVS( 555) = W( 114 )
  JVS( 556) = W( 165 )
  JVS( 557) = W( 167 )
  IF ( ABS( JVS( 558 )) < TINY(a) ) THEN
         IER = 115
         RETURN
  END IF
   W( 115 ) = JVS( 558 )
   W( 165 ) = JVS( 559 )
  JVS( 558) = W( 115 )
  JVS( 559) = W( 165 )
  IF ( ABS( JVS( 560 )) < TINY(a) ) THEN
         IER = 116
         RETURN
  END IF
   W( 116 ) = JVS( 560 )
   W( 161 ) = JVS( 561 )
   W( 165 ) = JVS( 562 )
   W( 168 ) = JVS( 563 )
  JVS( 560) = W( 116 )
  JVS( 561) = W( 161 )
  JVS( 562) = W( 165 )
  JVS( 563) = W( 168 )
  IF ( ABS( JVS( 564 )) < TINY(a) ) THEN
         IER = 117
         RETURN
  END IF
   W( 117 ) = JVS( 564 )
   W( 150 ) = JVS( 565 )
   W( 160 ) = JVS( 566 )
   W( 161 ) = JVS( 567 )
  JVS( 564) = W( 117 )
  JVS( 565) = W( 150 )
  JVS( 566) = W( 160 )
  JVS( 567) = W( 161 )
  IF ( ABS( JVS( 568 )) < TINY(a) ) THEN
         IER = 118
         RETURN
  END IF
   W( 118 ) = JVS( 568 )
   W( 128 ) = JVS( 569 )
   W( 157 ) = JVS( 570 )
   W( 161 ) = JVS( 571 )
   W( 167 ) = JVS( 572 )
  JVS( 568) = W( 118 )
  JVS( 569) = W( 128 )
  JVS( 570) = W( 157 )
  JVS( 571) = W( 161 )
  JVS( 572) = W( 167 )
  IF ( ABS( JVS( 573 )) < TINY(a) ) THEN
         IER = 119
         RETURN
  END IF
   W( 119 ) = JVS( 573 )
   W( 161 ) = JVS( 574 )
   W( 165 ) = JVS( 575 )
   W( 167 ) = JVS( 576 )
  JVS( 573) = W( 119 )
  JVS( 574) = W( 161 )
  JVS( 575) = W( 165 )
  JVS( 576) = W( 167 )
  IF ( ABS( JVS( 577 )) < TINY(a) ) THEN
         IER = 120
         RETURN
  END IF
   W( 120 ) = JVS( 577 )
   W( 158 ) = JVS( 578 )
   W( 165 ) = JVS( 579 )
   W( 166 ) = JVS( 580 )
   W( 168 ) = JVS( 581 )
  JVS( 577) = W( 120 )
  JVS( 578) = W( 158 )
  JVS( 579) = W( 165 )
  JVS( 580) = W( 166 )
  JVS( 581) = W( 168 )
  IF ( ABS( JVS( 582 )) < TINY(a) ) THEN
         IER = 121
         RETURN
  END IF
   W( 121 ) = JVS( 582 )
   W( 165 ) = JVS( 583 )
  JVS( 582) = W( 121 )
  JVS( 583) = W( 165 )
  IF ( ABS( JVS( 586 )) < TINY(a) ) THEN
         IER = 122
         RETURN
  END IF
   W( 115 ) = JVS( 584 )
   W( 121 ) = JVS( 585 )
   W( 122 ) = JVS( 586 )
   W( 165 ) = JVS( 587 )
  a = -W( 115 ) / JVS( 558 )
  W( 115 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 559 )
  a = -W( 121 ) / JVS( 582 )
  W( 121 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 583 )
  JVS( 584) = W( 115 )
  JVS( 585) = W( 121 )
  JVS( 586) = W( 122 )
  JVS( 587) = W( 165 )
  IF ( ABS( JVS( 590 )) < TINY(a) ) THEN
         IER = 123
         RETURN
  END IF
   W( 115 ) = JVS( 588 )
   W( 121 ) = JVS( 589 )
   W( 123 ) = JVS( 590 )
   W( 165 ) = JVS( 591 )
  a = -W( 115 ) / JVS( 558 )
  W( 115 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 559 )
  a = -W( 121 ) / JVS( 582 )
  W( 121 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 583 )
  JVS( 588) = W( 115 )
  JVS( 589) = W( 121 )
  JVS( 590) = W( 123 )
  JVS( 591) = W( 165 )
  IF ( ABS( JVS( 594 )) < TINY(a) ) THEN
         IER = 124
         RETURN
  END IF
   W( 115 ) = JVS( 592 )
   W( 121 ) = JVS( 593 )
   W( 124 ) = JVS( 594 )
   W( 157 ) = JVS( 595 )
   W( 165 ) = JVS( 596 )
  a = -W( 115 ) / JVS( 558 )
  W( 115 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 559 )
  a = -W( 121 ) / JVS( 582 )
  W( 121 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 583 )
  JVS( 592) = W( 115 )
  JVS( 593) = W( 121 )
  JVS( 594) = W( 124 )
  JVS( 595) = W( 157 )
  JVS( 596) = W( 165 )
  IF ( ABS( JVS( 597 )) < TINY(a) ) THEN
         IER = 125
         RETURN
  END IF
   W( 125 ) = JVS( 597 )
   W( 165 ) = JVS( 598 )
  JVS( 597) = W( 125 )
  JVS( 598) = W( 165 )
  IF ( ABS( JVS( 599 )) < TINY(a) ) THEN
         IER = 126
         RETURN
  END IF
   W( 126 ) = JVS( 599 )
   W( 156 ) = JVS( 600 )
   W( 165 ) = JVS( 601 )
  JVS( 599) = W( 126 )
  JVS( 600) = W( 156 )
  JVS( 601) = W( 165 )
  IF ( ABS( JVS( 604 )) < TINY(a) ) THEN
         IER = 127
         RETURN
  END IF
   W( 115 ) = JVS( 602 )
   W( 121 ) = JVS( 603 )
   W( 127 ) = JVS( 604 )
   W( 156 ) = JVS( 605 )
   W( 165 ) = JVS( 606 )
  a = -W( 115 ) / JVS( 558 )
  W( 115 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 559 )
  a = -W( 121 ) / JVS( 582 )
  W( 121 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 583 )
  JVS( 602) = W( 115 )
  JVS( 603) = W( 121 )
  JVS( 604) = W( 127 )
  JVS( 605) = W( 156 )
  JVS( 606) = W( 165 )
  IF ( ABS( JVS( 608 )) < TINY(a) ) THEN
         IER = 128
         RETURN
  END IF
   W( 118 ) = JVS( 607 )
   W( 128 ) = JVS( 608 )
   W( 138 ) = JVS( 609 )
   W( 157 ) = JVS( 610 )
   W( 161 ) = JVS( 611 )
   W( 167 ) = JVS( 612 )
  a = -W( 118 ) / JVS( 568 )
  W( 118 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 569 )
  W( 157 ) = W( 157 ) + a*JVS( 570 )
  W( 161 ) = W( 161 ) + a*JVS( 571 )
  W( 167 ) = W( 167 ) + a*JVS( 572 )
  JVS( 607) = W( 118 )
  JVS( 608) = W( 128 )
  JVS( 609) = W( 138 )
  JVS( 610) = W( 157 )
  JVS( 611) = W( 161 )
  JVS( 612) = W( 167 )
  IF ( ABS( JVS( 613 )) < TINY(a) ) THEN
         IER = 129
         RETURN
  END IF
   W( 129 ) = JVS( 613 )
   W( 158 ) = JVS( 614 )
   W( 161 ) = JVS( 615 )
   W( 165 ) = JVS( 616 )
   W( 166 ) = JVS( 617 )
  JVS( 613) = W( 129 )
  JVS( 614) = W( 158 )
  JVS( 615) = W( 161 )
  JVS( 616) = W( 165 )
  JVS( 617) = W( 166 )
  IF ( ABS( JVS( 620 )) < TINY(a) ) THEN
         IER = 130
         RETURN
  END IF
   W( 115 ) = JVS( 618 )
   W( 121 ) = JVS( 619 )
   W( 130 ) = JVS( 620 )
   W( 147 ) = JVS( 621 )
   W( 156 ) = JVS( 622 )
   W( 157 ) = JVS( 623 )
   W( 165 ) = JVS( 624 )
  a = -W( 115 ) / JVS( 558 )
  W( 115 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 559 )
  a = -W( 121 ) / JVS( 582 )
  W( 121 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 583 )
  JVS( 618) = W( 115 )
  JVS( 619) = W( 121 )
  JVS( 620) = W( 130 )
  JVS( 621) = W( 147 )
  JVS( 622) = W( 156 )
  JVS( 623) = W( 157 )
  JVS( 624) = W( 165 )
  IF ( ABS( JVS( 626 )) < TINY(a) ) THEN
         IER = 131
         RETURN
  END IF
   W( 121 ) = JVS( 625 )
   W( 131 ) = JVS( 626 )
   W( 138 ) = JVS( 627 )
   W( 157 ) = JVS( 628 )
   W( 161 ) = JVS( 629 )
   W( 165 ) = JVS( 630 )
  a = -W( 121 ) / JVS( 582 )
  W( 121 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 583 )
  JVS( 625) = W( 121 )
  JVS( 626) = W( 131 )
  JVS( 627) = W( 138 )
  JVS( 628) = W( 157 )
  JVS( 629) = W( 161 )
  JVS( 630) = W( 165 )
  IF ( ABS( JVS( 636 )) < TINY(a) ) THEN
         IER = 132
         RETURN
  END IF
   W( 115 ) = JVS( 631 )
   W( 121 ) = JVS( 632 )
   W( 122 ) = JVS( 633 )
   W( 123 ) = JVS( 634 )
   W( 124 ) = JVS( 635 )
   W( 132 ) = JVS( 636 )
   W( 142 ) = JVS( 637 )
   W( 146 ) = JVS( 638 )
   W( 148 ) = JVS( 639 )
   W( 156 ) = JVS( 640 )
   W( 157 ) = JVS( 641 )
   W( 165 ) = JVS( 642 )
  a = -W( 115 ) / JVS( 558 )
  W( 115 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 559 )
  a = -W( 121 ) / JVS( 582 )
  W( 121 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 583 )
  a = -W( 122 ) / JVS( 586 )
  W( 122 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 587 )
  a = -W( 123 ) / JVS( 590 )
  W( 123 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 591 )
  a = -W( 124 ) / JVS( 594 )
  W( 124 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 595 )
  W( 165 ) = W( 165 ) + a*JVS( 596 )
  JVS( 631) = W( 115 )
  JVS( 632) = W( 121 )
  JVS( 633) = W( 122 )
  JVS( 634) = W( 123 )
  JVS( 635) = W( 124 )
  JVS( 636) = W( 132 )
  JVS( 637) = W( 142 )
  JVS( 638) = W( 146 )
  JVS( 639) = W( 148 )
  JVS( 640) = W( 156 )
  JVS( 641) = W( 157 )
  JVS( 642) = W( 165 )
  IF ( ABS( JVS( 649 )) < TINY(a) ) THEN
         IER = 133
         RETURN
  END IF
   W( 122 ) = JVS( 643 )
   W( 123 ) = JVS( 644 )
   W( 125 ) = JVS( 645 )
   W( 126 ) = JVS( 646 )
   W( 127 ) = JVS( 647 )
   W( 132 ) = JVS( 648 )
   W( 133 ) = JVS( 649 )
   W( 135 ) = JVS( 650 )
   W( 137 ) = JVS( 651 )
   W( 139 ) = JVS( 652 )
   W( 140 ) = JVS( 653 )
   W( 142 ) = JVS( 654 )
   W( 143 ) = JVS( 655 )
   W( 144 ) = JVS( 656 )
   W( 145 ) = JVS( 657 )
   W( 146 ) = JVS( 658 )
   W( 147 ) = JVS( 659 )
   W( 148 ) = JVS( 660 )
   W( 149 ) = JVS( 661 )
   W( 150 ) = JVS( 662 )
   W( 152 ) = JVS( 663 )
   W( 153 ) = JVS( 664 )
   W( 156 ) = JVS( 665 )
   W( 157 ) = JVS( 666 )
   W( 165 ) = JVS( 667 )
  a = -W( 122 ) / JVS( 586 )
  W( 122 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 587 )
  a = -W( 123 ) / JVS( 590 )
  W( 123 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 591 )
  a = -W( 125 ) / JVS( 597 )
  W( 125 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 598 )
  a = -W( 126 ) / JVS( 599 )
  W( 126 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 600 )
  W( 165 ) = W( 165 ) + a*JVS( 601 )
  a = -W( 127 ) / JVS( 604 )
  W( 127 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 605 )
  W( 165 ) = W( 165 ) + a*JVS( 606 )
  a = -W( 132 ) / JVS( 636 )
  W( 132 ) = -a
  W( 142 ) = W( 142 ) + a*JVS( 637 )
  W( 146 ) = W( 146 ) + a*JVS( 638 )
  W( 148 ) = W( 148 ) + a*JVS( 639 )
  W( 156 ) = W( 156 ) + a*JVS( 640 )
  W( 157 ) = W( 157 ) + a*JVS( 641 )
  W( 165 ) = W( 165 ) + a*JVS( 642 )
  JVS( 643) = W( 122 )
  JVS( 644) = W( 123 )
  JVS( 645) = W( 125 )
  JVS( 646) = W( 126 )
  JVS( 647) = W( 127 )
  JVS( 648) = W( 132 )
  JVS( 649) = W( 133 )
  JVS( 650) = W( 135 )
  JVS( 651) = W( 137 )
  JVS( 652) = W( 139 )
  JVS( 653) = W( 140 )
  JVS( 654) = W( 142 )
  JVS( 655) = W( 143 )
  JVS( 656) = W( 144 )
  JVS( 657) = W( 145 )
  JVS( 658) = W( 146 )
  JVS( 659) = W( 147 )
  JVS( 660) = W( 148 )
  JVS( 661) = W( 149 )
  JVS( 662) = W( 150 )
  JVS( 663) = W( 152 )
  JVS( 664) = W( 153 )
  JVS( 665) = W( 156 )
  JVS( 666) = W( 157 )
  JVS( 667) = W( 165 )
  IF ( ABS( JVS( 674 )) < TINY(a) ) THEN
         IER = 134
         RETURN
  END IF
   W( 110 ) = JVS( 668 )
   W( 124 ) = JVS( 669 )
   W( 128 ) = JVS( 670 )
   W( 130 ) = JVS( 671 )
   W( 131 ) = JVS( 672 )
   W( 132 ) = JVS( 673 )
   W( 134 ) = JVS( 674 )
   W( 138 ) = JVS( 675 )
   W( 142 ) = JVS( 676 )
   W( 145 ) = JVS( 677 )
   W( 146 ) = JVS( 678 )
   W( 147 ) = JVS( 679 )
   W( 148 ) = JVS( 680 )
   W( 149 ) = JVS( 681 )
   W( 150 ) = JVS( 682 )
   W( 153 ) = JVS( 683 )
   W( 156 ) = JVS( 684 )
   W( 157 ) = JVS( 685 )
   W( 161 ) = JVS( 686 )
   W( 165 ) = JVS( 687 )
   W( 167 ) = JVS( 688 )
  a = -W( 110 ) / JVS( 544 )
  W( 110 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 545 )
  W( 167 ) = W( 167 ) + a*JVS( 546 )
  a = -W( 124 ) / JVS( 594 )
  W( 124 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 595 )
  W( 165 ) = W( 165 ) + a*JVS( 596 )
  a = -W( 128 ) / JVS( 608 )
  W( 128 ) = -a
  W( 138 ) = W( 138 ) + a*JVS( 609 )
  W( 157 ) = W( 157 ) + a*JVS( 610 )
  W( 161 ) = W( 161 ) + a*JVS( 611 )
  W( 167 ) = W( 167 ) + a*JVS( 612 )
  a = -W( 130 ) / JVS( 620 )
  W( 130 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 621 )
  W( 156 ) = W( 156 ) + a*JVS( 622 )
  W( 157 ) = W( 157 ) + a*JVS( 623 )
  W( 165 ) = W( 165 ) + a*JVS( 624 )
  a = -W( 131 ) / JVS( 626 )
  W( 131 ) = -a
  W( 138 ) = W( 138 ) + a*JVS( 627 )
  W( 157 ) = W( 157 ) + a*JVS( 628 )
  W( 161 ) = W( 161 ) + a*JVS( 629 )
  W( 165 ) = W( 165 ) + a*JVS( 630 )
  a = -W( 132 ) / JVS( 636 )
  W( 132 ) = -a
  W( 142 ) = W( 142 ) + a*JVS( 637 )
  W( 146 ) = W( 146 ) + a*JVS( 638 )
  W( 148 ) = W( 148 ) + a*JVS( 639 )
  W( 156 ) = W( 156 ) + a*JVS( 640 )
  W( 157 ) = W( 157 ) + a*JVS( 641 )
  W( 165 ) = W( 165 ) + a*JVS( 642 )
  JVS( 668) = W( 110 )
  JVS( 669) = W( 124 )
  JVS( 670) = W( 128 )
  JVS( 671) = W( 130 )
  JVS( 672) = W( 131 )
  JVS( 673) = W( 132 )
  JVS( 674) = W( 134 )
  JVS( 675) = W( 138 )
  JVS( 676) = W( 142 )
  JVS( 677) = W( 145 )
  JVS( 678) = W( 146 )
  JVS( 679) = W( 147 )
  JVS( 680) = W( 148 )
  JVS( 681) = W( 149 )
  JVS( 682) = W( 150 )
  JVS( 683) = W( 153 )
  JVS( 684) = W( 156 )
  JVS( 685) = W( 157 )
  JVS( 686) = W( 161 )
  JVS( 687) = W( 165 )
  JVS( 688) = W( 167 )
  IF ( ABS( JVS( 689 )) < TINY(a) ) THEN
         IER = 135
         RETURN
  END IF
   W( 135 ) = JVS( 689 )
   W( 152 ) = JVS( 690 )
   W( 156 ) = JVS( 691 )
   W( 157 ) = JVS( 692 )
   W( 165 ) = JVS( 693 )
  JVS( 689) = W( 135 )
  JVS( 690) = W( 152 )
  JVS( 691) = W( 156 )
  JVS( 692) = W( 157 )
  JVS( 693) = W( 165 )
  IF ( ABS( JVS( 699 )) < TINY(a) ) THEN
         IER = 136
         RETURN
  END IF
   W( 107 ) = JVS( 694 )
   W( 112 ) = JVS( 695 )
   W( 113 ) = JVS( 696 )
   W( 114 ) = JVS( 697 )
   W( 125 ) = JVS( 698 )
   W( 136 ) = JVS( 699 )
   W( 139 ) = JVS( 700 )
   W( 143 ) = JVS( 701 )
   W( 144 ) = JVS( 702 )
   W( 147 ) = JVS( 703 )
   W( 151 ) = JVS( 704 )
   W( 156 ) = JVS( 705 )
   W( 157 ) = JVS( 706 )
   W( 165 ) = JVS( 707 )
   W( 167 ) = JVS( 708 )
  a = -W( 107 ) / JVS( 534 )
  W( 107 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 535 )
  a = -W( 112 ) / JVS( 550 )
  W( 112 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 551 )
  a = -W( 113 ) / JVS( 552 )
  W( 113 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 553 )
  a = -W( 114 ) / JVS( 555 )
  W( 114 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 556 )
  W( 167 ) = W( 167 ) + a*JVS( 557 )
  a = -W( 125 ) / JVS( 597 )
  W( 125 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 598 )
  JVS( 694) = W( 107 )
  JVS( 695) = W( 112 )
  JVS( 696) = W( 113 )
  JVS( 697) = W( 114 )
  JVS( 698) = W( 125 )
  JVS( 699) = W( 136 )
  JVS( 700) = W( 139 )
  JVS( 701) = W( 143 )
  JVS( 702) = W( 144 )
  JVS( 703) = W( 147 )
  JVS( 704) = W( 151 )
  JVS( 705) = W( 156 )
  JVS( 706) = W( 157 )
  JVS( 707) = W( 165 )
  JVS( 708) = W( 167 )
  IF ( ABS( JVS( 709 )) < TINY(a) ) THEN
         IER = 137
         RETURN
  END IF
   W( 137 ) = JVS( 709 )
   W( 152 ) = JVS( 710 )
   W( 156 ) = JVS( 711 )
   W( 157 ) = JVS( 712 )
   W( 165 ) = JVS( 713 )
  JVS( 709) = W( 137 )
  JVS( 710) = W( 152 )
  JVS( 711) = W( 156 )
  JVS( 712) = W( 157 )
  JVS( 713) = W( 165 )
  IF ( ABS( JVS( 716 )) < TINY(a) ) THEN
         IER = 138
         RETURN
  END IF
   W( 124 ) = JVS( 714 )
   W( 131 ) = JVS( 715 )
   W( 138 ) = JVS( 716 )
   W( 157 ) = JVS( 717 )
   W( 159 ) = JVS( 718 )
   W( 160 ) = JVS( 719 )
   W( 161 ) = JVS( 720 )
   W( 162 ) = JVS( 721 )
   W( 163 ) = JVS( 722 )
   W( 164 ) = JVS( 723 )
   W( 165 ) = JVS( 724 )
   W( 167 ) = JVS( 725 )
  a = -W( 124 ) / JVS( 594 )
  W( 124 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 595 )
  W( 165 ) = W( 165 ) + a*JVS( 596 )
  a = -W( 131 ) / JVS( 626 )
  W( 131 ) = -a
  W( 138 ) = W( 138 ) + a*JVS( 627 )
  W( 157 ) = W( 157 ) + a*JVS( 628 )
  W( 161 ) = W( 161 ) + a*JVS( 629 )
  W( 165 ) = W( 165 ) + a*JVS( 630 )
  JVS( 714) = W( 124 )
  JVS( 715) = W( 131 )
  JVS( 716) = W( 138 )
  JVS( 717) = W( 157 )
  JVS( 718) = W( 159 )
  JVS( 719) = W( 160 )
  JVS( 720) = W( 161 )
  JVS( 721) = W( 162 )
  JVS( 722) = W( 163 )
  JVS( 723) = W( 164 )
  JVS( 724) = W( 165 )
  JVS( 725) = W( 167 )
  IF ( ABS( JVS( 726 )) < TINY(a) ) THEN
         IER = 139
         RETURN
  END IF
   W( 139 ) = JVS( 726 )
   W( 152 ) = JVS( 727 )
   W( 156 ) = JVS( 728 )
   W( 157 ) = JVS( 729 )
   W( 165 ) = JVS( 730 )
  JVS( 726) = W( 139 )
  JVS( 727) = W( 152 )
  JVS( 728) = W( 156 )
  JVS( 729) = W( 157 )
  JVS( 730) = W( 165 )
  IF ( ABS( JVS( 731 )) < TINY(a) ) THEN
         IER = 140
         RETURN
  END IF
   W( 140 ) = JVS( 731 )
   W( 152 ) = JVS( 732 )
   W( 156 ) = JVS( 733 )
   W( 157 ) = JVS( 734 )
   W( 165 ) = JVS( 735 )
  JVS( 731) = W( 140 )
  JVS( 732) = W( 152 )
  JVS( 733) = W( 156 )
  JVS( 734) = W( 157 )
  JVS( 735) = W( 165 )
  IF ( ABS( JVS( 744 )) < TINY(a) ) THEN
         IER = 141
         RETURN
  END IF
   W( 112 ) = JVS( 736 )
   W( 113 ) = JVS( 737 )
   W( 122 ) = JVS( 738 )
   W( 123 ) = JVS( 739 )
   W( 125 ) = JVS( 740 )
   W( 136 ) = JVS( 741 )
   W( 139 ) = JVS( 742 )
   W( 140 ) = JVS( 743 )
   W( 141 ) = JVS( 744 )
   W( 143 ) = JVS( 745 )
   W( 144 ) = JVS( 746 )
   W( 147 ) = JVS( 747 )
   W( 148 ) = JVS( 748 )
   W( 151 ) = JVS( 749 )
   W( 152 ) = JVS( 750 )
   W( 154 ) = JVS( 751 )
   W( 155 ) = JVS( 752 )
   W( 156 ) = JVS( 753 )
   W( 157 ) = JVS( 754 )
   W( 158 ) = JVS( 755 )
   W( 159 ) = JVS( 756 )
   W( 160 ) = JVS( 757 )
   W( 161 ) = JVS( 758 )
   W( 162 ) = JVS( 759 )
   W( 163 ) = JVS( 760 )
   W( 164 ) = JVS( 761 )
   W( 165 ) = JVS( 762 )
   W( 166 ) = JVS( 763 )
   W( 167 ) = JVS( 764 )
   W( 168 ) = JVS( 765 )
  a = -W( 112 ) / JVS( 550 )
  W( 112 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 551 )
  a = -W( 113 ) / JVS( 552 )
  W( 113 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 553 )
  a = -W( 122 ) / JVS( 586 )
  W( 122 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 587 )
  a = -W( 123 ) / JVS( 590 )
  W( 123 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 591 )
  a = -W( 125 ) / JVS( 597 )
  W( 125 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 598 )
  a = -W( 136 ) / JVS( 699 )
  W( 136 ) = -a
  W( 139 ) = W( 139 ) + a*JVS( 700 )
  W( 143 ) = W( 143 ) + a*JVS( 701 )
  W( 144 ) = W( 144 ) + a*JVS( 702 )
  W( 147 ) = W( 147 ) + a*JVS( 703 )
  W( 151 ) = W( 151 ) + a*JVS( 704 )
  W( 156 ) = W( 156 ) + a*JVS( 705 )
  W( 157 ) = W( 157 ) + a*JVS( 706 )
  W( 165 ) = W( 165 ) + a*JVS( 707 )
  W( 167 ) = W( 167 ) + a*JVS( 708 )
  a = -W( 139 ) / JVS( 726 )
  W( 139 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 727 )
  W( 156 ) = W( 156 ) + a*JVS( 728 )
  W( 157 ) = W( 157 ) + a*JVS( 729 )
  W( 165 ) = W( 165 ) + a*JVS( 730 )
  a = -W( 140 ) / JVS( 731 )
  W( 140 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 732 )
  W( 156 ) = W( 156 ) + a*JVS( 733 )
  W( 157 ) = W( 157 ) + a*JVS( 734 )
  W( 165 ) = W( 165 ) + a*JVS( 735 )
  JVS( 736) = W( 112 )
  JVS( 737) = W( 113 )
  JVS( 738) = W( 122 )
  JVS( 739) = W( 123 )
  JVS( 740) = W( 125 )
  JVS( 741) = W( 136 )
  JVS( 742) = W( 139 )
  JVS( 743) = W( 140 )
  JVS( 744) = W( 141 )
  JVS( 745) = W( 143 )
  JVS( 746) = W( 144 )
  JVS( 747) = W( 147 )
  JVS( 748) = W( 148 )
  JVS( 749) = W( 151 )
  JVS( 750) = W( 152 )
  JVS( 751) = W( 154 )
  JVS( 752) = W( 155 )
  JVS( 753) = W( 156 )
  JVS( 754) = W( 157 )
  JVS( 755) = W( 158 )
  JVS( 756) = W( 159 )
  JVS( 757) = W( 160 )
  JVS( 758) = W( 161 )
  JVS( 759) = W( 162 )
  JVS( 760) = W( 163 )
  JVS( 761) = W( 164 )
  JVS( 762) = W( 165 )
  JVS( 763) = W( 166 )
  JVS( 764) = W( 167 )
  JVS( 765) = W( 168 )
  IF ( ABS( JVS( 767 )) < TINY(a) ) THEN
         IER = 142
         RETURN
  END IF
   W( 140 ) = JVS( 766 )
   W( 142 ) = JVS( 767 )
   W( 147 ) = JVS( 768 )
   W( 152 ) = JVS( 769 )
   W( 156 ) = JVS( 770 )
   W( 157 ) = JVS( 771 )
   W( 165 ) = JVS( 772 )
  a = -W( 140 ) / JVS( 731 )
  W( 140 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 732 )
  W( 156 ) = W( 156 ) + a*JVS( 733 )
  W( 157 ) = W( 157 ) + a*JVS( 734 )
  W( 165 ) = W( 165 ) + a*JVS( 735 )
  JVS( 766) = W( 140 )
  JVS( 767) = W( 142 )
  JVS( 768) = W( 147 )
  JVS( 769) = W( 152 )
  JVS( 770) = W( 156 )
  JVS( 771) = W( 157 )
  JVS( 772) = W( 165 )
  IF ( ABS( JVS( 773 )) < TINY(a) ) THEN
         IER = 143
         RETURN
  END IF
   W( 143 ) = JVS( 773 )
   W( 152 ) = JVS( 774 )
   W( 156 ) = JVS( 775 )
   W( 157 ) = JVS( 776 )
   W( 165 ) = JVS( 777 )
  JVS( 773) = W( 143 )
  JVS( 774) = W( 152 )
  JVS( 775) = W( 156 )
  JVS( 776) = W( 157 )
  JVS( 777) = W( 165 )
  IF ( ABS( JVS( 778 )) < TINY(a) ) THEN
         IER = 144
         RETURN
  END IF
   W( 144 ) = JVS( 778 )
   W( 152 ) = JVS( 779 )
   W( 156 ) = JVS( 780 )
   W( 157 ) = JVS( 781 )
   W( 165 ) = JVS( 782 )
  JVS( 778) = W( 144 )
  JVS( 779) = W( 152 )
  JVS( 780) = W( 156 )
  JVS( 781) = W( 157 )
  JVS( 782) = W( 165 )
  IF ( ABS( JVS( 794 )) < TINY(a) ) THEN
         IER = 145
         RETURN
  END IF
   W( 115 ) = JVS( 783 )
   W( 121 ) = JVS( 784 )
   W( 122 ) = JVS( 785 )
   W( 123 ) = JVS( 786 )
   W( 126 ) = JVS( 787 )
   W( 127 ) = JVS( 788 )
   W( 131 ) = JVS( 789 )
   W( 135 ) = JVS( 790 )
   W( 138 ) = JVS( 791 )
   W( 143 ) = JVS( 792 )
   W( 144 ) = JVS( 793 )
   W( 145 ) = JVS( 794 )
   W( 146 ) = JVS( 795 )
   W( 152 ) = JVS( 796 )
   W( 156 ) = JVS( 797 )
   W( 157 ) = JVS( 798 )
   W( 159 ) = JVS( 799 )
   W( 160 ) = JVS( 800 )
   W( 161 ) = JVS( 801 )
   W( 162 ) = JVS( 802 )
   W( 163 ) = JVS( 803 )
   W( 164 ) = JVS( 804 )
   W( 165 ) = JVS( 805 )
   W( 167 ) = JVS( 806 )
  a = -W( 115 ) / JVS( 558 )
  W( 115 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 559 )
  a = -W( 121 ) / JVS( 582 )
  W( 121 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 583 )
  a = -W( 122 ) / JVS( 586 )
  W( 122 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 587 )
  a = -W( 123 ) / JVS( 590 )
  W( 123 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 591 )
  a = -W( 126 ) / JVS( 599 )
  W( 126 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 600 )
  W( 165 ) = W( 165 ) + a*JVS( 601 )
  a = -W( 127 ) / JVS( 604 )
  W( 127 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 605 )
  W( 165 ) = W( 165 ) + a*JVS( 606 )
  a = -W( 131 ) / JVS( 626 )
  W( 131 ) = -a
  W( 138 ) = W( 138 ) + a*JVS( 627 )
  W( 157 ) = W( 157 ) + a*JVS( 628 )
  W( 161 ) = W( 161 ) + a*JVS( 629 )
  W( 165 ) = W( 165 ) + a*JVS( 630 )
  a = -W( 135 ) / JVS( 689 )
  W( 135 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 690 )
  W( 156 ) = W( 156 ) + a*JVS( 691 )
  W( 157 ) = W( 157 ) + a*JVS( 692 )
  W( 165 ) = W( 165 ) + a*JVS( 693 )
  a = -W( 138 ) / JVS( 716 )
  W( 138 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 717 )
  W( 159 ) = W( 159 ) + a*JVS( 718 )
  W( 160 ) = W( 160 ) + a*JVS( 719 )
  W( 161 ) = W( 161 ) + a*JVS( 720 )
  W( 162 ) = W( 162 ) + a*JVS( 721 )
  W( 163 ) = W( 163 ) + a*JVS( 722 )
  W( 164 ) = W( 164 ) + a*JVS( 723 )
  W( 165 ) = W( 165 ) + a*JVS( 724 )
  W( 167 ) = W( 167 ) + a*JVS( 725 )
  a = -W( 143 ) / JVS( 773 )
  W( 143 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 774 )
  W( 156 ) = W( 156 ) + a*JVS( 775 )
  W( 157 ) = W( 157 ) + a*JVS( 776 )
  W( 165 ) = W( 165 ) + a*JVS( 777 )
  a = -W( 144 ) / JVS( 778 )
  W( 144 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 779 )
  W( 156 ) = W( 156 ) + a*JVS( 780 )
  W( 157 ) = W( 157 ) + a*JVS( 781 )
  W( 165 ) = W( 165 ) + a*JVS( 782 )
  JVS( 783) = W( 115 )
  JVS( 784) = W( 121 )
  JVS( 785) = W( 122 )
  JVS( 786) = W( 123 )
  JVS( 787) = W( 126 )
  JVS( 788) = W( 127 )
  JVS( 789) = W( 131 )
  JVS( 790) = W( 135 )
  JVS( 791) = W( 138 )
  JVS( 792) = W( 143 )
  JVS( 793) = W( 144 )
  JVS( 794) = W( 145 )
  JVS( 795) = W( 146 )
  JVS( 796) = W( 152 )
  JVS( 797) = W( 156 )
  JVS( 798) = W( 157 )
  JVS( 799) = W( 159 )
  JVS( 800) = W( 160 )
  JVS( 801) = W( 161 )
  JVS( 802) = W( 162 )
  JVS( 803) = W( 163 )
  JVS( 804) = W( 164 )
  JVS( 805) = W( 165 )
  JVS( 806) = W( 167 )
  IF ( ABS( JVS( 808 )) < TINY(a) ) THEN
         IER = 146
         RETURN
  END IF
   W( 140 ) = JVS( 807 )
   W( 146 ) = JVS( 808 )
   W( 147 ) = JVS( 809 )
   W( 152 ) = JVS( 810 )
   W( 156 ) = JVS( 811 )
   W( 157 ) = JVS( 812 )
   W( 165 ) = JVS( 813 )
  a = -W( 140 ) / JVS( 731 )
  W( 140 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 732 )
  W( 156 ) = W( 156 ) + a*JVS( 733 )
  W( 157 ) = W( 157 ) + a*JVS( 734 )
  W( 165 ) = W( 165 ) + a*JVS( 735 )
  JVS( 807) = W( 140 )
  JVS( 808) = W( 146 )
  JVS( 809) = W( 147 )
  JVS( 810) = W( 152 )
  JVS( 811) = W( 156 )
  JVS( 812) = W( 157 )
  JVS( 813) = W( 165 )
  IF ( ABS( JVS( 814 )) < TINY(a) ) THEN
         IER = 147
         RETURN
  END IF
   W( 147 ) = JVS( 814 )
   W( 152 ) = JVS( 815 )
   W( 156 ) = JVS( 816 )
   W( 157 ) = JVS( 817 )
   W( 165 ) = JVS( 818 )
  JVS( 814) = W( 147 )
  JVS( 815) = W( 152 )
  JVS( 816) = W( 156 )
  JVS( 817) = W( 157 )
  JVS( 818) = W( 165 )
  IF ( ABS( JVS( 821 )) < TINY(a) ) THEN
         IER = 148
         RETURN
  END IF
   W( 140 ) = JVS( 819 )
   W( 147 ) = JVS( 820 )
   W( 148 ) = JVS( 821 )
   W( 152 ) = JVS( 822 )
   W( 156 ) = JVS( 823 )
   W( 157 ) = JVS( 824 )
   W( 165 ) = JVS( 825 )
  a = -W( 140 ) / JVS( 731 )
  W( 140 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 732 )
  W( 156 ) = W( 156 ) + a*JVS( 733 )
  W( 157 ) = W( 157 ) + a*JVS( 734 )
  W( 165 ) = W( 165 ) + a*JVS( 735 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  JVS( 819) = W( 140 )
  JVS( 820) = W( 147 )
  JVS( 821) = W( 148 )
  JVS( 822) = W( 152 )
  JVS( 823) = W( 156 )
  JVS( 824) = W( 157 )
  JVS( 825) = W( 165 )
  IF ( ABS( JVS( 836 )) < TINY(a) ) THEN
         IER = 149
         RETURN
  END IF
   W( 101 ) = JVS( 826 )
   W( 108 ) = JVS( 827 )
   W( 112 ) = JVS( 828 )
   W( 113 ) = JVS( 829 )
   W( 125 ) = JVS( 830 )
   W( 135 ) = JVS( 831 )
   W( 137 ) = JVS( 832 )
   W( 139 ) = JVS( 833 )
   W( 146 ) = JVS( 834 )
   W( 147 ) = JVS( 835 )
   W( 149 ) = JVS( 836 )
   W( 151 ) = JVS( 837 )
   W( 152 ) = JVS( 838 )
   W( 153 ) = JVS( 839 )
   W( 154 ) = JVS( 840 )
   W( 155 ) = JVS( 841 )
   W( 156 ) = JVS( 842 )
   W( 157 ) = JVS( 843 )
   W( 159 ) = JVS( 844 )
   W( 160 ) = JVS( 845 )
   W( 162 ) = JVS( 846 )
   W( 163 ) = JVS( 847 )
   W( 164 ) = JVS( 848 )
   W( 165 ) = JVS( 849 )
  a = -W( 101 ) / JVS( 517 )
  W( 101 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 518 )
  a = -W( 108 ) / JVS( 536 )
  W( 108 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 537 )
  a = -W( 112 ) / JVS( 550 )
  W( 112 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 551 )
  a = -W( 113 ) / JVS( 552 )
  W( 113 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 553 )
  a = -W( 125 ) / JVS( 597 )
  W( 125 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 598 )
  a = -W( 135 ) / JVS( 689 )
  W( 135 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 690 )
  W( 156 ) = W( 156 ) + a*JVS( 691 )
  W( 157 ) = W( 157 ) + a*JVS( 692 )
  W( 165 ) = W( 165 ) + a*JVS( 693 )
  a = -W( 137 ) / JVS( 709 )
  W( 137 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 710 )
  W( 156 ) = W( 156 ) + a*JVS( 711 )
  W( 157 ) = W( 157 ) + a*JVS( 712 )
  W( 165 ) = W( 165 ) + a*JVS( 713 )
  a = -W( 139 ) / JVS( 726 )
  W( 139 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 727 )
  W( 156 ) = W( 156 ) + a*JVS( 728 )
  W( 157 ) = W( 157 ) + a*JVS( 729 )
  W( 165 ) = W( 165 ) + a*JVS( 730 )
  a = -W( 146 ) / JVS( 808 )
  W( 146 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 809 )
  W( 152 ) = W( 152 ) + a*JVS( 810 )
  W( 156 ) = W( 156 ) + a*JVS( 811 )
  W( 157 ) = W( 157 ) + a*JVS( 812 )
  W( 165 ) = W( 165 ) + a*JVS( 813 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  JVS( 826) = W( 101 )
  JVS( 827) = W( 108 )
  JVS( 828) = W( 112 )
  JVS( 829) = W( 113 )
  JVS( 830) = W( 125 )
  JVS( 831) = W( 135 )
  JVS( 832) = W( 137 )
  JVS( 833) = W( 139 )
  JVS( 834) = W( 146 )
  JVS( 835) = W( 147 )
  JVS( 836) = W( 149 )
  JVS( 837) = W( 151 )
  JVS( 838) = W( 152 )
  JVS( 839) = W( 153 )
  JVS( 840) = W( 154 )
  JVS( 841) = W( 155 )
  JVS( 842) = W( 156 )
  JVS( 843) = W( 157 )
  JVS( 844) = W( 159 )
  JVS( 845) = W( 160 )
  JVS( 846) = W( 162 )
  JVS( 847) = W( 163 )
  JVS( 848) = W( 164 )
  JVS( 849) = W( 165 )
  IF ( ABS( JVS( 870 )) < TINY(a) ) THEN
         IER = 150
         RETURN
  END IF
   W( 108 ) = JVS( 850 )
   W( 112 ) = JVS( 851 )
   W( 113 ) = JVS( 852 )
   W( 116 ) = JVS( 853 )
   W( 117 ) = JVS( 854 )
   W( 120 ) = JVS( 855 )
   W( 125 ) = JVS( 856 )
   W( 126 ) = JVS( 857 )
   W( 135 ) = JVS( 858 )
   W( 136 ) = JVS( 859 )
   W( 137 ) = JVS( 860 )
   W( 139 ) = JVS( 861 )
   W( 140 ) = JVS( 862 )
   W( 142 ) = JVS( 863 )
   W( 143 ) = JVS( 864 )
   W( 144 ) = JVS( 865 )
   W( 145 ) = JVS( 866 )
   W( 146 ) = JVS( 867 )
   W( 147 ) = JVS( 868 )
   W( 148 ) = JVS( 869 )
   W( 150 ) = JVS( 870 )
   W( 151 ) = JVS( 871 )
   W( 152 ) = JVS( 872 )
   W( 154 ) = JVS( 873 )
   W( 155 ) = JVS( 874 )
   W( 156 ) = JVS( 875 )
   W( 157 ) = JVS( 876 )
   W( 158 ) = JVS( 877 )
   W( 159 ) = JVS( 878 )
   W( 160 ) = JVS( 879 )
   W( 161 ) = JVS( 880 )
   W( 162 ) = JVS( 881 )
   W( 163 ) = JVS( 882 )
   W( 164 ) = JVS( 883 )
   W( 165 ) = JVS( 884 )
   W( 166 ) = JVS( 885 )
   W( 167 ) = JVS( 886 )
   W( 168 ) = JVS( 887 )
  a = -W( 108 ) / JVS( 536 )
  W( 108 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 537 )
  a = -W( 112 ) / JVS( 550 )
  W( 112 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 551 )
  a = -W( 113 ) / JVS( 552 )
  W( 113 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 553 )
  a = -W( 116 ) / JVS( 560 )
  W( 116 ) = -a
  W( 161 ) = W( 161 ) + a*JVS( 561 )
  W( 165 ) = W( 165 ) + a*JVS( 562 )
  W( 168 ) = W( 168 ) + a*JVS( 563 )
  a = -W( 117 ) / JVS( 564 )
  W( 117 ) = -a
  W( 150 ) = W( 150 ) + a*JVS( 565 )
  W( 160 ) = W( 160 ) + a*JVS( 566 )
  W( 161 ) = W( 161 ) + a*JVS( 567 )
  a = -W( 120 ) / JVS( 577 )
  W( 120 ) = -a
  W( 158 ) = W( 158 ) + a*JVS( 578 )
  W( 165 ) = W( 165 ) + a*JVS( 579 )
  W( 166 ) = W( 166 ) + a*JVS( 580 )
  W( 168 ) = W( 168 ) + a*JVS( 581 )
  a = -W( 125 ) / JVS( 597 )
  W( 125 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 598 )
  a = -W( 126 ) / JVS( 599 )
  W( 126 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 600 )
  W( 165 ) = W( 165 ) + a*JVS( 601 )
  a = -W( 135 ) / JVS( 689 )
  W( 135 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 690 )
  W( 156 ) = W( 156 ) + a*JVS( 691 )
  W( 157 ) = W( 157 ) + a*JVS( 692 )
  W( 165 ) = W( 165 ) + a*JVS( 693 )
  a = -W( 136 ) / JVS( 699 )
  W( 136 ) = -a
  W( 139 ) = W( 139 ) + a*JVS( 700 )
  W( 143 ) = W( 143 ) + a*JVS( 701 )
  W( 144 ) = W( 144 ) + a*JVS( 702 )
  W( 147 ) = W( 147 ) + a*JVS( 703 )
  W( 151 ) = W( 151 ) + a*JVS( 704 )
  W( 156 ) = W( 156 ) + a*JVS( 705 )
  W( 157 ) = W( 157 ) + a*JVS( 706 )
  W( 165 ) = W( 165 ) + a*JVS( 707 )
  W( 167 ) = W( 167 ) + a*JVS( 708 )
  a = -W( 137 ) / JVS( 709 )
  W( 137 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 710 )
  W( 156 ) = W( 156 ) + a*JVS( 711 )
  W( 157 ) = W( 157 ) + a*JVS( 712 )
  W( 165 ) = W( 165 ) + a*JVS( 713 )
  a = -W( 139 ) / JVS( 726 )
  W( 139 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 727 )
  W( 156 ) = W( 156 ) + a*JVS( 728 )
  W( 157 ) = W( 157 ) + a*JVS( 729 )
  W( 165 ) = W( 165 ) + a*JVS( 730 )
  a = -W( 140 ) / JVS( 731 )
  W( 140 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 732 )
  W( 156 ) = W( 156 ) + a*JVS( 733 )
  W( 157 ) = W( 157 ) + a*JVS( 734 )
  W( 165 ) = W( 165 ) + a*JVS( 735 )
  a = -W( 142 ) / JVS( 767 )
  W( 142 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 768 )
  W( 152 ) = W( 152 ) + a*JVS( 769 )
  W( 156 ) = W( 156 ) + a*JVS( 770 )
  W( 157 ) = W( 157 ) + a*JVS( 771 )
  W( 165 ) = W( 165 ) + a*JVS( 772 )
  a = -W( 143 ) / JVS( 773 )
  W( 143 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 774 )
  W( 156 ) = W( 156 ) + a*JVS( 775 )
  W( 157 ) = W( 157 ) + a*JVS( 776 )
  W( 165 ) = W( 165 ) + a*JVS( 777 )
  a = -W( 144 ) / JVS( 778 )
  W( 144 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 779 )
  W( 156 ) = W( 156 ) + a*JVS( 780 )
  W( 157 ) = W( 157 ) + a*JVS( 781 )
  W( 165 ) = W( 165 ) + a*JVS( 782 )
  a = -W( 145 ) / JVS( 794 )
  W( 145 ) = -a
  W( 146 ) = W( 146 ) + a*JVS( 795 )
  W( 152 ) = W( 152 ) + a*JVS( 796 )
  W( 156 ) = W( 156 ) + a*JVS( 797 )
  W( 157 ) = W( 157 ) + a*JVS( 798 )
  W( 159 ) = W( 159 ) + a*JVS( 799 )
  W( 160 ) = W( 160 ) + a*JVS( 800 )
  W( 161 ) = W( 161 ) + a*JVS( 801 )
  W( 162 ) = W( 162 ) + a*JVS( 802 )
  W( 163 ) = W( 163 ) + a*JVS( 803 )
  W( 164 ) = W( 164 ) + a*JVS( 804 )
  W( 165 ) = W( 165 ) + a*JVS( 805 )
  W( 167 ) = W( 167 ) + a*JVS( 806 )
  a = -W( 146 ) / JVS( 808 )
  W( 146 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 809 )
  W( 152 ) = W( 152 ) + a*JVS( 810 )
  W( 156 ) = W( 156 ) + a*JVS( 811 )
  W( 157 ) = W( 157 ) + a*JVS( 812 )
  W( 165 ) = W( 165 ) + a*JVS( 813 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  a = -W( 148 ) / JVS( 821 )
  W( 148 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 822 )
  W( 156 ) = W( 156 ) + a*JVS( 823 )
  W( 157 ) = W( 157 ) + a*JVS( 824 )
  W( 165 ) = W( 165 ) + a*JVS( 825 )
  JVS( 850) = W( 108 )
  JVS( 851) = W( 112 )
  JVS( 852) = W( 113 )
  JVS( 853) = W( 116 )
  JVS( 854) = W( 117 )
  JVS( 855) = W( 120 )
  JVS( 856) = W( 125 )
  JVS( 857) = W( 126 )
  JVS( 858) = W( 135 )
  JVS( 859) = W( 136 )
  JVS( 860) = W( 137 )
  JVS( 861) = W( 139 )
  JVS( 862) = W( 140 )
  JVS( 863) = W( 142 )
  JVS( 864) = W( 143 )
  JVS( 865) = W( 144 )
  JVS( 866) = W( 145 )
  JVS( 867) = W( 146 )
  JVS( 868) = W( 147 )
  JVS( 869) = W( 148 )
  JVS( 870) = W( 150 )
  JVS( 871) = W( 151 )
  JVS( 872) = W( 152 )
  JVS( 873) = W( 154 )
  JVS( 874) = W( 155 )
  JVS( 875) = W( 156 )
  JVS( 876) = W( 157 )
  JVS( 877) = W( 158 )
  JVS( 878) = W( 159 )
  JVS( 879) = W( 160 )
  JVS( 880) = W( 161 )
  JVS( 881) = W( 162 )
  JVS( 882) = W( 163 )
  JVS( 883) = W( 164 )
  JVS( 884) = W( 165 )
  JVS( 885) = W( 166 )
  JVS( 886) = W( 167 )
  JVS( 887) = W( 168 )
  IF ( ABS( JVS( 894 )) < TINY(a) ) THEN
         IER = 151
         RETURN
  END IF
   W( 114 ) = JVS( 888 )
   W( 139 ) = JVS( 889 )
   W( 143 ) = JVS( 890 )
   W( 144 ) = JVS( 891 )
   W( 146 ) = JVS( 892 )
   W( 147 ) = JVS( 893 )
   W( 151 ) = JVS( 894 )
   W( 152 ) = JVS( 895 )
   W( 156 ) = JVS( 896 )
   W( 157 ) = JVS( 897 )
   W( 158 ) = JVS( 898 )
   W( 160 ) = JVS( 899 )
   W( 165 ) = JVS( 900 )
   W( 167 ) = JVS( 901 )
  a = -W( 114 ) / JVS( 555 )
  W( 114 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 556 )
  W( 167 ) = W( 167 ) + a*JVS( 557 )
  a = -W( 139 ) / JVS( 726 )
  W( 139 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 727 )
  W( 156 ) = W( 156 ) + a*JVS( 728 )
  W( 157 ) = W( 157 ) + a*JVS( 729 )
  W( 165 ) = W( 165 ) + a*JVS( 730 )
  a = -W( 143 ) / JVS( 773 )
  W( 143 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 774 )
  W( 156 ) = W( 156 ) + a*JVS( 775 )
  W( 157 ) = W( 157 ) + a*JVS( 776 )
  W( 165 ) = W( 165 ) + a*JVS( 777 )
  a = -W( 144 ) / JVS( 778 )
  W( 144 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 779 )
  W( 156 ) = W( 156 ) + a*JVS( 780 )
  W( 157 ) = W( 157 ) + a*JVS( 781 )
  W( 165 ) = W( 165 ) + a*JVS( 782 )
  a = -W( 146 ) / JVS( 808 )
  W( 146 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 809 )
  W( 152 ) = W( 152 ) + a*JVS( 810 )
  W( 156 ) = W( 156 ) + a*JVS( 811 )
  W( 157 ) = W( 157 ) + a*JVS( 812 )
  W( 165 ) = W( 165 ) + a*JVS( 813 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  JVS( 888) = W( 114 )
  JVS( 889) = W( 139 )
  JVS( 890) = W( 143 )
  JVS( 891) = W( 144 )
  JVS( 892) = W( 146 )
  JVS( 893) = W( 147 )
  JVS( 894) = W( 151 )
  JVS( 895) = W( 152 )
  JVS( 896) = W( 156 )
  JVS( 897) = W( 157 )
  JVS( 898) = W( 158 )
  JVS( 899) = W( 160 )
  JVS( 900) = W( 165 )
  JVS( 901) = W( 167 )
  IF ( ABS( JVS( 912 )) < TINY(a) ) THEN
         IER = 152
         RETURN
  END IF
   W( 90 ) = JVS( 902 )
   W( 135 ) = JVS( 903 )
   W( 137 ) = JVS( 904 )
   W( 139 ) = JVS( 905 )
   W( 140 ) = JVS( 906 )
   W( 142 ) = JVS( 907 )
   W( 143 ) = JVS( 908 )
   W( 144 ) = JVS( 909 )
   W( 147 ) = JVS( 910 )
   W( 148 ) = JVS( 911 )
   W( 152 ) = JVS( 912 )
   W( 156 ) = JVS( 913 )
   W( 157 ) = JVS( 914 )
   W( 160 ) = JVS( 915 )
   W( 165 ) = JVS( 916 )
   W( 167 ) = JVS( 917 )
  a = -W( 90 ) / JVS( 480 )
  W( 90 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 481 )
  a = -W( 135 ) / JVS( 689 )
  W( 135 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 690 )
  W( 156 ) = W( 156 ) + a*JVS( 691 )
  W( 157 ) = W( 157 ) + a*JVS( 692 )
  W( 165 ) = W( 165 ) + a*JVS( 693 )
  a = -W( 137 ) / JVS( 709 )
  W( 137 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 710 )
  W( 156 ) = W( 156 ) + a*JVS( 711 )
  W( 157 ) = W( 157 ) + a*JVS( 712 )
  W( 165 ) = W( 165 ) + a*JVS( 713 )
  a = -W( 139 ) / JVS( 726 )
  W( 139 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 727 )
  W( 156 ) = W( 156 ) + a*JVS( 728 )
  W( 157 ) = W( 157 ) + a*JVS( 729 )
  W( 165 ) = W( 165 ) + a*JVS( 730 )
  a = -W( 140 ) / JVS( 731 )
  W( 140 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 732 )
  W( 156 ) = W( 156 ) + a*JVS( 733 )
  W( 157 ) = W( 157 ) + a*JVS( 734 )
  W( 165 ) = W( 165 ) + a*JVS( 735 )
  a = -W( 142 ) / JVS( 767 )
  W( 142 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 768 )
  W( 152 ) = W( 152 ) + a*JVS( 769 )
  W( 156 ) = W( 156 ) + a*JVS( 770 )
  W( 157 ) = W( 157 ) + a*JVS( 771 )
  W( 165 ) = W( 165 ) + a*JVS( 772 )
  a = -W( 143 ) / JVS( 773 )
  W( 143 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 774 )
  W( 156 ) = W( 156 ) + a*JVS( 775 )
  W( 157 ) = W( 157 ) + a*JVS( 776 )
  W( 165 ) = W( 165 ) + a*JVS( 777 )
  a = -W( 144 ) / JVS( 778 )
  W( 144 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 779 )
  W( 156 ) = W( 156 ) + a*JVS( 780 )
  W( 157 ) = W( 157 ) + a*JVS( 781 )
  W( 165 ) = W( 165 ) + a*JVS( 782 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  a = -W( 148 ) / JVS( 821 )
  W( 148 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 822 )
  W( 156 ) = W( 156 ) + a*JVS( 823 )
  W( 157 ) = W( 157 ) + a*JVS( 824 )
  W( 165 ) = W( 165 ) + a*JVS( 825 )
  JVS( 902) = W( 90 )
  JVS( 903) = W( 135 )
  JVS( 904) = W( 137 )
  JVS( 905) = W( 139 )
  JVS( 906) = W( 140 )
  JVS( 907) = W( 142 )
  JVS( 908) = W( 143 )
  JVS( 909) = W( 144 )
  JVS( 910) = W( 147 )
  JVS( 911) = W( 148 )
  JVS( 912) = W( 152 )
  JVS( 913) = W( 156 )
  JVS( 914) = W( 157 )
  JVS( 915) = W( 160 )
  JVS( 916) = W( 165 )
  JVS( 917) = W( 167 )
  IF ( ABS( JVS( 937 )) < TINY(a) ) THEN
         IER = 153
         RETURN
  END IF
   W( 107 ) = JVS( 918 )
   W( 112 ) = JVS( 919 )
   W( 113 ) = JVS( 920 )
   W( 122 ) = JVS( 921 )
   W( 123 ) = JVS( 922 )
   W( 125 ) = JVS( 923 )
   W( 127 ) = JVS( 924 )
   W( 129 ) = JVS( 925 )
   W( 135 ) = JVS( 926 )
   W( 137 ) = JVS( 927 )
   W( 139 ) = JVS( 928 )
   W( 142 ) = JVS( 929 )
   W( 143 ) = JVS( 930 )
   W( 144 ) = JVS( 931 )
   W( 146 ) = JVS( 932 )
   W( 147 ) = JVS( 933 )
   W( 148 ) = JVS( 934 )
   W( 151 ) = JVS( 935 )
   W( 152 ) = JVS( 936 )
   W( 153 ) = JVS( 937 )
   W( 154 ) = JVS( 938 )
   W( 155 ) = JVS( 939 )
   W( 156 ) = JVS( 940 )
   W( 157 ) = JVS( 941 )
   W( 158 ) = JVS( 942 )
   W( 160 ) = JVS( 943 )
   W( 161 ) = JVS( 944 )
   W( 165 ) = JVS( 945 )
   W( 166 ) = JVS( 946 )
   W( 167 ) = JVS( 947 )
  a = -W( 107 ) / JVS( 534 )
  W( 107 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 535 )
  a = -W( 112 ) / JVS( 550 )
  W( 112 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 551 )
  a = -W( 113 ) / JVS( 552 )
  W( 113 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 553 )
  a = -W( 122 ) / JVS( 586 )
  W( 122 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 587 )
  a = -W( 123 ) / JVS( 590 )
  W( 123 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 591 )
  a = -W( 125 ) / JVS( 597 )
  W( 125 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 598 )
  a = -W( 127 ) / JVS( 604 )
  W( 127 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 605 )
  W( 165 ) = W( 165 ) + a*JVS( 606 )
  a = -W( 129 ) / JVS( 613 )
  W( 129 ) = -a
  W( 158 ) = W( 158 ) + a*JVS( 614 )
  W( 161 ) = W( 161 ) + a*JVS( 615 )
  W( 165 ) = W( 165 ) + a*JVS( 616 )
  W( 166 ) = W( 166 ) + a*JVS( 617 )
  a = -W( 135 ) / JVS( 689 )
  W( 135 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 690 )
  W( 156 ) = W( 156 ) + a*JVS( 691 )
  W( 157 ) = W( 157 ) + a*JVS( 692 )
  W( 165 ) = W( 165 ) + a*JVS( 693 )
  a = -W( 137 ) / JVS( 709 )
  W( 137 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 710 )
  W( 156 ) = W( 156 ) + a*JVS( 711 )
  W( 157 ) = W( 157 ) + a*JVS( 712 )
  W( 165 ) = W( 165 ) + a*JVS( 713 )
  a = -W( 139 ) / JVS( 726 )
  W( 139 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 727 )
  W( 156 ) = W( 156 ) + a*JVS( 728 )
  W( 157 ) = W( 157 ) + a*JVS( 729 )
  W( 165 ) = W( 165 ) + a*JVS( 730 )
  a = -W( 142 ) / JVS( 767 )
  W( 142 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 768 )
  W( 152 ) = W( 152 ) + a*JVS( 769 )
  W( 156 ) = W( 156 ) + a*JVS( 770 )
  W( 157 ) = W( 157 ) + a*JVS( 771 )
  W( 165 ) = W( 165 ) + a*JVS( 772 )
  a = -W( 143 ) / JVS( 773 )
  W( 143 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 774 )
  W( 156 ) = W( 156 ) + a*JVS( 775 )
  W( 157 ) = W( 157 ) + a*JVS( 776 )
  W( 165 ) = W( 165 ) + a*JVS( 777 )
  a = -W( 144 ) / JVS( 778 )
  W( 144 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 779 )
  W( 156 ) = W( 156 ) + a*JVS( 780 )
  W( 157 ) = W( 157 ) + a*JVS( 781 )
  W( 165 ) = W( 165 ) + a*JVS( 782 )
  a = -W( 146 ) / JVS( 808 )
  W( 146 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 809 )
  W( 152 ) = W( 152 ) + a*JVS( 810 )
  W( 156 ) = W( 156 ) + a*JVS( 811 )
  W( 157 ) = W( 157 ) + a*JVS( 812 )
  W( 165 ) = W( 165 ) + a*JVS( 813 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  a = -W( 148 ) / JVS( 821 )
  W( 148 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 822 )
  W( 156 ) = W( 156 ) + a*JVS( 823 )
  W( 157 ) = W( 157 ) + a*JVS( 824 )
  W( 165 ) = W( 165 ) + a*JVS( 825 )
  a = -W( 151 ) / JVS( 894 )
  W( 151 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 895 )
  W( 156 ) = W( 156 ) + a*JVS( 896 )
  W( 157 ) = W( 157 ) + a*JVS( 897 )
  W( 158 ) = W( 158 ) + a*JVS( 898 )
  W( 160 ) = W( 160 ) + a*JVS( 899 )
  W( 165 ) = W( 165 ) + a*JVS( 900 )
  W( 167 ) = W( 167 ) + a*JVS( 901 )
  a = -W( 152 ) / JVS( 912 )
  W( 152 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 913 )
  W( 157 ) = W( 157 ) + a*JVS( 914 )
  W( 160 ) = W( 160 ) + a*JVS( 915 )
  W( 165 ) = W( 165 ) + a*JVS( 916 )
  W( 167 ) = W( 167 ) + a*JVS( 917 )
  JVS( 918) = W( 107 )
  JVS( 919) = W( 112 )
  JVS( 920) = W( 113 )
  JVS( 921) = W( 122 )
  JVS( 922) = W( 123 )
  JVS( 923) = W( 125 )
  JVS( 924) = W( 127 )
  JVS( 925) = W( 129 )
  JVS( 926) = W( 135 )
  JVS( 927) = W( 137 )
  JVS( 928) = W( 139 )
  JVS( 929) = W( 142 )
  JVS( 930) = W( 143 )
  JVS( 931) = W( 144 )
  JVS( 932) = W( 146 )
  JVS( 933) = W( 147 )
  JVS( 934) = W( 148 )
  JVS( 935) = W( 151 )
  JVS( 936) = W( 152 )
  JVS( 937) = W( 153 )
  JVS( 938) = W( 154 )
  JVS( 939) = W( 155 )
  JVS( 940) = W( 156 )
  JVS( 941) = W( 157 )
  JVS( 942) = W( 158 )
  JVS( 943) = W( 160 )
  JVS( 944) = W( 161 )
  JVS( 945) = W( 165 )
  JVS( 946) = W( 166 )
  JVS( 947) = W( 167 )
  IF ( ABS( JVS( 959 )) < TINY(a) ) THEN
         IER = 154
         RETURN
  END IF
   W( 112 ) = JVS( 948 )
   W( 113 ) = JVS( 949 )
   W( 125 ) = JVS( 950 )
   W( 137 ) = JVS( 951 )
   W( 139 ) = JVS( 952 )
   W( 142 ) = JVS( 953 )
   W( 146 ) = JVS( 954 )
   W( 147 ) = JVS( 955 )
   W( 148 ) = JVS( 956 )
   W( 151 ) = JVS( 957 )
   W( 152 ) = JVS( 958 )
   W( 154 ) = JVS( 959 )
   W( 155 ) = JVS( 960 )
   W( 156 ) = JVS( 961 )
   W( 157 ) = JVS( 962 )
   W( 158 ) = JVS( 963 )
   W( 160 ) = JVS( 964 )
   W( 165 ) = JVS( 965 )
   W( 166 ) = JVS( 966 )
   W( 167 ) = JVS( 967 )
   W( 168 ) = JVS( 968 )
  a = -W( 112 ) / JVS( 550 )
  W( 112 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 551 )
  a = -W( 113 ) / JVS( 552 )
  W( 113 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 553 )
  a = -W( 125 ) / JVS( 597 )
  W( 125 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 598 )
  a = -W( 137 ) / JVS( 709 )
  W( 137 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 710 )
  W( 156 ) = W( 156 ) + a*JVS( 711 )
  W( 157 ) = W( 157 ) + a*JVS( 712 )
  W( 165 ) = W( 165 ) + a*JVS( 713 )
  a = -W( 139 ) / JVS( 726 )
  W( 139 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 727 )
  W( 156 ) = W( 156 ) + a*JVS( 728 )
  W( 157 ) = W( 157 ) + a*JVS( 729 )
  W( 165 ) = W( 165 ) + a*JVS( 730 )
  a = -W( 142 ) / JVS( 767 )
  W( 142 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 768 )
  W( 152 ) = W( 152 ) + a*JVS( 769 )
  W( 156 ) = W( 156 ) + a*JVS( 770 )
  W( 157 ) = W( 157 ) + a*JVS( 771 )
  W( 165 ) = W( 165 ) + a*JVS( 772 )
  a = -W( 146 ) / JVS( 808 )
  W( 146 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 809 )
  W( 152 ) = W( 152 ) + a*JVS( 810 )
  W( 156 ) = W( 156 ) + a*JVS( 811 )
  W( 157 ) = W( 157 ) + a*JVS( 812 )
  W( 165 ) = W( 165 ) + a*JVS( 813 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  a = -W( 148 ) / JVS( 821 )
  W( 148 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 822 )
  W( 156 ) = W( 156 ) + a*JVS( 823 )
  W( 157 ) = W( 157 ) + a*JVS( 824 )
  W( 165 ) = W( 165 ) + a*JVS( 825 )
  a = -W( 151 ) / JVS( 894 )
  W( 151 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 895 )
  W( 156 ) = W( 156 ) + a*JVS( 896 )
  W( 157 ) = W( 157 ) + a*JVS( 897 )
  W( 158 ) = W( 158 ) + a*JVS( 898 )
  W( 160 ) = W( 160 ) + a*JVS( 899 )
  W( 165 ) = W( 165 ) + a*JVS( 900 )
  W( 167 ) = W( 167 ) + a*JVS( 901 )
  a = -W( 152 ) / JVS( 912 )
  W( 152 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 913 )
  W( 157 ) = W( 157 ) + a*JVS( 914 )
  W( 160 ) = W( 160 ) + a*JVS( 915 )
  W( 165 ) = W( 165 ) + a*JVS( 916 )
  W( 167 ) = W( 167 ) + a*JVS( 917 )
  JVS( 948) = W( 112 )
  JVS( 949) = W( 113 )
  JVS( 950) = W( 125 )
  JVS( 951) = W( 137 )
  JVS( 952) = W( 139 )
  JVS( 953) = W( 142 )
  JVS( 954) = W( 146 )
  JVS( 955) = W( 147 )
  JVS( 956) = W( 148 )
  JVS( 957) = W( 151 )
  JVS( 958) = W( 152 )
  JVS( 959) = W( 154 )
  JVS( 960) = W( 155 )
  JVS( 961) = W( 156 )
  JVS( 962) = W( 157 )
  JVS( 963) = W( 158 )
  JVS( 964) = W( 160 )
  JVS( 965) = W( 165 )
  JVS( 966) = W( 166 )
  JVS( 967) = W( 167 )
  JVS( 968) = W( 168 )
  IF ( ABS( JVS( 981 )) < TINY(a) ) THEN
         IER = 155
         RETURN
  END IF
   W( 112 ) = JVS( 969 )
   W( 121 ) = JVS( 970 )
   W( 125 ) = JVS( 971 )
   W( 139 ) = JVS( 972 )
   W( 140 ) = JVS( 973 )
   W( 143 ) = JVS( 974 )
   W( 144 ) = JVS( 975 )
   W( 146 ) = JVS( 976 )
   W( 147 ) = JVS( 977 )
   W( 148 ) = JVS( 978 )
   W( 151 ) = JVS( 979 )
   W( 152 ) = JVS( 980 )
   W( 155 ) = JVS( 981 )
   W( 156 ) = JVS( 982 )
   W( 157 ) = JVS( 983 )
   W( 158 ) = JVS( 984 )
   W( 160 ) = JVS( 985 )
   W( 162 ) = JVS( 986 )
   W( 163 ) = JVS( 987 )
   W( 164 ) = JVS( 988 )
   W( 165 ) = JVS( 989 )
   W( 166 ) = JVS( 990 )
   W( 167 ) = JVS( 991 )
   W( 168 ) = JVS( 992 )
  a = -W( 112 ) / JVS( 550 )
  W( 112 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 551 )
  a = -W( 121 ) / JVS( 582 )
  W( 121 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 583 )
  a = -W( 125 ) / JVS( 597 )
  W( 125 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 598 )
  a = -W( 139 ) / JVS( 726 )
  W( 139 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 727 )
  W( 156 ) = W( 156 ) + a*JVS( 728 )
  W( 157 ) = W( 157 ) + a*JVS( 729 )
  W( 165 ) = W( 165 ) + a*JVS( 730 )
  a = -W( 140 ) / JVS( 731 )
  W( 140 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 732 )
  W( 156 ) = W( 156 ) + a*JVS( 733 )
  W( 157 ) = W( 157 ) + a*JVS( 734 )
  W( 165 ) = W( 165 ) + a*JVS( 735 )
  a = -W( 143 ) / JVS( 773 )
  W( 143 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 774 )
  W( 156 ) = W( 156 ) + a*JVS( 775 )
  W( 157 ) = W( 157 ) + a*JVS( 776 )
  W( 165 ) = W( 165 ) + a*JVS( 777 )
  a = -W( 144 ) / JVS( 778 )
  W( 144 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 779 )
  W( 156 ) = W( 156 ) + a*JVS( 780 )
  W( 157 ) = W( 157 ) + a*JVS( 781 )
  W( 165 ) = W( 165 ) + a*JVS( 782 )
  a = -W( 146 ) / JVS( 808 )
  W( 146 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 809 )
  W( 152 ) = W( 152 ) + a*JVS( 810 )
  W( 156 ) = W( 156 ) + a*JVS( 811 )
  W( 157 ) = W( 157 ) + a*JVS( 812 )
  W( 165 ) = W( 165 ) + a*JVS( 813 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  a = -W( 148 ) / JVS( 821 )
  W( 148 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 822 )
  W( 156 ) = W( 156 ) + a*JVS( 823 )
  W( 157 ) = W( 157 ) + a*JVS( 824 )
  W( 165 ) = W( 165 ) + a*JVS( 825 )
  a = -W( 151 ) / JVS( 894 )
  W( 151 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 895 )
  W( 156 ) = W( 156 ) + a*JVS( 896 )
  W( 157 ) = W( 157 ) + a*JVS( 897 )
  W( 158 ) = W( 158 ) + a*JVS( 898 )
  W( 160 ) = W( 160 ) + a*JVS( 899 )
  W( 165 ) = W( 165 ) + a*JVS( 900 )
  W( 167 ) = W( 167 ) + a*JVS( 901 )
  a = -W( 152 ) / JVS( 912 )
  W( 152 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 913 )
  W( 157 ) = W( 157 ) + a*JVS( 914 )
  W( 160 ) = W( 160 ) + a*JVS( 915 )
  W( 165 ) = W( 165 ) + a*JVS( 916 )
  W( 167 ) = W( 167 ) + a*JVS( 917 )
  JVS( 969) = W( 112 )
  JVS( 970) = W( 121 )
  JVS( 971) = W( 125 )
  JVS( 972) = W( 139 )
  JVS( 973) = W( 140 )
  JVS( 974) = W( 143 )
  JVS( 975) = W( 144 )
  JVS( 976) = W( 146 )
  JVS( 977) = W( 147 )
  JVS( 978) = W( 148 )
  JVS( 979) = W( 151 )
  JVS( 980) = W( 152 )
  JVS( 981) = W( 155 )
  JVS( 982) = W( 156 )
  JVS( 983) = W( 157 )
  JVS( 984) = W( 158 )
  JVS( 985) = W( 160 )
  JVS( 986) = W( 162 )
  JVS( 987) = W( 163 )
  JVS( 988) = W( 164 )
  JVS( 989) = W( 165 )
  JVS( 990) = W( 166 )
  JVS( 991) = W( 167 )
  JVS( 992) = W( 168 )
  IF ( ABS( JVS( 1006 )) < TINY(a) ) THEN
         IER = 156
         RETURN
  END IF
   W( 126 ) = JVS( 993 )
   W( 127 ) = JVS( 994 )
   W( 135 ) = JVS( 995 )
   W( 137 ) = JVS( 996 )
   W( 139 ) = JVS( 997 )
   W( 140 ) = JVS( 998 )
   W( 142 ) = JVS( 999 )
   W( 143 ) = JVS( 1000 )
   W( 144 ) = JVS( 1001 )
   W( 146 ) = JVS( 1002 )
   W( 147 ) = JVS( 1003 )
   W( 148 ) = JVS( 1004 )
   W( 152 ) = JVS( 1005 )
   W( 156 ) = JVS( 1006 )
   W( 157 ) = JVS( 1007 )
   W( 159 ) = JVS( 1008 )
   W( 160 ) = JVS( 1009 )
   W( 161 ) = JVS( 1010 )
   W( 162 ) = JVS( 1011 )
   W( 163 ) = JVS( 1012 )
   W( 164 ) = JVS( 1013 )
   W( 165 ) = JVS( 1014 )
   W( 167 ) = JVS( 1015 )
  a = -W( 126 ) / JVS( 599 )
  W( 126 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 600 )
  W( 165 ) = W( 165 ) + a*JVS( 601 )
  a = -W( 127 ) / JVS( 604 )
  W( 127 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 605 )
  W( 165 ) = W( 165 ) + a*JVS( 606 )
  a = -W( 135 ) / JVS( 689 )
  W( 135 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 690 )
  W( 156 ) = W( 156 ) + a*JVS( 691 )
  W( 157 ) = W( 157 ) + a*JVS( 692 )
  W( 165 ) = W( 165 ) + a*JVS( 693 )
  a = -W( 137 ) / JVS( 709 )
  W( 137 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 710 )
  W( 156 ) = W( 156 ) + a*JVS( 711 )
  W( 157 ) = W( 157 ) + a*JVS( 712 )
  W( 165 ) = W( 165 ) + a*JVS( 713 )
  a = -W( 139 ) / JVS( 726 )
  W( 139 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 727 )
  W( 156 ) = W( 156 ) + a*JVS( 728 )
  W( 157 ) = W( 157 ) + a*JVS( 729 )
  W( 165 ) = W( 165 ) + a*JVS( 730 )
  a = -W( 140 ) / JVS( 731 )
  W( 140 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 732 )
  W( 156 ) = W( 156 ) + a*JVS( 733 )
  W( 157 ) = W( 157 ) + a*JVS( 734 )
  W( 165 ) = W( 165 ) + a*JVS( 735 )
  a = -W( 142 ) / JVS( 767 )
  W( 142 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 768 )
  W( 152 ) = W( 152 ) + a*JVS( 769 )
  W( 156 ) = W( 156 ) + a*JVS( 770 )
  W( 157 ) = W( 157 ) + a*JVS( 771 )
  W( 165 ) = W( 165 ) + a*JVS( 772 )
  a = -W( 143 ) / JVS( 773 )
  W( 143 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 774 )
  W( 156 ) = W( 156 ) + a*JVS( 775 )
  W( 157 ) = W( 157 ) + a*JVS( 776 )
  W( 165 ) = W( 165 ) + a*JVS( 777 )
  a = -W( 144 ) / JVS( 778 )
  W( 144 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 779 )
  W( 156 ) = W( 156 ) + a*JVS( 780 )
  W( 157 ) = W( 157 ) + a*JVS( 781 )
  W( 165 ) = W( 165 ) + a*JVS( 782 )
  a = -W( 146 ) / JVS( 808 )
  W( 146 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 809 )
  W( 152 ) = W( 152 ) + a*JVS( 810 )
  W( 156 ) = W( 156 ) + a*JVS( 811 )
  W( 157 ) = W( 157 ) + a*JVS( 812 )
  W( 165 ) = W( 165 ) + a*JVS( 813 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  a = -W( 148 ) / JVS( 821 )
  W( 148 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 822 )
  W( 156 ) = W( 156 ) + a*JVS( 823 )
  W( 157 ) = W( 157 ) + a*JVS( 824 )
  W( 165 ) = W( 165 ) + a*JVS( 825 )
  a = -W( 152 ) / JVS( 912 )
  W( 152 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 913 )
  W( 157 ) = W( 157 ) + a*JVS( 914 )
  W( 160 ) = W( 160 ) + a*JVS( 915 )
  W( 165 ) = W( 165 ) + a*JVS( 916 )
  W( 167 ) = W( 167 ) + a*JVS( 917 )
  JVS( 993) = W( 126 )
  JVS( 994) = W( 127 )
  JVS( 995) = W( 135 )
  JVS( 996) = W( 137 )
  JVS( 997) = W( 139 )
  JVS( 998) = W( 140 )
  JVS( 999) = W( 142 )
  JVS( 1000) = W( 143 )
  JVS( 1001) = W( 144 )
  JVS( 1002) = W( 146 )
  JVS( 1003) = W( 147 )
  JVS( 1004) = W( 148 )
  JVS( 1005) = W( 152 )
  JVS( 1006) = W( 156 )
  JVS( 1007) = W( 157 )
  JVS( 1008) = W( 159 )
  JVS( 1009) = W( 160 )
  JVS( 1010) = W( 161 )
  JVS( 1011) = W( 162 )
  JVS( 1012) = W( 163 )
  JVS( 1013) = W( 164 )
  JVS( 1014) = W( 165 )
  JVS( 1015) = W( 167 )
  IF ( ABS( JVS( 1045 )) < TINY(a) ) THEN
         IER = 157
         RETURN
  END IF
   W( 110 ) = JVS( 1016 )
   W( 119 ) = JVS( 1017 )
   W( 124 ) = JVS( 1018 )
   W( 128 ) = JVS( 1019 )
   W( 130 ) = JVS( 1020 )
   W( 131 ) = JVS( 1021 )
   W( 132 ) = JVS( 1022 )
   W( 134 ) = JVS( 1023 )
   W( 135 ) = JVS( 1024 )
   W( 137 ) = JVS( 1025 )
   W( 138 ) = JVS( 1026 )
   W( 139 ) = JVS( 1027 )
   W( 140 ) = JVS( 1028 )
   W( 141 ) = JVS( 1029 )
   W( 142 ) = JVS( 1030 )
   W( 143 ) = JVS( 1031 )
   W( 144 ) = JVS( 1032 )
   W( 145 ) = JVS( 1033 )
   W( 146 ) = JVS( 1034 )
   W( 147 ) = JVS( 1035 )
   W( 148 ) = JVS( 1036 )
   W( 149 ) = JVS( 1037 )
   W( 150 ) = JVS( 1038 )
   W( 151 ) = JVS( 1039 )
   W( 152 ) = JVS( 1040 )
   W( 153 ) = JVS( 1041 )
   W( 154 ) = JVS( 1042 )
   W( 155 ) = JVS( 1043 )
   W( 156 ) = JVS( 1044 )
   W( 157 ) = JVS( 1045 )
   W( 158 ) = JVS( 1046 )
   W( 159 ) = JVS( 1047 )
   W( 160 ) = JVS( 1048 )
   W( 161 ) = JVS( 1049 )
   W( 162 ) = JVS( 1050 )
   W( 163 ) = JVS( 1051 )
   W( 164 ) = JVS( 1052 )
   W( 165 ) = JVS( 1053 )
   W( 166 ) = JVS( 1054 )
   W( 167 ) = JVS( 1055 )
   W( 168 ) = JVS( 1056 )
  a = -W( 110 ) / JVS( 544 )
  W( 110 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 545 )
  W( 167 ) = W( 167 ) + a*JVS( 546 )
  a = -W( 119 ) / JVS( 573 )
  W( 119 ) = -a
  W( 161 ) = W( 161 ) + a*JVS( 574 )
  W( 165 ) = W( 165 ) + a*JVS( 575 )
  W( 167 ) = W( 167 ) + a*JVS( 576 )
  a = -W( 124 ) / JVS( 594 )
  W( 124 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 595 )
  W( 165 ) = W( 165 ) + a*JVS( 596 )
  a = -W( 128 ) / JVS( 608 )
  W( 128 ) = -a
  W( 138 ) = W( 138 ) + a*JVS( 609 )
  W( 157 ) = W( 157 ) + a*JVS( 610 )
  W( 161 ) = W( 161 ) + a*JVS( 611 )
  W( 167 ) = W( 167 ) + a*JVS( 612 )
  a = -W( 130 ) / JVS( 620 )
  W( 130 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 621 )
  W( 156 ) = W( 156 ) + a*JVS( 622 )
  W( 157 ) = W( 157 ) + a*JVS( 623 )
  W( 165 ) = W( 165 ) + a*JVS( 624 )
  a = -W( 131 ) / JVS( 626 )
  W( 131 ) = -a
  W( 138 ) = W( 138 ) + a*JVS( 627 )
  W( 157 ) = W( 157 ) + a*JVS( 628 )
  W( 161 ) = W( 161 ) + a*JVS( 629 )
  W( 165 ) = W( 165 ) + a*JVS( 630 )
  a = -W( 132 ) / JVS( 636 )
  W( 132 ) = -a
  W( 142 ) = W( 142 ) + a*JVS( 637 )
  W( 146 ) = W( 146 ) + a*JVS( 638 )
  W( 148 ) = W( 148 ) + a*JVS( 639 )
  W( 156 ) = W( 156 ) + a*JVS( 640 )
  W( 157 ) = W( 157 ) + a*JVS( 641 )
  W( 165 ) = W( 165 ) + a*JVS( 642 )
  a = -W( 134 ) / JVS( 674 )
  W( 134 ) = -a
  W( 138 ) = W( 138 ) + a*JVS( 675 )
  W( 142 ) = W( 142 ) + a*JVS( 676 )
  W( 145 ) = W( 145 ) + a*JVS( 677 )
  W( 146 ) = W( 146 ) + a*JVS( 678 )
  W( 147 ) = W( 147 ) + a*JVS( 679 )
  W( 148 ) = W( 148 ) + a*JVS( 680 )
  W( 149 ) = W( 149 ) + a*JVS( 681 )
  W( 150 ) = W( 150 ) + a*JVS( 682 )
  W( 153 ) = W( 153 ) + a*JVS( 683 )
  W( 156 ) = W( 156 ) + a*JVS( 684 )
  W( 157 ) = W( 157 ) + a*JVS( 685 )
  W( 161 ) = W( 161 ) + a*JVS( 686 )
  W( 165 ) = W( 165 ) + a*JVS( 687 )
  W( 167 ) = W( 167 ) + a*JVS( 688 )
  a = -W( 135 ) / JVS( 689 )
  W( 135 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 690 )
  W( 156 ) = W( 156 ) + a*JVS( 691 )
  W( 157 ) = W( 157 ) + a*JVS( 692 )
  W( 165 ) = W( 165 ) + a*JVS( 693 )
  a = -W( 137 ) / JVS( 709 )
  W( 137 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 710 )
  W( 156 ) = W( 156 ) + a*JVS( 711 )
  W( 157 ) = W( 157 ) + a*JVS( 712 )
  W( 165 ) = W( 165 ) + a*JVS( 713 )
  a = -W( 138 ) / JVS( 716 )
  W( 138 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 717 )
  W( 159 ) = W( 159 ) + a*JVS( 718 )
  W( 160 ) = W( 160 ) + a*JVS( 719 )
  W( 161 ) = W( 161 ) + a*JVS( 720 )
  W( 162 ) = W( 162 ) + a*JVS( 721 )
  W( 163 ) = W( 163 ) + a*JVS( 722 )
  W( 164 ) = W( 164 ) + a*JVS( 723 )
  W( 165 ) = W( 165 ) + a*JVS( 724 )
  W( 167 ) = W( 167 ) + a*JVS( 725 )
  a = -W( 139 ) / JVS( 726 )
  W( 139 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 727 )
  W( 156 ) = W( 156 ) + a*JVS( 728 )
  W( 157 ) = W( 157 ) + a*JVS( 729 )
  W( 165 ) = W( 165 ) + a*JVS( 730 )
  a = -W( 140 ) / JVS( 731 )
  W( 140 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 732 )
  W( 156 ) = W( 156 ) + a*JVS( 733 )
  W( 157 ) = W( 157 ) + a*JVS( 734 )
  W( 165 ) = W( 165 ) + a*JVS( 735 )
  a = -W( 141 ) / JVS( 744 )
  W( 141 ) = -a
  W( 143 ) = W( 143 ) + a*JVS( 745 )
  W( 144 ) = W( 144 ) + a*JVS( 746 )
  W( 147 ) = W( 147 ) + a*JVS( 747 )
  W( 148 ) = W( 148 ) + a*JVS( 748 )
  W( 151 ) = W( 151 ) + a*JVS( 749 )
  W( 152 ) = W( 152 ) + a*JVS( 750 )
  W( 154 ) = W( 154 ) + a*JVS( 751 )
  W( 155 ) = W( 155 ) + a*JVS( 752 )
  W( 156 ) = W( 156 ) + a*JVS( 753 )
  W( 157 ) = W( 157 ) + a*JVS( 754 )
  W( 158 ) = W( 158 ) + a*JVS( 755 )
  W( 159 ) = W( 159 ) + a*JVS( 756 )
  W( 160 ) = W( 160 ) + a*JVS( 757 )
  W( 161 ) = W( 161 ) + a*JVS( 758 )
  W( 162 ) = W( 162 ) + a*JVS( 759 )
  W( 163 ) = W( 163 ) + a*JVS( 760 )
  W( 164 ) = W( 164 ) + a*JVS( 761 )
  W( 165 ) = W( 165 ) + a*JVS( 762 )
  W( 166 ) = W( 166 ) + a*JVS( 763 )
  W( 167 ) = W( 167 ) + a*JVS( 764 )
  W( 168 ) = W( 168 ) + a*JVS( 765 )
  a = -W( 142 ) / JVS( 767 )
  W( 142 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 768 )
  W( 152 ) = W( 152 ) + a*JVS( 769 )
  W( 156 ) = W( 156 ) + a*JVS( 770 )
  W( 157 ) = W( 157 ) + a*JVS( 771 )
  W( 165 ) = W( 165 ) + a*JVS( 772 )
  a = -W( 143 ) / JVS( 773 )
  W( 143 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 774 )
  W( 156 ) = W( 156 ) + a*JVS( 775 )
  W( 157 ) = W( 157 ) + a*JVS( 776 )
  W( 165 ) = W( 165 ) + a*JVS( 777 )
  a = -W( 144 ) / JVS( 778 )
  W( 144 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 779 )
  W( 156 ) = W( 156 ) + a*JVS( 780 )
  W( 157 ) = W( 157 ) + a*JVS( 781 )
  W( 165 ) = W( 165 ) + a*JVS( 782 )
  a = -W( 145 ) / JVS( 794 )
  W( 145 ) = -a
  W( 146 ) = W( 146 ) + a*JVS( 795 )
  W( 152 ) = W( 152 ) + a*JVS( 796 )
  W( 156 ) = W( 156 ) + a*JVS( 797 )
  W( 157 ) = W( 157 ) + a*JVS( 798 )
  W( 159 ) = W( 159 ) + a*JVS( 799 )
  W( 160 ) = W( 160 ) + a*JVS( 800 )
  W( 161 ) = W( 161 ) + a*JVS( 801 )
  W( 162 ) = W( 162 ) + a*JVS( 802 )
  W( 163 ) = W( 163 ) + a*JVS( 803 )
  W( 164 ) = W( 164 ) + a*JVS( 804 )
  W( 165 ) = W( 165 ) + a*JVS( 805 )
  W( 167 ) = W( 167 ) + a*JVS( 806 )
  a = -W( 146 ) / JVS( 808 )
  W( 146 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 809 )
  W( 152 ) = W( 152 ) + a*JVS( 810 )
  W( 156 ) = W( 156 ) + a*JVS( 811 )
  W( 157 ) = W( 157 ) + a*JVS( 812 )
  W( 165 ) = W( 165 ) + a*JVS( 813 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  a = -W( 148 ) / JVS( 821 )
  W( 148 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 822 )
  W( 156 ) = W( 156 ) + a*JVS( 823 )
  W( 157 ) = W( 157 ) + a*JVS( 824 )
  W( 165 ) = W( 165 ) + a*JVS( 825 )
  a = -W( 149 ) / JVS( 836 )
  W( 149 ) = -a
  W( 151 ) = W( 151 ) + a*JVS( 837 )
  W( 152 ) = W( 152 ) + a*JVS( 838 )
  W( 153 ) = W( 153 ) + a*JVS( 839 )
  W( 154 ) = W( 154 ) + a*JVS( 840 )
  W( 155 ) = W( 155 ) + a*JVS( 841 )
  W( 156 ) = W( 156 ) + a*JVS( 842 )
  W( 157 ) = W( 157 ) + a*JVS( 843 )
  W( 159 ) = W( 159 ) + a*JVS( 844 )
  W( 160 ) = W( 160 ) + a*JVS( 845 )
  W( 162 ) = W( 162 ) + a*JVS( 846 )
  W( 163 ) = W( 163 ) + a*JVS( 847 )
  W( 164 ) = W( 164 ) + a*JVS( 848 )
  W( 165 ) = W( 165 ) + a*JVS( 849 )
  a = -W( 150 ) / JVS( 870 )
  W( 150 ) = -a
  W( 151 ) = W( 151 ) + a*JVS( 871 )
  W( 152 ) = W( 152 ) + a*JVS( 872 )
  W( 154 ) = W( 154 ) + a*JVS( 873 )
  W( 155 ) = W( 155 ) + a*JVS( 874 )
  W( 156 ) = W( 156 ) + a*JVS( 875 )
  W( 157 ) = W( 157 ) + a*JVS( 876 )
  W( 158 ) = W( 158 ) + a*JVS( 877 )
  W( 159 ) = W( 159 ) + a*JVS( 878 )
  W( 160 ) = W( 160 ) + a*JVS( 879 )
  W( 161 ) = W( 161 ) + a*JVS( 880 )
  W( 162 ) = W( 162 ) + a*JVS( 881 )
  W( 163 ) = W( 163 ) + a*JVS( 882 )
  W( 164 ) = W( 164 ) + a*JVS( 883 )
  W( 165 ) = W( 165 ) + a*JVS( 884 )
  W( 166 ) = W( 166 ) + a*JVS( 885 )
  W( 167 ) = W( 167 ) + a*JVS( 886 )
  W( 168 ) = W( 168 ) + a*JVS( 887 )
  a = -W( 151 ) / JVS( 894 )
  W( 151 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 895 )
  W( 156 ) = W( 156 ) + a*JVS( 896 )
  W( 157 ) = W( 157 ) + a*JVS( 897 )
  W( 158 ) = W( 158 ) + a*JVS( 898 )
  W( 160 ) = W( 160 ) + a*JVS( 899 )
  W( 165 ) = W( 165 ) + a*JVS( 900 )
  W( 167 ) = W( 167 ) + a*JVS( 901 )
  a = -W( 152 ) / JVS( 912 )
  W( 152 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 913 )
  W( 157 ) = W( 157 ) + a*JVS( 914 )
  W( 160 ) = W( 160 ) + a*JVS( 915 )
  W( 165 ) = W( 165 ) + a*JVS( 916 )
  W( 167 ) = W( 167 ) + a*JVS( 917 )
  a = -W( 153 ) / JVS( 937 )
  W( 153 ) = -a
  W( 154 ) = W( 154 ) + a*JVS( 938 )
  W( 155 ) = W( 155 ) + a*JVS( 939 )
  W( 156 ) = W( 156 ) + a*JVS( 940 )
  W( 157 ) = W( 157 ) + a*JVS( 941 )
  W( 158 ) = W( 158 ) + a*JVS( 942 )
  W( 160 ) = W( 160 ) + a*JVS( 943 )
  W( 161 ) = W( 161 ) + a*JVS( 944 )
  W( 165 ) = W( 165 ) + a*JVS( 945 )
  W( 166 ) = W( 166 ) + a*JVS( 946 )
  W( 167 ) = W( 167 ) + a*JVS( 947 )
  a = -W( 154 ) / JVS( 959 )
  W( 154 ) = -a
  W( 155 ) = W( 155 ) + a*JVS( 960 )
  W( 156 ) = W( 156 ) + a*JVS( 961 )
  W( 157 ) = W( 157 ) + a*JVS( 962 )
  W( 158 ) = W( 158 ) + a*JVS( 963 )
  W( 160 ) = W( 160 ) + a*JVS( 964 )
  W( 165 ) = W( 165 ) + a*JVS( 965 )
  W( 166 ) = W( 166 ) + a*JVS( 966 )
  W( 167 ) = W( 167 ) + a*JVS( 967 )
  W( 168 ) = W( 168 ) + a*JVS( 968 )
  a = -W( 155 ) / JVS( 981 )
  W( 155 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 982 )
  W( 157 ) = W( 157 ) + a*JVS( 983 )
  W( 158 ) = W( 158 ) + a*JVS( 984 )
  W( 160 ) = W( 160 ) + a*JVS( 985 )
  W( 162 ) = W( 162 ) + a*JVS( 986 )
  W( 163 ) = W( 163 ) + a*JVS( 987 )
  W( 164 ) = W( 164 ) + a*JVS( 988 )
  W( 165 ) = W( 165 ) + a*JVS( 989 )
  W( 166 ) = W( 166 ) + a*JVS( 990 )
  W( 167 ) = W( 167 ) + a*JVS( 991 )
  W( 168 ) = W( 168 ) + a*JVS( 992 )
  a = -W( 156 ) / JVS( 1006 )
  W( 156 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 1007 )
  W( 159 ) = W( 159 ) + a*JVS( 1008 )
  W( 160 ) = W( 160 ) + a*JVS( 1009 )
  W( 161 ) = W( 161 ) + a*JVS( 1010 )
  W( 162 ) = W( 162 ) + a*JVS( 1011 )
  W( 163 ) = W( 163 ) + a*JVS( 1012 )
  W( 164 ) = W( 164 ) + a*JVS( 1013 )
  W( 165 ) = W( 165 ) + a*JVS( 1014 )
  W( 167 ) = W( 167 ) + a*JVS( 1015 )
  JVS( 1016) = W( 110 )
  JVS( 1017) = W( 119 )
  JVS( 1018) = W( 124 )
  JVS( 1019) = W( 128 )
  JVS( 1020) = W( 130 )
  JVS( 1021) = W( 131 )
  JVS( 1022) = W( 132 )
  JVS( 1023) = W( 134 )
  JVS( 1024) = W( 135 )
  JVS( 1025) = W( 137 )
  JVS( 1026) = W( 138 )
  JVS( 1027) = W( 139 )
  JVS( 1028) = W( 140 )
  JVS( 1029) = W( 141 )
  JVS( 1030) = W( 142 )
  JVS( 1031) = W( 143 )
  JVS( 1032) = W( 144 )
  JVS( 1033) = W( 145 )
  JVS( 1034) = W( 146 )
  JVS( 1035) = W( 147 )
  JVS( 1036) = W( 148 )
  JVS( 1037) = W( 149 )
  JVS( 1038) = W( 150 )
  JVS( 1039) = W( 151 )
  JVS( 1040) = W( 152 )
  JVS( 1041) = W( 153 )
  JVS( 1042) = W( 154 )
  JVS( 1043) = W( 155 )
  JVS( 1044) = W( 156 )
  JVS( 1045) = W( 157 )
  JVS( 1046) = W( 158 )
  JVS( 1047) = W( 159 )
  JVS( 1048) = W( 160 )
  JVS( 1049) = W( 161 )
  JVS( 1050) = W( 162 )
  JVS( 1051) = W( 163 )
  JVS( 1052) = W( 164 )
  JVS( 1053) = W( 165 )
  JVS( 1054) = W( 166 )
  JVS( 1055) = W( 167 )
  JVS( 1056) = W( 168 )
  IF ( ABS( JVS( 1078 )) < TINY(a) ) THEN
         IER = 158
         RETURN
  END IF
   W( 107 ) = JVS( 1057 )
   W( 112 ) = JVS( 1058 )
   W( 113 ) = JVS( 1059 )
   W( 115 ) = JVS( 1060 )
   W( 121 ) = JVS( 1061 )
   W( 125 ) = JVS( 1062 )
   W( 137 ) = JVS( 1063 )
   W( 139 ) = JVS( 1064 )
   W( 140 ) = JVS( 1065 )
   W( 143 ) = JVS( 1066 )
   W( 144 ) = JVS( 1067 )
   W( 146 ) = JVS( 1068 )
   W( 147 ) = JVS( 1069 )
   W( 148 ) = JVS( 1070 )
   W( 151 ) = JVS( 1071 )
   W( 152 ) = JVS( 1072 )
   W( 153 ) = JVS( 1073 )
   W( 154 ) = JVS( 1074 )
   W( 155 ) = JVS( 1075 )
   W( 156 ) = JVS( 1076 )
   W( 157 ) = JVS( 1077 )
   W( 158 ) = JVS( 1078 )
   W( 159 ) = JVS( 1079 )
   W( 160 ) = JVS( 1080 )
   W( 161 ) = JVS( 1081 )
   W( 162 ) = JVS( 1082 )
   W( 163 ) = JVS( 1083 )
   W( 164 ) = JVS( 1084 )
   W( 165 ) = JVS( 1085 )
   W( 166 ) = JVS( 1086 )
   W( 167 ) = JVS( 1087 )
   W( 168 ) = JVS( 1088 )
  a = -W( 107 ) / JVS( 534 )
  W( 107 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 535 )
  a = -W( 112 ) / JVS( 550 )
  W( 112 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 551 )
  a = -W( 113 ) / JVS( 552 )
  W( 113 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 553 )
  a = -W( 115 ) / JVS( 558 )
  W( 115 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 559 )
  a = -W( 121 ) / JVS( 582 )
  W( 121 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 583 )
  a = -W( 125 ) / JVS( 597 )
  W( 125 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 598 )
  a = -W( 137 ) / JVS( 709 )
  W( 137 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 710 )
  W( 156 ) = W( 156 ) + a*JVS( 711 )
  W( 157 ) = W( 157 ) + a*JVS( 712 )
  W( 165 ) = W( 165 ) + a*JVS( 713 )
  a = -W( 139 ) / JVS( 726 )
  W( 139 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 727 )
  W( 156 ) = W( 156 ) + a*JVS( 728 )
  W( 157 ) = W( 157 ) + a*JVS( 729 )
  W( 165 ) = W( 165 ) + a*JVS( 730 )
  a = -W( 140 ) / JVS( 731 )
  W( 140 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 732 )
  W( 156 ) = W( 156 ) + a*JVS( 733 )
  W( 157 ) = W( 157 ) + a*JVS( 734 )
  W( 165 ) = W( 165 ) + a*JVS( 735 )
  a = -W( 143 ) / JVS( 773 )
  W( 143 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 774 )
  W( 156 ) = W( 156 ) + a*JVS( 775 )
  W( 157 ) = W( 157 ) + a*JVS( 776 )
  W( 165 ) = W( 165 ) + a*JVS( 777 )
  a = -W( 144 ) / JVS( 778 )
  W( 144 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 779 )
  W( 156 ) = W( 156 ) + a*JVS( 780 )
  W( 157 ) = W( 157 ) + a*JVS( 781 )
  W( 165 ) = W( 165 ) + a*JVS( 782 )
  a = -W( 146 ) / JVS( 808 )
  W( 146 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 809 )
  W( 152 ) = W( 152 ) + a*JVS( 810 )
  W( 156 ) = W( 156 ) + a*JVS( 811 )
  W( 157 ) = W( 157 ) + a*JVS( 812 )
  W( 165 ) = W( 165 ) + a*JVS( 813 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  a = -W( 148 ) / JVS( 821 )
  W( 148 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 822 )
  W( 156 ) = W( 156 ) + a*JVS( 823 )
  W( 157 ) = W( 157 ) + a*JVS( 824 )
  W( 165 ) = W( 165 ) + a*JVS( 825 )
  a = -W( 151 ) / JVS( 894 )
  W( 151 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 895 )
  W( 156 ) = W( 156 ) + a*JVS( 896 )
  W( 157 ) = W( 157 ) + a*JVS( 897 )
  W( 158 ) = W( 158 ) + a*JVS( 898 )
  W( 160 ) = W( 160 ) + a*JVS( 899 )
  W( 165 ) = W( 165 ) + a*JVS( 900 )
  W( 167 ) = W( 167 ) + a*JVS( 901 )
  a = -W( 152 ) / JVS( 912 )
  W( 152 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 913 )
  W( 157 ) = W( 157 ) + a*JVS( 914 )
  W( 160 ) = W( 160 ) + a*JVS( 915 )
  W( 165 ) = W( 165 ) + a*JVS( 916 )
  W( 167 ) = W( 167 ) + a*JVS( 917 )
  a = -W( 153 ) / JVS( 937 )
  W( 153 ) = -a
  W( 154 ) = W( 154 ) + a*JVS( 938 )
  W( 155 ) = W( 155 ) + a*JVS( 939 )
  W( 156 ) = W( 156 ) + a*JVS( 940 )
  W( 157 ) = W( 157 ) + a*JVS( 941 )
  W( 158 ) = W( 158 ) + a*JVS( 942 )
  W( 160 ) = W( 160 ) + a*JVS( 943 )
  W( 161 ) = W( 161 ) + a*JVS( 944 )
  W( 165 ) = W( 165 ) + a*JVS( 945 )
  W( 166 ) = W( 166 ) + a*JVS( 946 )
  W( 167 ) = W( 167 ) + a*JVS( 947 )
  a = -W( 154 ) / JVS( 959 )
  W( 154 ) = -a
  W( 155 ) = W( 155 ) + a*JVS( 960 )
  W( 156 ) = W( 156 ) + a*JVS( 961 )
  W( 157 ) = W( 157 ) + a*JVS( 962 )
  W( 158 ) = W( 158 ) + a*JVS( 963 )
  W( 160 ) = W( 160 ) + a*JVS( 964 )
  W( 165 ) = W( 165 ) + a*JVS( 965 )
  W( 166 ) = W( 166 ) + a*JVS( 966 )
  W( 167 ) = W( 167 ) + a*JVS( 967 )
  W( 168 ) = W( 168 ) + a*JVS( 968 )
  a = -W( 155 ) / JVS( 981 )
  W( 155 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 982 )
  W( 157 ) = W( 157 ) + a*JVS( 983 )
  W( 158 ) = W( 158 ) + a*JVS( 984 )
  W( 160 ) = W( 160 ) + a*JVS( 985 )
  W( 162 ) = W( 162 ) + a*JVS( 986 )
  W( 163 ) = W( 163 ) + a*JVS( 987 )
  W( 164 ) = W( 164 ) + a*JVS( 988 )
  W( 165 ) = W( 165 ) + a*JVS( 989 )
  W( 166 ) = W( 166 ) + a*JVS( 990 )
  W( 167 ) = W( 167 ) + a*JVS( 991 )
  W( 168 ) = W( 168 ) + a*JVS( 992 )
  a = -W( 156 ) / JVS( 1006 )
  W( 156 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 1007 )
  W( 159 ) = W( 159 ) + a*JVS( 1008 )
  W( 160 ) = W( 160 ) + a*JVS( 1009 )
  W( 161 ) = W( 161 ) + a*JVS( 1010 )
  W( 162 ) = W( 162 ) + a*JVS( 1011 )
  W( 163 ) = W( 163 ) + a*JVS( 1012 )
  W( 164 ) = W( 164 ) + a*JVS( 1013 )
  W( 165 ) = W( 165 ) + a*JVS( 1014 )
  W( 167 ) = W( 167 ) + a*JVS( 1015 )
  a = -W( 157 ) / JVS( 1045 )
  W( 157 ) = -a
  W( 158 ) = W( 158 ) + a*JVS( 1046 )
  W( 159 ) = W( 159 ) + a*JVS( 1047 )
  W( 160 ) = W( 160 ) + a*JVS( 1048 )
  W( 161 ) = W( 161 ) + a*JVS( 1049 )
  W( 162 ) = W( 162 ) + a*JVS( 1050 )
  W( 163 ) = W( 163 ) + a*JVS( 1051 )
  W( 164 ) = W( 164 ) + a*JVS( 1052 )
  W( 165 ) = W( 165 ) + a*JVS( 1053 )
  W( 166 ) = W( 166 ) + a*JVS( 1054 )
  W( 167 ) = W( 167 ) + a*JVS( 1055 )
  W( 168 ) = W( 168 ) + a*JVS( 1056 )
  JVS( 1057) = W( 107 )
  JVS( 1058) = W( 112 )
  JVS( 1059) = W( 113 )
  JVS( 1060) = W( 115 )
  JVS( 1061) = W( 121 )
  JVS( 1062) = W( 125 )
  JVS( 1063) = W( 137 )
  JVS( 1064) = W( 139 )
  JVS( 1065) = W( 140 )
  JVS( 1066) = W( 143 )
  JVS( 1067) = W( 144 )
  JVS( 1068) = W( 146 )
  JVS( 1069) = W( 147 )
  JVS( 1070) = W( 148 )
  JVS( 1071) = W( 151 )
  JVS( 1072) = W( 152 )
  JVS( 1073) = W( 153 )
  JVS( 1074) = W( 154 )
  JVS( 1075) = W( 155 )
  JVS( 1076) = W( 156 )
  JVS( 1077) = W( 157 )
  JVS( 1078) = W( 158 )
  JVS( 1079) = W( 159 )
  JVS( 1080) = W( 160 )
  JVS( 1081) = W( 161 )
  JVS( 1082) = W( 162 )
  JVS( 1083) = W( 163 )
  JVS( 1084) = W( 164 )
  JVS( 1085) = W( 165 )
  JVS( 1086) = W( 166 )
  JVS( 1087) = W( 167 )
  JVS( 1088) = W( 168 )
  IF ( ABS( JVS( 1099 )) < TINY(a) ) THEN
         IER = 159
         RETURN
  END IF
   W( 103 ) = JVS( 1089 )
   W( 140 ) = JVS( 1090 )
   W( 142 ) = JVS( 1091 )
   W( 146 ) = JVS( 1092 )
   W( 147 ) = JVS( 1093 )
   W( 148 ) = JVS( 1094 )
   W( 152 ) = JVS( 1095 )
   W( 156 ) = JVS( 1096 )
   W( 157 ) = JVS( 1097 )
   W( 158 ) = JVS( 1098 )
   W( 159 ) = JVS( 1099 )
   W( 160 ) = JVS( 1100 )
   W( 161 ) = JVS( 1101 )
   W( 162 ) = JVS( 1102 )
   W( 163 ) = JVS( 1103 )
   W( 164 ) = JVS( 1104 )
   W( 165 ) = JVS( 1105 )
   W( 166 ) = JVS( 1106 )
   W( 167 ) = JVS( 1107 )
   W( 168 ) = JVS( 1108 )
  a = -W( 103 ) / JVS( 522 )
  W( 103 ) = -a
  W( 159 ) = W( 159 ) + a*JVS( 523 )
  W( 167 ) = W( 167 ) + a*JVS( 524 )
  a = -W( 140 ) / JVS( 731 )
  W( 140 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 732 )
  W( 156 ) = W( 156 ) + a*JVS( 733 )
  W( 157 ) = W( 157 ) + a*JVS( 734 )
  W( 165 ) = W( 165 ) + a*JVS( 735 )
  a = -W( 142 ) / JVS( 767 )
  W( 142 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 768 )
  W( 152 ) = W( 152 ) + a*JVS( 769 )
  W( 156 ) = W( 156 ) + a*JVS( 770 )
  W( 157 ) = W( 157 ) + a*JVS( 771 )
  W( 165 ) = W( 165 ) + a*JVS( 772 )
  a = -W( 146 ) / JVS( 808 )
  W( 146 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 809 )
  W( 152 ) = W( 152 ) + a*JVS( 810 )
  W( 156 ) = W( 156 ) + a*JVS( 811 )
  W( 157 ) = W( 157 ) + a*JVS( 812 )
  W( 165 ) = W( 165 ) + a*JVS( 813 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  a = -W( 148 ) / JVS( 821 )
  W( 148 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 822 )
  W( 156 ) = W( 156 ) + a*JVS( 823 )
  W( 157 ) = W( 157 ) + a*JVS( 824 )
  W( 165 ) = W( 165 ) + a*JVS( 825 )
  a = -W( 152 ) / JVS( 912 )
  W( 152 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 913 )
  W( 157 ) = W( 157 ) + a*JVS( 914 )
  W( 160 ) = W( 160 ) + a*JVS( 915 )
  W( 165 ) = W( 165 ) + a*JVS( 916 )
  W( 167 ) = W( 167 ) + a*JVS( 917 )
  a = -W( 156 ) / JVS( 1006 )
  W( 156 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 1007 )
  W( 159 ) = W( 159 ) + a*JVS( 1008 )
  W( 160 ) = W( 160 ) + a*JVS( 1009 )
  W( 161 ) = W( 161 ) + a*JVS( 1010 )
  W( 162 ) = W( 162 ) + a*JVS( 1011 )
  W( 163 ) = W( 163 ) + a*JVS( 1012 )
  W( 164 ) = W( 164 ) + a*JVS( 1013 )
  W( 165 ) = W( 165 ) + a*JVS( 1014 )
  W( 167 ) = W( 167 ) + a*JVS( 1015 )
  a = -W( 157 ) / JVS( 1045 )
  W( 157 ) = -a
  W( 158 ) = W( 158 ) + a*JVS( 1046 )
  W( 159 ) = W( 159 ) + a*JVS( 1047 )
  W( 160 ) = W( 160 ) + a*JVS( 1048 )
  W( 161 ) = W( 161 ) + a*JVS( 1049 )
  W( 162 ) = W( 162 ) + a*JVS( 1050 )
  W( 163 ) = W( 163 ) + a*JVS( 1051 )
  W( 164 ) = W( 164 ) + a*JVS( 1052 )
  W( 165 ) = W( 165 ) + a*JVS( 1053 )
  W( 166 ) = W( 166 ) + a*JVS( 1054 )
  W( 167 ) = W( 167 ) + a*JVS( 1055 )
  W( 168 ) = W( 168 ) + a*JVS( 1056 )
  a = -W( 158 ) / JVS( 1078 )
  W( 158 ) = -a
  W( 159 ) = W( 159 ) + a*JVS( 1079 )
  W( 160 ) = W( 160 ) + a*JVS( 1080 )
  W( 161 ) = W( 161 ) + a*JVS( 1081 )
  W( 162 ) = W( 162 ) + a*JVS( 1082 )
  W( 163 ) = W( 163 ) + a*JVS( 1083 )
  W( 164 ) = W( 164 ) + a*JVS( 1084 )
  W( 165 ) = W( 165 ) + a*JVS( 1085 )
  W( 166 ) = W( 166 ) + a*JVS( 1086 )
  W( 167 ) = W( 167 ) + a*JVS( 1087 )
  W( 168 ) = W( 168 ) + a*JVS( 1088 )
  JVS( 1089) = W( 103 )
  JVS( 1090) = W( 140 )
  JVS( 1091) = W( 142 )
  JVS( 1092) = W( 146 )
  JVS( 1093) = W( 147 )
  JVS( 1094) = W( 148 )
  JVS( 1095) = W( 152 )
  JVS( 1096) = W( 156 )
  JVS( 1097) = W( 157 )
  JVS( 1098) = W( 158 )
  JVS( 1099) = W( 159 )
  JVS( 1100) = W( 160 )
  JVS( 1101) = W( 161 )
  JVS( 1102) = W( 162 )
  JVS( 1103) = W( 163 )
  JVS( 1104) = W( 164 )
  JVS( 1105) = W( 165 )
  JVS( 1106) = W( 166 )
  JVS( 1107) = W( 167 )
  JVS( 1108) = W( 168 )
  IF ( ABS( JVS( 1125 )) < TINY(a) ) THEN
         IER = 160
         RETURN
  END IF
   W( 111 ) = JVS( 1109 )
   W( 117 ) = JVS( 1110 )
   W( 141 ) = JVS( 1111 )
   W( 143 ) = JVS( 1112 )
   W( 144 ) = JVS( 1113 )
   W( 147 ) = JVS( 1114 )
   W( 148 ) = JVS( 1115 )
   W( 150 ) = JVS( 1116 )
   W( 151 ) = JVS( 1117 )
   W( 152 ) = JVS( 1118 )
   W( 154 ) = JVS( 1119 )
   W( 155 ) = JVS( 1120 )
   W( 156 ) = JVS( 1121 )
   W( 157 ) = JVS( 1122 )
   W( 158 ) = JVS( 1123 )
   W( 159 ) = JVS( 1124 )
   W( 160 ) = JVS( 1125 )
   W( 161 ) = JVS( 1126 )
   W( 162 ) = JVS( 1127 )
   W( 163 ) = JVS( 1128 )
   W( 164 ) = JVS( 1129 )
   W( 165 ) = JVS( 1130 )
   W( 166 ) = JVS( 1131 )
   W( 167 ) = JVS( 1132 )
   W( 168 ) = JVS( 1133 )
  a = -W( 111 ) / JVS( 547 )
  W( 111 ) = -a
  W( 160 ) = W( 160 ) + a*JVS( 548 )
  W( 165 ) = W( 165 ) + a*JVS( 549 )
  a = -W( 117 ) / JVS( 564 )
  W( 117 ) = -a
  W( 150 ) = W( 150 ) + a*JVS( 565 )
  W( 160 ) = W( 160 ) + a*JVS( 566 )
  W( 161 ) = W( 161 ) + a*JVS( 567 )
  a = -W( 141 ) / JVS( 744 )
  W( 141 ) = -a
  W( 143 ) = W( 143 ) + a*JVS( 745 )
  W( 144 ) = W( 144 ) + a*JVS( 746 )
  W( 147 ) = W( 147 ) + a*JVS( 747 )
  W( 148 ) = W( 148 ) + a*JVS( 748 )
  W( 151 ) = W( 151 ) + a*JVS( 749 )
  W( 152 ) = W( 152 ) + a*JVS( 750 )
  W( 154 ) = W( 154 ) + a*JVS( 751 )
  W( 155 ) = W( 155 ) + a*JVS( 752 )
  W( 156 ) = W( 156 ) + a*JVS( 753 )
  W( 157 ) = W( 157 ) + a*JVS( 754 )
  W( 158 ) = W( 158 ) + a*JVS( 755 )
  W( 159 ) = W( 159 ) + a*JVS( 756 )
  W( 160 ) = W( 160 ) + a*JVS( 757 )
  W( 161 ) = W( 161 ) + a*JVS( 758 )
  W( 162 ) = W( 162 ) + a*JVS( 759 )
  W( 163 ) = W( 163 ) + a*JVS( 760 )
  W( 164 ) = W( 164 ) + a*JVS( 761 )
  W( 165 ) = W( 165 ) + a*JVS( 762 )
  W( 166 ) = W( 166 ) + a*JVS( 763 )
  W( 167 ) = W( 167 ) + a*JVS( 764 )
  W( 168 ) = W( 168 ) + a*JVS( 765 )
  a = -W( 143 ) / JVS( 773 )
  W( 143 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 774 )
  W( 156 ) = W( 156 ) + a*JVS( 775 )
  W( 157 ) = W( 157 ) + a*JVS( 776 )
  W( 165 ) = W( 165 ) + a*JVS( 777 )
  a = -W( 144 ) / JVS( 778 )
  W( 144 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 779 )
  W( 156 ) = W( 156 ) + a*JVS( 780 )
  W( 157 ) = W( 157 ) + a*JVS( 781 )
  W( 165 ) = W( 165 ) + a*JVS( 782 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  a = -W( 148 ) / JVS( 821 )
  W( 148 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 822 )
  W( 156 ) = W( 156 ) + a*JVS( 823 )
  W( 157 ) = W( 157 ) + a*JVS( 824 )
  W( 165 ) = W( 165 ) + a*JVS( 825 )
  a = -W( 150 ) / JVS( 870 )
  W( 150 ) = -a
  W( 151 ) = W( 151 ) + a*JVS( 871 )
  W( 152 ) = W( 152 ) + a*JVS( 872 )
  W( 154 ) = W( 154 ) + a*JVS( 873 )
  W( 155 ) = W( 155 ) + a*JVS( 874 )
  W( 156 ) = W( 156 ) + a*JVS( 875 )
  W( 157 ) = W( 157 ) + a*JVS( 876 )
  W( 158 ) = W( 158 ) + a*JVS( 877 )
  W( 159 ) = W( 159 ) + a*JVS( 878 )
  W( 160 ) = W( 160 ) + a*JVS( 879 )
  W( 161 ) = W( 161 ) + a*JVS( 880 )
  W( 162 ) = W( 162 ) + a*JVS( 881 )
  W( 163 ) = W( 163 ) + a*JVS( 882 )
  W( 164 ) = W( 164 ) + a*JVS( 883 )
  W( 165 ) = W( 165 ) + a*JVS( 884 )
  W( 166 ) = W( 166 ) + a*JVS( 885 )
  W( 167 ) = W( 167 ) + a*JVS( 886 )
  W( 168 ) = W( 168 ) + a*JVS( 887 )
  a = -W( 151 ) / JVS( 894 )
  W( 151 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 895 )
  W( 156 ) = W( 156 ) + a*JVS( 896 )
  W( 157 ) = W( 157 ) + a*JVS( 897 )
  W( 158 ) = W( 158 ) + a*JVS( 898 )
  W( 160 ) = W( 160 ) + a*JVS( 899 )
  W( 165 ) = W( 165 ) + a*JVS( 900 )
  W( 167 ) = W( 167 ) + a*JVS( 901 )
  a = -W( 152 ) / JVS( 912 )
  W( 152 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 913 )
  W( 157 ) = W( 157 ) + a*JVS( 914 )
  W( 160 ) = W( 160 ) + a*JVS( 915 )
  W( 165 ) = W( 165 ) + a*JVS( 916 )
  W( 167 ) = W( 167 ) + a*JVS( 917 )
  a = -W( 154 ) / JVS( 959 )
  W( 154 ) = -a
  W( 155 ) = W( 155 ) + a*JVS( 960 )
  W( 156 ) = W( 156 ) + a*JVS( 961 )
  W( 157 ) = W( 157 ) + a*JVS( 962 )
  W( 158 ) = W( 158 ) + a*JVS( 963 )
  W( 160 ) = W( 160 ) + a*JVS( 964 )
  W( 165 ) = W( 165 ) + a*JVS( 965 )
  W( 166 ) = W( 166 ) + a*JVS( 966 )
  W( 167 ) = W( 167 ) + a*JVS( 967 )
  W( 168 ) = W( 168 ) + a*JVS( 968 )
  a = -W( 155 ) / JVS( 981 )
  W( 155 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 982 )
  W( 157 ) = W( 157 ) + a*JVS( 983 )
  W( 158 ) = W( 158 ) + a*JVS( 984 )
  W( 160 ) = W( 160 ) + a*JVS( 985 )
  W( 162 ) = W( 162 ) + a*JVS( 986 )
  W( 163 ) = W( 163 ) + a*JVS( 987 )
  W( 164 ) = W( 164 ) + a*JVS( 988 )
  W( 165 ) = W( 165 ) + a*JVS( 989 )
  W( 166 ) = W( 166 ) + a*JVS( 990 )
  W( 167 ) = W( 167 ) + a*JVS( 991 )
  W( 168 ) = W( 168 ) + a*JVS( 992 )
  a = -W( 156 ) / JVS( 1006 )
  W( 156 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 1007 )
  W( 159 ) = W( 159 ) + a*JVS( 1008 )
  W( 160 ) = W( 160 ) + a*JVS( 1009 )
  W( 161 ) = W( 161 ) + a*JVS( 1010 )
  W( 162 ) = W( 162 ) + a*JVS( 1011 )
  W( 163 ) = W( 163 ) + a*JVS( 1012 )
  W( 164 ) = W( 164 ) + a*JVS( 1013 )
  W( 165 ) = W( 165 ) + a*JVS( 1014 )
  W( 167 ) = W( 167 ) + a*JVS( 1015 )
  a = -W( 157 ) / JVS( 1045 )
  W( 157 ) = -a
  W( 158 ) = W( 158 ) + a*JVS( 1046 )
  W( 159 ) = W( 159 ) + a*JVS( 1047 )
  W( 160 ) = W( 160 ) + a*JVS( 1048 )
  W( 161 ) = W( 161 ) + a*JVS( 1049 )
  W( 162 ) = W( 162 ) + a*JVS( 1050 )
  W( 163 ) = W( 163 ) + a*JVS( 1051 )
  W( 164 ) = W( 164 ) + a*JVS( 1052 )
  W( 165 ) = W( 165 ) + a*JVS( 1053 )
  W( 166 ) = W( 166 ) + a*JVS( 1054 )
  W( 167 ) = W( 167 ) + a*JVS( 1055 )
  W( 168 ) = W( 168 ) + a*JVS( 1056 )
  a = -W( 158 ) / JVS( 1078 )
  W( 158 ) = -a
  W( 159 ) = W( 159 ) + a*JVS( 1079 )
  W( 160 ) = W( 160 ) + a*JVS( 1080 )
  W( 161 ) = W( 161 ) + a*JVS( 1081 )
  W( 162 ) = W( 162 ) + a*JVS( 1082 )
  W( 163 ) = W( 163 ) + a*JVS( 1083 )
  W( 164 ) = W( 164 ) + a*JVS( 1084 )
  W( 165 ) = W( 165 ) + a*JVS( 1085 )
  W( 166 ) = W( 166 ) + a*JVS( 1086 )
  W( 167 ) = W( 167 ) + a*JVS( 1087 )
  W( 168 ) = W( 168 ) + a*JVS( 1088 )
  a = -W( 159 ) / JVS( 1099 )
  W( 159 ) = -a
  W( 160 ) = W( 160 ) + a*JVS( 1100 )
  W( 161 ) = W( 161 ) + a*JVS( 1101 )
  W( 162 ) = W( 162 ) + a*JVS( 1102 )
  W( 163 ) = W( 163 ) + a*JVS( 1103 )
  W( 164 ) = W( 164 ) + a*JVS( 1104 )
  W( 165 ) = W( 165 ) + a*JVS( 1105 )
  W( 166 ) = W( 166 ) + a*JVS( 1106 )
  W( 167 ) = W( 167 ) + a*JVS( 1107 )
  W( 168 ) = W( 168 ) + a*JVS( 1108 )
  JVS( 1109) = W( 111 )
  JVS( 1110) = W( 117 )
  JVS( 1111) = W( 141 )
  JVS( 1112) = W( 143 )
  JVS( 1113) = W( 144 )
  JVS( 1114) = W( 147 )
  JVS( 1115) = W( 148 )
  JVS( 1116) = W( 150 )
  JVS( 1117) = W( 151 )
  JVS( 1118) = W( 152 )
  JVS( 1119) = W( 154 )
  JVS( 1120) = W( 155 )
  JVS( 1121) = W( 156 )
  JVS( 1122) = W( 157 )
  JVS( 1123) = W( 158 )
  JVS( 1124) = W( 159 )
  JVS( 1125) = W( 160 )
  JVS( 1126) = W( 161 )
  JVS( 1127) = W( 162 )
  JVS( 1128) = W( 163 )
  JVS( 1129) = W( 164 )
  JVS( 1130) = W( 165 )
  JVS( 1131) = W( 166 )
  JVS( 1132) = W( 167 )
  JVS( 1133) = W( 168 )
  IF ( ABS( JVS( 1177 )) < TINY(a) ) THEN
         IER = 161
         RETURN
  END IF
   W( 100 ) = JVS( 1134 )
   W( 106 ) = JVS( 1135 )
   W( 108 ) = JVS( 1136 )
   W( 111 ) = JVS( 1137 )
   W( 115 ) = JVS( 1138 )
   W( 116 ) = JVS( 1139 )
   W( 117 ) = JVS( 1140 )
   W( 118 ) = JVS( 1141 )
   W( 119 ) = JVS( 1142 )
   W( 120 ) = JVS( 1143 )
   W( 121 ) = JVS( 1144 )
   W( 122 ) = JVS( 1145 )
   W( 123 ) = JVS( 1146 )
   W( 126 ) = JVS( 1147 )
   W( 127 ) = JVS( 1148 )
   W( 128 ) = JVS( 1149 )
   W( 129 ) = JVS( 1150 )
   W( 132 ) = JVS( 1151 )
   W( 133 ) = JVS( 1152 )
   W( 135 ) = JVS( 1153 )
   W( 137 ) = JVS( 1154 )
   W( 138 ) = JVS( 1155 )
   W( 139 ) = JVS( 1156 )
   W( 140 ) = JVS( 1157 )
   W( 142 ) = JVS( 1158 )
   W( 143 ) = JVS( 1159 )
   W( 144 ) = JVS( 1160 )
   W( 145 ) = JVS( 1161 )
   W( 146 ) = JVS( 1162 )
   W( 147 ) = JVS( 1163 )
   W( 148 ) = JVS( 1164 )
   W( 149 ) = JVS( 1165 )
   W( 150 ) = JVS( 1166 )
   W( 151 ) = JVS( 1167 )
   W( 152 ) = JVS( 1168 )
   W( 153 ) = JVS( 1169 )
   W( 154 ) = JVS( 1170 )
   W( 155 ) = JVS( 1171 )
   W( 156 ) = JVS( 1172 )
   W( 157 ) = JVS( 1173 )
   W( 158 ) = JVS( 1174 )
   W( 159 ) = JVS( 1175 )
   W( 160 ) = JVS( 1176 )
   W( 161 ) = JVS( 1177 )
   W( 162 ) = JVS( 1178 )
   W( 163 ) = JVS( 1179 )
   W( 164 ) = JVS( 1180 )
   W( 165 ) = JVS( 1181 )
   W( 166 ) = JVS( 1182 )
   W( 167 ) = JVS( 1183 )
   W( 168 ) = JVS( 1184 )
  a = -W( 100 ) / JVS( 515 )
  W( 100 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 516 )
  a = -W( 106 ) / JVS( 531 )
  W( 106 ) = -a
  W( 161 ) = W( 161 ) + a*JVS( 532 )
  W( 165 ) = W( 165 ) + a*JVS( 533 )
  a = -W( 108 ) / JVS( 536 )
  W( 108 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 537 )
  a = -W( 111 ) / JVS( 547 )
  W( 111 ) = -a
  W( 160 ) = W( 160 ) + a*JVS( 548 )
  W( 165 ) = W( 165 ) + a*JVS( 549 )
  a = -W( 115 ) / JVS( 558 )
  W( 115 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 559 )
  a = -W( 116 ) / JVS( 560 )
  W( 116 ) = -a
  W( 161 ) = W( 161 ) + a*JVS( 561 )
  W( 165 ) = W( 165 ) + a*JVS( 562 )
  W( 168 ) = W( 168 ) + a*JVS( 563 )
  a = -W( 117 ) / JVS( 564 )
  W( 117 ) = -a
  W( 150 ) = W( 150 ) + a*JVS( 565 )
  W( 160 ) = W( 160 ) + a*JVS( 566 )
  W( 161 ) = W( 161 ) + a*JVS( 567 )
  a = -W( 118 ) / JVS( 568 )
  W( 118 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 569 )
  W( 157 ) = W( 157 ) + a*JVS( 570 )
  W( 161 ) = W( 161 ) + a*JVS( 571 )
  W( 167 ) = W( 167 ) + a*JVS( 572 )
  a = -W( 119 ) / JVS( 573 )
  W( 119 ) = -a
  W( 161 ) = W( 161 ) + a*JVS( 574 )
  W( 165 ) = W( 165 ) + a*JVS( 575 )
  W( 167 ) = W( 167 ) + a*JVS( 576 )
  a = -W( 120 ) / JVS( 577 )
  W( 120 ) = -a
  W( 158 ) = W( 158 ) + a*JVS( 578 )
  W( 165 ) = W( 165 ) + a*JVS( 579 )
  W( 166 ) = W( 166 ) + a*JVS( 580 )
  W( 168 ) = W( 168 ) + a*JVS( 581 )
  a = -W( 121 ) / JVS( 582 )
  W( 121 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 583 )
  a = -W( 122 ) / JVS( 586 )
  W( 122 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 587 )
  a = -W( 123 ) / JVS( 590 )
  W( 123 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 591 )
  a = -W( 126 ) / JVS( 599 )
  W( 126 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 600 )
  W( 165 ) = W( 165 ) + a*JVS( 601 )
  a = -W( 127 ) / JVS( 604 )
  W( 127 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 605 )
  W( 165 ) = W( 165 ) + a*JVS( 606 )
  a = -W( 128 ) / JVS( 608 )
  W( 128 ) = -a
  W( 138 ) = W( 138 ) + a*JVS( 609 )
  W( 157 ) = W( 157 ) + a*JVS( 610 )
  W( 161 ) = W( 161 ) + a*JVS( 611 )
  W( 167 ) = W( 167 ) + a*JVS( 612 )
  a = -W( 129 ) / JVS( 613 )
  W( 129 ) = -a
  W( 158 ) = W( 158 ) + a*JVS( 614 )
  W( 161 ) = W( 161 ) + a*JVS( 615 )
  W( 165 ) = W( 165 ) + a*JVS( 616 )
  W( 166 ) = W( 166 ) + a*JVS( 617 )
  a = -W( 132 ) / JVS( 636 )
  W( 132 ) = -a
  W( 142 ) = W( 142 ) + a*JVS( 637 )
  W( 146 ) = W( 146 ) + a*JVS( 638 )
  W( 148 ) = W( 148 ) + a*JVS( 639 )
  W( 156 ) = W( 156 ) + a*JVS( 640 )
  W( 157 ) = W( 157 ) + a*JVS( 641 )
  W( 165 ) = W( 165 ) + a*JVS( 642 )
  a = -W( 133 ) / JVS( 649 )
  W( 133 ) = -a
  W( 135 ) = W( 135 ) + a*JVS( 650 )
  W( 137 ) = W( 137 ) + a*JVS( 651 )
  W( 139 ) = W( 139 ) + a*JVS( 652 )
  W( 140 ) = W( 140 ) + a*JVS( 653 )
  W( 142 ) = W( 142 ) + a*JVS( 654 )
  W( 143 ) = W( 143 ) + a*JVS( 655 )
  W( 144 ) = W( 144 ) + a*JVS( 656 )
  W( 145 ) = W( 145 ) + a*JVS( 657 )
  W( 146 ) = W( 146 ) + a*JVS( 658 )
  W( 147 ) = W( 147 ) + a*JVS( 659 )
  W( 148 ) = W( 148 ) + a*JVS( 660 )
  W( 149 ) = W( 149 ) + a*JVS( 661 )
  W( 150 ) = W( 150 ) + a*JVS( 662 )
  W( 152 ) = W( 152 ) + a*JVS( 663 )
  W( 153 ) = W( 153 ) + a*JVS( 664 )
  W( 156 ) = W( 156 ) + a*JVS( 665 )
  W( 157 ) = W( 157 ) + a*JVS( 666 )
  W( 165 ) = W( 165 ) + a*JVS( 667 )
  a = -W( 135 ) / JVS( 689 )
  W( 135 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 690 )
  W( 156 ) = W( 156 ) + a*JVS( 691 )
  W( 157 ) = W( 157 ) + a*JVS( 692 )
  W( 165 ) = W( 165 ) + a*JVS( 693 )
  a = -W( 137 ) / JVS( 709 )
  W( 137 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 710 )
  W( 156 ) = W( 156 ) + a*JVS( 711 )
  W( 157 ) = W( 157 ) + a*JVS( 712 )
  W( 165 ) = W( 165 ) + a*JVS( 713 )
  a = -W( 138 ) / JVS( 716 )
  W( 138 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 717 )
  W( 159 ) = W( 159 ) + a*JVS( 718 )
  W( 160 ) = W( 160 ) + a*JVS( 719 )
  W( 161 ) = W( 161 ) + a*JVS( 720 )
  W( 162 ) = W( 162 ) + a*JVS( 721 )
  W( 163 ) = W( 163 ) + a*JVS( 722 )
  W( 164 ) = W( 164 ) + a*JVS( 723 )
  W( 165 ) = W( 165 ) + a*JVS( 724 )
  W( 167 ) = W( 167 ) + a*JVS( 725 )
  a = -W( 139 ) / JVS( 726 )
  W( 139 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 727 )
  W( 156 ) = W( 156 ) + a*JVS( 728 )
  W( 157 ) = W( 157 ) + a*JVS( 729 )
  W( 165 ) = W( 165 ) + a*JVS( 730 )
  a = -W( 140 ) / JVS( 731 )
  W( 140 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 732 )
  W( 156 ) = W( 156 ) + a*JVS( 733 )
  W( 157 ) = W( 157 ) + a*JVS( 734 )
  W( 165 ) = W( 165 ) + a*JVS( 735 )
  a = -W( 142 ) / JVS( 767 )
  W( 142 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 768 )
  W( 152 ) = W( 152 ) + a*JVS( 769 )
  W( 156 ) = W( 156 ) + a*JVS( 770 )
  W( 157 ) = W( 157 ) + a*JVS( 771 )
  W( 165 ) = W( 165 ) + a*JVS( 772 )
  a = -W( 143 ) / JVS( 773 )
  W( 143 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 774 )
  W( 156 ) = W( 156 ) + a*JVS( 775 )
  W( 157 ) = W( 157 ) + a*JVS( 776 )
  W( 165 ) = W( 165 ) + a*JVS( 777 )
  a = -W( 144 ) / JVS( 778 )
  W( 144 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 779 )
  W( 156 ) = W( 156 ) + a*JVS( 780 )
  W( 157 ) = W( 157 ) + a*JVS( 781 )
  W( 165 ) = W( 165 ) + a*JVS( 782 )
  a = -W( 145 ) / JVS( 794 )
  W( 145 ) = -a
  W( 146 ) = W( 146 ) + a*JVS( 795 )
  W( 152 ) = W( 152 ) + a*JVS( 796 )
  W( 156 ) = W( 156 ) + a*JVS( 797 )
  W( 157 ) = W( 157 ) + a*JVS( 798 )
  W( 159 ) = W( 159 ) + a*JVS( 799 )
  W( 160 ) = W( 160 ) + a*JVS( 800 )
  W( 161 ) = W( 161 ) + a*JVS( 801 )
  W( 162 ) = W( 162 ) + a*JVS( 802 )
  W( 163 ) = W( 163 ) + a*JVS( 803 )
  W( 164 ) = W( 164 ) + a*JVS( 804 )
  W( 165 ) = W( 165 ) + a*JVS( 805 )
  W( 167 ) = W( 167 ) + a*JVS( 806 )
  a = -W( 146 ) / JVS( 808 )
  W( 146 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 809 )
  W( 152 ) = W( 152 ) + a*JVS( 810 )
  W( 156 ) = W( 156 ) + a*JVS( 811 )
  W( 157 ) = W( 157 ) + a*JVS( 812 )
  W( 165 ) = W( 165 ) + a*JVS( 813 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  a = -W( 148 ) / JVS( 821 )
  W( 148 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 822 )
  W( 156 ) = W( 156 ) + a*JVS( 823 )
  W( 157 ) = W( 157 ) + a*JVS( 824 )
  W( 165 ) = W( 165 ) + a*JVS( 825 )
  a = -W( 149 ) / JVS( 836 )
  W( 149 ) = -a
  W( 151 ) = W( 151 ) + a*JVS( 837 )
  W( 152 ) = W( 152 ) + a*JVS( 838 )
  W( 153 ) = W( 153 ) + a*JVS( 839 )
  W( 154 ) = W( 154 ) + a*JVS( 840 )
  W( 155 ) = W( 155 ) + a*JVS( 841 )
  W( 156 ) = W( 156 ) + a*JVS( 842 )
  W( 157 ) = W( 157 ) + a*JVS( 843 )
  W( 159 ) = W( 159 ) + a*JVS( 844 )
  W( 160 ) = W( 160 ) + a*JVS( 845 )
  W( 162 ) = W( 162 ) + a*JVS( 846 )
  W( 163 ) = W( 163 ) + a*JVS( 847 )
  W( 164 ) = W( 164 ) + a*JVS( 848 )
  W( 165 ) = W( 165 ) + a*JVS( 849 )
  a = -W( 150 ) / JVS( 870 )
  W( 150 ) = -a
  W( 151 ) = W( 151 ) + a*JVS( 871 )
  W( 152 ) = W( 152 ) + a*JVS( 872 )
  W( 154 ) = W( 154 ) + a*JVS( 873 )
  W( 155 ) = W( 155 ) + a*JVS( 874 )
  W( 156 ) = W( 156 ) + a*JVS( 875 )
  W( 157 ) = W( 157 ) + a*JVS( 876 )
  W( 158 ) = W( 158 ) + a*JVS( 877 )
  W( 159 ) = W( 159 ) + a*JVS( 878 )
  W( 160 ) = W( 160 ) + a*JVS( 879 )
  W( 161 ) = W( 161 ) + a*JVS( 880 )
  W( 162 ) = W( 162 ) + a*JVS( 881 )
  W( 163 ) = W( 163 ) + a*JVS( 882 )
  W( 164 ) = W( 164 ) + a*JVS( 883 )
  W( 165 ) = W( 165 ) + a*JVS( 884 )
  W( 166 ) = W( 166 ) + a*JVS( 885 )
  W( 167 ) = W( 167 ) + a*JVS( 886 )
  W( 168 ) = W( 168 ) + a*JVS( 887 )
  a = -W( 151 ) / JVS( 894 )
  W( 151 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 895 )
  W( 156 ) = W( 156 ) + a*JVS( 896 )
  W( 157 ) = W( 157 ) + a*JVS( 897 )
  W( 158 ) = W( 158 ) + a*JVS( 898 )
  W( 160 ) = W( 160 ) + a*JVS( 899 )
  W( 165 ) = W( 165 ) + a*JVS( 900 )
  W( 167 ) = W( 167 ) + a*JVS( 901 )
  a = -W( 152 ) / JVS( 912 )
  W( 152 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 913 )
  W( 157 ) = W( 157 ) + a*JVS( 914 )
  W( 160 ) = W( 160 ) + a*JVS( 915 )
  W( 165 ) = W( 165 ) + a*JVS( 916 )
  W( 167 ) = W( 167 ) + a*JVS( 917 )
  a = -W( 153 ) / JVS( 937 )
  W( 153 ) = -a
  W( 154 ) = W( 154 ) + a*JVS( 938 )
  W( 155 ) = W( 155 ) + a*JVS( 939 )
  W( 156 ) = W( 156 ) + a*JVS( 940 )
  W( 157 ) = W( 157 ) + a*JVS( 941 )
  W( 158 ) = W( 158 ) + a*JVS( 942 )
  W( 160 ) = W( 160 ) + a*JVS( 943 )
  W( 161 ) = W( 161 ) + a*JVS( 944 )
  W( 165 ) = W( 165 ) + a*JVS( 945 )
  W( 166 ) = W( 166 ) + a*JVS( 946 )
  W( 167 ) = W( 167 ) + a*JVS( 947 )
  a = -W( 154 ) / JVS( 959 )
  W( 154 ) = -a
  W( 155 ) = W( 155 ) + a*JVS( 960 )
  W( 156 ) = W( 156 ) + a*JVS( 961 )
  W( 157 ) = W( 157 ) + a*JVS( 962 )
  W( 158 ) = W( 158 ) + a*JVS( 963 )
  W( 160 ) = W( 160 ) + a*JVS( 964 )
  W( 165 ) = W( 165 ) + a*JVS( 965 )
  W( 166 ) = W( 166 ) + a*JVS( 966 )
  W( 167 ) = W( 167 ) + a*JVS( 967 )
  W( 168 ) = W( 168 ) + a*JVS( 968 )
  a = -W( 155 ) / JVS( 981 )
  W( 155 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 982 )
  W( 157 ) = W( 157 ) + a*JVS( 983 )
  W( 158 ) = W( 158 ) + a*JVS( 984 )
  W( 160 ) = W( 160 ) + a*JVS( 985 )
  W( 162 ) = W( 162 ) + a*JVS( 986 )
  W( 163 ) = W( 163 ) + a*JVS( 987 )
  W( 164 ) = W( 164 ) + a*JVS( 988 )
  W( 165 ) = W( 165 ) + a*JVS( 989 )
  W( 166 ) = W( 166 ) + a*JVS( 990 )
  W( 167 ) = W( 167 ) + a*JVS( 991 )
  W( 168 ) = W( 168 ) + a*JVS( 992 )
  a = -W( 156 ) / JVS( 1006 )
  W( 156 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 1007 )
  W( 159 ) = W( 159 ) + a*JVS( 1008 )
  W( 160 ) = W( 160 ) + a*JVS( 1009 )
  W( 161 ) = W( 161 ) + a*JVS( 1010 )
  W( 162 ) = W( 162 ) + a*JVS( 1011 )
  W( 163 ) = W( 163 ) + a*JVS( 1012 )
  W( 164 ) = W( 164 ) + a*JVS( 1013 )
  W( 165 ) = W( 165 ) + a*JVS( 1014 )
  W( 167 ) = W( 167 ) + a*JVS( 1015 )
  a = -W( 157 ) / JVS( 1045 )
  W( 157 ) = -a
  W( 158 ) = W( 158 ) + a*JVS( 1046 )
  W( 159 ) = W( 159 ) + a*JVS( 1047 )
  W( 160 ) = W( 160 ) + a*JVS( 1048 )
  W( 161 ) = W( 161 ) + a*JVS( 1049 )
  W( 162 ) = W( 162 ) + a*JVS( 1050 )
  W( 163 ) = W( 163 ) + a*JVS( 1051 )
  W( 164 ) = W( 164 ) + a*JVS( 1052 )
  W( 165 ) = W( 165 ) + a*JVS( 1053 )
  W( 166 ) = W( 166 ) + a*JVS( 1054 )
  W( 167 ) = W( 167 ) + a*JVS( 1055 )
  W( 168 ) = W( 168 ) + a*JVS( 1056 )
  a = -W( 158 ) / JVS( 1078 )
  W( 158 ) = -a
  W( 159 ) = W( 159 ) + a*JVS( 1079 )
  W( 160 ) = W( 160 ) + a*JVS( 1080 )
  W( 161 ) = W( 161 ) + a*JVS( 1081 )
  W( 162 ) = W( 162 ) + a*JVS( 1082 )
  W( 163 ) = W( 163 ) + a*JVS( 1083 )
  W( 164 ) = W( 164 ) + a*JVS( 1084 )
  W( 165 ) = W( 165 ) + a*JVS( 1085 )
  W( 166 ) = W( 166 ) + a*JVS( 1086 )
  W( 167 ) = W( 167 ) + a*JVS( 1087 )
  W( 168 ) = W( 168 ) + a*JVS( 1088 )
  a = -W( 159 ) / JVS( 1099 )
  W( 159 ) = -a
  W( 160 ) = W( 160 ) + a*JVS( 1100 )
  W( 161 ) = W( 161 ) + a*JVS( 1101 )
  W( 162 ) = W( 162 ) + a*JVS( 1102 )
  W( 163 ) = W( 163 ) + a*JVS( 1103 )
  W( 164 ) = W( 164 ) + a*JVS( 1104 )
  W( 165 ) = W( 165 ) + a*JVS( 1105 )
  W( 166 ) = W( 166 ) + a*JVS( 1106 )
  W( 167 ) = W( 167 ) + a*JVS( 1107 )
  W( 168 ) = W( 168 ) + a*JVS( 1108 )
  a = -W( 160 ) / JVS( 1125 )
  W( 160 ) = -a
  W( 161 ) = W( 161 ) + a*JVS( 1126 )
  W( 162 ) = W( 162 ) + a*JVS( 1127 )
  W( 163 ) = W( 163 ) + a*JVS( 1128 )
  W( 164 ) = W( 164 ) + a*JVS( 1129 )
  W( 165 ) = W( 165 ) + a*JVS( 1130 )
  W( 166 ) = W( 166 ) + a*JVS( 1131 )
  W( 167 ) = W( 167 ) + a*JVS( 1132 )
  W( 168 ) = W( 168 ) + a*JVS( 1133 )
  JVS( 1134) = W( 100 )
  JVS( 1135) = W( 106 )
  JVS( 1136) = W( 108 )
  JVS( 1137) = W( 111 )
  JVS( 1138) = W( 115 )
  JVS( 1139) = W( 116 )
  JVS( 1140) = W( 117 )
  JVS( 1141) = W( 118 )
  JVS( 1142) = W( 119 )
  JVS( 1143) = W( 120 )
  JVS( 1144) = W( 121 )
  JVS( 1145) = W( 122 )
  JVS( 1146) = W( 123 )
  JVS( 1147) = W( 126 )
  JVS( 1148) = W( 127 )
  JVS( 1149) = W( 128 )
  JVS( 1150) = W( 129 )
  JVS( 1151) = W( 132 )
  JVS( 1152) = W( 133 )
  JVS( 1153) = W( 135 )
  JVS( 1154) = W( 137 )
  JVS( 1155) = W( 138 )
  JVS( 1156) = W( 139 )
  JVS( 1157) = W( 140 )
  JVS( 1158) = W( 142 )
  JVS( 1159) = W( 143 )
  JVS( 1160) = W( 144 )
  JVS( 1161) = W( 145 )
  JVS( 1162) = W( 146 )
  JVS( 1163) = W( 147 )
  JVS( 1164) = W( 148 )
  JVS( 1165) = W( 149 )
  JVS( 1166) = W( 150 )
  JVS( 1167) = W( 151 )
  JVS( 1168) = W( 152 )
  JVS( 1169) = W( 153 )
  JVS( 1170) = W( 154 )
  JVS( 1171) = W( 155 )
  JVS( 1172) = W( 156 )
  JVS( 1173) = W( 157 )
  JVS( 1174) = W( 158 )
  JVS( 1175) = W( 159 )
  JVS( 1176) = W( 160 )
  JVS( 1177) = W( 161 )
  JVS( 1178) = W( 162 )
  JVS( 1179) = W( 163 )
  JVS( 1180) = W( 164 )
  JVS( 1181) = W( 165 )
  JVS( 1182) = W( 166 )
  JVS( 1183) = W( 167 )
  JVS( 1184) = W( 168 )
  IF ( ABS( JVS( 1203 )) < TINY(a) ) THEN
         IER = 162
         RETURN
  END IF
   W( 102 ) = JVS( 1185 )
   W( 142 ) = JVS( 1186 )
   W( 143 ) = JVS( 1187 )
   W( 144 ) = JVS( 1188 )
   W( 145 ) = JVS( 1189 )
   W( 146 ) = JVS( 1190 )
   W( 147 ) = JVS( 1191 )
   W( 148 ) = JVS( 1192 )
   W( 152 ) = JVS( 1193 )
   W( 153 ) = JVS( 1194 )
   W( 154 ) = JVS( 1195 )
   W( 155 ) = JVS( 1196 )
   W( 156 ) = JVS( 1197 )
   W( 157 ) = JVS( 1198 )
   W( 158 ) = JVS( 1199 )
   W( 159 ) = JVS( 1200 )
   W( 160 ) = JVS( 1201 )
   W( 161 ) = JVS( 1202 )
   W( 162 ) = JVS( 1203 )
   W( 163 ) = JVS( 1204 )
   W( 164 ) = JVS( 1205 )
   W( 165 ) = JVS( 1206 )
   W( 166 ) = JVS( 1207 )
   W( 167 ) = JVS( 1208 )
   W( 168 ) = JVS( 1209 )
  a = -W( 102 ) / JVS( 519 )
  W( 102 ) = -a
  W( 162 ) = W( 162 ) + a*JVS( 520 )
  W( 167 ) = W( 167 ) + a*JVS( 521 )
  a = -W( 142 ) / JVS( 767 )
  W( 142 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 768 )
  W( 152 ) = W( 152 ) + a*JVS( 769 )
  W( 156 ) = W( 156 ) + a*JVS( 770 )
  W( 157 ) = W( 157 ) + a*JVS( 771 )
  W( 165 ) = W( 165 ) + a*JVS( 772 )
  a = -W( 143 ) / JVS( 773 )
  W( 143 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 774 )
  W( 156 ) = W( 156 ) + a*JVS( 775 )
  W( 157 ) = W( 157 ) + a*JVS( 776 )
  W( 165 ) = W( 165 ) + a*JVS( 777 )
  a = -W( 144 ) / JVS( 778 )
  W( 144 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 779 )
  W( 156 ) = W( 156 ) + a*JVS( 780 )
  W( 157 ) = W( 157 ) + a*JVS( 781 )
  W( 165 ) = W( 165 ) + a*JVS( 782 )
  a = -W( 145 ) / JVS( 794 )
  W( 145 ) = -a
  W( 146 ) = W( 146 ) + a*JVS( 795 )
  W( 152 ) = W( 152 ) + a*JVS( 796 )
  W( 156 ) = W( 156 ) + a*JVS( 797 )
  W( 157 ) = W( 157 ) + a*JVS( 798 )
  W( 159 ) = W( 159 ) + a*JVS( 799 )
  W( 160 ) = W( 160 ) + a*JVS( 800 )
  W( 161 ) = W( 161 ) + a*JVS( 801 )
  W( 162 ) = W( 162 ) + a*JVS( 802 )
  W( 163 ) = W( 163 ) + a*JVS( 803 )
  W( 164 ) = W( 164 ) + a*JVS( 804 )
  W( 165 ) = W( 165 ) + a*JVS( 805 )
  W( 167 ) = W( 167 ) + a*JVS( 806 )
  a = -W( 146 ) / JVS( 808 )
  W( 146 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 809 )
  W( 152 ) = W( 152 ) + a*JVS( 810 )
  W( 156 ) = W( 156 ) + a*JVS( 811 )
  W( 157 ) = W( 157 ) + a*JVS( 812 )
  W( 165 ) = W( 165 ) + a*JVS( 813 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  a = -W( 148 ) / JVS( 821 )
  W( 148 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 822 )
  W( 156 ) = W( 156 ) + a*JVS( 823 )
  W( 157 ) = W( 157 ) + a*JVS( 824 )
  W( 165 ) = W( 165 ) + a*JVS( 825 )
  a = -W( 152 ) / JVS( 912 )
  W( 152 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 913 )
  W( 157 ) = W( 157 ) + a*JVS( 914 )
  W( 160 ) = W( 160 ) + a*JVS( 915 )
  W( 165 ) = W( 165 ) + a*JVS( 916 )
  W( 167 ) = W( 167 ) + a*JVS( 917 )
  a = -W( 153 ) / JVS( 937 )
  W( 153 ) = -a
  W( 154 ) = W( 154 ) + a*JVS( 938 )
  W( 155 ) = W( 155 ) + a*JVS( 939 )
  W( 156 ) = W( 156 ) + a*JVS( 940 )
  W( 157 ) = W( 157 ) + a*JVS( 941 )
  W( 158 ) = W( 158 ) + a*JVS( 942 )
  W( 160 ) = W( 160 ) + a*JVS( 943 )
  W( 161 ) = W( 161 ) + a*JVS( 944 )
  W( 165 ) = W( 165 ) + a*JVS( 945 )
  W( 166 ) = W( 166 ) + a*JVS( 946 )
  W( 167 ) = W( 167 ) + a*JVS( 947 )
  a = -W( 154 ) / JVS( 959 )
  W( 154 ) = -a
  W( 155 ) = W( 155 ) + a*JVS( 960 )
  W( 156 ) = W( 156 ) + a*JVS( 961 )
  W( 157 ) = W( 157 ) + a*JVS( 962 )
  W( 158 ) = W( 158 ) + a*JVS( 963 )
  W( 160 ) = W( 160 ) + a*JVS( 964 )
  W( 165 ) = W( 165 ) + a*JVS( 965 )
  W( 166 ) = W( 166 ) + a*JVS( 966 )
  W( 167 ) = W( 167 ) + a*JVS( 967 )
  W( 168 ) = W( 168 ) + a*JVS( 968 )
  a = -W( 155 ) / JVS( 981 )
  W( 155 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 982 )
  W( 157 ) = W( 157 ) + a*JVS( 983 )
  W( 158 ) = W( 158 ) + a*JVS( 984 )
  W( 160 ) = W( 160 ) + a*JVS( 985 )
  W( 162 ) = W( 162 ) + a*JVS( 986 )
  W( 163 ) = W( 163 ) + a*JVS( 987 )
  W( 164 ) = W( 164 ) + a*JVS( 988 )
  W( 165 ) = W( 165 ) + a*JVS( 989 )
  W( 166 ) = W( 166 ) + a*JVS( 990 )
  W( 167 ) = W( 167 ) + a*JVS( 991 )
  W( 168 ) = W( 168 ) + a*JVS( 992 )
  a = -W( 156 ) / JVS( 1006 )
  W( 156 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 1007 )
  W( 159 ) = W( 159 ) + a*JVS( 1008 )
  W( 160 ) = W( 160 ) + a*JVS( 1009 )
  W( 161 ) = W( 161 ) + a*JVS( 1010 )
  W( 162 ) = W( 162 ) + a*JVS( 1011 )
  W( 163 ) = W( 163 ) + a*JVS( 1012 )
  W( 164 ) = W( 164 ) + a*JVS( 1013 )
  W( 165 ) = W( 165 ) + a*JVS( 1014 )
  W( 167 ) = W( 167 ) + a*JVS( 1015 )
  a = -W( 157 ) / JVS( 1045 )
  W( 157 ) = -a
  W( 158 ) = W( 158 ) + a*JVS( 1046 )
  W( 159 ) = W( 159 ) + a*JVS( 1047 )
  W( 160 ) = W( 160 ) + a*JVS( 1048 )
  W( 161 ) = W( 161 ) + a*JVS( 1049 )
  W( 162 ) = W( 162 ) + a*JVS( 1050 )
  W( 163 ) = W( 163 ) + a*JVS( 1051 )
  W( 164 ) = W( 164 ) + a*JVS( 1052 )
  W( 165 ) = W( 165 ) + a*JVS( 1053 )
  W( 166 ) = W( 166 ) + a*JVS( 1054 )
  W( 167 ) = W( 167 ) + a*JVS( 1055 )
  W( 168 ) = W( 168 ) + a*JVS( 1056 )
  a = -W( 158 ) / JVS( 1078 )
  W( 158 ) = -a
  W( 159 ) = W( 159 ) + a*JVS( 1079 )
  W( 160 ) = W( 160 ) + a*JVS( 1080 )
  W( 161 ) = W( 161 ) + a*JVS( 1081 )
  W( 162 ) = W( 162 ) + a*JVS( 1082 )
  W( 163 ) = W( 163 ) + a*JVS( 1083 )
  W( 164 ) = W( 164 ) + a*JVS( 1084 )
  W( 165 ) = W( 165 ) + a*JVS( 1085 )
  W( 166 ) = W( 166 ) + a*JVS( 1086 )
  W( 167 ) = W( 167 ) + a*JVS( 1087 )
  W( 168 ) = W( 168 ) + a*JVS( 1088 )
  a = -W( 159 ) / JVS( 1099 )
  W( 159 ) = -a
  W( 160 ) = W( 160 ) + a*JVS( 1100 )
  W( 161 ) = W( 161 ) + a*JVS( 1101 )
  W( 162 ) = W( 162 ) + a*JVS( 1102 )
  W( 163 ) = W( 163 ) + a*JVS( 1103 )
  W( 164 ) = W( 164 ) + a*JVS( 1104 )
  W( 165 ) = W( 165 ) + a*JVS( 1105 )
  W( 166 ) = W( 166 ) + a*JVS( 1106 )
  W( 167 ) = W( 167 ) + a*JVS( 1107 )
  W( 168 ) = W( 168 ) + a*JVS( 1108 )
  a = -W( 160 ) / JVS( 1125 )
  W( 160 ) = -a
  W( 161 ) = W( 161 ) + a*JVS( 1126 )
  W( 162 ) = W( 162 ) + a*JVS( 1127 )
  W( 163 ) = W( 163 ) + a*JVS( 1128 )
  W( 164 ) = W( 164 ) + a*JVS( 1129 )
  W( 165 ) = W( 165 ) + a*JVS( 1130 )
  W( 166 ) = W( 166 ) + a*JVS( 1131 )
  W( 167 ) = W( 167 ) + a*JVS( 1132 )
  W( 168 ) = W( 168 ) + a*JVS( 1133 )
  a = -W( 161 ) / JVS( 1177 )
  W( 161 ) = -a
  W( 162 ) = W( 162 ) + a*JVS( 1178 )
  W( 163 ) = W( 163 ) + a*JVS( 1179 )
  W( 164 ) = W( 164 ) + a*JVS( 1180 )
  W( 165 ) = W( 165 ) + a*JVS( 1181 )
  W( 166 ) = W( 166 ) + a*JVS( 1182 )
  W( 167 ) = W( 167 ) + a*JVS( 1183 )
  W( 168 ) = W( 168 ) + a*JVS( 1184 )
  JVS( 1185) = W( 102 )
  JVS( 1186) = W( 142 )
  JVS( 1187) = W( 143 )
  JVS( 1188) = W( 144 )
  JVS( 1189) = W( 145 )
  JVS( 1190) = W( 146 )
  JVS( 1191) = W( 147 )
  JVS( 1192) = W( 148 )
  JVS( 1193) = W( 152 )
  JVS( 1194) = W( 153 )
  JVS( 1195) = W( 154 )
  JVS( 1196) = W( 155 )
  JVS( 1197) = W( 156 )
  JVS( 1198) = W( 157 )
  JVS( 1199) = W( 158 )
  JVS( 1200) = W( 159 )
  JVS( 1201) = W( 160 )
  JVS( 1202) = W( 161 )
  JVS( 1203) = W( 162 )
  JVS( 1204) = W( 163 )
  JVS( 1205) = W( 164 )
  JVS( 1206) = W( 165 )
  JVS( 1207) = W( 166 )
  JVS( 1208) = W( 167 )
  JVS( 1209) = W( 168 )
  IF ( ABS( JVS( 1238 )) < TINY(a) ) THEN
         IER = 163
         RETURN
  END IF
   W( 104 ) = JVS( 1210 )
   W( 109 ) = JVS( 1211 )
   W( 115 ) = JVS( 1212 )
   W( 122 ) = JVS( 1213 )
   W( 123 ) = JVS( 1214 )
   W( 125 ) = JVS( 1215 )
   W( 132 ) = JVS( 1216 )
   W( 136 ) = JVS( 1217 )
   W( 139 ) = JVS( 1218 )
   W( 142 ) = JVS( 1219 )
   W( 143 ) = JVS( 1220 )
   W( 144 ) = JVS( 1221 )
   W( 146 ) = JVS( 1222 )
   W( 147 ) = JVS( 1223 )
   W( 148 ) = JVS( 1224 )
   W( 149 ) = JVS( 1225 )
   W( 151 ) = JVS( 1226 )
   W( 152 ) = JVS( 1227 )
   W( 153 ) = JVS( 1228 )
   W( 154 ) = JVS( 1229 )
   W( 155 ) = JVS( 1230 )
   W( 156 ) = JVS( 1231 )
   W( 157 ) = JVS( 1232 )
   W( 158 ) = JVS( 1233 )
   W( 159 ) = JVS( 1234 )
   W( 160 ) = JVS( 1235 )
   W( 161 ) = JVS( 1236 )
   W( 162 ) = JVS( 1237 )
   W( 163 ) = JVS( 1238 )
   W( 164 ) = JVS( 1239 )
   W( 165 ) = JVS( 1240 )
   W( 166 ) = JVS( 1241 )
   W( 167 ) = JVS( 1242 )
   W( 168 ) = JVS( 1243 )
  a = -W( 104 ) / JVS( 525 )
  W( 104 ) = -a
  W( 163 ) = W( 163 ) + a*JVS( 526 )
  W( 167 ) = W( 167 ) + a*JVS( 527 )
  a = -W( 109 ) / JVS( 538 )
  W( 109 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 539 )
  W( 143 ) = W( 143 ) + a*JVS( 540 )
  W( 144 ) = W( 144 ) + a*JVS( 541 )
  W( 156 ) = W( 156 ) + a*JVS( 542 )
  W( 165 ) = W( 165 ) + a*JVS( 543 )
  a = -W( 115 ) / JVS( 558 )
  W( 115 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 559 )
  a = -W( 122 ) / JVS( 586 )
  W( 122 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 587 )
  a = -W( 123 ) / JVS( 590 )
  W( 123 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 591 )
  a = -W( 125 ) / JVS( 597 )
  W( 125 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 598 )
  a = -W( 132 ) / JVS( 636 )
  W( 132 ) = -a
  W( 142 ) = W( 142 ) + a*JVS( 637 )
  W( 146 ) = W( 146 ) + a*JVS( 638 )
  W( 148 ) = W( 148 ) + a*JVS( 639 )
  W( 156 ) = W( 156 ) + a*JVS( 640 )
  W( 157 ) = W( 157 ) + a*JVS( 641 )
  W( 165 ) = W( 165 ) + a*JVS( 642 )
  a = -W( 136 ) / JVS( 699 )
  W( 136 ) = -a
  W( 139 ) = W( 139 ) + a*JVS( 700 )
  W( 143 ) = W( 143 ) + a*JVS( 701 )
  W( 144 ) = W( 144 ) + a*JVS( 702 )
  W( 147 ) = W( 147 ) + a*JVS( 703 )
  W( 151 ) = W( 151 ) + a*JVS( 704 )
  W( 156 ) = W( 156 ) + a*JVS( 705 )
  W( 157 ) = W( 157 ) + a*JVS( 706 )
  W( 165 ) = W( 165 ) + a*JVS( 707 )
  W( 167 ) = W( 167 ) + a*JVS( 708 )
  a = -W( 139 ) / JVS( 726 )
  W( 139 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 727 )
  W( 156 ) = W( 156 ) + a*JVS( 728 )
  W( 157 ) = W( 157 ) + a*JVS( 729 )
  W( 165 ) = W( 165 ) + a*JVS( 730 )
  a = -W( 142 ) / JVS( 767 )
  W( 142 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 768 )
  W( 152 ) = W( 152 ) + a*JVS( 769 )
  W( 156 ) = W( 156 ) + a*JVS( 770 )
  W( 157 ) = W( 157 ) + a*JVS( 771 )
  W( 165 ) = W( 165 ) + a*JVS( 772 )
  a = -W( 143 ) / JVS( 773 )
  W( 143 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 774 )
  W( 156 ) = W( 156 ) + a*JVS( 775 )
  W( 157 ) = W( 157 ) + a*JVS( 776 )
  W( 165 ) = W( 165 ) + a*JVS( 777 )
  a = -W( 144 ) / JVS( 778 )
  W( 144 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 779 )
  W( 156 ) = W( 156 ) + a*JVS( 780 )
  W( 157 ) = W( 157 ) + a*JVS( 781 )
  W( 165 ) = W( 165 ) + a*JVS( 782 )
  a = -W( 146 ) / JVS( 808 )
  W( 146 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 809 )
  W( 152 ) = W( 152 ) + a*JVS( 810 )
  W( 156 ) = W( 156 ) + a*JVS( 811 )
  W( 157 ) = W( 157 ) + a*JVS( 812 )
  W( 165 ) = W( 165 ) + a*JVS( 813 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  a = -W( 148 ) / JVS( 821 )
  W( 148 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 822 )
  W( 156 ) = W( 156 ) + a*JVS( 823 )
  W( 157 ) = W( 157 ) + a*JVS( 824 )
  W( 165 ) = W( 165 ) + a*JVS( 825 )
  a = -W( 149 ) / JVS( 836 )
  W( 149 ) = -a
  W( 151 ) = W( 151 ) + a*JVS( 837 )
  W( 152 ) = W( 152 ) + a*JVS( 838 )
  W( 153 ) = W( 153 ) + a*JVS( 839 )
  W( 154 ) = W( 154 ) + a*JVS( 840 )
  W( 155 ) = W( 155 ) + a*JVS( 841 )
  W( 156 ) = W( 156 ) + a*JVS( 842 )
  W( 157 ) = W( 157 ) + a*JVS( 843 )
  W( 159 ) = W( 159 ) + a*JVS( 844 )
  W( 160 ) = W( 160 ) + a*JVS( 845 )
  W( 162 ) = W( 162 ) + a*JVS( 846 )
  W( 163 ) = W( 163 ) + a*JVS( 847 )
  W( 164 ) = W( 164 ) + a*JVS( 848 )
  W( 165 ) = W( 165 ) + a*JVS( 849 )
  a = -W( 151 ) / JVS( 894 )
  W( 151 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 895 )
  W( 156 ) = W( 156 ) + a*JVS( 896 )
  W( 157 ) = W( 157 ) + a*JVS( 897 )
  W( 158 ) = W( 158 ) + a*JVS( 898 )
  W( 160 ) = W( 160 ) + a*JVS( 899 )
  W( 165 ) = W( 165 ) + a*JVS( 900 )
  W( 167 ) = W( 167 ) + a*JVS( 901 )
  a = -W( 152 ) / JVS( 912 )
  W( 152 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 913 )
  W( 157 ) = W( 157 ) + a*JVS( 914 )
  W( 160 ) = W( 160 ) + a*JVS( 915 )
  W( 165 ) = W( 165 ) + a*JVS( 916 )
  W( 167 ) = W( 167 ) + a*JVS( 917 )
  a = -W( 153 ) / JVS( 937 )
  W( 153 ) = -a
  W( 154 ) = W( 154 ) + a*JVS( 938 )
  W( 155 ) = W( 155 ) + a*JVS( 939 )
  W( 156 ) = W( 156 ) + a*JVS( 940 )
  W( 157 ) = W( 157 ) + a*JVS( 941 )
  W( 158 ) = W( 158 ) + a*JVS( 942 )
  W( 160 ) = W( 160 ) + a*JVS( 943 )
  W( 161 ) = W( 161 ) + a*JVS( 944 )
  W( 165 ) = W( 165 ) + a*JVS( 945 )
  W( 166 ) = W( 166 ) + a*JVS( 946 )
  W( 167 ) = W( 167 ) + a*JVS( 947 )
  a = -W( 154 ) / JVS( 959 )
  W( 154 ) = -a
  W( 155 ) = W( 155 ) + a*JVS( 960 )
  W( 156 ) = W( 156 ) + a*JVS( 961 )
  W( 157 ) = W( 157 ) + a*JVS( 962 )
  W( 158 ) = W( 158 ) + a*JVS( 963 )
  W( 160 ) = W( 160 ) + a*JVS( 964 )
  W( 165 ) = W( 165 ) + a*JVS( 965 )
  W( 166 ) = W( 166 ) + a*JVS( 966 )
  W( 167 ) = W( 167 ) + a*JVS( 967 )
  W( 168 ) = W( 168 ) + a*JVS( 968 )
  a = -W( 155 ) / JVS( 981 )
  W( 155 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 982 )
  W( 157 ) = W( 157 ) + a*JVS( 983 )
  W( 158 ) = W( 158 ) + a*JVS( 984 )
  W( 160 ) = W( 160 ) + a*JVS( 985 )
  W( 162 ) = W( 162 ) + a*JVS( 986 )
  W( 163 ) = W( 163 ) + a*JVS( 987 )
  W( 164 ) = W( 164 ) + a*JVS( 988 )
  W( 165 ) = W( 165 ) + a*JVS( 989 )
  W( 166 ) = W( 166 ) + a*JVS( 990 )
  W( 167 ) = W( 167 ) + a*JVS( 991 )
  W( 168 ) = W( 168 ) + a*JVS( 992 )
  a = -W( 156 ) / JVS( 1006 )
  W( 156 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 1007 )
  W( 159 ) = W( 159 ) + a*JVS( 1008 )
  W( 160 ) = W( 160 ) + a*JVS( 1009 )
  W( 161 ) = W( 161 ) + a*JVS( 1010 )
  W( 162 ) = W( 162 ) + a*JVS( 1011 )
  W( 163 ) = W( 163 ) + a*JVS( 1012 )
  W( 164 ) = W( 164 ) + a*JVS( 1013 )
  W( 165 ) = W( 165 ) + a*JVS( 1014 )
  W( 167 ) = W( 167 ) + a*JVS( 1015 )
  a = -W( 157 ) / JVS( 1045 )
  W( 157 ) = -a
  W( 158 ) = W( 158 ) + a*JVS( 1046 )
  W( 159 ) = W( 159 ) + a*JVS( 1047 )
  W( 160 ) = W( 160 ) + a*JVS( 1048 )
  W( 161 ) = W( 161 ) + a*JVS( 1049 )
  W( 162 ) = W( 162 ) + a*JVS( 1050 )
  W( 163 ) = W( 163 ) + a*JVS( 1051 )
  W( 164 ) = W( 164 ) + a*JVS( 1052 )
  W( 165 ) = W( 165 ) + a*JVS( 1053 )
  W( 166 ) = W( 166 ) + a*JVS( 1054 )
  W( 167 ) = W( 167 ) + a*JVS( 1055 )
  W( 168 ) = W( 168 ) + a*JVS( 1056 )
  a = -W( 158 ) / JVS( 1078 )
  W( 158 ) = -a
  W( 159 ) = W( 159 ) + a*JVS( 1079 )
  W( 160 ) = W( 160 ) + a*JVS( 1080 )
  W( 161 ) = W( 161 ) + a*JVS( 1081 )
  W( 162 ) = W( 162 ) + a*JVS( 1082 )
  W( 163 ) = W( 163 ) + a*JVS( 1083 )
  W( 164 ) = W( 164 ) + a*JVS( 1084 )
  W( 165 ) = W( 165 ) + a*JVS( 1085 )
  W( 166 ) = W( 166 ) + a*JVS( 1086 )
  W( 167 ) = W( 167 ) + a*JVS( 1087 )
  W( 168 ) = W( 168 ) + a*JVS( 1088 )
  a = -W( 159 ) / JVS( 1099 )
  W( 159 ) = -a
  W( 160 ) = W( 160 ) + a*JVS( 1100 )
  W( 161 ) = W( 161 ) + a*JVS( 1101 )
  W( 162 ) = W( 162 ) + a*JVS( 1102 )
  W( 163 ) = W( 163 ) + a*JVS( 1103 )
  W( 164 ) = W( 164 ) + a*JVS( 1104 )
  W( 165 ) = W( 165 ) + a*JVS( 1105 )
  W( 166 ) = W( 166 ) + a*JVS( 1106 )
  W( 167 ) = W( 167 ) + a*JVS( 1107 )
  W( 168 ) = W( 168 ) + a*JVS( 1108 )
  a = -W( 160 ) / JVS( 1125 )
  W( 160 ) = -a
  W( 161 ) = W( 161 ) + a*JVS( 1126 )
  W( 162 ) = W( 162 ) + a*JVS( 1127 )
  W( 163 ) = W( 163 ) + a*JVS( 1128 )
  W( 164 ) = W( 164 ) + a*JVS( 1129 )
  W( 165 ) = W( 165 ) + a*JVS( 1130 )
  W( 166 ) = W( 166 ) + a*JVS( 1131 )
  W( 167 ) = W( 167 ) + a*JVS( 1132 )
  W( 168 ) = W( 168 ) + a*JVS( 1133 )
  a = -W( 161 ) / JVS( 1177 )
  W( 161 ) = -a
  W( 162 ) = W( 162 ) + a*JVS( 1178 )
  W( 163 ) = W( 163 ) + a*JVS( 1179 )
  W( 164 ) = W( 164 ) + a*JVS( 1180 )
  W( 165 ) = W( 165 ) + a*JVS( 1181 )
  W( 166 ) = W( 166 ) + a*JVS( 1182 )
  W( 167 ) = W( 167 ) + a*JVS( 1183 )
  W( 168 ) = W( 168 ) + a*JVS( 1184 )
  a = -W( 162 ) / JVS( 1203 )
  W( 162 ) = -a
  W( 163 ) = W( 163 ) + a*JVS( 1204 )
  W( 164 ) = W( 164 ) + a*JVS( 1205 )
  W( 165 ) = W( 165 ) + a*JVS( 1206 )
  W( 166 ) = W( 166 ) + a*JVS( 1207 )
  W( 167 ) = W( 167 ) + a*JVS( 1208 )
  W( 168 ) = W( 168 ) + a*JVS( 1209 )
  JVS( 1210) = W( 104 )
  JVS( 1211) = W( 109 )
  JVS( 1212) = W( 115 )
  JVS( 1213) = W( 122 )
  JVS( 1214) = W( 123 )
  JVS( 1215) = W( 125 )
  JVS( 1216) = W( 132 )
  JVS( 1217) = W( 136 )
  JVS( 1218) = W( 139 )
  JVS( 1219) = W( 142 )
  JVS( 1220) = W( 143 )
  JVS( 1221) = W( 144 )
  JVS( 1222) = W( 146 )
  JVS( 1223) = W( 147 )
  JVS( 1224) = W( 148 )
  JVS( 1225) = W( 149 )
  JVS( 1226) = W( 151 )
  JVS( 1227) = W( 152 )
  JVS( 1228) = W( 153 )
  JVS( 1229) = W( 154 )
  JVS( 1230) = W( 155 )
  JVS( 1231) = W( 156 )
  JVS( 1232) = W( 157 )
  JVS( 1233) = W( 158 )
  JVS( 1234) = W( 159 )
  JVS( 1235) = W( 160 )
  JVS( 1236) = W( 161 )
  JVS( 1237) = W( 162 )
  JVS( 1238) = W( 163 )
  JVS( 1239) = W( 164 )
  JVS( 1240) = W( 165 )
  JVS( 1241) = W( 166 )
  JVS( 1242) = W( 167 )
  JVS( 1243) = W( 168 )
  IF ( ABS( JVS( 1256 )) < TINY(a) ) THEN
         IER = 164
         RETURN
  END IF
   W( 105 ) = JVS( 1244 )
   W( 130 ) = JVS( 1245 )
   W( 147 ) = JVS( 1246 )
   W( 152 ) = JVS( 1247 )
   W( 156 ) = JVS( 1248 )
   W( 157 ) = JVS( 1249 )
   W( 158 ) = JVS( 1250 )
   W( 159 ) = JVS( 1251 )
   W( 160 ) = JVS( 1252 )
   W( 161 ) = JVS( 1253 )
   W( 162 ) = JVS( 1254 )
   W( 163 ) = JVS( 1255 )
   W( 164 ) = JVS( 1256 )
   W( 165 ) = JVS( 1257 )
   W( 166 ) = JVS( 1258 )
   W( 167 ) = JVS( 1259 )
   W( 168 ) = JVS( 1260 )
  a = -W( 105 ) / JVS( 528 )
  W( 105 ) = -a
  W( 164 ) = W( 164 ) + a*JVS( 529 )
  W( 167 ) = W( 167 ) + a*JVS( 530 )
  a = -W( 130 ) / JVS( 620 )
  W( 130 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 621 )
  W( 156 ) = W( 156 ) + a*JVS( 622 )
  W( 157 ) = W( 157 ) + a*JVS( 623 )
  W( 165 ) = W( 165 ) + a*JVS( 624 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  a = -W( 152 ) / JVS( 912 )
  W( 152 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 913 )
  W( 157 ) = W( 157 ) + a*JVS( 914 )
  W( 160 ) = W( 160 ) + a*JVS( 915 )
  W( 165 ) = W( 165 ) + a*JVS( 916 )
  W( 167 ) = W( 167 ) + a*JVS( 917 )
  a = -W( 156 ) / JVS( 1006 )
  W( 156 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 1007 )
  W( 159 ) = W( 159 ) + a*JVS( 1008 )
  W( 160 ) = W( 160 ) + a*JVS( 1009 )
  W( 161 ) = W( 161 ) + a*JVS( 1010 )
  W( 162 ) = W( 162 ) + a*JVS( 1011 )
  W( 163 ) = W( 163 ) + a*JVS( 1012 )
  W( 164 ) = W( 164 ) + a*JVS( 1013 )
  W( 165 ) = W( 165 ) + a*JVS( 1014 )
  W( 167 ) = W( 167 ) + a*JVS( 1015 )
  a = -W( 157 ) / JVS( 1045 )
  W( 157 ) = -a
  W( 158 ) = W( 158 ) + a*JVS( 1046 )
  W( 159 ) = W( 159 ) + a*JVS( 1047 )
  W( 160 ) = W( 160 ) + a*JVS( 1048 )
  W( 161 ) = W( 161 ) + a*JVS( 1049 )
  W( 162 ) = W( 162 ) + a*JVS( 1050 )
  W( 163 ) = W( 163 ) + a*JVS( 1051 )
  W( 164 ) = W( 164 ) + a*JVS( 1052 )
  W( 165 ) = W( 165 ) + a*JVS( 1053 )
  W( 166 ) = W( 166 ) + a*JVS( 1054 )
  W( 167 ) = W( 167 ) + a*JVS( 1055 )
  W( 168 ) = W( 168 ) + a*JVS( 1056 )
  a = -W( 158 ) / JVS( 1078 )
  W( 158 ) = -a
  W( 159 ) = W( 159 ) + a*JVS( 1079 )
  W( 160 ) = W( 160 ) + a*JVS( 1080 )
  W( 161 ) = W( 161 ) + a*JVS( 1081 )
  W( 162 ) = W( 162 ) + a*JVS( 1082 )
  W( 163 ) = W( 163 ) + a*JVS( 1083 )
  W( 164 ) = W( 164 ) + a*JVS( 1084 )
  W( 165 ) = W( 165 ) + a*JVS( 1085 )
  W( 166 ) = W( 166 ) + a*JVS( 1086 )
  W( 167 ) = W( 167 ) + a*JVS( 1087 )
  W( 168 ) = W( 168 ) + a*JVS( 1088 )
  a = -W( 159 ) / JVS( 1099 )
  W( 159 ) = -a
  W( 160 ) = W( 160 ) + a*JVS( 1100 )
  W( 161 ) = W( 161 ) + a*JVS( 1101 )
  W( 162 ) = W( 162 ) + a*JVS( 1102 )
  W( 163 ) = W( 163 ) + a*JVS( 1103 )
  W( 164 ) = W( 164 ) + a*JVS( 1104 )
  W( 165 ) = W( 165 ) + a*JVS( 1105 )
  W( 166 ) = W( 166 ) + a*JVS( 1106 )
  W( 167 ) = W( 167 ) + a*JVS( 1107 )
  W( 168 ) = W( 168 ) + a*JVS( 1108 )
  a = -W( 160 ) / JVS( 1125 )
  W( 160 ) = -a
  W( 161 ) = W( 161 ) + a*JVS( 1126 )
  W( 162 ) = W( 162 ) + a*JVS( 1127 )
  W( 163 ) = W( 163 ) + a*JVS( 1128 )
  W( 164 ) = W( 164 ) + a*JVS( 1129 )
  W( 165 ) = W( 165 ) + a*JVS( 1130 )
  W( 166 ) = W( 166 ) + a*JVS( 1131 )
  W( 167 ) = W( 167 ) + a*JVS( 1132 )
  W( 168 ) = W( 168 ) + a*JVS( 1133 )
  a = -W( 161 ) / JVS( 1177 )
  W( 161 ) = -a
  W( 162 ) = W( 162 ) + a*JVS( 1178 )
  W( 163 ) = W( 163 ) + a*JVS( 1179 )
  W( 164 ) = W( 164 ) + a*JVS( 1180 )
  W( 165 ) = W( 165 ) + a*JVS( 1181 )
  W( 166 ) = W( 166 ) + a*JVS( 1182 )
  W( 167 ) = W( 167 ) + a*JVS( 1183 )
  W( 168 ) = W( 168 ) + a*JVS( 1184 )
  a = -W( 162 ) / JVS( 1203 )
  W( 162 ) = -a
  W( 163 ) = W( 163 ) + a*JVS( 1204 )
  W( 164 ) = W( 164 ) + a*JVS( 1205 )
  W( 165 ) = W( 165 ) + a*JVS( 1206 )
  W( 166 ) = W( 166 ) + a*JVS( 1207 )
  W( 167 ) = W( 167 ) + a*JVS( 1208 )
  W( 168 ) = W( 168 ) + a*JVS( 1209 )
  a = -W( 163 ) / JVS( 1238 )
  W( 163 ) = -a
  W( 164 ) = W( 164 ) + a*JVS( 1239 )
  W( 165 ) = W( 165 ) + a*JVS( 1240 )
  W( 166 ) = W( 166 ) + a*JVS( 1241 )
  W( 167 ) = W( 167 ) + a*JVS( 1242 )
  W( 168 ) = W( 168 ) + a*JVS( 1243 )
  JVS( 1244) = W( 105 )
  JVS( 1245) = W( 130 )
  JVS( 1246) = W( 147 )
  JVS( 1247) = W( 152 )
  JVS( 1248) = W( 156 )
  JVS( 1249) = W( 157 )
  JVS( 1250) = W( 158 )
  JVS( 1251) = W( 159 )
  JVS( 1252) = W( 160 )
  JVS( 1253) = W( 161 )
  JVS( 1254) = W( 162 )
  JVS( 1255) = W( 163 )
  JVS( 1256) = W( 164 )
  JVS( 1257) = W( 165 )
  JVS( 1258) = W( 166 )
  JVS( 1259) = W( 167 )
  JVS( 1260) = W( 168 )
  IF ( ABS( JVS( 1351 )) < TINY(a) ) THEN
         IER = 165
         RETURN
  END IF
   W( 3 ) = JVS( 1261 )
   W( 5 ) = JVS( 1262 )
   W( 66 ) = JVS( 1263 )
   W( 67 ) = JVS( 1264 )
   W( 68 ) = JVS( 1265 )
   W( 69 ) = JVS( 1266 )
   W( 70 ) = JVS( 1267 )
   W( 71 ) = JVS( 1268 )
   W( 72 ) = JVS( 1269 )
   W( 73 ) = JVS( 1270 )
   W( 74 ) = JVS( 1271 )
   W( 75 ) = JVS( 1272 )
   W( 76 ) = JVS( 1273 )
   W( 77 ) = JVS( 1274 )
   W( 78 ) = JVS( 1275 )
   W( 79 ) = JVS( 1276 )
   W( 80 ) = JVS( 1277 )
   W( 81 ) = JVS( 1278 )
   W( 82 ) = JVS( 1279 )
   W( 83 ) = JVS( 1280 )
   W( 84 ) = JVS( 1281 )
   W( 85 ) = JVS( 1282 )
   W( 86 ) = JVS( 1283 )
   W( 87 ) = JVS( 1284 )
   W( 88 ) = JVS( 1285 )
   W( 89 ) = JVS( 1286 )
   W( 90 ) = JVS( 1287 )
   W( 91 ) = JVS( 1288 )
   W( 92 ) = JVS( 1289 )
   W( 93 ) = JVS( 1290 )
   W( 94 ) = JVS( 1291 )
   W( 95 ) = JVS( 1292 )
   W( 96 ) = JVS( 1293 )
   W( 97 ) = JVS( 1294 )
   W( 98 ) = JVS( 1295 )
   W( 99 ) = JVS( 1296 )
   W( 100 ) = JVS( 1297 )
   W( 101 ) = JVS( 1298 )
   W( 106 ) = JVS( 1299 )
   W( 107 ) = JVS( 1300 )
   W( 108 ) = JVS( 1301 )
   W( 111 ) = JVS( 1302 )
   W( 112 ) = JVS( 1303 )
   W( 113 ) = JVS( 1304 )
   W( 115 ) = JVS( 1305 )
   W( 116 ) = JVS( 1306 )
   W( 119 ) = JVS( 1307 )
   W( 120 ) = JVS( 1308 )
   W( 121 ) = JVS( 1309 )
   W( 122 ) = JVS( 1310 )
   W( 123 ) = JVS( 1311 )
   W( 124 ) = JVS( 1312 )
   W( 125 ) = JVS( 1313 )
   W( 126 ) = JVS( 1314 )
   W( 127 ) = JVS( 1315 )
   W( 129 ) = JVS( 1316 )
   W( 130 ) = JVS( 1317 )
   W( 131 ) = JVS( 1318 )
   W( 132 ) = JVS( 1319 )
   W( 133 ) = JVS( 1320 )
   W( 134 ) = JVS( 1321 )
   W( 135 ) = JVS( 1322 )
   W( 136 ) = JVS( 1323 )
   W( 137 ) = JVS( 1324 )
   W( 138 ) = JVS( 1325 )
   W( 139 ) = JVS( 1326 )
   W( 140 ) = JVS( 1327 )
   W( 142 ) = JVS( 1328 )
   W( 143 ) = JVS( 1329 )
   W( 144 ) = JVS( 1330 )
   W( 145 ) = JVS( 1331 )
   W( 146 ) = JVS( 1332 )
   W( 147 ) = JVS( 1333 )
   W( 148 ) = JVS( 1334 )
   W( 149 ) = JVS( 1335 )
   W( 150 ) = JVS( 1336 )
   W( 151 ) = JVS( 1337 )
   W( 152 ) = JVS( 1338 )
   W( 153 ) = JVS( 1339 )
   W( 154 ) = JVS( 1340 )
   W( 155 ) = JVS( 1341 )
   W( 156 ) = JVS( 1342 )
   W( 157 ) = JVS( 1343 )
   W( 158 ) = JVS( 1344 )
   W( 159 ) = JVS( 1345 )
   W( 160 ) = JVS( 1346 )
   W( 161 ) = JVS( 1347 )
   W( 162 ) = JVS( 1348 )
   W( 163 ) = JVS( 1349 )
   W( 164 ) = JVS( 1350 )
   W( 165 ) = JVS( 1351 )
   W( 166 ) = JVS( 1352 )
   W( 167 ) = JVS( 1353 )
   W( 168 ) = JVS( 1354 )
  a = -W( 3 ) / JVS( 8 )
  W( 3 ) = -a
  a = -W( 5 ) / JVS( 10 )
  W( 5 ) = -a
  a = -W( 66 ) / JVS( 417 )
  W( 66 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 418 )
  a = -W( 67 ) / JVS( 419 )
  W( 67 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 420 )
  a = -W( 68 ) / JVS( 421 )
  W( 68 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 422 )
  a = -W( 69 ) / JVS( 423 )
  W( 69 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 424 )
  a = -W( 70 ) / JVS( 425 )
  W( 70 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 426 )
  a = -W( 71 ) / JVS( 427 )
  W( 71 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 428 )
  a = -W( 72 ) / JVS( 429 )
  W( 72 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 430 )
  a = -W( 73 ) / JVS( 431 )
  W( 73 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 432 )
  a = -W( 74 ) / JVS( 434 )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 435 )
  W( 165 ) = W( 165 ) + a*JVS( 436 )
  a = -W( 75 ) / JVS( 438 )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 439 )
  W( 165 ) = W( 165 ) + a*JVS( 440 )
  a = -W( 76 ) / JVS( 442 )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 443 )
  W( 165 ) = W( 165 ) + a*JVS( 444 )
  a = -W( 77 ) / JVS( 446 )
  W( 77 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 447 )
  W( 165 ) = W( 165 ) + a*JVS( 448 )
  a = -W( 78 ) / JVS( 450 )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 451 )
  W( 165 ) = W( 165 ) + a*JVS( 452 )
  a = -W( 79 ) / JVS( 454 )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 455 )
  W( 165 ) = W( 165 ) + a*JVS( 456 )
  a = -W( 80 ) / JVS( 458 )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 459 )
  W( 165 ) = W( 165 ) + a*JVS( 460 )
  a = -W( 81 ) / JVS( 462 )
  W( 81 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 463 )
  a = -W( 82 ) / JVS( 464 )
  W( 82 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 465 )
  a = -W( 83 ) / JVS( 466 )
  W( 83 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 467 )
  a = -W( 84 ) / JVS( 468 )
  W( 84 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 469 )
  a = -W( 85 ) / JVS( 470 )
  W( 85 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 471 )
  a = -W( 86 ) / JVS( 472 )
  W( 86 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 473 )
  a = -W( 87 ) / JVS( 474 )
  W( 87 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 475 )
  a = -W( 88 ) / JVS( 476 )
  W( 88 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 477 )
  a = -W( 89 ) / JVS( 478 )
  W( 89 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 479 )
  a = -W( 90 ) / JVS( 480 )
  W( 90 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 481 )
  a = -W( 91 ) / JVS( 483 )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 484 )
  W( 165 ) = W( 165 ) + a*JVS( 485 )
  a = -W( 92 ) / JVS( 487 )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 488 )
  W( 165 ) = W( 165 ) + a*JVS( 489 )
  a = -W( 93 ) / JVS( 491 )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 492 )
  W( 165 ) = W( 165 ) + a*JVS( 493 )
  a = -W( 94 ) / JVS( 495 )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 496 )
  W( 165 ) = W( 165 ) + a*JVS( 497 )
  a = -W( 95 ) / JVS( 499 )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 500 )
  W( 165 ) = W( 165 ) + a*JVS( 501 )
  a = -W( 96 ) / JVS( 503 )
  W( 96 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 504 )
  W( 165 ) = W( 165 ) + a*JVS( 505 )
  a = -W( 97 ) / JVS( 507 )
  W( 97 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 508 )
  W( 165 ) = W( 165 ) + a*JVS( 509 )
  a = -W( 98 ) / JVS( 511 )
  W( 98 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 512 )
  a = -W( 99 ) / JVS( 513 )
  W( 99 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 514 )
  a = -W( 100 ) / JVS( 515 )
  W( 100 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 516 )
  a = -W( 101 ) / JVS( 517 )
  W( 101 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 518 )
  a = -W( 106 ) / JVS( 531 )
  W( 106 ) = -a
  W( 161 ) = W( 161 ) + a*JVS( 532 )
  W( 165 ) = W( 165 ) + a*JVS( 533 )
  a = -W( 107 ) / JVS( 534 )
  W( 107 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 535 )
  a = -W( 108 ) / JVS( 536 )
  W( 108 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 537 )
  a = -W( 111 ) / JVS( 547 )
  W( 111 ) = -a
  W( 160 ) = W( 160 ) + a*JVS( 548 )
  W( 165 ) = W( 165 ) + a*JVS( 549 )
  a = -W( 112 ) / JVS( 550 )
  W( 112 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 551 )
  a = -W( 113 ) / JVS( 552 )
  W( 113 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 553 )
  a = -W( 115 ) / JVS( 558 )
  W( 115 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 559 )
  a = -W( 116 ) / JVS( 560 )
  W( 116 ) = -a
  W( 161 ) = W( 161 ) + a*JVS( 561 )
  W( 165 ) = W( 165 ) + a*JVS( 562 )
  W( 168 ) = W( 168 ) + a*JVS( 563 )
  a = -W( 119 ) / JVS( 573 )
  W( 119 ) = -a
  W( 161 ) = W( 161 ) + a*JVS( 574 )
  W( 165 ) = W( 165 ) + a*JVS( 575 )
  W( 167 ) = W( 167 ) + a*JVS( 576 )
  a = -W( 120 ) / JVS( 577 )
  W( 120 ) = -a
  W( 158 ) = W( 158 ) + a*JVS( 578 )
  W( 165 ) = W( 165 ) + a*JVS( 579 )
  W( 166 ) = W( 166 ) + a*JVS( 580 )
  W( 168 ) = W( 168 ) + a*JVS( 581 )
  a = -W( 121 ) / JVS( 582 )
  W( 121 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 583 )
  a = -W( 122 ) / JVS( 586 )
  W( 122 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 587 )
  a = -W( 123 ) / JVS( 590 )
  W( 123 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 591 )
  a = -W( 124 ) / JVS( 594 )
  W( 124 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 595 )
  W( 165 ) = W( 165 ) + a*JVS( 596 )
  a = -W( 125 ) / JVS( 597 )
  W( 125 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 598 )
  a = -W( 126 ) / JVS( 599 )
  W( 126 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 600 )
  W( 165 ) = W( 165 ) + a*JVS( 601 )
  a = -W( 127 ) / JVS( 604 )
  W( 127 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 605 )
  W( 165 ) = W( 165 ) + a*JVS( 606 )
  a = -W( 129 ) / JVS( 613 )
  W( 129 ) = -a
  W( 158 ) = W( 158 ) + a*JVS( 614 )
  W( 161 ) = W( 161 ) + a*JVS( 615 )
  W( 165 ) = W( 165 ) + a*JVS( 616 )
  W( 166 ) = W( 166 ) + a*JVS( 617 )
  a = -W( 130 ) / JVS( 620 )
  W( 130 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 621 )
  W( 156 ) = W( 156 ) + a*JVS( 622 )
  W( 157 ) = W( 157 ) + a*JVS( 623 )
  W( 165 ) = W( 165 ) + a*JVS( 624 )
  a = -W( 131 ) / JVS( 626 )
  W( 131 ) = -a
  W( 138 ) = W( 138 ) + a*JVS( 627 )
  W( 157 ) = W( 157 ) + a*JVS( 628 )
  W( 161 ) = W( 161 ) + a*JVS( 629 )
  W( 165 ) = W( 165 ) + a*JVS( 630 )
  a = -W( 132 ) / JVS( 636 )
  W( 132 ) = -a
  W( 142 ) = W( 142 ) + a*JVS( 637 )
  W( 146 ) = W( 146 ) + a*JVS( 638 )
  W( 148 ) = W( 148 ) + a*JVS( 639 )
  W( 156 ) = W( 156 ) + a*JVS( 640 )
  W( 157 ) = W( 157 ) + a*JVS( 641 )
  W( 165 ) = W( 165 ) + a*JVS( 642 )
  a = -W( 133 ) / JVS( 649 )
  W( 133 ) = -a
  W( 135 ) = W( 135 ) + a*JVS( 650 )
  W( 137 ) = W( 137 ) + a*JVS( 651 )
  W( 139 ) = W( 139 ) + a*JVS( 652 )
  W( 140 ) = W( 140 ) + a*JVS( 653 )
  W( 142 ) = W( 142 ) + a*JVS( 654 )
  W( 143 ) = W( 143 ) + a*JVS( 655 )
  W( 144 ) = W( 144 ) + a*JVS( 656 )
  W( 145 ) = W( 145 ) + a*JVS( 657 )
  W( 146 ) = W( 146 ) + a*JVS( 658 )
  W( 147 ) = W( 147 ) + a*JVS( 659 )
  W( 148 ) = W( 148 ) + a*JVS( 660 )
  W( 149 ) = W( 149 ) + a*JVS( 661 )
  W( 150 ) = W( 150 ) + a*JVS( 662 )
  W( 152 ) = W( 152 ) + a*JVS( 663 )
  W( 153 ) = W( 153 ) + a*JVS( 664 )
  W( 156 ) = W( 156 ) + a*JVS( 665 )
  W( 157 ) = W( 157 ) + a*JVS( 666 )
  W( 165 ) = W( 165 ) + a*JVS( 667 )
  a = -W( 134 ) / JVS( 674 )
  W( 134 ) = -a
  W( 138 ) = W( 138 ) + a*JVS( 675 )
  W( 142 ) = W( 142 ) + a*JVS( 676 )
  W( 145 ) = W( 145 ) + a*JVS( 677 )
  W( 146 ) = W( 146 ) + a*JVS( 678 )
  W( 147 ) = W( 147 ) + a*JVS( 679 )
  W( 148 ) = W( 148 ) + a*JVS( 680 )
  W( 149 ) = W( 149 ) + a*JVS( 681 )
  W( 150 ) = W( 150 ) + a*JVS( 682 )
  W( 153 ) = W( 153 ) + a*JVS( 683 )
  W( 156 ) = W( 156 ) + a*JVS( 684 )
  W( 157 ) = W( 157 ) + a*JVS( 685 )
  W( 161 ) = W( 161 ) + a*JVS( 686 )
  W( 165 ) = W( 165 ) + a*JVS( 687 )
  W( 167 ) = W( 167 ) + a*JVS( 688 )
  a = -W( 135 ) / JVS( 689 )
  W( 135 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 690 )
  W( 156 ) = W( 156 ) + a*JVS( 691 )
  W( 157 ) = W( 157 ) + a*JVS( 692 )
  W( 165 ) = W( 165 ) + a*JVS( 693 )
  a = -W( 136 ) / JVS( 699 )
  W( 136 ) = -a
  W( 139 ) = W( 139 ) + a*JVS( 700 )
  W( 143 ) = W( 143 ) + a*JVS( 701 )
  W( 144 ) = W( 144 ) + a*JVS( 702 )
  W( 147 ) = W( 147 ) + a*JVS( 703 )
  W( 151 ) = W( 151 ) + a*JVS( 704 )
  W( 156 ) = W( 156 ) + a*JVS( 705 )
  W( 157 ) = W( 157 ) + a*JVS( 706 )
  W( 165 ) = W( 165 ) + a*JVS( 707 )
  W( 167 ) = W( 167 ) + a*JVS( 708 )
  a = -W( 137 ) / JVS( 709 )
  W( 137 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 710 )
  W( 156 ) = W( 156 ) + a*JVS( 711 )
  W( 157 ) = W( 157 ) + a*JVS( 712 )
  W( 165 ) = W( 165 ) + a*JVS( 713 )
  a = -W( 138 ) / JVS( 716 )
  W( 138 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 717 )
  W( 159 ) = W( 159 ) + a*JVS( 718 )
  W( 160 ) = W( 160 ) + a*JVS( 719 )
  W( 161 ) = W( 161 ) + a*JVS( 720 )
  W( 162 ) = W( 162 ) + a*JVS( 721 )
  W( 163 ) = W( 163 ) + a*JVS( 722 )
  W( 164 ) = W( 164 ) + a*JVS( 723 )
  W( 165 ) = W( 165 ) + a*JVS( 724 )
  W( 167 ) = W( 167 ) + a*JVS( 725 )
  a = -W( 139 ) / JVS( 726 )
  W( 139 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 727 )
  W( 156 ) = W( 156 ) + a*JVS( 728 )
  W( 157 ) = W( 157 ) + a*JVS( 729 )
  W( 165 ) = W( 165 ) + a*JVS( 730 )
  a = -W( 140 ) / JVS( 731 )
  W( 140 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 732 )
  W( 156 ) = W( 156 ) + a*JVS( 733 )
  W( 157 ) = W( 157 ) + a*JVS( 734 )
  W( 165 ) = W( 165 ) + a*JVS( 735 )
  a = -W( 142 ) / JVS( 767 )
  W( 142 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 768 )
  W( 152 ) = W( 152 ) + a*JVS( 769 )
  W( 156 ) = W( 156 ) + a*JVS( 770 )
  W( 157 ) = W( 157 ) + a*JVS( 771 )
  W( 165 ) = W( 165 ) + a*JVS( 772 )
  a = -W( 143 ) / JVS( 773 )
  W( 143 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 774 )
  W( 156 ) = W( 156 ) + a*JVS( 775 )
  W( 157 ) = W( 157 ) + a*JVS( 776 )
  W( 165 ) = W( 165 ) + a*JVS( 777 )
  a = -W( 144 ) / JVS( 778 )
  W( 144 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 779 )
  W( 156 ) = W( 156 ) + a*JVS( 780 )
  W( 157 ) = W( 157 ) + a*JVS( 781 )
  W( 165 ) = W( 165 ) + a*JVS( 782 )
  a = -W( 145 ) / JVS( 794 )
  W( 145 ) = -a
  W( 146 ) = W( 146 ) + a*JVS( 795 )
  W( 152 ) = W( 152 ) + a*JVS( 796 )
  W( 156 ) = W( 156 ) + a*JVS( 797 )
  W( 157 ) = W( 157 ) + a*JVS( 798 )
  W( 159 ) = W( 159 ) + a*JVS( 799 )
  W( 160 ) = W( 160 ) + a*JVS( 800 )
  W( 161 ) = W( 161 ) + a*JVS( 801 )
  W( 162 ) = W( 162 ) + a*JVS( 802 )
  W( 163 ) = W( 163 ) + a*JVS( 803 )
  W( 164 ) = W( 164 ) + a*JVS( 804 )
  W( 165 ) = W( 165 ) + a*JVS( 805 )
  W( 167 ) = W( 167 ) + a*JVS( 806 )
  a = -W( 146 ) / JVS( 808 )
  W( 146 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 809 )
  W( 152 ) = W( 152 ) + a*JVS( 810 )
  W( 156 ) = W( 156 ) + a*JVS( 811 )
  W( 157 ) = W( 157 ) + a*JVS( 812 )
  W( 165 ) = W( 165 ) + a*JVS( 813 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  a = -W( 148 ) / JVS( 821 )
  W( 148 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 822 )
  W( 156 ) = W( 156 ) + a*JVS( 823 )
  W( 157 ) = W( 157 ) + a*JVS( 824 )
  W( 165 ) = W( 165 ) + a*JVS( 825 )
  a = -W( 149 ) / JVS( 836 )
  W( 149 ) = -a
  W( 151 ) = W( 151 ) + a*JVS( 837 )
  W( 152 ) = W( 152 ) + a*JVS( 838 )
  W( 153 ) = W( 153 ) + a*JVS( 839 )
  W( 154 ) = W( 154 ) + a*JVS( 840 )
  W( 155 ) = W( 155 ) + a*JVS( 841 )
  W( 156 ) = W( 156 ) + a*JVS( 842 )
  W( 157 ) = W( 157 ) + a*JVS( 843 )
  W( 159 ) = W( 159 ) + a*JVS( 844 )
  W( 160 ) = W( 160 ) + a*JVS( 845 )
  W( 162 ) = W( 162 ) + a*JVS( 846 )
  W( 163 ) = W( 163 ) + a*JVS( 847 )
  W( 164 ) = W( 164 ) + a*JVS( 848 )
  W( 165 ) = W( 165 ) + a*JVS( 849 )
  a = -W( 150 ) / JVS( 870 )
  W( 150 ) = -a
  W( 151 ) = W( 151 ) + a*JVS( 871 )
  W( 152 ) = W( 152 ) + a*JVS( 872 )
  W( 154 ) = W( 154 ) + a*JVS( 873 )
  W( 155 ) = W( 155 ) + a*JVS( 874 )
  W( 156 ) = W( 156 ) + a*JVS( 875 )
  W( 157 ) = W( 157 ) + a*JVS( 876 )
  W( 158 ) = W( 158 ) + a*JVS( 877 )
  W( 159 ) = W( 159 ) + a*JVS( 878 )
  W( 160 ) = W( 160 ) + a*JVS( 879 )
  W( 161 ) = W( 161 ) + a*JVS( 880 )
  W( 162 ) = W( 162 ) + a*JVS( 881 )
  W( 163 ) = W( 163 ) + a*JVS( 882 )
  W( 164 ) = W( 164 ) + a*JVS( 883 )
  W( 165 ) = W( 165 ) + a*JVS( 884 )
  W( 166 ) = W( 166 ) + a*JVS( 885 )
  W( 167 ) = W( 167 ) + a*JVS( 886 )
  W( 168 ) = W( 168 ) + a*JVS( 887 )
  a = -W( 151 ) / JVS( 894 )
  W( 151 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 895 )
  W( 156 ) = W( 156 ) + a*JVS( 896 )
  W( 157 ) = W( 157 ) + a*JVS( 897 )
  W( 158 ) = W( 158 ) + a*JVS( 898 )
  W( 160 ) = W( 160 ) + a*JVS( 899 )
  W( 165 ) = W( 165 ) + a*JVS( 900 )
  W( 167 ) = W( 167 ) + a*JVS( 901 )
  a = -W( 152 ) / JVS( 912 )
  W( 152 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 913 )
  W( 157 ) = W( 157 ) + a*JVS( 914 )
  W( 160 ) = W( 160 ) + a*JVS( 915 )
  W( 165 ) = W( 165 ) + a*JVS( 916 )
  W( 167 ) = W( 167 ) + a*JVS( 917 )
  a = -W( 153 ) / JVS( 937 )
  W( 153 ) = -a
  W( 154 ) = W( 154 ) + a*JVS( 938 )
  W( 155 ) = W( 155 ) + a*JVS( 939 )
  W( 156 ) = W( 156 ) + a*JVS( 940 )
  W( 157 ) = W( 157 ) + a*JVS( 941 )
  W( 158 ) = W( 158 ) + a*JVS( 942 )
  W( 160 ) = W( 160 ) + a*JVS( 943 )
  W( 161 ) = W( 161 ) + a*JVS( 944 )
  W( 165 ) = W( 165 ) + a*JVS( 945 )
  W( 166 ) = W( 166 ) + a*JVS( 946 )
  W( 167 ) = W( 167 ) + a*JVS( 947 )
  a = -W( 154 ) / JVS( 959 )
  W( 154 ) = -a
  W( 155 ) = W( 155 ) + a*JVS( 960 )
  W( 156 ) = W( 156 ) + a*JVS( 961 )
  W( 157 ) = W( 157 ) + a*JVS( 962 )
  W( 158 ) = W( 158 ) + a*JVS( 963 )
  W( 160 ) = W( 160 ) + a*JVS( 964 )
  W( 165 ) = W( 165 ) + a*JVS( 965 )
  W( 166 ) = W( 166 ) + a*JVS( 966 )
  W( 167 ) = W( 167 ) + a*JVS( 967 )
  W( 168 ) = W( 168 ) + a*JVS( 968 )
  a = -W( 155 ) / JVS( 981 )
  W( 155 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 982 )
  W( 157 ) = W( 157 ) + a*JVS( 983 )
  W( 158 ) = W( 158 ) + a*JVS( 984 )
  W( 160 ) = W( 160 ) + a*JVS( 985 )
  W( 162 ) = W( 162 ) + a*JVS( 986 )
  W( 163 ) = W( 163 ) + a*JVS( 987 )
  W( 164 ) = W( 164 ) + a*JVS( 988 )
  W( 165 ) = W( 165 ) + a*JVS( 989 )
  W( 166 ) = W( 166 ) + a*JVS( 990 )
  W( 167 ) = W( 167 ) + a*JVS( 991 )
  W( 168 ) = W( 168 ) + a*JVS( 992 )
  a = -W( 156 ) / JVS( 1006 )
  W( 156 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 1007 )
  W( 159 ) = W( 159 ) + a*JVS( 1008 )
  W( 160 ) = W( 160 ) + a*JVS( 1009 )
  W( 161 ) = W( 161 ) + a*JVS( 1010 )
  W( 162 ) = W( 162 ) + a*JVS( 1011 )
  W( 163 ) = W( 163 ) + a*JVS( 1012 )
  W( 164 ) = W( 164 ) + a*JVS( 1013 )
  W( 165 ) = W( 165 ) + a*JVS( 1014 )
  W( 167 ) = W( 167 ) + a*JVS( 1015 )
  a = -W( 157 ) / JVS( 1045 )
  W( 157 ) = -a
  W( 158 ) = W( 158 ) + a*JVS( 1046 )
  W( 159 ) = W( 159 ) + a*JVS( 1047 )
  W( 160 ) = W( 160 ) + a*JVS( 1048 )
  W( 161 ) = W( 161 ) + a*JVS( 1049 )
  W( 162 ) = W( 162 ) + a*JVS( 1050 )
  W( 163 ) = W( 163 ) + a*JVS( 1051 )
  W( 164 ) = W( 164 ) + a*JVS( 1052 )
  W( 165 ) = W( 165 ) + a*JVS( 1053 )
  W( 166 ) = W( 166 ) + a*JVS( 1054 )
  W( 167 ) = W( 167 ) + a*JVS( 1055 )
  W( 168 ) = W( 168 ) + a*JVS( 1056 )
  a = -W( 158 ) / JVS( 1078 )
  W( 158 ) = -a
  W( 159 ) = W( 159 ) + a*JVS( 1079 )
  W( 160 ) = W( 160 ) + a*JVS( 1080 )
  W( 161 ) = W( 161 ) + a*JVS( 1081 )
  W( 162 ) = W( 162 ) + a*JVS( 1082 )
  W( 163 ) = W( 163 ) + a*JVS( 1083 )
  W( 164 ) = W( 164 ) + a*JVS( 1084 )
  W( 165 ) = W( 165 ) + a*JVS( 1085 )
  W( 166 ) = W( 166 ) + a*JVS( 1086 )
  W( 167 ) = W( 167 ) + a*JVS( 1087 )
  W( 168 ) = W( 168 ) + a*JVS( 1088 )
  a = -W( 159 ) / JVS( 1099 )
  W( 159 ) = -a
  W( 160 ) = W( 160 ) + a*JVS( 1100 )
  W( 161 ) = W( 161 ) + a*JVS( 1101 )
  W( 162 ) = W( 162 ) + a*JVS( 1102 )
  W( 163 ) = W( 163 ) + a*JVS( 1103 )
  W( 164 ) = W( 164 ) + a*JVS( 1104 )
  W( 165 ) = W( 165 ) + a*JVS( 1105 )
  W( 166 ) = W( 166 ) + a*JVS( 1106 )
  W( 167 ) = W( 167 ) + a*JVS( 1107 )
  W( 168 ) = W( 168 ) + a*JVS( 1108 )
  a = -W( 160 ) / JVS( 1125 )
  W( 160 ) = -a
  W( 161 ) = W( 161 ) + a*JVS( 1126 )
  W( 162 ) = W( 162 ) + a*JVS( 1127 )
  W( 163 ) = W( 163 ) + a*JVS( 1128 )
  W( 164 ) = W( 164 ) + a*JVS( 1129 )
  W( 165 ) = W( 165 ) + a*JVS( 1130 )
  W( 166 ) = W( 166 ) + a*JVS( 1131 )
  W( 167 ) = W( 167 ) + a*JVS( 1132 )
  W( 168 ) = W( 168 ) + a*JVS( 1133 )
  a = -W( 161 ) / JVS( 1177 )
  W( 161 ) = -a
  W( 162 ) = W( 162 ) + a*JVS( 1178 )
  W( 163 ) = W( 163 ) + a*JVS( 1179 )
  W( 164 ) = W( 164 ) + a*JVS( 1180 )
  W( 165 ) = W( 165 ) + a*JVS( 1181 )
  W( 166 ) = W( 166 ) + a*JVS( 1182 )
  W( 167 ) = W( 167 ) + a*JVS( 1183 )
  W( 168 ) = W( 168 ) + a*JVS( 1184 )
  a = -W( 162 ) / JVS( 1203 )
  W( 162 ) = -a
  W( 163 ) = W( 163 ) + a*JVS( 1204 )
  W( 164 ) = W( 164 ) + a*JVS( 1205 )
  W( 165 ) = W( 165 ) + a*JVS( 1206 )
  W( 166 ) = W( 166 ) + a*JVS( 1207 )
  W( 167 ) = W( 167 ) + a*JVS( 1208 )
  W( 168 ) = W( 168 ) + a*JVS( 1209 )
  a = -W( 163 ) / JVS( 1238 )
  W( 163 ) = -a
  W( 164 ) = W( 164 ) + a*JVS( 1239 )
  W( 165 ) = W( 165 ) + a*JVS( 1240 )
  W( 166 ) = W( 166 ) + a*JVS( 1241 )
  W( 167 ) = W( 167 ) + a*JVS( 1242 )
  W( 168 ) = W( 168 ) + a*JVS( 1243 )
  a = -W( 164 ) / JVS( 1256 )
  W( 164 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 1257 )
  W( 166 ) = W( 166 ) + a*JVS( 1258 )
  W( 167 ) = W( 167 ) + a*JVS( 1259 )
  W( 168 ) = W( 168 ) + a*JVS( 1260 )
  JVS( 1261) = W( 3 )
  JVS( 1262) = W( 5 )
  JVS( 1263) = W( 66 )
  JVS( 1264) = W( 67 )
  JVS( 1265) = W( 68 )
  JVS( 1266) = W( 69 )
  JVS( 1267) = W( 70 )
  JVS( 1268) = W( 71 )
  JVS( 1269) = W( 72 )
  JVS( 1270) = W( 73 )
  JVS( 1271) = W( 74 )
  JVS( 1272) = W( 75 )
  JVS( 1273) = W( 76 )
  JVS( 1274) = W( 77 )
  JVS( 1275) = W( 78 )
  JVS( 1276) = W( 79 )
  JVS( 1277) = W( 80 )
  JVS( 1278) = W( 81 )
  JVS( 1279) = W( 82 )
  JVS( 1280) = W( 83 )
  JVS( 1281) = W( 84 )
  JVS( 1282) = W( 85 )
  JVS( 1283) = W( 86 )
  JVS( 1284) = W( 87 )
  JVS( 1285) = W( 88 )
  JVS( 1286) = W( 89 )
  JVS( 1287) = W( 90 )
  JVS( 1288) = W( 91 )
  JVS( 1289) = W( 92 )
  JVS( 1290) = W( 93 )
  JVS( 1291) = W( 94 )
  JVS( 1292) = W( 95 )
  JVS( 1293) = W( 96 )
  JVS( 1294) = W( 97 )
  JVS( 1295) = W( 98 )
  JVS( 1296) = W( 99 )
  JVS( 1297) = W( 100 )
  JVS( 1298) = W( 101 )
  JVS( 1299) = W( 106 )
  JVS( 1300) = W( 107 )
  JVS( 1301) = W( 108 )
  JVS( 1302) = W( 111 )
  JVS( 1303) = W( 112 )
  JVS( 1304) = W( 113 )
  JVS( 1305) = W( 115 )
  JVS( 1306) = W( 116 )
  JVS( 1307) = W( 119 )
  JVS( 1308) = W( 120 )
  JVS( 1309) = W( 121 )
  JVS( 1310) = W( 122 )
  JVS( 1311) = W( 123 )
  JVS( 1312) = W( 124 )
  JVS( 1313) = W( 125 )
  JVS( 1314) = W( 126 )
  JVS( 1315) = W( 127 )
  JVS( 1316) = W( 129 )
  JVS( 1317) = W( 130 )
  JVS( 1318) = W( 131 )
  JVS( 1319) = W( 132 )
  JVS( 1320) = W( 133 )
  JVS( 1321) = W( 134 )
  JVS( 1322) = W( 135 )
  JVS( 1323) = W( 136 )
  JVS( 1324) = W( 137 )
  JVS( 1325) = W( 138 )
  JVS( 1326) = W( 139 )
  JVS( 1327) = W( 140 )
  JVS( 1328) = W( 142 )
  JVS( 1329) = W( 143 )
  JVS( 1330) = W( 144 )
  JVS( 1331) = W( 145 )
  JVS( 1332) = W( 146 )
  JVS( 1333) = W( 147 )
  JVS( 1334) = W( 148 )
  JVS( 1335) = W( 149 )
  JVS( 1336) = W( 150 )
  JVS( 1337) = W( 151 )
  JVS( 1338) = W( 152 )
  JVS( 1339) = W( 153 )
  JVS( 1340) = W( 154 )
  JVS( 1341) = W( 155 )
  JVS( 1342) = W( 156 )
  JVS( 1343) = W( 157 )
  JVS( 1344) = W( 158 )
  JVS( 1345) = W( 159 )
  JVS( 1346) = W( 160 )
  JVS( 1347) = W( 161 )
  JVS( 1348) = W( 162 )
  JVS( 1349) = W( 163 )
  JVS( 1350) = W( 164 )
  JVS( 1351) = W( 165 )
  JVS( 1352) = W( 166 )
  JVS( 1353) = W( 167 )
  JVS( 1354) = W( 168 )
  IF ( ABS( JVS( 1396 )) < TINY(a) ) THEN
         IER = 166
         RETURN
  END IF
   W( 101 ) = JVS( 1355 )
   W( 107 ) = JVS( 1356 )
   W( 108 ) = JVS( 1357 )
   W( 112 ) = JVS( 1358 )
   W( 113 ) = JVS( 1359 )
   W( 115 ) = JVS( 1360 )
   W( 121 ) = JVS( 1361 )
   W( 122 ) = JVS( 1362 )
   W( 123 ) = JVS( 1363 )
   W( 124 ) = JVS( 1364 )
   W( 125 ) = JVS( 1365 )
   W( 126 ) = JVS( 1366 )
   W( 127 ) = JVS( 1367 )
   W( 129 ) = JVS( 1368 )
   W( 131 ) = JVS( 1369 )
   W( 135 ) = JVS( 1370 )
   W( 137 ) = JVS( 1371 )
   W( 138 ) = JVS( 1372 )
   W( 139 ) = JVS( 1373 )
   W( 140 ) = JVS( 1374 )
   W( 142 ) = JVS( 1375 )
   W( 143 ) = JVS( 1376 )
   W( 144 ) = JVS( 1377 )
   W( 146 ) = JVS( 1378 )
   W( 147 ) = JVS( 1379 )
   W( 148 ) = JVS( 1380 )
   W( 151 ) = JVS( 1381 )
   W( 152 ) = JVS( 1382 )
   W( 153 ) = JVS( 1383 )
   W( 154 ) = JVS( 1384 )
   W( 155 ) = JVS( 1385 )
   W( 156 ) = JVS( 1386 )
   W( 157 ) = JVS( 1387 )
   W( 158 ) = JVS( 1388 )
   W( 159 ) = JVS( 1389 )
   W( 160 ) = JVS( 1390 )
   W( 161 ) = JVS( 1391 )
   W( 162 ) = JVS( 1392 )
   W( 163 ) = JVS( 1393 )
   W( 164 ) = JVS( 1394 )
   W( 165 ) = JVS( 1395 )
   W( 166 ) = JVS( 1396 )
   W( 167 ) = JVS( 1397 )
   W( 168 ) = JVS( 1398 )
  a = -W( 101 ) / JVS( 517 )
  W( 101 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 518 )
  a = -W( 107 ) / JVS( 534 )
  W( 107 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 535 )
  a = -W( 108 ) / JVS( 536 )
  W( 108 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 537 )
  a = -W( 112 ) / JVS( 550 )
  W( 112 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 551 )
  a = -W( 113 ) / JVS( 552 )
  W( 113 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 553 )
  a = -W( 115 ) / JVS( 558 )
  W( 115 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 559 )
  a = -W( 121 ) / JVS( 582 )
  W( 121 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 583 )
  a = -W( 122 ) / JVS( 586 )
  W( 122 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 587 )
  a = -W( 123 ) / JVS( 590 )
  W( 123 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 591 )
  a = -W( 124 ) / JVS( 594 )
  W( 124 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 595 )
  W( 165 ) = W( 165 ) + a*JVS( 596 )
  a = -W( 125 ) / JVS( 597 )
  W( 125 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 598 )
  a = -W( 126 ) / JVS( 599 )
  W( 126 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 600 )
  W( 165 ) = W( 165 ) + a*JVS( 601 )
  a = -W( 127 ) / JVS( 604 )
  W( 127 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 605 )
  W( 165 ) = W( 165 ) + a*JVS( 606 )
  a = -W( 129 ) / JVS( 613 )
  W( 129 ) = -a
  W( 158 ) = W( 158 ) + a*JVS( 614 )
  W( 161 ) = W( 161 ) + a*JVS( 615 )
  W( 165 ) = W( 165 ) + a*JVS( 616 )
  W( 166 ) = W( 166 ) + a*JVS( 617 )
  a = -W( 131 ) / JVS( 626 )
  W( 131 ) = -a
  W( 138 ) = W( 138 ) + a*JVS( 627 )
  W( 157 ) = W( 157 ) + a*JVS( 628 )
  W( 161 ) = W( 161 ) + a*JVS( 629 )
  W( 165 ) = W( 165 ) + a*JVS( 630 )
  a = -W( 135 ) / JVS( 689 )
  W( 135 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 690 )
  W( 156 ) = W( 156 ) + a*JVS( 691 )
  W( 157 ) = W( 157 ) + a*JVS( 692 )
  W( 165 ) = W( 165 ) + a*JVS( 693 )
  a = -W( 137 ) / JVS( 709 )
  W( 137 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 710 )
  W( 156 ) = W( 156 ) + a*JVS( 711 )
  W( 157 ) = W( 157 ) + a*JVS( 712 )
  W( 165 ) = W( 165 ) + a*JVS( 713 )
  a = -W( 138 ) / JVS( 716 )
  W( 138 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 717 )
  W( 159 ) = W( 159 ) + a*JVS( 718 )
  W( 160 ) = W( 160 ) + a*JVS( 719 )
  W( 161 ) = W( 161 ) + a*JVS( 720 )
  W( 162 ) = W( 162 ) + a*JVS( 721 )
  W( 163 ) = W( 163 ) + a*JVS( 722 )
  W( 164 ) = W( 164 ) + a*JVS( 723 )
  W( 165 ) = W( 165 ) + a*JVS( 724 )
  W( 167 ) = W( 167 ) + a*JVS( 725 )
  a = -W( 139 ) / JVS( 726 )
  W( 139 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 727 )
  W( 156 ) = W( 156 ) + a*JVS( 728 )
  W( 157 ) = W( 157 ) + a*JVS( 729 )
  W( 165 ) = W( 165 ) + a*JVS( 730 )
  a = -W( 140 ) / JVS( 731 )
  W( 140 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 732 )
  W( 156 ) = W( 156 ) + a*JVS( 733 )
  W( 157 ) = W( 157 ) + a*JVS( 734 )
  W( 165 ) = W( 165 ) + a*JVS( 735 )
  a = -W( 142 ) / JVS( 767 )
  W( 142 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 768 )
  W( 152 ) = W( 152 ) + a*JVS( 769 )
  W( 156 ) = W( 156 ) + a*JVS( 770 )
  W( 157 ) = W( 157 ) + a*JVS( 771 )
  W( 165 ) = W( 165 ) + a*JVS( 772 )
  a = -W( 143 ) / JVS( 773 )
  W( 143 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 774 )
  W( 156 ) = W( 156 ) + a*JVS( 775 )
  W( 157 ) = W( 157 ) + a*JVS( 776 )
  W( 165 ) = W( 165 ) + a*JVS( 777 )
  a = -W( 144 ) / JVS( 778 )
  W( 144 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 779 )
  W( 156 ) = W( 156 ) + a*JVS( 780 )
  W( 157 ) = W( 157 ) + a*JVS( 781 )
  W( 165 ) = W( 165 ) + a*JVS( 782 )
  a = -W( 146 ) / JVS( 808 )
  W( 146 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 809 )
  W( 152 ) = W( 152 ) + a*JVS( 810 )
  W( 156 ) = W( 156 ) + a*JVS( 811 )
  W( 157 ) = W( 157 ) + a*JVS( 812 )
  W( 165 ) = W( 165 ) + a*JVS( 813 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  a = -W( 148 ) / JVS( 821 )
  W( 148 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 822 )
  W( 156 ) = W( 156 ) + a*JVS( 823 )
  W( 157 ) = W( 157 ) + a*JVS( 824 )
  W( 165 ) = W( 165 ) + a*JVS( 825 )
  a = -W( 151 ) / JVS( 894 )
  W( 151 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 895 )
  W( 156 ) = W( 156 ) + a*JVS( 896 )
  W( 157 ) = W( 157 ) + a*JVS( 897 )
  W( 158 ) = W( 158 ) + a*JVS( 898 )
  W( 160 ) = W( 160 ) + a*JVS( 899 )
  W( 165 ) = W( 165 ) + a*JVS( 900 )
  W( 167 ) = W( 167 ) + a*JVS( 901 )
  a = -W( 152 ) / JVS( 912 )
  W( 152 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 913 )
  W( 157 ) = W( 157 ) + a*JVS( 914 )
  W( 160 ) = W( 160 ) + a*JVS( 915 )
  W( 165 ) = W( 165 ) + a*JVS( 916 )
  W( 167 ) = W( 167 ) + a*JVS( 917 )
  a = -W( 153 ) / JVS( 937 )
  W( 153 ) = -a
  W( 154 ) = W( 154 ) + a*JVS( 938 )
  W( 155 ) = W( 155 ) + a*JVS( 939 )
  W( 156 ) = W( 156 ) + a*JVS( 940 )
  W( 157 ) = W( 157 ) + a*JVS( 941 )
  W( 158 ) = W( 158 ) + a*JVS( 942 )
  W( 160 ) = W( 160 ) + a*JVS( 943 )
  W( 161 ) = W( 161 ) + a*JVS( 944 )
  W( 165 ) = W( 165 ) + a*JVS( 945 )
  W( 166 ) = W( 166 ) + a*JVS( 946 )
  W( 167 ) = W( 167 ) + a*JVS( 947 )
  a = -W( 154 ) / JVS( 959 )
  W( 154 ) = -a
  W( 155 ) = W( 155 ) + a*JVS( 960 )
  W( 156 ) = W( 156 ) + a*JVS( 961 )
  W( 157 ) = W( 157 ) + a*JVS( 962 )
  W( 158 ) = W( 158 ) + a*JVS( 963 )
  W( 160 ) = W( 160 ) + a*JVS( 964 )
  W( 165 ) = W( 165 ) + a*JVS( 965 )
  W( 166 ) = W( 166 ) + a*JVS( 966 )
  W( 167 ) = W( 167 ) + a*JVS( 967 )
  W( 168 ) = W( 168 ) + a*JVS( 968 )
  a = -W( 155 ) / JVS( 981 )
  W( 155 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 982 )
  W( 157 ) = W( 157 ) + a*JVS( 983 )
  W( 158 ) = W( 158 ) + a*JVS( 984 )
  W( 160 ) = W( 160 ) + a*JVS( 985 )
  W( 162 ) = W( 162 ) + a*JVS( 986 )
  W( 163 ) = W( 163 ) + a*JVS( 987 )
  W( 164 ) = W( 164 ) + a*JVS( 988 )
  W( 165 ) = W( 165 ) + a*JVS( 989 )
  W( 166 ) = W( 166 ) + a*JVS( 990 )
  W( 167 ) = W( 167 ) + a*JVS( 991 )
  W( 168 ) = W( 168 ) + a*JVS( 992 )
  a = -W( 156 ) / JVS( 1006 )
  W( 156 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 1007 )
  W( 159 ) = W( 159 ) + a*JVS( 1008 )
  W( 160 ) = W( 160 ) + a*JVS( 1009 )
  W( 161 ) = W( 161 ) + a*JVS( 1010 )
  W( 162 ) = W( 162 ) + a*JVS( 1011 )
  W( 163 ) = W( 163 ) + a*JVS( 1012 )
  W( 164 ) = W( 164 ) + a*JVS( 1013 )
  W( 165 ) = W( 165 ) + a*JVS( 1014 )
  W( 167 ) = W( 167 ) + a*JVS( 1015 )
  a = -W( 157 ) / JVS( 1045 )
  W( 157 ) = -a
  W( 158 ) = W( 158 ) + a*JVS( 1046 )
  W( 159 ) = W( 159 ) + a*JVS( 1047 )
  W( 160 ) = W( 160 ) + a*JVS( 1048 )
  W( 161 ) = W( 161 ) + a*JVS( 1049 )
  W( 162 ) = W( 162 ) + a*JVS( 1050 )
  W( 163 ) = W( 163 ) + a*JVS( 1051 )
  W( 164 ) = W( 164 ) + a*JVS( 1052 )
  W( 165 ) = W( 165 ) + a*JVS( 1053 )
  W( 166 ) = W( 166 ) + a*JVS( 1054 )
  W( 167 ) = W( 167 ) + a*JVS( 1055 )
  W( 168 ) = W( 168 ) + a*JVS( 1056 )
  a = -W( 158 ) / JVS( 1078 )
  W( 158 ) = -a
  W( 159 ) = W( 159 ) + a*JVS( 1079 )
  W( 160 ) = W( 160 ) + a*JVS( 1080 )
  W( 161 ) = W( 161 ) + a*JVS( 1081 )
  W( 162 ) = W( 162 ) + a*JVS( 1082 )
  W( 163 ) = W( 163 ) + a*JVS( 1083 )
  W( 164 ) = W( 164 ) + a*JVS( 1084 )
  W( 165 ) = W( 165 ) + a*JVS( 1085 )
  W( 166 ) = W( 166 ) + a*JVS( 1086 )
  W( 167 ) = W( 167 ) + a*JVS( 1087 )
  W( 168 ) = W( 168 ) + a*JVS( 1088 )
  a = -W( 159 ) / JVS( 1099 )
  W( 159 ) = -a
  W( 160 ) = W( 160 ) + a*JVS( 1100 )
  W( 161 ) = W( 161 ) + a*JVS( 1101 )
  W( 162 ) = W( 162 ) + a*JVS( 1102 )
  W( 163 ) = W( 163 ) + a*JVS( 1103 )
  W( 164 ) = W( 164 ) + a*JVS( 1104 )
  W( 165 ) = W( 165 ) + a*JVS( 1105 )
  W( 166 ) = W( 166 ) + a*JVS( 1106 )
  W( 167 ) = W( 167 ) + a*JVS( 1107 )
  W( 168 ) = W( 168 ) + a*JVS( 1108 )
  a = -W( 160 ) / JVS( 1125 )
  W( 160 ) = -a
  W( 161 ) = W( 161 ) + a*JVS( 1126 )
  W( 162 ) = W( 162 ) + a*JVS( 1127 )
  W( 163 ) = W( 163 ) + a*JVS( 1128 )
  W( 164 ) = W( 164 ) + a*JVS( 1129 )
  W( 165 ) = W( 165 ) + a*JVS( 1130 )
  W( 166 ) = W( 166 ) + a*JVS( 1131 )
  W( 167 ) = W( 167 ) + a*JVS( 1132 )
  W( 168 ) = W( 168 ) + a*JVS( 1133 )
  a = -W( 161 ) / JVS( 1177 )
  W( 161 ) = -a
  W( 162 ) = W( 162 ) + a*JVS( 1178 )
  W( 163 ) = W( 163 ) + a*JVS( 1179 )
  W( 164 ) = W( 164 ) + a*JVS( 1180 )
  W( 165 ) = W( 165 ) + a*JVS( 1181 )
  W( 166 ) = W( 166 ) + a*JVS( 1182 )
  W( 167 ) = W( 167 ) + a*JVS( 1183 )
  W( 168 ) = W( 168 ) + a*JVS( 1184 )
  a = -W( 162 ) / JVS( 1203 )
  W( 162 ) = -a
  W( 163 ) = W( 163 ) + a*JVS( 1204 )
  W( 164 ) = W( 164 ) + a*JVS( 1205 )
  W( 165 ) = W( 165 ) + a*JVS( 1206 )
  W( 166 ) = W( 166 ) + a*JVS( 1207 )
  W( 167 ) = W( 167 ) + a*JVS( 1208 )
  W( 168 ) = W( 168 ) + a*JVS( 1209 )
  a = -W( 163 ) / JVS( 1238 )
  W( 163 ) = -a
  W( 164 ) = W( 164 ) + a*JVS( 1239 )
  W( 165 ) = W( 165 ) + a*JVS( 1240 )
  W( 166 ) = W( 166 ) + a*JVS( 1241 )
  W( 167 ) = W( 167 ) + a*JVS( 1242 )
  W( 168 ) = W( 168 ) + a*JVS( 1243 )
  a = -W( 164 ) / JVS( 1256 )
  W( 164 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 1257 )
  W( 166 ) = W( 166 ) + a*JVS( 1258 )
  W( 167 ) = W( 167 ) + a*JVS( 1259 )
  W( 168 ) = W( 168 ) + a*JVS( 1260 )
  a = -W( 165 ) / JVS( 1351 )
  W( 165 ) = -a
  W( 166 ) = W( 166 ) + a*JVS( 1352 )
  W( 167 ) = W( 167 ) + a*JVS( 1353 )
  W( 168 ) = W( 168 ) + a*JVS( 1354 )
  JVS( 1355) = W( 101 )
  JVS( 1356) = W( 107 )
  JVS( 1357) = W( 108 )
  JVS( 1358) = W( 112 )
  JVS( 1359) = W( 113 )
  JVS( 1360) = W( 115 )
  JVS( 1361) = W( 121 )
  JVS( 1362) = W( 122 )
  JVS( 1363) = W( 123 )
  JVS( 1364) = W( 124 )
  JVS( 1365) = W( 125 )
  JVS( 1366) = W( 126 )
  JVS( 1367) = W( 127 )
  JVS( 1368) = W( 129 )
  JVS( 1369) = W( 131 )
  JVS( 1370) = W( 135 )
  JVS( 1371) = W( 137 )
  JVS( 1372) = W( 138 )
  JVS( 1373) = W( 139 )
  JVS( 1374) = W( 140 )
  JVS( 1375) = W( 142 )
  JVS( 1376) = W( 143 )
  JVS( 1377) = W( 144 )
  JVS( 1378) = W( 146 )
  JVS( 1379) = W( 147 )
  JVS( 1380) = W( 148 )
  JVS( 1381) = W( 151 )
  JVS( 1382) = W( 152 )
  JVS( 1383) = W( 153 )
  JVS( 1384) = W( 154 )
  JVS( 1385) = W( 155 )
  JVS( 1386) = W( 156 )
  JVS( 1387) = W( 157 )
  JVS( 1388) = W( 158 )
  JVS( 1389) = W( 159 )
  JVS( 1390) = W( 160 )
  JVS( 1391) = W( 161 )
  JVS( 1392) = W( 162 )
  JVS( 1393) = W( 163 )
  JVS( 1394) = W( 164 )
  JVS( 1395) = W( 165 )
  JVS( 1396) = W( 166 )
  JVS( 1397) = W( 167 )
  JVS( 1398) = W( 168 )
  IF ( ABS( JVS( 1439 )) < TINY(a) ) THEN
         IER = 167
         RETURN
  END IF
   W( 102 ) = JVS( 1399 )
   W( 103 ) = JVS( 1400 )
   W( 104 ) = JVS( 1401 )
   W( 105 ) = JVS( 1402 )
   W( 110 ) = JVS( 1403 )
   W( 111 ) = JVS( 1404 )
   W( 114 ) = JVS( 1405 )
   W( 117 ) = JVS( 1406 )
   W( 118 ) = JVS( 1407 )
   W( 119 ) = JVS( 1408 )
   W( 128 ) = JVS( 1409 )
   W( 134 ) = JVS( 1410 )
   W( 138 ) = JVS( 1411 )
   W( 140 ) = JVS( 1412 )
   W( 141 ) = JVS( 1413 )
   W( 142 ) = JVS( 1414 )
   W( 143 ) = JVS( 1415 )
   W( 144 ) = JVS( 1416 )
   W( 145 ) = JVS( 1417 )
   W( 146 ) = JVS( 1418 )
   W( 147 ) = JVS( 1419 )
   W( 148 ) = JVS( 1420 )
   W( 149 ) = JVS( 1421 )
   W( 150 ) = JVS( 1422 )
   W( 151 ) = JVS( 1423 )
   W( 152 ) = JVS( 1424 )
   W( 153 ) = JVS( 1425 )
   W( 154 ) = JVS( 1426 )
   W( 155 ) = JVS( 1427 )
   W( 156 ) = JVS( 1428 )
   W( 157 ) = JVS( 1429 )
   W( 158 ) = JVS( 1430 )
   W( 159 ) = JVS( 1431 )
   W( 160 ) = JVS( 1432 )
   W( 161 ) = JVS( 1433 )
   W( 162 ) = JVS( 1434 )
   W( 163 ) = JVS( 1435 )
   W( 164 ) = JVS( 1436 )
   W( 165 ) = JVS( 1437 )
   W( 166 ) = JVS( 1438 )
   W( 167 ) = JVS( 1439 )
   W( 168 ) = JVS( 1440 )
  a = -W( 102 ) / JVS( 519 )
  W( 102 ) = -a
  W( 162 ) = W( 162 ) + a*JVS( 520 )
  W( 167 ) = W( 167 ) + a*JVS( 521 )
  a = -W( 103 ) / JVS( 522 )
  W( 103 ) = -a
  W( 159 ) = W( 159 ) + a*JVS( 523 )
  W( 167 ) = W( 167 ) + a*JVS( 524 )
  a = -W( 104 ) / JVS( 525 )
  W( 104 ) = -a
  W( 163 ) = W( 163 ) + a*JVS( 526 )
  W( 167 ) = W( 167 ) + a*JVS( 527 )
  a = -W( 105 ) / JVS( 528 )
  W( 105 ) = -a
  W( 164 ) = W( 164 ) + a*JVS( 529 )
  W( 167 ) = W( 167 ) + a*JVS( 530 )
  a = -W( 110 ) / JVS( 544 )
  W( 110 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 545 )
  W( 167 ) = W( 167 ) + a*JVS( 546 )
  a = -W( 111 ) / JVS( 547 )
  W( 111 ) = -a
  W( 160 ) = W( 160 ) + a*JVS( 548 )
  W( 165 ) = W( 165 ) + a*JVS( 549 )
  a = -W( 114 ) / JVS( 555 )
  W( 114 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 556 )
  W( 167 ) = W( 167 ) + a*JVS( 557 )
  a = -W( 117 ) / JVS( 564 )
  W( 117 ) = -a
  W( 150 ) = W( 150 ) + a*JVS( 565 )
  W( 160 ) = W( 160 ) + a*JVS( 566 )
  W( 161 ) = W( 161 ) + a*JVS( 567 )
  a = -W( 118 ) / JVS( 568 )
  W( 118 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 569 )
  W( 157 ) = W( 157 ) + a*JVS( 570 )
  W( 161 ) = W( 161 ) + a*JVS( 571 )
  W( 167 ) = W( 167 ) + a*JVS( 572 )
  a = -W( 119 ) / JVS( 573 )
  W( 119 ) = -a
  W( 161 ) = W( 161 ) + a*JVS( 574 )
  W( 165 ) = W( 165 ) + a*JVS( 575 )
  W( 167 ) = W( 167 ) + a*JVS( 576 )
  a = -W( 128 ) / JVS( 608 )
  W( 128 ) = -a
  W( 138 ) = W( 138 ) + a*JVS( 609 )
  W( 157 ) = W( 157 ) + a*JVS( 610 )
  W( 161 ) = W( 161 ) + a*JVS( 611 )
  W( 167 ) = W( 167 ) + a*JVS( 612 )
  a = -W( 134 ) / JVS( 674 )
  W( 134 ) = -a
  W( 138 ) = W( 138 ) + a*JVS( 675 )
  W( 142 ) = W( 142 ) + a*JVS( 676 )
  W( 145 ) = W( 145 ) + a*JVS( 677 )
  W( 146 ) = W( 146 ) + a*JVS( 678 )
  W( 147 ) = W( 147 ) + a*JVS( 679 )
  W( 148 ) = W( 148 ) + a*JVS( 680 )
  W( 149 ) = W( 149 ) + a*JVS( 681 )
  W( 150 ) = W( 150 ) + a*JVS( 682 )
  W( 153 ) = W( 153 ) + a*JVS( 683 )
  W( 156 ) = W( 156 ) + a*JVS( 684 )
  W( 157 ) = W( 157 ) + a*JVS( 685 )
  W( 161 ) = W( 161 ) + a*JVS( 686 )
  W( 165 ) = W( 165 ) + a*JVS( 687 )
  W( 167 ) = W( 167 ) + a*JVS( 688 )
  a = -W( 138 ) / JVS( 716 )
  W( 138 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 717 )
  W( 159 ) = W( 159 ) + a*JVS( 718 )
  W( 160 ) = W( 160 ) + a*JVS( 719 )
  W( 161 ) = W( 161 ) + a*JVS( 720 )
  W( 162 ) = W( 162 ) + a*JVS( 721 )
  W( 163 ) = W( 163 ) + a*JVS( 722 )
  W( 164 ) = W( 164 ) + a*JVS( 723 )
  W( 165 ) = W( 165 ) + a*JVS( 724 )
  W( 167 ) = W( 167 ) + a*JVS( 725 )
  a = -W( 140 ) / JVS( 731 )
  W( 140 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 732 )
  W( 156 ) = W( 156 ) + a*JVS( 733 )
  W( 157 ) = W( 157 ) + a*JVS( 734 )
  W( 165 ) = W( 165 ) + a*JVS( 735 )
  a = -W( 141 ) / JVS( 744 )
  W( 141 ) = -a
  W( 143 ) = W( 143 ) + a*JVS( 745 )
  W( 144 ) = W( 144 ) + a*JVS( 746 )
  W( 147 ) = W( 147 ) + a*JVS( 747 )
  W( 148 ) = W( 148 ) + a*JVS( 748 )
  W( 151 ) = W( 151 ) + a*JVS( 749 )
  W( 152 ) = W( 152 ) + a*JVS( 750 )
  W( 154 ) = W( 154 ) + a*JVS( 751 )
  W( 155 ) = W( 155 ) + a*JVS( 752 )
  W( 156 ) = W( 156 ) + a*JVS( 753 )
  W( 157 ) = W( 157 ) + a*JVS( 754 )
  W( 158 ) = W( 158 ) + a*JVS( 755 )
  W( 159 ) = W( 159 ) + a*JVS( 756 )
  W( 160 ) = W( 160 ) + a*JVS( 757 )
  W( 161 ) = W( 161 ) + a*JVS( 758 )
  W( 162 ) = W( 162 ) + a*JVS( 759 )
  W( 163 ) = W( 163 ) + a*JVS( 760 )
  W( 164 ) = W( 164 ) + a*JVS( 761 )
  W( 165 ) = W( 165 ) + a*JVS( 762 )
  W( 166 ) = W( 166 ) + a*JVS( 763 )
  W( 167 ) = W( 167 ) + a*JVS( 764 )
  W( 168 ) = W( 168 ) + a*JVS( 765 )
  a = -W( 142 ) / JVS( 767 )
  W( 142 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 768 )
  W( 152 ) = W( 152 ) + a*JVS( 769 )
  W( 156 ) = W( 156 ) + a*JVS( 770 )
  W( 157 ) = W( 157 ) + a*JVS( 771 )
  W( 165 ) = W( 165 ) + a*JVS( 772 )
  a = -W( 143 ) / JVS( 773 )
  W( 143 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 774 )
  W( 156 ) = W( 156 ) + a*JVS( 775 )
  W( 157 ) = W( 157 ) + a*JVS( 776 )
  W( 165 ) = W( 165 ) + a*JVS( 777 )
  a = -W( 144 ) / JVS( 778 )
  W( 144 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 779 )
  W( 156 ) = W( 156 ) + a*JVS( 780 )
  W( 157 ) = W( 157 ) + a*JVS( 781 )
  W( 165 ) = W( 165 ) + a*JVS( 782 )
  a = -W( 145 ) / JVS( 794 )
  W( 145 ) = -a
  W( 146 ) = W( 146 ) + a*JVS( 795 )
  W( 152 ) = W( 152 ) + a*JVS( 796 )
  W( 156 ) = W( 156 ) + a*JVS( 797 )
  W( 157 ) = W( 157 ) + a*JVS( 798 )
  W( 159 ) = W( 159 ) + a*JVS( 799 )
  W( 160 ) = W( 160 ) + a*JVS( 800 )
  W( 161 ) = W( 161 ) + a*JVS( 801 )
  W( 162 ) = W( 162 ) + a*JVS( 802 )
  W( 163 ) = W( 163 ) + a*JVS( 803 )
  W( 164 ) = W( 164 ) + a*JVS( 804 )
  W( 165 ) = W( 165 ) + a*JVS( 805 )
  W( 167 ) = W( 167 ) + a*JVS( 806 )
  a = -W( 146 ) / JVS( 808 )
  W( 146 ) = -a
  W( 147 ) = W( 147 ) + a*JVS( 809 )
  W( 152 ) = W( 152 ) + a*JVS( 810 )
  W( 156 ) = W( 156 ) + a*JVS( 811 )
  W( 157 ) = W( 157 ) + a*JVS( 812 )
  W( 165 ) = W( 165 ) + a*JVS( 813 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  a = -W( 148 ) / JVS( 821 )
  W( 148 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 822 )
  W( 156 ) = W( 156 ) + a*JVS( 823 )
  W( 157 ) = W( 157 ) + a*JVS( 824 )
  W( 165 ) = W( 165 ) + a*JVS( 825 )
  a = -W( 149 ) / JVS( 836 )
  W( 149 ) = -a
  W( 151 ) = W( 151 ) + a*JVS( 837 )
  W( 152 ) = W( 152 ) + a*JVS( 838 )
  W( 153 ) = W( 153 ) + a*JVS( 839 )
  W( 154 ) = W( 154 ) + a*JVS( 840 )
  W( 155 ) = W( 155 ) + a*JVS( 841 )
  W( 156 ) = W( 156 ) + a*JVS( 842 )
  W( 157 ) = W( 157 ) + a*JVS( 843 )
  W( 159 ) = W( 159 ) + a*JVS( 844 )
  W( 160 ) = W( 160 ) + a*JVS( 845 )
  W( 162 ) = W( 162 ) + a*JVS( 846 )
  W( 163 ) = W( 163 ) + a*JVS( 847 )
  W( 164 ) = W( 164 ) + a*JVS( 848 )
  W( 165 ) = W( 165 ) + a*JVS( 849 )
  a = -W( 150 ) / JVS( 870 )
  W( 150 ) = -a
  W( 151 ) = W( 151 ) + a*JVS( 871 )
  W( 152 ) = W( 152 ) + a*JVS( 872 )
  W( 154 ) = W( 154 ) + a*JVS( 873 )
  W( 155 ) = W( 155 ) + a*JVS( 874 )
  W( 156 ) = W( 156 ) + a*JVS( 875 )
  W( 157 ) = W( 157 ) + a*JVS( 876 )
  W( 158 ) = W( 158 ) + a*JVS( 877 )
  W( 159 ) = W( 159 ) + a*JVS( 878 )
  W( 160 ) = W( 160 ) + a*JVS( 879 )
  W( 161 ) = W( 161 ) + a*JVS( 880 )
  W( 162 ) = W( 162 ) + a*JVS( 881 )
  W( 163 ) = W( 163 ) + a*JVS( 882 )
  W( 164 ) = W( 164 ) + a*JVS( 883 )
  W( 165 ) = W( 165 ) + a*JVS( 884 )
  W( 166 ) = W( 166 ) + a*JVS( 885 )
  W( 167 ) = W( 167 ) + a*JVS( 886 )
  W( 168 ) = W( 168 ) + a*JVS( 887 )
  a = -W( 151 ) / JVS( 894 )
  W( 151 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 895 )
  W( 156 ) = W( 156 ) + a*JVS( 896 )
  W( 157 ) = W( 157 ) + a*JVS( 897 )
  W( 158 ) = W( 158 ) + a*JVS( 898 )
  W( 160 ) = W( 160 ) + a*JVS( 899 )
  W( 165 ) = W( 165 ) + a*JVS( 900 )
  W( 167 ) = W( 167 ) + a*JVS( 901 )
  a = -W( 152 ) / JVS( 912 )
  W( 152 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 913 )
  W( 157 ) = W( 157 ) + a*JVS( 914 )
  W( 160 ) = W( 160 ) + a*JVS( 915 )
  W( 165 ) = W( 165 ) + a*JVS( 916 )
  W( 167 ) = W( 167 ) + a*JVS( 917 )
  a = -W( 153 ) / JVS( 937 )
  W( 153 ) = -a
  W( 154 ) = W( 154 ) + a*JVS( 938 )
  W( 155 ) = W( 155 ) + a*JVS( 939 )
  W( 156 ) = W( 156 ) + a*JVS( 940 )
  W( 157 ) = W( 157 ) + a*JVS( 941 )
  W( 158 ) = W( 158 ) + a*JVS( 942 )
  W( 160 ) = W( 160 ) + a*JVS( 943 )
  W( 161 ) = W( 161 ) + a*JVS( 944 )
  W( 165 ) = W( 165 ) + a*JVS( 945 )
  W( 166 ) = W( 166 ) + a*JVS( 946 )
  W( 167 ) = W( 167 ) + a*JVS( 947 )
  a = -W( 154 ) / JVS( 959 )
  W( 154 ) = -a
  W( 155 ) = W( 155 ) + a*JVS( 960 )
  W( 156 ) = W( 156 ) + a*JVS( 961 )
  W( 157 ) = W( 157 ) + a*JVS( 962 )
  W( 158 ) = W( 158 ) + a*JVS( 963 )
  W( 160 ) = W( 160 ) + a*JVS( 964 )
  W( 165 ) = W( 165 ) + a*JVS( 965 )
  W( 166 ) = W( 166 ) + a*JVS( 966 )
  W( 167 ) = W( 167 ) + a*JVS( 967 )
  W( 168 ) = W( 168 ) + a*JVS( 968 )
  a = -W( 155 ) / JVS( 981 )
  W( 155 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 982 )
  W( 157 ) = W( 157 ) + a*JVS( 983 )
  W( 158 ) = W( 158 ) + a*JVS( 984 )
  W( 160 ) = W( 160 ) + a*JVS( 985 )
  W( 162 ) = W( 162 ) + a*JVS( 986 )
  W( 163 ) = W( 163 ) + a*JVS( 987 )
  W( 164 ) = W( 164 ) + a*JVS( 988 )
  W( 165 ) = W( 165 ) + a*JVS( 989 )
  W( 166 ) = W( 166 ) + a*JVS( 990 )
  W( 167 ) = W( 167 ) + a*JVS( 991 )
  W( 168 ) = W( 168 ) + a*JVS( 992 )
  a = -W( 156 ) / JVS( 1006 )
  W( 156 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 1007 )
  W( 159 ) = W( 159 ) + a*JVS( 1008 )
  W( 160 ) = W( 160 ) + a*JVS( 1009 )
  W( 161 ) = W( 161 ) + a*JVS( 1010 )
  W( 162 ) = W( 162 ) + a*JVS( 1011 )
  W( 163 ) = W( 163 ) + a*JVS( 1012 )
  W( 164 ) = W( 164 ) + a*JVS( 1013 )
  W( 165 ) = W( 165 ) + a*JVS( 1014 )
  W( 167 ) = W( 167 ) + a*JVS( 1015 )
  a = -W( 157 ) / JVS( 1045 )
  W( 157 ) = -a
  W( 158 ) = W( 158 ) + a*JVS( 1046 )
  W( 159 ) = W( 159 ) + a*JVS( 1047 )
  W( 160 ) = W( 160 ) + a*JVS( 1048 )
  W( 161 ) = W( 161 ) + a*JVS( 1049 )
  W( 162 ) = W( 162 ) + a*JVS( 1050 )
  W( 163 ) = W( 163 ) + a*JVS( 1051 )
  W( 164 ) = W( 164 ) + a*JVS( 1052 )
  W( 165 ) = W( 165 ) + a*JVS( 1053 )
  W( 166 ) = W( 166 ) + a*JVS( 1054 )
  W( 167 ) = W( 167 ) + a*JVS( 1055 )
  W( 168 ) = W( 168 ) + a*JVS( 1056 )
  a = -W( 158 ) / JVS( 1078 )
  W( 158 ) = -a
  W( 159 ) = W( 159 ) + a*JVS( 1079 )
  W( 160 ) = W( 160 ) + a*JVS( 1080 )
  W( 161 ) = W( 161 ) + a*JVS( 1081 )
  W( 162 ) = W( 162 ) + a*JVS( 1082 )
  W( 163 ) = W( 163 ) + a*JVS( 1083 )
  W( 164 ) = W( 164 ) + a*JVS( 1084 )
  W( 165 ) = W( 165 ) + a*JVS( 1085 )
  W( 166 ) = W( 166 ) + a*JVS( 1086 )
  W( 167 ) = W( 167 ) + a*JVS( 1087 )
  W( 168 ) = W( 168 ) + a*JVS( 1088 )
  a = -W( 159 ) / JVS( 1099 )
  W( 159 ) = -a
  W( 160 ) = W( 160 ) + a*JVS( 1100 )
  W( 161 ) = W( 161 ) + a*JVS( 1101 )
  W( 162 ) = W( 162 ) + a*JVS( 1102 )
  W( 163 ) = W( 163 ) + a*JVS( 1103 )
  W( 164 ) = W( 164 ) + a*JVS( 1104 )
  W( 165 ) = W( 165 ) + a*JVS( 1105 )
  W( 166 ) = W( 166 ) + a*JVS( 1106 )
  W( 167 ) = W( 167 ) + a*JVS( 1107 )
  W( 168 ) = W( 168 ) + a*JVS( 1108 )
  a = -W( 160 ) / JVS( 1125 )
  W( 160 ) = -a
  W( 161 ) = W( 161 ) + a*JVS( 1126 )
  W( 162 ) = W( 162 ) + a*JVS( 1127 )
  W( 163 ) = W( 163 ) + a*JVS( 1128 )
  W( 164 ) = W( 164 ) + a*JVS( 1129 )
  W( 165 ) = W( 165 ) + a*JVS( 1130 )
  W( 166 ) = W( 166 ) + a*JVS( 1131 )
  W( 167 ) = W( 167 ) + a*JVS( 1132 )
  W( 168 ) = W( 168 ) + a*JVS( 1133 )
  a = -W( 161 ) / JVS( 1177 )
  W( 161 ) = -a
  W( 162 ) = W( 162 ) + a*JVS( 1178 )
  W( 163 ) = W( 163 ) + a*JVS( 1179 )
  W( 164 ) = W( 164 ) + a*JVS( 1180 )
  W( 165 ) = W( 165 ) + a*JVS( 1181 )
  W( 166 ) = W( 166 ) + a*JVS( 1182 )
  W( 167 ) = W( 167 ) + a*JVS( 1183 )
  W( 168 ) = W( 168 ) + a*JVS( 1184 )
  a = -W( 162 ) / JVS( 1203 )
  W( 162 ) = -a
  W( 163 ) = W( 163 ) + a*JVS( 1204 )
  W( 164 ) = W( 164 ) + a*JVS( 1205 )
  W( 165 ) = W( 165 ) + a*JVS( 1206 )
  W( 166 ) = W( 166 ) + a*JVS( 1207 )
  W( 167 ) = W( 167 ) + a*JVS( 1208 )
  W( 168 ) = W( 168 ) + a*JVS( 1209 )
  a = -W( 163 ) / JVS( 1238 )
  W( 163 ) = -a
  W( 164 ) = W( 164 ) + a*JVS( 1239 )
  W( 165 ) = W( 165 ) + a*JVS( 1240 )
  W( 166 ) = W( 166 ) + a*JVS( 1241 )
  W( 167 ) = W( 167 ) + a*JVS( 1242 )
  W( 168 ) = W( 168 ) + a*JVS( 1243 )
  a = -W( 164 ) / JVS( 1256 )
  W( 164 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 1257 )
  W( 166 ) = W( 166 ) + a*JVS( 1258 )
  W( 167 ) = W( 167 ) + a*JVS( 1259 )
  W( 168 ) = W( 168 ) + a*JVS( 1260 )
  a = -W( 165 ) / JVS( 1351 )
  W( 165 ) = -a
  W( 166 ) = W( 166 ) + a*JVS( 1352 )
  W( 167 ) = W( 167 ) + a*JVS( 1353 )
  W( 168 ) = W( 168 ) + a*JVS( 1354 )
  a = -W( 166 ) / JVS( 1396 )
  W( 166 ) = -a
  W( 167 ) = W( 167 ) + a*JVS( 1397 )
  W( 168 ) = W( 168 ) + a*JVS( 1398 )
  JVS( 1399) = W( 102 )
  JVS( 1400) = W( 103 )
  JVS( 1401) = W( 104 )
  JVS( 1402) = W( 105 )
  JVS( 1403) = W( 110 )
  JVS( 1404) = W( 111 )
  JVS( 1405) = W( 114 )
  JVS( 1406) = W( 117 )
  JVS( 1407) = W( 118 )
  JVS( 1408) = W( 119 )
  JVS( 1409) = W( 128 )
  JVS( 1410) = W( 134 )
  JVS( 1411) = W( 138 )
  JVS( 1412) = W( 140 )
  JVS( 1413) = W( 141 )
  JVS( 1414) = W( 142 )
  JVS( 1415) = W( 143 )
  JVS( 1416) = W( 144 )
  JVS( 1417) = W( 145 )
  JVS( 1418) = W( 146 )
  JVS( 1419) = W( 147 )
  JVS( 1420) = W( 148 )
  JVS( 1421) = W( 149 )
  JVS( 1422) = W( 150 )
  JVS( 1423) = W( 151 )
  JVS( 1424) = W( 152 )
  JVS( 1425) = W( 153 )
  JVS( 1426) = W( 154 )
  JVS( 1427) = W( 155 )
  JVS( 1428) = W( 156 )
  JVS( 1429) = W( 157 )
  JVS( 1430) = W( 158 )
  JVS( 1431) = W( 159 )
  JVS( 1432) = W( 160 )
  JVS( 1433) = W( 161 )
  JVS( 1434) = W( 162 )
  JVS( 1435) = W( 163 )
  JVS( 1436) = W( 164 )
  JVS( 1437) = W( 165 )
  JVS( 1438) = W( 166 )
  JVS( 1439) = W( 167 )
  JVS( 1440) = W( 168 )
  IF ( ABS( JVS( 1472 )) < TINY(a) ) THEN
         IER = 168
         RETURN
  END IF
   W( 99 ) = JVS( 1441 )
   W( 114 ) = JVS( 1442 )
   W( 116 ) = JVS( 1443 )
   W( 125 ) = JVS( 1444 )
   W( 135 ) = JVS( 1445 )
   W( 136 ) = JVS( 1446 )
   W( 137 ) = JVS( 1447 )
   W( 139 ) = JVS( 1448 )
   W( 140 ) = JVS( 1449 )
   W( 143 ) = JVS( 1450 )
   W( 144 ) = JVS( 1451 )
   W( 147 ) = JVS( 1452 )
   W( 148 ) = JVS( 1453 )
   W( 149 ) = JVS( 1454 )
   W( 151 ) = JVS( 1455 )
   W( 152 ) = JVS( 1456 )
   W( 153 ) = JVS( 1457 )
   W( 154 ) = JVS( 1458 )
   W( 155 ) = JVS( 1459 )
   W( 156 ) = JVS( 1460 )
   W( 157 ) = JVS( 1461 )
   W( 158 ) = JVS( 1462 )
   W( 159 ) = JVS( 1463 )
   W( 160 ) = JVS( 1464 )
   W( 161 ) = JVS( 1465 )
   W( 162 ) = JVS( 1466 )
   W( 163 ) = JVS( 1467 )
   W( 164 ) = JVS( 1468 )
   W( 165 ) = JVS( 1469 )
   W( 166 ) = JVS( 1470 )
   W( 167 ) = JVS( 1471 )
   W( 168 ) = JVS( 1472 )
  a = -W( 99 ) / JVS( 513 )
  W( 99 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 514 )
  a = -W( 114 ) / JVS( 555 )
  W( 114 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 556 )
  W( 167 ) = W( 167 ) + a*JVS( 557 )
  a = -W( 116 ) / JVS( 560 )
  W( 116 ) = -a
  W( 161 ) = W( 161 ) + a*JVS( 561 )
  W( 165 ) = W( 165 ) + a*JVS( 562 )
  W( 168 ) = W( 168 ) + a*JVS( 563 )
  a = -W( 125 ) / JVS( 597 )
  W( 125 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 598 )
  a = -W( 135 ) / JVS( 689 )
  W( 135 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 690 )
  W( 156 ) = W( 156 ) + a*JVS( 691 )
  W( 157 ) = W( 157 ) + a*JVS( 692 )
  W( 165 ) = W( 165 ) + a*JVS( 693 )
  a = -W( 136 ) / JVS( 699 )
  W( 136 ) = -a
  W( 139 ) = W( 139 ) + a*JVS( 700 )
  W( 143 ) = W( 143 ) + a*JVS( 701 )
  W( 144 ) = W( 144 ) + a*JVS( 702 )
  W( 147 ) = W( 147 ) + a*JVS( 703 )
  W( 151 ) = W( 151 ) + a*JVS( 704 )
  W( 156 ) = W( 156 ) + a*JVS( 705 )
  W( 157 ) = W( 157 ) + a*JVS( 706 )
  W( 165 ) = W( 165 ) + a*JVS( 707 )
  W( 167 ) = W( 167 ) + a*JVS( 708 )
  a = -W( 137 ) / JVS( 709 )
  W( 137 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 710 )
  W( 156 ) = W( 156 ) + a*JVS( 711 )
  W( 157 ) = W( 157 ) + a*JVS( 712 )
  W( 165 ) = W( 165 ) + a*JVS( 713 )
  a = -W( 139 ) / JVS( 726 )
  W( 139 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 727 )
  W( 156 ) = W( 156 ) + a*JVS( 728 )
  W( 157 ) = W( 157 ) + a*JVS( 729 )
  W( 165 ) = W( 165 ) + a*JVS( 730 )
  a = -W( 140 ) / JVS( 731 )
  W( 140 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 732 )
  W( 156 ) = W( 156 ) + a*JVS( 733 )
  W( 157 ) = W( 157 ) + a*JVS( 734 )
  W( 165 ) = W( 165 ) + a*JVS( 735 )
  a = -W( 143 ) / JVS( 773 )
  W( 143 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 774 )
  W( 156 ) = W( 156 ) + a*JVS( 775 )
  W( 157 ) = W( 157 ) + a*JVS( 776 )
  W( 165 ) = W( 165 ) + a*JVS( 777 )
  a = -W( 144 ) / JVS( 778 )
  W( 144 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 779 )
  W( 156 ) = W( 156 ) + a*JVS( 780 )
  W( 157 ) = W( 157 ) + a*JVS( 781 )
  W( 165 ) = W( 165 ) + a*JVS( 782 )
  a = -W( 147 ) / JVS( 814 )
  W( 147 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 815 )
  W( 156 ) = W( 156 ) + a*JVS( 816 )
  W( 157 ) = W( 157 ) + a*JVS( 817 )
  W( 165 ) = W( 165 ) + a*JVS( 818 )
  a = -W( 148 ) / JVS( 821 )
  W( 148 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 822 )
  W( 156 ) = W( 156 ) + a*JVS( 823 )
  W( 157 ) = W( 157 ) + a*JVS( 824 )
  W( 165 ) = W( 165 ) + a*JVS( 825 )
  a = -W( 149 ) / JVS( 836 )
  W( 149 ) = -a
  W( 151 ) = W( 151 ) + a*JVS( 837 )
  W( 152 ) = W( 152 ) + a*JVS( 838 )
  W( 153 ) = W( 153 ) + a*JVS( 839 )
  W( 154 ) = W( 154 ) + a*JVS( 840 )
  W( 155 ) = W( 155 ) + a*JVS( 841 )
  W( 156 ) = W( 156 ) + a*JVS( 842 )
  W( 157 ) = W( 157 ) + a*JVS( 843 )
  W( 159 ) = W( 159 ) + a*JVS( 844 )
  W( 160 ) = W( 160 ) + a*JVS( 845 )
  W( 162 ) = W( 162 ) + a*JVS( 846 )
  W( 163 ) = W( 163 ) + a*JVS( 847 )
  W( 164 ) = W( 164 ) + a*JVS( 848 )
  W( 165 ) = W( 165 ) + a*JVS( 849 )
  a = -W( 151 ) / JVS( 894 )
  W( 151 ) = -a
  W( 152 ) = W( 152 ) + a*JVS( 895 )
  W( 156 ) = W( 156 ) + a*JVS( 896 )
  W( 157 ) = W( 157 ) + a*JVS( 897 )
  W( 158 ) = W( 158 ) + a*JVS( 898 )
  W( 160 ) = W( 160 ) + a*JVS( 899 )
  W( 165 ) = W( 165 ) + a*JVS( 900 )
  W( 167 ) = W( 167 ) + a*JVS( 901 )
  a = -W( 152 ) / JVS( 912 )
  W( 152 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 913 )
  W( 157 ) = W( 157 ) + a*JVS( 914 )
  W( 160 ) = W( 160 ) + a*JVS( 915 )
  W( 165 ) = W( 165 ) + a*JVS( 916 )
  W( 167 ) = W( 167 ) + a*JVS( 917 )
  a = -W( 153 ) / JVS( 937 )
  W( 153 ) = -a
  W( 154 ) = W( 154 ) + a*JVS( 938 )
  W( 155 ) = W( 155 ) + a*JVS( 939 )
  W( 156 ) = W( 156 ) + a*JVS( 940 )
  W( 157 ) = W( 157 ) + a*JVS( 941 )
  W( 158 ) = W( 158 ) + a*JVS( 942 )
  W( 160 ) = W( 160 ) + a*JVS( 943 )
  W( 161 ) = W( 161 ) + a*JVS( 944 )
  W( 165 ) = W( 165 ) + a*JVS( 945 )
  W( 166 ) = W( 166 ) + a*JVS( 946 )
  W( 167 ) = W( 167 ) + a*JVS( 947 )
  a = -W( 154 ) / JVS( 959 )
  W( 154 ) = -a
  W( 155 ) = W( 155 ) + a*JVS( 960 )
  W( 156 ) = W( 156 ) + a*JVS( 961 )
  W( 157 ) = W( 157 ) + a*JVS( 962 )
  W( 158 ) = W( 158 ) + a*JVS( 963 )
  W( 160 ) = W( 160 ) + a*JVS( 964 )
  W( 165 ) = W( 165 ) + a*JVS( 965 )
  W( 166 ) = W( 166 ) + a*JVS( 966 )
  W( 167 ) = W( 167 ) + a*JVS( 967 )
  W( 168 ) = W( 168 ) + a*JVS( 968 )
  a = -W( 155 ) / JVS( 981 )
  W( 155 ) = -a
  W( 156 ) = W( 156 ) + a*JVS( 982 )
  W( 157 ) = W( 157 ) + a*JVS( 983 )
  W( 158 ) = W( 158 ) + a*JVS( 984 )
  W( 160 ) = W( 160 ) + a*JVS( 985 )
  W( 162 ) = W( 162 ) + a*JVS( 986 )
  W( 163 ) = W( 163 ) + a*JVS( 987 )
  W( 164 ) = W( 164 ) + a*JVS( 988 )
  W( 165 ) = W( 165 ) + a*JVS( 989 )
  W( 166 ) = W( 166 ) + a*JVS( 990 )
  W( 167 ) = W( 167 ) + a*JVS( 991 )
  W( 168 ) = W( 168 ) + a*JVS( 992 )
  a = -W( 156 ) / JVS( 1006 )
  W( 156 ) = -a
  W( 157 ) = W( 157 ) + a*JVS( 1007 )
  W( 159 ) = W( 159 ) + a*JVS( 1008 )
  W( 160 ) = W( 160 ) + a*JVS( 1009 )
  W( 161 ) = W( 161 ) + a*JVS( 1010 )
  W( 162 ) = W( 162 ) + a*JVS( 1011 )
  W( 163 ) = W( 163 ) + a*JVS( 1012 )
  W( 164 ) = W( 164 ) + a*JVS( 1013 )
  W( 165 ) = W( 165 ) + a*JVS( 1014 )
  W( 167 ) = W( 167 ) + a*JVS( 1015 )
  a = -W( 157 ) / JVS( 1045 )
  W( 157 ) = -a
  W( 158 ) = W( 158 ) + a*JVS( 1046 )
  W( 159 ) = W( 159 ) + a*JVS( 1047 )
  W( 160 ) = W( 160 ) + a*JVS( 1048 )
  W( 161 ) = W( 161 ) + a*JVS( 1049 )
  W( 162 ) = W( 162 ) + a*JVS( 1050 )
  W( 163 ) = W( 163 ) + a*JVS( 1051 )
  W( 164 ) = W( 164 ) + a*JVS( 1052 )
  W( 165 ) = W( 165 ) + a*JVS( 1053 )
  W( 166 ) = W( 166 ) + a*JVS( 1054 )
  W( 167 ) = W( 167 ) + a*JVS( 1055 )
  W( 168 ) = W( 168 ) + a*JVS( 1056 )
  a = -W( 158 ) / JVS( 1078 )
  W( 158 ) = -a
  W( 159 ) = W( 159 ) + a*JVS( 1079 )
  W( 160 ) = W( 160 ) + a*JVS( 1080 )
  W( 161 ) = W( 161 ) + a*JVS( 1081 )
  W( 162 ) = W( 162 ) + a*JVS( 1082 )
  W( 163 ) = W( 163 ) + a*JVS( 1083 )
  W( 164 ) = W( 164 ) + a*JVS( 1084 )
  W( 165 ) = W( 165 ) + a*JVS( 1085 )
  W( 166 ) = W( 166 ) + a*JVS( 1086 )
  W( 167 ) = W( 167 ) + a*JVS( 1087 )
  W( 168 ) = W( 168 ) + a*JVS( 1088 )
  a = -W( 159 ) / JVS( 1099 )
  W( 159 ) = -a
  W( 160 ) = W( 160 ) + a*JVS( 1100 )
  W( 161 ) = W( 161 ) + a*JVS( 1101 )
  W( 162 ) = W( 162 ) + a*JVS( 1102 )
  W( 163 ) = W( 163 ) + a*JVS( 1103 )
  W( 164 ) = W( 164 ) + a*JVS( 1104 )
  W( 165 ) = W( 165 ) + a*JVS( 1105 )
  W( 166 ) = W( 166 ) + a*JVS( 1106 )
  W( 167 ) = W( 167 ) + a*JVS( 1107 )
  W( 168 ) = W( 168 ) + a*JVS( 1108 )
  a = -W( 160 ) / JVS( 1125 )
  W( 160 ) = -a
  W( 161 ) = W( 161 ) + a*JVS( 1126 )
  W( 162 ) = W( 162 ) + a*JVS( 1127 )
  W( 163 ) = W( 163 ) + a*JVS( 1128 )
  W( 164 ) = W( 164 ) + a*JVS( 1129 )
  W( 165 ) = W( 165 ) + a*JVS( 1130 )
  W( 166 ) = W( 166 ) + a*JVS( 1131 )
  W( 167 ) = W( 167 ) + a*JVS( 1132 )
  W( 168 ) = W( 168 ) + a*JVS( 1133 )
  a = -W( 161 ) / JVS( 1177 )
  W( 161 ) = -a
  W( 162 ) = W( 162 ) + a*JVS( 1178 )
  W( 163 ) = W( 163 ) + a*JVS( 1179 )
  W( 164 ) = W( 164 ) + a*JVS( 1180 )
  W( 165 ) = W( 165 ) + a*JVS( 1181 )
  W( 166 ) = W( 166 ) + a*JVS( 1182 )
  W( 167 ) = W( 167 ) + a*JVS( 1183 )
  W( 168 ) = W( 168 ) + a*JVS( 1184 )
  a = -W( 162 ) / JVS( 1203 )
  W( 162 ) = -a
  W( 163 ) = W( 163 ) + a*JVS( 1204 )
  W( 164 ) = W( 164 ) + a*JVS( 1205 )
  W( 165 ) = W( 165 ) + a*JVS( 1206 )
  W( 166 ) = W( 166 ) + a*JVS( 1207 )
  W( 167 ) = W( 167 ) + a*JVS( 1208 )
  W( 168 ) = W( 168 ) + a*JVS( 1209 )
  a = -W( 163 ) / JVS( 1238 )
  W( 163 ) = -a
  W( 164 ) = W( 164 ) + a*JVS( 1239 )
  W( 165 ) = W( 165 ) + a*JVS( 1240 )
  W( 166 ) = W( 166 ) + a*JVS( 1241 )
  W( 167 ) = W( 167 ) + a*JVS( 1242 )
  W( 168 ) = W( 168 ) + a*JVS( 1243 )
  a = -W( 164 ) / JVS( 1256 )
  W( 164 ) = -a
  W( 165 ) = W( 165 ) + a*JVS( 1257 )
  W( 166 ) = W( 166 ) + a*JVS( 1258 )
  W( 167 ) = W( 167 ) + a*JVS( 1259 )
  W( 168 ) = W( 168 ) + a*JVS( 1260 )
  a = -W( 165 ) / JVS( 1351 )
  W( 165 ) = -a
  W( 166 ) = W( 166 ) + a*JVS( 1352 )
  W( 167 ) = W( 167 ) + a*JVS( 1353 )
  W( 168 ) = W( 168 ) + a*JVS( 1354 )
  a = -W( 166 ) / JVS( 1396 )
  W( 166 ) = -a
  W( 167 ) = W( 167 ) + a*JVS( 1397 )
  W( 168 ) = W( 168 ) + a*JVS( 1398 )
  a = -W( 167 ) / JVS( 1439 )
  W( 167 ) = -a
  W( 168 ) = W( 168 ) + a*JVS( 1440 )
  JVS( 1441) = W( 99 )
  JVS( 1442) = W( 114 )
  JVS( 1443) = W( 116 )
  JVS( 1444) = W( 125 )
  JVS( 1445) = W( 135 )
  JVS( 1446) = W( 136 )
  JVS( 1447) = W( 137 )
  JVS( 1448) = W( 139 )
  JVS( 1449) = W( 140 )
  JVS( 1450) = W( 143 )
  JVS( 1451) = W( 144 )
  JVS( 1452) = W( 147 )
  JVS( 1453) = W( 148 )
  JVS( 1454) = W( 149 )
  JVS( 1455) = W( 151 )
  JVS( 1456) = W( 152 )
  JVS( 1457) = W( 153 )
  JVS( 1458) = W( 154 )
  JVS( 1459) = W( 155 )
  JVS( 1460) = W( 156 )
  JVS( 1461) = W( 157 )
  JVS( 1462) = W( 158 )
  JVS( 1463) = W( 159 )
  JVS( 1464) = W( 160 )
  JVS( 1465) = W( 161 )
  JVS( 1466) = W( 162 )
  JVS( 1467) = W( 163 )
  JVS( 1468) = W( 164 )
  JVS( 1469) = W( 165 )
  JVS( 1470) = W( 166 )
  JVS( 1471) = W( 167 )
  JVS( 1472) = W( 168 )
   END SUBROUTINE decomp_saprc99_mosaic_4bin_vbs9
END MODULE saprc99_mosaic_4bin_vbs9_Integrator
