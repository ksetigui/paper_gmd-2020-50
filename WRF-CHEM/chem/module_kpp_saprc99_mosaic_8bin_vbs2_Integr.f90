MODULE saprc99_mosaic_8bin_vbs2_Integrator
 USE saprc99_mosaic_8bin_vbs2_Parameters
 USE saprc99_mosaic_8bin_vbs2_Precision
 USE saprc99_mosaic_8bin_vbs2_JacobianSP
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
SUBROUTINE saprc99_mosaic_8bin_vbs2_INTEGRATE( TIN, TOUT, &
  FIX, VAR, RCONST, ATOL, RTOL, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )
   USE saprc99_mosaic_8bin_vbs2_Parameters
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
   CALL saprc99_mosaic_8bin_vbs2_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT, &
         ATOL,RTOL, &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
   STEPMIN = RCNTRL(ihexit)
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U)) IERR_U = IERR
END SUBROUTINE saprc99_mosaic_8bin_vbs2_INTEGRATE
SUBROUTINE saprc99_mosaic_8bin_vbs2_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol, &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
  USE saprc99_mosaic_8bin_vbs2_Parameters
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
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF
   Roundoff = saprc99_mosaic_8bin_vbs2_WLAMCH('E')
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO
   SELECT CASE (Method)
     CASE (1)
       CALL saprc99_mosaic_8bin_vbs2_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL saprc99_mosaic_8bin_vbs2_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL saprc99_mosaic_8bin_vbs2_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL saprc99_mosaic_8bin_vbs2_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL saprc99_mosaic_8bin_vbs2_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT
   CALL saprc99_mosaic_8bin_vbs2_ros_Integrator(Y,Tstart,Tend,Texit, &
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
 SUBROUTINE saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(Code,T,H,IERR)
   USE saprc99_mosaic_8bin_vbs2_Precision
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
 END SUBROUTINE saprc99_mosaic_8bin_vbs2_ros_ErrorMsg
 SUBROUTINE saprc99_mosaic_8bin_vbs2_ros_Integrator (Y, Tstart, Tend, T, &
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
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN
      CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF
   Hexit = H
   H = MIN(H,ABS(Tend-T))
   CALL saprc99_mosaic_8bin_vbs2_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF (.NOT.Autonomous) THEN
      CALL saprc99_mosaic_8bin_vbs2_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF
   CALL saprc99_mosaic_8bin_vbs2_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)
UntilAccepted: DO
   CALL saprc99_mosaic_8bin_vbs2_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec, Nsng )
   IF (Singular) THEN
       CALL saprc99_mosaic_8bin_vbs2_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF
Stage: DO istage = 1, ros_S
       ioffset = NVAR*(istage-1)
       IF ( istage == 1 ) THEN
         CALL saprc99_mosaic_8bin_vbs2_WCOPY(NVAR,Fcn0,1,Fcn,1)
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL saprc99_mosaic_8bin_vbs2_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL saprc99_mosaic_8bin_vbs2_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL saprc99_mosaic_8bin_vbs2_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF
       CALL saprc99_mosaic_8bin_vbs2_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL saprc99_mosaic_8bin_vbs2_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL saprc99_mosaic_8bin_vbs2_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL saprc99_mosaic_8bin_vbs2_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)
   END DO Stage
   CALL saprc99_mosaic_8bin_vbs2_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL saprc99_mosaic_8bin_vbs2_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO
   CALL saprc99_mosaic_8bin_vbs2_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL saprc99_mosaic_8bin_vbs2_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = saprc99_mosaic_8bin_vbs2_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )
   Fac = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac
   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN
      Nacc = Nacc+1
      CALL saprc99_mosaic_8bin_vbs2_WCOPY(NVAR,Ynew,1,Y,1)
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
  END SUBROUTINE saprc99_mosaic_8bin_vbs2_ros_Integrator
  REAL(kind=dp) FUNCTION saprc99_mosaic_8bin_vbs2_ros_ErrorNorm ( Y, Ynew, Yerr, &
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
    saprc99_mosaic_8bin_vbs2_ros_ErrorNorm = Err
  END FUNCTION saprc99_mosaic_8bin_vbs2_ros_ErrorNorm
  SUBROUTINE saprc99_mosaic_8bin_vbs2_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   IMPLICIT NONE
   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)
   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)
   INTEGER, INTENT(INOUT) ::Nfun
   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp
   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL saprc99_mosaic_8bin_vbs2_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL saprc99_mosaic_8bin_vbs2_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL saprc99_mosaic_8bin_vbs2_WSCAL(NVAR,(ONE/Delta),dFdT,1)
  END SUBROUTINE saprc99_mosaic_8bin_vbs2_ros_FunTimeDeriv
  SUBROUTINE saprc99_mosaic_8bin_vbs2_ros_PrepareMatrix ( H, Direction, gam, &
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
     CALL saprc99_mosaic_8bin_vbs2_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL saprc99_mosaic_8bin_vbs2_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO
     CALL saprc99_mosaic_8bin_vbs2_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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
  END SUBROUTINE saprc99_mosaic_8bin_vbs2_ros_PrepareMatrix
  SUBROUTINE saprc99_mosaic_8bin_vbs2_ros_Decomp( A, Pivot, ising, Ndec )
   IMPLICIT NONE
   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)
   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec
CALL decomp_saprc99_mosaic_8bin_vbs2 ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1
  END SUBROUTINE saprc99_mosaic_8bin_vbs2_ros_Decomp
  SUBROUTINE saprc99_mosaic_8bin_vbs2_ros_Solve( A, Pivot, b, Nsol )
   IMPLICIT NONE
   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)
   INTEGER, INTENT(INOUT) :: nsol
   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)
   CALL saprc99_mosaic_8bin_vbs2_KppSolve( A, b )
   Nsol = Nsol+1
  END SUBROUTINE saprc99_mosaic_8bin_vbs2_ros_Solve
  SUBROUTINE saprc99_mosaic_8bin_vbs2_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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
 END SUBROUTINE saprc99_mosaic_8bin_vbs2_Ros2
  SUBROUTINE saprc99_mosaic_8bin_vbs2_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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
  END SUBROUTINE saprc99_mosaic_8bin_vbs2_Ros3
  SUBROUTINE saprc99_mosaic_8bin_vbs2_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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
  END SUBROUTINE saprc99_mosaic_8bin_vbs2_Ros4
  SUBROUTINE saprc99_mosaic_8bin_vbs2_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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
  END SUBROUTINE saprc99_mosaic_8bin_vbs2_Rodas3
  SUBROUTINE saprc99_mosaic_8bin_vbs2_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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
  END SUBROUTINE saprc99_mosaic_8bin_vbs2_Rodas4
END SUBROUTINE saprc99_mosaic_8bin_vbs2_Rosenbrock
SUBROUTINE saprc99_mosaic_8bin_vbs2_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )
   USE saprc99_mosaic_8bin_vbs2_Parameters
   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)
   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun
   CALL saprc99_mosaic_8bin_vbs2_Fun( Y, FIX, RCONST, Ydot )
   Nfun = Nfun+1
END SUBROUTINE saprc99_mosaic_8bin_vbs2_FunTemplate
SUBROUTINE saprc99_mosaic_8bin_vbs2_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )
 USE saprc99_mosaic_8bin_vbs2_Parameters
 USE saprc99_mosaic_8bin_vbs2_Jacobian
    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)
    INTEGER :: Njac
    REAL(kind=dp) :: Jcb(LU_NONZERO)
    REAL(kind=dp) :: Told
    CALL saprc99_mosaic_8bin_vbs2_Jac_SP( Y, FIX, RCONST, Jcb )
    Njac = Njac+1
END SUBROUTINE saprc99_mosaic_8bin_vbs2_JacTemplate
SUBROUTINE saprc99_mosaic_8bin_vbs2_Fun ( V, F, RCT, Vdot )
  REAL(kind=dp) :: V(NVAR)
  REAL(kind=dp) :: F(NFIX)
  REAL(kind=dp) :: RCT(NREACT)
  REAL(kind=dp) :: Vdot(NVAR)
  REAL(kind=dp) :: A(NREACT)
  A(1) = RCT(1)*V(101)
  A(2) = RCT(2)*V(88)*F(2)
  A(3) = RCT(3)*V(88)*V(92)
  A(4) = RCT(4)*V(88)*V(93)*F(2)
  A(5) = RCT(5)*V(88)*V(101)
  A(6) = RCT(6)*V(88)*V(101)
  A(7) = RCT(7)*V(92)*V(93)
  A(8) = RCT(8)*V(92)*V(101)
  A(9) = RCT(9)*V(93)*V(99)
  A(10) = RCT(10)*V(93)*V(93)*F(2)
  A(11) = RCT(11)*V(99)*V(101)
  A(12) = RCT(12)*V(46)
  A(13) = RCT(13)*V(46)*F(1)
  A(14) = RCT(14)*V(99)*V(101)
  A(15) = RCT(15)*V(99)
  A(16) = RCT(16)*V(99)
  A(17) = RCT(17)*V(92)
  A(18) = RCT(18)*V(92)
  A(19) = RCT(19)*V(30)*F(1)
  A(20) = RCT(20)*V(30)*F(2)
  A(21) = RCT(21)*V(93)*V(103)
  A(22) = RCT(22)*V(47)
  A(23) = RCT(23)*V(47)
  A(24) = RCT(24)*V(47)*V(103)
  A(25) = RCT(25)*V(101)*V(103)
  A(26) = RCT(26)*V(99)*V(103)
  A(27) = RCT(27)*V(70)*V(103)
  A(28) = RCT(28)*V(70)
  A(29) = RCT(29)*V(69)*V(103)
  A(30) = RCT(30)*V(92)*V(103)
  A(31) = RCT(31)*V(93)*V(96)
  A(32) = RCT(32)*V(96)*V(101)
  A(33) = RCT(33)*V(55)
  A(34) = RCT(34)*V(55)
  A(35) = RCT(35)*V(55)*V(103)
  A(36) = RCT(36)*V(92)*V(96)
  A(37) = RCT(37)*V(96)*V(96)
  A(38) = RCT(38)*V(96)*V(96)*F(1)
  A(39) = RCT(39)*V(96)*V(99)
  A(40) = RCT(40)*V(99)*V(99)
  A(41) = RCT(41)*V(42)
  A(42) = RCT(42)*V(42)*V(103)
  A(43) = RCT(43)*V(96)*V(103)
  A(44) = RCT(44)*V(36)*V(103)
  A(45) = RCT(45)*V(103)*F(2)
  A(46) = RCT(46)*V(93)*V(102)
  A(47) = RCT(47)*V(96)*V(102)
  A(48) = RCT(48)*V(99)*V(102)
  A(49) = RCT(49)*V(102)*V(102)
  A(50) = RCT(50)*V(102)*V(102)
  A(51) = RCT(51)*V(93)*V(94)
  A(52) = RCT(52)*V(94)*V(96)
  A(53) = RCT(53)*V(94)*V(99)
  A(54) = RCT(54)*V(94)*V(102)
  A(55) = RCT(55)*V(94)*V(94)
  A(56) = RCT(56)*V(77)*V(93)
  A(57) = RCT(57)*V(77)*V(96)
  A(58) = RCT(58)*V(77)*V(99)
  A(59) = RCT(59)*V(77)*V(102)
  A(60) = RCT(60)*V(77)*V(94)
  A(61) = RCT(61)*V(77)*V(77)
  A(62) = RCT(62)*V(93)*V(104)
  A(63) = RCT(63)*V(96)*V(104)
  A(64) = RCT(64)*V(102)*V(104)
  A(65) = RCT(65)*V(99)*V(104)
  A(66) = RCT(66)*V(94)*V(104)
  A(67) = RCT(67)*V(77)*V(104)
  A(68) = RCT(68)*V(104)*V(104)
  A(69) = RCT(69)*V(95)*V(101)
  A(70) = RCT(70)*V(38)
  A(71) = RCT(71)*V(93)*V(95)
  A(72) = RCT(72)*V(95)*V(96)
  A(73) = RCT(73)*V(95)*V(99)
  A(74) = RCT(74)*V(95)*V(102)
  A(75) = RCT(75)*V(94)*V(95)
  A(76) = RCT(76)*V(77)*V(95)
  A(77) = RCT(77)*V(95)*V(104)
  A(78) = RCT(78)*V(95)*V(95)
  A(79) = RCT(79)*V(97)*V(101)
  A(80) = RCT(80)*V(39)
  A(81) = RCT(81)*V(93)*V(97)
  A(82) = RCT(82)*V(96)*V(97)
  A(83) = RCT(83)*V(97)*V(99)
  A(84) = RCT(84)*V(97)*V(102)
  A(85) = RCT(85)*V(94)*V(97)
  A(86) = RCT(86)*V(77)*V(97)
  A(87) = RCT(87)*V(97)*V(104)
  A(88) = RCT(88)*V(95)*V(97)
  A(89) = RCT(89)*V(97)*V(97)
  A(90) = RCT(90)*V(100)*V(101)
  A(91) = RCT(91)*V(40)
  A(92) = RCT(92)*V(93)*V(100)
  A(93) = RCT(93)*V(96)*V(100)
  A(94) = RCT(94)*V(99)*V(100)
  A(95) = RCT(95)*V(100)*V(102)
  A(96) = RCT(96)*V(94)*V(100)
  A(97) = RCT(97)*V(77)*V(100)
  A(98) = RCT(98)*V(100)*V(104)
  A(99) = RCT(99)*V(95)*V(100)
  A(100) = RCT(100)*V(97)*V(100)
  A(101) = RCT(101)*V(100)*V(100)
  A(102) = RCT(102)*V(98)*V(101)
  A(103) = RCT(103)*V(41)
  A(104) = RCT(104)*V(93)*V(98)
  A(105) = RCT(105)*V(96)*V(98)
  A(106) = RCT(106)*V(98)*V(99)
  A(107) = RCT(107)*V(98)*V(102)
  A(108) = RCT(108)*V(94)*V(98)
  A(109) = RCT(109)*V(77)*V(98)
  A(110) = RCT(110)*V(98)*V(104)
  A(111) = RCT(111)*V(95)*V(98)
  A(112) = RCT(112)*V(97)*V(98)
  A(113) = RCT(113)*V(98)*V(100)
  A(114) = RCT(114)*V(98)*V(98)
  A(115) = RCT(115)*V(49)*V(101)
  A(116) = RCT(116)*V(49)
  A(117) = RCT(117)*V(74)*V(101)
  A(118) = RCT(118)*V(74)*V(96)
  A(119) = RCT(119)*V(74)
  A(120) = RCT(120)*V(54)*V(101)
  A(121) = RCT(121)*V(54)*V(96)
  A(122) = RCT(122)*V(54)
  A(123) = RCT(123)*V(86)
  A(124) = RCT(124)*V(86)
  A(125) = RCT(125)*V(86)*V(103)
  A(126) = RCT(126)*V(86)*V(96)
  A(127) = RCT(127)*V(53)
  A(128) = RCT(128)*V(53)*V(93)
  A(129) = RCT(129)*V(86)*V(99)
  A(130) = RCT(130)*V(83)*V(103)
  A(131) = RCT(131)*V(83)
  A(132) = RCT(132)*V(83)*V(99)
  A(133) = RCT(133)*V(89)*V(103)
  A(134) = RCT(134)*V(89)
  A(135) = RCT(135)*V(89)*V(99)
  A(136) = RCT(136)*V(72)*V(103)
  A(137) = RCT(137)*V(72)
  A(138) = RCT(138)*V(90)*V(103)
  A(139) = RCT(139)*V(90)
  A(140) = RCT(140)*V(56)*V(103)
  A(141) = RCT(141)*V(44)*V(103)
  A(142) = RCT(142)*V(52)*V(103)
  A(143) = RCT(143)*V(52)
  A(144) = RCT(144)*V(65)*V(103)
  A(145) = RCT(145)*V(65)
  A(146) = RCT(146)*V(81)
  A(147) = RCT(147)*V(81)
  A(148) = RCT(148)*V(81)*V(103)
  A(149) = RCT(149)*V(81)*V(99)
  A(150) = RCT(150)*V(68)
  A(151) = RCT(151)*V(68)*V(103)
  A(152) = RCT(152)*V(68)*V(99)
  A(153) = RCT(153)*V(45)
  A(154) = RCT(154)*V(67)*V(103)
  A(155) = RCT(155)*V(67)*V(99)
  A(156) = RCT(156)*V(60)*V(103)
  A(157) = RCT(157)*V(60)*V(99)
  A(158) = RCT(158)*V(64)*V(99)
  A(159) = RCT(159)*V(66)*V(103)
  A(160) = RCT(160)*V(66)
  A(161) = RCT(161)*V(66)*V(99)
  A(162) = RCT(162)*V(78)*V(103)
  A(163) = RCT(163)*V(78)*V(92)
  A(164) = RCT(164)*V(78)*V(99)
  A(165) = RCT(165)*V(78)*V(88)
  A(166) = RCT(166)*V(78)
  A(167) = RCT(167)*V(85)*V(103)
  A(168) = RCT(168)*V(85)*V(92)
  A(169) = RCT(169)*V(85)*V(88)
  A(170) = RCT(170)*V(85)
  A(171) = RCT(171)*V(82)*V(103)
  A(172) = RCT(172)*V(82)*V(92)
  A(173) = RCT(173)*V(82)*V(99)
  A(174) = RCT(174)*V(82)
  A(175) = RCT(175)*V(91)*V(103)
  A(176) = RCT(176)*V(91)
  A(177) = RCT(177)*V(87)*V(103)
  A(178) = RCT(178)*V(87)
  A(179) = RCT(179)*V(62)*V(103)
  A(180) = RCT(180)*V(62)*V(92)
  A(181) = RCT(181)*V(58)*V(103)
  A(182) = RCT(182)*V(58)
  A(183) = RCT(183)*V(59)*V(103)
  A(184) = RCT(184)*V(59)
  A(185) = RCT(185)*V(35)*V(103)
  A(186) = RCT(186)*V(71)*V(103)
  A(187) = RCT(187)*V(71)*V(92)
  A(188) = RCT(188)*V(71)*V(99)
  A(189) = RCT(189)*V(71)*V(88)
  A(190) = RCT(190)*V(75)*V(103)
  A(191) = RCT(191)*V(75)*V(92)
  A(192) = RCT(192)*V(75)*V(99)
  A(193) = RCT(193)*V(75)*V(88)
  A(194) = RCT(194)*V(80)*V(103)
  A(195) = RCT(195)*V(80)*V(92)
  A(196) = RCT(196)*V(80)*V(99)
  A(197) = RCT(197)*V(80)*V(88)
  A(198) = RCT(198)*V(79)*V(103)
  A(199) = RCT(199)*V(79)*V(92)
  A(200) = RCT(200)*V(79)*V(99)
  A(201) = RCT(201)*V(79)*V(88)
  A(202) = RCT(202)*V(37)*V(103)
  A(203) = RCT(203)*V(43)*V(103)
  A(204) = RCT(204)*V(63)*V(103)
  A(205) = RCT(205)*V(48)*V(103)
  A(206) = RCT(206)*V(61)*V(103)
  A(207) = RCT(207)*V(50)*V(103)
  A(208) = RCT(208)*V(57)*V(103)
  A(209) = RCT(209)*V(51)*V(103)
  A(210) = RCT(210)*V(76)*V(103)
  A(211) = RCT(211)*V(76)*V(92)
  A(212) = RCT(212)*V(76)*V(99)
  A(213) = RCT(213)*V(76)*V(88)
  A(214) = RCT(214)*V(84)*V(103)
  A(215) = RCT(215)*V(84)*V(92)
  A(216) = RCT(216)*V(84)*V(99)
  A(217) = RCT(217)*V(84)*V(88)
  A(218) = RCT(218)*V(63)*V(92)
  A(219) = RCT(219)*V(73)*V(103)
  A(220) = RCT(220)*V(73)*V(92)
  A(221) = RCT(221)*V(73)*V(99)
  A(222) = RCT(222)*V(73)*V(88)
  A(223) = RCT(223)*V(36)
  A(224) = RCT(224)*V(96)
  A(225) = RCT(225)*V(36)
  A(226) = RCT(226)*V(1)
  A(227) = RCT(227)*V(70)
  A(228) = RCT(228)*V(42)
  A(229) = RCT(229)*V(2)
  A(230) = RCT(230)*V(57)*V(103)
  A(231) = RCT(231)*V(57)*V(103)
  A(232) = RCT(232)*V(57)*V(103)
  A(233) = RCT(233)*V(57)*V(103)
  A(234) = RCT(234)*V(51)*V(103)
  A(235) = RCT(235)*V(51)*V(103)
  A(236) = RCT(236)*V(51)*V(103)
  A(237) = RCT(237)*V(51)*V(103)
  A(238) = RCT(238)*V(75)*V(103)
  A(239) = RCT(239)*V(75)*V(103)
  A(240) = RCT(240)*V(75)*V(103)
  A(241) = RCT(241)*V(80)*V(103)
  A(242) = RCT(242)*V(80)*V(103)
  A(243) = RCT(243)*V(80)*V(103)
  A(244) = RCT(244)*V(79)*V(103)
  A(245) = RCT(245)*V(79)*V(103)
  A(246) = RCT(246)*V(75)*V(92)
  A(247) = RCT(247)*V(80)*V(92)
  A(248) = RCT(248)*V(80)*V(92)
  A(249) = RCT(249)*V(80)*V(92)
  A(250) = RCT(250)*V(79)*V(92)
  A(251) = RCT(251)*V(79)*V(92)
  A(252) = RCT(252)*V(79)*V(92)
  A(253) = RCT(253)*V(75)*V(99)
  A(254) = RCT(254)*V(75)*V(99)
  A(255) = RCT(255)*V(75)*V(99)
  A(256) = RCT(256)*V(80)*V(99)
  A(257) = RCT(257)*V(80)*V(99)
  A(258) = RCT(258)*V(79)*V(99)
  A(259) = RCT(259)*V(79)*V(99)
  A(260) = RCT(260)*V(3)*V(103)
  A(261) = RCT(261)*V(4)*V(103)
  A(262) = RCT(262)*V(32)*V(103)
  A(263) = RCT(263)*V(26)*V(103)
  A(264) = RCT(264)*V(5)*V(103)
  A(265) = RCT(265)*V(6)*V(103)
  A(266) = RCT(266)*V(34)*V(103)
  A(267) = RCT(267)*V(28)*V(103)
  A(268) = RCT(268)*V(31)*V(103)
  A(269) = RCT(269)*V(27)*V(103)
  A(270) = RCT(270)*V(33)*V(103)
  A(271) = RCT(271)*V(29)*V(103)
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
  Vdot(16) = A(230)+A(234)
  Vdot(17) = A(231)+A(235)
  Vdot(18) = A(232)+A(236)
  Vdot(19) = A(233)+A(237)
  Vdot(20) = A(238)
  Vdot(21) = A(239)+A(241)+A(244)+A(247)+A(250)+A(253)
  Vdot(22) = A(240)+A(242)+A(245)+A(246)+A(248)+A(251)+A(254)+A(256)+A(258)
  Vdot(23) = A(243)+A(249)+A(252)+A(255)+A(257)+A(259)
  Vdot(24) = A(21)+A(24)+A(25)+A(26)+A(27)+A(29)+A(30)+A(43)+A(44)+A(45)+A(125)+A(130)+A(133)+A(136)+A(138)+A(140)&
               &+A(141)+A(142)+A(144)+A(148)+A(151)+A(154)+A(156)+A(159)+A(162)+A(167)+A(171)+A(175)+A(177)+A(179)+A(181)&
               &+A(183)+A(185)+A(186)+A(190)+A(194)+A(198)+A(202)+A(203)+A(204)+A(205)+A(206)+A(207)+A(208)+A(209)+A(210)&
               &+A(214)+A(219)
  Vdot(25) = A(260)+A(261)+A(262)+A(263)+A(264)+A(265)+A(266)+A(267)+A(268)+A(269)+A(270)+A(271)
  Vdot(26) = 0.5*A(268)+A(269)
  Vdot(27) = -A(269)
  Vdot(28) = 0.5*A(270)+A(271)
  Vdot(29) = -A(271)
  Vdot(30) = A(18)-A(19)-A(20)
  Vdot(31) = -A(268)
  Vdot(32) = A(268)
  Vdot(33) = -A(270)
  Vdot(34) = A(270)
  Vdot(35) = -A(185)
  Vdot(36) = -A(44)-A(223)-A(225)
  Vdot(37) = -A(202)
  Vdot(38) = A(69)-A(70)
  Vdot(39) = A(79)-A(80)
  Vdot(40) = A(90)-A(91)
  Vdot(41) = A(102)-A(103)
  Vdot(42) = A(37)+A(38)-A(41)-A(42)-A(228)
  Vdot(43) = -A(203)
  Vdot(44) = -A(141)
  Vdot(45) = -A(153)+0.031*A(195)+0.031*A(199)+0.087*A(209)
  Vdot(46) = A(11)-A(12)-A(13)
  Vdot(47) = A(21)-A(22)-A(23)-A(24)
  Vdot(48) = -A(205)
  Vdot(49) = -A(115)-A(116)+0.236*A(205)
  Vdot(50) = -A(207)
  Vdot(51) = -A(209)
  Vdot(52) = A(47)-A(142)-A(143)
  Vdot(53) = A(126)-A(127)-A(128)
  Vdot(54) = -A(120)-A(121)-A(122)+A(158)
  Vdot(55) = A(32)-A(33)-A(34)-A(35)
  Vdot(56) = A(49)+0.25*A(54)+0.25*A(64)-A(140)
  Vdot(57) = -A(208)
  Vdot(58) = -A(181)-A(182)+0.108*A(208)+0.099*A(209)
  Vdot(59) = -A(183)-A(184)+0.051*A(208)+0.093*A(209)
  Vdot(60) = -A(156)-A(157)+0.207*A(208)+0.187*A(209)
  Vdot(61) = -A(206)
  Vdot(62) = -A(179)-A(180)+0.491*A(208)+0.561*A(209)
  Vdot(63) = -A(204)-A(218)
  Vdot(64) = A(117)+A(121)+A(122)-A(158)
  Vdot(65) = A(52)+A(63)-A(144)-A(145)
  Vdot(66) = -A(159)-A(160)-A(161)+0.059*A(208)+0.05*A(209)+0.061*A(214)+0.042*A(215)+0.015*A(216)
  Vdot(67) = A(118)+A(119)-A(154)-A(155)+0.017*A(208)
  Vdot(68) = -A(150)-A(151)-A(152)+0.23*A(156)+0.084*A(162)+0.9*A(163)+0.3*A(167)+0.95*A(168)+0.174*A(171)+0.742*A(172)&
               &+0.008*A(173)+0.5*A(182)+0.5*A(184)+0.119*A(208)+0.287*A(209)
  Vdot(69) = -A(29)+A(123)+A(124)+A(125)+A(129)+A(131)+0.034*A(133)+A(134)+2*A(146)+A(147)+1.26*A(148)+1.26*A(149)&
               &+A(150)+A(151)+A(152)+0.416*A(162)+0.45*A(163)+0.5*A(164)+0.67*A(166)+0.475*A(168)+0.7*A(170)+0.336*A(171)&
               &+0.498*A(172)+0.572*A(173)+1.233*A(174)+A(179)+1.5*A(180)+A(182)+A(184)+0.5*A(187)+0.491*A(189)+0.275*A(191)&
               &+0.157*A(195)+0.157*A(199)+0.393*A(204)+0.002*A(206)+0.345*A(211)+0.265*A(215)+0.012*A(217)+1.5*A(218)+0.51&
               &*A(220)
  Vdot(70) = 2*A(13)+A(25)-A(27)-A(28)+0.2*A(39)+A(129)+A(132)+A(135)+A(149)+A(152)+A(155)+A(157)+A(158)+A(161)+0.5&
               &*A(164)+0.15*A(173)-A(227)
  Vdot(71) = -A(186)-A(187)-A(188)-A(189)
  Vdot(72) = A(116)-A(136)-A(137)+0.006*A(177)+0.02*A(178)+0.13*A(195)+0.13*A(199)+0.704*A(203)+0.024*A(205)+0.452&
               &*A(206)+0.072*A(207)+0.005*A(210)+0.001*A(211)+0.024*A(212)+0.127*A(214)+0.045*A(215)+0.102*A(216)
  Vdot(73) = -A(219)-A(220)-A(221)-A(222)
  Vdot(74) = A(92)+A(94)+A(99)+A(100)+2*A(101)+A(113)-A(117)-A(118)-A(119)+0.24*A(154)+A(155)+0.24*A(156)+A(157)
  Vdot(75) = -A(190)-A(191)-A(192)-A(193)
  Vdot(76) = -A(210)-A(211)-A(212)-A(213)
  Vdot(77) = -A(56)-A(57)-A(58)-A(59)-A(60)-A(67)-A(76)-A(86)+A(92)+A(94)-A(97)+A(99)+A(100)+2*A(101)-A(109)+A(113)&
               &+A(136)+0.616*A(138)+0.675*A(167)+0.515*A(176)+0.596*A(177)+0.152*A(178)+A(181)+A(182)+A(183)+A(184)+0.079&
               &*A(190)+0.126*A(191)+0.187*A(192)+0.24*A(193)+0.5*A(194)+0.729*A(195)+0.75*A(196)+0.5*A(198)+0.729*A(199)&
               &+0.75*A(200)+0.559*A(205)+0.936*A(206)+0.948*A(207)+0.205*A(210)+0.488*A(212)+0.001*A(214)+0.137*A(215)&
               &+0.711*A(216)
  Vdot(78) = -A(162)-A(163)-A(164)-A(165)-A(166)+0.23*A(190)+0.39*A(191)+0.025*A(214)+0.026*A(215)+0.012*A(217)
  Vdot(79) = -A(198)-A(199)-A(200)-A(201)
  Vdot(80) = -A(194)-A(195)-A(196)-A(197)
  Vdot(81) = -A(146)-A(147)-A(148)-A(149)+0.23*A(154)+0.15*A(171)+0.023*A(172)+A(180)+0.5*A(182)+0.5*A(184)+0.009*A(189)&
               &+0.001*A(195)+0.001*A(199)+0.607*A(204)+0.118*A(208)+0.097*A(209)
  Vdot(82) = -A(171)-A(172)-A(173)-A(174)+0.357*A(190)+0.936*A(192)+0.025*A(214)
  Vdot(83) = A(81)+A(83)+A(88)+2*A(89)+A(100)+A(112)-A(130)-A(131)-A(132)+0.034*A(133)+A(134)+0.482*A(138)+A(139)+0.96&
               &*A(141)+0.129*A(171)+0.047*A(172)+0.467*A(174)+0.084*A(175)+0.246*A(176)+0.439*A(177)+0.431*A(178)+0.195&
               &*A(186)+0.25*A(189)+A(202)+0.445*A(205)+0.455*A(206)+0.099*A(207)+0.294*A(210)+0.154*A(211)+0.009*A(212)&
               &+0.732*A(214)+0.456*A(215)+0.507*A(216)+0.984*A(219)+0.5*A(220)
  Vdot(84) = -A(214)-A(215)-A(216)-A(217)
  Vdot(85) = -A(167)-A(168)-A(169)-A(170)+0.32*A(190)+0.16*A(191)+0.019*A(215)+0.048*A(216)
  Vdot(86) = A(46)+A(48)+A(49)+2*A(50)+0.75*A(54)+0.75*A(64)+A(74)+A(84)+A(95)+A(104)+A(106)+A(107)+A(111)+A(112)+A(113)&
               &+2*A(114)-A(123)-A(124)-A(125)-A(126)+A(127)-A(129)+A(136)+0.115*A(138)+A(140)+0.081*A(141)+0.35*A(142)&
               &+A(143)+A(147)+0.084*A(162)+0.2*A(163)+0.67*A(166)+0.3*A(167)+0.1*A(168)+0.055*A(171)+0.125*A(172)+0.227&
               &*A(173)+0.3*A(174)+0.213*A(175)+0.506*A(176)+0.01*A(177)+0.134*A(178)+1.61*A(186)+A(187)+0.191*A(189)+0.624&
               &*A(190)+0.592*A(191)+0.24*A(193)+0.276*A(194)+0.235*A(195)+0.276*A(198)+0.235*A(199)+0.096*A(204)+0.026&
               &*A(205)+0.024*A(206)+0.026*A(207)+0.732*A(210)+0.5*A(211)+0.244*A(214)+0.269*A(215)+0.079*A(216)+0.984&
               &*A(219)+0.5*A(220)
  Vdot(87) = A(62)+A(115)+0.572*A(173)-0.69*A(177)-A(178)+0.276*A(196)+0.276*A(200)+0.511*A(212)+0.321*A(216)
  Vdot(88) = A(1)-A(2)-A(3)-A(4)-A(5)-A(6)+A(16)+A(17)+A(20)-A(165)-A(169)-A(189)-A(193)-A(197)-A(201)-A(213)-A(217)&
               &-A(222)
  Vdot(89) = -A(133)-A(134)-A(135)+0.37*A(138)+A(144)+A(145)+A(165)+0.675*A(167)+0.45*A(169)+0.013*A(171)+0.218*A(173)&
               &+0.558*A(175)+0.71*A(176)+0.213*A(177)+0.147*A(178)+A(179)+A(181)+A(183)+A(188)+0.474*A(194)+0.205*A(195)&
               &+0.474*A(196)+0.147*A(197)+0.474*A(198)+0.205*A(199)+0.474*A(200)+0.147*A(201)+0.261*A(203)+0.122*A(205)&
               &+0.244*A(206)+0.204*A(207)+0.497*A(210)+0.363*A(211)+0.037*A(212)+0.45*A(213)+0.511*A(214)+0.305*A(215)&
               &+0.151*A(216)+0.069*A(217)+0.45*A(222)
  Vdot(90) = 0.5*A(64)+A(65)+0.5*A(66)+A(68)-A(138)-A(139)+0.416*A(162)+0.55*A(169)+0.15*A(171)+0.21*A(172)+0.233*A(174)&
               &+0.115*A(175)+0.177*A(177)+0.243*A(178)+0.332*A(205)+0.11*A(206)+0.089*A(207)+0.437*A(213)+0.072*A(214)&
               &+0.026*A(215)+0.001*A(216)+0.659*A(217)+0.55*A(222)
  Vdot(91) = 0.5*A(64)+0.5*A(66)+A(68)+A(77)+A(87)+A(98)+0.7*A(170)+0.332*A(171)-0.671*A(175)-A(176)+0.048*A(177)+0.435&
               &*A(178)+0.1*A(191)+0.75*A(193)+0.276*A(194)+0.276*A(195)+0.853*A(197)+0.276*A(198)+0.276*A(199)+0.853*A(201)&
               &+0.125*A(206)+0.417*A(207)+0.055*A(208)+0.119*A(210)+0.215*A(211)+0.113*A(213)+0.043*A(215)+0.259*A(217)
  Vdot(92) = A(2)-A(3)-A(7)-A(8)-A(17)-A(18)-A(30)-A(36)+0.25*A(72)+0.25*A(82)+0.25*A(93)+0.25*A(105)-A(163)-A(168)&
               &-A(172)-A(180)-A(187)-A(191)-A(195)-A(199)-A(211)-A(215)-A(218)-A(220)
  Vdot(93) = A(1)-A(4)+A(5)-A(7)-A(9)-2*A(10)+A(14)+A(15)-A(21)+A(22)-A(31)-A(46)-A(51)-A(56)-A(62)-A(71)-A(81)-A(92)&
               &-A(104)-A(128)
  Vdot(94) = -A(51)-A(52)-A(53)-A(54)-2*A(55)-A(66)-A(75)+A(81)+A(83)-A(85)+A(88)+2*A(89)-A(96)+A(100)-A(108)+A(112)&
               &+0.034*A(133)+A(134)+0.37*A(138)+A(139)+0.05*A(141)+0.34*A(144)+0.76*A(154)+0.76*A(156)+0.5*A(162)+0.1&
               &*A(163)+0.5*A(164)+0.33*A(166)+0.3*A(167)+0.05*A(168)+0.67*A(171)+0.048*A(172)+0.799*A(173)+0.473*A(175)&
               &+0.96*A(176)+0.376*A(177)+0.564*A(178)+A(179)+A(182)+A(184)+A(186)+A(188)+0.2*A(189)+0.907*A(190)+0.066&
               &*A(191)+0.749*A(192)+0.75*A(194)+0.031*A(195)+0.276*A(196)+0.75*A(198)+0.031*A(199)+0.276*A(200)+A(202)&
               &+0.965*A(203)+0.1*A(204)+0.695*A(205)+0.835*A(206)+0.653*A(207)+0.765*A(208)+0.804*A(209)+0.91*A(210)+0.022&
               &*A(211)+0.824*A(212)+0.918*A(214)+0.033*A(215)+0.442*A(216)+0.012*A(217)+0.984*A(219)+0.949*A(221)
  Vdot(95) = -A(69)+A(70)-A(71)-A(72)-A(73)-A(74)-A(75)-A(77)-2*A(78)-A(88)-A(99)+A(104)+A(106)+A(112)+A(113)+2*A(114)&
               &+A(130)+A(132)+A(136)+A(137)+0.492*A(138)+A(139)+A(150)+A(151)+A(152)+2*A(153)+0.67*A(166)+0.675*A(167)&
               &+0.467*A(174)+0.029*A(175)+0.667*A(176)+A(181)+0.5*A(182)+A(183)+0.5*A(184)+0.123*A(195)+0.123*A(199)+0.011&
               &*A(206)+0.137*A(215)
  Vdot(96) = A(23)+A(26)+A(29)+A(30)-A(31)-A(32)+A(33)+0.61*A(34)-A(36)-2*A(37)-2*A(38)-A(39)+A(42)-A(43)+A(44)+A(45)&
               &+A(46)-A(47)+A(48)+2*A(50)+A(51)-A(52)+A(53)+A(54)+A(55)-A(63)+A(64)+A(65)+A(66)+A(68)-A(72)-A(82)-A(93)&
               &-A(105)-A(118)-A(121)+2*A(123)+A(125)-A(126)+A(127)+A(128)+A(129)+A(131)+A(134)+A(140)+0.95*A(141)+A(143)&
               &+A(145)+2*A(146)+0.63*A(148)+0.63*A(149)+A(150)+0.008*A(163)+0.34*A(166)+0.064*A(168)+0.4*A(172)+1.233&
               &*A(174)+0.379*A(175)+0.113*A(177)+0.341*A(178)+1.5*A(180)+0.5*A(182)+0.5*A(184)+0.12*A(187)+0.5*A(189)+0.033&
               &*A(195)+0.033*A(199)+0.297*A(204)+0.224*A(208)+0.187*A(209)+0.056*A(211)+0.003*A(215)+0.013*A(217)+1.5&
               &*A(218)+0.06*A(220)-A(224)
  Vdot(97) = -A(79)+A(80)-A(81)-A(82)-A(83)-A(84)-A(85)-A(87)-A(88)-2*A(89)-A(100)-A(112)+0.965*A(133)+A(135)+0.096&
               &*A(138)+0.37*A(148)+0.37*A(149)+0.1*A(163)+0.05*A(168)+0.048*A(172)+0.3*A(174)+0.049*A(175)+0.333*A(176)&
               &+0.201*A(195)+0.201*A(199)+0.006*A(215)
  Vdot(98) = -A(102)+A(103)-A(104)-A(105)-A(106)-A(107)-A(108)-A(110)-A(111)-A(112)-A(113)-2*A(114)+0.5*A(162)+0.5&
               &*A(164)+0.33*A(166)+0.3*A(170)+0.289*A(171)+0.15*A(173)+0.192*A(191)+0.24*A(193)
  Vdot(99) = A(6)+A(8)-A(9)-A(11)+A(12)-A(14)-A(15)-A(16)-A(26)+A(27)+0.39*A(34)-A(39)-2*A(40)-A(48)-A(53)-A(58)-A(65)&
               &-A(73)-A(83)-A(94)-A(106)-A(129)-A(132)-A(135)-A(149)-A(152)-A(155)-A(157)-A(158)-A(161)-A(164)-A(173)&
               &-A(188)-A(192)-A(196)-A(200)-A(212)-A(216)-A(221)
  Vdot(100) = -A(90)+A(91)-A(92)-A(93)-A(94)-A(95)-A(96)-A(98)-A(99)-A(100)-2*A(101)-A(113)+A(159)+A(161)
  Vdot(101) = -A(1)+A(4)-A(5)-A(6)+A(7)-A(8)+2*A(9)+2*A(10)-A(11)+A(12)+A(16)+A(23)+A(24)-A(25)+A(26)+A(28)+A(31)-A(32)&
                &+A(33)+0.61*A(34)+A(35)+0.8*A(39)+2*A(40)+A(46)+A(48)+A(51)+A(53)+A(56)+A(58)+A(65)-A(69)+A(70)+A(71)+A(73)&
                &-A(79)+A(80)+A(81)+A(83)-A(90)+A(91)+A(92)+A(94)-A(102)+A(103)+A(104)+A(106)-A(115)-A(117)-A(120)+A(128)&
                &+0.338*A(177)+A(178)+0.187*A(192)+0.474*A(196)+0.474*A(200)+0.391*A(216)
  Vdot(102) = -A(46)-A(47)-A(48)-2*A(49)-2*A(50)-A(54)-A(64)+A(71)+A(73)-A(74)+2*A(78)-A(84)+A(88)-A(95)+A(99)-A(107)&
                &+A(111)+A(116)+A(131)+A(137)+0.65*A(142)+0.3*A(170)+A(185)+0.3*A(189)+0.25*A(193)+0.011*A(206)+0.076*A(211)&
                &+0.197*A(215)+0.03*A(216)+0.26*A(220)
  Vdot(103) = 2*A(19)-A(21)+A(22)-A(24)-A(25)-A(26)-A(27)+A(28)-A(29)-A(30)+A(31)+0.39*A(34)-A(35)+A(36)+0.8*A(39)+2&
                &*A(41)-A(42)-A(43)-A(44)-A(45)-A(125)-A(130)-A(133)-A(136)-A(138)-A(140)-A(141)-0.65*A(142)+A(143)-0.34&
                &*A(144)+A(145)-A(148)-A(151)-A(154)-A(156)-A(159)-A(162)+0.208*A(163)+0.33*A(166)-A(167)+0.164*A(168)&
                &-A(171)+0.285*A(172)-A(175)-A(177)-A(179)+0.5*A(180)-A(181)-A(183)-A(185)-A(186)+0.12*A(187)-A(190)+0.266&
                &*A(191)-A(194)+0.567*A(195)-A(198)+0.567*A(199)-A(202)-A(203)-0.397*A(204)-A(205)-A(206)-A(207)-A(208)&
                &-A(209)-A(210)+0.155*A(211)-A(214)+0.378*A(215)+0.5*A(218)-A(219)+0.32*A(220)-A(260)-A(262)-A(264)-A(266)&
                &-A(268)-A(270)
  Vdot(104) = -A(62)-A(63)-A(64)-A(65)-A(66)-2*A(68)-A(77)-A(87)-A(98)-A(110)+0.001*A(133)+0.042*A(138)+0.025*A(167)&
                &+0.041*A(171)+0.051*A(173)+0.07*A(175)+0.04*A(176)+0.173*A(177)+0.095*A(178)+0.093*A(190)+0.008*A(191)&
                &+0.064*A(192)+0.01*A(193)+0.25*A(194)+0.18*A(195)+0.25*A(196)+0.25*A(198)+0.18*A(199)+0.25*A(200)+0.035&
                &*A(203)+0.07*A(205)+0.143*A(206)+0.347*A(207)+0.011*A(208)+0.009*A(209)+0.09*A(210)+0.001*A(211)+0.176&
                &*A(212)+0.082*A(214)+0.002*A(215)+0.136*A(216)+0.001*A(217)+0.016*A(219)+0.051*A(221)
END SUBROUTINE saprc99_mosaic_8bin_vbs2_Fun
SUBROUTINE saprc99_mosaic_8bin_vbs2_Jac_SP ( V, F, RCT, JVS )
  REAL(kind=dp) :: V(NVAR)
  REAL(kind=dp) :: F(NFIX)
  REAL(kind=dp) :: RCT(NREACT)
  REAL(kind=dp) :: JVS(LU_NONZERO)
  REAL(kind=dp) :: B(485)
  B(1) = RCT(1)
  B(2) = RCT(2)*F(2)
  B(4) = RCT(3)*V(92)
  B(5) = RCT(3)*V(88)
  B(6) = RCT(4)*V(93)*F(2)
  B(7) = RCT(4)*V(88)*F(2)
  B(9) = RCT(5)*V(101)
  B(10) = RCT(5)*V(88)
  B(11) = RCT(6)*V(101)
  B(12) = RCT(6)*V(88)
  B(13) = RCT(7)*V(93)
  B(14) = RCT(7)*V(92)
  B(15) = RCT(8)*V(101)
  B(16) = RCT(8)*V(92)
  B(17) = RCT(9)*V(99)
  B(18) = RCT(9)*V(93)
  B(19) = RCT(10)*2*V(93)*F(2)
  B(21) = RCT(11)*V(101)
  B(22) = RCT(11)*V(99)
  B(23) = RCT(12)
  B(24) = RCT(13)*F(1)
  B(26) = RCT(14)*V(101)
  B(27) = RCT(14)*V(99)
  B(28) = RCT(15)
  B(29) = RCT(16)
  B(30) = RCT(17)
  B(31) = RCT(18)
  B(32) = RCT(19)*F(1)
  B(34) = RCT(20)*F(2)
  B(36) = RCT(21)*V(103)
  B(37) = RCT(21)*V(93)
  B(38) = RCT(22)
  B(39) = RCT(23)
  B(40) = RCT(24)*V(103)
  B(41) = RCT(24)*V(47)
  B(42) = RCT(25)*V(103)
  B(43) = RCT(25)*V(101)
  B(44) = RCT(26)*V(103)
  B(45) = RCT(26)*V(99)
  B(46) = RCT(27)*V(103)
  B(47) = RCT(27)*V(70)
  B(48) = RCT(28)
  B(49) = RCT(29)*V(103)
  B(50) = RCT(29)*V(69)
  B(51) = RCT(30)*V(103)
  B(52) = RCT(30)*V(92)
  B(53) = RCT(31)*V(96)
  B(54) = RCT(31)*V(93)
  B(55) = RCT(32)*V(101)
  B(56) = RCT(32)*V(96)
  B(57) = RCT(33)
  B(58) = RCT(34)
  B(59) = RCT(35)*V(103)
  B(60) = RCT(35)*V(55)
  B(61) = RCT(36)*V(96)
  B(62) = RCT(36)*V(92)
  B(63) = RCT(37)*2*V(96)
  B(64) = RCT(38)*2*V(96)*F(1)
  B(66) = RCT(39)*V(99)
  B(67) = RCT(39)*V(96)
  B(68) = RCT(40)*2*V(99)
  B(69) = RCT(41)
  B(70) = RCT(42)*V(103)
  B(71) = RCT(42)*V(42)
  B(72) = RCT(43)*V(103)
  B(73) = RCT(43)*V(96)
  B(74) = RCT(44)*V(103)
  B(75) = RCT(44)*V(36)
  B(76) = RCT(45)*F(2)
  B(78) = RCT(46)*V(102)
  B(79) = RCT(46)*V(93)
  B(80) = RCT(47)*V(102)
  B(81) = RCT(47)*V(96)
  B(82) = RCT(48)*V(102)
  B(83) = RCT(48)*V(99)
  B(84) = RCT(49)*2*V(102)
  B(85) = RCT(50)*2*V(102)
  B(86) = RCT(51)*V(94)
  B(87) = RCT(51)*V(93)
  B(88) = RCT(52)*V(96)
  B(89) = RCT(52)*V(94)
  B(90) = RCT(53)*V(99)
  B(91) = RCT(53)*V(94)
  B(92) = RCT(54)*V(102)
  B(93) = RCT(54)*V(94)
  B(94) = RCT(55)*2*V(94)
  B(95) = RCT(56)*V(93)
  B(96) = RCT(56)*V(77)
  B(97) = RCT(57)*V(96)
  B(98) = RCT(57)*V(77)
  B(99) = RCT(58)*V(99)
  B(100) = RCT(58)*V(77)
  B(101) = RCT(59)*V(102)
  B(102) = RCT(59)*V(77)
  B(103) = RCT(60)*V(94)
  B(104) = RCT(60)*V(77)
  B(105) = RCT(61)*2*V(77)
  B(106) = RCT(62)*V(104)
  B(107) = RCT(62)*V(93)
  B(108) = RCT(63)*V(104)
  B(109) = RCT(63)*V(96)
  B(110) = RCT(64)*V(104)
  B(111) = RCT(64)*V(102)
  B(112) = RCT(65)*V(104)
  B(113) = RCT(65)*V(99)
  B(114) = RCT(66)*V(104)
  B(115) = RCT(66)*V(94)
  B(116) = RCT(67)*V(104)
  B(117) = RCT(67)*V(77)
  B(118) = RCT(68)*2*V(104)
  B(119) = RCT(69)*V(101)
  B(120) = RCT(69)*V(95)
  B(121) = RCT(70)
  B(122) = RCT(71)*V(95)
  B(123) = RCT(71)*V(93)
  B(124) = RCT(72)*V(96)
  B(125) = RCT(72)*V(95)
  B(126) = RCT(73)*V(99)
  B(127) = RCT(73)*V(95)
  B(128) = RCT(74)*V(102)
  B(129) = RCT(74)*V(95)
  B(130) = RCT(75)*V(95)
  B(131) = RCT(75)*V(94)
  B(132) = RCT(76)*V(95)
  B(133) = RCT(76)*V(77)
  B(134) = RCT(77)*V(104)
  B(135) = RCT(77)*V(95)
  B(136) = RCT(78)*2*V(95)
  B(137) = RCT(79)*V(101)
  B(138) = RCT(79)*V(97)
  B(139) = RCT(80)
  B(140) = RCT(81)*V(97)
  B(141) = RCT(81)*V(93)
  B(142) = RCT(82)*V(97)
  B(143) = RCT(82)*V(96)
  B(144) = RCT(83)*V(99)
  B(145) = RCT(83)*V(97)
  B(146) = RCT(84)*V(102)
  B(147) = RCT(84)*V(97)
  B(148) = RCT(85)*V(97)
  B(149) = RCT(85)*V(94)
  B(150) = RCT(86)*V(97)
  B(151) = RCT(86)*V(77)
  B(152) = RCT(87)*V(104)
  B(153) = RCT(87)*V(97)
  B(154) = RCT(88)*V(97)
  B(155) = RCT(88)*V(95)
  B(156) = RCT(89)*2*V(97)
  B(157) = RCT(90)*V(101)
  B(158) = RCT(90)*V(100)
  B(159) = RCT(91)
  B(160) = RCT(92)*V(100)
  B(161) = RCT(92)*V(93)
  B(162) = RCT(93)*V(100)
  B(163) = RCT(93)*V(96)
  B(164) = RCT(94)*V(100)
  B(165) = RCT(94)*V(99)
  B(166) = RCT(95)*V(102)
  B(167) = RCT(95)*V(100)
  B(168) = RCT(96)*V(100)
  B(169) = RCT(96)*V(94)
  B(170) = RCT(97)*V(100)
  B(171) = RCT(97)*V(77)
  B(172) = RCT(98)*V(104)
  B(173) = RCT(98)*V(100)
  B(174) = RCT(99)*V(100)
  B(175) = RCT(99)*V(95)
  B(176) = RCT(100)*V(100)
  B(177) = RCT(100)*V(97)
  B(178) = RCT(101)*2*V(100)
  B(179) = RCT(102)*V(101)
  B(180) = RCT(102)*V(98)
  B(181) = RCT(103)
  B(182) = RCT(104)*V(98)
  B(183) = RCT(104)*V(93)
  B(184) = RCT(105)*V(98)
  B(185) = RCT(105)*V(96)
  B(186) = RCT(106)*V(99)
  B(187) = RCT(106)*V(98)
  B(188) = RCT(107)*V(102)
  B(189) = RCT(107)*V(98)
  B(190) = RCT(108)*V(98)
  B(191) = RCT(108)*V(94)
  B(192) = RCT(109)*V(98)
  B(193) = RCT(109)*V(77)
  B(194) = RCT(110)*V(104)
  B(195) = RCT(110)*V(98)
  B(196) = RCT(111)*V(98)
  B(197) = RCT(111)*V(95)
  B(198) = RCT(112)*V(98)
  B(199) = RCT(112)*V(97)
  B(200) = RCT(113)*V(100)
  B(201) = RCT(113)*V(98)
  B(202) = RCT(114)*2*V(98)
  B(203) = RCT(115)*V(101)
  B(204) = RCT(115)*V(49)
  B(205) = RCT(116)
  B(206) = RCT(117)*V(101)
  B(207) = RCT(117)*V(74)
  B(208) = RCT(118)*V(96)
  B(209) = RCT(118)*V(74)
  B(210) = RCT(119)
  B(211) = RCT(120)*V(101)
  B(212) = RCT(120)*V(54)
  B(213) = RCT(121)*V(96)
  B(214) = RCT(121)*V(54)
  B(215) = RCT(122)
  B(216) = RCT(123)
  B(217) = RCT(124)
  B(218) = RCT(125)*V(103)
  B(219) = RCT(125)*V(86)
  B(220) = RCT(126)*V(96)
  B(221) = RCT(126)*V(86)
  B(222) = RCT(127)
  B(223) = RCT(128)*V(93)
  B(224) = RCT(128)*V(53)
  B(225) = RCT(129)*V(99)
  B(226) = RCT(129)*V(86)
  B(227) = RCT(130)*V(103)
  B(228) = RCT(130)*V(83)
  B(229) = RCT(131)
  B(230) = RCT(132)*V(99)
  B(231) = RCT(132)*V(83)
  B(232) = RCT(133)*V(103)
  B(233) = RCT(133)*V(89)
  B(234) = RCT(134)
  B(235) = RCT(135)*V(99)
  B(236) = RCT(135)*V(89)
  B(237) = RCT(136)*V(103)
  B(238) = RCT(136)*V(72)
  B(239) = RCT(137)
  B(240) = RCT(138)*V(103)
  B(241) = RCT(138)*V(90)
  B(242) = RCT(139)
  B(243) = RCT(140)*V(103)
  B(244) = RCT(140)*V(56)
  B(245) = RCT(141)*V(103)
  B(246) = RCT(141)*V(44)
  B(247) = RCT(142)*V(103)
  B(248) = RCT(142)*V(52)
  B(249) = RCT(143)
  B(250) = RCT(144)*V(103)
  B(251) = RCT(144)*V(65)
  B(252) = RCT(145)
  B(253) = RCT(146)
  B(254) = RCT(147)
  B(255) = RCT(148)*V(103)
  B(256) = RCT(148)*V(81)
  B(257) = RCT(149)*V(99)
  B(258) = RCT(149)*V(81)
  B(259) = RCT(150)
  B(260) = RCT(151)*V(103)
  B(261) = RCT(151)*V(68)
  B(262) = RCT(152)*V(99)
  B(263) = RCT(152)*V(68)
  B(264) = RCT(153)
  B(265) = RCT(154)*V(103)
  B(266) = RCT(154)*V(67)
  B(267) = RCT(155)*V(99)
  B(268) = RCT(155)*V(67)
  B(269) = RCT(156)*V(103)
  B(270) = RCT(156)*V(60)
  B(271) = RCT(157)*V(99)
  B(272) = RCT(157)*V(60)
  B(273) = RCT(158)*V(99)
  B(274) = RCT(158)*V(64)
  B(275) = RCT(159)*V(103)
  B(276) = RCT(159)*V(66)
  B(277) = RCT(160)
  B(278) = RCT(161)*V(99)
  B(279) = RCT(161)*V(66)
  B(280) = RCT(162)*V(103)
  B(281) = RCT(162)*V(78)
  B(282) = RCT(163)*V(92)
  B(283) = RCT(163)*V(78)
  B(284) = RCT(164)*V(99)
  B(285) = RCT(164)*V(78)
  B(286) = RCT(165)*V(88)
  B(287) = RCT(165)*V(78)
  B(288) = RCT(166)
  B(289) = RCT(167)*V(103)
  B(290) = RCT(167)*V(85)
  B(291) = RCT(168)*V(92)
  B(292) = RCT(168)*V(85)
  B(293) = RCT(169)*V(88)
  B(294) = RCT(169)*V(85)
  B(295) = RCT(170)
  B(296) = RCT(171)*V(103)
  B(297) = RCT(171)*V(82)
  B(298) = RCT(172)*V(92)
  B(299) = RCT(172)*V(82)
  B(300) = RCT(173)*V(99)
  B(301) = RCT(173)*V(82)
  B(302) = RCT(174)
  B(303) = RCT(175)*V(103)
  B(304) = RCT(175)*V(91)
  B(305) = RCT(176)
  B(306) = RCT(177)*V(103)
  B(307) = RCT(177)*V(87)
  B(308) = RCT(178)
  B(309) = RCT(179)*V(103)
  B(310) = RCT(179)*V(62)
  B(311) = RCT(180)*V(92)
  B(312) = RCT(180)*V(62)
  B(313) = RCT(181)*V(103)
  B(314) = RCT(181)*V(58)
  B(315) = RCT(182)
  B(316) = RCT(183)*V(103)
  B(317) = RCT(183)*V(59)
  B(318) = RCT(184)
  B(319) = RCT(185)*V(103)
  B(320) = RCT(185)*V(35)
  B(321) = RCT(186)*V(103)
  B(322) = RCT(186)*V(71)
  B(323) = RCT(187)*V(92)
  B(324) = RCT(187)*V(71)
  B(325) = RCT(188)*V(99)
  B(326) = RCT(188)*V(71)
  B(327) = RCT(189)*V(88)
  B(328) = RCT(189)*V(71)
  B(329) = RCT(190)*V(103)
  B(330) = RCT(190)*V(75)
  B(331) = RCT(191)*V(92)
  B(332) = RCT(191)*V(75)
  B(333) = RCT(192)*V(99)
  B(334) = RCT(192)*V(75)
  B(335) = RCT(193)*V(88)
  B(336) = RCT(193)*V(75)
  B(337) = RCT(194)*V(103)
  B(338) = RCT(194)*V(80)
  B(339) = RCT(195)*V(92)
  B(340) = RCT(195)*V(80)
  B(341) = RCT(196)*V(99)
  B(342) = RCT(196)*V(80)
  B(343) = RCT(197)*V(88)
  B(344) = RCT(197)*V(80)
  B(345) = RCT(198)*V(103)
  B(346) = RCT(198)*V(79)
  B(347) = RCT(199)*V(92)
  B(348) = RCT(199)*V(79)
  B(349) = RCT(200)*V(99)
  B(350) = RCT(200)*V(79)
  B(351) = RCT(201)*V(88)
  B(352) = RCT(201)*V(79)
  B(353) = RCT(202)*V(103)
  B(354) = RCT(202)*V(37)
  B(355) = RCT(203)*V(103)
  B(356) = RCT(203)*V(43)
  B(357) = RCT(204)*V(103)
  B(358) = RCT(204)*V(63)
  B(359) = RCT(205)*V(103)
  B(360) = RCT(205)*V(48)
  B(361) = RCT(206)*V(103)
  B(362) = RCT(206)*V(61)
  B(363) = RCT(207)*V(103)
  B(364) = RCT(207)*V(50)
  B(365) = RCT(208)*V(103)
  B(366) = RCT(208)*V(57)
  B(367) = RCT(209)*V(103)
  B(368) = RCT(209)*V(51)
  B(369) = RCT(210)*V(103)
  B(370) = RCT(210)*V(76)
  B(371) = RCT(211)*V(92)
  B(372) = RCT(211)*V(76)
  B(373) = RCT(212)*V(99)
  B(374) = RCT(212)*V(76)
  B(375) = RCT(213)*V(88)
  B(376) = RCT(213)*V(76)
  B(377) = RCT(214)*V(103)
  B(378) = RCT(214)*V(84)
  B(379) = RCT(215)*V(92)
  B(380) = RCT(215)*V(84)
  B(381) = RCT(216)*V(99)
  B(382) = RCT(216)*V(84)
  B(383) = RCT(217)*V(88)
  B(384) = RCT(217)*V(84)
  B(385) = RCT(218)*V(92)
  B(386) = RCT(218)*V(63)
  B(387) = RCT(219)*V(103)
  B(388) = RCT(219)*V(73)
  B(389) = RCT(220)*V(92)
  B(390) = RCT(220)*V(73)
  B(391) = RCT(221)*V(99)
  B(392) = RCT(221)*V(73)
  B(393) = RCT(222)*V(88)
  B(394) = RCT(222)*V(73)
  B(395) = RCT(223)
  B(396) = RCT(224)
  B(397) = RCT(225)
  B(398) = RCT(226)
  B(399) = RCT(227)
  B(400) = RCT(228)
  B(401) = RCT(229)
  B(402) = RCT(230)*V(103)
  B(403) = RCT(230)*V(57)
  B(404) = RCT(231)*V(103)
  B(405) = RCT(231)*V(57)
  B(406) = RCT(232)*V(103)
  B(407) = RCT(232)*V(57)
  B(408) = RCT(233)*V(103)
  B(409) = RCT(233)*V(57)
  B(410) = RCT(234)*V(103)
  B(411) = RCT(234)*V(51)
  B(412) = RCT(235)*V(103)
  B(413) = RCT(235)*V(51)
  B(414) = RCT(236)*V(103)
  B(415) = RCT(236)*V(51)
  B(416) = RCT(237)*V(103)
  B(417) = RCT(237)*V(51)
  B(418) = RCT(238)*V(103)
  B(419) = RCT(238)*V(75)
  B(420) = RCT(239)*V(103)
  B(421) = RCT(239)*V(75)
  B(422) = RCT(240)*V(103)
  B(423) = RCT(240)*V(75)
  B(424) = RCT(241)*V(103)
  B(425) = RCT(241)*V(80)
  B(426) = RCT(242)*V(103)
  B(427) = RCT(242)*V(80)
  B(428) = RCT(243)*V(103)
  B(429) = RCT(243)*V(80)
  B(430) = RCT(244)*V(103)
  B(431) = RCT(244)*V(79)
  B(432) = RCT(245)*V(103)
  B(433) = RCT(245)*V(79)
  B(434) = RCT(246)*V(92)
  B(435) = RCT(246)*V(75)
  B(436) = RCT(247)*V(92)
  B(437) = RCT(247)*V(80)
  B(438) = RCT(248)*V(92)
  B(439) = RCT(248)*V(80)
  B(440) = RCT(249)*V(92)
  B(441) = RCT(249)*V(80)
  B(442) = RCT(250)*V(92)
  B(443) = RCT(250)*V(79)
  B(444) = RCT(251)*V(92)
  B(445) = RCT(251)*V(79)
  B(446) = RCT(252)*V(92)
  B(447) = RCT(252)*V(79)
  B(448) = RCT(253)*V(99)
  B(449) = RCT(253)*V(75)
  B(450) = RCT(254)*V(99)
  B(451) = RCT(254)*V(75)
  B(452) = RCT(255)*V(99)
  B(453) = RCT(255)*V(75)
  B(454) = RCT(256)*V(99)
  B(455) = RCT(256)*V(80)
  B(456) = RCT(257)*V(99)
  B(457) = RCT(257)*V(80)
  B(458) = RCT(258)*V(99)
  B(459) = RCT(258)*V(79)
  B(460) = RCT(259)*V(99)
  B(461) = RCT(259)*V(79)
  B(462) = RCT(260)*V(103)
  B(463) = RCT(260)*V(3)
  B(464) = RCT(261)*V(103)
  B(465) = RCT(261)*V(4)
  B(466) = RCT(262)*V(103)
  B(467) = RCT(262)*V(32)
  B(468) = RCT(263)*V(103)
  B(469) = RCT(263)*V(26)
  B(470) = RCT(264)*V(103)
  B(471) = RCT(264)*V(5)
  B(472) = RCT(265)*V(103)
  B(473) = RCT(265)*V(6)
  B(474) = RCT(266)*V(103)
  B(475) = RCT(266)*V(34)
  B(476) = RCT(267)*V(103)
  B(477) = RCT(267)*V(28)
  B(478) = RCT(268)*V(103)
  B(479) = RCT(268)*V(31)
  B(480) = RCT(269)*V(103)
  B(481) = RCT(269)*V(27)
  B(482) = RCT(270)*V(103)
  B(483) = RCT(270)*V(33)
  B(484) = RCT(271)*V(103)
  B(485) = RCT(271)*V(29)
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
  JVS(17) = 0.204*B(331)
  JVS(18) = 0.185*B(371)
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
  JVS(33) = B(130)
  JVS(34) = 0.25*B(124)+B(128)+B(131)+B(134)
  JVS(35) = 0.25*B(125)
  JVS(36) = B(129)
  JVS(37) = B(135)
  JVS(38) = 0
  JVS(39) = 0.15*B(331)
  JVS(40) = 0.119*B(371)
  JVS(41) = 0.189*B(347)
  JVS(42) = 0.189*B(339)
  JVS(43) = 0.372*B(298)
  JVS(44) = 0.247*B(379)
  JVS(45) = 0.372*B(299)+0.15*B(332)+0.189*B(340)+0.189*B(348)+0.119*B(372)+0.247*B(380)
  JVS(46) = B(148)+B(168)+B(190)
  JVS(47) = 0.25*B(142)+0.25*B(162)+0.25*B(184)
  JVS(48) = 0.25*B(143)+B(146)+B(149)+B(152)
  JVS(49) = 0.25*B(185)+B(188)+B(191)+2*B(194)
  JVS(50) = 0.25*B(163)+B(166)+B(169)+B(172)
  JVS(51) = B(147)+B(167)+B(189)
  JVS(52) = B(153)+B(173)+2*B(195)
  JVS(53) = 0
  JVS(54) = 0.75*B(124)
  JVS(55) = 0.75*B(125)
  JVS(56) = 0
  JVS(57) = 0.75*B(142)+0.75*B(162)+0.75*B(184)
  JVS(58) = 0.75*B(143)
  JVS(59) = 0.75*B(185)
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
  JVS(73) = 6*B(212)
  JVS(74) = 0.048*B(388)
  JVS(75) = 0
  JVS(76) = B(95)+B(99)
  JVS(77) = B(78)+B(86)+B(96)+B(106)
  JVS(78) = B(87)+B(90)
  JVS(79) = B(82)+B(91)+B(100)+B(112)
  JVS(80) = B(79)+B(83)
  JVS(81) = B(107)+B(113)
  JVS(82) = 0
  JVS(83) = B(97)+B(101)+B(103)+B(105)+B(116)
  JVS(84) = B(88)+B(92)+B(94)+B(104)+B(114)
  JVS(85) = B(80)+B(89)+B(98)+B(108)
  JVS(86) = B(81)+B(84)+B(85)+B(93)+B(102)+B(110)
  JVS(87) = B(109)+B(111)+B(115)+B(117)+B(118)
  JVS(88) = 0
  JVS(89) = B(410)
  JVS(90) = B(402)
  JVS(91) = B(403)+B(411)
  JVS(92) = 0
  JVS(93) = B(412)
  JVS(94) = B(404)
  JVS(95) = B(405)+B(413)
  JVS(96) = 0
  JVS(97) = B(414)
  JVS(98) = B(406)
  JVS(99) = B(407)+B(415)
  JVS(100) = 0
  JVS(101) = B(416)
  JVS(102) = B(408)
  JVS(103) = B(409)+B(417)
  JVS(104) = 0
  JVS(105) = B(418)
  JVS(106) = B(419)
  JVS(107) = 0
  JVS(108) = B(420)+B(448)
  JVS(109) = B(430)+B(442)
  JVS(110) = B(424)+B(436)
  JVS(111) = B(437)+B(443)
  JVS(112) = B(449)
  JVS(113) = B(421)+B(425)+B(431)
  JVS(114) = 0
  JVS(115) = B(422)+B(434)+B(450)
  JVS(116) = B(432)+B(444)+B(458)
  JVS(117) = B(426)+B(438)+B(454)
  JVS(118) = B(435)+B(439)+B(445)
  JVS(119) = B(451)+B(455)+B(459)
  JVS(120) = B(423)+B(427)+B(433)
  JVS(121) = 0
  JVS(122) = B(452)
  JVS(123) = B(446)+B(460)
  JVS(124) = B(428)+B(440)+B(456)
  JVS(125) = B(441)+B(447)
  JVS(126) = B(453)+B(457)+B(461)
  JVS(127) = B(429)
  JVS(128) = 0
  JVS(129) = B(319)
  JVS(130) = B(74)
  JVS(131) = B(353)
  JVS(132) = B(355)
  JVS(133) = B(245)
  JVS(134) = B(40)
  JVS(135) = B(359)
  JVS(136) = B(363)
  JVS(137) = B(367)
  JVS(138) = B(247)
  JVS(139) = B(243)
  JVS(140) = B(365)
  JVS(141) = B(313)
  JVS(142) = B(316)
  JVS(143) = B(269)
  JVS(144) = B(361)
  JVS(145) = B(309)
  JVS(146) = B(357)
  JVS(147) = B(250)
  JVS(148) = B(275)
  JVS(149) = B(265)
  JVS(150) = B(260)
  JVS(151) = B(49)
  JVS(152) = B(46)
  JVS(153) = B(321)
  JVS(154) = B(237)
  JVS(155) = B(387)
  JVS(156) = B(329)
  JVS(157) = B(369)
  JVS(158) = B(280)
  JVS(159) = B(345)
  JVS(160) = B(337)
  JVS(161) = B(255)
  JVS(162) = B(296)
  JVS(163) = B(227)
  JVS(164) = B(377)
  JVS(165) = B(289)
  JVS(166) = B(218)
  JVS(167) = B(306)
  JVS(168) = B(232)
  JVS(169) = B(240)
  JVS(170) = B(303)
  JVS(171) = B(51)
  JVS(172) = B(36)
  JVS(173) = B(72)
  JVS(174) = B(44)
  JVS(175) = B(42)
  JVS(176) = B(37)+B(41)+B(43)+B(45)+B(47)+B(50)+B(52)+B(73)+B(75)+B(76)+B(219)+B(228)+B(233)+B(238)+B(241)+B(244)&
               &+B(246)+B(248)+B(251)+B(256)+B(261)+B(266)+B(270)+B(276)+B(281)+B(290)+B(297)+B(304)+B(307)+B(310)+B(314)&
               &+B(317)+B(320)+B(322)+B(330)+B(338)+B(346)+B(354)+B(356)+B(358)+B(360)+B(362)+B(364)+B(366)+B(368)+B(370)&
               &+B(378)+B(388)
  JVS(177) = B(462)
  JVS(178) = B(464)
  JVS(179) = B(470)
  JVS(180) = B(472)
  JVS(181) = 0
  JVS(182) = B(468)
  JVS(183) = B(480)
  JVS(184) = B(476)
  JVS(185) = B(484)
  JVS(186) = B(478)
  JVS(187) = B(466)
  JVS(188) = B(482)
  JVS(189) = B(474)
  JVS(190) = B(463)+B(465)+B(467)+B(469)+B(471)+B(473)+B(475)+B(477)+B(479)+B(481)+B(483)+B(485)
  JVS(191) = 0
  JVS(192) = B(480)
  JVS(193) = 0.5*B(478)
  JVS(194) = 0.5*B(479)+B(481)
  JVS(195) = -B(480)
  JVS(196) = -B(481)
  JVS(197) = 0
  JVS(198) = B(484)
  JVS(199) = 0.5*B(482)
  JVS(200) = 0.5*B(483)+B(485)
  JVS(201) = -B(484)
  JVS(202) = -B(485)
  JVS(203) = -B(32)-B(34)
  JVS(204) = B(31)
  JVS(205) = -B(478)
  JVS(206) = -B(479)
  JVS(207) = B(478)
  JVS(208) = 0
  JVS(209) = B(479)
  JVS(210) = -B(482)
  JVS(211) = -B(483)
  JVS(212) = B(482)
  JVS(213) = 0
  JVS(214) = B(483)
  JVS(215) = -B(319)
  JVS(216) = -B(320)
  JVS(217) = -B(74)-B(395)-B(397)
  JVS(218) = -B(75)
  JVS(219) = -B(353)
  JVS(220) = -B(354)
  JVS(221) = -B(121)
  JVS(222) = B(119)
  JVS(223) = B(120)
  JVS(224) = -B(139)
  JVS(225) = B(137)
  JVS(226) = B(138)
  JVS(227) = -B(159)
  JVS(228) = B(157)
  JVS(229) = B(158)
  JVS(230) = -B(181)
  JVS(231) = B(179)
  JVS(232) = B(180)
  JVS(233) = -B(69)-B(70)-B(400)
  JVS(234) = B(63)+B(64)
  JVS(235) = -B(71)
  JVS(236) = -B(355)
  JVS(237) = -B(356)
  JVS(238) = -B(245)
  JVS(239) = -B(246)
  JVS(240) = -B(264)
  JVS(241) = 0.087*B(367)
  JVS(242) = 0.031*B(347)
  JVS(243) = 0.031*B(339)
  JVS(244) = 0.031*B(340)+0.031*B(348)
  JVS(245) = 0.087*B(368)
  JVS(246) = -B(23)-B(24)
  JVS(247) = B(21)
  JVS(248) = B(22)
  JVS(249) = -B(38)-B(39)-B(40)
  JVS(250) = B(36)
  JVS(251) = B(37)-B(41)
  JVS(252) = -B(359)
  JVS(253) = -B(360)
  JVS(254) = 0.236*B(359)
  JVS(255) = -B(203)-B(205)
  JVS(256) = -B(204)
  JVS(257) = 0.236*B(360)
  JVS(258) = -B(363)
  JVS(259) = -B(364)
  JVS(260) = -B(367)
  JVS(261) = -B(368)
  JVS(262) = -B(247)-B(249)
  JVS(263) = B(80)
  JVS(264) = B(81)
  JVS(265) = -B(248)
  JVS(266) = -B(222)-B(223)
  JVS(267) = B(220)
  JVS(268) = -B(224)
  JVS(269) = B(221)
  JVS(270) = -B(211)-B(213)-B(215)
  JVS(271) = B(273)
  JVS(272) = -B(214)
  JVS(273) = B(274)
  JVS(274) = -B(212)
  JVS(275) = -B(57)-B(58)-B(59)
  JVS(276) = B(55)
  JVS(277) = B(56)
  JVS(278) = -B(60)
  JVS(279) = -B(243)
  JVS(280) = 0.25*B(92)
  JVS(281) = B(84)+0.25*B(93)+0.25*B(110)
  JVS(282) = -B(244)
  JVS(283) = 0.25*B(111)
  JVS(284) = -B(365)
  JVS(285) = -B(366)
  JVS(286) = 0.099*B(367)
  JVS(287) = 0.108*B(365)
  JVS(288) = -B(313)-B(315)
  JVS(289) = -B(314)+0.108*B(366)+0.099*B(368)
  JVS(290) = 0.093*B(367)
  JVS(291) = 0.051*B(365)
  JVS(292) = -B(316)-B(318)
  JVS(293) = -B(317)+0.051*B(366)+0.093*B(368)
  JVS(294) = 0.187*B(367)
  JVS(295) = 0.207*B(365)
  JVS(296) = -B(269)-B(271)
  JVS(297) = -B(272)
  JVS(298) = -B(270)+0.207*B(366)+0.187*B(368)
  JVS(299) = -B(361)
  JVS(300) = -B(362)
  JVS(301) = 0.561*B(367)
  JVS(302) = 0.491*B(365)
  JVS(303) = -B(309)-B(311)
  JVS(304) = -B(312)
  JVS(305) = -B(310)+0.491*B(366)+0.561*B(368)
  JVS(306) = -B(357)-B(385)
  JVS(307) = -B(386)
  JVS(308) = -B(358)
  JVS(309) = B(213)+B(215)
  JVS(310) = -B(273)
  JVS(311) = B(206)
  JVS(312) = B(214)
  JVS(313) = -B(274)
  JVS(314) = B(207)
  JVS(315) = -B(250)-B(252)
  JVS(316) = B(88)
  JVS(317) = B(89)+B(108)
  JVS(318) = -B(251)
  JVS(319) = B(109)
  JVS(320) = 0.05*B(367)
  JVS(321) = 0.059*B(365)
  JVS(322) = -B(275)-B(277)-B(278)
  JVS(323) = 0.061*B(377)+0.042*B(379)+0.015*B(381)
  JVS(324) = 0.042*B(380)
  JVS(325) = -B(279)+0.015*B(382)
  JVS(326) = -B(276)+0.059*B(366)+0.05*B(368)+0.061*B(378)
  JVS(327) = 0.017*B(365)
  JVS(328) = -B(265)-B(267)
  JVS(329) = B(208)+B(210)
  JVS(330) = B(209)
  JVS(331) = -B(268)
  JVS(332) = -B(266)+0.017*B(366)
  JVS(333) = 0.287*B(367)
  JVS(334) = 0.119*B(365)
  JVS(335) = 0.5*B(315)
  JVS(336) = 0.5*B(318)
  JVS(337) = 0.23*B(269)
  JVS(338) = -B(259)-B(260)-B(262)
  JVS(339) = 0.084*B(280)+0.9*B(282)
  JVS(340) = 0.174*B(296)+0.742*B(298)+0.008*B(300)
  JVS(341) = 0.3*B(289)+0.95*B(291)
  JVS(342) = 0.9*B(283)+0.95*B(292)+0.742*B(299)
  JVS(343) = -B(263)+0.008*B(301)
  JVS(344) = -B(261)+0.23*B(270)+0.084*B(281)+0.3*B(290)+0.174*B(297)+0.119*B(366)+0.287*B(368)
  JVS(345) = B(315)
  JVS(346) = B(318)
  JVS(347) = 0.002*B(361)
  JVS(348) = B(309)+1.5*B(311)
  JVS(349) = 0.393*B(357)+1.5*B(385)
  JVS(350) = B(259)+B(260)+B(262)
  JVS(351) = -B(49)
  JVS(352) = 0.5*B(323)+0.491*B(327)
  JVS(353) = 0.51*B(389)
  JVS(354) = 0.275*B(331)
  JVS(355) = 0.345*B(371)
  JVS(356) = 0.416*B(280)+0.45*B(282)+0.5*B(284)+0.67*B(288)
  JVS(357) = 0.157*B(347)
  JVS(358) = 0.157*B(339)
  JVS(359) = 2*B(253)+B(254)+1.26*B(255)+1.26*B(257)
  JVS(360) = 0.336*B(296)+0.498*B(298)+0.572*B(300)+1.233*B(302)
  JVS(361) = B(229)
  JVS(362) = 0.265*B(379)+0.012*B(383)
  JVS(363) = 0.475*B(291)+0.7*B(295)
  JVS(364) = B(216)+B(217)+B(218)+B(225)
  JVS(365) = 0.491*B(328)+0.012*B(384)
  JVS(366) = 0.034*B(232)+B(234)
  JVS(367) = 0.45*B(283)+0.475*B(292)+0.498*B(299)+1.5*B(312)+0.5*B(324)+0.275*B(332)+0.157*B(340)+0.157*B(348)+0.345&
               &*B(372)+0.265*B(380)+1.5*B(386)+0.51*B(390)
  JVS(368) = B(226)+1.26*B(258)+B(263)+0.5*B(285)+0.572*B(301)
  JVS(369) = -B(50)+B(219)+0.034*B(233)+1.26*B(256)+B(261)+0.416*B(281)+0.336*B(297)+B(310)+0.393*B(358)+0.002*B(362)
  JVS(370) = 2*B(24)
  JVS(371) = B(271)
  JVS(372) = B(273)
  JVS(373) = B(278)
  JVS(374) = B(267)
  JVS(375) = B(262)
  JVS(376) = -B(46)-B(48)-B(399)
  JVS(377) = 0
  JVS(378) = 0.5*B(284)
  JVS(379) = B(257)
  JVS(380) = 0.15*B(300)
  JVS(381) = B(230)
  JVS(382) = 0
  JVS(383) = 0
  JVS(384) = B(225)
  JVS(385) = B(235)
  JVS(386) = 0
  JVS(387) = 0.2*B(66)
  JVS(388) = 0.2*B(67)+B(226)+B(231)+B(236)+B(258)+B(263)+B(268)+B(272)+B(274)+B(279)+0.5*B(285)+0.15*B(301)
  JVS(389) = B(42)
  JVS(390) = B(43)-B(47)
  JVS(391) = -B(321)-B(323)-B(325)-B(327)
  JVS(392) = -B(328)
  JVS(393) = -B(324)
  JVS(394) = -B(326)
  JVS(395) = -B(322)
  JVS(396) = 0.704*B(355)
  JVS(397) = 0.024*B(359)
  JVS(398) = B(205)
  JVS(399) = 0.072*B(363)
  JVS(400) = 0.452*B(361)
  JVS(401) = -B(237)-B(239)
  JVS(402) = 0.005*B(369)+0.001*B(371)+0.024*B(373)
  JVS(403) = 0.13*B(347)
  JVS(404) = 0.13*B(339)
  JVS(405) = 0.127*B(377)+0.045*B(379)+0.102*B(381)
  JVS(406) = 0.006*B(306)+0.02*B(308)
  JVS(407) = 0.13*B(340)+0.13*B(348)+0.001*B(372)+0.045*B(380)
  JVS(408) = 0.024*B(374)+0.102*B(382)
  JVS(409) = 0
  JVS(410) = -B(238)+0.006*B(307)+0.704*B(356)+0.024*B(360)+0.452*B(362)+0.072*B(364)+0.005*B(370)+0.127*B(378)
  JVS(411) = -B(387)-B(389)-B(391)-B(393)
  JVS(412) = -B(394)
  JVS(413) = -B(390)
  JVS(414) = -B(392)
  JVS(415) = -B(388)
  JVS(416) = 0.24*B(269)+B(271)
  JVS(417) = 0.24*B(265)+B(267)
  JVS(418) = -B(206)-B(208)-B(210)
  JVS(419) = B(160)
  JVS(420) = B(174)
  JVS(421) = -B(209)
  JVS(422) = B(176)
  JVS(423) = B(200)
  JVS(424) = B(164)+B(268)+B(272)
  JVS(425) = B(161)+B(165)+B(175)+B(177)+2*B(178)+B(201)
  JVS(426) = -B(207)
  JVS(427) = 0.24*B(266)+0.24*B(270)
  JVS(428) = -B(329)-B(331)-B(333)-B(335)
  JVS(429) = -B(336)
  JVS(430) = -B(332)
  JVS(431) = -B(334)
  JVS(432) = -B(330)
  JVS(433) = -B(369)-B(371)-B(373)-B(375)
  JVS(434) = -B(376)
  JVS(435) = -B(372)
  JVS(436) = -B(374)
  JVS(437) = -B(370)
  JVS(438) = 0.559*B(359)
  JVS(439) = 0.948*B(363)
  JVS(440) = B(313)+B(315)
  JVS(441) = B(316)+B(318)
  JVS(442) = 0.936*B(361)
  JVS(443) = B(237)
  JVS(444) = 0.079*B(329)+0.126*B(331)+0.187*B(333)+0.24*B(335)
  JVS(445) = 0.205*B(369)+0.488*B(373)
  JVS(446) = -B(95)-B(97)-B(99)-B(101)-B(103)-B(116)-B(132)-B(150)-B(170)-B(192)
  JVS(447) = 0.5*B(345)+0.729*B(347)+0.75*B(349)
  JVS(448) = 0.5*B(337)+0.729*B(339)+0.75*B(341)
  JVS(449) = 0.001*B(377)+0.137*B(379)+0.711*B(381)
  JVS(450) = 0.675*B(289)
  JVS(451) = 0.596*B(306)+0.152*B(308)
  JVS(452) = 0.24*B(336)
  JVS(453) = 0.616*B(240)
  JVS(454) = 0.515*B(305)
  JVS(455) = 0.126*B(332)+0.729*B(340)+0.729*B(348)+0.137*B(380)
  JVS(456) = -B(96)+B(160)
  JVS(457) = -B(104)
  JVS(458) = -B(133)+B(174)
  JVS(459) = -B(98)
  JVS(460) = -B(151)+B(176)
  JVS(461) = -B(193)+B(200)
  JVS(462) = -B(100)+B(164)+0.187*B(334)+0.75*B(342)+0.75*B(350)+0.488*B(374)+0.711*B(382)
  JVS(463) = B(161)+B(165)-B(171)+B(175)+B(177)+2*B(178)+B(201)
  JVS(464) = 0
  JVS(465) = -B(102)
  JVS(466) = B(238)+0.616*B(241)+0.675*B(290)+0.596*B(307)+B(314)+B(317)+0.079*B(330)+0.5*B(338)+0.5*B(346)+0.559*B(360)&
               &+0.936*B(362)+0.948*B(364)+0.205*B(370)+0.001*B(378)
  JVS(467) = -B(117)
  JVS(468) = 0.23*B(329)+0.39*B(331)
  JVS(469) = -B(280)-B(282)-B(284)-B(286)-B(288)
  JVS(470) = 0.025*B(377)+0.026*B(379)+0.012*B(383)
  JVS(471) = -B(287)+0.012*B(384)
  JVS(472) = -B(283)+0.39*B(332)+0.026*B(380)
  JVS(473) = -B(285)
  JVS(474) = -B(281)+0.23*B(330)+0.025*B(378)
  JVS(475) = -B(345)-B(347)-B(349)-B(351)
  JVS(476) = -B(352)
  JVS(477) = -B(348)
  JVS(478) = -B(350)
  JVS(479) = -B(346)
  JVS(480) = -B(337)-B(339)-B(341)-B(343)
  JVS(481) = -B(344)
  JVS(482) = -B(340)
  JVS(483) = -B(342)
  JVS(484) = -B(338)
  JVS(485) = 0.097*B(367)
  JVS(486) = 0.118*B(365)
  JVS(487) = 0.5*B(315)
  JVS(488) = 0.5*B(318)
  JVS(489) = B(311)
  JVS(490) = 0.607*B(357)
  JVS(491) = 0.23*B(265)
  JVS(492) = 0.009*B(327)
  JVS(493) = 0
  JVS(494) = 0.001*B(347)
  JVS(495) = 0.001*B(339)
  JVS(496) = -B(253)-B(254)-B(255)-B(257)
  JVS(497) = 0.15*B(296)+0.023*B(298)
  JVS(498) = 0.009*B(328)
  JVS(499) = 0.023*B(299)+B(312)+0.001*B(340)+0.001*B(348)
  JVS(500) = 0
  JVS(501) = 0
  JVS(502) = 0
  JVS(503) = 0
  JVS(504) = 0
  JVS(505) = -B(258)
  JVS(506) = 0
  JVS(507) = 0
  JVS(508) = -B(256)+0.23*B(266)+0.15*B(297)+0.607*B(358)+0.118*B(366)+0.097*B(368)
  JVS(509) = 0.357*B(329)+0.936*B(333)
  JVS(510) = -B(296)-B(298)-B(300)-B(302)
  JVS(511) = 0.025*B(377)
  JVS(512) = 0
  JVS(513) = -B(299)
  JVS(514) = -B(301)+0.936*B(334)
  JVS(515) = -B(297)+0.357*B(330)+0.025*B(378)
  JVS(516) = B(353)
  JVS(517) = 0.96*B(245)
  JVS(518) = 0.445*B(359)
  JVS(519) = 0.099*B(363)
  JVS(520) = 0.455*B(361)
  JVS(521) = 0.195*B(321)+0.25*B(327)
  JVS(522) = 0.984*B(387)+0.5*B(389)
  JVS(523) = 0.294*B(369)+0.154*B(371)+0.009*B(373)
  JVS(524) = 0.129*B(296)+0.047*B(298)+0.467*B(302)
  JVS(525) = -B(227)-B(229)-B(230)
  JVS(526) = 0.732*B(377)+0.456*B(379)+0.507*B(381)
  JVS(527) = 0.439*B(306)+0.431*B(308)
  JVS(528) = 0.25*B(328)
  JVS(529) = 0.034*B(232)+B(234)
  JVS(530) = 0.482*B(240)+B(242)
  JVS(531) = 0.084*B(303)+0.246*B(305)
  JVS(532) = 0.047*B(299)+0.154*B(372)+0.456*B(380)+0.5*B(390)
  JVS(533) = B(140)
  JVS(534) = B(154)
  JVS(535) = B(141)+B(144)+B(155)+2*B(156)+B(176)+B(198)
  JVS(536) = B(199)
  JVS(537) = B(145)-B(231)+0.009*B(374)+0.507*B(382)
  JVS(538) = B(177)
  JVS(539) = -B(228)+0.034*B(233)+0.482*B(241)+0.96*B(246)+0.129*B(297)+0.084*B(304)+0.439*B(307)+0.195*B(322)+B(354)&
               &+0.445*B(360)+0.455*B(362)+0.099*B(364)+0.294*B(370)+0.732*B(378)+0.984*B(388)
  JVS(540) = -B(377)-B(379)-B(381)-B(383)
  JVS(541) = -B(384)
  JVS(542) = -B(380)
  JVS(543) = -B(382)
  JVS(544) = -B(378)
  JVS(545) = 0.32*B(329)+0.16*B(331)
  JVS(546) = 0.019*B(379)+0.048*B(381)
  JVS(547) = -B(289)-B(291)-B(293)-B(295)
  JVS(548) = -B(294)
  JVS(549) = -B(292)+0.16*B(332)+0.019*B(380)
  JVS(550) = 0.048*B(382)
  JVS(551) = -B(290)+0.32*B(330)
  JVS(552) = 0.081*B(245)
  JVS(553) = 0.026*B(359)
  JVS(554) = 0.026*B(363)
  JVS(555) = 0.35*B(247)+B(249)
  JVS(556) = B(222)
  JVS(557) = B(243)
  JVS(558) = 0.024*B(361)
  JVS(559) = 0.096*B(357)
  JVS(560) = 1.61*B(321)+B(323)+0.191*B(327)
  JVS(561) = B(237)
  JVS(562) = 0.984*B(387)+0.5*B(389)
  JVS(563) = 0.624*B(329)+0.592*B(331)+0.24*B(335)
  JVS(564) = 0.732*B(369)+0.5*B(371)
  JVS(565) = 0.084*B(280)+0.2*B(282)+0.67*B(288)
  JVS(566) = 0.276*B(345)+0.235*B(347)
  JVS(567) = 0.276*B(337)+0.235*B(339)
  JVS(568) = B(254)
  JVS(569) = 0.055*B(296)+0.125*B(298)+0.227*B(300)+0.3*B(302)
  JVS(570) = 0.244*B(377)+0.269*B(379)+0.079*B(381)
  JVS(571) = 0.3*B(289)+0.1*B(291)
  JVS(572) = -B(216)-B(217)-B(218)-B(220)-B(225)
  JVS(573) = 0.01*B(306)+0.134*B(308)
  JVS(574) = 0.191*B(328)+0.24*B(336)
  JVS(575) = 0.115*B(240)
  JVS(576) = 0.213*B(303)+0.506*B(305)
  JVS(577) = 0.2*B(283)+0.1*B(292)+0.125*B(299)+B(324)+0.592*B(332)+0.235*B(340)+0.235*B(348)+0.5*B(372)+0.269*B(380)&
               &+0.5*B(390)
  JVS(578) = B(78)+B(182)
  JVS(579) = 0.75*B(92)
  JVS(580) = B(128)+B(196)
  JVS(581) = -B(221)
  JVS(582) = B(146)+B(198)
  JVS(583) = B(183)+B(186)+B(188)+B(197)+B(199)+B(200)+2*B(202)
  JVS(584) = B(82)+B(187)-B(226)+0.227*B(301)+0.079*B(382)
  JVS(585) = B(166)+B(201)
  JVS(586) = 0
  JVS(587) = B(79)+B(83)+B(84)+2*B(85)+0.75*B(93)+0.75*B(110)+B(129)+B(147)+B(167)+B(189)
  JVS(588) = -B(219)+B(238)+0.115*B(241)+B(244)+0.081*B(246)+0.35*B(248)+0.084*B(281)+0.3*B(290)+0.055*B(297)+0.213&
               &*B(304)+0.01*B(307)+1.61*B(322)+0.624*B(330)+0.276*B(338)+0.276*B(346)+0.096*B(358)+0.026*B(360)+0.024&
               &*B(362)+0.026*B(364)+0.732*B(370)+0.244*B(378)+0.984*B(388)
  JVS(589) = 0.75*B(111)
  JVS(590) = B(203)
  JVS(591) = 0.511*B(373)
  JVS(592) = 0.276*B(349)
  JVS(593) = 0.276*B(341)
  JVS(594) = 0.572*B(300)
  JVS(595) = 0.321*B(381)
  JVS(596) = -0.69*B(306)-B(308)
  JVS(597) = 0
  JVS(598) = 0
  JVS(599) = B(106)
  JVS(600) = 0.572*B(301)+0.276*B(342)+0.276*B(350)+0.511*B(374)+0.321*B(382)
  JVS(601) = B(204)
  JVS(602) = -0.69*B(307)
  JVS(603) = B(107)
  JVS(604) = B(34)
  JVS(605) = -B(327)
  JVS(606) = -B(393)
  JVS(607) = -B(335)
  JVS(608) = -B(375)
  JVS(609) = -B(286)
  JVS(610) = -B(351)
  JVS(611) = -B(343)
  JVS(612) = -B(383)
  JVS(613) = -B(293)
  JVS(614) = -B(2)-B(4)-B(6)-B(9)-B(11)-B(287)-B(294)-B(328)-B(336)-B(344)-B(352)-B(376)-B(384)-B(394)
  JVS(615) = -B(5)+B(30)
  JVS(616) = -B(7)
  JVS(617) = B(29)
  JVS(618) = B(1)-B(10)-B(12)
  JVS(619) = 0
  JVS(620) = 0.261*B(355)
  JVS(621) = 0.122*B(359)
  JVS(622) = 0.204*B(363)
  JVS(623) = B(313)
  JVS(624) = B(316)
  JVS(625) = 0.244*B(361)
  JVS(626) = B(309)
  JVS(627) = B(250)+B(252)
  JVS(628) = B(325)
  JVS(629) = 0.45*B(393)
  JVS(630) = 0.497*B(369)+0.363*B(371)+0.037*B(373)+0.45*B(375)
  JVS(631) = B(286)
  JVS(632) = 0.474*B(345)+0.205*B(347)+0.474*B(349)+0.147*B(351)
  JVS(633) = 0.474*B(337)+0.205*B(339)+0.474*B(341)+0.147*B(343)
  JVS(634) = 0.013*B(296)+0.218*B(300)
  JVS(635) = 0.511*B(377)+0.305*B(379)+0.151*B(381)+0.069*B(383)
  JVS(636) = 0.675*B(289)+0.45*B(293)
  JVS(637) = 0.213*B(306)+0.147*B(308)
  JVS(638) = B(287)+0.45*B(294)+0.147*B(344)+0.147*B(352)+0.45*B(376)+0.069*B(384)+0.45*B(394)
  JVS(639) = -B(232)-B(234)-B(235)
  JVS(640) = 0.37*B(240)
  JVS(641) = 0.558*B(303)+0.71*B(305)
  JVS(642) = 0.205*B(340)+0.205*B(348)+0.363*B(372)+0.305*B(380)
  JVS(643) = 0
  JVS(644) = 0
  JVS(645) = 0
  JVS(646) = -B(236)+0.218*B(301)+B(326)+0.474*B(342)+0.474*B(350)+0.037*B(374)+0.151*B(382)
  JVS(647) = 0
  JVS(648) = -B(233)+0.37*B(241)+B(251)+0.675*B(290)+0.013*B(297)+0.558*B(304)+0.213*B(307)+B(310)+B(314)+B(317)+0.474&
               &*B(338)+0.474*B(346)+0.261*B(356)+0.122*B(360)+0.244*B(362)+0.204*B(364)+0.497*B(370)+0.511*B(378)
  JVS(649) = 0
  JVS(650) = 0.332*B(359)
  JVS(651) = 0.089*B(363)
  JVS(652) = 0.11*B(361)
  JVS(653) = 0.55*B(393)
  JVS(654) = 0.437*B(375)
  JVS(655) = 0.416*B(280)
  JVS(656) = 0.15*B(296)+0.21*B(298)+0.233*B(302)
  JVS(657) = 0.072*B(377)+0.026*B(379)+0.001*B(381)+0.659*B(383)
  JVS(658) = 0.55*B(293)
  JVS(659) = 0.177*B(306)+0.243*B(308)
  JVS(660) = 0.55*B(294)+0.437*B(376)+0.659*B(384)+0.55*B(394)
  JVS(661) = -B(240)-B(242)
  JVS(662) = 0.115*B(303)
  JVS(663) = 0.21*B(299)+0.026*B(380)
  JVS(664) = 0
  JVS(665) = 0.5*B(114)
  JVS(666) = B(112)+0.001*B(382)
  JVS(667) = 0
  JVS(668) = 0.5*B(110)
  JVS(669) = -B(241)+0.416*B(281)+0.15*B(297)+0.115*B(304)+0.177*B(307)+0.332*B(360)+0.11*B(362)+0.089*B(364)+0.072&
               &*B(378)
  JVS(670) = 0.5*B(111)+B(113)+0.5*B(115)+B(118)
  JVS(671) = 0.417*B(363)
  JVS(672) = 0.055*B(365)
  JVS(673) = 0.125*B(361)
  JVS(674) = 0.1*B(331)+0.75*B(335)
  JVS(675) = 0.119*B(369)+0.215*B(371)+0.113*B(375)
  JVS(676) = 0.276*B(345)+0.276*B(347)+0.853*B(351)
  JVS(677) = 0.276*B(337)+0.276*B(339)+0.853*B(343)
  JVS(678) = 0.332*B(296)
  JVS(679) = 0.043*B(379)+0.259*B(383)
  JVS(680) = 0.7*B(295)
  JVS(681) = 0.048*B(306)+0.435*B(308)
  JVS(682) = 0.75*B(336)+0.853*B(344)+0.853*B(352)+0.113*B(376)+0.259*B(384)
  JVS(683) = -0.671*B(303)-B(305)
  JVS(684) = 0.1*B(332)+0.276*B(340)+0.276*B(348)+0.215*B(372)+0.043*B(380)
  JVS(685) = 0
  JVS(686) = 0.5*B(114)
  JVS(687) = B(134)
  JVS(688) = B(152)
  JVS(689) = 0
  JVS(690) = B(172)
  JVS(691) = 0
  JVS(692) = 0.5*B(110)
  JVS(693) = 0.332*B(297)-0.671*B(304)+0.048*B(307)+0.276*B(338)+0.276*B(346)+0.125*B(362)+0.417*B(364)+0.055*B(366)&
               &+0.119*B(370)
  JVS(694) = 0.5*B(111)+0.5*B(115)+B(118)+B(135)+B(153)+B(173)
  JVS(695) = -B(311)
  JVS(696) = -B(385)
  JVS(697) = -B(323)
  JVS(698) = -B(389)
  JVS(699) = -B(331)
  JVS(700) = -B(371)
  JVS(701) = -B(282)
  JVS(702) = -B(347)
  JVS(703) = -B(339)
  JVS(704) = -B(298)
  JVS(705) = -B(379)
  JVS(706) = -B(291)
  JVS(707) = B(2)-B(4)
  JVS(708) = -B(5)-B(13)-B(15)-B(30)-B(31)-B(51)-B(61)-B(283)-B(292)-B(299)-B(312)-B(324)-B(332)-B(340)-B(348)-B(372)&
               &-B(380)-B(386)-B(390)
  JVS(709) = -B(14)
  JVS(710) = 0.25*B(124)
  JVS(711) = -B(62)+0.25*B(125)+0.25*B(142)+0.25*B(162)+0.25*B(184)
  JVS(712) = 0.25*B(143)
  JVS(713) = 0.25*B(185)
  JVS(714) = 0
  JVS(715) = 0.25*B(163)
  JVS(716) = -B(16)
  JVS(717) = -B(52)
  JVS(718) = B(38)
  JVS(719) = -B(223)
  JVS(720) = -B(95)
  JVS(721) = 0
  JVS(722) = 0
  JVS(723) = 0
  JVS(724) = 0
  JVS(725) = 0
  JVS(726) = 0
  JVS(727) = -B(6)+B(9)
  JVS(728) = 0
  JVS(729) = 0
  JVS(730) = -B(13)
  JVS(731) = -B(7)-B(14)-B(17)-2*B(19)-B(36)-B(53)-B(78)-B(86)-B(96)-B(106)-B(122)-B(140)-B(160)-B(182)-B(224)
  JVS(732) = -B(87)
  JVS(733) = -B(123)
  JVS(734) = -B(54)
  JVS(735) = -B(141)
  JVS(736) = -B(183)
  JVS(737) = -B(18)+B(26)+B(28)
  JVS(738) = -B(161)
  JVS(739) = B(1)+B(10)+B(27)
  JVS(740) = -B(79)
  JVS(741) = -B(37)
  JVS(742) = -B(107)
  JVS(743) = B(353)
  JVS(744) = 0.965*B(355)
  JVS(745) = 0.05*B(245)
  JVS(746) = 0.695*B(359)
  JVS(747) = 0.653*B(363)
  JVS(748) = 0.804*B(367)
  JVS(749) = 0.765*B(365)
  JVS(750) = B(315)
  JVS(751) = B(318)
  JVS(752) = 0.76*B(269)
  JVS(753) = 0.835*B(361)
  JVS(754) = B(309)
  JVS(755) = 0.1*B(357)
  JVS(756) = 0.34*B(250)
  JVS(757) = 0.76*B(265)
  JVS(758) = B(321)+B(325)+0.2*B(327)
  JVS(759) = 0.984*B(387)+0.949*B(391)
  JVS(760) = 0
  JVS(761) = 0.907*B(329)+0.066*B(331)+0.749*B(333)
  JVS(762) = 0.91*B(369)+0.022*B(371)+0.824*B(373)
  JVS(763) = 0.5*B(280)+0.1*B(282)+0.5*B(284)+0.33*B(288)
  JVS(764) = 0.75*B(345)+0.031*B(347)+0.276*B(349)
  JVS(765) = 0.75*B(337)+0.031*B(339)+0.276*B(341)
  JVS(766) = 0.67*B(296)+0.048*B(298)+0.799*B(300)
  JVS(767) = 0.918*B(377)+0.033*B(379)+0.442*B(381)+0.012*B(383)
  JVS(768) = 0.3*B(289)+0.05*B(291)
  JVS(769) = 0.376*B(306)+0.564*B(308)
  JVS(770) = 0.2*B(328)+0.012*B(384)
  JVS(771) = 0.034*B(232)+B(234)
  JVS(772) = 0.37*B(240)+B(242)
  JVS(773) = 0.473*B(303)+0.96*B(305)
  JVS(774) = 0.1*B(283)+0.05*B(292)+0.048*B(299)+0.066*B(332)+0.031*B(340)+0.031*B(348)+0.022*B(372)+0.033*B(380)
  JVS(775) = -B(86)+B(140)
  JVS(776) = -B(87)-B(88)-B(90)-B(92)-2*B(94)-B(114)-B(130)-B(148)-B(168)-B(190)
  JVS(777) = -B(131)+B(154)
  JVS(778) = -B(89)
  JVS(779) = B(141)+B(144)-B(149)+B(155)+2*B(156)+B(176)+B(198)
  JVS(780) = -B(191)+B(199)
  JVS(781) = -B(91)+B(145)+0.5*B(285)+0.799*B(301)+B(326)+0.749*B(334)+0.276*B(342)+0.276*B(350)+0.824*B(374)+0.442&
               &*B(382)+0.949*B(392)
  JVS(782) = -B(169)+B(177)
  JVS(783) = 0
  JVS(784) = -B(93)
  JVS(785) = 0.034*B(233)+0.37*B(241)+0.05*B(246)+0.34*B(251)+0.76*B(266)+0.76*B(270)+0.5*B(281)+0.3*B(290)+0.67*B(297)&
               &+0.473*B(304)+0.376*B(307)+B(310)+B(322)+0.907*B(330)+0.75*B(338)+0.75*B(346)+B(354)+0.965*B(356)+0.1*B(358)&
               &+0.695*B(360)+0.835*B(362)+0.653*B(364)+0.765*B(366)+0.804*B(368)+0.91*B(370)+0.918*B(378)+0.984*B(388)
  JVS(786) = -B(115)
  JVS(787) = B(121)
  JVS(788) = 2*B(264)
  JVS(789) = 0
  JVS(790) = B(313)+0.5*B(315)
  JVS(791) = B(316)+0.5*B(318)
  JVS(792) = 0.011*B(361)
  JVS(793) = B(259)+B(260)+B(262)
  JVS(794) = B(237)+B(239)
  JVS(795) = 0
  JVS(796) = 0.67*B(288)
  JVS(797) = 0.123*B(347)
  JVS(798) = 0.123*B(339)
  JVS(799) = 0.467*B(302)
  JVS(800) = B(227)+B(230)
  JVS(801) = 0.137*B(379)
  JVS(802) = 0.675*B(289)
  JVS(803) = 0
  JVS(804) = 0
  JVS(805) = 0
  JVS(806) = 0.492*B(240)+B(242)
  JVS(807) = 0.029*B(303)+0.667*B(305)
  JVS(808) = 0.123*B(340)+0.123*B(348)+0.137*B(380)
  JVS(809) = -B(122)+B(182)
  JVS(810) = -B(130)
  JVS(811) = -B(119)-B(123)-B(124)-B(126)-B(128)-B(131)-B(134)-2*B(136)-B(154)-B(174)
  JVS(812) = -B(125)
  JVS(813) = -B(155)+B(198)
  JVS(814) = B(183)+B(186)+B(199)+B(200)+2*B(202)
  JVS(815) = -B(127)+B(187)+B(231)+B(263)
  JVS(816) = -B(175)+B(201)
  JVS(817) = -B(120)
  JVS(818) = -B(129)
  JVS(819) = B(228)+B(238)+0.492*B(241)+B(261)+0.675*B(290)+0.029*B(304)+B(314)+B(317)+0.011*B(362)
  JVS(820) = -B(135)
  JVS(821) = B(74)
  JVS(822) = B(70)
  JVS(823) = 0.95*B(245)
  JVS(824) = B(39)
  JVS(825) = 0.187*B(367)
  JVS(826) = B(249)
  JVS(827) = B(222)+B(223)
  JVS(828) = -B(213)
  JVS(829) = B(57)+0.61*B(58)
  JVS(830) = B(243)
  JVS(831) = 0.224*B(365)
  JVS(832) = 0.5*B(315)
  JVS(833) = 0.5*B(318)
  JVS(834) = 1.5*B(311)
  JVS(835) = 0.297*B(357)+1.5*B(385)
  JVS(836) = 0
  JVS(837) = B(252)
  JVS(838) = B(259)
  JVS(839) = B(49)
  JVS(840) = 0.12*B(323)+0.5*B(327)
  JVS(841) = 0.06*B(389)
  JVS(842) = -B(208)
  JVS(843) = 0
  JVS(844) = 0.056*B(371)
  JVS(845) = 0.008*B(282)+0.34*B(288)
  JVS(846) = 0.033*B(347)
  JVS(847) = 0.033*B(339)
  JVS(848) = 2*B(253)+0.63*B(255)+0.63*B(257)
  JVS(849) = 0.4*B(298)+1.233*B(302)
  JVS(850) = B(229)
  JVS(851) = 0.003*B(379)+0.013*B(383)
  JVS(852) = 0.064*B(291)
  JVS(853) = 2*B(216)+B(218)-B(220)+B(225)
  JVS(854) = 0.113*B(306)+0.341*B(308)
  JVS(855) = 0.5*B(328)+0.013*B(384)
  JVS(856) = B(234)
  JVS(857) = 0
  JVS(858) = 0.379*B(303)
  JVS(859) = B(51)-B(61)+0.008*B(283)+0.064*B(292)+0.4*B(299)+1.5*B(312)+0.12*B(324)+0.033*B(340)+0.033*B(348)+0.056&
               &*B(372)+0.003*B(380)+1.5*B(386)+0.06*B(390)
  JVS(860) = -B(53)+B(78)+B(86)+B(224)
  JVS(861) = B(87)-B(88)+B(90)+B(92)+B(94)+B(114)
  JVS(862) = -B(124)
  JVS(863) = -B(54)-B(55)-B(62)-2*B(63)-2*B(64)-B(66)-B(72)-B(80)-B(89)-B(108)-B(125)-B(142)-B(162)-B(184)-B(209)-B(214)&
               &-B(221)-B(396)
  JVS(864) = -B(143)
  JVS(865) = -B(185)
  JVS(866) = B(44)-B(67)+B(82)+B(91)+B(112)+B(226)+0.63*B(258)
  JVS(867) = -B(163)
  JVS(868) = -B(56)
  JVS(869) = B(79)-B(81)+B(83)+2*B(85)+B(93)+B(110)
  JVS(870) = B(45)+B(50)+B(52)+B(71)-B(73)+B(75)+B(76)+B(219)+B(244)+0.95*B(246)+0.63*B(256)+0.379*B(304)+0.113*B(307)&
               &+0.297*B(358)+0.224*B(366)+0.187*B(368)
  JVS(871) = -B(109)+B(111)+B(113)+B(115)+B(118)
  JVS(872) = B(139)
  JVS(873) = 0.1*B(282)
  JVS(874) = 0.201*B(347)
  JVS(875) = 0.201*B(339)
  JVS(876) = 0.37*B(255)+0.37*B(257)
  JVS(877) = 0.048*B(298)+0.3*B(302)
  JVS(878) = 0.006*B(379)
  JVS(879) = 0.05*B(291)
  JVS(880) = 0
  JVS(881) = 0.965*B(232)+B(235)
  JVS(882) = 0.096*B(240)
  JVS(883) = 0.049*B(303)+0.333*B(305)
  JVS(884) = 0.1*B(283)+0.05*B(292)+0.048*B(299)+0.201*B(340)+0.201*B(348)+0.006*B(380)
  JVS(885) = -B(140)
  JVS(886) = -B(148)
  JVS(887) = -B(154)
  JVS(888) = -B(142)
  JVS(889) = -B(137)-B(141)-B(143)-B(144)-B(146)-B(149)-B(152)-B(155)-2*B(156)-B(176)-B(198)
  JVS(890) = -B(199)
  JVS(891) = -B(145)+B(236)+0.37*B(258)
  JVS(892) = -B(177)
  JVS(893) = -B(138)
  JVS(894) = -B(147)
  JVS(895) = 0.965*B(233)+0.096*B(241)+0.37*B(256)+0.049*B(304)
  JVS(896) = -B(153)
  JVS(897) = B(181)
  JVS(898) = 0.192*B(331)+0.24*B(335)
  JVS(899) = 0.5*B(280)+0.5*B(284)+0.33*B(288)
  JVS(900) = 0.289*B(296)+0.15*B(300)
  JVS(901) = 0
  JVS(902) = 0.3*B(295)
  JVS(903) = 0.24*B(336)
  JVS(904) = 0.192*B(332)
  JVS(905) = -B(182)
  JVS(906) = -B(190)
  JVS(907) = -B(196)
  JVS(908) = -B(184)
  JVS(909) = -B(198)
  JVS(910) = -B(179)-B(183)-B(185)-B(186)-B(188)-B(191)-B(194)-B(197)-B(199)-B(200)-2*B(202)
  JVS(911) = -B(187)+0.5*B(285)+0.15*B(301)
  JVS(912) = -B(201)
  JVS(913) = -B(180)
  JVS(914) = -B(189)
  JVS(915) = 0.5*B(281)+0.289*B(297)
  JVS(916) = -B(195)
  JVS(917) = B(23)
  JVS(918) = 0.39*B(58)
  JVS(919) = -B(271)
  JVS(920) = -B(273)
  JVS(921) = -B(278)
  JVS(922) = -B(267)
  JVS(923) = -B(262)
  JVS(924) = B(46)
  JVS(925) = -B(325)
  JVS(926) = -B(391)
  JVS(927) = 0
  JVS(928) = -B(333)
  JVS(929) = -B(373)
  JVS(930) = -B(99)
  JVS(931) = -B(284)
  JVS(932) = -B(349)
  JVS(933) = -B(341)
  JVS(934) = -B(257)
  JVS(935) = -B(300)
  JVS(936) = -B(230)
  JVS(937) = -B(381)
  JVS(938) = 0
  JVS(939) = -B(225)
  JVS(940) = 0
  JVS(941) = B(11)
  JVS(942) = -B(235)
  JVS(943) = 0
  JVS(944) = 0
  JVS(945) = B(15)
  JVS(946) = -B(17)
  JVS(947) = -B(90)
  JVS(948) = -B(126)
  JVS(949) = -B(66)
  JVS(950) = -B(144)
  JVS(951) = -B(186)
  JVS(952) = -B(18)-B(21)-B(26)-B(28)-B(29)-B(44)-B(67)-2*B(68)-B(82)-B(91)-B(100)-B(112)-B(127)-B(145)-B(164)-B(187)&
               &-B(226)-B(231)-B(236)-B(258)-B(263)-B(268)-B(272)-B(274)-B(279)-B(285)-B(301)-B(326)-B(334)-B(342)-B(350)&
               &-B(374)-B(382)-B(392)
  JVS(953) = -B(165)
  JVS(954) = B(12)+B(16)-B(22)-B(27)
  JVS(955) = -B(83)
  JVS(956) = -B(45)+B(47)
  JVS(957) = -B(113)
  JVS(958) = B(159)
  JVS(959) = B(275)+B(278)
  JVS(960) = 0
  JVS(961) = 0
  JVS(962) = 0
  JVS(963) = -B(160)
  JVS(964) = -B(168)
  JVS(965) = -B(174)
  JVS(966) = -B(162)
  JVS(967) = -B(176)
  JVS(968) = -B(200)
  JVS(969) = -B(164)+B(279)
  JVS(970) = -B(157)-B(161)-B(163)-B(165)-B(166)-B(169)-B(172)-B(175)-B(177)-2*B(178)-B(201)
  JVS(971) = -B(158)
  JVS(972) = -B(167)
  JVS(973) = B(276)
  JVS(974) = -B(173)
  JVS(975) = B(121)
  JVS(976) = B(139)
  JVS(977) = B(159)
  JVS(978) = B(181)
  JVS(979) = B(23)
  JVS(980) = B(39)+B(40)
  JVS(981) = -B(203)
  JVS(982) = B(223)
  JVS(983) = -B(211)
  JVS(984) = B(57)+0.61*B(58)+B(59)
  JVS(985) = 0
  JVS(986) = B(48)
  JVS(987) = -B(206)
  JVS(988) = 0.187*B(333)
  JVS(989) = B(95)+B(99)
  JVS(990) = 0
  JVS(991) = 0.474*B(349)
  JVS(992) = 0.474*B(341)
  JVS(993) = 0
  JVS(994) = 0
  JVS(995) = 0
  JVS(996) = 0.391*B(381)
  JVS(997) = 0
  JVS(998) = 0
  JVS(999) = 0.338*B(306)+B(308)
  JVS(1000) = B(6)-B(9)-B(11)
  JVS(1001) = 0
  JVS(1002) = 0
  JVS(1003) = 0
  JVS(1004) = B(13)-B(15)
  JVS(1005) = B(7)+B(14)+2*B(17)+2*B(19)+B(53)+B(78)+B(86)+B(96)+B(122)+B(140)+B(160)+B(182)+B(224)
  JVS(1006) = B(87)+B(90)
  JVS(1007) = -B(119)+B(123)+B(126)
  JVS(1008) = B(54)-B(55)+0.8*B(66)
  JVS(1009) = -B(137)+B(141)+B(144)
  JVS(1010) = -B(179)+B(183)+B(186)
  JVS(1011) = 2*B(18)-B(21)+B(29)+B(44)+0.8*B(67)+2*B(68)+B(82)+B(91)+B(100)+B(112)+B(127)+B(145)+B(164)+B(187)+0.187&
                &*B(334)+0.474*B(342)+0.474*B(350)+0.391*B(382)
  JVS(1012) = -B(157)+B(161)+B(165)
  JVS(1013) = -B(1)-B(10)-B(12)-B(16)-B(22)-B(42)-B(56)-B(120)-B(138)-B(158)-B(180)-B(204)-B(207)-B(212)
  JVS(1014) = B(79)+B(83)
  JVS(1015) = B(41)-B(43)+B(45)+B(60)+0.338*B(307)
  JVS(1016) = B(113)
  JVS(1017) = B(319)
  JVS(1018) = B(205)
  JVS(1019) = 0.65*B(247)
  JVS(1020) = 0.011*B(361)
  JVS(1021) = 0.3*B(327)
  JVS(1022) = B(239)
  JVS(1023) = 0.26*B(389)
  JVS(1024) = 0.25*B(335)
  JVS(1025) = 0.076*B(371)
  JVS(1026) = 0
  JVS(1027) = 0
  JVS(1028) = B(229)
  JVS(1029) = 0.197*B(379)+0.03*B(381)
  JVS(1030) = 0.3*B(295)
  JVS(1031) = 0
  JVS(1032) = 0.3*B(328)+0.25*B(336)
  JVS(1033) = 0
  JVS(1034) = 0
  JVS(1035) = 0
  JVS(1036) = 0.076*B(372)+0.197*B(380)+0.26*B(390)
  JVS(1037) = -B(78)+B(122)
  JVS(1038) = -B(92)
  JVS(1039) = B(123)+B(126)-B(128)+2*B(136)+B(154)+B(174)+B(196)
  JVS(1040) = -B(80)
  JVS(1041) = -B(146)+B(155)
  JVS(1042) = -B(188)+B(197)
  JVS(1043) = -B(82)+B(127)+0.03*B(382)
  JVS(1044) = -B(166)+B(175)
  JVS(1045) = 0
  JVS(1046) = -B(79)-B(81)-B(83)-2*B(84)-2*B(85)-B(93)-B(110)-B(129)-B(147)-B(167)-B(189)
  JVS(1047) = 0.65*B(248)+B(320)+0.011*B(362)
  JVS(1048) = -B(111)
  JVS(1049) = -B(462)
  JVS(1050) = -B(470)
  JVS(1051) = 2*B(32)
  JVS(1052) = -B(478)
  JVS(1053) = -B(466)
  JVS(1054) = -B(482)
  JVS(1055) = -B(474)
  JVS(1056) = -B(319)
  JVS(1057) = -B(74)
  JVS(1058) = -B(353)
  JVS(1059) = 2*B(69)-B(70)
  JVS(1060) = -B(355)
  JVS(1061) = -B(245)
  JVS(1062) = B(38)-B(40)
  JVS(1063) = -B(359)
  JVS(1064) = -B(363)
  JVS(1065) = -B(367)
  JVS(1066) = -0.65*B(247)+B(249)
  JVS(1067) = 0.39*B(58)-B(59)
  JVS(1068) = -B(243)
  JVS(1069) = -B(365)
  JVS(1070) = -B(313)
  JVS(1071) = -B(316)
  JVS(1072) = -B(269)
  JVS(1073) = -B(361)
  JVS(1074) = -B(309)+0.5*B(311)
  JVS(1075) = -0.397*B(357)+0.5*B(385)
  JVS(1076) = -0.34*B(250)+B(252)
  JVS(1077) = -B(275)
  JVS(1078) = -B(265)
  JVS(1079) = -B(260)
  JVS(1080) = -B(49)
  JVS(1081) = -B(46)+B(48)
  JVS(1082) = -B(321)+0.12*B(323)
  JVS(1083) = -B(237)
  JVS(1084) = -B(387)+0.32*B(389)
  JVS(1085) = 0
  JVS(1086) = -B(329)+0.266*B(331)
  JVS(1087) = -B(369)+0.155*B(371)
  JVS(1088) = -B(280)+0.208*B(282)+0.33*B(288)
  JVS(1089) = -B(345)+0.567*B(347)
  JVS(1090) = -B(337)+0.567*B(339)
  JVS(1091) = -B(255)
  JVS(1092) = -B(296)+0.285*B(298)
  JVS(1093) = -B(227)
  JVS(1094) = -B(377)+0.378*B(379)
  JVS(1095) = -B(289)+0.164*B(291)
  JVS(1096) = -B(218)
  JVS(1097) = -B(306)
  JVS(1098) = 0
  JVS(1099) = -B(232)
  JVS(1100) = -B(240)
  JVS(1101) = -B(303)
  JVS(1102) = -B(51)+B(61)+0.208*B(283)+0.164*B(292)+0.285*B(299)+0.5*B(312)+0.12*B(324)+0.266*B(332)+0.567*B(340)+0.567&
                &*B(348)+0.155*B(372)+0.378*B(380)+0.5*B(386)+0.32*B(390)
  JVS(1103) = -B(36)+B(53)
  JVS(1104) = 0
  JVS(1105) = 0
  JVS(1106) = B(54)+B(62)+0.8*B(66)-B(72)
  JVS(1107) = 0
  JVS(1108) = 0
  JVS(1109) = -B(44)+0.8*B(67)
  JVS(1110) = 0
  JVS(1111) = -B(42)
  JVS(1112) = 0
  JVS(1113) = -B(37)-B(41)-B(43)-B(45)-B(47)-B(50)-B(52)-B(60)-B(71)-B(73)-B(75)-B(76)-B(219)-B(228)-B(233)-B(238)&
                &-B(241)-B(244)-B(246)-0.65*B(248)-0.34*B(251)-B(256)-B(261)-B(266)-B(270)-B(276)-B(281)-B(290)-B(297)&
                &-B(304)-B(307)-B(310)-B(314)-B(317)-B(320)-B(322)-B(330)-B(338)-B(346)-B(354)-B(356)-0.397*B(358)-B(360)&
                &-B(362)-B(364)-B(366)-B(368)-B(370)-B(378)-B(388)-B(463)-B(467)-B(471)-B(475)-B(479)-B(483)
  JVS(1114) = 0
  JVS(1115) = 0.035*B(355)
  JVS(1116) = 0.07*B(359)
  JVS(1117) = 0.347*B(363)
  JVS(1118) = 0.009*B(367)
  JVS(1119) = 0.011*B(365)
  JVS(1120) = 0.143*B(361)
  JVS(1121) = 0.016*B(387)+0.051*B(391)
  JVS(1122) = 0.093*B(329)+0.008*B(331)+0.064*B(333)+0.01*B(335)
  JVS(1123) = 0.09*B(369)+0.001*B(371)+0.176*B(373)
  JVS(1124) = 0.25*B(345)+0.18*B(347)+0.25*B(349)
  JVS(1125) = 0.25*B(337)+0.18*B(339)+0.25*B(341)
  JVS(1126) = 0.041*B(296)+0.051*B(300)
  JVS(1127) = 0.082*B(377)+0.002*B(379)+0.136*B(381)+0.001*B(383)
  JVS(1128) = 0.025*B(289)
  JVS(1129) = 0.173*B(306)+0.095*B(308)
  JVS(1130) = 0.01*B(336)+0.001*B(384)
  JVS(1131) = 0.001*B(232)
  JVS(1132) = 0.042*B(240)
  JVS(1133) = 0.07*B(303)+0.04*B(305)
  JVS(1134) = 0.008*B(332)+0.18*B(340)+0.18*B(348)+0.001*B(372)+0.002*B(380)
  JVS(1135) = -B(106)
  JVS(1136) = -B(114)
  JVS(1137) = -B(134)
  JVS(1138) = -B(108)
  JVS(1139) = -B(152)
  JVS(1140) = -B(194)
  JVS(1141) = -B(112)+0.051*B(301)+0.064*B(334)+0.25*B(342)+0.25*B(350)+0.176*B(374)+0.136*B(382)+0.051*B(392)
  JVS(1142) = -B(172)
  JVS(1143) = 0
  JVS(1144) = -B(110)
  JVS(1145) = 0.001*B(233)+0.042*B(241)+0.025*B(290)+0.041*B(297)+0.07*B(304)+0.173*B(307)+0.093*B(330)+0.25*B(338)+0.25&
                &*B(346)+0.035*B(356)+0.07*B(360)+0.143*B(362)+0.347*B(364)+0.011*B(366)+0.009*B(368)+0.09*B(370)+0.082&
                &*B(378)+0.016*B(388)
  JVS(1146) = -B(107)-B(109)-B(111)-B(113)-B(115)-2*B(118)-B(135)-B(153)-B(173)-B(195)
END SUBROUTINE saprc99_mosaic_8bin_vbs2_Jac_SP
SUBROUTINE saprc99_mosaic_8bin_vbs2_KppDecomp( JVS, IER )
      INTEGER :: IER
      REAL(kind=dp) :: JVS(1146), W(104), a
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
END SUBROUTINE saprc99_mosaic_8bin_vbs2_KppDecomp
SUBROUTINE saprc99_mosaic_8bin_vbs2_KppDecompCmplx( JVS, IER )
      INTEGER :: IER
      DOUBLE COMPLEX :: JVS(1146), W(104), a
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
END SUBROUTINE saprc99_mosaic_8bin_vbs2_KppDecompCmplx
SUBROUTINE saprc99_mosaic_8bin_vbs2_KppSolveIndirect( JVS, X )
      INTEGER i, j
      REAL(kind=dp) JVS(1146), X(104), sum
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
END SUBROUTINE saprc99_mosaic_8bin_vbs2_KppSolveIndirect
SUBROUTINE saprc99_mosaic_8bin_vbs2_KppSolveCmplx( JVS, X )
      INTEGER i, j
      DOUBLE COMPLEX JVS(1146), X(104), sum
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
END SUBROUTINE saprc99_mosaic_8bin_vbs2_KppSolveCmplx
SUBROUTINE saprc99_mosaic_8bin_vbs2_KppSolve ( JVS, X )
  REAL(kind=dp) :: JVS(LU_NONZERO)
  REAL(kind=dp) :: X(NVAR)
  X(25) = X(25)-JVS(177)*X(3)-JVS(178)*X(4)-JVS(179)*X(5)-JVS(180)*X(6)
  X(32) = X(32)-JVS(207)*X(31)
  X(34) = X(34)-JVS(212)*X(33)
  X(49) = X(49)-JVS(254)*X(48)
  X(58) = X(58)-JVS(286)*X(51)-JVS(287)*X(57)
  X(59) = X(59)-JVS(290)*X(51)-JVS(291)*X(57)
  X(60) = X(60)-JVS(294)*X(51)-JVS(295)*X(57)
  X(62) = X(62)-JVS(301)*X(51)-JVS(302)*X(57)
  X(64) = X(64)-JVS(309)*X(54)
  X(66) = X(66)-JVS(320)*X(51)-JVS(321)*X(57)
  X(67) = X(67)-JVS(327)*X(57)
  X(68) = X(68)-JVS(333)*X(51)-JVS(334)*X(57)-JVS(335)*X(58)-JVS(336)*X(59)-JVS(337)*X(60)
  X(69) = X(69)-JVS(345)*X(58)-JVS(346)*X(59)-JVS(347)*X(61)-JVS(348)*X(62)-JVS(349)*X(63)-JVS(350)*X(68)
  X(70) = X(70)-JVS(370)*X(46)-JVS(371)*X(60)-JVS(372)*X(64)-JVS(373)*X(66)-JVS(374)*X(67)-JVS(375)*X(68)
  X(72) = X(72)-JVS(396)*X(43)-JVS(397)*X(48)-JVS(398)*X(49)-JVS(399)*X(50)-JVS(400)*X(61)
  X(74) = X(74)-JVS(416)*X(60)-JVS(417)*X(67)
  X(77) = X(77)-JVS(438)*X(48)-JVS(439)*X(50)-JVS(440)*X(58)-JVS(441)*X(59)-JVS(442)*X(61)-JVS(443)*X(72)-JVS(444)*X(75)&
            &-JVS(445)*X(76)
  X(78) = X(78)-JVS(468)*X(75)
  X(81) = X(81)-JVS(485)*X(51)-JVS(486)*X(57)-JVS(487)*X(58)-JVS(488)*X(59)-JVS(489)*X(62)-JVS(490)*X(63)-JVS(491)*X(67)&
            &-JVS(492)*X(71)-JVS(493)*X(74)-JVS(494)*X(79)-JVS(495)*X(80)
  X(82) = X(82)-JVS(509)*X(75)
  X(83) = X(83)-JVS(516)*X(37)-JVS(517)*X(44)-JVS(518)*X(48)-JVS(519)*X(50)-JVS(520)*X(61)-JVS(521)*X(71)-JVS(522)*X(73)&
            &-JVS(523)*X(76)-JVS(524)*X(82)
  X(85) = X(85)-JVS(545)*X(75)-JVS(546)*X(84)
  X(86) = X(86)-JVS(552)*X(44)-JVS(553)*X(48)-JVS(554)*X(50)-JVS(555)*X(52)-JVS(556)*X(53)-JVS(557)*X(56)-JVS(558)*X(61)&
            &-JVS(559)*X(63)-JVS(560)*X(71)-JVS(561)*X(72)-JVS(562)*X(73)-JVS(563)*X(75)-JVS(564)*X(76)-JVS(565)*X(78)&
            &-JVS(566)*X(79)-JVS(567)*X(80)-JVS(568)*X(81)-JVS(569)*X(82)-JVS(570)*X(84)-JVS(571)*X(85)
  X(87) = X(87)-JVS(590)*X(49)-JVS(591)*X(76)-JVS(592)*X(79)-JVS(593)*X(80)-JVS(594)*X(82)-JVS(595)*X(84)
  X(88) = X(88)-JVS(604)*X(30)-JVS(605)*X(71)-JVS(606)*X(73)-JVS(607)*X(75)-JVS(608)*X(76)-JVS(609)*X(78)-JVS(610)*X(79)&
            &-JVS(611)*X(80)-JVS(612)*X(84)-JVS(613)*X(85)
  X(89) = X(89)-JVS(620)*X(43)-JVS(621)*X(48)-JVS(622)*X(50)-JVS(623)*X(58)-JVS(624)*X(59)-JVS(625)*X(61)-JVS(626)*X(62)&
            &-JVS(627)*X(65)-JVS(628)*X(71)-JVS(629)*X(73)-JVS(630)*X(76)-JVS(631)*X(78)-JVS(632)*X(79)-JVS(633)*X(80)&
            &-JVS(634)*X(82)-JVS(635)*X(84)-JVS(636)*X(85)-JVS(637)*X(87)-JVS(638)*X(88)
  X(90) = X(90)-JVS(650)*X(48)-JVS(651)*X(50)-JVS(652)*X(61)-JVS(653)*X(73)-JVS(654)*X(76)-JVS(655)*X(78)-JVS(656)*X(82)&
            &-JVS(657)*X(84)-JVS(658)*X(85)-JVS(659)*X(87)-JVS(660)*X(88)
  X(91) = X(91)-JVS(671)*X(50)-JVS(672)*X(57)-JVS(673)*X(61)-JVS(674)*X(75)-JVS(675)*X(76)-JVS(676)*X(79)-JVS(677)*X(80)&
            &-JVS(678)*X(82)-JVS(679)*X(84)-JVS(680)*X(85)-JVS(681)*X(87)-JVS(682)*X(88)
  X(92) = X(92)-JVS(695)*X(62)-JVS(696)*X(63)-JVS(697)*X(71)-JVS(698)*X(73)-JVS(699)*X(75)-JVS(700)*X(76)-JVS(701)*X(78)&
            &-JVS(702)*X(79)-JVS(703)*X(80)-JVS(704)*X(82)-JVS(705)*X(84)-JVS(706)*X(85)-JVS(707)*X(88)
  X(93) = X(93)-JVS(718)*X(47)-JVS(719)*X(53)-JVS(720)*X(77)-JVS(721)*X(79)-JVS(722)*X(80)-JVS(723)*X(84)-JVS(724)*X(85)&
            &-JVS(725)*X(86)-JVS(726)*X(87)-JVS(727)*X(88)-JVS(728)*X(90)-JVS(729)*X(91)-JVS(730)*X(92)
  X(94) = X(94)-JVS(743)*X(37)-JVS(744)*X(43)-JVS(745)*X(44)-JVS(746)*X(48)-JVS(747)*X(50)-JVS(748)*X(51)-JVS(749)*X(57)&
            &-JVS(750)*X(58)-JVS(751)*X(59)-JVS(752)*X(60)-JVS(753)*X(61)-JVS(754)*X(62)-JVS(755)*X(63)-JVS(756)*X(65)&
            &-JVS(757)*X(67)-JVS(758)*X(71)-JVS(759)*X(73)-JVS(760)*X(74)-JVS(761)*X(75)-JVS(762)*X(76)-JVS(763)*X(78)&
            &-JVS(764)*X(79)-JVS(765)*X(80)-JVS(766)*X(82)-JVS(767)*X(84)-JVS(768)*X(85)-JVS(769)*X(87)-JVS(770)*X(88)&
            &-JVS(771)*X(89)-JVS(772)*X(90)-JVS(773)*X(91)-JVS(774)*X(92)-JVS(775)*X(93)
  X(95) = X(95)-JVS(787)*X(38)-JVS(788)*X(45)-JVS(789)*X(51)-JVS(790)*X(58)-JVS(791)*X(59)-JVS(792)*X(61)-JVS(793)*X(68)&
            &-JVS(794)*X(72)-JVS(795)*X(76)-JVS(796)*X(78)-JVS(797)*X(79)-JVS(798)*X(80)-JVS(799)*X(82)-JVS(800)*X(83)&
            &-JVS(801)*X(84)-JVS(802)*X(85)-JVS(803)*X(87)-JVS(804)*X(88)-JVS(805)*X(89)-JVS(806)*X(90)-JVS(807)*X(91)&
            &-JVS(808)*X(92)-JVS(809)*X(93)-JVS(810)*X(94)
  X(96) = X(96)-JVS(821)*X(36)-JVS(822)*X(42)-JVS(823)*X(44)-JVS(824)*X(47)-JVS(825)*X(51)-JVS(826)*X(52)-JVS(827)*X(53)&
            &-JVS(828)*X(54)-JVS(829)*X(55)-JVS(830)*X(56)-JVS(831)*X(57)-JVS(832)*X(58)-JVS(833)*X(59)-JVS(834)*X(62)&
            &-JVS(835)*X(63)-JVS(836)*X(64)-JVS(837)*X(65)-JVS(838)*X(68)-JVS(839)*X(69)-JVS(840)*X(71)-JVS(841)*X(73)&
            &-JVS(842)*X(74)-JVS(843)*X(75)-JVS(844)*X(76)-JVS(845)*X(78)-JVS(846)*X(79)-JVS(847)*X(80)-JVS(848)*X(81)&
            &-JVS(849)*X(82)-JVS(850)*X(83)-JVS(851)*X(84)-JVS(852)*X(85)-JVS(853)*X(86)-JVS(854)*X(87)-JVS(855)*X(88)&
            &-JVS(856)*X(89)-JVS(857)*X(90)-JVS(858)*X(91)-JVS(859)*X(92)-JVS(860)*X(93)-JVS(861)*X(94)-JVS(862)*X(95)
  X(97) = X(97)-JVS(872)*X(39)-JVS(873)*X(78)-JVS(874)*X(79)-JVS(875)*X(80)-JVS(876)*X(81)-JVS(877)*X(82)-JVS(878)*X(84)&
            &-JVS(879)*X(85)-JVS(880)*X(88)-JVS(881)*X(89)-JVS(882)*X(90)-JVS(883)*X(91)-JVS(884)*X(92)-JVS(885)*X(93)&
            &-JVS(886)*X(94)-JVS(887)*X(95)-JVS(888)*X(96)
  X(98) = X(98)-JVS(897)*X(41)-JVS(898)*X(75)-JVS(899)*X(78)-JVS(900)*X(82)-JVS(901)*X(84)-JVS(902)*X(85)-JVS(903)*X(88)&
            &-JVS(904)*X(92)-JVS(905)*X(93)-JVS(906)*X(94)-JVS(907)*X(95)-JVS(908)*X(96)-JVS(909)*X(97)
  X(99) = X(99)-JVS(917)*X(46)-JVS(918)*X(55)-JVS(919)*X(60)-JVS(920)*X(64)-JVS(921)*X(66)-JVS(922)*X(67)-JVS(923)*X(68)&
            &-JVS(924)*X(70)-JVS(925)*X(71)-JVS(926)*X(73)-JVS(927)*X(74)-JVS(928)*X(75)-JVS(929)*X(76)-JVS(930)*X(77)&
            &-JVS(931)*X(78)-JVS(932)*X(79)-JVS(933)*X(80)-JVS(934)*X(81)-JVS(935)*X(82)-JVS(936)*X(83)-JVS(937)*X(84)&
            &-JVS(938)*X(85)-JVS(939)*X(86)-JVS(940)*X(87)-JVS(941)*X(88)-JVS(942)*X(89)-JVS(943)*X(90)-JVS(944)*X(91)&
            &-JVS(945)*X(92)-JVS(946)*X(93)-JVS(947)*X(94)-JVS(948)*X(95)-JVS(949)*X(96)-JVS(950)*X(97)-JVS(951)*X(98)
  X(100) = X(100)-JVS(958)*X(40)-JVS(959)*X(66)-JVS(960)*X(84)-JVS(961)*X(88)-JVS(962)*X(92)-JVS(963)*X(93)-JVS(964)&
             &*X(94)-JVS(965)*X(95)-JVS(966)*X(96)-JVS(967)*X(97)-JVS(968)*X(98)-JVS(969)*X(99)
  X(101) = X(101)-JVS(975)*X(38)-JVS(976)*X(39)-JVS(977)*X(40)-JVS(978)*X(41)-JVS(979)*X(46)-JVS(980)*X(47)-JVS(981)&
             &*X(49)-JVS(982)*X(53)-JVS(983)*X(54)-JVS(984)*X(55)-JVS(985)*X(64)-JVS(986)*X(70)-JVS(987)*X(74)-JVS(988)&
             &*X(75)-JVS(989)*X(77)-JVS(990)*X(78)-JVS(991)*X(79)-JVS(992)*X(80)-JVS(993)*X(81)-JVS(994)*X(82)-JVS(995)&
             &*X(83)-JVS(996)*X(84)-JVS(997)*X(85)-JVS(998)*X(86)-JVS(999)*X(87)-JVS(1000)*X(88)-JVS(1001)*X(89)-JVS(1002)&
             &*X(90)-JVS(1003)*X(91)-JVS(1004)*X(92)-JVS(1005)*X(93)-JVS(1006)*X(94)-JVS(1007)*X(95)-JVS(1008)*X(96)&
             &-JVS(1009)*X(97)-JVS(1010)*X(98)-JVS(1011)*X(99)-JVS(1012)*X(100)
  X(102) = X(102)-JVS(1017)*X(35)-JVS(1018)*X(49)-JVS(1019)*X(52)-JVS(1020)*X(61)-JVS(1021)*X(71)-JVS(1022)*X(72)&
             &-JVS(1023)*X(73)-JVS(1024)*X(75)-JVS(1025)*X(76)-JVS(1026)*X(79)-JVS(1027)*X(80)-JVS(1028)*X(83)-JVS(1029)&
             &*X(84)-JVS(1030)*X(85)-JVS(1031)*X(87)-JVS(1032)*X(88)-JVS(1033)*X(89)-JVS(1034)*X(90)-JVS(1035)*X(91)&
             &-JVS(1036)*X(92)-JVS(1037)*X(93)-JVS(1038)*X(94)-JVS(1039)*X(95)-JVS(1040)*X(96)-JVS(1041)*X(97)-JVS(1042)&
             &*X(98)-JVS(1043)*X(99)-JVS(1044)*X(100)-JVS(1045)*X(101)
  X(103) = X(103)-JVS(1049)*X(3)-JVS(1050)*X(5)-JVS(1051)*X(30)-JVS(1052)*X(31)-JVS(1053)*X(32)-JVS(1054)*X(33)&
             &-JVS(1055)*X(34)-JVS(1056)*X(35)-JVS(1057)*X(36)-JVS(1058)*X(37)-JVS(1059)*X(42)-JVS(1060)*X(43)-JVS(1061)&
             &*X(44)-JVS(1062)*X(47)-JVS(1063)*X(48)-JVS(1064)*X(50)-JVS(1065)*X(51)-JVS(1066)*X(52)-JVS(1067)*X(55)&
             &-JVS(1068)*X(56)-JVS(1069)*X(57)-JVS(1070)*X(58)-JVS(1071)*X(59)-JVS(1072)*X(60)-JVS(1073)*X(61)-JVS(1074)&
             &*X(62)-JVS(1075)*X(63)-JVS(1076)*X(65)-JVS(1077)*X(66)-JVS(1078)*X(67)-JVS(1079)*X(68)-JVS(1080)*X(69)&
             &-JVS(1081)*X(70)-JVS(1082)*X(71)-JVS(1083)*X(72)-JVS(1084)*X(73)-JVS(1085)*X(74)-JVS(1086)*X(75)-JVS(1087)&
             &*X(76)-JVS(1088)*X(78)-JVS(1089)*X(79)-JVS(1090)*X(80)-JVS(1091)*X(81)-JVS(1092)*X(82)-JVS(1093)*X(83)&
             &-JVS(1094)*X(84)-JVS(1095)*X(85)-JVS(1096)*X(86)-JVS(1097)*X(87)-JVS(1098)*X(88)-JVS(1099)*X(89)-JVS(1100)&
             &*X(90)-JVS(1101)*X(91)-JVS(1102)*X(92)-JVS(1103)*X(93)-JVS(1104)*X(94)-JVS(1105)*X(95)-JVS(1106)*X(96)&
             &-JVS(1107)*X(97)-JVS(1108)*X(98)-JVS(1109)*X(99)-JVS(1110)*X(100)-JVS(1111)*X(101)-JVS(1112)*X(102)
  X(104) = X(104)-JVS(1115)*X(43)-JVS(1116)*X(48)-JVS(1117)*X(50)-JVS(1118)*X(51)-JVS(1119)*X(57)-JVS(1120)*X(61)&
             &-JVS(1121)*X(73)-JVS(1122)*X(75)-JVS(1123)*X(76)-JVS(1124)*X(79)-JVS(1125)*X(80)-JVS(1126)*X(82)-JVS(1127)&
             &*X(84)-JVS(1128)*X(85)-JVS(1129)*X(87)-JVS(1130)*X(88)-JVS(1131)*X(89)-JVS(1132)*X(90)-JVS(1133)*X(91)&
             &-JVS(1134)*X(92)-JVS(1135)*X(93)-JVS(1136)*X(94)-JVS(1137)*X(95)-JVS(1138)*X(96)-JVS(1139)*X(97)-JVS(1140)&
             &*X(98)-JVS(1141)*X(99)-JVS(1142)*X(100)-JVS(1143)*X(101)-JVS(1144)*X(102)-JVS(1145)*X(103)
  X(104) = X(104)/JVS(1146)
  X(103) = (X(103)-JVS(1114)*X(104))/(JVS(1113))
  X(102) = (X(102)-JVS(1047)*X(103)-JVS(1048)*X(104))/(JVS(1046))
  X(101) = (X(101)-JVS(1014)*X(102)-JVS(1015)*X(103)-JVS(1016)*X(104))/(JVS(1013))
  X(100) = (X(100)-JVS(971)*X(101)-JVS(972)*X(102)-JVS(973)*X(103)-JVS(974)*X(104))/(JVS(970))
  X(99) = (X(99)-JVS(953)*X(100)-JVS(954)*X(101)-JVS(955)*X(102)-JVS(956)*X(103)-JVS(957)*X(104))/(JVS(952))
  X(98) = (X(98)-JVS(911)*X(99)-JVS(912)*X(100)-JVS(913)*X(101)-JVS(914)*X(102)-JVS(915)*X(103)-JVS(916)*X(104))&
            &/(JVS(910))
  X(97) = (X(97)-JVS(890)*X(98)-JVS(891)*X(99)-JVS(892)*X(100)-JVS(893)*X(101)-JVS(894)*X(102)-JVS(895)*X(103)-JVS(896)&
            &*X(104))/(JVS(889))
  X(96) = (X(96)-JVS(864)*X(97)-JVS(865)*X(98)-JVS(866)*X(99)-JVS(867)*X(100)-JVS(868)*X(101)-JVS(869)*X(102)-JVS(870)&
            &*X(103)-JVS(871)*X(104))/(JVS(863))
  X(95) = (X(95)-JVS(812)*X(96)-JVS(813)*X(97)-JVS(814)*X(98)-JVS(815)*X(99)-JVS(816)*X(100)-JVS(817)*X(101)-JVS(818)&
            &*X(102)-JVS(819)*X(103)-JVS(820)*X(104))/(JVS(811))
  X(94) = (X(94)-JVS(777)*X(95)-JVS(778)*X(96)-JVS(779)*X(97)-JVS(780)*X(98)-JVS(781)*X(99)-JVS(782)*X(100)-JVS(783)&
            &*X(101)-JVS(784)*X(102)-JVS(785)*X(103)-JVS(786)*X(104))/(JVS(776))
  X(93) = (X(93)-JVS(732)*X(94)-JVS(733)*X(95)-JVS(734)*X(96)-JVS(735)*X(97)-JVS(736)*X(98)-JVS(737)*X(99)-JVS(738)&
            &*X(100)-JVS(739)*X(101)-JVS(740)*X(102)-JVS(741)*X(103)-JVS(742)*X(104))/(JVS(731))
  X(92) = (X(92)-JVS(709)*X(93)-JVS(710)*X(95)-JVS(711)*X(96)-JVS(712)*X(97)-JVS(713)*X(98)-JVS(714)*X(99)-JVS(715)&
            &*X(100)-JVS(716)*X(101)-JVS(717)*X(103))/(JVS(708))
  X(91) = (X(91)-JVS(684)*X(92)-JVS(685)*X(93)-JVS(686)*X(94)-JVS(687)*X(95)-JVS(688)*X(97)-JVS(689)*X(99)-JVS(690)&
            &*X(100)-JVS(691)*X(101)-JVS(692)*X(102)-JVS(693)*X(103)-JVS(694)*X(104))/(JVS(683))
  X(90) = (X(90)-JVS(662)*X(91)-JVS(663)*X(92)-JVS(664)*X(93)-JVS(665)*X(94)-JVS(666)*X(99)-JVS(667)*X(101)-JVS(668)&
            &*X(102)-JVS(669)*X(103)-JVS(670)*X(104))/(JVS(661))
  X(89) = (X(89)-JVS(640)*X(90)-JVS(641)*X(91)-JVS(642)*X(92)-JVS(643)*X(93)-JVS(644)*X(94)-JVS(645)*X(96)-JVS(646)&
            &*X(99)-JVS(647)*X(101)-JVS(648)*X(103)-JVS(649)*X(104))/(JVS(639))
  X(88) = (X(88)-JVS(615)*X(92)-JVS(616)*X(93)-JVS(617)*X(99)-JVS(618)*X(101)-JVS(619)*X(103))/(JVS(614))
  X(87) = (X(87)-JVS(597)*X(88)-JVS(598)*X(92)-JVS(599)*X(93)-JVS(600)*X(99)-JVS(601)*X(101)-JVS(602)*X(103)-JVS(603)&
            &*X(104))/(JVS(596))
  X(86) = (X(86)-JVS(573)*X(87)-JVS(574)*X(88)-JVS(575)*X(90)-JVS(576)*X(91)-JVS(577)*X(92)-JVS(578)*X(93)-JVS(579)&
            &*X(94)-JVS(580)*X(95)-JVS(581)*X(96)-JVS(582)*X(97)-JVS(583)*X(98)-JVS(584)*X(99)-JVS(585)*X(100)-JVS(586)&
            &*X(101)-JVS(587)*X(102)-JVS(588)*X(103)-JVS(589)*X(104))/(JVS(572))
  X(85) = (X(85)-JVS(548)*X(88)-JVS(549)*X(92)-JVS(550)*X(99)-JVS(551)*X(103))/(JVS(547))
  X(84) = (X(84)-JVS(541)*X(88)-JVS(542)*X(92)-JVS(543)*X(99)-JVS(544)*X(103))/(JVS(540))
  X(83) = (X(83)-JVS(526)*X(84)-JVS(527)*X(87)-JVS(528)*X(88)-JVS(529)*X(89)-JVS(530)*X(90)-JVS(531)*X(91)-JVS(532)&
            &*X(92)-JVS(533)*X(93)-JVS(534)*X(95)-JVS(535)*X(97)-JVS(536)*X(98)-JVS(537)*X(99)-JVS(538)*X(100)-JVS(539)&
            &*X(103))/(JVS(525))
  X(82) = (X(82)-JVS(511)*X(84)-JVS(512)*X(88)-JVS(513)*X(92)-JVS(514)*X(99)-JVS(515)*X(103))/(JVS(510))
  X(81) = (X(81)-JVS(497)*X(82)-JVS(498)*X(88)-JVS(499)*X(92)-JVS(500)*X(93)-JVS(501)*X(95)-JVS(502)*X(96)-JVS(503)&
            &*X(97)-JVS(504)*X(98)-JVS(505)*X(99)-JVS(506)*X(100)-JVS(507)*X(101)-JVS(508)*X(103))/(JVS(496))
  X(80) = (X(80)-JVS(481)*X(88)-JVS(482)*X(92)-JVS(483)*X(99)-JVS(484)*X(103))/(JVS(480))
  X(79) = (X(79)-JVS(476)*X(88)-JVS(477)*X(92)-JVS(478)*X(99)-JVS(479)*X(103))/(JVS(475))
  X(78) = (X(78)-JVS(470)*X(84)-JVS(471)*X(88)-JVS(472)*X(92)-JVS(473)*X(99)-JVS(474)*X(103))/(JVS(469))
  X(77) = (X(77)-JVS(447)*X(79)-JVS(448)*X(80)-JVS(449)*X(84)-JVS(450)*X(85)-JVS(451)*X(87)-JVS(452)*X(88)-JVS(453)&
            &*X(90)-JVS(454)*X(91)-JVS(455)*X(92)-JVS(456)*X(93)-JVS(457)*X(94)-JVS(458)*X(95)-JVS(459)*X(96)-JVS(460)*X(97)&
            &-JVS(461)*X(98)-JVS(462)*X(99)-JVS(463)*X(100)-JVS(464)*X(101)-JVS(465)*X(102)-JVS(466)*X(103)-JVS(467)*X(104))&
            &/(JVS(446))
  X(76) = (X(76)-JVS(434)*X(88)-JVS(435)*X(92)-JVS(436)*X(99)-JVS(437)*X(103))/(JVS(433))
  X(75) = (X(75)-JVS(429)*X(88)-JVS(430)*X(92)-JVS(431)*X(99)-JVS(432)*X(103))/(JVS(428))
  X(74) = (X(74)-JVS(419)*X(93)-JVS(420)*X(95)-JVS(421)*X(96)-JVS(422)*X(97)-JVS(423)*X(98)-JVS(424)*X(99)-JVS(425)&
            &*X(100)-JVS(426)*X(101)-JVS(427)*X(103))/(JVS(418))
  X(73) = (X(73)-JVS(412)*X(88)-JVS(413)*X(92)-JVS(414)*X(99)-JVS(415)*X(103))/(JVS(411))
  X(72) = (X(72)-JVS(402)*X(76)-JVS(403)*X(79)-JVS(404)*X(80)-JVS(405)*X(84)-JVS(406)*X(87)-JVS(407)*X(92)-JVS(408)&
            &*X(99)-JVS(409)*X(101)-JVS(410)*X(103))/(JVS(401))
  X(71) = (X(71)-JVS(392)*X(88)-JVS(393)*X(92)-JVS(394)*X(99)-JVS(395)*X(103))/(JVS(391))
  X(70) = (X(70)-JVS(377)*X(74)-JVS(378)*X(78)-JVS(379)*X(81)-JVS(380)*X(82)-JVS(381)*X(83)-JVS(382)*X(84)-JVS(383)&
            &*X(85)-JVS(384)*X(86)-JVS(385)*X(89)-JVS(386)*X(92)-JVS(387)*X(96)-JVS(388)*X(99)-JVS(389)*X(101)-JVS(390)&
            &*X(103))/(JVS(376))
  X(69) = (X(69)-JVS(352)*X(71)-JVS(353)*X(73)-JVS(354)*X(75)-JVS(355)*X(76)-JVS(356)*X(78)-JVS(357)*X(79)-JVS(358)&
            &*X(80)-JVS(359)*X(81)-JVS(360)*X(82)-JVS(361)*X(83)-JVS(362)*X(84)-JVS(363)*X(85)-JVS(364)*X(86)-JVS(365)*X(88)&
            &-JVS(366)*X(89)-JVS(367)*X(92)-JVS(368)*X(99)-JVS(369)*X(103))/(JVS(351))
  X(68) = (X(68)-JVS(339)*X(78)-JVS(340)*X(82)-JVS(341)*X(85)-JVS(342)*X(92)-JVS(343)*X(99)-JVS(344)*X(103))/(JVS(338))
  X(67) = (X(67)-JVS(329)*X(74)-JVS(330)*X(96)-JVS(331)*X(99)-JVS(332)*X(103))/(JVS(328))
  X(66) = (X(66)-JVS(323)*X(84)-JVS(324)*X(92)-JVS(325)*X(99)-JVS(326)*X(103))/(JVS(322))
  X(65) = (X(65)-JVS(316)*X(94)-JVS(317)*X(96)-JVS(318)*X(103)-JVS(319)*X(104))/(JVS(315))
  X(64) = (X(64)-JVS(311)*X(74)-JVS(312)*X(96)-JVS(313)*X(99)-JVS(314)*X(101))/(JVS(310))
  X(63) = (X(63)-JVS(307)*X(92)-JVS(308)*X(103))/(JVS(306))
  X(62) = (X(62)-JVS(304)*X(92)-JVS(305)*X(103))/(JVS(303))
  X(61) = (X(61)-JVS(300)*X(103))/(JVS(299))
  X(60) = (X(60)-JVS(297)*X(99)-JVS(298)*X(103))/(JVS(296))
  X(59) = (X(59)-JVS(293)*X(103))/(JVS(292))
  X(58) = (X(58)-JVS(289)*X(103))/(JVS(288))
  X(57) = (X(57)-JVS(285)*X(103))/(JVS(284))
  X(56) = (X(56)-JVS(280)*X(94)-JVS(281)*X(102)-JVS(282)*X(103)-JVS(283)*X(104))/(JVS(279))
  X(55) = (X(55)-JVS(276)*X(96)-JVS(277)*X(101)-JVS(278)*X(103))/(JVS(275))
  X(54) = (X(54)-JVS(271)*X(64)-JVS(272)*X(96)-JVS(273)*X(99)-JVS(274)*X(101))/(JVS(270))
  X(53) = (X(53)-JVS(267)*X(86)-JVS(268)*X(93)-JVS(269)*X(96))/(JVS(266))
  X(52) = (X(52)-JVS(263)*X(96)-JVS(264)*X(102)-JVS(265)*X(103))/(JVS(262))
  X(51) = (X(51)-JVS(261)*X(103))/(JVS(260))
  X(50) = (X(50)-JVS(259)*X(103))/(JVS(258))
  X(49) = (X(49)-JVS(256)*X(101)-JVS(257)*X(103))/(JVS(255))
  X(48) = (X(48)-JVS(253)*X(103))/(JVS(252))
  X(47) = (X(47)-JVS(250)*X(93)-JVS(251)*X(103))/(JVS(249))
  X(46) = (X(46)-JVS(247)*X(99)-JVS(248)*X(101))/(JVS(246))
  X(45) = (X(45)-JVS(241)*X(51)-JVS(242)*X(79)-JVS(243)*X(80)-JVS(244)*X(92)-JVS(245)*X(103))/(JVS(240))
  X(44) = (X(44)-JVS(239)*X(103))/(JVS(238))
  X(43) = (X(43)-JVS(237)*X(103))/(JVS(236))
  X(42) = (X(42)-JVS(234)*X(96)-JVS(235)*X(103))/(JVS(233))
  X(41) = (X(41)-JVS(231)*X(98)-JVS(232)*X(101))/(JVS(230))
  X(40) = (X(40)-JVS(228)*X(100)-JVS(229)*X(101))/(JVS(227))
  X(39) = (X(39)-JVS(225)*X(97)-JVS(226)*X(101))/(JVS(224))
  X(38) = (X(38)-JVS(222)*X(95)-JVS(223)*X(101))/(JVS(221))
  X(37) = (X(37)-JVS(220)*X(103))/(JVS(219))
  X(36) = (X(36)-JVS(218)*X(103))/(JVS(217))
  X(35) = (X(35)-JVS(216)*X(103))/(JVS(215))
  X(34) = (X(34)-JVS(214)*X(103))/(JVS(213))
  X(33) = (X(33)-JVS(211)*X(103))/(JVS(210))
  X(32) = (X(32)-JVS(209)*X(103))/(JVS(208))
  X(31) = (X(31)-JVS(206)*X(103))/(JVS(205))
  X(30) = (X(30)-JVS(204)*X(92))/(JVS(203))
  X(29) = (X(29)-JVS(202)*X(103))/(JVS(201))
  X(28) = (X(28)-JVS(198)*X(29)-JVS(199)*X(33)-JVS(200)*X(103))/(JVS(197))
  X(27) = (X(27)-JVS(196)*X(103))/(JVS(195))
  X(26) = (X(26)-JVS(192)*X(27)-JVS(193)*X(31)-JVS(194)*X(103))/(JVS(191))
  X(25) = (X(25)-JVS(182)*X(26)-JVS(183)*X(27)-JVS(184)*X(28)-JVS(185)*X(29)-JVS(186)*X(31)-JVS(187)*X(32)-JVS(188)&
            &*X(33)-JVS(189)*X(34)-JVS(190)*X(103))/(JVS(181))
  X(24) = (X(24)-JVS(129)*X(35)-JVS(130)*X(36)-JVS(131)*X(37)-JVS(132)*X(43)-JVS(133)*X(44)-JVS(134)*X(47)-JVS(135)&
            &*X(48)-JVS(136)*X(50)-JVS(137)*X(51)-JVS(138)*X(52)-JVS(139)*X(56)-JVS(140)*X(57)-JVS(141)*X(58)-JVS(142)*X(59)&
            &-JVS(143)*X(60)-JVS(144)*X(61)-JVS(145)*X(62)-JVS(146)*X(63)-JVS(147)*X(65)-JVS(148)*X(66)-JVS(149)*X(67)&
            &-JVS(150)*X(68)-JVS(151)*X(69)-JVS(152)*X(70)-JVS(153)*X(71)-JVS(154)*X(72)-JVS(155)*X(73)-JVS(156)*X(75)&
            &-JVS(157)*X(76)-JVS(158)*X(78)-JVS(159)*X(79)-JVS(160)*X(80)-JVS(161)*X(81)-JVS(162)*X(82)-JVS(163)*X(83)&
            &-JVS(164)*X(84)-JVS(165)*X(85)-JVS(166)*X(86)-JVS(167)*X(87)-JVS(168)*X(89)-JVS(169)*X(90)-JVS(170)*X(91)&
            &-JVS(171)*X(92)-JVS(172)*X(93)-JVS(173)*X(96)-JVS(174)*X(99)-JVS(175)*X(101)-JVS(176)*X(103))/(JVS(128))
  X(23) = (X(23)-JVS(122)*X(75)-JVS(123)*X(79)-JVS(124)*X(80)-JVS(125)*X(92)-JVS(126)*X(99)-JVS(127)*X(103))/(JVS(121))
  X(22) = (X(22)-JVS(115)*X(75)-JVS(116)*X(79)-JVS(117)*X(80)-JVS(118)*X(92)-JVS(119)*X(99)-JVS(120)*X(103))/(JVS(114))
  X(21) = (X(21)-JVS(108)*X(75)-JVS(109)*X(79)-JVS(110)*X(80)-JVS(111)*X(92)-JVS(112)*X(99)-JVS(113)*X(103))/(JVS(107))
  X(20) = (X(20)-JVS(105)*X(75)-JVS(106)*X(103))/(JVS(104))
  X(19) = (X(19)-JVS(101)*X(51)-JVS(102)*X(57)-JVS(103)*X(103))/(JVS(100))
  X(18) = (X(18)-JVS(97)*X(51)-JVS(98)*X(57)-JVS(99)*X(103))/(JVS(96))
  X(17) = (X(17)-JVS(93)*X(51)-JVS(94)*X(57)-JVS(95)*X(103))/(JVS(92))
  X(16) = (X(16)-JVS(89)*X(51)-JVS(90)*X(57)-JVS(91)*X(103))/(JVS(88))
  X(15) = (X(15)-JVS(83)*X(77)-JVS(84)*X(94)-JVS(85)*X(96)-JVS(86)*X(102)-JVS(87)*X(104))/(JVS(82))
  X(14) = (X(14)-JVS(76)*X(77)-JVS(77)*X(93)-JVS(78)*X(94)-JVS(79)*X(99)-JVS(80)*X(102)-JVS(81)*X(104))/(JVS(75))
  X(13) = (X(13)-JVS(67)*X(54)-JVS(68)*X(66)-JVS(69)*X(73)-JVS(70)*X(88)-JVS(71)*X(92)-JVS(72)*X(99)-JVS(73)*X(101)&
            &-JVS(74)*X(103))/(JVS(66))
  X(12) = (X(12)-JVS(62)*X(54)-JVS(63)*X(73)-JVS(64)*X(99)-JVS(65)*X(101))/(JVS(61))
  X(11) = (X(11)-JVS(57)*X(96)-JVS(58)*X(97)-JVS(59)*X(98)-JVS(60)*X(100))/(JVS(56))
  X(10) = (X(10)-JVS(54)*X(95)-JVS(55)*X(96))/(JVS(53))
  X(9) = (X(9)-JVS(39)*X(75)-JVS(40)*X(76)-JVS(41)*X(79)-JVS(42)*X(80)-JVS(43)*X(82)-JVS(44)*X(84)-JVS(45)*X(92)-JVS(46)&
           &*X(94)-JVS(47)*X(96)-JVS(48)*X(97)-JVS(49)*X(98)-JVS(50)*X(100)-JVS(51)*X(102)-JVS(52)*X(104))/(JVS(38))
  X(8) = (X(8)-JVS(29)*X(73)-JVS(30)*X(76)-JVS(31)*X(84)-JVS(32)*X(92)-JVS(33)*X(94)-JVS(34)*X(95)-JVS(35)*X(96)-JVS(36)&
           &*X(102)-JVS(37)*X(104))/(JVS(28))
  X(7) = (X(7)-JVS(13)*X(53)-JVS(14)*X(63)-JVS(15)*X(71)-JVS(16)*X(73)-JVS(17)*X(75)-JVS(18)*X(76)-JVS(19)*X(78)-JVS(20)&
           &*X(79)-JVS(21)*X(80)-JVS(22)*X(82)-JVS(23)*X(84)-JVS(24)*X(85)-JVS(25)*X(92)-JVS(26)*X(93)-JVS(27)*X(103))&
           &/(JVS(12))
  X(6) = X(6)/JVS(11)
  X(5) = X(5)/JVS(10)
  X(4) = X(4)/JVS(9)
  X(3) = X(3)/JVS(8)
  X(2) = (X(2)-JVS(5)*X(63)-JVS(6)*X(73)-JVS(7)*X(92))/(JVS(4))
  X(1) = (X(1)-JVS(2)*X(36)-JVS(3)*X(103))/(JVS(1))
END SUBROUTINE saprc99_mosaic_8bin_vbs2_KppSolve
      SUBROUTINE saprc99_mosaic_8bin_vbs2_WCOPY(N,X,incX,Y,incY)
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
      END SUBROUTINE saprc99_mosaic_8bin_vbs2_WCOPY
      SUBROUTINE saprc99_mosaic_8bin_vbs2_WAXPY(N,Alpha,X,incX,Y,incY)
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
      END SUBROUTINE saprc99_mosaic_8bin_vbs2_WAXPY
      SUBROUTINE saprc99_mosaic_8bin_vbs2_WSCAL(N,Alpha,X,incX)
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
      END SUBROUTINE saprc99_mosaic_8bin_vbs2_WSCAL
      REAL(kind=dp) FUNCTION saprc99_mosaic_8bin_vbs2_WLAMCH( C )
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
          CALL saprc99_mosaic_8bin_vbs2_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10 Eps = Eps*2
        i = i-1
      END IF
      saprc99_mosaic_8bin_vbs2_WLAMCH = Eps
      END FUNCTION saprc99_mosaic_8bin_vbs2_WLAMCH
      SUBROUTINE saprc99_mosaic_8bin_vbs2_WLAMCH_ADD( A, B, Sum )
      REAL(kind=dp) A, B, Sum
      Sum = A + B
      END SUBROUTINE saprc99_mosaic_8bin_vbs2_WLAMCH_ADD
      SUBROUTINE saprc99_mosaic_8bin_vbs2_SET2ZERO(N,Y)
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
      END SUBROUTINE saprc99_mosaic_8bin_vbs2_SET2ZERO
      REAL(kind=dp) FUNCTION saprc99_mosaic_8bin_vbs2_WDOT (N, DX, incX, DY, incY)
      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N)
      INTEGER :: i, IX, IY, M, MP1, NS
      saprc99_mosaic_8bin_vbs2_WDOT = 0.0D0
      IF (N .LE. 0) RETURN
      IF (incX .EQ. incY) IF (incX-1) 5,20,60
    5 IX = 1
      IY = 1
      IF (incX .LT. 0) IX = (-N+1)*incX + 1
      IF (incY .LT. 0) IY = (-N+1)*incY + 1
      DO i = 1,N
        saprc99_mosaic_8bin_vbs2_WDOT = saprc99_mosaic_8bin_vbs2_WDOT + DX(IX)*DY(IY)
        IX = IX + incX
        IY = IY + incY
      END DO
      RETURN
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO i = 1,M
         saprc99_mosaic_8bin_vbs2_WDOT = saprc99_mosaic_8bin_vbs2_WDOT + DX(i)*DY(i)
      END DO
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO i = MP1,N,5
          saprc99_mosaic_8bin_vbs2_WDOT = saprc99_mosaic_8bin_vbs2_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) + &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)
      END DO
      RETURN
   60 NS = N*incX
      DO i = 1,NS,incX
        saprc99_mosaic_8bin_vbs2_WDOT = saprc99_mosaic_8bin_vbs2_WDOT + DX(i)*DY(i)
      END DO
      END FUNCTION saprc99_mosaic_8bin_vbs2_WDOT
   SUBROUTINE decomp_saprc99_mosaic_8bin_vbs2( JVS, IER )
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
   W( 36 ) = JVS( 2 )
   W( 103 ) = JVS( 3 )
  JVS( 1) = W( 1 )
  JVS( 2) = W( 36 )
  JVS( 3) = W( 103 )
  IF ( ABS( JVS( 4 )) < TINY(a) ) THEN
         IER = 2
         RETURN
  END IF
   W( 2 ) = JVS( 4 )
   W( 63 ) = JVS( 5 )
   W( 73 ) = JVS( 6 )
   W( 92 ) = JVS( 7 )
  JVS( 4) = W( 2 )
  JVS( 5) = W( 63 )
  JVS( 6) = W( 73 )
  JVS( 7) = W( 92 )
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
   W( 53 ) = JVS( 13 )
   W( 63 ) = JVS( 14 )
   W( 71 ) = JVS( 15 )
   W( 73 ) = JVS( 16 )
   W( 75 ) = JVS( 17 )
   W( 76 ) = JVS( 18 )
   W( 78 ) = JVS( 19 )
   W( 79 ) = JVS( 20 )
   W( 80 ) = JVS( 21 )
   W( 82 ) = JVS( 22 )
   W( 84 ) = JVS( 23 )
   W( 85 ) = JVS( 24 )
   W( 92 ) = JVS( 25 )
   W( 93 ) = JVS( 26 )
   W( 103 ) = JVS( 27 )
  JVS( 12) = W( 7 )
  JVS( 13) = W( 53 )
  JVS( 14) = W( 63 )
  JVS( 15) = W( 71 )
  JVS( 16) = W( 73 )
  JVS( 17) = W( 75 )
  JVS( 18) = W( 76 )
  JVS( 19) = W( 78 )
  JVS( 20) = W( 79 )
  JVS( 21) = W( 80 )
  JVS( 22) = W( 82 )
  JVS( 23) = W( 84 )
  JVS( 24) = W( 85 )
  JVS( 25) = W( 92 )
  JVS( 26) = W( 93 )
  JVS( 27) = W( 103 )
  IF ( ABS( JVS( 28 )) < TINY(a) ) THEN
         IER = 8
         RETURN
  END IF
   W( 8 ) = JVS( 28 )
   W( 73 ) = JVS( 29 )
   W( 76 ) = JVS( 30 )
   W( 84 ) = JVS( 31 )
   W( 92 ) = JVS( 32 )
   W( 94 ) = JVS( 33 )
   W( 95 ) = JVS( 34 )
   W( 96 ) = JVS( 35 )
   W( 102 ) = JVS( 36 )
   W( 104 ) = JVS( 37 )
  JVS( 28) = W( 8 )
  JVS( 29) = W( 73 )
  JVS( 30) = W( 76 )
  JVS( 31) = W( 84 )
  JVS( 32) = W( 92 )
  JVS( 33) = W( 94 )
  JVS( 34) = W( 95 )
  JVS( 35) = W( 96 )
  JVS( 36) = W( 102 )
  JVS( 37) = W( 104 )
  IF ( ABS( JVS( 38 )) < TINY(a) ) THEN
         IER = 9
         RETURN
  END IF
   W( 9 ) = JVS( 38 )
   W( 75 ) = JVS( 39 )
   W( 76 ) = JVS( 40 )
   W( 79 ) = JVS( 41 )
   W( 80 ) = JVS( 42 )
   W( 82 ) = JVS( 43 )
   W( 84 ) = JVS( 44 )
   W( 92 ) = JVS( 45 )
   W( 94 ) = JVS( 46 )
   W( 96 ) = JVS( 47 )
   W( 97 ) = JVS( 48 )
   W( 98 ) = JVS( 49 )
   W( 100 ) = JVS( 50 )
   W( 102 ) = JVS( 51 )
   W( 104 ) = JVS( 52 )
  JVS( 38) = W( 9 )
  JVS( 39) = W( 75 )
  JVS( 40) = W( 76 )
  JVS( 41) = W( 79 )
  JVS( 42) = W( 80 )
  JVS( 43) = W( 82 )
  JVS( 44) = W( 84 )
  JVS( 45) = W( 92 )
  JVS( 46) = W( 94 )
  JVS( 47) = W( 96 )
  JVS( 48) = W( 97 )
  JVS( 49) = W( 98 )
  JVS( 50) = W( 100 )
  JVS( 51) = W( 102 )
  JVS( 52) = W( 104 )
  IF ( ABS( JVS( 53 )) < TINY(a) ) THEN
         IER = 10
         RETURN
  END IF
   W( 10 ) = JVS( 53 )
   W( 95 ) = JVS( 54 )
   W( 96 ) = JVS( 55 )
  JVS( 53) = W( 10 )
  JVS( 54) = W( 95 )
  JVS( 55) = W( 96 )
  IF ( ABS( JVS( 56 )) < TINY(a) ) THEN
         IER = 11
         RETURN
  END IF
   W( 11 ) = JVS( 56 )
   W( 96 ) = JVS( 57 )
   W( 97 ) = JVS( 58 )
   W( 98 ) = JVS( 59 )
   W( 100 ) = JVS( 60 )
  JVS( 56) = W( 11 )
  JVS( 57) = W( 96 )
  JVS( 58) = W( 97 )
  JVS( 59) = W( 98 )
  JVS( 60) = W( 100 )
  IF ( ABS( JVS( 61 )) < TINY(a) ) THEN
         IER = 12
         RETURN
  END IF
   W( 12 ) = JVS( 61 )
   W( 54 ) = JVS( 62 )
   W( 73 ) = JVS( 63 )
   W( 99 ) = JVS( 64 )
   W( 101 ) = JVS( 65 )
  JVS( 61) = W( 12 )
  JVS( 62) = W( 54 )
  JVS( 63) = W( 73 )
  JVS( 64) = W( 99 )
  JVS( 65) = W( 101 )
  IF ( ABS( JVS( 66 )) < TINY(a) ) THEN
         IER = 13
         RETURN
  END IF
   W( 13 ) = JVS( 66 )
   W( 54 ) = JVS( 67 )
   W( 66 ) = JVS( 68 )
   W( 73 ) = JVS( 69 )
   W( 88 ) = JVS( 70 )
   W( 92 ) = JVS( 71 )
   W( 99 ) = JVS( 72 )
   W( 101 ) = JVS( 73 )
   W( 103 ) = JVS( 74 )
  JVS( 66) = W( 13 )
  JVS( 67) = W( 54 )
  JVS( 68) = W( 66 )
  JVS( 69) = W( 73 )
  JVS( 70) = W( 88 )
  JVS( 71) = W( 92 )
  JVS( 72) = W( 99 )
  JVS( 73) = W( 101 )
  JVS( 74) = W( 103 )
  IF ( ABS( JVS( 75 )) < TINY(a) ) THEN
         IER = 14
         RETURN
  END IF
   W( 14 ) = JVS( 75 )
   W( 77 ) = JVS( 76 )
   W( 93 ) = JVS( 77 )
   W( 94 ) = JVS( 78 )
   W( 99 ) = JVS( 79 )
   W( 102 ) = JVS( 80 )
   W( 104 ) = JVS( 81 )
  JVS( 75) = W( 14 )
  JVS( 76) = W( 77 )
  JVS( 77) = W( 93 )
  JVS( 78) = W( 94 )
  JVS( 79) = W( 99 )
  JVS( 80) = W( 102 )
  JVS( 81) = W( 104 )
  IF ( ABS( JVS( 82 )) < TINY(a) ) THEN
         IER = 15
         RETURN
  END IF
   W( 15 ) = JVS( 82 )
   W( 77 ) = JVS( 83 )
   W( 94 ) = JVS( 84 )
   W( 96 ) = JVS( 85 )
   W( 102 ) = JVS( 86 )
   W( 104 ) = JVS( 87 )
  JVS( 82) = W( 15 )
  JVS( 83) = W( 77 )
  JVS( 84) = W( 94 )
  JVS( 85) = W( 96 )
  JVS( 86) = W( 102 )
  JVS( 87) = W( 104 )
  IF ( ABS( JVS( 88 )) < TINY(a) ) THEN
         IER = 16
         RETURN
  END IF
   W( 16 ) = JVS( 88 )
   W( 51 ) = JVS( 89 )
   W( 57 ) = JVS( 90 )
   W( 103 ) = JVS( 91 )
  JVS( 88) = W( 16 )
  JVS( 89) = W( 51 )
  JVS( 90) = W( 57 )
  JVS( 91) = W( 103 )
  IF ( ABS( JVS( 92 )) < TINY(a) ) THEN
         IER = 17
         RETURN
  END IF
   W( 17 ) = JVS( 92 )
   W( 51 ) = JVS( 93 )
   W( 57 ) = JVS( 94 )
   W( 103 ) = JVS( 95 )
  JVS( 92) = W( 17 )
  JVS( 93) = W( 51 )
  JVS( 94) = W( 57 )
  JVS( 95) = W( 103 )
  IF ( ABS( JVS( 96 )) < TINY(a) ) THEN
         IER = 18
         RETURN
  END IF
   W( 18 ) = JVS( 96 )
   W( 51 ) = JVS( 97 )
   W( 57 ) = JVS( 98 )
   W( 103 ) = JVS( 99 )
  JVS( 96) = W( 18 )
  JVS( 97) = W( 51 )
  JVS( 98) = W( 57 )
  JVS( 99) = W( 103 )
  IF ( ABS( JVS( 100 )) < TINY(a) ) THEN
         IER = 19
         RETURN
  END IF
   W( 19 ) = JVS( 100 )
   W( 51 ) = JVS( 101 )
   W( 57 ) = JVS( 102 )
   W( 103 ) = JVS( 103 )
  JVS( 100) = W( 19 )
  JVS( 101) = W( 51 )
  JVS( 102) = W( 57 )
  JVS( 103) = W( 103 )
  IF ( ABS( JVS( 104 )) < TINY(a) ) THEN
         IER = 20
         RETURN
  END IF
   W( 20 ) = JVS( 104 )
   W( 75 ) = JVS( 105 )
   W( 103 ) = JVS( 106 )
  JVS( 104) = W( 20 )
  JVS( 105) = W( 75 )
  JVS( 106) = W( 103 )
  IF ( ABS( JVS( 107 )) < TINY(a) ) THEN
         IER = 21
         RETURN
  END IF
   W( 21 ) = JVS( 107 )
   W( 75 ) = JVS( 108 )
   W( 79 ) = JVS( 109 )
   W( 80 ) = JVS( 110 )
   W( 92 ) = JVS( 111 )
   W( 99 ) = JVS( 112 )
   W( 103 ) = JVS( 113 )
  JVS( 107) = W( 21 )
  JVS( 108) = W( 75 )
  JVS( 109) = W( 79 )
  JVS( 110) = W( 80 )
  JVS( 111) = W( 92 )
  JVS( 112) = W( 99 )
  JVS( 113) = W( 103 )
  IF ( ABS( JVS( 114 )) < TINY(a) ) THEN
         IER = 22
         RETURN
  END IF
   W( 22 ) = JVS( 114 )
   W( 75 ) = JVS( 115 )
   W( 79 ) = JVS( 116 )
   W( 80 ) = JVS( 117 )
   W( 92 ) = JVS( 118 )
   W( 99 ) = JVS( 119 )
   W( 103 ) = JVS( 120 )
  JVS( 114) = W( 22 )
  JVS( 115) = W( 75 )
  JVS( 116) = W( 79 )
  JVS( 117) = W( 80 )
  JVS( 118) = W( 92 )
  JVS( 119) = W( 99 )
  JVS( 120) = W( 103 )
  IF ( ABS( JVS( 121 )) < TINY(a) ) THEN
         IER = 23
         RETURN
  END IF
   W( 23 ) = JVS( 121 )
   W( 75 ) = JVS( 122 )
   W( 79 ) = JVS( 123 )
   W( 80 ) = JVS( 124 )
   W( 92 ) = JVS( 125 )
   W( 99 ) = JVS( 126 )
   W( 103 ) = JVS( 127 )
  JVS( 121) = W( 23 )
  JVS( 122) = W( 75 )
  JVS( 123) = W( 79 )
  JVS( 124) = W( 80 )
  JVS( 125) = W( 92 )
  JVS( 126) = W( 99 )
  JVS( 127) = W( 103 )
  IF ( ABS( JVS( 128 )) < TINY(a) ) THEN
         IER = 24
         RETURN
  END IF
   W( 24 ) = JVS( 128 )
   W( 35 ) = JVS( 129 )
   W( 36 ) = JVS( 130 )
   W( 37 ) = JVS( 131 )
   W( 43 ) = JVS( 132 )
   W( 44 ) = JVS( 133 )
   W( 47 ) = JVS( 134 )
   W( 48 ) = JVS( 135 )
   W( 50 ) = JVS( 136 )
   W( 51 ) = JVS( 137 )
   W( 52 ) = JVS( 138 )
   W( 56 ) = JVS( 139 )
   W( 57 ) = JVS( 140 )
   W( 58 ) = JVS( 141 )
   W( 59 ) = JVS( 142 )
   W( 60 ) = JVS( 143 )
   W( 61 ) = JVS( 144 )
   W( 62 ) = JVS( 145 )
   W( 63 ) = JVS( 146 )
   W( 65 ) = JVS( 147 )
   W( 66 ) = JVS( 148 )
   W( 67 ) = JVS( 149 )
   W( 68 ) = JVS( 150 )
   W( 69 ) = JVS( 151 )
   W( 70 ) = JVS( 152 )
   W( 71 ) = JVS( 153 )
   W( 72 ) = JVS( 154 )
   W( 73 ) = JVS( 155 )
   W( 75 ) = JVS( 156 )
   W( 76 ) = JVS( 157 )
   W( 78 ) = JVS( 158 )
   W( 79 ) = JVS( 159 )
   W( 80 ) = JVS( 160 )
   W( 81 ) = JVS( 161 )
   W( 82 ) = JVS( 162 )
   W( 83 ) = JVS( 163 )
   W( 84 ) = JVS( 164 )
   W( 85 ) = JVS( 165 )
   W( 86 ) = JVS( 166 )
   W( 87 ) = JVS( 167 )
   W( 89 ) = JVS( 168 )
   W( 90 ) = JVS( 169 )
   W( 91 ) = JVS( 170 )
   W( 92 ) = JVS( 171 )
   W( 93 ) = JVS( 172 )
   W( 96 ) = JVS( 173 )
   W( 99 ) = JVS( 174 )
   W( 101 ) = JVS( 175 )
   W( 103 ) = JVS( 176 )
  JVS( 128) = W( 24 )
  JVS( 129) = W( 35 )
  JVS( 130) = W( 36 )
  JVS( 131) = W( 37 )
  JVS( 132) = W( 43 )
  JVS( 133) = W( 44 )
  JVS( 134) = W( 47 )
  JVS( 135) = W( 48 )
  JVS( 136) = W( 50 )
  JVS( 137) = W( 51 )
  JVS( 138) = W( 52 )
  JVS( 139) = W( 56 )
  JVS( 140) = W( 57 )
  JVS( 141) = W( 58 )
  JVS( 142) = W( 59 )
  JVS( 143) = W( 60 )
  JVS( 144) = W( 61 )
  JVS( 145) = W( 62 )
  JVS( 146) = W( 63 )
  JVS( 147) = W( 65 )
  JVS( 148) = W( 66 )
  JVS( 149) = W( 67 )
  JVS( 150) = W( 68 )
  JVS( 151) = W( 69 )
  JVS( 152) = W( 70 )
  JVS( 153) = W( 71 )
  JVS( 154) = W( 72 )
  JVS( 155) = W( 73 )
  JVS( 156) = W( 75 )
  JVS( 157) = W( 76 )
  JVS( 158) = W( 78 )
  JVS( 159) = W( 79 )
  JVS( 160) = W( 80 )
  JVS( 161) = W( 81 )
  JVS( 162) = W( 82 )
  JVS( 163) = W( 83 )
  JVS( 164) = W( 84 )
  JVS( 165) = W( 85 )
  JVS( 166) = W( 86 )
  JVS( 167) = W( 87 )
  JVS( 168) = W( 89 )
  JVS( 169) = W( 90 )
  JVS( 170) = W( 91 )
  JVS( 171) = W( 92 )
  JVS( 172) = W( 93 )
  JVS( 173) = W( 96 )
  JVS( 174) = W( 99 )
  JVS( 175) = W( 101 )
  JVS( 176) = W( 103 )
  IF ( ABS( JVS( 181 )) < TINY(a) ) THEN
         IER = 25
         RETURN
  END IF
   W( 3 ) = JVS( 177 )
   W( 4 ) = JVS( 178 )
   W( 5 ) = JVS( 179 )
   W( 6 ) = JVS( 180 )
   W( 25 ) = JVS( 181 )
   W( 26 ) = JVS( 182 )
   W( 27 ) = JVS( 183 )
   W( 28 ) = JVS( 184 )
   W( 29 ) = JVS( 185 )
   W( 31 ) = JVS( 186 )
   W( 32 ) = JVS( 187 )
   W( 33 ) = JVS( 188 )
   W( 34 ) = JVS( 189 )
   W( 103 ) = JVS( 190 )
  a = -W( 3 ) / JVS( 8 )
  W( 3 ) = -a
  a = -W( 4 ) / JVS( 9 )
  W( 4 ) = -a
  a = -W( 5 ) / JVS( 10 )
  W( 5 ) = -a
  a = -W( 6 ) / JVS( 11 )
  W( 6 ) = -a
  JVS( 177) = W( 3 )
  JVS( 178) = W( 4 )
  JVS( 179) = W( 5 )
  JVS( 180) = W( 6 )
  JVS( 181) = W( 25 )
  JVS( 182) = W( 26 )
  JVS( 183) = W( 27 )
  JVS( 184) = W( 28 )
  JVS( 185) = W( 29 )
  JVS( 186) = W( 31 )
  JVS( 187) = W( 32 )
  JVS( 188) = W( 33 )
  JVS( 189) = W( 34 )
  JVS( 190) = W( 103 )
  IF ( ABS( JVS( 191 )) < TINY(a) ) THEN
         IER = 26
         RETURN
  END IF
   W( 26 ) = JVS( 191 )
   W( 27 ) = JVS( 192 )
   W( 31 ) = JVS( 193 )
   W( 103 ) = JVS( 194 )
  JVS( 191) = W( 26 )
  JVS( 192) = W( 27 )
  JVS( 193) = W( 31 )
  JVS( 194) = W( 103 )
  IF ( ABS( JVS( 195 )) < TINY(a) ) THEN
         IER = 27
         RETURN
  END IF
   W( 27 ) = JVS( 195 )
   W( 103 ) = JVS( 196 )
  JVS( 195) = W( 27 )
  JVS( 196) = W( 103 )
  IF ( ABS( JVS( 197 )) < TINY(a) ) THEN
         IER = 28
         RETURN
  END IF
   W( 28 ) = JVS( 197 )
   W( 29 ) = JVS( 198 )
   W( 33 ) = JVS( 199 )
   W( 103 ) = JVS( 200 )
  JVS( 197) = W( 28 )
  JVS( 198) = W( 29 )
  JVS( 199) = W( 33 )
  JVS( 200) = W( 103 )
  IF ( ABS( JVS( 201 )) < TINY(a) ) THEN
         IER = 29
         RETURN
  END IF
   W( 29 ) = JVS( 201 )
   W( 103 ) = JVS( 202 )
  JVS( 201) = W( 29 )
  JVS( 202) = W( 103 )
  IF ( ABS( JVS( 203 )) < TINY(a) ) THEN
         IER = 30
         RETURN
  END IF
   W( 30 ) = JVS( 203 )
   W( 92 ) = JVS( 204 )
  JVS( 203) = W( 30 )
  JVS( 204) = W( 92 )
  IF ( ABS( JVS( 205 )) < TINY(a) ) THEN
         IER = 31
         RETURN
  END IF
   W( 31 ) = JVS( 205 )
   W( 103 ) = JVS( 206 )
  JVS( 205) = W( 31 )
  JVS( 206) = W( 103 )
  IF ( ABS( JVS( 208 )) < TINY(a) ) THEN
         IER = 32
         RETURN
  END IF
   W( 31 ) = JVS( 207 )
   W( 32 ) = JVS( 208 )
   W( 103 ) = JVS( 209 )
  a = -W( 31 ) / JVS( 205 )
  W( 31 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 206 )
  JVS( 207) = W( 31 )
  JVS( 208) = W( 32 )
  JVS( 209) = W( 103 )
  IF ( ABS( JVS( 210 )) < TINY(a) ) THEN
         IER = 33
         RETURN
  END IF
   W( 33 ) = JVS( 210 )
   W( 103 ) = JVS( 211 )
  JVS( 210) = W( 33 )
  JVS( 211) = W( 103 )
  IF ( ABS( JVS( 213 )) < TINY(a) ) THEN
         IER = 34
         RETURN
  END IF
   W( 33 ) = JVS( 212 )
   W( 34 ) = JVS( 213 )
   W( 103 ) = JVS( 214 )
  a = -W( 33 ) / JVS( 210 )
  W( 33 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 211 )
  JVS( 212) = W( 33 )
  JVS( 213) = W( 34 )
  JVS( 214) = W( 103 )
  IF ( ABS( JVS( 215 )) < TINY(a) ) THEN
         IER = 35
         RETURN
  END IF
   W( 35 ) = JVS( 215 )
   W( 103 ) = JVS( 216 )
  JVS( 215) = W( 35 )
  JVS( 216) = W( 103 )
  IF ( ABS( JVS( 217 )) < TINY(a) ) THEN
         IER = 36
         RETURN
  END IF
   W( 36 ) = JVS( 217 )
   W( 103 ) = JVS( 218 )
  JVS( 217) = W( 36 )
  JVS( 218) = W( 103 )
  IF ( ABS( JVS( 219 )) < TINY(a) ) THEN
         IER = 37
         RETURN
  END IF
   W( 37 ) = JVS( 219 )
   W( 103 ) = JVS( 220 )
  JVS( 219) = W( 37 )
  JVS( 220) = W( 103 )
  IF ( ABS( JVS( 221 )) < TINY(a) ) THEN
         IER = 38
         RETURN
  END IF
   W( 38 ) = JVS( 221 )
   W( 95 ) = JVS( 222 )
   W( 101 ) = JVS( 223 )
  JVS( 221) = W( 38 )
  JVS( 222) = W( 95 )
  JVS( 223) = W( 101 )
  IF ( ABS( JVS( 224 )) < TINY(a) ) THEN
         IER = 39
         RETURN
  END IF
   W( 39 ) = JVS( 224 )
   W( 97 ) = JVS( 225 )
   W( 101 ) = JVS( 226 )
  JVS( 224) = W( 39 )
  JVS( 225) = W( 97 )
  JVS( 226) = W( 101 )
  IF ( ABS( JVS( 227 )) < TINY(a) ) THEN
         IER = 40
         RETURN
  END IF
   W( 40 ) = JVS( 227 )
   W( 100 ) = JVS( 228 )
   W( 101 ) = JVS( 229 )
  JVS( 227) = W( 40 )
  JVS( 228) = W( 100 )
  JVS( 229) = W( 101 )
  IF ( ABS( JVS( 230 )) < TINY(a) ) THEN
         IER = 41
         RETURN
  END IF
   W( 41 ) = JVS( 230 )
   W( 98 ) = JVS( 231 )
   W( 101 ) = JVS( 232 )
  JVS( 230) = W( 41 )
  JVS( 231) = W( 98 )
  JVS( 232) = W( 101 )
  IF ( ABS( JVS( 233 )) < TINY(a) ) THEN
         IER = 42
         RETURN
  END IF
   W( 42 ) = JVS( 233 )
   W( 96 ) = JVS( 234 )
   W( 103 ) = JVS( 235 )
  JVS( 233) = W( 42 )
  JVS( 234) = W( 96 )
  JVS( 235) = W( 103 )
  IF ( ABS( JVS( 236 )) < TINY(a) ) THEN
         IER = 43
         RETURN
  END IF
   W( 43 ) = JVS( 236 )
   W( 103 ) = JVS( 237 )
  JVS( 236) = W( 43 )
  JVS( 237) = W( 103 )
  IF ( ABS( JVS( 238 )) < TINY(a) ) THEN
         IER = 44
         RETURN
  END IF
   W( 44 ) = JVS( 238 )
   W( 103 ) = JVS( 239 )
  JVS( 238) = W( 44 )
  JVS( 239) = W( 103 )
  IF ( ABS( JVS( 240 )) < TINY(a) ) THEN
         IER = 45
         RETURN
  END IF
   W( 45 ) = JVS( 240 )
   W( 51 ) = JVS( 241 )
   W( 79 ) = JVS( 242 )
   W( 80 ) = JVS( 243 )
   W( 92 ) = JVS( 244 )
   W( 103 ) = JVS( 245 )
  JVS( 240) = W( 45 )
  JVS( 241) = W( 51 )
  JVS( 242) = W( 79 )
  JVS( 243) = W( 80 )
  JVS( 244) = W( 92 )
  JVS( 245) = W( 103 )
  IF ( ABS( JVS( 246 )) < TINY(a) ) THEN
         IER = 46
         RETURN
  END IF
   W( 46 ) = JVS( 246 )
   W( 99 ) = JVS( 247 )
   W( 101 ) = JVS( 248 )
  JVS( 246) = W( 46 )
  JVS( 247) = W( 99 )
  JVS( 248) = W( 101 )
  IF ( ABS( JVS( 249 )) < TINY(a) ) THEN
         IER = 47
         RETURN
  END IF
   W( 47 ) = JVS( 249 )
   W( 93 ) = JVS( 250 )
   W( 103 ) = JVS( 251 )
  JVS( 249) = W( 47 )
  JVS( 250) = W( 93 )
  JVS( 251) = W( 103 )
  IF ( ABS( JVS( 252 )) < TINY(a) ) THEN
         IER = 48
         RETURN
  END IF
   W( 48 ) = JVS( 252 )
   W( 103 ) = JVS( 253 )
  JVS( 252) = W( 48 )
  JVS( 253) = W( 103 )
  IF ( ABS( JVS( 255 )) < TINY(a) ) THEN
         IER = 49
         RETURN
  END IF
   W( 48 ) = JVS( 254 )
   W( 49 ) = JVS( 255 )
   W( 101 ) = JVS( 256 )
   W( 103 ) = JVS( 257 )
  a = -W( 48 ) / JVS( 252 )
  W( 48 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 253 )
  JVS( 254) = W( 48 )
  JVS( 255) = W( 49 )
  JVS( 256) = W( 101 )
  JVS( 257) = W( 103 )
  IF ( ABS( JVS( 258 )) < TINY(a) ) THEN
         IER = 50
         RETURN
  END IF
   W( 50 ) = JVS( 258 )
   W( 103 ) = JVS( 259 )
  JVS( 258) = W( 50 )
  JVS( 259) = W( 103 )
  IF ( ABS( JVS( 260 )) < TINY(a) ) THEN
         IER = 51
         RETURN
  END IF
   W( 51 ) = JVS( 260 )
   W( 103 ) = JVS( 261 )
  JVS( 260) = W( 51 )
  JVS( 261) = W( 103 )
  IF ( ABS( JVS( 262 )) < TINY(a) ) THEN
         IER = 52
         RETURN
  END IF
   W( 52 ) = JVS( 262 )
   W( 96 ) = JVS( 263 )
   W( 102 ) = JVS( 264 )
   W( 103 ) = JVS( 265 )
  JVS( 262) = W( 52 )
  JVS( 263) = W( 96 )
  JVS( 264) = W( 102 )
  JVS( 265) = W( 103 )
  IF ( ABS( JVS( 266 )) < TINY(a) ) THEN
         IER = 53
         RETURN
  END IF
   W( 53 ) = JVS( 266 )
   W( 86 ) = JVS( 267 )
   W( 93 ) = JVS( 268 )
   W( 96 ) = JVS( 269 )
  JVS( 266) = W( 53 )
  JVS( 267) = W( 86 )
  JVS( 268) = W( 93 )
  JVS( 269) = W( 96 )
  IF ( ABS( JVS( 270 )) < TINY(a) ) THEN
         IER = 54
         RETURN
  END IF
   W( 54 ) = JVS( 270 )
   W( 64 ) = JVS( 271 )
   W( 96 ) = JVS( 272 )
   W( 99 ) = JVS( 273 )
   W( 101 ) = JVS( 274 )
  JVS( 270) = W( 54 )
  JVS( 271) = W( 64 )
  JVS( 272) = W( 96 )
  JVS( 273) = W( 99 )
  JVS( 274) = W( 101 )
  IF ( ABS( JVS( 275 )) < TINY(a) ) THEN
         IER = 55
         RETURN
  END IF
   W( 55 ) = JVS( 275 )
   W( 96 ) = JVS( 276 )
   W( 101 ) = JVS( 277 )
   W( 103 ) = JVS( 278 )
  JVS( 275) = W( 55 )
  JVS( 276) = W( 96 )
  JVS( 277) = W( 101 )
  JVS( 278) = W( 103 )
  IF ( ABS( JVS( 279 )) < TINY(a) ) THEN
         IER = 56
         RETURN
  END IF
   W( 56 ) = JVS( 279 )
   W( 94 ) = JVS( 280 )
   W( 102 ) = JVS( 281 )
   W( 103 ) = JVS( 282 )
   W( 104 ) = JVS( 283 )
  JVS( 279) = W( 56 )
  JVS( 280) = W( 94 )
  JVS( 281) = W( 102 )
  JVS( 282) = W( 103 )
  JVS( 283) = W( 104 )
  IF ( ABS( JVS( 284 )) < TINY(a) ) THEN
         IER = 57
         RETURN
  END IF
   W( 57 ) = JVS( 284 )
   W( 103 ) = JVS( 285 )
  JVS( 284) = W( 57 )
  JVS( 285) = W( 103 )
  IF ( ABS( JVS( 288 )) < TINY(a) ) THEN
         IER = 58
         RETURN
  END IF
   W( 51 ) = JVS( 286 )
   W( 57 ) = JVS( 287 )
   W( 58 ) = JVS( 288 )
   W( 103 ) = JVS( 289 )
  a = -W( 51 ) / JVS( 260 )
  W( 51 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 261 )
  a = -W( 57 ) / JVS( 284 )
  W( 57 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 285 )
  JVS( 286) = W( 51 )
  JVS( 287) = W( 57 )
  JVS( 288) = W( 58 )
  JVS( 289) = W( 103 )
  IF ( ABS( JVS( 292 )) < TINY(a) ) THEN
         IER = 59
         RETURN
  END IF
   W( 51 ) = JVS( 290 )
   W( 57 ) = JVS( 291 )
   W( 59 ) = JVS( 292 )
   W( 103 ) = JVS( 293 )
  a = -W( 51 ) / JVS( 260 )
  W( 51 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 261 )
  a = -W( 57 ) / JVS( 284 )
  W( 57 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 285 )
  JVS( 290) = W( 51 )
  JVS( 291) = W( 57 )
  JVS( 292) = W( 59 )
  JVS( 293) = W( 103 )
  IF ( ABS( JVS( 296 )) < TINY(a) ) THEN
         IER = 60
         RETURN
  END IF
   W( 51 ) = JVS( 294 )
   W( 57 ) = JVS( 295 )
   W( 60 ) = JVS( 296 )
   W( 99 ) = JVS( 297 )
   W( 103 ) = JVS( 298 )
  a = -W( 51 ) / JVS( 260 )
  W( 51 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 261 )
  a = -W( 57 ) / JVS( 284 )
  W( 57 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 285 )
  JVS( 294) = W( 51 )
  JVS( 295) = W( 57 )
  JVS( 296) = W( 60 )
  JVS( 297) = W( 99 )
  JVS( 298) = W( 103 )
  IF ( ABS( JVS( 299 )) < TINY(a) ) THEN
         IER = 61
         RETURN
  END IF
   W( 61 ) = JVS( 299 )
   W( 103 ) = JVS( 300 )
  JVS( 299) = W( 61 )
  JVS( 300) = W( 103 )
  IF ( ABS( JVS( 303 )) < TINY(a) ) THEN
         IER = 62
         RETURN
  END IF
   W( 51 ) = JVS( 301 )
   W( 57 ) = JVS( 302 )
   W( 62 ) = JVS( 303 )
   W( 92 ) = JVS( 304 )
   W( 103 ) = JVS( 305 )
  a = -W( 51 ) / JVS( 260 )
  W( 51 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 261 )
  a = -W( 57 ) / JVS( 284 )
  W( 57 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 285 )
  JVS( 301) = W( 51 )
  JVS( 302) = W( 57 )
  JVS( 303) = W( 62 )
  JVS( 304) = W( 92 )
  JVS( 305) = W( 103 )
  IF ( ABS( JVS( 306 )) < TINY(a) ) THEN
         IER = 63
         RETURN
  END IF
   W( 63 ) = JVS( 306 )
   W( 92 ) = JVS( 307 )
   W( 103 ) = JVS( 308 )
  JVS( 306) = W( 63 )
  JVS( 307) = W( 92 )
  JVS( 308) = W( 103 )
  IF ( ABS( JVS( 310 )) < TINY(a) ) THEN
         IER = 64
         RETURN
  END IF
   W( 54 ) = JVS( 309 )
   W( 64 ) = JVS( 310 )
   W( 74 ) = JVS( 311 )
   W( 96 ) = JVS( 312 )
   W( 99 ) = JVS( 313 )
   W( 101 ) = JVS( 314 )
  a = -W( 54 ) / JVS( 270 )
  W( 54 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 271 )
  W( 96 ) = W( 96 ) + a*JVS( 272 )
  W( 99 ) = W( 99 ) + a*JVS( 273 )
  W( 101 ) = W( 101 ) + a*JVS( 274 )
  JVS( 309) = W( 54 )
  JVS( 310) = W( 64 )
  JVS( 311) = W( 74 )
  JVS( 312) = W( 96 )
  JVS( 313) = W( 99 )
  JVS( 314) = W( 101 )
  IF ( ABS( JVS( 315 )) < TINY(a) ) THEN
         IER = 65
         RETURN
  END IF
   W( 65 ) = JVS( 315 )
   W( 94 ) = JVS( 316 )
   W( 96 ) = JVS( 317 )
   W( 103 ) = JVS( 318 )
   W( 104 ) = JVS( 319 )
  JVS( 315) = W( 65 )
  JVS( 316) = W( 94 )
  JVS( 317) = W( 96 )
  JVS( 318) = W( 103 )
  JVS( 319) = W( 104 )
  IF ( ABS( JVS( 322 )) < TINY(a) ) THEN
         IER = 66
         RETURN
  END IF
   W( 51 ) = JVS( 320 )
   W( 57 ) = JVS( 321 )
   W( 66 ) = JVS( 322 )
   W( 84 ) = JVS( 323 )
   W( 92 ) = JVS( 324 )
   W( 99 ) = JVS( 325 )
   W( 103 ) = JVS( 326 )
  a = -W( 51 ) / JVS( 260 )
  W( 51 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 261 )
  a = -W( 57 ) / JVS( 284 )
  W( 57 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 285 )
  JVS( 320) = W( 51 )
  JVS( 321) = W( 57 )
  JVS( 322) = W( 66 )
  JVS( 323) = W( 84 )
  JVS( 324) = W( 92 )
  JVS( 325) = W( 99 )
  JVS( 326) = W( 103 )
  IF ( ABS( JVS( 328 )) < TINY(a) ) THEN
         IER = 67
         RETURN
  END IF
   W( 57 ) = JVS( 327 )
   W( 67 ) = JVS( 328 )
   W( 74 ) = JVS( 329 )
   W( 96 ) = JVS( 330 )
   W( 99 ) = JVS( 331 )
   W( 103 ) = JVS( 332 )
  a = -W( 57 ) / JVS( 284 )
  W( 57 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 285 )
  JVS( 327) = W( 57 )
  JVS( 328) = W( 67 )
  JVS( 329) = W( 74 )
  JVS( 330) = W( 96 )
  JVS( 331) = W( 99 )
  JVS( 332) = W( 103 )
  IF ( ABS( JVS( 338 )) < TINY(a) ) THEN
         IER = 68
         RETURN
  END IF
   W( 51 ) = JVS( 333 )
   W( 57 ) = JVS( 334 )
   W( 58 ) = JVS( 335 )
   W( 59 ) = JVS( 336 )
   W( 60 ) = JVS( 337 )
   W( 68 ) = JVS( 338 )
   W( 78 ) = JVS( 339 )
   W( 82 ) = JVS( 340 )
   W( 85 ) = JVS( 341 )
   W( 92 ) = JVS( 342 )
   W( 99 ) = JVS( 343 )
   W( 103 ) = JVS( 344 )
  a = -W( 51 ) / JVS( 260 )
  W( 51 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 261 )
  a = -W( 57 ) / JVS( 284 )
  W( 57 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 285 )
  a = -W( 58 ) / JVS( 288 )
  W( 58 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 289 )
  a = -W( 59 ) / JVS( 292 )
  W( 59 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 293 )
  a = -W( 60 ) / JVS( 296 )
  W( 60 ) = -a
  W( 99 ) = W( 99 ) + a*JVS( 297 )
  W( 103 ) = W( 103 ) + a*JVS( 298 )
  JVS( 333) = W( 51 )
  JVS( 334) = W( 57 )
  JVS( 335) = W( 58 )
  JVS( 336) = W( 59 )
  JVS( 337) = W( 60 )
  JVS( 338) = W( 68 )
  JVS( 339) = W( 78 )
  JVS( 340) = W( 82 )
  JVS( 341) = W( 85 )
  JVS( 342) = W( 92 )
  JVS( 343) = W( 99 )
  JVS( 344) = W( 103 )
  IF ( ABS( JVS( 351 )) < TINY(a) ) THEN
         IER = 69
         RETURN
  END IF
   W( 58 ) = JVS( 345 )
   W( 59 ) = JVS( 346 )
   W( 61 ) = JVS( 347 )
   W( 62 ) = JVS( 348 )
   W( 63 ) = JVS( 349 )
   W( 68 ) = JVS( 350 )
   W( 69 ) = JVS( 351 )
   W( 71 ) = JVS( 352 )
   W( 73 ) = JVS( 353 )
   W( 75 ) = JVS( 354 )
   W( 76 ) = JVS( 355 )
   W( 78 ) = JVS( 356 )
   W( 79 ) = JVS( 357 )
   W( 80 ) = JVS( 358 )
   W( 81 ) = JVS( 359 )
   W( 82 ) = JVS( 360 )
   W( 83 ) = JVS( 361 )
   W( 84 ) = JVS( 362 )
   W( 85 ) = JVS( 363 )
   W( 86 ) = JVS( 364 )
   W( 88 ) = JVS( 365 )
   W( 89 ) = JVS( 366 )
   W( 92 ) = JVS( 367 )
   W( 99 ) = JVS( 368 )
   W( 103 ) = JVS( 369 )
  a = -W( 58 ) / JVS( 288 )
  W( 58 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 289 )
  a = -W( 59 ) / JVS( 292 )
  W( 59 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 293 )
  a = -W( 61 ) / JVS( 299 )
  W( 61 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 300 )
  a = -W( 62 ) / JVS( 303 )
  W( 62 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 304 )
  W( 103 ) = W( 103 ) + a*JVS( 305 )
  a = -W( 63 ) / JVS( 306 )
  W( 63 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 307 )
  W( 103 ) = W( 103 ) + a*JVS( 308 )
  a = -W( 68 ) / JVS( 338 )
  W( 68 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 339 )
  W( 82 ) = W( 82 ) + a*JVS( 340 )
  W( 85 ) = W( 85 ) + a*JVS( 341 )
  W( 92 ) = W( 92 ) + a*JVS( 342 )
  W( 99 ) = W( 99 ) + a*JVS( 343 )
  W( 103 ) = W( 103 ) + a*JVS( 344 )
  JVS( 345) = W( 58 )
  JVS( 346) = W( 59 )
  JVS( 347) = W( 61 )
  JVS( 348) = W( 62 )
  JVS( 349) = W( 63 )
  JVS( 350) = W( 68 )
  JVS( 351) = W( 69 )
  JVS( 352) = W( 71 )
  JVS( 353) = W( 73 )
  JVS( 354) = W( 75 )
  JVS( 355) = W( 76 )
  JVS( 356) = W( 78 )
  JVS( 357) = W( 79 )
  JVS( 358) = W( 80 )
  JVS( 359) = W( 81 )
  JVS( 360) = W( 82 )
  JVS( 361) = W( 83 )
  JVS( 362) = W( 84 )
  JVS( 363) = W( 85 )
  JVS( 364) = W( 86 )
  JVS( 365) = W( 88 )
  JVS( 366) = W( 89 )
  JVS( 367) = W( 92 )
  JVS( 368) = W( 99 )
  JVS( 369) = W( 103 )
  IF ( ABS( JVS( 376 )) < TINY(a) ) THEN
         IER = 70
         RETURN
  END IF
   W( 46 ) = JVS( 370 )
   W( 60 ) = JVS( 371 )
   W( 64 ) = JVS( 372 )
   W( 66 ) = JVS( 373 )
   W( 67 ) = JVS( 374 )
   W( 68 ) = JVS( 375 )
   W( 70 ) = JVS( 376 )
   W( 74 ) = JVS( 377 )
   W( 78 ) = JVS( 378 )
   W( 81 ) = JVS( 379 )
   W( 82 ) = JVS( 380 )
   W( 83 ) = JVS( 381 )
   W( 84 ) = JVS( 382 )
   W( 85 ) = JVS( 383 )
   W( 86 ) = JVS( 384 )
   W( 89 ) = JVS( 385 )
   W( 92 ) = JVS( 386 )
   W( 96 ) = JVS( 387 )
   W( 99 ) = JVS( 388 )
   W( 101 ) = JVS( 389 )
   W( 103 ) = JVS( 390 )
  a = -W( 46 ) / JVS( 246 )
  W( 46 ) = -a
  W( 99 ) = W( 99 ) + a*JVS( 247 )
  W( 101 ) = W( 101 ) + a*JVS( 248 )
  a = -W( 60 ) / JVS( 296 )
  W( 60 ) = -a
  W( 99 ) = W( 99 ) + a*JVS( 297 )
  W( 103 ) = W( 103 ) + a*JVS( 298 )
  a = -W( 64 ) / JVS( 310 )
  W( 64 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 311 )
  W( 96 ) = W( 96 ) + a*JVS( 312 )
  W( 99 ) = W( 99 ) + a*JVS( 313 )
  W( 101 ) = W( 101 ) + a*JVS( 314 )
  a = -W( 66 ) / JVS( 322 )
  W( 66 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 323 )
  W( 92 ) = W( 92 ) + a*JVS( 324 )
  W( 99 ) = W( 99 ) + a*JVS( 325 )
  W( 103 ) = W( 103 ) + a*JVS( 326 )
  a = -W( 67 ) / JVS( 328 )
  W( 67 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 329 )
  W( 96 ) = W( 96 ) + a*JVS( 330 )
  W( 99 ) = W( 99 ) + a*JVS( 331 )
  W( 103 ) = W( 103 ) + a*JVS( 332 )
  a = -W( 68 ) / JVS( 338 )
  W( 68 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 339 )
  W( 82 ) = W( 82 ) + a*JVS( 340 )
  W( 85 ) = W( 85 ) + a*JVS( 341 )
  W( 92 ) = W( 92 ) + a*JVS( 342 )
  W( 99 ) = W( 99 ) + a*JVS( 343 )
  W( 103 ) = W( 103 ) + a*JVS( 344 )
  JVS( 370) = W( 46 )
  JVS( 371) = W( 60 )
  JVS( 372) = W( 64 )
  JVS( 373) = W( 66 )
  JVS( 374) = W( 67 )
  JVS( 375) = W( 68 )
  JVS( 376) = W( 70 )
  JVS( 377) = W( 74 )
  JVS( 378) = W( 78 )
  JVS( 379) = W( 81 )
  JVS( 380) = W( 82 )
  JVS( 381) = W( 83 )
  JVS( 382) = W( 84 )
  JVS( 383) = W( 85 )
  JVS( 384) = W( 86 )
  JVS( 385) = W( 89 )
  JVS( 386) = W( 92 )
  JVS( 387) = W( 96 )
  JVS( 388) = W( 99 )
  JVS( 389) = W( 101 )
  JVS( 390) = W( 103 )
  IF ( ABS( JVS( 391 )) < TINY(a) ) THEN
         IER = 71
         RETURN
  END IF
   W( 71 ) = JVS( 391 )
   W( 88 ) = JVS( 392 )
   W( 92 ) = JVS( 393 )
   W( 99 ) = JVS( 394 )
   W( 103 ) = JVS( 395 )
  JVS( 391) = W( 71 )
  JVS( 392) = W( 88 )
  JVS( 393) = W( 92 )
  JVS( 394) = W( 99 )
  JVS( 395) = W( 103 )
  IF ( ABS( JVS( 401 )) < TINY(a) ) THEN
         IER = 72
         RETURN
  END IF
   W( 43 ) = JVS( 396 )
   W( 48 ) = JVS( 397 )
   W( 49 ) = JVS( 398 )
   W( 50 ) = JVS( 399 )
   W( 61 ) = JVS( 400 )
   W( 72 ) = JVS( 401 )
   W( 76 ) = JVS( 402 )
   W( 79 ) = JVS( 403 )
   W( 80 ) = JVS( 404 )
   W( 84 ) = JVS( 405 )
   W( 87 ) = JVS( 406 )
   W( 92 ) = JVS( 407 )
   W( 99 ) = JVS( 408 )
   W( 101 ) = JVS( 409 )
   W( 103 ) = JVS( 410 )
  a = -W( 43 ) / JVS( 236 )
  W( 43 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 237 )
  a = -W( 48 ) / JVS( 252 )
  W( 48 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 253 )
  a = -W( 49 ) / JVS( 255 )
  W( 49 ) = -a
  W( 101 ) = W( 101 ) + a*JVS( 256 )
  W( 103 ) = W( 103 ) + a*JVS( 257 )
  a = -W( 50 ) / JVS( 258 )
  W( 50 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 259 )
  a = -W( 61 ) / JVS( 299 )
  W( 61 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 300 )
  JVS( 396) = W( 43 )
  JVS( 397) = W( 48 )
  JVS( 398) = W( 49 )
  JVS( 399) = W( 50 )
  JVS( 400) = W( 61 )
  JVS( 401) = W( 72 )
  JVS( 402) = W( 76 )
  JVS( 403) = W( 79 )
  JVS( 404) = W( 80 )
  JVS( 405) = W( 84 )
  JVS( 406) = W( 87 )
  JVS( 407) = W( 92 )
  JVS( 408) = W( 99 )
  JVS( 409) = W( 101 )
  JVS( 410) = W( 103 )
  IF ( ABS( JVS( 411 )) < TINY(a) ) THEN
         IER = 73
         RETURN
  END IF
   W( 73 ) = JVS( 411 )
   W( 88 ) = JVS( 412 )
   W( 92 ) = JVS( 413 )
   W( 99 ) = JVS( 414 )
   W( 103 ) = JVS( 415 )
  JVS( 411) = W( 73 )
  JVS( 412) = W( 88 )
  JVS( 413) = W( 92 )
  JVS( 414) = W( 99 )
  JVS( 415) = W( 103 )
  IF ( ABS( JVS( 418 )) < TINY(a) ) THEN
         IER = 74
         RETURN
  END IF
   W( 60 ) = JVS( 416 )
   W( 67 ) = JVS( 417 )
   W( 74 ) = JVS( 418 )
   W( 93 ) = JVS( 419 )
   W( 95 ) = JVS( 420 )
   W( 96 ) = JVS( 421 )
   W( 97 ) = JVS( 422 )
   W( 98 ) = JVS( 423 )
   W( 99 ) = JVS( 424 )
   W( 100 ) = JVS( 425 )
   W( 101 ) = JVS( 426 )
   W( 103 ) = JVS( 427 )
  a = -W( 60 ) / JVS( 296 )
  W( 60 ) = -a
  W( 99 ) = W( 99 ) + a*JVS( 297 )
  W( 103 ) = W( 103 ) + a*JVS( 298 )
  a = -W( 67 ) / JVS( 328 )
  W( 67 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 329 )
  W( 96 ) = W( 96 ) + a*JVS( 330 )
  W( 99 ) = W( 99 ) + a*JVS( 331 )
  W( 103 ) = W( 103 ) + a*JVS( 332 )
  JVS( 416) = W( 60 )
  JVS( 417) = W( 67 )
  JVS( 418) = W( 74 )
  JVS( 419) = W( 93 )
  JVS( 420) = W( 95 )
  JVS( 421) = W( 96 )
  JVS( 422) = W( 97 )
  JVS( 423) = W( 98 )
  JVS( 424) = W( 99 )
  JVS( 425) = W( 100 )
  JVS( 426) = W( 101 )
  JVS( 427) = W( 103 )
  IF ( ABS( JVS( 428 )) < TINY(a) ) THEN
         IER = 75
         RETURN
  END IF
   W( 75 ) = JVS( 428 )
   W( 88 ) = JVS( 429 )
   W( 92 ) = JVS( 430 )
   W( 99 ) = JVS( 431 )
   W( 103 ) = JVS( 432 )
  JVS( 428) = W( 75 )
  JVS( 429) = W( 88 )
  JVS( 430) = W( 92 )
  JVS( 431) = W( 99 )
  JVS( 432) = W( 103 )
  IF ( ABS( JVS( 433 )) < TINY(a) ) THEN
         IER = 76
         RETURN
  END IF
   W( 76 ) = JVS( 433 )
   W( 88 ) = JVS( 434 )
   W( 92 ) = JVS( 435 )
   W( 99 ) = JVS( 436 )
   W( 103 ) = JVS( 437 )
  JVS( 433) = W( 76 )
  JVS( 434) = W( 88 )
  JVS( 435) = W( 92 )
  JVS( 436) = W( 99 )
  JVS( 437) = W( 103 )
  IF ( ABS( JVS( 446 )) < TINY(a) ) THEN
         IER = 77
         RETURN
  END IF
   W( 48 ) = JVS( 438 )
   W( 50 ) = JVS( 439 )
   W( 58 ) = JVS( 440 )
   W( 59 ) = JVS( 441 )
   W( 61 ) = JVS( 442 )
   W( 72 ) = JVS( 443 )
   W( 75 ) = JVS( 444 )
   W( 76 ) = JVS( 445 )
   W( 77 ) = JVS( 446 )
   W( 79 ) = JVS( 447 )
   W( 80 ) = JVS( 448 )
   W( 84 ) = JVS( 449 )
   W( 85 ) = JVS( 450 )
   W( 87 ) = JVS( 451 )
   W( 88 ) = JVS( 452 )
   W( 90 ) = JVS( 453 )
   W( 91 ) = JVS( 454 )
   W( 92 ) = JVS( 455 )
   W( 93 ) = JVS( 456 )
   W( 94 ) = JVS( 457 )
   W( 95 ) = JVS( 458 )
   W( 96 ) = JVS( 459 )
   W( 97 ) = JVS( 460 )
   W( 98 ) = JVS( 461 )
   W( 99 ) = JVS( 462 )
   W( 100 ) = JVS( 463 )
   W( 101 ) = JVS( 464 )
   W( 102 ) = JVS( 465 )
   W( 103 ) = JVS( 466 )
   W( 104 ) = JVS( 467 )
  a = -W( 48 ) / JVS( 252 )
  W( 48 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 253 )
  a = -W( 50 ) / JVS( 258 )
  W( 50 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 259 )
  a = -W( 58 ) / JVS( 288 )
  W( 58 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 289 )
  a = -W( 59 ) / JVS( 292 )
  W( 59 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 293 )
  a = -W( 61 ) / JVS( 299 )
  W( 61 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 300 )
  a = -W( 72 ) / JVS( 401 )
  W( 72 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 402 )
  W( 79 ) = W( 79 ) + a*JVS( 403 )
  W( 80 ) = W( 80 ) + a*JVS( 404 )
  W( 84 ) = W( 84 ) + a*JVS( 405 )
  W( 87 ) = W( 87 ) + a*JVS( 406 )
  W( 92 ) = W( 92 ) + a*JVS( 407 )
  W( 99 ) = W( 99 ) + a*JVS( 408 )
  W( 101 ) = W( 101 ) + a*JVS( 409 )
  W( 103 ) = W( 103 ) + a*JVS( 410 )
  a = -W( 75 ) / JVS( 428 )
  W( 75 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 429 )
  W( 92 ) = W( 92 ) + a*JVS( 430 )
  W( 99 ) = W( 99 ) + a*JVS( 431 )
  W( 103 ) = W( 103 ) + a*JVS( 432 )
  a = -W( 76 ) / JVS( 433 )
  W( 76 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 434 )
  W( 92 ) = W( 92 ) + a*JVS( 435 )
  W( 99 ) = W( 99 ) + a*JVS( 436 )
  W( 103 ) = W( 103 ) + a*JVS( 437 )
  JVS( 438) = W( 48 )
  JVS( 439) = W( 50 )
  JVS( 440) = W( 58 )
  JVS( 441) = W( 59 )
  JVS( 442) = W( 61 )
  JVS( 443) = W( 72 )
  JVS( 444) = W( 75 )
  JVS( 445) = W( 76 )
  JVS( 446) = W( 77 )
  JVS( 447) = W( 79 )
  JVS( 448) = W( 80 )
  JVS( 449) = W( 84 )
  JVS( 450) = W( 85 )
  JVS( 451) = W( 87 )
  JVS( 452) = W( 88 )
  JVS( 453) = W( 90 )
  JVS( 454) = W( 91 )
  JVS( 455) = W( 92 )
  JVS( 456) = W( 93 )
  JVS( 457) = W( 94 )
  JVS( 458) = W( 95 )
  JVS( 459) = W( 96 )
  JVS( 460) = W( 97 )
  JVS( 461) = W( 98 )
  JVS( 462) = W( 99 )
  JVS( 463) = W( 100 )
  JVS( 464) = W( 101 )
  JVS( 465) = W( 102 )
  JVS( 466) = W( 103 )
  JVS( 467) = W( 104 )
  IF ( ABS( JVS( 469 )) < TINY(a) ) THEN
         IER = 78
         RETURN
  END IF
   W( 75 ) = JVS( 468 )
   W( 78 ) = JVS( 469 )
   W( 84 ) = JVS( 470 )
   W( 88 ) = JVS( 471 )
   W( 92 ) = JVS( 472 )
   W( 99 ) = JVS( 473 )
   W( 103 ) = JVS( 474 )
  a = -W( 75 ) / JVS( 428 )
  W( 75 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 429 )
  W( 92 ) = W( 92 ) + a*JVS( 430 )
  W( 99 ) = W( 99 ) + a*JVS( 431 )
  W( 103 ) = W( 103 ) + a*JVS( 432 )
  JVS( 468) = W( 75 )
  JVS( 469) = W( 78 )
  JVS( 470) = W( 84 )
  JVS( 471) = W( 88 )
  JVS( 472) = W( 92 )
  JVS( 473) = W( 99 )
  JVS( 474) = W( 103 )
  IF ( ABS( JVS( 475 )) < TINY(a) ) THEN
         IER = 79
         RETURN
  END IF
   W( 79 ) = JVS( 475 )
   W( 88 ) = JVS( 476 )
   W( 92 ) = JVS( 477 )
   W( 99 ) = JVS( 478 )
   W( 103 ) = JVS( 479 )
  JVS( 475) = W( 79 )
  JVS( 476) = W( 88 )
  JVS( 477) = W( 92 )
  JVS( 478) = W( 99 )
  JVS( 479) = W( 103 )
  IF ( ABS( JVS( 480 )) < TINY(a) ) THEN
         IER = 80
         RETURN
  END IF
   W( 80 ) = JVS( 480 )
   W( 88 ) = JVS( 481 )
   W( 92 ) = JVS( 482 )
   W( 99 ) = JVS( 483 )
   W( 103 ) = JVS( 484 )
  JVS( 480) = W( 80 )
  JVS( 481) = W( 88 )
  JVS( 482) = W( 92 )
  JVS( 483) = W( 99 )
  JVS( 484) = W( 103 )
  IF ( ABS( JVS( 496 )) < TINY(a) ) THEN
         IER = 81
         RETURN
  END IF
   W( 51 ) = JVS( 485 )
   W( 57 ) = JVS( 486 )
   W( 58 ) = JVS( 487 )
   W( 59 ) = JVS( 488 )
   W( 62 ) = JVS( 489 )
   W( 63 ) = JVS( 490 )
   W( 67 ) = JVS( 491 )
   W( 71 ) = JVS( 492 )
   W( 74 ) = JVS( 493 )
   W( 79 ) = JVS( 494 )
   W( 80 ) = JVS( 495 )
   W( 81 ) = JVS( 496 )
   W( 82 ) = JVS( 497 )
   W( 88 ) = JVS( 498 )
   W( 92 ) = JVS( 499 )
   W( 93 ) = JVS( 500 )
   W( 95 ) = JVS( 501 )
   W( 96 ) = JVS( 502 )
   W( 97 ) = JVS( 503 )
   W( 98 ) = JVS( 504 )
   W( 99 ) = JVS( 505 )
   W( 100 ) = JVS( 506 )
   W( 101 ) = JVS( 507 )
   W( 103 ) = JVS( 508 )
  a = -W( 51 ) / JVS( 260 )
  W( 51 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 261 )
  a = -W( 57 ) / JVS( 284 )
  W( 57 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 285 )
  a = -W( 58 ) / JVS( 288 )
  W( 58 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 289 )
  a = -W( 59 ) / JVS( 292 )
  W( 59 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 293 )
  a = -W( 62 ) / JVS( 303 )
  W( 62 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 304 )
  W( 103 ) = W( 103 ) + a*JVS( 305 )
  a = -W( 63 ) / JVS( 306 )
  W( 63 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 307 )
  W( 103 ) = W( 103 ) + a*JVS( 308 )
  a = -W( 67 ) / JVS( 328 )
  W( 67 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 329 )
  W( 96 ) = W( 96 ) + a*JVS( 330 )
  W( 99 ) = W( 99 ) + a*JVS( 331 )
  W( 103 ) = W( 103 ) + a*JVS( 332 )
  a = -W( 71 ) / JVS( 391 )
  W( 71 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 392 )
  W( 92 ) = W( 92 ) + a*JVS( 393 )
  W( 99 ) = W( 99 ) + a*JVS( 394 )
  W( 103 ) = W( 103 ) + a*JVS( 395 )
  a = -W( 74 ) / JVS( 418 )
  W( 74 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 419 )
  W( 95 ) = W( 95 ) + a*JVS( 420 )
  W( 96 ) = W( 96 ) + a*JVS( 421 )
  W( 97 ) = W( 97 ) + a*JVS( 422 )
  W( 98 ) = W( 98 ) + a*JVS( 423 )
  W( 99 ) = W( 99 ) + a*JVS( 424 )
  W( 100 ) = W( 100 ) + a*JVS( 425 )
  W( 101 ) = W( 101 ) + a*JVS( 426 )
  W( 103 ) = W( 103 ) + a*JVS( 427 )
  a = -W( 79 ) / JVS( 475 )
  W( 79 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 92 ) = W( 92 ) + a*JVS( 477 )
  W( 99 ) = W( 99 ) + a*JVS( 478 )
  W( 103 ) = W( 103 ) + a*JVS( 479 )
  a = -W( 80 ) / JVS( 480 )
  W( 80 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 92 ) = W( 92 ) + a*JVS( 482 )
  W( 99 ) = W( 99 ) + a*JVS( 483 )
  W( 103 ) = W( 103 ) + a*JVS( 484 )
  JVS( 485) = W( 51 )
  JVS( 486) = W( 57 )
  JVS( 487) = W( 58 )
  JVS( 488) = W( 59 )
  JVS( 489) = W( 62 )
  JVS( 490) = W( 63 )
  JVS( 491) = W( 67 )
  JVS( 492) = W( 71 )
  JVS( 493) = W( 74 )
  JVS( 494) = W( 79 )
  JVS( 495) = W( 80 )
  JVS( 496) = W( 81 )
  JVS( 497) = W( 82 )
  JVS( 498) = W( 88 )
  JVS( 499) = W( 92 )
  JVS( 500) = W( 93 )
  JVS( 501) = W( 95 )
  JVS( 502) = W( 96 )
  JVS( 503) = W( 97 )
  JVS( 504) = W( 98 )
  JVS( 505) = W( 99 )
  JVS( 506) = W( 100 )
  JVS( 507) = W( 101 )
  JVS( 508) = W( 103 )
  IF ( ABS( JVS( 510 )) < TINY(a) ) THEN
         IER = 82
         RETURN
  END IF
   W( 75 ) = JVS( 509 )
   W( 82 ) = JVS( 510 )
   W( 84 ) = JVS( 511 )
   W( 88 ) = JVS( 512 )
   W( 92 ) = JVS( 513 )
   W( 99 ) = JVS( 514 )
   W( 103 ) = JVS( 515 )
  a = -W( 75 ) / JVS( 428 )
  W( 75 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 429 )
  W( 92 ) = W( 92 ) + a*JVS( 430 )
  W( 99 ) = W( 99 ) + a*JVS( 431 )
  W( 103 ) = W( 103 ) + a*JVS( 432 )
  JVS( 509) = W( 75 )
  JVS( 510) = W( 82 )
  JVS( 511) = W( 84 )
  JVS( 512) = W( 88 )
  JVS( 513) = W( 92 )
  JVS( 514) = W( 99 )
  JVS( 515) = W( 103 )
  IF ( ABS( JVS( 525 )) < TINY(a) ) THEN
         IER = 83
         RETURN
  END IF
   W( 37 ) = JVS( 516 )
   W( 44 ) = JVS( 517 )
   W( 48 ) = JVS( 518 )
   W( 50 ) = JVS( 519 )
   W( 61 ) = JVS( 520 )
   W( 71 ) = JVS( 521 )
   W( 73 ) = JVS( 522 )
   W( 76 ) = JVS( 523 )
   W( 82 ) = JVS( 524 )
   W( 83 ) = JVS( 525 )
   W( 84 ) = JVS( 526 )
   W( 87 ) = JVS( 527 )
   W( 88 ) = JVS( 528 )
   W( 89 ) = JVS( 529 )
   W( 90 ) = JVS( 530 )
   W( 91 ) = JVS( 531 )
   W( 92 ) = JVS( 532 )
   W( 93 ) = JVS( 533 )
   W( 95 ) = JVS( 534 )
   W( 97 ) = JVS( 535 )
   W( 98 ) = JVS( 536 )
   W( 99 ) = JVS( 537 )
   W( 100 ) = JVS( 538 )
   W( 103 ) = JVS( 539 )
  a = -W( 37 ) / JVS( 219 )
  W( 37 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 220 )
  a = -W( 44 ) / JVS( 238 )
  W( 44 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 239 )
  a = -W( 48 ) / JVS( 252 )
  W( 48 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 253 )
  a = -W( 50 ) / JVS( 258 )
  W( 50 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 259 )
  a = -W( 61 ) / JVS( 299 )
  W( 61 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 300 )
  a = -W( 71 ) / JVS( 391 )
  W( 71 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 392 )
  W( 92 ) = W( 92 ) + a*JVS( 393 )
  W( 99 ) = W( 99 ) + a*JVS( 394 )
  W( 103 ) = W( 103 ) + a*JVS( 395 )
  a = -W( 73 ) / JVS( 411 )
  W( 73 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 92 ) = W( 92 ) + a*JVS( 413 )
  W( 99 ) = W( 99 ) + a*JVS( 414 )
  W( 103 ) = W( 103 ) + a*JVS( 415 )
  a = -W( 76 ) / JVS( 433 )
  W( 76 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 434 )
  W( 92 ) = W( 92 ) + a*JVS( 435 )
  W( 99 ) = W( 99 ) + a*JVS( 436 )
  W( 103 ) = W( 103 ) + a*JVS( 437 )
  a = -W( 82 ) / JVS( 510 )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 511 )
  W( 88 ) = W( 88 ) + a*JVS( 512 )
  W( 92 ) = W( 92 ) + a*JVS( 513 )
  W( 99 ) = W( 99 ) + a*JVS( 514 )
  W( 103 ) = W( 103 ) + a*JVS( 515 )
  JVS( 516) = W( 37 )
  JVS( 517) = W( 44 )
  JVS( 518) = W( 48 )
  JVS( 519) = W( 50 )
  JVS( 520) = W( 61 )
  JVS( 521) = W( 71 )
  JVS( 522) = W( 73 )
  JVS( 523) = W( 76 )
  JVS( 524) = W( 82 )
  JVS( 525) = W( 83 )
  JVS( 526) = W( 84 )
  JVS( 527) = W( 87 )
  JVS( 528) = W( 88 )
  JVS( 529) = W( 89 )
  JVS( 530) = W( 90 )
  JVS( 531) = W( 91 )
  JVS( 532) = W( 92 )
  JVS( 533) = W( 93 )
  JVS( 534) = W( 95 )
  JVS( 535) = W( 97 )
  JVS( 536) = W( 98 )
  JVS( 537) = W( 99 )
  JVS( 538) = W( 100 )
  JVS( 539) = W( 103 )
  IF ( ABS( JVS( 540 )) < TINY(a) ) THEN
         IER = 84
         RETURN
  END IF
   W( 84 ) = JVS( 540 )
   W( 88 ) = JVS( 541 )
   W( 92 ) = JVS( 542 )
   W( 99 ) = JVS( 543 )
   W( 103 ) = JVS( 544 )
  JVS( 540) = W( 84 )
  JVS( 541) = W( 88 )
  JVS( 542) = W( 92 )
  JVS( 543) = W( 99 )
  JVS( 544) = W( 103 )
  IF ( ABS( JVS( 547 )) < TINY(a) ) THEN
         IER = 85
         RETURN
  END IF
   W( 75 ) = JVS( 545 )
   W( 84 ) = JVS( 546 )
   W( 85 ) = JVS( 547 )
   W( 88 ) = JVS( 548 )
   W( 92 ) = JVS( 549 )
   W( 99 ) = JVS( 550 )
   W( 103 ) = JVS( 551 )
  a = -W( 75 ) / JVS( 428 )
  W( 75 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 429 )
  W( 92 ) = W( 92 ) + a*JVS( 430 )
  W( 99 ) = W( 99 ) + a*JVS( 431 )
  W( 103 ) = W( 103 ) + a*JVS( 432 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  JVS( 545) = W( 75 )
  JVS( 546) = W( 84 )
  JVS( 547) = W( 85 )
  JVS( 548) = W( 88 )
  JVS( 549) = W( 92 )
  JVS( 550) = W( 99 )
  JVS( 551) = W( 103 )
  IF ( ABS( JVS( 572 )) < TINY(a) ) THEN
         IER = 86
         RETURN
  END IF
   W( 44 ) = JVS( 552 )
   W( 48 ) = JVS( 553 )
   W( 50 ) = JVS( 554 )
   W( 52 ) = JVS( 555 )
   W( 53 ) = JVS( 556 )
   W( 56 ) = JVS( 557 )
   W( 61 ) = JVS( 558 )
   W( 63 ) = JVS( 559 )
   W( 71 ) = JVS( 560 )
   W( 72 ) = JVS( 561 )
   W( 73 ) = JVS( 562 )
   W( 75 ) = JVS( 563 )
   W( 76 ) = JVS( 564 )
   W( 78 ) = JVS( 565 )
   W( 79 ) = JVS( 566 )
   W( 80 ) = JVS( 567 )
   W( 81 ) = JVS( 568 )
   W( 82 ) = JVS( 569 )
   W( 84 ) = JVS( 570 )
   W( 85 ) = JVS( 571 )
   W( 86 ) = JVS( 572 )
   W( 87 ) = JVS( 573 )
   W( 88 ) = JVS( 574 )
   W( 90 ) = JVS( 575 )
   W( 91 ) = JVS( 576 )
   W( 92 ) = JVS( 577 )
   W( 93 ) = JVS( 578 )
   W( 94 ) = JVS( 579 )
   W( 95 ) = JVS( 580 )
   W( 96 ) = JVS( 581 )
   W( 97 ) = JVS( 582 )
   W( 98 ) = JVS( 583 )
   W( 99 ) = JVS( 584 )
   W( 100 ) = JVS( 585 )
   W( 101 ) = JVS( 586 )
   W( 102 ) = JVS( 587 )
   W( 103 ) = JVS( 588 )
   W( 104 ) = JVS( 589 )
  a = -W( 44 ) / JVS( 238 )
  W( 44 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 239 )
  a = -W( 48 ) / JVS( 252 )
  W( 48 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 253 )
  a = -W( 50 ) / JVS( 258 )
  W( 50 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 259 )
  a = -W( 52 ) / JVS( 262 )
  W( 52 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 263 )
  W( 102 ) = W( 102 ) + a*JVS( 264 )
  W( 103 ) = W( 103 ) + a*JVS( 265 )
  a = -W( 53 ) / JVS( 266 )
  W( 53 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 267 )
  W( 93 ) = W( 93 ) + a*JVS( 268 )
  W( 96 ) = W( 96 ) + a*JVS( 269 )
  a = -W( 56 ) / JVS( 279 )
  W( 56 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 280 )
  W( 102 ) = W( 102 ) + a*JVS( 281 )
  W( 103 ) = W( 103 ) + a*JVS( 282 )
  W( 104 ) = W( 104 ) + a*JVS( 283 )
  a = -W( 61 ) / JVS( 299 )
  W( 61 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 300 )
  a = -W( 63 ) / JVS( 306 )
  W( 63 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 307 )
  W( 103 ) = W( 103 ) + a*JVS( 308 )
  a = -W( 71 ) / JVS( 391 )
  W( 71 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 392 )
  W( 92 ) = W( 92 ) + a*JVS( 393 )
  W( 99 ) = W( 99 ) + a*JVS( 394 )
  W( 103 ) = W( 103 ) + a*JVS( 395 )
  a = -W( 72 ) / JVS( 401 )
  W( 72 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 402 )
  W( 79 ) = W( 79 ) + a*JVS( 403 )
  W( 80 ) = W( 80 ) + a*JVS( 404 )
  W( 84 ) = W( 84 ) + a*JVS( 405 )
  W( 87 ) = W( 87 ) + a*JVS( 406 )
  W( 92 ) = W( 92 ) + a*JVS( 407 )
  W( 99 ) = W( 99 ) + a*JVS( 408 )
  W( 101 ) = W( 101 ) + a*JVS( 409 )
  W( 103 ) = W( 103 ) + a*JVS( 410 )
  a = -W( 73 ) / JVS( 411 )
  W( 73 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 92 ) = W( 92 ) + a*JVS( 413 )
  W( 99 ) = W( 99 ) + a*JVS( 414 )
  W( 103 ) = W( 103 ) + a*JVS( 415 )
  a = -W( 75 ) / JVS( 428 )
  W( 75 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 429 )
  W( 92 ) = W( 92 ) + a*JVS( 430 )
  W( 99 ) = W( 99 ) + a*JVS( 431 )
  W( 103 ) = W( 103 ) + a*JVS( 432 )
  a = -W( 76 ) / JVS( 433 )
  W( 76 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 434 )
  W( 92 ) = W( 92 ) + a*JVS( 435 )
  W( 99 ) = W( 99 ) + a*JVS( 436 )
  W( 103 ) = W( 103 ) + a*JVS( 437 )
  a = -W( 78 ) / JVS( 469 )
  W( 78 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 470 )
  W( 88 ) = W( 88 ) + a*JVS( 471 )
  W( 92 ) = W( 92 ) + a*JVS( 472 )
  W( 99 ) = W( 99 ) + a*JVS( 473 )
  W( 103 ) = W( 103 ) + a*JVS( 474 )
  a = -W( 79 ) / JVS( 475 )
  W( 79 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 92 ) = W( 92 ) + a*JVS( 477 )
  W( 99 ) = W( 99 ) + a*JVS( 478 )
  W( 103 ) = W( 103 ) + a*JVS( 479 )
  a = -W( 80 ) / JVS( 480 )
  W( 80 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 92 ) = W( 92 ) + a*JVS( 482 )
  W( 99 ) = W( 99 ) + a*JVS( 483 )
  W( 103 ) = W( 103 ) + a*JVS( 484 )
  a = -W( 81 ) / JVS( 496 )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 497 )
  W( 88 ) = W( 88 ) + a*JVS( 498 )
  W( 92 ) = W( 92 ) + a*JVS( 499 )
  W( 93 ) = W( 93 ) + a*JVS( 500 )
  W( 95 ) = W( 95 ) + a*JVS( 501 )
  W( 96 ) = W( 96 ) + a*JVS( 502 )
  W( 97 ) = W( 97 ) + a*JVS( 503 )
  W( 98 ) = W( 98 ) + a*JVS( 504 )
  W( 99 ) = W( 99 ) + a*JVS( 505 )
  W( 100 ) = W( 100 ) + a*JVS( 506 )
  W( 101 ) = W( 101 ) + a*JVS( 507 )
  W( 103 ) = W( 103 ) + a*JVS( 508 )
  a = -W( 82 ) / JVS( 510 )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 511 )
  W( 88 ) = W( 88 ) + a*JVS( 512 )
  W( 92 ) = W( 92 ) + a*JVS( 513 )
  W( 99 ) = W( 99 ) + a*JVS( 514 )
  W( 103 ) = W( 103 ) + a*JVS( 515 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  a = -W( 85 ) / JVS( 547 )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 548 )
  W( 92 ) = W( 92 ) + a*JVS( 549 )
  W( 99 ) = W( 99 ) + a*JVS( 550 )
  W( 103 ) = W( 103 ) + a*JVS( 551 )
  JVS( 552) = W( 44 )
  JVS( 553) = W( 48 )
  JVS( 554) = W( 50 )
  JVS( 555) = W( 52 )
  JVS( 556) = W( 53 )
  JVS( 557) = W( 56 )
  JVS( 558) = W( 61 )
  JVS( 559) = W( 63 )
  JVS( 560) = W( 71 )
  JVS( 561) = W( 72 )
  JVS( 562) = W( 73 )
  JVS( 563) = W( 75 )
  JVS( 564) = W( 76 )
  JVS( 565) = W( 78 )
  JVS( 566) = W( 79 )
  JVS( 567) = W( 80 )
  JVS( 568) = W( 81 )
  JVS( 569) = W( 82 )
  JVS( 570) = W( 84 )
  JVS( 571) = W( 85 )
  JVS( 572) = W( 86 )
  JVS( 573) = W( 87 )
  JVS( 574) = W( 88 )
  JVS( 575) = W( 90 )
  JVS( 576) = W( 91 )
  JVS( 577) = W( 92 )
  JVS( 578) = W( 93 )
  JVS( 579) = W( 94 )
  JVS( 580) = W( 95 )
  JVS( 581) = W( 96 )
  JVS( 582) = W( 97 )
  JVS( 583) = W( 98 )
  JVS( 584) = W( 99 )
  JVS( 585) = W( 100 )
  JVS( 586) = W( 101 )
  JVS( 587) = W( 102 )
  JVS( 588) = W( 103 )
  JVS( 589) = W( 104 )
  IF ( ABS( JVS( 596 )) < TINY(a) ) THEN
         IER = 87
         RETURN
  END IF
   W( 49 ) = JVS( 590 )
   W( 76 ) = JVS( 591 )
   W( 79 ) = JVS( 592 )
   W( 80 ) = JVS( 593 )
   W( 82 ) = JVS( 594 )
   W( 84 ) = JVS( 595 )
   W( 87 ) = JVS( 596 )
   W( 88 ) = JVS( 597 )
   W( 92 ) = JVS( 598 )
   W( 93 ) = JVS( 599 )
   W( 99 ) = JVS( 600 )
   W( 101 ) = JVS( 601 )
   W( 103 ) = JVS( 602 )
   W( 104 ) = JVS( 603 )
  a = -W( 49 ) / JVS( 255 )
  W( 49 ) = -a
  W( 101 ) = W( 101 ) + a*JVS( 256 )
  W( 103 ) = W( 103 ) + a*JVS( 257 )
  a = -W( 76 ) / JVS( 433 )
  W( 76 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 434 )
  W( 92 ) = W( 92 ) + a*JVS( 435 )
  W( 99 ) = W( 99 ) + a*JVS( 436 )
  W( 103 ) = W( 103 ) + a*JVS( 437 )
  a = -W( 79 ) / JVS( 475 )
  W( 79 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 92 ) = W( 92 ) + a*JVS( 477 )
  W( 99 ) = W( 99 ) + a*JVS( 478 )
  W( 103 ) = W( 103 ) + a*JVS( 479 )
  a = -W( 80 ) / JVS( 480 )
  W( 80 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 92 ) = W( 92 ) + a*JVS( 482 )
  W( 99 ) = W( 99 ) + a*JVS( 483 )
  W( 103 ) = W( 103 ) + a*JVS( 484 )
  a = -W( 82 ) / JVS( 510 )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 511 )
  W( 88 ) = W( 88 ) + a*JVS( 512 )
  W( 92 ) = W( 92 ) + a*JVS( 513 )
  W( 99 ) = W( 99 ) + a*JVS( 514 )
  W( 103 ) = W( 103 ) + a*JVS( 515 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  JVS( 590) = W( 49 )
  JVS( 591) = W( 76 )
  JVS( 592) = W( 79 )
  JVS( 593) = W( 80 )
  JVS( 594) = W( 82 )
  JVS( 595) = W( 84 )
  JVS( 596) = W( 87 )
  JVS( 597) = W( 88 )
  JVS( 598) = W( 92 )
  JVS( 599) = W( 93 )
  JVS( 600) = W( 99 )
  JVS( 601) = W( 101 )
  JVS( 602) = W( 103 )
  JVS( 603) = W( 104 )
  IF ( ABS( JVS( 614 )) < TINY(a) ) THEN
         IER = 88
         RETURN
  END IF
   W( 30 ) = JVS( 604 )
   W( 71 ) = JVS( 605 )
   W( 73 ) = JVS( 606 )
   W( 75 ) = JVS( 607 )
   W( 76 ) = JVS( 608 )
   W( 78 ) = JVS( 609 )
   W( 79 ) = JVS( 610 )
   W( 80 ) = JVS( 611 )
   W( 84 ) = JVS( 612 )
   W( 85 ) = JVS( 613 )
   W( 88 ) = JVS( 614 )
   W( 92 ) = JVS( 615 )
   W( 93 ) = JVS( 616 )
   W( 99 ) = JVS( 617 )
   W( 101 ) = JVS( 618 )
   W( 103 ) = JVS( 619 )
  a = -W( 30 ) / JVS( 203 )
  W( 30 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 204 )
  a = -W( 71 ) / JVS( 391 )
  W( 71 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 392 )
  W( 92 ) = W( 92 ) + a*JVS( 393 )
  W( 99 ) = W( 99 ) + a*JVS( 394 )
  W( 103 ) = W( 103 ) + a*JVS( 395 )
  a = -W( 73 ) / JVS( 411 )
  W( 73 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 92 ) = W( 92 ) + a*JVS( 413 )
  W( 99 ) = W( 99 ) + a*JVS( 414 )
  W( 103 ) = W( 103 ) + a*JVS( 415 )
  a = -W( 75 ) / JVS( 428 )
  W( 75 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 429 )
  W( 92 ) = W( 92 ) + a*JVS( 430 )
  W( 99 ) = W( 99 ) + a*JVS( 431 )
  W( 103 ) = W( 103 ) + a*JVS( 432 )
  a = -W( 76 ) / JVS( 433 )
  W( 76 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 434 )
  W( 92 ) = W( 92 ) + a*JVS( 435 )
  W( 99 ) = W( 99 ) + a*JVS( 436 )
  W( 103 ) = W( 103 ) + a*JVS( 437 )
  a = -W( 78 ) / JVS( 469 )
  W( 78 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 470 )
  W( 88 ) = W( 88 ) + a*JVS( 471 )
  W( 92 ) = W( 92 ) + a*JVS( 472 )
  W( 99 ) = W( 99 ) + a*JVS( 473 )
  W( 103 ) = W( 103 ) + a*JVS( 474 )
  a = -W( 79 ) / JVS( 475 )
  W( 79 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 92 ) = W( 92 ) + a*JVS( 477 )
  W( 99 ) = W( 99 ) + a*JVS( 478 )
  W( 103 ) = W( 103 ) + a*JVS( 479 )
  a = -W( 80 ) / JVS( 480 )
  W( 80 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 92 ) = W( 92 ) + a*JVS( 482 )
  W( 99 ) = W( 99 ) + a*JVS( 483 )
  W( 103 ) = W( 103 ) + a*JVS( 484 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  a = -W( 85 ) / JVS( 547 )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 548 )
  W( 92 ) = W( 92 ) + a*JVS( 549 )
  W( 99 ) = W( 99 ) + a*JVS( 550 )
  W( 103 ) = W( 103 ) + a*JVS( 551 )
  JVS( 604) = W( 30 )
  JVS( 605) = W( 71 )
  JVS( 606) = W( 73 )
  JVS( 607) = W( 75 )
  JVS( 608) = W( 76 )
  JVS( 609) = W( 78 )
  JVS( 610) = W( 79 )
  JVS( 611) = W( 80 )
  JVS( 612) = W( 84 )
  JVS( 613) = W( 85 )
  JVS( 614) = W( 88 )
  JVS( 615) = W( 92 )
  JVS( 616) = W( 93 )
  JVS( 617) = W( 99 )
  JVS( 618) = W( 101 )
  JVS( 619) = W( 103 )
  IF ( ABS( JVS( 639 )) < TINY(a) ) THEN
         IER = 89
         RETURN
  END IF
   W( 43 ) = JVS( 620 )
   W( 48 ) = JVS( 621 )
   W( 50 ) = JVS( 622 )
   W( 58 ) = JVS( 623 )
   W( 59 ) = JVS( 624 )
   W( 61 ) = JVS( 625 )
   W( 62 ) = JVS( 626 )
   W( 65 ) = JVS( 627 )
   W( 71 ) = JVS( 628 )
   W( 73 ) = JVS( 629 )
   W( 76 ) = JVS( 630 )
   W( 78 ) = JVS( 631 )
   W( 79 ) = JVS( 632 )
   W( 80 ) = JVS( 633 )
   W( 82 ) = JVS( 634 )
   W( 84 ) = JVS( 635 )
   W( 85 ) = JVS( 636 )
   W( 87 ) = JVS( 637 )
   W( 88 ) = JVS( 638 )
   W( 89 ) = JVS( 639 )
   W( 90 ) = JVS( 640 )
   W( 91 ) = JVS( 641 )
   W( 92 ) = JVS( 642 )
   W( 93 ) = JVS( 643 )
   W( 94 ) = JVS( 644 )
   W( 96 ) = JVS( 645 )
   W( 99 ) = JVS( 646 )
   W( 101 ) = JVS( 647 )
   W( 103 ) = JVS( 648 )
   W( 104 ) = JVS( 649 )
  a = -W( 43 ) / JVS( 236 )
  W( 43 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 237 )
  a = -W( 48 ) / JVS( 252 )
  W( 48 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 253 )
  a = -W( 50 ) / JVS( 258 )
  W( 50 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 259 )
  a = -W( 58 ) / JVS( 288 )
  W( 58 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 289 )
  a = -W( 59 ) / JVS( 292 )
  W( 59 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 293 )
  a = -W( 61 ) / JVS( 299 )
  W( 61 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 300 )
  a = -W( 62 ) / JVS( 303 )
  W( 62 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 304 )
  W( 103 ) = W( 103 ) + a*JVS( 305 )
  a = -W( 65 ) / JVS( 315 )
  W( 65 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 316 )
  W( 96 ) = W( 96 ) + a*JVS( 317 )
  W( 103 ) = W( 103 ) + a*JVS( 318 )
  W( 104 ) = W( 104 ) + a*JVS( 319 )
  a = -W( 71 ) / JVS( 391 )
  W( 71 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 392 )
  W( 92 ) = W( 92 ) + a*JVS( 393 )
  W( 99 ) = W( 99 ) + a*JVS( 394 )
  W( 103 ) = W( 103 ) + a*JVS( 395 )
  a = -W( 73 ) / JVS( 411 )
  W( 73 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 92 ) = W( 92 ) + a*JVS( 413 )
  W( 99 ) = W( 99 ) + a*JVS( 414 )
  W( 103 ) = W( 103 ) + a*JVS( 415 )
  a = -W( 76 ) / JVS( 433 )
  W( 76 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 434 )
  W( 92 ) = W( 92 ) + a*JVS( 435 )
  W( 99 ) = W( 99 ) + a*JVS( 436 )
  W( 103 ) = W( 103 ) + a*JVS( 437 )
  a = -W( 78 ) / JVS( 469 )
  W( 78 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 470 )
  W( 88 ) = W( 88 ) + a*JVS( 471 )
  W( 92 ) = W( 92 ) + a*JVS( 472 )
  W( 99 ) = W( 99 ) + a*JVS( 473 )
  W( 103 ) = W( 103 ) + a*JVS( 474 )
  a = -W( 79 ) / JVS( 475 )
  W( 79 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 92 ) = W( 92 ) + a*JVS( 477 )
  W( 99 ) = W( 99 ) + a*JVS( 478 )
  W( 103 ) = W( 103 ) + a*JVS( 479 )
  a = -W( 80 ) / JVS( 480 )
  W( 80 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 92 ) = W( 92 ) + a*JVS( 482 )
  W( 99 ) = W( 99 ) + a*JVS( 483 )
  W( 103 ) = W( 103 ) + a*JVS( 484 )
  a = -W( 82 ) / JVS( 510 )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 511 )
  W( 88 ) = W( 88 ) + a*JVS( 512 )
  W( 92 ) = W( 92 ) + a*JVS( 513 )
  W( 99 ) = W( 99 ) + a*JVS( 514 )
  W( 103 ) = W( 103 ) + a*JVS( 515 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  a = -W( 85 ) / JVS( 547 )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 548 )
  W( 92 ) = W( 92 ) + a*JVS( 549 )
  W( 99 ) = W( 99 ) + a*JVS( 550 )
  W( 103 ) = W( 103 ) + a*JVS( 551 )
  a = -W( 87 ) / JVS( 596 )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 597 )
  W( 92 ) = W( 92 ) + a*JVS( 598 )
  W( 93 ) = W( 93 ) + a*JVS( 599 )
  W( 99 ) = W( 99 ) + a*JVS( 600 )
  W( 101 ) = W( 101 ) + a*JVS( 601 )
  W( 103 ) = W( 103 ) + a*JVS( 602 )
  W( 104 ) = W( 104 ) + a*JVS( 603 )
  a = -W( 88 ) / JVS( 614 )
  W( 88 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 615 )
  W( 93 ) = W( 93 ) + a*JVS( 616 )
  W( 99 ) = W( 99 ) + a*JVS( 617 )
  W( 101 ) = W( 101 ) + a*JVS( 618 )
  W( 103 ) = W( 103 ) + a*JVS( 619 )
  JVS( 620) = W( 43 )
  JVS( 621) = W( 48 )
  JVS( 622) = W( 50 )
  JVS( 623) = W( 58 )
  JVS( 624) = W( 59 )
  JVS( 625) = W( 61 )
  JVS( 626) = W( 62 )
  JVS( 627) = W( 65 )
  JVS( 628) = W( 71 )
  JVS( 629) = W( 73 )
  JVS( 630) = W( 76 )
  JVS( 631) = W( 78 )
  JVS( 632) = W( 79 )
  JVS( 633) = W( 80 )
  JVS( 634) = W( 82 )
  JVS( 635) = W( 84 )
  JVS( 636) = W( 85 )
  JVS( 637) = W( 87 )
  JVS( 638) = W( 88 )
  JVS( 639) = W( 89 )
  JVS( 640) = W( 90 )
  JVS( 641) = W( 91 )
  JVS( 642) = W( 92 )
  JVS( 643) = W( 93 )
  JVS( 644) = W( 94 )
  JVS( 645) = W( 96 )
  JVS( 646) = W( 99 )
  JVS( 647) = W( 101 )
  JVS( 648) = W( 103 )
  JVS( 649) = W( 104 )
  IF ( ABS( JVS( 661 )) < TINY(a) ) THEN
         IER = 90
         RETURN
  END IF
   W( 48 ) = JVS( 650 )
   W( 50 ) = JVS( 651 )
   W( 61 ) = JVS( 652 )
   W( 73 ) = JVS( 653 )
   W( 76 ) = JVS( 654 )
   W( 78 ) = JVS( 655 )
   W( 82 ) = JVS( 656 )
   W( 84 ) = JVS( 657 )
   W( 85 ) = JVS( 658 )
   W( 87 ) = JVS( 659 )
   W( 88 ) = JVS( 660 )
   W( 90 ) = JVS( 661 )
   W( 91 ) = JVS( 662 )
   W( 92 ) = JVS( 663 )
   W( 93 ) = JVS( 664 )
   W( 94 ) = JVS( 665 )
   W( 99 ) = JVS( 666 )
   W( 101 ) = JVS( 667 )
   W( 102 ) = JVS( 668 )
   W( 103 ) = JVS( 669 )
   W( 104 ) = JVS( 670 )
  a = -W( 48 ) / JVS( 252 )
  W( 48 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 253 )
  a = -W( 50 ) / JVS( 258 )
  W( 50 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 259 )
  a = -W( 61 ) / JVS( 299 )
  W( 61 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 300 )
  a = -W( 73 ) / JVS( 411 )
  W( 73 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 92 ) = W( 92 ) + a*JVS( 413 )
  W( 99 ) = W( 99 ) + a*JVS( 414 )
  W( 103 ) = W( 103 ) + a*JVS( 415 )
  a = -W( 76 ) / JVS( 433 )
  W( 76 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 434 )
  W( 92 ) = W( 92 ) + a*JVS( 435 )
  W( 99 ) = W( 99 ) + a*JVS( 436 )
  W( 103 ) = W( 103 ) + a*JVS( 437 )
  a = -W( 78 ) / JVS( 469 )
  W( 78 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 470 )
  W( 88 ) = W( 88 ) + a*JVS( 471 )
  W( 92 ) = W( 92 ) + a*JVS( 472 )
  W( 99 ) = W( 99 ) + a*JVS( 473 )
  W( 103 ) = W( 103 ) + a*JVS( 474 )
  a = -W( 82 ) / JVS( 510 )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 511 )
  W( 88 ) = W( 88 ) + a*JVS( 512 )
  W( 92 ) = W( 92 ) + a*JVS( 513 )
  W( 99 ) = W( 99 ) + a*JVS( 514 )
  W( 103 ) = W( 103 ) + a*JVS( 515 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  a = -W( 85 ) / JVS( 547 )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 548 )
  W( 92 ) = W( 92 ) + a*JVS( 549 )
  W( 99 ) = W( 99 ) + a*JVS( 550 )
  W( 103 ) = W( 103 ) + a*JVS( 551 )
  a = -W( 87 ) / JVS( 596 )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 597 )
  W( 92 ) = W( 92 ) + a*JVS( 598 )
  W( 93 ) = W( 93 ) + a*JVS( 599 )
  W( 99 ) = W( 99 ) + a*JVS( 600 )
  W( 101 ) = W( 101 ) + a*JVS( 601 )
  W( 103 ) = W( 103 ) + a*JVS( 602 )
  W( 104 ) = W( 104 ) + a*JVS( 603 )
  a = -W( 88 ) / JVS( 614 )
  W( 88 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 615 )
  W( 93 ) = W( 93 ) + a*JVS( 616 )
  W( 99 ) = W( 99 ) + a*JVS( 617 )
  W( 101 ) = W( 101 ) + a*JVS( 618 )
  W( 103 ) = W( 103 ) + a*JVS( 619 )
  JVS( 650) = W( 48 )
  JVS( 651) = W( 50 )
  JVS( 652) = W( 61 )
  JVS( 653) = W( 73 )
  JVS( 654) = W( 76 )
  JVS( 655) = W( 78 )
  JVS( 656) = W( 82 )
  JVS( 657) = W( 84 )
  JVS( 658) = W( 85 )
  JVS( 659) = W( 87 )
  JVS( 660) = W( 88 )
  JVS( 661) = W( 90 )
  JVS( 662) = W( 91 )
  JVS( 663) = W( 92 )
  JVS( 664) = W( 93 )
  JVS( 665) = W( 94 )
  JVS( 666) = W( 99 )
  JVS( 667) = W( 101 )
  JVS( 668) = W( 102 )
  JVS( 669) = W( 103 )
  JVS( 670) = W( 104 )
  IF ( ABS( JVS( 683 )) < TINY(a) ) THEN
         IER = 91
         RETURN
  END IF
   W( 50 ) = JVS( 671 )
   W( 57 ) = JVS( 672 )
   W( 61 ) = JVS( 673 )
   W( 75 ) = JVS( 674 )
   W( 76 ) = JVS( 675 )
   W( 79 ) = JVS( 676 )
   W( 80 ) = JVS( 677 )
   W( 82 ) = JVS( 678 )
   W( 84 ) = JVS( 679 )
   W( 85 ) = JVS( 680 )
   W( 87 ) = JVS( 681 )
   W( 88 ) = JVS( 682 )
   W( 91 ) = JVS( 683 )
   W( 92 ) = JVS( 684 )
   W( 93 ) = JVS( 685 )
   W( 94 ) = JVS( 686 )
   W( 95 ) = JVS( 687 )
   W( 97 ) = JVS( 688 )
   W( 99 ) = JVS( 689 )
   W( 100 ) = JVS( 690 )
   W( 101 ) = JVS( 691 )
   W( 102 ) = JVS( 692 )
   W( 103 ) = JVS( 693 )
   W( 104 ) = JVS( 694 )
  a = -W( 50 ) / JVS( 258 )
  W( 50 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 259 )
  a = -W( 57 ) / JVS( 284 )
  W( 57 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 285 )
  a = -W( 61 ) / JVS( 299 )
  W( 61 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 300 )
  a = -W( 75 ) / JVS( 428 )
  W( 75 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 429 )
  W( 92 ) = W( 92 ) + a*JVS( 430 )
  W( 99 ) = W( 99 ) + a*JVS( 431 )
  W( 103 ) = W( 103 ) + a*JVS( 432 )
  a = -W( 76 ) / JVS( 433 )
  W( 76 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 434 )
  W( 92 ) = W( 92 ) + a*JVS( 435 )
  W( 99 ) = W( 99 ) + a*JVS( 436 )
  W( 103 ) = W( 103 ) + a*JVS( 437 )
  a = -W( 79 ) / JVS( 475 )
  W( 79 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 92 ) = W( 92 ) + a*JVS( 477 )
  W( 99 ) = W( 99 ) + a*JVS( 478 )
  W( 103 ) = W( 103 ) + a*JVS( 479 )
  a = -W( 80 ) / JVS( 480 )
  W( 80 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 92 ) = W( 92 ) + a*JVS( 482 )
  W( 99 ) = W( 99 ) + a*JVS( 483 )
  W( 103 ) = W( 103 ) + a*JVS( 484 )
  a = -W( 82 ) / JVS( 510 )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 511 )
  W( 88 ) = W( 88 ) + a*JVS( 512 )
  W( 92 ) = W( 92 ) + a*JVS( 513 )
  W( 99 ) = W( 99 ) + a*JVS( 514 )
  W( 103 ) = W( 103 ) + a*JVS( 515 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  a = -W( 85 ) / JVS( 547 )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 548 )
  W( 92 ) = W( 92 ) + a*JVS( 549 )
  W( 99 ) = W( 99 ) + a*JVS( 550 )
  W( 103 ) = W( 103 ) + a*JVS( 551 )
  a = -W( 87 ) / JVS( 596 )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 597 )
  W( 92 ) = W( 92 ) + a*JVS( 598 )
  W( 93 ) = W( 93 ) + a*JVS( 599 )
  W( 99 ) = W( 99 ) + a*JVS( 600 )
  W( 101 ) = W( 101 ) + a*JVS( 601 )
  W( 103 ) = W( 103 ) + a*JVS( 602 )
  W( 104 ) = W( 104 ) + a*JVS( 603 )
  a = -W( 88 ) / JVS( 614 )
  W( 88 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 615 )
  W( 93 ) = W( 93 ) + a*JVS( 616 )
  W( 99 ) = W( 99 ) + a*JVS( 617 )
  W( 101 ) = W( 101 ) + a*JVS( 618 )
  W( 103 ) = W( 103 ) + a*JVS( 619 )
  JVS( 671) = W( 50 )
  JVS( 672) = W( 57 )
  JVS( 673) = W( 61 )
  JVS( 674) = W( 75 )
  JVS( 675) = W( 76 )
  JVS( 676) = W( 79 )
  JVS( 677) = W( 80 )
  JVS( 678) = W( 82 )
  JVS( 679) = W( 84 )
  JVS( 680) = W( 85 )
  JVS( 681) = W( 87 )
  JVS( 682) = W( 88 )
  JVS( 683) = W( 91 )
  JVS( 684) = W( 92 )
  JVS( 685) = W( 93 )
  JVS( 686) = W( 94 )
  JVS( 687) = W( 95 )
  JVS( 688) = W( 97 )
  JVS( 689) = W( 99 )
  JVS( 690) = W( 100 )
  JVS( 691) = W( 101 )
  JVS( 692) = W( 102 )
  JVS( 693) = W( 103 )
  JVS( 694) = W( 104 )
  IF ( ABS( JVS( 708 )) < TINY(a) ) THEN
         IER = 92
         RETURN
  END IF
   W( 62 ) = JVS( 695 )
   W( 63 ) = JVS( 696 )
   W( 71 ) = JVS( 697 )
   W( 73 ) = JVS( 698 )
   W( 75 ) = JVS( 699 )
   W( 76 ) = JVS( 700 )
   W( 78 ) = JVS( 701 )
   W( 79 ) = JVS( 702 )
   W( 80 ) = JVS( 703 )
   W( 82 ) = JVS( 704 )
   W( 84 ) = JVS( 705 )
   W( 85 ) = JVS( 706 )
   W( 88 ) = JVS( 707 )
   W( 92 ) = JVS( 708 )
   W( 93 ) = JVS( 709 )
   W( 95 ) = JVS( 710 )
   W( 96 ) = JVS( 711 )
   W( 97 ) = JVS( 712 )
   W( 98 ) = JVS( 713 )
   W( 99 ) = JVS( 714 )
   W( 100 ) = JVS( 715 )
   W( 101 ) = JVS( 716 )
   W( 103 ) = JVS( 717 )
  a = -W( 62 ) / JVS( 303 )
  W( 62 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 304 )
  W( 103 ) = W( 103 ) + a*JVS( 305 )
  a = -W( 63 ) / JVS( 306 )
  W( 63 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 307 )
  W( 103 ) = W( 103 ) + a*JVS( 308 )
  a = -W( 71 ) / JVS( 391 )
  W( 71 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 392 )
  W( 92 ) = W( 92 ) + a*JVS( 393 )
  W( 99 ) = W( 99 ) + a*JVS( 394 )
  W( 103 ) = W( 103 ) + a*JVS( 395 )
  a = -W( 73 ) / JVS( 411 )
  W( 73 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 92 ) = W( 92 ) + a*JVS( 413 )
  W( 99 ) = W( 99 ) + a*JVS( 414 )
  W( 103 ) = W( 103 ) + a*JVS( 415 )
  a = -W( 75 ) / JVS( 428 )
  W( 75 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 429 )
  W( 92 ) = W( 92 ) + a*JVS( 430 )
  W( 99 ) = W( 99 ) + a*JVS( 431 )
  W( 103 ) = W( 103 ) + a*JVS( 432 )
  a = -W( 76 ) / JVS( 433 )
  W( 76 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 434 )
  W( 92 ) = W( 92 ) + a*JVS( 435 )
  W( 99 ) = W( 99 ) + a*JVS( 436 )
  W( 103 ) = W( 103 ) + a*JVS( 437 )
  a = -W( 78 ) / JVS( 469 )
  W( 78 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 470 )
  W( 88 ) = W( 88 ) + a*JVS( 471 )
  W( 92 ) = W( 92 ) + a*JVS( 472 )
  W( 99 ) = W( 99 ) + a*JVS( 473 )
  W( 103 ) = W( 103 ) + a*JVS( 474 )
  a = -W( 79 ) / JVS( 475 )
  W( 79 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 92 ) = W( 92 ) + a*JVS( 477 )
  W( 99 ) = W( 99 ) + a*JVS( 478 )
  W( 103 ) = W( 103 ) + a*JVS( 479 )
  a = -W( 80 ) / JVS( 480 )
  W( 80 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 92 ) = W( 92 ) + a*JVS( 482 )
  W( 99 ) = W( 99 ) + a*JVS( 483 )
  W( 103 ) = W( 103 ) + a*JVS( 484 )
  a = -W( 82 ) / JVS( 510 )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 511 )
  W( 88 ) = W( 88 ) + a*JVS( 512 )
  W( 92 ) = W( 92 ) + a*JVS( 513 )
  W( 99 ) = W( 99 ) + a*JVS( 514 )
  W( 103 ) = W( 103 ) + a*JVS( 515 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  a = -W( 85 ) / JVS( 547 )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 548 )
  W( 92 ) = W( 92 ) + a*JVS( 549 )
  W( 99 ) = W( 99 ) + a*JVS( 550 )
  W( 103 ) = W( 103 ) + a*JVS( 551 )
  a = -W( 88 ) / JVS( 614 )
  W( 88 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 615 )
  W( 93 ) = W( 93 ) + a*JVS( 616 )
  W( 99 ) = W( 99 ) + a*JVS( 617 )
  W( 101 ) = W( 101 ) + a*JVS( 618 )
  W( 103 ) = W( 103 ) + a*JVS( 619 )
  JVS( 695) = W( 62 )
  JVS( 696) = W( 63 )
  JVS( 697) = W( 71 )
  JVS( 698) = W( 73 )
  JVS( 699) = W( 75 )
  JVS( 700) = W( 76 )
  JVS( 701) = W( 78 )
  JVS( 702) = W( 79 )
  JVS( 703) = W( 80 )
  JVS( 704) = W( 82 )
  JVS( 705) = W( 84 )
  JVS( 706) = W( 85 )
  JVS( 707) = W( 88 )
  JVS( 708) = W( 92 )
  JVS( 709) = W( 93 )
  JVS( 710) = W( 95 )
  JVS( 711) = W( 96 )
  JVS( 712) = W( 97 )
  JVS( 713) = W( 98 )
  JVS( 714) = W( 99 )
  JVS( 715) = W( 100 )
  JVS( 716) = W( 101 )
  JVS( 717) = W( 103 )
  IF ( ABS( JVS( 731 )) < TINY(a) ) THEN
         IER = 93
         RETURN
  END IF
   W( 47 ) = JVS( 718 )
   W( 53 ) = JVS( 719 )
   W( 77 ) = JVS( 720 )
   W( 79 ) = JVS( 721 )
   W( 80 ) = JVS( 722 )
   W( 84 ) = JVS( 723 )
   W( 85 ) = JVS( 724 )
   W( 86 ) = JVS( 725 )
   W( 87 ) = JVS( 726 )
   W( 88 ) = JVS( 727 )
   W( 90 ) = JVS( 728 )
   W( 91 ) = JVS( 729 )
   W( 92 ) = JVS( 730 )
   W( 93 ) = JVS( 731 )
   W( 94 ) = JVS( 732 )
   W( 95 ) = JVS( 733 )
   W( 96 ) = JVS( 734 )
   W( 97 ) = JVS( 735 )
   W( 98 ) = JVS( 736 )
   W( 99 ) = JVS( 737 )
   W( 100 ) = JVS( 738 )
   W( 101 ) = JVS( 739 )
   W( 102 ) = JVS( 740 )
   W( 103 ) = JVS( 741 )
   W( 104 ) = JVS( 742 )
  a = -W( 47 ) / JVS( 249 )
  W( 47 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 250 )
  W( 103 ) = W( 103 ) + a*JVS( 251 )
  a = -W( 53 ) / JVS( 266 )
  W( 53 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 267 )
  W( 93 ) = W( 93 ) + a*JVS( 268 )
  W( 96 ) = W( 96 ) + a*JVS( 269 )
  a = -W( 77 ) / JVS( 446 )
  W( 77 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 447 )
  W( 80 ) = W( 80 ) + a*JVS( 448 )
  W( 84 ) = W( 84 ) + a*JVS( 449 )
  W( 85 ) = W( 85 ) + a*JVS( 450 )
  W( 87 ) = W( 87 ) + a*JVS( 451 )
  W( 88 ) = W( 88 ) + a*JVS( 452 )
  W( 90 ) = W( 90 ) + a*JVS( 453 )
  W( 91 ) = W( 91 ) + a*JVS( 454 )
  W( 92 ) = W( 92 ) + a*JVS( 455 )
  W( 93 ) = W( 93 ) + a*JVS( 456 )
  W( 94 ) = W( 94 ) + a*JVS( 457 )
  W( 95 ) = W( 95 ) + a*JVS( 458 )
  W( 96 ) = W( 96 ) + a*JVS( 459 )
  W( 97 ) = W( 97 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  W( 99 ) = W( 99 ) + a*JVS( 462 )
  W( 100 ) = W( 100 ) + a*JVS( 463 )
  W( 101 ) = W( 101 ) + a*JVS( 464 )
  W( 102 ) = W( 102 ) + a*JVS( 465 )
  W( 103 ) = W( 103 ) + a*JVS( 466 )
  W( 104 ) = W( 104 ) + a*JVS( 467 )
  a = -W( 79 ) / JVS( 475 )
  W( 79 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 92 ) = W( 92 ) + a*JVS( 477 )
  W( 99 ) = W( 99 ) + a*JVS( 478 )
  W( 103 ) = W( 103 ) + a*JVS( 479 )
  a = -W( 80 ) / JVS( 480 )
  W( 80 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 92 ) = W( 92 ) + a*JVS( 482 )
  W( 99 ) = W( 99 ) + a*JVS( 483 )
  W( 103 ) = W( 103 ) + a*JVS( 484 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  a = -W( 85 ) / JVS( 547 )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 548 )
  W( 92 ) = W( 92 ) + a*JVS( 549 )
  W( 99 ) = W( 99 ) + a*JVS( 550 )
  W( 103 ) = W( 103 ) + a*JVS( 551 )
  a = -W( 86 ) / JVS( 572 )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 573 )
  W( 88 ) = W( 88 ) + a*JVS( 574 )
  W( 90 ) = W( 90 ) + a*JVS( 575 )
  W( 91 ) = W( 91 ) + a*JVS( 576 )
  W( 92 ) = W( 92 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 94 ) = W( 94 ) + a*JVS( 579 )
  W( 95 ) = W( 95 ) + a*JVS( 580 )
  W( 96 ) = W( 96 ) + a*JVS( 581 )
  W( 97 ) = W( 97 ) + a*JVS( 582 )
  W( 98 ) = W( 98 ) + a*JVS( 583 )
  W( 99 ) = W( 99 ) + a*JVS( 584 )
  W( 100 ) = W( 100 ) + a*JVS( 585 )
  W( 101 ) = W( 101 ) + a*JVS( 586 )
  W( 102 ) = W( 102 ) + a*JVS( 587 )
  W( 103 ) = W( 103 ) + a*JVS( 588 )
  W( 104 ) = W( 104 ) + a*JVS( 589 )
  a = -W( 87 ) / JVS( 596 )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 597 )
  W( 92 ) = W( 92 ) + a*JVS( 598 )
  W( 93 ) = W( 93 ) + a*JVS( 599 )
  W( 99 ) = W( 99 ) + a*JVS( 600 )
  W( 101 ) = W( 101 ) + a*JVS( 601 )
  W( 103 ) = W( 103 ) + a*JVS( 602 )
  W( 104 ) = W( 104 ) + a*JVS( 603 )
  a = -W( 88 ) / JVS( 614 )
  W( 88 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 615 )
  W( 93 ) = W( 93 ) + a*JVS( 616 )
  W( 99 ) = W( 99 ) + a*JVS( 617 )
  W( 101 ) = W( 101 ) + a*JVS( 618 )
  W( 103 ) = W( 103 ) + a*JVS( 619 )
  a = -W( 90 ) / JVS( 661 )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 662 )
  W( 92 ) = W( 92 ) + a*JVS( 663 )
  W( 93 ) = W( 93 ) + a*JVS( 664 )
  W( 94 ) = W( 94 ) + a*JVS( 665 )
  W( 99 ) = W( 99 ) + a*JVS( 666 )
  W( 101 ) = W( 101 ) + a*JVS( 667 )
  W( 102 ) = W( 102 ) + a*JVS( 668 )
  W( 103 ) = W( 103 ) + a*JVS( 669 )
  W( 104 ) = W( 104 ) + a*JVS( 670 )
  a = -W( 91 ) / JVS( 683 )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 684 )
  W( 93 ) = W( 93 ) + a*JVS( 685 )
  W( 94 ) = W( 94 ) + a*JVS( 686 )
  W( 95 ) = W( 95 ) + a*JVS( 687 )
  W( 97 ) = W( 97 ) + a*JVS( 688 )
  W( 99 ) = W( 99 ) + a*JVS( 689 )
  W( 100 ) = W( 100 ) + a*JVS( 690 )
  W( 101 ) = W( 101 ) + a*JVS( 691 )
  W( 102 ) = W( 102 ) + a*JVS( 692 )
  W( 103 ) = W( 103 ) + a*JVS( 693 )
  W( 104 ) = W( 104 ) + a*JVS( 694 )
  a = -W( 92 ) / JVS( 708 )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 709 )
  W( 95 ) = W( 95 ) + a*JVS( 710 )
  W( 96 ) = W( 96 ) + a*JVS( 711 )
  W( 97 ) = W( 97 ) + a*JVS( 712 )
  W( 98 ) = W( 98 ) + a*JVS( 713 )
  W( 99 ) = W( 99 ) + a*JVS( 714 )
  W( 100 ) = W( 100 ) + a*JVS( 715 )
  W( 101 ) = W( 101 ) + a*JVS( 716 )
  W( 103 ) = W( 103 ) + a*JVS( 717 )
  JVS( 718) = W( 47 )
  JVS( 719) = W( 53 )
  JVS( 720) = W( 77 )
  JVS( 721) = W( 79 )
  JVS( 722) = W( 80 )
  JVS( 723) = W( 84 )
  JVS( 724) = W( 85 )
  JVS( 725) = W( 86 )
  JVS( 726) = W( 87 )
  JVS( 727) = W( 88 )
  JVS( 728) = W( 90 )
  JVS( 729) = W( 91 )
  JVS( 730) = W( 92 )
  JVS( 731) = W( 93 )
  JVS( 732) = W( 94 )
  JVS( 733) = W( 95 )
  JVS( 734) = W( 96 )
  JVS( 735) = W( 97 )
  JVS( 736) = W( 98 )
  JVS( 737) = W( 99 )
  JVS( 738) = W( 100 )
  JVS( 739) = W( 101 )
  JVS( 740) = W( 102 )
  JVS( 741) = W( 103 )
  JVS( 742) = W( 104 )
  IF ( ABS( JVS( 776 )) < TINY(a) ) THEN
         IER = 94
         RETURN
  END IF
   W( 37 ) = JVS( 743 )
   W( 43 ) = JVS( 744 )
   W( 44 ) = JVS( 745 )
   W( 48 ) = JVS( 746 )
   W( 50 ) = JVS( 747 )
   W( 51 ) = JVS( 748 )
   W( 57 ) = JVS( 749 )
   W( 58 ) = JVS( 750 )
   W( 59 ) = JVS( 751 )
   W( 60 ) = JVS( 752 )
   W( 61 ) = JVS( 753 )
   W( 62 ) = JVS( 754 )
   W( 63 ) = JVS( 755 )
   W( 65 ) = JVS( 756 )
   W( 67 ) = JVS( 757 )
   W( 71 ) = JVS( 758 )
   W( 73 ) = JVS( 759 )
   W( 74 ) = JVS( 760 )
   W( 75 ) = JVS( 761 )
   W( 76 ) = JVS( 762 )
   W( 78 ) = JVS( 763 )
   W( 79 ) = JVS( 764 )
   W( 80 ) = JVS( 765 )
   W( 82 ) = JVS( 766 )
   W( 84 ) = JVS( 767 )
   W( 85 ) = JVS( 768 )
   W( 87 ) = JVS( 769 )
   W( 88 ) = JVS( 770 )
   W( 89 ) = JVS( 771 )
   W( 90 ) = JVS( 772 )
   W( 91 ) = JVS( 773 )
   W( 92 ) = JVS( 774 )
   W( 93 ) = JVS( 775 )
   W( 94 ) = JVS( 776 )
   W( 95 ) = JVS( 777 )
   W( 96 ) = JVS( 778 )
   W( 97 ) = JVS( 779 )
   W( 98 ) = JVS( 780 )
   W( 99 ) = JVS( 781 )
   W( 100 ) = JVS( 782 )
   W( 101 ) = JVS( 783 )
   W( 102 ) = JVS( 784 )
   W( 103 ) = JVS( 785 )
   W( 104 ) = JVS( 786 )
  a = -W( 37 ) / JVS( 219 )
  W( 37 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 220 )
  a = -W( 43 ) / JVS( 236 )
  W( 43 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 237 )
  a = -W( 44 ) / JVS( 238 )
  W( 44 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 239 )
  a = -W( 48 ) / JVS( 252 )
  W( 48 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 253 )
  a = -W( 50 ) / JVS( 258 )
  W( 50 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 259 )
  a = -W( 51 ) / JVS( 260 )
  W( 51 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 261 )
  a = -W( 57 ) / JVS( 284 )
  W( 57 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 285 )
  a = -W( 58 ) / JVS( 288 )
  W( 58 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 289 )
  a = -W( 59 ) / JVS( 292 )
  W( 59 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 293 )
  a = -W( 60 ) / JVS( 296 )
  W( 60 ) = -a
  W( 99 ) = W( 99 ) + a*JVS( 297 )
  W( 103 ) = W( 103 ) + a*JVS( 298 )
  a = -W( 61 ) / JVS( 299 )
  W( 61 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 300 )
  a = -W( 62 ) / JVS( 303 )
  W( 62 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 304 )
  W( 103 ) = W( 103 ) + a*JVS( 305 )
  a = -W( 63 ) / JVS( 306 )
  W( 63 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 307 )
  W( 103 ) = W( 103 ) + a*JVS( 308 )
  a = -W( 65 ) / JVS( 315 )
  W( 65 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 316 )
  W( 96 ) = W( 96 ) + a*JVS( 317 )
  W( 103 ) = W( 103 ) + a*JVS( 318 )
  W( 104 ) = W( 104 ) + a*JVS( 319 )
  a = -W( 67 ) / JVS( 328 )
  W( 67 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 329 )
  W( 96 ) = W( 96 ) + a*JVS( 330 )
  W( 99 ) = W( 99 ) + a*JVS( 331 )
  W( 103 ) = W( 103 ) + a*JVS( 332 )
  a = -W( 71 ) / JVS( 391 )
  W( 71 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 392 )
  W( 92 ) = W( 92 ) + a*JVS( 393 )
  W( 99 ) = W( 99 ) + a*JVS( 394 )
  W( 103 ) = W( 103 ) + a*JVS( 395 )
  a = -W( 73 ) / JVS( 411 )
  W( 73 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 92 ) = W( 92 ) + a*JVS( 413 )
  W( 99 ) = W( 99 ) + a*JVS( 414 )
  W( 103 ) = W( 103 ) + a*JVS( 415 )
  a = -W( 74 ) / JVS( 418 )
  W( 74 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 419 )
  W( 95 ) = W( 95 ) + a*JVS( 420 )
  W( 96 ) = W( 96 ) + a*JVS( 421 )
  W( 97 ) = W( 97 ) + a*JVS( 422 )
  W( 98 ) = W( 98 ) + a*JVS( 423 )
  W( 99 ) = W( 99 ) + a*JVS( 424 )
  W( 100 ) = W( 100 ) + a*JVS( 425 )
  W( 101 ) = W( 101 ) + a*JVS( 426 )
  W( 103 ) = W( 103 ) + a*JVS( 427 )
  a = -W( 75 ) / JVS( 428 )
  W( 75 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 429 )
  W( 92 ) = W( 92 ) + a*JVS( 430 )
  W( 99 ) = W( 99 ) + a*JVS( 431 )
  W( 103 ) = W( 103 ) + a*JVS( 432 )
  a = -W( 76 ) / JVS( 433 )
  W( 76 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 434 )
  W( 92 ) = W( 92 ) + a*JVS( 435 )
  W( 99 ) = W( 99 ) + a*JVS( 436 )
  W( 103 ) = W( 103 ) + a*JVS( 437 )
  a = -W( 78 ) / JVS( 469 )
  W( 78 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 470 )
  W( 88 ) = W( 88 ) + a*JVS( 471 )
  W( 92 ) = W( 92 ) + a*JVS( 472 )
  W( 99 ) = W( 99 ) + a*JVS( 473 )
  W( 103 ) = W( 103 ) + a*JVS( 474 )
  a = -W( 79 ) / JVS( 475 )
  W( 79 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 92 ) = W( 92 ) + a*JVS( 477 )
  W( 99 ) = W( 99 ) + a*JVS( 478 )
  W( 103 ) = W( 103 ) + a*JVS( 479 )
  a = -W( 80 ) / JVS( 480 )
  W( 80 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 92 ) = W( 92 ) + a*JVS( 482 )
  W( 99 ) = W( 99 ) + a*JVS( 483 )
  W( 103 ) = W( 103 ) + a*JVS( 484 )
  a = -W( 82 ) / JVS( 510 )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 511 )
  W( 88 ) = W( 88 ) + a*JVS( 512 )
  W( 92 ) = W( 92 ) + a*JVS( 513 )
  W( 99 ) = W( 99 ) + a*JVS( 514 )
  W( 103 ) = W( 103 ) + a*JVS( 515 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  a = -W( 85 ) / JVS( 547 )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 548 )
  W( 92 ) = W( 92 ) + a*JVS( 549 )
  W( 99 ) = W( 99 ) + a*JVS( 550 )
  W( 103 ) = W( 103 ) + a*JVS( 551 )
  a = -W( 87 ) / JVS( 596 )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 597 )
  W( 92 ) = W( 92 ) + a*JVS( 598 )
  W( 93 ) = W( 93 ) + a*JVS( 599 )
  W( 99 ) = W( 99 ) + a*JVS( 600 )
  W( 101 ) = W( 101 ) + a*JVS( 601 )
  W( 103 ) = W( 103 ) + a*JVS( 602 )
  W( 104 ) = W( 104 ) + a*JVS( 603 )
  a = -W( 88 ) / JVS( 614 )
  W( 88 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 615 )
  W( 93 ) = W( 93 ) + a*JVS( 616 )
  W( 99 ) = W( 99 ) + a*JVS( 617 )
  W( 101 ) = W( 101 ) + a*JVS( 618 )
  W( 103 ) = W( 103 ) + a*JVS( 619 )
  a = -W( 89 ) / JVS( 639 )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 640 )
  W( 91 ) = W( 91 ) + a*JVS( 641 )
  W( 92 ) = W( 92 ) + a*JVS( 642 )
  W( 93 ) = W( 93 ) + a*JVS( 643 )
  W( 94 ) = W( 94 ) + a*JVS( 644 )
  W( 96 ) = W( 96 ) + a*JVS( 645 )
  W( 99 ) = W( 99 ) + a*JVS( 646 )
  W( 101 ) = W( 101 ) + a*JVS( 647 )
  W( 103 ) = W( 103 ) + a*JVS( 648 )
  W( 104 ) = W( 104 ) + a*JVS( 649 )
  a = -W( 90 ) / JVS( 661 )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 662 )
  W( 92 ) = W( 92 ) + a*JVS( 663 )
  W( 93 ) = W( 93 ) + a*JVS( 664 )
  W( 94 ) = W( 94 ) + a*JVS( 665 )
  W( 99 ) = W( 99 ) + a*JVS( 666 )
  W( 101 ) = W( 101 ) + a*JVS( 667 )
  W( 102 ) = W( 102 ) + a*JVS( 668 )
  W( 103 ) = W( 103 ) + a*JVS( 669 )
  W( 104 ) = W( 104 ) + a*JVS( 670 )
  a = -W( 91 ) / JVS( 683 )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 684 )
  W( 93 ) = W( 93 ) + a*JVS( 685 )
  W( 94 ) = W( 94 ) + a*JVS( 686 )
  W( 95 ) = W( 95 ) + a*JVS( 687 )
  W( 97 ) = W( 97 ) + a*JVS( 688 )
  W( 99 ) = W( 99 ) + a*JVS( 689 )
  W( 100 ) = W( 100 ) + a*JVS( 690 )
  W( 101 ) = W( 101 ) + a*JVS( 691 )
  W( 102 ) = W( 102 ) + a*JVS( 692 )
  W( 103 ) = W( 103 ) + a*JVS( 693 )
  W( 104 ) = W( 104 ) + a*JVS( 694 )
  a = -W( 92 ) / JVS( 708 )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 709 )
  W( 95 ) = W( 95 ) + a*JVS( 710 )
  W( 96 ) = W( 96 ) + a*JVS( 711 )
  W( 97 ) = W( 97 ) + a*JVS( 712 )
  W( 98 ) = W( 98 ) + a*JVS( 713 )
  W( 99 ) = W( 99 ) + a*JVS( 714 )
  W( 100 ) = W( 100 ) + a*JVS( 715 )
  W( 101 ) = W( 101 ) + a*JVS( 716 )
  W( 103 ) = W( 103 ) + a*JVS( 717 )
  a = -W( 93 ) / JVS( 731 )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 732 )
  W( 95 ) = W( 95 ) + a*JVS( 733 )
  W( 96 ) = W( 96 ) + a*JVS( 734 )
  W( 97 ) = W( 97 ) + a*JVS( 735 )
  W( 98 ) = W( 98 ) + a*JVS( 736 )
  W( 99 ) = W( 99 ) + a*JVS( 737 )
  W( 100 ) = W( 100 ) + a*JVS( 738 )
  W( 101 ) = W( 101 ) + a*JVS( 739 )
  W( 102 ) = W( 102 ) + a*JVS( 740 )
  W( 103 ) = W( 103 ) + a*JVS( 741 )
  W( 104 ) = W( 104 ) + a*JVS( 742 )
  JVS( 743) = W( 37 )
  JVS( 744) = W( 43 )
  JVS( 745) = W( 44 )
  JVS( 746) = W( 48 )
  JVS( 747) = W( 50 )
  JVS( 748) = W( 51 )
  JVS( 749) = W( 57 )
  JVS( 750) = W( 58 )
  JVS( 751) = W( 59 )
  JVS( 752) = W( 60 )
  JVS( 753) = W( 61 )
  JVS( 754) = W( 62 )
  JVS( 755) = W( 63 )
  JVS( 756) = W( 65 )
  JVS( 757) = W( 67 )
  JVS( 758) = W( 71 )
  JVS( 759) = W( 73 )
  JVS( 760) = W( 74 )
  JVS( 761) = W( 75 )
  JVS( 762) = W( 76 )
  JVS( 763) = W( 78 )
  JVS( 764) = W( 79 )
  JVS( 765) = W( 80 )
  JVS( 766) = W( 82 )
  JVS( 767) = W( 84 )
  JVS( 768) = W( 85 )
  JVS( 769) = W( 87 )
  JVS( 770) = W( 88 )
  JVS( 771) = W( 89 )
  JVS( 772) = W( 90 )
  JVS( 773) = W( 91 )
  JVS( 774) = W( 92 )
  JVS( 775) = W( 93 )
  JVS( 776) = W( 94 )
  JVS( 777) = W( 95 )
  JVS( 778) = W( 96 )
  JVS( 779) = W( 97 )
  JVS( 780) = W( 98 )
  JVS( 781) = W( 99 )
  JVS( 782) = W( 100 )
  JVS( 783) = W( 101 )
  JVS( 784) = W( 102 )
  JVS( 785) = W( 103 )
  JVS( 786) = W( 104 )
  IF ( ABS( JVS( 811 )) < TINY(a) ) THEN
         IER = 95
         RETURN
  END IF
   W( 38 ) = JVS( 787 )
   W( 45 ) = JVS( 788 )
   W( 51 ) = JVS( 789 )
   W( 58 ) = JVS( 790 )
   W( 59 ) = JVS( 791 )
   W( 61 ) = JVS( 792 )
   W( 68 ) = JVS( 793 )
   W( 72 ) = JVS( 794 )
   W( 76 ) = JVS( 795 )
   W( 78 ) = JVS( 796 )
   W( 79 ) = JVS( 797 )
   W( 80 ) = JVS( 798 )
   W( 82 ) = JVS( 799 )
   W( 83 ) = JVS( 800 )
   W( 84 ) = JVS( 801 )
   W( 85 ) = JVS( 802 )
   W( 87 ) = JVS( 803 )
   W( 88 ) = JVS( 804 )
   W( 89 ) = JVS( 805 )
   W( 90 ) = JVS( 806 )
   W( 91 ) = JVS( 807 )
   W( 92 ) = JVS( 808 )
   W( 93 ) = JVS( 809 )
   W( 94 ) = JVS( 810 )
   W( 95 ) = JVS( 811 )
   W( 96 ) = JVS( 812 )
   W( 97 ) = JVS( 813 )
   W( 98 ) = JVS( 814 )
   W( 99 ) = JVS( 815 )
   W( 100 ) = JVS( 816 )
   W( 101 ) = JVS( 817 )
   W( 102 ) = JVS( 818 )
   W( 103 ) = JVS( 819 )
   W( 104 ) = JVS( 820 )
  a = -W( 38 ) / JVS( 221 )
  W( 38 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 222 )
  W( 101 ) = W( 101 ) + a*JVS( 223 )
  a = -W( 45 ) / JVS( 240 )
  W( 45 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 241 )
  W( 79 ) = W( 79 ) + a*JVS( 242 )
  W( 80 ) = W( 80 ) + a*JVS( 243 )
  W( 92 ) = W( 92 ) + a*JVS( 244 )
  W( 103 ) = W( 103 ) + a*JVS( 245 )
  a = -W( 51 ) / JVS( 260 )
  W( 51 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 261 )
  a = -W( 58 ) / JVS( 288 )
  W( 58 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 289 )
  a = -W( 59 ) / JVS( 292 )
  W( 59 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 293 )
  a = -W( 61 ) / JVS( 299 )
  W( 61 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 300 )
  a = -W( 68 ) / JVS( 338 )
  W( 68 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 339 )
  W( 82 ) = W( 82 ) + a*JVS( 340 )
  W( 85 ) = W( 85 ) + a*JVS( 341 )
  W( 92 ) = W( 92 ) + a*JVS( 342 )
  W( 99 ) = W( 99 ) + a*JVS( 343 )
  W( 103 ) = W( 103 ) + a*JVS( 344 )
  a = -W( 72 ) / JVS( 401 )
  W( 72 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 402 )
  W( 79 ) = W( 79 ) + a*JVS( 403 )
  W( 80 ) = W( 80 ) + a*JVS( 404 )
  W( 84 ) = W( 84 ) + a*JVS( 405 )
  W( 87 ) = W( 87 ) + a*JVS( 406 )
  W( 92 ) = W( 92 ) + a*JVS( 407 )
  W( 99 ) = W( 99 ) + a*JVS( 408 )
  W( 101 ) = W( 101 ) + a*JVS( 409 )
  W( 103 ) = W( 103 ) + a*JVS( 410 )
  a = -W( 76 ) / JVS( 433 )
  W( 76 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 434 )
  W( 92 ) = W( 92 ) + a*JVS( 435 )
  W( 99 ) = W( 99 ) + a*JVS( 436 )
  W( 103 ) = W( 103 ) + a*JVS( 437 )
  a = -W( 78 ) / JVS( 469 )
  W( 78 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 470 )
  W( 88 ) = W( 88 ) + a*JVS( 471 )
  W( 92 ) = W( 92 ) + a*JVS( 472 )
  W( 99 ) = W( 99 ) + a*JVS( 473 )
  W( 103 ) = W( 103 ) + a*JVS( 474 )
  a = -W( 79 ) / JVS( 475 )
  W( 79 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 92 ) = W( 92 ) + a*JVS( 477 )
  W( 99 ) = W( 99 ) + a*JVS( 478 )
  W( 103 ) = W( 103 ) + a*JVS( 479 )
  a = -W( 80 ) / JVS( 480 )
  W( 80 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 92 ) = W( 92 ) + a*JVS( 482 )
  W( 99 ) = W( 99 ) + a*JVS( 483 )
  W( 103 ) = W( 103 ) + a*JVS( 484 )
  a = -W( 82 ) / JVS( 510 )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 511 )
  W( 88 ) = W( 88 ) + a*JVS( 512 )
  W( 92 ) = W( 92 ) + a*JVS( 513 )
  W( 99 ) = W( 99 ) + a*JVS( 514 )
  W( 103 ) = W( 103 ) + a*JVS( 515 )
  a = -W( 83 ) / JVS( 525 )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 526 )
  W( 87 ) = W( 87 ) + a*JVS( 527 )
  W( 88 ) = W( 88 ) + a*JVS( 528 )
  W( 89 ) = W( 89 ) + a*JVS( 529 )
  W( 90 ) = W( 90 ) + a*JVS( 530 )
  W( 91 ) = W( 91 ) + a*JVS( 531 )
  W( 92 ) = W( 92 ) + a*JVS( 532 )
  W( 93 ) = W( 93 ) + a*JVS( 533 )
  W( 95 ) = W( 95 ) + a*JVS( 534 )
  W( 97 ) = W( 97 ) + a*JVS( 535 )
  W( 98 ) = W( 98 ) + a*JVS( 536 )
  W( 99 ) = W( 99 ) + a*JVS( 537 )
  W( 100 ) = W( 100 ) + a*JVS( 538 )
  W( 103 ) = W( 103 ) + a*JVS( 539 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  a = -W( 85 ) / JVS( 547 )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 548 )
  W( 92 ) = W( 92 ) + a*JVS( 549 )
  W( 99 ) = W( 99 ) + a*JVS( 550 )
  W( 103 ) = W( 103 ) + a*JVS( 551 )
  a = -W( 87 ) / JVS( 596 )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 597 )
  W( 92 ) = W( 92 ) + a*JVS( 598 )
  W( 93 ) = W( 93 ) + a*JVS( 599 )
  W( 99 ) = W( 99 ) + a*JVS( 600 )
  W( 101 ) = W( 101 ) + a*JVS( 601 )
  W( 103 ) = W( 103 ) + a*JVS( 602 )
  W( 104 ) = W( 104 ) + a*JVS( 603 )
  a = -W( 88 ) / JVS( 614 )
  W( 88 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 615 )
  W( 93 ) = W( 93 ) + a*JVS( 616 )
  W( 99 ) = W( 99 ) + a*JVS( 617 )
  W( 101 ) = W( 101 ) + a*JVS( 618 )
  W( 103 ) = W( 103 ) + a*JVS( 619 )
  a = -W( 89 ) / JVS( 639 )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 640 )
  W( 91 ) = W( 91 ) + a*JVS( 641 )
  W( 92 ) = W( 92 ) + a*JVS( 642 )
  W( 93 ) = W( 93 ) + a*JVS( 643 )
  W( 94 ) = W( 94 ) + a*JVS( 644 )
  W( 96 ) = W( 96 ) + a*JVS( 645 )
  W( 99 ) = W( 99 ) + a*JVS( 646 )
  W( 101 ) = W( 101 ) + a*JVS( 647 )
  W( 103 ) = W( 103 ) + a*JVS( 648 )
  W( 104 ) = W( 104 ) + a*JVS( 649 )
  a = -W( 90 ) / JVS( 661 )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 662 )
  W( 92 ) = W( 92 ) + a*JVS( 663 )
  W( 93 ) = W( 93 ) + a*JVS( 664 )
  W( 94 ) = W( 94 ) + a*JVS( 665 )
  W( 99 ) = W( 99 ) + a*JVS( 666 )
  W( 101 ) = W( 101 ) + a*JVS( 667 )
  W( 102 ) = W( 102 ) + a*JVS( 668 )
  W( 103 ) = W( 103 ) + a*JVS( 669 )
  W( 104 ) = W( 104 ) + a*JVS( 670 )
  a = -W( 91 ) / JVS( 683 )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 684 )
  W( 93 ) = W( 93 ) + a*JVS( 685 )
  W( 94 ) = W( 94 ) + a*JVS( 686 )
  W( 95 ) = W( 95 ) + a*JVS( 687 )
  W( 97 ) = W( 97 ) + a*JVS( 688 )
  W( 99 ) = W( 99 ) + a*JVS( 689 )
  W( 100 ) = W( 100 ) + a*JVS( 690 )
  W( 101 ) = W( 101 ) + a*JVS( 691 )
  W( 102 ) = W( 102 ) + a*JVS( 692 )
  W( 103 ) = W( 103 ) + a*JVS( 693 )
  W( 104 ) = W( 104 ) + a*JVS( 694 )
  a = -W( 92 ) / JVS( 708 )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 709 )
  W( 95 ) = W( 95 ) + a*JVS( 710 )
  W( 96 ) = W( 96 ) + a*JVS( 711 )
  W( 97 ) = W( 97 ) + a*JVS( 712 )
  W( 98 ) = W( 98 ) + a*JVS( 713 )
  W( 99 ) = W( 99 ) + a*JVS( 714 )
  W( 100 ) = W( 100 ) + a*JVS( 715 )
  W( 101 ) = W( 101 ) + a*JVS( 716 )
  W( 103 ) = W( 103 ) + a*JVS( 717 )
  a = -W( 93 ) / JVS( 731 )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 732 )
  W( 95 ) = W( 95 ) + a*JVS( 733 )
  W( 96 ) = W( 96 ) + a*JVS( 734 )
  W( 97 ) = W( 97 ) + a*JVS( 735 )
  W( 98 ) = W( 98 ) + a*JVS( 736 )
  W( 99 ) = W( 99 ) + a*JVS( 737 )
  W( 100 ) = W( 100 ) + a*JVS( 738 )
  W( 101 ) = W( 101 ) + a*JVS( 739 )
  W( 102 ) = W( 102 ) + a*JVS( 740 )
  W( 103 ) = W( 103 ) + a*JVS( 741 )
  W( 104 ) = W( 104 ) + a*JVS( 742 )
  a = -W( 94 ) / JVS( 776 )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 777 )
  W( 96 ) = W( 96 ) + a*JVS( 778 )
  W( 97 ) = W( 97 ) + a*JVS( 779 )
  W( 98 ) = W( 98 ) + a*JVS( 780 )
  W( 99 ) = W( 99 ) + a*JVS( 781 )
  W( 100 ) = W( 100 ) + a*JVS( 782 )
  W( 101 ) = W( 101 ) + a*JVS( 783 )
  W( 102 ) = W( 102 ) + a*JVS( 784 )
  W( 103 ) = W( 103 ) + a*JVS( 785 )
  W( 104 ) = W( 104 ) + a*JVS( 786 )
  JVS( 787) = W( 38 )
  JVS( 788) = W( 45 )
  JVS( 789) = W( 51 )
  JVS( 790) = W( 58 )
  JVS( 791) = W( 59 )
  JVS( 792) = W( 61 )
  JVS( 793) = W( 68 )
  JVS( 794) = W( 72 )
  JVS( 795) = W( 76 )
  JVS( 796) = W( 78 )
  JVS( 797) = W( 79 )
  JVS( 798) = W( 80 )
  JVS( 799) = W( 82 )
  JVS( 800) = W( 83 )
  JVS( 801) = W( 84 )
  JVS( 802) = W( 85 )
  JVS( 803) = W( 87 )
  JVS( 804) = W( 88 )
  JVS( 805) = W( 89 )
  JVS( 806) = W( 90 )
  JVS( 807) = W( 91 )
  JVS( 808) = W( 92 )
  JVS( 809) = W( 93 )
  JVS( 810) = W( 94 )
  JVS( 811) = W( 95 )
  JVS( 812) = W( 96 )
  JVS( 813) = W( 97 )
  JVS( 814) = W( 98 )
  JVS( 815) = W( 99 )
  JVS( 816) = W( 100 )
  JVS( 817) = W( 101 )
  JVS( 818) = W( 102 )
  JVS( 819) = W( 103 )
  JVS( 820) = W( 104 )
  IF ( ABS( JVS( 863 )) < TINY(a) ) THEN
         IER = 96
         RETURN
  END IF
   W( 36 ) = JVS( 821 )
   W( 42 ) = JVS( 822 )
   W( 44 ) = JVS( 823 )
   W( 47 ) = JVS( 824 )
   W( 51 ) = JVS( 825 )
   W( 52 ) = JVS( 826 )
   W( 53 ) = JVS( 827 )
   W( 54 ) = JVS( 828 )
   W( 55 ) = JVS( 829 )
   W( 56 ) = JVS( 830 )
   W( 57 ) = JVS( 831 )
   W( 58 ) = JVS( 832 )
   W( 59 ) = JVS( 833 )
   W( 62 ) = JVS( 834 )
   W( 63 ) = JVS( 835 )
   W( 64 ) = JVS( 836 )
   W( 65 ) = JVS( 837 )
   W( 68 ) = JVS( 838 )
   W( 69 ) = JVS( 839 )
   W( 71 ) = JVS( 840 )
   W( 73 ) = JVS( 841 )
   W( 74 ) = JVS( 842 )
   W( 75 ) = JVS( 843 )
   W( 76 ) = JVS( 844 )
   W( 78 ) = JVS( 845 )
   W( 79 ) = JVS( 846 )
   W( 80 ) = JVS( 847 )
   W( 81 ) = JVS( 848 )
   W( 82 ) = JVS( 849 )
   W( 83 ) = JVS( 850 )
   W( 84 ) = JVS( 851 )
   W( 85 ) = JVS( 852 )
   W( 86 ) = JVS( 853 )
   W( 87 ) = JVS( 854 )
   W( 88 ) = JVS( 855 )
   W( 89 ) = JVS( 856 )
   W( 90 ) = JVS( 857 )
   W( 91 ) = JVS( 858 )
   W( 92 ) = JVS( 859 )
   W( 93 ) = JVS( 860 )
   W( 94 ) = JVS( 861 )
   W( 95 ) = JVS( 862 )
   W( 96 ) = JVS( 863 )
   W( 97 ) = JVS( 864 )
   W( 98 ) = JVS( 865 )
   W( 99 ) = JVS( 866 )
   W( 100 ) = JVS( 867 )
   W( 101 ) = JVS( 868 )
   W( 102 ) = JVS( 869 )
   W( 103 ) = JVS( 870 )
   W( 104 ) = JVS( 871 )
  a = -W( 36 ) / JVS( 217 )
  W( 36 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 218 )
  a = -W( 42 ) / JVS( 233 )
  W( 42 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 234 )
  W( 103 ) = W( 103 ) + a*JVS( 235 )
  a = -W( 44 ) / JVS( 238 )
  W( 44 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 239 )
  a = -W( 47 ) / JVS( 249 )
  W( 47 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 250 )
  W( 103 ) = W( 103 ) + a*JVS( 251 )
  a = -W( 51 ) / JVS( 260 )
  W( 51 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 261 )
  a = -W( 52 ) / JVS( 262 )
  W( 52 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 263 )
  W( 102 ) = W( 102 ) + a*JVS( 264 )
  W( 103 ) = W( 103 ) + a*JVS( 265 )
  a = -W( 53 ) / JVS( 266 )
  W( 53 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 267 )
  W( 93 ) = W( 93 ) + a*JVS( 268 )
  W( 96 ) = W( 96 ) + a*JVS( 269 )
  a = -W( 54 ) / JVS( 270 )
  W( 54 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 271 )
  W( 96 ) = W( 96 ) + a*JVS( 272 )
  W( 99 ) = W( 99 ) + a*JVS( 273 )
  W( 101 ) = W( 101 ) + a*JVS( 274 )
  a = -W( 55 ) / JVS( 275 )
  W( 55 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 276 )
  W( 101 ) = W( 101 ) + a*JVS( 277 )
  W( 103 ) = W( 103 ) + a*JVS( 278 )
  a = -W( 56 ) / JVS( 279 )
  W( 56 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 280 )
  W( 102 ) = W( 102 ) + a*JVS( 281 )
  W( 103 ) = W( 103 ) + a*JVS( 282 )
  W( 104 ) = W( 104 ) + a*JVS( 283 )
  a = -W( 57 ) / JVS( 284 )
  W( 57 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 285 )
  a = -W( 58 ) / JVS( 288 )
  W( 58 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 289 )
  a = -W( 59 ) / JVS( 292 )
  W( 59 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 293 )
  a = -W( 62 ) / JVS( 303 )
  W( 62 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 304 )
  W( 103 ) = W( 103 ) + a*JVS( 305 )
  a = -W( 63 ) / JVS( 306 )
  W( 63 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 307 )
  W( 103 ) = W( 103 ) + a*JVS( 308 )
  a = -W( 64 ) / JVS( 310 )
  W( 64 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 311 )
  W( 96 ) = W( 96 ) + a*JVS( 312 )
  W( 99 ) = W( 99 ) + a*JVS( 313 )
  W( 101 ) = W( 101 ) + a*JVS( 314 )
  a = -W( 65 ) / JVS( 315 )
  W( 65 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 316 )
  W( 96 ) = W( 96 ) + a*JVS( 317 )
  W( 103 ) = W( 103 ) + a*JVS( 318 )
  W( 104 ) = W( 104 ) + a*JVS( 319 )
  a = -W( 68 ) / JVS( 338 )
  W( 68 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 339 )
  W( 82 ) = W( 82 ) + a*JVS( 340 )
  W( 85 ) = W( 85 ) + a*JVS( 341 )
  W( 92 ) = W( 92 ) + a*JVS( 342 )
  W( 99 ) = W( 99 ) + a*JVS( 343 )
  W( 103 ) = W( 103 ) + a*JVS( 344 )
  a = -W( 69 ) / JVS( 351 )
  W( 69 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 352 )
  W( 73 ) = W( 73 ) + a*JVS( 353 )
  W( 75 ) = W( 75 ) + a*JVS( 354 )
  W( 76 ) = W( 76 ) + a*JVS( 355 )
  W( 78 ) = W( 78 ) + a*JVS( 356 )
  W( 79 ) = W( 79 ) + a*JVS( 357 )
  W( 80 ) = W( 80 ) + a*JVS( 358 )
  W( 81 ) = W( 81 ) + a*JVS( 359 )
  W( 82 ) = W( 82 ) + a*JVS( 360 )
  W( 83 ) = W( 83 ) + a*JVS( 361 )
  W( 84 ) = W( 84 ) + a*JVS( 362 )
  W( 85 ) = W( 85 ) + a*JVS( 363 )
  W( 86 ) = W( 86 ) + a*JVS( 364 )
  W( 88 ) = W( 88 ) + a*JVS( 365 )
  W( 89 ) = W( 89 ) + a*JVS( 366 )
  W( 92 ) = W( 92 ) + a*JVS( 367 )
  W( 99 ) = W( 99 ) + a*JVS( 368 )
  W( 103 ) = W( 103 ) + a*JVS( 369 )
  a = -W( 71 ) / JVS( 391 )
  W( 71 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 392 )
  W( 92 ) = W( 92 ) + a*JVS( 393 )
  W( 99 ) = W( 99 ) + a*JVS( 394 )
  W( 103 ) = W( 103 ) + a*JVS( 395 )
  a = -W( 73 ) / JVS( 411 )
  W( 73 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 92 ) = W( 92 ) + a*JVS( 413 )
  W( 99 ) = W( 99 ) + a*JVS( 414 )
  W( 103 ) = W( 103 ) + a*JVS( 415 )
  a = -W( 74 ) / JVS( 418 )
  W( 74 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 419 )
  W( 95 ) = W( 95 ) + a*JVS( 420 )
  W( 96 ) = W( 96 ) + a*JVS( 421 )
  W( 97 ) = W( 97 ) + a*JVS( 422 )
  W( 98 ) = W( 98 ) + a*JVS( 423 )
  W( 99 ) = W( 99 ) + a*JVS( 424 )
  W( 100 ) = W( 100 ) + a*JVS( 425 )
  W( 101 ) = W( 101 ) + a*JVS( 426 )
  W( 103 ) = W( 103 ) + a*JVS( 427 )
  a = -W( 75 ) / JVS( 428 )
  W( 75 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 429 )
  W( 92 ) = W( 92 ) + a*JVS( 430 )
  W( 99 ) = W( 99 ) + a*JVS( 431 )
  W( 103 ) = W( 103 ) + a*JVS( 432 )
  a = -W( 76 ) / JVS( 433 )
  W( 76 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 434 )
  W( 92 ) = W( 92 ) + a*JVS( 435 )
  W( 99 ) = W( 99 ) + a*JVS( 436 )
  W( 103 ) = W( 103 ) + a*JVS( 437 )
  a = -W( 78 ) / JVS( 469 )
  W( 78 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 470 )
  W( 88 ) = W( 88 ) + a*JVS( 471 )
  W( 92 ) = W( 92 ) + a*JVS( 472 )
  W( 99 ) = W( 99 ) + a*JVS( 473 )
  W( 103 ) = W( 103 ) + a*JVS( 474 )
  a = -W( 79 ) / JVS( 475 )
  W( 79 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 92 ) = W( 92 ) + a*JVS( 477 )
  W( 99 ) = W( 99 ) + a*JVS( 478 )
  W( 103 ) = W( 103 ) + a*JVS( 479 )
  a = -W( 80 ) / JVS( 480 )
  W( 80 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 92 ) = W( 92 ) + a*JVS( 482 )
  W( 99 ) = W( 99 ) + a*JVS( 483 )
  W( 103 ) = W( 103 ) + a*JVS( 484 )
  a = -W( 81 ) / JVS( 496 )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 497 )
  W( 88 ) = W( 88 ) + a*JVS( 498 )
  W( 92 ) = W( 92 ) + a*JVS( 499 )
  W( 93 ) = W( 93 ) + a*JVS( 500 )
  W( 95 ) = W( 95 ) + a*JVS( 501 )
  W( 96 ) = W( 96 ) + a*JVS( 502 )
  W( 97 ) = W( 97 ) + a*JVS( 503 )
  W( 98 ) = W( 98 ) + a*JVS( 504 )
  W( 99 ) = W( 99 ) + a*JVS( 505 )
  W( 100 ) = W( 100 ) + a*JVS( 506 )
  W( 101 ) = W( 101 ) + a*JVS( 507 )
  W( 103 ) = W( 103 ) + a*JVS( 508 )
  a = -W( 82 ) / JVS( 510 )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 511 )
  W( 88 ) = W( 88 ) + a*JVS( 512 )
  W( 92 ) = W( 92 ) + a*JVS( 513 )
  W( 99 ) = W( 99 ) + a*JVS( 514 )
  W( 103 ) = W( 103 ) + a*JVS( 515 )
  a = -W( 83 ) / JVS( 525 )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 526 )
  W( 87 ) = W( 87 ) + a*JVS( 527 )
  W( 88 ) = W( 88 ) + a*JVS( 528 )
  W( 89 ) = W( 89 ) + a*JVS( 529 )
  W( 90 ) = W( 90 ) + a*JVS( 530 )
  W( 91 ) = W( 91 ) + a*JVS( 531 )
  W( 92 ) = W( 92 ) + a*JVS( 532 )
  W( 93 ) = W( 93 ) + a*JVS( 533 )
  W( 95 ) = W( 95 ) + a*JVS( 534 )
  W( 97 ) = W( 97 ) + a*JVS( 535 )
  W( 98 ) = W( 98 ) + a*JVS( 536 )
  W( 99 ) = W( 99 ) + a*JVS( 537 )
  W( 100 ) = W( 100 ) + a*JVS( 538 )
  W( 103 ) = W( 103 ) + a*JVS( 539 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  a = -W( 85 ) / JVS( 547 )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 548 )
  W( 92 ) = W( 92 ) + a*JVS( 549 )
  W( 99 ) = W( 99 ) + a*JVS( 550 )
  W( 103 ) = W( 103 ) + a*JVS( 551 )
  a = -W( 86 ) / JVS( 572 )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 573 )
  W( 88 ) = W( 88 ) + a*JVS( 574 )
  W( 90 ) = W( 90 ) + a*JVS( 575 )
  W( 91 ) = W( 91 ) + a*JVS( 576 )
  W( 92 ) = W( 92 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 94 ) = W( 94 ) + a*JVS( 579 )
  W( 95 ) = W( 95 ) + a*JVS( 580 )
  W( 96 ) = W( 96 ) + a*JVS( 581 )
  W( 97 ) = W( 97 ) + a*JVS( 582 )
  W( 98 ) = W( 98 ) + a*JVS( 583 )
  W( 99 ) = W( 99 ) + a*JVS( 584 )
  W( 100 ) = W( 100 ) + a*JVS( 585 )
  W( 101 ) = W( 101 ) + a*JVS( 586 )
  W( 102 ) = W( 102 ) + a*JVS( 587 )
  W( 103 ) = W( 103 ) + a*JVS( 588 )
  W( 104 ) = W( 104 ) + a*JVS( 589 )
  a = -W( 87 ) / JVS( 596 )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 597 )
  W( 92 ) = W( 92 ) + a*JVS( 598 )
  W( 93 ) = W( 93 ) + a*JVS( 599 )
  W( 99 ) = W( 99 ) + a*JVS( 600 )
  W( 101 ) = W( 101 ) + a*JVS( 601 )
  W( 103 ) = W( 103 ) + a*JVS( 602 )
  W( 104 ) = W( 104 ) + a*JVS( 603 )
  a = -W( 88 ) / JVS( 614 )
  W( 88 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 615 )
  W( 93 ) = W( 93 ) + a*JVS( 616 )
  W( 99 ) = W( 99 ) + a*JVS( 617 )
  W( 101 ) = W( 101 ) + a*JVS( 618 )
  W( 103 ) = W( 103 ) + a*JVS( 619 )
  a = -W( 89 ) / JVS( 639 )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 640 )
  W( 91 ) = W( 91 ) + a*JVS( 641 )
  W( 92 ) = W( 92 ) + a*JVS( 642 )
  W( 93 ) = W( 93 ) + a*JVS( 643 )
  W( 94 ) = W( 94 ) + a*JVS( 644 )
  W( 96 ) = W( 96 ) + a*JVS( 645 )
  W( 99 ) = W( 99 ) + a*JVS( 646 )
  W( 101 ) = W( 101 ) + a*JVS( 647 )
  W( 103 ) = W( 103 ) + a*JVS( 648 )
  W( 104 ) = W( 104 ) + a*JVS( 649 )
  a = -W( 90 ) / JVS( 661 )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 662 )
  W( 92 ) = W( 92 ) + a*JVS( 663 )
  W( 93 ) = W( 93 ) + a*JVS( 664 )
  W( 94 ) = W( 94 ) + a*JVS( 665 )
  W( 99 ) = W( 99 ) + a*JVS( 666 )
  W( 101 ) = W( 101 ) + a*JVS( 667 )
  W( 102 ) = W( 102 ) + a*JVS( 668 )
  W( 103 ) = W( 103 ) + a*JVS( 669 )
  W( 104 ) = W( 104 ) + a*JVS( 670 )
  a = -W( 91 ) / JVS( 683 )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 684 )
  W( 93 ) = W( 93 ) + a*JVS( 685 )
  W( 94 ) = W( 94 ) + a*JVS( 686 )
  W( 95 ) = W( 95 ) + a*JVS( 687 )
  W( 97 ) = W( 97 ) + a*JVS( 688 )
  W( 99 ) = W( 99 ) + a*JVS( 689 )
  W( 100 ) = W( 100 ) + a*JVS( 690 )
  W( 101 ) = W( 101 ) + a*JVS( 691 )
  W( 102 ) = W( 102 ) + a*JVS( 692 )
  W( 103 ) = W( 103 ) + a*JVS( 693 )
  W( 104 ) = W( 104 ) + a*JVS( 694 )
  a = -W( 92 ) / JVS( 708 )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 709 )
  W( 95 ) = W( 95 ) + a*JVS( 710 )
  W( 96 ) = W( 96 ) + a*JVS( 711 )
  W( 97 ) = W( 97 ) + a*JVS( 712 )
  W( 98 ) = W( 98 ) + a*JVS( 713 )
  W( 99 ) = W( 99 ) + a*JVS( 714 )
  W( 100 ) = W( 100 ) + a*JVS( 715 )
  W( 101 ) = W( 101 ) + a*JVS( 716 )
  W( 103 ) = W( 103 ) + a*JVS( 717 )
  a = -W( 93 ) / JVS( 731 )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 732 )
  W( 95 ) = W( 95 ) + a*JVS( 733 )
  W( 96 ) = W( 96 ) + a*JVS( 734 )
  W( 97 ) = W( 97 ) + a*JVS( 735 )
  W( 98 ) = W( 98 ) + a*JVS( 736 )
  W( 99 ) = W( 99 ) + a*JVS( 737 )
  W( 100 ) = W( 100 ) + a*JVS( 738 )
  W( 101 ) = W( 101 ) + a*JVS( 739 )
  W( 102 ) = W( 102 ) + a*JVS( 740 )
  W( 103 ) = W( 103 ) + a*JVS( 741 )
  W( 104 ) = W( 104 ) + a*JVS( 742 )
  a = -W( 94 ) / JVS( 776 )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 777 )
  W( 96 ) = W( 96 ) + a*JVS( 778 )
  W( 97 ) = W( 97 ) + a*JVS( 779 )
  W( 98 ) = W( 98 ) + a*JVS( 780 )
  W( 99 ) = W( 99 ) + a*JVS( 781 )
  W( 100 ) = W( 100 ) + a*JVS( 782 )
  W( 101 ) = W( 101 ) + a*JVS( 783 )
  W( 102 ) = W( 102 ) + a*JVS( 784 )
  W( 103 ) = W( 103 ) + a*JVS( 785 )
  W( 104 ) = W( 104 ) + a*JVS( 786 )
  a = -W( 95 ) / JVS( 811 )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 812 )
  W( 97 ) = W( 97 ) + a*JVS( 813 )
  W( 98 ) = W( 98 ) + a*JVS( 814 )
  W( 99 ) = W( 99 ) + a*JVS( 815 )
  W( 100 ) = W( 100 ) + a*JVS( 816 )
  W( 101 ) = W( 101 ) + a*JVS( 817 )
  W( 102 ) = W( 102 ) + a*JVS( 818 )
  W( 103 ) = W( 103 ) + a*JVS( 819 )
  W( 104 ) = W( 104 ) + a*JVS( 820 )
  JVS( 821) = W( 36 )
  JVS( 822) = W( 42 )
  JVS( 823) = W( 44 )
  JVS( 824) = W( 47 )
  JVS( 825) = W( 51 )
  JVS( 826) = W( 52 )
  JVS( 827) = W( 53 )
  JVS( 828) = W( 54 )
  JVS( 829) = W( 55 )
  JVS( 830) = W( 56 )
  JVS( 831) = W( 57 )
  JVS( 832) = W( 58 )
  JVS( 833) = W( 59 )
  JVS( 834) = W( 62 )
  JVS( 835) = W( 63 )
  JVS( 836) = W( 64 )
  JVS( 837) = W( 65 )
  JVS( 838) = W( 68 )
  JVS( 839) = W( 69 )
  JVS( 840) = W( 71 )
  JVS( 841) = W( 73 )
  JVS( 842) = W( 74 )
  JVS( 843) = W( 75 )
  JVS( 844) = W( 76 )
  JVS( 845) = W( 78 )
  JVS( 846) = W( 79 )
  JVS( 847) = W( 80 )
  JVS( 848) = W( 81 )
  JVS( 849) = W( 82 )
  JVS( 850) = W( 83 )
  JVS( 851) = W( 84 )
  JVS( 852) = W( 85 )
  JVS( 853) = W( 86 )
  JVS( 854) = W( 87 )
  JVS( 855) = W( 88 )
  JVS( 856) = W( 89 )
  JVS( 857) = W( 90 )
  JVS( 858) = W( 91 )
  JVS( 859) = W( 92 )
  JVS( 860) = W( 93 )
  JVS( 861) = W( 94 )
  JVS( 862) = W( 95 )
  JVS( 863) = W( 96 )
  JVS( 864) = W( 97 )
  JVS( 865) = W( 98 )
  JVS( 866) = W( 99 )
  JVS( 867) = W( 100 )
  JVS( 868) = W( 101 )
  JVS( 869) = W( 102 )
  JVS( 870) = W( 103 )
  JVS( 871) = W( 104 )
  IF ( ABS( JVS( 889 )) < TINY(a) ) THEN
         IER = 97
         RETURN
  END IF
   W( 39 ) = JVS( 872 )
   W( 78 ) = JVS( 873 )
   W( 79 ) = JVS( 874 )
   W( 80 ) = JVS( 875 )
   W( 81 ) = JVS( 876 )
   W( 82 ) = JVS( 877 )
   W( 84 ) = JVS( 878 )
   W( 85 ) = JVS( 879 )
   W( 88 ) = JVS( 880 )
   W( 89 ) = JVS( 881 )
   W( 90 ) = JVS( 882 )
   W( 91 ) = JVS( 883 )
   W( 92 ) = JVS( 884 )
   W( 93 ) = JVS( 885 )
   W( 94 ) = JVS( 886 )
   W( 95 ) = JVS( 887 )
   W( 96 ) = JVS( 888 )
   W( 97 ) = JVS( 889 )
   W( 98 ) = JVS( 890 )
   W( 99 ) = JVS( 891 )
   W( 100 ) = JVS( 892 )
   W( 101 ) = JVS( 893 )
   W( 102 ) = JVS( 894 )
   W( 103 ) = JVS( 895 )
   W( 104 ) = JVS( 896 )
  a = -W( 39 ) / JVS( 224 )
  W( 39 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 225 )
  W( 101 ) = W( 101 ) + a*JVS( 226 )
  a = -W( 78 ) / JVS( 469 )
  W( 78 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 470 )
  W( 88 ) = W( 88 ) + a*JVS( 471 )
  W( 92 ) = W( 92 ) + a*JVS( 472 )
  W( 99 ) = W( 99 ) + a*JVS( 473 )
  W( 103 ) = W( 103 ) + a*JVS( 474 )
  a = -W( 79 ) / JVS( 475 )
  W( 79 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 92 ) = W( 92 ) + a*JVS( 477 )
  W( 99 ) = W( 99 ) + a*JVS( 478 )
  W( 103 ) = W( 103 ) + a*JVS( 479 )
  a = -W( 80 ) / JVS( 480 )
  W( 80 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 92 ) = W( 92 ) + a*JVS( 482 )
  W( 99 ) = W( 99 ) + a*JVS( 483 )
  W( 103 ) = W( 103 ) + a*JVS( 484 )
  a = -W( 81 ) / JVS( 496 )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 497 )
  W( 88 ) = W( 88 ) + a*JVS( 498 )
  W( 92 ) = W( 92 ) + a*JVS( 499 )
  W( 93 ) = W( 93 ) + a*JVS( 500 )
  W( 95 ) = W( 95 ) + a*JVS( 501 )
  W( 96 ) = W( 96 ) + a*JVS( 502 )
  W( 97 ) = W( 97 ) + a*JVS( 503 )
  W( 98 ) = W( 98 ) + a*JVS( 504 )
  W( 99 ) = W( 99 ) + a*JVS( 505 )
  W( 100 ) = W( 100 ) + a*JVS( 506 )
  W( 101 ) = W( 101 ) + a*JVS( 507 )
  W( 103 ) = W( 103 ) + a*JVS( 508 )
  a = -W( 82 ) / JVS( 510 )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 511 )
  W( 88 ) = W( 88 ) + a*JVS( 512 )
  W( 92 ) = W( 92 ) + a*JVS( 513 )
  W( 99 ) = W( 99 ) + a*JVS( 514 )
  W( 103 ) = W( 103 ) + a*JVS( 515 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  a = -W( 85 ) / JVS( 547 )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 548 )
  W( 92 ) = W( 92 ) + a*JVS( 549 )
  W( 99 ) = W( 99 ) + a*JVS( 550 )
  W( 103 ) = W( 103 ) + a*JVS( 551 )
  a = -W( 88 ) / JVS( 614 )
  W( 88 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 615 )
  W( 93 ) = W( 93 ) + a*JVS( 616 )
  W( 99 ) = W( 99 ) + a*JVS( 617 )
  W( 101 ) = W( 101 ) + a*JVS( 618 )
  W( 103 ) = W( 103 ) + a*JVS( 619 )
  a = -W( 89 ) / JVS( 639 )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 640 )
  W( 91 ) = W( 91 ) + a*JVS( 641 )
  W( 92 ) = W( 92 ) + a*JVS( 642 )
  W( 93 ) = W( 93 ) + a*JVS( 643 )
  W( 94 ) = W( 94 ) + a*JVS( 644 )
  W( 96 ) = W( 96 ) + a*JVS( 645 )
  W( 99 ) = W( 99 ) + a*JVS( 646 )
  W( 101 ) = W( 101 ) + a*JVS( 647 )
  W( 103 ) = W( 103 ) + a*JVS( 648 )
  W( 104 ) = W( 104 ) + a*JVS( 649 )
  a = -W( 90 ) / JVS( 661 )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 662 )
  W( 92 ) = W( 92 ) + a*JVS( 663 )
  W( 93 ) = W( 93 ) + a*JVS( 664 )
  W( 94 ) = W( 94 ) + a*JVS( 665 )
  W( 99 ) = W( 99 ) + a*JVS( 666 )
  W( 101 ) = W( 101 ) + a*JVS( 667 )
  W( 102 ) = W( 102 ) + a*JVS( 668 )
  W( 103 ) = W( 103 ) + a*JVS( 669 )
  W( 104 ) = W( 104 ) + a*JVS( 670 )
  a = -W( 91 ) / JVS( 683 )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 684 )
  W( 93 ) = W( 93 ) + a*JVS( 685 )
  W( 94 ) = W( 94 ) + a*JVS( 686 )
  W( 95 ) = W( 95 ) + a*JVS( 687 )
  W( 97 ) = W( 97 ) + a*JVS( 688 )
  W( 99 ) = W( 99 ) + a*JVS( 689 )
  W( 100 ) = W( 100 ) + a*JVS( 690 )
  W( 101 ) = W( 101 ) + a*JVS( 691 )
  W( 102 ) = W( 102 ) + a*JVS( 692 )
  W( 103 ) = W( 103 ) + a*JVS( 693 )
  W( 104 ) = W( 104 ) + a*JVS( 694 )
  a = -W( 92 ) / JVS( 708 )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 709 )
  W( 95 ) = W( 95 ) + a*JVS( 710 )
  W( 96 ) = W( 96 ) + a*JVS( 711 )
  W( 97 ) = W( 97 ) + a*JVS( 712 )
  W( 98 ) = W( 98 ) + a*JVS( 713 )
  W( 99 ) = W( 99 ) + a*JVS( 714 )
  W( 100 ) = W( 100 ) + a*JVS( 715 )
  W( 101 ) = W( 101 ) + a*JVS( 716 )
  W( 103 ) = W( 103 ) + a*JVS( 717 )
  a = -W( 93 ) / JVS( 731 )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 732 )
  W( 95 ) = W( 95 ) + a*JVS( 733 )
  W( 96 ) = W( 96 ) + a*JVS( 734 )
  W( 97 ) = W( 97 ) + a*JVS( 735 )
  W( 98 ) = W( 98 ) + a*JVS( 736 )
  W( 99 ) = W( 99 ) + a*JVS( 737 )
  W( 100 ) = W( 100 ) + a*JVS( 738 )
  W( 101 ) = W( 101 ) + a*JVS( 739 )
  W( 102 ) = W( 102 ) + a*JVS( 740 )
  W( 103 ) = W( 103 ) + a*JVS( 741 )
  W( 104 ) = W( 104 ) + a*JVS( 742 )
  a = -W( 94 ) / JVS( 776 )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 777 )
  W( 96 ) = W( 96 ) + a*JVS( 778 )
  W( 97 ) = W( 97 ) + a*JVS( 779 )
  W( 98 ) = W( 98 ) + a*JVS( 780 )
  W( 99 ) = W( 99 ) + a*JVS( 781 )
  W( 100 ) = W( 100 ) + a*JVS( 782 )
  W( 101 ) = W( 101 ) + a*JVS( 783 )
  W( 102 ) = W( 102 ) + a*JVS( 784 )
  W( 103 ) = W( 103 ) + a*JVS( 785 )
  W( 104 ) = W( 104 ) + a*JVS( 786 )
  a = -W( 95 ) / JVS( 811 )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 812 )
  W( 97 ) = W( 97 ) + a*JVS( 813 )
  W( 98 ) = W( 98 ) + a*JVS( 814 )
  W( 99 ) = W( 99 ) + a*JVS( 815 )
  W( 100 ) = W( 100 ) + a*JVS( 816 )
  W( 101 ) = W( 101 ) + a*JVS( 817 )
  W( 102 ) = W( 102 ) + a*JVS( 818 )
  W( 103 ) = W( 103 ) + a*JVS( 819 )
  W( 104 ) = W( 104 ) + a*JVS( 820 )
  a = -W( 96 ) / JVS( 863 )
  W( 96 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 864 )
  W( 98 ) = W( 98 ) + a*JVS( 865 )
  W( 99 ) = W( 99 ) + a*JVS( 866 )
  W( 100 ) = W( 100 ) + a*JVS( 867 )
  W( 101 ) = W( 101 ) + a*JVS( 868 )
  W( 102 ) = W( 102 ) + a*JVS( 869 )
  W( 103 ) = W( 103 ) + a*JVS( 870 )
  W( 104 ) = W( 104 ) + a*JVS( 871 )
  JVS( 872) = W( 39 )
  JVS( 873) = W( 78 )
  JVS( 874) = W( 79 )
  JVS( 875) = W( 80 )
  JVS( 876) = W( 81 )
  JVS( 877) = W( 82 )
  JVS( 878) = W( 84 )
  JVS( 879) = W( 85 )
  JVS( 880) = W( 88 )
  JVS( 881) = W( 89 )
  JVS( 882) = W( 90 )
  JVS( 883) = W( 91 )
  JVS( 884) = W( 92 )
  JVS( 885) = W( 93 )
  JVS( 886) = W( 94 )
  JVS( 887) = W( 95 )
  JVS( 888) = W( 96 )
  JVS( 889) = W( 97 )
  JVS( 890) = W( 98 )
  JVS( 891) = W( 99 )
  JVS( 892) = W( 100 )
  JVS( 893) = W( 101 )
  JVS( 894) = W( 102 )
  JVS( 895) = W( 103 )
  JVS( 896) = W( 104 )
  IF ( ABS( JVS( 910 )) < TINY(a) ) THEN
         IER = 98
         RETURN
  END IF
   W( 41 ) = JVS( 897 )
   W( 75 ) = JVS( 898 )
   W( 78 ) = JVS( 899 )
   W( 82 ) = JVS( 900 )
   W( 84 ) = JVS( 901 )
   W( 85 ) = JVS( 902 )
   W( 88 ) = JVS( 903 )
   W( 92 ) = JVS( 904 )
   W( 93 ) = JVS( 905 )
   W( 94 ) = JVS( 906 )
   W( 95 ) = JVS( 907 )
   W( 96 ) = JVS( 908 )
   W( 97 ) = JVS( 909 )
   W( 98 ) = JVS( 910 )
   W( 99 ) = JVS( 911 )
   W( 100 ) = JVS( 912 )
   W( 101 ) = JVS( 913 )
   W( 102 ) = JVS( 914 )
   W( 103 ) = JVS( 915 )
   W( 104 ) = JVS( 916 )
  a = -W( 41 ) / JVS( 230 )
  W( 41 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 231 )
  W( 101 ) = W( 101 ) + a*JVS( 232 )
  a = -W( 75 ) / JVS( 428 )
  W( 75 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 429 )
  W( 92 ) = W( 92 ) + a*JVS( 430 )
  W( 99 ) = W( 99 ) + a*JVS( 431 )
  W( 103 ) = W( 103 ) + a*JVS( 432 )
  a = -W( 78 ) / JVS( 469 )
  W( 78 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 470 )
  W( 88 ) = W( 88 ) + a*JVS( 471 )
  W( 92 ) = W( 92 ) + a*JVS( 472 )
  W( 99 ) = W( 99 ) + a*JVS( 473 )
  W( 103 ) = W( 103 ) + a*JVS( 474 )
  a = -W( 82 ) / JVS( 510 )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 511 )
  W( 88 ) = W( 88 ) + a*JVS( 512 )
  W( 92 ) = W( 92 ) + a*JVS( 513 )
  W( 99 ) = W( 99 ) + a*JVS( 514 )
  W( 103 ) = W( 103 ) + a*JVS( 515 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  a = -W( 85 ) / JVS( 547 )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 548 )
  W( 92 ) = W( 92 ) + a*JVS( 549 )
  W( 99 ) = W( 99 ) + a*JVS( 550 )
  W( 103 ) = W( 103 ) + a*JVS( 551 )
  a = -W( 88 ) / JVS( 614 )
  W( 88 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 615 )
  W( 93 ) = W( 93 ) + a*JVS( 616 )
  W( 99 ) = W( 99 ) + a*JVS( 617 )
  W( 101 ) = W( 101 ) + a*JVS( 618 )
  W( 103 ) = W( 103 ) + a*JVS( 619 )
  a = -W( 92 ) / JVS( 708 )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 709 )
  W( 95 ) = W( 95 ) + a*JVS( 710 )
  W( 96 ) = W( 96 ) + a*JVS( 711 )
  W( 97 ) = W( 97 ) + a*JVS( 712 )
  W( 98 ) = W( 98 ) + a*JVS( 713 )
  W( 99 ) = W( 99 ) + a*JVS( 714 )
  W( 100 ) = W( 100 ) + a*JVS( 715 )
  W( 101 ) = W( 101 ) + a*JVS( 716 )
  W( 103 ) = W( 103 ) + a*JVS( 717 )
  a = -W( 93 ) / JVS( 731 )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 732 )
  W( 95 ) = W( 95 ) + a*JVS( 733 )
  W( 96 ) = W( 96 ) + a*JVS( 734 )
  W( 97 ) = W( 97 ) + a*JVS( 735 )
  W( 98 ) = W( 98 ) + a*JVS( 736 )
  W( 99 ) = W( 99 ) + a*JVS( 737 )
  W( 100 ) = W( 100 ) + a*JVS( 738 )
  W( 101 ) = W( 101 ) + a*JVS( 739 )
  W( 102 ) = W( 102 ) + a*JVS( 740 )
  W( 103 ) = W( 103 ) + a*JVS( 741 )
  W( 104 ) = W( 104 ) + a*JVS( 742 )
  a = -W( 94 ) / JVS( 776 )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 777 )
  W( 96 ) = W( 96 ) + a*JVS( 778 )
  W( 97 ) = W( 97 ) + a*JVS( 779 )
  W( 98 ) = W( 98 ) + a*JVS( 780 )
  W( 99 ) = W( 99 ) + a*JVS( 781 )
  W( 100 ) = W( 100 ) + a*JVS( 782 )
  W( 101 ) = W( 101 ) + a*JVS( 783 )
  W( 102 ) = W( 102 ) + a*JVS( 784 )
  W( 103 ) = W( 103 ) + a*JVS( 785 )
  W( 104 ) = W( 104 ) + a*JVS( 786 )
  a = -W( 95 ) / JVS( 811 )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 812 )
  W( 97 ) = W( 97 ) + a*JVS( 813 )
  W( 98 ) = W( 98 ) + a*JVS( 814 )
  W( 99 ) = W( 99 ) + a*JVS( 815 )
  W( 100 ) = W( 100 ) + a*JVS( 816 )
  W( 101 ) = W( 101 ) + a*JVS( 817 )
  W( 102 ) = W( 102 ) + a*JVS( 818 )
  W( 103 ) = W( 103 ) + a*JVS( 819 )
  W( 104 ) = W( 104 ) + a*JVS( 820 )
  a = -W( 96 ) / JVS( 863 )
  W( 96 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 864 )
  W( 98 ) = W( 98 ) + a*JVS( 865 )
  W( 99 ) = W( 99 ) + a*JVS( 866 )
  W( 100 ) = W( 100 ) + a*JVS( 867 )
  W( 101 ) = W( 101 ) + a*JVS( 868 )
  W( 102 ) = W( 102 ) + a*JVS( 869 )
  W( 103 ) = W( 103 ) + a*JVS( 870 )
  W( 104 ) = W( 104 ) + a*JVS( 871 )
  a = -W( 97 ) / JVS( 889 )
  W( 97 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 890 )
  W( 99 ) = W( 99 ) + a*JVS( 891 )
  W( 100 ) = W( 100 ) + a*JVS( 892 )
  W( 101 ) = W( 101 ) + a*JVS( 893 )
  W( 102 ) = W( 102 ) + a*JVS( 894 )
  W( 103 ) = W( 103 ) + a*JVS( 895 )
  W( 104 ) = W( 104 ) + a*JVS( 896 )
  JVS( 897) = W( 41 )
  JVS( 898) = W( 75 )
  JVS( 899) = W( 78 )
  JVS( 900) = W( 82 )
  JVS( 901) = W( 84 )
  JVS( 902) = W( 85 )
  JVS( 903) = W( 88 )
  JVS( 904) = W( 92 )
  JVS( 905) = W( 93 )
  JVS( 906) = W( 94 )
  JVS( 907) = W( 95 )
  JVS( 908) = W( 96 )
  JVS( 909) = W( 97 )
  JVS( 910) = W( 98 )
  JVS( 911) = W( 99 )
  JVS( 912) = W( 100 )
  JVS( 913) = W( 101 )
  JVS( 914) = W( 102 )
  JVS( 915) = W( 103 )
  JVS( 916) = W( 104 )
  IF ( ABS( JVS( 952 )) < TINY(a) ) THEN
         IER = 99
         RETURN
  END IF
   W( 46 ) = JVS( 917 )
   W( 55 ) = JVS( 918 )
   W( 60 ) = JVS( 919 )
   W( 64 ) = JVS( 920 )
   W( 66 ) = JVS( 921 )
   W( 67 ) = JVS( 922 )
   W( 68 ) = JVS( 923 )
   W( 70 ) = JVS( 924 )
   W( 71 ) = JVS( 925 )
   W( 73 ) = JVS( 926 )
   W( 74 ) = JVS( 927 )
   W( 75 ) = JVS( 928 )
   W( 76 ) = JVS( 929 )
   W( 77 ) = JVS( 930 )
   W( 78 ) = JVS( 931 )
   W( 79 ) = JVS( 932 )
   W( 80 ) = JVS( 933 )
   W( 81 ) = JVS( 934 )
   W( 82 ) = JVS( 935 )
   W( 83 ) = JVS( 936 )
   W( 84 ) = JVS( 937 )
   W( 85 ) = JVS( 938 )
   W( 86 ) = JVS( 939 )
   W( 87 ) = JVS( 940 )
   W( 88 ) = JVS( 941 )
   W( 89 ) = JVS( 942 )
   W( 90 ) = JVS( 943 )
   W( 91 ) = JVS( 944 )
   W( 92 ) = JVS( 945 )
   W( 93 ) = JVS( 946 )
   W( 94 ) = JVS( 947 )
   W( 95 ) = JVS( 948 )
   W( 96 ) = JVS( 949 )
   W( 97 ) = JVS( 950 )
   W( 98 ) = JVS( 951 )
   W( 99 ) = JVS( 952 )
   W( 100 ) = JVS( 953 )
   W( 101 ) = JVS( 954 )
   W( 102 ) = JVS( 955 )
   W( 103 ) = JVS( 956 )
   W( 104 ) = JVS( 957 )
  a = -W( 46 ) / JVS( 246 )
  W( 46 ) = -a
  W( 99 ) = W( 99 ) + a*JVS( 247 )
  W( 101 ) = W( 101 ) + a*JVS( 248 )
  a = -W( 55 ) / JVS( 275 )
  W( 55 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 276 )
  W( 101 ) = W( 101 ) + a*JVS( 277 )
  W( 103 ) = W( 103 ) + a*JVS( 278 )
  a = -W( 60 ) / JVS( 296 )
  W( 60 ) = -a
  W( 99 ) = W( 99 ) + a*JVS( 297 )
  W( 103 ) = W( 103 ) + a*JVS( 298 )
  a = -W( 64 ) / JVS( 310 )
  W( 64 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 311 )
  W( 96 ) = W( 96 ) + a*JVS( 312 )
  W( 99 ) = W( 99 ) + a*JVS( 313 )
  W( 101 ) = W( 101 ) + a*JVS( 314 )
  a = -W( 66 ) / JVS( 322 )
  W( 66 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 323 )
  W( 92 ) = W( 92 ) + a*JVS( 324 )
  W( 99 ) = W( 99 ) + a*JVS( 325 )
  W( 103 ) = W( 103 ) + a*JVS( 326 )
  a = -W( 67 ) / JVS( 328 )
  W( 67 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 329 )
  W( 96 ) = W( 96 ) + a*JVS( 330 )
  W( 99 ) = W( 99 ) + a*JVS( 331 )
  W( 103 ) = W( 103 ) + a*JVS( 332 )
  a = -W( 68 ) / JVS( 338 )
  W( 68 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 339 )
  W( 82 ) = W( 82 ) + a*JVS( 340 )
  W( 85 ) = W( 85 ) + a*JVS( 341 )
  W( 92 ) = W( 92 ) + a*JVS( 342 )
  W( 99 ) = W( 99 ) + a*JVS( 343 )
  W( 103 ) = W( 103 ) + a*JVS( 344 )
  a = -W( 70 ) / JVS( 376 )
  W( 70 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 377 )
  W( 78 ) = W( 78 ) + a*JVS( 378 )
  W( 81 ) = W( 81 ) + a*JVS( 379 )
  W( 82 ) = W( 82 ) + a*JVS( 380 )
  W( 83 ) = W( 83 ) + a*JVS( 381 )
  W( 84 ) = W( 84 ) + a*JVS( 382 )
  W( 85 ) = W( 85 ) + a*JVS( 383 )
  W( 86 ) = W( 86 ) + a*JVS( 384 )
  W( 89 ) = W( 89 ) + a*JVS( 385 )
  W( 92 ) = W( 92 ) + a*JVS( 386 )
  W( 96 ) = W( 96 ) + a*JVS( 387 )
  W( 99 ) = W( 99 ) + a*JVS( 388 )
  W( 101 ) = W( 101 ) + a*JVS( 389 )
  W( 103 ) = W( 103 ) + a*JVS( 390 )
  a = -W( 71 ) / JVS( 391 )
  W( 71 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 392 )
  W( 92 ) = W( 92 ) + a*JVS( 393 )
  W( 99 ) = W( 99 ) + a*JVS( 394 )
  W( 103 ) = W( 103 ) + a*JVS( 395 )
  a = -W( 73 ) / JVS( 411 )
  W( 73 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 92 ) = W( 92 ) + a*JVS( 413 )
  W( 99 ) = W( 99 ) + a*JVS( 414 )
  W( 103 ) = W( 103 ) + a*JVS( 415 )
  a = -W( 74 ) / JVS( 418 )
  W( 74 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 419 )
  W( 95 ) = W( 95 ) + a*JVS( 420 )
  W( 96 ) = W( 96 ) + a*JVS( 421 )
  W( 97 ) = W( 97 ) + a*JVS( 422 )
  W( 98 ) = W( 98 ) + a*JVS( 423 )
  W( 99 ) = W( 99 ) + a*JVS( 424 )
  W( 100 ) = W( 100 ) + a*JVS( 425 )
  W( 101 ) = W( 101 ) + a*JVS( 426 )
  W( 103 ) = W( 103 ) + a*JVS( 427 )
  a = -W( 75 ) / JVS( 428 )
  W( 75 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 429 )
  W( 92 ) = W( 92 ) + a*JVS( 430 )
  W( 99 ) = W( 99 ) + a*JVS( 431 )
  W( 103 ) = W( 103 ) + a*JVS( 432 )
  a = -W( 76 ) / JVS( 433 )
  W( 76 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 434 )
  W( 92 ) = W( 92 ) + a*JVS( 435 )
  W( 99 ) = W( 99 ) + a*JVS( 436 )
  W( 103 ) = W( 103 ) + a*JVS( 437 )
  a = -W( 77 ) / JVS( 446 )
  W( 77 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 447 )
  W( 80 ) = W( 80 ) + a*JVS( 448 )
  W( 84 ) = W( 84 ) + a*JVS( 449 )
  W( 85 ) = W( 85 ) + a*JVS( 450 )
  W( 87 ) = W( 87 ) + a*JVS( 451 )
  W( 88 ) = W( 88 ) + a*JVS( 452 )
  W( 90 ) = W( 90 ) + a*JVS( 453 )
  W( 91 ) = W( 91 ) + a*JVS( 454 )
  W( 92 ) = W( 92 ) + a*JVS( 455 )
  W( 93 ) = W( 93 ) + a*JVS( 456 )
  W( 94 ) = W( 94 ) + a*JVS( 457 )
  W( 95 ) = W( 95 ) + a*JVS( 458 )
  W( 96 ) = W( 96 ) + a*JVS( 459 )
  W( 97 ) = W( 97 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  W( 99 ) = W( 99 ) + a*JVS( 462 )
  W( 100 ) = W( 100 ) + a*JVS( 463 )
  W( 101 ) = W( 101 ) + a*JVS( 464 )
  W( 102 ) = W( 102 ) + a*JVS( 465 )
  W( 103 ) = W( 103 ) + a*JVS( 466 )
  W( 104 ) = W( 104 ) + a*JVS( 467 )
  a = -W( 78 ) / JVS( 469 )
  W( 78 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 470 )
  W( 88 ) = W( 88 ) + a*JVS( 471 )
  W( 92 ) = W( 92 ) + a*JVS( 472 )
  W( 99 ) = W( 99 ) + a*JVS( 473 )
  W( 103 ) = W( 103 ) + a*JVS( 474 )
  a = -W( 79 ) / JVS( 475 )
  W( 79 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 92 ) = W( 92 ) + a*JVS( 477 )
  W( 99 ) = W( 99 ) + a*JVS( 478 )
  W( 103 ) = W( 103 ) + a*JVS( 479 )
  a = -W( 80 ) / JVS( 480 )
  W( 80 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 92 ) = W( 92 ) + a*JVS( 482 )
  W( 99 ) = W( 99 ) + a*JVS( 483 )
  W( 103 ) = W( 103 ) + a*JVS( 484 )
  a = -W( 81 ) / JVS( 496 )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 497 )
  W( 88 ) = W( 88 ) + a*JVS( 498 )
  W( 92 ) = W( 92 ) + a*JVS( 499 )
  W( 93 ) = W( 93 ) + a*JVS( 500 )
  W( 95 ) = W( 95 ) + a*JVS( 501 )
  W( 96 ) = W( 96 ) + a*JVS( 502 )
  W( 97 ) = W( 97 ) + a*JVS( 503 )
  W( 98 ) = W( 98 ) + a*JVS( 504 )
  W( 99 ) = W( 99 ) + a*JVS( 505 )
  W( 100 ) = W( 100 ) + a*JVS( 506 )
  W( 101 ) = W( 101 ) + a*JVS( 507 )
  W( 103 ) = W( 103 ) + a*JVS( 508 )
  a = -W( 82 ) / JVS( 510 )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 511 )
  W( 88 ) = W( 88 ) + a*JVS( 512 )
  W( 92 ) = W( 92 ) + a*JVS( 513 )
  W( 99 ) = W( 99 ) + a*JVS( 514 )
  W( 103 ) = W( 103 ) + a*JVS( 515 )
  a = -W( 83 ) / JVS( 525 )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 526 )
  W( 87 ) = W( 87 ) + a*JVS( 527 )
  W( 88 ) = W( 88 ) + a*JVS( 528 )
  W( 89 ) = W( 89 ) + a*JVS( 529 )
  W( 90 ) = W( 90 ) + a*JVS( 530 )
  W( 91 ) = W( 91 ) + a*JVS( 531 )
  W( 92 ) = W( 92 ) + a*JVS( 532 )
  W( 93 ) = W( 93 ) + a*JVS( 533 )
  W( 95 ) = W( 95 ) + a*JVS( 534 )
  W( 97 ) = W( 97 ) + a*JVS( 535 )
  W( 98 ) = W( 98 ) + a*JVS( 536 )
  W( 99 ) = W( 99 ) + a*JVS( 537 )
  W( 100 ) = W( 100 ) + a*JVS( 538 )
  W( 103 ) = W( 103 ) + a*JVS( 539 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  a = -W( 85 ) / JVS( 547 )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 548 )
  W( 92 ) = W( 92 ) + a*JVS( 549 )
  W( 99 ) = W( 99 ) + a*JVS( 550 )
  W( 103 ) = W( 103 ) + a*JVS( 551 )
  a = -W( 86 ) / JVS( 572 )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 573 )
  W( 88 ) = W( 88 ) + a*JVS( 574 )
  W( 90 ) = W( 90 ) + a*JVS( 575 )
  W( 91 ) = W( 91 ) + a*JVS( 576 )
  W( 92 ) = W( 92 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 94 ) = W( 94 ) + a*JVS( 579 )
  W( 95 ) = W( 95 ) + a*JVS( 580 )
  W( 96 ) = W( 96 ) + a*JVS( 581 )
  W( 97 ) = W( 97 ) + a*JVS( 582 )
  W( 98 ) = W( 98 ) + a*JVS( 583 )
  W( 99 ) = W( 99 ) + a*JVS( 584 )
  W( 100 ) = W( 100 ) + a*JVS( 585 )
  W( 101 ) = W( 101 ) + a*JVS( 586 )
  W( 102 ) = W( 102 ) + a*JVS( 587 )
  W( 103 ) = W( 103 ) + a*JVS( 588 )
  W( 104 ) = W( 104 ) + a*JVS( 589 )
  a = -W( 87 ) / JVS( 596 )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 597 )
  W( 92 ) = W( 92 ) + a*JVS( 598 )
  W( 93 ) = W( 93 ) + a*JVS( 599 )
  W( 99 ) = W( 99 ) + a*JVS( 600 )
  W( 101 ) = W( 101 ) + a*JVS( 601 )
  W( 103 ) = W( 103 ) + a*JVS( 602 )
  W( 104 ) = W( 104 ) + a*JVS( 603 )
  a = -W( 88 ) / JVS( 614 )
  W( 88 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 615 )
  W( 93 ) = W( 93 ) + a*JVS( 616 )
  W( 99 ) = W( 99 ) + a*JVS( 617 )
  W( 101 ) = W( 101 ) + a*JVS( 618 )
  W( 103 ) = W( 103 ) + a*JVS( 619 )
  a = -W( 89 ) / JVS( 639 )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 640 )
  W( 91 ) = W( 91 ) + a*JVS( 641 )
  W( 92 ) = W( 92 ) + a*JVS( 642 )
  W( 93 ) = W( 93 ) + a*JVS( 643 )
  W( 94 ) = W( 94 ) + a*JVS( 644 )
  W( 96 ) = W( 96 ) + a*JVS( 645 )
  W( 99 ) = W( 99 ) + a*JVS( 646 )
  W( 101 ) = W( 101 ) + a*JVS( 647 )
  W( 103 ) = W( 103 ) + a*JVS( 648 )
  W( 104 ) = W( 104 ) + a*JVS( 649 )
  a = -W( 90 ) / JVS( 661 )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 662 )
  W( 92 ) = W( 92 ) + a*JVS( 663 )
  W( 93 ) = W( 93 ) + a*JVS( 664 )
  W( 94 ) = W( 94 ) + a*JVS( 665 )
  W( 99 ) = W( 99 ) + a*JVS( 666 )
  W( 101 ) = W( 101 ) + a*JVS( 667 )
  W( 102 ) = W( 102 ) + a*JVS( 668 )
  W( 103 ) = W( 103 ) + a*JVS( 669 )
  W( 104 ) = W( 104 ) + a*JVS( 670 )
  a = -W( 91 ) / JVS( 683 )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 684 )
  W( 93 ) = W( 93 ) + a*JVS( 685 )
  W( 94 ) = W( 94 ) + a*JVS( 686 )
  W( 95 ) = W( 95 ) + a*JVS( 687 )
  W( 97 ) = W( 97 ) + a*JVS( 688 )
  W( 99 ) = W( 99 ) + a*JVS( 689 )
  W( 100 ) = W( 100 ) + a*JVS( 690 )
  W( 101 ) = W( 101 ) + a*JVS( 691 )
  W( 102 ) = W( 102 ) + a*JVS( 692 )
  W( 103 ) = W( 103 ) + a*JVS( 693 )
  W( 104 ) = W( 104 ) + a*JVS( 694 )
  a = -W( 92 ) / JVS( 708 )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 709 )
  W( 95 ) = W( 95 ) + a*JVS( 710 )
  W( 96 ) = W( 96 ) + a*JVS( 711 )
  W( 97 ) = W( 97 ) + a*JVS( 712 )
  W( 98 ) = W( 98 ) + a*JVS( 713 )
  W( 99 ) = W( 99 ) + a*JVS( 714 )
  W( 100 ) = W( 100 ) + a*JVS( 715 )
  W( 101 ) = W( 101 ) + a*JVS( 716 )
  W( 103 ) = W( 103 ) + a*JVS( 717 )
  a = -W( 93 ) / JVS( 731 )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 732 )
  W( 95 ) = W( 95 ) + a*JVS( 733 )
  W( 96 ) = W( 96 ) + a*JVS( 734 )
  W( 97 ) = W( 97 ) + a*JVS( 735 )
  W( 98 ) = W( 98 ) + a*JVS( 736 )
  W( 99 ) = W( 99 ) + a*JVS( 737 )
  W( 100 ) = W( 100 ) + a*JVS( 738 )
  W( 101 ) = W( 101 ) + a*JVS( 739 )
  W( 102 ) = W( 102 ) + a*JVS( 740 )
  W( 103 ) = W( 103 ) + a*JVS( 741 )
  W( 104 ) = W( 104 ) + a*JVS( 742 )
  a = -W( 94 ) / JVS( 776 )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 777 )
  W( 96 ) = W( 96 ) + a*JVS( 778 )
  W( 97 ) = W( 97 ) + a*JVS( 779 )
  W( 98 ) = W( 98 ) + a*JVS( 780 )
  W( 99 ) = W( 99 ) + a*JVS( 781 )
  W( 100 ) = W( 100 ) + a*JVS( 782 )
  W( 101 ) = W( 101 ) + a*JVS( 783 )
  W( 102 ) = W( 102 ) + a*JVS( 784 )
  W( 103 ) = W( 103 ) + a*JVS( 785 )
  W( 104 ) = W( 104 ) + a*JVS( 786 )
  a = -W( 95 ) / JVS( 811 )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 812 )
  W( 97 ) = W( 97 ) + a*JVS( 813 )
  W( 98 ) = W( 98 ) + a*JVS( 814 )
  W( 99 ) = W( 99 ) + a*JVS( 815 )
  W( 100 ) = W( 100 ) + a*JVS( 816 )
  W( 101 ) = W( 101 ) + a*JVS( 817 )
  W( 102 ) = W( 102 ) + a*JVS( 818 )
  W( 103 ) = W( 103 ) + a*JVS( 819 )
  W( 104 ) = W( 104 ) + a*JVS( 820 )
  a = -W( 96 ) / JVS( 863 )
  W( 96 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 864 )
  W( 98 ) = W( 98 ) + a*JVS( 865 )
  W( 99 ) = W( 99 ) + a*JVS( 866 )
  W( 100 ) = W( 100 ) + a*JVS( 867 )
  W( 101 ) = W( 101 ) + a*JVS( 868 )
  W( 102 ) = W( 102 ) + a*JVS( 869 )
  W( 103 ) = W( 103 ) + a*JVS( 870 )
  W( 104 ) = W( 104 ) + a*JVS( 871 )
  a = -W( 97 ) / JVS( 889 )
  W( 97 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 890 )
  W( 99 ) = W( 99 ) + a*JVS( 891 )
  W( 100 ) = W( 100 ) + a*JVS( 892 )
  W( 101 ) = W( 101 ) + a*JVS( 893 )
  W( 102 ) = W( 102 ) + a*JVS( 894 )
  W( 103 ) = W( 103 ) + a*JVS( 895 )
  W( 104 ) = W( 104 ) + a*JVS( 896 )
  a = -W( 98 ) / JVS( 910 )
  W( 98 ) = -a
  W( 99 ) = W( 99 ) + a*JVS( 911 )
  W( 100 ) = W( 100 ) + a*JVS( 912 )
  W( 101 ) = W( 101 ) + a*JVS( 913 )
  W( 102 ) = W( 102 ) + a*JVS( 914 )
  W( 103 ) = W( 103 ) + a*JVS( 915 )
  W( 104 ) = W( 104 ) + a*JVS( 916 )
  JVS( 917) = W( 46 )
  JVS( 918) = W( 55 )
  JVS( 919) = W( 60 )
  JVS( 920) = W( 64 )
  JVS( 921) = W( 66 )
  JVS( 922) = W( 67 )
  JVS( 923) = W( 68 )
  JVS( 924) = W( 70 )
  JVS( 925) = W( 71 )
  JVS( 926) = W( 73 )
  JVS( 927) = W( 74 )
  JVS( 928) = W( 75 )
  JVS( 929) = W( 76 )
  JVS( 930) = W( 77 )
  JVS( 931) = W( 78 )
  JVS( 932) = W( 79 )
  JVS( 933) = W( 80 )
  JVS( 934) = W( 81 )
  JVS( 935) = W( 82 )
  JVS( 936) = W( 83 )
  JVS( 937) = W( 84 )
  JVS( 938) = W( 85 )
  JVS( 939) = W( 86 )
  JVS( 940) = W( 87 )
  JVS( 941) = W( 88 )
  JVS( 942) = W( 89 )
  JVS( 943) = W( 90 )
  JVS( 944) = W( 91 )
  JVS( 945) = W( 92 )
  JVS( 946) = W( 93 )
  JVS( 947) = W( 94 )
  JVS( 948) = W( 95 )
  JVS( 949) = W( 96 )
  JVS( 950) = W( 97 )
  JVS( 951) = W( 98 )
  JVS( 952) = W( 99 )
  JVS( 953) = W( 100 )
  JVS( 954) = W( 101 )
  JVS( 955) = W( 102 )
  JVS( 956) = W( 103 )
  JVS( 957) = W( 104 )
  IF ( ABS( JVS( 970 )) < TINY(a) ) THEN
         IER = 100
         RETURN
  END IF
   W( 40 ) = JVS( 958 )
   W( 66 ) = JVS( 959 )
   W( 84 ) = JVS( 960 )
   W( 88 ) = JVS( 961 )
   W( 92 ) = JVS( 962 )
   W( 93 ) = JVS( 963 )
   W( 94 ) = JVS( 964 )
   W( 95 ) = JVS( 965 )
   W( 96 ) = JVS( 966 )
   W( 97 ) = JVS( 967 )
   W( 98 ) = JVS( 968 )
   W( 99 ) = JVS( 969 )
   W( 100 ) = JVS( 970 )
   W( 101 ) = JVS( 971 )
   W( 102 ) = JVS( 972 )
   W( 103 ) = JVS( 973 )
   W( 104 ) = JVS( 974 )
  a = -W( 40 ) / JVS( 227 )
  W( 40 ) = -a
  W( 100 ) = W( 100 ) + a*JVS( 228 )
  W( 101 ) = W( 101 ) + a*JVS( 229 )
  a = -W( 66 ) / JVS( 322 )
  W( 66 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 323 )
  W( 92 ) = W( 92 ) + a*JVS( 324 )
  W( 99 ) = W( 99 ) + a*JVS( 325 )
  W( 103 ) = W( 103 ) + a*JVS( 326 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  a = -W( 88 ) / JVS( 614 )
  W( 88 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 615 )
  W( 93 ) = W( 93 ) + a*JVS( 616 )
  W( 99 ) = W( 99 ) + a*JVS( 617 )
  W( 101 ) = W( 101 ) + a*JVS( 618 )
  W( 103 ) = W( 103 ) + a*JVS( 619 )
  a = -W( 92 ) / JVS( 708 )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 709 )
  W( 95 ) = W( 95 ) + a*JVS( 710 )
  W( 96 ) = W( 96 ) + a*JVS( 711 )
  W( 97 ) = W( 97 ) + a*JVS( 712 )
  W( 98 ) = W( 98 ) + a*JVS( 713 )
  W( 99 ) = W( 99 ) + a*JVS( 714 )
  W( 100 ) = W( 100 ) + a*JVS( 715 )
  W( 101 ) = W( 101 ) + a*JVS( 716 )
  W( 103 ) = W( 103 ) + a*JVS( 717 )
  a = -W( 93 ) / JVS( 731 )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 732 )
  W( 95 ) = W( 95 ) + a*JVS( 733 )
  W( 96 ) = W( 96 ) + a*JVS( 734 )
  W( 97 ) = W( 97 ) + a*JVS( 735 )
  W( 98 ) = W( 98 ) + a*JVS( 736 )
  W( 99 ) = W( 99 ) + a*JVS( 737 )
  W( 100 ) = W( 100 ) + a*JVS( 738 )
  W( 101 ) = W( 101 ) + a*JVS( 739 )
  W( 102 ) = W( 102 ) + a*JVS( 740 )
  W( 103 ) = W( 103 ) + a*JVS( 741 )
  W( 104 ) = W( 104 ) + a*JVS( 742 )
  a = -W( 94 ) / JVS( 776 )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 777 )
  W( 96 ) = W( 96 ) + a*JVS( 778 )
  W( 97 ) = W( 97 ) + a*JVS( 779 )
  W( 98 ) = W( 98 ) + a*JVS( 780 )
  W( 99 ) = W( 99 ) + a*JVS( 781 )
  W( 100 ) = W( 100 ) + a*JVS( 782 )
  W( 101 ) = W( 101 ) + a*JVS( 783 )
  W( 102 ) = W( 102 ) + a*JVS( 784 )
  W( 103 ) = W( 103 ) + a*JVS( 785 )
  W( 104 ) = W( 104 ) + a*JVS( 786 )
  a = -W( 95 ) / JVS( 811 )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 812 )
  W( 97 ) = W( 97 ) + a*JVS( 813 )
  W( 98 ) = W( 98 ) + a*JVS( 814 )
  W( 99 ) = W( 99 ) + a*JVS( 815 )
  W( 100 ) = W( 100 ) + a*JVS( 816 )
  W( 101 ) = W( 101 ) + a*JVS( 817 )
  W( 102 ) = W( 102 ) + a*JVS( 818 )
  W( 103 ) = W( 103 ) + a*JVS( 819 )
  W( 104 ) = W( 104 ) + a*JVS( 820 )
  a = -W( 96 ) / JVS( 863 )
  W( 96 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 864 )
  W( 98 ) = W( 98 ) + a*JVS( 865 )
  W( 99 ) = W( 99 ) + a*JVS( 866 )
  W( 100 ) = W( 100 ) + a*JVS( 867 )
  W( 101 ) = W( 101 ) + a*JVS( 868 )
  W( 102 ) = W( 102 ) + a*JVS( 869 )
  W( 103 ) = W( 103 ) + a*JVS( 870 )
  W( 104 ) = W( 104 ) + a*JVS( 871 )
  a = -W( 97 ) / JVS( 889 )
  W( 97 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 890 )
  W( 99 ) = W( 99 ) + a*JVS( 891 )
  W( 100 ) = W( 100 ) + a*JVS( 892 )
  W( 101 ) = W( 101 ) + a*JVS( 893 )
  W( 102 ) = W( 102 ) + a*JVS( 894 )
  W( 103 ) = W( 103 ) + a*JVS( 895 )
  W( 104 ) = W( 104 ) + a*JVS( 896 )
  a = -W( 98 ) / JVS( 910 )
  W( 98 ) = -a
  W( 99 ) = W( 99 ) + a*JVS( 911 )
  W( 100 ) = W( 100 ) + a*JVS( 912 )
  W( 101 ) = W( 101 ) + a*JVS( 913 )
  W( 102 ) = W( 102 ) + a*JVS( 914 )
  W( 103 ) = W( 103 ) + a*JVS( 915 )
  W( 104 ) = W( 104 ) + a*JVS( 916 )
  a = -W( 99 ) / JVS( 952 )
  W( 99 ) = -a
  W( 100 ) = W( 100 ) + a*JVS( 953 )
  W( 101 ) = W( 101 ) + a*JVS( 954 )
  W( 102 ) = W( 102 ) + a*JVS( 955 )
  W( 103 ) = W( 103 ) + a*JVS( 956 )
  W( 104 ) = W( 104 ) + a*JVS( 957 )
  JVS( 958) = W( 40 )
  JVS( 959) = W( 66 )
  JVS( 960) = W( 84 )
  JVS( 961) = W( 88 )
  JVS( 962) = W( 92 )
  JVS( 963) = W( 93 )
  JVS( 964) = W( 94 )
  JVS( 965) = W( 95 )
  JVS( 966) = W( 96 )
  JVS( 967) = W( 97 )
  JVS( 968) = W( 98 )
  JVS( 969) = W( 99 )
  JVS( 970) = W( 100 )
  JVS( 971) = W( 101 )
  JVS( 972) = W( 102 )
  JVS( 973) = W( 103 )
  JVS( 974) = W( 104 )
  IF ( ABS( JVS( 1013 )) < TINY(a) ) THEN
         IER = 101
         RETURN
  END IF
   W( 38 ) = JVS( 975 )
   W( 39 ) = JVS( 976 )
   W( 40 ) = JVS( 977 )
   W( 41 ) = JVS( 978 )
   W( 46 ) = JVS( 979 )
   W( 47 ) = JVS( 980 )
   W( 49 ) = JVS( 981 )
   W( 53 ) = JVS( 982 )
   W( 54 ) = JVS( 983 )
   W( 55 ) = JVS( 984 )
   W( 64 ) = JVS( 985 )
   W( 70 ) = JVS( 986 )
   W( 74 ) = JVS( 987 )
   W( 75 ) = JVS( 988 )
   W( 77 ) = JVS( 989 )
   W( 78 ) = JVS( 990 )
   W( 79 ) = JVS( 991 )
   W( 80 ) = JVS( 992 )
   W( 81 ) = JVS( 993 )
   W( 82 ) = JVS( 994 )
   W( 83 ) = JVS( 995 )
   W( 84 ) = JVS( 996 )
   W( 85 ) = JVS( 997 )
   W( 86 ) = JVS( 998 )
   W( 87 ) = JVS( 999 )
   W( 88 ) = JVS( 1000 )
   W( 89 ) = JVS( 1001 )
   W( 90 ) = JVS( 1002 )
   W( 91 ) = JVS( 1003 )
   W( 92 ) = JVS( 1004 )
   W( 93 ) = JVS( 1005 )
   W( 94 ) = JVS( 1006 )
   W( 95 ) = JVS( 1007 )
   W( 96 ) = JVS( 1008 )
   W( 97 ) = JVS( 1009 )
   W( 98 ) = JVS( 1010 )
   W( 99 ) = JVS( 1011 )
   W( 100 ) = JVS( 1012 )
   W( 101 ) = JVS( 1013 )
   W( 102 ) = JVS( 1014 )
   W( 103 ) = JVS( 1015 )
   W( 104 ) = JVS( 1016 )
  a = -W( 38 ) / JVS( 221 )
  W( 38 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 222 )
  W( 101 ) = W( 101 ) + a*JVS( 223 )
  a = -W( 39 ) / JVS( 224 )
  W( 39 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 225 )
  W( 101 ) = W( 101 ) + a*JVS( 226 )
  a = -W( 40 ) / JVS( 227 )
  W( 40 ) = -a
  W( 100 ) = W( 100 ) + a*JVS( 228 )
  W( 101 ) = W( 101 ) + a*JVS( 229 )
  a = -W( 41 ) / JVS( 230 )
  W( 41 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 231 )
  W( 101 ) = W( 101 ) + a*JVS( 232 )
  a = -W( 46 ) / JVS( 246 )
  W( 46 ) = -a
  W( 99 ) = W( 99 ) + a*JVS( 247 )
  W( 101 ) = W( 101 ) + a*JVS( 248 )
  a = -W( 47 ) / JVS( 249 )
  W( 47 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 250 )
  W( 103 ) = W( 103 ) + a*JVS( 251 )
  a = -W( 49 ) / JVS( 255 )
  W( 49 ) = -a
  W( 101 ) = W( 101 ) + a*JVS( 256 )
  W( 103 ) = W( 103 ) + a*JVS( 257 )
  a = -W( 53 ) / JVS( 266 )
  W( 53 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 267 )
  W( 93 ) = W( 93 ) + a*JVS( 268 )
  W( 96 ) = W( 96 ) + a*JVS( 269 )
  a = -W( 54 ) / JVS( 270 )
  W( 54 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 271 )
  W( 96 ) = W( 96 ) + a*JVS( 272 )
  W( 99 ) = W( 99 ) + a*JVS( 273 )
  W( 101 ) = W( 101 ) + a*JVS( 274 )
  a = -W( 55 ) / JVS( 275 )
  W( 55 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 276 )
  W( 101 ) = W( 101 ) + a*JVS( 277 )
  W( 103 ) = W( 103 ) + a*JVS( 278 )
  a = -W( 64 ) / JVS( 310 )
  W( 64 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 311 )
  W( 96 ) = W( 96 ) + a*JVS( 312 )
  W( 99 ) = W( 99 ) + a*JVS( 313 )
  W( 101 ) = W( 101 ) + a*JVS( 314 )
  a = -W( 70 ) / JVS( 376 )
  W( 70 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 377 )
  W( 78 ) = W( 78 ) + a*JVS( 378 )
  W( 81 ) = W( 81 ) + a*JVS( 379 )
  W( 82 ) = W( 82 ) + a*JVS( 380 )
  W( 83 ) = W( 83 ) + a*JVS( 381 )
  W( 84 ) = W( 84 ) + a*JVS( 382 )
  W( 85 ) = W( 85 ) + a*JVS( 383 )
  W( 86 ) = W( 86 ) + a*JVS( 384 )
  W( 89 ) = W( 89 ) + a*JVS( 385 )
  W( 92 ) = W( 92 ) + a*JVS( 386 )
  W( 96 ) = W( 96 ) + a*JVS( 387 )
  W( 99 ) = W( 99 ) + a*JVS( 388 )
  W( 101 ) = W( 101 ) + a*JVS( 389 )
  W( 103 ) = W( 103 ) + a*JVS( 390 )
  a = -W( 74 ) / JVS( 418 )
  W( 74 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 419 )
  W( 95 ) = W( 95 ) + a*JVS( 420 )
  W( 96 ) = W( 96 ) + a*JVS( 421 )
  W( 97 ) = W( 97 ) + a*JVS( 422 )
  W( 98 ) = W( 98 ) + a*JVS( 423 )
  W( 99 ) = W( 99 ) + a*JVS( 424 )
  W( 100 ) = W( 100 ) + a*JVS( 425 )
  W( 101 ) = W( 101 ) + a*JVS( 426 )
  W( 103 ) = W( 103 ) + a*JVS( 427 )
  a = -W( 75 ) / JVS( 428 )
  W( 75 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 429 )
  W( 92 ) = W( 92 ) + a*JVS( 430 )
  W( 99 ) = W( 99 ) + a*JVS( 431 )
  W( 103 ) = W( 103 ) + a*JVS( 432 )
  a = -W( 77 ) / JVS( 446 )
  W( 77 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 447 )
  W( 80 ) = W( 80 ) + a*JVS( 448 )
  W( 84 ) = W( 84 ) + a*JVS( 449 )
  W( 85 ) = W( 85 ) + a*JVS( 450 )
  W( 87 ) = W( 87 ) + a*JVS( 451 )
  W( 88 ) = W( 88 ) + a*JVS( 452 )
  W( 90 ) = W( 90 ) + a*JVS( 453 )
  W( 91 ) = W( 91 ) + a*JVS( 454 )
  W( 92 ) = W( 92 ) + a*JVS( 455 )
  W( 93 ) = W( 93 ) + a*JVS( 456 )
  W( 94 ) = W( 94 ) + a*JVS( 457 )
  W( 95 ) = W( 95 ) + a*JVS( 458 )
  W( 96 ) = W( 96 ) + a*JVS( 459 )
  W( 97 ) = W( 97 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  W( 99 ) = W( 99 ) + a*JVS( 462 )
  W( 100 ) = W( 100 ) + a*JVS( 463 )
  W( 101 ) = W( 101 ) + a*JVS( 464 )
  W( 102 ) = W( 102 ) + a*JVS( 465 )
  W( 103 ) = W( 103 ) + a*JVS( 466 )
  W( 104 ) = W( 104 ) + a*JVS( 467 )
  a = -W( 78 ) / JVS( 469 )
  W( 78 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 470 )
  W( 88 ) = W( 88 ) + a*JVS( 471 )
  W( 92 ) = W( 92 ) + a*JVS( 472 )
  W( 99 ) = W( 99 ) + a*JVS( 473 )
  W( 103 ) = W( 103 ) + a*JVS( 474 )
  a = -W( 79 ) / JVS( 475 )
  W( 79 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 92 ) = W( 92 ) + a*JVS( 477 )
  W( 99 ) = W( 99 ) + a*JVS( 478 )
  W( 103 ) = W( 103 ) + a*JVS( 479 )
  a = -W( 80 ) / JVS( 480 )
  W( 80 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 92 ) = W( 92 ) + a*JVS( 482 )
  W( 99 ) = W( 99 ) + a*JVS( 483 )
  W( 103 ) = W( 103 ) + a*JVS( 484 )
  a = -W( 81 ) / JVS( 496 )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 497 )
  W( 88 ) = W( 88 ) + a*JVS( 498 )
  W( 92 ) = W( 92 ) + a*JVS( 499 )
  W( 93 ) = W( 93 ) + a*JVS( 500 )
  W( 95 ) = W( 95 ) + a*JVS( 501 )
  W( 96 ) = W( 96 ) + a*JVS( 502 )
  W( 97 ) = W( 97 ) + a*JVS( 503 )
  W( 98 ) = W( 98 ) + a*JVS( 504 )
  W( 99 ) = W( 99 ) + a*JVS( 505 )
  W( 100 ) = W( 100 ) + a*JVS( 506 )
  W( 101 ) = W( 101 ) + a*JVS( 507 )
  W( 103 ) = W( 103 ) + a*JVS( 508 )
  a = -W( 82 ) / JVS( 510 )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 511 )
  W( 88 ) = W( 88 ) + a*JVS( 512 )
  W( 92 ) = W( 92 ) + a*JVS( 513 )
  W( 99 ) = W( 99 ) + a*JVS( 514 )
  W( 103 ) = W( 103 ) + a*JVS( 515 )
  a = -W( 83 ) / JVS( 525 )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 526 )
  W( 87 ) = W( 87 ) + a*JVS( 527 )
  W( 88 ) = W( 88 ) + a*JVS( 528 )
  W( 89 ) = W( 89 ) + a*JVS( 529 )
  W( 90 ) = W( 90 ) + a*JVS( 530 )
  W( 91 ) = W( 91 ) + a*JVS( 531 )
  W( 92 ) = W( 92 ) + a*JVS( 532 )
  W( 93 ) = W( 93 ) + a*JVS( 533 )
  W( 95 ) = W( 95 ) + a*JVS( 534 )
  W( 97 ) = W( 97 ) + a*JVS( 535 )
  W( 98 ) = W( 98 ) + a*JVS( 536 )
  W( 99 ) = W( 99 ) + a*JVS( 537 )
  W( 100 ) = W( 100 ) + a*JVS( 538 )
  W( 103 ) = W( 103 ) + a*JVS( 539 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  a = -W( 85 ) / JVS( 547 )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 548 )
  W( 92 ) = W( 92 ) + a*JVS( 549 )
  W( 99 ) = W( 99 ) + a*JVS( 550 )
  W( 103 ) = W( 103 ) + a*JVS( 551 )
  a = -W( 86 ) / JVS( 572 )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 573 )
  W( 88 ) = W( 88 ) + a*JVS( 574 )
  W( 90 ) = W( 90 ) + a*JVS( 575 )
  W( 91 ) = W( 91 ) + a*JVS( 576 )
  W( 92 ) = W( 92 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 94 ) = W( 94 ) + a*JVS( 579 )
  W( 95 ) = W( 95 ) + a*JVS( 580 )
  W( 96 ) = W( 96 ) + a*JVS( 581 )
  W( 97 ) = W( 97 ) + a*JVS( 582 )
  W( 98 ) = W( 98 ) + a*JVS( 583 )
  W( 99 ) = W( 99 ) + a*JVS( 584 )
  W( 100 ) = W( 100 ) + a*JVS( 585 )
  W( 101 ) = W( 101 ) + a*JVS( 586 )
  W( 102 ) = W( 102 ) + a*JVS( 587 )
  W( 103 ) = W( 103 ) + a*JVS( 588 )
  W( 104 ) = W( 104 ) + a*JVS( 589 )
  a = -W( 87 ) / JVS( 596 )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 597 )
  W( 92 ) = W( 92 ) + a*JVS( 598 )
  W( 93 ) = W( 93 ) + a*JVS( 599 )
  W( 99 ) = W( 99 ) + a*JVS( 600 )
  W( 101 ) = W( 101 ) + a*JVS( 601 )
  W( 103 ) = W( 103 ) + a*JVS( 602 )
  W( 104 ) = W( 104 ) + a*JVS( 603 )
  a = -W( 88 ) / JVS( 614 )
  W( 88 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 615 )
  W( 93 ) = W( 93 ) + a*JVS( 616 )
  W( 99 ) = W( 99 ) + a*JVS( 617 )
  W( 101 ) = W( 101 ) + a*JVS( 618 )
  W( 103 ) = W( 103 ) + a*JVS( 619 )
  a = -W( 89 ) / JVS( 639 )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 640 )
  W( 91 ) = W( 91 ) + a*JVS( 641 )
  W( 92 ) = W( 92 ) + a*JVS( 642 )
  W( 93 ) = W( 93 ) + a*JVS( 643 )
  W( 94 ) = W( 94 ) + a*JVS( 644 )
  W( 96 ) = W( 96 ) + a*JVS( 645 )
  W( 99 ) = W( 99 ) + a*JVS( 646 )
  W( 101 ) = W( 101 ) + a*JVS( 647 )
  W( 103 ) = W( 103 ) + a*JVS( 648 )
  W( 104 ) = W( 104 ) + a*JVS( 649 )
  a = -W( 90 ) / JVS( 661 )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 662 )
  W( 92 ) = W( 92 ) + a*JVS( 663 )
  W( 93 ) = W( 93 ) + a*JVS( 664 )
  W( 94 ) = W( 94 ) + a*JVS( 665 )
  W( 99 ) = W( 99 ) + a*JVS( 666 )
  W( 101 ) = W( 101 ) + a*JVS( 667 )
  W( 102 ) = W( 102 ) + a*JVS( 668 )
  W( 103 ) = W( 103 ) + a*JVS( 669 )
  W( 104 ) = W( 104 ) + a*JVS( 670 )
  a = -W( 91 ) / JVS( 683 )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 684 )
  W( 93 ) = W( 93 ) + a*JVS( 685 )
  W( 94 ) = W( 94 ) + a*JVS( 686 )
  W( 95 ) = W( 95 ) + a*JVS( 687 )
  W( 97 ) = W( 97 ) + a*JVS( 688 )
  W( 99 ) = W( 99 ) + a*JVS( 689 )
  W( 100 ) = W( 100 ) + a*JVS( 690 )
  W( 101 ) = W( 101 ) + a*JVS( 691 )
  W( 102 ) = W( 102 ) + a*JVS( 692 )
  W( 103 ) = W( 103 ) + a*JVS( 693 )
  W( 104 ) = W( 104 ) + a*JVS( 694 )
  a = -W( 92 ) / JVS( 708 )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 709 )
  W( 95 ) = W( 95 ) + a*JVS( 710 )
  W( 96 ) = W( 96 ) + a*JVS( 711 )
  W( 97 ) = W( 97 ) + a*JVS( 712 )
  W( 98 ) = W( 98 ) + a*JVS( 713 )
  W( 99 ) = W( 99 ) + a*JVS( 714 )
  W( 100 ) = W( 100 ) + a*JVS( 715 )
  W( 101 ) = W( 101 ) + a*JVS( 716 )
  W( 103 ) = W( 103 ) + a*JVS( 717 )
  a = -W( 93 ) / JVS( 731 )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 732 )
  W( 95 ) = W( 95 ) + a*JVS( 733 )
  W( 96 ) = W( 96 ) + a*JVS( 734 )
  W( 97 ) = W( 97 ) + a*JVS( 735 )
  W( 98 ) = W( 98 ) + a*JVS( 736 )
  W( 99 ) = W( 99 ) + a*JVS( 737 )
  W( 100 ) = W( 100 ) + a*JVS( 738 )
  W( 101 ) = W( 101 ) + a*JVS( 739 )
  W( 102 ) = W( 102 ) + a*JVS( 740 )
  W( 103 ) = W( 103 ) + a*JVS( 741 )
  W( 104 ) = W( 104 ) + a*JVS( 742 )
  a = -W( 94 ) / JVS( 776 )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 777 )
  W( 96 ) = W( 96 ) + a*JVS( 778 )
  W( 97 ) = W( 97 ) + a*JVS( 779 )
  W( 98 ) = W( 98 ) + a*JVS( 780 )
  W( 99 ) = W( 99 ) + a*JVS( 781 )
  W( 100 ) = W( 100 ) + a*JVS( 782 )
  W( 101 ) = W( 101 ) + a*JVS( 783 )
  W( 102 ) = W( 102 ) + a*JVS( 784 )
  W( 103 ) = W( 103 ) + a*JVS( 785 )
  W( 104 ) = W( 104 ) + a*JVS( 786 )
  a = -W( 95 ) / JVS( 811 )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 812 )
  W( 97 ) = W( 97 ) + a*JVS( 813 )
  W( 98 ) = W( 98 ) + a*JVS( 814 )
  W( 99 ) = W( 99 ) + a*JVS( 815 )
  W( 100 ) = W( 100 ) + a*JVS( 816 )
  W( 101 ) = W( 101 ) + a*JVS( 817 )
  W( 102 ) = W( 102 ) + a*JVS( 818 )
  W( 103 ) = W( 103 ) + a*JVS( 819 )
  W( 104 ) = W( 104 ) + a*JVS( 820 )
  a = -W( 96 ) / JVS( 863 )
  W( 96 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 864 )
  W( 98 ) = W( 98 ) + a*JVS( 865 )
  W( 99 ) = W( 99 ) + a*JVS( 866 )
  W( 100 ) = W( 100 ) + a*JVS( 867 )
  W( 101 ) = W( 101 ) + a*JVS( 868 )
  W( 102 ) = W( 102 ) + a*JVS( 869 )
  W( 103 ) = W( 103 ) + a*JVS( 870 )
  W( 104 ) = W( 104 ) + a*JVS( 871 )
  a = -W( 97 ) / JVS( 889 )
  W( 97 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 890 )
  W( 99 ) = W( 99 ) + a*JVS( 891 )
  W( 100 ) = W( 100 ) + a*JVS( 892 )
  W( 101 ) = W( 101 ) + a*JVS( 893 )
  W( 102 ) = W( 102 ) + a*JVS( 894 )
  W( 103 ) = W( 103 ) + a*JVS( 895 )
  W( 104 ) = W( 104 ) + a*JVS( 896 )
  a = -W( 98 ) / JVS( 910 )
  W( 98 ) = -a
  W( 99 ) = W( 99 ) + a*JVS( 911 )
  W( 100 ) = W( 100 ) + a*JVS( 912 )
  W( 101 ) = W( 101 ) + a*JVS( 913 )
  W( 102 ) = W( 102 ) + a*JVS( 914 )
  W( 103 ) = W( 103 ) + a*JVS( 915 )
  W( 104 ) = W( 104 ) + a*JVS( 916 )
  a = -W( 99 ) / JVS( 952 )
  W( 99 ) = -a
  W( 100 ) = W( 100 ) + a*JVS( 953 )
  W( 101 ) = W( 101 ) + a*JVS( 954 )
  W( 102 ) = W( 102 ) + a*JVS( 955 )
  W( 103 ) = W( 103 ) + a*JVS( 956 )
  W( 104 ) = W( 104 ) + a*JVS( 957 )
  a = -W( 100 ) / JVS( 970 )
  W( 100 ) = -a
  W( 101 ) = W( 101 ) + a*JVS( 971 )
  W( 102 ) = W( 102 ) + a*JVS( 972 )
  W( 103 ) = W( 103 ) + a*JVS( 973 )
  W( 104 ) = W( 104 ) + a*JVS( 974 )
  JVS( 975) = W( 38 )
  JVS( 976) = W( 39 )
  JVS( 977) = W( 40 )
  JVS( 978) = W( 41 )
  JVS( 979) = W( 46 )
  JVS( 980) = W( 47 )
  JVS( 981) = W( 49 )
  JVS( 982) = W( 53 )
  JVS( 983) = W( 54 )
  JVS( 984) = W( 55 )
  JVS( 985) = W( 64 )
  JVS( 986) = W( 70 )
  JVS( 987) = W( 74 )
  JVS( 988) = W( 75 )
  JVS( 989) = W( 77 )
  JVS( 990) = W( 78 )
  JVS( 991) = W( 79 )
  JVS( 992) = W( 80 )
  JVS( 993) = W( 81 )
  JVS( 994) = W( 82 )
  JVS( 995) = W( 83 )
  JVS( 996) = W( 84 )
  JVS( 997) = W( 85 )
  JVS( 998) = W( 86 )
  JVS( 999) = W( 87 )
  JVS( 1000) = W( 88 )
  JVS( 1001) = W( 89 )
  JVS( 1002) = W( 90 )
  JVS( 1003) = W( 91 )
  JVS( 1004) = W( 92 )
  JVS( 1005) = W( 93 )
  JVS( 1006) = W( 94 )
  JVS( 1007) = W( 95 )
  JVS( 1008) = W( 96 )
  JVS( 1009) = W( 97 )
  JVS( 1010) = W( 98 )
  JVS( 1011) = W( 99 )
  JVS( 1012) = W( 100 )
  JVS( 1013) = W( 101 )
  JVS( 1014) = W( 102 )
  JVS( 1015) = W( 103 )
  JVS( 1016) = W( 104 )
  IF ( ABS( JVS( 1046 )) < TINY(a) ) THEN
         IER = 102
         RETURN
  END IF
   W( 35 ) = JVS( 1017 )
   W( 49 ) = JVS( 1018 )
   W( 52 ) = JVS( 1019 )
   W( 61 ) = JVS( 1020 )
   W( 71 ) = JVS( 1021 )
   W( 72 ) = JVS( 1022 )
   W( 73 ) = JVS( 1023 )
   W( 75 ) = JVS( 1024 )
   W( 76 ) = JVS( 1025 )
   W( 79 ) = JVS( 1026 )
   W( 80 ) = JVS( 1027 )
   W( 83 ) = JVS( 1028 )
   W( 84 ) = JVS( 1029 )
   W( 85 ) = JVS( 1030 )
   W( 87 ) = JVS( 1031 )
   W( 88 ) = JVS( 1032 )
   W( 89 ) = JVS( 1033 )
   W( 90 ) = JVS( 1034 )
   W( 91 ) = JVS( 1035 )
   W( 92 ) = JVS( 1036 )
   W( 93 ) = JVS( 1037 )
   W( 94 ) = JVS( 1038 )
   W( 95 ) = JVS( 1039 )
   W( 96 ) = JVS( 1040 )
   W( 97 ) = JVS( 1041 )
   W( 98 ) = JVS( 1042 )
   W( 99 ) = JVS( 1043 )
   W( 100 ) = JVS( 1044 )
   W( 101 ) = JVS( 1045 )
   W( 102 ) = JVS( 1046 )
   W( 103 ) = JVS( 1047 )
   W( 104 ) = JVS( 1048 )
  a = -W( 35 ) / JVS( 215 )
  W( 35 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 216 )
  a = -W( 49 ) / JVS( 255 )
  W( 49 ) = -a
  W( 101 ) = W( 101 ) + a*JVS( 256 )
  W( 103 ) = W( 103 ) + a*JVS( 257 )
  a = -W( 52 ) / JVS( 262 )
  W( 52 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 263 )
  W( 102 ) = W( 102 ) + a*JVS( 264 )
  W( 103 ) = W( 103 ) + a*JVS( 265 )
  a = -W( 61 ) / JVS( 299 )
  W( 61 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 300 )
  a = -W( 71 ) / JVS( 391 )
  W( 71 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 392 )
  W( 92 ) = W( 92 ) + a*JVS( 393 )
  W( 99 ) = W( 99 ) + a*JVS( 394 )
  W( 103 ) = W( 103 ) + a*JVS( 395 )
  a = -W( 72 ) / JVS( 401 )
  W( 72 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 402 )
  W( 79 ) = W( 79 ) + a*JVS( 403 )
  W( 80 ) = W( 80 ) + a*JVS( 404 )
  W( 84 ) = W( 84 ) + a*JVS( 405 )
  W( 87 ) = W( 87 ) + a*JVS( 406 )
  W( 92 ) = W( 92 ) + a*JVS( 407 )
  W( 99 ) = W( 99 ) + a*JVS( 408 )
  W( 101 ) = W( 101 ) + a*JVS( 409 )
  W( 103 ) = W( 103 ) + a*JVS( 410 )
  a = -W( 73 ) / JVS( 411 )
  W( 73 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 92 ) = W( 92 ) + a*JVS( 413 )
  W( 99 ) = W( 99 ) + a*JVS( 414 )
  W( 103 ) = W( 103 ) + a*JVS( 415 )
  a = -W( 75 ) / JVS( 428 )
  W( 75 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 429 )
  W( 92 ) = W( 92 ) + a*JVS( 430 )
  W( 99 ) = W( 99 ) + a*JVS( 431 )
  W( 103 ) = W( 103 ) + a*JVS( 432 )
  a = -W( 76 ) / JVS( 433 )
  W( 76 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 434 )
  W( 92 ) = W( 92 ) + a*JVS( 435 )
  W( 99 ) = W( 99 ) + a*JVS( 436 )
  W( 103 ) = W( 103 ) + a*JVS( 437 )
  a = -W( 79 ) / JVS( 475 )
  W( 79 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 92 ) = W( 92 ) + a*JVS( 477 )
  W( 99 ) = W( 99 ) + a*JVS( 478 )
  W( 103 ) = W( 103 ) + a*JVS( 479 )
  a = -W( 80 ) / JVS( 480 )
  W( 80 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 92 ) = W( 92 ) + a*JVS( 482 )
  W( 99 ) = W( 99 ) + a*JVS( 483 )
  W( 103 ) = W( 103 ) + a*JVS( 484 )
  a = -W( 83 ) / JVS( 525 )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 526 )
  W( 87 ) = W( 87 ) + a*JVS( 527 )
  W( 88 ) = W( 88 ) + a*JVS( 528 )
  W( 89 ) = W( 89 ) + a*JVS( 529 )
  W( 90 ) = W( 90 ) + a*JVS( 530 )
  W( 91 ) = W( 91 ) + a*JVS( 531 )
  W( 92 ) = W( 92 ) + a*JVS( 532 )
  W( 93 ) = W( 93 ) + a*JVS( 533 )
  W( 95 ) = W( 95 ) + a*JVS( 534 )
  W( 97 ) = W( 97 ) + a*JVS( 535 )
  W( 98 ) = W( 98 ) + a*JVS( 536 )
  W( 99 ) = W( 99 ) + a*JVS( 537 )
  W( 100 ) = W( 100 ) + a*JVS( 538 )
  W( 103 ) = W( 103 ) + a*JVS( 539 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  a = -W( 85 ) / JVS( 547 )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 548 )
  W( 92 ) = W( 92 ) + a*JVS( 549 )
  W( 99 ) = W( 99 ) + a*JVS( 550 )
  W( 103 ) = W( 103 ) + a*JVS( 551 )
  a = -W( 87 ) / JVS( 596 )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 597 )
  W( 92 ) = W( 92 ) + a*JVS( 598 )
  W( 93 ) = W( 93 ) + a*JVS( 599 )
  W( 99 ) = W( 99 ) + a*JVS( 600 )
  W( 101 ) = W( 101 ) + a*JVS( 601 )
  W( 103 ) = W( 103 ) + a*JVS( 602 )
  W( 104 ) = W( 104 ) + a*JVS( 603 )
  a = -W( 88 ) / JVS( 614 )
  W( 88 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 615 )
  W( 93 ) = W( 93 ) + a*JVS( 616 )
  W( 99 ) = W( 99 ) + a*JVS( 617 )
  W( 101 ) = W( 101 ) + a*JVS( 618 )
  W( 103 ) = W( 103 ) + a*JVS( 619 )
  a = -W( 89 ) / JVS( 639 )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 640 )
  W( 91 ) = W( 91 ) + a*JVS( 641 )
  W( 92 ) = W( 92 ) + a*JVS( 642 )
  W( 93 ) = W( 93 ) + a*JVS( 643 )
  W( 94 ) = W( 94 ) + a*JVS( 644 )
  W( 96 ) = W( 96 ) + a*JVS( 645 )
  W( 99 ) = W( 99 ) + a*JVS( 646 )
  W( 101 ) = W( 101 ) + a*JVS( 647 )
  W( 103 ) = W( 103 ) + a*JVS( 648 )
  W( 104 ) = W( 104 ) + a*JVS( 649 )
  a = -W( 90 ) / JVS( 661 )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 662 )
  W( 92 ) = W( 92 ) + a*JVS( 663 )
  W( 93 ) = W( 93 ) + a*JVS( 664 )
  W( 94 ) = W( 94 ) + a*JVS( 665 )
  W( 99 ) = W( 99 ) + a*JVS( 666 )
  W( 101 ) = W( 101 ) + a*JVS( 667 )
  W( 102 ) = W( 102 ) + a*JVS( 668 )
  W( 103 ) = W( 103 ) + a*JVS( 669 )
  W( 104 ) = W( 104 ) + a*JVS( 670 )
  a = -W( 91 ) / JVS( 683 )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 684 )
  W( 93 ) = W( 93 ) + a*JVS( 685 )
  W( 94 ) = W( 94 ) + a*JVS( 686 )
  W( 95 ) = W( 95 ) + a*JVS( 687 )
  W( 97 ) = W( 97 ) + a*JVS( 688 )
  W( 99 ) = W( 99 ) + a*JVS( 689 )
  W( 100 ) = W( 100 ) + a*JVS( 690 )
  W( 101 ) = W( 101 ) + a*JVS( 691 )
  W( 102 ) = W( 102 ) + a*JVS( 692 )
  W( 103 ) = W( 103 ) + a*JVS( 693 )
  W( 104 ) = W( 104 ) + a*JVS( 694 )
  a = -W( 92 ) / JVS( 708 )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 709 )
  W( 95 ) = W( 95 ) + a*JVS( 710 )
  W( 96 ) = W( 96 ) + a*JVS( 711 )
  W( 97 ) = W( 97 ) + a*JVS( 712 )
  W( 98 ) = W( 98 ) + a*JVS( 713 )
  W( 99 ) = W( 99 ) + a*JVS( 714 )
  W( 100 ) = W( 100 ) + a*JVS( 715 )
  W( 101 ) = W( 101 ) + a*JVS( 716 )
  W( 103 ) = W( 103 ) + a*JVS( 717 )
  a = -W( 93 ) / JVS( 731 )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 732 )
  W( 95 ) = W( 95 ) + a*JVS( 733 )
  W( 96 ) = W( 96 ) + a*JVS( 734 )
  W( 97 ) = W( 97 ) + a*JVS( 735 )
  W( 98 ) = W( 98 ) + a*JVS( 736 )
  W( 99 ) = W( 99 ) + a*JVS( 737 )
  W( 100 ) = W( 100 ) + a*JVS( 738 )
  W( 101 ) = W( 101 ) + a*JVS( 739 )
  W( 102 ) = W( 102 ) + a*JVS( 740 )
  W( 103 ) = W( 103 ) + a*JVS( 741 )
  W( 104 ) = W( 104 ) + a*JVS( 742 )
  a = -W( 94 ) / JVS( 776 )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 777 )
  W( 96 ) = W( 96 ) + a*JVS( 778 )
  W( 97 ) = W( 97 ) + a*JVS( 779 )
  W( 98 ) = W( 98 ) + a*JVS( 780 )
  W( 99 ) = W( 99 ) + a*JVS( 781 )
  W( 100 ) = W( 100 ) + a*JVS( 782 )
  W( 101 ) = W( 101 ) + a*JVS( 783 )
  W( 102 ) = W( 102 ) + a*JVS( 784 )
  W( 103 ) = W( 103 ) + a*JVS( 785 )
  W( 104 ) = W( 104 ) + a*JVS( 786 )
  a = -W( 95 ) / JVS( 811 )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 812 )
  W( 97 ) = W( 97 ) + a*JVS( 813 )
  W( 98 ) = W( 98 ) + a*JVS( 814 )
  W( 99 ) = W( 99 ) + a*JVS( 815 )
  W( 100 ) = W( 100 ) + a*JVS( 816 )
  W( 101 ) = W( 101 ) + a*JVS( 817 )
  W( 102 ) = W( 102 ) + a*JVS( 818 )
  W( 103 ) = W( 103 ) + a*JVS( 819 )
  W( 104 ) = W( 104 ) + a*JVS( 820 )
  a = -W( 96 ) / JVS( 863 )
  W( 96 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 864 )
  W( 98 ) = W( 98 ) + a*JVS( 865 )
  W( 99 ) = W( 99 ) + a*JVS( 866 )
  W( 100 ) = W( 100 ) + a*JVS( 867 )
  W( 101 ) = W( 101 ) + a*JVS( 868 )
  W( 102 ) = W( 102 ) + a*JVS( 869 )
  W( 103 ) = W( 103 ) + a*JVS( 870 )
  W( 104 ) = W( 104 ) + a*JVS( 871 )
  a = -W( 97 ) / JVS( 889 )
  W( 97 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 890 )
  W( 99 ) = W( 99 ) + a*JVS( 891 )
  W( 100 ) = W( 100 ) + a*JVS( 892 )
  W( 101 ) = W( 101 ) + a*JVS( 893 )
  W( 102 ) = W( 102 ) + a*JVS( 894 )
  W( 103 ) = W( 103 ) + a*JVS( 895 )
  W( 104 ) = W( 104 ) + a*JVS( 896 )
  a = -W( 98 ) / JVS( 910 )
  W( 98 ) = -a
  W( 99 ) = W( 99 ) + a*JVS( 911 )
  W( 100 ) = W( 100 ) + a*JVS( 912 )
  W( 101 ) = W( 101 ) + a*JVS( 913 )
  W( 102 ) = W( 102 ) + a*JVS( 914 )
  W( 103 ) = W( 103 ) + a*JVS( 915 )
  W( 104 ) = W( 104 ) + a*JVS( 916 )
  a = -W( 99 ) / JVS( 952 )
  W( 99 ) = -a
  W( 100 ) = W( 100 ) + a*JVS( 953 )
  W( 101 ) = W( 101 ) + a*JVS( 954 )
  W( 102 ) = W( 102 ) + a*JVS( 955 )
  W( 103 ) = W( 103 ) + a*JVS( 956 )
  W( 104 ) = W( 104 ) + a*JVS( 957 )
  a = -W( 100 ) / JVS( 970 )
  W( 100 ) = -a
  W( 101 ) = W( 101 ) + a*JVS( 971 )
  W( 102 ) = W( 102 ) + a*JVS( 972 )
  W( 103 ) = W( 103 ) + a*JVS( 973 )
  W( 104 ) = W( 104 ) + a*JVS( 974 )
  a = -W( 101 ) / JVS( 1013 )
  W( 101 ) = -a
  W( 102 ) = W( 102 ) + a*JVS( 1014 )
  W( 103 ) = W( 103 ) + a*JVS( 1015 )
  W( 104 ) = W( 104 ) + a*JVS( 1016 )
  JVS( 1017) = W( 35 )
  JVS( 1018) = W( 49 )
  JVS( 1019) = W( 52 )
  JVS( 1020) = W( 61 )
  JVS( 1021) = W( 71 )
  JVS( 1022) = W( 72 )
  JVS( 1023) = W( 73 )
  JVS( 1024) = W( 75 )
  JVS( 1025) = W( 76 )
  JVS( 1026) = W( 79 )
  JVS( 1027) = W( 80 )
  JVS( 1028) = W( 83 )
  JVS( 1029) = W( 84 )
  JVS( 1030) = W( 85 )
  JVS( 1031) = W( 87 )
  JVS( 1032) = W( 88 )
  JVS( 1033) = W( 89 )
  JVS( 1034) = W( 90 )
  JVS( 1035) = W( 91 )
  JVS( 1036) = W( 92 )
  JVS( 1037) = W( 93 )
  JVS( 1038) = W( 94 )
  JVS( 1039) = W( 95 )
  JVS( 1040) = W( 96 )
  JVS( 1041) = W( 97 )
  JVS( 1042) = W( 98 )
  JVS( 1043) = W( 99 )
  JVS( 1044) = W( 100 )
  JVS( 1045) = W( 101 )
  JVS( 1046) = W( 102 )
  JVS( 1047) = W( 103 )
  JVS( 1048) = W( 104 )
  IF ( ABS( JVS( 1113 )) < TINY(a) ) THEN
         IER = 103
         RETURN
  END IF
   W( 3 ) = JVS( 1049 )
   W( 5 ) = JVS( 1050 )
   W( 30 ) = JVS( 1051 )
   W( 31 ) = JVS( 1052 )
   W( 32 ) = JVS( 1053 )
   W( 33 ) = JVS( 1054 )
   W( 34 ) = JVS( 1055 )
   W( 35 ) = JVS( 1056 )
   W( 36 ) = JVS( 1057 )
   W( 37 ) = JVS( 1058 )
   W( 42 ) = JVS( 1059 )
   W( 43 ) = JVS( 1060 )
   W( 44 ) = JVS( 1061 )
   W( 47 ) = JVS( 1062 )
   W( 48 ) = JVS( 1063 )
   W( 50 ) = JVS( 1064 )
   W( 51 ) = JVS( 1065 )
   W( 52 ) = JVS( 1066 )
   W( 55 ) = JVS( 1067 )
   W( 56 ) = JVS( 1068 )
   W( 57 ) = JVS( 1069 )
   W( 58 ) = JVS( 1070 )
   W( 59 ) = JVS( 1071 )
   W( 60 ) = JVS( 1072 )
   W( 61 ) = JVS( 1073 )
   W( 62 ) = JVS( 1074 )
   W( 63 ) = JVS( 1075 )
   W( 65 ) = JVS( 1076 )
   W( 66 ) = JVS( 1077 )
   W( 67 ) = JVS( 1078 )
   W( 68 ) = JVS( 1079 )
   W( 69 ) = JVS( 1080 )
   W( 70 ) = JVS( 1081 )
   W( 71 ) = JVS( 1082 )
   W( 72 ) = JVS( 1083 )
   W( 73 ) = JVS( 1084 )
   W( 74 ) = JVS( 1085 )
   W( 75 ) = JVS( 1086 )
   W( 76 ) = JVS( 1087 )
   W( 78 ) = JVS( 1088 )
   W( 79 ) = JVS( 1089 )
   W( 80 ) = JVS( 1090 )
   W( 81 ) = JVS( 1091 )
   W( 82 ) = JVS( 1092 )
   W( 83 ) = JVS( 1093 )
   W( 84 ) = JVS( 1094 )
   W( 85 ) = JVS( 1095 )
   W( 86 ) = JVS( 1096 )
   W( 87 ) = JVS( 1097 )
   W( 88 ) = JVS( 1098 )
   W( 89 ) = JVS( 1099 )
   W( 90 ) = JVS( 1100 )
   W( 91 ) = JVS( 1101 )
   W( 92 ) = JVS( 1102 )
   W( 93 ) = JVS( 1103 )
   W( 94 ) = JVS( 1104 )
   W( 95 ) = JVS( 1105 )
   W( 96 ) = JVS( 1106 )
   W( 97 ) = JVS( 1107 )
   W( 98 ) = JVS( 1108 )
   W( 99 ) = JVS( 1109 )
   W( 100 ) = JVS( 1110 )
   W( 101 ) = JVS( 1111 )
   W( 102 ) = JVS( 1112 )
   W( 103 ) = JVS( 1113 )
   W( 104 ) = JVS( 1114 )
  a = -W( 3 ) / JVS( 8 )
  W( 3 ) = -a
  a = -W( 5 ) / JVS( 10 )
  W( 5 ) = -a
  a = -W( 30 ) / JVS( 203 )
  W( 30 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 204 )
  a = -W( 31 ) / JVS( 205 )
  W( 31 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 206 )
  a = -W( 32 ) / JVS( 208 )
  W( 32 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 209 )
  a = -W( 33 ) / JVS( 210 )
  W( 33 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 211 )
  a = -W( 34 ) / JVS( 213 )
  W( 34 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 214 )
  a = -W( 35 ) / JVS( 215 )
  W( 35 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 216 )
  a = -W( 36 ) / JVS( 217 )
  W( 36 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 218 )
  a = -W( 37 ) / JVS( 219 )
  W( 37 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 220 )
  a = -W( 42 ) / JVS( 233 )
  W( 42 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 234 )
  W( 103 ) = W( 103 ) + a*JVS( 235 )
  a = -W( 43 ) / JVS( 236 )
  W( 43 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 237 )
  a = -W( 44 ) / JVS( 238 )
  W( 44 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 239 )
  a = -W( 47 ) / JVS( 249 )
  W( 47 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 250 )
  W( 103 ) = W( 103 ) + a*JVS( 251 )
  a = -W( 48 ) / JVS( 252 )
  W( 48 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 253 )
  a = -W( 50 ) / JVS( 258 )
  W( 50 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 259 )
  a = -W( 51 ) / JVS( 260 )
  W( 51 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 261 )
  a = -W( 52 ) / JVS( 262 )
  W( 52 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 263 )
  W( 102 ) = W( 102 ) + a*JVS( 264 )
  W( 103 ) = W( 103 ) + a*JVS( 265 )
  a = -W( 55 ) / JVS( 275 )
  W( 55 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 276 )
  W( 101 ) = W( 101 ) + a*JVS( 277 )
  W( 103 ) = W( 103 ) + a*JVS( 278 )
  a = -W( 56 ) / JVS( 279 )
  W( 56 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 280 )
  W( 102 ) = W( 102 ) + a*JVS( 281 )
  W( 103 ) = W( 103 ) + a*JVS( 282 )
  W( 104 ) = W( 104 ) + a*JVS( 283 )
  a = -W( 57 ) / JVS( 284 )
  W( 57 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 285 )
  a = -W( 58 ) / JVS( 288 )
  W( 58 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 289 )
  a = -W( 59 ) / JVS( 292 )
  W( 59 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 293 )
  a = -W( 60 ) / JVS( 296 )
  W( 60 ) = -a
  W( 99 ) = W( 99 ) + a*JVS( 297 )
  W( 103 ) = W( 103 ) + a*JVS( 298 )
  a = -W( 61 ) / JVS( 299 )
  W( 61 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 300 )
  a = -W( 62 ) / JVS( 303 )
  W( 62 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 304 )
  W( 103 ) = W( 103 ) + a*JVS( 305 )
  a = -W( 63 ) / JVS( 306 )
  W( 63 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 307 )
  W( 103 ) = W( 103 ) + a*JVS( 308 )
  a = -W( 65 ) / JVS( 315 )
  W( 65 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 316 )
  W( 96 ) = W( 96 ) + a*JVS( 317 )
  W( 103 ) = W( 103 ) + a*JVS( 318 )
  W( 104 ) = W( 104 ) + a*JVS( 319 )
  a = -W( 66 ) / JVS( 322 )
  W( 66 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 323 )
  W( 92 ) = W( 92 ) + a*JVS( 324 )
  W( 99 ) = W( 99 ) + a*JVS( 325 )
  W( 103 ) = W( 103 ) + a*JVS( 326 )
  a = -W( 67 ) / JVS( 328 )
  W( 67 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 329 )
  W( 96 ) = W( 96 ) + a*JVS( 330 )
  W( 99 ) = W( 99 ) + a*JVS( 331 )
  W( 103 ) = W( 103 ) + a*JVS( 332 )
  a = -W( 68 ) / JVS( 338 )
  W( 68 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 339 )
  W( 82 ) = W( 82 ) + a*JVS( 340 )
  W( 85 ) = W( 85 ) + a*JVS( 341 )
  W( 92 ) = W( 92 ) + a*JVS( 342 )
  W( 99 ) = W( 99 ) + a*JVS( 343 )
  W( 103 ) = W( 103 ) + a*JVS( 344 )
  a = -W( 69 ) / JVS( 351 )
  W( 69 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 352 )
  W( 73 ) = W( 73 ) + a*JVS( 353 )
  W( 75 ) = W( 75 ) + a*JVS( 354 )
  W( 76 ) = W( 76 ) + a*JVS( 355 )
  W( 78 ) = W( 78 ) + a*JVS( 356 )
  W( 79 ) = W( 79 ) + a*JVS( 357 )
  W( 80 ) = W( 80 ) + a*JVS( 358 )
  W( 81 ) = W( 81 ) + a*JVS( 359 )
  W( 82 ) = W( 82 ) + a*JVS( 360 )
  W( 83 ) = W( 83 ) + a*JVS( 361 )
  W( 84 ) = W( 84 ) + a*JVS( 362 )
  W( 85 ) = W( 85 ) + a*JVS( 363 )
  W( 86 ) = W( 86 ) + a*JVS( 364 )
  W( 88 ) = W( 88 ) + a*JVS( 365 )
  W( 89 ) = W( 89 ) + a*JVS( 366 )
  W( 92 ) = W( 92 ) + a*JVS( 367 )
  W( 99 ) = W( 99 ) + a*JVS( 368 )
  W( 103 ) = W( 103 ) + a*JVS( 369 )
  a = -W( 70 ) / JVS( 376 )
  W( 70 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 377 )
  W( 78 ) = W( 78 ) + a*JVS( 378 )
  W( 81 ) = W( 81 ) + a*JVS( 379 )
  W( 82 ) = W( 82 ) + a*JVS( 380 )
  W( 83 ) = W( 83 ) + a*JVS( 381 )
  W( 84 ) = W( 84 ) + a*JVS( 382 )
  W( 85 ) = W( 85 ) + a*JVS( 383 )
  W( 86 ) = W( 86 ) + a*JVS( 384 )
  W( 89 ) = W( 89 ) + a*JVS( 385 )
  W( 92 ) = W( 92 ) + a*JVS( 386 )
  W( 96 ) = W( 96 ) + a*JVS( 387 )
  W( 99 ) = W( 99 ) + a*JVS( 388 )
  W( 101 ) = W( 101 ) + a*JVS( 389 )
  W( 103 ) = W( 103 ) + a*JVS( 390 )
  a = -W( 71 ) / JVS( 391 )
  W( 71 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 392 )
  W( 92 ) = W( 92 ) + a*JVS( 393 )
  W( 99 ) = W( 99 ) + a*JVS( 394 )
  W( 103 ) = W( 103 ) + a*JVS( 395 )
  a = -W( 72 ) / JVS( 401 )
  W( 72 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 402 )
  W( 79 ) = W( 79 ) + a*JVS( 403 )
  W( 80 ) = W( 80 ) + a*JVS( 404 )
  W( 84 ) = W( 84 ) + a*JVS( 405 )
  W( 87 ) = W( 87 ) + a*JVS( 406 )
  W( 92 ) = W( 92 ) + a*JVS( 407 )
  W( 99 ) = W( 99 ) + a*JVS( 408 )
  W( 101 ) = W( 101 ) + a*JVS( 409 )
  W( 103 ) = W( 103 ) + a*JVS( 410 )
  a = -W( 73 ) / JVS( 411 )
  W( 73 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 92 ) = W( 92 ) + a*JVS( 413 )
  W( 99 ) = W( 99 ) + a*JVS( 414 )
  W( 103 ) = W( 103 ) + a*JVS( 415 )
  a = -W( 74 ) / JVS( 418 )
  W( 74 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 419 )
  W( 95 ) = W( 95 ) + a*JVS( 420 )
  W( 96 ) = W( 96 ) + a*JVS( 421 )
  W( 97 ) = W( 97 ) + a*JVS( 422 )
  W( 98 ) = W( 98 ) + a*JVS( 423 )
  W( 99 ) = W( 99 ) + a*JVS( 424 )
  W( 100 ) = W( 100 ) + a*JVS( 425 )
  W( 101 ) = W( 101 ) + a*JVS( 426 )
  W( 103 ) = W( 103 ) + a*JVS( 427 )
  a = -W( 75 ) / JVS( 428 )
  W( 75 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 429 )
  W( 92 ) = W( 92 ) + a*JVS( 430 )
  W( 99 ) = W( 99 ) + a*JVS( 431 )
  W( 103 ) = W( 103 ) + a*JVS( 432 )
  a = -W( 76 ) / JVS( 433 )
  W( 76 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 434 )
  W( 92 ) = W( 92 ) + a*JVS( 435 )
  W( 99 ) = W( 99 ) + a*JVS( 436 )
  W( 103 ) = W( 103 ) + a*JVS( 437 )
  a = -W( 78 ) / JVS( 469 )
  W( 78 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 470 )
  W( 88 ) = W( 88 ) + a*JVS( 471 )
  W( 92 ) = W( 92 ) + a*JVS( 472 )
  W( 99 ) = W( 99 ) + a*JVS( 473 )
  W( 103 ) = W( 103 ) + a*JVS( 474 )
  a = -W( 79 ) / JVS( 475 )
  W( 79 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 92 ) = W( 92 ) + a*JVS( 477 )
  W( 99 ) = W( 99 ) + a*JVS( 478 )
  W( 103 ) = W( 103 ) + a*JVS( 479 )
  a = -W( 80 ) / JVS( 480 )
  W( 80 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 92 ) = W( 92 ) + a*JVS( 482 )
  W( 99 ) = W( 99 ) + a*JVS( 483 )
  W( 103 ) = W( 103 ) + a*JVS( 484 )
  a = -W( 81 ) / JVS( 496 )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 497 )
  W( 88 ) = W( 88 ) + a*JVS( 498 )
  W( 92 ) = W( 92 ) + a*JVS( 499 )
  W( 93 ) = W( 93 ) + a*JVS( 500 )
  W( 95 ) = W( 95 ) + a*JVS( 501 )
  W( 96 ) = W( 96 ) + a*JVS( 502 )
  W( 97 ) = W( 97 ) + a*JVS( 503 )
  W( 98 ) = W( 98 ) + a*JVS( 504 )
  W( 99 ) = W( 99 ) + a*JVS( 505 )
  W( 100 ) = W( 100 ) + a*JVS( 506 )
  W( 101 ) = W( 101 ) + a*JVS( 507 )
  W( 103 ) = W( 103 ) + a*JVS( 508 )
  a = -W( 82 ) / JVS( 510 )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 511 )
  W( 88 ) = W( 88 ) + a*JVS( 512 )
  W( 92 ) = W( 92 ) + a*JVS( 513 )
  W( 99 ) = W( 99 ) + a*JVS( 514 )
  W( 103 ) = W( 103 ) + a*JVS( 515 )
  a = -W( 83 ) / JVS( 525 )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 526 )
  W( 87 ) = W( 87 ) + a*JVS( 527 )
  W( 88 ) = W( 88 ) + a*JVS( 528 )
  W( 89 ) = W( 89 ) + a*JVS( 529 )
  W( 90 ) = W( 90 ) + a*JVS( 530 )
  W( 91 ) = W( 91 ) + a*JVS( 531 )
  W( 92 ) = W( 92 ) + a*JVS( 532 )
  W( 93 ) = W( 93 ) + a*JVS( 533 )
  W( 95 ) = W( 95 ) + a*JVS( 534 )
  W( 97 ) = W( 97 ) + a*JVS( 535 )
  W( 98 ) = W( 98 ) + a*JVS( 536 )
  W( 99 ) = W( 99 ) + a*JVS( 537 )
  W( 100 ) = W( 100 ) + a*JVS( 538 )
  W( 103 ) = W( 103 ) + a*JVS( 539 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  a = -W( 85 ) / JVS( 547 )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 548 )
  W( 92 ) = W( 92 ) + a*JVS( 549 )
  W( 99 ) = W( 99 ) + a*JVS( 550 )
  W( 103 ) = W( 103 ) + a*JVS( 551 )
  a = -W( 86 ) / JVS( 572 )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 573 )
  W( 88 ) = W( 88 ) + a*JVS( 574 )
  W( 90 ) = W( 90 ) + a*JVS( 575 )
  W( 91 ) = W( 91 ) + a*JVS( 576 )
  W( 92 ) = W( 92 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 94 ) = W( 94 ) + a*JVS( 579 )
  W( 95 ) = W( 95 ) + a*JVS( 580 )
  W( 96 ) = W( 96 ) + a*JVS( 581 )
  W( 97 ) = W( 97 ) + a*JVS( 582 )
  W( 98 ) = W( 98 ) + a*JVS( 583 )
  W( 99 ) = W( 99 ) + a*JVS( 584 )
  W( 100 ) = W( 100 ) + a*JVS( 585 )
  W( 101 ) = W( 101 ) + a*JVS( 586 )
  W( 102 ) = W( 102 ) + a*JVS( 587 )
  W( 103 ) = W( 103 ) + a*JVS( 588 )
  W( 104 ) = W( 104 ) + a*JVS( 589 )
  a = -W( 87 ) / JVS( 596 )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 597 )
  W( 92 ) = W( 92 ) + a*JVS( 598 )
  W( 93 ) = W( 93 ) + a*JVS( 599 )
  W( 99 ) = W( 99 ) + a*JVS( 600 )
  W( 101 ) = W( 101 ) + a*JVS( 601 )
  W( 103 ) = W( 103 ) + a*JVS( 602 )
  W( 104 ) = W( 104 ) + a*JVS( 603 )
  a = -W( 88 ) / JVS( 614 )
  W( 88 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 615 )
  W( 93 ) = W( 93 ) + a*JVS( 616 )
  W( 99 ) = W( 99 ) + a*JVS( 617 )
  W( 101 ) = W( 101 ) + a*JVS( 618 )
  W( 103 ) = W( 103 ) + a*JVS( 619 )
  a = -W( 89 ) / JVS( 639 )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 640 )
  W( 91 ) = W( 91 ) + a*JVS( 641 )
  W( 92 ) = W( 92 ) + a*JVS( 642 )
  W( 93 ) = W( 93 ) + a*JVS( 643 )
  W( 94 ) = W( 94 ) + a*JVS( 644 )
  W( 96 ) = W( 96 ) + a*JVS( 645 )
  W( 99 ) = W( 99 ) + a*JVS( 646 )
  W( 101 ) = W( 101 ) + a*JVS( 647 )
  W( 103 ) = W( 103 ) + a*JVS( 648 )
  W( 104 ) = W( 104 ) + a*JVS( 649 )
  a = -W( 90 ) / JVS( 661 )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 662 )
  W( 92 ) = W( 92 ) + a*JVS( 663 )
  W( 93 ) = W( 93 ) + a*JVS( 664 )
  W( 94 ) = W( 94 ) + a*JVS( 665 )
  W( 99 ) = W( 99 ) + a*JVS( 666 )
  W( 101 ) = W( 101 ) + a*JVS( 667 )
  W( 102 ) = W( 102 ) + a*JVS( 668 )
  W( 103 ) = W( 103 ) + a*JVS( 669 )
  W( 104 ) = W( 104 ) + a*JVS( 670 )
  a = -W( 91 ) / JVS( 683 )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 684 )
  W( 93 ) = W( 93 ) + a*JVS( 685 )
  W( 94 ) = W( 94 ) + a*JVS( 686 )
  W( 95 ) = W( 95 ) + a*JVS( 687 )
  W( 97 ) = W( 97 ) + a*JVS( 688 )
  W( 99 ) = W( 99 ) + a*JVS( 689 )
  W( 100 ) = W( 100 ) + a*JVS( 690 )
  W( 101 ) = W( 101 ) + a*JVS( 691 )
  W( 102 ) = W( 102 ) + a*JVS( 692 )
  W( 103 ) = W( 103 ) + a*JVS( 693 )
  W( 104 ) = W( 104 ) + a*JVS( 694 )
  a = -W( 92 ) / JVS( 708 )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 709 )
  W( 95 ) = W( 95 ) + a*JVS( 710 )
  W( 96 ) = W( 96 ) + a*JVS( 711 )
  W( 97 ) = W( 97 ) + a*JVS( 712 )
  W( 98 ) = W( 98 ) + a*JVS( 713 )
  W( 99 ) = W( 99 ) + a*JVS( 714 )
  W( 100 ) = W( 100 ) + a*JVS( 715 )
  W( 101 ) = W( 101 ) + a*JVS( 716 )
  W( 103 ) = W( 103 ) + a*JVS( 717 )
  a = -W( 93 ) / JVS( 731 )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 732 )
  W( 95 ) = W( 95 ) + a*JVS( 733 )
  W( 96 ) = W( 96 ) + a*JVS( 734 )
  W( 97 ) = W( 97 ) + a*JVS( 735 )
  W( 98 ) = W( 98 ) + a*JVS( 736 )
  W( 99 ) = W( 99 ) + a*JVS( 737 )
  W( 100 ) = W( 100 ) + a*JVS( 738 )
  W( 101 ) = W( 101 ) + a*JVS( 739 )
  W( 102 ) = W( 102 ) + a*JVS( 740 )
  W( 103 ) = W( 103 ) + a*JVS( 741 )
  W( 104 ) = W( 104 ) + a*JVS( 742 )
  a = -W( 94 ) / JVS( 776 )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 777 )
  W( 96 ) = W( 96 ) + a*JVS( 778 )
  W( 97 ) = W( 97 ) + a*JVS( 779 )
  W( 98 ) = W( 98 ) + a*JVS( 780 )
  W( 99 ) = W( 99 ) + a*JVS( 781 )
  W( 100 ) = W( 100 ) + a*JVS( 782 )
  W( 101 ) = W( 101 ) + a*JVS( 783 )
  W( 102 ) = W( 102 ) + a*JVS( 784 )
  W( 103 ) = W( 103 ) + a*JVS( 785 )
  W( 104 ) = W( 104 ) + a*JVS( 786 )
  a = -W( 95 ) / JVS( 811 )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 812 )
  W( 97 ) = W( 97 ) + a*JVS( 813 )
  W( 98 ) = W( 98 ) + a*JVS( 814 )
  W( 99 ) = W( 99 ) + a*JVS( 815 )
  W( 100 ) = W( 100 ) + a*JVS( 816 )
  W( 101 ) = W( 101 ) + a*JVS( 817 )
  W( 102 ) = W( 102 ) + a*JVS( 818 )
  W( 103 ) = W( 103 ) + a*JVS( 819 )
  W( 104 ) = W( 104 ) + a*JVS( 820 )
  a = -W( 96 ) / JVS( 863 )
  W( 96 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 864 )
  W( 98 ) = W( 98 ) + a*JVS( 865 )
  W( 99 ) = W( 99 ) + a*JVS( 866 )
  W( 100 ) = W( 100 ) + a*JVS( 867 )
  W( 101 ) = W( 101 ) + a*JVS( 868 )
  W( 102 ) = W( 102 ) + a*JVS( 869 )
  W( 103 ) = W( 103 ) + a*JVS( 870 )
  W( 104 ) = W( 104 ) + a*JVS( 871 )
  a = -W( 97 ) / JVS( 889 )
  W( 97 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 890 )
  W( 99 ) = W( 99 ) + a*JVS( 891 )
  W( 100 ) = W( 100 ) + a*JVS( 892 )
  W( 101 ) = W( 101 ) + a*JVS( 893 )
  W( 102 ) = W( 102 ) + a*JVS( 894 )
  W( 103 ) = W( 103 ) + a*JVS( 895 )
  W( 104 ) = W( 104 ) + a*JVS( 896 )
  a = -W( 98 ) / JVS( 910 )
  W( 98 ) = -a
  W( 99 ) = W( 99 ) + a*JVS( 911 )
  W( 100 ) = W( 100 ) + a*JVS( 912 )
  W( 101 ) = W( 101 ) + a*JVS( 913 )
  W( 102 ) = W( 102 ) + a*JVS( 914 )
  W( 103 ) = W( 103 ) + a*JVS( 915 )
  W( 104 ) = W( 104 ) + a*JVS( 916 )
  a = -W( 99 ) / JVS( 952 )
  W( 99 ) = -a
  W( 100 ) = W( 100 ) + a*JVS( 953 )
  W( 101 ) = W( 101 ) + a*JVS( 954 )
  W( 102 ) = W( 102 ) + a*JVS( 955 )
  W( 103 ) = W( 103 ) + a*JVS( 956 )
  W( 104 ) = W( 104 ) + a*JVS( 957 )
  a = -W( 100 ) / JVS( 970 )
  W( 100 ) = -a
  W( 101 ) = W( 101 ) + a*JVS( 971 )
  W( 102 ) = W( 102 ) + a*JVS( 972 )
  W( 103 ) = W( 103 ) + a*JVS( 973 )
  W( 104 ) = W( 104 ) + a*JVS( 974 )
  a = -W( 101 ) / JVS( 1013 )
  W( 101 ) = -a
  W( 102 ) = W( 102 ) + a*JVS( 1014 )
  W( 103 ) = W( 103 ) + a*JVS( 1015 )
  W( 104 ) = W( 104 ) + a*JVS( 1016 )
  a = -W( 102 ) / JVS( 1046 )
  W( 102 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 1047 )
  W( 104 ) = W( 104 ) + a*JVS( 1048 )
  JVS( 1049) = W( 3 )
  JVS( 1050) = W( 5 )
  JVS( 1051) = W( 30 )
  JVS( 1052) = W( 31 )
  JVS( 1053) = W( 32 )
  JVS( 1054) = W( 33 )
  JVS( 1055) = W( 34 )
  JVS( 1056) = W( 35 )
  JVS( 1057) = W( 36 )
  JVS( 1058) = W( 37 )
  JVS( 1059) = W( 42 )
  JVS( 1060) = W( 43 )
  JVS( 1061) = W( 44 )
  JVS( 1062) = W( 47 )
  JVS( 1063) = W( 48 )
  JVS( 1064) = W( 50 )
  JVS( 1065) = W( 51 )
  JVS( 1066) = W( 52 )
  JVS( 1067) = W( 55 )
  JVS( 1068) = W( 56 )
  JVS( 1069) = W( 57 )
  JVS( 1070) = W( 58 )
  JVS( 1071) = W( 59 )
  JVS( 1072) = W( 60 )
  JVS( 1073) = W( 61 )
  JVS( 1074) = W( 62 )
  JVS( 1075) = W( 63 )
  JVS( 1076) = W( 65 )
  JVS( 1077) = W( 66 )
  JVS( 1078) = W( 67 )
  JVS( 1079) = W( 68 )
  JVS( 1080) = W( 69 )
  JVS( 1081) = W( 70 )
  JVS( 1082) = W( 71 )
  JVS( 1083) = W( 72 )
  JVS( 1084) = W( 73 )
  JVS( 1085) = W( 74 )
  JVS( 1086) = W( 75 )
  JVS( 1087) = W( 76 )
  JVS( 1088) = W( 78 )
  JVS( 1089) = W( 79 )
  JVS( 1090) = W( 80 )
  JVS( 1091) = W( 81 )
  JVS( 1092) = W( 82 )
  JVS( 1093) = W( 83 )
  JVS( 1094) = W( 84 )
  JVS( 1095) = W( 85 )
  JVS( 1096) = W( 86 )
  JVS( 1097) = W( 87 )
  JVS( 1098) = W( 88 )
  JVS( 1099) = W( 89 )
  JVS( 1100) = W( 90 )
  JVS( 1101) = W( 91 )
  JVS( 1102) = W( 92 )
  JVS( 1103) = W( 93 )
  JVS( 1104) = W( 94 )
  JVS( 1105) = W( 95 )
  JVS( 1106) = W( 96 )
  JVS( 1107) = W( 97 )
  JVS( 1108) = W( 98 )
  JVS( 1109) = W( 99 )
  JVS( 1110) = W( 100 )
  JVS( 1111) = W( 101 )
  JVS( 1112) = W( 102 )
  JVS( 1113) = W( 103 )
  JVS( 1114) = W( 104 )
  IF ( ABS( JVS( 1146 )) < TINY(a) ) THEN
         IER = 104
         RETURN
  END IF
   W( 43 ) = JVS( 1115 )
   W( 48 ) = JVS( 1116 )
   W( 50 ) = JVS( 1117 )
   W( 51 ) = JVS( 1118 )
   W( 57 ) = JVS( 1119 )
   W( 61 ) = JVS( 1120 )
   W( 73 ) = JVS( 1121 )
   W( 75 ) = JVS( 1122 )
   W( 76 ) = JVS( 1123 )
   W( 79 ) = JVS( 1124 )
   W( 80 ) = JVS( 1125 )
   W( 82 ) = JVS( 1126 )
   W( 84 ) = JVS( 1127 )
   W( 85 ) = JVS( 1128 )
   W( 87 ) = JVS( 1129 )
   W( 88 ) = JVS( 1130 )
   W( 89 ) = JVS( 1131 )
   W( 90 ) = JVS( 1132 )
   W( 91 ) = JVS( 1133 )
   W( 92 ) = JVS( 1134 )
   W( 93 ) = JVS( 1135 )
   W( 94 ) = JVS( 1136 )
   W( 95 ) = JVS( 1137 )
   W( 96 ) = JVS( 1138 )
   W( 97 ) = JVS( 1139 )
   W( 98 ) = JVS( 1140 )
   W( 99 ) = JVS( 1141 )
   W( 100 ) = JVS( 1142 )
   W( 101 ) = JVS( 1143 )
   W( 102 ) = JVS( 1144 )
   W( 103 ) = JVS( 1145 )
   W( 104 ) = JVS( 1146 )
  a = -W( 43 ) / JVS( 236 )
  W( 43 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 237 )
  a = -W( 48 ) / JVS( 252 )
  W( 48 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 253 )
  a = -W( 50 ) / JVS( 258 )
  W( 50 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 259 )
  a = -W( 51 ) / JVS( 260 )
  W( 51 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 261 )
  a = -W( 57 ) / JVS( 284 )
  W( 57 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 285 )
  a = -W( 61 ) / JVS( 299 )
  W( 61 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 300 )
  a = -W( 73 ) / JVS( 411 )
  W( 73 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 92 ) = W( 92 ) + a*JVS( 413 )
  W( 99 ) = W( 99 ) + a*JVS( 414 )
  W( 103 ) = W( 103 ) + a*JVS( 415 )
  a = -W( 75 ) / JVS( 428 )
  W( 75 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 429 )
  W( 92 ) = W( 92 ) + a*JVS( 430 )
  W( 99 ) = W( 99 ) + a*JVS( 431 )
  W( 103 ) = W( 103 ) + a*JVS( 432 )
  a = -W( 76 ) / JVS( 433 )
  W( 76 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 434 )
  W( 92 ) = W( 92 ) + a*JVS( 435 )
  W( 99 ) = W( 99 ) + a*JVS( 436 )
  W( 103 ) = W( 103 ) + a*JVS( 437 )
  a = -W( 79 ) / JVS( 475 )
  W( 79 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 92 ) = W( 92 ) + a*JVS( 477 )
  W( 99 ) = W( 99 ) + a*JVS( 478 )
  W( 103 ) = W( 103 ) + a*JVS( 479 )
  a = -W( 80 ) / JVS( 480 )
  W( 80 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 92 ) = W( 92 ) + a*JVS( 482 )
  W( 99 ) = W( 99 ) + a*JVS( 483 )
  W( 103 ) = W( 103 ) + a*JVS( 484 )
  a = -W( 82 ) / JVS( 510 )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 511 )
  W( 88 ) = W( 88 ) + a*JVS( 512 )
  W( 92 ) = W( 92 ) + a*JVS( 513 )
  W( 99 ) = W( 99 ) + a*JVS( 514 )
  W( 103 ) = W( 103 ) + a*JVS( 515 )
  a = -W( 84 ) / JVS( 540 )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 541 )
  W( 92 ) = W( 92 ) + a*JVS( 542 )
  W( 99 ) = W( 99 ) + a*JVS( 543 )
  W( 103 ) = W( 103 ) + a*JVS( 544 )
  a = -W( 85 ) / JVS( 547 )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 548 )
  W( 92 ) = W( 92 ) + a*JVS( 549 )
  W( 99 ) = W( 99 ) + a*JVS( 550 )
  W( 103 ) = W( 103 ) + a*JVS( 551 )
  a = -W( 87 ) / JVS( 596 )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 597 )
  W( 92 ) = W( 92 ) + a*JVS( 598 )
  W( 93 ) = W( 93 ) + a*JVS( 599 )
  W( 99 ) = W( 99 ) + a*JVS( 600 )
  W( 101 ) = W( 101 ) + a*JVS( 601 )
  W( 103 ) = W( 103 ) + a*JVS( 602 )
  W( 104 ) = W( 104 ) + a*JVS( 603 )
  a = -W( 88 ) / JVS( 614 )
  W( 88 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 615 )
  W( 93 ) = W( 93 ) + a*JVS( 616 )
  W( 99 ) = W( 99 ) + a*JVS( 617 )
  W( 101 ) = W( 101 ) + a*JVS( 618 )
  W( 103 ) = W( 103 ) + a*JVS( 619 )
  a = -W( 89 ) / JVS( 639 )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 640 )
  W( 91 ) = W( 91 ) + a*JVS( 641 )
  W( 92 ) = W( 92 ) + a*JVS( 642 )
  W( 93 ) = W( 93 ) + a*JVS( 643 )
  W( 94 ) = W( 94 ) + a*JVS( 644 )
  W( 96 ) = W( 96 ) + a*JVS( 645 )
  W( 99 ) = W( 99 ) + a*JVS( 646 )
  W( 101 ) = W( 101 ) + a*JVS( 647 )
  W( 103 ) = W( 103 ) + a*JVS( 648 )
  W( 104 ) = W( 104 ) + a*JVS( 649 )
  a = -W( 90 ) / JVS( 661 )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 662 )
  W( 92 ) = W( 92 ) + a*JVS( 663 )
  W( 93 ) = W( 93 ) + a*JVS( 664 )
  W( 94 ) = W( 94 ) + a*JVS( 665 )
  W( 99 ) = W( 99 ) + a*JVS( 666 )
  W( 101 ) = W( 101 ) + a*JVS( 667 )
  W( 102 ) = W( 102 ) + a*JVS( 668 )
  W( 103 ) = W( 103 ) + a*JVS( 669 )
  W( 104 ) = W( 104 ) + a*JVS( 670 )
  a = -W( 91 ) / JVS( 683 )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 684 )
  W( 93 ) = W( 93 ) + a*JVS( 685 )
  W( 94 ) = W( 94 ) + a*JVS( 686 )
  W( 95 ) = W( 95 ) + a*JVS( 687 )
  W( 97 ) = W( 97 ) + a*JVS( 688 )
  W( 99 ) = W( 99 ) + a*JVS( 689 )
  W( 100 ) = W( 100 ) + a*JVS( 690 )
  W( 101 ) = W( 101 ) + a*JVS( 691 )
  W( 102 ) = W( 102 ) + a*JVS( 692 )
  W( 103 ) = W( 103 ) + a*JVS( 693 )
  W( 104 ) = W( 104 ) + a*JVS( 694 )
  a = -W( 92 ) / JVS( 708 )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 709 )
  W( 95 ) = W( 95 ) + a*JVS( 710 )
  W( 96 ) = W( 96 ) + a*JVS( 711 )
  W( 97 ) = W( 97 ) + a*JVS( 712 )
  W( 98 ) = W( 98 ) + a*JVS( 713 )
  W( 99 ) = W( 99 ) + a*JVS( 714 )
  W( 100 ) = W( 100 ) + a*JVS( 715 )
  W( 101 ) = W( 101 ) + a*JVS( 716 )
  W( 103 ) = W( 103 ) + a*JVS( 717 )
  a = -W( 93 ) / JVS( 731 )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 732 )
  W( 95 ) = W( 95 ) + a*JVS( 733 )
  W( 96 ) = W( 96 ) + a*JVS( 734 )
  W( 97 ) = W( 97 ) + a*JVS( 735 )
  W( 98 ) = W( 98 ) + a*JVS( 736 )
  W( 99 ) = W( 99 ) + a*JVS( 737 )
  W( 100 ) = W( 100 ) + a*JVS( 738 )
  W( 101 ) = W( 101 ) + a*JVS( 739 )
  W( 102 ) = W( 102 ) + a*JVS( 740 )
  W( 103 ) = W( 103 ) + a*JVS( 741 )
  W( 104 ) = W( 104 ) + a*JVS( 742 )
  a = -W( 94 ) / JVS( 776 )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 777 )
  W( 96 ) = W( 96 ) + a*JVS( 778 )
  W( 97 ) = W( 97 ) + a*JVS( 779 )
  W( 98 ) = W( 98 ) + a*JVS( 780 )
  W( 99 ) = W( 99 ) + a*JVS( 781 )
  W( 100 ) = W( 100 ) + a*JVS( 782 )
  W( 101 ) = W( 101 ) + a*JVS( 783 )
  W( 102 ) = W( 102 ) + a*JVS( 784 )
  W( 103 ) = W( 103 ) + a*JVS( 785 )
  W( 104 ) = W( 104 ) + a*JVS( 786 )
  a = -W( 95 ) / JVS( 811 )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 812 )
  W( 97 ) = W( 97 ) + a*JVS( 813 )
  W( 98 ) = W( 98 ) + a*JVS( 814 )
  W( 99 ) = W( 99 ) + a*JVS( 815 )
  W( 100 ) = W( 100 ) + a*JVS( 816 )
  W( 101 ) = W( 101 ) + a*JVS( 817 )
  W( 102 ) = W( 102 ) + a*JVS( 818 )
  W( 103 ) = W( 103 ) + a*JVS( 819 )
  W( 104 ) = W( 104 ) + a*JVS( 820 )
  a = -W( 96 ) / JVS( 863 )
  W( 96 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 864 )
  W( 98 ) = W( 98 ) + a*JVS( 865 )
  W( 99 ) = W( 99 ) + a*JVS( 866 )
  W( 100 ) = W( 100 ) + a*JVS( 867 )
  W( 101 ) = W( 101 ) + a*JVS( 868 )
  W( 102 ) = W( 102 ) + a*JVS( 869 )
  W( 103 ) = W( 103 ) + a*JVS( 870 )
  W( 104 ) = W( 104 ) + a*JVS( 871 )
  a = -W( 97 ) / JVS( 889 )
  W( 97 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 890 )
  W( 99 ) = W( 99 ) + a*JVS( 891 )
  W( 100 ) = W( 100 ) + a*JVS( 892 )
  W( 101 ) = W( 101 ) + a*JVS( 893 )
  W( 102 ) = W( 102 ) + a*JVS( 894 )
  W( 103 ) = W( 103 ) + a*JVS( 895 )
  W( 104 ) = W( 104 ) + a*JVS( 896 )
  a = -W( 98 ) / JVS( 910 )
  W( 98 ) = -a
  W( 99 ) = W( 99 ) + a*JVS( 911 )
  W( 100 ) = W( 100 ) + a*JVS( 912 )
  W( 101 ) = W( 101 ) + a*JVS( 913 )
  W( 102 ) = W( 102 ) + a*JVS( 914 )
  W( 103 ) = W( 103 ) + a*JVS( 915 )
  W( 104 ) = W( 104 ) + a*JVS( 916 )
  a = -W( 99 ) / JVS( 952 )
  W( 99 ) = -a
  W( 100 ) = W( 100 ) + a*JVS( 953 )
  W( 101 ) = W( 101 ) + a*JVS( 954 )
  W( 102 ) = W( 102 ) + a*JVS( 955 )
  W( 103 ) = W( 103 ) + a*JVS( 956 )
  W( 104 ) = W( 104 ) + a*JVS( 957 )
  a = -W( 100 ) / JVS( 970 )
  W( 100 ) = -a
  W( 101 ) = W( 101 ) + a*JVS( 971 )
  W( 102 ) = W( 102 ) + a*JVS( 972 )
  W( 103 ) = W( 103 ) + a*JVS( 973 )
  W( 104 ) = W( 104 ) + a*JVS( 974 )
  a = -W( 101 ) / JVS( 1013 )
  W( 101 ) = -a
  W( 102 ) = W( 102 ) + a*JVS( 1014 )
  W( 103 ) = W( 103 ) + a*JVS( 1015 )
  W( 104 ) = W( 104 ) + a*JVS( 1016 )
  a = -W( 102 ) / JVS( 1046 )
  W( 102 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 1047 )
  W( 104 ) = W( 104 ) + a*JVS( 1048 )
  a = -W( 103 ) / JVS( 1113 )
  W( 103 ) = -a
  W( 104 ) = W( 104 ) + a*JVS( 1114 )
  JVS( 1115) = W( 43 )
  JVS( 1116) = W( 48 )
  JVS( 1117) = W( 50 )
  JVS( 1118) = W( 51 )
  JVS( 1119) = W( 57 )
  JVS( 1120) = W( 61 )
  JVS( 1121) = W( 73 )
  JVS( 1122) = W( 75 )
  JVS( 1123) = W( 76 )
  JVS( 1124) = W( 79 )
  JVS( 1125) = W( 80 )
  JVS( 1126) = W( 82 )
  JVS( 1127) = W( 84 )
  JVS( 1128) = W( 85 )
  JVS( 1129) = W( 87 )
  JVS( 1130) = W( 88 )
  JVS( 1131) = W( 89 )
  JVS( 1132) = W( 90 )
  JVS( 1133) = W( 91 )
  JVS( 1134) = W( 92 )
  JVS( 1135) = W( 93 )
  JVS( 1136) = W( 94 )
  JVS( 1137) = W( 95 )
  JVS( 1138) = W( 96 )
  JVS( 1139) = W( 97 )
  JVS( 1140) = W( 98 )
  JVS( 1141) = W( 99 )
  JVS( 1142) = W( 100 )
  JVS( 1143) = W( 101 )
  JVS( 1144) = W( 102 )
  JVS( 1145) = W( 103 )
  JVS( 1146) = W( 104 )
   END SUBROUTINE decomp_saprc99_mosaic_8bin_vbs2
END MODULE saprc99_mosaic_8bin_vbs2_Integrator
