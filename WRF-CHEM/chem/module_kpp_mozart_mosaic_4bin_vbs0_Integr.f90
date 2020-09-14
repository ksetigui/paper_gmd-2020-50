MODULE mozart_mosaic_4bin_vbs0_Integrator
 USE mozart_mosaic_4bin_vbs0_Parameters
 USE mozart_mosaic_4bin_vbs0_Precision
 USE mozart_mosaic_4bin_vbs0_JacobianSP
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
SUBROUTINE mozart_mosaic_4bin_vbs0_INTEGRATE( TIN, TOUT, &
  FIX, VAR, RCONST, ATOL, RTOL, &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )
   USE mozart_mosaic_4bin_vbs0_Parameters
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
   CALL mozart_mosaic_4bin_vbs0_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT, &
         ATOL,RTOL, &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
   STEPMIN = RCNTRL(ihexit)
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U)) IERR_U = IERR
END SUBROUTINE mozart_mosaic_4bin_vbs0_INTEGRATE
SUBROUTINE mozart_mosaic_4bin_vbs0_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol, &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IERR)
  USE mozart_mosaic_4bin_vbs0_Parameters
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
      CALL mozart_mosaic_4bin_vbs0_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL mozart_mosaic_4bin_vbs0_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF
   Roundoff = mozart_mosaic_4bin_vbs0_WLAMCH('E')
   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL mozart_mosaic_4bin_vbs0_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL mozart_mosaic_4bin_vbs0_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL mozart_mosaic_4bin_vbs0_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL mozart_mosaic_4bin_vbs0_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL mozart_mosaic_4bin_vbs0_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL mozart_mosaic_4bin_vbs0_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL mozart_mosaic_4bin_vbs0_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL mozart_mosaic_4bin_vbs0_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO
   SELECT CASE (Method)
     CASE (1)
       CALL mozart_mosaic_4bin_vbs0_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL mozart_mosaic_4bin_vbs0_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL mozart_mosaic_4bin_vbs0_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL mozart_mosaic_4bin_vbs0_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL mozart_mosaic_4bin_vbs0_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL mozart_mosaic_4bin_vbs0_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT
   CALL mozart_mosaic_4bin_vbs0_ros_Integrator(Y,Tstart,Tend,Texit, &
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
 SUBROUTINE mozart_mosaic_4bin_vbs0_ros_ErrorMsg(Code,T,H,IERR)
   USE mozart_mosaic_4bin_vbs0_Precision
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
 END SUBROUTINE mozart_mosaic_4bin_vbs0_ros_ErrorMsg
 SUBROUTINE mozart_mosaic_4bin_vbs0_ros_Integrator (Y, Tstart, Tend, T, &
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
      CALL mozart_mosaic_4bin_vbs0_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN
      CALL mozart_mosaic_4bin_vbs0_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF
   Hexit = H
   H = MIN(H,ABS(Tend-T))
   CALL mozart_mosaic_4bin_vbs0_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF (.NOT.Autonomous) THEN
      CALL mozart_mosaic_4bin_vbs0_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF
   CALL mozart_mosaic_4bin_vbs0_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)
UntilAccepted: DO
   CALL mozart_mosaic_4bin_vbs0_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec, Nsng )
   IF (Singular) THEN
       CALL mozart_mosaic_4bin_vbs0_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF
Stage: DO istage = 1, ros_S
       ioffset = NVAR*(istage-1)
       IF ( istage == 1 ) THEN
         CALL mozart_mosaic_4bin_vbs0_WCOPY(NVAR,Fcn0,1,Fcn,1)
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL mozart_mosaic_4bin_vbs0_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL mozart_mosaic_4bin_vbs0_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL mozart_mosaic_4bin_vbs0_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF
       CALL mozart_mosaic_4bin_vbs0_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL mozart_mosaic_4bin_vbs0_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL mozart_mosaic_4bin_vbs0_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL mozart_mosaic_4bin_vbs0_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)
   END DO Stage
   CALL mozart_mosaic_4bin_vbs0_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL mozart_mosaic_4bin_vbs0_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO
   CALL mozart_mosaic_4bin_vbs0_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL mozart_mosaic_4bin_vbs0_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = mozart_mosaic_4bin_vbs0_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )
   Fac = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac
   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN
      Nacc = Nacc+1
      CALL mozart_mosaic_4bin_vbs0_WCOPY(NVAR,Ynew,1,Y,1)
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
  END SUBROUTINE mozart_mosaic_4bin_vbs0_ros_Integrator
  REAL(kind=dp) FUNCTION mozart_mosaic_4bin_vbs0_ros_ErrorNorm ( Y, Ynew, Yerr, &
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
    mozart_mosaic_4bin_vbs0_ros_ErrorNorm = Err
  END FUNCTION mozart_mosaic_4bin_vbs0_ros_ErrorNorm
  SUBROUTINE mozart_mosaic_4bin_vbs0_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   IMPLICIT NONE
   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)
   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)
   INTEGER, INTENT(INOUT) ::Nfun
   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp
   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL mozart_mosaic_4bin_vbs0_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL mozart_mosaic_4bin_vbs0_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL mozart_mosaic_4bin_vbs0_WSCAL(NVAR,(ONE/Delta),dFdT,1)
  END SUBROUTINE mozart_mosaic_4bin_vbs0_ros_FunTimeDeriv
  SUBROUTINE mozart_mosaic_4bin_vbs0_ros_PrepareMatrix ( H, Direction, gam, &
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
     CALL mozart_mosaic_4bin_vbs0_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL mozart_mosaic_4bin_vbs0_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO
     CALL mozart_mosaic_4bin_vbs0_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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
  END SUBROUTINE mozart_mosaic_4bin_vbs0_ros_PrepareMatrix
  SUBROUTINE mozart_mosaic_4bin_vbs0_ros_Decomp( A, Pivot, ising, Ndec )
   IMPLICIT NONE
   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)
   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec
CALL decomp_mozart_mosaic_4bin_vbs0 ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1
  END SUBROUTINE mozart_mosaic_4bin_vbs0_ros_Decomp
  SUBROUTINE mozart_mosaic_4bin_vbs0_ros_Solve( A, Pivot, b, Nsol )
   IMPLICIT NONE
   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)
   INTEGER, INTENT(INOUT) :: nsol
   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)
   CALL mozart_mosaic_4bin_vbs0_KppSolve( A, b )
   Nsol = Nsol+1
  END SUBROUTINE mozart_mosaic_4bin_vbs0_ros_Solve
  SUBROUTINE mozart_mosaic_4bin_vbs0_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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
 END SUBROUTINE mozart_mosaic_4bin_vbs0_Ros2
  SUBROUTINE mozart_mosaic_4bin_vbs0_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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
  END SUBROUTINE mozart_mosaic_4bin_vbs0_Ros3
  SUBROUTINE mozart_mosaic_4bin_vbs0_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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
  END SUBROUTINE mozart_mosaic_4bin_vbs0_Ros4
  SUBROUTINE mozart_mosaic_4bin_vbs0_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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
  END SUBROUTINE mozart_mosaic_4bin_vbs0_Rodas3
  SUBROUTINE mozart_mosaic_4bin_vbs0_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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
  END SUBROUTINE mozart_mosaic_4bin_vbs0_Rodas4
END SUBROUTINE mozart_mosaic_4bin_vbs0_Rosenbrock
SUBROUTINE mozart_mosaic_4bin_vbs0_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )
   USE mozart_mosaic_4bin_vbs0_Parameters
   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)
   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun
   CALL mozart_mosaic_4bin_vbs0_Fun( Y, FIX, RCONST, Ydot )
   Nfun = Nfun+1
END SUBROUTINE mozart_mosaic_4bin_vbs0_FunTemplate
SUBROUTINE mozart_mosaic_4bin_vbs0_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )
 USE mozart_mosaic_4bin_vbs0_Parameters
 USE mozart_mosaic_4bin_vbs0_Jacobian
    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)
    INTEGER :: Njac
    REAL(kind=dp) :: Jcb(LU_NONZERO)
    REAL(kind=dp) :: Told
    CALL mozart_mosaic_4bin_vbs0_Jac_SP( Y, FIX, RCONST, Jcb )
    Njac = Njac+1
END SUBROUTINE mozart_mosaic_4bin_vbs0_JacTemplate
SUBROUTINE mozart_mosaic_4bin_vbs0_Fun ( V, F, RCT, Vdot )
  REAL(kind=dp) :: V(NVAR)
  REAL(kind=dp) :: F(NFIX)
  REAL(kind=dp) :: RCT(NREACT)
  REAL(kind=dp) :: Vdot(NVAR)
  REAL(kind=dp) :: A(NREACT)
  A(1) = RCT(1)*F(2)
  A(2) = RCT(2)*V(83)
  A(3) = RCT(3)*V(83)
  A(4) = RCT(4)*V(15)
  A(5) = RCT(5)*V(88)
  A(6) = RCT(6)*V(21)
  A(7) = RCT(7)*V(45)
  A(8) = RCT(8)*V(81)
  A(9) = RCT(9)*V(34)
  A(10) = RCT(10)*V(33)
  A(11) = RCT(11)*V(79)
  A(12) = RCT(12)*V(79)
  A(13) = RCT(13)*V(19)
  A(14) = RCT(14)*V(70)
  A(15) = RCT(15)*V(46)
  A(16) = RCT(16)*V(40)
  A(17) = RCT(17)*V(47)
  A(18) = RCT(18)*V(50)
  A(19) = RCT(19)*V(74)
  A(20) = RCT(20)*V(77)
  A(21) = RCT(21)*V(28)
  A(22) = RCT(22)*V(29)
  A(23) = RCT(23)*V(30)
  A(24) = RCT(24)*V(58)
  A(25) = RCT(25)*V(75)
  A(26) = RCT(26)*V(23)
  A(27) = RCT(27)*V(68)
  A(28) = RCT(28)*V(49)
  A(29) = RCT(29)*V(73)
  A(30) = RCT(30)*V(62)
  A(31) = RCT(31)*V(41)
  A(32) = RCT(32)*V(53)
  A(33) = RCT(33)*V(43)
  A(34) = RCT(34)*V(56)
  A(35) = RCT(35)*V(32)
  A(36) = RCT(36)*V(38)
  A(37) = RCT(37)*V(44)
  A(38) = RCT(38)*V(60)*F(2)
  A(39) = RCT(39)*V(60)*V(83)
  A(40) = RCT(40)*V(55)*F(2)
  A(41) = RCT(41)*V(55)*F(1)
  A(42) = RCT(42)*V(27)*V(55)
  A(43) = RCT(43)*V(27)*V(84)
  A(44) = RCT(44)*V(60)*V(84)
  A(45) = RCT(45)*V(60)*V(86)
  A(46) = RCT(46)*V(83)*V(84)
  A(47) = RCT(47)*V(83)*V(86)
  A(48) = RCT(48)*V(86)*V(86)*F(1)
  A(49) = RCT(49)*V(19)*V(84)
  A(50) = RCT(50)*V(84)*V(86)
  A(51) = RCT(51)*V(84)*V(84)
  A(52) = RCT(52)*V(84)*V(84)
  A(53) = RCT(53)*V(15)*V(55)
  A(54) = RCT(54)*V(15)*V(55)
  A(55) = RCT(55)*V(86)*V(87)
  A(56) = RCT(56)*V(83)*V(87)
  A(57) = RCT(57)*V(60)*V(88)
  A(58) = RCT(58)*V(83)*V(88)
  A(59) = RCT(59)*V(81)*V(86)
  A(60) = RCT(60)*V(81)*V(88)
  A(61) = RCT(61)*V(21)
  A(62) = RCT(62)*V(84)*V(88)
  A(63) = RCT(63)*V(45)*V(84)
  A(64) = RCT(64)*V(81)*V(87)
  A(65) = RCT(65)*V(86)*V(88)
  A(66) = RCT(66)*V(34)*V(84)
  A(67) = RCT(67)*V(34)
  A(68) = RCT(68)*V(21)*F(2)
  A(69) = RCT(69)*V(81)
  A(70) = RCT(70)*V(88)
  A(71) = RCT(71)*V(48)*V(84)
  A(72) = RCT(72)*V(48)*V(55)
  A(73) = RCT(73)*V(82)*V(87)
  A(74) = RCT(74)*V(82)*V(82)
  A(75) = RCT(75)*V(82)*V(82)
  A(76) = RCT(76)*V(82)*V(86)
  A(77) = RCT(77)*V(33)*V(84)
  A(78) = RCT(78)*V(79)*V(81)
  A(79) = RCT(79)*V(79)*V(84)
  A(80) = RCT(80)*V(65)*V(84)
  A(81) = RCT(81)*V(36)*V(84)
  A(82) = RCT(82)*V(36)*V(83)
  A(83) = RCT(83)*V(16)*V(84)
  A(84) = RCT(84)*V(43)*V(84)
  A(85) = RCT(85)*V(39)*V(87)
  A(86) = RCT(86)*V(20)*F(2)
  A(87) = RCT(87)*V(20)
  A(88) = RCT(88)*V(11)*V(84)
  A(89) = RCT(89)*V(64)*V(87)
  A(90) = RCT(90)*V(64)*V(86)
  A(91) = RCT(91)*V(64)*V(82)
  A(92) = RCT(92)*V(28)*V(84)
  A(93) = RCT(93)*V(69)*V(84)
  A(94) = RCT(94)*V(69)*V(83)
  A(95) = RCT(95)*V(69)*V(81)
  A(96) = RCT(96)*V(59)*V(87)
  A(97) = RCT(97)*V(59)*V(86)
  A(98) = RCT(98)*V(46)*V(84)
  A(99) = RCT(99)*V(70)*V(84)
  A(100) = RCT(100)*V(70)*V(81)
  A(101) = RCT(101)*V(87)*V(89)
  A(102) = RCT(102)*V(88)*V(89)
  A(103) = RCT(103)*V(86)*V(89)
  A(104) = RCT(104)*V(82)*V(89)
  A(105) = RCT(105)*V(40)*V(84)
  A(106) = RCT(106)*V(47)
  A(107) = RCT(107)*V(89)*V(89)
  A(108) = RCT(108)*V(12)*V(84)
  A(109) = RCT(109)*V(66)*V(87)
  A(110) = RCT(110)*V(66)*V(86)
  A(111) = RCT(111)*V(66)*V(82)
  A(112) = RCT(112)*V(29)*V(84)
  A(113) = RCT(113)*V(58)*V(84)
  A(114) = RCT(114)*V(72)*V(87)
  A(115) = RCT(115)*V(72)*V(86)
  A(116) = RCT(116)*V(72)*V(82)
  A(117) = RCT(117)*V(30)*V(84)
  A(118) = RCT(118)*V(13)*V(84)
  A(119) = RCT(119)*V(31)*V(87)
  A(120) = RCT(120)*V(14)*V(84)
  A(121) = RCT(121)*V(63)*V(87)
  A(122) = RCT(122)*V(63)*V(86)
  A(123) = RCT(123)*V(56)*V(84)
  A(124) = RCT(124)*V(37)*V(84)
  A(125) = RCT(125)*V(41)*V(84)
  A(126) = RCT(126)*V(57)*V(87)
  A(127) = RCT(127)*V(57)*V(86)
  A(128) = RCT(128)*V(32)*V(84)
  A(129) = RCT(129)*V(17)*V(84)
  A(130) = RCT(130)*V(18)*V(84)
  A(131) = RCT(131)*V(22)*V(88)
  A(132) = RCT(132)*V(52)*V(87)
  A(133) = RCT(133)*V(52)*V(86)
  A(134) = RCT(134)*V(38)*V(84)
  A(135) = RCT(135)*V(67)*V(84)
  A(136) = RCT(136)*V(67)*V(83)
  A(137) = RCT(137)*V(80)*V(87)
  A(138) = RCT(138)*V(80)*V(81)
  A(139) = RCT(139)*V(80)*V(86)
  A(140) = RCT(140)*V(49)*V(84)
  A(141) = RCT(141)*V(80)*V(82)
  A(142) = RCT(142)*V(80)*V(89)
  A(143) = RCT(143)*V(77)*V(84)
  A(144) = RCT(144)*V(77)*V(83)
  A(145) = RCT(145)*V(74)*V(84)
  A(146) = RCT(146)*V(74)*V(83)
  A(147) = RCT(147)*V(78)*V(87)
  A(148) = RCT(148)*V(78)*V(87)
  A(149) = RCT(149)*V(78)*V(81)
  A(150) = RCT(150)*V(78)*V(86)
  A(151) = RCT(151)*V(78)*V(82)
  A(152) = RCT(152)*V(78)*V(89)
  A(153) = RCT(153)*V(26)*V(84)
  A(154) = RCT(154)*V(85)*V(87)
  A(155) = RCT(155)*V(81)*V(85)
  A(156) = RCT(156)*V(85)*V(86)
  A(157) = RCT(157)*V(82)*V(85)
  A(158) = RCT(158)*V(85)*V(89)
  A(159) = RCT(159)*V(85)*V(85)
  A(160) = RCT(160)*V(85)*V(88)*F(2)
  A(161) = RCT(161)*V(50)*F(2)
  A(162) = RCT(162)*V(54)*V(84)
  A(163) = RCT(163)*V(54)*V(83)
  A(164) = RCT(164)*V(54)*V(81)
  A(165) = RCT(165)*V(71)*V(87)
  A(166) = RCT(166)*V(71)*V(86)
  A(167) = RCT(167)*V(44)*V(84)
  A(168) = RCT(168)*V(42)*V(84)
  A(169) = RCT(169)*V(67)*V(81)
  A(170) = RCT(170)*V(61)*V(87)
  A(171) = RCT(171)*V(61)*V(81)
  A(172) = RCT(172)*V(61)*V(86)
  A(173) = RCT(173)*V(75)*V(84)
  A(174) = RCT(174)*V(75)*V(81)
  A(175) = RCT(175)*V(68)*V(84)
  A(176) = RCT(176)*V(68)*V(81)
  A(177) = RCT(177)*V(35)*V(84)
  A(178) = RCT(178)*V(76)*V(87)
  A(179) = RCT(179)*V(76)*V(81)
  A(180) = RCT(180)*V(76)*V(86)
  A(181) = RCT(181)*V(76)*V(82)
  A(182) = RCT(182)*V(76)*V(89)
  A(183) = RCT(183)*V(23)*V(84)
  A(184) = RCT(184)*V(23)*V(84)
  A(185) = RCT(185)*V(51)*V(84)
  A(186) = RCT(186)*V(25)*V(84)
  A(187) = RCT(187)*V(50)*V(84)
  A(188) = RCT(188)*V(47)*V(84)
  A(189) = RCT(189)*V(73)*V(84)
  A(190) = RCT(190)*V(62)*V(84)
  A(191) = RCT(191)*V(24)*V(84)
  A(192) = RCT(192)*V(24)*V(84)
  A(193) = RCT(193)*V(24)*V(81)
  A(194) = RCT(194)*V(10)*V(84)
  A(195) = RCT(195)*V(86)
  A(196) = RCT(196)*V(64)*V(64)
  A(197) = RCT(197)*V(7)*V(84)
  A(198) = RCT(198)*V(9)*V(84)
  A(199) = RCT(199)*V(67)*V(84)
  A(200) = RCT(200)*V(54)*V(84)
  Vdot(1) = A(83)
  Vdot(2) = A(114)
  Vdot(3) = A(115)
  Vdot(4) = A(199)
  Vdot(5) = A(200)
  Vdot(6) = A(197)
  Vdot(7) = -A(197)
  Vdot(8) = A(198)
  Vdot(9) = -A(198)
  Vdot(10) = -A(194)
  Vdot(11) = -A(88)
  Vdot(12) = -A(108)
  Vdot(13) = -A(118)
  Vdot(14) = -A(120)
  Vdot(15) = -A(4)-A(53)-A(54)
  Vdot(16) = -A(83)+A(191)+0.5*A(192)+A(193)
  Vdot(17) = -A(129)
  Vdot(18) = 0.25*A(129)-A(130)
  Vdot(19) = -A(13)+A(48)-A(49)+A(52)+0.5*A(195)
  Vdot(20) = A(85)-A(86)-A(87)
  Vdot(21) = -A(6)+A(60)-A(61)-A(68)
  Vdot(22) = A(130)-A(131)
  Vdot(23) = -A(26)+A(180)-A(183)-A(184)
  Vdot(24) = -A(191)-A(192)-A(193)
  Vdot(25) = 0.2*A(91)-A(186)+0.4*A(196)
  Vdot(26) = A(150)-A(153)
  Vdot(27) = A(12)-A(42)-A(43)+0.05*A(72)
  Vdot(28) = -A(21)+A(90)-A(92)
  Vdot(29) = -A(22)+A(110)-A(112)
  Vdot(30) = -A(23)+A(115)-A(117)
  Vdot(31) = A(118)-A(119)
  Vdot(32) = -A(35)+A(127)-A(128)
  Vdot(33) = -A(10)+A(76)-A(77)
  Vdot(34) = -A(9)+A(65)-A(66)-A(67)
  Vdot(35) = 0.37*A(137)+0.4*A(138)+0.3*A(141)+0.4*A(142)+A(175)+A(176)-A(177)
  Vdot(36) = -A(81)-A(82)
  Vdot(37) = A(95)+0.1*A(121)-A(124)
  Vdot(38) = -A(36)+A(133)-A(134)
  Vdot(39) = 0.75*A(81)-A(85)
  Vdot(40) = -A(16)+0.75*A(103)-A(105)+0.75*A(156)
  Vdot(41) = -A(31)+0.8*A(34)+0.75*A(121)-A(125)
  Vdot(42) = 0.25*A(82)+0.25*A(94)+0.25*A(103)+0.1*A(104)+0.2*A(136)+0.25*A(156)-A(168)
  Vdot(43) = 0.13*A(32)-A(33)+0.45*A(36)-A(84)+0.45*A(132)+0.2*A(190)
  Vdot(44) = -A(37)+A(166)-A(167)
  Vdot(45) = -A(7)+A(62)-A(63)+2*A(68)+A(69)+0.5*A(70)+A(78)+A(100)+A(174)+A(193)
  Vdot(46) = -A(15)+A(97)-A(98)
  Vdot(47) = -A(17)+A(102)-A(106)-A(188)
  Vdot(48) = -A(71)-A(72)+0.08*A(94)
  Vdot(49) = -A(28)+A(139)-A(140)
  Vdot(50) = -A(18)+A(160)-A(161)-A(187)
  Vdot(51) = A(75)+0.3*A(91)+0.5*A(116)+0.25*A(141)+0.25*A(151)+0.3*A(181)-A(185)
  Vdot(52) = 0.7*A(129)-A(132)-A(133)+A(134)
  Vdot(53) = -A(32)+0.9*A(36)+0.7*A(131)+0.9*A(132)
  Vdot(54) = -A(162)-A(163)-A(164)
  Vdot(55) = A(2)+A(4)-A(40)-A(41)-A(42)-A(53)-A(54)-A(72)
  Vdot(56) = -A(34)+A(122)-A(123)
  Vdot(57) = A(125)-A(126)-A(127)+A(128)
  Vdot(58) = 0.82*A(22)-A(24)+0.25*A(34)+0.1*A(37)+0.82*A(109)+0.82*A(111)-A(113)+0.5*A(119)+0.25*A(121)+0.1*A(165)
  Vdot(59) = A(93)-A(96)-A(97)+0.5*A(98)
  Vdot(60) = 2*A(1)+A(3)+A(5)-A(38)-A(39)+A(40)-A(44)-A(45)+A(51)-A(57)
  Vdot(61) = A(169)-A(170)-A(171)-A(172)
  Vdot(62) = -A(30)+A(86)+0.53*A(147)+0.53*A(149)+0.26*A(151)+0.53*A(152)+0.25*A(178)+0.25*A(179)+0.1*A(181)+0.25*A(182)&
               &-A(190)
  Vdot(63) = A(120)-A(121)-A(122)+A(123)
  Vdot(64) = A(31)+A(88)-A(89)-A(90)-A(91)+0.5*A(92)-2*A(196)
  Vdot(65) = A(11)+A(12)+A(14)+0.67*A(19)+0.7*A(20)+A(25)+A(27)+A(30)+0.45*A(32)+2*A(33)+A(78)+A(79)-A(80)+0.5*A(82)&
               &+A(84)+0.56*A(94)+0.3*A(136)+0.05*A(144)+0.2*A(146)+0.22*A(147)+0.22*A(149)+0.11*A(151)+0.22*A(152)+A(173)&
               &+A(174)+A(178)+A(179)+0.4*A(181)+A(182)
  Vdot(66) = A(108)-A(109)-A(110)-A(111)+A(112)
  Vdot(67) = -A(135)-A(136)-A(169)
  Vdot(68) = -A(27)+0.08*A(137)+0.8*A(148)+0.794*A(170)+0.794*A(171)+0.794*A(172)-A(175)-A(176)
  Vdot(69) = 0.7*A(20)-A(93)-A(94)-A(95)+0.07*A(136)
  Vdot(70) = -A(14)+A(15)+A(21)+0.4*A(34)+A(35)+A(89)+0.8*A(91)+0.5*A(92)+0.5*A(94)+A(96)-A(99)-A(100)+0.27*A(109)&
               &+A(119)+0.4*A(121)+A(126)+0.04*A(144)+A(186)+1.6*A(196)
  Vdot(71) = A(162)+A(164)-A(165)-A(166)+A(167)
  Vdot(72) = A(113)-A(114)-A(115)-A(116)+A(117)
  Vdot(73) = -A(29)+0.5*A(98)+0.2*A(116)+0.22*A(147)+0.22*A(149)+0.23*A(151)+0.22*A(152)+0.25*A(178)+0.25*A(179)+0.1&
               &*A(181)+0.25*A(182)+0.5*A(187)-A(189)
  Vdot(74) = -A(19)+0.288*A(28)+A(37)+0.4*A(136)+0.23*A(137)+0.25*A(138)+0.19*A(141)+0.25*A(142)-A(145)-A(146)+A(163)&
               &+A(165)+0.167*A(170)+0.167*A(171)+0.167*A(172)
  Vdot(75) = -A(25)+0.18*A(32)+0.45*A(36)+0.5*A(116)+A(124)+0.45*A(132)+0.95*A(144)+0.8*A(146)+0.25*A(147)+0.25*A(149)&
               &+0.24*A(151)+0.25*A(152)-A(173)-A(174)+0.25*A(178)+0.25*A(179)+0.1*A(181)+0.25*A(182)+A(189)
  Vdot(76) = 0.5*A(140)+A(177)-A(178)-A(179)-A(180)-A(181)-A(182)+A(183)
  Vdot(77) = -A(20)+0.402*A(28)+A(37)+0.2*A(136)+0.32*A(137)+0.35*A(138)+0.26*A(141)+0.35*A(142)-A(143)-A(144)+A(163)&
               &+A(165)+0.039*A(170)+0.039*A(171)+0.039*A(172)
  Vdot(78) = A(143)+0.5*A(145)-A(147)-A(148)-A(149)-A(150)-A(151)-A(152)+0.2*A(153)
  Vdot(79) = A(10)-A(11)-A(12)+A(15)+0.67*A(19)+A(23)+A(27)+0.69*A(28)+A(29)+A(30)+0.1*A(34)+0.25*A(72)+A(73)+2*A(74)&
               &+A(75)+0.3*A(77)-A(78)-A(79)+0.5*A(81)+A(82)+2*A(87)+0.7*A(91)+0.54*A(94)+A(96)+A(104)+0.5*A(105)+A(111)&
               &+A(114)+0.8*A(116)+0.5*A(119)+0.1*A(121)+0.6*A(136)+0.55*A(137)+0.6*A(138)+1.2*A(141)+0.6*A(142)+0.8*A(144)&
               &+0.7*A(146)+0.25*A(147)+0.25*A(149)+0.88*A(151)+0.25*A(152)+A(154)+A(155)+2*A(157)+A(158)+2*A(159)+0.072&
               &*A(170)+0.072*A(171)+0.008*A(172)+0.7*A(181)+A(185)+0.5*A(187)+A(188)+0.8*A(190)
  Vdot(80) = A(135)-A(137)-A(138)-A(139)+0.5*A(140)-A(141)-A(142)
  Vdot(81) = A(6)-A(8)+0.33*A(9)+0.4*A(17)+A(58)-A(59)-A(60)+A(61)+A(63)-A(64)-A(69)-A(78)-A(95)-A(100)-A(138)-A(149)&
               &-A(155)-A(164)-A(169)-A(171)-A(174)-A(176)-A(179)+0.5*A(187)+A(188)-A(193)
  Vdot(82) = A(14)+A(16)+0.4*A(17)+0.3*A(20)+A(24)+A(71)+0.75*A(72)-A(73)-2*A(74)-2*A(75)-A(76)+0.7*A(77)-A(91)+0.31&
               &*A(94)+A(101)-0.1*A(104)+2*A(107)-A(111)-A(116)-A(141)+A(142)-A(151)+A(152)-A(157)+A(158)+A(168)-A(181)&
               &+A(182)
  Vdot(83) = -A(2)-A(3)+0.89*A(8)+A(38)-A(39)-A(46)-A(47)-A(56)-A(58)-A(82)-A(94)+0.25*A(103)-0.9*A(136)-0.8*A(144)-0.8&
               &*A(146)+0.25*A(156)-A(163)
  Vdot(84) = A(7)+0.33*A(9)+A(10)+2*A(13)+A(15)+A(16)+0.33*A(19)+A(21)+A(22)+A(23)+A(26)+A(34)+A(35)+A(36)+A(37)+2*A(41)&
               &+A(42)-A(43)-A(44)+A(45)-A(46)+A(47)-A(49)-A(50)-2*A(51)-2*A(52)+A(55)+A(59)-A(62)-A(63)-A(66)+0.5*A(70)&
               &-A(71)+0.75*A(72)-0.7*A(77)-A(79)-A(80)-A(81)+0.12*A(82)-A(83)-A(84)-A(88)-0.5*A(92)-A(93)+0.33*A(94)-0.5&
               &*A(98)-A(99)-A(105)-A(108)-A(112)-A(113)-A(117)-A(118)-A(120)-A(123)-A(124)-A(125)-A(128)-A(129)-A(130)&
               &-A(134)-A(135)+0.27*A(136)-A(140)-A(143)+0.08*A(144)-A(145)+0.215*A(146)-0.9*A(153)-A(162)+0.7*A(163)-A(167)&
               &-A(168)-A(173)-A(175)-A(177)-A(183)-A(185)-A(186)-A(187)-A(188)-A(189)-A(190)-A(191)-A(192)-A(194)
  Vdot(85) = A(18)+0.33*A(19)+0.2*A(136)+0.5*A(145)+0.5*A(153)-A(154)-A(155)-A(156)-A(157)-A(158)-2*A(159)-A(160)+A(161)
  Vdot(86) = 0.66*A(9)+A(10)+2*A(11)+A(14)+A(15)+0.67*A(19)+A(21)+A(22)+A(25)+A(27)+A(28)+A(29)+2*A(30)+0.56*A(32)+2&
               &*A(33)+0.9*A(34)+A(37)+A(42)+A(43)+A(44)-A(45)+A(46)-A(47)-2*A(48)+A(49)-A(50)-A(55)-A(59)-A(65)+A(67)+0.4&
               &*A(72)+A(73)+2*A(74)-A(76)+A(78)+A(79)+A(80)+0.25*A(81)+0.12*A(82)+A(84)+A(86)+A(87)+A(89)-A(90)+A(91)+0.19&
               &*A(94)+A(96)-A(97)-A(103)+0.9*A(104)+A(109)-A(110)+A(111)-A(115)+0.3*A(116)+A(119)+0.9*A(121)-A(122)-A(127)&
               &+0.25*A(129)+0.7*A(131)+0.9*A(132)-A(133)+0.06*A(136)+A(137)+A(138)-A(139)+A(141)+A(142)+0.06*A(144)+0.275&
               &*A(146)+0.47*A(147)+0.47*A(149)-A(150)+0.73*A(151)+0.47*A(152)+0.2*A(153)-A(156)+A(157)+A(163)+A(165)-A(166)&
               &+0.794*A(170)+0.794*A(171)-0.206*A(172)+A(175)+A(176)+1.5*A(178)+1.5*A(179)-A(180)+A(181)+1.5*A(182)+A(185)&
               &+A(186)+0.5*A(187)+A(189)+A(190)+0.5*A(192)-A(195)+1.2*A(196)
  Vdot(87) = A(5)+0.11*A(8)+2*A(53)-A(55)-A(56)+A(57)-A(64)+0.5*A(70)-A(73)-A(85)-A(89)-A(96)-A(101)-A(109)-A(114)&
               &-A(119)-A(121)-A(126)-A(132)-A(137)-A(147)-A(148)-A(154)-A(165)-A(170)-A(178)
  Vdot(88) = -A(5)+A(6)+A(7)+0.89*A(8)+0.66*A(9)+0.6*A(17)+A(18)+A(27)+A(55)+A(56)-A(57)-A(58)+A(59)-A(60)+A(61)-A(62)+2&
               &*A(64)-A(65)+A(66)+A(67)-A(70)+A(73)+A(85)+A(89)+A(96)+A(101)-A(102)+A(106)+A(109)+A(114)+A(119)+0.9*A(121)&
               &+A(124)+A(126)-0.3*A(131)+0.9*A(132)+0.92*A(137)+A(138)+A(147)+A(149)+A(154)+A(155)-A(160)+A(161)+A(164)&
               &+A(165)+1.206*A(170)+1.206*A(171)+0.206*A(172)+0.4*A(175)+A(176)+A(178)+A(179)
  Vdot(89) = 0.6*A(17)+0.67*A(19)+0.3*A(20)+A(23)+A(24)+A(25)+A(29)+A(31)+0.13*A(32)+A(35)+A(99)+A(100)-A(101)-A(102)&
               &-A(103)-A(104)+0.5*A(105)+A(106)-2*A(107)+A(114)+0.3*A(116)+A(126)-A(142)+0.53*A(147)+0.53*A(149)+0.26&
               &*A(151)-0.47*A(152)+A(154)+A(155)+A(157)+2*A(159)+A(173)+A(174)-A(182)
END SUBROUTINE mozart_mosaic_4bin_vbs0_Fun
SUBROUTINE mozart_mosaic_4bin_vbs0_Jac_SP ( V, F, RCT, JVS )
  REAL(kind=dp) :: V(NVAR)
  REAL(kind=dp) :: F(NFIX)
  REAL(kind=dp) :: RCT(NREACT)
  REAL(kind=dp) :: JVS(LU_NONZERO)
  REAL(kind=dp) :: B(350)
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
  B(21) = RCT(21)
  B(22) = RCT(22)
  B(23) = RCT(23)
  B(24) = RCT(24)
  B(25) = RCT(25)
  B(26) = RCT(26)
  B(27) = RCT(27)
  B(28) = RCT(28)
  B(29) = RCT(29)
  B(30) = RCT(30)
  B(31) = RCT(31)
  B(32) = RCT(32)
  B(33) = RCT(33)
  B(34) = RCT(34)
  B(35) = RCT(35)
  B(36) = RCT(36)
  B(37) = RCT(37)
  B(38) = RCT(38)*F(2)
  B(40) = RCT(39)*V(83)
  B(41) = RCT(39)*V(60)
  B(42) = RCT(40)*F(2)
  B(44) = RCT(41)*F(1)
  B(46) = RCT(42)*V(55)
  B(47) = RCT(42)*V(27)
  B(48) = RCT(43)*V(84)
  B(49) = RCT(43)*V(27)
  B(50) = RCT(44)*V(84)
  B(51) = RCT(44)*V(60)
  B(52) = RCT(45)*V(86)
  B(53) = RCT(45)*V(60)
  B(54) = RCT(46)*V(84)
  B(55) = RCT(46)*V(83)
  B(56) = RCT(47)*V(86)
  B(57) = RCT(47)*V(83)
  B(58) = RCT(48)*2*V(86)*F(1)
  B(60) = RCT(49)*V(84)
  B(61) = RCT(49)*V(19)
  B(62) = RCT(50)*V(86)
  B(63) = RCT(50)*V(84)
  B(64) = RCT(51)*2*V(84)
  B(65) = RCT(52)*2*V(84)
  B(66) = RCT(53)*V(55)
  B(67) = RCT(53)*V(15)
  B(68) = RCT(54)*V(55)
  B(69) = RCT(54)*V(15)
  B(70) = RCT(55)*V(87)
  B(71) = RCT(55)*V(86)
  B(72) = RCT(56)*V(87)
  B(73) = RCT(56)*V(83)
  B(74) = RCT(57)*V(88)
  B(75) = RCT(57)*V(60)
  B(76) = RCT(58)*V(88)
  B(77) = RCT(58)*V(83)
  B(78) = RCT(59)*V(86)
  B(79) = RCT(59)*V(81)
  B(80) = RCT(60)*V(88)
  B(81) = RCT(60)*V(81)
  B(82) = RCT(61)
  B(83) = RCT(62)*V(88)
  B(84) = RCT(62)*V(84)
  B(85) = RCT(63)*V(84)
  B(86) = RCT(63)*V(45)
  B(87) = RCT(64)*V(87)
  B(88) = RCT(64)*V(81)
  B(89) = RCT(65)*V(88)
  B(90) = RCT(65)*V(86)
  B(91) = RCT(66)*V(84)
  B(92) = RCT(66)*V(34)
  B(93) = RCT(67)
  B(94) = RCT(68)*F(2)
  B(96) = RCT(69)
  B(97) = RCT(70)
  B(98) = RCT(71)*V(84)
  B(99) = RCT(71)*V(48)
  B(100) = RCT(72)*V(55)
  B(101) = RCT(72)*V(48)
  B(102) = RCT(73)*V(87)
  B(103) = RCT(73)*V(82)
  B(104) = RCT(74)*2*V(82)
  B(105) = RCT(75)*2*V(82)
  B(106) = RCT(76)*V(86)
  B(107) = RCT(76)*V(82)
  B(108) = RCT(77)*V(84)
  B(109) = RCT(77)*V(33)
  B(110) = RCT(78)*V(81)
  B(111) = RCT(78)*V(79)
  B(112) = RCT(79)*V(84)
  B(113) = RCT(79)*V(79)
  B(114) = RCT(80)*V(84)
  B(115) = RCT(80)*V(65)
  B(116) = RCT(81)*V(84)
  B(117) = RCT(81)*V(36)
  B(118) = RCT(82)*V(83)
  B(119) = RCT(82)*V(36)
  B(120) = RCT(83)*V(84)
  B(121) = RCT(83)*V(16)
  B(122) = RCT(84)*V(84)
  B(123) = RCT(84)*V(43)
  B(124) = RCT(85)*V(87)
  B(125) = RCT(85)*V(39)
  B(126) = RCT(86)*F(2)
  B(128) = RCT(87)
  B(129) = RCT(88)*V(84)
  B(130) = RCT(88)*V(11)
  B(131) = RCT(89)*V(87)
  B(132) = RCT(89)*V(64)
  B(133) = RCT(90)*V(86)
  B(134) = RCT(90)*V(64)
  B(135) = RCT(91)*V(82)
  B(136) = RCT(91)*V(64)
  B(137) = RCT(92)*V(84)
  B(138) = RCT(92)*V(28)
  B(139) = RCT(93)*V(84)
  B(140) = RCT(93)*V(69)
  B(141) = RCT(94)*V(83)
  B(142) = RCT(94)*V(69)
  B(143) = RCT(95)*V(81)
  B(144) = RCT(95)*V(69)
  B(145) = RCT(96)*V(87)
  B(146) = RCT(96)*V(59)
  B(147) = RCT(97)*V(86)
  B(148) = RCT(97)*V(59)
  B(149) = RCT(98)*V(84)
  B(150) = RCT(98)*V(46)
  B(151) = RCT(99)*V(84)
  B(152) = RCT(99)*V(70)
  B(153) = RCT(100)*V(81)
  B(154) = RCT(100)*V(70)
  B(155) = RCT(101)*V(89)
  B(156) = RCT(101)*V(87)
  B(157) = RCT(102)*V(89)
  B(158) = RCT(102)*V(88)
  B(159) = RCT(103)*V(89)
  B(160) = RCT(103)*V(86)
  B(161) = RCT(104)*V(89)
  B(162) = RCT(104)*V(82)
  B(163) = RCT(105)*V(84)
  B(164) = RCT(105)*V(40)
  B(165) = RCT(106)
  B(166) = RCT(107)*2*V(89)
  B(167) = RCT(108)*V(84)
  B(168) = RCT(108)*V(12)
  B(169) = RCT(109)*V(87)
  B(170) = RCT(109)*V(66)
  B(171) = RCT(110)*V(86)
  B(172) = RCT(110)*V(66)
  B(173) = RCT(111)*V(82)
  B(174) = RCT(111)*V(66)
  B(175) = RCT(112)*V(84)
  B(176) = RCT(112)*V(29)
  B(177) = RCT(113)*V(84)
  B(178) = RCT(113)*V(58)
  B(179) = RCT(114)*V(87)
  B(180) = RCT(114)*V(72)
  B(181) = RCT(115)*V(86)
  B(182) = RCT(115)*V(72)
  B(183) = RCT(116)*V(82)
  B(184) = RCT(116)*V(72)
  B(185) = RCT(117)*V(84)
  B(186) = RCT(117)*V(30)
  B(187) = RCT(118)*V(84)
  B(188) = RCT(118)*V(13)
  B(189) = RCT(119)*V(87)
  B(190) = RCT(119)*V(31)
  B(191) = RCT(120)*V(84)
  B(192) = RCT(120)*V(14)
  B(193) = RCT(121)*V(87)
  B(194) = RCT(121)*V(63)
  B(195) = RCT(122)*V(86)
  B(196) = RCT(122)*V(63)
  B(197) = RCT(123)*V(84)
  B(198) = RCT(123)*V(56)
  B(199) = RCT(124)*V(84)
  B(200) = RCT(124)*V(37)
  B(201) = RCT(125)*V(84)
  B(202) = RCT(125)*V(41)
  B(203) = RCT(126)*V(87)
  B(204) = RCT(126)*V(57)
  B(205) = RCT(127)*V(86)
  B(206) = RCT(127)*V(57)
  B(207) = RCT(128)*V(84)
  B(208) = RCT(128)*V(32)
  B(209) = RCT(129)*V(84)
  B(210) = RCT(129)*V(17)
  B(211) = RCT(130)*V(84)
  B(212) = RCT(130)*V(18)
  B(213) = RCT(131)*V(88)
  B(214) = RCT(131)*V(22)
  B(215) = RCT(132)*V(87)
  B(216) = RCT(132)*V(52)
  B(217) = RCT(133)*V(86)
  B(218) = RCT(133)*V(52)
  B(219) = RCT(134)*V(84)
  B(220) = RCT(134)*V(38)
  B(221) = RCT(135)*V(84)
  B(222) = RCT(135)*V(67)
  B(223) = RCT(136)*V(83)
  B(224) = RCT(136)*V(67)
  B(225) = RCT(137)*V(87)
  B(226) = RCT(137)*V(80)
  B(227) = RCT(138)*V(81)
  B(228) = RCT(138)*V(80)
  B(229) = RCT(139)*V(86)
  B(230) = RCT(139)*V(80)
  B(231) = RCT(140)*V(84)
  B(232) = RCT(140)*V(49)
  B(233) = RCT(141)*V(82)
  B(234) = RCT(141)*V(80)
  B(235) = RCT(142)*V(89)
  B(236) = RCT(142)*V(80)
  B(237) = RCT(143)*V(84)
  B(238) = RCT(143)*V(77)
  B(239) = RCT(144)*V(83)
  B(240) = RCT(144)*V(77)
  B(241) = RCT(145)*V(84)
  B(242) = RCT(145)*V(74)
  B(243) = RCT(146)*V(83)
  B(244) = RCT(146)*V(74)
  B(245) = RCT(147)*V(87)
  B(246) = RCT(147)*V(78)
  B(247) = RCT(148)*V(87)
  B(248) = RCT(148)*V(78)
  B(249) = RCT(149)*V(81)
  B(250) = RCT(149)*V(78)
  B(251) = RCT(150)*V(86)
  B(252) = RCT(150)*V(78)
  B(253) = RCT(151)*V(82)
  B(254) = RCT(151)*V(78)
  B(255) = RCT(152)*V(89)
  B(256) = RCT(152)*V(78)
  B(257) = RCT(153)*V(84)
  B(258) = RCT(153)*V(26)
  B(259) = RCT(154)*V(87)
  B(260) = RCT(154)*V(85)
  B(261) = RCT(155)*V(85)
  B(262) = RCT(155)*V(81)
  B(263) = RCT(156)*V(86)
  B(264) = RCT(156)*V(85)
  B(265) = RCT(157)*V(85)
  B(266) = RCT(157)*V(82)
  B(267) = RCT(158)*V(89)
  B(268) = RCT(158)*V(85)
  B(269) = RCT(159)*2*V(85)
  B(270) = RCT(160)*V(88)*F(2)
  B(271) = RCT(160)*V(85)*F(2)
  B(273) = RCT(161)*F(2)
  B(275) = RCT(162)*V(84)
  B(276) = RCT(162)*V(54)
  B(277) = RCT(163)*V(83)
  B(278) = RCT(163)*V(54)
  B(279) = RCT(164)*V(81)
  B(280) = RCT(164)*V(54)
  B(281) = RCT(165)*V(87)
  B(282) = RCT(165)*V(71)
  B(283) = RCT(166)*V(86)
  B(284) = RCT(166)*V(71)
  B(285) = RCT(167)*V(84)
  B(286) = RCT(167)*V(44)
  B(287) = RCT(168)*V(84)
  B(288) = RCT(168)*V(42)
  B(289) = RCT(169)*V(81)
  B(290) = RCT(169)*V(67)
  B(291) = RCT(170)*V(87)
  B(292) = RCT(170)*V(61)
  B(293) = RCT(171)*V(81)
  B(294) = RCT(171)*V(61)
  B(295) = RCT(172)*V(86)
  B(296) = RCT(172)*V(61)
  B(297) = RCT(173)*V(84)
  B(298) = RCT(173)*V(75)
  B(299) = RCT(174)*V(81)
  B(300) = RCT(174)*V(75)
  B(301) = RCT(175)*V(84)
  B(302) = RCT(175)*V(68)
  B(303) = RCT(176)*V(81)
  B(304) = RCT(176)*V(68)
  B(305) = RCT(177)*V(84)
  B(306) = RCT(177)*V(35)
  B(307) = RCT(178)*V(87)
  B(308) = RCT(178)*V(76)
  B(309) = RCT(179)*V(81)
  B(310) = RCT(179)*V(76)
  B(311) = RCT(180)*V(86)
  B(312) = RCT(180)*V(76)
  B(313) = RCT(181)*V(82)
  B(314) = RCT(181)*V(76)
  B(315) = RCT(182)*V(89)
  B(316) = RCT(182)*V(76)
  B(317) = RCT(183)*V(84)
  B(318) = RCT(183)*V(23)
  B(319) = RCT(184)*V(84)
  B(320) = RCT(184)*V(23)
  B(321) = RCT(185)*V(84)
  B(322) = RCT(185)*V(51)
  B(323) = RCT(186)*V(84)
  B(324) = RCT(186)*V(25)
  B(325) = RCT(187)*V(84)
  B(326) = RCT(187)*V(50)
  B(327) = RCT(188)*V(84)
  B(328) = RCT(188)*V(47)
  B(329) = RCT(189)*V(84)
  B(330) = RCT(189)*V(73)
  B(331) = RCT(190)*V(84)
  B(332) = RCT(190)*V(62)
  B(333) = RCT(191)*V(84)
  B(334) = RCT(191)*V(24)
  B(335) = RCT(192)*V(84)
  B(336) = RCT(192)*V(24)
  B(337) = RCT(193)*V(81)
  B(338) = RCT(193)*V(24)
  B(339) = RCT(194)*V(84)
  B(340) = RCT(194)*V(10)
  B(341) = RCT(195)
  B(342) = RCT(196)*2*V(64)
  B(343) = RCT(197)*V(84)
  B(344) = RCT(197)*V(7)
  B(345) = RCT(198)*V(84)
  B(346) = RCT(198)*V(9)
  B(347) = RCT(199)*V(84)
  B(348) = RCT(199)*V(67)
  B(349) = RCT(200)*V(84)
  B(350) = RCT(200)*V(54)
  JVS(1) = 0
  JVS(2) = B(120)
  JVS(3) = B(121)
  JVS(4) = 0
  JVS(5) = B(179)
  JVS(6) = B(180)
  JVS(7) = 0
  JVS(8) = B(181)
  JVS(9) = B(182)
  JVS(10) = 0
  JVS(11) = B(347)
  JVS(12) = B(348)
  JVS(13) = 0
  JVS(14) = B(349)
  JVS(15) = B(350)
  JVS(16) = 0
  JVS(17) = B(343)
  JVS(18) = B(344)
  JVS(19) = -B(343)
  JVS(20) = -B(344)
  JVS(21) = 0
  JVS(22) = B(345)
  JVS(23) = B(346)
  JVS(24) = -B(345)
  JVS(25) = -B(346)
  JVS(26) = -B(339)
  JVS(27) = -B(340)
  JVS(28) = -B(129)
  JVS(29) = -B(130)
  JVS(30) = -B(167)
  JVS(31) = -B(168)
  JVS(32) = -B(187)
  JVS(33) = -B(188)
  JVS(34) = -B(191)
  JVS(35) = -B(192)
  JVS(36) = -B(4)-B(66)-B(68)
  JVS(37) = -B(67)-B(69)
  JVS(38) = -B(120)
  JVS(39) = B(333)+0.5*B(335)+B(337)
  JVS(40) = B(338)
  JVS(41) = -B(121)+B(334)+0.5*B(336)
  JVS(42) = -B(209)
  JVS(43) = -B(210)
  JVS(44) = 0.25*B(209)
  JVS(45) = -B(211)
  JVS(46) = 0.25*B(210)-B(212)
  JVS(47) = -B(13)-B(60)
  JVS(48) = -B(61)+B(65)
  JVS(49) = B(58)+0.5*B(341)
  JVS(50) = -B(126)-B(128)
  JVS(51) = B(124)
  JVS(52) = B(125)
  JVS(53) = -B(6)-B(82)-B(94)
  JVS(54) = B(80)
  JVS(55) = B(81)
  JVS(56) = B(211)
  JVS(57) = -B(213)
  JVS(58) = B(212)
  JVS(59) = -B(214)
  JVS(60) = -B(26)-B(317)-B(319)
  JVS(61) = B(311)
  JVS(62) = -B(318)-B(320)
  JVS(63) = B(312)
  JVS(64) = -B(333)-B(335)-B(337)
  JVS(65) = -B(338)
  JVS(66) = -B(334)-B(336)
  JVS(67) = -B(323)
  JVS(68) = 0.2*B(135)+0.4*B(342)
  JVS(69) = 0.2*B(136)
  JVS(70) = -B(324)
  JVS(71) = -B(257)
  JVS(72) = B(251)
  JVS(73) = -B(258)
  JVS(74) = B(252)
  JVS(75) = -B(46)-B(48)
  JVS(76) = 0.05*B(100)
  JVS(77) = -B(47)+0.05*B(101)
  JVS(78) = B(12)
  JVS(79) = -B(49)
  JVS(80) = -B(21)-B(137)
  JVS(81) = B(133)
  JVS(82) = -B(138)
  JVS(83) = B(134)
  JVS(84) = -B(22)-B(175)
  JVS(85) = B(171)
  JVS(86) = -B(176)
  JVS(87) = B(172)
  JVS(88) = -B(23)-B(185)
  JVS(89) = B(181)
  JVS(90) = -B(186)
  JVS(91) = B(182)
  JVS(92) = B(187)
  JVS(93) = -B(189)
  JVS(94) = B(188)
  JVS(95) = -B(190)
  JVS(96) = -B(35)-B(207)
  JVS(97) = B(205)
  JVS(98) = -B(208)
  JVS(99) = B(206)
  JVS(100) = -B(10)-B(108)
  JVS(101) = B(106)
  JVS(102) = -B(109)
  JVS(103) = B(107)
  JVS(104) = -B(9)-B(91)-B(93)
  JVS(105) = -B(92)
  JVS(106) = B(89)
  JVS(107) = B(90)
  JVS(108) = -B(305)
  JVS(109) = B(301)+B(303)
  JVS(110) = 0.37*B(225)+0.4*B(227)+0.3*B(233)+0.4*B(235)
  JVS(111) = 0.4*B(228)+B(304)
  JVS(112) = 0.3*B(234)
  JVS(113) = B(302)-B(306)
  JVS(114) = 0.37*B(226)
  JVS(115) = 0.4*B(236)
  JVS(116) = -B(116)-B(118)
  JVS(117) = -B(119)
  JVS(118) = -B(117)
  JVS(119) = -B(199)
  JVS(120) = 0.1*B(193)
  JVS(121) = B(143)
  JVS(122) = B(144)
  JVS(123) = -B(200)
  JVS(124) = 0.1*B(194)
  JVS(125) = -B(36)-B(219)
  JVS(126) = B(217)
  JVS(127) = -B(220)
  JVS(128) = B(218)
  JVS(129) = 0.75*B(116)
  JVS(130) = -B(124)
  JVS(131) = 0
  JVS(132) = 0.75*B(117)
  JVS(133) = -B(125)
  JVS(134) = -B(16)-B(163)
  JVS(135) = -B(164)
  JVS(136) = 0.75*B(263)
  JVS(137) = 0.75*B(159)+0.75*B(264)
  JVS(138) = 0.75*B(160)
  JVS(139) = -B(31)-B(201)
  JVS(140) = 0.8*B(34)
  JVS(141) = 0.75*B(193)
  JVS(142) = -B(202)
  JVS(143) = 0.75*B(194)
  JVS(144) = 0.25*B(118)
  JVS(145) = -B(287)
  JVS(146) = 0.2*B(223)
  JVS(147) = 0.25*B(141)
  JVS(148) = 0.1*B(161)
  JVS(149) = 0.25*B(119)+0.25*B(142)+0.2*B(224)
  JVS(150) = -B(288)
  JVS(151) = 0.25*B(263)
  JVS(152) = 0.25*B(159)+0.25*B(264)
  JVS(153) = 0.25*B(160)+0.1*B(162)
  JVS(154) = 0.45*B(36)
  JVS(155) = -B(33)-B(122)
  JVS(156) = 0.45*B(215)
  JVS(157) = 0.13*B(32)
  JVS(158) = 0.2*B(331)
  JVS(159) = -B(123)+0.2*B(332)
  JVS(160) = 0
  JVS(161) = 0.45*B(216)
  JVS(162) = -B(37)-B(285)
  JVS(163) = B(283)
  JVS(164) = -B(286)
  JVS(165) = B(284)
  JVS(166) = 2*B(94)
  JVS(167) = B(337)
  JVS(168) = -B(7)-B(85)
  JVS(169) = B(153)
  JVS(170) = B(299)
  JVS(171) = B(110)
  JVS(172) = B(96)+B(111)+B(154)+B(300)+B(338)
  JVS(173) = B(83)-B(86)
  JVS(174) = B(84)+0.5*B(97)
  JVS(175) = -B(15)-B(149)
  JVS(176) = B(147)
  JVS(177) = -B(150)
  JVS(178) = B(148)
  JVS(179) = -B(17)-B(165)-B(327)
  JVS(180) = -B(328)
  JVS(181) = B(157)
  JVS(182) = B(158)
  JVS(183) = -B(98)-B(100)
  JVS(184) = -B(101)
  JVS(185) = 0.08*B(141)
  JVS(186) = 0.08*B(142)
  JVS(187) = -B(99)
  JVS(188) = -B(28)-B(231)
  JVS(189) = B(229)
  JVS(190) = -B(232)
  JVS(191) = B(230)
  JVS(192) = -B(18)-B(273)-B(325)
  JVS(193) = -B(326)
  JVS(194) = B(270)
  JVS(195) = B(271)
  JVS(196) = -B(321)
  JVS(197) = 0.3*B(135)
  JVS(198) = 0.5*B(183)
  JVS(199) = 0.3*B(313)
  JVS(200) = 0.25*B(253)
  JVS(201) = 0.25*B(233)
  JVS(202) = B(105)+0.3*B(136)+0.5*B(184)+0.25*B(234)+0.25*B(254)+0.3*B(314)
  JVS(203) = -B(322)
  JVS(204) = 0.7*B(209)
  JVS(205) = B(219)
  JVS(206) = -B(215)-B(217)
  JVS(207) = 0.7*B(210)+B(220)
  JVS(208) = -B(218)
  JVS(209) = -B(216)
  JVS(210) = 0.7*B(213)
  JVS(211) = 0.9*B(36)
  JVS(212) = 0.9*B(215)
  JVS(213) = -B(32)
  JVS(214) = 0
  JVS(215) = 0
  JVS(216) = 0.9*B(216)
  JVS(217) = 0.7*B(214)
  JVS(218) = -B(275)-B(277)-B(279)
  JVS(219) = -B(280)
  JVS(220) = -B(278)
  JVS(221) = -B(276)
  JVS(222) = B(4)-B(66)-B(68)
  JVS(223) = -B(46)
  JVS(224) = -B(100)
  JVS(225) = -B(42)-B(44)-B(47)-B(67)-B(69)-B(101)
  JVS(226) = 0
  JVS(227) = 0
  JVS(228) = B(2)
  JVS(229) = 0
  JVS(230) = -B(34)-B(197)
  JVS(231) = B(195)
  JVS(232) = -B(198)
  JVS(233) = B(196)
  JVS(234) = B(207)
  JVS(235) = B(201)
  JVS(236) = 0
  JVS(237) = -B(203)-B(205)
  JVS(238) = 0
  JVS(239) = B(202)+B(208)
  JVS(240) = -B(206)
  JVS(241) = -B(204)
  JVS(242) = 0.82*B(22)
  JVS(243) = 0.5*B(189)
  JVS(244) = 0.1*B(37)
  JVS(245) = 0.25*B(34)
  JVS(246) = -B(24)-B(177)
  JVS(247) = 0.25*B(193)
  JVS(248) = 0.82*B(169)+0.82*B(173)
  JVS(249) = 0.1*B(281)
  JVS(250) = 0.82*B(174)
  JVS(251) = -B(178)
  JVS(252) = 0
  JVS(253) = 0.82*B(170)+0.5*B(190)+0.25*B(194)+0.1*B(282)
  JVS(254) = 0.5*B(149)
  JVS(255) = -B(145)-B(147)
  JVS(256) = B(139)
  JVS(257) = B(140)+0.5*B(150)
  JVS(258) = -B(148)
  JVS(259) = -B(146)
  JVS(260) = B(42)
  JVS(261) = -B(38)-B(40)-B(50)-B(52)-B(74)
  JVS(262) = 0
  JVS(263) = 0
  JVS(264) = B(3)-B(41)
  JVS(265) = -B(51)+B(64)
  JVS(266) = -B(53)
  JVS(267) = B(5)-B(75)
  JVS(268) = -B(291)-B(293)-B(295)
  JVS(269) = B(289)
  JVS(270) = B(290)-B(294)
  JVS(271) = -B(296)
  JVS(272) = -B(292)
  JVS(273) = B(126)
  JVS(274) = 0
  JVS(275) = -B(30)-B(331)
  JVS(276) = 0.25*B(307)+0.25*B(309)+0.1*B(313)+0.25*B(315)
  JVS(277) = 0.53*B(245)+0.53*B(249)+0.26*B(253)+0.53*B(255)
  JVS(278) = 0.53*B(250)+0.25*B(310)
  JVS(279) = 0.26*B(254)+0.1*B(314)
  JVS(280) = 0
  JVS(281) = -B(332)
  JVS(282) = 0.53*B(246)+0.25*B(308)
  JVS(283) = 0.53*B(256)+0.25*B(316)
  JVS(284) = B(191)
  JVS(285) = B(197)
  JVS(286) = -B(193)-B(195)
  JVS(287) = B(192)+B(198)
  JVS(288) = -B(196)
  JVS(289) = -B(194)
  JVS(290) = B(129)
  JVS(291) = 0.5*B(137)
  JVS(292) = B(31)
  JVS(293) = 0
  JVS(294) = 0
  JVS(295) = -B(131)-B(133)-B(135)-2*B(342)
  JVS(296) = -B(136)
  JVS(297) = B(130)+0.5*B(138)
  JVS(298) = -B(134)
  JVS(299) = -B(132)
  JVS(300) = 0.5*B(118)
  JVS(301) = 2*B(33)+B(122)
  JVS(302) = 0
  JVS(303) = 0.45*B(32)
  JVS(304) = B(30)
  JVS(305) = -B(114)
  JVS(306) = 0.3*B(223)
  JVS(307) = B(27)
  JVS(308) = 0.56*B(141)
  JVS(309) = B(14)
  JVS(310) = 0.67*B(19)+0.2*B(243)
  JVS(311) = B(25)+B(297)+B(299)
  JVS(312) = B(307)+B(309)+0.4*B(313)+B(315)
  JVS(313) = 0.7*B(20)+0.05*B(239)
  JVS(314) = 0.22*B(245)+0.22*B(249)+0.11*B(253)+0.22*B(255)
  JVS(315) = B(11)+B(12)+B(110)+B(112)
  JVS(316) = B(111)+0.22*B(250)+B(300)+B(310)
  JVS(317) = 0.11*B(254)+0.4*B(314)
  JVS(318) = 0.5*B(119)+0.56*B(142)+0.3*B(224)+0.05*B(240)+0.2*B(244)
  JVS(319) = B(113)-B(115)+B(123)+B(298)
  JVS(320) = 0
  JVS(321) = 0.22*B(246)+B(308)
  JVS(322) = 0
  JVS(323) = 0.22*B(256)+B(316)
  JVS(324) = B(167)
  JVS(325) = B(175)
  JVS(326) = -B(169)-B(171)-B(173)
  JVS(327) = -B(174)
  JVS(328) = B(168)+B(176)
  JVS(329) = -B(172)
  JVS(330) = -B(170)
  JVS(331) = -B(221)-B(223)-B(289)
  JVS(332) = -B(290)
  JVS(333) = -B(224)
  JVS(334) = -B(222)
  JVS(335) = 0.794*B(291)+0.794*B(293)+0.794*B(295)
  JVS(336) = 0
  JVS(337) = -B(27)-B(301)-B(303)
  JVS(338) = 0.8*B(247)
  JVS(339) = 0.08*B(225)
  JVS(340) = 0.794*B(294)-B(304)
  JVS(341) = 0
  JVS(342) = -B(302)
  JVS(343) = 0.794*B(296)
  JVS(344) = 0.08*B(226)+0.8*B(248)+0.794*B(292)
  JVS(345) = 0.07*B(223)
  JVS(346) = -B(139)-B(141)-B(143)
  JVS(347) = 0.7*B(20)
  JVS(348) = -B(144)
  JVS(349) = -B(142)+0.07*B(224)
  JVS(350) = -B(140)
  JVS(351) = B(323)
  JVS(352) = B(21)+0.5*B(137)
  JVS(353) = B(189)
  JVS(354) = B(35)
  JVS(355) = B(15)
  JVS(356) = 0.4*B(34)
  JVS(357) = B(203)
  JVS(358) = B(145)
  JVS(359) = 0.4*B(193)
  JVS(360) = B(131)+0.8*B(135)+1.6*B(342)
  JVS(361) = 0.27*B(169)
  JVS(362) = 0.5*B(141)
  JVS(363) = -B(14)-B(151)-B(153)
  JVS(364) = 0.04*B(239)
  JVS(365) = -B(154)
  JVS(366) = 0.8*B(136)
  JVS(367) = 0.5*B(142)+0.04*B(240)
  JVS(368) = 0.5*B(138)-B(152)+B(324)
  JVS(369) = 0
  JVS(370) = B(132)+B(146)+0.27*B(170)+B(190)+0.4*B(194)+B(204)
  JVS(371) = B(285)
  JVS(372) = B(275)+B(279)
  JVS(373) = -B(281)-B(283)
  JVS(374) = B(280)
  JVS(375) = 0
  JVS(376) = B(276)+B(286)
  JVS(377) = -B(284)
  JVS(378) = -B(282)
  JVS(379) = B(185)
  JVS(380) = B(177)
  JVS(381) = 0
  JVS(382) = 0
  JVS(383) = 0
  JVS(384) = -B(179)-B(181)-B(183)
  JVS(385) = 0
  JVS(386) = -B(184)
  JVS(387) = 0
  JVS(388) = B(178)+B(186)
  JVS(389) = -B(182)
  JVS(390) = -B(180)
  JVS(391) = 0.5*B(149)
  JVS(392) = 0.5*B(325)
  JVS(393) = 0
  JVS(394) = 0
  JVS(395) = 0.2*B(183)
  JVS(396) = -B(29)-B(329)
  JVS(397) = 0.25*B(307)+0.25*B(309)+0.1*B(313)+0.25*B(315)
  JVS(398) = 0
  JVS(399) = 0.22*B(245)+0.22*B(249)+0.23*B(253)+0.22*B(255)
  JVS(400) = 0.22*B(250)+0.25*B(310)
  JVS(401) = 0.2*B(184)+0.23*B(254)+0.1*B(314)
  JVS(402) = 0
  JVS(403) = 0.5*B(150)+0.5*B(326)-B(330)
  JVS(404) = 0
  JVS(405) = 0
  JVS(406) = 0.22*B(246)+0.25*B(308)
  JVS(407) = 0
  JVS(408) = 0.22*B(256)+0.25*B(316)
  JVS(409) = B(37)
  JVS(410) = 0.288*B(28)
  JVS(411) = B(277)
  JVS(412) = 0.167*B(291)+0.167*B(293)+0.167*B(295)
  JVS(413) = 0.4*B(223)
  JVS(414) = B(281)
  JVS(415) = -B(19)-B(241)-B(243)
  JVS(416) = 0.23*B(225)+0.25*B(227)+0.19*B(233)+0.25*B(235)
  JVS(417) = 0.25*B(228)+0.167*B(294)
  JVS(418) = 0.19*B(234)
  JVS(419) = 0.4*B(224)-B(244)+B(278)
  JVS(420) = -B(242)
  JVS(421) = 0.167*B(296)
  JVS(422) = 0.23*B(226)+B(282)+0.167*B(292)
  JVS(423) = 0.25*B(236)
  JVS(424) = B(199)
  JVS(425) = 0.45*B(36)
  JVS(426) = 0.45*B(215)
  JVS(427) = 0.18*B(32)
  JVS(428) = 0
  JVS(429) = 0
  JVS(430) = 0.5*B(183)
  JVS(431) = B(329)
  JVS(432) = 0.8*B(243)
  JVS(433) = -B(25)-B(297)-B(299)
  JVS(434) = 0.25*B(307)+0.25*B(309)+0.1*B(313)+0.25*B(315)
  JVS(435) = 0.95*B(239)
  JVS(436) = 0.25*B(245)+0.25*B(249)+0.24*B(253)+0.25*B(255)
  JVS(437) = 0
  JVS(438) = 0.25*B(250)-B(300)+0.25*B(310)
  JVS(439) = 0.5*B(184)+0.24*B(254)+0.1*B(314)
  JVS(440) = 0.95*B(240)+0.8*B(244)
  JVS(441) = B(200)-B(298)+B(330)
  JVS(442) = 0
  JVS(443) = 0
  JVS(444) = 0.45*B(216)+0.25*B(246)+0.25*B(308)
  JVS(445) = 0
  JVS(446) = 0.25*B(256)+0.25*B(316)
  JVS(447) = B(317)
  JVS(448) = B(305)
  JVS(449) = 0.5*B(231)
  JVS(450) = 0
  JVS(451) = -B(307)-B(309)-B(311)-B(313)-B(315)
  JVS(452) = 0
  JVS(453) = 0
  JVS(454) = -B(310)
  JVS(455) = -B(314)
  JVS(456) = 0
  JVS(457) = 0.5*B(232)+B(306)+B(318)
  JVS(458) = -B(312)
  JVS(459) = -B(308)
  JVS(460) = -B(316)
  JVS(461) = B(37)
  JVS(462) = 0.402*B(28)
  JVS(463) = B(277)
  JVS(464) = 0.039*B(291)+0.039*B(293)+0.039*B(295)
  JVS(465) = 0.2*B(223)
  JVS(466) = B(281)
  JVS(467) = -B(20)-B(237)-B(239)
  JVS(468) = 0.32*B(225)+0.35*B(227)+0.26*B(233)+0.35*B(235)
  JVS(469) = 0.35*B(228)+0.039*B(294)
  JVS(470) = 0.26*B(234)
  JVS(471) = 0.2*B(224)-B(240)+B(278)
  JVS(472) = -B(238)
  JVS(473) = 0.039*B(296)
  JVS(474) = 0.32*B(226)+B(282)+0.039*B(292)
  JVS(475) = 0.35*B(236)
  JVS(476) = 0.2*B(257)
  JVS(477) = 0.5*B(241)
  JVS(478) = B(237)
  JVS(479) = -B(245)-B(247)-B(249)-B(251)-B(253)-B(255)
  JVS(480) = 0
  JVS(481) = -B(250)
  JVS(482) = -B(254)
  JVS(483) = 0
  JVS(484) = B(238)+0.5*B(242)+0.2*B(258)
  JVS(485) = -B(252)
  JVS(486) = -B(246)-B(248)
  JVS(487) = -B(256)
  JVS(488) = 2*B(128)
  JVS(489) = B(23)
  JVS(490) = 0.5*B(189)
  JVS(491) = B(10)+0.3*B(108)
  JVS(492) = 0.5*B(116)+B(118)
  JVS(493) = 0
  JVS(494) = 0.5*B(163)
  JVS(495) = B(15)
  JVS(496) = B(327)
  JVS(497) = 0.25*B(100)
  JVS(498) = 0.69*B(28)
  JVS(499) = 0.5*B(325)
  JVS(500) = B(321)
  JVS(501) = 0.25*B(101)
  JVS(502) = 0.1*B(34)
  JVS(503) = B(145)
  JVS(504) = 0.072*B(291)+0.072*B(293)+0.008*B(295)
  JVS(505) = B(30)+0.8*B(331)
  JVS(506) = 0.1*B(193)
  JVS(507) = 0.7*B(135)
  JVS(508) = B(173)
  JVS(509) = 0.6*B(223)
  JVS(510) = B(27)
  JVS(511) = 0.54*B(141)
  JVS(512) = B(179)+0.8*B(183)
  JVS(513) = B(29)
  JVS(514) = 0.67*B(19)+0.7*B(243)
  JVS(515) = 0.7*B(313)
  JVS(516) = 0.8*B(239)
  JVS(517) = 0.25*B(245)+0.25*B(249)+0.88*B(253)+0.25*B(255)
  JVS(518) = -B(11)-B(12)-B(110)-B(112)
  JVS(519) = 0.55*B(225)+0.6*B(227)+1.2*B(233)+0.6*B(235)
  JVS(520) = -B(111)+0.6*B(228)+0.25*B(250)+B(261)+0.072*B(294)
  JVS(521) = B(102)+2*B(104)+B(105)+0.7*B(136)+B(161)+B(174)+0.8*B(184)+1.2*B(234)+0.88*B(254)+2*B(265)+0.7*B(314)
  JVS(522) = B(119)+0.54*B(142)+0.6*B(224)+0.8*B(240)+0.7*B(244)
  JVS(523) = 0.3*B(109)-B(113)+0.5*B(117)+0.5*B(164)+B(322)+0.5*B(326)+B(328)+0.8*B(332)
  JVS(524) = B(259)+B(262)+2*B(266)+B(267)+2*B(269)
  JVS(525) = 0.008*B(296)
  JVS(526) = B(103)+B(146)+B(180)+0.5*B(190)+0.1*B(194)+0.55*B(226)+0.25*B(246)+B(260)+0.072*B(292)
  JVS(527) = 0
  JVS(528) = B(162)+0.6*B(236)+0.25*B(256)+B(268)
  JVS(529) = 0.5*B(231)
  JVS(530) = B(221)
  JVS(531) = -B(225)-B(227)-B(229)-B(233)-B(235)
  JVS(532) = -B(228)
  JVS(533) = -B(234)
  JVS(534) = 0
  JVS(535) = B(222)+0.5*B(232)
  JVS(536) = -B(230)
  JVS(537) = -B(226)
  JVS(538) = -B(236)
  JVS(539) = B(6)+B(82)
  JVS(540) = -B(337)
  JVS(541) = 0.33*B(9)
  JVS(542) = B(85)
  JVS(543) = 0.4*B(17)+B(327)
  JVS(544) = 0.5*B(325)
  JVS(545) = -B(279)
  JVS(546) = -B(293)
  JVS(547) = -B(289)
  JVS(548) = -B(303)
  JVS(549) = -B(143)
  JVS(550) = -B(153)
  JVS(551) = -B(299)
  JVS(552) = -B(309)
  JVS(553) = 0
  JVS(554) = -B(249)
  JVS(555) = -B(110)
  JVS(556) = -B(227)
  JVS(557) = -B(8)-B(78)-B(80)-B(87)-B(96)-B(111)-B(144)-B(154)-B(228)-B(250)-B(261)-B(280)-B(290)-B(294)-B(300)-B(304)&
               &-B(310)-B(338)
  JVS(558) = 0
  JVS(559) = B(76)
  JVS(560) = B(86)+0.5*B(326)+B(328)
  JVS(561) = -B(262)
  JVS(562) = -B(79)
  JVS(563) = -B(88)
  JVS(564) = B(77)-B(81)
  JVS(565) = 0
  JVS(566) = 0.7*B(108)
  JVS(567) = B(16)
  JVS(568) = B(287)
  JVS(569) = 0.4*B(17)
  JVS(570) = B(98)+0.75*B(100)
  JVS(571) = 0.75*B(101)
  JVS(572) = B(24)
  JVS(573) = 0
  JVS(574) = -B(135)
  JVS(575) = -B(173)
  JVS(576) = 0
  JVS(577) = 0.31*B(141)
  JVS(578) = B(14)
  JVS(579) = 0
  JVS(580) = -B(183)
  JVS(581) = -B(313)+B(315)
  JVS(582) = 0.3*B(20)
  JVS(583) = -B(253)+B(255)
  JVS(584) = 0
  JVS(585) = -B(233)+B(235)
  JVS(586) = 0
  JVS(587) = -B(102)-2*B(104)-2*B(105)-B(106)-B(136)-0.1*B(161)-B(174)-B(184)-B(234)-B(254)-B(265)-B(314)
  JVS(588) = 0.31*B(142)
  JVS(589) = B(99)+0.7*B(109)+B(288)
  JVS(590) = -B(266)+B(267)
  JVS(591) = -B(107)
  JVS(592) = -B(103)+B(155)
  JVS(593) = 0
  JVS(594) = B(156)-0.1*B(162)+2*B(166)+B(236)+B(256)+B(268)+B(316)
  JVS(595) = -B(118)
  JVS(596) = -B(277)
  JVS(597) = B(38)-B(40)
  JVS(598) = -0.9*B(223)
  JVS(599) = -B(141)
  JVS(600) = -0.8*B(243)
  JVS(601) = -0.8*B(239)
  JVS(602) = 0
  JVS(603) = 0
  JVS(604) = 0.89*B(8)
  JVS(605) = 0
  JVS(606) = -B(2)-B(3)-B(41)-B(54)-B(56)-B(72)-B(76)-B(119)-B(142)-0.9*B(224)-0.8*B(240)-0.8*B(244)-B(278)
  JVS(607) = -B(55)
  JVS(608) = 0.25*B(263)
  JVS(609) = -B(57)+0.25*B(159)+0.25*B(264)
  JVS(610) = -B(73)
  JVS(611) = -B(77)
  JVS(612) = 0.25*B(160)
  JVS(613) = -B(339)
  JVS(614) = -B(129)
  JVS(615) = -B(167)
  JVS(616) = -B(187)
  JVS(617) = -B(191)
  JVS(618) = -B(120)
  JVS(619) = -B(209)
  JVS(620) = -B(211)
  JVS(621) = 2*B(13)-B(60)
  JVS(622) = B(26)-B(317)
  JVS(623) = -B(333)-B(335)
  JVS(624) = -B(323)
  JVS(625) = -0.9*B(257)
  JVS(626) = B(46)-B(48)
  JVS(627) = B(21)-0.5*B(137)
  JVS(628) = B(22)-B(175)
  JVS(629) = B(23)-B(185)
  JVS(630) = B(35)-B(207)
  JVS(631) = B(10)-0.7*B(108)
  JVS(632) = 0.33*B(9)-B(91)
  JVS(633) = -B(305)
  JVS(634) = -B(116)+0.12*B(118)
  JVS(635) = -B(199)
  JVS(636) = B(36)-B(219)
  JVS(637) = B(16)-B(163)
  JVS(638) = -B(201)
  JVS(639) = -B(287)
  JVS(640) = -B(122)
  JVS(641) = B(37)-B(285)
  JVS(642) = B(7)-B(85)
  JVS(643) = B(15)-0.5*B(149)
  JVS(644) = -B(327)
  JVS(645) = -B(98)+0.75*B(100)
  JVS(646) = -B(231)
  JVS(647) = -B(325)
  JVS(648) = -B(321)
  JVS(649) = 0
  JVS(650) = 0
  JVS(651) = -B(275)+0.7*B(277)
  JVS(652) = 2*B(44)+B(47)+0.75*B(101)
  JVS(653) = B(34)-B(197)
  JVS(654) = 0
  JVS(655) = -B(177)
  JVS(656) = 0
  JVS(657) = -B(50)+B(52)
  JVS(658) = -B(331)
  JVS(659) = 0
  JVS(660) = 0
  JVS(661) = -B(114)
  JVS(662) = 0
  JVS(663) = -B(221)+0.27*B(223)
  JVS(664) = -B(301)
  JVS(665) = -B(139)+0.33*B(141)
  JVS(666) = -B(151)
  JVS(667) = 0
  JVS(668) = 0
  JVS(669) = -B(329)
  JVS(670) = 0.33*B(19)-B(241)+0.215*B(243)
  JVS(671) = -B(297)
  JVS(672) = 0
  JVS(673) = -B(237)+0.08*B(239)
  JVS(674) = 0
  JVS(675) = -B(112)
  JVS(676) = 0
  JVS(677) = B(78)
  JVS(678) = 0
  JVS(679) = -B(54)+B(56)+0.12*B(119)+0.33*B(142)+0.27*B(224)+0.08*B(240)+0.215*B(244)+0.7*B(278)
  JVS(680) = -B(49)-B(51)-B(55)-B(61)-B(62)-2*B(64)-2*B(65)-B(83)-B(86)-B(92)-B(99)-0.7*B(109)-B(113)-B(115)-B(117)&
               &-B(121)-B(123)-B(130)-0.5*B(138)-B(140)-0.5*B(150)-B(152)-B(164)-B(168)-B(176)-B(178)-B(186)-B(188)-B(192)&
               &-B(198)-B(200)-B(202)-B(208)-B(210)-B(212)-B(220)-B(222)-B(232)-B(238)-B(242)-0.9*B(258)-B(276)-B(286)&
               &-B(288)-B(298)-B(302)-B(306)-B(318)-B(322)-B(324)-B(326)-B(328)-B(330)-B(332)-B(334)-B(336)-B(340)
  JVS(681) = 0
  JVS(682) = B(53)+B(57)-B(63)+B(70)+B(79)
  JVS(683) = B(71)
  JVS(684) = -B(84)+0.5*B(97)
  JVS(685) = 0
  JVS(686) = 0.5*B(257)
  JVS(687) = B(18)+B(273)
  JVS(688) = 0.2*B(223)
  JVS(689) = 0.33*B(19)+0.5*B(241)
  JVS(690) = 0
  JVS(691) = 0
  JVS(692) = -B(261)
  JVS(693) = -B(265)
  JVS(694) = 0.2*B(224)
  JVS(695) = 0.5*B(242)+0.5*B(258)
  JVS(696) = -B(259)-B(262)-B(263)-B(266)-B(267)-2*B(269)-B(270)
  JVS(697) = -B(264)
  JVS(698) = -B(260)
  JVS(699) = -B(271)
  JVS(700) = -B(268)
  JVS(701) = 0.25*B(209)
  JVS(702) = B(60)
  JVS(703) = B(126)+B(128)
  JVS(704) = 0.7*B(213)
  JVS(705) = 0.5*B(335)
  JVS(706) = B(323)
  JVS(707) = 0.2*B(257)
  JVS(708) = B(46)+B(48)
  JVS(709) = B(21)
  JVS(710) = B(22)
  JVS(711) = B(189)
  JVS(712) = B(10)
  JVS(713) = 0.66*B(9)+B(93)
  JVS(714) = 0.25*B(116)+0.12*B(118)
  JVS(715) = 0
  JVS(716) = 2*B(33)+B(122)
  JVS(717) = B(37)
  JVS(718) = B(15)
  JVS(719) = 0.4*B(100)
  JVS(720) = B(28)
  JVS(721) = 0.5*B(325)
  JVS(722) = B(321)
  JVS(723) = 0.9*B(215)-B(217)
  JVS(724) = 0.56*B(32)
  JVS(725) = B(277)
  JVS(726) = B(47)+0.4*B(101)
  JVS(727) = 0.9*B(34)
  JVS(728) = -B(205)
  JVS(729) = B(145)-B(147)
  JVS(730) = B(50)-B(52)
  JVS(731) = 0.794*B(291)+0.794*B(293)-0.206*B(295)
  JVS(732) = 2*B(30)+B(331)
  JVS(733) = 0.9*B(193)-B(195)
  JVS(734) = B(131)-B(133)+B(135)+1.2*B(342)
  JVS(735) = B(114)
  JVS(736) = B(169)-B(171)+B(173)
  JVS(737) = 0.06*B(223)
  JVS(738) = B(27)+B(301)+B(303)
  JVS(739) = 0.19*B(141)
  JVS(740) = B(14)
  JVS(741) = B(281)-B(283)
  JVS(742) = -B(181)+0.3*B(183)
  JVS(743) = B(29)+B(329)
  JVS(744) = 0.67*B(19)+0.275*B(243)
  JVS(745) = B(25)
  JVS(746) = 1.5*B(307)+1.5*B(309)-B(311)+B(313)+1.5*B(315)
  JVS(747) = 0.06*B(239)
  JVS(748) = 0.47*B(245)+0.47*B(249)-B(251)+0.73*B(253)+0.47*B(255)
  JVS(749) = 2*B(11)+B(110)+B(112)
  JVS(750) = B(225)+B(227)-B(229)+B(233)+B(235)
  JVS(751) = -B(78)+B(111)+B(228)+0.47*B(250)+0.794*B(294)+B(304)+1.5*B(310)
  JVS(752) = B(102)+2*B(104)-B(106)+B(136)+0.9*B(161)+B(174)+0.3*B(184)+B(234)+0.73*B(254)+B(265)+B(314)
  JVS(753) = B(54)-B(56)+0.12*B(119)+0.19*B(142)+0.06*B(224)+0.06*B(240)+0.275*B(244)+B(278)
  JVS(754) = B(49)+B(51)+B(55)+B(61)-B(62)+B(113)+B(115)+0.25*B(117)+B(123)+0.25*B(210)+0.2*B(258)+B(302)+B(322)+B(324)&
               &+0.5*B(326)+B(330)+B(332)+0.5*B(336)
  JVS(755) = -B(263)+B(266)
  JVS(756) = -B(53)-B(57)-2*B(58)-B(63)-B(70)-B(79)-B(89)-B(107)-B(134)-B(148)-B(159)-B(172)-B(182)-B(196)-B(206)-B(218)&
               &-B(230)-B(252)-B(264)-B(284)-0.206*B(296)-B(312)-B(341)
  JVS(757) = -B(71)+B(103)+B(132)+B(146)+B(170)+B(190)+0.9*B(194)+0.9*B(216)+B(226)+0.47*B(246)+B(282)+0.794*B(292)+1.5&
               &*B(308)
  JVS(758) = -B(90)+0.7*B(214)
  JVS(759) = -B(160)+0.9*B(162)+B(236)+0.47*B(256)+1.5*B(316)
  JVS(760) = 2*B(66)
  JVS(761) = -B(189)
  JVS(762) = -B(124)
  JVS(763) = -B(215)
  JVS(764) = 2*B(67)
  JVS(765) = -B(203)
  JVS(766) = -B(145)
  JVS(767) = B(74)
  JVS(768) = -B(291)
  JVS(769) = -B(193)
  JVS(770) = -B(131)
  JVS(771) = -B(169)
  JVS(772) = 0
  JVS(773) = 0
  JVS(774) = -B(281)
  JVS(775) = -B(179)
  JVS(776) = -B(307)
  JVS(777) = 0
  JVS(778) = -B(245)-B(247)
  JVS(779) = 0
  JVS(780) = -B(225)
  JVS(781) = 0.11*B(8)-B(87)
  JVS(782) = -B(102)
  JVS(783) = -B(72)
  JVS(784) = 0
  JVS(785) = -B(259)
  JVS(786) = -B(70)
  JVS(787) = -B(71)-B(73)-B(88)-B(103)-B(125)-B(132)-B(146)-B(155)-B(170)-B(180)-B(190)-B(194)-B(204)-B(216)-B(226)&
               &-B(246)-B(248)-B(260)-B(282)-B(292)-B(308)
  JVS(788) = B(5)+B(75)+0.5*B(97)
  JVS(789) = -B(156)
  JVS(790) = B(6)+B(82)
  JVS(791) = -0.3*B(213)
  JVS(792) = B(189)
  JVS(793) = 0.66*B(9)+B(91)+B(93)
  JVS(794) = B(199)
  JVS(795) = B(124)
  JVS(796) = B(7)
  JVS(797) = 0.6*B(17)+B(165)
  JVS(798) = B(18)+B(273)
  JVS(799) = 0.9*B(215)
  JVS(800) = B(279)
  JVS(801) = B(203)
  JVS(802) = B(145)
  JVS(803) = -B(74)
  JVS(804) = 1.206*B(291)+1.206*B(293)+0.206*B(295)
  JVS(805) = 0.9*B(193)
  JVS(806) = B(131)
  JVS(807) = B(169)
  JVS(808) = 0
  JVS(809) = B(27)+0.4*B(301)+B(303)
  JVS(810) = 0
  JVS(811) = 0
  JVS(812) = B(281)
  JVS(813) = B(179)
  JVS(814) = 0
  JVS(815) = B(307)+B(309)
  JVS(816) = 0
  JVS(817) = B(245)+B(249)
  JVS(818) = 0
  JVS(819) = 0.92*B(225)+B(227)
  JVS(820) = 0.89*B(8)+B(78)-B(80)+2*B(87)+B(228)+B(250)+B(261)+B(280)+1.206*B(294)+B(304)+B(310)
  JVS(821) = B(102)
  JVS(822) = B(72)-B(76)
  JVS(823) = -B(83)+B(92)+B(200)+0.4*B(302)
  JVS(824) = B(259)+B(262)-B(270)
  JVS(825) = B(70)+B(79)-B(89)+0.206*B(296)
  JVS(826) = B(71)+B(73)+2*B(88)+B(103)+B(125)+B(132)+B(146)+B(155)+B(170)+B(180)+B(190)+0.9*B(194)+B(204)+0.9*B(216)&
               &+0.92*B(226)+B(246)+B(260)+B(282)+1.206*B(292)+B(308)
  JVS(827) = -B(5)-B(75)-B(77)-B(81)-B(84)-B(90)-B(97)-B(157)-0.3*B(214)-B(271)
  JVS(828) = B(156)-B(158)
  JVS(829) = B(23)
  JVS(830) = B(35)
  JVS(831) = 0.5*B(163)
  JVS(832) = B(31)
  JVS(833) = 0.6*B(17)+B(165)
  JVS(834) = 0.13*B(32)
  JVS(835) = 0
  JVS(836) = B(203)
  JVS(837) = B(24)
  JVS(838) = 0
  JVS(839) = 0
  JVS(840) = B(151)+B(153)
  JVS(841) = 0
  JVS(842) = B(179)+0.3*B(183)
  JVS(843) = B(29)
  JVS(844) = 0.67*B(19)
  JVS(845) = B(25)+B(297)+B(299)
  JVS(846) = -B(315)
  JVS(847) = 0.3*B(20)
  JVS(848) = 0.53*B(245)+0.53*B(249)+0.26*B(253)-0.47*B(255)
  JVS(849) = -B(235)
  JVS(850) = B(154)+0.53*B(250)+B(261)+B(300)
  JVS(851) = -B(161)+0.3*B(184)+0.26*B(254)+B(265)
  JVS(852) = 0
  JVS(853) = B(152)+0.5*B(164)+B(298)
  JVS(854) = B(259)+B(262)+B(266)+2*B(269)
  JVS(855) = -B(159)
  JVS(856) = -B(155)+B(180)+B(204)+0.53*B(246)+B(260)
  JVS(857) = -B(157)
  JVS(858) = -B(156)-B(158)-B(160)-B(162)-2*B(166)-B(236)-0.47*B(256)-B(316)
END SUBROUTINE mozart_mosaic_4bin_vbs0_Jac_SP
SUBROUTINE mozart_mosaic_4bin_vbs0_KppDecomp( JVS, IER )
      INTEGER :: IER
      REAL(kind=dp) :: JVS(858), W(89), a
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
END SUBROUTINE mozart_mosaic_4bin_vbs0_KppDecomp
SUBROUTINE mozart_mosaic_4bin_vbs0_KppDecompCmplx( JVS, IER )
      INTEGER :: IER
      DOUBLE COMPLEX :: JVS(858), W(89), a
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
END SUBROUTINE mozart_mosaic_4bin_vbs0_KppDecompCmplx
SUBROUTINE mozart_mosaic_4bin_vbs0_KppSolveIndirect( JVS, X )
      INTEGER i, j
      REAL(kind=dp) JVS(858), X(89), sum
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
END SUBROUTINE mozart_mosaic_4bin_vbs0_KppSolveIndirect
SUBROUTINE mozart_mosaic_4bin_vbs0_KppSolveCmplx( JVS, X )
      INTEGER i, j
      DOUBLE COMPLEX JVS(858), X(89), sum
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
END SUBROUTINE mozart_mosaic_4bin_vbs0_KppSolveCmplx
SUBROUTINE mozart_mosaic_4bin_vbs0_KppSolve ( JVS, X )
  REAL(kind=dp) :: JVS(LU_NONZERO)
  REAL(kind=dp) :: X(NVAR)
  X(18) = X(18)-JVS(44)*X(17)
  X(22) = X(22)-JVS(56)*X(18)
  X(31) = X(31)-JVS(92)*X(13)
  X(39) = X(39)-JVS(129)*X(36)
  X(42) = X(42)-JVS(144)*X(36)
  X(43) = X(43)-JVS(154)*X(38)
  X(45) = X(45)-JVS(166)*X(21)-JVS(167)*X(24)
  X(52) = X(52)-JVS(204)*X(17)-JVS(205)*X(38)
  X(53) = X(53)-JVS(210)*X(22)-JVS(211)*X(38)-JVS(212)*X(52)
  X(55) = X(55)-JVS(222)*X(15)-JVS(223)*X(27)-JVS(224)*X(48)
  X(57) = X(57)-JVS(234)*X(32)-JVS(235)*X(41)-JVS(236)*X(56)
  X(58) = X(58)-JVS(242)*X(29)-JVS(243)*X(31)-JVS(244)*X(44)-JVS(245)*X(56)
  X(59) = X(59)-JVS(254)*X(46)
  X(60) = X(60)-JVS(260)*X(55)
  X(62) = X(62)-JVS(273)*X(20)-JVS(274)*X(39)
  X(63) = X(63)-JVS(284)*X(14)-JVS(285)*X(56)
  X(64) = X(64)-JVS(290)*X(11)-JVS(291)*X(28)-JVS(292)*X(41)-JVS(293)*X(56)-JVS(294)*X(63)
  X(65) = X(65)-JVS(300)*X(36)-JVS(301)*X(43)-JVS(302)*X(52)-JVS(303)*X(53)-JVS(304)*X(62)
  X(66) = X(66)-JVS(324)*X(12)-JVS(325)*X(29)
  X(68) = X(68)-JVS(335)*X(61)-JVS(336)*X(67)
  X(69) = X(69)-JVS(345)*X(67)
  X(70) = X(70)-JVS(351)*X(25)-JVS(352)*X(28)-JVS(353)*X(31)-JVS(354)*X(32)-JVS(355)*X(46)-JVS(356)*X(56)-JVS(357)*X(57)&
            &-JVS(358)*X(59)-JVS(359)*X(63)-JVS(360)*X(64)-JVS(361)*X(66)-JVS(362)*X(69)
  X(71) = X(71)-JVS(371)*X(44)-JVS(372)*X(54)
  X(72) = X(72)-JVS(379)*X(30)-JVS(380)*X(58)-JVS(381)*X(63)-JVS(382)*X(66)-JVS(383)*X(71)
  X(73) = X(73)-JVS(391)*X(46)-JVS(392)*X(50)-JVS(393)*X(59)-JVS(394)*X(69)-JVS(395)*X(72)
  X(74) = X(74)-JVS(409)*X(44)-JVS(410)*X(49)-JVS(411)*X(54)-JVS(412)*X(61)-JVS(413)*X(67)-JVS(414)*X(71)
  X(75) = X(75)-JVS(424)*X(37)-JVS(425)*X(38)-JVS(426)*X(52)-JVS(427)*X(53)-JVS(428)*X(63)-JVS(429)*X(69)-JVS(430)*X(72)&
            &-JVS(431)*X(73)-JVS(432)*X(74)
  X(76) = X(76)-JVS(447)*X(23)-JVS(448)*X(35)-JVS(449)*X(49)-JVS(450)*X(68)
  X(77) = X(77)-JVS(461)*X(44)-JVS(462)*X(49)-JVS(463)*X(54)-JVS(464)*X(61)-JVS(465)*X(67)-JVS(466)*X(71)
  X(78) = X(78)-JVS(476)*X(26)-JVS(477)*X(74)-JVS(478)*X(77)
  X(79) = X(79)-JVS(488)*X(20)-JVS(489)*X(30)-JVS(490)*X(31)-JVS(491)*X(33)-JVS(492)*X(36)-JVS(493)*X(39)-JVS(494)*X(40)&
            &-JVS(495)*X(46)-JVS(496)*X(47)-JVS(497)*X(48)-JVS(498)*X(49)-JVS(499)*X(50)-JVS(500)*X(51)-JVS(501)*X(55)&
            &-JVS(502)*X(56)-JVS(503)*X(59)-JVS(504)*X(61)-JVS(505)*X(62)-JVS(506)*X(63)-JVS(507)*X(64)-JVS(508)*X(66)&
            &-JVS(509)*X(67)-JVS(510)*X(68)-JVS(511)*X(69)-JVS(512)*X(72)-JVS(513)*X(73)-JVS(514)*X(74)-JVS(515)*X(76)&
            &-JVS(516)*X(77)-JVS(517)*X(78)
  X(80) = X(80)-JVS(529)*X(49)-JVS(530)*X(67)
  X(81) = X(81)-JVS(539)*X(21)-JVS(540)*X(24)-JVS(541)*X(34)-JVS(542)*X(45)-JVS(543)*X(47)-JVS(544)*X(50)-JVS(545)*X(54)&
            &-JVS(546)*X(61)-JVS(547)*X(67)-JVS(548)*X(68)-JVS(549)*X(69)-JVS(550)*X(70)-JVS(551)*X(75)-JVS(552)*X(76)&
            &-JVS(553)*X(77)-JVS(554)*X(78)-JVS(555)*X(79)-JVS(556)*X(80)
  X(82) = X(82)-JVS(566)*X(33)-JVS(567)*X(40)-JVS(568)*X(42)-JVS(569)*X(47)-JVS(570)*X(48)-JVS(571)*X(55)-JVS(572)*X(58)&
            &-JVS(573)*X(63)-JVS(574)*X(64)-JVS(575)*X(66)-JVS(576)*X(67)-JVS(577)*X(69)-JVS(578)*X(70)-JVS(579)*X(71)&
            &-JVS(580)*X(72)-JVS(581)*X(76)-JVS(582)*X(77)-JVS(583)*X(78)-JVS(584)*X(79)-JVS(585)*X(80)-JVS(586)*X(81)
  X(83) = X(83)-JVS(595)*X(36)-JVS(596)*X(54)-JVS(597)*X(60)-JVS(598)*X(67)-JVS(599)*X(69)-JVS(600)*X(74)-JVS(601)*X(77)&
            &-JVS(602)*X(79)-JVS(603)*X(80)-JVS(604)*X(81)-JVS(605)*X(82)
  X(84) = X(84)-JVS(613)*X(10)-JVS(614)*X(11)-JVS(615)*X(12)-JVS(616)*X(13)-JVS(617)*X(14)-JVS(618)*X(16)-JVS(619)*X(17)&
            &-JVS(620)*X(18)-JVS(621)*X(19)-JVS(622)*X(23)-JVS(623)*X(24)-JVS(624)*X(25)-JVS(625)*X(26)-JVS(626)*X(27)&
            &-JVS(627)*X(28)-JVS(628)*X(29)-JVS(629)*X(30)-JVS(630)*X(32)-JVS(631)*X(33)-JVS(632)*X(34)-JVS(633)*X(35)&
            &-JVS(634)*X(36)-JVS(635)*X(37)-JVS(636)*X(38)-JVS(637)*X(40)-JVS(638)*X(41)-JVS(639)*X(42)-JVS(640)*X(43)&
            &-JVS(641)*X(44)-JVS(642)*X(45)-JVS(643)*X(46)-JVS(644)*X(47)-JVS(645)*X(48)-JVS(646)*X(49)-JVS(647)*X(50)&
            &-JVS(648)*X(51)-JVS(649)*X(52)-JVS(650)*X(53)-JVS(651)*X(54)-JVS(652)*X(55)-JVS(653)*X(56)-JVS(654)*X(57)&
            &-JVS(655)*X(58)-JVS(656)*X(59)-JVS(657)*X(60)-JVS(658)*X(62)-JVS(659)*X(63)-JVS(660)*X(64)-JVS(661)*X(65)&
            &-JVS(662)*X(66)-JVS(663)*X(67)-JVS(664)*X(68)-JVS(665)*X(69)-JVS(666)*X(70)-JVS(667)*X(71)-JVS(668)*X(72)&
            &-JVS(669)*X(73)-JVS(670)*X(74)-JVS(671)*X(75)-JVS(672)*X(76)-JVS(673)*X(77)-JVS(674)*X(78)-JVS(675)*X(79)&
            &-JVS(676)*X(80)-JVS(677)*X(81)-JVS(678)*X(82)-JVS(679)*X(83)
  X(85) = X(85)-JVS(686)*X(26)-JVS(687)*X(50)-JVS(688)*X(67)-JVS(689)*X(74)-JVS(690)*X(78)-JVS(691)*X(80)-JVS(692)*X(81)&
            &-JVS(693)*X(82)-JVS(694)*X(83)-JVS(695)*X(84)
  X(86) = X(86)-JVS(701)*X(17)-JVS(702)*X(19)-JVS(703)*X(20)-JVS(704)*X(22)-JVS(705)*X(24)-JVS(706)*X(25)-JVS(707)*X(26)&
            &-JVS(708)*X(27)-JVS(709)*X(28)-JVS(710)*X(29)-JVS(711)*X(31)-JVS(712)*X(33)-JVS(713)*X(34)-JVS(714)*X(36)&
            &-JVS(715)*X(39)-JVS(716)*X(43)-JVS(717)*X(44)-JVS(718)*X(46)-JVS(719)*X(48)-JVS(720)*X(49)-JVS(721)*X(50)&
            &-JVS(722)*X(51)-JVS(723)*X(52)-JVS(724)*X(53)-JVS(725)*X(54)-JVS(726)*X(55)-JVS(727)*X(56)-JVS(728)*X(57)&
            &-JVS(729)*X(59)-JVS(730)*X(60)-JVS(731)*X(61)-JVS(732)*X(62)-JVS(733)*X(63)-JVS(734)*X(64)-JVS(735)*X(65)&
            &-JVS(736)*X(66)-JVS(737)*X(67)-JVS(738)*X(68)-JVS(739)*X(69)-JVS(740)*X(70)-JVS(741)*X(71)-JVS(742)*X(72)&
            &-JVS(743)*X(73)-JVS(744)*X(74)-JVS(745)*X(75)-JVS(746)*X(76)-JVS(747)*X(77)-JVS(748)*X(78)-JVS(749)*X(79)&
            &-JVS(750)*X(80)-JVS(751)*X(81)-JVS(752)*X(82)-JVS(753)*X(83)-JVS(754)*X(84)-JVS(755)*X(85)
  X(87) = X(87)-JVS(760)*X(15)-JVS(761)*X(31)-JVS(762)*X(39)-JVS(763)*X(52)-JVS(764)*X(55)-JVS(765)*X(57)-JVS(766)*X(59)&
            &-JVS(767)*X(60)-JVS(768)*X(61)-JVS(769)*X(63)-JVS(770)*X(64)-JVS(771)*X(66)-JVS(772)*X(67)-JVS(773)*X(69)&
            &-JVS(774)*X(71)-JVS(775)*X(72)-JVS(776)*X(76)-JVS(777)*X(77)-JVS(778)*X(78)-JVS(779)*X(79)-JVS(780)*X(80)&
            &-JVS(781)*X(81)-JVS(782)*X(82)-JVS(783)*X(83)-JVS(784)*X(84)-JVS(785)*X(85)-JVS(786)*X(86)
  X(88) = X(88)-JVS(790)*X(21)-JVS(791)*X(22)-JVS(792)*X(31)-JVS(793)*X(34)-JVS(794)*X(37)-JVS(795)*X(39)-JVS(796)*X(45)&
            &-JVS(797)*X(47)-JVS(798)*X(50)-JVS(799)*X(52)-JVS(800)*X(54)-JVS(801)*X(57)-JVS(802)*X(59)-JVS(803)*X(60)&
            &-JVS(804)*X(61)-JVS(805)*X(63)-JVS(806)*X(64)-JVS(807)*X(66)-JVS(808)*X(67)-JVS(809)*X(68)-JVS(810)*X(69)&
            &-JVS(811)*X(70)-JVS(812)*X(71)-JVS(813)*X(72)-JVS(814)*X(75)-JVS(815)*X(76)-JVS(816)*X(77)-JVS(817)*X(78)&
            &-JVS(818)*X(79)-JVS(819)*X(80)-JVS(820)*X(81)-JVS(821)*X(82)-JVS(822)*X(83)-JVS(823)*X(84)-JVS(824)*X(85)&
            &-JVS(825)*X(86)-JVS(826)*X(87)
  X(89) = X(89)-JVS(829)*X(30)-JVS(830)*X(32)-JVS(831)*X(40)-JVS(832)*X(41)-JVS(833)*X(47)-JVS(834)*X(53)-JVS(835)*X(56)&
            &-JVS(836)*X(57)-JVS(837)*X(58)-JVS(838)*X(63)-JVS(839)*X(66)-JVS(840)*X(70)-JVS(841)*X(71)-JVS(842)*X(72)&
            &-JVS(843)*X(73)-JVS(844)*X(74)-JVS(845)*X(75)-JVS(846)*X(76)-JVS(847)*X(77)-JVS(848)*X(78)-JVS(849)*X(80)&
            &-JVS(850)*X(81)-JVS(851)*X(82)-JVS(852)*X(83)-JVS(853)*X(84)-JVS(854)*X(85)-JVS(855)*X(86)-JVS(856)*X(87)&
            &-JVS(857)*X(88)
  X(89) = X(89)/JVS(858)
  X(88) = (X(88)-JVS(828)*X(89))/(JVS(827))
  X(87) = (X(87)-JVS(788)*X(88)-JVS(789)*X(89))/(JVS(787))
  X(86) = (X(86)-JVS(757)*X(87)-JVS(758)*X(88)-JVS(759)*X(89))/(JVS(756))
  X(85) = (X(85)-JVS(697)*X(86)-JVS(698)*X(87)-JVS(699)*X(88)-JVS(700)*X(89))/(JVS(696))
  X(84) = (X(84)-JVS(681)*X(85)-JVS(682)*X(86)-JVS(683)*X(87)-JVS(684)*X(88)-JVS(685)*X(89))/(JVS(680))
  X(83) = (X(83)-JVS(607)*X(84)-JVS(608)*X(85)-JVS(609)*X(86)-JVS(610)*X(87)-JVS(611)*X(88)-JVS(612)*X(89))/(JVS(606))
  X(82) = (X(82)-JVS(588)*X(83)-JVS(589)*X(84)-JVS(590)*X(85)-JVS(591)*X(86)-JVS(592)*X(87)-JVS(593)*X(88)-JVS(594)&
            &*X(89))/(JVS(587))
  X(81) = (X(81)-JVS(558)*X(82)-JVS(559)*X(83)-JVS(560)*X(84)-JVS(561)*X(85)-JVS(562)*X(86)-JVS(563)*X(87)-JVS(564)&
            &*X(88)-JVS(565)*X(89))/(JVS(557))
  X(80) = (X(80)-JVS(532)*X(81)-JVS(533)*X(82)-JVS(534)*X(83)-JVS(535)*X(84)-JVS(536)*X(86)-JVS(537)*X(87)-JVS(538)&
            &*X(89))/(JVS(531))
  X(79) = (X(79)-JVS(519)*X(80)-JVS(520)*X(81)-JVS(521)*X(82)-JVS(522)*X(83)-JVS(523)*X(84)-JVS(524)*X(85)-JVS(525)&
            &*X(86)-JVS(526)*X(87)-JVS(527)*X(88)-JVS(528)*X(89))/(JVS(518))
  X(78) = (X(78)-JVS(480)*X(80)-JVS(481)*X(81)-JVS(482)*X(82)-JVS(483)*X(83)-JVS(484)*X(84)-JVS(485)*X(86)-JVS(486)&
            &*X(87)-JVS(487)*X(89))/(JVS(479))
  X(77) = (X(77)-JVS(468)*X(80)-JVS(469)*X(81)-JVS(470)*X(82)-JVS(471)*X(83)-JVS(472)*X(84)-JVS(473)*X(86)-JVS(474)&
            &*X(87)-JVS(475)*X(89))/(JVS(467))
  X(76) = (X(76)-JVS(452)*X(78)-JVS(453)*X(80)-JVS(454)*X(81)-JVS(455)*X(82)-JVS(456)*X(83)-JVS(457)*X(84)-JVS(458)&
            &*X(86)-JVS(459)*X(87)-JVS(460)*X(89))/(JVS(451))
  X(75) = (X(75)-JVS(434)*X(76)-JVS(435)*X(77)-JVS(436)*X(78)-JVS(437)*X(80)-JVS(438)*X(81)-JVS(439)*X(82)-JVS(440)&
            &*X(83)-JVS(441)*X(84)-JVS(442)*X(85)-JVS(443)*X(86)-JVS(444)*X(87)-JVS(445)*X(88)-JVS(446)*X(89))/(JVS(433))
  X(74) = (X(74)-JVS(416)*X(80)-JVS(417)*X(81)-JVS(418)*X(82)-JVS(419)*X(83)-JVS(420)*X(84)-JVS(421)*X(86)-JVS(422)&
            &*X(87)-JVS(423)*X(89))/(JVS(415))
  X(73) = (X(73)-JVS(397)*X(76)-JVS(398)*X(77)-JVS(399)*X(78)-JVS(400)*X(81)-JVS(401)*X(82)-JVS(402)*X(83)-JVS(403)&
            &*X(84)-JVS(404)*X(85)-JVS(405)*X(86)-JVS(406)*X(87)-JVS(407)*X(88)-JVS(408)*X(89))/(JVS(396))
  X(72) = (X(72)-JVS(385)*X(81)-JVS(386)*X(82)-JVS(387)*X(83)-JVS(388)*X(84)-JVS(389)*X(86)-JVS(390)*X(87))/(JVS(384))
  X(71) = (X(71)-JVS(374)*X(81)-JVS(375)*X(83)-JVS(376)*X(84)-JVS(377)*X(86)-JVS(378)*X(87))/(JVS(373))
  X(70) = (X(70)-JVS(364)*X(77)-JVS(365)*X(81)-JVS(366)*X(82)-JVS(367)*X(83)-JVS(368)*X(84)-JVS(369)*X(86)-JVS(370)&
            &*X(87))/(JVS(363))
  X(69) = (X(69)-JVS(347)*X(77)-JVS(348)*X(81)-JVS(349)*X(83)-JVS(350)*X(84))/(JVS(346))
  X(68) = (X(68)-JVS(338)*X(78)-JVS(339)*X(80)-JVS(340)*X(81)-JVS(341)*X(83)-JVS(342)*X(84)-JVS(343)*X(86)-JVS(344)&
            &*X(87))/(JVS(337))
  X(67) = (X(67)-JVS(332)*X(81)-JVS(333)*X(83)-JVS(334)*X(84))/(JVS(331))
  X(66) = (X(66)-JVS(327)*X(82)-JVS(328)*X(84)-JVS(329)*X(86)-JVS(330)*X(87))/(JVS(326))
  X(65) = (X(65)-JVS(306)*X(67)-JVS(307)*X(68)-JVS(308)*X(69)-JVS(309)*X(70)-JVS(310)*X(74)-JVS(311)*X(75)-JVS(312)&
            &*X(76)-JVS(313)*X(77)-JVS(314)*X(78)-JVS(315)*X(79)-JVS(316)*X(81)-JVS(317)*X(82)-JVS(318)*X(83)-JVS(319)*X(84)&
            &-JVS(320)*X(86)-JVS(321)*X(87)-JVS(322)*X(88)-JVS(323)*X(89))/(JVS(305))
  X(64) = (X(64)-JVS(296)*X(82)-JVS(297)*X(84)-JVS(298)*X(86)-JVS(299)*X(87))/(JVS(295))
  X(63) = (X(63)-JVS(287)*X(84)-JVS(288)*X(86)-JVS(289)*X(87))/(JVS(286))
  X(62) = (X(62)-JVS(276)*X(76)-JVS(277)*X(78)-JVS(278)*X(81)-JVS(279)*X(82)-JVS(280)*X(83)-JVS(281)*X(84)-JVS(282)&
            &*X(87)-JVS(283)*X(89))/(JVS(275))
  X(61) = (X(61)-JVS(269)*X(67)-JVS(270)*X(81)-JVS(271)*X(86)-JVS(272)*X(87))/(JVS(268))
  X(60) = (X(60)-JVS(262)*X(69)-JVS(263)*X(79)-JVS(264)*X(83)-JVS(265)*X(84)-JVS(266)*X(86)-JVS(267)*X(88))/(JVS(261))
  X(59) = (X(59)-JVS(256)*X(69)-JVS(257)*X(84)-JVS(258)*X(86)-JVS(259)*X(87))/(JVS(255))
  X(58) = (X(58)-JVS(247)*X(63)-JVS(248)*X(66)-JVS(249)*X(71)-JVS(250)*X(82)-JVS(251)*X(84)-JVS(252)*X(86)-JVS(253)&
            &*X(87))/(JVS(246))
  X(57) = (X(57)-JVS(238)*X(63)-JVS(239)*X(84)-JVS(240)*X(86)-JVS(241)*X(87))/(JVS(237))
  X(56) = (X(56)-JVS(231)*X(63)-JVS(232)*X(84)-JVS(233)*X(86))/(JVS(230))
  X(55) = (X(55)-JVS(226)*X(69)-JVS(227)*X(79)-JVS(228)*X(83)-JVS(229)*X(84))/(JVS(225))
  X(54) = (X(54)-JVS(219)*X(81)-JVS(220)*X(83)-JVS(221)*X(84))/(JVS(218))
  X(53) = (X(53)-JVS(214)*X(84)-JVS(215)*X(86)-JVS(216)*X(87)-JVS(217)*X(88))/(JVS(213))
  X(52) = (X(52)-JVS(207)*X(84)-JVS(208)*X(86)-JVS(209)*X(87))/(JVS(206))
  X(51) = (X(51)-JVS(197)*X(64)-JVS(198)*X(72)-JVS(199)*X(76)-JVS(200)*X(78)-JVS(201)*X(80)-JVS(202)*X(82)-JVS(203)&
            &*X(84))/(JVS(196))
  X(50) = (X(50)-JVS(193)*X(84)-JVS(194)*X(85)-JVS(195)*X(88))/(JVS(192))
  X(49) = (X(49)-JVS(189)*X(80)-JVS(190)*X(84)-JVS(191)*X(86))/(JVS(188))
  X(48) = (X(48)-JVS(184)*X(55)-JVS(185)*X(69)-JVS(186)*X(83)-JVS(187)*X(84))/(JVS(183))
  X(47) = (X(47)-JVS(180)*X(84)-JVS(181)*X(88)-JVS(182)*X(89))/(JVS(179))
  X(46) = (X(46)-JVS(176)*X(59)-JVS(177)*X(84)-JVS(178)*X(86))/(JVS(175))
  X(45) = (X(45)-JVS(169)*X(70)-JVS(170)*X(75)-JVS(171)*X(79)-JVS(172)*X(81)-JVS(173)*X(84)-JVS(174)*X(88))/(JVS(168))
  X(44) = (X(44)-JVS(163)*X(71)-JVS(164)*X(84)-JVS(165)*X(86))/(JVS(162))
  X(43) = (X(43)-JVS(156)*X(52)-JVS(157)*X(53)-JVS(158)*X(62)-JVS(159)*X(84)-JVS(160)*X(86)-JVS(161)*X(87))/(JVS(155))
  X(42) = (X(42)-JVS(146)*X(67)-JVS(147)*X(69)-JVS(148)*X(82)-JVS(149)*X(83)-JVS(150)*X(84)-JVS(151)*X(85)-JVS(152)&
            &*X(86)-JVS(153)*X(89))/(JVS(145))
  X(41) = (X(41)-JVS(140)*X(56)-JVS(141)*X(63)-JVS(142)*X(84)-JVS(143)*X(87))/(JVS(139))
  X(40) = (X(40)-JVS(135)*X(84)-JVS(136)*X(85)-JVS(137)*X(86)-JVS(138)*X(89))/(JVS(134))
  X(39) = (X(39)-JVS(131)*X(83)-JVS(132)*X(84)-JVS(133)*X(87))/(JVS(130))
  X(38) = (X(38)-JVS(126)*X(52)-JVS(127)*X(84)-JVS(128)*X(86))/(JVS(125))
  X(37) = (X(37)-JVS(120)*X(63)-JVS(121)*X(69)-JVS(122)*X(81)-JVS(123)*X(84)-JVS(124)*X(87))/(JVS(119))
  X(36) = (X(36)-JVS(117)*X(83)-JVS(118)*X(84))/(JVS(116))
  X(35) = (X(35)-JVS(109)*X(68)-JVS(110)*X(80)-JVS(111)*X(81)-JVS(112)*X(82)-JVS(113)*X(84)-JVS(114)*X(87)-JVS(115)&
            &*X(89))/(JVS(108))
  X(34) = (X(34)-JVS(105)*X(84)-JVS(106)*X(86)-JVS(107)*X(88))/(JVS(104))
  X(33) = (X(33)-JVS(101)*X(82)-JVS(102)*X(84)-JVS(103)*X(86))/(JVS(100))
  X(32) = (X(32)-JVS(97)*X(57)-JVS(98)*X(84)-JVS(99)*X(86))/(JVS(96))
  X(31) = (X(31)-JVS(94)*X(84)-JVS(95)*X(87))/(JVS(93))
  X(30) = (X(30)-JVS(89)*X(72)-JVS(90)*X(84)-JVS(91)*X(86))/(JVS(88))
  X(29) = (X(29)-JVS(85)*X(66)-JVS(86)*X(84)-JVS(87)*X(86))/(JVS(84))
  X(28) = (X(28)-JVS(81)*X(64)-JVS(82)*X(84)-JVS(83)*X(86))/(JVS(80))
  X(27) = (X(27)-JVS(76)*X(48)-JVS(77)*X(55)-JVS(78)*X(79)-JVS(79)*X(84))/(JVS(75))
  X(26) = (X(26)-JVS(72)*X(78)-JVS(73)*X(84)-JVS(74)*X(86))/(JVS(71))
  X(25) = (X(25)-JVS(68)*X(64)-JVS(69)*X(82)-JVS(70)*X(84))/(JVS(67))
  X(24) = (X(24)-JVS(65)*X(81)-JVS(66)*X(84))/(JVS(64))
  X(23) = (X(23)-JVS(61)*X(76)-JVS(62)*X(84)-JVS(63)*X(86))/(JVS(60))
  X(22) = (X(22)-JVS(58)*X(84)-JVS(59)*X(88))/(JVS(57))
  X(21) = (X(21)-JVS(54)*X(81)-JVS(55)*X(88))/(JVS(53))
  X(20) = (X(20)-JVS(51)*X(39)-JVS(52)*X(87))/(JVS(50))
  X(19) = (X(19)-JVS(48)*X(84)-JVS(49)*X(86))/(JVS(47))
  X(18) = (X(18)-JVS(46)*X(84))/(JVS(45))
  X(17) = (X(17)-JVS(43)*X(84))/(JVS(42))
  X(16) = (X(16)-JVS(39)*X(24)-JVS(40)*X(81)-JVS(41)*X(84))/(JVS(38))
  X(15) = (X(15)-JVS(37)*X(55))/(JVS(36))
  X(14) = (X(14)-JVS(35)*X(84))/(JVS(34))
  X(13) = (X(13)-JVS(33)*X(84))/(JVS(32))
  X(12) = (X(12)-JVS(31)*X(84))/(JVS(30))
  X(11) = (X(11)-JVS(29)*X(84))/(JVS(28))
  X(10) = (X(10)-JVS(27)*X(84))/(JVS(26))
  X(9) = (X(9)-JVS(25)*X(84))/(JVS(24))
  X(8) = (X(8)-JVS(22)*X(9)-JVS(23)*X(84))/(JVS(21))
  X(7) = (X(7)-JVS(20)*X(84))/(JVS(19))
  X(6) = (X(6)-JVS(17)*X(7)-JVS(18)*X(84))/(JVS(16))
  X(5) = (X(5)-JVS(14)*X(54)-JVS(15)*X(84))/(JVS(13))
  X(4) = (X(4)-JVS(11)*X(67)-JVS(12)*X(84))/(JVS(10))
  X(3) = (X(3)-JVS(8)*X(72)-JVS(9)*X(86))/(JVS(7))
  X(2) = (X(2)-JVS(5)*X(72)-JVS(6)*X(87))/(JVS(4))
  X(1) = (X(1)-JVS(2)*X(16)-JVS(3)*X(84))/(JVS(1))
END SUBROUTINE mozart_mosaic_4bin_vbs0_KppSolve
      SUBROUTINE mozart_mosaic_4bin_vbs0_WCOPY(N,X,incX,Y,incY)
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
      END SUBROUTINE mozart_mosaic_4bin_vbs0_WCOPY
      SUBROUTINE mozart_mosaic_4bin_vbs0_WAXPY(N,Alpha,X,incX,Y,incY)
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
      END SUBROUTINE mozart_mosaic_4bin_vbs0_WAXPY
      SUBROUTINE mozart_mosaic_4bin_vbs0_WSCAL(N,Alpha,X,incX)
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
      END SUBROUTINE mozart_mosaic_4bin_vbs0_WSCAL
      REAL(kind=dp) FUNCTION mozart_mosaic_4bin_vbs0_WLAMCH( C )
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
          CALL mozart_mosaic_4bin_vbs0_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10 Eps = Eps*2
        i = i-1
      END IF
      mozart_mosaic_4bin_vbs0_WLAMCH = Eps
      END FUNCTION mozart_mosaic_4bin_vbs0_WLAMCH
      SUBROUTINE mozart_mosaic_4bin_vbs0_WLAMCH_ADD( A, B, Sum )
      REAL(kind=dp) A, B, Sum
      Sum = A + B
      END SUBROUTINE mozart_mosaic_4bin_vbs0_WLAMCH_ADD
      SUBROUTINE mozart_mosaic_4bin_vbs0_SET2ZERO(N,Y)
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
      END SUBROUTINE mozart_mosaic_4bin_vbs0_SET2ZERO
      REAL(kind=dp) FUNCTION mozart_mosaic_4bin_vbs0_WDOT (N, DX, incX, DY, incY)
      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N)
      INTEGER :: i, IX, IY, M, MP1, NS
      mozart_mosaic_4bin_vbs0_WDOT = 0.0D0
      IF (N .LE. 0) RETURN
      IF (incX .EQ. incY) IF (incX-1) 5,20,60
    5 IX = 1
      IY = 1
      IF (incX .LT. 0) IX = (-N+1)*incX + 1
      IF (incY .LT. 0) IY = (-N+1)*incY + 1
      DO i = 1,N
        mozart_mosaic_4bin_vbs0_WDOT = mozart_mosaic_4bin_vbs0_WDOT + DX(IX)*DY(IY)
        IX = IX + incX
        IY = IY + incY
      END DO
      RETURN
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO i = 1,M
         mozart_mosaic_4bin_vbs0_WDOT = mozart_mosaic_4bin_vbs0_WDOT + DX(i)*DY(i)
      END DO
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO i = MP1,N,5
          mozart_mosaic_4bin_vbs0_WDOT = mozart_mosaic_4bin_vbs0_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) + &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)
      END DO
      RETURN
   60 NS = N*incX
      DO i = 1,NS,incX
        mozart_mosaic_4bin_vbs0_WDOT = mozart_mosaic_4bin_vbs0_WDOT + DX(i)*DY(i)
      END DO
      END FUNCTION mozart_mosaic_4bin_vbs0_WDOT
   SUBROUTINE decomp_mozart_mosaic_4bin_vbs0( JVS, IER )
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
   W( 16 ) = JVS( 2 )
   W( 84 ) = JVS( 3 )
  JVS( 1) = W( 1 )
  JVS( 2) = W( 16 )
  JVS( 3) = W( 84 )
  IF ( ABS( JVS( 4 )) < TINY(a) ) THEN
         IER = 2
         RETURN
  END IF
   W( 2 ) = JVS( 4 )
   W( 72 ) = JVS( 5 )
   W( 87 ) = JVS( 6 )
  JVS( 4) = W( 2 )
  JVS( 5) = W( 72 )
  JVS( 6) = W( 87 )
  IF ( ABS( JVS( 7 )) < TINY(a) ) THEN
         IER = 3
         RETURN
  END IF
   W( 3 ) = JVS( 7 )
   W( 72 ) = JVS( 8 )
   W( 86 ) = JVS( 9 )
  JVS( 7) = W( 3 )
  JVS( 8) = W( 72 )
  JVS( 9) = W( 86 )
  IF ( ABS( JVS( 10 )) < TINY(a) ) THEN
         IER = 4
         RETURN
  END IF
   W( 4 ) = JVS( 10 )
   W( 67 ) = JVS( 11 )
   W( 84 ) = JVS( 12 )
  JVS( 10) = W( 4 )
  JVS( 11) = W( 67 )
  JVS( 12) = W( 84 )
  IF ( ABS( JVS( 13 )) < TINY(a) ) THEN
         IER = 5
         RETURN
  END IF
   W( 5 ) = JVS( 13 )
   W( 54 ) = JVS( 14 )
   W( 84 ) = JVS( 15 )
  JVS( 13) = W( 5 )
  JVS( 14) = W( 54 )
  JVS( 15) = W( 84 )
  IF ( ABS( JVS( 16 )) < TINY(a) ) THEN
         IER = 6
         RETURN
  END IF
   W( 6 ) = JVS( 16 )
   W( 7 ) = JVS( 17 )
   W( 84 ) = JVS( 18 )
  JVS( 16) = W( 6 )
  JVS( 17) = W( 7 )
  JVS( 18) = W( 84 )
  IF ( ABS( JVS( 19 )) < TINY(a) ) THEN
         IER = 7
         RETURN
  END IF
   W( 7 ) = JVS( 19 )
   W( 84 ) = JVS( 20 )
  JVS( 19) = W( 7 )
  JVS( 20) = W( 84 )
  IF ( ABS( JVS( 21 )) < TINY(a) ) THEN
         IER = 8
         RETURN
  END IF
   W( 8 ) = JVS( 21 )
   W( 9 ) = JVS( 22 )
   W( 84 ) = JVS( 23 )
  JVS( 21) = W( 8 )
  JVS( 22) = W( 9 )
  JVS( 23) = W( 84 )
  IF ( ABS( JVS( 24 )) < TINY(a) ) THEN
         IER = 9
         RETURN
  END IF
   W( 9 ) = JVS( 24 )
   W( 84 ) = JVS( 25 )
  JVS( 24) = W( 9 )
  JVS( 25) = W( 84 )
  IF ( ABS( JVS( 26 )) < TINY(a) ) THEN
         IER = 10
         RETURN
  END IF
   W( 10 ) = JVS( 26 )
   W( 84 ) = JVS( 27 )
  JVS( 26) = W( 10 )
  JVS( 27) = W( 84 )
  IF ( ABS( JVS( 28 )) < TINY(a) ) THEN
         IER = 11
         RETURN
  END IF
   W( 11 ) = JVS( 28 )
   W( 84 ) = JVS( 29 )
  JVS( 28) = W( 11 )
  JVS( 29) = W( 84 )
  IF ( ABS( JVS( 30 )) < TINY(a) ) THEN
         IER = 12
         RETURN
  END IF
   W( 12 ) = JVS( 30 )
   W( 84 ) = JVS( 31 )
  JVS( 30) = W( 12 )
  JVS( 31) = W( 84 )
  IF ( ABS( JVS( 32 )) < TINY(a) ) THEN
         IER = 13
         RETURN
  END IF
   W( 13 ) = JVS( 32 )
   W( 84 ) = JVS( 33 )
  JVS( 32) = W( 13 )
  JVS( 33) = W( 84 )
  IF ( ABS( JVS( 34 )) < TINY(a) ) THEN
         IER = 14
         RETURN
  END IF
   W( 14 ) = JVS( 34 )
   W( 84 ) = JVS( 35 )
  JVS( 34) = W( 14 )
  JVS( 35) = W( 84 )
  IF ( ABS( JVS( 36 )) < TINY(a) ) THEN
         IER = 15
         RETURN
  END IF
   W( 15 ) = JVS( 36 )
   W( 55 ) = JVS( 37 )
  JVS( 36) = W( 15 )
  JVS( 37) = W( 55 )
  IF ( ABS( JVS( 38 )) < TINY(a) ) THEN
         IER = 16
         RETURN
  END IF
   W( 16 ) = JVS( 38 )
   W( 24 ) = JVS( 39 )
   W( 81 ) = JVS( 40 )
   W( 84 ) = JVS( 41 )
  JVS( 38) = W( 16 )
  JVS( 39) = W( 24 )
  JVS( 40) = W( 81 )
  JVS( 41) = W( 84 )
  IF ( ABS( JVS( 42 )) < TINY(a) ) THEN
         IER = 17
         RETURN
  END IF
   W( 17 ) = JVS( 42 )
   W( 84 ) = JVS( 43 )
  JVS( 42) = W( 17 )
  JVS( 43) = W( 84 )
  IF ( ABS( JVS( 45 )) < TINY(a) ) THEN
         IER = 18
         RETURN
  END IF
   W( 17 ) = JVS( 44 )
   W( 18 ) = JVS( 45 )
   W( 84 ) = JVS( 46 )
  a = -W( 17 ) / JVS( 42 )
  W( 17 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 43 )
  JVS( 44) = W( 17 )
  JVS( 45) = W( 18 )
  JVS( 46) = W( 84 )
  IF ( ABS( JVS( 47 )) < TINY(a) ) THEN
         IER = 19
         RETURN
  END IF
   W( 19 ) = JVS( 47 )
   W( 84 ) = JVS( 48 )
   W( 86 ) = JVS( 49 )
  JVS( 47) = W( 19 )
  JVS( 48) = W( 84 )
  JVS( 49) = W( 86 )
  IF ( ABS( JVS( 50 )) < TINY(a) ) THEN
         IER = 20
         RETURN
  END IF
   W( 20 ) = JVS( 50 )
   W( 39 ) = JVS( 51 )
   W( 87 ) = JVS( 52 )
  JVS( 50) = W( 20 )
  JVS( 51) = W( 39 )
  JVS( 52) = W( 87 )
  IF ( ABS( JVS( 53 )) < TINY(a) ) THEN
         IER = 21
         RETURN
  END IF
   W( 21 ) = JVS( 53 )
   W( 81 ) = JVS( 54 )
   W( 88 ) = JVS( 55 )
  JVS( 53) = W( 21 )
  JVS( 54) = W( 81 )
  JVS( 55) = W( 88 )
  IF ( ABS( JVS( 57 )) < TINY(a) ) THEN
         IER = 22
         RETURN
  END IF
   W( 18 ) = JVS( 56 )
   W( 22 ) = JVS( 57 )
   W( 84 ) = JVS( 58 )
   W( 88 ) = JVS( 59 )
  a = -W( 18 ) / JVS( 45 )
  W( 18 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 46 )
  JVS( 56) = W( 18 )
  JVS( 57) = W( 22 )
  JVS( 58) = W( 84 )
  JVS( 59) = W( 88 )
  IF ( ABS( JVS( 60 )) < TINY(a) ) THEN
         IER = 23
         RETURN
  END IF
   W( 23 ) = JVS( 60 )
   W( 76 ) = JVS( 61 )
   W( 84 ) = JVS( 62 )
   W( 86 ) = JVS( 63 )
  JVS( 60) = W( 23 )
  JVS( 61) = W( 76 )
  JVS( 62) = W( 84 )
  JVS( 63) = W( 86 )
  IF ( ABS( JVS( 64 )) < TINY(a) ) THEN
         IER = 24
         RETURN
  END IF
   W( 24 ) = JVS( 64 )
   W( 81 ) = JVS( 65 )
   W( 84 ) = JVS( 66 )
  JVS( 64) = W( 24 )
  JVS( 65) = W( 81 )
  JVS( 66) = W( 84 )
  IF ( ABS( JVS( 67 )) < TINY(a) ) THEN
         IER = 25
         RETURN
  END IF
   W( 25 ) = JVS( 67 )
   W( 64 ) = JVS( 68 )
   W( 82 ) = JVS( 69 )
   W( 84 ) = JVS( 70 )
  JVS( 67) = W( 25 )
  JVS( 68) = W( 64 )
  JVS( 69) = W( 82 )
  JVS( 70) = W( 84 )
  IF ( ABS( JVS( 71 )) < TINY(a) ) THEN
         IER = 26
         RETURN
  END IF
   W( 26 ) = JVS( 71 )
   W( 78 ) = JVS( 72 )
   W( 84 ) = JVS( 73 )
   W( 86 ) = JVS( 74 )
  JVS( 71) = W( 26 )
  JVS( 72) = W( 78 )
  JVS( 73) = W( 84 )
  JVS( 74) = W( 86 )
  IF ( ABS( JVS( 75 )) < TINY(a) ) THEN
         IER = 27
         RETURN
  END IF
   W( 27 ) = JVS( 75 )
   W( 48 ) = JVS( 76 )
   W( 55 ) = JVS( 77 )
   W( 79 ) = JVS( 78 )
   W( 84 ) = JVS( 79 )
  JVS( 75) = W( 27 )
  JVS( 76) = W( 48 )
  JVS( 77) = W( 55 )
  JVS( 78) = W( 79 )
  JVS( 79) = W( 84 )
  IF ( ABS( JVS( 80 )) < TINY(a) ) THEN
         IER = 28
         RETURN
  END IF
   W( 28 ) = JVS( 80 )
   W( 64 ) = JVS( 81 )
   W( 84 ) = JVS( 82 )
   W( 86 ) = JVS( 83 )
  JVS( 80) = W( 28 )
  JVS( 81) = W( 64 )
  JVS( 82) = W( 84 )
  JVS( 83) = W( 86 )
  IF ( ABS( JVS( 84 )) < TINY(a) ) THEN
         IER = 29
         RETURN
  END IF
   W( 29 ) = JVS( 84 )
   W( 66 ) = JVS( 85 )
   W( 84 ) = JVS( 86 )
   W( 86 ) = JVS( 87 )
  JVS( 84) = W( 29 )
  JVS( 85) = W( 66 )
  JVS( 86) = W( 84 )
  JVS( 87) = W( 86 )
  IF ( ABS( JVS( 88 )) < TINY(a) ) THEN
         IER = 30
         RETURN
  END IF
   W( 30 ) = JVS( 88 )
   W( 72 ) = JVS( 89 )
   W( 84 ) = JVS( 90 )
   W( 86 ) = JVS( 91 )
  JVS( 88) = W( 30 )
  JVS( 89) = W( 72 )
  JVS( 90) = W( 84 )
  JVS( 91) = W( 86 )
  IF ( ABS( JVS( 93 )) < TINY(a) ) THEN
         IER = 31
         RETURN
  END IF
   W( 13 ) = JVS( 92 )
   W( 31 ) = JVS( 93 )
   W( 84 ) = JVS( 94 )
   W( 87 ) = JVS( 95 )
  a = -W( 13 ) / JVS( 32 )
  W( 13 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 33 )
  JVS( 92) = W( 13 )
  JVS( 93) = W( 31 )
  JVS( 94) = W( 84 )
  JVS( 95) = W( 87 )
  IF ( ABS( JVS( 96 )) < TINY(a) ) THEN
         IER = 32
         RETURN
  END IF
   W( 32 ) = JVS( 96 )
   W( 57 ) = JVS( 97 )
   W( 84 ) = JVS( 98 )
   W( 86 ) = JVS( 99 )
  JVS( 96) = W( 32 )
  JVS( 97) = W( 57 )
  JVS( 98) = W( 84 )
  JVS( 99) = W( 86 )
  IF ( ABS( JVS( 100 )) < TINY(a) ) THEN
         IER = 33
         RETURN
  END IF
   W( 33 ) = JVS( 100 )
   W( 82 ) = JVS( 101 )
   W( 84 ) = JVS( 102 )
   W( 86 ) = JVS( 103 )
  JVS( 100) = W( 33 )
  JVS( 101) = W( 82 )
  JVS( 102) = W( 84 )
  JVS( 103) = W( 86 )
  IF ( ABS( JVS( 104 )) < TINY(a) ) THEN
         IER = 34
         RETURN
  END IF
   W( 34 ) = JVS( 104 )
   W( 84 ) = JVS( 105 )
   W( 86 ) = JVS( 106 )
   W( 88 ) = JVS( 107 )
  JVS( 104) = W( 34 )
  JVS( 105) = W( 84 )
  JVS( 106) = W( 86 )
  JVS( 107) = W( 88 )
  IF ( ABS( JVS( 108 )) < TINY(a) ) THEN
         IER = 35
         RETURN
  END IF
   W( 35 ) = JVS( 108 )
   W( 68 ) = JVS( 109 )
   W( 80 ) = JVS( 110 )
   W( 81 ) = JVS( 111 )
   W( 82 ) = JVS( 112 )
   W( 84 ) = JVS( 113 )
   W( 87 ) = JVS( 114 )
   W( 89 ) = JVS( 115 )
  JVS( 108) = W( 35 )
  JVS( 109) = W( 68 )
  JVS( 110) = W( 80 )
  JVS( 111) = W( 81 )
  JVS( 112) = W( 82 )
  JVS( 113) = W( 84 )
  JVS( 114) = W( 87 )
  JVS( 115) = W( 89 )
  IF ( ABS( JVS( 116 )) < TINY(a) ) THEN
         IER = 36
         RETURN
  END IF
   W( 36 ) = JVS( 116 )
   W( 83 ) = JVS( 117 )
   W( 84 ) = JVS( 118 )
  JVS( 116) = W( 36 )
  JVS( 117) = W( 83 )
  JVS( 118) = W( 84 )
  IF ( ABS( JVS( 119 )) < TINY(a) ) THEN
         IER = 37
         RETURN
  END IF
   W( 37 ) = JVS( 119 )
   W( 63 ) = JVS( 120 )
   W( 69 ) = JVS( 121 )
   W( 81 ) = JVS( 122 )
   W( 84 ) = JVS( 123 )
   W( 87 ) = JVS( 124 )
  JVS( 119) = W( 37 )
  JVS( 120) = W( 63 )
  JVS( 121) = W( 69 )
  JVS( 122) = W( 81 )
  JVS( 123) = W( 84 )
  JVS( 124) = W( 87 )
  IF ( ABS( JVS( 125 )) < TINY(a) ) THEN
         IER = 38
         RETURN
  END IF
   W( 38 ) = JVS( 125 )
   W( 52 ) = JVS( 126 )
   W( 84 ) = JVS( 127 )
   W( 86 ) = JVS( 128 )
  JVS( 125) = W( 38 )
  JVS( 126) = W( 52 )
  JVS( 127) = W( 84 )
  JVS( 128) = W( 86 )
  IF ( ABS( JVS( 130 )) < TINY(a) ) THEN
         IER = 39
         RETURN
  END IF
   W( 36 ) = JVS( 129 )
   W( 39 ) = JVS( 130 )
   W( 83 ) = JVS( 131 )
   W( 84 ) = JVS( 132 )
   W( 87 ) = JVS( 133 )
  a = -W( 36 ) / JVS( 116 )
  W( 36 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 117 )
  W( 84 ) = W( 84 ) + a*JVS( 118 )
  JVS( 129) = W( 36 )
  JVS( 130) = W( 39 )
  JVS( 131) = W( 83 )
  JVS( 132) = W( 84 )
  JVS( 133) = W( 87 )
  IF ( ABS( JVS( 134 )) < TINY(a) ) THEN
         IER = 40
         RETURN
  END IF
   W( 40 ) = JVS( 134 )
   W( 84 ) = JVS( 135 )
   W( 85 ) = JVS( 136 )
   W( 86 ) = JVS( 137 )
   W( 89 ) = JVS( 138 )
  JVS( 134) = W( 40 )
  JVS( 135) = W( 84 )
  JVS( 136) = W( 85 )
  JVS( 137) = W( 86 )
  JVS( 138) = W( 89 )
  IF ( ABS( JVS( 139 )) < TINY(a) ) THEN
         IER = 41
         RETURN
  END IF
   W( 41 ) = JVS( 139 )
   W( 56 ) = JVS( 140 )
   W( 63 ) = JVS( 141 )
   W( 84 ) = JVS( 142 )
   W( 87 ) = JVS( 143 )
  JVS( 139) = W( 41 )
  JVS( 140) = W( 56 )
  JVS( 141) = W( 63 )
  JVS( 142) = W( 84 )
  JVS( 143) = W( 87 )
  IF ( ABS( JVS( 145 )) < TINY(a) ) THEN
         IER = 42
         RETURN
  END IF
   W( 36 ) = JVS( 144 )
   W( 42 ) = JVS( 145 )
   W( 67 ) = JVS( 146 )
   W( 69 ) = JVS( 147 )
   W( 82 ) = JVS( 148 )
   W( 83 ) = JVS( 149 )
   W( 84 ) = JVS( 150 )
   W( 85 ) = JVS( 151 )
   W( 86 ) = JVS( 152 )
   W( 89 ) = JVS( 153 )
  a = -W( 36 ) / JVS( 116 )
  W( 36 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 117 )
  W( 84 ) = W( 84 ) + a*JVS( 118 )
  JVS( 144) = W( 36 )
  JVS( 145) = W( 42 )
  JVS( 146) = W( 67 )
  JVS( 147) = W( 69 )
  JVS( 148) = W( 82 )
  JVS( 149) = W( 83 )
  JVS( 150) = W( 84 )
  JVS( 151) = W( 85 )
  JVS( 152) = W( 86 )
  JVS( 153) = W( 89 )
  IF ( ABS( JVS( 155 )) < TINY(a) ) THEN
         IER = 43
         RETURN
  END IF
   W( 38 ) = JVS( 154 )
   W( 43 ) = JVS( 155 )
   W( 52 ) = JVS( 156 )
   W( 53 ) = JVS( 157 )
   W( 62 ) = JVS( 158 )
   W( 84 ) = JVS( 159 )
   W( 86 ) = JVS( 160 )
   W( 87 ) = JVS( 161 )
  a = -W( 38 ) / JVS( 125 )
  W( 38 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 126 )
  W( 84 ) = W( 84 ) + a*JVS( 127 )
  W( 86 ) = W( 86 ) + a*JVS( 128 )
  JVS( 154) = W( 38 )
  JVS( 155) = W( 43 )
  JVS( 156) = W( 52 )
  JVS( 157) = W( 53 )
  JVS( 158) = W( 62 )
  JVS( 159) = W( 84 )
  JVS( 160) = W( 86 )
  JVS( 161) = W( 87 )
  IF ( ABS( JVS( 162 )) < TINY(a) ) THEN
         IER = 44
         RETURN
  END IF
   W( 44 ) = JVS( 162 )
   W( 71 ) = JVS( 163 )
   W( 84 ) = JVS( 164 )
   W( 86 ) = JVS( 165 )
  JVS( 162) = W( 44 )
  JVS( 163) = W( 71 )
  JVS( 164) = W( 84 )
  JVS( 165) = W( 86 )
  IF ( ABS( JVS( 168 )) < TINY(a) ) THEN
         IER = 45
         RETURN
  END IF
   W( 21 ) = JVS( 166 )
   W( 24 ) = JVS( 167 )
   W( 45 ) = JVS( 168 )
   W( 70 ) = JVS( 169 )
   W( 75 ) = JVS( 170 )
   W( 79 ) = JVS( 171 )
   W( 81 ) = JVS( 172 )
   W( 84 ) = JVS( 173 )
   W( 88 ) = JVS( 174 )
  a = -W( 21 ) / JVS( 53 )
  W( 21 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 54 )
  W( 88 ) = W( 88 ) + a*JVS( 55 )
  a = -W( 24 ) / JVS( 64 )
  W( 24 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 65 )
  W( 84 ) = W( 84 ) + a*JVS( 66 )
  JVS( 166) = W( 21 )
  JVS( 167) = W( 24 )
  JVS( 168) = W( 45 )
  JVS( 169) = W( 70 )
  JVS( 170) = W( 75 )
  JVS( 171) = W( 79 )
  JVS( 172) = W( 81 )
  JVS( 173) = W( 84 )
  JVS( 174) = W( 88 )
  IF ( ABS( JVS( 175 )) < TINY(a) ) THEN
         IER = 46
         RETURN
  END IF
   W( 46 ) = JVS( 175 )
   W( 59 ) = JVS( 176 )
   W( 84 ) = JVS( 177 )
   W( 86 ) = JVS( 178 )
  JVS( 175) = W( 46 )
  JVS( 176) = W( 59 )
  JVS( 177) = W( 84 )
  JVS( 178) = W( 86 )
  IF ( ABS( JVS( 179 )) < TINY(a) ) THEN
         IER = 47
         RETURN
  END IF
   W( 47 ) = JVS( 179 )
   W( 84 ) = JVS( 180 )
   W( 88 ) = JVS( 181 )
   W( 89 ) = JVS( 182 )
  JVS( 179) = W( 47 )
  JVS( 180) = W( 84 )
  JVS( 181) = W( 88 )
  JVS( 182) = W( 89 )
  IF ( ABS( JVS( 183 )) < TINY(a) ) THEN
         IER = 48
         RETURN
  END IF
   W( 48 ) = JVS( 183 )
   W( 55 ) = JVS( 184 )
   W( 69 ) = JVS( 185 )
   W( 83 ) = JVS( 186 )
   W( 84 ) = JVS( 187 )
  JVS( 183) = W( 48 )
  JVS( 184) = W( 55 )
  JVS( 185) = W( 69 )
  JVS( 186) = W( 83 )
  JVS( 187) = W( 84 )
  IF ( ABS( JVS( 188 )) < TINY(a) ) THEN
         IER = 49
         RETURN
  END IF
   W( 49 ) = JVS( 188 )
   W( 80 ) = JVS( 189 )
   W( 84 ) = JVS( 190 )
   W( 86 ) = JVS( 191 )
  JVS( 188) = W( 49 )
  JVS( 189) = W( 80 )
  JVS( 190) = W( 84 )
  JVS( 191) = W( 86 )
  IF ( ABS( JVS( 192 )) < TINY(a) ) THEN
         IER = 50
         RETURN
  END IF
   W( 50 ) = JVS( 192 )
   W( 84 ) = JVS( 193 )
   W( 85 ) = JVS( 194 )
   W( 88 ) = JVS( 195 )
  JVS( 192) = W( 50 )
  JVS( 193) = W( 84 )
  JVS( 194) = W( 85 )
  JVS( 195) = W( 88 )
  IF ( ABS( JVS( 196 )) < TINY(a) ) THEN
         IER = 51
         RETURN
  END IF
   W( 51 ) = JVS( 196 )
   W( 64 ) = JVS( 197 )
   W( 72 ) = JVS( 198 )
   W( 76 ) = JVS( 199 )
   W( 78 ) = JVS( 200 )
   W( 80 ) = JVS( 201 )
   W( 82 ) = JVS( 202 )
   W( 84 ) = JVS( 203 )
  JVS( 196) = W( 51 )
  JVS( 197) = W( 64 )
  JVS( 198) = W( 72 )
  JVS( 199) = W( 76 )
  JVS( 200) = W( 78 )
  JVS( 201) = W( 80 )
  JVS( 202) = W( 82 )
  JVS( 203) = W( 84 )
  IF ( ABS( JVS( 206 )) < TINY(a) ) THEN
         IER = 52
         RETURN
  END IF
   W( 17 ) = JVS( 204 )
   W( 38 ) = JVS( 205 )
   W( 52 ) = JVS( 206 )
   W( 84 ) = JVS( 207 )
   W( 86 ) = JVS( 208 )
   W( 87 ) = JVS( 209 )
  a = -W( 17 ) / JVS( 42 )
  W( 17 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 43 )
  a = -W( 38 ) / JVS( 125 )
  W( 38 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 126 )
  W( 84 ) = W( 84 ) + a*JVS( 127 )
  W( 86 ) = W( 86 ) + a*JVS( 128 )
  JVS( 204) = W( 17 )
  JVS( 205) = W( 38 )
  JVS( 206) = W( 52 )
  JVS( 207) = W( 84 )
  JVS( 208) = W( 86 )
  JVS( 209) = W( 87 )
  IF ( ABS( JVS( 213 )) < TINY(a) ) THEN
         IER = 53
         RETURN
  END IF
   W( 22 ) = JVS( 210 )
   W( 38 ) = JVS( 211 )
   W( 52 ) = JVS( 212 )
   W( 53 ) = JVS( 213 )
   W( 84 ) = JVS( 214 )
   W( 86 ) = JVS( 215 )
   W( 87 ) = JVS( 216 )
   W( 88 ) = JVS( 217 )
  a = -W( 22 ) / JVS( 57 )
  W( 22 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 58 )
  W( 88 ) = W( 88 ) + a*JVS( 59 )
  a = -W( 38 ) / JVS( 125 )
  W( 38 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 126 )
  W( 84 ) = W( 84 ) + a*JVS( 127 )
  W( 86 ) = W( 86 ) + a*JVS( 128 )
  a = -W( 52 ) / JVS( 206 )
  W( 52 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 207 )
  W( 86 ) = W( 86 ) + a*JVS( 208 )
  W( 87 ) = W( 87 ) + a*JVS( 209 )
  JVS( 210) = W( 22 )
  JVS( 211) = W( 38 )
  JVS( 212) = W( 52 )
  JVS( 213) = W( 53 )
  JVS( 214) = W( 84 )
  JVS( 215) = W( 86 )
  JVS( 216) = W( 87 )
  JVS( 217) = W( 88 )
  IF ( ABS( JVS( 218 )) < TINY(a) ) THEN
         IER = 54
         RETURN
  END IF
   W( 54 ) = JVS( 218 )
   W( 81 ) = JVS( 219 )
   W( 83 ) = JVS( 220 )
   W( 84 ) = JVS( 221 )
  JVS( 218) = W( 54 )
  JVS( 219) = W( 81 )
  JVS( 220) = W( 83 )
  JVS( 221) = W( 84 )
  IF ( ABS( JVS( 225 )) < TINY(a) ) THEN
         IER = 55
         RETURN
  END IF
   W( 15 ) = JVS( 222 )
   W( 27 ) = JVS( 223 )
   W( 48 ) = JVS( 224 )
   W( 55 ) = JVS( 225 )
   W( 69 ) = JVS( 226 )
   W( 79 ) = JVS( 227 )
   W( 83 ) = JVS( 228 )
   W( 84 ) = JVS( 229 )
  a = -W( 15 ) / JVS( 36 )
  W( 15 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 37 )
  a = -W( 27 ) / JVS( 75 )
  W( 27 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 76 )
  W( 55 ) = W( 55 ) + a*JVS( 77 )
  W( 79 ) = W( 79 ) + a*JVS( 78 )
  W( 84 ) = W( 84 ) + a*JVS( 79 )
  a = -W( 48 ) / JVS( 183 )
  W( 48 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 184 )
  W( 69 ) = W( 69 ) + a*JVS( 185 )
  W( 83 ) = W( 83 ) + a*JVS( 186 )
  W( 84 ) = W( 84 ) + a*JVS( 187 )
  JVS( 222) = W( 15 )
  JVS( 223) = W( 27 )
  JVS( 224) = W( 48 )
  JVS( 225) = W( 55 )
  JVS( 226) = W( 69 )
  JVS( 227) = W( 79 )
  JVS( 228) = W( 83 )
  JVS( 229) = W( 84 )
  IF ( ABS( JVS( 230 )) < TINY(a) ) THEN
         IER = 56
         RETURN
  END IF
   W( 56 ) = JVS( 230 )
   W( 63 ) = JVS( 231 )
   W( 84 ) = JVS( 232 )
   W( 86 ) = JVS( 233 )
  JVS( 230) = W( 56 )
  JVS( 231) = W( 63 )
  JVS( 232) = W( 84 )
  JVS( 233) = W( 86 )
  IF ( ABS( JVS( 237 )) < TINY(a) ) THEN
         IER = 57
         RETURN
  END IF
   W( 32 ) = JVS( 234 )
   W( 41 ) = JVS( 235 )
   W( 56 ) = JVS( 236 )
   W( 57 ) = JVS( 237 )
   W( 63 ) = JVS( 238 )
   W( 84 ) = JVS( 239 )
   W( 86 ) = JVS( 240 )
   W( 87 ) = JVS( 241 )
  a = -W( 32 ) / JVS( 96 )
  W( 32 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 97 )
  W( 84 ) = W( 84 ) + a*JVS( 98 )
  W( 86 ) = W( 86 ) + a*JVS( 99 )
  a = -W( 41 ) / JVS( 139 )
  W( 41 ) = -a
  W( 56 ) = W( 56 ) + a*JVS( 140 )
  W( 63 ) = W( 63 ) + a*JVS( 141 )
  W( 84 ) = W( 84 ) + a*JVS( 142 )
  W( 87 ) = W( 87 ) + a*JVS( 143 )
  a = -W( 56 ) / JVS( 230 )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 231 )
  W( 84 ) = W( 84 ) + a*JVS( 232 )
  W( 86 ) = W( 86 ) + a*JVS( 233 )
  JVS( 234) = W( 32 )
  JVS( 235) = W( 41 )
  JVS( 236) = W( 56 )
  JVS( 237) = W( 57 )
  JVS( 238) = W( 63 )
  JVS( 239) = W( 84 )
  JVS( 240) = W( 86 )
  JVS( 241) = W( 87 )
  IF ( ABS( JVS( 246 )) < TINY(a) ) THEN
         IER = 58
         RETURN
  END IF
   W( 29 ) = JVS( 242 )
   W( 31 ) = JVS( 243 )
   W( 44 ) = JVS( 244 )
   W( 56 ) = JVS( 245 )
   W( 58 ) = JVS( 246 )
   W( 63 ) = JVS( 247 )
   W( 66 ) = JVS( 248 )
   W( 71 ) = JVS( 249 )
   W( 82 ) = JVS( 250 )
   W( 84 ) = JVS( 251 )
   W( 86 ) = JVS( 252 )
   W( 87 ) = JVS( 253 )
  a = -W( 29 ) / JVS( 84 )
  W( 29 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 85 )
  W( 84 ) = W( 84 ) + a*JVS( 86 )
  W( 86 ) = W( 86 ) + a*JVS( 87 )
  a = -W( 31 ) / JVS( 93 )
  W( 31 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 94 )
  W( 87 ) = W( 87 ) + a*JVS( 95 )
  a = -W( 44 ) / JVS( 162 )
  W( 44 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 163 )
  W( 84 ) = W( 84 ) + a*JVS( 164 )
  W( 86 ) = W( 86 ) + a*JVS( 165 )
  a = -W( 56 ) / JVS( 230 )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 231 )
  W( 84 ) = W( 84 ) + a*JVS( 232 )
  W( 86 ) = W( 86 ) + a*JVS( 233 )
  JVS( 242) = W( 29 )
  JVS( 243) = W( 31 )
  JVS( 244) = W( 44 )
  JVS( 245) = W( 56 )
  JVS( 246) = W( 58 )
  JVS( 247) = W( 63 )
  JVS( 248) = W( 66 )
  JVS( 249) = W( 71 )
  JVS( 250) = W( 82 )
  JVS( 251) = W( 84 )
  JVS( 252) = W( 86 )
  JVS( 253) = W( 87 )
  IF ( ABS( JVS( 255 )) < TINY(a) ) THEN
         IER = 59
         RETURN
  END IF
   W( 46 ) = JVS( 254 )
   W( 59 ) = JVS( 255 )
   W( 69 ) = JVS( 256 )
   W( 84 ) = JVS( 257 )
   W( 86 ) = JVS( 258 )
   W( 87 ) = JVS( 259 )
  a = -W( 46 ) / JVS( 175 )
  W( 46 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 176 )
  W( 84 ) = W( 84 ) + a*JVS( 177 )
  W( 86 ) = W( 86 ) + a*JVS( 178 )
  JVS( 254) = W( 46 )
  JVS( 255) = W( 59 )
  JVS( 256) = W( 69 )
  JVS( 257) = W( 84 )
  JVS( 258) = W( 86 )
  JVS( 259) = W( 87 )
  IF ( ABS( JVS( 261 )) < TINY(a) ) THEN
         IER = 60
         RETURN
  END IF
   W( 55 ) = JVS( 260 )
   W( 60 ) = JVS( 261 )
   W( 69 ) = JVS( 262 )
   W( 79 ) = JVS( 263 )
   W( 83 ) = JVS( 264 )
   W( 84 ) = JVS( 265 )
   W( 86 ) = JVS( 266 )
   W( 88 ) = JVS( 267 )
  a = -W( 55 ) / JVS( 225 )
  W( 55 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 226 )
  W( 79 ) = W( 79 ) + a*JVS( 227 )
  W( 83 ) = W( 83 ) + a*JVS( 228 )
  W( 84 ) = W( 84 ) + a*JVS( 229 )
  JVS( 260) = W( 55 )
  JVS( 261) = W( 60 )
  JVS( 262) = W( 69 )
  JVS( 263) = W( 79 )
  JVS( 264) = W( 83 )
  JVS( 265) = W( 84 )
  JVS( 266) = W( 86 )
  JVS( 267) = W( 88 )
  IF ( ABS( JVS( 268 )) < TINY(a) ) THEN
         IER = 61
         RETURN
  END IF
   W( 61 ) = JVS( 268 )
   W( 67 ) = JVS( 269 )
   W( 81 ) = JVS( 270 )
   W( 86 ) = JVS( 271 )
   W( 87 ) = JVS( 272 )
  JVS( 268) = W( 61 )
  JVS( 269) = W( 67 )
  JVS( 270) = W( 81 )
  JVS( 271) = W( 86 )
  JVS( 272) = W( 87 )
  IF ( ABS( JVS( 275 )) < TINY(a) ) THEN
         IER = 62
         RETURN
  END IF
   W( 20 ) = JVS( 273 )
   W( 39 ) = JVS( 274 )
   W( 62 ) = JVS( 275 )
   W( 76 ) = JVS( 276 )
   W( 78 ) = JVS( 277 )
   W( 81 ) = JVS( 278 )
   W( 82 ) = JVS( 279 )
   W( 83 ) = JVS( 280 )
   W( 84 ) = JVS( 281 )
   W( 87 ) = JVS( 282 )
   W( 89 ) = JVS( 283 )
  a = -W( 20 ) / JVS( 50 )
  W( 20 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 51 )
  W( 87 ) = W( 87 ) + a*JVS( 52 )
  a = -W( 39 ) / JVS( 130 )
  W( 39 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 131 )
  W( 84 ) = W( 84 ) + a*JVS( 132 )
  W( 87 ) = W( 87 ) + a*JVS( 133 )
  JVS( 273) = W( 20 )
  JVS( 274) = W( 39 )
  JVS( 275) = W( 62 )
  JVS( 276) = W( 76 )
  JVS( 277) = W( 78 )
  JVS( 278) = W( 81 )
  JVS( 279) = W( 82 )
  JVS( 280) = W( 83 )
  JVS( 281) = W( 84 )
  JVS( 282) = W( 87 )
  JVS( 283) = W( 89 )
  IF ( ABS( JVS( 286 )) < TINY(a) ) THEN
         IER = 63
         RETURN
  END IF
   W( 14 ) = JVS( 284 )
   W( 56 ) = JVS( 285 )
   W( 63 ) = JVS( 286 )
   W( 84 ) = JVS( 287 )
   W( 86 ) = JVS( 288 )
   W( 87 ) = JVS( 289 )
  a = -W( 14 ) / JVS( 34 )
  W( 14 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 35 )
  a = -W( 56 ) / JVS( 230 )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 231 )
  W( 84 ) = W( 84 ) + a*JVS( 232 )
  W( 86 ) = W( 86 ) + a*JVS( 233 )
  JVS( 284) = W( 14 )
  JVS( 285) = W( 56 )
  JVS( 286) = W( 63 )
  JVS( 287) = W( 84 )
  JVS( 288) = W( 86 )
  JVS( 289) = W( 87 )
  IF ( ABS( JVS( 295 )) < TINY(a) ) THEN
         IER = 64
         RETURN
  END IF
   W( 11 ) = JVS( 290 )
   W( 28 ) = JVS( 291 )
   W( 41 ) = JVS( 292 )
   W( 56 ) = JVS( 293 )
   W( 63 ) = JVS( 294 )
   W( 64 ) = JVS( 295 )
   W( 82 ) = JVS( 296 )
   W( 84 ) = JVS( 297 )
   W( 86 ) = JVS( 298 )
   W( 87 ) = JVS( 299 )
  a = -W( 11 ) / JVS( 28 )
  W( 11 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 29 )
  a = -W( 28 ) / JVS( 80 )
  W( 28 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 81 )
  W( 84 ) = W( 84 ) + a*JVS( 82 )
  W( 86 ) = W( 86 ) + a*JVS( 83 )
  a = -W( 41 ) / JVS( 139 )
  W( 41 ) = -a
  W( 56 ) = W( 56 ) + a*JVS( 140 )
  W( 63 ) = W( 63 ) + a*JVS( 141 )
  W( 84 ) = W( 84 ) + a*JVS( 142 )
  W( 87 ) = W( 87 ) + a*JVS( 143 )
  a = -W( 56 ) / JVS( 230 )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 231 )
  W( 84 ) = W( 84 ) + a*JVS( 232 )
  W( 86 ) = W( 86 ) + a*JVS( 233 )
  a = -W( 63 ) / JVS( 286 )
  W( 63 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 287 )
  W( 86 ) = W( 86 ) + a*JVS( 288 )
  W( 87 ) = W( 87 ) + a*JVS( 289 )
  JVS( 290) = W( 11 )
  JVS( 291) = W( 28 )
  JVS( 292) = W( 41 )
  JVS( 293) = W( 56 )
  JVS( 294) = W( 63 )
  JVS( 295) = W( 64 )
  JVS( 296) = W( 82 )
  JVS( 297) = W( 84 )
  JVS( 298) = W( 86 )
  JVS( 299) = W( 87 )
  IF ( ABS( JVS( 305 )) < TINY(a) ) THEN
         IER = 65
         RETURN
  END IF
   W( 36 ) = JVS( 300 )
   W( 43 ) = JVS( 301 )
   W( 52 ) = JVS( 302 )
   W( 53 ) = JVS( 303 )
   W( 62 ) = JVS( 304 )
   W( 65 ) = JVS( 305 )
   W( 67 ) = JVS( 306 )
   W( 68 ) = JVS( 307 )
   W( 69 ) = JVS( 308 )
   W( 70 ) = JVS( 309 )
   W( 74 ) = JVS( 310 )
   W( 75 ) = JVS( 311 )
   W( 76 ) = JVS( 312 )
   W( 77 ) = JVS( 313 )
   W( 78 ) = JVS( 314 )
   W( 79 ) = JVS( 315 )
   W( 81 ) = JVS( 316 )
   W( 82 ) = JVS( 317 )
   W( 83 ) = JVS( 318 )
   W( 84 ) = JVS( 319 )
   W( 86 ) = JVS( 320 )
   W( 87 ) = JVS( 321 )
   W( 88 ) = JVS( 322 )
   W( 89 ) = JVS( 323 )
  a = -W( 36 ) / JVS( 116 )
  W( 36 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 117 )
  W( 84 ) = W( 84 ) + a*JVS( 118 )
  a = -W( 43 ) / JVS( 155 )
  W( 43 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 156 )
  W( 53 ) = W( 53 ) + a*JVS( 157 )
  W( 62 ) = W( 62 ) + a*JVS( 158 )
  W( 84 ) = W( 84 ) + a*JVS( 159 )
  W( 86 ) = W( 86 ) + a*JVS( 160 )
  W( 87 ) = W( 87 ) + a*JVS( 161 )
  a = -W( 52 ) / JVS( 206 )
  W( 52 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 207 )
  W( 86 ) = W( 86 ) + a*JVS( 208 )
  W( 87 ) = W( 87 ) + a*JVS( 209 )
  a = -W( 53 ) / JVS( 213 )
  W( 53 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 214 )
  W( 86 ) = W( 86 ) + a*JVS( 215 )
  W( 87 ) = W( 87 ) + a*JVS( 216 )
  W( 88 ) = W( 88 ) + a*JVS( 217 )
  a = -W( 62 ) / JVS( 275 )
  W( 62 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 276 )
  W( 78 ) = W( 78 ) + a*JVS( 277 )
  W( 81 ) = W( 81 ) + a*JVS( 278 )
  W( 82 ) = W( 82 ) + a*JVS( 279 )
  W( 83 ) = W( 83 ) + a*JVS( 280 )
  W( 84 ) = W( 84 ) + a*JVS( 281 )
  W( 87 ) = W( 87 ) + a*JVS( 282 )
  W( 89 ) = W( 89 ) + a*JVS( 283 )
  JVS( 300) = W( 36 )
  JVS( 301) = W( 43 )
  JVS( 302) = W( 52 )
  JVS( 303) = W( 53 )
  JVS( 304) = W( 62 )
  JVS( 305) = W( 65 )
  JVS( 306) = W( 67 )
  JVS( 307) = W( 68 )
  JVS( 308) = W( 69 )
  JVS( 309) = W( 70 )
  JVS( 310) = W( 74 )
  JVS( 311) = W( 75 )
  JVS( 312) = W( 76 )
  JVS( 313) = W( 77 )
  JVS( 314) = W( 78 )
  JVS( 315) = W( 79 )
  JVS( 316) = W( 81 )
  JVS( 317) = W( 82 )
  JVS( 318) = W( 83 )
  JVS( 319) = W( 84 )
  JVS( 320) = W( 86 )
  JVS( 321) = W( 87 )
  JVS( 322) = W( 88 )
  JVS( 323) = W( 89 )
  IF ( ABS( JVS( 326 )) < TINY(a) ) THEN
         IER = 66
         RETURN
  END IF
   W( 12 ) = JVS( 324 )
   W( 29 ) = JVS( 325 )
   W( 66 ) = JVS( 326 )
   W( 82 ) = JVS( 327 )
   W( 84 ) = JVS( 328 )
   W( 86 ) = JVS( 329 )
   W( 87 ) = JVS( 330 )
  a = -W( 12 ) / JVS( 30 )
  W( 12 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 31 )
  a = -W( 29 ) / JVS( 84 )
  W( 29 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 85 )
  W( 84 ) = W( 84 ) + a*JVS( 86 )
  W( 86 ) = W( 86 ) + a*JVS( 87 )
  JVS( 324) = W( 12 )
  JVS( 325) = W( 29 )
  JVS( 326) = W( 66 )
  JVS( 327) = W( 82 )
  JVS( 328) = W( 84 )
  JVS( 329) = W( 86 )
  JVS( 330) = W( 87 )
  IF ( ABS( JVS( 331 )) < TINY(a) ) THEN
         IER = 67
         RETURN
  END IF
   W( 67 ) = JVS( 331 )
   W( 81 ) = JVS( 332 )
   W( 83 ) = JVS( 333 )
   W( 84 ) = JVS( 334 )
  JVS( 331) = W( 67 )
  JVS( 332) = W( 81 )
  JVS( 333) = W( 83 )
  JVS( 334) = W( 84 )
  IF ( ABS( JVS( 337 )) < TINY(a) ) THEN
         IER = 68
         RETURN
  END IF
   W( 61 ) = JVS( 335 )
   W( 67 ) = JVS( 336 )
   W( 68 ) = JVS( 337 )
   W( 78 ) = JVS( 338 )
   W( 80 ) = JVS( 339 )
   W( 81 ) = JVS( 340 )
   W( 83 ) = JVS( 341 )
   W( 84 ) = JVS( 342 )
   W( 86 ) = JVS( 343 )
   W( 87 ) = JVS( 344 )
  a = -W( 61 ) / JVS( 268 )
  W( 61 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 269 )
  W( 81 ) = W( 81 ) + a*JVS( 270 )
  W( 86 ) = W( 86 ) + a*JVS( 271 )
  W( 87 ) = W( 87 ) + a*JVS( 272 )
  a = -W( 67 ) / JVS( 331 )
  W( 67 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 332 )
  W( 83 ) = W( 83 ) + a*JVS( 333 )
  W( 84 ) = W( 84 ) + a*JVS( 334 )
  JVS( 335) = W( 61 )
  JVS( 336) = W( 67 )
  JVS( 337) = W( 68 )
  JVS( 338) = W( 78 )
  JVS( 339) = W( 80 )
  JVS( 340) = W( 81 )
  JVS( 341) = W( 83 )
  JVS( 342) = W( 84 )
  JVS( 343) = W( 86 )
  JVS( 344) = W( 87 )
  IF ( ABS( JVS( 346 )) < TINY(a) ) THEN
         IER = 69
         RETURN
  END IF
   W( 67 ) = JVS( 345 )
   W( 69 ) = JVS( 346 )
   W( 77 ) = JVS( 347 )
   W( 81 ) = JVS( 348 )
   W( 83 ) = JVS( 349 )
   W( 84 ) = JVS( 350 )
  a = -W( 67 ) / JVS( 331 )
  W( 67 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 332 )
  W( 83 ) = W( 83 ) + a*JVS( 333 )
  W( 84 ) = W( 84 ) + a*JVS( 334 )
  JVS( 345) = W( 67 )
  JVS( 346) = W( 69 )
  JVS( 347) = W( 77 )
  JVS( 348) = W( 81 )
  JVS( 349) = W( 83 )
  JVS( 350) = W( 84 )
  IF ( ABS( JVS( 363 )) < TINY(a) ) THEN
         IER = 70
         RETURN
  END IF
   W( 25 ) = JVS( 351 )
   W( 28 ) = JVS( 352 )
   W( 31 ) = JVS( 353 )
   W( 32 ) = JVS( 354 )
   W( 46 ) = JVS( 355 )
   W( 56 ) = JVS( 356 )
   W( 57 ) = JVS( 357 )
   W( 59 ) = JVS( 358 )
   W( 63 ) = JVS( 359 )
   W( 64 ) = JVS( 360 )
   W( 66 ) = JVS( 361 )
   W( 69 ) = JVS( 362 )
   W( 70 ) = JVS( 363 )
   W( 77 ) = JVS( 364 )
   W( 81 ) = JVS( 365 )
   W( 82 ) = JVS( 366 )
   W( 83 ) = JVS( 367 )
   W( 84 ) = JVS( 368 )
   W( 86 ) = JVS( 369 )
   W( 87 ) = JVS( 370 )
  a = -W( 25 ) / JVS( 67 )
  W( 25 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 68 )
  W( 82 ) = W( 82 ) + a*JVS( 69 )
  W( 84 ) = W( 84 ) + a*JVS( 70 )
  a = -W( 28 ) / JVS( 80 )
  W( 28 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 81 )
  W( 84 ) = W( 84 ) + a*JVS( 82 )
  W( 86 ) = W( 86 ) + a*JVS( 83 )
  a = -W( 31 ) / JVS( 93 )
  W( 31 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 94 )
  W( 87 ) = W( 87 ) + a*JVS( 95 )
  a = -W( 32 ) / JVS( 96 )
  W( 32 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 97 )
  W( 84 ) = W( 84 ) + a*JVS( 98 )
  W( 86 ) = W( 86 ) + a*JVS( 99 )
  a = -W( 46 ) / JVS( 175 )
  W( 46 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 176 )
  W( 84 ) = W( 84 ) + a*JVS( 177 )
  W( 86 ) = W( 86 ) + a*JVS( 178 )
  a = -W( 56 ) / JVS( 230 )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 231 )
  W( 84 ) = W( 84 ) + a*JVS( 232 )
  W( 86 ) = W( 86 ) + a*JVS( 233 )
  a = -W( 57 ) / JVS( 237 )
  W( 57 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 238 )
  W( 84 ) = W( 84 ) + a*JVS( 239 )
  W( 86 ) = W( 86 ) + a*JVS( 240 )
  W( 87 ) = W( 87 ) + a*JVS( 241 )
  a = -W( 59 ) / JVS( 255 )
  W( 59 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 256 )
  W( 84 ) = W( 84 ) + a*JVS( 257 )
  W( 86 ) = W( 86 ) + a*JVS( 258 )
  W( 87 ) = W( 87 ) + a*JVS( 259 )
  a = -W( 63 ) / JVS( 286 )
  W( 63 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 287 )
  W( 86 ) = W( 86 ) + a*JVS( 288 )
  W( 87 ) = W( 87 ) + a*JVS( 289 )
  a = -W( 64 ) / JVS( 295 )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 296 )
  W( 84 ) = W( 84 ) + a*JVS( 297 )
  W( 86 ) = W( 86 ) + a*JVS( 298 )
  W( 87 ) = W( 87 ) + a*JVS( 299 )
  a = -W( 66 ) / JVS( 326 )
  W( 66 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 327 )
  W( 84 ) = W( 84 ) + a*JVS( 328 )
  W( 86 ) = W( 86 ) + a*JVS( 329 )
  W( 87 ) = W( 87 ) + a*JVS( 330 )
  a = -W( 69 ) / JVS( 346 )
  W( 69 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  W( 83 ) = W( 83 ) + a*JVS( 349 )
  W( 84 ) = W( 84 ) + a*JVS( 350 )
  JVS( 351) = W( 25 )
  JVS( 352) = W( 28 )
  JVS( 353) = W( 31 )
  JVS( 354) = W( 32 )
  JVS( 355) = W( 46 )
  JVS( 356) = W( 56 )
  JVS( 357) = W( 57 )
  JVS( 358) = W( 59 )
  JVS( 359) = W( 63 )
  JVS( 360) = W( 64 )
  JVS( 361) = W( 66 )
  JVS( 362) = W( 69 )
  JVS( 363) = W( 70 )
  JVS( 364) = W( 77 )
  JVS( 365) = W( 81 )
  JVS( 366) = W( 82 )
  JVS( 367) = W( 83 )
  JVS( 368) = W( 84 )
  JVS( 369) = W( 86 )
  JVS( 370) = W( 87 )
  IF ( ABS( JVS( 373 )) < TINY(a) ) THEN
         IER = 71
         RETURN
  END IF
   W( 44 ) = JVS( 371 )
   W( 54 ) = JVS( 372 )
   W( 71 ) = JVS( 373 )
   W( 81 ) = JVS( 374 )
   W( 83 ) = JVS( 375 )
   W( 84 ) = JVS( 376 )
   W( 86 ) = JVS( 377 )
   W( 87 ) = JVS( 378 )
  a = -W( 44 ) / JVS( 162 )
  W( 44 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 163 )
  W( 84 ) = W( 84 ) + a*JVS( 164 )
  W( 86 ) = W( 86 ) + a*JVS( 165 )
  a = -W( 54 ) / JVS( 218 )
  W( 54 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 219 )
  W( 83 ) = W( 83 ) + a*JVS( 220 )
  W( 84 ) = W( 84 ) + a*JVS( 221 )
  JVS( 371) = W( 44 )
  JVS( 372) = W( 54 )
  JVS( 373) = W( 71 )
  JVS( 374) = W( 81 )
  JVS( 375) = W( 83 )
  JVS( 376) = W( 84 )
  JVS( 377) = W( 86 )
  JVS( 378) = W( 87 )
  IF ( ABS( JVS( 384 )) < TINY(a) ) THEN
         IER = 72
         RETURN
  END IF
   W( 30 ) = JVS( 379 )
   W( 58 ) = JVS( 380 )
   W( 63 ) = JVS( 381 )
   W( 66 ) = JVS( 382 )
   W( 71 ) = JVS( 383 )
   W( 72 ) = JVS( 384 )
   W( 81 ) = JVS( 385 )
   W( 82 ) = JVS( 386 )
   W( 83 ) = JVS( 387 )
   W( 84 ) = JVS( 388 )
   W( 86 ) = JVS( 389 )
   W( 87 ) = JVS( 390 )
  a = -W( 30 ) / JVS( 88 )
  W( 30 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 89 )
  W( 84 ) = W( 84 ) + a*JVS( 90 )
  W( 86 ) = W( 86 ) + a*JVS( 91 )
  a = -W( 58 ) / JVS( 246 )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 247 )
  W( 66 ) = W( 66 ) + a*JVS( 248 )
  W( 71 ) = W( 71 ) + a*JVS( 249 )
  W( 82 ) = W( 82 ) + a*JVS( 250 )
  W( 84 ) = W( 84 ) + a*JVS( 251 )
  W( 86 ) = W( 86 ) + a*JVS( 252 )
  W( 87 ) = W( 87 ) + a*JVS( 253 )
  a = -W( 63 ) / JVS( 286 )
  W( 63 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 287 )
  W( 86 ) = W( 86 ) + a*JVS( 288 )
  W( 87 ) = W( 87 ) + a*JVS( 289 )
  a = -W( 66 ) / JVS( 326 )
  W( 66 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 327 )
  W( 84 ) = W( 84 ) + a*JVS( 328 )
  W( 86 ) = W( 86 ) + a*JVS( 329 )
  W( 87 ) = W( 87 ) + a*JVS( 330 )
  a = -W( 71 ) / JVS( 373 )
  W( 71 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 86 ) = W( 86 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  JVS( 379) = W( 30 )
  JVS( 380) = W( 58 )
  JVS( 381) = W( 63 )
  JVS( 382) = W( 66 )
  JVS( 383) = W( 71 )
  JVS( 384) = W( 72 )
  JVS( 385) = W( 81 )
  JVS( 386) = W( 82 )
  JVS( 387) = W( 83 )
  JVS( 388) = W( 84 )
  JVS( 389) = W( 86 )
  JVS( 390) = W( 87 )
  IF ( ABS( JVS( 396 )) < TINY(a) ) THEN
         IER = 73
         RETURN
  END IF
   W( 46 ) = JVS( 391 )
   W( 50 ) = JVS( 392 )
   W( 59 ) = JVS( 393 )
   W( 69 ) = JVS( 394 )
   W( 72 ) = JVS( 395 )
   W( 73 ) = JVS( 396 )
   W( 76 ) = JVS( 397 )
   W( 77 ) = JVS( 398 )
   W( 78 ) = JVS( 399 )
   W( 81 ) = JVS( 400 )
   W( 82 ) = JVS( 401 )
   W( 83 ) = JVS( 402 )
   W( 84 ) = JVS( 403 )
   W( 85 ) = JVS( 404 )
   W( 86 ) = JVS( 405 )
   W( 87 ) = JVS( 406 )
   W( 88 ) = JVS( 407 )
   W( 89 ) = JVS( 408 )
  a = -W( 46 ) / JVS( 175 )
  W( 46 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 176 )
  W( 84 ) = W( 84 ) + a*JVS( 177 )
  W( 86 ) = W( 86 ) + a*JVS( 178 )
  a = -W( 50 ) / JVS( 192 )
  W( 50 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 193 )
  W( 85 ) = W( 85 ) + a*JVS( 194 )
  W( 88 ) = W( 88 ) + a*JVS( 195 )
  a = -W( 59 ) / JVS( 255 )
  W( 59 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 256 )
  W( 84 ) = W( 84 ) + a*JVS( 257 )
  W( 86 ) = W( 86 ) + a*JVS( 258 )
  W( 87 ) = W( 87 ) + a*JVS( 259 )
  a = -W( 69 ) / JVS( 346 )
  W( 69 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  W( 83 ) = W( 83 ) + a*JVS( 349 )
  W( 84 ) = W( 84 ) + a*JVS( 350 )
  a = -W( 72 ) / JVS( 384 )
  W( 72 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 385 )
  W( 82 ) = W( 82 ) + a*JVS( 386 )
  W( 83 ) = W( 83 ) + a*JVS( 387 )
  W( 84 ) = W( 84 ) + a*JVS( 388 )
  W( 86 ) = W( 86 ) + a*JVS( 389 )
  W( 87 ) = W( 87 ) + a*JVS( 390 )
  JVS( 391) = W( 46 )
  JVS( 392) = W( 50 )
  JVS( 393) = W( 59 )
  JVS( 394) = W( 69 )
  JVS( 395) = W( 72 )
  JVS( 396) = W( 73 )
  JVS( 397) = W( 76 )
  JVS( 398) = W( 77 )
  JVS( 399) = W( 78 )
  JVS( 400) = W( 81 )
  JVS( 401) = W( 82 )
  JVS( 402) = W( 83 )
  JVS( 403) = W( 84 )
  JVS( 404) = W( 85 )
  JVS( 405) = W( 86 )
  JVS( 406) = W( 87 )
  JVS( 407) = W( 88 )
  JVS( 408) = W( 89 )
  IF ( ABS( JVS( 415 )) < TINY(a) ) THEN
         IER = 74
         RETURN
  END IF
   W( 44 ) = JVS( 409 )
   W( 49 ) = JVS( 410 )
   W( 54 ) = JVS( 411 )
   W( 61 ) = JVS( 412 )
   W( 67 ) = JVS( 413 )
   W( 71 ) = JVS( 414 )
   W( 74 ) = JVS( 415 )
   W( 80 ) = JVS( 416 )
   W( 81 ) = JVS( 417 )
   W( 82 ) = JVS( 418 )
   W( 83 ) = JVS( 419 )
   W( 84 ) = JVS( 420 )
   W( 86 ) = JVS( 421 )
   W( 87 ) = JVS( 422 )
   W( 89 ) = JVS( 423 )
  a = -W( 44 ) / JVS( 162 )
  W( 44 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 163 )
  W( 84 ) = W( 84 ) + a*JVS( 164 )
  W( 86 ) = W( 86 ) + a*JVS( 165 )
  a = -W( 49 ) / JVS( 188 )
  W( 49 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 189 )
  W( 84 ) = W( 84 ) + a*JVS( 190 )
  W( 86 ) = W( 86 ) + a*JVS( 191 )
  a = -W( 54 ) / JVS( 218 )
  W( 54 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 219 )
  W( 83 ) = W( 83 ) + a*JVS( 220 )
  W( 84 ) = W( 84 ) + a*JVS( 221 )
  a = -W( 61 ) / JVS( 268 )
  W( 61 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 269 )
  W( 81 ) = W( 81 ) + a*JVS( 270 )
  W( 86 ) = W( 86 ) + a*JVS( 271 )
  W( 87 ) = W( 87 ) + a*JVS( 272 )
  a = -W( 67 ) / JVS( 331 )
  W( 67 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 332 )
  W( 83 ) = W( 83 ) + a*JVS( 333 )
  W( 84 ) = W( 84 ) + a*JVS( 334 )
  a = -W( 71 ) / JVS( 373 )
  W( 71 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 86 ) = W( 86 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  JVS( 409) = W( 44 )
  JVS( 410) = W( 49 )
  JVS( 411) = W( 54 )
  JVS( 412) = W( 61 )
  JVS( 413) = W( 67 )
  JVS( 414) = W( 71 )
  JVS( 415) = W( 74 )
  JVS( 416) = W( 80 )
  JVS( 417) = W( 81 )
  JVS( 418) = W( 82 )
  JVS( 419) = W( 83 )
  JVS( 420) = W( 84 )
  JVS( 421) = W( 86 )
  JVS( 422) = W( 87 )
  JVS( 423) = W( 89 )
  IF ( ABS( JVS( 433 )) < TINY(a) ) THEN
         IER = 75
         RETURN
  END IF
   W( 37 ) = JVS( 424 )
   W( 38 ) = JVS( 425 )
   W( 52 ) = JVS( 426 )
   W( 53 ) = JVS( 427 )
   W( 63 ) = JVS( 428 )
   W( 69 ) = JVS( 429 )
   W( 72 ) = JVS( 430 )
   W( 73 ) = JVS( 431 )
   W( 74 ) = JVS( 432 )
   W( 75 ) = JVS( 433 )
   W( 76 ) = JVS( 434 )
   W( 77 ) = JVS( 435 )
   W( 78 ) = JVS( 436 )
   W( 80 ) = JVS( 437 )
   W( 81 ) = JVS( 438 )
   W( 82 ) = JVS( 439 )
   W( 83 ) = JVS( 440 )
   W( 84 ) = JVS( 441 )
   W( 85 ) = JVS( 442 )
   W( 86 ) = JVS( 443 )
   W( 87 ) = JVS( 444 )
   W( 88 ) = JVS( 445 )
   W( 89 ) = JVS( 446 )
  a = -W( 37 ) / JVS( 119 )
  W( 37 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 120 )
  W( 69 ) = W( 69 ) + a*JVS( 121 )
  W( 81 ) = W( 81 ) + a*JVS( 122 )
  W( 84 ) = W( 84 ) + a*JVS( 123 )
  W( 87 ) = W( 87 ) + a*JVS( 124 )
  a = -W( 38 ) / JVS( 125 )
  W( 38 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 126 )
  W( 84 ) = W( 84 ) + a*JVS( 127 )
  W( 86 ) = W( 86 ) + a*JVS( 128 )
  a = -W( 52 ) / JVS( 206 )
  W( 52 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 207 )
  W( 86 ) = W( 86 ) + a*JVS( 208 )
  W( 87 ) = W( 87 ) + a*JVS( 209 )
  a = -W( 53 ) / JVS( 213 )
  W( 53 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 214 )
  W( 86 ) = W( 86 ) + a*JVS( 215 )
  W( 87 ) = W( 87 ) + a*JVS( 216 )
  W( 88 ) = W( 88 ) + a*JVS( 217 )
  a = -W( 63 ) / JVS( 286 )
  W( 63 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 287 )
  W( 86 ) = W( 86 ) + a*JVS( 288 )
  W( 87 ) = W( 87 ) + a*JVS( 289 )
  a = -W( 69 ) / JVS( 346 )
  W( 69 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  W( 83 ) = W( 83 ) + a*JVS( 349 )
  W( 84 ) = W( 84 ) + a*JVS( 350 )
  a = -W( 72 ) / JVS( 384 )
  W( 72 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 385 )
  W( 82 ) = W( 82 ) + a*JVS( 386 )
  W( 83 ) = W( 83 ) + a*JVS( 387 )
  W( 84 ) = W( 84 ) + a*JVS( 388 )
  W( 86 ) = W( 86 ) + a*JVS( 389 )
  W( 87 ) = W( 87 ) + a*JVS( 390 )
  a = -W( 73 ) / JVS( 396 )
  W( 73 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 397 )
  W( 77 ) = W( 77 ) + a*JVS( 398 )
  W( 78 ) = W( 78 ) + a*JVS( 399 )
  W( 81 ) = W( 81 ) + a*JVS( 400 )
  W( 82 ) = W( 82 ) + a*JVS( 401 )
  W( 83 ) = W( 83 ) + a*JVS( 402 )
  W( 84 ) = W( 84 ) + a*JVS( 403 )
  W( 85 ) = W( 85 ) + a*JVS( 404 )
  W( 86 ) = W( 86 ) + a*JVS( 405 )
  W( 87 ) = W( 87 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 89 ) = W( 89 ) + a*JVS( 408 )
  a = -W( 74 ) / JVS( 415 )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 416 )
  W( 81 ) = W( 81 ) + a*JVS( 417 )
  W( 82 ) = W( 82 ) + a*JVS( 418 )
  W( 83 ) = W( 83 ) + a*JVS( 419 )
  W( 84 ) = W( 84 ) + a*JVS( 420 )
  W( 86 ) = W( 86 ) + a*JVS( 421 )
  W( 87 ) = W( 87 ) + a*JVS( 422 )
  W( 89 ) = W( 89 ) + a*JVS( 423 )
  JVS( 424) = W( 37 )
  JVS( 425) = W( 38 )
  JVS( 426) = W( 52 )
  JVS( 427) = W( 53 )
  JVS( 428) = W( 63 )
  JVS( 429) = W( 69 )
  JVS( 430) = W( 72 )
  JVS( 431) = W( 73 )
  JVS( 432) = W( 74 )
  JVS( 433) = W( 75 )
  JVS( 434) = W( 76 )
  JVS( 435) = W( 77 )
  JVS( 436) = W( 78 )
  JVS( 437) = W( 80 )
  JVS( 438) = W( 81 )
  JVS( 439) = W( 82 )
  JVS( 440) = W( 83 )
  JVS( 441) = W( 84 )
  JVS( 442) = W( 85 )
  JVS( 443) = W( 86 )
  JVS( 444) = W( 87 )
  JVS( 445) = W( 88 )
  JVS( 446) = W( 89 )
  IF ( ABS( JVS( 451 )) < TINY(a) ) THEN
         IER = 76
         RETURN
  END IF
   W( 23 ) = JVS( 447 )
   W( 35 ) = JVS( 448 )
   W( 49 ) = JVS( 449 )
   W( 68 ) = JVS( 450 )
   W( 76 ) = JVS( 451 )
   W( 78 ) = JVS( 452 )
   W( 80 ) = JVS( 453 )
   W( 81 ) = JVS( 454 )
   W( 82 ) = JVS( 455 )
   W( 83 ) = JVS( 456 )
   W( 84 ) = JVS( 457 )
   W( 86 ) = JVS( 458 )
   W( 87 ) = JVS( 459 )
   W( 89 ) = JVS( 460 )
  a = -W( 23 ) / JVS( 60 )
  W( 23 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 61 )
  W( 84 ) = W( 84 ) + a*JVS( 62 )
  W( 86 ) = W( 86 ) + a*JVS( 63 )
  a = -W( 35 ) / JVS( 108 )
  W( 35 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 109 )
  W( 80 ) = W( 80 ) + a*JVS( 110 )
  W( 81 ) = W( 81 ) + a*JVS( 111 )
  W( 82 ) = W( 82 ) + a*JVS( 112 )
  W( 84 ) = W( 84 ) + a*JVS( 113 )
  W( 87 ) = W( 87 ) + a*JVS( 114 )
  W( 89 ) = W( 89 ) + a*JVS( 115 )
  a = -W( 49 ) / JVS( 188 )
  W( 49 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 189 )
  W( 84 ) = W( 84 ) + a*JVS( 190 )
  W( 86 ) = W( 86 ) + a*JVS( 191 )
  a = -W( 68 ) / JVS( 337 )
  W( 68 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 338 )
  W( 80 ) = W( 80 ) + a*JVS( 339 )
  W( 81 ) = W( 81 ) + a*JVS( 340 )
  W( 83 ) = W( 83 ) + a*JVS( 341 )
  W( 84 ) = W( 84 ) + a*JVS( 342 )
  W( 86 ) = W( 86 ) + a*JVS( 343 )
  W( 87 ) = W( 87 ) + a*JVS( 344 )
  JVS( 447) = W( 23 )
  JVS( 448) = W( 35 )
  JVS( 449) = W( 49 )
  JVS( 450) = W( 68 )
  JVS( 451) = W( 76 )
  JVS( 452) = W( 78 )
  JVS( 453) = W( 80 )
  JVS( 454) = W( 81 )
  JVS( 455) = W( 82 )
  JVS( 456) = W( 83 )
  JVS( 457) = W( 84 )
  JVS( 458) = W( 86 )
  JVS( 459) = W( 87 )
  JVS( 460) = W( 89 )
  IF ( ABS( JVS( 467 )) < TINY(a) ) THEN
         IER = 77
         RETURN
  END IF
   W( 44 ) = JVS( 461 )
   W( 49 ) = JVS( 462 )
   W( 54 ) = JVS( 463 )
   W( 61 ) = JVS( 464 )
   W( 67 ) = JVS( 465 )
   W( 71 ) = JVS( 466 )
   W( 77 ) = JVS( 467 )
   W( 80 ) = JVS( 468 )
   W( 81 ) = JVS( 469 )
   W( 82 ) = JVS( 470 )
   W( 83 ) = JVS( 471 )
   W( 84 ) = JVS( 472 )
   W( 86 ) = JVS( 473 )
   W( 87 ) = JVS( 474 )
   W( 89 ) = JVS( 475 )
  a = -W( 44 ) / JVS( 162 )
  W( 44 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 163 )
  W( 84 ) = W( 84 ) + a*JVS( 164 )
  W( 86 ) = W( 86 ) + a*JVS( 165 )
  a = -W( 49 ) / JVS( 188 )
  W( 49 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 189 )
  W( 84 ) = W( 84 ) + a*JVS( 190 )
  W( 86 ) = W( 86 ) + a*JVS( 191 )
  a = -W( 54 ) / JVS( 218 )
  W( 54 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 219 )
  W( 83 ) = W( 83 ) + a*JVS( 220 )
  W( 84 ) = W( 84 ) + a*JVS( 221 )
  a = -W( 61 ) / JVS( 268 )
  W( 61 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 269 )
  W( 81 ) = W( 81 ) + a*JVS( 270 )
  W( 86 ) = W( 86 ) + a*JVS( 271 )
  W( 87 ) = W( 87 ) + a*JVS( 272 )
  a = -W( 67 ) / JVS( 331 )
  W( 67 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 332 )
  W( 83 ) = W( 83 ) + a*JVS( 333 )
  W( 84 ) = W( 84 ) + a*JVS( 334 )
  a = -W( 71 ) / JVS( 373 )
  W( 71 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 86 ) = W( 86 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  JVS( 461) = W( 44 )
  JVS( 462) = W( 49 )
  JVS( 463) = W( 54 )
  JVS( 464) = W( 61 )
  JVS( 465) = W( 67 )
  JVS( 466) = W( 71 )
  JVS( 467) = W( 77 )
  JVS( 468) = W( 80 )
  JVS( 469) = W( 81 )
  JVS( 470) = W( 82 )
  JVS( 471) = W( 83 )
  JVS( 472) = W( 84 )
  JVS( 473) = W( 86 )
  JVS( 474) = W( 87 )
  JVS( 475) = W( 89 )
  IF ( ABS( JVS( 479 )) < TINY(a) ) THEN
         IER = 78
         RETURN
  END IF
   W( 26 ) = JVS( 476 )
   W( 74 ) = JVS( 477 )
   W( 77 ) = JVS( 478 )
   W( 78 ) = JVS( 479 )
   W( 80 ) = JVS( 480 )
   W( 81 ) = JVS( 481 )
   W( 82 ) = JVS( 482 )
   W( 83 ) = JVS( 483 )
   W( 84 ) = JVS( 484 )
   W( 86 ) = JVS( 485 )
   W( 87 ) = JVS( 486 )
   W( 89 ) = JVS( 487 )
  a = -W( 26 ) / JVS( 71 )
  W( 26 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 72 )
  W( 84 ) = W( 84 ) + a*JVS( 73 )
  W( 86 ) = W( 86 ) + a*JVS( 74 )
  a = -W( 74 ) / JVS( 415 )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 416 )
  W( 81 ) = W( 81 ) + a*JVS( 417 )
  W( 82 ) = W( 82 ) + a*JVS( 418 )
  W( 83 ) = W( 83 ) + a*JVS( 419 )
  W( 84 ) = W( 84 ) + a*JVS( 420 )
  W( 86 ) = W( 86 ) + a*JVS( 421 )
  W( 87 ) = W( 87 ) + a*JVS( 422 )
  W( 89 ) = W( 89 ) + a*JVS( 423 )
  a = -W( 77 ) / JVS( 467 )
  W( 77 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 468 )
  W( 81 ) = W( 81 ) + a*JVS( 469 )
  W( 82 ) = W( 82 ) + a*JVS( 470 )
  W( 83 ) = W( 83 ) + a*JVS( 471 )
  W( 84 ) = W( 84 ) + a*JVS( 472 )
  W( 86 ) = W( 86 ) + a*JVS( 473 )
  W( 87 ) = W( 87 ) + a*JVS( 474 )
  W( 89 ) = W( 89 ) + a*JVS( 475 )
  JVS( 476) = W( 26 )
  JVS( 477) = W( 74 )
  JVS( 478) = W( 77 )
  JVS( 479) = W( 78 )
  JVS( 480) = W( 80 )
  JVS( 481) = W( 81 )
  JVS( 482) = W( 82 )
  JVS( 483) = W( 83 )
  JVS( 484) = W( 84 )
  JVS( 485) = W( 86 )
  JVS( 486) = W( 87 )
  JVS( 487) = W( 89 )
  IF ( ABS( JVS( 518 )) < TINY(a) ) THEN
         IER = 79
         RETURN
  END IF
   W( 20 ) = JVS( 488 )
   W( 30 ) = JVS( 489 )
   W( 31 ) = JVS( 490 )
   W( 33 ) = JVS( 491 )
   W( 36 ) = JVS( 492 )
   W( 39 ) = JVS( 493 )
   W( 40 ) = JVS( 494 )
   W( 46 ) = JVS( 495 )
   W( 47 ) = JVS( 496 )
   W( 48 ) = JVS( 497 )
   W( 49 ) = JVS( 498 )
   W( 50 ) = JVS( 499 )
   W( 51 ) = JVS( 500 )
   W( 55 ) = JVS( 501 )
   W( 56 ) = JVS( 502 )
   W( 59 ) = JVS( 503 )
   W( 61 ) = JVS( 504 )
   W( 62 ) = JVS( 505 )
   W( 63 ) = JVS( 506 )
   W( 64 ) = JVS( 507 )
   W( 66 ) = JVS( 508 )
   W( 67 ) = JVS( 509 )
   W( 68 ) = JVS( 510 )
   W( 69 ) = JVS( 511 )
   W( 72 ) = JVS( 512 )
   W( 73 ) = JVS( 513 )
   W( 74 ) = JVS( 514 )
   W( 76 ) = JVS( 515 )
   W( 77 ) = JVS( 516 )
   W( 78 ) = JVS( 517 )
   W( 79 ) = JVS( 518 )
   W( 80 ) = JVS( 519 )
   W( 81 ) = JVS( 520 )
   W( 82 ) = JVS( 521 )
   W( 83 ) = JVS( 522 )
   W( 84 ) = JVS( 523 )
   W( 85 ) = JVS( 524 )
   W( 86 ) = JVS( 525 )
   W( 87 ) = JVS( 526 )
   W( 88 ) = JVS( 527 )
   W( 89 ) = JVS( 528 )
  a = -W( 20 ) / JVS( 50 )
  W( 20 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 51 )
  W( 87 ) = W( 87 ) + a*JVS( 52 )
  a = -W( 30 ) / JVS( 88 )
  W( 30 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 89 )
  W( 84 ) = W( 84 ) + a*JVS( 90 )
  W( 86 ) = W( 86 ) + a*JVS( 91 )
  a = -W( 31 ) / JVS( 93 )
  W( 31 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 94 )
  W( 87 ) = W( 87 ) + a*JVS( 95 )
  a = -W( 33 ) / JVS( 100 )
  W( 33 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 101 )
  W( 84 ) = W( 84 ) + a*JVS( 102 )
  W( 86 ) = W( 86 ) + a*JVS( 103 )
  a = -W( 36 ) / JVS( 116 )
  W( 36 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 117 )
  W( 84 ) = W( 84 ) + a*JVS( 118 )
  a = -W( 39 ) / JVS( 130 )
  W( 39 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 131 )
  W( 84 ) = W( 84 ) + a*JVS( 132 )
  W( 87 ) = W( 87 ) + a*JVS( 133 )
  a = -W( 40 ) / JVS( 134 )
  W( 40 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 135 )
  W( 85 ) = W( 85 ) + a*JVS( 136 )
  W( 86 ) = W( 86 ) + a*JVS( 137 )
  W( 89 ) = W( 89 ) + a*JVS( 138 )
  a = -W( 46 ) / JVS( 175 )
  W( 46 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 176 )
  W( 84 ) = W( 84 ) + a*JVS( 177 )
  W( 86 ) = W( 86 ) + a*JVS( 178 )
  a = -W( 47 ) / JVS( 179 )
  W( 47 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 180 )
  W( 88 ) = W( 88 ) + a*JVS( 181 )
  W( 89 ) = W( 89 ) + a*JVS( 182 )
  a = -W( 48 ) / JVS( 183 )
  W( 48 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 184 )
  W( 69 ) = W( 69 ) + a*JVS( 185 )
  W( 83 ) = W( 83 ) + a*JVS( 186 )
  W( 84 ) = W( 84 ) + a*JVS( 187 )
  a = -W( 49 ) / JVS( 188 )
  W( 49 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 189 )
  W( 84 ) = W( 84 ) + a*JVS( 190 )
  W( 86 ) = W( 86 ) + a*JVS( 191 )
  a = -W( 50 ) / JVS( 192 )
  W( 50 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 193 )
  W( 85 ) = W( 85 ) + a*JVS( 194 )
  W( 88 ) = W( 88 ) + a*JVS( 195 )
  a = -W( 51 ) / JVS( 196 )
  W( 51 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 197 )
  W( 72 ) = W( 72 ) + a*JVS( 198 )
  W( 76 ) = W( 76 ) + a*JVS( 199 )
  W( 78 ) = W( 78 ) + a*JVS( 200 )
  W( 80 ) = W( 80 ) + a*JVS( 201 )
  W( 82 ) = W( 82 ) + a*JVS( 202 )
  W( 84 ) = W( 84 ) + a*JVS( 203 )
  a = -W( 55 ) / JVS( 225 )
  W( 55 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 226 )
  W( 79 ) = W( 79 ) + a*JVS( 227 )
  W( 83 ) = W( 83 ) + a*JVS( 228 )
  W( 84 ) = W( 84 ) + a*JVS( 229 )
  a = -W( 56 ) / JVS( 230 )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 231 )
  W( 84 ) = W( 84 ) + a*JVS( 232 )
  W( 86 ) = W( 86 ) + a*JVS( 233 )
  a = -W( 59 ) / JVS( 255 )
  W( 59 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 256 )
  W( 84 ) = W( 84 ) + a*JVS( 257 )
  W( 86 ) = W( 86 ) + a*JVS( 258 )
  W( 87 ) = W( 87 ) + a*JVS( 259 )
  a = -W( 61 ) / JVS( 268 )
  W( 61 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 269 )
  W( 81 ) = W( 81 ) + a*JVS( 270 )
  W( 86 ) = W( 86 ) + a*JVS( 271 )
  W( 87 ) = W( 87 ) + a*JVS( 272 )
  a = -W( 62 ) / JVS( 275 )
  W( 62 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 276 )
  W( 78 ) = W( 78 ) + a*JVS( 277 )
  W( 81 ) = W( 81 ) + a*JVS( 278 )
  W( 82 ) = W( 82 ) + a*JVS( 279 )
  W( 83 ) = W( 83 ) + a*JVS( 280 )
  W( 84 ) = W( 84 ) + a*JVS( 281 )
  W( 87 ) = W( 87 ) + a*JVS( 282 )
  W( 89 ) = W( 89 ) + a*JVS( 283 )
  a = -W( 63 ) / JVS( 286 )
  W( 63 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 287 )
  W( 86 ) = W( 86 ) + a*JVS( 288 )
  W( 87 ) = W( 87 ) + a*JVS( 289 )
  a = -W( 64 ) / JVS( 295 )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 296 )
  W( 84 ) = W( 84 ) + a*JVS( 297 )
  W( 86 ) = W( 86 ) + a*JVS( 298 )
  W( 87 ) = W( 87 ) + a*JVS( 299 )
  a = -W( 66 ) / JVS( 326 )
  W( 66 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 327 )
  W( 84 ) = W( 84 ) + a*JVS( 328 )
  W( 86 ) = W( 86 ) + a*JVS( 329 )
  W( 87 ) = W( 87 ) + a*JVS( 330 )
  a = -W( 67 ) / JVS( 331 )
  W( 67 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 332 )
  W( 83 ) = W( 83 ) + a*JVS( 333 )
  W( 84 ) = W( 84 ) + a*JVS( 334 )
  a = -W( 68 ) / JVS( 337 )
  W( 68 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 338 )
  W( 80 ) = W( 80 ) + a*JVS( 339 )
  W( 81 ) = W( 81 ) + a*JVS( 340 )
  W( 83 ) = W( 83 ) + a*JVS( 341 )
  W( 84 ) = W( 84 ) + a*JVS( 342 )
  W( 86 ) = W( 86 ) + a*JVS( 343 )
  W( 87 ) = W( 87 ) + a*JVS( 344 )
  a = -W( 69 ) / JVS( 346 )
  W( 69 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  W( 83 ) = W( 83 ) + a*JVS( 349 )
  W( 84 ) = W( 84 ) + a*JVS( 350 )
  a = -W( 72 ) / JVS( 384 )
  W( 72 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 385 )
  W( 82 ) = W( 82 ) + a*JVS( 386 )
  W( 83 ) = W( 83 ) + a*JVS( 387 )
  W( 84 ) = W( 84 ) + a*JVS( 388 )
  W( 86 ) = W( 86 ) + a*JVS( 389 )
  W( 87 ) = W( 87 ) + a*JVS( 390 )
  a = -W( 73 ) / JVS( 396 )
  W( 73 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 397 )
  W( 77 ) = W( 77 ) + a*JVS( 398 )
  W( 78 ) = W( 78 ) + a*JVS( 399 )
  W( 81 ) = W( 81 ) + a*JVS( 400 )
  W( 82 ) = W( 82 ) + a*JVS( 401 )
  W( 83 ) = W( 83 ) + a*JVS( 402 )
  W( 84 ) = W( 84 ) + a*JVS( 403 )
  W( 85 ) = W( 85 ) + a*JVS( 404 )
  W( 86 ) = W( 86 ) + a*JVS( 405 )
  W( 87 ) = W( 87 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 89 ) = W( 89 ) + a*JVS( 408 )
  a = -W( 74 ) / JVS( 415 )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 416 )
  W( 81 ) = W( 81 ) + a*JVS( 417 )
  W( 82 ) = W( 82 ) + a*JVS( 418 )
  W( 83 ) = W( 83 ) + a*JVS( 419 )
  W( 84 ) = W( 84 ) + a*JVS( 420 )
  W( 86 ) = W( 86 ) + a*JVS( 421 )
  W( 87 ) = W( 87 ) + a*JVS( 422 )
  W( 89 ) = W( 89 ) + a*JVS( 423 )
  a = -W( 76 ) / JVS( 451 )
  W( 76 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 452 )
  W( 80 ) = W( 80 ) + a*JVS( 453 )
  W( 81 ) = W( 81 ) + a*JVS( 454 )
  W( 82 ) = W( 82 ) + a*JVS( 455 )
  W( 83 ) = W( 83 ) + a*JVS( 456 )
  W( 84 ) = W( 84 ) + a*JVS( 457 )
  W( 86 ) = W( 86 ) + a*JVS( 458 )
  W( 87 ) = W( 87 ) + a*JVS( 459 )
  W( 89 ) = W( 89 ) + a*JVS( 460 )
  a = -W( 77 ) / JVS( 467 )
  W( 77 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 468 )
  W( 81 ) = W( 81 ) + a*JVS( 469 )
  W( 82 ) = W( 82 ) + a*JVS( 470 )
  W( 83 ) = W( 83 ) + a*JVS( 471 )
  W( 84 ) = W( 84 ) + a*JVS( 472 )
  W( 86 ) = W( 86 ) + a*JVS( 473 )
  W( 87 ) = W( 87 ) + a*JVS( 474 )
  W( 89 ) = W( 89 ) + a*JVS( 475 )
  a = -W( 78 ) / JVS( 479 )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 480 )
  W( 81 ) = W( 81 ) + a*JVS( 481 )
  W( 82 ) = W( 82 ) + a*JVS( 482 )
  W( 83 ) = W( 83 ) + a*JVS( 483 )
  W( 84 ) = W( 84 ) + a*JVS( 484 )
  W( 86 ) = W( 86 ) + a*JVS( 485 )
  W( 87 ) = W( 87 ) + a*JVS( 486 )
  W( 89 ) = W( 89 ) + a*JVS( 487 )
  JVS( 488) = W( 20 )
  JVS( 489) = W( 30 )
  JVS( 490) = W( 31 )
  JVS( 491) = W( 33 )
  JVS( 492) = W( 36 )
  JVS( 493) = W( 39 )
  JVS( 494) = W( 40 )
  JVS( 495) = W( 46 )
  JVS( 496) = W( 47 )
  JVS( 497) = W( 48 )
  JVS( 498) = W( 49 )
  JVS( 499) = W( 50 )
  JVS( 500) = W( 51 )
  JVS( 501) = W( 55 )
  JVS( 502) = W( 56 )
  JVS( 503) = W( 59 )
  JVS( 504) = W( 61 )
  JVS( 505) = W( 62 )
  JVS( 506) = W( 63 )
  JVS( 507) = W( 64 )
  JVS( 508) = W( 66 )
  JVS( 509) = W( 67 )
  JVS( 510) = W( 68 )
  JVS( 511) = W( 69 )
  JVS( 512) = W( 72 )
  JVS( 513) = W( 73 )
  JVS( 514) = W( 74 )
  JVS( 515) = W( 76 )
  JVS( 516) = W( 77 )
  JVS( 517) = W( 78 )
  JVS( 518) = W( 79 )
  JVS( 519) = W( 80 )
  JVS( 520) = W( 81 )
  JVS( 521) = W( 82 )
  JVS( 522) = W( 83 )
  JVS( 523) = W( 84 )
  JVS( 524) = W( 85 )
  JVS( 525) = W( 86 )
  JVS( 526) = W( 87 )
  JVS( 527) = W( 88 )
  JVS( 528) = W( 89 )
  IF ( ABS( JVS( 531 )) < TINY(a) ) THEN
         IER = 80
         RETURN
  END IF
   W( 49 ) = JVS( 529 )
   W( 67 ) = JVS( 530 )
   W( 80 ) = JVS( 531 )
   W( 81 ) = JVS( 532 )
   W( 82 ) = JVS( 533 )
   W( 83 ) = JVS( 534 )
   W( 84 ) = JVS( 535 )
   W( 86 ) = JVS( 536 )
   W( 87 ) = JVS( 537 )
   W( 89 ) = JVS( 538 )
  a = -W( 49 ) / JVS( 188 )
  W( 49 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 189 )
  W( 84 ) = W( 84 ) + a*JVS( 190 )
  W( 86 ) = W( 86 ) + a*JVS( 191 )
  a = -W( 67 ) / JVS( 331 )
  W( 67 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 332 )
  W( 83 ) = W( 83 ) + a*JVS( 333 )
  W( 84 ) = W( 84 ) + a*JVS( 334 )
  JVS( 529) = W( 49 )
  JVS( 530) = W( 67 )
  JVS( 531) = W( 80 )
  JVS( 532) = W( 81 )
  JVS( 533) = W( 82 )
  JVS( 534) = W( 83 )
  JVS( 535) = W( 84 )
  JVS( 536) = W( 86 )
  JVS( 537) = W( 87 )
  JVS( 538) = W( 89 )
  IF ( ABS( JVS( 557 )) < TINY(a) ) THEN
         IER = 81
         RETURN
  END IF
   W( 21 ) = JVS( 539 )
   W( 24 ) = JVS( 540 )
   W( 34 ) = JVS( 541 )
   W( 45 ) = JVS( 542 )
   W( 47 ) = JVS( 543 )
   W( 50 ) = JVS( 544 )
   W( 54 ) = JVS( 545 )
   W( 61 ) = JVS( 546 )
   W( 67 ) = JVS( 547 )
   W( 68 ) = JVS( 548 )
   W( 69 ) = JVS( 549 )
   W( 70 ) = JVS( 550 )
   W( 75 ) = JVS( 551 )
   W( 76 ) = JVS( 552 )
   W( 77 ) = JVS( 553 )
   W( 78 ) = JVS( 554 )
   W( 79 ) = JVS( 555 )
   W( 80 ) = JVS( 556 )
   W( 81 ) = JVS( 557 )
   W( 82 ) = JVS( 558 )
   W( 83 ) = JVS( 559 )
   W( 84 ) = JVS( 560 )
   W( 85 ) = JVS( 561 )
   W( 86 ) = JVS( 562 )
   W( 87 ) = JVS( 563 )
   W( 88 ) = JVS( 564 )
   W( 89 ) = JVS( 565 )
  a = -W( 21 ) / JVS( 53 )
  W( 21 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 54 )
  W( 88 ) = W( 88 ) + a*JVS( 55 )
  a = -W( 24 ) / JVS( 64 )
  W( 24 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 65 )
  W( 84 ) = W( 84 ) + a*JVS( 66 )
  a = -W( 34 ) / JVS( 104 )
  W( 34 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 105 )
  W( 86 ) = W( 86 ) + a*JVS( 106 )
  W( 88 ) = W( 88 ) + a*JVS( 107 )
  a = -W( 45 ) / JVS( 168 )
  W( 45 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 169 )
  W( 75 ) = W( 75 ) + a*JVS( 170 )
  W( 79 ) = W( 79 ) + a*JVS( 171 )
  W( 81 ) = W( 81 ) + a*JVS( 172 )
  W( 84 ) = W( 84 ) + a*JVS( 173 )
  W( 88 ) = W( 88 ) + a*JVS( 174 )
  a = -W( 47 ) / JVS( 179 )
  W( 47 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 180 )
  W( 88 ) = W( 88 ) + a*JVS( 181 )
  W( 89 ) = W( 89 ) + a*JVS( 182 )
  a = -W( 50 ) / JVS( 192 )
  W( 50 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 193 )
  W( 85 ) = W( 85 ) + a*JVS( 194 )
  W( 88 ) = W( 88 ) + a*JVS( 195 )
  a = -W( 54 ) / JVS( 218 )
  W( 54 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 219 )
  W( 83 ) = W( 83 ) + a*JVS( 220 )
  W( 84 ) = W( 84 ) + a*JVS( 221 )
  a = -W( 61 ) / JVS( 268 )
  W( 61 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 269 )
  W( 81 ) = W( 81 ) + a*JVS( 270 )
  W( 86 ) = W( 86 ) + a*JVS( 271 )
  W( 87 ) = W( 87 ) + a*JVS( 272 )
  a = -W( 67 ) / JVS( 331 )
  W( 67 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 332 )
  W( 83 ) = W( 83 ) + a*JVS( 333 )
  W( 84 ) = W( 84 ) + a*JVS( 334 )
  a = -W( 68 ) / JVS( 337 )
  W( 68 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 338 )
  W( 80 ) = W( 80 ) + a*JVS( 339 )
  W( 81 ) = W( 81 ) + a*JVS( 340 )
  W( 83 ) = W( 83 ) + a*JVS( 341 )
  W( 84 ) = W( 84 ) + a*JVS( 342 )
  W( 86 ) = W( 86 ) + a*JVS( 343 )
  W( 87 ) = W( 87 ) + a*JVS( 344 )
  a = -W( 69 ) / JVS( 346 )
  W( 69 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  W( 83 ) = W( 83 ) + a*JVS( 349 )
  W( 84 ) = W( 84 ) + a*JVS( 350 )
  a = -W( 70 ) / JVS( 363 )
  W( 70 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 364 )
  W( 81 ) = W( 81 ) + a*JVS( 365 )
  W( 82 ) = W( 82 ) + a*JVS( 366 )
  W( 83 ) = W( 83 ) + a*JVS( 367 )
  W( 84 ) = W( 84 ) + a*JVS( 368 )
  W( 86 ) = W( 86 ) + a*JVS( 369 )
  W( 87 ) = W( 87 ) + a*JVS( 370 )
  a = -W( 75 ) / JVS( 433 )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 434 )
  W( 77 ) = W( 77 ) + a*JVS( 435 )
  W( 78 ) = W( 78 ) + a*JVS( 436 )
  W( 80 ) = W( 80 ) + a*JVS( 437 )
  W( 81 ) = W( 81 ) + a*JVS( 438 )
  W( 82 ) = W( 82 ) + a*JVS( 439 )
  W( 83 ) = W( 83 ) + a*JVS( 440 )
  W( 84 ) = W( 84 ) + a*JVS( 441 )
  W( 85 ) = W( 85 ) + a*JVS( 442 )
  W( 86 ) = W( 86 ) + a*JVS( 443 )
  W( 87 ) = W( 87 ) + a*JVS( 444 )
  W( 88 ) = W( 88 ) + a*JVS( 445 )
  W( 89 ) = W( 89 ) + a*JVS( 446 )
  a = -W( 76 ) / JVS( 451 )
  W( 76 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 452 )
  W( 80 ) = W( 80 ) + a*JVS( 453 )
  W( 81 ) = W( 81 ) + a*JVS( 454 )
  W( 82 ) = W( 82 ) + a*JVS( 455 )
  W( 83 ) = W( 83 ) + a*JVS( 456 )
  W( 84 ) = W( 84 ) + a*JVS( 457 )
  W( 86 ) = W( 86 ) + a*JVS( 458 )
  W( 87 ) = W( 87 ) + a*JVS( 459 )
  W( 89 ) = W( 89 ) + a*JVS( 460 )
  a = -W( 77 ) / JVS( 467 )
  W( 77 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 468 )
  W( 81 ) = W( 81 ) + a*JVS( 469 )
  W( 82 ) = W( 82 ) + a*JVS( 470 )
  W( 83 ) = W( 83 ) + a*JVS( 471 )
  W( 84 ) = W( 84 ) + a*JVS( 472 )
  W( 86 ) = W( 86 ) + a*JVS( 473 )
  W( 87 ) = W( 87 ) + a*JVS( 474 )
  W( 89 ) = W( 89 ) + a*JVS( 475 )
  a = -W( 78 ) / JVS( 479 )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 480 )
  W( 81 ) = W( 81 ) + a*JVS( 481 )
  W( 82 ) = W( 82 ) + a*JVS( 482 )
  W( 83 ) = W( 83 ) + a*JVS( 483 )
  W( 84 ) = W( 84 ) + a*JVS( 484 )
  W( 86 ) = W( 86 ) + a*JVS( 485 )
  W( 87 ) = W( 87 ) + a*JVS( 486 )
  W( 89 ) = W( 89 ) + a*JVS( 487 )
  a = -W( 79 ) / JVS( 518 )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 519 )
  W( 81 ) = W( 81 ) + a*JVS( 520 )
  W( 82 ) = W( 82 ) + a*JVS( 521 )
  W( 83 ) = W( 83 ) + a*JVS( 522 )
  W( 84 ) = W( 84 ) + a*JVS( 523 )
  W( 85 ) = W( 85 ) + a*JVS( 524 )
  W( 86 ) = W( 86 ) + a*JVS( 525 )
  W( 87 ) = W( 87 ) + a*JVS( 526 )
  W( 88 ) = W( 88 ) + a*JVS( 527 )
  W( 89 ) = W( 89 ) + a*JVS( 528 )
  a = -W( 80 ) / JVS( 531 )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 532 )
  W( 82 ) = W( 82 ) + a*JVS( 533 )
  W( 83 ) = W( 83 ) + a*JVS( 534 )
  W( 84 ) = W( 84 ) + a*JVS( 535 )
  W( 86 ) = W( 86 ) + a*JVS( 536 )
  W( 87 ) = W( 87 ) + a*JVS( 537 )
  W( 89 ) = W( 89 ) + a*JVS( 538 )
  JVS( 539) = W( 21 )
  JVS( 540) = W( 24 )
  JVS( 541) = W( 34 )
  JVS( 542) = W( 45 )
  JVS( 543) = W( 47 )
  JVS( 544) = W( 50 )
  JVS( 545) = W( 54 )
  JVS( 546) = W( 61 )
  JVS( 547) = W( 67 )
  JVS( 548) = W( 68 )
  JVS( 549) = W( 69 )
  JVS( 550) = W( 70 )
  JVS( 551) = W( 75 )
  JVS( 552) = W( 76 )
  JVS( 553) = W( 77 )
  JVS( 554) = W( 78 )
  JVS( 555) = W( 79 )
  JVS( 556) = W( 80 )
  JVS( 557) = W( 81 )
  JVS( 558) = W( 82 )
  JVS( 559) = W( 83 )
  JVS( 560) = W( 84 )
  JVS( 561) = W( 85 )
  JVS( 562) = W( 86 )
  JVS( 563) = W( 87 )
  JVS( 564) = W( 88 )
  JVS( 565) = W( 89 )
  IF ( ABS( JVS( 587 )) < TINY(a) ) THEN
         IER = 82
         RETURN
  END IF
   W( 33 ) = JVS( 566 )
   W( 40 ) = JVS( 567 )
   W( 42 ) = JVS( 568 )
   W( 47 ) = JVS( 569 )
   W( 48 ) = JVS( 570 )
   W( 55 ) = JVS( 571 )
   W( 58 ) = JVS( 572 )
   W( 63 ) = JVS( 573 )
   W( 64 ) = JVS( 574 )
   W( 66 ) = JVS( 575 )
   W( 67 ) = JVS( 576 )
   W( 69 ) = JVS( 577 )
   W( 70 ) = JVS( 578 )
   W( 71 ) = JVS( 579 )
   W( 72 ) = JVS( 580 )
   W( 76 ) = JVS( 581 )
   W( 77 ) = JVS( 582 )
   W( 78 ) = JVS( 583 )
   W( 79 ) = JVS( 584 )
   W( 80 ) = JVS( 585 )
   W( 81 ) = JVS( 586 )
   W( 82 ) = JVS( 587 )
   W( 83 ) = JVS( 588 )
   W( 84 ) = JVS( 589 )
   W( 85 ) = JVS( 590 )
   W( 86 ) = JVS( 591 )
   W( 87 ) = JVS( 592 )
   W( 88 ) = JVS( 593 )
   W( 89 ) = JVS( 594 )
  a = -W( 33 ) / JVS( 100 )
  W( 33 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 101 )
  W( 84 ) = W( 84 ) + a*JVS( 102 )
  W( 86 ) = W( 86 ) + a*JVS( 103 )
  a = -W( 40 ) / JVS( 134 )
  W( 40 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 135 )
  W( 85 ) = W( 85 ) + a*JVS( 136 )
  W( 86 ) = W( 86 ) + a*JVS( 137 )
  W( 89 ) = W( 89 ) + a*JVS( 138 )
  a = -W( 42 ) / JVS( 145 )
  W( 42 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 146 )
  W( 69 ) = W( 69 ) + a*JVS( 147 )
  W( 82 ) = W( 82 ) + a*JVS( 148 )
  W( 83 ) = W( 83 ) + a*JVS( 149 )
  W( 84 ) = W( 84 ) + a*JVS( 150 )
  W( 85 ) = W( 85 ) + a*JVS( 151 )
  W( 86 ) = W( 86 ) + a*JVS( 152 )
  W( 89 ) = W( 89 ) + a*JVS( 153 )
  a = -W( 47 ) / JVS( 179 )
  W( 47 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 180 )
  W( 88 ) = W( 88 ) + a*JVS( 181 )
  W( 89 ) = W( 89 ) + a*JVS( 182 )
  a = -W( 48 ) / JVS( 183 )
  W( 48 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 184 )
  W( 69 ) = W( 69 ) + a*JVS( 185 )
  W( 83 ) = W( 83 ) + a*JVS( 186 )
  W( 84 ) = W( 84 ) + a*JVS( 187 )
  a = -W( 55 ) / JVS( 225 )
  W( 55 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 226 )
  W( 79 ) = W( 79 ) + a*JVS( 227 )
  W( 83 ) = W( 83 ) + a*JVS( 228 )
  W( 84 ) = W( 84 ) + a*JVS( 229 )
  a = -W( 58 ) / JVS( 246 )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 247 )
  W( 66 ) = W( 66 ) + a*JVS( 248 )
  W( 71 ) = W( 71 ) + a*JVS( 249 )
  W( 82 ) = W( 82 ) + a*JVS( 250 )
  W( 84 ) = W( 84 ) + a*JVS( 251 )
  W( 86 ) = W( 86 ) + a*JVS( 252 )
  W( 87 ) = W( 87 ) + a*JVS( 253 )
  a = -W( 63 ) / JVS( 286 )
  W( 63 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 287 )
  W( 86 ) = W( 86 ) + a*JVS( 288 )
  W( 87 ) = W( 87 ) + a*JVS( 289 )
  a = -W( 64 ) / JVS( 295 )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 296 )
  W( 84 ) = W( 84 ) + a*JVS( 297 )
  W( 86 ) = W( 86 ) + a*JVS( 298 )
  W( 87 ) = W( 87 ) + a*JVS( 299 )
  a = -W( 66 ) / JVS( 326 )
  W( 66 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 327 )
  W( 84 ) = W( 84 ) + a*JVS( 328 )
  W( 86 ) = W( 86 ) + a*JVS( 329 )
  W( 87 ) = W( 87 ) + a*JVS( 330 )
  a = -W( 67 ) / JVS( 331 )
  W( 67 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 332 )
  W( 83 ) = W( 83 ) + a*JVS( 333 )
  W( 84 ) = W( 84 ) + a*JVS( 334 )
  a = -W( 69 ) / JVS( 346 )
  W( 69 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  W( 83 ) = W( 83 ) + a*JVS( 349 )
  W( 84 ) = W( 84 ) + a*JVS( 350 )
  a = -W( 70 ) / JVS( 363 )
  W( 70 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 364 )
  W( 81 ) = W( 81 ) + a*JVS( 365 )
  W( 82 ) = W( 82 ) + a*JVS( 366 )
  W( 83 ) = W( 83 ) + a*JVS( 367 )
  W( 84 ) = W( 84 ) + a*JVS( 368 )
  W( 86 ) = W( 86 ) + a*JVS( 369 )
  W( 87 ) = W( 87 ) + a*JVS( 370 )
  a = -W( 71 ) / JVS( 373 )
  W( 71 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 86 ) = W( 86 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  a = -W( 72 ) / JVS( 384 )
  W( 72 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 385 )
  W( 82 ) = W( 82 ) + a*JVS( 386 )
  W( 83 ) = W( 83 ) + a*JVS( 387 )
  W( 84 ) = W( 84 ) + a*JVS( 388 )
  W( 86 ) = W( 86 ) + a*JVS( 389 )
  W( 87 ) = W( 87 ) + a*JVS( 390 )
  a = -W( 76 ) / JVS( 451 )
  W( 76 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 452 )
  W( 80 ) = W( 80 ) + a*JVS( 453 )
  W( 81 ) = W( 81 ) + a*JVS( 454 )
  W( 82 ) = W( 82 ) + a*JVS( 455 )
  W( 83 ) = W( 83 ) + a*JVS( 456 )
  W( 84 ) = W( 84 ) + a*JVS( 457 )
  W( 86 ) = W( 86 ) + a*JVS( 458 )
  W( 87 ) = W( 87 ) + a*JVS( 459 )
  W( 89 ) = W( 89 ) + a*JVS( 460 )
  a = -W( 77 ) / JVS( 467 )
  W( 77 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 468 )
  W( 81 ) = W( 81 ) + a*JVS( 469 )
  W( 82 ) = W( 82 ) + a*JVS( 470 )
  W( 83 ) = W( 83 ) + a*JVS( 471 )
  W( 84 ) = W( 84 ) + a*JVS( 472 )
  W( 86 ) = W( 86 ) + a*JVS( 473 )
  W( 87 ) = W( 87 ) + a*JVS( 474 )
  W( 89 ) = W( 89 ) + a*JVS( 475 )
  a = -W( 78 ) / JVS( 479 )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 480 )
  W( 81 ) = W( 81 ) + a*JVS( 481 )
  W( 82 ) = W( 82 ) + a*JVS( 482 )
  W( 83 ) = W( 83 ) + a*JVS( 483 )
  W( 84 ) = W( 84 ) + a*JVS( 484 )
  W( 86 ) = W( 86 ) + a*JVS( 485 )
  W( 87 ) = W( 87 ) + a*JVS( 486 )
  W( 89 ) = W( 89 ) + a*JVS( 487 )
  a = -W( 79 ) / JVS( 518 )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 519 )
  W( 81 ) = W( 81 ) + a*JVS( 520 )
  W( 82 ) = W( 82 ) + a*JVS( 521 )
  W( 83 ) = W( 83 ) + a*JVS( 522 )
  W( 84 ) = W( 84 ) + a*JVS( 523 )
  W( 85 ) = W( 85 ) + a*JVS( 524 )
  W( 86 ) = W( 86 ) + a*JVS( 525 )
  W( 87 ) = W( 87 ) + a*JVS( 526 )
  W( 88 ) = W( 88 ) + a*JVS( 527 )
  W( 89 ) = W( 89 ) + a*JVS( 528 )
  a = -W( 80 ) / JVS( 531 )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 532 )
  W( 82 ) = W( 82 ) + a*JVS( 533 )
  W( 83 ) = W( 83 ) + a*JVS( 534 )
  W( 84 ) = W( 84 ) + a*JVS( 535 )
  W( 86 ) = W( 86 ) + a*JVS( 536 )
  W( 87 ) = W( 87 ) + a*JVS( 537 )
  W( 89 ) = W( 89 ) + a*JVS( 538 )
  a = -W( 81 ) / JVS( 557 )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 558 )
  W( 83 ) = W( 83 ) + a*JVS( 559 )
  W( 84 ) = W( 84 ) + a*JVS( 560 )
  W( 85 ) = W( 85 ) + a*JVS( 561 )
  W( 86 ) = W( 86 ) + a*JVS( 562 )
  W( 87 ) = W( 87 ) + a*JVS( 563 )
  W( 88 ) = W( 88 ) + a*JVS( 564 )
  W( 89 ) = W( 89 ) + a*JVS( 565 )
  JVS( 566) = W( 33 )
  JVS( 567) = W( 40 )
  JVS( 568) = W( 42 )
  JVS( 569) = W( 47 )
  JVS( 570) = W( 48 )
  JVS( 571) = W( 55 )
  JVS( 572) = W( 58 )
  JVS( 573) = W( 63 )
  JVS( 574) = W( 64 )
  JVS( 575) = W( 66 )
  JVS( 576) = W( 67 )
  JVS( 577) = W( 69 )
  JVS( 578) = W( 70 )
  JVS( 579) = W( 71 )
  JVS( 580) = W( 72 )
  JVS( 581) = W( 76 )
  JVS( 582) = W( 77 )
  JVS( 583) = W( 78 )
  JVS( 584) = W( 79 )
  JVS( 585) = W( 80 )
  JVS( 586) = W( 81 )
  JVS( 587) = W( 82 )
  JVS( 588) = W( 83 )
  JVS( 589) = W( 84 )
  JVS( 590) = W( 85 )
  JVS( 591) = W( 86 )
  JVS( 592) = W( 87 )
  JVS( 593) = W( 88 )
  JVS( 594) = W( 89 )
  IF ( ABS( JVS( 606 )) < TINY(a) ) THEN
         IER = 83
         RETURN
  END IF
   W( 36 ) = JVS( 595 )
   W( 54 ) = JVS( 596 )
   W( 60 ) = JVS( 597 )
   W( 67 ) = JVS( 598 )
   W( 69 ) = JVS( 599 )
   W( 74 ) = JVS( 600 )
   W( 77 ) = JVS( 601 )
   W( 79 ) = JVS( 602 )
   W( 80 ) = JVS( 603 )
   W( 81 ) = JVS( 604 )
   W( 82 ) = JVS( 605 )
   W( 83 ) = JVS( 606 )
   W( 84 ) = JVS( 607 )
   W( 85 ) = JVS( 608 )
   W( 86 ) = JVS( 609 )
   W( 87 ) = JVS( 610 )
   W( 88 ) = JVS( 611 )
   W( 89 ) = JVS( 612 )
  a = -W( 36 ) / JVS( 116 )
  W( 36 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 117 )
  W( 84 ) = W( 84 ) + a*JVS( 118 )
  a = -W( 54 ) / JVS( 218 )
  W( 54 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 219 )
  W( 83 ) = W( 83 ) + a*JVS( 220 )
  W( 84 ) = W( 84 ) + a*JVS( 221 )
  a = -W( 60 ) / JVS( 261 )
  W( 60 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 262 )
  W( 79 ) = W( 79 ) + a*JVS( 263 )
  W( 83 ) = W( 83 ) + a*JVS( 264 )
  W( 84 ) = W( 84 ) + a*JVS( 265 )
  W( 86 ) = W( 86 ) + a*JVS( 266 )
  W( 88 ) = W( 88 ) + a*JVS( 267 )
  a = -W( 67 ) / JVS( 331 )
  W( 67 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 332 )
  W( 83 ) = W( 83 ) + a*JVS( 333 )
  W( 84 ) = W( 84 ) + a*JVS( 334 )
  a = -W( 69 ) / JVS( 346 )
  W( 69 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  W( 83 ) = W( 83 ) + a*JVS( 349 )
  W( 84 ) = W( 84 ) + a*JVS( 350 )
  a = -W( 74 ) / JVS( 415 )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 416 )
  W( 81 ) = W( 81 ) + a*JVS( 417 )
  W( 82 ) = W( 82 ) + a*JVS( 418 )
  W( 83 ) = W( 83 ) + a*JVS( 419 )
  W( 84 ) = W( 84 ) + a*JVS( 420 )
  W( 86 ) = W( 86 ) + a*JVS( 421 )
  W( 87 ) = W( 87 ) + a*JVS( 422 )
  W( 89 ) = W( 89 ) + a*JVS( 423 )
  a = -W( 77 ) / JVS( 467 )
  W( 77 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 468 )
  W( 81 ) = W( 81 ) + a*JVS( 469 )
  W( 82 ) = W( 82 ) + a*JVS( 470 )
  W( 83 ) = W( 83 ) + a*JVS( 471 )
  W( 84 ) = W( 84 ) + a*JVS( 472 )
  W( 86 ) = W( 86 ) + a*JVS( 473 )
  W( 87 ) = W( 87 ) + a*JVS( 474 )
  W( 89 ) = W( 89 ) + a*JVS( 475 )
  a = -W( 79 ) / JVS( 518 )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 519 )
  W( 81 ) = W( 81 ) + a*JVS( 520 )
  W( 82 ) = W( 82 ) + a*JVS( 521 )
  W( 83 ) = W( 83 ) + a*JVS( 522 )
  W( 84 ) = W( 84 ) + a*JVS( 523 )
  W( 85 ) = W( 85 ) + a*JVS( 524 )
  W( 86 ) = W( 86 ) + a*JVS( 525 )
  W( 87 ) = W( 87 ) + a*JVS( 526 )
  W( 88 ) = W( 88 ) + a*JVS( 527 )
  W( 89 ) = W( 89 ) + a*JVS( 528 )
  a = -W( 80 ) / JVS( 531 )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 532 )
  W( 82 ) = W( 82 ) + a*JVS( 533 )
  W( 83 ) = W( 83 ) + a*JVS( 534 )
  W( 84 ) = W( 84 ) + a*JVS( 535 )
  W( 86 ) = W( 86 ) + a*JVS( 536 )
  W( 87 ) = W( 87 ) + a*JVS( 537 )
  W( 89 ) = W( 89 ) + a*JVS( 538 )
  a = -W( 81 ) / JVS( 557 )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 558 )
  W( 83 ) = W( 83 ) + a*JVS( 559 )
  W( 84 ) = W( 84 ) + a*JVS( 560 )
  W( 85 ) = W( 85 ) + a*JVS( 561 )
  W( 86 ) = W( 86 ) + a*JVS( 562 )
  W( 87 ) = W( 87 ) + a*JVS( 563 )
  W( 88 ) = W( 88 ) + a*JVS( 564 )
  W( 89 ) = W( 89 ) + a*JVS( 565 )
  a = -W( 82 ) / JVS( 587 )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 588 )
  W( 84 ) = W( 84 ) + a*JVS( 589 )
  W( 85 ) = W( 85 ) + a*JVS( 590 )
  W( 86 ) = W( 86 ) + a*JVS( 591 )
  W( 87 ) = W( 87 ) + a*JVS( 592 )
  W( 88 ) = W( 88 ) + a*JVS( 593 )
  W( 89 ) = W( 89 ) + a*JVS( 594 )
  JVS( 595) = W( 36 )
  JVS( 596) = W( 54 )
  JVS( 597) = W( 60 )
  JVS( 598) = W( 67 )
  JVS( 599) = W( 69 )
  JVS( 600) = W( 74 )
  JVS( 601) = W( 77 )
  JVS( 602) = W( 79 )
  JVS( 603) = W( 80 )
  JVS( 604) = W( 81 )
  JVS( 605) = W( 82 )
  JVS( 606) = W( 83 )
  JVS( 607) = W( 84 )
  JVS( 608) = W( 85 )
  JVS( 609) = W( 86 )
  JVS( 610) = W( 87 )
  JVS( 611) = W( 88 )
  JVS( 612) = W( 89 )
  IF ( ABS( JVS( 680 )) < TINY(a) ) THEN
         IER = 84
         RETURN
  END IF
   W( 10 ) = JVS( 613 )
   W( 11 ) = JVS( 614 )
   W( 12 ) = JVS( 615 )
   W( 13 ) = JVS( 616 )
   W( 14 ) = JVS( 617 )
   W( 16 ) = JVS( 618 )
   W( 17 ) = JVS( 619 )
   W( 18 ) = JVS( 620 )
   W( 19 ) = JVS( 621 )
   W( 23 ) = JVS( 622 )
   W( 24 ) = JVS( 623 )
   W( 25 ) = JVS( 624 )
   W( 26 ) = JVS( 625 )
   W( 27 ) = JVS( 626 )
   W( 28 ) = JVS( 627 )
   W( 29 ) = JVS( 628 )
   W( 30 ) = JVS( 629 )
   W( 32 ) = JVS( 630 )
   W( 33 ) = JVS( 631 )
   W( 34 ) = JVS( 632 )
   W( 35 ) = JVS( 633 )
   W( 36 ) = JVS( 634 )
   W( 37 ) = JVS( 635 )
   W( 38 ) = JVS( 636 )
   W( 40 ) = JVS( 637 )
   W( 41 ) = JVS( 638 )
   W( 42 ) = JVS( 639 )
   W( 43 ) = JVS( 640 )
   W( 44 ) = JVS( 641 )
   W( 45 ) = JVS( 642 )
   W( 46 ) = JVS( 643 )
   W( 47 ) = JVS( 644 )
   W( 48 ) = JVS( 645 )
   W( 49 ) = JVS( 646 )
   W( 50 ) = JVS( 647 )
   W( 51 ) = JVS( 648 )
   W( 52 ) = JVS( 649 )
   W( 53 ) = JVS( 650 )
   W( 54 ) = JVS( 651 )
   W( 55 ) = JVS( 652 )
   W( 56 ) = JVS( 653 )
   W( 57 ) = JVS( 654 )
   W( 58 ) = JVS( 655 )
   W( 59 ) = JVS( 656 )
   W( 60 ) = JVS( 657 )
   W( 62 ) = JVS( 658 )
   W( 63 ) = JVS( 659 )
   W( 64 ) = JVS( 660 )
   W( 65 ) = JVS( 661 )
   W( 66 ) = JVS( 662 )
   W( 67 ) = JVS( 663 )
   W( 68 ) = JVS( 664 )
   W( 69 ) = JVS( 665 )
   W( 70 ) = JVS( 666 )
   W( 71 ) = JVS( 667 )
   W( 72 ) = JVS( 668 )
   W( 73 ) = JVS( 669 )
   W( 74 ) = JVS( 670 )
   W( 75 ) = JVS( 671 )
   W( 76 ) = JVS( 672 )
   W( 77 ) = JVS( 673 )
   W( 78 ) = JVS( 674 )
   W( 79 ) = JVS( 675 )
   W( 80 ) = JVS( 676 )
   W( 81 ) = JVS( 677 )
   W( 82 ) = JVS( 678 )
   W( 83 ) = JVS( 679 )
   W( 84 ) = JVS( 680 )
   W( 85 ) = JVS( 681 )
   W( 86 ) = JVS( 682 )
   W( 87 ) = JVS( 683 )
   W( 88 ) = JVS( 684 )
   W( 89 ) = JVS( 685 )
  a = -W( 10 ) / JVS( 26 )
  W( 10 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 27 )
  a = -W( 11 ) / JVS( 28 )
  W( 11 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 29 )
  a = -W( 12 ) / JVS( 30 )
  W( 12 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 31 )
  a = -W( 13 ) / JVS( 32 )
  W( 13 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 33 )
  a = -W( 14 ) / JVS( 34 )
  W( 14 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 35 )
  a = -W( 16 ) / JVS( 38 )
  W( 16 ) = -a
  W( 24 ) = W( 24 ) + a*JVS( 39 )
  W( 81 ) = W( 81 ) + a*JVS( 40 )
  W( 84 ) = W( 84 ) + a*JVS( 41 )
  a = -W( 17 ) / JVS( 42 )
  W( 17 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 43 )
  a = -W( 18 ) / JVS( 45 )
  W( 18 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 46 )
  a = -W( 19 ) / JVS( 47 )
  W( 19 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 48 )
  W( 86 ) = W( 86 ) + a*JVS( 49 )
  a = -W( 23 ) / JVS( 60 )
  W( 23 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 61 )
  W( 84 ) = W( 84 ) + a*JVS( 62 )
  W( 86 ) = W( 86 ) + a*JVS( 63 )
  a = -W( 24 ) / JVS( 64 )
  W( 24 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 65 )
  W( 84 ) = W( 84 ) + a*JVS( 66 )
  a = -W( 25 ) / JVS( 67 )
  W( 25 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 68 )
  W( 82 ) = W( 82 ) + a*JVS( 69 )
  W( 84 ) = W( 84 ) + a*JVS( 70 )
  a = -W( 26 ) / JVS( 71 )
  W( 26 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 72 )
  W( 84 ) = W( 84 ) + a*JVS( 73 )
  W( 86 ) = W( 86 ) + a*JVS( 74 )
  a = -W( 27 ) / JVS( 75 )
  W( 27 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 76 )
  W( 55 ) = W( 55 ) + a*JVS( 77 )
  W( 79 ) = W( 79 ) + a*JVS( 78 )
  W( 84 ) = W( 84 ) + a*JVS( 79 )
  a = -W( 28 ) / JVS( 80 )
  W( 28 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 81 )
  W( 84 ) = W( 84 ) + a*JVS( 82 )
  W( 86 ) = W( 86 ) + a*JVS( 83 )
  a = -W( 29 ) / JVS( 84 )
  W( 29 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 85 )
  W( 84 ) = W( 84 ) + a*JVS( 86 )
  W( 86 ) = W( 86 ) + a*JVS( 87 )
  a = -W( 30 ) / JVS( 88 )
  W( 30 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 89 )
  W( 84 ) = W( 84 ) + a*JVS( 90 )
  W( 86 ) = W( 86 ) + a*JVS( 91 )
  a = -W( 32 ) / JVS( 96 )
  W( 32 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 97 )
  W( 84 ) = W( 84 ) + a*JVS( 98 )
  W( 86 ) = W( 86 ) + a*JVS( 99 )
  a = -W( 33 ) / JVS( 100 )
  W( 33 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 101 )
  W( 84 ) = W( 84 ) + a*JVS( 102 )
  W( 86 ) = W( 86 ) + a*JVS( 103 )
  a = -W( 34 ) / JVS( 104 )
  W( 34 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 105 )
  W( 86 ) = W( 86 ) + a*JVS( 106 )
  W( 88 ) = W( 88 ) + a*JVS( 107 )
  a = -W( 35 ) / JVS( 108 )
  W( 35 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 109 )
  W( 80 ) = W( 80 ) + a*JVS( 110 )
  W( 81 ) = W( 81 ) + a*JVS( 111 )
  W( 82 ) = W( 82 ) + a*JVS( 112 )
  W( 84 ) = W( 84 ) + a*JVS( 113 )
  W( 87 ) = W( 87 ) + a*JVS( 114 )
  W( 89 ) = W( 89 ) + a*JVS( 115 )
  a = -W( 36 ) / JVS( 116 )
  W( 36 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 117 )
  W( 84 ) = W( 84 ) + a*JVS( 118 )
  a = -W( 37 ) / JVS( 119 )
  W( 37 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 120 )
  W( 69 ) = W( 69 ) + a*JVS( 121 )
  W( 81 ) = W( 81 ) + a*JVS( 122 )
  W( 84 ) = W( 84 ) + a*JVS( 123 )
  W( 87 ) = W( 87 ) + a*JVS( 124 )
  a = -W( 38 ) / JVS( 125 )
  W( 38 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 126 )
  W( 84 ) = W( 84 ) + a*JVS( 127 )
  W( 86 ) = W( 86 ) + a*JVS( 128 )
  a = -W( 40 ) / JVS( 134 )
  W( 40 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 135 )
  W( 85 ) = W( 85 ) + a*JVS( 136 )
  W( 86 ) = W( 86 ) + a*JVS( 137 )
  W( 89 ) = W( 89 ) + a*JVS( 138 )
  a = -W( 41 ) / JVS( 139 )
  W( 41 ) = -a
  W( 56 ) = W( 56 ) + a*JVS( 140 )
  W( 63 ) = W( 63 ) + a*JVS( 141 )
  W( 84 ) = W( 84 ) + a*JVS( 142 )
  W( 87 ) = W( 87 ) + a*JVS( 143 )
  a = -W( 42 ) / JVS( 145 )
  W( 42 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 146 )
  W( 69 ) = W( 69 ) + a*JVS( 147 )
  W( 82 ) = W( 82 ) + a*JVS( 148 )
  W( 83 ) = W( 83 ) + a*JVS( 149 )
  W( 84 ) = W( 84 ) + a*JVS( 150 )
  W( 85 ) = W( 85 ) + a*JVS( 151 )
  W( 86 ) = W( 86 ) + a*JVS( 152 )
  W( 89 ) = W( 89 ) + a*JVS( 153 )
  a = -W( 43 ) / JVS( 155 )
  W( 43 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 156 )
  W( 53 ) = W( 53 ) + a*JVS( 157 )
  W( 62 ) = W( 62 ) + a*JVS( 158 )
  W( 84 ) = W( 84 ) + a*JVS( 159 )
  W( 86 ) = W( 86 ) + a*JVS( 160 )
  W( 87 ) = W( 87 ) + a*JVS( 161 )
  a = -W( 44 ) / JVS( 162 )
  W( 44 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 163 )
  W( 84 ) = W( 84 ) + a*JVS( 164 )
  W( 86 ) = W( 86 ) + a*JVS( 165 )
  a = -W( 45 ) / JVS( 168 )
  W( 45 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 169 )
  W( 75 ) = W( 75 ) + a*JVS( 170 )
  W( 79 ) = W( 79 ) + a*JVS( 171 )
  W( 81 ) = W( 81 ) + a*JVS( 172 )
  W( 84 ) = W( 84 ) + a*JVS( 173 )
  W( 88 ) = W( 88 ) + a*JVS( 174 )
  a = -W( 46 ) / JVS( 175 )
  W( 46 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 176 )
  W( 84 ) = W( 84 ) + a*JVS( 177 )
  W( 86 ) = W( 86 ) + a*JVS( 178 )
  a = -W( 47 ) / JVS( 179 )
  W( 47 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 180 )
  W( 88 ) = W( 88 ) + a*JVS( 181 )
  W( 89 ) = W( 89 ) + a*JVS( 182 )
  a = -W( 48 ) / JVS( 183 )
  W( 48 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 184 )
  W( 69 ) = W( 69 ) + a*JVS( 185 )
  W( 83 ) = W( 83 ) + a*JVS( 186 )
  W( 84 ) = W( 84 ) + a*JVS( 187 )
  a = -W( 49 ) / JVS( 188 )
  W( 49 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 189 )
  W( 84 ) = W( 84 ) + a*JVS( 190 )
  W( 86 ) = W( 86 ) + a*JVS( 191 )
  a = -W( 50 ) / JVS( 192 )
  W( 50 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 193 )
  W( 85 ) = W( 85 ) + a*JVS( 194 )
  W( 88 ) = W( 88 ) + a*JVS( 195 )
  a = -W( 51 ) / JVS( 196 )
  W( 51 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 197 )
  W( 72 ) = W( 72 ) + a*JVS( 198 )
  W( 76 ) = W( 76 ) + a*JVS( 199 )
  W( 78 ) = W( 78 ) + a*JVS( 200 )
  W( 80 ) = W( 80 ) + a*JVS( 201 )
  W( 82 ) = W( 82 ) + a*JVS( 202 )
  W( 84 ) = W( 84 ) + a*JVS( 203 )
  a = -W( 52 ) / JVS( 206 )
  W( 52 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 207 )
  W( 86 ) = W( 86 ) + a*JVS( 208 )
  W( 87 ) = W( 87 ) + a*JVS( 209 )
  a = -W( 53 ) / JVS( 213 )
  W( 53 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 214 )
  W( 86 ) = W( 86 ) + a*JVS( 215 )
  W( 87 ) = W( 87 ) + a*JVS( 216 )
  W( 88 ) = W( 88 ) + a*JVS( 217 )
  a = -W( 54 ) / JVS( 218 )
  W( 54 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 219 )
  W( 83 ) = W( 83 ) + a*JVS( 220 )
  W( 84 ) = W( 84 ) + a*JVS( 221 )
  a = -W( 55 ) / JVS( 225 )
  W( 55 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 226 )
  W( 79 ) = W( 79 ) + a*JVS( 227 )
  W( 83 ) = W( 83 ) + a*JVS( 228 )
  W( 84 ) = W( 84 ) + a*JVS( 229 )
  a = -W( 56 ) / JVS( 230 )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 231 )
  W( 84 ) = W( 84 ) + a*JVS( 232 )
  W( 86 ) = W( 86 ) + a*JVS( 233 )
  a = -W( 57 ) / JVS( 237 )
  W( 57 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 238 )
  W( 84 ) = W( 84 ) + a*JVS( 239 )
  W( 86 ) = W( 86 ) + a*JVS( 240 )
  W( 87 ) = W( 87 ) + a*JVS( 241 )
  a = -W( 58 ) / JVS( 246 )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 247 )
  W( 66 ) = W( 66 ) + a*JVS( 248 )
  W( 71 ) = W( 71 ) + a*JVS( 249 )
  W( 82 ) = W( 82 ) + a*JVS( 250 )
  W( 84 ) = W( 84 ) + a*JVS( 251 )
  W( 86 ) = W( 86 ) + a*JVS( 252 )
  W( 87 ) = W( 87 ) + a*JVS( 253 )
  a = -W( 59 ) / JVS( 255 )
  W( 59 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 256 )
  W( 84 ) = W( 84 ) + a*JVS( 257 )
  W( 86 ) = W( 86 ) + a*JVS( 258 )
  W( 87 ) = W( 87 ) + a*JVS( 259 )
  a = -W( 60 ) / JVS( 261 )
  W( 60 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 262 )
  W( 79 ) = W( 79 ) + a*JVS( 263 )
  W( 83 ) = W( 83 ) + a*JVS( 264 )
  W( 84 ) = W( 84 ) + a*JVS( 265 )
  W( 86 ) = W( 86 ) + a*JVS( 266 )
  W( 88 ) = W( 88 ) + a*JVS( 267 )
  a = -W( 62 ) / JVS( 275 )
  W( 62 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 276 )
  W( 78 ) = W( 78 ) + a*JVS( 277 )
  W( 81 ) = W( 81 ) + a*JVS( 278 )
  W( 82 ) = W( 82 ) + a*JVS( 279 )
  W( 83 ) = W( 83 ) + a*JVS( 280 )
  W( 84 ) = W( 84 ) + a*JVS( 281 )
  W( 87 ) = W( 87 ) + a*JVS( 282 )
  W( 89 ) = W( 89 ) + a*JVS( 283 )
  a = -W( 63 ) / JVS( 286 )
  W( 63 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 287 )
  W( 86 ) = W( 86 ) + a*JVS( 288 )
  W( 87 ) = W( 87 ) + a*JVS( 289 )
  a = -W( 64 ) / JVS( 295 )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 296 )
  W( 84 ) = W( 84 ) + a*JVS( 297 )
  W( 86 ) = W( 86 ) + a*JVS( 298 )
  W( 87 ) = W( 87 ) + a*JVS( 299 )
  a = -W( 65 ) / JVS( 305 )
  W( 65 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 306 )
  W( 68 ) = W( 68 ) + a*JVS( 307 )
  W( 69 ) = W( 69 ) + a*JVS( 308 )
  W( 70 ) = W( 70 ) + a*JVS( 309 )
  W( 74 ) = W( 74 ) + a*JVS( 310 )
  W( 75 ) = W( 75 ) + a*JVS( 311 )
  W( 76 ) = W( 76 ) + a*JVS( 312 )
  W( 77 ) = W( 77 ) + a*JVS( 313 )
  W( 78 ) = W( 78 ) + a*JVS( 314 )
  W( 79 ) = W( 79 ) + a*JVS( 315 )
  W( 81 ) = W( 81 ) + a*JVS( 316 )
  W( 82 ) = W( 82 ) + a*JVS( 317 )
  W( 83 ) = W( 83 ) + a*JVS( 318 )
  W( 84 ) = W( 84 ) + a*JVS( 319 )
  W( 86 ) = W( 86 ) + a*JVS( 320 )
  W( 87 ) = W( 87 ) + a*JVS( 321 )
  W( 88 ) = W( 88 ) + a*JVS( 322 )
  W( 89 ) = W( 89 ) + a*JVS( 323 )
  a = -W( 66 ) / JVS( 326 )
  W( 66 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 327 )
  W( 84 ) = W( 84 ) + a*JVS( 328 )
  W( 86 ) = W( 86 ) + a*JVS( 329 )
  W( 87 ) = W( 87 ) + a*JVS( 330 )
  a = -W( 67 ) / JVS( 331 )
  W( 67 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 332 )
  W( 83 ) = W( 83 ) + a*JVS( 333 )
  W( 84 ) = W( 84 ) + a*JVS( 334 )
  a = -W( 68 ) / JVS( 337 )
  W( 68 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 338 )
  W( 80 ) = W( 80 ) + a*JVS( 339 )
  W( 81 ) = W( 81 ) + a*JVS( 340 )
  W( 83 ) = W( 83 ) + a*JVS( 341 )
  W( 84 ) = W( 84 ) + a*JVS( 342 )
  W( 86 ) = W( 86 ) + a*JVS( 343 )
  W( 87 ) = W( 87 ) + a*JVS( 344 )
  a = -W( 69 ) / JVS( 346 )
  W( 69 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  W( 83 ) = W( 83 ) + a*JVS( 349 )
  W( 84 ) = W( 84 ) + a*JVS( 350 )
  a = -W( 70 ) / JVS( 363 )
  W( 70 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 364 )
  W( 81 ) = W( 81 ) + a*JVS( 365 )
  W( 82 ) = W( 82 ) + a*JVS( 366 )
  W( 83 ) = W( 83 ) + a*JVS( 367 )
  W( 84 ) = W( 84 ) + a*JVS( 368 )
  W( 86 ) = W( 86 ) + a*JVS( 369 )
  W( 87 ) = W( 87 ) + a*JVS( 370 )
  a = -W( 71 ) / JVS( 373 )
  W( 71 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 86 ) = W( 86 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  a = -W( 72 ) / JVS( 384 )
  W( 72 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 385 )
  W( 82 ) = W( 82 ) + a*JVS( 386 )
  W( 83 ) = W( 83 ) + a*JVS( 387 )
  W( 84 ) = W( 84 ) + a*JVS( 388 )
  W( 86 ) = W( 86 ) + a*JVS( 389 )
  W( 87 ) = W( 87 ) + a*JVS( 390 )
  a = -W( 73 ) / JVS( 396 )
  W( 73 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 397 )
  W( 77 ) = W( 77 ) + a*JVS( 398 )
  W( 78 ) = W( 78 ) + a*JVS( 399 )
  W( 81 ) = W( 81 ) + a*JVS( 400 )
  W( 82 ) = W( 82 ) + a*JVS( 401 )
  W( 83 ) = W( 83 ) + a*JVS( 402 )
  W( 84 ) = W( 84 ) + a*JVS( 403 )
  W( 85 ) = W( 85 ) + a*JVS( 404 )
  W( 86 ) = W( 86 ) + a*JVS( 405 )
  W( 87 ) = W( 87 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 89 ) = W( 89 ) + a*JVS( 408 )
  a = -W( 74 ) / JVS( 415 )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 416 )
  W( 81 ) = W( 81 ) + a*JVS( 417 )
  W( 82 ) = W( 82 ) + a*JVS( 418 )
  W( 83 ) = W( 83 ) + a*JVS( 419 )
  W( 84 ) = W( 84 ) + a*JVS( 420 )
  W( 86 ) = W( 86 ) + a*JVS( 421 )
  W( 87 ) = W( 87 ) + a*JVS( 422 )
  W( 89 ) = W( 89 ) + a*JVS( 423 )
  a = -W( 75 ) / JVS( 433 )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 434 )
  W( 77 ) = W( 77 ) + a*JVS( 435 )
  W( 78 ) = W( 78 ) + a*JVS( 436 )
  W( 80 ) = W( 80 ) + a*JVS( 437 )
  W( 81 ) = W( 81 ) + a*JVS( 438 )
  W( 82 ) = W( 82 ) + a*JVS( 439 )
  W( 83 ) = W( 83 ) + a*JVS( 440 )
  W( 84 ) = W( 84 ) + a*JVS( 441 )
  W( 85 ) = W( 85 ) + a*JVS( 442 )
  W( 86 ) = W( 86 ) + a*JVS( 443 )
  W( 87 ) = W( 87 ) + a*JVS( 444 )
  W( 88 ) = W( 88 ) + a*JVS( 445 )
  W( 89 ) = W( 89 ) + a*JVS( 446 )
  a = -W( 76 ) / JVS( 451 )
  W( 76 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 452 )
  W( 80 ) = W( 80 ) + a*JVS( 453 )
  W( 81 ) = W( 81 ) + a*JVS( 454 )
  W( 82 ) = W( 82 ) + a*JVS( 455 )
  W( 83 ) = W( 83 ) + a*JVS( 456 )
  W( 84 ) = W( 84 ) + a*JVS( 457 )
  W( 86 ) = W( 86 ) + a*JVS( 458 )
  W( 87 ) = W( 87 ) + a*JVS( 459 )
  W( 89 ) = W( 89 ) + a*JVS( 460 )
  a = -W( 77 ) / JVS( 467 )
  W( 77 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 468 )
  W( 81 ) = W( 81 ) + a*JVS( 469 )
  W( 82 ) = W( 82 ) + a*JVS( 470 )
  W( 83 ) = W( 83 ) + a*JVS( 471 )
  W( 84 ) = W( 84 ) + a*JVS( 472 )
  W( 86 ) = W( 86 ) + a*JVS( 473 )
  W( 87 ) = W( 87 ) + a*JVS( 474 )
  W( 89 ) = W( 89 ) + a*JVS( 475 )
  a = -W( 78 ) / JVS( 479 )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 480 )
  W( 81 ) = W( 81 ) + a*JVS( 481 )
  W( 82 ) = W( 82 ) + a*JVS( 482 )
  W( 83 ) = W( 83 ) + a*JVS( 483 )
  W( 84 ) = W( 84 ) + a*JVS( 484 )
  W( 86 ) = W( 86 ) + a*JVS( 485 )
  W( 87 ) = W( 87 ) + a*JVS( 486 )
  W( 89 ) = W( 89 ) + a*JVS( 487 )
  a = -W( 79 ) / JVS( 518 )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 519 )
  W( 81 ) = W( 81 ) + a*JVS( 520 )
  W( 82 ) = W( 82 ) + a*JVS( 521 )
  W( 83 ) = W( 83 ) + a*JVS( 522 )
  W( 84 ) = W( 84 ) + a*JVS( 523 )
  W( 85 ) = W( 85 ) + a*JVS( 524 )
  W( 86 ) = W( 86 ) + a*JVS( 525 )
  W( 87 ) = W( 87 ) + a*JVS( 526 )
  W( 88 ) = W( 88 ) + a*JVS( 527 )
  W( 89 ) = W( 89 ) + a*JVS( 528 )
  a = -W( 80 ) / JVS( 531 )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 532 )
  W( 82 ) = W( 82 ) + a*JVS( 533 )
  W( 83 ) = W( 83 ) + a*JVS( 534 )
  W( 84 ) = W( 84 ) + a*JVS( 535 )
  W( 86 ) = W( 86 ) + a*JVS( 536 )
  W( 87 ) = W( 87 ) + a*JVS( 537 )
  W( 89 ) = W( 89 ) + a*JVS( 538 )
  a = -W( 81 ) / JVS( 557 )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 558 )
  W( 83 ) = W( 83 ) + a*JVS( 559 )
  W( 84 ) = W( 84 ) + a*JVS( 560 )
  W( 85 ) = W( 85 ) + a*JVS( 561 )
  W( 86 ) = W( 86 ) + a*JVS( 562 )
  W( 87 ) = W( 87 ) + a*JVS( 563 )
  W( 88 ) = W( 88 ) + a*JVS( 564 )
  W( 89 ) = W( 89 ) + a*JVS( 565 )
  a = -W( 82 ) / JVS( 587 )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 588 )
  W( 84 ) = W( 84 ) + a*JVS( 589 )
  W( 85 ) = W( 85 ) + a*JVS( 590 )
  W( 86 ) = W( 86 ) + a*JVS( 591 )
  W( 87 ) = W( 87 ) + a*JVS( 592 )
  W( 88 ) = W( 88 ) + a*JVS( 593 )
  W( 89 ) = W( 89 ) + a*JVS( 594 )
  a = -W( 83 ) / JVS( 606 )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 607 )
  W( 85 ) = W( 85 ) + a*JVS( 608 )
  W( 86 ) = W( 86 ) + a*JVS( 609 )
  W( 87 ) = W( 87 ) + a*JVS( 610 )
  W( 88 ) = W( 88 ) + a*JVS( 611 )
  W( 89 ) = W( 89 ) + a*JVS( 612 )
  JVS( 613) = W( 10 )
  JVS( 614) = W( 11 )
  JVS( 615) = W( 12 )
  JVS( 616) = W( 13 )
  JVS( 617) = W( 14 )
  JVS( 618) = W( 16 )
  JVS( 619) = W( 17 )
  JVS( 620) = W( 18 )
  JVS( 621) = W( 19 )
  JVS( 622) = W( 23 )
  JVS( 623) = W( 24 )
  JVS( 624) = W( 25 )
  JVS( 625) = W( 26 )
  JVS( 626) = W( 27 )
  JVS( 627) = W( 28 )
  JVS( 628) = W( 29 )
  JVS( 629) = W( 30 )
  JVS( 630) = W( 32 )
  JVS( 631) = W( 33 )
  JVS( 632) = W( 34 )
  JVS( 633) = W( 35 )
  JVS( 634) = W( 36 )
  JVS( 635) = W( 37 )
  JVS( 636) = W( 38 )
  JVS( 637) = W( 40 )
  JVS( 638) = W( 41 )
  JVS( 639) = W( 42 )
  JVS( 640) = W( 43 )
  JVS( 641) = W( 44 )
  JVS( 642) = W( 45 )
  JVS( 643) = W( 46 )
  JVS( 644) = W( 47 )
  JVS( 645) = W( 48 )
  JVS( 646) = W( 49 )
  JVS( 647) = W( 50 )
  JVS( 648) = W( 51 )
  JVS( 649) = W( 52 )
  JVS( 650) = W( 53 )
  JVS( 651) = W( 54 )
  JVS( 652) = W( 55 )
  JVS( 653) = W( 56 )
  JVS( 654) = W( 57 )
  JVS( 655) = W( 58 )
  JVS( 656) = W( 59 )
  JVS( 657) = W( 60 )
  JVS( 658) = W( 62 )
  JVS( 659) = W( 63 )
  JVS( 660) = W( 64 )
  JVS( 661) = W( 65 )
  JVS( 662) = W( 66 )
  JVS( 663) = W( 67 )
  JVS( 664) = W( 68 )
  JVS( 665) = W( 69 )
  JVS( 666) = W( 70 )
  JVS( 667) = W( 71 )
  JVS( 668) = W( 72 )
  JVS( 669) = W( 73 )
  JVS( 670) = W( 74 )
  JVS( 671) = W( 75 )
  JVS( 672) = W( 76 )
  JVS( 673) = W( 77 )
  JVS( 674) = W( 78 )
  JVS( 675) = W( 79 )
  JVS( 676) = W( 80 )
  JVS( 677) = W( 81 )
  JVS( 678) = W( 82 )
  JVS( 679) = W( 83 )
  JVS( 680) = W( 84 )
  JVS( 681) = W( 85 )
  JVS( 682) = W( 86 )
  JVS( 683) = W( 87 )
  JVS( 684) = W( 88 )
  JVS( 685) = W( 89 )
  IF ( ABS( JVS( 696 )) < TINY(a) ) THEN
         IER = 85
         RETURN
  END IF
   W( 26 ) = JVS( 686 )
   W( 50 ) = JVS( 687 )
   W( 67 ) = JVS( 688 )
   W( 74 ) = JVS( 689 )
   W( 78 ) = JVS( 690 )
   W( 80 ) = JVS( 691 )
   W( 81 ) = JVS( 692 )
   W( 82 ) = JVS( 693 )
   W( 83 ) = JVS( 694 )
   W( 84 ) = JVS( 695 )
   W( 85 ) = JVS( 696 )
   W( 86 ) = JVS( 697 )
   W( 87 ) = JVS( 698 )
   W( 88 ) = JVS( 699 )
   W( 89 ) = JVS( 700 )
  a = -W( 26 ) / JVS( 71 )
  W( 26 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 72 )
  W( 84 ) = W( 84 ) + a*JVS( 73 )
  W( 86 ) = W( 86 ) + a*JVS( 74 )
  a = -W( 50 ) / JVS( 192 )
  W( 50 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 193 )
  W( 85 ) = W( 85 ) + a*JVS( 194 )
  W( 88 ) = W( 88 ) + a*JVS( 195 )
  a = -W( 67 ) / JVS( 331 )
  W( 67 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 332 )
  W( 83 ) = W( 83 ) + a*JVS( 333 )
  W( 84 ) = W( 84 ) + a*JVS( 334 )
  a = -W( 74 ) / JVS( 415 )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 416 )
  W( 81 ) = W( 81 ) + a*JVS( 417 )
  W( 82 ) = W( 82 ) + a*JVS( 418 )
  W( 83 ) = W( 83 ) + a*JVS( 419 )
  W( 84 ) = W( 84 ) + a*JVS( 420 )
  W( 86 ) = W( 86 ) + a*JVS( 421 )
  W( 87 ) = W( 87 ) + a*JVS( 422 )
  W( 89 ) = W( 89 ) + a*JVS( 423 )
  a = -W( 78 ) / JVS( 479 )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 480 )
  W( 81 ) = W( 81 ) + a*JVS( 481 )
  W( 82 ) = W( 82 ) + a*JVS( 482 )
  W( 83 ) = W( 83 ) + a*JVS( 483 )
  W( 84 ) = W( 84 ) + a*JVS( 484 )
  W( 86 ) = W( 86 ) + a*JVS( 485 )
  W( 87 ) = W( 87 ) + a*JVS( 486 )
  W( 89 ) = W( 89 ) + a*JVS( 487 )
  a = -W( 80 ) / JVS( 531 )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 532 )
  W( 82 ) = W( 82 ) + a*JVS( 533 )
  W( 83 ) = W( 83 ) + a*JVS( 534 )
  W( 84 ) = W( 84 ) + a*JVS( 535 )
  W( 86 ) = W( 86 ) + a*JVS( 536 )
  W( 87 ) = W( 87 ) + a*JVS( 537 )
  W( 89 ) = W( 89 ) + a*JVS( 538 )
  a = -W( 81 ) / JVS( 557 )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 558 )
  W( 83 ) = W( 83 ) + a*JVS( 559 )
  W( 84 ) = W( 84 ) + a*JVS( 560 )
  W( 85 ) = W( 85 ) + a*JVS( 561 )
  W( 86 ) = W( 86 ) + a*JVS( 562 )
  W( 87 ) = W( 87 ) + a*JVS( 563 )
  W( 88 ) = W( 88 ) + a*JVS( 564 )
  W( 89 ) = W( 89 ) + a*JVS( 565 )
  a = -W( 82 ) / JVS( 587 )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 588 )
  W( 84 ) = W( 84 ) + a*JVS( 589 )
  W( 85 ) = W( 85 ) + a*JVS( 590 )
  W( 86 ) = W( 86 ) + a*JVS( 591 )
  W( 87 ) = W( 87 ) + a*JVS( 592 )
  W( 88 ) = W( 88 ) + a*JVS( 593 )
  W( 89 ) = W( 89 ) + a*JVS( 594 )
  a = -W( 83 ) / JVS( 606 )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 607 )
  W( 85 ) = W( 85 ) + a*JVS( 608 )
  W( 86 ) = W( 86 ) + a*JVS( 609 )
  W( 87 ) = W( 87 ) + a*JVS( 610 )
  W( 88 ) = W( 88 ) + a*JVS( 611 )
  W( 89 ) = W( 89 ) + a*JVS( 612 )
  a = -W( 84 ) / JVS( 680 )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 681 )
  W( 86 ) = W( 86 ) + a*JVS( 682 )
  W( 87 ) = W( 87 ) + a*JVS( 683 )
  W( 88 ) = W( 88 ) + a*JVS( 684 )
  W( 89 ) = W( 89 ) + a*JVS( 685 )
  JVS( 686) = W( 26 )
  JVS( 687) = W( 50 )
  JVS( 688) = W( 67 )
  JVS( 689) = W( 74 )
  JVS( 690) = W( 78 )
  JVS( 691) = W( 80 )
  JVS( 692) = W( 81 )
  JVS( 693) = W( 82 )
  JVS( 694) = W( 83 )
  JVS( 695) = W( 84 )
  JVS( 696) = W( 85 )
  JVS( 697) = W( 86 )
  JVS( 698) = W( 87 )
  JVS( 699) = W( 88 )
  JVS( 700) = W( 89 )
  IF ( ABS( JVS( 756 )) < TINY(a) ) THEN
         IER = 86
         RETURN
  END IF
   W( 17 ) = JVS( 701 )
   W( 19 ) = JVS( 702 )
   W( 20 ) = JVS( 703 )
   W( 22 ) = JVS( 704 )
   W( 24 ) = JVS( 705 )
   W( 25 ) = JVS( 706 )
   W( 26 ) = JVS( 707 )
   W( 27 ) = JVS( 708 )
   W( 28 ) = JVS( 709 )
   W( 29 ) = JVS( 710 )
   W( 31 ) = JVS( 711 )
   W( 33 ) = JVS( 712 )
   W( 34 ) = JVS( 713 )
   W( 36 ) = JVS( 714 )
   W( 39 ) = JVS( 715 )
   W( 43 ) = JVS( 716 )
   W( 44 ) = JVS( 717 )
   W( 46 ) = JVS( 718 )
   W( 48 ) = JVS( 719 )
   W( 49 ) = JVS( 720 )
   W( 50 ) = JVS( 721 )
   W( 51 ) = JVS( 722 )
   W( 52 ) = JVS( 723 )
   W( 53 ) = JVS( 724 )
   W( 54 ) = JVS( 725 )
   W( 55 ) = JVS( 726 )
   W( 56 ) = JVS( 727 )
   W( 57 ) = JVS( 728 )
   W( 59 ) = JVS( 729 )
   W( 60 ) = JVS( 730 )
   W( 61 ) = JVS( 731 )
   W( 62 ) = JVS( 732 )
   W( 63 ) = JVS( 733 )
   W( 64 ) = JVS( 734 )
   W( 65 ) = JVS( 735 )
   W( 66 ) = JVS( 736 )
   W( 67 ) = JVS( 737 )
   W( 68 ) = JVS( 738 )
   W( 69 ) = JVS( 739 )
   W( 70 ) = JVS( 740 )
   W( 71 ) = JVS( 741 )
   W( 72 ) = JVS( 742 )
   W( 73 ) = JVS( 743 )
   W( 74 ) = JVS( 744 )
   W( 75 ) = JVS( 745 )
   W( 76 ) = JVS( 746 )
   W( 77 ) = JVS( 747 )
   W( 78 ) = JVS( 748 )
   W( 79 ) = JVS( 749 )
   W( 80 ) = JVS( 750 )
   W( 81 ) = JVS( 751 )
   W( 82 ) = JVS( 752 )
   W( 83 ) = JVS( 753 )
   W( 84 ) = JVS( 754 )
   W( 85 ) = JVS( 755 )
   W( 86 ) = JVS( 756 )
   W( 87 ) = JVS( 757 )
   W( 88 ) = JVS( 758 )
   W( 89 ) = JVS( 759 )
  a = -W( 17 ) / JVS( 42 )
  W( 17 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 43 )
  a = -W( 19 ) / JVS( 47 )
  W( 19 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 48 )
  W( 86 ) = W( 86 ) + a*JVS( 49 )
  a = -W( 20 ) / JVS( 50 )
  W( 20 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 51 )
  W( 87 ) = W( 87 ) + a*JVS( 52 )
  a = -W( 22 ) / JVS( 57 )
  W( 22 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 58 )
  W( 88 ) = W( 88 ) + a*JVS( 59 )
  a = -W( 24 ) / JVS( 64 )
  W( 24 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 65 )
  W( 84 ) = W( 84 ) + a*JVS( 66 )
  a = -W( 25 ) / JVS( 67 )
  W( 25 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 68 )
  W( 82 ) = W( 82 ) + a*JVS( 69 )
  W( 84 ) = W( 84 ) + a*JVS( 70 )
  a = -W( 26 ) / JVS( 71 )
  W( 26 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 72 )
  W( 84 ) = W( 84 ) + a*JVS( 73 )
  W( 86 ) = W( 86 ) + a*JVS( 74 )
  a = -W( 27 ) / JVS( 75 )
  W( 27 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 76 )
  W( 55 ) = W( 55 ) + a*JVS( 77 )
  W( 79 ) = W( 79 ) + a*JVS( 78 )
  W( 84 ) = W( 84 ) + a*JVS( 79 )
  a = -W( 28 ) / JVS( 80 )
  W( 28 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 81 )
  W( 84 ) = W( 84 ) + a*JVS( 82 )
  W( 86 ) = W( 86 ) + a*JVS( 83 )
  a = -W( 29 ) / JVS( 84 )
  W( 29 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 85 )
  W( 84 ) = W( 84 ) + a*JVS( 86 )
  W( 86 ) = W( 86 ) + a*JVS( 87 )
  a = -W( 31 ) / JVS( 93 )
  W( 31 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 94 )
  W( 87 ) = W( 87 ) + a*JVS( 95 )
  a = -W( 33 ) / JVS( 100 )
  W( 33 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 101 )
  W( 84 ) = W( 84 ) + a*JVS( 102 )
  W( 86 ) = W( 86 ) + a*JVS( 103 )
  a = -W( 34 ) / JVS( 104 )
  W( 34 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 105 )
  W( 86 ) = W( 86 ) + a*JVS( 106 )
  W( 88 ) = W( 88 ) + a*JVS( 107 )
  a = -W( 36 ) / JVS( 116 )
  W( 36 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 117 )
  W( 84 ) = W( 84 ) + a*JVS( 118 )
  a = -W( 39 ) / JVS( 130 )
  W( 39 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 131 )
  W( 84 ) = W( 84 ) + a*JVS( 132 )
  W( 87 ) = W( 87 ) + a*JVS( 133 )
  a = -W( 43 ) / JVS( 155 )
  W( 43 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 156 )
  W( 53 ) = W( 53 ) + a*JVS( 157 )
  W( 62 ) = W( 62 ) + a*JVS( 158 )
  W( 84 ) = W( 84 ) + a*JVS( 159 )
  W( 86 ) = W( 86 ) + a*JVS( 160 )
  W( 87 ) = W( 87 ) + a*JVS( 161 )
  a = -W( 44 ) / JVS( 162 )
  W( 44 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 163 )
  W( 84 ) = W( 84 ) + a*JVS( 164 )
  W( 86 ) = W( 86 ) + a*JVS( 165 )
  a = -W( 46 ) / JVS( 175 )
  W( 46 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 176 )
  W( 84 ) = W( 84 ) + a*JVS( 177 )
  W( 86 ) = W( 86 ) + a*JVS( 178 )
  a = -W( 48 ) / JVS( 183 )
  W( 48 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 184 )
  W( 69 ) = W( 69 ) + a*JVS( 185 )
  W( 83 ) = W( 83 ) + a*JVS( 186 )
  W( 84 ) = W( 84 ) + a*JVS( 187 )
  a = -W( 49 ) / JVS( 188 )
  W( 49 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 189 )
  W( 84 ) = W( 84 ) + a*JVS( 190 )
  W( 86 ) = W( 86 ) + a*JVS( 191 )
  a = -W( 50 ) / JVS( 192 )
  W( 50 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 193 )
  W( 85 ) = W( 85 ) + a*JVS( 194 )
  W( 88 ) = W( 88 ) + a*JVS( 195 )
  a = -W( 51 ) / JVS( 196 )
  W( 51 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 197 )
  W( 72 ) = W( 72 ) + a*JVS( 198 )
  W( 76 ) = W( 76 ) + a*JVS( 199 )
  W( 78 ) = W( 78 ) + a*JVS( 200 )
  W( 80 ) = W( 80 ) + a*JVS( 201 )
  W( 82 ) = W( 82 ) + a*JVS( 202 )
  W( 84 ) = W( 84 ) + a*JVS( 203 )
  a = -W( 52 ) / JVS( 206 )
  W( 52 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 207 )
  W( 86 ) = W( 86 ) + a*JVS( 208 )
  W( 87 ) = W( 87 ) + a*JVS( 209 )
  a = -W( 53 ) / JVS( 213 )
  W( 53 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 214 )
  W( 86 ) = W( 86 ) + a*JVS( 215 )
  W( 87 ) = W( 87 ) + a*JVS( 216 )
  W( 88 ) = W( 88 ) + a*JVS( 217 )
  a = -W( 54 ) / JVS( 218 )
  W( 54 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 219 )
  W( 83 ) = W( 83 ) + a*JVS( 220 )
  W( 84 ) = W( 84 ) + a*JVS( 221 )
  a = -W( 55 ) / JVS( 225 )
  W( 55 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 226 )
  W( 79 ) = W( 79 ) + a*JVS( 227 )
  W( 83 ) = W( 83 ) + a*JVS( 228 )
  W( 84 ) = W( 84 ) + a*JVS( 229 )
  a = -W( 56 ) / JVS( 230 )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 231 )
  W( 84 ) = W( 84 ) + a*JVS( 232 )
  W( 86 ) = W( 86 ) + a*JVS( 233 )
  a = -W( 57 ) / JVS( 237 )
  W( 57 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 238 )
  W( 84 ) = W( 84 ) + a*JVS( 239 )
  W( 86 ) = W( 86 ) + a*JVS( 240 )
  W( 87 ) = W( 87 ) + a*JVS( 241 )
  a = -W( 59 ) / JVS( 255 )
  W( 59 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 256 )
  W( 84 ) = W( 84 ) + a*JVS( 257 )
  W( 86 ) = W( 86 ) + a*JVS( 258 )
  W( 87 ) = W( 87 ) + a*JVS( 259 )
  a = -W( 60 ) / JVS( 261 )
  W( 60 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 262 )
  W( 79 ) = W( 79 ) + a*JVS( 263 )
  W( 83 ) = W( 83 ) + a*JVS( 264 )
  W( 84 ) = W( 84 ) + a*JVS( 265 )
  W( 86 ) = W( 86 ) + a*JVS( 266 )
  W( 88 ) = W( 88 ) + a*JVS( 267 )
  a = -W( 61 ) / JVS( 268 )
  W( 61 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 269 )
  W( 81 ) = W( 81 ) + a*JVS( 270 )
  W( 86 ) = W( 86 ) + a*JVS( 271 )
  W( 87 ) = W( 87 ) + a*JVS( 272 )
  a = -W( 62 ) / JVS( 275 )
  W( 62 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 276 )
  W( 78 ) = W( 78 ) + a*JVS( 277 )
  W( 81 ) = W( 81 ) + a*JVS( 278 )
  W( 82 ) = W( 82 ) + a*JVS( 279 )
  W( 83 ) = W( 83 ) + a*JVS( 280 )
  W( 84 ) = W( 84 ) + a*JVS( 281 )
  W( 87 ) = W( 87 ) + a*JVS( 282 )
  W( 89 ) = W( 89 ) + a*JVS( 283 )
  a = -W( 63 ) / JVS( 286 )
  W( 63 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 287 )
  W( 86 ) = W( 86 ) + a*JVS( 288 )
  W( 87 ) = W( 87 ) + a*JVS( 289 )
  a = -W( 64 ) / JVS( 295 )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 296 )
  W( 84 ) = W( 84 ) + a*JVS( 297 )
  W( 86 ) = W( 86 ) + a*JVS( 298 )
  W( 87 ) = W( 87 ) + a*JVS( 299 )
  a = -W( 65 ) / JVS( 305 )
  W( 65 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 306 )
  W( 68 ) = W( 68 ) + a*JVS( 307 )
  W( 69 ) = W( 69 ) + a*JVS( 308 )
  W( 70 ) = W( 70 ) + a*JVS( 309 )
  W( 74 ) = W( 74 ) + a*JVS( 310 )
  W( 75 ) = W( 75 ) + a*JVS( 311 )
  W( 76 ) = W( 76 ) + a*JVS( 312 )
  W( 77 ) = W( 77 ) + a*JVS( 313 )
  W( 78 ) = W( 78 ) + a*JVS( 314 )
  W( 79 ) = W( 79 ) + a*JVS( 315 )
  W( 81 ) = W( 81 ) + a*JVS( 316 )
  W( 82 ) = W( 82 ) + a*JVS( 317 )
  W( 83 ) = W( 83 ) + a*JVS( 318 )
  W( 84 ) = W( 84 ) + a*JVS( 319 )
  W( 86 ) = W( 86 ) + a*JVS( 320 )
  W( 87 ) = W( 87 ) + a*JVS( 321 )
  W( 88 ) = W( 88 ) + a*JVS( 322 )
  W( 89 ) = W( 89 ) + a*JVS( 323 )
  a = -W( 66 ) / JVS( 326 )
  W( 66 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 327 )
  W( 84 ) = W( 84 ) + a*JVS( 328 )
  W( 86 ) = W( 86 ) + a*JVS( 329 )
  W( 87 ) = W( 87 ) + a*JVS( 330 )
  a = -W( 67 ) / JVS( 331 )
  W( 67 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 332 )
  W( 83 ) = W( 83 ) + a*JVS( 333 )
  W( 84 ) = W( 84 ) + a*JVS( 334 )
  a = -W( 68 ) / JVS( 337 )
  W( 68 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 338 )
  W( 80 ) = W( 80 ) + a*JVS( 339 )
  W( 81 ) = W( 81 ) + a*JVS( 340 )
  W( 83 ) = W( 83 ) + a*JVS( 341 )
  W( 84 ) = W( 84 ) + a*JVS( 342 )
  W( 86 ) = W( 86 ) + a*JVS( 343 )
  W( 87 ) = W( 87 ) + a*JVS( 344 )
  a = -W( 69 ) / JVS( 346 )
  W( 69 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  W( 83 ) = W( 83 ) + a*JVS( 349 )
  W( 84 ) = W( 84 ) + a*JVS( 350 )
  a = -W( 70 ) / JVS( 363 )
  W( 70 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 364 )
  W( 81 ) = W( 81 ) + a*JVS( 365 )
  W( 82 ) = W( 82 ) + a*JVS( 366 )
  W( 83 ) = W( 83 ) + a*JVS( 367 )
  W( 84 ) = W( 84 ) + a*JVS( 368 )
  W( 86 ) = W( 86 ) + a*JVS( 369 )
  W( 87 ) = W( 87 ) + a*JVS( 370 )
  a = -W( 71 ) / JVS( 373 )
  W( 71 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 86 ) = W( 86 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  a = -W( 72 ) / JVS( 384 )
  W( 72 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 385 )
  W( 82 ) = W( 82 ) + a*JVS( 386 )
  W( 83 ) = W( 83 ) + a*JVS( 387 )
  W( 84 ) = W( 84 ) + a*JVS( 388 )
  W( 86 ) = W( 86 ) + a*JVS( 389 )
  W( 87 ) = W( 87 ) + a*JVS( 390 )
  a = -W( 73 ) / JVS( 396 )
  W( 73 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 397 )
  W( 77 ) = W( 77 ) + a*JVS( 398 )
  W( 78 ) = W( 78 ) + a*JVS( 399 )
  W( 81 ) = W( 81 ) + a*JVS( 400 )
  W( 82 ) = W( 82 ) + a*JVS( 401 )
  W( 83 ) = W( 83 ) + a*JVS( 402 )
  W( 84 ) = W( 84 ) + a*JVS( 403 )
  W( 85 ) = W( 85 ) + a*JVS( 404 )
  W( 86 ) = W( 86 ) + a*JVS( 405 )
  W( 87 ) = W( 87 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 89 ) = W( 89 ) + a*JVS( 408 )
  a = -W( 74 ) / JVS( 415 )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 416 )
  W( 81 ) = W( 81 ) + a*JVS( 417 )
  W( 82 ) = W( 82 ) + a*JVS( 418 )
  W( 83 ) = W( 83 ) + a*JVS( 419 )
  W( 84 ) = W( 84 ) + a*JVS( 420 )
  W( 86 ) = W( 86 ) + a*JVS( 421 )
  W( 87 ) = W( 87 ) + a*JVS( 422 )
  W( 89 ) = W( 89 ) + a*JVS( 423 )
  a = -W( 75 ) / JVS( 433 )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 434 )
  W( 77 ) = W( 77 ) + a*JVS( 435 )
  W( 78 ) = W( 78 ) + a*JVS( 436 )
  W( 80 ) = W( 80 ) + a*JVS( 437 )
  W( 81 ) = W( 81 ) + a*JVS( 438 )
  W( 82 ) = W( 82 ) + a*JVS( 439 )
  W( 83 ) = W( 83 ) + a*JVS( 440 )
  W( 84 ) = W( 84 ) + a*JVS( 441 )
  W( 85 ) = W( 85 ) + a*JVS( 442 )
  W( 86 ) = W( 86 ) + a*JVS( 443 )
  W( 87 ) = W( 87 ) + a*JVS( 444 )
  W( 88 ) = W( 88 ) + a*JVS( 445 )
  W( 89 ) = W( 89 ) + a*JVS( 446 )
  a = -W( 76 ) / JVS( 451 )
  W( 76 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 452 )
  W( 80 ) = W( 80 ) + a*JVS( 453 )
  W( 81 ) = W( 81 ) + a*JVS( 454 )
  W( 82 ) = W( 82 ) + a*JVS( 455 )
  W( 83 ) = W( 83 ) + a*JVS( 456 )
  W( 84 ) = W( 84 ) + a*JVS( 457 )
  W( 86 ) = W( 86 ) + a*JVS( 458 )
  W( 87 ) = W( 87 ) + a*JVS( 459 )
  W( 89 ) = W( 89 ) + a*JVS( 460 )
  a = -W( 77 ) / JVS( 467 )
  W( 77 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 468 )
  W( 81 ) = W( 81 ) + a*JVS( 469 )
  W( 82 ) = W( 82 ) + a*JVS( 470 )
  W( 83 ) = W( 83 ) + a*JVS( 471 )
  W( 84 ) = W( 84 ) + a*JVS( 472 )
  W( 86 ) = W( 86 ) + a*JVS( 473 )
  W( 87 ) = W( 87 ) + a*JVS( 474 )
  W( 89 ) = W( 89 ) + a*JVS( 475 )
  a = -W( 78 ) / JVS( 479 )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 480 )
  W( 81 ) = W( 81 ) + a*JVS( 481 )
  W( 82 ) = W( 82 ) + a*JVS( 482 )
  W( 83 ) = W( 83 ) + a*JVS( 483 )
  W( 84 ) = W( 84 ) + a*JVS( 484 )
  W( 86 ) = W( 86 ) + a*JVS( 485 )
  W( 87 ) = W( 87 ) + a*JVS( 486 )
  W( 89 ) = W( 89 ) + a*JVS( 487 )
  a = -W( 79 ) / JVS( 518 )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 519 )
  W( 81 ) = W( 81 ) + a*JVS( 520 )
  W( 82 ) = W( 82 ) + a*JVS( 521 )
  W( 83 ) = W( 83 ) + a*JVS( 522 )
  W( 84 ) = W( 84 ) + a*JVS( 523 )
  W( 85 ) = W( 85 ) + a*JVS( 524 )
  W( 86 ) = W( 86 ) + a*JVS( 525 )
  W( 87 ) = W( 87 ) + a*JVS( 526 )
  W( 88 ) = W( 88 ) + a*JVS( 527 )
  W( 89 ) = W( 89 ) + a*JVS( 528 )
  a = -W( 80 ) / JVS( 531 )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 532 )
  W( 82 ) = W( 82 ) + a*JVS( 533 )
  W( 83 ) = W( 83 ) + a*JVS( 534 )
  W( 84 ) = W( 84 ) + a*JVS( 535 )
  W( 86 ) = W( 86 ) + a*JVS( 536 )
  W( 87 ) = W( 87 ) + a*JVS( 537 )
  W( 89 ) = W( 89 ) + a*JVS( 538 )
  a = -W( 81 ) / JVS( 557 )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 558 )
  W( 83 ) = W( 83 ) + a*JVS( 559 )
  W( 84 ) = W( 84 ) + a*JVS( 560 )
  W( 85 ) = W( 85 ) + a*JVS( 561 )
  W( 86 ) = W( 86 ) + a*JVS( 562 )
  W( 87 ) = W( 87 ) + a*JVS( 563 )
  W( 88 ) = W( 88 ) + a*JVS( 564 )
  W( 89 ) = W( 89 ) + a*JVS( 565 )
  a = -W( 82 ) / JVS( 587 )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 588 )
  W( 84 ) = W( 84 ) + a*JVS( 589 )
  W( 85 ) = W( 85 ) + a*JVS( 590 )
  W( 86 ) = W( 86 ) + a*JVS( 591 )
  W( 87 ) = W( 87 ) + a*JVS( 592 )
  W( 88 ) = W( 88 ) + a*JVS( 593 )
  W( 89 ) = W( 89 ) + a*JVS( 594 )
  a = -W( 83 ) / JVS( 606 )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 607 )
  W( 85 ) = W( 85 ) + a*JVS( 608 )
  W( 86 ) = W( 86 ) + a*JVS( 609 )
  W( 87 ) = W( 87 ) + a*JVS( 610 )
  W( 88 ) = W( 88 ) + a*JVS( 611 )
  W( 89 ) = W( 89 ) + a*JVS( 612 )
  a = -W( 84 ) / JVS( 680 )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 681 )
  W( 86 ) = W( 86 ) + a*JVS( 682 )
  W( 87 ) = W( 87 ) + a*JVS( 683 )
  W( 88 ) = W( 88 ) + a*JVS( 684 )
  W( 89 ) = W( 89 ) + a*JVS( 685 )
  a = -W( 85 ) / JVS( 696 )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 697 )
  W( 87 ) = W( 87 ) + a*JVS( 698 )
  W( 88 ) = W( 88 ) + a*JVS( 699 )
  W( 89 ) = W( 89 ) + a*JVS( 700 )
  JVS( 701) = W( 17 )
  JVS( 702) = W( 19 )
  JVS( 703) = W( 20 )
  JVS( 704) = W( 22 )
  JVS( 705) = W( 24 )
  JVS( 706) = W( 25 )
  JVS( 707) = W( 26 )
  JVS( 708) = W( 27 )
  JVS( 709) = W( 28 )
  JVS( 710) = W( 29 )
  JVS( 711) = W( 31 )
  JVS( 712) = W( 33 )
  JVS( 713) = W( 34 )
  JVS( 714) = W( 36 )
  JVS( 715) = W( 39 )
  JVS( 716) = W( 43 )
  JVS( 717) = W( 44 )
  JVS( 718) = W( 46 )
  JVS( 719) = W( 48 )
  JVS( 720) = W( 49 )
  JVS( 721) = W( 50 )
  JVS( 722) = W( 51 )
  JVS( 723) = W( 52 )
  JVS( 724) = W( 53 )
  JVS( 725) = W( 54 )
  JVS( 726) = W( 55 )
  JVS( 727) = W( 56 )
  JVS( 728) = W( 57 )
  JVS( 729) = W( 59 )
  JVS( 730) = W( 60 )
  JVS( 731) = W( 61 )
  JVS( 732) = W( 62 )
  JVS( 733) = W( 63 )
  JVS( 734) = W( 64 )
  JVS( 735) = W( 65 )
  JVS( 736) = W( 66 )
  JVS( 737) = W( 67 )
  JVS( 738) = W( 68 )
  JVS( 739) = W( 69 )
  JVS( 740) = W( 70 )
  JVS( 741) = W( 71 )
  JVS( 742) = W( 72 )
  JVS( 743) = W( 73 )
  JVS( 744) = W( 74 )
  JVS( 745) = W( 75 )
  JVS( 746) = W( 76 )
  JVS( 747) = W( 77 )
  JVS( 748) = W( 78 )
  JVS( 749) = W( 79 )
  JVS( 750) = W( 80 )
  JVS( 751) = W( 81 )
  JVS( 752) = W( 82 )
  JVS( 753) = W( 83 )
  JVS( 754) = W( 84 )
  JVS( 755) = W( 85 )
  JVS( 756) = W( 86 )
  JVS( 757) = W( 87 )
  JVS( 758) = W( 88 )
  JVS( 759) = W( 89 )
  IF ( ABS( JVS( 787 )) < TINY(a) ) THEN
         IER = 87
         RETURN
  END IF
   W( 15 ) = JVS( 760 )
   W( 31 ) = JVS( 761 )
   W( 39 ) = JVS( 762 )
   W( 52 ) = JVS( 763 )
   W( 55 ) = JVS( 764 )
   W( 57 ) = JVS( 765 )
   W( 59 ) = JVS( 766 )
   W( 60 ) = JVS( 767 )
   W( 61 ) = JVS( 768 )
   W( 63 ) = JVS( 769 )
   W( 64 ) = JVS( 770 )
   W( 66 ) = JVS( 771 )
   W( 67 ) = JVS( 772 )
   W( 69 ) = JVS( 773 )
   W( 71 ) = JVS( 774 )
   W( 72 ) = JVS( 775 )
   W( 76 ) = JVS( 776 )
   W( 77 ) = JVS( 777 )
   W( 78 ) = JVS( 778 )
   W( 79 ) = JVS( 779 )
   W( 80 ) = JVS( 780 )
   W( 81 ) = JVS( 781 )
   W( 82 ) = JVS( 782 )
   W( 83 ) = JVS( 783 )
   W( 84 ) = JVS( 784 )
   W( 85 ) = JVS( 785 )
   W( 86 ) = JVS( 786 )
   W( 87 ) = JVS( 787 )
   W( 88 ) = JVS( 788 )
   W( 89 ) = JVS( 789 )
  a = -W( 15 ) / JVS( 36 )
  W( 15 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 37 )
  a = -W( 31 ) / JVS( 93 )
  W( 31 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 94 )
  W( 87 ) = W( 87 ) + a*JVS( 95 )
  a = -W( 39 ) / JVS( 130 )
  W( 39 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 131 )
  W( 84 ) = W( 84 ) + a*JVS( 132 )
  W( 87 ) = W( 87 ) + a*JVS( 133 )
  a = -W( 52 ) / JVS( 206 )
  W( 52 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 207 )
  W( 86 ) = W( 86 ) + a*JVS( 208 )
  W( 87 ) = W( 87 ) + a*JVS( 209 )
  a = -W( 55 ) / JVS( 225 )
  W( 55 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 226 )
  W( 79 ) = W( 79 ) + a*JVS( 227 )
  W( 83 ) = W( 83 ) + a*JVS( 228 )
  W( 84 ) = W( 84 ) + a*JVS( 229 )
  a = -W( 57 ) / JVS( 237 )
  W( 57 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 238 )
  W( 84 ) = W( 84 ) + a*JVS( 239 )
  W( 86 ) = W( 86 ) + a*JVS( 240 )
  W( 87 ) = W( 87 ) + a*JVS( 241 )
  a = -W( 59 ) / JVS( 255 )
  W( 59 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 256 )
  W( 84 ) = W( 84 ) + a*JVS( 257 )
  W( 86 ) = W( 86 ) + a*JVS( 258 )
  W( 87 ) = W( 87 ) + a*JVS( 259 )
  a = -W( 60 ) / JVS( 261 )
  W( 60 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 262 )
  W( 79 ) = W( 79 ) + a*JVS( 263 )
  W( 83 ) = W( 83 ) + a*JVS( 264 )
  W( 84 ) = W( 84 ) + a*JVS( 265 )
  W( 86 ) = W( 86 ) + a*JVS( 266 )
  W( 88 ) = W( 88 ) + a*JVS( 267 )
  a = -W( 61 ) / JVS( 268 )
  W( 61 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 269 )
  W( 81 ) = W( 81 ) + a*JVS( 270 )
  W( 86 ) = W( 86 ) + a*JVS( 271 )
  W( 87 ) = W( 87 ) + a*JVS( 272 )
  a = -W( 63 ) / JVS( 286 )
  W( 63 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 287 )
  W( 86 ) = W( 86 ) + a*JVS( 288 )
  W( 87 ) = W( 87 ) + a*JVS( 289 )
  a = -W( 64 ) / JVS( 295 )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 296 )
  W( 84 ) = W( 84 ) + a*JVS( 297 )
  W( 86 ) = W( 86 ) + a*JVS( 298 )
  W( 87 ) = W( 87 ) + a*JVS( 299 )
  a = -W( 66 ) / JVS( 326 )
  W( 66 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 327 )
  W( 84 ) = W( 84 ) + a*JVS( 328 )
  W( 86 ) = W( 86 ) + a*JVS( 329 )
  W( 87 ) = W( 87 ) + a*JVS( 330 )
  a = -W( 67 ) / JVS( 331 )
  W( 67 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 332 )
  W( 83 ) = W( 83 ) + a*JVS( 333 )
  W( 84 ) = W( 84 ) + a*JVS( 334 )
  a = -W( 69 ) / JVS( 346 )
  W( 69 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  W( 83 ) = W( 83 ) + a*JVS( 349 )
  W( 84 ) = W( 84 ) + a*JVS( 350 )
  a = -W( 71 ) / JVS( 373 )
  W( 71 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 86 ) = W( 86 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  a = -W( 72 ) / JVS( 384 )
  W( 72 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 385 )
  W( 82 ) = W( 82 ) + a*JVS( 386 )
  W( 83 ) = W( 83 ) + a*JVS( 387 )
  W( 84 ) = W( 84 ) + a*JVS( 388 )
  W( 86 ) = W( 86 ) + a*JVS( 389 )
  W( 87 ) = W( 87 ) + a*JVS( 390 )
  a = -W( 76 ) / JVS( 451 )
  W( 76 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 452 )
  W( 80 ) = W( 80 ) + a*JVS( 453 )
  W( 81 ) = W( 81 ) + a*JVS( 454 )
  W( 82 ) = W( 82 ) + a*JVS( 455 )
  W( 83 ) = W( 83 ) + a*JVS( 456 )
  W( 84 ) = W( 84 ) + a*JVS( 457 )
  W( 86 ) = W( 86 ) + a*JVS( 458 )
  W( 87 ) = W( 87 ) + a*JVS( 459 )
  W( 89 ) = W( 89 ) + a*JVS( 460 )
  a = -W( 77 ) / JVS( 467 )
  W( 77 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 468 )
  W( 81 ) = W( 81 ) + a*JVS( 469 )
  W( 82 ) = W( 82 ) + a*JVS( 470 )
  W( 83 ) = W( 83 ) + a*JVS( 471 )
  W( 84 ) = W( 84 ) + a*JVS( 472 )
  W( 86 ) = W( 86 ) + a*JVS( 473 )
  W( 87 ) = W( 87 ) + a*JVS( 474 )
  W( 89 ) = W( 89 ) + a*JVS( 475 )
  a = -W( 78 ) / JVS( 479 )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 480 )
  W( 81 ) = W( 81 ) + a*JVS( 481 )
  W( 82 ) = W( 82 ) + a*JVS( 482 )
  W( 83 ) = W( 83 ) + a*JVS( 483 )
  W( 84 ) = W( 84 ) + a*JVS( 484 )
  W( 86 ) = W( 86 ) + a*JVS( 485 )
  W( 87 ) = W( 87 ) + a*JVS( 486 )
  W( 89 ) = W( 89 ) + a*JVS( 487 )
  a = -W( 79 ) / JVS( 518 )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 519 )
  W( 81 ) = W( 81 ) + a*JVS( 520 )
  W( 82 ) = W( 82 ) + a*JVS( 521 )
  W( 83 ) = W( 83 ) + a*JVS( 522 )
  W( 84 ) = W( 84 ) + a*JVS( 523 )
  W( 85 ) = W( 85 ) + a*JVS( 524 )
  W( 86 ) = W( 86 ) + a*JVS( 525 )
  W( 87 ) = W( 87 ) + a*JVS( 526 )
  W( 88 ) = W( 88 ) + a*JVS( 527 )
  W( 89 ) = W( 89 ) + a*JVS( 528 )
  a = -W( 80 ) / JVS( 531 )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 532 )
  W( 82 ) = W( 82 ) + a*JVS( 533 )
  W( 83 ) = W( 83 ) + a*JVS( 534 )
  W( 84 ) = W( 84 ) + a*JVS( 535 )
  W( 86 ) = W( 86 ) + a*JVS( 536 )
  W( 87 ) = W( 87 ) + a*JVS( 537 )
  W( 89 ) = W( 89 ) + a*JVS( 538 )
  a = -W( 81 ) / JVS( 557 )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 558 )
  W( 83 ) = W( 83 ) + a*JVS( 559 )
  W( 84 ) = W( 84 ) + a*JVS( 560 )
  W( 85 ) = W( 85 ) + a*JVS( 561 )
  W( 86 ) = W( 86 ) + a*JVS( 562 )
  W( 87 ) = W( 87 ) + a*JVS( 563 )
  W( 88 ) = W( 88 ) + a*JVS( 564 )
  W( 89 ) = W( 89 ) + a*JVS( 565 )
  a = -W( 82 ) / JVS( 587 )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 588 )
  W( 84 ) = W( 84 ) + a*JVS( 589 )
  W( 85 ) = W( 85 ) + a*JVS( 590 )
  W( 86 ) = W( 86 ) + a*JVS( 591 )
  W( 87 ) = W( 87 ) + a*JVS( 592 )
  W( 88 ) = W( 88 ) + a*JVS( 593 )
  W( 89 ) = W( 89 ) + a*JVS( 594 )
  a = -W( 83 ) / JVS( 606 )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 607 )
  W( 85 ) = W( 85 ) + a*JVS( 608 )
  W( 86 ) = W( 86 ) + a*JVS( 609 )
  W( 87 ) = W( 87 ) + a*JVS( 610 )
  W( 88 ) = W( 88 ) + a*JVS( 611 )
  W( 89 ) = W( 89 ) + a*JVS( 612 )
  a = -W( 84 ) / JVS( 680 )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 681 )
  W( 86 ) = W( 86 ) + a*JVS( 682 )
  W( 87 ) = W( 87 ) + a*JVS( 683 )
  W( 88 ) = W( 88 ) + a*JVS( 684 )
  W( 89 ) = W( 89 ) + a*JVS( 685 )
  a = -W( 85 ) / JVS( 696 )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 697 )
  W( 87 ) = W( 87 ) + a*JVS( 698 )
  W( 88 ) = W( 88 ) + a*JVS( 699 )
  W( 89 ) = W( 89 ) + a*JVS( 700 )
  a = -W( 86 ) / JVS( 756 )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 757 )
  W( 88 ) = W( 88 ) + a*JVS( 758 )
  W( 89 ) = W( 89 ) + a*JVS( 759 )
  JVS( 760) = W( 15 )
  JVS( 761) = W( 31 )
  JVS( 762) = W( 39 )
  JVS( 763) = W( 52 )
  JVS( 764) = W( 55 )
  JVS( 765) = W( 57 )
  JVS( 766) = W( 59 )
  JVS( 767) = W( 60 )
  JVS( 768) = W( 61 )
  JVS( 769) = W( 63 )
  JVS( 770) = W( 64 )
  JVS( 771) = W( 66 )
  JVS( 772) = W( 67 )
  JVS( 773) = W( 69 )
  JVS( 774) = W( 71 )
  JVS( 775) = W( 72 )
  JVS( 776) = W( 76 )
  JVS( 777) = W( 77 )
  JVS( 778) = W( 78 )
  JVS( 779) = W( 79 )
  JVS( 780) = W( 80 )
  JVS( 781) = W( 81 )
  JVS( 782) = W( 82 )
  JVS( 783) = W( 83 )
  JVS( 784) = W( 84 )
  JVS( 785) = W( 85 )
  JVS( 786) = W( 86 )
  JVS( 787) = W( 87 )
  JVS( 788) = W( 88 )
  JVS( 789) = W( 89 )
  IF ( ABS( JVS( 827 )) < TINY(a) ) THEN
         IER = 88
         RETURN
  END IF
   W( 21 ) = JVS( 790 )
   W( 22 ) = JVS( 791 )
   W( 31 ) = JVS( 792 )
   W( 34 ) = JVS( 793 )
   W( 37 ) = JVS( 794 )
   W( 39 ) = JVS( 795 )
   W( 45 ) = JVS( 796 )
   W( 47 ) = JVS( 797 )
   W( 50 ) = JVS( 798 )
   W( 52 ) = JVS( 799 )
   W( 54 ) = JVS( 800 )
   W( 57 ) = JVS( 801 )
   W( 59 ) = JVS( 802 )
   W( 60 ) = JVS( 803 )
   W( 61 ) = JVS( 804 )
   W( 63 ) = JVS( 805 )
   W( 64 ) = JVS( 806 )
   W( 66 ) = JVS( 807 )
   W( 67 ) = JVS( 808 )
   W( 68 ) = JVS( 809 )
   W( 69 ) = JVS( 810 )
   W( 70 ) = JVS( 811 )
   W( 71 ) = JVS( 812 )
   W( 72 ) = JVS( 813 )
   W( 75 ) = JVS( 814 )
   W( 76 ) = JVS( 815 )
   W( 77 ) = JVS( 816 )
   W( 78 ) = JVS( 817 )
   W( 79 ) = JVS( 818 )
   W( 80 ) = JVS( 819 )
   W( 81 ) = JVS( 820 )
   W( 82 ) = JVS( 821 )
   W( 83 ) = JVS( 822 )
   W( 84 ) = JVS( 823 )
   W( 85 ) = JVS( 824 )
   W( 86 ) = JVS( 825 )
   W( 87 ) = JVS( 826 )
   W( 88 ) = JVS( 827 )
   W( 89 ) = JVS( 828 )
  a = -W( 21 ) / JVS( 53 )
  W( 21 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 54 )
  W( 88 ) = W( 88 ) + a*JVS( 55 )
  a = -W( 22 ) / JVS( 57 )
  W( 22 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 58 )
  W( 88 ) = W( 88 ) + a*JVS( 59 )
  a = -W( 31 ) / JVS( 93 )
  W( 31 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 94 )
  W( 87 ) = W( 87 ) + a*JVS( 95 )
  a = -W( 34 ) / JVS( 104 )
  W( 34 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 105 )
  W( 86 ) = W( 86 ) + a*JVS( 106 )
  W( 88 ) = W( 88 ) + a*JVS( 107 )
  a = -W( 37 ) / JVS( 119 )
  W( 37 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 120 )
  W( 69 ) = W( 69 ) + a*JVS( 121 )
  W( 81 ) = W( 81 ) + a*JVS( 122 )
  W( 84 ) = W( 84 ) + a*JVS( 123 )
  W( 87 ) = W( 87 ) + a*JVS( 124 )
  a = -W( 39 ) / JVS( 130 )
  W( 39 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 131 )
  W( 84 ) = W( 84 ) + a*JVS( 132 )
  W( 87 ) = W( 87 ) + a*JVS( 133 )
  a = -W( 45 ) / JVS( 168 )
  W( 45 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 169 )
  W( 75 ) = W( 75 ) + a*JVS( 170 )
  W( 79 ) = W( 79 ) + a*JVS( 171 )
  W( 81 ) = W( 81 ) + a*JVS( 172 )
  W( 84 ) = W( 84 ) + a*JVS( 173 )
  W( 88 ) = W( 88 ) + a*JVS( 174 )
  a = -W( 47 ) / JVS( 179 )
  W( 47 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 180 )
  W( 88 ) = W( 88 ) + a*JVS( 181 )
  W( 89 ) = W( 89 ) + a*JVS( 182 )
  a = -W( 50 ) / JVS( 192 )
  W( 50 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 193 )
  W( 85 ) = W( 85 ) + a*JVS( 194 )
  W( 88 ) = W( 88 ) + a*JVS( 195 )
  a = -W( 52 ) / JVS( 206 )
  W( 52 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 207 )
  W( 86 ) = W( 86 ) + a*JVS( 208 )
  W( 87 ) = W( 87 ) + a*JVS( 209 )
  a = -W( 54 ) / JVS( 218 )
  W( 54 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 219 )
  W( 83 ) = W( 83 ) + a*JVS( 220 )
  W( 84 ) = W( 84 ) + a*JVS( 221 )
  a = -W( 57 ) / JVS( 237 )
  W( 57 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 238 )
  W( 84 ) = W( 84 ) + a*JVS( 239 )
  W( 86 ) = W( 86 ) + a*JVS( 240 )
  W( 87 ) = W( 87 ) + a*JVS( 241 )
  a = -W( 59 ) / JVS( 255 )
  W( 59 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 256 )
  W( 84 ) = W( 84 ) + a*JVS( 257 )
  W( 86 ) = W( 86 ) + a*JVS( 258 )
  W( 87 ) = W( 87 ) + a*JVS( 259 )
  a = -W( 60 ) / JVS( 261 )
  W( 60 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 262 )
  W( 79 ) = W( 79 ) + a*JVS( 263 )
  W( 83 ) = W( 83 ) + a*JVS( 264 )
  W( 84 ) = W( 84 ) + a*JVS( 265 )
  W( 86 ) = W( 86 ) + a*JVS( 266 )
  W( 88 ) = W( 88 ) + a*JVS( 267 )
  a = -W( 61 ) / JVS( 268 )
  W( 61 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 269 )
  W( 81 ) = W( 81 ) + a*JVS( 270 )
  W( 86 ) = W( 86 ) + a*JVS( 271 )
  W( 87 ) = W( 87 ) + a*JVS( 272 )
  a = -W( 63 ) / JVS( 286 )
  W( 63 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 287 )
  W( 86 ) = W( 86 ) + a*JVS( 288 )
  W( 87 ) = W( 87 ) + a*JVS( 289 )
  a = -W( 64 ) / JVS( 295 )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 296 )
  W( 84 ) = W( 84 ) + a*JVS( 297 )
  W( 86 ) = W( 86 ) + a*JVS( 298 )
  W( 87 ) = W( 87 ) + a*JVS( 299 )
  a = -W( 66 ) / JVS( 326 )
  W( 66 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 327 )
  W( 84 ) = W( 84 ) + a*JVS( 328 )
  W( 86 ) = W( 86 ) + a*JVS( 329 )
  W( 87 ) = W( 87 ) + a*JVS( 330 )
  a = -W( 67 ) / JVS( 331 )
  W( 67 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 332 )
  W( 83 ) = W( 83 ) + a*JVS( 333 )
  W( 84 ) = W( 84 ) + a*JVS( 334 )
  a = -W( 68 ) / JVS( 337 )
  W( 68 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 338 )
  W( 80 ) = W( 80 ) + a*JVS( 339 )
  W( 81 ) = W( 81 ) + a*JVS( 340 )
  W( 83 ) = W( 83 ) + a*JVS( 341 )
  W( 84 ) = W( 84 ) + a*JVS( 342 )
  W( 86 ) = W( 86 ) + a*JVS( 343 )
  W( 87 ) = W( 87 ) + a*JVS( 344 )
  a = -W( 69 ) / JVS( 346 )
  W( 69 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  W( 83 ) = W( 83 ) + a*JVS( 349 )
  W( 84 ) = W( 84 ) + a*JVS( 350 )
  a = -W( 70 ) / JVS( 363 )
  W( 70 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 364 )
  W( 81 ) = W( 81 ) + a*JVS( 365 )
  W( 82 ) = W( 82 ) + a*JVS( 366 )
  W( 83 ) = W( 83 ) + a*JVS( 367 )
  W( 84 ) = W( 84 ) + a*JVS( 368 )
  W( 86 ) = W( 86 ) + a*JVS( 369 )
  W( 87 ) = W( 87 ) + a*JVS( 370 )
  a = -W( 71 ) / JVS( 373 )
  W( 71 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 86 ) = W( 86 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  a = -W( 72 ) / JVS( 384 )
  W( 72 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 385 )
  W( 82 ) = W( 82 ) + a*JVS( 386 )
  W( 83 ) = W( 83 ) + a*JVS( 387 )
  W( 84 ) = W( 84 ) + a*JVS( 388 )
  W( 86 ) = W( 86 ) + a*JVS( 389 )
  W( 87 ) = W( 87 ) + a*JVS( 390 )
  a = -W( 75 ) / JVS( 433 )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 434 )
  W( 77 ) = W( 77 ) + a*JVS( 435 )
  W( 78 ) = W( 78 ) + a*JVS( 436 )
  W( 80 ) = W( 80 ) + a*JVS( 437 )
  W( 81 ) = W( 81 ) + a*JVS( 438 )
  W( 82 ) = W( 82 ) + a*JVS( 439 )
  W( 83 ) = W( 83 ) + a*JVS( 440 )
  W( 84 ) = W( 84 ) + a*JVS( 441 )
  W( 85 ) = W( 85 ) + a*JVS( 442 )
  W( 86 ) = W( 86 ) + a*JVS( 443 )
  W( 87 ) = W( 87 ) + a*JVS( 444 )
  W( 88 ) = W( 88 ) + a*JVS( 445 )
  W( 89 ) = W( 89 ) + a*JVS( 446 )
  a = -W( 76 ) / JVS( 451 )
  W( 76 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 452 )
  W( 80 ) = W( 80 ) + a*JVS( 453 )
  W( 81 ) = W( 81 ) + a*JVS( 454 )
  W( 82 ) = W( 82 ) + a*JVS( 455 )
  W( 83 ) = W( 83 ) + a*JVS( 456 )
  W( 84 ) = W( 84 ) + a*JVS( 457 )
  W( 86 ) = W( 86 ) + a*JVS( 458 )
  W( 87 ) = W( 87 ) + a*JVS( 459 )
  W( 89 ) = W( 89 ) + a*JVS( 460 )
  a = -W( 77 ) / JVS( 467 )
  W( 77 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 468 )
  W( 81 ) = W( 81 ) + a*JVS( 469 )
  W( 82 ) = W( 82 ) + a*JVS( 470 )
  W( 83 ) = W( 83 ) + a*JVS( 471 )
  W( 84 ) = W( 84 ) + a*JVS( 472 )
  W( 86 ) = W( 86 ) + a*JVS( 473 )
  W( 87 ) = W( 87 ) + a*JVS( 474 )
  W( 89 ) = W( 89 ) + a*JVS( 475 )
  a = -W( 78 ) / JVS( 479 )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 480 )
  W( 81 ) = W( 81 ) + a*JVS( 481 )
  W( 82 ) = W( 82 ) + a*JVS( 482 )
  W( 83 ) = W( 83 ) + a*JVS( 483 )
  W( 84 ) = W( 84 ) + a*JVS( 484 )
  W( 86 ) = W( 86 ) + a*JVS( 485 )
  W( 87 ) = W( 87 ) + a*JVS( 486 )
  W( 89 ) = W( 89 ) + a*JVS( 487 )
  a = -W( 79 ) / JVS( 518 )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 519 )
  W( 81 ) = W( 81 ) + a*JVS( 520 )
  W( 82 ) = W( 82 ) + a*JVS( 521 )
  W( 83 ) = W( 83 ) + a*JVS( 522 )
  W( 84 ) = W( 84 ) + a*JVS( 523 )
  W( 85 ) = W( 85 ) + a*JVS( 524 )
  W( 86 ) = W( 86 ) + a*JVS( 525 )
  W( 87 ) = W( 87 ) + a*JVS( 526 )
  W( 88 ) = W( 88 ) + a*JVS( 527 )
  W( 89 ) = W( 89 ) + a*JVS( 528 )
  a = -W( 80 ) / JVS( 531 )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 532 )
  W( 82 ) = W( 82 ) + a*JVS( 533 )
  W( 83 ) = W( 83 ) + a*JVS( 534 )
  W( 84 ) = W( 84 ) + a*JVS( 535 )
  W( 86 ) = W( 86 ) + a*JVS( 536 )
  W( 87 ) = W( 87 ) + a*JVS( 537 )
  W( 89 ) = W( 89 ) + a*JVS( 538 )
  a = -W( 81 ) / JVS( 557 )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 558 )
  W( 83 ) = W( 83 ) + a*JVS( 559 )
  W( 84 ) = W( 84 ) + a*JVS( 560 )
  W( 85 ) = W( 85 ) + a*JVS( 561 )
  W( 86 ) = W( 86 ) + a*JVS( 562 )
  W( 87 ) = W( 87 ) + a*JVS( 563 )
  W( 88 ) = W( 88 ) + a*JVS( 564 )
  W( 89 ) = W( 89 ) + a*JVS( 565 )
  a = -W( 82 ) / JVS( 587 )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 588 )
  W( 84 ) = W( 84 ) + a*JVS( 589 )
  W( 85 ) = W( 85 ) + a*JVS( 590 )
  W( 86 ) = W( 86 ) + a*JVS( 591 )
  W( 87 ) = W( 87 ) + a*JVS( 592 )
  W( 88 ) = W( 88 ) + a*JVS( 593 )
  W( 89 ) = W( 89 ) + a*JVS( 594 )
  a = -W( 83 ) / JVS( 606 )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 607 )
  W( 85 ) = W( 85 ) + a*JVS( 608 )
  W( 86 ) = W( 86 ) + a*JVS( 609 )
  W( 87 ) = W( 87 ) + a*JVS( 610 )
  W( 88 ) = W( 88 ) + a*JVS( 611 )
  W( 89 ) = W( 89 ) + a*JVS( 612 )
  a = -W( 84 ) / JVS( 680 )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 681 )
  W( 86 ) = W( 86 ) + a*JVS( 682 )
  W( 87 ) = W( 87 ) + a*JVS( 683 )
  W( 88 ) = W( 88 ) + a*JVS( 684 )
  W( 89 ) = W( 89 ) + a*JVS( 685 )
  a = -W( 85 ) / JVS( 696 )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 697 )
  W( 87 ) = W( 87 ) + a*JVS( 698 )
  W( 88 ) = W( 88 ) + a*JVS( 699 )
  W( 89 ) = W( 89 ) + a*JVS( 700 )
  a = -W( 86 ) / JVS( 756 )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 757 )
  W( 88 ) = W( 88 ) + a*JVS( 758 )
  W( 89 ) = W( 89 ) + a*JVS( 759 )
  a = -W( 87 ) / JVS( 787 )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 788 )
  W( 89 ) = W( 89 ) + a*JVS( 789 )
  JVS( 790) = W( 21 )
  JVS( 791) = W( 22 )
  JVS( 792) = W( 31 )
  JVS( 793) = W( 34 )
  JVS( 794) = W( 37 )
  JVS( 795) = W( 39 )
  JVS( 796) = W( 45 )
  JVS( 797) = W( 47 )
  JVS( 798) = W( 50 )
  JVS( 799) = W( 52 )
  JVS( 800) = W( 54 )
  JVS( 801) = W( 57 )
  JVS( 802) = W( 59 )
  JVS( 803) = W( 60 )
  JVS( 804) = W( 61 )
  JVS( 805) = W( 63 )
  JVS( 806) = W( 64 )
  JVS( 807) = W( 66 )
  JVS( 808) = W( 67 )
  JVS( 809) = W( 68 )
  JVS( 810) = W( 69 )
  JVS( 811) = W( 70 )
  JVS( 812) = W( 71 )
  JVS( 813) = W( 72 )
  JVS( 814) = W( 75 )
  JVS( 815) = W( 76 )
  JVS( 816) = W( 77 )
  JVS( 817) = W( 78 )
  JVS( 818) = W( 79 )
  JVS( 819) = W( 80 )
  JVS( 820) = W( 81 )
  JVS( 821) = W( 82 )
  JVS( 822) = W( 83 )
  JVS( 823) = W( 84 )
  JVS( 824) = W( 85 )
  JVS( 825) = W( 86 )
  JVS( 826) = W( 87 )
  JVS( 827) = W( 88 )
  JVS( 828) = W( 89 )
  IF ( ABS( JVS( 858 )) < TINY(a) ) THEN
         IER = 89
         RETURN
  END IF
   W( 30 ) = JVS( 829 )
   W( 32 ) = JVS( 830 )
   W( 40 ) = JVS( 831 )
   W( 41 ) = JVS( 832 )
   W( 47 ) = JVS( 833 )
   W( 53 ) = JVS( 834 )
   W( 56 ) = JVS( 835 )
   W( 57 ) = JVS( 836 )
   W( 58 ) = JVS( 837 )
   W( 63 ) = JVS( 838 )
   W( 66 ) = JVS( 839 )
   W( 70 ) = JVS( 840 )
   W( 71 ) = JVS( 841 )
   W( 72 ) = JVS( 842 )
   W( 73 ) = JVS( 843 )
   W( 74 ) = JVS( 844 )
   W( 75 ) = JVS( 845 )
   W( 76 ) = JVS( 846 )
   W( 77 ) = JVS( 847 )
   W( 78 ) = JVS( 848 )
   W( 80 ) = JVS( 849 )
   W( 81 ) = JVS( 850 )
   W( 82 ) = JVS( 851 )
   W( 83 ) = JVS( 852 )
   W( 84 ) = JVS( 853 )
   W( 85 ) = JVS( 854 )
   W( 86 ) = JVS( 855 )
   W( 87 ) = JVS( 856 )
   W( 88 ) = JVS( 857 )
   W( 89 ) = JVS( 858 )
  a = -W( 30 ) / JVS( 88 )
  W( 30 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 89 )
  W( 84 ) = W( 84 ) + a*JVS( 90 )
  W( 86 ) = W( 86 ) + a*JVS( 91 )
  a = -W( 32 ) / JVS( 96 )
  W( 32 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 97 )
  W( 84 ) = W( 84 ) + a*JVS( 98 )
  W( 86 ) = W( 86 ) + a*JVS( 99 )
  a = -W( 40 ) / JVS( 134 )
  W( 40 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 135 )
  W( 85 ) = W( 85 ) + a*JVS( 136 )
  W( 86 ) = W( 86 ) + a*JVS( 137 )
  W( 89 ) = W( 89 ) + a*JVS( 138 )
  a = -W( 41 ) / JVS( 139 )
  W( 41 ) = -a
  W( 56 ) = W( 56 ) + a*JVS( 140 )
  W( 63 ) = W( 63 ) + a*JVS( 141 )
  W( 84 ) = W( 84 ) + a*JVS( 142 )
  W( 87 ) = W( 87 ) + a*JVS( 143 )
  a = -W( 47 ) / JVS( 179 )
  W( 47 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 180 )
  W( 88 ) = W( 88 ) + a*JVS( 181 )
  W( 89 ) = W( 89 ) + a*JVS( 182 )
  a = -W( 53 ) / JVS( 213 )
  W( 53 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 214 )
  W( 86 ) = W( 86 ) + a*JVS( 215 )
  W( 87 ) = W( 87 ) + a*JVS( 216 )
  W( 88 ) = W( 88 ) + a*JVS( 217 )
  a = -W( 56 ) / JVS( 230 )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 231 )
  W( 84 ) = W( 84 ) + a*JVS( 232 )
  W( 86 ) = W( 86 ) + a*JVS( 233 )
  a = -W( 57 ) / JVS( 237 )
  W( 57 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 238 )
  W( 84 ) = W( 84 ) + a*JVS( 239 )
  W( 86 ) = W( 86 ) + a*JVS( 240 )
  W( 87 ) = W( 87 ) + a*JVS( 241 )
  a = -W( 58 ) / JVS( 246 )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 247 )
  W( 66 ) = W( 66 ) + a*JVS( 248 )
  W( 71 ) = W( 71 ) + a*JVS( 249 )
  W( 82 ) = W( 82 ) + a*JVS( 250 )
  W( 84 ) = W( 84 ) + a*JVS( 251 )
  W( 86 ) = W( 86 ) + a*JVS( 252 )
  W( 87 ) = W( 87 ) + a*JVS( 253 )
  a = -W( 63 ) / JVS( 286 )
  W( 63 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 287 )
  W( 86 ) = W( 86 ) + a*JVS( 288 )
  W( 87 ) = W( 87 ) + a*JVS( 289 )
  a = -W( 66 ) / JVS( 326 )
  W( 66 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 327 )
  W( 84 ) = W( 84 ) + a*JVS( 328 )
  W( 86 ) = W( 86 ) + a*JVS( 329 )
  W( 87 ) = W( 87 ) + a*JVS( 330 )
  a = -W( 70 ) / JVS( 363 )
  W( 70 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 364 )
  W( 81 ) = W( 81 ) + a*JVS( 365 )
  W( 82 ) = W( 82 ) + a*JVS( 366 )
  W( 83 ) = W( 83 ) + a*JVS( 367 )
  W( 84 ) = W( 84 ) + a*JVS( 368 )
  W( 86 ) = W( 86 ) + a*JVS( 369 )
  W( 87 ) = W( 87 ) + a*JVS( 370 )
  a = -W( 71 ) / JVS( 373 )
  W( 71 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 86 ) = W( 86 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  a = -W( 72 ) / JVS( 384 )
  W( 72 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 385 )
  W( 82 ) = W( 82 ) + a*JVS( 386 )
  W( 83 ) = W( 83 ) + a*JVS( 387 )
  W( 84 ) = W( 84 ) + a*JVS( 388 )
  W( 86 ) = W( 86 ) + a*JVS( 389 )
  W( 87 ) = W( 87 ) + a*JVS( 390 )
  a = -W( 73 ) / JVS( 396 )
  W( 73 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 397 )
  W( 77 ) = W( 77 ) + a*JVS( 398 )
  W( 78 ) = W( 78 ) + a*JVS( 399 )
  W( 81 ) = W( 81 ) + a*JVS( 400 )
  W( 82 ) = W( 82 ) + a*JVS( 401 )
  W( 83 ) = W( 83 ) + a*JVS( 402 )
  W( 84 ) = W( 84 ) + a*JVS( 403 )
  W( 85 ) = W( 85 ) + a*JVS( 404 )
  W( 86 ) = W( 86 ) + a*JVS( 405 )
  W( 87 ) = W( 87 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 89 ) = W( 89 ) + a*JVS( 408 )
  a = -W( 74 ) / JVS( 415 )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 416 )
  W( 81 ) = W( 81 ) + a*JVS( 417 )
  W( 82 ) = W( 82 ) + a*JVS( 418 )
  W( 83 ) = W( 83 ) + a*JVS( 419 )
  W( 84 ) = W( 84 ) + a*JVS( 420 )
  W( 86 ) = W( 86 ) + a*JVS( 421 )
  W( 87 ) = W( 87 ) + a*JVS( 422 )
  W( 89 ) = W( 89 ) + a*JVS( 423 )
  a = -W( 75 ) / JVS( 433 )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 434 )
  W( 77 ) = W( 77 ) + a*JVS( 435 )
  W( 78 ) = W( 78 ) + a*JVS( 436 )
  W( 80 ) = W( 80 ) + a*JVS( 437 )
  W( 81 ) = W( 81 ) + a*JVS( 438 )
  W( 82 ) = W( 82 ) + a*JVS( 439 )
  W( 83 ) = W( 83 ) + a*JVS( 440 )
  W( 84 ) = W( 84 ) + a*JVS( 441 )
  W( 85 ) = W( 85 ) + a*JVS( 442 )
  W( 86 ) = W( 86 ) + a*JVS( 443 )
  W( 87 ) = W( 87 ) + a*JVS( 444 )
  W( 88 ) = W( 88 ) + a*JVS( 445 )
  W( 89 ) = W( 89 ) + a*JVS( 446 )
  a = -W( 76 ) / JVS( 451 )
  W( 76 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 452 )
  W( 80 ) = W( 80 ) + a*JVS( 453 )
  W( 81 ) = W( 81 ) + a*JVS( 454 )
  W( 82 ) = W( 82 ) + a*JVS( 455 )
  W( 83 ) = W( 83 ) + a*JVS( 456 )
  W( 84 ) = W( 84 ) + a*JVS( 457 )
  W( 86 ) = W( 86 ) + a*JVS( 458 )
  W( 87 ) = W( 87 ) + a*JVS( 459 )
  W( 89 ) = W( 89 ) + a*JVS( 460 )
  a = -W( 77 ) / JVS( 467 )
  W( 77 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 468 )
  W( 81 ) = W( 81 ) + a*JVS( 469 )
  W( 82 ) = W( 82 ) + a*JVS( 470 )
  W( 83 ) = W( 83 ) + a*JVS( 471 )
  W( 84 ) = W( 84 ) + a*JVS( 472 )
  W( 86 ) = W( 86 ) + a*JVS( 473 )
  W( 87 ) = W( 87 ) + a*JVS( 474 )
  W( 89 ) = W( 89 ) + a*JVS( 475 )
  a = -W( 78 ) / JVS( 479 )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 480 )
  W( 81 ) = W( 81 ) + a*JVS( 481 )
  W( 82 ) = W( 82 ) + a*JVS( 482 )
  W( 83 ) = W( 83 ) + a*JVS( 483 )
  W( 84 ) = W( 84 ) + a*JVS( 484 )
  W( 86 ) = W( 86 ) + a*JVS( 485 )
  W( 87 ) = W( 87 ) + a*JVS( 486 )
  W( 89 ) = W( 89 ) + a*JVS( 487 )
  a = -W( 80 ) / JVS( 531 )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 532 )
  W( 82 ) = W( 82 ) + a*JVS( 533 )
  W( 83 ) = W( 83 ) + a*JVS( 534 )
  W( 84 ) = W( 84 ) + a*JVS( 535 )
  W( 86 ) = W( 86 ) + a*JVS( 536 )
  W( 87 ) = W( 87 ) + a*JVS( 537 )
  W( 89 ) = W( 89 ) + a*JVS( 538 )
  a = -W( 81 ) / JVS( 557 )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 558 )
  W( 83 ) = W( 83 ) + a*JVS( 559 )
  W( 84 ) = W( 84 ) + a*JVS( 560 )
  W( 85 ) = W( 85 ) + a*JVS( 561 )
  W( 86 ) = W( 86 ) + a*JVS( 562 )
  W( 87 ) = W( 87 ) + a*JVS( 563 )
  W( 88 ) = W( 88 ) + a*JVS( 564 )
  W( 89 ) = W( 89 ) + a*JVS( 565 )
  a = -W( 82 ) / JVS( 587 )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 588 )
  W( 84 ) = W( 84 ) + a*JVS( 589 )
  W( 85 ) = W( 85 ) + a*JVS( 590 )
  W( 86 ) = W( 86 ) + a*JVS( 591 )
  W( 87 ) = W( 87 ) + a*JVS( 592 )
  W( 88 ) = W( 88 ) + a*JVS( 593 )
  W( 89 ) = W( 89 ) + a*JVS( 594 )
  a = -W( 83 ) / JVS( 606 )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 607 )
  W( 85 ) = W( 85 ) + a*JVS( 608 )
  W( 86 ) = W( 86 ) + a*JVS( 609 )
  W( 87 ) = W( 87 ) + a*JVS( 610 )
  W( 88 ) = W( 88 ) + a*JVS( 611 )
  W( 89 ) = W( 89 ) + a*JVS( 612 )
  a = -W( 84 ) / JVS( 680 )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 681 )
  W( 86 ) = W( 86 ) + a*JVS( 682 )
  W( 87 ) = W( 87 ) + a*JVS( 683 )
  W( 88 ) = W( 88 ) + a*JVS( 684 )
  W( 89 ) = W( 89 ) + a*JVS( 685 )
  a = -W( 85 ) / JVS( 696 )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 697 )
  W( 87 ) = W( 87 ) + a*JVS( 698 )
  W( 88 ) = W( 88 ) + a*JVS( 699 )
  W( 89 ) = W( 89 ) + a*JVS( 700 )
  a = -W( 86 ) / JVS( 756 )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 757 )
  W( 88 ) = W( 88 ) + a*JVS( 758 )
  W( 89 ) = W( 89 ) + a*JVS( 759 )
  a = -W( 87 ) / JVS( 787 )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 788 )
  W( 89 ) = W( 89 ) + a*JVS( 789 )
  a = -W( 88 ) / JVS( 827 )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 828 )
  JVS( 829) = W( 30 )
  JVS( 830) = W( 32 )
  JVS( 831) = W( 40 )
  JVS( 832) = W( 41 )
  JVS( 833) = W( 47 )
  JVS( 834) = W( 53 )
  JVS( 835) = W( 56 )
  JVS( 836) = W( 57 )
  JVS( 837) = W( 58 )
  JVS( 838) = W( 63 )
  JVS( 839) = W( 66 )
  JVS( 840) = W( 70 )
  JVS( 841) = W( 71 )
  JVS( 842) = W( 72 )
  JVS( 843) = W( 73 )
  JVS( 844) = W( 74 )
  JVS( 845) = W( 75 )
  JVS( 846) = W( 76 )
  JVS( 847) = W( 77 )
  JVS( 848) = W( 78 )
  JVS( 849) = W( 80 )
  JVS( 850) = W( 81 )
  JVS( 851) = W( 82 )
  JVS( 852) = W( 83 )
  JVS( 853) = W( 84 )
  JVS( 854) = W( 85 )
  JVS( 855) = W( 86 )
  JVS( 856) = W( 87 )
  JVS( 857) = W( 88 )
  JVS( 858) = W( 89 )
   END SUBROUTINE decomp_mozart_mosaic_4bin_vbs0
END MODULE mozart_mosaic_4bin_vbs0_Integrator
