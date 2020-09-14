 MODULE module_lightning_nox_decaria_cup
 IMPLICIT NONE
 CONTAINS
 SUBROUTINE lightning_nox_decaria_cup ( &
                            dx, dy, xland, ht, t_phy, rho, z, p, &
                            ic_flashrate, cg_flashrate, &
                            cldfra_cup, htop, &
                            N_IC, N_CG, &
                            ltng_temp_upper,ltng_temp_lower, &
                            ids, ide, jds, jde, kds, kde, &
                            ims, ime, jms, jme, kms, kme, &
                            ips, ipe, jps, jpe, kps, kpe, &
                            lnox_ic_tend, lnox_cg_tend &
                          )
 USE module_state_description
 USE module_model_constants
 USE module_wrf_error
 USE module_dm, only: wrf_dm_max_real, wrf_dm_min_real, wrf_dm_sum_real
 IMPLICIT NONE
 REAL, INTENT(IN ) :: dx, dy
 REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN ) :: xland, ht, htop
 REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN ) :: t_phy, rho, z, p
 REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN ) :: cldfra_cup
 REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN ) :: ic_flashrate , cg_flashrate
 REAL, INTENT(IN ) :: N_IC, N_CG
 REAL, INTENT(IN ) :: ltng_temp_lower, ltng_temp_upper
 INTEGER, INTENT(IN ) :: ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN ) :: ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN ) :: ips,ipe, jps,jpe, kps,kpe
 REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( OUT) :: lnox_ic_tend,lnox_cg_tend
 INTEGER :: i,k,j
 INTEGER :: ktop,kbtm,kupper,klower
 REAL :: ic_fr, cg_fr, delta
 CHARACTER (LEN=250) :: message
 REAL, DIMENSION( kps:kpe ) :: molesofair, cldfra_cup_1d, t_1d, z_1d
 REAL :: htop_ij
 REAL, DIMENSION( kps:kpe ) :: fd, fd2, dz
 lnox_ic_tend (ips:ipe,kps:kpe,jps:jpe ) = 0.
 lnox_cg_tend (ips:ipe,kps:kpe,jps:jpe ) = 0.
 DO i=ips,ipe
   DO j=jps,jpe
   IF ((ic_flashrate(i,j) .gt. 0) .OR. (cg_flashrate(i,j) .gt. 0)) THEN
   cldfra_cup_1d(kps:kpe) = cldfra_cup(i,kps:kpe,j)
   htop_ij = htop(i,j)
   t_1d(kps:kpe) = t_phy(i,kps:kpe,j)
   z_1d(kps:kpe) = z(i,kps:kpe,j)
   ic_fr = ic_flashrate(i,j)
   cg_fr = cg_flashrate(i,j)
   CALL kfind ( cldfra_cup_1d, htop_ij, t_1d, &
                 ltng_temp_upper,ltng_temp_lower, &
                 i, j, &
                 ips, ipe, jps, jpe, kps, kpe, &
                 ktop,kbtm,kupper,klower )
   molesofair(kps:kpe) = rho(i,kps:kpe,j) * 1E3 * dx * dy / .02897
   IF (( ic_fr .gt. 0 ) .and. (( ktop .gt. klower ) .and. (kbtm .lt. ktop) ) )THEN
     call bellcurve(kbtm,ktop,klower,z_1d, kps,kpe, fd, dz)
     if (ktop .gt. kupper) then
       call bellcurve(kbtm,ktop,kupper,z_1d, kps,kpe, fd2, dz)
       fd(kbtm:ktop) = 0.5*( fd(kbtm:ktop) + fd2(kbtm:ktop) )
     endif
     DO k=kbtm,ktop
       IF ( cldfra_cup_1d(k) .gt. 0.01 ) THEN
         delta = (ic_fr * N_IC) * fd(k) / (molesofair(k)*dz(k)) * 1E6
         WRITE(message, * ) ' LNOx_driver: k, delta, fd = ', k, delta, fd(k)
         CALL wrf_debug ( 100, message )
         lnox_ic_tend(i,k,j) = delta
       ENDIF
     ENDDO
  ENDIF
  IF ((cg_fr .gt. 0 ) .and. (( ktop .gt. klower ) .and. (kbtm .lt. ktop) ) ) THEN
    call bellcurve(kps,ktop,klower,z_1d, kps,kpe, fd, dz)
    delta = (cg_fr * N_CG ) * fd(k) / (molesofair(k)*dz(k)) * 1E6
    IF ( cldfra_cup_1d(kbtm) .gt. 0.01 ) THEN
      lnox_cg_tend(i,kps:(kbtm-1),j) = delta
    ENDIF
    do k = kbtm,ktop
      delta = (cg_fr * N_CG ) * fd(k) / (molesofair(k)*dz(k)) * 1E6
      lnox_cg_tend(i,k,j) = delta
    enddo
  ENDIF
  ENDIF
  ENDDO
ENDDO
 END SUBROUTINE lightning_nox_decaria_cup
 SUBROUTINE bellcurve ( k_min, k_max, k_mu, z, kps,kpe, f, dz )
 IMPLICIT NONE
 INTEGER, INTENT(IN ) :: k_min, k_max, k_mu
 REAL, DIMENSION( kps:kpe ), INTENT(IN ) :: z
 INTEGER, INTENT(IN ) :: kps,kpe
 REAL, DIMENSION( kps:kpe ), INTENT( OUT) :: f, dz
 INTEGER :: i,j,k
 REAL, DIMENSION( kps:kpe ) :: ex
 REAL :: sigma, z_mu, cuml_f_dist
 REAL, PARAMETER :: two_pi = 6.2831854
 f(kps:kpe) = 0.
 z_mu = z(k_mu)
 sigma = AMIN1(z(k_max)-z_mu,z_mu-z(k_min))/3.0
 ex(k_min:k_max) = (z(k_min:k_max)-z_mu)/sigma
 f(k_min:k_max) = (1.0/(sqrt(two_pi)*sigma))*exp(-ex(k_min:k_max)*ex(k_min:k_max)/2.0)
 dz(kps) = (z(kps+1) - z(kps))*.5
 dz(kpe) = (z(kpe) - z(kpe-1))*.5
 DO k=kps+1,kpe-1
   dz(k) = (z(k+1) - z(k-1))*.5
 ENDDO
 cuml_f_dist = sum(f(k_min:k_max) * dz(k_min:k_max))
 f(k_min:k_max) = f(k_min:k_max)*dz(k_min:k_max)/cuml_f_dist
 END SUBROUTINE bellcurve
 SUBROUTINE kfind ( &
                cldfra_cup_1d, htop_ij, t_1d, &
                ltng_temp_upper,ltng_temp_lower, &
                i , j, &
                ips, ipe, jps, jpe, kps, kpe, &
                ktop,kbtm,kupper,klower &
            )
 USE module_state_description
 USE module_model_constants
 USE module_dm, only: wrf_dm_max_real, wrf_dm_min_real, wrf_dm_sum_real
 IMPLICIT NONE
 REAL, DIMENSION( kps:kpe ), INTENT(IN ) :: cldfra_cup_1d
 REAL, DIMENSION( kps:kpe ), INTENT(IN ) :: t_1d
 REAL, INTENT(IN ) :: htop_ij
 REAL, INTENT(IN ) :: ltng_temp_lower, ltng_temp_upper
 INTEGER, INTENT(IN ) :: ips,ipe, jps,jpe, kps,kpe
 INTEGER, INTENT(IN ) :: i, j
 INTEGER, INTENT( OUT) :: ktop,kbtm,kupper,klower
 CHARACTER (LEN=250) :: message
 REAL :: ktop_r, kbtm_r, kupper_r, klower_r
 INTEGER :: k
 ktop = kps
 kbtm = kps
 kupper = kps
 klower = kps
 k = kpe
 DO WHILE ( cldfra_cup_1d(k) .le. 0.01 .and. k .gt. kps)
  k = k-1
 ENDDO
 ktop = k
 k = kps
 DO WHILE( cldfra_cup_1d(k) .le. 0.01 .and. k .le. ktop )
  k = k+1
 ENDDO
 kbtm = k
 k = kps
 DO WHILE ( t_1d(k) .gt. ltng_temp_lower + 273.15 .and. k .lt. kpe )
   k = k + 1
 ENDDO
 klower = k
 DO WHILE ( t_1d(k) .gt. ltng_temp_upper + 273.15 .and. k .lt. kpe )
   k = k + 1
 ENDDO
 kupper = k
 ktop = nint(htop_ij)
 WRITE(message, * ) ' LNOx_driver: kbtm, ktop, klower, kupper = ', kbtm, ktop, klower, kupper
 CALL wrf_debug ( 100, message )
 END SUBROUTINE kfind
 END MODULE module_lightning_nox_decaria_cup
