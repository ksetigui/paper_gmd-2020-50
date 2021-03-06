


MODULE module_diagnostics
CONTAINS
   SUBROUTINE diagnostic_output_calc( &
                      ids,ide, jds,jde, kds,kde, &
                      ims,ime, jms,jme, kms,kme, &
                      ips,ipe, jps,jpe, kps,kpe, &
                      i_start,i_end,j_start,j_end,kts,kte,num_tiles &
                     ,dpsdt,dmudt &
                     ,p8w,pk1m,mu_2,mu_2m &
                     ,u,v &
                     ,raincv,rainncv,rainc,rainnc &
                     ,i_rainc,i_rainnc &
                     ,hfx,sfcevp,lh &
                     ,ACSWUPT,ACSWUPTC,ACSWDNT,ACSWDNTC &
                     ,ACSWUPB,ACSWUPBC,ACSWDNB,ACSWDNBC &
                     ,ACLWUPT,ACLWUPTC,ACLWDNT,ACLWDNTC &
                     ,ACLWUPB,ACLWUPBC,ACLWDNB,ACLWDNBC &
                     ,I_ACSWUPT,I_ACSWUPTC,I_ACSWDNT,I_ACSWDNTC &
                     ,I_ACSWUPB,I_ACSWUPBC,I_ACSWDNB,I_ACSWDNBC &
                     ,I_ACLWUPT,I_ACLWUPTC,I_ACLWDNT,I_ACLWDNTC &
                     ,I_ACLWUPB,I_ACLWUPBC,I_ACLWDNB,I_ACLWDNBC &
                     ,dt,xtime,sbw,t2 &
                     ,diag_print &
                     ,bucket_mm, bucket_J &
                     ,prec_acc_c, prec_acc_nc, snow_acc_nc &
                     ,snowncv, prec_acc_dt, curr_secs &
                     ,nwp_diagnostics, diagflag &
                     ,history_interval &
                     ,itimestep &
                     ,u10,v10,w &
                     ,wspd10max &
                     ,up_heli_max &
                     ,w_up_max,w_dn_max &
                     ,znw,w_colmean &
                     ,numcolpts,w_mean &
                     ,grpl_max,grpl_colint,refd_max,refl_10cm &
                     ,qg_curr &
                     ,rho,ph,phb,g &
                                                                      )


  USE module_dm, ONLY: wrf_dm_sum_real, wrf_dm_maxval

   IMPLICIT NONE
   INTEGER, INTENT(IN ) :: &
                                      ids,ide, jds,jde, kds,kde, &
                                      ims,ime, jms,jme, kms,kme, &
                                      ips,ipe, jps,jpe, kps,kpe, &
                                                        kts,kte, &
                                                      num_tiles
   INTEGER, DIMENSION(num_tiles), INTENT(IN) :: &
     & i_start,i_end,j_start,j_end
   INTEGER, INTENT(IN ) :: diag_print
   REAL, INTENT(IN ) :: bucket_mm, bucket_J
   INTEGER, INTENT(IN ) :: nwp_diagnostics
   LOGICAL, INTENT(IN ) :: diagflag
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), &
         INTENT(IN ) :: u &
                                                    , v &
                                                    , p8w
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: &
                                                           MU_2 &
                                                    , RAINNCV &
                                                    , RAINCV &
                                                    , SNOWNCV &
                                                    , HFX &
                                                    , LH &
                                                    , SFCEVP &
                                                    , T2
   REAL, DIMENSION( ims:ime , jms:jme ), &
          INTENT(INOUT) :: DPSDT &
                                                    , DMUDT &
                                                    , RAINNC &
                                                    , RAINC &
                                                    , MU_2M &
                                                    , PK1M
   REAL, INTENT(IN ) :: DT, XTIME
   INTEGER, INTENT(IN ) :: SBW
   INTEGER, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT) :: &
                                                       I_RAINC, &
                                                       I_RAINNC
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT) ::&
                      ACSWUPT,ACSWUPTC,ACSWDNT,ACSWDNTC, &
                      ACSWUPB,ACSWUPBC,ACSWDNB,ACSWDNBC, &
                      ACLWUPT,ACLWUPTC,ACLWDNT,ACLWDNTC, &
                      ACLWUPB,ACLWUPBC,ACLWDNB,ACLWDNBC
   INTEGER, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT) ::&
                      I_ACSWUPT,I_ACSWUPTC,I_ACSWDNT,I_ACSWDNTC, &
                      I_ACSWUPB,I_ACSWUPBC,I_ACSWDNB,I_ACSWDNBC, &
                      I_ACLWUPT,I_ACLWUPTC,I_ACLWDNT,I_ACLWDNTC, &
                      I_ACLWUPB,I_ACLWUPBC,I_ACLWDNB,I_ACLWDNBC
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT) ::&
                      PREC_ACC_C, PREC_ACC_NC, SNOW_ACC_NC
   REAL, OPTIONAL, INTENT(IN):: PREC_ACC_DT, CURR_SECS
   INTEGER :: i,j,k,its,ite,jts,jte,ij
   INTEGER :: idp,jdp,irc,jrc,irnc,jrnc,isnh,jsnh
   INTEGER :: prfreq
   REAL :: no_points
   REAL :: dpsdt_sum, dmudt_sum, dardt_sum, drcdt_sum, drndt_sum
   REAL :: hfx_sum, lh_sum, sfcevp_sum, rainc_sum, rainnc_sum, raint_sum
   REAL :: dmumax, raincmax, rainncmax, snowhmax
   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
   CHARACTER*256 :: outstring
   CHARACTER*6 :: grid_str
   INTEGER, INTENT(IN) :: &
                                     history_interval,itimestep
   REAL, DIMENSION( kms:kme ), INTENT(IN) :: &
                                                            znw
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN) :: &
                                                              w &
                                                       ,qg_curr &
                                                           ,rho &
                                                     ,refl_10cm &
                                                        ,ph,phb
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN) :: &
                                                            u10 &
                                                           ,v10
   REAL, INTENT(IN) :: g
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT) :: &
                                                      wspd10max &
                                                   ,up_heli_max &
                                             ,w_up_max,w_dn_max &
                                    ,w_colmean,numcolpts,w_mean &
                                          ,grpl_max,grpl_colint &
                                                      ,refd_max
   INTEGER :: idump
   REAL :: wind_vel
   REAL :: depth
   IF(bucket_mm .gt. 0. .AND. MOD(NINT(XTIME),360) .EQ. 0)THEN
! !$OMP PARALLEL DO &
! !$OMP PRIVATE ( ij )
   DO ij = 1 , num_tiles
      IF (xtime .eq. 0.0)THEN
        DO j=j_start(ij),j_end(ij)
        DO i=i_start(ij),i_end(ij)
          i_rainnc(i,j) = 0
          i_rainc(i,j) = 0
        ENDDO
        ENDDO
      ENDIF
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
        IF(rainnc(i,j) .gt. bucket_mm)THEN
          rainnc(i,j) = rainnc(i,j) - bucket_mm
          i_rainnc(i,j) = i_rainnc(i,j) + 1
        ENDIF
        IF(rainc(i,j) .gt. bucket_mm)THEN
          rainc(i,j) = rainc(i,j) - bucket_mm
          i_rainc(i,j) = i_rainc(i,j) + 1
        ENDIF
      ENDDO
      ENDDO
      IF (xtime .eq. 0.0 .and. bucket_J .gt. 0.0 .and. PRESENT(ACSWUPT))THEN
        DO j=j_start(ij),j_end(ij)
        DO i=i_start(ij),i_end(ij)
          i_acswupt(i,j) = 0
          i_acswuptc(i,j) = 0
          i_acswdnt(i,j) = 0
          i_acswdntc(i,j) = 0
          i_acswupb(i,j) = 0
          i_acswupbc(i,j) = 0
          i_acswdnb(i,j) = 0
          i_acswdnbc(i,j) = 0
        ENDDO
        ENDDO
      ENDIF
      IF (xtime .eq. 0.0 .and. bucket_J .gt. 0.0 .and. PRESENT(ACLWUPT))THEN
        DO j=j_start(ij),j_end(ij)
        DO i=i_start(ij),i_end(ij)
          i_aclwupt(i,j) = 0
          i_aclwuptc(i,j) = 0
          i_aclwdnt(i,j) = 0
          i_aclwdntc(i,j) = 0
          i_aclwupb(i,j) = 0
          i_aclwupbc(i,j) = 0
          i_aclwdnb(i,j) = 0
          i_aclwdnbc(i,j) = 0
        ENDDO
        ENDDO
      ENDIF
      IF (PRESENT(ACSWUPT) .and. bucket_J .gt. 0.0)THEN
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
        IF(acswupt(i,j) .gt. bucket_J)THEN
          acswupt(i,j) = acswupt(i,j) - bucket_J
          i_acswupt(i,j) = i_acswupt(i,j) + 1
        ENDIF
        IF(acswuptc(i,j) .gt. bucket_J)THEN
          acswuptc(i,j) = acswuptc(i,j) - bucket_J
          i_acswuptc(i,j) = i_acswuptc(i,j) + 1
        ENDIF
        IF(acswdnt(i,j) .gt. bucket_J)THEN
          acswdnt(i,j) = acswdnt(i,j) - bucket_J
          i_acswdnt(i,j) = i_acswdnt(i,j) + 1
        ENDIF
        IF(acswdntc(i,j) .gt. bucket_J)THEN
          acswdntc(i,j) = acswdntc(i,j) - bucket_J
          i_acswdntc(i,j) = i_acswdntc(i,j) + 1
        ENDIF
        IF(acswupb(i,j) .gt. bucket_J)THEN
          acswupb(i,j) = acswupb(i,j) - bucket_J
          i_acswupb(i,j) = i_acswupb(i,j) + 1
        ENDIF
        IF(acswupbc(i,j) .gt. bucket_J)THEN
          acswupbc(i,j) = acswupbc(i,j) - bucket_J
          i_acswupbc(i,j) = i_acswupbc(i,j) + 1
        ENDIF
        IF(acswdnb(i,j) .gt. bucket_J)THEN
          acswdnb(i,j) = acswdnb(i,j) - bucket_J
          i_acswdnb(i,j) = i_acswdnb(i,j) + 1
        ENDIF
        IF(acswdnbc(i,j) .gt. bucket_J)THEN
          acswdnbc(i,j) = acswdnbc(i,j) - bucket_J
          i_acswdnbc(i,j) = i_acswdnbc(i,j) + 1
        ENDIF
      ENDDO
      ENDDO
      ENDIF
      IF (PRESENT(ACLWUPT) .and. bucket_J .gt. 0.0)THEN
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
        IF(aclwupt(i,j) .gt. bucket_J)THEN
          aclwupt(i,j) = aclwupt(i,j) - bucket_J
          i_aclwupt(i,j) = i_aclwupt(i,j) + 1
        ENDIF
        IF(aclwuptc(i,j) .gt. bucket_J)THEN
          aclwuptc(i,j) = aclwuptc(i,j) - bucket_J
          i_aclwuptc(i,j) = i_aclwuptc(i,j) + 1
        ENDIF
        IF(aclwdnt(i,j) .gt. bucket_J)THEN
          aclwdnt(i,j) = aclwdnt(i,j) - bucket_J
          i_aclwdnt(i,j) = i_aclwdnt(i,j) + 1
        ENDIF
        IF(aclwdntc(i,j) .gt. bucket_J)THEN
          aclwdntc(i,j) = aclwdntc(i,j) - bucket_J
          i_aclwdntc(i,j) = i_aclwdntc(i,j) + 1
        ENDIF
        IF(aclwupb(i,j) .gt. bucket_J)THEN
          aclwupb(i,j) = aclwupb(i,j) - bucket_J
          i_aclwupb(i,j) = i_aclwupb(i,j) + 1
        ENDIF
        IF(aclwupbc(i,j) .gt. bucket_J)THEN
          aclwupbc(i,j) = aclwupbc(i,j) - bucket_J
          i_aclwupbc(i,j) = i_aclwupbc(i,j) + 1
        ENDIF
        IF(aclwdnb(i,j) .gt. bucket_J)THEN
          aclwdnb(i,j) = aclwdnb(i,j) - bucket_J
          i_aclwdnb(i,j) = i_aclwdnb(i,j) + 1
        ENDIF
        IF(aclwdnbc(i,j) .gt. bucket_J)THEN
          aclwdnbc(i,j) = aclwdnbc(i,j) - bucket_J
          i_aclwdnbc(i,j) = i_aclwdnbc(i,j) + 1
        ENDIF
      ENDDO
      ENDDO
      ENDIF
   ENDDO
! !$OMP END PARALLEL DO
   ENDIF
   IF (prec_acc_dt .gt. 0.) THEN
! !$OMP PARALLEL DO &
! !$OMP PRIVATE ( ij )
   DO ij = 1 , num_tiles
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
         IF (mod(curr_secs, 60.* prec_acc_dt) == 0.) THEN
            prec_acc_c(i,j) = 0.
            prec_acc_nc(i,j) = 0.
            snow_acc_nc(i,j) = 0.
         ENDIF
         prec_acc_c(i,j) = prec_acc_c(i,j) + RAINCV(i,j)
         prec_acc_nc(i,j) = prec_acc_nc(i,j) + RAINNCV(i,j)
         prec_acc_c(i,j) = MAX (prec_acc_c(i,j), 0.0)
         prec_acc_nc(i,j) = MAX (prec_acc_nc(i,j), 0.0)
         snow_acc_nc(i,j) = snow_acc_nc(i,j) + SNOWNCV(I,J)
         IF ( t2(i,j) .lt. 273.15 ) THEN
         snow_acc_nc(i,j) = snow_acc_nc(i,j) + RAINCV(i,j)
         snow_acc_nc(i,j) = MAX (snow_acc_nc(i,j), 0.0)
         ENDIF
      ENDDO
      ENDDO
   ENDDO
! !$OMP END PARALLEL DO
   ENDIF
   IF ( nwp_diagnostics .EQ. 1 ) THEN
     idump = (history_interval * 60.) / dt
   IF ( MOD((itimestep - 1), idump) .eq. 0 ) THEN
     WRITE(outstring,*) 'NSSL Diagnostics: Resetting max arrays for domain with dt = ', dt
     CALL wrf_debug ( 10,TRIM(outstring) )
! !$OMP PARALLEL DO &
! !$OMP PRIVATE ( ij )
     DO ij = 1 , num_tiles
       DO j=j_start(ij),j_end(ij)
       DO i=i_start(ij),i_end(ij)
         wspd10max(i,j) = 0.
         up_heli_max(i,j) = 0.
         w_up_max(i,j) = 0.
         w_dn_max(i,j) = 0.
         w_mean(i,j) = 0.
         grpl_max(i,j) = 0.
         refd_max(i,j) = 0.
       ENDDO
       ENDDO
     ENDDO
! !$OMP END PARALLEL DO
   ENDIF
! !$OMP PARALLEL DO &
! !$OMP PRIVATE ( ij )
   DO ij = 1 , num_tiles
     DO j=j_start(ij),j_end(ij)
     DO i=i_start(ij),i_end(ij)
       w_colmean(i,j) = 0.
       numcolpts(i,j) = 0.
       grpl_colint(i,j) = 0.
     ENDDO
     ENDDO
   ENDDO
! !$OMP END PARALLEL DO
! !$OMP PARALLEL DO &
! !$OMP PRIVATE ( ij )
   DO ij = 1 , num_tiles
     DO j=j_start(ij),j_end(ij)
     DO k=kms,kme
     DO i=i_start(ij),i_end(ij)
       IF ( p8w(i,k,j) .GT. 40000. .AND. w(i,k,j) .GT. w_up_max(i,j) ) THEN
         w_up_max(i,j) = w(i,k,j)
       ENDIF
       IF ( p8w(i,k,j) .GT. 40000. .AND. w(i,k,j) .LT. w_dn_max(i,j) ) THEN
         w_dn_max(i,j) = w(i,k,j)
       ENDIF
       IF ( znw(k) .GE. 0.5 .AND. znw(k) .LE. 0.8 ) THEN
         w_colmean(i,j) = w_colmean(i,j) + w(i,k,j)
         numcolpts(i,j) = numcolpts(i,j) + 1
       ENDIF
     ENDDO
     ENDDO
     ENDDO
   ENDDO
! !$OMP END PARALLEL DO
! !$OMP PARALLEL DO &
! !$OMP PRIVATE ( ij )
   DO ij = 1 , num_tiles
     DO j=j_start(ij),j_end(ij)
     DO k=kms,kme-1
     DO i=i_start(ij),i_end(ij)
       depth = ( ( ph(i,k+1,j) + phb(i,k+1,j) ) / g ) - &
               ( ( ph(i,k ,j) + phb(i,k ,j) ) / g )
       grpl_colint(i,j) = grpl_colint(i,j) + qg_curr(i,k,j) * depth * rho(i,k,j)
     ENDDO
     ENDDO
     ENDDO
   ENDDO
! !$OMP END PARALLEL DO
! !$OMP PARALLEL DO &
! !$OMP PRIVATE ( ij )
   DO ij = 1 , num_tiles
     DO j=j_start(ij),j_end(ij)
     DO i=i_start(ij),i_end(ij)
       wind_vel = sqrt ( u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j) )
       IF ( wind_vel .GT. wspd10max(i,j) ) THEN
         wspd10max(i,j) = wind_vel
       ENDIF
       w_mean(i,j) = w_mean(i,j) + w_colmean(i,j) / numcolpts(i,j)
       IF ( MOD(itimestep, idump) .eq. 0 ) THEN
         w_mean(i,j) = w_mean(i,j) / idump
       ENDIF
       IF ( grpl_colint(i,j) .gt. grpl_max(i,j) ) THEN
          grpl_max(i,j) = grpl_colint(i,j)
       ENDIF
       IF ( refl_10cm(i,kms,j) .GT. refd_max(i,j) ) THEN
         refd_max(i,j) = refl_10cm(i,kms,j)
       ENDIF
     ENDDO
     ENDDO
   ENDDO
! !$OMP END PARALLEL DO
   ENDIF
   if (diag_print .eq. 0 ) return
   IF ( xtime .ne. 0. ) THEN
    if(diag_print.eq.1) then
       prfreq = dt
    else
       prfreq=10
    endif
    IF (MOD(nint(dt),prfreq) == 0) THEN
   no_points = float((ide-ids)*(jde-jds))
! !$OMP PARALLEL DO &
! !$OMP PRIVATE ( ij )
   dmumax = 0.
   DO ij = 1 , num_tiles
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
         dpsdt(i,j)=(p8w(i,kms,j)-pk1m(i,j))/dt
         dmudt(i,j)=(mu_2(i,j)-mu_2m(i,j))/dt
         if(abs(dmudt(i,j)*dt).gt.dmumax)then
           dmumax=abs(dmudt(i,j)*dt)
           idp=i
           jdp=j
         endif
      ENDDO
      ENDDO
   ENDDO
! !$OMP END PARALLEL DO
   dmumax = dmumax*1.e-5
   CALL wrf_dm_maxval ( dmumax, idp, jdp )
   dpsdt_sum = 0.
   dmudt_sum = 0.
   DO j = jps, min(jpe,jde-1)
     DO i = ips, min(ipe,ide-1)
       dpsdt_sum = dpsdt_sum + abs(dpsdt(i,j))
       dmudt_sum = dmudt_sum + abs(dmudt(i,j))
     ENDDO
   ENDDO
   dpsdt_sum = wrf_dm_sum_real ( dpsdt_sum )
   dmudt_sum = wrf_dm_sum_real ( dmudt_sum )
   IF ( diag_print .eq. 2 ) THEN
   dardt_sum = 0.
   drcdt_sum = 0.
   drndt_sum = 0.
   rainc_sum = 0.
   raint_sum = 0.
   rainnc_sum = 0.
   sfcevp_sum = 0.
   hfx_sum = 0.
   lh_sum = 0.
   raincmax = 0.
   rainncmax = 0.
   DO j = jps, min(jpe,jde-1)
     DO i = ips, min(ipe,ide-1)
       drcdt_sum = drcdt_sum + abs(raincv(i,j))
       drndt_sum = drndt_sum + abs(rainncv(i,j))
       dardt_sum = dardt_sum + abs(raincv(i,j)) + abs(rainncv(i,j))
       rainc_sum = rainc_sum + abs(rainc(i,j))
       IF(rainc(i,j).gt.raincmax)then
          raincmax=rainc(i,j)
          irc=i
          jrc=j
       ENDIF
       rainnc_sum = rainnc_sum + abs(rainnc(i,j))
       IF(rainnc(i,j).gt.rainncmax)then
          rainncmax=rainnc(i,j)
          irnc=i
          jrnc=j
       ENDIF
       raint_sum = raint_sum + abs(rainc(i,j)) + abs(rainnc(i,j))
       sfcevp_sum = sfcevp_sum + abs(sfcevp(i,j))
       hfx_sum = hfx_sum + abs(hfx(i,j))
       lh_sum = lh_sum + abs(lh(i,j))
     ENDDO
   ENDDO
   CALL wrf_dm_maxval ( raincmax, irc, jrc )
   CALL wrf_dm_maxval ( rainncmax, irnc, jrnc )
   drcdt_sum = wrf_dm_sum_real ( drcdt_sum )
   drndt_sum = wrf_dm_sum_real ( drndt_sum )
   dardt_sum = wrf_dm_sum_real ( dardt_sum )
   rainc_sum = wrf_dm_sum_real ( rainc_sum )
   rainnc_sum = wrf_dm_sum_real ( rainnc_sum )
   raint_sum = wrf_dm_sum_real ( raint_sum )
   sfcevp_sum = wrf_dm_sum_real ( sfcevp_sum )
   hfx_sum = wrf_dm_sum_real ( hfx_sum )
   lh_sum = wrf_dm_sum_real ( lh_sum )
   ENDIF
   CALL get_current_grid_name( grid_str )
   IF ( wrf_dm_on_monitor() ) THEN
     WRITE(outstring,*) grid_str,'Domain average of dpsdt, dmudt (mb/3h): ', xtime, &
           dpsdt_sum/no_points*108., &
           dmudt_sum/no_points*108.
     CALL wrf_message ( TRIM(outstring) )
     WRITE(outstring,*) grid_str,'Max mu change time step: ', idp,jdp,dmumax
     CALL wrf_message ( TRIM(outstring) )
     IF ( diag_print .eq. 2) THEN
     WRITE(outstring,*) grid_str,'Domain average of dardt, drcdt, drndt (mm/sec): ', xtime, &
           dardt_sum/dt/no_points, &
           drcdt_sum/dt/no_points, &
           drndt_sum/dt/no_points
     CALL wrf_message ( TRIM(outstring) )
     WRITE(outstring,*) grid_str,'Domain average of rt_sum, rc_sum, rnc_sum (mm): ', xtime, &
           raint_sum/no_points, &
           rainc_sum/no_points, &
           rainnc_sum/no_points
     CALL wrf_message ( TRIM(outstring) )
     WRITE(outstring,*) grid_str,'Max Accum Resolved Precip,   I,J  (mm): ' ,&
           rainncmax,irnc,jrnc
     CALL wrf_message ( TRIM(outstring) )
     WRITE(outstring,*) grid_str,'Max Accum Convective Precip,   I,J  (mm): ' ,&
           raincmax,irc,jrc
     CALL wrf_message ( TRIM(outstring) )
     WRITE(outstring,*) grid_str,'Domain average of sfcevp, hfx, lh: ', xtime, &
           sfcevp_sum/no_points, &
           hfx_sum/no_points, &
           lh_sum/no_points
     CALL wrf_message ( TRIM(outstring) )
     ENDIF
   ENDIF
    ENDIF
   ENDIF
   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij,i,j )
   DO ij = 1 , num_tiles
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
         pk1m(i,j)=p8w(i,kms,j)
         mu_2m(i,j)=mu_2(i,j)
      ENDDO
      ENDDO
      IF ( xtime .lt. 0.0001 ) THEN
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
         dpsdt(i,j)=0.
         dmudt(i,j)=0.
      ENDDO
      ENDDO
      ENDIF
   ENDDO
   !$OMP END PARALLEL DO
   END SUBROUTINE diagnostic_output_calc
   SUBROUTINE clwrf_output_calc( &
                      ids,ide, jds,jde, kds,kde, &
                      ims,ime, jms,jme, kms,kme, &
                      ips,ipe, jps,jpe, kps,kpe, &
                      i_start,i_end,j_start,j_end,kts,kte,num_tiles &
                     ,dpsdt,dmudt &
                     ,p8w,pk1m,mu_2,mu_2m &
                     ,u,v &
                     ,is_restart &
                     ,clwrfH,t2,q2,u10,v10, skintemp &
                     ,t2clmin,t2clmax,tt2clmin,tt2clmax &
                     ,t2clmean,t2clstd &
                     ,q2clmin,q2clmax,tq2clmin,tq2clmax &
                     ,q2clmean,q2clstd &
                     ,u10clmax,v10clmax,spduv10clmax,tspduv10clmax &
                     ,u10clmean,v10clmean,spduv10clmean &
                     ,u10clstd,v10clstd,spduv10clstd &
                     ,raincclmax,rainncclmax,traincclmax,trainncclmax &
                     ,raincclmean,rainncclmean,raincclstd,rainncclstd &
                     ,skintempclmin,skintempclmax &
                     ,tskintempclmin,tskintempclmax &
                     ,skintempclmean,skintempclstd &
                     ,raincv,rainncv,rainc,rainnc &
                     ,i_rainc,i_rainnc &
                     ,hfx,sfcevp,lh &
                     ,dt,xtime,sbw &
                     ,diag_print &
                     ,bucket_mm, bucket_J &
                                                                      )
  USE module_dm, ONLY: wrf_dm_sum_real, wrf_dm_maxval
  USE module_configure
   IMPLICIT NONE
   INTEGER, INTENT(IN ) :: &
                                      ids,ide, jds,jde, kds,kde, &
                                      ims,ime, jms,jme, kms,kme, &
                                      ips,ipe, jps,jpe, kps,kpe, &
                                                        kts,kte, &
                                                      num_tiles
   INTEGER, DIMENSION(num_tiles), INTENT(IN) :: i_start, &
                                      i_end,j_start,j_end
   INTEGER, INTENT(IN ) :: diag_print
   REAL, INTENT(IN ) :: bucket_mm, &
                                      bucket_J
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), &
                                       INTENT(IN ) :: u,v,p8w
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN) :: MU_2, &
                                      RAINNCV, RAINCV, HFX, &
                                      SFCEVP, LH, SKINTEMP
   REAL, DIMENSION( ims:ime , jms:jme ), &
                                     INTENT(INOUT) :: DPSDT, &
                                     DMUDT, RAINNC, RAINC, &
                                     MU_2M, PK1M
   REAL, INTENT(IN ) :: DT, XTIME
   INTEGER, INTENT(IN ) :: SBW
   INTEGER, DIMENSION( ims:ime , jms:jme ), &
                                     INTENT(INOUT) :: I_RAINC, &
                                     I_RAINNC
   INTEGER :: i,j,k,its,ite,jts,jte,ij
   INTEGER :: idp,jdp,irc,jrc,irnc,jrnc,isnh,jsnh
   INTEGER :: prfreq
   REAL :: dpsdt_sum, dmudt_sum, dardt_sum, &
                          drcdt_sum, drndt_sum
   REAL :: hfx_sum, lh_sum, sfcevp_sum, &
                          rainc_sum, rainnc_sum, raint_sum
   REAL :: dmumax, raincmax, rainncmax, &
                          snowhmax
   REAL :: xtimep
   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
   CHARACTER*256 :: outstring
   CHARACTER*6 :: grid_str
   CHARACTER (LEN=80) :: timestr
   REAL, DIMENSION( ims:ime , jms:jme ), &
                          INTENT(IN) :: t2, q2, u10, v10
   REAL, DIMENSION( ims:ime , jms:jme ), &
                          INTENT(OUT) :: t2clmin, t2clmax, tt2clmin, &
                          tt2clmax, t2clmean, t2clstd, &
                          q2clmin, q2clmax, tq2clmin, tq2clmax, q2clmean, q2clstd,&
                          u10clmax, v10clmax, spduv10clmax, tspduv10clmax, &
                          u10clmean, v10clmean, spduv10clmean, &
                          u10clstd, v10clstd, spduv10clstd, skintempclmin, &
                          skintempclmax, tskintempclmin, tskintempclmax, &
                          skintempclmean, skintempclstd
   REAL, DIMENSION( ims:ime , jms:jme ), &
                          INTENT(OUT) :: raincclmax, rainncclmax, &
                          traincclmax, trainncclmax, raincclmean, rainncclmean, &
                          raincclstd, rainncclstd
   REAL, PARAMETER :: minimum0= 1000000., &
                          maximum0= -1000000.
   REAL :: value
   INTEGER, INTENT(IN) :: clwrfH
   CHARACTER (LEN=1024) :: message
   REAL, SAVE :: nsteps
   LOGICAL :: is_restart
! !$OMP PARALLEL DO &
! !$OMP PRIVATE ( ij )
  IF (( MOD(NINT(XTIME/dt*60.),NINT(clwrfH/dt*60.)) == 0) .AND. (.NOT.is_restart)) THEN
    DO ij = 1 , num_tiles
      IF ( wrf_dm_on_monitor() ) THEN
        WRITE(message, *)'CLWRFdiag - T2; tile: ',ij,' T2clmin:', &
          t2clmin(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' T2clmax:', &
          t2clmax(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' TT2clmin:', &
          tt2clmin(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' TT2clmax:', &
          tt2clmax(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' T2clmean:', &
          t2clmean(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' T2clstd:', &
          t2clstd(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2)
        CALL wrf_debug(75, message)
        WRITE(message, *)'CLWRFdiag - Q2; tile: ',ij,' Q2clmin:', &
          q2clmin(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' Q2clmax:', &
          q2clmax(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' TQ2clmin:', &
          tq2clmin(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' TQ2clmax:', &
          tq2clmax(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' Q2clmean:', &
          q2clmean(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' Q2clstd:', &
          q2clstd(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2)
        CALL wrf_debug(75, message)
        WRITE(message, *)'CLWRFdiag - WINDSPEED; tile: ',ij,' U10clmax:', &
          u10clmax(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' V10clmax:', &
          v10clmax(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' SPDUV10clmax:', &
          spduv10clmax(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' TSPDUV10clmax:', &
          tspduv10clmax(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' U10clmean:', &
          u10clmean(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' V10clmean:', &
          v10clmean(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' SPDUV10clmean:', &
          spduv10clmean(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' U10clstd:', &
          u10clstd(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' V10clstd:', &
          v10clstd(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' SPDUV10clstd:', &
          spduv10clstd(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2)
        CALL wrf_debug(75, message)
        WRITE(message, *)'CLWRFdiag - RAIN; tile: ',ij,' RAINCclmax:', &
          raincclmax(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' RAINNCclmax:', &
          rainncclmax(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' TRAINCclmax:', &
          traincclmax(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' TRAINNCclmax:', &
          trainncclmax(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' RAINCclmean:', &
          raincclmean(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' RAINNCclmean:', &
          rainncclmean(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' RAINCclstd:', &
          raincclstd(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' RAINNCclstd:', &
          rainncclstd(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2)
        CALL wrf_debug(75, message)
        WRITE(message,*)'CLWRFdiag - SKINTEMP; tile: ',ij,' SKINTEMPclmin:',&
          skintempclmin(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' SKINTEMPclmax:', &
          skintempclmax(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' TSKINTEMPclmin:', &
          tskintempclmin(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' TSKINTEMPclmax:', &
          tskintempclmax(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' SKINTEMPclmean:', &
          skintempclmean(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2),' SKINTEMPclstd:', &
          skintempclstd(i_start(ij)+(i_end(ij)-i_start(ij))/2, &
          j_start(ij)+(j_end(ij)-j_start(ij))/2)
        CALL wrf_debug(75, message)
      ENDIF
      DO j = j_start(ij), j_end(ij)
        DO i = i_start(ij), i_end(ij)
          t2clmin(i,j)=t2(i,j)
          t2clmax(i,j)=t2(i,j)
          t2clmean(i,j)=t2(i,j)
          t2clstd(i,j)=t2(i,j)*t2(i,j)
          q2clmin(i,j)=q2(i,j)
          q2clmax(i,j)=q2(i,j)
          q2clmean(i,j)=q2(i,j)
          q2clstd(i,j)=q2(i,j)*q2(i,j)
          spduv10clmax(i,j)=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
          u10clmean(i,j)=u10(i,j)
          v10clmean(i,j)=v10(i,j)
          spduv10clmean(i,j)=sqrt(u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j))
          u10clstd(i,j)=u10(i,j)*u10(i,j)
          v10clstd(i,j)=v10(i,j)*v10(i,j)
          spduv10clstd(i,j)=u10(i,j)*u10(i,j)+v10(i,j)*v10(i,j)
          raincclmax(i,j)=raincv(i,j)/dt
          rainncclmax(i,j)=rainncv(i,j)/dt
          raincclmean(i,j)=raincv(i,j)/dt
          rainncclmean(i,j)=rainncv(i,j)/dt
          raincclstd(i,j)=(raincv(i,j)/dt)*(raincv(i,j)/dt)
          rainncclstd(i,j)=(rainncv(i,j)/dt)*(rainncv(i,j)/dt)
          skintempclmin(i,j)=skintemp(i,j)
          skintempclmax(i,j)=skintemp(i,j)
          skintempclmean(i,j)=skintemp(i,j)
          skintempclstd(i,j)=skintemp(i,j)*skintemp(i,j)
        ENDDO
      ENDDO
  ENDDO
    nsteps=clwrfH*60./dt
  ELSE
    xtimep = xtime + dt/60.
    nsteps=clwrfH*60./dt
          CALL varstatistics(t2,xtimep,ime-ims+1,jme-jms+1,t2clmin,t2clmax, &
            tt2clmin,tt2clmax,t2clmean,t2clstd)
          CALL varstatistics(q2,xtimep,ime-ims+1,jme-jms+1,q2clmin,q2clmax, &
            tq2clmin,tq2clmax,q2clmean,q2clstd)
          CALL varstatisticsWIND(u10,v10,xtimep,ime-ims+1,jme-jms+1,u10clmax, &
            v10clmax,spduv10clmax,tspduv10clmax,u10clmean,v10clmean, &
            spduv10clmean,u10clstd,v10clstd,spduv10clstd)
          CALL varstatisticsMAX(raincv/dt,xtimep,ime-ims+1,jme-jms+1, &
            raincclmax,traincclmax,raincclmean,raincclstd)
          CALL varstatisticsMAX(rainncv/dt,xtimep,ime-ims+1,jme-jms+1, &
            rainncclmax,trainncclmax,rainncclmean,rainncclstd)
          CALL varstatistics(skintemp,xtimep,ime-ims+1,jme-jms+1,skintempclmin,&
            skintempclmax, tskintempclmin,tskintempclmax,skintempclmean, &
            skintempclstd)
           IF ((MOD(NINT((XTIME+dt/60.)*60./dt),NINT(clwrfH*60./dt)) == 0)) THEN
             IF ( wrf_dm_on_monitor() ) PRINT *,'nsteps=',nsteps,' xtime:', &
               xtime,' clwrfH:',clwrfH
               t2clmean=t2clmean/nsteps
               t2clstd=SQRT(t2clstd/nsteps-t2clmean**2.)
               q2clmean=q2clmean/nsteps
               q2clstd=SQRT(q2clstd/nsteps-q2clmean**2.)
               u10clmean=u10clmean/nsteps
               v10clmean=v10clmean/nsteps
               spduv10clmean=spduv10clmean/nsteps
               u10clstd=SQRT(u10clstd/nsteps-u10clmean**2.)
               v10clstd=SQRT(v10clstd/nsteps-v10clmean**2.)
               spduv10clstd=SQRT(spduv10clstd/nsteps- &
                 spduv10clmean**2)
               raincclmean=raincclmean/nsteps
               rainncclmean=rainncclmean/nsteps
               raincclstd=SQRT(raincclstd/nsteps-raincclmean**2.)
               rainncclstd=SQRT(rainncclstd/nsteps-rainncclmean**2.)
               skintempclmean=skintempclmean/nsteps
              skintempclstd=SQRT(skintempclstd/nsteps-skintempclmean**2.)
            END IF
  ENDIF
! !$OMP END PARALLEL DO
   END SUBROUTINE clwrf_output_calc
SUBROUTINE varstatisticsWIND(varu, varv, tt, dx, dy, varumax, varvmax, &
  varuvmax, tvaruvmax, varumean, varvmean, varuvmean, varustd, varvstd, &
  varuvstd)
IMPLICIT NONE
INTEGER :: i, j
INTEGER, INTENT(IN) :: dx, dy
REAL, DIMENSION(dx,dy), INTENT(IN) :: varu, varv
REAL, INTENT(IN) :: tt
REAL, DIMENSION(dx,dy), INTENT(INOUT) :: varumax, &
  varvmax, varuvmax, tvaruvmax, varumean, varvmean, varuvmean, varustd, &
  varvstd, varuvstd
REAL :: varuv
DO i=1,dx
  DO j=1,dy
    varuv=sqrt(varu(i,j)*varu(i,j)+varv(i,j)*varv(i,j))
      IF (varuv > varuvmax(i,j)) THEN
        varumax(i,j)=varu(i,j)
        varvmax(i,j)=varv(i,j)
        varuvmax(i,j)=varuv
        tvaruvmax(i,j)=tt
      END IF
    varuvmean(i,j)=varuvmean(i,j)+varuv
    varuvstd(i,j)=varuvstd(i,j)+varuv**2
  END DO
END DO
varumean=varumean+varu
varvmean=varvmean+varv
varustd=varustd+varu**2
varvstd=varvstd+varv**2
END SUBROUTINE varstatisticsWIND
SUBROUTINE varstatisticsMAX(var, tt, dx, dy, varmax, tvarmax, varmean, &
   varstd)
IMPLICIT NONE
INTEGER :: i,j
INTEGER, INTENT(IN) :: dx, dy
REAL, DIMENSION(dx,dy), INTENT(IN) :: var
REAL, INTENT(IN) :: tt
REAL, DIMENSION(dx,dy), INTENT(INOUT) :: varmax, &
  tvarmax, varmean, varstd
DO i=1,dx
  DO j=1,dy
    IF (var(i,j) > varmax(i,j)) THEN
      varmax(i,j)=var(i,j)
      tvarmax(i,j)=tt
    END IF
  END DO
END DO
varmean=varmean+var
varstd=varstd+var**2
END SUBROUTINE varstatisticsMAX
SUBROUTINE varstatistics(var, tt, dx, dy, varmin, varmax, tvarmin, tvarmax, &
  varmean, varstd)
IMPLICIT NONE
INTEGER :: i,j
INTEGER, INTENT(IN) :: dx, dy
REAL, DIMENSION(dx,dy), INTENT(IN) :: var
REAL, INTENT(IN) :: tt
REAL, DIMENSION(dx,dy), INTENT(INOUT) :: varmin, &
  varmax, tvarmin, tvarmax, varmean, varstd
DO i=1,dx
  DO j=1,dy
    IF (var(i,j) < varmin(i,j)) THEN
      varmin(i,j)=var(i,j)
      tvarmin(i,j)=tt
    END IF
    IF (var(i,j) > varmax(i,j)) THEN
      varmax(i,j)=var(i,j)
      tvarmax(i,j)=tt
    END IF
  END DO
END DO
varmean=varmean+var
varstd=varstd+var**2
END SUBROUTINE varstatistics
SUBROUTINE pld ( u,v,w,t,qv,zp,zb,pp,pb,p,pw, &
                 msfux,msfuy,msfvx,msfvy,msftx,msfty, &
                 f,e, &
                 use_tot_or_hyd_p,missing, &
                 num_press_levels,max_press_levels,press_levels, &
                 p_pl,u_pl,v_pl,t_pl,rh_pl,ght_pl,s_pl,td_pl, &
                 ids,ide, jds,jde, kds,kde, &
                 ims,ime, jms,jme, kms,kme, &
                 its,ite, jts,jte, kts,kte )
   USE module_model_constants
   IMPLICIT NONE
   INTEGER, INTENT(IN ) :: ids,ide, jds,jde, kds,kde, &
                                                                      ims,ime, jms,jme, kms,kme, &
                                                                      its,ite, jts,jte, kts,kte
   REAL , INTENT(IN ) , DIMENSION(ims:ime , jms:jme) :: msfux,msfuy,msfvx,msfvy,msftx,msfty, &
                                                                      f,e
   INTEGER, INTENT(IN ) :: use_tot_or_hyd_p
   REAL , INTENT(IN ) :: missing
   REAL , INTENT(IN ) , DIMENSION(ims:ime , kms:kme , jms:jme) :: u,v,w,t,qv,zp,zb,pp,pb,p,pw
   INTEGER, INTENT(IN ) :: num_press_levels, max_press_levels
   REAL , INTENT(IN ) , DIMENSION(max_press_levels) :: press_levels
   REAL , INTENT( OUT) , DIMENSION(num_press_levels) :: p_pl
   REAL , INTENT( OUT) , DIMENSION(ims:ime , num_press_levels , jms:jme) :: u_pl,v_pl,t_pl,rh_pl,ght_pl,s_pl,td_pl
   REAL, PARAMETER :: eps = 0.622, t_kelvin = svpt0 , s1 = 243.5, s2 = svp2 , s3 = svp1*10., s4 = 611.0, s5 = 5418.12
   INTEGER :: i, j, ke, kp, ke_h, ke_f
   REAL :: pu, pd, pm , &
              tu, td , &
              su, sd , &
              uu, ud , &
              vu, vd , &
              zu, zd , &
              qu, qd, qm , &
              eu, ed, em , &
              du, dd
   REAL :: es, qs
   DO kp = 1 , num_press_levels
      p_pl(kp) = press_levels(kp)
   END DO
   DO j = jts , jte
      DO kp = 1 , num_press_levels
         DO i = its , ite
            u_pl (i,kp,j) = missing
            v_pl (i,kp,j) = missing
            t_pl (i,kp,j) = missing
            rh_pl (i,kp,j) = missing
            ght_pl(i,kp,j) = missing
            s_pl (i,kp,j) = missing
            td_pl (i,kp,j) = missing
         END DO
      END DO
   END DO
   j_loop : DO j = jts , MIN(jte,jde-1)
      i_loop : DO i = its , MIN(ite,ide-1)
         ke_h = kts
         ke_f = kts
         kp_loop : DO kp = 1 , num_press_levels
            ke_loop_half : DO ke = ke_h , kte-2
               IF ( use_tot_or_hyd_p .EQ. 1 ) THEN
                  pu = pp(i,ke+1,j)+pb(i,ke+1,j)
                  pd = pp(i,ke ,j)+pb(i,ke ,j)
               ELSE IF ( use_tot_or_hyd_p .EQ. 2 ) THEN
                  pu = p(i,ke+1,j)
                  pd = p(i,ke ,j)
               END IF
               pm = p_pl(kp)
               IF ( ( pd .GE. pm ) .AND. &
                    ( pu .LT. pm ) ) THEN
                  tu = (t(i,ke+1,j)+t0)*(pu/p1000mb)**rcp
                  td = (t(i,ke ,j)+t0)*(pd/p1000mb)**rcp
                  t_pl(i,kp,j) = ( tu * (pm-pd) + td * (pu-pm) ) / (pu-pd)
                  su = 0.5 * SQRT ( ( u(i,ke+1,j)+u(i+1,ke+1,j) )**2 + ( v(i,ke+1,j)+v(i,ke+1,j+1) )**2 )
                  sd = 0.5 * SQRT ( ( u(i,ke ,j)+u(i+1,ke ,j) )**2 + ( v(i,ke ,j)+v(i,ke ,j+1) )**2 )
                  s_pl(i,kp,j) = ( su * (pm-pd) + sd * (pu-pm) ) / (pu-pd)
                  uu = 0.5 * ( u(i,ke+1,j)+u(i+1,ke+1,j) )
                  ud = 0.5 * ( u(i,ke ,j)+u(i+1,ke ,j) )
                  u_pl(i,kp,j) = ( uu * (pm-pd) + ud * (pu-pm) ) / (pu-pd)
                  vu = 0.5 * ( v(i,ke+1,j)+v(i,ke+1,j+1) )
                  vd = 0.5 * ( v(i,ke ,j)+v(i,ke ,j+1) )
                  v_pl(i,kp,j) = ( vu * (pm-pd) + vd * (pu-pm) ) / (pu-pd)
                  qu = MAX(qv(i,ke+1,j),0.)
                  qd = MAX(qv(i,ke ,j),0.)
                  eu = qu * pu * 0.01 / ( eps + qu )
                  ed = qd * pd * 0.01 / ( eps + qd )
                  eu = max(eu, 0.001)
                  ed = max(ed, 0.001)
                  du = t_kelvin + ( s1 / ((s2 / log(eu/s3)) - 1.0) )
                  dd = t_kelvin + ( s1 / ((s2 / log(ed/s3)) - 1.0) )
                  td_pl(i,kp,j) = ( du * (pm-pd) + dd * (pu-pm) ) / (pu-pd)
                  qm = ( qu * (pm-pd) + qd * (pu-pm) ) / (pu-pd)
                  es = s4 * exp(s5 * (1.0 / 273.0 - 1.0 / t_pl(i,kp,j)))
                  qs = eps * es / (pm - es)
                  rh_pl(i,kp,j) = qm / qs * 100.
                  ke_h = ke
                  EXIT ke_loop_half
               END IF
            END DO ke_loop_half
            ke_loop_full : DO ke = ke_f , kte-1
               IF ( ( pw(i,ke ,j) .GE. p_pl(kp) ) .AND. &
                    ( pw(i,ke+1,j) .LT. p_pl(kp) ) ) THEN
                  pu = LOG(pw(i,ke+1,j))
                  pm = LOG(p_pl(kp))
                  pd = LOG(pw(i,ke ,j))
                  zu = ( zp(i,ke+1,j)+zb(i,ke+1,j) ) / g
                  zd = ( zp(i,ke ,j)+zb(i,ke ,j) ) / g
                  ght_pl(i,kp,j) = ( zu * (pm-pd) + zd * (pu-pm) ) / (pu-pd)
                  ke_f = ke
                  EXIT ke_loop_full
               END IF
            END DO ke_loop_full
         END DO kp_loop
      END DO i_loop
   END DO j_loop
END SUBROUTINE pld
END MODULE module_diagnostics
