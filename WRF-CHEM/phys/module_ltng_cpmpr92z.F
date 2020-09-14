! WRF:MODEL_LAYER:PHYSICS
!
! Lightning flash rate prediction based on cloud-top height. Implemented
! for models using convective parameterization. Assume use of sub-grid LNB.
!
! Price, C., and D. Rind (1992), A Simple Lightning Parameterization for Calculating
!   Global Lightning Distributions, J. Geophys. Res., 97(D9), 9919–9933, doi:10.1029/92JD00719.
!
! Wong, J., M. Barth, and D. Noone (2012), Evaluating a Lightning Parameterization
!   at Resolutions with Partially-Resolved Convection, GMDD, in preparation.
!
! Contact: J. Wong <johnwong@ucar.edu>
!
!**********************************************************************
!
!lm 2015/09: The predicted lightning rates are now scaled by cldfra_cup and by
!            the square of the map factor on the mass grid
!lm 2015/09: I rewrote this for KFCuP !! WARNING !! This is no longer working
!with Grell convection.


 MODULE module_ltng_cpmpr92z
 CONTAINS

 SUBROUTINE ltng_cpmpr92z ( &
                          ! Frequently used prognostics
                            dx, dy, xland, ht, z, t,              &
                          ! Scheme specific prognostics
                            kLNB,                                 &
                          ! Scheme specific namelist inputs
                            cldtop_adjustment,                    &
!lm
                            msft, &
                            cldfra_cup, &
                            htop, &
                            shall, &
!lm
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate                       &
                          )
!-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t

! Scheme specific prognostics
 INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: kLNB     ! model LNB from cu_physics
!lm
 REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: cldfra_cup
 REAL, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: htop
 REAL, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: shall
 REAL, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: msft
!lm

! Scheme specific namelist inputs
 REAL,    INTENT(IN   )    ::       cldtop_adjustment

! Order dependent args for domain, mem, and tile dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Local variables
 REAL :: dA              ! grid area dx*dy in km2
 REAL :: zkm             ! AGL z in km
!lm
 REAL :: cldfra_scaling_f    !average cloud fraction between cloud bottom and cloud top
 REAL :: mapfac_scaling_f
!lm

 REAL, PARAMETER:: baseArea=1296. ! base-case area, dx = 36 km

 INTEGER :: i,k,j
!lm
 INTEGER :: kbtm, ktop
 CHARACTER (LEN=250) :: message
!lm

!-----------------------------------------------------------------

 dA = dx*dy/1E6
!lm
 cldfra_scaling_f=0.
 mapfac_scaling_f=1.
!lm
 total_flashrate( ips:ipe,jps:jpe ) = 0.

! Compute AGL heights in km
 jloop: DO j=jps,jpe
    iloop: DO i=ips,ipe
        mapfac_scaling_f = 1.0 / (msft(i,j) * msft(i,j))
!lm This is super ugly, fix someday (this won't work with the Grell routines anymore)
        ktop = kLNB(i,j)
        IF (shall(i,j)>0.5) THEN !no deep convection
          ktop = -1
        ELSE !deep convection
          ktop=nint(htop(i,j))
          ! Look for kbtm
          k = kps
          DO WHILE( cldfra_cup(i,k,j) .le. 0.01 .and. k .le. ktop )
            k = k+1
          ENDDO
          kbtm = k
          IF ( (ktop - kbtm) > 0 ) THEN
            cldfra_scaling_f=SUM(cldfra_cup(i,kbtm:ktop,j))/REAL(ktop-kbtm+1)
            WRITE(message, * ) ' lightning_cpm: scaling by average cloud fraction (kbotom, ktop) = ', cldfra_scaling_f, kbtm, ktop
            CALL wrf_debug ( 15, message )
            IF ( t(i,ktop,j) .lt. 273.15 .and. &
              ktop .ge. kps .and. ktop .le. kpe ) THEN              ! Cloud top must freeze
              zkm = ( z(i,ktop,j) - ht(i,j) )/1E3 + cldtop_adjustment    ! Compute AGL height in km
              IF ( zkm .gt. 0. ) THEN                                         ! Must be above ground
                IF ( xland(i,j) .lt. 1.5 ) THEN
                  total_flashrate(i,j) = 3.44E-5 * (zkm**4.9) /60.            ! Continental equation
                ELSE
                  total_flashrate(i,j) = 6.57E-6 * (zkm**4.9) /60.            ! Marine equation (Michalon 99)
                ENDIF
!lm scale by cloud fraction and map factor
                total_flashrate(i,j) = total_flashrate(i,j) * cldfra_scaling_f * mapfac_scaling_f
!                total_flashrate(i,j) = total_flashrate(i,j) * cldfra_scaling_f
!lm
              ENDIF !AGL height is above ground
            ENDIF !cloud top is freezing
          ENDIF !cloud depth is > 0
        ENDIF !deep convection is on
!lm
    ENDDO iloop
 ENDDO jloop

! Scale by areal ratio
 total_flashrate(ips:ipe,jps:jpe) = total_flashrate(ips:ipe,jps:jpe) * dA/baseArea

 END SUBROUTINE ltng_cpmpr92z

 END MODULE module_ltng_cpmpr92z
