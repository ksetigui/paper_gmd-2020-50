 MODULE module_ltng_cpmpr92z
 CONTAINS
 SUBROUTINE ltng_cpmpr92z ( &
                            dx, dy, xland, ht, z, t, &
                            kLNB, &
                            cldtop_adjustment, &
                            msft, &
                            cldfra_cup, &
                            htop, &
                            shall, &
                            ids, ide, jds, jde, kds, kde, &
                            ims, ime, jms, jme, kms, kme, &
                            ips, ipe, jps, jpe, kps, kpe, &
                            total_flashrate &
                          )
 USE module_state_description
 USE module_model_constants
 USE module_wrf_error
 IMPLICIT NONE
 REAL, INTENT(IN ) :: dx, dy
 REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN ) :: xland, ht
 REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN ) :: z, t
 INTEGER, DIMENSION( ims:ime, jms:jme ), INTENT(IN ) :: kLNB
 REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN ) :: cldfra_cup
 REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN ) :: htop
 REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN ) :: shall
 REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN ) :: msft
 REAL, INTENT(IN ) :: cldtop_adjustment
 INTEGER, INTENT(IN ) :: ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN ) :: ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN ) :: ips,ipe, jps,jpe, kps,kpe
 REAL, DIMENSION( ims:ime, jms:jme ), INTENT( OUT) :: total_flashrate
 REAL :: dA
 REAL :: zkm
 REAL :: cldfra_scaling_f
 REAL :: mapfac_scaling_f
 REAL, PARAMETER:: baseArea=1296.
 INTEGER :: i,k,j
 INTEGER :: kbtm, ktop
 CHARACTER (LEN=250) :: message
 dA = dx*dy/1E6
 cldfra_scaling_f=0.
 mapfac_scaling_f=1.
 total_flashrate( ips:ipe,jps:jpe ) = 0.
 jloop: DO j=jps,jpe
    iloop: DO i=ips,ipe
        mapfac_scaling_f = 1.0 / (msft(i,j) * msft(i,j))
        ktop = kLNB(i,j)
        IF (shall(i,j)>0.5) THEN
          ktop = -1
        ELSE
          ktop=nint(htop(i,j))
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
              ktop .ge. kps .and. ktop .le. kpe ) THEN
              zkm = ( z(i,ktop,j) - ht(i,j) )/1E3 + cldtop_adjustment
              IF ( zkm .gt. 0. ) THEN
                IF ( xland(i,j) .lt. 1.5 ) THEN
                  total_flashrate(i,j) = 3.44E-5 * (zkm**4.9) /60.
                ELSE
                  total_flashrate(i,j) = 6.57E-6 * (zkm**4.9) /60.
                ENDIF
                total_flashrate(i,j) = total_flashrate(i,j) * cldfra_scaling_f * mapfac_scaling_f
              ENDIF
            ENDIF
          ENDIF
        ENDIF
    ENDDO iloop
 ENDDO jloop
 total_flashrate(ips:ipe,jps:jpe) = total_flashrate(ips:ipe,jps:jpe) * dA/baseArea
 END SUBROUTINE ltng_cpmpr92z
 END MODULE module_ltng_cpmpr92z
