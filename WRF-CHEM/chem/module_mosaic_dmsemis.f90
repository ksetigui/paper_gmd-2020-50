MODULE module_mosaic_dmsemis
CONTAINS
   subroutine mosaic_dmsemis( id, dtstep, dz8w, config_flags, &
               u_phy, v_phy, rho_phy, chem, emis_ant, alt, &
               u10, v10, xland, tsk, &
               ids,ide, jds,jde, kds,kde, &
               ims,ime, jms,jme, kms,kme, &
               its,ite, jts,jte, kts,kte )
  USE module_configure
  USE module_state_description
  USE module_data_radm2
  IMPLICIT NONE
   TYPE(grid_config_rec_type), INTENT(IN ) :: config_flags
   INTEGER, INTENT(IN ) :: id, &
                                  ids,ide, jds,jde, kds,kde, &
                                  ims,ime, jms,jme, kms,kme, &
                                  its,ite, jts,jte, kts,kte
   REAL, INTENT(IN ) :: dtstep
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ), &
         INTENT(INOUT ) :: chem
   REAL, DIMENSION( ims:ime, kms:config_flags%kemit, jms:jme,num_emis_ant),&
         INTENT(IN ) :: &
                         emis_ant
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ) , &
          INTENT(IN ) :: dz8w, rho_phy,alt, u_phy, v_phy
   REAL, DIMENSION( ims:ime , jms:jme ) , &
          INTENT(IN ) :: u10, v10, xland, tsk
   integer :: i,j,k
   real :: conv
   real :: w10, dms_emi, k_dms, dms_ocean_sfc = 0.0, sc_co2 = 600.0, sc_dms, &
                sst_ij_cels
   do j=jts,jte
      do i=its,ite
         conv = 4.828e-4/rho_phy(i,kts,j)*dtstep/(dz8w(i,kts,j)*60.)
             if (xland(i,j) .gt. 1.5) then
                 w10 = sqrt( u10(i,j)** 2.0 + v10(i,j)** 2.0 )
                 if(dz8w(i,kts,j).lt.12.)w10=sqrt(u_phy(i,kts,j)*u_phy(i,kts,j)+v_phy(i,kts,j)*v_phy(i,kts,j))
                 sst_ij_cels = ( tsk(i,j) - 273.15 )
                 if (sst_ij_cels < 5.0) sst_ij_cels = 5.0
                 if (sst_ij_cels > 30.0) sst_ij_cels = 30.0
                 sc_dms = 2674.0 - 147.12 * sst_ij_cels + 3.726 * sst_ij_cels ** 2.0 &
                             - 0.038 * sst_ij_cels ** 3
                 k_dms = ( 0.222 * w10**2 + 0.333 * w10) * (sc_dms / sc_co2)**(-0.5) / 3600. / 100.
                 dms_ocean_sfc=emis_ant(i,1,j,p_e_dms_oc)
                 dms_emi = k_dms * dms_ocean_sfc
                 dms_emi = dms_emi * 1.0E6 * 3600.0
                 chem(i,kts,j,p_dms) = chem(i,kts,j,p_dms) &
                             +dms_emi*conv
             end if
      end do
   end do
END subroutine mosaic_dmsemis
END MODULE module_mosaic_dmsemis
