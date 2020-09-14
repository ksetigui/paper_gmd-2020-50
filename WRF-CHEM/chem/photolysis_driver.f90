




      SUBROUTINE photolysis_driver (id,curr_secs,ktau,dtstep, &
               config_flags,haveaer, &
               gmt,julday,t_phy,moist,aerwrf,p8w,t8w,p_phy, &
               chem,rho_phy,dz8w,xlat,xlong,z_at_w,gd_cloud,gd_cloud2, &
               ph_macr,ph_o31d,ph_o33p,ph_no2,ph_no3o2, &
               ph_no3o,ph_hno2,ph_hno3,ph_hno4,ph_h2o2, &
               ph_ch2or,ph_ch2om,ph_ch3cho,ph_ch3coch3, &
               ph_ch3coc2h5,ph_hcocho,ph_ch3cocho, &
               ph_hcochest,ph_ch3o2h,ph_ch3coo2h,ph_ch3ono2, &
               ph_hcochob, ph_n2o5,ph_o2,ph_n2o, &
               ph_pan,ph_mpan,ph_acetol,ph_gly, &
               ph_bigald,ph_mek,ph_c2h5ooh,ph_c3h7ooh,ph_pooh, &
               ph_rooh,ph_xooh,ph_isopooh,ph_alkooh, &
               ph_mekooh,ph_tolooh,ph_terpooh,ph_mvk, &
               ph_glyald,ph_hyac, &
               nref0, nw0, tuv_jmax0, &
               ph_radfld, ph_adjcoe, ph_prate, &
               wc, zref, &

               snowh, snowc, xice, lu_index, &

               tauaer1,tauaer2,tauaer3,tauaer4, &
               gaer1,gaer2,gaer3,gaer4, &
               waer1,waer2,waer3,waer4, &
               bscoef1,bscoef2,bscoef3,bscoef4, &
               l2aer,l3aer,l4aer,l5aer,l6aer,l7aer, &
               pm2_5_dry,pm2_5_water,uvrad,ivgtyp, &
               ids,ide, jds,jde, kds,kde, &
               ims,ime, jms,jme, kms,kme, &
               its,ite, jts,jte, kts,kte )


   USE module_configure
   USE module_state_description
   USE module_model_constants
   USE module_phot_mad
   USE module_phot_fastj
   USE module_ftuv_driver

   INTEGER, INTENT(IN ) :: id,julday, &
                                  ids,ide, jds,jde, kds,kde, &
                                  ims,ime, jms,jme, kms,kme, &
                                  its,ite, jts,jte, kts,kte
   INTEGER, INTENT(IN ) :: ktau
   REAL(KIND=8), INTENT(IN ) :: curr_secs
   REAL, INTENT(IN ) :: dtstep,gmt



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ), &
         INTENT(IN ) :: moist



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), &
         INTENT(INOUT ) :: &
               pm2_5_dry,pm2_5_water, aerwrf




   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), &
         INTENT(INOUT ) :: &
           ph_macr,ph_o31d,ph_o33p,ph_no2,ph_no3o2,ph_no3o,ph_hno2, &
           ph_hno3,ph_hno4,ph_h2o2,ph_ch2or,ph_ch2om,ph_ch3cho, &
           ph_ch3coch3,ph_ch3coc2h5,ph_hcocho,ph_ch3cocho, &
           ph_hcochest,ph_ch3o2h,ph_ch3coo2h,ph_ch3ono2,ph_hcochob, &
           ph_n2o5,ph_o2,ph_n2o,ph_pan,ph_mpan,ph_acetol,ph_gly, &
           ph_bigald,ph_mek,ph_c2h5ooh,ph_c3h7ooh,ph_pooh,ph_rooh, &
           ph_xooh,ph_isopooh,ph_alkooh,ph_mekooh,ph_tolooh, &
           ph_terpooh,ph_mvk,ph_glyald,ph_hyac

   INTEGER, INTENT(IN ) :: nref0, nw0, tuv_jmax0
   real, dimension( ims:ime, nref0, jms:jme, nw0 ), &
            intent(out ) :: ph_radfld
   real, dimension( ims:ime, nref0, jms:jme, tuv_jmax0 ), &
            intent(out ) :: ph_adjcoe
   real, dimension( ims:ime, nref0, jms:jme, tuv_jmax0 ), &
            intent(out ) :: ph_prate
   real, dimension(nw0), &
            intent(out ) :: wc
   real, dimension(nref0), &
            intent(out ) :: zref
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), &
         OPTIONAL, &
         INTENT(INOUT ) :: &
           gd_cloud,gd_cloud2



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), &
         INTENT(IN ) :: &
           tauaer1,tauaer2,tauaer3,tauaer4, &
           gaer1,gaer2,gaer3,gaer4, &
           waer1,waer2,waer3,waer4, &
           bscoef1,bscoef2,bscoef3,bscoef4
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, 1:4 ), &
         INTENT(IN ) :: &
           l2aer,l3aer,l4aer,l5aer,l6aer,l7aer



   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ), &
         INTENT(INOUT ) :: chem



   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ) , &
          INTENT(IN ) :: &
                                                      t_phy, &
                                                      p_phy, &
                                                      dz8w, &
                                              t8w,p8w,z_at_w , &
                                                    rho_phy
   REAL, DIMENSION( ims:ime , jms:jme ) , &
          INTENT(INOUT ) :: uvrad
   REAL, DIMENSION( ims:ime , jms:jme ) , &
          INTENT(IN ) :: &
                                                     xlat, &
                                                     xlong

   REAL, DIMENSION( ims:ime , jms:jme ) , &
          INTENT(IN ) :: snowh, snowc, xice, lu_index

   TYPE(grid_config_rec_type), INTENT(IN ) :: config_flags

   LOGICAL, INTENT(IN) :: haveaer
   integer, INTENT(IN ) :: ivgtyp(ims:ime,jms:jme)
   character(len=132) :: dbg_msg







   IF (config_flags%phot_opt .eq. 0) return

   write(dbg_msg,*) 'photolysis_driver: called for domain ',id
   CALL wrf_message( trim(dbg_msg) )



   chem_phot_select: SELECT CASE(config_flags%phot_opt)

     CASE (PHOTMAD)
       CALL wrf_debug(15,'calling madronich1_driver')
       call madronich1_driver(id,curr_secs,ktau,config_flags,haveaer, &
               gmt,julday,t_phy,moist,aerwrf,p8w,t8w,p_phy, &
               chem,rho_phy,dz8w,xlat,xlong,z_at_w,gd_cloud,gd_cloud2, &
               ph_macr,ph_o31d,ph_o33p,ph_no2,ph_no3o2,ph_no3o,ph_hno2,&
               ph_hno3,ph_hno4,ph_h2o2,ph_ch2or,ph_ch2om,ph_ch3cho, &
               ph_ch3coch3,ph_ch3coc2h5,ph_hcocho,ph_ch3cocho, &
               ph_hcochest,ph_ch3o2h,ph_ch3coo2h,ph_ch3ono2,ph_hcochob,&
               pm2_5_dry,pm2_5_water,uvrad, &
               ids,ide, jds,jde, kds,kde, &
               ims,ime, jms,jme, kms,kme, &
               its,ite, jts,jte, kts,kte )
     CASE (PHOTFASTJ)

       call wrf_debug(15,'calling fastj_driver')
       call fastj_driver(id,curr_secs,dtstep,config_flags, &
               gmt,julday,t_phy,moist,p8w,p_phy, &
               chem,rho_phy,dz8w,xlat,xlong,z_at_w, &
               ph_o2,ph_o31d,ph_o33p,ph_no2,ph_no3o2,ph_no3o,ph_hno2, &
               ph_hno3,ph_hno4,ph_h2o2,ph_ch2or,ph_ch2om,ph_ch3cho, &
               ph_ch3coch3,ph_ch3coc2h5,ph_hcocho,ph_ch3cocho, &
               ph_hcochest,ph_ch3o2h,ph_ch3coo2h,ph_ch3ono2,ph_hcochob,&
               ph_n2o5, &
               ivgtyp, &

               snowh, snowc, xice, lu_index, &

               tauaer1,tauaer2,tauaer3,tauaer4, &
               gaer1,gaer2,gaer3,gaer4, &
               waer1,waer2,waer3,waer4, &
               bscoef1,bscoef2,bscoef3,bscoef4, &
               l2aer,l3aer,l4aer,l5aer,l6aer,l7aer, &
               ids,ide, jds,jde, kds,kde, &
               ims,ime, jms,jme, kms,kme, &
               its,ite, jts,jte, kts,kte )
     CASE (FTUV)
       call wrf_debug(15,'calling ftuv_driver')
       call ftuv_driver( id, curr_secs, dtstep, config_flags, &
                              gmt, julday, &
                              p_phy, t_phy, rho_phy, p8w, t8w, &
                              xlat, xlong, z_at_w, &
                              moist, chem, gd_cloud,gd_cloud2, &
                              ph_no2,ph_o31d,ph_o33p,ph_hno2, &
                              ph_hno3,ph_hno4,ph_no3o2,ph_no3o, &
                              ph_h2o2,ph_ch2om,ph_ch2or,ph_ch3cho, &
                              ph_ch3o2h,ph_o2,ph_ch3coo2h,ph_ch3coch3,ph_hcocho, &
                              ph_hcochob,ph_ch3cocho,ph_hcochest,ph_ch3ono2, &
                              ph_macr,ph_ch3coc2h5,ph_n2o, &
                              ph_pan,ph_mpan,ph_acetol,ph_gly, &
                              ph_bigald,ph_mek,ph_c2h5ooh,ph_c3h7ooh, &
                              ph_pooh,ph_rooh,ph_xooh,ph_isopooh, &
                              ph_alkooh,ph_mekooh,ph_tolooh,ph_terpooh,&
                              ph_n2o5,ph_mvk,ph_glyald,ph_hyac, &
                              ivgtyp, &
                              ph_radfld, ph_adjcoe, ph_prate, &
                              wc,zref, &
                              tauaer1, tauaer2, tauaer3, tauaer4, &
                              waer1, waer2, waer3, waer4, &
                              gaer1, gaer2, gaer3, gaer4, &
                              ids,ide, jds,jde, kds,kde, &
                              ims,ime, jms,jme, kms,kme, &
                              its,ite, jts,jte, kts,kte )
     CASE DEFAULT

   END SELECT chem_phot_select

END SUBROUTINE photolysis_driver
