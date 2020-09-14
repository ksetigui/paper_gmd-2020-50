    module module_mosaic_settling
    contains
    subroutine mosaic_settling_driver( &
        t_phy, rho_phy, dz8w, chem, dtstep, &
        ids,ide, jds,jde, kds,kde, &
        ims,ime, jms,jme, kms,kme, &
        its,ite, jts,jte, kts,kte )
    USE module_data_mosaic_asect, ONLY: &
        dens_water_aer, mw_water_aer, dens_aer, &
        dlo_sect, dcen_sect, dhi_sect, volumlo_sect, &
        volumhi_sect, ai_phase, ncomp_aer, nsize_aer, ntype_aer, &
        numptr_aer, waterptr_aer, massptr_aer, lptr_so4_aer
    USE module_configure, ONLY: num_chem
    USE module_data_mosaic_other, ONLY: pi
    USE module_wrf_error
    implicit none
    integer, intent(in) :: &
        ids, ide, jds, jde, kds, kde, &
        ims, ime, jms, jme, kms, kme, &
        its, ite, jts, jte, kts, kte
    real, intent(in), &
        dimension( ims:ime, kms:kme, jms:jme ) :: &
        t_phy, rho_phy, dz8w
    real, intent(inout), &
        dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
        chem
    real, intent(in) :: dtstep
    integer it, jt, kt, idiagaa, ntdt, ndt_settl, nsettl
    integer iphase, itype, n, ll, l1
    integer, parameter :: r8=8
    real(r8) :: temp
    real(r8) :: airdens
    real(r8) :: dumfact, dummass, rnum
    real(r8) :: dtmax, aer_mass_above, aer_mass, aer_water_above, aer_water, aer_num
    real(r8) :: dt_settl, air_volume_vert_ratio
    real(r8) :: aer_mass_change, aer_vol_change, aer_num_change
    real(r8) :: wetdp, wetdens, wetmass, wetvol, &
         drydp, drydens, drymass, dryvol
    real(r8) :: vsettl
    real(r8), dimension(kts:kte) :: vsettl_col, dz_ij
    real(r8), parameter :: densdefault = 2.0
    real(r8), parameter :: smallmassaa = 1.0e-19
    real(r8), parameter :: smallmassbb = 1.0e-29
    real(r8), parameter :: piover6 = pi/6.0
    real(r8), parameter :: onethird = 1.0/3.0
    CHARACTER (LEN=250) :: message
    real(r8), dimension(kts:kte) :: check_mass_after
    real(r8), dimension(kts:kte) :: check_mass_before
    real(r8) :: aer_mass_before
    do jt = jts, jte
    do it = its, ite
    iphase = ai_phase
    do itype = 1, ntype_aer
    do n = 1, nsize_aer(itype)
    do kt = kts, kte
      dryvol = 0.0_r8
      drymass = 0.0_r8
      do ll = 1, ncomp_aer(itype)
          l1 = massptr_aer(ll,n,itype,iphase)
          dummass = chem(it,kt,jt,l1) * 1.0E-6_r8
          drymass = drymass + dummass
          dryvol = dryvol + dummass/dens_aer(ll,itype)
      end do
      l1 = waterptr_aer(n,itype)
      dummass = chem(it,kt,jt,l1) * 1.0E-6_r8
      wetmass = drymass + dummass
      wetvol = dryvol + dummass/dens_water_aer
      l1 = numptr_aer(n,itype,iphase)
      rnum = chem(it,kt,jt,l1)
      if (drymass .le. smallmassbb) then
          drydp = dcen_sect(n,itype)
          drydens = densdefault
          wetdp = drydp
          wetdens = drydens
      else
          if (drymass .le. smallmassaa) then
              wetmass = drymass
              wetvol = dryvol
          end if
          drydens = drymass/dryvol
          wetdens = wetmass/wetvol
          if (rnum .ge. dryvol/volumlo_sect(n,itype)) then
              drydp = dlo_sect(n,itype)
          else if (rnum .le. dryvol/volumhi_sect(n,itype)) then
              drydp = dhi_sect(n,itype)
          else
              drydp = (dryvol/(piover6*rnum))**onethird
          end if
          if(abs(wetvol).gt.(1000.*abs(dryvol))) then
              dumfact=10.0
          else
              dumfact=abs(wetvol/dryvol)**onethird
              dumfact=max(1.0,min(dumfact,10.0))
          endif
          wetdp = drydp*dumfact
      endif
      temp = t_phy(it,kt,jt)
      airdens = rho_phy(it,kt,jt)*1.0e-3_r8
      call aerosol_depvel( &
           wetdp, wetdens,temp, airdens, vsettl)
      vsettl_col(kt)=vsettl
    end do
    dz_ij(kts:kte) = dz8w(it,:,jt)*1.0E2_r8
    l1 = lptr_so4_aer(n,itype,iphase)
    check_mass_before(kts:kte)=(dz8w(it,kts:kte,jt)*chem(it,kts:kte,jt,l1)*rho_phy(it,kts:kte,jt))
    do kt = kts, kte
      ntdt=INT(dtstep)
      dtmax = dz_ij(kt+1) / vsettl_col(kt+1)
      dtmax = max( dtmax, dz_ij(kt) / vsettl_col(kt) )
      ndt_settl = MAX( 1, INT( ntdt /dtmax) )
      dt_settl = REAL(ntdt) / REAL(ndt_settl)
      IF (ndt_settl > 12) ndt_settl = 12
      air_volume_vert_ratio = ( rho_phy(it,kt+1,jt) * dz_ij(kt+1) ) / (rho_phy(it,kt,jt) * dz_ij(kt) )
      l1 = numptr_aer(n,itype,iphase)
      aer_num = chem(it,kt,jt,l1)
      do ll = 1, ncomp_aer(itype)
        l1 = massptr_aer(ll,n,itype,iphase)
        aer_mass = chem(it,kt,jt,l1)
        aer_mass_above = chem(it,kt+1,jt,l1)
        aer_mass_before = aer_mass
        do nsettl = 1, ndt_settl
           if (kt .eq. kts) then
             aer_mass = aer_mass + aer_mass_above * ( vsettl_col(kt+1) * dt_settl/dz_ij(kt+1) ) * &
                         air_volume_vert_ratio
             aer_mass_change = aer_mass_above * ( vsettl_col(kt+1) * dt_settl/dz_ij(kt+1) ) * &
                               air_volume_vert_ratio * 1.0E-6_r8
             aer_vol_change = aer_mass_change / dens_aer(ll,itype)
             aer_num_change = 6.0_r8 * aer_vol_change * 1.0_r8/pi * 1.0_r8 / dcen_sect(n,itype) ** 3.0_r8
             aer_num = aer_num + aer_num_change
           elseif (kt .eq. kte) then
             aer_mass = aer_mass * (1.0_r8 - vsettl_col(kt)*dt_settl/dz_ij(kt))
             aer_mass_change = - aer_mass * vsettl_col(kt)*dt_settl/dz_ij(kt) * 1.0E-6_r8
             aer_vol_change = aer_mass_change / dens_aer(ll,itype)
             aer_num_change = 6.0_r8 * aer_vol_change * 1.0_r8/pi * 1.0_r8 / dcen_sect(n,itype) ** 3.0_r8
             aer_num = aer_num + aer_num_change
           else
             aer_mass = aer_mass * ( 1.0_r8 - vsettl_col(kt)*dt_settl / dz_ij(kt) ) + &
                         aer_mass_above * ( vsettl_col(kt+1) * dt_settl / dz_ij(kt+1) ) * &
                         air_volume_vert_ratio
             aer_mass_change = (- aer_mass * vsettl_col(kt) * dt_settl / dz_ij(kt) + &
                               aer_mass_above * ( vsettl_col(kt+1) * dt_settl / dz_ij(kt+1) ) * &
                               air_volume_vert_ratio) * 1.0E-6_r8
             aer_vol_change = aer_mass_change / dens_aer(ll,itype)
             aer_num_change = 6.0_r8 * aer_vol_change * 1.0_r8/pi * 1.0_r8 / dcen_sect(n,itype) ** 3.0_r8
             aer_num = aer_num + aer_num_change
           end if
        end do
          idiagaa = 0
          if (idiagaa>0) print 9310, it, jt, kt, n, ll, &
               vsettl_col(kt), dt_settl, dz_ij(kt), aer_mass, aer_mass_above, aer_mass_change, aer_mass-aer_mass_before
               9310 format( 'aerdep', 4i4, i3, 1p, 7e10.2, &
               2x, 0p )
        if (aer_mass < 0.0_r8 ) aer_mass = 0.0_r8
        l1 = massptr_aer(ll,n,itype,iphase)
        chem(it,kt,jt,l1) = aer_mass
      end do
      if (aer_num < 0.0_r8) aer_num = 0.0_r8
      l1 = numptr_aer(n,itype,iphase)
      chem(it,kt,jt,l1) = aer_num
      l1 = waterptr_aer(n,itype)
      aer_water = chem(it,kt,jt,l1)
      aer_water_above = chem(it,kt+1,jt,l1)
      do nsettl = 1, ndt_settl
         if (kt==kts) then
           aer_water = aer_water + aer_water_above * ( vsettl_col(kt+1) * dt_settl/dz_ij(kt+1) ) * &
                        air_volume_vert_ratio
         elseif (kt==kte) then
           aer_water = aer_water * (1.0_r8 - vsettl_col(kt)*dt_settl/dz_ij(kt))
         else
           aer_water = aer_water * ( 1.0_r8 - vsettl_col(kt)*dt_settl / dz_ij(kt) )+ &
                        aer_water_above * ( vsettl_col(kt+1) * dt_settl /dz_ij(kt+1) ) * &
                        air_volume_vert_ratio
         end if
      end do
      if (aer_water < 0.0_r8) aer_water = 0.0_r8
      l1 = waterptr_aer(n,itype)
      chem(it,kt,jt,l1)=aer_water
    end do
    l1 = lptr_so4_aer(n,itype,iphase)
    check_mass_after(kts:kte)=dz8w(it,kts:kte,jt)*chem(it,kts:kte,jt,l1)*rho_phy(it,kts:kte,jt)
    WRITE(message, * ) ' aer_settling it,jt,masscol_before,masscol_after = ',it, jt, sum(check_mass_before), sum(check_mass_after)
    IF (idiagaa >0) CALL wrf_debug ( 15, message )
    end do
    end do
    end do
    end do
    end subroutine mosaic_settling_driver
    subroutine aerosol_depvel( &
              dgnum, aerodens, &
              temp, airdens, &
              vsettl )
    implicit none
    integer, parameter :: r8=8
    real(r8), intent(in) :: temp, airdens
    real(r8), intent(in) :: dgnum, aerodens
    real(r8), intent(out) :: vsettl
    real airkinvisc
    real freepath
    real vsettl_dgnum
    real cunningham_fact
    real pi
    parameter (pi = 3.1415926536)
    real gravity
    parameter (gravity = 980.616)
    airkinvisc = ( 1.8325e-4 * (416.16/(temp+120.0)) * &
                      ((temp/296.16)**1.5) ) / airdens
    freepath = 7.39758e-4 * airkinvisc / sqrt(temp)
    cunningham_fact = 1.0 + 2.0 * freepath/dgnum * ( 1.257+ 0.400 * exp( -0.55 * dgnum/freepath ) )
    vsettl_dgnum = (gravity*aerodens*dgnum*dgnum)/ &
              (18.*airkinvisc*airdens)
    vsettl = vsettl_dgnum * cunningham_fact
    return
    end subroutine aerosol_depvel
    end module module_mosaic_settling
