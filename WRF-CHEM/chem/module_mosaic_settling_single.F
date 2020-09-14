!*******************************************************************************
! Louis Marelle, June 2015: calculate aerosol sedimentation velocities and do 
! the sedimentation for MOSAIC aerosols (mass, number and aerosol water). This
! routine only performs sedimentation  above the first level. Deposition by
! sedimentation in the first level is treated in mixactivate (done in the
! released code). This module is  based on module_mosaic_drydep and 
! module_vash_settling
!
! For each i j, this routine calculates the column of deposition velocities for 
! each aerosol size, taking into account predicted aerosol wet radius (based on 
! water_* aerosols) and aerosol density. The method used is based on the Stokes 
! velocity + cunningham correction factor.
! The routine then updates aerosol mixing ratios at level k to account for the 
! loss of mass from settling in level k to level k-1 and increase of mass from 
! settling in level k+1 to level k (using an explicit method)
!
! For now, this routine only takes into accound the settling in the 
! interstitial (*_a01 etc.) phase. Activated aerosols do not settle here, 
! although cloud droplet settling might already be treated in another module.
!
! This routine should conserve mass (except in case of numerical errrors,
! roundings etc.)
!
!Things to look into are indicated by !lm flags
!lm no sedimentation for cw phase and hysteresis water ?
!lm this would be more efficient (one less k loop) if the settling velocities
!were calculated from top to bottom
!lm how is aerosol number removed in mixactivate, is this consistent ?
!lm I should add a mass check calculation (calculate column number and mass
!before and after the settling is done)
!**********************************************************************************  
    module module_mosaic_settling


    contains

!-----------------------------------------------------------------------
    subroutine mosaic_settling_driver(                                 &
        t_phy, rho_phy, dz8w, chem, dtstep,                            &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte                                      )

! Import mosaic data (aerosol density by chem composition, molecular weights,
! bin diameters, bin volumes, etc.). dens_* is in (g/cm3) d*_sect etc. are 
! in (cm), vol_* in cm3
! Also import pointers to number, water and aerosol arrays
    use module_data_mosaic_asect, only:  &
        dens_water_aer, mw_water_aer, dens_aer, & 
        dlo_sect, dcen_sect, dhi_sect, volumlo_sect, &
        volumhi_sect, ai_phase, ncomp_aer, nsize_aer, ntype_aer, &
        numptr_aer, waterptr_aer, massptr_aer
    use module_configure, only:  num_chem
    use module_data_mosaic_other, only: pi !lm make sure that this still works with the 'only' statement

    implicit none

!   subroutine arguments
    integer, intent(in) ::   &
        ids, ide, jds, jde, kds, kde,   &
        ims, ime, jms, jme, kms, kme,   &
        its, ite, jts, jte, kts, kte

    real, intent(in),   &
        dimension( ims:ime, kms:kme, jms:jme ) :: &
        t_phy, rho_phy, dz8w

    real, intent(inout),   &
        dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
        chem

    real, intent(in) :: dtstep

!   local variables
    integer it, jt, kt, idiagaa, ntdt, ndt_settl, nsettl
    integer iphase, itype, n, ll, l1 
    real temp
    real airdens
    real dumfact, dummass, rnum
    real dtmax, aer_mass_above, aer_mass, aer_water_above, aer_water, aer_num
    real dt_settl, air_volume_vert_ratio
    real aer_mass_change, aer_vol_change, aer_num_change
    real wetdp, wetdens, wetmass, wetvol, &
         drydp, drydens, drymass, dryvol 
    real vsettl
    real,  dimension(kts:kte) :: vsettl_col, dz_ij
!   diameters are in (cm), volumes are (cm^3/kg air) masses are (g/kg air)
!   parameters
    real, parameter :: densdefault = 2.0 ! g/cm3
!    real, parameter :: smallmassaa = 1.0e-20 !g/mole air
!    real, parameter :: smallmassbb = 1.0e-30
    real, parameter :: smallmassaa = 1.0e-19 ! g/kg air
    real, parameter :: smallmassbb = 1.0e-29
    real, parameter :: piover6 = pi/6.0
    real, parameter :: onethird = 1.0/3.0

! ************************ Variables description *******************************
! *** Subr arguments ***
!    t_phy: 3D air temperature (K)
!    rho_phy: 3D air density (kg m-3)
!    dz8w: 3D grid height (m)
!    chem: advected chemical/aer species (ppm or µg kg-1)
!    dtstep: time step (s)
!
! *** Local variables ***
!    it, jt, kt, iphase, itype, n, ll, : loop indices for wrf i dimension, wrf j 
!    dimension, wrf k dimensuon, aerosol phase, aerosol type, aerosol size, 
!    aerosol chemical composition
! 
!    idiagaa : switch for diagnostic output
!
!    temp : air temperature (K)
!    airdens: air density (g cm-3)
!
!    ndt_settl : number of sub time loops with time step dt_settl 
!                to use if the model time step is too high to respect the CFL condition
!    dt_settl : time step used in this routine to comply to the CFL condition.
!               In most cases, this will be equal to model dtstep (s)
! 
!    vsettl : settling velocity (cm s-1) at i,j,k,n
!    vesttl_col : settling velocity column (cm s-1) at i,j,n
!    dz_ij : height of vertical levels at i,j (cm ; column) 
!
!    aer_mass_above, aer_mass: aerosol mass mixing ratios (µg kg-1) at (i,j,k+1) and (i,j,k)
!    aer_water_above, aer_water: aerosol water mass mixing ratios (µg kg-1) at
!                                (i,j,k+1) and (i,j,k)
!    aer_num: aerosol number concentration (#/kg) at (i,j,k)
!    aer_mass_change, aer_vol_change, aer_num_change : mass, volume, number
!                     change in grid i,j,k from sedimentation
!    air_volume_vert_ratio : ratio [air_mass in grid (i,j,k+1)]/[air mass in grid i,j,k]
!  
!    rnum : aerosol number (# kg-1 air)
!    wetdp, wetdens, wetmass, wetvol: aerosol wet diameter (cm), density (g cm-3),
!                          mass (g kg-1 air), volume (cm-3 kg-1 air)
!    drydp, drydens, drymass, dryvol: same for dry aerosol
!    ***************************************************************************

!   loop over i and j
    do jt = jts, jte
    do it = its, ite

!   Get aerosol wet density and wet diameter

!   *** for now, just calc vdep_aer for iphase = ai_phase !lm, check in mixactivate how cw aerosols are removed
    iphase = ai_phase

    do itype = 1, ntype_aer
    do n = 1, nsize_aer(itype)

    do kt = kts, kte

      dryvol = 0.0
      drymass = 0.0 
  
      do ll = 1, ncomp_aer(itype) 
! This loop is "summing" chemical compounds; e.g. for bin 1, drymass = bc_a01+so4_a01+... etc and 
! dryvol = vol_bc_a01+... etc.  to compute bin diameters ((dryvol/(piover6*rnum))**onethird) 
          l1 = massptr_aer(ll,n,itype,iphase)
          dummass = chem(it,kt,jt,l1) * 1.0E-6 ! chem is in ug/kg, convert to g/kg air !lm check the units, check that the values are not too low
          drymass = drymass + dummass
          dryvol = dryvol + dummass/dens_aer(ll,itype) ! cm3/kg air
      end do
  
      l1 = waterptr_aer(n,itype)
      dummass = chem(it,kt,jt,l1) * 1.0E-6 ! chem is in ug/kg, convert to g/kg air
      wetmass = drymass + dummass
      wetvol = dryvol + dummass/dens_water_aer ! cm3/kg air, dens is in g/cm3
  
      l1 = numptr_aer(n,itype,iphase)
      rnum = chem(it,kt,jt,l1) ! #/kg
  
      if (drymass .le. smallmassbb) then
          drydp = dcen_sect(n,itype)
          drydens = densdefault
          wetdp = drydp
          wetdens = drydens
  
      else ! (drymass > smallmassbb)
          if (drymass .le. smallmassaa) then 
              wetmass = drymass
              wetvol = dryvol
          end if
          drydens = drymass/dryvol ! g/cm3
          wetdens = wetmass/wetvol ! g/cm3
  
          if (rnum .ge. dryvol/volumlo_sect(n,itype)) then 
          ! This is capping the dp diameters to [dlo - dhi],
          ! This routine is estimating one sedimentation velocity per bin
          ! (assuming internal mixing in the bin)
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
  
      endif !if (drymass .le. smallmassbb)
  
!   Get air temperature & density at it, kt, jt
      temp = t_phy(it,kt,jt)
      airdens = rho_phy(it,kt,jt)*1.0e-3 !kg/m3 -> g/cm3
  
!   Calculate settling velocity for mass
      call aerosol_depvel(   &
           wetdp, wetdens,temp, airdens, vsettl)
  
!lm Found this in the dry_dep driver, I'm not sure if this is needed here so I
!left it out for now 
! limit aerosol ddvels to <= 0.5 m/s
! drydep routines occasionally produce unrealistically-large
! particle diameter leading to unrealistically-large sedimentation
! velocity


!     Keep the vsettl column
      vsettl_col(kt)=vsettl

    end do !kt

! Sedimentation

!   vertical level depths in cm at i,j
    dz_ij(kts:kte) = dz8w(it,:,jt)*1.0E2 !m -> cm 

    do kt = kts, kte 

      ! Determine the maximum time-step satisying the CFL condition:
      ! dt <= (dz)_min / v_settl
      ntdt=INT(dtstep)
      dtmax = dz_ij(kt+1) / vsettl_col(kt+1) 
      dtmax = max( dtmax, dz_ij(kt) / vsettl_col(kt) ) 
      !lm It might be simpler to just use one dtmax for the whole column like in settl vash
      ndt_settl = MAX( 1, INT( ntdt /dtmax) )
      dt_settl = REAL(ntdt) / REAL(ndt_settl)
      ! limit maximum number of iterations. In this case, settling will be
      ! underestimated (this can only happen with very big time steps, big
      ! particles and very thin vertical levels) !lm add a warning, or call to
      ! wrf debug
      IF (ndt_settl > 12) ndt_settl = 12

      ! This ratio takes into account the change of air mass in grids
      ! between levels kt+1 and kt. If this is not taken into account, mass will
      ! not be conserved during settling (we need to transfer mass, not mixing
      ! ratios, between levels kt+1 and kt)
      air_volume_vert_ratio = ( rho_phy(it,kt+1,jt) * dz_ij(kt+1) ) / (rho_phy(it,kt,jt) * dz_ij(kt) )

      ! aerosol number at kt
      l1 = numptr_aer(n,itype,iphase)
      aer_num = chem(it,kt,jt,l1)

      ! loop on aerosol chemical composition (does not contain aerosol water)
      do ll = 1, ncomp_aer(itype)

        ! aerosol mass at kt and kt+1 for species ll
        l1 = massptr_aer(ll,n,itype,iphase)
        aer_mass = chem(it,kt,jt,l1)
        aer_mass_above = chem(it,kt+1,jt,l1)

!----------------------------------------------------------------------
! Do sedimentation for aerosol mass and number using an explicit scheme
!----------------------------------------------------------------------
!lm Delta values might be low enough to be approximated to 0 (or truncated) by
!   the floating point sum/substraction when doing the sedimentation. I'm not 
!   sure if this is an issue or not (I do not care if epsilon ~ 0 is approx. to 0).

!lm Maybe rewrite in smaller calculations or add parenthesis to limit numerical
!   approximations. Another possibility would be to do some calculations in
!   double precision
        do nsettl = 1, ndt_settl
           if (kt==kts) then !Surface, no loss to the level below (deposition already done elsewhere)
             aer_mass =  aer_mass + aer_mass_above * ( vsettl_col(kt+1) * dt_settl/dz_ij(kt+1) ) * &
                         air_volume_vert_ratio 
             aer_mass_change = aer_mass_above * ( vsettl_col(kt+1) * dt_settl/dz_ij(kt+1) ) * & 
                               air_volume_vert_ratio * 1.0E-6 ! 1E-6 = µg/kg to g/kg
             aer_vol_change = aer_mass_change / dens_aer(ll,itype) ! cm3/kg
!            lm: I calculate the number change using the center diameter, to be consistent with the 
!            number concentration update in mosaic_addemiss
             aer_num_change = 6 * aer_vol_change * 1./pi * 1. / dcen_sect(n,itype) ** 3  ! #/kg 
             aer_num = aer_num + aer_num_change
           elseif (kt==kte) then !Top level, no gain from the level above
             aer_mass =  aer_mass * (1. - vsettl_col(kt)*dt_settl/dz_ij(kt)) 
             aer_mass_change = - aer_mass * vsettl_col(kt)*dt_settl/dz_ij(kt) * 1.0E-6 
             aer_vol_change = aer_mass_change / dens_aer(ll,itype)
             aer_num_change = 6 * aer_vol_change * 1./pi * 1. / dcen_sect(n,itype) ** 3
             aer_num = aer_num + aer_num_change
           else !Default case, aerosols enter from above and exit below
             aer_mass =  aer_mass * ( 1. - vsettl_col(kt)*dt_settl / dz_ij(kt) ) + &
                         aer_mass_above * ( vsettl_col(kt+1) * dt_settl / dz_ij(kt+1) ) * &
                         air_volume_vert_ratio
             aer_mass_change = (- aer_mass * vsettl_col(kt) * dt_settl / dz_ij(kt) + & 
                               aer_mass_above * ( vsettl_col(kt+1) * dt_settl / dz_ij(kt+1) ) * &
                               air_volume_vert_ratio) * 1.0E-6
             aer_vol_change = aer_mass_change / dens_aer(ll,itype)
             aer_num_change = 6 * aer_vol_change * 1./pi * 1. / dcen_sect(n,itype) ** 3 
             aer_num = aer_num + aer_num_change
           end if
             
!         Diagnostic output
          idiagaa = 0
          if (idiagaa>0) print 9310, it, jt, kt, ndt_settl, n, itype, iphase,   &
               dcen_sect(n,itype), vsettl_col(kt), dt_settl, aer_mass_change
               9310  format( 'aerdep', 4i4, 3i3, 1p, 4e10.2,   &
               2x, 0p )

        end do !dtsettl 

!       Update aerosol mass mixing ratios
        if (aer_mass < 0.0 ) aer_mass = 0.0
        l1 = massptr_aer(ll,n,itype,iphase) 
        chem(it,kt,jt,l1) = aer_mass
 
      end do !ll, aer chem comp

!     Update number concentrations
      if (aer_num < 0.0) aer_num = 0.0
      l1 = numptr_aer(n,itype,iphase)
      chem(it,kt,jt,l1) = aer_num

!     Now do sedimentation for aerosol water (using a subroutine for this and
!     the mass change might be cleaner) 
      l1 = waterptr_aer(n,itype)
      aer_water = chem(it,kt,jt,l1)
      aer_water_above = chem(it,kt+1,jt,l1)
      do nsettl = 1, ndt_settl
         if (kt==kts) then 
           aer_water =  aer_water + aer_water_above * ( vsettl_col(kt+1) * dt_settl/dz_ij(kt+1) ) * &
                        air_volume_vert_ratio
         elseif (kt==kte) then
           aer_water =  aer_water * (1. - vsettl_col(kt)*dt_settl/dz_ij(kt))
         else
           aer_water =  aer_water * ( 1. - vsettl_col(kt)*dt_settl / dz_ij(kt) )+ &
                        aer_water_above * ( vsettl_col(kt+1) * dt_settl /dz_ij(kt+1) ) * &
                        air_volume_vert_ratio
         end if
      end do !dtsettl
      !Update the aerosol water mixing ratio
      if (aer_water < 0.0) aer_water = 0.0
      l1 = waterptr_aer(n,itype)
      chem(it,kt,jt,l1)=aer_water

    end do !kt

    end do !n, size
    end do !itype, type

    end do !it
    end do !jt


    end subroutine mosaic_settling_driver




!------------------------------------------------------------------------
    subroutine aerosol_depvel(   &
              dgnum, aerodens,   &
              temp, airdens,            &
              vsettl )

! 2015/06/25, Louis Marelle: This routine comes from module_mosaic_drydep.F
!
!   Computes the gravitational settling velocity (stokes velocity + cunningham
!   correction factor)
!
!   input parameters
!    dgnum - geometric mean diameter for aerosol number (cm)
!    aerodens - aerosol density (dgnum and aerodens are for the
!        actual wet distribution)
!    temp - temperature (K)
!    airdens - air density (g/cm^3)
!   output parameters
!    vsettl - gravitational settling velocity (cm/s)
!    
    implicit none

!   subroutine arguments
    real, intent(in) :: temp, airdens
    real, intent(in) :: dgnum, aerodens
    real, intent(out) :: vsettl

!   local variables
!    airkinvisc - air kinematic viscosity (cm^2/s)
    real airkinvisc
!   air molecular freepath (cm)
    real freepath
    real vsettl_dgnum
!   cunningham factor (adimensional)
    real cunningham_fact
    real pi
    parameter (pi = 3.1415926536)
!   gravity = gravitational acceleration in cm/s^2
    real gravity
    parameter (gravity = 980.616)
    
!   Calculate the kinetic viscosity and molecular freepath of air
    airkinvisc = ( 1.8325e-4 * (416.16/(temp+120.0)) *   &
                      ((temp/296.16)**1.5) ) / airdens
    freepath = 7.39758e-4 * airkinvisc / sqrt(temp)

!   cunningham correction factor for the stokes velocity
    cunningham_fact = 1.0 + 2.0 * freepath/dgnum * ( 1.257+ 0.400 * exp( -0.55 * dgnum/freepath ) )

!   Calculate the sedimentation velocity
    vsettl_dgnum = (gravity*aerodens*dgnum*dgnum)/   &
              (18.*airkinvisc*airdens)

    vsettl = vsettl_dgnum * cunningham_fact

    return
    end subroutine aerosol_depvel

!------------------------------------------------------------------------
    end module module_mosaic_settling

