MODULE GOCART_DUST_AFWA




  USE module_data_gocart_dust

  IMPLICIT NONE

  INTRINSIC max, min

CONTAINS
  subroutine gocart_dust_afwa_driver(ktau,dt,config_flags,julday,alt,t_phy,moist,u_phy, &
         v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,dustin,snowh,zs, &
         ivgtyp,isltyp,vegfra,xland,xlat,xlong,gsw,dx,g,emis_dust, &
         ust,znt,clay,sand,alpha,gamma,smtune, &
         ids,ide, jds,jde, kds,kde, &
         ims,ime, jms,jme, kms,kme, &
         its,ite, jts,jte, kts,kte )
  USE module_configure
  USE module_state_description

   TYPE(grid_config_rec_type), INTENT(IN ) :: config_flags

   INTEGER, INTENT(IN ) :: julday, ktau, &
                                  ids,ide, jds,jde, kds,kde, &
                                  ims,ime, jms,jme, kms,kme, &
                                  its,ite, jts,jte, kts,kte
   INTEGER,DIMENSION( ims:ime , jms:jme ) , &
          INTENT(IN ) :: &
                                                     ivgtyp, &
                                                     isltyp
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ), &
         INTENT(IN ) :: moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ), &
         INTENT(INOUT ) :: chem
   REAL, DIMENSION( ims:ime, 1, jms:jme,num_emis_dust),OPTIONAL, &
         INTENT(INOUT ) :: &
         emis_dust
   REAL, DIMENSION( ims:ime, config_flags%num_soil_layers, jms:jme ) , &
         INTENT(IN ) :: smois
   REAL, DIMENSION( config_flags%num_soil_layers ) , &
         INTENT(IN ) :: zs
   REAL, DIMENSION( ims:ime , jms:jme, ndcls ) , &
         INTENT(IN ) :: erod
   REAL, DIMENSION( ims:ime , jms:jme, 5 ) , &
         INTENT(INOUT) :: dustin
   REAL, DIMENSION( ims:ime , jms:jme ) , &
         INTENT(IN ) :: &
                                                     u10, &
                                                     v10, &
                                                     gsw, &
                                                     vegfra, &
                                                     xland, &
                                                     xlat, &
                                                     xlong, &
                                                     ust, &
                                                     znt, &
                                                     clay, &
                                                     sand, &
                                                     snowh
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ), &
         INTENT(IN ) :: &
                                                     alt, &
                                                     t_phy, &
                                                     dz8w,p8w, &
                                                     u_phy,v_phy,rho_phy
  REAL, INTENT(IN ) :: dt,dx,g



  INTEGER :: nmx,smx,i,j,k,imx,jmx,lmx,lhave
  INTEGER,DIMENSION (1,1) :: ilwi
  REAL*8, DIMENSION (1,1) :: erodtot
  REAL*8, DIMENSION (1,1) :: gravsm
  REAL*8, DIMENSION (1,1) :: drylimit
  REAL*8, DIMENSION (5) :: tc,bems
  REAL*8, DIMENSION (1,1) :: airden,airmas,ustar
  REAL*8, DIMENSION (1) :: dxy
  REAL*8, DIMENSION (3) :: massfrac
  REAL :: conver,converi,volsm
  REAL*8 :: zwant
  REAL, INTENT(IN ) :: alpha, gamma, smtune

  conver=1.e-9
  converi=1.e9



  imx=1
  jmx=1
  lmx=1
  nmx=ndust
  smx=ngsalt

  k=kts
  DO j=jts,jte
  DO i=its,ite



    IF (xland(i,j) .lt. 1.5) THEN
      ilwi(1,1)=1



      tc(1)=chem(i,kts,j,p_dust_1)*conver
      tc(2)=chem(i,kts,j,p_dust_2)*conver
      tc(3)=chem(i,kts,j,p_dust_3)*conver
      tc(4)=chem(i,kts,j,p_dust_4)*conver
      tc(5)=chem(i,kts,j,p_dust_5)*conver



      airmas(1,1)=-(p8w(i,kts+1,j)-p8w(i,kts,j))*dx*dx/g
      airden(1,1)=rho_phy(i,kts,j)
      ustar(1,1)=ust(i,j)
      dxy(1)=dx*dx



      erodtot(1,1)=SUM(erod(i,j,:))



      massfrac(1)=clay(i,j)
      massfrac(2)=1-(clay(i,j)+sand(i,j))
      massfrac(3)=sand(i,j)
      IF (znt(i,j) .gt. 0.2) then
        ilwi(1,1)=0
      ENDIF
      IF (isltyp(i,j) .eq. 15 .or. isltyp(i,j) .eq. 16. .or. &
          isltyp(i,j) .eq. 18) then
        ilwi(1,1)=0
      ENDIF
      IF (snowh(i,j) .gt. 0.01) then
        ilwi(1,1)=0
      ENDIF
      sfc_select: SELECT CASE(config_flags%sf_surface_physics)
        CASE (RUCLSMSCHEME)
          volsm=max((smois(i,1,j)+drypoint(isltyp(i,j)))*smtune,0.)
        CASE DEFAULT
          volsm=max(smois(i,1,j)*smtune,0.)
      END SELECT sfc_select
      gravsm(1,1)=100*volsm/((1.-porosity(isltyp(i,j)))*(2.65*(1-clay(i,j))+2.50*clay(i,j)))
      drylimit(1,1)=14.0*clay(i,j)*clay(i,j)+17.0*clay(i,j)
      call source_dust(imx, jmx, lmx, nmx, smx, dt, tc, ustar, massfrac, &
                       erodtot, ilwi, dxy, gravsm, airden, airmas, &
                       bems, g, drylimit, alpha, gamma)
      IF(config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11 ) then
        dustin(i,j,1:5)=tc(1:5)*converi
      ELSE
        chem(i,kts,j,p_dust_1)=tc(1)*converi
        chem(i,kts,j,p_dust_2)=tc(2)*converi
        chem(i,kts,j,p_dust_3)=tc(3)*converi
        chem(i,kts,j,p_dust_4)=tc(4)*converi
        chem(i,kts,j,p_dust_5)=tc(5)*converi
      ENDIF
      emis_dust(i,1,j,p_edust1)=bems(1)
      emis_dust(i,1,j,p_edust2)=bems(2)
      emis_dust(i,1,j,p_edust3)=bems(3)
      emis_dust(i,1,j,p_edust4)=bems(4)
      emis_dust(i,1,j,p_edust5)=bems(5)
    ENDIF
  ENDDO
  ENDDO
end subroutine gocart_dust_afwa_driver
  SUBROUTINE source_dust(imx, jmx, lmx, nmx, smx, dt1, tc, ustar, massfrac,&
                         erod, ilwi, dxy, gravsm, airden, airmas, &
                         bems, g0, drylimit, alpha, gamma)
  INTEGER, INTENT(IN) :: nmx,imx,jmx,lmx,smx
  INTEGER, INTENT(IN) :: ilwi(imx,jmx)
  REAL*8, INTENT(IN) :: erod(imx,jmx)
  REAL*8, INTENT(IN) :: ustar(imx,jmx)
  REAL*8, INTENT(IN) :: gravsm(imx,jmx)
  REAL*8, INTENT(IN) :: drylimit(imx,jmx)
  REAL*8, INTENT(IN) :: dxy(jmx)
  REAL*8, INTENT(IN) :: airden(imx,jmx,lmx), airmas(imx,jmx,lmx)
  REAL*8, INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)
  REAL*8, INTENT(OUT) :: bems(imx,jmx,nmx)
  REAL, INTENT(IN) :: g0,dt1
  REAL*8 :: den(smx), diam(smx)
  REAL*8 :: dvol(nmx), dlndp(nmx)
  REAL*8 :: dsurface(smx), ds_rel(smx)
  REAL*8 :: massfrac(3)
  REAL*8 :: u_ts0, u_ts, dsrc, srce, dmass, dvol_tot
  REAL*8 :: emit, emit_vol
  REAL :: rhoa, g
  REAL*8 :: salt, stotal
  INTEGER :: i, j, m, n, s
  REAL, INTENT(IN) :: alpha
  REAL, PARAMETER :: betamax=5.25E-4
  REAL*8 :: beta
  REAL, INTENT(IN) :: gamma
  REAL, PARAMETER :: cmb=1.0
  DO n=1,smx
    dmass=massfrac(spoint(n))*frac_salt(n)
    dsurface(n)=0.75*dmass/(den_salt(n)*reff_salt(n))
  ENDDO
  stotal=SUM(dsurface(:))
  DO n=1,smx
    ds_rel(n)=dsurface(n)/stotal
  ENDDO
 g = g0*1.0E2
 emit=0.0
 DO n = 1, smx
   den(n) = den_salt(n)*1.0D-3
   diam(n) = 2.0*reff_salt(n)*1.0D2
   DO i = 1,imx
     DO j = 1,jmx
       rhoa = airden(i,j,1)*1.0D-3
       u_ts0 = 0.13*1.0D-2*SQRT(den(n)*g*diam(n)/rhoa)* &
               SQRT(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
               SQRT(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0)
       IF (gravsm(i,j) > drylimit(i,j)) THEN
         u_ts = MAX(0.0D+0,u_ts0*(sqrt(1.0+1.21*(gravsm(i,j)-drylimit(i,j))**0.68)))
       ELSE
         u_ts = u_ts0
       END IF
       IF (ustar(i,j) .gt. u_ts .and. erod(i,j) .gt. 0.0 .and. ilwi(i,j) == 1) THEN
         salt = cmb*ds_rel(n)*(airden(i,j,1)/g0)*(ustar(i,j)**3)* &
                (1. + u_ts/ustar(i,j))*(1. - (u_ts**2)/(ustar(i,j)**2))
       ELSE
         salt = 0.D0
       ENDIF
       beta=10**(13.6*massfrac(1)-6.0)
       if (beta .gt. betamax) then
         beta=betamax
       endif
       emit=emit+salt*(erod(i,j)**gamma)*alpha*beta
     END DO
   END DO
 END DO
 DO n=1,nmx
   DO i=1,imx
     DO j=1,jmx
        dsrc = emit*distr_dust(n)*dxy(j)*dt1
        IF (dsrc < 0.0) dsrc = 0.0
        tc(i,j,1,n) = tc(i,j,1,n) + dsrc / airmas(i,j,1)
        bems(i,j,n) = 1000.*dsrc/(dxy(j)*dt1)
     END DO
   END DO
 END DO
END SUBROUTINE source_dust
END MODULE GOCART_DUST_AFWA
