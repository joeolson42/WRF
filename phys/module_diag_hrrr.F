!WRF:MEDIATION_LAYER:PHYSICS
!
MODULE module_diag_hrrr
      
      USE module_model_constants
      PRIVATE :: WGAMMA
      PRIVATE :: GAMMLN
CONTAINS

   SUBROUTINE hrrr_refl_first_time_step(                               &    
                     ims,ime, jms,jme, kms,kme, kts,kte                &    
                    ,i_start,i_end,j_start,j_end,num_tiles             &    
                    ,mp_physics, th, pii, p, rho_curr, qv_curr         &    
                    ,qc_curr, qnr_curr, qr_curr, qs_curr, qg_curr      &    
                    ,refl_10cm                                         )    
!----------------------------------------------------------------------

  USE module_state_description, ONLY : THOMPSON,THOMPSONAERO,TEMPO
  USE module_mp_thompson, ONLY : calc_refl10cm

   IMPLICIT NONE 
!======================================================================             
!  dcd  subroutine to compute reflectivity and fall speed for wrfout at t=0
!       currently works for Thompson microphysics scheme only
!
! Definitions
!-----------
! Input:
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- kts
!-- kte
!-- i_start       start indices for i in tile
!-- i_end         end indices for i in tile
!-- j_start       start indices for j in tile
!-- j_end         end indices for j in tile
!-- num_tiles     number of tiles
!-- mp_physics    index of precipitation microphysics scheme
!-- th
!-- pii
!-- p
!-- rho_curr      air density (kg m-3)
!-- qv_curr       water vapor mixing ratio (kg kg-1)
!-- qc_curr       cloud water mixing ratio (kg kg-1)
!-- qnr_curr      rain number concentration (kg-1)
!-- qr_curr       rain water mixing ratio (kg kg-1)
!-- qs_curr       snow mixing ratio (kg kg-1)
!-- qg_curr       graupel mixing ratio (kg kg-1)
!
! Output:
!-- refl_10cm             radar reflectivity factor for 10-cm radar (dBZ)
!-- vt_dbz_wt             reflectivity-weighted precipitation fall speed (m s-1)
!======================================================================

   integer,      intent(in   )    ::                             &    
                             ims,ime, jms,jme, kms,kme, kts,kte, &
                                                      num_tiles

   integer, dimension(num_tiles), intent(in) ::                  &    
     &           i_start,i_end,j_start,j_end

   integer, intent(in) :: mp_physics

   real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) ::   &
                                qv_curr,qr_curr,qs_curr,qg_curr  &
                                              ,qc_curr,qnr_curr  &
                                             ,rho_curr,th,pii,p

   real, dimension( ims:ime, kms:kme, jms:jme ), intent(inout) :: refl_10cm
 !  real, dimension( ims:ime, kms:kme, jms:jme ), intent(out) :: vt_dbz_wt

   real    :: rand1
   integer :: i,j,k,ij,m,n

   real, dimension(kts:kte) :: qv1d, qc1d, qr1d, nr1d, qs1d, qg1d, &
            ng1d, qb1d, t1d, p1d, dbz, fallspd, rho
   logical :: qnr_initialized

   !thompson parameters:
   real :: ygra1,zans1,n0_exp,lam_exp,lamg
   real,    parameter:: r1    = 1.e-12
   real,    parameter:: r2    = 1.e-6
   real,    parameter:: pi    = 3.1415926536
!..densities of rain, snow, graupel, and cloud ice.
   real,    parameter:: rho_w = 1000.0
   real,    parameter:: rho_s = 100.0
   real,    parameter:: rho_i = 890.0
   integer, parameter:: nrhg  = 9
   integer, parameter:: nrhg1 = 1
   integer :: dimnrhg
   real, dimension(nrhg), parameter:: &
         rho_g = (/50., 100., 200., 300., 400., 500., 600., 700., 800./)
   integer, parameter :: idx_bg1 = 5 ! index for rhog when mp=8 or 28
   real, dimension(nrhg), parameter:: am_g = (/                       &
   &            pi*rho_g(1)/6.0, pi*rho_g(2)/6.0, pi*rho_g(3)/6.0,    &
   &            pi*rho_g(4)/6.0, pi*rho_g(5)/6.0, pi*rho_g(6)/6.0,    &
   &            pi*rho_g(7)/6.0, pi*rho_g(8)/6.0, pi*rho_g(9)/6.0 /)
   real, parameter:: bm_g = 3.0, mu_g = 0.0
   real, dimension(12,nrhg):: cge, cgg
   real:: oge1, ogg1, ogg2, ogg3, obmg
   real, dimension(nrhg):: & ! computed from a. heymsfield: best - reynolds relation
     &    bv_g = (/0.640961647, 0.640961647, 0.640961647, 0.640961647,  &
     &             0.640961647, 0.640961647, 0.640961647, 0.640961647,  &
     &             0.640961647/)


   if ( mp_physics .eq. thompson   .or. &
      mp_physics .eq. thompsonaero .or. &
      mp_physics .eq. tempo) then 

      cge(1,:) = bm_g + 1.
      cge(2,:) = mu_g + 1.
      cge(3,:) = bm_g + mu_g + 1.
      cge(4,:) = bm_g*2. + mu_g + 1.
      cge(10,:) = mu_g + 2.
      cge(12,:) = bm_g*0.5 + mu_g + 1.
      do m = 1, NRHG
         cge(5,m) = bm_g*2. + mu_g + bv_g(m) + 1.
         cge(6,m) = bm_g + mu_g + bv_g(m) + 1.
         cge(7,m) = bm_g*0.5 + mu_g + bv_g(m) + 1.
         cge(8,m) = mu_g + bv_g(m) + 1.                                  ! not used                                                                                                      
         cge(9,m) = mu_g + bv_g(m) + 3.
         cge(11,m) = 0.5*(bv_g(m) + 5. + 2.*mu_g)
      enddo
      do m = 1, NRHG
      do n = 1, 12
         cgg(n,m) = WGAMMA(cge(n,m))
      enddo
      enddo
      obmg = 1./bm_g
      oge1 = 1./cge(1,1)
      ogg1 = 1./cgg(1,1)
      ogg2 = 1./cgg(2,1)
      ogg3 = 1./cgg(3,1)


!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )

     DO ij = 1 , num_tiles

       qnr_initialized = .true.
       rand1 = 0.0

       DO j=j_start(ij),j_end(ij)
         DO i=i_start(ij),i_end(ij)
           DO k = kts, kte
             rho(k)=rho_curr(i,k,j)
             IF ( (qnr_curr(i,k,j) .le. 0.0) .and. (qr_curr(i,k,j) .ge. 1.e-12) ) THEN 
!              Rain number concentration was not initialized properly, but we would still
!              like a reflectivity estimate.  Let's arbitrarily assume the drops have a
!              Marshall-Palmer distribution with intercept parameter 8000000 m**-4.
!              The desired units of qnr are kg**-1.
               nr1d(k) =                                                                &    
                 (8.e6 / rho(k)) *                                                      &    
                 (rho(k) * qr_curr(i,k,j) / (3.14159 * 1000.0 * 8.e6) )**0.25
               qnr_initialized = .false.
             ELSE 
               nr1d(k) = qnr_curr(i,k,j)
             ENDIF

             t1d(k) = th(i,k,j)*pii(i,k,j)
             p1d(k) = p(i,k,j)
             qv1d(k) = qv_curr(i,k,j)
             qc1d(k) = qc_curr(i,k,j)
             qr1d(k) = qr_curr(i,k,j)
             qs1d(k) = qs_curr(i,k,j)
             qg1d(k) = qg_curr(i,k,j)
             if (qg1d(k).gt.R1) then
                ygra1 = alog10(max(1.E-9, qg1d(k)*rho(k)))
                zans1 = 3.0 + 2./7.*(ygra1+8.)
                zans1 = MAX(2., MIN(zans1, 6.))
                N0_exp = 10.**(zans1)
                lam_exp = (N0_exp*am_g(idx_bg1)*cgg(1,1)/(rho(k)*qg1d(k)))**oge1
                lamg = lam_exp * (cgg(3,1)*ogg2*ogg1)**obmg
                ng1d(k) = cgg(2,1)*ogg3*rho(k)*qg1d(k)*lamg**bm_g / am_g(idx_bg1)
                ng1d(k) = MAX(R2, ng1d(k)/rho(k))
                qb1d(k) = qg1d(k)/rho_g(idx_bg1)
             else
                ng1d(k) = 0
                qb1d(k) = 0
             endif 
           ENDDO

          call calc_refl10cm (qv1d, qc1d, qr1d, nr1d, qs1d, qg1d,       &
                      ng1d, qb1d, t1d, p1d, dBZ, kts, kte, i, j, kte)

           DO k = kts, kte
             refl_10cm(i,k,j) = MAX(-35., dbz(k))
             !vt_dbz_wt(i,k,j) = fallspd(k)
           ENDDO

         ENDDO
       ENDDO

       IF (.not. qnr_initialized) THEN 
         print *, 'QNR not initialized.  Marshall-Palmer distribution will be assumed for reflectivity calculation.'
       ENDIF

     ENDDO  ! ij 
!  !$OMP END PARALLEL DO

   ELSE 

     print*, 'skipping this procedure for mp_physics = ', mp_physics
     !CALL wrf_error_fatal('This microphysics scheme is not supported in refl_first_time_step.')

   ENDIF   ! IF ( mp_physics .EQ. THOMPSON )

   END SUBROUTINE hrrr_refl_first_time_step

   SUBROUTINE hrrr_diagnostic_driver(                                 &
                      ids,ide, jds,jde, kds,kde,                      &
                      ims,ime, jms,jme, kms,kme,                      &
                      ips,ipe, jps,jpe, kps,kpe,                      & ! patch  dims
                      i_start,i_end,j_start,j_end,kts,kte,num_tiles   &
                     ,hrrr_diagnostics                                &
                     ,dpsdt,dmudt                                     &
                     ,p8w,pk1m,mu_2,mu_2m                             &
                     ,u,v, temp                                       &
                     ,pblh                                            &
                     ,raincv,rainncv,rainc,rainnc                     &
                     ,i_rainc,i_rainnc                                & ! Optional - bucketr_opt =1  from bucket_mm
                     ,hfx,sfcevp,lh                                   &
                     ,ACSWUPT,ACSWUPTC,ACSWDNT,ACSWDNTC               & ! Optional - SW radiation package
                     ,ACSWUPB,ACSWUPBC,ACSWDNB,ACSWDNBC               & ! Optional - SW radiation package
                     ,ACLWUPT,ACLWUPTC,ACLWDNT,ACLWDNTC               & ! Optional - LW radiation package
                     ,ACLWUPB,ACLWUPBC,ACLWDNB,ACLWDNBC               & ! Optional - LW radiation package
                     ,I_ACSWUPT,I_ACSWUPTC,I_ACSWDNT,I_ACSWDNTC       & ! Optional - bucketf_opt =1  from bucket_J
                     ,I_ACSWUPB,I_ACSWUPBC,I_ACSWDNB,I_ACSWDNBC       & ! Optional - bucketf_opt =1  from bucket_J
                     ,I_ACLWUPT,I_ACLWUPTC,I_ACLWDNT,I_ACLWDNTC       & ! Optional - bucketf_opt =1  from bucket_J
                     ,I_ACLWUPB,I_ACLWUPBC,I_ACLWDNB,I_ACLWDNBC       & ! Optional - bucketf_opt =1  from bucket_J
                     ,dt,xtime,sbw,t2                                 &
                     ,diag_print                                      & ! Option to print out time-series
                     ,dfi_stage                                       &
                     ,bucket_mm, bucket_J                             &
                     ,mphysics_opt                                    &
                     ,gsfcgce_hail, gsfcgce_2ice                      &
                     ,mpuse_hail                                      &
                     ,nssl_cnoh, nssl_cnohl                           & !namelist parameters
                     ,nssl_rho_qh, nssl_rho_qhl                       & !namelist parameters
                     ,nssl_alphah, nssl_alphahl                       & !namelist parameters 
                     ,prec_acc_c, prec_acc_nc, snow_acc_nc            & ! Optional - prec_acc_opt =1
                     ,graupelncv                                      &
                     ,acsnow, acrunoff                                &
                     ,sfcrunoff, udrunoff                             &
                     ,snowncv, prec_acc_dt, curr_secs2                &
                     ,diagflag                                        &
                     ,history_interval                                &
!                     ,reset_interval1                                 &
                     ,itimestep                                       &
!                     ,ntimesteps                                      &
                     ,t,u10,v10,w,u_phy,v_phy,p,pb                    & !JOE
                     ,wspd10max                                       &
                     ,wspd10,wspd80,wspd160,LWP,IWP,SWP               &
                     ,maxcldfra,maxcldfra_pbl                         & !JOE: Optional - hrrr_diagnostics =1
                     ,cloud_regime                                    &
                     ,up_heli_max                                     & ! Optional - hrrr_diagnostics =1 
                     ,w_up_max,w_dn_max                               & ! Optional - hrrr_diagnostics =1
                     ,znw,w_colmean                                   &
                     ,w_mean                                          &
                     ,grpl_max,grpl_colint,refd_max,refl_10cm         &
                     ,msft                                            &
                     ,refl_10cm_1km,refl_10cm_4km,composite_refl_10cm & ! dcd
                     ,ht                                              & ! dcd 
                     ,hail_maxk1,hail_max2d                           & ! Optional - hrrr_diagnostics =1
                     ,qg_curr,qs_curr                                 &
                     ,ng_curr,qh_curr,nh_curr,qr_curr,nr_curr         & ! Optional (gthompsn)
                     ,qc_curr,qi_curr                                 & !JOE
                     ,rho,ph,phb,grav                                 &
                     ,cldfra,icloud_bl,cldfra_bl,qc_bl,qi_bl          & !JOE
                                                                      )
!----------------------------------------------------------------------

  USE module_dm, ONLY: wrf_dm_sum_real, wrf_dm_maxval
  USE module_state_description, ONLY :                                  &
      KESSLERSCHEME, LINSCHEME, SBU_YLINSCHEME, WSM3SCHEME, WSM5SCHEME, &
      WSM6SCHEME, ETAMPNEW, THOMPSON, THOMPSONAERO,                     &
      MORR_TWO_MOMENT, GSFCGCESCHEME, WDM5SCHEME, WDM6SCHEME,           &
      NSSL_2MOM, MILBRANDT2MOM , CAMMGMPSCHEME, TEMPO    !,MILBRANDT3MOM, NSSL_3MOM

   IMPLICIT NONE
!======================================================================
! Definitions
!-----------
!-- DIAG_PRINT    print control: 0 - no diagnostics; 1 - dmudt only; 2 - all
!-- DT            time step (second)
!-- XTIME         forecast time
!-- SBW           specified boundary width - used later
!
!-- P8W           3D pressure array at full eta levels
!-- MU            dry column hydrostatic pressure
!-- RAINC         cumulus scheme precipitation since hour 0
!-- RAINCV        cumulus scheme precipitation in one time step (mm)
!-- RAINNC        explicit scheme precipitation since hour 0
!-- RAINNCV       explicit scheme precipitation in one time step (mm)
!-- SNOWNCV       explicit scheme snow in one time step (mm)
!-- HFX           surface sensible heat flux
!-- LH            surface latent heat flux
!-- SFCEVP        total surface evaporation
!-- U             u component of wind - to be used later to compute k.e.
!-- V             v component of wind - to be used later to compute k.e.
!-- PREC_ACC_C    accumulated convective precip over accumulation time prec_acc_dt
!-- PREC_ACC_NC   accumulated explicit precip over accumulation time prec_acc_dt
!-- SNOW_ACC_NC   accumulated explicit snow precip over accumulation time prec_acc_dt
!-- GRAUP_ACC_NC  accumulated explicit graupel precip over accumulation time prec_acc_dt
!-- PREC_ACC_C1   accumulated convective precip over accumulation time prec_acc_dt1
!-- PREC_ACC_NC1  accumulated explicit precip over accumulation time prec_acc_dt1
!-- SNOW_ACC_NC1  accumulated explicit snow precip over accumulation time prec_acc_dt1
!-- GRAUP_ACC_NC1 accumulated explicit graupel precip over accumulation time prec_acc_dt1
!-- PREC_ACC_DT   precip accumulation time, default is 60 min
!-- PREC_ACC_DT1  precip accumulation time, default is 15 min
!-- ACSNOW        run-total accumulated snow
!-- ACGRAUP       run-total accumulated graupel
!-- ACFRAIN       run-total accumulated freezing rain (T<=273)
!-- CURR_SECS2    Time (s) since the beginning of the restart
!-- NWP_DIAGNOSTICS  = 1, compute hourly maximum fields *** JOE: removed - compute always
!-- DIAGFLAG      logical flag to indicate if this is a history output time
!-- U10, V10      10 m wind components
!-- WSPD10MAX     10 m max wind speed
!-- UP_HELI_MAX   max updraft helicity
!-- W_UP_MAX      max updraft vertical velocity
!-- W_DN_MAX      max downdraft vertical velocity
!-- W_COLMEAN     column mean vertical velocity
!-- NUMCOLPTS     no of column points
!-- GRPL_MAX      max column-integrated graupel
!-- GRPL_COLINT   column-integrated graupel
!-- REF_MAX       max derived radar reflectivity
!-- REFL_10CM     model computed 3D reflectivity
!-- REFL_10CM_1KM         reflectivity at 1 km AGL (dBZ)
!-- REFL_10CM_4KM         reflectivity at 4 km AGL (dBZ)
!-- COMPOSITE_REFL_10CM   maximum reflectivity in column (dBZ)
!-- HT            ht            terrain height (m)
!
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- ips           start index for i in patch
!-- ipe           end index for i in patch
!-- jps           start index for j in patch
!-- jpe           end index for j in patch
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- i_start       start indices for i in tile
!-- i_end         end indices for i in tile
!-- j_start       start indices for j in tile
!-- j_end         end indices for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!-- num_tiles     number of tiles
!
!======================================================================

   integer,   intent(in)    ::                                   &
                                      ids,ide, jds,jde, kds,kde, &
                                      ims,ime, jms,jme, kms,kme, &
                                      ips,ipe, jps,jpe, kps,kpe, &
                                                        kts,kte, &
                                                      num_tiles

   integer, dimension(num_tiles), intent(in) ::                  &
     &           i_start,i_end,j_start,j_end

   integer,   intent(in)    :: diag_print
   integer,   intent(in)    :: dfi_stage
   real,      intent(in)    :: bucket_mm, bucket_j
   integer,   intent(in)    :: mphysics_opt
   integer,   intent(in)    :: gsfcgce_hail, gsfcgce_2ice, mpuse_hail
   real,      intent(in)    ::    nssl_cnoh, nssl_cnohl          &
                               ,nssl_rho_qh, nssl_rho_qhl        &
                               ,nssl_alphah, nssl_alphahl
   integer,   intent(in)    ::   icloud_bl
   integer,   intent(in)    ::   hrrr_diagnostics
   logical,   intent(in)    ::   diagflag

   real, dimension( ims:ime, kms:kme, jms:jme ),                 &
         intent(in ) ::                                       u  &
                                                    ,         v  &
                                                    ,       p8w

   real, dimension( ims:ime , jms:jme ), intent(in) ::           &
                                                           mu_2  &
                                                    ,   rainncv  &
                                                    ,    raincv  &
                                                    ,   snowncv  &
                                                    ,graupelncv  &
                                                    ,       hfx  &
                                                    ,        lh  &
                                                    ,    sfcevp  &  
                                                    ,        t2  &
                                                    ,      pblh

   real, dimension( ims:ime, kms:kme, jms:jme ), optional,       &
         intent(in ) ::         cldfra, cldfra_bl, qc_bl, qi_bl     

   real, dimension( ims:ime , jms:jme ),                         &
          intent(inout) ::                                dpsdt  &
                                                    ,     dmudt  &
                                                    ,    rainnc  &
                                                    ,     rainc  &
                                                    ,     mu_2m  &
                                                    ,      pk1m
 
   real,  intent(in   ) :: dt, xtime
   integer,  intent(in   ) :: sbw
   integer, dimension( ims:ime , jms:jme ), intent(inout) ::     &
                                                       i_rainc,  &
                                                       i_rainnc

! these optional radiation variables should only be used when
! hrrr_diagnostics == 2
   real, dimension( ims:ime, jms:jme ), optional, intent(inout) ::&
                      acswupt,acswuptc,acswdnt,acswdntc,          &
                      acswupb,acswupbc,acswdnb,acswdnbc,          &
                      aclwupt,aclwuptc,aclwdnt,aclwdntc,          &
                      aclwupb,aclwupbc,aclwdnb,aclwdnbc
   integer, dimension( ims:ime, jms:jme ), optional, intent(inout) ::&
                      i_acswupt,i_acswuptc,i_acswdnt,i_acswdntc,  &
                      i_acswupb,i_acswupbc,i_acswdnb,i_acswdnbc,  &
                      i_aclwupt,i_aclwuptc,i_aclwdnt,i_aclwdntc,  &
                      i_aclwupb,i_aclwupbc,i_aclwdnb,i_aclwdnbc

   real, dimension( ims:ime, jms:jme ), optional, intent(inout) ::&
                      prec_acc_c, prec_acc_nc, snow_acc_nc        &
                                                   , acsnow       &
                                                   , acrunoff     &
                                                   , sfcrunoff    &
                                                   , udrunoff


   real, optional, intent(in)::  prec_acc_dt, curr_secs2

   integer           :: i,j,k,its,ite,jts,jte,ij
   integer           :: idp,jdp,irc,jrc,irnc,jrnc,isnh,jsnh
   integer           :: prfreq
   real, parameter   :: dbzmin=-10.0          ! dcd
   integer           :: k1, k2                ! dcd
   real              :: h_agl                 ! dcd
   real              :: cuprate               ! dcd
   real              :: ze, ze_conv, dbz_sum  ! dcd
   logical           :: found_1km, found_4km  ! dcd
   real              :: no_points
   real              :: dpsdt_sum, dmudt_sum, dardt_sum, drcdt_sum, drndt_sum
   real              :: hfx_sum, lh_sum, sfcevp_sum, rainc_sum, rainnc_sum, raint_sum
   real              :: dmumax, raincmax, rainncmax, snowhmax
   logical, external :: wrf_dm_on_monitor

   character*256     :: outstring
   character*6       :: grid_str

   integer, intent(in) ::                                        &
                     history_interval,itimestep

   real, dimension( kms:kme ), intent(in) ::                     &
                                                            znw

   real, dimension( ims:ime, kms:kme, jms:jme ), intent(in) ::   &
                                                              w  &
                                                          ,temp  &
                                               ,qs_curr,qg_curr  &
                                               ,qc_curr,qi_curr  &
                                                           ,rho  &
                                                        ,ph,phb  &
                                                          ,p,pb  &
                                                 ,t,u_phy,v_phy

   real, dimension( ims:ime, kms:kme, jms:jme ), intent(inout) :: & !dcd
                                                      refl_10cm

   real, dimension(ims:ime,kms:kme,jms:jme), optional, intent(in) ::    &
                                       ng_curr, qh_curr, nh_curr        &
                                               ,qr_curr, nr_curr

   real, dimension( ims:ime, jms:jme ), intent(in) ::            &
                                                            u10  &
                                                           ,v10  &
                                                           ,msft

   real, intent(in) :: grav

   real, dimension( ims:ime , jms:jme ), intent(inout) ::        &
                                                      wspd10max  &
                                                   ,up_heli_max  &
                                             ,w_up_max,w_dn_max  &
                                              ,w_colmean,w_mean  &
                                          ,grpl_max,grpl_colint  &
                                         ,hail_maxk1,hail_max2d  &
                                                      ,refd_max  &
               ,refl_10cm_1km,refl_10cm_4km,composite_refl_10cm  &  ! dcd
                                                            ,ht  &  ! dcd
                             ,wspd10,wspd80,wspd160,lwp,iwp,swp  &
                                       ,maxcldfra,maxcldfra_pbl  &
                                                 ,cloud_regime

   integer ::           ntimesteps
   real, dimension(ims:ime,kms:kme,jms:jme):: temp_qg, temp_ng,   &
                                              temp_qr, temp_nr

!   integer, dimension( ims:ime, jms:jme ), intent(inout) ::      &    
!                                        ltg1_calc

   real, dimension(ims:ime,jms:jme):: w_mean_sum

   integer :: idump

   real :: time_minutes
   real :: afrac
   real :: uh, vh, wind_vel
   real :: depth, zm
   real :: tval2,tval1
   real :: qr_val, qs_val, qg_val, rho_val
   real :: thr1mx, thr2mx, thr1, thr2
   real :: dbz,dbzr,dbzs,dbzg
   real, parameter :: tm10c=263.15

   integer :: ithr1mx,jthr1mx,ithr2mx,jthr2mx
   integer :: ithr1,jthr1,ithr2,jthr2

! derived composite radar reflectivity every 5 minutes

   real, parameter :: refddt=5.

! constants used in ltg threat calculations

   real, parameter    :: coef1=0.042*1000.*1.22
   real, parameter    :: coef2=0.20*1.22
   real, parameter    :: clim1=1.50
   real, parameter    :: clim2=0.40*1.22
   real, parameter    :: clim3=0.02*1.22
   real, parameter    :: capa=0.28589641e0

   double precision:: hail_max
   real:: hail_max_sp
   double precision, parameter:: thresh_conc = 0.0005d0                 ! number conc. of graupel/hail per cubic meter
   logical:: scheme_has_graupel
   integer, parameter:: ngbins=50
   double precision, dimension(ngbins+1):: xxdx
   double precision, dimension(ngbins):: xxdg, xdtg
   real:: xrho_g, xam_g, xbm_g, xmu_g
   real, dimension(3):: cge, cgg
   double precision:: f_d, sum_ng, sum_t, lamg, ilamg, n0_g, lam_exp, n0exp
   double precision:: lamr, n0min
   real:: mvd_r, xslw1, ygra1, zans1
   integer:: ng, n
   real :: dp,sum1,sum2,sum3,sum4,sum5,zm0,wt1,ter,u100,v100,t100, &
           qctotal,qitotal,qstotal

!-----------------------------------------------------------------
! Handle accumulations with buckets to prevent round-off truncation in long runs
! This is done every 360 minutes assuming time step fits exactly into 360 minutes
 IF (bucket_mm .gt. 0. .AND. MOD(NINT(XTIME),360) .EQ. 0)THEN
! SET START AND END POINTS FOR TILES
!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )

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
             i_rainnc(i,j) =  i_rainnc(i,j) + 1
          ENDIF
          IF(rainc(i,j) .gt. bucket_mm)THEN
             rainc(i,j) = rainc(i,j) - bucket_mm
             i_rainc(i,j) =  i_rainc(i,j) + 1
          ENDIF
       ENDDO      
       ENDDO

       !LIMIT THESE RADIATION DIAGS TO NON-REALTIME CONFIGURATION
       IF (hrrr_diagnostics .eq. 1) THEN
          !only set them == 0 for realtime application
          !and avoid using them again
          IF (xtime .eq. 0.0 .and. bucket_J .gt. 0.0 .and. PRESENT(ACSWUPT))THEN
             DO j=j_start(ij),j_end(ij)
             DO i=i_start(ij),i_end(ij)
                i_acswupt(i,j) = 1
                i_acswuptc(i,j) = 1
                i_acswdnt(i,j) = 1
                i_acswdntc(i,j) = 1
                i_acswupb(i,j) = 1
                i_acswupbc(i,j) = 1
                i_acswdnb(i,j) = 1
                i_acswdnbc(i,j) = 1
                acswupt(i,j) = 0.0
                acswuptc(i,j) = 0.0
                acswdnt(i,j) = 0.0
                acswdntc(i,j) = 0.0
                acswupb(i,j) = 0.0
                acswupbc(i,j) = 0.0
                acswdnb(i,j) = 0.0
                acswdnbc(i,j) = 0.0
             ENDDO
             ENDDO
          ENDIF
          IF (xtime .eq. 0.0  .and. bucket_J .gt. 0.0 .and. PRESENT(ACLWUPT))THEN
             DO j=j_start(ij),j_end(ij)
             DO i=i_start(ij),i_end(ij)
                i_aclwupt(i,j) = 1
                i_aclwuptc(i,j) = 1
                i_aclwdnt(i,j) = 1
                i_aclwdntc(i,j) = 1
                i_aclwupb(i,j) = 1
                i_aclwupbc(i,j) = 1
                i_aclwdnb(i,j) = 1
                i_aclwdnbc(i,j) = 1
                aclwupt(i,j) = 0.0
                aclwuptc(i,j) = 0.0
                aclwdnt(i,j) = 0.0
                aclwdntc(i,j) = 0.0
                aclwupb(i,j) = 0.0
                aclwupbc(i,j) = 0.0
                aclwdnb(i,j) = 0.0
                aclwdnbc(i,j) = 0.0
             ENDDO
             ENDDO
          ENDIF
       ELSEIF (hrrr_diagnostics .eq. 2) THEN
          !For research (hrrr_diagnostics=2), calculate these fields
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
          IF (xtime .eq. 0.0  .and. bucket_J .gt. 0.0 .and. PRESENT(ACLWUPT))THEN
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
                   i_acswupt(i,j) =  i_acswupt(i,j) + 1
                ENDIF
                IF(acswuptc(i,j) .gt. bucket_J)THEN
                   acswuptc(i,j) = acswuptc(i,j) - bucket_J
                   i_acswuptc(i,j) =  i_acswuptc(i,j) + 1
                ENDIF
                IF(acswdnt(i,j) .gt. bucket_J)THEN
                   acswdnt(i,j) = acswdnt(i,j) - bucket_J
                   i_acswdnt(i,j) =  i_acswdnt(i,j) + 1
                ENDIF
                IF(acswdntc(i,j) .gt. bucket_J)THEN
                   acswdntc(i,j) = acswdntc(i,j) - bucket_J
                   i_acswdntc(i,j) =  i_acswdntc(i,j) + 1
                ENDIF
                IF(acswupb(i,j) .gt. bucket_J)THEN
                   acswupb(i,j) = acswupb(i,j) - bucket_J
                   i_acswupb(i,j) =  i_acswupb(i,j) + 1
                ENDIF
                IF(acswupbc(i,j) .gt. bucket_J)THEN
                   acswupbc(i,j) = acswupbc(i,j) - bucket_J
                   i_acswupbc(i,j) =  i_acswupbc(i,j) + 1
                ENDIF
                IF(acswdnb(i,j) .gt. bucket_J)THEN
                   acswdnb(i,j) = acswdnb(i,j) - bucket_J
                   i_acswdnb(i,j) =  i_acswdnb(i,j) + 1
                ENDIF
                IF(acswdnbc(i,j) .gt. bucket_J)THEN
                   acswdnbc(i,j) = acswdnbc(i,j) - bucket_J
                   i_acswdnbc(i,j) =  i_acswdnbc(i,j) + 1
                ENDIF
             ENDDO      
             ENDDO
          ENDIF
          IF (PRESENT(ACLWUPT) .and. bucket_J .gt. 0.0)THEN
             DO j=j_start(ij),j_end(ij)
             DO i=i_start(ij),i_end(ij)
                IF(aclwupt(i,j) .gt. bucket_J)THEN
                aclwupt(i,j) = aclwupt(i,j) - bucket_J
                i_aclwupt(i,j) =  i_aclwupt(i,j) + 1
             ENDIF
             IF(aclwuptc(i,j) .gt. bucket_J)THEN
                aclwuptc(i,j) = aclwuptc(i,j) - bucket_J
                i_aclwuptc(i,j) =  i_aclwuptc(i,j) + 1
             ENDIF
             IF(aclwdnt(i,j) .gt. bucket_J)THEN
                aclwdnt(i,j) = aclwdnt(i,j) - bucket_J
                i_aclwdnt(i,j) =  i_aclwdnt(i,j) + 1
             ENDIF
             IF(aclwdntc(i,j) .gt. bucket_J)THEN
                aclwdntc(i,j) = aclwdntc(i,j) - bucket_J
                i_aclwdntc(i,j) =  i_aclwdntc(i,j) + 1
             ENDIF
             IF(aclwupb(i,j) .gt. bucket_J)THEN
                aclwupb(i,j) = aclwupb(i,j) - bucket_J
                i_aclwupb(i,j) =  i_aclwupb(i,j) + 1
             ENDIF
             IF(aclwupbc(i,j) .gt. bucket_J)THEN
                aclwupbc(i,j) = aclwupbc(i,j) - bucket_J
                i_aclwupbc(i,j) =  i_aclwupbc(i,j) + 1
             ENDIF
             IF(aclwdnb(i,j) .gt. bucket_J)THEN
                aclwdnb(i,j) = aclwdnb(i,j) - bucket_J
                i_aclwdnb(i,j) =  i_aclwdnb(i,j) + 1
             ENDIF
             IF(aclwdnbc(i,j) .gt. bucket_J)THEN
                aclwdnbc(i,j) = aclwdnbc(i,j) - bucket_J
                i_aclwdnbc(i,j) =  i_aclwdnbc(i,j) + 1
             ENDIF
          ENDDO      
          ENDDO
       ENDIF
    ENDIF !hrrr_diagnostic check
   ENDDO  !tile loop
!  !$OMP END PARALLEL DO
ENDIF !check:MOD(NINT(XTIME),360)

!
IF (prec_acc_dt .gt. 0.) THEN

!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )

   DO ij = 1 , num_tiles

      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
!         IF (mod(curr_secs2, 60.* prec_acc_dt) == 0.) THEN
         IF (curr_secs2 == 0 .or. mod(curr_secs2, 60.* prec_acc_dt) == 0.) THEN
!    print *,'zero out buckets', curr_secs2
            prec_acc_c(i,j)  = 0.
            prec_acc_nc(i,j) = 0.
            snow_acc_nc(i,j)  = 0.
 !           graup_acc_nc(i,j) = 0.
            sfcrunoff  (i,j) = 0.
            udrunoff   (i,j) = 0.
         ENDIF
         prec_acc_c(i,j)  = prec_acc_c(i,j)  +  RAINCV(i,j)
         prec_acc_nc(i,j) = prec_acc_nc(i,j) + RAINNCV(i,j)
         prec_acc_c(i,j)  = MAX (prec_acc_c(i,j), 0.0)
         prec_acc_nc(i,j) = MAX (prec_acc_nc(i,j), 0.0)
         snow_acc_nc(i,j)   = snow_acc_nc(i,j) + SNOWNCV(I,J)
         acsnow(i,j) = acsnow(i,j) + SNOWNCV(I,J)
!         graup_acc_nc(i,j)  = graup_acc_nc(i,j) + graupelncv(i,j)
!         acgraup(i,j) = acgraup(i,j) + graupelncv(i,j)
!         acfrain(i,j) = acfrain(i,j) + frain(i,j)
! add convective precip to snow bucket if t2 < 273.15
         IF ( t2(i,j) .lt. 273.15 ) THEN
            snow_acc_nc(i,j)   = snow_acc_nc(i,j) +  RAINCV(i,j)
            acsnow(i,j) = acsnow(i,j) + RAINCV(i,j)
         ENDIF
         snow_acc_nc(i,j)   = MAX (snow_acc_nc(i,j), 0.0)
!         graup_acc_nc(i,j)  = MAX (graup_acc_nc(i,j), 0.0)
         acsnow(i,j) = MAX (acsnow(i,j), 0.0)
!         acgraup(i,j)= MAX (acgraup(i,j), 0.0)

      ENDDO     
      ENDDO     

   ENDDO     

!  !$OMP END PARALLEL DO
 ENDIF



!JOE: removed   IF ( nwp_diagnostics .EQ. 1 ) THEN

 idump = (history_interval * 60.) / dt
     time_minutes = MOD(xtime,60.)

   !print *,' in HRRR diagnostics:'
   !print *,' diagflag         = ', diagflag
   !print *,' history_interval = ', history_interval
   !print *,' itimestep        = ', itimestep
   !print *,' idump            = ', idump
   !print *,' xtime            = ', xtime


! IF ( MOD(itimestep, idump) .eq. 0 ) THEN
!    WRITE(outstring,*) 'Computing PH0 for this domain with curr_secs2 = ', curr_secs2
!    CALL wrf_message ( TRIM(outstring) )

!   IF ( MOD((itimestep - 1), idump) .eq. 0 ) THEN
 IF ( MOD(curr_secs2,60.*history_interval) .eq. 0. ) THEN
    WRITE(outstring,*) 'Resetting max arrays for domain with dt = ', dt
    WRITE(outstring,*) 'Diagnostics: curr_sec2 = ', curr_secs2
    CALL wrf_debug ( 10,TRIM(outstring) )

    ntimesteps = 0 

!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )
    do ij = 1 , num_tiles
       do j=j_start(ij),j_end(ij)
       do i=i_start(ij),i_end(ij)
          wspd10max(i,j)    = 0. 
          up_heli_max(i,j)  = 0. 
          w_up_max(i,j)     = 0. 
          w_dn_max(i,j)     = 0. 
          w_mean_sum(i,j)   = 0. 
          w_mean(i,j)       = 0. 
          grpl_max(i,j)     = 0. 
          refd_max(i,j)     = 0. 
          hail_maxk1(i,j)   = 0. 
          hail_max2d(i,j)   = 0. 
          wspd10(i,j)       =-99.
          wspd80(i,j)       =-99.
          wspd160(i,j)      =-99.
          lwp(i,j)          = 0.
          iwp(i,j)          = 0.
          swp(i,j)          = 0.
          maxcldfra(i,j)    = 0.
          maxcldfra_pbl(i,j)= 0.
          cloud_regime(i,j) = 0.
       enddo
       enddo
    enddo
!  !$OMP END PARALLEL DO
 ENDIF

!   ntimesteps = ntimesteps + 1
!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )
 do ij = 1 , num_tiles
    do j=j_start(ij),j_end(ij)
    do i=i_start(ij),i_end(ij)
       ! zero some accounting arrays that will be used below
       w_colmean(i,j)   = 0.
       grpl_colint(i,j) = 0.
    enddo
    enddo
 enddo
!  !$OMP END PARALLEL DO

!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )
 do ij = 1 , num_tiles
    do j=j_start(ij),j_end(ij)
    do k=kms,kme
    do i=i_start(ij),i_end(ij)
       ! find vertical velocity max (up and down) below 100 mb
       if ( p8w(i,k,j) .gt. 10000. .and. w(i,k,j) .gt. w_up_max(i,j) ) then
          w_up_max(i,j) = w(i,k,j)
       endif

       if ( p8w(i,k,j) .gt. 10000. .and. w(i,k,j) .lt. w_dn_max(i,j) ) then
          w_dn_max(i,j) = w(i,k,j)
       endif

       ! for the column mean vertical velocity calculation, first
       ! total the vertical velocity between sigma levels 0.5 and 0.8
       if ( znw(k) .ge. 0.5 .and. znw(k) .le. 0.8 ) then
          w_colmean(i,j) = w_colmean(i,j) + w(i,k,j)
       endif
    enddo
    enddo
    enddo
 enddo
!  !$OMP END PARALLEL DO

!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )
 do ij = 1 , num_tiles
    do j=j_start(ij),j_end(ij)
    do k=kms,kme-1
    do i=i_start(ij),i_end(ij)
       ! calculate the column integrated graupel
       depth = ( ( ph(i,k+1,j) + phb(i,k+1,j) ) / grav ) - &
               ( ( ph(i,k  ,j) + phb(i,k  ,j) ) / grav )
       grpl_colint(i,j) = grpl_colint(i,j) + qg_curr(i,k,j) * depth * rho(i,k,j)
    enddo
    enddo
    enddo
 enddo
!  !$OMP END PARALLEL DO

!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )
 do ij = 1 , num_tiles
    do j=j_start(ij),j_end(ij)
    do i=i_start(ij),i_end(ij)
       ! calculate the max 10 m wind speed between output times
       wind_vel = sqrt ( u10(i,j)*u10(i,j) + v10(i,j)*v10(i,j) )
       if ( wind_vel .gt. wspd10max(i,j) ) then
         wspd10max(i,j) = wind_vel
       endif

       ! calculate the max column integrated graupel between output times
       if ( grpl_colint(i,j) .gt. grpl_max(i,j) ) then
          grpl_max(i,j) = grpl_colint(i,j)
       endif

       ! calculate the max radar reflectivity between output times
       !if ( refl_10cm(i,kms,j) .gt. refd_max(i,j) ) then
       !  refd_max(i,j) = refl_10cm(i,kms,j)
       !endif

       ! use k=kms for now, but can be reset to whatever level seems best
       dbzr = 0.
       dbzs = 0.
       dbzg = 0.
       if ( qr_curr(i,kms,j) .gt. 0. ) then
         dbzr = (( qr_curr(i,kms,j) * rho(i,kms,j) ) ** 1.75 )*     &
                3.630803e-9 * 1.e18
       endif

       if ( qs_curr(i,kms,j) .gt. 0. ) then
         dbzs = (( qs_curr(i,kms,j) * rho(i,kms,j) ) ** 1.75 )*     &
                2.18500e-10 * 1.e18
       endif

       if ( qg_curr(i,kms,j) .gt. 0. ) then
         dbzg = (( qg_curr(i,kms,j) * rho(i,kms,j) ) ** 1.75 )*     &
                1.033267e-9 * 1.e18
       endif

       dbz = dbzr + dbzs + dbzg
       if ( dbz  .gt. 0. ) dbz  = 10.0 * log10 ( dbz  )
       dbz  = max ( dbzmin, dbz  )

       if ( dbz .gt. refd_max(i,j) ) then
         refd_max(i,j) = dbz
       endif
    enddo
    enddo
 enddo
!  !$OMP END PARALLEL DO
! Calculate the derived radar reflectivity at the -10 deg C
! level between output times -- used in CI algorithm

!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )
 DO ij = 1 , num_tiles
    DO j=j_start(ij),j_end(ij)
    DO k=kme-1,kms+1,-1
    DO i=i_start(ij),i_end(ij)
       tval1 = ( t(i,k,j) + t0 ) *  &
               (( p(i,k,j) + pb(i,k,j) ) * 1.e-5 ) ** capa
       tval2 = ( t(i,k-1,j) + t0 ) *  &
               (( p(i,k-1,j) + pb(i,k-1,j) ) * 1.e-5 ) ** capa

!       if ( j.eq.j_start(ij) .and. i.eq.i_start(ij) )then
!         print *,' k, tval1, tval2 = ', k,tval1,tval2,refdm10c_calc(i,j)
!       endif

       IF ( tval1 .LE. tm10c .AND. tval2 .GE. tm10c ) THEN

         afrac = ( tm10c - tval2 ) / ( tval1 - tval2 )

         qr_val = qr_curr(i,k-1,j) + afrac * ( qr_curr(i,k,j) - qr_curr(i,k-1,j) )
         qs_val = qs_curr(i,k-1,j) + afrac * ( qs_curr(i,k,j) - qs_curr(i,k-1,j) )
         qg_val = qg_curr(i,k-1,j) + afrac * ( qg_curr(i,k,j) - qg_curr(i,k-1,j) )

         rho_val = rho(i,k-1,j) + afrac * ( rho(i,k,j) - rho(i,k-1,j) )

! Derived radar reflectivity between output times

! Use k=kms for now, but can be reset to whatever level seems best

         dbzr = 0.
         dbzs = 0.
         dbzg = 0.

         if ( qr_val .gt. 0. ) then
           dbzr = (( qr_val * rho_val ) ** 1.75 )*     &
                  3.630803e-9 * 1.e18
         endif

         if ( qs_val .gt. 0. ) then
           dbzs = (( qs_val * rho_val ) ** 1.75 )*     &
                  2.18500e-10 * 1.e18
         endif

         if ( qg_val .gt. 0. ) then
           dbzg = (( qg_val * rho_val ) ** 1.75 )*     &
                  1.033267e-9 * 1.e18
         endif

         dbz = dbzr + dbzs + dbzg
         if ( dbz  .gt. 0. ) dbz  = 10.0 * log10 ( dbz  )
         dbz = max ( dbzmin, dbz  )

       endif
    enddo
    enddo
    enddo
 enddo
!  !$OMP END PARALLEL DO

! dcd  reflectivity diagnostics for RAP and HRRR                                                                                                      
 IF (diagflag) THEN                                                                                                                                 
    k1 = kms                                                                                                                                         
    k2 = kme-1                                                                                                                                       
!  !$OMP PARALLEL DO   &                                                                                                                              
!  !$OMP PRIVATE ( ij )                                                                                                                               
    do ij = 1 , num_tiles                                                                                                                            
                                                                                                                                                      
       do j=j_start(ij),j_end(ij)                                                                                                                     
       do i=i_start(ij),i_end(ij)
          composite_refl_10cm(i,j) = dbzmin                                                                                                            
          refl_10cm_1km(i,j) = dbzmin                                                                                                                  
          refl_10cm_4km(i,j) = dbzmin                                                                                                                  
          found_1km = .false.                                                                                                                          
          found_4km = .false.                                                                                                                          
                                                                                                                                                      
          !reflectivity parameterization for parameterized convection (reference:  unipost mdlfld.f)                                                    
          ze_conv = 0.0                                                                                                                                
          cuprate = raincv(i,j) * 3600.0 / dt          ! cu precip rate (mm/h)                                                                         
          ze_conv = 300.0 * cuprate**1.4                                                                                                               
                                                                                                                                                      
          do k=k1,k2                                                                                                                                   
                                                                                                                                                      
             !approximation:  use base-state geopotential to determine height
             h_agl = phb(i,k,j) / grav - ht(i,j)
             ze = 10.0 ** (0.1 * refl_10cm(i,k,j))                                                                                                      
             dbz_sum = max(dbzmin, 10.0 * log10(ze + ze_conv))                                                                                          
             refl_10cm(i,k,j) = dbz_sum                                                                                                                 
                                                                                                                                                      
             if ( (.not. found_1km) .and. (h_agl.gt.1000.0) ) then                                                                                      
                refl_10cm_1km(i,j) = dbz_sum                                                                                                             
                found_1km = .true.                                                                                                                       
             endif                                                                                                                                      
             if ( (.not. found_4km) .and. (h_agl.gt.4000.0) ) then                                                                                      
                refl_10cm_4km(i,j) = dbz_sum                                                                                                             
                found_4km = .true.                                                                                                                       
             endif                                                                                                                                      
             if ( dbz_sum .gt. composite_refl_10cm(i,j) ) then                                                                                          
                composite_refl_10cm(i,j) = dbz_sum                                                                                                       
             endif                                                                                                                                      
                                                                                                                                                      
          enddo                                                                                                                                        
                                                                                                                                                      
       enddo                                                                                                                                          
       enddo                                                                                                                                          
                                                                                                                                                      
    enddo                                                                                                                                            
!  !$OMP END PARALLEL DO                                                                                                                              
 ENDIF                                                                                                                                              


!Diagnostics for PBL/Cloud Development
 IF (diagflag) THEN
!   print*," doing 10, 80 m winds & lwp, iwp, swp..."
!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )
   do ij = 1 , num_tiles

     do j=j_start(ij),j_end(ij)
     do i=i_start(ij),i_end(ij)

       ter = ( ph(i,1,j) + phb(i,1,j) ) / grav
       zm=0.
       sum1=0.
       sum2=0.
       sum3=0.
       sum4=0.
       sum5=0.

       do k=kms,kme-1
          !note: ph & phb are defined at the layer interface, so we average to
          !      get a height (amsl) at the middle of the model layers
          zm0=zm
          zm = 0.5*((ph(i,k,j) + phb(i,k,j) + ph(i,k+1,j) + phb(i,k+1,j)) / grav )
          depth = zm-zm0

          !10 m wind speed
          if (zm-ter > 10. .and. wspd10(i,j)==-99.) then
             if (k==kms) then
                wspd10(i,j)=sqrt(u_phy(i,k,j)**2 + v_phy(i,k,j)**2)
             else
                wt1=(zm-ter - 10.)/depth    ! portion of layer above 10 m
                wt1=min(max(0.,wt1),1.0)
                wspd10(i,j)=(1.-wt1)*sqrt(u_phy(i,k,j)**2 + v_phy(i,k,j)**2) + &
                             wt1*sqrt(u_phy(i,k-1,j)**2 + v_phy(i,k-1,j)**2)
                !print*,"z-ter=",zm-ter," wt1=",wt1
                !print*,i,j," wspd10=",wspd10(i,j)
             endif
          endif

          !80 m wind speed
          if (zm-ter > 80. .and. wspd80(i,j)==-99.) then
             if (k==kms) then
                wspd80(i,j)=sqrt(u_phy(i,k,j)**2 + v_phy(i,k,j)**2)
             else
                wt1=(zm-ter - 80.)/depth ! portion of layer above 80 m
                wt1=min(max(0.,wt1),1.0)
                wspd80(i,j)=(1.-wt1)*sqrt(u_phy(i,k,j)**2 + v_phy(i,k,j)**2) + &
                                wt1*sqrt(u_phy(i,k-1,j)**2 + v_phy(i,k-1,j)**2)
                !print*,"z-ter=",zm-ter," wt1=",wt1
                !print*,i,j," wspd80=",wspd80(i,j)
             endif
          endif

          !160 m wind speed
          if (zm-ter > 160. .and. wspd160(i,j)==-99.) then
             if (k==kms) then
                wspd160(i,j)=sqrt(u_phy(i,k,j)**2 + v_phy(i,k,j)**2)
             else
                wt1=(zm-ter - 160.)/depth ! portion of layer above 160 m
                wt1=min(max(0.,wt1),1.0)
                wspd160(i,j)=(1.-wt1)*sqrt(u_phy(i,k,j)**2 + v_phy(i,k,j)**2) + &
                                  wt1*sqrt(u_phy(i,k-1,j)**2 + v_phy(i,k-1,j)**2)
             endif
          endif
          
          !joe -add lwp, iwp & swp & maxcldfra; use temp (abs temperature in kelvin)
          dp=p8w(i,k,j)-p8w(i,k+1,j)
          if (icloud_bl>0 .and. present(qc_bl) .and. &
             present(cldfra_bl) .and. qc_curr(i,k,j)<1e-6 .and. &
             cldfra_bl(i,k,j)>0.01) then
             !already temp-dependent sorted, and now grid mean, so dont multiply by cf 
             qctotal = qc_bl(i,k,j)
          else
             qctotal = qc_curr(i,k,j)
          endif
          if (icloud_bl>0 .and. present(qi_bl) .and. &
             present(cldfra_bl) .and. qi_curr(i,k,j)<1e-9 .and. &
             cldfra_bl(i,k,j)>0.01) then
             !already temp-dependent sorted, and now grid mean, so dont multiply by cf
             qitotal = qi_bl(i,k,j)
          else
             qitotal = qi_curr(i,k,j)
          endif
          qstotal = qs_curr(i,k,j) !ther ecurrently is no sgs snow...
          sum1 = sum1 + max((dp/grav)*qctotal,0.0)
          sum2 = sum2 + max((dp/grav)*(qitotal+qstotal),0.0) !actually frozen water path
          sum3 = sum3 + max((dp/grav)*qstotal,0.0)
          sum4 = max(sum4,cldfra_bl(i,k,j))
          if (zm < (pblh(i,j)+1000.)) sum5 = max(sum5,cldfra_bl(i,k,j))
       enddo

       lwp(i,j)=sum1*1000. ! kg m-2  --> g m-2
       iwp(i,j)=sum2*1000.
       swp(i,j)=sum3*1000.
       maxcldfra(i,j)=sum4
       maxcldfra_pbl(i,j)=sum5
       
     enddo
     enddo

   enddo
!  !$OMP END PARALLEL DO
ENDIF !check for diagflag
!End extra diags for cloud/PBL development

!+---+-----------------------------------------------------------------+ 
!+---+-----------------------------------------------------------------+ 
!..Calculate a maximum hail diameter from the characteristics of the
!.. graupel category mixing ratio and number concentration (or hail, if
!.. available).  This diagnostic uses the actual spectral distribution
!.. assumptions, calculated by breaking the distribution into 50 bins
!.. from 0.5mm to 7.5cm.  Once a minimum number concentration of 0.01
!.. particle per cubic meter of air is reached, from the upper size
!.. limit, then this bin is considered the max size.
!+---+-----------------------------------------------------------------+ 

      WRITE(outstring,*) 'GT-Diagnostics, computing max-hail diameter'
      CALL wrf_debug (100, TRIM(outstring))


   IF (PRESENT(qh_curr)) THEN
   WRITE(outstring,*) 'GT-Debug, this mp scheme, ', mphysics_opt, ' has hail mixing ratio'
   CALL wrf_debug (150, TRIM(outstring))
!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )
   DO ij = 1 , num_tiles
     DO j=j_start(ij),j_end(ij)
     DO k=kms,kme-1
     DO i=i_start(ij),i_end(ij)
        temp_qg(i,k,j) = MAX(1.E-12, qh_curr(i,k,j)*rho(i,k,j))
     ENDDO
     ENDDO
     ENDDO
   ENDDO
!  !$OMP END PARALLEL DO
   ELSE
!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )
   DO ij = 1 , num_tiles
     DO j=j_start(ij),j_end(ij)
     DO k=kms,kme-1
     DO i=i_start(ij),i_end(ij)
        temp_qg(i,k,j) = MAX(1.E-12, qg_curr(i,k,j)*rho(i,k,j))
     ENDDO
     ENDDO
     ENDDO
   ENDDO
!  !$OMP END PARALLEL DO
   ENDIF

   IF (PRESENT(nh_curr)) THEN
   WRITE(outstring,*) 'GT-Debug, this mp scheme, ', mphysics_opt, ' has hail number concentration'
   CALL wrf_debug (150, TRIM(outstring))
!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )
   DO ij = 1 , num_tiles
     DO j=j_start(ij),j_end(ij)
     DO k=kms,kme-1
     DO i=i_start(ij),i_end(ij)
        temp_ng(i,k,j) = MAX(1.E-8, nh_curr(i,k,j)*rho(i,k,j))
     ENDDO
     ENDDO
     ENDDO
   ENDDO
!  !$OMP END PARALLEL DO
   ELSEIF (PRESENT(ng_curr)) THEN
      WRITE(outstring,*) 'This mp scheme, ', mphysics_opt, &
      ' has graupel num conc'
!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )
   DO ij = 1 , num_tiles
     DO j=j_start(ij),j_end(ij)
     DO k=kms,kme-1
     DO i=i_start(ij),i_end(ij)
        temp_ng(i,k,j) = MAX(1.E-8, ng_curr(i,k,j)*rho(i,k,j))
     ENDDO
     ENDDO
     ENDDO
   ENDDO
!  !$OMP END PARALLEL DO
   ELSE
!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )
   DO ij = 1 , num_tiles
     DO j=j_start(ij),j_end(ij)
     DO k=kms,kme-1
     DO i=i_start(ij),i_end(ij)
        temp_ng(i,k,j) = 1.E-8
     ENDDO
     ENDDO
     ENDDO
   ENDDO
!  !$OMP END PARALLEL DO
   ENDIF


      ! Calculate the max hail size from graupel/hail parameters in microphysics scheme (gthompsn 12Mar2015)
      ! First, we do know multiple schemes have assumed inverse-exponential distribution with constant
      ! intercept parameter and particle density.  Please leave next 4 settings alone for common
      ! use and copy these and cge constants to re-assign per scheme if needed (e.g. NSSL).

      xrho_g = 500.
      xam_g = 3.1415926536/6.0*xrho_g     ! Assumed m(D) = a*D**b, where a=PI/6*rho_g and b=3
      xbm_g = 3.                          ! in other words, spherical graupel/hail
      xmu_g = 0.
      scheme_has_graupel = .false.

      !..Some constants. These *MUST* get changed below per scheme
      !.. *IF* either xbm_g or xmu_g value is changed from 3 and zero, respectively.

      cge(1) = xbm_g + 1.
      cge(2) = xmu_g + 1.
      cge(3) = xbm_g + xmu_g + 1.
      do n = 1, 3
         cgg(n) = WGAMMA(cge(n))
      enddo

   mp_select: SELECT CASE(mphysics_opt)

     CASE (WSM6SCHEME)
       scheme_has_graupel = .true.
       xrho_g = 500.
       xam_g = 3.1415926536/6.0*xrho_g
       N0exp = 4.e6
!      !$OMP PARALLEL DO   &
!      !$OMP PRIVATE ( ij )
       DO ij = 1 , num_tiles
         DO j=j_start(ij),j_end(ij)
         DO k=kme-1, kms, -1
         DO i=i_start(ij),i_end(ij)
           if (temp_qg(i,k,j) .LT. 1.E-6) CYCLE
           lam_exp = (N0exp*xam_g*cgg(1)/temp_qg(i,k,j))**(1./cge(1))
           temp_ng(i,k,j) = N0exp*cgg(2)*lam_exp**(-cge(2))
         ENDDO
         ENDDO
         ENDDO
       ENDDO
!      !$OMP END PARALLEL DO

     CASE (WDM6SCHEME)
       scheme_has_graupel = .true.
       xrho_g = 500.
       N0exp = 4.e6
       if (mpuse_hail .eq. 1) then
         xrho_g = 700.
         N0exp = 4.e4
       endif
       xam_g = 3.1415926536/6.0*xrho_g
!      !$OMP PARALLEL DO   &
!      !$OMP PRIVATE ( ij )
       DO ij = 1 , num_tiles
         DO j=j_start(ij),j_end(ij)
         DO k=kme-1, kms, -1
         DO i=i_start(ij),i_end(ij)
           if (temp_qg(i,k,j) .LT. 1.E-6) CYCLE
           lam_exp = (N0exp*xam_g*cgg(1)/temp_qg(i,k,j))**(1./cge(1))
           temp_ng(i,k,j) = N0exp*cgg(2)*lam_exp**(-cge(2))
         ENDDO
         ENDDO
         ENDDO
       ENDDO
!      !$OMP END PARALLEL DO

     CASE (THOMPSON, THOMPSONAERO, TEMPO)

       scheme_has_graupel = .true.
       xmu_g = 1.
       cge(1) = xbm_g + 1.
       cge(2) = xmu_g + 1.
       cge(3) = xbm_g + xmu_g + 1.
       do n = 1, 3
          cgg(n) = WGAMMA(cge(n))
       enddo

!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )
      DO ij = 1 , num_tiles
       DO j=j_start(ij),j_end(ij)
       DO i=i_start(ij),i_end(ij)
        DO k=kme-1, kms, -1
         if (temp_qg(i,k,j) .LT. 1.E-6) CYCLE
         zans1 = (2.5 + 2./7. * (ALOG10(temp_qg(i,k,j))+7.))
         zans1 = MAX(2., MIN(zans1, 7.))
         N0exp = 10.**zans1
         lam_exp = (N0exp*xam_g*cgg(1)/temp_qg(i,k,j))**(1./cge(1))
         lamg = lam_exp * (cgg(3)/cgg(2)/cgg(1))**(1./xbm_g)
         N0_g = N0exp/(cgg(2)*lam_exp) * lamg**cge(2)
         temp_ng(i,k,j) = N0_g*cgg(2)*lamg**(-cge(2))
        ENDDO
       ENDDO
       ENDDO
      ENDDO
!  !$OMP END PARALLEL DO


     CASE (NSSL_2MOM)

       scheme_has_graupel = .true.
       xrho_g = nssl_rho_qh
       N0exp = nssl_cnoh
       if (PRESENT(qh_curr)) then
          xrho_g = nssl_rho_qhl
          N0exp = nssl_cnohl
       endif
       xam_g = 3.1415926536/6.0*xrho_g

       if (PRESENT(ng_curr)) xmu_g = nssl_alphah
       if (PRESENT(nh_curr)) xmu_g = nssl_alphahl

       if (xmu_g .NE. 0.) then
          cge(1) = xbm_g + 1.
          cge(2) = xmu_g + 1.
          cge(3) = xbm_g + xmu_g + 1.
          do n = 1, 3
             cgg(n) = WGAMMA(cge(n))
          enddo
       endif

       ! NSSL scheme has many options, but, if single-moment, just fill
       ! in the number array for graupel from built-in assumptions.

       if (.NOT.(PRESENT(nh_curr).OR.PRESENT(ng_curr)) ) then
!      !$OMP PARALLEL DO   &
!      !$OMP PRIVATE ( ij )
       DO ij = 1 , num_tiles
         DO j=j_start(ij),j_end(ij)
         DO k=kme-1, kms, -1
         DO i=i_start(ij),i_end(ij)
           if (temp_qg(i,k,j) .LT. 1.E-6) CYCLE
           lam_exp = (N0exp*xam_g*cgg(1)/temp_qg(i,k,j))**(1./cge(1))
           temp_ng(i,k,j) = N0exp*cgg(2)*lam_exp**(-cge(2))
         ENDDO
         ENDDO
         ENDDO
       ENDDO
!      !$OMP END PARALLEL DO
       endif

     CASE DEFAULT

   END SELECT mp_select


   if (scheme_has_graupel) then

!..Create bins of graupel/hail (from 500 microns up to 7.5 cm).
      xxDx(1) = 500.D-6
      xxDx(ngbins+1) = 0.075d0
      do n = 2, ngbins
         xxDx(n) = DEXP(DFLOAT(n-1)/DFLOAT(ngbins) &
                  *DLOG(xxDx(ngbins+1)/xxDx(1)) +DLOG(xxDx(1)))
      enddo
      do n = 1, ngbins
         xxDg(n) = DSQRT(xxDx(n)*xxDx(n+1))
         xdtg(n) = xxDx(n+1) - xxDx(n)
      enddo


!  !$OMP PARALLEL DO   &
!  !$OMP PRIVATE ( ij )
      do ij = 1 , num_tiles
        do j=j_start(ij),j_end(ij)
        do k=kms,kme-1
        do i=i_start(ij),i_end(ij)
           if (temp_qg(i,k,j) .lt. 1.e-6) cycle
           lamg = (xam_g*cgg(3)/cgg(2)*temp_ng(i,k,j)/temp_qg(i,k,j))**(1./xbm_g)
           n0_g = temp_ng(i,k,j)/cgg(2)*lamg**cge(2)
           sum_ng = 0.0d0
           sum_t  = 0.0d0
           do ng = ngbins, 1, -1
              f_d = n0_g*xxdg(ng)**xmu_g * dexp(-lamg*xxdg(ng))*xdtg(ng)
              sum_ng = sum_ng + f_d
              if (sum_ng .gt. thresh_conc) then
                 exit
              endif
              sum_t = sum_ng
           enddo
           if (ng .ge. ngbins) then
              hail_max = xxdg(ngbins)
           elseif (xxdg(ng+1) .gt. 1.e-3) then
              hail_max = xxdg(ng) - (sum_ng-thresh_conc)/(sum_ng-sum_t) &
     &                              * (xxdg(ng)-xxdg(ng+1))
           else
              hail_max = 1.e-4
           endif
           if (hail_max .gt. 1e-2) then
            write(outstring,*) 'gt-debug-hail, ', hail_max*1000.
            call wrf_debug (350, trim(outstring))
           endif
           hail_max_sp = hail_max
           if (k.eq.kms) then
              hail_maxk1(i,j) = max(hail_maxk1(i,j), hail_max_sp)
           endif
           hail_max2d(i,j) = max(hail_max2d(i,j), hail_max_sp)
        enddo
        enddo
        enddo
      enddo
!  !$OMP END PARALLEL DO

   endif


! NSSL

   if (diag_print .eq. 0 ) return

   IF ( xtime .ne. 0. ) THEN

      if(diag_print.eq.1) then
         prfreq = dt
         !prfreq = max(2,int(dt/60.))   ! in min
      else
         prfreq=10                   ! in min
      endif
   
      IF (MOD(nint(dt),prfreq) == 0) THEN

        ! COMPUTE THE NUMBER OF MASS GRID POINTS
        no_points = float((ide-ids)*(jde-jds))

        ! SET START AND END POINTS FOR TILES
        !  !$OMP PARALLEL DO   &
        !  !$OMP PRIVATE ( ij )

        dmumax = 0.
        do ij = 1 , num_tiles
          do j=j_start(ij),j_end(ij)
          do i=i_start(ij),i_end(ij)
             dpsdt(i,j)=(p8w(i,kms,j)-pk1m(i,j))/dt
             dmudt(i,j)=(mu_2(i,j)-mu_2m(i,j))/dt
             if(abs(dmudt(i,j)*dt).gt.dmumax)then
               dmumax=abs(dmudt(i,j)*dt)
               idp=i
               jdp=j
             endif
          enddo      
          enddo
        enddo
        !  !$OMP END PARALLEL DO

        ! convert DMUMAX from (PA) to (bars) per time step
        dmumax = dmumax*1.e-5
        ! compute global MAX
        CALL wrf_dm_maxval ( dmumax,  idp, jdp )

        dpsdt_sum = 0.
        dmudt_sum = 0.

        do j = jps, min(jpe,jde-1)
        do i = ips, min(ipe,ide-1)
           dpsdt_sum = dpsdt_sum + abs(dpsdt(i,j))
           dmudt_sum = dmudt_sum + abs(dmudt(i,j))
        enddo
        enddo

        ! compute global sum
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
              ! MAX for accumulated conv precip
              IF(rainc(i,j).gt.raincmax)then
                raincmax=rainc(i,j)
                irc=i
                jrc=j
              ENDIF
              rainnc_sum = rainnc_sum + abs(rainnc(i,j))
              ! MAX for accumulated resolved precip
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

           ! compute global MAX
           CALL wrf_dm_maxval ( raincmax, irc, jrc )
           CALL wrf_dm_maxval ( rainncmax, irnc, jrnc )

           ! compute global sum
           drcdt_sum = wrf_dm_sum_real ( drcdt_sum )
           drndt_sum = wrf_dm_sum_real ( drndt_sum )
           dardt_sum = wrf_dm_sum_real ( dardt_sum )
           rainc_sum = wrf_dm_sum_real ( rainc_sum )
           rainnc_sum = wrf_dm_sum_real ( rainnc_sum )
           raint_sum = wrf_dm_sum_real ( raint_sum )
           sfcevp_sum = wrf_dm_sum_real ( sfcevp_sum )
           hfx_sum = wrf_dm_sum_real ( hfx_sum )
           lh_sum = wrf_dm_sum_real ( lh_sum )

         ENDIF !diag_print .eq. 2

         ! print out the average values
          CALL get_current_grid_name( grid_str )

#ifdef DM_PARALLEL
          IF ( wrf_dm_on_monitor() ) THEN
#endif
          WRITE(outstring,*) grid_str,&
          'Domain average of dpsdt, dmudt (mb/3h): ', xtime, &
          dpsdt_sum/no_points*108., &
          dmudt_sum/no_points*108.
          CALL wrf_message ( TRIM(outstring) )

          WRITE(outstring,*) grid_str,'Max mu change time step: ', idp,jdp,dmumax
          CALL wrf_message ( TRIM(outstring) )

         IF ( diag_print .eq. 2) THEN
            WRITE(outstring,*) grid_str,&
            'Domain average of dardt, drcdt, drndt (mm/sec): ', xtime, &
            dardt_sum/dt/no_points, &
            drcdt_sum/dt/no_points, &
            drndt_sum/dt/no_points
            CALL wrf_message ( TRIM(outstring) )
            WRITE(outstring,*) grid_str,&
            'Domain average of rt_sum, rc_sum, rnc_sum (mm): ', xtime, &
            raint_sum/no_points, &
            rainc_sum/no_points, &
            rainnc_sum/no_points
            CALL wrf_message ( TRIM(outstring) )
            WRITE(outstring,*) grid_str,&
            'Max Accum Resolved Precip,   I,J  (mm): '               ,&
            rainncmax,irnc,jrnc
            CALL wrf_message ( TRIM(outstring) )
            WRITE(outstring,*) grid_str,&
            'Max Accum Convective Precip,   I,J  (mm): '             ,&
            raincmax,irc,jrc
            CALL wrf_message ( TRIM(outstring) )
            WRITE(outstring,*) grid_str,&
            'Domain average of sfcevp, hfx, lh: ', xtime, &
            sfcevp_sum/no_points, &
            hfx_sum/no_points, &
            lh_sum/no_points
            CALL wrf_message ( TRIM(outstring) )
         ENDIF
#ifdef DM_PARALLEL
         ENDIF
#endif

    ENDIF  ! print frequency
   ENDIF   ! xtime .ne. 0.

! save values at this time step
   !$OMP PARALLEL DO   &
   !$OMP PRIVATE ( ij,i,j )
   do ij = 1 , num_tiles

      do j=j_start(ij),j_end(ij)
      do i=i_start(ij),i_end(ij)
         pk1m(i,j)=p8w(i,kms,j)
         mu_2m(i,j)=mu_2(i,j)
      enddo
      enddo

      if ( xtime .lt. 0.0001 ) then
      do j=j_start(ij),j_end(ij)
      do i=i_start(ij),i_end(ij)
         dpsdt(i,j)=0.
         dmudt(i,j)=0.
      enddo
      enddo
      endif

   enddo
   !$omp end parallel do

   end subroutine hrrr_diagnostic_driver


!+---+-----------------------------------------------------------------+ 
      real function gammln(xx)
!     --- returns the value ln(gamma(xx)) for xx > 0.
      implicit none
      real, intent(in):: xx
      double precision, parameter:: stp = 2.5066282746310005d0
      double precision, dimension(6), parameter:: &
               cof = (/76.18009172947146d0, -86.50532032941677d0, &
                       24.01409824083091d0, -1.231739572450155d0, &
                      .1208650973866179d-2, -.5395239384953d-5/)
      double precision:: ser,tmp,x,y
      integer:: j

      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      end function gammln
!  (C) Copr. 1986-92 Numerical Recipes Software 2.02
!+---+-----------------------------------------------------------------+ 
      real function wgamma(y)

      implicit none
      real, intent(in):: y

      wgamma = exp(gammln(y))

      end function wgamma
!+---+-----------------------------------------------------------------+ 

END MODULE module_diag_hrrr

