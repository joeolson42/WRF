!> \file module_bl_mynn_wrapper.F90
!!  This serves as the interface between the WRF PBL driver and the MYNN 
!!  eddy-diffusivity mass-flux scheme in module_bl_mynn.F. 

!>\ingroup gsd_mynn_edmf
!> The following references best describe the code within
!!    Olson et al. (2019, NOAA Technical Memorandum)
!!    Nakanishi and Niino (2009) \cite NAKANISHI_2009
      MODULE module_bl_mynn_wrapper

      use module_bl_mynn_common

      contains

!> \section arg_table_mynnedmf_wrapper_init Argument Table
!! \htmlinclude mynnedmf_wrapper_init.html
!!
      subroutine mynnedmf_wrapper_init (              &
        &  RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,RQCBLTEN,&
        &  RQIBLTEN,QKE,                              &
        &  restart,allowed_to_read,                   &
        &  P_QC,P_QI,PARAM_FIRST_SCALAR,              &
        &  IDS,IDE,JDS,JDE,KDS,KDE,                   &
        &  IMS,IME,JMS,JME,KMS,KME,                   &
        &  ITS,ITE,JTS,JTE,KTS,KTE                    )

        implicit none
        
        LOGICAL,INTENT(IN) :: ALLOWED_TO_READ,RESTART

        INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE,  &
             &                IMS,IME,JMS,JME,KMS,KME,  &
             &                ITS,ITE,JTS,JTE,KTS,KTE


        REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(INOUT) :: &
             &RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,                 &
             &RQCBLTEN,RQIBLTEN,QKE

        INTEGER,  intent(in) :: P_QC,P_QI,PARAM_FIRST_SCALAR

        INTEGER :: I,J,K,ITF,JTF,KTF

        JTF=MIN0(JTE,JDE-1)
        KTF=MIN0(KTE,KDE-1)
        ITF=MIN0(ITE,IDE-1)

        IF (.NOT.RESTART) THEN
         DO J=JTS,JTF
          DO K=KTS,KTF
           DO I=ITS,ITF
              RUBLTEN(i,k,j)=0.
              RVBLTEN(i,k,j)=0.
              RTHBLTEN(i,k,j)=0.
              RQVBLTEN(i,k,j)=0.
              if( p_qc >= param_first_scalar ) RQCBLTEN(i,k,j)=0.
              if( p_qi >= param_first_scalar ) RQIBLTEN(i,k,j)=0.
           ENDDO
          ENDDO
         ENDDO
        ENDIF

      end subroutine mynnedmf_wrapper_init

      subroutine mynnedmf_wrapper_finalize ()
      end subroutine mynnedmf_wrapper_finalize

! \brief This scheme (1) performs pre-mynnedmf work, (2) runs the mynnedmf, and (3) performs post-mynnedmf work
!> \section arg_table_mynnedmf_wrapper_run Argument Table
!! \htmlinclude mynnedmf_wrapper_run.html
!!
SUBROUTINE mynnedmf_wrapper_run(                 &
     &  initflag,restart,cycling,                &
     &  delt,dz,dxc,znt,                         &
     &  u,v,w,th,                                &
     &  qv,qc,qi,qs,qnc,qni,qnwfa,qnifa,qnbca,   &
!     &  ozone,                                   &
     &  p,exner,rho,t3d,                         &
     &  xland,ts,qsfc,ps,                        &
     &  ust,ch,hfx,qfx,rmol,wspd,                &
     &  uoce,voce,                               &
     &  qke,qke_adv,sh3d,sm3d,                   &
!--- chem/smoke
#if (WRF_CHEM == 1)
     &  mix_chem,chem3d,vd3d,nchem,kdvel,        &
     &  ndvel,num_vert_mix,                      &
!     &  frp_mean,emis_ant_no,enh_mix,            & !to be included soon
#endif
!--- end chem/smoke 
     &  tsq,qsq,cov,                             &
     &  rublten,rvblten,rthblten,                &
     &  rqvblten,rqcblten,rqiblten,rqsblten,     &
     &  rqncblten,rqniblten,                     &
     &  rqnwfablten,rqnifablten,rqnbcablten,     &
!     &  ro3blten,                                &
     &  exch_h,exch_m,pblh,kpbl,el_pbl,          &
     &  dqke,qwt,qshear,qbuoy,qdiss,             &
     &  qc_bl,qi_bl,cldfra_bl,                   &
     &  edmf_a,edmf_w,edmf_qt,                   &
     &  edmf_thl,edmf_ent,edmf_qc,               &
     &  sub_thl3d,sub_sqv3d,                     &
     &  det_thl3d,det_sqv3d,                     &
     &  maxwidth,maxMF,ztop_plume,ktop_plume,    &
     &  rthraten,                                &
     &  tke_budget,         bl_mynn_tkeadvect,   &
     &  bl_mynn_cloudpdf,   bl_mynn_mixlength,   &
     &  icloud_bl,          bl_mynn_edmf,        &
     &  bl_mynn_edmf_mom,   bl_mynn_edmf_tke,    &
     &  bl_mynn_cloudmix,   bl_mynn_mixqt,       &
     &  bl_mynn_output,     bl_mynn_closure,     &
     &  bl_mynn_mixscalars,                      &
     &  spp_pbl,pattern_spp_pbl,                 &
     &  flag_qc,flag_qi,flag_qs,                 &
     &  flag_qnc,flag_qni,                       &
     &  flag_qnwfa,flag_qnifa,flag_qnbca,        &
     &  ids,ide,jds,jde,kds,kde,                 &
     &  ims,ime,jms,jme,kms,kme,                 &
     &  its,ite,jts,jte,kts,kte                  )

     use module_bl_mynn, only: mynn_bl_driver

!------------------------------------------------------------------- 
     implicit none
!------------------------------------------------------------------- 

     !smoke/chem: disclaimer: all smoke-related variables are still
     !considered under development in CCPP. Until that work is
     !completed, these flags/arrays must be kept hard-coded as is.
#if (WRF_CHEM == 1)
     logical, intent(in) :: mix_chem
     integer, intent(in) :: nchem, ndvel, kdvel, num_vert_mix
     logical, parameter ::                                  &
     &       rrfs_sd    =.false.,                           &
     &       smoke_dbg  =.false.,                           &
     &       enh_mix    =.false.
#else
     logical, parameter ::                                  &
     &       mix_chem   =.false.,                           &
     &       enh_mix    =.false.,                           &
     &       rrfs_sd    =.false.,                           &
     &       smoke_dbg  =.false.
     integer, parameter :: nchem=2, ndvel=2, kdvel=1,       &
             num_vert_mix = 1
#endif

! NAMELIST OPTIONS (INPUT):
     logical, intent(in) ::                                 &
     &       bl_mynn_tkeadvect,                             &
     &       cycling
      integer, intent(in) ::                                &
     &       bl_mynn_cloudpdf,                              &
     &       bl_mynn_mixlength,                             &
     &       icloud_bl,                                     &
     &       bl_mynn_edmf,                                  &
     &       bl_mynn_edmf_mom,                              &
     &       bl_mynn_edmf_tke,                              &
     &       bl_mynn_cloudmix,                              &
     &       bl_mynn_mixqt,                                 &
     &       bl_mynn_output,                                &
     &       bl_mynn_mixscalars,                            &
     &       spp_pbl,                                       &
     &       tke_budget
      real(kind_phys), intent(in) ::                        &
     &       bl_mynn_closure

      logical, intent(in) ::                                &
     &       FLAG_QI, FLAG_QNI, FLAG_QC, FLAG_QNC,          &
     &       FLAG_QS, FLAG_QNWFA, FLAG_QNIFA, FLAG_QNBCA
      logical, parameter :: FLAG_OZONE = .false.

!MYNN-1D
      REAL(kind_phys),    intent(in) :: delt, dxc
      LOGICAL, intent(in) :: restart
      INTEGER :: i, j, k, itf, jtf, n
      INTEGER, intent(in) :: initflag,                      &
     &           IDS,IDE,JDS,JDE,KDS,KDE,                   &
     &           IMS,IME,JMS,JME,KMS,KME,                   &
     &           ITS,ITE,JTS,JTE,KTS,KTE

!MYNN-3D
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(in) ::               &
     &       u,v,w,t3d,th,rho,exner,p,dz
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(inout) ::            &
     &       rublten,rvblten,rthblten,                                                 &
     &       rqvblten,rqcblten,rqiblten,rqsblten,                                      &
     &       rqncblten,rqniblten,                                                      &
     &       rqnwfablten,rqnifablten,rqnbcablten !,ro3blten
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(inout) ::            &
     &       qke, qke_adv, el_pbl, sh3d, sm3d, tsq, qsq, cov
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(inout) ::            &
     &       exch_h, exch_m
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), intent(in) :: rthraten

!optional 3D arrays
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), optional, intent(in)  ::    &
     &       pattern_spp_pbl
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), optional, intent(inout) ::  &
     &       qc_bl, qi_bl, cldfra_bl
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), optional, intent(inout) ::  &
     &       edmf_a,edmf_w,edmf_qt,                                                    &
     &       edmf_thl,edmf_ent,edmf_qc,                                                &
     &       sub_thl3d,sub_sqv3d,det_thl3d,det_sqv3d
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), optional, intent(inout) ::  &
     &       dqke,qWT,qSHEAR,qBUOY,qDISS
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme), optional, intent(inout) ::  &
     &       qv,qc,qi,qs,qnc,qni,qnwfa,qnifa,qnbca!,o3

!(non-optional) 1D arrays
      real(kind_phys), dimension(kts:kte)  ::                                          &
     &       edmf_a1,edmf_w1,edmf_qt1,edmf_thl1,edmf_ent1,edmf_qc1,                    &
     &       sub_thl1,sub_sqv1,det_thl1,det_sqv1
      real(kind_phys), dimension(kts:kte) ::                                           &
     &       qc_bl1, qi_bl1, cldfra_bl1, pattern_spp_pbl1
      real(kind_phys), dimension(kts:kte) ::                                           &
     &       dqke1,qWT1,qSHEAR1,qBUOY1,qDISS1
      real(kind_phys), dimension(kts:kte) ::                                           &
     &       qc1,qi1,qs1,qnc1,qni1,qnwfa1,qnifa1,qnbca1,ozone1

!smoke/chem arrays - no if-defs in module_bl_mynn.F, so must define arrays
#if (WRF_CHEM == 1)
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme,nchem), intent(in) :: chem3d
      real(kind_phys), dimension(ims:ime,kdvel,jms:jme, ndvel), intent(in)  :: vd3d
      real(kind_phys), dimension(ims:ime,kms:kme,nchem) :: chem
      real(kind_phys), dimension(ims:ime,ndvel)         :: vd
      real(kind_phys), dimension(ims:ime)               :: frp_mean, emis_ant_no
#else
      real(kind_phys), dimension(ims:ime,kms:kme,nchem) :: chem
      real(kind_phys), dimension(ims:ime,ndvel)         :: vd
      real(kind_phys), dimension(ims:ime) :: frp_mean, emis_ant_no
#endif

!MYNN-2D
      real(kind_phys), dimension(ims:ime,jms:jme), intent(in) ::                       &
     &        xland,ts,qsfc,ps,ch
      real(kind_phys), dimension(ims:ime,jms:jme), intent(inout) ::                    &
     &        znt,pblh,maxwidth,maxmf,ztop_plume,rmol,hfx,qfx,ust,wspd,                &
     &        uoce,voce
      integer, dimension(ims:ime,jms:jme), intent(inout) ::                            &
     &        kpbl,ktop_plume

!Local
      real(kind_phys), dimension(kts:kte)                 :: delp,sqv1,sqc1,sqi1,sqs1,ikzero
      real(kind_phys), dimension(ims:ime)                 :: dx
      logical, parameter                                  :: debug = .false.
      real(kind_phys), dimension(ims:ime,kms:kme,jms:jme) :: ozone,rO3blten

   if (debug) then
      write(0,*)"=============================================="
      write(0,*)"in mynn wrapper..."
      write(0,*)"initflag=",initflag," restart =",restart
   endif

   jtf=MIN0(JTE,JDE-1)
   itf=MIN0(ITE,IDE-1)

   !For now, initialized bogus array
   ozone=0.0
   rO3blten=0.0
   ikzero=0.0

   !initialize 1d arrays:
   !subgrid clouds:
   qc_bl1           =0.0
   qi_bl1           =0.0
   cldfra_bl1       =0.0
   !spp
   pattern_spp_pbl1 =0.0
   !tke budget
   dqke1            =0.0
   qWT1             =0.0
   qSHEAR1          =0.0
   qBUOY1           =0.0
   qDISS1           =0.0
   !mass-flux arrays
   edmf_a1          =0.0
   edmf_w1          =0.0
   edmf_qt1         =0.0
   edmf_thl1        =0.0
   edmf_ent1        =0.0
   edmf_qc1         =0.0
   sub_thl1         =0.0
   sub_sqv1         =0.0
   det_thl1         =0.0
   det_sqv1         =0.0
   !moist species
   qc1              =0.0
   qi1              =0.0
   qs1              =0.0
   qnc1             =0.0
   qni1             =0.0
   qnwfa1           =0.0
   qnifa1           =0.0
   qnbca1           =0.0

   !---------------------------------
   !Begin looping in the i- and j-direction
   !---------------------------------
   do j = jts,jtf
   do i = its,itf

      !need sgs cloud info input for diagnostic-decay
      if (icloud_bl > 0) then
         do k=kts,kte
            qc_bl1(k)     = qc_bl(i,k,j)
            qi_bl1(k)     = qi_bl(i,k,j)
            cldfra_bl1(k) = cldfra_bl(i,k,j)
         enddo
      endif

      !spp input
      if (spp_pbl > 0) then
         do k=kts,kte
            pattern_spp_pbl1(k) = pattern_spp_pbl(i,k,j)
         enddo
      endif

      !intialize moist species
      if (flag_qc) then
         do k=kts,kte
            qc1(k) = qc(i,k,j)
         enddo
      endif
      if (flag_qi) then
         do k=kts,kte
            qi1(k) = qi(i,k,j)
         enddo
      endif
      if (flag_qs) then
         do k=kts,kte
            qs1(k) = qs(i,k,j)
         enddo
      endif
      if (flag_qnc) then
         do k=kts,kte
            qnc1(k) = qnc(i,k,j)
         enddo
      endif
      if (flag_qni) then
         do k=kts,kte
            qni1(k) = qni(i,k,j)
         enddo
      endif
      if (flag_qnwfa) then
         do k=kts,kte
            qnwfa1(k) = qnwfa(i,k,j)
         enddo
      endif
      if (flag_qnifa) then
         do k=kts,kte
            qnifa1(k) = qnifa(i,k,j)
         enddo
      endif
      if (flag_qnbca) then
         do k=kts,kte
            qnbca1(k) = qnbca(i,k,j)
         enddo
      endif

#if (WRF_CHEM == 1)
      if (mix_chem) then
         do n=1,nchem
         do k=kts,kte
            chem(k,n)=chem3d(i,k,j,n)
         enddo
         enddo

         !set kdvel =1
         do n=1,ndvel
            vd(n)  =vd3d(i,1,j,n)
         enddo
      endif
      frp_mean = 0.0
      emis_ant_no = 0.0
#else
      chem     = 0.0
      vd       = 0.0
      frp_mean = 0.0
      emis_ant_no = 0.0
#endif


      !find/fix negative mixing ratios
!      call moisture_check2(kte, delt,                 &
!                           delp(i,:), exner(i,:,j),   &
!                           qv(i,:,j), qc(i,:,j),      &
!                           qi(i,:,j), t3d(i,:,j)      )

      !In WRF, mixing ratio is incoming. Convert to specific humidity:
      do k=kts,kte
         sqv1(k)=qv(i,k,j)/(1.0 + qv(i,k,j))
         sqc1(k)=qc1(k)/(1.0 + qv(i,k,j))
         sqi1(k)=qi1(k)/(1.0 + qv(i,k,j))
      enddo
      if (flag_qs) then
         do k=kts,kte
            sqs1(k)=max(qs1(k)/(1.0 + qv(i,k,j)), 0.0)
         enddo
      else
         sqs1=0.0
      endif

      if (debug) then
         print*
         write(0,*)"===CALLING mynn_bl_driver; input:"
         print*,"tke_budget=",tke_budget
         print*,"bl_mynn_tkeadvect=",bl_mynn_tkeadvect
         print*,"bl_mynn_cloudpdf=",bl_mynn_cloudpdf
         print*,"bl_mynn_mixlength=",bl_mynn_mixlength
         print*,"bl_mynn_edmf=",bl_mynn_edmf
         print*,"bl_mynn_edmf_mom=",bl_mynn_edmf_mom
         print*,"bl_mynn_edmf_tke=",bl_mynn_edmf_tke
         print*,"bl_mynn_cloudmix=",bl_mynn_cloudmix
         print*,"bl_mynn_mixqt=",bl_mynn_mixqt
         print*,"icloud_bl=",icloud_bl
         print*,"T:",t3d(i,1,j),t3d(i,2,j),t3d(i,kte,j)
         print*,"TH:",th(i,1,j),th(i,2,j),th(i,kte,j)
         print*,"rho:",rho(i,1,j),rho(i,2,j),rho(i,kte,j)
         print*,"exner:",exner(i,1,j),exner(i,2,j),exner(i,kte,j)
         print*,"p:",p(i,1,j),p(i,2,j),p(i,kte,j)
         print*,"dz:",dz(i,1,j),dz(i,2,j),dz(i,kte,j)
         print*,"u:",u(i,1,j),u(i,2,j),u(i,kte,j)
         print*,"v:",v(i,1,j),v(i,2,j),v(i,kte,j)
         print*,"sqv:",sqv1(1),sqv1(2),sqv1(kte)
         print*,"sqc:",sqc1(1),sqc1(2),sqc1(kte)
         print*,"sqi:",sqi1(1),sqi1(2),sqi1(kte)
         print*,"rmol:",rmol(i,j)," ust:",ust(i,j)
         print*,"dx=",dx(i),"initflag=",initflag
         print*,"Thetasurf:",ts(i,j)
         print*,"HFX:",hfx(i,j)," qfx",qfx(i,j)
         print*,"qsfc:",qsfc(i,j)," ps:",ps(i,j)
         print*,"wspd:",wspd(i,j)
         print*,"znt:",znt(i,j)," delt=",delt
         print*,"ite=",ite," kte=",kte
         print*,"PBLH=",pblh(i,j)," KPBL=",KPBL(i,j)," xland=",xland(i,j)
         print*,"ch=",ch(i,j)
         print*,"qke:",qke(i,1,j),qke(i,2,j),qke(i,kte,j)
         print*,"el_pbl:",el_pbl(i,1,j),el_pbl(i,2,j),el_pbl(i,kte,j)
         print*,"Sh3d:",Sh3d(i,1,j),sh3d(i,2,j),sh3d(i,kte,j)
         print*,"max cf_bl:",maxval(cldfra_bl(i,:,j))
      endif

!print*,"In mynn wrapper, calling mynn_bl_driver"
              CALL  mynn_bl_driver(                                    &
     &        initflag=initflag,         restart=restart,              &
     &        cycling=cycling,           delt=delt,                    &
     &        dz1=dz(i,kts:kte,j),       dx=dxc,                       &
     &        znt=znt(i,j),               &
     &        u1=u(i,kts:kte,j),         v1=v(i,kts:kte,j),            &
     &        w1=w(i,kts:kte,j),         th1=th(i,kts:kte,j),          &
     &        sqv1=sqv1,                 sqc1=sqc1,                    &
     &        sqi1=sqi1,                 sqs1=sqs1,                    &
     &        qnc1=qnc1,                 qni1=qni1,                    &
     &        qnwfa1=qnwfa1,             qnifa1=qnifa1,                &
     &        qnbca1=qnbca1,             ozone1=ozone(i,kts:kte,j),    &
     &        p1=p(i,kts:kte,j),         ex1=exner(i,kts:kte,j),       &
     &        rho1=rho(i,kts:kte,j),     tk1=t3d(i,kts:kte,j),         &
              !2d variables
     &        xland=xland(i,j),          ts=ts(i,j),                   &
     &        qsfc=qsfc(i,j),            ps=ps(i,j),                   &
     &        ust=ust(i,j),              ch=ch(i,j),                   &
     &        hfx=hfx(i,j),              qfx=qfx(i,j),                 &
     &        rmol=rmol(i,j),            wspd=wspd(i,j),               &
     &        pblh=pblh(i,j),            KPBL=KPBL(i,j),               & !output
     &        uoce=uoce(i,j),            voce=voce(i,j),               & !input
              !smoke/chem
     &        nchem=nchem,kdvel=kdvel,ndvel=ndvel,                     &
     &        Chem3d=chem,Vdep=vd,                                     &
     &        FRP=frp_mean(i),EMIS_ANT_NO=emis_ant_no(i),              &
     &        mix_chem=mix_chem,enh_mix=enh_mix,                       &
     &        rrfs_sd=rrfs_sd,smoke_dbg=smoke_dbg,                     &
              !higher-order variables and subgrid clouds
     &        tsq1=tsq(i,kts:kte,j),     qsq1=qsq(i,kts:kte,j),        & !output
     &        cov1=cov(i,kts:kte,j),     cldfra_bl1=cldfra_bl1,        & !output
     &        qc_bl1=qc_bl1,             qi_bl1=qi_bl1,                & !output
              !eddy diffusivity/viscosity and components
     &        kh1=exch_h(i,kts:kte+1,j), km1=exch_m(i,kts:kte+1,j),    & !output
     &        qke1=QKE(i,kts:kte,j),     qke_adv1=qke_adv(i,kts:kte,j),& !output
     &        sh1=Sh3d(i,kts:kte,j),     sm1=Sm3d(i,kts:kte,j),        & !output
     &        el1=el_pbl(i,kts:kte,j),                                 &
              !tendencies
     &        du1=RUBLTEN(i,kts:kte,j),  dv1=RVBLTEN(i,kts:kte,j),     & !output
     &        dth1=RTHBLTEN(i,kts:kte,j),dqv1=RQVBLTEN(i,kts:kte,j),   & !output
     &        dqc1=rqcblten(i,kts:kte,j),dqi1=rqiblten(i,kts:kte,j),   & !output
     &        dqnc1=rqncblten(i,kts:kte,j),dqni1=rqniblten(i,kts:kte,j),& !output
     &        dqs1=ikzero,                                             & !there is no RQSBLTEN, so use dummy arary
     &        dqnwfa1=RQNWFABLTEN(i,kts:kte,j),dqnifa1=RQNIFABLTEN(i,kts:kte,j),   & !output
     &        dqnbca1=RQNBCABLTEN(i,kts:kte,j),dozone1=rO3blten(i,kts:kte,j),      & !output
              !tke budget
     &        dqke1=dqke1,                                             & !output
     &        qWT1=qWT1,                 qSHEAR1=qSHEAR1,              & !output
     &        qBUOY1=qBUOY1,             qDISS1=qDISS1,                & !output
              !input namelist parameters
     &        bl_mynn_tkeadvect  = bl_mynn_tkeadvect,                  & !input parameter
     &        tke_budget         = tke_budget,                         & !input parameter
     &        bl_mynn_cloudpdf   = bl_mynn_cloudpdf,                   & !input parameter
     &        bl_mynn_mixlength  = bl_mynn_mixlength,                  & !input parameter
     &        icloud_bl          = icloud_bl,                          & !input parameter
     &        closure            = bl_mynn_closure,                    & !input parameter
     &        bl_mynn_edmf       = bl_mynn_edmf,                       & !input parameter
     &        bl_mynn_edmf_mom   = bl_mynn_edmf_mom,                   & !input parameter
     &        bl_mynn_edmf_tke   = bl_mynn_edmf_tke,                   & !input parameter
     &        bl_mynn_mixscalars = bl_mynn_mixscalars,                 & !input parameter
     &        bl_mynn_output     = bl_mynn_output,                     & !input parameter
     &        bl_mynn_cloudmix   = bl_mynn_cloudmix,                   & !input parameter
     &        bl_mynn_mixqt      = bl_mynn_mixqt,                      & !input parameter
     &        edmf_a1=edmf_a1,           edmf_w1=edmf_w1,              & !output
     &        edmf_qt1=edmf_qt1,         edmf_thl1=edmf_thl1,          & !output
     &        edmf_ent1=edmf_ent1,       edmf_qc1=edmf_qc1,            & !output
     &        sub_thl1=sub_thl1,         sub_sqv1=sub_sqv1,            & !output
     &        det_thl1=det_thl1,         det_sqv1=det_sqv1,            & !output
     &        maxwidth=maxwidth(i,j),    maxmf=maxMF(i,j),             & !output
     &        ztop_plume=ztop_plume(i,j),ktop_plume=ktop_plume(i,j),   & !output
     &        spp_pbl=spp_pbl, pattern_spp_pbl1=pattern_spp_pbl1,      & !input
     &        RTHRATEN1=rthraten(i,kts:kte,j),                         & !input
     &        FLAG_QI=flag_qi,FLAG_QNI=flag_qni,FLAG_QS=flag_qs,       & !input
     &        FLAG_QC=flag_qc,FLAG_QNC=flag_qnc,                       & !input
     &        FLAG_QNWFA=FLAG_QNWFA,FLAG_QNIFA=FLAG_QNIFA,             & !input
     &        FLAG_QNBCA=FLAG_QNBCA,FLAG_OZONE=flag_ozone,             & !input
     &        IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde,         & !input
     &        IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme,         & !input
     &        ITS=its,ITE=itf,JTS=jts,JTE=jtf,KTS=kts,KTE=kte          ) !input
!print*,"In mynn wrapper, after bl_mynn_driver"

!SHOULD CONVERT BACK WITH THE UPDATED QV??????????
     !- Convert spec hum to mixing ratio:
      do k=kts,kte
         RQVBLTEN(i,k,j) = RQVBLTEN(i,k,j)/(1.0 - sqv1(k))
         RQCBLTEN(i,k,j) = RQCBLTEN(i,k,j)/(1.0 - sqv1(k))
         RQIBLTEN(i,k,j) = RQIBLTEN(i,k,j)/(1.0 - sqv1(k))
      enddo
      if (.false.) then !as of now, there is no RQSBLTEN in WRF
        do k=kts,kte
           RQSBLTEN(i,k,j) = RQSBLTEN(i,k,j)/(1.0 - sqv1(k))
        enddo
      endif

     !- Collect 3D ouput:
      if (icloud_bl > 0) then
         do k=kts,kte
            qc_bl(i,k,j)     = qc_bl1(k)/(1.0 - sqv1(k))
            qi_bl(i,k,j)     = qi_bl1(k)/(1.0 - sqv1(k))
            cldfra_bl(i,k,j) = cldfra_bl1(k)
         enddo
      endif

      if (tke_budget .eq. 1) then
         do k=kts,kte
            dqke(i,k,j)     = dqke1(k)
            qwt(i,k,j)      = qwt1(k)
            qshear(i,k,j)   = qshear1(k)
            qbuoy(i,k,j)    = qbuoy1(k)
            qdiss(i,k,j)    = qdiss1(k)
         enddo
      endif

      if (bl_mynn_output > 0) then
         do k=kts,kte
            edmf_a(i,k,j)     = edmf_a1(k)
            edmf_w(i,k,j)     = edmf_w1(k)
            edmf_qt(i,k,j)    = edmf_qt1(k)
            edmf_thl(i,k,j)   = edmf_thl1(k)
            edmf_ent(i,k,j)   = edmf_ent1(k)
            edmf_qc(i,k,j)    = edmf_qc1(k)
            sub_thl3d(i,k,j)  = sub_thl1(k)
            sub_sqv3d(i,k,j)  = sub_sqv1(k)
            det_thl3d(i,k,j)  = det_thl1(k)
            det_sqv3d(i,k,j)  = det_sqv1(k)
         enddo
      endif

      if (debug) then
          print*
          print*,"===Finished with mynn_bl_driver; output:"
          print*,"T:",t3d(i,1,j),t3d(i,2,j),t3d(i,kte,j)
          print*,"TH:",th(i,1,j),th(i,2,j),th(i,kte,j)
          print*,"rho:",rho(i,1,j),rho(i,2,j),rho(i,kte,j)
          print*,"exner:",exner(i,1,j),exner(i,2,j),exner(i,kte,j)
          print*,"p:",p(i,1,j),p(i,2,j),p(i,kte,j)
          print*,"dz:",dz(i,1,j),dz(i,2,j),dz(i,kte,j)
          print*,"u:",u(i,1,j),u(i,2,j),u(i,kte,j)
          print*,"v:",v(i,1,j),v(i,2,j),v(i,kte,j)
          print*,"sqv:",sqv1(1),sqv1(2),sqv1(kte)
          print*,"sqc:",sqc1(1),sqc1(2),sqc1(kte)
          print*,"sqi:",sqi1(1),sqi1(2),sqi1(kte)
          print*,"rmol:",rmol(i,j)," ust:",ust(i,j)
          print*,"dx(i,j)=",dx(i),"initflag=",initflag
          print*,"Thetasurf:",ts(i,j)
          print*,"HFX:",hfx(i,j)," qfx",qfx(i,j)
          print*,"qsfc:",qsfc(i,j)," ps:",ps(i,j)
          print*,"wspd:",wspd(i,j)
          print*,"znt:",znt(i,j)," delt=",delt
          print*,"im=",ite," kte=",kte
          print*,"PBLH=",pblh(i,j)," KPBL=",KPBL(i,j)," xland=",xland(i,j)
          print*,"ch=",ch(i,j)
          print*,"qke:",qke(i,1,j),qke(i,2,j),qke(i,kte,j)
          print*,"el_pbl:",el_pbl(i,1,j),el_pbl(i,2,j),el_pbl(i,kte,j)
          print*,"Sh3d:",Sh3d(i,1,j),sh3d(i,2,j),sh3d(i,kte,j)
          print*,"exch_h:",exch_h(i,1,j),exch_h(i,2,j),exch_h(i,kte,j)
          print*,"exch_m:",exch_m(i,1,j),exch_m(i,2,j),exch_m(i,kte,j)
          print*,"max cf_bl:",maxval(cldfra_bl(i,:,j))
          print*,"max qc_bl:",maxval(qc_bl(i,:,j))
          print*,"dtdt:",rthblten(i,1,j),rthblten(i,2,j),rthblten(i,kte,j)
          print*,"dudt:",rublten(i,1,j),rublten(i,2,j),rublten(i,kte,j)
          print*,"dvdt:",rvblten(i,1,j),rvblten(i,2,j),rvblten(i,kte,j)
          print*,"dqdt:",rqvblten(i,1,j),rqvblten(i,2,j),rqvblten(i,kte,j)
          print*,"ztop_plume:",ztop_plume(i,j)," maxmf:",maxmf(i,j)
          print*
       endif

   enddo  !end j-loop
   enddo  !end i-loop

   if (debug) then
      print*,"In mynn wrapper, at end"
   endif

  CONTAINS

! ==================================================================
  SUBROUTINE moisture_check2(kte, delt, dp, exner, &
                             qv, qc, qi, th        )
  !
  ! If qc < qcmin, qi < qimin, or qv < qvmin happens in any layer,
  ! force them to be larger than minimum value by (1) condensating 
  ! water vapor into liquid or ice, and (2) by transporting water vapor 
  ! from the very lower layer.
  ! 
  ! We then update the final state variables and tendencies associated
  ! with this correction. If any condensation happens, update theta/temperature too.
  ! Note that (qv,qc,qi,th) are the final state variables after
  ! applying corresponding input tendencies and corrective tendencies.

    implicit none
    integer,  intent(in)     :: kte
    real, intent(in)     :: delt
    real, dimension(kte), intent(in)     :: dp
    real, dimension(kte), intent(in)     :: exner
    real, dimension(kte), intent(inout)  :: qv, qc, qi, th
    integer   k
    real ::  dqc2, dqi2, dqv2, sum, aa, dum
    real, parameter :: qvmin1= 1e-8,    & !min at k=1
                       qvmin = 1e-20,   & !min above k=1
                       qcmin = 0.0,     &
                       qimin = 0.0

    do k = kte, 1, -1  ! From the top to the surface
       dqc2 = max(0.0, qcmin-qc(k)) !qc deficit (>=0)
       dqi2 = max(0.0, qimin-qi(k)) !qi deficit (>=0)

       !update species
       qc(k)  = qc(k)  +  dqc2
       qi(k)  = qi(k)  +  dqi2
       qv(k)  = qv(k)  -  dqc2 - dqi2
       !for theta
       !th(k)  = th(k)  +  xlvcp/exner(k)*dqc2 + &
       !                   xlscp/exner(k)*dqi2
       !for temperature
       th(k)  = th(k)  +  xlvcp*dqc2 + &
                          xlscp*dqi2

       !then fix qv if lending qv made it negative
       if (k .eq. 1) then
          dqv2   = max(0.0, qvmin1-qv(k)) !qv deficit (>=0)
          qv(k)  = qv(k)  + dqv2
          qv(k)  = max(qv(k),qvmin1)
          dqv2   = 0.0
       else
          dqv2   = max(0.0, qvmin-qv(k)) !qv deficit (>=0)
          qv(k)  = qv(k)  + dqv2
          qv(k-1)= qv(k-1)  - dqv2*dp(k)/dp(k-1)
          qv(k)  = max(qv(k),qvmin)
       endif
       qc(k) = max(qc(k),qcmin)
       qi(k) = max(qi(k),qimin)
    end do

    ! Extra moisture used to satisfy 'qv(1)>=qvmin' is proportionally
    ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
    ! preserves column moisture.
    if( dqv2 .gt. 1.e-20 ) then
        sum = 0.0
        do k = 1, kte
           if( qv(k) .gt. 2.0*qvmin ) sum = sum + qv(k)*dp(k)
        enddo
        aa = dqv2*dp(1)/max(1.e-20,sum)
        if( aa .lt. 0.5 ) then
            do k = 1, kte
               if( qv(k) .gt. 2.0*qvmin ) then
                   dum    = aa*qv(k)
                   qv(k)  = qv(k) - dum
               endif
            enddo
        else
        ! For testing purposes only (not yet found in any output):
        !    write(*,*) 'Full moisture conservation is impossible'
        endif
    endif

    return

  END SUBROUTINE moisture_check2

  END SUBROUTINE mynnedmf_wrapper_run

!###=================================================================

END MODULE module_bl_mynn_wrapper
