   MODULE global_parameters
   
    IMPLICIT NONE 
 
    INTEGER, PARAMETER :: INT_K  = 1
    INTEGER, PARAMETER :: CGREAL = SELECTED_REAL_KIND(P=14,R=30) 

    ! flow_type   
    INTEGER, PARAMETER :: VISCOUS_FLOW             = 1, &
                          POTENTIAL_FLOW           = 2
    ! [xyz]grid_unif
    INTEGER, PARAMETER :: UNIFORM_GRID             = 1, &
                          NONUNIFORM_GRID          = 2
    ! bcx1, bcx2 ...     
    INTEGER, PARAMETER :: BC_TYPE_DIRICHLET        = 1, & 
                          BC_TYPE_ZERO_GRADIENT    = 2, &
                          BC_TYPE_PULSATILE_INFLOW = 3, &  
                          BC_TYPE_SYMMETRY         = 4, &
                          BC_TYPE_PERIODIC         = 5, &
                          BC_TYPE_USER_SPECIFIED   = 6, &
                          BC_TYPE_SHEAR            = 7
      
    INTEGER, PARAMETER :: IT_SOLVER_TYPE_LSOR      = 1, & 
                          IT_SOLVER_TYPE_PETSC     = 2, &
                          IT_SOLVER_TYPE_MG        = 3, &
                          IT_SOLVER_TYPE_AZ        = 4
              
    INTEGER, PARAMETER :: FIXED_BOUNDARY           = 1, & 
                          MOVING_BOUNDARY          = 2

    ! boundary_motion_type              
    INTEGER, PARAMETER :: STATIONARY               = 0, &
                          PRESCRIBED_RIGID         = 1, & 
                          INDUCED_RIGID            = 2, & 
                          PRESCRIBED_ARBITRARY     = 3, &
                          ELASTIC_MOTION           = 4

    ! forced_motion_spec
    INTEGER, PARAMETER :: SPECIFY_VELOCITY         = 0, &
                          SPECIFY_POSITION         = 1
    ! pbcx1, pbcx2 ...
    INTEGER, PARAMETER :: PBC_DIRICHLET            = 1, &
                          PBC_NEUMANN              = 2

    ! advec_scheme
    INTEGER, PARAMETER :: ADAMS_BASHFORTH2         = 1, &
                          RUNGE_KUTTA3             = 2, &
                          CRANK_NICOLSON2          = 3   
    ! wall_type
    INTEGER, PARAMETER :: NONPOROUS_AND_NONSLIP    = 0, &
                          POROUS_OR_SLIP           = 1
              
    INTEGER, PARAMETER :: NONE                     = 0, &
                          GENERAL                  = 1, &
                          CANONICAL                = 2
    ! canonical_body_type
    INTEGER, PARAMETER :: ELLIPTIC_CYLINDER        = 1, &
                          GENERAL_CYLINDER         = 2, &
                          ELLIPSOID                = 3, &
                          UNSTRUCTURED_SURFACE     = 4

    ! unstruc_surface_type
    INTEGER, PARAMETER :: SOLID_BODY               = 1, &
                          MEMBRANE                 = 2
   
    INTEGER, PARAMETER :: BODY_DIM2                = 2, &
                          BODY_DIM3                = 3
    
    INTEGER, PARAMETER :: INTR_BOUND_NONE          = 0, &
                          INTR_BOUND_PRESENT       = 1
    ! frac_step                                   
    INTEGER, PARAMETER :: NO_VAN_KAN               = 0, &
                          VAN_KAN                  = 1
    ! format_dump                                   
    INTEGER, PARAMETER :: TECPLOT                  = 1, & 
                          FIELDVIEW                = 2, & 
                          MATLAB                   = 3 
              
    INTEGER, PARAMETER :: IBLANK_READ              = 1

    INTEGER, PARAMETER :: DIM_2D                   = 2, &
                          DIM_3D                   = 3

    INTEGER, PARAMETER :: IBLANK_USED              = 1
    
    INTEGER, PARAMETER :: NO_INTERNAL_BOUNDARY     = 0, &
                          SSM_METHOD               = 1, &
                          GCM_METHOD               = 2
    
    INTEGER, PARAMETER :: INVISCID                 = 1                  

    INTEGER, PARAMETER :: ICOORD                   = 1, &
                          JCOORD                   = 2, &
                          KCOORD                   = 3

    INTEGER, PARAMETER :: STATS_NONE               = 0

    INTEGER, PARAMETER :: STDOUT                   = 6

    INTEGER, PARAMETER :: INACTIVE                 = 0, &
                          ACTIVE                   = 1, &
                          ERR_NONE                 = 0            

    INTEGER, PARAMETER :: YES                      = 1, &
                          NO                       = 0

   END MODULE global_parameters
!------------------------------------------------------

   MODULE flow_parameters 

    USE global_parameters
    
    IMPLICIT NONE
 
    REAL(KIND=CGREAL), PARAMETER :: sidw  = 2.0_CGREAL ! parameter used in IDW interpolation
    REAL(KIND=CGREAL)            :: pi

    INTEGER  :: nread
    INTEGER  :: ndim
    INTEGER  :: flow_type
    INTEGER  :: nx,ny,nz
    INTEGER  :: xgrid_unif,ygrid_unif,zgrid_unif
    INTEGER  :: bcx1,bcx2,bcy1,bcy2,bcz1,bcz2
    INTEGER  :: no_tsteps,nmonitor,ndump,nrestart,nstat,nmonitor_probe_liftdrag
    INTEGER  :: format_dump
    INTEGER  :: ntime,ntime_start
    INTEGER  :: it_solver_type, advec_scheme
    
    INTEGER  :: boundary_motion, nbody, body_type 
    INTEGER  :: nPtsMax, nGhostMax, nDeadMax
    INTEGER  :: boundary_formulation=GCM_METHOD
    INTEGER  :: mix_GC_form, elastic_present
    INTEGER  :: itermax_pson, itermax_ad, itermax_gc
    
    INTEGER  :: ifuBodyIn, ifuInput, ifuIblankIn, ifuMarkerIn, ifuProbeIn, &
                ifuRstrtFlowIn, ifuRstrtBodyIn, ifuUnstrucSurfIn, ifuUnstrucSurfOut, &
                ifuBodyOut
    INTEGER  :: ifuDragOut, ifuFreshCellOut, ifuMarkerOut, ifuProbeOut,     &
                ifuRstrtFlowOut , ifuRstrtBodyOut, ifuStatOut, ifuStatPlot
    
    INTEGER  :: internal_boundary_present, frac_step_type
    
    INTEGER, DIMENSION(:), ALLOCATABLE  :: nPtsBodyMarkerOrig,canonical_body_type, &
                                           body_dim,boundary_motion_type,wall_type
    INTEGER, DIMENSION(:), ALLOCATABLE  :: unstruc_surface_type

    INTEGER, DIMENSION(:), ALLOCATABLE  :: nPtsBodyMarker
    
    REAL(KIND=CGREAL)  :: xOrigin,yOrigin,zOrigin
    REAL(KIND=CGREAL)  :: xout,yout,zout 
    REAL(KIND=CGREAL)  :: uinit,vinit,winit
    REAL(KIND=CGREAL)  :: ux1,ux2,vx1,vx2,wx1,wx2
    REAL(KIND=CGREAL)  :: uy1,uy2,vy1,vy2,wy1,wy2
    REAL(KIND=CGREAL)  :: uz1,uz2,vz1,vz2,wz1,wz2
    REAL(KIND=CGREAL)  :: freq_ux1,freq_vx1,freq_wx1
    REAL(KIND=CGREAL)  :: freq_ux2,freq_vx2,freq_wx2
    REAL(KIND=CGREAL)  :: freq_uy1,freq_vy1,freq_wy1
    REAL(KIND=CGREAL)  :: freq_uy2,freq_vy2,freq_wy2
    REAL(KIND=CGREAL)  :: freq_uz1,freq_vz1,freq_wz1
    REAL(KIND=CGREAL)  :: freq_uz2,freq_vz2,freq_wz2
  
    REAL(KIND=CGREAL)  :: re,dt,reinv,dtinv
    
    REAL(KIND=CGREAL)  :: restol_ad, restol_gcp, restol_gcv, omega_pson, omega_ad
    REAL(KIND=CGREAL)  :: restol_pson, restol_phiPrime, restol_div, resDiv
    
    REAL(KIND=CGREAL)  :: time, time_start
    REAL(KIND=CGREAL)  :: bodyInterceptWeight, imagePointWeight, probeLengthNormalized, &
                          probeLengthNormalizedD, gcmFlag

    REAL(KIND=CGREAL)  :: membrane_tkns_factor, membrane_tkns

    REAL(KIND=CGREAL)  :: area_left,area_right,  &
                          area_bot,area_top,     &
                          area_back,area_front,  &
                          outflow_area
    
    ! variable for rigid body kinematics                             
    
    INTEGER, DIMENSION(:), ALLOCATABLE  :: n_theta,n_phi

    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: radiusx,radiusy,radiusz
     
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: vxcent,vycent,vzcent,  &
                                                   angvx,angvy,angvz,     &
                                                   xcent,ycent,zcent

    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: angvx_old,angvy_old,angvz_old
     
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: xcentinit,ycentinit,zcentinit,        &
                                                   vxcentTrans,vycentTrans,vzcentTrans,  &
                                                   ampx,ampy,ampz,freqx,freqy,freqz
                                                        
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: angvxinit,angvyinit,angvzinit,  &
                                                   ampangx,ampangy,ampangz,        &
                                                   freqangx,freqangy,freqangz

    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: angx ,angy , angz,              &
                                                   angxinit,angyinit,angzinit
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: xPivot,yPivot,zPivot

    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: vxphase,vyphase,vzphase
    REAL(KIND=CGREAL),DIMENSION(:), ALLOCATABLE :: angxphase,angyphase,angzphase
    INTEGER,          DIMENSION(:), ALLOCATABLE :: forced_motion_spec

    ! For Flow-Induced Motion.
    REAL(KIND=CGREAL) ,DIMENSION(:),ALLOCATABLE  :: xcentConstr,ycentConstr,zcentConstr ! Centroid Constraint Flag
    
    REAL(KIND=CGREAL) :: density_fluid 
    
    REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE  :: density_solid,zoom_factor

    INTEGER            :: is_PBC_homogeneous
    INTEGER            :: pbcx1,pbcx2, pbcy1,pbcy2, pbcz1,pbcz2  !
    REAL(KIND=CGREAL)  :: pppx1,pppx2, pppy1,pppy2, pppz1,pppz2  !
    LOGICAL            :: smooth_divergence                      !
    REAL(KIND=CGREAL)  :: alfa_upwind, x1_upwind                 !

    REAL(KIND=CGREAL)  :: alfa        ! Weighting factor for hybrid scheme - added by Rupesh
    REAL(KIND=CGREAL)  :: vper        ! for 2D-3D initial random perturbations
  
    ! probe parameters
    INTEGER                           :: Imonitor_probe, nProbeFlow
    INTEGER, DIMENSION(:),ALLOCATABLE :: iProbe, jProbe, kProbe

    INTEGER                           :: Imonitor_drag

    INTEGER, PARAMETER                :: nIntpDataMax = 20, &
                                         iRowMax      = 8
                                         
    ! parameters for Runge-Kutta scheme
    INTEGER                         :: rk_step
    REAL(KIND=CGREAL), DIMENSION(3) :: rk_bet, rk_gam, rk_zet


    REAL(KIND=CGREAL) :: I_XX,I_YY,I_ZZ,I_XY,I_YZ,I_XZ

    REAL(KIND=CGREAL) :: moment_x,moment_y,moment_z
    REAL(KIND=CGREAL) :: force_x,force_y,force_z

    REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE  :: cxt,cyt,czt  ! total force on a body
    REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE  :: cxs,cys,czs  ! total shear on a body
    REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE  :: cmxt,cmyt,cmzt ! total torque on a body
    REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE  :: cmxs,cmys,cmzs ! total shear torque on a body

    REAL(KIND=CGREAL) :: vq1_i, vq1
    REAL(KIND=CGREAL) :: vq2_i, vq2
    REAL(KIND=CGREAL) :: vq3_i, vq3
    REAL(KIND=CGREAL) :: vq4_i, vq4

    !Ye,goa
    REAL(KIND=CGREAL) :: goa_area_i, goa_area

   END MODULE flow_parameters
!------------------------------------------------------

   MODULE grid_arrays

    USE global_parameters
    
    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: x,y,z,xc,yc,zc
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dx,dy,dz
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dxinv,dyinv,dzinv
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dxc,dyc,dzc
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: dxcinv,dycinv,dzcinv
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: fx,fy,fz
   
    REAL(KIND=CGREAL) :: dxMin,dyMin,dzMin

   END MODULE grid_arrays
!------------------------------------------------------

   MODULE flow_arrays

    USE global_parameters
    
    IMPLICIT NONE
   
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: u,v,w
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: uast,vast,wast

    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: face_u,face_v,face_w
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bcxu,bcxv,bcxw
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bcyu,bcyv,bcyw
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bczu,bczv,bczw
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bcxp,bcyp,bczp
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: viscTot
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: bcxvisc,bcyvisc,bczvisc

    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: p, pPrime, p_temp
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: dpdx,dpdy,dpdz

    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: uOld, vOld, wOld    
   END MODULE flow_arrays
!------------------------------------------------------

   MODULE multiuse_arrays

    USE global_parameters
    
    IMPLICIT NONE
   
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: nlu,nlv,nlw
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: nluold,nlvold,nlwold
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE :: div, lmd

   END MODULE multiuse_arrays
!------------------------------------------------------

   MODULE boundary_arrays

    USE global_parameters
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER       :: NSIZE = 1000
    INTEGER, PARAMETER       :: MSIZE = 20
    
    INTEGER  :: nFresh, nDead, iblank_reset 
    INTEGER, DIMENSION(:,:),  ALLOCATABLE :: marker2dID

    INTEGER :: nProbeBody
    INTEGER, DIMENSION(:,:),  ALLOCATABLE :: markerProbe

    INTEGER(KIND=INT_K), DIMENSION(:,:,:), ALLOCATABLE :: igmark,iblank,dead_cell, &
                                                          bodyNum, ghostCellMark, iblankTemp, &
                                                          iblankUndecided
    
    !Ye,test
    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE  :: iblank_real

    !Ye, pure solid node for DC
    INTEGER(KIND=INT_K), DIMENSION(:,:,:), ALLOCATABLE :: ps1, ps2

    !Ye,goa
    INTEGER(KIND=INT_K), DIMENSION(:,:,:), ALLOCATABLE :: goa
    REAL(KIND=CGREAL), DIMENSION(:,:,:,:), ALLOCATABLE :: xyzDCIP

    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ghostCellIndex

    REAL(KIND=CGREAL), DIMENSION(:,:,:), ALLOCATABLE  :: xBItable, yBItable, zBItable
    INTEGER,           DIMENSION(:,:,:), ALLOCATABLE  :: cELtable
    INTEGER(KIND=INT_K),DIMENSION(:,:,:),ALLOCATABLE  :: is_BIset
    !INTEGER,           DIMENSION(:),ALLOCATABLE       :: NeighElemInd
    !REAL(KIND=CGREAL), DIMENSION(:),ALLOCATABLE       :: distMarker
    
    REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: xBodyMarker,yBodyMarker,zBodyMarker,    &
                                                        uBodyMarker,vBodyMarker,wBodyMarker,    &
                                                        axBodyMarker,ayBodyMarker,azBodyMarker, &
                                                        sBodyMarker,dsBodyMarker,               &
                                                        xNormBodyMarker,                        &
                                                        yNormBodyMarker,zNormBodyMarker,        &
                                                        dpdnBodyMarker

    !Ye,monitor the true BM velocity before filtered from the solid
    REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: u1BodyMarker,v1BodyMarker,w1BodyMarker

    REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: pgradx1,pgradx2
    REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: pgrady1,pgrady2 
    REAL(KIND=CGREAL), DIMENSION(:,:),   ALLOCATABLE :: pgradz1,pgradz2 

    REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: xBodyMarkerOld,yBodyMarkerOld,zBodyMarkerOld

    REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: uBodyMarkerOld,vBodyMarkerOld,wBodyMarkerOld, &
                                                        uBodyMarkerTmp,vBodyMarkerTmp,wBodyMarkerTmp

    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: triElemNormX,triElemNormY, &
                                                      triElemNormZ,triElemArea
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: triElemCentX,triElemCentY, &
                                                      triElemCentZ
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: triElemTang1X,triElemTang1Y, &
                                                      triElemTang1Z
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: triElemTang2X,triElemTang2Y, &
                                                      triElemTang2Z
    REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: pointOutsideBodyX,pointOutsideBodyY, &
                                                      pointOutsideBodyZ,surfArea

    INTEGER, DIMENSION(:,:,:),         ALLOCATABLE :: triElemNeig
    INTEGER, DIMENSION(:),             ALLOCATABLE :: totNumTriElem
    INTEGER, DIMENSION(:,:),           ALLOCATABLE :: Flag_outside_Marker

    REAL (KIND=CGREAL)                             :: normDirFlag

    !INTEGER, DIMENSION(:,:), ALLOCATABLE :: iCellIndexM,jCellIndexM,kCellIndexM
    !REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: xTang1BodyMarker, yTang1BodyMarker,     &
    !                                                    zTang1BodyMarker,                       &
    !                                                    xTang2BodyMarker, yTang2BodyMarker,     &
    !                                                    zTang2BodyMarker

    REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: MarkerArea
    REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: xMarkerForce,yMarkerForce,zMarkerForce, &
                                                        xMarkerShear,yMarkerShear,zMarkerShear, &
                                                        xMarkerStress,yMarkerStress,zMarkerStress
    REAL(KIND=CGREAL), DIMENSION(:)    , ALLOCATABLE :: xyzMarkerForceTmp,xyzMarkerForceTmp_i

    REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: xMarkerForceOld,yMarkerForceOld, &
                                                        zMarkerForceOld
    
    !Ye, store the marker force before filtering
    REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: xMarkerForce0,yMarkerForce0,zMarkerForce0

   END MODULE boundary_arrays
!------------------------------------------------------

   MODULE solver_arrays
   
    USE global_parameters
    
    IMPLICIT NONE
   
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: amx,apx,acx
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: amy,apy,acy
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: amz,apz,acz
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: rhs,dummy
    !REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: face1, face2

    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: amx_ad,apx_ad
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: amy_ad,apy_ad
    REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE :: amz_ad,apz_ad  
     
   END MODULE solver_arrays
!------------------------------------------------------

   MODULE GCM_arrays
   
    USE global_parameters
   
    IMPLICIT NONE
   
    INTEGER                                :: nGhost
    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: incI, incJ, incK, iPvt
    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: closestMarker,            &
                                              closestElementGC,         &
                                              iGhost,jGhost,kGhost,     &
                                              iCellIndex,jCellIndex,kCellIndex
                                              
    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: iDead,jDead,kDead      &
                                             ,iDeadCellIndex,jDeadCellIndex,kDeadCellIndex
                                             
    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: closestMarkerDead
    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: closestElementDead

    REAL(KIND=CGREAL), DIMENSION(2)                :: det  
    REAL(KIND=CGREAL), DIMENSION(:)  , ALLOCATABLE :: work
    REAL(KIND=CGREAL), DIMENSION(:)  , ALLOCATABLE :: closestMarkerRatio, &
                                                      probeLength,        &
                                                      xBodyInterceptTang, &
                                                      yBodyInterceptTang, &
                                                      zBodyInterceptTang, &
                                                      xBodyInterceptNorm, &
                                                      yBodyInterceptNorm, &
                                                      zBodyInterceptNorm, &
                                                      xBodyIntercept,     &
                                                      yBodyIntercept,     &
                                                      zBodyIntercept,     &
                                                      uBodyIntercept,     &
                                                      vBodyIntercept,     &
                                                      wBodyIntercept,     &
                                                      pBodyIntercept,     &
                                                      dpdnBodyIntercept,  &
                                                      dpdtBodyIntercept,  &
                                                      xImagePoint,        &
                                                      yImagePoint,        &
                                                      zImagePoint,        &
                                                      res_laplacian,      &
                                                      alphaGhost

    REAL(KIND=CGREAL), DIMENSION(:),   ALLOCATABLE :: closestMarkerRatioDead, &
                                                      xBodyInterceptDead,     &
                                                      yBodyInterceptDead,     &
                                                      zBodyInterceptDead,     &
                                                      uBodyInterceptDead,     &
                                                      vBodyInterceptDead,     &
                                                      wBodyInterceptDead,     &   
                                                      xNormInterceptDead,     &
                                                      yNormInterceptDead,     &
                                                      zNormInterceptDead,     &
                                                       dpdnInterceptDead,     &
                                                         probeLengthDead,     &
                                                         alphaDead
    
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: ghostNScoeff

    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: coeffGCMD, coeffGCMN,   &
                                                      coeffGCMDeadD,          &
                                                      coeffGCMDeadN
    
    REAL(KIND=CGREAL), DIMENSION(:)  , ALLOCATABLE :: coeffGCMDirc, coeffGCMNeum
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: vanMatrixD, vanMatrixN
                                                              
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: xBodyCentroid,yBodyCentroid,  &
                                                      sBodyCentroid,                &
                                                      xCentroidTang,yCentroidTang, &
                                                      xCentroidNorm,yCentroidNorm

    INTEGER, DIMENSION(:)  , ALLOCATABLE   :: intI, intJ, intK

   END MODULE GCM_arrays
   
!------------------------------------------------------
   MODULE stat_arrays

    USE global_parameters

    IMPLICIT NONE

    INTEGER                                        :: statCtr
    REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: uAv,vAv,wAv,pAv
    REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: uvAv,vwAv,uwAv
    REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: uuAv,vvAv,wwAv

   END MODULE stat_arrays

!------------------------------------------------------
   MODULE stat_vort_arrays
 
    USE global_parameters
 
    IMPLICIT NONE
 
    INTEGER                                        :: statCtrv
    REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: oxAv,oyAv,ozAv
    REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:,:,:) :: oxoxAv,oyoyAv,ozozAv
 
   END MODULE stat_vort_arrays
!------------------------------------------------------

   MODULE blasius_profile

    USE global_parameters

    IMPLICIT NONE

    REAL(KIND=CGREAL)                          :: cavity_H,slot_H,ddratio,d,delta,uinf
    REAL(KIND=CGREAL),ALLOCATABLE,DIMENSION(:) :: eta,u_blasius
    INTEGER                                    :: l,junk,i_start

   END MODULE blasius_profile
!------------------------------------------------------
   MODULE mg_module   ! Multigrid module
 
    USE global_parameters
    USE flow_parameters
 
    IMPLICIT NONE
     
    INTEGER, PARAMETER :: max_grid_level = 10
    INTEGER  :: mgrids_x, mgrids_y, mgrids_z, mLevel
    INTEGER  :: mgcyclex, mgcycley, mgcyclez, infoconv, incrlev
    INTEGER  :: iterFinest, iterInter, iterCoarest
    INTEGER  :: ittt1

    INTEGER :: iRedBlack, TNcolorX, TNcolorY, TNcolorZ, iStep, jStep, kStep  
               !new for Redblack LSOR

    REAL(KIND=CGREAL) :: pmgx1, pmgx2, pmgy1, pmgy2, pmgz1, pmgz2

    INTEGER, DIMENSION(:), ALLOCATABLE :: mgrid_I, mgrid_J, mgrid_K
    INTEGER, DIMENSION(:), ALLOCATABLE :: zc_start_mg, zc_end_mg, z_start_mg, z_end_mg, &
                                          zb1_mg, zb2_mg, num_slicesZ_mg
    INTEGER, DIMENSION(:), ALLOCATABLE :: l_proc_mg, r_proc_mg, leftmost_mg, rightmost_mg

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: all_slicesY_mg, all_slicesZ_mg

    TYPE mg_array_type
       REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: dxcinv,dycinv,dzcinv
       REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: dxinv,dyinv,dzinv
       REAL(KIND=CGREAL),DIMENSION(:),ALLOCATABLE :: x, xc, y, yc, &
                                                     z, zc

       REAL(KIND=CGREAL),  DIMENSION(:,:,:),ALLOCATABLE :: rhs,phi,res
       INTEGER(KIND=INT_K),DIMENSION(:,:,:),ALLOCATABLE :: ibk,igk
    ENDTYPE mg_array_type

    TYPE(mg_array_type),DIMENSION(:),ALLOCATABLE :: mg_array

    INTEGER, DIMENSION(:), ALLOCATABLE :: yc_start_mg, yc_end_mg, y_start_mg, y_end_mg, &
                                          yb1_mg, yb2_mg, num_slicesY_mg

    INTEGER, DIMENSION(:), ALLOCATABLE :: b_proc_mg, t_proc_mg, btommost_mg, topmost_mg

   END MODULE mg_module
!------------------------------------------------------
   MODULE usr_module
    
    USE global_parameters
     
    REAL(KIND=CGREAL) :: density_ratio
    REAL(KIND=CGREAL) :: lScale,vScale
    REAL(KIND=CGREAL) :: non_dim_volume,volume

    INTEGER           :: ifuMarkerTrace
    REAL(KIND=CGREAL) :: solidGap, minSolidGap, xflux_rate

    !!! Added by Paulo Ferreira de Sousa on 02/27/2009
    REAL(KIND=CGREAL), DIMENSION(:,:)  , ALLOCATABLE :: xBMPosition,yBMPosition,zBMPosition,    &
                                                        uBMPosition,vBMPosition,wBMPosition

    ! eigenmode specification
    REAL(KIND=CGREAL),DIMENSION(:,:),ALLOCATABLE  :: eigenmodeX, &
                                                     eigenmodeY, &
                                                     eigenmodeZ, &
                                                     xBodyMarkerOrg, &
                                                     yBodyMarkerOrg, &
                                                     zBodyMarkerOrg

    END MODULE usr_module 

!------------------------------------------------------
   MODULE fsi_module      ! Flow/structure interaction module.

     USE global_parameters  
     IMPLICIT NONE

     ! iteration parameters for FSI
     INTEGER  :: itermax_FSI, iter_FSI,  is_FSI_cvg, nsub

     !Ye,0:fsi not activated; 1:fsi activated
     INTEGER  :: fsi_on

     ! torlerance on the residuals of pressure, velocity and displacement
     ! of the boundary
     REAL(KIND=CGREAL) :: restol_FSIp, restol_FSIv, restol_FSIdsp
     REAL(KIND=CGREAL) :: resMax_FSIp, resMax_FSIv, resMax_FSIdsp

     !Ye, adaptive relaxtion factor for velocity
     REAL(KIND=CGREAL) :: resMax_FSIv_old, relax_FSIp0, relax_FSIv0

     !Ye, adaptive relaxtion factor for force
     REAL(KIND=CGREAL) :: resMax_FSIp_old

     ! relaxation factors for boundary pressure, velocity, and displacement
     REAL(KIND=CGREAL) :: relax_FSIp, relax_FSIv, relax_FSIdsp

     REAL(KIND=CGREAL) :: FSI_delay          ! time delay constant for gradually enabling FSI
     REAL(KIND=CGREAL) :: fsi_filter_f, fsi_filter_f0, fsi_filter_v  ! fsi filter constant
  
     INTEGER           :: ifuElasId
     INTEGER, DIMENSION(:),   ALLOCATABLE :: bodyInFlowID
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: elasticMarkerID

     ! parameters in the solid solver
     INTEGER           :: nPtsMax_s, nElasticBody, ndata_s
     INTEGER           :: nread_s, nrestart_s, nstep_s
     REAL(KIND=CGREAL) :: dt_s


     ! arrays to store data from the solid solver
     REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: xdata_s, ydata_s, zdata_s

     ! stress multiplied by the element area
     !REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE :: xf_elem, yf_elem, zf_elem
     REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE  :: cnt_fx,cnt_fy,cnt_fz
     REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE  :: cnt_force
     REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE  :: cfx,cfy,cfz

     REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE  :: cnt_d12,cnt_d21,cnt_d13,  &
                                                        cnt_d31,cnt_d23,cnt_d32
     
     REAL(KIND=CGREAL), DIMENSION(:), ALLOCATABLE  :: cntd12,cntd21,cntd13,  &
                                                      cntd31,cntd23,cntd32

     INTEGER, DIMENSION(:,:), ALLOCATABLE :: cnt_EM12,cnt_EM21,cnt_EM13,  &
                                             cnt_EM31,cnt_EM23,cnt_EM32

     INTEGER, DIMENSION(:), ALLOCATABLE :: ElemNum12,ElemNum21,ElemNum13,  &
                                           ElemNum31,ElemNum23,ElemNum32

     !Ye, test MarkerForce Interpolation
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: dotP, dotP2

     !INTEGER, DIMENSION(:), ALLOCATABLE :: v_num
      
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: set_old 

     INTEGER, DIMENSION(:,:), ALLOCATABLE :: cnt_fg12, cntflg12
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: cnt_fg21, cntflg21
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: cnt_fg13, cntflg13
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: cnt_fg31, cntflg31
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: cnt_fg23, cntflg23
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: cnt_fg32, cntflg32

   END MODULE fsi_module
!------------------------------------------------------
   MODULE MPI_module

    use global_parameters

    IMPLICIT NONE

    include "mpif.h"

    INTEGER  ::  nProcs_g, nProcs, lProc_g, lProc
    INTEGER  ::  ierr
    integer  ::  istatus(MPI_STATUS_SIZE)
    integer  ::  PROC_M, proc_s(8), nElastic
    INTEGER  ::  WORLD_GROUP, FLOW_GROUP, FLOW_COMM

    INTEGER  ::  z_slab_unif, kSlices, zc_start, zc_end, z_start, z_end, &
                 zb1, zb2
    INTEGER  ::  lProc_leftmost, lProc_rightmost
    INTEGER  ::  lProc_left, lProc_right

    INTEGER, DIMENSION(:), ALLOCATABLE :: numSlices
    
    ! z_slab_unif
    INTEGER, PARAMETER :: UNIFORM_SLAB             = 1, &
                          NONUNIFORM_SLAB          = 2

    integer :: toleft=1, torght=2
    integer :: isend_rq_toleft,isend_rq_torght
    integer :: irecv_rq_fromleft,irecv_rq_fromrght
    integer, dimension(MPI_STATUS_SIZE) :: isend_stat_tl,isend_stat_tr
    integer, dimension(MPI_STATUS_SIZE) :: irecv_stat_fl,irecv_stat_fr

    INTEGER(KIND=INT_K),DIMENSION(:,:,:),ALLOCATABLE :: i_lbuff, i_rbuff
    
    !Added by Ye
    INTEGER  :: y_slab_unif,jProc,kProc,jSlices
    INTEGER, DIMENSION(:), ALLOCATABLE :: numSlicesY,numSlicesZ
    INTEGER  ::  yc_start, yc_end
    INTEGER  ::  y_start, y_end
    INTEGER  ::  yb1, yb2
    INTEGER :: nProcY
    INTEGER :: nProcZ
    INTEGER  ::  lProc_btommost, lProc_topmost
    INTEGER  ::  lProc_top, lProc_bottom
    integer :: tobtom=3, totop=4
    integer :: isend_rq_tobtom,isend_rq_totop
    integer :: irecv_rq_frombtom,irecv_rq_fromtop
    integer, dimension(MPI_STATUS_SIZE) :: isend_stat_tb,isend_stat_tt
    integer, dimension(MPI_STATUS_SIZE) :: irecv_stat_fb,irecv_stat_ft

   END MODULE MPI_module
!------------------------------------------------------
