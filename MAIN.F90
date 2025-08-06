!---------------------------------------
!  SUBROUTINE open_files()
!  SUBROUTINE read_inputs()
!  SUBROUTINE check_inputs()
!  SUBROUTINE read_inputs_bound()
!  SUBROUTINE check_inputs_bound()
!  SUBROUTINE init_simulation()
!  SUBROUTINE stop_program()
!  SUBROUTINE init_flow()
!  SUBROUTINE read_restart_flow()
!  SUBROUTINE read_restart_body()
!  SUBROUTINE write_restart()
!  SUBROUTINE set_arrays_zero()
!  SUBROUTINE initialize_flowfield()
!---------------------------------------
!
!-------------------------------------------------
! Code to simulate 3-D incompressible, unsteady, viscous 
! flow in domain with multiple, complex embedded obstacles
! including zero-thickness membranes 
! 
! Only Cartesian grids (uniform or non-uniform) are allowed
!
!
! spatial discretization -  2nd order central difference
! temporal discretization-  Crank-Nicoloson for non-linear term and viscous terms 
!
!     |-------|-------|-------|-------|-------
!     |       |       |       |       |
!     |   o   |   o   |   o   |   o   |   o
!     |       |       |       |       |
!     |-------|-------|-------|-------|-------
!     |       |*******|*******|*******|
!     |   o   +***1***|***1***|***1***|   o
!     |       |*******|*******|*******|
!     |-------|-------|-------|-------|-------
!     |       |*******|*******|       |
!     |   o   |***1***|***1***|   o   |   o
!     |       |*******|*******|       |
!     |-------|-------|-------|-------|-------
!     |       |*******|*******|       |
!     |   o   |***1***|***1***|   o   |   o
!     |       |*******|*******|       |
!     |-------|-------|-------|-------|-------
!     |       |       |       |       |
!     |   o   |   o   |   o   |   o   |   o
!     |       |       |       |       |
!     |-------|---+---|-------|-------|-------
!

!---------------------------------------------
! Main Program
!---------------------------------------------

   PROGRAM FEICAR3D

   USE global_parameters
   USE flow_parameters
   USE MPI_module
   USE boundary_arrays

    implicit none
    integer       ::i
!***********************

     CALL MPI_Initialize()
     
     !if(lProc_g .lt. nProcs_g-nElastic) then    ! if the processor is for flow
     if (lProc_g .gt. nElastic-1) then             ! for flow processes

       if(lProc .eq. proc_m) PRINT*,'ENTERING INIT_PROGRAM'
       CALL init_program()

       if(lProc .eq. proc_m) PRINT*,'ENTERING READ_INPUTS'
       CALL read_inputs()  !read from 'input.dat'
       CALL slice_assignment()
 
       if(lProc.eq.proc_m) then
          PRINT*,'ENTERING CHECK_INPUTS'          
          CALL check_inputs()
          CALL check_inputs_body()

          !PRINT*,'ENTERING TURB_CHECK_INPUTS'
          !CALL turb_check_inputs() 
       endif

       if(lProc.eq.proc_m) then
         PRINT*,'ENTERING init_simulation'
       endif
       !PRINT*,'ENTERING INIT_SIMULATION', lProc
       CALL init_simulation()

       if(lProc .eq. PROC_M) CALL initialize_fsi()
 
       !Ye,segmentation test
       !call compute_marker_stress()
       !if (lProc.eq.proc_m)then
       !   do i=1,10482
       !      write(400,*)xMarkerforce(i,2),yMarkerforce(i,2),zMarkerforce(i,2)
       !   enddo
       !endif
       !call mpi_barrier(MPI_COMM_WORLD,ierr)
       !stop

       if(lProc.eq.proc_m) PRINT*,'ENTERING TIME_STEP'

       CALL time_step_viscous()

       CALL stop_program()

    else  ! lProc is for the solid solver
       print*,'lProc_g',lProc_g,' calling the solid solver...'
       CALL nonstad_main(lProc_g)
    endif  ! lProc
    
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call MPI_Finalize(ierr)
 
    STOP   
    
  END PROGRAM FEICAR3D

!---------------------------------------------
! File numbers 50 to 100 reserved for input/output files
!---------------------------------------------

   SUBROUTINE init_program()

    USE global_parameters
    USE flow_parameters 

    IMPLICIT NONE

    pi = 4.0_CGREAL*ATAN(1.0_CGREAL)

    ! set File Units numbers
 
    ifuInput         = 50
    ifuIblankIn      = 52
    ifuRstrtFlowIn   = 53
    ifuRstrtFlowOut  = 54
    ifuRstrtBodyIn   = 55
    ifuRstrtBodyOut  = 56
    ifuBodyIn        = 57
    ifuProbeIn       = 58
    ifuMarkerIn      = 59
    ifuUnstrucSurfIn = 60 
    ifuUnstrucSurfOut= 61 
    ifuBodyOut       = 62

    ifuDragOut      = 100
    ifuFreshCellOut = 150
    ifuMarkerOut    = 160
    ifuStatOut      = 170
    ifuStatPlot     = 171
    ifuProbeOut     = 180

    ! Open files
    !OPEN(ifuIblankIn,      FILE='iblank_in.dat')
    !OPEN(ifuFreshCellOut,  FILE='freshcell_out.dat',     STATUS='UNKNOWN')
    !OPEN(ifuStatOut,       FILE='stat_out.dat'          ,FORM='UNFORMATTED') 

   END SUBROUTINE init_program
!---------------------------------------------
   SUBROUTINE read_inputs()

    USE global_parameters
    USE flow_parameters
    USE mg_module
    USE fsi_module
    USE MPI_module

    IMPLICIT NONE

    INTEGER :: n

    OPEN(ifuInput,         FILE='input.dat',             STATUS='OLD')

    READ(ifuInput,*)
    READ(ifuInput,*) nread             
    READ(ifuInput,*) 
    READ(ifuInput,*) ndim
    READ(ifuInput,*)
    READ(ifuInput,*) nx,ny,nz, y_slab_unif, z_slab_unif
    READ(ifuInput,*)
    READ(ifuInput,*) xgrid_unif,xOrigin,xout
    READ(ifuInput,*)
    READ(ifuInput,*) ygrid_unif,yOrigin,yout
    READ(ifuInput,*)
    READ(ifuInput,*) zgrid_unif,zOrigin,zout
    READ(ifuInput,*)
    READ(ifuInput,*) uinit,vinit,winit, vper !new change for 2D-3D perturbations
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) bcx1,ux1,vx1,wx1,freq_ux1,freq_vx1,freq_wx1
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) bcx2,ux2,vx2,wx2,freq_ux2,freq_vx2,freq_wx2
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) bcy1,uy1,vy1,wy1,freq_uy1,freq_vy1,freq_wy1
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) bcy2,uy2,vy2,wy2,freq_uy2,freq_vy2,freq_wy2
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) bcz1,uz1,vz1,wz1,freq_uz1,freq_vz1,freq_wz1
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) bcz2,uz2,vz2,wz2,freq_uz2,freq_vz2,freq_wz2
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) pbcx1, pbcx2, pbcy1, pbcy2, pbcz1, pbcz2
    READ(ifuInput,*)
    READ(ifuInput,*) pppx1, pppx2, pppy1, pppy2, pppz1, pppz2
    READ(ifuInput,*)
    READ(ifuInput,*)            
    READ(ifuInput,*) no_tsteps,nmonitor,ndump,nrestart,nstat,nmonitor_probe_liftdrag
    READ(ifuInput,*)
    READ(ifuInput,*) Imonitor_drag,Imonitor_probe,format_dump
    READ(ifuInput,*)
    READ(ifuInput,*) re,dt
    READ(ifuInput,*)
    READ(ifuInput,*) frac_step_type, advec_scheme
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) alfa_upwind, x1_upwind
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) internal_boundary_present
    READ(ifuInput,*)
    READ(ifuInput,*) probeLengthNormalized, itermax_gc, restol_gcv, restol_gcp 
    READ(ifuInput,*)
    READ(ifuInput,*) mix_GC_form, probeLengthNormalizedD, membrane_tkns_factor
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) itermax_ad,restol_ad,omega_ad
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) it_solver_type, iRedBlack
    READ(ifuInput,*)
    READ(ifuInput,*) omega_pson
    READ(ifuInput,*)
    READ(ifuInput,*) itermax_pson,restol_pson
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) mgcyclex, mgcycley, mgcyclez
    READ(ifuInput,*)
    READ(ifuInput,*) iterFinest, iterInter, iterCoarest
    READ(ifuInput,*)
    READ(ifuInput,*) mgrids_x, mgrids_y, mgrids_z 
    READ(ifuInput,*)
    READ(ifuInput,*) infoconv
    READ(ifuInput,*)
    READ(ifuInput,*)
    READ(ifuInput,*) itermax_FSI, FSI_delay, fsi_filter_f0, fsi_filter_v
    READ(ifuInput,*)
    READ(ifuInput,*) restol_FSIp, restol_FSIv, restol_FSIdsp
    READ(ifuInput,*)
    READ(ifuInput,*) relax_FSIp, relax_FSIv, relax_FSIdsp

    close(ifuInput)

    !Ye,if there is no solid body involved, need this line to define density since
    !density_fluid will not be read from canonical_body_in.dat
    density_fluid = 1.0d0
 
    ! scale the boundary pressure by the fluid density as p/rho is solved in VICAR3d.
    ! This applies only to cases where pressure Dirichlet B.C. is present.
    pppx1 = pppx1/density_fluid
    pppx2 = pppx2/density_fluid
    pppy1 = pppy1/density_fluid
    pppy2 = pppy2/density_fluid
    pppz1 = pppz1/density_fluid
    pppz2 = pppz2/density_fluid

    relax_FSIv0 = relax_FSIv
    relax_FSIp0 = relax_FSIp

    ! print input
    ! if(lProc .eq. PROC_M) CALL print_inputs()
    
    ! new for Red-Black LSOR
    IF (iRedBlack == 0) THEN
         TNcolorX = 1
         TNcolorY = 1
         TNcolorZ = 1
         iStep = 1
         jStep = 1
         kStep = 1
    ELSE IF(iRedBlack == 1) THEN
         TNcolorX = 2
         TNcolorY = 1
         TNcolorZ = 1

         iStep = 1
         jStep = 2
         kStep = (ndim-DIM_2D) + 1
    ELSE 
         TNcolorX = 2
         TNcolorY = 2
         TNcolorZ = 2

         iStep = 2
         jStep = 2
         kStep = (ndim-DIM_2D) + 1 
    END IF

   END SUBROUTINE read_inputs
!---------------------------------------------

!---------------------------------------------
   SUBROUTINE print_inputs()

    USE global_parameters
    USE flow_parameters

    IMPLICIT NONE

    PRINT*,' NREAD in Input = ',nread
    PRINT*,' NDIM in Input  = ',ndim
        
        PRINT*,' NX in Input = ',nx
    PRINT*,' NY in Input = ',ny
    PRINT*,' NZ in Input = ',nz
 
    PRINT*,' XGRID_UNIF in Input = ',xgrid_unif
    PRINT*,' YGRID_UNIF in Input = ',ygrid_unif
    PRINT*,' ZGRID_UNIF in Input = ',zgrid_unif

    PRINT*,' advec_scheme   in Input = ',advec_scheme
    PRINT*,' IT_SOLVER_TYPE in Input = ',it_solver_type

    PRINT*,' advec_scheme   in Input = ',advec_scheme
    PRINT*,' Itermax_ad     in Input = ',itermax_ad
    PRINT*,' Itermax_pson   in Input = ',itermax_pson
   
    PRINT*,' restol_ad       in Input = ',restol_ad
    PRINT*,' restol_pson     in Input = ',restol_pson

    PRINT*,' nGhostMax       in Input = ',nGhostMax

   END SUBROUTINE print_inputs
!---------------------------------------------
   SUBROUTINE check_inputs()

    USE global_parameters
    USE flow_parameters
    !USE MG_parameters

    IMPLICIT NONE
    
    INTEGER :: iBody,n

    IF ( nread > 1 .OR. nread < 0 ) THEN
      PRINT*,'Incorrect Value for NREAD in Input',nread
      PRINT*,'Either use 0 or 1'
      STOP
    END IF ! nread
        
    IF ( ndim > DIM_3D .OR. ndim < DIM_2D ) THEN
      PRINT*,'Incorrect Value for NDIM in Input',ndim
      PRINT*,'Either use 2 or 3'
      STOP
    END IF ! ndim
        
    IF ( nx < 3 ) THEN
      PRINT*,'Incorrect Value for NX in Input',nx
      PRINT*,'Use at Least 3'
      STOP
    END IF ! nx
    
    IF ( ny < 3 ) THEN
      PRINT*,'Incorrect Value for NY in Input',ny
      PRINT*,'Use at Least 3'
      STOP
    END IF ! ny
    
    IF ( nz < 3 ) THEN
      PRINT*,'Incorrect Value for NZ in Input',nz
      PRINT*,'Use at Least 3'
      STOP
    END IF ! nz

    IF ( xgrid_unif < UNIFORM_GRID .AND. xgrid_unif > NONUNIFORM_GRID ) THEN
      PRINT*,'Incorrect Value for XGRID_UNIF in Input',xgrid_unif
      PRINT*,'Use either 1 (Uniform) or 2 (Non-Uniform)'
      STOP
    END IF ! xgrid_unif

    IF ( ygrid_unif < UNIFORM_GRID .AND. ygrid_unif > NONUNIFORM_GRID ) THEN
      PRINT*,'Incorrect Value for YGRID_UNIF in Input',ygrid_unif
      PRINT*,'Use either 1 (Uniform) or 2 (Non-Uniform)'
      STOP
    END IF ! ygrid_unif

    IF ( zgrid_unif < UNIFORM_GRID .AND. zgrid_unif > NONUNIFORM_GRID ) THEN
      PRINT*,'Incorrect Value for ZGRID_UNIF in Input',zgrid_unif
      PRINT*,'Use either 1 (Uniform) or 2 (Non-Uniform)'
      STOP
    END IF ! zgrid_unif

    IF ( nDim == DIM_2D .AND. zgrid_unif /= UNIFORM_GRID ) THEN
      PRINT*,'Incorrect Value for ZGRID_UNIF in Input',zgrid_unif
      PRINT*,'Use a value of 1 (Uniform) for 2D simulations'
      STOP
    END IF ! nDim

    IF ( it_solver_type <= 0 .OR. it_solver_type > 4 ) THEN
      PRINT*,'Incorrect Value for IT_SOLVER_TYPE in Input',it_solver_type
      PRINT*,'Either use 1 (LSOR) or 2 (PETSC) or 3 (MG)'
      STOP
    END IF ! it_solver_type

   END SUBROUTINE check_inputs
!---------------------------------------------
!---------------------------------------------
   SUBROUTINE read_inputs_body()

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE MPI_module

    IMPLICIT NONE


    REAL(KIND=CGREAL) :: readDummyR, Aratio
    INTEGER           :: n,readDummyInt

    !IF ( internal_boundary_present == INTR_BOUND_NONE) GOTO 999 

    OPEN(ifuBodyIn,FILE='canonical_body_in.dat',STATUS='UNKNOWN')
    READ(ifuBodyIn,*)nbody
    !PRINT*,'   nBody = ',nBody

   ! allocate memory -------------------------------------------------------------

    ALLOCATE(canonical_body_type(nBody))
    ALLOCATE(body_dim(nBody))
    ALLOCATE(boundary_motion_type(nBody))
    ALLOCATE(unstruc_surface_type(nBody))
    ALLOCATE(wall_type(nBody))
    
    ALLOCATE(nPtsBodyMarkerOrig(nbody))
    ALLOCATE(nPtsBodyMarker(nbody))
    ALLOCATE(totNumTriElem(nBody))
    
    ALLOCATE(n_phi(nbody))
    ALLOCATE(n_theta(nbody))
    
    ALLOCATE(radiusx(nbody))
    ALLOCATE(radiusy(nbody))
    ALLOCATE(radiusz(nbody))

    ALLOCATE(xcent(nbody))
    ALLOCATE(ycent(nbody))
    ALLOCATE(zcent(nbody))
    ALLOCATE(xcentinit(nbody))
    ALLOCATE(ycentinit(nbody))
    ALLOCATE(zcentinit(nbody))

    ALLOCATE(vxcent(nbody))
    ALLOCATE(vycent(nbody))
    ALLOCATE(vzcent(nbody))

    ALLOCATE(angvx(nbody))
    ALLOCATE(angvy(nbody))
    ALLOCATE(angvz(nbody))

    ALLOCATE(angvx_old(nbody))
    ALLOCATE(angvy_old(nbody))
    ALLOCATE(angvz_old(nbody))

    ALLOCATE(vxcentTrans(nbody))
    ALLOCATE(vycentTrans(nbody))
    ALLOCATE(vzcentTrans(nbody))
    ALLOCATE(ampx(nbody))
    ALLOCATE(ampy(nbody))
    ALLOCATE(ampz(nbody))
    ALLOCATE(freqx(nbody))
    ALLOCATE(freqy(nbody))
    ALLOCATE(freqz(nbody))
       
    ALLOCATE(vxphase(nbody))
    ALLOCATE(vyphase(nbody))
    ALLOCATE(vzphase(nbody))

    ALLOCATE(angvxinit(nbody))
    ALLOCATE(angvyinit(nbody))
    ALLOCATE(angvzinit(nbody))
    
    ALLOCATE(angx(nbody))
    ALLOCATE(angy(nbody))
    ALLOCATE(angz(nbody))

    ALLOCATE(angxinit(nbody))
    ALLOCATE(angyinit(nbody))
    ALLOCATE(angzinit(nbody))
    
    ALLOCATE(xPivot(nbody))
    ALLOCATE(yPivot(nbody))
    ALLOCATE(zPivot(nbody))

    ALLOCATE(ampangx(nbody))
    ALLOCATE(ampangy(nbody))
    ALLOCATE(ampangz(nbody))
    ALLOCATE(freqangx(nbody))
    ALLOCATE(freqangy(nbody))
    ALLOCATE(freqangz(nbody))

    ALLOCATE(angxphase(nbody))
    ALLOCATE(angyphase(nbody))
    ALLOCATE(angzphase(nbody))

    ALLOCATE(xcentConstr(nBody),ycentConstr(nBody),zcentConstr(nBody))
    ALLOCATE(forced_motion_spec(nbody))
    ALLOCATE(density_solid(nBody),zoom_factor(nBody))

   ! initialize values and arrays ------------------------------------------------
 
    nPtsBodyMarker = 0

    n_phi       = 0
    n_theta     = 0
    
    radiusx     = 0.0_CGREAL
    radiusy     = 0.0_CGREAL
    radiusz     = 0.0_CGREAL

    xcent       = 0.0_CGREAL 
    ycent       = 0.0_CGREAL 
    zcent       = 0.0_CGREAL 

    vxcent      = 0.0_CGREAL 
    vycent      = 0.0_CGREAL 
    vzcent      = 0.0_CGREAL 

    angvx       = 0.0_CGREAL 
    angvy       = 0.0_CGREAL 
    angvz       = 0.0_CGREAL 

    !new arrays for rotation
    angvx_old   = 0.0_CGREAL
    angvy_old   = 0.0_CGREAL
    angvz_old   = 0.0_CGREAL

    xcentinit   = 0.0_CGREAL 
    ycentinit   = 0.0_CGREAL 
    zcentinit   = 0.0_CGREAL 

    vxcentTrans = 0.0_CGREAL 
    vycentTrans = 0.0_CGREAL 
    vzcentTrans = 0.0_CGREAL 

    ampx        = 0.0_CGREAL 
    ampy        = 0.0_CGREAL
    ampz        = 0.0_CGREAL

    freqx       = 0.0_CGREAL 
    freqy       = 0.0_CGREAL
    freqz       = 0.0_CGREAL

    angvxinit   = 0.0_CGREAL 
    angvyinit   = 0.0_CGREAL 
    angvzinit   = 0.0_CGREAL 

    angxinit    = 0.0_CGREAL 
    angyinit    = 0.0_CGREAL 
    angzinit    = 0.0_CGREAL 

    ampangx     = 0.0_CGREAL 
    ampangy     = 0.0_CGREAL
    ampangz     = 0.0_CGREAL

    freqangx    = 0.0_CGREAL 
    freqangy    = 0.0_CGREAL
    freqangz    = 0.0_CGREAL        

   ! read input file pertinent to internal boundary ------------------------------

    boundary_motion = FIXED_BOUNDARY
    elastic_present = 0
    nElastic = 0

    if(lProc_g .eq. 0) PRINT*,'Reading canonical_body_in.dat ...'
    DO n=1,nbody
        READ(ifuBodyIn,*)
        READ(ifuBodyIn,*)
        READ(ifuBodyIn,*)canonical_body_type(n), body_dim(n)
        READ(ifuBodyIn,*)unstruc_surface_type(n),boundary_motion_type(n)        
        READ(ifuBodyIn,*)wall_type(n), forced_motion_spec(n)

        IF ( boundary_motion_type(n) /= STATIONARY )  boundary_motion = MOVING_BOUNDARY

        IF ( boundary_motion_type(n) == ELASTIC_MOTION ) THEN
           elastic_present = 1
           nElastic = nElastic + 1
        ENDIF

        SELECT CASE (canonical_body_type(n))

          CASE(ELLIPTIC_CYLINDER:GENERAL_CYLINDER)

             READ(ifuBodyIn,*)nPtsBodyMarker(n),readDummyInt
             READ(ifuBodyIn,*)radiusx(n),radiusy(n)
             READ(ifuBodyIn,*)xcent(n),ycent(n)
             READ(ifuBodyIn,*)angzinit(n)
             READ(ifuBodyIn,*)vxcentTrans(n),vycentTrans(n),vzcentTrans(n)
             READ(ifuBodyIn,*)vxphase(n),vyphase(n),vzphase(n)
             READ(ifuBodyIn,*)ampx(n),ampy(n),ampz(n)
             READ(ifuBodyIn,*)freqx(n),freqy(n),freqz(n)
             READ(ifuBodyIn,*)angvx(n),angvy(n),angvz(n)
             READ(ifuBodyIn,*)angzphase(n)            ! phase angle advance of angular velocity 
                                                      ! over translational velocity oscillation
             READ(ifuBodyIn,*)ampangx(n),ampangy(n),ampangz(n)
             READ(ifuBodyIn,*)freqangx(n),freqangy(n),freqangz(n)
             READ(ifuBodyIn,*)density_fluid, density_solid(n)
             READ(ifuBodyIn,*)xcentConstr(n),ycentConstr(n),zcentConstr(n)
            
          CASE(ELLIPSOID)

             READ(ifuBodyIn,*)n_theta(n),n_phi(n)  
             READ(ifuBodyIn,*)radiusx(n),radiusy(n),radiusz(n)
             READ(ifuBodyIn,*)xcent(n),ycent(n),zcent(n)
             READ(ifuBodyIn,*)angzinit(n)
             READ(ifuBodyIn,*)vxcentTrans(n),vycentTrans(n),vzcentTrans(n)
             READ(ifuBodyIn,*)vxphase(n),vyphase(n),vzphase(n)
             READ(ifuBodyIn,*)ampx(n),ampy(n),ampz(n)
             READ(ifuBodyIn,*)freqx(n),freqy(n),freqz(n)
             READ(ifuBodyIn,*)angvx(n),angvy(n),angvz(n)

             READ(ifuBodyIn,*)angzphase(n)
             READ(ifuBodyIn,*)ampangx(n),ampangy(n),ampangz(n)            
             READ(ifuBodyIn,*)freqangx(n),freqangy(n),freqangz(n)
             READ(ifuBodyIn,*)density_fluid, density_solid(n)
             READ(ifuBodyIn,*)xcentConstr(n),ycentConstr(n),zcentConstr(n)

             Aratio = radiusz(n)/radiusx(n)

             !IF ( MOD(n_theta(n),2) == 0 ) THEN
             !  n_phi(n) = IDNINT(0.5_CGREAL*REAL(n_theta(n),KIND=CGREAL)*Aratio)
             !ELSE
             !  n_phi(n) = IDNINT(0.5_CGREAL*REAL(n_theta(n)-1,KIND=CGREAL)*Aratio)+1
             !ENDIF ! n_theta

             nPtsBodyMarker(n) = n_theta(n)*n_phi(n)

           CASE(UNSTRUCTURED_SURFACE)

             ! Note: for unstruc surface, nPtsBodyMarker = total number of nodes             
         
             READ(ifuBodyIn,*)nPtsBodyMarker(n),totNumTriElem(n)
             READ(ifuBodyIn,*)xPivot(n),yPivot(n),zPivot(n)
             READ(ifuBodyIn,*)xcent(n),ycent(n),zcent(n)
             READ(ifuBodyIn,*)angxinit(n),angyinit(n),angzinit(n)
             READ(ifuBodyIn,*)vxcentTrans(n),vycentTrans(n),vzcentTrans(n)
             READ(ifuBodyIn,*)vxphase(n),vyphase(n),vzphase(n)
             READ(ifuBodyIn,*)ampx(n),ampy(n),ampz(n)
             READ(ifuBodyIn,*)freqx(n),freqy(n),freqz(n)
             READ(ifuBodyIn,*)angvx(n),angvy(n),angvz(n)

             READ(ifuBodyIn,*)angxphase(n),angyphase(n),angzphase(n)
             READ(ifuBodyIn,*)ampangx(n),ampangy(n),ampangz(n)
             READ(ifuBodyIn,*)freqangx(n),freqangy(n),freqangz(n)
             READ(ifuBodyIn,*)density_fluid, density_solid(n),zoom_factor(n)
             READ(ifuBodyIn,*)xcentConstr(n),ycentConstr(n),zcentConstr(n)

             !  alpha(n)     = 0.0_CGREAL
          
        END SELECT ! canonical_body_type
        ! initialize variables
               
        xcentinit(n) = xcent(n)
        ycentinit(n) = ycent(n)
        zcentinit(n) = zcent(n)
               
        vxcent(n)    = vxcentTrans(n)
        vycent(n)    = vycentTrans(n)
        vzcent(n)    = vzcentTrans(n)
               
        angvxinit(n) = angvx(n)
        angvyinit(n) = angvy(n)
        angvzinit(n) = angvz(n)

     ENDDO ! n

999 CONTINUE

   END SUBROUTINE read_inputs_body
!---------------------------------------------

!---------------------------------------------
   SUBROUTINE check_inputs_body()

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays

    IMPLICIT NONE
    INTEGER           :: n

    DO n=1,nbody
       ! write variables
       PRINT*,'  iBody = ', n
       PRINT*,'  canonical_body_type   = ',canonical_body_type(n)
       PRINT*,'  body_dim              = ',body_dim(n)
       PRINT*,'  unstruc_surface_type  = ',unstruc_surface_type(n)
       PRINT*,'  boundary_motion_type  = ',boundary_motion_type(n)
       PRINT*,'  wall_type             = ',wall_type(n)
       PRINT*,'  forced_motion_spec    = ',forced_motion_spec(n)


       PRINT*,'   nPtsBodyMarker = ',nPtsBodyMarker(n)
       PRINT*,'   Canonical Body Type = ', canonical_body_type(n)

       SELECT CASE(canonical_body_type(n))

       CASE(ELLIPTIC_CYLINDER:GENERAL_CYLINDER)
          write(*,'(a,3E15.7)')'   RadiusX, RadiusY           = ',radiusx(n),radiusy(n)
          write(*,'(a,3E15.7)')'   XCent, YCent               = ',xcent(n),ycent(n)

       CASE(ELLIPSOID) 
          write(*,'(a,3E15.7)')'   RadiusX, RadiusY, RadiusZ  = ',radiusx(n),radiusy(n),radiusz(n)
          write(*,'(a,3E15.7)')'   XCent, YCent, ZCent        = ',xcent(n),ycent(n),zcent(n)

       CASE(UNSTRUCTURED_SURFACE)
          write(*,'(a,I15)')    '   totNumTriElem =',totNumTriElem(n) 
          write(*,'(a,3E15.7)')'   XCent, YCent, ZCent        = ',xcent(n),ycent(n),zcent(n) 
       END SELECT ! canonical_body_type

       write(*,'(a,3E15.7)')'   Alpha                         = ',angzinit(n)
       write(*,'(a,3E15.7)')'   VXCentTrans, VYCentTrans, VZCentTrans     = ',&
            vxcentTrans(n),vycentTrans(n),vzcentTrans(n)
       write(*,'(a,3E15.7)')'   AmpX, AmpY, AmpZ             = ',ampx(n),ampy(n),ampz(n)
       write(*,'(a,3E15.7)')'   FreqX, FreqY, FreqZ          = ',freqx(n),freqy(n),freqz(n)
       write(*,'(a,3E15.7)')'   AngVx, AngVy, AngVz          = ',angvxinit(n),angvyinit(n),angvzinit(n)
       write(*,'(a,3E15.7)')'   AmpAngX, AmpAngY, AmpAngZ    = ',ampangx(n),ampangy(n),ampangz(n)
       write(*,'(a,3E15.7)')'   FreqAngX, FreqAngY, FreqAngZ = ',freqangx(n),freqangy(n),freqangz(n)
       write(*,'(a,3E15.7)')'   Density fluid, Density Solid = ',density_fluid, density_solid(n)         

    ENDDO ! n

   END SUBROUTINE check_inputs_body
!---------------------------------------------

!---------------------------------------------
   SUBROUTINE init_simulation()

    USE global_parameters
    USE flow_parameters
    USE fsi_module
    USE flow_arrays
    USE MPI_module

    use grid_arrays

    use boundary_arrays

    IMPLICIT NONE

    INTEGER              :: i, j,   k

    CALL allocate_memory()
  
    CALL BOUNDARY_allocate_memory()
    
    CALL make_grid()

    CALL metrics()

    if(it_solver_type == IT_SOLVER_TYPE_MG)  call mg_initial() 

    !Ye,debug
    if (lProc.eq.proc_m)then
       write(*,*)'Entering init_flow...'
    endif

    CALL init_flow()

    !Testing send and receive subroutine
    !call test_send_receive()
    !stop

    ! initialize markers and iblank
    !PRINT*,'CALLING initialize_marker from MAIN.F90'
    !print*, boundary_motion,moving_boundary

    !Ye,debug
    if (lProc.eq.proc_m)then
       write(*,*)'Entering initialize_marker...'
    endif
 
    call MPI_BARRIER(flow_comm,ierr)
   
    !read from 'unstruc_surf_in'
    CALL initialize_marker()    ! set initial location and velocity of markers

    call MPI_BARRIER(flow_comm,ierr)

    ! compute new iBlank, define ghostcells, 
    CALL set_iblank_body_fast()

            CALL write_dump()
            !call sleep(15)
            !STOP

    call MPI_BARRIER(flow_comm,ierr)

    !Ye,debug
    if (lProc.eq.proc_m)then
       write(*,*)'Exiting set_iblank_body_fast...'
    endif

    ! define ghostcells, compute projection, and prepare interpolation
    CALL GCM_set_internal_boundary()

    !Ye,debug
    if (lProc.eq.proc_m)then
       write(*,*)'Exiting GCM_set_internal_boundary...'
    endif

    ! Internal boundaries conditions at the body-intercept
    CALL GCM_SetBodyInterceptValues()
    
    !Ye,debug
    if (lProc.eq.proc_m)then
       write(*,*)'Exiting GCM_SetBodyInterceptValues...'
    endif
    !Ye,test
    call MPI_BARRIER(flow_comm,ierr)
    !STOP

    CALL face_vel(u,v,w)

    if(Imonitor_drag  == 1 .and. lProc == PROC_M) CALL open_drag_files()    
    if(Imonitor_probe == 1 .and. lProc == PROC_M) then
       CALL read_probe_inputs()
       CALL open_probe_files()
    endif

    if(elastic_present == 1) then
       !PRINT*,'ENTERING initialize_fsi()'    
       !CALL initialize_fsi() 
    else
       itermax_FSI = 1       ! only 1 iteration step for solving the NS equation.
    endif

    !Prepare iblank and ghostcell marks for multigrids
    !if(it_solver_type == IT_SOLVER_TYPE_MG)  call mg_prepare()
 
    rk_bet(1) = 8.0_CGREAL/15.0_CGREAL
    rk_bet(2) = 2.0_CGREAL/15.0_CGREAL
    rk_bet(3) = 1.0_CGREAL/ 3.0_CGREAL

    rk_gam(1) = 8.0_CGREAL/15.0_CGREAL
    rk_gam(2) = 5.0_CGREAL/12.0_CGREAL
    rk_gam(3) = 3.0_CGREAL/ 4.0_CGREAL

    rk_zet(1) =  0.0_CGREAL
    rk_zet(2) =-17.0_CGREAL/60.0_CGREAL
    rk_zet(3) = -5.0_CGREAL/12.0_CGREAL

    is_PBC_homogeneous = 0             

    !Ye,test constant fluid force
    !if (lProc.eq.proc_m)then
    !   open(unit=1,file='ventrical_num.dat')
    !   do i=1,10482
    !      read(1,*)j,v_num(i)
    !   enddo
    !   close(1)
    !   write(*,*)'finish reading v_num.dat!!',lProc
    !endif
    !Broadcast it to each processor
    !call MPI_BCAST(v_num, 10482, MPI_INTEGER,proc_m,flow_comm, ierr)
    !write(*,*)'finish broadcasting v_num.dat!!',lProc

   END SUBROUTINE init_simulation

!---------------------------------------------
   SUBROUTINE stop_program()

   USE MPI_module

    IMPLICIT NONE

    !CALL finalize_solid()
    CALL deallocate_memory()
   ! call MPI_Finalize(ierr)

   END SUBROUTINE stop_program
!-------------------------------------------  
   
!------------------------------------------- 
   SUBROUTINE init_flow()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE MPI_module
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER                   :: i
    INTEGER                   :: k,j    

    reinv   = 1.0_CGREAL/re
    dtinv   = 1.0_CGREAL/dt
    
    time    = 0.0_CGREAL
    time_start = 0.0_CGREAL

    ntime_start = 0
    ntime   = 0

    u(0:,yb1:,zb1:)       = uinit
    v(0:,yb1:,zb1:)       = vinit
    w(0:,yb1:,zb1:)       = winit
    p(0:,yb1:,zb1:)       = 0.0_CGREAL
    pPrime (0:,yb1:,zb1:) = 0.0_CGREAL
    p_temp (0:,yb1:,zb1:) = 0.0_CGREAL
    viscTot(0:,yb1:,zb1:) = reinv

    !Ye, pressure initialization using linear distribution===
    !do i=1,nx-1
    !   p(i:,0:,zb1:)=((pppx2-pppx1)*xc(i)+(pppx1*xout-pppx2*xOrigin))/(xout-xOrigin)
    !enddo
    !pPrime = p
    !===============================

    IF (nread == 1) then
       !PRINT*,'ENTERING READ_RESTART_FLOW'
       CALL read_restart()
       ntime = ntime_start
       time_start = time
    ELSE

       !Ye
       if (lProc.eq.proc_m)then
          write(*,*)'Entering set_outer_velocity_bc'
       endif

       ! set the outer boundary condition and the ghost velocity
       CALL set_outer_velocity_bc()
       CALL set_outer_ghost_velocity()

       !Ye
       if (lProc.eq.proc_m)then
          write(*,*)'Exiting set_outer_ghost_velocity'
       endif
 
    !Ye,test
    !write(600+lProc,*)ntime, lProc
    !write(600+lProc,*)yb1,yb2
    !write(600+lProc,*)zb1,zb2
    !write(600+lProc,*)'==='

       CALL set_outer_pressure_bc   (pPrime,nx,yb1,yb2,zb1,zb2,1)  ! inhomogeneous B.C.
       CALL set_outer_ghost_pressure(pPrime,nx,yb1,yb2,zb1,zb2,1)

    !Ye,DEBUG
    !write(700+lProc,*)'VARIABLES="X","Y","Z","P"'
    !write(700+lProc,*)'ZONE F=POINT, I=',nx+1,', J=',jSlices+2,' ,K=',kSlices+2
    !do k=zc_start-1,zc_end+1   !1,nz-1
    !do j=yc_start-1,yc_end+1
    !do i=0,nx
    !   write(700+lProc,'(4(2X,F12.5))')xc(i),yc(j),zc(k),pPrime(i,j,k)
    !enddo
    !enddo
    !enddo
    !close(700+lProc)
    !call mpi_barrier(MPI_COMM_WORLD,ierr)
    !stop

    ENDIF

   END SUBROUTINE init_flow
!-------------------------------------------     
   SUBROUTINE read_restart()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE mg_module
    use grid_arrays
    use multiuse_arrays
    use MPI_module

    implicit none
    integer::iBody, m, npts,nelm
    real*8:: dummy
    CHARACTER*20              :: fname1

!**********************************

    write(fname1,'(I4.4)') lProc
    fname1 = 'restart_'//trim(fname1)//'.dat'
    print*,'Read from file: ',fname1

    OPEN(UNIT=ifuRstrtFlowIn,FILE=fname1,STATUS='old',FORM='unformatted')
    read(ifuRstrtFlowIn)  time,ntime_start

    read(ifuRstrtFlowIn)                          &
               u(:,:,:), v(:,:,:), w(:,:,:), p(:,:,:),   &
               face_u(:,:,:),face_v(:,:,:),face_w(:,:,:),&
               pPrime(:,:,:), &   
               bcxu, bcxv, bcxw, bcxp,   &
               bcyu, bcyv, bcyw, bcyp,   &
               bczu, bczv, bczw, bczp,   &
               pgradx1, pgradx2, pgrady1,&
               pgrady2, pgradz1, pgradz2,&               
               iblank(:,:,:), ghostcellmark(:,:,:)
    
    !IF ( boundary_motion == MOVING_BOUNDARY ) THEN
    DO iBody = 1, nBody
       npts = nPtsBodyMarker(iBody)
       nelm = totNumTriElem(iBody)
       read(ifuRstrtFlowIn)                   &
                   xBodyMarker(1:npts,iBody), &
                   yBodyMarker(1:npts,iBody), &
                   zBodyMarker(1:npts,iBody), &
                   uBodyMarker(1:npts,iBody), &
                   vBodyMarker(1:npts,iBody), &
                   wBodyMarker(1:npts,iBody), &
                   xMarkerForce(1:npts,iBody), &
                   yMarkerForce(1:npts,iBody), &
                   zMarkerForce(1:npts,iBody), &
                   triElemNeig(1:3,1:nelm,iBody),               &
                   xcent(ibody), ycent(ibody), zcent(iBody),    &
                   angx (ibody), angy (ibody), angz (ibody),    &
                   vxcent(ibody), vycent(ibody), vzcent(ibody), &
                   angvx (ibody), angvy (ibody), angvz (ibody)
    ENDDO
    ! ENDIF

    xMarkerForceOld = xMarkerForce
    yMarkerForceOld = yMarkerForce
    zMarkerForceOld = zMarkerForce

   CLOSE(ifuRstrtFlowIn)

   END SUBROUTINE read_restart
!-------------------------------------------     
   SUBROUTINE write_restart()
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE mg_module
    use grid_arrays
    use multiuse_arrays
    use MPI_module

    implicit none
    integer       ::i,j,k,iBody,npts,nelm
    CHARACTER*20  :: fname1

!**********************************

    !A nice way of writing file name
    write(fname1,'(I4.4)') lProc
    fname1 = 'restart_'//trim(fname1)//'.dat'
    print*,'Write to file: ',fname1

    !IF ( lProc >= 0 .AND.  lProc .le. 9  )           &
    !    write(fname1,340) lProc 
    !IF ( lProc >= 10 .AND. lProc .le. 99 )           &
    !    write(fname1,341) lProc
    !IF ( lProc >= 100 .AND. lProc .le. 999 )         &
    !     write(fname1,342) lProc                    
    !IF ( lProc >= 1000 .AND. lProc .le. 9999 )       &
    !     write(fname1,343) lProc  


340 format('restart_000',i1,'.dat')
341 format('restart_00', i2,'.dat')
342 format('restart_0',  i3,'.dat')
343 format('restart_',   i4,'.dat')

    OPEN(UNIT=ifuRstrtFlowOut,FILE=fname1,FORM='unformatted')

    write(ifuRstrtFlowOut)  time,ntime
    write(ifuRstrtFlowOut)                   &
               u(:,:,:), v(:,:,:), w(:,:,:), p(:,:,:),   &
               face_u(:,:,:),face_v(:,:,:),face_w(:,:,:),&
               pPrime(:,:,:), &   
               bcxu, bcxv, bcxw, bcxp,   &
               bcyu, bcyv, bcyw, bcyp,   &
               bczu, bczv, bczw, bczp,   &
               pgradx1, pgradx2, pgrady1,&
               pgrady2, pgradz1, pgradz2,&               
               iblank(:,:,:), ghostcellmark(:,:,:)
    
    !IF ( boundary_motion == MOVING_BOUNDARY ) THEN
    DO iBody = 1, nBody
       npts = nPtsBodyMarker(iBody)
       nelm = totNumTriElem(iBody)
       write(ifuRstrtFlowOut)                 &
                   xBodyMarker(1:npts,iBody), &
                   yBodyMarker(1:npts,iBody), &
                   zBodyMarker(1:npts,iBody), &
                   uBodyMarker(1:npts,iBody), &
                   vBodyMarker(1:npts,iBody), &
                   wBodyMarker(1:npts,iBody), &
                   xMarkerForce(1:npts,iBody), &
                   yMarkerForce(1:npts,iBody), &
                   zMarkerForce(1:npts,iBody), &
                   triElemNeig(1:3,1:nelm,iBody),               &
                   xcent(ibody), ycent(ibody), zcent(iBody),    &
                   angx (ibody), angy (ibody), angz (ibody),    &
                   vxcent(ibody), vycent(ibody), vzcent(ibody), &
                   angvx (ibody), angvy (ibody), angvz (ibody)
   ENDDO
   ! ENDIF
   CLOSE(ifuRstrtFlowOut)

   END SUBROUTINE  write_restart
!-------------------------------------------
   SUBROUTINE spanwise_pert()
 
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE MPI_module

    IMPLICIT NONE
 
    INTEGER :: I, J, K
    REAL(KIND=CGREAL) :: har1, har2
 
    DO K = zb1,zb2   ! 0,nz
    DO J = 0, ny
    DO I = 0, nx
       CALL RANDOM_NUMBER(har1)
       CALL RANDOM_NUMBER(har2)
       u(i,j,k) = u(i,j,k) + vper*(har1-har2)
       v(i,j,k) = v(i,j,k) + vper*(har1-har2)
       w(i,j,k) = w(i,j,k) + vper*(har1-har2)
       !face_u(i,j,k) = face_u(i,j,k) + vper*(har1-har2)
       !face_v(i,j,k) = face_v(i,j,k) + vper*(har1-har2)
       !face_w(i,j,k) = face_w(i,j,k) + vper*(har1-har2)
    END DO
    END DO
    END DO
     
   END SUBROUTINE spanwise_pert
!-------------------------------------------

!-------------------------------------------     
   SUBROUTINE write_res()
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE mg_module
    use grid_arrays
    use multiuse_arrays
    use MPI_module

    implicit none
    integer       ::i,j,k,iBody,npts,nelm
    CHARACTER*20  :: fname1

!**********************************

    !A nice way of writing file name
    write(fname1,'(I4.4)') lProc
    fname1 = 'res_'//trim(fname1)//'.dat'
    print*,'Write to file: ',fname1

    !IF ( lProc >= 0 .AND.  lProc .le. 9  )           &
    !    write(fname1,340) lProc 
    !IF ( lProc >= 10 .AND. lProc .le. 99 )           &
    !    write(fname1,341) lProc
    !IF ( lProc >= 100 .AND. lProc .le. 999 )         &
    !     write(fname1,342) lProc                    
    !IF ( lProc >= 1000 .AND. lProc .le. 9999 )       &
    !     write(fname1,343) lProc  


340 format('res_000',i1,'.dat')
341 format('res_00', i2,'.dat')
342 format('res_0',  i3,'.dat')
343 format('res_',   i4,'.dat')

    OPEN(UNIT=ifuRstrtFlowOut,FILE=fname1,STATUS='UNKNOWN')

    write(ifuRstrtFlowOut,'(2(1X,1PE16.8))')  time,ntime

    do i=0,nx+1
    do j=0,ny+1
    do k=zb1,zb2
       write(ifuRstrtFlowOut,'(8(1X,1PE16.8),2(1X,I6))')u(i,j,k),v(i,j,k),w(i,j,k), &
            p(i,j,k),face_u(i,j,k),face_v(i,j,k),face_w(i,j,k),pPrime(i,j,k), &
            iblank(i,j,k),ghostcellmark(i,j,k)
    enddo
    enddo
    enddo

    write(ifuRstrtFlowOut,'(18(1X,1PE16.8))')bcxu, bcxv, bcxw, bcxp,   &
                                            bcyu, bcyv, bcyw, bcyp,   &
                                            bczu, bczv, bczw, bczp,   &
                                            pgradx1, pgradx2, pgrady1,&
                                            pgrady2, pgradz1, pgradz2
    
    !IF ( boundary_motion == MOVING_BOUNDARY ) THEN
    DO iBody = 1, nBody
       npts = nPtsBodyMarker(iBody)
       nelm = totNumTriElem(iBody)
       do i=1,npts
          write(ifuRstrtFlowOut,'(9(1X,1PE16.8))')xBodyMarker(i,iBody), &
                                  yBodyMarker(i,iBody), &
                                  zBodyMarker(i,iBody),uBodyMarker(i,iBody), &
                                  vBodyMarker(i,iBody),wBodyMarker(i,iBody), &
                                  xMarkerForce(i,iBody),yMarkerForce(i,iBody), &
                                  zMarkerForce(i,iBody)
       enddo

       do i=1,nelm
          write(ifuRstrtFlowOut,'(3(1X,I8))')triElemNeig(1,i,iBody), &
                                             triElemNeig(2,i,iBody), &
                                             triElemNeig(3,i,iBody)
       enddo

        write(ifuRstrtFlowOut,'(12(1X,1PE16.8))')xcent(ibody), ycent(ibody), &
                                                 zcent(iBody), &
                                                 angx(ibody), angy(ibody), &
                                                 angz(ibody), &
                                                 vxcent(ibody), vycent(ibody), &
                                                 vzcent(ibody), &
                                                 angvx(ibody), angvy(ibody), &
                                                 angvz(ibody)

   ENDDO
   ! ENDIF
   CLOSE(ifuRstrtFlowOut)

   END SUBROUTINE  write_res

!-------------------------------------------     
   SUBROUTINE read_res()
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE mg_module
    use grid_arrays
    use multiuse_arrays
    use MPI_module

    implicit none
    integer       ::i,j,k,iBody,npts,nelm
    CHARACTER*20  :: fname1

!**********************************

    !A nice way of writing file name
    write(fname1,'(I4.4)') lProc
    fname1 = 'res_'//trim(fname1)//'.dat'
    print*,'Read from file: ',fname1

    OPEN(UNIT=ifuRstrtFlowIn,FILE=fname1,STATUS='UNKNOWN')
    READ(ifuRstrtFlowIn,'(2(1X,1PE16.8))')  time,ntime_start

    do i=0,nx+1
    do j=0,ny+1
    do k=zb1,zb2
       READ(ifuRstrtFlowIn,'(8(1X,1PE16.8),2(1X,I6))')u(i,j,k),v(i,j,k),w(i,j,k), &
           p(i,j,k), &
           face_u(i,j,k),face_v(i,j,k),face_w(i,j,k),pPrime(i,j,k), &
           iblank(i,j,k),ghostcellmark(i,j,k)
    enddo
    enddo
    enddo

    READ(ifuRstrtFlowIn,'(18(1X,1PE16.8))')bcxu, bcxv, bcxw, bcxp,   &
                                            bcyu, bcyv, bcyw, bcyp,   &
                                            bczu, bczv, bczw, bczp,   &
                                            pgradx1, pgradx2, pgrady1,&
                                            pgrady2, pgradz1, pgradz2
    
    DO iBody = 1, nBody
       npts = nPtsBodyMarker(iBody)
       nelm = totNumTriElem(iBody)
       do i=1,npts
          READ(ifuRstrtFlowIn,'(9(1X,1PE16.8))')xBodyMarker(i,iBody), &
                                  yBodyMarker(i,iBody), &
                                  zBodyMarker(i,iBody),uBodyMarker(i,iBody), &
                                  vBodyMarker(i,iBody),wBodyMarker(i,iBody), &
                                  xMarkerForce(i,iBody),yMarkerForce(i,iBody), &
                                  zMarkerForce(i,iBody)
       enddo

       do i=1,nelm
          READ(ifuRstrtFlowIn,'(3(1X,I8))')triElemNeig(1,i,iBody), &
                                           triElemNeig(2,i,iBody), &
                                           triElemNeig(3,i,iBody)
       enddo

       READ(ifuRstrtFlowIn,'(12(1X,1PE16.8))')xcent(ibody), ycent(ibody), &
                                              zcent(iBody), &
                                              angx(ibody), angy(ibody), &
                                              angz(ibody), &
                                              vxcent(ibody), vycent(ibody), &
                                              vzcent(ibody), &
                                              angvx(ibody), angvy(ibody), &
                                              angvz(ibody)

   ENDDO
   ! ENDIF

    xMarkerForceOld = xMarkerForce
    yMarkerForceOld = yMarkerForce
    zMarkerForceOld = zMarkerForce

   CLOSE(ifuRstrtFlowIn)

   END SUBROUTINE  read_res

