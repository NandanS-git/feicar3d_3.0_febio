!---------------------------------------------
! SUBROUTINE initialize_fsi()
! SUBROUTINE deallocate_fsi()
!
!---------------------
!
! Author: Haoxiang Luo
! September, 2009
!---------------------------------------------
!-------------------------------------------------------------------------------------
   SUBROUTINE initialize_fsi()

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE fsi_module
    USE MPI_module

!*******************************************************************
        
    IMPLICIT NONE

    INTEGER :: i,j,k,n,m, idfile, iBody,ie
    REAL(KIND=CGREAL)  :: dtemp, xShift, yShift, zShift
    REAL(KIND=CGREAL)  :: dmax,temp, temp_id

    allocate(bodyInFlowID(nBody))

    if (elastic_present .eq. 0) then
       print*, 'No elastic body is present.'
       goto 1000
    else
        write(*,*) "lProc_g:",lProc_g," Elastic body is present: nElastic=",nElastic
    endif

    nElasticBody    = 0

    DO n = 1,nBody
       IF (boundary_motion_type(n) == ELASTIC_MOTION ) THEN
          !elastic_present = 1
          
          nElasticBody = nElasticBody + 1
          bodyInFlowID(nElasticBody) = n
       ENDIF
    ENDDO

    if(nElasticbody /= nElastic) then
       print*,' Inconsistent nElastic: '
       print*, 'nElasticBody = ',  nElasticBody, nElastic
    endif

    allocate(elasticMarkerID(nPtsMax,1))
    !allocate(xf_elem        (nPtsMax,1))
    !allocate(yf_elem        (nPtsMax,1))
    !allocate(zf_elem        (nPtsMax,1))

    ifuElasId = 200
    !open(unit=ifuElasId, file='elastic_nodeID.dat')
    !rewind ifuElasId

!    DO n = 1, nElastic  
!       iBody = bodyInFlowID(n)      ! get the corresponding body id in flow
!       DO m=1, nPtsBodyMarker(iBody)
!          elasticMarkerID(m, n) = m  !set the corresponding node id in solid
!       ENDDO
!    ENDDO

   do ie = 1, nElastic
    !check if the parameters from two solvers are the same.
    call MPI_RECV(dt_s, 1, MPI_DOUBLE_PRECISION, proc_s(ie), 1,  &
                  MPI_COMM_WORLD,istatus,ierr)

      IF( lProc.eq.PROC_M) THEN
      print*, "dt, dt_s", dt, dt_s
      ENDIF

 !    nsub = NINT(dt/dt_s)
     if(dt .ne. dt_s) then
 !      if(abs(dt/dt_s-NINT(dt/dt_s)).lt.1e-7) then
 !        print*, 'Nonstad substep: ', nsub
 !      else
         print*, "Error: step sizes in the two solvers do not match!"
!		 print*, "Step sizes in fluid is forced to be the same with that in solid solver"
!		 dt=dt_s 
         stop
!       endif
     endif

    call MPI_RECV(nread_s, 1, MPI_INTEGER, proc_s(ie), 1,  &
                  MPI_COMM_WORLD,istatus,ierr)
    if(nread .ne. nread_s) then
       print*, "ERROR: nread in the two solvers do not match!"
       stop
    endif

    call MPI_RECV(nrestart_s, 1, MPI_INTEGER, proc_s(ie), 1,  &
                  MPI_COMM_WORLD,istatus,ierr)
!    if(nrestart .ne. nrestart_s/nsub) then
     if(nrestart .ne. nrestart_s) then
       print*, "ERROR: nrestart in the two solvers do not match!"
       stop
    endif

    call MPI_RECV(nstep_s, 1, MPI_INTEGER, proc_s(ie), 1,  &
                  MPI_COMM_WORLD,istatus,ierr)
!    if(no_tsteps .ne. (nstep_s-1)/nsub) then
     if(no_tsteps .ne. (nstep_s-1)) then
       print*, no_tsteps, nstep_s
       print*, "Stopped. Set structure steps to be ", 1+no_tsteps
       stop
    endif


    ! receive the initial marker point position and check 
    ! if matches the flow solver
    call MPI_RECV(nPtsMax_s, 1, MPI_INTEGER, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

    print*, 'nPtsMax_s = ', nPtsMax_s
    if(ie.eq.1) then
    allocate(xdata_s(nPtsMax_s, 1))
    allocate(ydata_s(nPtsMax_s, 1))
    allocate(zdata_s(nPtsMax_s, 1))
    endif
    xdata_s(:, :) = 0.0_CGREAL
    ydata_s(:, :) = 0.0_CGREAL
    zdata_s(:, :) = 0.0_CGREAL

    ndata_s = nPtsMax_s * 1

    call MPI_RECV(xdata_s, ndata_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1,  &
                  MPI_COMM_WORLD,istatus,ierr)


    call MPI_RECV(ydata_s, ndata_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1,  &
                  MPI_COMM_WORLD,istatus,ierr)

    call MPI_RECV(zdata_s, ndata_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1,  &
                  MPI_COMM_WORLD,istatus,ierr)

    dmax = 0.0d0
    iBody = bodyInFlowID(ie)      ! get the corresponding body id in flow
   
    !print*, ndata_s, ie, iBody, nPtsBodyMarker(iBody)

    DO m=1, nPtsBodyMarker(iBody)
       dtemp = abs(xBodyMarker(m,iBody) - xdata_s(m,1)) &
             + abs(yBodyMarker(m,iBody) - ydata_s(m,1)) &
             + abs(zBodyMarker(m,iBody) - zdata_s(m,1))
       !write(333,'(7F12.5)') xBodyMarker(m,iBody),yBodyMarker(m,iBody),zBodyMarker(m,iBody) &
       !                     ,xdata_s(m,1),ydata_s(m,1),zdata_s(m,1),dtemp
       if(dtemp .gt. dmax) dmax=dtemp
    ENDDO ! m
 
    print*,'Marker difference from two solver: ', dmax
    if (dmax .ge. restol_FSIdsp) then  !restol_FSIdsp=0.0001 from input.dat
       print*, 'Initialize_fsi: markers in two solvers do not match!'
       print*, ibody, dmax
       print*, xBodyMarker(1:3,iBody)
       print*, xdata_s(1:3,1)
       stop
       ! write(*,'(4I6, 1PE12.5)') iBody, m, n, i, dtemp
       ! write(*,'(4F15.7)')xBodyMarker(m,iBody), xdata_s(i,n),  &
       !                    yBodyMarker(m,iBody), ydata_s(i,n),  &
       !                    zBodyMarker(m,iBody), zdata_s(i,n)
    endif
                
    !call MPI_send(dt, 1, MPI_double_precision, proc_s(ie), 1, &
    !              MPI_COMM_WORLD,istatus,ierr)


    ENDDO    

 1000 CONTINUE

   END SUBROUTINE initialize_fsi
!-------------------------------------------------------------------------------------
   SUBROUTINE interp_force_to_elastic()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE fsi_module
    USE MPI_module

    IMPLICIT NONE

    INTEGER           :: i,j,k,m,n,iBody,kk,ie
    INTEGER           :: node(3)
    INTEGER           :: itime_s

    REAL(KIND=CGREAL) :: fac, xfmax, yfmax, zfmax
    INTEGER   :: maxforcexloc, maxforceyloc, maxforcezloc
    REAL(KIND=CGREAL) :: xMF, yMF, zMF

    ! receive the time level from the solid solver and check
    ! if matches the flow solver
    do ie= 1, nElastic

    call MPI_RECV(itime_s, 1, MPI_INTEGER, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

!    if(itime_s .ne. ntime-ntime_start) then
     if(itime_s .ne. ntime) then
       print*, proc_s(ie), 'WARNING: itime_s not equal to ntime.  Two solvers are at   &
       different time levels!'
       print*,'itime s = ', itime_s+1, 'ntime =', ntime
       stop
    endif

    call MPI_SEND(relax_FSIp, 1, MPI_DOUBLE_PRECISION, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

! compute the nodal force by multiplying the distributed force by 
! the effective area
!--------

    xdata_s(:, :) = 0.0_CGREAL
    ydata_s(:, :) = 0.0_CGREAL
    zdata_s(:, :) = 0.0_CGREAL

    fsi_on = 1

    fac = density_fluid
    !if( (time-time_start) < FSI_delay) then
    if ( (nread.eq.0).and.(time-time_start) < FSI_delay )then
       fsi_on = 0
       fac = 1.0d0 - exp(-5.0d0*(time-time_start)/FSI_delay)
       fac =  density_fluid * fac
       print*,'FSI disabled because time < fsi_delay.' 
    endif
 
    iBody = bodyInFlowID(ie)      ! get the corresponding body id in flow

    ! integrate the filtered force over all the elements
    DO m = 1, nPtsBodyMarker(iBody) 
       xdata_s(m, 1) = xMarkerForce(m,iBody) * fac
       ydata_s(m, 1) = yMarkerForce(m,iBody) * fac
       zdata_s(m, 1) = zMarkerForce(m,iBody) * fac
    ENDDO     ! end m 

    xfmax = maxval(abs(xdata_s(1:nPtsMax_s, 1)))
    yfmax = maxval(abs(ydata_s(1:nPtsMax_s, 1)))
    zfmax = maxval(abs(zdata_s(1:nPtsMax_s, 1)))
    
    maxforcexloc = 1
    do m = 1, nPtsBodyMarker(iBody)
       if (abs(xdata_s(m,1)) >= abs(xdata_s(maxforcexloc,1))) then
          maxforcexloc = m
       endif
    enddo

    maxforceyloc = 1
    do m = 1, nPtsBodyMarker(iBody)
       if (abs(ydata_s(m,1)) >= abs(ydata_s(maxforceyloc,1))) then
          maxforceyloc = m
       endif
    enddo

    maxforcezloc = 1
    do m = 1, nPtsBodyMarker(iBody)
       if (abs(zdata_s(m,1)) >= abs(zdata_s(maxforcezloc,1))) then
          maxforcezloc = m
       endif
    enddo

    ! send the marker point force to the solid solver
    print*, 'Flow sending force to elastic body #', ie
    write(*,'(a, 3E15.5)') ' Maximum force = :', xfmax, yfmax, zfmax
    write(*,*)'Maxumum force @ :', maxforcexloc, maxforceyloc,  &
                                   maxforcezloc

    !call MPI_SEND(fsi_on, 1, MPI_INTEGER, proc_s(ie), 1, &
    !              MPI_COMM_WORLD,istatus,ierr)

    !call MPI_SEND(iter_FSI, 1, MPI_INTEGER,            &
    !                       proc_s(ie), 1, MPI_COMM_WORLD,istatus,ierr)

    !call MPI_SEND(density_fluid, 1, MPI_DOUBLE_PRECISION, proc_s(ie), 1, &
    !              MPI_COMM_WORLD,istatus,ierr)

    call MPI_SEND(xdata_s, ndata_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1,  &
                  MPI_COMM_WORLD,istatus,ierr)

    call MPI_SEND(ydata_s, ndata_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1,  &
                  MPI_COMM_WORLD,istatus,ierr)

    call MPI_SEND(zdata_s, ndata_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1,  &
                  MPI_COMM_WORLD,istatus,ierr)

    !write(*,100) 'Interp_force:', xdata_s(1:10,1)
    !write(*,100) 'Interp_force:', ydata_s(1:10,1)
    !write(*,100) 'Interp_force:', zdata_s(1:10,1)
    call MPI_recv(resMax_FSIp, 1, MPI_double_precision,            &
                  proc_s(ie), 1, MPI_COMM_WORLD,istatus,ierr)

    enddo !enddo ie

 100  format(a20, 10f10.5)

   END SUBROUTINE interp_force_to_elastic
!-------------------------------------------------------------------------------------
   SUBROUTINE get_elastic_marker_pos_and_vel()

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE fsi_module
    USE MPI_module

    IMPLICIT NONE

    INTEGER            :: i,j,k,m,n,iBody,ie
    REAL(KIND=CGREAL)  :: disp_max, dvel_max
    REAL(KIND=CGREAL)  :: resDispX, resDispY, resDispZ,  &
                          resVelX,  resVelY,  resVelZ
    INTEGER            :: max_disp_loc, max_vel_loc
    INTEGER            :: res_v_flag

    REAL(KIND=CGREAL)  :: beta
!****************************************************

    !Ye,velocity filter coeff
    !if (ntime.le.3000)then!.and.(ntime.le.1800))then
       beta = fsi_filter_v!From input.dat
    !elseif ((ntime.gt.3000).and.(ntime.le.4000))then
    !   beta = 0.225 * time -2.6
    !elseif (ntime.gt.4000)then
    !   beta = 1.0d0
    !endif

    !if (ntime.gt.3545)then
    !   beta = 0.01d0
    !endif

    resMax_FSIdsp = 0.0_CGREAL  ! maximum residual
    resMax_FSIv   = 0.0_CGREAL

  !Ye, adaptive relaxation factor for velocity
  if (iter_FSI.eq.0) then
     resMax_FSIv_old = 0.0_CGREAL
  endif

  do ie=1, nelastic
    ! Get the node position first
    print*, 'Flow receiving from elastic body #', ie

    xdata_s(:, :) = 0.0_CGREAL
    ydata_s(:, :) = 0.0_CGREAL
    zdata_s(:, :) = 0.0_CGREAL
    call MPI_RECV(xdata_s, ndata_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1,  &
                  MPI_COMM_WORLD,istatus,ierr)

    call MPI_RECV(ydata_s, ndata_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1,  &
                  MPI_COMM_WORLD,istatus,ierr)

    call MPI_RECV(zdata_s, ndata_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1,  &
                  MPI_COMM_WORLD,istatus,ierr)

    ! Update the node position with a relaxation factor
    iBody = bodyInFlowID(ie)      ! get the corresponding body id in flow

    DO m=1, nPtsBodyMarker(iBody)
       i = m 

       !Ye,debug
       !if( (m.eq.4180).and.(iBody.eq.2) )then
       !   write(55,'(8(E16.8,1X))')time,xdata_s(i,1),ydata_s(i,1),zdata_s(i,1),  &
       !         xBodyMarker(m,iBody),yBodyMarker(m,iBody),zBodyMarker(m,iBody)
       !endif

       resDispX  = xdata_s(i,1) - xBodyMarker(m,iBody)
       resDispY  = ydata_s(i,1) - yBodyMarker(m,iBody)
       resDispZ  = zdata_s(i,1) - zBodyMarker(m,iBody)

       xBodyMarker(m,iBody) = xBodyMarker(m,iBody) + relax_FSIdsp*resDispX
       yBodyMarker(m,iBody) = yBodyMarker(m,iBody) + relax_FSIdsp*resDispY
       zBodyMarker(m,iBody) = zBodyMarker(m,iBody) + relax_FSIdsp*resDispZ

       !Ye,add parts for the location of resMax
       IF (abs(resDispX).ge.resMax_FSIdsp) THEN
          resMax_FSIdsp= abs(resDispX)
          max_disp_loc=m
       ENDIF
       IF (abs(resDispY).ge.resMax_FSIdsp) THEN
          resMax_FSIdsp= abs(resDispY)
          max_disp_loc=m
       ENDIF
       IF (abs(resDispZ).ge.resMax_FSIdsp) THEN
          resMax_FSIdsp= abs(resDispZ)
          max_disp_loc=m
       ENDIF
       
       !resMax_FSIdsp= max(resMax_FSIdsp, abs(resDispX))
       !resMax_FSIdsp= max(resMax_FSIdsp, abs(resDispY))
       !resMax_FSIdsp= max(resMax_FSIdsp, abs(resDispZ))
    ENDDO ! m

    !write(*,100) 'Get_pos: ', xBodyMarker(1:10,1)
    !write(*,100) 'Get_pos: ', yBodyMarker(1:10,1)
    !write(*,100) 'Get_pos: ', zBodyMarker(1:10,1)

    ! Get the node velocity
    call MPI_RECV(xdata_s, ndata_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1,  &
                  MPI_COMM_WORLD,istatus,ierr)

    call MPI_RECV(ydata_s, ndata_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1,  &
                  MPI_COMM_WORLD,istatus,ierr)

    call MPI_RECV(zdata_s, ndata_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1,  &
                  MPI_COMM_WORLD,istatus,ierr)

    if ( mod(ntime, nmonitor) == 0 ) then
       write(*,'(a,2E15.5)') ' Min-Max of xMarkerVel:',minval(xdata_s(:,1)),maxval(xdata_s(:,1))
       write(*,'(a,2E15.5)') ' Min-Max of yMarkerVel:',minval(ydata_s(:,1)),maxval(ydata_s(:,1))
       write(*,'(a,2E15.5)') ' Min-Max of zMarkerVel:',minval(zdata_s(:,1)),maxval(zdata_s(:,1))
    endif

    ! Update the node velocity with a relaxation factor
    ! iBody = bodyInFlowID(ie)      ! get the corresponding body id in flow

!xdata_s=0.0d0
!ydata_s=0.0d0
!zdata_s=0.0d0

    !Ye, store true velocity from solid before filtered
    DO m=1, nPtsBodyMarker(iBody)
       i = m
       u1BodyMarker(m,iBody) = xdata_s(i,1)
       v1BodyMarker(m,iBody) = ydata_s(i,1)
       w1BodyMarker(m,iBody) = zdata_s(i,1)
    ENDDO

    !Ye, velocity filter
    DO m=1, nPtsBodyMarker(iBody)
       i = m
       xdata_s(i,1) = (1.0d0-beta)*uBodyMarkerOld(m,iBody) + beta*xdata_s(i,1)
       ydata_s(i,1) = (1.0d0-beta)*vBodyMarkerOld(m,iBody) + beta*ydata_s(i,1)
       zdata_s(i,1) = (1.0d0-beta)*wBodyMarkerOld(m,iBody) + beta*zdata_s(i,1)
    ENDDO

    res_v_flag = 0
    DO m=1, nPtsBodyMarker(iBody)
       i = m

       !Ye,debug
       !if( (m.eq.4180).and.(iBody.eq.2) )then
       !   write(77,'(8(E16.8,1X))')time,xdata_s(i,1),ydata_s(i,1),zdata_s(i,1),  &
       !         uBodyMarker(m,iBody),vBodyMarker(m,iBody),wBodyMarker(m,iBody)
       !endif

       resVelX  = xdata_s(i,1) - uBodyMarker(m,iBody)
       resVelY  = ydata_s(i,1) - vBodyMarker(m,iBody)
       resVelZ  = zdata_s(i,1) - wBodyMarker(m,iBody)
       
       uBodyMarker(m,iBody) = uBodyMarker(m,iBody) + relax_FSIv*resVelX
       vBodyMarker(m,iBody) = vBodyMarker(m,iBody) + relax_FSIv*resVelY
       wBodyMarker(m,iBody) = wBodyMarker(m,iBody) + relax_FSIv*resVelZ

       !Ye,add parts for the location of resMax
       IF (abs(resVelX).ge.resMax_FSIv) THEN
          resMax_FSIv= abs(resVelX)
          max_vel_loc=m
          res_v_flag = 1
       ENDIF
       IF (abs(resVelY).ge.resMax_FSIv) THEN
          resMax_FSIv= abs(resVelY)
          max_vel_loc=m
          res_v_flag = 2
       ENDIF
       IF (abs(resVelZ).ge.resMax_FSIv) THEN
          resMax_FSIv= abs(resVelZ)
          max_vel_loc=m
          res_v_flag = 3
       ENDIF

       !resMax_FSIv  = max(resMax_FSIv,   abs(resVelX ))
       !resMax_FSIv  = max(resMax_FSIv,   abs(resVelY ))
       !resMax_FSIv  = max(resMax_FSIv,   abs(resVelZ ))

    ENDDO ! m

   enddo ! ie,body#

!   if (iter_FSI.eq.1) then
      relax_FSIv = relax_FSIv0
      relax_FSIp = relax_FSIp0
!   endif

!   if( (iter_FSI.ge.2).and.(resMax_FSIv.gt.resMax_FSIv_old) ) then
!     relax_FSIv = relax_FSIv / 2.0d0
!     relax_FSIp = relax_FSIp / 1.3d0
!   endif

!   if (relax_FSIv.lt.0.50d0*relax_FSIv0) then
!      relax_FSIv = relax_FSIv0 / 2.0d0
!   endif

   resMax_FSIv_old = resMax_FSIv

    !write(*,100) 'Get_vel: ', uBodyMarker(1:10,1)
    !write(*,100) 'Get_vel: ', vBodyMarker(1:10,1)
    !write(*,100) 'Get_vel: ', wBodyMarker(1:10,1)

100 format(a20, 10f10.5)

   if ( mod(ntime, nmonitor) == 0 ) then

      !DO m = 1,nPtsBodyMarker(iBody)
      !   write(*,'(4F12.5)')xBodyMarker(m,iBody), yBodyMarker(m,iBody)
      !ENDDO

      ! we can compute the maxmum elastic deformation here
      
      !write(*,*) ' Maximum displacement, velocity, and residual of Body #', iBody
      !write(*,'(a,2E15.5)') ' Min-Max of uBodyMarker:',minval(uBodyMarker(:,iBody)),maxval(uBodyMarker(:,iBody))
      !write(*,'(a,2E15.5)') ' Min-Max of vBodyMarker:',minval(vBodyMarker(:,iBody)),maxval(vBodyMarker(:,iBody))
      !write(*,'(a,2E15.5)') ' Min-Max of wBodyMarker:',minval(wBodyMarker(:,iBody)),maxval(wBodyMarker(:,iBody))
      !write(*,'(a,2F15.5)') 'Min-Max of uBodyMOld:  ',minval(uBodyMarkerOld(:,iBody)),maxval(uBodyMarkerOld(:,iBody))
      !write(*,'(a,2F15.5)') 'Min-Max of vBodyMOld:  ',minval(vBodyMarkerOld(:,iBody)),maxval(vBodyMarkerOld(:,iBody))
      !write(*,'(a,2F15.5)') 'Min-Max of wBodyMOld:  ',minval(wBodyMarkerOld(:,iBody)),maxval(wBodyMarkerOld(:,iBody))
      write(*,'(a, 2I6,5(E12.4,1X), 3I6)')'FSI residual:',ntime,iter_FSI,resMax_FSIdsp,resMax_FSIv,resMax_FSIp &
                                                         ,relax_FSIv, relax_FSIp, res_v_flag
      write(*,'(a,4(E12.4,1X))')'force filter choose:', fsi_filter_f0, fsi_filter_f, resMax_FSIp_old, resMax_FSIp
      write(*,'(a, 3(E12.4,1X))')'U_BM:',uBodyMarker(max_vel_loc,2),vBodyMarker(max_vel_loc,2),wBodyMarker(max_vel_loc,2)
      write(*,*)'max FSI residual at:',max_disp_loc,max_vel_loc
      print*,'---------F-------S-------I-------- '

   endif
   
!   ! debug data transfer
!   iBody = 1
!   do m=1, nPtsBodyMarker(iBody)
!      write(412,'(I6, 6f10.5)')m,xBodyMarker(m,iBody),yBodyMarker(m,iBody),zBodyMarker(m,iBody)  &
!              ,uBodyMarker(m,iBody),vBodyMarker(m,iBody),wBodyMarker(m,iBody)
!   enddo
!   close(412)

   END SUBROUTINE get_elastic_marker_pos_and_vel

!-------------------------------------------------------------------------------------
   SUBROUTINE deallocate_fsi()
     USE flow_parameters
     USE boundary_arrays
     USE fsi_module

     IMPLICIT NONE

     deallocate(bodyInFlowID)
     deallocate(elasticMarkerID)
     !deallocate(xf_elem, yf_elem, zf_elem)
     deallocate(xdata_s, ydata_s, zdata_s)

   END SUBROUTINE deallocate_fsi
