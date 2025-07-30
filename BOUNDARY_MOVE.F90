! Move body markers and centroids and set new velocity of marker points and centroid 

  SUBROUTINE move_boundary()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE flow_arrays
    USE usr_module 
    USE fsi_module  
    USE mpi_module

    IMPLICIT NONE

    INTEGER           :: ie,ndata,i,n,ifortCent,iBody,iFort,m,j
    REAL(KIND=CGREAL) :: dx_temp, dy_temp, dz_temp, ds_max

    ! Update body velocity

    if(lProc .eq. PROC_M) then

       ! Update body velocity
       DO iBody=1,nBody
          SELECT CASE (boundary_motion_type(iBody))
          CASE (PRESCRIBED_RIGID)
             CALL forced_motion(iBody)
             CALL compute_marker_vel(iBody)
             CALL compute_marker_poisition(iBody)

          CASE (INDUCED_RIGID)
             CALL induced_rigid_motion(iBody)

          CASE (PRESCRIBED_ARBITRARY)
             CALL user_marker_vel(iBody)
             CALL compute_marker_poisition(iBody)

          END SELECT
       ENDDO

       if (elastic_present .eq. 1) then
          CALL get_elastic_marker_pos_and_vel() !get resMax_FSIdsp/v
       endif

       if (elastic_present .eq. 1) then
          is_fsi_cvg    = 0   ! set to zero first

          if ((iter_FSI > 1 .and. resMax_FSIdsp <= restol_FSIdsp .and.        &
               resMax_FSIv <= restol_FSIv .and. resMax_FSIp <= restol_FSIp ) &
               .or. iter_FSI >= itermax_FSI)                                  &
               is_fsi_cvg = 1

          !send the convergence info to the solid solver
          do ie=1,nElastic
             call MPI_SEND(is_fsi_cvg, 1, MPI_INTEGER,            &
                           proc_s(ie), 1, MPI_COMM_WORLD,istatus,ierr)
          enddo
       else
          is_fsi_cvg    = 1       ! for non-FSI, this is set to 1, so that the FSI
                                  ! loop is executed only once.
       endif  ! end elastic_present
    endif  ! proc_m

    ! broadcast the boundary location and velocity to all processors from PROC_M
    ndata = nPtsMax*nBody

    call MPI_BCAST(xBodyMarker, ndata, MPI_DOUBLE_PRECISION,proc_m,flow_comm, ierr)
    call MPI_BCAST(yBodyMarker, ndata, MPI_DOUBLE_PRECISION,proc_m,flow_comm, ierr)
    call MPI_BCAST(zBodyMarker, ndata, MPI_DOUBLE_PRECISION,proc_m,flow_comm, ierr)

    call MPI_BCAST(uBodyMarker, ndata, MPI_DOUBLE_PRECISION,proc_m,flow_comm, ierr)
    call MPI_BCAST(vBodyMarker, ndata, MPI_DOUBLE_PRECISION,proc_m,flow_comm, ierr)
    call MPI_BCAST(wBodyMarker, ndata, MPI_DOUBLE_PRECISION,proc_m,flow_comm, ierr)

    ! Broadcast the residule.  This is important.  Otherwise the processes won't
    ! terminate the FSI loop simultaneously!  --Haoxiang
    call MPI_BCAST(resMax_FSIdsp, 1, MPI_DOUBLE_PRECISION,proc_m,flow_comm, ierr)
    call MPI_BCAST(resMax_FSIv  , 1, MPI_DOUBLE_PRECISION,proc_m,flow_comm, ierr)
    call MPI_BCAST(resMax_FSIp  , 1, MPI_DOUBLE_PRECISION,proc_m,flow_comm, ierr)
    call MPI_BCAST(is_fsi_cvg   , 1, MPI_INTEGER ,        proc_m,flow_comm, ierr)

    !All processors calculate boundary variables
    call MPI_BARRIER(flow_comm,ierr)
    CALL calculate_arclength_norm_ds()

    if(ntime > 1) then
      ! dp/dn = -D u_n/Dt

      DO n=1,nBody
      DO i=1,nPtsBodyMarker(n)
        dpdnBodyMarker(i,n) = (uBodyMarker(i,n) - uBodyMarkerOld(i,n))/dt*xNormBodyMarker(i,n)  &
                            + (vBodyMarker(i,n) - vBodyMarkerOld(i,n))/dt*yNormBodyMarker(i,n)  &
                            + (wBodyMarker(i,n) - wBodyMarkerOld(i,n))/dt*zNormBodyMarker(i,n)

        dpdnBodyMarker(i,n) = -dpdnBodyMarker(i,n)
        !write(*,'(I6,3F12.5)')i, xNormBodyMarker(i,n), yNormBodyMarker(i,n), dpdnBodyMarker(i,n)
      ENDDO ! i
      ENDDO ! n
    else
      dpdnBodyMarker = 0.0_CGREAL
    endif

    if ( mod(ntime, nmonitor) == 0 .and. lProc==PROC_M) then
       do iBody=1,nBody
          write(*,*) 'Body #', iBody
          write(*,'(a,2E15.5)') ' Min-Max of uBodyMarker:',minval(uBodyMarker(:,iBody)),maxval(uBodyMarker(:,iBody))
          write(*,'(a,2E15.5)') ' Min-Max of vBodyMarker:',minval(vBodyMarker(:,iBody)),maxval(vBodyMarker(:,iBody))
          write(*,'(a,2E15.5)') ' Min-Max of wBodyMarker:',minval(wBodyMarker(:,iBody)),maxval(wBodyMarker(:,iBody))
       enddo
    endif

 !   if(elastic_present==1 .and. iter_FSI>0 .and. is_fsi_cvg==0 ) then  !reset the velocity 
 !      u = uOld
 !      v = vOld
 !      w = wOld
 !   endif

    ! Calculate the max. displacement of the boundary w.r.t. last snapshot;
    ! determine if iblank needs to be re-done
    ds_max = 0.0d0
    DO n=1,nBody
    DO i=1,nPtsBodyMarker(n)
       dx_temp = xBodyMarker(i,n) - xBodyMarkerOld(i,n)
       dy_temp = yBodyMarker(i,n) - yBodyMarkerOld(i,n)
       dz_temp = zBodyMarker(i,n) - zBodyMarkerOld(i,n)

       dsBodyMarker(i,n) = sqrt(dx_temp**2 + dy_temp**2 + dz_temp**2)
       if(dsBodyMarker(i,n) > ds_max) ds_max = dsBodyMarker(i,n)
    ENDDO ! i
    ENDDO ! n
    
    IF(ds_max > 0.05d0*min(dxMin,dyMin,dzMin)) THEN
       iblank_reset = 1

       xBodyMarkerOld = xBodyMarker
       yBodyMarkerOld = yBodyMarker
       zBodyMarkerOld = zBodyMarker
       if(lProc==PROC_M) write(*,'(a,E15.5,a,E15.5)')'ds_max = ', ds_max,     &
                         ' ; ds_max/min(dxyz)=',ds_max/min(dxMin,dyMin,dzMin)
    ELSE
       iblank_reset = 0
    ENDIF

  END SUBROUTINE move_boundary
!
!---------------------------------------------------------------------------
!---------------------------------------------------------------                                                      
! integrate the marker point velocity to obtain the                                                                   
! marker point position                                                                                               

  SUBROUTINE compute_marker_poisition(iBody)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

!                                                                                                                     
!      n+1     n                                                                                                      
!    x    -  x                                                                                                        
!    - i     - i        n+1                                                                                           
!   ------------   =  v                                                                                               
!         dt          - i                                                                                             

!...Update location of boundary points                                                                                
!...Update centroid location                                                                                          

    INTEGER,INTENT (IN) :: iBody
    INTEGER             :: i,n

    n = iBody

    DO i=1,nPtsBodyMarker(n)
       xBodyMarker(i,n) = xBodyMarker(i,n) + dt*uBodyMarker(i,n)
       yBodyMarker(i,n) = yBodyMarker(i,n) + dt*vBodyMarker(i,n)
       zBodyMarker(i,n) = zBodyMarker(i,n) + dt*wBodyMarker(i,n)
    ENDDO ! i                                                                                                         

    xcent(n) = xcent(n) + dt*vxcent(n)
    ycent(n) = ycent(n) + dt*vycent(n)
    zcent(n) = zcent(n) + dt*vzcent(n)

  END SUBROUTINE compute_marker_poisition
!---------------------------------------------------------------------------                                          

  SUBROUTINE save_lastStep()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE flow_arrays
    Use usr_module
    USE fsi_module
    use MPI_module

    IMPLICIT NONE

    INTEGER           :: i

    angvx_old(1:nbody) = angvx(1:nbody)
    angvy_old(1:nbody) = angvy(1:nbody)
    angvz_old(1:nbody) = angvz(1:nbody)

    uBodyMarkerOld = uBodyMarker  !
    vBodyMarkerOld = vBodyMarker  !
    wBodyMarkerOld = wBodyMarker  !

    uOld = u
    vOld = v
    wOld = w

    !Ye,use constant force at specific timestep for test
    !v_num=1:constant force
    !if (time.ge.10000.0d0)then
    !   do i=1,10590
    !      if (v_num(i).eq.0)then
    !         xMarkerForceOld(i,2) = xMarkerForce(i,2)
    !         yMarkerForceOld(i,2) = yMarkerForce(i,2)
    !         zMarkerForceOld(i,2) = zMarkerForce(i,2)
    !      else
             !Const fluid force
    !         xMarkerForce(i,2) = xMarkerForceOld(i,2)
    !         yMarkerForce(i,2) = yMarkerForceOld(i,2)
    !         zMarkerForce(i,2) = zMarkerForceOld(i,2)
    !      endif
    !   enddo
    !else
       xMarkerForceOld = xMarkerForce
       yMarkerForceOld = yMarkerForce
       zMarkerForceOld = zMarkerForce
    !endif

333 continue

  END SUBROUTINE save_lastStep
!---------------------------------------------------------------------------                                          
