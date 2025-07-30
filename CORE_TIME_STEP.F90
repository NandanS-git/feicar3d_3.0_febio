!-----------------------------------------------------------------------------
   SUBROUTINE time_step_viscous()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE multiuse_arrays
    USE stat_arrays
    USE stat_vort_arrays
    USE fsi_module
    USE MPI_module
    USE GCM_arrays
    USE solver_arrays

    IMPLICIT NONE

!... Local variables 

    REAL(KIND=CGREAL)    :: sumd
    INTEGER              :: clock0, clock1, clock2, clock3, clock4, clock5, clock_rate
    INTEGER              :: ie, i, m, iBody, ndata

    INTEGER              :: j,k
!******************************************************************************

    ! Start time stepping loop
    DO ntime = ntime_start+1,ntime_start+no_tsteps        ! move solution from n --> n+1
        
       call MPI_BARRIER(flow_comm,ierr)
       IF(lProc.eq.PROC_M) call system_clock(clock0)  

       call calc_flow_rate()
       call calc_goa()

       time = time + dt
       IF( lProc.eq.PROC_M) THEN
          WRITE(*,*)'---------------------------------------------------------'
          WRITE(*,*)'===============TIME STEP', ntime, time, '================'
          WRITE(*,*)'---------------------------------------------------------'
       END IF


       IF( MOD(ntime,nmonitor) == 0 .and. lProc.eq.PROC_M) THEN
          WRITE(*,*)'============================================================'
          WRITE(*,'(A,I6,A,F15.7,A,I6)') 'NTIME = ',ntime,',  Time = ',time,  &
               ',  NTIME TO GO = ',no_tsteps+ntime_start-ntime
          WRITE(4,*)'======================================================='
          WRITE(4,'(A,I6,A,F15.7,A,I6)') 'NTIME = ',ntime,',  Time = ',time,  &
                ',  NTIME TO GO = ',no_tsteps+ntime_start-ntime
          !Ye,write out the flow rate
          WRITE(*,*)'Q-1:',vq1,'Q-2:',vq2
          WRITE(*,*)'Q-3:',vq3,'Q-4:',vq4
          WRITE(*,*)'GOA:',goa_area
          call write_flowrate()
       ENDIF

       iter_FSI      = 0
       is_fsi_cvg    = 0
       resMax_FSIdsp = 1.0E10_CGREAL

       CALL save_lastStep() !Save all old velocity

       !----------
       ! FSI loop
       !----------
       DO WHILE(is_fsi_cvg .eq. 0 .and. iter_FSI < itermax_FSI)

          IF( lProc.eq.PROC_M) THEN
          WRITE(*,*)'===============FSI START====',iter_FSI,'====================='
          END IF

          rk_step = 0   ! necessary to set this to zero even for non-RK schemes

          IF (elastic_present .eqv. .true. ) THEN
             !Send force to solid nodes
             if(lProc .eq. PROC_M)  CALL interp_force_to_elastic() 
     
             !Ye, for ALL processors
             !reduce the force filter factor for convergence enhancement
             !from the 3rd iteration (before iter_fsi=iter_fsi+1, 3)
             !laterest resMax_FSIp from 'interp_force_to_elastic'
             !fsi_filter_f will be used in 'compute_marker_stress'
             !USE itermax_FSI TO DISABLE THIS IF BLOCK
             if ( (iter_FSI.GE.itermax_FSI).and.(resMax_FSIp.GE.resMax_FSIp_old) )then
                fsi_filter_f = fsi_filter_f0!/10.0
             else
                fsi_filter_f = fsi_filter_f0
             endif
          ENDIF

          call system_clock(clock3) !start the clock for the current FSI iteration
!------------------------------------------------------------------------------
!     Move Body Boundary and compute coefficients
!------------------------------------------------------------------------------

          IF ( ntime .eq. ntime_start+1) THEN
             CALL set_solve_ad()
          ENDIF

          IF ( boundary_motion == MOVING_BOUNDARY .AND. ntime >= 1) THEN

             if(iblank_reset == 1) then 
                call system_clock(clock1) 
                ! reset iblank and immersed-boundary variables
                CALL set_iblank_body_fast()
                call system_clock(clock2, clock_rate)
                IF (lProc.eq.PROC_M .and. MOD(ntime,nmonitor)==0 )       &
                   WRITE(*,*) ' Time for setting iblank is:',REAL(clock2-clock1)/REAL(clock_rate)

                call system_clock(clock1)
                CALL GCM_set_internal_boundary()
                call system_clock(clock2, clock_rate)
                IF (lProc.eq.PROC_M .and. MOD(ntime,nmonitor)==0 )       &
                   WRITE(*,*) ' Time for setting boundary is:', REAL(clock2-clock1)/REAL(clock_rate)
             endif
 
             CALL GCM_SetBodyInterceptValues()
             CALL GCM_ghostcell_velocity(1, uOld, vOld, wOld)
             CALL face_vel(uOld,vOld,wOld)
             CALL set_solve_ad()
          END IF ! boundary_motion

500       CONTINUE       ! if RK scheme is used, the substep iteration starts here 

          call MPI_BARRIER(flow_comm,ierr)

          !update rk_step. Necessary for a non-RK scheme, since rk_step is used
          rk_step = rk_step + 1
          !if(advec_scheme == RUNGE_KUTTA3)  CALL INTERPOLATE_RK_BCs(rk_step)

!------------------------------------------------------------------------------
!     Compute advection-diffusion terms
!------------------------------------------------------------------------------

          !call MPI_BARRIER(flow_comm,ierr)
          IF(lProc.eq.PROC_M) call system_clock(clock1)  
                                        
          if(lProc==proc_m) WRITE(*,*) ' SOLVE_AD ... '      !            n         n
          CALL rhs_advec_diff()                              ! compute NLU   and VIS  

          if(lProc==proc_m) WRITE(*,*) ' rhs_advec_diff OK... '
          CALL solve_ad()
          if(lProc==proc_m) WRITE(*,*) ' solve_ad OK... '

          !Ye,debug
          !call write_flowfield_test()
          !if(lProc==proc_m) WRITE(*,*) ' solve_ad OOK... '
          !call MPI_BARRIER(flow_comm,ierr)
          !stop

          !call MPI_BARRIER(flow_comm,ierr)
          CALL face_vel(u,v,w)
          if(lProc==proc_m) WRITE(*,*) ' face_vel OK... '

          !  if(advec_scheme == RUNGE_KUTTA3 .and. rk_step<3) go to 500

          IF(lProc.eq.PROC_M) then
             call system_clock(clock2, clock_rate)

             IF ( MOD(ntime,nmonitor) == 0 ) WRITE(*,*) ' Total time for solve_ad is:', &
                                                         REAL(clock2-clock1)/REAL(clock_rate)
             !WRITE(4,*) ' Total time for solve_ad is:', &
             !            REAL(clock2-clock1)/REAL(clock_rate)
          ENDIF
      
!------------------------------------------------------------------------------
!     Compute RHS for the Poisson Pressure Equation (PPE)
!------------------------------------------------------------------------------
          IF(lProc.eq.PROC_M) call system_clock(clock1)
          CALL rhs_poisson(sumd)

!------------------------------------------------------------------------------
!    Solve the Poisson Pressure Equation (PPE)
!------------------------------------------------------------------------------

          IF(lProc.eq.PROC_M) write(*,*) ' SOLVE_POISSON ... '
          CALL solve_poisson()
    
          IF(lProc.eq.PROC_M) then
             call system_clock(clock2, clock_rate)
             IF ( MOD(ntime,nmonitor) == 0 ) & 
                  WRITE(*,*) ' Total time for solve_poisson is:', &
                              REAL(clock2-clock1)/REAL(clock_rate)           
          ENDIF

    !call MPI_BARRIER(flow_comm,ierr)
    !!Ye,debug
    !write(800+lProc,*)'VARIABLES="X","Y","Z","u","v","w","P","xbi","ybi","zbi","IBK","GC","DC","BN"'
    !write(800+lProc,*)'ZONE F=POINT, I=',nx+1,', J=',jSlices+2,', K=',kSlices+2
    !do k=zc_start-1,zc_end+1   !1,nz-1
    !do j=yc_start-1,yc_end+1
    !do i=0,nx
    ! write(800+lProc,'(10(2X,F12.5),4(3X,I10))')xc(i),yc(j),zc(k),u(i,j,k),v(i,j,k), &
    !                                           w(i,j,k),pPrime(i,j,k), &
    !                                           xBItable(i,j,k),yBItable(i,j,k), &
    !                                           zBItable(i,j,k), &
    !                                           iblank(i,j,k),ghostCellMark(i,j,k), &
    !                                           dead_cell(i,j,k),bodyNum(i,j,k)
    !enddo
    !enddo
    !enddo
    !close(800+lProc)
    !call MPI_BARRIER(flow_comm,ierr)
    !STOP

!------------------------------------------------------------------------------
!    Correct velocity field and update pressure
!------------------------------------------------------------------------------

          CALL correct_vel()
          CALL correct_ghostcell_velocity(pPrime)
          !Update buffer slices for each subdomain
          call send_receive_slices_real_y(u, 0,nx+1,yb1,yb2,zb1,zb2,2)
          call send_receive_slices_real_z(u, 0,nx+1,yb1,yb2,zb1,zb2,2)
          call send_receive_slices_real_y(v, 0,nx+1,yb1,yb2,zb1,zb2,2)
          call send_receive_slices_real_z(v, 0,nx+1,yb1,yb2,zb1,zb2,2)
          call send_receive_slices_real_y(w, 0,nx+1,yb1,yb2,zb1,zb2,2)
          call send_receive_slices_real_z(w, 0,nx+1,yb1,yb2,zb1,zb2,2)

          CALL face_vel(u,v,w)
          CALL update_pressure()

          CALL compute_marker_stress()

          call MPI_BARRIER(flow_comm,ierr)

          call system_clock(clock4, clock_rate)
          IF (lProc==PROC_M .and. MOD(ntime,nmonitor)==0 )       &
             WRITE(*,*) ' Total time for flow (within FSI) is:', REAL(clock4-clock3)/REAL(clock_rate)

          iter_FSI = iter_FSI + 1
!------------------------------------------------------------------------------
!    Monitor output
!------------------------------------------------------------------------------

          IF ( MOD(ntime,nmonitor) == 0 ) THEN
             CALL write_monitor()
          ENDIF

          IF ( boundary_motion == MOVING_BOUNDARY .AND. ntime >= 1) THEN
             IF(lProc==PROC_M) PRINT*,'CALL move_boundary()'   
             call system_clock(clock1)
             CALL move_boundary()    ! Update the boundary position and velocity
             !Ye, for ALL processors
             !save the current resMax_FSIp for use in next fsi iteration 
             resMax_FSIp_old = resMax_FSIp
             call system_clock(clock2, clock_rate)
             IF (lProc==PROC_M .and. MOD(ntime,nmonitor)==0 )       &
                WRITE(*,*) ' Time for move_boundary is:', REAL(clock2-clock1)/REAL(clock_rate)
          ENDIF

          call MPI_BARRIER(flow_comm,ierr)

         !Ye,debug
         !IF ((ntime.eq.9401).and.(iter_FSI.eq.2))THEN
         !   CALL write_dump()
         !ENDIF

          !Ye,debug use
          if(lProc.eq.PROC_M .and. Imonitor_probe == 1) CALL write_probe_files()


!----------------------------
       ENDDO   ! end FSI loop
!---------------------------- 

      IF( lProc.eq.PROC_M) THEN
         WRITE(*,*)'---------------------------------------------------------'
         WRITE(*,*)'===============FSI CONVERGED============================='
         WRITE(*,*)'---------------------------------------------------------'
      END IF

     !======================================================================
     !Ye,receive contact force from solid solver
     if(lProc .eq. PROC_M) then

        cnt_d12 = 1.0E+2
        cnt_d21 = 1.0E+2
        cnt_d13 = 1.0E+2
        cnt_d31 = 1.0E+2
        cnt_d23 = 1.0E+2
        cnt_d32 = 1.0E+2

        cnt_EM12 = -1
        cnt_EM21 = -1
        cnt_EM13 = -1
        cnt_EM31 = -1
        cnt_EM23 = -1
        cnt_EM32 = -1

     do ie=1, nelastic

     CALL MPI_RECV(nPtsMax_s, 1, MPI_INTEGER, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     if (ie.eq.1) then
        allocate(cntflg12(nPtsMax_s,1))
        allocate(cntflg21(nPtsMax_s,1))
        allocate(cntflg13(nPtsMax_s,1))
        allocate(cntflg31(nPtsMax_s,1))
        allocate(cntflg23(nPtsMax_s,1))
        allocate(cntflg32(nPtsMax_s,1))

        allocate(cfx(nPtsMax_s,1))
        allocate(cfy(nPtsMax_s,1))
        allocate(cfz(nPtsMax_s,1))

        allocate(cntd12(nPtsMax_s))
        allocate(cntd21(nPtsMax_s))
        allocate(cntd13(nPtsMax_s))
        allocate(cntd31(nPtsMax_s))
        allocate(cntd23(nPtsMax_s))
        allocate(cntd32(nPtsMax_s))

        allocate(ElemNum12(nPtsMax_s))
        allocate(ElemNum21(nPtsMax_s))
        allocate(ElemNum13(nPtsMax_s))
        allocate(ElemNum31(nPtsMax_s))
        allocate(ElemNum23(nPtsMax_s))
        allocate(ElemNum32(nPtsMax_s))
     endif

     cntflg12 = 0
     cntflg21 = 0
     cntflg13 = 0
     cntflg31 = 0
     cntflg23 = 0
     cntflg32 = 0

     cfx = 0.0_CGREAL
     cfy = 0.0_CGREAL
     cfz = 0.0_CGREAL

     cntd12 = 1.0E+2
     cntd21 = 1.0E+2
     cntd13 = 1.0E+2
     cntd31 = 1.0E+2
     cntd23 = 1.0E+2
     cntd32 = 1.0E+2

     ElemNum12 = -1
     ElemNum21 = -1
     ElemNum13 = -1
     ElemNum31 = -1
     ElemNum23 = -1
     ElemNum32 = -1

     CALL MPI_RECV(cntflg12, nPtsMax_s, MPI_INTEGER, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(cntflg21, nPtsMax_s, MPI_INTEGER, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(cntflg13, nPtsMax_s, MPI_INTEGER, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(cntflg31, nPtsMax_s, MPI_INTEGER, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(cntflg23, nPtsMax_s, MPI_INTEGER, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(cntflg32, nPtsMax_s, MPI_INTEGER, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(cfx, nPtsMax_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(cfy, nPtsMax_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(cfz, nPtsMax_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(cntd12, nPtsMax_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(cntd21, nPtsMax_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(cntd13, nPtsMax_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(cntd31, nPtsMax_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(cntd23, nPtsMax_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(cntd32, nPtsMax_s, MPI_DOUBLE_PRECISION, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(ElemNum12, nPtsMax_s, MPI_INTEGER, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(ElemNum21, nPtsMax_s, MPI_INTEGER, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(ElemNum13, nPtsMax_s, MPI_INTEGER, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(ElemNum31, nPtsMax_s, MPI_INTEGER, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(ElemNum23, nPtsMax_s, MPI_INTEGER, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     CALL MPI_RECV(ElemNum32, nPtsMax_s, MPI_INTEGER, proc_s(ie), 1, &
                  MPI_COMM_WORLD,istatus,ierr)

     iBody = bodyInFlowID(ie)

     DO m=1, nPtsBodyMarker(iBody)
       i = m 
       cnt_fg12(m,iBody)  = cntflg12(i,1)
       cnt_fg21(m,iBody)  = cntflg21(i,1)
       cnt_fg13(m,iBody)  = cntflg13(i,1)
       cnt_fg31(m,iBody)  = cntflg31(i,1)
       cnt_fg23(m,iBody)  = cntflg23(i,1)
       cnt_fg32(m,iBody)  = cntflg32(i,1)

       cnt_fx(m,iBody)  = cfx(i,1)
       cnt_fy(m,iBody)  = cfy(i,1)
       cnt_fz(m,iBody)  = cfz(i,1)

       cnt_d12(m,iBody) = cntd12(i)
       cnt_d21(m,iBody) = cntd21(i)
       cnt_d13(m,iBody) = cntd13(i)
       cnt_d31(m,iBody) = cntd31(i)
       cnt_d23(m,iBody) = cntd23(i)
       cnt_d32(m,iBody) = cntd32(i)

       cnt_EM12(m,iBody) = ElemNum12(i)
       cnt_EM21(m,iBody) = ElemNum21(i)
       cnt_EM13(m,iBody) = ElemNum13(i)
       cnt_EM31(m,iBody) = ElemNum31(i)
       cnt_EM23(m,iBody) = ElemNum23(i)
       cnt_EM32(m,iBody) = ElemNum32(i)
     ENDDO ! m

     enddo!ie
     !write(*,*),'Instad: Contact force received !'
     endif!lProc

     ndata = nPtsMax*nBody
     call MPI_BCAST(cnt_fx, ndata, MPI_DOUBLE_PRECISION,proc_m,flow_comm, ierr)
     call MPI_BCAST(cnt_fy, ndata, MPI_DOUBLE_PRECISION,proc_m,flow_comm, ierr)
     call MPI_BCAST(cnt_fz, ndata, MPI_DOUBLE_PRECISION,proc_m,flow_comm, ierr)

     IF (lProc.eq.PROC_M) THEN
        DO iBody=1, nBody
           Do m=1, nPtsBodyMarker(iBody)
              cnt_force(m,ibody)=sqrt(cnt_fx(m,iBody)**2+cnt_fy(m,iBody)**2+  &
                                      cnt_fz(m,iBody)**2)
           ENDDO
        ENDDO
     ENDIF
     
     call MPI_BCAST(cnt_force, ndata, MPI_DOUBLE_PRECISION,proc_m,flow_comm, ierr)
     !write(*,*),'Instad: Total contact force calculated !'
     !======================================================================

       IF ( MOD(ntime,nmonitor_probe_liftdrag) == 0 ) THEN
          if(lProc.eq.PROC_M .and. Imonitor_drag  == 1) CALL write_total_force()
          !if(lProc.eq.PROC_M .and. Imonitor_probe == 1) CALL write_probe_files()
       ENDIF
       IF ( MOD(ntime,ndump) == 0)    THEN
          CALL write_dump()
       ENDIF

       IF ( MOD(ntime,nrestart) == 0 .OR. ntime==ntime_start+no_tsteps) THEN
          CALL write_restart()
          !CALL write_res()
       ENDIF
       call MPI_BARRIER(flow_comm,ierr)

       IF(lProc.eq.PROC_M) then
          call system_clock(clock5, clock_rate)

          IF ( MOD(ntime,nmonitor) == 0 ) & 
               WRITE(*,*) ' Total time for entire step is:', &
                           REAL(clock5-clock0)/REAL(clock_rate)
               !WRITE(4,*) ' Total time for timestep is:', &
               !          REAL(clock5-clock0)/REAL(clock_rate)
               !WRITE(4,121) NTIME,TIME,  REAL(clock2-clock1)/REAL(clock_rate), &
               !                     REAL(clock3-clock2)/REAL(clock_rate), &
               !                     REAL(clock5-clock0)/REAL(clock_rate)
121       FORMAT(1X,I10,5(3x,1PE12.5))

       ENDIF

       if( (lProc .eq. PROC_M).and.(elastic_present .eqv. .true.) )then
         deallocate(cntflg12)
         deallocate(cntflg21)
         deallocate(cntflg13)
         deallocate(cntflg31)
         deallocate(cntflg23)
         deallocate(cntflg32)

         deallocate(cfx)
         deallocate(cfy)
         deallocate(cfz)

         deallocate(cntd12)
         deallocate(cntd21)
         deallocate(cntd13)
         deallocate(cntd31)
         deallocate(cntd23)
         deallocate(cntd32)

         deallocate(ElemNum12)
         deallocate(ElemNum21)
         deallocate(ElemNum13)
         deallocate(ElemNum31)
         deallocate(ElemNum23)
         deallocate(ElemNum32)
       endif
!--------
    ENDDO ! end of time loop
!--------
      IF( lProc.eq.PROC_M) THEN
         WRITE(*,*)'---------------------------------------------------------'
         WRITE(*,*)'===============TIME STEP FINISHED ======================='
         WRITE(*,*)'---------------------------------------------------------'
      END IF

    ntime = ntime_start+no_tsteps 
            
    !CALL write_dump()                     
    
   END SUBROUTINE time_step_viscous


   !Ye,calculate the volumetric flow rate=======================
   SUBROUTINE calc_flow_rate()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE GCM_arrays
    USE MPI_module

    IMPLICIT NONE

    INTEGER             :: i,j,k
    REAL(KIND=CGREAL)   :: flag

    vq1_i   = 0.0_CGREAL
    vq1     = 0.0_CGREAL

    vq2_i   = 0.0_CGREAL
    vq2     = 0.0_CGREAL

    vq3_i   = 0.0_CGREAL
    vq3     = 0.0_CGREAL

    vq4_i   = 0.0_CGREAL
    vq4     = 0.0_CGREAL

    i = 5
    DO k=zc_start,zc_end
    DO j=yc_start,yc_end
       flag = (1.0_CGREAL-REAL(iblank       (i,j,k),KIND=CGREAL))  
       !     * (1.0_CGREAL-REAL(ghostcellMark(i,j,k),KIND=CGREAL))

       vq1_i = vq1_i + flag * face_u(i+1,j,k)*dy(j)*dz(k)
    ENDDO
    ENDDO

    CALL MPI_ALLREDUCE(vq1_i, vq1, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, flow_comm, ierr)

    i = 100
    DO k=zc_start,zc_end
    DO j=yc_start,yc_end
       flag = (1.0_CGREAL-REAL(iblank       (i,j,k),KIND=CGREAL))  
       !     * (1.0_CGREAL-REAL(ghostcellMark(i,j,k),KIND=CGREAL))

       vq2_i = vq2_i + flag * face_u(i+1,j,k)*dy(j)*dz(k)
    ENDDO
    ENDDO

    CALL MPI_ALLREDUCE(vq2_i, vq2, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, flow_comm, ierr)

    i = 170
    DO k=zc_start,zc_end
    DO j=yc_start,yc_end
       flag = (1.0_CGREAL-REAL(iblank       (i,j,k),KIND=CGREAL))
       !     * (1.0_CGREAL-REAL(ghostcellMark(i,j,k),KIND=CGREAL))

       vq3_i = vq3_i + flag * face_u(i+1,j,k)*dy(j)*dz(k)
    ENDDO
    ENDDO

    CALL MPI_ALLREDUCE(vq3_i, vq3, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, flow_comm, ierr)

    i = nx-4
    DO k=zc_start,zc_end
    DO j=yc_start,yc_end
       flag = (1.0_CGREAL-REAL(iblank       (i,j,k),KIND=CGREAL))
       !     * (1.0_CGREAL-REAL(ghostcellMark(i,j,k),KIND=CGREAL))

       vq4_i = vq4_i + flag * face_u(i+1,j,k)*dy(j)*dz(k)
    ENDDO
    ENDDO

    CALL MPI_ALLREDUCE(vq4_i, vq4, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, flow_comm, ierr)

   END SUBROUTINE calc_flow_rate


   !Ye,Calculate GOA=================================
   SUBROUTINE calc_goa()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE GCM_arrays
    USE MPI_module

    IMPLICIT NONE

    INTEGER             :: i,j,k

    goa = 0
    goa_area_i = 0.0d0
    goa_area = 0.0d0

    DO k=zc_start,zc_end
    DO j=yc_start,yc_end
       DO i=2,nx-2
          if (iblank(i,j,k).eq.1)then
             goa(1,j,k) = 1
          endif
       ENDDO
    ENDDO
    ENDDO

    DO k=zc_start,zc_end
       DO j=yc_start,yc_end
          if (goa(1,j,k).eq.0)then
             goa_area_i = goa_area_i + dyc(j)*dzc(k)
          endif
       ENDDO
    ENDDO

    CALL MPI_ALLREDUCE(goa_area_i, goa_area, 1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, flow_comm, ierr)

   END SUBROUTINE calc_goa

   !Ye, write flowrate and GOA to file==================================
   SUBROUTINE write_flowrate()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE multiuse_arrays
    USE stat_arrays
    USE stat_vort_arrays
    USE fsi_module
    USE MPI_module
    USE GCM_arrays
    USE solver_arrays

    IMPLICIT NONE

    write(9999,'(6F15.8)') time, vq1, vq2, vq3, vq4, goa_area

   END
