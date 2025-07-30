   SUBROUTINE GCM_ghostcell_velocity(iOpt, u_vel,v_vel,w_vel) 
!
!  update u*, v* and w* values for ghost cell
!
!
!  Ghost cell velocity satisfies the following equations
!
!  Integral [ DEL.(u*) dv ] =- Uo.nDS  with coeff of "dead" faces = 0
!
! [                                                               ]
! [ U  +  (imagePointWeight) U   =   U  * (bodyInterceptWeight)   ] . tau
! [  gp                       ip      b                           ]   ---
!
! and
!        (          ) 
!   U  = (coeffGCMD ) U
!    ip  (         i)  i


    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE flow_arrays
    USE grid_arrays
    USE GCM_arrays
    USE multiuse_arrays
    use mpi_module

    use fsi_module

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iOpt   ! iOpt = 1, update the dead cell only
                                  ! iOpt = 2, update both dead and hybrid cells
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2), INTENT(INOUT) :: u_vel
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2), INTENT(INOUT) :: v_vel
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2), INTENT(INOUT) :: w_vel

    INTEGER           :: iBody, iRow, n, iterGC, loc(3)
    INTEGER           :: i,j,k,ii,jj,kk
    REAL(KIND=CGREAL) :: resVel,resVelMax,uPrev,vPrev,wPrev,res_i
    REAL(KIND=CGREAL) :: uIP,vIP,wIP
    REAL(KIND=CGREAL) :: bmx,bpx,bcx,bmy,bpy,bcy,bmz,bcz,bpz
    REAL(KIND=CGREAL) :: uG,vG,wG,flag,uTemp,vTemp,wTemp
    REAL(KIND=CGREAL) :: rhsu, rhsv, rhsw    
    REAL(KIND=CGREAL) :: up, vp, wp

    INTEGER           :: degenerate_dc
    INTEGER           :: degenerate_gc 

    REAL(KIND=CGREAL) :: u_inv_dist_int, v_inv_dist_int, w_inv_dist_int
    REAL(KIND=CGREAL) :: dist, dist_int 
    REAL(KIND=CGREAL) :: xIP,yIP,zIP

    INTEGER           :: sumk
    REAL(KIND=CGREAL) :: res_ip(2),res_loc(2),res2
    INTEGER           :: ip_loc

! Iterate to correct interdependent ghost points

    ! Initialize values
    iterGC    = 0
    resVelMax = 1.0E10_CGREAL
    loc  =  1

    !bodyInterceptWeight = (probeLengthNormalized-1.0_CGREAL) / probeLengthNormalized
    !imagePointWeight    = 1.0_CGREAL/probeLengthNormalized   

    bodyInterceptWeight = probeLengthNormalizedD / (probeLengthNormalizedD-1.0_CGREAL)
    imagePointWeight    = 1.0_CGREAL/(probeLengthNormalizedD-1.0_CGREAL)

    DO WHILE ((iterGC .LT. itermax_gc) .AND. (resVelMax .GT. restol_gcv))

       resVelMax = 0.0_CGREAL

       !--------------------
       ! update the dead cells
       !--------------------
       DO n = 1, nDead

          i=iDead(n)
          j=jDead(n)
          k=kDead(n)

          xIP = xyzDCIP(i,j,k,0)
          yIP = xyzDCIP(i,j,k,1)
          zIP = xyzDCIP(i,j,k,2)

          !skip the dead cells in the buffer layer
          if(j<yc_start .or. j>yc_end) cycle
          if(k<zc_start .or. k>zc_end) cycle

          uIP = 0.0_CGREAL
          vIP = 0.0_CGREAL
          wIP = 0.0_CGREAL

          uPrev = u_vel(i,j,k)
          vPrev = v_vel(i,j,k)
          wPrev = w_vel(i,j,k)

          !Ye,inverse distance sqr weighted inetrpolation for DC
          u_inv_dist_int = 0.0_CGREAL
          v_inv_dist_int = 0.0_CGREAL
          w_inv_dist_int = 0.0_CGREAL
          dist_int       = 0.0_CGREAL

          sumk = 0
          !sum_valid(i,j,k) = 0
          DO iRow = 1, iRowMax
             ii = iDeadCellIndex(n) + incI(iRow)
             jj = jDeadCellIndex(n) + incJ(iRow)
             kk = kDeadCellIndex(n) + incK(iRow)

             !Ye,debug
             !if(ntime==8964 .and. lProc==87 .and. i==147 .and. j==94 .and.k==65)then
             !   write(4000+lProc,*)iRow,ii,jj,kk,iblank(ii,jj,kk),dead_cell(ii,jj,kk),&
             !            ghostCellMark(ii,jj,kk)
             !endif

             if ( (ii.eq.i).and.(jj.eq.j).and.(kk.eq.k) ) cycle

             dist = (xc(ii)-xIP)**2 &
                  + (yc(jj)-yIP)**2 &
                  + (zc(kk)-zIP)**2

             ! include in the interpolation if the cell is a fluid cell or dead
             ! cell
             IF ( iblank(ii,jj,kk) == 0 .OR. dead_cell(ii,jj,kk) == -1 ) THEN
                sumk = sumk + 1
                !sum_valid(i,j,k) = sum_valid(i,j,k) + 1
                dist_int       = dist_int + (1.0_CGREAL/dist)
                u_inv_dist_int = u_inv_dist_int + (1.0_CGREAL/dist)*u(ii,jj,kk)
                v_inv_dist_int = v_inv_dist_int + (1.0_CGREAL/dist)*v(ii,jj,kk)
                w_inv_dist_int = w_inv_dist_int + (1.0_CGREAL/dist)*w(ii,jj,kk)
             ENDIF
          ENDDO

        !Ye,debug
        !if(ntime==8964 .and. lProc==87 .and. i==147 .and. j==94 .and.k==65)then
        !  write(3000,*)n,sumk,dist_int
        !  write(3000,*)xIP,yIP,zIP
        !  write(3000,*)iDeadCellIndex(n),jDeadCellIndex(n),kDeadCellIndex(n)
        !  write(3000,*)'---------'
        !endif 

          DO iRow = 1, iRowMax
             ii = iDeadCellIndex(n) + incI(iRow)
             jj = jDeadCellIndex(n) + incJ(iRow)
             kk = kDeadCellIndex(n) + incK(iRow)

             !If the stencil involves iblank cells where flow variables are not
             !defined,
             !we use inverse-distance averaged velocity.
             IF( iblank(ii,jj,kk) == 1 .AND. dead_cell(ii,jj,kk)== 0) THEN
               if (dist_int < 0.000001d0 .or. sumk==1) then ! no/1 point is found for the interpolation
                  !write(*,'(6F12.5)') xc(ii), yc(jj), zc(kk), dist_int
                  uIP = 0.0_CGREAL
                  vIP = 0.0_CGREAL
                  wIP = 0.0_CGREAL
               else
                  uIP = u_inv_dist_int/dist_int
                  vIP = v_inv_dist_int/dist_int
                  wIP = w_inv_dist_int/dist_int
                  GOTO 566
               endif
             ENDIF
          ENDDO

          DO iRow = 1, iRowMax

             ii = iDeadCellIndex(n) + incI(iRow)
             jj = jDeadCellIndex(n) + incJ(iRow)
             kk = kDeadCellIndex(n) + incK(iRow)

             IF ( ii /= i .OR. jj /= j .OR. kk /= k) THEN
                uIP = uIP + coeffGCMDeadD(iRow,n)* u_vel(ii,jj,kk)
                vIP = vIP + coeffGCMDeadD(iRow,n)* v_vel(ii,jj,kk)
                wIP = wIP + coeffGCMDeadD(iRow,n)* w_vel(ii,jj,kk)
             ELSE
                uIP = uIP + coeffGCMDeadD(iRow,n)* uBodyInterceptDead(n)
                vIP = vIP + coeffGCMDeadD(iRow,n)* vBodyInterceptDead(n)
                wIP = wIP + coeffGCMDeadD(iRow,n)* wBodyInterceptDead(n)
             ENDIF ! ii

          ENDDO ! iRow

566   CONTINUE

          u_vel(i,j,k) =  uBodyInterceptDead(n)*bodyInterceptWeight - uIP*imagePointWeight
          v_vel(i,j,k) =  vBodyInterceptDead(n)*bodyInterceptWeight - vIP*imagePointWeight
          w_vel(i,j,k) =  wBodyInterceptDead(n)*bodyInterceptWeight - wIP*imagePointWeight

          ! Compute residual
          resVel = ABS( u_vel(i,j,k)-uPrev ) + ABS( v_vel(i,j,k)-vPrev )  &
                 + ABS( w_vel(i,j,k)-wPrev )

          if (resVel > resVelMax ) then
             resVelMax = resVel
             loc(1) = i; loc(2) = j; loc(3) = k
             up = uPrev; vp = vPrev; wp = wPrev;
          endif

       ENDDO ! n

       if(iOpt == 1) goto 100   ! skip the hybrid node update

       !--------------------
       ! update the ghost cells
       !--------------------
       DO n = 1, nGhost

          i=iGhost(n)
          j=jGhost(n)
          k=kGhost(n)

          uG = 0.0_CGREAL
          vG = 0.0_CGREAL
          wG = 0.0_CGREAL

          uPrev = u(i,j,k)
          vPrev = v(i,j,k)
          wPrev = w(i,j,k)

          DO iRow = 1, iRowMax
            ii = iCellIndex(n) + incI(iRow)
            jj = jCellIndex(n) + incJ(iRow)
            kk = kCellIndex(n) + incK(iRow)

             IF ( ii /= i .OR. jj /= j .OR. kk /= k) THEN
                uG = uG + coeffGCMD(iRow,n)* u(ii,jj,kk)
                vG = vG + coeffGCMD(iRow,n)* v(ii,jj,kk)
                wG = wG + coeffGCMD(iRow,n)* w(ii,jj,kk)
             ELSE
                uG = uG + coeffGCMD(iRow,n)* uBodyIntercept(n)
                vG = vG + coeffGCMD(iRow,n)* vBodyIntercept(n)
                wG = wG + coeffGCMD(iRow,n)* wBodyIntercept(n)
             ENDIF ! ii
          ENDDO ! iRow

          !-------------------
          ! Use a hybrid formulation to update ghost-node.
          ! This is done during time-stepping of the NS equation.
          !-------------------
          flag  = alphaGhost(n)

          if(mix_GC_form == 1) then    
             bcx = ghostNScoeff(1,n)
             bmx = ghostNScoeff(2,n)
             bpx = ghostNScoeff(3,n)
             bmy = ghostNScoeff(4,n)
             bpy = ghostNScoeff(5,n)
             bmz = ghostNScoeff(6,n)
             bpz = ghostNScoeff(7,n)

             rhsu= ghostNScoeff(8 ,n)
             rhsv= ghostNScoeff(9 ,n)
             rhsw= ghostNScoeff(10,n)

             uTemp = (rhsu - bmx*u(i-1,j,k) - bpx*u(i+1,j,k)  &
                           - bmy*u(i,j-1,k) - bpy*u(i,j+1,k)  &
                           - bmz*u(i,j,k-1) - bpz*u(i,j,k+1) ) / bcx

             vTemp = (rhsv - bmx*v(i-1,j,k) - bpx*v(i+1,j,k)    &
                           - bmy*v(i,j-1,k) - bpy*v(i,j+1,k)    & 
                           - bmz*v(i,j,k-1) - bpz*v(i,j,k+1) ) / bcx

             wTemp = (rhsw - bmx*w(i-1,j,k) - bpx*w(i+1,j,k)    &
                           - bmy*w(i,j-1,k) - bpy*w(i,j+1,k)    & 
                           - bmz*w(i,j,k-1) - bpz*w(i,j,k+1) ) / bcx

          else  ! mix_GC_form
             ! otherwise, use the interpolation value
             uTemp = uG
             vTemp = vG
             wTemp = wG
          endif ! mix_GC_form

          u(i,j,k) = (1.0_CGREAL-flag)*uG + flag*uTemp
          v(i,j,k) = (1.0_CGREAL-flag)*vG + flag*vTemp
          w(i,j,k) = (1.0_CGREAL-flag)*wG + flag*wTemp
          !write(*,'(I5,8F12.5)') n, uG, uTemp, u(i,j,k), flag

          ! Compute residual
          resVel = ABS( u(i,j,k)-uPrev ) + ABS( v(i,j,k)-vPrev )  &
                 + ABS( w(i,j,k)-wPrev )

          if (resVel > resVelMax ) then
             resVelMax = resVel
             loc(1) = i; loc(2) = j; loc(3) = k
             up = uPrev; vp = vPrev; wp = wPrev;
          endif
       ENDDO ! n

100    continue

       !write(*,'(I5,E12.5,3I6,6F12.5)')lProc,resVelMax,loc(1:3),u(loc(1),loc(2),loc(3)),up,vp,wp
       !write(*,'(I5,E12.5,3I6,3F12.5)')lProc,resVelMax,loc(1:3),xc(loc(1)),yc(loc(2)),zc(loc(3))
       !Data exchange between subdomains for 2 slices
       call send_receive_slices_real_y(u, 0,nx+1,yb1,yb2,zb1,zb2,2)
       call send_receive_slices_real_z(u, 0,nx+1,yb1,yb2,zb1,zb2,2)
       call send_receive_slices_real_y(v, 0,nx+1,yb1,yb2,zb1,zb2,2)
       call send_receive_slices_real_z(v, 0,nx+1,yb1,yb2,zb1,zb2,2)
       call send_receive_slices_real_y(w, 0,nx+1,yb1,yb2,zb1,zb2,2)
       call send_receive_slices_real_z(w, 0,nx+1,yb1,yb2,zb1,zb2,2)

       !Find the global max. residual so that all procs. enter the do-while loop
       res_i = resVelMax
       !Ye,debug
       !write(6000+lProc,'(I6,I6,I6,E12.5)')ntime,iter_FSI,iterGC,res_i
       !if(ntime == 8964 .and. abs(res_i-0.2169).le.1.0E-4 )then
       !  write(7000+lProc,'(I6,I6,I6,E12.5)')ntime,lProc,iter_FSI,res_i
       !endif
       !if(ntime == 8964 .and. abs(res_i-0.2169).le.1.0E-4 )then
       !  write(8000+lProc,'(I6,I6,I6,I6)')ntime,loc(1:3)
       !endif
       call MPI_ALLREDUCE(res_i, resVelMax, 1, MPI_DOUBLE_PRECISION,  &
                          MPI_MAX, flow_comm, ierr)

       !Pack up residual and processor ID for comparison
       res_ip(1) = res_i
       res_ip(2) = lProc
       !Find max error among all subdomains and propagate it to all
       !processors.
       CALL MPI_ALLREDUCE(res_ip, res_loc, 1, MPI_2DOUBLE_PRECISION, &
                          MPI_MAXLOC, FLOW_COMM, ierr)
       res2   = res_loc(1)
       ip_loc = res_loc(2)  ! coerce into integer
     
       if(ip_loc /= PROC_M) then  !send the location to PROC_M
          if(lProc==ip_loc) &
               call MPI_SEND(loc, 3, MPI_INTEGER, PROC_M,1,FLOW_COMM,istatus,ierr)
          if(lProc==PROC_M) &
               call MPI_RECV(loc, 3, MPI_INTEGER, ip_loc,1,FLOW_COMM,istatus,ierr)
       endif

       iterGC = iterGC + 1

       IF ( MOD(ntime,nmonitor) == 0 .and. lProc==PROC_M .and. (resVelMax.le.restol_gcv)) then
          write(*,'(a,I5,a,I5,E12.5)')'lProc',lProc,' Ghostcell Velocity Convergence: ',iterGC,resVelMax
       ENDIF

    ENDDO ! iterGC

    IF ( iterGC .EQ. itermax_gc .AND. resVelMax .GT. restol_gcv ) THEN
       PRINT*,'GCM_vel_set_bc_internal for iBody :', iBody
       PRINT*,'   GhostCell u* did not converge in ',itermax_gc,' iterations'
       PRINT*,'   Final residual = ',resVelMax
       IF(lProc==proc_m .and. MOD(ntime,nmonitor) == 0) THEN
         write(6,'(a,I5,a,PE15.5)'),'iter=',iterGC,' GhostVel max residual:',res2
         write(6,'(a,4I5)'),'  Location of max residual:',ip_loc,loc(1),loc(2),loc(3)
       ENDIF
       call write_subdomain()
       if(lProc==PROC_M) call write_markers()
       call MPI_BARRIER(flow_comm,ierr)
       call sleep(15)
       STOP
    ENDIF

  END SUBROUTINE  GCM_GhostCell_Velocity
!----------------------------------------------------------------------
  SUBROUTINE GCM_ghostcell_pressure(pres,r)

!----------------------------------------------------------------------
! Compute Ghost point values at internal boundary for the velocity field
!----------------------------------------------------------------------
!
!     
!   P  =   P   
!    gp     ip 
!
! and
!        (          ) 
!   P  = (coeffGCMN ) P
!    ip  (         i)  i


    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE MPI_module

    USE fsi_module

    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2), INTENT(INOUT) :: pres
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2), INTENT(IN   ) :: r

    INTEGER :: iBody,iRow,n,loc(3)
    INTEGER :: iG,ii,ind,jG,jj,kG,kk,iterGC,i2,j2,k2
    
    REAL(KIND=CGREAL) :: pIP, resPres, resPresMax, pPrev, dpdn, res_i

    REAL(KIND=CGREAL) :: dpdxE, dpdxW, dpdyN, dpdyS, dpdzF, dpdzB
    REAL(KIND=CGREAL) :: coeffG, pG, pG2, flag

    REAL(KIND=CGREAL) :: dist, dist_int, p_inv_dist_int, pIP2
    REAL(KIND=CGREAL) :: xIP, yIP, zIP

    INTEGER :: degenerate_dc, is_change_dc
    INTEGER :: degenerate_gc, is_change_gc

    INTEGER           :: sumk
    REAL(KIND=CGREAL) :: res_ip(2),res_loc(2),res2
    INTEGER           :: ip_loc
!*****************************************************************************************

! Loop over all immersed bodies

    iterGC    = 0
    resPresMax = 1.0E10_CGREAL
    loc  =  1

    ! Iterate to correct interdependent ghost points
    DO WHILE ((iterGC .LT. itermax_gc) .AND. (resPresMax .GT. restol_gcp))

       resPresMax = 0.0_CGREAL

       !--------------
       !Update dead cells
       !--------------
       DO n = 1,nDead

          iG = iDead(n)
          jG = jDead(n)
          kG = kDead(n)

          xIP = xyzDCIP(iG,jG,kG,0)
          yIP = xyzDCIP(iG,jG,kG,1)
          zIP = xyzDCIP(iG,jG,kG,2)

          !skip the dead cells in the buffer layer
          if(jG<yc_start .or. jG>yc_end) cycle
          if(kG<zc_start .or. kG>zc_end) cycle

          ! Save previous values
          pPrev = pres(iG,jG,kG)

          ! Initialize values
          pIP = 0.0_CGREAL

          if(is_PBC_homogeneous == 1) then
             dpdn = 0.0_CGREAL
          else
             dpdn = dpdnInterceptDead(n)
          endif

          !Ye,inverse distance sqr weighted inetrpolation
          p_inv_dist_int = 0.0_CGREAL
          dist_int       = 0.0_CGREAL

          sumk = 0
          !sum_valid(iG,jG,kG) = 0
          DO iRow = 1, iRowMax
             ii = iDeadCellIndex(n) + incI(iRow) 
             jj = jDeadCellIndex(n) + incJ(iRow)
             kk = kDeadCellIndex(n) + incK(iRow)

             !Ye,debug
             IF(ntime==9618 .and. iG==140 .and. jG==94 .and. kG==65)THEN
               write(3000+lProc,*)iRow,ii,jj,kk,iblank(ii,jj,kk),dead_cell(ii,jj,kk),&
                                  ghostCellMark(ii,jj,kk)
             ENDIF

             if ( (ii.eq.iG).and.(jj.eq.jG).and.(kk.eq.kG) ) cycle

             dist = (xc(ii)-xIP)**2 &
                  + (yc(jj)-yIP)**2 &
                  + (zc(kk)-zIP)**2
 
             ! include in the interpolation if the cell is a fluid cell or dead cell
             IF ( iblank(ii,jj,kk) == 0 .OR. dead_cell(ii,jj,kk) == -1 ) THEN
                sumk = sumk + 1
                !sum_valid(iG,jG,kG) = sum_valid(iG,jG,kG) + 1
                dist_int       = dist_int + (1.0_CGREAL/dist)
                p_inv_dist_int = p_inv_dist_int + (1.0_CGREAL/dist)*p(ii,jj,kk)
             ENDIF
          ENDDO

          !Ye,debug
          !IF(ntime==9618 .and. iG==140 .and. jG==94 .and. kG==65)THEN
          !  write(4000+lProc,*)sumK
          !ENDIF

          DO iRow = 1, iRowMax
             ii = iDeadCellIndex(n) + incI(iRow) 
             jj = jDeadCellIndex(n) + incJ(iRow)
             kk = kDeadCellIndex(n) + incK(iRow)

             !If the stencil involves iblank cells where flow variables are not defined,
             !we use inverse-distance averaged pressure.
             IF( iblank(ii,jj,kk) == 1 .AND. dead_cell(ii,jj,kk)== 0) THEN
               if (dist_int < 0.000001d0 .or. sumk==1) then ! no/1 point is found for the interpolation
               !if (dist_int < 0.000001d0) then ! no point is found
                  !write(*,'(6F12.5)') xc(ii), yc(jj), zc(kk), dist_int
                  pIP = 0.0_CGREAL
               else
                  pIP = p_inv_dist_int/dist_int
                  GOTO 567
               endif
             ENDIF
          ENDDO

          ! Compute pressure field at Image points
          DO iRow = 1, iRowMax
             ii = iDeadCellIndex(n) + incI(iRow) 
             jj = jDeadCellIndex(n) + incJ(iRow)
             kk = kDeadCellIndex(n) + incK(iRow)

             IF ( ii == iG .AND. jj == jG .AND. kk == kG ) THEN
                pIP = pIP + coeffGCMDeadN(iRow,n)*dpdn
             ELSE
                pIP = pIP + coeffGCMDeadN(iRow,n)* pres(ii,jj,kk)
             ENDIF ! ii          
          ENDDO ! iRow

567   CONTINUE

          pres(iG,jG,kG) = pIP + dpdn * probeLengthDead(n)

          !Ye,debug
          !if ((iG.eq.144).and.(jG.eq.102).and.(kG.eq.68))then
          !   write(700+lProc,*)ntime,iter_FSI,iterGC,n
          !   write(700+lProc,*)pIP,dpdn,probeLengthDead(n)
          !   write(700+lProc,*)pres(iG,jG,kG)
          !   write(700+lProc,*)'----'
          !endif

          ! Compute residual
          resPres = ABS( pres(iG,jG,kG) - pPrev )
          IF (resPres > resPresMax )THEN
             resPresMax = resPres
             loc(1) = iG; loc(2) = jG; loc(3) = kG
          ENDIF

          DO iRow = 1, iRowMax
             ii = iDeadCellIndex(n) + incI(iRow)
             jj = jDeadCellIndex(n) + incJ(iRow)
             kk = kDeadCellIndex(n) + incK(iRow)

             !Ye,debug===================
             !To see if it contains pure solid cubic nodes or not
             !IF( (iblank(ii,jj,kk).eq.1).AND.  &
             !    (dead_cell(ii,jj,kk).eq.0).AND. &
             !    (resPresMax.le.restol_gcp) )THEN
             !  WRITE(700,*)ntime,iter_FSI,lProc,iterGC,n
             !  WRITE(700,*)iG,jG,kG,ii,jj,kk
             !  WRITE(700,*)iDeadCellIndex(n),jDeadCellIndex(n),  &
             !              kDeadCellIndex(n)
             !  WRITE(700,*)'=============='
             !ENDIF
             !=========================
          ENDDO ! iRow
 
       ENDDO ! n 

       !------------------
       !update ghost cells
       !------------------
       DO n = 1,nGhost

          iG = iGhost(n)
          jG = jGhost(n)
          kG = kGhost(n)

          ! Save previous values
          pPrev = pres(iG,jG,kG)

          ! Initialize values
          pG = 0.0_CGREAL

          if(is_PBC_homogeneous == 1) then
             dpdn = 0.0_CGREAL
          else
             dpdn = dpdnBodyIntercept(n)
          endif

          ! Compute pressure field at Image points
          DO iRow = 1, iRowMax
             ii = iCellIndex(n) + incI(iRow) 
             jj = jCellIndex(n) + incJ(iRow)
             kk = kCellIndex(n) + incK(iRow)

          !Ye,debug,check ghostcell projection
          IF( (ntime.eq.9).AND.(iG.eq.204).AND.(jG.eq.217).AND.  &
              (kG.eq.131).AND.(iterGC.le.109) )THEN
              write(4000+lProc,'(9(I6,1X),3(F12.5,1X))')iter_FSI,iterGC,&
                   iRow,ii,jj,kk, &
                   iblank(ii,jj,kk),dead_cell(ii,jj,kk),ghostCellMark(ii,jj,kk),&
                   coeffGCMN(iRow,n),pres(ii,jj,kk),pG
          ENDIF

             IF ( ii == iG .AND. jj == jG .AND. kk == kG ) THEN
                pG = pG + coeffGCMN(iRow,n)*dpdn
             ELSE
                pG = pG + coeffGCMN(iRow,n)* pres(ii,jj,kk)
             ENDIF ! ii          
          ENDDO ! iRow

       if(ntime == 9 .and. iG==204 .and. jG==217 .and. kG==131)then
         write(6000+lProc,'(4(I6,1X),F12.5)')ntime,n,iter_FSI,iterGC,pG
       endif

          !-------------------
          ! Use a hybrid formulation to update ghost-node.
          ! This is done during time-stepping of the NS equation.
          !-------------------

          if(mix_GC_form == 1) then  

            dpdxW = (pres(iG,jG,kG) - pres(iG-1,jG,kG)) * dxcinv(iG)
            dpdxE = (pres(iG+1,jG,kG) - pres(iG,jG,kG)) * dxcinv(iG+1)
            dpdyS = (pres(iG,jG,kG) - pres(iG,jG-1,kG)) * dycinv(jG)
            dpdyN = (pres(iG,jG+1,kG) - pres(iG,jG,kG)) * dycinv(jG+1) 

            dpdzB = (pres(iG,jG,kG) - pres(iG,jG,kG-1)) * dzcinv(kG)
            dpdzF = (pres(iG,jG,kG+1) - pres(iG,jG,kG)) * dzcinv(kG+1)

            flag  = alphaGhost(n)

            coeffG = -(dxcinv(iG) + dxcinv(iG+1)) * dxinv(iG)   &
                     -(dycinv(jG) + dycinv(jG+1)) * dyinv(jG)   &
                     -(dzcinv(kG) + dzcinv(kG+1)) * dzinv(kG)

            resPres =  (dpdxE - dpdxW) * dxinv(iG)     &
                    +  (dpdyN - dpdyS) * dyinv(jG)     &
                    +  (dpdzF - dpdzB) * dzinv(kG) - r(iG,jG,kG)

            pG2     = pres(iG,jG,kG) + resPres / (-coeffG)
             
            pres(iG,jG,kG) = (1.0_CGREAL - flag)*pG + flag*pG2

            ! residual of the Poisson equation at the ghost node
            res_laplacian(n) = (-coeffG)*(pG2 - pG)*(1.0_CGREAL  - flag)

          ELSE

             pres(iG,jG,kG) = pG          ! use the interpolation only
             res_laplacian(n) = 0.0_CGREAL
          ENDIF
          !write(*,'(I5,8F12.5)') n, pG, pG2, pres(iG,jG,kG), flag

       if(ntime == 9 .and. iG==204 .and. jG==217 .and. kG==131)then
         write(5000+lProc,'(4(I6,1X),2(F12.5,1X))')ntime,lProc,iter_FSI,iterGC,pG2,flag
       endif

          ! Compute residual
          resPres = ABS( pres(iG,jG,kG) - pPrev )
          IF (resPres > resPresMax )THEN
             resPresMax = resPres
             loc(1) = iG; loc(2) = jG; loc(3) = kG
          ENDIF

          !Ye,debug===================
          !To see if it contains pure solid cubic nodes or not
          !DO iRow = 1, iRowMax
          !   ii = iCellIndex(n) + incI(iRow)
          !   jj = jCellIndex(n) + incJ(iRow)
          !   kk = kCellIndex(n) + incK(iRow)
             !IF( (iblank(ii,jj,kk).eq.1).AND.  &
             !    (dead_cell(ii,jj,kk).eq.0) .AND.  &
             !    (resPresMax.le.restol_gcp) )THEN
             !  WRITE(800,*)ntime,iter_FSI,lProc,iterGC,n
             !  WRITE(800,*)iG,jG,kG,ii,jj,kk
             !  WRITE(800,*)iCellIndex(n),jCellIndex(n),  &
             !              kCellIndex(n)
             !  WRITE(800,*)'=============='
             !ENDIF
          !ENDDO ! iRow
          !=========================

          !Ye,debug,check ghostcell projection
          !IF( (ntime.eq.4804).AND.(iG.eq.81).AND.(jG.eq.81).AND.  &
          !    (kG.eq.50).AND.(resPresMax.le.restol_gcp).AND.  &
          !    (iter_FSI.eq.2).AND.(lProc.eq.4) )THEN
          !  WRITE(807,*)ntime,iter_FSI,lProc,n,iterGC
          !  WRITE(807,*)iCellIndex(n),jCellIndex(n),  &
          !              kCellIndex(n)
          !  WRITE(807,*)pres(iG,jG,kG),dpdn
          !  WRITE(807,*)pG, pG2, flag
          !  WRITE(807,*)coeffGCMN(1:8,n)
          !  WRITE(807,*)'========'
          !ENDIF

       ENDDO ! n    

       !Data exchange between subdomains for 2 slices
       call send_receive_slices_real_y(pres, 0,nx+1,yb1,yb2,zb1,zb2,2)
       call send_receive_slices_real_z(pres, 0,nx+1,yb1,yb2,zb1,zb2,2)

       !Find the global max. residual so that all procs. enter the do-while loop
       res_i = resPresMax
       call MPI_ALLREDUCE(res_i, resPresMax, 1, MPI_DOUBLE_PRECISION,  &
                          MPI_MAX, flow_comm, ierr)

       !Pack up residual and processor ID for comparison
       res_ip(1) = res_i
       res_ip(2) = lProc
       !Find max error among all subdomains and propagate it to all
       !processors.
       CALL MPI_ALLREDUCE(res_ip, res_loc, 1, MPI_2DOUBLE_PRECISION, &
                          MPI_MAXLOC, FLOW_COMM, ierr)
       res2   = res_loc(1)
       ip_loc = res_loc(2)  ! coerce into integer

       if(ip_loc /= PROC_M) then  !send the location to PROC_M
          if(lProc==ip_loc) &
               call MPI_SEND(loc, 3, MPI_INTEGER, PROC_M,1,FLOW_COMM,istatus,ierr)
          if(lProc==PROC_M) &
               call MPI_RECV(loc, 3, MPI_INTEGER, ip_loc,1,FLOW_COMM,istatus,ierr)
       endif

       !if(ntime == 9 )then
       !  write(7000+lProc,'(I6,I6,I6,E12.5)')ntime,lProc,iter_FSI,res_i
       !endif
       !if(ntime == 9 )then
       !  write(8000+lProc,'(I6,I6,I6,I6)')ntime,loc(1:3)
       !endif

       iterGC = iterGC + 1

       !IF (MOD(ntime,nmonitor) == 0 .and. lProc==PROC_M)  &
       !     write(6,'(I5,a,I6,3x,E12.5)')lProc,' GhostCell pressure convergence :',iterGC,resPresMax

    ENDDO ! iterGC

    IF ( iterGC .EQ. itermax_gc .AND. resPresMax .GT. restol_gcp ) THEN
       PRINT*,'GhostCell Pressure did not converge in ',itermax_gc,' iterations'
       PRINT*,'Final residual = ',resPresMax
       IF(lProc==proc_m .and. MOD(ntime,nmonitor) == 0) THEN
         write(6,'(a,I5,a,PE15.5)'),'iter=',iterGC,' GhostPressure max residual:',res2
         write(6,'(a,4I5)'),'  Location of max residual:',ip_loc,loc(1),loc(2),loc(3)
       ENDIF
    ENDIF
    
  END SUBROUTINE GCM_ghostcell_pressure
!----------------------------------------------------------------------

  SUBROUTINE correct_ghost_massflux(massflux)

    USE global_parameters
    USE flow_parameters
    USE boundary_arrays
    USE flow_arrays
    USE grid_arrays
    USE GCM_arrays
    USE multiuse_arrays

    IMPLICIT NONE

    REAL(KIND=CGREAL), INTENT(INOUT) :: massflux

!... Local variables

    INTEGER           :: i,j,k, n
    REAL(KIND=CGREAL) :: flag
!****************************************************************

    DO n = 1, nGhost

       i=iGhost(n)
       j=jGhost(n)
       k=kGhost(n)

       flag  = alphaGhost(n)
       flag  = cos(PI/2.0d0*(1.0d0 - flag))

       massflux = massflux +  flag *                                  &
                  ( (-face_u(i,j,k) + face_u(i+1,j,k))*dy(j)*dz(k) +  &
                    (-face_v(i,j,k) + face_v(i,j+1,k))*dx(i)*dz(k) +  &
                    (-face_w(i,j,k) + face_w(i,j,k+1))*dx(i)*dy(j)    &
                  )

       !Calculate the divergence at the ghost cells
       div(i,j,k)    = (dtinv)*                                       &
                      ( ( face_u(i+1,j,k) - face_u(i,j,k) )*dxinv(i)  &
                       +( face_v(i,j+1,k) - face_v(i,j,k) )*dyinv(j)  &
                       +( face_w(i,j,k+1) - face_w(i,j,k) )*dzinv(k)  )

       div(i,j,k)    = div(i,j,k)*flag

    ENDDO

  END SUBROUTINE correct_ghost_massflux
!-------------------------------------------------------------------------------    
  SUBROUTINE correct_ghostcell_velocity(pVar)
!--------------------------------------------
! This subroutine uses the interpolated pressure gradient
! to correct the velocity at the ghost nodes.
!--------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE GCM_arrays
    USE MPI_module

    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2),  INTENT (IN) :: pVar

    INTEGER :: n,iBody
    INTEGER :: info
    INTEGER :: i,j,k

    REAL(KIND=CGREAL) :: dt_fac, flag

    REAL(KIND=CGREAL)    :: pe,pw,pn,ps,pf,pb
    REAL(KIND=CGREAL)    :: pgx,pgy,pgz
!************************************************************

    if(advec_scheme==RUNGE_KUTTA3) then
       dt_fac = 2.0_CGREAL * rk_bet(rk_step) * dt
    else
       dt_fac = dt
    endif    

    DO n = 1, nGhost

       i = iGhost(n) 
       j = jGhost(n)
       k = kGhost(n)
       
       if(mix_GC_form == 1) then
          flag  = alphaGhost(n)
       else
          flag  = 0.0_CGREAL  ! do not update the hybrid "ghost" nodes
       endif
       !flag  = 1.0_CGREAL

       pe = ( fx(i+1)*pVar(i+1,j,k) + (1.0_CGREAL-fx(i+1))*pVar(i,j,k)   )
       pw = ( fx(i)  *pVar(i,j,k)   + (1.0_CGREAL-fx(i))  *pVar(i-1,j,k) ) 

       pn = ( fy(j+1)*pVar(i,j+1,k) + (1.0_CGREAL-fy(j+1))*pVar(i,j,k)   ) 
       ps = ( fy(j)  *pVar(i,j,k)   + (1.0_CGREAL-fy(j))  *pVar(i,j-1,k) )

       pf = ( fz(k+1)*pVar(i,j,k+1) + (1.0_CGREAL-fz(k+1))*pVar(i,j,k)   ) 
       pb = ( fz(k)  *pVar(i,j,k)   + (1.0_CGREAL-fz(k))  *pVar(i,j,k-1) ) 

       pgx= (pe-pw)*dxinv(i)
       pgy= (pn-ps)*dyinv(j)
       pgz= (pf-pb)*dzinv(k)

       u(i,j,k) = u(i,j,k) - dt_fac*pgx*flag
       v(i,j,k) = v(i,j,k) - dt_fac*pgy*flag
       w(i,j,k) = w(i,j,k) - dt_fac*pgz*flag

       !u(i,j,k) = u(i,j,k) - dt_fac*dpdx(i,j,k)*flag
       !v(i,j,k) = v(i,j,k) - dt_fac*dpdy(i,j,k)*flag
       !w(i,j,k) = w(i,j,k) - dt_fac*dpdz(i,j,k)*flag

    ENDDO

    END SUBROUTINE correct_ghostcell_velocity
!----------------------------------------------
! compute mass flux at all boundaries and adjust outflow BC so as to satisfy
! global mass conservation.

   SUBROUTINE GCM_enforce_p_compatibility(pres)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE GCM_arrays
    USE MPI_module

    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2), INTENT(IN) :: pres

    INTEGER             :: i,j,k,n
    REAL(KIND=CGREAL)   :: pgradflux,correction_pgrad,flag
    REAL(KIND=CGREAL)   :: pgxw,pgxe,pgys,pgyn,pgzb,pgzf
    REAL(KIND=CGREAL)   :: flux_i

    IF (pbcx1 == PBC_DIRICHLET .OR. pbcx2 == PBC_DIRICHLET .OR.  & 
        pbcy1 == PBC_DIRICHLET .OR. pbcy2 == PBC_DIRICHLET .OR.  & 
        pbcz1 == PBC_DIRICHLET .OR. pbcz2 == PBC_DIRICHLET ) THEN    
       ! If pressure Dirichlet condition present, then there is no need to
       ! enforce the compatibility condition.
       goto 1000                        ! Do nothing
    END IF

    pgradflux = 0.0_CGREAL

    ! add pressure gradient to face velocities
    DO k = zc_start,zc_end  !1,nz-1
    DO j = yc_start,yc_end
    DO i = 1,nx-1
      pgxw     = (pres(i,j,k)  -pres(i-1,j,k))*dxcinv(i)
      pgxe     = (pres(i+1,j,k)-pres(i,j,k)  )*dxcinv(i+1)

      pgys     = (pres(i,j,k)  -pres(i,j-1,k))*dycinv(j)
      pgyn     = (pres(i,j+1,k)-pres(i,j,k)  )*dycinv(j+1)

      pgzb     = (pres(i,j,k)  -pres(i,j,k-1))*dzcinv(k)
      pgzf     = (pres(i,j,k+1)-pres(i,j,k)  )*dzcinv(k+1)


      flag     = (1.0_CGREAL-REAL(iblank(i,j,k),       KIND=CGREAL))   &
               * (1.0_CGREAL-REAL(ghostCellMark(i,j,k),KIND=CGREAL))

      !if(i==1)    pgxw = pgradx1(j,k)
      !if(i==nx-1) pgxe = pgradx2(j,k)
      !if(j==1)    pgys = pgrady1(i,k)
      !if(j==ny-1) pgyn = pgrady2(i,k)
      !if(k==1)    pgzb = pgradz1(i,j)
      !if(k==nz-1) pgzf = pgradz2(i,j)
      
      pgradflux= pgradflux +  flag *                                 &
                 ( (-pgxw + pgxe)*dy(j)*dz(k) + (-pgys + pgyn)*dx(i)*dz(k) + &
                   (-pgzb + pgzf)*dx(i)*dy(j)      ) 
    ENDDO
    ENDDO
    ENDDO

    IF(mix_GC_form == 1) then
      DO n = 1,nGhost

        i = iGhost(n)
        j = jGhost(n)
        k = kGhost(n)

        pgxw = (pres(i,j,k)  -pres(i-1,j,k))*dxcinv(i)
        pgxe = (pres(i+1,j,k)-pres(i,j,k)  )*dxcinv(i+1)
        
        pgys = (pres(i,j,k)  -pres(i,j-1,k))*dycinv(j)
        pgyn = (pres(i,j+1,k)-pres(i,j,k)  )*dycinv(j+1)

        pgzb = (pres(i,j,k)  -pres(i,j,k-1))*dzcinv(k)
        pgzf = (pres(i,j,k+1)-pres(i,j,k)  )*dzcinv(k+1)

        pgradflux= pgradflux +                                                 &
                   ( (-pgxw + pgxe)*dy(j)*dz(k) + (-pgys + pgyn)*dx(i)*dz(k) + &
                     (-pgzb + pgzf)*dx(i)*dy(j)                                &
                   ) - res_laplacian(n)*dx(i)*dy(j)*dz(k)
      ENDDO
    ENDIF
    flux_i = pgradflux

    !sum up flux_i from all subdomains, then all processors get the result.
    call MPI_ALLREDUCE(flux_i, pgradflux, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, FLOW_COMM, ierr)

    ! adjust BC at outflow to satisfy global mass conservation
    correction_pgrad =-pgradflux/outflow_area  


    !IF ( MOD(ntime,nmonitor) == 0 ) THEN
    !   PRINT*,' SET_BC:pgrad flux       = ',pgradflux
    !   PRINT*,' SET_BC:outflow_area     = ',outflow_area
    !   PRINT*,' SET_BC:correction_pgrad = ',correction_pgrad
    !END IF ! ntime

    Do k=zc_start-1,zc_end+1   !0,nz
    DO j=yc_start-1,yc_end+1
       IF (bcx1 == BC_TYPE_ZERO_GRADIENT) pgradx1(j,k) = pgradx1(j,k) - correction_pgrad
       IF (bcx2 == BC_TYPE_ZERO_GRADIENT) pgradx2(j,k) = pgradx2(j,k) + correction_pgrad
    ENDDO ! j
    ENDDO 

    if (jProc==0)then
        DO k=zc_start-1,zc_end+1   !0,nz
        DO i=0,nx
           IF (bcy1 == BC_TYPE_ZERO_GRADIENT) pgrady1(i,k) = pgrady1(i,k) - correction_pgrad 
        ENDDO ! i
        ENDDO
    endif

    if (jProc==nProcY-1)then
        DO k=zc_start-1,zc_end+1   !0,nz
        DO i=0,nx
           IF (bcy2 == BC_TYPE_ZERO_GRADIENT) pgrady2(i,k) = pgrady2(i,k) + correction_pgrad 
        ENDDO ! i
        ENDDO
    endif

    if(kProc==0) then
       DO j=yc_start-1,yc_end+1
       DO i=0,nx
          IF (bcz1 == BC_TYPE_ZERO_GRADIENT) pgradz1(i,j) = pgradz1(i,j) - correction_pgrad 
       ENDDO ! i
       ENDDO
    endif
    if(kProc==nProcZ-1) then
       DO j=yc_start-1,yc_end+1
       DO i=0,nx
          IF (bcz2 == BC_TYPE_ZERO_GRADIENT) pgradz2(i,j) = pgradz2(i,j) + correction_pgrad 
       ENDDO ! i
       ENDDO
    endif

1000  continue 

   END SUBROUTINE GCM_enforce_p_compatibility
!-------------------------------------------------------------------------
