!-------------------------------------------------
!  SUBROUTINE solve_ad()
!  SUBROUTINE itsolv_ad(var,r)
!  SUBROUTINE itsolv_ad_x(var,r)
!  SUBROUTINE itsolv_ad_y(var,r)
!  SUBROUTINE calc_residual_ad(var,r,resm)
!-------------------------------------------------

!-------------------------------------------------
!-------------------------------------------------
   SUBROUTINE solve_ad()

    USE global_parameters
    USE flow_parameters
    USE flow_parameters
    USE flow_arrays
    USE multiuse_arrays
    USE boundary_arrays
    USE MPI_module

    USE grid_arrays

    IMPLICIT NONE

!... Local Variables

    INTEGER           :: iter,i,j,k,loc(3),locu(3),locv(3),locw(3)
    REAL(KIND=CGREAL) :: res2u,res2v,res2w
    REAL(KIND=CGREAL) :: res2_i,res_ip(2),res_loc(2),res2
    INTEGER           :: ip_loc, iter2

!******************************************************************************

    iter      = 0
    res2      = 1.E10_CGREAL 

    DO WHILE ((iter .LT. itermax_ad) .AND. (res2 .GT. restol_ad))

       ! set the outer boundary condition and the ghost velocity
       ! solve the advection-diffusion equations
       CALL itsolv_ad(u,nlu)
       CALL itsolv_ad(v,nlv)
       CALL itsolv_ad(w,nlw)

       !Update buffer slices for each subdomain
       call send_receive_slices_real_y(u, 0,nx+1,yb1,yb2,zb1,zb2,2)
       call send_receive_slices_real_z(u, 0,nx+1,yb1,yb2,zb1,zb2,2)
       call send_receive_slices_real_y(v, 0,nx+1,yb1,yb2,zb1,zb2,2)
       call send_receive_slices_real_z(v, 0,nx+1,yb1,yb2,zb1,zb2,2)
       call send_receive_slices_real_y(w, 0,nx+1,yb1,yb2,zb1,zb2,2)
       call send_receive_slices_real_z(w, 0,nx+1,yb1,yb2,zb1,zb2,2)

       ! adjust the outer B.C. to enforce the global mass conservation
       CALL GCM_ghostcell_velocity(2,u,v,w)
       !Update buffer slices for each subdomain
       !call send_receive_slices_real(u, 0,nx+1,0,ny+1,zb1,zb2,1)
       !call send_receive_slices_real(v, 0,nx+1,0,ny+1,zb1,zb2,1)
       !call send_receive_slices_real(w, 0,nx+1,0,ny+1,zb1,zb2,1)
       !call write_subdomain()

       CALL face_vel(u, v, w)

       CALL enforce_global_mass_consv()
       CALL set_outer_velocity_bc()
       CALL set_outer_ghost_velocity()

       CALL calc_residual_ad(u,nlu,res2u,locu)
       CALL calc_residual_ad(v,nlv,res2v,locv)
       CALL calc_residual_ad(w,nlw,res2w,locw)
       !write(6,'(a,I5,3E15.5)'),'lProc,res2u,res2v,res2w', &
       !                           lProc,res2u,res2v,res2w

       !Find the max residual from 3 velocity components
       res2_i = abs(res2u)
       loc    = locu
       IF (abs(res2v) > res2_i) THEN
          res2_i = abs(res2v)
          loc    = locv
       ENDIF
       IF (abs(res2w) > res2_i) THEN
          res2_i = abs(res2w)
          loc    = locw
       ENDIF
       !if(lProc==2) res2_i = 0.0001   ! test MPI_MAXLOC
       !write(6,'(a,I5,a,1E15.5,3I5)'),'lProc',lProc,' res2_i =',res2_i,loc(1:3)

       !Pack up residual and processor ID for comparison
       res_ip(1) = res2_i
       res_ip(2) = lProc

       !Find max error among all subdomains and propagate it to all
       !processors.
       CALL MPI_ALLREDUCE(res_ip, res_loc, 1, MPI_2DOUBLE_PRECISION, &
                          MPI_MAXLOC, FLOW_COMM, ierr)

       !if(lProc==PROC_M) then
       !    write(6,'(a,E15.5,F5.0)'),'Max velocity residual and lProc',res_loc(1),res_loc(2)
       !endif
       res2   = res_loc(1)
       ip_loc = res_loc(2)  ! coerce into integer

       if(ip_loc /= PROC_M) then  !send the location to PROC_M
          if(lProc==ip_loc) &
               call MPI_SEND(loc, 3, MPI_INTEGER, PROC_M, 1,FLOW_COMM,istatus,ierr)
          if(lProc==PROC_M) &
               call MPI_RECV(loc, 3, MPI_INTEGER, ip_loc, 1,FLOW_COMM,istatus,ierr)
       endif

       iter = iter + 1

       !call mpi_barrier(FLOW_COMM,ierr)

       IF(lProc==proc_m) THEN
          IF ( iter .EQ. itermax_ad .AND. res2 .GT. restol_ad ) THEN
            write(6,'(a,I5,a)'),'Velocity did not converge in ',itermax_ad,' iterations'
            write(6,'(a,E15.5)'),'Final residual = ',res2
            write(6,'(a,3I5)'),'  Location of MAX residual:',loc(1),loc(2),loc(3)
          ELSE
           IF (MOD(ntime,nmonitor) == 0) THEN
             write(6,'(a,I5,a,PE15.5)'),'iter=',iter,' Velocity max residual:',res2
             write(6,'(a,3I5)'),'  Location of max residual:',loc(1),loc(2),loc(3)
           END IF ! ntime
          ENDIF
       ENDIF
    ENDDO ! do while

   END SUBROUTINE solve_ad
!----------------------------------------------------------
! currently coded as Line SOR with Gauss Siedel as smoother
!----------------------------------------------------------
   SUBROUTINE itsolv_ad(var,r)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_arrays
    USE MPI_module

    IMPLICIT NONE

!... Parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2),  INTENT (IN OUT) ::var
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2),  INTENT (IN)     ::r

    !Ye,test
    INTEGER           :: i,j,k

    ! Line solve in the x-direction
    CALL itsolv_ad_x(var,r)

    ! Line solve in the y-direction    
    CALL itsolv_ad_y(var,r)

    ! Line solve in the z-direction
    CALL itsolv_ad_z(var,r)

   END SUBROUTINE itsolv_ad
!----------------------------------------------

   SUBROUTINE itsolv_ad_x(var,r)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_arrays
    USE GCM_arrays      
    USE flow_arrays
    USE MPI_module

    IMPLICIT NONE

!... Parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2),  INTENT (IN OUT) ::var
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2),  INTENT (IN)     ::r

!... Loop Variables

    INTEGER :: i,j,k

!... Local Variables

    REAL(KIND=CGREAL) :: oned, twod, half, half_dt
    REAL(KIND=CGREAL) :: alpha_p, beta_p, tmp1, tmp2
    REAL(KIND=CGREAL) :: riblank

    INTEGER :: IP,IM,JP,JM,KP,KM

    REAL(KIND=CGREAL) :: hyb1, hyb2
!******************************************************************************

    oned     = 1.0_CGREAL
    twod     = 2.0_CGREAL
    half     = 0.5_CGREAL

    IF(advec_scheme == CRANK_NICOLSON2)        THEN
       half_dt  = 0.5_CGREAL * dt
    ELSE
       half_dt  = 0.0_CGREAL
    ENDIF

! Line solve in the x-direction

  DO k=zc_start,zc_end     !1,nz-1
    kp   = k+1
    km   = k-1    

    DO j=yc_start,yc_end
     JP   = j+1 
     JM   = j-1 

      DO i=1,nx-1

       amx(i) = amx_ad(i)
       apx(i) = apx_ad(i)
       acx(i) =- ( amx(i) + apx(i) )      

       amy(j) = amy_ad(j)
       apy(j) = apy_ad(j)
       acy(j) =- ( amy(j) + apy(j) )     

       amz(k) = amz_ad(k)
       apz(k) = apz_ad(k)
       acz(k) =- ( amz(k) + apz(k) )     

       ! Note that this part only takes effect when the adection terms were
       ! treated with the CN scheme, since otherwise half_dt = 0.
       !--------------------------------------------------------------
       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_u(i,j,k)))
       hyb2 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_u(i+1,j,k)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   = oned - fx(i  );   
       tmp2   = fx(i+1)
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          amx(i) = amx(i) - half_dt * face_u(i  ,j,k)* tmp1 *dxinv(i)
          apx(i) = apx(i) + half_dt * face_u(i+1,j,k)* tmp2 *dxinv(i) 
       ELSE
          amx(i) = amx(i)-half_dt*hyb1*face_u(i,j,k)*dxinv(i)
          apx(i) = apx(i)+half_dt*hyb2*face_u(i+1,j,k)*dxinv(i)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_u(i,j,k)))
       hyb2 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_u(i+1,j,k)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   =         fx(i  )
       tmp2   = (oned - fx(i+1))
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          acx(i) = acx(i) + half_dt * (face_u(i+1,j,k)*tmp2 - face_u(i,j,k)*tmp1)*dxinv(i)
       ELSE
          acx(i) = acx(i)+half_dt*(face_u(i+1,j,k)*hyb2-face_u(i,j,k)*hyb1)*dxinv(i)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_v(i,j,k)))
       hyb2 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_v(i,j+1,k)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   = oned - fy(j  )
       tmp2   =        fy(j+1)
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          amy(j) = amy(j) - half_dt * face_v(i,j  ,k)* tmp1 *dyinv(j)
          apy(j) = apy(j) + half_dt * face_v(i,j+1,k)* tmp2 *dyinv(j)
       ELSE
          amy(j) = amy(j)-half_dt*hyb1*face_v(i,j,k)*dyinv(j)
          apy(j) = apy(j)+half_dt*hyb2*face_v(i,j+1,k)*dyinv(j)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_v(i,j,k)))
       hyb2 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_v(i,j+1,k)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   =         fy(j  )
       tmp2   = (oned - fy(j+1))
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          acy(j) = acy(j) + half_dt * (face_v(i,j+1,k)*tmp2 - face_v(i,j,k)*tmp1)*dyinv(j)
       ELSE
          acy(j) = acy(j)+half_dt*(face_v(i,j+1,k)*hyb2-face_v(i,j,k)*hyb1)*dyinv(j)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_w(i,j,k)))
       hyb2 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_w(i,j,k+1)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   = oned - fz(k  )  
       tmp2   =        fz(k+1)
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          amz(k) = amz(k) - half_dt * face_w(i,j,k  )* tmp1 *dzinv(k)
          apz(k) = apz(k) + half_dt * face_w(i,j,k+1)* tmp2 *dzinv(k)
       ELSE
          amz(k) = amz(k)-half_dt*hyb1*face_w(i,j,k)*dzinv(k)
          apz(k) = apz(k)+half_dt*hyb2*face_w(i,j,k+1)*dzinv(k)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_w(i,j,k)))
       hyb2 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_w(i,j,k+1)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   =         fz(k  )
       tmp2   = (oned - fz(k+1))
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          acz(k) = acz(k) + half_dt * (face_w(i,j,k+1)*tmp2 - face_w(i,j,k)*tmp1)*dzinv(k)
       ELSE
          acz(k) = acz(k)+half_dt*(face_w(i,j,k+1)*hyb2-face_w(i,j,k)*hyb1)*dzinv(k)
       ENDIF
       !--------------------------------------------------------------

       rhs(i) = r(i,j,k) - ( var(i,JM,k)*amy(j) + var(i,JP,k)*apy(j) &
                           + var(i,j,km)*amz(k) + var(i,j,kp)*apz(k))

       riblank   = (1.0_CGREAL-REAL(iblank(i,j,k),KIND=CGREAL))       &
                 * (1.0_CGREAL-REAL(abs(ghostCellMark(i,j,k)),KIND=CGREAL))

       amx(i) = amx(i)*riblank
       apx(i) = apx(i)*riblank
       acx(i) = 1.0_CGREAL + ( acx(i)+acy(j)+acz(k) ) *riblank

       rhs(i) = rhs(i)*riblank + var(i,j,k)*(1.0_CGREAL-riblank) 
       
      ENDDO ! i

      ! set up trivial equations at outer boundaries
      amx(0) = 0.0_CGREAL
      apx(0) = 0.0_CGREAL
      acx(0) = 1.0_CGREAL
      rhs(0) = var(0,j,k)

      amx(nx)= 0.0_CGREAL
      apx(nx)= 0.0_CGREAL
      acx(nx)= 1.0_CGREAL
      rhs(nx)= var(nx,j,k)

      CALL tdma(amx,acx,apx,rhs,dummy,0,nx)

      DO i=1,nx-1
         var(i,j,k) = var(i,j,k) + omega_ad*(dummy(i)-var(i,j,k))
         !var(i,j,k) = dummy(i)
      ENDDO

    ENDDO ! k
  ENDDO ! j
 
  call send_receive_slices_real_y(var, 0,nx+1,yb1,yb2,zb1,zb2,1)
  call send_receive_slices_real_z(var, 0,nx+1,yb1,yb2,zb1,zb2,1)
 
END SUBROUTINE itsolv_ad_x
!----------------------------------------------
!
!----------------------------------------------
   SUBROUTINE itsolv_ad_y(var,r)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_arrays
    USE GCM_arrays
    USE flow_arrays
    USE MPI_module

    IMPLICIT NONE

!... Parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2),  INTENT (IN OUT) ::var
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2),  INTENT (IN)     ::r

!... Loop Variables

    INTEGER :: i,j,k

!... Local Variables

    REAL(KIND=CGREAL) :: oned, twod, half, half_dt
    REAL(KIND=CGREAL) :: riblank
    REAL(KIND=CGREAL) :: alpha_p, beta_p, tmp1, tmp2

    INTEGER :: IP,IM,JP,JM,KP,KM
    INTEGER :: j0,j1

    REAL(KIND=CGREAL) :: hyb1, hyb2
!******************************************************************************

    !-------------------
    oned     = 1.0_CGREAL
    twod     = 2.0_CGREAL
    half     = 0.5_CGREAL

    IF(advec_scheme == CRANK_NICOLSON2) THEN
       half_dt  = 0.5_CGREAL * dt
    ELSE
       half_dt  = 0.0_CGREAL
    ENDIF
    !---------------------

! Line solve in the y-direction
  j0 = yc_start-1
  j1 = yc_end  +1

  DO k=zc_start,zc_end   !1,nz-1

    kP   = k + 1 
    kM   = k - 1 

    DO i=1,nx-1

     IP   = i + 1 
     IM   = i - 1 

      DO j=yc_start,yc_end

       amx(i) = amx_ad(i)
       apx(i) = apx_ad(i)
       acx(i) =- ( amx(i) + apx(i) )      

       amy(j) = amy_ad(j)
       apy(j) = apy_ad(j)
       acy(j) =- ( amy(j) + apy(j) )     

    
       amz(k) = amz_ad(k)
       apz(k) = apz_ad(k)
       acz(k) =- ( amz(k) + apz(k) )     

       ! Note that this part only takes effect when the adection terms were
       ! treated with the CN scheme, since otherwise half_dt = 0.
       !----------------------------------------------------------
       !Ye, 1st-order upwind for implicit term
       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_u(i,j,k)))
       hyb2 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_u(i+1,j,k)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   = oned - fx(i  );   
       tmp2   = fx(i+1)
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          amx(i) = amx(i) - half_dt * face_u(i  ,j,k)* tmp1 *dxinv(i)
          apx(i) = apx(i) + half_dt * face_u(i+1,j,k)* tmp2 *dxinv(i) 
       ELSE
          amx(i) = amx(i)-half_dt*hyb1*face_u(i,j,k)*dxinv(i)
          apx(i) = apx(i)+half_dt*hyb2*face_u(i+1,j,k)*dxinv(i)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_u(i,j,k)))
       hyb2 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_u(i+1,j,k)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   =         fx(i  )
       tmp2   = (oned - fx(i+1))
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          acx(i) = acx(i) + half_dt * (face_u(i+1,j,k)*tmp2 - face_u(i,j,k)*tmp1)*dxinv(i)
       ELSE
          acx(i) = acx(i)+half_dt*(face_u(i+1,j,k)*hyb2-face_u(i,j,k)*hyb1)*dxinv(i)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_v(i,j,k)))
       hyb2 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_v(i,j+1,k)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   = oned - fy(j  )
       tmp2   =        fy(j+1)
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          amy(j) = amy(j) - half_dt * face_v(i,j  ,k)* tmp1 *dyinv(j)
          apy(j) = apy(j) + half_dt * face_v(i,j+1,k)* tmp2 *dyinv(j)
       ELSE
          amy(j) = amy(j)-half_dt*hyb1*face_v(i,j,k)*dyinv(j)
          apy(j) = apy(j)+half_dt*hyb2*face_v(i,j+1,k)*dyinv(j)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_v(i,j,k)))
       hyb2 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_v(i,j+1,k)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   =         fy(j  )
       tmp2   = (oned - fy(j+1))
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          acy(j) = acy(j) + half_dt * (face_v(i,j+1,k)*tmp2 - face_v(i,j,k)*tmp1)*dyinv(j)
       ELSE
          acy(j) = acy(j)+half_dt*(face_v(i,j+1,k)*hyb2-face_v(i,j,k)*hyb1)*dyinv(j)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_w(i,j,k)))
       hyb2 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_w(i,j,k+1)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   = oned - fz(k  )  
       tmp2   =        fz(k+1)
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          amz(k) = amz(k) - half_dt * face_w(i,j,k  )* tmp1 *dzinv(k)
          apz(k) = apz(k) + half_dt * face_w(i,j,k+1)* tmp2 *dzinv(k)
       ELSE
          amz(k) = amz(k)-half_dt*hyb1*face_w(i,j,k)*dzinv(k)
          apz(k) = apz(k)+half_dt*hyb2*face_w(i,j,k+1)*dzinv(k)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_w(i,j,k)))
       hyb2 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_w(i,j,k+1)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   =         fz(k  )
       tmp2   = (oned - fz(k+1))
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          acz(k) = acz(k) + half_dt * (face_w(i,j,k+1)*tmp2 - face_w(i,j,k)*tmp1)*dzinv(k)
       ELSE
          acz(k) = acz(k)+half_dt*(face_w(i,j,k+1)*hyb2-face_w(i,j,k)*hyb1)*dzinv(k)
       ENDIF
       !----------------------------------------------------------

       rhs(j) = r(i,j,k) -( var(IM,j,k)*amx(i)  + var(IP,j,k)*apx(i)  &
                          + var(i,j,km)*amz(k)  + var(i,j,kp)*apz(k) )

       riblank   = (1.0_CGREAL-REAL(iblank(i,j,k),KIND=CGREAL))       &
                 * (1.0_CGREAL-REAL(abs(ghostCellMark(i,j,k)),KIND=CGREAL))

       amy(j) = amy(j)*riblank
       apy(j) = apy(j)*riblank

       acy(j) = 1.0_CGREAL + ( acx(i)+acy(j)+acz(k) ) * riblank
       rhs(j) = rhs(j)*riblank + var(i,j,k)*(1.0_CGREAL-riblank)        
      
      ENDDO ! j

      ! set up trivial equations at outer boundaries
      amy(j0)  = 0.0_CGREAL
      apy(j0)  = 0.0_CGREAL
      acy(j0)  = 1.0_CGREAL
      rhs(j0)  = var(i,yc_start-1,k)

      amy(j1) = 0.0_CGREAL
      apy(j1) = 0.0_CGREAL
      acy(j1) = 1.0_CGREAL
      rhs(j1) = var(i,yc_end+1,k)

      CALL tdma(amy,acy,apy,rhs,dummy,j0,j1)

      DO j=yc_start,yc_end
         var(i,j,k) = var(i,j,k) + omega_ad*(dummy(j)-var(i,j,k))
         !var(i,j,k) = dummy(j)
      ENDDO ! j

    ENDDO ! i
   ENDDO ! k

  call send_receive_slices_real_y(var, 0,nx+1,yb1,yb2,zb1,zb2,1)
  call send_receive_slices_real_z(var, 0,nx+1,yb1,yb2,zb1,zb2,1)

   END SUBROUTINE itsolv_ad_y
!----------------------------------------------------------

!----------------------------------------------------------
   SUBROUTINE itsolv_ad_z(var,r)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_arrays
    USE GCM_arrays
    USE flow_arrays 
    USE MPI_module

    IMPLICIT NONE

!... Parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2),  INTENT (IN OUT) ::var
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2),  INTENT (IN)     ::r

!... Loop Variables

    INTEGER :: i,j,k

!... Local Variables

    REAL(KIND=CGREAL) :: oned, twod, half, half_dt
    REAL(KIND=CGREAL) :: riblank
    REAL(KIND=CGREAL) :: alpha_p, beta_p, tmp1, tmp2

    INTEGER :: n
    INTEGER :: IP,IM,JP,JM,KP,KM
    INTEGER :: k0,k1

    REAL(KIND=CGREAL) :: hyb1, hyb2
!******************************************************************************

    !-------------------
    oned     = 1.0_CGREAL
    twod     = 2.0_CGREAL
    half     = 0.5_CGREAL

    IF( advec_scheme == CRANK_NICOLSON2) THEN
       half_dt  = 0.5_CGREAL * dt
    ELSE
       half_dt  = 0.0_CGREAL
    ENDIF
    !---------------------
   
! Line solve in the z-direction
   k0    = zc_start-1
   k1    = zc_end  +1

   DO j=yc_start,yc_end

    jP   = j + 1 
    jM   = j - 1 

    DO i=1,nx-1

      IP   = i + 1 
      IM   = i - 1 

      DO k=zc_start,zc_end  !1,nz-1   

       amx(i) = amx_ad(i)
       apx(i) = apx_ad(i)
       acx(i) =- ( amx(i) + apx(i) )      

       amy(j) = amy_ad(j)
       apy(j) = apy_ad(j)
       acy(j) =- ( amy(j) + apy(j) )     

    
       amz(k) = amz_ad(k)
       apz(k) = apz_ad(k)
       acz(k) =- ( amz(k) + apz(k) )     

   
       ! Note that this part only takes effect when the adection terms were
       ! treated with the CN scheme, since otherwise half_dt = 0.
       !------------------------------------------------------------
       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_u(i,j,k)))
       hyb2 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_u(i+1,j,k)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   = oned - fx(i  );   
       tmp2   = fx(i+1)
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          amx(i) = amx(i) - half_dt * face_u(i  ,j,k)* tmp1 *dxinv(i)
          apx(i) = apx(i) + half_dt * face_u(i+1,j,k)* tmp2 *dxinv(i) 
       ELSE
          amx(i) = amx(i)-half_dt*hyb1*face_u(i,j,k)*dxinv(i)
          apx(i) = apx(i)+half_dt*hyb2*face_u(i+1,j,k)*dxinv(i)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_u(i,j,k)))
       hyb2 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_u(i+1,j,k)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   =         fx(i  )
       tmp2   = (oned - fx(i+1))
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          acx(i) = acx(i) + half_dt * (face_u(i+1,j,k)*tmp2 - face_u(i,j,k)*tmp1)*dxinv(i)
       ELSE
          acx(i) = acx(i)+half_dt*(face_u(i+1,j,k)*hyb2-face_u(i,j,k)*hyb1)*dxinv(i)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_v(i,j,k)))
       hyb2 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_v(i,j+1,k)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   = oned - fy(j  )
       tmp2   =        fy(j+1)
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          amy(j) = amy(j) - half_dt * face_v(i,j  ,k)* tmp1 *dyinv(j)
          apy(j) = apy(j) + half_dt * face_v(i,j+1,k)* tmp2 *dyinv(j)
       ELSE
          amy(j) = amy(j)-half_dt*hyb1*face_v(i,j,k)*dyinv(j)
          apy(j) = apy(j)+half_dt*hyb2*face_v(i,j+1,k)*dyinv(j)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_v(i,j,k)))
       hyb2 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_v(i,j+1,k)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   =         fy(j  )
       tmp2   = (oned - fy(j+1))
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          acy(j) = acy(j) + half_dt * (face_v(i,j+1,k)*tmp2 - face_v(i,j,k)*tmp1)*dyinv(j)
       ELSE
          acy(j) = acy(j)+half_dt*(face_v(i,j+1,k)*hyb2-face_v(i,j,k)*hyb1)*dyinv(j)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_w(i,j,k)))
       hyb2 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_w(i,j,k+1)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   = oned - fz(k  )  
       tmp2   =        fz(k+1)
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          amz(k) = amz(k) - half_dt * face_w(i,j,k  )* tmp1 *dzinv(k)
          apz(k) = apz(k) + half_dt * face_w(i,j,k+1)* tmp2 *dzinv(k)
       ELSE
          amz(k) = amz(k)-half_dt*hyb1*face_w(i,j,k)*dzinv(k)
          apz(k) = apz(k)+half_dt*hyb2*face_w(i,j,k+1)*dzinv(k)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_w(i,j,k)))
       hyb2 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_w(i,j,k+1)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   =         fz(k  )
       tmp2   = (oned - fz(k+1))
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          acz(k) = acz(k) + half_dt * (face_w(i,j,k+1)*tmp2 - face_w(i,j,k)*tmp1)*dzinv(k)
       ELSE
          acz(k) = acz(k)+half_dt*(face_w(i,j,k+1)*hyb2-face_w(i,j,k)*hyb1)*dzinv(k)
       ENDIF
       !------------------------------------------------------------

       rhs(k) = r(i,j,k) -( var(IM,j,k)*amx(i)  + var(IP,j,k)*apx(i)  &
                          + var(i,jm,k)*amy(j)  + var(i,jp,k)*apy(j) )

       if(abs(ghostCellMark(i,j,k)) == 1 ) then
          n = ghostCellIndex(i,j,k)
          ghostNScoeff(1,n) = 1.0_CGREAL + ( acx(i)+acy(j) + acz(k))
          ghostNScoeff(2,n) = amx(i)
          ghostNScoeff(3,n) = apx(i)
          ghostNScoeff(4,n) = amy(j)
          ghostNScoeff(5,n) = apy(j)
          ghostNScoeff(6,n) = amz(k)
          ghostNScoeff(7,n) = apz(k)
       endif 

       riblank   = (1.0_CGREAL-REAL(iblank(i,j,k),KIND=CGREAL))       &
                 * (1.0_CGREAL-REAL(abs(ghostCellMark(i,j,k)),KIND=CGREAL))

       amz(k) = amz(k)*riblank
       apz(k) = apz(k)*riblank

       acz(k) = 1.0_CGREAL + ( acx(i)+acy(j)+acz(k) ) * riblank
       rhs(k) = rhs(k)*riblank + var(i,j,k)*(1.0_CGREAL-riblank)       
      
      ENDDO ! k

      !! set up trivial equations at outer boundaries
      !amz(0)  = 0.0_CGREAL
      !apz(0)  = 0.0_CGREAL
      !acz(0)  = 1.0_CGREAL
      !rhs(0)  = var(i,j,0)

      !amz(nz) = 0.0_CGREAL
      !apz(nz) = 0.0_CGREAL
      !acz(nz) = 1.0_CGREAL
      !rhs(nz) = var(i,j,nz)

      ! set up trivial equations at buffer slices
      amz(k0)  = 0.0_CGREAL
      apz(k0)  = 0.0_CGREAL
      acz(k0)  = 1.0_CGREAL
      rhs(k0)  = var(i,j,zc_start-1)

      amz(k1) = 0.0_CGREAL
      apz(k1) = 0.0_CGREAL
      acz(k1) = 1.0_CGREAL
      rhs(k1) = var(i,j,zc_end+1)

      ! Note that these arrays have actual size of 0:nz+1, but only 
      ! zc_start-1:zc_end+1 are being used. It is very important that
      ! the arrays in TDMA have the same starting index, i.e., 0, even
      ! though the ending index is may be shorter.
      CALL tdma(amz,acz,apz,rhs,dummy,k0,k1)

      DO k=zc_start,zc_end  !1,nz-1
         var(i,j,k) = var(i,j,k) + omega_ad*(dummy(k)-var(i,j,k))
         !var(i,j,k) = dummy(k)
      ENDDO ! j

    ENDDO ! i
    ENDDO ! j

  call send_receive_slices_real_y(var, 0,nx+1,yb1,yb2,zb1,zb2,1)
  call send_receive_slices_real_z(var, 0,nx+1,yb1,yb2,zb1,zb2,1)

   END SUBROUTINE itsolv_ad_z
!--------------------------------

!----------------------------------------------------------
   SUBROUTINE calc_residual_ad(var,r,resm,loc)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE solver_arrays
    USE GCM_arrays
    USE flow_arrays
    USE multiuse_arrays
    USE MPI_module

    IMPLICIT NONE

!... Parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2), INTENT (IN)  ::var,r
    REAL(KIND=CGREAL),                            INTENT (OUT) ::resm
    INTEGER           :: loc(3)

!... Loop Variables

    INTEGER            :: i,j,k

!... Local Variables
    
    INTEGER            :: IP,IM,JP,JM,KP,KM

    REAL(KIND=CGREAL)  :: res, res2
    REAL(KIND=CGREAL)  :: oned, twod, half, half_dt, tmp1, tmp2
    INTEGER            :: iLoc, jLoc,kLoc                    

    REAL(KIND=CGREAL) :: hyb1, hyb2
!******************************************************************************
    oned     = 1.0_CGREAL
    twod     = 2.0_CGREAL
    half     = 0.5_CGREAL

    IF(advec_scheme == CRANK_NICOLSON2) THEN
       half_dt  = 0.5_CGREAL * dt
    ELSE
       half_dt  = 0.0_CGREAL
    ENDIF

    resm = 0.0_CGREAL
    res2 = 0.0_cgreal

    DO k=zc_start,zc_end    !1,nz-1  
    DO j=yc_start,yc_end
    DO i=1,nx-1

       IP   = i + 1 
       IM   = i - 1 
 
       JP   = j + 1 
       JM   = j - 1 
       kP   = k + 1 
       kM   = k - 1 
 
       amx(i) = amx_ad(i)
       apx(i) = apx_ad(i)
       acx(i) = - ( amx(i) + apx(i) )

       amy(j) = amy_ad(j)
       apy(j) = apy_ad(j)
       acy(j) = - ( amy(j) + apy(j) )
  
       amz(k) = amz_ad(k)
       apz(k) = apz_ad(k)
       acz(k) = - ( amz(k) + apz(k) )

       ! Note that this part only takes effect when the adection terms were
       ! treated with the CN scheme, since otherwise half_dt = 0.
       !-------------------------------------------------------
       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_u(i,j,k)))
       hyb2 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_u(i+1,j,k)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   = oned - fx(i  );   
       tmp2   = fx(i+1)
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          amx(i) = amx(i) - half_dt * face_u(i  ,j,k)* tmp1 *dxinv(i)
          apx(i) = apx(i) + half_dt * face_u(i+1,j,k)* tmp2 *dxinv(i) 
       ELSE
          amx(i) = amx(i)-half_dt*hyb1*face_u(i,j,k)*dxinv(i)
          apx(i) = apx(i)+half_dt*hyb2*face_u(i+1,j,k)*dxinv(i)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_u(i,j,k)))
       hyb2 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_u(i+1,j,k)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   =         fx(i  )
       tmp2   = (oned - fx(i+1))
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          acx(i) = acx(i) + half_dt * (face_u(i+1,j,k)*tmp2 - face_u(i,j,k)*tmp1)*dxinv(i)
       ELSE
          acx(i) = acx(i)+half_dt*(face_u(i+1,j,k)*hyb2-face_u(i,j,k)*hyb1)*dxinv(i)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_v(i,j,k)))
       hyb2 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_v(i,j+1,k)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   = oned - fy(j  )
       tmp2   =        fy(j+1)
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          amy(j) = amy(j) - half_dt * face_v(i,j  ,k)* tmp1 *dyinv(j)
          apy(j) = apy(j) + half_dt * face_v(i,j+1,k)* tmp2 *dyinv(j)
       ELSE
          amy(j) = amy(j)-half_dt*hyb1*face_v(i,j,k)*dyinv(j)
          apy(j) = apy(j)+half_dt*hyb2*face_v(i,j+1,k)*dyinv(j)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_v(i,j,k)))
       hyb2 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_v(i,j+1,k)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   =         fy(j  )
       tmp2   = (oned - fy(j+1))
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          acy(j) = acy(j) + half_dt * (face_v(i,j+1,k)*tmp2 - face_v(i,j,k)*tmp1)*dyinv(j)
       ELSE
          acy(j) = acy(j)+half_dt*(face_v(i,j+1,k)*hyb2-face_v(i,j,k)*hyb1)*dyinv(j)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_w(i,j,k)))
       hyb2 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_w(i,j,k+1)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   = oned - fz(k  )  
       tmp2   =        fz(k+1)
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          amz(k) = amz(k) - half_dt * face_w(i,j,k  )* tmp1 *dzinv(k)
          apz(k) = apz(k) + half_dt * face_w(i,j,k+1)* tmp2 *dzinv(k)
       ELSE
          amz(k) = amz(k)-half_dt*hyb1*face_w(i,j,k)*dzinv(k)
          apz(k) = apz(k)+half_dt*hyb2*face_w(i,j,k+1)*dzinv(k)
       ENDIF

       !Ye, 1st-order upwind scheme for implicit term
       hyb1 = 0.5d0*(1.0d0-SIGN(1.0_CGREAL,face_w(i,j,k)))
       hyb2 = 0.5d0*(1.0d0+SIGN(1.0_CGREAL,face_w(i,j,k+1)))
       !Ye, 2nd-order central scheme for implicit term
       tmp1   =         fz(k  )
       tmp2   = (oned - fz(k+1))
       IF ( (i.ge.25).and.(i.le.375) ) THEN
          acz(k) = acz(k) + half_dt * (face_w(i,j,k+1)*tmp2 - face_w(i,j,k)*tmp1)*dzinv(k)
       ELSE
          acz(k) = acz(k)+half_dt*(face_w(i,j,k+1)*hyb2-face_w(i,j,k)*hyb1)*dzinv(k)
       ENDIF
       !-------------------------------------------------------


       res    = r(i,j,k)   - var(i,j,k) * (1.0_CGREAL + acx(i) + acy(j)+acz(k)  )          &
                         - var(IM,j,k)*amx(i)                                     &
                         - var(IP,j,k)*apx(i)                                     &
                         - var(i,JM,k)*amy(j)                                     &
                         - var(i,JP,k)*apy(j)                                     &
                         - var(i,J,km)*amz(k)                                     &
                         - var(i,J,kP)*apz(k)             
      
       !res    = res*(1.0_CGREAL-REAL(iblank(i,j,k),KIND=CGREAL))
       res    = res*(1.0_CGREAL-REAL(iblank(i,j,k),KIND=CGREAL))       &
                   *(1.0_CGREAL-REAL(ghostCellMark(i,j,k),KIND=CGREAL))

       IF (ABS(res) > resm ) THEN
          resm = ABS(res)
          loc(1) = i
          loc(2) = j
          loc(3) = k
       ENDIF
       !if(ghostCellMark(i,j,k)==1) write(*,'(4I5,3X,F10.5)')lProc,i,j,k,res
       res2 = res2+res**2   ! L2-norm
    ENDDO
    ENDDO
    ENDDO

    !!Find the L-2 norm
    !write(*,'(a,I6,a,1PE15.5)')'lProc:', lProc, ' L-2 residual: ', res2
    !call MPI_ALLREDUCE(res2, resm, 1, MPI_DOUBLE_PRECISION, &
    !                   MPI_SUM, flow_comm, ierr)
    !resm = sqrt(resm)/dfloat(nx*ny*nz)

  END SUBROUTINE calc_residual_ad  
!----------------------------------------------------------
  SUBROUTINE vel_adjust2D() 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays

    IMPLICIT NONE

    INTEGER              :: i,j,k

! Set Velocity field for 2D calculations

! copy k=1 plane to other planes
   
    DO k = 2,nz-1
    DO j = 1,ny-1
    DO i = 1,nx-1
      u(i,j,k) = u(i,j,1)
      v(i,j,k) = v(i,j,1)
      face_u(i,j,k) = face_u(i,j,1)
      face_v(i,j,k) = face_v(i,j,1)
      face_u(i+1,j,k) = face_u(i+1,j,1)
      face_v(i,j+1,k) = face_v(i,j+1,1)
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k
    
! zero w-component
   
    DO k = 1,nz-1
    DO j = 1,ny-1
    DO i = 1,nx-1
      w(i,j,k) = 0.0_CGREAL
      face_w(i,j,k) = 0.0_CGREAL
      face_w(i,j,k+1) = 0.0_CGREAL
    ENDDO ! i
    ENDDO ! j
    ENDDO ! k

  END SUBROUTINE  vel_adjust2D
!-------------------------------------------------------------------------------
