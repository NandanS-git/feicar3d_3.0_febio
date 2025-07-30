!-------------------------------------------------
! Subroutines :  face_vel() ;  correct_vel()
!-------------------------------------------------

!-------------------------------------------------------------------------------
   SUBROUTINE face_vel(u_vel,v_vel,w_vel)
!  compute face velocities

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE GCM_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE MPI_module

    IMPLICIT NONE

    INTEGER :: i,j,k,iBody,n
    REAL(KIND=CGREAL)    :: pe,pw,pn,ps,pf,pb
    REAL(KIND=CGREAL)    :: pgx,pgy,pgz
    REAL(KIND=CGREAL)    :: pgxe,pgxw,pgyn,pgys,pgzf,pgzb
    REAL(KIND=CGREAL)    :: ghostY, ghostN
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2), INTENT(IN) :: u_vel
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2), INTENT(IN) :: v_vel
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2), INTENT(IN) :: w_vel
    INTEGER :: k1,k2,j1,j2
 
    !face_u
    DO k = zc_start,zc_end  !1,nz-1
    DO j = yc_start,yc_end
    DO i = 2,nx-1  ! start from the 2nd grid point

       face_u(i,j,k)   = (              fx(i) *u_vel(i,j,k)           &
                          + (1.0_CGREAL-fx(i))*u_vel(i-1,j,k)   )
    ENDDO
    ENDDO
    ENDDO

    !face_v
    if (jProc==0)then
       j1 = y_start+1
    else
       j1 = y_start
    endif
    if (jProc==nProcY-1)then
       j2 = y_end-1
    else
       j2 = y_end
    endif

    DO j = j1,j2
    DO k = zc_start,zc_end
    DO i = 1,nx-1
       face_v(i,j,k)   = (              fy(j) *v_vel(i,j,k)           &
                          + (1.0_CGREAL-fy(j))*v_vel(i,j-1,k)   )
    ENDDO
    ENDDO
    ENDDO

    !face_w
    if(kProc==0) then
       k1 = z_start+1
    else
       k1 = z_start
    endif
    if(kProc==nProcZ-1) then
       k2 = z_end-1   
    else
       k2 = z_end
    endif

    DO k = k1,k2
    DO j = yc_start,yc_end
    DO i = 1,nx-1
       face_w(i,j,k)   = (              fz(k) *w_vel(i,j,k)           &
                          + (1.0_CGREAL-fz(k))*w_vel(i,j,k-1)   )
    ENDDO
    ENDDO
    ENDDO

    !Apply the outer boundary conditions for face_u
    DO k = zc_start,zc_end   !1,nz-1
    DO j = yc_start,yc_end
       face_u(1, j,k) = bcxu(0,j,k)  
       face_u(nx,j,k) = bcxu(1,j,k)
    ENDDO
    ENDDO

    !Apply the outer boundary conditions for face_v....
    if (jProc==0)then
       DO k = zc_start,zc_end   !1,nz-1
       DO i = 1,nx-1
          face_v(i,1 ,k) = bcyv(i,0,k)
       ENDDO
       ENDDO
    endif
    if (jProc==nProcZ-1)then
       DO k = zc_start,zc_end   !1,nz-1
       DO i = 1,nx-1
          face_v(i,ny,k) = bcyv(i,1,k)
       ENDDO
       ENDDO
    endif

    !Apply the outer boundary conditions for face_w
    !DO j = 1,ny-1
    !DO i = 1,nx-1
    !   face_w(i,j,1 ) = bczw(i,j,0)
    !   face_w(i,j,nz) = bczw(i,j,1)
    !ENDDO
    !ENDDO
    if(kProc==0) then
       DO j = yc_start,yc_end
       DO i = 1,nx-1
          face_w(i,j,1 ) = bczw(i,j,0)
       ENDDO
       ENDDO
    endif
    if(kProc==nProcZ-1) then
       DO j = yc_start,yc_end
       DO i = 1,nx-1
          face_w(i,j,nz) = bczw(i,j,1)
       ENDDO
       ENDDO
    endif

   END SUBROUTINE face_vel   
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
   SUBROUTINE correct_vel()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE MPI_module

    IMPLICIT NONE
    
    INTEGER :: i,j,k,nmax
    REAL(KIND=CGREAL)    :: pe,pw,pn,ps,pf,pb
    REAL(KIND=CGREAL)    :: pgx,pgy,pgz
    REAL(KIND=CGREAL)    :: pgxe,pgxw,pgyn,pgys,pgzf,pgzb
    REAL(KIND=CGREAL)    :: flag, dt_fac

! correct nodal velocities
! need to account for boundaries...

!    if(advec_scheme==RUNGE_KUTTA3) then
!       dt_fac = 2.0_CGREAL * rk_bet(rk_step) * dt
!    else
       dt_fac = dt
!    endif    

    DO k = zc_start,zc_end   !1,nz-1
    DO j = yc_start,yc_end
    DO i = 1,nx-1

      pe = ( fx(i+1)*pPrime(i+1,j,k) + (1.0_CGREAL-fx(i+1))*pPrime(i,j,k)   )

      pw = ( fx(i)  *pPrime(i,j,k)   + (1.0_CGREAL-fx(i))  *pPrime(i-1,j,k) ) 

      pn = ( fy(j+1)*pPrime(i,j+1,k) + (1.0_CGREAL-fy(j+1))*pPrime(i,j,k)   ) 

      ps = ( fy(j)  *pPrime(i,j,k)   + (1.0_CGREAL-fy(j))  *pPrime(i,j-1,k) )

      pf = ( fz(k+1)*pPrime(i,j,k+1) + (1.0_CGREAL-fz(k+1))*pPrime(i,j,k)   ) 

      pb = ( fz(k)  *pPrime(i,j,k)   + (1.0_CGREAL-fz(k))  *pPrime(i,j,k-1) ) 

      pgx= (pe-pw)*dxinv(i)
      pgy= (pn-ps)*dyinv(j)
      pgz= (pf-pb)*dzinv(k)

      ! exclude the solid cells, but also update the ghost cells
      flag   = (1.0_CGREAL-REAL(iblank       (i,j,k),KIND=CGREAL))  !&
             !* (1.0_CGREAL-REAL(ghostCellMark(i,j,k),KIND=CGREAL))

      u(i,j,k) = u(i,j,k) - dt_fac*pgx*flag
      v(i,j,k) = v(i,j,k) - dt_fac*pgy*flag
      w(i,j,k) = w(i,j,k) - dt_fac*pgz*flag
    ENDDO
    ENDDO
    ENDDO

    !Update buffer slices for each subdomain
    !call send_receive_slices_real(u, 0,nx+1,0,ny+1,zb1,zb2,1)
    !call send_receive_slices_real(v, 0,nx+1,0,ny+1,zb1,zb2,1)
    !call send_receive_slices_real(w, 0,nx+1,0,ny+1,zb1,zb2,1)

!    if(MOD(ntime,nmonitor) == 0) then
!    DO k = zb1,zb2   !1,nz-1
!    DO j = 0,ny+1
!    DO i = 0,nx+1
!      if(pPrime(i,j,k) .lt. -1.5) then
!        write(*,'(5I4,7F12.5)')i,j,k,iblank(i,j,k),ghostCellMark(i,j,k),pPrime(i,j,k),pgx,pgy,pgz,u(i,j,k),v(i,j,k),w(i,j,k)
!      endif
!    ENDDO
!    ENDDO
!    ENDDO
!    endif

   END SUBROUTINE correct_vel
!-------------------------------------------------------------------------------
   SUBROUTINE update_pressure()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE gcm_arrays
    USE MPI_module
    USE boundary_arrays

    IMPLICIT NONE
    INTEGER :: i,j,k,j1,k1

    !Set pressure at trivious nodes (pure solid nodes) to zero
    DO k = zc_start-2,zc_end+2 
       !Ye, kProc=0
       if (k.eq.-1) cycle
    DO j = yc_start-2,yc_end+2
       !Ye, jProc=0
       if (j.eq.-1) cycle
    DO i = 0,nx
       if((iblank(i,j,k).eq.0).or.(dead_cell(i,j,k).eq.-1)) cycle
       pPrime(i,j,k) = 0.0d0          
    ENDDO
    ENDDO
    ENDDO

    if (jProc.eq.0)then
       j1 = yc_start-1
    else
       j1 = yc_start-2
    endif
    if (kProc.eq.0)then
       k1 = zc_start-1
    else
       k1 = zc_start-2
    endif

    ! Also zero out the domain edge points that are trivial
    if ((jProc.eq.0).and.(kProc.eq.0))then
       pPrime(0:nx, 0,  k1        :zc_start-1) = 0.0d0!01
       pPrime(0:nx, 0,  zc_end+1  :zc_end+2  ) = 0.0d0!01
    endif

    if ((jProc.eq.nProcY-1).and.(kProc.eq.0))then
       pPrime(0:nx, ny, k1        :zc_start-1) = 0.0d0!02
       pPrime(0:nx, ny, zc_end+1  :zc_end+2  ) = 0.0d0!02
    endif

    if ((jProc.eq.0).and.(kProc.eq.nProcZ-1))then
       pPrime(0:nx, 0,  k1        :zc_start-1  ) = 0.0d0!03
       pPrime(0:nx, 0,  zc_end+1  :zc_end+2    ) = 0.0d0!03      
    endif

    if ((jProc.eq.nProcY-1).and.(kProc.eq.nProcZ-1))then
       pPrime(0:nx, ny, k1        :zc_start-1   ) = 0.0d0!04
       pPrime(0:nx, ny, zc_end+1  :zc_end+2     ) = 0.0d0!04
    endif

    if (kProc.eq.0)then
       pPrime( 0, j1:yc_end+2,  zc_start-1:zc_start-1) = 0.0d0!05
       pPrime( 0, j1:yc_end+2,  zc_end+1  :zc_end+2  ) = 0.0d0!05

       pPrime(nx, j1:yc_end+2,  zc_start-1:zc_start-1) = 0.0d0!06
       pPrime(nx, j1:yc_end+2,  zc_end+1  :zc_end+2  ) = 0.0d0!06
    endif

    if (kProc.eq.nProcZ-1)then
       pPrime( 0, j1:yc_end+2,  zc_start-2:zc_start-1) = 0.0d0!07
       pPrime( 0, j1:yc_end+2,  zc_end+1  :zc_end+2  ) = 0.0d0!07

       pPrime(nx, j1:yc_end+2,  zc_start-2:zc_start-1) = 0.0d0!08
       pPrime(nx, j1:yc_end+2,  zc_end+1  :zc_end+2  ) = 0.0d0!08
    endif

    if (jProc.eq.0)then
       pPrime( 0,   0,    k1:zc_end+2  ) = 0.0d0!09
       pPrime(nx,   0,    k1:zc_end+2  ) = 0.0d0!10
    endif

    if (jProc.eq.nProcY-1)then
       pPrime( 0,   ny,   k1:zc_end+2  ) = 0.0d0!11
       pPrime(nx,   ny,   k1:zc_end+2  ) = 0.0d0!12
    endif

    ! Copy pPrime to p
    p = pPrime

  END SUBROUTINE update_pressure
!-------------------------------------------------------------------------------
   SUBROUTINE calc_pres_grad()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE MPI_module

    IMPLICIT NONE
    
    INTEGER :: i,j,k,nmax
    REAL(KIND=CGREAL)    :: pe,pw,pn,ps,pf,pb
    REAL(KIND=CGREAL)    :: pgx,pgy,pgz
    REAL(KIND=CGREAL)    :: pgxe,pgxw,pgyn,pgys,pgzf,pgzb
    REAL(KIND=CGREAL)    :: flag, dt_fac

! calculate pressure gradient
! need to account for boundaries...   

    DO k = zc_start,zc_end  !1,nz-1
    DO j = 1,ny-1
    DO i = 1,nx-1

      pe = ( fx(i+1)*pPrime(i+1,j,k) + (1.0_CGREAL-fx(i+1))*pPrime(i,j,k)   )

      pw = ( fx(i)  *pPrime(i,j,k)   + (1.0_CGREAL-fx(i))  *pPrime(i-1,j,k) ) 

      pn = ( fy(j+1)*pPrime(i,j+1,k) + (1.0_CGREAL-fy(j+1))*pPrime(i,j,k)   ) 

      ps = ( fy(j)  *pPrime(i,j,k)   + (1.0_CGREAL-fy(j))  *pPrime(i,j-1,k) )

      pf = ( fz(k+1)*pPrime(i,j,k+1) + (1.0_CGREAL-fz(k+1))*pPrime(i,j,k)   ) 

      pb = ( fz(k)  *pPrime(i,j,k)   + (1.0_CGREAL-fz(k))  *pPrime(i,j,k-1) ) 

      ! exclude the solid cells, but also update the ghost cells
      flag   = (1.0_CGREAL-REAL(iblank       (i,j,k),KIND=CGREAL))  !&
             !* (1.0_CGREAL-REAL(ghostCellMark(i,j,k),KIND=CGREAL))

      dpdx(i,j,k) = (pe-pw)*dxinv(i)*flag
      dpdy(i,j,k) = (pn-ps)*dyinv(j)*flag
      dpdz(i,j,k) = (pf-pb)*dzinv(k)*flag

    ENDDO
    ENDDO
    ENDDO


   END SUBROUTINE calc_pres_grad
!-------------------------------------------------------------------------------
   SUBROUTINE calc_pres_gradx()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE multiuse_arrays
    USE solver_arrays

    IMPLICIT NONE
    
    INTEGER :: i,j,k,nmax
    REAL(KIND=CGREAL)    :: pe,pw,pn,ps,pf,pb
    REAL(KIND=CGREAL)    :: pgx,pgy,pgz
    REAL(KIND=CGREAL)    :: pgxe,pgxw,pgyn,pgys,pgzf,pgzb
    REAL(KIND=CGREAL)    :: flag, dt_fac

! calculate pressure gradient
! need to account for boundaries...   

    DO k = 1,nz-1
    DO j = 1,ny-1
    DO i = 1,nx-1

      pe = ( fx(i+1)*pPrime(i+1,j,k) + (1.0_CGREAL-fx(i+1))*pPrime(i,j,k)   )

      pw = ( fx(i)  *pPrime(i,j,k)   + (1.0_CGREAL-fx(i))  *pPrime(i-1,j,k) ) 

      pn = ( fy(j+1)*pPrime(i,j+1,k) + (1.0_CGREAL-fy(j+1))*pPrime(i,j,k)   ) 

      ps = ( fy(j)  *pPrime(i,j,k)   + (1.0_CGREAL-fy(j))  *pPrime(i,j-1,k) )

      pf = ( fz(k+1)*pPrime(i,j,k+1) + (1.0_CGREAL-fz(k+1))*pPrime(i,j,k)   ) 

      pb = ( fz(k)  *pPrime(i,j,k)   + (1.0_CGREAL-fz(k))  *pPrime(i,j,k-1) ) 

      ! exclude the solid cells, but also update the ghost cells
      flag   = (1.0_CGREAL-REAL(iblank       (i,j,k),KIND=CGREAL))  !&
             !* (1.0_CGREAL-REAL(ghostCellMark(i,j,k),KIND=CGREAL))

      dpdx(i,j,k) = (pe-pw)*dxinv(i)*flag
!      dpdy(i,j,k) = (pn-ps)*dyinv(j)*flag
!      dpdz(i,j,k) = (pf-pb)*dzinv(k)*flag

    ENDDO
    ENDDO
    ENDDO


   END SUBROUTINE calc_pres_gradx
!-------------------------------------------------------------------------------
   SUBROUTINE calc_pres_grady()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE multiuse_arrays
    USE solver_arrays

    IMPLICIT NONE
    
    INTEGER :: i,j,k,nmax
    REAL(KIND=CGREAL)    :: pe,pw,pn,ps,pf,pb
    REAL(KIND=CGREAL)    :: pgx,pgy,pgz
    REAL(KIND=CGREAL)    :: pgxe,pgxw,pgyn,pgys,pgzf,pgzb
    REAL(KIND=CGREAL)    :: flag, dt_fac

! calculate pressure gradient
! need to account for boundaries...   

    DO k = 1,nz-1
    DO j = 1,ny-1
    DO i = 1,nx-1

      pe = ( fx(i+1)*pPrime(i+1,j,k) + (1.0_CGREAL-fx(i+1))*pPrime(i,j,k)   )

      pw = ( fx(i)  *pPrime(i,j,k)   + (1.0_CGREAL-fx(i))  *pPrime(i-1,j,k) ) 

      pn = ( fy(j+1)*pPrime(i,j+1,k) + (1.0_CGREAL-fy(j+1))*pPrime(i,j,k)   ) 

      ps = ( fy(j)  *pPrime(i,j,k)   + (1.0_CGREAL-fy(j))  *pPrime(i,j-1,k) )

      pf = ( fz(k+1)*pPrime(i,j,k+1) + (1.0_CGREAL-fz(k+1))*pPrime(i,j,k)   ) 

      pb = ( fz(k)  *pPrime(i,j,k)   + (1.0_CGREAL-fz(k))  *pPrime(i,j,k-1) ) 

      ! exclude the solid cells, but also update the ghost cells
      flag   = (1.0_CGREAL-REAL(iblank       (i,j,k),KIND=CGREAL))  !&
             !* (1.0_CGREAL-REAL(ghostCellMark(i,j,k),KIND=CGREAL))

!      dpdx(i,j,k) = (pe-pw)*dxinv(i)*flag
      dpdy(i,j,k) = (pn-ps)*dyinv(j)*flag
!      dpdz(i,j,k) = (pf-pb)*dzinv(k)*flag

    ENDDO
    ENDDO
    ENDDO


   END SUBROUTINE calc_pres_grady
!-------------------------------------------------------------------------------
   SUBROUTINE calc_pres_gradz()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE multiuse_arrays
    USE solver_arrays

    IMPLICIT NONE
    
    INTEGER :: i,j,k,nmax
    REAL(KIND=CGREAL)    :: pe,pw,pn,ps,pf,pb
    REAL(KIND=CGREAL)    :: pgx,pgy,pgz
    REAL(KIND=CGREAL)    :: pgxe,pgxw,pgyn,pgys,pgzf,pgzb
    REAL(KIND=CGREAL)    :: flag, dt_fac

! calculate pressure gradient
! need to account for boundaries...   

    DO k = 1,nz-1
    DO j = 1,ny-1
    DO i = 1,nx-1

      pe = ( fx(i+1)*pPrime(i+1,j,k) + (1.0_CGREAL-fx(i+1))*pPrime(i,j,k)   )

      pw = ( fx(i)  *pPrime(i,j,k)   + (1.0_CGREAL-fx(i))  *pPrime(i-1,j,k) ) 

      pn = ( fy(j+1)*pPrime(i,j+1,k) + (1.0_CGREAL-fy(j+1))*pPrime(i,j,k)   ) 

      ps = ( fy(j)  *pPrime(i,j,k)   + (1.0_CGREAL-fy(j))  *pPrime(i,j-1,k) )

      pf = ( fz(k+1)*pPrime(i,j,k+1) + (1.0_CGREAL-fz(k+1))*pPrime(i,j,k)   ) 

      pb = ( fz(k)  *pPrime(i,j,k)   + (1.0_CGREAL-fz(k))  *pPrime(i,j,k-1) ) 

      ! exclude the solid cells, but also update the ghost cells
      flag   = (1.0_CGREAL-REAL(iblank       (i,j,k),KIND=CGREAL))  !&
             !* (1.0_CGREAL-REAL(ghostCellMark(i,j,k),KIND=CGREAL))

!      dpdx(i,j,k) = (pe-pw)*dxinv(i)*flag
!      dpdy(i,j,k) = (pn-ps)*dyinv(j)*flag
      dpdz(i,j,k) = (pf-pb)*dzinv(k)*flag

    ENDDO
    ENDDO
    ENDDO


   END SUBROUTINE calc_pres_gradz
