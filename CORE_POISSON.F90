!--------------------------------------------
!  SUBROUTINE rhs_poisson(sum) 
!  SUBROUTINE solve_poisson()
!  SUBROUTINE itsolv(var,r)
!  SUBROUTINE calc_residual(var,r,resm)
!--------------------------------------------

!-------------------------------------------------------------------------------
   SUBROUTINE rhs_poisson(sumd) 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE MPI_module

    IMPLICIT NONE
    
    REAL(KIND=CGREAL), INTENT(OUT)    :: sumd
    
    REAL(KIND=CGREAL) :: flag, sumd_i, resDiv_i
    INTEGER :: i,j,k

!******************************************************************************

    ! rhs = [ d(U)/dx + d(V)/dy + d(W)/dz ] / dt

    CALL enforce_global_mass_consv()

    sumd       = 0.0_CGREAL
    resDiv     = 0.0_CGREAL
    div(:,:,:) = 0.0_CGREAL

    DO k = zc_start,zc_end   !1,nz-1
    DO j = yc_start,yc_end
    DO i = 1,nx-1
      flag     = (1.0_CGREAL-REAL(iblank(i,j,k),       KIND=CGREAL))   &
               * (1.0_CGREAL-REAL(ghostCellMark(i,j,k),KIND=CGREAL))

      div(i,j,k)    = (dtinv)*                                       &
                     ( ( face_u(i+1,j,k) - face_u(i,j,k) )*dxinv(i)  &
                      +( face_v(i,j+1,k) - face_v(i,j,k) )*dyinv(j)  &
                      +( face_w(i,j,k+1) - face_w(i,j,k) )*dzinv(k)  )
      div(i,j,k)    = div(i,j,k)*flag

      sumd = sumd + div(i,j,k)*dx(i)*dy(j)*dz(k)

      resDiv = max(resDiv, abs(div(i,j,k)/dtinv) )
    ENDDO
    ENDDO
    ENDDO
    sumd  = sumd / dtinv  !Remove dtinv factor from massflux calculation

    !Correct massflux at ghostcells. 
    IF(mix_GC_form == 1) THEN
       CALL correct_ghost_massflux(sumd)
    ENDIF
    resDiv_i = resDiv

    ! Sum up massflux from all subdomains, and all processors will get the final sum
    sumd_i = sumd
    call MPI_ALLREDUCE(sumd_i, sumd, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, FLOW_COMM, ierr)       
    !print*, 'sumd_i and sumd of lProc ', lProc, 'are ', sumd_i, sumd

    ! Find the max Div among all subdomains.
    call MPI_ALLREDUCE(resDiv_i, resDiv, 1, MPI_DOUBLE_PRECISION, &
                          MPI_MAX, FLOW_COMM, ierr)
    !write(6,'(a,I5,a,E12.5,a,E12.5,a,E12.5,a,E12.5)'),'lProc=',lProc   &
    !       ,' sumd_i=',sumd_i,' sumd=',sumd,' resDiv_i=',resDiv_i,' resDiv=',resDiv

    IF (lProc==PROC_M .and. MOD(ntime,nmonitor) == 0 )  then
       print*,'RHS of Poisson equation: ', sumd
       print*,'Max of div = ', resDiv
    ENDIF

   END SUBROUTINE rhs_poisson
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
   SUBROUTINE solve_poisson()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE MPI_module

    IMPLICIT NONE

!... Local variables

    INTEGER              :: iter,i,j,loc(3),k
    REAL(KIND=CGREAL)    :: maxres, res_i, maxres2
    REAL(KIND=CGREAL)    :: pTotal, volTotal, pAverage
    !integer :: itr, clock1, clock2, clock_rate
!------------------

    !-- debugging the Poisson solver
    !call test_mg(0)
    !call sleep(50)
    !stop

   iter   = 0
   maxres = 1.0E10_CGREAL 
   maxres2 = 1.0E10_CGREAL

   IF (MOD(ntime,nmonitor) == 0) THEN

      CALL calc_residual(pPrime,div,maxres,loc)
      if(lProc==PROC_M) write(STDOUT,'(A,2(2X,1PE12.5),3(2X,I3))')  &
                       'Initial Residual Entering Poisson Solver:', &
                       maxres,restol_pson       !,loc(1),loc(2),loc(3)
   END IF ! ntime

   !For multigrid method, set up igmark for each level.
   if(it_solver_type == IT_SOLVER_TYPE_MG .and. boundary_motion == MOVING_BOUNDARY ) then
      call mg_prepare() 
   endif

   !DO WHILE ((iter .LT. itermax_pson) .AND. (maxres .GT. restol_pson))
   DO WHILE ((iter .LT. itermax_pson) .AND. (maxres2 .GT. restol_gcp))

      !------- subtract average pressure ----------
      IF(pbcx1==PBC_NEUMANN .and. pbcx2==PBC_NEUMANN .and.  &
         pbcy1==PBC_NEUMANN .and. pbcy2==PBC_NEUMANN .and.  &
         pbcz1==PBC_NEUMANN .and. pbcz2==PBC_NEUMANN        ) THEN
             
         if(lProc==0) then
            ! use the lower left corner pressure as the reference.
            pAverage = pPrime(1,1,1)     
         endif
         ! Broadcast pAverage to all processors from the 1st processor.
         ! Note that all processors need to execuate this call.
         call MPI_BCAST(pAverage, 1, MPI_DOUBLE_PRECISION, 0, FLOW_COMM, ierr)

         call MPI_BARRIER(FLOW_COMM,ierr)
         !All subdomains use the same reference pressure
         do k=zc_start-1,zc_end+1
         do j=yc_start-1,yc_end+1
         do i=0,nx
            pPrime(i,j,k) = pPrime(i,j,k) - pAverage
         enddo
         enddo
         enddo
      ENDIF ! if(pbcx1 ==

      CALL set_outer_pressure_bc(pPrime,nx,yb1,yb2,zb1,zb2,1)   ! inhomogeneous B.C.
      CALL set_outer_ghost_pressure(pPrime,nx,yb1,yb2,zb1,zb2,1)
      !CALL set_outer_ghost_pressure_2nd(pPrime)
      CALL GCM_enforce_p_compatibility(pPrime)

      call send_receive_slices_real_y(pPrime, 0,nx+1,yb1,yb2,zb1,zb2,2) ! exchange 2 slices      
      call send_receive_slices_real_z(pPrime, 0,nx+1,yb1,yb2,zb1,zb2,2)

      ! Save the pressure from the previous iteration.
      ! This will be used for residual calculation.
      p_temp = pPrime

      CALL GCM_ghostcell_pressure(pPrime,div)

      SELECT CASE(it_solver_type)

      CASE (IT_SOLVER_TYPE_LSOR)
         CALL itsolv(pPrime,div)

      CASE (IT_SOLVER_TYPE_MG) 
         CALL mg_solve(pPrime,div)   

      CASE (IT_SOLVER_TYPE_AZ)
         print*, 'AZ solve not implemented. STOP'
         stop
      END SELECT

      ! Residual based on satisfaction of the Poisson Eq.
      CALL calc_residual(pPrime,div,maxres,loc)

      !Residual based on the difference between iterations
      CALL calc_residual2(pPrime,div,maxres2,loc)

      iter = iter + 1

      IF(MOD(ntime,nmonitor) == 0 .and. lProc==PROC_M) THEN
         write(STDOUT,'(A,I4,2X,1PE15.7,2X,1PE15.7)')  &  !4(2X,I4)
              'Pressure Convergence : ',       &
               iter,maxres,maxres2 !loc(1),loc(2),loc(3),ghostcellMark(loc(1),loc(2),loc(3))
      ENDIF
   ENDDO  ! end do while

   !IF (lProc==PROC_M .and. iter.EQ.itermax_pson .AND. maxres .GT. restol_pson ) THEN
   IF (lProc==PROC_M .and. iter.EQ.itermax_pson .AND. maxres2 .GT. restol_gcp) THEN
      PRINT*,'Pressure did not converge in ',itermax_pson,' iterations'
      PRINT*,'Final residual = ',maxres
   ENDIF

   END SUBROUTINE solve_poisson

!----------------------------------------------------------
! Line SOR with Gauss Siedel as smoother
!----------------------------------------------------------
   SUBROUTINE itsolv(var,r)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MPI_module

    IMPLICIT NONE

!... parameters

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2),  INTENT (IN OUT) ::var
    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2),  INTENT (IN)     ::r

!... Local variables

    REAL(KIND=CGREAL), DIMENSION(nx-1) :: FF

    REAL(KIND=CGREAL) :: VolCell, flag
    INTEGER :: i,j,k
    INTEGER :: iBody, iRow, iG,jG,kG, iNode, jNode,kNode, n

    INTEGER :: iinit, istep, jinit, jstep, Ncolor ,kinit, kstep
    INTEGER :: kend,k0,k1,jend,j0,j1
!******************************************************************************

! Line solve in the x-direction
  DO Ncolor = 1, 1!2

    !Ye,for 2d_decomposition
    !jinit = mod(Ncolor+1,2) + 1  ! color scheme in y
    !jstep = 2
    jinit = yc_start
    jstep = 1
    jend  = yc_end
    kinit = zc_start
    kstep = 1
    kend  = zc_end

    !CALL set_outer_pressure_bc   (pPrime,nx,ny,zb1,zb2,1)   ! inhomogeneous B.C.
    !CALL set_outer_ghost_pressure(pPrime,nx,ny,zb1,zb2,1)

    DO k = kinit, kend, kstep 
    DO j = jinit, jend, jstep

      DO i = 1,nx-1
       amx(i) =   dxcinv(i)  *dxinv(i)
       apx(i) =   dxcinv(i+1)*dxinv(i)
       acx(i) = - ( amx(i) + apx(i) )

       amy(j) =   dycinv(j)  *dyinv(j)
       apy(j) =   dycinv(j+1)*dyinv(j)
       acy(j) = - ( amy(j) + apy(j) )

       amz(k) =   dzcinv(k)  *dzinv(k)
       apz(k) =   dzcinv(k+1)*dzinv(k)
       acz(k) = - ( amz(k) + apz(k) )

       rhs(i) = r(i,j,k) - var(i,j-1,k)*amy(j)  &
                         - var(i,j+1,k)*apy(j)  &               
                         - var(i,j,k-1)*amz(k)  & 
                         - var(i,j,k+1)*apz(k)

       flag   = (1.0_CGREAL-REAL(iblank(i,j,k),       KIND=CGREAL))   &
              * (1.0_CGREAL-REAL(ghostCellMark(i,j,k),KIND=CGREAL))   

       amx(i) = amx(i)*flag
       apx(i) = apx(i)*flag
       acx(i) = (acx(i)+acy(j)+acz(k))*flag + (1.0_CGREAL - flag)
         
       rhs(i) = rhs(i)*flag + (1.0_CGREAL - flag)*var(i,j,k)
      ENDDO ! i 

      ! set up trivial equations for the outer boundaries
      amx(0)  = 0.0_CGREAL
      apx(0)  = 0.0_CGREAL
      acx(0)  = 1.0_CGREAL
      rhs(0)  = var(0, j,k)
      amx(nx) = 0.0_CGREAL
      apx(nx) = 0.0_CGREAL
      acx(nx) = 1.0_CGREAL
      rhs(nx) = var(nx,j,k)                                 

      CALL tdma(amx,acx,apx,rhs,dummy,0,nx)

      DO i=1,nx-1
         var(i,j,k) = var(i,j,k) + omega_pson*(dummy(i)-var(i,j,k))
      ENDDO
    ENDDO ! j 
    ENDDO ! k
  ENDDO ! Ncolor

  !Exchange 1 slice of buffer between subdomains
  call send_receive_slices_real_y(var, 0,nx+1,yb1,yb2,zb1,zb2,1)
  call send_receive_slices_real_z(var, 0,nx+1,yb1,yb2,zb1,zb2,1)

! Line solve in the y-direction
  !Ye, for 2d-decomposition, keep the Ncolor
  DO Ncolor = 1, 2             
    
    iinit = mod(Ncolor+1,2) + 1   ! color scheme in x
    istep = 2
    kinit = zc_start  
    kstep = 1
    kend  = zc_end
    j0    = yc_start-1
    j1    = yc_end  +1

    !CALL set_outer_pressure_bc   (pPrime,nx,ny,zb1,zb2,1)   ! inhomogeneous B.C.
    !CALL set_outer_ghost_pressure(pPrime,nx,ny,zb1,zb2,1)

    DO k = kinit,kend,kstep
    DO i = iinit,nx-1,istep
      !DO j = 0,ny  
      DO j = yc_start,yc_end 
       amx(i) =   dxcinv(i)  *dxinv(i)
       apx(i) =   dxcinv(i+1)*dxinv(i)
       acx(i) = - ( amx(i) + apx(i) )

       amy(j) =   dycinv(j)  *dyinv(j)
       apy(j) =   dycinv(j+1)*dyinv(j)
       acy(j) = - ( amy(j) + apy(j) )

       amz(k) =   dzcinv(k)  *dzinv(k)
       apz(k) =   dzcinv(k+1)*dzinv(k)
       acz(k) = - ( amz(k) + apz(k) )

       rhs(j) = r(i,j,k) - var(i-1,j,k)*amx(i)   &  
                         - var(i+1,j,k)*apx(i)   &
                         - var(i,j,k-1)*amz(k)   &
                         - var(i,j,k+1)*apz(k)   

!-- Modify rhs and coefficients for Ghost cells in GCM
       flag   = (1.0_CGREAL-REAL(iblank(i,j,k),       KIND=CGREAL))   &
              * (1.0_CGREAL-REAL(ghostCellMark(i,j,k),KIND=CGREAL))

       amy(j) =  amy(j)*flag
       apy(j) =  apy(j)*flag
       acy(j) = (acx(i)+acy(j)+acz(k))*flag + (1.0_CGREAL - flag)      
       rhs(j) =  rhs(j)*flag + (1.0_CGREAL - flag)*var(i,j,k) 

      ENDDO

      ! set up trivial equations for the outer boundaries
      amy(j0)  = 0.0_CGREAL
      apy(j0)  = 0.0_CGREAL
      acy(j0)  = 1.0_CGREAL
      rhs(j0)  = var(i, yc_start-1,k)                                          
      amy(j1) = 0.0_CGREAL
      apy(j1) = 0.0_CGREAL
      acy(j1) = 1.0_CGREAL
      rhs(j1) = var(i,yc_end+1,k)                                          
                                         
      CALL tdma(amy,acy,apy,rhs,dummy,j0,j1)

      DO j=yc_start,yc_end
         var(i,j,k) = var(i,j,k) + omega_pson*(dummy(j)-var(i,j,k))
      ENDDO
    ENDDO
    ENDDO
  ENDDO ! Ncolor

  !Exchange 1 slice of buffer between subdomains
  call send_receive_slices_real_y(var, 0,nx+1,yb1,yb2,zb1,zb2,1)
  call send_receive_slices_real_z(var, 0,nx+1,yb1,yb2,zb1,zb2,1)

! Line solve in the z-direction
  DO Ncolor = 1, 1!2        ! color scheme in y     
    
    iinit = 1  !mod(Ncolor+1,2) + 1
    istep = 1
    !jinit = mod(Ncolor+1,2) + 1
    !jstep = 2
    !Ye,2d-decomposition
    jinit = yc_start
    jstep = 1
    jend  = yc_end
    k0    = zc_start-1
    k1    = zc_end  +1 

    !CALL set_outer_pressure_bc   (pPrime,nx,ny,zb1,zb2,1)   ! inhomogeneous B.C.
    !CALL set_outer_ghost_pressure(pPrime,nx,ny,zb1,zb2,1)

    DO j = jinit,jend,jstep
    DO i = iinit,nx-1,istep

      DO k = zc_start,zc_end

       amx(i) =   dxcinv(i)  *dxinv(i)
       apx(i) =   dxcinv(i+1)*dxinv(i)
       acx(i) = - ( amx(i) + apx(i) )

       amy(j) =   dycinv(j)  *dyinv(j)
       apy(j) =   dycinv(j+1)*dyinv(j)
       acy(j) = - ( amy(j) + apy(j) )

       !kk     = k-zc_start+1   !shift the 1D arrays, but NOT the others
       amz(k) =   dzcinv(k)  *dzinv(k)
       apz(k) =   dzcinv(k+1)*dzinv(k)
       acz(k) = - ( amz(k) + apz(k) )

       rhs(k) = r(i,j,k) - var(i-1,j,k)*amx(i) &  
                         - var(i+1,j,k)*apx(i)   &
                         - var(i,j-1,k)*amy(j)   &
                         - var(i,j+1,k)*apy(j)   

!-- Modify rhs and coefficients for Ghost cells in GCM
       flag   = (1.0_CGREAL-REAL(iblank(i,j,k),       KIND=CGREAL))   &
              * (1.0_CGREAL-REAL(ghostCellMark(i,j,k),KIND=CGREAL))

       amz(k) = amz(k)*flag
       apz(k) = apz(k)*flag
       acz(k) = (acx(i)+acy(j)+acz(k))*flag + (1.0_CGREAL - flag)      
       rhs(k) = rhs(k)*flag + (1.0_CGREAL - flag)*var(i,j,k) 

      ENDDO

      ! set up trivial equations for the outer boundaries of the subdomain
      amz(k0)  = 0.0_CGREAL
      apz(k0)  = 0.0_CGREAL
      acz(k0)  = 1.0_CGREAL
      rhs(k0)  = var(i, j,zc_start-1)                                          
      amz(k1) = 0.0_CGREAL
      apz(k1) = 0.0_CGREAL
      acz(k1) = 1.0_CGREAL
      rhs(k1) = var(i,j,zc_end+1)                                          

      ! Note that these arrays have actual size of 0:nz+1, but only 
      ! zc_start-1:zc_end+1 are being used. It is very important that
      ! the arrays in TDMA have the same starting index, i.e., 0, even
      ! though the ending index is may be shorter.
      CALL tdma(amz,acz,apz,rhs,dummy,k0,k1) 

      DO k=zc_start,zc_end   !1,nz-1
         !kk     = k-zc_start+1
         var(i,j,k) = var(i,j,k) + omega_pson*(dummy(k)-var(i,j,k))
      ENDDO

    ENDDO
    ENDDO
  ENDDO ! Ncolor

  !Exchange 1 slice of buffer between subdomains
  call send_receive_slices_real_y(var, 0,nx+1,yb1,yb2,zb1,zb2,1)
  call send_receive_slices_real_z(var, 0,nx+1,yb1,yb2,zb1,zb2,1)

  END SUBROUTINE itsolv
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
   SUBROUTINE calc_residual(var,r,resm,loc)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MPI_module

    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2),INTENT (IN)  ::var,r
    REAL(KIND=CGREAL),                            INTENT (OUT) ::resm
    INTEGER,           DIMENSION(3),              INTENT (OUT) ::loc

    INTEGER              :: i,j,k
    INTEGER              :: iG,jG,kG,iBody,iRow,n
    REAL(KIND=CGREAL)    :: res, res2, volCell,flag
    REAL(KIND=CGREAL) :: bmx,bpx,bcx,bc
    REAL(KIND=CGREAL) :: bmy,bpy,bcy
    REAL(KIND=CGREAL) :: bmz,bpz,bcz

!******************************************************************************

    resm = 0.0_CGREAL
    res2 = 0.0_CGREAL

    DO k=zc_start,zc_end  !1,nz-1
    DO j=yc_start,yc_end
    DO i=1,nx-1

       bmx =   dxcinv(i)  *dxinv(i)
       bpx =   dxcinv(i+1)*dxinv(i)
       bcx = - ( bmx + bpx )
  
       bmy =   dycinv(j)  *dyinv(j)
       bpy =   dycinv(j+1)*dyinv(j)
       bcy = - ( bmy + bpy )

       bmz =   dzcinv(k)  *dzinv(k)
       bpz =   dzcinv(k+1)*dzinv(k)
       bcz = - ( bmz + bpz )


       ! Skip the blank-cells and ghost-cells. Note that dead_cell = -1
       flag     = (1.0_CGREAL-REAL(iblank(i,j,k),       KIND=CGREAL))   &
                * (1.0_CGREAL-REAL(ghostCellMark(i,j,k),KIND=CGREAL))  !&
               !* (1.0_CGREAL-REAL(dead_cell(i,j,k),    KIND=CGREAL))
     

       bc  = (bcx+bcy+bcz)

       res    = r(i,j,k) - var(i,j,k)*bc                 &
                       - var(i-1,j,k)*bmx                &
                       - var(i+1,j,k)*bpx                &
                       - var(i,j-1,k)*bmy                &
                       - var(i,j+1,k)*bpy                &
                       - var(i,j,k-1)*bmz                &
                       - var(i,j,k+1)*bpz                

       ! skip the residual check at ghost cells and solid cells
       res = res * flag

       IF (ABS(res) > resm ) THEN
         resm = ABS(res)     !L-inf norm
         loc(1) = i
         loc(2) = j
         loc(3) = k
       
       ENDIF 
       res2 = res2 + res**2  ! L2-norm

    ENDDO
    ENDDO
    ENDDO

    ! Find the infinity norm
    !res2 = resm
    !call MPI_ALLREDUCE(res2, resm, 1, MPI_DOUBLE_PRECISION, &   ! change from res to res2
    !                   MPI_MAX, flow_comm, ierr)
    
    !Find the L-2 norm
    call MPI_ALLREDUCE(res2, resm, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, flow_comm, ierr)
    resm = sqrt( resm/dfloat(nx*ny*nz) )
    
    !print*,lProc, res2, resm

   END SUBROUTINE calc_residual
!-------------------------------------------------------------------------------
!
! Use the following subroutine to test the MG Poisson solver using a 
! user specified problem.
!
!  d^2     d^2    d^2   
! (---- + ---- + ----)P = sum[sin(k*x)*sin(k*y)*sin(k*z)]_{k=1,2,3,4}
!  dx^2   dy^2   dz^2
!
! To run this test, simply call this subroutine somewhere in
! the code, e.g., within solve_poisson. Enjoy watching the fast convergence!
!
!                               --H. Luo, July 2015
   SUBROUTINE test_mg(iOpt)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MPI_module

    IMPLICIT NONE
    INTEGER :: iOpt

!... Local variables
    INTEGER :: I,J,K,L,LL
    INTEGER :: iinit, istep, jinit, jstep, Ncolor ,kinit, kstep
    REAL(KIND=CGREAL) :: PI2, xx, yy, zz

    INTEGER              :: iter,loc(3)
    REAL(KIND=CGREAL)    :: maxres, res_i, pAverage
!**************************************************************

    PI2 = 8 * atan (1.0_8)

    !r.h.s. forcing term
    DO I = 1,nx-1
    DO J = yc_start, yc_end
    DO L = zc_start, zc_end

       div(I,J,L) = 0.0d0
       
       xx = (xc(I)-xOrigin)/(xOut-xOrigin)*PI2 
       yy = (yc(J)-yOrigin)/(yOut-yOrigin)*PI2
       zz = (zc(L)-zOrigin)/(zOut-zOrigin)*PI2

       DO K = 1,4
          div(I,J,L)=div(I,J,L) + 10*cos(xx*K)*cos(yy*K)*cos(zz*K)
       END DO
    END DO
    END DO
    END DO
    
    pPrime = 0.0d0

    !comment out the following 3 lines 
    !to test immersed boundary

    iblank(:,:,:) = 0           ! clear iblank and ghostcells
    ghostcellmark(:,:,:) = 0    
    nGhost = 0; nDead = 0

    iter   = 0
    maxres = 1.0E10_CGREAL 

   !For multigrid method, set up igmark for each level.
   if(it_solver_type == IT_SOLVER_TYPE_MG)  call mg_prepare() 

   DO WHILE ((iter .LT. itermax_pson) .AND. (maxres .GT. restol_pson))
      if(iOpt == 1) then
        CALL set_outer_pressure_bc   (pPrime,nx,yb1,yb2,zb1,zb2,1)   ! inhomogeneous B.C.
        CALL set_outer_ghost_pressure(pPrime,nx,yb1,yb2,zb1,zb2,1)
      endif

      SELECT CASE(it_solver_type)

      CASE (IT_SOLVER_TYPE_LSOR)
         CALL itsolv(pPrime,div)

      CASE (IT_SOLVER_TYPE_MG) 
         CALL mg_solve(pPrime,div)   

      END SELECT

      call send_receive_slices_real_y(pPrime, 0,nx+1,yb1,yb2,zb1,zb2,2) ! exchange 2 slices      
      call send_receive_slices_real_z(pPrime, 0,nx+1,yb1,yb2,zb1,zb2,2)
      CALL GCM_ghostcell_pressure(pPrime,div)
      CALL calc_residual(pPrime,div,maxres,loc)
      iter = iter + 1

      IF(lProc==PROC_M) THEN
         write(STDOUT,'(A,I4,2X,1PE15.7,4(2X,I4))')  &
              'Pressure Convergence : ',iter,maxres
      ENDIF
   ENDDO  ! end do while

         if(lProc==0) then
            pAverage = pPrime(1,1,1)     ! use the lower left corner pressure as the reference.
         endif
         ! Broadcast pAverage to all processors from the 1st processor.
         ! Note that all processors need to execuate this call.
         call MPI_BCAST(pAverage, 1, MPI_DOUBLE_PRECISION, 0, FLOW_COMM, ierr)

         call MPI_BARRIER(FLOW_COMM,ierr)
         !All subdomains use the same reference pressure
         do k=zc_start-1,zc_end+1
         do j=yc_start-1,yc_end+1
         do i=0,nx
            pPrime(i,j,k) = pPrime(i,j,k) - pAverage
         enddo
         enddo
         enddo

   p = pPrime
   lmd = div

   call write_subdomain()

   END SUBROUTINE test_mg
!-------------------------------------------------------------------------------
! This subroutine calculate the residual of pressure using difference
! of two consecutive iterations; L-inf norm is used.

   SUBROUTINE calc_residual2(var,r,resm,loc)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE MPI_module

    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(0:nx+1,yb1:yb2,zb1:zb2),INTENT (IN)  ::var,r
    REAL(KIND=CGREAL),                            INTENT (OUT) ::resm
    INTEGER,           DIMENSION(3),              INTENT (OUT) ::loc

    INTEGER              :: i, j, k
    REAL(KIND=CGREAL)    :: res, res1, flag

    resm = 0.0_CGREAL
    res1 = 0.0_CGREAL

    DO k=zc_start,zc_end
    DO j=yc_start,yc_end
    DO i=1,nx-1
       flag = (1.0_CGREAL-REAL(iblank(i,j,k),       KIND=CGREAL))
       res  = ABS(var(i,j,k)-p_temp(i,j,k))*flag

       if ( res > res1 )then
          res1 = res
          loc(1) = i
          loc(2) = j
          loc(3) = k
       endif

    ENDDO
    ENDDO
    ENDDO

    !Call MPI reduction to get the max residual from all processors
    call MPI_ALLREDUCE(res1, resm, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, flow_comm, ierr)

    END SUBROUTINE calc_residual2
