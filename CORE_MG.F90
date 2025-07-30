!--------------------------------------------------------------------
! This files contains all subroutines for multi-grid method.
! It includes:
!  SUBROUTINE MG_initial() 
!  SUBROUTINE itsolv_mg()
!  SUBROUTINE residual_MG()
!  SUBROUTINE mg_restrict()
!  SUBROUTINE mg_prolong()
!  SUBROUTINE mg_solver() 
!  SUBROUTINE mg_solv_x, mg_solv_y, mg_solv_z
!  SUBROUTINE mg_prepare()
!---------------------------------------------------------------------
!
! Rewritten by H. Luo, July 2015
!
!---------------------------------------------------------------------


!--------------------------------------------------------
!Purpose:  Initializing necessary arrays for MG method.
!--------------------------------------------------------
   SUBROUTINE MG_initial() 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE mg_module
    USE MPI_module

    IMPLICIT NONE

    INTEGER :: I, J, K, N, M, ii, iim, jj, jjm, kk, kkm, ncoargrids

!******************************************************************************

!---------
! Parameters for alternating sweep scheme 
! for the line-SOR iteration
!---------

    IF (iRedBlack == 0) THEN
       TNcolorX = 1
       TNcolorY = 1
       TNcolorz = 1
       iStep = 1
       jStep = 1
       kStep = 1
    ELSE IF(iRedBlack == 1) THEN
       TNcolorX = 2
       TNcolorY = 2
       TNcolorZ = 1
       
       iStep = 1
       jStep = 2
       kStep = 2
    ELSE 
       TNcolorX = 2
       TNcolorY = 2
       TNcolorZ = 2
       
       iStep = 2
       jStep = 2
       kStep = 2
    END IF

    ncoargrids = 4    ! minimum points

    IF (mgrids_x .EQ. 0 .or. mgrids_x .gt. max_grid_level) THEN
       II = NX-1
       DO N = 1, max_grid_level
          II = II/2
          IF (II.GE.ncoargrids) CYCLE   ! keep going
          EXIT 
       END DO
       mgrids_x = N
    END IF
    if(lProc==PROC_M) WRITE(*,*) 'X-multigrid level is:',mgrids_x 
    
    IF (mgrids_y .EQ. 0 .or. mgrids_y .gt. max_grid_level) THEN
       JJ = NY-1
       DO N = 1, max_grid_level
          JJ = JJ/2
          IF (JJ.GE.ncoargrids) CYCLE 
          EXIT  
       END DO
       mgrids_y = N
    END IF
    if(lProc==PROC_M) WRITE(*,*) 'Y-multigrid level is:', mgrids_y

    IF (mgrids_z .EQ. 0 .or. mgrids_z .gt. max_grid_level) THEN
       KK = NZ-1
       DO N = 1, max_grid_level
          KK = KK/2
          IF (KK.GE.ncoargrids) CYCLE
          EXIT
       END DO
       mgrids_z = N
    END IF
    if(lProc==PROC_M) WRITE(*,*) 'Z-multigrid level is:', mgrids_z

    !Number of levels for simultaneous coarsening
    mLevel = max(mgrids_x,mgrids_y,mgrids_z)

    !Number of points at each level (whole domain)
    ALLOCATE( mgrid_I(mLevel), mgrid_J(mLevel),mgrid_K(mLevel))

    !start/end and number of slices at each level (subdomain)
    !Z-direction
    ALLOCATE( z_start_mg(mLevel),  z_end_mg(mLevel))
    ALLOCATE( zc_start_mg(mLevel), zc_end_mg(mLevel))
    ALLOCATE( zb1_mg(mLevel), zb2_mg(mLevel) )
    ALLOCATE( num_slicesZ_mg(mLevel), all_slicesZ_mg(0:nProcs-1,mLevel))

    !left and right processors IDs (some processors have no slices)
    ALLOCATE( l_proc_mg(mLevel),   r_proc_mg(mLevel) )
    ALLOCATE( leftmost_mg(mLevel), rightmost_mg(mLevel) )

    !start/end and number of slices at each level (subdomain)
    !Y-direction
    ALLOCATE( y_start_mg(mLevel),  y_end_mg(mLevel))
    ALLOCATE( yc_start_mg(mLevel), yc_end_mg(mLevel))
    ALLOCATE( yb1_mg(mLevel), yb2_mg(mLevel) )
    ALLOCATE( num_slicesY_mg(mLevel), all_slicesY_mg(0:nProcs-1,mLevel))

    !bottom and top processors IDs (some processors have no slices)
    ALLOCATE( b_proc_mg(mLevel),   t_proc_mg(mLevel) )
    ALLOCATE( btommost_mg(mLevel), topmost_mg(mLevel) )

    !data arrays at each level
    ALLOCATE( mg_array(mLevel))
!---------------- Calculate the number of grids at each level -------------------

    ! Finest grid, or the original grid
    mgrid_I(1) = nx
    mgrid_J(1) = ny
    mgrid_K(1) = nz
    
    ! Other grids
    DO n = 2, mLevel 
       if(n <= mgrids_x) then
          ! Coarsening by half
          IF(mod(mgrid_I(n-1),2) .EQ. 0) THEN
             mgrid_I(n) = mgrid_I(n-1)/2 + 1
          ELSE
             mgrid_I(n) = (mgrid_I(n-1) + 1)/2
          END IF
       else
          mgrid_I(n) = mgrid_I(n-1)  ! no further coarsening; use the same grid
       endif

       if(n <= mgrids_y) then
          ! Coarsening by half
          IF(mod(mgrid_J(n-1),2) .EQ. 0) THEN
             mgrid_J(n) = mgrid_J(n-1)/2 + 1
          ELSE
             mgrid_J(n) = (mgrid_J(n-1) + 1)/2
          END IF
       else
          mgrid_J(n) = mgrid_J(n-1)  ! no further coarsening; use the same grid
       endif

       if(n <= mgrids_z) then
          ! Coarsening by half
          IF(mod(mgrid_K(n-1),2) .EQ. 0) THEN
             mgrid_K(n) = mgrid_K(n-1)/2 + 1
          ELSE
             mgrid_K(n) = (mgrid_K(n-1) + 1)/2
          END IF
       else
          mgrid_K(n) = mgrid_K(n-1)  ! no further coarsening; use the same grid
       endif
    END DO ! end n

    ! Grid array allocation
    DO n = 1, mLevel
       ! allocate grid arrays for each level
       ALLOCATE( mg_array(n)%x( 0:mgrid_I(n)+1 ) )
       ALLOCATE( mg_array(n)%y( 0:mgrid_J(n)+1 ) )
       ALLOCATE( mg_array(n)%z( 0:mgrid_K(n)+1 ) )

       ALLOCATE( mg_array(n)%xc( 0:mgrid_I(n)+1 ) )
       ALLOCATE( mg_array(n)%yc( 0:mgrid_J(n)+1 ) )
       ALLOCATE( mg_array(n)%zc( 0:mgrid_K(n)+1 ) )

       ALLOCATE( mg_array(n)%dxinv( 0:mgrid_I(n)+1 ) )
       ALLOCATE( mg_array(n)%dyinv( 0:mgrid_J(n)+1 ) )
       ALLOCATE( mg_array(n)%dzinv( 0:mgrid_K(n)+1 ) )

       ALLOCATE( mg_array(n)%dxcinv( 0:mgrid_I(n)+1 ) )
       ALLOCATE( mg_array(n)%dycinv( 0:mgrid_J(n)+1 ) )
       ALLOCATE( mg_array(n)%dzcinv( 0:mgrid_K(n)+1 ) )
 
    END DO ! end n

    !----------- Calculate xc, yc, zc, dxcinv, dycinv, dzcinv etc. -------
    
    ! Finest grid: copying from the original grid.
    DO I = 0, nx
       mg_array(1)%dxcinv(i) = dxcinv(i)
       mg_array(1)%dxinv(i) = dxinv(i)
       mg_array(1)%xc(i) = xc(i)
    END DO

    DO J = 0, ny 
       mg_array(1)%dycinv(j) = dycinv(j)
       mg_array(1)%dyinv(j) = dyinv(j)
       mg_array(1)%yc(j) = yc(j) 
    END DO

    DO k = 0, nz 
       mg_array(1)%dzcinv(k) = dzcinv(k)
       mg_array(1)%dzinv(k) = dzinv(k)
       mg_array(1)%zc(k) = zc(k) 
    END DO

    DO I = 0, nx+1
       mg_array(1)%x(i) = x(i)
    END DO
    DO J = 0, ny+1
       mg_array(1)%y(j) = y(j) 
    END DO
    DO k = 0, nz+1
       mg_array(1)%z(k) = z(k) 
    END DO


    !==== Other grids ========
    DO n = 2, mLevel

       ! x-grid
       DO i=1,mgrid_I(n)-1
          if( mgrid_I(n) .eq. mgrid_I(n-1)) then
             iim = i
          else
             iim = 2*i-1
          endif
          mg_array(n)%x(i) = mg_array(n-1)%x(iim)
       ENDDO
       ! The last node is located at the outer boundary
       mg_array(n)%x( mgrid_I(n) ) = mg_array(n-1)%x( mgrid_I(n-1) )
       if(lProc==PROC_M) then
          write(111, '(a,I5)') 'x-grid, n = ', n, 'mgrid_I(n) = ',mgrid_I(n)
          write(111, '(500F12.5)') mg_array(n)%x( 1:mgrid_I(n) )
       endif

       DO i=1,mgrid_I(n)-1
          mg_array(n)%dxinv(i) = 1.0_CGREAL/(mg_array(n)%x(i+1)-mg_array(n)%x(i))

          mg_array(n)%xc(i) = 0.5_CGREAL*(mg_array(n)%x(i) + mg_array(n)%x(i+1))
       ENDDO

       mg_array(n)%xc(0)         = 2.0_CGREAL * x(1 ) - mg_array(n)%xc(1)
       mg_array(n)%xc(mgrid_I(n))= 2.0_CGREAL * x(nx) - mg_array(n)%xc(mgrid_I(n)-1)
       if(lProc==PROC_M) then
          write(111, '(a,I5)') 'xc-grid, n = ', n, '0:mgrid_I(n) = ',mgrid_I(n)+1
          write(111, '(500F12.5)') mg_array(n)%xc( 0:mgrid_I(n))
       endif

       DO i= 1,mgrid_I(n)
          mg_array(n)%dxcinv(i) = 1.0_CGREAL/(mg_array(n)%xc(i)-mg_array(n)%xc(i-1))
       ENDDO

       ! y-grid
       DO j=1,mgrid_J(n)-1
          if( mgrid_J(n) .eq. mgrid_J(n-1)) then
             jjm = j
          else
             jjm = 2*j-1
          endif
          mg_array(n)%y(j) = mg_array(n-1)%y(jjm)
       ENDDO
       ! The last node is located at the outer boundary
       mg_array(n)%y( mgrid_J(n) ) = mg_array(n-1)%y( mgrid_J(n-1) )
       if(lProc==PROC_M) then
          write(111, '(a,I5)') 'y-grid, n = ', n, 'mgrid_J(n) = ',mgrid_J(n)
          write(111, '(500F12.5)') mg_array(n)%y( 1:mgrid_J(n) )
       endif

       DO j=1,mgrid_J(n)-1
          mg_array(n)%dyinv(j) = 1.0_CGREAL/(mg_array(n)%y(j+1)-mg_array(n)%y(j)) 

          mg_array(n)%yc(j) = 0.5_CGREAL*(mg_array(n)%y(j)+mg_array(n)%y(j+1))
       ENDDO

       mg_array(n)%yc(0)         = 2.0_CGREAL * y(1 ) - mg_array(n)%yc(1)
       mg_array(n)%yc(mgrid_J(n))= 2.0_CGREAL * y(ny) - mg_array(n)%yc(mgrid_J(n)-1)
       if(lProc==PROC_M) then
          write(111, '(a,I5)') 'yc-grid, n = ', n, '0:mgrid_J(n) = ',mgrid_J(n)+1
          write(111, '(500F12.5)') mg_array(n)%yc( 0:mgrid_J(n))
       endif

       DO j= 1,mgrid_J(n)
           mg_array(n)%dycinv(j) = 1.0_CGREAL/( mg_array(n)%yc(j)- mg_array(n)%yc(j-1))
       ENDDO

       ! z-grid
       DO k=1,mgrid_K(n)-1
          if( mgrid_K(n) .eq. mgrid_K(n-1)) then
             kkm = k
          else
             kkm = 2*k-1
          endif
          mg_array(n)%z(k) = mg_array(n-1)%z(kkm)
       ENDDO
       ! The last node is located at the outer boundary
       mg_array(n)%z( mgrid_K(n) ) = mg_array(n-1)%z( mgrid_K(n-1) )
       if(lProc==PROC_M) then
          write(111, '(a,I5)') 'z-grid, n = ', n, 'mgrid_K(n) = ',mgrid_K(n)
          write(111, '(500F12.5)') mg_array(n)%z( 1:mgrid_K(n) )
       endif

       DO k=1,mgrid_K(n)-1
          mg_array(n)%dzinv(k) = 1.0_CGREAL/(mg_array(n)%z(k+1)-mg_array(n)%z(k)) 

          mg_array(n)%zc(k) = 0.5_CGREAL*(mg_array(n)%z(k)+mg_array(n)%z(k+1))
       ENDDO
      
       mg_array(n)%zc(0)         = 2.0_CGREAL * z(1 ) - mg_array(n)%zc(1)
       mg_array(n)%zc(mgrid_K(n))= 2.0_CGREAL * z(nz) - mg_array(n)%zc(mgrid_K(n)-1)
       if(lProc==PROC_M) then
          write(111, '(a,I5)') 'zc-grid, n = ', n, '0:mgrid_K(n) = ',mgrid_K(n)+1
          write(111, '(500F12.5)') mg_array(n)%zc( 0:mgrid_K(n))
       endif

       DO k= 1,mgrid_K(n)
          mg_array(n)%dzcinv(k) = 1.0_CGREAL/(mg_array(n)%zc(k)-mg_array(n)%zc(k-1))
       ENDDO

    END DO ! enddo for n

    !----
    !Determine slice assignment on multigrid for MPI parallelization
    !Note that some processors will have 0 slices.
    !----

    !--- at level 1
    !y-direction
    y_start_mg (1) = y_start
    yc_start_mg(1) = yc_start

    y_end_mg (1)  = y_end
    yc_end_mg(1)  = yc_end

    yb1_mg(1)     = yb1
    yb2_mg(1)     = yb2

    !z-direction
    z_start_mg (1) = z_start
    zc_start_mg(1) = zc_start

    z_end_mg (1)  = z_end
    zc_end_mg(1)  = zc_end

    zb1_mg(1)     = zb1
    zb2_mg(1)     = zb2

    !y-direction
    l_proc_mg(1) = lProc_left
    r_proc_mg(1) = lProc_right
    leftmost_mg (1) = 0         !nProcY*kProc
    rightmost_mg(1) = nProcY-1  !nProcY*kProc+(nProcY-1)
    num_slicesY_mg(1)= jSlices

    !z-direction
    b_proc_mg(1) = lProc_bottom
    t_proc_mg(1) = lProc_top
    btommost_mg (1) = 0         !jProc
    topmost_mg(1)   = nProcZ-1  !nProcY*(nProcZ-1)+jProc
    num_slicesZ_mg(1)= kSlices

    !--- other levels
    DO n = 2, mLevel
 
       !Y-direction
       !If no coarsening in y, copy all slice info from the finer level
       if( mgrid_J(n) .eq. mgrid_J(n-1)) then
          y_start_mg (n) = y_start_mg (n-1)
          y_end_mg   (n) = y_end_mg   (n-1)
          yc_start_mg(n) = yc_start_mg(n-1)
          yc_end_mg  (n) = yc_end_mg  (n-1)

          yb1_mg     (n) = yb1_mg     (n-1)
          yb2_mg     (n) = yb2_mg     (n-1)

          num_slicesY_mg(n) = num_slicesY_mg(n-1)

          l_proc_mg(n)    = l_proc_mg(n-1)
          r_proc_mg(n)    = r_proc_mg(n-1)
          leftmost_mg (n) = leftmost_mg (n-1)
          rightmost_mg(n) = rightmost_mg(n-1)

          cycle  ! go to next level
       endif

      !Otherwise, the grid is coarsened. Redo the slice assignment
      if( mod(y_start_mg(n-1), 2)==0 )  then
         y_start_mg(n) = y_start_mg(n-1)/2 + 1
      else
         y_start_mg(n) = (y_start_mg(n-1)+1)/2
      endif

      if( mod(y_end_mg(n-1), 2)==0 )  then
         y_end_mg(n) = y_end_mg(n-1)/2 + 1
      else
         y_end_mg(n) = (y_end_mg(n-1)+1)/2
      endif
      num_slicesY_mg(n) = y_end_mg(n) - y_start_mg(n)

      yc_start_mg(n) = y_start_mg(n)
      yc_end_mg(n)   = y_end_mg(n) - 1

      yb1_mg(n)      = yc_start_mg(n) - 2
      yb2_mg(n)      = yc_end_mg(n) + 2

      write(*,'(a,I2,a,I4,a,7I4)')'Lev:',n,', id:',lProc,', Ystart/end/slices:',&
               yb1_mg(n),y_start_mg(n),yc_start_mg(n),yc_end_mg(n),y_end_mg(n),&
               yb2_mg(n),num_slicesY_mg(n)

      !Each processor gathers slice info from all others
      call MPI_ALLGather(num_slicesY_mg(n),1,MPI_INTEGER,   &
                         all_slicesY_mg(0,n),1,MPI_INTEGER,flow_comm,ierr)

      if (ierr /= 0) print*,'Warning: MPI_allgather Y says', ierr
      !write(*,'(a,I2,a,I4,a,100I4)')'Grid level=',n,', lProc=',lProc, &
      !        ', all slices:',all_slices_mg(0:nProcs-1,n)
      call MPI_BARRIER(flow_comm,ierr)

       !Z-direction
       !If no coarsening in z, copy all slice info from the finer level
       if( mgrid_K(n) .eq. mgrid_K(n-1)) then
          z_start_mg (n) = z_start_mg (n-1)
          z_end_mg   (n) = z_end_mg   (n-1)
          zc_start_mg(n) = zc_start_mg(n-1)
          zc_end_mg  (n) = zc_end_mg  (n-1)

          zb1_mg     (n) = zb1_mg     (n-1)
          zb2_mg     (n) = zb2_mg     (n-1)

          num_slicesZ_mg(n) = num_slicesZ_mg(n-1)

          b_proc_mg(n)    = b_proc_mg(n-1)
          t_proc_mg(n)    = t_proc_mg(n-1)
          btommost_mg (n) = btommost_mg (n-1)
          topmost_mg(n) = topmost_mg(n-1)

          cycle  ! go to next level
       endif

      !Otherwise, the grid is coarsened. Redo the slice assignment
      if( mod(z_start_mg(n-1), 2)==0 )  then
         z_start_mg(n) = z_start_mg(n-1)/2 + 1
      else
         z_start_mg(n) = (z_start_mg(n-1)+1)/2  
      endif

      if( mod(z_end_mg(n-1), 2)==0 )  then
         z_end_mg(n) = z_end_mg(n-1)/2 + 1
      else
         z_end_mg(n) = (z_end_mg(n-1)+1)/2  
      endif
      num_slicesZ_mg(n) = z_end_mg(n) - z_start_mg(n)

      zc_start_mg(n) = z_start_mg(n)
      zc_end_mg(n)   = z_end_mg(n) - 1

      zb1_mg(n)      = zc_start_mg(n) - 2
      zb2_mg(n)      = zc_end_mg(n) + 2      

      write(*,'(a,I2,a,I4,a,7I4)')'Lev:',n,', id:',lProc,', Z start/end/slices:', &
               zb1_mg(n),z_start_mg(n),zc_start_mg(n),zc_end_mg(n),z_end_mg(n), &
               zb2_mg(n),num_slicesZ_mg(n)

       !Each processor gathers slice info from all others
       call MPI_ALLGather(num_slicesZ_mg(n),1,MPI_INTEGER,   &
                          all_slicesZ_mg(0,n),1,MPI_INTEGER,flow_comm,ierr)

       if (ierr /= 0) print*,'Warning: MPI_allgather Z says', ierr
       !write(*,'(a,I2,a,I4,a,100I4)')'Grid level=',n,', lProc=',lProc, &
       !        ', all slices:',all_slices_mg(0:nProcs-1,n)
       call MPI_BARRIER(flow_comm,ierr)

       !set up left and right neighbors for communication at each level
       ii = jProc - 1
       do while(ii>=0)
          if(all_slicesY_mg(nProcY*kProc+ii,n) == 0) then
            ii = ii - 1
          else
            exit
          endif
       enddo
       l_proc_mg(n) = nProcY*kProc+ii

       ii = jProc + 1
       do while(ii<=nProcY-1)
          if(all_slicesY_mg(nProcY*kProc+ii,n) == 0) then
            ii = ii + 1
          else
            exit
          endif
       enddo
       r_proc_mg(n) = nProcY*kProc+ii    
       
       if(num_slicesY_mg(n)==0) then
          !For empty subdomains, use itself as left/right neighbours
          l_proc_mg(n) = lProc
          r_proc_mg(n) = lProc
       endif

       !Determine the leftmost and rightmost processors, in directional ID
       leftmost_mg( n) = 0        !kProc*nProcY
       rightmost_mg(n) = nProcY-1 !kProc*nProcY+nProcY-1

       if(l_proc_mg(n)-nProcY*kProc == -1  .and. num_slicesY_mg(n)>0) then
          leftmost_mg(n) = 0
          yb1_mg(n) = 0         ! start from 0 for its subdomain; consistent with the original grid
       endif

       if(r_proc_mg(n)-nProcY*kProc == nProcY .and. num_slicesY_mg(n)>0)  &
          rightmost_mg(n) = jProc  !Only THIS processor needs to know if it is at the boundary.
                                   !For the other processors,  rightmost_mg(n)
                                   !is arbitarily set to nProcY-1 initially.

       !write(*,'(a,I2,a,I4,a,2I4,a,2I4)')'Grid level=',n,', lProc=',lProc, &
       !        ' L/R proc =',l_proc_mg(n),r_proc_mg(n),' L/R most=',leftmost_mg(n),rightmost_mg(n)

       !set up bottom and top neighbors for communication at each level
       ii = kProc - 1
       do while(ii>=0)
          if(all_slicesZ_mg(nProcY*ii+jProc,n) == 0) then
            ii = ii - 1
          else
            exit
          endif
       enddo
       b_proc_mg(n) = nProcY*ii+jProc

       ii = kProc + 1
       do while(ii<=nProcZ-1)
          if(all_slicesZ_mg(nProcY*ii+jProc,n) == 0) then
            ii = ii + 1
          else
            exit
          endif
       enddo
       t_proc_mg(n) = nProcY*ii+jProc

       if(num_slicesZ_mg(n)==0) then
          !For empty subdomains, use itself as left/right neighbours
          b_proc_mg(n) = lProc
          t_proc_mg(n) = lProc
       endif
 
       !Determine the btommost and topmost processors, using directional ID
       btommost_mg( n) = 0
       topmost_mg(n)   = nProcZ-1  !nProcY*(nProcZ-1)+jProc

       !If the curent processor is on the bottom row
       !if(b_proc_mg(n)+(kProc+1)*nProcY == lProc .and. num_slicesZ_mg(n)>0) then
       if(b_proc_mg(n)<0 .and. num_slicesZ_mg(n)>0) then
          btommost_mg(n) = 0
          zb1_mg(n) = 0         ! start from 0 for its subdomain; consistent with the original grid
       endif

       if(t_proc_mg(n) == nProcY*nProcZ+jProc .and. num_slicesZ_mg(n)>0)  &
          topmost_mg(n) = kProc

    END DO ! end n

    !Ye,check multi-grid
    !do n=1,mLevel
    !   write(300+lProc,*)'nlev = ', n, num_slicesY_mg(n), num_slicesZ_mg(n)
    !   write(300+lProc,*)'mgrid_I/J/K = ',mgrid_I(n),mgrid_J(n),mgrid_K(n)
    !   write(300+lProc,*)l_proc_mg(n),r_proc_mg(n),b_proc_mg(n),t_proc_mg(n)
    !   write(300+lProc,*)leftmost_mg(n),rightmost_mg(n),btommost_mg(n),topmost_mg(n)
    !   write(300+lProc,*)yc_start_mg(n),yc_end_mg(n),y_start_mg(n),y_end_mg(n),yb1_mg(n),yb2_mg(n)
    !   write(300+lProc,*)zc_start_mg(n),zc_end_mg(n),z_start_mg(n),z_end_mg(n),zb1_mg(n),zb2_mg(n)
    !   write(300+lProc,*)'=========='
    !enddo

    ! Field array allocation for coarse grids
    DO n = 1, mLevel
       !write(*,'(a,5I7)')'MG array alloc:',lProc,n,zb1_mg(n),zb2_mg(n),num_slices_mg(n)

       if(num_slicesY_mg(n)>0 .and. num_slicesZ_mg(n)>0) then
          ! allocate field variable arrays; the start/end indices are consistent with
          ! those on the finest grid, except that in the z-direction, only 1 ghost slice
          ! on each side is added. Note that there are 2 ghost slices on the original grid.
          ALLOCATE( mg_array(n)%rhs(0:mgrid_I(n)+1, yb1_mg(n):yb2_mg(n), zb1_mg(n):zb2_mg(n)) )
          ALLOCATE( mg_array(n)%phi(0:mgrid_I(n)+1, yb1_mg(n):yb2_mg(n), zb1_mg(n):zb2_mg(n)) )
          ALLOCATE( mg_array(n)%res(0:mgrid_I(n)+1, yb1_mg(n):yb2_mg(n), zb1_mg(n):zb2_mg(n)) )
          ALLOCATE( mg_array(n)%ibk(0:mgrid_I(n)+1, yb1_mg(n):yb2_mg(n), zb1_mg(n):zb2_mg(n)) )
          ALLOCATE( mg_array(n)%igk(0:mgrid_I(n)+1, yb1_mg(n):yb2_mg(n), zb1_mg(n):zb2_mg(n)) )
       endif
    END DO ! end n
 
   END SUBROUTINE MG_initial

!------------------------------------------------------------------------------
!
! Purpose: Line SOR with Gauss Siedel method are used as smoother at each level. 
!
! Input: 
!        var  -- initial guesses
!         r   -- values at the right-hand side
!      igmark_mg -- iblank and ghostcell mark at the corrent grid
!       nlev  -- current grid level
!   IL,JL,[zb1g,zb2g] -- grid points in 3 directions (start/end in z)
!       
! Output: var -- storing the approximation of the solution.
!
!------------------------------------------------------------------------------
   SUBROUTINE itsolv_mg(var, r, igmark_mg, nlev, nx_mg, yb1g, yb2g, zb1g, zb2g)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE mg_module
    USE MPI_module

    IMPLICIT NONE

!... parameters
    INTEGER      :: nlev, nx_mg, ny_mg, zb1g, zb2g, yb1g, yb2g
    REAL(KIND=CGREAL), DIMENSION(0:nx_mg+1,yb1g:yb2g,zb1g:zb2g), INTENT (IN)     ::r
    REAL(KIND=CGREAL), DIMENSION(0:nx_mg+1,yb1g:yb2g,zb1g:zb2g), INTENT (INOUT)  ::var
    INTEGER(KIND=INT_K),DIMENSION(0:nx_mg+1,yb1g:yb2g,zb1g:zb2g),INTENT (IN)     ::igmark_mg

!... Local variables
    INTEGER :: i,j,k,ihomo
    INTEGER :: iBody, iRow, iG,jG,kG, iNode, jNode, n,knode
    INTEGER :: iinit, jinit,kinit,Ncolor
    INTEGER :: k0, k1, KL1, KL2, JL1, JL2, j0, j1
    REAL(KIND=CGREAL) :: temp,flag

    !start/end slices in y- and z-direction
    JL1 = yc_start_mg(nlev)
    JL2 = yc_end_mg(nlev)
    KL1 = zc_start_mg(nlev)
    KL2 = zc_end_mg(nlev)

!-----  Line solve in the x-direction -------------------
    DO Ncolor = 1, 1!2, Ye, disable color scheme in y

       iinit = 1;      istep = 1
       !jinit = Ncolor; jstep = 2
       jinit = JL1;    jstep = 1
       kinit = KL1;    kstep = 1

       CALL set_outer_pressure_bc   (var,nx_mg,yb1g,yb2g,zb1g,zb2g,nlev)
       CALL set_outer_ghost_pressure(var,nx_mg,yb1g,yb2g,zb1g,zb2g,nlev)

       !Ye,debug
       !write(*,*)'lProc and nlev = ',nlev, lProc

       DO k = kinit, KL2, kStep
          amz(k) =   mg_array(nlev)%dzcinv(k)  * mg_array(nlev)%dzinv(k)
          apz(k) =   mg_array(nlev)%dzcinv(k+1)* mg_array(nlev)%dzinv(k)
          acz(k) = - ( amz(k) + apz(k) )

          DO j = jinit, JL2, jStep
             amy(j) =   mg_array(nlev)%dycinv(j)  * mg_array(nlev)%dyinv(j)
             apy(j) =   mg_array(nlev)%dycinv(j+1)* mg_array(nlev)%dyinv(j)
             acy(j) = - ( amy(j) + apy(j) )

             DO i=1, nx_mg-1
                amx(i) =   mg_array(nlev)%dxcinv(i)  * mg_array(nlev)%dxinv(i)
                apx(i) =   mg_array(nlev)%dxcinv(i+1)* mg_array(nlev)%dxinv(i)
                acx(i) = - ( amx(i) + apx(i) )

                rhs(i) =    r(i,j,k)  - var(i,j-1,k)*amy(j) &
                                      - var(i,j+1,k)*apy(j) &
                                      - var(i,j,k-1)*amz(k) &
                                      - var(i,j,k+1)*apz(k) 

                flag   = (1 - igmark_mg(i,j,k))  
                amx(i) = amx(i)*flag
                apx(i) = apx(i)*flag
                acx(i) = (acx(i)+acy(j)+acz(k))*flag + (1.0_CGREAL - flag)

                rhs(i) = rhs(i)*flag + (1.0_CGREAL - flag)*var(i,j,k)

                IF(abs(acx(i)).LT.1.0e-12) THEN
                   amx(i) = 0.0_CGREAL
                   apx(i) = 0.0_CGREAL
                   acx(i) = 1.0_CGREAL
                   rhs(i) = var(i,j,k)
                ENDIF
             ENDDO ! i 

             ! Include Derichelt BC implicitly
             !if(nlev>0) then
             !  amx(0)  = 0.0_CGREAL
             !  acx(0)  = 1.0_CGREAL
             !  apx(0)  = 1.0_CGREAL
             !  rhs(0)  = bcxp(0, j,k)*2.0d0

             !  amx(nx_mg) = 1.0_CGREAL
             !  acx(nx_mg) = 1.0_CGREAL
             !  apx(nx_mg) = 0.0_CGREAL
             !  rhs(nx_mg) = bcxp(1,j,k)*2.0d0
             !else

             amx(0)  = 0.0_CGREAL
             apx(0)  = 0.0_CGREAL
             acx(0)  = 1.0_CGREAL
             rhs(0)  = var(0, j,k)

             amx(nx_mg) = 0.0_CGREAL
             apx(nx_mg) = 0.0_CGREAL
             acx(nx_mg) = 1.0_CGREAL
             rhs(nx_mg) = var(nx_mg,j,k)

             !endif

             CALL tdma(amx,acx,apx,rhs,dummy,0,nx_mg)

             !DO i=1,nx_mg-1
             DO i=0,nx_mg
                var(i,j,k) = var(i,j,k)*(1.0_CGREAL-omega_pson) + omega_pson*dummy(i) 
             ENDDO

          ENDDO ! j 
       ENDDO !k
  
    ENDDO ! NcolorX

!    if(nlev==1) then
       !Exchange 2 slices of data and
       !update ghostcell pressure on the finest grid.
!       call send_receive_slices_real(var, 0,nx+1,0,ny+1,zb1,zb2,2) ! exchange 2 slices      
!       call GCM_ghostcell_pressure(var, div)
!    else
       !Exchange 1 slice of buffer between subdomains
       call send_receive_slices_real_y(var, 0,nx_mg+1,yb1g,yb2g,zb1g,zb2g,1)
       call send_receive_slices_real_z(var, 0,nx_mg+1,yb1g,yb2g,zb1g,zb2g,1)
!    endif

!------- Line solve in the y-direction ------------------
    j0 = JL1 - 1
    j1 = JL2 + 1

    DO Ncolor = 1, 2   !color scheme in x

      iinit = Ncolor; istep = 2
      jinit = JL1;    jstep = 1
      kinit = KL1;    kstep = 1

      CALL set_outer_pressure_bc   (var,nx_mg,yb1g,yb2g,zb1g,zb2g,nlev)
      CALL set_outer_ghost_pressure(var,nx_mg,yb1g,yb2g,zb1g,zb2g,nlev)

        DO k= kinit, KL2, kstep
            amz(k) =   mg_array(nlev)%dzcinv(k)  * mg_array(nlev)%dzinv(k)
            apz(k) =   mg_array(nlev)%dzcinv(k+1)* mg_array(nlev)%dzinv(k)
            acz(k) = - ( amz(k) + apz(k) )

            DO i= iinit, nx_mg-1, istep
               amx(i) =   mg_array(nlev)%dxcinv(i)  * mg_array(nlev)%dxinv(i)
               apx(i) =   mg_array(nlev)%dxcinv(i+1)* mg_array(nlev)%dxinv(i)
               acx(i) = - ( amx(i) + apx(i) )

              DO j=jinit, JL2, jstep
                 amy(j) =   mg_array(nlev)%dycinv(j)  * mg_array(nlev)%dyinv(j)
                 apy(j) =   mg_array(nlev)%dycinv(j+1)* mg_array(nlev)%dyinv(j)
                 acy(j) = - ( amy(j) + apy(j) )
      
                 rhs(j) = r(i,j,k) - var(i-1,j,k)*amx(i)  &  
                                   - var(i+1,j,k)*apx(i)  &
                                   - var(i,j,k-1)*amz(k)  &
                                   - var(i,j,k+1)*apz(k) 
 
                 flag   = (1 - igmark_mg(i,j,k))   

                 amy(j) = amy(j)*flag
                 apy(j) = apy(j)*flag
                 acy(j) = (acx(i)+acy(j)+acz(k))*flag + (1.0_CGREAL - flag)
                 rhs(j) = rhs(j)*flag + (1.0_CGREAL - flag)*var(i,j,k)

                 if(abs(acy(j)).LT.1.0e-12) THEN 
                    amy(j) = 0.0_CGREAL 
                    apy(j) = 0.0_CGREAL 
                    acy(j) = 1.0_CGREAL 
                    rhs(j) = var(i,j,k) 
                 end if
              ENDDO

              !if(nlev==0) then
              !   amy(0)  = 0.0_CGREAL
              !   acy(0)  =-1.0_CGREAL
              !   apy(0)  = 1.0_CGREAL
              !   rhs(0)  = 0.0d0

              !   amy(ny_mg) =-1.0_CGREAL
              !   acy(ny_mg) = 1.0_CGREAL
              !   apy(ny_mg) = 0.0_CGREAL
              !   rhs(ny_mg) = 0.0d0
              !else

              amy(j0)  = 0.0_CGREAL
              apy(j0)  = 0.0_CGREAL
              acy(j0)  = 1.0_CGREAL
              rhs(j0)  = var(i,JL1-1,k)
              amy(j1) = 0.0_CGREAL
              apy(j1) = 0.0_CGREAL
              acy(j1) = 1.0_CGREAL
              rhs(j1) = var(i,JL2+1,k)

              !endif

              CALL tdma(amy,acy,apy,rhs,dummy,j0,j1)

              !DO j=1,ny_mg-1
              !DO j=JL1,JL2
              DO j=j0,j1
                 var(i,j,k) = var(i,j,k)*(1.0_CGREAL-omega_pson) + omega_pson*dummy(j) 
              ENDDO
      
           ENDDO
        ENDDO
     ENDDO ! Ncolor

!     if(nlev==1) then
        !Exchange 2 slices of data and
        !update ghostcell pressure on the finest grid.
!        call send_receive_slices_real(var, 0,nx+1,0,ny+1,zb1,zb2,2) ! exchange 2 slices      
!        call GCM_ghostcell_pressure(var, div)
!     else
        !Exchange 1 slice of buffer between subdomains
        call send_receive_slices_real_y(var, 0,nx_mg+1,yb1g,yb2g,zb1g,zb2g,1)
        call send_receive_slices_real_z(var, 0,nx_mg+1,yb1g,yb2g,zb1g,zb2g,1)
!     endif

!------- Line solve in the z-direction ------------------
     k0 = KL1 - 1
     k1 = KL2 + 1

     DO Ncolor = 1, 1!2   !color scheme in y
        CALL set_outer_pressure_bc   (var,nx_mg,yb1g,yb2g,zb1g,zb2g,nlev)
        CALL set_outer_ghost_pressure(var,nx_mg,yb1g,yb2g,zb1g,zb2g,nlev)

        iinit = 1;      istep = 1
        !jinit = Ncolor; jstep = 2
        jinit = JL1;    jstep = 1
        kinit = KL1;    kstep = 1

        DO j= jinit, JL2, jstep
           amy(j) =   mg_array(nlev)%dycinv(j)  * mg_array(nlev)%dyinv(j)
           apy(j) =   mg_array(nlev)%dycinv(j+1)* mg_array(nlev)%dyinv(j)
           acy(j) = - ( amy(j) + apy(j) )
 
           DO i= iinit, nx_mg-1, istep
              amx(i) =   mg_array(nlev)%dxcinv(i)  * mg_array(nlev)%dxinv(i)
              apx(i) =   mg_array(nlev)%dxcinv(i+1)* mg_array(nlev)%dxinv(i)
              acx(i) = - ( amx(i) + apx(i) )
 
              DO k=KL1, KL2

                amz(k) =   mg_array(nlev)%dzcinv(k)  * mg_array(nlev)%dzinv(k)
                apz(k) =   mg_array(nlev)%dzcinv(k+1)* mg_array(nlev)%dzinv(k)
                acz(k) = - ( amz(k) + apz(k) )

                rhs(k) = r(i,j,k)  - var(i-1,j,k)*amx(i)  &  
                                   - var(i+1,j,k)*apx(i)  &
                                   - var(i,j-1,k)*amy(j)  &
                                   - var(i,j+1,k)*apy(j) 

                flag   = (1 - igmark_mg(i,j,k))   

                amz(k) = amz(k)*flag
                apz(k) = apz(k)*flag
                acz(k) = (acx(i)+acy(j)+acz(k))*flag + (1.0_CGREAL - flag)
                rhs(k) = rhs(k)*flag + (1.0_CGREAL - flag)*var(i,j,k)

                if(abs(acz(k)).LT.1.0e-12) THEN 
                   amz(k) = 0.0_CGREAL 
                   apz(k) = 0.0_CGREAL 
                   acz(k) = 1.0_CGREAL 
                   rhs(k) = var(i,j,k) 
                end if
             ENDDO

             !if(nlev==0) then
             !   amz(k0)  = 0.0_CGREAL
             !   acz(k0)  =-1.0_CGREAL
             !   apz(k0)  = 1.0_CGREAL
             !   rhs(k0)  = 0.0d0

             !   amz(k1) =-1.0_CGREAL
             !   acz(k1) = 1.0_CGREAL
             !   apz(k1) = 0.0_CGREAL
             !   rhs(k1) = 0.0d0
             !else
   
             amz(k0)  = 0.0_CGREAL
             apz(k0)  = 0.0_CGREAL
             acz(k0)  = 1.0_CGREAL
             rhs(k0)  = var(i,j,k0)
             amz(k1) = 0.0_CGREAL 
             apz(k1) = 0.0_CGREAL
             acz(k1) = 1.0_CGREAL
             rhs(k1) = var(i,j,k1)

             !endif

             ! Note that these arrays have actual size of 0:nz+1, but only 
             ! zc_start_mg-1:zc_end_mg+1 are being used. It is very important that
             ! the arrays in TDMA have the same starting index, i.e., 0, even
             ! though the ending index is may be shorter.              
             CALL tdma(amz,acz,apz,rhs,dummy,k0,k1)

             !DO k=KL1,KL2
              DO k=k0,k1
                var(i,j,k) = var(i,j,k)*(1.0_CGREAL-omega_pson) + omega_pson*dummy(k) 
             ENDDO
      
          ENDDO
       ENDDO
    ENDDO ! Ncolor

!     if(nlev==1) then
        !Exchange 2 slices of data and
        !update ghostcell pressure on the finest grid.
!        call send_receive_slices_real(var, 0,nx+1,0,ny+1,zb1,zb2,2) ! exchange 2 slices      
!        call GCM_ghostcell_pressure(var, div)
!     else
        !Exchange 1 slice of buffer between subdomains
        call send_receive_slices_real_y(var, 0,nx_mg+1,yb1g,yb2g,zb1g,zb2g,1)
        call send_receive_slices_real_z(var, 0,nx_mg+1,yb1g,yb2g,zb1g,zb2g,1)
!     endif

   END SUBROUTINE itsolv_mg
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! Purpose: Calculate the residues at each level. 
!          Prepare for the RESTRICTION step in MG method
!
! Input: 
!        nx_mg,ny_mg,[zb1g,zb2g] --  the number of grid points in x, y, z direction respectively
!                                    at current level
!        var  -- initial guesses
!        rrr   -- values at the right-hand side
!  
! Output: rrr -- storing the residues at current level
!
!------------------------------------------------------------------------------

   SUBROUTINE residual_mg(var, rrr, igmark_mg, nlev, nx_mg, yb1g, yb2g, zb1g, zb2g)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE GCM_arrays
    USE mg_module
    USE MPI_module

    IMPLICIT NONE

!... parameters
    INTEGER      :: nlev, nx_mg, ny_mg, zb1g, zb2g, yb1g, yb2g
    REAL(KIND=CGREAL), DIMENSION(0:nx_mg+1,yb1g:yb2g,zb1g:zb2g), INTENT (INOUT) ::rrr
    REAL(KIND=CGREAL), DIMENSION(0:nx_mg+1,yb1g:yb2g,zb1g:zb2g), INTENT (IN)    ::var
    INTEGER(KIND=INT_K),DIMENSION(0:nx_mg+1,yb1g:yb2g,zb1g:zb2g), INTENT (IN)   ::igmark_mg

!... local variables
    INTEGER              :: i,j,k, KL1, KL2, JL1, JL2
    REAL(KIND=CGREAL)    :: res, res2, resmax=0.0_CGREAL,flag
    REAL(KIND=CGREAL) :: bmx,bpx,bcx,bc
    REAL(KIND=CGREAL) :: bmy,bpy,bcy
    REAL(KIND=CGREAL) :: bmz,bpz,bcz

!-----------------------------------------------------------------------------------

    !start/end slices in y- and z-direction
    JL1 = yc_start_mg(nlev)
    JL2 = yc_end_mg(nlev)
    KL1 = zc_start_mg(nlev)
    KL2 = zc_end_mg(nlev)
    !print*,'KL1,KL2',KL1,KL2

    res2 = 0.0d0
    DO K=KL1, KL2
      DO J = JL1, JL2
        DO I = 1, nx_mg-1
           bmz =   mg_array(nlev)%dzcinv(k)  * mg_array(nlev)%dzinv(k)
           bpz =   mg_array(nlev)%dzcinv(k+1)* mg_array(nlev)%dzinv(k)
           bcz = - ( bmz + bpz )

           bmy =   mg_array(nlev)%dycinv(j)  * mg_array(nlev)%dyinv(j)
           bpy =   mg_array(nlev)%dycinv(j+1)* mg_array(nlev)%dyinv(j)  
           bcy = - ( bmy + bpy )

           bmx =   mg_array(nlev)%dxcinv(i)  * mg_array(nlev)%dxinv(i)
           bpx =   mg_array(nlev)%dxcinv(i+1)* mg_array(nlev)%dxinv(i) 
           bcx = - ( bmx + bpx )
                   
           flag     = (1.0_CGREAL-REAL(igmark_mg(i,j,k),       KIND=CGREAL))  

           bc  = (bcx+bcy+bcz)

           res    = rrr(i,j,k) - var(i,j,k)*bc              &
                          - var(i-1,j,k)*bmx                &
                          - var(i+1,j,k)*bpx                &
                          - var(i,j-1,k)*bmy                &
                          - var(i,j+1,k)*bpy                &
                          - var(i,j,k-1)*bmz                &
                          - var(i,j,k+1)*bpz                                      

           rrr(i,j,k) = res*flag

           res2 = res2 + res**2*flag
        ENDDO
      ENDDO
    ENDDO

!---------- will be exectued if the user wants to look at the convergence history --------
!---------- should be turned off for normal running case (set infoconv == 0 in input.dat)
    IF (infoconv .EQ. 1) THEN
        !res = MAXVAL(ABS(rrr(1:nx_mg-1,1:ny_mg-1,KL1:KL2)))  ! L-inf norm
        res = sqrt(res2 / dfloat(nx_mg*ny_mg*(KL2-KL1+1)))    ! L-2 norm

        ! Find the max residual among all subdomains. DON't call this
        ! if some subdomains are empty (the corresponding proc won't execute this subroutine)
        !call MPI_ALLREDUCE(res, resmax, 1, MPI_DOUBLE_PRECISION, &
        !                   MPI_MAX, flow_comm, ierr)
        write(*,100) nlev, res  !resmax
    END IF

100 FORMAT('MG: residual check : ',1x,I4,2x,PE12.5)

   END SUBROUTINE residual_mg


!------------------------------------------------------------------------------
!
! Purpose: Inject/extract the residual from a fine grid (level N)
!          to a coarse grid (level N+1)
!
! Input: 
!        nx_mg,ny_mg,[zb1g,zb2g] --  the number of grid points in x, y, z direction respectively
!                                    
!        rr1  -- residual on the fine grid
!        nlev -- level of the fine grid
!  
! Output: rr2 -- residual on the coarse grid
!
!------------------------------------------------------------------------------

  SUBROUTINE mg_inject(N, rr1, nx_mg1, yb1g1, yb2g1, zb1g1, zb2g1,  &
                          rr2, nx_mg2, yb1g2, yb2g2, zb1g2, zb2g2)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE mg_module
    USE MPI_module

    IMPLICIT NONE

!... parameters
    INTEGER      :: N, nx_mg1, ny_mg1, zb1g1, zb2g1, yb1g1, yb2g1
    INTEGER      ::    nx_mg2, ny_mg2, zb1g2, zb2g2, yb1g2, yb2g2

    REAL(KIND=CGREAL), DIMENSION(0:nx_mg1+1,yb1g1:yb2g1,zb1g1:zb2g1), INTENT (IN) ::rr1
    REAL(KIND=CGREAL), DIMENSION(0:nx_mg2+1,yb1g2:yb2g2,zb1g2:zb2g2), INTENT (OUT)::rr2

!...
    INTEGER           :: i, j, k, ii, jj, kk
    INTEGER           :: i1, i2, j1, j2, k1, k2, KL1,KL2, JL1, JL2
    INTEGER           :: i_coarse, j_coarse, k_coarse, ibk_inside
    REAL(KIND=CGREAL) :: alphaX, alphaY, alphaZ
    REAL(KIND=CGREAL) :: xc0,yc0,zc0, r_tmp, r_avg
    REAL(KIND=CGREAL) :: dltX, dltY, dltZ, vol0, vol_i, vol_c, sum0
!--------------------
    !print*,shape(rr1),nx_mg1, ny_mg1, zb1g1, zb2g1
    !print*,shape(rr2),nx_mg2, ny_mg2, zb1g2, zb2g2

    rr2(:,:,:) = 0.0d0

    ! check if the grid is coarsened in each direction
    i_coarse = 1
    j_coarse = 1
    k_coarse = 1
    if(mgrid_I(N) == mgrid_I(N+1)) i_coarse = 0
    if(mgrid_J(N) == mgrid_J(N+1)) j_coarse = 0
    if(mgrid_K(N) == mgrid_K(N+1)) k_coarse = 0
    !print*,N,i_coarse,j_coarse,k_coarse

    JL1 = yc_start_mg(N+1)
    JL2 = yc_end_mg(N+1)
    KL1 = zc_start_mg(N+1)
    KL2 = zc_end_mg(N+1)

    DO k = KL1, KL2
    DO j = JL1, JL2
    DO i = 1, nx_mg2-1

       !location of the node on level N+1
       xc0 = mg_array(N+1)%xc(i)
       yc0 = mg_array(N+1)%yc(j)
       zc0 = mg_array(N+1)%zc(k)

       !surrounding nodes on the level N
       i1 = (1-i_coarse)*i + i_coarse*(2*i - 1)
       i2 = (1-i_coarse)*i + i_coarse*(2*i    )

       j1 = (1-j_coarse)*j + j_coarse*(2*j - 1)
       j2 = (1-j_coarse)*j + j_coarse*(2*j    )

       k1 = (1-k_coarse)*k + k_coarse*(2*k - 1)
       k2 = (1-k_coarse)*k + k_coarse*(2*k    )

       ! trilinear interpolation at the node.
       alphaX = (xc0 - mg_array(N)%xc(i1)) * mg_array(N)%dxcinv(i1+1)
       alphaY = (yc0 - mg_array(N)%yc(j1)) * mg_array(N)%dycinv(j1+1)
       alphaZ = (zc0 - mg_array(N)%zc(k1)) * mg_array(N)%dzcinv(k1+1)

       r_tmp  = rr1(i1,j1,k1)*(1.0_CGREAL - alphaX)*(1.0_CGREAL -alphaY)*(1.0_CGREAL - alphaZ) &
           + rr1(i1+1,j1,k1)*alphaX*(1.0_CGREAL - alphaY)*(1.0_CGREAL - alphaZ)             &
           + rr1(i1,j1+1,k1)*(1.0_CGREAL - alphaX)*alphaY*(1.0_CGREAL - alphaZ)             &
           + rr1(i1,j1,k1+1)*(1.0_CGREAL - alphaX)*(1.0_CGREAL - alphaY)*alphaZ             &
           + rr1(i1+1,j1+1,k1)*alphaX*alphaY*(1.0_CGREAL - alphaZ)             &
           + rr1(i1+1,j1,k1+1)*alphaX*(1.0_CGREAL - alphaY)*alphaZ             &
           + rr1(i1,j1+1,k1+1)*(1.0_CGREAL - alphaX)*alphaY*alphaZ             &
           + rr1(i1+1,j1+1,k1+1)*alphaX*alphaY*alphaZ

       !alternative average based on cell volumes if iblank cells are present
       vol0 = 0.0d0
       sum0 = 0.0d0
       ibk_inside = 0
       do kk=k1,k2  ! note k2,j2,and i2 are the same as k1,j1,i1
       do jj=j1,j2  ! respectively, if there is no coarsening.
       do ii=i1,i2
         !If the stencil involves iblank cells where flow variables are not defined,
         if( mg_array(N)%ibk(ii,jj,kk) == 1) then
            ibk_inside = 1
         else
            dltX  = mg_array(N)%x(ii+1)-mg_array(N)%x(ii)
            dltY  = mg_array(N)%y(jj+1)-mg_array(N)%y(jj)
            dltZ  = mg_array(N)%z(kk+1)-mg_array(N)%z(kk)
            vol_i = dltX * dltY * dltZ

            vol0 = vol0 + vol_i
            sum0 = rr1(ii,jj,kk) * vol_i 
         endif
       enddo
       enddo
       enddo

       vol_c = 1.0d0 / mg_array(N)%dxinv(i1) &
                     / mg_array(N)%dyinv(j1) &
                     / mg_array(N)%dzinv(k1)   !reference volume
       if(vol0 > 0.01 * vol_c) then
          r_avg = sum0 / vol0
       else
          !No cell data was used in the average (all iblank cells)
          r_avg = 0.0d0 
       endif

       !Use either trilinear or area-based interpolations
       rr2(i,j,k) = r_tmp * (1-ibk_inside) + r_avg * ibk_inside
    ENDDO
    ENDDO
    ENDDO
 
    !What about the slices in the buffer zone?

  END SUBROUTINE mg_inject

!------------------------------------------------------------------------------
!
! Purpose: Interpolate the solution from a coarse grid
!          onto a fine grid, as an addition
!
! Input: 
!        nx_mg,ny_mg,[zb1g,zb2g] --  the number of grid points in x, y, z direction respectively
!                                    
!        phi1  -- solution on the coarse grid
!        nlev -- level of the coarse grid
!  
! Output: phi2 -- solution on the fine grid
!         dphi -- solution increment (interpolated from phi1)
!------------------------------------------------------------------------------

  SUBROUTINE mg_interp(N, phi1, nx_mg1,  yb1g1, yb2g1, zb1g1, zb2g1,  &
                          phi2, nx_mg2, ny_mg2,nz_mg2, yb1g2, yb2g2, zb1g2, zb2g2,  &
                          dphi)

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE multiuse_arrays
    USE solver_arrays
    USE mg_module
    USE MPI_module

    IMPLICIT NONE

!... parameters
    INTEGER      :: N, nx_mg1, ny_mg1, zb1g1, zb2g1, yb1g1, yb2g1
    INTEGER      ::    nx_mg2, ny_mg2, zb1g2, zb2g2, yb1g2, yb2g2
    INTEGER      ::    nz_mg2

    REAL(KIND=CGREAL), DIMENSION(0:nx_mg1+1,yb1g1:yb2g1,zb1g1:zb2g1), INTENT (IN) ::phi1
    REAL(KIND=CGREAL), DIMENSION(0:nx_mg2+1,yb1g2:yb2g2,zb1g2:zb2g2), INTENT (INOUT)::phi2
    REAL(KIND=CGREAL), DIMENSION(0:nx_mg2+1,yb1g2:yb2g2,zb1g2:zb2g2), INTENT (INOUT)::dphi

!...
    INTEGER           :: i, j, k, ii, jj, kk
    INTEGER           :: i1, i2, j1, j2, k1, k2, KL1, KL2, JL1, JL2
    INTEGER           :: i_coarse, j_coarse, k_coarse, ibk_inside
    REAL(KIND=CGREAL) :: alphaX, alphaY, alphaZ
    REAL(KIND=CGREAL) :: xc0,yc0,zc0, r_tmp, r_avg
    REAL(KIND=CGREAL) :: dltX, dltY, dltZ, vol0, vol_i, vol_c, sum0
!--------------------

    dphi(:,:,:) = 0.0d0

    ! check if the grid is coarsened in each direction
    i_coarse = 1
    j_coarse = 1
    k_coarse = 1
    if(mgrid_I(N) == mgrid_I(N-1)) i_coarse = 0
    if(mgrid_J(N) == mgrid_J(N-1)) j_coarse = 0
    if(mgrid_K(N) == mgrid_K(N-1)) k_coarse = 0
    !print*,N,i_coarse,j_coarse,k_coarse

    JL1 = yc_start_mg(N-1)
    JL2 = yc_end_mg(N-1)
    KL1 = zc_start_mg(N-1)
    KL2 = zc_end_mg(N-1)

    DO k = KL1, KL2
    DO j = JL1, JL2
    DO i = 1, nx_mg2-1

       ! Skip the edge points which would damage the interpolation
       if((i-1       ) + (j-1       ) == 0 ) cycle
       if((i-1       ) + (ny_mg2-1-j) == 0 ) cycle
       if((nx_mg2-1-i) + (j-1       ) == 0 ) cycle
       if((nx_mg2-1-i) + (ny_mg2-1-j) == 0 ) cycle

       if((i-1       ) + (k-KL1     ) == 0 ) cycle
       if((i-1       ) + (KL2-k     ) == 0 ) cycle
       if((nx_mg2-1-i) + (k-KL1     ) == 0 ) cycle
       if((nx_mg2-1-i) + (KL2-k     ) == 0 ) cycle

       if((j-1       ) + (k-KL1     ) == 0 ) cycle
       if((j-1       ) + (KL2-k     ) == 0 ) cycle
       if((ny_mg2-1-j) + (k-KL1     ) == 0 ) cycle
       if((ny_mg2-1-j) + (KL2-k     ) == 0 ) cycle

       !Ye, btom & top
       if((k-1       ) + (j-JL1     ) == 0 ) cycle
       if((k-1       ) + (JL2-j     ) == 0 ) cycle
       if((nz_mg2-1-k) + (j-JL1     ) == 0 ) cycle
       if((nz_mg2-1-k) + (JL2-j     ) == 0 ) cycle

       if(mg_array(N-1)%ibk(i,j,k) == 1) cycle  ! skip for iblank cell

       !location of the node on level N-1
       xc0 = mg_array(N-1)%xc(i)
       yc0 = mg_array(N-1)%yc(j)
       zc0 = mg_array(N-1)%zc(k)

       !surrounding nodes on the level N
       ii = int(i/2)  ! get the floor integer
       i1 = (1-i_coarse)*i + i_coarse* ii
       i2 = (1-i_coarse)*i + i_coarse*(ii+1)

       jj = int(j/2)
       j1 = (1-j_coarse)*j + j_coarse* jj
       j2 = (1-j_coarse)*j + j_coarse*(jj+1)

       kk = int(k/2)
       k1 = (1-k_coarse)*k + k_coarse* kk
       k2 = (1-k_coarse)*k + k_coarse*(kk+1)

       ! trilinear interpolation at the node.
       alphaX = (xc0 - mg_array(N)%xc(i1)) * mg_array(N)%dxcinv(i1+1)
       alphaY = (yc0 - mg_array(N)%yc(j1)) * mg_array(N)%dycinv(j1+1)
       alphaZ = (zc0 - mg_array(N)%zc(k1)) * mg_array(N)%dzcinv(k1+1)

       r_tmp  = phi1(i1,j1,k1)*(1.0_CGREAL - alphaX)*(1.0_CGREAL -alphaY)*(1.0_CGREAL - alphaZ) &
           + phi1(i1+1,j1,k1)*alphaX*(1.0_CGREAL - alphaY)*(1.0_CGREAL - alphaZ)             &
           + phi1(i1,j1+1,k1)*(1.0_CGREAL - alphaX)*alphaY*(1.0_CGREAL - alphaZ)             &
           + phi1(i1,j1,k1+1)*(1.0_CGREAL - alphaX)*(1.0_CGREAL - alphaY)*alphaZ             &
           + phi1(i1+1,j1+1,k1)*alphaX*alphaY*(1.0_CGREAL - alphaZ)             &
           + phi1(i1+1,j1,k1+1)*alphaX*(1.0_CGREAL - alphaY)*alphaZ             &
           + phi1(i1,j1+1,k1+1)*(1.0_CGREAL - alphaX)*alphaY*alphaZ             &
           + phi1(i1+1,j1+1,k1+1)*alphaX*alphaY*alphaZ

       !alternative average based on cell volumes if iblank cells are present
       vol0 = 0.0d0
       sum0 = 0.0d0
       ibk_inside = 0
       do kk=k1,k2  ! note k2,j2,and i2 are the same as k1,j1,i1
       do jj=j1,j2  ! respectively, if there is no coarsening.
       do ii=i1,i2
         !skip iblank cells,
         if( mg_array(N)%ibk(ii,jj,kk) == 1) then
            ibk_inside = 1
         else
            dltX  = mg_array(N)%x(ii+1)-mg_array(N)%x(ii)
            dltY  = mg_array(N)%y(jj+1)-mg_array(N)%y(jj)
            dltZ  = mg_array(N)%z(kk+1)-mg_array(N)%z(kk)
            vol_i = dltX * dltY * dltZ

            vol0 = vol0 + vol_i
            sum0 = phi1(ii,jj,kk) * vol_i
         endif
       enddo
       enddo
       enddo

       vol_c = 1.0d0 / mg_array(N)%dxinv(i1) &
                     / mg_array(N)%dyinv(j1) &
                     / mg_array(N)%dzinv(k1)   !reference volume
       if(vol0 > 0.01 * vol_c) then
          r_avg = sum0 / vol0
       else
          !No cell data was used in the average (all iblank cells)
          r_avg = 0.0d0
       endif

       !Use either trilinear or area-based interpolations
       dphi(i,j,k) = r_tmp * (1-ibk_inside) + r_avg * ibk_inside
       phi2(i,j,k) = phi2(i,j,k) + dphi(i,j,k)
    ENDDO
    ENDDO
    ENDDO

    !What about the slices in the buffer zone?

  END SUBROUTINE mg_interp

!------------------------------------------------------------------------------
!
! Purpose: Prepare MG solver for each time step: set igmark, etc.
!
! Input:
!         var  -- input pressure values
!         r   -- values at the right-hand side
!
! Output: var -- storing the approximation of the solution.
!
!------------------------------------------------------------------------------

    SUBROUTINE mg_prepare()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE mg_module
    USE MPI_module
 
    IMPLICIT NONE
    !REAL(KIND=CGREAL),DIMENSION(0:nx+1,0:ny+1,zb1:zb2),INTENT (IN) :: rr 
    !REAL(KIND=CGREAL),DIMENSION(0:nx+1,0:ny+1,zb1:zb2),INTENT (IN) :: var

    INTEGER           :: i, j, k, ii, jj, kk
    INTEGER           :: i1, i2, j1, j2, k1, k2
    INTEGER           :: N, nx_mg, ny_mg, KL1, KL2, zb1g, zb2g
    INTEGER           :: i_coarse, j_coarse, k_coarse, sumIblank_mg
    REAL(KIND=CGREAL) :: dltX, dltY, dltZ, vol0, vol_ibk, vol_i

    INTEGER           :: JL1, JL2, yb1g, yb2g

    igmark(:,:,:) = 0
    DO k=zb1,zb2
    DO j=yb1,yb2
    DO i=0,nx
       if(iblank(i,j,k)==1 .or. ghostcellmark(i,j,k)==1)  igmark(i,j,k) = 1
    ENDDO
    ENDDO
    ENDDO

!-----------  
! preparing for igmark_mg for each level 
!-----------

    !Populate the arrays at the finest level

    !zb1g = zb1_mg(1)
    !zb2g = zb2_mg(1)
    !write(*,'(a,10I5)')'Shape:',lProc,shape(mg_array(1)%igk), &
    !        shape(mg_array(1)%rhs),shape(mg_array(1)%phi)
    DO k = zb1,zb2
    DO j = yb1,yb2
    DO i = 0, nx
       mg_array(1)%ibk(i,j,k) = iblank(i,j,k)
       mg_array(1)%igk(i,j,k) = igmark(i,j,k)
       !mg_array(1)%rhs(i,j,k) = rr    (i,j,k)
       !mg_array(1)%phi(i,j,k) = var   (i,j,k)
    ENDDO
    ENDDO
    ENDDO

    DO N = 2, mLevel

       if((num_slicesY_mg(N)==0).or.(num_slicesZ_mg(N)==0) )cycle  ! Skip if this is an empty subdomain

       nx_mg = mgrid_I(N)

       JL1 = yc_start_mg(N)
       JL2 = yc_end_mg(N)

       yb1g= yb1_mg(N)
       yb2g= yb2_mg(N)

       KL1 = zc_start_mg(N)
       KL2 = zc_end_mg(N)

       zb1g= zb1_mg(N)
       zb2g= zb2_mg(N)
       !write(*,'(a,I4,a,I5,2I5)')'Lev:',N,', id:',lProc,zb1g,zb2g

       lProc_leftmost  = leftmost_mg (N)
       lProc_rightmost = rightmost_mg(N)
       lProc_left      = l_proc_mg   (N)
       lProc_right     = r_proc_mg   (N)

       lProc_btommost  = btommost_mg (N)
       lProc_topmost   = topmost_mg(N)
       lProc_bottom    = b_proc_mg   (N)
       lProc_top       = t_proc_mg   (N)

       ! check if the grid is coarsened in each direction
       i_coarse = 1
       j_coarse = 1
       k_coarse = 1
       if(mgrid_I(N) == mgrid_I(N-1)) i_coarse = 0
       if(mgrid_J(N) == mgrid_J(N-1)) j_coarse = 0
       if(mgrid_K(N) == mgrid_K(N-1)) k_coarse = 0
       !write(*,'(a,I4,a,3I4)')'Lev:',N,', ijk coarse:',i_coarse,j_coarse,k_coarse

       mg_array(N)%igk(:, :, :) = 0   ! initialize
       mg_array(N)%ibk(:, :, :) = 0
       sumIblank_mg = 0

       DO k = KL1, KL2
       DO j = JL1, JL2
       DO i = 1, nx_mg-1

          i1 = (1-i_coarse)*i + i_coarse*(2*i - 1)
          i2 = (1-i_coarse)*i + i_coarse*(2*i    )

          j1 = (1-j_coarse)*j + j_coarse*(2*j - 1)
          j2 = (1-j_coarse)*j + j_coarse*(2*j    )

          k1 = (1-k_coarse)*k + k_coarse*(2*k - 1)
          k2 = (1-k_coarse)*k + k_coarse*(2*k    )

          vol0    = 0.0d0
          vol_ibk = 0.0d0
          loopK: do kk=k1,k2
          do jj=j1,j2
          do ii=i1,i2
             
             !if any one fine cell is iblank or ghost cell, the entire coarse cell 
             ! is marked
             if(mg_array(N-1)%igk(ii,jj,kk) == 1) then
                mg_array(N)%igk(i, j, k) = 1
                !exit loopK
             endif

             !Determine iblank based on the volume fraction
             dltX  = mg_array(N-1)%x(ii+1)-mg_array(N-1)%x(ii)
             dltY  = mg_array(N-1)%y(jj+1)-mg_array(N-1)%y(jj)
             dltZ  = mg_array(N-1)%z(kk+1)-mg_array(N-1)%z(kk)
             vol_i = dltX * dltY * dltZ 
             vol0  = vol0 + vol_i
             if(mg_array(N-1)%ibk(ii,jj,kk) == 1) then
                vol_ibk = vol_ibk + vol_i
             endif
          enddo
          enddo           
          enddo loopK

          !Set as iblank cell if there is more iblank volume
          !Otherwise set as non-iblank cell, even though it could be an igmark cell.
          if (vol_ibk / vol0 >= 0.5d0) then
             mg_array(N)%ibk(i,j,k) = 1
             sumIblank_mg = sumIblank_mg + 1
          endif
       ENDDO
       ENDDO
       ENDDO
       !IF ( MOD(ntime,nmonitor) == 0) &
       !write(*,'(a,I5,a,I4,a,I8)')'lProc',lProc,', N lev =',N,', sumIblank = ',sumIblank_mg

       !Exchange iblank and igmark between subdomains
       call send_receive_slices_integer_y(mg_array(N)%ibk,0,nx_mg+1,yb1g,yb2g,zb1g,zb2g,1)
       call send_receive_slices_integer_z(mg_array(N)%ibk,0,nx_mg+1,yb1g,yb2g,zb1g,zb2g,1)
    
       call send_receive_slices_integer_y(mg_array(N)%igk,0,nx_mg+1,yb1g,yb2g,zb1g,zb2g,1)
       call send_receive_slices_integer_z(mg_array(N)%igk,0,nx_mg+1,yb1g,yb2g,zb1g,zb2g,1)

!----
!DEBUG
!----
!    IF(N==6) THEN
!     OPEN(UNIT=900+lProc,STATUS='UNKNOWN')
    
!      write(900+lProc,*)'VARIABLES="X","Y","Z","IGK"'
!      write(900+lProc,*)'ZONE F=POINT, I=',nx_mg-1,', J=',num_slicesY_mg(N),' K=',num_slicesZ_mg(N)
!      do k=KL1, KL2
!      do j=JL1, JL2
!      do i=1,nx_mg-1
!         write(900+lProc,123)mg_array(N)%xc(i),mg_array(N)%yc(j),mg_array(N)%zc(k), &
!                             mg_array(N)%igk(i,j,k)
!      enddo
!      enddo
!      enddo
!      close(900+lProc)
!123   format(3(2X,F14.7),3(2X,I4))
!    ENDIF

       
       !initialize the arrays
       mg_array(N)%rhs(:,:,:)    = 0.0_CGREAL
       mg_array(N)%phi(:,:,:)    = 0.0_CGREAL
 
    END DO ! end for N=2, mLevel      

    !reset left- and right lProcs
    lProc_leftmost  = leftmost_mg (1)
    lProc_rightmost = rightmost_mg(1)
    lProc_btommost  = btommost_mg(1)
    lProc_topmost   = topmost_mg(1)

    lProc_left  =  kProc*nProcY+(jProc-1)
    lProc_right =  kProc*nProcY+(jProc+1)
    lProc_top   = (kProc+1)*nProcY+jProc
    lProc_bottom= (kProc-1)*nProcY+jProc

  END SUBROUTINE mg_prepare

!------------------------------------------------------------------------------
!========================================

!------------------------------------------------------------------------------
!
! Purpose: apply MG method with coarsening for all 3 directions
!
! Input:
!         var  -- input pressure values
!         r   -- values at the right-hand side
!
! Output: var -- storing the approximation of the solution.
!
!------------------------------------------------------------------------------

    SUBROUTINE mg_solve(var, rr)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
!    USE GCM_arrays
    USE mg_module
    USE MPI_module
 
    IMPLICIT NONE

    REAL(KIND=CGREAL),DIMENSION(0:nx+1,yb1:yb2,zb1:zb2),INTENT (IN) :: rr 
    REAL(KIND=CGREAL),DIMENSION(0:nx+1,yb1:yb2,zb1:zb2),INTENT (INOUT) :: var

    INTEGER           :: i, j, k, ii, jj, kk
    INTEGER           :: N, nx_mg, ny_mg, KL1, KL2, zb1g, zb2g
    INTEGER           :: i_coarse, j_coarse, k_coarse
    INTEGER           :: nx_mg2, ny_mg2, nz_mg2, zb1g2, zb2g2
    INTEGER           :: iter_mg, loc(3)
    REAL(KIND=CGREAL) :: resm, resm_i
   
    INTEGER           :: JL1, JL2, yb1g, yb2g, yb1g2, yb2g2

!---------
    !Initialize the solution on finest grid
    DO k = zb1, zb2
    DO j = yb1, yb2
    DO i = 0, nx
       mg_array(1)%phi(i,j,k) = var(i,j,k)
       mg_array(1)%rhs(i,j,k) = rr (i,j,k)
    ENDDO
    ENDDO
    ENDDO

!---
! V-loop iterations
!---
   iter_mg = 0
   resm = 1.0E6
   DO WHILE ((iter_mg<=iterInter) .and. (resm .GT. restol_pson))

!------
! From fine levels to coarse levels
!------
    DO N = 1, mLevel
       !Test the MG solver on the finest grid
       !CALL itsolv_mg(var, rr, igmark, 1, nx, ny, zb1, zb2)
       !CALL residual_mg(var, nlu, igmark, 1, nx, ny, zb1, zb2)

       if( (num_slicesY_mg(N)==0).or.(num_slicesZ_mg(N)==0) )  cycle  ! skip if subdomain is empty at this level

       !Leftmost and rightmost processors on each grid.
       ! For lProc not left- or right-most, the values were set to -1 or nProcs.
       lProc_leftmost  = leftmost_mg (N)
       lProc_rightmost = rightmost_mg(N)
       lProc_left      = l_proc_mg   (N)
       lProc_right     = r_proc_mg   (N)

       lProc_btommost = btommost_mg(N)
       lProc_topmost  =  topmost_mg(N)
       lProc_bottom   =   b_proc_mg(N)
       lProc_top      =   t_proc_mg(N)
       !write(*,'(a,I3,a,I4,a,I5,a,2I4,a,2I4)') &
       !     'Lev: ',N,', id:',lProc, ', slices:', num_slices_mg(N), &
       !     ', L/R: ',lProc_left,lProc_right,   &
       !     ', L/R most:',lProc_leftmost,lProc_rightmost

       nx_mg = mgrid_I(N)
       !ny_mg = mgrid_J(N)

       !Y-direction
       JL1 = yc_start_mg(N)
       JL2 = yc_end_mg(N)

       yb1g = yb1_mg(N)
       yb2g = yb2_mg(N)

       !Z-direction
       KL1 = zc_start_mg(N)
       KL2 = zc_end_mg(N)

       zb1g= zb1_mg(N)
       zb2g= zb2_mg(N)

       !Initialize solution on the coarse grid
       if(N>1)   mg_array(N)%phi(:,:,:)    = 0.0_CGREAL

       !write(*,'(a,I4,a,I3,a,5I5)')'id:',lProc,' solving on lev:',N,', grid:',nx_mg,JL1,JL2,KL1,KL2
       CALL itsolv_mg(mg_array(N)%phi, mg_array(N)%rhs, mg_array(N)%igk, &
                      N, nx_mg, yb1g, yb2g, zb1g, zb2g)

       !Update ghostcell pressure on the finest grid.
       !if(N==1) then
       !   call send_receive_slices_real(mg_array(N)%phi, 0,nx+1,0,ny+1,zb1,zb2,2) ! exchange 2 slices      
       !   CALL GCM_ghostcell_pressure(mg_array(N)%phi,div)
       !endif

       !Copy over the r.h.s. (the rhs array should remain unchanged!!!)
       DO k = KL1, KL2
       DO j = JL1, JL2
       DO i = 1, nx_mg
          mg_array(N)%res(i,j,k) = mg_array(N)%rhs(i,j,k)
       ENDDO
       ENDDO
       ENDDO

       ! Calculate residual on the current grid; exchange among subdomains.
       CALL residual_mg(mg_array(N)%phi, mg_array(N)%res, mg_array(N)%igk,  &
                        N, nx_mg, yb1g, yb2g, zb1g, zb2g)
       call send_receive_slices_real_y(mg_array(N)%res,0,nx_mg+1,yb1g,yb2g,zb1g,zb2g,1)
       call send_receive_slices_real_z(mg_array(N)%res,0,nx_mg+1,yb1g,yb2g,zb1g,zb2g,1)

       !Restrict the residual to the coarser, non-empty, grid
       IF(N /= mLevel) THEN
          nx_mg2 = mgrid_I(N+1)
          !ny_mg2 = mgrid_J(N+1)

          yb1g2 = yb1_mg(N+1)
          yb2g2 = yb2_mg(N+1)

          zb1g2 = zb1_mg(N+1)
          zb2g2 = zb2_mg(N+1)          

          IF( num_slicesZ_mg(N+1)>0 .and.num_slicesY_mg(N+1)>0 ) then
             !Note that data is restricted from res at N to rhs at N+1
             CALL mg_inject(N, mg_array(N)%res,   nx_mg,  yb1g, yb2g,  zb1g,  zb2g,  &
                              mg_array(N+1)%rhs, nx_mg2, yb1g2, yb2g2, zb1g2, zb2g2)
          ENDIF

       ENDIF
    ENDDO ! enddo N

    !Synchronize the processors (some have empty subdomain 
    !and thus are much faster)
    !call MPI_BARRIER(flow_comm,ierr)
    !if (lProc.eq.proc_m)then
    !   write(*,*)'MG:from Fine to Coarse...'
    !endif
  
!------
! From coarse levels to fine levels
!------
    DO N = mLevel, 2, -1

       if( num_slicesY_mg(N)==0 .or. num_slicesZ_mg(N)==0) cycle  ! skip if subdomain is empty at this level

       !Leftmost and rightmost processors on each grid.
       ! For lProc not left- or right-most, the values were set to -1 or nProcs.
       lProc_leftmost  = leftmost_mg (N)
       lProc_rightmost = rightmost_mg(N)
       lProc_left      = l_proc_mg   (N)
       lProc_right     = r_proc_mg   (N)
       !write(*,'(a,4I5)')'Left-/right procs',N,lProc,lProc_leftmost,lProc_rightmost

       lProc_btommost = btommost_mg(N)
       lProc_topmost  = topmost_mg(N)
       lProc_bottom   = b_proc_mg(N)
       lProc_top      = t_proc_mg(N)
 
       nx_mg = mgrid_I(N)
       !ny_mg = mgrid_J(N)

       yb1g = yb1_mg(N)
       yb2g = yb2_mg(N)

       zb1g= zb1_mg(N)
       zb2g= zb2_mg(N)

       IF(N /= mLevel) THEN  ! if not the coarsest level
          !write(*,'(a,I7,a,I3,a,5I5)')'id:',lProc,' solving on lev:',N,', grid:',nx_mg,YL1,YL2,KL1,KL2
          CALL itsolv_mg(mg_array(N)%phi, mg_array(N)%rhs, mg_array(N)%igk, &
                         N, nx_mg, yb1g, yb2g, zb1g, zb2g)
       ENDIF 

       IF(N > 1) THEN  ! if not the finest level

          nx_mg2 = mgrid_I(N-1)
          ny_mg2 = mgrid_J(N-1)
          nz_mg2 = mgrid_K(N-1)

          yb1g2 = yb1_mg(N-1)
          yb2g2 = yb2_mg(N-1)

          zb1g2 = zb1_mg(N-1)
          zb2g2 = zb2_mg(N-1)          
          
          !write(*,'(a,I3,a,4I5)')'Interp to lev ',N-1,' grid:',nx_mg2,ny_mg2,zb1g2,zb2g2
          CALL mg_interp(N, mg_array(N)%phi,   nx_mg,  yb1g,  yb2g,  zb1g,  zb2g,  &
                           mg_array(N-1)%phi, nx_mg2, ny_mg2,nz_mg2, yb1g2, yb2g2, zb1g2, zb2g2, &
                           mg_array(N-1)%res)
       ENDIF

    ENDDO ! enddo N

!DEBUG
!!$if(iter_mg==iterInter-1) then
!!$  DO N=1,3
!!$    if(lProc==lProc_leftmost) then
!!$      KL1 = zc_start_mg(N) - 1
!!$    else
!!$      KL1 = zc_start_mg(N) 
!!$    endif
!!$
!!$    if(lProc==lProc_rightmost) then
!!$      KL2 = zc_end_mg(N) + 1
!!$    else
!!$      KL2 = zc_end_mg(N) 
!!$    endif
!!$    nx_mg = mgrid_I(N)
!!$    ny_mg = mgrid_J(N)
!!$
!!$    write(600+lProc+10*N,*)'VARIABLES="X","Y","Z","phi","res"'
!!$    write(600+lProc+10*N,*)'ZONE F=POINT, I=',nx_mg+1,', J=',ny_mg+1,' ,K=',KL2-KL1+1
!!$    do k=KL1,KL2
!!$    do j=0,ny_mg
!!$    do i=0,nx_mg
!!$       write(600+lProc+10*N,'(6(2X,F12.5))')mg_array(N)%xc(i),mg_array(N)%yc(j),mg_array(N)%zc(k), &
!!$                                            mg_array(N)%phi(i,j,k), mg_array(N)%rhs(i,j,k)
!!$    enddo
!!$    enddo
!!$    enddo
!!$    close(600+lProc+10*N)
!!$  ENDDO ! enddo N
!!$endif
!DEBUG

1000 continue

    !CALL set_outer_pressure_bc(mg_array(1)%phi,nx,ny,zb1,zb2,1)   ! inhomogeneous B.C.
    !CALL set_outer_ghost_pressure(mg_array(1)%phi,nx,ny,zb1,zb2,1)

    call send_receive_slices_real_y(mg_array(1)%phi, 0,nx+1,yb1,yb2,zb1,zb2,2)
    call send_receive_slices_real_z(mg_array(1)%phi, 0,nx+1,yb1,yb2,zb1,zb2,2)

    !Calculate the residual on the finest grid
    call calc_residual(mg_array(1)%phi, mg_array(1)%rhs, resm,loc)

    IF (infoconv .EQ. 1) THEN
       write(*,110) iter_mg, resm
    END IF

110 FORMAT('MG: residual check at lev 1: ',1x,I4,1(2X,1PE12.5))

    iter_mg = iter_mg + 1

    IF (mod(ntime,nmonitor)==0 .and. lProc==PROC_M .and. (resm.le.restol_pson)) THEN
       write(*,'(a,1x,I4,1(2X,1PE12.5))')'MG V-loop residual:',iter_mg, resm
    END IF

!------
     ENDDO  ! Finishing V-loop
!------


    !Return the solution at the finest level
    DO k = zb1, zb2
    DO j = yb1, yb2
    DO i = 0, nx
       var   (i,j,k) = mg_array(1)%phi(i,j,k)
    ENDDO
    ENDDO
    ENDDO

    !reset left- and right lProcs
    lProc_leftmost  = leftmost_mg (1)
    lProc_rightmost = rightmost_mg(1)
    lProc_btommost  = btommost_mg(1)
    lProc_topmost   = topmost_mg(1)

    lProc_left  =  kProc*nProcY+(jProc-1)
    lProc_right =  kProc*nProcY+(jProc+1)
    lProc_top   = (kProc+1)*nProcY+jProc
    lProc_bottom= (kProc-1)*nProcY+jProc

  END SUBROUTINE mg_solve

!------------------------------------------------------------------------------
!========================================
!
! Issues to think about: 
! 1) outer BC at coarse levels
! 2) interp. and inject involving iblank or igmark cells
! 3) ... involving outer cells
