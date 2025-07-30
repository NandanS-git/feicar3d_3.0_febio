!------------------------------------------    
!  SUBROUTINE write_monitor() 
!  SUBROUTINE write_dump()
!  SUBROUTINE write_minmaxvals
!  SUBROUTINE vorticity()
!  SUBROUTINE divergence()
!------------------------------------------    

!------------------------------------------    
    SUBROUTINE write_monitor() 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE MPI_module

    IMPLICIT NONE
    INTEGER                 :: i,j,k
    INTEGER, DIMENSION(3)   :: maxcflloc,minddfloc
    REAL(KIND=CGREAL)       :: maxcfl,minddf,rnDim,cfl_i

! cfl = dt*{ |u|/dx  + |v|/dy + |w|/dz }
! diagonal dominance factor of advection diffusion equation
!  = (1 + rx + ry + rz)/(rx + ry + rz)

    nlv = 0.0_CGREAL
    nlw = 0.0_CGREAL

    rnDim = REAL((nDim-DIM_2D),KIND=CGREAL)

    DO k = zc_start,zc_end  !1,nz-1
    DO j = yc_start,yc_end
    DO i = 1,nx-1
      nlv(i,j,k) = dt*( abs(u(i,j,k))*dxinv(i)   &
                       +abs(v(i,j,k))*dyinv(j)   &
                       +abs(w(i,j,k))*dzinv(k) ) &
                 *( 1.0_CGREAL-(REAL(iblank(i,j,k),KIND=CGREAL)) )

      nlw(i,j,k) =  dt*reinv* ( dxinv(i)**2       &
                               +dyinv(j)**2       &
                               +rnDim*dzinv(k)**2 )
      nlw(i,j,k) = 1.0_CGREAL + 1.0_CGREAL/nlw(i,j,k)

    ENDDO
    ENDDO
    ENDDO

    maxcfl    = MAXVAL(nlv(1:nx-1,yc_start:yc_end,zc_start:zc_end))
    maxcflloc = MAXLOC(nlv(1:nx-1,yc_start:yc_end,zc_start:zc_end))

    cfl_i = maxcfl
    call MPI_ALLREDUCE(cfl_i, maxcfl, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, flow_comm, ierr)

    if(lProc==PROC_M) WRITE(STDOUT,'(A,E15.7,A,3(2X,I4))')' Max CFL              = ',maxcfl

    minddf    = MINVAL(nlw(1:nx-1,yc_start:yc_end,zc_start:zc_end))
    minddfloc = MINLOC(nlw(1:nx-1,yc_start:yc_end,zc_start:zc_end))

    !WRITE(STDOUT,'(A,E15.7,A,3(2X,I4))')' Min Diagonal Dom. Fac AD = ',minddf,' at ',minddfloc
    !WRITE(998,'(F15.6,1X,F12.5,3(2X,I4), F12.5)') time, maxcfl, maxcflloc, minddf

    CALL divergence()

    CALL write_minmaxvals()

    END SUBROUTINE write_monitor 
!------------------------------------------
!    
!------------------------------------------    
    SUBROUTINE write_minmaxvals() 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE MPI_module
    
    IMPLICIT NONE
    INTEGER                 :: i,j,k
    REAL(KIND=CGREAL)       :: maxp
    INTEGER                 :: loc(3)
    REAL(KIND=CGREAL)       :: uMin,   uMax,   vMin,   vMax,   wMin,   wMax,   pMin,   pMax
    REAL(KIND=CGREAL)       :: uMin_i, uMax_i, vMin_i, vMax_i, wMin_i, wMax_i, pMin_i, pMax_i

    uMin = minval(u(0:nx,yc_start:yc_end,zc_start:zc_end))
    uMax = maxval(u(0:nx,yc_start:yc_end,zc_start:zc_end))
    vMin = minval(v(0:nx,yc_start:yc_end,zc_start:zc_end))
    vMax = maxval(v(0:nx,yc_start:yc_end,zc_start:zc_end))
    wMin = minval(w(0:nx,yc_start:yc_end,zc_start:zc_end))
    wMax = maxval(w(0:nx,yc_start:yc_end,zc_start:zc_end))
    pMin = minval(p(0:nx,yc_start:yc_end,zc_start:zc_end))
    pMax = maxval(p(0:nx,yc_start:yc_end,zc_start:zc_end))

    uMin_i = uMin 
    call MPI_ALLREDUCE(uMin_i, uMin, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MIN, flow_comm, ierr)
    uMax_i = uMax
    call MPI_ALLREDUCE(uMax_i, uMax, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, flow_comm, ierr)
    vMin_i = vMin
    call MPI_ALLREDUCE(vMin_i, vMin, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MIN, flow_comm, ierr)
    vMax_i = vMax  
    call MPI_ALLREDUCE(vMax_i, vMax, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, flow_comm, ierr)
    wMin_i = wMin
    call MPI_ALLREDUCE(wMin_i, wMin, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MIN, flow_comm, ierr)
    wMax_i = wMax  
    call MPI_ALLREDUCE(wMax_i, wMax, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, flow_comm, ierr)
    pMin_i = pMin
    call MPI_ALLREDUCE(pMin_i, pMin, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MIN, flow_comm, ierr)
    pMax_i = pMax  
    call MPI_ALLREDUCE(pMax_i, pMax, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, flow_comm, ierr)

    if(lProc==PROC_M) then 
        write(6,'(a,I4,a,2E15.5)')'lProc:',lProc,' Min-Max of U(0: =', uMin, uMax
        write(6,'(a,I4,a,2E15.5)')'lProc:',lProc,' Min-Max of V(0: =', vMin, vMax
        write(6,'(a,I4,a,2E15.5)')'lProc:',lProc,' Min-Max of W(0: =', wMin, wMax
        write(6,'(a,I4,a,2E15.5)')'lProc:',lProc,' Min-Max of P(0: =', pMin, pMax
     endif

     !Ye, find the location of max Pressure
     !do i=0, nx
     !   do j=0, ny
     !      do k=zc_start, zc_end
     !         if (abs(p(i,j,k)-pMax).lt.1.0E-8) then
     !            write(*,*)'pMax at:', lProc, P(i,j,k), i, j, k
     !         endif
     !      enddo
     !   enddo
     !enddo

    END SUBROUTINE write_minmaxvals 
!------------------------------------------
!
!------------------------------------------    
    SUBROUTINE write_dump() 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE MPI_module

    IMPLICIT NONE

    !PRINT*,'Writing out dump files'
    CALL write_flowfield()
    if(lProc == PROC_M) CALL write_markers()

    END SUBROUTINE write_dump
!------------------------------------------
!
!------------------------------------------    
    SUBROUTINE write_flowfield() 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE MPI_module

    IMPLICIT NONE
   
    REAL                      :: fsmach,alp,density,xhat,yhat
    INTEGER                   :: i,j,k,iBody,n,m,nCylinder
    CHARACTER*20              :: fname1, fname2
    CHARACTER*100             :: command


    !PRINT*,'Writing out flowfield file'

    write(fname1,'(I7.7)')  ntime
    write(fname2,'(I4.4)')  lProc
    fname1 = 'q_'//trim(fname1)//'.'//fname2
    print*,'Write to file: ',fname1

    OPEN(UNIT=70,FILE=fname1,STATUS='UNKNOWN')

    call vorticity()

    !-------------------------------
    IF (format_dump == TECPLOT  ) THEN
    !-------------------------------
    IF (nDim == DIM_2D) THEN
      write(70,*)'VARIABLES="X","Y","U","V","P","OZ","IBLANK","GHOST","DEAD"'
      write(70,*)'ZONE F=POINT, I=',nx-1,', J=',ny-1
      k = 1
      do j=1,ny-1
      do i=1,nx-1
         write(70,122)xc(i),yc(j),u(i,j,k),v(i,j,k),p(i,j,k)       &
                     ,nlw(i,j,k), iblank(i,j,k),ghostcellmark(i,j,k) &
                     ,dead_cell(i,j,k)
      enddo
      enddo
    ELSE
      write(70,'(a)')'VARIABLES="X","Y","Z","U","V","W","P","OX","OY","OZ","LMD","IBLANK","GHOST","DEAD","B#","goa"'
      write(70,*)'ZONE F=POINT, I=',nx+1,', J=',jSlices+2,' K=',kSlices+2
      do k=zc_start-1,zc_end+1
      do j=yc_start-1,yc_end+1
      do i=0,nx
        write(70,123)xc(i),yc(j),zc(k) &
                    ,u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k)   &
                    ,nlu(i,j,k),nlv(i,j,k),nlw(i,j,k),lmd(i,j,k)             &
                    ,iblank(i,j,k),ghostcellmark(i,j,k),dead_cell(i,j,k)     &
                    ,bodyNum(i,j,k),goa(i,j,k)
!       write(70,124)u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k)   &
!                   ,iblank(i,j,k),ghostcellmark(i,j,k),dead_cell(i,j,k)

      enddo
      enddo
      enddo
    ENDIF
    CLOSE(70)

    !System call to compress the data file
    command = 'gzip '//trim(fname1)    
    call system(command)
    print*, command

127 FORMAT(7(2x,e14.7),2(2x,i4))
      
    !--------------
    ENDIF ! tecplot
    !--------------

122 format(6(2x, pe12.5),3(2x,i2)) 
123 format(11(2x,pe12.5),5(2x,i2))
124 format( 4(2x,pe12.5),3(2x,i2))

    END SUBROUTINE write_flowfield
!------------------------------------------
!
!------------------------------------------    
    SUBROUTINE write_markers() 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE fsi_module

    IMPLICIT NONE
   
    REAL                      :: fsmach,alp,density,xhat,yhat
    INTEGER                   :: i,j,k,iBody,n,m,nCylinder
    CHARACTER*13              :: fname1
    CHARACTER*20              :: fname2

    !print*,'Write out marker file...'

    !   Marker points viz file
    
    write(fname2,'(I7.7)')  ntime
    fname2 = 'marker_'//trim(fname2)
    print*,'Write to file: ',fname2

    !IF ( boundary_motion == MOVING_BOUNDARY ) THEN

      OPEN(UNIT=234,FILE=fname2,STATUS='UNKNOWN')

      WRITE(234,*)'TITLE="3D TRIANGULAR SURFACE DATA"'
      WRITE(234,'(a)')'VARIABLES="X","Y","Z","u","v","w","sx","sy","sz","fx","fy","fz",&
                             "u1","v1","w1", &
                             "cx","cy","cz","ct","d12","d21","d13","d31","d23","d32", &
                             "en12","en21","en13","en31","en23","en32","cfg12","cfg21", &
                             "cfg13","cfg31","cfg23","cfg32"'
 
      DO iBody = 1, nBody
         WRITE(234,*)'ZONE T="unstruc"','N=',nPtsBodyMarker(iBody),  &
                     'E=',totNumTriElem(iBody),'F=FEPOINT  ET=TRIANGLE'
         DO i=1,nPtsBodyMarker(iBody)
           write(234,'(25(1X,1PE13.5),12(1X,I6))')xBodyMarker(i,iBody),yBodyMarker(i,iBody),zBodyMarker(i,iBody),&
                                  uBodyMarker(i,iBody),vBodyMarker(i,iBody),wBodyMarker(i,iBody), &
                                  xMarkerStress(i,iBody),yMarkerStress(i,iBody),zMarkerStress(i,iBody), &
                                  xMarkerForce(i,iBody),yMarkerForce(i,iBody),zMarkerForce(i,iBody), &
                                  u1BodyMarker(i,iBody),v1BodyMarker(i,iBody),w1BodyMarker(i,iBody), &
                                  cnt_fx(i,iBody), cnt_fy(i,iBody), cnt_fz(i,iBody), &
                                  cnt_force(i,iBody), cnt_d12(i,iBody),cnt_d21(i,iBody),cnt_d13(i,iBody), &
                                  cnt_d31(i,iBody),cnt_d23(i,iBody),cnt_d32(i,iBody), &
                                  cnt_EM12(i,iBody),cnt_EM21(i,iBody),cnt_EM13(i,iBody), &
                                  cnt_EM31(i,iBody),cnt_EM23(i,iBody),cnt_EM32(i,iBody),cnt_fg12(i,iBody), &
                                  cnt_fg21(i,iBody),cnt_fg13(i,iBody),cnt_fg31(i,iBody),cnt_fg23(i,iBody), &
                                  cnt_fg32(i,iBody)
         ENDDO ! i
 
         DO  j=1,totNumTriElem(iBody)
            WRITE(234,'(3I10)') triElemNeig(1,j,iBody),triElemNeig(2,j,iBody),triElemNeig(3,j,iBody)
         ENDDO
 
      ENDDO
      CLOSE(234)

    !END IF ! boundary_motion

    END SUBROUTINE write_markers

!------------------------------------------------    
! three components of vorticity in nlu, nlv, nlw
!------------------------------------------------    
    SUBROUTINE vorticity()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE MPI_module

    IMPLICIT NONE
    
    INTEGER           :: i,j,k,ok
    REAL(KIND=CGREAL) :: un,us,uf,ub
    REAL(KIND=CGREAL) :: ve,vw,vf,vb
    REAL(KIND=CGREAL) :: wn,ws,we,ww
    REAL(KIND=CGREAL) :: flag,dmy,LMD_MAX
    REAL(KIND=CGREAL) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
    COMPLEX*16        :: gradV(3,3), eigV(3), DUMMY(1,1), WORK2(6)

    nlu(0:,:,:) = 0.0_CGREAL
    nlv(0:,:,:) = 0.0_CGREAL
    nlw(0:,:,:) = 0.0_CGREAL
    lmd(0:,:,:) = 0.0_CGREAL    

    DO k = zc_start,zc_end   !1,nz-1
    DO j = yc_start,yc_end
    DO i = 1,nx-1
      wn = ( fy(j+1)*w(i,j+1,k) + (1.0_CGREAL-fy(j+1))*w(i,j,k)   ) 

      ws = ( fy(j)  *w(i,j,k)   + (1.0_CGREAL-fy(j)  )*w(i,j-1,k) )

      vf = ( fz(k+1)*v(i,j,k+1) + (1.0_CGREAL-fz(k+1))*v(i,j,k)   )

      vb = ( fz(k)  *v(i,j,k)   + (1.0_CGREAL-fz(k)  )*v(i,j,k-1) )

      we = ( fx(i+1)*w(i+1,j,k) + (1.0_CGREAL-fx(i+1))*w(i,j,k)   )

      ww = ( fx(i)  *w(i,j,k)   + (1.0_CGREAL-fx(i)  )*w(i-1,j,k) )

      uf = ( fz(k+1)*u(i,j,k+1) + (1.0_CGREAL-fz(k+1))*u(i,j,k)   )

      ub = ( fz(k)  *u(i,j,k)   + (1.0_CGREAL-fz(k)  )*u(i,j,k-1) )

      ve = ( fx(i+1)*v(i+1,j,k) + (1.0_CGREAL-fx(i+1))*v(i,j,k)   )

      vw = ( fx(i)  *v(i,j,k)   + (1.0_CGREAL-fx(i)  )*v(i-1,j,k) )

      un = ( fy(j+1)*u(i,j+1,k) + (1.0_CGREAL-fy(j+1))*u(i,j,k)   )

      us = ( fy(j)  *u(i,j,k)   + (1.0_CGREAL-fy(j)  )*u(i,j-1,k) )

 
      flag = (1.0_CGREAL-iblank(i,j,k)) !*(1.0_CGREAL-ghostCellMark(i,j,k))

      dudx = (face_u(i+1,j,k)-face_u(i,j,k))*dxinv(i)
      dvdy = (face_v(i,j+1,k)-face_v(i,j,k))*dyinv(j)
      dwdz = (face_w(i,j,k+1)-face_w(i,j,k))*dzinv(k)

      dudy = ( un - us )*dyinv(j)
      dudz = ( uf - ub )*dzinv(k)
      dvdz = ( vf - vb )*dzinv(k) 
      dvdx = ( ve - vw )*dxinv(i)
      dwdy = ( wn - ws )*dyinv(j) 
      dwdx = ( we - ww )*dxinv(i)

      gradV(1,1) = dudx
      gradV(2,2) = dvdy
      gradV(3,3) = dwdz 

      gradV(1,2) = dudy
      gradV(1,3) = dudz

      gradV(2,1) = dvdx
      gradV(2,3) = dvdz

      gradV(3,1) = dwdx
      gradV(3,2) = dwdy


      !nlu(i,j,k) =(  ( wn - ws )*dyinv(j)       &
      !              -( vf - vb )*dzinv(k) )*flag
      nlu(i,j,k) =( dwdy - dvdz) *flag
      nlv(i,j,k) =( dudz - dwdx) *flag
      nlw(i,j,k) =( dvdx - dudy) *flag

      !nlv(i,j,k) =(  ( uf - ub )*dzinv(k)  &
      !              -( we - ww )*dxinv(i) )*flag


      !nlw(i,j,k) =(  ( ve - vw )*dxinv(i)  & 
      !              -( un - us )*dyinv(j) )*flag

      !...Find the solution using the LAPACK routine ZGEEV
      CALL ZGEEV('N', 'N', 3, gradV, 3, eigV, DUMMY, 1, DUMMY,      &
                  1, WORK2, 6, WORK2, ok)
       
      !... Output the Eigenvalues
      IF (ok .eq. 0) THEN
         LMD_MAX = MAX ( AIMAG(eigV(1)),AIMAG(eigV(2)),AIMAG(eigV(3)) )
      ELSE
         WRITE (*,*) "An error occured"
      ENDIF

      lmd(i,j,k) = LMD_MAX
    ENDDO
    ENDDO
    ENDDO

    !Update buffer slices for each subdomain
    call send_receive_slices_real_y(nlu, 0,nx+1,yb1,yb2,zb1,zb2,1)
    call send_receive_slices_real_z(nlu, 0,nx+1,yb1,yb2,zb1,zb2,1)
    call send_receive_slices_real_y(nlv, 0,nx+1,yb1,yb2,zb1,zb2,1)
    call send_receive_slices_real_z(nlv, 0,nx+1,yb1,yb2,zb1,zb2,1)
    call send_receive_slices_real_y(nlw, 0,nx+1,yb1,yb2,zb1,zb2,1)
    call send_receive_slices_real_z(nlw, 0,nx+1,yb1,yb2,zb1,zb2,1)
    call send_receive_slices_real_y(lmd, 0,nx+1,yb1,yb2,zb1,zb2,1)
    call send_receive_slices_real_z(lmd, 0,nx+1,yb1,yb2,zb1,zb2,1)

    END SUBROUTINE vorticity

!-------------------------------------------------------------------------------
   SUBROUTINE divergence()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE MPI_module

    IMPLICIT NONE

    INTEGER              :: i,j,k
    REAL(KIND=CGREAL)    :: dsum,divmax,divmin,flag
    INTEGER , DIMENSION(1:3):: maxl,minl

! div = [ d(U)/dx + d(V)/dy + d(W)/dz ] 

    dsum          = 0.0_CGREAL
    div(0:,:,:) = 0.0_CGREAL

    DO k = zc_start,zc_end   !1,nz-1
    DO j = yc_start,yc_end
    DO i = 1,nx-1
      div(i,j,k)    = ( ( face_u(i+1,j,k) - face_u(i,j,k) )*dxinv(i)  &
                       +( face_v(i,j+1,k) - face_v(i,j,k) )*dyinv(j)  &
                       +( face_w(i,j,k+1) - face_w(i,j,k) )*dzinv(k) )

      flag = (1.0_CGREAL-iblank(i,j,k))*(1.0_CGREAL-ghostCellMark(i,j,k))

      div(i,j,k)    = div(i,j,k) * flag
      dsum = dsum + div(i,j,k)
    ENDDO
    ENDDO
    ENDDO

    !divmax = MAXVAL(abs(div(1:nx-1,1:ny-1,1:nz-1)))
    !maxl   = MAXLOC(div(1:nx-1,1:ny-1,1:nz-1))

    !WRITE(STDOUT,'(A,E15.7,A,3(2X,I4))') 'Max Divergence       = ',divmax,' at ',maxl

   END SUBROUTINE divergence
!-------------------------------------------------------------------------------
!------------------------------------------    
    SUBROUTINE write_flowfield_test() 

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE multiuse_arrays
    USE MPI_module

    IMPLICIT NONE
   
    REAL                      :: fsmach,alp,density,xhat,yhat
    INTEGER                   :: i,j,k,iBody,n,m,nCylinder
    CHARACTER*20              :: fname1, fname2
    CHARACTER*100             :: command


    !PRINT*,'Writing out flowfield file'

    write(fname1,'(I7.7)')  ntime
    write(fname2,'(I4.4)')  lProc
    fname1 = 'q_'//trim(fname1)//'.'//fname2
    print*,'Write to file: ',fname1

    OPEN(UNIT=70,FILE=fname1,STATUS='UNKNOWN')

    !-------------------------------
    IF (format_dump == TECPLOT  ) THEN
    !-------------------------------
    IF (nDim == DIM_2D) THEN
      write(70,*)'VARIABLES="X","Y","U","V","P","OZ","IBLANK","GHOST","DEAD"'
      write(70,*)'ZONE F=POINT, I=',nx-1,', J=',ny-1
      k = 1
      do j=1,ny-1
      do i=1,nx-1
         write(70,222)xc(i),yc(j),u(i,j,k),v(i,j,k),p(i,j,k)       &
                     ,nlw(i,j,k), iblank(i,j,k),ghostcellmark(i,j,k) &
                     ,dead_cell(i,j,k)
      enddo
      enddo
    ELSE
      write(70,'(a)')'VARIABLES="X","Y","Z","U","V","W","P","IBLANK","GHOST","DEAD","B#", &
                      "xbi","ybi","zbi"'
      write(70,*)'ZONE F=POINT, I=',nx+1,', J=',jSlices+2,' K=',kSlices+2
      do k=zc_start-1,zc_end+1
      do j=yc_start-1,yc_end+1
      do i=0,nx
        write(70,223)xc(i),yc(j),zc(k),u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k)   &
                    ,iblank(i,j,k),ghostcellmark(i,j,k),dead_cell(i,j,k)     &
                    ,bodyNum(i,j,k),xBItable(i,j,k),yBItable(i,j,k),zBItable(i,j,k)
!       write(70,124)u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k)   &
!                   ,iblank(i,j,k),ghostcellmark(i,j,k),dead_cell(i,j,k)

      enddo
      enddo
      enddo
    ENDIF
    CLOSE(70)

    !System call to compress the data file
    command = 'gzip '//trim(fname1)    
    call system(command)
    print*, command

    !--------------
    ENDIF ! tecplot
    !--------------

222 format(6(2x, pe12.5),3(2x,i2)) 
223 format(7(2x,pe12.5),4(2x,i2),3(2x,pe12.5))
224 format( 4(2x,pe12.5),3(2x,i2))

    END SUBROUTINE write_flowfield_test
!------------------------------------------
