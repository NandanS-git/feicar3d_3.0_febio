!---------------------------------------------
!   SUBROUTINE make_grid()
!   SUBROUTINE metrics()
!---------------------------------------------



!---------------------------------------------
! grid convention
!
!      ny+1---------------------------------------
!          |  ||  |  |  |  |   |  |  |  |  |  ||  |
!       ny +==++==+==+==+==+===+==+==+==+==+==++==+
!          |  ||  |  |  |  |   |  |  |  |  |  ||  |
!     ^    +--++--+--+--+--+---+--+--+--+--+--++--+
!   dy|    |  ||  |  |  |  | * |  |  |  |  |  ||  |
!     -  j +--++--+--+--+--+---+--+--+--+--+--++--+
!          |  ||  |  |  |  |   |  |  |  |  |  ||  |
!        2 +--++--+--+--+--+---+--+--+--+--+--++--+
!          |  ||  |  |  |  |   |  |  |  |  |  ||  |
!        1 +==++==+==+==+==+===+==+==+==+==+==++==+
!          |  ||  |  |  |  |   |  |  |  |  |  ||  |
!        0 ---------------------------------------
!          0  1   2        i                  nx  nx+1
!                          <--->
!                           dx
!
!---------------------------------------------
   SUBROUTINE make_grid()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays

    IMPLICIT NONE

    REAL(KIND=CGREAL)    :: delta
    INTEGER :: i,j,k,junk
 
    IF (xgrid_unif == UNIFORM_GRID) then
      delta  = (xout-xOrigin)/REAL(nx-1,KIND=CGREAL)
      DO i=1,nx
        x(i) = xOrigin + REAL(i-1,KIND=CGREAL)*delta
      ENDDO
    ELSE
      OPEN(UNIT=10,FILE='xgrid.dat')
      DO i=1,nx
        read(10,*)junk,x(i)
      ENDDO

      xOrigin=x(1 )
      xout   =x(nx) 
    ENDIF
    x(0)    = xOrigin*2.0_CGREAL - x(2)
    x(nx+1) = x(nx) + x(nx)-x(nx-1)
  
    IF (ygrid_unif == UNIFORM_GRID) then
      delta  = (yout-yOrigin)/REAL(ny-1,KIND=CGREAL)
      DO i=1,ny
        y(i) = yOrigin + REAL(i-1,KIND=CGREAL)*delta
      ENDDO
    ELSE
      OPEN(UNIT=11,FILE='ygrid.dat')
      DO i=1,ny
        read(11,*)junk,y(i)
      ENDDO

      yOrigin=y(1 )
      yout   =y(ny)
    ENDIF
    y(0)    = yOrigin*2.0_CGREAL - y(2)
    y(ny+1) = y(ny) + y(ny)-y(ny-1)
  
    IF (zgrid_unif == UNIFORM_GRID) then
      delta  = (zout-zOrigin)/REAL(nz-1,KIND=CGREAL)
      DO i=1,nz
        z(i) = zOrigin + REAL(i-1,KIND=CGREAL)*delta
      ENDDO
    ELSE
      OPEN(UNIT=12,FILE='zgrid.dat')
      DO i=1,nz
        read(12,*)junk,z(i)
      ENDDO

      zOrigin=z(1 )
      zout   =z(nz)
    ENDIF
    z(0)    = zOrigin*2.0_CGREAL - z(2)
    z(nz+1) = z(nz) + z(nz)-z(nz-1)

! write out grid file for plotting
!    OPEN(UNIT=9,FILE='grid_plot.dat')

!    WRITE(9,*)'VARIABLES="X","Y","Z"'
!    WRITE(9,*)'ZONE F=POINT, I=',nx,', J=',ny,' K=',nz
!    DO k=1,nz
!    DO j=1,ny
!    DO i=1,nx
!       WRITE(9,123)x(i),y(j),z(k)!
!    ENDDO
!    ENDDO
!    ENDDO
!    CLOSE(9)
123 FORMAT(3(2x,e14.7))

    dx = 0.0_CGREAL
    dy = 0.0_CGREAL
    dz = 0.0_CGREAL
    dxc= 0.0_CGREAL
    dyc= 0.0_CGREAL
    dzc= 0.0_CGREAL
 
    dxinv = 0.0_CGREAL
    dyinv = 0.0_CGREAL
    dzinv = 0.0_CGREAL
    dxcinv= 0.0_CGREAL
    dycinv= 0.0_CGREAL
    dzcinv= 0.0_CGREAL
 
    DO i=0,nx
      xc(i)    = 0.5_CGREAL*(x(i+1)+x(i)) 
      dx(i)    =             x(i+1)-x(i)  
      dxinv(i) = 1.0_CGREAL/dx(i)
    ENDDO
    xc(nx+1)= x(nx+1)
   
    DO i=0,ny
      yc(i)    = 0.5_CGREAL*(y(i+1)+y(i)) 
      dy(i)    =             y(i+1)-y(i)  
      dyinv(i) = 1.0_CGREAL/dy(i)
    ENDDO
    yc(ny+1)= y(ny+1)
   
    DO i=0,nz
      zc(i)    = 0.5_CGREAL*(z(i+1)+z(i)) 
      dz(i)    =             z(i+1)-z(i)  
      dzinv(i) = 1.0_CGREAL/dz(i)
    ENDDO
    zc(nz+1)= z(nz+1)

    DO i=1,nx
      dxc(i)    = xc(i)-xc(i-1)  
      dxcinv(i) = 1.0_CGREAL/dxc(i)
    ENDDO
   
    DO i=1,ny
      dyc(i)    = yc(i)-yc(i-1)  
      dycinv(i) = 1.0_CGREAL/dyc(i)
    ENDDO
   
    DO i=1,nz
      dzc(i)    = zc(i)-zc(i-1)  
      dzcinv(i) = 1.0_CGREAL/dzc(i)
    ENDDO

    write(100,*) 'x-grid'
    write(100,'(300F12.6)') (x(1:nx))
    write(100,*) 'xc-grid'
    write(100,'(300F12.6)') (xc(1:nx))

    write(100,*) 'y-grid'
    write(100,'(300F12.6)') (y(1:ny))
    write(100,*) 'yc-grid'
    write(100,'(300F12.6)') (yc(1:ny))

    write(100,*) 'z-grid'
    write(100,'(300F12.6)') (z(1:nz))
    write(100,*) 'zc-grid'
    write(100,'(300F12.6)') (zc(1:nz))
    close(100)

   END SUBROUTINE make_grid
!------------------------------------------- 
    
!-------------------------------------------     
   SUBROUTINE metrics()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays

    IMPLICIT NONE

    INTEGER :: i,j,k

!            |--*--|--*--|--*--|-----|
!              i-1 w  i  e i+1 
!
!        g(w) = fx(i)g(i) + (1-f(i))*g(i-1)
!
! Interpolation metrics
    DO i=1,nx
      fx(i) = ( x(i) - xc(i-1) )/( xc(i) - xc(i-1) )
    ENDDO

    DO i=1,ny
      fy(i) = ( y(i) - yc(i-1) )/( yc(i) - yc(i-1) )
    ENDDO

    DO i=1,nz
      fz(i) = ( z(i) - zc(i-1) )/( zc(i) - zc(i-1) )
    ENDDO

   END SUBROUTINE metrics
!------------------------------------------------------------------------------
