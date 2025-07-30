!---------------------------------------------------------------
! subroutine enforces general mean + sinusoidal component on translational
! and angular velocity

   SUBROUTINE forced_motion(iBody)

    USE global_parameters
    USE flow_parameters

    IMPLICIT NONE

    INTEGER,INTENT (IN) ::iBody

    vxcent(iBody) = vxcentTrans(iBody)  &
                        + ampx(iBody)*SIN(2.0_CGREAL*PI*freqx(iBody)*time)
    vycent(iBody) = vycentTrans(iBody)  &
                        + ampy(iBody)*SIN(2.0_CGREAL*PI*freqy(iBody)*time)
    vzcent(iBody) = vzcentTrans(iBody)  &
                        + ampz(iBody)*SIN(2.0_CGREAL*PI*freqz(iBody)*time)

    angvx(iBody) = angvxinit(iBody)   &
                        + ampangx(iBody)   &
                        * SIN(2.0_CGREAL*PI*freqangx(iBody)*time + angxphase(iBody) )
    angvy(iBody) = angvyinit(iBody)   &
                        + ampangy(iBody)    &
                        * SIN(2.0_CGREAL*PI*freqangy(iBody)*time + angyphase(iBody) )
    angvz(iBody) = angvzinit(iBody)   &
                        + ampangz(iBody)    &
                        * SIN(2.0_CGREAL*PI*freqangz(iBody)*time + angzphase(iBody) )  

    !print*,ibody, angvx(ibody),angvy(ibody),angvz(ibody)

    CALL compute_marker_vel(iBody)

   END SUBROUTINE forced_motion
!------------------------------------------------

   SUBROUTINE forced_motion_heaving(iBody)

    USE global_parameters
    USE flow_parameters

    IMPLICIT NONE

    INTEGER,INTENT (IN) ::iBody

    vxcent(iBody) = vxcentTrans(iBody)  &
                        + ampx(iBody)*SIN(2.0_CGREAL*PI*freqx(iBody)*time)
    vycent(iBody) = vycentTrans(iBody)  &
                        + ampy(iBody)*SIN(2.0_CGREAL*PI*freqy(iBody)*time)
    vzcent(iBody) = vzcentTrans(iBody)  &
                        + ampz(iBody)*SIN(2.0_CGREAL*PI*freqz(iBody)*time)

    !print*,ibody, angvx(ibody),angvy(ibody),angvz(ibody)

    CALL compute_marker_vel(iBody)

   END SUBROUTINE forced_motion_heaving
!------------------------------------------------
  SUBROUTINE compute_marker_vel(n)
! this is a second-order accurate algorithm for motion developed by
! R. Mittal that preserves length during rotation.

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER,INTENT(IN)  :: n

    INTEGER             :: i,m,j
    REAL(KIND=CGREAL)   :: uB,vB,wB
    REAL(KIND=CGREAL)   :: temp_uB, temp_vB, temp_wB, temp_angvx, temp_angvy, temp_angvz
    REAL(KIND=CGREAL)   :: total_omega

      temp_angvx  = 0.5_CGREAL*(angvx_old(n)+angvx(n))
      temp_angvy  = 0.5_CGREAL*(angvy_old(n)+angvy(n))
      temp_angvz  = 0.5_CGREAL*(angvz_old(n)+angvz(n))
      total_omega = 1.0_CGREAL/(1.0_CGREAL + &
                    0.25_CGREAL*dt**2*(angvx(n)**2 + angvy(n)**2 + angvz(n)**2))

!write(6,'(I4,100f7.3)'),n,temp_angvx,temp_angvy,temp_angvz,total_omega
!call sleep(1)
      DO i=1,nPtsBodyMarker(n)

        temp_uB =  ( temp_angvy*(zBodyMarker(i,n)-zcent(n)) &
                   - temp_angvz*(yBodyMarker(i,n)-ycent(n)) )

        temp_vB = -( temp_angvx*(zBodyMarker(i,n)-zcent(n)) &
                   - temp_angvz*(xBodyMarker(i,n)-xcent(n)) )

        temp_wB =  ( temp_angvx*(yBodyMarker(i,n)-ycent(n)) &
                   - temp_angvy*(xBodyMarker(i,n)-xcent(n)) )

        uB =   vxcent(n)                            &
             + total_omega*( temp_uB*(1.0_CGREAL + 0.25_CGREAL*dt**2*angvx(n)**2) &
             + temp_vB*( -0.5_CGREAL*dt*angvz(n) - 0.25_CGREAL*dt**2*angvx(n)*angvy(n)) &
             + temp_wB*( -0.5_CGREAL*dt*angvy(n) + 0.25_CGREAL*dt**2*angvx(n)*angvz(n)))

        vB =   vycent(n)                            &
             + total_omega*( temp_uB*(0.5_CGREAL*dt*angvz(n) - 0.25_CGREAL*dt**2*angvx(n)*angvy(n)) &
             + temp_vB*(1.0_CGREAL + 0.25_CGREAL*dt**2*angvy(n)**2) &
             + temp_wB*( -0.5_CGREAL*dt*angvx(n) - 0.25_CGREAL*dt**2*angvy(n)*angvz(n)))

        wB =   vzcent(n)                            &
             + total_omega*( temp_uB*(0.5_CGREAL*dt*angvy(n) + 0.25_CGREAL*dt**2*angvx(n)*angvz(n)) &
             + temp_vB*( 0.5_CGREAL*dt*angvx(n) - 0.25_CGREAL*dt**2*angvy(n)*angvz(n)) &
             + temp_wB*( 1.0_CGREAL + 0.25_CGREAL*dt**2*angvz(n)**2))

        uBodyMarker(i,n) =  uB
        vBodyMarker(i,n) =  vB
        wBodyMarker(i,n) =  wB
        !write(6,'(I4,100f7.3)'),n,ub,vb,wb

        !if(abs(wb)>5) print*,vzcent(n),total_omega,angvx,angvy,angvz,temp_ub,temp_vb,temp_wb  
      ENDDO ! i
  END SUBROUTINE compute_marker_vel

!-----------------------------------------------------------
! Flow induced rigid motion.  Rigid bodies are assumed to be 
! supported by springs and dashpots.
!
! Author: Haoxiang Luo
!-----------------------------------------------------------
SUBROUTINE induced_rigid_motion(iBody)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays,     ONLY : xc,yc,zc,dx,dy,dz
    USE boundary_arrays, ONLY : iblank, uBodyMarker,vBodyMarker,wBodyMarker
    Use usr_module
    
    IMPLICIT NONE
    
    INTEGER :: i,j,k, n
    INTEGER,INTENT (IN) :: iBody
    REAL(KIND=CGREAL)   :: spring_Kx, spring_Ky, damping, mass   
   
    n = iBody
    density_ratio = density_solid(n) / density_fluid
         
!!$    ! Use variables [xyz]centConstr for spring stiffness
!!$    ! and damping coeff.
!!$    spring_Kx = xcentConstr(n) 
!!$    spring_Ky = ycentConstr(n)
!!$    damping   = zcentConstr(n)
!!$	
!!$    IF ( body_dim(n) == 2 ) THEN  
!!$       force_x = sCx(n) / zout
!!$       force_y = sCy(n) / zout
!!$    ELSE
!!$       print*, 'induced_rigid_motion: funtion not implemented!'
!!$       STOP
!!$    ENDIF
!!$
!!$    ! compute the volume
!!$    volume = 0.0_CGREAL
!!$    k = 1
!!$    DO j =2, ny-1
!!$    DO i =2, nx-1
!!$       volume = volume + dx(i)*dy(j)*iblank(i,j,k)
!!$    ENDDO
!!$    ENDDO
!!$	 
!!$    mass = density_solid(n) * volume
!!$
!!$    vxcent_prev	= vxcent(n)
!!$    vycent_prev	= vycent(n)
!!$    vzcent_prev	= vzcent(n)
!!$
!!$! Linear Momentum Equation:  m du/dt = -k*x - c*u + F 
!!$
!!$    force_x = force_x - spring_Kx * (xcent(n) - xcentinit(n))  &
!!$                      - damping   * vxcent(n)
!!$
!!$    force_y = force_y - spring_Ky * (ycent(n) - ycentinit(n))  &
!!$                      - damping   * vycent(n)
!!$	 
!!$    vxcent(n) = vxcent_prev + dt * force_x / mass
!!$    vycent(n) = vycent_prev + dt * force_y / mass
!!$    vzcent(n) = 0.0_CGREAL
!!$	
!!$    DO i=1,nPtsBodyMarker(n)
!!$       uBodyMarker(i,n) = vxcent(n)
!!$       vBodyMarker(i,n) = vycent(n)
!!$       wBodyMarker(i,n) = vzcent(n)
!!$    ENDDO
!!$
!!$    write(*,*)'Rigid Motion: time cx, cy, vx, vy, fx, fy' 
!!$    write(*,  '(7(F12.6,1X))')time,xcent(n),ycent(n),vxcent(n),vycent(n),force_x,force_y
!!$    write(ifuBodyTrace,'(7(F12.6,1X))')  &
!!$                              time,xcent(n),ycent(n),vxcent(n),vycent(n),force_x,force_y

   END SUBROUTINE induced_rigid_motion
!----------------------------------------------------------   
   SUBROUTINE user_marker_vel(iBody)

    USE global_parameters
    USE flow_parameters

    IMPLICIT NONE

    INTEGER,INTENT (IN) ::iBody

  END SUBROUTINE user_marker_vel
