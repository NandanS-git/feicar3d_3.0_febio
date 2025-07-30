!******************************************************************************
!
! Purpose: Compute velocity components and normal pressure gradient at internal 
!          body intercept points. These are needed as BC in solvers.
!
! Description: none.
!
! Input: [u,v,w]BodyMarker  = velocity components of BM points
!        closestMarker      = closest Marker for Ghost point,
!        closestMarkerRatio = closest Marker ratio for Ghost point
!
! Output: [u,v,w]BodyIntercept = velocity components of BI point,
!         dpdnBodyIntercept    = normal pressure gradient of BI point,
!         dpdtBodyIntercept    = tangential pressure gradient of BI point.
!
! Notes: Currently dpdnBodyIntercept and dpdtBodyIntercept are set to zero.
!
!******************************************************************************
!
!******************************************************************************

  SUBROUTINE GCM_SetBodyInterceptValues()

!----------------------------------------------------------------------
! Compute velocity components and normal pressure gradient at internal 
! body intercept points. These are needed as BC in solvers
!----------------------------------------------------------------------
!

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE grid_arrays
    USE GCM_arrays
    
    USE fsi_module
    USE MPI_module 
   
    IMPLICIT NONE

!... Loop variables

    INTEGER :: n, k

!... Local variables

    INTEGER           :: iG,jG,kG,iBody
    REAL(KIND=CGREAL) :: xGC, yGC, zGC, xBI, yBI, zBI
    REAL(KIND=CGREAL) :: uBIold, vBIold, wBIOld
    !INTEGER           ::llProc  
!******************************************************************************

!---------------------------
! Loop over all ghost points in the subdomain
! and compute the velocity and dp/dn 
! at the body-intercepts
!---------------------------

    DO n=1,nGhost
      iG      = iGhost(n)
      jG      = jGhost(n)
      kG      = kGhost(n)

      !xGC     = xc(iG)
      !yGC     = yc(jG)
      !zGC     = zc(kG)

      xBI = xBodyIntercept(n)
      yBI = yBodyIntercept(n)
      zBI = zBodyIntercept(n)
      
      CALL GCM_calc_BIVelocity_Unstruc( iG, jG, kG, xBI, yBI, zBI, closestElementGC(n),      &
                                       uBodyIntercept(n),vBodyIntercept(n),wBodyIntercept(n),&
                                       uBIOld,vBIOld,wBIOld)

      !Ye,debug
      !IF((ntime.eq.2504).AND.(lProc.eq.5).AND.(n.eq.1880))THEN
      !  write(77,*)iter_FSI, iG, jG, kG
      !  write(77,'(3F10.4)')uBodyIntercept(n),vBodyIntercept(n),wBodyIntercept(n)
      !  write(77,'(3F10.4)')uBIOld,vBIOld,wBIOld
      !ENDIF
 
      if(ntime > 1) then
        ! dp/dn = -D u_n /Dt

        dpdnBodyIntercept(n) = (uBodyIntercept(n) - uBIold)/dt*xBodyInterceptNorm(n)  &
                             + (vBodyIntercept(n) - vBIold)/dt*yBodyInterceptNorm(n)  &
                             + (wBodyIntercept(n) - wBIold)/dt*zBodyInterceptNorm(n)

        dpdnBodyIntercept(n) = -dpdnBodyIntercept(n)
      else
        dpdnBodyIntercept(n) = 0.0_CGREAL
      endif

      !Ye,debug
      !IF( (lProc.eq.5).AND.(n.eq.1880) )THEN
      !   write(75,*)ntime,iter_FSI
      !   write(75,*)xBI,yBI,zBI
      !   write(75,*)dpdnBodyIntercept(n)
      !ENDIF

    ENDDO ! n

    IF (MOD(ntime,nmonitor) == 0 .and. nGhost.gt.0) THEN

      !write(*,'(a,2E20.5)')'Min-Max uBI-Actual= ',MINVAL(uBodyIntercept(1:nGhost)),&
      !                                             MAXVAL(uBodyIntercept(1:nGhost))
      !write(*,'(a,2E20.5)')'Min-Max vBI-Actual= ',MINVAL(vBodyIntercept(1:nGhost)),&
      !                                             MAXVAL(vBodyIntercept(1:nGhost))
      !write(*,'(a,2E20.5)')'Min-Max wBI-Actual= ',MINVAL(wBodyIntercept(1:nGhost)),&
      !                                             MAXVAL(wBodyIntercept(1:nGhost))
      !write(*,'(a,2E20.5)')'Min-Max ds        = ',MINVAL(probeLength(1:nGhost)),   &
      !                                             MAXVAL(probeLength(1:nGhost))
      !PRINT*, 'Min-Max xBINorm-Actual= ',MINVAL(xBodyInterceptNorm(1:nGhost)),&
      !                                   MAXVAL(xBodyInterceptNorm(1:nGhost))
      !PRINT*, 'Min-Max yBINorm-Actual= ',MINVAL(yBodyInterceptNorm(1:nGhost)),&
      !                                   MAXVAL(yBodyInterceptNorm(1:nGhost))
    ENDIF ! ntime
 
!---------------------------
! Loop over all dead cells
! and compute the velocity and dp/dn 
! at the body-intercepts
!---------------------------

    DO n=1,nDead
      iG      = iDead(n)
      jG      = jDead(n)
      kG      = kDead(n)

      !xGC     = xc(iG)
      !yGC     = yc(jG)
      !zGC     = zc(kG)

      xBI = xBodyInterceptDead(n)
      yBI = yBodyInterceptDead(n)
      zBI = zBodyInterceptDead(n)
      
      CALL GCM_calc_BIVelocity_Unstruc( iG, jG, kG, xBI, yBI, zBI, closestElementDead(n),  &
                                       uBodyInterceptDead(n),vBodyInterceptDead(n),        &
                                       wBodyInterceptDead(n),                              &
                                       uBIOld,vBIOld,wBIOld)

      if(ntime > 1) then
        ! dp/dn = -D u_n /Dt

        dpdnInterceptDead(n) = (uBodyInterceptDead(n) - uBIold)/dt*xNormInterceptDead(n)  &
                             + (vBodyInterceptDead(n) - vBIold)/dt*yNormInterceptDead(n)  &
                             + (wBodyInterceptDead(n) - wBIold)/dt*zNormInterceptDead(n)

        dpdnInterceptDead(n) = -dpdnInterceptDead(n)
      else
        dpdnInterceptDead(n) = 0.0_CGREAL
      endif

      !Ye,debug
      !if ((lProc.eq.5).AND.(iG.eq.100).AND.(jG.eq.81).AND.(kG.eq.53)) then
      !   write(83,*)ntime,iter_FSI
      !   write(83,*)n,dpdnInterceptDead(n)
      !   write(83,*)uBodyMarker(3384,2),vBodyMarker(3384,2),wBodyMarker(3384,2)
      !   write(83,*)uBodyMarker(3382,2),vBodyMarker(3382,2),wBodyMarker(3382,2)
      !   write(83,*)uBodyMarker(3329,2),vBodyMarker(3329,2),wBodyMarker(3329,2)
      !   write(83,*)uBodyInterceptDead(n),vBodyInterceptDead(n),wBodyInterceptDead(n)
      !   write(83,*)uBIold,vBIold,wBIold
      !   write(83,*)xNormBodyMarker(3384,2),yNormBodyMarker(3384,2), &
      !              zNormBodyMarker(3384,2)
      !   write(83,*)xNormBodyMarker(3382,2),yNormBodyMarker(3382,2), &
      !              zNormBodyMarker(3382,2)
      !   write(83,*)xNormBodyMarker(3329,2),yNormBodyMarker(3329,2), &
      !              zNormBodyMarker(3329,2)
      !   write(83,*)xNormInterceptDead(n),yNormInterceptDead(n),zNormInterceptDead(n)
      !   write(83,*)'======================'
      !endif

    ENDDO ! n
    
  END SUBROUTINE GCM_SetBodyInterceptValues      
!------------------------------------------------------------------------------

!******************************************************************************
!
! Purpose: generalized kernel to compute the velocity components of
!          the intercept points for any generic point onto the body markers.
!
! Description: none.
!
! Input: iGP, jGP, kGP      = indices of Generic Point,
!        closestMarkerGP = closest Marker for Generic point,
!
! Output: uGP, vGP, wGP = velocity components of Generic Point.
!
! Notes: none.
!
!
!******************************************************************************

  SUBROUTINE GCM_calc_BIVelocity_Unstruc( iGBI, jGBI, kGBI, xGBI, yGBI, zGBI, closestElementGBI,          &
                                          uGBI, vGBI, wGBI, uGBIOld, vGBIOld, wGBIOld )

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    USE gcm_arrays

    IMPLICIT NONE

!... parameters

    INTEGER,           INTENT(IN)  :: iGBI, jGBI, kGBI, closestElementGBI 
    REAL(KIND=CGREAL), INTENT(IN)  :: xGBI, yGBI,zGBI
    REAL(KIND=CGREAL), INTENT(OUT) :: uGBI, vGBI, wGBI
    REAL(KIND=CGREAL), INTENT(OUT) :: uGBIold, vGBIold, wGBIOld
!   
!... loop variables

    INTEGER :: i

!... local variables

    INTEGER                           :: iBody,node1,node2,node3,nMarker
    INTEGER                           :: info
    REAL(KIND=CGREAL)                 :: cX, cY, cZ, cC, rCond
    REAL(KIND=CGREAL), DIMENSION(4,4) :: vanTri
    REAL(KIND=CGREAL), DIMENSION(4)   :: rhsTri

! use the following approach
!
!  u = a x + b y + c z + d 
!  (a,b,c,d) determined by using four conditions
!   u = u(i) at ith node for i=1,3
!   GRAD(u) . n = 0  where n is normal to plane of triangle.
!  
!******************************************************************************

    iBody   = bodyNum(iGBI,jGBI,kGBI)
    nMarker = nPtsBodyMarker(iBody)
!  
!   assume linear variation of velocity across element and then compute value at intercept.

    node1   = triElemNeig(1,closestElementGBI,iBody)
    node2   = triElemNeig(2,closestElementGBI,iBody)
    node3   = triElemNeig(3,closestElementGBI,iBody)
!
    vanTri(1,1) = xBodyMarker(node1,iBody) 
    vanTri(1,2) = yBodyMarker(node1,iBody) 
    vanTri(1,3) = zBodyMarker(node1,iBody) 
    vanTri(1,4) = 1.0_CGREAL

    vanTri(2,1) = xBodyMarker(node2,iBody) 
    vanTri(2,2) = yBodyMarker(node2,iBody) 
    vanTri(2,3) = zBodyMarker(node2,iBody) 
    vanTri(2,4) = 1.0_CGREAL

    vanTri(3,1) = xBodyMarker(node3,iBody) 
    vanTri(3,2) = yBodyMarker(node3,iBody) 
    vanTri(3,3) = zBodyMarker(node3,iBody) 
    vanTri(3,4) = 1.0_CGREAL

    vanTri(4,1) = triElemNormx(closestElementGBI,iBody) 
    vanTri(4,2) = triElemNormy(closestElementGBI,iBody) 
    vanTri(4,3) = triElemNormz(closestElementGBI,iBody) 
    vanTri(4,4) = 0.0_CGREAL

    CALL DGETRF(4, 4, vanTri,4,iPvt, info)
    CALL DGETRI(4, vanTri,4,iPvt,work, 4, info)

! compute uGBI
    rhsTri(1) = uBodyMarker(node1,iBody)
    rhsTri(2) = uBodyMarker(node2,iBody)
    rhsTri(3) = uBodyMarker(node3,iBody)
    rhsTri(4) = 0.0_CGREAL

    cX = 0.0_CGREAL  
    cY = 0.0_CGREAL  
    cZ = 0.0_CGREAL  
    cC = 0.0_CGREAL  
    DO i = 1,4
      cX  = cX + vanTri(1,i)*rhsTri(i) 
      cY  = cY + vanTri(2,i)*rhsTri(i) 
      cZ  = cZ + vanTri(3,i)*rhsTri(i) 
      cC  = cC + vanTri(4,i)*rhsTri(i) 
    ENDDO

    uGBI  = cX * xGBI + cY * yGBI + cZ * zGBI + cC

! compute vGBI
    rhsTri(1) = vBodyMarker(node1,iBody)
    rhsTri(2) = vBodyMarker(node2,iBody)
    rhsTri(3) = vBodyMarker(node3,iBody)
    rhsTri(4) = 0.0_CGREAL

    cX = 0.0_CGREAL  
    cY = 0.0_CGREAL  
    cZ = 0.0_CGREAL  
    cC = 0.0_CGREAL  
    DO i = 1,4
      cX  = cX + vanTri(1,i)*rhsTri(i) 
      cY  = cY + vanTri(2,i)*rhsTri(i) 
      cZ  = cZ + vanTri(3,i)*rhsTri(i) 
      cC  = cC + vanTri(4,i)*rhsTri(i) 
    ENDDO

    vGBI  = cX * xGBI + cY * yGBI + cZ * zGBI + cC

! compute wGBI
    rhsTri(1) = wBodyMarker(node1,iBody)
    rhsTri(2) = wBodyMarker(node2,iBody)
    rhsTri(3) = wBodyMarker(node3,iBody)
    rhsTri(4) = 0.0_CGREAL

    cX = 0.0_CGREAL  
    cY = 0.0_CGREAL  
    cZ = 0.0_CGREAL  
    cC = 0.0_CGREAL  
    DO i = 1,4
      cX  = cX + vanTri(1,i)*rhsTri(i) 
      cY  = cY + vanTri(2,i)*rhsTri(i) 
      cZ  = cZ + vanTri(3,i)*rhsTri(i) 
      cC  = cC + vanTri(4,i)*rhsTri(i) 
    ENDDO
  
    wGBI  = cX * xGBI + cY * yGBI + cZ * zGBI + cC

!-----------------
! Interpolate from the previous velocity 
!-----------------

! compute uGBIOld
    rhsTri(1) = uBodyMarkerOld(node1,iBody)
    rhsTri(2) = uBodyMarkerOld(node2,iBody)
    rhsTri(3) = uBodyMarkerOld(node3,iBody)
    rhsTri(4) = 0.0_CGREAL

    cX = 0.0_CGREAL  
    cY = 0.0_CGREAL  
    cZ = 0.0_CGREAL  
    cC = 0.0_CGREAL  
    DO i = 1,4
      cX  = cX + vanTri(1,i)*rhsTri(i) 
      cY  = cY + vanTri(2,i)*rhsTri(i) 
      cZ  = cZ + vanTri(3,i)*rhsTri(i) 
      cC  = cC + vanTri(4,i)*rhsTri(i) 
    ENDDO

    uGBIOld  = cX * xGBI + cY * yGBI + cZ * zGBI + cC

! compute vGBIOld
    rhsTri(1) = vBodyMarkerOld(node1,iBody)
    rhsTri(2) = vBodyMarkerOld(node2,iBody)
    rhsTri(3) = vBodyMarkerOld(node3,iBody)
    rhsTri(4) = 0.0_CGREAL

    cX = 0.0_CGREAL  
    cY = 0.0_CGREAL  
    cZ = 0.0_CGREAL  
    cC = 0.0_CGREAL  
    DO i = 1,4
      cX  = cX + vanTri(1,i)*rhsTri(i) 
      cY  = cY + vanTri(2,i)*rhsTri(i) 
      cZ  = cZ + vanTri(3,i)*rhsTri(i) 
      cC  = cC + vanTri(4,i)*rhsTri(i) 
    ENDDO

    vGBIOld  = cX * xGBI + cY * yGBI + cZ * zGBI + cC

! compute wGBIOld
    rhsTri(1) = wBodyMarkerOld(node1,iBody)
    rhsTri(2) = wBodyMarkerOld(node2,iBody)
    rhsTri(3) = wBodyMarkerOld(node3,iBody)
    rhsTri(4) = 0.0_CGREAL

    cX = 0.0_CGREAL  
    cY = 0.0_CGREAL  
    cZ = 0.0_CGREAL  
    cC = 0.0_CGREAL  
    DO i = 1,4
      cX  = cX + vanTri(1,i)*rhsTri(i) 
      cY  = cY + vanTri(2,i)*rhsTri(i) 
      cZ  = cZ + vanTri(3,i)*rhsTri(i) 
      cC  = cC + vanTri(4,i)*rhsTri(i) 
    ENDDO
  
    wGBIOld  = cX * xGBI + cY * yGBI + cZ * zGBI + cC


  END SUBROUTINE GCM_calc_BIVelocity_Unstruc
