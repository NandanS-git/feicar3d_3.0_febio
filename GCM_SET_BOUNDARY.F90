!------------------------------------------------------------
  SUBROUTINE GCM_set_internal_boundary()
!
! Give a set of marker points, this subroutine computes
!   1) normal intercepts and associated geometrical info. from ghost nodes to body
!   2) Image point location
!   3) Weights in stencil for computing values at the image points
!

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE gcm_arrays
    USE MPI_module
    USE fsi_module
        
    IMPLICIT NONE

!... loop variables

    INTEGER :: i,iBody,iRow,j,k,m,n

!... local variables

    INTEGER :: iG, jG, kG, nbdr, iCIndx, jCIndx, kCIndx , iCIndxS, jCIndxS
    INTEGER :: iMin, iMax, jMin, jMax, kMin, kMax
    INTEGER :: iM, iP, jM, jP, kM, kP
    INTEGER :: iRange, jRange, kRange

    REAL(KIND=CGREAL) :: cosTheta,dsIntercept,sinTheta,               &
                         xBI,xBIN,xGC,xIP,xIPS,yBI,yBIN,yGC,yIP,yIPS,zBI,zBIN,zGC,zIP, &
                         slopeX, slopeY, slopeZ, maxDelta
    REAL(KIND=CGREAL) :: distBIElem, dotNorm

    INTEGER :: sum_ib, kk, jj, ii, i2, j2, k2, nloop, iloop
    INTEGER :: Bloop, jloop
    !===============
    !CHARACTER*6          :: gcfile
    !CHARACTER*9          :: gcpfile
    !CHARACTER*18         :: gc
    !CHARACTER*21         :: gcp
    !===============

!*****************************************************************************************

    !print*,'GCM_set_internal_boundary: lProc=',lProc
    ! Define ghostcells.  Need to exchange iblank slices before defining ghost cells
    ! for each subdomain

    !Ye,test block================================================
    !Ye, add the block here to test send_receive_real_y/z====
    !do i=0,nx+1
    !do j=yb1,yb2
    !do k=zb1,zb2
    !   iblank_real(i,j,k) = iblank(i,j,k)/1.0
    !enddo
    !enddo
    !enddo

    !write(400+lProc,*)'VARIABLES="X","Y","Z","ibreal","IBLANK", "BodyNum"'
    !write(400+lProc,*)'ZONE F=POINT, I=',nx-1,', J=',jSlices+2,', K=',kSlices+2
    !do k=zc_start-1,zc_end+1   !1,nz-1
    !do j=yc_start-1,yc_end+1
    !do i=1,nx-1
    ! write(400+lProc,'(4(2X,F12.5),2(3X,I10))')xc(i),yc(j),zc(k),iblank_real(i,j,k),iblank(i,j,k),bodyNum(i,j,k)
    !enddo
    !enddo
    !enddo
    !close(400+lProc)

    call send_receive_slices_integer_y(iblank, 0,nx+1,yb1,yb2,zb1,zb2,2)
    call send_receive_slices_integer_y(bodyNum,0,nx+1,yb1,yb2,zb1,zb2,2)

    !Ye,test========
    !call send_receive_slices_real_y(iblank_real, 0,nx+1,yb1,yb2,zb1,zb2,2)

    !Ye,check iblank after 3/2D smoothing
    !IF (ntime.eq.2502) THEN
    !write(500+lProc,*)'VARIABLES="X","Y","Z","ibreal","IBLANK", "BodyNum"'
    !write(500+lProc,*)'ZONE F=POINT, I=',nx-1,', J=',jSlices+2,', K=',kSlices+2
    !do k=zc_start-1,zc_end+1   !1,nz-1
    !do j=yc_start-1,yc_end+1
    !do i=1,nx-1
    ! write(500+lProc,'(4(2X,F12.5),2(3X,I10))')xc(i),yc(j),zc(k),iblank_real(i,j,k),iblank(i,j,k),bodyNum(i,j,k)
    !enddo
    !enddo
    !enddo
    !close(500+lProc)
    !ENDIF
 
    !call MPI_BARRIER(flow_comm,ierr)
    !stop

    call send_receive_slices_integer_z(iblank, 0,nx+1,yb1,yb2,zb1,zb2,2)
    call send_receive_slices_integer_z(bodyNum,0,nx+1,yb1,yb2,zb1,zb2,2)

    !write(600+lProc,*)'VARIABLES="X","Y","Z","IBUNDEC","IBLANK", "BodyNum"'
    !write(600+lProc,*)'ZONE F=POINT, I=',nx-1,', J=',jSlices+2,', K=',kSlices+2
    !do k=zc_start-1,zc_end+1   !1,nz-1
    !do j=yc_start-1,yc_end+1
    !do i=1,nx-1
    ! write(600+lProc,'(3(2X,F12.5),3(3X,I10))')xc(i),yc(j),zc(k),iblankUndecided(i,j,k),iblank(i,j,k),bodyNum(i,j,k)
    !enddo
    !enddo
    !enddo
    !close(600+lProc)

    !===============BIG SMOOTH LOOP==============
    !Bloop = 2 

    !DO jloop=1, Bloop

    call smooth_iblank()

    !Ye,check iblank after 3/2D smoothing
    !IF (ntime.eq.2502) THEN
    !write(700+lProc,*)'VARIABLES="X","Y","Z","IBUNDEC","IBLANK", "BodyNum"'
    !write(700+lProc,*)'ZONE F=POINT, I=',nx-1,', J=',jSlices+2,', K=',kSlices+2
    !do k=zc_start-1,zc_end+1   !1,nz-1
    !do j=yc_start-1,yc_end+1
    !do i=1,nx-1
    ! write(700+lProc,'(3(2X,F12.5),3(3X,I10))')xc(i),yc(j),zc(k),iblankUndecided(i,j,k),iblank(i,j,k),bodyNum(i,j,k)
    !enddo
    !enddo
    !enddo
    !close(700+lProc)
    !ENDIF

    !Ye,2D smooth_iblank along x direction
    nloop = 5
    nbdr  = 0
    !x-direction
    DO iloop=1, nloop
    DO i = 1, nx-1
       DO j = yc_start, yc_end
          DO k = zc_start,zc_end
             IF (iblank(i,j,k).eq.1) CYCLE
             sum_ib = 0
             DO kk = -1,1
             DO jj = -1,1
                j2 = j + jj
                k2 = k + kk
                if(abs(jj)+abs(kk) /= 1) cycle ! skip corner points
                if(iblank(i,j2,k2) == 1) then
                  sum_ib = sum_ib + 1
                  if(bodyNum(i,j2,k2) > 0) ibody = bodyNum(i,j2,k2)
                endif
             ENDDO!jj
             ENDDO!kk

         if(sum_ib >= 3) then
            iblank(i,j,k) = 1 ! mark as iblank cell
            bodyNum(i,j,k) = ibody
            nbdr = nbdr + 1
         endif

         ENDDO!k
       ENDDO!j
    ENDDO!i

    call send_receive_slices_integer_y(iblank, 0,nx+1,yb1,yb2,zb1,zb2,2)
    call send_receive_slices_integer_y(bodyNum,0,nx+1,yb1,yb2,zb1,zb2,2)

    call send_receive_slices_integer_z(iblank, 0,nx+1,yb1,yb2,zb1,zb2,2)
    call send_receive_slices_integer_z(bodyNum,0,nx+1,yb1,yb2,zb1,zb2,2)

    ENDDO!iloop

    !call MPI_BARRIER(flow_comm,ierr)
    !if (lProc.eq.0)then
    !   write(*,*)'2d x-direction smooth OK...'
    !endif

    !Ye,check iblank after 3/2D smoothing
    !write(800+lProc,*)'VARIABLES="X","Y","Z","IBUNDEC","IBLANK", "BodyNum"'
    !write(800+lProc,*)'ZONE F=POINT, I=',nx-1,', J=',jSlices+2,', K=',kSlices+2
    !do k=zc_start-1,zc_end+1   !1,nz-1
    !do j=yc_start-1,yc_end+1
    !do i=1,nx-1
    ! write(800+lProc,'(3(2X,F12.5),3(3X,I10))')xc(i),yc(j),zc(k),iblankUndecided(i,j,k),iblank(i,j,k),bodyNum(i,j,k)
    !enddo
    !enddo
    !enddo
    !close(800+lProc)

    !y-direction
    DO iloop=1, nloop
    DO j = yc_start, yc_end
       DO i = 1, nx-1
          DO k = zc_start,zc_end
             IF (iblank(i,j,k).eq.1) CYCLE
             sum_ib = 0
             DO kk = -1,1
             DO ii = -1,1
                k2 = k + kk
                i2 = i + ii
                if(abs(ii)+abs(kk) /= 1) cycle ! skip corner points
                if(iblank(i2,j,k2) == 1) then
                  sum_ib = sum_ib + 1
                  if(bodyNum(i2,j,k2) > 0) ibody = bodyNum(i2,j,k2)
                endif
             ENDDO!kk
             ENDDO!ii

          if(sum_ib >= 3) then
            iblank(i,j,k) = 1 ! mark as iblank cell
            bodyNum(i,j,k) = ibody
            nbdr = nbdr + 1
         endif

          ENDDO!k
       ENDDO!i
    ENDDO!j

    call send_receive_slices_integer_y(iblank, 0,nx+1,yb1,yb2,zb1,zb2,2)
    call send_receive_slices_integer_z(iblank, 0,nx+1,yb1,yb2,zb1,zb2,2)

    call send_receive_slices_integer_y(bodyNum,0,nx+1,yb1,yb2,zb1,zb2,2)
    call send_receive_slices_integer_z(bodyNum,0,nx+1,yb1,yb2,zb1,zb2,2)

    ENDDO!iloop

    !Ye,debug
    !call MPI_BARRIER(flow_comm,ierr)
    !if (lProc.eq.0)then
    !   write(*,*)'2d y-direction smooth OK...'
    !endif

    IF (MOD(ntime,nmonitor) == 0 .and. nbdr>0) &
      print*,'2D-Smooth_iblank: number of isolated cells removed:',nbdr

    !ENDDO !BIG SMOOTH LOOP

    !Ye,check iblank after 3/2D smoothing
    !IF (ntime.eq.836) THEN
    !write(7000+lProc,*)'VARIABLES="X","Y","Z","xbi","ybi","zbi","IBLANK","G","D","BNum","isBI"'
    !write(7000+lProc,*)'ZONE F=POINT, I=',nx-1,', J=',jSlices+2,', K=',kSlices+2
    !do k=zc_start-1,zc_end+1   !1,nz-1
    !do j=yc_start-1,yc_end+1
    !do i=1,nx-1
    ! write(7000+lProc,'(6(2X,F12.5),3(3X,I10))')xc(i),yc(j),zc(k), &
    !      xBItable(i,j,k),yBItable(i,j,k),zBItable(i,j,k), &
    !      iblank(i,j,k),bodyNum(i,j,k),is_BIset(i,j,k)
    !enddo
    !enddo
    !enddo
    !close(7000+lProc)
    !call MPI_BARRIER(flow_comm,ierr)
    !STOP
    !ENDIF
!========================================================================

    call GCM_define_ghostcells()

    !write(800+lProc,*)'VARIABLES="X","Y","Z","gc","IBLANK", "BodyNum"'
    !write(800+lProc,*)'ZONE F=POINT, I=',nx-1,', J=',jSlices+2,', K=',kSlices+2
    !do k=zc_start-1,zc_end+1   !1,nz-1
    !do j=yc_start-1,yc_end+1
    !do i=1,nx-1
    ! write(800+lProc,'(3(2X,F12.5),3(3X,I10))')xc(i),yc(j),zc(k),ghostCellMark(i,j,k),iblank(i,j,k),bodyNum(i,j,k)
    !enddo
    !enddo
    !enddo
    !close(800+lProc)
    !call MPI_BARRIER(flow_comm,ierr)
    !stop

    ! count the number of ghostcells
    nbdr = 0
    DO k = zc_start,zc_end    !1, nz-1
    DO j = yc_start,yc_end
    DO i = 1, nx-1
       IF ( ghostCellMark(i,j,k) == 1 )   THEN
          nbdr = nbdr + 1
       ENDIF ! ghostCellMark    
    ENDDO ! i 
    ENDDO ! j
    ENDDO ! k
      
    nGhost = nbdr         ! total number of ghost points

      IF ( MOD(ntime,nmonitor) == 0 .OR. ntime == 1) &
         PRINT*,'lProc',lProc,'GCM_set_internal_boundary: nGhost = ',nGhost 

    ! Deallocate arrays pertinent to Ghost Cells 
    IF ( ntime >= ntime_start+1 ) THEN
       CALL GCM_DeallocateGhostCellArrays()
    ENDIF ! ntime

    ! Allocate Arrays pertinent to Ghost Cells
    CALL GCM_AllocateGhostCellArrays()

    ! Set appropriate values for iGhost and jGhost by doing search
    nbdr = 0
    DO k = zc_start,zc_end    !1, nz-1
    DO j = yc_start,yc_end
    DO i = 1, nx-1
       IF ( ghostCellMark(i,j,k) == 1 )   THEN
          nbdr         = nbdr + 1
          iGhost(nbdr) = i
          jGhost(nbdr) = j
          kGhost(nbdr) = k

          ghostCellIndex(i,j,k) = nbdr
       ENDIF ! ghostCellMark
    ENDDO ! i 
    ENDDO ! j
    ENDDO ! k

    !call write_dump()
    !call MPI_BARRIER(flow_comm,ierr)
!    call sleep(20)
    !stop

!---------------------------------------------	
! project each ghostnode onto the boundary  
!---------------------------------------------

    nbdr=0
    DO n = 1, nGhost

        iG = iGhost(n) 
        jG = jGhost(n)
        kG = kGhost(n)
                
        iCellIndex(n) = -1
        jCellIndex(n) = -1
        kCellIndex(n) = -1

        iBody = bodyNum(iG,jG,kG)
 
        xGC = xc(iG)
        yGC = yc(jG)
        zGC = zc(kG)
   
        !if( (is_BIset(iG,jG,kG)==1).and.(iblank(iG,jG,kG)==1) )then
        !   !Use the BI information from the saved table
        !   xBI  = xBItable(iG,jG,kG)
        !   yBI  = yBItable(iG,jG,kG)
        !   zBI  = zBItable(iG,jG,kG)
        !   closestElementGC(n)  = cELtable(iG,jG,kG)

        !else      
           !write(*,'(4I4,I8,6F10.5)')iG,jG,kG,ghostcellMark(iG,jG,kG), &
           !     iblankUndecided(iG,jG,kG),xGC,yGC,zGC,xBI,yBI,zBI
        !   nbdr=nbdr+1

           !Use the Body intercept subroutine to determine BI.
           !Note this might cause problem if IBLANK is not updated during
           !multiple time steps.
           CALL calc_bodyIntercept_Unstruc(iBody, iG,jG,kG, xBI,yBI,zBI, &
                                           closestElementGC(n), distBIElem, dotNorm)

           !print*, 'iG, jG, kG, closestElement = ', iG, jG, kG, closestElementGC(n)
     
           IF( body_dim(iBody) == 2 )  zBI = zGC ! Enforce for 2D Body
        !endif ! is_BIset 

        ! Extract coordinates of Body Intercept.
        ! Note that, for a membrane structure, the true BI point 
          
        ! is saved.
        xBodyIntercept(n) = xBI
        yBodyIntercept(n) = yBI
        zBodyIntercept(n) = zBI !store BI info for each GhostCell

        !-- Get length of normal
        dsIntercept = SQRT( (xGC-xBI)**2 + (yGC-yBI)**2 + (zGC-zBI)**2 )

        slopeX = (xGC - xBI)/dsIntercept
        slopeY = (yGC - yBI)/dsIntercept
        slopeZ = (zGC - zBI)/dsIntercept

        ! set up virtual intercept for membrane-type structure
        if((unstruc_surface_type(iBody)==MEMBRANE)) then
           xBI = xBI + membrane_tkns/2.0_CGREAL * slopeX
           yBI = yBI + membrane_tkns/2.0_CGREAL * slopeY
           zBI = zBI + membrane_tkns/2.0_CGREAL * slopeZ

           dsIntercept = SQRT( (xGC-xBI)**2 + (yGC-yBI)**2 + (zGC-zBI)**2 )
        endif

        ! intercept normal (pointing into the body)
        xBodyInterceptNorm(n) = -slopeX
        yBodyInterceptNorm(n) = -slopeY
        zBodyInterceptNorm(n) = -slopeZ

        !-- Now compute location of probe-tip (Image Point)
        !-- Equation of 3D line  (parametric form)
        !-- x = xo + slopeX . s
        !-- y = yo + slopeY . s
        !-- z = zo + slopeZ . s

         probeLength(n) = dsIntercept*probeLengthNormalized

         !xImagePoint(n) = xBI + slopeX*probeLength(n)
         !yImagePoint(n) = yBI + slopeY*probeLength(n)
         !zImagePoint(n) = zBI + slopeZ*probeLength(n)

         xImagePoint(n) = xBI - (xBI-xGC)*probeLengthNormalized !1.1
         yImagePoint(n) = yBI - (yBI-yGC)*probeLengthNormalized
         zImagePoint(n) = zBI - (zBI-zGC)*probeLengthNormalized

         xIP = xImagePoint(n)
         yIP = yImagePoint(n)
         zIP = zImagePoint(n)

         !Ye,check the IP by using xyzDCIP, shared with DCs
         xyzDCIP(iG,jG,kG,0) = xIP
         xyzDCIP(iG,jG,kG,1) = yIP
         xyzDCIP(iG,jG,kG,2) = zIP

         ! Find the lower left grid point to Image Point in Physical domain
         maxDelta = MAX(dx(iG),dy(jG),dz(kG))

         !if(dsIntercept > maxDelta)  &
         !write(113,'(I5,12F10.5)')n,xGC,yGC,zGC,xBI,yBI,zBI,xIP,yIP,zIP,dsIntercept,probeLength(n)

         iRange  = 2*NINT(maxDelta/dx(iG) + 0.5_CGREAL)
         jRange  = 2*NINT(maxDelta/dy(jG) + 0.5_CGREAL)
         kRange  = 2*NINT(maxDelta/dz(kG) + 0.5_CGREAL)

         iMin = iG-iRange
         iMax = iG+iRange
         iMin = MAX(iMin,0)    ! note image point is allowed to be between xc(0) and x(1)
         iMax = MIN(iMax,nx)   ! note image point is allowed to be between x(nx) and xc(nx)

         DO i = iMin,iMax 
           IF ( ( xc(i) <= xIP .AND. xc(i+1) > xIP ) ) THEN
             iCellIndex(n) = i
           ENDIF ! xc
         ENDDO ! i

         jMin = jG-jRange
         jMax = jG+jRange
         jMin = MAX(jMin,yc_start-2)
         jMax = MIN(jMax,yc_end  +1)

         DO j = jMin,jMax
           IF ( ( yc(j) <= yIP .AND. yc(j+1) > yIP ) ) THEN
             jCellIndex(n) = j
           ENDIF ! xc
         ENDDO ! j

         kMin = kG-kRange
         kMax = kG+kRange

         kMin = MAX(kMin,zc_start-2)
         kMax = MIN(kMax,zc_end  +1)

!        IF (n==1) THEN
!          print*,iMin,iMax
!          print*,jMin,jMax
!          print*,kMin,kMax
!        ENDIF

         DO k = kMin,kMax
           IF ( ( zc(k) <= zIP .AND. zc(k+1) > zIP ) ) THEN
             kCellIndex(n) = k
           ENDIF ! xc
         ENDDO ! k

         !Ye,debug
         !IF( (ntime.eq.4803).AND.(lProc.eq.4).AND.(iG.eq.81).AND.  &
         !    (jG.eq.81).AND.(kG.eq.50) ) THEN
         !   write(74,*)'INFO:',ntime, iter_FSI, n
         !   write(74,*)iG, jG, kG, ghostCellIndex(iG,jG,kG)
         !   write(74,*)xBodyIntercept(n),yBodyIntercept(n),zBodyIntercept(n)
         !   write(74,*)xBI,yBI,zBI
         !   write(74,*)xIP,yIP,zIP
         !   write(74,*)iCellIndex(n),jCellIndex(n),kCellIndex(n)
         !ENDIF

         IF ( ( iCellIndex(n) ==-1 .OR. jCellIndex(n) ==-1 .OR. kCellIndex(n) ==-1) ) THEN
            PRINT*,'Failed to find the hosting box for an image point of a ghost node'           
            PRINT*,ntime,lProc,n,iBody,iG, jG, kG
            PRINT*,iG, jG, kG, iBody, closestElementGC(n)
            PRINT*,iCellIndex(n),jCellIndex(n),kCellIndex(n)
            PRINT*,xgc,ygc,zgc
            PRINT*,'BIt',xBodyIntercept(n),yBodyIntercept(n),zBodyIntercept(n)
            PRINT*,'BIv',xBI,yBI,zBI
            PRINT*,xIP,yIP,zIP
            PRINT*, 'Nandan: xc(iMin), xc(iMax)', xc(iMin), xc(iMax)
            PRINT*, 'Nandan: iMin, iMax, iRange', iMin, iMax, iRange
            PRINT*,'Aborting Run; write dump...'
            CALL write_dump()
            call sleep(15)
            !call MPI_BARRIER(flow_comm,ierr)
            !STOP
         ENDIF

!!$         ! all 8 nodes should be fluid nodes
!!$         IF (sum(iblank(iCellIndex(n):iCellIndex(n)+1,  &
!!$                        jCellIndex(n):jCellIndex(n)+1,  &
!!$                        kCellIndex(n):kCellIndex(n)+1)) ) THEN
!!$            PRINT*,'The interpolation stencil for the image point includes solid cell!'
!!$            PRINT*,n
!!$            PRINT*,iG, jG, kG
!!$            PRINT*,xgc,ygc,zgc
!!$            PRINT*,xBI,yBI,zBI
!!$            PRINT*,xIP,yIP,zIP
!!$            PRINT*,'Aborting Run'
!!$            CALL write_dump()
!!$            STOP
!!$         ENDIF

        ! Perform bilinear interpolation
        iCIndx = iCellIndex(n)
        jCIndx = jCellIndex(n)
        kCIndx = kCellIndex(n)

        !xBIN = triElemNormx(closestElementGC(n),iBody)
        !yBIN = triElemNormy(closestElementGC(n),iBody)
        !zBIN = triElemNormz(closestElementGC(n),iBody)
        xBIN = xBodyInterceptNorm(n)
        yBIN = yBodyInterceptNorm(n)
        zBIN = zBodyInterceptNorm(n)

        CALL GCM_Calc_vanMatrixDN( iG, jG, kG, iCIndx, jCIndx, kCIndx,             &
                                   xGC, yGC, zGC, xBI, yBI, zBI, xBIN, yBIN, zBIN, &
                                   coeffGCMDirc, coeffGCMNeum      )

        ! For Van Kan to work for 2D bodies, I found that coeffGCMD(5:8,:) 
        ! and coeffGCMN(5:8,:) have to be exactly zeroed out.
        ! Also, coeffGCMD and coeffGCMN have to be exactly the same for z layers
        ! -- that is, the geometry is strictly 2D.
        !
        ! H. Luo, Aug. 16, 2006

        coeffGCMD(1:iRowMax,n) = coeffGCMDirc(1:iRowMax)
        coeffGCMN(1:iRowMax,n) = coeffGCMNeum(1:iRowMax)
        
        !write(200,'(I8,1X,12I5)')n,iG,jG,kG,iCIndx,jCIndx,kCIndx
        !write(200,'(I8,1X,12F15.5)')n,xBI,yBI,zBI,xIP,yIP,zIP
      ENDDO ! n

      !print*,lProc,'BItable not set for ghostcells: nbdr/nGhost=',nbdr,nGhost

      ! Set up dead-cells
      Call GCM_set_deadcells()

!!$      do n=1,nGhost
!!$         iG = iGhost(n) 
!!$         jG = jGhost(n)
!!$         kG = kGhost(n)
!!$
!!$         xGC = xc(iG)
!!$         yGC = yc(jG)
!!$         zGC = zc(kG)
!!$
!!$         iCIndx = iCellIndex(n)
!!$         jCIndx = jCellIndex(n)
!!$         kCIndx = kCellIndex(n)
!!$         
!!$         xBI = xBodyIntercept(n)
!!$         yBI = yBodyIntercept(n)
!!$         zBI = zBodyIntercept(n)
!!$         
!!$         xIP = xImagePoint(n)
!!$         yIP = yImagePoint(n)
!!$         zIP = zImagePoint(n)
!!$
!!$         do iRow= 1, iRowMax
!!$            i  = iCIndx + incI(iRow)
!!$            j  = jCIndx + incJ(iRow)
!!$            k  = kCIndx + incK(iRow)
!!$            if(dead_cell(i,j,k)==-1) then
!!$               write(200,'(I8,1X,12I5)')n,iG,jG,kG,iCIndx,jCIndx,kCIndx,i,j,k
!!$               write(200,'(I8,1X,12F15.5)')n,xGC,yGC,zGC,xBI,yBI,zBI,xIP,yIP,zIP    
!!$            endif
!!$            if(iblank(i,j,k)==1 .and. dead_cell(i,j,k)==0) then
!!$               write(200,'(a)')'====================='
!!$               write(200,'(I8,1X,12I5)')n,iG,jG,kG,iCIndx,jCIndx,kCIndx,i,j,k
!!$               write(200,'(I8,1X,12F15.5)')n,xGC,yGC,zGC,xBI,yBI,zBI,xIP,yIP,zIP                
!!$               write(200,'(a)')'====================='
!!$            endif
!!$         enddo
!!$      enddo


    !Ye,check iblank after 3/2D smoothing
    !IF (ntime.eq.836) THEN
    !write(7000+lProc,*)'VARIABLES="X","Y","Z","xb","yb","zb","xp","yp","zp","U","B","G","D","BN","is"'
    !write(7000+lProc,*)'ZONE F=POINT, I=',nx-1,', J=',jSlices+2,', K=',kSlices+2
    !do k=zc_start-1,zc_end+1   !1,nz-1
    !do j=yc_start-1,yc_end+1
    !do i=1,nx-1
    ! write(7000+lProc,'(9(2X,F12.5),6(3X,I10))')xc(i),yc(j),zc(k), &
    !      xBItable(i,j,k),yBItable(i,j,k),zBItable(i,j,k), &
    !      xyzDCIP(i,j,k,0),xyzDCIP(i,j,k,1),xyzDCIP(i,j,k,2), &
    !      iblankUndecided(i,j,k), &
    !      iblank(i,j,k),ghostcellMark(i,j,k),dead_cell(i,j,k), &
    !      bodyNum(i,j,k),is_BIset(i,j,k)
    !enddo
    !enddo
    !enddo
    !close(7000+lProc)
    !call MPI_BARRIER(flow_comm,ierr)
    !STOP
    !ENDIF

  END SUBROUTINE GCM_set_internal_boundary 
!----------------------------------------------------------------------

  SUBROUTINE GCM_Calc_VanMatrixDN( iG, jG, kG, iCIndex,jCIndex, kCIndex,      &
                                  xIP, yIP, zIP, xBI, yBI, zBI, xBIN, yBIN, zBIN, &
                                  coeffDirc, coeffNeum      )
    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE gcm_arrays
    USE boundary_arrays
    USE MPI_module

    IMPLICIT NONE

!... parameters variables

    INTEGER, INTENT(IN)           :: iG, jG, kG, iCIndex,jCIndex, kCIndex
    REAL(KIND=CGREAL), INTENT(IN) :: xIP, yIP, zIP, xBI, yBI, zBI, xBIN, yBIN, zBIN
    REAL(KIND=CGREAL), DIMENSION(iRowMax), INTENT(OUT) :: coeffDirc, &
                                                          coeffNeum

!... loop variables

    INTEGER :: i,j,k,iRow
    INTEGER :: info

!... local variables
    
    REAL(KIND=CGREAL) :: rCond, xC1,xC2,xC3, xN1,xN2,xN3
  
!*****************************************************************************************
  
!   |-------|-------|---/---|-------|--         N : Nth ghost point
!   |   ii  |  iii  |  *    |       |           * : markers
!   |   0...|...O   | / .   |   .   |           O : other nodes used in bilinear interpolation
!   |   .   |   .   |*      |       |           + : probe tip (Image Point) 
!   |---.---|--+.---/-------|-------|--
!   |   .   |   .  *|       |       |
!   |   0...|. .O / |   N   |   .   |
!   |   i   |  iv*  |       |       |
!   |-------| --/ --|-------|-------|--

! interpolant      U = a X X X  + b X X  + c X X  + d X X
!                         1 2 3      1 2      1 3      2 3
!
!                    + e X  + f X + g X  + h
!                         1      2     3
!
!
!         [  X X     X     X   1  ]  [   ]     [     ] 
!      i  [   1 2     1     2     ]  [ a ]     [ U   ]
!         [                       ]  [   ]     [  i  ]
!         [  X X     X     X   1  ]  [   ]     [     ]
!      ii [   1 2     1     2     ]  [ b ]     [ U   ]
!         [                       ]  [   ]  =  [  ii ]
!     iii [  X X     X     X   1  ]  [   ]     [     ]
!         [   1 2     1     2     ]  [ c ]     [ U   ]
!         [                       ]  [   ]     [  iii]
!     iv  [  X X     X     X   1  ]  [   ]     [     ]
!         [   1 2     1     2     ]  [ d ]     [ U   ]
!         [                       ]  [   ]     [  iv ]
!
!
!   Van Matrix For Dirichlet conditions at Intersection Point (N)
!
!         [  X X     X     X   1  ]  [   ]     [     ] 
!      i  [   1 2     1     2     ]  [ a ]     [ U   ]
!         [                       ]  [   ]     [  i  ]
!         [  X X     X     X   1  ]  [   ]     [     ]
!      ii [   1 2     1     2     ]  [ b ]     [ U   ]
!         [                       ]  [   ]  =  [  ii ]
!     iii [  X X     X     X   1  ]  [   ]     [     ]
!         [   1 2     1     2     ]  [ c ]     [ U   ]
!         [                       ]  [   ]     [  iii]
!      N  [  X X     X     X   1  ]  [   ]     [     ]
!         [   1 2     1     2     ]  [ d ]     [ U   ]
!         [                       ]  [   ]     [  N  ]
!
!   Van Matrix For Neumann conditions at Intersection point (N)
!    B1 = n_x, B2 = n_y (components of normal vectors)
!    F_m = value of normal derivative 
!
!         [  X X           X     X   1  ]  [   ]     [     ] 
!      i  [   1 2           1     2     ]  [ a ]     [ U   ]
!         [                             ]  [   ]     [  i  ]
!         [  X X           X     X   1  ]  [   ]     [     ]
!      ii [   1 2           1     2     ]  [ b ]     [ U   ]
!         [                             ]  [   ]  =  [  ii ]
!     iii [  X X           X     X   1  ]  [   ]     [     ]
!         [   1 2           1     2     ]  [ c ]     [ U   ]
!         [                             ]  [   ]     [  iii]
!      N  [  B X  + B X    B     B   0  ]  [   ]     [     ]
!         [   1 2    2  1   1     2     ]  [ d ]     [ F   ]
!         [                             ]  [   ]     [  N  ]
!


    DO iRow= 1, iRowMax
      i  = iCIndex + incI(iRow)
      j  = jCIndex + incJ(iRow)
      k  = kCIndex + incK(iRow)

      xC1 = xc(i)
      xC2 = yc(j)
      xC3 = zc(k)

!-- Construct Vandermonde Matrices

!--- Dirichlet conditions for velocity field

      vanMatrixD(iRow,1) = xC1*xC2*xC3
      vanMatrixD(iRow,2) = xC1*xC2
      vanMatrixD(iRow,3) = xC1*xC3
      vanMatrixD(iRow,4) = xC2*xC3
      vanMatrixD(iRow,5) = xC1
      vanMatrixD(iRow,6) = xC2
      vanMatrixD(iRow,7) = xC3
      vanMatrixD(iRow,8) = 1.0_CGREAL

!--- Neumann conditions for pressure field


      vanMatrixN(iRow,1) = xC1*xC2*xC3
      vanMatrixN(iRow,2) = xC1*xC2
      vanMatrixN(iRow,3) = xC1*xC3
      vanMatrixN(iRow,4) = xC2*xC3
      vanMatrixN(iRow,5) = xC1
      vanMatrixN(iRow,6) = xC2
      vanMatrixN(iRow,7) = xC3
      vanMatrixN(iRow,8) = 1.0_CGREAL

!-- Correct For Ghost node part of cell formation, switch to Body Intercept point

      IF ( i==iG .AND. j == jG  .AND. k== kG) THEN
         xC1 = xBI
         xC2 = yBI
         xC3 = zBI
         xN1 = xBIN
         xN2 = yBIN
         xN3 = zBIN

         vanMatrixD(iRow,1) = xC1*xC2*xC3
         vanMatrixD(iRow,2) = xC1*xC2
         vanMatrixD(iRow,3) = xC1*xC3
         vanMatrixD(iRow,4) = xC2*xC3
         vanMatrixD(iRow,5) = xC1
         vanMatrixD(iRow,6) = xC2
         vanMatrixD(iRow,7) = xC3
         vanMatrixD(iRow,8) = 1.0_CGREAL

         vanMatrixN(iRow,1) = xN1*xC2*XC3 + xN2*xC1*XC3 + xN3*XC1*XC2
         vanMatrixN(iRow,2) = xN1*xC2 + xN2*xC1
         vanMatrixN(iRow,3) = xN1*xC3 + xN3*xC1
         vanMatrixN(iRow,4) = xN2*xC3 + xN3*xC2
         vanMatrixN(iRow,5) = xN1
         vanMatrixN(iRow,6) = xN2
         vanMatrixN(iRow,7) = xN3
         vanMatrixN(iRow,8) = 0.0_CGREAL

      ENDIF ! i
    ENDDO ! iRow		

! Compute inverse of Vandermonde Matrices

    CALL DGETRF(8, 8, vanMatrixD,8,iPvt, info) 
    CALL DGETRI(8, vanMatrixD,8,iPvt,work, 8, info) 

    CALL DGETRF(8, 8, vanMatrixN,8,iPvt, info)
    CALL DGETRI(8, vanMatrixN,8,iPvt,work, 8, info)

! Load Coeff-Matrices

    DO iRow = 1, iRowMax
       coeffDirc(iRow) = vanMatrixD(1,iRow)*xIP*yIP*zIP  &
                       + vanMatrixD(2,iRow)*xIP*yIP      &
                       + vanMatrixD(3,iRow)*xIP*zIP      &
                       + vanMatrixD(4,iRow)*yIP*zIP      &
                       + vanMatrixD(5,iRow)*xIP          &
                       + vanMatrixD(6,iRow)*yIP          &
                       + vanMatrixD(7,iRow)*zIP          &
                       + vanMatrixD(8,iRow)

       coeffNeum(iRow) = vanMatrixN(1,iRow)*xIP*yIP*zIP  &
                       + vanMatrixN(2,iRow)*xIP*yIP      &
                       + vanMatrixN(3,iRow)*xIP*zIP      &
                       + vanMatrixN(4,iRow)*yIP*zIP      &
                       + vanMatrixN(5,iRow)*xIP          &
                       + vanMatrixN(6,iRow)*yIP          &
                       + vanMatrixN(7,iRow)*zIP          &
                       + vanMatrixN(8,iRow)
    ENDDO ! iRow 

  END SUBROUTINE GCM_Calc_VanMatrixDN
!----------------------------------------------------------------------
!--------------------------------------------------------------
  SUBROUTINE GCM_define_ghostcells()
!
! This subroutine finds the ghostcells around solid structures,
! 
    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE gcm_arrays
    USE MPI_module
    
    IMPLICIT NONE

    INTEGER :: iBody,i,j,k,m,n, ii, jj, kk, i2,j2,k2

      ghostcellMark(:,:,:) = 0

      !Include one extra slice in z on each side of the subdomain
      !to incorporate iblank cells from adjacent subdomains.
      DO k = zc_start-1,zc_end+1    !1, nz
      !DO j = 1, ny
      !DO i = 1, nx
      !Ye,test
      DO j = yc_start-1,yc_end+1
      DO i = 1, nx-1

         IF (iblank(i,j,k) == 1) THEN  !iblank cell
            DO kk = -1,1
            DO jj = -1,1
            DO ii = -1,1

               i2 = i + ii
               j2 = j + jj
               k2 = k + kk

!!$               ! The following min-max operation may include additional ghost cells if
!!$               ! iblank is 1 at the 1st/last element.
!!$               i2 = MAX(i2,1) 
!!$               i2 = MIN(i2,nx)
!!$               
!!$               j2 = MAX(j2,1) 
!!$               j2 = MIN(j2,ny)
!!$               
!!$               k2 = MAX(k2,zc_start)  !MAX(k2,1) 
!!$               k2 = MIN(k2,zc_end)    !MIN(k2,nz)

               ! Use the following exclusion statements to replace the above min-max
               ! operations.
               !if(i2.lt.1 .or. i2.gt.nx) cycle
               !if(j2.lt.1 .or. j2.gt.ny) cycle
               !Ye,remove ghost cell located at outer boundary to prevent a
               !ghostcell index
               if(i2.lt.1 .or. i2.gt.nx-1) cycle
               if(j2.lt.yc_start .or. j2.gt.yc_end) cycle
               if(k2.lt.zc_start .or. k2.gt.zc_end) cycle

               if(abs(ii)+abs(jj)+abs(kk) /= 1) cycle  ! skip edge and corner points

               IF(iblank(i2, j2, k2) == 0) THEN
                  ! check if ghostcell is set
                  if (ghostCellMark(i2,j2,k2) == 0) then
                     ! set ghost cell
                     ghostCellMark(i2,j2,k2) = 1
                     bodyNum      (i2,j2,k2) = bodyNum(i,j,k)
                  endif
               ENDIF
            ENDDO     !ii
            ENDDO     !jj
            ENDDO     !kk
         ENDIF    !iblank(i,j,k)

      ENDDO ! i
      ENDDO ! j
      ENDDO ! k

    END SUBROUTINE GCM_define_ghostcells
!-----------------------------------------------------
!----------------------------------------------------------------------

  SUBROUTINE GCM_set_deadcells()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE gcm_arrays
    USE MPI_module
    USE fsi_module

    IMPLICIT NONE

    INTEGER :: i,j,k,iBody,iRow,m,n,nG

    INTEGER :: iG, jG, kG, nbdr, iCIndx, jCIndx, kCIndx
    INTEGER :: iMin, iMax, jMin, jMax, kMin, kMax
    INTEGER :: ii, jj, kk, iM, iP, jM, jP, kM, kP
    INTEGER :: iRange, jRange, kRange, nMarker

    REAL(KIND=CGREAL) :: xBI,xBIN,xIP,yBI,yBIN,yIP,zBI,zBIN,zIP, &
                         slopeX, slopeY, slopeZ, maxDelta
    REAL(KIND=CGREAL) :: xDC, yDC, zDC, dsIntercept, temp
    REAL(KIND=CGREAL) :: ds1
    REAL(KIND=CGREAL) :: distBIElem, dotNorm

!*****************************************************************************************

    dead_cell(:,:,:) = 0

    !Ye, pure solid node for DC
    ps1 = 0
    ps2 = 0

!! The following alogrithm to determine the deadcell doesn't work any more for subdomains, as
!! there are ghostcells in buffer slices that are not counted in the current subdomain. Instead,
!! we need to sweep through the entire subdomain, including the buffer slices.
!
!    DO n = 1, nGhost
!
!       i = iGhost(n) 
!       j = jGhost(n)
!       k = kGhost(n)
!
!       iBody = bodyNum(i,j,k)
!
!       ! -1 means that the dead-cell is associated with a ghost-cell
!       ! on the negative side of the surface ('fluid' side for a solid body)
!       if(iblank(i+1,j,k)==1)  dead_cell(i+1,j,k) = -1
!       if(iblank(i-1,j,k)==1)  dead_cell(i-1,j,k) = -1
!       if(iblank(i,j+1,k)==1)  dead_cell(i,j+1,k) = -1
!       if(iblank(i,j-1,k)==1)  dead_cell(i,j-1,k) = -1
!       if(iblank(i,j,k+1)==1)  dead_cell(i,j,k+1) = -1
!       if(iblank(i,j,k-1)==1)  dead_cell(i,j,k-1) = -1
!    ENDDO

    !If we want to use ghostcells in the buffer slices at zc_start-1 and zc_end+1. 
    !
    !call send_receive_slices_integer(ghostCellMark,0,nx+1,0,ny+1,zb1,zb2,2)
    !  then use the following range for the loop
    !  do k = zc_start-1,zc_end+1 
    !
    ! Find deadcells in the subdomain. The extra deadcells found in these two
    ! buffer slices will be used to determine alphaGhost for associated ghostcell in
    ! the same subdomain.  Do not exchange deadcells between subdomains, as the deadcells
    ! associated with same subdomain could be erased.
    DO k = zc_start,zc_end
    DO j = yc_start,yc_end
    DO i = 1, nx-1
       
       if(ghostCellMark(i,j,k) /= 1) cycle   ! skip cells that are not ghost cell

       !iBody = bodyNum(i,j,k)

       ! -1 means that the dead-cell is associated with a ghost-cell
       ! on the negative side of the surface ('fluid' side for a solid body)
       if(iblank(i+1,j,k)==1)  dead_cell(i+1,j,k) = -1
       if(iblank(i-1,j,k)==1)  dead_cell(i-1,j,k) = -1
       if(iblank(i,j+1,k)==1)  dead_cell(i,j+1,k) = -1
       if(iblank(i,j-1,k)==1)  dead_cell(i,j-1,k) = -1
       if(iblank(i,j,k+1)==1)  dead_cell(i,j,k+1) = -1
       if(iblank(i,j,k-1)==1)  dead_cell(i,j,k-1) = -1       
    ENDDO ! i 
    ENDDO ! j
    ENDDO ! k

    !Ye,check dc before exchange
    !write(500+lProc,*)'VARIABLES="X","Y","Z","IBLANK", "BodyNum","deadcell"'
    !write(500+lProc,*)'ZONE F=POINT, I=',nx-1,', J=',jSlices+2,', K=',kSlices+2
    !do k=zc_start-1,zc_end+1   !1,nz-1
    !do j=yc_start-1,yc_end+1
    !do i=1,nx-1
    ! write(500+lProc,'(3(2X,F12.5),3(3X,I10))')xc(i),yc(j),zc(k),iblank(i,j,k),bodyNum(i,j,k),dead_cell(i,j,k)
    !enddo
    !enddo
    !enddo
    !close(500+lProc)

    !There are dead cells in the buffer layers. Set them in the actual subdomains.
    call exchange_buffer_layer_deadcells_z()
    call exchange_buffer_layer_deadcells_y(dead_cell,0,nx+1,yb1,yb2,zb1,zb2,2)
 
    !Ye,check dc before exchange
    !write(700+lProc,*)'VARIABLES="X","Y","Z","IBLANK", "BodyNum","deadcell"'
    !write(700+lProc,*)'ZONE F=POINT, I=',nx-1,', J=',jSlices+2,', K=',kSlices+2
    !do k=zc_start-1,zc_end+1   !1,nz-1
    !do j=yc_start-1,yc_end+1
    !do i=1,nx-1
    ! write(700+lProc,'(3(2X,F12.5),3(3X,I10))')xc(i),yc(j),zc(k),iblank(i,j,k),bodyNum(i,j,k),dead_cell(i,j,k)
    !enddo
    !enddo
    !enddo
    !close(700+lProc)
    !Ye,debug
    !call MPI_BARRIER(flow_comm,ierr)
    !stop

    ! count the number of dead-cells, including those in the buffer slices
    nbdr = 0
    DO k = zc_start-1,zc_end+1    !1, nz-1
    DO j = yc_start-1,yc_end+1
    DO i = 1, nx-1

       ! unmark the deadcells right at outer boundary
       !if(i.le.1.or.i.ge.nx-1 .or. j.le.1.or.j.ge.ny-1 .or. k.le.1.or.k.ge.nz-1) then
       !   dead_cell(i,j,k) = 0
       !endif

       IF ( abs(dead_cell(i,j,k)) == 1 )   THEN
          nbdr = nbdr + 1
       ENDIF
    ENDDO ! i 
    ENDDO ! j
    ENDDO ! k

    nDead = nbdr         ! total number of fresh points

    !IF ( MOD(ntime,nmonitor) == 0 .OR. ntime == 1) &
    !    PRINT*,'lProc',lProc,'GCM_set_internal_boundary: nDead = ',nDead
    !call write_subdomain()

    !if(nDead > nDeadMax) then
    !   print*,'WARNING: More dead cells encountered. Needs to increase nDeadMax in input.dat!'
    !   stop 
    !endif

    ! Deallocate arrays pertinent to dead Cells 
    !IF (nread==0 .and. ntime >= ntime_start+1 ) THEN
    IF ( ntime >= ntime_start+1 ) THEN
       CALL GCM_DeallocateDeadCellArrays
    ENDIF ! ntime
    ! Allocate Arrays pertinent to Dead Cells
    CALL GCM_AllocateDeadCellArrays()

   ! Set appropriate values for iGhost and jGhost by doing search      
    nbdr = 0
    DO k = zc_start-1,zc_end+1    !1, nz-1
    DO j = yc_start-1,yc_end+1
    DO i = 1, nx-1
       IF ( abs(dead_cell(i,j,k)) == 1 )   THEN
          nbdr         = nbdr + 1
          iDead(nbdr) = i
          jDead(nbdr) = j
          kDead(nbdr) = k
       ENDIF 
    ENDDO ! i 
    ENDDO ! j
    ENDDO ! k

    alphaGhost(:) = 0.0_CGREAL

    nbdr = 0
    ! project each fresh node onto the boundary       
    DO n = 1, nDead

       iG = iDead(n) 
       jG = jDead(n)
       kG = kDead(n)

       iBody = bodyNum(iG,jG,kG)

       if(iBody==0) then
          print*,'GCM_set_deadcells: iBody = 0!'
          print*,iG,jG,kG,xc(iG),yc(jG),zc(kG)
          call write_subdomain()
          stop
       endif

       iDeadCellIndex(n) = -1
       jDeadCellIndex(n) = -1
       kDeadCellIndex(n) = -1

       xDC = xc(iG)
       yDC = yc(jG)
       zDC = zc(kG)

       !if(is_BIset(iG,jG,kG)==1) then
          !Use the BI information from the saved table
       !   xBI  = xBItable(iG,jG,kG)
       !   yBI  = yBItable(iG,jG,kG)
       !   zBI  = zBItable(iG,jG,kG)
       !   closestElementDead(n)  = cELtable(iG,jG,kG)

       !else  
       !   nbdr = nbdr + 1

          !call the body intercept subroutine
          CALL calc_bodyIntercept_Unstruc(iBody, iG,jG,kG, xBI,yBI,zBI, &
                                          closestElementDead(n), distBIElem, dotNorm)

          IF( body_dim(iBody) == 2 )  zBI = zDC ! Enforce for 2D Body
       !endif ! is_BIset

       ! Extract coordinates of Body Intercept.
       ! Note that, for a membrane structure, the true BI point 
       ! is saved.
       xBodyInterceptDead(n) = xBI
       yBodyInterceptDead(n) = yBI
       zBodyInterceptDead(n) = zBI

       !-- Get length of normal
       dsIntercept = SQRT( (xDC-xBI)**2 + (yDC-yBI)**2 + (zDC-zBI)**2 )

       slopeX = (xBI - xDC)/dsIntercept
       slopeY = (yBI - yDC)/dsIntercept
       slopeZ = (zBI - zDC)/dsIntercept

       ! intercept normal (solid: pointing out of the body)
       ! intercept normal (membrane: pointing into the body)
       xNormInterceptDead(n) = -slopeX
       yNormInterceptDead(n) = -slopeY
       zNormInterceptDead(n) = -slopeZ

       ! set up virtual intercept for membrane-type structure
       ! by reversing the direction of the intercept
       if((unstruc_surface_type(iBody)==MEMBRANE)) then

          !Ye,debug=====================================
          !xNormInterceptDead(n) = -xNormInterceptDead(n)
          !yNormInterceptDead(n) = -yNormInterceptDead(n)
          !zNormInterceptDead(n) = -zNormInterceptDead(n)
          !=============================================

          slopeX = - slopeX
          slopeY = - slopeY
          slopeZ = - slopeZ

          !Ye,debug
          !if ( (iG.eq.144).and.(jG.eq.102).and. &
          !     (kG.eq.68) )then
          !   write(500+lProc,*)ntime,iter_FSI,lProc,n
          !endif

          !Ye,debug,ture BI position
          !if ( (ntime.eq.6625).and.(lProc.eq.2).and.(n.eq.10508) )then
          !   write(504,'(4(i8,1X))')ntime,n,iG,jG,kG
          !   write(504,'(2(F12.7,1X))')dsIntercept,membrane_tkns/2.0_CGREAL
          !   write(504,'(3(F12.7,1X))')xBI,yBI,zBI
          !endif

          !Ye.debug,add the situation when dc cell locates outside
          !of real membrane due to iblank_smooth
          !=================
          IF (dsIntercept .ge.membrane_tkns/2.0_CGREAL) THEN
             xBI = xBI + (dsIntercept + 0.1 * (dsIntercept - &
                   membrane_tkns/2.0_CGREAL)) * slopeX
             yBI = yBI + (dsIntercept + 0.1 * (dsIntercept - &
                   membrane_tkns/2.0_CGREAL)) * slopeY
             zBI = zBI + (dsIntercept + 0.1 * (dsIntercept - &
                   membrane_tkns/2.0_CGREAL)) * slopeZ
          ELSE
             xBI = xBI + membrane_tkns/2.0_CGREAL * slopeX
             yBI = yBI + membrane_tkns/2.0_CGREAL * slopeY
             zBI = zBI + membrane_tkns/2.0_CGREAL * slopeZ
          ENDIF

          ! check if the image point will be outside the subdomain.  If so,
          ! scale the BI along the normal and reduce the distance.
          zIP = zDC + (zBI-zDC)*probeLengthNormalizedD
          !Ye,debug,virtual BI
          !if ( (ntime.eq.6625).and.(lProc.eq.2).and.(n.eq.10508) )then
          !   write(505,'(4(i8,1X))')ntime,n,iG,jG,kG
          !   write(505,'(2(i8,1X))')zb1,zb2
          !   write(505,'(3(F12.7,1X))')zIP,zc(zb1),zc(zb2)
          !   write(505,'(3(F12.7,1X))')xBI,yBI,zBI
          !endif

          !Ye,debug, distance btw DC & virtual BI
          !ds1 = SQRT( (xDC-xBI)**2 + (yDC-yBI)**2 + (zDC-zBI)**2 )
 
          !if(zIP .lt. zc(zb1) .or. zIP .gt. zc(zb2)) then 
             !Cancel the previous shift and add a new/shorter distance.
             !Ye,debug
          !   if( (ntime.eq.6625).and.(lProc.eq.2).and.(n.eq.10508) )then
          !     temp = -0.5d0 * ds1
          !   else
          !     temp = -membrane_tkns/2.0_CGREAL + (0.5_CGREAL-0.01d0)*dz(kG)
          !   endif
 
          !   xBI = xBI + temp * slopeX
          !   yBI = yBI + temp * slopeY
          !   zBI = zBI + temp * slopeZ
          !endif

          !Ye,NEW way to shift BI
          DO WHILE ((yIP .lt. yc(yb1)) .or. (yIP .gt. yc(yb2)))
             !Ye, distance btw DC & virtual BI
             ds1 = SQRT( (xDC-xBI)**2 + (yDC-yBI)**2 + (zDC-zBI)**2 )
             temp = -0.5d0 * ds1

             xBI = xBI + temp * slopeX
             yBI = yBI + temp * slopeY
             zBI = zBI + temp * slopeZ

             yIP = yDC + (yBI-yDC)*probeLengthNormalizedD
          ENDDO!DO WHILE
          DO WHILE ((zIP .lt. zc(zb1)) .or. (zIP .gt. zc(zb2)))
             !Ye, distance btw DC & virtual BI
             ds1 = SQRT( (xDC-xBI)**2 + (yDC-yBI)**2 + (zDC-zBI)**2 )
             temp = -0.5d0 * ds1

             xBI = xBI + temp * slopeX
             yBI = yBI + temp * slopeY
             zBI = zBI + temp * slopeZ         

             zIP = zDC + (zBI-zDC)*probeLengthNormalizedD
          ENDDO!DO WHILE

          dsIntercept = SQRT( (xDC-xBI)**2 + (yDC-yBI)**2 + (zDC-zBI)**2 )
          !Ye,debug,shift BI because zIp outside the subdomain
          !if ( (ntime.eq.6625).and.(lProc.eq.2).and.(n.eq.10508) )then
          !   write(506,'(4(i8,1X))')ntime,n,iG,jG,kG
          !   write(506,'(3(F12.7,1X))')slopeX,slopeY,slopeZ
          !   write(506,'(1(F12.7,1X))')temp,dz(kG)
          !   write(506,'(3(F12.7,1X))')xBI,yBI,zBI
          !endif
       endif!MEMBRANE

       ! set up alphaGhost for associated ghost-cell
       ! Only handle ghostcells that are in the same subdomain
       if(abs(ghostcellMark(iG+1,jG,kG))==1) then
          if(kG.lt.zc_start .or. kG.gt.zc_end) cycle  !skip ghost cells in the buffer slices

          temp = dsIntercept / dxc(iG+1)
          nG   = ghostCellIndex(iG+1,jG,kG)
          if(nG==0) print*,iG,jG,kG,iG+1,jG,kG,ghostcellMark(iG+1,jG,kG)
          !alphaGhost(nG) = max(alphaGhost(nG), temp)
          alphaGhost(nG) = alphaGhost(nG) + temp**2
       endif
       if(abs(ghostcellMark(iG-1,jG,kG))==1) then 
          if(kG.lt.zc_start .or. kG.gt.zc_end) cycle  !skip ghost cells in the buffer slices

          temp = dsIntercept / dxc(iG)
          nG   = ghostCellIndex(iG-1,jG,kG)
          if(nG==0) print*,iG,jG,kG,iG-1,jG,kG,ghostcellMark(iG-1,jG,kG)
          !alphaGhost(nG) = max(alphaGhost(nG), temp)
          alphaGhost(nG) = alphaGhost(nG) + temp**2
       endif

       if(abs(ghostcellMark(iG,jG+1,kG))==1) then
          if(kG.lt.zc_start .or. kG.gt.zc_end) cycle  !skip ghost cells in the buffer slices

          temp = dsIntercept / dyc(jG+1)
          nG   = ghostCellIndex(iG,jG+1,kG)
          if(nG==0) print*,iG,jG,kG,iG,jG+1,kG,ghostcellMark(iG,jG+1,kG)
          !alphaGhost(nG) = max(alphaGhost(nG), temp)
          alphaGhost(nG) = alphaGhost(nG) + temp**2
       endif
       if(abs(ghostcellMark(iG,jG-1,kG))==1) then
          if(kG.lt.zc_start .or. kG.gt.zc_end) cycle  !skip ghost cells in the buffer slices

          temp = dsIntercept / dyc(jG)
          nG   = ghostCellIndex(iG,jG-1,kG)
          if(nG==0) print*,iG,jG,kG,iG,jG-1,kG,ghostcellMark(iG,jG-1,kG)
          !alphaGhost(nG) = max(alphaGhost(nG), temp)
          alphaGhost(nG) = alphaGhost(nG) + temp**2
       endif

       ! Only handle ghostcells that are in the same subdomain
       if(abs(ghostcellMark(iG,jG,kG+1))==1) then          
          if(kG+1.lt.zc_start .or. kG+1.gt.zc_end) cycle  !skip ghost cells in the buffer slices

          temp = dsIntercept / dzc(kG+1)
          nG   = ghostCellIndex(iG,jG,kG+1)
          if(nG==0) print*,iG,jG,kG,iG,jG,kG+1,ghostcellMark(iG,jG,kG+1)
          !alphaGhost(nG) = max(alphaGhost(nG), temp)
          alphaGhost(nG) = alphaGhost(nG) + temp**2          
       endif

       if(abs(ghostcellMark(iG,jG,kG-1))==1) then   
          if(kG-1.lt.zc_start .or. kG-1.gt.zc_end) cycle  !skip ghost cells in the buffer slices

          temp = dsIntercept / dzc(kG)
          nG   = ghostCellIndex(iG,jG,kG-1)
          if(nG==0) print*,iG,jG,kG,iG,jG,kG-1,ghostcellMark(iG,jG,kG-1)
          !alphaGhost(nG) = max(alphaGhost(nG), temp)
          alphaGhost(nG) = alphaGhost(nG) + temp**2
       endif

        ! compute the image point
        !xIP = xDC + slopeX*dsIntercept*probeLengthNormalizedD
        !yIP = yDC + slopeY*dsIntercept*probeLengthNormalizedD
        !zIP = zDC + slopeZ*dsIntercept*probeLengthNormalizedD

        ! in case dsIntercept is 0, the following lines are better
        xIP = xDC + (xBI-xDC)*probeLengthNormalizedD
        yIP = yDC + (yBI-yDC)*probeLengthNormalizedD
        zIP = zDC + (zBI-zDC)*probeLengthNormalizedD

        probeLengthDead(n) = dsIntercept*probeLengthNormalizedD

        !Ye,debug, shifr IP since IP outside the subdomain
        !if ( (ntime.eq.6625).and.(lProc.eq.2).and.(n.eq.10508) )then
        !   write(507,'(4(i8,1X))')ntime,n,iG,jG,kG
        !   write(507,'(3(F12.7,1X))')xIP,yIP,zIP
        !   write(507,'(2(F12.7,1X))')zc(zb1),zc(zb2)
        !endif

        !if the IP is outside the subdomain, adjust it back
        if(yIP .lt. yc(yb1) .or. yIP .gt. yc(yb2)) then
           if(abs(yIP-yc(yb1)) <  abs(yIP-yc(yb2))) then
              yIP = yc(yb1)
           else
              yIP = yc(yb2)
           endif
           probeLengthDead(n) = sqrt( (xDC-xIP)**2 + (yDC-yIP)**2 + (zDC-zIP)**2 )
        endif
        if(zIP .lt. zc(zb1) .or. zIP .gt. zc(zb2)) then 
           if(abs(zIP-zc(zb1)) <  abs(zIP-zc(zb2))) then
              zIP = zc(zb1)
           else
              zIP = zc(zb2)
           endif
           probeLengthDead(n) = sqrt( (xDC-xIP)**2 + (yDC-yIP)**2 + (zDC-zIP)**2 )
        endif

        xyzDCIP(iG,jG,kG,0) = xIP
        xyzDCIP(iG,jG,kG,1) = yIP
        xyzDCIP(iG,jG,kG,2) = zIP

        ! Find the lower left grid point to Image Point in Physical domai
        maxDelta = MAX(dx(iG),dy(jG),dz(kG))
 
        iRange  = 2*NINT(maxDelta/dx(iG) + 0.5_CGREAL)
        jRange  = 2*NINT(maxDelta/dy(jG) + 0.5_CGREAL)
        kRange  = 2*NINT(maxDelta/dz(kG) + 0.5_CGREAL)

        iMin = iG-iRange
        iMax = iG+iRange
        iMin = MAX(iMin,0)    ! note image point is allowed to be between xc(0) and x(1)
        iMax = MIN(iMax,nx)   ! note image point is allowed to be between x(nx) and xc(nx)

        DO i = iMin,iMax 
           IF ( ( xc(i) <= xIP .AND. xc(i+1) >= xIP ) ) THEN
              iDeadCellIndex(n) = i
           ENDIF ! xc
        ENDDO ! i

        jMin = jG-jRange
        jMax = jG+jRange
        jMin = MAX(jMin,yc_start-2)
        jMax = MIN(jMax,yc_end  +1)

        DO j = jMin,jMax
           IF ( ( yc(j) <= yIP .AND. yc(j+1) >= yIP ) ) THEN
              jDeadCellIndex(n) = j
           ENDIF ! xc
        ENDDO ! j

        kMin = kG-kRange
        kMax = kG+kRange

        kMin = MAX(kMin,zc_start-2)
        kMax = MIN(kMax,zc_end  +1)

        DO k = kMin,kMax
           IF ( ( zc(k) <= zIP .AND. zc(k+1) >= zIP ) ) THEN
              kDeadCellIndex(n) = k
           ENDIF ! xc
        ENDDO ! k

     !IF((lProc.eq.5).AND.(iG.eq.100).AND.(jG.eq.81).AND.(kG.eq.53))THEN
     !  write(82,*)'INFO:',ntime,iter_FSI
     !  write(82,*)dead_cell(100,81,53)
     !  write(82,'(3F12.4)')xDC, yDC, zDC
     !  write(82,'(3F12.4)')xBI, yBI, zBI
     !  write(82,'(3F12.4)')xIP, yIP, zIP
     !  write(82,*)iDeadCellIndex(n),jDeadCellIndex(n),kDeadCellIndex(n)
     !ENDIF

        IF ( iDeadCellIndex(n) == -1 .OR. jDeadCellIndex(n) ==-1  &
                                     .OR. kDeadCellIndex(n) ==-1) THEN
           PRINT*,'Failed to find the hosting box for the image point of a dead cell'
           !PRINT*,ntime
           !PRINT*,n
           PRINT*,iG, jG, kG
           PRINT*,iRange, jRange, kRange
           PRINT*,iBody,closestElementDead(n)
           PRINT*,iDeadCellIndex(n),jDeadCellIndex(n),kDeadCellIndex(n)
           PRINT*,xdc,ydc,zdc
           PRINT*,xBI,yBI,zBI
           PRINT*,xIP,yIP,zIP
           PRINT*,lproc,kMin,xc(iMin),yc(jMin),zc(kMin)
           PRINT*,lProc,kMax,xc(iMax),yc(jMax),zc(kMax)
           PRINT*,dsIntercept
    
           !PRINT*,jMin, jMax, yc(jMin), yc(jMax)
           PRINT*,'Aborting Run; write dump ...'
           CALL write_dump()
           call sleep(15)
           call MPI_BARRIER(flow_comm,ierr)
           STOP
        ENDIF

        !Ye,Mar==========================
        !Find DCs with pure solid node in interpolation
        DO iRow = 1, iRowMax
           ii = iDeadCellIndex(n) + incI(iRow)
           jj = jDeadCellIndex(n) + incJ(iRow)
           kk = kDeadCellIndex(n) + incK(iRow)
           IF( (iblank(ii,jj,kk).eq.1).AND.  &
                 (dead_cell(ii,jj,kk).ne.-1) )THEN
             if(jG<yc_start-1 .or. jG>yc_end+1) cycle
             if(kG<zc_start-1 .or. kG>zc_end+1) cycle
             ps1(iG,jG,kG) = 10
             xyzDCIP(iG,jG,kG,0) = xIP
             xyzDCIP(iG,jG,kG,1) = yIP
             xyzDCIP(iG,jG,kG,2) = zIP
           ENDIF
        ENDDO

        !Mark the nearby GCs
        IF (ps1(iG,jG,kG).eq.10)THEN
           if ( (iG-1.ge.0).and.(iG-1.le.nx).and. &
                (ghostCellMark(iG-1,jG,kG).eq.1) )then
              ps2(iG-1,jG,kG)=10
           endif
           if ( (iG+1.ge.0).and.(iG+1.le.nx).and. &
                (ghostCellMark(iG+1,jG,kG).eq.1) )then
              ps2(iG+1,jG,kG)=10
           endif

           if ( (jG-1.ge.0).and.(jG-1.le.ny).and. &
                (ghostCellMark(iG,jG-1,kG).eq.1) )then
              ps2(iG,jG-1,kG)=10
           endif
           if ( (jG+1.ge.0).and.(jG+1.le.ny).and. &
                (ghostCellMark(iG,jG+1,kG).eq.1) )then
              ps2(iG,jG+1,kG)=10
           endif

           if ( (kG-1.ge.zc_start).and.(kG-1.le.zc_end).and. &
                (ghostCellMark(iG,jG,kG-1).eq.1) )then
              ps2(iG,jG,kG-1)=10
           endif
           if ( (kG+1.ge.zc_start).and.(kG+1.le.zc_end).and. &
                (ghostCellMark(iG,jG,kG+1).eq.1) )then
              ps2(iG,jG,kG+1)=10
           endif
        ENDIF
        !================================

        ! Perform bilinear interpolation
        iCIndx = iDeadCellIndex(n)
        jCIndx = jDeadCellIndex(n)
        kCIndx = kDeadCellIndex(n)

        xBIN = xNormInterceptDead(n)
        yBIN = yNormInterceptDead(n)
        zBIN = zNormInterceptDead(n)

        CALL GCM_Calc_vanMatrixDN( iG, jG, kG, iCIndx, jCIndx, kCIndx,             &
                                   xIP, yIP, zIP, xBI, yBI, zBI, xBIN, yBIN, zBIN, &
                                   coeffGCMDirc, coeffGCMNeum      )

        coeffGCMDeadD(1:iRowMax,n) = coeffGCMDirc(1:iRowMax)
        coeffGCMDeadN(1:iRowMax,n) = coeffGCMNeum(1:iRowMax)

        !write(200,'(I8,1X,12I5)')n,iG,jG,kG,iCIndx,jCIndx,kCIndx
        !write(200,'(I8,1X,12F15.5)')n,xBI,yBI,zBI,xIP,yIP,zIP

        !Ye,debug, MARKER1/2/3/4
        !IF((ntime.eq.4803).AND.(iG.eq.81).AND.(jG.eq.82).AND.(kG.eq.51))THEN
        !   WRITE(75,*)ntime,iter_FSI,lProc,n,closestElementDead(n)
        !   WRITE(75,*)bodyNum(iG,jG,kG)
        !   WRITE(75,*)xBodyInterceptDead(n),yBodyInterceptDead(n),  &
        !              zBodyInterceptDead(n)
        !   WRITE(75,*)xBI,yBI,zBI
        !   WRITE(75,*)xIP,yIP,zIP
        !   WRITE(75,*)iCIndx,jCIndx,kCIndx
        !   WRITE(75,*)'========'
        !ENDIF
        
     ENDDO ! n 
     !print*,lProc,'BItable not set for deadcells: nbdr/nDead=',nbdr,nDead

     ! set up alphaGhost for each ghost-cell
     DO n = 1, nGhost
        iG = iGhost(n) 
        jG = jGhost(n)
        kG = kGhost(n)

        !alphaGhost(n) = 1.0_CGREAL - alphaGhost(n)
        alphaGhost(n) = 1.0_CGREAL - sqrt(alphaGhost(n))
        !alphaGhost(n) = 1.0_CGREAL  !Either completely from interpolation or from NS Eq., 
        !                            !depending on value of mix_gc_form

        !write(114+lProc,'(I5,6F12.5)')n,xc(iG),yc(jG),zc(iG),probeLength(n)/probeLengthNormalized,alphaGhost(n)

     ENDDO

  END SUBROUTINE GCM_set_deadcells
!-----------------------------------------------------------------------
  SUBROUTINE exchange_buffer_layer_deadcells_z()

  !set additional dead cells identified by the neighbour subdomains.

    USE flow_parameters
    USE boundary_arrays
    USE MPI_module

    IMPLICIT NONE
    
    INTEGER :: i,j,k,ndata, l1, l2

    i_lbuff = 0
    i_rbuff = 0

    ! number of data to be exchanged.
    ndata = (nx+2)*(yb2-yb1+1)*2
    !print*, 'ndata = ',ndata

    !Ye,debug
    !if (lProc.eq.13)then
    !   write(80,*)kProc,lProc_bottom,ndata
    !   write(80,*)yb1,yb2,zb1
    !   write(80,*)'---'
    !endif

    !Send to the left; provide the 1st element location and number of data
    l1 = zb1+1
    l2 = l1+1
    if (kProc > 0)  &
         call MPI_ISEND(dead_cell(0,yb1,l1), ndata, MPI_INTEGER1, lProc_bottom, &
                       toleft, FLOW_COMM, isend_rq_toleft, ierr)

    !Ye,test
    !call MPI_BARRIER(flow_comm,ierr)
    !if (lProc.eq.proc_m)then
    !   write(*,*)'DC send to bottom OK!'
    !endif
    !if (lProc.eq.13)then
    !   write(90,*)lProc_bottom
    !   write(90,*)yb1,yb2,l1,l2
    !   write(90,*)'---'
    !endif

    !Send to the right
    l2 = zb2-1
    l1 = l2-1    
    if (kProc < nProcZ-1)  &
         call MPI_ISEND(dead_cell(0,yb1,l1), ndata, MPI_INTEGER1, lProc_top, &
                       torght, FLOW_COMM, isend_rq_torght, ierr)

    !Ye,test
    !call MPI_BARRIER(flow_comm,ierr)
    !if (lProc.eq.proc_m)then
    !   write(*,*)'DC send to top OK!'
    !endif

    !Receiving from right; provide the 1st element location and number of data      
    if (kProc < nProcZ-1)  &
         call MPI_IRECV(i_rbuff(0,yb1,1), ndata, MPI_INTEGER1, lProc_top, &
                       toleft, FLOW_COMM, irecv_rq_fromrght, ierr)

    !Ye,test
    !call MPI_BARRIER(flow_comm,ierr)
    !if (lProc.eq.proc_m)then
    !   write(*,*)'DC receive from top OK!'
    !endif
    !if (lProc.eq.13)then
    !   write(91,*)lProc_top
    !   write(91,*)yb1,yb2,l1,l2
    !   write(91,*)'---'
    !endif

    !Receiving from left 
    if (kProc > 0) &
         call MPI_IRECV(i_lbuff(0,yb1,1), ndata, MPI_INTEGER1, lProc_bottom, &
                       torght, FLOW_COMM, irecv_rq_fromleft, ierr)
   
    !Ye,test
    !call MPI_BARRIER(flow_comm,ierr)
    !if (lProc.eq.proc_m)then
    !   write(*,*)'DC receive from bottom OK!'
    !endif

    !Wait for sending and receiving to complete
    if (kProc > 0) &
         call MPI_WAIT(isend_rq_toleft, isend_stat_tl, ierr)
    if (ierr /= 0) print*,'Warning: Wait toleft says', ierr

    if (kProc < nProcZ-1)  &
         call MPI_WAIT(isend_rq_torght, isend_stat_tr, ierr)
    if (ierr /= 0) print*,'Warning: Wait torght says', ierr

    if (kProc < nProcZ-1)  &
         call MPI_WAIT(irecv_rq_fromrght, irecv_stat_fr, ierr)
    if (ierr /= 0) print*,'Warning: Wait toleft (fromrght) says', ierr
    
    if (kProc > 0) &
         call MPI_WAIT(irecv_rq_fromleft, irecv_stat_fl, ierr)
    if (ierr /= 0) print*,'Warning: Wait torght (fromleft) says', ierr

    !print*,'Min  dc lbuff: ',MINVAL(dead_cell(1:nx-1,1:ny-1,zb1:zb1+1))
    !print*,'Min  dc rbuff: ',MINVAL(dead_cell(1:nx-1,1:ny-1,zb2-1:zb2))
    !print*,'Min i_lbuff: ',MINVAL(i_lbuff(1:nx-1,1:ny-1,1:2))
    !print*,'Min i_rbuff: ',MINVAL(i_rbuff(1:nx-1,1:ny-1,1:2))

    !Compare set dead cells on both sides
    do k=1,2
    do j=yb1,yb2
    do i=0,nx
       if(i_lbuff(i,j,k) == -1) then
          !print*,'i_lbuff: ', i,j,k
          dead_cell(i,j,zc_start+k-2) = -1
       endif
       if(i_rbuff(i,j,k) == -1) then
          !print*,'i_rbuff: ', i,j,k
          dead_cell(i,j,zc_end  -1+k) = -1
       endif
    enddo
    enddo
    enddo

  END SUBROUTINE exchange_buffer_layer_deadcells_z
!--------------------------------------------------------
!--------------------------------------------------------------
  SUBROUTINE smooth_iblank()
!
! This subroutine identify the isolated or nearly isolated fluid cells,
! then mark them as iblank cells

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE grid_arrays
    USE boundary_arrays
    USE gcm_arrays
    USE MPI_module
    
    IMPLICIT NONE

    INTEGER :: i,j,k, ii, jj, kk, i2,j2,k2
    INTEGER :: iloop, nloop, ibk_sum, ibody, nbdr

    nloop = 5
    nbdr  = 0

    DO iloop = 1,nloop

      DO k = zc_start, zc_end
      DO j = yc_start, yc_end
      DO i = 1, nx-1

         IF (iblank(i,j,k) == 1) cycle !skip iblank cells
            
         ibk_sum = 0

         DO kk = -1,1
         DO jj = -1,1
         DO ii = -1,1

            i2 = i + ii
            j2 = j + jj
            k2 = k + kk
            if(abs(ii)+abs(jj)+abs(kk) /= 1) cycle  ! skip edge and corner points
            
            if(iblank(i2,j2,k2) == 1) then 
               ibk_sum = ibk_sum + 1
               if(bodyNum(i2,j2,k2) > 0) ibody = bodyNum(i2,j2,k2)
            endif
         ENDDO     !ii
         ENDDO     !jj
         ENDDO     !kk

         if(ibk_sum >= 5) then
            iblank(i,j,k) = 1 ! mark as iblank cell
            bodyNum(i,j,k) = ibody
            nbdr = nbdr + 1
            !Ye
            !write(88,*)'INFO:',lProc, ntime, i, j, k
         endif
      ENDDO ! i
      ENDDO ! j
      ENDDO ! k

      !Ye
      !IF (lProc.eq.13) THEN
      !   write(89,*)'INFO:', yb1,yb2,zb1,zb2
      !   write(89,*)'INFO:', lProc_left,lProc_right
      !   write(89,*)'INFO:', lProc_bottom,lProc_top
      !   write(89,*)'INFO:', jProc,kProc
      !   write(89,*)'-----'
      !ENDIF

      call send_receive_slices_integer_y(iblank, 0,nx+1,yb1,yb2,zb1,zb2,1) 
      call send_receive_slices_integer_y(bodyNum,0,nx+1,yb1,yb2,zb1,zb2,1)

      call send_receive_slices_integer_z(iblank, 0,nx+1,yb1,yb2,zb1,zb2,1)
      call send_receive_slices_integer_z(bodyNum,0,nx+1,yb1,yb2,zb1,zb2,1)

   ENDDO ! end iloop

   IF (MOD(ntime,nmonitor) == 0 .and. nbdr>0) &
      print*,'3D-Smooth_iblank: number of isolated cells removed:',nbdr

   END SUBROUTINE smooth_iblank
!-----------------------------------------------------
!-----------------------------------------------------------------------
  SUBROUTINE exchange_buffer_layer_deadcells_y(f,i1,i2,j1,j2,k1,k2,ks)

  !set additional dead cells identified by the neighbour subdomains.

    USE flow_parameters
    USE boundary_arrays
    USE MPI_module

    IMPLICIT NONE
    
    INTEGER,INTENT (IN) :: i1,i2,j1,j2,k1,k2,ks
    INTEGER(KIND=INT_K), DIMENSION(i1:i2,j1:j2,k1:k2),INTENT (IN OUT) :: f

    INTEGER :: i,j,k,ndatay,l1,l2
       
    !Data needs to be continuous in memory for MPI passing.
    !Use a temporary array to hold the matrix transpose; only ghost slices are
    !needed.
    INTEGER(KIND=INT_K), DIMENSION(i1:i2,k1:k2,0:ks-1) :: g_tr1, g_tr2
    INTEGER(KIND=INT_K), DIMENSION(i1:i2,k1:k2,0:ks-1) :: g_tr3, g_tr4

    g_tr1 = 0
    g_tr2 = 0
    g_tr3 = 0
    g_tr4 = 0

    !left buffer transpose
    if (jProc > lProc_leftmost)then
    l1 = j1+1
    l2 = l1+ks-1
    do i=i1,i2
    do j=l1,l2
    do k=k1,k2
       g_tr1(i,k,j-l1) = f(i,j,k)
    enddo
    enddo
    enddo
    endif

    !right buffer transpose
    if (jProc < lProc_rightmost)then
    l2 = j2-1
    l1 = l2-ks+1
    do i=i1,i2
    do j=l1,l2
    do k=k1,k2
       g_tr2(i,k,j-l1) = f(i,j,k)
    enddo
    enddo
    enddo
    endif

    ! number of data to be exchanged, where ks is either 1 or 2.
    ndatay = (i2-i1+1)*ks*(k2-k1+1)
    !print*, 'ndata = ',ndata

    !Send to the left; provide the 1st element location and number of data
    if (jProc > 0)  &
         call MPI_ISEND(g_tr1(i1,k1,0), ndatay, MPI_INTEGER1, lProc_left, &
                       toleft, FLOW_COMM, isend_rq_toleft, ierr)
    if (ierr /= 0) print*,'Warning: Isend toleft says', ierr

    !Ye,test
    !call MPI_BARRIER(flow_comm,ierr)
    !if (lProc.eq.proc_m)then
    !   write(*,*)'DC send to LEFT OK!'
    !endif
    !if (lProc.eq.13)then
    !   write(90,*)lProc_left,ndatay
    !   write(90,*)j1,j2,l1,l2
    !   write(90,*)'---'
    !endif

    !Send to the right    
    if (jProc < nProcY-1)  &
         call MPI_ISEND(g_tr2(i1,k1,0), ndatay, MPI_INTEGER1, lProc_right, &
                       torght, FLOW_COMM, isend_rq_torght, ierr)
    if (ierr /= 0) print*,'Warning: Isend torght says', ierr

    !Ye,test
    !call MPI_BARRIER(flow_comm,ierr)
    !if (lProc.eq.proc_m)then
    !   write(*,*)'DC send to RIGHT OK!'
    !endif
    !if (lProc.eq.13)then
    !   write(91,*)lProc_right
    !   write(91,*)yb1,yb2,l1,l2
    !   write(91,*)'---'
    !endif

    !Receiving from right; provide the 1st element location and number of data
    if (jProc < nProcY-1)  &
         call MPI_IRECV(g_tr3(i1,k1,0), ndatay, MPI_INTEGER1, lProc_right, &
                       toleft, FLOW_COMM, irecv_rq_fromrght, ierr)
    if (ierr /= 0) print*,'Warning: Irecv toleft (fromrght) says', ierr

    !Ye,test
    !call MPI_BARRIER(flow_comm,ierr)
    !if (lProc.eq.proc_m)then
    !   write(*,*)'DC receive from RIGHT OK!'
    !endif
    !if (lProc.eq.13)then
    !   write(92,*)lProc_right
    !   write(92,*)j1,j2,l1,l2
    !   write(92,*)'---'
    !endif

    !Receiving from left
    if (jProc > 0) &
         call MPI_IRECV(g_tr4(i1,k1,0), ndatay, MPI_INTEGER1, lProc_left, &
                       torght, FLOW_COMM, irecv_rq_fromleft, ierr)
    if (ierr /= 0) print*,'Warning: Irecv torght (fromleft) says', ierr

    !Ye,test
    !call MPI_BARRIER(flow_comm,ierr)
    !if (lProc.eq.proc_m)then
    !   write(*,*)'DC receive from LEFT OK!'
    !endif

    !Wait for sending and receiving to complete
    if (jProc > 0) &
         call MPI_WAIT(isend_rq_toleft, isend_stat_tl, ierr)
    if (ierr /= 0) print*,'Warning: Wait toleft says', ierr

    if (jProc < nProcY-1)  &
         call MPI_WAIT(isend_rq_torght, isend_stat_tr, ierr)
    if (ierr /= 0) print*,'Warning: Wait torght says', ierr

    if (jProc < nProcY-1)  &
         call MPI_WAIT(irecv_rq_fromrght, irecv_stat_fr, ierr)
    if (ierr /= 0) print*,'Warning: Wait toleft (fromrght) says', ierr
    
    if (jProc > 0) &
         call MPI_WAIT(irecv_rq_fromleft, irecv_stat_fl, ierr)
    if (ierr /= 0) print*,'Warning: Wait torght (fromleft) says', ierr

    !print*,'Min  dc lbuff: ',MINVAL(dead_cell(1:nx-1,1:ny-1,zb1:zb1+1))
    !print*,'Min  dc rbuff: ',MINVAL(dead_cell(1:nx-1,1:ny-1,zb2-1:zb2))
    !print*,'Min i_lbuff: ',MINVAL(i_lbuff(1:nx-1,1:ny-1,1:2))
    !print*,'Min i_rbuff: ',MINVAL(i_rbuff(1:nx-1,1:ny-1,1:2))

    !Compare and set dead cells on left and right side
    do i=i1,i2
    do k=k1,k2
    do j=0,ks-1
       !Left
       if (jProc > 0)then
       if (g_tr4(i,k,j) == -1)then
          f(i,j1+3-ks+j,k) = -1
       endif
       endif 
       !Right
       if (jProc < nProcY-1)then
       if (g_tr3(i,k,j) == -1)then
          f(i,j2-2+j,k) = -1
       endif
       endif
    enddo
    enddo
    enddo

  END SUBROUTINE exchange_buffer_layer_deadcells_y
!--------------------------------------------------------
