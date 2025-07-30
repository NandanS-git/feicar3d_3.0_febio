!------------------------------------------------------------------------------
   SUBROUTINE set_iblank_body_fast()
!------------------------------------------------------------------------------
    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE MPI_module
    
    IMPLICIT NONE

    INTEGER, PARAMETER :: UNDECIDED_VALUE = 100

    INTEGER           :: i,j,k,m,iBody,m_theta,m_phi,np,num_fresh,num_dead,inside
    INTEGER           :: cMarker,cElement,nBodySelect
    INTEGER           :: ii,jj,kk,k1,k2,j1,j2
    INTEGER           :: iMin,iMax,jMin,jMax,kMin,kMax
    INTEGER           :: iCellMin,iCellMax,jCellMin,jCellMax,kCellMin,kCellMax
    INTEGER           :: scanRange,sumIblank,tmpIblank
    INTEGER           :: iBeg,jBeg,kBeg
    INTEGER           :: nLevelRefine,nSubDiv
    INTEGER           :: iclock1,iclock2,iclock3,iclock4,iclock5,iclock6,clock_rate
    INTEGER           :: iVert,mVert1,mVert2,mVert3,nVertTot
    INTEGER           :: iTri,jTri,kTri
    INTEGER           :: cElement2,nBodySelect2,tmpIblank2,Nm_tkns
    !INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: iblankUndecided, iblankTemp
    REAL(KIND=CGREAL) :: rad,theta,phi,x_prime,y_prime,z_prime
    REAL(KIND=CGREAL) :: distElement,dMin,dotNorm,dMinUnstruc
    REAL(KIND=CGREAL) :: xBM,yBM,zBM,xP,yP,zP,xM,yM,zM
    REAL(KIND=CGREAL) :: xTemp,yTemp,zTemp
    REAL(KIND=CGREAL) :: lenGridMin,lenElemMax,ratioElemGrid,checkDist
    REAL(KIND=CGREAL) :: rSubDivInv,uTri,vTri,wTri
    REAL(KIND=CGREAL) :: dMinUnstruc0,dMinUnstruc1,dMinUnstruc2,dMinUnstruc3
    REAL(KIND=CGREAL) :: dotNorm1,dotNorm2,dotNorm3
    REAL(KIND=CGREAL) :: xBoundMax,xBoundMin,yBoundMax,yBoundMin,&
                         zBoundMax,zBoundMin
    REAL(KIND=CGREAL) :: xElem,yElem,zElem
    REAL(KIND=CGREAL) :: distMin1,distMin2
    REAL(KIND=CGREAL) :: distxBoundMin,distxBoundMax
    REAL(KIND=CGREAL) :: distyBoundMin,distyBoundMax
    REAL(KIND=CGREAL) :: distzBoundMin,distzBoundMax
    REAL(KIND=CGREAL), DIMENSION(3) :: delta,lenElem,lenGrid 
    REAL(KIND=CGREAL), DIMENSION(3) :: xVert,yVert,zVert
    REAL(KIND=CGREAL) :: rtemp
! ============================================================================
!   Initialize values
! ============================================================================

    num_fresh  = 0
    num_dead   = 0


    !fresh_cell = 0 

    ! Initialize the iblank variables
    bodyNum        (:,:,:) = 0
    iblankUndecided(:,:,:) = 0
    iblankTemp     (:,:,:) = 0

    ! zero out BI tables
    xBItable(:,:,:) = 0.0_CGREAL
    yBItable(:,:,:) = 0.0_CGREAL
    zBItable(:,:,:) = 0.0_CGREAL
    cELtable(:,:,:) = 0
    is_BIset(:,:,:) = 0
    ! Set the number extension of the marker of the detecting region

! ============================================================================
!   Initialize iblank to undecided value over all bodies in domain
! ============================================================================

!    iblank = 0
! DEBUG
!   FMN: To test is all iblanks are being filled.
!     iblank = -20*UNDECIDED_VALUE
! END DEBUG

    SELECT CASE (nDim)
      CASE (DIM_2D)
        k=1
        DO j=0,ny
        DO i=0,nx
          iblankUndecided(i,j,k) = UNDECIDED_VALUE
        END DO ! i
        END DO ! j
      CASE (DIM_3D)
        DO k=zb1,zb2  !0,nz
        DO j=yb1,yb2
        DO i=0,nx
          iblankUndecided(i,j,k) = UNDECIDED_VALUE
        END DO ! i
        END DO ! j
        END DO ! k
    END SELECT ! nDim

! ============================================================================
!   Loop over all bodies in domain and body surface elements
!   Find the points outside the domain
!   by Fang-Bao Tian; later modified by Ye Chen
! ============================================================================

    Flag_outside_Marker=0

    ! Limit of projection for each subdomain
    if (jProc.gt.0)then
       j1 = yc_start-1
    else
       j1 = 1
    endif
    if (jProc.lt.nProcY-1)then
       j2 = yc_end+1
    else
       j2 = ny-1
    endif

    if(kProc .gt. 0) then  !changed by song
       k1 = zc_start-1
    else
       k1 = 1       ! first slab
    endif
    if(kProc .lt. nProcZ-1) then 
       k2 = zc_end + 1
    else
       k2 = nz-1    ! last slab
    endif

    DO iBody=1,nBody
       DO m=1,nPtsBodyMarker(iBody)
         xM = xBodyMarker(m,iBody)
         yM = yBodyMarker(m,iBody)
         zM = zBodyMarker(m,iBody)
 
         !IF (( xM <= xc(1) .or. xM >= xc(nx-1) ) .OR.  &
         !    ( yM <= yc(1) .or. yM >= yc(ny-1) ) .OR.  &
         !    ( zM <= zc(1) .or. zM >= zc(nz-1) )) THEN
         
         ! Only do this for subdomain
         IF (( xM <= xc(1) .or. xM >= xc(nx-1) ) .OR.  &
             ( yM <= yc(j1) .or. yM >= yc(j2) ) .OR.  &
             ( zM <= zc(k1) .or. zM >= zc(k2) )) THEN         
            Flag_outside_Marker(m,iBody)=1
         ENDIF
       Enddo
    Enddo !end do iBody


! ============================================================================
!   Loop over all bodies in domain and body surface elements
! ============================================================================

    !CALL system_clock(iclock1)
    DO iBody=nBody,1,-1
      ! skip the MEMBRANE-type bodies
      ! if(unstruc_surface_type(iBody) == MEMBRANE) cycle
      
      DO m=1,totNumTriElem(iBody)

! ******************************************************************************
!       Extract cell indices and vertices of triangular element
! ******************************************************************************

         mVert1    = triElemNeig(1,m,iBody)
         mVert2    = triElemNeig(2,m,iBody)
         mVert3    = triElemNeig(3,m,iBody)

         xElem     = triElemCentx(m,iBody)
         yElem     = triElemCenty(m,iBody)
         zElem     = triElemCentz(m,iBody)

         xVert(1)  = xBodyMarker(mVert1,iBody)
         yVert(1)  = yBodyMarker(mVert1,iBody)
         zVert(1)  = zBodyMarker(mVert1,iBody)

         xVert(2)  = xBodyMarker(mVert2,iBody)
         yVert(2)  = yBodyMarker(mVert2,iBody)
         zVert(2)  = zBodyMarker(mVert2,iBody)

         xVert(3)  = xBodyMarker(mVert3,iBody)
         yVert(3)  = yBodyMarker(mVert3,iBody)
         zVert(3)  = zBodyMarker(mVert3,iBody)

! ******************************************************************************
! Skip if the entire element is outside computational domain. Note if the element
! is between the boundary and the 1st cell center, then the 1st cell center will
! be included and checked.
! For the membrane, larger domain need to be searched for the iblankUndecided, 
! edit by song
! ******************************************************************************
         IF(unstruc_surface_type(iBody) /= MEMBRANE) THEN
             Nm_tkns = 0
         ELSE
             !Use 2 instead of 1 to enlarge the Bounding Box for correct iblank
             !identification!
             Nm_tkns = 3!2!1
         ENDIF

         !For membrane structures (e.g., a hummingbird wing), include the body 
         ! elements outside the domain or subdomain.
         !  IF ( ( xVert(1) <= x(1)  .AND. xVert(2) <= x(1)          &
         !                           .AND. xVert(3) <= x(1)  ) .OR.  &
         !       ( xVert(1) >= x(nx) .AND. xVert(2) >= x(nx)         &
         !                           .AND. xVert(3) >= x(nx) ) .OR.  &
         !       ( yVert(1) <= y(1)  .AND. yVert(2) <= y(1)          &
         !                           .AND. yVert(3) <= y(1)  ) .OR.  &
         !       ( yVert(1) >= y(ny) .AND. yVert(2) >= y(ny)         &
         !                           .AND. yVert(3) >= y(ny) ) .OR.  &
         !       ( zVert(1) <= z(z_start)  .AND. zVert(2) <= z(z_start)       &
         !                                 .AND. zVert(3) <= z(z_start)) .OR. &
         !       ( zVert(1) >= z(z_end)    .AND. zVert(2) >= z(z_end)         &
         !                                 .AND. zVert(3) >= z(z_end) )     ) &
         !        CYCLE
         !  ENDIF            

        !  IF(unstruc_surface_type(iBody) == MEMBRANE) THEN

            !xIF ( ( xVert(1) <= x(1)  .AND. xVert(2) <= x(1)          &
            !x                        .AND. xVert(3) <= x(1)  ) .OR.   &
            !x    ( xVert(1) >= x(nx) .AND. xVert(2) >= x(nx)          &
            !x                        .AND. xVert(3) >= x(nx) ) .OR.   &
            !x    ( yVert(1) <= y(y_start)  .AND. yVert(2) <= y(y_start)     &
            !x                      .AND. yVert(3) <= y(y_start)  ) .OR.   &
            !x  ( yVert(1) >= y(y_end) .AND. yVert(2) >= y(y_end)          &
            !x                      .AND. yVert(3) >= y(y_end) ) .OR.   &
            !x  ( zVert(1) + membrane_tkns*Nm_tkns <= z(z_start)  .AND.       &
            !x    zVert(2) + membrane_tkns*Nm_tkns <= z(z_start)  .AND.       &
            !x    zVert(3) + membrane_tkns*Nm_tkns <= z(z_start)) .OR.        &
            !x  ( zVert(1) - membrane_tkns*Nm_tkns >= z(z_end)    .AND.       &
            !x    zVert(2) - membrane_tkns*Nm_tkns >= z(z_end)     .AND.      &
            !x    zVert(3) - membrane_tkns*Nm_tkns >= z(z_end) ))             &
            !x  CYCLE
         !  ENDIF

            IF ( ( xVert(1) <= x(1)  .AND. xVert(2) <= x(1)          &
                                    .AND. xVert(3) <= x(1)  ) .OR.   &
                ( xVert(1) >= x(nx) .AND. xVert(2) >= x(nx)          &
                                    .AND. xVert(3) >= x(nx) ) .OR.   &
                ( yVert(1) + membrane_tkns*Nm_tkns <= y(y_start)  .AND.       &
                  yVert(2) + membrane_tkns*Nm_tkns <= y(y_start)  .AND.       &
                  yVert(3) + membrane_tkns*Nm_tkns <= y(y_start)) .OR.        &
                ( yVert(1) - membrane_tkns*Nm_tkns >= y(y_end)    .AND.       &
                  yVert(2) - membrane_tkns*Nm_tkns >= y(y_end)    .AND.       &
                  yVert(3) - membrane_tkns*Nm_tkns >= y(y_end) )  .OR.   &
                ( zVert(1) + membrane_tkns*Nm_tkns <= z(z_start)  .AND.       &
                  zVert(2) + membrane_tkns*Nm_tkns <= z(z_start)  .AND.       &
                  zVert(3) + membrane_tkns*Nm_tkns <= z(z_start)) .OR.        &
                ( zVert(1) - membrane_tkns*Nm_tkns >= z(z_end)    .AND.       &
                  zVert(2) - membrane_tkns*Nm_tkns >= z(z_end)     .AND.      &
                  zVert(3) - membrane_tkns*Nm_tkns >= z(z_end) ))             &
                CYCLE


! ******************************************************************************
!        Find all cells within the bounding box of the vertices  
!        for each element
! ******************************************************************************

!DEBUG
!   IF (m == 1993) THEN
!      write(455,*)  iBody,m
!      write(455,*)  xVert(1), xVert(2) , xVert(3)
!      write(455,*)  yVert(1), yVert(2) , yVert(3)
!      write(455,*)  zVert(1), zVert(2) , zVert(3)
!   ENDIF
! END DEBUG

         xBoundMin = MINVAL(xVert(1:3))
         xBoundMax = MAXVAL(xVert(1:3))

         yBoundMin = MINVAL(yVert(1:3))
         yBoundMax = MAXVAL(yVert(1:3))

         zBoundMin = MINVAL(zVert(1:3))
         zBoundMax = MAXVAL(zVert(1:3))


         ! extend the search box for membrane-type structure, edited by song
         if(unstruc_surface_type(iBody) == MEMBRANE) then
            xBoundMin = xBoundMin - membrane_tkns*Nm_tkns !  /2.0_CGREAL
            xBoundMax = xBoundMax + membrane_tkns*Nm_tkns !  /2.0_CGREAL
            yBoundMin = yBoundMin - membrane_tkns*Nm_tkns !  /2.0_CGREAL
            yBoundMax = yBoundMax + membrane_tkns*Nm_tkns !  /2.0_CGREAL
            zBoundMin = zBoundMin - membrane_tkns*Nm_tkns !  /2.0_CGREAL
            zBoundMax = zBoundMax + membrane_tkns*Nm_tkns !   /2.0_CGREAL
         endif

!DEBUG
 !        WRITE(*,*)'iBody,m = ',iBody,m
 !        WRITE(*,*)'xBoundMin-Max =',xBoundMin,xBoundMax
 !        WRITE(*,*)'yBoundMin-Max =',yBoundMin,yBoundMax
 !        WRITE(*,*)'zBoundMin-Max =',zBoundMin,zBoundMax
 !        WRITE(*,*)'xMin-Max     =',MINVAL(x(1:nx)),MAXVAL(x(1:nx))
 !        WRITE(*,*)'yMin-Max     =',MINVAL(y(1:ny)),MAXVAL(y(1:ny))
 !        WRITE(*,*)'zMin-Max     =',MINVAL(z(1:nz)),MAXVAL(z(1:nz))
!END DEBUG

         iCellMin = 0
         iCellMax = 0
         jCellMin = yc_start-1
         jCellMax = yc_start-1
         kCellMin = zc_start-1  !0   song
         kCellMax = zc_start-1  !0

! ******************************************************************************
!        i-direction 
! ******************************************************************************

         iMin = 0
         iMax = nx
         DO i = iMin, iMax-1
           IF ( (xc(i)-xBoundMin)*(xc(i+1)-xBoundMin) <= 0.0_CGREAL ) &
              iCellMin = i
           IF ( (xc(i)-xBoundMax)*(xc(i+1)-xBoundMax) <= 0.0_CGREAL ) &
              iCellMax = i+1
         ENDDO
         IF (iCellMin == 0 )                    iCellMin = 1
         IF (iCellMax == 0 .OR. iCellMax == nx) iCellMax = nx-1

! ******************************************************************************
!        j-direction 
! ******************************************************************************

         jMin = yc_start-1
         jMax = yc_end+1
         DO j = jMin, jMax-1
           IF ( (yc(j)-yBoundMin)*(yc(j+1)-yBoundMin) <= 0.0_CGREAL ) &
              jCellMin = j
           IF ( (yc(j)-yBoundMax)*(yc(j+1)-yBoundMax) <= 0.0_CGREAL ) &
              jCellMax = j+1
         ENDDO
         IF (jCellMin == yc_start-1 )                          jCellMin = yc_start
         IF (jCellMax == yc_start-1 .OR. jCellMax == yc_end+1) jCellMax = yc_end

! ******************************************************************************
!        k-direction 
! ******************************************************************************

         SELECT CASE (nDim)
         CASE (DIM_2D)
            kCellMin = 1
            kCellMax = nz-1

         CASE (DIM_3D)

            !kMin = 0
            !kMax = nz
            !DO k = kMin, kMax-1
            !  IF ( (zc(k)-zBoundMin)*(zc(k+1)-zBoundMin) < 0.0_CGREAL ) &
            !  kCellMin = k
            !  IF ( (zc(k)-zBoundMax)*(zc(k+1)-zBoundMax) < 0.0_CGREAL ) &
            !  kCellMax = k+1
            !ENDDO
            !IF (kCellMin == 0 ) kCellMin = 1
            !IF (kCellMax == 0 .OR. kCellMax == nz) kCellMax = nz-1
            
            kMin = zc_start-1
            kMax = zc_end+1
            DO k = kMin, kMax-1
               IF ( (zc(k)-zBoundMin)*(zc(k+1)-zBoundMin) <= 0.0_CGREAL ) &
                    kCellMin = k
               IF ( (zc(k)-zBoundMax)*(zc(k+1)-zBoundMax) <= 0.0_CGREAL ) &
                    kCellMax = k+1
            ENDDO
            IF (kCellMin == zc_start-1 )                          kCellMin = zc_start
            IF (kCellMax == zc_start-1 .OR. kCellMax == zc_end+1) kCellMax = zc_end

         END SELECT ! nDim

! DEBUG
!if(kCellMin < 2) then
!      WRITE(*,*) ' m   = ',m
!      WRITE(*,*) 'iCellMin-Max = ',iCellMin,iCellMax
!      WRITE(*,*) 'jCellMin-Max = ',jCellMin,jCellMax
!      WRITE(*,*) 'kCellMin-Max = ',kCellMin,kCellMax
!      write(*,'(3F12.5)')xVert(1:3)
!      write(*,'(3F12.5)')yVert(1:3)
!      write(*,'(3F12.5)')zVert(1:3)
!endif
!END DEBUG

! ******************************************************************************
!        set iblankUndecided to NEGATIVE undecided value for the ring 
!            of cells extracted. Also save bodyNum
! ******************************************************************************

         !if(lProc==1) print*,kCellMin, kCellMax
         DO k = kCellMin, kCellMax
         DO j = jCellMin, jCellMax
         DO i = iCellMin, iCellMax
              iblankUndecided(i,j,k)  = -UNDECIDED_VALUE
         END DO ! i 
         END DO ! j 
         END DO ! k 

      END DO ! m
    END DO ! iBody
    !CALL system_clock(iclock3)

! ============================================================================
!   Loop over cells whose iblank is -UNDECIDED_VALUE for unstructured surfaces
!    invoking dot normal algorithm
! ============================================================================

    !CALL system_clock(iclock4)
    SELECT CASE(nDim)
    CASE (DIM_2D)
       kMin=1
       kMax=1

    CASE (DIM_3D)
       jMin = yc_start
       jMax = yc_end
       kMin = zc_start
       kMax = zc_end
    END SELECT ! nDim

    DO k=kMin,kMax
    DO j=jMin,jMax
    DO i=1,nx-1

      tmpiblank    = 0
      xP = xc(i)
      yP = yc(j)
      zP = zc(k)

      IF ( iblankUndecided(i,j,k) == -UNDECIDED_VALUE ) THEN
! ******************************************************************************
!       Search for vertex closest to cells in band
!       Extract the elements that share that vertex
!       Drop the normal and find BI point
!       Check if normal intercept lies inside closest triangular element
!       Authors: Fady and Rajat Oct 18,2005; modified by Haoxiang later
! ******************************************************************************
        DO iBody = nBody,1,-1

          ! set a finite thickness for the MEMBRANE-type bodies
          IF(unstruc_surface_type(iBody) == MEMBRANE) THEN

            CALL calc_bodyIntercept_Unstruc(iBody,i,j,k,xBM,yBM,zBM,cElement,distElement,dotNorm)

            distElement = SQRT( (xP-xBM)**2 + (yP-yBM)**2 + (zP-zBM)**2)
            !write(*,'(I5,8F12.5)')n,xP,yP,zP,xBM,yBM,zBM,distElement
            !if(i==73.and.j==31.and.k==40) then
            !  write(*,'(I5,8F12.5)')iBody,xP,yP,zP,xBM,yBM,zBM,distElement,membrane_tkns
            !endif

            IF (distElement < membrane_tkns/2.0_CGREAL) THEN 
               iblankTemp(i,j,k)  = 1
               bodyNum(i,j,k)     = iBody
               !write(*,'(I5,8F12.5)')n,xP,yP,zP,xBM,yBM,zBM,distElement
               GOTO 888
            ENDIF
          ELSE  !if(unstruc_ ...
            ! for other solid bodies, check the dot product with the surface
            ! normal.
            CALL calc_bodyIntercept_Unstruc(iBody,i,j,k,xBM,yBM,zBM,cElement,distElement,dotNorm)

            IF (dotNorm >= 0.0_CGREAL) THEN
              iblankTemp(i,j,k)  = 1
              bodyNum(i,j,k)     = iBody
              GOTO 888
            ENDIF
          ENDIF  !if(unstruc_ ...
        ENDDO

        iblankTemp(i,j,k)  = 0

888     CONTINUE

      END IF ! iblankUndecided

    ENDDO ! i
    ENDDO ! j
    ENDDO ! k
    !CALL system_clock(iclock5)

!DEBUG
    !IF (ntime.eq.836) THEN
    !write(6000+lProc,*)'VARIABLES="X","Y","Z","xb","yb","zb","IBUNDEC","IBLANK","BdNum","is"'
    !write(6000+lProc,*)'ZONE F=POINT, I=',nx-1,', J=',jSlices+2,' ,K=',kSlices+2   !nz-1
    !do k=zc_start-1,zc_end+1   !1,nz-1
    !do j=yc_start-1,yc_end+1
    !do i=1,nx-1
    !   write(6000+lProc,'(6(2X,F12.5),4(3X,I10))')xc(i),yc(j),zc(k), &
    !               xBItable(i,j,k),yBItable(i,j,k),zBItable(i,j,k), &
    !               iblankUndecided(i,j,k), &
    !               iblankTemp(i,j,k),bodyNum(i,j,k),is_BIset(i,j,k)
    !enddo
    !enddo
    !enddo
    !close(6000+lProc)
    !call MPI_BARRIER(flow_comm,ierr)
    !stop
    !ENDIF
!DEBUG

! ============================================================================
!   Set undecided iblank values outside the ring of cells for body
!    by searching horizontal, similar to a ray tracing routine. 
!    Set iblank value at grid cell by searching for first DECIDED VALUE
!    Move along i-direction, j-direction, then k-direction
! ============================================================================

! ******************************************************************************
!   i-direction 
! ******************************************************************************

    !CALL system_clock(iclock6)
    iBeg = 0
    DO k=kMin,kMax
    DO j=jMin,jMax
      initLoopI: DO ii=1,nx-1  !nx
         !if iblank is not set, skip the cell.
         IF ( iblankUndecided(ii,j,k) == UNDECIDED_VALUE ) CYCLE
         
         !Otherwise, set the 1st cell of the line. Note that if 
         !the entire line is skipped, then the 1st cell by default is a fluid cell
         iblankTemp(iBeg,j,k) = iblankTemp(ii,j,k)
            bodyNum(iBeg,j,k) =    bodyNum(ii,j,k)   
         EXIT initLoopI  ! terminate the loop
      END DO  initLoopI
    END DO ! j
    END DO ! k

    DO k=kMin,kMax
    DO j=jMin,jMax
    DO i=1,nx-1   !nx
       ! if iblank is set, skip the cell
       IF ( iblankUndecided(i,j,k) /= UNDECIDED_VALUE ) CYCLE

       ! otherwise, set the iblank to be the same as the preceding cell.
       iblankTemp(i,j,k)  = iblankTemp(i-1,j,k)
       IF (iblankTemp(i,j,k)==1) bodyNum(i,j,k)     = bodyNum(i-1,j,k)
    END DO ! i
    END DO ! j
    END DO ! k

!DEBUG
    !IF (ntime.eq.806) THEN
    !write(7000+lProc,*)'VARIABLES="X","Y","Z","IBUNDEC","IBLANK", "BodyNum"'
    !write(7000+lProc,*)'ZONE F=POINT, I=',nx-1,', J=',jSlices+2,' ,K=',kSlices+2   !nz-1
    !do k=zc_start-1,zc_end+1   !1,nz-1
    !do j=yc_start-1,yc_end+1
    !do i=1,nx-1
!!$      !write(700+lProc,'(3(3X,  I10  ),3(3X,I10))')i,j,k, iblankUndecided(i,j,k),iblankTemp(i,j,k),bodyNum(i,j,k)
    !   write(7000+lProc,'(3(2X,F12.5),3(3X,I10))')xc(i),yc(j),zc(k),iblankUndecided(i,j,k),iblankTemp(i,j,k),bodyNum(i,j,k)
    !enddo
    !enddo
    !enddo
    !close(7000+lProc)
    !call MPI_BARRIER(flow_comm,ierr)
    !stop
    !ENDIF
!DEBUG


! ******************************************************************************
!   j-direction 
! ******************************************************************************

    !jBeg = 0
    !DO k=kMin,kMax
    !DO i=1,nx-1
    !  initLoopJ: DO jj=yb1,yb2  !ny
    !     !if iblank is not set, skip the cell.
    !     IF ( iblankUndecided(i,jj,k) == UNDECIDED_VALUE ) CYCLE

    !     !Otherwise, set the 1st cell of the line. Note that if 
    !     !the entire line is skipped, then the 1st cell by default is a fluid cell
    !     iblankTemp(i,jBeg,k) = iblankTemp(i,jj,k)
    !     bodyNum(i,jBeg,k)    = bodyNum(i,jj,k) 
    !     EXIT initLoopJ  ! terminate the loop
    !  END DO  initLoopJ
    !END DO ! i
    !END DO ! k

    !DO k=kMin,kMax
    !DO j=yb1,yb2
    !DO i=1,nx-1
    !   ! if iblank is set, skip the cell
    !   IF ( iblankUndecided(i,j,k) /= UNDECIDED_VALUE ) CYCLE

    !   ! otherwise, set the iblank to be the same as the preceding cell.
    !   iblankTemp(i,j,k) = iblankTemp(i,j-1,k)
    !   IF (iblankTemp(i,j,k)==1) bodyNum(i,j,k)    = bodyNum(i,j-1,k)
    !END DO ! i
    !END DO ! k
    !END DO ! j

    !IF( (ntime.eq.6012).AND.(lProc.eq.1) )THEN
    !  WRITE(*,*)'jBODY#',bodyNum(65,66,12)
    !ENDIF

!DEBUG
    !write(800+lProc,*)'VARIABLES="X","Y","Z","IBUNDEC","IBLANK", "BodyNum"'
    !write(800+lProc,*)'ZONE F=POINT, I=',nx-1,', J=',jSlices+2,' ,K=',kSlices+2   !nz-1
    !do k=zc_start-1,zc_end+1   !1,nz-1
    !do j=yc_start-1,yc_end+1
    !do i=1,nx-1
!!$       !write(800+lProc,'(3(3X,  I10  ),3(3X,I10))')i,j,k, iblankUndecided(i,j,k),iblankTemp(i,j,k),bodyNum(i,j,k)
    !   write(800+lProc,'(3(2X,F12.5),3(3X,I10))')xc(i),yc(j),zc(k),iblankUndecided(i,j,k),iblankTemp(i,j,k),bodyNum(i,j,k)
    !enddo
    !enddo
    !enddo
    !close(800+lProc)
!DEBUG

! ******************************************************************************
!   k-direction 
!
! Currently we disable the z-direction tracing, as the entire z-line within the subdomain 
! could be inside a solid body and iblank could be mis-determined if we assume that for such a 'clean' line,
! the 1st cell is always a fluid cell. Tracing in the x- and y- directions is sufficient. 
!                                                       -- Haoxiang Luo, April 2015
! ******************************************************************************

!    IF (nDim == DIM_3D) THEN
!      kBeg = zc_start-1  !0
!      DO j=1,ny-1
!      DO i=1,nx-1
!        initLoopK: DO kk=kMin,kMax  !kMax+1
!           !if iblank is not set, skip the cell.
!           IF ( iblankUndecided(i,j,kk) == UNDECIDED_VALUE ) CYCLE
!
!           !Otherwise, set the 1st cell of the line.
!           iblankTemp(i,j,kBeg)  = iblankTemp(i,j,kk)
!           bodyNum(i,j,kBeg)     = bodyNum(i,j,kk)
!           EXIT initLoopK   ! terminate the loop
!        END DO  initLoopK
!      END DO ! i
!      END DO ! j
!
!      DO k=kMin,kMax+1
!      DO j=1,ny-1
!      DO i=1,nx-1
!         ! if iblank is set, skip the cell
!         IF ( iblankUndecided(i,j,k) /= UNDECIDED_VALUE ) CYCLE
!
!         ! otherwise, set the iblank to be the same as the preceding cell.
!         iblankTemp(i,j,k)  = iblankTemp(i,j,k-1)
!         IF (iblankTemp(i,j,k)==1) bodyNum(i,j,k)     = bodyNum(i,j,k-1)         
!      END DO ! i
!      END DO ! j
!      END DO   ! k
!    ENDIF !nDim

    !Ye,ad-hoc for realistic epitrochoid aorta model=====
    !do k=kMin,kMax
    !do j=1,ny-1
    !do i=1,nx-1
    !   if ((xc(i)>=-0.2).and.(xc(i)<=0.02))then
    !      rtemp=0.428+(0.02-xc(i))*(0.478-0.428)/0.22
    !   elseif ((xc(i)>0.02).and.(xc(i)<=0.6))then
    !      rtemp=0.428+(xc(i)-0.02)*(0.478-0.428)/0.58
    !   else
    !      rtemp=0.478
    !   endif

    !   if(yc(j)**2+zc(k)**2-rtemp**2>0)then
    !      iblankTemp(i,j,k) = 1
    !      bodyNum(i,j,k) = 1
    !   endif
    !enddo
    !enddo
    !enddo
    !=====================================================

     !Ye,ad-hoc for straight cylinder aorta model==========
!     do k=kMin,kMax
!        do j=1,ny-1
!         do i=1,nx-1
!            if(yc(j)**2+zc(k)**2-0.428**2>0)then
!               iblankTemp(i,j,k) = 1
!               bodyNum(i,j,k) = 1
!            endif
!         enddo
!        enddo
!     enddo

!      !Ye,ad-hoc for epitrochoid-type aorta model=====
!      DO j=jMin,jMax
!         DO k=kMin,kMax
!            IF (yc(j)**2+zc(k)**2-1.08712**2>=0) THEN!0.428*2.54
!               iblankTemp(0,j,k) = 1
!               bodyNum(0,j,k) = 1
!               iblankTemp(nx,j,k) = 1
!               bodyNum(nx,j,k) = 1
!            ENDIF
!         ENDDO
!      ENDDO
     
!      DO j=jMin,jMax
!         DO k=kMin,kMax
!            IF (yc(j)**2+zc(k)**2-1.08712**2<0) CYCLE
!            initLoopI1: DO ii=1,nx-1
!               IF( (iblankTemp(ii,j,k).eq.0).OR. &
!                   (iblankUndecided(ii,j,k).EQ.UNDECIDED_VALUE) ) THEN
!                   iblankTemp(ii,j,k)  = iblankTemp(ii-1,j,k)
!                   IF (iblankTemp(ii,j,k)==1) bodyNum(ii,j,k) = bodyNum(ii-1,j,k)
!                   GOTO 444
!               ENDIF
!               EXIT initLoopI1
! 444        CONTINUE
!            ENDDO initLoopI1
!         ENDDO!i
!      ENDDO!j

!      DO j=jMin,jMax
!         DO k=kMin,kMax
!            IF (yc(j)**2+zc(k)**2-1.08712**2<0) CYCLE
!            initLoopI2: DO ii=nx-1,1,-1
!               IF( (iblankTemp(ii,j,k).eq.0).OR. &
!                   (iblankUndecided(ii,j,k).EQ.UNDECIDED_VALUE) ) THEN
!                   iblankTemp(ii,j,k)  = iblankTemp(ii+1,j,k)
!                   IF (iblankTemp(ii,j,k)==1) bodyNum(ii,j,k) = bodyNum(ii+1,j,k)
!                   GOTO 555
!               ENDIF
!               EXIT initLoopI2
! 555        CONTINUE
!            ENDDO initLoopI2
!         ENDDO!i
!      ENDDO!j

     !Ye, fix the iblankTemp at the outflow, xc=nx
     DO j=jMin,jMax
        DO k=kMin,kMax
           iblankTemp(nx,j,k)  = iblankTemp(nx-1,j,k)
           bodyNum(nx,j,k) = bodyNum(nx-1,j,k)
        ENDDO
     ENDDO
     !=====================================================

    !DEBUG
    !IF (ntime.eq.0) THEN
    !write(900+lProc,*)'VARIABLES="X","Y","Z","IBUNDEC","IBLANK", "BodyNum"'
    !write(900+lProc,*)'ZONE F=POINT, I=',nx+1,', J=',jSlices+2,' ,K=',kSlices+2   !nz-1
    !do k=zc_start-1,zc_end+1   !1,nz-1
    !do j=yc_start-1,yc_end+1
    !do i=0,nx
       !write(900+lProc,'(3(3X,  I10  ),3(3X,I10))')i,j,k, iblankUndecided(i,j,k),iblankTemp(i,j,k),bodyNum(i,j,k)
    !   write(900+lProc,'(3(2X,F12.5),3(3X,I10))')xc(i),yc(j),zc(k),iblankUndecided(i,j,k),iblankTemp(i,j,k),bodyNum(i,j,k)
    !enddo
    !enddo
    !enddo
    !close(900+lProc)
    !ENDIF

! ============================================================================
!   Extend iblank for 2D simulations 
! ============================================================================

    IF (nDim == DIM_2D) THEN
      DO k=kMin+1,nz
      DO j=1,ny-1
      DO i=1,nx-1
        iblankTemp(i,j,k)  = iblankTemp(i,j,1)
        IF ( iblankTemp(i,j,k) == 1 ) bodyNum(i,j,k) = bodyNum(i,j,1)
      END DO !i
      END DO !j
      END DO !k
    END IF ! nDim

    DO k=zc_start,zc_end  !1,nz-1
    !DO j=1,ny-1
    !DO i=1,nx-1
    !Ye,set iblank=1 at outer boundary
    DO j=yc_start,yc_end
    DO i=0,nx
       IF  ( iblankTemp(i,j,k) == 1 ) THEN
          IF ( iblank(i,j,k) == 0  .AND. ntime > ntime_start+1 ) THEN
             !WRITE(ifuFreshCellOut,*)ntime,i,j,k,'   --- dead cell'
             num_dead     = num_dead+1
          ENDIF
       ENDIF
       IF  ( iblankTemp(i,j,k) == 0 ) THEN
          IF ( iblank(i,j,k) == 1 .AND. ntime > 1 ) THEN
             !fresh_cell(i,j,k)=  1
             !bodyNum(i,j,k)   = bodyNum(i,j,k)  ! this is done to show that this array has correct value
             num_fresh        = num_fresh+1

             !WRITE(ifuFreshCellOut,*)ntime,i,j,k,'   --- fresh cell'
          ENDIF
       ENDIF
       iblank(i,j,k)  = iblankTemp(i,j,k)
    END DO ! i
    END DO ! j
    END DO ! k

    !Need to call MPI_ALLReduce to sum up num_dead and num_fresh
    
    nDead  = num_dead                                
    !if(boundary_formulation == SSM_METHOD ) nFresh = num_fresh  !

!    IF (MOD(ntime,nmonitor)==0) THEN
!      PRINT*,'Number of Dead  Cells = ',num_dead
!      PRINT*,'Number of Fresh Cells = ',num_fresh
!    ENDIF ! ntime

!    CALL system_clock(iclock2,clock_rate)
!    IF ( ntime==1           .OR. &
!         ntime==ntime_start .OR. &
!         MOD(ntime,nmonitor) == 0 )  THEN
!      WRITE(*,*) 'CPU Time for Fast Initial Iblank = ',&
!      REAL(iclock2-iclock1,KIND=CGREAL)/REAL(clock_rate,KIND=CGREAL)
!      WRITE(*,*) '    CPU Time for Initial Setup = ',&
!      REAL(iclock3-iclock1,KIND=CGREAL)/REAL(clock_rate,KIND=CGREAL)
!      WRITE(*,*) '    CPU Time for dotNorm = ',&
!      REAL(iclock5-iclock4,KIND=CGREAL)/REAL(clock_rate,KIND=CGREAL)
!      WRITE(*,*) '    CPU Time for Filling Iblank = ',&
!      REAL(iclock2-iclock6,KIND=CGREAL)/REAL(clock_rate,KIND=CGREAL)

       !WRITE(*,*) 'SUM of IBLANK = ', sum(IBLANK(1:nx-1,1:ny-1,1:nz-1))
!    END IF ! ntime

    ! count the number of iblank cells
    sumIblank = 0
    DO k = zc_start,zc_end    !1, nz-1
    DO j = yc_start,yc_end
    DO i = 1, nx-1
       IF ( iblank(i,j,k) == 1 )   THEN
          sumIblank = sumIblank + 1
       ENDIF
    ENDDO ! i 
    ENDDO ! j
    ENDDO ! k

    !IF ( MOD(ntime,nmonitor) == 0) &
    !    PRINT*,'lProc',lProc,'set_iblank_body_fast: sumIblank = ',sumIblank

    !Ye
    !IF (ntime.eq.836) THEN
    !write(2000+lProc,*)'VARIABLES="X","Y","Z","IBLANK", "BodyNum"'
    !write(2000+lProc,*)'ZONE F=POINT, I=',nx-1,', J=',jSlices+2,' ,K=',kSlices+2   !nz-1
    !do k=zc_start-1,zc_end+1   !1,nz-1
    !do j=yc_start-1,yc_end+1
    !do i=1,nx-1
    !   write(2000+lProc,'(3(2X,F12.5),2(3X,I10))')xc(i),yc(j),zc(k),iblank(i,j,k),bodyNum(i,j,k)
    !enddo
    !enddo
    !enddo
    !close(2000+lProc)
    !call MPI_BARRIER(flow_comm,ierr)
    !stop
    !ENDIF

    !Ye
    !IF ((ntime.eq.2502).AND.(lProc.eq.1)) THEN
    !   write(88,*)iblank(79,36,12)
    !   write(88,*)'BF',iblank(78,36,12),iblank(80,36,12)
    !   write(88,*)'LR',iblank(79,35,12),iblank(79,37,12)
    !   write(88,*)'DU',iblank(79,36,11),iblank(79,36,13)
    !ENDIF

   END SUBROUTINE set_iblank_body_fast

!------------------------------------------------------------------------------
   SUBROUTINE calc_bodyIntercept_Unstruc(iBody,iCell,jCell,kCell,xBI,yBI,zBI, &
                                         closestElement, distBIElem, dotNorm)
!------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays

    USE MPI_module
    USE fsi_module

    IMPLICIT NONE

!... parameters

    INTEGER, INTENT(IN)             :: iBody,iCell,jCell,kCell
    REAL(KIND=CGREAL), INTENT(OUT)  :: xBI,yBI,zBI
    INTEGER, INTENT(OUT)            :: closestElement
    REAL(KIND=CGREAL) , INTENT(OUT) :: distBIElem, dotNorm
!... loop variables

    INTEGER :: m,n,nc
    
!... local variables

    INTEGER,           DIMENSION(:),ALLOCATABLE   :: NeighElemInd
    REAL(KIND=CGREAL), DIMENSION(:),ALLOCATABLE   :: distMarker

    INTEGER                  :: cElement_in,cElement_out
    INTEGER                  :: iErr0,nBodySelect,numNeighElement
    INTEGER                  :: elemInd,node1,node2,node3,nMarker
    INTEGER                  :: nCheck
    INTEGER,DIMENSION(1)     :: iDummy(1)
    INTEGER                  :: shortestProbe
    INTEGER,DIMENSION(1:MSIZE) :: cElement,closestVert

    REAL(KIND=CGREAL)        :: inside_flag
    REAL(KIND=CGREAL)        :: dist_in, distTemp
    REAL(KIND=CGREAL)        :: xCell,yCell,zCell
    REAL(KIND=CGREAL)        :: dMin,dsIntercept,xM,yM,zM
    REAL(KIND=CGREAL)        :: areaDiffMin
    REAL(KIND=CGREAL)        :: distBIElemMin
    REAL(KIND=CGREAL)        :: planeConst,distanceToPlane,distPointToPlane
    REAL(KIND=CGREAL)        :: side12,side23,side31,side14,side24,side34
    REAL(KIND=CGREAL)        :: area123,area124,area234,area314
    REAL(KIND=CGREAL)        :: semiPerimeter123,semiPerimeter124,semiPerimeter234,semiPerimeter314
    REAL(KIND=CGREAL)        :: epsiArea,areaDiff
    REAL(KIND=CGREAL)        :: xBI_in,yBI_in,zBI_in
    REAL(KIND=CGREAL)        :: xBI_out,yBI_out,zBI_out
    REAL(KIND=CGREAL)        :: xBITemp, yBITemp, zBITemp
    REAL(KIND=CGREAL),DIMENSION(3) :: xVert, yVert, zVert
    REAL(KIND=CGREAL),DIMENSION(3) :: xVertN, yVertN, zVertN
    REAL(KIND=CGREAL),DIMENSION(1:MSIZE) :: xBIT,yBIT,zBIT,dist
    REAL(KIND=CGREAL),DIMENSION(1:MSIZE) :: xBIN,yBIN,zBIN
    REAL(KIND=CGREAL)        :: xNorm,yNorm,zNorm,xCG,yCG,zCG
    REAL(KIND=CGREAL)        :: xBInorm, yBInorm, zBInorm
    REAL(KIND=CGREAL)        :: xBInormTemp, yBInormTemp, zBInormTemp
!******************************************************************************

    xCell = xc(iCell)
    yCell = yc(jCell)
    zCell = zc(kCell)

    dotNorm = -1.0_CGREAL

    nMarker = nPtsBodyMarker(iBody)

! NCheck:  Number of closesest nodes to check
! and also CPU time for finding body intercept.

    nCheck = 1

    IF (nCheck > MSIZE) THEN
       PRINT*,'nCheck in GCM_calc_bodyIntercept_Unstruc is limited to', MSIZE
       PRINT*,'Increase array size'
       STOP
    ENDIF

! ============================================================================
!   Allocate local array 
! ============================================================================

    ALLOCATE(distMarker(nMarker),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'calc_bodyIntercept_Unstruc: Memory Allocation Error for distMarker'
      STOP
    ENDIF ! ierr

    ALLOCATE(NeighElemInd(NSIZE),STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'calc_bodyIntercept_Unstruc: Memory Allocation Error for NeighElemInd'
      STOP
    ENDIF ! ierr

! ============================================================================
!   Get closestMarker for cell 
! ============================================================================

    dMin = 1.0E+16_CGREAL

    DO m = 1, nMarker
      xM = xBodyMarker(m,iBody)
      yM = yBodyMarker(m,iBody)
      zM = zBodyMarker(m,iBody)
      distMarker(m) = (xM-xCell)**2 + (yM-yCell)**2 + (zM-zCell)**2

      !if the marker is outside the subdomain, set a large distance
      !if(Flag_outside_Marker(m,iBody)==1) distMarker(m)=1.0E3_CGREAL
    ENDDO

    DO nc = 1,NCheck
       iDummy                        = MINLOC(distMarker(1:nMarker))
       closestVert(nc)               = iDummy(1)
       distMarker(closestVert(nc))   = 1.0E20_CGREAL
    ENDDO

    !Ye,debug
    !if ( (iCell.eq.1).and.(jCell.eq.56).and.(kCell.eq.41) )then
    !   write(301,*)ntime,nCheck,iBody,nMarker
    !   write(301,*)closestVert(1:NCheck)
    !   write(301,*)distMarker(closestVert(1:NCheck))
    !   write(301,*)'---'
    !endif

! ============================================================================
!   Find elements that share closest node/marker
! ============================================================================
    DO nc = 1, NCheck

    numNeighElement = 0
    DO m=1,totNumTriElem(iBody)
       IF ( triElemNeig(1,m,iBody) == closestVert(nc) .OR. &
            triElemNeig(2,m,iBody) == closestVert(nc) .OR. &
            triElemNeig(3,m,iBody) == closestVert(nc)      ) THEN
          numNeighElement               = numNeighElement + 1
          NeighElemInd(numNeighElement) = m
       ENDIF
    ENDDO ! m

! ============================================================================
!   Trap error if array NeighElemenInd overflows 
! ============================================================================

    IF ( numNeighElement > NSIZE ) THEN
      WRITE(STDOUT,*) &
       'search_vertex_dotNorm: Memory Overflow Error for NeighElemInd'
      WRITE(STDOUT,*) ' Allocated size = ',NSIZE
      WRITE(STDOUT,*) ' Current size   = ',numNeighElement
      WRITE(STDOUT,*) ' Aborting Run'
      STOP
    ENDIF ! NeighElemInd

! ============================================================================
!   Check which element contains normal intercept
! ============================================================================

    distBIElemMin = 1.0E+16_CGREAL
    areaDiffMin   = 1.0E+16_CGREAL
    epsiArea      = 1.0E-6_CGREAL

    closestElement= 0
!===========================
    DO n = 1,numNeighElement
     elemInd = NeighElemInd(n)

     node1   = triElemNeig(1,elemInd,iBody)
     node2   = triElemNeig(2,elemInd,iBody)
     node3   = triElemNeig(3,elemInd,iBody)

     xVert(1)  = xBodyMarker(node1,iBody)
     yVert(1)  = yBodyMarker(node1,iBody)
     zVert(1)  = zBodyMarker(node1,iBody)

     xVert(2)  = xBodyMarker(node2,iBody)
     yVert(2)  = yBodyMarker(node2,iBody)
     zVert(2)  = zBodyMarker(node2,iBody)

     xVert(3)  = xBodyMarker(node3,iBody)
     yVert(3)  = yBodyMarker(node3,iBody)
     zVert(3)  = zBodyMarker(node3,iBody)

     !extract the emement normal
     xNorm     = triElemNormx(elemInd,iBody)
     yNorm     = triElemNormy(elemInd,iBody)
     zNorm     = triElemNormz(elemInd,iBody)

     !extract the bodymarker normal
     xVertN(1)  = xNormBodyMarker(node1,iBody)
     yVertN(1)  = yNormBodyMarker(node1,iBody)
     zVertN(1)  = zNormBodyMarker(node1,iBody)

     xVertN(2)  = xNormBodyMarker(node2,iBody)
     yVertN(2)  = yNormBodyMarker(node2,iBody)
     zVertN(2)  = zNormBodyMarker(node2,iBody)

     xVertN(3)  = xNormBodyMarker(node3,iBody)
     yVertN(3)  = yNormBodyMarker(node3,iBody)
     zVertN(3)  = zNormBodyMarker(node3,iBody)

! ******************************************************************************
! equation of plane (note our normals are unit normals)
!     
!  n  x + n  y + n  z + planeConst = 0
!   x         y      z
! ******************************************************************************

     planeConst=- xNorm*xVert(1) &
                - yNorm*yVert(1) &
                - zNorm*zVert(1)

! ******************************************************************************
! Compute coordinates of normal intercept
!      
! Consider point Po = (xo,yo,zo)  not on plane
!   and  point   P1 = (x1,y1,z1)  on the plane
!                                                               ^
! equation of line through Po normal to plane is  P(s) = Po + s n
!      
! normal distance from Po to Plane is given by
!            ^                               ^         ^ ^
!            n. ( P(s) - P1 ) = 0   => so = -n.(Po-P1)/n.n                     ^ ^
!                                         = -(n xo + n yo + n zo + planeConst)/n.n
!                                              x      y      z
!
!                                                 ^
!   subsequently normal intersection point = Po + so n  
!                   
! ******************************************************************************

     distanceToPlane = -(  xNorm*xCell  &
                         + yNorm*yCell  &
                         + zNorm*zCell  &
                         + planeConst )

     xBITemp = xCell + xNorm*distanceToPlane
     yBITemp = yCell + yNorm*distanceToPlane
     zBITemp = zCell + zNorm*distanceToPlane
                
! ******************************************************************************
!    Check if BI inside the triangle of the surface element
!     through area differences
! ******************************************************************************
    
     CALL check_BIInsideTriangle(xVert,yVert,zVert,xNorm,yNorm,zNorm,       &
                                 xBITemp,yBITemp,zBITemp,area123,areaDiff)

!!!!!!!!!!!!!!!!!!!!
!    IF (iCell == 1 .and. jCell == 1 .and. kCell == 1) then
!      WRITE(355,*)'ZONE'
!      WRITE(355,*)xBodyMarker(node1,iBody),yBodyMarker(node1,iBody),zBodyMarker(node1,iBody)
!      WRITE(355,*)xBodyMarker(node2,iBody),yBodyMarker(node2,iBody),zBodyMarker(node2,iBody)
!      WRITE(355,*)xBodyMarker(node3,iBody),yBodyMarker(node3,iBody),zBodyMarker(node3,iBody)
!      WRITE(355,*)xBodyMarker(node1,iBody),yBodyMarker(node1,iBody),zBodyMarker(node1,iBody)
!      WRITE(355,*)'ZONE'
!      WRITE(355,*)xCell,yCell,zCell
!      WRITE(355,*)xBItemp,yBItemp,zBItemp
!      WRITE(356,*) 'n,elemInd,areaDiff,area123 = ',n,elemInd,areadiff,area123
!   ENDIF
!!!!!!!!!!!!!!!!!!!!!

     distBIElem  =  SQRT( (xCell - xBITemp)**2 &
                        + (yCell - yBITemp)**2 &
                        + (zCell - zBITemp)**2 )

     IF ( ABS(areaDiff) < epsiArea*area123) THEN

         xBInormTemp = xNorm
         yBInormTemp = yNorm
         zBInormTemp = zNorm

     ELSE

        !The projection is outside the triangle; find the closest point on the
        !triangle
        CALL calc_BIOutsideTriangle(xVert,yVert,zVert,xCell,yCell,zCell,     &
                                  xVertN,yVertN,zVertN,                      &
                                  xBITemp,yBITemp,zBITemp,distBIElem,        &
                                  xBInormTemp,yBInormTemp,zBInormTemp)

     ENDIF ! areaDiff

     IF (distBIElem <= distBIElemMin) THEN

          distBIElemMin = distBIElem
          closestElement = elemInd
          xBI = xBITemp
          yBI = yBITemp
          zBI = zBITemp

          xBInorm = xBInormTemp
          yBInorm = yBInormTemp
          zBInorm = zBInormTemp

     ENDIF ! distBIElem

    ENDDO ! DO n = 1,numNeighElement
!==========================================================

    xBIT(nc) = xBI
    yBIT(nc) = yBI
    zBIT(nc) = zBI
    cElement(nc) = closestElement

    xBIN(nc) = xBInorm
    yBIN(nc) = yBInorm
    zBIN(nc) = zBInorm

  ENDDO  ! nc; end of nCheck

  DO nc = 1,nCheck
     dist(nc) = SQRT( (xBIT(nc)-xCell)**2 &
                    + (yBIT(nc)-yCell)**2 &
                    + (zBIT(nc)-zCell)**2 )
  ENDDO
  iDummy           = MINLOC(dist(1:nCheck))
  shortestProbe    = iDummy(1)

  xBI              = xBIT(shortestProbe)
  yBI              = yBIT(shortestProbe)
  zBI              = zBIT(shortestProbe)
  closestElement   = cElement(shortestProbe)
  distBIElem       = dist(shortestProbe)

  dotNorm = (xCell - xBI) * xBIN(shortestProbe) &
          + (yCell - yBI) * yBIN(shortestProbe) &
          + (zCell - zBI) * zBIN(shortestProbe)

  xBItable(iCell,jCell,kCell) = xBI
  yBItable(iCell,jCell,kCell) = yBI
  zBItable(iCell,jCell,kCell) = zBI
  cELtable(iCell,jCell,kCell) = closestElement
  is_BIset(iCell,jCell,kCell) = 1

    !Ye,debug
    !if ( (iCell.eq.1).and.(jCell.eq.56).and.(kCell.eq.41) )then
    !   write(302,*)ntime,nCheck,iBody
    !   write(302,*)xBItable(iCell,jCell,kCell),yBItable(iCell,jCell,kCell), &
    !               zBItable(iCell,jCell,kCell)
    !   write(302,*)cELtable(iCell,jCell,kCell)
    !   write(302,*)'---'
    !endif

! ============================================================================
!   Deallocate local array 
! ============================================================================

    DEALLOCATE(NeighElemInd,STAT=iErr)
    DEALLOCATE(distMarker, STAT=iErr)
    IF ( iErr /= ERR_NONE ) THEN
      WRITE(STDOUT,*) &
       'calc_bodyIntercept_Unstruc: Memory Deallocation Error for NeighElemInd'
      STOP
    ENDIF ! ierr

   END SUBROUTINE calc_bodyIntercept_Unstruc
!------------------------------------------------------------------------------
                            
!------------------------------------------------------------------------------
   SUBROUTINE check_BIInsideTriangle(xVert,yVert,zVert,xNorm,yNorm,zNorm,      &
                                     xBITemp,yBITemp,zBITemp,area123,areaDiff)
!------------------------------------------------------------------------------
           
    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays

    IMPLICIT NONE

!... parameters

    REAL(KIND=CGREAL), DIMENSION(3), INTENT(IN) :: xVert,yVert,zVert
    REAL(KIND=CGREAL), INTENT(IN)               :: xNorm,yNorm,zNorm

    REAL(KIND=CGREAL), INTENT(IN)               :: xBITemp,yBITemp,zBITemp
    REAL(KIND=CGREAL), INTENT(OUT)              :: area123,areaDiff

!... loop variables

    INTEGER :: iside

!... local variables

    REAL(KIND=CGREAL)        :: distanceToPlane,distPointToPlane
    REAL(KIND=CGREAL)        :: side12,side23,side31,side14,side24,side34
    REAL(KIND=CGREAL)        :: area124,area234,area314
    REAL(KIND=CGREAL)        :: semiPerimeter123,semiPerimeter124
    REAL(KIND=CGREAL)        :: semiPerimeter234,semiPerimeter314

! ******************************************************************************
! Check to see if normal intercept lies inside the closest trianglular element
!               3 
!               *  .
!              /  \   .
!             /    \    .
!            /      \    * 4=BI
!           /        \  .
!         1*__________*2
!         
! Basic Idea :  IF [ AREA(124) + AREA(234) + AREA(314) ] > AREA(123) THEN  POINT(4) is
! outside triangle (123)
!
! using Heron formula for area of triangle
! AREA(123) = SQRT[ S * ( S - S12) * (S - S23) * (S - S31) ]
! S = 0.5*(S12 + S23 + S31) 
! ******************************************************************************

     side12 =  SQRT( (xVert(2)-xVert(1))**2  &
                    +(yVert(2)-yVert(1))**2  &
                    +(zVert(2)-zVert(1))**2  )
     side23 =  SQRT( (xVert(3)-xVert(2))**2  &
                    +(yVert(3)-yVert(2))**2  &
                    +(zVert(3)-zVert(2))**2  )
     side31 =  SQRT( (xVert(1)-xVert(3))**2  &
                    +(yVert(1)-yVert(3))**2  &
                    +(zVert(1)-zVert(3))**2  )
     side14 =  SQRT( (xBITemp-xVert(1))**2  &
                    +(yBITemp-yVert(1))**2  &
                    +(zBITemp-zVert(1))**2  )
     side24 =  SQRT( (xBITemp-xVert(2))**2  &
                    +(yBITemp-yVert(2))**2  &
                    +(zBITemp-zVert(2))**2  )
     side34 =  SQRT( (xBITemp-xVert(3))**2  &
                    +(yBITemp-yVert(3))**2  &
                    +(zBITemp-zVert(3))**2  )

     semiPerimeter123 = 0.5_CGREAL*(side12 + side23 + side31)
     semiPerimeter124 = 0.5_CGREAL*(side12 + side24 + side14)
     semiPerimeter234 = 0.5_CGREAL*(side23 + side24 + side34)
     semiPerimeter314 = 0.5_CGREAL*(side31 + side34 + side14)

     area123       = SQRT( semiPerimeter123*abs(semiPerimeter123-side12) &
                                           *abs(semiPerimeter123-side23) &
                                           *abs(semiPerimeter123-side31)   )
   
     area124       = SQRT( semiPerimeter124*abs(semiPerimeter124-side12) &
                                           *abs(semiPerimeter124-side24) &
                                           *abs(semiPerimeter124-side14)   )
    
     area234       = SQRT( semiPerimeter234*abs(semiPerimeter234-side23) &
                                           *abs(semiPerimeter234-side24) &
                                           *abs(semiPerimeter234-side34)   )

     area314       = SQRT( semiPerimeter314*abs(semiPerimeter314-side31) &
                                           *abs(semiPerimeter314-side34) &
                                           *abs(semiPerimeter314-side14)   )

     areaDiff  = area124 + area234 + area314 - area123

   END SUBROUTINE check_BIInsideTriangle 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
   SUBROUTINE calc_BIOutsideTriangle(xVert,yVert,zVert,xCG,yCG,zCG,     &
                                    xVertN,yVertN,zVertN,               &
                                    xBITemp,yBITemp,zBITemp,distBIElem,xBInorm,yBInorm,zBInorm)
!------------------------------------------------------------------------------

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE MPI_module
    
    IMPLICIT NONE

!... parameters

    REAL(KIND=CGREAL), DIMENSION(3), INTENT(IN) :: xVert,yVert,zVert
    REAL(KIND=CGREAL), DIMENSION(3), INTENT(IN) :: xVertN,yVertN,zVertN
    REAL(KIND=CGREAL), INTENT(IN)               :: xCG,yCG,zCG
    REAL(KIND=CGREAL), INTENT(INOUT)  :: xBITemp,yBITemp,zBITemp
    REAL(KIND=CGREAL), INTENT(OUT) :: distBIElem
    REAL(KIND=CGREAL), INTENT(OUT) :: xBInorm,yBInorm,zBInorm

!... loop variables

    INTEGER :: iside

!... local variables

    INTEGER :: isideSelect,mside,kside
    INTEGER :: nodeSelect1,nodeSelect2

    REAL(KIND=CGREAL) :: aCrossbVectMagn,dotVal
    REAL(KIND=CGREAL) :: distIntBI,distIntCG,distIntMin
    REAL(KIND=CGREAL) :: distNorm,distNode1BINorm,distNode2BINorm
    REAL(KIND=CGREAL) :: distVert1BI,distVert2BI
    REAL(KIND=CGREAL) :: magnitude12,magnitudeBICG,magnitude,projectedLength
    REAL(KIND=CGREAL) :: node12x,node12y,node12z
    REAL(KIND=CGREAL) :: vec01x,vec01y,vec01z
    REAL(KIND=CGREAL) :: vec12x,vec12y,vec12z
    REAL(KIND=CGREAL), DIMENSION(3) :: aVect,bVect,cVect,aCrossbVect,cCrossbVect
    REAL(KIND=CGREAL), DIMENSION(3) :: xInt,yInt,zInt
    REAL(KIND=CGREAL), DIMENSION(3) :: vect1,vect2,vect3,vect4
    REAL(KIND=CGREAL), DIMENSION(3,3) :: vectInt
    REAL(kind=CGREAL) :: epstol
    REAL(kind=CGREAL) :: xBI, yBI, zBI
!******************************************************************************
!
! If the projection, (xBITemp, yBITemp, zBITemp), falls outside a triangle,
! find the closest point on the edge (including the vertex) to the projection.   
!
!             x2
!             *
!             |
!             |
!    x3       | Int    x4  (BITemp)
!    *--------*-------*
!             |
!             | iside = 1  
!             | 
!             * x1
!
! -- Haoxiang Luo

    vect4(1:3) = (/xBITemp,yBITemp,zBITemp/)

    DO iside = 1,3
      mside = iside +1
      IF(iside == 3) mside = 1

      if(iside == 1) kside = 3
      if(iside == 2) kside = 1
      if(iside == 3) kside = 2

      vect1(1:3) =(/xVert(iside),yVert(iside),zVert(iside)/)
      vect2(1:3) =(/xVert(mside),yVert(mside),zVert(mside)/)
      vect3(1:3) =(/xVert(kside),yVert(kside),zVert(kside)/)

      aVect(1:3) = vect3(1:3)-vect1(1:3)  !vec13      
      bVect(1:3) = vect2(1:3)-vect1(1:3)  !vec12
      cVect(1:3) = vect4(1:3)-vect1(1:3)  !vec14

      call calc_crossProduct(aVect,bVect,aCrossbVect)
      call calc_crossProduct(cVect,bVect,cCrossbVect)

      dotVal = DOT_PRODUCT(cCrossbVect,aCrossbVect)
      if(dotVal <= 0.0d0) then ! the BItemp and a vertex are two different sides
        isideSelect = iside
        go to 111
     endif

    END DO ! iside 
      
111 continue

! ============================================================================
!   Trap error for isideSelect 
! ============================================================================
   
     IF ( isideSelect < 0 ) THEN
      WRITE(STDOUT,*) &
       'calc_BIdistMin: Incorrect selection of iside (Should be either 1, 2 or 3'
      WRITE(STDOUT,*) &
       '                default value selected = ',isideSelect
      WRITE(STDOUT,*) xVert(1),yVert(1),zVert(1)
      WRITE(STDOUT,*) xVert(2),yVert(2),zVert(2)
      WRITE(STDOUT,*) xVert(3),yVert(3),zVert(3)
      WRITE(STDOUT,*) vect1(1:3)
      WRITE(STDOUT,*) vect2(1:3)
      call write_subdomain()

      STOP
     END IF ! isideSelect 

! ============================================================================
!   Select appropriate vertices from isideSelect 
! ============================================================================
   
     SELECT CASE(isideSelect)
       CASE(1)
         nodeSelect1 = 1
         nodeSelect2 = 2
       CASE(2)
         nodeSelect1 = 2
         nodeSelect2 = 3
       CASE(3)
         nodeSelect1 = 3
         nodeSelect2 = 1
     END SELECT ! isideSelect

! ============================================================================
!   Drop normals from BI to selected side 
!    and find coordinates of intersection point
!
!   unit vector from node 1 to node 2
!   use formula for distance between point and line from Mathworld
!   http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
! 
!    x1                  x2
!    *-------------------*
!             |                              |(x2-x1) x (x1-x0)|
!             | d                       d =  --------------------
!             |                                   |x2-x1|
!             * x0
!   x0: BI
! ============================================================================


! ============================================================================
!    Project vector BI-node1 onto node12 to find body intercept point
! ============================================================================

     node12x = xVert(nodeSelect2) - xVert(nodeSelect1)
     node12y = yVert(nodeSelect2) - yVert(nodeSelect1)
     node12z = zVert(nodeSelect2) - zVert(nodeSelect1)

     magnitude = SQRT(node12x**2 + node12y**2 + node12z**2)

     node12x = node12x/magnitude
     node12y = node12y/magnitude
     node12z = node12z/magnitude

     projectedLength = (xBITemp - xVert(nodeSelect1))*node12x  &
                      +(yBITemp - yVert(nodeSelect1))*node12y  &
                      +(zBITemp - zVert(nodeSelect1))*node12z

     xBI = xVert(nodeSelect1) + projectedLength*node12x
     yBI = yVert(nodeSelect1) + projectedLength*node12y
     zBI = zVert(nodeSelect1) + projectedLength*node12z

! ============================================================================
!    Determine whether the BI is in between the two nodes or is outside of the
!    segment.
! ============================================================================

      distNode1BINorm = SQRT( (xVert(nodeSelect1)-xBI)**2 &
                            + (yVert(nodeSelect1)-yBI)**2 &
                            + (zVert(nodeSelect1)-zBI)**2 )/magnitude

      distNode2BINorm = SQRT( (xVert(nodeSelect2)-xBI)**2 &
                            + (yVert(nodeSelect2)-yBI)**2 &
                            + (zVert(nodeSelect2)-zBI)**2 )/magnitude

      IF ( distNode1BINorm <= 1.0_CGREAL .AND. distNode2BINorm <= 1.0_CGREAL)THEN
         !The edge projection point is between two vertices of the side

         xBITemp = xBI
         yBITemp = yBI
         zBITemp = zBI

         xBInorm = 0.5_CGREAL*(xVertN(nodeSelect1) + xVertN(nodeSelect2))
         yBInorm = 0.5_CGREAL*(yVertN(nodeSelect1) + yVertN(nodeSelect2))
         zBInorm = 0.5_CGREAL*(zVertN(nodeSelect1) + zVertN(nodeSelect2))

      ELSEIF(distNode2BINorm > 1.0_CGREAL ) THEN
         !The edge projection is closest to node 1
         xBITemp = xVert(nodeSelect1)
         yBITemp = yVert(nodeSelect1)
         zBITemp = zVert(nodeSelect1)

         xBInorm = xVertN(nodeSelect1)
         yBInorm = yVertN(nodeSelect1)
         zBInorm = zVertN(nodeSelect1)

      ELSE
         !The edge projection is closest to node 2
         xBITemp = xVert(nodeSelect2)
         yBITemp = yVert(nodeSelect2)
         zBITemp = zVert(nodeSelect2)

         xBInorm = xVertN(nodeSelect2)
         yBInorm = yVertN(nodeSelect2)
         zBInorm = zVertN(nodeSelect2)
      END IF ! distNode1BINorm

      distBIElem  =  SQRT( (xCG - xBITemp)**2 &
                         + (yCG - yBITemp)**2 &
                         + (zCG - zBITemp)**2 )

   END SUBROUTINE calc_BIOutsideTriangle 
!------------------------------------------------------------------------------

   SUBROUTINE calc_crossProduct(r,s,cross_product)
   
    USE global_parameters
    IMPLICIT NONE

    REAL(KIND=CGREAL), DIMENSION(3), INTENT(IN)  :: r,s
    REAL(KIND=CGREAL), DIMENSION(3), INTENT(OUT) :: cross_product 

    INTEGER :: component,i,j

    DO component = 1,3
      i = MODULO(component,3) + 1
      j = MODULO(i,3) + 1
      cross_product(component) = r(i)*s(j) - s(i)*r(j)
    END DO ! component 

   END SUBROUTINE calc_crossProduct
!------------------------------------------------------------------------------

