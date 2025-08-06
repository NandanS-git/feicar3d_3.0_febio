!------------------------------------------------------------------------------
   SUBROUTINE initialize_marker()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays   
    USE usr_module
    USE MPI_module
    
    IMPLICIT NONE

!...Loop variables
    INTEGER           :: i,j,iBody
    INTEGER           :: k, nnt

!...Local variables
    INTEGER           :: m, nBodyMarkerIn,nBodyRstrt,n, &
                         nPtsBodyMarkerIn, totNumTriElemIn
    INTEGER, DIMENSION(nBody) :: body_type_orig

    REAL(KIND=CGREAL) :: dMin,                                      &
                         bodyResolution, bodyResolutionNormalized,  &
                         rad,theta,phi, xTemp, yTemp, zTemp
    REAL(KIND=CGREAL) :: cosalpha,sinalpha

    ! by Nandan for reading unstruc mesh from Gmsh
    INTEGER                :: numTags, i1
    CHARACTER (LEN = 72)   :: cLine
    INTEGER, allocatable   :: ibElTagVal(:)
!------------------------------------------------ 
    
    ! Set read Marker flag

    ! PRINT*,'SETTING UP CANONICAL BODIES in INITIALIZE_MARKER '
     
    IF ( nread == 1) GOTO 1000
    
    ! Save copy of canonical body type
    body_type_orig(1:nBody) = canonical_body_type(1:nBody)

    ! Initialize for  no restart
    uBodyMarker = 0.0_CGREAL
    vBodyMarker = 0.0_CGREAL
    wBodyMarker = 0.0_CGREAL

    xBodyMarker = 0.0_CGREAL
    yBodyMarker = 0.0_CGREAL
    zBodyMarker = 0.0_CGREAL

    !OPEN(ifuUnstrucSurfIn, FILE='unstruc_surface_in.dat',STATUS='UNKNOWN')
    OPEN(ifuUnstrucSurfIn, FILE='cylinderPlateTransfiniteSurface.msh',STATUS='UNKNOWN')
    OPEN(ifuMarkerIn,      FILE='marker2D_in.dat',       STATUS='UNKNOWN')

    ! Select appropriate body type
    ! Setting up marker locations

    DO iBody = 1, nBody

      cosalpha = cos(angz(iBody))
      sinalpha = sin(angz(iBody))

      SELECT CASE (canonical_body_type(iBody))

      CASE(ELLIPTIC_CYLINDER)
         PRINT*,'  SETTING UP ELLIPTIC CYLINDER'

         !-- Test resolution for GCM
         IF (boundary_formulation == GCM_METHOD) THEN
             dxMin          = MINVAL(dx(1:nx-1))
             dyMin          = MINVAL(dy(1:ny-1))
             dMin           = MIN(dxMin,dyMin)
             bodyResolution = PI*( radiusx(iBody)+radiusy(iBody) ) /   &
                              REAL(nPtsBodyMarker(iBody),KIND=CGREAL)
             bodyResolutionNormalized = dMin/bodyResolution

             PRINT*,'    dx Min ',dxMin
             PRINT*,'    dy Min ',dyMin
             PRINT*,'    d  Min ',dMin
             PRINT*,'    Body Resolution',bodyResolution
             PRINT*,'    Current Normalized Resolution for Body ',bodyResolutionNormalized
             IF (bodyResolutionNormalized < 2.0_CGREAL ) THEN
               PRINT*,'    Ideal Normalized Resolution Should be at LEAST 2 .. aborting'
               STOP
             ENDIF
          ENDIF

          DO m = 1, nPtsBodyMarkerOrig(iBody)
             theta  = (REAL(m,KIND=CGREAL)-1.0_CGREAL)*2.0_CGREAL*PI/ &
                       REAL(nPtsBodyMarkerOrig(iBody),KIND=CGREAL)        
             !xTemp                = radiusx(iBody)*COS(theta)
             !yTemp                = radiusy(iBody)*SIN(theta)
             xBodyMarker(m,iBody)  = radiusx(iBody)*COS(theta)
             yBodyMarker(m,iBody)  = radiusy(iBody)*SIN(theta)
            ! xBodyMarker(m,iBody) = xcent(iBody) + xTemp*cosalpha(iBody) - yTemp*sinalpha(iBody) 
            ! yBodyMarker(m,iBody) = ycent(iBody) + xTemp*sinalpha(iBody) + yTemp*cosalpha(iBody) 
             zBodyMarker(m,iBody)  = 0.0_CGREAL  !z(1)
          ENDDO ! m

          CALL extend_cylinder_3D(iBody)

       CASE(GENERAL_CYLINDER)
          PRINT*,'  SETTING UP GENERAL CYLINDER'
          READ(ifuMarkerIn,*) nBodyMarkerIn
          IF ( nBodyMarkerIn /= nPtsBodyMarkerOrig(iBody) ) THEN
             PRINT*,'Init_Marker: Inconsistent body_in.dat and marker_in.dat files for body = ', iBody
             PRINT*,'             Reading in body_in.dat nPtsBodyMarker = ',nPtsBodyMarkerOrig(iBody)
             PRINT*,'             Reading from marker_in.dat nBodyMarkerIn = ', nBodyMarkerIn
             STOP
          ENDIF ! nBodyMarkerIn
            
          PRINT*,'  iBody ', iBody
          PRINT*,'     nBodyMarkerIn ', nBodyMarkerIn
          PRINT*,'     nPtsBodyMarker ', nPtsBodyMarker(iBody)
          DO m = 1, nPtsBodyMarkerOrig(iBody)
             READ(ifuMarkerIn,*)k, xBodyMarker(m,iBody),yBodyMarker(m,iBody)
             zBodyMarker(m,iBody) = z(1)
          ENDDO

          CALL extend_cylinder_3D(iBody)
           
       CASE(ELLIPSOID)

          PRINT*,'  SETTING UP ELLIPSOID'
          m = 0 
          DO i = 1, n_phi(iBody)
             phi  = ( REAL(i,KIND=CGREAL)-1.0_CGREAL )*PI/ &
                      (REAL(n_phi(iBody),KIND=CGREAL)-1.0_CGREAL)
             IF (i .EQ. 1 .OR. i .EQ. n_phi(iBody)) THEN
                  theta = 0.0_CGREAL
                  m = m + 1
 
                  xTemp                = radiusx(iBody)*COS(theta)*SIN(phi)
                  yTemp                = radiusy(iBody)*SIN(theta)*SIN(phi)
                  xBodyMarker(m,iBody) = xcent(iBody) + xTemp*cosalpha - yTemp*sinalpha
                  yBodyMarker(m,iBody) = ycent(iBody) + xTemp*sinalpha + yTemp*cosalpha
                  zBodyMarker(m,iBody) = zcent(iBody) + radiusz(iBody)*COS(phi)
             ELSE
               DO j = 1,n_theta(iBody)
                 theta = (REAL(j,KIND=CGREAL)-1.0_CGREAL)*2.0_CGREAL*PI/ &
                          REAL(n_theta(iBody),KIND=CGREAL)
                  m = m + 1

                  xTemp                = radiusx(iBody)*COS(theta)*SIN(phi)
                  yTemp                = radiusy(iBody)*SIN(theta)*SIN(phi)
                  xBodyMarker(m,iBody) = xcent(iBody) + xTemp*cosalpha - yTemp*sinalpha
                  yBodyMarker(m,iBody) = ycent(iBody) + xTemp*sinalpha + yTemp*cosalpha
                  zBodyMarker(m,iBody) = zcent(iBody) + radiusz(iBody)*COS(phi)
               ENDDO ! j
             END IF
           ENDDO ! i

           nPtsBodyMarker(iBody) = m
           !write(*,*) 'Total number of marker points are: ', m
           j = 0
           i = 1

           DO k = 1, n_theta(iBody)
                 j = j + 1
                 triElemNeig(1,j,iBody) = i
                 triElemNeig(2,j,iBody) = i + k
                 IF (i+k .EQ. n_theta(ibody) + 1) THEN
                    triElemNeig(3,j,iBody) = i + 1
                 ELSE 
                    triElemNeig(3,j,iBody) = i + k + 1
                 END IF
           END DO ! end k

           nnt = 1

           DO i = 2, m-n_theta(iBody)-1
                 j = j + 1
                 triElemNeig(1,j,iBody) = i
                 triElemNeig(2,j,iBody) = i + n_theta(iBody)
                 !write (*,*) 'i, n_theta =', i, nnt*n_theta(ibody)+1

                 IF ( i .EQ. nnt*n_theta(ibody)+1 ) THEN
                   triElemNeig(3,j,iBody) = i + 1                 
                 ELSE
                   triElemNeig(3,j,iBody) = i + n_theta(iBody) + 1
                 END IF

                 j = j + 1 
                 triElemNeig(1,j,iBody) = i
                 IF ( i .EQ. nnt*n_theta(ibody) + 1 ) THEN
                   triElemNeig(2,j,iBody) = i + 1
                   triElemNeig(3,j,iBody) = i + 1 - n_theta(iBody)
                    nnt = nnt + 1                   
                 ELSE
                   triElemNeig(2,j,iBody) = i + n_theta(iBody) + 1
                   triElemNeig(3,j,iBody) = i + 1
                 END IF     
            ENDDO ! i

           DO i = m-n_theta(iBody), m-1
 
                 j = j + 1
                 triElemNeig(1,j,iBody) = i
                 triElemNeig(2,j,iBody) = m
                 IF (i+1 .EQ. m) THEN
                   triElemNeig(3,j,iBody) = m-n_theta(iBody) 
                 ELSE
                   triElemNeig(3,j,iBody) = i + 1
                 END IF
           END DO ! end i

          totNumTriElem(iBody) = j
          !write(*,*) 'Total # of elements: ', j

       CASE(UNSTRUCTURED_SURFACE)
          !PRINT*,'  SETTING UP UNSTRUCTURED SURFACE for Body #',iBody

         !  READ(ifuUnstrucSurfIn,*)
         !  READ(ifuUnstrucSurfIn,*)
         !  READ(ifuUnstrucSurfIn,*) nPtsBodyMarkerIn, totNumTriElemIn
         !  READ(ifuUnstrucSurfIn,*)
         !  IF ( nPtsBodyMarkerIn /= nPtsBodyMarker(iBody) ) THEN
         !     PRINT*,'Init_Marker: Inconsistent canonical_body_in.dat and unstruc_surface_in.dat files for body = ', iBody
         !     PRINT*,'             Reading in canonical_body_in.dat    nPtsBodyMarker   = ',nPtsBodyMarker(iBody)
         !     PRINT*,'             Reading from unstruc_surface_in.dat nPtsBodyMarkerIn = ',nPtsBodyMarkerIn
         !     STOP
         !  ENDIF ! nPtsBodyMarkerIn
             
         !  IF ( totNumTriElemIn /= totNumTriElem(iBody) ) THEN
         !     PRINT*,'Init_Marker: Inconsistent canonical_body_in.dat and unstruc_surface_in.dat files for body = ', iBody
         !     PRINT*,'             Reading in canonical_body_in.dat     totNumTriElem   = ', totNumTriElem(iBody)
         !     PRINT*,'             Reading from unstruc_surface_in.dat  totNumTriElemIn = ', totNumTriElemIn
         !  ENDIF !  totNumTriElemIn
 
         !  DO m=1,nPtsBodyMarker(iBody)
         !     !READ(ifuUnstrucSurfIn,*) i,xBodyMarker(m,iBody),yBodyMarker(m,iBody),zBodyMarker(m,iBody)

         !     READ(ifuUnstrucSurfIn,*) i, xTemp, yTemp, zTemp

         !     xBodyMarker(m,iBody) = xTemp*zoom_factor(iBody)
         !     yBodyMarker(m,iBody) = yTemp*zoom_factor(iBody)
         !     zBodyMarker(m,iBody) = zTemp*zoom_factor(iBody)
         !  ENDDO

         !  READ(ifuUnstrucSurfIn,*)
         !  DO  j=1,totNumTriElem(iBody)
         !     READ(ifuUnstrucSurfIn,*) i,triElemNeig(1,j,iBody),triElemNeig(2,j,iBody),triElemNeig(3,j,iBody)
         !  ENDDO
         !  READ(ifuUnstrucSurfIn,*)
         !  READ(ifuUnstrucSurfIn,*)pointOutsideBodyX(iBody),pointOutsideBodyY(iBody),pointOutsideBodyZ(iBody)

         ! by Nandan
         DO i=1,4
            READ(ifuUnstrucSurfIn,*)
         END DO

         READ(ifuUnstrucSurfIn,*) numTags
         print*, 'numTags', numTags
         ALLOCATE ( ibElTagVal(numTags) )

         DO i = 1, numTags
           READ (ifuUnstrucSurfIn,*) i1, ibElTagVal(n), cLine
           WRITE (6,*) i1, ibElTagVal(n), cLine
         END DO

         DO n = 1, 2
           READ (ifuUnstrucSurfIn,*) cLine
         END DO

         READ(ifuUnstrucSurfIn,*) nPtsBodyMarkerIn

         IF ( nPtsBodyMarkerIn /= nPtsBodyMarker(iBody) ) THEN
            PRINT*,'Init_Marker: Inconsistent canonical_body_in.dat and unstruc_surface_in.dat files for body = ', iBody
            PRINT*,'             Reading in canonical_body_in.dat    nPtsBodyMarker   = ',nPtsBodyMarker(iBody)
            PRINT*,'             Reading from unstruc_surface_in.dat nPtsBodyMarkerIn = ',nPtsBodyMarkerIn
            STOP
         ENDIF ! nPtsBodyMarkerIn

         DO n = 1, nPtsBodyMarkerIn
            READ (ifuUnstrucSurfIn,*) i1, xBodyMarker(n,iBody), yBodyMarker(n,iBody), zBodyMarker(n,iBody)
         END DO

         DO n = 1, 2
            READ (ifuUnstrucSurfIn,*)
         END DO  

         READ(ifuUnstrucSurfIn,*) totNumTriElemIn

         IF ( totNumTriElemIn /= totNumTriElem(iBody) ) THEN
            PRINT*,'Init_Marker: Inconsistent canonical_body_in.dat and unstruc_surface_in.dat files for body = ', iBody
            PRINT*,'             Reading in canonical_body_in.dat     totNumTriElem   = ', totNumTriElem(iBody)
            PRINT*,'             Reading from unstruc_surface_in.dat  totNumTriElemIn = ', totNumTriElemIn
         ENDIF !  totNumTriElemIn

         DO n = 1, totNumTriElemIn
            READ (ifuUnstrucSurfIn,*) i1, i1, i1, i1, i1, triElemNeig(1,n,iBody), triElemNeig(2,n,iBody), triElemNeig(3,n,iBody)
         END DO

         open (101, file='pointOutsideBody.dat')
            READ(101,*) pointOutsideBodyX(iBody),pointOutsideBodyY(iBody),pointOutsideBodyZ(iBody)
         close(101)

       END SELECT ! canonical_body_type
    ENDDO ! iBody

    !--Write surface mesh data in Tecplot Format to check

    if(lProc .ne. PROC_M) goto 900  ! skip output for other processors
    
    OPEN(ifuMarkerOut,     FILE='marker_unstruc_out.dat', STATUS='UNKNOWN')   
    OPEN(ifuBodyOut,       FILE='canonical_body_out.dat', STATUS='UNKNOWN')
    OPEN(ifuUnstrucSurfOut,FILE='unstruc_surface_out.dat',STATUS='UNKNOWN')

    DO iBody = 1, nBody
       WRITE(ifuMarkerOut,*) 'TITLE="3D TRIANGULAR SURFACE DATA"'
       WRITE(ifuMarkerOut,*) 'VARIABLES="X","Y","Z"'
       WRITE(ifuMarkerOut,*) 'ZONE N=',nPtsBodyMarker(iBody),',E=',totNumTriElem(iBody),'F=FEPOINT, ET=TRIANGLE'
       DO m=1,nPtsBodyMarker(iBody)
          WRITE(ifuMarkerOut,'(3F15.7)') xBodyMarker(m,iBody),yBodyMarker(m,iBody),zBodyMarker(m,iBody)
       ENDDO
       DO j=1,totNumTriElem(iBody)
          WRITE(ifuMarkerOut,'(3I10)') triElemNeig(1,j,iBody),triElemNeig(2,j,iBody),triElemNeig(3,j,iBody)
       ENDDO
    ENDDO ! iBody

    !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    ! set up initial marker velocities
    DO iBody = 1, nBody
       ! Zero for everything. Will be changed in the future
       DO m = 1, nPtsBodyMarkerOrig(iBody)
          uBodyMarker(m,iBody) =  0.0_CGREAL
          vBodyMarker(m,iBody) =  0.0_CGREAL
          wBodyMarker(m,iBody) =  0.0_CGREAL
       ENDDO
    ENDDO ! iBody
  
!----------------------
! Write out the body information
! and marker points
!----------------------

    WRITE(ifuBodyOut,*)nbody

    ! write body information to Output file
    DO iBody = 1, nBody

      ! writing out body parameter file
      IF(canonical_body_type(iBody) <= GENERAL_CYLINDER) THEN
         canonical_body_type(iBody) = 4
         body_dim(iBody)            = 2
         radiusz(iBody)             = 0.0_CGREAL
         zcent(iBody)             = 0.0_CGREAL
     
      ELSE IF (canonical_body_type(iBody) == ELLIPSOID .AND. &
              boundary_formulation == GCM_METHOD) THEN
         canonical_body_type(iBody) = 4
      ENDIF

      WRITE(ifuBodyOut,'(4I10)')canonical_body_type(iBody),body_dim(iBody),unstruc_surface_type(iBody), &
                                boundary_motion_type(iBody)
      WRITE(ifuBodyOut,*)
      WRITE(ifuBodyOut,'(2I10)')wall_type(iBody), forced_motion_spec(iBody)
      WRITE(ifuBodyOut,'(2I20)')nPtsBodyMarker(iBody),totNumTriElem(iBody)
      WRITE(ifuBodyOut,'(3E15.7)')radiusx(iBody),radiusy(iBody),radiusz(iBody)
      WRITE(ifuBodyOut,'(3E15.7)')xcent(iBody),ycent(iBody),zcent(iBody)
      WRITE(ifuBodyOut,'(3E15.7)')angzinit(iBody)
      WRITE(ifuBodyOut,'(3E15.7)')vxcentTrans(iBody),vycentTrans(iBody),vzcentTrans(iBody)
      WRITE(ifuBodyOut,'(3E15.7)')ampx(iBody),ampy(iBody),ampz(iBody)
      WRITE(ifuBodyOut,'(3E15.7)')freqx(iBody),freqy(iBody),freqz(iBody)

      WRITE(ifuBodyOut,'(3E15.7)')angvx(iBody),angvy(iBody),angvz(iBody)
      WRITE(ifuBodyOut,'(3E15.7)')angxphase(iBody),angyphase(iBody),angzphase(iBody)
      WRITE(ifuBodyOut,'(3E15.7)')ampangx(iBody),ampangy(iBody),ampangz(iBody)
      WRITE(ifuBodyOut,'(3E15.7)')freqangx(iBody),freqangy(iBody),freqangz(iBody)
      WRITE(ifuBodyOut,'(3E15.7)')density_fluid, density_solid(iBody)
      WRITE(ifuBodyOut,'(3E15.7)')xcentConstr(iBody),ycentConstr(iBody),zcentConstr(iBody) 

      ! writing out unstructured surface

      WRITE(ifuUnstrucSurfOut,*)
      WRITE(ifuUnstrucSurfOut,*) 'Body #', iBody
      WRITE(ifuUnstrucSurfOut,*) nPtsBodyMarker(iBody), totNumTriElem(iBody)
      WRITE(ifuUnstrucSurfOut,*)
      DO m=1,nPtsBodyMarker(iBody)
         WRITE(ifuUnstrucSurfOut,'(I6,1X,3F15.7)') m,xBodyMarker(m,iBody),yBodyMarker(m,iBody),zBodyMarker(m,iBody)
      ENDDO
      WRITE(ifuUnstrucSurfOut,*) 'Triangular mesh'
      DO  j=1,totNumTriElem(iBody)
        WRITE(ifuUnstrucSurfOut,'(4I10)') j,triElemNeig(1,j,iBody),triElemNeig(2,j,iBody),triElemNeig(3,j,iBody)
      ENDDO
      WRITE(ifuUnstrucSurfOut,*) 'A point outside'
      ! WRITE(ifuUnstrucSurfOut,'(3F15.7)')pointOutsideBodyX(iBody),pointOutsideBodyY(iBody),pointOutsideBodyZ(iBody)
      WRITE(ifuUnstrucSurfOut,'(3F15.7)') minval(xBodyMarker(:,iBody))-1.0 ,minval(yBodyMarker(:,iBody))-1.0,  &
                                          minval(zBodyMarker(:,iBody))-1.0

    ENDDO  ! end iBody

    DO iBody = 1, nBody
      IF(body_type_orig(iBody) <= GENERAL_CYLINDER .OR. &
         body_type_orig(iBody) == ELLIPSOID ) THEN

    PRINT*,'######################CODE NEEDS TO BE RERUN ##########################################'
    PRINT*,'Body parametrs have been written out in       :        canonical_body_out.dat'
    PRINT*,'Cylinder Surface(s) have been written out in  :  unstructured_surface_out.dat'
    PRINT*,'Following needs to be done:'
    PRINT*,'(1) Rename    unstruc_surface_out.dat    TO  unstruc_surface_in.dat'
    PRINT*,'(2) Rename    canonical_body_out.dat     TO  canonical_body_in.dat'
    PRINT*,'(3) Run code again'
        STOP
      ENDIF
    ENDDO

    CLOSE(ifuMarkerOut)
    CLOSE(ifuUnstrucSurfOut)
    CLOSE(ifuBodyOut)

900 CONTINUE
    
    CLOSE(ifuMarkerIn)
    CLOSE(ifuUnstrucSurfIn)

1000 CONTINUE

    if(lProc==PROC_M) PRINT*, 'CALL calculate_arclength_norm_ds() in initialize_marker'
    CALL calculate_arclength_norm_ds()
    if(lProc==PROC_M) PRINT*, 'calculate_arclength_norm_ds() OK!'

    ! set up the membrane thickness threshold
    dxMin = minval(dx(1:nx-1))
    dyMin = minval(dy(1:ny-1))
    dzMin = minval(dz(1:nz-1))
    membrane_tkns = membrane_tkns_factor * sqrt(dxMin**2 + dyMin**2 + dzMin**2)

    if(lProc==PROC_M)then
    write(*,*),'Membrane thickness threshold is: ', membrane_tkns
    endif

    !Initialize [xyz]BodyMarkerOld
    xBodyMarkerOld = xBodyMarker
    yBodyMarkerOld = yBodyMarker
    zBodyMarkerOld = zBodyMarker

   END SUBROUTINE initialize_marker 
!-----------------------------------------------------------------

   ! extend cylinder to make pseudo-3D body

   SUBROUTINE extend_cylinder_3D(iBody)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER, INTENT(IN)           :: iBody

!...Loop variables
    INTEGER           :: i,j,k,m,mInd,mp, &
                         nPtsBodyMarkerIn, totNumTriElemIn

    DO k=2,nz
    DO m=1,nPtsBodymarkerOrig(iBody)
       mInd = (k-1)*nPtsBodymarkerOrig(iBody)+m
       xBodyMarker(mInd,iBody) = xBodyMarker(m,iBody)
       yBodyMarker(mInd,iBody) = yBodyMarker(m,iBody)
       zBodyMarker(mInd,iBody) = z(k)
    ENDDO
    ENDDO

      j = 0
      DO k = 1, nz-1
      DO m = 1, nPtsBodyMarkerOrig(iBody) - 1   ! 
                                                ! for non-closing surface 
        j = j+1

        mp = m+1
        IF (m == nPtsBodyMarkerOrig(iBody)) mp = 1
        triElemNeig(1,j,iBody) = (k-1)*nPtsBodymarkerOrig(iBody) + m
        triElemNeig(2,j,iBody) =  k*nPtsBodymarkerOrig(iBody)    + m
        triElemNeig(3,j,iBody) = (k-1)*nPtsBodymarkerOrig(iBody) + mp

        j = j+1
        triElemNeig(1,j,iBody) = (k-1)*nPtsBodymarkerOrig(iBody) + mp
        triElemNeig(2,j,iBody) =  k*nPtsBodymarkerOrig(iBody)    + m
        triElemNeig(3,j,iBody) =  k*nPtsBodymarkerOrig(iBody)    + mp

      ENDDO
      ENDDO
      totNumTriElem(iBody) = j   ! H. Luo

      WRITE(ifuMarkerOut,*) 'TITLE="3D TRIANGULAR SURFACE DATA"'
      WRITE(ifuMarkerOut,*) 'VARIABLES="X","Y","Z"'
      WRITE(ifuMarkerOut,*) 'ZONE N=',nPtsBodyMarker(iBody),',E=',totNumTriElem(iBody),'F=FEPOINT, ET=TRIANGLE'
      DO m=1,nPtsBodyMarker(iBody)
         WRITE(ifuMarkerOut,'(3F15.7)') xBodyMarker(m,iBody),yBodyMarker(m,iBody),zBodyMarker(m,iBody)
      ENDDO
      DO j=1,totNumTriElem(iBody)
        WRITE(ifuMarkerOut,'(3I10)') triElemNeig(1,j,iBody),triElemNeig(2,j,iBody),triElemNeig(3,j,iBody)
      ENDDO

      pointOutsideBodyX(iBody) =-10.0_CGREAL
      pointOutsideBodyy(iBody) =-10.0_CGREAL
      pointOutsideBodyz(iBody) =-10.0_CGREAL

  END SUBROUTINE extend_cylinder_3D

!------------------------------------------------------------------------------

! extend cylinder to make pseudo-3D body

   SUBROUTINE extend_cylinder_vel_3D(iBody)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    IMPLICIT NONE

    INTEGER, INTENT(IN)           :: iBody

!...Loop variables
    INTEGER           :: i,j,k,m,mInd

      DO k=2,nz
      DO m=1,nPtsBodymarkerOrig(iBody)
        mInd = (k-1)*nPtsBodymarkerOrig(iBody)+m
        uBodyMarker(mInd,iBody) = uBodyMarker(m,iBody)
        vBodyMarker(mInd,iBody) = vBodyMarker(m,iBody)
        wBodyMarker(mInd,iBody) = wBodyMarker(m,iBody) 
      ENDDO
      ENDDO

  END SUBROUTINE extend_cylinder_vel_3D

!------------------------------------------------------------------------------
