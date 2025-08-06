!---------------------------------------
!  SUBROUTINE open_probe_files()
!  SUBROUTINE read_probe_inputs()
!  SUBROUTINE write_probe_files()
!---------------------------------------
!
!---------------------------------------------

   SUBROUTINE open_probe_files()
    
    USE flow_parameters
    USE boundary_arrays
    USE grid_arrays
    
    IMPLICIT NONE

    CHARACTER*20  :: probeFile
    CHARACTER*30  :: inProbeFile

    INTEGER :: m

    !Open files for body probe
    DO m = 1, nProbeBody
      probeFile = TRIM("probe_marker_out")
      WRITE(inProbeFile,101) trim(probeFile),m

      IF (nread==0) THEN
         OPEN(UNIT=ifuMarkerOut+m-1,FILE=inProbeFile,FORM='formatted',ACTION="WRITE")

      ELSE ! append for restart simulations
         OPEN(UNIT=ifuMarkerOut+m-1,FILE=inProbeFile,FORM='formatted',POSITION='append',ACTION="WRITE")
      ENDIF
    END DO ! m

101  FORMAT(a,'_',i3.3,'.dat')

   END SUBROUTINE open_probe_files

!---------------------------------------------
   SUBROUTINE read_probe_inputs()

    USE flow_parameters
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER :: m

    OPEN(ifuProbeIn,FILE='probe_body_in.dat',STATUS='UNKNOWN')
    READ(ifuProbeIn,*) nProbeBody
    PRINT*,'   nProbeBody = ',nProbeBody

    ALLOCATE(markerProbe(nProbeBody,2))
    
    PRINT*,'Reading probe_flow_in.dat ...'
    DO m= 1, nProbeBody
       !read in body ID and marker ID
       READ(ifuProbeIn,*) markerProbe(m,1), markerProbe(m,2)   
    ENDDO ! m
    CLOSE(ifuProbeIn)

   END SUBROUTINE read_probe_inputs

!---------------------------------------------
   SUBROUTINE write_probe_files()

    USE global_parameters
    USE flow_parameters
    USE flow_arrays
    USE boundary_arrays
    
    IMPLICIT NONE
    
    INTEGER :: i,j,k,m,iBody,n
    REAL(KIND=CGREAL) :: xP,yP,zP,uP,vP,wP,fxP,fyP,fzP,sxP,syP,szP
    REAL(KIND=CGREAL) :: fxP1,fyP1,fzP1,uP1,vP1,wP1    

    print*, "in write_probe_files", nProbeBody
    DO m= 1, nProbeBody

      iBody = markerProbe(m,1)
      n     = markerProbe(m,2)

      xP = xBodyMarker(n,iBody)
      yP = yBodyMarker(n,iBody)
      zP = zBodyMarker(n,iBody)

      !Ye, filtered BM velocity
      uP = uBodyMarker(n,iBody)
      vP = vBodyMarker(n,iBody)
      wP = wBodyMarker(n,iBody)

      !Ye, unfiltered BM velocity
      uP1 = u1BodyMarker(n,iBody)
      vP1 = v1BodyMarker(n,iBody)
      wP1 = w1BodyMarker(n,iBody) 

      !Ye, filtered BM force
      fxP= xMarkerForce(n,iBody)
      fyP= yMarkerForce(n,iBody)
      fzP= zMarkerForce(n,iBody)

      !Ye, unfiltered BM force
      fxP1= xMarkerForce0(n,iBody)
      fyP1= yMarkerForce0(n,iBody)
      fzP1= zMarkerForce0(n,iBody)

      sxP= xMarkerStress(n,iBody)
      syP= yMarkerStress(n,iBody)
      szP= zMarkerStress(n,iBody)

      WRITE(ifuMarkerOut+m-1,1001) time,xP,yP,zP,uP1,vP1,wP1,uP,vP,wP, &
                                   sxP,syP,szP,fxP1,fyP1,fzP1,fxP,fyP,fzP

    ENDDO ! m  

1001 FORMAT(F16.5, 19E15.5)

   END SUBROUTINE write_probe_files
!-------------------------------------------   
