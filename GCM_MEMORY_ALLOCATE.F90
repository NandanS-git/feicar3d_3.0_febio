!------------------------------------------------------------------------
   SUBROUTINE BOUNDARY_allocate_memory()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE MPI_module
    
    IMPLICIT NONE

    INTEGER ::iBody

    ! arrays required by all internal boundaries
    ALLOCATE(       iblank(0:nx+1,yb1:yb2,zb1:zb2))
    ALLOCATE(    dead_cell(0:nx+1,yb1:yb2,zb1:zb2))
    ALLOCATE(ghostCellMark(0:nx+1,yb1:yb2,zb1:zb2))
    ALLOCATE(      bodyNum(0:nx+1,yb1:yb2,zb1:zb2))
    allocate(       igmark(0:nx+1,yb1:yb2,zb1:zb2))
    ALLOCATE(ghostCellIndex(0:nx+1,yb1:yb2,zb1:zb2))

    ALLOCATE(iblankUndecided(0:nx+1,yb1:yb2,zb1:zb2))
    ALLOCATE(iblankTemp     (0:nx+1,yb1:yb2,zb1:zb2))

    !YE,TEST
    ALLOCATE(    iblank_real(0:nx+1,yb1:yb2,zb1:zb2))
    iblank_real = 0.0

    !Ye,pure solid node for DC
    allocate(ps1(0:nx+1,yb1:yb2,zb1:zb2))
    allocate(ps2(0:nx+1,yb1:yb2,zb1:zb2))
    allocate(xyzDCIP(0:nx+1,yb1:yb2,zb1:zb2,0:2))
    ps1 = 0
    ps2 = 0
    xyzDCIP = 0.0

    !Ye,goa
    ALLOCATE(goa(0:nx+1,yb1:yb2,zb1:zb2))
    goa = 0

    igmark = 0
    ghostCellMark = 0
    bodyNum       = 0
    ghostCellIndex= 0

    iblank      = 0
    dead_cell   = 0

    iblankUndecided = 0
    iblankTemp      = 0

    ALLOCATE(    i_lbuff(0:nx+1,yb1:yb2,1:2))
    ALLOCATE(    i_rbuff(0:nx+1,yb1:yb2,1:2))
    i_lbuff = 0
    i_rbuff = 0

    !ALLOCATE(    j_lbuff(0:nx+1,zb1:zb2,1:2))
    !ALLOCATE(    j_rbuff(0:nx+1,zb1:zb2,1:2))
    !j_lbuff = 0
    !j_rbuff = 0
    
    !Arrays associated with projection
    ALLOCATE( xBItable(0:nx+1,yb1:yb2,zb1:zb2) )
    ALLOCATE( yBItable(0:nx+1,yb1:yb2,zb1:zb2) )
    ALLOCATE( zBItable(0:nx+1,yb1:yb2,zb1:zb2) )
    ALLOCATE( cELtable(0:nx+1,yb1:yb2,zb1:zb2) )
    ALLOCATE( is_BIset(0:nx+1,yb1:yb2,zb1:zb2) )

    xBItable = 0.0_CGREAL
    yBItable = 0.0_CGREAL
    zBItable = 0.0_CGREAL
    cELtable = 0
    is_BIset = 0    
    
    ! since elliptic & general cylinder will be converted into
    ! 3D unstruc surfaces, we need to determine memory requirement for these

    DO iBody = 1, nBody

      SELECT CASE (canonical_body_type(iBody))

      CASE(ELLIPTIC_CYLINDER:GENERAL_CYLINDER)
         nPtsBodyMarkerOrig(iBody)= nPtsBodyMarker(iBody)
         nPtsBodyMarker(iBody)    = nPtsBodyMarkerOrig(iBody)*nz
         totNumTriElem(iBody)     = 2*nPtsBodyMarkerOrig(iBody)*(nz-1)

      CASE(ELLIPSOID)
         
         IF (boundary_formulation == GCM_METHOD) THEN
            nPtsBodyMarkerOrig(iBody)= nPtsBodyMarker(iBody)
            totNumTriElem(iBody)     = 2*nPtsBodyMarker(iBody) + 5
         ELSE
            nPtsBodyMarkerOrig(iBody)= nPtsBodyMarker(iBody)
            totNumTriElem(iBody)     = 2*nPtsBodyMarker(iBody) + 5
         END IF

      CASE(UNSTRUCTURED_SURFACE)
         nPtsBodyMarkerOrig(iBody)= nPtsBodyMarker(iBody)
      END SELECT ! canonical_body_type
    ENDDO        ! iBody

    ! Allocate arrays associated with boundary markers and interpolation setup.
    ! Note that each processor has a full copy of boundary marker arrays, but each
    ! has only a subdomain copy of the ghost or dead cell arrays 
    CALL MARKER_allocate_memory()
    CALL UNSTRUC_allocate_memory()

    CALL GCM_allocate_static_arrays()
    !print*, 'lProc ', lProc, ' allocated boundary arrays.'

   END SUBROUTINE BOUNDARY_allocate_memory
!------------------------------------------------------------------------
   SUBROUTINE GCM_AllocateGhostCellArrays()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    
    IMPLICIT NONE

    INTEGER :: nFaceMax
    
    !iRowMax  = 8

! Allocate Memory for various arrays pertinent to GCM Ghost Cell arrays

    ALLOCATE(iGhost(1:nGhost))
    ALLOCATE(jGhost(1:nGhost)) 
    ALLOCATE(kGhost(1:nGhost)) 

    ALLOCATE(coeffGCMD(iRowMax,nGhost))
    ALLOCATE(coeffGCMN(iRowMax,nGhost))
    
    ALLOCATE(xBodyInterceptTang(nGhost))
    ALLOCATE(yBodyInterceptTang(nGhost))
    ALLOCATE(zBodyInterceptTang(nGhost))

    ALLOCATE(xBodyInterceptNorm(nGhost))
    ALLOCATE(yBodyInterceptNorm(nGhost))
    ALLOCATE(zBodyInterceptNorm(nGhost))
    
    ALLOCATE(xBodyIntercept(nGhost))
    ALLOCATE(yBodyIntercept(nGhost))
    ALLOCATE(zBodyIntercept(nGhost))

    ALLOCATE(uBodyIntercept(nGhost))
    ALLOCATE(vBodyIntercept(nGhost)) 
    ALLOCATE(wBodyIntercept(nGhost)) 
    ALLOCATE(pBodyIntercept(nGhost)) 
    ALLOCATE(dpdnBodyIntercept(nGhost)) 
    ALLOCATE(dpdtBodyIntercept(nGhost)) 
                        
    ALLOCATE(closestMarker(nGhost)) 
    ALLOCATE(closestMarkerRatio(nGhost)) 
    ALLOCATE(closestElementGC(nGhost))
          
    ALLOCATE(iCellIndex(nGhost))
    ALLOCATE(jCellIndex(nGhost))
    ALLOCATE(kCellIndex(nGhost))
 
    ALLOCATE(xImagePoint(nGhost))
    ALLOCATE(yImagePoint(nGhost)) 
    ALLOCATE(zImagePoint(nGhost)) 
    ALLOCATE(probeLength(nGhost))

    ALLOCATE(res_laplacian(nGhost))
    ALLOCATE(alphaGhost   (nGhost))

    ALLOCATE(ghostNScoeff(10,nGhost))

    !   Initialize arrays

    iGhost             = 0
    jGhost             = 0 
    kGhost             = 0

    coeffGCMD          = 0.0_CGREAL
    coeffGCMN          = 0.0_CGREAL

    xBodyInterceptTang = 0.0_CGREAL
    yBodyInterceptTang = 0.0_CGREAL
    zBodyInterceptTang = 0.0_CGREAL

    xBodyInterceptNorm = 0.0_CGREAL
    yBodyInterceptNorm = 0.0_CGREAL
    zBodyInterceptNorm = 0.0_CGREAL
      
    xBodyIntercept     = 0.0_CGREAL
    yBodyIntercept     = 0.0_CGREAL 
    zBodyIntercept     = 0.0_CGREAL 
      
    uBodyIntercept     = 0.0_CGREAL
    vBodyIntercept     = 0.0_CGREAL
    wBodyIntercept     = 0.0_CGREAL 

    pBodyIntercept     = 0.0_CGREAL 
    dpdnBodyIntercept  = 0.0_CGREAL 
    dpdtBodyIntercept  = 0.0_CGREAL 
        
    xImagePoint        = 0.0_CGREAL
    yImagePoint        = 0.0_CGREAL
    zImagePoint        = 0.0_CGREAL
    probeLength        = 0.0_CGREAL

    closestMarker      = 0
    closestMarkerRatio = 0.0_CGREAL
    closestElementGC   = 0.0_CGREAL

    res_laplacian = 0.0_CGREAL
    alphaGhost    = 0.0_CGREAL

    ghostNScoeff(:,:)  = 0.0_CGREAL
    ghostNScoeff(1,:)  = 1.0_CGREAL

   END SUBROUTINE GCM_AllocateGhostCellArrays
!------------------------------------------------------------------------
   SUBROUTINE GCM_AllocateDeadCellArrays()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    
    IMPLICIT NONE

! Allocate Memory for various arrays pertinent to GCM Dead Cells

    ALLOCATE(iDead(nDead))
    ALLOCATE(jDead(nDead))
    ALLOCATE(kDead(nDead))

    ALLOCATE(closestMarkerDead(nDead))
    ALLOCATE(closestElementDead(nDead))
    ALLOCATE(closestMarkerRatioDead(nDead))

    ALLOCATE(xBodyInterceptDead(nDead))
    ALLOCATE(yBodyInterceptDead(nDead))
    ALLOCATE(zBodyInterceptDead(nDead))

    ALLOCATE(iDeadCellIndex(nDead))
    ALLOCATE(jDeadCellIndex(nDead))
    ALLOCATE(kDeadCellIndex(nDead))

    ALLOCATE(coeffGCMDeadN(iRowMax,nDead))
    ALLOCATE(coeffGCMDeadD(iRowMax,nDead))

    ALLOCATE(uBodyInterceptDead(nDead))
    ALLOCATE(vBodyInterceptDead(nDead))
    ALLOCATE(wBodyInterceptDead(nDead))

    ALLOCATE(alphaDead(nDead))

    ALLOCATE(dpdnInterceptDead(nDead))

    ALLOCATE(xNormInterceptDead(nDead))
    ALLOCATE(yNormInterceptDead(nDead))
    ALLOCATE(zNormInterceptDead(nDead))

    ALLOCATE(probeLengthDead   (nDead))

    !   Initialize arrays

    iDead                   = 0
    jDead                   = 0
    kDead                   = 0

    closestMarkerDead       = 0
    closestElementDead      = 0
    closestMarkerRatioDead  = 0.0_CGREAL

    xBodyInterceptDead      = 0.0_CGREAL
    yBodyInterceptDead      = 0.0_CGREAL
    zBodyInterceptDead      = 0.0_CGREAL

    iDeadCellIndex          = 0
    jDeadCellIndex          = 0
    kDeadCellIndex          = 0

    coeffGCMDeadN           = 0.0_CGREAL
    coeffGCMDeadD           = 0.0_CGREAL

    uBodyInterceptDead      = 0.0_CGREAL
    vBodyInterceptDead      = 0.0_CGREAL
    wBodyInterceptDead      = 0.0_CGREAL

    xNormInterceptDead      = 0.0_CGREAL
    yNormInterceptDead      = 0.0_CGREAL
    zNormInterceptDead      = 0.0_CGREAL
    dpdnInterceptDead       = 0.0_CGREAL

   END SUBROUTINE GCM_AllocateDeadCellArrays

!------------------------------------------------------------------------
   SUBROUTINE MARKER_allocate_memory()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    USE usr_module
    USE fsi_module

    IMPLICIT NONE

    INTEGER :: nBodyPtsMax,iBody

    ! Marker point arrays
    nPtsMax = MAXVAL(nPtsBodyMarker(:))

    ALLOCATE(xBodyMarker(nPtsMax,nBody))
    ALLOCATE(yBodyMarker(nPtsMax,nBody))
    ALLOCATE(zBodyMarker(nPtsMax,nBody))
    ALLOCATE(Flag_outside_Marker(nPtsMax,nBody)) 
    !ALLOCATE(distMarker(nPtsMax))
    !ALLOCATE(NeighElemInd(NSIZE))

    ALLOCATE(cnt_fg12(nPtsMax,nBody))
    ALLOCATE(cnt_fg21(nPtsMax,nBody))
    ALLOCATE(cnt_fg13(nPtsMax,nBody))
    ALLOCATE(cnt_fg31(nPtsMax,nBody))
    ALLOCATE(cnt_fg23(nPtsMax,nBody))
    ALLOCATE(cnt_fg32(nPtsMax,nBody))

    ALLOCATE(cnt_fx(nPtsMax,nBody))
    ALLOCATE(cnt_fy(nPtsMax,nBody))
    ALLOCATE(cnt_fz(nPtsMax,nBody))
    ALLOCATE(cnt_force(nPtsMax,nBody))
    
    ALLOCATE(cnt_d12(nPtsMax,nBody))
    ALLOCATE(cnt_d21(nPtsMax,nBody))
    ALLOCATE(cnt_d13(nPtsMax,nBody))
    ALLOCATE(cnt_d31(nPtsMax,nBody))
    ALLOCATE(cnt_d23(nPtsMax,nBody))
    ALLOCATE(cnt_d32(nPtsMax,nBody))

    ALLOCATE(cnt_EM12(nPtsMax,nBody))
    ALLOCATE(cnt_EM21(nPtsMax,nBody))
    ALLOCATE(cnt_EM13(nPtsMax,nBody))
    ALLOCATE(cnt_EM31(nPtsMax,nBody))
    ALLOCATE(cnt_EM23(nPtsMax,nBody))
    ALLOCATE(cnt_EM32(nPtsMax,nBody))

    ALLOCATE(dotP(nPtsMax,nBody))
    ALLOCATE(dotP2(nPtsMax,nBody))
    !ALLOCATE(v_num(10590))
    ALLOCATE(set_old(nPtsMax,nBody))

    xBodyMarker     = 0.0_CGREAL
    yBodyMarker     = 0.0_CGREAL
    zBodyMarker     = 0.0_CGREAL
    !distMarker      = 0.0_CGREAL

    cnt_fg12 = 0
    cnt_fg21 = 0
    cnt_fg13 = 0
    cnt_fg31 = 0
    cnt_fg23 = 0
    cnt_fg32 = 0

    cnt_fx = 0.0_CGREAL
    cnt_fy = 0.0_CGREAL
    cnt_fz = 0.0_CGREAL 
    cnt_force = 0.0_CGREAL
 
    dotP = 0
    dotP2 = 0
    !v_num = -1
    set_old = 0
  
!!$    !!! Added by Paulo Ferreira de Sousa on 02/27/2009
!!$    ALLOCATE(xBMPosition(1:nPtsMax,nBody))
!!$    ALLOCATE(yBMPosition(1:nPtsMax,nBody))
!!$    ALLOCATE(zBMPosition(1:nPtsMax,nBody))
!!$    xBMPosition     = 0.0_CGREAL
!!$    yBMPosition     = 0.0_CGREAL
!!$    zBMPosition     = 0.0_CGREAL

    ALLOCATE(xBodyMarkerOld(nPtsMax,nBody))
    ALLOCATE(yBodyMarkerOld(nPtsMax,nBody))
    ALLOCATE(zBodyMarkerOld(nPtsMax,nBody))    
    xBodyMarkerOld  = 0.0_CGREAL
    yBodyMarkerOld  = 0.0_CGREAL
    zBodyMarkerOld  = 0.0_CGREAL

    ALLOCATE(sBodyMarker(nPtsMax,nBody))
    ALLOCATE(dsBodyMarker(nPtsMax,nBody))
    ALLOCATE(xNormBodyMarker(nPtsMax,nBody))
    ALLOCATE(yNormBodyMarker(nPtsMax,nBody))
    ALLOCATE(zNormBodyMarker(nPtsMax,nBody))
    sBodyMarker     = 0.0_CGREAL
    dsBodyMarker    = 0.0_CGREAL
    xNormBodyMarker = 0.0_CGREAL
    yNormBodyMarker = 0.0_CGREAL
    zNormBodyMarker = 0.0_CGREAL

    ALLOCATE(uBodyMarker(nPtsMax,nBody))
    ALLOCATE(vBodyMarker(nPtsMax,nBody)) 
    ALLOCATE(wBodyMarker(nPtsMax,nBody)) 
    uBodyMarker = 0.0_CGREAL
    vBodyMarker = 0.0_CGREAL
    wBodyMarker = 0.0_CGREAL

    !Ye, store the true velocity from solid before filtered
    ALLOCATE(u1BodyMarker(nPtsMax,nBody))
    ALLOCATE(v1BodyMarker(nPtsMax,nBody))
    ALLOCATE(w1BodyMarker(nPtsMax,nBody))
    u1BodyMarker = 0.0_CGREAL
    v1BodyMarker = 0.0_CGREAL
    w1BodyMarker = 0.0_CGREAL

    ALLOCATE(axBodyMarker(nPtsMax,nBody))
    ALLOCATE(ayBodyMarker(nPtsMax,nBody)) 
    ALLOCATE(azBodyMarker(nPtsMax,nBody)) 
    axBodyMarker = 0.0_CGREAL
    ayBodyMarker = 0.0_CGREAL
    azBodyMarker = 0.0_CGREAL

    ALLOCATE(uBodyMarkerOld(nPtsMax,nBody))
    ALLOCATE(vBodyMarkerOld(nPtsMax,nBody))
    ALLOCATE(wBodyMarkerOld(nPtsMax,nBody))
    uBodyMarkerOld = 0.0_CGREAL
    vBodyMarkerOld = 0.0_CGREAL
    wBodyMarkerOld = 0.0_CGREAL

    ALLOCATE(uBodyMarkerTmp(nPtsMax,nBody))
    ALLOCATE(vBodyMarkerTmp(nPtsMax,nBody))
    ALLOCATE(wBodyMarkerTmp(nPtsMax,nBody))
    uBodyMarkerTmp = 0.0_CGREAL
    vBodyMarkerTmp = 0.0_CGREAL
    wBodyMarkerTmp = 0.0_CGREAL

    ALLOCATE(xMarkerShear(nPtsMax,nBody))
    ALLOCATE(yMarkerShear(nPtsMax,nBody))
    ALLOCATE(zMarkerShear(nPtsMax,nBody))

    ALLOCATE(xMarkerForce(nPtsMax,nBody))
    ALLOCATE(yMarkerForce(nPtsMax,nBody))
    ALLOCATE(zMarkerForce(nPtsMax,nBody))

    !Ye, store the unfiltered BM force
    ALLOCATE(xMarkerForce0(nPtsMax,nBody))
    ALLOCATE(yMarkerForce0(nPtsMax,nBody))
    ALLOCATE(zMarkerForce0(nPtsMax,nBody))

    ALLOCATE(xMarkerForceOld(nPtsMax,nBody))
    ALLOCATE(yMarkerForceOld(nPtsMax,nBody))
    ALLOCATE(zMarkerForceOld(nPtsMax,nBody))

    ALLOCATE(xMarkerStress(nPtsMax,nBody))
    ALLOCATE(yMarkerStress(nPtsMax,nBody))
    ALLOCATE(zMarkerStress(nPtsMax,nBody))
    ALLOCATE(MarkerArea  (nPtsMax,nBody))
    ALLOCATE(xyzMarkerForceTmp  (nPtsMax*6))
    ALLOCATE(xyzMarkerForceTmp_i(nPtsMax*6))

    xMarkerShear = 0.0_CGREAL
    yMarkerShear = 0.0_CGREAL
    zMarkerShear = 0.0_CGREAL

    xMarkerForce = 0.0_CGREAL
    yMarkerForce = 0.0_CGREAL
    zMarkerForce = 0.0_CGREAL

    xMarkerForce0 = 0.0_CGREAL
    yMarkerForce0 = 0.0_CGREAL
    zMarkerForce0 = 0.0_CGREAL

    xMarkerStress = 0.0_CGREAL
    yMarkerStress = 0.0_CGREAL
    zMarkerStress = 0.0_CGREAL
    MarkerArea         = 0.0_CGREAL
    xyzMarkerForceTmp  = 0.0_CGREAL
    xyzMarkerForceTmp_i= 0.0_CGREAL

    ALLOCATE(dpdnBodyMarker(nPtsMax,nBody))
    dpdnBodyMarker = 0.0_CGREAL

    !ALLOCATE(xTang1BodyMarker(nPtsMax,nBody))
    !ALLOCATE(yTang1BodyMarker(nPtsMax,nBody))
    !ALLOCATE(zTang1BodyMarker(nPtsMax,nBody))

    !xTang1BodyMarker = 0.0_CGREAL
    !yTang1BodyMarker = 0.0_CGREAL
    !zTang1BodyMarker = 0.0_CGREAL

    !ALLOCATE(xTang2BodyMarker(nPtsMax,nBody))
    !ALLOCATE(yTang2BodyMarker(nPtsMax,nBody))
    !ALLOCATE(zTang2BodyMarker(nPtsMax,nBody))

    !xTang2BodyMarker = 0.0_CGREAL
    !yTang2BodyMarker = 0.0_CGREAL
    !zTang2BodyMarker = 0.0_CGREAL

    ALLOCATE(cxt(nBody))  
    ALLOCATE(cyt(nBody))  
    ALLOCATE(czt(nBody))  

    ALLOCATE(cxs(nBody))  
    ALLOCATE(cys(nBody))  
    ALLOCATE(czs(nBody))  

    ALLOCATE(cmxt(nBody)) 
    ALLOCATE(cmyt(nBody)) 
    ALLOCATE(cmzt(nBody)) 

    ALLOCATE(cmxs(nBody)) 
    ALLOCATE(cmys(nBody)) 
    ALLOCATE(cmzs(nBody)) 

  END SUBROUTINE MARKER_allocate_memory

!-------------------------------------------------------------------
! Allocate Memory for various arrays pertinent to unstructured surface

   SUBROUTINE UNSTRUC_allocate_memory()
    
    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    
    IMPLICIT NONE

    INTEGER :: nBodyPtsMax,nTriElemMax,iBody
    LOGICAL :: unstruc


    nTriElemMax = MAXVAL(totNumTriElem(:))
    ALLOCATE(triElemNeig(3,nTriElemMax,nBody))

    ALLOCATE(triElemtang1X(nTriElemMax,nBody))
    ALLOCATE(triElemtang1Y(nTriElemMax,nBody))
    ALLOCATE(triElemtang1Z(nTriElemMax,nBody))

    ALLOCATE(triElemtang2X(nTriElemMax,nBody))
    ALLOCATE(triElemtang2Y(nTriElemMax,nBody))
    ALLOCATE(triElemtang2Z(nTriElemMax,nBody))
    
    ALLOCATE(triElemNormX(nTriElemMax,nBody))
    ALLOCATE(triElemNormY(nTriElemMax,nBody))
    ALLOCATE(triElemNormZ(nTriElemMax,nBody))
    
    ALLOCATE(triElemCentX(nTriElemMax,nBody))
    ALLOCATE(triElemCentY(nTriElemMax,nBody))
    ALLOCATE(triElemCentZ(nTriElemMax,nBody))
    
    ALLOCATE(triElemArea(nTriElemMax,nBody))
    
    ALLOCATE(pointOutsideBodyX(nBody))
    ALLOCATE(pointOutsideBodyY(nBody))
    ALLOCATE(pointOutsideBodyZ(nBody))
    ALLOCATE(surfArea(nBody))

   END SUBROUTINE UNSTRUC_allocate_memory

!------------------------------------------------------------------------
! static arrays for GCM that need to be declared only once

   SUBROUTINE GCM_allocate_static_arrays()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays

    IMPLICIT NONE

    ALLOCATE(incI(iRowMax))
    ALLOCATE(incJ(iRowMax))
    ALLOCATE(incK(iRowMax))
    ALLOCATE(iPvt(iRowMax))
    ALLOCATE(work(iRowMax))
    ALLOCATE(vanMatrixD(iRowMax,iRowMax))
    ALLOCATE(vanMatrixN(iRowMax,iRowMax))

    ALLOCATE(coeffGCMDirc(iRowMax))
    ALLOCATE(coeffGCMNeum(iRowMax))

    ALLOCATE(intI(6))
    ALLOCATE(intJ(6))
    ALLOCATE(intK(6))

    ! These allow us to define stencil image point
    ! Assumed clockwise from lower left corner.

      incI(1) = 0
      incJ(1) = 0
      incK(1) = 0

      incI(2) = 0
      incJ(2) = 1
      incK(2) = 0

      incI(3) = 1
      incJ(3) = 1
      incK(3) = 0

      incI(4) = 1
      incJ(4) = 0
      incK(4) = 0

      incI(5) = 0 
      incJ(5) = 0
      incK(5) = 1

      incI(6) = 0
      incJ(6) = 1
      incK(6) = 1

      incI(7) = 1
      incJ(7) = 1
      incK(7) = 1

      incI(8) = 1
      incJ(8) = 0
      incK(8) = 1 

      !======
      intI(1) = 1
      intJ(1) = 0
      intK(1) = 0

      intI(2) = -1
      intJ(2) = 0
      intK(2) = 0

      intI(3) = 0
      intJ(3) = 1
      intK(3) = 0

      intI(4) = 0
      intJ(4) = -1
      intK(4) = 0

      intI(5) = 0
      intJ(5) = 0
      intK(5) = 1

      intI(6) = 0
      intJ(6) = 0
      intK(6) = -1
      !======

    END SUBROUTINE GCM_allocate_static_arrays
!------------------------------------------------------------------------
   SUBROUTINE GCM_DeallocateGhostCellArrays()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays

    IMPLICIT NONE

! Deallocate Memory for various arrays pertinent to GCM Body Markers

    DEALLOCATE(iGhost)
    DEALLOCATE(jGhost) 
    DEALLOCATE(kGhost) 

    DEALLOCATE(coeffGCMD)
    DEALLOCATE(coeffGCMN)
    
    DEALLOCATE(xBodyInterceptTang)
    DEALLOCATE(yBodyInterceptTang)
    DEALLOCATE(zBodyInterceptTang)

    DEALLOCATE(xBodyInterceptNorm)
    DEALLOCATE(yBodyInterceptNorm)
    DEALLOCATE(zBodyInterceptNorm)

    DEALLOCATE(xBodyIntercept)
    DEALLOCATE(yBodyIntercept)
    DEALLOCATE(zBodyIntercept)

    DEALLOCATE(uBodyIntercept)
    DEALLOCATE(vBodyIntercept) 
    DEALLOCATE(wBodyIntercept) 
    DEALLOCATE(pBodyIntercept) 
    DEALLOCATE(dpdnBodyIntercept) 
    DEALLOCATE(dpdtBodyIntercept) 

    DEALLOCATE(closestMarker) 
    DEALLOCATE(closestMarkerRatio) 
    DEALLOCATE(closestElementGC)
          
    DEALLOCATE(iCellIndex)
    DEALLOCATE(jCellIndex)
    DEALLOCATE(kCellIndex)

    DEALLOCATE(xImagePoint)
    DEALLOCATE(yImagePoint) 
    DEALLOCATE(zImagePoint) 
    DEALLOCATE(probeLength)

    DEALLOCATE(res_laplacian)
    DEALLOCATE(alphaGhost)

    DEALLOCATE(ghostNScoeff)

   END SUBROUTINE GCM_DeallocateGhostCellArrays
!------------------------------------------------------------------------
   SUBROUTINE GCM_DeallocateDeadCellArrays()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays
    USE GCM_arrays
    
    IMPLICIT NONE

! Deallocate Memory for various arrays pertinent to GCM Dead Cells

    DEALLOCATE(iDead)
    DEALLOCATE(jDead)
    DEALLOCATE(kDead)

    DEALLOCATE(closestMarkerDead)
    DEALLOCATE(closestElementDead)
    DEALLOCATE(closestMarkerRatioDead)

    DEALLOCATE(xBodyInterceptDead)
    DEALLOCATE(yBodyInterceptDead)
    DEALLOCATE(zBodyInterceptDead)

    DEALLOCATE(iDeadCellIndex)
    DEALLOCATE(jDeadCellIndex)
    DEALLOCATE(kDeadCellIndex)

    DEALLOCATE(coeffGCMDeadN)
    DEALLOCATE(coeffGCMDeadD)

    DEALLOCATE(uBodyInterceptDead)
    DEALLOCATE(vBodyInterceptDead)
    DEALLOCATE(wBodyInterceptDead)

    DEALLOCATE(xNormInterceptDead)
    DEALLOCATE(yNormInterceptDead)
    DEALLOCATE(zNormInterceptDead)
    DEALLOCATE(dpdnInterceptDead)

    DEALLOCATE(alphaDead)
    DEALLOCATE(probeLengthDead)

   END SUBROUTINE GCM_DeallocateDeadCellArrays
!------------------------------------------------------------------------
