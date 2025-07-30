!******************************************************************************
!
! Purpose: generalized kernel to compute the coordinates of the intercept 
!          points for any generic point onto the unstructured surface mesh
!
! Description: none.
!
! Input: iGP, jGP, kGP = indices of Generic Point
!
! Output: xBI, yBI , zBI         = coordinates of Body Intercept,
!         closestElementGP       = closest Element for Generic point,
!         xBITang1GP, yBITang1GP = coordinates of tangent vector 1
!         xBITang2GP, yBITang2GP = coordinates of tangent vector 2
!         xBINormGP, yBINormGP   = coordinates of normal vector
!
! Notes: none.
!
!
!******************************************************************************

      SUBROUTINE calc_bodyIntercept_Unstruc_cnt(iPt,fg,L1,L2,
     &           cnt_d,ElemNum)

       USE contact

       IMPLICIT NONE

!... parameters
       INTEGER, INTENT(IN)  :: iPt, fg, L1, L2
       REAL*8               :: xGP, yGP, zGP
       REAL*8               :: xBI, yBI, zBI
       REAL*8               :: cnt_d
       INTEGER  ::  ElemNum
    
!... loop variables

       INTEGER :: iEdge,m,n,nc,ii

!... local variables
 
       INTEGER, PARAMETER       :: NSIZE = 10
       INTEGER, PARAMETER       :: MSIZE = 10
       INTEGER   :: NeighElemInd(1:NSIZE)
       REAL*8    :: distMarker(1:nPtsBM)

      INTEGER                :: nCheck
      INTEGER,DIMENSION(1)   :: iDummy(1)
      INTEGER                :: numNeighElement
      INTEGER                :: elemInd,node1,node2,node3,nMarker,iErr
      INTEGER                :: nEdges,nodeCV,nodeA,nodeB
      INTEGER           :: cElementGP(1:MSIZE),closestNodeGP(1:MSIZE)
      INTEGER                  :: shortestProbe
      INTEGER                  :: closestElementGP

      REAL*8        :: cNormal,dMin,dsIntercept,xM,yM,zM
      REAL*8        :: area123,areaDiff,distBIElem,distBIElemMin
      REAL*8        :: epsiArea,distGPEI,distGPEIMin
      REAL*8        :: xBITemp, yBITemp, zBITemp
      REAL*8        :: xCV, yCV, zCV
      REAL*8        :: xEI, yEI, zEI
      REAL*8        :: xEITemp, yEITemp, zEITemp
      REAL*8        :: vec01x, vec01y, vec01z
      REAL*8        :: vec12x, vec12y, vec12z
      REAL*8        :: magnitude12,magnitude12Inv,projectedLength
      REAL*8 :: xBIT(1:MSIZE),yBIT(1:MSIZE),zBIT(1:MSIZE),dist(1:MSIZE)
      REAL*8        :: xNorm,yNorm,zNorm,xCG,yCG,zCG
      REAL*8        :: planeConst,distanceToPlane
      REAL*8, DIMENSION(3) :: xVert,yVert,zVert

      INTEGER                :: i, is_in, is_in2

      REAL*8        :: NormMag
      REAL*8        :: NormXtemp0,NormYtemp0,NormZtemp0
      REAL*8        :: NormXtemp1,NormYtemp1,NormZtemp1
      REAL*8        :: NormXtemp2,NormYtemp2,NormZtemp2
      REAL*8        :: normX(1:MSIZE),normY(1:MSIZE),normZ(1:MSIZE)
      REAL*8        :: rto

!******************************************************************************
! NCheck:  Number of closesest nodes to check
! high values of this variable increases robustness of procedure 
! and also CPU time for finding body intercept.

       !Ye,debug
       !if (iPt.eq.45) then
       !   write(560,*)L1,L2
       !endif

       BI_normX(iPt) = 0.0d0
       BI_normY(iPt) = 0.0d0
       BI_normZ(iPt) = 0.0d0

       NormXtemp0 = 0.0d0
       NormYtemp0 = 0.0d0
       NormZtemp0 = 0.0d0
     
       NormXtemp1 = 0.0d0
       NormYtemp1 = 0.0d0
       NormZtemp1 = 0.0d0

       NormXtemp2 = 0.0d0
       NormYtemp2 = 0.0d0
       NormZtemp2 = 0.0d0

       NormMag = 0.0d0

       normX = 0.0d0
       normY = 0.0d0
       normZ = 0.0d0

       xGP = BM_x(iPt)
       yGP = BM_y(iPt)
       zGP = BM_z(iPt)

       nCheck = 3 

! ============================================================================
! Get closestMarker for generic point 
! ============================================================================

      dMin = 1.0E+5

      !Ye,debug
      !if (iPt.eq.45)then
      !   write(559,*)xGP,yGP,zGP
      !endif
      !=======

      DO m = 1, nPtsBM
       if( (m.lt.L1).OR.(m.gt.L2).OR.(mark(m).eq.0) )THEN
         distMarker(m)=1.0E+20
       else
         distMarker(m)=(BM_x(m)-xGP)**2+(BM_y(m)-yGP)**2
     &                +(BM_z(m)-zGP)**2
       endif
      ENDDO!m

      !Hsu: abs(x1-x2)<c,default cnt_d, causing 0 cnt force
      if (minval(distmarker(1:nPtsBM)).ge.cnt_c**2) then
         closestElementGP = 1
         goto 918
      endif

      DO nc = 1,NCheck
         iDummy  = MINLOC(distmarker(1:nPtsBM))
         closestNodeGP(nc)             = iDummy(1)
         distmarker(closestNodeGP(nc)) = 1.0E+20
      ENDDO

! ============================================================================
! Find elements that share closest node/marker
! ============================================================================
      is_in = 0     

      DO nc = 1,nCheck

      is_in2 = 0

      numNeighElement = 0
!      DO m=1,totNumElem
!         IF ( (ElemNeig(1,m) == closestNodeGP(nc) .or.
!     &         ElemNeig(2,m) == closestNodeGP(nc) .or.
!     &         ElemNeig(3,m) == closestNodeGP(nc)) .AND.
!     &        (mark(ElemNeig(1,m)).ne.0 .and.
!     &         mark(ElemNeig(2,m)).ne.0 .and.
!     &         mark(ElemNeig(3,m)).ne.0) )THEN
!              numNeighElement               = numNeighElement + 1
!              NeighElemInd(numNeighElement) = m
!         ENDIF
!      ENDDO

       do m=1,8
            ii=ElemC(closestNodeGP(nc),m)
            IF (ii.eq.0) CYCLE
            IF ( (mark(ElemNeig(1,ii)).eq.0).or.
     &           (mark(ElemNeig(2,ii)).eq.0).or.
     &           (mark(ElemNeig(3,ii)).eq.0) ) CYCLE
            numNeighElement = numNeighElement + 1
            NeighElemInd(numNeighElement) = ii
         enddo

         IF (numNeighElement.eq.0) THEN
            WRITE(*,*)'ERROR!'
         ENDIF
! ============================================================================
!   Trap error if array NeighElemenInd overflows
! ============================================================================

      IF ( numNeighElement > NSIZE ) THEN
         WRITE(*,*)'CNT_calc_bodyIntercept_Unstruc: ,
     &             Memory Overflow Error for NeighElemInd'
        WRITE(*,*) ' Allocated size = ',NSIZE
        WRITE(*,*) ' Current size   = ',numNeighElement
        WRITE(*,*) ' Aborting Run'

        STOP
      ENDIF ! NeighElemInd

! ============================================================================
!   Determine which element contains normal intercept
! ============================================================================

      closestElementGP = 0

      epsiArea = 1.0E-04
      distBIElemMin = 1.0E+16

      DO n = 1,numNeighElement

         elemInd = NeighElemInd(n)

         node1   = ElemNeig(1,elemInd)
         node2   = ElemNeig(2,elemInd)
         node3   = ElemNeig(3,elemInd)

         IF ( ((mark(node1)).eq.0).AND.((mark(node2)).eq.0).AND.
     &        ((mark(node3)).eq.0) ) CYCLE

         xVert(1)  = BM_x(node1)
         yVert(1)  = BM_y(node1)
         zVert(1)  = BM_z(node1)

         xVert(2)  = BM_x(node2)
         yVert(2)  = BM_y(node2)
         zVert(2)  = BM_z(node2)

         xVert(3)  = BM_x(node3)
         yVert(3)  = BM_y(node3)
         zVert(3)  = BM_z(node3)

         xNorm     = ElemNorm_x(elemInd)
         yNorm     = ElemNorm_y(elemInd)
         zNorm     = ElemNorm_z(elemInd)

! ******************************************************************************
! equation of plane (note our normals are unit normals)
     
!  n  x + n  y + n  z + planeConst = 0
!   x         y      z
! ******************************************************************************

         planeConst=- xNorm*xVert(1)
     &              - yNorm*yVert(1)
     &              - zNorm*zVert(1)

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

          distanceToPlane = -(  xNorm*xGP
     &                        + yNorm*yGP    
     &                        + zNorm*zGP   
     &                        + planeConst )  

         xBITemp = xGP + xNorm*distanceToPlane
         yBITemp = yGP + yNorm*distanceToPlane
         zBITemp = zGP + zNorm*distanceToPlane
                
! ******************************************************************************
!    Check if BI inside the triangle of the surface element
!     through area differences
! ******************************************************************************

        CALL check_BIInside_cnt(xVert,yVert,zVert,xNorm,yNorm,zNorm,  
     &                     xBITemp,yBITemp,zBITemp,area123,areaDiff)

!DEBUG
!     IF (iGP == 38 .and. jGP == 33 .and. kGP == 19) then
!        WRITE(*,*)elemInd
!        WRITE(*,'(6F12.5)')xGP,yGP,zGP
!        WRITE(*,'(6F12.5)')xVert(1),yVert(1),zVert(1)
!        WRITE(*,'(6F12.5)')xVert(2),yVert(2),zVert(2)
!        WRITE(*,'(6F12.5)')xVert(3),yVert(3),zVert(3)
!        WRITE(*,'(6F12.5)')xNorm,yNorm,zNorm,PlaneConst
!        WRITE(*,'(6F12.5)')xBITemp,yBITemp,zBITemp,areadiff
!        pause
!     ENDIF

!--------------------------------------------------------------------------
!    Select closest Elem and BI coordinates:
!     If BI falls inside the element use that
!     Else Base the selection on the minimum distance
!       between BI and either the norm to closest side or vertices of side
!--------------------------------------------------------------------------
       !Ye,debug
       !if ((iPt.eq.2549).and.(fg.eq.12))then
       !   write(554,*)iPt,nc,closestNodeGP(nc),elemInd
       !   write(554,*)xBITemp,yBITemp,zBITemp
       !   write(554,*)ABS(areaDiff),epsiArea*area123   
       !   write(554,*)'----'
       !endif

       IF ( ABS(areaDiff) < epsiArea*area123 ) THEN
          xBI = xBITemp
          yBI = yBITemp
          zBI = zBITemp
          closestElementGP = elemInd
          is_in  = 1
          is_in2 = 1
          GOTO 999
       ELSE
          xCG       = ElemCentx(elemInd)
          yCG       = ElemCenty(elemInd)
          zCG       = ElemCentz(elemInd)
        
         CALL calc_BIOutside_cnt(xVert,yVert,zVert,xCG,yCG,zCG,    
     &                           xBITemp,yBITemp,zBITemp,distBIElem)
        
         IF (distBIElem <= distBIElemMin) THEN
            distBIElemMin = distBIElem
            closestElementGP = elemInd
            xBI = xBITemp
            yBI = yBITemp
            zBI = zBITemp
          ENDIF ! distBIElem
       ENDIF ! areaDiff

!DEBUG
!    IF (iGP == 49 .and. jGP == 2 .and. kGP == 2) then
!       WRITE(*,*)elemInd
!       WRITE(*,*)distBIElemMin
!       WRITE(*,'(6F12.5)')xVert(1),yVert(1),zVert(1)
!       WRITE(*,'(6F12.5)')xVert(2),yVert(2),zVert(2)
!       WRITE(*,'(6F12.5)')xVert(3),yVert(3),zVert(3)
!       WRITE(*,'(6F12.5)')xBItemp,yBItemp,zBItemp
!       WRITE(*,'(1F12.5,1I6)')areadiff,closestElementGP
!    ENDIF
!DEBUG

       ENDDO ! n

! ============================================================================
!   Compute coordinates of Body Intercept in a robust manner
!    for the case where the temporary BI is located outside 
!    all the surface elements
!    1. Load coordinates of closest vertex (CV)
! ============================================================================

         xCV = BM_x(closestNodeGP(nc))
         yCV = BM_y(closestNodeGP(nc))
         zCV = BM_z(closestNodeGP(nc))

         distGPEIMin = 1.0E+16

! ============================================================================
!    2. Detemine the indices of the 2 vertices connected to CV
! ============================================================================

         node1   = ElemNeig(1,closestElementGP)
         node2   = ElemNeig(2,closestElementGP)
         node3   = ElemNeig(3,closestElementGP)

         IF ( node1 == closestNodeGP(nc) ) THEN
            nodeCV = node1
            nodeA  = node2
            nodeB  = node3
         ELSEIF ( node2 == closestNodeGP(nc) ) THEN
            nodeCV = node2
            nodeA  = node3
            nodeB  = node1
         ELSEIF ( node3 == closestNodeGP(nc) ) THEN
            nodeCV = node3
            nodeA  = node1
            nodeB  = node2
         END IF ! node1 

       !Ye,debug
!       if ((iPt.eq.5003).and.(fg.eq.21))then
!          write(700,*)'nc:',nc
!          write(700,*)closestNodeGP(nc),closestElementGP
!          write(700,*)nodeCV, nodeA, nodeB
!          write(700,*)ElemNorm_x(closestElementGP),
!     &                ElemNorm_y(closestElementGP),
!     &                ElemNorm_z(closestElementGP)
!          write(700,*)'---'
!       endif

! ============================================================================
!   3. Compute edge01 (CV-->GP), edge12 (CV-->A), edge13 (CV-->B) vectors
!      Project vector GP-CV onto CV-A or CV-B to find temporary edge intercept
!      If the projectedLength is < 0 or > edgeLength, EI is outside
!      Else EI is inside then compute its Location and distance to GP
! ============================================================================

       vec01x = xGP - xCV
       vec01y = yGP - yCV
       vec01z = zGP - zCV

       nEdges = 2

       DO iEdge = 1, nEdges
        SELECT CASE(iEdge)
        CASE(1)
         node1 = nodeA
        CASE(2)
         node1 = nodeB
        END SELECT ! iEdge

       vec12x = BM_x(node1) - xCV
       vec12y = BM_y(node1) - yCV
       vec12z = BM_z(node1) - zCV

       magnitude12 = SQRT(vec12x**2 + vec12y**2 + vec12z**2)

       magnitude12Inv = 1.0d0/magnitude12
     
       vec12x = vec12x*magnitude12Inv
       vec12y = vec12y*magnitude12Inv
       vec12z = vec12z*magnitude12Inv

       projectedLength = vec01x*vec12x +vec01y*vec12y +vec01z*vec12z 
 
       !Ye,debug
!       if ((iPt.eq.5003).and.(fg.eq.21))then
!          write(701,*)'nc:',nc
!          write(701,*)'Edge:',iEdge
!          write(701,*)closestNodeGP(nc),closestElementGP
!          write(701,*)nodeCV, nodeA, nodeB
!          write(701,*)vec12x,vec12y,vec12z
!          write(701,*)projectedLength,magnitude12
!          write(701,*)'---'
!       endif

!----------------------------------------------------------------------------
!     Edge-Intercept (EI) point is outside Edge if pL < 0 or pL>magnitude12
!      else EI is inside and compute coordinates and distance
!      No need to take SQRT for distGPEI to save computations
!      Load EI into temporary value
!----------------------------------------------------------------------------

      IF ( projectedLength < 0.0d0 .OR.
     &     projectedLength > magnitude12     ) THEN
        xEITemp = xCV
        yEITemp = yCV
        zEITemp = zCV
        NormXtemp1 = xNormBM(nodeCV)
        NormYtemp1 = yNormBM(nodeCV)
        NormZtemp1 = zNormBM(nodeCV)
      ELSE
        xEITemp = xCV + projectedLength*vec12x
        yEITemp = yCV + projectedLength*vec12y
        zEITemp = zCV + projectedLength*vec12z
        !Simple Average=====
        !NormXtemp0=0.5*(xNormBM(nodeCV)+xNormBM(node1))
        !NormYtemp0=0.5*(yNormBM(nodeCV)+yNormBM(node1))
        !NormZtemp0=0.5*(zNormBM(nodeCV)+zNormBM(node1))
        !Dist-based Average=====
        rto = projectedLength/magnitude12
        NormXtemp0=rto*xNormBM(nodeCV)+(1.0-rto)*xNormBM(node1)
        NormYtemp0=rto*yNormBM(nodeCV)+(1.0-rto)*yNormBM(node1)
        NormZtemp0=rto*zNormBM(nodeCV)+(1.0-rto)*zNormBM(node1)

        NormMag = SQRT(NormXtemp0**2+NormYtemp0**2+NormZtemp0**2)
        if (NormMag.eq.0.0d0)then
           write(*,*)'Error! NormMag = 0.0'
        endif
        NormXtemp1 = NormXtemp0/NormMag
        NormYtemp1 = NormYtemp0/NormMag
        NormZtemp1 = NormZtemp0/NormMag
      END IF ! projectedLength

       !Ye,debug
!       if ((iPt.eq.5003).and.(fg.eq.21))then
!             write(702,*)'nc:',nc
!             write(702,*)'Edge:',iEdge,node1
!             write(702,*)rto,projectedLength,magnitude12
!             write(702,*)nodeCV,node1
!             write(702,*)xNormBM(nodeCV),
!     &                   yNormBM(nodeCV),
!     &                   zNormBM(nodeCV)
!             write(702,*)xNormBM(node1),yNormBM(node1),
!     &                   zNormBM(node1)
!             write(702,*)NormXtemp0,NormYtemp0,NormZtemp0
!             write(702,*)NormXtemp1,NormYtemp1,NormZtemp1
!             write(702,*)'---'
!       endif
!----------------------------------------------------------------------------
!     Find mininum value of |GP-EI| and corresponding EI
!----------------------------------------------------------------------------

      distGPEI = (xEITemp-xGP)**2 
     &         + (yEITemp-yGP)**2 
     &         + (zEITemp-zGP)**2 
  
      IF ( distGPEI < distGPEIMin ) THEN
        distGPEIMin = distGPEI
        xEI = xEITemp
        yEI = yEITemp
        zEI = zEITemp
        NormXtemp2 = NormXtemp1
        NormYtemp2 = NormYtemp1
        NormZtemp2 = NormZtemp1
      END IF ! distGPEI

      END DO ! iEdge
 
      xBI = xEI 
      yBI = yEI
      zBI = zEI

      !Ye,debug
!       if ((iPt.eq.5003).and.(fg.eq.21))then
!          if (nc.eq.1)then
!             write(703,*)NormXtemp2,NormYtemp2,
!     &                   NormZtemp2
!             write(703,*)'---'
!          endif
!       endif

!DEBUG
!    IF (iGP == 49 .and. jGP == 2 .and. kGP == 2) then
!       WRITE(355,*)distGPEIMin
!       WRITE(355,*)'ZONE'
!       WRITE(355,*)xGP,yGP,zGP
!       WRITE(355,*)xBI,yBI,zBI
!    ENDIF
!DEBUG

! TEMPORARY
!    WRITE(365,*)xBI,yBI,zBI
! END TEMPORARY

999    CONTINUE 

       xBIT(nc) = xBI
       yBIT(nc) = yBI
       zBIT(nc) = zBI
       cElementGP(nc) = closestElementGP

       if (is_in2.eq.1)then
          normX(nc) = ElemNorm_x(closestElementGP)
          normY(nc) = ElemNorm_y(closestElementGP)
          normZ(nc) = ElemNorm_z(closestElementGP)
       else
          normX(nc) = NormXtemp2
          normY(nc) = NormYtemp2
          normZ(nc) = NormZtemp2
       endif

       ENDDO  ! nc

       !Ye,find the outside projection BodyMarkers===
       IF (fg.eq.12)THEN
       if (is_in.eq.1)then
          cnt_flg12(iPt) =  1
       else
          cnt_flg12(iPt) = -1
       endif
       ENDIF
       IF (fg.eq.21)THEN
       if (is_in.eq.1)then
          cnt_flg21(iPt) =  1
       else
          cnt_flg21(iPt) = -1
       endif
       ENDIF
       IF (fg.eq.13)THEN
       if (is_in.eq.1)then
          cnt_flg13(iPt) =  1
       else
          cnt_flg13(iPt) = -1
       endif
       ENDIF
       IF (fg.eq.31)THEN
       if (is_in.eq.1)then
          cnt_flg31(iPt) =  1
       else
          cnt_flg31(iPt) = -1
       endif
       ENDIF
       IF (fg.eq.23)THEN
       if (is_in.eq.1)then
          cnt_flg23(iPt) =  1
       else
          cnt_flg23(iPt) = -1
       endif
       ENDIF
       IF (fg.eq.32)THEN
       if (is_in.eq.1)then
          cnt_flg32(iPt) =  1
       else
          cnt_flg32(iPt) = -1
       endif
       ENDIF
       !===========================================
      
       !Ye,debug
       !if ((iPt.eq.2549).and.(fg.eq.12))then
       !   write(601,*)normX(1),normY(1),normZ(1)
       !   write(601,*)normX(2),normY(2),normZ(2)
       !   write(601,*)normX(3),normY(3),normZ(3)
       !   write(601,*)'---'
       !endif

       DO nc = 1,nCheck
          dist(nc) =  (xBIT(nc)-xGP)**2 
     &              + (yBIT(nc)-yGP)**2 
     &              + (zBIT(nc)-zGP)**2
       ENDDO

       iDummy           = MINLOC(dist(1:nCheck))
       shortestProbe    = iDummy(1)
       xBI              = xBIT(shortestProbe)
       yBI              = yBIT(shortestProbe)
       zBI              = zBIT(shortestProbe)
       closestElementGP = cElementGP(shortestProbe)

       !Ye,debug
       !if ((iPt.eq.1415))then
       !   write(500,*)fg
       !   write(500,*)xGP,yGP,zGP
       !   write(500,*)xBIT(1:nCheck)
       !   write(500,*)yBIT(1:nCheck)
       !   write(500,*)zBIT(1:nCheck)
       !   write(500,*)dist(1:nCheck)
       !   write(500,*)iDummy(1)
       !   write(500,*)'-----'
       !endif

       !Ye,debug
       !if ((iPt.eq.2549).and.(fg.eq.12))then
       !   write(602,*)shortestProbe
       !   write(602,*)normX(1),normY(1),normZ(1)
       !   write(602,*)normX(2),normY(2),normZ(2)
       !   write(602,*)normX(3),normY(3),normZ(3)
       !   write(602,*)'---'
       !endif

       BI_normX(iPt) = normX(shortestProbe)
       BI_normY(iPt) = normY(shortestProbe)
       BI_normZ(iPt) = normZ(shortestProbe)

       !Ye,debug
       !if ((iPt.eq.2549).and.(fg.eq.12))then
       !   write(603,*)shortestProbe
       !   write(603,*)BI_normX(iPt),BI_normY(iPt),BI_normZ(iPt)
       !   write(603,*)'---'
       !endif

       if (is_in.eq.1)then
          cnt_d =  (xBI-xGP)*ElemNorm_x(closestElementGP)
     &            +(yBI-yGP)*ElemNorm_y(closestElementGP)
     &            +(zBI-zGP)*ElemNorm_z(closestElementGP)
       else
          cnt_d = SQRT((xBI-xGP)**2 +
     &                 (yBI-yGP)**2 +
     &                 (zBI-zGP)**2 )
       endif

 
918    CONTINUE

       ElemNum  = closestElementGP

!      Ye,debug
!       if ((iPt.eq.5003).and.(fg.eq.21))then
!          write(555,*)is_in,closestElementGP,closestNodeGP(1:NCheck)
!          write(555,*)xBI,yBI,zBI
!          write(555,*)xGP,yGP,zGP
!          write(555,*)'shortestProbe:',shortestProbe
!          write(555,*)dist(1:nCheck)
!          write(555,*)ElemNorm_x(closestElementGP),
!     &                ElemNorm_y(closestElementGP),
!     &                ElemNorm_z(closestElementGP)
!          write(555,*)BI_normX(iPt),BI_normY(iPt),BI_normZ(iPt)
!          write(555,*)cnt_d
!          write(555,*)'----'
!       endif

       !Ye,test
       !IF (ABS(xBI-xGP).GE.0.05)THEN
       !   cnt_d = 1.0E+20
       !ENDIF

       END SUBROUTINE calc_bodyIntercept_Unstruc_cnt
!------------------------------------------------------------------------------

      SUBROUTINE check_BIInside_cnt(xVert,yVert,zVert,xNorm,yNorm,zNorm,
     &           xBITemp,yBITemp,zBITemp,area123,areaDiff)
!------------------------------------------------------------------------------

      USE contact
 
      IMPLICIT NONE

!... parameters

      REAL*8, DIMENSION(3), INTENT(IN) :: xVert,yVert,zVert
      REAL*8, INTENT(IN)               :: xNorm,yNorm,zNorm
      REAL*8, INTENT(IN)               :: xBITemp,yBITemp,zBITemp
      REAL*8, INTENT(OUT)              :: area123,areaDiff

!... loop variables

      INTEGER :: iside

!... local variables

      REAL*8        :: distanceToPlane,distPointToPlane
      REAL*8        :: side12,side23,side31,side14,side24,side34
      REAL*8        :: area124,area234,area314
      REAL*8        :: semiPerimeter123,semiPerimeter124
      REAL*8        :: semiPerimeter234,semiPerimeter314

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

      side12 =  SQRT( (xVert(2)-xVert(1))**2  
     &               +(yVert(2)-yVert(1))**2  
     &               +(zVert(2)-zVert(1))**2  )
      side23 =  SQRT( (xVert(3)-xVert(2))**2  
     &               +(yVert(3)-yVert(2))**2  
     &               +(zVert(3)-zVert(2))**2  )
      side31 =  SQRT( (xVert(1)-xVert(3))**2  
     &               +(yVert(1)-yVert(3))**2  
     &               +(zVert(1)-zVert(3))**2  )
      side14 =  SQRT( (xBITemp-xVert(1))**2  
     &               +(yBITemp-yVert(1))**2  
     &               +(zBITemp-zVert(1))**2  )
      side24 =  SQRT( (xBITemp-xVert(2))**2  
     &               +(yBITemp-yVert(2))**2  
     &               +(zBITemp-zVert(2))**2  )
      side34 =  SQRT( (xBITemp-xVert(3))**2  
     &               +(yBITemp-yVert(3))**2  
     &               +(zBITemp-zVert(3))**2  )

      semiPerimeter123 = 0.5d0*(side12 + side23 + side31)
      semiPerimeter124 = 0.5d0*(side12 + side24 + side14)
      semiPerimeter234 = 0.5d0*(side23 + side24 + side34)
      semiPerimeter314 = 0.5d0*(side31 + side34 + side14)
 
      area123  = SQRT( semiPerimeter123*abs(semiPerimeter123-side12)
     &                 *abs(semiPerimeter123-side23)
     &                 *abs(semiPerimeter123-side31)   )

      area124  = SQRT( semiPerimeter124*abs(semiPerimeter124-side12)
     &                 *abs(semiPerimeter124-side24)
     &                 *abs(semiPerimeter124-side14)   )

      area234  = SQRT( semiPerimeter234*abs(semiPerimeter234-side23) 
     &                 *abs(semiPerimeter234-side24) 
     &                 *abs(semiPerimeter234-side34)   )

      area314  = SQRT( semiPerimeter314*abs(semiPerimeter314-side31) 
     &                 *abs(semiPerimeter314-side34) 
     &                 *abs(semiPerimeter314-side14)   )

      areaDiff  = area124 + area234 + area314 - area123

      END SUBROUTINE check_BIInside_cnt
!===================================================================
      SUBROUTINE calc_BIOutside_cnt(xVert,yVert,zVert,xCG,yCG,zCG,  
     &                              xBITemp,yBITemp,zBITemp,distBIElem)
!------------------------------------------------------------------------------

      USE contact

      IMPLICIT NONE

!... parameters

      REAL*8, DIMENSION(3), INTENT(IN) :: xVert,yVert,zVert
      REAL*8, INTENT(IN)               :: xCG,yCG,zCG
      REAL*8, INTENT(IN)  :: xBITemp,yBITemp,zBITemp
      REAL*8, INTENT(OUT) :: distBIElem

!... loop variables

      INTEGER :: iside

!... local variables

      INTEGER :: isideSelect,mside
      INTEGER :: nodeSelect1,nodeSelect2

      REAL*8 :: aCrossbVectMagn,dotVal
      REAL*8 :: distIntBI,distIntCG,distIntMin
      REAL*8 :: distNorm,distNode1BINorm,distNode2BINorm
      REAL*8 :: distVert1BI,distVert2BI
      REAL*8 :: magnitude12,magnitudeBICG,magnitude,projectedLength
      REAL*8 :: node12x,node12y,node12z
      REAL*8 :: vec01x,vec01y,vec01z
      REAL*8 :: vec12x,vec12y,vec12z
      REAL*8 :: xBINorm,yBINorm,zBINorm
      REAL*8, DIMENSION(3) :: aVect,bVect,cVect,aCrossbVect,cCrossbVect
      REAL*8, DIMENSION(3) :: xInt,yInt,zInt
      REAL*8, DIMENSION(3) :: vect1,vect2,vect3,vect4
      REAL*8, DIMENSION(3,3) :: vectInt
      REAL*8 :: epstol
!******************************************************************************

! ============================================================================
!   Construct Intersection points between 
!     line linking BI and Centroid of Surface Element and the triangle sides
!     L1: BI-CG, L2:Sides of Vertices
!
!   use formula for intersection point between 2 co-planar lines from Mathworld
!   http://mathworld.wolfram.com/Line-LineIntersect.html
!
!             x4
!             *
!             |
!             |
!    x1       | Int      x2
!    *--------*----------*
!             |
!             |  
!             | 
!             * x3
!
! ============================================================================

        vect1(1:3) = (/xBITemp,yBITemp,zBITemp/)
        vect2(1:3) = (/xCG,yCG,zCG/)

       aVect(1:3) = vect2(1:3)-vect1(1:3)

      DO iside = 1,3
         mside = iside +1
         IF(iside == 3) mside = 1

        vect3(1:3) =(/xVert(iside),yVert(iside),zVert(iside)/)
        vect4(1:3) =(/xVert(mside),yVert(mside),zVert(mside)/)

      bVect(1:3)  = vect4(1:3) -vect3(1:3)
      cVect(1:3)  = vect3(1:3) -vect1(1:3)

      call calc_crossProduct_cnt(aVect,bVect,aCrossbVect)
      call calc_crossProduct_cnt(cVect,bVect,cCrossbVect)

      aCrossbVectMagn = aCrossbVect(1)**2 + aCrossbVect(2)**2 +
     &                  aCrossbVect(3)**2

      dotVal = DOT_PRODUCT(cCrossbVect,aCrossbVect)

      vectInt(1:3,iside)=vect1(1:3)+aVect(1:3)*dotVal/aCrossbVectMagn
      END DO ! iside 

! ============================================================================
!   Choose closest intersection point lying between BI and CG
!     Normalsize value with L1
! ============================================================================
      magnitudeBICG = SQRT( (vect1(1)-vect2(1))**2 
     &                    + (vect1(2)-vect2(2))**2 
     &                    + (vect1(3)-vect2(3))**2 )

      distIntMin = 1.0E+16
      isideSelect = -1000
      epstol = 1.0e-3 !1.0E-6_cgreal

      DO iside = 1,3
      distIntBI = SQRT( (vect1(1)-vectInt(1,iside))**2 
     &                + (vect1(2)-vectInt(2,iside))**2 
     &                + (vect1(3)-vectInt(3,iside))**2 )/magnitudeBICG

      distIntCG = SQRT( (vect2(1)-vectInt(1,iside))**2 
     &                + (vect2(2)-vectInt(2,iside))**2 
     &                + (vect2(3)-vectInt(3,iside))**2 )/magnitudeBICG

      IF(distIntBI <= 1.0d0+epstol .AND. distIntCG <= 1.0d0+epstol)THEN
        distIntMin  = DMIN1(distIntBI,distIntCG)
        isideSelect = iside
      END IF ! distIntBI

      END DO ! iside 

! ============================================================================
!   Trap error for isideSelect 
! ============================================================================

      IF ( isideSelect < 0 ) THEN
      WRITE(*,*)
     & 'CNT_BIdistMin: Incorrect selection of iside (Should be either 1,
     & 2 or 3'
      WRITE(*,*)
     &  'default value selected = ',isideSelect
      WRITE(*,*) xVert(1),yVert(1),zVert(1)
      WRITE(*,*) xVert(2),yVert(2),zVert(2)
      WRITE(*,*) xVert(3),yVert(3),zVert(3)
      WRITE(*,*) vect1(1:3)
      WRITE(*,*) vect2(1:3)

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

      vec12x = xVert(nodeSelect2) - xVert(nodeSelect1)
      vec12y = yVert(nodeSelect2) - yVert(nodeSelect1)
      vec12z = zVert(nodeSelect2) - zVert(nodeSelect1)

      magnitude12 = SQRT(vec12x**2 + vec12y**2 + vec12z**2)

      vec01x = xVert(nodeSelect1) - xBITemp
      vec01y = yVert(nodeSelect1) - yBITemp
      vec01z = zVert(nodeSelect1) - zBITemp

      distNorm = SQRT(  (vec12y*vec01z - vec12z*vec01y)**2  
     &           + (vec12z*vec01x - vec01z*vec12x)**2  
     &           + (vec12x*vec01y - vec01x*vec12y)**2  )/magnitude12

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

      projectedLength = (xBITemp - xVert(nodeSelect1))*node12x 
     &                 +(yBITemp - yVert(nodeSelect1))*node12y 
     &                 +(zBITemp - zVert(nodeSelect1))*node12z

      xBINorm = xVert(nodeSelect1) + projectedLength*node12x
      yBINorm = yVert(nodeSelect1) + projectedLength*node12y
      zBINorm = zVert(nodeSelect1) + projectedLength*node12z

! ============================================================================
!    Determine distance between BINorm and vertices of selected side.
!     If normal point lies inside the side, select that distance.
!     If it lies outside, find the minimum distance with vertices
!     Use normalized length
! ============================================================================

      distNode1BINorm = SQRT( (xVert(nodeSelect1)-xBINorm)**2 
     &                  + (yVert(nodeSelect1)-yBINorm)**2 
     &                  + (zVert(nodeSelect1)-zBINorm)**2 )/magnitude

! H. Luo changed distNode1BINorm to distNode2BINorm in the following line.

      distNode2BINorm = SQRT( (xVert(nodeSelect2)-xBINorm)**2 
     &                  + (yVert(nodeSelect2)-yBINorm)**2 
     &                  + (zVert(nodeSelect2)-zBINorm)**2 )/magnitude

      IF(distNode1BINorm <= 1.0d0 .AND. distNode2BINorm <= 1.0d0)THEN
        distBIElem  =  distNorm
      ELSE
        distVert1BI = SQRT( (xVert(nodeSelect1)-xBITemp)**2 
     &                    + (yVert(nodeSelect1)-yBITemp)**2 
     &                    + (zVert(nodeSelect1)-zBITemp)**2 )

        distVert2BI = SQRT( (xVert(nodeSelect2)-xBITemp)**2 
     &                    + (yVert(nodeSelect2)-yBITemp)**2 
     &                    + (zVert(nodeSelect2)-zBITemp)**2 )

        distBIElem = DMIN1(distVert1BI,distVert2BI)
      END IF ! distNode1BINorm

      END SUBROUTINE calc_BIOutside_cnt
!==============================================================
      SUBROUTINE calc_crossProduct_cnt(r,s,cross_product)

      IMPLICIT NONE

      REAL*8, DIMENSION(3), INTENT(IN)  :: r,s
      REAL*8, DIMENSION(3), INTENT(OUT) :: cross_product

      INTEGER :: component,i,j

      DO component = 1,3
         i = MODULO(component,3) + 1
         j = MODULO(i,3) + 1
         cross_product(component) = r(i)*s(j) - s(i)*r(j)
      END DO ! component 

      END SUBROUTINE calc_crossProduct_cnt
