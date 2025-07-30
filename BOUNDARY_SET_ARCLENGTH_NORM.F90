!------------------------------------------------------------------------------
   SUBROUTINE calculate_arclength_norm_ds()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

!... Loop variables
    INTEGER           :: iBody,m, mMax, mMin, cNode, cElement
    INTEGER           :: node1, node2, node3
    REAL(KIND=CGREAL) :: tangX, tangY, distNode, dMinUnstruc
    REAL(KIND=CGREAL) :: normMag, tangMag
    REAL(KIND=CGREAL) :: vectorx,vectory,vectorz,vectProduct
    REAL(KIND=CGREAL) :: triElemVectorX(2),triElemVectorY(2),triElemVectorZ(2)
    !REAL(KIND=CGREAL) :: dxMin,dyMin,dzMin
    REAL(KIND=CGREAL) :: cellAreaMin,triElemAreaAv,bodyResolutionNormalized


!... Local variables

! check direction of surface normal and adjust ordering of triangle if normal does not point into body
! this is done only during fresh start
    IF( ntime == 0 ) THEN
       DO iBody = 1, nBody
!          IF(unstruc_surface_type(iBody) == SOLID_BODY) THEN
             normDirFlag = 1.0_CGREAL
!             PRINT*,'CALL determine_norm_dir_unstruc(iBody)'
             call check_orientation_unstruc(iBody)
             CALL determine_norm_dir_unstruc(iBody)
!          ENDIF
       ENDDO ! iBody
    ENDIF
 
!    PRINT*,'Computing Surface Quantities for Unstructured surface'
    sBodyMarker = 0.0_CGREAL
    surfArea    = 0.0_CGREAL

    ! initialize the surface normal at the marker points
    xNormBodyMarker = 0.0_CGREAL
    yNormBodyMarker = 0.0_CGREAL
    zNormBodyMarker = 0.0_CGREAL

    DO iBody = nBody,1,-1

       !--Sweep surface of each triangular element
       surfArea(iBody)  = 0.0_CGREAL
                                    
       DO m=1,totNumTriElem(iBody)

          node1 = triElemNeig(1,m,iBody)
          node2 = triElemNeig(2,m,iBody)
          node3 = triElemNeig(3,m,iBody)

          !--First vector of each element
           triElemVectorx(1)= xBodyMarker(node2,iBody) &
                             -xBodyMarker(node1,iBody)
           triElemVectory(1)= yBodyMarker(node2,iBody) &
                             -yBodyMarker(node1,iBody)
           triElemVectorz(1)= zBodyMarker(node2,iBody) &
                             -zBodyMarker(node1,iBody)

           !--Second vector of each element
           triElemVectorx(2)= xBodyMarker(node3,iBody) &
                             -xBodyMarker(node2,iBody)
           triElemVectory(2)= yBodyMarker(node3,iBody) &
                             -yBodyMarker(node2,iBody)
           triElemVectorz(2)= zBodyMarker(node3,iBody) &
                             -zBodyMarker(node2,iBody)

           !--Normal vector of the element
           triElemNormx(m,iBody)= (  triElemVectory(1)*triElemVectorz(2)  &
                                   - triElemVectorz(1)*triElemVectory(2)  )
           triElemNormy(m,iBody)= (  triElemVectorz(1)*triElemVectorx(2)  &
                                   - triElemVectorx(1)*triElemVectorz(2)  )
           triElemNormz(m,iBody)= (  triElemVectorx(1)*triElemVectory(2)  &
                                   - triElemVectory(1)*triElemVectorx(2)  )

           normMag  = SQRT(  triElemNormx(m,iBody)**2   &
                           + triElemNormy(m,iBody)**2   &
                           + triElemNormz(m,iBody)**2     )

           ! Area of element
           triElemArea(m,iBody)  = 0.5_CGREAL*normMag
           surfArea(iBody)       = surfArea(iBody) + triElemArea(m,iBody)

           ! Unit Normal vector
           triElemNormx(m,iBody) = triElemNormx(m,iBody)/normMag
           triElemNormy(m,iBody) = triElemNormy(m,iBody)/normMag
           triElemNormz(m,iBody) = triElemNormz(m,iBody)/normMag

           xNormBodyMarker(node1,iBody) = xNormBodyMarker(node1,iBody) + triElemNormx(m,iBody) 
           xNormBodyMarker(node2,iBody) = xNormBodyMarker(node2,iBody) + triElemNormx(m,iBody)
           xNormBodyMarker(node3,iBody) = xNormBodyMarker(node3,iBody) + triElemNormx(m,iBody)

           yNormBodyMarker(node1,iBody) = yNormBodyMarker(node1,iBody) + triElemNormy(m,iBody) 
           yNormBodyMarker(node2,iBody) = yNormBodyMarker(node2,iBody) + triElemNormy(m,iBody)
           yNormBodyMarker(node3,iBody) = yNormBodyMarker(node3,iBody) + triElemNormy(m,iBody)

           zNormBodyMarker(node1,iBody) = zNormBodyMarker(node1,iBody) + triElemNormz(m,iBody) 
           zNormBodyMarker(node2,iBody) = zNormBodyMarker(node2,iBody) + triElemNormz(m,iBody)
           zNormBodyMarker(node3,iBody) = zNormBodyMarker(node3,iBody) + triElemNormz(m,iBody)

           ! Unit Tangents
           ! Tangent-2 defined parallel to vector from vertex-1 to vertex-2

           triElemTang2x(m,iBody)  = xBodyMarker(node2,iBody) &
                                    -xBodyMarker(node1,iBody)
           triElemTang2y(m,iBody)  = yBodyMarker(node2,iBody) &
                                    -yBodyMarker(node1,iBody)
           triElemTang2z(m,iBody)  = zBodyMarker(node2,iBody) &
                                    -zBodyMarker(node1,iBody)

           tangMag     = SQRT(  triElemTang2x(m,iBody)**2   &
                              + triElemTang2y(m,iBody)**2   &
                              + triElemTang2z(m,iBody)**2 )

           triElemTang2x(m,iBody)  = triElemTang2x(m,iBody)/tangMag
           triElemTang2y(m,iBody)  = triElemTang2y(m,iBody)/tangMag
           triElemTang2z(m,iBody)  = triElemTang2z(m,iBody)/tangMag

           ! t1 = t2 x n

           triElemTang1x(m,iBody)  =  triElemTang2y(m,iBody)*triElemNormz(m,iBody)  &
                                    - triElemTang2z(m,iBody)*triElemNormy(m,iBody)
           triElemTang1y(m,iBody)  =- triElemTang2x(m,iBody)*triElemNormz(m,iBody)  &
                                    + triElemTang2z(m,iBody)*triElemNormx(m,iBody)
           triElemTang1z(m,iBody)  =  triElemTang2x(m,iBody)*triElemNormy(m,iBody)  &
                                    - triElemTang2y(m,iBody)*triElemNormx(m,iBody)

           !--Centroid of the each element
           triElemCentx(m,iBody)=(xBodyMarker(node1,iBody) &
                                 +xBodyMarker(node2,iBody) &
                                 +xBodyMarker(node3,iBody))/3.0_CGREAL
           triElemCenty(m,iBody)=(yBodyMarker(node1,iBody) &
                                 +yBodyMarker(node2,iBody) &
                                 +yBodyMarker(node3,iBody))/3.0_CGREAL
           triElemCentz(m,iBody)=(zBodyMarker(node1,iBody) &
                                 +zBodyMarker(node2,iBody) &
                                 +zBodyMarker(node3,iBody))/3.0_CGREAL

           !WRITE(800+iBody,'(I6,1X,6F12.5)') m,triElemCentx(m,iBody),triElemCenty(m,iBody),triElemCentz(m,iBody), &
           !                                  triElemNormx(m,iBody),triElemNormy(m,iBody),triElemNormz(m,iBody)
        ENDDO ! m
        !stop
 
        !dxMin                    = MINVAL(dx(1:nx-1))
        !dyMin                    = MINVAL(dy(1:ny-1))
        !dzMin                    = MINVAL(dz(1:nz-1))
        cellAreaMin              = (dxMin*dyMIn*dzMin)**(2.0_CGREAL/3.0_CGREAL)
        triElemAreaAv            = surfArea(iBody)/totNumTriElem(iBody)
        bodyResolutionNormalized = cellAreaMin/triElemAreaAv
!        print*,'Min Cell Area                             = ',cellAreaMin
!        print*,'Surface area of body',ibody,            ' = ',surfArea(iBody)
!        print*,'Average Element Area for',ibody,        ' = ',triElemAreaAv
!        print*,'Norm. Surface resolution of body',ibody,' = ',bodyResolutionNormalized

        DO m=1,nPtsBodyMarker(iBody)
           normMag  = SQRT(xNormBodyMarker(m,iBody)**2   &
                         + yNormBodyMarker(m,iBody)**2   &
                         + zNormBodyMarker(m,iBody)**2 )

           ! Unit Normal vector
           xNormBodyMarker(m,iBody) = xNormBodyMarker(m,iBody)/normMag
           yNormBodyMarker(m,iBody) = yNormBodyMarker(m,iBody)/normMag
           zNormBodyMarker(m,iBody) = zNormBodyMarker(m,iBody)/normMag
           !WRITE(900+iBody,'(I6,1X,8F12.5)') m,                                     &
           !         xBodyMarker(m,iBody),yBodyMarker(m,iBody),zBodyMarker(m,iBody), &
           !         xNormBodyMarker(m,iBody),yNormBodyMarker(m,iBody),zNormBodyMarker(m,iBody)
        ENDDO
        !close(900)

    ENDDO ! iBody

    !!output BodyMarker norm info
    !IF (ntime==ntime_start+1)THEN
    !DO iBody = nBody,1,-1
    !   DO m=1,nPtsBodyMarker(iBody)
    !      write(513,'(6F12.5)')xBodyMarker(m,iBody),yBodyMarker(m,iBody),   &
    !           zBodyMarker(m,iBody),xNormBodyMarker(m,iBody),        &
    !           yNormBodyMarker(m,iBody),zNormBodyMarker(m,iBody)
    !   ENDDO
    !ENDDO
    !ENDIF


   END SUBROUTINE calculate_arclength_norm_ds
!
!
!------------------------------------------------------------------------------
! Subroutine determines if for a given unstructured surface mesh,
! the norm points in or out and also reorders vertex numbers

! Note: our convention is that positive normal vector points into body
! Per convention of solid geometry, if three vertices (1,2,3) of a triangle are
! ordered clockwise when viewed from one side, then (p2-p1)x(p3-p1) produces a vector
! that points out towards the direction from where it is being viewed. If we force
! the normal to point inwards, then we should also change the vertex ordering
!
!------------------------------------------------------------------------------

   SUBROUTINE determine_norm_dir_unstruc(iBody)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER,INTENT(IN)::iBody

!... Loop variables
    INTEGER           :: m,cNode, cElement

!... Local variables
    INTEGER           :: node1,node2,node3
    REAL(KIND=CGREAL) :: distNode, dMinUnstruc
    REAL(KIND=CGREAL) :: vectorx,vectory,vectorz,vectProduct
    REAL(KIND=CGREAL) :: triElemVectorX(2),triElemVectorY(2),triElemVectorZ(2)
    REAL(KIND=CGREAL) :: cTriElemNormx,cTriElemNormy,cTriElemNormz
    REAL(KIND=CGREAL) :: cTriElemCentx,cTriElemCenty,cTriElemCentz


    cNode = 0

    !---Find node closest to outside point
    dMinUnstruc = 1.0E8_CGREAL
    DO m=1,nPtsBodyMarker(iBody)
       distNode= SQRT( (xBodyMarker(m,iBody)-pointOutsideBodyX(iBody))**2  &
                      +(yBodyMarker(m,iBody)-pointOutsideBodyY(iBody))**2  &
                      +(zBodyMarker(m,iBody)-pointOutsideBodyZ(iBody))**2   )
       IF(distNode <= dMinUnstruc) THEN
         dMinUnstruc  = distNode
         cNode        = m
       ENDIF
    ENDDO

    !---Find element corresponing to closest node
    DO m=1,totNumTriElem(iBody)
       IF ( triElemNeig(1,m,iBody) == cNode .OR. &
            triElemNeig(2,m,iBody) == cNode .OR. &
            triElemNeig(3,m,iBody) == cNode       ) cElement = m
    ENDDO

!RRRRRRRRRRRRRRRRRRRRRRRRRRR
!         print*,'Closest Node to outside point             = ',cNode
!         print*,'Element Corresponding to closest node is = ',cElement
!RRRRRRRRRRRRRRRRRRRRRRRRRRR


!     1 *-------------* 3
!        \           /
!         \         /
!          \       /
!           \     /
!            \   /
!             \ /
!            2 *
!
!--Sweep surface of CLOSEST element to determine normal direction of triangular elements

       node1 = triElemNeig(1,cElement,iBody)
       node2 = triElemNeig(2,cElement,iBody)
       node3 = triElemNeig(3,cElement,iBody)


       triElemVectorx(1)= xBodyMarker(node2,iBody) &
                         -xBodyMarker(node1,iBody)
       triElemVectory(1)= yBodyMarker(node2,iBody) &
                         -yBodyMarker(node1,iBody)
       triElemVectorz(1)= zBodyMarker(node2,iBody) &
                         -zBodyMarker(node1,iBody)

       triElemVectorx(2)= xBodyMarker(node3,iBody) &
                         -xBodyMarker(node2,iBody)
       triElemVectory(2)= yBodyMarker(node3,iBody) &
                         -yBodyMarker(node2,iBody)
       triElemVectorz(2)= zBodyMarker(node3,iBody) &
                         -zBodyMarker(node2,iBody)

!-- Normal of closest element
       cTriElemNormx = ( triElemVectory(1)*triElemVectorz(2)  &
                        -triElemVectorz(1)*triElemVectory(2) )
       cTriElemNormy = ( triElemVectorz(1)*triElemVectorx(2)  &
                        -triElemVectorx(1)*triElemVectorz(2) )
       cTriElemNormz = ( triElemVectorx(1)*triElemVectory(2)  &
                        -triElemVectory(1)*triElemVectorx(2) )

!--Centroid of the closest element
       cTriElemCentx = (xBodyMarker(node1,iBody) &
                       +xBodyMarker(node2,iBody) &
                       +xBodyMarker(node3,iBody))/3.0_CGREAL
       cTriElemCenty = (yBodyMarker(node1,iBody) &
                       +yBodyMarker(node2,iBody) &
                       +yBodyMarker(node3,iBody))/3.0_CGREAL
       cTriElemCentz = (zBodyMarker(node1,iBody) &
                       +zBodyMarker(node2,iBody) &
                       +zBodyMarker(node3,iBody))/3.0_CGREAL

       vectorx = cTriElemCentX - pointOutsideBodyX(iBody)
       vectory = cTriElemCentY - pointOutsideBodyY(iBody)
       vectorz = cTriElemCentZ - pointOutsideBodyZ(iBody)

       vectProduct = vectorx*cTriElemNormx  &
                   + vectory*cTriElemNormy  &
                   + vectorz*cTriElemNormz

!RRRRRRRRRRRRRRRRRRRRRRRRRRRRr
!        print*,'vectProduct = ',vectproduct
!RRRRRRRRRRRRRRRRRRRRRRRRRRRRr

       normDirFlag = 1.0_CGREAL

       IF (vectProduct < 0.0_CGREAL)  THEN
         normDirFlag =-1.0_CGREAL
!         print*,'Reordering triangle vertices since normal points out of the solid body'
         DO m=1,totNumTriElem(iBody)         ! changing vertex ordering
           node2 = triElemNeig(2,m,iBody)
           node3 = triElemNeig(3,m,iBody)
           triElemNeig(2,m,iBody) = node3
           triElemNeig(3,m,iBody) = node2
         ENDDO
       ENDIF

   END SUBROUTINE determine_norm_dir_unstruc
!------------------------------------------------------------------------------
   SUBROUTINE check_orientation_unstruc(iBody)

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE boundary_arrays

    IMPLICIT NONE

    INTEGER,INTENT(IN)::iBody

!... Loop variables
    INTEGER           :: m,n,i,j,info

!... Local variables
    INTEGER           :: node1(3),node2(3),v1(2),v2(2),i1,d1,d2,nodetemp2,nodetemp3,eci
    INTEGER, DIMENSION(:), ALLOCATABLE  :: queue,iNorm
    INTEGER, DIMENSION(:,:), ALLOCATABLE  :: ElemConnect 
    REAL(KIND=CGREAL) :: distNode, dMinUnstruc
    REAL(KIND=CGREAL) :: vectorx,vectory,vectorz,vectProduct
    REAL(KIND=CGREAL) :: triElemVectorX(2),triElemVectorY(2),triElemVectorZ(2)
    REAL(KIND=CGREAL) :: cTriElemNormx,cTriElemNormy,cTriElemNormz
    REAL(KIND=CGREAL) :: cTriElemCentx,cTriElemCenty,cTriElemCentz



!     1 *-------------* 3
!        \           /
!         \         /
!          \       /
!           \     /
!            \   /
!             \ /
!            2 *
!
!
    allocate(ElemConnect(3,totNumTriElem(iBody)))
    allocate(queue(totNumTriElem(iBody)))
    allocate(iNorm(totNumTriElem(iBody)))

    iNorm=0


    DO m=1,totNumTriElem(iBody)

       node1(1) = triElemNeig(1,m,iBody)
       node1(2) = triElemNeig(2,m,iBody)
       node1(3) = triElemNeig(3,m,iBody)

        eci=0
        do n=1,totNumTriElem(iBody)

            if( (m.ne.n) )then !.and.(iswitch(m).eq.1 ) 
                
                node2(1) = triElemNeig(1,n,iBody)
                node2(2) = triElemNeig(2,n,iBody)
                node2(3) = triElemNeig(3,n,iBody)
    
                i1=0
                do i=1,3
                do j=1,3                
                    if(node1(i).eq.node2(j)) then
                        i1=i1+1
                        v1(i1)=i    
                        v2(i1)=j                      
                    endif
                enddo
                enddo

                if (i1.gt.2) then
                    print *,m,n,i1
                    stop
                endif

                if(i1.eq.2) then
!                    print *,m,n                   
                    eci=eci+1
                    ElemConnect(eci,m)=n
                    if(eci.eq.3) goto 888

                endif

            endif


        enddo

888     continue

    enddo


    call FIFO_INIT(m,queue)
    
    m=1
    call FIFO_ADD(m,totNumTriElem(iBody),queue)
    iNorm(m)=1

    call FIFO_EMPTY(totNumTriElem(iBody),queue,info)          
    
    do while (info.eq.0)

        call FIFO_REMOVE(m,totNumTriElem(iBody),queue)

        node1(1) = triElemNeig(1,m,iBody)
        node1(2) = triElemNeig(2,m,iBody)
        node1(3) = triElemNeig(3,m,iBody)

        do eci=1,3
            n=ElemConnect(eci,m)

            !By Ye Chen
            !n=0 if the element vertices have correct ordering.
            !Need to skip the following block to avoid segmentation error.
            if (n.eq.0) goto 777

            if( (iNorm(n).eq.0) .and. (n.gt.0) ) then
                node2(1) = triElemNeig(1,n,iBody)
                node2(2) = triElemNeig(2,n,iBody)
                node2(3) = triElemNeig(3,n,iBody)

                i1=0
                do i=1,3
                do j=1,3                
                    if(node1(i).eq.node2(j)) then
                        i1=i1+1
                        v1(i1)=i    
                        v2(i1)=j                      
                    endif
                enddo
                enddo


                if(i1.eq.2) then
!                    print *,m,n
                    d1=0
                    if(((v1(1).eq.1).and.(v1(2).eq.2) ).or. &
                       ((v1(1).eq.2).and.(v1(2).eq.3) ).or. &
                       ((v1(1).eq.3).and.(v1(2).eq.1) ) )then
                            d1=1
                    endif
                    d2=0
                    if(((v2(1).eq.1).and.(v2(2).eq.2) ).or. &
                       ((v2(1).eq.2).and.(v2(2).eq.3) ).or. &
                       ((v2(1).eq.3).and.(v2(2).eq.1) ) )then
                            d2=1
                    endif

                    if(d1.eq.d2) then
!                       print*,'Reordering triangle vertices',m,n
                       nodetemp2 = triElemNeig(2,n,iBody)
                       nodetemp3 = triElemNeig(3,n,iBody)
                       triElemNeig(2,n,iBody) = nodetemp3
                       triElemNeig(3,n,iBody) = nodetemp2
                    endif
                endif

                iNorm(n)=1
                call FIFO_ADD(n,totNumTriElem(iBody),queue)

            endif

777     continue
        enddo

        call FIFO_EMPTY(totNumTriElem(iBody),queue,info)

        if (info.eq.1) then
        DO m=1,totNumTriElem(iBody)
            if( iNorm(m).eq.0 ) then
                print *,'Element', m,' unchanged'
                call FIFO_ADD(m,totNumTriElem(iBody),queue)
                info=0
                goto 999
            endif
        enddo
999     continue
        endif

    enddo      

    DO m=1,totNumTriElem(iBody)
        if( iNorm(m).eq.0 ) print *,'Element', m,'remained unchanged!!!'
    enddo

   END SUBROUTINE check_orientation_unstruc

    
!---------------------------------------------------------------------
      SUBROUTINE FIFO_INIT(NP,QUEUE)
 
      IMPLICIT none
 
        INTEGER NP
        INTEGER QUEUE(NP)
 
 
        QUEUE(1)=3
        QUEUE(2)=3
        RETURN
        END
 
!---------------------------------------------------------------------
      SUBROUTINE FIFO_ADD(IPOS,NP,QUEUE)
 
      IMPLICIT none
 
        INTEGER NP,IPOS,QFIRST,QLAST
        INTEGER QUEUE(NP)
 
        QFIRST=QUEUE(1)
        QLAST=QUEUE(2)
 

        QUEUE(QLAST)=IPOS
        IF(QLAST.EQ.NP) THEN
           QLAST=3
        ELSE
           QLAST=QLAST+1
        ENDIF
        QUEUE(2)=QLAST
        RETURN
        END
 
!---------------------------------------------------------------------
      SUBROUTINE FIFO_REMOVE(IPOS,NP,QUEUE)
 
      IMPLICIT none
 
        INTEGER NP,IPOS,QFIRST,QLAST
        INTEGER QUEUE(NP)
 
        QFIRST=QUEUE(1)
        QLAST=QUEUE(2)
 
        IPOS=QUEUE(QFIRST)
 
        IF(QFIRST.EQ.NP) THEN
           QFIRST=3
        ELSE
           QFIRST=QFIRST+1
        ENDIF
 
        QUEUE(1)=QFIRST
 
        RETURN
        END
 
!---------------------------------------------------------------------
      SUBROUTINE FIFO_EMPTY(NP,QUEUE,INFO)
 
      IMPLICIT none
 
        INTEGER NP,INFO,QFIRST,QLAST
        INTEGER QUEUE(NP*NP)
 
        QFIRST=QUEUE(1)
        QLAST=QUEUE(2)
 
        INFO=0
        IF(QFIRST.EQ.QLAST) INFO=1
 
        RETURN
        END

