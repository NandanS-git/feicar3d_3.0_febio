!------------------------------------------------------------------------------
      subroutine calculate_arclength_norm_cnt(itime1,istart1,ksub1)

        USE contact

        IMPLICIT NONE

!... Loop variables
        INTEGER           :: iBody,m, mMax, mMin, cNode, cElement
        INTEGER           :: node1, node2, node3
        REAL*8 :: tangX, tangY, distNode, dMinUnstruc
        REAL*8 :: normMag, tangMag
        REAL*8 :: vectorx,vectory,vectorz,vectProduct
        REAL*8 :: triElemVectorX(2),triElemVectorY(2),triElemVectorZ(2)

        INTEGER :: itime1,istart1,ksub1

! check direction of surface normal and adjust ordering of triangle if normal does not point into body
! this is done only during fresh start

        IF( (itime1.eq.istart1).AND.(ksub1.eq.1) ) THEN
          normDirFlag = 1.0d0
          call check_orientation_unstruc_cnt()
          CALL determine_norm_dir_unstruc_cnt()
        ENDIF

        !xNormBM = 0.0d0
        !yNormBM = 0.0d0
        !zNormBM = 0.0d0
 
       !--Sweep surface of each triangular element
        DO m=1,totNumElem

           node1 = ElemNeig(1,m)
           node2 = ElemNeig(2,m)
           node3 = ElemNeig(3,m)

          !--First vector of each element
           triElemVectorx(1)= BM_x(node2)
     &                       -BM_x(node1)
           triElemVectory(1)= BM_y(node2)
     &                       -BM_y(node1)
           triElemVectorz(1)= BM_z(node2)
     &                       -BM_z(node1)

           !--Second vector of each element
           triElemVectorx(2)= BM_x(node3)
     &                       -BM_x(node2)
           triElemVectory(2)= BM_y(node3)
     &                       -BM_y(node2)
           triElemVectorz(2)= BM_z(node3)
     &                       -BM_z(node2)

           !--Normal vector of the element
           ElemNorm_x(m)= (  triElemVectory(1)*triElemVectorz(2)
     &                     - triElemVectorz(1)*triElemVectory(2)  )
           ElemNorm_y(m)= (  triElemVectorz(1)*triElemVectorx(2)
     &                     - triElemVectorx(1)*triElemVectorz(2)  )
           ElemNorm_z(m)= (  triElemVectorx(1)*triElemVectory(2)
     &                     - triElemVectory(1)*triElemVectorx(2)  )

           normMag  = SQRT(  ElemNorm_x(m)**2
     &                     + ElemNorm_y(m)**2
     &                     + ElemNorm_z(m)**2    )

           ! Unit Normal vector
           ElemNorm_x(m) = ElemNorm_x(m)/normMag
           ElemNorm_y(m) = ElemNorm_y(m)/normMag
           ElemNorm_z(m) = ElemNorm_z(m)/normMag

           xNormBM(node1) = xNormBM(node1) + ElemNorm_x(m)
           xNormBM(node2) = xNormBM(node2) + ElemNorm_x(m)
           xNormBM(node3) = xNormBM(node3) + ElemNorm_x(m)

           yNormBM(node1) = yNormBM(node1) + ElemNorm_y(m)
           yNormBM(node2) = yNormBM(node2) + ElemNorm_y(m)
           yNormBM(node3) = yNormBM(node3) + ElemNorm_y(m)

           zNormBM(node1) = zNormBM(node1) + ElemNorm_z(m)
           zNormBM(node2) = zNormBM(node2) + ElemNorm_z(m)
           zNormBM(node3) = zNormBM(node3) + ElemNorm_z(m)

           !--Centroid of the each element
           ElemCentx(m)=(BM_x(node1)+BM_x(node2)+BM_x(node3))/3.0d0
           ElemCenty(m)=(BM_y(node1)+BM_y(node2)+BM_y(node3))/3.0d0
           ElemCentz(m)=(BM_z(node1)+BM_z(node2)+BM_z(node3))/3.0d0

        ENDDO ! m

        DO m=1,nPtsBM
           normMag = SQRT( xNormBM(m)**2  
     &                   + yNormBM(m)**2
     &                   + zNormBM(m)**2 )

           xNormBM(m) = xNormBM(m) /normMag
           yNormBM(m) = yNormBM(m) /normMag
           zNormBM(m) = zNormBM(m) /normMag
        ENDDO

        !IF( (itime1.eq.istart1).AND.(ksub1.eq.1) ) THEN
        !DO m=1,totNumElem
        !write(514,'(3F12.5)')ElemNorm_x(m),ElemNorm_y(m),ElemNorm_z(m)
        !ENDDO
        !ENDIF

!        IF( (itime1.eq.istart1).AND.(ksub1.eq.1) ) THEN
!        DO m=1,nPtsBM
!        write(515,'(6F12.5)')BM_x(m),BM_y(m),BM_z(m),
!     &       xNormBM(m),yNormBM(m),zNormBM(m)
!        ENDDO
!        ENDIF

       END SUBROUTINE calculate_arclength_norm_cnt
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

      SUBROUTINE determine_norm_dir_unstruc_cnt()

        USE contact

        IMPLICIT NONE

!... Loop variables
        INTEGER           :: m,cNode, cElement

!... Local variables
        INTEGER           :: node1,node2,node3
        REAL*8 :: distNode, dMinUnstruc
        REAL*8 :: vectorx,vectory,vectorz,vectProduct
        REAL*8 :: triElemVectorX(2),triElemVectorY(2),
     &                       triElemVectorZ(2)
        REAL*8 :: cTriElemNormx,cTriElemNormy,cTriElemNormz
        REAL*8 :: cTriElemCentx,cTriElemCenty,cTriElemCentz


        cNode = 0

!---Find node closest to outside point
        dMinUnstruc = 1.0E8
        DO m=1,nPtsBM
           distNode= SQRT( (BM_x(m)-pointOutsideX)**2  
     &                    +(BM_y(m)-pointOutsideY)**2  
     &                    +(BM_z(m)-pointOutsideZ)**2   )
          IF(distNode <= dMinUnstruc) THEN
            dMinUnstruc  = distNode
            cNode        = m
          ENDIF
        ENDDO

!---Find element corresponing to closest node
        DO m=1,totNumElem
           IF ( (ElemNeig(1,m) == cNode) .OR.
     &          (ElemNeig(2,m) == cNode) .OR.
     &          (ElemNeig(3,m) == cNode)       ) cElement = m
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

       node1 = ElemNeig(1,cElement)
       node2 = ElemNeig(2,cElement)
       node3 = ElemNeig(3,cElement)


        triElemVectorx(1)= BM_x(node2)
     &                    -BM_x(node1)
        triElemVectory(1)= BM_y(node2) 
     &                    -BM_y(node1)
        triElemVectorz(1)= BM_z(node2) 
     &                    -BM_z(node1)

        triElemVectorx(2)= BM_x(node3)
     &                    -BM_x(node2)
        triElemVectory(2)= BM_y(node3) 
     &                    -BM_y(node2)
        triElemVectorz(2)= BM_z(node3)
     &                    -BM_z(node2)

!-- Normal of closest element
       cTriElemNormx = ( triElemVectory(1)*triElemVectorz(2)  
     &                  -triElemVectorz(1)*triElemVectory(2) )
       cTriElemNormy = ( triElemVectorz(1)*triElemVectorx(2)  
     &                  -triElemVectorx(1)*triElemVectorz(2) )
       cTriElemNormz = ( triElemVectorx(1)*triElemVectory(2)  
     &                  -triElemVectory(1)*triElemVectorx(2) )

!--Centroid of the closest element
        cTriElemCentx = ( BM_x(node1) 
     &                  + BM_x(node2) 
     &                  + BM_x(node3) )/3.0d0
        cTriElemCenty = ( BM_y(node1) 
     &                  + BM_y(node2) 
     &                  + BM_y(node3) )/3.0d0
        cTriElemCentz = ( BM_z(node1) 
     &                  + BM_z(node2) 
     &                  + BM_z(node3) )/3.0d0

       vectorx = cTriElemCentX - pointOutsideX
       vectory = cTriElemCentY - pointOutsideY
       vectorz = cTriElemCentZ - pointOutsideZ

        vectProduct = vectorx*cTriElemNormx
     &              + vectory*cTriElemNormy
     &              + vectorz*cTriElemNormz

!RRRRRRRRRRRRRRRRRRRRRRRRRRRRr
!        print*,'vectProduct = ',vectproduct
!RRRRRRRRRRRRRRRRRRRRRRRRRRRRr

       normDirFlag = 1.0d0

       IF (vectProduct < 0.0d0)  THEN
         normDirFlag =-1.0d0
!         print*,'Reordering triangle vertices since normal points out of the solid body'
         DO m=1,totNumElem         ! changing vertex ordering
           node2 = ElemNeig(2,m)
           node3 = ElemNeig(3,m)
           ElemNeig(2,m) = node3
           ElemNeig(3,m) = node2
         ENDDO
       ENDIF

      END SUBROUTINE determine_norm_dir_unstruc_cnt
!------------------------------------------------------------------------------
      SUBROUTINE check_orientation_unstruc_cnt()

       USE contact

       IMPLICIT NONE

!... Loop variables
       INTEGER           :: m,n,i,j,info

!... Local variables
      INTEGER           :: node1(3),node2(3),v1(2),v2(2),i1,d1,d2,
     &                     nodetemp2,nodetemp3,eci
      INTEGER, DIMENSION(:), ALLOCATABLE  :: queue,iNorm
      INTEGER, DIMENSION(:,:), ALLOCATABLE  :: ElemConnect 
      REAL :: distNode, dMinUnstruc
      REAL :: vectorx,vectory,vectorz,vectProduct
      REAL :: triElemVectorX(2),triElemVectorY(2),triElemVectorZ(2)
      REAL :: cTriElemNormx,cTriElemNormy,cTriElemNormz
      REAL :: cTriElemCentx,cTriElemCenty,cTriElemCentz


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
      allocate(ElemConnect(3,totNumElem))
      allocate(queue(totNumElem))
      allocate(iNorm(totNumElem))

      iNorm=0

      DO m=1,totNumElem

         node1(1) = ElemNeig(1,m)
         node1(2) = ElemNeig(2,m)
         node1(3) = ElemNeig(3,m)

        eci=0
        do n=1,totNumElem

            if( (m.ne.n) )then !.and.(iswitch(m).eq.1 ) 
                
                node2(1) = ElemNeig(1,n)
                node2(2) = ElemNeig(2,n)
                node2(3) = ElemNeig(3,n)
    
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



      call FIFO_INIT_cnt(m,queue)
    
      m=1
      call FIFO_ADD_cnt(m,totNumElem,queue)
      iNorm(m)=1
      call FIFO_EMPTY_cnt(totNumElem,queue,info)          
    
      do while (info.eq.0)

        call FIFO_REMOVE_cnt(m,totNumElem,queue)

        node1(1) = ElemNeig(1,m)
        node1(2) = ElemNeig(2,m)
        node1(3) = ElemNeig(3,m)

        do eci=1,3
            n=ElemConnect(eci,m)

            if( (iNorm(n).eq.0) .and. (n.gt.0) ) then
                node2(1) = ElemNeig(1,n)
                node2(2) = ElemNeig(2,n)
                node2(3) = ElemNeig(3,n)

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
                    if(((v1(1).eq.1).and.(v1(2).eq.2) ).or. 
     &                 ((v1(1).eq.2).and.(v1(2).eq.3) ).or. 
     &                 ((v1(1).eq.3).and.(v1(2).eq.1) ) )then
                            d1=1
                    endif
                    d2=0
                    if(((v2(1).eq.1).and.(v2(2).eq.2) ).or. 
     &                 ((v2(1).eq.2).and.(v2(2).eq.3) ).or. 
     &                 ((v2(1).eq.3).and.(v2(2).eq.1) ) )then
                            d2=1
                    endif

                    if(d1.eq.d2) then
!                       print*,'Reordering triangle vertices',m,n
                       nodetemp2 = ElemNeig(2,n)
                       nodetemp3 = ElemNeig(3,n)
                       ElemNeig(2,n) = nodetemp3
                       ElemNeig(3,n) = nodetemp2
                    endif
                endif

                iNorm(n)=1
                call FIFO_ADD_cnt(n,totNumElem,queue)

            endif


        enddo

        call FIFO_EMPTY_cnt(totNumElem,queue,info)

        if (info.eq.1) then
        DO m=1,totNumElem
            if( iNorm(m).eq.0 ) then
                !Ye
                !print *,'Element', m,' unchanged'
                call FIFO_ADD_cnt(m,totNumElem,queue)
                info=0
                goto 999
            endif
        enddo
999     continue
        endif

      enddo      

      DO m=1,totNumElem
         !Ye
         !if( iNorm(m).eq.0 ) print *,'Element', m,'remained unchanged!!'
      enddo

      END SUBROUTINE check_orientation_unstruc_cnt

   
!---------------------------------------------------------------------
      SUBROUTINE FIFO_INIT_cnt(NP,QUEUE)
 
      IMPLICIT none
 
        INTEGER NP
        INTEGER QUEUE(NP)
 
 
        QUEUE(1)=3
        QUEUE(2)=3
        RETURN
      END
 
!---------------------------------------------------------------------
      SUBROUTINE FIFO_ADD_cnt(IPOS,NP,QUEUE)
 
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
      SUBROUTINE FIFO_REMOVE_cnt(IPOS,NP,QUEUE)
 
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
      SUBROUTINE FIFO_EMPTY_cnt(NP,QUEUE,INFO)
 
      IMPLICIT none
 
        INTEGER NP,INFO,QFIRST,QLAST
        INTEGER QUEUE(NP*NP)
 
        QFIRST=QUEUE(1)
        QLAST=QUEUE(2)
 
        INFO=0
        IF(QFIRST.EQ.QLAST) INFO=1
 
        RETURN
        END

