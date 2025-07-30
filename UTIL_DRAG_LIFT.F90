!---------------------------------------------

   SUBROUTINE open_drag_files()
    
    USE global_parameters
    USE flow_parameters
    
    IMPLICIT NONE
    INTEGER :: ibody
    
    CHARACTER*9          :: dragfile
    CHARACTER*25         :: indragfile
    
    !Open drag and lift file

    DO ibody = 1, nbody
      dragfile = TRIM("drag_lift")
      WRITE(indragfile,101) dragfile,ibody
101   FORMAT(a,'_body_',i3.3,'.dat')
      IF (nread==0) THEN
         OPEN(UNIT=ifuDragOut+ibody-1,FILE=indragfile,FORM='formatted',ACTION="write")
      ELSE  
         OPEN(UNIT=ifuDragOut+ibody-1,FILE=indragfile,FORM='formatted',POSITION='append',ACTION="write")
      ENDIF
    END DO ! ibody

   END SUBROUTINE open_drag_files
!-------------------------------------------------------------------------------
   SUBROUTINE compute_marker_stress()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE boundary_arrays
    USE MPI_module
    USE fsi_module

    IMPLICIT NONE

    INTEGER           :: i,j,k,m,iBody,nElems
    INTEGER           :: node(3),n,nMarker

    REAL(KIND=CGREAL) :: xNorm, yNorm, zNorm, sign
    REAL(KIND=CGREAL) :: xEc, yEc, zEc
    REAL(KIND=CGREAL) :: dUtdn(3), dUtdn2(3), pEc, pEc2 
    REAL(KIND=CGREAL) :: fdotN
    REAL(KIND=CGREAL) :: xpf,ypf,zpf,xsf,ysf,zsf,elemArea
    REAL(KIND=CGREAL) :: force_t(9),force_i(9)
    REAL(KIND=CGREAL) :: beta
!*****************************************************


   ! total forces 
   cxt(1:nBody) = 0.0_CGREAL
   cyt(1:nBody) = 0.0_CGREAL
   czt(1:nBody) = 0.0_CGREAL

   cxs(1:nBody) = 0.0_CGREAL
   cys(1:nBody) = 0.0_CGREAL
   czs(1:nBody) = 0.0_CGREAL

   cmxt(1:nBody)= 0.0_CGREAL
   cmyt(1:nBody)= 0.0_CGREAL
   cmzt(1:nBody)= 0.0_CGREAL

   cmxs(1:nBody)= 0.0_CGREAL
   cmys(1:nBody)= 0.0_CGREAL
   cmzs(1:nBody)= 0.0_CGREAL

   !zero out marker forces
   xMarkerShear = 0.0d0
   yMarkerShear = 0.0d0
   zMarkerShear = 0.0d0

   xMarkerForce = 0.0d0
   yMarkerForce = 0.0d0
   zMarkerForce = 0.0d0

   xMarkerStress = 0.0_CGREAL
   yMarkerStress = 0.0_CGREAL
   zMarkerStress = 0.0_CGREAL
   MarkerArea    = 0.0_CGREAL
 
   dotP = 0
   dotP2 = 0
   
   DO iBody = 1, nBody   
     xyzMarkerForceTmp  (:) = 0.0d0
     xyzMarkerForceTmp_i(:) = 0.0d0

     nElems = totNumTriElem(iBody)

     DO m = 1, nElems
        xEc = triElemCentX(m,iBody)
        yEc = triElemCentY(m,iBody)
        zEc = triElemCentZ(m,iBody) ! coordinates of element centroid

        xNorm = triElemNormx(m,iBody)
        yNorm = triElemNormy(m,iBody)
        zNorm = triElemNormz(m,iBody)

        elemArea= triElemArea(m,iBody)

        node(1) = triElemNeig(1,m,iBody)
        node(2) = triElemNeig(2,m,iBody)
        node(3) = triElemNeig(3,m,iBody)

        sign = -1.0_CGREAL
        !print*, 'Entering marker_interpolation ... '

        !IF (Flag_outside_Marker(node(1),iBody)*Flag_outside_Marker(node(2),iBody) &
        !     *Flag_outside_Marker(node(3),iBody)==1)then
        !   dUtdn(1:3) = 0.0_CGREAL
        !   pEc        = 0.0_CGREAL
        !   dUtdn2(1:3) = 0.0_CGREAL
        !   pEc2        = 0.0_CGREAL
        !   goto 2010
        !endif
        
        ! if the element center is within the currect subdomain, compute stresses.
        if( zEC .le. z(z_end) .and. zEC .gt. z(z_start) .and. &
            yEC .le. y(y_end) .and. yEC .gt. y(y_start) .and. &
            xEC .le. xc(nx-1) .and. xEC .gt. xc(1)      ) then
           !Ye,debug
           !if (  (ntime.eq.16645).and.(iter_FSI.eq.0).and.(lProc.eq.7).and.  &
           !      (iBody.eq.2).and.(m.eq.973)  ) then
           !   write(401,*)xEc,yEc,zEc
           !   write(401,*)z(z_start),z(z_end)
           !   write(401,*)xNorm,yNorm,zNorm
           !endif

           CALL marker_interpolation(iBody,m,xEc,yEc,zEc,xNorm,yNorm,zNorm,  &
                                     sign,dUtdn,pEc)        
           !write(*,'(a,I5,3F10.5,4E12.5)')'Marker stress:', lProc,xEc,yEc,zEC,dUtdn(1:3),pEc

           !For membrane structure, compute stresses on the other side as well.
!           if(unstruc_surface_type(iBody)==MEMBRANE) then
!              ! compute the stresses on the other side
!              sign = 1.0_CGREAL
!              CALL marker_interpolation(iBody,xEc,yEc,zEc,xNorm,yNorm,zNorm,  &
!                                        sign,dUtdn2,pEc2)
!           else
              dUtdn2(1:3) = 0.0_CGREAL
              pEc2        = 0.0_CGREAL
!           endif
        else
           dUtdn(1:3) = 0.0_CGREAL
           pEc        = 0.0_CGREAL
           dUtdn2(1:3)= 0.0_CGREAL
           pEc2       = 0.0_CGREAL           
        endif ! zEC

2010    continue   !Fang-Bao Tian

        !write(*,'(I5,20F10.5)')m,pEc,pEc2,dUtdn(1:3),dUtdn2(1:3)

        ! compute difference between two sides of the surface
        dUtdn(1:3) = dUtdn(1:3) - dUtdn2(1:3)
        pEc        = pEc - pEc2
        
        ! Remove the normal component in the shear stress
        ! It doesn't matter whether the normal points to inside or outside of the body
        fdotN = dUtdn(1)*xnorm + dUtdn(2)*ynorm + dUtdn(3)*znorm
        dUtdn(1) = dUtdn(1) - fdotN * xnorm
        dUtdn(2) = dUtdn(2) - fdotN * ynorm
        dUtdn(3) = dUtdn(3) - fdotN * znorm

        xpf      = pEc * xnorm*elemArea
        ypf      = pEc * ynorm*elemArea
        zpf      = pEc * znorm*elemArea

        xsf      = dUtdn(1)*reinv*elemArea
        ysf      = dUtdn(2)*reinv*elemArea
        zsf      = dUtdn(3)*reinv*elemArea

        !Ye,debug
        !if (  (iBody.eq.2).and.( (node(1).eq.581).or.  &
        !     (node(2).eq.581).or.(node(3).eq.581) )  )then
        !  write(300+lProc,*)ntime,time,iter_FSI,lProc,m
        !  write(300+lProc,'(1(E16.8,1X))')pEc
        !  write(300+lProc,'(3(E16.8,1X))')dUtdn 
        !  write(300+lProc,'(3(E16.8,1X))')xpf, ypf, zpf
        !  write(300+lProc,'(3(E16.8,1X))')xsf, ysf, zsf
        !  write(300+lProc,*)'======'
        !endif

        ! total forces
        cxt(iBody) = cxt(iBody) + xpf + xsf
        cyt(iBody) = cyt(iBody) + ypf + ysf
        czt(iBody) = czt(iBody) + zpf + zsf

        cmxt(iBody) = cmxt(iBody) + yEc*(zsf+zpf) - zEc*(ysf+ypf)
        cmyt(iBody) = cmyt(iBody) + zEc*(xsf+xpf) - xEc*(zsf+zpf)
        cmzt(iBody) = cmzt(iBody) + xEc*(ysf+ypf) - yEc*(xsf+xpf)

        ! total shear forces
        cxs(iBody) = cxs(iBody) + xsf
        cys(iBody) = cys(iBody) + ysf
        czs(iBody) = czs(iBody) + zsf
      
        ! shear and total forces at the marker point
        do i=1,3
           n = node(i)
           MarkerArea  (n,iBody) = MarkerArea  (n,iBody) + elemArea/3.0_CGREAL

           xMarkerShear(n,iBody) = xMarkerShear(n,iBody) + xsf/3.0_CGREAL
           yMarkerShear(n,iBody) = yMarkerShear(n,iBody) + ysf/3.0_CGREAL
           zMarkerShear(n,iBody) = zMarkerShear(n,iBody) + zsf/3.0_CGREAL

           xMarkerForce(n,iBody) = xMarkerForce(n,iBody) + (xpf+xsf)/3.0_CGREAL
           yMarkerForce(n,iBody) = yMarkerForce(n,iBody) + (ypf+ysf)/3.0_CGREAL
           zMarkerForce(n,iBody) = zMarkerForce(n,iBody) + (zpf+zsf)/3.0_CGREAL
        enddo

     ENDDO! end m  nElems

     !Ye,debug, erroneous marker force interpolation
     call MPI_ALLREDUCE(dotP(1,iBody), dotP2(1,iBody), nPtsMax, &
          MPI_INTEGER, MPI_SUM, flow_comm, ierr)
     !if ( (ntime.eq.2).and.(iBody.eq.2) )then
     !   write(600+lProc,'(4(I8,1X))')ntime,iter_fsi,lProc,m
     !   write(600+lProc,'(3(I8,1X))')dotP2(4999,iBody), &
     !                                dotP2(3583,iBody), &
     !                                dotP2(4996,iBody)
     !   write(600+lProc,*)'======='
     !endif
     do i=1,nPtsMax
        if (dotP2(i,iBody).lt.0)then
           dotP2(i,iBody) = -1
        endif
     enddo

     !Ye,debug
     !if (iBody.eq.2)then
     !   write(600+lProc,'(i6,1X,1(E16.8,1X),i4,1X,3(E16.8,1X))')ntime,time,iter_FSI, &
     !        xMarkerForce(4401,iBody),yMarkerForce(4401,iBody),zMarkerForce(4401,iBody)
     !endif

     !if(abs(cxt(iBody)) > 0.0d0) write(*,'(a,I5,F15.8)')'Before MPI sum:',lProc,cxt(iBody)

     !Combine data from all processors using MPI_SUM. Use a temporary variable to
     !catch the sum, as sometimes using the same variable is not allowed.
     force_i(1) = cxt (iBody);  force_i(2) = cyt (iBody); force_i(3) = czt (iBody);
     force_i(4) = cmxt(iBody);  force_i(5) = cmyt(iBody); force_i(6) = cmzt(iBody);
     force_i(7) = cxs (iBody);  force_i(8) = cys (iBody); force_i(9) = czs (iBody);

     call MPI_ALLREDUCE(force_i, force_t, 9, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, FLOW_COMM, ierr)

     cxt (iBody) = force_t(1);  cyt (iBody) = force_t(2);  czt (iBody) = force_t(3);
     cmxt(iBody) = force_t(4);  cmyt(iBody) = force_t(5);  cmzt(iBody) = force_t(6);
     cxs (iBody) = force_t(7);  cys (iBody) = force_t(8);  czs (iBody) = force_t(9);
     !if(abs(cxt(iBody)) > 0.0d0) write(*,'(a,I5,F15.8)')'After MPI sum:',lProc,cxt(iBody)

     !pack the marker stress data
     nMarker = nPtsBodyMarker(iBody)
     DO m=1,nMarker
        xyzMarkerForceTmp_i(m)           = xMarkerForce(m,iBody)
        xyzMarkerForceTmp_i(m+nMarker  ) = yMarkerForce(m,iBody)
        xyzMarkerForceTmp_i(m+nMarker*2) = zMarkerForce(m,iBody)
        xyzMarkerForceTmp_i(m+nMarker*3) = xMarkerShear(m,iBody)
        xyzMarkerForceTmp_i(m+nMarker*4) = yMarkerShear(m,iBody)
        xyzMarkerForceTmp_i(m+nMarker*5) = zMarkerShear(m,iBody)
     ENDDO
     call MPI_ALLREDUCE(xyzMarkerForceTmp_i, xyzMarkerForceTmp, nMarker*6, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, FLOW_COMM, ierr)

     !unpack the data
     DO m=1,nMarker
        xMarkerForce(m,iBody) = xyzMarkerForceTmp(m)
        yMarkerForce(m,iBody) = xyzMarkerForceTmp(m+nMarker  )
        zMarkerForce(m,iBody) = xyzMarkerForceTmp(m+nMarker*2)

        xMarkerShear(m,iBody) = xyzMarkerForceTmp(m+nMarker*3)
        yMarkerShear(m,iBody) = xyzMarkerForceTmp(m+nMarker*4)
        zMarkerShear(m,iBody) = xyzMarkerForceTmp(m+nMarker*5)

        !Distributed stress at the marker points
        xMarkerStress(m,iBody)= xMarkerForce(m,iBody) / MarkerArea(m,iBody)
        yMarkerStress(m,iBody)= yMarkerForce(m,iBody) / MarkerArea(m,iBody)
        zMarkerStress(m,iBody)= zMarkerForce(m,iBody) / MarkerArea(m,iBody)
     ENDDO

   ENDDO  ! end iBody

   !Ye,debug pressure force and shear force
   !if (lProc==PROC_M)then
   !   write(405,*)ntime,iter_FSI
   !   write(405,*)xMarkerForce( 581,2),yMarkerForce( 581,2),zMarkerForce( 581,2)
      !write(405,*)xMarkerShear(2083,2),yMarkerShear(2083,2),zMarkerShear(2083,2)
   !   write(405,*)xMarkerForce(1542,2),yMarkerForce(1542,2),zMarkerForce(1542,2)
      !write(405,*)xMarkerShear(5252,2),yMarkerShear(5252,2),zMarkerShear(5252,2)
   !   write(405,*)'======='
   !endif

   beta = fsi_filter_f

   DO iBody = 1, nBody
      nElems = nPtsBodyMarker(iBody)
      DO m = 1, nElems
         !Ye,store the unfiltered marker force
         xMarkerForce0(m,iBody) = xMarkerForce(m,iBody)
         yMarkerForce0(m,iBody) = yMarkerForce(m,iBody)
         zMarkerForce0(m,iBody) = zMarkerForce(m,iBody)

         !Ye, filter the marker force
         xMarkerForce(m,iBody) = (1.0d0 - beta)*xMarkerForceOld(m,iBody) + beta*xMarkerForce(m,iBody)
         yMarkerForce(m,iBody) = (1.0d0 - beta)*yMarkerForceOld(m,iBody) + beta*yMarkerForce(m,iBody)
         zMarkerForce(m,iBody) = (1.0d0 - beta)*zMarkerForceOld(m,iBody) + beta*zMarkerForce(m,iBody)
     ENDDO
   ENDDO

!   IF ( MOD(ntime,nmonitor) == 0 ) THEN
!      DO iBody = 1, nBody
!         write(*,'(a)') ' Max local force and shear:'
!         write(*,'(I5,6F12.5)')iBody,maxval(abs(xMarkerForce(:,:))),maxval(abs(yMarkerForce(:,:))),&
!                                     maxval(abs(zMarkerForce(:,:))),maxval(abs(xMarkerShear(:,:))),&
!                                     maxval(abs(yMarkerShear(:,:))),maxval(abs(zMarkerShear(:,:)))

         !write(*,'(a)') 'Drag and lift:'
         !WRITE(*,121)time,cxt(iBody),cyt(iBody),czt(iBody), &
         !                 cxs(iBody),cys(iBody),czs(iBody)
!      ENDDO  ! end iBody    
!   ENDIF

    !Ye,use constant force at specific timestep for test, FOR OUTPUT
    !For iBody=2 only(valve),update the force on Aortic side only
    !if (time.ge.10000.0d0)then
    !   do i=1,10482
    !      if (v_num(i).eq.1)then
    !         !v_num=1:constant force
    !         xMarkerForce(i,2) = xMarkerForceOld(i,2)
    !         yMarkerForce(i,2) = yMarkerForceOld(i,2)
    !         zMarkerForce(i,2) = zMarkerForceOld(i,2)
    !      endif
    !   enddo
    !endif

121 FORMAT(20(3x,1PE12.5))

   END SUBROUTINE  compute_marker_stress
!-------------------------------------------------------------------------------
   SUBROUTINE write_total_force()

    USE global_parameters
    USE flow_parameters
    USE grid_arrays
    USE flow_arrays
    USE boundary_arrays
    use fsi_module

    IMPLICIT NONE

    INTEGER           :: iBody
!*****************************************************
   DO iBody = 1, nBody
     write(*,'(a)') ' Max local force and shear:'
     write(*,'(I5,6E15.5)')iBody,maxval(abs(xMarkerForce(:,iBody))),maxval(abs(yMarkerForce(:,iBody))),&
                                 maxval(abs(zMarkerForce(:,iBody))),maxval(abs(xMarkerShear(:,iBody))),&
                                 maxval(abs(yMarkerShear(:,iBody))),maxval(abs(zMarkerShear(:,iBody)))

!     write(*,'(a)') 'Drag and lift:'
!     WRITE(*,121)time,cxt(iBody),cyt(iBody),czt(iBody), &
!                      cxs(iBody),cys(iBody),czs(iBody)

     if(is_fsi_cvg .eq.1)WRITE(ifuDragOut+iBody-1,121)time,cxt(iBody),cyt(iBody),czt(iBody), &
                              cxs(iBody),cys(iBody),czs(iBody),cmxt(ibody),cmyt(ibody),cmzt(ibody)
   ENDDO  ! end iBody    

121 FORMAT(20(3x,1PE12.5))

   END SUBROUTINE write_total_force
!-------------------------------------------------------------------------------
   SUBROUTINE marker_interpolation(iBody,m,xEc0,yEc0,zEc0,xNorm,yNorm,zNorm,  &
                                   sign,dUdn,pEc)

     USE global_parameters
     USE flow_parameters
     USE grid_arrays
     USE flow_arrays
     USE boundary_arrays
     USE GCM_arrays
     USE MPI_module
     USE fsi_module

     IMPLICIT NONE

     INTEGER          , INTENT (IN) :: iBody,m
     REAL(KIND=CGREAL), INTENT (IN) :: xEc0,yEc0,zEc0
     REAL(KIND=CGREAL), INTENT (IN) :: xNorm,yNorm,zNorm
     REAL(KIND=CGREAL), INTENT (IN) :: sign
     REAL(KIND=CGREAL), INTENT (OUT):: dUdn(3), pEc

     INTEGER :: i,j,k,ii,jj,kk,iRow, info
     INTEGER :: n, nG, iG, jG, kG
     INTEGER :: j1, j2, k1, k2, iMin, iMax, jMin, jMax, kMin, kMax

     REAL(KIND=CGREAL) :: xEc, yEc, zEc
     REAL(KIND=CGREAL) :: xIP, yIP, zIP, uIP, vIP, wIP
     REAL(KIND=CGREAL) :: alphaX,alphaY,alphaZ,minDelta
     REAL(KIND=CGREAL) :: dist, distMin, dist_int, p_inv_dist_int

     REAL(KIND=CGREAL) :: dot_temp
     !**********************************************************
     ii = -1
     jj = -1
     kk = -1

     pEc = 0
     dUdn = 0

     !Y-direction
     if(jProc .gt. 0) then
        j1 = yc_start-2
     else
        j1 = 1       ! first slab
     endif
     if(jProc .lt. nProcY-1) then       
        j2 = yc_end + 2
     else
        j2 = ny-1    ! last slab
     endif

     !Z-direction
     if(kProc .gt. 0) then
        k1 = zc_start-2
     else
        k1 = 1       ! first slab
     endif
     if(kProc .lt. nProcZ-1) then 
        k2 = zc_end + 2
     else
        k2 = nz-1    ! last slab
     endif

     ! for membrane-type structure, offset the marker point
     ! position by the virtual thickness
     if(unstruc_surface_type(iBody)==MEMBRANE) then
        xEc = xEc0 + sign*xNorm*membrane_tkns/2.0_CGREAL
        yEc = yEc0 + sign*yNorm*membrane_tkns/2.0_CGREAL
        zEc = zEc0 + sign*zNorm*membrane_tkns/2.0_CGREAL
     else
        xEc = xEc0
        yEc = yEc0 
        zEc = zEc0 
     endif

     ! find closest dead cell
     distMin = 1.0E8_CGREAL
     DO n=1,nDead
        i = iDead(n)
        j = jDead(n)
        k = kDead(n)

        dist = (xc(i)-xEc)**2 + (yc(j)-yEc)**2  + (zc(k)-zEc)**2

        IF ( dist <= distMin ) THEN
           distMin = dist
           nG  = n
           iG  = i
           jG  = j
           kG  = k
        ENDIF
     ENDDO  ! end n

     !Ye,debug
     !if (  (iBody.eq.2).and. (m.eq.754)  )then
     !   write(502,*)ntime,iter_fsi,m
     !   write(502,'(3F12.7)')xEc,yEc,zEc
     !   write(502,*)nG, iG, jG, kG
     !   write(502,*)iDeadCellIndex(nG),jDeadCellIndex(nG),kDeadCellIndex(nG)
     !   write(502,*)'===='
     !endif

     ! search vicinity of the dead cell to find nodes surrounding element
     iMin = Max(1,iG -20); iMax = Min(nx-1,iG+20)
     jMin = Max(j1,jG-20); jMax = Min(j2,jG+20)
     kMin = Max(k1,kG-20); kMax = Min(k2,kG+20)

     !Ye,debug
     !if (  (ntime.eq.16645).and.(iter_FSI.eq.0).and.(lProc.eq.7).and.  &
     !      (iBody.eq.2).and.(m.eq.973)  ) then
     !   write(404,*)iMin,iMax,xc(iMin),xc(iMax)
     !   write(404,*)jMin,jMax,yc(jMin),yc(jMax)
     !   write(404,*)kMin,kMax,zc(kMin),zc(kMax)
     !endif

     DO i = iMin, iMax

        IF ( ( xc(i)-xEc)*(xc(i+1)-xEc ) <= 0.0_CGREAL ) ii = i
     ENDDO
     DO j = jMin, jMax
        IF ( ( yc(j)-yEc)*(yc(j+1)-yEc ) <= 0.0_CGREAL ) jj = j
     ENDDO
     DO k = kMin, kMax
        IF ( ( zc(k)-zEc)*(zc(k+1)-zEc ) <= 0.0_CGREAL ) kk = k
     ENDDO

     if ( ((ii.eq.-1).or.(jj.eq.-1).or.(kk.eq.-1)).and.(iBody.eq.2) )then
     write(*,*)'BodyMarkerForce Warning!',lProc,m,ii,jj,kk,iG,jG,kG,xEc,yEc,zEc
     endif

     !Ye,debug
     !dot_temp=(xc(iG)-xEc0)*xNorm+(yc(jG)-yEc0)*yNorm+(zc(kG)-zEc0)*zNorm
     !if (  (dot_temp.gt.0).or.(-dot_temp.gt.membrane_tkns) )then
     !   dotP(triElemNeig(1,m,iBody),iBody) = -1
     !   dotP(triElemNeig(2,m,iBody),iBody) = -1
     !   dotP(triElemNeig(3,m,iBody),iBody) = -1
        !Ye,debug
     !   if ( (ntime.eq.100).and.(iBody.eq.2) )then
     !   write(500+lProc,'(4(I8,1X))')ntime,iter_fsi,lProc,m
     !   write(500+lProc,'(3(I8,1X))')iG,jG,kG
     !   write(500+lProc,'(3F12.7)')xNorm,yNorm,zNorm
     !   write(500+lProc,'(1F12.7)')dot_temp
     !   write(500+lProc,'(3(I8,1X))')triElemNeig(1,m,iBody), &
     !                                triElemNeig(2,m,iBody), &
     !                                triElemNeig(3,m,iBody)
     !   write(500+lProc,'(3(I8,1X))')dotP(triElemNeig(1,m,iBody),iBody), &
     !                                dotP(triElemNeig(2,m,iBody),iBody), &
     !                                dotP(triElemNeig(3,m,iBody),iBody)
     !   write(500+lProc,*)'======='
     !   endif
     !endif

     !Ye,debug, m: elem #
     !if ( ((ntime.ge.3618).or.(ntime.le.3622)).and.(iBody.eq.2).and.(m.eq.9385) )then
     !   write(500+lProc,'(3(I8,1X))')ntime,iter_fsi,lProc
     !   write(500+lProc,'(3(I8,1X))')iG,jG,kG
     !   write(500+lProc,'(3F12.7)')xNorm,yNorm,zNorm
     !   write(500+lProc,'(1F12.7)')dot_temp
     !   write(500+lProc,'(3(I8,1X))')triElemNeig(1,m,iBody), &
     !                                triElemNeig(2,m,iBody), &
     !                                triElemNeig(3,m,iBody)
     !   write(500+lProc,'(3(I8,1X))')dotP(triElemNeig(1,m,iBody),iBody), &
     !                                dotP(triElemNeig(2,m,iBody),iBody), &
     !                                dotP(triElemNeig(3,m,iBody),iBody)
     !endif

     !Ye,debug
     !if (  (lProc.eq.32).and.  &
     !      (iBody.eq. 2).and.(m.eq.754)  ) then
     !   write(406,*)ntime,iter_FSI
     !   write(406,*)nG,iG,jG,kG
     !   write(406,*)p(iG,jG,kG)
     !   write(406,*)ii,jj,kk
     !   DO k=0,1
     !   DO j=0,1
     !   DO i=0,1
     !      write(406,*)iblank(ii+i,jj+j,kk+k),dead_cell(ii+i,jj+j,kk+k), &
     !                  p(ii+i,jj+j,kk+k)     
     !   ENDDO
     !   ENDDO
     !   ENDDO
     !   write(406,*)'-----'
     !endif

! inverse distance sqr weighted interpolation
     p_inv_dist_int = 0.0_CGREAL
     dist_int       = 0.0_CGREAL

     DO k=0,1
     DO j=0,1
     DO i=0,1
        if ( (ii.eq.-1).or.(jj.eq.-1).or.(kk.eq.-1) ) cycle

        dist = (xc(ii+i)-xEc)**2 &
              +(yc(jj+j)-yEc)**2 &
              +(zc(kk+k)-zEc)**2

        ! include in the interpolation if the cell is a fluid cell or dead cell
        IF ( iblank(ii+i,jj+j,kk+k) == 0 .OR. dead_cell(ii+i,jj+j,kk+k) == -1 ) THEN
           dist_int       = dist_int + (1.0_CGREAL/dist)
           p_inv_dist_int = p_inv_dist_int + (1.0_CGREAL/dist)*p(ii+i,jj+j,kk+k)
        ENDIF
      ENDDO
      ENDDO
      ENDDO

      DO k=0,1
      DO j=0,1
      DO i=0,1

         if ( (ii.eq.-1).or.(jj.eq.-1).or.(kk.eq.-1) ) cycle

         !If the stencil involves iblank cells where flow variables are not defined,
         !we use inverse-distance averaged pressure.
         IF( iblank(ii+i,jj+j,kk+k) == 1 .AND. dead_cell(ii+i,jj+j,kk+k)== 0) THEN
            if (dist_int < 0.000001d0) then   ! no point is found for the interpolation
               !write(*,'(6F12.5)') xc(ii), yc(jj), zc(kk), dist_int
               pEc = 0.0_CGREAL
            else
               pEc = p_inv_dist_int/dist_int
               GOTO 567
            endif
         ENDIF
      ENDDO
      ENDDO
      ENDDO


      if ( (ii.eq.-1).or.(jj.eq.-1).or.(kk.eq.-1) ) goto 567

      !if the stencil only involves valid cells, we switch to 
      ! trilinear interpolation at the centroid.
      alphaX = (xEc-xc(ii))*dxcinv(ii+1)
      alphaY = (yEc-yc(jj))*dycinv(jj+1)
      alphaZ = (zEc-zc(kk))*dzcinv(kk+1)

      pEc  = p(ii,jj,kk)*(1.0_CGREAL - alphaX)*(1.0_CGREAL -alphaY)*(1.0_CGREAL - alphaZ) &
           + p(ii+1,jj,kk)*alphaX*(1.0_CGREAL - alphaY)*(1.0_CGREAL - alphaZ)             &
           + p(ii,jj+1,kk)*(1.0_CGREAL - alphaX)*alphaY*(1.0_CGREAL - alphaZ)             &
           + p(ii,jj,kk+1)*(1.0_CGREAL - alphaX)*(1.0_CGREAL - alphaY)*alphaZ             &
           + p(ii+1,jj+1,kk)*alphaX*alphaY*(1.0_CGREAL - alphaZ)             &
           + p(ii+1,jj,kk+1)*alphaX*(1.0_CGREAL - alphaY)*alphaZ             &
           + p(ii,jj+1,kk+1)*(1.0_CGREAL - alphaX)*alphaY*alphaZ             &
           + p(ii+1,jj+1,kk+1)*alphaX*alphaY*alphaZ

567   CONTINUE

      ! compute shear stress at element centroid
      ! find velocity at image point corresponding to ghost point

     !Ye,debug
     !if (  (ntime.eq.16645).and.(iter_FSI.eq.0).and.(lProc.eq.7).and.  &
     !      (iBody.eq.2).and.(m.eq.973)  ) then
     !   write(406,*)nG
     !   write(406,*)iDeadCellIndex(nG),jDeadCellIndex(nG),kDeadCellIndex(nG)
     !   write(406,*)coeffGCMDeadD(1:8,nG)
     !   write(406,*)probeLengthDead(nG)
     !   write(406,*)uBodyInterceptDead(nG),vBodyInterceptDead(nG),  &
     !               wBodyInterceptDead(nG)
     !endif

      if ( (ii.eq.-1).or.(jj.eq.-1).or.(kk.eq.-1) ) goto 568

      uIP = 0.0_CGREAL
      vIP = 0.0_CGREAL
      wIP = 0.0_CGREAL

      DO iRow = 1,iRowMax
         ii = iDeadCellIndex(nG) + incI(iRow)
         jj = jDeadCellIndex(nG) + incJ(iRow)
         kk = kDeadCellIndex(nG) + incK(iRow)

         IF ( ii /= iG .OR. jj /= jG .OR. kk /= kG) THEN
            uIP = uIP + coeffGCMDeadD(iRow,nG)* u(ii,jj,kk)
            vIP = vIP + coeffGCMDeadD(iRow,nG)* v(ii,jj,kk)
            wIP = wIP + coeffGCMDeadD(iRow,nG)* w(ii,jj,kk)
         ELSE
            uIP = uIP + coeffGCMDeadD(iRow,nG)* uBodyInterceptDead(nG)
            vIP = vIP + coeffGCMDeadD(iRow,nG)* vBodyInterceptDead(nG)
            wIP = wIP + coeffGCMDeadD(iRow,nG)* wBodyInterceptDead(nG)
         ENDIF ! ii

      ENDDO ! iRow
     !---- dUdn -------
     dUdn(1)= (uIP - u(iG,jG,kG)) / probeLengthDead(nG) * (-sign)
     dUdn(2)= (vIP - v(iG,jG,kG)) / probeLengthDead(nG) * (-sign)
     dUdn(3)= (wIP - w(iG,jG,kG)) / probeLengthDead(nG) * (-sign)

568   CONTINUE

     !Ye,debug
     !if ( (ntime.eq.6650).and.llProc.eq.5).and.(iBody.eq.2).and.  &
     !     (m.eq.3447) )then
     !   write(503,*)ntime,iter_fsi,m
     !   write(503,'(3F12.7)')uIP,vIP,wIP
     !   write(503,'(3F12.7)')u(iG,jG,kG),v(iG,jG,kG),w(iG,jG,kG)
     !   write(503,'(3F12.7)')dUdn
     !   write(503,'(1F12.7)')probeLengthDead(nG)
     !   write(503,*)'===='
     !endif

     ! compute the tangential component, dU/dn = (1 - nn) dU/dn
     !dotNorm = dUdn(1)*xNorm + dUdn(2)*yNorm + dUdn(3)*zNorm

     !dUdn(1)= dUdn(1) - dotNorm * xNorm
     !dUdn(2)= dUdn(2) - dotNorm * yNorm
     !dUdn(3)= dUdn(3) - dotNorm * zNorm

     !write(*,'(I5,10F12.5)')m, coeffU(1:9)
     !write(*,'(9F12.5)') xM,yM,zM,uM,vM,wM,uIP,vIP,wIP
     !write(*,'(8F12.5)') dUdn(1:3), pM

   END SUBROUTINE marker_interpolation
!---------------------------------------------------------
  SUBROUTINE pseudo_inverse(A, M, N)

! this routine computes the left pseudo-inverse of A
! of dimension (M,N).
    USE global_parameters
    USE flow_parameters

    IMPLICIT NONE   
   
    INTEGER          ,INTENT(IN)     :: M, N
    REAL(KIND=CGREAL),INTENT(INOUT)  :: A (iRowMax+1, iRowMax)

    !...local variable
    INTEGER            :: i,j,k
    INTEGER            :: N2, LDA, LDU, LDVT, LWORK, INFO
    REAL(KIND=CGREAL)  :: S(N), U(M,M), VT(N,N)  !Note: dummy arguments can be
                                                 ! to define an automatic array.
    REAL(KIND=CGREAL), DIMENSION(:,:), ALLOCATABLE  :: WORK
!**********************************************************
     
    LDA   = iRowMax+1
    LDU   = M
    LDVT  = N

! according the LAPACK manual, the following line for LWORK 
! should be enough. But on some computers I found LWORK has to be
! much larger.
!   LWORK = MAX(3*MIN(M,N)+MAX(M,N), 5*MIN(M,N)-4) 

    LWORK = MAX(3*MIN(M,N)+MAX(M,N), 5*MIN(M,N)-4) * 2

    ALLOCATE(WORK(LWORK, LWORK))
    
!   print*, 'LWORK = ', LWORK

! call the LAPACK routine to compute the SVD decomposition
! A = U * S * V^T. Only first N columns of U are returned.
 
    CALL DGESVD('S', 'A', M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, &
                LWORK, INFO )

    if (INFO < 0) then
       write(STDOUT,*)'pseudo_inverse: INFO < 0, the SVD argument had an illegal value.'
    elseif (INFO > 0) then   
       write(STDOUT,*)'pseudo_inverse: DBDSQR  did  not  converge!' 
    endif

    N2 = N
    DO j=1,N
       IF (S(j) < 1e-8*S(1) ) THEN
          N2 = j-1
          EXIT  ! discard the small singular values
       ENDIF
       DO i=1,M
          U(i,j) = U(i,j) / S(j)
       ENDDO
    ENDDO

    A = 0.0_CGREAL

    DO j=1,M
    DO i=1,N
      DO k=1,N2
        A(j,i) = A(j,i) + VT(k,i) * U(j,k)
      ENDDO  
    ENDDO
    ENDDO

    DEALLOCATE(WORK)       
  END SUBROUTINE pseudo_inverse
!---------------------------------------------------------
